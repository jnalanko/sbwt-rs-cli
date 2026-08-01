// Code by Martin Kostadinov.

/// Prefer to use the [define_variants] macro in order to ensure all of the variants have a
/// corresponding value in the enum defined by [define_variants_enum].
macro_rules! define_scan {
    ($name:ident, $word:ty, $element:ty) => {
        #[derive(Clone, Debug, PartialEq, Eq)]
        pub struct $name {
            pub words: Vec<$word>,
            pub n: usize,
        }

        impl $name {
            const ZERO: [$element; Self::LANES] = [0; Self::LANES];
        }

        impl Scan for $name {
            type Word = $word;
            type Element = $element;
            const LANES: usize = Self::Word::LANES as usize;
            const BYTES_PER_ELEMENT: usize = (Self::Word::BITS as u32 / u8::BITS) as usize / Self::LANES;

            fn from_iterator<T, I>(input: I, n: usize, k: usize) -> Self
            where
                T: Into<usize>,
                I: Iterator<Item = T>,
            {
                assert!(crate::vodbg::util::bit_width(k) <= Self::Element::BITS as usize);

                let ten_percent = n / 10;
                let mut border = ten_percent;
                let mut percent_count = 0;

                #[allow(clippy::manual_div_ceil)]
                let mut words = Vec::with_capacity((n + Self::LANES - 1) / Self::LANES);
                let mut array = Self::ZERO;
                for (i, item) in input.enumerate() {
                    if i != 0 && i % Self::LANES == 0 {
                        words.push(Self::Word::new(array));
                        array = Self::ZERO;
                    }
                    array[i % Self::LANES] = item.into() as $element;

                    if i >= border {
                        border += ten_percent;
                        percent_count += 1;
                        log::info!("[{}] {}0%", stringify!($name), percent_count);
                    }
                }
                words.push(Self::Word::new(array));

                Self { words, n }
            }

            fn scan_left(&self, index: usize, target_length: usize) -> usize {
                if index >= self.n {
                    return 0;
                }

                let target_length = target_length as Self::Element;
                let word_index = index / Self::LANES;
                let index_in_word = index % Self::LANES;

                // Scan the values in the SIMD word the index is located in individually.
                let near_word = self.words[word_index].as_array();
                for i in (0..=index_in_word).rev() {
                    if near_word[i] < target_length {
                        return word_index * Self::LANES + i;
                    }
                }

                // Scan the rest of the words using SIMD operations.
                for w in (0..word_index).rev() {
                    let comparison_result = self.words[w].simd_lt(target_length);
                    if !comparison_result.any() {
                        continue;
                    }
                    let bitmask = comparison_result.to_bitmask();
                    let rightmost_smaller_element = u32::BITS as usize - 1 - bitmask.leading_zeros() as usize;
                    return w * Self::LANES + rightmost_smaller_element;
                }

                0
            }

            fn scan_right(&self, index: usize, target_length: usize) -> usize {
                if index >= self.n {
                    return self.n;
                }

                let target_length = target_length as Self::Element;
                let word_index = index / Self::LANES;
                let index_in_word = index % Self::LANES;

                // Similarly to the scan_left procedure, first scan values in the word the index is located
                // in individually.
                let near_word = self.words[word_index].as_array();
                for i in index_in_word..Self::LANES {
                    if near_word[i] < target_length {
                        return word_index * Self::LANES + i;
                    }
                }

                // Then scan the rest of the words using SIMD.
                for w in (word_index + 1)..self.words.len() {
                    let comparison_result = self.words[w].simd_lt(target_length);
                    if !comparison_result.any() {
                        continue;
                    }
                    let bitmask = comparison_result.to_bitmask();
                    let leftmost_smaller_element = bitmask.trailing_zeros() as usize;
                    let result = w * Self::LANES + leftmost_smaller_element;
                    return result.min(self.n);
                }

                self.n
            }

            /// Bound gives the maximum number of words to be scanned in addition to the one the index is located.
            fn scan_left_bounded(&self, index: usize, target_length: usize, bound: usize) -> Result<usize, usize> {
                if index >= self.n {
                    return Err(self.n);
                }

                let target_length = target_length as Self::Element;
                let word_index = index / Self::LANES;
                let index_in_word = index % Self::LANES;

                // Scan the values in the SIMD word the index is located in individually.
                let near_word = self.words[word_index].as_array();
                for i in (0..=index_in_word).rev() {
                    if near_word[i] < target_length {
                        return Ok(word_index * Self::LANES + i);
                    }
                }

                let lower_bound_word_index = word_index.saturating_sub(bound);

                // Scan the rest of the words using SIMD operations.
                for w in (lower_bound_word_index..word_index).rev() {
                    let comparison_result = self.words[w].simd_lt(target_length);
                    if !comparison_result.any() {
                        continue;
                    }
                    let bitmask = comparison_result.to_bitmask();
                    let rightmost_smaller_element = u32::BITS as usize - 1 - bitmask.leading_zeros() as usize;
                    return Ok(w * Self::LANES + rightmost_smaller_element);
                }

                Err(lower_bound_word_index * Self::LANES)
            }

            /// Bound gives the maximum number of words to be scanned in addition to the one the index is located.
            fn scan_right_bounded(&self, index: usize, target_length: usize, bound: usize) -> Result<usize, usize> {
                if index >= self.n {
                    return Err(self.n);
                }

                let target_length = target_length as Self::Element;
                let word_index = index / Self::LANES;
                let index_in_word = index % Self::LANES;

                // Similarly to the scan_left procedure, first scan values in the word the index is located
                // in individually.
                let near_word = self.words[word_index].as_array();
                for i in index_in_word..Self::LANES {
                    if near_word[i] < target_length {
                        return Ok(word_index * Self::LANES + i);
                    }
                }

                let mut upper_bound_word_index = (word_index + bound).min(self.words.len());
                if upper_bound_word_index < self.words.len() {
                    upper_bound_word_index += 1;
                }

                // Then scan the rest of the words using SIMD.
                for w in (word_index + 1)..upper_bound_word_index {
                    let comparison_result = self.words[w].simd_lt(target_length);
                    if !comparison_result.any() {
                        continue;
                    }
                    let bitmask = comparison_result.to_bitmask();
                    let leftmost_smaller_element = bitmask.trailing_zeros() as usize;
                    let result = w * Self::LANES + leftmost_smaller_element;
                    return Ok(result.min(self.n));
                }

                Err(upper_bound_word_index * Self::LANES - 1)
            }

            fn serialize<W: std::io::Write>(&self, out: &mut W) -> std::io::Result<usize> {
                log::info!("[{}::serialize] serializing...", stringify!($name));
                let mut written: usize = 0;
                out.write_all(&(self.n as u64).to_le_bytes())?;
                written += size_of::<u64>();
                for word in &self.words {
                    for element in word.to_array() {
                        out.write_all(&element.to_le_bytes())?;
                    }
                }
                written += self.words.len() * Self::LANES * Self::BYTES_PER_ELEMENT;
                log::info!("[{}::serialize] done.", stringify!($name));
                Ok(written)
            }

            fn load<R: std::io::Read>(input: &mut R) -> std::io::Result<Self> {
                log::info!("[{}::load] loading...", stringify!($name));
                use byteorder::{ReadBytesExt, LittleEndian};
                let n = input.read_u64::<LittleEndian>()? as usize;
                #[allow(clippy::manual_div_ceil)]
                let word_count = (n + Self::LANES - 1) / Self::LANES;
                let mut words: Vec<Self::Word> = Vec::with_capacity(word_count);
                let mut array = Self::ZERO;
                let mut bytes = [0; Self::BYTES_PER_ELEMENT];
                for _ in 0..word_count {
                    for i in 0..Self::LANES {
                        input.read_exact(&mut bytes)?;
                        array[i] = <$element>::from_le_bytes(bytes);
                    }
                    words.push(Self::Word::new(array));
                }
                let result = Self {
                    words,
                    n
                };
                log::info!("[{}::load] done.", stringify!($name));
                Ok(result)
            }

            #[inline]
            fn word_count(&self) -> usize {
                self.words.len()
            }

            #[inline]
            fn len(&self) -> usize {
                self.n
            }

            #[inline]
            fn is_empty(&self) -> bool {
                self.n == 0
            }
        }
    };
}

/// Prefer to use the [define_variants] macro in order to ensure all of the variants have a
/// corresponding value in the enum defined by [define_variants_enum].
macro_rules! define_variants_enum {
    ($name:ident, $($variant:ident),+ $(,)?) => {
        #[derive(Clone, Debug, PartialEq, Eq)]
        pub enum $name {
            $($variant($variant),)+
        }

        impl $name {
        }

        $(
            impl From<$variant> for $name {
                fn from(variant: $variant) -> Self {
                    Self::$variant(variant)
                }
            }
        )+

        #[allow(unused)]
        impl Scan for $name {
            type Word = ();
            type Element = usize;
            const LANES: usize = 0;
            const BYTES_PER_ELEMENT: usize = 0;

            #[allow(unused)]
            fn from_iterator<T, I>(input: I, n: usize, k: usize) -> Self
            where T: Into<usize>, I: Iterator<Item = T> {
                let bit_width = crate::vodbg::util::bit_width(k);
                $(
                    if bit_width <= <$variant as Scan>::Element::BITS as usize {
                        return $variant::from_iterator(input, n, k).into();
                    }
                )+
                panic!("No variants' elements can store the maximum value of k.");
            }

            fn scan_left(&self, index: usize, target_length: usize) -> usize {
                match self {
                    $(Self::$variant(variant) => variant.scan_left(index, target_length),)+
                }
            }

            fn scan_right(&self, index: usize, target_length: usize) -> usize {
                match self {
                    $(Self::$variant(variant) => variant.scan_right(index, target_length),)+
                }
            }

            fn scan_left_bounded(&self, index: usize, target_length: usize, bound: usize) -> Result<usize, usize> {
                match self {
                    $(Self::$variant(variant) => variant.scan_left_bounded(index, target_length, bound),)+
                }
            }

            fn scan_right_bounded(&self, index: usize, target_length: usize, bound: usize) -> Result<usize, usize> {
                match self {
                    $(Self::$variant(variant) => variant.scan_right_bounded(index, target_length, bound),)+
                }
            }

            fn serialize<W: std::io::Write>(&self, out: &mut W) -> std::io::Result<usize> {
                let mut written = 0;
                match self {
                    $(
                        Self::$variant(variant) => {
                            out.write_all(&($variant::LANES as u64).to_le_bytes())?;
                            out.write_all(&($variant::BYTES_PER_ELEMENT as u64).to_le_bytes())?;
                            written += 2 * size_of::<u64>();
                            written += variant.serialize(out)?;
                            Ok(written)
                        }
                    )+
                }
            }

            fn load<R: std::io::Read>(input: &mut R) -> std::io::Result<Self> {
                use byteorder::{ReadBytesExt, LittleEndian};
                let lanes = input.read_u64::<LittleEndian>()? as usize;
                let bytes_per_element = input.read_u64::<LittleEndian>()? as usize;
                match (lanes, bytes_per_element) {
                    $(
                        ($variant::LANES, $variant::BYTES_PER_ELEMENT) => {
                            let variant = $variant::load(input)?;
                            Ok(variant.into())
                        }
                    )+
                    _ => {
                        use std::io::{Error, ErrorKind};
                        Err(Error::new(
                            ErrorKind::Other,
                            "None of the variants match the number of lanes and the bytes per element."
                        ))
                    }
                }
            }

            #[inline]
            fn lanes(&self) -> usize {
                match self {
                    $(Self::$variant(_) => $variant::LANES,)+
                }
            }

            #[inline]
            fn bytes_per_element(&self) -> usize {
                match self {
                    $(Self::$variant(_) => $variant::BYTES_PER_ELEMENT,)+
                }
            }

            #[inline]
            fn word_count(&self) -> usize {
                match self {
                    $(Self::$variant(variant) => variant.words.len(),)+
                }
            }

            #[inline]
            fn len(&self) -> usize {
                match self {
                    $(Self::$variant(variant) => variant.len(),)+
                }
            }

            #[inline]
            fn is_empty(&self) -> bool {
                match self {
                    $(Self::$variant(variant) => variant.is_empty(),)+
                }
            }
        }
    };
}

macro_rules! define_variants {
    ($enum_name:ident; $($variant:ident, $word:ty, $element:ty $(,)?);+ $(;)?) => {
        $(
            macros::define_scan!($variant, $word, $element);
        )+

        macros::define_variants_enum!($enum_name, $($variant,)+);
    };
}

pub(super) use define_scan;
pub(super) use define_variants_enum;
pub(super) use define_variants;


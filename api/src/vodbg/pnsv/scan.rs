// Code by Martin Kostadinov.

use super::bp::nearest_neighbor_dictionary::NearestNeighbourDictionary as NND;
use super::Pnsv;

pub mod augmented_bounded_scan;
mod macros;
use macros::define_variants;

pub use augmented_bounded_scan as abs;

pub trait Scan: Sized {
    type Word;
    type Element;
    const LANES: usize;
    const BYTES_PER_ELEMENT: usize;
    /// Constructs the structure from an iterator.
    fn from_iterator<T, I>(input: I, n: usize, k: usize) -> Self
    where T: Into<usize>, I: Iterator<Item = T>;
    /// Look at [Pnsv::previous].
    fn scan_left(&self, index: usize, target_length: usize) -> usize;
    /// Look at [Pnsv::next].
    fn scan_right(&self, index: usize, target_length: usize) -> usize;
    /// The same as [Scan::scan_left], however, the number of words which are scanned beyond the
    /// first word (in which the index is located which is scanned element by element) is bounded.
    fn scan_left_bounded(&self, index: usize, target_length: usize, bound: usize) -> Result<usize, usize>;
    /// The same as [Scan::scan_right], however, the number of words which are scanned beyond the
    /// first word (in which the index is located which is scanned element by element) is bounded.
    fn scan_right_bounded(&self, index: usize, target_length: usize, bound: usize) -> Result<usize, usize>;
    /// Serializes this data structure in a binary format and returns the number of bytes written.
    fn serialize<W: std::io::Write>(&self, out: &mut W) -> std::io::Result<usize>;
    /// Loads this data structure from binary.
    fn load<R: std::io::Read>(input: &mut R) -> std::io::Result<Self>;
    /// The number of lanes (elements) in the SIMD word this structure is using.
    fn lanes(&self) -> usize { Self::LANES }
    /// How many bytes each element in the SIMD word has.
    fn bytes_per_element(&self) -> usize { Self::BYTES_PER_ELEMENT }
    /// The number of words in the array. The word count times the number of lanes does not
    /// necessarily equal the number of elements in the LCS array. The LCS array is padded with 0s
    /// to make the number of elements 
    fn word_count(&self) -> usize;
    /// The number of elements in the LCS array without the additional 0s.
    fn len(&self) -> usize;
    fn is_empty(&self) -> bool;
}

impl<T: Scan> Pnsv for T {
    fn previous(&self, index: usize, target_length: usize) -> usize {
        self.scan_left(index, target_length)
    }

    fn next(&self, index: usize, target_length: usize) -> usize {
        self.scan_right(index, target_length)
    }
}

define_variants! {
    LcsSimd;
    LcsSimd8x32, wide::u8x32, u8;
    LcsSimd16x32, wide::u16x32, u16;
    LcsSimd32x16, wide::u32x16, u32;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn definition_order_is_correct() {
        let items: &[u8] = &[0, 1, 2, 3];
        let lcs_simd_1 = LcsSimd::from_iterator(items.iter().cloned(), items.len(), 255);
        let lcs_simd_2 = LcsSimd::from_iterator(items.iter().cloned(), items.len(), 1 << 8);
        let lcs_simd_4 = LcsSimd::from_iterator(items.iter().cloned(), items.len(), 1 << 16);
        assert_eq!(lcs_simd_1.bytes_per_element(), 1);
        assert_eq!(lcs_simd_2.bytes_per_element(), 2);
        assert_eq!(lcs_simd_4.bytes_per_element(), 4);
    }

    macro_rules! serialize_and_load_body {
        ($structure:ty, $max_k:expr) => {{
            let items: &[u8] = &[
                0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
                24, 25, 26, 27, 28, 29, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
            ];
            let lcs_simd = <$structure>::from_iterator(items.iter().cloned(), items.len(), $max_k);
            let mut buffer = Vec::<u8>::new();
            let written = lcs_simd.serialize(&mut buffer).unwrap();
            assert_eq!(buffer.len(), written);
            let lcs_simd_loaded = <$structure>::load(&mut buffer.as_slice()).unwrap();
            assert_eq!(lcs_simd, lcs_simd_loaded);
        }};
    }

    #[test]
    fn serialize_and_load() {
        serialize_and_load_body!(LcsSimd8x32, 31);
        serialize_and_load_body!(LcsSimd16x32, 1<<8);
        serialize_and_load_body!(LcsSimd32x16, 1<<16);
        serialize_and_load_body!(LcsSimd, 31);
        serialize_and_load_body!(LcsSimd, 1<<8);
        serialize_and_load_body!(LcsSimd, 1<<16);
    }

    macro_rules! words_are_correct_body {
        ($structure:ty, $max_k:expr, $element:ty) => {{
            let items: &[u8] = &[
                0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
                24, 25, 26, 27, 28, 29, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
            ];
            let lcs_simd = <$structure>::from_iterator(items.iter().cloned(), items.len(), 31);
            let lanes = <$structure as Scan>::LANES;
            for i in 0..items.len() {
                let item_in_word = lcs_simd.words[i / lanes].as_array()[i % lanes];
                assert_eq!(items[i] as $element, item_in_word);
            }
        }};
    }
    #[test]
    fn words_are_correct() {
        words_are_correct_body!(LcsSimd8x32, 31, u8);
        words_are_correct_body!(LcsSimd16x32, 1<<8, u16);
        words_are_correct_body!(LcsSimd32x16, 1<<16, u32);
    }

    macro_rules! scan_left_body {
        ($structure:ty, $max_k:expr) => {{
            let items: &[u8] = &[
                5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
                5, 5, 5, 5, 5, 5, 5, 5, 5, 1, 5, 5, 5, 5, 5, 5, // 25: 1

                5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
                5, 5, 5, 5, 5, 5, 5, 5, 5, 2, 5, 5, 5, 5, 5, 5, // 57: 2

                5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
                5, 5, 5, 5, 5, 5, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5, // 89: 3

                5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
                5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 5, 5, 5, 5, 5, 5, // 121: 4
            ];

            let lcs_simd = <$structure>::from_iterator(items.iter().cloned(), items.len(), $max_k);
            assert_eq!(
                lcs_simd.previous(125, 5),
                121,
                "Result in the same SIMD word failed."
            );
            assert_eq!(
                lcs_simd.previous(100, 5),
                89,
                "Result in the previous SIMD word failed."
            );
            assert_eq!(lcs_simd.previous(100, 4), 89, "Smaller than 4 failed.");
            assert_eq!(lcs_simd.previous(100, 3), 57, "Smaller than 3 failed.");
            assert_eq!(lcs_simd.previous(100, 2), 25, "Smaller than 2 failed.");
            assert_eq!(lcs_simd.previous(100, 1), 0, "Scan to beginning failed.");
        }};
    }

    #[test]
    fn scan_left() {
        scan_left_body!(LcsSimd8x32, 31);
        scan_left_body!(LcsSimd16x32, 1<<8);
        scan_left_body!(LcsSimd32x16, 1<<16);
        scan_left_body!(LcsSimd, 31);
        scan_left_body!(LcsSimd, 1<<8);
        scan_left_body!(LcsSimd, 1<<16);
    }

    macro_rules! scan_right_body {
        ($structure:ty, $max_k:expr) => {{
            let items: &[u8] = &[
                5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
                5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 5, 5, 5, 5, 5, 5, // 25: 4

                5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
                5, 5, 5, 5, 5, 5, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5, // 57: 3

                5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
                5, 5, 5, 5, 5, 5, 5, 5, 5, 2, 5, 5, 5, 5, 5, 5, // 89: 2

                5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
                5, 5, 5, 5, 5, 5, 5, 5, 5, 1, 5, 5, 5, 5, 5, 5, // 121: 1
            ];

            let lcs_simd = <$structure>::from_iterator(items.iter().cloned(), items.len(), $max_k);
            assert_eq!(
                lcs_simd.next(20, 5),
                25,
                "Result in the same SIMD word failed."
            );
            assert_eq!(
                lcs_simd.next(30, 5),
                57,
                "Result in the next SIMD word failed."
            );
            assert_eq!(lcs_simd.next(30, 4), 57, "Smaller than 4 failed.");
            assert_eq!(lcs_simd.next(30, 3), 89, "Smaller than 3 failed.");
            assert_eq!(lcs_simd.next(30, 2), 121, "Smaller than 2 failed.");
            assert_eq!(lcs_simd.next(30, 1), items.len(), "Scan to the end failed.");
        }};
    }

    #[test]
    fn scan_right() {
        scan_right_body!(LcsSimd8x32, 31);
        scan_right_body!(LcsSimd16x32, 1<<8);
        scan_right_body!(LcsSimd32x16, 1<<16);
        scan_right_body!(LcsSimd, 31);
        scan_right_body!(LcsSimd, 1<<8);
        scan_right_body!(LcsSimd, 1<<16);
    }

    macro_rules! scan_left_bounded_body {
        ($structure:ty, $max_k:expr) => {{
            let items: &[u8] = &[
                5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
                5, 5, 5, 5, 5, 5, 5, 5, 5, 1, 5, 5, 5, 5, 5, 5, // 25: 1

                5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
                5, 5, 5, 5, 5, 5, 5, 5, 5, 2, 5, 5, 5, 5, 5, 5, // 57: 2

                5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
                5, 5, 5, 5, 5, 5, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5, // 89: 3

                5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
                5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 5, 5, 5, 5, 5, 5, // 121: 4
            ];

            let lcs_simd = <$structure>::from_iterator(items.iter().cloned(), items.len(), $max_k);
            let max_bound = lcs_simd.word_count();
            assert_eq!(
                lcs_simd.scan_left_bounded(125, 5, max_bound),
                Ok(121),
                "Result in the same SIMD word failed."
            );
            assert_eq!(
                lcs_simd.scan_left_bounded(100, 5, max_bound),
                Ok(89),
                "Result in the previous SIMD word failed."
            );
            assert_eq!(
                lcs_simd.scan_left_bounded(100, 4, max_bound),
                Ok(89),
                "Smaller than 4 failed."
            );
            assert_eq!(
                lcs_simd.scan_left_bounded(100, 3, max_bound),
                Ok(57),
                "Smaller than 3 failed."
            );
            assert_eq!(
                lcs_simd.scan_left_bounded(100, 2, max_bound),
                Ok(25),
                "Smaller than 2 failed."
            );
            assert!(
                lcs_simd.scan_left_bounded(100, 1, max_bound).is_err(),
                "Scan to beginning failed."
            );

            let bound_multiplier = 32 / lcs_simd.lanes(); 
            assert!(lcs_simd.scan_left_bounded(100, 5, 0 * bound_multiplier).is_err());
            assert!(lcs_simd.scan_left_bounded(100, 4, 0 * bound_multiplier).is_err());
            assert!(lcs_simd.scan_left_bounded(100, 3, 1 * bound_multiplier).is_err());
            assert!(lcs_simd.scan_left_bounded(100, 2, 2 * bound_multiplier).is_err());

            assert_eq!(lcs_simd.scan_left_bounded(125, 5, 0 * bound_multiplier), Ok(121));
            assert_eq!(lcs_simd.scan_left_bounded(100, 5, 1 * bound_multiplier), Ok(89));
            assert_eq!(lcs_simd.scan_left_bounded(100, 4, 1 * bound_multiplier), Ok(89));
            assert_eq!(lcs_simd.scan_left_bounded(100, 3, 2 * bound_multiplier), Ok(57));
            assert_eq!(lcs_simd.scan_left_bounded(100, 2, 3 * bound_multiplier), Ok(25));
        }};
    }

    #[test]
    fn scan_left_bounded() {
        scan_left_bounded_body!(LcsSimd8x32, 31);
        scan_left_bounded_body!(LcsSimd16x32, 1<<8);
        scan_left_bounded_body!(LcsSimd32x16, 1<<16);
        scan_left_bounded_body!(LcsSimd, 31);
        scan_left_bounded_body!(LcsSimd, 1<<8);
        scan_left_bounded_body!(LcsSimd, 1<<16);
    }


    macro_rules! scan_right_bounded_body {
        ($structure:ty, $max_k:expr) => {{
            let items: &[u8] = &[
                5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
                5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 5, 5, 5, 5, 5, 5, // 25: 4

                5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
                5, 5, 5, 5, 5, 5, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5, // 57: 3

                5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
                5, 5, 5, 5, 5, 5, 5, 5, 5, 2, 5, 5, 5, 5, 5, 5, // 89: 2

                5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
                5, 5, 5, 5, 5, 5, 5, 5, 5, 1, 5, 5, 5, 5, 5, 5, // 121: 1
            ];

            let lcs_simd = <$structure>::from_iterator(items.iter().cloned(), items.len(), $max_k);
            let max_bound = lcs_simd.word_count();
            assert_eq!(
                lcs_simd.scan_right_bounded(20, 5, max_bound),
                Ok(25),
                "Result in the same SIMD word failed."
            );
            assert_eq!(
                lcs_simd.scan_right_bounded(30, 5, max_bound),
                Ok(57),
                "Result in the next SIMD word failed."
            );
            assert_eq!(
                lcs_simd.scan_right_bounded(30, 4, max_bound),
                Ok(57),
                "Smaller than 4 failed."
            );
            assert_eq!(
                lcs_simd.scan_right_bounded(30, 3, max_bound),
                Ok(89),
                "Smaller than 3 failed."
            );
            assert_eq!(
                lcs_simd.scan_right_bounded(30, 2, max_bound),
                Ok(121),
                "Smaller than 2 failed."
            );
            assert!(
                lcs_simd.scan_right_bounded(30, 1, max_bound).is_err(),
                "Scan to the end failed."
            );

            let bound_multiplier = 32 / lcs_simd.lanes();
            assert!(lcs_simd.scan_right_bounded(30, 5, 0 * bound_multiplier).is_err());
            assert!(lcs_simd.scan_right_bounded(30, 4, 0 * bound_multiplier).is_err());
            assert!(lcs_simd.scan_right_bounded(30, 3, 1 * bound_multiplier).is_err());
            assert!(lcs_simd.scan_right_bounded(30, 2, 2 * bound_multiplier).is_err());

            assert_eq!(lcs_simd.scan_right_bounded(20, 5, 0 * bound_multiplier), Ok(25));
            assert_eq!(lcs_simd.scan_right_bounded(30, 5, 1 * bound_multiplier), Ok(57));
            assert_eq!(lcs_simd.scan_right_bounded(30, 4, 1 * bound_multiplier), Ok(57));
            assert_eq!(lcs_simd.scan_right_bounded(30, 3, 2 * bound_multiplier), Ok(89));
            assert_eq!(lcs_simd.scan_right_bounded(30, 2, 3 * bound_multiplier), Ok(121));
        }};
    }

    #[test]
    fn scan_right_bounded() {
        scan_right_bounded_body!(LcsSimd8x32, 31);
        scan_right_bounded_body!(LcsSimd16x32, 1<<8);
        scan_right_bounded_body!(LcsSimd32x16, 1<<16);
        scan_right_bounded_body!(LcsSimd, 31);
        scan_right_bounded_body!(LcsSimd, 1<<8);
        scan_right_bounded_body!(LcsSimd, 1<<16);
    }
}

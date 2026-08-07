use std::io::{self, Read, Write};

const LOG_B: usize = 10;
const B: usize = 1 << LOG_B;
const PREFIX_WORDS: usize = 2;
const WORDS_PER_BLOCK: usize = PREFIX_WORDS + 2 * (B / 64);

#[derive(Clone, Eq, PartialEq, Debug)]
pub struct Base4RankVector {
    n: usize,
    n_bits: usize,
    bits: Vec<u64>,
}

impl Base4RankVector {
    #[inline(always)]
    unsafe fn write_prefix(bits: &mut [u64], word: usize, psums: &[u32; 4]) {
        let p = bits.as_mut_ptr().add(word) as *mut u32;
        *p.add(0) = psums[0];
        *p.add(1) = psums[1];
        *p.add(2) = psums[2];
        *p.add(3) = psums[3];
    }

    #[inline(always)]
    unsafe fn read_prefix(bits: &[u64], word: usize, sym: usize) -> u32 {
        *((bits.as_ptr().add(word) as *const u32).add(sym))
    }

    pub fn from_symbols(seq: &[u8]) -> Self {
        let n = seq.len();
        let nblocks = (n + B - 1) / B;
        let n_bits = nblocks * (2 * B + 128);

        let mut bits = vec![0u64; n_bits / 64];
        let mut psums = [0u32; 4];

        let mut bi = 0usize;
        let mut i = 0usize;

        while i < n {
            if i % B == 0 {
                unsafe {
                    Self::write_prefix(&mut bits, bi, &psums);
                }
                bi += PREFIX_WORDS;
            }

            let mut upper_w = 0u64;
            let mut lower_w = 0u64;
            let mut j = 0usize;

            while j < 64 && i + j < n {
                let sym = seq[i + j] as usize;

                debug_assert!(sym < 4);

                psums[sym] += 1;

                upper_w |= (((sym & 0x2) != 0) as u64) << j;
                lower_w |= (((sym & 0x1) != 0) as u64) << j;

                j += 1;
            }

            bits[bi] = upper_w;
            bi += 1;

            bits[bi] = lower_w;
            bi += 1;

            i += j;
        }

        //
        // Final block only contains prefix sums when n is a multiple of B.
        //
        if n % B == 0 && n > 0 {
            dbg!(n, nblocks);
            
            let word = (nblocks - 1) * WORDS_PER_BLOCK;
            unsafe {
                Self::write_prefix(&mut bits, word, &psums);
            }
        }

        Self { n, n_bits, bits }
    }

    #[inline(always)]
    pub fn len(&self) -> usize {
        self.n
    }

    #[inline(always)]
    pub fn size_in_bytes(&self) -> usize {
        std::mem::size_of::<Self>() + self.bits.len() * std::mem::size_of::<u64>()
    }

    /// Rank of `sym` in the half-open interval [0,pos).
    /// Returns (rank, contains_at_pos).
    #[inline(always)]
    pub fn rank_with_contains(&self, pos: usize, sym: u8) -> (usize, bool) {
        debug_assert!(sym < 4);
        debug_assert!(pos <= self.n);

        let upper_set = (sym & 0x2) != 0;
        let lower_set = (sym & 0x1) != 0;

        let blockstart = (pos >> LOG_B) * WORDS_PER_BLOCK;

        let pre_block_rank =
            unsafe { Self::read_prefix(&self.bits, blockstart, sym as usize) } as usize;

        let blockwords = &self.bits[blockstart + PREFIX_WORDS..];

        let blocki = ((pos & (B - 1)) >> 6) << 1;

        let mut whole_word_rank = 0usize;

        unsafe {
            let mut p = blockwords.as_ptr();
            let end = p.add(blocki);

            while p < end {
                let mut upper_w = *p;
                let mut lower_w = *p.add(1);

                p = p.add(2);

                upper_w = if upper_set { upper_w } else { !upper_w };
                lower_w = if lower_set { lower_w } else { !lower_w };

                whole_word_rank += (upper_w & lower_w).count_ones() as usize;
            }

            let mut leftover_rank = 0usize;

            let contains;
            if (pos & 63) != 0 {
                let mut upper_w = *p;
                let mut lower_w = *p.add(1);

                upper_w = if upper_set { upper_w } else { !upper_w };
                lower_w = if lower_set { lower_w } else { !lower_w };

                let shift = 64 - (pos & 63);

                leftover_rank = ((upper_w & lower_w) << shift).count_ones() as usize;

                // check if the symbol at pos equals sym
                contains = pos < self.n && ((((upper_w & lower_w) >> (pos & 63)) & 1) != 0) as bool;
            } else {
                // we should inspect the first element of blocki
                let mut upper_w = *p;
                let mut lower_w = *p.add(1);
                upper_w = if upper_set { upper_w } else { !upper_w };
                lower_w = if lower_set { lower_w } else { !lower_w };

                let upper_bit = ((upper_w) & 1) != 0;
                let lower_bit = ((lower_w) & 1) != 0;

                contains = pos < self.n && (upper_bit & lower_bit);
            }

            let rank = pre_block_rank + whole_word_rank + leftover_rank;

            (rank, contains)
        }
    }

    #[inline(always)]
    pub fn window256_starting_from_pos_with_sym_word_aligned(
        &self,
        pos: usize,
        sym: u8,
    ) -> [u64; 4] {
        let mut words = [0u64; 4];

        if pos >= self.n {
            return words;
        }

        debug_assert!(sym < 4);

        let upper_set = (sym & 0x2) != 0;
        let lower_set = (sym & 0x1) != 0;

        let global_word = pos >> 6;

        let words_per_block = B / 64;

        let mut block = global_word / words_per_block;
        let mut word_in_block = global_word % words_per_block;

        unsafe {
            for i in 0..4 {
                let word_ptr = self
                    .bits
                    .as_ptr()
                    .add(block * WORDS_PER_BLOCK + PREFIX_WORDS + word_in_block * 2);

                let mut upper = *word_ptr;
                let mut lower = *word_ptr.add(1);

                upper = if upper_set { upper } else { !upper };
                lower = if lower_set { lower } else { !lower };

                let value = upper & lower;

                words[i] = value;

                word_in_block += 1;

                if word_in_block == words_per_block {
                    block += 1;
                    word_in_block = 0;
                }
            }
        }

        words
    }

    #[inline(always)]
    pub fn get_words_in_range_for_sym(&self, range: std::ops::Range<usize>, sym: u8) -> Vec<u64> {
        // length of vector with first position offset to be aligned to a 64 bit word as in the data layout
        let len = (range.start & 63) + range.end.min(self.n - 1) - range.start;
        let num_words = (len + 64 - 1) / 64;
        let mut words: Vec<u64> = vec![0; num_words];
        let pos = range.start;
        if pos >= self.n {
            return words;
        }

        debug_assert!(sym < 4);

        let upper_set = (sym & 0x2) != 0;
        let lower_set = (sym & 0x1) != 0;

        let global_word = pos >> 6;

        let words_per_block = B / 64;

        let mut block = global_word / words_per_block;
        let mut word_in_block = global_word % words_per_block;

        unsafe {
            for i in 0..num_words {
                let word_ptr = self
                    .bits
                    .as_ptr()
                    .add(block * WORDS_PER_BLOCK + PREFIX_WORDS + word_in_block * 2);

                let mut upper = *word_ptr;
                let mut lower = *word_ptr.add(1);

                upper = if upper_set { upper } else { !upper };
                lower = if lower_set { lower } else { !lower };

                let value = upper & lower;

                words[i] = value;

                word_in_block += 1;

                if word_in_block == words_per_block {
                    block += 1;
                    word_in_block = 0;
                }
            }
        }

        words
    }

    /// concat[pos] == sym
    #[inline(always)]
    pub fn contains(&self, pos: usize, sym: u8) -> bool {
        debug_assert!(sym < 4);

        if pos >= self.n {
            return false;
        }

        self.sym_at_pos(pos) == sym
    }

    /// returns the sym at concat[pos]
    #[inline(always)]
    pub fn sym_at_pos(&self, pos: usize) -> u8 {
        debug_assert!(pos < self.n);

        let blockstart = (pos >> LOG_B) * WORDS_PER_BLOCK;
        let blockwords = &self.bits[blockstart + PREFIX_WORDS..];
        let blocki = ((pos & (B - 1)) >> 6) << 1;

        unsafe {
            let p = blockwords.as_ptr().add(blocki);

            let word = (pos & 63) as u64;

            let upper_bit = ((*p >> word) & 1) != 0;
            let lower_bit = ((*p.add(1) >> word) & 1) != 0;

            // value is the symbol written at pos
            let value = ((upper_bit as u8) << 1) | (lower_bit as u8);

            value
        }
    }

    pub fn serialize<W: Write>(&self, mut out: W) -> io::Result<usize> {
        debug_assert_eq!(self.bits.len(), (self.n_bits + 63) / 64);

        let mut written = 0;

        out.write_all(&(self.n as u64).to_le_bytes())?;
        written += 8;

        out.write_all(&(self.n_bits as u64).to_le_bytes())?;
        written += 8;

        out.write_all(&(self.bits.len() as u64).to_le_bytes())?;
        written += 8;

        unsafe {
            let bytes = std::slice::from_raw_parts(
                self.bits.as_ptr() as *const u8,
                self.bits.len() * std::mem::size_of::<u64>(),
            );
            out.write_all(bytes)?;
        }

        written += self.bits.len() * std::mem::size_of::<u64>();

        Ok(written)
    }

    pub fn load<R: Read>(mut input: R) -> io::Result<Self> {
        fn read_u64<R: Read>(r: &mut R) -> io::Result<u64> {
            let mut buf = [0u8; 8];
            r.read_exact(&mut buf)?;
            Ok(u64::from_le_bytes(buf))
        }

        let n = read_u64(&mut input)? as usize;
        let n_bits = read_u64(&mut input)? as usize;
        let bits_len = read_u64(&mut input)? as usize;

        debug_assert_eq!(bits_len, (n_bits + 63) / 64);

        let mut bits = vec![0u64; bits_len];

        unsafe {
            let bytes = std::slice::from_raw_parts_mut(
                bits.as_mut_ptr() as *mut u8,
                bits.len() * std::mem::size_of::<u64>(),
            );
            input.read_exact(bytes)?;
        }

        Ok(Self { n, n_bits, bits })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn build_test_sequence(len: usize) -> Vec<u8> {
        (0..len).map(|i| (i % 4) as u8).collect()
    }

    #[test]
    fn sym_at_pos_roundtrip() {
        let seq = build_test_sequence(5000);
        let rv = Base4RankVector::from_symbols(&seq);

        assert_eq!(rv.len(), seq.len());

        for i in 0..seq.len() {
            assert_eq!(rv.sym_at_pos(i), seq[i]);
        }
    }

    #[test]
    fn contains_matches_sequence() {
        let seq = build_test_sequence(3000);
        let rv = Base4RankVector::from_symbols(&seq);

        for i in 0..seq.len() {
            for sym in 0..4 {
                assert_eq!(rv.contains(i, sym), seq[i] == sym);
            }
        }

        for sym in 0..4 {
            assert!(!rv.contains(seq.len(), sym));
            assert!(!rv.contains(seq.len() + 100, sym));
        }
    }

    #[test]
    fn rank_matches_naive() {
        let seq = build_test_sequence(5000);
        let rv = Base4RankVector::from_symbols(&seq);

        for pos in 0..=seq.len() {
            for sym in 0..4 {
                let expected_rank = seq[..pos].iter().filter(|&&x| x == sym).count();

                let expected_contains = pos < seq.len() && seq[pos] == sym;

                let (rank, contains) = rv.rank_with_contains(pos, sym);

                assert_eq!(rank, expected_rank);
                assert_eq!(contains, expected_contains);
            }
        }
    }

    #[test]
    fn rank_across_word_boundaries() {
        let seq = build_test_sequence(300);
        let rv = Base4RankVector::from_symbols(&seq);

        let positions = [
            0, 1, 63, 64, 65, 127, 128, 129, 191, 192, 193, 255, 256, 257, 299, 300,
        ];

        for &pos in &positions {
            for sym in 0..4 {
                let expected_rank = seq[..pos].iter().filter(|&&x| x == sym).count();

                let expected_contains = pos < seq.len() && seq[pos] == sym;

                let (rank, contains) = rv.rank_with_contains(pos, sym);

                assert_eq!(rank, expected_rank);
                assert_eq!(contains, expected_contains);
            }
        }
    }

    #[test]
    fn rank_across_block_boundaries() {
        let seq = build_test_sequence(3000);
        let rv = Base4RankVector::from_symbols(&seq);

        let positions = [1023, 1024, 1025, 2047, 2048, 2049, 2999, 3000];

        for &pos in &positions {
            for sym in 0..4 {
                let expected_rank = seq[..pos].iter().filter(|&&x| x == sym).count();

                let expected_contains = pos < seq.len() && seq[pos] == sym;

                let (rank, contains) = rv.rank_with_contains(pos, sym);

                assert_eq!(rank, expected_rank);
                assert_eq!(contains, expected_contains);
            }
        }
    }

    #[test]
    fn serialize_roundtrip() {
        let seq = build_test_sequence(4096);

        let rv = Base4RankVector::from_symbols(&seq);

        let mut bytes = Vec::new();
        rv.serialize(&mut bytes).unwrap();

        let loaded = Base4RankVector::load(bytes.as_slice()).unwrap();

        assert_eq!(rv, loaded);

        for i in 0..seq.len() {
            assert_eq!(loaded.sym_at_pos(i), seq[i]);
        }
    }

    #[test]
    fn window256_matches_naive() {
        let seq = build_test_sequence(5000);
        let rv = Base4RankVector::from_symbols(&seq);

        let starts = [0, 64, 128, 512, 960, 1024, 2048, 3000, 490];

        for &start in &starts {
            for sym in 0..4 {
                let words = rv.window256_starting_from_pos_with_sym_word_aligned(start, sym);

                for w in 0..4 {
                    for bit in 0..64 {
                        let pos = start.saturating_sub(start & 63) + w * 64 + bit;

                        let expected = pos < seq.len() && seq[pos] == sym;

                        let actual = ((words[w] >> bit) & 1) != 0;

                        assert_eq!(
                            actual, expected,
                            "start={}, word={}, bit={}, sym={}",
                            start, w, bit, sym
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn get_words_matches_naive() {
        let seq = build_test_sequence(3000);
        let rv = Base4RankVector::from_symbols(&seq);

        let ranges = [0..64, 5..120, 63..140, 64..192, 127..320, 1023..1200];

        for range in ranges {
            for sym in 0..4 {
                let words = rv.get_words_in_range_for_sym(range.clone(), sym);

                let offset = range.start & 63;

                for (wi, word) in words.iter().enumerate() {
                    for bit in 0..64 {
                        let absolute = range.start + wi * 64 + bit - offset;

                        let expected = seq[absolute] == sym;
                        let actual = ((*word >> bit) & 1) != 0;

                        assert_eq!(
                            actual, expected,
                            "range={:?}, sym={}, absolute={}",
                            range, sym, absolute
                        );
                    }
                }
            }
        }
    }
}

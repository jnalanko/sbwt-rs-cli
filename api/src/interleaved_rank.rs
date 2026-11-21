//! The interleaved rank: bit vectors whose rank counts are interleaved with the bits (as in
//! SDSL's `bit_vector_il`), so that a rank query reads one cache line.
//!
//! The bits are stored in 64-byte, cache-line-aligned blocks of a 96-bit header (the number
//! of ones before the block and cumulative counts per word) and [BITS_PER_LINE] = 416 bits of
//! the vector. A rank query reads only the block and does one popcount. Space overhead
//! 96/416 = 23.1%, about that of rank9 (25%), but with one cache line per query instead of two.
//!
//! [InterleavedRank] is one such bit vector, and [SubsetInterleavedMatrix] is a [SubsetSeq] with one
//! of them per character. It has its own file format; an index of [SubsetMatrix] is converted
//! once from its rows (`transform-index --type interleaved_matrix`, which builds it with
//! [SubsetSeq::new_from_bit_vectors] from [SubsetSeq::into_bitvectors] of the matrix).

use std::ops::Range;

use crate::subsetseq::{SubsetMatrix, SubsetSeq};
use bitvec::order::Lsb0;

/// Number of bits of the vector stored in one block.
pub const BITS_PER_LINE: usize = 416;

/// [InterleavedRank::select] keeps the block of every `SELECT_SAMPLE`-th one, so that it searches
/// only the blocks between two samples (computed when the vector is built or loaded, not stored).
const SELECT_SAMPLE: usize = 4096;

// Bits of a block before the data: the header.
const HEADER_BITS: usize = 96;

// One cache line: words[0] and the low half of words[1] are the header, the high half of
// words[1] and words[2..8] hold 416 bits of the vector.
#[repr(C, align(64))]
#[derive(Clone, Copy, Default, Debug, PartialEq, Eq)]
struct Line {
    words: [u64; 8],
}

// SAFETY: Line is 8 u64s with alignment 64 and size 64: no padding, and every bit pattern is valid.
unsafe impl bytemuck::Zeroable for Line {}
unsafe impl bytemuck::Pod for Line {}

// Header of a Line. words[0] bits 0..40: the number of ones before the block (so rows can
// have up to 2^40 bits). Then the cumulative number of ones in the data before each data
// word, each field just wide enough for its maximum (data before word w: 32 + 64 (w - 2)
// bits). No field crosses a word boundary. Indexed by the data word 1..8; the data before
// word 1 is empty, so its count is 0 and not stored.
const ABSOLUTE_BITS: u32 = 40;
const COUNT_WORD: [usize; 8] = [0, 0, 0, 0, 0, 1, 1, 1];
const COUNT_SHIFT: [u32; 8] = [0, 0, 40, 46, 53, 0, 8, 17];
const COUNT_MASK: [u64; 8] = [0, 0, 0x3F, 0x7F, 0xFF, 0xFF, 0x1FF, 0x1FF];
// Bits of each word that hold data: the low half of words[1] is header.
const DATA_MASK: [u64; 8] = [0, 0xFFFF_FFFF_0000_0000, !0, !0, !0, !0, !0, !0];

// A word with the lowest n bits set, for n <= 64.
fn low_bits(n: usize) -> u64 {
    if n == 64 {
        u64::MAX
    } else {
        (1_u64 << n) - 1
    }
}

// Returns n <= 64 bits of the bit vector `words` of length `len` starting at bit `start`,
// least significant bit first. Bits at or past `len` are zero.
fn get_bits(words: &[u64], len: usize, start: usize, n: usize) -> u64 {
    if start >= len {
        return 0;
    }
    let (w, shift) = (start / 64, start % 64);
    let mut bits = words[w] >> shift;
    if shift > 0 && w + 1 < words.len() {
        bits |= words[w + 1] << (64 - shift);
    }
    bits & low_bits(std::cmp::min(n, len - start))
}

impl Line {
    // Number of ones before the block.
    #[inline]
    fn absolute(&self) -> usize {
        (self.words[0] & low_bits(ABSOLUTE_BITS as usize)) as usize
    }

    // Number of ones in the first `offset` data bits of the block plus the ones before it,
    // with offset < BITS_PER_LINE.
    #[inline]
    fn rank(&self, offset: usize) -> usize {
        let position = offset + HEADER_BITS; // Bit position within the block
        let word = position / 64; // 1..8
        let bit = position % 64;
        let before_word = (self.words[COUNT_WORD[word]] >> COUNT_SHIFT[word]) & COUNT_MASK[word];
        let within_word = (self.words[word] & ((1_u64 << bit) - 1) & DATA_MASK[word]).count_ones() as usize;
        self.absolute() + before_word as usize + within_word
    }

    // The data bits of word w (1..8), shifted so that bit 0 is the first data bit of the word.
    #[inline]
    fn data_word(&self, w: usize) -> u64 {
        if w == 1 { self.words[1] >> 32 } else { self.words[w] }
    }

    // Offset in the data of the first bit of data word w (1..8).
    #[inline]
    fn data_start(w: usize) -> usize {
        if w == 1 { 0 } else { 64 * w - HEADER_BITS }
    }
}

/// A bit vector in interleaved blocks: rank reads one cache line and does one popcount, with no
/// other table. See the [module documentation](self).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct InterleavedRank {
    lines: Vec<Line>,
    len: usize,
    ones: usize,
    // The block of the one of rank j * SELECT_SAMPLE, for every j.
    select_samples: Vec<u32>,
}

// The block of the one of rank j * SELECT_SAMPLE for every j (see InterleavedRank::select).
fn select_samples(lines: &[Line], ones: usize) -> Vec<u32> {
    assert!(lines.len() <= u32::MAX as usize, "too many blocks for the select samples");
    let mut samples = Vec::with_capacity(ones.div_ceil(SELECT_SAMPLE));
    let mut l = 0;
    for r in (0..ones).step_by(SELECT_SAMPLE) {
        while l + 1 < lines.len() && lines[l + 1].absolute() <= r {
            l += 1;
        }
        samples.push(l as u32);
    }
    samples
}

impl InterleavedRank {
    /// Builds the vector from `len` bits packed into `words`, least significant bit first:
    /// bit `i` is bit `i % 64` of `words[i / 64]`. Bits of `words` past `len` are ignored.
    pub fn new(words: &[u64], len: usize) -> Self {
        assert!(words.len() * 64 >= len, "not enough words for {} bits", len);
        assert!((len as u64) < 1 << ABSOLUTE_BITS, "too many bits for InterleavedRank");
        // One more block than needed for the bits, so that rank(len) also has a block.
        let mut lines = vec![Line::default(); len / BITS_PER_LINE + 1];
        let mut ones = 0_usize;
        for (l, line) in lines.iter_mut().enumerate() {
            let start = l * BITS_PER_LINE;
            // Data: 32 bits in the high half of words[1], then 64 bits in each of words[2..8].
            line.words[1] = get_bits(words, len, start, 32) << 32;
            for w in 2..8 {
                line.words[w] = get_bits(words, len, start + Line::data_start(w), 64);
            }
            // Header: absolute count, then cumulative counts before data words 2..8.
            line.words[0] |= ones as u64;
            let mut before_word = 0_u64;
            for w in 1..8 {
                debug_assert!(before_word <= COUNT_MASK[w]); // 0 for w == 1
                line.words[COUNT_WORD[w]] |= before_word << COUNT_SHIFT[w];
                before_word += (line.words[w] & DATA_MASK[w]).count_ones() as u64;
            }
            ones += before_word as usize;
        }
        let select_samples = select_samples(&lines, ones);
        InterleavedRank { lines, len, ones, select_samples }
    }

    /// Number of bits in the vector.
    pub fn len(&self) -> usize {
        self.len
    }

    /// Number of ones in the vector.
    #[cfg(test)]
    pub fn count_ones(&self) -> usize {
        self.ones
    }

    /// Number of ones in positions `[0..i)`. Requires `i <= len`.
    #[inline]
    pub fn rank(&self, i: usize) -> usize {
        self.lines[i / BITS_PER_LINE].rank(i % BITS_PER_LINE)
    }

    /// Returns `(self.rank(i), self.rank(j))`. Requires `i, j <= len`.
    #[inline]
    pub fn rank_pair(&self, i: usize, j: usize) -> (usize, usize) {
        (self.rank(i), self.rank(j))
    }

    /// Like [InterleavedRank::rank_pair], but if `i` and `j` are in the same block, reads the
    /// block only once. Branches on whether they are in the same block.
    #[inline]
    pub fn rank_pair_shared(&self, i: usize, j: usize) -> (usize, usize) {
        let l = i / BITS_PER_LINE;
        if l == j / BITS_PER_LINE {
            let line = &self.lines[l];
            (line.rank(i % BITS_PER_LINE), line.rank(j % BITS_PER_LINE))
        } else {
            self.rank_pair(i, j)
        }
    }

    /// Hints the CPU to start loading the block that [InterleavedRank::rank] reads for `i`.
    #[inline]
    pub fn prefetch(&self, i: usize) {
        // Plain address arithmetic, no bounds check: a prefetch of any address is harmless.
        crate::util::prefetch_read(self.lines.as_ptr().wrapping_add(i / BITS_PER_LINE));
    }

    /// The bit at position `i < len`.
    #[inline]
    pub fn get(&self, i: usize) -> bool {
        let position = i % BITS_PER_LINE + HEADER_BITS;
        (self.lines[i / BITS_PER_LINE].words[position / 64] >> (position % 64)) & 1 == 1
    }

    /// Position of the one of rank `r` (the (r+1)-th one), or None if there are at most `r` ones:
    /// a binary search over the counts in the headers of the blocks between the two select
    /// samples around `r`, then a scan of the block.
    pub fn select(&self, r: usize) -> Option<usize> {
        if r >= self.ones {
            return None;
        }
        // The last block with at most r ones before it, between the samples around r.
        let j = r / SELECT_SAMPLE;
        let lo = self.select_samples[j] as usize;
        let hi = self.select_samples.get(j + 1).map_or(self.lines.len(), |&s| s as usize + 1);
        let l = lo + self.lines[lo..hi].partition_point(|line| line.absolute() <= r) - 1;
        let line = &self.lines[l];
        let mut left = r - line.absolute();
        for w in 1..8 {
            let mut bits = line.data_word(w);
            let n = bits.count_ones() as usize;
            if left < n {
                for _ in 0..left {
                    bits &= bits - 1;
                }
                return Some(l * BITS_PER_LINE + Line::data_start(w) + bits.trailing_zeros() as usize);
            }
            left -= n;
        }
        unreachable!("the block headers do not match the data")
    }

    /// Position of the first one at or after `i`, or None if there is none: a forward scan of
    /// the data words, skipping the block headers.
    pub fn next_one(&self, i: usize) -> Option<usize> {
        if i >= self.len {
            return None;
        }
        let (mut l, offset) = (i / BITS_PER_LINE, i % BITS_PER_LINE);
        // The data word that holds offset, and the bits of it before offset.
        let mut w = if offset < 32 { 1 } else { (offset + HEADER_BITS) / 64 };
        let mut skip = offset - Line::data_start(w);
        while l < self.lines.len() {
            let line = &self.lines[l];
            while w < 8 {
                let mut bits = line.data_word(w);
                if skip > 0 {
                    bits &= !low_bits(skip);
                    skip = 0;
                }
                if bits != 0 {
                    let position = l * BITS_PER_LINE + Line::data_start(w) + bits.trailing_zeros() as usize;
                    return (position < self.len).then_some(position);
                }
                w += 1;
            }
            (l, w) = (l + 1, 1);
        }
        None
    }

    /// Calls `f` on the position of every one in `range`, in increasing order.
    pub fn for_each_one<F: FnMut(usize)>(&self, range: Range<usize>, mut f: F) {
        if range.is_empty() {
            return;
        }
        for l in range.start / BITS_PER_LINE..=(range.end - 1) / BITS_PER_LINE {
            let line = &self.lines[l];
            let line_start = l * BITS_PER_LINE;
            for w in 1..8 {
                let start = line_start + Line::data_start(w); // Position of the word's first bit
                let n = if w == 1 { 32 } else { 64 };
                if start + n <= range.start || start >= range.end {
                    continue;
                }
                let mut bits = line.data_word(w);
                if range.start > start {
                    bits &= !low_bits(range.start - start);
                }
                if range.end < start + n {
                    bits &= low_bits(range.end - start);
                }
                while bits != 0 {
                    f(start + bits.trailing_zeros() as usize);
                    bits &= bits - 1;
                }
            }
        }
    }

    /// Size of the structure in memory in bytes: the blocks and the select samples (the samples
    /// are not stored in the file).
    pub fn size_in_bytes(&self) -> usize {
        self.lines.len() * std::mem::size_of::<Line>() + self.select_samples.len() * 4
    }

    fn serialize<W: std::io::Write>(&self, out: &mut W) -> std::io::Result<usize> {
        out.write_all(&(self.len as u64).to_le_bytes())?;
        out.write_all(&(self.ones as u64).to_le_bytes())?;
        out.write_all(&(self.lines.len() as u64).to_le_bytes())?;
        let blocks: &[u8] = bytemuck::cast_slice(self.lines.as_slice());
        out.write_all(blocks)?;
        Ok(24 + blocks.len())
    }

    fn load<R: std::io::Read>(input: &mut R) -> std::io::Result<Self> {
        let len = read_u64(input)? as usize;
        let ones = read_u64(input)? as usize;
        let n_lines = read_u64(input)? as usize;
        if n_lines != len / BITS_PER_LINE + 1 {
            return Err(std::io::Error::new(std::io::ErrorKind::InvalidData, "interleaved rank: wrong number of blocks"));
        }
        // Read straight into the aligned blocks, without a second copy.
        let mut lines = vec![Line::default(); n_lines];
        input.read_exact(bytemuck::cast_slice_mut(lines.as_mut_slice()))?;
        let select_samples = select_samples(&lines, ones);
        Ok(InterleavedRank { lines, len, ones, select_samples })
    }
}

fn read_u64<R: std::io::Read>(input: &mut R) -> std::io::Result<u64> {
    let mut buf = [0u8; 8];
    input.read_exact(&mut buf)?;
    Ok(u64::from_le_bytes(buf))
}

//-----------------------------------------------------------------------------

// Version of the file format of SubsetInterleavedMatrix.
const FORMAT: u64 = 1;

/// A [SubsetSeq] with one [InterleavedRank] bit vector per character: the i-th bit of the c-th
/// vector is 1 if and only if the i-th subset contains c, as in [SubsetMatrix], but each rank
/// query reads one cache line. Rank and select support are always there. See the
/// [module documentation](self).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SubsetInterleavedMatrix {
    rows: Vec<InterleavedRank>,
    len: usize,
}

impl SubsetInterleavedMatrix {
    /// The alphabet size (the number of rows).
    pub fn sigma(&self) -> usize {
        self.rows.len()
    }

    /// Size of the subset sequence in memory, in bytes.
    pub fn size_in_bytes(&self) -> usize {
        self.rows.iter().map(|row| row.size_in_bytes()).sum()
    }
}

impl SubsetSeq for SubsetInterleavedMatrix {
    fn new(subset_seq: Vec<Vec<u8>>, sigma: usize) -> Self {
        Self::new_from_bit_vectors(SubsetMatrix::new(subset_seq, sigma).into_bitvectors())
    }

    fn new_from_bit_vectors(vecs: Vec<bitvec::vec::BitVec<u64, Lsb0>>) -> Self {
        let len = vecs.first().map_or(0, |v| v.len());
        let rows = vecs.iter().map(|v| InterleavedRank::new(v.as_raw_slice(), v.len())).collect();
        SubsetInterleavedMatrix { rows, len }
    }

    fn len(&self) -> usize {
        self.len
    }

    // Rank and select support are always there.
    fn build_rank(&mut self) {}

    fn build_select(&mut self) {}

    fn has_select_support(&self) -> bool {
        true
    }

    fn has_rank_support(&self) -> bool {
        true
    }

    #[inline]
    fn rank(&self, c: u8, i: usize) -> usize {
        self.rows[c as usize].rank(i)
    }

    #[inline]
    fn prefetch_rank(&self, c: u8, i: usize) {
        self.rows[c as usize].prefetch(i)
    }

    #[inline]
    fn rank_pair(&self, c: u8, i: usize, j: usize) -> (usize, usize) {
        self.rows[c as usize].rank_pair(i, j)
    }

    #[inline]
    fn rank_pair_shared(&self, c: u8, i: usize, j: usize) -> (usize, usize) {
        self.rows[c as usize].rank_pair_shared(i, j)
    }

    #[inline]
    fn prefetch_rank_pair(&self, c: u8, i: usize, j: usize) {
        let row = &self.rows[c as usize];
        row.prefetch(i);
        row.prefetch(j);
    }

    fn select(&self, c: u8, i: usize) -> Option<usize> {
        self.rows[c as usize].select(i)
    }

    fn append_set_to_buf(&self, i: usize, buf: &mut Vec<u8>) {
        for (c, row) in self.rows.iter().enumerate() {
            if row.get(i) {
                buf.push(c as u8);
            }
        }
    }

    fn subset_size(&self, i: usize) -> usize {
        self.rows.iter().filter(|row| row.get(i)).count()
    }

    fn set_contains(&self, set_idx: usize, character: u8) -> bool {
        self.rows[character as usize].get(set_idx)
    }

    fn next_set_with_char(&self, set_idx: usize, c: u8) -> Option<usize> {
        self.rows[c as usize].next_one(set_idx)
    }

    fn serialize<W: std::io::Write>(&self, out: &mut W) -> std::io::Result<usize> {
        out.write_all(&FORMAT.to_le_bytes())?;
        out.write_all(&(self.len as u64).to_le_bytes())?;
        out.write_all(&(self.rows.len() as u64).to_le_bytes())?;
        let mut written = 24;
        for row in &self.rows {
            written += row.serialize(out)?;
        }
        Ok(written)
    }

    fn load<R: std::io::Read>(input: &mut R) -> std::io::Result<Self> {
        let format = read_u64(input)?;
        if format != FORMAT {
            return Err(std::io::Error::new(std::io::ErrorKind::InvalidData,
                "interleaved matrix of another format version: convert the index again (transform-index --type interleaved_matrix)"));
        }
        let len = read_u64(input)? as usize;
        let n_rows = read_u64(input)? as usize;
        let rows = (0..n_rows).map(|_| InterleavedRank::load(input)).collect::<std::io::Result<Vec<_>>>()?;
        if rows.iter().any(|row| row.len() != len) {
            return Err(std::io::Error::new(std::io::ErrorKind::InvalidData, "interleaved matrix: rows of different lengths"));
        }
        Ok(SubsetInterleavedMatrix { rows, len })
    }

    fn call_on_char_occurrences<F: FnMut(usize)>(&self, range: Range<usize>, c: u8, callback: F) {
        self.rows[c as usize].for_each_one(range, callback)
    }

    fn push_labels_forward(&self, labels: &[u8], sbwt_input_range: Range<usize>, output_ranges: Vec<&mut [u8]>) {
        assert_eq!(labels.len(), sbwt_input_range.len());
        assert_eq!(output_ranges.len(), self.rows.len());
        for (row, output_slice) in self.rows.iter().zip(output_ranges) {
            let mut output_offset = 0_usize;
            row.for_each_one(sbwt_input_range.clone(), |i| {
                output_slice[output_offset] = labels[i - sbwt_input_range.start];
                output_offset += 1;
            });
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{Rng, SeedableRng};

    // Rank by counting bit by bit.
    fn naive_ranks(words: &[u64], len: usize) -> Vec<usize> {
        let mut ranks = vec![0];
        for i in 0..len {
            let bit = (words[i / 64] >> (i % 64)) & 1;
            ranks.push(ranks[i] + bit as usize);
        }
        ranks
    }

    fn check(words: &[u64], len: usize) {
        let v = InterleavedRank::new(words, len);
        let ranks = naive_ranks(words, len);
        let bit = |i: usize| (words[i / 64] >> (i % 64)) & 1 == 1;
        let ones: Vec<usize> = (0..len).filter(|&i| bit(i)).collect();
        assert_eq!(v.len(), len);
        assert_eq!(v.count_ones(), ranks[len]);
        for i in 0..=len {
            assert_eq!(v.rank(i), ranks[i], "rank({}) with len {}", i, len);
        }
        for i in 0..len {
            assert_eq!(v.get(i), bit(i), "get({}) with len {}", i, len);
        }
        for (r, &i) in ones.iter().enumerate() {
            assert_eq!(v.select(r), Some(i), "select({}) with len {}", r, len);
        }
        assert_eq!(v.select(ones.len()), None);
        // The next one from every position, past the end included.
        for i in 0..=len + 1 {
            let expected = ones.iter().copied().find(|&o| o >= i);
            assert_eq!(v.next_one(i), expected, "next_one({}) with len {}", i, len);
        }
        // Ranges starting and ending inside and at the borders of words and blocks.
        for start in (0..=len).step_by(13) {
            for end in [start, (start + 1).min(len), (start + 31).min(len), (start + 100).min(len), (start + 900).min(len), len] {
                let mut got = vec![];
                v.for_each_one(start..end, |i| got.push(i));
                let expected: Vec<usize> = ones.iter().copied().filter(|&i| i >= start && i < end).collect();
                assert_eq!(got, expected, "range {}..{} with len {}", start, end, len);
            }
        }
        // Pairs: neighbors, within a block, across blocks, and (i, len); prefetching first
        // (it must not change anything).
        for i in (0..=len).step_by(7) {
            for j in [i, (i + 1).min(len), (i + 63).min(len), (i + 500).min(len), (i + 15000).min(len), len] {
                v.prefetch(i);
                v.prefetch(j);
                assert_eq!(v.rank_pair(i, j), (ranks[i], ranks[j]));
                assert_eq!(v.rank_pair_shared(i, j), (ranks[i], ranks[j]), "pair ({}, {})", i, j);
            }
        }
        // Serialization round trip.
        let mut buf = vec![];
        let written = v.serialize(&mut buf).unwrap();
        assert_eq!(written, buf.len());
        assert_eq!(InterleavedRank::load(&mut buf.as_slice()).unwrap(), v);
    }

    #[test]
    fn matches_naive() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(1);
        // Lengths around word and block boundaries.
        for len in [0, 1, 31, 32, 33, 63, 64, 65, 95, 96, 97, 415, 416, 417, 447, 448, 449, 831, 832, 833, 4480, 30000] {
            // Density 1.0 reaches the largest values of all header fields.
            for density in [0.0, 0.1, 0.5, 1.0] {
                let words: Vec<u64> = (0..(len + 63) / 64 + 1)
                    .map(|_| (0..64).fold(0_u64, |w, b| w | ((rng.gen::<f64>() < density) as u64) << b))
                    .collect();
                check(&words, len);
            }
        }
    }

    // select over many select samples, with long stretches of zeros between dense parts (so that
    // samples are many blocks apart), and after a serialization round trip (the samples are
    // computed again when loading).
    #[test]
    fn select_with_samples() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(3);
        let len = 2_000_000;
        let mut words = vec![0_u64; len / 64 + 1];
        for i in 0..len {
            let dense = (i / 100_000) % 3 != 1; // every third stretch of 100 000 bits is empty
            if dense && rng.gen::<f64>() < 0.3 {
                words[i / 64] |= 1 << (i % 64);
            }
        }
        let v = InterleavedRank::new(&words, len);
        let ones: Vec<usize> = (0..len).filter(|&i| (words[i / 64] >> (i % 64)) & 1 == 1).collect();
        assert!(ones.len() > 50 * SELECT_SAMPLE);
        let mut buf = vec![];
        v.serialize(&mut buf).unwrap();
        let loaded = InterleavedRank::load(&mut buf.as_slice()).unwrap();
        assert_eq!(loaded, v);
        for (r, &i) in ones.iter().enumerate() {
            assert_eq!(v.select(r), Some(i), "select({})", r);
        }
        assert_eq!(v.select(ones.len()), None);
    }

    #[test]
    fn bits_past_len_are_ignored() {
        let v = InterleavedRank::new(&[u64::MAX, u64::MAX], 70);
        assert_eq!((v.count_ones(), v.rank(70)), (70, 70));
    }

    #[test]
    fn matches_subset_matrix() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(2);
        let sets: Vec<Vec<u8>> = (0..20000).map(|_| (0..4).filter(|_| rng.gen::<bool>()).collect()).collect();
        let mut plain = SubsetMatrix::new(sets.clone(), 4);
        plain.build_rank();
        plain.build_select();
        let interleaved = SubsetInterleavedMatrix::new(sets, 4);
        let n = plain.len();
        assert_eq!(interleaved.len(), n);
        for c in 0..4 {
            for i in 0..=n {
                assert_eq!(interleaved.rank(c, i), plain.rank(c, i));
                assert_eq!(interleaved.next_set_with_char(i, c), plain.next_set_with_char(i, c), "next_set_with_char({}, {})", i, c);
            }
            for r in 0..=plain.rank(c, n) {
                assert_eq!(interleaved.select(c, r), plain.select(c, r));
            }
            let (mut a, mut b) = (vec![], vec![]);
            interleaved.call_on_char_occurrences(123..15000, c, |i| a.push(i));
            plain.call_on_char_occurrences(123..15000, c, |i| b.push(i));
            assert_eq!(a, b);
        }
        for i in 0..n {
            assert_eq!(interleaved.access(i), plain.access(i));
            assert_eq!(interleaved.subset_size(i), plain.subset_size(i));
        }
        assert_eq!(interleaved.get_C_array(), plain.get_C_array());
        // push_labels_forward on a range, against the plain matrix.
        let range = 1000..9000;
        let labels: Vec<u8> = (0..range.len()).map(|i| (i % 251) as u8).collect();
        let counts: Vec<usize> = (0..4).map(|c| plain.rank(c, range.end) - plain.rank(c, range.start)).collect();
        let mut out_a: Vec<Vec<u8>> = counts.iter().map(|&m| vec![0; m]).collect();
        let mut out_b = out_a.clone();
        interleaved.push_labels_forward(&labels, range.clone(), out_a.iter_mut().map(|v| v.as_mut_slice()).collect());
        plain.push_labels_forward(&labels, range, out_b.iter_mut().map(|v| v.as_mut_slice()).collect());
        assert_eq!(out_a, out_b);
        // Bit vectors and serialization round trip.
        assert_eq!(SubsetInterleavedMatrix::new_from_bit_vectors(plain.clone().into_bitvectors()), interleaved);
        let mut buf = vec![];
        let written = interleaved.serialize(&mut buf).unwrap();
        assert_eq!(written, buf.len());
        assert_eq!(SubsetInterleavedMatrix::load(&mut buf.as_slice()).unwrap(), interleaved);
    }

    #[test]
    fn header_fields_fit() {
        assert_eq!(std::mem::size_of::<Line>(), 64);
        assert_eq!(std::mem::align_of::<Line>(), 64);
        assert_eq!(HEADER_BITS + BITS_PER_LINE, 512);
        // Each count field fits its maximum (the data before word w), lies in the header part
        // of its word, and does not overlap the absolute count or other fields.
        let mut used = [0_u64; 2];
        used[0] = low_bits(ABSOLUTE_BITS as usize);
        for w in 2..8 {
            let max_before = 32 + 64 * (w as u64 - 2);
            assert!(COUNT_MASK[w] >= max_before, "field of word {}", w);
            let field = COUNT_MASK[w] << COUNT_SHIFT[w];
            assert_eq!(field >> COUNT_SHIFT[w], COUNT_MASK[w], "field of word {} crosses the word", w);
            let hw = COUNT_WORD[w];
            assert_eq!(used[hw] & field, 0, "field of word {} overlaps", w);
            used[hw] |= field;
        }
        assert_eq!(used[1] & DATA_MASK[1], 0, "header fields overlap the data in word 1");
    }
}

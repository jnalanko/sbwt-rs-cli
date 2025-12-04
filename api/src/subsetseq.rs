//! A module for representing a sequence of subsets of an alphabet, with support for rank and select queries
//! on the elements of the subsets.

use std::ops::Range;

use bitvec::order::Lsb0;
use simple_sds_sbwt::bit_vector::*;
use simple_sds_sbwt::raw_vector::*;
use simple_sds_sbwt::ops::*;
use simple_sds_sbwt::serialize::*;

use crate::util::bitvec_to_simple_sds_raw_bitvec;

/// This trait represents a sequence of subsets from alphabet {0, 1, ..., sigma-1}, where sigma is the alphabet size.
/// The trait provides access to the subsets and rank and select queries for the elements inside the subsets.
#[allow(clippy::len_without_is_empty)]
pub trait SubsetSeq{

    // Todo: make char type into a generic unsigned integer.
    // Issues with that: it can't seem to figure out how to use such
    // generic integers as array indexes. Lol.

    /// Creates a new subset sequence from a vector of subsets of {0, 1, ..., sigma-1}. Order
    /// of characters inside a subset does not matter. The rank and select supports are not automatically
    /// initialized, so if the you need those functions, you need to call [`SubsetSeq::build_rank`] and [`SubsetSeq::build_select`], respectively.
    fn new(subset_seq: Vec<Vec<u8>>, sigma: usize) -> Self;

    /// Create a new subset sequence from indicator [bit vectors](bitvec::vec::BitVec), where the i-th bit of the j-th bit vector
    /// is 1 if and only if the i-th subset contains the j-th character. The resulting subset sequence has
    /// rank and select support if the provided bit vectors have rank and select support enabled. Otherwise, those
    /// supports need to be initialized by calling [`SubsetSeq::build_rank`] and [`SubsetSeq::build_select`], respectively.
    fn new_from_bit_vectors(vecs: Vec<bitvec::vec::BitVec::<u64, Lsb0>>) -> Self;

    /// Number of sets in the sequence (**not** the total length of the sets).
    fn len(&self) -> usize;

    /// Initialize rank support for the elements of the subsets.
    fn build_rank(&mut self);

    /// Initialize select support for the elements of the subsets.
    fn build_select(&mut self);

    /// Returns true if the select support has been initialized.
    fn has_select_support(&self) -> bool;

    /// Returns true if the rank support has been initialized.
    fn has_rank_support(&self) -> bool;

    /// Returns the number of sets in the range [0..i) that have character `c`.
    fn rank(&self, c: u8, i: usize) -> usize;

    /// Returns the index of the set that contains the i-th occurence of character c, if exists.
    fn select(&self, c: u8, i: usize) -> Option<usize>;

    /// Appends the elements of the i-th set to the buffer.
    fn append_set_to_buf(&self, i: usize, buf: &mut Vec<u8>);

    /// Returns the elements in the i-th set.
    fn access(&self, i: usize) -> Vec<u8>{
        let mut v = Vec::new();
        self.append_set_to_buf(i, &mut v);
        v
    }

    /// Returns the size of the i-th subset.
    fn subset_size(&self, i: usize) -> usize;

    /// Returns true if the set at index `set_idx` contains the character.
    fn set_contains(&self, set_idx: usize, character: u8) -> bool;

    /// Returns the index of the next set from set_idx (possibly set_idx itself) 
    /// that contains character c, or None if does not exist.
    /// The index set_idx can be from 0 to self.len().
    fn next_set_with_char(&self, set_idx: usize, c: u8) -> Option<usize>;

    /// Writes the subset sequence to the given writer.
    /// The sequence can be loaded later back with [`SubsetSeq::load`].
    /// Returns the number of bytes written.
    fn serialize<W: std::io::Write>(&self, out: &mut W) -> std::io::Result<usize>;

    /// Loads a subset sequence that was previously written with [`SubsetSeq::serialize`].
    fn load<R: std::io::Read>(input: &mut R) -> std::io::Result<Self> where Self: Sized;

    // Calls the callback on the index of every set that contains character c in the range.
    // The callback is called in increasing order of set index.
    fn call_on_char_occurrences<F: FnMut(usize)>(&self, range: Range<usize>, c: u8, callback: F);

    /// Build the cumulative sum array C required in [`crate::sbwt::SbwtIndex`].
    fn get_C_array(&self) -> Vec<usize> {
        let sigma = 4; // TODO
        let n = self.len();

        let mut C: Vec<usize> = vec![0; sigma];
        for i in 0..n {
            for c in 0..(sigma as u8) {
                if self.set_contains(i, c) {
                    for d in (c + 1)..(sigma as u8) {
                        C[d as usize] += 1;
                    }
                }
            }
        }

        // Plus one for the ghost dollar
        #[allow(clippy::needless_range_loop)] // Is perfectly clear this way
        for c in 0..sigma {
            C[c] += 1;
        }

        C
    }

    // This a key subroutine for SbwtIndex::push_labels_forward. We want to put it here in the SubsetSeq
    // trait because then we have have a really optimized version for particular implementations
    // like SubsetMatrix. See the code for what it does.
    fn push_labels_forward(&self, labels: &[u8], sbwt_input_range: std::ops::Range<usize>, output_ranges: Vec<&mut[u8]>);
}

/// An implementation of [SubsetSeq] with a matrix of sigma indicator bit vectors: the i-th bit of the j-th bit vector
/// is 1 if and only if the i-th subset contains the j-th character. Rank and select queries are reduced to
/// bit vector rank and select queries on the indicator bit vectors.
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct SubsetMatrix{
    rows: Vec<BitVector>
}

impl SubsetSeq for SubsetMatrix{
    
    fn serialize<W: std::io::Write>(&self, out: &mut W) -> std::io::Result<usize>{
        let mut n_written = 0_usize;

        let n_rows = self.rows.len() as u64;
        out.write_all(&n_rows.to_le_bytes())?;
        for row in self.rows.iter(){
            row.serialize(out)?;
            n_written += row.size_in_bytes();
        }
        Ok(n_written)
    }

    fn load<R: std::io::Read>(input: &mut R) -> std::io::Result<Self>{

        let n_rows = u64::load(input)? as usize;

        let mut rows = Vec::<BitVector>::new();
        for _ in 0..n_rows{
            rows.push(BitVector::load(input)?);
        }
        Ok(Self{rows})
    }

    fn new_from_bit_vectors(rows: Vec<bitvec::vec::BitVec<u64, Lsb0>>) -> Self{
        Self{rows: rows.into_iter().map(|row| 
            simple_sds_sbwt::bit_vector::BitVector::from(
                bitvec_to_simple_sds_raw_bitvec(row))
            ).collect()
        }
    }

    fn new(subset_seq: Vec<Vec<u8>>, sigma: usize) -> Self{
        let n = subset_seq.len();
        let mut rawrows = Vec::<RawVector>::new();
        for _ in 0..sigma{
            rawrows.push(RawVector::with_len(n, false));
        }
        for (i, set) in subset_seq.iter().enumerate(){
            for c in set.iter(){
                rawrows[*c as usize].set_bit(i, true)
            }
        }

        let rows: Vec<BitVector> = rawrows.into_iter().map(BitVector::from).collect();
        Self{rows}
    }

    fn build_rank(&mut self) {
        for row in self.rows.iter_mut(){
            row.enable_rank();
        }
    }

    fn build_select(&mut self) {
        for row in self.rows.iter_mut(){
            row.enable_select();
        }
    }

    fn has_select_support(&self) -> bool {
        self.rows[0].supports_select()
    }

    fn has_rank_support(&self) -> bool {
        self.rows[0].supports_rank()
    }

    fn rank(&self, c: u8, i: usize) -> usize{
        self.rows[c as usize].rank(i)
    }

    fn select(&self, c: u8, i: usize) -> Option<usize>{
        assert!(self.rows[c as usize].supports_select());
        self.rows[c as usize].select(i)
    }

    fn len(&self) -> usize{
        if self.rows.is_empty() { 0 }
        else { self.rows[0].len() }
    }

    fn subset_size(&self, i: usize) -> usize {
        self.rows.iter().filter(|&row| row.get(i)).count()
    }

    fn set_contains(&self, set_idx: usize, character: u8) -> bool {
        self.rows[character as usize].get(set_idx)
    }

    fn append_set_to_buf(&self, i: usize, buf: &mut Vec<u8>){
        for c in 0..self.rows.len(){
            if self.rows[c].get(i){
                buf.push(c as u8);
            }
        }
    }

    fn next_set_with_char(&self, set_idx: usize, c: u8) -> Option<usize> {
        // Use the bitvec crate on the raw data of the row
        let bv = bitvec::slice::BitSlice::<u64, Lsb0>::from_slice(self.rows[c as usize].as_ref().as_ref());
        let off = bv[set_idx..].first_one()?;
        Some(set_idx + off)
    }

    fn call_on_char_occurrences<F: FnMut(usize)>(&self, range: Range<usize>, c: u8, mut callback: F) {
        let bv = bitvec::slice::BitSlice::<u64, Lsb0>::from_slice(self.rows[c as usize].as_ref().as_ref());
        let slice = &bv[range.clone()];
        for i in slice.iter_ones() {
            callback(range.start + i);
        }
    }

    fn push_labels_forward(&self, labels: &[u8], sbwt_input_range: std::ops::Range<usize>, output_ranges: Vec<&mut[u8]>) {
        assert_eq!(labels.len(), sbwt_input_range.len());
        assert_eq!(output_ranges.len(), self.rows.len());
        if sbwt_input_range.is_empty() { return }

        #[allow(clippy::needless_range_loop)]
        for (char_idx, output_slice) in output_ranges.into_iter().enumerate() {
            let mut output_offset = 0_usize;

            let words = self.rows[char_idx].as_ref().get_words();

            // We roll our own one-bit iterator because the one in the bitvec crate
            // has a lot of overhead. This one is something like 6x faster!
            for_each_one_bit(words, sbwt_input_range.start, sbwt_input_range.end, |i| {
                output_slice[output_offset] = labels[i - sbwt_input_range.start];
                output_offset += 1;
            });
        }
    }

}

fn for_each_one_bit<F>(words: &[u64], s: usize, e: usize, mut cb: F)
where
    F: FnMut(usize),
{
    if words.is_empty() || s >= e {
        return;
    }

    let bit_len = words.len() * 64;
    if s >= bit_len {
        return;
    }

    let start_word = s / 64;
    let end_word = (e - 1) / 64; // inclusive

    // Helper: iterate 1-bits in a single word, adding `base_bit` to positions.
    #[inline]
    fn scan_word<F>(mut word: u64, base_bit: usize, cb: &mut F)
    where
        F: FnMut(usize),
    {
        while word != 0 {
            let tz = word.trailing_zeros() as usize;
            cb(base_bit + tz);
            word &= word - 1; // clear lowest set bit
        }
    }

    // Single-word range: mask both ends in the same word.
    if start_word == end_word {
        let start_bit = s % 64;
        let end_bit = e % 64; // exclusive

        let mut word = words[start_word];

        // Mask out bits before s
        word &= !0u64 << start_bit;

        // Mask out bits at/after e
        if end_bit != 0 {
            word &= (1u64 << end_bit) - 1;
        }

        scan_word(word, start_word * 64, &mut cb);
        return;
    }

    // --- First word: mask only the low bits ---
    {
        let start_bit = s % 64;
        let mut word = words[start_word];

        word &= !0u64 << start_bit;
        scan_word(word, start_word * 64, &mut cb);
    }

    // --- Middle words: no masking at all ---
    for (rel_idx, &word) in words[start_word + 1..end_word].iter().enumerate() {
        scan_word(word, (start_word + 1 + rel_idx) * 64, &mut cb);
    }

    // --- Last word: mask only the high bits ---
    {
        let end_bit = e % 64; // exclusive; 0 means "all 64 bits are in range"
        let mut word = words[end_word];

        if end_bit != 0 {
            word &= (1u64 << end_bit) - 1;
        }
        scan_word(word, end_word * 64, &mut cb);
    }
}

/// Formats the subset matrix as an ASCII bit matrix of 0s and 1s, where row i
/// is the indicator bit vector for the i-th character.
impl std::fmt::Display for SubsetMatrix{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for row in self.rows.iter(){
            for i in 0..row.len(){
                if row.get(i){
                    write!(f, "1")?;
                } else{
                    write!(f, "0")?;
                }
            }
            writeln!(f)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn serialize_and_load(){
        let sets: Vec<Vec<u8>> = vec![vec![1,2,3], vec![0,2], vec![0,1,3,4], vec![], vec![0,1,2]];
        let mut sm = SubsetMatrix::new(sets, 5);
        sm.build_rank();
        let mut buf = Vec::<u8>::new();
        sm.serialize(&mut buf).unwrap();
        let sm2 = SubsetMatrix::load(&mut buf.as_slice()).unwrap();
        assert_eq!(sm, sm2);
    }

    #[test]
    fn next_set_with_char(){
        let sets: Vec<Vec<u8>> = vec![vec![1,2,3], vec![0,2], vec![0,1,3,4], vec![], vec![0,1,2]];
        let sm = SubsetMatrix::new(sets, 5);
        assert_eq!(sm.next_set_with_char(0,0), Some(1));
        assert_eq!(sm.next_set_with_char(0,1), Some(0));
        assert_eq!(sm.next_set_with_char(0,2), Some(0));
        assert_eq!(sm.next_set_with_char(0,3), Some(0));
        assert_eq!(sm.next_set_with_char(0,4), Some(2));

        assert_eq!(sm.next_set_with_char(3,0), Some(4));
        assert_eq!(sm.next_set_with_char(3,1), Some(4));
        assert_eq!(sm.next_set_with_char(3,2), Some(4));
        assert_eq!(sm.next_set_with_char(3,3), None);
        assert_eq!(sm.next_set_with_char(3,4), None);

        assert_eq!(sm.next_set_with_char(5,0), None); // One past the end
    }
}
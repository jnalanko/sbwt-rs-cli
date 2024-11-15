//! Miscellaneous utility functions and constants used in the crate.

use bitvec::prelude::*;
use simple_sds_sbwt::bit_vector::BitVector;
use simple_sds_sbwt::ops::Select;
use rayon::iter::IntoParallelRefIterator;
use rayon::iter::IndexedParallelIterator;
use rayon::iter::ParallelIterator;

// Returns the number of bytes written
pub(crate) fn write_bytes<W: std::io::Write>(out: &mut W, bytes: &[u8]) -> std::io::Result<usize>{
    out.write_all(bytes)?;
    Ok(bytes.len() + 8)
}

// Searcher the range [0..n)
// Return the index of the answer, or n if does not exist
pub(crate) fn binary_search_leftmost_that_fulfills_pred<T, Access: Fn(usize) -> T, Pred: Fn(T) -> bool>(access: Access, pred: Pred, n: usize) -> usize {
    let mut ans = n;
    let mut step = n;
    while step > 0 {
        while ans as isize - step as isize >= 0 && pred(access(ans-step)) {
            ans -= step;
        }
        step /= 2;
    }
    ans
}

pub const DNA_ALPHABET: [u8; 4] = [b'A', b'C', b'G', b'T'];

// This bit vector of length 256 marks the ascii values of these characters: acgtACGT
const IS_DNA: BitArray<[u32; 8]> = bitarr![const u32, Lsb0; 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];

pub(crate) fn is_dna(c: u8) -> bool {
    IS_DNA[c as usize]
}

// A table mapping mapping ascii A -> 0, C -> 1, G -> 2, T -> 3. Same for lower case.
// All other charcters map to 255. Other code depends on this choice: don't touch it. 
pub const ACGT_TO_0123: [u8; 256] = [255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 0, 255, 1, 255, 255, 255, 2, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 3, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 0, 255, 1, 255, 255, 255, 2, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 3, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255];

#[allow(non_snake_case)] // C-array is an established convention in BWT indexes
pub(crate) fn get_C_array(rawrows: &[simple_sds_sbwt::raw_vector::RawVector]) -> Vec<usize> {
    let sigma = rawrows.len();
    assert!(sigma > 0);

    let mut C: Vec<usize> = rawrows.par_iter().enumerate().map(|(c, rawrow)| {
        let bv = BitVector::from(rawrow.clone());
        let mut C: Vec<usize> = vec![0; sigma];
        bv.one_iter().for_each(|_| {
            for d in (c + 1)..(sigma) {
                C[d as usize] += 1;
            }
        });
        C
    }).reduce(|| vec![0; sigma], |a, b| a.iter().zip(b.iter()).map(|(x, y)| x + y).collect());

    // Plus one for the ghost dollar
    #[allow(clippy::needless_range_loop)] // Is perfectly clear this way
    for c in 0..sigma {
        C[c] += 1;
    }

    C
}

/// Reverses the given ASCII DNA sequence and replaces each nucleotide with its complement.
pub fn reverse_complement_in_place(seq: &mut [u8]){
    jseqio::reverse_complement_in_place(seq);
}

#[allow(dead_code)]
pub(crate) struct FastXReader{
    inner: jseqio::reader::DynamicFastXReader
}

impl crate::SeqStream for FastXReader{
    fn stream_next(&mut self) -> Option<&[u8]> {
        self.inner.read_next().unwrap().map(|x| x.seq)
    }
}

/// Creates a [crate::SeqStream] out of a slice of ASCII sequences.
pub struct SliceSeqStream<'a>{
    slices: &'a [ &'a [u8]],
    cur_slice_idx: usize,
}

impl<'a> SliceSeqStream<'a> {
    /// Creates a [crate::SeqStream] out of a slice of ASCII sequences.
    pub fn new(slices: &'a [&'a [u8]]) -> Self {
        Self {slices, cur_slice_idx: 0}
    }
}

impl<'a> crate::SeqStream for SliceSeqStream<'a> {
    fn stream_next(&mut self) -> Option<&[u8]> {
        if self.cur_slice_idx == self.slices.len() {
            None
        } else {
            let s = self.slices[self.cur_slice_idx];
            self.cur_slice_idx += 1;
            Some(s)
        }
    } 
}

/// Creates a [crate::SeqStream] out of a slice of ascii vectors.
pub struct VecSeqStream<'a> {
    seqs: &'a [Vec<u8>],
    cur_seq_idx: usize,
}

impl<'a> VecSeqStream<'a> {
    /// Creates a [crate::SeqStream] out of a slice of ascii vectors.
    pub fn new(seqs: &'a [Vec<u8>]) -> Self {
        Self {seqs, cur_seq_idx: 0}
    }
}

impl<'a> crate::SeqStream for VecSeqStream<'a> {
    fn stream_next(&mut self) -> Option<&[u8]> {
        if self.cur_seq_idx == self.seqs.len() {
            None
        } else {
            let s = &self.seqs[self.cur_seq_idx];
            self.cur_seq_idx += 1;
            Some(s)
        }
    } 
}

//! Miscellaneous utility functions and constants used in the crate.

use std::cmp::min;
use std::io::Cursor;
use std::io::Read;
use std::ops::Range;

use bitvec::prelude::*;
use rayon::iter::ParallelBridge;
use simple_sds_sbwt::bit_vector::BitVector;
use simple_sds_sbwt::ops::Select;
use rayon::iter::IntoParallelRefIterator;
use rayon::iter::IndexedParallelIterator;
use rayon::iter::ParallelIterator;
use simple_sds_sbwt::raw_vector::AccessRaw;
use simple_sds_sbwt::serialize::Serialize;

type BitVec = bitvec::vec::BitVec<u64, Lsb0>;
type BitSlice = bitvec::slice::BitSlice<u64, Lsb0>;

// Returns the number of bytes written
pub(crate) fn write_bytes<W: std::io::Write>(out: &mut W, bytes: &[u8]) -> std::io::Result<usize>{
    out.write_all(bytes)?;
    Ok(bytes.len() + 8)
}

// Searcher the range [0..n)
// Return the index of the answer, or n if does not exist
pub(crate) fn binary_search_leftmost_that_fulfills_pred<T, Access: Fn(usize) -> T, Pred: Fn(T) -> bool>(access: Access, pred: Pred, n: usize) -> usize {
    let mut ans = n;
    let mut step = n.next_power_of_two();
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
fn popcounts_to_C_array(popcounts: &[usize]) -> Vec<usize> {
    let sigma = popcounts.len();
    let mut C: Vec<usize> = vec![0; sigma];
    for c in 0..sigma {
        for d in (c + 1)..sigma {
            C[d] += popcounts[c];
        }
    }

    // Plus one for the ghost dollar
    #[allow(clippy::needless_range_loop)] // Is perfectly clear this way
    for c in 0..sigma {
        C[c] += 1;
    }

    C

}

#[allow(non_snake_case)] // C-array is an established convention in BWT indexes
pub(crate) fn get_C_array(rawrows: &[simple_sds_sbwt::raw_vector::RawVector]) -> Vec<usize> {
    let popcounts: Vec<usize> = rawrows.iter().map(|row| row.count_ones()).collect(); 
    popcounts_to_C_array(&popcounts)
}

#[allow(non_snake_case, dead_code)] // C-array is an established convention in BWT indexes
pub(crate) fn get_C_array_parallel(rawrows: &[simple_sds_sbwt::raw_vector::RawVector], n_threads: usize) -> Vec<usize> {
    let sigma = rawrows.len();
    assert!(sigma > 0);

    let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();
    thread_pool.install(||{
        let popcounts: Vec<usize> = rawrows.par_iter().map(|row| row.count_ones()).collect(); 
        popcounts_to_C_array(&popcounts)
    })
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

impl crate::SeqStream for SliceSeqStream<'_> {
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

impl crate::SeqStream for VecSeqStream<'_> {
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


// This function runs in parallel, so a rayon thread pool must be initialized.
pub(crate) fn parallel_bitvec_concat(bitvecs: Vec<BitVec>) -> BitVec {
    if bitvecs.len() == 1 {
        return bitvecs.into_iter().next().unwrap(); // Nothing to do
    }
    let total_length = bitvecs.iter().fold(0_usize, |acc, s| acc + s.len());
    let n_words = total_length.div_ceil(64);
    let mut output_data = vec![0_u64; n_words];

    let mut exclusive_input_bitslices = Vec::<&BitSlice>::new(); // One per nonempty input bitvec
    let mut exclusive_output_word_ranges = Vec::<Range<usize>>::new(); // One per nonempty input bitvec

    // Figure out which words will be potentially shared between input slices in the concatenation,
    // and set those bits. Push the exclusive input and output regions to vecs.
    let mut bits_so_far = 0_usize;
    for s in bitvecs.iter() {
        if s.is_empty() {
            continue // Nothing to concatenate
        }

        // Copy the part that falls in the first output word
        let first_word = bits_so_far / 64; // The first word that will be written to
        let first_word_bit_offset = bits_so_far % 64; // Index of the first bit that is written in the first word
        let n_bits_written_to_first_word = min(s.len(), 64 - first_word_bit_offset);
        let first_out = BitSlice::from_element_mut(&mut output_data[first_word]);
        let first_out_slice = &mut first_out[first_word_bit_offset..first_word_bit_offset + n_bits_written_to_first_word];
        first_out_slice.copy_from_bitslice(&s[0..n_bits_written_to_first_word]);

        // Copy the part that falls in the last output word (can be the same
        // as the first output word but that's okay).
        bits_so_far += s.len();
        let last_bit = bits_so_far - 1; // >= 0 because s.len() > 0
        let last_word = last_bit / 64;
        let last_word_bit_offset = last_bit % 64; // Index of the last bit that is written in the last word
        let n_bits_written_to_last_word = min(s.len(), last_word_bit_offset + 1);
        let last_out_word = BitSlice::from_element_mut(&mut output_data[last_word]);
        let last_out_slice = &mut last_out_word[(1 + last_word_bit_offset - n_bits_written_to_last_word)..=last_word_bit_offset];
        let last_in_slice = &s[(s.len() - n_bits_written_to_last_word)..];
        last_out_slice.copy_from_bitslice(last_in_slice);

        let first_exclusive_word = (first_word + 1) as isize;
        let last_exclusive_word = last_word as isize - 1;

        if first_exclusive_word <= last_exclusive_word {
            // Nonempty range of exclusive words
            exclusive_input_bitslices.push(&s[n_bits_written_to_first_word..(s.len()-n_bits_written_to_last_word)]);
            exclusive_output_word_ranges.push(first_exclusive_word as usize .. last_exclusive_word as usize +1);
        }

    }

    let exclusive_output_word_ranges = split_to_mut_regions(&mut output_data, exclusive_output_word_ranges);

    // Copy non-overlapping parts in parallel
    assert_eq!(exclusive_input_bitslices.len(), exclusive_output_word_ranges.len());
    exclusive_output_word_ranges.into_iter().enumerate().par_bridge().for_each(|(i, out_word_range)| {
        let out = BitSlice::from_slice_mut(out_word_range);
        out.copy_from_bitslice(exclusive_input_bitslices[i]);
    });

   let mut concat = BitVec::from_vec(output_data); // Reinterpret as BitVec
   concat.truncate(total_length); // Get rid of the tail in the last word
   concat

}

/// Returns a mutable slice for every region in `regions`.
///
/// ## Preconditions  (guaranteed by the caller)
/// * All ranges are inside `0..v.len()`.
/// * The ranges are sorted by `start` and do **not** overlap.
pub(crate) fn split_to_mut_regions(
    v: &mut Vec<u64>,
    regions: Vec<Range<usize>>,
) -> Vec<&mut [u64]> {
    let mut result = Vec::with_capacity(regions.len());
    let mut tail: &mut [u64] = v.as_mut_slice();
    let mut consumed = 0; // absolute index we have reached in `v`

    for r in regions {
        // translate the absolute `start` to an index inside `tail`
        let rel_start = r.start - consumed;
        let len       = r.end   - r.start;

        // First split off everything before the wanted region …
        let (_, after_start)    = tail.split_at_mut(rel_start);
        // ... then split that remainder into the desired region and the rest.
        let (region, after_end) = after_start.split_at_mut(len);

        result.push(region); // keep the region
        tail = after_end; // keep working with the suffix
        consumed = r.end; // advance the absolute cursor
    }

    result
}

pub(crate) fn bitvec_to_simple_sds_raw_bitvec(bv: bitvec::vec::BitVec<u64, Lsb0>) -> simple_sds_sbwt::raw_vector::RawVector {
    // Let's use the deserialization function in simple_sds_sbwt for a raw bitvector.
    // It requires the following header:
    let mut header = [0u64, 0u64]; // bits, words
    header[0] = bv.len() as u64; // Assumes little-endian byte order
    header[1] = bv.len().div_ceil(64) as u64;

    let header_bytes = bytemuck::cast_slice(&header);
    let raw_data = bytemuck::cast_slice(&bv.as_raw_slice());
    let mut data_with_header = Cursor::new(header_bytes).chain(Cursor::new(raw_data));

    simple_sds_sbwt::raw_vector::RawVector::load(&mut data_with_header).unwrap()

}

pub(crate) fn segment_range(range: Range<usize>, n_pieces: usize) -> Vec<Range<usize>> {
        let segment_len = range.len().div_ceil(n_pieces);
        let mut pieces: Vec<Range<usize>> = vec![];
        for t in 0..n_pieces{
            pieces.push(range.start + t*segment_len .. 
                        range.start + min((t+1)*segment_len, range.len())); // Final segments may be empty. Is ok.
        }
        pieces
}

#[cfg(test)]
pub(crate) fn gen_random_dna_string(len: usize, seed: u64) -> Vec<u8> {
    use rand_chacha::rand_core::{RngCore, SeedableRng};

    let mut rng = rand_chacha::ChaCha20Rng::seed_from_u64(1234);
    (0..len).map(|_| { 
        match rng.next_u64() % 4 {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            3 => b'T',
            _ => panic!("Impossible")
        }
    }).collect()
}


#[cfg(test)]
mod tests {

    use rand::{Rng, SeedableRng};

    use super::*;

    #[test]
    fn parallel_concatenation() {
        // Generate 1001 pseudorandom bit vectors with lengths between 0 and 256 bits
        let seed = 42;
        let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
        let mut bitvecs: Vec<BitVec> = Vec::new();
        for i in 0..1000 {
            let len = rng.gen_range(0, 256);
            if i == 500 {
                // Add an empty bit vector
                bitvecs.push(BitVec::new());
            }
            let mut bits = BitVec::with_capacity(len);
            for _ in 0..len {
                bits.push(rng.gen_bool(0.5));
            }
            bitvecs.push(bits);
        }

        let true_concat = bitvecs.iter().fold(BitVec::new(), |mut acc, v| {
            acc.extend_from_bitslice(v);
            acc
        });

        let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(3).build().unwrap();
        let our_concat = thread_pool.install(|| {
            parallel_bitvec_concat(bitvecs)
        });

        assert_eq!(true_concat.len(), our_concat.len());
        assert_eq!(true_concat, our_concat);

    }

}
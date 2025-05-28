use std::cmp::min;
use std::ops::Range;

use bitvec::vec::BitVec;
use bitvec::prelude::*;
use rayon::iter::IntoParallelIterator;
use rayon::iter::IntoParallelRefIterator;
use rayon::iter::ParallelIterator;
use crate::subsetseq::*;
use crate::sbwt::*;

pub struct MergeInterleaving {
    // Has one bit per colex position in the merged SBWT
    // s1[i] is 1 iff this k-mer is in the first SBWT
    pub s1: BitVec, 

    // Has one bit per colex position in the merged SBWT
    // s2[i] is 1 iff this k-mer is in the second SBWT
    pub s2: BitVec,

    // Has one bit per colex position in the merged SBWT, marking
    // the dummy nodes.
    pub is_dummy: BitVec, 
    
    // Has one bit per colex position in the merged SBWT, marking
    // the positions whose k-mer has a different (k-1)-length suffix
    // than the previous k-mer in colex order. 
    // Edge case: is_leader[0] = 1.
    pub is_leader: BitVec,
}

impl MergeInterleaving {

    // Parallel version of split_to_pieces
    fn split_to_pieces_par(s: &BitSlice, n_pieces: usize, n_threads: usize) -> Vec<(usize, Range<usize>)> {
        const BLOCK_SIZE: usize = 1024;

        let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();
        thread_pool.install(||{
            // Strategy: compute block popcounts in parallel, then locate the blocks where the
            // pieces start, and do bit-by-bit counting to find the precise start points of pieces.

            assert!(n_pieces > 0);
            if s.len() > 0 {
                // The last bit should always be 1
                assert!(s.last().unwrap() == true);
            }

            let blocks = s.chunks(BLOCK_SIZE).collect::<Vec::<&BitSlice>>();
            let block_popcounts: Vec<usize> = blocks.par_iter().map(|block| block.count_ones()).collect();
            let total_popcount: usize = block_popcounts.iter().sum();
            let ones_per_piece = total_popcount.div_ceil(n_pieces); // Last piece may have fewer

            let mut starts : Vec<usize> = vec![0]; // Init with the start of the first piece
            let mut n_zeros_before_piece: Vec<usize> = vec![0]; // Init with #zeros before the first piece

            let mut n_ones = 0_usize;
            let mut n_zeros = 0_usize;

            // Let N be the number of ones in a piece. The i-th block piece just after
            // the one-bit with zero-based rank N*i - 1. That is, if N = 10, then the fifth
            // piece starts at the one-bit with rank 49 (= the 50th 1-bit)
            for (block_idx, block) in blocks.iter().enumerate() {
                while starts.len() < n_pieces && n_ones + block_popcounts[block_idx] >= starts.len() * ones_per_piece {
                    // The check for starts.len() < n_pieces is to avoid creating an empty piece after
                    // the last one in case the total popcount is divisible by n_pieces.

                    // the 1-bit just before the start the next piece is in this block.
                    // Find where it is
                    let mut n_ones_precise = n_ones;
                    let mut n_zeros_precise = n_zeros;
                    let target = starts.len() * ones_per_piece;
                    let mut i = 0_usize;
                    loop {
                        n_ones_precise += block[i] as usize;
                        n_zeros_precise += 1 - (block[i] as usize);
                        if n_ones_precise == target {
                            break;
                        } else {
                            i += 1;
                        }
                    }
                    assert_eq!(n_ones_precise, target);
                    starts.push(n_ones + n_zeros + i + 1);
                    n_zeros_before_piece.push(n_zeros_precise);
                }
                n_ones += block_popcounts[block_idx];
                n_zeros += block.len() - block_popcounts[block_idx];
            }
            assert_eq!(starts.len(), n_zeros_before_piece.len());
            assert!(starts.len() > 0);
            while starts.len() < n_pieces {
                // Add empty pieces to the end
                starts.push(s.len());
                n_zeros_before_piece.push(s.len() - total_popcount);
            }

            let mut pieces: Vec<(usize, Range<usize>)> = vec![];
            assert_eq!(starts.len(), n_pieces);
            starts.push(s.len()); // End sentinel for the end of the last range
            for i in 0..n_pieces {
                pieces.push((n_zeros_before_piece[i], starts[i]..starts[i+1]));
            }

            pieces
        })
    }

    // Returns a segmentation s[l_1..r_1), s[l_2..r_2], ... such that
    // each segment ends in a 1-bit and has an approximately equal number of 1-bits.
    // Also returns the number of 0-bits before each segment.
    fn split_to_pieces(s: &BitSlice, n_pieces: usize) -> Vec<(usize, Range<usize>)> {
        assert!(n_pieces > 0);
        if s.len() > 0 {
            // The last bit should always be 1
            assert!(s.last().unwrap() == true);
        }
        let s_popcount = s.count_ones();
        let ones_per_piece = s_popcount.div_ceil(n_pieces); // Last pieces may have fewer

        let mut s_ranges: Vec<Range<usize>> = vec![];
        let mut zeros_before: Vec<usize> = vec![];

        let mut cur_popcount = 0_usize;
        let mut n_zeros = 0_usize;
        for s_idx in s.iter_ones(){
            cur_popcount += 1;
            if cur_popcount == ones_per_piece || s_idx == s.len() - 1 {
                let prev_end = match s_ranges.last() {
                    Some(r) => r.end,
                    None => 0,
                };
                let new_end = s_idx + 1;
                s_ranges.push(prev_end..new_end);

                zeros_before.push(n_zeros);
                n_zeros += (new_end - prev_end) - cur_popcount;

                cur_popcount = 0;
            }
        }
        
        while s_ranges.len() < n_pieces {
            s_ranges.push(s.len()..s.len()); // Empty range
            zeros_before.push(n_zeros);
        }

        // Pair of vectors into vector of pairs
        zeros_before.into_iter().zip(s_ranges).collect()
        
    }

    pub fn new<SS: SubsetSeq + Send + Sync>(index1: &SbwtIndex::<SS>, index2: &SbwtIndex<SS>, n_threads: usize) -> MergeInterleaving {

        let k = index1.k();
        assert_eq!(k, index2.k());

        // We invert the SBWTs column by column and maintain ranges
        // in both SBWTs that are so far equal. The ranges are in increasing
        // order and the partition the SBWTs, to it's enough to just store their
        // sizes in order. The sizes are stored in concatenated unary representations.
        // Empty ranges are allowed.

        let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();
        thread_pool.install(||{
            let mut leader_bits = None;

            // Initialize unary concatenations with empty ranges
            let mut s1 = bitvec![0; index1.n_sets()];
            let mut s2 = bitvec![0; index2.n_sets()];
            s1.push(true);
            s2.push(true);

            let mut chars1 = index1.build_last_column();
            let mut chars2 = index2.build_last_column();

            let mut temp_char_buf_1 = vec![0u8; chars1.len()];
            let mut temp_char_buf_2 = vec![0u8; chars2.len()];

            for round in 0..k {
                log::info!("Round {}/{}", round+1, k);

                // Split work into pieces for different threads
                log::info!("Splitting work");
                let p1 = Self::split_to_pieces_par(&s1, n_threads, n_threads);
                let p2 = Self::split_to_pieces_par(&s2, n_threads, n_threads);
                assert_eq!(p1.len(), n_threads);
                assert_eq!(p2.len(), n_threads);
                // Zip pairs of tuples into 4-tuples
                let pieces = (0..n_threads).map(|i| (p1[i].0, p2[i].0, p1[i].1.clone(), p2[i].1.clone())).collect();

                log::info!("Refining segmentation");
                let new_arrays = Self::refine_segmentation(s1, s2, &chars1, &chars2, pieces, round == k-1);
                (s1, s2, leader_bits) = new_arrays;

                if round != k-1 {
                    log::info!("Pushing labels forward in the SBWT graph");
                    index1.push_all_labels_forward(&chars1, &mut temp_char_buf_1, n_threads);
                    std::mem::swap(&mut chars1, &mut temp_char_buf_1);

                    index2.push_all_labels_forward(&chars2, &mut temp_char_buf_2, n_threads);
                    std::mem::swap(&mut chars2, &mut temp_char_buf_2);
                }
            }

            let leader_bits = leader_bits.unwrap(); // Compute in the last round
            log::info!("Number of suffix groups: {}", leader_bits.count_ones());

            Self::compress_in_place(&mut s1);
            Self::compress_in_place(&mut s2);

            assert_eq!(s1.len(), s2.len());
            let merged_len = s1.len();

            // Identify dummies in the merged SBWT
            let mut is_dummy = bitvec::vec::BitVec::new();
            is_dummy.resize(merged_len, false);
            let mut c1_idx = 0_usize;
            let mut c2_idx = 0_usize;
            for colex in 0..merged_len {
                let d1 = s1[colex] && c1_idx < chars1.len() && chars1[c1_idx] == b'$';
                let d2 = s2[colex] && c2_idx < chars2.len() && chars2[c2_idx] == b'$';
                is_dummy.set(colex, d1 || d2); 

                c1_idx += s1[colex] as usize;
                c2_idx += s2[colex] as usize;
            }

            log::info!("Number of dummies: {}", is_dummy.count_ones());


            MergeInterleaving { s1, s2, is_dummy, is_leader: leader_bits}
        })
    }

    pub fn intersection_size(&self) -> usize {
        assert_eq!(self.s1.len(), self.s2.len());
        let mut ans = 0_usize;
        for i in 0..self.s1.len() {
            // Do not count dummy nodes
            ans += (!self.is_dummy[i] && self.s1[i] && self.s2[i]) as usize;
        }
        ans
    }

    pub fn union_size(&self) -> usize {
        assert_eq!(self.s1.len(), self.s2.len());
        let mut ans = 0_usize;
        for i in 0..self.s1.len() {
            // Do not count dummy nodes
            ans += (!self.is_dummy[i] && (self.s1[i] || self.s2[i])) as usize;
        }
        ans
    }

    // Helper for construction. Takes a unary concatenation of binary numbers 0^b 1, like:
    // 101011011101 (= 0,1,1,0,1,0,0,1)
    // And produces a bit vector encoding the bits:
    // 01101001
    fn compress_in_place(s1: &mut bitvec::vec::BitVec) {

            let mut s1_i = 0_usize; // Index in s1

            let mut new_idx = 0_usize;
            while s1_i < s1.len() {
                let len1 = Self::leading_zeros(&s1[s1_i..]);
                assert!(len1 <= 1); // This is the colex range of a k-mer, so it should be empty or singleton

                s1_i += len1 + 1; // Length of the unary number we just parsed

                s1.set(new_idx, len1 != 0);
                new_idx += 1;
            }
            assert_eq!(s1_i, s1.len());
            s1.resize(new_idx, false);
            s1.shrink_to_fit();

    }

    // Assumes there is at least one 1-bit
    fn leading_zeros(s: &BitSlice) -> usize {
        s.first_one().unwrap()
    }

    // Assumes that seq is a string with some number of c in the beginning (possibly none)
    // followed by a tail of non-c characters which are all larger than c.
    // Returns the number of c in the beginning.
    // Falls back to binary search for long sequences
    #[inline]
    fn run_length_in_sorted_seq(seq: &[u8], c: u8) -> usize {
        if seq.len() == 0 {
            return 0;
        }
        if seq.len() > 200 {
            // Binary search the first element that is larger than c
            let end = seq.binary_search_by(|&x| {
                if x > c {
                    std::cmp::Ordering::Greater
                } else {
                    std::cmp::Ordering::Less
                }
            }).unwrap_err(); // Is always Err because we never return Ordering::Equal
            end
        } else {
            // Linear scan
            let mut i = 0;
            while i < seq.len() && seq[i] == c {
                i += 1;
            }
            i
        }
    }

    #[inline]
    fn zero_extend(v: &mut BitVec, howmany: usize) {
        if howmany >= 32 { // Faster extension in bulk (up to 64 bits at a time)
            v.resize(v.len() + howmany, false);
        } else {
            for _ in 0..howmany {
                v.push(false)
            }
        }
    }


    // Helper of a helper function
    fn refine_piece(s1: &BitVec, s2: &BitVec, chars1: &[u8], chars2: &[u8], mut c1_i: usize, mut c2_i: usize, s1_range: Range<usize>, s2_range: Range<usize>, compute_leaders: bool) 
    -> (BitVec, BitVec, Option<BitVec>) {

        let mut out1 = bitvec::bitvec![];
        let mut out2 = bitvec::bitvec![];

        let mut leader_bits = bitvec::bitvec![]; // Last round only

        // c1_i and c2_i are current indices in chars1 and chars2 respectively
        // s1_i and s2_i are current indices in s1 and s2 respectively
        let mut s1_i = s1_range.start;
        let mut s2_i = s2_range.start;

        while s1_i < s1_range.end {

            assert!(s2_i < s2_range.end);

            let len1 = Self::leading_zeros(&s1[s1_i..]);
            let len2 = Self::leading_zeros(&s2[s2_i..]);

            let c1_end = c1_i + len1; // One past the end
            let c2_end = c2_i + len2; // One past the end

            let mut is_leader = true;
            while c1_i < c1_end || c2_i < c2_end {
                let c1 = if c1_i == c1_end { u8::MAX } else { chars1[c1_i] };
                let c2 = if c2_i == c2_end { u8::MAX } else { chars2[c2_i] };
                let c = min(c1,c2);

                let r1 = Self::run_length_in_sorted_seq(&chars1[c1_i..c1_end], c);
                let r2 = Self::run_length_in_sorted_seq(&chars2[c2_i..c2_end], c);

                // Write r1 and r2 in unary
                Self::zero_extend(&mut out1, r1);
                Self::zero_extend(&mut out2, r2);
                // Terminate unary representations
                out1.push(true);
                out2.push(true);

                // Advance indexes in chars
                c1_i += r1;
                c2_i += r2;

                if compute_leaders {
                    leader_bits.push(is_leader);
                }

                is_leader = false;
            }

            assert_eq!(c1_i, c1_end);
            assert_eq!(c2_i, c2_end);

            s1_i += len1 + 1;
            s2_i += len2 + 1;
        }

        assert_eq!(s1_i, s1_range.end);
        assert_eq!(s2_i, s2_range.end);

        (out1, out2, if compute_leaders { Some(leader_bits) } else { None })
    }

    // Helper function for construction.
    // Input pieces are pairs (start in chars1, start in chars2, range in s1, range in s2)
    // Input piece are for parallelism: one piece to work on for each thread.
    // Returns the new segmentation.
    // On the last round also returns the leader bit vector, which marks
    // the smallest k-mer in each group of k-mers with the same suffix of length (k-1).
    // This function runs in parallel, so a rayon thread pool must be initialized.
    fn refine_segmentation(s1: BitVec, s2: BitVec, chars1: &[u8], chars2: &[u8], input_pieces: Vec<(usize, usize, Range<usize>, Range<usize>)>, last_round: bool) -> (BitVec, BitVec, Option<BitVec>) {
        log::info!("Refining pieces in parallel");
        let output_pieces: Vec<(BitVec, BitVec, Option<BitVec>)> = input_pieces.par_iter().map(|piece| {
            let (c1_i, c2_1, s1_range, s2_range) = piece;
            Self::refine_piece(&s1, &s2, chars1, chars2, *c1_i, *c2_1, s1_range.clone(), s2_range.clone(), last_round)
        }).collect();

        log::info!("Freeing memory");
        // Free memory
        drop(s1);
        drop(s2);

        // Concatenate pieces
        let new_s1_len = output_pieces.iter().fold(0_usize, |acc, v| acc + v.0.len());
        let new_s2_len = output_pieces.iter().fold(0_usize, |acc, v| acc + v.1.len());
        log::info!("Concatenating pieces (total length {})", new_s1_len + new_s2_len);
        let mut new_s1 = bitvec![];
        let mut new_s2 = bitvec![];
        new_s1.reserve_exact(new_s1_len);
        new_s2.reserve_exact(new_s2_len);
        let mut new_leader_bits = if last_round { Some(bitvec![]) } else { None }; // Todo: reserve up front
        for (s1_piece, s2_piece, leader_piece) in output_pieces {
            new_s1.extend_from_bitslice(&s1_piece);
            new_s2.extend_from_bitslice(&s2_piece);
            if let Some(v) = new_leader_bits.as_mut() {
                v.extend_from_bitslice(&leader_piece.unwrap());
            }
        }

        (new_s1, new_s2, new_leader_bits)
    }

}

fn parallel_bitslice_concat(slices: &[&BitSlice::<u64, Lsb0>]) -> BitVec<u64, Lsb0> {
    let total_length = slices.iter().fold(0_usize, |acc, s| acc + s.len());
    let n_words = total_length.div_ceil(64);
    let mut output_data = vec![0_u64; n_words];

    // Figure out which words will be potentially shared between input slices in the concatenation,
    // and set those bits
    let mut last_word_indices = Vec::<usize>::new();
    let mut bits_so_far = 0_usize;
    for s in slices.iter() {
        // Copy the part that falls in the first output word
        let first_word = bits_so_far / 64;
        let first_word_bit_offset = bits_so_far % 64;
        let first_out = BitSlice::<u64, Lsb0>::from_element_mut(&mut output_data[first_word]);
        first_out[first_word_bit_offset..].copy_from_bitslice(&s[0..(64-first_word_bit_offset)]);


        // Copy the part that falls in the last output word (can be the same
        // as the first output word but that's okay).
        bits_so_far += s.len();
        if s.len() > 0 {
            let last_bit = bits_so_far - 1; // >= 0 because s.len() > 0
            let last_word = last_bit / 64;
            let last_word_bit_offset = last_bit % 64; // Index of the last bit that is written in the last word
            let n_bits_written_to_last_word = min(s.len(), last_word_bit_offset + 1);
            let last_out_word = BitSlice::<u64, Lsb0>::from_element_mut(&mut output_data[last_word]);
            let last_out_slice = &mut last_out_word[(1 + last_word_bit_offset - n_bits_written_to_last_word)..=last_word_bit_offset];
            let last_in_slice = &s[(s.len() - n_bits_written_to_last_word)..];
            last_out_slice.copy_from_bitslice(last_in_slice);

            last_word_indices.push(last_word);
        }

    }

    // Split output data into independent regions. Let w0, w1, w2... be the last
    // word indices. Split into [0..w0), (w0..w1), (w1..w2)..
    let mut remaining_output_slice = output_data.as_mut_slice();

    // Tuples of input bitslice, output u64 slice slice, bit index of first output bit
    let mut exclusive_output_regions = Vec::<(&BitSlice::<u64, Lsb0>, &mut[u64], usize)>::new();

    // TODO: build exclusive_output_regions here

    // Copy non-overlapping parts in parallel
    for (in_bitslice, out_wordslice, first_out_bit_idx) in exclusive_output_regions {
        //let out_bitslice = BitSlice::<u64, Lsb0>::from_slice_mut(out_wordslice);
        //out_bitslice[first_out_bit_idx..].copy_from_bitslice(in_bitslice);
    }

    todo!();

}

#[cfg(test)]
mod tests {

    use super::*;

    /*
    #[test]
    fn quick(){
        let v: Vec<usize> = vec![];
        let s = &v[1..1];
        eprintln!("{:?}", s);
    }
    */

    #[test]
    fn split_to_pieces() {
        // 7 one-bits -> ceil(7/3) = 3 per piece
        //              e              e              e     e
        let s = bitvec![0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1];
        let pieces = MergeInterleaving::split_to_pieces(&s, 3);
        assert_eq!(pieces.len(), 3);
        assert_eq!(pieces[0], (0, 0..5));
        assert_eq!(pieces[1], (2, 5..10));
        assert_eq!(pieces[2], (4, 10..12));
    }

    #[test]
    fn split_to_pieces_empty() {
        let s = bitvec![];
        let pieces = MergeInterleaving::split_to_pieces(&s, 3);
        assert_eq!(pieces.len(), 3);
        assert_eq!(pieces[0], (0, 0..0));
        assert_eq!(pieces[1], (0, 0..0));
        assert_eq!(pieces[2], (0, 0..0));
    }

    #[test]
    fn split_to_pieces_run_out_of_pieces() {
        //              0  1  2  3  4  5  6  7  8  9  10 11
        let s = bitvec![0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1];
        let pieces = MergeInterleaving::split_to_pieces(&s, 20);
        assert_eq!(pieces.len(), 20);
        assert_eq!(pieces[0], (0, 0..2));
        assert_eq!(pieces[1], (1, 2..3));
        assert_eq!(pieces[2], (1, 3..5));
        assert_eq!(pieces[3], (2, 5..8));
        assert_eq!(pieces[4], (4, 8..9));
        assert_eq!(pieces[5], (4, 9..10));
        assert_eq!(pieces[6], (4, 10..12));
        for i in 7..pieces.len(){
            assert_eq!(pieces[i], (5, 12..12));
        }
    }

    #[test]
    fn split_to_pieces_par() {
        // 7 one-bits -> ceil(7/3) = 3 per piece
        //              e              e              e     e
        let s = bitvec![0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1];
        let pieces = MergeInterleaving::split_to_pieces(&s, 3);
        let pieces_par = MergeInterleaving::split_to_pieces_par(&s, 3, 4);
        assert_eq!(pieces, pieces_par);

        let s = bitvec![];
        let pieces = MergeInterleaving::split_to_pieces(&s, 3);
        let pieces_par = MergeInterleaving::split_to_pieces_par(&s, 3, 4);
        assert_eq!(pieces, pieces_par);

        //              0  1  2  3  4  5  6  7  8  9  10 11
        let s = bitvec![0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1];
        let pieces = MergeInterleaving::split_to_pieces(&s, 20);
        let pieces_par = MergeInterleaving::split_to_pieces_par(&s, 20, 4);
        assert_eq!(pieces, pieces_par);
    }

}
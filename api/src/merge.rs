use std::cmp::min;
use std::ops::Range;

use bitvec::vec::BitVec;
use bitvec::prelude::*;
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

    pub fn new<SS: SubsetSeq + Send + Sync>(index1: &SbwtIndex::<SS>, index2: &SbwtIndex<SS>, n_threads: usize) -> MergeInterleaving {

        let k = index1.k();
        assert_eq!(k, index2.k());

        // We invert the SBWTs column by column and maintain ranges
        // in both SBWTs that are so far equal. The ranges are in increasing
        // order and the partition the SBWTs, to it's enough to just store their
        // sizes in order. The sizes are stored in concatenated unary representations.
        // Empty ranges are allowed.

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

            let new_arrays = Self::refine_segmentation(s1, s2, &chars1, &chars2, round == k-1);
            (s1, s2, leader_bits) = new_arrays;

            if round != k-1 {
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


    // Helper of a helper function
    fn refine_piece(s1: BitVec, s2: BitVec, chars1: &[u8], chars2: &[u8], mut c1_i: usize, mut c2_i: usize, s1_range: Range<usize>, s2_range: Range<usize>, compute_leaders: bool) 
    -> (BitVec, BitVec, Option<BitVec>) {

        let mut out1 = bitvec::bitvec![];
        let mut out2 = bitvec::bitvec![];

        let mut leader_bits = bitvec::bitvec![]; // Last round only

        // c1_i and c2_i are current indices in chars1 and chars2 respectively
        // s1_i and s2_i are current indices in s1 and s2 respectively
        let mut s1_i = 0_usize;
        let mut s2_i = 0_usize;

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

                while c1_i < c1_end && chars1[c1_i] == c {
                    c1_i += 1;
                    out1.push(false); // Unary bit
                }

                while c2_i < c2_end && chars2[c2_i] == c {
                    c2_i += 1;
                    out2.push(false); // Unary bit
                }

                // Terminate unary representations
                out1.push(true);
                out2.push(true);

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
    fn refine_segmentation(s1: BitVec, s2: BitVec, chars1: &[u8], chars2: &[u8], /*input_pieces: Vec<(usize, usize, Range<usize>, Range<usize>)>,*/ last_round: bool) -> (BitVec, BitVec, Option<BitVec>) {

        let s1_range = 0..s1.len();
        let s2_range = 0..s2.len();
        Self::refine_piece(s1, s2, chars1, chars2, 0, 0, s1_range, s2_range, last_round)

    }

}
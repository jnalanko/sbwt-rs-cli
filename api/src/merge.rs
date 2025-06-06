use std::cmp::min;
use std::ops::Range;

use bitvec::prelude::*;
use rayon::iter::IntoParallelIterator;
use rayon::iter::IntoParallelRefIterator;
use rayon::iter::ParallelIterator;
use crate::compact_int_vector::CompactIntVector;
use crate::subsetseq::*;
use crate::sbwt::*;

type BitVec = bitvec::vec::BitVec<u64, Lsb0>;
type BitSlice = bitvec::slice::BitSlice<u64, Lsb0>;

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

enum CharVector {
    ByteAlphabet(Vec<u8>),
    Compact(CompactIntVector<3>) // 3 bits per symbol for alphabet $,A,C,G,T
}

impl CharVector{
    fn len(&self) -> usize {
        match self {
            Self::ByteAlphabet(v) => v.len(),
            Self::Compact(v) => v.len(),
        }
    }

    fn init_with_byte_alphabet(len: usize) -> CharVector {
        CharVector::ByteAlphabet(vec![0; len])
    }

    fn init_compact(len: usize) -> CharVector {
        CharVector::Compact(CompactIntVector::<3>::new(len))
    }
}

impl MergeInterleaving {


    /// optimize_peak_ram enables optimizations to reduce the RAM peak at the expense of running time.
    pub fn new<SS: SubsetSeq + Send + Sync>(index1: &SbwtIndex::<SS>, index2: &SbwtIndex<SS>, optimize_peak_ram: bool, n_threads: usize) -> MergeInterleaving {

        use CharVector::*;

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
            let mut s1 = bitvec![u64, Lsb0; 0; index1.n_sets()];
            let mut s2 = bitvec![u64, Lsb0; 0; index2.n_sets()];
            s1.push(true);
            s2.push(true);

            let (mut chars1, mut chars2, mut temp_char_buf_1, mut temp_char_buf_2 ) = if optimize_peak_ram {
                (Compact(index1.build_last_column_compact()), 
                 Compact(index2.build_last_column_compact()),
                 None, // Temp buffers allocated on demand
                 None  // Temp buffers allocated on demand
                )
            } else {
                (ByteAlphabet(index1.build_last_column()), 
                 ByteAlphabet(index2.build_last_column()),
                 Some(CharVector::init_with_byte_alphabet(index1.n_sets())),
                 Some(CharVector::init_with_byte_alphabet(index2.n_sets()))
                )
            };

            for round in 0..k {
                log::info!("Round {}/{}", round+1, k);

                // Split work into pieces for different threads
                log::debug!("Splitting work");
                let p1 = split_to_pieces_par(&s1, n_threads, n_threads);
                let p2 = split_to_pieces_par(&s2, n_threads, n_threads);
                assert_eq!(p1.len(), n_threads);
                assert_eq!(p2.len(), n_threads);
                // Zip pairs of tuples into 4-tuples
                let pieces = (0..n_threads).map(|i| (p1[i].0, p2[i].0, p1[i].1.clone(), p2[i].1.clone())).collect();

                log::debug!("Refining segmentation");
                let new_arrays = refine_segmentation(s1, s2, &chars1, &chars2, pieces, round == k-1);
                (s1, s2, leader_bits) = new_arrays;

                if round != k-1 {
                    log::debug!("Pushing labels forward in the SBWT graph");
                    if let (ByteAlphabet(c1), ByteAlphabet(c2), Some(ByteAlphabet(temp1)), Some(ByteAlphabet(temp2))) = (&mut chars1, &mut chars2, &mut temp_char_buf_1, &mut temp_char_buf_2) {
                    if let (Some(buf1), Some(buf2)) = (&mut temp_char_buf_1, &mut temp_char_buf_2) {
                        // High memory mode
                        index1.push_all_labels_forward(&c1, temp1, n_threads);
                        std::mem::swap(&mut c1, &mut temp1);

                        index2.push_all_labels_forward(&c2, temp2, n_threads);
                        std::mem::swap(&mut c2, &mut temp2);
                    } else if let (Compact(c1), Compact(c2), None, None) = (&mut chars1, &mut chars2, &mut temp_char_buf_1, &mut temp_char_buf_2) {
                        // Low memory mode
                        let mut temp1 = CharVector::init_compact(chars1.len()) ;
                        index1.push_all_labels_forward_compact(&chars1, &mut buf1, n_threads);
                        std::mem::swap(&mut chars1, &mut buf1);
                        drop(buf1);

                        let mut temp2 = CharVector::init_compact(chars2.len()) ;
                        index2.push_all_labels_forward(&chars2, &mut buf2, n_threads);
                        std::mem::swap(&mut chars2, &mut buf2);
                        drop(buf2)
                    } else {
                        panic!("Programmer messed up");
                    }
                }
            }

            drop(temp_char_buf_1);
            drop(temp_char_buf_2);

            let leader_bits = leader_bits.unwrap(); // Computed in the last round
            log::debug!("Number of suffix groups: {}", leader_bits.count_ones());

            assert_eq!(s1.len(), s2.len());
            let merged_len = s1.len();

            // Identify dummies in the merged SBWT

            log::debug!("Marking dummy nodes");
            let piece_len = merged_len.div_ceil(n_threads);
            let merged_piece_ranges: Vec<Range<usize>> = (0..n_threads).map(|t| t*piece_len..min((t+1)*piece_len, merged_len)).collect();
            let s1_piece_popcounts: Vec<usize> = merged_piece_ranges.par_iter().map(|range| s1[range.clone()].count_ones()).collect();
            let s2_piece_popcounts: Vec<usize> = merged_piece_ranges.par_iter().map(|range| s2[range.clone()].count_ones()).collect();
            let is_dummy_pieces = (0..n_threads).into_par_iter().map(|thread_idx| {
                let colex_range = &merged_piece_ranges[thread_idx];
                let mut is_dummy_bits = bitvec![u64, Lsb0 ;];
                is_dummy_bits.resize(colex_range.len(), false);

                let mut c1_idx: usize = s1_piece_popcounts[..thread_idx].iter().sum(); // Skip over previous pieces
                let mut c2_idx: usize = s2_piece_popcounts[..thread_idx].iter().sum(); // Skip over previous pieces
                for colex in colex_range.clone() {
                    let d1 = s1[colex] && c1_idx < chars1.len() && chars1[c1_idx] == b'$';
                    let d2 = s2[colex] && c2_idx < chars2.len() && chars2[c2_idx] == b'$';

                    let rel_colex = colex - colex_range.start;
                    is_dummy_bits.set(rel_colex, d1 || d2); 

                    c1_idx += s1[colex] as usize;
                    c2_idx += s2[colex] as usize;
                }
                is_dummy_bits 
            }).collect();

            let is_dummy = crate::util::parallel_bitvec_concat(is_dummy_pieces);

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
    if seq.is_empty() {
        return 0;
    }
    if seq.len() > 200 {
        // Binary search the first element that is larger than c
        seq.binary_search_by(|&x| {
            if x > c {
                std::cmp::Ordering::Greater
            } else {
                std::cmp::Ordering::Less
            }
        }).unwrap_err() // Is always Err because we never return Ordering::Equal
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


// Helper of a helper function.
// Yeah I know it's got a lot of arguments but it's an internal helper and collecting the arguments
// into bigger structs would not really make it much clearer and would just increase the number of lines of code. 
#[allow(clippy::too_many_arguments)] 
fn refine_piece(s1: &BitVec, s2: &BitVec, chars1: &[u8], chars2: &[u8], mut c1_i: usize, mut c2_i: usize, s1_range: Range<usize>, s2_range: Range<usize>, last_round: bool) 
-> (BitVec, BitVec, Option<BitVec>) {

    let mut out1 = bitvec::bitvec![u64, Lsb0;];
    let mut out2 = bitvec::bitvec![u64, Lsb0;];

    let mut leader_bits = bitvec::bitvec![u64, Lsb0;]; // Last round only

    // c1_i and c2_i are current indices in chars1 and chars2 respectively
    // s1_i and s2_i are current indices in s1 and s2 respectively
    let mut s1_i = s1_range.start;
    let mut s2_i = s2_range.start;

    while s1_i < s1_range.end {

        assert!(s2_i < s2_range.end);

        let len1 = leading_zeros(&s1[s1_i..]);
        let len2 = leading_zeros(&s2[s2_i..]);

        let c1_end = c1_i + len1; // One past the end
        let c2_end = c2_i + len2; // One past the end

        let mut is_leader = true;
        while c1_i < c1_end || c2_i < c2_end {
            let c1 = if c1_i == c1_end { u8::MAX } else { chars1[c1_i] };
            let c2 = if c2_i == c2_end { u8::MAX } else { chars2[c2_i] };
            let c = min(c1,c2);

            let r1 = run_length_in_sorted_seq(&chars1[c1_i..c1_end], c);
            let r2 = run_length_in_sorted_seq(&chars2[c2_i..c2_end], c);

            if last_round {
                // We know that r1 and r2 are at most 1 -> no need to have a full unary code.
                // Let's not assert this here because this is a tight inner loop.
                out1.push(r1 > 0);
                out2.push(r2 > 0);
            } else {
                // Write r1 and r2 in unary
                zero_extend(&mut out1, r1);
                zero_extend(&mut out2, r2);
                // Terminate unary representations
                out1.push(true);
                out2.push(true);
            }

            // Advance indexes in chars
            c1_i += r1;
            c2_i += r2;

            if last_round {
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

    (out1, out2, if last_round { Some(leader_bits) } else { None })
}

// Helper function for construction.
// Input pieces are pairs (start in chars1, start in chars2, range in s1, range in s2)
// Input piece are for parallelism: one piece to work on for each thread.
// Returns the new segmentation in concatenate unary form.
// One the last round the output is different: now the two output bit vectors would have
// only 0 and 1 encoded in unary ("1" and "01"). Instead we write just the bits 0 and 1.
// On the also round the function also returns the leader bit vector, which marks
// the smallest k-mer in each group of k-mers with the same suffix of length (k-1).
// This function runs in parallel, so a rayon thread pool must be initialized.
fn refine_segmentation(s1: BitVec, s2: BitVec, chars1: &[u8], chars2: &[u8], input_pieces: Vec<(usize, usize, Range<usize>, Range<usize>)>, last_round: bool) -> (BitVec, BitVec, Option<BitVec>) {
    let output_pieces: Vec<(BitVec, BitVec, Option<BitVec>)> = input_pieces.par_iter().map(|piece| {
        let (c1_i, c2_1, s1_range, s2_range) = piece;
        refine_piece(&s1, &s2, chars1, chars2, *c1_i, *c2_1, s1_range.clone(), s2_range.clone(), last_round)
    }).collect();

    // Free memory
    drop(s1);
    drop(s2);

    // Turn output pieces (a vec of triples) into triple of vecs 
    let mut new_s1_pieces = vec![];
    let mut new_s2_pieces = vec![];
    let mut leader_pieces = if last_round {Some(vec![])} else {None};
    for (a,b,c) in output_pieces.into_iter() {
        new_s1_pieces.push(a);
        new_s2_pieces.push(b);
        if last_round {
            leader_pieces.as_mut().unwrap().push(c.unwrap());
        }
    }
    log::debug!("Concatenating pieces");
    let new_s1 = crate::util::parallel_bitvec_concat(new_s1_pieces);
    let new_s2 = crate::util::parallel_bitvec_concat(new_s2_pieces);
    let new_leader_bits = leader_pieces.map(crate::util::parallel_bitvec_concat);

    (new_s1, new_s2, new_leader_bits)
}

// Parallel version of split_to_pieces (see mod tests for a single-threaded reference implementation)
// Returns a segmentation s[l_1..r_1), s[l_2..r_2], ... such that
// each segment ends in a 1-bit and has an approximately equal number of 1-bits.
// Also returns the number of 0-bits before each segment.
fn split_to_pieces_par(s: &BitSlice, n_pieces: usize, n_threads: usize) -> Vec<(usize, Range<usize>)> {
    const BLOCK_SIZE: usize = 1024;

    let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();
    thread_pool.install(||{
        // Strategy: compute block popcounts in parallel, then locate the blocks where the
        // pieces start, and do bit-by-bit counting to find the precise start points of pieces.

        assert!(n_pieces > 0);
        if !s.is_empty() {
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
        assert!(!starts.is_empty());
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



#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn split_to_pieces() {
        // 7 one-bits -> ceil(7/3) = 3 per piece
        //              e              e              e     e
        let s = bitvec![u64, Lsb0; 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1];
        let pieces = split_to_pieces_ref_implementation(&s, 3);
        assert_eq!(pieces.len(), 3);
        assert_eq!(pieces[0], (0, 0..5));
        assert_eq!(pieces[1], (2, 5..10));
        assert_eq!(pieces[2], (4, 10..12));
    }

    #[test]
    fn split_to_pieces_empty() {
        let s = bitvec![u64, Lsb0;];
        let pieces = split_to_pieces_ref_implementation(&s, 3);
        assert_eq!(pieces.len(), 3);
        assert_eq!(pieces[0], (0, 0..0));
        assert_eq!(pieces[1], (0, 0..0));
        assert_eq!(pieces[2], (0, 0..0));
    }

    #[test]
    fn split_to_pieces_run_out_of_pieces() {
        //              0  1  2  3  4  5  6  7  8  9  10 11
        let s = bitvec![u64, Lsb0; 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1];
        let pieces = split_to_pieces_ref_implementation(&s, 20);
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
    fn test_split_to_pieces_par() {
        // 7 one-bits -> ceil(7/3) = 3 per piece
        //              e              e              e     e
        let s = bitvec![u64, Lsb0; 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1];
        let pieces = split_to_pieces_ref_implementation(&s, 3);
        let pieces_par = split_to_pieces_par(&s, 3, 4);
        assert_eq!(pieces, pieces_par);

        let s = bitvec![u64, Lsb0;];
        let pieces = split_to_pieces_ref_implementation(&s, 3);
        let pieces_par = split_to_pieces_par(&s, 3, 4);
        assert_eq!(pieces, pieces_par);

        //              0  1  2  3  4  5  6  7  8  9  10 11
        let s = bitvec![u64, Lsb0; 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1];
        let pieces = split_to_pieces_ref_implementation(&s, 20);
        let pieces_par = split_to_pieces_par(&s, 20, 4);
        assert_eq!(pieces, pieces_par);
    }

    // This is not a test. This is a single-threaded reference implementation.
    // Returns a segmentation s[l_1..r_1), s[l_2..r_2], ... such that
    // each segment ends in a 1-bit and has an approximately equal number of 1-bits.
    // Also returns the number of 0-bits before each segment.
    fn split_to_pieces_ref_implementation(s: &BitSlice, n_pieces: usize) -> Vec<(usize, Range<usize>)> {
        assert!(n_pieces > 0);
        if !s.is_empty() {
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


}
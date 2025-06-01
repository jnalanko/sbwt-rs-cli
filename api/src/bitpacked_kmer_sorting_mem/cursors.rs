use bitvec::order::Lsb0;
use rand::AsByteSliceMut;
use simple_sds_sbwt::ops::Push;
use simple_sds_sbwt::serialize::Serialize;
use std::cmp::min;
use std::io::Cursor;
use std::io::Read;
use std::ops::Range;
use crate::kmer::LongKmer;
use crate::util::binary_search_leftmost_that_fulfills_pred;

use rayon::iter::IntoParallelIterator;
use rayon::iter::ParallelIterator;

use super::dummies::KmersWithLengths;

pub struct KmerDummyMergeSlice<'a, const B: usize> {
    all_dummies: &'a KmersWithLengths<B>,
    all_kmers: &'a [LongKmer<B>],
    dummy_idx: usize,
    kmer_idx: usize,
    dummy_range: Range<usize>,
    kmer_range: Range<usize>,
    k: usize,
}

impl<'a, const B:usize> KmerDummyMergeSlice<'a, B> {

    pub fn len(&self) -> usize {
        self.dummy_range.len() + self.kmer_range.len()
    }

    pub fn cur_merged_index(&self) -> usize {
        self.dummy_idx + self.kmer_idx
    }

    fn advance_dummy(&mut self) -> (LongKmer<B>, u8) {
        let x = self.all_dummies.get(self.dummy_idx);
        self.dummy_idx += 1;
        x
    }

    fn advance_kmers(&mut self) -> (LongKmer<B>, u8) {
        let x = self.all_kmers[self.kmer_idx];
        self.kmer_idx += 1;
        (x, self.k as u8)
    }

    pub fn peek(&mut self) -> Option<(LongKmer<B>, u8)> {
        // Todo: refactor this logic together with next()
        if self.dummy_idx == self.dummy_range.end && self.kmer_idx == self.kmer_range.end {
            None
        } else if self.dummy_idx == self.dummy_range.end {
            Some((self.all_kmers[self.kmer_idx], self.k as u8))
        } else if self.kmer_idx == self.kmer_range.end {
            Some(self.all_dummies.get(self.dummy_idx))
        } else if (self.all_kmers[self.kmer_idx], self.k as u8) < self.all_dummies.get(self.dummy_idx) {
            Some((self.all_kmers[self.kmer_idx], self.k as u8))
        } else {
            Some(self.all_dummies.get(self.dummy_idx))
        }
    }

    pub fn next(&mut self) -> Option<(LongKmer<B>, u8)> {
        if self.dummy_idx == self.dummy_range.end && self.kmer_idx == self.kmer_range.end {
            None
        } else if self.dummy_idx == self.dummy_range.end {
            Some(self.advance_kmers())
        } else if self.kmer_idx == self.kmer_range.end {
            Some(self.advance_dummy())
        } else if (self.all_kmers[self.kmer_idx], self.k as u8) < self.all_dummies.get(self.dummy_idx) {
            Some(self.advance_kmers())
        } else {
            Some(self.advance_dummy())
        }
    }

    pub fn new(all_dummies: &'a KmersWithLengths<B>, all_kmers: &'a [LongKmer<B>], merged_range: Range<usize>, k: usize) -> Self {
        assert!(merged_range.end <= all_dummies.len() + all_kmers.len());

        // Binary search the starts and ends in kmers and dummies.

        let (kmer_start, dummy_start, _) = binary_search_position_in_merged_list(
            |i| (all_kmers[i], k as u8), 
            |j| all_dummies.get(j), 
            merged_range.start, 
            all_kmers.len(), 
            all_dummies.len()
        );

        let (kmer_end, dummy_end, _) = binary_search_position_in_merged_list(
            |i| (all_kmers[i], k as u8), 
            |j| all_dummies.get(j), 
            merged_range.end, 
            all_kmers.len(), 
            all_dummies.len()
        );

        Self {
            all_dummies, all_kmers, 
            dummy_idx: dummy_start,
            kmer_idx: kmer_start,
            dummy_range: dummy_start..dummy_end,
            kmer_range: kmer_start..kmer_end,
            k
        }
    }
}

// log^2(n) time complexity 
pub fn get_ith_merged_kmer<const B: usize>(all_kmers: &Vec<LongKmer<B>>, all_dummies: &KmersWithLengths<B>, i: usize, k: usize) -> (LongKmer<B>, u8) {
    let (kmers_idx, dummy_idx, in_kmers) = binary_search_position_in_merged_list(
        |i| (all_kmers[i], k as u8), 
        |j| all_dummies.get(j), 
        i, all_kmers.len(), all_dummies.len());

    if in_kmers {
        (all_kmers[kmers_idx], k as u8)
    } else {
        all_dummies.get(dummy_idx)
    }
}

// Assumes the two sorted lists have no duplicates in them or between them.
// Returns (pos_a, pos_b, take_from_a)
// Takes O(log^2(n)) time
pub fn binary_search_position_in_merged_list<T: PartialOrd + Eq, Access1: Fn(usize) -> T, Access2: Fn(usize) -> T>(access_a: Access1, access_b: Access2, target_pos: usize, len_a: usize, len_b: usize) -> (usize, usize, bool) {

    assert!(target_pos <= len_a + len_b); // One-past the end allowed 

    let pos_in_merged_list = |a_idx: usize| {
        // Returns the index of a[a_idx] in the merged list of a and b
        // That is equal to to number of element in a that are smaller than a[a_idx],
        // and the number of elements in b that are smaller than a[a_idx].
        // We imagine that a[a_len] is infinity to make that corner case work.

        if a_idx == len_a {
            // All of a and b are smaller than infinity
            len_a + len_b 
        } else {
            let b_count = binary_search_leftmost_that_fulfills_pred(|j| j, |b_idx| access_b(b_idx) > access_a(a_idx), len_b);
            a_idx + b_count
        }

    };

    // Find the smallest a_idx that is at or after the target position
    let a_idx = binary_search_leftmost_that_fulfills_pred(|i| i, |a_idx| pos_in_merged_list(a_idx) >= target_pos, len_a);
    if pos_in_merged_list(a_idx) == target_pos {
        // Bingo. Find the b_idx of the next element in b, that is
        // the smallest element in b that is larger than a[a_idx]
        if a_idx == len_a {
            (len_a, len_b, true) // End
        } else {
            let x = access_a(a_idx); 
            let b_idx = binary_search_leftmost_that_fulfills_pred(|j| j, |b_idx| access_b(b_idx) > x, len_b);
            assert_eq!(a_idx + b_idx, target_pos);
            (a_idx, b_idx, true)
        }
    } else {
        // a[a_idx] is after the target position, and a[a_idx-1] is before it (if exists)
        // This means that there are a_idx elements from a before the target position.
        (a_idx, target_pos - a_idx, false)
    }
}

pub fn find_first_starting_with<const B: usize>(
    v: &[LongKmer<B>],
    c: u8,
) -> usize {
    binary_search_leftmost_that_fulfills_pred(
        |pos| { v[pos] },
        |kmer: LongKmer::<B>| { kmer.get_from_left(0) >= c },
        v.len()
    )
}


// Returns the LCS array.
// This first builds the uncompressed LCS array of &[u16], and then compressed it to log(k) bits per element.
// This means that the peak memory is high, but we have the k-mers in memory anyway, so it's not so bad
// TODO: compute the compressed LCS directly. 
pub fn build_lcs_array<const B: usize>(
    mut merged_slice: KmerDummyMergeSlice<B>,
    k: usize,
) -> simple_sds_sbwt::int_vector::IntVector {
    // LCS values are between 0 and k-1
    assert!(k > 0);
    assert!(k < u16::MAX as usize);
    //assert!(kmers.len() > 0);

    let bitwidth = 64 - (k as u64 - 1).leading_zeros();

    let mut non_compressed_lcs = Vec::<u16>::with_capacity(merged_slice.len());
    non_compressed_lcs.push(0);
    let (mut prev_kmer, mut prev_len) = merged_slice.next().unwrap();

    while let Some((kmer, len)) = merged_slice.next() {
        let lcp_value = LongKmer::<B>::lcp_with_different_lengths((&prev_kmer, prev_len), (&kmer, len));
        non_compressed_lcs.push(lcp_value as u16);
        (prev_kmer, prev_len) = (kmer,len);
    }

    // Compress into log(k) bits per element
    let mut compressed_lcs = simple_sds_sbwt::int_vector::IntVector::with_capacity(merged_slice.len(), bitwidth as usize).unwrap();
    for x in non_compressed_lcs {
        compressed_lcs.push(x as u64);
    }
    compressed_lcs
}

/// Assumes both input vector have no duplicates, and any
/// pair (kmer, len) exists in only one of the inputs.
pub fn merge_kmers_and_dummies<const B: usize>(
    mut kmers: Vec<LongKmer<B>>,
    dummies: Vec<(LongKmer<B>, u8)>,
    k: usize) -> KmersWithLengths<B> {

    let n_merged = kmers.len() + dummies.len();
    let n_non_dummies = kmers.len();

    // Merge from back to front, putting the merged kmers to the same vector as the input kmers.
    // Make more space
    kmers.resize(n_merged, LongKmer::from_u64_data([0; B]));

    // Allocate space for the lengths
    let mut lengths: Vec<u8> = vec![0; n_merged];

    let mut nondummy_in = n_non_dummies as isize - 1; // Input pointer in non-dummies
    let mut dummy_in = dummies.len() as isize - 1; // Input pointer in dummies
    let mut merged_out = n_merged as isize - 1; // Output pointer

    // The merge algorithm is easy because there are no duplicates
    while merged_out >= 0 {
        if dummy_in < 0 || (nondummy_in >= 0 && (kmers[nondummy_in as usize], k as u8) > dummies[dummy_in as usize]) {
            // Take an element from nondummies
            kmers[merged_out as usize] = kmers[nondummy_in as usize];
            lengths[merged_out as usize] = k as u8;
            merged_out -= 1;
            nondummy_in -= 1;
        } else {
            // Take an element from dummies
            kmers[merged_out as usize] = dummies[dummy_in as usize].0;
            lengths[merged_out as usize] = dummies[dummy_in as usize].1;
            merged_out -= 1;
            dummy_in -= 1;
        }
    }
    KmersWithLengths{kmers, lengths}
}

// Here c is from 0..3
fn prepend_c<const B: usize>(kmer: (LongKmer<B>, u8), k: usize, c: usize) -> (LongKmer<B>, u8) {
    if kmer.1 as usize == k {
        (
            kmer.0
                .copy_set_from_left(k - 1, 0)
                .right_shifted(1)
                .copy_set_from_left(0, c as u8),
            k as u8,
        )
    } else {
        (kmer.0.right_shifted(1).copy_set_from_left(0, c as u8), kmer.1 + 1) // Dummy
    }
}

// Returns the SBWT bit vectors and optionally the LCS array
pub fn build_sbwt_bit_vectors<const B: usize>(
    kmers: Vec<LongKmer<B>>,
    dummies: KmersWithLengths<B>,
    k: usize,
    sigma: usize,
    build_lcs: bool,
    n_threads: usize,
) -> (Vec<simple_sds_sbwt::raw_vector::RawVector>, Option<simple_sds_sbwt::int_vector::IntVector>) {

    let n = kmers.len() + dummies.len(); // Merged lengths

    let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();
    thread_pool.install(||{
        // Split the merged range 0..n into segments for threads
        let segment_len = n.div_ceil(n_threads);
        let mut input_ranges: Vec<Range<usize>> = vec![];
        for t in 0..n_threads {
            input_ranges.push(t*segment_len .. min((t+1)*segment_len,n)); // Final segments may be empty. Is ok.
        }

        let mut rows = vec![];
        for c in 0..sigma {
            let row_pieces = input_ranges.clone().into_par_iter().map(|input_range|{
                let mut row_piece: bitvec::vec::BitVec::<u64, Lsb0> = bitvec::vec::BitVec::with_capacity(input_range.len());
                row_piece.resize(input_range.len(), false);

                if !input_range.is_empty() {
                    // Let x be the first input kmer in this range. The destinations of the c-edges start
                    // from the smallest k-mer that is larger or equal to cx.
                    let x: (LongKmer<B>, u8) = get_ith_merged_kmer(&kmers, &dummies, input_range.start, k);
                    let (cx, cx_len) = prepend_c(x, k, c);
                    let cx_kmer_insertion_index = binary_search_leftmost_that_fulfills_pred(|i| kmers[i], |y| (y, k as u8) >= (cx, cx_len), kmers.len()); 
                    let cx_dummy_insertion_index = binary_search_leftmost_that_fulfills_pred(|i| dummies.get(i), |y| y >= (cx, cx_len), dummies.len());
                    let dest_start_idx = cx_kmer_insertion_index + cx_dummy_insertion_index;

                    let mut src_pointer = KmerDummyMergeSlice::new(&dummies, &kmers, input_range.clone(), k); // Origin of edge
                    let mut dest_pointer = KmerDummyMergeSlice::new(&dummies, &kmers, dest_start_idx..n, k); // Destination of edge

                    // We might be starting at a k-mer that is not a suffix group leader. We must not
                    // add any edges for it. Rewind forward until we are at a suffix group leader.
                    let mut cur = x;
                    if input_range.start > 0 {
                        let mut prev = get_ith_merged_kmer(&kmers, &dummies, input_range.start - 1, k);
                        while LongKmer::<B>::lcp_with_different_lengths((&prev.0, prev.1), (&cur.0, cur.1)) == k-1 {
                            src_pointer.next(); // Advance
                            prev = cur;
                            if let Some(next) = src_pointer.peek() {
                                cur = next;
                            } else {
                                break;
                            }
                        }
                    }

                    // Now src_pointer should point to the first leader of a suffix group
                    // in the input range.
                    while let Some((kmer, len)) = src_pointer.next() {
                        let kmer_c = prepend_c((kmer,len), k, c);

                        while dest_pointer.peek().is_some_and(|y| y < kmer_c) {
                            dest_pointer.next();
                        }

                        if dest_pointer.peek().is_some_and(|y| y == kmer_c) {
                            row_piece.set(src_pointer.cur_merged_index() - 1 - input_range.start, true); // -1 because we have advanced past the current k-mer
                            dest_pointer.next().unwrap(); // Don't point to this k-mer anymore because it now has an edge
                        }
                    };
                }
                row_piece
            }).collect();
            rows.push(crate::util::parallel_bitvec_concat(row_pieces));
        }


        // Convert rows into simple_sds_sbwt vectors
        let rows: Vec<simple_sds_sbwt::raw_vector::RawVector> = rows.into_par_iter().map(|mut row| {
            // Let's use the deserialization function in simple_sds_sbwt for a raw bitvector.
            // It requires the following header: (TODO: this code is repeated in sbwt.rs... refactor)
            let mut header = [0u64, 0u64]; // bits, words
            header[0] = row.len() as u64; // Assumes little-endian byte order
            header[1] = row.len().div_ceil(64) as u64;

            let header_bytes = header.as_byte_slice_mut();
            let raw_data = row.as_raw_mut_slice().as_byte_slice_mut();
            let mut data_with_header = Cursor::new(header_bytes).chain(Cursor::new(raw_data));

            simple_sds_sbwt::raw_vector::RawVector::load(&mut data_with_header).unwrap()
        }).collect();

        let lcs = if build_lcs {
            // LCS values are between 0 and k-1
            Some(build_lcs_array(KmerDummyMergeSlice::new(&dummies, &kmers, 0..n, k), k))
        } else {
            None
        };

        (rows, lcs)
    })
}


#[cfg(test)]
mod tests {
    use rand::{seq::SliceRandom, Rng, SeedableRng};

    #[test]
    fn test_binary_search_merged_list() {
        // Test empty vs empty
        assert_eq!((0,0,true), super::binary_search_position_in_merged_list(|i| i, |j| j, 0, 0, 0));

        // Random testcases
        let seed = 1234;
        let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
        for rep in 0..100 { // Stress test with 100 random runs
            let mut v: Vec<usize> = (0..20).collect();
            v.shuffle(&mut rng);
            let split_point = if rep == 0 {
                0 // Test empty v1    
            } else if rep == 1 {
                v.len() // Test empty v2
            } else {
                // Random split point
                rng.gen_range(0,v.len()+1)
            };

            let mut v1 = v[0..split_point].to_vec();
            let mut v2 = v[split_point..].to_vec();
            v1.sort();
            v2.sort();

            let mut merged: Vec<(usize, usize, bool)> = vec![]; // (i_v1, i_v2, b). b tells which vector it's from
            merged.extend(v1.iter().enumerate().map(|x| (*x.1, x.0, true)));
            merged.extend(v2.iter().enumerate().map(|x| (*x.1, x.0, false)));
            merged.sort();
            let n_merged = merged.len();
            merged.push((v1.len(), v2.len(), true)); // One past the end. The bit is true to match how the search treats this.

            for query in 0..=n_merged {
                let mut true_i = 0;
                let mut true_j = 0;
                for (_,_,from_a) in &merged[0..query] {
                    true_i += *from_a as usize;
                    true_j += !(*from_a) as usize;
                }
                let true_from_a = merged[query].2;
                let (i,j,from_a) = super::binary_search_position_in_merged_list(|i| v1[i], |j| v2[j], query, v1.len(), v2.len());
                assert_eq!((i,j,from_a), (true_i, true_j, true_from_a));
            }
        }
    }
}
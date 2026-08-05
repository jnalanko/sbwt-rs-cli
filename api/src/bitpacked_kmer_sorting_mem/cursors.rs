use bitvec::order::Lsb0;
use simple_sds_sbwt::ops::Push;
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
enum NextInfo<const B: usize> {
    None,
    Kmer((LongKmer<B>, u8)),
    Dummy((LongKmer<B>, u8)),
}

impl<'a, const B:usize> KmerDummyMergeSlice<'a, B> {

    pub fn len(&self) -> usize {
        self.dummy_range.len() + self.kmer_range.len()
    }

    pub fn cur_merged_index(&self) -> usize {
        self.dummy_idx + self.kmer_idx
    }

    fn determine_next(&self) -> NextInfo<B> {
        let dummies_done = self.dummy_idx == self.dummy_range.end;
        let kmers_done = self.kmer_idx == self.kmer_range.end;
        if kmers_done && dummies_done {
            NextInfo::None
        } else if dummies_done || (!kmers_done && (self.all_kmers[self.kmer_idx], self.k as u8) < self.all_dummies.get(self.dummy_idx)) {
            NextInfo::Kmer((self.all_kmers[self.kmer_idx], self.k as u8))
        } else { 
            NextInfo::Dummy(self.all_dummies.get(self.dummy_idx))
        }

    }

    pub fn peek(&mut self) -> Option<(LongKmer<B>, u8)> {
        match self.determine_next() {
            NextInfo::None => None,
            NextInfo::Kmer(x) => Some(x),
            NextInfo::Dummy(x) => Some(x),
        }
    }

    pub fn next(&mut self) -> Option<(LongKmer<B>, u8)> {
        match self.determine_next() {
            NextInfo::None => None,
            NextInfo::Kmer(x) => {
                self.kmer_idx += 1;
                Some(x)
            },
            NextInfo::Dummy(x) => {
                self.dummy_idx += 1;
                Some(x)
            }
        }
    }

    pub fn new(all_dummies: &'a KmersWithLengths<B>, all_kmers: &'a [LongKmer<B>], merged_range: Range<usize>, k: usize) -> Self {
        assert!(merged_range.end <= all_dummies.len() + all_kmers.len());

        // Binary search the starts and ends in kmers and dummies.

        let (kmer_start, dummy_start, _) = crate::util::binary_search_position_in_merged_list(
            |i| (all_kmers[i], k as u8),
            |j| all_dummies.get(j),
            merged_range.start,
            all_kmers.len(),
            all_dummies.len()
        );

        let (kmer_end, dummy_end, _) = crate::util::binary_search_position_in_merged_list(
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
pub fn get_ith_merged_kmer<const B: usize>(all_kmers: &[LongKmer<B>], all_dummies: &KmersWithLengths<B>, i: usize, k: usize) -> (LongKmer<B>, u8) {
    let (kmers_idx, dummy_idx, in_kmers) = crate::util::binary_search_position_in_merged_list(
        |i| (all_kmers[i], k as u8),
        |j| all_dummies.get(j),
        i, all_kmers.len(), all_dummies.len());

    if in_kmers {
        (all_kmers[kmers_idx], k as u8)
    } else {
        all_dummies.get(dummy_idx)
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
    kmers: &[LongKmer<B>],
    dummies: &KmersWithLengths<B>,
    k: usize, n_threads: usize
) -> simple_sds_sbwt::int_vector::IntVector {
    // LCS values are between 0 and k-1
    assert!(k > 0);
    assert!(k < u16::MAX as usize);
    //assert!(kmers.len() > 0);

    log::info!("Computing LCS values");

    let full_slice = KmerDummyMergeSlice::new(dummies, kmers, 0..(kmers.len()+dummies.len()), k);

    // We start the segmentation from 1 so that we always have a previous k-mer to compare against
    let segments = crate::util::segment_range(1..full_slice.len(), n_threads);
    let lcs_pieces: Vec<Vec<u16>> = segments.into_par_iter().map(|range| {
        let (mut prev_kmer, mut prev_len) = get_ith_merged_kmer(kmers, dummies, range.start - 1, k); // range.start >= 1 so this is ok
        let mut subslice = KmerDummyMergeSlice::new(dummies, kmers, range.clone(), k);
        let mut lcs_piece = Vec::<u16>::with_capacity(range.len());
        while let Some((kmer, len)) = subslice.next() {
            let lcp_value = LongKmer::<B>::lcp_with_different_lengths((&prev_kmer, prev_len), (&kmer, len));
            lcs_piece.push(lcp_value as u16);
            (prev_kmer, prev_len) = (kmer,len);
        }
        lcs_piece
    }).collect();


    // Compress into log(k) bits per element
    log::info!("Compressing LCS array to log(k) bits per element");
    let bitwidth = 64 - (k as u64 - 1).leading_zeros();
    let mut compressed_lcs = simple_sds_sbwt::int_vector::IntVector::with_capacity(full_slice.len(), bitwidth as usize).unwrap();
    compressed_lcs.push(0); // lcs[0] = 0 by definition
    for piece in lcs_pieces { // Concatenate pieces. Todo: could this be done in parallel?
        compressed_lcs.extend(piece);
    }

    compressed_lcs
}

/// Assumes both input vector have no duplicates, and any
/// pair (kmer, len) exists in only one of the inputs.
#[allow(dead_code)] 
// Dead code: Could be useful later because if we can spend some memory to do the concenation
// then we might it might speed up the SBWT rows construction because then the merge logic happens only
// once here and then it's easy to scan for the SBWT rows.
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
) -> (Vec<bitvec::vec::BitVec::<u64, Lsb0>>, Option<simple_sds_sbwt::int_vector::IntVector>) {

    let n = kmers.len() + dummies.len(); // Merged lengths

    let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();
    thread_pool.install(||{
        // Split the merged range 0..n into segments for threads
        let input_ranges = crate::util::segment_range(0..n, n_threads);

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

        let lcs = if build_lcs {
            // LCS values are between 0 and k-1
            Some(build_lcs_array(&kmers, &dummies, k, n_threads))
        } else {
            None
        };

        (rows, lcs)
    })
}
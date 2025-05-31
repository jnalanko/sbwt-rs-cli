use bitvec::order::Lsb0;
use bitvec::prelude;
use rand::AsByteSliceMut;
use simple_sds_sbwt::ops::Push;
use simple_sds_sbwt::raw_vector::AccessRaw;
use simple_sds_sbwt::serialize::Serialize;
use std::cmp::min;
use std::io::Cursor;
use std::io::Read;
use std::ops::Range;
use crate::kmer::Kmer;
use crate::kmer::LongKmer;
use crate::util::binary_search_leftmost_that_fulfills_pred;
use crate::util::parallel_bitvec_concat;

use rayon::iter::IntoParallelIterator;
use rayon::iter::ParallelIterator;
use rayon::prelude::ParallelSlice;

use super::dummies::get_ith_merged_kmer;
use super::dummies::KmerDummyMergeSlice;
use super::dummies::KmersWithLengths;


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

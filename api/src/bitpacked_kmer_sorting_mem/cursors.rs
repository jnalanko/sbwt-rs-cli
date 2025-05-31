use simple_sds_sbwt::ops::Push;
use simple_sds_sbwt::raw_vector::AccessRaw;
use std::cmp::min;
use crate::kmer::Kmer;
use crate::kmer::LongKmer;

use rayon::iter::IntoParallelIterator;
use rayon::iter::ParallelIterator;
use rayon::prelude::ParallelSlice;

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
        let mut lcp_value = LongKmer::<B>::lcp(&prev_kmer, &kmer);
        lcp_value = min(lcp_value, min(prev_len as usize, len as usize));
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

// Returns the SBWT bit vectors and optionally the LCS array
pub fn build_sbwt_bit_vectors<const B: usize>(
    kmers: Vec<LongKmer<B>>,
    dummies: KmersWithLengths<B>,
    k: usize,
    sigma: usize,
    build_lcs: bool,
) -> (Vec<simple_sds_sbwt::raw_vector::RawVector>, Option<simple_sds_sbwt::int_vector::IntVector>) {

    let n = kmers.len() + dummies.len(); // Merged lengths

    let rawrows = (0..sigma).collect::<Vec<usize>>().into_iter().map(|c|{
        let mut rawrow = simple_sds_sbwt::raw_vector::RawVector::with_len(n, false);
        let mut src_pointer = KmerDummyMergeSlice::new(&dummies, &kmers, 0..n, k); // Origin of edge
        let mut dest_pointer = KmerDummyMergeSlice::new(&dummies, &kmers, 0..n, k); // Destination of edge
        //let mut pointed_idx = merged.first_that_starts_with(c as u8); // TODO: from C array
        //let end = merged.first_that_starts_with(c as u8 + 1); // TODO: from C array

        while let Some((kmer, len)) = src_pointer.next() {
            let kmer_c = if len as usize == k {
                (
                    kmer
                        .copy_set_from_left(k - 1, 0)
                        .right_shifted(1)
                        .copy_set_from_left(0, c as u8),
                    k as u8,
                )
            } else {
                (kmer.right_shifted(1).copy_set_from_left(0, c as u8), len + 1) // Dummy
            };

            while dest_pointer.peek().is_some_and(|x| x < kmer_c) {
                dest_pointer.next();
            }

            if dest_pointer.peek().is_some_and(|x| x == kmer_c) {
                rawrow.set_bit(src_pointer.cur_merged_index()-1, true); // -1 because we have advanced past the current k-mer
                dest_pointer.next().unwrap();
            }
        };
        rawrow
    }).collect::<Vec<simple_sds_sbwt::raw_vector::RawVector>>();

    let lcs = if build_lcs {
        // LCS values are between 0 and k-1
        Some(build_lcs_array(KmerDummyMergeSlice::new(&dummies, &kmers, 0..n, k), k))
    } else {
        None
    };

    (rawrows, lcs)
}

use simple_sds_sbwt::ops::Push;
use simple_sds_sbwt::raw_vector::AccessRaw;
use std::cmp::min;
use crate::kmer::LongKmer;
use crate::util::binary_search_leftmost_that_fulfills_pred;

use rayon::iter::IntoParallelIterator;
use rayon::iter::ParallelIterator;
use rayon::prelude::ParallelSlice;

pub fn find_in_dummy<const B: usize>(
    dummy_file: &[(LongKmer<B>, u8)],
    c: u8,
) -> u64 {
    let dummy_file_len = dummy_file.len();

    let access_fn = |pos| {
        dummy_file[pos]
    };

    let pred_fn = |kmer: (LongKmer::<B>,u8)| {
        kmer.1 > 0 && kmer.0.get_from_left(0) >= c
    };

    let start = binary_search_leftmost_that_fulfills_pred(
        access_fn,
        pred_fn,
        dummy_file_len);

    start as u64
}

pub fn find_in_nondummy<const B: usize>(
    nondummy_file: &[LongKmer<B>],
    c: u8,
) -> u64 {
    let nondummy_file_len = nondummy_file.len();

    let access_fn = |pos| {
        nondummy_file[pos]
    };

    let pred_fn = |kmer: LongKmer::<B>| {
        kmer.get_from_left(0) >= c
    };

    let start = binary_search_leftmost_that_fulfills_pred(
        access_fn,
        pred_fn,
        nondummy_file_len);

    start as u64
}

// Returns the LCS array.
// This first builds the uncompressed LCS array, and then compressed it to log(k) bits per element.
// This means that the peak memory is high, but we have the k-mers in memory anyway, so it's only
// twice that for B = 8 (k <= 32). TODO: compute the compressed LCS directly. 
pub fn build_lcs_array<const B: usize>(
    kmers: &[(LongKmer::<B>, u8)],
    k: usize,
) -> simple_sds_sbwt::int_vector::IntVector {
    // LCS values are between 0 and k-1
    assert!(k > 0);

    let bitwidth = 64 - (k as u64 - 1).leading_zeros();

    let non_compressed_lcs = std::iter::once(0_usize).chain(kmers.par_windows(2).map(|window| {
        // The longest common suffix is the longest common prefix of reversed k-mers
        let prev_kmer = &window[0].0;
        let prev_len = &window[0].1;
        let kmer = &window[1].0;
        let len = &window[1].1;
        let mut lcs_value = LongKmer::<B>::lcp(prev_kmer, kmer);
        lcs_value = min(lcs_value, min(*prev_len as usize, *len as usize));
        lcs_value
    }).collect::<Vec<usize>>());

    // Compress into log(k) bits per element
    let mut compressed_lcs = simple_sds_sbwt::int_vector::IntVector::with_capacity(kmers.len(), bitwidth as usize).unwrap();
    for x in non_compressed_lcs {
        compressed_lcs.push(x as u64);
    }
    compressed_lcs
}

// Read the next kmer or dummy from data stored separately in memory
pub fn read_kmer_or_dummy<const B: usize>(
    kmers: &[LongKmer<B>], kmers_pos: &mut usize,
    dummies: &[(LongKmer<B>, u8)], dummies_pos: &mut usize,
    k: usize,
) -> (LongKmer<B>, u8) {
    let n_kmers = kmers.len();
    let n_dummies = dummies.len();

    if *kmers_pos == n_kmers {
        *dummies_pos += 1;
        dummies[*dummies_pos - 1]
    } else if *dummies_pos == n_dummies {
        *kmers_pos += 1;
        (kmers[*kmers_pos - 1], k as u8)
    } else if dummies[*dummies_pos] < (kmers[*kmers_pos], k as u8) {
        *dummies_pos += 1;
        dummies[*dummies_pos - 1]
    } else {
        *kmers_pos += 1;
        (kmers[*kmers_pos - 1], k as u8)
    }
}


pub fn merge_kmers_and_dummies<const B: usize>(
    kmers: &[LongKmer<B>],
    dummies: &[(LongKmer<B>, u8)],
    k: usize,
) -> Vec<(LongKmer<B>, u8)> {

    let n_merged = kmers.len() + dummies.len();

    let mut kmers_pos = 0;
    let mut dummies_pos = 0;

    // TODO: reduce the space overhead
    let mut merged = Vec::<(LongKmer<B>, u8)>::with_capacity(n_merged);
    for _ in 0..n_merged {
        merged.push(read_kmer_or_dummy(kmers, &mut kmers_pos, dummies, &mut dummies_pos, k))
    }
    merged
}

// Returns the SBWT bit vectors and optionally the LCS array
pub fn build_sbwt_bit_vectors<const B: usize>(
    merged: &[(LongKmer<B>, u8)],
    k: usize,
    sigma: usize,
    build_lcs: bool,
) -> (Vec<simple_sds_sbwt::raw_vector::RawVector>, Option<simple_sds_sbwt::int_vector::IntVector>) {

    let n = merged.len();

    let rawrows = (0..sigma).collect::<Vec<usize>>().into_par_iter().map(|c|{
        let mut rawrow = simple_sds_sbwt::raw_vector::RawVector::with_len(n, false);
        let mut pointed_idx = find_in_dummy(merged, c as u8) as usize;
        let end = if c < sigma - 1 { find_in_dummy(merged, c as u8 + 1) as usize } else { n };

        merged.iter().enumerate().for_each(|(kmer_idx, (kmer, len))| {
            let kmer_c = if *len as usize == k {
                (
                    kmer
                        .set_from_left(k - 1, 0)
                        .right_shift(1)
                        .set_from_left(0, c as u8),
                    k as u8,
                )
            } else {
                (kmer.right_shift(1).set_from_left(0, c as u8), len + 1) // Dummy
            };

            while pointed_idx < end && merged[pointed_idx] < kmer_c {
                pointed_idx += 1;
            }

            if pointed_idx < end && merged[pointed_idx] == kmer_c {
                rawrow.set_bit(kmer_idx, true);
                pointed_idx += 1;
            }
        });
        rawrow
    }).collect::<Vec<simple_sds_sbwt::raw_vector::RawVector>>();

    let lcs = if build_lcs {
        // LCS values are between 0 and k-1
        Some(build_lcs_array(merged, k))
    } else {
        None
    };

    (rawrows, lcs)
}

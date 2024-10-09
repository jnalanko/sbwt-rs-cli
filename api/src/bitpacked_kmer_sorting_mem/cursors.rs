use simple_sds_sbwt::{raw_vector::AccessRaw};
use std::cmp::min;
use crate::kmer::LongKmer;
use crate::util::binary_search_leftmost_that_fulfills_pred;

use rayon::iter::IntoParallelIterator;
use rayon::iter::ParallelIterator;
use rayon::prelude::ParallelSlice;

pub fn find_in_dummy<const B: usize>(
    dummy_file: &std::io::Cursor<Vec<(LongKmer<B>, u8)>>,
    c: u8,
) -> u64 {
    let dummy_file_len = dummy_file.get_ref().len();

    let access_fn = |pos| {
        let record = dummy_file.get_ref()[pos];
        record
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
    nondummy_file: &std::io::Cursor<Vec<LongKmer<B>>>,
    c: u8,
) -> u64 {
    let nondummy_file_len = nondummy_file.get_ref().len() as usize;

    let access_fn = |pos| {
        nondummy_file.get_ref()[pos as usize]
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

// Returns the LCS array
pub fn build_lcs_array<const B: usize>(
    kmers: &Vec<(LongKmer::<B>, u8)>,
    k: usize,
) -> simple_sds_sbwt::int_vector::IntVector {
    // LCS values are between 0 and k-1
    assert!(k > 0);

    std::iter::once(0_usize).chain(kmers.par_windows(2).map(|window| {
        // The longest common suffix is the longest common prefix of reversed k-mers
        let prev_kmer = &window[0].0;
        let prev_len = &window[0].1;
        let kmer = &window[1].0;
        let len = &window[1].1;
        let mut lcs_value = LongKmer::<B>::lcp(&prev_kmer, &kmer);
        lcs_value = min(lcs_value, min(*prev_len as usize, *len as usize));
        lcs_value
    }).collect::<Vec<usize>>().into_iter()).collect::<simple_sds_sbwt::int_vector::IntVector>()
}

// Read the next kmer or dummy from data stored separately in memory
pub fn read_kmer_or_dummy<const B: usize>(
    kmers: &mut std::io::Cursor<Vec<LongKmer<B>>>,
    dummies: &mut std::io::Cursor<Vec<(LongKmer<B>, u8)>>,
    k: usize,
) -> (LongKmer<B>, u8) {
    let n_kmers = kmers.get_ref().len();
    let n_dummies = dummies.get_ref().len();
    let kmers_pos = kmers.position() as usize;
    let dummies_pos = dummies.position() as usize;

    if kmers_pos == n_kmers {
        dummies.set_position(dummies_pos as u64 + 1);
        dummies.get_ref()[dummies_pos]
    } else if dummies_pos == n_dummies {
        kmers.set_position(kmers_pos as u64 + 1);
        (kmers.get_ref()[kmers_pos], k as u8)
    } else {
        if dummies.get_ref()[dummies_pos] < (kmers.get_ref()[kmers_pos], k as u8) {
            dummies.set_position(dummies_pos as u64 + 1);
            dummies.get_ref()[dummies_pos]
        } else {
            kmers.set_position(kmers_pos as u64 + 1);
            (kmers.get_ref()[kmers_pos], k as u8)
        }
    }
}


pub fn merge_kmers_and_dummies<const B: usize>(
    kmers: &mut std::io::Cursor<Vec<LongKmer<B>>>,
    dummies: &mut std::io::Cursor<Vec<(LongKmer<B>, u8)>>,
    k: usize,
) -> std::io::Cursor<Vec<(LongKmer<B>, u8)>> {

    let n_merged = kmers.get_ref().len() + dummies.get_ref().len();

    std::io::Cursor::new(Vec::from((0..n_merged).map(|_| {
        read_kmer_or_dummy(kmers, dummies, k)
    }).collect::<Vec<(LongKmer::<B>, u8)>>()))
}

// Returns the SBWT bit vectors and optionally the LCS array
pub fn build_sbwt_bit_vectors<const B: usize>(
    merged: &std::io::Cursor<Vec<(LongKmer<B>, u8)>>,
    k: usize,
    sigma: usize,
    build_lcs: bool,
) -> (Vec<simple_sds_sbwt::raw_vector::RawVector>, Option<simple_sds_sbwt::int_vector::IntVector>) {

    let n = merged.get_ref().len();

    let rawrows = (0..sigma).collect::<Vec<usize>>().into_par_iter().map(|c|{
        let mut rawrow = simple_sds_sbwt::raw_vector::RawVector::with_len(n, false);
        let mut pointed_idx = find_in_dummy(&merged, c as u8) as usize;
        let end = if c < sigma - 1 { find_in_dummy(&merged, c as u8 + 1) as usize } else { n };

        merged.get_ref().iter().enumerate().for_each(|(kmer_idx, (kmer, len))| {
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

            while pointed_idx < end && merged.get_ref()[pointed_idx] < kmer_c {
                pointed_idx += 1;
            }

            if pointed_idx < end && merged.get_ref()[pointed_idx] == kmer_c {
                rawrow.set_bit(kmer_idx, true);
                pointed_idx += 1;
            }
        });
        rawrow
    }).collect::<Vec<simple_sds_sbwt::raw_vector::RawVector>>();

    let lcs = if build_lcs {
        // LCS values are between 0 and k-1
        Some(build_lcs_array(merged.get_ref(), k))
    } else {
        None
    };

    (rawrows, lcs)
}

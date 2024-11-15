//! Build an [SbwtIndex] fully in memory using an algorithm based on bitpacked k-mer sorting.

mod dummies;
mod kmer_splitter;
mod cursors;

use crate::{sbwt::{PrefixLookupTable, SbwtIndex}, streaming_index::LcsArray, subsetseq::SubsetSeq, util::DNA_ALPHABET};
/// Build SBWT and optionally the LCS array fully in memory using bitpacked k-mer sorting.
///
/// See [SbwtIndexBuilder](crate::builder::SbwtIndexBuilder) for a wrapper with a more
/// user-friendly interface. B is the number u64 words in a k-mer.
///
/// Unused arguments for signature compatibility with the disk based building algorithm.
pub fn build_with_bitpacked_kmer_sorting<const B: usize, IN: crate::SeqStream + Send, SS: SubsetSeq + Send>(
    seqs: IN,
    k: usize,
    dedup_batches: bool,
    build_lcs: bool,
) -> (SbwtIndex::<SS>, Option<LcsArray>) {

    let sigma = DNA_ALPHABET.len();

    log::info!("Splitting k-mers into bins");
    let mut kmer_bins = kmer_splitter::split_to_bins::<B, IN>(seqs, k, dedup_batches);

    log::info!("Sorting and deduplicating bins");
    kmer_splitter::par_sort_and_dedup_bin_files::<B>(&mut kmer_bins);

    let (merged, n_kmers) = {
        let mut kmers = kmer_splitter::concat_files(&mut kmer_bins);
        let n_kmers = kmers.get_ref().len();

        log::info!("{} distinct k-mers found", n_kmers);
        let mut dummies = dummies::get_sorted_dummies::<B>(&mut kmers, sigma, k);

        (cursors::merge_kmers_and_dummies(&mut kmers, &mut dummies, k), n_kmers)
    };

    log::info!("Constructing the sbwt subset sequence");

    let (rawrows, lcs) = cursors::build_sbwt_bit_vectors::<B>(&merged, k, sigma, build_lcs);

    // Create the C array
    #[allow(non_snake_case)] // C-array is an established convention in BWT indexes
    let C: Vec<usize> = crate::util::get_C_array(&rawrows);

    log::info!("Building the subset rank structure");
    let mut subsetseq = SS::new_from_bit_vectors(rawrows.into_iter().map(simple_sds_sbwt::bit_vector::BitVector::from).collect());
    subsetseq.build_rank();
    let n_sets = subsetseq.len();
    let (mut index, lcs) = (SbwtIndex::<SS>::from_components(
        subsetseq,
        n_kmers,
        k,
        C,
        PrefixLookupTable::new_empty(n_sets))
                            , lcs.map(LcsArray::new));

    let lut = PrefixLookupTable::new(&index, 8);
    index.set_lookup_table(lut);
    (index, lcs)

}

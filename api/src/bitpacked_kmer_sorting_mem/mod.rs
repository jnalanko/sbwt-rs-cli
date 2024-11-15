//! Build an [SbwtIndex] fully in memory using an algorithm based on bitpacked k-mer sorting.

mod dummies;
mod kmer_splitter;
mod cursors;

use kmer_splitter::get_bitpacked_sorted_distinct_kmers;

use crate::{sbwt::{PrefixLookupTable, SbwtIndex}, streaming_index::LcsArray, subsetseq::SubsetSeq, util::DNA_ALPHABET};
/// Build SBWT and optionally the LCS array fully in memory using bitpacked k-mer sorting.
///
/// See [SbwtIndexBuilder](crate::builder::SbwtIndexBuilder) for a wrapper with a more
/// user-friendly interface. B is the number u64 words in a k-mer.
///
/// Unused arguments for signature compatibility with the disk based building algorithm.
pub fn build_with_bitpacked_kmer_sorting<const B: usize, IN: crate::SeqStream + Send, SS: SubsetSeq + Send>(
    mut seqs: IN,
    k: usize,
    n_threads: usize,
    dedup_batches: bool,
    build_lcs: bool,
) -> (SbwtIndex::<SS>, Option<LcsArray>) {

    let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();

    thread_pool.install(||{

        let sigma = DNA_ALPHABET.len();

        log::info!("Bit-packing and sorting all k-mers of all input sequences (dedup batches: {}).", dedup_batches);
        let mut rev_kmers = get_bitpacked_sorted_distinct_kmers::<B, IN>(seqs, k, n_threads, dedup_batches);

        let (merged, n_kmers) = {
            let n_kmers = rev_kmers.len();
            log::info!("{} distinct k-mers found", n_kmers);
            let dummies = dummies::get_sorted_dummies::<B>(&rev_kmers, sigma, k);
            (cursors::merge_kmers_and_dummies(&rev_kmers, &dummies, k), n_kmers)
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

        })

}

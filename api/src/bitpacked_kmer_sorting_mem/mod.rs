//! Build an [SbwtIndex] fully in memory using an algorithm based on bitpacked k-mer sorting.

mod dummies;
mod kmer_splitter;
mod cursors;

use dummies::KmersWithLengths;
use kmer_splitter::get_bitpacked_sorted_distinct_kmers;

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
    n_threads: usize,
    approx_mem_gb: usize,
    dedup_batches: bool,
    build_lcs: bool,
    add_all_dummy_paths: bool,
) -> (SbwtIndex::<SS>, Option<LcsArray>) {

    let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();

    thread_pool.install(||{

        let sigma = DNA_ALPHABET.len();

        log::info!("Bit-packing and sorting all k-mers of all input sequences (dedup batches: {}).", dedup_batches);
        let (rev_kmers, rev_first_mers) = get_bitpacked_sorted_distinct_kmers::<B, IN>(seqs, k, n_threads, dedup_batches, add_all_dummy_paths, approx_mem_gb);

        if let Some(rev_first_mers) = rev_first_mers {
            println!("TODO: DO SOMETHING WITH THESE -MERS");
            for (rev_mer, len) in rev_first_mers {
                let mut s = rev_mer.to_string()[0..len].as_bytes().to_owned();
                s.reverse();
                eprintln!("{}", String::from_utf8_lossy(&s));
            }
        }

        let n_kmers = rev_kmers.len();
        log::info!("{} distinct k-mers found", n_kmers);

        log::info!("Constructing dummy k-mers");
        let dummies = dummies::get_sorted_dummies::<B>(&rev_kmers, sigma, k, n_threads);

        log::info!("{} dummy nodes constructed", dummies.len());

        log::info!("Compacting dummies and lengths");
        let dummies_and_lengths = KmersWithLengths{ // This fixes word alignment overhead
            kmers: dummies.iter().map(|pair| pair.0).collect(),
            lengths: dummies.into_iter().map(|pair| pair.1).collect()
        };

        log::info!("Constructing the sbwt subset sequence");

        let (rawrows, lcs) = cursors::build_sbwt_bit_vectors::<B>(rev_kmers, dummies_and_lengths, k, sigma, build_lcs, n_threads);

        // Create the C array
        #[allow(non_snake_case)] // C-array is an established convention in BWT indexes
        let C: Vec<usize> = crate::util::get_C_array(&rawrows);

        log::info!("Building the subset rank structure");
        let mut subsetseq = SS::new_from_bit_vectors(rawrows);
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

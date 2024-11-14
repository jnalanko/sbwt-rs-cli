//! Build an [SbwtIndex] fully in memory using an algorithm based on bitpacked k-mer sorting.

mod dummies;
mod kmer_splitter;
mod cursors;

use rayon::iter::IntoParallelRefIterator;
use rayon::iter::ParallelIterator;
use rayon::prelude::ParallelSliceMut;
use crate::kmer::LongKmer;

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

        log::info!("Reading input sequences into memory");

        let mut input_sequences = Vec::<Vec<u8>>::new();
        while let Some(seq) = seqs.stream_next() {
            input_sequences.push(seq.to_owned());
        }


        log::info!("Bit-packing all k-mers of all input sequences.");
        let mut kmers = input_sequences.par_iter().map(|seq| {
            seq.windows(k).filter_map(|bytes| {
                LongKmer::<B>::from_ascii(bytes).ok()
            }).collect::<Vec<LongKmer::<B>>>()
        }).flatten().collect::<Vec<LongKmer::<B>>>();


        log::info!("Sorting all k-mers");
        kmers.par_sort_unstable_by_key(|kmer| {
            kmer.get_from_left(0) as usize * 16 + kmer.get_from_left(1) as usize * 4 + kmer.get_from_left(2) as usize
        });

        log::info!("Deduplicating k-mers");
        kmers.dedup();

        let (merged, n_kmers) = {
            let n_kmers = kmers.len();
            log::info!("{} distinct k-mers found", n_kmers);
            let dummies = dummies::get_sorted_dummies::<B>(&kmers, sigma, k);
            (cursors::merge_kmers_and_dummies(&kmers, &dummies, k), n_kmers)
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

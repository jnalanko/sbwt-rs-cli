//! Build an [SbwtIndex] fully in memory using an algorithm based on bitpacked k-mer sorting.

mod dummies;
mod kmer_splitter;
mod cursors;

use std::cmp::max;

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
) -> (SbwtIndex::<SS>, Option<LcsArray>) {

    let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();

    thread_pool.install(||{

        let sigma = DNA_ALPHABET.len();
        let kmer_bytes = std::mem::size_of::<crate::kmer::LongKmer<B>>();

        log::info!("Bit-packing and sorting all k-mers of all input sequences (dedup batches: {}).", dedup_batches);
        // Buffer sizes (todo: calculate from memory budged)
        let producer_buf_cap = 1_usize << 20;
        let thread_local_buf_caps = 1_usize << 16; 
        let shared_buf_caps = if dedup_batches { 
            let mem_remaining = approx_mem_gb as isize * (1 << 30) - producer_buf_cap as isize - thread_local_buf_caps as isize;
            let mem_remaining = max(mem_remaining, 1_isize << 30) as usize; // At least 1 GB

            // The total memory in the buffers will be:
            // 64 * buf_cap * sizeof(LongKmer<B>>
            // Set this equal to mem_remaining and solve for buf_cap:
            let cap = mem_remaining / 64 / std::mem::size_of::<crate::kmer::LongKmer<B>>();
            max(cap, (1_usize << 30) / 64) // At least 1 GB total
        } else {
            1 << 20
        };

        // One bin per thread per 3-mer
        let thread_local_total = thread_local_buf_caps * n_threads * 64 * kmer_bytes; 
        let shared_total = shared_buf_caps * 64 * kmer_bytes; // One bin per 3-mer
        log::info!("Producer buffer capacity: {}", human_bytes::human_bytes(producer_buf_cap as f64));
        log::info!("Thread-local buffer capacity: {} ({} total)", human_bytes::human_bytes(thread_local_buf_caps as f64), human_bytes::human_bytes(thread_local_total as f64)); // 64 local bins per thread
        log::info!("Shared bin buffer capacity: {} ({} total)", human_bytes::human_bytes(shared_buf_caps as f64), human_bytes::human_bytes(shared_total as f64));
        if producer_buf_cap + thread_local_total + shared_total > approx_mem_gb * (1_usize << 30) {
            log::warn!("Exceeding memory budget");
        }

        let rev_kmers = get_bitpacked_sorted_distinct_kmers::<B, IN>(seqs, k, n_threads, dedup_batches, producer_buf_cap, thread_local_buf_caps, shared_buf_caps);

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

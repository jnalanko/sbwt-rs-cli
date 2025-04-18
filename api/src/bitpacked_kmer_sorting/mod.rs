//! Build an [SbwtIndex] using an algorithm based on bitpacked k-mer sorting.

mod dummies;
mod kmer_splitter;
mod cursors;
mod kmer_chunk;

use crate::kmer::LongKmer;
use human_bytes::human_bytes;

use std::{cmp::max, io::Seek, path::Path};

use crate::{sbwt::{PrefixLookupTable, SbwtIndex}, streaming_index::LcsArray, subsetseq::SubsetSeq, tempfile::TempFileManager, util::DNA_ALPHABET};

fn file_size(path: &Path) -> usize {
    std::fs::metadata(path).unwrap().len() as usize
}

/// Build using bitpacked k-mer sorting. See [SbwtIndexBuilder](crate::builder::SbwtIndexBuilder) for a wrapper with a more 
/// user-friendly interface. B is the number u64 words in a k-mer.
pub fn build_with_bitpacked_kmer_sorting<const B: usize, IN: crate::SeqStream + Send, SS: SubsetSeq + Send>(seqs: IN, k: usize, mem_gb: usize, n_threads: usize, dedup_batches: bool, build_lcs: bool, temp_file_manager: &mut TempFileManager) -> (SbwtIndex::<SS>, Option<LcsArray>) {

    let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();

    thread_pool.install(||{
        let sigma = DNA_ALPHABET.len(); 

        log::info!("Splitting k-mers into bins");
        let (bin_files, n_bytes_in_bins) = kmer_splitter::split_to_bins::<B, IN>(seqs, k, mem_gb, n_threads, dedup_batches, temp_file_manager);

        log::info!("Total size of k-mer bins: {} bytes ({})", n_bytes_in_bins, human_bytes(n_bytes_in_bins as f64));

        log::info!("Sorting and deduplicating bins");
        let (bin_files, n_bytes_after_dedup) = kmer_splitter::par_sort_and_dedup_bin_files::<B>(bin_files, mem_gb, n_threads);

        log::info!("Total size of deduplicated k-mer bins: {} bytes ({})", n_bytes_after_dedup, human_bytes(n_bytes_after_dedup as f64));

        let mut kmers_file = temp_file_manager.create_new_file("kmers-", 10, ".bin");
        let concat_space_overhead = kmer_splitter::concat_files(bin_files, &mut kmers_file.file);

        let concat_space_peak = n_bytes_after_dedup + concat_space_overhead;
        log::info!("Disk peak space during concatenation: {} bytes ({})", concat_space_peak, human_bytes(concat_space_peak as f64));
        kmers_file.file.seek(std::io::SeekFrom::Start(0)).unwrap();

        let n_kmers = file_size(&kmers_file.path) / LongKmer::<B>::byte_size();

        log::info!("{} distinct k-mers found", n_kmers);

        let required_dummies = dummies::get_sorted_dummies::<B>(&kmers_file.path, sigma, k, temp_file_manager);

        log::info!("{} dummy nodes needed", required_dummies.len());

        let n = n_kmers + required_dummies.len();

        // Write dummies to disk
        log::info!("Writing dummies to disk");
        let mut dummy_file = temp_file_manager.create_new_file("dummies-", 10, ".bin");
        dummies::write_to_disk(required_dummies, &mut dummy_file.file);

        let dummy_merge_peak = file_size(&kmers_file.path) + file_size(&dummy_file.path);
        let disk_peak_total = max(max(dummy_merge_peak, concat_space_peak), n_bytes_in_bins);
        log::info!("Temporary disk space peak: {} bytes ({})", disk_peak_total, human_bytes(disk_peak_total as f64));

        log::info!("Constructing the sbwt subset sequence");

        let char_cursors = cursors::init_char_cursors::<B>(&dummy_file.path, &kmers_file.path, k, sigma);

        let global_cursor = cursors::DummyNodeMerger::new(
            std::io::BufReader::new(std::fs::File::open(&dummy_file.path).unwrap()),
            std::io::BufReader::new(std::fs::File::open(&kmers_file.path).unwrap()),
            k,
        );

        let (rawrows, lcs) = cursors::build_sbwt_bit_vectors(global_cursor, char_cursors, n, k, sigma, build_lcs);

        drop(kmers_file); // Free disk space
        drop(dummy_file); // Free disk space

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
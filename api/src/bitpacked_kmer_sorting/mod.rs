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
pub fn build_with_bitpacked_kmer_sorting<const B: usize, IN: crate::SeqStream + Send, SS: SubsetSeq + Send>(seqs: IN, k: usize, mem_gb: usize, n_threads: usize, dedup_batches: bool, build_lcs: bool, add_all_dummy_paths: bool, temp_file_manager: &mut TempFileManager) -> (SbwtIndex::<SS>, Option<LcsArray>) {

    let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();

    thread_pool.install(||{
        let sigma = DNA_ALPHABET.len(); 

        log::info!("Splitting k-mers into bins");
        let (bin_files, n_bytes_in_bins, first_mers) = kmer_splitter::split_to_bins::<B, IN>(seqs, k, mem_gb, n_threads, dedup_batches, add_all_dummy_paths, temp_file_manager);

        log::info!("Total size of k-mer bins: {} bytes ({})", n_bytes_in_bins, human_bytes(n_bytes_in_bins as f64));

        log::info!("Sorting and deduplicating bins");
        let (bin_files, n_bytes_after_dedup) = kmer_splitter::par_sort_and_dedup_bin_files::<B>(bin_files, mem_gb, n_threads);

        log::info!("Total size of deduplicated k-mer bins: {} bytes ({})", n_bytes_after_dedup, human_bytes(n_bytes_after_dedup as f64));

        let mut kmers_file = temp_file_manager.create_new_file("kmers-", 10, ".bin");

        log::info!("Concatenating deduplicated bins");
        let concat_space_overhead = kmer_splitter::concat_files(bin_files, &mut kmers_file.file);

        let concat_space_peak = n_bytes_after_dedup + concat_space_overhead;
        log::info!("Disk peak space during concatenation: {} bytes ({})", concat_space_peak, human_bytes(concat_space_peak as f64));
        kmers_file.file.seek(std::io::SeekFrom::Start(0)).unwrap();

        let n_kmers = file_size(&kmers_file.path) / LongKmer::<B>::byte_size();

        log::info!("{} distinct k-mers found", n_kmers);

        let required_dummies = dummies::get_sorted_dummies::<B>(&kmers_file.path, sigma, k, temp_file_manager, first_mers);

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

#[cfg(test)]
mod tests {
    use std::cmp::min;
    use crate::util::is_dna;

    #[test]
    fn test_add_all_dummy_paths() {
        let k = 3;
        let seq1 = b"ACACTGANNNNNGATANNNNTGNNNAN";
        let seq2 = b"NNNNNNNNNN";
        let seq3 = b"ACCCCCCCNNNNTACNNN";
        let seqs: Vec<&[u8]> = vec![seq1, seq2, seq3];
        let mut temp_file_manager = crate::tempfile::TempFileManager::new(std::path::Path::new("/tmp"));
        let (mut sbwt, _lcs) = super::build_with_bitpacked_kmer_sorting::<1, _, crate::subsetseq::SubsetMatrix>(
            crate::SliceSeqStream::new(seqs.as_slice()),
            k,
            1,     // mem_gb
            1,     // n_threads
            false, // dedup_batches
            true,  // build_lcs
            true,  // add_all_dummy_paths
            &mut temp_file_manager,
        );
        sbwt.build_select();

        for seq in [seq1.as_ref(), seq2.as_ref(), seq3.as_ref()] {
            crate::util::for_each_run_with_key(seq, |c| crate::util::is_dna(*c), |run_range| {
                if crate::util::is_dna(seq[run_range.start]) {
                    let s = &seq[run_range];
                    // All dummy prefixes up to length k-1 should be in the index
                    for l in 1..min(k, s.len()) {
                        let colex_range = sbwt.search(&s[0..l]).unwrap();
                        assert!(colex_range.len() > 0);
                        let smallest = sbwt.access_kmer(colex_range.start);
                        let n_dna = smallest.iter().filter(|c| is_dna(**c)).count();
                        assert_eq!(n_dna, l);
                    }
                }
            });
        }

        // Collect true k-mers from the index (rows where all k characters are DNA)
        let index_kmers: std::collections::HashSet<Vec<u8>> = (0..sbwt.n_sets())
            .map(|i| sbwt.access_kmer(i))
            .filter(|kmer| kmer.iter().all(|&c| is_dna(c)))
            .collect();

        // Collect k-mers from the input sequences
        let mut input_kmers: std::collections::HashSet<Vec<u8>> = std::collections::HashSet::new();
        for seq in [seq1.as_ref(), seq2.as_ref(), seq3.as_ref()] {
            crate::util::for_each_run_with_key(seq, |c| is_dna(*c), |run_range| {
                if is_dna(seq[run_range.start]) {
                    let s = &seq[run_range];
                    for i in 0..s.len().saturating_sub(k - 1) {
                        input_kmers.insert(s[i..i + k].to_vec());
                    }
                }
            });
        }

        assert_eq!(index_kmers, input_kmers);
    }
}
//! Build an [SbwtIndex] using an algorithm based on bitpacked k-mer sorting.

mod dummies;
mod kmer_splitter;
mod cursors;
mod kmer_chunk;
mod disk_access;

use human_bytes::human_bytes;

use std::{cmp::max, io::{BufWriter, Seek}, path::Path};

use crate::{sbwt::{PrefixLookupTable, SbwtIndex}, streaming_index::LcsArray, subsetseq::SubsetSeq, tempfile::TempFileManager, util::DNA_ALPHABET};

fn file_size(path: &Path) -> usize {
    std::fs::metadata(path).unwrap().len() as usize
}

// Returns:
// * A file with all the reverse k-mers in lexicogprahic order
// * if `add_all_dummy_paths` is  given, also a file with the reversed up to k characters of each input string.
pub fn sort_and_dedup_kmers_into_file<const B: usize, IN: crate::SeqStream + Send, SS: SubsetSeq + Send>(seqs: IN, k: usize, mem_gb: usize, n_threads: usize, dedup_batches: bool, add_all_dummy_paths: bool, temp_file_manager: &mut TempFileManager) -> (crate::tempfile::TempFile, Option<crate::tempfile::TempFile>){
    let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();

    thread_pool.install(||{

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

        let n_kmers = disk_access::n_kmer_records::<B>(&kmers_file.path);

        log::info!("{} distinct k-mers found", n_kmers);

        if let Some(first_mers) = first_mers {
            let mut first_mers_file = temp_file_manager.create_new_file("sbwt-temp-first-mers", 8, ".bin");
            log::info!("Writing first -mers to {}", first_mers_file.path.display());
            dummies::write_to_disk(first_mers, &mut first_mers_file.file);
            (kmers_file, Some(first_mers_file))
        } else {
            (kmers_file, None)
        }

    })
}

pub fn build_from_kmers_on_disk<const B: usize, SS: SubsetSeq + Send>(k: usize, n_threads: usize, build_lcs: bool, temp_file_manager: &mut TempFileManager, kmers_file: &Path, first_mers_file: Option<&Path>) -> (SbwtIndex::<SS>, Option<LcsArray>) {

    let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();

    thread_pool.install(||{

        let sigma = DNA_ALPHABET.len(); 
        let n_kmers = disk_access::n_kmer_records::<B>(&kmers_file);
        let required_dummies = dummies::get_sorted_dummies::<B>(kmers_file, n_kmers, sigma, k, n_threads, first_mers_file);

        log::info!("{} dummy nodes needed", required_dummies.len());

        let n_dummies = required_dummies.len();

        // Write dummies to disk
        log::info!("Writing dummies to disk");
        let mut dummy_file = temp_file_manager.create_new_file("dummies-", 10, ".bin");
        dummies::write_to_disk(required_dummies, &mut dummy_file.file);

        let mut dummy_merge_peak = file_size(&kmers_file) + file_size(&dummy_file.path);
        if let Some(first_mer_file) = first_mers_file { dummy_merge_peak += file_size(first_mer_file) }
        log::info!("Temporary disk space peak during dummy construction: {} bytes ({})", dummy_merge_peak, human_bytes(dummy_merge_peak as f64));

        log::info!("Constructing the sbwt subset sequence");
        let (rawrows, lcs) = cursors::build_sbwt_bit_vectors::<B>(&dummy_file.path, &kmers_file, n_dummies, n_kmers, k, sigma, build_lcs, n_threads);

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

/// Build using bitpacked k-mer sorting. See [SbwtIndexBuilder](crate::builder::SbwtIndexBuilder) for a wrapper with a more 
/// user-friendly interface. B is the number u64 words in a k-mer.
pub fn build_with_bitpacked_kmer_sorting<const B: usize, IN: crate::SeqStream + Send, SS: SubsetSeq + Send>(seqs: IN, k: usize, mem_gb: usize, n_threads: usize, dedup_batches: bool, build_lcs: bool, add_all_dummy_paths: bool, temp_file_manager: &mut TempFileManager) -> (SbwtIndex::<SS>, Option<LcsArray>) {
    let (kmers_file, first_mers_file) = sort_and_dedup_kmers_into_file::<B, IN, SS>(seqs, k, mem_gb, n_threads, dedup_batches, add_all_dummy_paths, temp_file_manager);

    let first_mers_file_path_option: Option<&Path> = first_mers_file.as_ref().map(|f| f.path.as_path());
    build_from_kmers_on_disk::<B, SS>(k, n_threads, build_lcs, temp_file_manager, &kmers_file.path, first_mers_file_path_option)
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

    #[test]
    fn disk_vs_mem_parallel_differential() {
        use crate::builder::{BitPackedKmerSortingDisk, BitPackedKmerSortingMem, SbwtIndexBuilder};
        use crate::util::gen_random_dna_string;

        // k must be >= BIN_PREFIX_LEN (3), which the mem-based bin splitter requires.
        for &k in &[3_usize, 8, 31] {
            let seq = gen_random_dna_string(500, k as u64);
            for &n_threads in &[1_usize, 2, 3, 8] {
                let (mem_sbwt, mem_lcs) = SbwtIndexBuilder::<BitPackedKmerSortingMem>::new()
                    .k(k).n_threads(n_threads).build_lcs(true)
                    .run_from_slices(&[seq.as_slice()]);

                let (disk_sbwt, disk_lcs) = SbwtIndexBuilder::<BitPackedKmerSortingDisk>::new()
                    .k(k).n_threads(n_threads).build_lcs(true)
                    .run_from_slices(&[seq.as_slice()]);

                assert_eq!(mem_sbwt, disk_sbwt, "SBWT mismatch at k={} n_threads={}", k, n_threads);
                assert_eq!(mem_lcs, disk_lcs, "LCS mismatch at k={} n_threads={}", k, n_threads);
            }
        }
    }

    // Stress test for the suffix-group-leader rewind logic: constructs many k-mers that share
    // the same (k-1)-length suffix (varying only the first character), so that with many
    // threads, segment boundaries are likely to land inside such a "suffix group".
    #[test]
    fn disk_vs_mem_suffix_group_boundary_stress() {
        use crate::builder::{BitPackedKmerSortingDisk, BitPackedKmerSortingMem, SbwtIndexBuilder};
        use crate::util::gen_random_dna_string;

        let k = 4;
        let mut seqs: Vec<Vec<u8>> = vec![];
        for seed in 0..40_u64 {
            let suffix = gen_random_dna_string(k - 1, seed);
            for &first in &[b'A', b'C', b'G', b'T'] {
                let mut kmer = vec![first];
                kmer.extend_from_slice(&suffix);
                seqs.push(kmer);
            }
        }
        let seq_refs: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();

        for &n_threads in &[1_usize, 2, 4, 8, 16] {
            let (mem_sbwt, mem_lcs) = SbwtIndexBuilder::<BitPackedKmerSortingMem>::new()
                .k(k).n_threads(n_threads).build_lcs(true)
                .run_from_slices(&seq_refs);

            let (disk_sbwt, disk_lcs) = SbwtIndexBuilder::<BitPackedKmerSortingDisk>::new()
                .k(k).n_threads(n_threads).build_lcs(true)
                .run_from_slices(&seq_refs);

            assert_eq!(mem_sbwt, disk_sbwt, "SBWT mismatch at n_threads={}", n_threads);
            assert_eq!(mem_lcs, disk_lcs, "LCS mismatch at n_threads={}", n_threads);
        }
    }
}
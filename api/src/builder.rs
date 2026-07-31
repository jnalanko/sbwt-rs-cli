//! A builder pattern interface for building an [SbwtIndex].

use std::str::FromStr;

use crate::JSeqIOSeqStreamWrapper;
use crate::{subsetseq::SubsetMatrix, SeqStream};
use crate::sbwt::{PrefixLookupTable, SbwtIndex};
use crate::streaming_index::LcsArray;

#[derive(Clone, Eq, PartialEq, Debug)]
struct SeqStreamWithPossiblyRevComp<SS: SeqStream + Send>{
    inner: SS, 
    rc_buf: Vec<u8>,
    parity: bool, // Every other sequence we return is a reverse complement of the previous. Initialize to false.
    enable_rev_comp: bool,
}

impl<SS: SeqStream + Send> crate::SeqStream for SeqStreamWithPossiblyRevComp<SS>{
    fn stream_next(&mut self) -> Option<&[u8]> {
        if !self.enable_rev_comp {
            return self.inner.stream_next();
        }

        self.parity = !self.parity;
        if self.parity {
            #[allow(clippy::question_mark)] // More space to write comments on the cases
            let new = match self.inner.stream_next() {
                None => return None, // End of stream
                Some(r) => r // Will return this at the end of the function
            };

            // Store the sequence to rc_buf for the next call
            self.rc_buf.clear();
            self.rc_buf.extend(new);

            Some(new)

        } else {
            jseqio::reverse_complement_in_place(&mut self.rc_buf);
            Some(&self.rc_buf)
        }
    }
}

/// Any struct implementing this interface can be used in [SbwtIndexBuilder] to construct the SBWT index and optionally the LCS array.
pub trait SbwtConstructionAlgorithm {
    fn run(self, k: usize, n_threads: usize, build_lcs: bool, add_all_dummy_paths: bool, add_rev_comp: bool) -> (SbwtIndex<SubsetMatrix>, Option<LcsArray>);
}

/// A construction algorithm based on sorting of bit-packed k-mers using temporary disk space.
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct BitPackedKmerSortingDisk<SS: SeqStream + Send> {
    mem_gb: usize,
    dedup_batches: bool,
    temp_dir: std::path::PathBuf,
    input: SS,
}

impl<SS: SeqStream + Send> BitPackedKmerSortingDisk<SS> {

    /// Initializes the algorithm with default settings:
    /// - 4 GB of memory.
    /// - do not deduplicate k-mer batches before sorting.
    /// - use the current directory as the temporary directory.
    pub fn new(input: SS) -> Self {
        Self{input, mem_gb: 4, dedup_batches: false, temp_dir: std::path::PathBuf::from_str(".").unwrap()}
    }

    /// Set the amount of memory to use in gigabytes. This is not strictly enforced, but the algorithm will try to stay within this limit.
    pub fn mem_gb(mut self, mem_gb: usize) -> Self {
        self.mem_gb = mem_gb;
        self
    }

    /// Whether to deduplicate k-mer batches before sorting. If the input has many duplicate k-mers, this will reduce the disk space required by the algorithm.
    pub fn dedup_batches(mut self, enable: bool) -> Self {
        self.dedup_batches = enable;
        self
    }

    /// Set the temporary directory where the algorithm can store temporary files.
    pub fn temp_dir(mut self, temp_dir: &std::path::Path) -> Self {
        self.temp_dir = temp_dir.to_path_buf();
        std::fs::create_dir_all(&self.temp_dir).unwrap();
        self
    }
}

/// Initialize the algorithm for the given ASCII nucleotide strings, and return the SBWT index and optionally the LCS array if [SbwtIndexBuilder::build_lcs] was set.
pub fn init_bitpacked_kmer_sorting_disk_from_slices<'a>(input: &'a[&'a [u8]]) -> BitPackedKmerSortingDisk<crate::SliceSeqStream<'a>> {
    let stream = crate::util::SliceSeqStream::new(input);
    BitPackedKmerSortingDisk::new(stream)
} 

/// Initialize the algorithm for the given ASCII nucleotide strings, and return the SBWT index and optionally the LCS array if [SbwtIndexBuilder::build_lcs] was set.
pub fn init_bitpacked_kmer_sorting_disk_from_vecs<'a>(input: &'a[Vec<u8>]) -> BitPackedKmerSortingDisk<crate::VecSeqStream<'a>> {
    let stream = crate::util::VecSeqStream::new(input);
    BitPackedKmerSortingDisk::new(stream)
} 

/// Initialize the algorithm from a FASTA-formatted stream of sequences, and return the SBWT index and optionally the LCS array if [SbwtIndexBuilder::build_lcs] was set.
pub fn init_bitpacked_kmer_sorting_disk_from_fasta<R: std::io::Read + Send + Sync + 'static>(input: R) -> BitPackedKmerSortingDisk<JSeqIOSeqStreamWrapper> {
    let input = std::io::BufReader::new(input);
    let input = crate::JSeqIOSeqStreamWrapper{inner: jseqio::reader::DynamicFastXReader::new(input).unwrap()};
    BitPackedKmerSortingDisk::new(input)
} 

/// Initialize the algorithm from a FASTQ-formatted stream of sequences, and return the SBWT index and optionally the LCS array if [SbwtIndexBuilder::build_lcs] was set.
pub fn init_bitpacked_kmer_sorting_disk_from_fastq<R: std::io::Read + Send + Sync + 'static>(input: R) -> BitPackedKmerSortingDisk<JSeqIOSeqStreamWrapper> {
    let input = std::io::BufReader::new(input);
    let input = crate::JSeqIOSeqStreamWrapper{inner: jseqio::reader::DynamicFastXReader::new(input).unwrap()};
    BitPackedKmerSortingDisk::new(input)
} 

impl<SS: SeqStream + Send> SbwtConstructionAlgorithm for BitPackedKmerSortingDisk<SS> {
    fn run(self, k: usize, n_threads: usize, build_lcs: bool, add_all_dummy_paths: bool, add_rev_comp: bool) -> (SbwtIndex<SubsetMatrix>, Option<LcsArray>) {
        let mem_gb = self.mem_gb;
        let dedup_batches = self.dedup_batches;
        let mut temp_file_manager = crate::tempfile::TempFileManager::new(&self.temp_dir);
        let input = SeqStreamWithPossiblyRevComp{ inner: self.input, rc_buf: Vec::<u8>::new(), parity: false, enable_rev_comp: add_rev_comp };
        match k {
            0..=32 => {
                crate::bitpacked_kmer_sorting::build_with_bitpacked_kmer_sorting::<1,_,SubsetMatrix>(input, k, mem_gb, n_threads, dedup_batches, build_lcs, add_all_dummy_paths, &mut temp_file_manager)
            }
            33..=64 => {
                crate::bitpacked_kmer_sorting::build_with_bitpacked_kmer_sorting::<2,_,SubsetMatrix>(input, k, mem_gb, n_threads, dedup_batches, build_lcs, add_all_dummy_paths, &mut temp_file_manager)
            }
            65..=96 => {
                crate::bitpacked_kmer_sorting::build_with_bitpacked_kmer_sorting::<3,_,SubsetMatrix>(input, k, mem_gb, n_threads, dedup_batches, build_lcs, add_all_dummy_paths, &mut temp_file_manager)
            }
            97..=128 => {
                crate::bitpacked_kmer_sorting::build_with_bitpacked_kmer_sorting::<4,_,SubsetMatrix>(input, k, mem_gb, n_threads, dedup_batches, build_lcs, add_all_dummy_paths, &mut temp_file_manager)
            }
            129..=256 => {
                crate::bitpacked_kmer_sorting::build_with_bitpacked_kmer_sorting::<8,_,SubsetMatrix>(input, k, mem_gb, n_threads, dedup_batches, build_lcs, add_all_dummy_paths, &mut temp_file_manager)
            }
            _ => {
                panic!("k > 256 not supported with bitpacked sorting algorithm.");
            }
        }
    }
}

/// A construction algorithm based on sorting of bit-packed k-mers in entirely in RAM.
/// Faster and scales better with parallelism than [BitPackedKmerSortingDisk], but takes
/// more RAM.
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct BitPackedKmerSortingMem<SS: SeqStream + Send>{
    dedup_batches: bool,
    mem_gb: usize,
    input: SS,
}

impl<SS: SeqStream + Send> BitPackedKmerSortingMem<SS> {

    /// Initializes the algorithm with default settings:
    /// - do not deduplicate k-mer batches before sorting.
    /// - 8GB memory
    pub fn new(input: SS) -> Self {
        Self{input, dedup_batches: false, mem_gb: 8}
    }

    /// Whether to deduplicate k-mer batches before sorting. If the input has many duplicate k-mers, this will reduce the disk space required by the algorithm.
    pub fn dedup_batches(mut self, enable: bool) -> Self {
        self.dedup_batches = enable;
        self
    }

    /// Set the amount of memory to use in gigabytes. The deduplication buffer sizes will be tuned to try to
    /// stay within this memory, but there are no guarantees. 
    pub fn mem_gb(mut self, mem_gb: usize) -> Self {
        self.mem_gb = mem_gb;
        self
    }
}

/// Initialize the algorithm for the given ASCII nucleotide strings, and return the SBWT index and optionally the LCS array if [SbwtIndexBuilder::build_lcs] was set.
pub fn init_bitpacked_kmer_sorting_from_slices<'a>(input: &'a[&'a [u8]]) -> BitPackedKmerSortingMem<crate::SliceSeqStream<'a>> {
    let stream = crate::util::SliceSeqStream::new(input);
    BitPackedKmerSortingMem::new(stream)
} 

/// Initialize the algorithm for the given ASCII nucleotide strings, and return the SBWT index and optionally the LCS array if [SbwtIndexBuilder::build_lcs] was set.
pub fn init_bitpacked_kmer_sorting_from_vecs<'a>(input: &'a[Vec<u8>]) -> BitPackedKmerSortingMem<crate::VecSeqStream<'a>> {
    let stream = crate::util::VecSeqStream::new(input);
    BitPackedKmerSortingMem::new(stream)
} 

/// Initialize the algorithm from a FASTA-formatted stream of sequences, and return the SBWT index and optionally the LCS array if [SbwtIndexBuilder::build_lcs] was set.
pub fn init_bitpacked_kmer_sorting_from_fasta<R: std::io::Read + Send + Sync + 'static>(input: R) -> BitPackedKmerSortingMem<JSeqIOSeqStreamWrapper> {
    let input = std::io::BufReader::new(input);
    let input = crate::JSeqIOSeqStreamWrapper{inner: jseqio::reader::DynamicFastXReader::new(input).unwrap()};
    BitPackedKmerSortingMem::new(input)
} 

/// Initialize the algorithm from a FASTQ-formatted stream of sequences, and return the SBWT index and optionally the LCS array if [SbwtIndexBuilder::build_lcs] was set.
pub fn init_bitpacked_kmer_sorting_from_fastq<R: std::io::Read + Send + Sync + 'static>(input: R) -> BitPackedKmerSortingMem<JSeqIOSeqStreamWrapper> {
    let input = std::io::BufReader::new(input);
    let input = crate::JSeqIOSeqStreamWrapper{inner: jseqio::reader::DynamicFastXReader::new(input).unwrap()};
    BitPackedKmerSortingMem::new(input)
} 

impl<SS: SeqStream + Send> SbwtConstructionAlgorithm for BitPackedKmerSortingMem<SS> {
    fn run(self, k: usize, n_threads: usize, build_lcs: bool, add_all_dummy_paths: bool, add_rev_comp: bool) -> (SbwtIndex<SubsetMatrix>, Option<LcsArray>) {
        let dedup_batches = self.dedup_batches;
        let mem_gb = self.mem_gb;
        let input = SeqStreamWithPossiblyRevComp{ inner: self.input, rc_buf: Vec::<u8>::new(), parity: false, enable_rev_comp: add_rev_comp };
        match k {
            0..=32 => {
                crate::bitpacked_kmer_sorting_mem::build_with_bitpacked_kmer_sorting::<1,_,SubsetMatrix>(input, k, n_threads, mem_gb, dedup_batches, build_lcs, add_all_dummy_paths)
            }
            33..=64 => {
                crate::bitpacked_kmer_sorting_mem::build_with_bitpacked_kmer_sorting::<2,_,SubsetMatrix>(input, k, n_threads, mem_gb, dedup_batches, build_lcs, add_all_dummy_paths)
            }
            65..=96 => {
                crate::bitpacked_kmer_sorting_mem::build_with_bitpacked_kmer_sorting::<3,_,SubsetMatrix>(input, k, n_threads, mem_gb, dedup_batches, build_lcs, add_all_dummy_paths)
            }
            97..=128 => {
                crate::bitpacked_kmer_sorting_mem::build_with_bitpacked_kmer_sorting::<4,_,SubsetMatrix>(input, k, n_threads, mem_gb, dedup_batches, build_lcs, add_all_dummy_paths)
            }
            129..=256 => {
                crate::bitpacked_kmer_sorting_mem::build_with_bitpacked_kmer_sorting::<8,_,SubsetMatrix>(input, k, n_threads, mem_gb, dedup_batches, build_lcs, add_all_dummy_paths)
            }
            _ => {
                panic!("k > 256 not supported with bitpacked sorting algorithm.");
            }
        }
    }
}

/// A builder for constructing an SBWT index.
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct SbwtIndexBuilder<A: SbwtConstructionAlgorithm> {
    k: usize,
    n_threads: usize,
    algorithm: A,
    build_lcs: bool,
    add_rev_comp: bool,
    build_select_support: bool,
    precalc_length: usize,
    add_all_dummy_paths: bool,
}

impl<A: SbwtConstructionAlgorithm> SbwtIndexBuilder<A> {

    /// Sets up the builder with default values:
    /// - k = 31.
    /// - n_threads = 4.
    /// - do not build the LCS array.
    /// - do not add the reverse complement of the input sequences.
    /// - do not build the select support.
    /// - precalc_length = 8.
    /// - default settings for the chosen algorithm.
    pub fn new(algorithm: A) -> Self {
        Self{k: 31, n_threads: 4, algorithm, build_lcs: false, add_rev_comp: false, build_select_support: false, precalc_length: 8, add_all_dummy_paths: false}
    }

    /// Sets the k-mer length.
    pub fn k(mut self, k: usize) -> Self {
        self.k = k;
        self
    }

    /// Set the algorithm to use for constructing the SBWT index.
    pub fn algorithm(mut self, alg: A) -> Self {
        self.algorithm = alg;
        self
    }

    /// Whether to build the LCS array.
    pub fn build_lcs(mut self, enable: bool) -> Self {
        self.build_lcs = enable;
        self
    }

    /// Whether to build the select support.
    pub fn build_select_support(mut self, enable: bool) -> Self {
        self.build_select_support = enable;
        self
    }

    /// Whether to add the reverse complement of the input sequences.
    pub fn add_rev_comp(mut self, enable: bool) -> Self {
        self.add_rev_comp = enable;
        self
    }

    /// Set the length of the prefix in the prefix lookup table.
    pub fn precalc_length(mut self, precalc_length: usize) -> Self {
        self.precalc_length = precalc_length;
        self
    }

    /// Set the number of threads to use.
    pub fn n_threads(mut self, n_threads: usize) -> Self {
        self.n_threads = n_threads;
        self
    }

    /// Include all dummy paths for every DNA run in the input, not only those strictly required by the SBWT structure.
    pub fn add_all_dummy_paths(mut self, enable: bool) -> Self {
        self.add_all_dummy_paths = enable;
        self
    }

    /// Run the algorithm and return the SBWT index and optionally the LCS array if [build_lcs](SbwtIndexBuilder::build_lcs) was set.
    /// See also [run_from_slices](SbwtIndexBuilder::run_from_slices), [run_from_fasta](SbwtIndexBuilder::run_from_fasta) and [run_from_fastq](SbwtIndexBuilder::run_from_fastq).
    pub fn run(self) -> (SbwtIndex<SubsetMatrix>, Option<LcsArray>) {
        let (mut sbwt, lcs) = self.algorithm.run(self.k, self.n_threads, self.build_lcs, self.add_all_dummy_paths, self.add_rev_comp);

        if self.build_select_support {
            sbwt.build_select();
        }

        if sbwt.get_lookup_table().prefix_length != self.precalc_length {
            let lut = PrefixLookupTable::new(&sbwt, self.precalc_length);
            sbwt.set_lookup_table(lut);
        }

        (sbwt, lcs)
    } 
}

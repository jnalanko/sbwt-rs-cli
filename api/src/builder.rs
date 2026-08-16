//! A builder pattern interface for building an [SbwtIndex].

use std::path::Path;
use std::str::FromStr;

use crate::tempfile::TempFileManager;
use crate::{subsetseq::SubsetMatrix, SeqStream};
use crate::sbwt::{PrefixLookupTable, SbwtIndex};
use crate::streaming_index::LcsArray;
use crate::util::SeqStreamWithPossiblyRevComp;
/// A construction algorithm based on sorting of bit-packed k-mers using temporary disk space.
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct BitPackedKmerSortingDisk<SS: SeqStream + Send> {
    mem_gb: usize,
    dedup_batches: bool,
    temp_dir: std::path::PathBuf,
    input: SS,

    k: usize,
    n_threads: usize,
    build_lcs: bool,
    add_rev_comp: bool,
    build_select_support: bool,
    precalc_length: usize,
    add_all_dummy_paths: bool,
}

impl<SS: SeqStream + Send> BitPackedKmerSortingDisk<SS> {

    /// Initializes the algorithm for a [crate::SeqStream] with a given k and default settings.
    /// - Use the current directory as the temporary directory.
    /// - 8 GiB memory
    /// - 4 threads
    /// - no LCS array
    /// - no reverse complements added
    /// - no select support
    /// - precalc length min(8,k)
    /// - add_all_dummy_paths = false
    /// - dedup_batches = false
    pub fn new(input: SS, k: usize) -> Self {
        Self{input, dedup_batches: false, mem_gb: 8, k, n_threads: 4, build_lcs: false, add_rev_comp: false, build_select_support: false, precalc_length: std::cmp::min(8, k), add_all_dummy_paths: false, temp_dir: std::path::PathBuf::from_str(".").unwrap()}
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

    /// Sets the k-mer length.
    pub fn k(mut self, k: usize) -> Self {
        self.k = k;
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

    /// Run the construction algorithm, consuming the input stream, and return the SBWT index and
    /// the LCS array if [build_lcs](BitPackedKmerSortingDisk::build_lcs) was set, otherwise `None`.
    ///
    /// Panics if k > 256.
    pub fn run(self) -> (SbwtIndex<SubsetMatrix>, Option<LcsArray>) {
        let mut temp_file_manager = crate::tempfile::TempFileManager::new(&self.temp_dir);
        let input = SeqStreamWithPossiblyRevComp::new(self.input, self.add_rev_comp);
        let (mut sbwt, lcs) = match self.k {
            0..=32 => {
                crate::bitpacked_kmer_sorting::build_with_bitpacked_kmer_sorting::<1,_,SubsetMatrix>(input, self.k, self.mem_gb, self.n_threads, self.dedup_batches, self.build_lcs, self.add_all_dummy_paths, &mut temp_file_manager)
            }
            33..=64 => {
                crate::bitpacked_kmer_sorting::build_with_bitpacked_kmer_sorting::<2,_,SubsetMatrix>(input, self.k, self.mem_gb, self.n_threads, self.dedup_batches, self.build_lcs, self.add_all_dummy_paths, &mut temp_file_manager)
            }
            65..=96 => {
                crate::bitpacked_kmer_sorting::build_with_bitpacked_kmer_sorting::<3,_,SubsetMatrix>(input, self.k, self.mem_gb, self.n_threads, self.dedup_batches, self.build_lcs, self.add_all_dummy_paths, &mut temp_file_manager)
            }
            97..=128 => {
                crate::bitpacked_kmer_sorting::build_with_bitpacked_kmer_sorting::<4,_,SubsetMatrix>(input, self.k, self.mem_gb, self.n_threads, self.dedup_batches, self.build_lcs, self.add_all_dummy_paths, &mut temp_file_manager)
            }
            129..=256 => {
                crate::bitpacked_kmer_sorting::build_with_bitpacked_kmer_sorting::<8,_,SubsetMatrix>(input, self.k, self.mem_gb, self.n_threads, self.dedup_batches, self.build_lcs, self.add_all_dummy_paths, &mut temp_file_manager)
            }
            _ => {
                panic!("k > 256 not supported with a bitpacked sorting algorithm.");
            }
        };

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

impl<'a> BitPackedKmerSortingDisk<crate::SliceSeqStream<'a>> {
    /// Initialize the algorithm for a slice of ASCII sequences, with the same defaults as [BitPackedKmerSortingDisk::new].
    pub fn new_from_slices(input: &'a[&'a [u8]], k: usize) -> Self {
        let stream = crate::util::SliceSeqStream::new(input);
        BitPackedKmerSortingDisk::new(stream, k)
    }
}

impl<'a> BitPackedKmerSortingDisk<crate::VecSeqStream<'a>> {
    /// Initialize the algorithm for a slice of owned ASCII sequences, with the same defaults as [BitPackedKmerSortingDisk::new].
    pub fn new_from_vecs(input: &'a[Vec<u8>], k: usize) -> Self {
        let stream = crate::util::VecSeqStream::new(input);
        BitPackedKmerSortingDisk::new(stream, k)
    }
}

impl BitPackedKmerSortingDisk<crate::util::FastXReader> {
    /// Initialize the algorithm for a fasta file, with the same defaults as [BitPackedKmerSortingDisk::new].
    pub fn new_from_fasta<R: std::io::Read + Send + Sync + 'static>(input: R, k: usize) -> Self {
        let input = std::io::BufReader::new(input);
        let input = crate::util::FastXReader{inner: jseqio::reader::DynamicFastXReader::new(input).unwrap()};
        BitPackedKmerSortingDisk::new(input, k)
    }

    /// Initialize the algorithm for a fastq file, with the same defaults as [BitPackedKmerSortingDisk::new].
    pub fn new_from_fastq<R: std::io::Read + Send + Sync + 'static>(input: R, k: usize) -> Self {
        let input = std::io::BufReader::new(input);
        let input = crate::util::FastXReader{inner: jseqio::reader::DynamicFastXReader::new(input).unwrap()};
        BitPackedKmerSortingDisk::new(input, k)
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

    k: usize,
    n_threads: usize,
    build_lcs: bool,
    add_rev_comp: bool,
    build_select_support: bool,
    precalc_length: usize,
    add_all_dummy_paths: bool,
}

impl<SS: SeqStream + Send> BitPackedKmerSortingMem<SS> {

    /// Initializes the algorithm for a [crate::SeqStream] with a given k and default settings.
    /// - 8 GiB memory
    /// - 4 threads
    /// - no LCS array
    /// - no reverse complements added
    /// - no select support
    /// - precalc length min(8,k)
    /// - add_all_dummy_paths = false
    /// - dedup_batches = false
    pub fn new(input: SS, k: usize) -> Self {
        Self{input, dedup_batches: false, mem_gb: 8, k, n_threads: 4, build_lcs: false, add_rev_comp: false, build_select_support: false, precalc_length: std::cmp::min(8, k), add_all_dummy_paths: false}
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

    /// Sets the k-mer length.
    pub fn k(mut self, k: usize) -> Self {
        self.k = k;
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

    /// Run the construction algorithm, consuming the input stream, and return the SBWT index and
    /// the LCS array if [build_lcs](BitPackedKmerSortingMem::build_lcs) was set, otherwise `None`.
    ///
    /// Panics if k > 256.
    pub fn run(self) -> (SbwtIndex<SubsetMatrix>, Option<LcsArray>) {
        let input = SeqStreamWithPossiblyRevComp::new(self.input, self.add_rev_comp);
        let (mut sbwt, lcs) = match self.k {
            0..=32 => {
                crate::bitpacked_kmer_sorting_mem::build_with_bitpacked_kmer_sorting::<1,_,SubsetMatrix>(input, self.k, self.n_threads, self.mem_gb, self.dedup_batches, self.build_lcs, self.add_all_dummy_paths)
            }
            33..=64 => {
                crate::bitpacked_kmer_sorting_mem::build_with_bitpacked_kmer_sorting::<2,_,SubsetMatrix>(input, self.k, self.n_threads, self.mem_gb, self.dedup_batches, self.build_lcs, self.add_all_dummy_paths)
            }
            65..=96 => {
                crate::bitpacked_kmer_sorting_mem::build_with_bitpacked_kmer_sorting::<3,_,SubsetMatrix>(input, self.k, self.n_threads, self.mem_gb, self.dedup_batches, self.build_lcs, self.add_all_dummy_paths)
            }
            97..=128 => {
                crate::bitpacked_kmer_sorting_mem::build_with_bitpacked_kmer_sorting::<4,_,SubsetMatrix>(input, self.k, self.n_threads, self.mem_gb, self.dedup_batches, self.build_lcs, self.add_all_dummy_paths)
            }
            129..=256 => {
                crate::bitpacked_kmer_sorting_mem::build_with_bitpacked_kmer_sorting::<8,_,SubsetMatrix>(input, self.k, self.n_threads, self.mem_gb, self.dedup_batches, self.build_lcs, self.add_all_dummy_paths)
            }
            _ => {
                panic!("k > 256 not supported with a bitpacked sorting algorithm.");
            }
        };

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

impl<'a> BitPackedKmerSortingMem<crate::SliceSeqStream<'a>> {
    /// Initialize the algorithm for a slice of ASCII sequences, with the same defaults as [BitPackedKmerSortingMem::new].
    pub fn new_from_slices(input: &'a[&'a [u8]], k: usize) -> Self {
        let stream = crate::util::SliceSeqStream::new(input);
        BitPackedKmerSortingMem::new(stream, k)
    }
}

impl<'a> BitPackedKmerSortingMem<crate::VecSeqStream<'a>> {
    /// Initialize the algorithm for a slice of owned ASCII sequences, with the same defaults as [BitPackedKmerSortingMem::new].
    pub fn new_from_vecs(input: &'a[Vec<u8>], k: usize) -> Self {
        let stream = crate::util::VecSeqStream::new(input);
        BitPackedKmerSortingMem::new(stream, k)
    }
}

impl BitPackedKmerSortingMem<crate::util::FastXReader> {
    /// Initialize the algorithm for a fasta file, with the same defaults as [BitPackedKmerSortingMem::new].
    pub fn new_from_fasta<R: std::io::Read + Send + Sync + 'static>(input: R, k: usize) -> Self {
        let input = std::io::BufReader::new(input);
        let input = crate::util::FastXReader{inner: jseqio::reader::DynamicFastXReader::new(input).unwrap()};
        BitPackedKmerSortingMem::new(input, k)
    }

    /// Initialize the algorithm for a fastq file, with the same defaults as [BitPackedKmerSortingMem::new].
    pub fn new_from_fastq<R: std::io::Read + Send + Sync + 'static>(input: R, k: usize) -> Self {
        let input = std::io::BufReader::new(input);
        let input = crate::util::FastXReader{inner: jseqio::reader::DynamicFastXReader::new(input).unwrap()};
        BitPackedKmerSortingMem::new(input, k)
    }
}

/// Wraps a [SeqStream] so that each sequence is sanitized (characters outside of the DNA
/// alphabet are replaced with the separator `$`) and reversed, which is the form the
/// bounded suffix sorting construction expects the concatenation to be in.
struct SanitizedReversedSeqStream<SS: SeqStream + Send> {
    inner: SS,
    buf: Vec<u8>,
}

impl<SS: SeqStream + Send> crate::SeqStream for SanitizedReversedSeqStream<SS> {
    fn stream_next(&mut self) -> Option<&[u8]> {
        // `self.inner` and `self.buf` are disjoint fields, so the borrow of the former
        // does not prevent us from writing to the latter.
        let seq = self.inner.stream_next()?;
        self.buf.clear();
        self.buf.extend(seq);
        crate::alternative_construction::preprocessing::sanitise(&mut self.buf);
        self.buf.reverse();
        Some(&self.buf)
    }
}

/// A construction algorithm based on sorting the k-bounded contexts of the input.
/// Unlike [BitPackedKmerSortingMem] and [BitPackedKmerSortingDisk], this algorithm holds the
/// concatenation of the input sequences in memory, and is not limited to k <= 256.
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct BuildByBoundedSuffixSort<SS: SeqStream + Send> {
    input: SS,
    prefix_length_for_bucket_sort: usize,

    k: usize,
    n_threads: usize,
    build_lcs: bool,
    add_rev_comp: bool,
    build_select_support: bool,
    precalc_length: usize,
    add_all_dummy_paths: bool,
}

impl<SS: SeqStream + Send> BuildByBoundedSuffixSort<SS> {

    /// Initializes the algorithm for a [crate::SeqStream] with a given k and default settings.
    /// - 4 threads
    /// - no LCS array
    /// - no reverse complements added
    /// - no select support
    /// - precalc length min(8,k)
    /// - add_all_dummy_paths = false
    /// - bucket sort prefix length 4
    pub fn new(input: SS, k: usize) -> Self {
        Self{input, prefix_length_for_bucket_sort: 4, k, n_threads: 4, build_lcs: false, add_rev_comp: false, build_select_support: false, precalc_length: std::cmp::min(8, k), add_all_dummy_paths: false}
    }

    /// Set the length of the prefix used to bucket the contexts before sorting them.
    pub fn prefix_length_for_bucket_sort(mut self, prefix_length: usize) -> Self {
        self.prefix_length_for_bucket_sort = prefix_length;
        self
    }

    /// Sets the k-mer length.
    pub fn k(mut self, k: usize) -> Self {
        self.k = k;
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

    /// Run the construction algorithm, consuming the input stream, and return the SBWT index and
    /// the LCS array if [build_lcs](BuildByBoundedSuffixSort::build_lcs) was set, otherwise `None`.
    ///
    /// The input sequences are first read into a single concatenation in memory, separated by
    /// the `$` character and preceded by the sentinel `#`.
    pub fn run(self) -> (SbwtIndex<SubsetMatrix>, Option<LcsArray>) {
        let input = SeqStreamWithPossiblyRevComp::new(self.input, self.add_rev_comp);
        let mut input = SanitizedReversedSeqStream{ inner: input, buf: Vec::<u8>::new() };

        let mut concatenation = Vec::<u8>::new();
        crate::alternative_construction::preprocessing::concatenate_sequences(&mut input, &mut concatenation).unwrap();

        let crate::alternative_construction::Output{mut sbwt, lcs, counts: _} =
            crate::alternative_construction::build_with_bounded_suffix_array::<SubsetMatrix>(
                self.n_threads,
                concatenation,
                self.k,
                self.prefix_length_for_bucket_sort,
                self.build_lcs,
                self.add_all_dummy_paths,
                false, // Counts are not part of the return value of this interface
            );

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

impl<'a> BuildByBoundedSuffixSort<crate::SliceSeqStream<'a>> {
    /// Initialize the algorithm for a slice of ASCII sequences, with the same defaults as [BuildByBoundedSuffixSort::new].
    pub fn new_from_slices(input: &'a[&'a [u8]], k: usize) -> Self {
        let stream = crate::util::SliceSeqStream::new(input);
        BuildByBoundedSuffixSort::new(stream, k)
    }
}

impl<'a> BuildByBoundedSuffixSort<crate::VecSeqStream<'a>> {
    /// Initialize the algorithm for a slice of owned ASCII sequences, with the same defaults as [BuildByBoundedSuffixSort::new].
    pub fn new_from_vecs(input: &'a[Vec<u8>], k: usize) -> Self {
        let stream = crate::util::VecSeqStream::new(input);
        BuildByBoundedSuffixSort::new(stream, k)
    }
}

impl BuildByBoundedSuffixSort<crate::util::FastXReader> {
    /// Initialize the algorithm for a fasta file, with the same defaults as [BuildByBoundedSuffixSort::new].
    pub fn new_from_fasta<R: std::io::Read + Send + Sync + 'static>(input: R, k: usize) -> Self {
        let input = std::io::BufReader::new(input);
        let input = crate::util::FastXReader{inner: jseqio::reader::DynamicFastXReader::new(input).unwrap()};
        BuildByBoundedSuffixSort::new(input, k)
    }

    /// Initialize the algorithm for a fastq file, with the same defaults as [BuildByBoundedSuffixSort::new].
    pub fn new_from_fastq<R: std::io::Read + Send + Sync + 'static>(input: R, k: usize) -> Self {
        let input = std::io::BufReader::new(input);
        let input = crate::util::FastXReader{inner: jseqio::reader::DynamicFastXReader::new(input).unwrap()};
        BuildByBoundedSuffixSort::new(input, k)
    }
}

/// A construction algorithm that uses the `libsais` crate to build a suffix array and an LCP
/// array of the concatenation of the (reversed) input sequences, and derives the BWT and LCP
/// that [crate::alternative_construction::build_from_input] expects from them. Like
/// [BuildByBoundedSuffixSort], this algorithm holds the concatenation of the input sequences in
/// memory, and is not limited to k <= 256.
///
/// Only available when the optional `libsais` feature is enabled, since it pulls in the
/// `libsais` crate as a dependency. The other construction algorithms are always available.
#[cfg(feature = "libsais")]
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct BuildByLibsais<SS: SeqStream + Send> {
    input: SS,

    k: usize,
    n_threads: usize,
    build_lcs: bool,
    add_rev_comp: bool,
    build_select_support: bool,
    precalc_length: usize,
    add_all_dummy_paths: bool,
}

#[cfg(feature = "libsais")]
impl<SS: SeqStream + Send> BuildByLibsais<SS> {

    /// Initializes the algorithm for a [crate::SeqStream] with a given k and default settings.
    /// - 4 threads
    /// - no LCS array
    /// - no reverse complements added
    /// - no select support
    /// - precalc length min(8,k)
    /// - add_all_dummy_paths = false
    pub fn new(input: SS, k: usize) -> Self {
        Self{input, k, n_threads: 4, build_lcs: false, add_rev_comp: false, build_select_support: false, precalc_length: std::cmp::min(8, k), add_all_dummy_paths: false}
    }

    /// Sets the k-mer length.
    pub fn k(mut self, k: usize) -> Self {
        self.k = k;
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

    /// Set the number of threads for `libsais` to use to construct the suffix array and the LCP array.
    pub fn n_threads(mut self, n_threads: usize) -> Self {
        self.n_threads = n_threads;
        self
    }

    /// Include all dummy paths for every DNA run in the input, not only those strictly required by the SBWT structure.
    pub fn add_all_dummy_paths(mut self, enable: bool) -> Self {
        self.add_all_dummy_paths = enable;
        self
    }

    /// Run the construction algorithm, consuming the input stream, and return the SBWT index and
    /// the LCS array if [build_lcs](BuildByLibsais::build_lcs) was set, otherwise `None`.
    ///
    /// The input sequences are first read into a single concatenation in memory, separated by
    /// the `$` character and preceded by the sentinel `#`. A suffix array and an LCP array of the
    /// concatenation are then built with `libsais`, from which the BWT and LCP that
    /// [crate::alternative_construction::build_from_input] expects are derived.
    pub fn run(self) -> (SbwtIndex<SubsetMatrix>, Option<LcsArray>) {
        let input = SeqStreamWithPossiblyRevComp::new(self.input, self.add_rev_comp);
        let mut input = SanitizedReversedSeqStream{ inner: input, buf: Vec::<u8>::new() };

        let mut concatenation = Vec::<u8>::new();
        crate::alternative_construction::preprocessing::concatenate_sequences(&mut input, &mut concatenation).unwrap();

        // let (bwt_bytes, lcp_bytes) = libsais_bwt_and_lcp(&concatenation, self.n_threads);
        let suffix_array = libsais_suffix_array(&concatenation, self.n_threads);

        let crate::alternative_construction::Output{mut sbwt, lcs, counts: _} =
            crate::alternative_construction::build::<SubsetMatrix>(
                self.n_threads,
                concatenation,
                suffix_array,
                // &mut bwt_bytes.as_slice(),
                // &mut lcp_bytes.as_slice(),
                // concatenation.len(),
                self.k,
                self.build_lcs,
                self.add_all_dummy_paths,
                false, // Counts are not part of the return value of this interface
            );

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

#[cfg(feature = "libsais")]
impl<'a> BuildByLibsais<crate::SliceSeqStream<'a>> {
    /// Initialize the algorithm for a slice of ASCII sequences, with the same defaults as [BuildByLibsais::new].
    pub fn new_from_slices(input: &'a[&'a [u8]], k: usize) -> Self {
        let stream = crate::util::SliceSeqStream::new(input);
        BuildByLibsais::new(stream, k)
    }
}

#[cfg(feature = "libsais")]
impl<'a> BuildByLibsais<crate::VecSeqStream<'a>> {
    /// Initialize the algorithm for a slice of owned ASCII sequences, with the same defaults as [BuildByLibsais::new].
    pub fn new_from_vecs(input: &'a[Vec<u8>], k: usize) -> Self {
        let stream = crate::util::VecSeqStream::new(input);
        BuildByLibsais::new(stream, k)
    }
}

#[cfg(feature = "libsais")]
impl BuildByLibsais<crate::util::FastXReader> {
    /// Initialize the algorithm for a fasta file, with the same defaults as [BuildByLibsais::new].
    pub fn new_from_fasta<R: std::io::Read + Send + Sync + 'static>(input: R, k: usize) -> Self {
        let input = std::io::BufReader::new(input);
        let input = crate::util::FastXReader{inner: jseqio::reader::DynamicFastXReader::new(input).unwrap()};
        BuildByLibsais::new(input, k)
    }

    /// Initialize the algorithm for a fastq file, with the same defaults as [BuildByLibsais::new].
    pub fn new_from_fastq<R: std::io::Read + Send + Sync + 'static>(input: R, k: usize) -> Self {
        let input = std::io::BufReader::new(input);
        let input = crate::util::FastXReader{inner: jseqio::reader::DynamicFastXReader::new(input).unwrap()};
        BuildByLibsais::new(input, k)
    }
}

/// Builds a suffix array of `concatenation` using `libsais`.
///
/// This relies on `'#'`, the sentinel character that
/// [crate::alternative_construction::preprocessing::concatenate_sequences] places at the start
/// of the concatenation, being byte-value-smaller than every other character used
/// (`'#' < '$' < A,C,G,T` in ASCII) and occurring exactly once: `libsais`'s suffix array of
/// `concatenation` (which behaves as if an even smaller, implicit sentinel were appended and
/// then omitted from the output) then coincides exactly, position by position, with the plain
/// byte-value suffix sort that [crate::alternative_construction] is built around, and the BWT can
/// be derived from it with the ordinary `bwt[i] = concatenation[(sa[i] + n - 1) % n]` formula.
#[cfg(feature = "libsais")]
fn libsais_suffix_array(concatenation: &[u8], n_threads: usize) -> Vec<usize> {
    if concatenation.len() <= libsais::LIBSAIS_I32_OUTPUT_MAXIMUM_SIZE {
        run_libsais_suffix_array::<i32>(concatenation, n_threads)
    } else {
        run_libsais_suffix_array::<i64>(concatenation, n_threads)
    }
}

#[cfg(feature = "libsais")]
fn run_libsais_suffix_array<O: libsais::OutputElement>(concatenation: &[u8], n_threads: usize) -> Vec<usize> {
    let sa_construction = libsais::SuffixArrayConstruction::for_text(concatenation).in_owned_buffer::<O>();
    let sa_result = if n_threads > 1 {
        sa_construction.multi_threaded(libsais::ThreadCount::fixed(n_threads as u16)).run()
    } else {
        sa_construction.single_threaded().run()
    }.expect("libsais suffix array construction failed");

    sa_result.into_vec().into_iter().map(|element| element.to_usize().unwrap()).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// All three construction algorithms must produce exactly the same index and LCS array.
    #[test]
    fn all_algorithms_agree() {
        // The input has 540 k-mers in total. With k = 5 there are 4^5 = 1024 distinct k-mers,
        // so repeated k-mers are a certainty, while only about 40% of all possible k-mers
        // occur, which keeps the k-mer set sparse.
        let k = 5;
        let seqs: Vec<Vec<u8>> = (0..8_u64).map(|i| {
            crate::util::gen_random_dna_string(50 + 7 * i as usize, 1234 + i)
        }).collect();

        for add_all_dummy_paths in [false, true] {
            let (mem_sbwt, mem_lcs) = BitPackedKmerSortingMem::new_from_vecs(&seqs, k)
                .build_lcs(true)
                .add_all_dummy_paths(add_all_dummy_paths)
                .run();

            let (disk_sbwt, disk_lcs) = BitPackedKmerSortingDisk::new_from_vecs(&seqs, k)
                .build_lcs(true)
                .add_all_dummy_paths(add_all_dummy_paths)
                .temp_dir(&std::env::temp_dir())
                .run();

            let (bounded_sbwt, bounded_lcs) = BuildByBoundedSuffixSort::new_from_vecs(&seqs, k)
                .build_lcs(true)
                .add_all_dummy_paths(add_all_dummy_paths)
                .run();

            assert_eq!(mem_sbwt, disk_sbwt);
            assert_eq!(mem_sbwt, bounded_sbwt);

            assert!(mem_lcs.is_some());
            assert_eq!(mem_lcs, disk_lcs);
            assert_eq!(mem_lcs, bounded_lcs);

            #[cfg(feature = "libsais")]
            {
                let (libsais_sbwt, libsais_lcs) = BuildByLibsais::new_from_vecs(&seqs, k)
                    .build_lcs(true)
                    .add_all_dummy_paths(add_all_dummy_paths)
                    .run();
                assert_eq!(mem_sbwt, libsais_sbwt);
                assert_eq!(mem_lcs, libsais_lcs);
            }
        }
    }
}

/// Sort and deduplicate the reverse k-mers of the input sequences into a file on disk, to be later
/// consumed by [build_from_kmers_on_disk]. Returns the path to the k-mers file, and if
/// `add_all_dummy_paths` is set, also the path to the file of first k-mers of each input sequence.
/// The returned files are NOT deleted: it is the caller's responsibility to delete them once they
/// are no longer needed.
pub fn sort_and_dedup_kmers_into_file<SS: crate::SubsetSeq + Send, IN: SeqStream + Send>(input: IN, k: usize, mem_gb: usize, n_threads: usize, dedup_batches: bool, add_all_dummy_paths: bool, temp_dir: &Path) -> (std::path::PathBuf, Option<std::path::PathBuf>) {
    let mut tfm = TempFileManager::new(temp_dir);
    let (mut kmers_file, mut first_mers_file) = match k {
        0..=32 => {
            crate::bitpacked_kmer_sorting::sort_and_dedup_kmers_into_file::<1, IN, SS>(input, k, mem_gb, n_threads, dedup_batches, add_all_dummy_paths, &mut tfm)
        }
        33..=64 => {
            crate::bitpacked_kmer_sorting::sort_and_dedup_kmers_into_file::<2, IN, SS>(input, k, mem_gb, n_threads, dedup_batches, add_all_dummy_paths, &mut tfm)
        }
        65..=96 => {
            crate::bitpacked_kmer_sorting::sort_and_dedup_kmers_into_file::<3, IN, SS>(input, k, mem_gb, n_threads, dedup_batches, add_all_dummy_paths, &mut tfm)
        }
        97..=128 => {
            crate::bitpacked_kmer_sorting::sort_and_dedup_kmers_into_file::<4, IN, SS>(input, k, mem_gb, n_threads, dedup_batches, add_all_dummy_paths, &mut tfm)
        }
        129..=256 => {
            crate::bitpacked_kmer_sorting::sort_and_dedup_kmers_into_file::<8, IN, SS>(input, k, mem_gb, n_threads, dedup_batches, add_all_dummy_paths, &mut tfm)
        }
        _ => {
            panic!("k > 256 not supported with bitpacked sorting algorithm.");
        }
    };

    // These files must survive past the end of this function, so tell them not to delete themselves on drop.
    kmers_file.set_delete_on_drop(false);
    if let Some(first_mers_file) = first_mers_file.as_mut() {
        first_mers_file.set_delete_on_drop(false);
    }

    (kmers_file.path.clone(), first_mers_file.map(|f| f.path.clone()))
}

pub fn build_from_kmers_on_disk<SS: crate::SubsetSeq + Send>(k: usize, n_threads: usize, build_lcs: bool, temp_dir: &Path, kmers_file: &Path, first_mers_file: Option<&Path>) -> (SbwtIndex::<SS>, Option<LcsArray>) {
    let mut tfm = TempFileManager::new(temp_dir);
    match k {
        0..=32 => {
            crate::bitpacked_kmer_sorting::build_from_kmers_on_disk::<1, SS>(k, n_threads, build_lcs, &mut tfm, kmers_file, first_mers_file)
        }
        33..=64 => {
            crate::bitpacked_kmer_sorting::build_from_kmers_on_disk::<2, SS>(k, n_threads, build_lcs, &mut tfm, kmers_file, first_mers_file)
        }
        65..=96 => {
            crate::bitpacked_kmer_sorting::build_from_kmers_on_disk::<3, SS>(k, n_threads, build_lcs, &mut tfm, kmers_file, first_mers_file)
        }
        97..=128 => {
            crate::bitpacked_kmer_sorting::build_from_kmers_on_disk::<4, SS>(k, n_threads, build_lcs, &mut tfm, kmers_file, first_mers_file)
        }
        129..=256 => {
            crate::bitpacked_kmer_sorting::build_from_kmers_on_disk::<8, SS>(k, n_threads, build_lcs, &mut tfm, kmers_file, first_mers_file)
        }
        _ => {
            panic!("k > 256 not supported with bitpacked sorting algorithm.");
        }
    }
}
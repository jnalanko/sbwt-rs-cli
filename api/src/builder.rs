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
        let input = SeqStreamWithPossiblyRevComp{ inner: self.input, rc_buf: Vec::<u8>::new(), parity: false, enable_rev_comp: self.add_rev_comp };
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

impl BitPackedKmerSortingDisk<JSeqIOSeqStreamWrapper> {
    /// Initialize the algorithm for a fasta file, with the same defaults as [BitPackedKmerSortingDisk::new].
    pub fn new_from_fasta<R: std::io::Read + Send + Sync + 'static>(input: R, k: usize) -> Self {
        let input = std::io::BufReader::new(input);
        let input = crate::JSeqIOSeqStreamWrapper{inner: jseqio::reader::DynamicFastXReader::new(input).unwrap()};
        BitPackedKmerSortingDisk::new(input, k)
    }

    /// Initialize the algorithm for a fastq file, with the same defaults as [BitPackedKmerSortingDisk::new].
    pub fn new_from_fastq<R: std::io::Read + Send + Sync + 'static>(input: R, k: usize) -> Self {
        let input = std::io::BufReader::new(input);
        let input = crate::JSeqIOSeqStreamWrapper{inner: jseqio::reader::DynamicFastXReader::new(input).unwrap()};
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
        let input = SeqStreamWithPossiblyRevComp{ inner: self.input, rc_buf: Vec::<u8>::new(), parity: false, enable_rev_comp: self.add_rev_comp };
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

impl BitPackedKmerSortingMem<JSeqIOSeqStreamWrapper> {
    /// Initialize the algorithm for a fasta file, with the same defaults as [BitPackedKmerSortingMem::new].
    pub fn new_from_fasta<R: std::io::Read + Send + Sync + 'static>(input: R, k: usize) -> Self {
        let input = std::io::BufReader::new(input);
        let input = crate::JSeqIOSeqStreamWrapper{inner: jseqio::reader::DynamicFastXReader::new(input).unwrap()};
        BitPackedKmerSortingMem::new(input, k)
    }

    /// Initialize the algorithm for a fastq file, with the same defaults as [BitPackedKmerSortingMem::new].
    pub fn new_from_fastq<R: std::io::Read + Send + Sync + 'static>(input: R, k: usize) -> Self {
        let input = std::io::BufReader::new(input);
        let input = crate::JSeqIOSeqStreamWrapper{inner: jseqio::reader::DynamicFastXReader::new(input).unwrap()};
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
        let input = SeqStreamWithPossiblyRevComp{ inner: self.input, rc_buf: Vec::<u8>::new(), parity: false, enable_rev_comp: self.add_rev_comp };
        let mut input = SanitizedReversedSeqStream{ inner: input, buf: Vec::<u8>::new() };

        let mut concatenation = Vec::<u8>::new();
        crate::alternative_construction::preprocessing::concatenate_sequences(&mut input, &mut concatenation).unwrap();

        let crate::alternative_construction::Output{mut sbwt, lcs, counts: _} =
            crate::bounded_alternative_construction::build::<SubsetMatrix>(
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

impl BuildByBoundedSuffixSort<JSeqIOSeqStreamWrapper> {
    /// Initialize the algorithm for a fasta file, with the same defaults as [BuildByBoundedSuffixSort::new].
    pub fn new_from_fasta<R: std::io::Read + Send + Sync + 'static>(input: R, k: usize) -> Self {
        let input = std::io::BufReader::new(input);
        let input = crate::JSeqIOSeqStreamWrapper{inner: jseqio::reader::DynamicFastXReader::new(input).unwrap()};
        BuildByBoundedSuffixSort::new(input, k)
    }

    /// Initialize the algorithm for a fastq file, with the same defaults as [BuildByBoundedSuffixSort::new].
    pub fn new_from_fastq<R: std::io::Read + Send + Sync + 'static>(input: R, k: usize) -> Self {
        let input = std::io::BufReader::new(input);
        let input = crate::JSeqIOSeqStreamWrapper{inner: jseqio::reader::DynamicFastXReader::new(input).unwrap()};
        BuildByBoundedSuffixSort::new(input, k)
    }
}

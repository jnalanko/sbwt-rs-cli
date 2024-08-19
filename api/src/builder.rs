//! A builder pattern interface for building an [SbwtIndex].

use std::str::FromStr;

use crate::{subsetseq::SubsetMatrix, SeqStream};
use crate::sbwt::{PrefixLookupTable, SbwtIndex};
use crate::streaming_index::LcsArray;

// 'a must be a higher-ranked lifetime because it's tied to the lifetime
// of the borrow in stream_next? That function should be callabed with
// *any* borrow length, so that's why a regular generic lifetime parameter
// does not compile?
#[derive(Clone, Eq, PartialEq, Debug)]
struct SeqStreamWithRevComp<SS: SeqStream + Send>{
    inner: SS, 
    rc_buf: Vec<u8>,
    parity: bool, // Every other sequence we return is a reverse complement of the previous. Initialize to false.
}

impl<SS: SeqStream + Send> crate::SeqStream for SeqStreamWithRevComp<SS>{
    fn stream_next(&mut self) -> Option<&[u8]> {
        self.parity = !self.parity;

        if self.parity {
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
    fn run<SS: SeqStream + Send>(self, input: SS, k: usize, n_threads: usize, build_lcs: bool) -> (SbwtIndex<SubsetMatrix>, Option<LcsArray>);
}

/// A construction algorithm based on sorting of bit-packed k-mers.
#[derive(Default)]
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct BitPackedKmerSorting{
    mem_gb: usize,
    dedup_batches: bool,
    temp_dir: std::path::PathBuf,
}

impl BitPackedKmerSorting {

    /// Initializes the algorithm with default settings:
    /// - 4 GB of memory.
    /// - do not deduplicate k-mer batches before sorting.
    /// - use the current directory as the temporary directory.
    pub fn new() -> Self {
        Self{mem_gb: 4, dedup_batches: false, temp_dir: std::path::PathBuf::from_str(".").unwrap()}
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

impl SbwtConstructionAlgorithm for BitPackedKmerSorting {
    fn run<SS: SeqStream + Send>(self, input: SS, k: usize, n_threads: usize, build_lcs: bool) -> (SbwtIndex<SubsetMatrix>, Option<LcsArray>) {
        let mem_gb = self.mem_gb;
        let dedup_batches = self.dedup_batches;
        let mut temp_file_manager = crate::tempfile::TempFileManager::new(&self.temp_dir);
        match k {
            0..=32 => {
                crate::bitpacked_kmer_sorting::build_with_bitpacked_kmer_sorting::<1,_,SubsetMatrix>(input, k, mem_gb, n_threads, dedup_batches, build_lcs, &mut temp_file_manager)
            }
            33..=64 => {
                crate::bitpacked_kmer_sorting::build_with_bitpacked_kmer_sorting::<2,_,SubsetMatrix>(input, k, mem_gb, n_threads, dedup_batches, build_lcs, &mut temp_file_manager)
            }
            65..=96 => {
                crate::bitpacked_kmer_sorting::build_with_bitpacked_kmer_sorting::<3,_,SubsetMatrix>(input, k, mem_gb, n_threads, dedup_batches, build_lcs, &mut temp_file_manager)
            }
            97..=128 => {
                crate::bitpacked_kmer_sorting::build_with_bitpacked_kmer_sorting::<4,_,SubsetMatrix>(input, k, mem_gb, n_threads, dedup_batches, build_lcs, &mut temp_file_manager)
            }
            129..=256 => {
                crate::bitpacked_kmer_sorting::build_with_bitpacked_kmer_sorting::<8,_,SubsetMatrix>(input, k, mem_gb, n_threads, dedup_batches, build_lcs, &mut temp_file_manager)
            }
            _ => {
                panic!("k > 256 not supported with bitpacked sorting algorithm.");
            }
        } 
    }
}

/// A builder for constructing an SBWT index.
#[derive(Default)]
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct SbwtIndexBuilder<A: SbwtConstructionAlgorithm + Default> {
    k: usize,
    n_threads: usize,
    algorithm: A,
    build_lcs: bool,
    add_rev_comp: bool,
    build_select_support: bool,
    precalc_length: usize,
}

impl<A: SbwtConstructionAlgorithm + Default> SbwtIndexBuilder<A> {

    /// Sets up the builder with default values:
    /// - k = 31.
    /// - n_threads = 4.
    /// - do not build the LCS array.
    /// - do not add the reverse complement of the input sequences.
    /// - do not build the select support.
    /// - precalc_length = 8.
    /// - default settings for the chosen algorithm.
    pub fn new() -> Self {
        Self{k: 31, n_threads: 4, algorithm: A::default(), build_lcs: false, add_rev_comp: false, build_select_support: false, precalc_length: 8}
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

    /// Run the algorithm and return the SBWT index and optionally the LCS array if [build_lcs](SbwtIndexBuilder::build_lcs) was set.
    /// See also [run_from_slices](SbwtIndexBuilder::run_from_slices), [run_from_fasta](SbwtIndexBuilder::run_from_fasta) and [run_from_fastq](SbwtIndexBuilder::run_from_fastq).
    pub fn run<SS: SeqStream + Send>(self, input: SS) -> (SbwtIndex<SubsetMatrix>, Option<LcsArray>) {
        let (mut sbwt, lcs) = if self.add_rev_comp {
            let input_with_rc = SeqStreamWithRevComp{inner: input, rc_buf: Vec::<u8>::new(), parity: false}; 
            self.algorithm.run(input_with_rc, self.k, self.n_threads, self.build_lcs)
        } else {
            self.algorithm.run(input, self.k, self.n_threads, self.build_lcs)
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


    /// Run the algorithm for the given ASCII nucleotide strings, and return the SBWT index and optionally the LCS array if [SbwtIndexBuilder::build_lcs] was set.
    pub fn run_from_slices(self, input: &[&[u8]]) -> (SbwtIndex<SubsetMatrix>, Option<LcsArray>) {
        let stream = crate::util::SliceSeqStream::new(input);
        self.run(stream)
    } 

    /// Run the algorithm for the given ASCII nucleotide strings, and return the SBWT index and optionally the LCS array if [SbwtIndexBuilder::build_lcs] was set.
    pub fn run_from_vecs(self, input: &[Vec<u8>]) -> (SbwtIndex<SubsetMatrix>, Option<LcsArray>) {
        let stream = crate::util::VecSeqStream::new(input);
        self.run(stream)
    } 

    /// Run the algorithm from a FASTA-formatted stream of sequences, and return the SBWT index and optionally the LCS array if [SbwtIndexBuilder::build_lcs] was set.
    pub fn run_from_fasta<R: std::io::Read + Send + 'static>(self, input: R) -> (SbwtIndex<SubsetMatrix>, Option<LcsArray>) {
        let input = std::io::BufReader::new(input);
        let input = crate::JSeqIOSeqStreamWrapper{inner: jseqio::reader::DynamicFastXReader::new(input).unwrap()};
        self.run(input)
    } 

    /// Run the algorithm from a FASTQ-formatted stream of sequences, and return the SBWT index and optionally the LCS array if [SbwtIndexBuilder::build_lcs] was set.
    pub fn run_from_fastq<R: std::io::Read + Send + 'static>(self, input: R) -> (SbwtIndex<SubsetMatrix>, Option<LcsArray>) {
        let input = std::io::BufReader::new(input);
        let input = crate::JSeqIOSeqStreamWrapper{inner: jseqio::reader::DynamicFastXReader::new(input).unwrap()};
        self.run(input)
    } 
}

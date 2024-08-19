//! # Introduction
//! 
//! This crate contains an implementation of the 
//! [Bit Matrix SBWT data structure](`SbwtIndex<SubsetMatrix>`), as described in 
//! [Small Searchable k-Spectra via Subset Rank Queries on the Spectral Burrows-Wheeler Transform](https://epubs.siam.org/doi/abs/10.1137/1.9781611977714.20),
//! for the DNA alphabet ACGT. 
//! The data structure uses a variant of the Burrows-Wheeler transform to compress a set of k-mers in way that allows fast lookup queries. 
//! If the input k-mers are consecutive k-mers from longer underlying sequences,
//! the index takes typically around 5 bits per distinct k-mer, supporting k-mer lookup queries at a speed of around
//! 1 μs / k-mer on modern hardware.
//! 
//! Queries can be further sped up by using the [Longest common suffix array](LcsArray), (see [here](https://link.springer.com/chapter/10.1007/978-3-031-43980-3_1))
//! taking roughly log(k) bits of space per k-mer. The LCS array also enables the 
//! computation of [*k-bounded matching statistics*][streaming_index::StreamingIndex::matching_statistics] 
//! and [*shortest frequency-bounded suffixes*][streaming_index::StreamingIndex::shortest_freq_bound_suffixes]
//! (both described [here](https://www.biorxiv.org/content/10.1101/2024.02.19.580943v1)) 
//! Finally, the crate provides an [interface][dbg::Dbg] for traversing the node-centric de Bruijn graph
//! of the k-mers.
//! 
//! # API Quick start 
//! 
//! ```
//! use sbwt::*;
//! use std::io::BufReader;
//! use std::io::BufWriter;
//! use std::fs::File;
//! use std::path::Path;
//! 
//! // Build the sbwt
//! let seqs: Vec<&[u8]> = vec![b"AACTGACTGATCGTCTTGACTCGTTTATCTACGGT", b"ACTGACAGCTCTGCGATGCGA"];
//! let seq_stream = sbwt::SliceSeqStream::new(seqs.as_slice());
//! let (sbwt, lcs) = SbwtIndexBuilder::new()
//!     .k(6).n_threads(4).build_lcs(true).add_rev_comp(true)
//!     .algorithm(BitPackedKmerSorting::new()
//!         .mem_gb(2)
//!         .dedup_batches(false)
//!         .temp_dir(Path::new("./temp")))
//!     .run(seq_stream);
//!
//! // Query a k-mer
//! let query_kmer = b"GACTCG";
//! if sbwt.search(query_kmer).is_some() {
//!     println!("{} is found", &String::from_utf8_lossy(query_kmer));
//! }
//!
//! // Query all k-mers of a longer query using k-bounded matching statistics
//! let lcs = lcs.unwrap(); // Ok because we used build_lcs(true)
//! let streaming_index = StreamingIndex::new(&sbwt, &lcs);
//! let long_query = b"TGATACGTCTTAGTGACTCGTTT";
//! for (i, (len, range)) in streaming_index.matching_statistics(long_query).iter().enumerate() {
//!     // Kmer ending at long_query[i] exists iff len == k
//!     println!("Longest match ending at {} has length {} and colex range {:?}", i, len, range);
//! }
//!
//! // Write index to disk for later use
//! sbwt.serialize(&mut BufWriter::new(File::create("index.sbwt").unwrap())).unwrap();
//! lcs.serialize(&mut BufWriter::new(File::create("index.lcs").unwrap())).unwrap();
//! 
//! ```
//! 
//! # Reverse complements and input preprocessing
//! 
//! If your input is very large, you may want to preprocess it to remove duplicate k-mers reduce construction time, memory and disk space.
//! We recommend using a specialized tool such as
//! [GGCAT](https://github.com/algbio/ggcat/), 
//! [Cuttlefish](https://github.com/COMBINE-lab/cuttlefish) or 
//! [BCALM2](https://github.com/GATB/bcalm) for this purpose.
//! 
//! This crate considers a k-mer distinct from
//! its reverse complement. In the use case where a k-mer is considered equal to its reverse
//! complement, you need to either feed both directions to the index, or feed
//! just one direction but query each k-mer both ways. The former approach can be easily implemented by
//! [enabling reverse complements](builder::SbwtIndexBuilder::add_rev_comp) in the builder.
//! For the latter approach, we recommend using one of the unitig construction tools mentioned above 
//! to turn the data into canonical unitigs containing each k-mer in only one orientation. After this, it is
//! **very important** to [reorient the unitigs](optimize_unitig_orientation) to maintain a consistent orientation for neighboring k-mers
//! since as explained [here](#details-on-the-space-usage-of-the-index), each
//! k-mer without a direct predecessor increases the index size. 
//! The payoff of this approach is up to two times smaller index size compared to indexing both orientations explicitly, with the drawback that now queries
//! need to be run in both orientations.
//! 
//! # De Bruijn graph operations
//! 
//! With the addition of two extra bitvectors, the data structure supports traversal on the node-centric de-Bruijn graph of the input k-mers. The set
//! of nodes is the set of distinct k-mers, and there is an edge from x to y iff x[1..k) = y[0..k-1).
//! The label of the edge is the last character of y. 
//! The graph is not aware of reverse complements, so if you want to traverse the bi-directed de Bruijn
//! graph that includes edges where one or both of the endpoints reverse complemented, you need to take care of
//! that logic yourself. See [Dbg][`dbg::Dbg`] for the API.
//! 
//! # Construction algorithms 
//! 
//! The crate currently provides one construction algorithm: [bitpacked k-mer sorting][`builder::BitPackedKmerSorting`]. 
//! This extracts all k-mers in the input, packs each k-mer to 2k bits, and sorts them in parallel. This requires O(nk) time and disk
//! space, so it is suitable only for small k. There is still a lot of room for optimization in the implementation. 
//! For larger k, a suffix-array-based construction algorithm is planned.
//! 
//! # Details on the space usage of the index
//! 
//! The index exploits overlaps between k-mers to encode them in small space.
//! We say that a k-mer x is a *source k-mer* if it has no incoming edges in the
//! node-centric de Bruijn graph of the input k-mers S, that is, there does not
//! exist a k-mer y ∈ S such that y[1..k) = x[0..k-1).
//! The number of bits in [`SbwtIndex<SubsetMatrix>`] is 5(n + n') plus a small constant, 
//! where n is the number
//! of distinct k-mers in the dataset, and n' is the number of nodes in the trie
//! of all prefixes of length k-1 of all source k-mers (See [here](SbwtIndex#sbwt-graph) for more
//! details on the inner workings of the SBWT to understand what is going on). When the k-mers are from
//! biological sequences, and unitigs are [oriented consistently](optimize_unitig_orientation), the term n' is typically negligible, but if the k-mers
//! are for example a randomly sampled subset, then the benefit of overlaps is
//! lost, and the term n' dominates. In the worst case, n' can be up to n(k-1) + 1.
//! 
//! # Limitations
//! 
//! The implementation only supports the DNA alphabet ACGT. For best compression, the input k-mers
//! should originate from a longer underlying sequence so that sbwt is able to exploit the
//! overlaps for better compression. For non-overlapping k-mer sets, a simple hash table is likely
//! a better choice.
//! 
//! 

// We're using standard upper-case names in string data structures, like SA and BWT
#![allow(non_snake_case)]

// String algorithms are often clearer with explicit indexing
#![allow(clippy::needless_range_loop)]

mod bitpacked_kmer_sorting;
mod tempfile;
mod util;

pub mod dbg;
pub mod benchmark;

mod sbwt;
pub use sbwt::*;

mod builder;
pub use builder::*;

mod streaming_index;
pub use streaming_index::*;

mod subsetseq;
pub use subsetseq::*;

pub use unitig_flipper::Orientation;
pub use util::reverse_complement_in_place;
pub use util::VecSeqStream;
pub use util::SliceSeqStream;

/// A stream of ASCII-encoded DNA-sequences. This is not necessarily a standard Rust iterator
/// because we want to support streaming sequences from disk, which is not possible
/// with a regular iterator due to lifetime constraints of the Iterator trait.
pub trait SeqStream{
    fn stream_next(&mut self) -> Option<&[u8]>;
}

pub(crate) struct UnitigFlipperSeqStreamWrapper<SS: SeqStream> {
    inner: SS,
}

impl<'a, SS: crate::SeqStream> unitig_flipper::SeqStream<'a> for UnitigFlipperSeqStreamWrapper<SS> {
    fn stream_next<'b>(&'b mut self) -> Option<&'b [u8]> where 'a : 'b {
        self.inner.stream_next() 
    }
}

pub(crate) struct JSeqIOSeqStreamWrapper {
    pub(crate) inner: jseqio::reader::DynamicFastXReader,
}

impl crate::SeqStream for JSeqIOSeqStreamWrapper {
    fn stream_next(&mut self) -> Option<&[u8]> {
        self.inner.read_next().unwrap().map(|rec| rec.seq)
    }
}


/// Given an iterator of canonical unitigs of the node-centric de Bruijn graph of order k, 
/// returns a vector of orientations, one for each sequence, aiming to
/// minimize the number of unitigs which do not have an incoming edge in the de Bruijn graph
/// of order k. Canonical unitigs means that each k-mer occur in only one orientation.
/// The input can be a [VecSeqStream], a [SliceSeqStream] or any struct implementing [SeqStream].
/// 
/// Usage example:
/// 
/// ```
/// let mut seqs: Vec<Vec<u8>> = vec![b"ACGTA".to_vec(), b"TAAGAT".to_vec(), b"CAGTC".to_vec()];
/// let mut stream = sbwt::VecSeqStream::new(&seqs);
/// let k = 4;
/// let orientations = sbwt::optimize_unitig_orientation(stream, k);
/// for (i, &ori) in orientations.iter().enumerate(){
///     if ori == sbwt::Orientation::Reverse {
///         sbwt::reverse_complement_in_place(&mut seqs[i]);
///     }
/// }
/// ```
/// 
pub fn optimize_unitig_orientation<SS: crate::SeqStream>(input: SS, k: usize) -> Vec<Orientation> {
    let stream = UnitigFlipperSeqStreamWrapper{inner: input};
    unitig_flipper::optimize_unitig_orientation(stream, k)
}

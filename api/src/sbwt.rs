//! The [SbwtIndex] data structure. Construct with [SbwtIndexBuilder](crate::SbwtIndexBuilder).

use std::cmp::min;
use std::io::Read;
use std::io::Write;
use bitvec::slice::BitSlice;
use bitvec::vec::BitVec;
use bitvec::prelude::*;
use byteorder::ReadBytesExt;

use byteorder::LittleEndian;
use num::traits::ToBytes;
use rayon::iter::IntoParallelIterator;
use rayon::iter::ParallelIterator;

use crate::sdsl_compatibility::load_known_width_sdsl_int_vector;
use crate::sdsl_compatibility::load_sdsl_bit_vector;
use crate::subsetseq::*;
use crate::util;
use crate::util::ACGT_TO_0123;
use crate::util::DNA_ALPHABET;

pub const CARGO_API_VERSION: &str = env!("CARGO_PKG_VERSION");

/// The SBWT index data structure. Construct with [SbwtIndexBuilder](crate::SbwtIndexBuilder). For the [SubsetSeq] trait implementation,  we recommend using the bit matrix implementation [SubsetMatrix]. 
///
/// # SBWT index 
/// 
/// The SBWT index is a compressed index for searching for k-mers in a set of k-mers. It can be seen a version of the
/// FM-index on sets of k-mers. Here we give a brief overview of the data structure and key concepts that are
/// helpful for understanding the API. For further details, see the [paper](https://doi.org/10.1137/1.9781611977714.20).
/// 
/// ## SBWT graph
/// 
/// To understand the SBWT index, it's helpful to first understand the SBWT graph.
/// The node-centric de Bruijn graph of a set of k-mers R (also called *a k-spectrum*) is a directed graph where the nodes are the k-mers in R,
/// and there is an edge from x to y if x[1..k) = y[0..k-1). The label of the edge is the last character of y.
/// The **SBWT graph** is a modified version of the node-centric de Bruijn graph. There are two modifications:
/// 1. We add a set R' of "dummy nodes". First, the **source set** R' ⊆ R is the subset of k-mers in R that 
///    do not have an incoming edge in the de Bruijn graph, that is x ∈ R' iff x[0..k-1) is not the suffix of 
///    any k-mer in R. The set of dummy nodes is the set of all proper prefixes of the k-mers in R'. We pad the 
///    dummy nodes with dollar symbols from the left so that all dummy nodes have length k. The set R ∪ R' is called 
///    the **padded k-spectrum**. We add all de Bruijn graph edges between overlaps of the padded k-spectrum, that is,
///    for every x ∈ R ∪ R' and y ∈ R ∪ R' such that x[1..k) = y[0..k-1), we add an edge from x to y in the SBWT graph.
///    As a special case, the empty string $^k is always included in R', and the incoming self-loop edge labeled with $ is not included.
/// 2. For every node in R ∪ R' that has more than one incoming edge, we delete all incoming edges except the 
///    one that comes from the colexicographically smallest k-mer. A k-mer x is colexicographically smaller than k-mer y iff the reverse of x is lexicographically smaller than the reverse of y.
/// 
/// 
/// For example, below is the SBWT graph of all 4-mers of strings ["TGTTTG", "TTGCTAT", "ACGTAGTATAT", "TGTAAA"].
/// The source node set is shown in blue, the dummy node set R' in orange, and the de Bruijn graph edges that 
/// have been removed are shown in red.
/// 
/// ![SBWT graph][sbwt_graph] 
/// 
/// ## SBWT definition
/// 
/// The SBWT is a sequence of subsets of {A,C,G,T} such that the i-th subset contains the edge labels of outgoing
/// edges from the i-th k-mer in the SBWT graph in colexicographic order of the k-mers. 
/// If we take the example from above and stack the nodes in colexicographic order, we get: 
///
/// ![SBWT graph][sbwt_sequence] 
/// 
/// The SBWT subset sequence in this example is the sequence of outgoing edge label sets read from top to bottom:
/// {A,T}, {C}, {}, {A}, {T}, {T}, {A,G,T}, {}, {}, {G}, {T}, {T}, {T}, {T}, {C}, {G}, {A}, {}, {}, {A}, {A}, {A}, {A,T}, {T}, {G}.
/// The key property that makes the SBWT index work is that the i-th outgoing edge with a given label on the right column
/// in the picture is the same edge as the i-th incoming edge with the same label on the left column.
/// Thanks to this property, once the subset sequnce has been constructed, the k-mers
/// and the graph can be discarded, with no loss of information. That is, the k-mers and the graph can be reconstructed from the
/// subset sequence alone.
/// 
/// Given the subset sequence, there is also an algorithm to search for a k-mer in O(k) time. The algorithm is essentially the same
/// as the search algorithm on the FM-index. Given a k-mer P[0..k], at iteration i, we have the range of rows in the picture whose
/// k-mers have P[0..i) as a suffix. When i = 0, we have the empty suffix P[0..0), which is a suffix of all k-mers.
/// To update the range for the next iteration, we follow the edges labeled with the next character of the pattern, using the key property
/// mentioned above, as shown in the figure below for query k-mer TATA. This update step can be implemented efficiently by preprocessing the SBWT set sequence for
/// *subset rank queries*. See [this paper](https://doi.org/10.1137/1.9781611977714.20) for more on how to implement these rank queries. We provide
/// an implementation based on a bit matrix at [SubsetMatrix]. Any struct implementing the [SubsetSeq] trait can be used to query the SBWT.
///  
/// ![SBWT graph][sbwt_search] 
/// 
/// 
#[embed_doc_image::embed_doc_image("sbwt_graph", "doc_images/sbwt_graph.svg")] 
#[embed_doc_image::embed_doc_image("sbwt_sequence", "doc_images/sbwt_figure.drawio.svg")] 
#[embed_doc_image::embed_doc_image("sbwt_search", "doc_images/sbwt_figure_with_search.drawio.png")] // This is as .png because there is a bug in vertical centering in the svg export of drawio.

#[derive(Clone, Eq, PartialEq, Debug)]
#[allow(non_snake_case)] // C-array is an established convention in BWT indexes
pub struct SbwtIndex<SS: SubsetSeq> {
    pub(crate) sbwt: SS, // pub(crate) for testing from submodules
    n_kmers: usize,
    k: usize,
    C: Vec<usize>, // Cumulative character counts (includes one ghost dollar)
    prefix_lookup_table: PrefixLookupTable,
}

/// An enum listing SbwtIndex types built on different subset rank implementations provided in this crate. 
/// Currently, only the [SubsetMatrix] structure is supported, but in the future there might be more.
pub enum SbwtIndexVariant {
    SubsetMatrix(SbwtIndex<SubsetMatrix>),
    //SubsetConcat(SbwtIndex<SubsetConcat>), // Remember to update the doc comment if a variant is added
}

/// Loads an index that is wrapped in an enum describing the used subset rank structure type. 
/// The format includes a type identifier so the correct variant can later be loaded with [load_sbwt_index_variant].
pub fn write_sbwt_index_variant(sbwt: &SbwtIndexVariant, out: &mut impl std::io::Write) -> std::io::Result<usize> {
    match sbwt {
        SbwtIndexVariant::SubsetMatrix(sbwt) => {
            out.write_all(&(b"SubsetMatrix".len() as u64).to_le_bytes())?;
            out.write_all(b"SubsetMatrix")?;
            sbwt.serialize(out)
        }
    }
}

/// Loads an index that was stored with [write_sbwt_index_variant]. This includes a type identifier
/// to load the correct subset rank variant.
pub fn load_sbwt_index_variant(input: &mut impl std::io::Read) -> Result<SbwtIndexVariant, Box<dyn std::error::Error>> {
    let type_id_len = input.read_u64::<LittleEndian>().unwrap();
    let mut type_id = vec![0_u8; type_id_len as usize];
    input.read_exact(&mut type_id)?;

    if type_id == b"SubsetMatrix" {
        Ok(SbwtIndexVariant::SubsetMatrix(SbwtIndex::<SubsetMatrix>::load(input)?))
    } else {
        Err("Unknown SBWT index type".into())
    }

}

/// Loads an index that was previously serialized either with the C++ API https://github.com/algbio/SBWT,
/// or the associated CLI. Supports only version v0.1.
#[allow(non_snake_case)] // For C-array
pub fn load_from_cpp_plain_matrix_format<R: std::io::Read>(input: &mut R) -> std::io::Result<SbwtIndexVariant> {

    // The file format with the CLI is:
    // [type id len] [type id] [version string id] [version string] [data]

    // For example:
    // 00000000  0c 00 00 00 00 00 00 00  70 6c 61 69 6e 2d 6d 61  |........plain-ma|
    // 00000010  74 72 69 78 04 00 00 00  00 00 00 00 76 30 2e 31  |trix........v0.1|
    // 00000020  b6 09 2c 0f 00 00 00 00  ff ff ff 8b 48 44 2a 95  |..,.........HD*.|
    
    // We only support type-id "plain-matrix" and version v0.1.

    // The file format of the raw API is:
    // [version string id] [version string] [data]

    // We can tell which format we have by checking what the first encoded string is:
    // - If it's the version string, we check it and proceed to load the sbwt.
    // - If it's plain-matrix, we ignore it, then read the next string, and check that it's the correct version.
    // - If it's a string other than plain-matrix or a version string, we throw an error

    let SUPPORTED_CPP_SBWT_VARIANT = b"plain-matrix"; // If you change this, remember to update all the comments above.
    let SUPPORTED_CPP_FILE_FORMAT = b"v0.1"; // If you change this, remember to update all the comments above, including the doc comment.

    let mut string_length = input.read_u64::<LittleEndian>().unwrap();
    if string_length > 100 { // The first string is never this long
        return Err(std::io::ErrorKind::InvalidData.into());
    }
    let mut string_buf = vec![0u8; string_length as usize];
    input.read_exact(&mut string_buf).unwrap();

    if string_buf == SUPPORTED_CPP_SBWT_VARIANT { 
        // Ok. Ignore and read the next string.

        string_length = input.read_u64::<LittleEndian>().unwrap();
        if string_length > 100 { // The second string is never this long
            return Err(std::io::ErrorKind::InvalidData.into());
        }
        string_buf = vec![0u8; string_length as usize];
        input.read_exact(&mut string_buf).unwrap();
    }

    if string_buf == SUPPORTED_CPP_FILE_FORMAT { 
        // Ok. Proceed to deserialize. 

        log::info!("Reading bit vectors");
        let A_bits = simple_sds_sbwt::bit_vector::BitVector::from(load_sdsl_bit_vector(input)?);
        let C_bits = simple_sds_sbwt::bit_vector::BitVector::from(load_sdsl_bit_vector(input)?);
        let G_bits = simple_sds_sbwt::bit_vector::BitVector::from(load_sdsl_bit_vector(input)?);
        let T_bits = simple_sds_sbwt::bit_vector::BitVector::from(load_sdsl_bit_vector(input)?);

        let _A_rank_support = load_known_width_sdsl_int_vector(input, 64)?; // Ignore: we build our own
        let _C_rank_support = load_known_width_sdsl_int_vector(input, 64)?; // Ignore: we build our own
        let _G_rank_support = load_known_width_sdsl_int_vector(input, 64)?; // Ignore: we build our own
        let _T_rank_support = load_known_width_sdsl_int_vector(input, 64)?; // Ignore: we build our own 

        let _suffix_group_starts = load_sdsl_bit_vector(input); // Ignore: we don't use this

        log::info!("Reading the C-array");
        let C_array_byte_length = input.read_u64::<LittleEndian>()? as usize;
        assert!(C_array_byte_length % 8 == 0);
        let C_array_len = C_array_byte_length / 8;
        let mut C_array = Vec::<usize>::with_capacity(C_array_len);
        for _ in 0..C_array_len {
            C_array.push(input.read_u64::<LittleEndian>()? as usize);
        }

        log::info!("Reading the precalc table");
        let precalc_array_byte_length = input.read_u64::<LittleEndian>()?;
        assert!(precalc_array_byte_length % 16 == 0);
        let mut precalc_array_bytes = vec![0u8; precalc_array_byte_length as usize];
        input.read_exact(precalc_array_bytes.as_mut_slice())?;

        log::info!("Reading constants");
        let precalc_k = input.read_u64::<LittleEndian>()? as usize;
        let n_nodes = input.read_u64::<LittleEndian>()? as usize;
        let n_kmers = input.read_u64::<LittleEndian>()? as usize;
        let k = input.read_u64::<LittleEndian>()? as usize;

        log::info!("Permuting the precalc table");
        // For some stupid reason the lookup table in the Rust code stores the ranges in the
        // lexicographic order of the prefixes, not colexicographic like in the C++ code, which
        // also would make more sense since then the intervals are non-decreasing. We could change
        // to colex order in the Rust code, but that would break index compatibility. So we permute
        // the lookup table here to the lexicographic order.
        let n_ranges = (precalc_array_byte_length/16) as usize;
        let mut ranges = vec![0..0; n_ranges];
        let mut buf = [0_u8; 8];
        for i in 0..n_ranges {

            buf.copy_from_slice(&precalc_array_bytes[(2*i)*8..(2*i+1)*8]);
            let left = i64::from_le_bytes(buf);
            buf.copy_from_slice(&precalc_array_bytes[(2*i+1)*8..(2*i+2)*8]);
            let right = i64::from_le_bytes(buf);

            let mut lut_idx = 0; // Index in the new lookup table
            for j in 0..precalc_k {
                lut_idx = (lut_idx << 2) | ((i >> (2*j)) & 0x3); // Trust me
            }

            ranges[lut_idx] = if left == -1 {
                0..0
            } else {
                left as usize .. (right+1) as usize // +1 to make the end exclusive
            }
        }
        if precalc_k == 0 {
            ranges.push(0..n_nodes); // Special case: range of the empty string
        }
        drop(precalc_array_bytes); // Free some memory

        log::info!("Building rank structures");
        let mut subset_rank = SubsetMatrix::new_from_bit_vectors(vec![A_bits, C_bits, G_bits, T_bits]);
        subset_rank.build_rank();

        let prefix_lut = PrefixLookupTable{ranges, prefix_length: precalc_k};

        Ok(SbwtIndexVariant::SubsetMatrix(SbwtIndex::from_components(subset_rank, n_kmers, k, C_array, prefix_lut)))

    } else {
        Err(std::io::ErrorKind::InvalidData.into())
    }

}


// When non-compatible changes to the serialization format occur, update the version number here to the current version 
const SERIALIZATION_MAGIC_STRING: &[u8] = b"sbwtfile-v0.3.0"; 

impl<SS: SubsetSeq> SbwtIndex<SS> {

    /// Number of k-mers in the index, not counting dummy k-mers.
    pub fn n_kmers(&self) -> usize {
        self.n_kmers
    }

    /// Number of sets in the SBWT.
    pub fn n_sets(&self) -> usize {
        self.sbwt.len()
    }

    /// Length of the k-mers in the index.
    pub fn k(&self) -> usize {
        self.k
    }

    /// Maps A,C,G,T to 0,1,2,3.
    /// Panics if c is not a DNA character.
    pub fn char_idx(&self, c: u8) -> usize {
        let idx = ACGT_TO_0123[c as usize] as usize;
        if idx > DNA_ALPHABET.len() {
            panic!("Invalid character: {}", char::from(c));
        }
        idx
    }

    /// Returns the alphabet ACGT of the index in ascii.
    #[allow(unused)]
    pub fn alphabet(&self) -> &[u8] {
        &DNA_ALPHABET
    }

    pub fn interval_of_empty_string(&self) -> std::ops::Range<usize> {
        0..self.n_sets()
    }

    /// Returns the C-array of the index.
    /// The array is such that C\[i\] is 1 plus the number of sets in the SBWT that
    /// contain a character that is smaller than the i-th character in the alphabet.
    #[allow(non_snake_case)]
    #[allow(unused)]
    fn C(&self) -> &[usize] {
        self.C.as_slice()
    }

    /// Writes the index to the writer and retuns the number of bytes written.
    /// The array can later be loaded with [SbwtIndex::load].
    pub fn serialize<W: std::io::Write>(&self, out: &mut W) -> std::io::Result<usize> {
        let mut n_written = 0_usize;

        let magic_string_length = [SERIALIZATION_MAGIC_STRING.len() as u8];
        n_written += util::write_bytes(out, &magic_string_length)?;
        n_written += util::write_bytes(out, SERIALIZATION_MAGIC_STRING)?;

        n_written += self.sbwt.serialize(out)?;

        // We're not using serde because we want full control over the bytes
        // in order to guarantee compatibility across languages

        n_written += util::write_bytes(out, &self.n_kmers.to_le_bytes())?;
        n_written += util::write_bytes(out, &self.k.to_le_bytes())?;
        n_written += util::write_bytes(out, &self.C.len().to_le_bytes())?; // TODO: check at build time that usize = u64, crash otherwise
        n_written += util::write_bytes(
            out,
            &self
                .C
                .iter()
                .flat_map(|x| x.to_le_bytes())
                .collect::<Vec<u8>>(),
        )?;

        n_written += self.prefix_lookup_table.serialize(out)?; 

        Ok(n_written)
    }

    /// Loads an index that was previously serialized with [SbwtIndex::serialize].
    #[allow(non_snake_case)] // For C-array
    pub fn load<R: std::io::Read>(input: &mut R) -> std::io::Result<Self> {

        let magic_string_length = input.read_u8().unwrap();
        let mut magic_string_buf = vec![0u8; magic_string_length as usize];
        input.read_exact(&mut magic_string_buf).unwrap();
        if magic_string_buf != SERIALIZATION_MAGIC_STRING {
            panic!("Error loading SBWT: incorrect version string: expected \"{}\", found \"{}\"", String::from_utf8(SERIALIZATION_MAGIC_STRING.to_vec()).unwrap(), String::from_utf8_lossy(&magic_string_buf));
        }

        let subset_rank = SS::load(input)?;

        let n_kmers = byteorder::ReadBytesExt::read_u64::<LittleEndian>(input).unwrap() as usize;
        let k = byteorder::ReadBytesExt::read_u64::<LittleEndian>(input).unwrap() as usize;
        let C_len = byteorder::ReadBytesExt::read_u64::<LittleEndian>(input).unwrap() as usize;

        let mut C_bytes = vec![0_u8; C_len * 8];
        input.read_exact(&mut C_bytes)?;

        let C = C_bytes
            .chunks_exact(8)
            .map(|chunk| {
                let mut arr = [0_u8; 8];
                arr.copy_from_slice(chunk);
                u64::from_le_bytes(arr) as usize
            })
            .collect::<Vec<usize>>();

        let prefix_lookup_table = PrefixLookupTable::load(input)?;

        let index = Self {
            sbwt: subset_rank,
            n_kmers,
            k,
            C,
            prefix_lookup_table
        };

        Ok(index)

    }

    /// Returns the edge label on the incoming edge to the i-th node in colexicographic order 
    /// in the SBWT graph (not the same graph as the de Bruijn graph!), or None if the node has 
    /// no incoming edge. This is well-defined since every node in the SBWT graph has at most 1 incoming edge.
    pub fn inlabel(&self, i: usize) -> Option<u8> {
        if i == 0 {
            return None; // Source node
        }

        for c_idx in 0..DNA_ALPHABET.len() {
            if self.C[c_idx] > i {
                return Some(DNA_ALPHABET[c_idx-1]); // c_idx > 0 because of the source node
            }
        }
        Some(*DNA_ALPHABET.last().unwrap()) // Last character

        // For example, if
        // C = 1, 234, 654, 954
        //
        // Then:
        // [0..1) $
        // [1..234) A
        // [234..654) C
        // [654..954) G
        // [954..n) T
    }

    /// Build select support for the SBWT subset sequence. This is required for
    /// [SbwtIndex::inverse_lf_step] and [SbwtIndex::access_kmer].
    pub fn build_select(&mut self) {
        self.sbwt.build_select();
    }

    /// A low-level function returning `C[char_idx] + SBWT.rank(char_idx, i)`.
    pub fn lf_step(&self, i: usize, char_idx: usize) -> usize {
        self.C[char_idx] + self.sbwt.rank(char_idx as u8, i)
    }

    /// The inverse function of [`SbwtIndex::lf_step`].
    /// Requires that the [select support](SbwtIndex::build_select) has been built.
    /// Returns None if called on the source node of the SBWT graph.
    pub fn inverse_lf_step(&self, i: usize) -> Option<usize> {
        match self.inlabel(i){
            None => None, // Source node
            Some(c) => {
                let char_idx = self.char_idx(c);
                let rank_within_c = i - self.C[char_idx];
                Some(self.sbwt.select(char_idx as u8, rank_within_c).unwrap())
            }
        }
    }

    /// Requires that the [select support](SbwtIndex::build_select) has been built.
    /// Accesses the k-mer associated with the node with colexicographic rank
    /// `colex_rank`. Note that this k-mer may contain dollar symbols if it's
    /// a dummy node. The k-mer string is pushed to `buf`.
    pub fn push_kmer_to_vec(&self, mut colex_rank: usize, buf: &mut Vec<u8>){
        let buf_start = buf.len();
        for _ in 0..self.k {
            match self.inlabel(colex_rank) {
                Some(c) => {
                    buf.push(c);
                    colex_rank = self.inverse_lf_step(colex_rank).unwrap(); // Can unwrap because c != None
                },
                None => {
                    buf.push(b'$');
                },
            }
        }
        buf[buf_start..buf_start+self.k].reverse();
    } 

    /// Requires that the [select support](SbwtIndex::build_select) has been built.
    /// Accesses the k-mer associated with the node with colexicographic rank
    /// `colex_rank`. Note that this k-mer may contain dollar symbols if it's
    /// a dummy node. The k-mer is retuned as a `Vec<u8>`. To push to an existing
    /// buffer for better performance, see [SbwtIndex::push_kmer_to_vec].
    pub fn access_kmer(&self, colex_rank: usize) -> Vec<u8> {
        let mut buf = Vec::with_capacity(self.k);
        self.push_kmer_to_vec(colex_rank, &mut buf);
        buf
    }

    /// Returns the colexicographic interval of the pattern, if found.
    /// The pattern is given in ascii characters.
    pub fn search(&self, pattern: &[u8]) -> Option<std::ops::Range<usize>> {
        let l = self.prefix_lookup_table.prefix_length;
        match pattern.len() >= l {
            true => {
                // Initialize search from lookup table
                let I = self.prefix_lookup_table.lookup(&pattern[0..l]);
                self.search_from(I, &pattern[l..])
            },
            false => {
                // Search from scratch 
                self.search_from(0_usize..self.sbwt.len(), pattern)
            }
        }
    }

    /// Searches the pattern starting from the given colexicographic interval.
    /// Returns the colexicographic interval after searching the pattern, if found.
    /// The pattern is given in ascii characters.
    pub fn search_from(&self, interval: std::ops::Range<usize>, pattern: &[u8]) -> Option<std::ops::Range<usize>> {
        let mut left = interval.start;
        let mut right = interval.end; 
        for chr in pattern.iter() {
            let c = ACGT_TO_0123[*chr as usize];
            if c as usize > DNA_ALPHABET.len() {
                return None; // Character does not exist in index
            }
            left = self.lf_step(left, c as usize);
            right = self.lf_step(right, c as usize);
            if left >= right {
                return None;
            }
        }
        Some(left..right)
    }

    /// Reconstruct all k-mers in the [padded k-spectrum](#sbwt-graph) in the data structure.
    /// The reconstructed k-mers concatenated into a single ascii string in colexicographic order.
    /// The i-th k-mer starts at position i*k in the string.
    pub fn reconstruct_padded_spectrum(&self, n_threads: usize) -> Vec<u8> 
    where Self: Sync {
        let n_nodes = self.n_sets();
        let k = self.k;

        // Build the last column of the SBWT matrix
        log::info!("Building column {} of the SBWT matrix", k-1);
        let mut labels_in = self.build_last_column();
        let mut labels_out = labels_in.clone();

        // Iterate columns from right to left
        let mut kmers_concat = vec![0u8; n_nodes * k];
        for round in 0..k {
            if round > 0 {
                log::info!("Building column {} of the SBWT matrix", k-1-round);
                self.push_all_labels_forward(&labels_in, &mut labels_out, n_threads);
            }

            for i in 0..n_nodes {
                let pos = k - 1 - round;
                kmers_concat[i * k + pos] = labels_out[i];
            }

            (labels_in, labels_out) = (labels_out, labels_in);
        }

        kmers_concat
    }

    /// Set the prefix lookup table of the data structure.
    pub fn set_lookup_table(&mut self, prefix_lookup_table: PrefixLookupTable){
        self.prefix_lookup_table = prefix_lookup_table;
    }

    /// Get the prefix lookup table of the data structure.
    pub fn get_lookup_table(&self) -> &PrefixLookupTable {
        &self.prefix_lookup_table
    }

    /// Get the SBWT set sequence of the data structure.
    pub fn sbwt(&self) -> &SS {
        &self.sbwt
    }

    /// Internal function: construct from parts.
    #[allow(non_snake_case)]
    pub(crate) fn from_components(subset_rank: SS, n_kmers: usize, k: usize, C: Vec<usize>, prefix_lookup_table: PrefixLookupTable) -> Self {
        Self {sbwt: subset_rank, n_kmers, k, C, prefix_lookup_table}
    }

    /// Construct from existing [`crate::subsetseq::SubsetSeq`].
    #[allow(non_snake_case)]
    pub fn from_subset_seq(subset_rank: SS, n_kmers: usize, k: usize, precalc_prefix_length: usize) -> Self {
        let C = subset_rank.get_C_array();
        let n = subset_rank.len();
        let mut index = Self{sbwt: subset_rank, n_kmers, k, C, prefix_lookup_table: PrefixLookupTable::new_empty(n)};
        index.prefix_lookup_table = PrefixLookupTable::new(&index, precalc_prefix_length);

        index
    }

    /// Private helper function for label propagation.
    /// The output range must have length self.n_sets().
    fn split_output_range_by_char<'a>(&self, output_range: &'a mut[u8]) -> Vec::<&'a mut[u8]>{
        let mut remaining_output_range = &mut output_range[1..]; // Drop index 0 because it's the root of the graph
        let mut output_ranges = Vec::<&mut[u8]>::new();

        for i in 0..self.alphabet().len() {
            let range_length = if i + 1 < self.alphabet().len() {
                self.C[i+1] - self.C[i]
            } else {
                self.n_sets() - self.C[i]
            };

            let (head, tail) = remaining_output_range.split_at_mut(range_length); 
            output_ranges.push(head);
            remaining_output_range = tail;
        }
        assert!(remaining_output_range.is_empty());
        output_ranges
    }

    /// Push labels forward along the SBWT graph.
    /// The effect is this: if `labels_in` is a list of labels, one for each node in the SBWT in colex order, then the output
    /// is those labels pushed forward along the edges in the SBWT graph. Those nodes that do not have a predecessor
    /// (just the root of the graph) get a dollar.
    pub fn push_all_labels_forward(&self, labels_in: &[u8], labels_out: &mut [u8], n_threads: usize) 
    where Self: Sync {

        labels_out[0] = b'$'; // The root has no predecessor so it always gets a dollar.
        let mut remaining_char_output_ranges = self.split_output_range_by_char(labels_out);

        let input_piece_len = self.n_sets().div_ceil(n_threads);

        // Subdivide char output ranges for threads.
        // For each thread one input sbwt range and four (= |ACGT|) output ranges. 
        let mut thread_input_output_ranges = Vec::<(std::ops::Range<usize>, Vec::<&mut[u8]>)>::new();

        for t in 0..n_threads {

            let input_range_start = t*input_piece_len;
            let input_range_end = std::cmp::min((t+1)*input_piece_len, self.n_sets());
            let input_subrange = input_range_start..input_range_end;

            let mut thread_output_subranges: Vec<&mut[u8]> = Vec::new();

            let mut range_tails = Vec::<&mut[u8]>::new();
            for (char_idx, output_range) in remaining_char_output_ranges.into_iter().enumerate() {

                // The length of the output range is equal to the number of occurrences of the current character in the input range
                let output_range_len = self.sbwt.rank(char_idx as u8, input_range_end) - self.sbwt.rank(char_idx as u8, input_range_start);
                let (output_subrange, tail) = output_range.split_at_mut(output_range_len); 
                range_tails.push(tail);

                thread_output_subranges.push(output_subrange);
            }
            remaining_char_output_ranges = range_tails;
            
            thread_input_output_ranges.push((input_subrange, thread_output_subranges));
        }

        let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();
        thread_pool.install(||{
            thread_input_output_ranges.into_par_iter().for_each(|(input_subrange, output_subranges)| {
                self.push_labels_forward(&labels_in[input_subrange.clone()], input_subrange.clone(), output_subranges);
            });
        });
    }

    /// Internal function. Takes a range in the SBWT, and a vector of labels (one for each position in the input range).
    /// Output: writes the input labels to the output ranges. Read the code for details. 
    /// The length of output_ranges[j] must be equal to the number of occurrences of character j in the input range.
    pub(crate) fn push_labels_forward(&self, labels: &[u8], input_range: std::ops::Range<usize>, mut output_ranges: Vec<&mut[u8]>){
        assert_eq!(labels.len(), input_range.len());
        assert_eq!(output_ranges.len(), self.alphabet().len());

        let mut range_offsets = vec![0_usize; output_ranges.len()];
        #[allow(clippy::needless_range_loop)]
        for i in input_range.clone() {
            for char_idx in 0..DNA_ALPHABET.len() { // Using the constant DNA_ALPHABET array here to allow loop unrolling
                if self.sbwt.set_contains(i, char_idx as u8) {
                    let input_char = labels[i - input_range.start];
                    output_ranges[char_idx][range_offsets[char_idx]] = input_char;
                    range_offsets[char_idx] += 1;
                }
            }
        }
    }

    /// Build the last column of the SBWT matrix.
    pub fn build_last_column(&self) -> Vec<u8> {
        let mut last = Vec::<u8>::with_capacity(self.n_sets());
        last.push(b'$');
        for c_i in 0..self.alphabet().len(){
            let count = self.sbwt().rank(c_i as u8, self.n_sets());
            let c = self.alphabet()[c_i];
            for _ in 0..count {
                last.push(c);
            }
        }
        assert_eq!(last.len(), self.n_sets());
        last
    } 

    // Internal function: marks the first k-mer of every (k-1)-suffix group.
    pub(crate) fn mark_k_minus_1_mers(&self, n_threads: usize) -> bitvec::vec::BitVec 
    where Self: Sync {
        let n_nodes = self.n_sets();
        let k = self.k();

        // Build the last column of the SBWT matrix
        log::info!("Building column {} of the SBWT matrix", k-1);
        let mut labels_in = self.build_last_column(); 
        let mut labels_out = labels_in.clone();

        let mut marks = bitvec::bitvec![0; self.n_sets()];
        marks.set(0, true);

        // Iterate columns from right to left
        for round in 0..k-1 {
            if round > 0 {
                log::info!("Building column {} of the SBWT matrix", k-1-round);

                self.push_all_labels_forward(&labels_in, &mut labels_out, n_threads);
            }

            for i in 1..n_nodes {
                if labels_out[i] != labels_out[i-1] {
                    marks.set(i, true);
                }
            }

            (labels_in, labels_out) = (labels_out, labels_in);
        }

        marks
    }

}

pub struct MergeInterleaving {
    // Has one bit per colex position in the merged SBWT
    // s1[i] is 1 iff this k-mer is in the first SBWT
    pub s1: BitVec, 

    // Has one bit per colex position in the merged SBWT
    // s2[i] is 1 iff this k-mer is in the second SBWT
    pub s2: BitVec,

    // Has one bit per colex position in the merged SBWT, marking
    // the dummy nodes.
    pub is_dummy: BitVec, 
}

impl MergeInterleaving {

    pub fn new<SS: SubsetSeq + Send + Sync>(index1: &SbwtIndex::<SS>, index2: &SbwtIndex<SS>, n_threads: usize) -> MergeInterleaving {

        let k = index1.k();
        assert_eq!(k, index2.k());

        // We invert the SBWTs column by column and maintain ranges
        // in both SBWTs that are so far equal. The ranges are in increasing
        // order and the partition the SBWTs, to it's enough to just store their
        // sizes in order. The sizes are stored in concatenated unary representations.
        // Empty ranges are allowed.

        // Initialize unary concatenations with empty ranges
        let mut s1 = bitvec![0; index1.n_sets()];
        let mut s2 = bitvec![0; index2.n_sets()];
        s1.push(true);
        s2.push(true);

        let mut chars1 = index1.build_last_column();
        let mut chars2 = index2.build_last_column();

        let mut temp_char_buf_1 = vec![0u8; chars1.len()];
        let mut temp_char_buf_2 = vec![0u8; chars2.len()];

        for round in 0..k {
            log::info!("Round {}/{}", round+1, k);

            let new_arrays = refine_segmentation(s1, s2, &chars1, &chars2);
            (s1, s2) = new_arrays;

            if round != k-1 {
                index1.push_all_labels_forward(&chars1, &mut temp_char_buf_1, n_threads);
                chars1.copy_from_slice(&temp_char_buf_1);

                index2.push_all_labels_forward(&chars2, &mut temp_char_buf_2, n_threads);
                chars2.copy_from_slice(&temp_char_buf_2);
            }
        }

        Self::compress_in_place(&mut s1);
        Self::compress_in_place(&mut s2);

        assert_eq!(s1.len(), s2.len());
        let merged_len = s1.len();

        // Identify dummies in the merged SBWT
        let mut is_dummy = bitvec::vec::BitVec::new();
        is_dummy.resize(merged_len, false);
        let mut c1_idx = 0_usize;
        let mut c2_idx = 0_usize;
        for colex in 0..merged_len {
            let d1 = s1[colex] && c1_idx < chars1.len() && chars1[c1_idx] == b'$';
            let d2 = s2[colex] && c2_idx < chars2.len() && chars2[c2_idx] == b'$';
            is_dummy.set(colex, d1 || d2); 

            c1_idx += s1[colex] as usize;
            c2_idx += s2[colex] as usize;
        }


        MergeInterleaving { s1, s2, is_dummy }
    }

    pub fn intersection_size(&self) -> usize {
        assert_eq!(self.s1.len(), self.s2.len());
        let mut ans = 0_usize;
        for i in 0..self.s1.len() {
            // Do not count dummy nodes
            ans += (!self.is_dummy[i] && self.s1[i] && self.s2[i]) as usize;
        }
        ans
    }

    pub fn union_size(&self) -> usize {
        assert_eq!(self.s1.len(), self.s2.len());
        let mut ans = 0_usize;
        for i in 0..self.s1.len() {
            // Do not count dummy nodes
            ans += (!self.is_dummy[i] && (self.s1[i] || self.s2[i])) as usize;
        }
        ans
    }

    // Helper for construction. Takes a unary concatenation of binary numbers 0^b 1, like:
    // 101011011101 (= 0,1,1,0,1,0,0,1)
    // And produces a bit vector encoding the bits:
    // 01101001
    fn compress_in_place(s1: &mut bitvec::vec::BitVec) {

            let mut s1_i = 0_usize; // Index in s1

            let mut new_idx = 0_usize;
            while s1_i < s1.len() {
                let len1 = leading_zeros(&s1[s1_i..]);
                assert!(len1 <= 1); // This is the colex range of a k-mer, so it should be empty or singleton

                s1_i += len1 + 1; // Length of the unary number we just parsed

                s1.set(new_idx, len1 != 0);
                new_idx += 1;
            }
            assert_eq!(s1_i, s1.len());
            s1.resize(new_idx, false);
            s1.shrink_to_fit();

    }

}

// Assumes there is at least one 1-bit
fn leading_zeros(s: &BitSlice) -> usize {
    s.first_one().unwrap()
}

// Utility function for SBWT merging
fn refine_segmentation(s1: BitVec, s2: BitVec, chars1: &[u8], chars2: &[u8]) -> (BitVec, BitVec) {
    let mut s1_i = 0_usize; // Index in s1
    let mut s2_i = 0_usize; // Index in s2

    let mut c1_i = 0_usize; // Index in chars1
    let mut c2_i = 0_usize; // Index in chars2

    let mut out1 = bitvec::bitvec![];
    let mut out2 = bitvec::bitvec![];

    while s1_i < s1.len() {
        assert!(s2_i < s2.len());

        let len1 = leading_zeros(&s1[s1_i..]);
        let len2 = leading_zeros(&s2[s2_i..]);

        let c1_end = c1_i + len1; // One past the end
        let c2_end = c2_i + len2; // One past the end

        while c1_i < c1_end || c2_i < c2_end {
            let c1 = if c1_i == c1_end { u8::MAX } else { chars1[c1_i] };
            let c2 = if c2_i == c2_end { u8::MAX } else { chars2[c2_i] };
            let c = min(c1,c2);

            while c1_i < c1_end && chars1[c1_i] == c {
                c1_i += 1;
                out1.push(false); // Unary bit
            }

            while c2_i < c2_end && chars2[c2_i] == c {
                c2_i += 1;
                out2.push(false); // Unary bit
            }

            // Terminate unary representations
            out1.push(true);
            out2.push(true);

        }

        assert_eq!(c1_i, c1_end);
        assert_eq!(c2_i, c2_end);

        s1_i += len1 + 1;
        s2_i += len2 + 1;

    }

    assert_eq!(s1_i, s1.len());
    assert_eq!(s2_i, s2.len());

    (out1, out2)
}

/// A table storing the SBWT intervals of all 4^p possible p-mers.
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct PrefixLookupTable {
    /// ranges\[i\] is the interval of the p-mer with colexicographic rank
    /// i in the sorted list of all possible p-mers.
    /// If the p-mer does not exist in the SBWT, the range is [0..0).
    pub ranges: Vec<std::ops::Range<usize>>,

    /// Prefix length p.
    pub prefix_length: usize, 
}

impl PrefixLookupTable {

    /// Create a new prefix lookup table containing only the interval of the
    /// empty string.
    #[allow(clippy::single_range_in_vec_init)] // Clippy false positive. It's actually intended like this
    pub fn new_empty(n_sets_in_sbwt: usize) -> PrefixLookupTable {
        Self{ranges: vec![0..n_sets_in_sbwt], prefix_length: 0}
    }

    /// Create a new prefix lookup table by searching all DNA strings of length `prefix_length`
    /// in the given sbwt.
    pub fn new<SS: SubsetSeq>(sbwt: &SbwtIndex<SS>, prefix_length: usize) -> PrefixLookupTable {
        let mut pmer = vec![0u8; prefix_length];
        let mut ranges = vec![0..0; num::pow(4_usize, prefix_length)];
        for x in 0..num::pow(4, prefix_length) as u64{
            pmer.clear();

            // Construct the p-mer string
            for i in 0..prefix_length {
                let char_idx = (x >> (2*(prefix_length - 1 - i))) & 0x3;
                let c = DNA_ALPHABET[char_idx as usize];
                pmer.push(c);
            }

            if let Some(range) = sbwt.search(&pmer) {
                ranges[x as usize] = range;
            } // Else left as 0..0
        }
        PrefixLookupTable{ranges, prefix_length}
    }

    /// Look up the colex interval of the prefix. 
    pub fn lookup(&self, prefix: &[u8]) -> std::ops::Range<usize> {
        assert!(prefix.len() == self.prefix_length);
        let mut table_idx = 0_usize;
        for (i, c) in prefix.iter().rev().enumerate() {
            let char_idx = ACGT_TO_0123[*c as usize];
            if char_idx == 255 {
                return 0..0; // Not a DNA character
            }
            table_idx |= ((char_idx as u64) << (2*i)) as usize;
        }
        self.ranges[table_idx].clone()
    }

    /// Write the lookup table to the given writer.
    /// The lookup table can be then later loaded with [PrefixLookupTable::load].
    /// Returns number of bytes written. 
    pub fn serialize<W: Write>(&self, out: &mut W) -> std::io::Result<usize> {
        let mut n_written = 0_usize;
        n_written += util::write_bytes(out, &(self.prefix_length as u64).to_le_bytes())?;
        n_written += util::write_bytes(out, &(self.ranges.len() as u64).to_le_bytes())?;
        for range in self.ranges.iter(){
            n_written += util::write_bytes(out, &(range.start as u64).to_le_bytes())?;
            n_written += util::write_bytes(out, &(range.end as u64).to_le_bytes())?;
        }
        Ok(n_written)
    }

    /// Loads a a prefix lookup table that was previosly serialized with [PrefixLookupTable::serialize].
    pub fn load<R: Read>(input: &mut R) -> std::io::Result<Self> {

        let prefix_length = byteorder::ReadBytesExt::read_u64::<LittleEndian>(input).unwrap() as usize;
        let ranges_len = byteorder::ReadBytesExt::read_u64::<LittleEndian>(input).unwrap() as usize;

        let mut ranges = vec![0..0; ranges_len];
        for range in ranges.iter_mut(){
            let start = byteorder::ReadBytesExt::read_u64::<LittleEndian>(input).unwrap() as usize; 
            let end = byteorder::ReadBytesExt::read_u64::<LittleEndian>(input).unwrap() as usize; 
            *range = start..end;
        }

        Ok(Self{ranges, prefix_length})
    }
}

#[cfg(test)]
mod tests {

    use crate::builder::{BitPackedKmerSorting, SbwtIndexBuilder};

    #[cfg(feature = "bpks-mem")]
    use crate::builder::BitPackedKmerSortingMem;

    use super::*;

    #[allow(non_snake_case)]
    fn ACGT_to_0123(c: u8) -> u8 {
        match c {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => panic!("Invalid DNA character {} ", c as char),
        }
    }

    fn encode(s: &str) -> Vec<u8> {
        s.as_bytes()
            .iter()
            .map(|c| ACGT_to_0123(*c))
            .collect::<Vec<u8>>()
    }
    
    #[test_log::test]
    #[allow(non_snake_case)]
    fn doc_example() {
        // The example used in the documentation page for SbwtIndex.
        let seqs: Vec<&[u8]> = vec![b"TGTTTG", b"TTGCTAT", b"ACGTAGTATAT", b"TGTAAA"]; 
        let (sbwt, _) = SbwtIndexBuilder::<BitPackedKmerSorting>::new().k(4).run_from_slices(seqs.as_slice());
        let mut doc_sbwt = vec![vec![b'A',b'T'], vec![b'C'], vec![], vec![b'A'], vec![b'T'], vec![b'T'], vec![b'A',b'G',b'T'], vec![], vec![], vec![b'G'], vec![b'T'], vec![b'T'], vec![b'T'], vec![b'T'], vec![b'C'], vec![b'G'], vec![b'A'], vec![], vec![], vec![b'A'], vec![b'A'], vec![b'A'], vec![b'A',b'T'], vec![b'T'], vec![b'G']];

        // ACGT to 0123
        for set in doc_sbwt.iter_mut() {
            for c in set.iter_mut(){
                *c = sbwt.char_idx(*c) as u8;
            }
        }
        let computed_sbwt: Vec<Vec<u8>> = (0..sbwt.n_sets()).map(|i| sbwt.sbwt.access(i)).collect();
        assert_eq!(doc_sbwt, computed_sbwt);
    }

    #[test_log::test]
    #[allow(non_snake_case)]
    fn LCS_paper_example() {

        let seqs: Vec<&[u8]> = vec![b"AGGTAAA", b"ACAGGTAGGAAAGGAAAGT"];

        let (mut sbwt, _) = SbwtIndexBuilder::<BitPackedKmerSorting>::new().k(4).run_from_slices(seqs.as_slice());

        assert_eq!(sbwt.sbwt.len(), 18);

        assert_eq!(sbwt.sbwt.access(0), encode("A")); //   $$$$
        assert_eq!(sbwt.sbwt.access(1), encode("C")); //   $$$A
        assert_eq!(sbwt.sbwt.access(2), encode("G")); //   GAAA
        assert_eq!(sbwt.sbwt.access(3), encode("")); //    TAAA
        assert_eq!(sbwt.sbwt.access(4), encode("A")); //   GGAA
        assert_eq!(sbwt.sbwt.access(5), encode("A")); //   GTAA
        assert_eq!(sbwt.sbwt.access(6), encode("G")); //   $ACA
        assert_eq!(sbwt.sbwt.access(7), encode("A")); //   AGGA
        assert_eq!(sbwt.sbwt.access(8), encode("AG")); //  GGTA
        assert_eq!(sbwt.sbwt.access(9), encode("A")); //   $$AC
        assert_eq!(sbwt.sbwt.access(10), encode("GT")); // AAAG
        assert_eq!(sbwt.sbwt.access(11), encode("G")); //  ACAG
        assert_eq!(sbwt.sbwt.access(12), encode("G")); //  GTAG
        assert_eq!(sbwt.sbwt.access(13), encode("AT")); // AAGG
        assert_eq!(sbwt.sbwt.access(14), encode("")); //   CAGG
        assert_eq!(sbwt.sbwt.access(15), encode("")); //   TAGG
        assert_eq!(sbwt.sbwt.access(16), encode("")); //   AAGT
        assert_eq!(sbwt.sbwt.access(17), encode("A")); //  AGGT

        let true_padded_spectrum: Vec<&[u8]> = vec![b"$$$$", b"$$$A", b"GAAA", b"TAAA", b"GGAA", b"GTAA", b"$ACA", b"AGGA", b"GGTA", b"$$AC", b"AAAG", b"ACAG", b"GTAG", b"AAGG", b"CAGG", b"TAGG", b"AAGT", b"AGGT"];

        let reconstructed_padded_spectrum = sbwt.reconstruct_padded_spectrum(3);
        sbwt.sbwt.build_select();
        #[allow(clippy::needless_range_loop)]
        for i in 0..18 {
            eprintln!("{} {}",i, String::from_utf8_lossy(&sbwt.access_kmer(i)));
            assert_eq!(sbwt.access_kmer(i), true_padded_spectrum[i]);
            assert_eq!(&reconstructed_padded_spectrum[i*sbwt.k..(i+1)*sbwt.k], true_padded_spectrum[i]);
        }


        assert_eq!(sbwt.search(b""), Some(0..18));
        assert_eq!(sbwt.search(b"AGG"), Some(13..16));
        assert_eq!(sbwt.search(b"AAGT"), Some(16..17));

        // Test prefix looup table
        let two_mers = [b"AA", b"AC", b"AG", b"AT", b"CA", b"CC", b"CG", b"CT", b"GA", b"GC", b"GG", b"GT", b"TA", b"TC", b"TG", b"TT"];
        let lut = PrefixLookupTable::new(&sbwt, 2);
        for two_mer in two_mers {
            let I1 = match sbwt.search(two_mer){
                Some(I) => I,
                None => 0..0
            };
            let I2 = lut.lookup(two_mer);
            assert_eq!(I1, I2);
        }
    }

    #[test]
    fn serialize_and_load() {
        let seqs: Vec<&[u8]> = vec![b"AGGTAAA", b"ACAGGTAGGAAAGGAAAGT"];

        let (sbwt, _) = crate::builder::SbwtIndexBuilder::<BitPackedKmerSorting>::new().k(4).run_from_slices(seqs.as_slice());

        let mut buf = Vec::<u8>::new();
        sbwt.serialize(&mut buf).unwrap();
        let sbwt2 = SbwtIndex::<SubsetMatrix>::load(&mut buf.as_slice()).unwrap();

        assert_eq!(sbwt, sbwt2);
    }

    #[test]
    #[allow(non_snake_case)]
    fn non_ACGT(){
        let seqs: Vec<&[u8]> = vec![b"AGGTAAA", b"ACAGGTAGGANAAGGAAAGT"];           
        //..................................................^...................

        let (sbwt, _) = crate::builder::SbwtIndexBuilder::<BitPackedKmerSorting>::new().k(4).run_from_slices(seqs.as_slice());

        let mut buf = Vec::<u8>::new();
        sbwt.serialize(&mut buf).unwrap();
        let sbwt2 = SbwtIndex::<SubsetMatrix>::load(&mut buf.as_slice()).unwrap();

        assert_eq!(sbwt, sbwt2);
    }

    #[test]
    fn from_subset_seq() {
        let seqs: Vec<&[u8]> = vec![b"AGGTAAA", b"ACAGGTAGGAAAGGAAAGT"];
        let (sbwt_index, _) = crate::builder::SbwtIndexBuilder::<BitPackedKmerSorting>::new().k(4).run_from_slices(seqs.as_slice());
        let ss = sbwt_index.sbwt().clone();
        let sbwt_index2 = SbwtIndex::<SubsetMatrix>::from_subset_seq(ss, sbwt_index.n_sets(), sbwt_index.k(), sbwt_index.prefix_lookup_table.prefix_length);

        let kmers1 = sbwt_index.reconstruct_padded_spectrum(3);
        let kmers2 = sbwt_index2.reconstruct_padded_spectrum(3);
        assert_eq!(kmers1, kmers2);

    }

    // https://stackoverflow.com/questions/52987181/how-can-i-convert-a-hex-string-to-a-u8-slice
    pub fn decode_hex(s: &str) -> Result<Vec<u8>, std::num::ParseIntError> {
        (0..s.len())
            .step_by(3)
            .map(|i| u8::from_str_radix(&s[i..i + 2], 16))
            .collect()
    }


    #[test]
    fn from_cpp_plain_matrix() {

        // C++ SBWT of {ACACTG, GCACTAA}
        // Built with ./build/bin/sbwt build -i small.fna -k 3 -p 2 -o small.sbwt

        let data_hex = concat!(
        "0c 00 00 00 00 00 00 00 70 6c 61 69 6e 2d 6d 61 ",
        "74 72 69 78 04 00 00 00 00 00 00 00 76 30 2e 31 ",
        "0a 00 00 00 00 00 00 00 70 02 00 00 00 00 00 00 ",
        "0a 00 00 00 00 00 00 00 84 00 00 00 00 00 00 00 ",
        "0a 00 00 00 00 00 00 00 01 02 00 00 00 00 00 00 ",
        "0a 00 00 00 00 00 00 00 20 00 00 00 00 00 00 00 ",
        "80 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 ",
        "00 00 00 00 00 00 00 00 80 00 00 00 00 00 00 00 ",
        "00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 ",
        "80 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 ",
        "00 00 00 00 00 00 00 00 80 00 00 00 00 00 00 00 ",
        "00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 ",
        "0a 00 00 00 00 00 00 00 f7 03 00 00 00 00 00 00 ",
        "20 00 00 00 00 00 00 00 01 00 00 00 00 00 00 00 ",
        "05 00 00 00 00 00 00 00 07 00 00 00 00 00 00 00 ",
        "09 00 00 00 00 00 00 00 00 01 00 00 00 00 00 00 ",
        "01 00 00 00 00 00 00 00 01 00 00 00 00 00 00 00 ",
        "02 00 00 00 00 00 00 00 03 00 00 00 00 00 00 00 ",
        "ff ff ff ff ff ff ff ff ff ff ff ff ff ff ff ff ",
        "04 00 00 00 00 00 00 00 04 00 00 00 00 00 00 00 ",
        "05 00 00 00 00 00 00 00 05 00 00 00 00 00 00 00 ",
        "ff ff ff ff ff ff ff ff ff ff ff ff ff ff ff ff ",
        "06 00 00 00 00 00 00 00 06 00 00 00 00 00 00 00 ",
        "ff ff ff ff ff ff ff ff ff ff ff ff ff ff ff ff ",
        "ff ff ff ff ff ff ff ff ff ff ff ff ff ff ff ff ",
        "ff ff ff ff ff ff ff ff ff ff ff ff ff ff ff ff ",
        "ff ff ff ff ff ff ff ff ff ff ff ff ff ff ff ff ",
        "08 00 00 00 00 00 00 00 08 00 00 00 00 00 00 00 ",
        "ff ff ff ff ff ff ff ff ff ff ff ff ff ff ff ff ",
        "09 00 00 00 00 00 00 00 09 00 00 00 00 00 00 00 ",
        "ff ff ff ff ff ff ff ff ff ff ff ff ff ff ff ff ",
        "ff ff ff ff ff ff ff ff ff ff ff ff ff ff ff ff ",
        "02 00 00 00 00 00 00 00 0a 00 00 00 00 00 00 00 ",
        "07 00 00 00 00 00 00 00 03 00 00 00 00 00 00 00");

        let data = decode_hex(data_hex).unwrap();

        let SbwtIndexVariant::SubsetMatrix(cpp_sbwt) = load_from_cpp_plain_matrix_format(&mut std::io::Cursor::new(&data)).unwrap();

        let (rust_sbwt, _lcs) = SbwtIndexBuilder::<BitPackedKmerSorting>::new().k(3).precalc_length(2).run_from_slices(&[b"ACACTG", b"GCACTAA"]);

        assert_eq!(cpp_sbwt, rust_sbwt);

        // Check also that the loading works when the plain-matrix type id is not present.

        let data_without_type_id = &data[20..]; // Skip 8 bytes of (u64 length of type id string) + 12 bytes of "plain-matrix"
        let SbwtIndexVariant::SubsetMatrix(cpp_sbwt2) = load_from_cpp_plain_matrix_format(&mut std::io::Cursor::new(data_without_type_id)).unwrap();

        assert_eq!(cpp_sbwt2, rust_sbwt);

    }

    ///////////////////////////////////////////////////
    //
    // Same tests for bitpacked k-mer sorting in memory
    //
    ///////////////////////////////////////////////////

    #[test_log::test]
    #[allow(non_snake_case)]
    #[cfg(feature = "bpks-mem")]
    fn doc_example_mem() {
        // The example used in the documentation page for SbwtIndex.
        let seqs: Vec<&[u8]> = vec![b"TGTTTG", b"TTGCTAT", b"ACGTAGTATAT", b"TGTAAA"];
        let (sbwt, _) = SbwtIndexBuilder::<BitPackedKmerSortingMem>::new().k(4).run_from_slices(seqs.as_slice());
        let mut doc_sbwt = vec![vec![b'A',b'T'], vec![b'C'], vec![], vec![b'A'], vec![b'T'], vec![b'T'], vec![b'A',b'G',b'T'], vec![], vec![], vec![b'G'], vec![b'T'], vec![b'T'], vec![b'T'], vec![b'T'], vec![b'C'], vec![b'G'], vec![b'A'], vec![], vec![], vec![b'A'], vec![b'A'], vec![b'A'], vec![b'A',b'T'], vec![b'T'], vec![b'G']];

        // ACGT to 0123
        for set in doc_sbwt.iter_mut() {
            for c in set.iter_mut(){
                *c = sbwt.char_idx(*c) as u8;
            }
        }
        let computed_sbwt: Vec<Vec<u8>> = (0..sbwt.n_sets()).map(|i| sbwt.sbwt.access(i)).collect();
        assert_eq!(doc_sbwt, computed_sbwt);
    }

    #[test_log::test]
    #[allow(non_snake_case)]
    #[cfg(feature = "bpks-mem")]
    fn LCS_paper_example_mem() {

        let seqs: Vec<&[u8]> = vec![b"AGGTAAA", b"ACAGGTAGGAAAGGAAAGT"];

        let (mut sbwt, _) = SbwtIndexBuilder::<BitPackedKmerSortingMem>::new().k(4).run_from_slices(seqs.as_slice());

        assert_eq!(sbwt.sbwt.len(), 18);

        assert_eq!(sbwt.sbwt.access(0), encode("A")); //   $$$$
        assert_eq!(sbwt.sbwt.access(1), encode("C")); //   $$$A
        assert_eq!(sbwt.sbwt.access(2), encode("G")); //   GAAA
        assert_eq!(sbwt.sbwt.access(3), encode("")); //    TAAA
        assert_eq!(sbwt.sbwt.access(4), encode("A")); //   GGAA
        assert_eq!(sbwt.sbwt.access(5), encode("A")); //   GTAA
        assert_eq!(sbwt.sbwt.access(6), encode("G")); //   $ACA
        assert_eq!(sbwt.sbwt.access(7), encode("A")); //   AGGA
        assert_eq!(sbwt.sbwt.access(8), encode("AG")); //  GGTA
        assert_eq!(sbwt.sbwt.access(9), encode("A")); //   $$AC
        assert_eq!(sbwt.sbwt.access(10), encode("GT")); // AAAG
        assert_eq!(sbwt.sbwt.access(11), encode("G")); //  ACAG
        assert_eq!(sbwt.sbwt.access(12), encode("G")); //  GTAG
        assert_eq!(sbwt.sbwt.access(13), encode("AT")); // AAGG
        assert_eq!(sbwt.sbwt.access(14), encode("")); //   CAGG
        assert_eq!(sbwt.sbwt.access(15), encode("")); //   TAGG
        assert_eq!(sbwt.sbwt.access(16), encode("")); //   AAGT
        assert_eq!(sbwt.sbwt.access(17), encode("A")); //  AGGT

        let true_padded_spectrum: Vec<&[u8]> = vec![b"$$$$", b"$$$A", b"GAAA", b"TAAA", b"GGAA", b"GTAA", b"$ACA", b"AGGA", b"GGTA", b"$$AC", b"AAAG", b"ACAG", b"GTAG", b"AAGG", b"CAGG", b"TAGG", b"AAGT", b"AGGT"];

        let reconstructed_padded_spectrum = sbwt.reconstruct_padded_spectrum();
        sbwt.sbwt.build_select();
        #[allow(clippy::needless_range_loop)]
        for i in 0..18 {
            eprintln!("{} {}",i, String::from_utf8_lossy(&sbwt.access_kmer(i)));
            assert_eq!(sbwt.access_kmer(i), true_padded_spectrum[i]);
            assert_eq!(&reconstructed_padded_spectrum[i*sbwt.k..(i+1)*sbwt.k], true_padded_spectrum[i]);
        }


        assert_eq!(sbwt.search(b""), Some(0..18));
        assert_eq!(sbwt.search(b"AGG"), Some(13..16));
        assert_eq!(sbwt.search(b"AAGT"), Some(16..17));

        // Test prefix looup table
        let two_mers = [b"AA", b"AC", b"AG", b"AT", b"CA", b"CC", b"CG", b"CT", b"GA", b"GC", b"GG", b"GT", b"TA", b"TC", b"TG", b"TT"];
        let lut = PrefixLookupTable::new(&sbwt, 2);
        for two_mer in two_mers {
            let I1 = match sbwt.search(two_mer){
                Some(I) => I,
                None => 0..0
            };
            let I2 = lut.lookup(two_mer);
            assert_eq!(I1, I2);
        }
    }

    #[test]
    #[cfg(feature = "bpks-mem")]
    fn serialize_and_load_mem() {
        let seqs: Vec<&[u8]> = vec![b"AGGTAAA", b"ACAGGTAGGAAAGGAAAGT"];

        let (sbwt, _) = crate::builder::SbwtIndexBuilder::<BitPackedKmerSortingMem>::new().k(4).run_from_slices(seqs.as_slice());

        let mut buf = Vec::<u8>::new();
        sbwt.serialize(&mut buf).unwrap();
        let sbwt2 = SbwtIndex::<SubsetMatrix>::load(&mut buf.as_slice()).unwrap();

        assert_eq!(sbwt, sbwt2);
    }

    #[test]
    #[allow(non_snake_case)]
    #[cfg(feature = "bpks-mem")]
    fn non_ACGT_mem(){
        let seqs: Vec<&[u8]> = vec![b"AGGTAAA", b"ACAGGTAGGANAAGGAAAGT"];
        //..................................................^...................

        let (sbwt, _) = crate::builder::SbwtIndexBuilder::<BitPackedKmerSortingMem>::new().k(4).run_from_slices(seqs.as_slice());

        let mut buf = Vec::<u8>::new();
        sbwt.serialize(&mut buf).unwrap();
        let sbwt2 = SbwtIndex::<SubsetMatrix>::load(&mut buf.as_slice()).unwrap();

        assert_eq!(sbwt, sbwt2);
    }

    #[test]
    #[cfg(feature = "bpks-mem")]
    fn from_subset_seq_mem() {
        let seqs: Vec<&[u8]> = vec![b"AGGTAAA", b"ACAGGTAGGAAAGGAAAGT"];
        let (sbwt_index, _) = crate::builder::SbwtIndexBuilder::<BitPackedKmerSortingMem>::new().k(4).run_from_slices(seqs.as_slice());
        let ss = sbwt_index.sbwt().clone();
        let sbwt_index2 = SbwtIndex::<SubsetMatrix>::from_subset_seq(ss, sbwt_index.n_sets(), sbwt_index.k(), sbwt_index.prefix_lookup_table.prefix_length);

        let kmers1 = sbwt_index.reconstruct_padded_spectrum();
        let kmers2 = sbwt_index2.reconstruct_padded_spectrum();
        assert_eq!(kmers1, kmers2);

    }
}

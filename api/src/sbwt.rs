//! The [SbwtIndex] data structure. Construct with [SbwtIndexBuilder](crate::SbwtIndexBuilder).

use std::io::Read;
use std::io::Write;
use byteorder::ReadBytesExt;

use byteorder::LittleEndian;
use num::traits::ToBytes;

use crate::subsetseq::*;
use crate::util;
use crate::util::ACGT_TO_0123;
use crate::util::DNA_ALPHABET;

/// The SBWT index data structure. Construct with [SbwtIndexBuilder](crate::SbwtIndexBuilder).
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

pub enum SbwtIndexVariant {
    SubsetMatrix(SbwtIndex<SubsetMatrix>),
    //SubsetConcat(SbwtIndex<SubsetConcat>),
}

pub fn write_sbwt_index_variant(sbwt: &SbwtIndexVariant, out: &mut impl std::io::Write) -> std::io::Result<usize> {
    match sbwt {
        SbwtIndexVariant::SubsetMatrix(sbwt) => {
            out.write(&(b"SubsetMatrix".len() as u64).to_le_bytes())?;
            out.write(b"SubsetMatrix")?;
            sbwt.serialize(out)
        }
    }
}

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

#[derive(Clone, Eq, PartialEq, Debug)]
#[allow(non_snake_case)] // C-array is an established convention in BWT indexes
pub struct SbwtIndex<SS: SubsetSeq> {
    pub(crate) sbwt: SS, // pub(crate) for testing from submodules
    n_kmers: usize,
    k: usize,
    C: Vec<usize>, // Cumulative character counts (includes one ghost dollar)
    prefix_lookup_table: PrefixLookupTable,
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
        n_written += util::write_bytes(out, &SERIALIZATION_MAGIC_STRING)?;

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
    pub fn reconstruct_padded_spectrum(&self) -> Vec<u8> {
        let n_nodes = self.n_sets();
        let k = self.k;

        // Build the last column of the SBWT matrix
        log::info!("Building column {} of the SBWT matrix", k-1);
        let mut last = self.build_last_column();

        // Iterate columns from right to left
        let mut kmers_concat = vec![0u8; n_nodes * k];
        for round in 0..k {
            if round > 0 {
                log::info!("Building column {} of the SBWT matrix", k-1-round);
                last = self.push_labels_forward(&last);
            }

            for i in 0..n_nodes {
                let pos = k - 1 - round;
                kmers_concat[i * k + pos] = last[i];
            }

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
    pub(crate) fn from_components(subset_rank:SS, n_kmers: usize, k: usize, C: Vec<usize>, prefix_lookup_table: PrefixLookupTable) -> Self {
        Self {sbwt: subset_rank, n_kmers, k, C, prefix_lookup_table}
    }

    /// Internal function: Returns vector v such that v[i] = labels[inverse_lf_step(i)].
    pub(crate) fn push_labels_forward(&self, labels: &[u8]) -> Vec<u8> {
        assert_eq!(labels.len(), self.n_sets());
        let mut propagated = vec![b'$'; self.n_sets()];
        let mut ptrs = self.C.clone();
        #[allow(clippy::needless_range_loop)]
        for i in 0..self.n_sets(){
            for char_idx in 0..DNA_ALPHABET.len() {
                if self.sbwt.set_contains(i, char_idx as u8) {
                    propagated[ptrs[char_idx]] = labels[i];
                    ptrs[char_idx] += 1;
                }
            }
        }
        propagated
    }

    /// Internal function: builds the last column of the SBWT matrix.
    pub(crate) fn build_last_column(&self) -> Vec<u8> {
        let mut last = Vec::<u8>::with_capacity(self.n_sets());
        last.push(b'$');
        for c_i in 0..self.alphabet().len(){
            let count = self.sbwt().rank(c_i as u8, self.n_sets());
            for _ in 0..count {
                last.push(self.alphabet()[c_i]);
            }
        }
        assert_eq!(last.len(), self.n_sets());
        last
    } 

    // Internal function: marks the first k-mer of every (k-1)-suffix group.
    pub(crate) fn mark_k_minus_1_mers(&self) -> bitvec::vec::BitVec {
        let n_nodes = self.n_sets();
        let k = self.k();

        // Build the last column of the SBWT matrix
        log::info!("Building column {} of the SBWT matrix", k-1);
        let mut last = self.build_last_column(); 

        let mut marks = bitvec::bitvec![0; self.n_sets()];
        marks.set(0, true);

        // Iterate columns from right to left
        for round in 0..k-1 {
            if round > 0 {
                log::info!("Building column {} of the SBWT matrix", k-1-round);

                last = self.push_labels_forward(&last);
            }

            for i in 1..n_nodes {
                if last[i] != last[i-1] {
                    marks.set(i, true);
                }
            }
        }

        marks
    }

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
}

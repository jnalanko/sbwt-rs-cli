//! An important observation for some operations on the VoDbg is that any correct range of k-mers
//! in the SBWT will contain at most one dummy k-mer which is comprised entirely of $ except for
//! the shared suffix for the given range. Another important note is that this k-mer node will
//! always be located at the beginning of the range.

// Module and submodule contributions by Martin Kostadinov.

#![allow(clippy::ptr_arg)]

use crate::{ContractLeft, ExtendRight, SbwtIndex};
use crate::subsetseq::SubsetSeq;
use pnsv::Pnsv;

use simple_sds_sbwt::bit_vector::BitVector;
use simple_sds_sbwt::int_vector::IntVector;
use simple_sds_sbwt::raw_vector::{AccessRaw, RawVector};
use simple_sds_sbwt::ops::{BitVec, Rank, Access};

pub mod count;
pub mod iter;
/// Module for Previous and Next Smaller Value support.
pub mod pnsv;
pub mod benchmark;
pub mod util {
    /// Given an integer, returns the minimum number of bits which are needed to store it.
    pub fn bit_width(value: usize) -> usize {
        64 - u64::leading_zeros(value as u64) as usize
    }
}

use count::Counts;
use simple_sds_sbwt::serialize::Serialize;

// A struct supporting de Bruijng graph operations on the k-mers stored in the SBWT. More notably
// it supports changing the order i.e. the length of the string corresponding to a given node.
#[derive(Clone, Debug)]
pub struct VoDbg<'a, SS: SubsetSeq + Send + Sync, P: Pnsv + Send + Sync> {
    pub sbwt: &'a SbwtIndex<SS>,
    pub pnsv: &'a P,
    /// A bitvector which marks the dummy nodes in the SBWT.
    dummy_marks: BitVector,
    /// A packed integer array with the maximum length of a suffix of a dummy node which does not
    /// contain $ characters. These lengths appear in the same order as the marks in the
    /// [VoDbg::dummy_marks] bitvector and together with this packed integer array supports random
    /// access to these values given the position of a dummy node.
    dummy_lengths: IntVector,
    counts: Option<Counts>,
}

#[derive(Copy, Clone, Eq, PartialEq, Debug, Hash)]
pub struct Node {
    pub start: usize,
    pub end: usize,
    pub k: usize,
}

#[inline]
pub fn new_node(start: usize, end: usize, k: usize) -> Node {
    Node {
        start,
        end,
        k,
    }
}

pub trait DummyInfo {
    fn is_dummy(&self, position: usize) -> bool;
    fn get_dummy_length(&self, position: usize) -> usize;
}

impl<'a, SS, P> VoDbg<'a, SS, P>
where
    SS: SubsetSeq + Send + Sync,
    P: Pnsv + Send + Sync
{
    /// Marks the dummy k-mers and records the number of $ each has.
    pub fn compute_auxiliary_data_about_dummies(sbwt: &SbwtIndex<SS>) -> (BitVector, IntVector) {
        // Node, depth.
        let mut dfs_stack = Vec::<(usize, usize)>::new(); 
        let mut outlabels = Vec::<u8>::new();

        let mut dummy_count = 0;
        let mut dummy_marks = RawVector::with_len(sbwt.n_sets(), false); // BitVec::new(sbwt.n_sets());

        // First pass to calculate the dummies.
        dfs_stack.push((0, 0)); // Colex rank of $, depth of $
        while let Some((node, depth)) = dfs_stack.pop() { 
            if !dummy_marks.bit(node) {
                dummy_count += 1;
            }

            dummy_marks.set_bit(node, true);

            if depth + 1 < sbwt.k() {
                outlabels.clear();
                sbwt.sbwt.append_set_to_buf(node, &mut outlabels);
                for &c_idx in outlabels.iter() {
                    let u = sbwt.lf_step(node, c_idx as usize);
                    dfs_stack.push((u, depth + 1));
                }
            }
        }
        
        let mut dummy_marks = BitVector::from(dummy_marks);// Rank9::new(dummy_marks);
        dummy_marks.enable_rank();

        let bit_width = util::bit_width(sbwt.k());
        let mut dummy_lengths = IntVector::with_len(dummy_count, bit_width, 0).unwrap();

        // Second pass to calculate their depths now that we have their positions.
        dfs_stack.push((0, 0)); // Colex rank of $, depth of $
        while let Some((node, depth)) = dfs_stack.pop() { 
            let dummy_index = dummy_marks.rank(node);
            dummy_lengths.set(dummy_index, depth as u64);

            if depth + 1 < sbwt.k() {
                outlabels.clear();
                sbwt.sbwt.append_set_to_buf(node, &mut outlabels);
                for &c_idx in outlabels.iter() {
                    let u = sbwt.lf_step(node, c_idx as usize);
                    dfs_stack.push((u, depth + 1));
                }
            }
        }

        (dummy_marks, dummy_lengths)
    }

    /// Initializes supports for de Bruijn graph operation based on the given [SbwtIndex].
    /// If the Lcs array of the SBWT is available, it can be given to significantly speed up construction.
    /// IMPORTANT: [select support][SbwtIndex::build_select()] must be built before calling this function. 
    pub fn new(sbwt: &'a SbwtIndex<SS>, pnsv: &'a P) -> Self
    {
        assert!(sbwt.sbwt.has_select_support());
        log::info!("[VoDbg::new] computing auxiliary data about dummy k-mers...");
        let (dummy_marks, dummy_lengths) = Self::compute_auxiliary_data_about_dummies(sbwt);
        Self {
            sbwt,
            pnsv,
            dummy_marks,
            dummy_lengths,
            counts: None,
        }
    }

    pub fn counts(&self) -> Option<&Counts> {
        self.counts.as_ref()
    }

    #[allow(clippy::result_unit_err)]
    pub fn build_counts<Stream>(
        &mut self,
        sequence_stream: Stream,
        use_hash_map: bool,
        sample_distance: usize,
        additional_memory_bound_gb: usize,
        thread_count: usize,
        batch_size: usize,
    ) -> Result<(), ()>
    where
        Stream: crate::SeqStream + Send + Clone,
    {
        if self.counts.is_some() {
            return Ok(());
        }
        let streaming_index = crate::StreamingIndex {
            extend_right: self.sbwt,
            contract_left: self.pnsv,
            n: self.sbwt.n_sets(),
            k: self.sbwt.k(),
        };
        let result = if use_hash_map {
            Counts::try_new_concurrent_with_hashmap(
                sequence_stream,
                &streaming_index,
                self,
                sample_distance,
                additional_memory_bound_gb,
                thread_count,
                batch_size
            )
        } else {
            Counts::try_new_concurrent_two_passes(
                sequence_stream,
                &streaming_index,
                self,
                sample_distance,
                additional_memory_bound_gb,
                thread_count,
                batch_size
            )
        };
        if let Some(counts) = result {
            self.counts = Some(counts);
            Ok(())
        } else {
            Err(())
        }
    }

    /// Push the k-mer string of the node to the given buffer.
    pub fn push_node_kmer(&self, node: Node, buf: &mut Vec<u8>) {
        // assert!(node.k < self.sbwt.k() || !self.is_dummy(node.start));
        let mut colex_rank = node.start;
        let buf_start = buf.len();
        for _ in 0..node.k {
            match self.sbwt.inlabel(colex_rank) {
                Some(c) => {
                    buf.push(c);
                    // Can unwrap because c != None.
                    colex_rank = self.sbwt.inverse_lf_step(colex_rank).unwrap(); 
                },
                None => {
                    buf.push(b'$');
                },
            }
        }
        buf[buf_start..buf_start+node.k].reverse();
    }

    /// Get the k-mer string label of a node. To avoid memory allocation, check
    /// [VoDbg::push_node_kmer].
    pub fn get_kmer(&self, node: Node) -> Vec<u8> {
        // assert!(node.k < self.sbwt.k() || !self.is_dummy(node.start));
        let mut buf = Vec::<u8>::with_capacity(node.k);
        self.push_node_kmer(node, &mut buf);
        buf
    }

    /// Get the number of occurrences of the k-mer corresponding to the given k-mer. Support for
    /// the [counts][VoDbg::build_counts] must have been called beforehand, otherwise this method
    /// returns 0.
    pub fn get_count(&self, node: Node) -> u64 {
        if let Some(counts) = self.counts.as_ref() {
            counts.range_sum(node.start, node.end)
        } else {
            0
        }
    }

    /// Get a handle to the node corresponding to the given k-mer, if exists in the graph.
    pub fn get_node(&self, kmer: &[u8]) -> Option<Node> {
        assert!(kmer.len() <= self.sbwt.k());
        self.sbwt.search(kmer).map(|range| new_node(
            range.start,
            range.end,
            kmer.len(),
        ))
    }

    /// Climb to a "lower" level of the graph where the strings in the nodes are shorter by
    /// removing characters from the left.
    pub fn contract_left(&self, node: Node, target_length: usize) -> Node {
        assert!(node.k < self.sbwt.k() || !self.is_dummy(node.start));
        assert!(node.k > target_length);
        let range = self.pnsv.contract_left(node.start..node.end, target_length);
        new_node(range.start, range.end, target_length)
    }

    /// Climb to a "lower" level of the graph where the strings in the nodes are shorter by
    /// removing characters from the right.
    pub fn contract_right(&self, node: Node, target_length: usize) -> Option<Node> {
        assert!(node.k < self.sbwt.k() || !self.is_dummy(node.start));
        assert!(node.k > target_length);
        let representative = self.get_representative(node);
        let mut current_length = node.k;
        let mut start = representative;
        while current_length > target_length {
            start = self.sbwt.inverse_lf_step(start)?;
            current_length -= 1;
        }
        let end = self.pnsv.next(start + 1, target_length);
        let start = self.pnsv.previous(start, target_length);
        let node = new_node(start, end, target_length);
        Some(node)
    }

    /// Climb to an "upper" level of the graph where the strings in the nodes are one longer. Given
    /// an index returns the corresponding node in colexicographic order if it exists.
    pub fn extend_left_with_index(&self, node: Node, index: usize) -> Option<Node> {
        assert!(node.k < self.sbwt.k() || !self.is_dummy(node.start));

        let mut start = node.start;
        // Skip the dummy node which has only $ except for the suffix of the range.
        if self.is_dummy(start) && self.get_dummy_length(start) <= node.k {
            start += 1;
        }

        let target_length = node.k + 1;
        let mut end;
        let mut current_index = 0;
        while start < node.end {
            end = self.pnsv.next(start + 1, target_length);
            if current_index == index {
                let extended_node = new_node(start, end, target_length);
                return Some(extended_node);
            }
            start = end;
            current_index += 1;
        }

        None
    }

    /// Climb to an "upper" level of the graph where the strings in the nodes are one longer.
    /// Returns a node only if it exists.
    pub fn extend_left_with_character(&self, node: Node, character: u8, kmer_buffer: &mut Vec<u8>) -> Option<Node> {
        assert!(node.k < self.sbwt.k() || !self.is_dummy(node.start));
        // let mut kmer_buffer = Vec::<u8>::with_capacity(node.k + 1);

        let in_left_half = {
            let half = self.sbwt.alphabet().len() >> 1;
            match self.sbwt.alphabet().iter().position(|&c| c == character) {
                Some(index) => index < half,
                None => return None,
            }
        };

        let mut start = node.start;
        // Skip the dummy node which has only $ except for the suffix of the range.
        if self.is_dummy(start) && self.get_dummy_length(start) <= node.k {
            start += 1;
        }

        let target_length = node.k + 1;
        let mut end; 

        if in_left_half {
            while start < node.end {
                end = self.pnsv.next(start + 1, target_length);
                let extended_node = new_node(start, end, target_length);
                kmer_buffer.clear();
                self.push_node_kmer(extended_node, kmer_buffer);
                
                // It is guaranteed that there is at least one character in the suffix.
                if kmer_buffer[0] == character {
                    return Some(extended_node);
                }
                start = end;
            }

            return None;
        }

        end = node.end;
        let mut new_start;
        while end > start {
            new_start = self.pnsv.previous(end - 1, target_length);
            let extended_node = new_node(new_start, end, target_length);
            kmer_buffer.clear();
            self.push_node_kmer(extended_node, kmer_buffer);

            if kmer_buffer[0] == character {
                return Some(extended_node);
            }

            end = new_start;
        }

        None
    }

    /// Climb to an "upper" level of the graph where the strings in the nodes are one longer.
    pub fn extend_right(&self, node: Node, character: u8) -> Option<Node> {
        let result = self.sbwt.extend_right(node.start..node.end, character);
        let length_increase = if node.k < self.sbwt.k() { 1 } else { 0 };
        if result.is_empty() {
            return None;
        }
        let node = new_node(result.start, result.end, node.k + length_increase);
        Some(node)
    }

    /// Returns the number of outgoing edges from the given node.
    pub fn outdegree(&self, node: Node) -> usize {
        assert!(node.k < self.sbwt.k() || !self.is_dummy(node.start));
        // note(mk): We have to do a contract left and an extend right to find whether the given
        // neighbour exists (node centric de Bruijn graph). This is a costly operation. This note
        // is also valid for every other location where information about the out neighbours is
        // needed.
        let mut outdegree = 0;
        let contracted = self.contract_left(node, node.k - 1);
        for &c in self.sbwt.alphabet() {
            if self.extend_right(contracted, c).is_some() {
                outdegree += 1;
            }
        }
        outdegree
    }

    /// Returns the number of incoming edges to the given node.
    pub fn indegree(&self, node: Node) -> usize {
        assert!(node.k < self.sbwt.k() || !self.is_dummy(node.start));

        //
        // The overall idea is to count the number of ranges for suffixes of length node.k in the
        // range of the suffix which is equal to the (node.k-1) prefix of the node's k-mer.
        //

        let in_neighbours_whole_range = self.contract_right(node, node.k - 1);
        let (mut in_neighbour_start, end) = if let Some(node) = in_neighbours_whole_range {
            (node.start, node.end)
        } else {
            return 0;
        };

        //
        // Let's assume the given node's k-mer is ---ACG i.e. the max k is 6, but the node
        // is part of a lower order de Bruijn graph. After the contract right operation the maximal
        // k-mers in the new range will be of the form ----AC. We want to find the number of slices
        // of the range which have the same (node.k)-suffix. For example, if the range contains:
        //
        // ---AAC
        // ---AAC
        // ---CAC
        // ---GAC
        //
        // then the indegree of the node should be 3. As per the observation about the ranges,
        // there might be a single k-mer at the beginning of the range whose (node.k)-suffix
        // contains at most 1 $ symbol. We should skip it.
        //

        if self.is_dummy(in_neighbour_start) {
            let length = self.get_dummy_length(in_neighbour_start);
            if length < node.k {
                in_neighbour_start += 1;
            }
        }

        let mut count = 0;
        let target_length = node.k;
        while in_neighbour_start < end {
            count += 1;
            in_neighbour_start = self.pnsv.next(in_neighbour_start + 1, target_length);
        }

        count
    }
 
    /// The k-mers in the SBWT which are colexicographically the smallest with a given (k-1) suffix
    /// will be referred to as "representatives" i.e. those k-mers in the SBWT which (should) have
    /// a non-empty set of outgoing edges (where the value k in (k-1) is equal to the k of the
    /// SBWT). This method returns the reprsentative for the k-mer in the SBWT at the start of the
    /// range of the given node.
    fn get_representative(&self, node: Node) -> usize {
        // For reference I will refer to nodes and k-mers whose length is equal to the k-mers
        // in the SBWT as "maximal" and shorter k-mers and their corresponding nodes as "non-maximal".
        //
        // If the node is non-maximal, then the start of its range is guaranteed to be a
        // representative k-mer. This is because if the start of the range is not the first element
        // with a given (k-1) suffix, then there exists a k-mer immediately before it (in the
        // colexicographic order) with the same (k-1) suffix. This means that this previous k-mer
        // shares any suffix with length less than or equal to (k-1) with the k-mer at the start of
        // the range which means that the range is not correct. Assuming that every node has a
        // correct range, this is a contradiction.
        //
        // If, however, the node is maximal, then its range contains only 1 k-mer which is not
        // guaranteed to be a representative and thus we must find that representative.
        if node.k == self.sbwt.k() {
            self.get_representative_of_maximal_node(node)
        } else {
            node.start
        }
    }

    /// Returns the colex rank of the smallest k-mer (possibly dummy) that has the same suffix of
    /// length (k-1) as the given colex position (possibly dummy). Serves a similar purpose as
    /// [super::dbg::Dbg::get_suffix_group_start()] from the Dbg structure.
    #[inline]
    fn get_representative_of_maximal_node(&self, node: Node) -> usize {
        self.pnsv.previous(node.start, self.sbwt.k() - 1)
    }

    /// For each outgoing edge from the given node to nodes in the same order de Bruijn graph,
    /// pushes to the output vector a pair (v, c), where v is the target node and c is the edge
    /// label.
    pub fn push_out_neighbors(&self, node: Node, output: &mut Vec<(Node, u8)>) {
        assert!(node.k < self.sbwt.k() || !self.is_dummy(node.start));
        let contracted = self.contract_left(node, node.k - 1);
        for &c in self.sbwt.alphabet() {
            if let Some(neighbor) = self.extend_right(contracted, c) {
                output.push((neighbor, c));
            }
        }
    }

    /// For each incoming edge to the given node in the same order de Bruijn graph, pushes to the
    /// output vector a pair (v, c), where v is the source node and c is the edge label. The edge
    /// label will be the same for all in-neighbors because it has to be equal to the last character
    /// of the destination k-mer.
    pub fn push_in_neighbors(&self, node: Node, output: &mut Vec<(Node, u8)>) {
        assert!(node.k < self.sbwt.k() || !self.is_dummy(node.start));
        let inlabel = self.get_last_character(node);
        let in_neighbours_whole_range = self.contract_right(node, node.k - 1);
        let (mut in_neighbour_start, end) = if let Some(node) = in_neighbours_whole_range {
            (node.start, node.end)
        } else {
            return;
        };

        if self.is_dummy(in_neighbour_start) {
            let length = self.get_dummy_length(in_neighbour_start);
            if length < node.k {
                in_neighbour_start += 1;
            }
        }

        let target_length = node.k;
        while in_neighbour_start < end {
            let in_neighbour_end = self.pnsv.next(in_neighbour_start + 1, target_length);
            let innode = new_node(in_neighbour_start, in_neighbour_end, target_length);
            output.push((innode, inlabel));
            in_neighbour_start = in_neighbour_end;
        }
    }

    /// Gets the last character of the k-mer string of the given node. Panics if the k-mer the node
    /// represents is empty.
    pub fn get_last_character(&self, node: Node) -> u8 {
        assert!(0 < node.k);
        assert!(node.k < self.sbwt.k() || !self.is_dummy(node.start));
        // If the length of the suffix of the node is greater than 0, then the last character of
        // the k-mer this node represents is not a $. That is, node.start == 0 if and only if
        // node.k == 0.
        self.sbwt.inlabel(node.start).unwrap() 
    }

    /// Returns whether the given node has an outgoing edge labeled with `edge_label` in the same
    /// order de Bruijn graph.
    pub fn has_outlabel(&self, node: Node, edge_label: u8) -> bool {
        assert!(node.k < self.sbwt.k() || !self.is_dummy(node.start));
        self.follow_outedge(node, edge_label).is_some()
    }

    /// Pushes the labels of all outgoing edges in the same order de Bruijn graph from the given
    /// node to the output vector.
    pub fn push_outlabels(&self, node: Node, output: &mut Vec<u8>) {
        assert!(node.k < self.sbwt.k() || !self.is_dummy(node.start));
        let contracted = self.contract_left(node, node.k - 1);
        for &c in self.sbwt.alphabet() {
            if self.extend_right(contracted, c).is_some() {
                output.push(c);
            }
        }
    }

    /// Follows the outgoing edge labeled with `edge_label` from the given node in the same order de
    /// Bruijn graph. Returns None if the edge does not exist.
    pub fn follow_outedge(&self, node: Node, edge_label: u8) -> Option<Node>{
        assert!(node.k < self.sbwt.k() || !self.is_dummy(node.start));
        // The SBWT is based on a node-centric de Bruijn graph and thus it is not guaranteed that
        // for each two nodes there is a corresponding (k+1)-mer.
        let contracted = self.contract_left(node, node.k - 1);
        self.extend_right(contracted, edge_label)
    }

    /// Follows backward the incoming edge that comes from the i-th smallest k-mer
    /// (i ∈ [0, indegree(node)) in colexicographic order that has an outgoing edge to `node` in the
    /// same order de Bruijn graph. Returns None if i ≥ indegree(node).
    pub fn follow_inedge(&self, node: Node, index: usize) -> Option<Node>{
        assert!(node.k < self.sbwt.k() || !self.is_dummy(node.start));
        let in_neighbours_whole_range = self.contract_right(node, node.k - 1);
        let (mut in_neighbour_start, end) = if let Some(node) = in_neighbours_whole_range {
            (node.start, node.end)
        } else {
            return None;
        };

        if self.is_dummy(in_neighbour_start) {
            let length = self.get_dummy_length(in_neighbour_start);
            if length < node.k {
                in_neighbour_start += 1;
            }
        }

        let target_length = node.k;
        let mut current_index = 0;
        while in_neighbour_start < end {
            let in_neighbour_end = self.pnsv.next(in_neighbour_start + 1, target_length);
            if current_index == index {
                let innode = new_node(in_neighbour_start, in_neighbour_end, target_length);
                return Some(innode);
            }
            if in_neighbour_end >= end {
                break;
            }
            current_index += 1;
            in_neighbour_start = in_neighbour_end;
        }
        None
    }

    const SOME_COUNT: u8 = 1;
    const NONE_COUNT: u8 = 1 - Self::SOME_COUNT;

    pub fn serialize<W: std::io::Write>(&self, out: &mut W) -> std::io::Result<usize> {
        use byteorder::WriteBytesExt;
        let mut written = 0;
        if let Some(count) = self.counts.as_ref() {
            out.write_u8(Self::SOME_COUNT)?;
            written += count.serialize(out)?;
        } else {
            out.write_u8(Self::NONE_COUNT)?;
        }
        written += 1;
        self.dummy_marks.serialize(out)?;
        self.dummy_lengths.serialize(out)?;
        written += self.dummy_marks.size_in_bytes();
        written += self.dummy_lengths.size_in_bytes();
        Ok(written)
    }

    pub fn load<R: std::io::Read>(input: &mut R, sbwt: &'a SbwtIndex<SS>, pnsv: &'a P) -> std::io::Result<Self> {
        use byteorder::ReadBytesExt;
        let count_is_some = input.read_u8()?;
        let counts = if count_is_some == Self::SOME_COUNT {
            Some(Counts::load(input)?)
        } else {
            None
        };
        let dummy_marks = BitVector::load(input)?;
        let dummy_lengths = IntVector::load(input)?;
        let result = Self {
            sbwt,
            pnsv,
            dummy_marks,
            dummy_lengths,
            counts,
        };
        Ok(result)
    }
}

impl<'a, SS, P> DummyInfo for VoDbg<'a, SS, P>
where
    SS: SubsetSeq + Send + Sync,
    P: Pnsv + Send + Sync
{
    #[inline]
    fn is_dummy(&self, position: usize) -> bool {
        self.dummy_marks.get(position)
    }

    fn get_dummy_length(&self, position: usize) -> usize {
        assert!(self.is_dummy(position));
        let dummy_index = self.dummy_marks.rank(position);
        self.dummy_lengths.get(dummy_index) as usize
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{BitPackedKmerSortingMem, LcsArray, SbwtIndexBuilder, SubsetMatrix};
    use crate::dbg::Dbg;
    use pnsv::PnsvTuned;

    #[test]
    fn serialize_and_load() {
        use rand_chacha::ChaCha20Rng;
        use rand_chacha::rand_core::SeedableRng;
        use rand_chacha::rand_core::RngCore;

        let max_k: usize = 16;
        let kmer_count = 128;
        let mut rng = ChaCha20Rng::from_seed([42; 32]);

        let mut seqs = Vec::<Vec<u8>>::new();
        for _ in 0..kmer_count {
            let kmer: Vec<u8> = (0..max_k).map(|_| match rng.next_u32() % 4 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            }).collect();
            seqs.push(kmer);
        }

        seqs.sort();
        seqs.dedup();

        let (sbwt, lcs) = SbwtIndexBuilder::<BitPackedKmerSortingMem>::new()
            .k(max_k).build_lcs(true)
            .add_all_dummy_paths(true)
            .build_select_support(true)
            .run_from_vecs(&seqs);
        let lcs = lcs.unwrap();
        let pnsv_tuned = PnsvTuned::new_default(&sbwt, &lcs, max_k);

        let mut vodbg = VoDbg::new(&sbwt, &pnsv_tuned);

        // Without Counts being built.
        let mut buffer = Vec::<u8>::new();
        let written = vodbg.serialize(&mut buffer).unwrap();
        assert_eq!(written, buffer.len());
        let vodbg_loaded = VoDbg::load(&mut buffer.as_slice(), &sbwt, &pnsv_tuned).unwrap();
        assert_eq!(vodbg.counts, vodbg_loaded.counts);
        assert_eq!(vodbg.dummy_marks, vodbg_loaded.dummy_marks);
        assert_eq!(vodbg.dummy_lengths, vodbg_loaded.dummy_lengths);

        // With Counts being built.
        let sequence_stream = crate::util::VecSeqStream::new(&seqs);
        vodbg.build_counts(
            sequence_stream,
            true,
            Counts::DEFAULT_SAMPLE_DISTANCE,
            1, 4,
            Counts::DEFAULT_BATCH_SIZE_IN_BYTES
        ).unwrap();
        buffer.clear();
        let written = vodbg.serialize(&mut buffer).unwrap();
        assert_eq!(written, buffer.len());
        let vodbg_loaded = VoDbg::load(&mut buffer.as_slice(), &sbwt, &pnsv_tuned).unwrap();
        assert_eq!(vodbg.counts, vodbg_loaded.counts);
        assert_eq!(vodbg.dummy_marks, vodbg_loaded.dummy_marks);
        assert_eq!(vodbg.dummy_lengths, vodbg_loaded.dummy_lengths);
    }

    #[test]
    fn randomised_kmers() {
        use rand_chacha::ChaCha20Rng;
        use rand_chacha::rand_core::SeedableRng;
        use rand_chacha::rand_core::RngCore;

        const MIN_K: usize = 3;
        let max_k: usize = 20;
        let kmer_count = 256;
        let mut rng = ChaCha20Rng::from_seed([42; 32]);

        let mut seqs = Vec::<Vec<u8>>::new();
        for _ in 0..kmer_count {
            let kmer: Vec<u8> = (0..max_k).map(|_| match rng.next_u32() % 4 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            }).collect();
            seqs.push(kmer);
        }

        seqs.sort();
        seqs.dedup();

        let mut sbwt_indices: Vec<(SbwtIndex<SubsetMatrix>, Option<LcsArray>)> = Vec::with_capacity(max_k);
        let mut graphs = Vec::with_capacity(max_k);

        for i in MIN_K..=max_k {
            let (sbwt, lcs) = SbwtIndexBuilder::<BitPackedKmerSortingMem>::new()
                .k(i).build_lcs(true)
                .build_select_support(true)
                .run_from_vecs(seqs.as_slice());
            sbwt_indices.push((sbwt, lcs));
        }

        for i in 0..sbwt_indices.len() {
            let dbg = Dbg::new(&sbwt_indices[i].0, sbwt_indices[i].1.as_ref(), 3);
            graphs.push(dbg);
        }

        let vodbg_sbwt = &sbwt_indices[max_k - MIN_K].0;
        let vodbg_lcs = sbwt_indices[max_k - MIN_K].1.as_ref().unwrap();
        let pnsv_tuned = PnsvTuned::new_default(vodbg_sbwt, vodbg_lcs, max_k);
        let vodbg = VoDbg::new(vodbg_sbwt, &pnsv_tuned);

        let alphabet = sbwt_indices[max_k - MIN_K].0.alphabet();

        let mut dbg_buffer = vec![];
        let mut vodbg_buffer = vec![];

        let mut dbg_outlabels = vec![];
        let mut vodbg_outlabels = vec![];

        let mut kmer_buffer = vec![];

        for current_k in MIN_K..=max_k {
            let dbg_index = current_k - MIN_K;
            let dbg = &graphs[dbg_index];
            for sequence_index in 0..seqs.len() {
                for sequence_start in 0..max_k-current_k {
                    let sequence = &seqs[sequence_index][sequence_start..sequence_start + current_k];

                    // get_node
                    let dbg_node = dbg.get_node(sequence).expect("Should exist.");
                    let vodbg_node = vodbg.get_node(sequence).expect("Should exist.");
 
                    // get_kmer
                    assert_eq!(dbg.get_kmer(dbg_node), sequence);
                    assert_eq!(vodbg.get_kmer(vodbg_node), sequence);

                    // indegree, outdegree
                    assert_eq!(dbg.indegree(dbg_node), vodbg.indegree(vodbg_node));
                    assert_eq!(dbg.outdegree(dbg_node), vodbg.outdegree(vodbg_node));

                    // push_in_neighbors
                    dbg_buffer.clear();
                    dbg.push_in_neighbors(dbg_node, &mut dbg_buffer);
                    vodbg_buffer.clear();
                    vodbg.push_in_neighbors(vodbg_node, &mut vodbg_buffer);

                    assert_eq!(dbg_buffer.len(), vodbg_buffer.len());
                    for i in 0..dbg_buffer.len() {
                        assert_eq!(dbg_buffer[i].1, vodbg_buffer[i].1);
                        let dbg_in_neighbor = dbg_buffer[i].0;
                        let vodbg_in_neighbor = vodbg_buffer[i].0;

                        let dbg_in_neighbour_kmer = dbg.get_kmer(dbg_in_neighbor);
                        let vodbg_in_neighbour_kmer = vodbg.get_kmer(vodbg_in_neighbor);
                        assert_eq!(dbg_in_neighbour_kmer, vodbg_in_neighbour_kmer);

                        // follow_inedge
                        let dbg_from_inedge = dbg.follow_inedge(dbg_node, i).expect("Should exist.");
                        let vodbg_from_inedge = vodbg.follow_inedge(vodbg_node, i).expect("Should exist.");
                        assert_eq!(dbg_in_neighbour_kmer, dbg.get_kmer(dbg_from_inedge));
                        assert_eq!(vodbg_in_neighbour_kmer, vodbg.get_kmer(vodbg_from_inedge));
                    }

                    // extend_left_with_character
                    let vodbg_contracted = vodbg.contract_left(vodbg_node, vodbg_node.k - 1);
                    let vodbg_extended = vodbg.extend_left_with_character(vodbg_contracted, sequence[0], &mut kmer_buffer)
                        .expect("Should exist.");
                    assert_eq!(vodbg_node, vodbg_extended);

                    // extend_left_with_index
                    {
                        let mut index = 0;
                        loop {
                            let extended_op = vodbg.extend_left_with_index(vodbg_node, index);
                            match extended_op {
                                Some(extended) => {
                                    let contracted = vodbg.contract_left(extended, vodbg_node.k);
                                    assert_eq!(contracted, vodbg_node);
                                },
                                None => {
                                    break;
                                }
                            }
                            index += 1;
                        }
                    }

                    // push_out_neighbors
                    dbg_buffer.clear();
                    dbg.push_out_neighbors(dbg_node, &mut dbg_buffer);
                    vodbg_buffer.clear();
                    vodbg.push_out_neighbors(vodbg_node, &mut vodbg_buffer);

                    assert_eq!(dbg_buffer.len(), vodbg_buffer.len());
                    for i in 0..dbg_buffer.len() {
                        assert_eq!(dbg_buffer[i].1, vodbg_buffer[i].1);
                        let dbg_out_neighbor = dbg_buffer[i].0;
                        let vodbg_out_neighbor = vodbg_buffer[i].0;
                        assert_eq!(dbg.get_kmer(dbg_out_neighbor), vodbg.get_kmer(vodbg_out_neighbor));
                    }

                    // get_last_character
                    assert_eq!(dbg.get_last_character(dbg_node), vodbg.get_last_character(vodbg_node));

                    // has_outlabel
                    for &c in alphabet {
                        assert_eq!(dbg.has_outlabel(dbg_node, c), vodbg.has_outlabel(vodbg_node, c));
                    }

                    // push_outlabels
                    dbg_outlabels.clear();
                    dbg.push_outlabels(dbg_node, &mut dbg_outlabels);
                    vodbg_outlabels.clear();
                    vodbg.push_outlabels(vodbg_node, &mut vodbg_outlabels);
                    assert_eq!(dbg_outlabels, vodbg_outlabels);
                }
            }
        }

        const EXTRA_KMERS_LOWER_BOUND_K: usize = 12;
        const EXTRA_SEQUENCE_COUNT: usize = 512;
        let mut sequence: Vec<u8> = Vec::with_capacity(max_k);
        for current_k in EXTRA_KMERS_LOWER_BOUND_K..=max_k {
            let dbg_index = current_k - MIN_K;
            let dbg = &graphs[dbg_index];
            for _ in 0..EXTRA_SEQUENCE_COUNT {
                let iterator = (0..current_k).map(|_| match rng.next_u32() % 4 {
                    0 => b'A',
                    1 => b'C',
                    2 => b'G',
                    _ => b'T',
                });
                sequence.clear();
                sequence.extend(iterator);

                let dbg_node = dbg.get_node(&sequence);
                let vodbg_node = vodbg.get_node(&sequence);

                assert_eq!(dbg_node.is_some(), vodbg_node.is_some());
            }
        }
    }

    #[test]
    fn smaller_values_of_k() {
        // All possible 2-mers.
        let seqs = vec![b"AACCGGTTAGCTGATCA".to_vec()];
        let (sbwt, lcs) = SbwtIndexBuilder::<BitPackedKmerSortingMem>::new()
            .k(3).build_lcs(true)
            .build_select_support(true)
            .run_from_vecs(seqs.as_slice());
        let lcs = lcs.unwrap();

        let pnsv_tuned = PnsvTuned::new_default(&sbwt, &lcs, sbwt.k());
        let vodbg = VoDbg::new(&sbwt, &pnsv_tuned);
        for kmer_end in 1..seqs[0].len() {
            let kmer = &seqs[0][kmer_end-1..=kmer_end];
            let node = vodbg.get_node(kmer).expect("Should exist.");
            assert_eq!(vodbg.indegree(node), 4);
            assert_eq!(vodbg.outdegree(node), 4);

            let one_mer = vodbg.contract_left(node, 1);
            assert_eq!(vodbg.indegree(one_mer), 4);
            assert_eq!(vodbg.outdegree(one_mer), 4);
        }
    }
}

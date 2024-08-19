//! de Bruijn graph operations using [SbwtIndex].

use crate::sbwt::SbwtIndex;
use crate::subsetseq::SubsetSeq;
use std::{io::Write, sync::{Arc, Mutex}};
use rayon::prelude::*;
use bitvec::bitvec;

/// A struct supporting de Bruijn graph operations on the spectrum of k-mers
/// encoded in an [SbwtIndex]. The graph is **node-centric**, meaning that the nodes
/// are the distinct k-mers in the index (not including the [dummy k-mers](crate::sbwt::SbwtIndex#sbwt-graph)) and there is an 
/// edge from x to y iff x[1..k) = y[0..k-1).
/// Edge (x,y) is labeled with the last character of y. Reverse complements are not modeled. 
/// The struct takes 2n bits of extra space on top of the SBWT, where n is the number of sets in the SBWT.
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct Dbg<'a, SS: SubsetSeq + Send + Sync> {
    sbwt: &'a SbwtIndex<SS>,
    k_minus_1_marks: bitvec::vec::BitVec,
    dummy_marks: bitvec::vec::BitVec,
}

/// A node in the de Bruijn graph.
#[derive(Copy, Clone, Eq, PartialEq, Debug, Hash)]
pub struct Node {
    pub id: usize,
}

/// An iterator over the nodes of the de Bruijn graph.
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct NodeIterator<'a> {
    colex: usize,
    dummy_marks: &'a bitvec::vec::BitVec,
}

impl<'a> Iterator for NodeIterator<'a> {
    type Item = Node;
    fn next(&mut self) -> Option<Self::Item> {
        while self.colex < self.dummy_marks.len() && self.dummy_marks[self.colex] {
            self.colex += 1;
        }
        let node = Node{id: self.colex};
        self.colex += 1; // Starting point for next iteration

        if node.id < self.dummy_marks.len(){
            Some(node)
        } else {
            None
        }
    }
}

impl<'a, SS: SubsetSeq + Send + Sync> Dbg<'a, SS> {

    /// Returns an iterator over all nodes of the de Bruijn graph.
    pub fn node_iterator(&'a self) -> NodeIterator<'a> {
        NodeIterator{colex: 0, dummy_marks: &self.dummy_marks}
    }

    /// An internal function for marking the dummy nodes in the SBWT.
    fn mark_dummies(sbwt: &SbwtIndex<SS>) -> bitvec::vec::BitVec {
        let mut dummy_marks = bitvec![0; sbwt.n_sets()];
        let mut dfs_stack = Vec::<(usize, usize)>::new(); // Node, depth
        dfs_stack.push((0,0)); // Colex rank of $, depth of $
        let mut outlabels = Vec::<u8>::new();
        while let Some((v, depth)) = dfs_stack.pop() { 
            outlabels.clear();
            sbwt.sbwt.append_set_to_buf(v, &mut outlabels);
            dummy_marks.set(v, true);

            if depth + 1 < sbwt.k() {
                for &c_idx in outlabels.iter() {
                    let u = sbwt.lf_step(v, c_idx as usize);
                    dfs_stack.push((u, depth + 1));
                }
            }
        }

        dummy_marks
    } 

    /// An internal function marking for each (k-1)-mer the smallest k-mer that has that (k-1)-mer as a suffix.
    fn mark_k_minus_1_mers(lcs: &crate::streaming_index::LcsArray, k: usize) -> bitvec::vec::BitVec {
        let mut k_minus_1_marks = bitvec![0; lcs.len()];
        for i in 0..lcs.len(){
            let len = lcs.access(i);
            if len < k - 1 {
                k_minus_1_marks.set(i,true);
            }
        }
        k_minus_1_marks
    }

    /// Returns next 1-bit to the right of i, or bv.len() if does not exist
    fn next_1_bit(bv: &bitvec::vec::BitVec, mut i: usize) -> usize {
        while i < bv.len() && !bv[i] {
            i += 1;
        }
        i
    }

    /// Initializes supports for de Bruijn graph operation based on the given [SbwtIndex].
    /// If the Lcs array of the SBWT is available, it can be given to significantly speed up construction.
    /// IMPORTANT: [select support][SbwtIndex::build_select()] must be built before calling this function. 
    pub fn new(sbwt: &'a SbwtIndex<SS>, lcs: Option<&crate::streaming_index::LcsArray>) -> Self {
        assert!(sbwt.sbwt.has_select_support());
        let k_minus_1_marks = match lcs {
            Some(lcs) => {
                log::info!("Building (k-1)-mer marks from LCS array");
                Self::mark_k_minus_1_mers(lcs, sbwt.k())
            }
            None => {
                log::info!("No LCS-array given. Building (k-1)-mer marks with column inversion.");
                sbwt.mark_k_minus_1_mers()
            }
        };
        let dummy_marks = Self::mark_dummies(sbwt);
        Self{sbwt, k_minus_1_marks, dummy_marks}
    }

    /// Push the k-mer string of the node to the given buffer.
    pub fn push_node_kmer(&self, node: Node, buf: &mut Vec<u8>) {
        assert!(!self.dummy_marks[node.id]);
        self.sbwt.push_kmer_to_vec(node.id, buf);
    }

    /// Get the k-mer string label of a node. To avoid memory allocation, check
    /// [Dbg::push_node_kmer].
    pub fn get_kmer(&self, node: Node) -> Vec<u8> {
        assert!(!self.dummy_marks[node.id]);
        let mut buf = Vec::<u8>::with_capacity(self.sbwt.k());
        self.push_node_kmer(node, &mut buf);
        buf
    }

    /// Get a handle to the node corresponding to the given k-mer, if exists in the graph.
    pub fn get_node(&self, kmer: &[u8]) -> Option<Node> {
        assert!(kmer.len() == self.sbwt.k());
        self.sbwt.search(kmer).map(|range| Node{id: range.start})
    }

    /// Returns the number of outgoing edges from the given node.
    pub fn outdegree(&self, node: Node) -> usize {
        assert!(!self.dummy_marks[node.id]);
        self.sbwt.sbwt.subset_size(self.get_representative_k_minus_1_mer(node).id)
    }

    /// Returns the number of incoming edges to the given node.
    pub fn indegree(&self, node: Node) -> usize {
        assert!(!self.dummy_marks[node.id]);
        match self.follow_inedge(node, 0) {
            Some(v) => self.k_minus_1_mer_freq(v),
            None => 0,
        }
    }
    
    fn get_representative_k_minus_1_mer(&self, node: Node) -> Node{
        assert!(!self.dummy_marks[node.id]);
        let mut v = Node{id: node.id};

        // Go to the smallest k-mer that has the same suffix as our node
        while !self.k_minus_1_marks[v.id] { // index 0 is always marked so we're good
            v.id -= 1;
        }
        v
    }

    /// For each outgoing edge from the given node, pushes to the output vector a pair
    /// (v, c), where v is the target node and c is the edge label.
    pub fn push_out_neighbors(&self, node: Node, output: &mut Vec<(Node, u8)>){
        assert!(!self.dummy_marks[node.id]);

        let rep = self.get_representative_k_minus_1_mer(node);

        for (i, &c) in self.sbwt.alphabet().iter().enumerate() {
            if self.sbwt.sbwt.set_contains(rep.id, i as u8) {
                let outnode = Node{id: self.sbwt.lf_step(rep.id, i)};
                output.push((outnode, c));
            }
        }
    }

    /// For each incoming edge to the given node, pushes to the output vector a pair
    /// (v, c), where v is the source node and c is the edge label.
    pub fn push_in_neighbors(&self, node: Node, output: &mut Vec<(Node, u8)>){
        assert!(!self.dummy_marks[node.id]);
        if let Some(v) = self.sbwt.inverse_lf_step(node.id) { // Predecessor
            if self.dummy_marks[v] { return; }
            let vrep = self.get_representative_k_minus_1_mer(Node{id: v}).id;
            let end = Self::next_1_bit(&self.k_minus_1_marks, vrep+1);
            let inlabel = self.get_last_character(node);
            (vrep..end).filter(|&i| !self.dummy_marks[i]).for_each(|i|{
                output.push((Node{id: i}, inlabel));
            });
        }
    }

    /// Gets the last character of the k-mer string of the given node.
    pub fn get_last_character(&self, node: Node) -> u8 {
        assert!(!self.dummy_marks[node.id]);
        self.sbwt.inlabel(node.id).unwrap() // Can unwrap because this is not a dummy node
    }

    /// Returns whether the given node has an outgoing edge labeled with `edge_label`.
    pub fn has_outlabel(&self, node: Node, edge_label: u8) -> bool {
        assert!(!self.dummy_marks[node.id]);
        let rep = self.get_representative_k_minus_1_mer(node);
        let c_idx = self.sbwt.char_idx(edge_label) as u8;
        self.sbwt.sbwt.set_contains(rep.id, c_idx)
    }

    /// Pushes the labels of all outgoing edges from the given node to the output vector.
    pub fn push_outlabels(&self, node: Node, output: &mut Vec<u8>) {
        assert!(!self.dummy_marks[node.id]);
        let rep = self.get_representative_k_minus_1_mer(node);
        self.sbwt.sbwt.append_set_to_buf(rep.id, output);
        for c in output.iter_mut() { // Map from 0123 to ACGT
            *c = self.sbwt.alphabet()[*c as usize];
        }
    }

    /// Follows the outgoing edge labeled with edge_label from the given node.
    /// Returns None if the edge does not exist.
    pub fn follow_outedge(&self, node: Node, edge_label: u8) -> Option<Node>{
        assert!(!self.dummy_marks[node.id]);
        if !self.has_outlabel(node, edge_label) {
            return None;
        }
        let rep = self.get_representative_k_minus_1_mer(node);
        Some(Node{id: self.sbwt.lf_step(rep.id, self.sbwt.char_idx(edge_label))})
    }

    // Takes in the representative (= colex smallest) node suffixed by the (k-1)-mer
    // Returns the number of full k-mers that have the same (k-1)-mer suffix as the given
    // representative (including itself if rep itself is a full k-mer).
    fn k_minus_1_mer_freq(&self, rep: Node) -> usize {
        assert!(!self.dummy_marks[rep.id]);
        assert!(self.k_minus_1_marks[rep.id]);
        let start = rep.id;
        let end = Self::next_1_bit(&self.k_minus_1_marks, start+1);
        self.dummy_marks[start..end].iter().filter(|b| !(**b)).count() // Number of 0-bits in range
    }

    /// Follows backward the incoming edge that comes from the i-th smallest k-mer
    /// (i ∈ [0, indegree(node)) in colexicographic order that has an outgoing edge to `node`. 
    /// Returns None if i ≥ indegree(node). 
    pub fn follow_inedge(&self, node: Node, i: usize) -> Option<Node>{
        assert!(!self.dummy_marks[node.id]);
        let v = self.sbwt.inverse_lf_step(node.id)?;
        if self.dummy_marks[v] { return None; }
        let v = Node {id: v};
        let vrep = self.get_representative_k_minus_1_mer(v).id;
        let end = Self::next_1_bit(&self.k_minus_1_marks, vrep+1);
        let mut non_dummies = 0_usize;

        // Return the position of the 0-bit with rank in_edge_number in the range
        for j in vrep..end {
            if !self.dummy_marks[j] {
                if non_dummies == i {
                    return Some(Node{id: j});
                }
                non_dummies += 1; 
            }
        }

        None
    }


    /// Writes the unitigs of the graph to the given writer in FASTA format.
    pub fn parallel_export_unitigs<W: Write + Send + Sync>(&self, fasta_out: W){
        let start_time = std::time::Instant::now();

        let unitig_id = 0_usize;
        let visited = bitvec![0; self.sbwt.n_sets()];

        let shared_data = Arc::new(Mutex::new((visited, fasta_out, unitig_id)));

        log::info!("Listing acyclic unitigs");
        self.node_iterator().filter(|&v| self.is_first_kmer_of_unitig(v)).par_bridge().for_each(|v| {
            let mut out_labels_buf = Vec::<u8>::new();
            let (nodes, unitig_string) = self.walk_unitig_from(v, &mut out_labels_buf);

            let (visited, fasta_out, unitig_id) = &mut *shared_data.lock().unwrap();
            for u in nodes {
                assert!(!visited[u.id]);
                visited.set(u.id, true);
            }

            write!(fasta_out, ">{}\n{}\n", unitig_id, String::from_utf8_lossy(&unitig_string)).unwrap();
            *unitig_id += 1;
        });

        log::info!("Listing cyclic unitigs");

        let mut out_labels_buf = Vec::<u8>::new();
        let (visited, fasta_out, unitig_id) = &mut *shared_data.lock().unwrap();

        // Only disjoint cyclic unitigs remain
        for v in self.node_iterator(){
            if visited[v.id] {
                continue;
            }

            out_labels_buf.clear();
            let (nodes, unitig_string) = self.walk_unitig_from(v, &mut out_labels_buf);
            for u in nodes {
                assert!(!visited[u.id]);
                visited.set(u.id, true);
            }

            write!(fasta_out, ">{}\n{}\n", unitig_id, String::from_utf8_lossy(&unitig_string)).unwrap();
            *unitig_id += 1;

        }

        let end_time = std::time::Instant::now();
        log::info!("Wrote all {} unitigs in {} seconds", unitig_id, (end_time - start_time).as_secs_f64());
    }

    fn is_first_kmer_of_unitig(&self, v: Node) -> bool {
        if self.indegree(v) > 1 {
            return true;
        }
        if let Some(u) = self.follow_inedge(v, 0) {
            self.outdegree(u) > 1
        } else {
            true
        }
    }

    // Returns the sequence of nodes and the label of the unitig
    // The out_labels_buf is working space for the function. Don't assume
    // anything about its contents when the function returns.
    fn walk_unitig_from(&self, mut v: Node, out_labels_buf: &mut Vec<u8>) -> (Vec<Node>, Vec<u8>){
        let v0 = v;
        let mut nodes = Vec::<Node>::new();
        nodes.push(v);
        
        let mut label = Vec::<u8>::new();
        self.push_node_kmer(v, &mut label); 

        while self.outdegree(v) == 1 {
            out_labels_buf.clear();
            self.push_outlabels(v, out_labels_buf);
            let c = out_labels_buf[0];
            v = self.follow_outedge(v, c).unwrap();
            if v != v0 && self.indegree(v) == 1 {
                label.push(c);
                nodes.push(v);
            } else { break; }
        }

        (nodes, label)
}

}


#[cfg(test)]
mod tests {
    use std::io::BufRead;

    use crate::{builder::{BitPackedKmerSorting, SbwtIndexBuilder}, util};
    use bitvec::prelude::*;
    use rand_chacha::rand_core::RngCore;
    use super::*;

    #[test]
    fn finimizer_paper_example_unitig_export(){
        // Note: this test does not cover the cyclic unitig case
        let seqs: Vec<&[u8]> = vec![b"GTAAGTCT", b"AGGAAA", b"ACAGG", b"GTAGG", b"AGGTA"];
        let mut seqs_copy: Vec<Vec<u8>> = seqs.iter().map(|x| x.to_vec()).collect(); // For verification in the end

        let (mut sbwt, lcs) = SbwtIndexBuilder::<BitPackedKmerSorting>::new().k(4).build_lcs(true).run_from_slices(seqs.as_slice());
        sbwt.build_select();
        let dbg = Dbg::new(&sbwt, lcs.as_ref());

        let mut output = Vec::<u8>::new();
        let out_cursor = std::io::Cursor::new(&mut output);
        dbg.parallel_export_unitigs(out_cursor);

        let mut unitigs: Vec<Vec<u8>> = vec![];
        for line in output.lines(){
            let line = line.unwrap();
            if !line.starts_with('>'){
                unitigs.push(line.as_bytes().to_vec());
            }
        }
        unitigs.sort();
        seqs_copy.sort();

        assert_eq!(unitigs, seqs_copy);

    }

    
    #[test]
    fn finimizer_paper_example_dbg_operations(){
        let seqs: Vec<&[u8]> = vec![b"GTAAGTCT", b"AGGAAA", b"ACAGG", b"GTAGG", b"AGGTA"];
        let (mut sbwt, lcs) = SbwtIndexBuilder::<BitPackedKmerSorting>::new().k(4).build_lcs(true).run_from_slices(seqs.as_slice());
        sbwt.build_select();

        let lcs = lcs.unwrap();
        let dbg = Dbg::new(&sbwt, Some(&lcs));
        let dbg_without_lcs = Dbg::new(&sbwt, None);

        let true_dummy_marks = bitvec![1,1,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0];
        assert_eq!(dbg.dummy_marks, true_dummy_marks);
        assert_eq!(dbg_without_lcs.dummy_marks, true_dummy_marks);

        let true_k_minus_1_marks = bitvec![1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1];
        assert_eq!(dbg.k_minus_1_marks, true_k_minus_1_marks);
        assert_eq!(dbg_without_lcs.k_minus_1_marks, true_k_minus_1_marks);

        // Get node 

        assert!(dbg.get_node(b"TTAT").is_none());
        let v = dbg.get_node(b"ACAG").unwrap();
        assert_eq!(v.id, 11);
        assert_eq!(dbg.outdegree(v), 1);

        // Push node kmer

        let mut buf = Vec::<u8>::new();
        dbg.push_node_kmer(v, &mut buf);
        assert_eq!(&buf, b"ACAG");

        // Has outlabel

        assert!(!dbg.has_outlabel(v , b'A'));
        assert!(!dbg.has_outlabel(v , b'C'));
        assert!(dbg.has_outlabel(v , b'G'));
        assert!(!dbg.has_outlabel(v , b'T'));

        // Follow outedge

        assert!(dbg.follow_outedge(v, b'A').is_none());
        let v = dbg.follow_outedge(v, b'G').unwrap();

        // Outdegree

        assert_eq!(dbg.outdegree(v), 2);

        // Out labels

        let mut outlabels = Vec::<u8>::new();
        dbg.push_outlabels(v, &mut outlabels);
        assert_eq!(outlabels, vec![b'A', b'T']);

        let v = dbg.get_node(b"TAGG").unwrap();
        let mut outlabels = Vec::<u8>::new();
        dbg.push_outlabels(v, &mut outlabels);
        assert_eq!(outlabels, vec![b'A', b'T']);

        // Out neighbors

        let v = dbg.get_node(b"CAGG").unwrap();
        let mut out_neighbors = Vec::<(Node, u8)>::new();
        dbg.push_out_neighbors(v, &mut out_neighbors);
        assert_eq!(out_neighbors, vec![(dbg.get_node(b"AGGA").unwrap(), b'A'), (dbg.get_node(b"AGGT").unwrap(), b'T')]);

        let mut out_neighbors = Vec::<(Node, u8)>::new();
        dbg.push_out_neighbors(dbg.get_node(b"GTCT").unwrap(), &mut out_neighbors);
        assert!(out_neighbors.is_empty());

        // Get inlabel

        assert_eq!(dbg.get_last_character(v), b'G');
        assert_eq!(dbg.get_last_character(Node{id:11}), b'G'); // ACAG

        // Indegree

        let v = dbg.get_node(b"AGGA").unwrap();
        assert_eq!(dbg.indegree(v), 2);

        let v = dbg.get_node(b"AGGT").unwrap();
        assert_eq!(dbg.indegree(v), 2);

        let v = dbg.get_node(b"TAGG").unwrap();
        assert_eq!(dbg.indegree(v), 1);

        assert_eq!(dbg.indegree(Node{id: 11}), 0); // ACAG

        // In neighbors 

        let v = dbg.get_node(b"AGGA").unwrap();
        let mut in_neighbors = Vec::<(Node, u8)>::new();
        dbg.push_in_neighbors(v, &mut in_neighbors);
        assert_eq!(in_neighbors, vec![(dbg.get_node(b"CAGG").unwrap(), b'A'), (dbg.get_node(b"TAGG").unwrap(), b'A')]);

        let v = dbg.get_node(b"AGGT").unwrap();
        let mut in_neighbors = Vec::<(Node, u8)>::new();
        dbg.push_in_neighbors(v, &mut in_neighbors);
        assert_eq!(in_neighbors, vec![(dbg.get_node(b"CAGG").unwrap(), b'T'), (dbg.get_node(b"TAGG").unwrap(), b'T')]);

        let mut in_neighbors = Vec::<(Node, u8)>::new();
        dbg.push_in_neighbors(Node{id: 11}, &mut in_neighbors); // ACAG
        assert!(in_neighbors.is_empty());

    }

    #[test]
    fn cyclic_unitigs_in_export(){
        use rand_chacha::ChaCha20Rng;
        use rand_chacha::rand_core::SeedableRng;

        let k = 10_usize;
        let mut rng = ChaCha20Rng::from_seed([123; 32]);
        let mut seqs = Vec::<Vec<u8>>::new();
        for _ in 0..10 {
            let seq: Vec<u8> = (0..3*k).map(|_| match rng.next_u32() % 4 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            }).collect();
            seqs.push(seq.clone());
        }

        let x0 = seqs[0][0..k].to_vec();
        seqs[0].extend(x0); // Make cyclic

        let x1 = seqs[1][0..k].to_vec();
        seqs[1].extend(x1); // Make cyclic

        let x2 = seqs[2][0..k].to_vec();
        seqs[2].extend(x2); // Make cyclic

        let (sbwt, lcs) = SbwtIndexBuilder::<BitPackedKmerSorting>::new().k(k).build_lcs(true).build_select_support(true).run_from_vecs(seqs.as_slice());
        let dbg = Dbg::new(&sbwt, lcs.as_ref());

        let mut unitig_ascii_out = Vec::<u8>::new();
        dbg.parallel_export_unitigs(std::io::Cursor::new(&mut unitig_ascii_out));
        let unitigs: Vec<Vec<u8>> = unitig_ascii_out.lines().map(|s| s.unwrap().as_bytes().to_owned()).filter(|s| s[0] != b'>').collect();

        let mut n_kmers = 0_usize;
        for unitig in unitigs.iter(){
            let self_overlap: bool = unitig[0..k-1] == unitig[unitig.len()-k+1..];
            eprintln!("{}", String::from_utf8_lossy(unitig));
            eprintln!("self_overlap: {}", self_overlap);
            for (i, kmer) in unitig.windows(k).enumerate() {
                assert!(dbg.get_node(kmer).is_some());
                let node = dbg.get_node(kmer).unwrap();
                let indeg = dbg.indegree(node);
                let outdeg = dbg.outdegree(node);

                if i == 0 { // Check that the unitig is maximal to the left
                    if self_overlap {
                        assert_eq!(indeg, 1);
                    } else if indeg == 1 { // Indeg 0 or >= 2 are always ok.
                        let pred = dbg.follow_inedge(node, 0).unwrap();
                        assert!(dbg.outdegree(pred) >= 2);
                    }
                }
                if i + k == unitig.len() { // Check that the unitig is maximal to the right
                    if self_overlap {
                        assert_eq!(outdeg, 1);
                    } else if outdeg == 1 { // Outdeg 0 or >= 2 are always ok.
                        let mut outlabels = Vec::<u8>::new();
                        dbg.push_outlabels(node, &mut outlabels);
                        let succ = dbg.follow_outedge(node, *outlabels.first().unwrap()).unwrap();
                        assert!(dbg.indegree(succ) >= 2);
                    }
                }
                if i > 0 {
                    assert_eq!(indeg, 1);
                }
                if i + k != unitig.len() {
                    assert_eq!(outdeg, 1);
                }
                n_kmers += 1;
            }
        }
        assert_eq!(n_kmers, sbwt.n_kmers());
    }

    #[test]
    fn randomized_test(){
        use rand_chacha::ChaCha20Rng;
        use rand_chacha::rand_core::SeedableRng;

        // Generate 1000 random k-mers using a seeded rng.
        let k = 5_usize;
        let mut rng = ChaCha20Rng::from_seed([123; 32]);

        let mut seqs = Vec::<Vec<u8>>::new();
        let mut seqs_hashset = std::collections::HashSet::<Vec<u8>>::new();
        for _ in 0..1000 {
            let kmer: Vec<u8> = (0..k).map(|_| match rng.next_u32() % 4 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            }).collect();
            seqs.push(kmer.clone());
            seqs_hashset.insert(kmer);
        }

        seqs.sort();
        seqs.dedup();

        let (sbwt, lcs) = SbwtIndexBuilder::<BitPackedKmerSorting>::new().k(k).build_lcs(true).build_select_support(true).run_from_vecs(seqs.as_slice());
        let dbg = Dbg::new(&sbwt, lcs.as_ref());

        
        { // Check that node iterator iterates all k-mers (tests node_iterator, get_kmer)
            let mut extracted_kmers: Vec<Vec<u8>> = dbg.node_iterator().map(|v| dbg.get_kmer(v)).collect();
            extracted_kmers.sort();
            eprintln!("{} {}", extracted_kmers.len(), seqs.len());
            assert_eq!(extracted_kmers, seqs);
        }

        { // Test get_node
            for kmer in seqs.iter() {
                assert_eq!(dbg.get_kmer(dbg.get_node(kmer).unwrap()), *kmer);
            }
            // Try to get a non-existent k-mer.
            assert!(dbg.get_node(b"XXXXX").is_none());
        }

        { // Test outdegree, has_outlabel, push_outlabels, follow_outedge, push_out_neighbors
            for v in dbg.node_iterator(){
                let kmer = dbg.get_kmer(v);
                eprintln!("Processing {} {}", v.id, String::from_utf8_lossy(&kmer));
                let mut true_outdegree = 0_usize;
                let mut true_outlabels = Vec::<u8>::new(); 
                let mut true_out_kmers = Vec::<Vec::<u8>>::new(); 
                for &c in util::DNA_ALPHABET.iter() {
                    let mut next = kmer[1..].to_vec();
                    next.push(c);
                    let has_c = seqs_hashset.contains(&next);
                    true_outdegree += has_c as usize;
                    assert_eq!(has_c, dbg.has_outlabel(v, c));
                    if has_c {
                        true_outlabels.push(c);
                        true_out_kmers.push(next.clone());
                        let next_node = dbg.follow_outedge(v, c).unwrap();
                        eprintln!("From {} to {}", v.id, next_node.id);
                        assert_eq!(dbg.get_kmer(next_node), next);
                    } else {
                        assert!(dbg.follow_outedge(v, c).is_none());
                    }
                }

                let mut outlabels = Vec::<u8>::new();
                dbg.push_outlabels(v, &mut outlabels);

                let mut outneighbors = Vec::<(Node, u8)>::new();
                dbg.push_out_neighbors(v, &mut outneighbors);

                assert_eq!(true_outdegree, dbg.outdegree(v));
                assert_eq!(outlabels, true_outlabels);

                assert_eq!(outneighbors.len(), true_outdegree);
                for i in 0..true_outdegree {
                    assert_eq!(dbg.get_kmer(outneighbors[i].0), true_out_kmers[i]);
                    assert_eq!(outneighbors[i].1, outlabels[i]);
                }

            }
        }

        { // Test indegree, get_last_character, follow_inedge, push_in_neighbors
            for v in dbg.node_iterator(){
                let kmer = dbg.get_kmer(v);
                eprintln!("Processing {} {}", v.id, String::from_utf8_lossy(&kmer));
                let mut true_indegree = 0_usize;
                let mut true_in_kmers = Vec::<Vec::<u8>>::new(); 
                let true_last_character = *kmer.last().unwrap();
                assert_eq!(true_last_character, dbg.get_last_character(v));
                for &c in util::DNA_ALPHABET.iter() {
                    let mut prev = vec![c];
                    prev.extend(&kmer[..k-1]);
                    let has_c = seqs_hashset.contains(&prev);
                    if has_c {
                        let prev_node = dbg.follow_inedge(v, true_indegree).unwrap();
                        assert_eq!(dbg.get_kmer(prev_node), prev);
                        true_in_kmers.push(prev);
                        true_indegree += 1;
                    }
                }
                assert!(dbg.follow_inedge(v, true_indegree).is_none()); // As specced in follow_inedge documentation comment
                assert_eq!(true_indegree, dbg.indegree(v));

                let mut in_neighbors = Vec::<(Node, u8)>::new();
                dbg.push_in_neighbors(v, &mut in_neighbors);

                for i in 0..true_indegree {
                    assert_eq!(dbg.get_kmer(in_neighbors[i].0), true_in_kmers[i]);
                    assert_eq!(in_neighbors[i].1, true_last_character);
                }
            }
        }

    }
}


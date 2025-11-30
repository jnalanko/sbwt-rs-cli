//! de Bruijn graph operations using [SbwtIndex].

use crate::{atomic_bitmap::AtomicBitmap, sbwt::SbwtIndex};
use crate::subsetseq::SubsetSeq;
use std::sync::atomic::AtomicUsize;
use std::sync::atomic::Ordering::{Acquire, Release};
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

impl Iterator for NodeIterator<'_> {
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
    pub fn new(sbwt: &'a SbwtIndex<SS>, lcs: Option<&crate::streaming_index::LcsArray>, n_threads: usize) -> Self
    where SS: Sync {
        assert!(sbwt.sbwt.has_select_support());
        let k_minus_1_marks = match lcs {
            Some(lcs) => {
                log::info!("Building (k-1)-mer marks from LCS array");
                Self::mark_k_minus_1_mers(lcs, sbwt.k())
            }
            None => {
                log::info!("No LCS-array given. Building (k-1)-mer marks with column inversion.");
                sbwt.mark_k_minus_1_mers(n_threads)
            }
        };

        log::info!("Marking dummy nodes");
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
        self.sbwt.sbwt.subset_size(self.get_suffix_group_start(node))
    }

    /// Returns the number of incoming edges to the given node.
    pub fn indegree(&self, node: Node) -> usize {
        assert!(!self.dummy_marks[node.id]);
        match self.follow_inedge(node, 0) {
            Some(v) => {
                let s = v.id; // This is the smallest non-dummy in-neighbor
                let e = Self::next_1_bit(&self.k_minus_1_marks, s + 1);
                e - s
            },
            None => 0,
        }
    }
    
    // Returns the colex rank of the smallest k-mer (possibly dummy)
    // that has the same suffix as the given node.
    fn get_suffix_group_start(&self, node: Node) -> usize {
        assert!(!self.dummy_marks[node.id]);
        let mut v = Node{id: node.id};

        // Go to the smallest k-mer that has the same suffix as our node
        while !self.k_minus_1_marks[v.id] { // index 0 is always marked so we're good
            v.id -= 1;
        }
        v.id
    }

    /// For each outgoing edge from the given node, pushes to the output vector a pair
    /// (v, c), where v is the target node and c is the edge label.
    pub fn push_out_neighbors(&self, node: Node, output: &mut Vec<(Node, u8)>){
        assert!(!self.dummy_marks[node.id]);

        let rep = self.get_suffix_group_start(node);

        for (i, &c) in self.sbwt.alphabet().iter().enumerate() {
            if self.sbwt.sbwt.set_contains(rep, i as u8) {
                let outnode = Node{id: self.sbwt.lf_step(rep, i)};
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
            let vrep = self.get_suffix_group_start(Node{id: v});
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
        let rep = self.get_suffix_group_start(node);
        let c_idx = self.sbwt.char_idx(edge_label) as u8;
        self.sbwt.sbwt.set_contains(rep, c_idx)
    }

    /// Pushes the labels of all outgoing edges from the given node to the output vector.
    pub fn push_outlabels(&self, node: Node, output: &mut Vec<u8>) {
        assert!(!self.dummy_marks[node.id]);
        let rep = self.get_suffix_group_start(node);
        self.sbwt.sbwt.append_set_to_buf(rep, output);
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
        let rep = self.get_suffix_group_start(node);
        Some(Node{id: self.sbwt.lf_step(rep, self.sbwt.char_idx(edge_label))})
    }

    /// Follows backward the incoming edge that comes from the i-th smallest k-mer
    /// (i ∈ [0, indegree(node)) in colexicographic order that has an outgoing edge to `node`. 
    /// Returns None if i ≥ indegree(node). 
    pub fn follow_inedge(&self, node: Node, i: usize) -> Option<Node>{
        assert!(!self.dummy_marks[node.id]);
        let v = self.sbwt.inverse_lf_step(node.id)?;
        let vrep = if self.dummy_marks[v] {
            // We can't use get_representative_k_minus_1_mer here because
            // that asserts that we're not calling it with a dummy.
            v
        } else {
            self.get_suffix_group_start(Node{id: v})
            // ^ This might be unnecessary since inverse_lf_step always takes
            // us to the smallest k-mer of a suffix group?
        };
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


    // Iterates unitigs in parallel and calls the given callback on each. The callback
    // is given the nodes in the unitig in order, and the string label of the unitig as &[u8].
    // The starting point of cyclic unitigs may vary nondeterministically between calls because
    // of parallelism.
    pub fn iter_unitigs_with_callback<F: Fn(&[Node], &[u8]) + Send + Sync>(&self, callback: F, n_threads: usize){

        let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();

        thread_pool.install(|| {
            rayon::scope(|s| {
                s.spawn(|_| {
                    let start_time = std::time::Instant::now();

                    let n_unitigs = AtomicUsize::new(0);
                    let visited = AtomicBitmap::new(self.sbwt.n_sets());

                    log::info!("Iterating acyclic unitigs");
                    (0..self.sbwt.n_sets())
                    .into_par_iter()
                    .filter(|&colex| !self.dummy_marks[colex])
                    .filter(|&colex| self.is_first_kmer_of_unitig(Node{id: colex}))
                    .for_each(|colex| {
                        let mut out_labels_buf = Vec::<u8>::new();
                        let (nodes, unitig_string) = self.walk_unitig_from(Node{id: colex}, &mut out_labels_buf);
                        callback(&nodes, &unitig_string);

                        // Mark the nodes as visited
                        for u in nodes {
                            assert!(!visited.get(u.id));
                            visited.set(u.id, true);
                        }

                        n_unitigs.fetch_add(1, Release);
                    });

                    log::info!("Iterating cyclic unitigs");

                    let mut visited = visited.into_bitvec(); // Non-atomic

                    let mut out_labels_buf = Vec::<u8>::new();

                    // Only disjoint cyclic unitigs remain. There are probably very few of these,
                    // so single-threaded iteration is ok.
                    let mut colex = 0_usize;
                    while colex < self.sbwt.n_sets() {
                        match visited[colex..].first_zero() {
                            Some(off) => colex += off,
                            None => break,
                        }
                        if !self.dummy_marks[colex] {
                            out_labels_buf.clear();
                            let (nodes, unitig_string) = self.walk_unitig_from(Node{id: colex}, &mut out_labels_buf);
                            callback(&nodes, &unitig_string);
                            for u in nodes {
                                assert!(!visited[u.id]);
                                visited.set(u.id, true);
                            }
                            n_unitigs.fetch_add(1, Release);
                        }
                        colex += 1;
                    }

                    let end_time = std::time::Instant::now();
                    log::info!("Iterated all {} unitigs in {} seconds", n_unitigs.load(Acquire), (end_time - start_time).as_secs_f64());
                });
            });
        });
    }

    /// Writes the unitigs of the graph to the given writer in FASTA format.
    /// Uses as many threads as are available with the Rayon initialization.
    pub fn parallel_export_unitigs<W: Write + Send + Sync>(&self, fasta_out: W, n_threads: usize){
        let unitig_id = 0_usize;
        let shared_data = Arc::new(Mutex::new((fasta_out, unitig_id)));

        let write_unitig_callback = move |_ : &[Node], unitig_string: &[u8]| {
            let (fasta_out, unitig_id) = &mut *shared_data.lock().unwrap();
            write!(fasta_out, ">{}\n{}\n", unitig_id, String::from_utf8_lossy(unitig_string)).unwrap();
            *unitig_id += 1;
        };

        self.iter_unitigs_with_callback(write_unitig_callback, n_threads);

    }

    pub fn is_first_kmer_of_unitig(&self, v: Node) -> bool {
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
    pub fn walk_unitig_from(&self, mut v: Node, out_labels_buf: &mut Vec<u8>) -> (Vec<Node>, Vec<u8>){
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

    pub fn is_dummy_colex_position(&self, pos: usize) -> bool {
        self.dummy_marks[pos]
    }

}


#[cfg(test)]
mod tests {
    use std::io::BufRead;

    use crate::{builder::{BitPackedKmerSortingDisk, SbwtIndexBuilder}, util};

    use crate::builder::BitPackedKmerSortingMem;

    use bitvec::prelude::*;
    use rand_chacha::rand_core::RngCore;
    use super::*;

    #[test]
    fn indegree_with_redundant_dummy_in_neighbor() {
        let k = 10;
        let s1: &[u8] =  b"TACGTACACA";
        let s2: &[u8] = b"ATACGTACACA";
        let seqs_1: &[&[u8]] = &[s1];
        let seqs_2: &[&[u8]] = &[s2];
        // We will merge SBWTs of these two sequences. The 10-mer TACGTACACA
        // requires a dummy path in seqs_1, but not in seqs_2. The merged SBWT
        // will have a dummy in-neighbor and a real in-neighbor. The dummy in-neighbor
        // should not be counted into the indegree.

        let (mut sbwt_1, lcs_1) = SbwtIndexBuilder::<BitPackedKmerSortingMem>::new().k(k).build_lcs(true).run_from_slices(seqs_1);
        let (mut sbwt_2, lcs_2) = SbwtIndexBuilder::<BitPackedKmerSortingMem>::new().k(k).build_lcs(true).run_from_slices(seqs_2);
        let merge_plan = crate::merge::MergeInterleaving::new(&sbwt_1, &sbwt_2, true, 3);
        let mut sbwt_merged = crate::merge::merge(sbwt_1, sbwt_2, merge_plan, 0, 3);
        sbwt_merged.build_select();

        let dbg = Dbg::new(&sbwt_merged, None, 3);
        let v = dbg.get_node(b"TACGTACACA".as_slice()).unwrap();
        assert_eq!(dbg.indegree(v), 1);
        assert_eq!(dbg.follow_inedge(v, 0).unwrap(), dbg.get_node(b"ATACGTACAC".as_slice()).unwrap());
    } 

    #[test]
    fn finimizer_paper_example_unitig_export(){
        // Note: this test does not cover the cyclic unitig case
        let seqs: Vec<&[u8]> = vec![b"GTAAGTCT", b"AGGAAA", b"ACAGG", b"GTAGG", b"AGGTA"];
        let mut seqs_copy: Vec<Vec<u8>> = seqs.iter().map(|x| x.to_vec()).collect(); // For verification in the end

        let (mut sbwt, lcs) = SbwtIndexBuilder::<BitPackedKmerSortingDisk>::new().k(4).build_lcs(true).run_from_slices(seqs.as_slice());
        sbwt.build_select();
        let dbg = Dbg::new(&sbwt, lcs.as_ref(), 3);

        let mut output = Vec::<u8>::new();
        let out_cursor = std::io::Cursor::new(&mut output);
        dbg.parallel_export_unitigs(out_cursor, 3);

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
        let (mut sbwt, lcs) = SbwtIndexBuilder::<BitPackedKmerSortingDisk>::new().k(4).build_lcs(true).run_from_slices(seqs.as_slice());
        sbwt.build_select();

        let lcs = lcs.unwrap();
        let dbg = Dbg::new(&sbwt, Some(&lcs), 3);
        let dbg_without_lcs = Dbg::new(&sbwt, None, 3);

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

        let (sbwt, lcs) = SbwtIndexBuilder::<BitPackedKmerSortingDisk>::new().k(k).build_lcs(true).build_select_support(true).run_from_vecs(seqs.as_slice());
        let dbg = Dbg::new(&sbwt, lcs.as_ref(), 3);

        let mut unitig_ascii_out = Vec::<u8>::new();
        dbg.parallel_export_unitigs(std::io::Cursor::new(&mut unitig_ascii_out), 3);
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

        let (sbwt, lcs) = SbwtIndexBuilder::<BitPackedKmerSortingDisk>::new().k(k).build_lcs(true).build_select_support(true).run_from_vecs(seqs.as_slice());
        let dbg = Dbg::new(&sbwt, lcs.as_ref(), 3);

        
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

    ///////////////////////////////////////////////////
    //
    // Same tests for bitpacked k-mer sorting in memory
    //
    ///////////////////////////////////////////////////


    #[test]
    fn finimizer_paper_example_unitig_export_mem(){
        // Note: this test does not cover the cyclic unitig case
        let seqs: Vec<&[u8]> = vec![b"GTAAGTCT", b"AGGAAA", b"ACAGG", b"GTAGG", b"AGGTA"];
        let mut seqs_copy: Vec<Vec<u8>> = seqs.iter().map(|x| x.to_vec()).collect(); // For verification in the end

        let (mut sbwt, lcs) = SbwtIndexBuilder::<BitPackedKmerSortingMem>::new().k(4).build_lcs(true).run_from_slices(seqs.as_slice());
        sbwt.build_select();
        let dbg = Dbg::new(&sbwt, lcs.as_ref(), 3);

        let mut output = Vec::<u8>::new();
        let out_cursor = std::io::Cursor::new(&mut output);
        dbg.parallel_export_unitigs(out_cursor, 3);

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
    fn finimizer_paper_example_dbg_operations_mem(){
        let seqs: Vec<&[u8]> = vec![b"GTAAGTCT", b"AGGAAA", b"ACAGG", b"GTAGG", b"AGGTA"];
        let (mut sbwt, lcs) = SbwtIndexBuilder::<BitPackedKmerSortingMem>::new().k(4).build_lcs(true).run_from_slices(seqs.as_slice());
        sbwt.build_select();

        let lcs = lcs.unwrap();
        let dbg = Dbg::new(&sbwt, Some(&lcs), 3);
        let dbg_without_lcs = Dbg::new(&sbwt, None, 3);

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
    fn cyclic_unitigs_in_export_mem(){
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

        let (sbwt, lcs) = SbwtIndexBuilder::<BitPackedKmerSortingMem>::new().k(k).build_lcs(true).build_select_support(true).run_from_vecs(seqs.as_slice());
        let dbg = Dbg::new(&sbwt, lcs.as_ref(), 3);

        let mut unitig_ascii_out = Vec::<u8>::new();
        dbg.parallel_export_unitigs(std::io::Cursor::new(&mut unitig_ascii_out), 3);
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
    fn randomized_test_mem(){
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

        let (sbwt, lcs) = SbwtIndexBuilder::<BitPackedKmerSortingMem>::new().k(k).build_lcs(true).build_select_support(true).run_from_vecs(seqs.as_slice());
        let dbg = Dbg::new(&sbwt, lcs.as_ref(), 3);


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

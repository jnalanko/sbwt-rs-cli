
use std::cmp::min;
use std::ops::Range;

use bitvec::prelude::*;
use rayon::iter::IntoParallelIterator;
use rayon::iter::IntoParallelRefIterator;
use rayon::iter::IndexedParallelIterator;
use rayon::iter::ParallelIterator;
use crate::compact_int_vector::CompactIntVector;
use crate::subsetseq::*;
use crate::sbwt::*;

type BitVec = bitvec::vec::BitVec<u64, Lsb0>;
type BitSlice = bitvec::slice::BitSlice<u64, Lsb0>;

/// An interleaving plan for combining two [SbwtIndex] structures with a set operation.
///
/// Two indexes of the same order k can be combined into an index of the [union](crate::merge), the
/// [intersection](crate::intersect) or the [difference](crate::difference) of their k-mer sets, without
/// going back to the underlying sequences. All three follow the same two steps:
/// 1. **Compute the interleaving plan.** Create a [MergeInterleaving] with [MergeInterleaving::new].
///    It records, for every node of the merged colex order, which of the two inputs it came from, and
///    serves as the blueprint of the operation.
/// 2. **Execute the operation.** Pass the interleaving and the two indexes to
///    [merge](crate::merge), [intersect](crate::intersect) or [difference](crate::difference).
///
/// The interleaving alone already determines the [size of the union](MergeInterleaving::union_size),
/// the [intersection](MergeInterleaving::intersection_size) and the
/// [difference](MergeInterleaving::difference_size), so these can be queried without materializing the
/// result index, which is what a Jaccard index computation needs.
///
/// The interleaving algorithm is an adaptation of the Wheeler graph merge algorithm described in
/// [*"Buffering updates enables efficient dynamic de Bruijn graphs"* (Alanko et al. 2021)](https://doi.org/10.1016/j.csbj.2021.06.047),
/// tailored specifically for the SBWT.
///
/// The intersection and the difference may turn a k-mer into a *source k-mer* (see
/// [here](crate#details-on-the-space-usage-of-the-index)) whose dummy predecessor chain is present in
/// neither input, and may leave dummy chains of the inputs leading nowhere. Both are repaired in
/// the interleaving of the inputs: the dummy nodes of the missing chains are found by searching the
/// prefixes of the source k-mers, and the new ones are inserted, so the result is always
/// structurally identical to an index built from scratch from the result k-mers.
/// For this, [intersect](crate::intersect) and [difference](crate::difference) require
/// [select support](SbwtIndex::build_select) on the first input, and return an error if it is missing.
///
/// The indexes and the interleaving are passed in as [Arc](std::sync::Arc)s so that the operation can
/// free them as early as possible when the caller holds no other reference, which lowers the memory
/// peak. Passing them as Arcs also allows the caller to keep holding on to them, in which case dropping
/// the Arcs inside the operation does not free the memory.
///
/// # Example
///
/// ```
/// use sbwt::*;
/// use std::sync::Arc;
///
/// let seqs1: Vec<&[u8]> = vec![b"AACTGACTGATCGTCTTGACTCGTTTATCTACGGT"];
/// let seqs2: Vec<&[u8]> = vec![b"ACTGACAGCTCTGCGATGCGACTCGTTTATCTAC"];
/// let k = 6;
///
/// // intersect and difference require select support on the first index
/// let (sbwt1, _) = BitPackedKmerSortingMem::new_from_slices(&seqs1, k)
///     .build_select_support(true).run();
/// let (sbwt2, _) = BitPackedKmerSortingMem::new_from_slices(&seqs2, k).run();
///
/// // The interleaving is the shared blueprint of all three set operations
/// let n_threads = 4;
/// let optimize_peak_ram = false; // Trades running time for a lower memory peak
/// let interleaving = MergeInterleaving::new(&sbwt1, &sbwt2, optimize_peak_ram, n_threads);
///
/// // Set sizes are available without building the result index
/// println!("|A ∪ B| = {}", interleaving.union_size());
/// println!("|A ∩ B| = {}", interleaving.intersection_size());
/// println!("|A \\ B| = {}", interleaving.difference_size());
///
/// let sbwt1 = Arc::new(sbwt1);
/// let sbwt2 = Arc::new(sbwt2);
/// let interleaving = Arc::new(interleaving);
/// let prefix_lookup_table_length = 4; // Of the result index
///
/// let union = merge(sbwt1.clone(), sbwt2.clone(), interleaving.clone(),
///     prefix_lookup_table_length, n_threads);
/// let isec = intersect(sbwt1.clone(), sbwt2.clone(), interleaving.clone(),
///     prefix_lookup_table_length, optimize_peak_ram, n_threads).unwrap();
/// let diff = difference(sbwt1, sbwt2, interleaving,
///     prefix_lookup_table_length, optimize_peak_ram, n_threads).unwrap();
///
/// // The results are ordinary SBWT indexes, ready to be queried or serialized
/// assert!(union.search(b"GACTCG").is_some());
/// ```
#[derive(Clone, Debug, Eq, PartialEq)]
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
    
    // Has one bit per colex position in the merged SBWT, marking
    // the positions whose k-mer has a different (k-1)-length suffix
    // than the previous k-mer in colex order. 
    // Edge case: is_leader[0] = 1.
    pub is_leader: BitVec,
}

trait ReadOnlyIntVector {
    fn get(&self, i: usize) -> usize;
    fn len(&self) -> usize;
    fn dollar_symbol() -> usize;
}

impl ReadOnlyIntVector for Vec<u8> {
    fn get(&self, i: usize) -> usize {
        self[i] as usize
    }

    fn len(&self) -> usize {
        self.len() 
    }

    fn dollar_symbol() -> usize {
        b'$' as usize
    }
}

impl ReadOnlyIntVector for CompactIntVector<3> {
    fn get(&self, i: usize) -> usize {
        self.get(i)
    }

    fn len(&self) -> usize {
        self.len() 
    }

    fn dollar_symbol() -> usize {
        0
    }
}

enum CharVector {
    ByteAlphabet(Vec<u8>),
    Compact(CompactIntVector<3>) // 3 bits per symbol for alphabet $,A,C,G,T
}

impl CharVector {
    fn init_with_byte_alphabet(len: usize) -> CharVector {
        CharVector::ByteAlphabet(vec![0; len])
    }
}

impl MergeInterleaving {

    /// Computes the merge interleaving between `index1` and `index2`. 
    /// The `optimize_peak_ram` flag enables optimizations to reduce the RAM peak at the expense of running time.
    /// `n_threads` is the number of parallel threads.
    pub fn new<SS: SubsetSeq + Send + Sync>(index1: &SbwtIndex::<SS>, index2: &SbwtIndex<SS>, optimize_peak_ram: bool, n_threads: usize) -> MergeInterleaving {

        use CharVector::*;

        let k = index1.k();
        assert_eq!(k, index2.k());

        log::info!("Computing merge interleaving: index1 has {} sets ({} k-mers), index2 has {} sets ({} k-mers), k={}",
            index1.n_sets(), index1.n_kmers(), index2.n_sets(), index2.n_kmers(), k);

        // Validate the SBWT invariant: sum of in-edge counts == n_sets - 1
        // (every non-root node has exactly one in-coming edge).
        // rank() takes 0-based alphabet indices (0..sigma), not ASCII character values.
        {
            let sigma = index1.alphabet().len();
            let in_edges1: usize = (0..sigma).map(|c| index1.sbwt().rank(c as u8, index1.n_sets())).sum();
            assert_eq!(in_edges1 + 1, index1.n_sets(),
                "index1 SBWT invariant violated: {} in-edge(s) + root != {} sets (diff: {}); \
                 {} nodes have no in-coming edge",
                in_edges1, index1.n_sets(),
                index1.n_sets() as isize - in_edges1 as isize - 1,
                index1.n_sets() - 1 - in_edges1);

            let in_edges2: usize = (0..sigma).map(|c| index2.sbwt().rank(c as u8, index2.n_sets())).sum();
            assert_eq!(in_edges2 + 1, index2.n_sets(),
                "index2 SBWT invariant violated: {} in-edge(s) + root != {} sets (diff: {}); \
                 {} nodes have no in-coming edge",
                in_edges2, index2.n_sets(),
                index2.n_sets() as isize - in_edges2 as isize - 1,
                index2.n_sets() - 1 - in_edges2);
        }

        // We invert the SBWTs column by column and maintain ranges
        // in both SBWTs that are so far equal. The ranges are in increasing
        // order and they partition the SBWTs, so it's enough to just store their
        // sizes in order. The sizes are stored in concatenated unary representations.
        // Empty ranges are allowed.

        let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();
        thread_pool.install(||{
            let mut leader_bits = None;

            // Initialize unary concatenations with empty ranges
            let mut s1 = bitvec![u64, Lsb0; 0; index1.n_sets()];
            let mut s2 = bitvec![u64, Lsb0; 0; index2.n_sets()];
            s1.push(true);
            s2.push(true);

            let (mut chars1, mut chars2, mut temp_char_buf_1, mut temp_char_buf_2 ) = if optimize_peak_ram {
                (Compact(index1.build_last_column_compact()), 
                 Compact(index2.build_last_column_compact()),
                 None, // Temp buffers allocated on demand
                 None  // Temp buffers allocated on demand
                )
            } else {
                (ByteAlphabet(index1.build_last_column()), 
                 ByteAlphabet(index2.build_last_column()),
                 Some(CharVector::init_with_byte_alphabet(index1.n_sets())),
                 Some(CharVector::init_with_byte_alphabet(index2.n_sets()))
                )
            };

            for round in 0..k {
                log::info!("Round {}/{}", round+1, k);

                // Split work into pieces for different threads
                log::debug!("Splitting work");
                let p1 = split_to_pieces_par(&s1, n_threads);
                let p2 = split_to_pieces_par(&s2, n_threads);
                assert_eq!(p1.len(), n_threads);
                assert_eq!(p2.len(), n_threads);
                // Zip pairs of tuples into 4-tuples
                let pieces = (0..n_threads).map(|i| (p1[i].0, p2[i].0, p1[i].1.clone(), p2[i].1.clone())).collect();

                log::debug!("Refining segmentation");
                let new_arrays = match (&chars1, &chars2) {
                    (ByteAlphabet(c1), ByteAlphabet(c2)) => {
                        refine_segmentation(s1, s2, c1, c2, pieces, round == k-1)
                    },
                    (Compact(c1), Compact(c2)) => {
                        refine_segmentation(s1, s2, c1, c2, pieces, round == k-1)
                    },
                    _ => panic!("Programmer messed up")
                };
                (s1, s2, leader_bits) = new_arrays;

                if round != k-1 {
                    log::debug!("Pushing labels forward in the SBWT graph");
                    if let (ByteAlphabet(ref mut c1), ByteAlphabet(ref mut c2), Some(ByteAlphabet(ref mut temp1)), Some(ByteAlphabet(ref mut temp2))) = (&mut chars1, &mut chars2, &mut temp_char_buf_1, &mut temp_char_buf_2) {
                        // High memory mode
                        index1.push_all_labels_forward(c1, temp1, n_threads);
                        std::mem::swap(c1, temp1);

                        index2.push_all_labels_forward(c2, temp2, n_threads);
                        std::mem::swap(c2, temp2);
                    } else if let (Compact(ref mut c1), Compact(ref mut c2), None, None) = (&mut chars1, &mut chars2, &mut temp_char_buf_1, &mut temp_char_buf_2) {
                        // Low memory mode
                        let mut temp1 = CompactIntVector::<3>::new(c1.len());
                        index1.push_all_labels_forward_compact(c1, &mut temp1, n_threads);
                        std::mem::swap(c1, &mut temp1);
                        drop(temp1);

                        let mut temp2 = CompactIntVector::<3>::new(c2.len());
                        index2.push_all_labels_forward_compact(c2, &mut temp2, n_threads);
                        std::mem::swap(c2, &mut temp2);
                        drop(temp2)
                    } else {
                        panic!("Programmer messed up");
                    }
                }
            }

            drop(temp_char_buf_1);
            drop(temp_char_buf_2);

            let leader_bits = leader_bits.unwrap(); // Computed in the last round
            log::debug!("Number of suffix groups: {}", leader_bits.count_ones());

            // Identify dummies in the merged SBWT

            log::debug!("Marking dummy nodes");
            let is_dummy = match (chars1, chars2) {
                (ByteAlphabet(c1), ByteAlphabet(c2)) => {
                    mark_dummy_nodes(&s1, &s2, &c1, &c2, n_threads)
                },
                (Compact(c1), Compact(c2)) => {
                    mark_dummy_nodes(&s1, &s2, &c1, &c2, n_threads)
                },
                _ => panic!("Programmer messed up")
            };

            log::info!("Number of dummies: {}", is_dummy.count_ones());

            MergeInterleaving { s1, s2, is_dummy, is_leader: leader_bits}
        })
    }

    /// Number of k-mers in the intersection of the two SBWTs.
    pub fn intersection_size(&self) -> usize {
        assert_eq!(self.s1.len(), self.s2.len());
        let s1w = self.s1.as_raw_slice();
        let s2w = self.s2.as_raw_slice();
        let dw = self.is_dummy.as_raw_slice();
        let mut ans = 0_usize;
        for i in 0..s1w.len() {
            // intersection = s1 & s2, excluding dummies
            ans += (s1w[i] & s2w[i] & !dw[i]).count_ones() as usize;
        }
        ans
    }

    /// Number of k-mers in `index1` that are **not** in `index2` (set difference index1 \ index2).
    pub fn difference_size(&self) -> usize {
        assert_eq!(self.s1.len(), self.s2.len());
        let s1w = self.s1.as_raw_slice();
        let s2w = self.s2.as_raw_slice();
        let dw = self.is_dummy.as_raw_slice();
        let mut ans = 0_usize;
        for i in 0..s1w.len() {
            // difference = s1 & !s2, excluding dummies
            ans += (s1w[i] & !s2w[i] & !dw[i]).count_ones() as usize;
        }
        ans
    }

    /// Number of k-mers in the union of the two SBWTs.
    pub fn union_size(&self) -> usize {
        assert_eq!(self.s1.len(), self.s2.len());
        let s1w = self.s1.as_raw_slice();
        let s2w = self.s2.as_raw_slice();
        let dw = self.is_dummy.as_raw_slice();
        let mut ans = 0_usize;
        for i in 0..s1w.len() {
            // union = s1 | s2, excluding dummies
            ans += ((s1w[i] | s2w[i]) & !dw[i]).count_ones() as usize;
        }
        ans
    }
}

/// An interleaving of `index1`, `index2` and a third set of nodes: the interleaving that the result
/// of an intersection, a difference or a filter with repaired dummy chains consists of, where the
/// third set (s3) is the dummy nodes of those chains.
///
/// For every colex position `i`:
/// * `s1[i]` — position belongs to `index1`,
/// * `s2[i]` — position belongs to `index2`,
/// * `s3[i]` — position belongs to the third set,
/// * `is_leader[i]` — position is the colex-minimum node of its (k-1)-suffix group.
///
/// For an intersection, a position is a **result node** iff `(s1[i] && s2[i]) || s3[i]` (except
/// for dead dummy chains), and the edge-bits at a result node are `s3_or[c] || (s1_or[c] &&
/// s2_or[c])`, where the ORs are accumulated over all positions in the same (k-1)-suffix group.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ThreeWayInterleaving {
    pub s1: BitVec,
    pub s2: BitVec,
    pub s3: BitVec,
    pub is_leader: BitVec,
}

/// The (k-1)-suffix groups of the result of a set operation, as [Dbg::new_with_k_minus_1_marks]
/// takes them: one bit per node of the result, set for the first node of each group. The result
/// consists of the positions of a merge plan ([MergeInterleaving] or [ThreeWayInterleaving]) that
/// are not `removed` (None if none are), in order, and `is_leader` are the starts of the groups of
/// that merge plan. A group of the result consists of its nodes in one group of the merge plan,
/// since both are ranges in colex order, so a node of the result starts a group if no node of the
/// result comes before it in its group of the merge plan.
///
/// [Dbg::new_with_k_minus_1_marks]: crate::dbg::Dbg::new_with_k_minus_1_marks
pub fn result_k_minus_1_marks(is_leader: &BitVec, removed: Option<&BitVec>, n_threads: usize) -> bitvec::vec::BitVec {
    let n = is_leader.len();
    assert!(removed.is_none_or(|removed| removed.len() == n));
    assert!(n == 0 || is_leader[0], "the root does not start a group");
    let in_result = |p: usize| removed.is_none_or(|removed| !removed[p]);
    let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();
    let parts: Vec<BitVec> = thread_pool.install(|| crate::util::segment_range(0..n, n_threads).into_par_iter().map(|range| {
        let mut marks = BitVec::new();
        if range.is_empty() {
            return marks;
        }
        // Whether a node of the result comes before the current position in its group: for the
        // first position, from the positions before it in its group (a group has a few nodes)
        let mut result_in_group = false;
        let mut p = range.start;
        while !is_leader[p] {
            p -= 1;
            if in_result(p) {
                result_in_group = true;
                break;
            }
        }
        for p in range {
            if is_leader[p] {
                result_in_group = false;
            }
            if in_result(p) {
                marks.push(!result_in_group);
                result_in_group = true;
            }
        }
        marks
    }).collect());
    let marks = crate::util::parallel_bitvec_concat(parts);
    // The words of a Dbg are usize
    let len = marks.len();
    let mut marks = bitvec::vec::BitVec::from_vec(marks.into_vec().into_iter().map(|w| w as usize).collect());
    marks.truncate(len);
    marks
}

fn mark_dummy_nodes<V: ReadOnlyIntVector + Send + Sync>(s1: &BitVec, s2: &BitVec, chars1: &V, chars2: &V, n_threads: usize) -> BitVec {
    assert_eq!(s1.len(), s2.len());
    let merged_len = s1.len();

    let merged_piece_ranges = crate::util::segment_range(0..merged_len, n_threads);
    let s1_piece_popcounts: Vec<usize> = merged_piece_ranges.par_iter().map(|range| s1[range.clone()].count_ones()).collect();
    let s2_piece_popcounts: Vec<usize> = merged_piece_ranges.par_iter().map(|range| s2[range.clone()].count_ones()).collect();
    let is_dummy_pieces = (0..merged_piece_ranges.len()).into_par_iter().map(|piece_idx| {
        let colex_range = &merged_piece_ranges[piece_idx];
        let mut is_dummy_bits = bitvec![u64, Lsb0 ;];
        is_dummy_bits.resize(colex_range.len(), false);

        let mut c1_idx: usize = s1_piece_popcounts[..piece_idx].iter().sum(); // Skip over previous pieces
        let mut c2_idx: usize = s2_piece_popcounts[..piece_idx].iter().sum(); // Skip over previous pieces
        for colex in colex_range.clone() {
            let d1 = s1[colex] && c1_idx < chars1.len() && chars1.get(c1_idx) == V::dollar_symbol();
            let d2 = s2[colex] && c2_idx < chars2.len() && chars2.get(c2_idx) == V::dollar_symbol();

            let rel_colex = colex - colex_range.start;
            is_dummy_bits.set(rel_colex, d1 || d2); 

            c1_idx += s1[colex] as usize;
            c2_idx += s2[colex] as usize;
        }
        is_dummy_bits 
    }).collect();

    crate::util::parallel_bitvec_concat(is_dummy_pieces)

}

// Computes per-piece ranges aligned to group leaders and per-piece popcounts for s1, s2, s3
// in a ThreeWayInterleaving. Returns at most n_threads pieces, all non-empty.
pub(super) fn compute_piece_ranges_three_way(
    n_threads: usize,
    interleaving: &ThreeWayInterleaving,
) -> (Vec<Range<usize>>, Vec<usize>, Vec<usize>, Vec<usize>) {
    let ranges = leader_aligned_ranges(&interleaving.is_leader, n_threads);
    let s1_pc: Vec<usize> = ranges.par_iter().map(|r| interleaving.s1[r.clone()].count_ones()).collect();
    let s2_pc: Vec<usize> = ranges.par_iter().map(|r| interleaving.s2[r.clone()].count_ones()).collect();
    let s3_pc: Vec<usize> = ranges.par_iter().map(|r| interleaving.s3[r.clone()].count_ones()).collect();
    (ranges, s1_pc, s2_pc, s3_pc)
}

// Whether `merged_colex` is a result node of the intersection (`difference` false) or the set
// difference (`difference` true) built from the three-way interleaving: `(s1 && s2 && !dead) || s3`
// or `(s1 && !s2 && !dead) || s3`. `s1_colex` is the number of s1 positions before `merged_colex`.
//
// `redundant_s1` (in index1 colex space) marks dead dummy-chain nodes: those are not result nodes
// through the s1-side predicate (they can still be result nodes through `s3` when a repaired dummy
// chain needs the same node). Every pass that needs the result nodes uses this function, so
// that they all agree.
#[inline]
pub(super) fn is_result_node(
    interleaving: &ThreeWayInterleaving,
    merged_colex: usize,
    s1_colex: usize,
    difference: bool,
    redundant_s1: Option<&BitVec>,
) -> bool {
    let s1_side = if difference {
        interleaving.s1[merged_colex] && !interleaving.s2[merged_colex]
    } else {
        interleaving.s1[merged_colex] && interleaving.s2[merged_colex]
    };
    (s1_side && redundant_s1.map_or(true, |red| !red[s1_colex])) || interleaving.s3[merged_colex]
}

// Counts the result nodes of each piece (see is_result_node).
pub(super) fn count_result_nodes_per_piece(
    merged_piece_ranges: &[Range<usize>],
    interleaving: &ThreeWayInterleaving,
    difference: bool,
    redundant_s1: Option<&BitVec>,
    s1_pc: &[usize],
) -> Vec<usize> {
    merged_piece_ranges.par_iter().enumerate().map(|(piece_idx, colex_range)| {
        let mut s1_colex: usize = s1_pc[..piece_idx].iter().sum();
        let mut count = 0usize;
        for merged_colex in colex_range.clone() {
            count += is_result_node(interleaving, merged_colex, s1_colex, difference, redundant_s1) as usize;
            s1_colex += interleaving.s1[merged_colex] as usize;
        }
        count
    }).collect()
}

// Assumes there is at least one 1-bit
fn leading_zeros(s: &BitSlice) -> usize {
    s.first_one().unwrap()
}

// Assumes that seq is a string with some number of c in the beginning (possibly none)
// followed by a tail of non-c characters which are all larger than c.
// Returns the number of c in the beginning.
// Falls back to binary search for long sequences
#[inline]
fn run_length_in_sorted_seq<V: ReadOnlyIntVector + Send + Sync>(seq: &V, range: Range<usize>, c: u8) -> usize {
    if range.is_empty() {
        return 0;
    }
    if range.len() > 200 {
        // Binary search the first element that is larger than c
        crate::util::binary_search_leftmost_that_fulfills_pred(|i| i, |i| seq.get(i+range.start) as u8 > c, range.len())
    } else {
        // Linear scan
        let mut i = 0;
        while i < range.len() && seq.get(range.start + i) as u8 == c {
            i += 1;
        }
        i
    }
}

// TODO: In a profile of themisto2's merge on E. coli assemblies, about 8% of the time was in the
// BitVec::resize of this function, called from refine_piece for every unary code. A larger
// pre-allocated capacity of the outputs did not change the time, so the cost is in writing the
// zeros, not in reallocating. Writing the unary codes a word at a time (e.g. only advancing the
// length over zeros that are already there, or appending whole words) might be faster.
#[inline]
fn zero_extend(v: &mut BitVec, howmany: usize) {
    v.resize(v.len() + howmany, false);
}


// Inner per-piece worker for refine_segmentation.
// Takes many arguments because it is a tight inner-loop helper; bundling them into
// structs would add boilerplate without improving clarity.
#[allow(clippy::too_many_arguments)] 
fn refine_piece<V: ReadOnlyIntVector + Send + Sync>(s1: &BitVec, s2: &BitVec, chars1: &V, chars2: &V, mut c1_i: usize, mut c2_i: usize, s1_range: Range<usize>, s2_range: Range<usize>, last_round: bool) 
-> (BitVec, BitVec, Option<BitVec>) {

    // Pre-allocate with estimated capacity based on input size
    let est_cap = s1_range.len();
    let mut out1 = BitVec::with_capacity(est_cap);
    let mut out2 = BitVec::with_capacity(est_cap);

    let mut leader_bits = if last_round { BitVec::with_capacity(est_cap) } else { BitVec::new() };

    // c1_i and c2_i are current indices in chars1 and chars2 respectively
    // s1_i and s2_i are current indices in s1 and s2 respectively
    let mut s1_i = s1_range.start;
    let mut s2_i = s2_range.start;

    while s1_i < s1_range.end {

        assert!(s2_i < s2_range.end);

        let len1 = leading_zeros(&s1[s1_i..]);
        let len2 = leading_zeros(&s2[s2_i..]);

        let c1_end = c1_i + len1; // One past the end
        let c2_end = c2_i + len2; // One past the end

        let mut is_leader = true;
        while c1_i < c1_end || c2_i < c2_end {
            let c1 = if c1_i == c1_end { u8::MAX } else { chars1.get(c1_i) as u8 };
            let c2 = if c2_i == c2_end { u8::MAX } else { chars2.get(c2_i) as u8 };
            let c = min(c1,c2);

            let r1 = run_length_in_sorted_seq(chars1, c1_i..c1_end, c);
            let r2 = run_length_in_sorted_seq(chars2, c2_i..c2_end, c);

            if last_round {
                // We know that r1 and r2 are at most 1 -> no need to have a full unary code.
                // Let's not assert this here because this is a tight inner loop.
                out1.push(r1 > 0);
                out2.push(r2 > 0);
            } else {
                // Write r1 and r2 in unary
                zero_extend(&mut out1, r1);
                zero_extend(&mut out2, r2);
                // Terminate unary representations
                out1.push(true);
                out2.push(true);
            }

            // Advance indexes in chars
            c1_i += r1;
            c2_i += r2;

            if last_round {
                leader_bits.push(is_leader);
            }

            is_leader = false;
        }

        assert_eq!(c1_i, c1_end);
        assert_eq!(c2_i, c2_end);

        s1_i += len1 + 1;
        s2_i += len2 + 1;
    }

    assert_eq!(s1_i, s1_range.end);
    assert_eq!(s2_i, s2_range.end);

    (out1, out2, if last_round { Some(leader_bits) } else { None })
}

// Helper function for construction.
// Input pieces are pairs (start in chars1, start in chars2, range in s1, range in s2)
// Input piece are for parallelism: one piece to work on for each thread.
// Returns the new segmentation in concatenate unary form.
// On the last round the output is different: now the two output bit vectors would have
// only 0 and 1 encoded in unary ("1" and "01"). Instead we write just the bits 0 and 1.
// On the last round the function also returns the leader bit vector, which marks
// the smallest k-mer in each group of k-mers with the same suffix of length (k-1).
// This function runs in parallel, so a rayon thread pool must be initialized.
fn refine_segmentation<V: ReadOnlyIntVector + Send + Sync>(s1: BitVec, s2: BitVec, chars1: &V, chars2: &V, input_pieces: Vec<(usize, usize, Range<usize>, Range<usize>)>, last_round: bool) -> (BitVec, BitVec, Option<BitVec>) {
    let output_pieces: Vec<(BitVec, BitVec, Option<BitVec>)> = input_pieces.par_iter().map(|piece| {
        let (c1_i, c2_1, s1_range, s2_range) = piece;
        refine_piece(&s1, &s2, chars1, chars2, *c1_i, *c2_1, s1_range.clone(), s2_range.clone(), last_round)
    }).collect();

    // Free memory
    drop(s1);
    drop(s2);

    // Turn output pieces (a vec of triples) into triple of vecs 
    let mut new_s1_pieces = vec![];
    let mut new_s2_pieces = vec![];
    let mut leader_pieces = if last_round {Some(vec![])} else {None};
    for (a,b,c) in output_pieces.into_iter() {
        new_s1_pieces.push(a);
        new_s2_pieces.push(b);
        if last_round {
            leader_pieces.as_mut().unwrap().push(c.unwrap());
        }
    }
    log::debug!("Concatenating pieces");
    let new_s1 = crate::util::parallel_bitvec_concat(new_s1_pieces);
    let new_s2 = crate::util::parallel_bitvec_concat(new_s2_pieces);
    let new_leader_bits = leader_pieces.map(crate::util::parallel_bitvec_concat);

    (new_s1, new_s2, new_leader_bits)
}

// Parallel version of split_to_pieces (see mod tests for a single-threaded reference implementation)
// Returns a segmentation s[l_1..r_1), s[l_2..r_2], ... such that
// each segment ends in a 1-bit and has an approximately equal number of 1-bits.
// Also returns the number of 0-bits before each segment.
// Runs on the caller's rayon thread pool — do not create a pool here.
pub(super) fn split_to_pieces_par(s: &BitSlice, n_pieces: usize) -> Vec<(usize, Range<usize>)> {
    const BLOCK_SIZE: usize = 1024;

    // Strategy: compute block popcounts in parallel, then locate the blocks where the
    // pieces start, and do bit-by-bit counting to find the precise start points of pieces.

    assert!(n_pieces > 0);
    if !s.is_empty() {
        // The last bit should always be 1
        assert!(s.last().unwrap() == true);
    }

    let blocks = s.chunks(BLOCK_SIZE).collect::<Vec::<&BitSlice>>();
    let block_popcounts: Vec<usize> = blocks.par_iter().map(|block| block.count_ones()).collect();
    let total_popcount: usize = block_popcounts.iter().sum();
    let ones_per_piece = total_popcount.div_ceil(n_pieces); // Last piece may have fewer

    let mut starts : Vec<usize> = vec![0]; // Init with the start of the first piece
    let mut n_zeros_before_piece: Vec<usize> = vec![0]; // Init with #zeros before the first piece

    let mut n_ones = 0_usize;
    let mut n_zeros = 0_usize;

    // Let N be the number of ones in a piece. The i-th block piece just after
    // the one-bit with zero-based rank N*i - 1. That is, if N = 10, then the fifth
    // piece starts at the one-bit with rank 49 (= the 50th 1-bit)
    for (block_idx, block) in blocks.iter().enumerate() {
        while starts.len() < n_pieces && n_ones + block_popcounts[block_idx] >= starts.len() * ones_per_piece {
            // The check for starts.len() < n_pieces is to avoid creating an empty piece after
            // the last one in case the total popcount is divisible by n_pieces.

            // the 1-bit just before the start the next piece is in this block.
            // Find where it is
            let mut n_ones_precise = n_ones;
            let mut n_zeros_precise = n_zeros;
            let target = starts.len() * ones_per_piece;
            let mut i = 0_usize;
            loop {
                n_ones_precise += block[i] as usize;
                n_zeros_precise += 1 - (block[i] as usize);
                if n_ones_precise == target {
                    break;
                } else {
                    i += 1;
                }
            }
            assert_eq!(n_ones_precise, target);
            starts.push(n_ones + n_zeros + i + 1);
            n_zeros_before_piece.push(n_zeros_precise);
        }
        n_ones += block_popcounts[block_idx];
        n_zeros += block.len() - block_popcounts[block_idx];
    }
    assert_eq!(starts.len(), n_zeros_before_piece.len());
    assert!(!starts.is_empty());
    while starts.len() < n_pieces {
        // Add empty pieces to the end
        starts.push(s.len());
        n_zeros_before_piece.push(s.len() - total_popcount);
    }

    let mut pieces: Vec<(usize, Range<usize>)> = vec![];
    assert_eq!(starts.len(), n_pieces);
    starts.push(s.len()); // End sentinel for the end of the last range
    for i in 0..n_pieces {
        pieces.push((n_zeros_before_piece[i], starts[i]..starts[i+1]));
    }

    pieces
}


// Splits 0..is_leader.len() into at most n_pieces non-empty ranges of roughly equal length,
// each starting at a group leader. There may be fewer than n_pieces ranges, e.g. when there
// are fewer positions than pieces, so callers must iterate over the returned ranges.
pub(super) fn leader_aligned_ranges(is_leader: &BitSlice, n_pieces: usize) -> Vec<Range<usize>> {
    let len = is_leader.len();
    let mut starts: Vec<usize> = crate::util::segment_range(0..len, n_pieces).into_iter().map(|r| r.start).collect();
    // Move each start forward to the next leader. This keeps the starts non-decreasing.
    // The first start is 0, which is always a leader.
    for s in starts.iter_mut().skip(1) {
        while *s < len && !is_leader[*s] {
            *s += 1;
        }
    }
    starts.push(len);
    starts.windows(2).map(|w| w[0]..w[1]).filter(|r| !r.is_empty()).collect()
}

// Aligns piece boundaries to group leaders and computes per-piece popcounts for s1 and s2.
// Returns at most n_threads pieces, all non-empty.
pub(super) fn compute_piece_ranges(n_threads: usize, interleaving: &MergeInterleaving) -> (Vec<Range<usize>>, Vec<usize>, Vec<usize>) {
    let merged_piece_ranges = leader_aligned_ranges(&interleaving.is_leader, n_threads);

    let s1_piece_popcounts: Vec<usize> = merged_piece_ranges.par_iter().map(|range| interleaving.s1[range.clone()].count_ones()).collect();
    let s2_piece_popcounts: Vec<usize> = merged_piece_ranges.par_iter().map(|range| interleaving.s2[range.clone()].count_ones()).collect();
    (merged_piece_ranges, s1_piece_popcounts, s2_piece_popcounts)
}


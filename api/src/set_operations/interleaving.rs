
use std::cmp::min;
use std::ops::Range;

use bitvec::prelude::*;
use rayon::iter::IntoParallelIterator;
use rayon::iter::IntoParallelRefIterator;
use rayon::iter::ParallelIterator;
use crate::compact_int_vector::CompactIntVector;
use crate::subsetseq::*;
use crate::sbwt::*;

type BitVec = bitvec::vec::BitVec<u64, Lsb0>;
type BitSlice = bitvec::slice::BitSlice<u64, Lsb0>;

/// An interleaving plan for [merging](merge) two [SbwtIndex] structures.
///
/// To merge two `SbwtIndex` structures, follow these steps:
/// 1. **Compute the interleaving plan.** Create a [MergeInterleaving] instance using [MergeInterleaving::new]. This interleaving serves as a blueprint for how the two SBWTs will be merged. It can also be queried to compute the size of the [intersection](MergeInterleaving::intersection_size) or the [union](MergeInterleaving::union_size) of the k-mer sets in the SBWTs.
/// 2. **Execute the merge.** Pass the interleaving and the two SBWTs to the [merge] function.
///
/// The merge algorithm is an adaptation of the Wheeler graph merge algorithm described in
/// [*"Buffering updates enables efficient dynamic de Bruijn graphs"* (Alanko et al. 2021)](https://doi.org/10.1016/j.csbj.2021.06.047),
/// tailored specifically for the SBWT.
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

/// A three-way interleaving of `index1`, `index2`, and a small auxiliary SBWT (`aux`).
///
/// Used when some result k-mers become new source nodes and need fresh dummy chains built via
/// an auxiliary SBWT. Replaces three separate 2-way interleaving computations (plus two
/// intermediate merges) with a single O((N1+N2+Naux)·k) pass.
///
/// For every colex position `i` in the three-way merged space:
/// * `s1[i]` — position belongs to `index1`,
/// * `s2[i]` — position belongs to `index2`,
/// * `s3[i]` — position belongs to `aux`,
/// * `is_leader[i]` — position is the colex-minimum node of its (k-1)-suffix group.
///
/// A position is a **result node** iff `(s1[i] && s2[i]) || s3[i]`.
/// Edge-bits at a result node are `(s1_or[c] || s3_or[c]) && (s2_or[c] || s3_or[c])`,
/// where the ORs are accumulated over all positions in the same (k-1)-suffix group.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ThreeWayInterleaving {
    pub s1: BitVec,
    pub s2: BitVec,
    pub s3: BitVec,
    pub is_leader: BitVec,
}

impl ThreeWayInterleaving {
    /// Compute the three-way interleaving of `index1`, `index2`, and `aux`.
    /// All three indices must share the same `k`.
    /// If `optimize_peak_ram` is true, uses compact 3-bit character encoding to reduce memory usage.
    pub fn new<SS: SubsetSeq + Send + Sync>(
        index1: &SbwtIndex<SS>,
        index2: &SbwtIndex<SS>,
        aux: &SbwtIndex<SS>,
        optimize_peak_ram: bool,
        n_threads: usize,
    ) -> ThreeWayInterleaving {


        let k = index1.k();
        assert_eq!(k, index2.k());
        assert_eq!(k, aux.k());

        log::info!(
            "Computing three-way intersection interleaving: index1={} sets, index2={} sets, aux={} sets, k={}",
            index1.n_sets(), index2.n_sets(), aux.n_sets(), k
        );

        let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();
        thread_pool.install(|| {
            let mut leader_bits: Option<BitVec> = None;

            // Start each segmentation with a single group containing all nodes of that index.
            let mut s1 = bitvec![u64, Lsb0; 0; index1.n_sets()];
            let mut s2 = bitvec![u64, Lsb0; 0; index2.n_sets()];
            let mut s3 = bitvec![u64, Lsb0; 0; aux.n_sets()];
            s1.push(true);
            s2.push(true);
            s3.push(true);

            use CharVector::*;
            let (mut chars1, mut chars2, mut chars3, mut temp_char_buf_1, mut temp_char_buf_2, mut temp_char_buf_3) = if optimize_peak_ram {
                (Compact(index1.build_last_column_compact()),
                 Compact(index2.build_last_column_compact()),
                 Compact(aux.build_last_column_compact()),
                 None,
                 None,
                 None,
                )
            } else {
                (ByteAlphabet(index1.build_last_column()),
                 ByteAlphabet(index2.build_last_column()),
                 ByteAlphabet(aux.build_last_column()),
                 Some(CharVector::init_with_byte_alphabet(index1.n_sets())),
                 Some(CharVector::init_with_byte_alphabet(index2.n_sets())),
                 Some(CharVector::init_with_byte_alphabet(aux.n_sets())),
                )
            };

            for round in 0..k {
                log::info!("Three-way interleaving round {}/{}", round + 1, k);

                let p1 = split_to_pieces_par(&s1, n_threads);
                let p2 = split_to_pieces_par(&s2, n_threads);
                let p3 = split_to_pieces_par(&s3, n_threads);

                let pieces: Vec<_> = (0..n_threads)
                    .map(|i| (p1[i].0, p2[i].0, p3[i].0,
                              p1[i].1.clone(), p2[i].1.clone(), p3[i].1.clone()))
                    .collect();

                let new_arrays = match (&chars1, &chars2, &chars3) {
                    (ByteAlphabet(c1), ByteAlphabet(c2), ByteAlphabet(c3)) => {
                        refine_segmentation_three_way(s1, s2, s3, c1, c2, c3, pieces, round == k - 1)
                    },
                    (Compact(c1), Compact(c2), Compact(c3)) => {
                        refine_segmentation_three_way(s1, s2, s3, c1, c2, c3, pieces, round == k - 1)
                    },
                    _ => panic!("Programmer messed up")
                };
                (s1, s2, s3, leader_bits) = new_arrays;

                if round != k - 1 {
                    if let (ByteAlphabet(ref mut c1), ByteAlphabet(ref mut c2), ByteAlphabet(ref mut c3),
                            Some(ByteAlphabet(ref mut t1)), Some(ByteAlphabet(ref mut t2)), Some(ByteAlphabet(ref mut t3)))
                        = (&mut chars1, &mut chars2, &mut chars3, &mut temp_char_buf_1, &mut temp_char_buf_2, &mut temp_char_buf_3)
                    {
                        index1.push_all_labels_forward(c1, t1, n_threads);
                        std::mem::swap(c1, t1);
                        index2.push_all_labels_forward(c2, t2, n_threads);
                        std::mem::swap(c2, t2);
                        aux.push_all_labels_forward(c3, t3, n_threads);
                        std::mem::swap(c3, t3);
                    } else if let (Compact(ref mut c1), Compact(ref mut c2), Compact(ref mut c3), None, None, None)
                        = (&mut chars1, &mut chars2, &mut chars3, &mut temp_char_buf_1, &mut temp_char_buf_2, &mut temp_char_buf_3)
                    {
                        let mut t1 = CompactIntVector::<3>::new(c1.len());
                        index1.push_all_labels_forward_compact(c1, &mut t1, n_threads);
                        std::mem::swap(c1, &mut t1);
                        drop(t1);
                        let mut t2 = CompactIntVector::<3>::new(c2.len());
                        index2.push_all_labels_forward_compact(c2, &mut t2, n_threads);
                        std::mem::swap(c2, &mut t2);
                        drop(t2);
                        let mut t3 = CompactIntVector::<3>::new(c3.len());
                        aux.push_all_labels_forward_compact(c3, &mut t3, n_threads);
                        std::mem::swap(c3, &mut t3);
                        drop(t3);
                    } else {
                        panic!("Programmer messed up");
                    }
                }
            }

            drop(temp_char_buf_1);
            drop(temp_char_buf_2);
            drop(temp_char_buf_3);

            let leader_bits = leader_bits.unwrap();

            ThreeWayInterleaving { s1, s2, s3, is_leader: leader_bits }
        })
    }
}

fn mark_dummy_nodes<V: ReadOnlyIntVector + Send + Sync>(s1: &BitVec, s2: &BitVec, chars1: &V, chars2: &V, n_threads: usize) -> BitVec {
    assert_eq!(s1.len(), s2.len());
    let merged_len = s1.len();

    let piece_len = merged_len.div_ceil(n_threads);
    let merged_piece_ranges: Vec<Range<usize>> = (0..n_threads).map(|t| t*piece_len..min((t+1)*piece_len, merged_len)).collect();
    let s1_piece_popcounts: Vec<usize> = merged_piece_ranges.par_iter().map(|range| s1[range.clone()].count_ones()).collect();
    let s2_piece_popcounts: Vec<usize> = merged_piece_ranges.par_iter().map(|range| s2[range.clone()].count_ones()).collect();
    let is_dummy_pieces = (0..n_threads).into_par_iter().map(|thread_idx| {
        let colex_range = &merged_piece_ranges[thread_idx];
        let mut is_dummy_bits = bitvec![u64, Lsb0 ;];
        is_dummy_bits.resize(colex_range.len(), false);

        let mut c1_idx: usize = s1_piece_popcounts[..thread_idx].iter().sum(); // Skip over previous pieces
        let mut c2_idx: usize = s2_piece_popcounts[..thread_idx].iter().sum(); // Skip over previous pieces
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

// Three-way version of refine_piece: merges three sorted character ranges into sub-groups.
#[allow(clippy::too_many_arguments)]
fn refine_piece_three_way<V: ReadOnlyIntVector + Send + Sync>(
    s1: &BitVec, s2: &BitVec, s3: &BitVec,
    chars1: &V, chars2: &V, chars3: &V,
    mut c1_i: usize, mut c2_i: usize, mut c3_i: usize,
    s1_range: Range<usize>, s2_range: Range<usize>, s3_range: Range<usize>,
    last_round: bool,
) -> (BitVec, BitVec, BitVec, Option<BitVec>) {
    let est_cap = s1_range.len() + s2_range.len() + s3_range.len();
    let mut out1 = BitVec::with_capacity(est_cap);
    let mut out2 = BitVec::with_capacity(est_cap);
    let mut out3 = BitVec::with_capacity(est_cap);
    let mut leader_bits = if last_round { BitVec::with_capacity(est_cap) } else { BitVec::new() };

    let mut s1_i = s1_range.start;
    let mut s2_i = s2_range.start;
    let mut s3_i = s3_range.start;

    while s1_i < s1_range.end {
        assert!(s2_i < s2_range.end);
        assert!(s3_i < s3_range.end);

        let len1 = leading_zeros(&s1[s1_i..]);
        let len2 = leading_zeros(&s2[s2_i..]);
        let len3 = leading_zeros(&s3[s3_i..]);

        let c1_end = c1_i + len1;
        let c2_end = c2_i + len2;
        let c3_end = c3_i + len3;

        let mut is_leader = true;
        while c1_i < c1_end || c2_i < c2_end || c3_i < c3_end {
            let cv1 = if c1_i == c1_end { u8::MAX } else { chars1.get(c1_i) as u8 };
            let cv2 = if c2_i == c2_end { u8::MAX } else { chars2.get(c2_i) as u8 };
            let cv3 = if c3_i == c3_end { u8::MAX } else { chars3.get(c3_i) as u8 };
            let c = min(cv1, min(cv2, cv3));

            let r1 = run_length_in_sorted_seq(chars1, c1_i..c1_end, c);
            let r2 = run_length_in_sorted_seq(chars2, c2_i..c2_end, c);
            let r3 = run_length_in_sorted_seq(chars3, c3_i..c3_end, c);

            if last_round {
                out1.push(r1 > 0);
                out2.push(r2 > 0);
                out3.push(r3 > 0);
            } else {
                zero_extend(&mut out1, r1); out1.push(true);
                zero_extend(&mut out2, r2); out2.push(true);
                zero_extend(&mut out3, r3); out3.push(true);
            }

            c1_i += r1;
            c2_i += r2;
            c3_i += r3;

            if last_round { leader_bits.push(is_leader); }
            is_leader = false;
        }

        s1_i += len1 + 1;
        s2_i += len2 + 1;
        s3_i += len3 + 1;
    }

    assert_eq!(s1_i, s1_range.end);
    assert_eq!(s2_i, s2_range.end);
    assert_eq!(s3_i, s3_range.end);

    (out1, out2, out3, if last_round { Some(leader_bits) } else { None })
}

// Three-way version of refine_segmentation.
// Input pieces are 6-tuples: (start_in_chars1, start_in_chars2, start_in_chars3,
//                              range_in_s1, range_in_s2, range_in_s3).
fn refine_segmentation_three_way<V: ReadOnlyIntVector + Send + Sync>(
    s1: BitVec, s2: BitVec, s3: BitVec,
    chars1: &V, chars2: &V, chars3: &V,
    input_pieces: Vec<(usize, usize, usize, Range<usize>, Range<usize>, Range<usize>)>,
    last_round: bool,
) -> (BitVec, BitVec, BitVec, Option<BitVec>) {
    let output_pieces: Vec<(BitVec, BitVec, BitVec, Option<BitVec>)> =
        input_pieces.par_iter().map(|(c1, c2, c3, r1, r2, r3)| {
            refine_piece_three_way(
                &s1, &s2, &s3, chars1, chars2, chars3,
                *c1, *c2, *c3, r1.clone(), r2.clone(), r3.clone(), last_round,
            )
        }).collect();

    drop(s1);
    drop(s2);
    drop(s3);

    let mut ps1 = vec![];
    let mut ps2 = vec![];
    let mut ps3 = vec![];
    let mut lp = if last_round { Some(vec![]) } else { None };
    for (a, b, c, d) in output_pieces {
        ps1.push(a);
        ps2.push(b);
        ps3.push(c);
        if last_round { lp.as_mut().unwrap().push(d.unwrap()); }
    }

    (
        crate::util::parallel_bitvec_concat(ps1),
        crate::util::parallel_bitvec_concat(ps2),
        crate::util::parallel_bitvec_concat(ps3),
        lp.map(crate::util::parallel_bitvec_concat),
    )
}

// Computes per-piece ranges aligned to group leaders and per-piece popcounts for s1, s2, s3
// in a ThreeWayInterleaving.
pub(super) fn compute_piece_ranges_three_way(
    merged_length: usize,
    n_threads: usize,
    interleaving: &ThreeWayInterleaving,
) -> (Vec<Range<usize>>, Vec<usize>, Vec<usize>, Vec<usize>) {
    let piece_len = merged_length.div_ceil(n_threads);
    let mut ranges: Vec<Range<usize>> =
        (0..n_threads).map(|t| t * piece_len..min((t + 1) * piece_len, merged_length)).collect();
    for piece_idx in 1..ranges.len() {
        let pair = &mut ranges[piece_idx - 1..=piece_idx];
        while !pair[1].is_empty() && !interleaving.is_leader[pair[1].start] {
            pair[1].start += 1;
            pair[0].end += 1;
        }
    }
    let s1_pc: Vec<usize> = ranges.par_iter().map(|r| interleaving.s1[r.clone()].count_ones()).collect();
    let s2_pc: Vec<usize> = ranges.par_iter().map(|r| interleaving.s2[r.clone()].count_ones()).collect();
    let s3_pc: Vec<usize> = ranges.par_iter().map(|r| interleaving.s3[r.clone()].count_ones()).collect();
    (ranges, s1_pc, s2_pc, s3_pc)
}

// Count result positions per piece. When `difference` is false: `(s1 && s2) || s3`
// (intersection with auxiliary). When `difference` is true: `(s1 && !s2) || s3`
// (set-difference with auxiliary).
pub(super) fn count_result_nodes_per_piece(
    merged_piece_ranges: &[Range<usize>],
    interleaving: &ThreeWayInterleaving,
    difference: bool,
) -> Vec<usize> {
    merged_piece_ranges.par_iter().map(|colex_range| {
        let mut count = 0usize;
        for merged_colex in colex_range.clone() {
            let is_result = if difference {
                (interleaving.s1[merged_colex] && !interleaving.s2[merged_colex])
                    || interleaving.s3[merged_colex]
            } else {
                (interleaving.s1[merged_colex] && interleaving.s2[merged_colex])
                    || interleaving.s3[merged_colex]
            };
            if is_result { count += 1; }
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


// Aligns piece boundaries to group leaders and computes per-piece popcounts for s1 and s2.
pub(super) fn compute_piece_ranges(merged_length: usize, n_threads: usize, interleaving: &MergeInterleaving) -> (Vec<Range<usize>>, Vec<usize>, Vec<usize>) {
    let piece_len = merged_length.div_ceil(n_threads);
    let mut merged_piece_ranges: Vec<Range<usize>> = (0..n_threads).map(|t| t*piece_len..min((t+1)*piece_len, merged_length)).collect();

    // Adjust the ranges so that they start with a leader
    for piece_idx in 1..merged_piece_ranges.len() {
        let pair = &mut merged_piece_ranges[piece_idx-1..=piece_idx]; // Borrow a pair of elements
        while !pair[1].is_empty() && !interleaving.is_leader[pair[1].start] {
            pair[1].start += 1;
            pair[0].end += 1;
        }
    }

    let s1_piece_popcounts: Vec<usize> = merged_piece_ranges.par_iter().map(|range| interleaving.s1[range.clone()].count_ones()).collect();
    let s2_piece_popcounts: Vec<usize> = merged_piece_ranges.par_iter().map(|range| interleaving.s2[range.clone()].count_ones()).collect();
    (merged_piece_ranges, s1_piece_popcounts, s2_piece_popcounts)
}


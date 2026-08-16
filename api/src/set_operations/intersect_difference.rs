use std::cmp::min;
use std::ops::Range;
use std::sync::Arc;
use bitvec::prelude::*;
use rayon::iter::IntoParallelIterator;
use rayon::iter::IntoParallelRefIterator;
use rayon::iter::IndexedParallelIterator;
use rayon::iter::ParallelIterator;
use crate::builder::BitPackedKmerSortingMem;
use crate::subsetseq::*;
use crate::sbwt::*;
use crate::atomic_bitmap::AtomicBitmap;
use super::interleaving::{MergeInterleaving, ThreeWayInterleaving,
    compute_piece_ranges, compute_piece_ranges_three_way, count_result_nodes_per_piece};
use super::common::{allocate_rows, transpose_and_concat_pieces, build_index,
    word_and_popcount_range, word_diff_popcount_range};

type BitVec = bitvec::vec::BitVec<u64, Lsb0>;

fn convert_index<SS: SubsetSeq>(idx: SbwtIndex<SubsetMatrix>, n_threads: usize, precalc_length: usize) -> SbwtIndex<SS> {
    let n_sets = idx.n_sets();
    let n_kmers = idx.n_kmers();
    let k = idx.k();
    let sigma = crate::util::DNA_ALPHABET.len();

    let pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();
    let new_rows: Vec<BitVec> = pool.install(|| {
        (0..sigma).into_par_iter().map(|c| {
            let mut row = BitVec::new();
            row.resize(n_sets, false);
            for i in 0..n_sets {
                if idx.sbwt().set_contains(i, c as u8) {
                    row.set(i, true);
                }
            }
            row
        }).collect()
    });

    build_index::<SS>(new_rows, n_kmers, k, n_threads, precalc_length)
}

/// Counts the number of result positions and dummy result positions using
/// parallel word-level operations over the interleaving's raw bit slices.
/// When `difference` is false, the predicate is `s1 & s2` (intersection);
/// when true, it is `s1 & !s2` (set difference).
fn compute_result_counts(interleaving: &MergeInterleaving, difference: bool) -> (usize, usize) {
    let s1w = interleaving.s1.as_raw_slice();
    let s2w = interleaving.s2.as_raw_slice();
    let dw  = interleaving.is_dummy.as_raw_slice();
    s1w.par_iter()
        .zip(s2w.par_iter())
        .zip(dw.par_iter())
        .map(|((a, b), d)| {
            let result: u64 = if difference { a & !b } else { a & b };
            (result.count_ones() as usize, (result & d).count_ones() as usize)
        })
        .reduce(|| (0, 0), |(i0, d0), (i1, d1)| (i0 + i1, d0 + d1))
}

/// Shared implementation for pass 1 of both intersection and set-difference.
///
/// For each (k-1)-suffix group that contributes a result node, marks every LF-step
/// target in `index1` colex space as having an incoming edge via a lock-free
/// [AtomicBitmap].  Each Rayon thread processes its own `piece_ranges` slice independently.
///
/// When `difference` is `false` (intersection), a group contributes iff it contains a
/// position where both `s1` and `s2` are set, and edge-bit `c` is kept when
/// `s1_or[c] && s2_or[c]`.  When `difference` is `true` (set difference), the predicate
/// becomes `s1 && !s2` and the edge condition becomes `s1_or[c] && !s2_or[c]`.
///
/// The intersection variant additionally pre-sets bit 0 (the shared root) so that the
/// root is always counted as covered, matching how `isec_length` is computed.
///
/// Also returns **dead-end dummies**: dummy result nodes with no surviving outgoing
/// edge, detected as a byproduct of the same `flush_group` scan (no extra pass, no
/// atomics needed since a group's successor status is thread-local). A (k-1)-suffix
/// group can have more than one `in_result` member (e.g. at a shared branch point), but
/// by the SBWT pruning invariant only the first one in colex order can retain a
/// non-empty row — every later member is unconditionally all-empty, hence a dead end
/// whenever it's a dummy. So each group buffers its `in_result` members (≤ `sigma`) and
/// resolves them at flush time.
fn pass1_has_incoming_impl<SS: SubsetSeq + Send + Sync>(
    index1: &SbwtIndex<SS>,
    index2: &SbwtIndex<SS>,
    interleaving: &MergeInterleaving,
    result_length: usize,
    sigma: usize,
    piece_ranges: &[Range<usize>],
    s1_pops: &[usize],
    s2_pops: &[usize],
    n_threads: usize,
    thread_pool: &rayon::ThreadPool,
    difference: bool,
) -> (AtomicBitmap, Vec<usize>) {
    let has_incoming = AtomicBitmap::new(index1.n_sets());
    // Intersection: pre-set bit 0 (shared root); difference: omit to avoid inflating n_incoming.
    if !difference && result_length > 0 { has_incoming.set(0, true); }

    let dead_end_dummies_per_piece: Vec<Vec<usize>> = thread_pool.install(|| {
        (0..n_threads).into_par_iter().map(|thread_idx| {
            let colex_range = &piece_ranges[thread_idx];
            let mut s1_colex: usize = s1_pops[..thread_idx].iter().sum();
            let mut s2_colex: usize = s2_pops[..thread_idx].iter().sum();

            assert!(interleaving.is_leader[colex_range.start]);

            let mut s1_group_or = vec![false; sigma];
            let mut s2_group_or = vec![false; sigma];
            let mut s1_group_read = false;
            let mut s2_group_read = false;
            // s1_colex of the first s1 element in the current group (used for lf_step).
            let mut s1_first_in_group = 0usize;
            // (s1_colex, is_dummy) of each `in_result` member seen so far in this group,
            // in colex order (bounded by `sigma`, see doc comment above).
            let mut group_result_members: Vec<(usize, bool)> = Vec::with_capacity(sigma);
            let mut local_dead_ends: Vec<usize> = Vec::new();

            // Marks LF-step targets for each surviving edge in the current group.
            // Returns whether at least one outgoing edge survived.
            let flush_group = |s1_or: &[bool], s2_or: &[bool], s1_first: usize| -> bool {
                let mut any_edge = false;
                for c in 0..sigma {
                    let edge = if difference {
                        s1_or[c] && !s2_or[c]
                    } else {
                        s1_or[c] && s2_or[c]
                    };
                    if edge {
                        any_edge = true;
                        has_incoming.set(index1.lf_step(s1_first, c), true);
                    }
                }
                any_edge
            };

            // Resolves dead-end status for the group's buffered `in_result` members: only
            // the first can have a surviving edge, every later one is unconditionally dead.
            // Colex position 0 is the root `$^k`, which is never removable, so it is never
            // reported even when the intersection leaves it without a single outgoing edge.
            let flush_dead_ends = |members: &[(usize, bool)], any_edge: bool, dead_ends: &mut Vec<usize>| {
                for (idx, &(pos, is_dummy)) in members.iter().enumerate() {
                    let dead = if idx == 0 { is_dummy && !any_edge } else { is_dummy };
                    if dead && pos != 0 { dead_ends.push(pos); }
                }
            };

            for merged_colex in colex_range.clone() {
                if interleaving.is_leader[merged_colex] && merged_colex > colex_range.start {
                    // Flush: if any result position exists, mark LF-step targets and
                    // resolve dead-end status for every buffered member.
                    if !group_result_members.is_empty() {
                        let any_edge = flush_group(&s1_group_or, &s2_group_or, s1_first_in_group);
                        flush_dead_ends(&group_result_members, any_edge, &mut local_dead_ends);
                    }
                    s1_group_read = false;
                    s2_group_read = false;
                    s1_group_or.fill(false);
                    s2_group_or.fill(false);
                    group_result_members.clear();
                }

                if !s1_group_read && interleaving.s1[merged_colex] {
                    for c in 0..sigma { s1_group_or[c] = index1.sbwt.set_contains(s1_colex, c as u8); }
                    s1_first_in_group = s1_colex;
                    s1_group_read = true;
                }
                if !s2_group_read && interleaving.s2[merged_colex] {
                    for c in 0..sigma { s2_group_or[c] = index2.sbwt.set_contains(s2_colex, c as u8); }
                    s2_group_read = true;
                }
                let in_result = if difference {
                    interleaving.s1[merged_colex] && !interleaving.s2[merged_colex]
                } else {
                    interleaving.s1[merged_colex] && interleaving.s2[merged_colex]
                };
                // `in_result` implies s1[merged_colex], so s1_colex (not yet incremented
                // this iteration) is this member's identifying position.
                if in_result {
                    group_result_members.push((s1_colex, interleaving.is_dummy[merged_colex]));
                }

                s1_colex += interleaving.s1[merged_colex] as usize;
                s2_colex += interleaving.s2[merged_colex] as usize;
            }

            // Flush the last group in this piece.
            if !group_result_members.is_empty() {
                let any_edge = flush_group(&s1_group_or, &s2_group_or, s1_first_in_group);
                flush_dead_ends(&group_result_members, any_edge, &mut local_dead_ends);
            }

            local_dead_ends
        }).collect()
    });

    // Pieces were scanned in increasing colex-range order and s1_colex is monotonically
    // increasing within and across pieces, so this concatenation is already sorted.
    let dead_end_dummies: Vec<usize> = dead_end_dummies_per_piece.into_iter().flatten().collect();

    (has_incoming, dead_end_dummies)
}

/// Pass 1 for intersection: delegates to [`pass1_has_incoming_impl`] with `difference = false`.
/// Returns the incoming-edge bitmap together with the `s1_colex` positions of dead-end
/// dummies (dummy result nodes with no surviving outgoing edge); see
/// [`pass1_has_incoming_impl`] for details.
fn pass1_has_incoming<SS: SubsetSeq + Send + Sync>(
    index1: &SbwtIndex<SS>,
    index2: &SbwtIndex<SS>,
    interleaving: &MergeInterleaving,
    isec_length: usize,
    sigma: usize,
    piece_ranges: &[Range<usize>],
    s1_pops: &[usize],
    s2_pops: &[usize],
    n_threads: usize,
    thread_pool: &rayon::ThreadPool,
) -> (AtomicBitmap, Vec<usize>) {
    pass1_has_incoming_impl(index1, index2, interleaving, isec_length, sigma,
        piece_ranges, s1_pops, s2_pops, n_threads, thread_pool, false)
}

/// Collects s1-colex positions of result nodes that have no incoming edge and therefore
/// need auxiliary dummy chains.  When `difference` is false the result predicate is `s1 & s2`
/// (intersection); when true it is `s1 & !s2` (set difference).
/// Positions marked in `redundant_s1` (dead dummy-chain nodes) are excluded: they are
/// dropped from the result, so rebuilding their chains would just resurrect them.
/// The scan is split into `n_threads` pieces; each piece derives its starting s1_colex
/// via a prefix popcount on `interleaving.s1`, so all pieces run fully in parallel.
fn collect_result_source_nodes(
    interleaving: &MergeInterleaving,
    has_incoming: &AtomicBitmap,
    redundant_s1: Option<&BitVec>,
    merged_length: usize,
    expected_count: usize,
    difference: bool,
    n_threads: usize,
    thread_pool: &rayon::ThreadPool,
) -> Vec<usize> {
    let piece_len = merged_length.div_ceil(n_threads);
    let pieces: Vec<Range<usize>> = (0..n_threads)
        .map(|t| t * piece_len..min((t + 1) * piece_len, merged_length))
        .filter(|r| !r.is_empty())
        .collect();

    // Pre-compute each piece's starting s1_colex via a parallel prefix popcount.
    let s1_starts: Vec<usize> = thread_pool.install(|| {
        pieces.par_iter()
            .map(|r| interleaving.s1[..r.start].count_ones())
            .collect()
    });

    // Scan each piece independently, collecting matching s1_colex values.
    let mut piece_results: Vec<Vec<usize>> = thread_pool.install(|| {
        pieces.par_iter()
            .zip(s1_starts.par_iter())
            .map(|(range, &s1_start)| {
                let mut local = Vec::new();
                let mut s1_colex = s1_start;
                for merged_colex in range.clone() {
                    if interleaving.s1[merged_colex] {
                        let in_result = if difference { !interleaving.s2[merged_colex] }
                                        else          {  interleaving.s2[merged_colex] };
                        if in_result
                            && !has_incoming.get(s1_colex)
                            && redundant_s1.map_or(true, |red| !red[s1_colex])
                        {
                            local.push(s1_colex);
                        }
                        s1_colex += 1;
                    }
                }
                local
            })
            .collect()
    });

    // Flatten in piece order; s1_colex values are globally increasing across pieces.
    let mut result = Vec::with_capacity(expected_count);
    for v in piece_results.drain(..) { result.extend(v); }
    result
}

/// Reconstructs k-mers from s1 colex positions via inverse-LF walks.
/// Each k-mer is independent, so the work is distributed across Rayon threads.
/// `index1` must already have select built before this function is called.
fn reconstruct_source_kmers<SS: SubsetSeq + Send + Sync>(
    index1: &SbwtIndex<SS>,
    source_colexes: &[usize],
    k: usize,
) -> Vec<Vec<u8>> {
    source_colexes.par_iter().map(|&s1_colex| {
        let mut buf = Vec::with_capacity(k);
        let mut pos = s1_colex;
        for _ in 0..k {
            match index1.inlabel(pos) {
                Some(c) => {
                    buf.push(c);
                    let char_idx  = index1.char_idx(c);
                    let c_start   = index1.lf_step(0, char_idx);
                    let rank_in_c = pos - c_start;
                    pos = index1.sbwt.select(char_idx as u8, rank_in_c).unwrap();
                }
                None => { buf.push(b'$'); }
            }
        }
        buf.reverse();
        buf
    }).collect()
}

/// Marks the dead-end dummies and, unconditionally, every ancestor on their dummy chains,
/// walking up to the root via inverse-LF in `index1`. The predecessor of a dummy is always
/// itself a dummy, so the walk never leaves the dummy subgraph; chains that merge stop at
/// the first already-marked node. Colex position 0 (the root) is never marked.
///
/// This deliberately over-kills: a marked ancestor may still feed a live sibling subtree.
/// Any such orphaned child is recognizable purely from `has_incoming`: a shared dummy's
/// child receives its unique incoming edge from that (now marked) parent, so an unmarked
/// child with its bit set has just lost its incoming edge. The sweep clears those bits, so
/// the orphans are collected as source nodes and the auxiliary SBWT rebuilds their padding
/// chains — re-creating exactly the marked ancestors they still need.
///
/// Returns the marked set in `index1`'s colex space.
fn mark_dead_chains<SS: SubsetSeq>(
    index1: &SbwtIndex<SS>,
    dead_ends: &[usize],
    has_incoming: &AtomicBitmap,
    sigma: usize,
) -> BitVec {
    let mut redundant = BitVec::new();
    redundant.resize(index1.n_sets(), false);

    let mut marked: Vec<usize> = Vec::new();
    for &d in dead_ends {
        let mut q = d;
        while !redundant[q] {
            redundant.set(q, true);
            marked.push(q);
            match index1.inverse_lf_step(q) {
                Some(p) if p != 0 => q = p,
                _ => break,
            }
        }
    }

    // Orphan sweep (after all marking, so `redundant` is final).
    let mut n_orphans = 0usize;
    for &p in &marked {
        for c in 0..sigma {
            if index1.sbwt.set_contains(p, c as u8) {
                let t = index1.lf_step(p, c);
                if !redundant[t] && has_incoming.get(t) {
                    has_incoming.set(t, false);
                    n_orphans += 1;
                }
            }
        }
    }

    log::info!("Dead-chain removal: {} node(s) marked from {} dead end(s); {} orphan(s) re-routed to the auxiliary index",
        marked.len(), dead_ends.len(), n_orphans);
    redundant
}

/// Direct pass 2: builds the intersection SBWT bit-rows when every intersection
/// position already has an incoming edge and no dummy chain repair is required.
/// For each group the AND of the two per-index ORs is written at the group's isec leader.
///
fn intersect_rows_direct<SS: SubsetSeq + Send + Sync>(
    index1: Arc<SbwtIndex<SS>>,
    index2: Arc<SbwtIndex<SS>>,
    interleaving: Arc<MergeInterleaving>,
    sigma: usize,
    piece_ranges: Vec<Range<usize>>,
    s1_pops: Vec<usize>,
    s2_pops: Vec<usize>,
    n_threads: usize,
    thread_pool: &rayon::ThreadPool,
) -> Vec<BitVec> {
    thread_pool.install(|| {
        let isec_piece_pops: Vec<usize> = piece_ranges.par_iter().map(|range| {
            word_and_popcount_range(
                interleaving.s1.as_raw_slice(),
                interleaving.s2.as_raw_slice(),
                range.clone(),
            )
        }).collect();

        let pieces_vecvec: Vec<Vec<BitVec>> = (0..n_threads).into_par_iter().map(|thread_idx| {
            let colex_range  = &piece_ranges[thread_idx];
            let mut new_rows = allocate_rows(sigma, isec_piece_pops[thread_idx]);

            let mut s1_colex: usize = s1_pops[..thread_idx].iter().sum();
            let mut s2_colex: usize = s2_pops[..thread_idx].iter().sum();

            assert!(interleaving.is_leader[colex_range.start]);

            // Per-group accumulators: OR each index's edge-bits, then AND at flush time.
            let mut s1_group_or = vec![false; sigma];
            let mut s2_group_or = vec![false; sigma];
            let mut s1_group_read = false;
            let mut s2_group_read = false;
            // Piece-relative index of the first shared position in the current group.
            let mut piece_rel_isec_leader: Option<usize> = None;
            let mut isec_colex_in_piece = 0usize;

            for merged_colex in colex_range.clone() {
                if interleaving.is_leader[merged_colex] && merged_colex > colex_range.start {
                    // Flush: write AND-of-ORs at the isec leader.
                    if let Some(l) = piece_rel_isec_leader {
                        for c in 0..sigma {
                            if s1_group_or[c] && s2_group_or[c] { new_rows[c].set(l, true); }
                        }
                    }
                    s1_group_read = false;
                    s2_group_read = false;
                    piece_rel_isec_leader = None;
                }

                if !s1_group_read && interleaving.s1[merged_colex] {
                    for c in 0..sigma { s1_group_or[c] = index1.sbwt.set_contains(s1_colex, c as u8); }
                    s1_group_read = true;
                }
                if !s2_group_read && interleaving.s2[merged_colex] {
                    for c in 0..sigma { s2_group_or[c] = index2.sbwt.set_contains(s2_colex, c as u8); }
                    s2_group_read = true;
                }
                if interleaving.s1[merged_colex] && interleaving.s2[merged_colex] {
                    if piece_rel_isec_leader.is_none() {
                        piece_rel_isec_leader = Some(isec_colex_in_piece);
                    }
                    isec_colex_in_piece += 1;
                }

                s1_colex += interleaving.s1[merged_colex] as usize;
                s2_colex += interleaving.s2[merged_colex] as usize;
            }

            // Flush the last group in this piece.
            if let Some(l) = piece_rel_isec_leader {
                for c in 0..sigma {
                    if s1_group_or[c] && s2_group_or[c] { new_rows[c].set(l, true); }
                }
            }

            assert_eq!(s1_colex, s1_pops[..=thread_idx].iter().sum());
            assert_eq!(s2_colex, s2_pops[..=thread_idx].iter().sum());
            assert_eq!(isec_colex_in_piece, isec_piece_pops[thread_idx]);

            new_rows
        }).collect();

        drop(index1);
        drop(index2);
        drop(interleaving);

        log::info!("[intersect] Transposing and concatenating pieces");
        transpose_and_concat_pieces(pieces_vecvec, sigma)
    })
}

/// Dummy-repair pass 2: some intersection k-mers are new source nodes (their predecessor k-mers
/// exist in each individual SBWT but not in the intersection), so their dummy predecessor chains
/// are absent. Remedy: reconstruct those source k-mers (in parallel), build a small auxiliary
/// SBWT from them, then use a single three-way interleaving of (index1, index2, aux) to produce
/// the complete bit-rows. Edge rule: c is set iff aux has it OR both index1 and index2 have it.
///
/// The same pass also leaves out the dead dummy chains: `redundant_s1` (index1 colex space,
/// precomputed by [`mark_dead_chains`]) marks every node of every chain whose real
/// k-mer did not survive the intersection. Marked nodes are not emitted as result nodes and
/// edges whose LF-step target is marked are not written. This is safe even when a marked
/// node lies on a new source's padding chain: the aux contains the very same node, so the
/// `s3` term of the predicate re-emits it and `s3_or` supplies both its row and its
/// incoming edge — dead is dead, and the aux brings back exactly what it needs.
fn intersect_rows_with_dummy_repair<SS: SubsetSeq + Send + Sync + Clone>(
    index1: SbwtIndex<SS>,
    index2: Arc<SbwtIndex<SS>>,
    source_colexes: Vec<usize>,
    redundant_s1: Option<&BitVec>,
    sigma: usize,
    k: usize,
    optimize_peak_ram: bool,
    n_threads: usize,
    thread_pool: &rayon::ThreadPool,
) -> Vec<BitVec> {
    log::info!("[intersect] Dummy repair: reconstructing {} source k-mer(s) via inverse-LF (parallel)", source_colexes.len());
    let source_kmers = thread_pool.install(|| reconstruct_source_kmers(&index1, &source_colexes, k));
    drop(source_colexes);
    let index1 = Arc::new(index1);

    log::info!("[intersect] Dummy repair: building auxiliary SBWT from source k-mers");
    let (aux_submatrix, _) = BitPackedKmerSortingMem::new_from_vecs(&source_kmers, k)
        .n_threads(n_threads)
        .add_all_dummy_paths(true)
        .run();
    drop(source_kmers);

    let aux_arc = Arc::new(convert_index::<SS>(aux_submatrix, n_threads, 0));

    // Build a single three-way interleaving of (index1, index2, aux).
    // This replaces three separate 2-way interleaving computations + two intermediate merges,
    // reducing dummy-repair overhead from O(3·(N+Naux)·k) to O((N+Naux)·k).
    log::info!("[intersect] Dummy repair: building three-way interleaving");
    let three_way = ThreeWayInterleaving::new(&index1, &index2, &aux_arc, optimize_peak_ram, n_threads);
    let three_way_len = three_way.s1.len();

    log::info!("[intersect] Dummy repair: building intersection SBWT bit-rows via three-way pass (parallel)");
    let (tw_ranges, s1_pc, s2_pc, s3_pc) =
        thread_pool.install(|| compute_piece_ranges_three_way(three_way_len, n_threads, &three_way));
    let result_piece_counts = thread_pool.install(||
        count_result_nodes_per_piece(&tw_ranges, &three_way, false, redundant_s1, &s1_pc));

    // For each (k-1)-suffix group: result node exists iff (s1 && s2 && !dead) || s3.
    // Edge-bit c is set iff s3_or[c] || (s1_or[c] && s2_or[c] && target not dead).
    thread_pool.install(|| {
        let pieces_vecvec: Vec<Vec<BitVec>> = (0..n_threads).into_par_iter().map(|thread_idx| {
            let colex_range = &tw_ranges[thread_idx];
            let mut piece_rows = allocate_rows(sigma, result_piece_counts[thread_idx]);

            let mut s1_colex: usize = s1_pc[..thread_idx].iter().sum();
            let mut s2_colex: usize = s2_pc[..thread_idx].iter().sum();
            let mut s3_colex: usize = s3_pc[..thread_idx].iter().sum();

            assert!(three_way.is_leader[colex_range.start]);

            // Per-group ORs; reset to false at each group boundary to prevent stale bleed-through
            // (stale s1/s2 could trigger spurious s3_or || (s1_or && s2_or) writes otherwise).
            let mut s1_or = vec![false; sigma];
            let mut s2_or = vec![false; sigma];
            let mut s3_or = vec![false; sigma];
            let mut s1_read = false;
            let mut s2_read = false;
            let mut s3_read = false;
            // s1_colex of the first s1 element in the current group (used for lf_step
            // when filtering edges into dead nodes).
            let mut s1_first_in_group = 0usize;
            let mut piece_rel_result_leader: Option<usize> = None;
            let mut result_colex_in_piece = 0usize;

            for merged_colex in colex_range.clone() {
                if three_way.is_leader[merged_colex] && merged_colex > colex_range.start {
                    // Flush: write edge bits at the result leader (if any). Edges into the
                    // dead set can only originate from the root's group or a dead
                    // (aux-revived) group — every dead node's parent is dead or the root —
                    // so all other groups skip the per-edge lf_step target check.
                    if let Some(l) = piece_rel_result_leader {
                        let check_dead_targets = redundant_s1
                            .filter(|red| s1_first_in_group == 0 || red[s1_first_in_group]);
                        for c in 0..sigma {
                            if s3_or[c] || (s1_or[c] && s2_or[c]
                                && check_dead_targets.map_or(true, |red| !red[index1.lf_step(s1_first_in_group, c)]))
                            { piece_rows[c].set(l, true); }
                        }
                    }
                    s1_read = false; s2_read = false; s3_read = false;
                    s1_or.fill(false); s2_or.fill(false); s3_or.fill(false);
                    piece_rel_result_leader = None;
                }

                if !s1_read && three_way.s1[merged_colex] {
                    for c in 0..sigma { s1_or[c] = index1.sbwt.set_contains(s1_colex, c as u8); }
                    s1_first_in_group = s1_colex;
                    s1_read = true;
                }
                if !s2_read && three_way.s2[merged_colex] {
                    for c in 0..sigma { s2_or[c] = index2.sbwt.set_contains(s2_colex, c as u8); }
                    s2_read = true;
                }
                if !s3_read && three_way.s3[merged_colex] {
                    for c in 0..sigma { s3_or[c] = aux_arc.sbwt.set_contains(s3_colex, c as u8); }
                    s3_read = true;
                }

                // s1&&s2 implies s1[merged_colex], so s1_colex (not yet incremented this
                // iteration) is this node's position in index1's colex space; dead nodes
                // are not result nodes unless the aux re-emits them via s3.
                let is_result = (three_way.s1[merged_colex] && three_way.s2[merged_colex]
                        && redundant_s1.map_or(true, |red| !red[s1_colex]))
                    || three_way.s3[merged_colex];
                if is_result {
                    if piece_rel_result_leader.is_none() {
                        piece_rel_result_leader = Some(result_colex_in_piece);
                    }
                    result_colex_in_piece += 1;
                }

                s1_colex += three_way.s1[merged_colex] as usize;
                s2_colex += three_way.s2[merged_colex] as usize;
                s3_colex += three_way.s3[merged_colex] as usize;
            }

            // Flush the last group (same dead-target gating as above).
            if let Some(l) = piece_rel_result_leader {
                let check_dead_targets = redundant_s1
                    .filter(|red| s1_first_in_group == 0 || red[s1_first_in_group]);
                for c in 0..sigma {
                    if s3_or[c] || (s1_or[c] && s2_or[c]
                        && check_dead_targets.map_or(true, |red| !red[index1.lf_step(s1_first_in_group, c)]))
                    { piece_rows[c].set(l, true); }
                }
            }

            assert_eq!(s1_colex, s1_pc[..=thread_idx].iter().sum());
            assert_eq!(s2_colex, s2_pc[..=thread_idx].iter().sum());
            assert_eq!(s3_colex, s3_pc[..=thread_idx].iter().sum());
            assert_eq!(result_colex_in_piece, result_piece_counts[thread_idx]);

            piece_rows
        }).collect();

        drop(index1);
        drop(index2);
        drop(aux_arc);

        log::info!("[intersect] Dummy repair: transposing and concatenating pieces");
        transpose_and_concat_pieces(pieces_vecvec, sigma)
    })
}

// ── public API ─────────────────────────────────────────────────────────────────

/// Computes the intersection of two [SbwtIndex] structures. The intersection k-mers are those
/// present in **both** `index1` and `index2`. This mirrors [merge] in structure but uses AND
/// logic for incoming edges instead of OR, and restricts output positions to those shared by
/// both SBWTs.
///
/// For each (k-1)-suffix group in the merged SBWT, incoming-edge bits are OR-ed separately for
/// each index across all positions in the group, and the two resulting ORs are then AND-ed. The
/// result is stored at the first shared position of the group (the group's leader in the
/// intersection SBWT). This differs from [merge] where every position writes OR-into-leader directly.
///
/// Only dummy nodes shared by both SBWTs are retained. After pass 1, a parallel count detects
/// whether any k-mers became source nodes in the intersection, and pass 1 also reports the
/// dead-end dummies (shared dummies whose chain no longer leads to a surviving k-mer). If there
/// is neither, pass 2 returns the result directly. Otherwise, an auxiliary SBWT is built from
/// those source k-mers (with parallel k-mer reconstruction) and merged via a three-way
/// interleaving that repairs the missing dummy chains and strips the dead-end ones.
pub fn intersect<SS: SubsetSeq + Send + Sync + Clone>(
    index1: Arc<SbwtIndex<SS>>,
    index2: Arc<SbwtIndex<SS>>,
    interleaving: Arc<MergeInterleaving>,
    new_prefix_lookup_table_length: usize,
    optimize_peak_ram: bool,
    n_threads: usize,
) -> SbwtIndex<SS> {
    let sigma = crate::util::DNA_ALPHABET.len();

    assert!(index1.k() == index2.k());
    let k = index1.k();
    let merged_length = interleaving.s1.len();

    let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();

    let (isec_length, n_dummies) = thread_pool.install(|| compute_result_counts(&interleaving, false));
    let n_kmers = isec_length - n_dummies;
    log::info!("Intersecting into {} distinct k-mers", n_kmers);
    log::info!("Number of sets in intersection SBWT: {}", isec_length);

    let (piece_ranges, s1_pops, s2_pops) =
        thread_pool.install(|| compute_piece_ranges(merged_length, n_threads, &interleaving));

    log::info!("[intersect] Pass 1: computing incoming-edge coverage (parallel)");
    let (has_incoming, dead_end_dummies) = pass1_has_incoming(
        &index1, &index2, &interleaving,
        isec_length, sigma, &piece_ranges, &s1_pops, &s2_pops,
        n_threads, &thread_pool,
    );
    log::info!("[intersect] Dead-end dummies (no surviving outgoing edge): {}", dead_end_dummies.len());

    // Count covered positions in parallel; bit 0 is pre-set (root), so complete == isec_length.
    let n_incoming: usize = thread_pool.install(|| {
        has_incoming.data.par_iter()
            .map(|w| w.load(std::sync::atomic::Ordering::Acquire).count_ones() as usize)
            .sum()
    });

    let new_rows = if dead_end_dummies.is_empty() && n_incoming == isec_length {
        log::info!("All {} intersection positions have incoming edges; no auxiliary dummy chains needed", isec_length);
        log::info!("[intersect] Pass 2: building intersection SBWT bit-rows (parallel)");
        intersect_rows_direct(
            index1, index2, interleaving,
            sigma, piece_ranges, s1_pops, s2_pops,
            n_threads, &thread_pool,
        )
    } else {
        log::info!("Intersection SBWT structurally incomplete ({} source(s) missing dummy chains); merging with auxiliary index",
            isec_length - n_incoming);
        // Unwrap the Arc (refcount is 1 here since nothing else cloned it) and build select
        // now: both the dead-chain walk and source k-mer reconstruction need it.
        let mut index1 = Arc::try_unwrap(index1).unwrap_or_else(|_| panic!("index1 Arc must be uniquely owned at this point"));
        index1.sbwt.build_select();

        // Mark the dead chains (unconditionally, up to the root) and clear the incoming-edge
        // bits of any orphaned children BEFORE collecting sources: orphans must be collected
        // so the auxiliary index rebuilds their padding chains. The three-way pass excludes
        // marked nodes and the edges into them while building the rows.
        let redundant_s1: Option<BitVec> = if dead_end_dummies.is_empty() { None } else {
            Some(mark_dead_chains(&index1, &dead_end_dummies, &has_incoming, sigma))
        };
        drop(dead_end_dummies);

        let source_colexes = collect_result_source_nodes(
            &interleaving, &has_incoming, redundant_s1.as_ref(),
            merged_length, isec_length - n_incoming,
            false, n_threads, &thread_pool,
        );
        drop(interleaving); // No longer needed; free before the aux-index build.
        intersect_rows_with_dummy_repair(
            index1, index2, source_colexes, redundant_s1.as_ref(),
            sigma, k, optimize_peak_ram, n_threads, &thread_pool,
        )
    };

    log::info!("[intersect] Building subset rank structure");
    build_index(new_rows, n_kmers, k, n_threads, new_prefix_lookup_table_length)
}

// ── difference helpers ──────────────────────────────────────────────────────

/// Pass 1 for set difference: delegates to [`pass1_has_incoming_impl`] with `difference = true`.
fn pass1_has_incoming_diff<SS: SubsetSeq + Send + Sync>(
    index1: &SbwtIndex<SS>,
    index2: &SbwtIndex<SS>,
    interleaving: &MergeInterleaving,
    diff_length: usize,
    sigma: usize,
    piece_ranges: &[Range<usize>],
    s1_pops: &[usize],
    s2_pops: &[usize],
    n_threads: usize,
    thread_pool: &rayon::ThreadPool,
) -> (AtomicBitmap, Vec<usize>) {
    pass1_has_incoming_impl(index1, index2, interleaving, diff_length, sigma,
        piece_ranges, s1_pops, s2_pops, n_threads, thread_pool, true)
}

fn difference_rows_direct<SS: SubsetSeq + Send + Sync>(
    index1: Arc<SbwtIndex<SS>>,
    index2: Arc<SbwtIndex<SS>>,
    interleaving: Arc<MergeInterleaving>,
    sigma: usize,
    piece_ranges: Vec<Range<usize>>,
    s1_pops: Vec<usize>,
    s2_pops: Vec<usize>,
    n_threads: usize,
    thread_pool: &rayon::ThreadPool,
) -> Vec<BitVec> {
    thread_pool.install(|| {
        // Count difference positions per piece: s1 & !s2
        let diff_piece_pops: Vec<usize> = piece_ranges.par_iter().map(|range| {
            word_diff_popcount_range(
                interleaving.s1.as_raw_slice(),
                interleaving.s2.as_raw_slice(),
                range.clone(),
            )
        }).collect();

        let pieces_vecvec: Vec<Vec<BitVec>> = (0..n_threads).into_par_iter().map(|thread_idx| {
            let colex_range  = &piece_ranges[thread_idx];
            let mut new_rows = allocate_rows(sigma, diff_piece_pops[thread_idx]);

            let mut s1_colex: usize = s1_pops[..thread_idx].iter().sum();
            let mut s2_colex: usize = s2_pops[..thread_idx].iter().sum();

            assert!(interleaving.is_leader[colex_range.start]);

            let mut s1_group_or = vec![false; sigma];
            let mut s2_group_or = vec![false; sigma];
            let mut s1_group_read = false;
            let mut s2_group_read = false;
            let mut piece_rel_diff_leader: Option<usize> = None;
            let mut diff_colex_in_piece = 0usize;

            for merged_colex in colex_range.clone() {
                if interleaving.is_leader[merged_colex] && merged_colex > colex_range.start {
                    if let Some(l) = piece_rel_diff_leader {
                        for c in 0..sigma {
                            if s1_group_or[c] && !s2_group_or[c] { new_rows[c].set(l, true); }
                        }
                    }
                    s1_group_read = false;
                    s2_group_read = false;
                    s1_group_or.fill(false);
                    s2_group_or.fill(false);
                    piece_rel_diff_leader = None;
                }

                if !s1_group_read && interleaving.s1[merged_colex] {
                    for c in 0..sigma { s1_group_or[c] = index1.sbwt.set_contains(s1_colex, c as u8); }
                    s1_group_read = true;
                }
                if !s2_group_read && interleaving.s2[merged_colex] {
                    for c in 0..sigma { s2_group_or[c] = index2.sbwt.set_contains(s2_colex, c as u8); }
                    s2_group_read = true;
                }
                if interleaving.s1[merged_colex] && !interleaving.s2[merged_colex] {
                    if piece_rel_diff_leader.is_none() {
                        piece_rel_diff_leader = Some(diff_colex_in_piece);
                    }
                    diff_colex_in_piece += 1;
                }

                s1_colex += interleaving.s1[merged_colex] as usize;
                s2_colex += interleaving.s2[merged_colex] as usize;
            }

            // Flush the last group in this piece.
            if let Some(l) = piece_rel_diff_leader {
                for c in 0..sigma {
                    if s1_group_or[c] && !s2_group_or[c] { new_rows[c].set(l, true); }
                }
            }

            assert_eq!(s1_colex, s1_pops[..=thread_idx].iter().sum());
            assert_eq!(s2_colex, s2_pops[..=thread_idx].iter().sum());
            assert_eq!(diff_colex_in_piece, diff_piece_pops[thread_idx]);

            new_rows
        }).collect();

        drop(interleaving);

        log::info!("[difference] Transposing and concatenating pieces");
        let mut rows = transpose_and_concat_pieces(pieces_vecvec, sigma);

        // The SBWT invariant requires position 0 to be the root node ($$$…$).  The root is
        // always shared (s1 = s2 = 1) so it never appears among the s1 & !s2 diff positions.
        // We prepend it here with the appropriate edge bits: edge c is set iff $$$$c is
        // exclusive to index1 (i.e. the immediately following dummy is a diff node).
        for c in 0..sigma {
            let root_edge = index1.sbwt.set_contains(0, c as u8)
                && !index2.sbwt.set_contains(0, c as u8);
            let mut new_row = std::mem::take(&mut rows[c]);
            new_row.insert(0, root_edge);
            rows[c] = new_row;
        }

        drop(index1);
        drop(index2);

        rows
    })
}

fn difference_rows_with_dummy_repair<SS: SubsetSeq + Send + Sync + Clone>(
    index1: SbwtIndex<SS>,
    index2: Arc<SbwtIndex<SS>>,
    source_colexes: Vec<usize>,
    redundant_s1: Option<&BitVec>,
    sigma: usize,
    k: usize,
    optimize_peak_ram: bool,
    n_threads: usize,
    thread_pool: &rayon::ThreadPool,
) -> Vec<BitVec> {
    log::info!("[difference] Dummy repair: reconstructing {} source k-mer(s) via inverse-LF (parallel)", source_colexes.len());
    let source_kmers = thread_pool.install(|| reconstruct_source_kmers(&index1, &source_colexes, k));
    drop(source_colexes);
    let index1 = Arc::new(index1);

    log::info!("[difference] Dummy repair: building auxiliary SBWT from source k-mers");
    let (aux_submatrix, _) = BitPackedKmerSortingMem::new_from_vecs(&source_kmers, k)
        .n_threads(n_threads)
        .add_all_dummy_paths(true)
        .run();
    drop(source_kmers);

    let aux_arc = Arc::new(convert_index::<SS>(aux_submatrix, n_threads, 0));

    log::info!("[difference] Dummy repair: building three-way interleaving");
    // Reuse the same three-way interleaving structure; only the result predicate changes below.
    let three_way = ThreeWayInterleaving::new(&index1, &index2, &aux_arc, optimize_peak_ram, n_threads);
    let three_way_len = three_way.s1.len();

    log::info!("[difference] Dummy repair: building difference SBWT bit-rows via three-way pass (parallel)");
    let (tw_ranges, s1_pc, s2_pc, s3_pc) =
        thread_pool.install(|| compute_piece_ranges_three_way(three_way_len, n_threads, &three_way));
    let result_piece_counts = thread_pool.install(||
        count_result_nodes_per_piece(&tw_ranges, &three_way, true, redundant_s1, &s1_pc));

    // For each (k-1)-suffix group: result node exists iff (s1 && !s2 && !dead) || s3.
    // Edge-bit c is set iff s3_or[c] || (s1_or[c] && !s2_or[c] && target not dead).
    thread_pool.install(|| {
        let pieces_vecvec: Vec<Vec<BitVec>> = (0..n_threads).into_par_iter().map(|thread_idx| {
            let colex_range = &tw_ranges[thread_idx];
            let mut piece_rows = allocate_rows(sigma, result_piece_counts[thread_idx]);

            let mut s1_colex: usize = s1_pc[..thread_idx].iter().sum();
            let mut s2_colex: usize = s2_pc[..thread_idx].iter().sum();
            let mut s3_colex: usize = s3_pc[..thread_idx].iter().sum();

            assert!(three_way.is_leader[colex_range.start]);

            let mut s1_or = vec![false; sigma];
            let mut s2_or = vec![false; sigma];
            let mut s3_or = vec![false; sigma];
            let mut s1_read = false;
            let mut s2_read = false;
            let mut s3_read = false;
            // s1_colex of the first s1 element in the current group (used for lf_step
            // when filtering edges into dead nodes).
            let mut s1_first_in_group = 0usize;
            let mut piece_rel_result_leader: Option<usize> = None;
            let mut result_colex_in_piece = 0usize;

            for merged_colex in colex_range.clone() {
                if three_way.is_leader[merged_colex] && merged_colex > colex_range.start {
                    // Flush: write edge bits at the result leader (if any); see the
                    // intersection variant for the dead-target gating rationale.
                    if let Some(l) = piece_rel_result_leader {
                        let check_dead_targets = redundant_s1
                            .filter(|red| s1_first_in_group == 0 || red[s1_first_in_group]);
                        for c in 0..sigma {
                            if s3_or[c] || (s1_or[c] && !s2_or[c]
                                && check_dead_targets.map_or(true, |red| !red[index1.lf_step(s1_first_in_group, c)]))
                            { piece_rows[c].set(l, true); }
                        }
                    }
                    s1_read = false; s2_read = false; s3_read = false;
                    s1_or.fill(false); s2_or.fill(false); s3_or.fill(false);
                    piece_rel_result_leader = None;
                }

                if !s1_read && three_way.s1[merged_colex] {
                    for c in 0..sigma { s1_or[c] = index1.sbwt.set_contains(s1_colex, c as u8); }
                    s1_first_in_group = s1_colex;
                    s1_read = true;
                }
                if !s2_read && three_way.s2[merged_colex] {
                    for c in 0..sigma { s2_or[c] = index2.sbwt.set_contains(s2_colex, c as u8); }
                    s2_read = true;
                }
                if !s3_read && three_way.s3[merged_colex] {
                    for c in 0..sigma { s3_or[c] = aux_arc.sbwt.set_contains(s3_colex, c as u8); }
                    s3_read = true;
                }

                // Dead nodes are not result nodes unless the aux re-emits them via s3;
                // s1_colex (not yet incremented) is this node's index1 colex position.
                let is_result = (three_way.s1[merged_colex] && !three_way.s2[merged_colex]
                        && redundant_s1.map_or(true, |red| !red[s1_colex]))
                    || three_way.s3[merged_colex];
                if is_result {
                    if piece_rel_result_leader.is_none() {
                        piece_rel_result_leader = Some(result_colex_in_piece);
                    }
                    result_colex_in_piece += 1;
                }

                s1_colex += three_way.s1[merged_colex] as usize;
                s2_colex += three_way.s2[merged_colex] as usize;
                s3_colex += three_way.s3[merged_colex] as usize;
            }

            // Flush the last group (same dead-target gating as above).
            if let Some(l) = piece_rel_result_leader {
                let check_dead_targets = redundant_s1
                    .filter(|red| s1_first_in_group == 0 || red[s1_first_in_group]);
                for c in 0..sigma {
                    if s3_or[c] || (s1_or[c] && !s2_or[c]
                        && check_dead_targets.map_or(true, |red| !red[index1.lf_step(s1_first_in_group, c)]))
                    { piece_rows[c].set(l, true); }
                }
            }

            assert_eq!(s1_colex, s1_pc[..=thread_idx].iter().sum());
            assert_eq!(s2_colex, s2_pc[..=thread_idx].iter().sum());
            assert_eq!(s3_colex, s3_pc[..=thread_idx].iter().sum());
            assert_eq!(result_colex_in_piece, result_piece_counts[thread_idx]);

            piece_rows
        }).collect();

        drop(index1);
        drop(index2);
        drop(aux_arc);

        log::info!("[difference] Dummy repair: transposing and concatenating pieces");
        transpose_and_concat_pieces(pieces_vecvec, sigma)
    })
}

/// Computes the set difference `index1 \ index2`: the [SbwtIndex] containing exactly the k-mers
/// that are present in `index1` but **not** in `index2`.
///
/// The algorithm mirrors [`intersect`] but uses the complement predicate: a position in the
/// merged SBWT is a result node iff `s1[i] && !s2[i]`. Incoming-edge bits follow the same
/// group-OR logic but the AND is replaced by `s1_or[c] && !s2_or[c]`: an outgoing character `c`
/// is kept iff the successor k-mer (first character `c` followed by the current (k-1)-suffix)
/// is in `index1` but not in `index2`.
///
/// Like [`intersect`], a direct pass is taken when no difference k-mer becomes a source node
/// and no dummy became a dead end. Otherwise, an auxiliary SBWT is built from the reconstructed
/// source k-mers and a three-way interleaving with edge rule `s3_or[c] || (s1_or[c] && !s2_or[c])`
/// is used to repair the missing dummy chains, while dead dummy chains (marked by
/// [`mark_dead_chains`]) are excluded from the result exactly as in the intersection.
pub fn difference<SS: SubsetSeq + Send + Sync + Clone>(
    index1: Arc<SbwtIndex<SS>>,
    index2: Arc<SbwtIndex<SS>>,
    interleaving: Arc<MergeInterleaving>,
    new_prefix_lookup_table_length: usize,
    optimize_peak_ram: bool,
    n_threads: usize,
) -> SbwtIndex<SS> {
    let sigma = crate::util::DNA_ALPHABET.len();

    assert!(index1.k() == index2.k());
    let k = index1.k();
    let merged_length = interleaving.s1.len();

    let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();

    let (diff_length, n_dummies) = thread_pool.install(|| compute_result_counts(&interleaving, true));
    let n_kmers = diff_length - n_dummies;
    log::info!("Differencing into {} distinct k-mers", n_kmers);
    log::info!("Number of sets in difference SBWT: {}", diff_length);

    let (piece_ranges, s1_pops, s2_pops) =
        thread_pool.install(|| compute_piece_ranges(merged_length, n_threads, &interleaving));

    log::info!("[difference] Pass 1: computing incoming-edge coverage (parallel)");
    let (has_incoming, dead_end_dummies) = pass1_has_incoming_diff(
        &index1, &index2, &interleaving,
        diff_length, sigma, &piece_ranges, &s1_pops, &s2_pops,
        n_threads, &thread_pool,
    );
    log::info!("[difference] Dead-end dummies (no surviving outgoing edge): {}", dead_end_dummies.len());

    // Count covered positions in parallel. Every covered position is a difference result
    // node (a surviving edge's target is never in index2), so full coverage
    // (n_incoming == diff_length) means no result k-mer became a source node.
    let n_incoming: usize = thread_pool.install(|| {
        has_incoming.data.par_iter()
            .map(|w| w.load(std::sync::atomic::Ordering::Acquire).count_ones() as usize)
            .sum()
    });

    let new_rows = if dead_end_dummies.is_empty() && n_incoming == diff_length {
        log::info!("[difference] All real difference k-mers have incoming edges; no auxiliary dummy chains needed");
        log::info!("[difference] Pass 2: building difference SBWT bit-rows (parallel)");
        difference_rows_direct(
            index1, index2, interleaving,
            sigma, piece_ranges, s1_pops, s2_pops,
            n_threads, &thread_pool,
        )
    } else {
        log::info!("[difference] {} source(s) missing dummy chains; merging with auxiliary index",
            diff_length - n_incoming);
        // Unwrap the Arc (refcount is 1 here since nothing else cloned it) and build select
        // now: both the dead-chain walk and source k-mer reconstruction need it.
        let mut index1 = Arc::try_unwrap(index1).unwrap_or_else(|_| panic!("index1 Arc must be uniquely owned at this point"));
        index1.sbwt.build_select();

        // Mark the dead chains and clear the incoming-edge bits of any orphaned children
        // BEFORE collecting sources, exactly as in [intersect]. The empty-source case is
        // sound: an aux built from zero k-mers is a root-only SBWT, and the difference
        // result receives its root via the s3 term (the root is never s1 && !s2).
        let redundant_s1: Option<BitVec> = if dead_end_dummies.is_empty() { None } else {
            Some(mark_dead_chains(&index1, &dead_end_dummies, &has_incoming, sigma))
        };
        drop(dead_end_dummies);

        let source_colexes = collect_result_source_nodes(
            &interleaving, &has_incoming, redundant_s1.as_ref(),
            merged_length,
            n_kmers, // upper-bound capacity hint
            true, n_threads, &thread_pool,
        );
        drop(interleaving);
        difference_rows_with_dummy_repair(
            index1, index2, source_colexes, redundant_s1.as_ref(),
            sigma, k, optimize_peak_ram, n_threads, &thread_pool,
        )
    };

    log::info!("[difference] Building subset rank structure");
    build_index(new_rows, n_kmers, k, n_threads, new_prefix_lookup_table_length)
}


use std::cmp::min;
use std::ops::Range;
use std::sync::Arc;
use bitvec::prelude::*;
use rayon::iter::IntoParallelIterator;
use rayon::iter::IntoParallelRefIterator;
use rayon::iter::IndexedParallelIterator;
use rayon::iter::ParallelIterator;
use crate::builder::{SbwtIndexBuilder, BitPackedKmerSortingMem};
use crate::subsetseq::*;
use crate::sbwt::*;
use crate::atomic_bitmap::AtomicBitmap;
use super::interleaving::{MergeInterleaving, ThreeWayInterleaving,
    compute_piece_ranges, compute_piece_ranges_three_way, count_result_nodes_per_piece,
    count_result_nodes_per_piece_masked};
use super::common::{allocate_rows, transpose_and_concat_pieces, build_index,
    word_and_popcount_range, word_diff_nondummy_popcount_range};

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

/// Pass 1: for each (k-1)-suffix group that contributes a node to the intersection,
/// marks every LF-step target in `index1` colex space as having an incoming edge.
/// Each thread processes its own piece; the AtomicBitmap is shared for lock-free writes.
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
) -> AtomicBitmap {
    // Bit p = 1 iff s1 colex p receives at least one incoming edge in the intersection.
    // Bit 0 (root/dollar) is pre-set: it is structural and never needs an incoming edge.
    let has_incoming = AtomicBitmap::new(index1.n_sets());
    if isec_length > 0 { has_incoming.set(0, true); }

    thread_pool.install(|| {
        (0..n_threads).into_par_iter().for_each(|thread_idx| {
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
            // True iff the current group contains at least one shared (s1 && s2) position.
            let mut group_has_shared = false;

            for merged_colex in colex_range.clone() {
                if interleaving.is_leader[merged_colex] && merged_colex > colex_range.start {
                    // Flush: if any shared position exists, mark LF-step targets.
                    if group_has_shared {
                        for c in 0..sigma {
                            if s1_group_or[c] && s2_group_or[c] {
                                has_incoming.set(index1.lf_step(s1_first_in_group, c), true);
                            }
                        }
                    }
                    s1_group_read = false;
                    s2_group_read = false;
                    s1_group_or.fill(false);
                    s2_group_or.fill(false);
                    group_has_shared = false;
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
                if interleaving.s1[merged_colex] && interleaving.s2[merged_colex] {
                    group_has_shared = true;
                }

                s1_colex += interleaving.s1[merged_colex] as usize;
                s2_colex += interleaving.s2[merged_colex] as usize;
            }

            // Flush the last group in this piece.
            if group_has_shared {
                for c in 0..sigma {
                    if s1_group_or[c] && s2_group_or[c] {
                        has_incoming.set(index1.lf_step(s1_first_in_group, c), true);
                    }
                }
            }
        });
    });

    has_incoming
}

/// Collects s1-colex positions of real result k-mers that have no incoming edge and therefore
/// need auxiliary dummy chains.  When `difference` is false the result predicate is `s1 & s2`
/// (intersection); when true it is `s1 & !s2` (set difference).
/// The scan is split into `n_threads` pieces; each piece derives its starting s1_colex
/// via a prefix popcount on `interleaving.s1`, so all pieces run fully in parallel.
fn collect_result_source_nodes(
    interleaving: &MergeInterleaving,
    has_incoming: &AtomicBitmap,
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
                            && !interleaving.is_dummy[merged_colex]
                            && !has_incoming.get(s1_colex)
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

/// Builds a bitmap in `index1` colex space marking the rows that are real (non-dummy) k-mers
/// of the result: bit `p` is set iff `index1` row `p` is a result k-mer. When `difference` is
/// false the predicate is `s1 & s2 & !is_dummy` (intersection); when true it is
/// `s1 & !s2 & !is_dummy` (set difference).
///
/// The dummy-repair pass tracks `s1` colex positions but has no dummy marking of its own
/// ([ThreeWayInterleaving] does not compute one), so this carries the verdict from the two-way
/// interleaving into that pass. It costs `index1.n_sets()` bits, versus one dummy bitmap per
/// input if the three-way interleaving marked dummies itself.
///
/// The scan is split into `n_threads` pieces; each piece writes the contiguous slice of `s1`
/// colex space that it covers, and the pieces are concatenated in order.
fn compute_result_kmer_mask_s1(
    interleaving: &MergeInterleaving,
    difference: bool,
    n_threads: usize,
    thread_pool: &rayon::ThreadPool,
) -> BitVec {
    let merged_length = interleaving.s1.len();
    let piece_len = merged_length.div_ceil(n_threads);
    let pieces: Vec<Range<usize>> = (0..n_threads)
        .map(|t| t * piece_len..min((t + 1) * piece_len, merged_length))
        .collect();

    let piece_masks: Vec<BitVec> = thread_pool.install(|| {
        pieces.par_iter().map(|range| {
            let mut mask = BitVec::new();
            mask.resize(interleaving.s1[range.clone()].count_ones(), false);
            let mut piece_rel_s1_colex = 0usize;
            for merged_colex in range.clone() {
                if interleaving.s1[merged_colex] {
                    let in_result = if difference { !interleaving.s2[merged_colex] }
                                    else          {  interleaving.s2[merged_colex] };
                    if in_result && !interleaving.is_dummy[merged_colex] {
                        mask.set(piece_rel_s1_colex, true);
                    }
                    piece_rel_s1_colex += 1;
                }
            }
            mask
        }).collect()
    });

    crate::util::parallel_bitvec_concat(piece_masks)
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

/// Direct pass 2: builds the intersection SBWT bit-rows when every intersection
/// position already has an incoming edge and no dummy chain repair is required.
/// For each group the AND of the two per-index ORs is written at the group's isec leader.
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
fn intersect_rows_with_dummy_repair<SS: SubsetSeq + Send + Sync + Clone>(
    mut index1: SbwtIndex<SS>,
    index2: Arc<SbwtIndex<SS>>,
    source_colexes: Vec<usize>,
    sigma: usize,
    k: usize,
    optimize_peak_ram: bool,
    n_threads: usize,
    thread_pool: &rayon::ThreadPool,
) -> Vec<BitVec> {
    log::info!("[intersect] Dummy repair: reconstructing {} source k-mer(s) via inverse-LF (parallel)", source_colexes.len());
    // Build select before wrapping in Arc so no clone of the SubsetSeq is needed.
    index1.sbwt.build_select();
    let source_kmers = thread_pool.install(|| reconstruct_source_kmers(&index1, &source_colexes, k));
    drop(source_colexes);
    let index1 = Arc::new(index1);

    log::info!("[intersect] Dummy repair: building auxiliary SBWT from source k-mers");
    let (aux_submatrix, _) = SbwtIndexBuilder::<BitPackedKmerSortingMem>::new()
        .k(k)
        .n_threads(n_threads)
        .run_from_vecs(&source_kmers);
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
    let result_piece_counts =
        thread_pool.install(|| count_result_nodes_per_piece(&tw_ranges, &three_way, false));

    // For each (k-1)-suffix group: result node exists iff (s1&&s2)||s3.
    // Edge-bit c is set iff s3_or[c] || (s1_or[c] && s2_or[c]).
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
            let mut piece_rel_result_leader: Option<usize> = None;
            let mut result_colex_in_piece = 0usize;

            for merged_colex in colex_range.clone() {
                if three_way.is_leader[merged_colex] && merged_colex > colex_range.start {
                    // Flush: write edge bits at the result leader (if any).
                    if let Some(l) = piece_rel_result_leader {
                        for c in 0..sigma {
                            if s3_or[c] || (s1_or[c] && s2_or[c]) { piece_rows[c].set(l, true); }
                        }
                    }
                    s1_read = false; s2_read = false; s3_read = false;
                    s1_or.fill(false); s2_or.fill(false); s3_or.fill(false);
                    piece_rel_result_leader = None;
                }

                if !s1_read && three_way.s1[merged_colex] {
                    for c in 0..sigma { s1_or[c] = index1.sbwt.set_contains(s1_colex, c as u8); }
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

                let is_result = (three_way.s1[merged_colex] && three_way.s2[merged_colex])
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

            // Flush the last group.
            if let Some(l) = piece_rel_result_leader {
                for c in 0..sigma {
                    if s3_or[c] || (s1_or[c] && s2_or[c]) { piece_rows[c].set(l, true); }
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
/// whether any k-mers became source nodes in the intersection. If none did, pass 2 returns
/// the result directly. Otherwise, an auxiliary SBWT is built from those source
/// k-mers (with parallel k-mer reconstruction) and merged via a three-way interleaving to repair
/// the missing dummy chains.
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
    let has_incoming = pass1_has_incoming(
        &index1, &index2, &interleaving,
        isec_length, sigma, &piece_ranges, &s1_pops, &s2_pops,
        n_threads, &thread_pool,
    );

    // Count covered positions in parallel; bit 0 is pre-set (root), so complete == isec_length.
    let n_incoming: usize = thread_pool.install(|| {
        has_incoming.data.par_iter()
            .map(|w| w.load(std::sync::atomic::Ordering::Acquire).count_ones() as usize)
            .sum()
    });

    let new_rows = if n_incoming == isec_length {
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
        let source_colexes = collect_result_source_nodes(
            &interleaving, &has_incoming, merged_length, isec_length - n_incoming,
            false, n_threads, &thread_pool,
        );
        drop(interleaving); // No longer needed; free before the aux-index build.
        // Unwrap the Arc; refcount is 1 here since nothing else cloned it.
        let index1 = Arc::try_unwrap(index1).unwrap_or_else(|_| panic!("index1 Arc must be uniquely owned at this point"));
        intersect_rows_with_dummy_repair(
            index1, index2, source_colexes,
            sigma, k, optimize_peak_ram, n_threads, &thread_pool,
        )
    };

    log::info!("[intersect] Building subset rank structure");
    build_index(new_rows, n_kmers, k, n_threads, new_prefix_lookup_table_length)
}

// ── difference helpers ──────────────────────────────────────────────────────

/// Pass 1 for the set difference: marks every real difference k-mer that receives an incoming
/// edge from another real difference k-mer, in `index1` colex space.
///
/// Only real k-mers may vouch for their successors. The dummies of `index1` are not carried into
/// the result — the result's dummies are rebuilt from the auxiliary index — so a dummy must not
/// mark its successor as covered, or a real k-mer whose only `index1` predecessor is a dummy
/// would be left without a dummy chain of its own.
///
/// The unmarked real difference k-mers are exactly the ones needing fresh dummy chains, i.e. the
/// set `R'` of Definition 3 of the paper for the k-mer set `index1 \ index2`.
fn pass1_has_incoming_diff<SS: SubsetSeq + Send + Sync>(
    index1: &SbwtIndex<SS>,
    index2: &SbwtIndex<SS>,
    interleaving: &MergeInterleaving,
    sigma: usize,
    piece_ranges: &[Range<usize>],
    s1_pops: &[usize],
    s2_pops: &[usize],
    n_threads: usize,
    thread_pool: &rayon::ThreadPool,
) -> AtomicBitmap {
    let has_incoming = AtomicBitmap::new(index1.n_sets());
    // Do NOT pre-set position 0: the shared root ($$...$) is never a real difference k-mer,
    // and the root of the result comes from the auxiliary index instead.

    thread_pool.install(|| {
        (0..n_threads).into_par_iter().for_each(|thread_idx| {
            let colex_range = &piece_ranges[thread_idx];
            let mut s1_colex: usize = s1_pops[..thread_idx].iter().sum();
            let mut s2_colex: usize = s2_pops[..thread_idx].iter().sum();

            assert!(interleaving.is_leader[colex_range.start]);

            let mut s1_group_or = vec![false; sigma];
            let mut s2_group_or = vec![false; sigma];
            let mut s1_group_read = false;
            let mut s2_group_read = false;
            let mut s1_first_in_group = 0usize;
            // True iff the current group contains at least one real difference k-mer.
            let mut group_has_real_diff = false;

            for merged_colex in colex_range.clone() {
                if interleaving.is_leader[merged_colex] && merged_colex > colex_range.start {
                    if group_has_real_diff {
                        for c in 0..sigma {
                            // Edge exists iff successor is in s1 but not s2. A group holding a
                            // real k-mer has a $-free (k-1)-suffix, so all of its successors are
                            // real, and "in s1 but not s2" means exactly "is a difference k-mer".
                            if s1_group_or[c] && !s2_group_or[c] {
                                has_incoming.set(index1.lf_step(s1_first_in_group, c), true);
                            }
                        }
                    }
                    s1_group_read = false;
                    s2_group_read = false;
                    s1_group_or.fill(false);
                    s2_group_or.fill(false);
                    group_has_real_diff = false;
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
                // Only real (non-dummy) difference k-mers may vouch for their successors.
                if interleaving.s1[merged_colex] && !interleaving.s2[merged_colex]
                    && !interleaving.is_dummy[merged_colex] {
                    group_has_real_diff = true;
                }

                s1_colex += interleaving.s1[merged_colex] as usize;
                s2_colex += interleaving.s2[merged_colex] as usize;
            }

            // Flush the last group in this piece.
            if group_has_real_diff {
                for c in 0..sigma {
                    if s1_group_or[c] && !s2_group_or[c] {
                        has_incoming.set(index1.lf_step(s1_first_in_group, c), true);
                    }
                }
            }
        });
    });

    has_incoming
}

/// Direct pass 2 for the set difference: builds the SBWT bit-rows when no difference k-mer needs
/// a fresh dummy chain, which happens only when every difference k-mer has a difference k-mer as
/// its predecessor, i.e. when the result's de Bruijn graph consists of cycles (or is empty). The
/// padded spectrum then needs no dummies at all beyond the root, so the only rows are the real
/// difference k-mers plus a prepended, edge-free root.
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
        // Count real difference positions per piece: s1 & !s2 & !is_dummy
        let diff_piece_pops: Vec<usize> = piece_ranges.par_iter().map(|range| {
            word_diff_nondummy_popcount_range(
                interleaving.s1.as_raw_slice(),
                interleaving.s2.as_raw_slice(),
                interleaving.is_dummy.as_raw_slice(),
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
                if interleaving.s1[merged_colex] && !interleaving.s2[merged_colex]
                    && !interleaving.is_dummy[merged_colex] {
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

        // The SBWT invariant requires position 0 to be the root node ($$$…$). The root is a
        // dummy, so it is never among the real difference positions; we prepend it here. Its
        // successors are the dummies $$$…$c, none of which exist on this path, so it gets no
        // outgoing edges.
        for c in 0..sigma {
            let mut new_row = std::mem::take(&mut rows[c]);
            new_row.insert(0, false);
            rows[c] = new_row;
        }

        drop(index1);
        drop(index2);

        rows
    })
}

/// Dummy-repair pass 2 for the set difference. The result's real k-mers are the real difference
/// k-mers of `index1` (marked in `result_kmer_mask_s1`); all of its dummies come from the
/// auxiliary index built from the k-mers in `source_colexes`, whose padded spectrum supplies
/// exactly the minimal dummy chains the result needs.
fn difference_rows_with_dummy_repair<SS: SubsetSeq + Send + Sync + Clone>(
    mut index1: SbwtIndex<SS>,
    index2: Arc<SbwtIndex<SS>>,
    source_colexes: Vec<usize>,
    result_kmer_mask_s1: BitVec,
    sigma: usize,
    k: usize,
    optimize_peak_ram: bool,
    n_threads: usize,
    thread_pool: &rayon::ThreadPool,
) -> Vec<BitVec> {
    log::info!("[difference] Dummy repair: reconstructing {} source k-mer(s) via inverse-LF (parallel)", source_colexes.len());
    // Build select before wrapping in Arc so no clone of the SubsetSeq is needed.
    index1.sbwt.build_select();
    let source_kmers = thread_pool.install(|| reconstruct_source_kmers(&index1, &source_colexes, k));
    drop(source_colexes);
    let index1 = Arc::new(index1);

    log::info!("[difference] Dummy repair: building auxiliary SBWT from source k-mers");
    let (aux_submatrix, _) = SbwtIndexBuilder::<BitPackedKmerSortingMem>::new()
        .k(k)
        .n_threads(n_threads)
        .run_from_vecs(&source_kmers);
    drop(source_kmers);

    let aux_arc = Arc::new(convert_index::<SS>(aux_submatrix, n_threads, 0));

    log::info!("[difference] Dummy repair: building three-way interleaving");
    // Reuse the same three-way interleaving structure; only the result predicate changes below.
    let three_way = ThreeWayInterleaving::new(&index1, &index2, &aux_arc, optimize_peak_ram, n_threads);
    let three_way_len = three_way.s1.len();

    log::info!("[difference] Dummy repair: building difference SBWT bit-rows via three-way pass (parallel)");
    let (tw_ranges, s1_pc, s2_pc, s3_pc) =
        thread_pool.install(|| compute_piece_ranges_three_way(three_way_len, n_threads, &three_way));
    let result_piece_counts = thread_pool.install(|| {
        count_result_nodes_per_piece_masked(&tw_ranges, &three_way, &s1_pc, &result_kmer_mask_s1)
    });

    // For each (k-1)-suffix group: a result node exists iff s3 || (s1 && it is a real difference
    // k-mer).  Edge bits come from whichever of the two supplies the group's successors:
    //
    //   * A group holding a real difference k-mer has a $-free (k-1)-suffix, so its successors
    //     are real and are result nodes exactly when they are in index1 but not index2.
    //   * Any other group's successors are dummies, which exist only if the auxiliary index
    //     has them. The (s1_or && !s2_or) rule must NOT be applied here: it would set edges to
    //     index1 dummy chains that the result dropped, pointing at rows that do not exist.
    //
    // The two cases are disjoint: if a group holds a real difference k-mer then none of its
    // successors needs a fresh chain, so the auxiliary index contributes nothing to it.
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
            let mut piece_rel_result_leader: Option<usize> = None;
            let mut result_colex_in_piece = 0usize;
            // True iff the current group contains at least one real difference k-mer.
            let mut group_has_real_diff = false;

            for merged_colex in colex_range.clone() {
                if three_way.is_leader[merged_colex] && merged_colex > colex_range.start {
                    if let Some(l) = piece_rel_result_leader {
                        for c in 0..sigma {
                            let edge = if group_has_real_diff { s1_or[c] && !s2_or[c] }
                                       else                   { s3_or[c] };
                            if edge { piece_rows[c].set(l, true); }
                        }
                    }
                    s1_read = false; s2_read = false; s3_read = false;
                    s1_or.fill(false); s2_or.fill(false); s3_or.fill(false);
                    piece_rel_result_leader = None;
                    group_has_real_diff = false;
                }

                if !s1_read && three_way.s1[merged_colex] {
                    for c in 0..sigma { s1_or[c] = index1.sbwt.set_contains(s1_colex, c as u8); }
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

                let is_real_diff = three_way.s1[merged_colex] && result_kmer_mask_s1[s1_colex];
                group_has_real_diff |= is_real_diff;

                let is_result = is_real_diff || three_way.s3[merged_colex];
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

            // Flush the last group.
            if let Some(l) = piece_rel_result_leader {
                for c in 0..sigma {
                    let edge = if group_has_real_diff { s1_or[c] && !s2_or[c] }
                               else                   { s3_or[c] };
                    if edge { piece_rows[c].set(l, true); }
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
/// The result's **real** k-mers are the positions of the merged SBWT with `s1[i] && !s2[i] &&
/// !is_dummy[i]`. Its **dummy** k-mers, however, cannot be obtained by differencing the dummies
/// of the inputs: a dummy present in both inputs may still be needed by the result, and a dummy
/// exclusive to `index1` may be dead in it. They are therefore always rebuilt from scratch.
///
/// Pass 1 marks each real difference k-mer that has another real difference k-mer as its
/// predecessor. The unmarked ones are exactly those needing a dummy chain (the set `R'` of
/// Definition 3 of the paper, for the k-mer set `index1 \ index2`); they are reconstructed and
/// fed to the SBWT builder, and the padded spectrum of that auxiliary index supplies the
/// result's minimal dummy set, root included. A three-way interleaving of `index1`, `index2` and
/// the auxiliary index then assembles the rows.
///
/// The direct pass is taken only when no dummy chain is needed at all, which means the result is
/// either empty or made up entirely of cycles in the de Bruijn graph.
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

    // The number of sets in the result is not known until the dummy chains are known, so only
    // the k-mer count is computed here.
    let (diff_length, n_dummies) = thread_pool.install(|| compute_result_counts(&interleaving, true));
    let n_kmers = diff_length - n_dummies;
    log::info!("Differencing into {} distinct k-mers", n_kmers);

    let (piece_ranges, s1_pops, s2_pops) =
        thread_pool.install(|| compute_piece_ranges(merged_length, n_threads, &interleaving));

    log::info!("[difference] Pass 1: computing incoming-edge coverage (parallel)");
    let has_incoming = pass1_has_incoming_diff(
        &index1, &index2, &interleaving,
        sigma, &piece_ranges, &s1_pops, &s2_pops,
        n_threads, &thread_pool,
    );

    let source_colexes = collect_result_source_nodes(
        &interleaving, &has_incoming, merged_length,
        n_kmers, // upper-bound capacity hint
        true, n_threads, &thread_pool,
    );

    let new_rows = if source_colexes.is_empty() {
        log::info!("[difference] No difference k-mer needs a dummy chain; the result is cyclic or empty");
        log::info!("[difference] Pass 2: building difference SBWT bit-rows (parallel)");
        difference_rows_direct(
            index1, index2, interleaving,
            sigma, piece_ranges, s1_pops, s2_pops,
            n_threads, &thread_pool,
        )
    } else {
        log::info!("[difference] {} difference k-mer(s) need fresh dummy chains; merging with auxiliary index",
            source_colexes.len());
        let result_kmer_mask_s1 =
            compute_result_kmer_mask_s1(&interleaving, true, n_threads, &thread_pool);
        drop(interleaving);
        // Unwrap the Arc; refcount is 1 here since nothing else cloned it.
        let index1 = Arc::try_unwrap(index1).unwrap_or_else(|_| panic!("index1 Arc must be uniquely owned at this point"));
        difference_rows_with_dummy_repair(
            index1, index2, source_colexes, result_kmer_mask_s1,
            sigma, k, optimize_peak_ram, n_threads, &thread_pool,
        )
    };

    log::info!("[difference] Building subset rank structure");
    build_index(new_rows, n_kmers, k, n_threads, new_prefix_lookup_table_length)
}


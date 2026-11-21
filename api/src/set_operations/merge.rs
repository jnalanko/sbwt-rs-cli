use std::ops::Range;
use std::sync::Arc;
use bitvec::prelude::*;
use rayon::iter::IndexedParallelIterator;
use rayon::iter::IntoParallelIterator;
use rayon::iter::IntoParallelRefIterator;
use rayon::iter::IntoParallelRefMutIterator;
use rayon::iter::ParallelIterator;
use rayon::slice::ParallelSlice;
use simple_sds_sbwt::bit_vector::BitVector as SdsBitVector;
use simple_sds_sbwt::ops::{BitVec as _, Rank, Select};
use simple_sds_sbwt::raw_vector::{AccessRaw, RawVector};
use crate::atomic_bitmap::AtomicBitmap;
use crate::subsetseq::*;
use crate::sbwt::*;
use super::interleaving::{MergeInterleaving, compute_piece_ranges, leader_aligned_ranges};
use super::common::{allocate_rows, transpose_and_concat_pieces, build_index};

type BitVec = bitvec::vec::BitVec<u64, Lsb0>;

// ── merge helpers ──────────────────────────────────────────────────────────────

/// Cleanup, step 1 (mark): marks the directly redundant dummy nodes, i.e. those whose (k-1)-suffix group also
/// contains a real k-mer. Returns `None` if there are none (nothing to strip).
///
/// A dummy node is always the leader of its group: the k-mers of a group differ only in their
/// first character, and the dummy's first character `$` sorts first. A group also has at most one
/// dummy, since two dummies with the same (k-1)-suffix are the same k-mer. So a dummy is directly
/// redundant iff the next position is in the same group, i.e. is not a leader.
fn mark_directly_redundant(
    interleaving: &MergeInterleaving,
    thread_pool: &rayon::ThreadPool,
) -> Option<BitVec> {
    let n = interleaving.is_dummy.len();
    debug_assert!(interleaving.is_dummy.iter_ones().all(|p| interleaving.is_leader[p]),
        "a dummy node is not the leader of its group");

    let is_dummy = interleaving.is_dummy.as_raw_slice();
    let is_leader = interleaving.is_leader.as_raw_slice();
    let words: Vec<u64> = thread_pool.install(|| {
        (0..is_dummy.len()).into_par_iter().map(|i| {
            // Bit j of next_is_leader is is_leader[64*i + j + 1]
            let next_word = is_leader.get(i + 1).copied().unwrap_or(0);
            let next_is_leader = (is_leader[i] >> 1) | (next_word << 63);
            is_dummy[i] & !next_is_leader
        }).collect()
    });
    let mut directly_redundant = BitVec::from_vec(words);
    directly_redundant.truncate(n);
    directly_redundant.set_uninitialized(false); // Clear the bits past the end in the last word
    if n > 0 {
        // The last position has no next position, so its group ends there
        directly_redundant.set(n - 1, false);
    }

    let n_redundant = directly_redundant.count_ones();
    if n_redundant > 0 {
        log::info!("{} directly redundant dummy node(s) identified; will propagate and strip", n_redundant);
        Some(directly_redundant)
    } else {
        log::info!("No redundant dummy nodes in merged SBWT; skipping strip pass");
        None
    }
}

/// Builds the merged SBWT bit-rows from `index1` and `index2` according to the
/// interleaving. Each position writes OR-into-leader: edge bits at position `p` are forward-
/// OR'd into the leader of `p`'s (k-1)-suffix group. Consumes the index Arcs to free them
/// as early as possible once all pieces are done.
///
/// # Correctness invariant
/// After this function returns, exactly `interleaving.s1.len() - 1` bits are set across all rows.
/// (Each non-root node has exactly one incoming edge; the root has none.)
fn merge_rows<SS: SubsetSeq + Send + Sync>(
    index1: Arc<SbwtIndex<SS>>,
    index2: Arc<SbwtIndex<SS>>,
    interleaving: Arc<MergeInterleaving>,
    sigma: usize,
    n_threads: usize,
    thread_pool: &rayon::ThreadPool,
) -> Vec<BitVec> {
    thread_pool.install(|| {
        let (piece_ranges, s1_pops, s2_pops) =
            compute_piece_ranges(n_threads, &interleaving);

        let pieces_vecvec: Vec<Vec<BitVec>> = (0..piece_ranges.len()).into_par_iter().map(|piece_idx| {
            let colex_range = &piece_ranges[piece_idx];
            let mut new_rows = allocate_rows(sigma, colex_range.len());

            let mut s1_colex: usize = s1_pops[..piece_idx].iter().sum();
            let mut s2_colex: usize = s2_pops[..piece_idx].iter().sum();

            assert!(interleaving.is_leader[colex_range.start]);

            // Track the piece-relative index of the current group leader so every position
            // in the group OR-writes into its leader's row slot.
            let mut piece_rel_current_leader = 0;

            for merged_colex in colex_range.clone() {
                if interleaving.is_leader[merged_colex] {
                    piece_rel_current_leader = merged_colex - colex_range.start;
                }
                let in_s1 = interleaving.s1[merged_colex];
                let in_s2 = interleaving.s2[merged_colex];
                if in_s1 {
                    for c in 0..sigma {
                        if index1.sbwt.set_contains(s1_colex, c as u8) {
                            new_rows[c].set(piece_rel_current_leader, true);
                        }
                    }
                }
                if in_s2 {
                    for c in 0..sigma {
                        if index2.sbwt.set_contains(s2_colex, c as u8) {
                            new_rows[c].set(piece_rel_current_leader, true);
                        }
                    }
                }
                s1_colex += in_s1 as usize;
                s2_colex += in_s2 as usize;
            }

            assert_eq!(s1_colex, s1_pops[..=piece_idx].iter().sum());
            assert_eq!(s2_colex, s2_pops[..=piece_idx].iter().sum());

            new_rows
        }).collect();

        drop(index1);
        drop(index2);
        drop(interleaving);

        log::info!("[merge] Transposing and concatenating pieces");
        transpose_and_concat_pieces(pieces_vecvec, sigma)
    })
}

// ── public API ─────────────────────────────────────────────────────────────────

/// Checks that the SBWT invariant holds on `index` via rank queries on the built
/// SubsetSeq rank structure: sum_c rank(c, n_sets) + 1 == n_sets.
/// Compiled away in release builds.
#[cfg(debug_assertions)]
fn debug_assert_sbwt_invariant<SS: SubsetSeq>(index: &SbwtIndex<SS>, label: &str) {
    let n = index.n_sets();
    let sigma = index.alphabet().len();
    // rank() takes 0-based alphabet indices (0..sigma), not ASCII character values.
    let in_edges: usize = (0..sigma).map(|c| index.sbwt().rank(c as u8, n)).sum();
    assert_eq!(in_edges + 1, n,
        "{label} SBWT invariant violated after build_index: \
         {in_edges} in-edge(s) + root != {n} sets; {} node(s) have no in-coming edge",
        n - 1 - in_edges);
}

/// Merge `index1` and `index2` according to `interleaving`. After the merge, a [PrefixLookupTable] with
/// prefix length `new_prefix_lookup_table_length` will be added to new index. The number of threads used
/// in the merge is `n_threads`. Indexes are passed in as Arcs because if those are the only existing references,
/// this function can free the input SBWTs early which lowers the memory peak. Passing as Arc also allows
/// for use cases where the caller still wants to hold onto the sbwts: in that case dropping the Arcs
/// here will not free the memory. Same goes for the interleaving.
///
/// Redundant dummy nodes are removed: a dummy node is redundant if the merged index has a real k-mer that
/// makes it unnecessary. This includes redundant dummy nodes that were already in the inputs, e.g. in the
/// output of [merge_without_cleanup], so the result is the same as an index built from scratch. See
/// [merge_without_cleanup] to keep them.
pub fn merge<SS: SubsetSeq + Send + Sync>(
    index1: Arc<SbwtIndex<SS>>,
    index2: Arc<SbwtIndex<SS>>,
    interleaving: Arc<MergeInterleaving>,
    new_prefix_lookup_table_length: usize,
    n_threads: usize,
) -> SbwtIndex<SS> {
    merge_impl(index1, index2, interleaving, new_prefix_lookup_table_length, true, n_threads).0
}

/// Like [merge], but keeps the redundant dummy nodes introduced by the merge. The result is still
/// a valid SBWT with the same k-mers, but may have more dummy nodes. Colex position `i` of the
/// merged index corresponds exactly to position `i` of `interleaving`.
pub fn merge_without_cleanup<SS: SubsetSeq + Send + Sync>(
    index1: Arc<SbwtIndex<SS>>,
    index2: Arc<SbwtIndex<SS>>,
    interleaving: Arc<MergeInterleaving>,
    new_prefix_lookup_table_length: usize,
    n_threads: usize,
) -> SbwtIndex<SS> {
    merge_impl(index1, index2, interleaving, new_prefix_lookup_table_length, false, n_threads).0
}

/// Like [merge], but also returns which positions of `interleaving` were removed: a bit vector with one
/// bit per position of `interleaving`, set for the removed positions. Only dummy nodes are ever removed,
/// and the remaining positions keep their order: colex position `j` of the merged index corresponds to
/// the `j`-th unset bit.
pub fn merge_with_removed_positions<SS: SubsetSeq + Send + Sync>(
    index1: Arc<SbwtIndex<SS>>,
    index2: Arc<SbwtIndex<SS>>,
    interleaving: Arc<MergeInterleaving>,
    new_prefix_lookup_table_length: usize,
    n_threads: usize,
) -> (SbwtIndex<SS>, BitVec) {
    let (index, removed) = merge_impl(index1, index2, interleaving, new_prefix_lookup_table_length, true, n_threads);
    (index, removed.unwrap()) // Some because stripping was requested
}

/// Implementation of [merge], [merge_without_cleanup] and [merge_with_removed_positions].
///
/// If `strip_redundant_dummies` is true, the second return value is `Some(removed)`, where `removed`
/// has one bit per position of `interleaving`, set for the positions that were removed. Only dummy
/// nodes are ever removed, and the remaining positions keep their order: merged colex position `j`
/// corresponds to the `j`-th unset bit of `removed`. If `strip_redundant_dummies` is false, the second
/// return value is `None`.
pub(crate) fn merge_impl<SS: SubsetSeq + Send + Sync>(
    index1: Arc<SbwtIndex<SS>>,
    index2: Arc<SbwtIndex<SS>>,
    interleaving: Arc<MergeInterleaving>,
    new_prefix_lookup_table_length: usize,
    strip_redundant_dummies: bool,
    n_threads: usize,
) -> (SbwtIndex<SS>, Option<BitVec>) {
    let sigma = crate::util::DNA_ALPHABET.len();

    assert!(index1.k() == index2.k());
    let k = index1.k();
    let merged_length = interleaving.s1.len();
    let n_kmers = merged_length - interleaving.is_dummy.count_ones();
    log::info!("Merging into {} distinct k-mers", n_kmers);
    log::info!("Number of sets in merged SBWT: {}", merged_length);

    let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();

    let directly_redundant = if strip_redundant_dummies {
        log::info!("[merge] Marking directly redundant dummy nodes");
        mark_directly_redundant(&interleaving, &thread_pool)
    } else {
        log::info!("[merge] Keeping redundant dummy nodes");
        None
    };

    log::info!("[merge] Building the merged rows");
    let new_rows = merge_rows(
        index1, index2, Arc::clone(&interleaving),
        sigma, n_threads, &thread_pool,
    );

    // Verify the SBWT invariant on the freshly-built rows before any strip pass:
    // sum_c popcount(row[c]) must equal merged_length - 1.
    #[cfg(debug_assertions)]
    {
        let total_in_edges: usize = new_rows.iter().map(|r| r.count_ones()).sum();
        assert_eq!(total_in_edges + 1, merged_length,
            "[merge] merged rows violate SBWT invariant after merge_rows: \
             {} in-edge(s) + root != {} sets; {} node(s) have no in-coming edge",
            total_in_edges, merged_length, merged_length - 1 - total_in_edges);
    }

    if let Some(dr) = directly_redundant {
        log::info!("[merge] Removing redundant dummy nodes");
        let (new_rows, removed) = strip_redundant_from_rows(
            new_rows, n_kmers,
            dr, &interleaving.is_dummy, &interleaving.is_leader, &thread_pool,
        );
        drop(interleaving);
        // Verify the SBWT invariant is preserved by the strip pass.
        #[cfg(debug_assertions)]
        {
            let stripped_len = new_rows[0].len();
            let total_in_edges: usize = new_rows.iter().map(|r| r.count_ones()).sum();
            assert_eq!(total_in_edges + 1, stripped_len,
                "[merge] strip pass violated SBWT invariant: \
                 {} in-edge(s) + root != {} sets; {} node(s) have no in-coming edge",
                total_in_edges, stripped_len, stripped_len - 1 - total_in_edges);
        }
        log::info!("Building the subset rank structure");
        let result = build_index(new_rows, n_kmers, k, n_threads, new_prefix_lookup_table_length);
        #[cfg(debug_assertions)]
        debug_assert_sbwt_invariant(&result, "[merge/strip]");
        (result, Some(removed))
    } else {
        drop(interleaving);
        log::info!("Building the subset rank structure");
        let result = build_index(new_rows, n_kmers, k, n_threads, new_prefix_lookup_table_length);
        #[cfg(debug_assertions)]
        debug_assert_sbwt_invariant(&result, "[merge]");
        // Nothing was removed, either because stripping was not requested or because there was nothing to strip
        let removed = strip_redundant_dummies.then(|| bitvec![u64, Lsb0; 0; merged_length]);
        (result, removed)
    }
}

// Removes the redundant dummy nodes from `new_rows`. Returns the filtered rows and a bit vector
// marking the removed positions.
//
// A dummy node d is redundant when its (k-1)-suffix group also contains a real k-mer — that
// real k-mer provides a concrete predecessor for the lower-level node d pointed to, so d's
// dummy chain is no longer needed.
//
// The cleanup has three steps:
// 1. Mark – mark "directly redundant" dummies: those whose group also contains a real k-mer
//           (done in mark_directly_redundant, before the rows are built).
// 2. Propagate – BFS from the directly-redundant set using FL (inverse of lf_step).
//           FL(q) = select(c, q-C[c]) where c is the in-label of q.  A dummy leader p is
//           also redundant when *all* of its lf_step targets are redundant.  A colex scan
//           is NOT sufficient because lf_step(p,c) < p is not guaranteed.
// 3. Filter – rebuild each bit-row keeping only non-redundant positions.  When the
//           first surviving position of a group is not the original leader (because the
//           leader was removed), the leader's bits are OR-transferred to the new leader so
//           no outgoing-edge information is lost.
fn strip_redundant_from_rows(
    new_rows: Vec<BitVec>,
    n_kmers: usize,
    directly_redundant: BitVec,   // precomputed in merge_impl before interleaving was dropped
    is_dummy: &BitVec,
    is_leader: &BitVec,
    pool: &rayon::ThreadPool,
) -> (Vec<BitVec>, BitVec) {
    let n_threads = pool.current_num_threads();
    let sigma = crate::util::DNA_ALPHABET.len();

    // Step 1 (mark) was done in merge_impl; start with those marks.
    let mut redundant = directly_redundant;

    log::info!("[strip] Propagating redundancy up the dummy chains (BFS)");
    // The BFS needs rank (lf_step) and select (fl_step) on the rows. The rows are converted to
    // simple-sds bit vectors with rank support (rank9, 25% of a row) and select support
    // (select_mcl, a few percent of a row). The conversion copies a row, so it is done one row
    // at a time; the original row is freed as soon as it is converted.
    let mut rows: Vec<SdsBitVector> = new_rows.into_iter()
        .map(|row| SdsBitVector::from(crate::util::bitvec_to_simple_sds_raw_bitvec(row)))
        .collect();
    pool.install(|| rows.par_iter_mut().for_each(|row| {
        row.enable_rank();
        row.enable_select();
    }));

    // Build C-array: C[c] = 1 + #{in-labels < c} = 1 + sum_{c'<c} popcount(rows[c'])
    let c_count: Vec<usize> = rows.iter().map(|row| row.count_ones()).collect();
    let mut C = vec![0usize; sigma];
    C[0] = 1; // position 0 is the root (in-label '$')
    for c in 1..sigma {
        C[c] = C[c - 1] + c_count[c - 1];
    }

    // --- Step 2 (propagate): up the dummy chains via FL (inverse of lf_step) ---
    //
    // BFS from all directly-redundant positions using FL (inverse lf):
    //   Position q has in-label c iff C[c] <= q < C[c] + c_count[c].
    //   FL(q) = select(c, q - C[c]).

    let n_direct = redundant.count_ones();
    log::info!("[strip] BFS: {} directly redundant node(s); propagating up dummy chains (level-parallel)", n_direct);

    // Level-synchronous parallel BFS.
    //
    // Within one level every node in the frontier is already recorded in `redundant`.
    // The check for a candidate predecessor p only *reads* `redundant` and `rows`,
    // so the whole frontier can be evaluated in parallel.  After the parallel step we
    // mark the newly-redundant p's and use them as the next frontier.
    //
    // The frontiers are bit vectors over the merged nodes, so they take n bits each however
    // many nodes they hold (a list of positions would take 64 bits per node, which is more
    // than n bits in total on fragmented inputs with many dummies).  Two nodes of the frontier
    // can share the same predecessor p (if p has several outgoing edges into the frontier);
    // p then simply gets set twice in the next frontier.
    const CHUNK_WORDS: usize = 1024; // Frontier words per parallel task (65536 nodes)

    // Returns the predecessor p = fl_step(q) if it becomes redundant now that q is redundant.
    let newly_redundant_predecessor = |q: usize, redundant: &BitVec| -> Option<usize> {
        // --- inline fl_step ---
        if q == 0 { return None; }
        let c = (0..sigma).find(|&c| C[c] <= q && q < C[c] + c_count[c])?;
        let rank_of_q = q - C[c];
        let p = rows[c].select(rank_of_q)?;

        // The root (colex position 0) is never redundant: it must remain even with no out-edges.
        if p == 0 || redundant[p] || !is_dummy[p] || !is_leader[p] { return None; }
        let has_any_edge = (0..sigma).any(|c| rows[c].get(p));
        if !has_any_edge { return None; }
        let all_targets_redundant = (0..sigma)
            .filter(|&c| rows[c].get(p))
            .all(|c| {
                // --- inline lf_step ---
                let target = C[c] + rows[c].rank(p);
                redundant[target]
            });
        if all_targets_redundant { Some(p) } else { None }
    };

    let n = redundant.len();
    let mut frontier = redundant.clone();
    let mut frontier_size = n_direct;
    let mut bfs_levels = 0usize;
    let mut n_newly_redundant = 0usize;
    pool.install(|| {
        while frontier_size > 0 {
            bfs_levels += 1;
            log::info!("[strip] BFS level {}: {} node(s) in frontier", bfs_levels, frontier_size);

            let next = AtomicBitmap::new(n);
            frontier.as_raw_slice().par_chunks(CHUNK_WORDS).enumerate().for_each(|(chunk_idx, chunk)| {
                let chunk_start = chunk_idx * CHUNK_WORDS * 64;
                crate::util::for_each_one_bit(chunk, 0..chunk.len() * 64, |i| {
                    if let Some(p) = newly_redundant_predecessor(chunk_start + i, &redundant) {
                        next.set(p, true);
                    }
                });
            });
            frontier = BitVec::from_vec(next.data.into_iter().map(|w| w.into_inner()).collect());
            frontier.truncate(n);

            // Mark the new frontier as redundant
            redundant.as_raw_mut_slice().par_iter_mut().zip(frontier.as_raw_slice().par_iter())
                .for_each(|(r, f)| *r |= f);
            frontier_size = frontier.count_ones();
            n_newly_redundant += frontier_size;
        }
    });
    drop(frontier);
    log::info!("[strip] BFS complete: {} level(s), {} node(s) newly marked redundant", bfs_levels, n_newly_redundant);

    let n_redundant = redundant.count_ones();

    // --- Step 3 (filter): parallel filter rows ---
    //
    // Divide 0..merged_length into at most n_threads pieces, aligning each boundary to the next
    // is_leader position so that `pending` is always zero at a piece start.
    // Pre-compute initial rank_so_far[c] for each piece via O(1) rank queries on the rows.
    let piece_ranges = leader_aligned_ranges(is_leader, n_threads);
    let actual_n_pieces = piece_ranges.len();

    // Pre-compute per-piece initial rank_so_far values.
    let initial_ranks: Vec<Vec<usize>> = piece_ranges.iter().map(|range| {
        (0..sigma).map(|c| rows[c].rank(range.start)).collect()
    }).collect();
    // Filtering only reads bits: drop the rank and select support (this does not copy the bits)
    let rows: Vec<RawVector> = rows.into_iter().map(RawVector::from).collect();

    log::info!("[strip] Filtering rows: removing {} redundant dummy node(s), keeping {} real k-mers ({} parallel pieces)",
        n_redundant, n_kmers, actual_n_pieces);

    let pieces_vecvec: Vec<Vec<BitVec>> = pool.install(|| {
        (0..actual_n_pieces).into_par_iter().map(|t| {
            let Range { start, end } = piece_ranges[t];
            let mut rank_so_far = initial_ranks[t].clone();
            let mut pending = vec![false; sigma];
            let piece_keep_count = (end - start) - redundant[start..end].count_ones();
            let mut piece_filtered: Vec<BitVec> = (0..sigma).map(|_| {
                let mut bv = BitVec::new();
                bv.reserve(piece_keep_count);
                bv
            }).collect();

            for p in start..end {
                // lf_step target for (p, c) is C[c] + rank_so_far[c] (rank exclusive of p).
                if is_leader[p] {
                    pending.fill(false);
                }
                if redundant[p] {
                    for c in 0..sigma {
                        if rows[c].bit(p) && !redundant[C[c] + rank_so_far[c]] {
                            pending[c] = true;
                        }
                    }
                } else {
                    for c in 0..sigma {
                        let local_kept = rows[c].bit(p) && !redundant[C[c] + rank_so_far[c]];
                        piece_filtered[c].push(pending[c] || local_kept);
                    }
                    pending.fill(false);
                }
                // Advance rank counters after using them (rank is exclusive of p).
                for c in 0..sigma {
                    rank_so_far[c] += rows[c].bit(p) as usize;
                }
            }
            piece_filtered
        }).collect()
    });
    drop(rows); // Free the unfiltered rows before the filtered ones are concatenated

    (transpose_and_concat_pieces(pieces_vecvec, sigma), redundant)
}

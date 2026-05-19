use std::sync::Arc;
use bitvec::prelude::*;
use rayon::iter::IntoParallelIterator;
use rayon::iter::IntoParallelRefIterator;
use rayon::iter::ParallelIterator;
use crate::subsetseq::*;
use crate::sbwt::*;
use super::interleaving::{MergeInterleaving, compute_piece_ranges};
use super::common::{allocate_rows, transpose_and_concat_pieces, build_index};

type BitVec = bitvec::vec::BitVec<u64, Lsb0>;

// ── merge helpers ──────────────────────────────────────────────────────────────

/// Pass 1: for each (k-1)-suffix group that mixes exclusive-s1 dummies with exclusive-s2
/// real k-mers (or vice versa), marks those dummies as directly redundant.
/// Returns `None` when no such groups exist (nothing to strip).
fn mark_directly_redundant(
    interleaving: &MergeInterleaving,
    merged_length: usize,
    n_threads: usize,
    thread_pool: &rayon::ThreadPool,
) -> Option<BitVec> {
    let (piece_ranges, _, _) = thread_pool.install(|| compute_piece_ranges(merged_length, n_threads, interleaving));

    let pieces: Vec<(BitVec, bool)> = thread_pool.install(|| {
        piece_ranges.par_iter().map(|colex_range| {
            let mut piece_redundant = bitvec![u64, Lsb0; 0; colex_range.len()];
            let mut has_mixed = false;

            let mut pos = colex_range.start;
            while pos < colex_range.end {
                if interleaving.is_leader[pos] {
                    // First pass: determine whether the group has real k-mers exclusive to each index.
                    let mut group_has_s1_only_real = false;
                    let mut group_has_s2_only_real = false;
                    let mut gpos = pos;
                    loop {
                        let s1    = interleaving.s1[gpos];
                        let s2    = interleaving.s2[gpos];
                        let dummy = interleaving.is_dummy[gpos];
                        if s1 && !s2 && !dummy { group_has_s1_only_real = true; }
                        if s2 && !s1 && !dummy { group_has_s2_only_real = true; }
                        gpos += 1;
                        if gpos >= colex_range.end || interleaving.is_leader[gpos] { break; }
                    }
                    // Second pass: mark redundant dummies (only if needed).
                    if group_has_s1_only_real || group_has_s2_only_real {
                        has_mixed = true;
                        let mut gpos2 = pos;
                        loop {
                            let s1    = interleaving.s1[gpos2];
                            let s2    = interleaving.s2[gpos2];
                            let dummy = interleaving.is_dummy[gpos2];
                            if dummy {
                                if s1 && !s2 && group_has_s2_only_real {
                                    piece_redundant.set(gpos2 - colex_range.start, true);
                                }
                                if s2 && !s1 && group_has_s1_only_real {
                                    piece_redundant.set(gpos2 - colex_range.start, true);
                                }
                            }
                            gpos2 += 1;
                            if gpos2 >= colex_range.end || interleaving.is_leader[gpos2] { break; }
                        }
                    }
                    pos = gpos;
                } else {
                    pos += 1;
                }
            }
            (piece_redundant, has_mixed)
        }).collect()
    });

    if pieces.iter().any(|(_, mixed)| *mixed) {
        let directly_redundant = crate::util::parallel_bitvec_concat(
            pieces.into_iter().map(|(bv, _)| bv).collect()
        );
        log::info!("{} directly redundant dummy node(s) identified; will propagate and strip",
            directly_redundant.count_ones());
        Some(directly_redundant)
    } else {
        log::info!("No redundant dummy nodes in merged SBWT; skipping strip pass");
        None
    }
}

/// Pass 2: builds the merged SBWT bit-rows from `index1` and `index2` according to the
/// interleaving. Each position writes OR-into-leader: edge bits at position `p` are forward-
/// OR'd into the leader of `p`'s (k-1)-suffix group. Consumes the index Arcs to free them
/// as early as possible once all pieces are done.
///
/// # Correctness invariant
/// After this function returns, exactly `merged_length - 1` bits are set across all rows.
/// (Each non-root node has exactly one incoming edge; the root has none.)
fn merge_rows<SS: SubsetSeq + Send + Sync>(
    index1: Arc<SbwtIndex<SS>>,
    index2: Arc<SbwtIndex<SS>>,
    interleaving: Arc<MergeInterleaving>,
    sigma: usize,
    merged_length: usize,
    n_threads: usize,
    thread_pool: &rayon::ThreadPool,
) -> Vec<BitVec> {
    thread_pool.install(|| {
        let (piece_ranges, s1_pops, s2_pops) =
            compute_piece_ranges(merged_length, n_threads, &interleaving);

        let pieces_vecvec: Vec<Vec<BitVec>> = (0..n_threads).into_par_iter().map(|thread_idx| {
            let colex_range = &piece_ranges[thread_idx];
            let mut new_rows = allocate_rows(sigma, colex_range.len());

            let mut s1_colex: usize = s1_pops[..thread_idx].iter().sum();
            let mut s2_colex: usize = s2_pops[..thread_idx].iter().sum();

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

            assert_eq!(s1_colex, s1_pops[..=thread_idx].iter().sum());
            assert_eq!(s2_colex, s2_pops[..=thread_idx].iter().sum());

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
pub fn merge<SS: SubsetSeq + Send + Sync>(
    index1: Arc<SbwtIndex<SS>>,
    index2: Arc<SbwtIndex<SS>>,
    interleaving: Arc<MergeInterleaving>,
    new_prefix_lookup_table_length: usize,
    n_threads: usize,
) -> SbwtIndex<SS> {
    let sigma = crate::util::DNA_ALPHABET.len();

    assert!(index1.k() == index2.k());
    let k = index1.k();
    let merged_length = interleaving.s1.len();
    let n_kmers = merged_length - interleaving.is_dummy.count_ones();
    log::info!("Merging into {} distinct k-mers", n_kmers);
    log::info!("Number of sets in merged SBWT: {}", merged_length);

    let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();

    log::info!("[merge] Pass 1: scanning for directly redundant dummy nodes");
    let directly_redundant = mark_directly_redundant(&interleaving, merged_length, n_threads, &thread_pool);

    log::info!("[merge] Pass 2: building merged SBWT bit-rows (parallel)");
    let new_rows = merge_rows(
        index1, index2, Arc::clone(&interleaving),
        sigma, merged_length, n_threads, &thread_pool,
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

    log::info!("Building the subset rank structure");

    if let Some(dr) = directly_redundant {
        log::info!("[merge] Stripping redundant dummy nodes");
        let new_rows = strip_redundant_from_rows(
            new_rows, merged_length, n_kmers,
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
        let result = build_index(new_rows, n_kmers, k, n_threads, new_prefix_lookup_table_length);
        #[cfg(debug_assertions)]
        debug_assert_sbwt_invariant(&result, "[merge/strip]");
        result
    } else {
        drop(interleaving);
        let result = build_index(new_rows, n_kmers, k, n_threads, new_prefix_lookup_table_length);
        #[cfg(debug_assertions)]
        debug_assert_sbwt_invariant(&result, "[merge]");
        result
    }
}

/// Build a prefix popcount array for O(1) rank queries on a BitVec.
/// `prefix_ranks[w]` = popcount(bv[..w*64]).
fn build_prefix_ranks(bv: &BitVec) -> Vec<u64> {
    let words = bv.as_raw_slice();
    let mut prefix = Vec::with_capacity(words.len() + 1);
    let mut cum = 0u64;
    prefix.push(0u64);
    for &w in words {
        cum += w.count_ones() as u64;
        prefix.push(cum);
    }
    prefix
}

/// rank(p) = popcount(bv[..p]) using the precomputed prefix array. O(1).
#[inline]
fn rank_with_prefix(bv: &BitVec, prefix: &[u64], p: usize) -> usize {
    let wi = p / 64;
    let bi = p % 64;
    let base = prefix[wi] as usize;
    if bi == 0 { return base; }
    let word = bv.as_raw_slice()[wi];
    base + (word & ((1u64 << bi) - 1)).count_ones() as usize
}

/// select(r) = position of the r-th one-bit (0-indexed). O(log(n/64)) via binary search
/// on the prefix array, then a linear scan within the target word.
fn select_with_prefix(bv: &BitVec, prefix: &[u64], r: usize) -> Option<usize> {
    let n_words = bv.as_raw_slice().len();
    if n_words == 0 || r >= prefix[n_words] as usize { return None; }
    // Binary search for the word containing the r-th one-bit:
    // find largest w such that prefix[w] <= r and prefix[w+1] > r
    let mut lo = 0usize;
    let mut hi = n_words;
    while lo < hi {
        let mid = lo + (hi - lo) / 2;
        if (prefix[mid + 1] as usize) <= r {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    let remaining = r - prefix[lo] as usize;
    let mut w = bv.as_raw_slice()[lo];
    // Clear the first `remaining` one-bits to find the target bit
    for _ in 0..remaining {
        w &= w - 1;
    }
    Some(lo * 64 + w.trailing_zeros() as usize)
}

// Filters `new_rows` in-place to remove redundant dummy nodes introduced by the merge.
//
// A dummy node d from one index is redundant when its (k-1)-suffix group also contains a real
// k-mer from the *other* index — that real k-mer provides a concrete predecessor for the
// lower-level node d pointed to, so d's dummy chain is no longer needed.
//
// Pass 1  – mark "directly redundant" dummies: scan every group; for any group that has a
//           wrong-index dummy alongside a real k-mer from the other index, mark those dummies.
// Pass 2  – propagate: BFS from the directly-redundant set using FL (inverse of lf_step).
//           FL(q) = select(c, q-C[c]) where c is the in-label of q.  A dummy leader p is
//           also redundant when *all* of its lf_step targets are redundant.  A colex scan
//           is NOT sufficient because lf_step(p,c) < p is not guaranteed.
// Pass 3  – filter: rebuild each bit-row keeping only non-redundant positions.  When the
//           first surviving position of a group is not the original leader (because the
//           leader was removed), the leader's bits are OR-transferred to the new leader so
//           no outgoing-edge information is lost.
fn strip_redundant_from_rows(
    new_rows: Vec<BitVec>,
    merged_length: usize,
    n_kmers: usize,
    directly_redundant: BitVec,   // precomputed in merge_impl before interleaving was dropped
    is_dummy: &BitVec,
    is_leader: &BitVec,
    pool: &rayon::ThreadPool,
) -> Vec<BitVec> {
    let n_threads = pool.current_num_threads();
    let sigma = crate::util::DNA_ALPHABET.len();

    // Pass 1 was done in merge_impl; start with those marks.
    let mut redundant = directly_redundant;

    log::info!("[strip] Pass 2: BFS propagation of redundancy up dummy chains");
    // --- Pass 2: propagate backward up dummy chains ---
    // Build C-array: C[c] = 1 + #{in-labels < c} = 1 + sum_{c'<c} popcount(new_rows[c'])
    let mut C = vec![0usize; sigma];
    C[0] = 1; // position 0 is the root (in-label '$')
    for c in 1..sigma {
        C[c] = C[c - 1] + new_rows[c - 1].count_ones();
    }

    // Build prefix rank arrays for O(1) rank and O(log n) select queries.
    let prefix_ranks: Vec<Vec<u64>> = (0..sigma).map(|c| build_prefix_ranks(&new_rows[c])).collect();

    // --- Pass 2: propagate up dummy chains via FL (inverse of lf_step) ---
    //
    // A right-to-left colex scan is NOT sufficient because lf_step(p, c) < p is not
    // guaranteed — the target of a dummy can have a *larger* colex rank than its source
    // (e.g. `$ACAT` < `$$CAT` in colex even though `$$CAT` lf-steps to `$ACAT`).
    //
    // BFS from all directly-redundant positions using FL (inverse lf):
    //   Position q has in-label c iff C[c] <= q < C[c] + count_ones(new_rows[c]).
    //   FL(q) = select(c, q - C[c]).
    let c_count: Vec<usize> = (0..sigma).map(|c| new_rows[c].count_ones()).collect();

    let n_direct = redundant.count_ones();
    log::info!("[strip] BFS: {} directly redundant node(s); propagating up dummy chains (level-parallel)", n_direct);

    // Level-synchronous parallel BFS.
    //
    // Within one level every node in the frontier is already recorded in `redundant`.
    // The check for a candidate predecessor p only *reads* `redundant` and `new_rows`,
    // so the whole frontier can be evaluated in parallel.  After the parallel step we
    // collect the newly-redundant p's, mark them, and use them as the next frontier.
    //
    // Note: two different q's in the same level CAN share the same predecessor p (if p
    // has multiple outgoing edges that both land in the current frontier).  `dedup` on
    // the sorted candidate list handles that case — p is marked at most once.
    let mut frontier: Vec<usize> = redundant.iter_ones().collect();
    let mut bfs_levels = 0usize;
    let mut n_newly_redundant = 0usize;
    pool.install(|| {
        while !frontier.is_empty() {
            bfs_levels += 1;
            log::info!("[strip] BFS level {}: {} node(s) in frontier", bfs_levels, frontier.len());

            // For each q in the frontier, find the unique predecessor p = fl_step(q) and
            // check whether p should now be marked redundant.  All reads are from the
            // immutable-for-this-level `redundant` snapshot.
            let mut candidates: Vec<usize> = frontier
                .par_iter()
                .filter_map(|&q| {
                    // --- inline fl_step ---
                    if q == 0 { return None; }
                    let c = (0..sigma).find(|&c| C[c] <= q && q < C[c] + c_count[c])?;
                    let rank_of_q = q - C[c];
                    let p = select_with_prefix(&new_rows[c], &prefix_ranks[c], rank_of_q)?;

                    if redundant[p] || !is_dummy[p] || !is_leader[p] { return None; }
                    let has_any_edge = (0..sigma).any(|c| new_rows[c][p]);
                    if !has_any_edge { return None; }
                    let all_targets_redundant = (0..sigma)
                        .filter(|&c| new_rows[c][p])
                        .all(|c| {
                            // --- inline lf_step ---
                            let target = C[c] + rank_with_prefix(&new_rows[c], &prefix_ranks[c], p);
                            redundant[target]
                        });
                    if all_targets_redundant { Some(p) } else { None }
                })
                .collect();

            // Deduplicate: multiple q's may share the same predecessor p.
            candidates.sort_unstable();
            candidates.dedup();

            for &p in &candidates {
                redundant.set(p, true);
            }
            n_newly_redundant += candidates.len();
            frontier = candidates;
        }
    });
    log::info!("[strip] BFS complete: {} level(s), {} node(s) newly marked redundant", bfs_levels, n_newly_redundant);

    let n_redundant = redundant.count_ones();

    // --- Pass 3: parallel filter rows ---
    //
    // Divide 0..merged_length into n_threads pieces, aligning each boundary to the next
    // is_leader position so that `pending` is always zero at a piece start.
    // Pre-compute initial rank_so_far[c] for each piece via O(1) rank queries on the
    // prefix-rank arrays built in pass 2.
    let piece_len = merged_length.div_ceil(n_threads);
    let mut piece_starts: Vec<usize> = (0..n_threads).map(|t| t * piece_len).collect();
    // Align each boundary to the next leader (first piece starts at 0 which is always a leader).
    for t in 1..n_threads {
        while piece_starts[t] < merged_length && !is_leader[piece_starts[t]] {
            piece_starts[t] += 1;
        }
    }
    piece_starts.dedup(); // in case merging squashed some pieces
    let actual_n_pieces = piece_starts.len();
    let piece_ends: Vec<usize> = piece_starts[1..].iter().copied()
        .chain(std::iter::once(merged_length))
        .collect();

    // Pre-compute per-piece initial rank_so_far values.
    let initial_ranks: Vec<Vec<usize>> = piece_starts.iter().map(|&start| {
        (0..sigma).map(|c| rank_with_prefix(&new_rows[c], &prefix_ranks[c], start)).collect()
    }).collect();

    log::info!("[strip] Pass 3: filtering rows ({} redundant node(s) to remove, {} parallel pieces)",
        n_redundant, actual_n_pieces);
    log::info!("Removing {} redundant dummy node(s) from merged SBWT ({} real k-mers kept)",
        n_redundant, n_kmers);

    let pieces_vecvec: Vec<Vec<BitVec>> = pool.install(|| {
        (0..actual_n_pieces).into_par_iter().map(|t| {
            let start = piece_starts[t];
            let end = piece_ends[t];
            let mut rank_so_far = initial_ranks[t].clone();
            let mut pending = vec![false; sigma];
            let piece_keep_count = (start..end).filter(|&p| !redundant[p]).count();
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
                        if new_rows[c][p] && !redundant[C[c] + rank_so_far[c]] {
                            pending[c] = true;
                        }
                    }
                } else {
                    for c in 0..sigma {
                        let local_kept = new_rows[c][p] && !redundant[C[c] + rank_so_far[c]];
                        piece_filtered[c].push(pending[c] || local_kept);
                    }
                    pending.fill(false);
                }
                // Advance rank counters after using them (rank is exclusive of p).
                for c in 0..sigma {
                    rank_so_far[c] += new_rows[c][p] as usize;
                }
            }
            piece_filtered
        }).collect()
    });

    transpose_and_concat_pieces(pieces_vecvec, sigma)
}

// ── intersect helpers ──────────────────────────────────────────────────────────


// Converts an SbwtIndex<SubsetMatrix> to SbwtIndex<SS> by re-extracting the bit matrix.
// Used in the dummy-repair path of intersect when we need to merge the partial intersection (type SS)
// with an auxiliary index built by SbwtIndexBuilder (which always returns SubsetMatrix).
// Since SubsetMatrix is the only SubsetSeq implementation in practice, SS = SubsetMatrix always,
// and this conversion is effectively a no-op at runtime, but is needed to satisfy the type checker.

use std::ops::Range;
use std::sync::Arc;
use bitvec::prelude::*;
use rayon::iter::IntoParallelIterator;
use rayon::iter::IntoParallelRefIterator;
use rayon::iter::IndexedParallelIterator;
use rayon::iter::ParallelIterator;
use rayon::slice::ParallelSlice;
use crate::subsetseq::*;
use crate::sbwt::*;
use crate::atomic_bitmap::AtomicBitmap;
use super::interleaving::{MergeInterleaving, ThreeWayInterleaving,
    compute_piece_ranges, compute_piece_ranges_three_way, count_result_nodes_per_piece, is_result_node};
use super::dummy_chains::{needed_dummies, ExistingDummy, NeededDummies, NewDummy};
use super::common::{allocate_rows, transpose_and_concat_pieces, build_index,
    word_and_popcount_range, word_diff_popcount_range};

type BitVec = bitvec::vec::BitVec<u64, Lsb0>;

#[derive(Debug, Clone)]
pub struct MissingSelectSupport; // An error type

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

/// Pass 1 of both intersection and set difference.
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
fn pass1_has_incoming<SS: SubsetSeq + Send + Sync>(
    index1: &SbwtIndex<SS>,
    second: &SecondInput<'_, SS>,
    interleaving: &MergeInterleaving,
    result_length: usize,
    sigma: usize,
    piece_ranges: &[Range<usize>],
    s1_pops: &[usize],
    s2_pops: &[usize],
    thread_pool: &rayon::ThreadPool,
    difference: bool,
) -> (AtomicBitmap, Vec<usize>) {
    let has_incoming = AtomicBitmap::new(index1.n_sets());
    // Intersection: pre-set bit 0 (shared root); difference: omit to avoid inflating n_incoming.
    if !difference && result_length > 0 { has_incoming.set(0, true); }

    let dead_end_dummies_per_piece: Vec<Vec<usize>> = thread_pool.install(|| {
        (0..piece_ranges.len()).into_par_iter().map(|piece_idx| {
            let colex_range = &piece_ranges[piece_idx];
            let mut s1_colex: usize = s1_pops[..piece_idx].iter().sum();
            let mut s2_colex: usize = s2_pops[..piece_idx].iter().sum();

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
            let flush_group = |s1_or: &[bool], s2_or: &mut [bool], s1_first: usize| -> bool {
                second.complete_group_edges(index1, s1_first, s1_or, s2_or);
                let mut any_edge = false;
                for c in 0..sigma {
                    if keeps(s1_or[c], s2_or[c], difference) {
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
                        let any_edge = flush_group(&s1_group_or, &mut s2_group_or, s1_first_in_group);
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
                    second.read_group_edges(s2_colex, &mut s2_group_or);
                    s2_group_read = true;
                }
                let in_result = keeps(interleaving.s1[merged_colex], interleaving.s2[merged_colex], difference);
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
                let any_edge = flush_group(&s1_group_or, &mut s2_group_or, s1_first_in_group);
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

/// Collects s1-colex positions of result nodes that have no incoming edge and therefore
/// need new dummy chains.  When `difference` is false the result predicate is `s1 & s2`
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
    let pieces: Vec<Range<usize>> = crate::util::segment_range(0..merged_length, n_threads)
        .into_iter()
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
                        if keeps(true, interleaving.s2[merged_colex], difference)
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
/// the orphans are collected as source nodes and the repair rebuilds their padding chains —
/// re-creating exactly the marked ancestors they still need.
///
/// Returns the marked set in `index1`'s colex space. The dead ends are walked in parallel in
/// `thread_pool`: a walk stops at a node that another walk already marked, which goes on from
/// there, so the marked set does not depend on the order.
fn mark_dead_chains<SS: SubsetSeq + Sync>(
    index1: &SbwtIndex<SS>,
    dead_ends: &[usize],
    has_incoming: &AtomicBitmap,
    sigma: usize,
    thread_pool: &rayon::ThreadPool,
) -> BitVec {
    let marks = AtomicBitmap::new(index1.n_sets());
    thread_pool.install(|| dead_ends.par_iter().for_each(|&d| {
        let mut q = d;
        while !marks.test_and_set(q) {
            match index1.inverse_lf_step(q) {
                Some(p) if p != 0 => q = p,
                _ => break,
            }
        }
    }));
    let redundant = marks.into_bitvec_u64();

    // Orphan sweep (after all marking, so `redundant` is final). Every node has one incoming edge,
    // so each orphan is found once.
    const CHUNK_WORDS: usize = 1024; // Words of the marks per parallel task (65536 nodes)
    let (n_marked, n_orphans) = thread_pool.install(|| redundant.as_raw_slice().par_chunks(CHUNK_WORDS).enumerate().map(|(chunk, words)| {
        let (mut n_marked, mut n_orphans) = (0usize, 0usize);
        for (w, &word) in words.iter().enumerate() {
            let mut word = word;
            while word != 0 {
                let p = (chunk * CHUNK_WORDS + w) * 64 + word.trailing_zeros() as usize;
                word &= word - 1; // Clear the lowest set bit
                n_marked += 1;
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
        }
        (n_marked, n_orphans)
    }).reduce(|| (0, 0), |a, b| (a.0 + b.0, a.1 + b.1)));

    log::info!("Dead-chain removal: {} node(s) marked from {} dead end(s); {} orphan(s) get new dummy chains",
        n_marked, dead_ends.len(), n_orphans);
    redundant
}

/// Whether an element (a node or an edge) that is in index1 if `s1` and in index2 if `s2` is in
/// the result: for the intersection if it is in both, for the set difference if it is in index1
/// and not in index2.
#[inline]
fn keeps(s1: bool, s2: bool, difference: bool) -> bool {
    if difference { s1 && !s2 } else { s1 && s2 }
}

/// The second input of an intersection or a difference: another index, or, for a filter, the
/// nodes of index1 to keep.
enum SecondInput<'a, SS: SubsetSeq> {
    Index(Arc<SbwtIndex<SS>>),
    /// One bit per node of index1, set for the nodes to keep, dummy nodes included. A filter is
    /// the intersection of index1 and this input, which acts as an index that has exactly the kept
    /// nodes of index1: its positions in the merge plan are those nodes (see [filter]), and it has
    /// the edges of index1 that lead to a kept node.
    Keep(&'a BitVec),
}

impl<SS: SubsetSeq> SecondInput<'_, SS> {
    /// Reads the edge-bits of the group of the second input whose first node is at `s2_colex`.
    /// A keep mask has no rows of its own: see [SecondInput::complete_group_edges].
    fn read_group_edges(&self, s2_colex: usize, s2_or: &mut [bool]) {
        if let SecondInput::Index(index2) = self {
            for (c, bit) in s2_or.iter_mut().enumerate() { *bit = index2.sbwt.set_contains(s2_colex, c as u8); }
        }
    }

    /// The prefix of the log messages, the name of the result in them, and the verb of the
    /// operation: a filter is computed as an intersection, but logged as a filter.
    fn log_names(&self, difference: bool) -> (&'static str, &'static str, &'static str) {
        match (self, difference) {
            (SecondInput::Keep(_), _) => ("[filter]", "filter", "Filtering"),
            (SecondInput::Index(_), true) => ("[difference]", "difference", "Differencing"),
            (SecondInput::Index(_), false) => ("[intersect]", "intersection", "Intersecting"),
        }
    }

    /// For a keep mask, sets the edge-bits of the second input for the group of index1 whose first
    /// node is at `s1_first`, given the edge-bits `s1_or` of index1 there: an edge is in the second
    /// input if index1 has it and its target is kept. Nothing to do for an index.
    fn complete_group_edges(&self, index1: &SbwtIndex<SS>, s1_first: usize, s1_or: &[bool], s2_or: &mut [bool]) {
        if let SecondInput::Keep(kept) = self {
            for c in 0..s1_or.len() {
                s2_or[c] = s1_or[c] && kept[index1.lf_step(s1_first, c)];
            }
        }
    }
}

/// The merge plan of the dummy repair: `merge_plan` with the new dummy nodes of `dummies` inserted
/// as positions that are only in s3, and s3 also set at the positions of the needed dummy nodes
/// that exist, the root included.
///
/// A new node `$^i α` goes right before the first node that ends with α, which starts a
/// (k-1)-suffix group. With two or more dollars, its own (k-1)-suffix starts with a dollar, so it
/// starts a group of its own. With one dollar, its (k-1)-suffix is α, so it is the first node of
/// the group of that node, which then no longer starts the group. The new nodes before a node are
/// sorted by decreasing number of dollars, so one with one dollar comes right before it. The edges
/// of the s3 nodes are not copied: see [S3Edges].
///
/// The columns are built one at a time, and if `merge_plan` has no other reference, each of its
/// columns is freed once its new column is built (and is_dummy right away): at most two columns
/// of the merge plan and one column in pieces are then held next to the new columns, instead of
/// the whole merge plan, the new plan in pieces and a concatenated column.
fn insert_dummies(merge_plan: Arc<MergeInterleaving>, dummies: &NeededDummies, n_threads: usize) -> ThreeWayInterleaving {
    let n = merge_plan.s1.len();
    let pieces = crate::util::segment_range(0..n, n_threads);
    let (mut own_s1, mut own_s2, own_is_leader, shared) = match Arc::try_unwrap(merge_plan) {
        Ok(MergeInterleaving { s1, s2, is_leader, .. }) => (Some(s1), Some(s2), Some(is_leader), None),
        Err(shared) => {
            log::warn!("The merge plan has other references: its columns are freed only afterwards, which raises the peak memory of the dummy repair");
            (None, None, None, Some(shared))
        },
    };
    let old_is_leader = own_is_leader.as_ref().unwrap_or_else(|| &shared.as_ref().unwrap().is_leader);

    let s1 = {
        let old = own_s1.as_ref().unwrap_or_else(|| &shared.as_ref().unwrap().s1);
        insert_dummies_column(&pieces, dummies, old_is_leader, false, Some(old), false)
    };
    drop(own_s1.take());
    let s2 = {
        let old = own_s2.as_ref().unwrap_or_else(|| &shared.as_ref().unwrap().s2);
        insert_dummies_column(&pieces, dummies, old_is_leader, false, Some(old), false)
    };
    drop(own_s2.take());
    // The new dummy nodes, and the needed nodes that exist, each after the new nodes inserted
    // before it or before an earlier position
    let mut s3 = insert_dummies_column(&pieces, dummies, old_is_leader, true, None, false);
    for e in &dummies.existing {
        s3.set(e.pos + dummies.new.partition_point(|d| d.pos <= e.pos), true);
    }
    let is_leader = insert_dummies_column(&pieces, dummies, old_is_leader, true, Some(old_is_leader), true);
    ThreeWayInterleaving { s1, s2, s3, is_leader }
}

/// One column of the merge plan of [insert_dummies], built in parallel over `pieces` (ranges of
/// positions of the merge plan): `new_bit` at the new dummy nodes, and at the other positions the
/// bits of `old` (the column in the merge plan), or zeros without `old`. With `clear_joining`, a
/// position before which a new node with one dollar was inserted gets a zero instead: that node
/// starts the group of the position (for is_leader). `is_leader` is that of the merge plan.
///
/// The positions between two insertion points are copied as whole ranges of the old column, not
/// bit by bit: there are few new nodes (thousands to millions) in billions of positions.
fn insert_dummies_column(pieces: &[Range<usize>], dummies: &NeededDummies, is_leader: &BitVec, new_bit: bool, old: Option<&BitVec>, clear_joining: bool) -> BitVec {
    // Appends the positions `range` of the old column
    let copy = |part: &mut BitVec, range: Range<usize>| match old {
        Some(old) => part.extend_from_bitslice(&old[range]),
        None => part.resize(part.len() + range.len(), false),
    };
    let parts: Vec<BitVec> = pieces.par_iter().map(|range| {
        // The new dummies of this piece, which come before their position
        let first = dummies.new.partition_point(|d| d.pos < range.start);
        let last = dummies.new.partition_point(|d| d.pos < range.end);
        let mut part = BitVec::with_capacity(range.len() + last - first);
        let mut copied = range.start; // The positions before this one are in part
        let mut i = first;
        while i < last {
            let pos = dummies.new[i].pos;
            assert!(is_leader[pos], "a new dummy node is not before the first node of a group");
            copy(&mut part, copied..pos);
            let mut joins_group = false; // Whether the last new node before pos has one dollar
            while i < last && dummies.new[i].pos == pos {
                part.push(new_bit);
                joins_group = dummies.new[i].n_dollars == 1;
                i += 1;
            }
            copied = pos;
            if clear_joining && joins_group {
                part.push(false);
                copied = pos + 1;
            }
        }
        copy(&mut part, copied..range.end);
        debug_assert_eq!(part.len(), range.len() + last - first);
        part
    }).collect();
    crate::util::parallel_bitvec_concat(parts)
}

/// The out-edges of the s3 nodes of the merge plan of the dummy repair (see [insert_dummies]), read
/// from the needed dummy nodes without copying them to the positions of that merge plan. A position
/// that is in neither input is a new dummy node, and the new dummy nodes are in the order of
/// `dummies.new`; any other position p is position p - (the number of new nodes before it) of the
/// merge plan that the needed dummy nodes refer to.
struct S3Edges<'a> {
    dummies: &'a NeededDummies,
    n_new_before: usize, // The number of new dummy nodes before the current position
    next_existing: usize, // The first existing needed node at or after the current position
}

impl<'a> S3Edges<'a> {
    /// For reading from position `start` of the merge plan of the dummy repair on, in order.
    fn new(dummies: &'a NeededDummies, start: usize) -> Self {
        // The i-th new node is at position pos + i, after the i new nodes before it, so these
        // positions increase with i
        let new = &dummies.new;
        let (mut lo, mut hi) = (0, new.len());
        while lo < hi {
            let mid = (lo + hi) / 2;
            if new[mid].pos + mid < start { lo = mid + 1 } else { hi = mid }
        }
        // The position in the merge plan that the needed dummy nodes refer to
        let old_start = start - lo;
        Self { dummies, n_new_before: lo, next_existing: dummies.existing.partition_point(|e| e.pos < old_start) }
    }

    /// Writes the out-edges of the s3 node at `pos` into `edges` if `read`, and moves past `pos`,
    /// which is a new dummy node if `is_new`. Called for every position in order.
    fn advance(&mut self, pos: usize, read: bool, is_new: bool, edges: &mut [bool]) {
        if read {
            if is_new {
                let d = &self.dummies.new[self.n_new_before];
                debug_assert_eq!(d.pos + self.n_new_before, pos);
                for (c, e) in edges.iter_mut().enumerate() { *e = (d.edges >> c) & 1 == 1; }
            } else {
                let old = pos - self.n_new_before;
                // The positions come in order, so the node is found by moving forward in the list
                let existing = &self.dummies.existing;
                while existing[self.next_existing].pos < old {
                    self.next_existing += 1;
                }
                let node = existing[self.next_existing];
                debug_assert_eq!(node.pos, old, "an s3 node is neither new nor an existing needed node");
                for (c, e) in edges.iter_mut().enumerate() { *e = (node.edges >> c) & 1 == 1; }
            }
        }
        self.n_new_before += is_new as usize;
    }
}

/// Direct pass 2: builds the result SBWT bit-rows when every result position already has an
/// incoming edge and no dummy chain repair is required. The result positions are the positions of
/// the interleaving that the operation keeps (see [keeps]), and the root (position 0), which is in
/// both inputs: the intersection keeps it anyway, and the difference needs it as its root. For each
/// group, the edges that the operation keeps, from the ORs of the edge-bits of each index over the
/// group, are written at the group's first result position.
#[allow(clippy::too_many_arguments)]
fn rows_direct<SS: SubsetSeq + Send + Sync>(
    index1: Arc<SbwtIndex<SS>>,
    second: SecondInput<'_, SS>,
    interleaving: Arc<MergeInterleaving>,
    sigma: usize,
    piece_ranges: Vec<Range<usize>>,
    s1_pops: Vec<usize>,
    s2_pops: Vec<usize>,
    difference: bool,
    thread_pool: &rayon::ThreadPool,
) -> Vec<BitVec> {
    let (log_prefix, _, _) = second.log_names(difference);
    thread_pool.install(|| {
        let result_piece_pops: Vec<usize> = piece_ranges.par_iter().map(|range| {
            let (s1, s2) = (interleaving.s1.as_raw_slice(), interleaving.s2.as_raw_slice());
            if difference {
                // The root is in both inputs, so it is not counted by the predicate
                word_diff_popcount_range(s1, s2, range.clone()) + range.contains(&0) as usize
            } else {
                word_and_popcount_range(s1, s2, range.clone())
            }
        }).collect();

        let pieces_vecvec: Vec<Vec<BitVec>> = (0..piece_ranges.len()).into_par_iter().map(|piece_idx| {
            let colex_range  = &piece_ranges[piece_idx];
            let mut new_rows = allocate_rows(sigma, result_piece_pops[piece_idx]);

            let mut s1_colex: usize = s1_pops[..piece_idx].iter().sum();
            let mut s2_colex: usize = s2_pops[..piece_idx].iter().sum();

            assert!(interleaving.is_leader[colex_range.start]);

            // Per-group ORs of each index's edge-bits; reset to false at each group boundary,
            // because a group that is not in both inputs reads only one of them.
            let mut s1_group_or = vec![false; sigma];
            let mut s2_group_or = vec![false; sigma];
            let mut s1_group_read = false;
            let mut s2_group_read = false;
            // s1_colex of the first s1 element in the current group (used for lf_step).
            let mut s1_first_in_group = 0usize;
            // Piece-relative index of the first result position in the current group.
            let mut piece_rel_result_leader: Option<usize> = None;
            let mut result_colex_in_piece = 0usize;

            // Writes the kept edges of the current group at its result leader `l`
            let flush = |l: usize, s1_or: &[bool], s2_or: &mut [bool], s1_first: usize, rows: &mut [BitVec]| {
                second.complete_group_edges(&index1, s1_first, s1_or, s2_or);
                for c in 0..sigma {
                    if keeps(s1_or[c], s2_or[c], difference) { rows[c].set(l, true); }
                }
            };

            for merged_colex in colex_range.clone() {
                if interleaving.is_leader[merged_colex] && merged_colex > colex_range.start {
                    if let Some(l) = piece_rel_result_leader {
                        flush(l, &s1_group_or, &mut s2_group_or, s1_first_in_group, &mut new_rows);
                    }
                    s1_group_read = false;
                    s2_group_read = false;
                    s1_group_or.fill(false);
                    s2_group_or.fill(false);
                    piece_rel_result_leader = None;
                }

                if !s1_group_read && interleaving.s1[merged_colex] {
                    for c in 0..sigma { s1_group_or[c] = index1.sbwt.set_contains(s1_colex, c as u8); }
                    s1_first_in_group = s1_colex;
                    s1_group_read = true;
                }
                if !s2_group_read && interleaving.s2[merged_colex] {
                    second.read_group_edges(s2_colex, &mut s2_group_or);
                    s2_group_read = true;
                }
                if keeps(interleaving.s1[merged_colex], interleaving.s2[merged_colex], difference) || merged_colex == 0 {
                    if piece_rel_result_leader.is_none() {
                        piece_rel_result_leader = Some(result_colex_in_piece);
                    }
                    result_colex_in_piece += 1;
                }

                s1_colex += interleaving.s1[merged_colex] as usize;
                s2_colex += interleaving.s2[merged_colex] as usize;
            }

            // Flush the last group in this piece.
            if let Some(l) = piece_rel_result_leader {
                flush(l, &s1_group_or, &mut s2_group_or, s1_first_in_group, &mut new_rows);
            }

            assert_eq!(s1_colex, s1_pops[..=piece_idx].iter().sum());
            assert_eq!(s2_colex, s2_pops[..=piece_idx].iter().sum());
            assert_eq!(result_colex_in_piece, result_piece_pops[piece_idx]);

            new_rows
        }).collect();

        drop(index1);
        drop(second);
        drop(interleaving);

        log::info!("{} Transposing and concatenating pieces", log_prefix);
        transpose_and_concat_pieces(pieces_vecvec, sigma)
    })
}

/// Dummy-repair pass 2: some result k-mers are new source nodes (their predecessor k-mers exist
/// in index1 but not in the result), so their dummy predecessor chains are absent. Remedy:
/// reconstruct those source k-mers (in parallel), find the dummy nodes of their chains in the merge
/// plan (see [needed_dummies]), insert the new ones into it, and produce the complete bit-rows from
/// that merge plan, whose s3 are the dummy nodes of the chains. Edge rule: c is set iff a chain
/// has it, or index1 and index2 keep it (see [keeps]).
///
/// The same pass also leaves out the dead dummy chains: `redundant_s1` (index1 colex space,
/// precomputed by [`mark_dead_chains`]) marks every node of every chain whose real
/// k-mer did not survive the operation. Marked nodes are not emitted as result nodes and
/// edges whose LF-step target is marked are not written. This is safe even when a marked
/// node lies on a new source's padding chain: that node is in s3, so the `s3` term of the
/// predicate re-emits it and `s3_or` supplies both its row and its incoming edge.
///
/// If `record_result_nodes` is true, also returns the three-way interleaving and the positions of
/// the result nodes in it (one bit per position, set for the result nodes).
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
fn rows_with_dummy_repair<SS: SubsetSeq + Send + Sync + Clone>(
    index1: Arc<SbwtIndex<SS>>,
    second: SecondInput<'_, SS>,
    merge_plan: Arc<MergeInterleaving>,
    source_colexes: Vec<usize>,
    redundant_s1: Option<&BitVec>,
    sigma: usize,
    k: usize,
    difference: bool,
    record_result_nodes: bool,
    n_threads: usize,
    thread_pool: &rayon::ThreadPool,
) -> (Vec<BitVec>, Option<(ThreeWayInterleaving, BitVec)>) {
    let (log_prefix, result_name, _) = second.log_names(difference);
    log::info!("{} Dummy repair: reconstructing {} source k-mer(s) via inverse-LF (parallel)", log_prefix, source_colexes.len());
    let source_kmers = thread_pool.install(|| reconstruct_source_kmers(&index1, &source_colexes, k));
    drop(source_colexes);

    log::info!("{} Dummy repair: finding the dummy chains of the sources in the merge plan", log_prefix);
    // The root is needed also without sources (when only dead chains are removed): the difference
    // needs it, because the root is in both inputs, so it is in the difference only through s3.
    let dummies = match &second {
        SecondInput::Index(index2) => thread_pool.install(|| needed_dummies(&index1, Some(index2), &merge_plan.s1, Some(&merge_plan.s2), &source_kmers)),
        // The positions of the merge plan are the nodes of index1
        SecondInput::Keep(_) => thread_pool.install(|| needed_dummies(&index1, None, &merge_plan.s1, None, &source_kmers)),
    };
    drop(source_kmers);
    // The memory of the dummy repair, the peak of filter, in the log: what it consists of
    let bytes = |bits: &BitVec| bits.capacity() / 8;
    let dummies_bytes = dummies.new.capacity() * std::mem::size_of::<NewDummy>() + dummies.existing.capacity() * std::mem::size_of::<ExistingDummy>();
    let dead_bytes = redundant_s1.map_or(0, bytes);
    crate::util::log_memory(log_prefix, "after finding the dummy chains", &[
        ("merge plan", bytes(&merge_plan.s1) + bytes(&merge_plan.s2) + bytes(&merge_plan.is_dummy) + bytes(&merge_plan.is_leader)),
        ("needed dummies", dummies_bytes), ("dead chains", dead_bytes)]);

    log::info!("{} Dummy repair: inserting {} new dummy node(s) into the merge plan", log_prefix, dummies.new.len());
    // Frees the merge plan as it goes, if this is its only reference
    let three_way = thread_pool.install(|| insert_dummies(merge_plan, &dummies, n_threads));
    let three_way_bytes = bytes(&three_way.s1) + bytes(&three_way.s2) + bytes(&three_way.s3) + bytes(&three_way.is_leader);
    crate::util::log_memory(log_prefix, "after inserting the dummies", &[
        ("three-way merge plan", three_way_bytes), ("needed dummies", dummies_bytes), ("dead chains", dead_bytes)]);

    log::info!("{} Dummy repair: building {} SBWT bit-rows via three-way pass (parallel)", log_prefix, result_name);
    let (tw_ranges, s1_pc, s2_pc, s3_pc) =
        thread_pool.install(|| compute_piece_ranges_three_way(n_threads, &three_way));
    let result_piece_counts = thread_pool.install(||
        count_result_nodes_per_piece(&tw_ranges, &three_way, difference, redundant_s1, &s1_pc));

    // For each (k-1)-suffix group: result node exists iff (index1 and index2 keep it && !dead) || s3
    // (see is_result_node). Edge-bit c is set iff s3_or[c] || (keeps && target not dead).
    // If record_result_nodes is true, each piece also marks its result nodes.
    let pieces: Vec<(Vec<BitVec>, BitVec)> = thread_pool.install(|| {
        (0..tw_ranges.len()).into_par_iter().map(|piece_idx| {
            let colex_range = &tw_ranges[piece_idx];
            let mut piece_rows = allocate_rows(sigma, result_piece_counts[piece_idx]);
            let mut piece_result_marks = if record_result_nodes { BitVec::with_capacity(colex_range.len()) } else { BitVec::new() };

            let mut s1_colex: usize = s1_pc[..piece_idx].iter().sum();
            let mut s2_colex: usize = s2_pc[..piece_idx].iter().sum();
            let mut s3_colex: usize = s3_pc[..piece_idx].iter().sum();

            assert!(three_way.is_leader[colex_range.start]);

            // Per-group ORs; reset to false at each group boundary to prevent stale bleed-through
            // (stale s1/s2 could trigger spurious s3_or || keeps writes otherwise).
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
            let mut s3_edges = S3Edges::new(&dummies, colex_range.start);

            // Writes the edge bits of the current group at its result leader `l`. Edges into the
            // dead set can only originate from the root's group or a dead (revived through s3) group —
            // every dead node's parent is dead or the root — so all other groups skip the
            // per-edge lf_step target check.
            let flush = |l: usize, s1_or: &[bool], s2_or: &mut [bool], s3_or: &[bool], s1_first: usize, rows: &mut [BitVec]| {
                second.complete_group_edges(&index1, s1_first, s1_or, s2_or);
                let check_dead_targets = redundant_s1.filter(|red| s1_first == 0 || red[s1_first]);
                for c in 0..sigma {
                    if s3_or[c] || (keeps(s1_or[c], s2_or[c], difference)
                        && check_dead_targets.map_or(true, |red| !red[index1.lf_step(s1_first, c)]))
                    { rows[c].set(l, true); }
                }
            };

            for merged_colex in colex_range.clone() {
                if three_way.is_leader[merged_colex] && merged_colex > colex_range.start {
                    if let Some(l) = piece_rel_result_leader {
                        flush(l, &s1_or, &mut s2_or, &s3_or, s1_first_in_group, &mut piece_rows);
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
                    second.read_group_edges(s2_colex, &mut s2_or);
                    s2_read = true;
                }
                // The positions of the merge plan are in s1 or s2, the new dummy nodes in neither
                let is_new = !three_way.s1[merged_colex] && !three_way.s2[merged_colex];
                let read_s3 = !s3_read && three_way.s3[merged_colex];
                s3_edges.advance(merged_colex, read_s3, is_new, &mut s3_or);
                s3_read |= read_s3;

                // s1_colex (not yet incremented) is this node's index1 colex position; dead nodes
                // are not result nodes unless a dummy chain needs them (s3).
                let is_result = is_result_node(&three_way, merged_colex, s1_colex, difference, redundant_s1);
                if record_result_nodes {
                    piece_result_marks.push(is_result);
                }
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

            // Flush the last group in this piece.
            if let Some(l) = piece_rel_result_leader {
                flush(l, &s1_or, &mut s2_or, &s3_or, s1_first_in_group, &mut piece_rows);
            }

            assert_eq!(s1_colex, s1_pc[..=piece_idx].iter().sum());
            assert_eq!(s2_colex, s2_pc[..=piece_idx].iter().sum());
            assert_eq!(s3_colex, s3_pc[..=piece_idx].iter().sum());
            assert_eq!(result_colex_in_piece, result_piece_counts[piece_idx]);

            (piece_rows, piece_result_marks)
        }).collect()
    });

    let piece_bytes: usize = pieces.iter().map(|(rows, marks)| rows.iter().map(bytes).sum::<usize>() + bytes(marks)).sum();
    crate::util::log_memory(log_prefix, "after the three-way pass", &[
        ("three-way merge plan", three_way_bytes), ("rows and result marks in pieces", piece_bytes), ("needed dummies", dummies_bytes), ("dead chains", dead_bytes)]);

    drop(index1);
    drop(second);
    drop(dummies);
    let (row_pieces, result_mark_pieces): (Vec<Vec<BitVec>>, Vec<BitVec>) = pieces.into_iter().unzip();
    // Without recording, the interleaving and the empty marks are freed before the rows are
    // concatenated, which allocates the rows of the whole result
    let positions = record_result_nodes.then(|| (three_way, crate::util::parallel_bitvec_concat(result_mark_pieces)));

    log::info!("{} Dummy repair: transposing and concatenating pieces", log_prefix);
    let rows = thread_pool.install(|| transpose_and_concat_pieces(row_pieces, sigma));
    crate::util::log_memory(log_prefix, "after concatenating the rows", &[
        ("rows", rows.iter().map(bytes).sum()), ("kept three-way merge plan and result marks", positions.as_ref().map_or(0, |(_, marks)| three_way_bytes + bytes(marks)))]);
    (rows, positions)
}

// ── public API ─────────────────────────────────────────────────────────────────

/// The result of a set operation with its positions, which were requested, as returned by the
/// variants with removed positions: the merge plan and the positions that are not in the result.
fn with_removed_positions<SS: SubsetSeq>((index, positions): (SbwtIndex<SS>, Option<(ThreeWayInterleaving, BitVec)>)) -> (SbwtIndex<SS>, ThreeWayInterleaving, BitVec) {
    let (merge_plan, in_result) = positions.unwrap(); // Some because they were requested
    (index, merge_plan, !in_result)
}

/// Computes the intersection of two [SbwtIndex] structures. The intersection k-mers are those
/// present in **both** `index1` and `index2`. This mirrors [merge](crate::merge) in structure but uses AND
/// logic for incoming edges instead of OR, and restricts output positions to those shared by
/// both SBWTs.
/// 
/// IMPORTANT: index1 must have select support, otherwise returns an error.
///
/// For each (k-1)-suffix group in the merged SBWT, incoming-edge bits are OR-ed separately for
/// each index across all positions in the group, and the two resulting ORs are then AND-ed. The
/// result is stored at the first shared position of the group (the group's leader in the
/// intersection SBWT). This differs from [merge](crate::merge) where every position writes OR-into-leader directly.
///
/// Only dummy nodes shared by both SBWTs are retained. After pass 1, a parallel count detects
/// whether any k-mers became source nodes in the intersection, and pass 1 also reports the
/// dead-end dummies (shared dummies whose chain no longer leads to a surviving k-mer). If there
/// is neither, pass 2 returns the result directly. Otherwise, the dummy nodes of the missing
/// chains are found in the given interleaving (by searching the prefixes of the source k-mers in
/// the inputs), the new ones are inserted into it, and the rows are built from that, without the
/// dead-end chains.
///
/// `optimize_peak_ram` has no effect: the interleaving is given, and the repair computes none.
pub fn intersect<SS: SubsetSeq + Send + Sync + Clone>(
    index1: Arc<SbwtIndex<SS>>,
    index2: Arc<SbwtIndex<SS>>,
    interleaving: Arc<MergeInterleaving>,
    new_prefix_lookup_table_length: usize,
    optimize_peak_ram: bool,
    n_threads: usize,
) -> Result<SbwtIndex<SS>, MissingSelectSupport> {
    let _ = optimize_peak_ram; // No effect: the merge plan is given, and the repair computes none
    intersect_or_difference(index1, SecondInput::Index(index2), interleaving, new_prefix_lookup_table_length, false, false, n_threads)
        .map(|(index, _)| index)
}

/// Like [intersect], but also returns the interleaving that the result consists of, and which of
/// its positions are not in the result: a bit vector with one bit per position of the returned
/// interleaving, set for the removed positions. The remaining positions keep their order: colex
/// position `j` of the result is the `j`-th unset bit. Every node of the result, dummy nodes
/// included, is at one of the positions, so its position tells in which inputs it is.
///
/// If no dummy chains need to be repaired, the returned interleaving is `interleaving` with an
/// empty s3. Its bit vectors are moved out of the Arc if this is the only reference to it, and
/// copied otherwise. If dummy chains are repaired, it is `interleaving` with the new dummy nodes of
/// the chains inserted, and s3 marks the dummy nodes of the chains: a position that is only in s3
/// is a new dummy node, in neither input.
pub fn intersect_with_removed_positions<SS: SubsetSeq + Send + Sync + Clone>(
    index1: Arc<SbwtIndex<SS>>,
    index2: Arc<SbwtIndex<SS>>,
    interleaving: Arc<MergeInterleaving>,
    new_prefix_lookup_table_length: usize,
    optimize_peak_ram: bool,
    n_threads: usize,
) -> Result<(SbwtIndex<SS>, ThreeWayInterleaving, BitVec), MissingSelectSupport> {
    let _ = optimize_peak_ram; // No effect: the merge plan is given, and the repair computes none
    intersect_or_difference(index1, SecondInput::Index(index2), interleaving, new_prefix_lookup_table_length, false, true, n_threads)
        .map(with_removed_positions)
}

/// Computes the set difference `index1 \ index2`: the [SbwtIndex] containing exactly the k-mers
/// that are present in `index1` but **not** in `index2`.
/// 
/// IMPORTANT: index1 must have select support, otherwise returns an error.
///
/// The algorithm mirrors [`intersect`] but uses the complement predicate: a position in the
/// merged SBWT is a result node iff `s1[i] && !s2[i]`. Incoming-edge bits follow the same
/// group-OR logic but the AND is replaced by `s1_or[c] && !s2_or[c]`: an outgoing character `c`
/// is kept iff the successor k-mer (first character `c` followed by the current (k-1)-suffix)
/// is in `index1` but not in `index2`.
///
/// Like [`intersect`], a direct pass is taken when no difference k-mer becomes a source node
/// and no dummy became a dead end. Otherwise, the missing dummy chains are repaired and the dead
/// ones excluded exactly as in the intersection, with the edge rule
/// `s3_or[c] || (s1_or[c] && !s2_or[c])`, where `s3_or` are the edges of the dummy chains.
///
/// `optimize_peak_ram` has no effect, as for [intersect].
pub fn difference<SS: SubsetSeq + Send + Sync + Clone>(
    index1: Arc<SbwtIndex<SS>>,
    index2: Arc<SbwtIndex<SS>>,
    interleaving: Arc<MergeInterleaving>,
    new_prefix_lookup_table_length: usize,
    optimize_peak_ram: bool,
    n_threads: usize,
) -> Result<SbwtIndex<SS>, MissingSelectSupport> {
    let _ = optimize_peak_ram; // No effect: the merge plan is given, and the repair computes none
    intersect_or_difference(index1, SecondInput::Index(index2), interleaving, new_prefix_lookup_table_length, true, false, n_threads)
        .map(|(index, _)| index)
}

/// Like [difference], but also returns the interleaving that the result consists of, and which of
/// its positions are not in the result, as [intersect_with_removed_positions] does for the
/// intersection. The root (position 0) is in the result although it is in both inputs.
pub fn difference_with_removed_positions<SS: SubsetSeq + Send + Sync + Clone>(
    index1: Arc<SbwtIndex<SS>>,
    index2: Arc<SbwtIndex<SS>>,
    interleaving: Arc<MergeInterleaving>,
    new_prefix_lookup_table_length: usize,
    optimize_peak_ram: bool,
    n_threads: usize,
) -> Result<(SbwtIndex<SS>, ThreeWayInterleaving, BitVec), MissingSelectSupport> {
    let _ = optimize_peak_ram; // No effect: the merge plan is given, and the repair computes none
    intersect_or_difference(index1, SecondInput::Index(index2), interleaving, new_prefix_lookup_table_length, true, true, n_threads)
        .map(with_removed_positions)
}

/// Keeps the nodes of `index` that are marked in `keep` (one bit per node; the bits of the dummy
/// nodes are ignored): the result has exactly the k-mers of `index` whose bit is set, with the
/// dummy nodes they need, the same as an index built from scratch from those k-mers.
///
/// IMPORTANT: index must have select support, otherwise returns an error.
///
/// The filter is computed as an [intersect] of `index` and an input that has exactly the kept
/// k-mers, without building that input: a merge plan of `index` with itself gives the positions,
/// and the edges of that input are those of `index` that lead to a kept k-mer. The (k-1)-suffix
/// groups of `index` are taken from `group_starts` (one bit per node, set for the first node of
/// each group, e.g. from an LCS array), or computed from `index` if it is None. `keep` and
/// `group_starts` are taken over, so that they are freed as soon as they are no longer needed.
///
/// `optimize_peak_ram` has no effect, as for [intersect].
pub fn filter<SS: SubsetSeq + Send + Sync + Clone>(
    index: Arc<SbwtIndex<SS>>,
    keep: BitVec,
    group_starts: Option<BitVec>,
    new_prefix_lookup_table_length: usize,
    optimize_peak_ram: bool,
    n_threads: usize,
) -> Result<SbwtIndex<SS>, MissingSelectSupport> {
    let _ = optimize_peak_ram; // No effect: the merge plan is given, and the repair computes none
    filter_impl(index, keep, group_starts, new_prefix_lookup_table_length, false, n_threads)
        .map(|(index, _)| index)
}

/// Like [filter], but also returns the merge plan that the result consists of, and which of its
/// positions are not in the result, as [intersect_with_removed_positions] does. Its s1 are the
/// nodes of `index`, and its s2 are the kept nodes, dummy nodes included.
pub fn filter_with_removed_positions<SS: SubsetSeq + Send + Sync + Clone>(
    index: Arc<SbwtIndex<SS>>,
    keep: BitVec,
    group_starts: Option<BitVec>,
    new_prefix_lookup_table_length: usize,
    optimize_peak_ram: bool,
    n_threads: usize,
) -> Result<(SbwtIndex<SS>, ThreeWayInterleaving, BitVec), MissingSelectSupport> {
    let _ = optimize_peak_ram; // No effect: the merge plan is given, and the repair computes none
    filter_impl(index, keep, group_starts, new_prefix_lookup_table_length, true, n_threads)
        .map(with_removed_positions)
}

/// A bit vector with usize words, as computed for an index, with the u64 words of this module.
fn to_u64_words(bits: bitvec::vec::BitVec) -> BitVec {
    let n = bits.len();
    let mut result = BitVec::from_vec(bits.into_vec().into_iter().map(|w| w as u64).collect());
    result.truncate(n);
    result
}

#[allow(clippy::type_complexity)]
fn filter_impl<SS: SubsetSeq + Send + Sync + Clone>(
    index: Arc<SbwtIndex<SS>>,
    keep: BitVec,
    group_starts: Option<BitVec>,
    new_prefix_lookup_table_length: usize,
    return_positions: bool,
    n_threads: usize,
) -> Result<(SbwtIndex<SS>, Option<(ThreeWayInterleaving, BitVec)>), MissingSelectSupport> {
    let n = index.n_sets();
    assert_eq!(keep.len(), n);
    let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();
    let is_leader = match group_starts {
        Some(group_starts) => {
            assert_eq!(group_starts.len(), n);
            group_starts
        },
        None => {
            log::info!("[filter] Computing the (k-1)-suffix groups");
            to_u64_words(index.mark_k_minus_1_mers(n_threads))
        },
    };
    let is_dummy = to_u64_words(thread_pool.install(|| index.compute_dummy_node_marks()));

    // All dummy nodes are kept: the dummy chains that lead to no kept k-mer are removed as dead
    // chains, as in an intersection.
    let mut kept = keep;
    kept |= &is_dummy;
    // The word-wise counts of the merge plan also see the bits after the end of the last word
    kept.set_uninitialized(false);
    let mut s1 = bitvec![u64, Lsb0; 1; n];
    s1.set_uninitialized(false);
    let merge_plan = MergeInterleaving { s1, s2: kept.clone(), is_dummy, is_leader };
    intersect_or_difference(index, SecondInput::Keep(&kept), Arc::new(merge_plan), new_prefix_lookup_table_length, false, return_positions, n_threads)
}

/// Implementation of [intersect], [difference] (if `difference` is true) and [filter], and of their
/// variants with removed positions. If `return_positions` is true, also returns the interleaving
/// that the result consists of, and the positions of the result in it (one bit per position, set
/// for the result nodes).
#[allow(clippy::type_complexity, clippy::too_many_arguments)]
fn intersect_or_difference<SS: SubsetSeq + Send + Sync + Clone>(
    index1: Arc<SbwtIndex<SS>>,
    second: SecondInput<'_, SS>,
    interleaving: Arc<MergeInterleaving>,
    new_prefix_lookup_table_length: usize,
    difference: bool,
    return_positions: bool,
    n_threads: usize,
) -> Result<(SbwtIndex<SS>, Option<(ThreeWayInterleaving, BitVec)>), MissingSelectSupport> {

    if !index1.sbwt().has_select_support() {
        return Err(MissingSelectSupport)
    }

    let sigma = crate::util::DNA_ALPHABET.len();

    if let SecondInput::Index(index2) = &second {
        assert!(index1.k() == index2.k());
    }
    let k = index1.k();
    let merged_length = interleaving.s1.len();
    let (log_prefix, result_name, verb) = second.log_names(difference);

    let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();

    // The number of result positions (without the root for the difference) and of dummies among them
    let (result_length, n_dummies) = thread_pool.install(|| compute_result_counts(&interleaving, difference));
    let n_kmers = result_length - n_dummies;
    log::info!("{} into {} distinct k-mers", verb, n_kmers);
    log::info!("Number of sets in {} SBWT: {}", result_name, result_length);

    let (piece_ranges, s1_pops, s2_pops) =
        thread_pool.install(|| compute_piece_ranges(n_threads, &interleaving));

    log::info!("{} Pass 1: computing incoming-edge coverage (parallel)", log_prefix);
    let (has_incoming, dead_end_dummies) = pass1_has_incoming(
        &index1, &second, &interleaving,
        result_length, sigma, &piece_ranges, &s1_pops, &s2_pops,
        &thread_pool, difference,
    );
    log::info!("{} Dead-end dummies (no surviving outgoing edge): {}", log_prefix, dead_end_dummies.len());

    // Count covered positions in parallel. For the intersection, bit 0 is pre-set (root), so
    // complete == result_length. For the difference, every covered position is a result node (a
    // surviving edge's target is never in index2), so full coverage means no result k-mer became
    // a source node.
    let n_incoming: usize = thread_pool.install(|| {
        has_incoming.data.par_iter()
            .map(|w| w.load(std::sync::atomic::Ordering::Acquire).count_ones() as usize)
            .sum()
    });

    let (new_rows, positions) = if dead_end_dummies.is_empty() && n_incoming == result_length {
        log::info!("{} All {} {} positions have incoming edges; no dummy chains to repair", log_prefix, result_length, result_name);
        log::info!("{} Pass 2: building {} SBWT bit-rows (parallel)", log_prefix, result_name);
        if return_positions {
            let rows = rows_direct(
                index1, second, Arc::clone(&interleaving),
                sigma, piece_ranges, s1_pops, s2_pops, difference,
                &thread_pool,
            );
            // The result is exactly the positions that the operation keeps, and the root (see
            // rows_direct). The interleaving is moved out of the Arc if this is the only
            // reference, and otherwise only the needed vectors are copied.
            let mut in_result = interleaving.s2.clone();
            if difference {
                in_result = !in_result;
            }
            in_result &= &interleaving.s1;
            if merged_length > 0 {
                in_result.set(0, true);
            }
            let (s1, s2, is_leader) = match Arc::try_unwrap(interleaving) {
                Ok(il) => (il.s1, il.s2, il.is_leader),
                Err(shared) => (shared.s1.clone(), shared.s2.clone(), shared.is_leader.clone()),
            };
            let s3 = bitvec![u64, Lsb0; 0; s1.len()];
            (rows, Some((ThreeWayInterleaving { s1, s2, s3, is_leader }, in_result)))
        } else {
            // The rows function frees the interleaving before concatenating the rows
            let rows = rows_direct(
                index1, second, interleaving,
                sigma, piece_ranges, s1_pops, s2_pops, difference,
                &thread_pool,
            );
            (rows, None)
        }
    } else {
        log::info!("{} {} source(s) missing dummy chains; repairing them",
            log_prefix, result_length - n_incoming);

        // Mark the dead chains (unconditionally, up to the root) and clear the incoming-edge
        // bits of any orphaned children BEFORE collecting sources: orphans must be collected
        // so the repair rebuilds their padding chains. The rows pass of the repair excludes
        // marked nodes and the edges into them.
        let redundant_s1: Option<BitVec> = if dead_end_dummies.is_empty() { None } else {
            Some(mark_dead_chains(&index1, &dead_end_dummies, &has_incoming, sigma, &thread_pool))
        };
        drop(dead_end_dummies);

        let source_colexes = collect_result_source_nodes(
            &interleaving, &has_incoming, redundant_s1.as_ref(),
            merged_length, result_length - n_incoming, // capacity hint
            difference, n_threads, &thread_pool,
        );
        rows_with_dummy_repair(
            index1, second, interleaving, source_colexes, redundant_s1.as_ref(),
            sigma, k, difference, return_positions, n_threads, &thread_pool,
        )
    };

    log::info!("{} Building subset rank structure", log_prefix);
    Ok((build_index(new_rows, n_kmers, k, n_threads, new_prefix_lookup_table_length), positions))
}

#[cfg(test)]
mod tests {
    use rand_chacha::rand_core::{RngCore, SeedableRng};
    use rayon::prelude::*;

    use super::{insert_dummies_column, BitVec, NeededDummies};
    use crate::set_operations::dummy_chains::NewDummy;

    /// insert_dummies_column as it was, bit by bit: the reference for the range copies.
    fn by_bit(pieces: &[std::ops::Range<usize>], dummies: &NeededDummies, new_bit: bool, bit: impl Fn(usize, bool) -> bool + Sync) -> BitVec {
        let parts: Vec<BitVec> = pieces.par_iter().map(|range| {
            let first = dummies.new.partition_point(|d| d.pos < range.start);
            let last = dummies.new.partition_point(|d| d.pos < range.end);
            let mut next = dummies.new[first..last].iter().peekable();
            let mut part = BitVec::new();
            for pos in range.clone() {
                let mut joins_group = false;
                while let Some(d) = next.next_if(|d| d.pos == pos) {
                    part.push(new_bit);
                    joins_group = d.n_dollars == 1;
                }
                part.push(bit(pos, joins_group));
            }
            part
        }).collect();
        crate::util::parallel_bitvec_concat(parts)
    }

    /// A merge plan of `n` positions with random columns and groups, and `n_new` new dummy nodes
    /// before random group starts, several before some (by decreasing number of dollars).
    fn random_plan(n: usize, n_new: usize, seed: u64) -> (BitVec, BitVec, NeededDummies) {
        let mut rng = rand_chacha::ChaCha8Rng::seed_from_u64(seed);
        let column: BitVec = (0..n).map(|_| rng.next_u32() % 2 == 0).collect();
        let mut is_leader: BitVec = (0..n).map(|_| rng.next_u32() % 4 == 0).collect();
        if n > 0 { is_leader.set(0, true); }
        let leaders: Vec<usize> = is_leader.iter_ones().collect();
        let mut new: Vec<NewDummy> = (0..n_new).filter(|_| !leaders.is_empty()).map(|_| NewDummy {
            pos: leaders[rng.next_u32() as usize % leaders.len()],
            n_dollars: 1 + rng.next_u32() % 3,
            edges: 1,
        }).collect();
        new.sort_by_key(|d| (d.pos, std::cmp::Reverse(d.n_dollars)));
        new.dedup_by_key(|d| (d.pos, d.n_dollars));
        (column, is_leader, NeededDummies { new, existing: Vec::new() })
    }

    /// The range copies give the same columns as bit by bit, for each kind of column, with pieces
    /// that start and end anywhere, also between new dummy nodes of one position.
    #[test]
    fn insert_dummies_column_matches_by_bit() {
        let pool = rayon::ThreadPoolBuilder::new().num_threads(3).build().unwrap();
        pool.install(|| {
            for (n, n_new) in [(0, 0), (1, 0), (1, 3), (100, 0), (100, 40), (5000, 3000), (70000, 200)] {
                let (column, is_leader, dummies) = random_plan(n, n_new, n as u64);
                for n_pieces in [1, 3, 7] {
                    let pieces = crate::util::segment_range(0..n, n_pieces);
                    let case = format!("n {n}, {} new, {n_pieces} piece(s)", dummies.new.len());
                    assert_eq!(insert_dummies_column(&pieces, &dummies, &is_leader, false, Some(&column), false),
                        by_bit(&pieces, &dummies, false, |pos, _| column[pos]), "s1/s2, {case}");
                    assert_eq!(insert_dummies_column(&pieces, &dummies, &is_leader, true, None, false),
                        by_bit(&pieces, &dummies, true, |_, _| false), "s3, {case}");
                    assert_eq!(insert_dummies_column(&pieces, &dummies, &is_leader, true, Some(&is_leader), true),
                        by_bit(&pieces, &dummies, true, |pos, joins_group| is_leader[pos] && !joins_group), "is_leader, {case}");
                }
            }
        });
    }

    /// Times the range copies against bit by bit on a merge plan of 2^32 positions with a million
    /// new dummy nodes, on 32 threads: cargo test --release insert_dummies_column_speed -- --ignored --nocapture
    #[test]
    #[ignore]
    fn insert_dummies_column_speed() {
        let n = 1usize << 32;
        let pool = rayon::ThreadPoolBuilder::new().num_threads(32).build().unwrap();
        let (column, is_leader, dummies) = random_plan(n, 1_000_000, 1);
        let pieces = crate::util::segment_range(0..n, 32);
        pool.install(|| {
            let t = std::time::Instant::now();
            let fast = insert_dummies_column(&pieces, &dummies, &is_leader, false, Some(&column), false);
            let fast_s = t.elapsed().as_secs_f64();
            let t = std::time::Instant::now();
            let slow = by_bit(&pieces, &dummies, false, |pos, _| column[pos]);
            let slow_s = t.elapsed().as_secs_f64();
            assert_eq!(fast, slow);
            println!("{} positions, {} new dummy nodes: range copies {:.2} s, bit by bit {:.2} s", n, dummies.new.len(), fast_s, slow_s);
        });
    }
}

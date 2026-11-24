//! A k-mer lookup index that uses the suffix group marks of the SBWT instead of the LCS array.

use std::ops::Range;

use crate::sbwt::SbwtIndex;
use crate::subsetseq::SubsetSeq;
use crate::util::ACGT_TO_0123;

/// Looks up the k-mers of a batch of reads in an SBWT, using its (k-1)-suffix group marks
/// ([Sufgroups]) instead of an LCS array. The batching and prefetching are done inside
/// [KmerLookupIndex::streaming_search]: give it many reads at once.
/// The Copy trait is implemented because this struct contains only references to components of the
/// index, so it can be copied around cheaply.
#[derive(Clone, Copy, Eq, PartialEq, Debug)]
pub struct KmerLookupIndex<'a, SS: SubsetSeq> {
    pub sbwt: &'a SbwtIndex<SS>,
    pub sufgroups: &'a Sufgroups,
    /// Whether the index has only one orientation of each k-mer, so that a k-mer of a read is
    /// looked up in both orientations. Not supported yet.
    pub canonical: bool,
    pub n: usize, // The empty string has colex range [0,n)
    pub k: usize,
}

impl<'a, SS: SubsetSeq> KmerLookupIndex<'a, SS> {

    /// Create a new k-mer lookup index from an SBWT index and its suffix group marks.
    /// `canonical`: see [KmerLookupIndex::canonical] (not supported yet).
    pub fn new(sbwt: &'a SbwtIndex<SS>, sufgroups: &'a Sufgroups, canonical: bool) -> Self {
        assert_eq!(sufgroups.len(), sbwt.n_sets(), "the suffix group marks do not belong to the SBWT");
        assert!(!canonical, "one-orientation (canonical) indexes are not supported yet");
        KmerLookupIndex{sbwt, sufgroups, canonical, n: sbwt.n_sets(), k: sbwt.k()}
    }

    pub fn k(&self) -> usize {
        self.k
    }

    pub fn sbwt_len(&self) -> usize {
        self.n
    }

    /// For each read of `queries`, the colex rank of each of its k-mers in order, or None if the
    /// k-mer is not in the index (including k-mers with characters other than ACGT). A read
    /// shorter than k has no k-mers.
    pub fn streaming_search(&self, queries: &[&[u8]]) -> Vec<Vec<Option<usize>>> {
        // 32 reads in flight, their searches from scratch 16 at a time
        const READS: usize = 32;
        const BATCH_SIZE: usize = 16;
        let k = self.k;
        let mut result: Vec<Vec<Option<usize>>> = queries.iter().map(|q| vec![None; (q.len() + 1).saturating_sub(k)]).collect();
        restart_kmers_queued_fskip::<SS, &[u8], READS, BATCH_SIZE>(self.sbwt, self.sufgroups, queries, |r, s, x| {
            result[r][s] = Some(x);
        });
        result
    }
}

/// The (k-1)-suffix group marks of an SBWT: one bit per node, set for the first node of each group
/// of nodes whose k-mers have the same last k-1 characters (the marks of
/// [Dbg::k_minus_1_marks](crate::dbg::Dbg::k_minus_1_marks)). Built with the construction by
/// e.g. [run_with_sufgroups](crate::BitPackedKmerSortingMem::run_with_sufgroups).
///
/// The file format is that of the suffix group marks of Themisto 2 (`kmers.sufgroups`): the number
/// of bits as a little-endian u64, then the bits in 64-bit words, little-endian.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Sufgroups {
    marks: bitvec::vec::BitVec,
}

impl Sufgroups {

    /// Construct from the marks, one bit per node of the SBWT.
    pub fn new(marks: bitvec::vec::BitVec) -> Self {
        Self{marks}
    }

    /// The marks from the LCS array of an SBWT with k-mers of length `k`: node i is marked iff
    /// lcs\[i\] < k - 1, and node 0 (the root dummy) always.
    pub fn from_lcs(lcs: &crate::LcsArray, k: usize) -> Self {
        let mut marks: bitvec::vec::BitVec = (0..lcs.len()).map(|i| lcs.access(i) + 1 < k).collect();
        if !marks.is_empty() {
            marks.set(0, true); // Node 0 (the root dummy) starts the first suffix group
        }
        Self{marks}
    }

    /// The marks computed from the SBWT alone, by inverting its columns (as
    /// [Dbg](crate::dbg::Dbg) does without an LCS array). If `low_ram`, the labels of the inversion
    /// take 3 bits per node instead of a byte, and it takes about twice as long.
    pub fn from_sbwt<SS: SubsetSeq + Sync>(sbwt: &SbwtIndex<SS>, n_threads: usize, low_ram: bool) -> Self {
        Self{marks: if low_ram { sbwt.mark_k_minus_1_mers_low_memory(n_threads) } else { sbwt.mark_k_minus_1_mers(n_threads) }}
    }

    /// The marks, one bit per node of the SBWT.
    pub fn marks(&self) -> &bitvec::vec::BitVec {
        &self.marks
    }

    /// The marks, one bit per node of the SBWT, e.g. for
    /// [Dbg::new_with_k_minus_1_marks](crate::dbg::Dbg::new_with_k_minus_1_marks).
    pub fn into_marks(self) -> bitvec::vec::BitVec {
        self.marks
    }

    /// Returns the number of marks (the number of nodes of the SBWT).
    pub fn len(&self) -> usize {
        self.marks.len()
    }

    /// Returns true iff there are no marks.
    pub fn is_empty(&self) -> bool {
        self.marks.is_empty()
    }

    /// Returns the number of bytes used by the marks.
    pub fn size_in_bytes(&self) -> usize {
        self.marks.as_raw_slice().len() * std::mem::size_of::<usize>()
    }

    /// The suffix group of node `x`: from the last mark at or before `x` to the next mark after
    /// `x` (or the end). Requires that node 0 is marked, as it always is.
    #[inline]
    pub fn group(&self, x: usize) -> Range<usize> {
        const W: usize = usize::BITS as usize;
        let words = self.marks.as_raw_slice();
        let n = self.marks.len();
        let (w, b) = (x / W, x % W);
        // Bits 0..=b of word w
        let at_or_before = words[w] & (usize::MAX >> (W - 1 - b));
        let start = if at_or_before != 0 {
            w * W + W - 1 - at_or_before.leading_zeros() as usize
        } else {
            let mut w = w;
            loop {
                w -= 1; // Node 0 is marked, so this stops
                if words[w] != 0 {
                    break w * W + W - 1 - words[w].leading_zeros() as usize;
                }
            }
        };
        // Bits b+1.. of word w, then the next words (the end is at most n: the bits of the last
        // word past n are not looked at)
        let after = if b == W - 1 { 0 } else { words[w] & (usize::MAX << (b + 1)) };
        let end = if after != 0 {
            w * W + after.trailing_zeros() as usize
        } else {
            let mut w = w + 1;
            loop {
                if w * W >= n {
                    break n;
                }
                if words[w] != 0 {
                    break w * W + words[w].trailing_zeros() as usize;
                }
                w += 1;
            }
        };
        start..std::cmp::min(end, n)
    }

    /// Hints the CPU to load the marks around node `x`. Does not affect any results.
    #[inline]
    pub fn prefetch(&self, x: usize) {
        crate::util::prefetch_read(self.marks.as_raw_slice().as_ptr().wrapping_add(x / usize::BITS as usize));
    }

    /// Serialize the marks to the given writer.
    /// They can then later be loaded with [Sufgroups::load].
    /// Returns number of bytes written.
    pub fn serialize<W: std::io::Write>(&self, out: &mut W) -> std::io::Result<usize> {
        out.write_all(&(self.marks.len() as u64).to_le_bytes())?;
        let n_words = self.marks.len().div_ceil(64);
        for word in self.marks.as_raw_slice()[..n_words].iter() {
            out.write_all(&(*word as u64).to_le_bytes())?;
        }
        Ok(8 + 8 * n_words)
    }

    /// Load marks that were previously serialized with [Sufgroups::serialize].
    pub fn load<R: std::io::Read>(input: &mut R) -> std::io::Result<Self> {
        let mut buf = [0_u8; 8];
        input.read_exact(&mut buf)?;
        let len = u64::from_le_bytes(buf) as usize;
        let n_words = len.div_ceil(64);
        let mut words = Vec::<usize>::with_capacity(n_words);
        for _ in 0..n_words {
            input.read_exact(&mut buf)?;
            words.push(u64::from_le_bytes(buf) as usize);
        }
        let mut marks = bitvec::vec::BitVec::from_vec(words);
        marks.truncate(len);
        Ok(Self{marks})
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum FskipPhase {
    // A continuation or search from scratch of window `s` (the forward run).
    Run,
    // The look-ahead search of window `t` over the undecided windows `lo..`.
    Probe,
    // The search of window `w`, from right to left down to `lo`.
    Back,
    // The left gallop: window `bf_w` was found; galloping to the left for the first window of
    // its run, probe at `w` ...
    BackGallop,
    // ... bisection between `left` and `run_s`, probe at `w` ...
    BackBisect,
    // ... and the continuations from `run_s` towards `bf_w`, window `s` now.
    BackCont,
}

// The restart search with the look-ahead skip (fskip), for an index with both orientations of its
// k-mers, over many reads at once: READS reads are in flight, and their searches from scratch are
// batched, BATCH_SIZE at a time, with SbwtIndex::search_mini_batch_from_len. A read continues from
// a found window to the next one with one LF step from its k-mer, or from its suffix group if that
// fails. After a window that is not found, it looks ahead (a search that fails after l characters
// shows that the windows w + l - k..=w containing them are absent), searches the windows the
// look-ahead leaves from right to left, and a window found there starts a gallop to the left end
// of its run of found windows, which is then covered by continuations. Calls `on_kmer(r, s, x)`
// for every window s of read r that is in the index (not in order), with its colex rank x.
// When fewer than BATCH_SIZE searches are waiting and no read can continue (few reads, or the end
// of the reads), two or more waiting searches go out as one mini batch padded with repeats of its
// first search (their results are not used), and a single one is searched alone.
// Requires `BATCH_SIZE <= READS <= 64`.
#[allow(non_snake_case)]
fn restart_kmers_queued_fskip<'a, SS: SubsetSeq, R: AsRef<[u8]>, const READS: usize, const BATCH_SIZE: usize>(
    sbwt: &SbwtIndex<SS>,
    marks: &Sufgroups,
    reads: &'a [R],
    mut on_kmer: impl FnMut(usize, usize, usize),
) {
    assert!(BATCH_SIZE <= READS && READS <= 64, "requires BATCH_SIZE <= READS <= 64");
    let (n, k) = (sbwt.n_sets(), sbwt.k());
    let table = sbwt.get_lookup_table();
    let p = table.prefix_length;
    assert!(p <= k, "the prefix lookup table is longer than k");
    let extend = |I: &Range<usize>, c: usize| sbwt.lf_step(I.start, c)..sbwt.lf_step(I.end, c);

    let mut read_of = [0_usize; READS];
    let mut slot_read: [&'a [u8]; READS] = [&[]; READS];
    let mut phase = [FskipPhase::Run; READS];
    let mut s_ = [0_usize; READS]; // Run: current window
    let mut lo_ = [0_usize; READS]; // Probe/Back: first undecided window
    let mut t_ = [0_usize; READS]; // Probe: look-ahead window
    let mut w_ = [0_usize; READS]; // Back: current window
    let mut resume_x: [Option<usize>; READS] = [None; READS]; // Back: continue from the found probe t, ...
    let mut resume_run = [false; READS]; // ... or search t + 1 from scratch (after a targeted probe), or look ahead from t + 1
    let mut targeted = [false; READS]; // Probe: of the last window over the base of a failed continuation
    let mut run_search = 0_u64; // Run: the pending step is a search from scratch (not a continuation)
    let mut l_est = [15.min(k); READS];
    let mut prev_kmer = [0_usize; READS];
    let mut active = 0_u64;
    let mut cont = 0_u64;
    let mut group = 0_u64;
    let mut waiting = 0_u64;
    let mut search_win = [0_usize; READS];
    let mut table_idx: [Option<usize>; READS] = [None; READS];
    let mut next_read = 0;
    // The left gallop: the found window whose run is looked for and its rank, the first window
    // of the run found so far and its rank, the first window not known absent left of it, the
    // rightmost failed probe and its length, and the ranges left by a gap in a run.
    let mut bf_w = [0_usize; READS];
    let mut bf_x = [0_usize; READS];
    let mut run_s = [0_usize; READS];
    let mut run_s_x = [0_usize; READS];
    let mut left = [0_usize; READS];
    let mut gd = [0_usize; READS];
    let mut absent: [Option<(usize, usize)>; READS] = [None; READS];
    let mut pending: [Vec<(usize, usize)>; READS] = std::array::from_fn(|_| Vec::new());

    macro_rules! enqueue {
        ($q:expr, $win:expr) => {
            let (q, win): (usize, usize) = ($q, $win);
            waiting |= 1 << q;
            run_search |= 1 << q;
            search_win[q] = win;
            if p > 0 {
                table_idx[q] = table_index(slot_read[q], win, p);
                if let Some(idx) = table_idx[q] {
                    crate::util::prefetch_read(table.ranges.as_ptr().wrapping_add(idx));
                }
            }
        };
    }
    macro_rules! load {
        ($q:expr) => {
            let q: usize = $q;
            while next_read < reads.len() && reads[next_read].as_ref().len() < k {
                next_read += 1;
            }
            if next_read < reads.len() {
                (read_of[q], slot_read[q], s_[q], phase[q], l_est[q]) = (next_read, reads[next_read].as_ref(), 0, FskipPhase::Run, 15.min(k));
                pending[q].clear();
                next_read += 1;
                active |= 1 << q;
                enqueue!(q, 0);
            } else {
                active &= !(1 << q);
            }
        };
    }
    // Slot q continues its run at window s from k-mer x (or loads the next read after the last window).
    macro_rules! run_from {
        ($q:expr, $s:expr, $x:expr) => {
            let (q, s, x): (usize, usize, usize) = ($q, $s, $x);
            if s == slot_read[q].len() - k + 1 {
                load!(q);
            } else {
                (phase[q], s_[q], prev_kmer[q]) = (FskipPhase::Run, s, x);
                cont |= 1 << q;
                run_search &= !(1 << q);
                let c = ACGT_TO_0123[slot_read[q][s + k - 1] as usize];
                if c < 4 {
                    sbwt.prefetch_lf_step_interval(&(x..x + 1), c as usize);
                }
            }
        };
    }
    // Slot q searches window s from scratch in its run (or loads the next read after the last window).
    macro_rules! run_search_at {
        ($q:expr, $s:expr) => {
            let (q, s): (usize, usize) = ($q, $s);
            if s == slot_read[q].len() - k + 1 {
                load!(q);
            } else {
                (phase[q], s_[q]) = (FskipPhase::Run, s);
                enqueue!(q, s);
            }
        };
    }
    // Slot q looks ahead from the undecided window lo (or loads the next read if there is none).
    macro_rules! probe_from {
        ($q:expr, $lo:expr) => {
            let (q, lo): (usize, usize) = ($q, $lo);
            let m = slot_read[q].len() - k + 1;
            if lo >= m {
                load!(q);
            } else {
                let t = (lo + k.saturating_sub(l_est[q]).max(1)).min(m - 1);
                (phase[q], lo_[q], t_[q], targeted[q]) = (FskipPhase::Probe, lo, t, false);
                enqueue!(q, t);
            }
        };
    }
    // Slot q ends a right-to-left pass: continues from the found probe, or looks ahead after it.
    macro_rules! end_back {
        ($q:expr) => {
            let q: usize = $q;
            match resume_x[q] {
                Some(x) => {
                    run_from!(q, t_[q] + 1, x);
                }
                None if resume_run[q] => {
                    run_search_at!(q, t_[q] + 1);
                }
                None => {
                    probe_from!(q, t_[q] + 1);
                }
            }
        };
    }
    // The range lo_..= of the right-to-left pass of slot q is decided: the next range left by a
    // gap, or the end of the pass.
    macro_rules! back_range_done {
        ($q:expr) => {
            let q: usize = $q;
            match pending[q].pop() {
                Some((lo, hi)) => {
                    (phase[q], lo_[q], w_[q]) = (FskipPhase::Back, lo, hi);
                    enqueue!(q, hi);
                }
                None => {
                    end_back!(q);
                }
            }
        };
    }
    // Slot q has decided the run ending at bf_w: the pass goes on left of the windows the failed
    // probe next to it showed absent, or the range is done.
    macro_rules! back_after_run {
        ($q:expr) => {
            let q: usize = $q;
            let mut next = None;
            if let Some((a, l)) = absent[q] {
                l_est[q] = l;
                let c = (a + l).saturating_sub(k);
                if c > lo_[q] {
                    next = Some(c - 1);
                }
            }
            match next {
                Some(w2) => {
                    (phase[q], w_[q]) = (FskipPhase::Back, w2);
                    enqueue!(q, w2);
                }
                None => {
                    back_range_done!(q);
                }
            }
        };
    }
    // The next bisection probe of slot q, or the continuations from the first window of the run.
    macro_rules! back_bisect {
        ($q:expr) => {
            let q: usize = $q;
            if left[q] < run_s[q] {
                let mid = left[q] + (run_s[q] - left[q]) / 2;
                (phase[q], w_[q]) = (FskipPhase::BackBisect, mid);
                enqueue!(q, mid);
            } else {
                on_kmer(read_of[q], run_s[q], run_s_x[q]);
                if run_s[q] + 1 < bf_w[q] {
                    let (s, x) = (run_s[q] + 1, run_s_x[q]);
                    (phase[q], s_[q], prev_kmer[q]) = (FskipPhase::BackCont, s, x);
                    cont |= 1 << q;
                    run_search &= !(1 << q);
                    let c = ACGT_TO_0123[slot_read[q][s + k - 1] as usize];
                    if c < 4 {
                        sbwt.prefetch_lf_step_interval(&(x..x + 1), c as usize);
                    }
                } else {
                    if run_s[q] < bf_w[q] {
                        on_kmer(read_of[q], bf_w[q], bf_x[q]);
                    }
                    back_after_run!(q);
                }
            }
        };
    }
    // The left gallop of slot q from the found window bf_w: the next probe, or the bisection.
    macro_rules! back_gallop {
        ($q:expr) => {
            let q: usize = $q;
            if bf_w[q] - lo_[q] >= gd[q] {
                let probe = bf_w[q] - gd[q];
                (phase[q], w_[q]) = (FskipPhase::BackGallop, probe);
                enqueue!(q, probe);
            } else {
                left[q] = lo_[q];
                back_bisect!(q);
            }
        };
    }
    // Slot q starts the left gallop from window w, found with rank x.
    macro_rules! start_gallop {
        ($q:expr, $w:expr, $x:expr) => {
            let (q, w, x): (usize, usize, usize) = ($q, $w, $x);
            (bf_w[q], bf_x[q], run_s[q], run_s_x[q], absent[q], gd[q]) = (w, x, w, x, None, 1);
            back_gallop!(q);
        };
    }
    // The result of the pending continuation or search of slot q (Err(l): the first l characters
    // are in no k-mer; l = k for a failed continuation).
    macro_rules! handle {
        ($q:expr, $res:expr) => {
            let (q, res): (usize, Result<usize, usize>) = ($q, $res);
            match phase[q] {
                FskipPhase::Run => match res {
                    Ok(x) => {
                        on_kmer(read_of[q], s_[q], x);
                        run_from!(q, s_[q] + 1, x);
                    }
                    Err(l) => {
                        let s = s_[q];
                        let m = slot_read[q].len() - k + 1;
                        if run_search & (1 << q) != 0 {
                            // A failed search from scratch: a stretch not in the index.
                            l_est[q] = l;
                            probe_from!(q, s + 1);
                        } else if s + 1 >= m {
                            load!(q);
                        } else {
                            // A failed continuation: look ahead to the last window over its base.
                            let t = (s + k - 1).min(m - 1);
                            (phase[q], lo_[q], t_[q], targeted[q]) = (FskipPhase::Probe, s + 1, t, true);
                            enqueue!(q, t);
                        }
                    }
                },
                FskipPhase::Probe => {
                    let (lo, t) = (lo_[q], t_[q]);
                    match res {
                        Err(l) => {
                            l_est[q] = l;
                            let c = (t + l).saturating_sub(k);
                            if c > lo {
                                (phase[q], w_[q], resume_x[q], resume_run[q]) = (FskipPhase::Back, c - 1, None, targeted[q]);
                                enqueue!(q, c - 1);
                            } else if targeted[q] {
                                run_search_at!(q, t + 1);
                            } else {
                                probe_from!(q, t + 1);
                            }
                        }
                        Ok(x) => {
                            if t > lo {
                                resume_x[q] = Some(x);
                                start_gallop!(q, t, x);
                            } else {
                                on_kmer(read_of[q], t, x);
                                run_from!(q, t + 1, x);
                            }
                        }
                    }
                }
                // A window found starts the gallop to the left end of its run.
                FskipPhase::Back if res.is_ok() => {
                    let (w, x) = (w_[q], res.unwrap());
                    start_gallop!(q, w, x);
                }
                FskipPhase::BackGallop => match res {
                    Ok(x) => {
                        (run_s[q], run_s_x[q]) = (w_[q], x);
                        gd[q] = 2 * gd[q] + 1;
                        back_gallop!(q);
                    }
                    Err(l) => {
                        (absent[q], left[q]) = (Some((w_[q], l)), w_[q] + 1);
                        back_bisect!(q);
                    }
                },
                FskipPhase::BackBisect => {
                    match res {
                        Ok(x) => (run_s[q], run_s_x[q]) = (w_[q], x),
                        Err(l) => (absent[q], left[q]) = (Some((w_[q], l)), w_[q] + 1),
                    }
                    back_bisect!(q);
                }
                FskipPhase::BackCont => {
                    let (s, w) = (s_[q], bf_w[q]);
                    match res {
                        Ok(x) => {
                            on_kmer(read_of[q], s, x);
                            if s + 1 == w {
                                on_kmer(read_of[q], w, bf_x[q]);
                                back_after_run!(q);
                            } else {
                                s_[q] = s + 1;
                                prev_kmer[q] = x;
                                cont |= 1 << q;
                                let c = ACGT_TO_0123[slot_read[q][s + 1 + k - 1] as usize];
                                if c < 4 {
                                    sbwt.prefetch_lf_step_interval(&(x..x + 1), c as usize);
                                }
                            }
                        }
                        Err(_) => {
                            // A gap the probes missed: windows s + 1..w are decided later.
                            if s + 1 < w {
                                pending[q].push((s + 1, w - 1));
                            }
                            on_kmer(read_of[q], w, bf_x[q]);
                            back_after_run!(q);
                        }
                    }
                }
                // Not found: windows w + l - k..=w contain the first l characters of window w.
                FskipPhase::Back => {
                    let (lo, w) = (lo_[q], w_[q]);
                    let l = res.unwrap_err();
                    l_est[q] = l;
                    let c = (w + l).saturating_sub(k);
                    if c > lo {
                        w_[q] = c - 1;
                        enqueue!(q, c - 1);
                    } else {
                        back_range_done!(q);
                    }
                }
            }
        };
    }

    fn start_of<'b>(seq: &'b [u8], win: usize, idx: Option<usize>, ranges: &[Range<usize>], n: usize, p: usize, k: usize) -> (Range<usize>, &'b [u8]) {
        let I = if p == 0 { 0..n } else { idx.map_or(0..0, |idx| ranges[idx].clone()) };
        (I, &seq[win + p..win + k])
    }

    for q in 0..READS {
        load!(q);
    }
    while active != 0 {
        // One continuation round.
        let mut to_visit = cont;
        while to_visit != 0 {
            let q = to_visit.trailing_zeros() as usize;
            to_visit &= to_visit - 1;
            let read = slot_read[q];
            let x = prev_kmer[q];
            let c = ACGT_TO_0123[read[s_[q] + k - 1] as usize] as usize;
            let J = if c >= 4 {
                0..0
            } else if group & (1 << q) != 0 {
                group &= !(1 << q);
                extend(&marks.group(x), c)
            } else {
                let J = extend(&(x..x + 1), c);
                if J.is_empty() {
                    group |= 1 << q;
                    marks.prefetch(x);
                    continue;
                }
                J
            };
            cont &= !(1 << q);
            handle!(q, if J.is_empty() { Err(k) } else { Ok(J.start) });
        }
        // Searches from scratch.
        if waiting.count_ones() as usize >= BATCH_SIZE {
            let mut slots = [0_usize; BATCH_SIZE];
            for slot in slots.iter_mut() {
                *slot = waiting.trailing_zeros() as usize;
                waiting &= waiting - 1;
            }
            let results = {
                let starts: [(Range<usize>, &[u8]); BATCH_SIZE] = std::array::from_fn(|i| {
                    let q = slots[i];
                    start_of(slot_read[q], search_win[q], table_idx[q], &table.ranges, n, p, k)
                });
                let patterns: [&[u8]; BATCH_SIZE] = std::array::from_fn(|i| starts[i].1);
                let intervals: [Range<usize>; BATCH_SIZE] = std::array::from_fn(|i| starts[i].0.clone());
                sbwt.search_mini_batch_from_len(intervals, &patterns)
            };
            for (i, result) in results.into_iter().enumerate() {
                handle!(slots[i], result.map(|I| I.start).map_err(|c| p + c));
            }
        } else if cont == 0 && waiting.count_ones() >= 2 {
            // Fewer than BATCH_SIZE waiting and nothing else to do: one mini batch, padded with
            // repeats of its first search
            let count = waiting.count_ones() as usize;
            let mut slots = [0_usize; BATCH_SIZE];
            for slot in slots.iter_mut().take(count) {
                *slot = waiting.trailing_zeros() as usize;
                waiting &= waiting - 1;
            }
            for i in count..BATCH_SIZE {
                slots[i] = slots[0];
            }
            let results = {
                let starts: [(Range<usize>, &[u8]); BATCH_SIZE] = std::array::from_fn(|i| {
                    let q = slots[i];
                    start_of(slot_read[q], search_win[q], table_idx[q], &table.ranges, n, p, k)
                });
                let patterns: [&[u8]; BATCH_SIZE] = std::array::from_fn(|i| starts[i].1);
                let intervals: [Range<usize>; BATCH_SIZE] = std::array::from_fn(|i| starts[i].0.clone());
                sbwt.search_mini_batch_from_len(intervals, &patterns)
            };
            for (i, result) in results.into_iter().enumerate().take(count) {
                handle!(slots[i], result.map(|I| I.start).map_err(|c| p + c));
            }
        } else if cont == 0 {
            // A single search waiting: alone
            let mut to_visit = waiting;
            waiting = 0;
            while to_visit != 0 {
                let q = to_visit.trailing_zeros() as usize;
                to_visit &= to_visit - 1;
                let result = {
                    let (I, rest) = start_of(slot_read[q], search_win[q], table_idx[q], &table.ranges, n, p, k);
                    search_from_len(sbwt, I, rest)
                };
                handle!(q, result.map(|I| I.start).map_err(|c| p + c));
            }
        }
    }
}

// Index of the p-mer `read[s..s + p]` in the prefix lookup table (as in
// PrefixLookupTable::lookup), or None if it contains a character that is not A, C, G or T.
#[inline]
pub(crate) fn table_index(read: &[u8], s: usize, p: usize) -> Option<usize> {
    let mut idx = 0_usize;
    for (i, &c) in read[s..s + p].iter().rev().enumerate() {
        let c = ACGT_TO_0123[c as usize];
        if c >= 4 {
            return None;
        }
        idx |= (c as usize) << (2 * i);
    }
    Some(idx)
}

// [SbwtIndex::search_from], also returning for a pattern that is not found how many of its
// characters were consumed (0 if the starting interval is empty).
#[allow(non_snake_case)]
pub(crate) fn search_from_len<SS: SubsetSeq>(sbwt: &SbwtIndex<SS>, mut I: Range<usize>, pattern: &[u8]) -> Result<Range<usize>, usize> {
    if I.is_empty() {
        return Err(0);
    }
    for (i, &chr) in pattern.iter().enumerate() {
        let c = ACGT_TO_0123[chr as usize] as usize;
        if c >= 4 {
            return Err(i + 1);
        }
        I = sbwt.lf_step(I.start, c)..sbwt.lf_step(I.end, c);
        if I.is_empty() {
            return Err(i + 1);
        }
    }
    Ok(I)
}
#[cfg(test)]
mod tests {
    use super::*;

    // streaming_search gives, for every k-mer of every read, the colex rank that SbwtIndex::search
    // gives: for more reads than are in flight at once, reads from the index with substitutions
    // (failed continuations and look-aheads), reads half from the index and half random (the
    // right-to-left searches and the gallop), random reads, reads with an N, short and empty reads;
    // on several index types and k.
    #[test]
    fn streaming_search_matches_search() {
        use rand::{Rng, SeedableRng};
        let seqs: Vec<Vec<u8>> = (0..6_u64).map(|i| crate::util::gen_random_dna_string(300 + 31 * i as usize, 55 + i)).collect();
        let mut rng = rand::rngs::StdRng::seed_from_u64(7);
        let mut reads: Vec<Vec<u8>> = vec![vec![], b"ACG".to_vec()];
        for i in 0..150 {
            let seq = &seqs[i % seqs.len()];
            let len = rng.gen_range(20, 150);
            let start = rng.gen_range(0, seq.len() - len);
            let mut read = seq[start..start + len].to_vec();
            match i % 5 {
                0 => {} // From the index
                1 | 2 => {
                    // Substitutions at a rate of 1-5%
                    let rate = [0.01, 0.05][i % 2];
                    for c in read.iter_mut() {
                        if rng.gen::<f64>() < rate {
                            *c = b"ACGT"[rng.gen_range(0, 4)];
                        }
                    }
                }
                3 => {
                    // Half from the index, half random, in either order
                    let random = crate::util::gen_random_dna_string(len / 2, 1000 + i as u64);
                    if i % 2 == 0 { read.truncate(len / 2); read.extend(random); } else { read = [random, read[len / 2..].to_vec()].concat(); }
                }
                _ => {
                    if i % 2 == 0 {
                        read = crate::util::gen_random_dna_string(len, 2000 + i as u64); // Random
                    } else {
                        read[len / 2] = b'N';
                    }
                }
            }
            reads.push(read);
        }
        let queries: Vec<&[u8]> = reads.iter().map(|r| r.as_slice()).collect();

        fn check<SS: SubsetSeq>(sbwt: &SbwtIndex<SS>, sufgroups: &Sufgroups, queries: &[&[u8]], name: &str) {
            let k = sbwt.k();
            let result = KmerLookupIndex::new(sbwt, sufgroups, false).streaming_search(queries);
            assert_eq!(result.len(), queries.len());
            for (read, ranks) in queries.iter().zip(result.iter()) {
                let expected: Vec<Option<usize>> = if read.len() < k { vec![] } else {
                    read.windows(k).map(|kmer| sbwt.search(kmer).map(|interval| interval.start)).collect()
                };
                assert_eq!(*ranks, expected, "{}, k = {}, read {}", name, k, String::from_utf8_lossy(read));
            }
        }
        for k in [5, 21, 31, 32] {
            let (matrix, _, sufgroups) = crate::BitPackedKmerSortingMem::new_from_vecs(&seqs, k).add_rev_comp(true).run_with_sufgroups();
            check(&matrix, &sufgroups, &queries, "bit matrix");
            // Batches of fewer reads than are in flight (32), fewer than a mini batch (16), and
            // just over 32
            for n_reads in [1, 2, 15, 16, 31, 32, 33, 50] {
                for start in [0, 2, 40] {
                    check(&matrix, &sufgroups, &queries[start..start + n_reads], &format!("bit matrix, {} reads from {}", n_reads, start));
                }
            }
            let interleaved = crate::BitPackedKmerSortingMem::new_from_vecs(&seqs, k).add_rev_comp(true)
                .subset_seq::<crate::SubsetInterleavedMatrix>().run().0;
            check(&interleaved, &sufgroups, &queries, "interleaved matrix");
            let correction_sets = crate::BitPackedKmerSortingMem::new_from_vecs(&seqs, k).add_rev_comp(true)
                .subset_seq::<crate::SubsetCorrectionSets>().run().0;
            check(&correction_sets, &sufgroups, &queries, "correction sets");
        }
    }

    // Sufgroups::from_lcs and Sufgroups::from_sbwt (with and without low_ram) give the marks the
    // builders give.
    #[test]
    fn sufgroups_from_lcs_and_sbwt() {
        let seqs: Vec<Vec<u8>> = (0..6_u64).map(|i| crate::util::gen_random_dna_string(200 + 31 * i as usize, 31 + i)).collect();
        for k in [3, 5, 31] {
            for add_rev_comp in [false, true] {
                let (sbwt, lcs, sufgroups) = crate::BitPackedKmerSortingMem::new_from_vecs(&seqs, k)
                    .add_rev_comp(add_rev_comp).build_lcs(true).run_with_sufgroups();
                assert_eq!(Sufgroups::from_lcs(&lcs.unwrap(), k), sufgroups, "k = {}", k);
                assert_eq!(Sufgroups::from_sbwt(&sbwt, 2, false), sufgroups, "k = {}", k);
                assert_eq!(Sufgroups::from_sbwt(&sbwt, 2, true), sufgroups, "k = {}, low ram", k);
            }
        }
    }

    // The file of Sufgroups is that of the suffix group marks of Themisto 2: the number of bits as a
    // little-endian u64, then the bits in 64-bit words, little-endian (bit i in bit i % 64 of word
    // i / 64); and it loads back to the same marks.
    #[test]
    fn sufgroups_file_format() {
        for len in [0_usize, 1, 63, 64, 65, 130] {
            let marks: bitvec::vec::BitVec = (0..len).map(|i| i % 3 == 0 || i % 7 == 1).collect();
            let sufgroups = Sufgroups::new(marks.clone());
            let mut expected = (len as u64).to_le_bytes().to_vec();
            for w in 0..len.div_ceil(64) {
                let word = (w * 64..std::cmp::min(len, (w + 1) * 64)).filter(|&i| marks[i]).fold(0_u64, |word, i| word | 1 << (i % 64));
                expected.extend(word.to_le_bytes());
            }
            let mut buf = vec![];
            assert_eq!(sufgroups.serialize(&mut buf).unwrap(), buf.len());
            assert_eq!(buf, expected, "len {}", len);
            let loaded = Sufgroups::load(&mut buf.as_slice()).unwrap();
            assert_eq!(loaded, sufgroups);
            assert_eq!(loaded.len(), len);
            assert_eq!(loaded.size_in_bytes(), 8 * len.div_ceil(64));
        }
    }

    // Sufgroups::group gives the range from the last mark at or before x to the next mark after
    // x, for marks across word boundaries, dense and sparse, and with bits set past the end in the
    // last word (as a loaded file may have).
    #[test]
    fn sufgroups_group() {
        for len in [1_usize, 2, 63, 64, 65, 127, 128, 129, 300, 1000] {
            for density in [1, 3, 70, 500] {
                let mut marks: bitvec::vec::BitVec = (0..len).map(|i| i == 0 || (i * 7919) % density == 0).collect();
                let sufgroups = Sufgroups::new(marks.clone());
                for x in 0..len {
                    let start = (0..=x).rev().find(|&i| marks[i]).unwrap();
                    let end = (x + 1..len).find(|&i| marks[i]).unwrap_or(len);
                    assert_eq!(sufgroups.group(x), start..end, "len {}, density {}, x {}", len, density, x);
                    sufgroups.prefetch(x);
                }
                // Bits past the end of the marks in the last word do not count
                let words = len.div_ceil(64);
                marks.resize(words * 64, true);
                marks.truncate(len);
                let sufgroups = Sufgroups::new(marks.clone());
                let last = len - 1;
                let start = (0..=last).rev().find(|&i| marks[i]).unwrap();
                assert_eq!(sufgroups.group(last), start..len);
            }
        }
    }
}

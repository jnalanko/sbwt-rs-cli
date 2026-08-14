pub mod interleaving;
mod common;
mod merge;
mod intersect_difference;

pub use interleaving::MergeInterleaving;
pub use merge::merge;
pub use intersect_difference::{intersect, difference};

#[cfg(test)]
mod tests {

    use crate::{BitPackedKmerSortingDisk, BitPackedKmerSortingMem, SbwtIndexBuilder, SbwtIndex, SubsetMatrix, SubsetSeq};

    use super::*;
    use super::interleaving::split_to_pieces_par;
    use std::ops::Range;
    use std::sync::Arc;
    use bitvec::prelude::*;
    type BitSlice = bitvec::slice::BitSlice<u64, Lsb0>;

    #[test]
    fn split_to_pieces() {
        // 7 one-bits -> ceil(7/3) = 3 per piece
        //              e              e              e     e
        let s = bitvec![u64, Lsb0; 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1];
        let pieces = split_to_pieces_ref_implementation(&s, 3);
        assert_eq!(pieces.len(), 3);
        assert_eq!(pieces[0], (0, 0..5));
        assert_eq!(pieces[1], (2, 5..10));
        assert_eq!(pieces[2], (4, 10..12));
    }

    #[test]
    fn split_to_pieces_empty() {
        let s = bitvec![u64, Lsb0;];
        let pieces = split_to_pieces_ref_implementation(&s, 3);
        assert_eq!(pieces.len(), 3);
        assert_eq!(pieces[0], (0, 0..0));
        assert_eq!(pieces[1], (0, 0..0));
        assert_eq!(pieces[2], (0, 0..0));
    }

    #[test]
    fn split_to_pieces_run_out_of_pieces() {
        //              0  1  2  3  4  5  6  7  8  9  10 11
        let s = bitvec![u64, Lsb0; 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1];
        let pieces = split_to_pieces_ref_implementation(&s, 20);
        assert_eq!(pieces.len(), 20);
        assert_eq!(pieces[0], (0, 0..2));
        assert_eq!(pieces[1], (1, 2..3));
        assert_eq!(pieces[2], (1, 3..5));
        assert_eq!(pieces[3], (2, 5..8));
        assert_eq!(pieces[4], (4, 8..9));
        assert_eq!(pieces[5], (4, 9..10));
        assert_eq!(pieces[6], (4, 10..12));
        for i in 7..pieces.len(){
            assert_eq!(pieces[i], (5, 12..12));
        }
    }

    #[test]
    fn test_split_to_pieces_par() {
        // 7 one-bits -> ceil(7/3) = 3 per piece
        //              e              e              e     e
        let s = bitvec![u64, Lsb0; 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1];
        let pieces = split_to_pieces_ref_implementation(&s, 3);
        let pieces_par = split_to_pieces_par(&s, 3);
        assert_eq!(pieces, pieces_par);

        let s = bitvec![u64, Lsb0;];
        let pieces = split_to_pieces_ref_implementation(&s, 3);
        let pieces_par = split_to_pieces_par(&s, 3);
        assert_eq!(pieces, pieces_par);

        //              0  1  2  3  4  5  6  7  8  9  10 11
        let s = bitvec![u64, Lsb0; 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1];
        let pieces = split_to_pieces_ref_implementation(&s, 20);
        let pieces_par = split_to_pieces_par(&s, 20);
        assert_eq!(pieces, pieces_par);
    }

    // This is not a test. This is a single-threaded reference implementation.
    // Returns a segmentation s[l_1..r_1), s[l_2..r_2], ... such that
    // each segment ends in a 1-bit and has an approximately equal number of 1-bits.
    // Also returns the number of 0-bits before each segment.
    fn split_to_pieces_ref_implementation(s: &BitSlice, n_pieces: usize) -> Vec<(usize, Range<usize>)> {
        assert!(n_pieces > 0);
        if !s.is_empty() {
            // The last bit should always be 1
            assert!(s.last().unwrap() == true);
        }
        let s_popcount = s.count_ones();
        let ones_per_piece = s_popcount.div_ceil(n_pieces); // Last pieces may have fewer

        let mut s_ranges: Vec<Range<usize>> = vec![];
        let mut zeros_before: Vec<usize> = vec![];

        let mut cur_popcount = 0_usize;
        let mut n_zeros = 0_usize;
        for s_idx in s.iter_ones(){
            cur_popcount += 1;
            if cur_popcount == ones_per_piece || s_idx == s.len() - 1 {
                let prev_end = match s_ranges.last() {
                    Some(r) => r.end,
                    None => 0,
                };
                let new_end = s_idx + 1;
                s_ranges.push(prev_end..new_end);

                zeros_before.push(n_zeros);
                n_zeros += (new_end - prev_end) - cur_popcount;

                cur_popcount = 0;
            }
        }
        
        while s_ranges.len() < n_pieces {
            s_ranges.push(s.len()..s.len()); // Empty range
            zeros_before.push(n_zeros);
        }

        // Pair of vectors into vector of pairs
        zeros_before.into_iter().zip(s_ranges).collect()
        
    }

    // ── Merge helpers ──────────────────────────────────────────────────

fn check_merge(seq1: &[u8], seq2: &[u8], k: usize, n_threads: usize) {
        // Original logic: Build individual indices
        let (sbwt1, _) = SbwtIndexBuilder::<BitPackedKmerSortingDisk>::new()
            .k(k)
            .run_from_slices(vec![seq1].as_slice());
        let (sbwt2, _) = SbwtIndexBuilder::<BitPackedKmerSortingDisk>::new()
            .k(k)
            .run_from_slices(vec![seq2].as_slice());

        // Original logic: Build ground truth index
        let (sbwt_both, _) = SbwtIndexBuilder::<BitPackedKmerSortingDisk>::new()
            .k(k)
            .run_from_slices(vec![seq1, seq2].as_slice());

        // Verify interleaving consistency
        let inter_high_ram = MergeInterleaving::new(&sbwt1, &sbwt2, false, n_threads); 
        let inter_low_ram = MergeInterleaving::new(&sbwt1, &sbwt2, true, n_threads); 
        let inter_high_ram_1_thread = MergeInterleaving::new(&sbwt1, &sbwt2, false, 1);
        let inter_low_ram_1_thread = MergeInterleaving::new(&sbwt1, &sbwt2, true, 1);

        assert_eq!(inter_high_ram, inter_low_ram);
        assert_eq!(inter_high_ram, inter_high_ram_1_thread);
        assert_eq!(inter_high_ram, inter_low_ram_1_thread);

        // Perform merge
        let sbwt_merged = merge(Arc::new(sbwt1.clone()), Arc::new(sbwt2.clone()), Arc::new(inter_high_ram.clone()), 4, n_threads);

        // Compare spectra
        let spectrum_true = sbwt_both.reconstruct_padded_spectrum(n_threads);
        let spectrum_merged = sbwt_merged.reconstruct_padded_spectrum(n_threads);

        // Convert to chunks and sort to ensure deterministic comparison
        let mut true_kmers: Vec<Vec<u8>> = spectrum_true.chunks(k).map(|s| s.to_vec()).collect();
        let mut merged_kmers: Vec<Vec<u8>> = spectrum_merged.chunks(k).map(|s| s.to_vec()).collect();

        true_kmers.sort();
        merged_kmers.sort();

        assert_eq!(true_kmers, merged_kmers);
    }

    #[test]
    fn test_merge_original_random() {
        let seq1 = crate::util::gen_random_dna_string(1000, 42);
        let seq2 = crate::util::gen_random_dna_string(1000, 5235);
        check_merge(seq1.as_slice(), seq2.as_slice(), 5, 3);
    }

    #[test]
    fn test_merge_small_overlap() {
        let seq1 = b"ATCGATCGATCGATCG";
        let seq2 = b"GGGGATCGGGGGATCG";
        check_merge(seq1, seq2, 4, 1);
    }

    #[test]
    fn test_merge_no_overlap() {
        let seq1 = b"AAAAAAAAAAAAAA";
        let seq2 = b"TTTTTTTTTTTTTT";
        check_merge(seq1, seq2, 5, 3);
    }

    #[test]
    fn test_merge_redundant_dummy() {
        let seq1 = b"ACTGA";
        let seq2 = b"CTGAC";
        check_merge(seq1, seq2, 5, 1);
    }

    #[test]
    fn test_merge_redundant_dummy_long() {
        // Same redundant-dummy trigger as test_merge_redundant_dummy ("ACTGA" / "CTGAC"),
        // but both sequences are extended with long random tails so that the merged index
        // is large enough to exercise multi-threaded code paths.
        let tail1 = crate::util::gen_random_dna_string(2000, 101);
        let tail2 = crate::util::gen_random_dna_string(2000, 202);
        let seq1: Vec<u8> = b"ACTGA".iter().chain(tail1.iter()).copied().collect();
        let seq2: Vec<u8> = b"CTGAC".iter().chain(tail2.iter()).copied().collect();
        for n_threads in [1, 2, 4] {
            check_merge(&seq1, &seq2, 5, n_threads);
        }
    }

    // ── Intersection helpers ──────────────────────────────────────────────────

/// All distinct k-mers found in `seq`.
    fn extract_kmers(seq: &[u8], k: usize) -> std::collections::BTreeSet<Vec<u8>> {
        if seq.len() < k { return std::collections::BTreeSet::new(); }
        seq.windows(k).map(|w| w.to_vec()).collect()
    }

    /// Core correctness check for `intersect`:
    ///   1. Build individual SBWTs for seq1 and seq2.
    ///   2. Compute the intersection SBWT via `intersect()`.
    ///   3. Manually intersect the k-mer sets and build a ground-truth SBWT.
    ///   4. Assert that both yield the same set of real k-mers AND dummy nodes.
    ///
    /// Tested for both optimize_peak_ram variants of the interleaving.
    fn check_intersect(seq1: &[u8], seq2: &[u8], k: usize, n_threads: usize) {
        // Ground truth: manual k-mer intersection.
        let kmers1 = extract_kmers(seq1, k);
        let kmers2 = extract_kmers(seq2, k);
        let mut isec_kmers: Vec<Vec<u8>> = kmers1.intersection(&kmers2).cloned().collect();
        
        // Ensure deterministic order for ground-truth SBWT construction
        isec_kmers.sort();

        let (sbwt_true, _) = SbwtIndexBuilder::<BitPackedKmerSortingMem>::new()
            .k(k)
            .n_threads(n_threads)
            .run_from_vecs(&isec_kmers);

        // Build the two individual SBWTs.
        let (sbwt1, _) = SbwtIndexBuilder::<BitPackedKmerSortingMem>::new()
            .k(k)
            .n_threads(n_threads)
            .run_from_slices(&[seq1]);
        let (sbwt2, _) = SbwtIndexBuilder::<BitPackedKmerSortingMem>::new()
            .k(k)
            .n_threads(n_threads)
            .run_from_slices(&[seq2]);

        let spectrum_true_raw = sbwt_true.reconstruct_padded_spectrum(n_threads);
        let mut kmers_true: Vec<Vec<u8>> = spectrum_true_raw
            .chunks(k)
            .map(|c| c.to_vec())
            .collect();
        kmers_true.sort();

        // Test both high-RAM and low-RAM interleaving variants.
        for optimize_peak_ram in [false, true] {
            let interleaving = MergeInterleaving::new(&sbwt1, &sbwt2, optimize_peak_ram, n_threads);
            let sbwt_isec = intersect(
                Arc::new(sbwt1.clone()),
                Arc::new(sbwt2.clone()),
                Arc::new(interleaving),
                0,
                /* optimize_peak_ram */
                true,
                n_threads,
            );

            let spectrum_isec_raw = sbwt_isec.reconstruct_padded_spectrum(n_threads);
            let mut kmers_isec: Vec<Vec<u8>> = spectrum_isec_raw
                .chunks(k)
                .map(|c| c.to_vec())
                .collect();
            kmers_isec.sort();

            assert_eq!(
                kmers_true, kmers_isec,
                "k={k}, n_threads={n_threads}, optimize_peak_ram={optimize_peak_ram}: \n\
                 Intersection spectrum (including dummies) differs from ground truth.\n\
                 Expected (first 5): {:?}\n\
                 Found    (first 5): {:?}",
                 &kmers_true[..std::cmp::min(5, kmers_true.len())],
                 &kmers_isec[..std::cmp::min(5, kmers_isec.len())]
            );
        }
    }

    // ── Intersection tests ────────────────────────────────────────────────────

    #[test]
    fn test_intersect_no_dummy_repair() {
        // k=5.  Both sequences share the prefix "AAGTACGAAAAAAAAAAAA", so the shared k-mers
        // {AAGTA, AGTAC, GTACG, TACGA, ACGAA, CGAAA, GAAAA, AAAAA} are all in the intersection.
        // AAGTA was already a source node in both individual indexes (no "?AAGT" k-mer exists
        // in either), so the intersection inherits its dummy chain without introducing any new
        // source nodes — no dummy chain repair needed.
        let seq1 = b"AAGTACGAAAAAAAAAAACCCC";
        let seq2 = b"AAGTACGAAAAAAAAAAATTTT";
        for n_threads in [1, 3] {
            check_intersect(seq1, seq2, 5, n_threads);
        }
    }

    #[test]
    fn test_intersect_with_dummy_repair() {
        // k=5.  seq1 = "A" + shared, seq2 = "C" + shared.
        // Both contain {GTACG, TACGA, ACGAA, CGAAA, GAAAA, AAAAA, …} from the shared suffix.
        // GTACG's predecessor is AGTAC in seq1 and CGTAC in seq2; neither is in the
        // intersection, so GTACG becomes a NEW source node in the intersection.
        // An auxiliary index must be built and merged into both inputs to supply
        // the missing dummy predecessor chain.
        let shared = b"GTACGAAAAAAAAAAAAAAA" as &[u8];
        let seq1: Vec<u8> = b"A".iter().chain(shared).copied().collect();
        let seq2: Vec<u8> = b"C".iter().chain(shared).copied().collect();
        for n_threads in [1, 3] {
            check_intersect(&seq1, &seq2, 5, n_threads);
        }
    }

    #[test]
    fn test_intersect_random() {
        // Larger randomised pairs exercise the full algorithm across many k-mer overlap
        // patterns, including both direct and dummy-repair scenarios.
        let k = 5;
        for (seed1, seed2) in [(42, 1042), (1337, 2337), (9999, 10999)] {
            let seq1 = crate::util::gen_random_dna_string(1000, seed1);
            let seq2 = crate::util::gen_random_dna_string(1000, seed2);
            for n_threads in [1, 3] {
                check_intersect(&seq1, &seq2, k, n_threads);
            }
        }
    }

    #[test]
    fn test_intersect_with_dummy_repair_long() {
        // Same dummy-repair trigger as test_intersect_with_dummy_repair: seq1 = "A" + shared,
        // seq2 = "C" + shared, so the first shared k-mer (GTACG) becomes a new source
        // node in the intersection (its predecessors AGTAC / CGTAC are not shared).
        //
        // Both sequences have long but disjoint tails (all-A vs all-C) to make the
        // indices large and exercise multi-threaded behaviour without accidentally
        // introducing extra shared k-mers that would provide GTACG an incoming edge
        // and suppress dummy chain repair (random tails of length ≥ ~1000 saturate all 4^5
        // 5-mers, so both tails would share AGTAC and cancel the dummy-repair trigger).
        let shared = b"GTACGTTTTTTTTTTTTTT" as &[u8];
        let tail1 = vec![b'A'; 2000];
        let tail2 = vec![b'C'; 2000];
        let seq1: Vec<u8> = b"A".iter().chain(shared).chain(tail1.iter()).copied().collect();
        let seq2: Vec<u8> = b"C".iter().chain(shared).chain(tail2.iter()).copied().collect();
        for n_threads in [1, 2, 4] {
            check_intersect(&seq1, &seq2, 5, n_threads);
        }
    }

    #[test]
    fn test_intersect_no_overlap() {
        let seq1 = b"AAAAAAAAAA";
        let seq2 = b"CCCCCCCCCC";
        // Should result in an empty k-mer set, not a panic.
        check_intersect(seq1, seq2, 5, 1);
    }

    #[test]
    fn test_intersect_single_disjoint_kmers() {
        // index1 = {ACCTT}, index2 = {ATGCG}: each sequence is exactly one k-mer long
        // (k=5) and the two k-mers share nothing but the "A" prefix. The intersection
        // should be empty.
        let seq1 = b"ACCTT";
        let seq2 = b"ATGCG";
        check_intersect(seq1, seq2, 5, 1);
    }

    #[test]
    fn test_intersect_dead_end_dummy_chain() {
        // k=5.  Both sequences start with "ACCT" and diverge only at the 5th character, so
        // index1 and index2 share the whole padding chain $$$$A, $$$AC, $$ACC, $ACCT while
        // the real k-mer it was built for (ACCTT vs ACCTG) is not shared.  In the
        // intersection $ACCT loses its only outgoing edge, so the entire four-node chain
        // above it becomes purposeless and must be stripped (a BFS, not a single node).
        //
        // The tails share "GGTACGTTTT…", but GGTAC's predecessor differs (GGGTA vs AGGTA),
        // so GGTAC is also a new source node — dead-end stripping and dummy repair run
        // together in the same pass.
        let seq1 = b"ACCTTGGGTACGTTTTTTTTTTTT";
        let seq2 = b"ACCTGAGGTACGTTTTTTTTTTTT";
        for n_threads in [1, 2, 4] {
            check_intersect(seq1, seq2, 5, n_threads);
        }
    }

    #[test]
    fn test_intersect_dead_end_dummy_chain_long() {
        // Same dead-end chain trigger as test_intersect_dead_end_dummy_chain, but with long
        // disjoint tails so the indexes are big enough for the multi-threaded piece splitting
        // in the three-way pass and in the strip pass.  The tails are all-A / all-C rather
        // than random so they cannot accidentally supply the missing predecessors.
        let head1 = b"ACCTTGGGTACGTTTTTTTTTTTT" as &[u8];
        let head2 = b"ACCTGAGGTACGTTTTTTTTTTTT" as &[u8];
        let seq1: Vec<u8> = head1.iter().chain(vec![b'A'; 2000].iter()).copied().collect();
        let seq2: Vec<u8> = head2.iter().chain(vec![b'C'; 2000].iter()).copied().collect();
        for n_threads in [1, 2, 4] {
            check_intersect(&seq1, &seq2, 5, n_threads);
        }
    }

    #[test]
    fn test_intersect_orphaned_branch_dummy() {
        // k=4. index1 = {ACAA, ACGT}, index2 = {ACAC, ACGT} (two sequences each, all four
        // k-mers are sources in their own index). The shared dummy $ACA is a dead end
        // (continuations A vs C, nothing survives), so its ancestors $$AC and $$$A are
        // marked dead unconditionally. But $$AC also feeds the live $ACG chain of the
        // shared k-mer ACGT — killing it orphans $ACG. The orphan's incoming-edge bit is
        // cleared, it is collected as a source node, and the auxiliary index rebuilds the
        // chain $$$A, $$AC, $ACG, re-creating exactly the killed ancestors still needed.
        // Expected result: padded spectrum of {ACGT}.
        let k = 4;
        let mut isec_kmers: Vec<Vec<u8>> = vec![b"ACGT".to_vec()];
        isec_kmers.sort();
        let (sbwt_true, _) = SbwtIndexBuilder::<BitPackedKmerSortingMem>::new()
            .k(k).run_from_vecs(&isec_kmers);
        let spectrum_true_raw = sbwt_true.reconstruct_padded_spectrum(1);
        let mut kmers_true: Vec<Vec<u8>> = spectrum_true_raw.chunks(k).map(|c| c.to_vec()).collect();
        kmers_true.sort();

        let (sbwt1, _) = SbwtIndexBuilder::<BitPackedKmerSortingMem>::new()
            .k(k).run_from_slices(&[b"ACAA" as &[u8], b"ACGT"]);
        let (sbwt2, _) = SbwtIndexBuilder::<BitPackedKmerSortingMem>::new()
            .k(k).run_from_slices(&[b"ACAC" as &[u8], b"ACGT"]);

        for n_threads in [1, 2] {
            for optimize_peak_ram in [false, true] {
                let interleaving = MergeInterleaving::new(&sbwt1, &sbwt2, optimize_peak_ram, n_threads);
                let sbwt_isec = intersect(
                    Arc::new(sbwt1.clone()), Arc::new(sbwt2.clone()),
                    Arc::new(interleaving), 0, true, n_threads,
                );
                let spectrum_isec_raw = sbwt_isec.reconstruct_padded_spectrum(n_threads);
                let mut kmers_isec: Vec<Vec<u8>> = spectrum_isec_raw.chunks(k).map(|c| c.to_vec()).collect();
                kmers_isec.sort();
                assert_eq!(kmers_true, kmers_isec,
                    "n_threads={n_threads}, optimize_peak_ram={optimize_peak_ram}: \
                     orphaned branch dummy not repaired correctly");
            }
        }
    }

    #[test]
    fn test_intersect_revived_dead_end_dummy() {
        // k=4. Intersection = {ACGT}, a new source (predecessors TACG / GACG not shared).
        // Both indexes contain dummy $$AC (via sources ACAA / ACCC) whose edges (A vs C)
        // AND to nothing -> pass 1 flags it dead. But ACGT's fresh chain $$$A,$$AC,$ACG
        // passes through $$AC, so the aux index must revive it.
        let seq1 = b"ACAATACGT";
        let seq2 = b"ACCCGACGT";
        check_intersect(seq1, seq2, 4, 1);
    }

    #[test]
    fn test_difference_single_disjoint_kmers() {
        let seq1 = b"ACCTT";
        let seq2 = b"ATGCG";
        check_difference(seq1, seq2, 5, 1);
    }

    // ── Difference helpers and tests ──────────────────────────────────────────

    fn check_difference(seq1: &[u8], seq2: &[u8], k: usize, n_threads: usize) {
        // Ground truth: manual k-mer difference (kmers in seq1 but not in seq2).
        let kmers1 = extract_kmers(seq1, k);
        let kmers2 = extract_kmers(seq2, k);
        let mut diff_kmers: Vec<Vec<u8>> = kmers1.difference(&kmers2).cloned().collect();
        diff_kmers.sort();

        let (sbwt_true, _) = SbwtIndexBuilder::<BitPackedKmerSortingMem>::new()
            .k(k)
            .n_threads(n_threads)
            .run_from_vecs(&diff_kmers);

        let (sbwt1, _) = SbwtIndexBuilder::<BitPackedKmerSortingMem>::new()
            .k(k)
            .n_threads(n_threads)
            .run_from_slices(&[seq1]);
        let (sbwt2, _) = SbwtIndexBuilder::<BitPackedKmerSortingMem>::new()
            .k(k)
            .n_threads(n_threads)
            .run_from_slices(&[seq2]);

        let spectrum_true_raw = sbwt_true.reconstruct_padded_spectrum(n_threads);
        let mut kmers_true: Vec<Vec<u8>> = spectrum_true_raw.chunks(k).map(|c| c.to_vec()).collect();
        kmers_true.sort();

        for optimize_peak_ram in [false, true] {
            let interleaving = MergeInterleaving::new(&sbwt1, &sbwt2, optimize_peak_ram, n_threads);

            assert_eq!(
                interleaving.difference_size(),
                diff_kmers.len(),
                "difference_size() mismatch: k={k}, n_threads={n_threads}"
            );

            let sbwt_diff = difference(
                Arc::new(sbwt1.clone()),
                Arc::new(sbwt2.clone()),
                Arc::new(interleaving),
                0,
                true,
                n_threads,
            );

            assert_eq!(sbwt_diff.n_kmers(), diff_kmers.len(),
                "n_kmers mismatch: k={k}, n_threads={n_threads}, optimize_peak_ram={optimize_peak_ram}");

            // When the difference is empty the SBWT has no real k-mers; spectrum
            // reconstruction is meaningless for a 0-node index so skip it.
            if diff_kmers.is_empty() { continue; }

            let spectrum_diff_raw = sbwt_diff.reconstruct_padded_spectrum(n_threads);
            let mut kmers_diff: Vec<Vec<u8>> = spectrum_diff_raw.chunks(k).map(|c| c.to_vec()).collect();
            kmers_diff.sort();

            assert_eq!(
                kmers_true, kmers_diff,
                "k={k}, n_threads={n_threads}, optimize_peak_ram={optimize_peak_ram}: \n\
                 Difference spectrum differs from ground truth.\n\
                 Expected (first 5): {:?}\n\
                 Found    (first 5): {:?}",
                &kmers_true[..std::cmp::min(5, kmers_true.len())],
                &kmers_diff[..std::cmp::min(5, kmers_diff.len())]
            );
        }
    }

    #[test]
    fn test_difference_dead_dummy_chain() {
        // k=5, index1 = {ACCTT}, index2 = S_5(CACCTT) = {CACCT, ACCTT}.  The difference is
        // empty, so the only correct output is the padded spectrum of the empty set, {$$$$$}.
        //
        // But index1's dummy chain $$$$A, $$$AC, $$ACC, $ACCT is exclusive to index1 (index2
        // pads CACCT instead, giving $$$$C, $$$CA, $$CAC, $CACC), so the s1 && !s2 predicate
        // keeps all four rows even though the real k-mer they were introduced to support is
        // gone.  The result is a dangling dummy chain leading nowhere.
        let seq1 = b"ACCTT" as &[u8];
        let seq2 = b"CACCTT" as &[u8];
        let n_threads = 1;

        let (sbwt1, _) = SbwtIndexBuilder::<BitPackedKmerSortingMem>::new()
            .k(5)
            .n_threads(n_threads)
            .run_from_slices(&[seq1]);
        let (sbwt2, _) = SbwtIndexBuilder::<BitPackedKmerSortingMem>::new()
            .k(5)
            .n_threads(n_threads)
            .run_from_slices(&[seq2]);

        for optimize_peak_ram in [false, true] {
            let interleaving = MergeInterleaving::new(&sbwt1, &sbwt2, optimize_peak_ram, n_threads);
            let sbwt_diff = difference(
                Arc::new(sbwt1.clone()),
                Arc::new(sbwt2.clone()),
                Arc::new(interleaving),
                0,
                true,
                n_threads,
            );

            assert_eq!(sbwt_diff.n_kmers(), 0,
                "optimize_peak_ram={optimize_peak_ram}");
            assert_eq!(
                sbwt_diff.n_sets(), 1,
                "optimize_peak_ram={optimize_peak_ram}: expected only the root row $$$$$, \
                 got {} rows (dangling dummy chain retained)",
                sbwt_diff.n_sets()
            );
        }
    }

    #[test]
    fn test_difference_with_dummy_repair() {
        // seq1 and seq2 share a common prefix; unique tails ensure seq1 has k-mers not in seq2.
        let seq1 = b"AAGTACGAAAAAAAAAAACCCC";
        let seq2 = b"AAGTACGAAAAAAAAAAATTTT";
        for n_threads in [1, 3] {
            check_difference(seq1, seq2, 5, n_threads);
        }
    }

    #[test]
    fn test_difference_no_dummy_repair() {
        // Same no-dummy-repair trigger as test_intersect_no_dummy_repair but for difference.
        // seq1 = "A" + shared, seq2 = "C" + shared.
        // The first k-mer of the shared block (GTACG) is in seq1 but NOT in seq2's unique part,
        // however its predecessor AGTAC is NOT in the difference — no dummy chain repair needed.
        let shared = b"GTACGAAAAAAAAAAAAAAA" as &[u8];
        let seq1: Vec<u8> = b"A".iter().chain(shared).copied().collect();
        let seq2: Vec<u8> = b"C".iter().chain(shared).copied().collect();
        for n_threads in [1, 3] {
            check_difference(&seq1, &seq2, 5, n_threads);
        }
    }

    #[test]
    fn test_difference_random() {
        let k = 5;
        for (seed1, seed2) in [(42, 1042), (1337, 2337), (9999, 10999)] {
            let seq1 = crate::util::gen_random_dna_string(1000, seed1);
            let seq2 = crate::util::gen_random_dna_string(1000, seed2);
            for n_threads in [1, 3] {
                check_difference(&seq1, &seq2, k, n_threads);
            }
        }
    }

    #[test]
    fn test_difference_no_overlap() {
        // When seq1 and seq2 are completely disjoint, difference == seq1's k-mers.
        let seq1 = b"AAAATTTTAA";
        let seq2 = b"CCCCGGGGCC";
        check_difference(seq1, seq2, 5, 1);
    }

    #[test]
    fn test_difference_full_overlap() {
        // When seq1 == seq2, the difference should be empty.
        let seq = b"AAGTACGAAAAAAAAAAACCCC";
        check_difference(seq, seq, 5, 1);
    }

    #[test]
    fn test_difference_random_wide_k() {
        // test_difference_random uses k=5 on 1000-base strings, which saturates all 4^5 5-mers
        // and so barely exercises dummy chains at all.  Larger k leaves the k-mer space sparse,
        // giving long dummy chains in both inputs, chains that index2 removes the real k-mers
        // from, and chains index1 has no counterpart for.
        for k in [3, 5, 8, 12, 20] {
            for (seed1, seed2) in [(1, 2), (31, 32)] {
                for (len1, len2) in [(60, 200), (500, 500)] {
                    let a = crate::util::gen_random_dna_string(len1, seed1);
                    let b = crate::util::gen_random_dna_string(len2, seed2);
                    // Partial overlap: b prefixed with the first half of a.
                    let shared: Vec<u8> = a[..len1 / 2].iter().chain(b.iter()).copied().collect();
                    for (x, y) in [(&a, &b), (&a, &shared), (&shared, &a)] {
                        check_difference(x, y, k, 1);
                    }
                }
            }
        }
    }

    #[test]
    fn test_difference_thread_and_low_ram_consistency() {
        // Thread count and low-RAM mode must not change the result.  Compared directly rather
        // than via check_difference, which reconstructs the ground-truth spectrum with
        // n_threads and panics when the difference is empty (a 1-row index cannot be split).
        for k in [5, 12] {
            let a = crate::util::gen_random_dna_string(2000, 1);
            let b = crate::util::gen_random_dna_string(2000, 2);
            let shared: Vec<u8> = a[..1000].iter().chain(b.iter()).copied().collect();
            for (x, y) in [(&a, &b), (&a, &shared), (&shared, &a)] {
                let (s1, _) = SbwtIndexBuilder::<BitPackedKmerSortingMem>::new()
                    .k(k).n_threads(1).run_from_slices(&[x.as_slice()]);
                let (s2, _) = SbwtIndexBuilder::<BitPackedKmerSortingMem>::new()
                    .k(k).n_threads(1).run_from_slices(&[y.as_slice()]);
                let mut reference = None;
                for n_threads in [1, 2, 4] {
                    for opt_ram in [false, true] {
                        let il = MergeInterleaving::new(&s1, &s2, opt_ram, n_threads);
                        let d = difference(Arc::new(s1.clone()), Arc::new(s2.clone()),
                            Arc::new(il), 0, opt_ram, n_threads);
                        let rows: Vec<Vec<bool>> = (0..d.n_sets()).map(|i| {
                            (0..4).map(|c| d.sbwt().set_contains(i, c)).collect()
                        }).collect();
                        let got = (d.n_kmers(), d.n_sets(), rows);
                        match &reference {
                            None => reference = Some(got),
                            Some(r) => assert_eq!(*r, got,
                                "k={k} n_threads={n_threads} opt_ram={opt_ram}"),
                        }
                    }
                }
            }
        }
    }


    // ── Difference stress-test helpers ────────────────────────────────────────
    //
    // These build indexes from explicit k-mer *sets* rather than from sequences.  Every SBWT is
    // the padded spectrum of some k-mer set, so this reaches dummy-chain configurations that are
    // awkward or impossible to reach by picking input sequences, which is where the difference's
    // dummy handling has been wrong.

    fn build_from_kmer_set(kmers: &[Vec<u8>], k: usize) -> SbwtIndex<SubsetMatrix> {
        let mut sorted = kmers.to_vec();
        sorted.sort();
        SbwtIndexBuilder::<BitPackedKmerSortingMem>::new()
            .k(k).n_threads(1).run_from_vecs(&sorted).0
    }

    /// Full identifying content of an index: k-mer count, row count, and every subset.  Two
    /// indexes with equal fingerprints are equal, since the C array is derived from the rows.
    /// Stronger than comparing reconstructed spectra, and it works for empty results, where
    /// reconstruct_padded_spectrum cannot be used.
    fn index_fingerprint(d: &SbwtIndex<SubsetMatrix>) -> (usize, usize, Vec<Vec<bool>>) {
        let rows = (0..d.n_sets())
            .map(|i| (0..4).map(|c| d.sbwt().set_contains(i, c)).collect())
            .collect();
        (d.n_kmers(), d.n_sets(), rows)
    }

    /// Asserts `difference(build(kmers1), build(kmers2))` equals a from-scratch build of the set
    /// difference, in both the default and low-RAM modes.  `case` labels the failure.
    fn assert_difference_correct(kmers1: &[Vec<u8>], kmers2: &[Vec<u8>], k: usize, case: &str) {
        // These are k-mer sets: drop repeats so the expected count is a set difference.
        let dedup = |v: &[Vec<u8>]| {
            let mut out: Vec<Vec<u8>> = Vec::new();
            for x in v { if !out.contains(x) { out.push(x.clone()); } }
            out
        };
        let kmers1 = &dedup(kmers1);
        let kmers2 = &dedup(kmers2);
        let expected: Vec<Vec<u8>> =
            kmers1.iter().filter(|x| !kmers2.contains(x)).cloned().collect();
        let truth = index_fingerprint(&build_from_kmer_set(&expected, k));

        for optimize_peak_ram in [false, true] {
            let a = build_from_kmer_set(kmers1, k);
            let b = build_from_kmer_set(kmers2, k);
            let il = MergeInterleaving::new(&a, &b, optimize_peak_ram, 1);
            assert_eq!(il.difference_size(), expected.len(),
                "{case}: difference_size() mismatch (k={k}, optimize_peak_ram={optimize_peak_ram})");
            let d = difference(Arc::new(a), Arc::new(b), Arc::new(il), 0, optimize_peak_ram, 1);
            assert_eq!(truth, index_fingerprint(&d),
                "{case}: differs from ground truth (k={k}, optimize_peak_ram={optimize_peak_ram})\n\
                 index1 = {kmers1:?}\nindex2 = {kmers2:?}\nexpected k-mers = {expected:?}");
        }
    }

    /// All k-mers of `seq`, in the order they first occur.
    fn kmers_of(seq: &[u8], k: usize) -> Vec<Vec<u8>> {
        let mut out: Vec<Vec<u8>> = Vec::new();
        for w in seq.windows(k) {
            if !out.contains(&w.to_vec()) { out.push(w.to_vec()); }
        }
        out
    }

    #[test]
    fn test_difference_shared_dummy_prefix_sweep() {
        // The result may need a dummy that both inputs have.  A and B are single k-mers sharing
        // exactly a j-base prefix, so their dummy chains coincide for the first j+1 rows and
        // diverge after: $^k and $^(k-1)A... are shared, the rest are exclusive.  Sweeping j
        // walks the shared/exclusive boundary along the whole chain.  j=1, k=5 is the originally
        // reported failure (ACCTT vs ATGCG): $$$$A is shared, so differencing the dummies dropped
        // it and broke the chain to ACCTT.
        // k >= 3 is a requirement of the bitpacked k-mer sorting builder (BIN_PREFIX_LEN).
        for k in [3, 5, 8] {
            for j in 0..k {
                let shared_prefix: Vec<u8> = vec![b'A'; j];
                let a: Vec<u8> = shared_prefix.iter().copied()
                    .chain(std::iter::repeat(b'C').take(k - j)).collect();
                let b: Vec<u8> = shared_prefix.iter().copied()
                    .chain(std::iter::repeat(b'G').take(k - j)).collect();
                assert_difference_correct(&[a.clone()], &[b.clone()], k,
                    &format!("shared prefix of length {j}"));
                // And the reverse direction, which differs because only index1's dummies are
                // candidates for reuse.
                assert_difference_correct(&[b], &[a], k,
                    &format!("shared prefix of length {j}, reversed"));
            }
        }
    }

    #[test]
    fn test_difference_dead_dummy_chain_sweep() {
        // The mirror case: index1's chain is exclusive to index1, but the real k-mer it was
        // introduced to support is removed by index2, leaving the chain dangling.  index2 holds
        // that k-mer plus a predecessor for it, so index2 pads the predecessor instead and has no
        // chain in common with index1.  Sweeping the shared suffix length varies how much of
        // index1's chain is exclusive.
        for k in [3, 5, 8] {
            for tail in 1..k {
                // x = A^(k-tail) C^tail, and index2 additionally holds G + x[..k-1], a
                // predecessor of x, so x needs no chain in index2.
                let x: Vec<u8> = std::iter::repeat(b'A').take(k - tail)
                    .chain(std::iter::repeat(b'C').take(tail)).collect();
                let pred: Vec<u8> = std::iter::once(b'G').chain(x[..k - 1].iter().copied()).collect();
                assert_difference_correct(&[x.clone()], &[x.clone(), pred], k,
                    &format!("dead chain, tail={tail}"));
            }
        }
    }

    #[test]
    fn test_difference_chain_broken_at_each_position() {
        // A path of real k-mers Y1 -> Y2 -> ... -> Ym in index1, with index2 removing exactly one
        // of them.  Every k-mer after the removed one loses its predecessor and needs a fresh
        // dummy chain, while the ones before it keep index1's.  Sweeping the removed position
        // covers the boundary between reused and rebuilt chains at every point in the path.
        for k in [3, 5, 8] {
            let seq = b"ACGTACGGTTACAGGCATTGCA";
            let path = kmers_of(seq, k);
            for i in 0..path.len() {
                let removed = vec![path[i].clone()];
                assert_difference_correct(&path, &removed, k,
                    &format!("path k-mer {i} of {} removed", path.len()));
            }
            // And removing every other one, so many chains are rebuilt at once.
            let every_other: Vec<Vec<u8>> =
                path.iter().step_by(2).cloned().collect();
            assert_difference_correct(&path, &every_other, k, "every other path k-mer removed");
        }
    }

    #[test]
    fn test_difference_branching_suffix_groups() {
        // Suffix groups holding several real k-mers, with index2 removing some but not all.  The
        // outgoing edges of a group are stored once, at its first row, and are computed by OR-ing
        // over the group's rows, so a group that is partly removed exercises that per-group
        // reduction and the "does this group hold a real difference k-mer" test that selects
        // between the index1 and auxiliary-index edge rules.
        for k in [3, 5] {
            // All four k-mers cA^(k-1) share the (k-1)-suffix A^(k-1), so they form one group.
            let group: Vec<Vec<u8>> = [b'A', b'C', b'G', b'T'].iter()
                .map(|&c| std::iter::once(c).chain(std::iter::repeat(b'A').take(k - 1)).collect())
                .collect();
            for keep_mask in 0..16u32 {
                let removed: Vec<Vec<u8>> = (0..4)
                    .filter(|i| (keep_mask >> i) & 1 == 1)
                    .map(|i| group[i].clone())
                    .collect();
                assert_difference_correct(&group, &removed, k,
                    &format!("branching group, removed mask {keep_mask:04b}"));
            }
        }
    }

    #[test]
    fn test_difference_phantom_edge_guard() {
        // Two chains leaving a common dummy prefix: index1 holds ACCTT and AGGGG, so its chains
        // $$$$A -> $$$AC -> ... and $$$$A -> $$$AG -> ... branch at $$$$A.  index2 removes AGGGG
        // and supplies a real predecessor for it (TAGGG), so index2 has no $$$AG chain of its own.
        //
        // The difference is {ACCTT}, whose chain must keep the C branch and drop the G branch.
        // Selecting edges by "in index1 but not index2" would set both, because index2 has
        // nothing in the $$$A group to cancel the G — a phantom edge into a row the result does
        // not contain.  The G branch must come from the auxiliary index or not at all.
        let a = vec![b"ACCTT".to_vec(), b"AGGGG".to_vec()];
        let b = vec![b"AGGGG".to_vec(), b"TAGGG".to_vec()];
        assert_difference_correct(&a, &b, 5, "phantom edge, C branch survives");

        // Same shape with the surviving branch on the other side.
        let a2 = vec![b"ACCTT".to_vec(), b"AGGGG".to_vec()];
        let b2 = vec![b"ACCTT".to_vec(), b"TACCT".to_vec()];
        assert_difference_correct(&a2, &b2, 5, "phantom edge, G branch survives");
    }

    #[test]
    fn test_difference_cyclic_results() {
        // Results whose de Bruijn graph is all cycles need no dummy chains beyond the root, which
        // is the only case that still takes the direct pass.  A self-loop (A^k) and longer cycles
        // from repeated periodic sequences both land there.
        for k in [3, 5, 8] {
            let self_loop = vec![vec![b'A'; k]];
            assert_difference_correct(&self_loop, &[vec![b'C'; k]], k, "self-loop preserved");
            assert_difference_correct(&self_loop, &self_loop, k, "self-loop removed");

            // ACGACGACG... is periodic, so its k-mers form one cycle with no source node.
            let periodic: Vec<u8> = b"ACG".iter().cycle().take(k + 32).copied().collect();
            let cycle = kmers_of(&periodic, k);
            assert_difference_correct(&cycle, &[vec![b'T'; k]], k, "cycle preserved");
            assert_difference_correct(&cycle, &cycle, k, "cycle removed");
            // A cycle plus a tail hanging off it: the tail needs a chain, the cycle does not.
            let with_tail: Vec<u8> = periodic.iter().copied().chain(b"TTTT".iter().copied()).collect();
            assert_difference_correct(&kmers_of(&with_tail, k), &cycle, k, "tail off a cycle");
        }
    }

    #[test]
    fn test_difference_degenerate_set_relations() {
        // Empty inputs, and the four containment regimes, at several k.
        for k in [3, 5, 8] {
            let a = kmers_of(b"ACGTTGCAACGT", k);
            let subset: Vec<Vec<u8>> = a.iter().take(a.len() / 2).cloned().collect();
            let disjoint = kmers_of(b"GGGGGGGGGGGG", k);
            let empty: Vec<Vec<u8>> = vec![];

            assert_difference_correct(&a, &a, k, "identical");
            assert_difference_correct(&a, &subset, k, "superset minus subset");
            assert_difference_correct(&subset, &a, k, "subset minus superset");
            assert_difference_correct(&a, &disjoint, k, "disjoint");
            assert_difference_correct(&a, &empty, k, "minus empty");
            assert_difference_correct(&empty, &a, k, "empty minus");
            assert_difference_correct(&empty, &empty, k, "empty minus empty");
        }
    }

    #[test]
    fn test_difference_random_kmer_sets() {
        // Random k-mer sets rather than random sequences.  A sequence's k-mers form a path, so
        // its padded spectrum has few dummy chains; an arbitrary set is fragmented into many
        // short unitigs and so has many, which is what the dummy handling operates on.  The
        // density sweep moves between "mostly isolated k-mers, one chain each" and "mostly
        // connected, few chains".
        for k in [3, 5, 8, 12] {
            for seed in [1u64, 7, 99] {
                for density_num in [1usize, 2, 4, 8] {
                    // Draw a pool of random k-mers by chunking a random string, then split it
                    // between the two inputs with an overlapping middle section.
                    let pool_len = 64 * k;
                    let pool_seq = crate::util::gen_random_dna_string(pool_len, seed);
                    let pool: Vec<Vec<u8>> = pool_seq.chunks_exact(k)
                        .map(|c| c.to_vec())
                        .step_by(density_num.max(1))
                        .collect();
                    if pool.len() < 4 { continue; }

                    let cut1 = pool.len() / 3;
                    let cut2 = 2 * pool.len() / 3;
                    let s1: Vec<Vec<u8>> = pool[..cut2].to_vec();
                    let s2: Vec<Vec<u8>> = pool[cut1..].to_vec();
                    assert_difference_correct(&s1, &s2, k,
                        &format!("random sets, seed={seed}, step={density_num}"));
                    assert_difference_correct(&s2, &s1, k,
                        &format!("random sets reversed, seed={seed}, step={density_num}"));
                }
            }
        }
    }

    #[test]
    fn test_difference_multi_sequence_inputs() {
        // Indexes built from several sequences at once, so each input has many separate dummy
        // chains, and the two inputs share some sequences but not others.
        let seqs: Vec<Vec<u8>> = (0..6)
            .map(|i| crate::util::gen_random_dna_string(40 + 7 * i, 1000 + i as u64))
            .collect();
        for k in [3, 5, 8, 12] {
            let mut a: Vec<Vec<u8>> = Vec::new();
            for s in &seqs[..4] { for x in kmers_of(s, k) { if !a.contains(&x) { a.push(x); } } }
            let mut b: Vec<Vec<u8>> = Vec::new();
            for s in &seqs[2..] { for x in kmers_of(s, k) { if !b.contains(&x) { b.push(x); } } }
            assert_difference_correct(&a, &b, k, "multi-sequence inputs");
            assert_difference_correct(&b, &a, k, "multi-sequence inputs reversed");
        }
    }
}

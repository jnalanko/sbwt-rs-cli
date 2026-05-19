pub mod interleaving;
mod common;
mod merge;
mod intersect_difference;

pub use interleaving::MergeInterleaving;
pub use merge::merge;
pub use intersect_difference::{intersect, difference};

#[cfg(test)]
mod tests {

    use crate::{BitPackedKmerSortingDisk, BitPackedKmerSortingMem, SbwtIndexBuilder};

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
}

//! `check` and `check-set-operation`: verifying that an SBWT index has exactly the k-mer set
//! it's supposed to have, without ever dumping an index's own (possibly huge) k-mer set as
//! text - only ever searching for k-mers we already know about from a bounded-size source
//! (a sequence file, or another index queried via a bounded-size sequence file).

use jseqio::reader::DynamicFastXReader;
use sbwt::sbwt_index_variant::SbwtIndexVariant;
use sbwt::*;

use crate::{default_lcs_path, load_index, load_lcs_if_exists};

/// Verification logic behind the `check` subcommand; see `main.rs::check_command` for the
/// clap argument parsing that calls this.
pub fn run_check(
    index_path: &std::path::Path,
    cpp_format: bool,
    invariant_check: bool,
    original_sequences: Option<&std::path::PathBuf>,
    add_revcomp: bool,
) {
    log::info!("Loading index from {:?}", index_path);
    let mut index = load_index(index_path, cpp_format);
    let lcs = load_lcs_if_exists(&default_lcs_path(index_path), &format!("LCS file not found at {} -> running without it", default_lcs_path(index_path).display()));

    let n = index.n_sets();
    let sigma = index.alphabet().len();
    log::info!("Index loaded: {} sets, {} k-mers, k={}, sigma={}", n, index.n_kmers(), index.k(), sigma);

    if invariant_check {
        log::info!("Checking SBWT structural invariant: sum_c rank(c, n_sets) + 1 == n_sets");
        let in_edges: usize = (0..sigma).map(|c| index.rank(c as u8, n)).sum();
        if in_edges + 1 == n {
            log::info!("OK: {} in-edge(s) + root == {} sets", in_edges, n);
            println!("PASS");
        } else {
            log::error!(
                "FAIL: invariant violated: {} in-edge(s) + root != {} sets; {} node(s) have no in-coming edge",
                in_edges, n, n - 1 - in_edges
            );
            println!("FAIL");
            std::process::exit(1);
        }
    }

    let streaming_index = lcs.as_ref().map(|lcs| {
        index.get_streaming_index(&lcs)
    });

    if let Some(original_sequences) = original_sequences {
        // Open reader early to fail early if there is a problem
        let reader = sbwt::FastXReader::new(original_sequences).unwrap();
        let mut reader = sbwt::SeqStreamWithPossiblyRevComp::new(reader, add_revcomp);

        log::info!("Marking dummy nodes");
        let mut visited_marks = index.compute_dummy_node_marks();
        log::info!("Searching all sequences in {}", original_sequences.display());
        let k = index.k();
        while let Some(seq) = reader.stream_next() {
            if let Some(streaming_index) = &streaming_index {
                let mut prev_non_acgt_char = -1_isize;
                for (pos, (len, colex)) in streaming_index.matching_statistics_iter(seq).enumerate() {
                    if !sbwt::is_dna(seq[pos]) {
                        prev_non_acgt_char = pos as isize;
                    }
                    let dna_run_len = (pos as isize - prev_non_acgt_char) as usize;
                    if dna_run_len >= k {
                        // Found in input data -> must be found in the SBWT
                        assert!(len == k); // is found in the SBWT
                        assert!(colex.len() == 1); // k-mer should have colex range of length 1

                        // Mark as found
                        visited_marks.set(colex.start, true);
                    }
                }
            } else {
                for kmer in seq.windows(k) {
                    if kmer.iter().all(|&c| is_dna(c)) {
                        // Found in input data -> must be found in the SBWT
                        let colex = index.search(kmer).expect(&format!("k-mer {} not found in the SBWT", String::from_utf8_lossy(kmer)));
                        assert!(colex.len() == 1); // k-mer should have colex range of length 1

                        // Mark as found
                        visited_marks.set(colex.start, true);

                    }
                }
            }
        }

        // Check that all the k-mers in the SBWT were in the input sequences

        if visited_marks.count_ones() != visited_marks.len() {
            log::info!("SBWT has a k-mer that was not in the input sequences. Fetching the k-mer...");
            index.build_select();
            let kmer = index.access_kmer(visited_marks.first_zero().unwrap());
            panic!("SBWT has k-mer {}, which is not in the input sequences.", String::from_utf8_lossy(&kmer));
        }
        log::info!("OK. The SBWT has exactly the same k-mers as the input sequences.");
    }

}

/// Calls `visit(kmer)` for every ACGT-only k-mer window across every record in the
/// fasta/fastq file at `path`.
fn for_each_valid_kmer(path: &std::path::Path, k: usize, visit: &mut impl FnMut(&[u8])) {
    let mut reader = DynamicFastXReader::from_file(&path).unwrap();
    while let Some(rec) = reader.read_next().unwrap() {
        for kmer in rec.seq.windows(k) {
            if kmer.iter().all(|&c| is_dna(c)) {
                visit(kmer);
            }
        }
    }
}

/// Calls `expect_in_result(kmer)` for every ACGT-only k-mer window of `seq1_path` that is also
/// present in `sbwt2` (used to verify `result = intersect(sbwt1, sbwt2)`). Uses `streaming2`
/// for a single streaming pass per record when available (much faster than searching from
/// scratch for every k-mer), or falls back to a plain per-k-mer search against `sbwt2`.
fn for_each_kmer_in_seq1_and_sbwt2(
    seq1_path: &std::path::Path,
    k: usize,
    sbwt2: &SbwtIndexVariant,
    streaming2: Option<&StreamingIndex<SbwtIndexVariant, LcsArray>>,
    expect_in_result: &mut impl FnMut(&[u8]),
) {
    let mut reader = DynamicFastXReader::from_file(&seq1_path).unwrap();
    while let Some(rec) = reader.read_next().unwrap() {
        let seq = rec.seq;
        if let Some(streaming2) = streaming2 {
            let mut prev_non_acgt_char = -1_isize;
            for (pos, (len, _)) in streaming2.matching_statistics_iter(seq).enumerate() {
                if !is_dna(seq[pos]) {
                    prev_non_acgt_char = pos as isize;
                }
                let dna_run_len = (pos as isize - prev_non_acgt_char) as usize;
                if dna_run_len >= k && len == k {
                    // Found in sbwt2 -> also expected in the intersection
                    expect_in_result(&seq[pos + 1 - k..=pos]);
                }
            }
        } else {
            for kmer in seq.windows(k) {
                if kmer.iter().all(|&c| is_dna(c)) && sbwt2.search(kmer).is_some() {
                    expect_in_result(kmer);
                }
            }
        }
    }
}

/// Calls `expect_in_result(kmer)` for every ACGT-only k-mer window of `seq1_path` that is NOT
/// present in `sbwt2` (used to verify `result = difference(sbwt1, sbwt2)`, i.e. seq1 \ seq2).
/// Uses `streaming2` for a single streaming pass per record when available, or falls back to a
/// plain per-k-mer search against `sbwt2`.
fn for_each_kmer_in_seq1_not_sbwt2(
    seq1_path: &std::path::Path,
    k: usize,
    sbwt2: &SbwtIndexVariant,
    streaming2: Option<&StreamingIndex<SbwtIndexVariant, LcsArray>>,
    expect_in_result: &mut impl FnMut(&[u8]),
) {
    let mut reader = DynamicFastXReader::from_file(&seq1_path).unwrap();
    while let Some(rec) = reader.read_next().unwrap() {
        let seq = rec.seq;
        if let Some(streaming2) = streaming2 {
            let mut prev_non_acgt_char = -1_isize;
            for (pos, (len, _)) in streaming2.matching_statistics_iter(seq).enumerate() {
                if !is_dna(seq[pos]) {
                    prev_non_acgt_char = pos as isize;
                }
                let dna_run_len = (pos as isize - prev_non_acgt_char) as usize;
                if dna_run_len >= k && len != k {
                    // Not found in sbwt2 -> expected in the difference
                    expect_in_result(&seq[pos + 1 - k..=pos]);
                }
            }
        } else {
            for kmer in seq.windows(k) {
                if kmer.iter().all(|&c| is_dna(c)) && sbwt2.search(kmer).is_none() {
                    expect_in_result(kmer);
                }
            }
        }
    }
}

/// Verification logic behind the `check-set-operation` subcommand; see
/// `main.rs::check_set_operation_command` for the clap argument parsing that calls this.
///
/// Verifies that `result` (claimed to be `sbwt1 <op> sbwt2`) contains exactly the expected
/// k-mer set, computed from `seq1`/`seq2` (the original sequences `sbwt1`/`sbwt2` were built
/// from) rather than by enumerating `result`'s own k-mers, since `result` may be far too large
/// to dump as text. Uses streaming search (via the `.lcs` file next to `sbwt2`, if present)
/// when checking k-mer membership in `sbwt2`, falling back to plain per-k-mer search otherwise.
pub fn run_check_set_operation(
    op: &str,
    seq1_path: &std::path::Path,
    seq2_path: &std::path::Path,
    sbwt1_path: &std::path::Path,
    sbwt2_path: &std::path::Path,
    result_path: &std::path::Path,
    cpp_format: bool,
) {
    log::info!("Loading sbwt1 from {}", sbwt1_path.display());
    let index1 = load_index(sbwt1_path, cpp_format);
    log::info!("Loading sbwt2 from {}", sbwt2_path.display());
    let index2 = load_index(sbwt2_path, cpp_format);
    log::info!("Loading result index from {}", result_path.display());
    let mut result = load_index(result_path, cpp_format);

    let k = index1.k();
    assert_eq!(k, index2.k(), "sbwt1 and sbwt2 must have the same k");
    assert_eq!(k, result.k(), "the result index must have the same k as sbwt1/sbwt2");

    let lcs2 = load_lcs_if_exists(
        &default_lcs_path(sbwt2_path),
        "falling back to non-streaming search against sbwt2",
    );
    let streaming2 = lcs2.as_ref().map(|lcs| index2.get_streaming_index(lcs));

    // Pre-marks dummy nodes as visited, same as `check_command`: a node is only "unexpected"
    // at the end if it's a real (non-dummy) k-mer that never got marked below.
    let mut visited_marks = result.compute_dummy_node_marks();
    let mut expect_in_result = |kmer: &[u8]| {
        let colex = result.search(kmer).unwrap_or_else(|| {
            panic!(
                "expected k-mer {} to be in the result index (op={op}), but it was not found",
                String::from_utf8_lossy(kmer)
            )
        });
        assert_eq!(colex.len(), 1, "k-mer {} should have a colex range of length 1", String::from_utf8_lossy(kmer));
        visited_marks.set(colex.start, true);
    };

    match op {
        "merge" => {
            log::info!("Checking that every k-mer of seq1 and of seq2 is in the result (union)");
            for_each_valid_kmer(seq1_path, k, &mut expect_in_result);
            for_each_valid_kmer(seq2_path, k, &mut expect_in_result);
        }
        "intersect" => {
            log::info!("Checking that every k-mer of seq1 that's also in sbwt2 is in the result (intersection)");
            for_each_kmer_in_seq1_and_sbwt2(seq1_path, k, &index2, streaming2.as_ref(), &mut expect_in_result);
        }
        "difference" => {
            log::info!("Checking that every k-mer of seq1 that's not in sbwt2 is in the result (seq1 \\ seq2)");
            for_each_kmer_in_seq1_not_sbwt2(seq1_path, k, &index2, streaming2.as_ref(), &mut expect_in_result);
        }
        other => {
            eprintln!("error: unknown --op '{other}', expected one of: merge, intersect, difference");
            std::process::exit(2);
        }
    }

    if visited_marks.count_ones() != visited_marks.len() {
        log::info!("Result index has an unexpected k-mer. Fetching it...");
        result.build_select();
        let unvisited = visited_marks.first_zero().unwrap();
        let kmer = result.access_kmer(unvisited);
        panic!(
            "result index has k-mer {}, which should not be in the {op} of the given sequences",
            String::from_utf8_lossy(&kmer)
        );
    }

    log::info!("OK. The result index has exactly the expected k-mer set for {op}.");
    println!("PASS");
}

use std::io::BufWriter;
use std::path::{Path, PathBuf};

use jseqio::reader::DynamicFastXReader;
use sbwt::dbg::Dbg;
use sbwt::sbwt_index_variant::SbwtIndexVariant;
use sbwt::*;

use crate::{default_lcs_path, load_index, load_lcs_if_exists};

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

fn lookup_each_valid_kmer_nonstreaming_and_callback (
    seq1_path: &std::path::Path,
    index2: &SbwtIndexVariant,
    callback: &mut impl FnMut((Option<usize>, &[u8])), // Calls back on the colex rank of the k-mer and the k-mer string (None if not found)
) {

    let k = index2.k();
    let mut reader = DynamicFastXReader::from_file(&seq1_path).unwrap();
    while let Some(rec) = reader.read_next().unwrap() {
        for kmer in rec.seq.windows(k) {
            if kmer.iter().all(|&c| is_dna(c)) {
                callback((index2.search(kmer).map(|range| range.start), &kmer));
            }
        }
    }
}

fn lookup_each_valid_kmer_and_callback (
    seq1_path: &std::path::Path,
    index2: &StreamingIndex<SbwtIndexVariant, LcsArray>,
    callback: &mut impl FnMut((Option<usize>, &[u8])), // Calls back on the colex rank of the k-mer and the k-mer string (None if not found)
) {

    let k = index2.k;
    let mut reader = DynamicFastXReader::from_file(&seq1_path).unwrap();
    while let Some(rec) = reader.read_next().unwrap() {
        let seq = rec.seq;
        let mut prev_non_acgt_char = -1_isize;
        for (pos, (len, colex)) in index2.matching_statistics_iter(seq).enumerate() {
            if !is_dna(seq[pos]) {
                prev_non_acgt_char = pos as isize;
            }
            let dna_run_len = (pos as isize - prev_non_acgt_char) as usize;
            if dna_run_len >= k { // Valid k-mer
                if len == k { // Found in index2
                    assert!(colex.len() == 1);
                    callback((Some(colex.start), &seq[pos+1-k..pos+1]));
                } else {
                    callback((None, &seq[pos+1-k..pos+1]))
                }
            }
        }    
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SetOperation {
    Intersection,
    Union,
    Difference,
}

pub struct InputTriple<'a> {
    pub source_seqs: PathBuf,
    pub index: &'a SbwtIndexVariant,
    pub lcs: &'a LcsArray,
}

pub fn run_check_set_operation<'a>(
    op: SetOperation,
    input1: InputTriple<'a>,
    input2: InputTriple<'a>,
    result_index: &mut SbwtIndexVariant,
) {

    let (index1, seq1_source, lcs1) = (input1.index, input1.source_seqs, input1.lcs);
    let (index2, seq2_source, lcs2) = (input2.index, input2.source_seqs, input2.lcs);
    let _streaming1 = index1.get_streaming_index(lcs1);
    let streaming2 = index2.get_streaming_index(lcs2);

    let k = index1.k();
    assert_eq!(k, index2.k(), "sbwt1 and sbwt2 must have the same k");
    assert_eq!(k, result_index.k(), "the result index must have the same k as sbwt1/sbwt2");

    let mut visited_marks = result_index.compute_dummy_node_marks();

    match op {
        SetOperation::Union => {
            log::info!("Checking that every k-mer of seq1 and of seq2 is in the result (union)");
            let mut check = |(colex_result, kmer): (Option<usize>, &[u8])| {
                let colex_result = colex_result.expect(&format!("k-mer {} not found in union", String::from_utf8_lossy(kmer)));
                visited_marks.set(colex_result, true);
            };

            lookup_each_valid_kmer_nonstreaming_and_callback(&seq1_source, &result_index, &mut check);
            lookup_each_valid_kmer_nonstreaming_and_callback(&seq2_source, &result_index, &mut check);
        },
        SetOperation::Intersection => {
            log::info!("Checking that every k-mer of seq1 that's also in sbwt2 is in the result (intersection)");
            let mut check = |(colex2, kmer): (Option<usize>, &[u8])| {
                if colex2.is_some() {
                    // This k-mer is in the intersection
                    let colex_result = result_index.search(kmer).expect(&format!("k-mer {} should be in the intersection", String::from_utf8_lossy(kmer)));
                    visited_marks.set(colex_result.start, true);
                }
            };

            lookup_each_valid_kmer_and_callback(&seq1_source, &streaming2, &mut check);
        },
        SetOperation::Difference => {
            log::info!("Checking that every k-mer of seq1 that's not in sbwt2 is in the result (seq1 \\ seq2)");
            let mut check = |(colex2, kmer): (Option<usize>, &[u8])| {
                if colex2.is_none() {
                    // This k-mer is in the difference
                    let colex_result = result_index.search(kmer).expect(&format!("k-mer {} should be in the difference", String::from_utf8_lossy(kmer)));
                    visited_marks.set(colex_result.start, true);
                }
            };

            lookup_each_valid_kmer_and_callback(&seq1_source, &streaming2, &mut check);
        }
    }

    if visited_marks.count_ones() != visited_marks.len() {
        log::info!("Result index has an unexpected k-mer. Fetching it...");
        result_index.build_select();
        let unvisited = visited_marks.first_zero().unwrap();
        let kmer = result_index.access_kmer(unvisited);
        panic!(
            "result index has k-mer {}, which should not be in the {op:?} of the given sequences",
            String::from_utf8_lossy(&kmer)
        );
    }

    log::info!("OK. The result index has exactly the expected k-mer set for {op:?}.");
    println!("PASS");
}

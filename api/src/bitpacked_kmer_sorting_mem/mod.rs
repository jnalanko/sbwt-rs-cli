//! Build an [SbwtIndex] fully in memory using an algorithm based on bitpacked k-mer sorting.

mod dummies;
mod kmer_splitter;
mod cursors;

use dummies::KmersWithLengths;
use kmer_splitter::get_bitpacked_sorted_distinct_kmers;

use crate::{sbwt::{PrefixLookupTable, SbwtIndex}, streaming_index::LcsArray, subsetseq::SubsetSeq, util::DNA_ALPHABET};
/// Build SBWT and optionally the LCS array fully in memory using bitpacked k-mer sorting.
///
/// See [SbwtIndexBuilder](crate::builder::SbwtIndexBuilder) for a wrapper with a more
/// user-friendly interface. B is the number u64 words in a k-mer.
///
/// Unused arguments for signature compatibility with the disk based building algorithm.
pub fn build_with_bitpacked_kmer_sorting<const B: usize, IN: crate::SeqStream + Send, SS: SubsetSeq + Send>(
    seqs: IN,
    k: usize,
    n_threads: usize,
    approx_mem_gb: usize,
    dedup_batches: bool,
    build_lcs: bool,
    add_all_dummy_paths: bool,
) -> (SbwtIndex::<SS>, Option<LcsArray>) {

    let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();

    thread_pool.install(||{

        let sigma = DNA_ALPHABET.len();

        log::info!("Bit-packing and sorting all k-mers of all input sequences (dedup batches: {}).", dedup_batches);
        let (rev_kmers, rev_first_mers) = get_bitpacked_sorted_distinct_kmers::<B, IN>(seqs, k, n_threads, dedup_batches, add_all_dummy_paths, approx_mem_gb);

        let n_kmers = rev_kmers.len();
        log::info!("{} distinct k-mers found", n_kmers);

        log::info!("Constructing dummy k-mers");
        let dummies = dummies::get_sorted_dummies::<B>(&rev_kmers, rev_first_mers, sigma, k, n_threads);

        log::info!("{} dummy nodes constructed", dummies.len());

        log::info!("Compacting dummies and lengths");
        let dummies_and_lengths = KmersWithLengths{ // This fixes word alignment overhead
            kmers: dummies.iter().map(|pair| pair.0).collect(),
            lengths: dummies.into_iter().map(|pair| pair.1).collect()
        };

        log::info!("Constructing the sbwt subset sequence");

        let (rawrows, lcs) = cursors::build_sbwt_bit_vectors::<B>(rev_kmers, dummies_and_lengths, k, sigma, build_lcs, n_threads);

        // Create the C array
        #[allow(non_snake_case)] // C-array is an established convention in BWT indexes
        let C: Vec<usize> = crate::util::get_C_array(&rawrows);

        log::info!("Building the subset rank structure");
        let mut subsetseq = SS::new_from_bit_vectors(rawrows);
        subsetseq.build_rank();
        let n_sets = subsetseq.len();
        let (mut index, lcs) = (SbwtIndex::<SS>::from_components(
            subsetseq,
            n_kmers,
            k,
            C,
            PrefixLookupTable::new_empty(n_sets))
                                , lcs.map(LcsArray::new));

        let lut = PrefixLookupTable::new(&index, 8);
        index.set_lookup_table(lut);
        (index, lcs)

        })

}

#[cfg(test)]
mod tests {
    use std::cmp::min;

    use crate::util::is_dna;


    #[test]
    fn test_add_all_dummy_paths() {
        let k = 3;
        let seq1 = b"ACACTGANNNNNGATANNNNTGNNNAN";
        let seq2 = b"NNNNNNNNNN";
        let seq3 = b"ACCCCCCCNNNNTACNNN";
        let seqs: Vec<&[u8]> = vec![seq1, seq2, seq3];
        let (mut sbwt, lcs) = super::build_with_bitpacked_kmer_sorting::<1, _, crate::subsetseq::SubsetMatrix>(
            crate::SliceSeqStream::new(seqs.as_slice()),
            k,
            3,     // n_threads
            1,     // approx_mem_gb
            false, // dedup_batches
            true, // build_lcs
            true,  // add_all_dummy_paths
        );
        sbwt.build_select();

        for seq in [seq1.as_ref(), seq2.as_ref(), seq3.as_ref()] {
            crate::util::for_each_run_with_key(seq, |c| crate::util::is_dna(*c), |run_range| {
                if crate::util::is_dna(seq[run_range.start]) {
                    let s = &seq[run_range];
                    // All dummy prefixes up to length k-1 should be in the index
                    for l in 1..min(k, s.len()) {
                        let colex_range = sbwt.search(&s[0..l]).unwrap();
                        assert!(colex_range.len() > 0);
                        let smallest = sbwt.access_kmer(colex_range.start);
                        let n_dna = smallest.iter().filter(|c| is_dna(**c)).count();
                        eprintln!("{}", String::from_utf8(smallest).unwrap());
                        assert_eq!(n_dna, l);
                    }
                }
            });
        }

        // Collect true k-mers from the index (rows where all k characters are DNA)
        let index_kmers: std::collections::HashSet<Vec<u8>> = (0..sbwt.n_sets())
            .map(|i| sbwt.access_kmer(i))
            .filter(|kmer| kmer.iter().all(|&c| is_dna(c)))
            .collect();

        // Collect k-mers from the input sequences
        let mut input_kmers: std::collections::HashSet<Vec<u8>> = std::collections::HashSet::new();
        for seq in [seq1.as_ref(), seq2.as_ref(), seq3.as_ref()] {
            crate::util::for_each_run_with_key(seq, |c| is_dna(*c), |run_range| {
                if is_dna(seq[run_range.start]) {
                    let s = &seq[run_range];
                    for i in 0..s.len().saturating_sub(k - 1) {
                        input_kmers.insert(s[i..i + k].to_vec());
                    }
                }
            });
        }

        assert_eq!(index_kmers, input_kmers);

    }
}
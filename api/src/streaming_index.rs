//! A streaming query index that uses right extensions and left contractions to find matches.

use crate::subsetseq::*;
use crate::sbwt::*;
use simple_sds_sbwt::ops::Access;
use simple_sds_sbwt::ops::Push;
use simple_sds_sbwt::ops::Vector;
use simple_sds_sbwt::serialize::Serialize;

/// Extending a search pattern to the right.
pub trait ExtendRight {
    /// Takes the interval I of a pattern P in an SBWT and returns the interval of
    /// pattern Pc.
    fn extend_right(&self, I: std::ops::Range<usize>, c: u8) -> std::ops::Range<usize>;
}

/// Contracting a search pattern from the left.
pub trait ContractLeft {
    /// Takes the interval I of a pattern cP in an SBWT, where c is a single character, and returns the interval of
    /// pattern P. If I is the empty string, then this should return back the interval
    /// of the empty string.
    fn contract_left(&self, I: std::ops::Range<usize>, target_len: usize) -> std::ops::Range<usize>;
}

/// An index that uses right extensions and left contractions to find matches in a streaming fashion.
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct StreamingIndex<'a, E: ExtendRight, C: ContractLeft> {
    extend_right: &'a E,
    contract_left: &'a C,
    n: usize, // The empty string has colex range [0,n)
    k: usize, // Maximum length of a match
}

impl<'a, SS: SubsetSeq> StreamingIndex<'a, SbwtIndex<SS>, LcsArray> {

    /// Create a new streaming index using an SBWT index for right extensions
    /// and an LCS array for left contractions.
    pub fn new(sbwt: &'a SbwtIndex<SS>, lcs: &'a LcsArray) -> Self {
        StreamingIndex{contract_left: lcs, extend_right: sbwt, n: sbwt.n_sets(), k: sbwt.k()}
    }
}

/// An array that stores the lengths of the longest common suffixes of consecutive k-mers in 
/// colexicographic order. 
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct LcsArray {
    lcs: simple_sds_sbwt::int_vector::IntVector,
}

impl ContractLeft for LcsArray {

    /// Expands the interval by scanning left and right until the LCS value is less than `target_len`.
    /// This has worst case time complexity equal to the length of LCS array, but if used on short intervals, the
    /// performance is typically quite good. 
    #[allow(non_snake_case)]
    fn contract_left(&self, I: std::ops::Range<usize>, target_len: usize) -> std::ops::Range<usize> {
        let mut new_start = I.start;
        let mut new_end = I.end;
        while new_start > 0 && self.lcs.get(new_start) >= target_len as u64 {
            new_start -= 1;
        }
        let n = self.lcs.len();
        while new_end < n && self.lcs.get(new_end) >= target_len as u64 {
            new_end += 1;
        }
        new_start..new_end
    }
}

// Any SbwtIndex can be used for right extensions.
impl<SS: SubsetSeq> ExtendRight for SbwtIndex<SS> {

    /// Right extensions implemented time O(t), where t is the time for a rank query in the
    /// subset rank query implementation of `SS`. This is O(1) for [SubsetMatrix]. 
    #[allow(non_snake_case)]
    fn extend_right(&self, I: std::ops::Range<usize>, c: u8) -> std::ops::Range<usize> {
        let c_idx = self.char_idx(c);
        self.lf_step(I.start, c_idx)..self.lf_step(I.end, c_idx)
    }

}

/// An LCS array implementation using ceil(log(k)) bits per element. 
impl LcsArray {

    /// Construct from a pre-built LCS array stored in an [IntVector](simple_sds_sbwt::int_vector::IntVector).
    pub fn new(v: simple_sds_sbwt::int_vector::IntVector) -> Self {
        Self{lcs: v}
    }
    
    /// Access the LCS value at position i.
    pub fn access(&self, i: usize) -> usize {
        self.lcs.get(i) as usize
    }

    /// Returns the number of elements in the LCS array.
    pub fn len(&self) -> usize {
        self.lcs.len()
    }

    /// Returns true iff the length of the array is zero.
    pub fn is_empty(&self) -> bool {
        self.lcs.len() == 0
    }

    /// Returns the number of bytes used by the LCS array.
    pub fn size_in_bytes(&self) -> usize {
        self.lcs.size_in_bytes()
    }

    /// Serialize the LCS array to the given writer.
    /// The array can then later be loaded with [LcsArray::load].
    /// Returns number of bytes written.
    pub fn serialize<W: std::io::Write>(&self, out: &mut W) -> std::io::Result<usize>{
        let mut n_written = 0_usize;
        self.lcs.serialize(out)?;
        n_written += self.lcs.size_in_bytes();
        Ok(n_written)
    }

    /// Load an LCS array that was previously serialized with [LcsArray::serialize].
    pub fn load<R: std::io::Read>(input: &mut R) -> std::io::Result<Self>{
        let lcs = simple_sds_sbwt::int_vector::IntVector::load(input)?;
        Ok(Self{lcs})
    }

    /// Build an LCS array using the given SBWT index. This runs an O(nk) time algorithm, where
    /// n is the number of k-mers in the padded spectrum encoded in the SBWT, and k is the length
    /// of the k-mers. The algorithm works by inverting the SBWT column by column, which takes
    /// 2n bytes of extra working space. The algorithm is described in the paper
    /// [Longest Common Prefix Arrays for Succinct k-spectra](https://doi.org/10.48550/arXiv.2306.04850).
    pub fn from_sbwt<SS: SubsetSeq>(sbwt: &SbwtIndex<SS>) -> Self {
        // The array stores values in the range [0..k-1], so we need
        // ceil(log2(k)) = bits per element.
        let k = sbwt.k();
        let bit_width = 64 - u64::leading_zeros((k as isize -1) as u64);
        let mut lcs = simple_sds_sbwt::int_vector::IntVector::with_capacity(sbwt.n_sets(), bit_width as usize).unwrap();
        let n_nodes = sbwt.n_sets();
        for _ in 0..n_nodes {
            lcs.push(0);
        }

        // Build the last column of the SBWT matrix
        log::info!("Building column {} of the SBWT matrix", k-1);
        let mut last = sbwt.build_last_column(); 

        let mut computed_values = bitvec::bitvec![0; sbwt.n_sets()];
        lcs.set(0, 0);
        computed_values.set(0, true);

        // Iterate columns from right to left
        for round in 0..k {
            if round > 0 {
                log::info!("Building column {} of the SBWT matrix", k-1-round);
                last = sbwt.push_labels_forward(&last);
            }

            for i in 1..n_nodes {
                if !computed_values[i] && last[i] != last[i-1] {
                    lcs.set(i, round as u64);
                    computed_values.set(i, true);
                }
            }
        }

        Self{lcs}

    }

}

impl<'a, E: ExtendRight, C: ContractLeft> StreamingIndex<'a, E, C>{

    /// Returns an array A of length `query.len()`, where A\[i\] is a pair (d, I) where d is the length of the shortest suffix
    /// of `query[0..=i]` that occurs at most `freq_bound` times in the de Bruijn graph of the k-mers in the index, and
    /// I is the colexicographic interval of that suffix. Such match may not exist, in which case A\[i\] is None.
    /// The algorithm is described in the paper [Finimizers: Variable-length bounded-frequency minimizers for k-mer sets](https://doi.org/10.1101/2024.02.19.580943).
    #[allow(non_snake_case)]
    pub fn shortest_freq_bound_suffixes(&self, query: &[u8], freq_bound: usize) -> Vec<Option<(usize, std::ops::Range<usize>)>> {
        let mut SFS = vec![None; query.len()];
        let mut I = 0..self.n;
        let mut d = 0_usize; // String length of interval I
        for i in 0..query.len() {
            if crate::util::is_dna(query[i]) {
                let mut Ic = self.extend_right.extend_right(I.clone(), query[i]);
                while d > 0 && Ic.is_empty() {
                    I = self.contract_left.contract_left(I, d-1);
                    d -= 1;
                    Ic = self.extend_right.extend_right(I.clone(), query[i]);
                }

                if !Ic.is_empty() {
                    I = Ic;
                    d = std::cmp::min(self.k, d+1);
                }
            } else {
                I = 0..self.n;
                d = 0;
            };

            // Contract I to the sortest frequency-bounded suffix
            while I.len() <= freq_bound && d > 0 {
                let J = self.contract_left.contract_left(I.clone(), d-1);
                if J.len() <= freq_bound {
                    I = J;
                    d -= 1;
                } else {
                    break
                }
            }

            if I.len() <= freq_bound{
                SFS[i] = Some((d, I.clone()));
            }
        }
        SFS
    }

    /// Returns an array A of length `query.len()`, where A\[i\] is a pair (d, I) where d is the length of the longest
    /// suffix of `query[i-k+1..=i]` that is found in the index, and I is the colexicographic interval of the match.
    /// The algorithm is described in the paper [Finimizers: Variable-length bounded-frequency minimizers for k-mer sets](https://doi.org/10.1101/2024.02.19.580943).
    #[allow(non_snake_case)]
    pub fn matching_statistics(&self, query: &[u8]) -> Vec<(usize, std::ops::Range<usize>)> {
        let mut d = 0_usize;
        let mut MS = vec![(0_usize, 0..0); query.len()];
        let mut I = 0..self.n;
        for i in 0..query.len() {
            if !crate::util::is_dna(query[i]) {
                // Invalid character. Attempting to right-extend would panic, and anyway
                // we would need to contract all the way to the empty string, and that would
                // take a long time with the current LCS-scanning implementation.
                // So, we just reset the interval to the interval of the empty string.
                I = 0..self.n;
                d = 0;
            } else {
                let mut Ic = self.extend_right.extend_right(I.clone(), query[i]);
                while d > 0 && Ic.is_empty() {
                    I = self.contract_left.contract_left(I, d-1);
                    d -= 1;
                    Ic = self.extend_right.extend_right(I.clone(), query[i]);
                }
                if !Ic.is_empty() {
                    I = Ic;
                    d = std::cmp::min(self.k, d+1);
                }
            }
            MS[i] = (d, I.clone());
        }
        MS
    }

}

#[cfg(test)]
mod tests {

    use crate::builder::BitPackedKmerSorting;

    use super::*;

    #[test_log::test]
    #[allow(non_snake_case)]
    fn LCS_paper_example() {
        let seqs: Vec<&[u8]> = vec![b"AGGTAAA", b"ACAGGTAGGAAAGGAAAGT"];

        let (sbwt, lcs) = crate::builder::SbwtIndexBuilder::<BitPackedKmerSorting>::new().k(4).build_lcs(true).run_from_slices(seqs.as_slice());
        let lcs = lcs.unwrap();
        let from_sbwt = LcsArray::from_sbwt(&sbwt);

        let true_lcs = [0,0,1,3,2,2,1,1,1,0,0,2,2,1,3,3,0,2];
        for i in 0..lcs.len() {
            println!("LCS {}", lcs.access(i));
            assert_eq!(true_lcs[i], lcs.access(i));
            assert_eq!(true_lcs[i], from_sbwt.access(i));
        }

        // Test streaming support
        let SS = StreamingIndex::<SbwtIndex::<SubsetMatrix>, LcsArray>{contract_left: &lcs, extend_right: &sbwt, n: sbwt.n_sets(), k: sbwt.k()};
        let mut I = 0..sbwt.sbwt.len();
        I = SS.extend_right.extend_right(I, b'A') ;
        I = SS.extend_right.extend_right(I, b'G') ;
        I = SS.extend_right.extend_right(I, b'G') ;
        assert_eq!(I, 13..16);

        let I3 = SS.contract_left.contract_left(I.clone(), 3);
        assert_eq!(I3, 13..16);

        let I2 = SS.contract_left.contract_left(I.clone(), 2);
        assert_eq!(I2, 13..16);

        let I1 = SS.contract_left.contract_left(I.clone(), 1);
        assert_eq!(I1, 10..16);

        let I0 = SS.contract_left.contract_left(I.clone(), 0);
        assert_eq!(I0, 0..18);
        dbg!(I3, I2, I1, I0);

    }

    #[test]
    #[allow(non_snake_case)]
    fn finimizer_paper_example(){

        // This is a pretty weak test but it's something.
        // Using the example from the finimizer paper.
        let seqs: Vec<&[u8]> = vec![b"GTAAGTCT", b"AGGAAA", b"ACAGG", b"GTAGG", b"AGGTA"];
        let k = 4;

        let (sbwt, lcs) = crate::builder::SbwtIndexBuilder::<BitPackedKmerSorting>::new().k(k).build_lcs(true).run_from_slices(seqs.as_slice());
        let lcs = lcs.unwrap();
        let SS = StreamingIndex::new(&sbwt, &lcs);
        let MS = SS.matching_statistics(b"AAGTAA");
        eprintln!("{:?}", MS);
        let true_MS_lengths: [usize; 6] = [1,2,3,4,3,4];
        let true_MS_freqs: [usize; 6] = [7,3,1,1,1,1];
        let true_MS_colex_starts: [usize; 6] = [2,3,11,17,8,5];
        assert_eq!(MS.len(), true_MS_lengths.len());
        for i in 0..MS.len(){
            assert_eq!(MS[i].0, true_MS_lengths[i]);
            assert_eq!(MS[i].1.len(), true_MS_freqs[i]);
            assert_eq!(MS[i].1.start + 1, true_MS_colex_starts[i]); // Paper has 1-indexing
        }

        let mut SFS = SS.shortest_freq_bound_suffixes(b"AAGTAA", 1);

        let true_SFS = [None, None, Some((3,11..12)), Some((3, 17..18)), Some((2, 8..9)), Some((3, 5..6))];

        // Put our SFS in 1-based indexing
        for item in SFS.iter_mut(){
            if let Some(X) = item{
                *item = Some((X.0, X.1.start + 1..X.1.end + 1));
            }

        }
        eprintln!("{:?}", SFS);

        assert_eq!(SFS, true_SFS);

        // Let's try a query that has an invalid character
        let mut SFS = SS.shortest_freq_bound_suffixes(b"AANAAGTAA", 1);

        let true_SFS = [None, None, None, None, None, Some((3,11..12)), Some((3, 17..18)), Some((2, 8..9)), Some((3, 5..6))];

        // Put our SFS in 1-based indexing
        for item in SFS.iter_mut(){
            if let Some(X) = item{
                *item = Some((X.0, X.1.start + 1..X.1.end + 1));
            }

        }
        assert_eq!(SFS, true_SFS);

        eprintln!("{:?}", SFS);
    }


}
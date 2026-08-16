//! A module for representing a sequence of subsets of an alphabet, with support for rank and select queries
//! on the elements of the subsets.

use std::ops::Range;

use bitvec::order::Lsb0;
use simple_sds_sbwt::bit_vector::*;
use simple_sds_sbwt::ops::*;
use simple_sds_sbwt::raw_vector::*;
use simple_sds_sbwt::serialize::*;
use bitvec::prelude::*;

// Import simple_sds_sbwt::ops::BitVec explicitly.
// It is already imported in the glob import above, but there is a name conflict:
// bitvec::prelude::* brings the struct bitvec::vec::BitVec. It compiles for now,
// but future Rust versions will not support this.
// See https://github.com/rust-lang/rust/issues/147992
// This import fixes the ambiguity.
use simple_sds_sbwt::ops::BitVec;

use crate::util::bitvec_to_simple_sds_raw_bitvec;

/// This trait represents a sequence of subsets from alphabet {0, 1, ..., sigma-1}, where sigma is the alphabet size.
/// The trait provides access to the subsets and rank and select queries for the elements inside the subsets.
#[allow(clippy::len_without_is_empty)]
pub trait SubsetSeq {
    // Todo: make char type into a generic unsigned integer.
    // Issues with that: it can't seem to figure out how to use such
    // generic integers as array indexes. Lol.

    /// Creates a new subset sequence from a vector of subsets of {0, 1, ..., sigma-1}. Order
    /// of characters inside a subset does not matter. The rank and select supports are not automatically
    /// initialized, so if the you need those functions, you need to call [`SubsetSeq::build_rank`] and [`SubsetSeq::build_select`], respectively.
    fn new(subset_seq: Vec<Vec<u8>>, sigma: usize) -> Self;

    /// Create a new subset sequence from indicator [bit vectors](bitvec::vec::BitVec), where the i-th bit of the j-th bit vector
    /// is 1 if and only if the i-th subset contains the j-th character. Rank and select
    /// supports need to be initialized by calling [`SubsetSeq::build_rank`] and [`SubsetSeq::build_select`], respectively.
    fn new_from_bit_vectors(vecs: Vec<bitvec::vec::BitVec<u64, Lsb0>>) -> Self;

    /// Number of sets in the sequence (**not** the total length of the sets).
    fn len(&self) -> usize;

    /// Initialize rank support for the elements of the subsets.
    fn build_rank(&mut self);

    /// Initialize select support for the elements of the subsets.
    fn build_select(&mut self);

    /// Returns true if the select support has been initialized.
    fn has_select_support(&self) -> bool;

    /// Returns true if the rank support has been initialized.
    fn has_rank_support(&self) -> bool;

    /// Returns the number of sets in the range [0..i) that have character `c`.
    fn rank(&self, c: u8, i: usize) -> usize;

    /// Returns the index of the set that contains the i-th occurence of character c, if exists.
    fn select(&self, c: u8, i: usize) -> Option<usize>;

    /// Appends the elements of the i-th set to the buffer.
    fn append_set_to_buf(&self, i: usize, buf: &mut Vec<u8>);

    /// Returns the elements in the i-th set.
    fn access(&self, i: usize) -> Vec<u8> {
        let mut v = Vec::new();
        self.append_set_to_buf(i, &mut v);
        v
    }

    /// Returns the size of the i-th subset.
    fn subset_size(&self, i: usize) -> usize;

    /// Returns true if the set at index `set_idx` contains the character.
    fn set_contains(&self, set_idx: usize, character: u8) -> bool;

    /// Returns the index of the next set from set_idx (possibly set_idx itself)
    /// that contains character c, or None if does not exist.
    /// The index set_idx can be from 0 to self.len().
    fn next_set_with_char(&self, set_idx: usize, c: u8) -> Option<usize>;

    /// Writes the subset sequence to the given writer.
    /// The sequence can be loaded later back with [`SubsetSeq::load`].
    /// Returns the number of bytes written.
    fn serialize<W: std::io::Write>(&self, out: &mut W) -> std::io::Result<usize>;

    /// Loads a subset sequence that was previously written with [`SubsetSeq::serialize`].
    fn load<R: std::io::Read>(input: &mut R) -> std::io::Result<Self>
    where
        Self: Sized;

    // Calls the callback on the index of every set that contains character c in the range.
    // The callback is called in increasing order of set index.
    fn call_on_char_occurrences<F: FnMut(usize)>(&self, range: Range<usize>, c: u8, callback: F);

    /// Build the cumulative sum array C required in [`crate::sbwt::SbwtIndex`].
    fn get_C_array(&self) -> Vec<usize> {
        let sigma = 4; // TODO
        let n = self.len();

        let mut C: Vec<usize> = vec![0; sigma];
        for i in 0..n {
            for c in 0..(sigma as u8) {
                if self.set_contains(i, c) {
                    for d in (c + 1)..(sigma as u8) {
                        C[d as usize] += 1;
                    }
                }
            }
        }

        // Plus one for the ghost dollar
        #[allow(clippy::needless_range_loop)] // Is perfectly clear this way
        for c in 0..sigma {
            C[c] += 1;
        }

        C
    }

    // This a key subroutine for SbwtIndex::push_labels_forward. We want to put it here in the SubsetSeq
    // trait because then we have have a really optimized version for particular implementations
    // like SubsetMatrix. See the code for what it does.
    fn push_labels_forward(
        &self,
        labels: &[u8],
        sbwt_input_range: std::ops::Range<usize>,
        output_ranges: Vec<&mut [u8]>,
    );

    /// Creates an indicator [bit vectors](bitvec::vec::BitVec) for each symbol of the alphabet, where the i-th bit of the j-th bit vector
    /// is 1 if and only if the i-th subset contains the j-th character.
    fn into_bitvectors(self) -> Vec<bitvec::vec::BitVec<u64, Lsb0>> where Self: Sized {
        let sigma = 4; // TODO
        (0..sigma).map(|c| {
            let mut row = bitvec![u64, Lsb0; 0; self.len()];
            self.call_on_char_occurrences(0..self.len(), c as u8, |i| {
                row.set(i, true);
            });
            row
        }).collect()
    }
}

/// An implementation of [SubsetSeq] with a matrix of sigma indicator bit vectors: the i-th bit of the j-th bit vector
/// is 1 if and only if the i-th subset contains the j-th character. Rank and select queries are reduced to
/// bit vector rank and select queries on the indicator bit vectors.
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct SubsetMatrix {
    rows: Vec<simple_sds_sbwt::bit_vector::BitVector>,
}

impl SubsetSeq for SubsetMatrix {
    fn serialize<W: std::io::Write>(&self, out: &mut W) -> std::io::Result<usize> {
        let mut n_written = 0_usize;

        let n_rows = self.rows.len() as u64;
        out.write_all(&n_rows.to_le_bytes())?;
        for row in self.rows.iter() {
            row.serialize(out)?;
            n_written += row.size_in_bytes();
        }
        Ok(n_written)
    }

    fn load<R: std::io::Read>(input: &mut R) -> std::io::Result<Self> {
        let n_rows = u64::load(input)? as usize;

        let mut rows = Vec::<BitVector>::new();
        for _ in 0..n_rows {
            rows.push(BitVector::load(input)?);
        }
        Ok(Self { rows })
    }

    fn new_from_bit_vectors(rows: Vec<bitvec::vec::BitVec<u64, Lsb0>>) -> Self {
        Self {
            rows: rows
                .into_iter()
                .map(|row| {
                    simple_sds_sbwt::bit_vector::BitVector::from(bitvec_to_simple_sds_raw_bitvec(
                        row,
                    ))
                })
                .collect(),
        }
    }

    fn new(subset_seq: Vec<Vec<u8>>, sigma: usize) -> Self {
        let n = subset_seq.len();
        let mut rawrows = Vec::<RawVector>::new();
        for _ in 0..sigma {
            rawrows.push(RawVector::with_len(n, false));
        }
        for (i, set) in subset_seq.iter().enumerate() {
            for c in set.iter() {
                rawrows[*c as usize].set_bit(i, true)
            }
        }

        let rows: Vec<BitVector> = rawrows.into_iter().map(BitVector::from).collect();
        Self { rows }
    }

    fn build_rank(&mut self) {
        for row in self.rows.iter_mut() {
            row.enable_rank();
        }
    }

    fn build_select(&mut self) {
        for row in self.rows.iter_mut() {
            row.enable_select();
        }
    }

    fn has_select_support(&self) -> bool {
        self.rows[0].supports_select()
    }

    fn has_rank_support(&self) -> bool {
        self.rows[0].supports_rank()
    }

    fn rank(&self, c: u8, i: usize) -> usize {
        self.rows[c as usize].rank(i)
    }

    fn select(&self, c: u8, i: usize) -> Option<usize> {
        assert!(self.rows[c as usize].supports_select());
        self.rows[c as usize].select(i)
    }

    fn len(&self) -> usize {
        if self.rows.is_empty() {
            0
        } else {
            self.rows[0].len()
        }
    }

    fn subset_size(&self, i: usize) -> usize {
        self.rows.iter().filter(|&row| row.get(i)).count()
    }

    fn set_contains(&self, set_idx: usize, character: u8) -> bool {
        self.rows[character as usize].get(set_idx)
    }

    fn append_set_to_buf(&self, i: usize, buf: &mut Vec<u8>) {
        for c in 0..self.rows.len() {
            if self.rows[c].get(i) {
                buf.push(c as u8);
            }
        }
    }

    fn next_set_with_char(&self, set_idx: usize, c: u8) -> Option<usize> {
        // Use the bitvec crate on the raw data of the row
        let bv = bitvec::slice::BitSlice::<u64, Lsb0>::from_slice(
            self.rows[c as usize].as_ref().as_ref(),
        );
        let off = bv[set_idx..].first_one()?;
        Some(set_idx + off)
    }

    fn call_on_char_occurrences<F: FnMut(usize)>(
        &self,
        range: Range<usize>,
        c: u8,
        mut callback: F,
    ) {
        let bv = bitvec::slice::BitSlice::<u64, Lsb0>::from_slice(
            self.rows[c as usize].as_ref().as_ref(),
        );
        let slice = &bv[range.clone()];
        for i in slice.iter_ones() {
            callback(range.start + i);
        }
    }

    fn push_labels_forward(
        &self,
        labels: &[u8],
        sbwt_input_range: std::ops::Range<usize>,
        output_ranges: Vec<&mut [u8]>,
    ) {
        assert_eq!(labels.len(), sbwt_input_range.len());
        assert_eq!(output_ranges.len(), self.rows.len());
        if sbwt_input_range.is_empty() {
            return;
        }

        #[allow(clippy::needless_range_loop)]
        for (char_idx, output_slice) in output_ranges.into_iter().enumerate() {
            let mut output_offset = 0_usize;

            let words = self.rows[char_idx].as_ref().get_words();

            // We roll our own one-bit iterator because the one in the bitvec crate
            // has a lot of overhead. This one is something like 6x faster!
            crate::util::for_each_one_bit(words, sbwt_input_range.clone(), |i| {
                output_slice[output_offset] = labels[i - sbwt_input_range.start];
                output_offset += 1;
            });
        }
    }

    fn into_bitvectors(self) -> Vec<bitvec::vec::BitVec<u64, Lsb0>> where Self: Sized {
        // This is zero-copy throughout
        self.rows.into_iter()
          .map(simple_sds_sbwt::raw_vector::RawVector::from)
          .map(crate::util::simple_sds_raw_bitvec_to_bitvec)
          .collect()
    }
}

/// Formats the subset matrix as an ASCII bit matrix of 0s and 1s, where row i
/// is the indicator bit vector for the i-th character.
impl std::fmt::Display for SubsetMatrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for row in self.rows.iter() {
            for i in 0..row.len() {
                if row.get(i) {
                    write!(f, "1")?;
                } else {
                    write!(f, "0")?;
                }
            }
            writeln!(f)?;
        }
        Ok(())
    }
}
#[cfg(test)]
mod tests {
    //use crate::SubsetCorrectionSets;

    use super::*;

    #[test]
    fn serialize_and_load() {
        let sets: Vec<Vec<u8>> = vec![
            vec![1, 2, 3],
            vec![0, 2],
            vec![0, 1, 3, 4],
            vec![],
            vec![0, 1, 2],
        ];
        let mut sm = SubsetMatrix::new(sets, 5);
        sm.build_rank();
        let mut buf = Vec::<u8>::new();
        sm.serialize(&mut buf).unwrap();
        let sm2 = SubsetMatrix::load(&mut buf.as_slice()).unwrap();
        assert_eq!(sm, sm2);
    }

    #[test]
    fn next_set_with_char() {
        let sets: Vec<Vec<u8>> = vec![
            vec![1, 2, 3],
            vec![0, 2],
            vec![0, 1, 3, 4],
            vec![],
            vec![0, 1, 2],
        ];
        let sm = SubsetMatrix::new(sets, 5);
        assert_eq!(sm.next_set_with_char(0, 0), Some(1));
        assert_eq!(sm.next_set_with_char(0, 1), Some(0));
        assert_eq!(sm.next_set_with_char(0, 2), Some(0));
        assert_eq!(sm.next_set_with_char(0, 3), Some(0));
        assert_eq!(sm.next_set_with_char(0, 4), Some(2));

        assert_eq!(sm.next_set_with_char(3, 0), Some(4));
        assert_eq!(sm.next_set_with_char(3, 1), Some(4));
        assert_eq!(sm.next_set_with_char(3, 2), Some(4));
        assert_eq!(sm.next_set_with_char(3, 3), None);
        assert_eq!(sm.next_set_with_char(3, 4), None);

        assert_eq!(sm.next_set_with_char(5, 0), None); // One past the end
    }

    /* fn test_access_generic<SS: SubsetSeq>(ss: SS) {
        assert_eq!(ss.access(0), vec![1, 2, 3]);
    }

    fn get_example<SS: SubsetSeq>() -> SS {
        let sets: Vec<Vec<u8>> = vec![
            vec![1, 2, 3],
            vec![0, 2],
            vec![0, 1, 2, 3],
            vec![],
            vec![0, 1, 2],
        ];
        let ss = SS::new(sets, 4);
        ss
    } */

    // unit tests for each function
    /* #[test]
       fn test_access() {
           test_access_generic(get_example::<SubsetMatrix>());
           test_access_generic(get_example::<SubsetCorrectionSets>());
       }
    */
    fn get_example<SS: SubsetSeq>() -> SS {
        let sets = vec![
            vec![1, 2, 3],
            vec![0, 2],
            vec![0, 1, 2, 3],
            vec![],
            vec![0, 1, 2],
        ];

        SS::new(sets, 4)
    }

    fn test_access_generic<SS: SubsetSeq>() {
        let expected = vec![
            vec![1, 2, 3],
            vec![0, 2],
            vec![0, 1, 2, 3],
            vec![],
            vec![0, 1, 2],
        ];

        let ss = get_example::<SS>();

        for (i, set) in expected.iter().enumerate() {
            assert_eq!(ss.access(i), *set);
        }
    }

    fn test_subset_size_generic<SS: SubsetSeq>() {
        let ss = get_example::<SS>();

        let expected = [3, 2, 4, 0, 3];

        for (i, &sz) in expected.iter().enumerate() {
            assert_eq!(ss.subset_size(i), sz);
        }
    }

    fn test_set_contains_generic<SS: SubsetSeq>() {
        let ss = get_example::<SS>();

        for i in 0..ss.len() {
            let set = ss.access(i);

            for c in 0..4 {
                assert_eq!(ss.set_contains(i, c), set.contains(&(c as u8)));
            }
        }
    }

    fn test_rank_generic<SS: SubsetSeq>() {
        let mut ss = get_example::<SS>();
        ss.build_rank();

        for c in 0..4 {
            for i in 0..=ss.len() {
                let expected = (0..i)
                    .filter(|&j| ss.access(j).contains(&(c as u8)))
                    .count();

                assert_eq!(ss.rank(c as u8, i), expected, "char={}, i={}", c, i);
            }
        }
    }

    fn test_select_generic<SS: SubsetSeq>() {
        let mut ss = get_example::<SS>();
        ss.build_select();

        for c in 0..4 {
            let occs: Vec<_> = (0..ss.len()).filter(|&i| ss.set_contains(i, c)).collect();

            for (k, &expected) in occs.iter().enumerate() {
                assert_eq!(ss.select(c, k), Some(expected));
            }

            assert_eq!(ss.select(c, occs.len()), None);
        }
    }

    fn test_next_set_with_char_generic<SS: SubsetSeq>() {
        let ss = get_example::<SS>();

        for c in 0..4 {
            for start in 0..=ss.len() {
                let expected = (start..ss.len()).find(|&i| ss.set_contains(i, c));

                assert_eq!(ss.next_set_with_char(start, c), expected);
            }
        }
    }

    fn test_call_on_occurrences_generic<SS: SubsetSeq>() {
        let ss = get_example::<SS>();

        for c in 0..4 {
            let mut got = Vec::new();

            ss.call_on_char_occurrences(1..4, c, |i| got.push(i));

            let expected: Vec<_> = (1..4).filter(|&i| ss.set_contains(i, c)).collect();

            assert_eq!(got, expected);
        }
    }

    fn test_serialize_generic<SS>()
    where
        SS: SubsetSeq + Eq + std::fmt::Debug,
    {
        let mut ss = get_example::<SS>();
        ss.build_rank();
        ss.build_select();

        let mut bytes = Vec::new();

        ss.serialize(&mut bytes).unwrap();

        let loaded = SS::load(&mut bytes.as_slice()).unwrap();

        assert_eq!(ss, loaded);
    }

    fn test_rank_select_consistency<SS: SubsetSeq>() {
        let mut ss = get_example::<SS>();
        ss.build_rank();
        ss.build_select();

        for c in 0..4 {
            let occs: Vec<_> = (0..ss.len()).filter(|&i| ss.set_contains(i, c)).collect();

            for (k, &pos) in occs.iter().enumerate() {
                assert_eq!(ss.select(c, k), Some(pos));
                assert_eq!(ss.rank(c, pos), k);
                assert_eq!(ss.rank(c, pos + 1), k + 1);
                assert!(ss.set_contains(pos, c));
            }
        }
    }

    fn test_access_consistency<SS: SubsetSeq>() {
        let ss = get_example::<SS>();

        for i in 0..ss.len() {
            let access = ss.access(i);

            assert_eq!(access.len(), ss.subset_size(i));

            for c in 0..4 {
                assert_eq!(access.contains(&(c as u8)), ss.set_contains(i, c));
            }
        }
    }

    fn run_subset_seq_tests<SS>()
    where
        SS: SubsetSeq + Eq + std::fmt::Debug,
    {
        test_access_generic::<SS>();
        test_subset_size_generic::<SS>();
        test_set_contains_generic::<SS>();
        test_rank_generic::<SS>();
        test_select_generic::<SS>();
        test_next_set_with_char_generic::<SS>();
        test_call_on_occurrences_generic::<SS>();
        test_serialize_generic::<SS>();
        test_rank_select_consistency::<SS>();
        test_access_consistency::<SS>();
    }

    #[test]
    fn subset_matrix_generic_tests() {
        run_subset_seq_tests::<SubsetMatrix>();
    }

    /*
    #[test]
    fn subset_correction_sets_generic_tests() {
        run_subset_seq_tests::<SubsetCorrectionSets>();
    }

    use rand::rngs::StdRng;
    use rand::{Rng, SeedableRng};

    fn build_or_dump(sets: Vec<Vec<u8>>) -> SubsetCorrectionSets {
        match std::panic::catch_unwind(|| SubsetCorrectionSets::new(sets.clone(), 4)) {
            Ok(ss) => ss,
            Err(e) => {
                eprintln!("Failing input:");
                eprintln!("{:#?}", sets.len());
                std::panic::resume_unwind(e);
            }
        }
    }

    #[test]
    fn differential_random() {
        let mut rng = StdRng::seed_from_u64(12345);

        for _case in 0..1000 {
            let n = rng.gen_range(1, 1000);

            let mut sets = Vec::with_capacity(n);

            for _ in 0..n {
                let mut set = Vec::new();
                for c in 0..4 {
                    if rng.gen_bool(0.5) {
                        set.push(c);
                    }
                }
                sets.push(set);
            }

            let mut matrix = SubsetMatrix::new(sets.clone(), 4);
            matrix.build_rank();
            matrix.build_select();

            let corr = build_or_dump(sets.clone());

            assert_eq!(matrix.len(), corr.len());
            // access / subset_size
            // dbg!("Access/subset size");
            for i in 0..matrix.len() {
                assert_eq!(matrix.access(i), corr.access(i));
                assert_eq!(matrix.subset_size(i), corr.subset_size(i));

                for c in 0..4 {
                    assert_eq!(
                        matrix.set_contains(i, c),
                        corr.set_contains(i, c),
                        "set_contains({}, {})",
                        i,
                        c
                    );
                }
            }
            // rank
            //dbg!("Rank");
            for c in 0..4 {
                for i in 0..=matrix.len() {
                    assert_eq!(matrix.rank(c, i), corr.rank(c, i), "rank({}, {})", c, i);
                }
            }
            // select
            //dbg!("Select");
            for c in 0..4 {
                let occs = matrix.rank(c, matrix.len());
                for k in 0..=occs {
                    assert_eq!(
                        matrix.select(c, k),
                        corr.select(c, k),
                        "select({}, {})",
                        c,
                        k
                    );
                }
            }

            // next_set_with_char
            //dbg!("next_set_with_char");
            for c in 0..4 {
                for start in 0..=matrix.len() {
                    assert_eq!(
                        matrix.next_set_with_char(start, c),
                        corr.next_set_with_char(start, c),
                        "next({}, {})",
                        start,
                        c
                    );
                }
            }
        }
    }

    #[test]
    fn differential_call_on_occurrences() {
        let mut rng = StdRng::seed_from_u64(999);

        for _ in 0..500 {
            let n = rng.gen_range(1, 300);

            let mut sets = Vec::new();

            for _ in 0..n {
                let mut set = Vec::new();
                for c in 0..4 {
                    if rng.gen_bool(0.5) {
                        set.push(c);
                    }
                }
                sets.push(set);
            }

            let matrix = SubsetMatrix::new(sets.clone(), 4);
            let corr = build_or_dump(sets);

            for _ in 0..50 {
                let a = rng.gen_range(0, n);
                let b = rng.gen_range(a, n);

                for c in 0..4 {
                    let mut v1 = Vec::new();
                    matrix.call_on_char_occurrences(a..b, c, |i| v1.push(i));

                    let mut v2 = Vec::new();
                    corr.call_on_char_occurrences(a..b, c, |i| v2.push(i));

                    assert_eq!(v1, v2);
                }
            }
        }
    }

    #[test]
    fn differential_push_labels_forward() {
        let mut rng = StdRng::seed_from_u64(7);

        for _ in 0..500 {
            let n = rng.gen_range(1, 300);

            let mut sets = Vec::new();

            for _ in 0..n {
                let mut set = Vec::new();
                for c in 0..4 {
                    if rng.gen_bool(0.5) {
                        set.push(c);
                    }
                }
                sets.push(set);
            }

            let matrix = SubsetMatrix::new(sets.clone(), 4);
            let corr = build_or_dump(sets.clone());
            let labels: Vec<u8> = (0..n).map(|_| rng.gen()).collect();

            let a = rng.gen_range(0, n);
            let b = rng.gen_range(a, n);

            let range = a..b;

            let mut out1 = vec![Vec::<u8>::new(); 4];
            let mut out2 = vec![Vec::<u8>::new(); 4];

            for c in 0..4 {
                let count = (range.clone())
                    .filter(|&i| matrix.set_contains(i, c as u8))
                    .count();

                out1[c].resize(count, 0);
                out2[c].resize(count, 0);
            }

            matrix.push_labels_forward(
                &labels[range.clone()],
                range.clone(),
                out1.iter_mut().map(Vec::as_mut_slice).collect(),
            );

            corr.push_labels_forward(
                &labels[range.clone()],
                range.clone(),
                out2.iter_mut().map(Vec::as_mut_slice).collect(),
            );

            assert_eq!(out1, out2);
        }
    }
    */
}

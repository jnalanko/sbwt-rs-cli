//! A module for representing a sequence of subsets of an alphabet, with support for rank and select queries
//! on the elements of the subsets.

use simple_sds_sbwt::bit_vector::*;
use simple_sds_sbwt::raw_vector::*;
use simple_sds_sbwt::ops::*;
use simple_sds_sbwt::serialize::*;

/// This trait represents a sequence of subsets from alphabet {0, 1, ..., sigma-1}, where sigma is the alphabet size.
/// The trait provides access to the subsets and rank and select queries for the elements inside the subsets.
#[allow(clippy::len_without_is_empty)]
pub trait SubsetSeq{

    // Todo: make char type into a generic unsigned integer.
    // Issues with that: it can't seem to figure out how to use such
    // generic integers as array indexes. Lol.

    /// Creates a new subset sequence from a vector of subsets of {0, 1, ..., sigma-1}. Order
    /// of characters inside a subset does not matter. The rank and select supports are not automatically
    /// initialized, so if the you need those functions, you need to call [`SubsetSeq::build_rank`] and [`SubsetSeq::build_select`], respectively.
    fn new(subset_seq: Vec<Vec<u8>>, sigma: usize) -> Self;

    /// Create a new subset sequence from indicator [bit vectors](simple_sds_sbwt::bit_vector::BitVector), where the i-th bit of the j-th bit vector
    /// is 1 if and only if the i-th subset contains the j-th character. The resulting subset sequence has
    /// rank and select support if the provided bit vectors have rank and select support enabled. Otherwise, those
    /// supports need to be initialized by calling [`SubsetSeq::build_rank`] and [`SubsetSeq::build_select`], respectively.
    fn new_from_bit_vectors(vecs: Vec<simple_sds_sbwt::bit_vector::BitVector>) -> Self;

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
    fn access(&self, i: usize) -> Vec<u8>{
        let mut v = Vec::new();
        self.append_set_to_buf(i, &mut v);
        v
    }

    /// Returns the size of the i-th subset.
    fn subset_size(&self, i: usize) -> usize;

    /// Returns true if the set at index `set_idx` contains the character.
    fn set_contains(&self, set_idx: usize, character: u8) -> bool;

    /// Writes the subset sequence to the given writer.
    /// The sequence can be loaded later back with [`SubsetSeq::load`].
    /// Returns the number of bytes written.
    fn serialize<W: std::io::Write>(&self, out: &mut W) -> std::io::Result<usize>;

    /// Loads a subset sequence that was previously written with [`SubsetSeq::serialize`].
    fn load<R: std::io::Read>(input: &mut R) -> std::io::Result<Self> where Self: Sized;
}

/// An implementation of [SubsetSeq] with a matrix of sigma indicator bit vectors: the i-th bit of the j-th bit vector
/// is 1 if and only if the i-th subset contains the j-th character. Rank and select queries are reduced to
/// bit vector rank and select queries on the indicator bit vectors.
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct SubsetMatrix{
    rows: Vec<BitVector>
}

impl SubsetSeq for SubsetMatrix{
    
    fn serialize<W: std::io::Write>(&self, out: &mut W) -> std::io::Result<usize>{
        let mut n_written = 0_usize;

        let n_rows = self.rows.len() as u64;
        out.write_all(&n_rows.to_le_bytes())?;
        for row in self.rows.iter(){
            row.serialize(out)?;
            n_written += row.size_in_bytes();
        }
        Ok(n_written)
    }

    fn load<R: std::io::Read>(input: &mut R) -> std::io::Result<Self>{

        let n_rows = u64::load(input)? as usize;

        let mut rows = Vec::<BitVector>::new();
        for _ in 0..n_rows{
            rows.push(BitVector::load(input)?);
        }
        Ok(Self{rows})
    }

    fn new_from_bit_vectors(rows: Vec<simple_sds_sbwt::bit_vector::BitVector>) -> Self{
        Self{rows}
    }

    fn new(subset_seq: Vec<Vec<u8>>, sigma: usize) -> Self{
        let n = subset_seq.len();
        let mut rawrows = Vec::<RawVector>::new();
        for _ in 0..sigma{
            rawrows.push(RawVector::with_len(n, false));
        }
        for (i, set) in subset_seq.iter().enumerate(){
            for c in set.iter(){
                rawrows[*c as usize].set_bit(i, true)
            }
        }

        let rows: Vec<BitVector> = rawrows.into_iter().map(BitVector::from).collect();
        Self::new_from_bit_vectors(rows)
    }

    fn build_rank(&mut self) {
        for row in self.rows.iter_mut(){
            row.enable_rank();
        }
    }

    fn build_select(&mut self) {
        for row in self.rows.iter_mut(){
            row.enable_select();
        }
    }

    fn has_select_support(&self) -> bool {
        self.rows[0].supports_select()
    }

    fn has_rank_support(&self) -> bool {
        self.rows[0].supports_rank()
    }

    fn rank(&self, c: u8, i: usize) -> usize{
        self.rows[c as usize].rank(i)
    }

    fn select(&self, c: u8, i: usize) -> Option<usize>{
        assert!(self.rows[c as usize].supports_select());
        self.rows[c as usize].select(i)
    }

    fn len(&self) -> usize{
        if self.rows.is_empty() { 0 }
        else { self.rows[0].len() }
    }

    fn subset_size(&self, i: usize) -> usize {
        self.rows.iter().filter(|&row| row.get(i)).count()
    }

    fn set_contains(&self, set_idx: usize, character: u8) -> bool {
        self.rows[character as usize].get(set_idx)
    }

    fn append_set_to_buf(&self, i: usize, buf: &mut Vec<u8>){
        for c in 0..self.rows.len(){
            if self.rows[c].get(i){
                buf.push(c as u8);
            }
        }
    }
}

/// Formats the subset matrix as an ASCII bit matrix of 0s and 1s, where row i
/// is the indicator bit vector for the i-th character.
impl std::fmt::Display for SubsetMatrix{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for row in self.rows.iter(){
            for i in 0..row.len(){
                if row.get(i){
                    write!(f, "1")?;
                } else{
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
    use super::*;

    #[test]
    fn serialize_and_load(){
        let sets: Vec<Vec<u8>> = vec![vec![1,2,3], vec![0,2], vec![0,1,3,4], vec![], vec![0,1,2]];
        let mut sm = SubsetMatrix::new(sets, 5);
        sm.build_rank();
        let mut buf = Vec::<u8>::new();
        sm.serialize(&mut buf).unwrap();
        let sm2 = SubsetMatrix::load(&mut buf.as_slice()).unwrap();
        assert_eq!(sm, sm2);
    }
}
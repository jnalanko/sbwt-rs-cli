use std::cmp::min;
use std::ops::Range;

use crate::kmer::Kmer;
use crate::kmer::LongKmer;
use crate::util::binary_search_leftmost_that_fulfills_pred;

use simple_sds_sbwt::ops::Select;
use simple_sds_sbwt::ops::SelectZero;
use simple_sds_sbwt::bit_vector::BitVector;
use simple_sds_sbwt::raw_vector::*;
use rayon::prelude::*;

pub struct KmersWithLengths<const B: usize> {
    pub kmers : Vec<LongKmer<B>>,
    pub lengths: Vec<u8>, // B can be at most 8, which means max length 256
}

impl<const B: usize> KmersWithLengths<B> {
    pub fn get(&self, i: usize) -> (LongKmer<B>, u8) {
        (self.kmers[i], self.lengths[i])
    }

    pub fn len(&self) -> usize {
        self.kmers.len()
    }

    // Returns the index of the first k-mer that starts with c,
    // self.len() if not found.
    pub fn first_that_starts_with(&self, c: u8) -> usize {
        binary_search_leftmost_that_fulfills_pred(
            |i| self.get(i), 
            |(kmer,len)| len > 0 && kmer.get_from_left(0) >= c, 
            self.len()
        )
    }
}

pub struct KmerDummyMergeSlice<'a, const B: usize> {
    all_dummies: &'a KmersWithLengths<B>,
    all_kmers: &'a [LongKmer<B>],
    dummy_idx: usize,
    kmer_idx: usize,
    dummy_end: usize,
    kmer_end: usize,
    k: usize,
}

impl<'a, const B:usize> KmerDummyMergeSlice<'a, B> {

    fn advance_dummy(&mut self) -> (LongKmer<B>, u8) {
        let x = self.all_dummies.get(self.dummy_idx);
        self.dummy_idx += 1;
        x
    }

    fn advance_kmers(&mut self) -> (LongKmer<B>, u8) {
        let x = self.all_kmers[self.dummy_idx];
        self.kmer_idx += 1;
        (x, self.k as u8)
    }

    pub fn next(&mut self) -> Option<(LongKmer<B>, u8)> {
        if self.dummy_idx == self.dummy_end && self.kmer_idx == self.kmer_end {
            None
        } else if self.dummy_idx == self.dummy_end {
            Some(self.advance_kmers())
        } else if self.kmer_idx == self.kmer_end {
            Some(self.advance_dummy())
        } else if (self.all_kmers[self.kmer_idx], self.k as u8) < self.all_dummies.get(self.dummy_idx) {
            Some(self.advance_kmers())
        } else {
            Some(self.advance_dummy())
        }
    }

    pub fn new(all_dummies: &'a KmersWithLengths<B>, all_kmers: &'a [LongKmer<B>], merged_range: Range<usize>, k: usize) -> Self {
        assert!(merged_range.end <= all_dummies.len() + all_kmers.len());

        todo!();
        // Binary search the starts and ends in kmers and dummies.



        /*Self {
            all_dummies, all_kmers, 
            dummy_idx: dummy_range.start, dummy_end: dummy_range.end,
            kmer_idx: kmer_range.start, kmer_end: kmer_range.end,
            k
        }*/
    }
}

// Assumes the two sorted lists have no duplicates in them or between them.
// Returns (pos_a, pos_b, take_from_a)
fn binary_search_position_in_merged_list<T: PartialOrd + Eq, Access1: Fn(usize) -> T, Access2: Fn(usize) -> T>(access_a: Access1, access_b: Access2, target_pos: usize, len_a: usize, len_b: usize) -> (usize, usize, bool) {

    assert!(target_pos < len_a + len_b); // Must be a valid index in the merged list
    assert!(len_a > 0); // TODO: don't assume this
    assert!(len_b > 0); // TODO: don't assume this

    // Guess an index i in a and find how many elements of b are smaller than that (index j in b)
    // The element a[i] is at index i+j in the merged list. Binary search with this info. If the sought-after
    // merged index was found in a, we're good. Otherwise, the merged index is in b. Then we m be the number
    // of elements we're short of the target. The answer is then j+m.
    let is_ge_target = |a_idx| {
        let b_idx = binary_search_leftmost_that_fulfills_pred(|j| j, |b_idx| access_b(b_idx) > access_a(min(a_idx, len_a-1)), len_b);
        dbg!(a_idx, b_idx);
        a_idx + b_idx >= target_pos
    };

    let a_idx = binary_search_leftmost_that_fulfills_pred(|i| i, is_ge_target, len_a);
    dbg!("Landed on", a_idx);

    // Find the corresponding b position (could remember this from the first search but whatever, it's just one more search).
    let b_idx = binary_search_leftmost_that_fulfills_pred(|j| j, |j| access_b(j) > access_a(min(a_idx, len_a-1)), len_b);
    dbg!("Found", target_pos, a_idx, b_idx);

    if a_idx + b_idx == target_pos {
        (a_idx, b_idx, true)
    } else { 
        // The target position is in b. a[a_idx] is larger than the answer in b (or a_idx == a_len).
        // That means we take a_idx from a. and the rest from b. 
        (a_idx, target_pos - a_idx, false)
    }
}

fn find_first_starting_with<const B: usize>(
    v: &[LongKmer<B>],
    c: u8,
) -> usize {
    binary_search_leftmost_that_fulfills_pred(
        |pos| { v[pos] },
        |kmer: LongKmer::<B>| { kmer.get_from_left(0) >= c },
        v.len()
    )
}

pub fn get_has_predecessor_marks<const B: usize>(
    kmers: &[LongKmer::<B>],
    char_slice: (&[LongKmer<B>], usize),
    k: usize,
    c: u8,
) -> simple_sds_sbwt::raw_vector::RawVector {
    let mut bits = simple_sds_sbwt::raw_vector::RawVector::new();
    bits.resize(kmers.len(), false); // Todo: only need the colex slice for c, not the whole range
    let mut pointed_idx = 0;
    kmers.iter().for_each(|x| {
        let xc = x.copy_set_from_left(k-1, 0).right_shifted(1).copy_set_from_left(0, c);

        while pointed_idx < char_slice.0.len() {
            match char_slice.0[pointed_idx].cmp(&xc) {
                std::cmp::Ordering::Greater => {
                    break
                },
                std::cmp::Ordering::Equal => {
                    bits.set_bit(pointed_idx + char_slice.1, true);
                    pointed_idx += 1; // Advance
                    break
                },
                std::cmp::Ordering::Less => {
                    pointed_idx += 1; // Advance
                    // no break
                }
            }
        }
    });

    bits
}

// We take in a path and not a file object because we need multiple readers to the same file
pub fn get_sorted_dummies<const B: usize>(
    sorted_kmers: &[LongKmer::<B>],
    sigma: usize, k: usize,
) -> Vec<(LongKmer<B>, u8)> {
    // Number of k-mers in file
    let n = sorted_kmers.len();

    let mut char_cursors: Vec<(&[LongKmer::<B>], usize)> = (0..sigma).map(|c|{
        let start = find_first_starting_with(sorted_kmers, c as u8);
        let end = find_first_starting_with(sorted_kmers, c as u8 + 1);
        (&sorted_kmers[start..end], start)
    }).collect();

    log::info!("Identifying k-mers without predecessors");
    let has_predecessor = char_cursors.par_iter_mut().enumerate().map(|(c, cursor)| {
        get_has_predecessor_marks(sorted_kmers, *cursor, k, c as u8)
    }).reduce(|| {
        let mut res = simple_sds_sbwt::raw_vector::RawVector::new();
        res.resize(n, false);
        res
    }, |mut a, b| {
        let bv = BitVector::from(b); // Todo: unnecessary copy
        bv.one_iter().for_each(|idx| a.set_bit(idx.1, true)); // Todo: 64 bits at a time
        // Both of the todos above become irrelevant if we just build
        // the colex slices for each char and concatenate the results.
        // But the best way to do this could be to split the colex space into
        // n_threads chunks.
        a
    });


    log::info!("Constructing dummy k-mers");
    let iterable = BitVector::from(has_predecessor);
    let mut required_dummies: Vec::<(LongKmer::<B>, u8)> = iterable.zero_iter().par_bridge().map(|x| {
        let mut prefix = sorted_kmers[x.1];
        (0..k).collect::<Vec<usize>>().iter().map(|i| {
            let len = k - i - 1;
            prefix = prefix.left_shifted(1);
            (prefix, len as u8)
        }).collect::<Vec<(LongKmer::<B>, u8)>>()
    }).flatten().collect();

    // We always assume that the empty k-mer exists. This assumption is reflected in the C-array
    // later, which adds one "ghost dollar" count to all counts.
    required_dummies.push((LongKmer::<B>::from_ascii(b"").unwrap(), 0));

    log::info!("Sorting dummy k-mers");
    required_dummies.par_sort_unstable();
    required_dummies.dedup();
    required_dummies.shrink_to_fit();

    required_dummies
}

#[cfg(test)]
mod tests {
    use super::binary_search_position_in_merged_list;

    #[test]
    fn test_binary_search_merged_list() {
        let v1: Vec<usize> = vec![0,1,2,5,7,8,10];
        let v2: Vec<usize> = vec![3,4,6,9,11];

        let mut merged: Vec<(usize, usize, bool)> = vec![]; // (i_v1, i_v2, b). b tells which vector it's from
        merged.extend(v1.iter().enumerate().map(|x| (*x.1, x.0, true)));
        merged.extend(v2.iter().enumerate().map(|x| (*x.1, x.0, false)));
        merged.sort();
        eprintln!("{:?}", &merged);

        for query in 0..12 {
            let mut true_i = 0;
            let mut true_j = 0;
            for (_,_,from_a) in &merged[0..query] {
                true_i += *from_a as usize;
                true_j += !(*from_a) as usize;
            }
            let true_from_a = merged[query].2;
            let (i,j,from_a) = binary_search_position_in_merged_list(|i| v1[i], |j| v2[j], query, v1.len(), v2.len());
            dbg!((i,j,from_a), (true_i, true_j, true_from_a));
            assert_eq!((i,j,from_a), (true_i, true_j, true_from_a));

        }

    }
//fn binary_search_position_in_merged_list<T: PartialOrd + Eq, Access: Fn(usize) -> T>(access_a: Access, access_b: Access, target_pos: usize, len_a: usize, len_b: usize) -> (usize, usize, bool) {
}
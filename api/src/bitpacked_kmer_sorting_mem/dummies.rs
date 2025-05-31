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
    dummy_range: Range<usize>,
    kmer_range: Range<usize>,
    k: usize,
}

impl<'a, const B:usize> KmerDummyMergeSlice<'a, B> {

    pub fn len(&self) -> usize {
        self.dummy_range.len() + self.kmer_range.len()
    }

    pub fn cur_merged_index(&self) -> usize {
        self.dummy_idx + self.kmer_idx
    }

    fn advance_dummy(&mut self) -> (LongKmer<B>, u8) {
        let x = self.all_dummies.get(self.dummy_idx);
        self.dummy_idx += 1;
        x
    }

    fn advance_kmers(&mut self) -> (LongKmer<B>, u8) {
        let x = self.all_kmers[self.kmer_idx];
        self.kmer_idx += 1;
        (x, self.k as u8)
    }

    pub fn peek(&mut self) -> Option<(LongKmer<B>, u8)> {
        // Todo: refactor this logic together with next()
        if self.dummy_idx == self.dummy_range.end && self.kmer_idx == self.kmer_range.end {
            None
        } else if self.dummy_idx == self.dummy_range.end {
            Some((self.all_kmers[self.kmer_idx], self.k as u8))
        } else if self.kmer_idx == self.kmer_range.end {
            Some(self.all_dummies.get(self.dummy_idx))
        } else if (self.all_kmers[self.kmer_idx], self.k as u8) < self.all_dummies.get(self.dummy_idx) {
            Some((self.all_kmers[self.kmer_idx], self.k as u8))
        } else {
            Some(self.all_dummies.get(self.dummy_idx))
        }
    }

    pub fn next(&mut self) -> Option<(LongKmer<B>, u8)> {
        if self.dummy_idx == self.dummy_range.end && self.kmer_idx == self.kmer_range.end {
            None
        } else if self.dummy_idx == self.dummy_range.end {
            Some(self.advance_kmers())
        } else if self.kmer_idx == self.kmer_range.end {
            Some(self.advance_dummy())
        } else if (self.all_kmers[self.kmer_idx], self.k as u8) < self.all_dummies.get(self.dummy_idx) {
            Some(self.advance_kmers())
        } else {
            Some(self.advance_dummy())
        }
    }

    pub fn new(all_dummies: &'a KmersWithLengths<B>, all_kmers: &'a [LongKmer<B>], merged_range: Range<usize>, k: usize) -> Self {
        assert!(merged_range.end <= all_dummies.len() + all_kmers.len());

        // Binary search the starts and ends in kmers and dummies.

        let (kmer_start, dummy_start, _) = binary_search_position_in_merged_list(
            |i| (all_kmers[i], k as u8), 
            |j| all_dummies.get(j), 
            merged_range.start, 
            all_kmers.len(), 
            all_dummies.len()
        );

        let (kmer_end, dummy_end, _) = binary_search_position_in_merged_list(
            |i| (all_kmers[i], k as u8), 
            |j| all_dummies.get(j), 
            merged_range.end, 
            all_kmers.len(), 
            all_dummies.len()
        );

        Self {
            all_dummies, all_kmers, 
            dummy_idx: dummy_start,
            kmer_idx: kmer_start,
            dummy_range: dummy_start..dummy_end,
            kmer_range: kmer_start..kmer_end,
            k
        }

    }
}

// Assumes the two sorted lists have no duplicates in them or between them.
// Returns (pos_a, pos_b, take_from_a)
fn binary_search_position_in_merged_list<T: PartialOrd + Eq, Access1: Fn(usize) -> T, Access2: Fn(usize) -> T>(access_a: Access1, access_b: Access2, target_pos: usize, len_a: usize, len_b: usize) -> (usize, usize, bool) {

    assert!(target_pos <= len_a + len_b); // One-past the end allowed 

    let pos_in_merged_list = |a_idx: usize| {
        // Returns the index of a[a_idx] in the merged list of a and b
        // That is equal to to number of element in a that are smaller than a[a_idx],
        // and the number of elements in b that are smaller than a[a_idx].
        // We imagine that a[a_len] is infinity to make that corner case work.

        if a_idx == len_a {
            // All of a and b are smaller than infinity
            len_a + len_b 
        } else {
            let b_count = binary_search_leftmost_that_fulfills_pred(|j| j, |b_idx| access_b(b_idx) > access_a(a_idx), len_b);
            a_idx + b_count
        }

    };

    // Find the smallest a_idx that is at or after the target position
    let a_idx = binary_search_leftmost_that_fulfills_pred(|i| i, |a_idx| pos_in_merged_list(a_idx) >= target_pos, len_a);
    if pos_in_merged_list(a_idx) == target_pos {
        // Bingo. Find the b_idx of the next element in b, that is
        // the smallest element in b that is larger than a[a_idx]
        if a_idx == len_a {
            (len_a, len_b, true) // End
        } else {
            let x = access_a(a_idx); 
            let b_idx = binary_search_leftmost_that_fulfills_pred(|j| j, |b_idx| access_b(b_idx) > x, len_b);
            assert_eq!(a_idx + b_idx, target_pos);
            (a_idx, b_idx, true)
        }
    } else {
        // a[a_idx] is after the target position, and a[a_idx-1] is before it (if exists)
        // This means that there are a_idx elements from a before the target position.
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
    use rand::{seq::SliceRandom, Rng, SeedableRng};

    use super::binary_search_position_in_merged_list;

    #[test]
    fn test_binary_search_merged_list() {
        // Test empty vs empty
        assert_eq!((0,0,true), binary_search_position_in_merged_list(|i| i, |j| j, 0, 0, 0));

        // Random testcases
        let seed = 1234;
        let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
        for rep in 0..100 { // Stress test with 100 random runs
            let mut v: Vec<usize> = (0..20).collect();
            v.shuffle(&mut rng);
            let split_point = if rep == 0 {
                0 // Test empty v1    
            } else if rep == 1 {
                v.len() // Test empty v2
            } else {
                // Random split point
                rng.gen_range(0,v.len()+1)
            };

            let mut v1 = v[0..split_point].to_vec();
            let mut v2 = v[split_point..].to_vec();
            v1.sort();
            v2.sort();

            let mut merged: Vec<(usize, usize, bool)> = vec![]; // (i_v1, i_v2, b). b tells which vector it's from
            merged.extend(v1.iter().enumerate().map(|x| (*x.1, x.0, true)));
            merged.extend(v2.iter().enumerate().map(|x| (*x.1, x.0, false)));
            merged.sort();
            let n_merged = merged.len();
            merged.push((v1.len(), v2.len(), true)); // One past the end. The bit is true to match how the search treats this.

            for query in 0..=n_merged {
                let mut true_i = 0;
                let mut true_j = 0;
                for (_,_,from_a) in &merged[0..query] {
                    true_i += *from_a as usize;
                    true_j += !(*from_a) as usize;
                }
                let true_from_a = merged[query].2;
                let (i,j,from_a) = binary_search_position_in_merged_list(|i| v1[i], |j| v2[j], query, v1.len(), v2.len());
                assert_eq!((i,j,from_a), (true_i, true_j, true_from_a));
            }
        }

    }
//fn binary_search_position_in_merged_list<T: PartialOrd + Eq, Access: Fn(usize) -> T>(access_a: Access, access_b: Access, target_pos: usize, len_a: usize, len_b: usize) -> (usize, usize, bool) {
}
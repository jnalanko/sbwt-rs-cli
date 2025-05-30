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
        let xc = x.set_from_left(k-1, 0).right_shift(1).set_from_left(0, c);

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
        let start = find_first_starting_with::<B>(sorted_kmers, c as u8);
        let end = find_first_starting_with::<B>(sorted_kmers, c as u8 + 1);
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


    log::info!("Contructing dummy k-mers");
    let iterable = BitVector::from(has_predecessor);
    let mut required_dummies: Vec::<(LongKmer::<B>, u8)> = iterable.zero_iter().par_bridge().map(|x| {
        let mut prefix = sorted_kmers[x.1];
        (0..k).collect::<Vec<usize>>().iter().map(|i| {
            let len = k - i - 1;
            prefix = prefix.left_shift(1);
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

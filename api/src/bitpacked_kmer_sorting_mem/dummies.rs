use crate::kmer::LongKmer;
use crate::util::binary_search_leftmost_that_fulfills_pred;

use bitvec::order::Lsb0;
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

pub fn get_has_predecessor_marks<const B: usize>(
    source_slice: &[LongKmer::<B>],
    dest_slice: &[LongKmer<B>],
    k: usize,
    c: u8,
) -> bitvec::vec::BitVec<u64, Lsb0> {
    let mut bits = bitvec::vec::BitVec::<u64,Lsb0>::new(); 
    bits.resize(dest_slice.len(), false);
    let mut pointed_idx = 0;
    source_slice.iter().for_each(|x| {
        let xc = x.copy_set_from_left(k-1, 0).right_shifted(1).copy_set_from_left(0, c);

        while pointed_idx < dest_slice.len() {
            match dest_slice[pointed_idx].cmp(&xc) {
                std::cmp::Ordering::Greater => {
                    break
                },
                std::cmp::Ordering::Equal => {
                    bits.set(pointed_idx, true);
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
    sigma: usize, k: usize, n_threads: usize
) -> Vec<(LongKmer<B>, u8)> {

    let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();
    thread_pool.install(||{

        let mut rows = Vec::<bitvec::vec::BitVec::<u64, Lsb0>>::new();

        log::info!("Identifying k-mers without predecessors");
        for c in 0..sigma {
            let input_ranges = crate::util::segment_range(0..sorted_kmers.len(), n_threads);
            let pieces = input_ranges.into_par_iter().map(|range| {
                if !range.is_empty() {
                    let source_slice = &sorted_kmers[range.clone()];
                    let x_start = sorted_kmers[range.start];
                    let cx_start = x_start.copy_set_from_left(k-1, 0).right_shifted(1).copy_set_from_left(0, c as u8);
                    let dest_slice_start = match sorted_kmers.binary_search_by(|probe| probe.cmp(&cx_start)) {
                        Ok(y) => y,
                        Err(y) => y,
                    };
                    let dest_slice_end = if range.end < sorted_kmers.len() {
                        let x_end = sorted_kmers[range.end]; 
                        let cx_end = x_end.copy_set_from_left(k-1, 0).right_shifted(1).copy_set_from_left(0, c as u8);
                        match sorted_kmers.binary_search_by(|probe| probe.cmp(&cx_end)) {
                            Ok(y) => y,
                            Err(y) => y,
                        }
                    } else {
                        sorted_kmers.len()
                    };

                    get_has_predecessor_marks(source_slice, &sorted_kmers[dest_slice_start..dest_slice_end], k, c as u8) 
                } else {
                    bitvec::vec::BitVec::<u64,Lsb0>::new()
                }
            }).collect();
            let row = crate::util::parallel_bitvec_concat(pieces);
            rows.push(row);
        }

        // Bitwise-or the rows together
        let mut rows_iter = rows.into_iter();
        let mut row = rows_iter.next().unwrap();
        for other in rows_iter {
            assert_eq!(row.len(), other.len());
            row |= other;
        }

        log::info!("Constructing dummy k-mers");
        let mut required_dummies: Vec::<(LongKmer::<B>, u8)> = row.iter_zeros().par_bridge().map(|x| {
            let mut prefix = sorted_kmers[x];
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
    })
}
use crate::bitpacked_kmer_sorting_mem::cursors::find_in_nondummy;
use crate::kmer::LongKmer;

use simple_sds_sbwt::ops::Select;
use simple_sds_sbwt::ops::SelectZero;
use simple_sds_sbwt::bit_vector::BitVector;
use simple_sds_sbwt::raw_vector::*;
use rayon::prelude::*;

pub fn get_set_bits<const B: usize>(
    kmers: &[LongKmer::<B>],
    cursor: &mut (std::io::Cursor<&[LongKmer<B>]>, usize),
    k: usize,
    c: u8,
) -> simple_sds_sbwt::raw_vector::RawVector {
    let mut bits = simple_sds_sbwt::raw_vector::RawVector::new();
    bits.resize(kmers.len(), false);
    // let mut set_bits: Vec<usize> = Vec::new();
    let mut pointed_idx = 0;
    kmers.iter().for_each(|x| {
        let xc = x.set_from_left(k-1, 0).right_shift(1).set_from_left(0, c as u8);

        while pointed_idx < cursor.0.get_ref().len() {
            match cursor.0.get_ref()[pointed_idx].cmp(&xc) {
                std::cmp::Ordering::Greater => {
                    break
                },
                std::cmp::Ordering::Equal => {
                    // set_bits.push(cursor.nondummy_position());
                    bits.set_bit(pointed_idx + cursor.1, true);
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

    let mut sorted_kmers_cursor = std::io::Cursor::new(sorted_kmers);

    let mut char_cursors: Vec<(std::io::Cursor::<&[LongKmer::<B>]>, usize)> = (0..sigma).map(|c|{
        let pos = find_in_nondummy::<B>(&mut sorted_kmers_cursor, c as u8);
        let start = pos as usize;
        let end = if c < sigma - 1 {
            find_in_nondummy::<B>(&mut sorted_kmers_cursor, c as u8 + 1)
        } else {
            sorted_kmers.len() as u64
        };
        (std::io::Cursor::new(&sorted_kmers[start..(end as usize)]), pos as usize)
    }).collect();

    let has_predecessor = char_cursors.par_iter_mut().enumerate().map(|(c, cursor)| {
        get_set_bits(sorted_kmers, cursor, k, c as u8)
    }).reduce(|| {
        let mut res = simple_sds_sbwt::raw_vector::RawVector::new();
        res.resize(n, false);
        res
    }, |mut a, b| {
        let bv = BitVector::from(b);
        bv.one_iter().for_each(|idx| a.set_bit(idx.1, true));
        a
    });

    let iterable = BitVector::from(has_predecessor);
    let mut required_dummies: Vec::<(LongKmer::<B>, u8)> = iterable.zero_iter().par_bridge().map(|x| {
        let mut prefix = sorted_kmers[x.1];
        (0..k).collect::<Vec<usize>>().iter().map(|i| {
            let len = k - i - 1;
            prefix = prefix.left_shift(1);
            (prefix, len as u8)
        }).collect::<Vec<(LongKmer::<B>, u8)>>()
    }).flatten().collect();

    // We always assume that the empty k-mer exists. This assumption is reflected in the C-arrya
    // later, which adds one "ghost dollar" count to all counts.
    required_dummies.push((LongKmer::<B>::from_ascii(b"").unwrap(), 0));

    required_dummies.par_sort_unstable();
    required_dummies.dedup();
    required_dummies.shrink_to_fit();

    required_dummies
}

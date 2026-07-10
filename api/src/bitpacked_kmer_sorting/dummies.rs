use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::ops::Range;
use std::path::Path;

use crate::bitpacked_kmer_sorting::disk_access::dummy_record_len;
use crate::kmer::LongKmer;
use crate::util::binary_search_leftmost_that_fulfills_pred;

use bitvec::order::Lsb0;
use byteorder::ReadBytesExt;
use rayon::prelude::*;

use super::disk_access;

/// Returns the index of the first k-mer in the sorted k-mer file that starts with character
/// `c`, or `n_kmers` if there is none. Since the file is globally sorted, all k-mers starting
/// with a given first character form one contiguous block.
fn find_first_kmer_starting_with<const B: usize>(path: &Path, n_kmers: usize, c: u8) -> usize {
    binary_search_leftmost_that_fulfills_pred(
        |i| disk_access::read_kmer_at::<B>(path, i),
        |kmer: LongKmer<B>| kmer.get_from_left(0) >= c,
        n_kmers,
    )
}

/// Disk-streaming analogue of the in-memory `get_has_predecessor_marks`: for every k-mer `x`
/// read (in order) from `source_range`, computes its `c`-successor `xc` and marks the
/// corresponding position (relative to `dest_range`) in the output bitvector if `xc` is found
/// there. Both ranges are read purely sequentially, seeked only once at the start.
fn get_has_predecessor_marks_disk<const B: usize>(
    path: &Path,
    source_range: Range<usize>,
    dest_range: Range<usize>,
    k: usize,
    c: u8,
) -> bitvec::vec::BitVec<u64, Lsb0> {
    let mut bits = bitvec::vec::BitVec::<u64, Lsb0>::new();
    bits.resize(dest_range.len(), false);

    if source_range.is_empty() {
        return bits;
    }

    let mut source_reader = disk_access::seek_kmer_reader::<B>(path, source_range.start);
    let mut dest_reader = disk_access::seek_kmer_reader::<B>(path, dest_range.start);
    let mut pointed_idx = 0_usize;
    let mut dest_lookahead: Option<LongKmer<B>> = if dest_range.is_empty() {
        None
    } else {
        Some(LongKmer::<B>::load(&mut dest_reader).unwrap().unwrap())
    };

    for _ in 0..source_range.len() {
        let x = LongKmer::<B>::load(&mut source_reader).unwrap().unwrap();
        let xc = x.copy_set_from_left(k - 1, 0).right_shifted(1).copy_set_from_left(0, c);

        while pointed_idx < dest_range.len() {
            let dest_kmer = dest_lookahead.unwrap();
            match dest_kmer.cmp(&xc) {
                std::cmp::Ordering::Greater => {
                    break
                },
                std::cmp::Ordering::Equal => {
                    bits.set(pointed_idx, true);
                    pointed_idx += 1;
                    dest_lookahead = if pointed_idx < dest_range.len() {
                        Some(LongKmer::<B>::load(&mut dest_reader).unwrap().unwrap())
                    } else {
                        None
                    };
                    break
                },
                std::cmp::Ordering::Less => {
                    pointed_idx += 1;
                    dest_lookahead = if pointed_idx < dest_range.len() {
                        Some(LongKmer::<B>::load(&mut dest_reader).unwrap().unwrap())
                    } else {
                        None
                    };
                    // no break
                }
            }
        }
    }

    bits
}

/// Phase A: identifies, for every k-mer in the sorted file, whether it has an incoming
/// predecessor edge. Parallelized by segmenting the k-mer index range per character (mirrors
/// `bitpacked_kmer_sorting_mem::dummies::get_sorted_dummies`'s Phase A).
fn get_disk_has_predecessor_bits<const B: usize>(
    path: &Path,
    n_kmers: usize,
    sigma: usize,
    k: usize,
    n_threads: usize,
) -> bitvec::vec::BitVec<u64, Lsb0> {
    let char_ranges: Vec<Range<usize>> = (0..sigma).map(|c| {
        let start = find_first_kmer_starting_with::<B>(path, n_kmers, c as u8);
        let end = find_first_kmer_starting_with::<B>(path, n_kmers, c as u8 + 1);
        start..end
    }).collect();

    let mut pieces = Vec::<bitvec::vec::BitVec::<u64, Lsb0>>::new();
    for c in 0..sigma {
        let input_ranges = crate::util::segment_range(0..n_kmers, n_threads);
        let char_pieces: Vec<bitvec::vec::BitVec::<u64, Lsb0>> = input_ranges.into_par_iter().map(|range| {
            if !range.is_empty() {
                let x_start = disk_access::read_kmer_at::<B>(path, range.start);
                let cx_start = x_start.copy_set_from_left(k - 1, 0).right_shifted(1).copy_set_from_left(0, c as u8);
                let dest_start = if range.start == 0 {
                    // Make sure the output bitvector slice covers the start even if the first k-mers
                    // do not have predecessors.
                    char_ranges[c].start
                } else {
                    binary_search_leftmost_that_fulfills_pred(
                        |i| disk_access::read_kmer_at::<B>(path, i),
                        |y: LongKmer<B>| y >= cx_start,
                        n_kmers,
                    )
                };
                let dest_end = if range.end < n_kmers {
                    let x_end = disk_access::read_kmer_at::<B>(path, range.end);
                    let cx_end = x_end.copy_set_from_left(k - 1, 0).right_shifted(1).copy_set_from_left(0, c as u8);
                    binary_search_leftmost_that_fulfills_pred(
                        |i| disk_access::read_kmer_at::<B>(path, i),
                        |y: LongKmer<B>| y >= cx_end,
                        n_kmers,
                    )
                } else {
                    char_ranges[c].end
                };

                get_has_predecessor_marks_disk::<B>(path, range, dest_start..dest_end, k, c as u8)
            } else {
                bitvec::vec::BitVec::<u64, Lsb0>::new()
            }
        }).collect();

        pieces.extend(char_pieces.into_iter());
    }

    crate::util::parallel_bitvec_concat(pieces)
}

/// Phase B: for every k-mer without a predecessor, emits the k dummy prefixes of decreasing
/// length. Parallelized by segmenting the k-mer index range and streaming sequentially through
/// each segment.
fn build_dummy_prefixes_disk<const B: usize>(
    path: &Path,
    has_predecessor: &bitvec::vec::BitVec<u64, Lsb0>,
    n_kmers: usize,
    k: usize,
    n_threads: usize,
) -> Vec<(LongKmer<B>, u8)> {
    let colex_pieces = crate::util::segment_range(0..n_kmers, n_threads);
    let required_dummies_pieces: Vec<Vec<(LongKmer<B>, u8)>> = colex_pieces.into_par_iter().map(|range| {
        let mut dummies = Vec::<(LongKmer<B>, u8)>::new();
        if range.is_empty() {
            return dummies;
        }

        let mut reader = disk_access::seek_kmer_reader::<B>(path, range.start);
        for rel in 0..range.len() {
            let kmer = LongKmer::<B>::load(&mut reader).unwrap().unwrap();
            if !has_predecessor[range.start + rel] {
                let mut prefix = kmer;
                for i in 0..k {
                    let len = k - i - 1;
                    prefix = prefix.left_shifted(1);
                    dummies.push((prefix, len as u8))
                }
            }
        }
        dummies
    }).collect();

    // Concatenate pieces
    required_dummies_pieces.into_iter().fold(vec![], |mut acc, v| { acc.extend(v); acc })
}

fn file_size(path: &Path) -> usize {
    std::fs::metadata(path).unwrap().len() as usize
}

pub fn get_sorted_dummies<const B: usize>(
    sorted_kmers_filepath: &Path,
    n_kmers: usize,
    sigma: usize,
    k: usize,
    n_threads: usize,
    given_mers: Option<&Path>, // All dummy suffixes of these will be included.
) -> Vec<(LongKmer::<B>, u8)>{

    let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();
    thread_pool.install(|| {

        log::info!("Identifying k-mers without predecessors");
        let has_predecessor = get_disk_has_predecessor_bits::<B>(sorted_kmers_filepath, n_kmers, sigma, k, n_threads);

        log::info!("Building dummy k-mer strings");
        let mut required_dummies = build_dummy_prefixes_disk::<B>(sorted_kmers_filepath, &has_predecessor, n_kmers, k, n_threads);

        if let Some(given_mers) = given_mers {
            log::info!("Adding extra dummy k-mer strings");
            let dummy_record_len = crate::bitpacked_kmer_sorting::disk_access::dummy_record_len::<B>();
            assert!(file_size(given_mers) % dummy_record_len == 0);
            let n_given_mers = file_size(given_mers) / dummy_record_len;
            let mut first_mers_reader = BufReader::new(File::open(given_mers).unwrap());
            for _ in 0..n_given_mers {
                let mut prefix = LongKmer::<B>::load(&mut first_mers_reader).unwrap().unwrap();
                let mut prefix_full_len = first_mers_reader.read_u8().unwrap() as usize;

                if prefix_full_len as usize == k {
                    prefix = prefix.left_shifted(1);
                    prefix_full_len -= 1;
                }
                for i in 0..prefix_full_len {
                    let len = prefix_full_len - i;
                    required_dummies.push((prefix, len as u8));
                    prefix = prefix.left_shifted(1);
                }
            }
        }

        // We always assume that the empty k-mer exists. This assumption is reflected in the C-array
        // later, which adds one "ghost dollar" count to all counts.
        required_dummies.push((LongKmer::<B>::from_ascii(b"").unwrap(), 0));

        log::info!("Sorting and deduplicating dummies");
        required_dummies.par_sort_unstable();
        required_dummies.dedup();
        required_dummies.shrink_to_fit();
        required_dummies
    })
}

// Returns the number of bytes written
pub fn write_to_disk<const B: usize>(dummies: Vec<(LongKmer::<B>, u8)>, writer: &mut std::fs::File){
    let mut bw = BufWriter::new(writer);
    for (kmer, len) in dummies.iter(){
        kmer.serialize(&mut bw).unwrap();
        bw.write_all(&[*len]).unwrap();
    }
    bw.flush().unwrap();
}

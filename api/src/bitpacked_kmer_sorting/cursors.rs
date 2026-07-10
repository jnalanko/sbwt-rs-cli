use std::{fs::File, io::BufReader, ops::Range, path::Path};

use bitvec::order::Lsb0;
use rayon::iter::IntoParallelIterator;
use rayon::iter::ParallelIterator;
use simple_sds_sbwt::ops::Push;

use crate::kmer::LongKmer;
use crate::util::binary_search_leftmost_that_fulfills_pred;

use super::disk_access;

pub struct DummyNodeMerger<R: std::io::Read, const B: usize> {
    dummy_reader: R, // Stream of k-mer objects
    nondummy_reader: R, // Stream of pairs (kmer, len)

    dummy_kmer: Option<(LongKmer::<B>, u8)>,
    nondummy_kmer: Option<(LongKmer::<B>, u8)>,

    k: usize,

    dummy_position: usize, // Position of the dummy cursor
    nondummy_position: usize, // Position of the nondummy cursor

    // The cursor treats a side as exhausted (returns None for it) once its position counter
    // reaches the corresponding range end, even if the underlying file has more data past that
    // point. Unbounded mergers (see `new`/`new_with_initial_positions`) use usize::MAX here, so
    // they only stop at the real end of file, matching the previous unbounded behavior exactly.
    dummy_range_end: usize,
    nondummy_range_end: usize,
}

impl <R: std::io::Read, const B: usize> DummyNodeMerger<R, B> {

    pub fn read_from_dummy_reader(dummy_reader: &mut R) -> Option<(LongKmer::<B>, u8)>{
        let kmer = match LongKmer::<B>::load(dummy_reader){
            Ok(kmer_opt) => {
                #[allow(clippy::question_mark)] // More space to write comments on the cases
                match kmer_opt{
                    Some(kmer) => kmer,
                    None => return None, // End of stream
                }
            },
            Err(e) => panic!("IO error while streaming sorted k-mers: {}", e),
        };

        // Read length
        let mut buf = [0_u8; 1];
        dummy_reader.read_exact(&mut buf).unwrap();
        let len = u8::from_le_bytes(buf);

        Some((kmer, len))
    }

    pub fn read_from_non_dummy_reader(nondummy_reader: &mut R, k: usize) -> Option<(LongKmer::<B>, u8)>{
        match LongKmer::<B>::load(nondummy_reader){
            Ok(kmer_opt) => {
                kmer_opt.map(|kmer| (kmer, k as u8))
            },
            Err(e) => panic!("IO error while streaming sorted k-mers: {}", e),
        }
    }

    #[allow(dead_code)] // Unbounded constructor; kept for API completeness and used in tests.
    pub fn new(mut dummy_reader: R, mut nondummy_reader: R, k: usize) -> Self {
        let dummy_kmer = Self::read_from_dummy_reader(&mut dummy_reader);
        let nondummy_kmer = Self::read_from_non_dummy_reader(&mut nondummy_reader, k);

        Self {
            dummy_reader,
            nondummy_reader,
            dummy_kmer,
            nondummy_kmer,
            k,
            dummy_position: 0,
            nondummy_position: 0,
            dummy_range_end: usize::MAX,
            nondummy_range_end: usize::MAX,
        }
    }

    // TODO: this is stupid. The given positions are just for bookkeeping that the caller might use. They don't affect the
    // cursor at all. Need to refactor.
    #[allow(dead_code)] // Unbounded constructor; kept for API completeness.
    pub fn new_with_initial_positions(mut dummy_reader: R, mut nondummy_reader: R, k: usize, dummy_position: usize, nondummy_position: usize) -> Self {
        let dummy_kmer = Self::read_from_dummy_reader(&mut dummy_reader);
        let nondummy_kmer = Self::read_from_non_dummy_reader(&mut nondummy_reader, k);

        Self {
            dummy_reader,
            nondummy_reader,
            dummy_kmer,
            nondummy_kmer,
            k,
            dummy_position,
            nondummy_position,
            dummy_range_end: usize::MAX,
            nondummy_range_end: usize::MAX,
        }
    }

    /// Constructs a merger whose two sides act exhausted once their position counters reach
    /// `dummy_range_end`/`nondummy_range_end`, even if the underlying file has more data past
    /// that point. The readers must already be seeked to `dummy_position`/`nondummy_position`.
    pub fn new_range_bounded(
        mut dummy_reader: R, mut nondummy_reader: R, k: usize,
        dummy_position: usize, dummy_range_end: usize,
        nondummy_position: usize, nondummy_range_end: usize,
    ) -> Self {
        let dummy_kmer = if dummy_position < dummy_range_end {
            Self::read_from_dummy_reader(&mut dummy_reader)
        } else {
            None
        };
        let nondummy_kmer = if nondummy_position < nondummy_range_end {
            Self::read_from_non_dummy_reader(&mut nondummy_reader, k)
        } else {
            None
        };

        Self {
            dummy_reader,
            nondummy_reader,
            dummy_kmer,
            nondummy_kmer,
            k,
            dummy_position,
            nondummy_position,
            dummy_range_end,
            nondummy_range_end,
        }
    }

    pub fn peek(&self) -> Option<(LongKmer::<B>, u8)>{
        match (self.dummy_kmer, self.nondummy_kmer){
            (None, None) => None,
            (Some(dummy_kmer), None) => {
                Some(dummy_kmer)
            },
            (None, Some(nondummy_kmer)) => {
                Some(nondummy_kmer)
            },
            (Some(dummy_kmer), Some(nondummy_kmer)) => {
                if dummy_kmer < nondummy_kmer {
                    Some(dummy_kmer)
                } else {
                    Some(nondummy_kmer)
                }
            }
        }
    }

    #[allow(dead_code)]
    pub fn dummy_position(&self) -> usize{
        self.dummy_position
    }

    #[allow(dead_code)]
    pub fn nondummy_position(&self)  -> usize{
        self.nondummy_position
    }

    /// The index of the next element to be returned, in the conceptual merged (dummy + nondummy) order.
    pub fn cur_merged_index(&self) -> usize {
        self.dummy_position + self.nondummy_position
    }
}

impl<const B: usize> Iterator for DummyNodeMerger<BufReader<File>, B> {
    type Item = (LongKmer<B>, u8);

    // Produces pairs (kmer, length)
    fn next(&mut self) -> Option<(LongKmer::<B>, u8)> {
        match (self.dummy_kmer, self.nondummy_kmer){
            (None, None) => None,
            (Some(dummy_kmer), None) => {
                self.dummy_position += 1;
                self.dummy_kmer = if self.dummy_position < self.dummy_range_end {
                    Self::read_from_dummy_reader(&mut self.dummy_reader)
                } else {
                    None
                };
                Some(dummy_kmer)
            },
            (None, Some(nondummy_kmer)) => {
                self.nondummy_position += 1;
                self.nondummy_kmer = if self.nondummy_position < self.nondummy_range_end {
                    Self::read_from_non_dummy_reader(&mut self.nondummy_reader, self.k)
                } else {
                    None
                };
                Some(nondummy_kmer)
            },
            (Some(dummy_kmer), Some(nondummy_kmer)) => {
                if dummy_kmer < nondummy_kmer {
                    self.dummy_position += 1;
                    self.dummy_kmer = if self.dummy_position < self.dummy_range_end {
                        Self::read_from_dummy_reader(&mut self.dummy_reader)
                    } else {
                        None
                    };
                    Some(dummy_kmer)
                } else {
                    self.nondummy_position += 1;
                    self.nondummy_kmer = if self.nondummy_position < self.nondummy_range_end {
                        Self::read_from_non_dummy_reader(&mut self.nondummy_reader, self.k)
                    } else {
                        None
                    };
                    Some(nondummy_kmer)
                }
            }
        }
    }

}

impl<const B: usize> DummyNodeMerger<BufReader<File>, B> {
    /// Builds a merger scoped to `merged_range` within the conceptual merged (dummy + nondummy)
    /// sorted sequence. Both endpoints of `merged_range` are located via binary search over the
    /// two sorted disk files (`n_dummies`/`n_kmers` records respectively), and the two readers
    /// are seeked to the resulting split point. `merged_range.end` may be `n_dummies + n_kmers`
    /// (one past the very end), in which case the dest side is effectively unbounded.
    pub fn new_bounded(
        dummy_path: &Path,
        nondummy_path: &Path,
        k: usize,
        merged_range: Range<usize>,
        n_dummies: usize,
        n_kmers: usize,
    ) -> Self {
        let (dummy_start, nondummy_start, _) = crate::util::binary_search_position_in_merged_list(
            |j| disk_access::read_dummy_at::<B>(dummy_path, j),
            |i| (disk_access::read_kmer_at::<B>(nondummy_path, i), k as u8),
            merged_range.start, n_dummies, n_kmers,
        );
        let (dummy_end, nondummy_end, _) = crate::util::binary_search_position_in_merged_list(
            |j| disk_access::read_dummy_at::<B>(dummy_path, j),
            |i| (disk_access::read_kmer_at::<B>(nondummy_path, i), k as u8),
            merged_range.end, n_dummies, n_kmers,
        );

        let dummy_reader = disk_access::seek_dummy_reader::<B>(dummy_path, dummy_start);
        let nondummy_reader = disk_access::seek_kmer_reader::<B>(nondummy_path, nondummy_start);

        Self::new_range_bounded(dummy_reader, nondummy_reader, k, dummy_start, dummy_end, nondummy_start, nondummy_end)
    }
}

/// Disk analogue of the in-memory `get_ith_merged_kmer`: returns the element at position `i`
/// (0-indexed) in the conceptual merged (dummy + nondummy) sorted sequence. O(log^2(n)) time.
pub fn get_ith_merged_kmer_disk<const B: usize>(
    dummy_path: &Path,
    nondummy_path: &Path,
    i: usize,
    k: usize,
    n_dummies: usize,
    n_kmers: usize,
) -> (LongKmer<B>, u8) {
    let (dummy_idx, nondummy_idx, in_dummies) = crate::util::binary_search_position_in_merged_list(
        |j| disk_access::read_dummy_at::<B>(dummy_path, j),
        |i| (disk_access::read_kmer_at::<B>(nondummy_path, i), k as u8),
        i, n_dummies, n_kmers,
    );

    if in_dummies {
        disk_access::read_dummy_at::<B>(dummy_path, dummy_idx)
    } else {
        (disk_access::read_kmer_at::<B>(nondummy_path, nondummy_idx), k as u8)
    }
}

// Here c is from 0..3. Returns the k-mer (or dummy) reached by walking backward from `kmer`
// along a `c`-labeled edge: `c` prepended to `kmer`'s first (len-1) characters.
fn prepend_c<const B: usize>(kmer: (LongKmer<B>, u8), k: usize, c: usize) -> (LongKmer<B>, u8) {
    if kmer.1 as usize == k {
        (
            kmer.0
                .copy_set_from_left(k - 1, 0)
                .right_shifted(1)
                .copy_set_from_left(0, c as u8),
            k as u8,
        )
    } else {
        (kmer.0.right_shifted(1).copy_set_from_left(0, c as u8), kmer.1 + 1) // Dummy
    }
}

/// Disk analogue of the in-memory `build_lcs_array`. Computes the (compressed) LCS array by
/// segmenting the merged index range and, per segment, streaming forward through a
/// range-bounded `DummyNodeMerger` while computing pairwise LCP values against the previous
/// merged element (looked up once per segment via `get_ith_merged_kmer_disk`).
pub fn build_lcs_array_disk<const B: usize>(
    dummy_filepath: &Path,
    nondummy_filepath: &Path,
    n_dummies: usize,
    n_kmers: usize,
    k: usize,
    n_threads: usize,
) -> simple_sds_sbwt::int_vector::IntVector {
    // LCS values are between 0 and k-1
    assert!(k > 0);
    assert!(k < u16::MAX as usize);

    let n = n_kmers + n_dummies;

    log::info!("Computing LCS values");

    // We start the segmentation from 1 so that we always have a previous k-mer to compare against
    let segments = crate::util::segment_range(1..n, n_threads);
    let lcs_pieces: Vec<Vec<u16>> = segments.into_par_iter().map(|range| {
        let (mut prev_kmer, mut prev_len) = get_ith_merged_kmer_disk::<B>(dummy_filepath, nondummy_filepath, range.start - 1, k, n_dummies, n_kmers); // range.start >= 1 so this is ok
        let mut subrange = DummyNodeMerger::<BufReader<File>, B>::new_bounded(dummy_filepath, nondummy_filepath, k, range.clone(), n_dummies, n_kmers);
        let mut lcs_piece = Vec::<u16>::with_capacity(range.len());
        while let Some((kmer, len)) = subrange.next() {
            let lcp_value = LongKmer::<B>::lcp_with_different_lengths((&prev_kmer, prev_len), (&kmer, len));
            lcs_piece.push(lcp_value as u16);
            (prev_kmer, prev_len) = (kmer, len);
        }
        lcs_piece
    }).collect();

    // Compress into log(k) bits per element
    log::info!("Compressing LCS array to log(k) bits per element");
    let bitwidth = 64 - (k as u64 - 1).leading_zeros();
    let mut compressed_lcs = simple_sds_sbwt::int_vector::IntVector::with_capacity(n, bitwidth as usize).unwrap();
    compressed_lcs.push(0); // lcs[0] = 0 by definition
    for piece in lcs_pieces { // Concatenate pieces. Todo: could this be done in parallel?
        compressed_lcs.extend(piece);
    }

    compressed_lcs
}

// Returns the SBWT bit vectors and optionally the LCS array
pub fn build_sbwt_bit_vectors<const B: usize>(
    dummy_filepath: &Path,
    nondummy_filepath: &Path,
    n_dummies: usize,
    n_kmers: usize,
    k: usize,
    sigma: usize,
    build_lcs: bool,
    n_threads: usize,
) -> (Vec<bitvec::vec::BitVec::<u64, Lsb0>>, Option<simple_sds_sbwt::int_vector::IntVector>)
{
    let n = n_kmers + n_dummies;

    let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();
    thread_pool.install(|| {
        // Split the merged range 0..n into segments for threads
        let input_ranges = crate::util::segment_range(0..n, n_threads);

        let mut rows = vec![];
        for c in 0..sigma {
            let row_pieces: Vec<bitvec::vec::BitVec<u64, Lsb0>> = input_ranges.clone().into_par_iter().map(|input_range| {
                let mut row_piece: bitvec::vec::BitVec::<u64, Lsb0> = bitvec::vec::BitVec::with_capacity(input_range.len());
                row_piece.resize(input_range.len(), false);

                if !input_range.is_empty() {
                    // Let x be the first input kmer in this range. The destinations of the c-edges start
                    // from the smallest k-mer that is larger or equal to cx.
                    let x: (LongKmer<B>, u8) = get_ith_merged_kmer_disk::<B>(dummy_filepath, nondummy_filepath, input_range.start, k, n_dummies, n_kmers);
                    let (cx, cx_len) = prepend_c(x, k, c);

                    let cx_dummy_insertion_index = binary_search_leftmost_that_fulfills_pred(
                        |j| disk_access::read_dummy_at::<B>(dummy_filepath, j),
                        |y: (LongKmer<B>, u8)| y >= (cx, cx_len),
                        n_dummies,
                    );
                    let cx_kmer_insertion_index = binary_search_leftmost_that_fulfills_pred(
                        |i| disk_access::read_kmer_at::<B>(nondummy_filepath, i),
                        |y: LongKmer<B>| (y, k as u8) >= (cx, cx_len),
                        n_kmers,
                    );
                    let dest_start_idx = cx_dummy_insertion_index + cx_kmer_insertion_index;

                    let mut src_pointer = DummyNodeMerger::<BufReader<File>, B>::new_bounded(dummy_filepath, nondummy_filepath, k, input_range.clone(), n_dummies, n_kmers); // Origin of edge
                    let mut dest_pointer = DummyNodeMerger::<BufReader<File>, B>::new_bounded(dummy_filepath, nondummy_filepath, k, dest_start_idx..n, n_dummies, n_kmers); // Destination of edge

                    // We might be starting at a k-mer that is not a suffix group leader. We must not
                    // add any edges for it. Rewind forward until we are at a suffix group leader.
                    let mut cur = x;
                    if input_range.start > 0 {
                        let mut prev = get_ith_merged_kmer_disk::<B>(dummy_filepath, nondummy_filepath, input_range.start - 1, k, n_dummies, n_kmers);
                        while LongKmer::<B>::lcp_with_different_lengths((&prev.0, prev.1), (&cur.0, cur.1)) == k-1 {
                            src_pointer.next(); // Advance
                            prev = cur;
                            if let Some(next) = src_pointer.peek() {
                                cur = next;
                            } else {
                                break;
                            }
                        }
                    }

                    // Now src_pointer should point to the first leader of a suffix group
                    // in the input range.
                    while let Some((kmer, len)) = src_pointer.next() {
                        let kmer_c = prepend_c((kmer,len), k, c);

                        while dest_pointer.peek().is_some_and(|y| y < kmer_c) {
                            dest_pointer.next();
                        }

                        if dest_pointer.peek().is_some_and(|y| y == kmer_c) {
                            row_piece.set(src_pointer.cur_merged_index() - 1 - input_range.start, true); // -1 because we have advanced past the current k-mer
                            dest_pointer.next().unwrap(); // Don't point to this k-mer anymore because it now has an edge
                        }
                    };
                }
                row_piece
            }).collect();
            rows.push(crate::util::parallel_bitvec_concat(row_pieces));
        }

        let lcs = if build_lcs {
            Some(build_lcs_array_disk::<B>(dummy_filepath, nondummy_filepath, n_dummies, n_kmers, k, n_threads))
        } else {
            None
        };

        (rows, lcs)
    })
}

#[cfg(test)]
mod tests{

    use std::io::Write;

    use super::*;

    #[test]
    fn test_new_bounded_matches_sequential_merge(){
        let nondummies = [
            LongKmer::<2>::from_ascii(b"ACGT").unwrap(),
            LongKmer::<2>::from_ascii(b"AGGT").unwrap(),
            LongKmer::<2>::from_ascii(b"GGAA").unwrap(),
            LongKmer::<2>::from_ascii(b"GGGT").unwrap()
        ];
        let dummies = [
            (LongKmer::<2>::from_ascii(b"AAAA").unwrap(),0), // This is actually the empty dummy
            (LongKmer::<2>::from_ascii(b"AAAA").unwrap(),1),
            (LongKmer::<2>::from_ascii(b"ACAA").unwrap(),2),
            (LongKmer::<2>::from_ascii(b"ACAA").unwrap(),3),
            (LongKmer::<2>::from_ascii(b"GGTT").unwrap(),3),
        ];

        let mut temp_file_manager = crate::tempfile::TempFileManager::new(std::path::Path::new("/tmp"));

        let mut nondummy_file = temp_file_manager.create_new_file("test-", 10, ".nondummy");
        let mut dummy_file = temp_file_manager.create_new_file("test-", 10, ".dummy");
        let nondummy_path = nondummy_file.path.clone();
        let dummy_path = dummy_file.path.clone();

        for kmer in nondummies.iter(){
            kmer.serialize(&mut nondummy_file).unwrap();
        }
        for (kmer, len) in dummies.iter(){
            kmer.serialize(&mut dummy_file).unwrap();
            let len_byte = *len as u8;
            dummy_file.write_all(&[len_byte]).unwrap();
        }

        // Flush
        dummy_file.flush().unwrap();
        nondummy_file.flush().unwrap();

        let n_dummies = dummies.len();
        let n_kmers = nondummies.len();
        let n = n_dummies + n_kmers;

        // Build the ground truth via a plain unbounded sequential merge.
        let full_merge = DummyNodeMerger::new(
            std::io::BufReader::new(std::fs::File::open(&dummy_path).unwrap()),
            std::io::BufReader::new(std::fs::File::open(&nondummy_path).unwrap()),
            4,
        );
        let expected: Vec<(LongKmer::<2>, u8)> = full_merge.collect();
        assert_eq!(expected.len(), n);

        // Check get_ith_merged_kmer_disk against the ground truth for every position.
        for i in 0..n {
            assert_eq!(get_ith_merged_kmer_disk::<2>(&dummy_path, &nondummy_path, i, 4, n_dummies, n_kmers), expected[i]);
        }

        // Check that new_bounded, restricted to every possible sub-range, reproduces the
        // corresponding slice of the ground truth merge.
        for start in 0..=n {
            for end in start..=n {
                let bounded = DummyNodeMerger::<BufReader<File>, 2>::new_bounded(&dummy_path, &nondummy_path, 4, start..end, n_dummies, n_kmers);
                let got: Vec<(LongKmer::<2>, u8)> = bounded.collect();
                assert_eq!(got, expected[start..end].to_vec());
            }
        }
    }
}

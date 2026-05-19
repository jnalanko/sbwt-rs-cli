use std::ops::Range;
use bitvec::prelude::*;
use crate::subsetseq::*;
use crate::sbwt::*;

type BitVec = bitvec::vec::BitVec<u64, Lsb0>;

pub(super) fn allocate_rows(sigma: usize, len: usize) -> Vec<BitVec> {
    (0..sigma).map(|_| {
        let mut row = BitVec::new();
        row.resize(len, false);
        row
    }).collect()
}

// Transposes pieces_vecvec[piece][char] into new_rows[char] and concatenates pieces for each char.
pub(super) fn transpose_and_concat_pieces(pieces_vecvec: Vec<Vec<bitvec::vec::BitVec<u64, Lsb0>>>, sigma: usize) -> Vec<bitvec::vec::BitVec<u64, Lsb0>> {
    // Collect pieces for each char ("transpose the Vec<Vec<...>>")
    let mut char_to_piece_list = Vec::<Vec::<bitvec::vec::BitVec::<u64, Lsb0>>>::new();
    for _ in 0..sigma {
        char_to_piece_list.push(vec![]); // Init new piece list for this char
    }
    for char_vecs_for_piece in pieces_vecvec {
        for (c, vec) in char_vecs_for_piece.into_iter().enumerate() {
            char_to_piece_list[c].push(vec);
        }
    }
    // Concatenate pieces for each char
    char_to_piece_list.into_iter().map(crate::util::parallel_bitvec_concat).collect()
}

// Computes the C array, builds the subset rank structure, and assembles the final SbwtIndex.
pub(super) fn build_index<SS: SubsetSeq>(new_rows: Vec<bitvec::vec::BitVec<u64, Lsb0>>, n_kmers: usize, k: usize, n_threads: usize, new_prefix_lookup_table_length: usize) -> SbwtIndex<SS> {
    #[allow(non_snake_case)] // C-array is an established convention in BWT indexes
    let C: Vec<usize> = crate::util::get_C_array_parallel(&new_rows, n_threads);

    let mut subsetseq = SS::new_from_bit_vectors(new_rows);
    subsetseq.build_rank();
    let n_sets = subsetseq.len();
    let mut index = SbwtIndex::<SS>::from_components(
        subsetseq, n_kmers, k, C,
        PrefixLookupTable::new_empty(n_sets));

    let lut = PrefixLookupTable::new(&index, new_prefix_lookup_table_length);
    index.set_lookup_table(lut);
    index
}

/// Popcount of `(a[range] & !b[range])` using word-level operations.
pub(super) fn word_diff_popcount_range(a: &[u64], b: &[u64], range: Range<usize>) -> usize {
    if range.is_empty() { return 0; }
    let start_word = range.start / 64;
    let end_word   = (range.end - 1) / 64;
    let start_bit  = range.start % 64;
    let end_bit    = range.end % 64; // exclusive; 0 means full last word

    if start_word == end_word {
        let mut w = a[start_word] & !b[start_word];
        w &= !0u64 << start_bit;
        if end_bit != 0 { w &= (1u64 << end_bit) - 1; }
        return w.count_ones() as usize;
    }

    let mut count = 0usize;
    count += ((a[start_word] & !b[start_word]) & (!0u64 << start_bit)).count_ones() as usize;
    for w in (start_word + 1)..end_word {
        count += (a[w] & !b[w]).count_ones() as usize;
    }
    let mut w = a[end_word] & !b[end_word];
    if end_bit != 0 { w &= (1u64 << end_bit) - 1; }
    count += w.count_ones() as usize;
    count
}

/// Popcount of `(a[range] & b[range])` using word-level operations.
pub(super) fn word_and_popcount_range(a: &[u64], b: &[u64], range: Range<usize>) -> usize {
    if range.is_empty() { return 0; }
    let start_word = range.start / 64;
    let end_word = (range.end - 1) / 64;
    let start_bit = range.start % 64;
    let end_bit = range.end % 64; // exclusive; 0 means full last word

    if start_word == end_word {
        let mut w = a[start_word] & b[start_word];
        w &= !0u64 << start_bit;
        if end_bit != 0 { w &= (1u64 << end_bit) - 1; }
        return w.count_ones() as usize;
    }

    let mut count = 0usize;
    // First (partial) word
    count += ((a[start_word] & b[start_word]) & (!0u64 << start_bit)).count_ones() as usize;
    // Full middle words
    for w in (start_word + 1)..end_word {
        count += (a[w] & b[w]).count_ones() as usize;
    }
    // Last (partial) word
    let mut w = a[end_word] & b[end_word];
    if end_bit != 0 { w &= (1u64 << end_bit) - 1; }
    count += w.count_ones() as usize;
    count
}

// Code by Martin Kostadinov.

use rayon::slice::ParallelSliceMut;
use simple_sds_sbwt::bit_vector::BitVector;
use simple_sds_sbwt::ops::{BitVec as BitVecTrait, Push, Rank, Select};
use simple_sds_sbwt::raw_vector::{AccessRaw, RawVector};

/// Bounded Burrows-Wheeler Transform.
use super::alternative_construction::input_structures::{
    Bwt, Lcp, CHAR_TO_INDEX
};

type Word = wide::u8x32;
const LANES: usize = Word::LANES as usize;

/// Needs to be executed in a rayon context.
pub fn par_bounded_context_suffix_array(input: &mut Vec<u8>, k: usize) -> Vec<usize> {
    let length = input.len();
    let word_count = k.div_ceil(LANES);

    // Pad the input with '#' characters so that the comparisons need not do bound checks.
    for _ in 1..(word_count * LANES) {
        input.push(b'#');
    }

    let mut suffix_array: Vec<usize> = Vec::<usize>::with_capacity(length);
    suffix_array.extend(0..length);

    suffix_array.par_sort_by(|suffix_a, suffix_b| {
        let mut cursor_a = *suffix_a;
        let mut cursor_b = *suffix_b;

        for _ in 0..word_count {
            let slice_a = &input[cursor_a..cursor_a + LANES];
            let simd_word_a = Word::new(slice_a.try_into().unwrap());
            cursor_a += LANES;

            let slice_b = &input[cursor_b..cursor_b + LANES];
            let simd_word_b = Word::new(slice_b.try_into().unwrap());
            cursor_b += LANES;

            let equal = simd_word_a.simd_eq(simd_word_b).to_bitmask();
            if equal != 0 {
                continue;
            }

            let less = simd_word_a.simd_lt(simd_word_b).to_bitmask();
            if less != 0 {
                return std::cmp::Ordering::Less;
            }
            return std::cmp::Ordering::Greater;
        }

        // Ensure the suffixes are stably sorted.
        suffix_a.cmp(suffix_b)
    });

    suffix_array
}

struct FullAuxiliaryData {
    kmer_count: usize,
    bwtk: Bwt,
    lcp: Lcp,
    shorter_than_k: BitVector,
    equal_to_k: RawVector,
    k_minus_one_ranges: BitVector,
    k_ranges: RawVector,
}

fn build_full_auxiliary_data(
    input: Vec<u8>,
    bounded_context_suffix_array: Vec<usize>,
    k: usize
) -> FullAuxiliaryData {
    // The input should still be padded with '#'.
    let length = bounded_context_suffix_array.len();
    let word_count = k.div_ceil(LANES);

    let mut bwt_raw_vectors = [
        RawVector::with_len(length, false), // $
        RawVector::with_len(length, false), // A
        RawVector::with_len(length, false), // C
        RawVector::with_len(length, false), // G
        RawVector::with_len(length, false), // T
    ];
    // The suffix which is equivalent to the entire string is always at the beginning since the
    // input string should begin with '#' which is lexicographically the smallest. The input
    // sequence must also end with '$' so the preceding character of '#' is always '$'
    bwt_raw_vectors[0].set_bit(0, true);

    let mut lcp = {
        let k_bit_width = usize::BITS - k.leading_zeros();
        let width = (k_bit_width.div_ceil(u8::BITS) as usize).next_power_of_two();
        let data = Vec::<u8>::with_capacity(length * width);
        Lcp::new_with_width(data, width)
    };
    lcp.push(0);

    let mut kmer_count: usize = 0;
    let mut shorter_than_k     = RawVector::with_len(length, false);
    let mut equal_to_k         = RawVector::with_len(length, false);
    let mut k_minus_one_ranges = RawVector::with_len(length, false);
    let mut k_ranges           = RawVector::with_len(length, false);

    for rank in 1..length {
        let current_suffix_index = bounded_context_suffix_array[rank];

        {
            // Set the bit in the bitvector corresponding to the previous character in the input
            // for the (bounded) BWT.
            let previous_character_index_in_input = (current_suffix_index + length - 1) % length;
            let previous_character = input[previous_character_index_in_input];
            let previous_character_index = CHAR_TO_INDEX[previous_character as usize];
            if previous_character_index < bwt_raw_vectors.len() {
                bwt_raw_vectors[previous_character_index].set_bit(rank, true);
            }
        }

        let previous_suffix_index = bounded_context_suffix_array[rank - 1];
        let (length, lcp_value) = find_length_and_lcp_value(
            &input,
            word_count,
            current_suffix_index,
            previous_suffix_index
        );
        lcp.push(lcp_value);

        if length < k {
            shorter_than_k.set_bit(rank, true);
        } else if length == k {
            equal_to_k.set_bit(rank, true);
        }

        if lcp_value < length {
            k_ranges.set_bit(rank, true);
            if length >= k {
                kmer_count += 1;
            }
            if length < k || lcp_value < k - 1 {
                k_minus_one_ranges.set_bit(rank, true);
            }
        }
    }

    let mut shorter_than_k = BitVector::from(shorter_than_k);
    shorter_than_k.enable_rank();

    let mut k_minus_one_ranges = BitVector::from(k_minus_one_ranges);
    k_minus_one_ranges.enable_rank();
    k_minus_one_ranges.enable_select();

    let mut bwt_bit_vectors = bwt_raw_vectors.into_iter().map(BitVector::from).collect::<Vec<_>>();
    bwt_bit_vectors.iter_mut().for_each(|vector| {
        vector.enable_rank();
        vector.enable_select();
    });
    let bwtk = Bwt::new(bwt_bit_vectors);

    FullAuxiliaryData {
        kmer_count,
        bwtk,
        lcp,
        shorter_than_k,
        equal_to_k,
        k_minus_one_ranges,
        k_ranges
    }
}

#[inline]
fn find_length_and_lcp_value(
    input: &[u8],
    word_count: usize,
    current_suffix_index: usize,
    previous_suffix_index: usize
) -> (usize, usize) {
    let mut length = 0;
    let mut stop_accumulating_length = false;

    let mut lcp_value = 0;
    let mut found_lcp = false;
    for _ in 0..word_count {
        let slice_for_current = &input[current_suffix_index..current_suffix_index + LANES];
        let simd_word_current = Word::new(slice_for_current.try_into().unwrap());
        let slice_for_previous = &input[previous_suffix_index..previous_suffix_index + LANES];
        let simd_word_previous = Word::new(slice_for_previous.try_into().unwrap());

        if !stop_accumulating_length {
            let bitmask = simd_word_current.simd_eq(b'$').to_bitmask();
            length += bitmask.trailing_zeros();
            if bitmask != 0 {
                stop_accumulating_length = true;
            }
        }

        if !found_lcp {
            let bitmask = simd_word_current.simd_eq(simd_word_previous).to_bitmask();
            lcp_value += bitmask.trailing_ones();
            if bitmask != u32::MAX {
                found_lcp = true;
            }
            if lcp_value >= length {
                found_lcp = true
            }
        }

        if found_lcp && stop_accumulating_length {
            break;
        }
    }

    (length as usize, lcp_value as usize)
}


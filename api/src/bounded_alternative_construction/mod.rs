// Code by Martin Kostadinov.

use crate::{LcsArray, SbwtIndex, SubsetSeq};

use super::alternative_construction::{
    Output,
    FULL_SET,
    push_set,
    input_structures::{
        Bwt, Lcp, CHAR_TO_INDEX
    }
};

use rayon::slice::ParallelSliceMut;
use bitvec::vec::BitVec;
use simple_sds_sbwt::{bit_vector::BitVector, int_vector::IntVector};
use simple_sds_sbwt::ops::{Access, BitVec as BitVecTrait, Push, Rank, Select};
use simple_sds_sbwt::raw_vector::{AccessRaw, RawVector};

type Word = wide::u8x32;
const LANES: usize = Word::LANES as usize;

pub fn build_without_redundant_dummies<SS: SubsetSeq + Send>(
    mut input: Vec<u8>, threads: usize, k: usize, build_lcs: bool
) -> Output<SS> {
    let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(threads).build().unwrap();
    let mut result: Option<Vec<usize>> = None;
    thread_pool.scope(|scope| {
        scope.spawn(|_| {
            result = Some(par_bounded_context_suffix_array(&mut input, k));
        });
    });
    let bounded_context_suffix_array = result.unwrap();
    println!("{:?}", bounded_context_suffix_array);

    let aux = build_full_auxiliary_data(input, bounded_context_suffix_array, k);
    let dummy_marks = build_dummy_marks(&aux, k);

    let FullAuxiliaryData {
        kmer_count,
        bwtk,
        lcp,
        shorter_than_k,
        equal_to_k,
        k_minus_one_ranges,
        k_ranges
    } = aux;
    drop(equal_to_k);

    let mut rows = Vec::<BitVec<u64>>::new();
    for _ in 0..4 {
        rows.push(BitVec::with_capacity(bwtk.len()));
    }

    let mut lcs: Option<IntVector> = if build_lcs {
        let bit_width = usize::BITS - (k.overflowing_sub(1).0).leading_zeros();
        let bit_width = bit_width as usize;
        let mut value = IntVector::with_capacity(bwtk.len(), bit_width).unwrap();
        value.push(0); // '$...$' dummy k-mer
        Some(value)
    } else {
        None
    };

    let mut current_set: u8 = 0;
    let separator_count = bwtk.counts[1];

    for index in 1..separator_count {
        if dummy_marks.keep_dummy.bit(index) {
            current_set |= (1 << dummy_marks.outedge.get(index)) as u8;
            if current_set & FULL_SET == FULL_SET {
                break;
            }
        }
    }
    push_set(&mut rows, current_set);

    current_set = 0;
    let mut current_lcs_value = k - 1;
    let mut include_dummy_kmer = false;
    let mut has_dummy_kmer     = false;
    let mut k_range_count = 0;

    for index in separator_count..bwtk.len() {
        if k_minus_one_ranges.get(index) {
            println!("---");
            if has_dummy_kmer && !include_dummy_kmer {
                k_range_count -= 1;
            }
            while k_range_count > 0 {
                println!("{:?}", current_set);
                push_set(&mut rows, current_set);
                current_set = 0;
                k_range_count -= 1;
            }

            current_set = 0;
            has_dummy_kmer = false;
            include_dummy_kmer = false;
            k_range_count = 0;
        }

        if build_lcs {
            current_lcs_value = current_lcs_value.min(lcp.get(index));
        }

        let is_start_of_k_range = k_ranges.bit(index);
        if is_start_of_k_range {
            k_range_count += 1;
        }

        if shorter_than_k.get(index) {
            has_dummy_kmer = true;
            let outedge = dummy_marks.outedge.get(index);
            if dummy_marks.keep_dummy.bit(index) {
                if !include_dummy_kmer && build_lcs {
                    lcs.as_mut().unwrap().push(current_lcs_value as u64);
                    current_lcs_value = k - 1;
                }
                include_dummy_kmer = true;
                current_set |= (1 << outedge) as u8;
            }
            if outedge != 0 {
                current_set |= (1 << outedge) as u8;
            }
        } else {
            if build_lcs && is_start_of_k_range {
                lcs.as_mut().unwrap().push(current_lcs_value as u64);
                current_lcs_value = k - 1;
            }
            current_set |= (1 << bwtk.get_char_index(index)) as u8;
        }
    }

    println!("---");
    if has_dummy_kmer && !include_dummy_kmer {
        k_range_count -= 1;
    }
    while k_range_count > 0 {
        println!("{:?}", current_set);
        push_set(&mut rows, current_set);
        current_set = 0;
        k_range_count -= 1;
    }

    let C: Vec<usize> = crate::util::get_C_array(&rows);
    let mut subset_rank = SS::new_from_bit_vectors(rows);
    subset_rank.build_rank();
    let n_sets  = subset_rank.len();
    let n_kmers = kmer_count;
    let mut index = SbwtIndex::<SS>::from_components(
        subset_rank,
        n_kmers,
        k,
        C,
        crate::PrefixLookupTable::new_empty(n_sets)
    );
    let prefix_lookup_table = crate::PrefixLookupTable::new(&index, 8);
    index.set_lookup_table(prefix_lookup_table);
    let lcs = lcs.map(LcsArray::new);

    Output {
        sbwt: index,
        lcs,
        counts: None,
    }
}

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
            if equal == u32::MAX {
                continue;
            }

            let less    = simd_word_a.simd_lt(simd_word_b).to_bitmask();
            let greater = simd_word_a.simd_gt(simd_word_b).to_bitmask();
            if less.trailing_zeros() < greater.trailing_zeros() {
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

const OUTEDGE_WIDTH: usize = 3;
struct DummyMarks {
    keep_dummy: RawVector,
    outedge: IntVector,
}

fn build_dummy_marks(aux: &FullAuxiliaryData, k: usize) -> DummyMarks {
    use super::alternative_construction::has_full_kmer_predecessor;

    let len = aux.bwtk.len();
    let mut keep_dummy = RawVector::with_len(len, false);
    let mut outedge    = IntVector::with_len(len, OUTEDGE_WIDTH, 0).unwrap();

    let start = aux.bwtk.counts[1];
    let mut predecessor_confirmed = false;
    for index in start..len {
        if aux.k_ranges.bit(index) {
            predecessor_confirmed = false;
        }

        if aux.equal_to_k.bit(index) {
            let (predecessor, char_index) = aux.bwtk.inverse_lf_step(index);
            if !predecessor_confirmed {
                predecessor_confirmed |= has_full_kmer_predecessor(
                    predecessor, &aux.bwtk, &aux.k_minus_one_ranges, &aux.shorter_than_k
                );
            }

            if predecessor_confirmed {
                outedge.set(predecessor, char_index as u64);
            } else {
                keep_predecessors(predecessor, char_index, &aux.bwtk, k, &mut keep_dummy, &mut outedge);
            }
        } else if !aux.shorter_than_k.get(index) {
            predecessor_confirmed = true;
        }
    }

    DummyMarks {
        keep_dummy,
        outedge,
    }
}

fn keep_predecessors(
    mut predecessor: usize,
    mut char_index: usize,
    bwtk: &Bwt,
    mut k: usize,
    keep_dummy: &mut RawVector,
    outedge: &mut IntVector,
) {
    while k > 0 {
        keep_dummy.set_bit(predecessor, true);
        outedge.set(predecessor, char_index as u64);
        let (next_predecessor, next_outedge) = bwtk.inverse_lf_step(predecessor);
        predecessor = next_predecessor;
        char_index = next_outedge;
        k -= 1;
    }
}


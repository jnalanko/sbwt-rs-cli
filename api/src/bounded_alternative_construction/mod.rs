// Code by Martin Kostadinov.

use crate::SubsetSeq;

use crate::vodbg::count::{Counts, Sample};
use crate::alternative_construction::{
    Output,
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

pub fn build<SS: SubsetSeq + Send>(
    threads: usize,
    input: Vec<u8>,
    k: usize,
    build_lcs: bool,
    add_all_dummies: bool,
    build_counts: bool
) -> Output<SS> {
    log::info!("[build] begin [length: {}]", input.len());
    let result = if !add_all_dummies {
        build_without_redundant_dummies(threads, input, k, build_lcs)
    } else {
        build_with_all_dummies(threads, input, k, build_lcs, build_counts)
    };
    log::info!("[build] done");
    result
}

pub fn build_without_redundant_dummies<SS: SubsetSeq + Send>(
    threads: usize, mut input: Vec<u8>, k: usize, build_lcs: bool
) -> Output<SS> {
    log::info!("[build_without_redundant_dummies] begin");
    let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(threads).build().unwrap();
    let mut result: Option<Vec<usize>> = None;
    thread_pool.scope(|scope| {
        scope.spawn(|_| {
            result = Some(par_bounded_context_suffix_array(&mut input, k));
        });
    });
    let bounded_context_suffix_array = result.unwrap();

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

    let separator_count = bwtk.counts[1];
    let (rows, lcs) = super::alternative_construction::_build_without_redundant_dummies(
        k,
        separator_count,
        bwtk.len(),
        &lcp,
        &k_minus_one_ranges,
        &k_ranges,
        &shorter_than_k,
        &dummy_marks.keep_dummy,
        build_lcs,
        |current_set, index| {
            let outedge = dummy_marks.outedge.get(index);
            *current_set |= (1 << outedge) as u8;
        },
        |current_set, index| {
            let outedge = bwtk.get_char_index(index);
            *current_set |= (1 << outedge) as u8;
        }
    );
    let result = super::alternative_construction::collect_output(k, rows, lcs, None, kmer_count);
    log::info!("[build_without_redundant_dummies] done");
    result
}

/// Needs to be executed in a rayon context.
pub fn par_bounded_context_suffix_array(input: &mut Vec<u8>, k: usize) -> Vec<usize> {
    log::info!("[par_bounded_context_suffix_array] begin");
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

    log::info!("[par_bounded_context_suffix_array] done");
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
    log::info!("[build_full_auxiliary_data] begin");
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
            k,
            &input,
            word_count,
            current_suffix_index,
            previous_suffix_index
        );
        let length = length.min(k);
        let lcp_value = lcp_value.min(length);
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

    log::info!("[build_full_auxiliary_data] shorter_than_k");
    let mut shorter_than_k = BitVector::from(shorter_than_k);
    shorter_than_k.enable_rank();

    log::info!("[build_full_auxiliary_data] k_minus_one_ranges");
    let mut k_minus_one_ranges = BitVector::from(k_minus_one_ranges);
    k_minus_one_ranges.enable_rank();
    k_minus_one_ranges.enable_select();

    log::info!("[build_full_auxiliary_data] bwt");
    let mut bwt_bit_vectors = bwt_raw_vectors.into_iter().map(BitVector::from).collect::<Vec<_>>();
    bwt_bit_vectors.iter_mut().for_each(|vector| {
        vector.enable_rank();
        vector.enable_select();
    });
    let bwtk = Bwt::new(bwt_bit_vectors);

    log::info!("[build_full_auxiliary_data] done");
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
    k: usize,
    input: &[u8],
    word_count: usize,
    current_suffix_index: usize,
    previous_suffix_index: usize
) -> (usize, usize) {
    let k = k as u32;

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
            if length >= k {
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
    log::info!("[build_dummy_marks] begin");
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

    log::info!("[build_dummy_marks] done");
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

pub fn build_with_all_dummies<SS: SubsetSeq + Send>(
    threads: usize, mut input: Vec<u8>, k: usize, build_lcs: bool, build_counts: bool
) -> Output<SS> {
    log::info!("[build_with_all_dummies] begin");
    let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(threads).build().unwrap();
    let mut result: Option<Vec<usize>> = None;
    thread_pool.scope(|scope| {
        scope.spawn(|_| {
            result = Some(par_bounded_context_suffix_array(&mut input, k));
        });
    });
    let bounded_context_suffix_array = result.unwrap();
    let word_count = k.div_ceil(LANES);
    let length = bounded_context_suffix_array.len();

    let mut rows = Vec::<BitVec<u64>>::new();
    for _ in 0..4 {
        rows.push(BitVec::with_capacity(length));
    }

    let mut lcs: Option<IntVector> = if build_lcs {
        let bit_width = usize::BITS - (k.overflowing_sub(1).0).leading_zeros();
        let bit_width = bit_width as usize;
        let mut value = IntVector::with_capacity(length, bit_width).unwrap();
        value.push(0); // '$...$' dummy k-mer
        Some(value)
    } else {
        None
    };

    let mut counts: Option<Counts> = if build_counts {
        let sample_capacity = length / Counts::DEFAULT_SAMPLE_DISTANCE;
        let mut value = Counts {
            individual_counts: Vec::with_capacity(length),
            sample_distance: Counts::DEFAULT_SAMPLE_DISTANCE,
            sample_information: Vec::with_capacity(sample_capacity + 2),
            large_counts: Vec::with_capacity(sample_capacity),
        };
        value.sample_information.push(Sample { count: 0, large_counts_up_to_sample: 0 });
        Some(value)
    } else {
        None
    };

    let mut current_set: u8 = 0;

    let mut kmer_count: usize = 0;
    let mut sets_count: usize = 1;

    let mut k_range_count   : u64 = 1; // Always output the $ set.
    let mut count           : u64 = 0;
    let mut individual_count: u64 = 0;
    let mut large_counts_up_to_sample: usize = 0;

    let mut start_of_k_range;
    let mut start_of_k_minus_one_range;

    for rank in 1..length {
        let current_suffix_index = bounded_context_suffix_array[rank];

        let previous_character_index_in_input = (current_suffix_index + length - 1) % length;
        let previous_character = input[previous_character_index_in_input];
        let outedge = CHAR_TO_INDEX[previous_character as usize];

        let first_character = input[current_suffix_index];
        let first_character_index = CHAR_TO_INDEX[first_character as usize];

        let (length, lcp_value) = if first_character_index == 0 {
            (0, 0)
        } else {
            let previous_suffix_index = bounded_context_suffix_array[rank - 1];
            let (length, lcp_value) = find_length_and_lcp_value(
                k,
                &input,
                word_count,
                current_suffix_index,
                previous_suffix_index
            );
            let length = length.min(k);
            let lcp_value = lcp_value.min(length);
            (length, lcp_value)
        };

        start_of_k_range           = lcp_value < length;
        start_of_k_minus_one_range = start_of_k_range && (length < k || lcp_value < k - 1);

        if start_of_k_minus_one_range {
            while k_range_count > 0 {
                push_set(&mut rows, current_set);
                current_set = 0;
                k_range_count -= 1;
            }

            current_set = 0;
        }

        if start_of_k_range {
            k_range_count += 1;
            if length >= k {
                kmer_count += 1;
            }

            if build_lcs {
                lcs.as_mut().unwrap().push(lcp_value as u64);
            }

            if build_counts {
                let counts = counts.as_mut().unwrap();
                if sets_count % Counts::DEFAULT_SAMPLE_DISTANCE == 0 {
                    counts.sample_information.push(Sample {
                        count,
                        large_counts_up_to_sample,
                    });
                }
                if individual_count >= u8::MAX as u64 {
                    counts.individual_counts.push(u8::MAX);
                    counts.large_counts.push(individual_count - u8::MAX as u64);
                } else {
                    counts.individual_counts.push(individual_count as u8);
                }
                individual_count = 0;
            }

            sets_count += 1;
        }

        if build_counts && length > 0 {
            count += 1;
            individual_count += 1;
            if  individual_count == u8::MAX as u64 {
                large_counts_up_to_sample += 1;
            }
        }

        current_set |= (1 << outedge) as u8;
    }

    while k_range_count > 0 {
        push_set(&mut rows, current_set);
        current_set = 0;
        k_range_count -= 1;
    }

    if build_counts {
        let counts = counts.as_mut().unwrap();
        counts.sample_information.push(Sample {
            count,
            large_counts_up_to_sample,
        });
        if individual_count >= u8::MAX as u64 {
            counts.individual_counts.push(u8::MAX);
            counts.large_counts.push(individual_count - u8::MAX as u64);
        } else {
            counts.individual_counts.push(individual_count as u8);
        }
    }

    let result = super::alternative_construction::collect_output(k, rows, lcs, counts, kmer_count);
    log::info!("[build_with_all_dummies] done");
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alternative_construction::make_concatenation;
    use crate::{BitPackedKmerSortingMem, SbwtIndexBuilder, SubsetMatrix, VecSeqStream};

    #[test]
    fn randomised_kmers() {
        use rand_chacha::ChaCha20Rng;
        use rand_chacha::rand_core::SeedableRng;
        use rand_chacha::rand_core::RngCore;

        let k: usize = 16;
        let kmer_length: usize = 48;
        let kmer_count = 256;
        let mut rng = ChaCha20Rng::from_seed([42; 32]);

        let mut seqs = Vec::<Vec<u8>>::new();
        for _ in 0..kmer_count {
            let kmer: Vec<u8> = (0..kmer_length).map(|_| match rng.next_u32() % 4 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            }).collect();
            seqs.push(kmer);
        }

        seqs.push(b"WITH_INCORRECT_CHARACTERS".to_vec());
        seqs.push(vec![b'A'; 1024]);

        seqs.sort();
        seqs.dedup();

        let concatenation = make_concatenation(&seqs);

        {
            // Without redundant dummies.
            let (correct_sbwt, correct_lcs) = SbwtIndexBuilder::<BitPackedKmerSortingMem>::new()
                .k(k).build_lcs(true)
                .run_from_vecs(&seqs);

            let Output {
                sbwt: constructed_sbwt,
                lcs: constructed_lcs,
                counts
            } = build::<SubsetMatrix>(2, concatenation.clone(), k, true, false, false);

            assert!(counts.is_none());

            let mut correct_buf = Vec::<u8>::new();
            let mut constructed_buf = Vec::<u8>::new();

            correct_sbwt.serialize(&mut correct_buf).unwrap();
            constructed_sbwt.serialize(&mut constructed_buf).unwrap();
            assert_eq!(correct_buf, constructed_buf);

            correct_buf.clear();
            constructed_buf.clear();

            let correct_lcs = correct_lcs.unwrap();
            let constructed_lcs = constructed_lcs.unwrap();

            correct_lcs    .serialize(&mut correct_buf).unwrap();
            constructed_lcs.serialize(&mut constructed_buf).unwrap();
            assert_eq!(correct_buf, constructed_buf);
        }

        {
            // With all dummies.
            let (mut correct_sbwt, correct_lcs) = SbwtIndexBuilder::<BitPackedKmerSortingMem>::new()
                .k(k).build_lcs(true)
                .add_all_dummy_paths(true)
                .run_from_vecs(&seqs);

            let Output {
                sbwt: constructed_sbwt,
                lcs: constructed_lcs,
                counts
            } = build::<SubsetMatrix>(2, concatenation, k, true, true, true);

            let mut correct_buf = Vec::<u8>::new();
            let mut constructed_buf = Vec::<u8>::new();

            correct_sbwt.serialize(&mut correct_buf).unwrap();
            constructed_sbwt.serialize(&mut constructed_buf).unwrap();
            assert_eq!(correct_buf, constructed_buf);

            correct_buf.clear();
            constructed_buf.clear();

            let correct_lcs = correct_lcs.unwrap();
            let constructed_lcs = constructed_lcs.unwrap();
            assert_eq!(correct_lcs.len(), constructed_lcs.len());

            correct_lcs    .serialize(&mut correct_buf).unwrap();
            constructed_lcs.serialize(&mut constructed_buf).unwrap();
            assert_eq!(correct_buf, constructed_buf);

            correct_sbwt.build_select();
            let pnsv = crate::vodbg::pnsv::PnsvTuned::new_default(&correct_sbwt, &correct_lcs, k);
            let mut vodbg = crate::vodbg::VoDbg::new(&correct_sbwt, &pnsv);
            vodbg.build_counts(
                VecSeqStream::new(&seqs),
                true,
                Counts::DEFAULT_SAMPLE_DISTANCE,
                1, 4, Counts::DEFAULT_BATCH_SIZE_IN_BYTES).unwrap();

            correct_buf.clear();
            constructed_buf.clear();

            vodbg.counts().unwrap().serialize(&mut correct_buf).unwrap();
            counts.unwrap().serialize(&mut constructed_buf).unwrap();
            assert_eq!(correct_buf, constructed_buf);
        }
    }
}

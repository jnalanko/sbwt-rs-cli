// Code by Martin Kostadinov.

pub mod preprocessing;
pub mod input_structures;

use input_structures::{Bwt, Lcp, CHAR_TO_INDEX};
use crate::atomic_bitmap::AtomicBitmap;
use crate::vodbg::count::{Counts, Sample};
use crate::{LcsArray, SbwtIndex, SubsetSeq};

use std::sync::{Arc, Mutex, atomic::AtomicUsize};

use bitvec::vec::BitVec;
use simple_sds_sbwt::{bit_vector::BitVector, int_vector::IntVector};
use simple_sds_sbwt::ops::{Access, BitVec as BitVecTrait, Push, Rank, Select};
use simple_sds_sbwt::raw_vector::{AccessRaw, RawVector};

type Word = wide::u8x16;
const LANES: usize = Word::LANES as usize;
type Bitmask = u16;
const _: () = assert!(Bitmask::BITS as usize == LANES);
const FULL_SET: u8 = 0b00011110;

pub struct Output<SS: SubsetSeq + Send> {
    pub sbwt: SbwtIndex<SS>,
    pub lcs: Option<LcsArray>,
    pub counts: Option<Counts>,
}

pub fn build<SS: SubsetSeq + Send>(
    threads: usize,
    input: Vec<u8>,
    suffix_array: Vec<usize>,
    k: usize,
    build_lcs: bool,
    add_all_dummies: bool,
    build_counts: bool
) -> Output<SS> {
    log::info!("[build] begin");
    let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(threads).build().unwrap();
    let mut result = None;
    thread_pool.scope(|scope| {
        scope.spawn(|_| {
            let output = par_build(threads, input, suffix_array, k, build_lcs, add_all_dummies, build_counts);
            result = Some(output);
        });
    });
    log::info!("[build] done");
    result.unwrap()
}

pub fn build_with_bounded_suffix_array<SS: SubsetSeq + Send>(
    threads: usize,
    mut input: Vec<u8>,
    k: usize,
    prefix_length_for_bucket_sort: usize,
    build_lcs: bool,
    add_all_dummies: bool,
    build_counts: bool
) -> Output<SS> {
    log::info!("[build_with_bounded_suffix_array] begin");
    let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(threads).build().unwrap();
    let mut result = None;
    thread_pool.scope(|scope| {
        scope.spawn(|_| {
            let suffix_array = par_bounded_context_suffix_array_bucket_sort(&mut input, k, prefix_length_for_bucket_sort);
            let output = par_build::<SS>(threads, input, suffix_array, k, build_lcs, add_all_dummies, build_counts);
            result = Some(output);
        });
    });
    log::info!("[build_with_bounded_suffix_array] done");
    result.unwrap()
}

pub fn par_build<SS: SubsetSeq + Send>(
    threads: usize,
    input: Vec<u8>,
    suffix_array: Vec<usize>,
    k: usize,
    build_lcs: bool,
    add_all_dummies: bool,
    build_counts: bool
) -> Output<SS> {
    log::info!(
        "[par_build] begin [length: {} | build_lcs: {} | add_all_dummies: {} | build_counts {}]",
        input.len(), build_lcs, add_all_dummies, build_counts
    );
    let output = if !add_all_dummies {
        par_build_without_redundant_dummies(threads, input, suffix_array, k, build_lcs)
    } else {
        build_with_all_dummies(input, suffix_array, k, build_lcs, build_counts)
    };
    log::info!("[par_build] done");
    output
}

pub fn par_build_without_redundant_dummies<SS: SubsetSeq + Send>(
    threads: usize, input: Vec<u8>, bounded_context_suffix_array: Vec<usize>, k: usize, build_lcs: bool
) -> Output<SS> {
    log::info!("[par_build_without_redundant_dummies] begin");
    let _ = threads;
    let aux = par_build_full_auxiliary_data(threads, input, bounded_context_suffix_array, k);
    let dummy_marks = par_build_dummy_marks(threads, k, &aux);
    let kmer_count = aux.kmer_count;

    let (rows, lcs) = _par_build_without_redundant_dummies(
        threads,
        k,
        aux,
        dummy_marks,
        build_lcs,
    );
    let result = collect_output(k, rows, lcs, None, kmer_count);
    log::info!("[par_build_without_redundant_dummies] done");
    result
}

#[allow(clippy::too_many_arguments)]
#[inline]
fn _par_build_without_redundant_dummies(
    threads: usize,
    k: usize,
    aux: FullAuxiliaryData,
    dummy_marks: DummyMarks,
    build_lcs: bool,
) -> (Vec<BitVec<u64>>, Option<IntVector>)
{
    log::info!("[_par_build_without_redundant_dummies] begin");

    let length = aux.bwt.len();
    let separator_count = aux.bwt.counts[1];

    let mut dummy_set: u8 = 0;
    for index in 1..separator_count {
        if dummy_marks.keep_dummy.bit(index) {
            dummy_set = include_letter(&aux.bwt, index, dummy_set);
        }
        if dummy_set & FULL_SET == FULL_SET {
            break;
        }
    }

    let results = Arc::new(Mutex::new(Vec::<SbwtRegionResult>::with_capacity(threads)));
    let threads_total_length = length - separator_count;
    let thread_region_length = threads_total_length.div_ceil(threads);

    log::info!("[_par_build_without_redundant_dummies] regions begin");
    rayon::scope(|s| {

        let aux         = &aux;
        let dummy_marks = &dummy_marks;
        let results     = &results;

        for thread_index in 0..threads {
            s.spawn(move |_| {
                let start = separator_count + thread_index * thread_region_length;
                let end = (start + thread_region_length).min(length);
                let is_last_thread = thread_index == threads - 1;
                let result = _build_sbwt_region(
                    k,
                    start,
                    end,
                    is_last_thread,
                    aux,
                    dummy_marks,
                    build_lcs
                );
                results.lock().unwrap().push(result);
            });
        }
    });
    log::info!("[_par_build_without_redundant_dummies] regions done");

    drop(aux);

    let mut results = {
        // No other thread should be holding the mutex now.
        let mutex = Arc::try_unwrap(results).unwrap();
        mutex.into_inner().unwrap()
    };
    results.sort_by_key(|result| result.start);

    log::info!("[_par_build_without_redundant_dummies] agglomerating result");
    let rows = Arc::new(Mutex::new(Vec::<(usize, BitVec<u64>)>::new()));
    let mut lcs = None;
    rayon::scope(|s| {
        let results = &results;
        let rows = &rows;
        let lcs = &mut lcs;

        if build_lcs {
            s.spawn(move |_| {
                let bit_width = usize::BITS - (k.overflowing_sub(1).0).leading_zeros();
                let bit_width = bit_width as usize;
                let mut value = IntVector::with_capacity(length, bit_width).unwrap();
                value.push(0);

                let mut last_lcs_value: usize = 0;
                for result in results {
                    if result.rows[0].is_empty() {
                        last_lcs_value = last_lcs_value.min(result.last_lcs_value);
                        continue;
                    }
                    let result_lcs = result.lcs.as_ref().unwrap();
                    let first_value = result_lcs.get(0) as usize;
                    let new_first_value =  first_value.min(last_lcs_value);
                    value.push(new_first_value as u64);
                    value.extend(result_lcs.iter().skip(1));
                    last_lcs_value = result.last_lcs_value;
                }

                *lcs = Some(value);
            });
        }

        for row_index in 0..4 {
            s.spawn(move |_| {
                let mut local_row = BitVec::with_capacity(length);
                local_row.push(false);
                for result in results {
                    local_row.extend_from_bitslice(result.rows[row_index].as_bitslice());
                }
                rows.lock().unwrap().push((row_index, local_row));
            });
        }
    });

    let mut rows_with_index = {
        let mutex = Arc::try_unwrap(rows).unwrap();
        mutex.into_inner().unwrap()
    };
    rows_with_index.sort_by_key(|(index, _)| *index);
    let mut rows = rows_with_index.into_iter().map(|(_, row)| row).collect::<Vec<_>>();
    for i in 1..=4 {
        rows[i - 1].set(0, dummy_set & (1 << i) != 0);
    }

    log::info!("[_par_build_without_redundant_dummies] done");
    (rows, lcs)
}

#[derive(Debug)]
struct SbwtRegionResult {
    start: usize,
    rows: Vec<BitVec<u64>>,
    lcs: Option<IntVector>,
    last_lcs_value: usize,
}

#[allow(clippy::too_many_arguments)]
fn _build_sbwt_region(
    k: usize,
    start: usize,
    end: usize,
    is_last_thread: bool,
    aux: &FullAuxiliaryData,
    dummy_marks: &DummyMarks,
    build_lcs: bool,
) -> SbwtRegionResult {
    let length = aux.bwt.len();
    let capacity = end - start;

    let mut rows = Vec::<BitVec<u64>>::new();
    for _ in 0..4 {
        rows.push(BitVec::with_capacity(capacity));
    }

    let mut lcs: Option<IntVector> = if build_lcs {
        let bit_width = usize::BITS - (k.overflowing_sub(1).0).leading_zeros();
        let bit_width = bit_width as usize;
        let value = IntVector::with_capacity(capacity, bit_width).unwrap();
        Some(value)
    } else {
        None
    };

    let mut index = start;
    let mut current_set = 0;
    let mut current_lcs_value  = k - 1;
    let mut include_dummy_kmer = false;
    let mut has_dummy_kmer     = false;
    let mut k_range_count = 0;

    while index < end {
        if build_lcs {
            current_lcs_value = current_lcs_value.min(aux.lcp.get(index));
        }
        if aux.k_minus_one_ranges.get(index) {
            break;
        }
        index += 1;
    }

    if index == end {
        return SbwtRegionResult {
            start,
            rows,
            lcs,
            last_lcs_value: current_lcs_value
        };
    }

    while index < length {
        if aux.k_minus_one_ranges.get(index) {
            if has_dummy_kmer && !include_dummy_kmer {
                k_range_count -= 1;
            }
            while k_range_count > 0 {
                push_set(&mut rows, current_set);
                current_set = 0;
                k_range_count -= 1;
            }

            current_set = 0;
            has_dummy_kmer = false;
            include_dummy_kmer = false;
            k_range_count = 0;

            if index >= end {
                break;
            }
        }

        if build_lcs {
            current_lcs_value = current_lcs_value.min(aux.lcp.get(index));
        }

        let is_start_of_k_range = aux.k_ranges.bit(index);
        if is_start_of_k_range {
            k_range_count += 1;
        }

        if aux.shorter_than_k.get(index) {
            has_dummy_kmer = true;
            if dummy_marks.keep_dummy.bit(index) {
                if build_lcs && !include_dummy_kmer {
                    lcs.as_mut().unwrap().push(current_lcs_value as u64);
                    current_lcs_value = k - 1;
                }
                include_dummy_kmer = true;
            }
            if dummy_marks.keep_dummy.bit(index) || aux.equal_to_k_minus_one_or_k.bit(index) {
                current_set = include_letter(&aux.bwt, index, current_set);
            }
        } else {
            if build_lcs && is_start_of_k_range {
                lcs.as_mut().unwrap().push(current_lcs_value as u64);
                current_lcs_value = k - 1;
            }
            current_set = include_letter(&aux.bwt, index, current_set);
        }

        index += 1;
    }

    if is_last_thread {
        if has_dummy_kmer && !include_dummy_kmer {
            k_range_count -= 1;
        }
        while k_range_count > 0 {
            push_set(&mut rows, current_set);
            current_set = 0;
            k_range_count -= 1;
        }
    }

    SbwtRegionResult {
        start,
        rows,
        lcs,
        last_lcs_value: current_lcs_value
    }
}

pub(crate) fn collect_output<SS: SubsetSeq + Send>(
    k: usize,
    rows: Vec<BitVec<u64>>,
    lcs: Option<IntVector>,
    counts: Option<Counts>,
    kmer_count: usize
) -> Output<SS> {
    log::info!("[collect_output] constructing sbwt");
    let C: Vec<usize> = crate::util::get_C_array(&rows);
    let mut subset_rank = SS::new_from_bit_vectors(rows);
    subset_rank.build_rank();
    let n_sets  = subset_rank.len();
    let n_kmers = kmer_count;
    let mut index = SbwtIndex::<SS>::from_parts(
        subset_rank,
        n_kmers,
        k,
        C,
        crate::PrefixLookupTable::new_empty(n_sets)
    );
    // TODO: add a parameter for prefix_length?
    let prefix_lookup_table = crate::PrefixLookupTable::new(&index, 8);
    index.set_lookup_table(prefix_lookup_table);
    let lcs = lcs.map(LcsArray::new);

    log::info!("[collect_output] constructing done");
    Output {
        sbwt: index,
        lcs,
        counts,
    }
}

const CHARACTER_COUNT: usize = 5;
const MAX_PREFIX_LENGTH_FOR_BUCKET_SORT: usize = 8;

/// The length is the length of the original input sequences. If the input buffer has length less
/// than that, ensure that it is padded.
fn pad_input(input: &mut Vec<u8>, k: usize, length: usize) {
    let word_count = k.div_ceil(LANES);
    if input.len() > length {
        return;
    }
    for _ in 1..(word_count * LANES) {
        input.push(b'$');
    }
}

pub fn par_bounded_context_suffix_array_bucket_sort(input: &mut Vec<u8>, k: usize, prefix_length_for_bucket_sort: usize) -> Vec<usize> {
    log::info!("[par_bounded_context_suffix_array_bucket_sort] begin");
    let length = input.len();
    let word_count = k.div_ceil(LANES);
    let prefix_length_for_bucket_sort = prefix_length_for_bucket_sort.min(MAX_PREFIX_LENGTH_FOR_BUCKET_SORT);
    let bucket_count: usize = CHARACTER_COUNT.pow(prefix_length_for_bucket_sort as u32);

    // Pad the input so that the comparisons need not do bound checks.
    pad_input(input, k, input.len());

    use rayon::prelude::*;
    rayon::scope(|s| {
        s.spawn(|_| {
            input.par_iter_mut().for_each(|character| {
                *character = CHAR_TO_INDEX[*character as usize] as u8;
            });
        });
    });

    let identify_bucket = |index: usize| -> usize {
        if input[index] == 0 {
            return 0;
        }
        let mut power = 1;
        let mut bucket = 0;
        for i in (0..prefix_length_for_bucket_sort).rev() {
            let digit = input[index + i] as usize;
            bucket += digit * power;
            power *= CHARACTER_COUNT;
        }
        bucket
    };

    // Count for the bucket sort step.
    use std::sync::atomic::{AtomicUsize, Ordering};
    let buckets = vec![0_usize; bucket_count].into_iter().map(AtomicUsize::new).collect::<Vec<_>>();
    rayon::scope(|s| {
        s.spawn(|_| {
            (1..length).into_par_iter().for_each(|index| {
                let bucket = identify_bucket(index);
                buckets[bucket].fetch_add(1, Ordering::AcqRel);
            });
        });
    });
    let buckets = buckets.into_iter().map(|value| value.load(Ordering::Relaxed)).collect::<Vec<_>>();

    let mut result = vec![0_usize; length];

    {
        // Sort the suffixes into buckets.
        let mut buckets_scratch = Vec::<usize>::with_capacity(buckets.len());
        let mut last = 1;
        buckets_scratch.push(last);
        for i in 1..buckets.len() {
            last += buckets[i - 1];
            buckets_scratch.push(last);
        }


        for item in 1..length {
            let bucket = identify_bucket(item);
            let index_in_swap = buckets_scratch[bucket];
            buckets_scratch[bucket] += 1;
            result[index_in_swap] = item;
        }
    }

    rayon::scope(|s| {
        use input_structures::INDEX_TO_CHAR;
        s.spawn(|_| {
            input.par_iter_mut().for_each(|character| *character = INDEX_TO_CHAR[*character as usize]);
        });
    });

    rayon::scope(|s| {
        // Sort the buckets.
        let start = 1 + buckets[0];
        let mut slice = &mut result[start..];
        #[allow(clippy::needless_range_loop)]
        for i in 1..buckets.len() {
            let (to_sort, residual) = slice.split_at_mut(buckets[i]);
            s.spawn(|_| {
                to_sort.sort_unstable_by(|a, b| {
                    compare_suffixes(input, word_count, *a, *b)
                });
            });
            slice = residual;
        }
    });

    log::info!("[par_bounded_context_suffix_array_bucket_sort] done");
    result
}

fn compare_suffixes(input: &[u8], word_count: usize, mut cursor_a: usize, mut cursor_b: usize) -> std::cmp::Ordering {
    for _ in 0..word_count {
        let slice_a = &input[cursor_a..cursor_a + LANES];
        let simd_word_a = Word::new(slice_a.try_into().unwrap());
        cursor_a += LANES;

        let slice_b = &input[cursor_b..cursor_b + LANES];
        let simd_word_b = Word::new(slice_b.try_into().unwrap());
        cursor_b += LANES;

        let equal = simd_word_a.simd_eq(simd_word_b).to_bitmask() as Bitmask;
        if equal == Bitmask::MAX {
            continue;
        }

        let less    = simd_word_a.simd_lt(simd_word_b).to_bitmask() as Bitmask;
        let greater = !(equal | less);
        if less.trailing_zeros() < greater.trailing_zeros() {
            return std::cmp::Ordering::Less;
        }
        return std::cmp::Ordering::Greater;
    }

    std::cmp::Ordering::Equal
}

pub(crate) struct FullAuxiliaryData {
    /// It is easier to calculate the true k-mers while creating the auxiliary bitvectors rather
    /// than trying to figure out their count later.
    pub(crate) kmer_count: usize,
    pub(crate) bwt: Bwt,
    pub(crate) lcp: Lcp,
    /// Used to figure out whether a k-mer at the beginning of a sequence has a true k-mer as a
    /// predecessor in order to figure out whether the dummy k-mer is necessary. In the final pass
    /// it is used to figure out if a given (k-1)-range contains a region of dummy k-mers.
    pub(crate) shorter_than_k: BitVector,
    /// Used in the pass which marks the dummy k-mers that need to be kept in order to identify the
    /// k-mers which are at the beginning of an input sequence.
    ///
    /// Also marks (k-1)-mers in order to keep their outedge as it is not guaranteed that the
    /// non-dummy k-mer that exists as a predecessor to a non-dummy k-mer at the beginning of a
    /// sequence has the needed label. This is needed as we don't want to include all letters of
    /// dummy k-mers whose "true" length is less than k-1 and there is no way to figure out which
    /// ones are those from the [FullAuxiliaryBitVectors::shorter_than_k] bitvector alone.
    pub(crate) equal_to_k_minus_one_or_k: RawVector,
    /// Used in order to figure out the bounds at which to check the [Self::shorter_than_k]
    /// bitvector whether a non-dummy k-mer at the beginning of an input sequence has a non-dummy
    /// k-mer as a predecessor. In addition, it is used to figure out the bounds at which to
    /// "collect" outedges for the first k-mer with a given (k-1) suffix in the SBWT.
    pub(crate) k_minus_one_ranges: BitVector,
    /// Used to enumerate the sets of the SBWT.
    pub(crate) k_ranges: RawVector,
}

pub(crate) fn par_build_full_auxiliary_data(
    threads: usize,
    mut input: Vec<u8>,
    bounded_context_suffix_array: Vec<usize>,
    k: usize
) -> FullAuxiliaryData {
    // The prefix of a suffix up to the first '$' will be referred to as the true prefix and
    // its length as the true length.
    //
    // An N-range is a contiguous range of suffixes which have the same true prefix after being
    // truncated from the right to a length of N (or if the true length of the suffix is less than
    // N, they are padded with imaginary '$').
    //
    // A (k-1)-range which contains suffixes with true lengths less than k-1 will be referred
    // to as a small range. A (k-1)-range which contains suffixes with true length equal to k
    // will be referred to as a big range.
    //
    // A big range can be further divided into k-ranges.

    log::info!("[par_build_full_auxiliary_data] begin");

    let length = bounded_context_suffix_array.len();
    pad_input(&mut input, k, length);
    let word_count = k.div_ceil(LANES);

    let bwt_vectors = [
        AtomicBitmap::new(length), // $
        AtomicBitmap::new(length), // A
        AtomicBitmap::new(length), // C
        AtomicBitmap::new(length), // G
        AtomicBitmap::new(length), // T
    ];
    bwt_vectors[0].set(0, true);

    let mut lcp_result_data = Vec::with_capacity(threads);
    for _ in 0..threads {
        lcp_result_data.push(None);
    }
    let k_bit_width = usize::BITS - k.leading_zeros();
    let lcp_width = (k_bit_width.div_ceil(u8::BITS) as usize).next_power_of_two();
    let lcp_result: Arc<Mutex<Vec<Option<Lcp>>>> = Arc::new(Mutex::new(
        lcp_result_data
    ));

    let kmer_count = AtomicUsize::new(0);
    let shorter_than_k     = AtomicBitmap::new(length);
    let equal_to_k         = AtomicBitmap::new(length);
    let k_minus_one_ranges = AtomicBitmap::new(length);
    let k_ranges           = AtomicBitmap::new(length);

    let region_size = length.div_ceil(threads);

    rayon::scope(|s| {
        let input              = &input;
        let bounded_context_suffix_array = &bounded_context_suffix_array;
        let bwt_vectors        = &bwt_vectors;
        let lcp_result         = &lcp_result;
        let kmer_count         = &kmer_count;
        let shorter_than_k     = &shorter_than_k;
        let equal_to_k         = &equal_to_k;
        let k_minus_one_ranges = &k_minus_one_ranges;
        let k_ranges           = &k_ranges;

        for thread_index in 0..threads {
            s.spawn(move |_| {

                let start = 1 + thread_index * region_size;
                let end = (start + region_size).min(length);
                let local_region_length = end - start;

                let mut local_lcp = {
                    let mut capacity = local_region_length * lcp_width;
                    if thread_index == 0 {
                        capacity += 1;
                    }
                    let data = Vec::<u8>::with_capacity(capacity);
                    Lcp::new_with_width(data, lcp_width)
                };

                if thread_index == 0 {
                    local_lcp.push(0);
                }

                for rank in start..end {
                    let current_suffix_index = bounded_context_suffix_array[rank];

                    // Set the bit in the bitvector corresponding to the previous character in the input
                    // for the (bounded) BWT.
                    let previous_character_index_in_input = (current_suffix_index + length - 1) % length;
                    let previous_character = input[previous_character_index_in_input];
                    let previous_character_index = CHAR_TO_INDEX[previous_character as usize];
                    if previous_character_index < bwt_vectors.len() {
                        bwt_vectors[previous_character_index].set(rank, true);
                    }

                    let previous_suffix_index = bounded_context_suffix_array[rank - 1];
                    let (length, lcp_value) = find_length_and_lcp_value(
                        k,
                        input,
                        word_count,
                        current_suffix_index,
                        previous_suffix_index
                    );
                    local_lcp.push(lcp_value);

                    if length < k {
                        shorter_than_k.set(rank, true);
                        if length == k - 1 {
                            equal_to_k.set(rank, true);
                        }
                    } else if length == k {
                        equal_to_k.set(rank, true);
                    }

                    if lcp_value < length {
                        k_ranges.set(rank, true);
                        if length >= k {
                            kmer_count.fetch_add(1, std::sync::atomic::Ordering::AcqRel);
                        }
                        if length < k || lcp_value < k - 1 {
                            k_minus_one_ranges.set(rank, true);
                        }
                    }
                }

                let mut lcp_ref = lcp_result.lock().unwrap();
                lcp_ref[thread_index] = Some(local_lcp);
            });
        }
    });

    // note(mk): Perhaps a bit more time can be reduced here if the conversion of the different
    // structures is done on separate threads.

    log::info!("[par_build_full_auxiliary_data] bwt");
    let bwt_bit_vectors = bwt_vectors.into_iter()
        .map(|value| {
            let intermediate = convert_atomic_bitmap(value);
            BitVector::from(intermediate)
        })
        .collect::<Vec<_>>();
    let bwtk = Bwt::new(bwt_bit_vectors);

    log::info!("[par_build_full_auxiliary_data] lcp");
    let mut collected_lcp_data = Vec::<u8>::with_capacity(length * lcp_width);
    for lcp in lcp_result.lock().unwrap().iter_mut() {
        let bytes: Vec<u8> = lcp.take().unwrap().into();
        collected_lcp_data.extend_from_slice(&bytes);
    }
    let lcp = Lcp::new_with_width(collected_lcp_data, lcp_width);

    fn convert_atomic_bitmap(value: AtomicBitmap) -> RawVector {
        let intermediate = value.into_bitvec_u64();
        crate::util::bitvec_to_simple_sds_raw_bitvec(intermediate)
    }

    log::info!("[par_build_full_auxiliary_data] shorter_than_k");
    let shorter_than_k     = convert_atomic_bitmap(shorter_than_k);
    let mut shorter_than_k = BitVector::from(shorter_than_k);
    shorter_than_k.enable_rank();

    let equal_to_k         = convert_atomic_bitmap(equal_to_k);

    log::info!("[par_build_full_auxiliary_data] k_minus_one_ranges");
    let k_minus_one_ranges = convert_atomic_bitmap(k_minus_one_ranges);
    let mut k_minus_one_ranges = BitVector::from(k_minus_one_ranges);
    k_minus_one_ranges.enable_rank();
    k_minus_one_ranges.enable_select();

    let k_ranges           = convert_atomic_bitmap(k_ranges);

    FullAuxiliaryData {
        kmer_count: kmer_count.load(std::sync::atomic::Ordering::Relaxed),
        bwt: bwtk,
        lcp,
        shorter_than_k,
        equal_to_k_minus_one_or_k: equal_to_k,
        k_minus_one_ranges,
        k_ranges
    }
}

#[inline]
fn find_length_and_lcp_value(
    k: usize,
    input: &[u8],
    word_count: usize,
    mut current_suffix_index: usize,
    mut previous_suffix_index: usize
) -> (usize, usize) {
    let k = k as u32;

    let mut length = 0;
    let mut stop_accumulating_length = false;

    let mut lcp_value = 0;
    let mut found_lcp = false;
    for _ in 0..word_count {
        let slice_for_current = &input[current_suffix_index..current_suffix_index + LANES];
        current_suffix_index += LANES;
        let simd_word_current = Word::new(slice_for_current.try_into().unwrap());

        let slice_for_previous = &input[previous_suffix_index..previous_suffix_index + LANES];
        previous_suffix_index += LANES;
        let simd_word_previous = Word::new(slice_for_previous.try_into().unwrap());

        if !stop_accumulating_length {
            let bitmask = simd_word_current.simd_eq(b'$').to_bitmask() as Bitmask;
            length += bitmask.trailing_zeros();
            if length >= k {
                stop_accumulating_length = true;
                length = k;
            }
            if bitmask != 0 {
                stop_accumulating_length = true;
            }
        }

        if !found_lcp {
            let bitmask = simd_word_current.simd_eq(simd_word_previous).to_bitmask() as Bitmask;
            lcp_value += bitmask.trailing_ones();
            if stop_accumulating_length && lcp_value >= length {
                found_lcp = true;
                lcp_value = length;
            }
            if bitmask != Bitmask::MAX {
                found_lcp = true;
            }
        }

        if found_lcp && stop_accumulating_length {
            break;
        }
    }

    (length as usize, lcp_value as usize)
}

struct DummyMarks {
    /// We keep a dummy (k-1)-mer and all of its predecessors if a k-mer at the beginning of an
    /// input sequence would not have a non-dummy k-mer as a predecessor in the SBWT graph.
    keep_dummy: RawVector,
}

fn par_build_dummy_marks(threads: usize, k: usize, aux: &FullAuxiliaryData) -> DummyMarks {
    log::info!("[par_build_dummy_marks] begin");

    let bwt = &aux.bwt;
    let scan_start = bwt.counts[1];
    let length = bwt.len();
    let keep_dummy  = AtomicBitmap::new(length);

    let scan_length = length - scan_start;
    let thread_region_length = scan_length.div_ceil(threads);

    rayon::scope(|s| {
        let keep_dummy  = &keep_dummy;
        for thread_index in 0..threads {
            s.spawn(move |_| {
                let start = scan_start + thread_index * thread_region_length;
                let end = (start + thread_region_length).min(length);

                let mut index = start;
                while index < end && !aux.k_ranges.bit(index) {
                    index += 1;
                }

                let mut predecessor_confirmed = false;
                while index < length {
                    if aux.k_ranges.bit(index) {
                        if index >= end {
                            break;
                        }
                        predecessor_confirmed = false;
                    }

                    if !aux.shorter_than_k.get(index) {
                        if aux.equal_to_k_minus_one_or_k.bit(index) {
                            // Equal to k.
                            let predecessor = bwt.inverse_lf_step(index);
                            if !predecessor_confirmed {
                                predecessor_confirmed |= has_full_kmer_predecessor(
                                    predecessor, bwt, &aux.k_minus_one_ranges, &aux.shorter_than_k
                                );
                            }

                            if !predecessor_confirmed {
                                keep_predecessors_atomic(predecessor, bwt, k, keep_dummy);
                            }
                        } else {
                            // Longer than k.
                            predecessor_confirmed = true;
                        }
                    }
                    index += 1;
                }
            });
        }
    });

    let keep_dummy = crate::util::bitvec_to_simple_sds_raw_bitvec(
        keep_dummy.into_bitvec_u64()
    );

    log::info!("[par_build_dummy_marks] done");
    DummyMarks {
        keep_dummy,
    }
}

/// Finds the end of the (k-1)-range and performs two rank queries on the shorter_than_k bitvector
/// in order to figure out if there is a non-dummy k-mer in this (k-1) range. If there is then the
/// dummy k-mer is not needed.
pub(crate) fn has_full_kmer_predecessor(
    predecessor: usize,
    bwt: &Bwt,
    k_minus_one_ranges: &BitVector, 
    shorter_than_k: &BitVector
) -> bool {
    let range_start = predecessor;
    let one_index = k_minus_one_ranges.rank(range_start + 1);
    let range_end = if one_index == k_minus_one_ranges.count_ones() {
        bwt.len()
    } else {
        // There is at least one 1 after the current position.
        k_minus_one_ranges.select(one_index).unwrap()
    };
    let range_length = range_end - range_start;
    let number_of_prefixes_with_true_length_smaller_than_k =
        shorter_than_k.rank(range_end) - shorter_than_k.rank(range_start);
    number_of_prefixes_with_true_length_smaller_than_k < range_length
}

fn keep_predecessors_atomic(mut predecessor: usize, bwt: &Bwt, mut k: usize, keep_dummy: &AtomicBitmap) {
    while k > 0 {
        keep_dummy.set(predecessor, true);
        let next = bwt.inverse_lf_step(predecessor);
        predecessor = next;
        k -= 1;
    }
}

pub fn build_with_all_dummies<SS: SubsetSeq + Send>(
    mut input: Vec<u8>, bounded_context_suffix_array: Vec<usize>, k: usize, build_lcs: bool, build_counts: bool
) -> Output<SS> {
    log::info!("[build_with_all_dummies] begin");
    let word_count = k.div_ceil(LANES);
    let length = bounded_context_suffix_array.len();
    pad_input(&mut input, k, length);

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
            find_length_and_lcp_value(
                k,
                &input,
                word_count,
                current_suffix_index,
                previous_suffix_index
            )
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

    let result = collect_output(k, rows, lcs, counts, kmer_count);
    log::info!("[build_with_all_dummies] done");
    result
}

#[inline]
fn push_set(rows: &mut [BitVec<u64>], set: u8) {
    for i in 1..=4 {
        rows[i - 1].push(set & (1 << i) != 0);
    }
}

#[inline]
fn include_letter(bwt: &Bwt, index: usize, current_set: u8) -> u8 {
    (1 << bwt.get_char_index(index)) as u8 | current_set
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{BitPackedKmerSortingMem, SubsetMatrix, VecSeqStream};

    struct RevVecSeqStream<'a> {
        seqs: &'a [Vec<u8>],
        copy: Vec<u8>,
        index: usize,
    }

    impl<'a> RevVecSeqStream<'a> {
        fn new(seqs: &'a [Vec<u8>]) -> Self {
            Self {
                seqs,
                copy: vec![],
                index: 0
            }
        }
    }

    impl<'a> crate::SeqStream for RevVecSeqStream<'a> {
        fn stream_next(&mut self) -> Option<&[u8]> {
            if self.index >= self.seqs.len() {
                return None;
            }
            self.copy.clear();
            self.copy.extend(&self.seqs[self.index]);
            preprocessing::sanitise(&mut self.copy);
            self.copy.reverse();
            let result = &self.copy;
            self.index += 1;
            Some(result)
        }
    }

    pub fn make_concatenation(seqs: &[Vec<u8>]) -> Vec<u8> {
        let mut stream = RevVecSeqStream::new(seqs);
        let capacity: usize = seqs.iter().map(|seq| seq.len() + 1).sum::<usize>() + 2;
        let mut concatenation = Vec::<u8>::with_capacity(capacity);
        preprocessing::concatenate_sequences(&mut stream, &mut concatenation).unwrap();
        concatenation
    }

    macro_rules! seqs {
        ($($seq:expr),* $(,)?) => {
            vec![$($seq.to_vec(),)*]
        };
    }

    #[test]
    fn concatenate_sequences() {
        let seqs = seqs![b"ACGT", b"ANOPT", b"TGCA"];
        let concatenation = make_concatenation(&seqs);
        assert_eq!(b"#$TGCA$T$$$A$ACGT$".as_slice(), &concatenation);
    }

    fn make_suffix_array(concatenation: &[u8]) -> Vec<usize> {
        let mut suffix_array = (0..concatenation.len()).collect::<Vec<_>>();
        suffix_array.sort_by_key(|index| &concatenation[*index..]);
        suffix_array
    }

    #[test]
    fn randomised_kmers() {
        use rand_chacha::ChaCha20Rng;
        use rand_chacha::rand_core::SeedableRng;
        use rand_chacha::rand_core::RngCore;

        let k: usize = 17;
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
        let suffix_array = make_suffix_array(&concatenation);

        {
            // Without redundant dummies.
            let (correct_sbwt, correct_lcs) = BitPackedKmerSortingMem::new_from_vecs(&seqs, k)
                .build_lcs(true)
                .run();

            // With bounded context suffix array.
            let Output {
                sbwt: constructed_sbwt,
                lcs: constructed_lcs,
                counts
            } = build_with_bounded_suffix_array::<SubsetMatrix>(4, concatenation.clone(), k, 4, true, false, false);

            let correct_lcs = correct_lcs.unwrap();
            let constructed_lcs = constructed_lcs.unwrap();

            assert!(counts.is_none());
            assert_eq!(correct_sbwt, constructed_sbwt);
            assert_eq!(correct_lcs, constructed_lcs);

            // With full context suffix array.
            let Output {
                sbwt: constructed_sbwt,
                lcs: constructed_lcs,
                counts
            } = build::<SubsetMatrix>(4, concatenation.clone(), suffix_array.clone(), k, true, false, false);

            let constructed_lcs = constructed_lcs.unwrap();
            assert!(counts.is_none());
            assert_eq!(correct_sbwt, constructed_sbwt);
            assert_eq!(correct_lcs, constructed_lcs);
        }

        {
            // With all dummies.
            let (mut correct_sbwt, correct_lcs) = BitPackedKmerSortingMem::new_from_vecs(&seqs, k)
                .build_lcs(true)
                .add_all_dummy_paths(true)
                .run();

            // With bounded context suffix array.
            let Output {
                sbwt: constructed_sbwt,
                lcs: constructed_lcs,
                counts
            } = build_with_bounded_suffix_array::<SubsetMatrix>(4, concatenation.clone(), k, 4, true, true, true);

            let correct_lcs = correct_lcs.unwrap();
            let constructed_lcs = constructed_lcs.unwrap();

            assert_eq!(correct_sbwt, constructed_sbwt);
            assert_eq!(correct_lcs.len(), constructed_lcs.len());
            assert_eq!(correct_lcs, constructed_lcs);

            correct_sbwt.build_select();
            let pnsv = crate::vodbg::pnsv::PnsvTuned::new_default(&correct_sbwt, &correct_lcs, k);
            let mut vodbg = crate::vodbg::VoDbg::new(&correct_sbwt, &pnsv);
            vodbg.build_counts(
                VecSeqStream::new(&seqs),
                true,
                Counts::DEFAULT_SAMPLE_DISTANCE,
                1, 4, Counts::DEFAULT_BATCH_SIZE_IN_BYTES).unwrap();
            assert_eq!(vodbg.counts().unwrap(), &counts.unwrap());

            // With full context suffix array.
            let Output {
                sbwt: mut constructed_sbwt,
                lcs: constructed_lcs,
                counts
            } = build::<SubsetMatrix>(4, concatenation, suffix_array, k, true, true, true);

            constructed_sbwt.build_select();
            let constructed_lcs = constructed_lcs.unwrap();

            assert_eq!(correct_sbwt, constructed_sbwt);
            assert_eq!(constructed_lcs, constructed_lcs);
            assert_eq!(vodbg.counts().unwrap(), &counts.unwrap());
        }
    }
}

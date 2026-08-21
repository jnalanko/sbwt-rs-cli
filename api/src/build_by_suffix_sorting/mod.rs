//! Construction algorithms for an [SbwtIndex] via suffix sorting. Good for large k.

// Module and submodule contributions by Martin Kostadinov.

//
// # Construction without redundant dummies
//
// There are a three steps to this version of the algorithm. Firstly some auxiliary arrays and
// bitvectors are constructed which are then used in the second step in which dummy k-mers found
// necessary are marked. The final step is the population of the sets in the SBWT index.
//
// ## Auxiliary data
//
// To describe the auxiliary data, here are a few definitions:
// 1. The "true prefix" of a given suffix is the prefix up to the first $. The "true length" of
//    the suffix is the length of the true prefix.
// 2. The "k-view prefix" (k-VP) of a given suffix is the true prefix of the suffix truncated from
//    the right to a length of k, if the true length is greater than k, or padded to the right
//    with $ if the true length is smaller than k. The k-length of a suffix is equal to the number
//    of non $ characters in the k-VP.
// 3. A dummy suffix is a suffix with true length less than k. A non-dummy suffix is a suffix with
//    true length greater than or equal to k.
// 4. A k-range is a sequence of suffixes in the order of the suffix array which share the same
//    k-VP. A (k-1)-range contains one or more k-ranges.
// 5. Given a k-mer A and B, if A[1..k] == B[0..k-1] then A is a predecessor to B and B is a
//    successor to A. Given the k-VP of two suffixes A and B, if A[0..k-1] == B[1..k], then the
//    suffix of A is a predecessor to the suffix of B and the suffix of B is a successor to the
//    suffix of A. This is because, the input sequences are reversed in the input concatenation.
//    The (k-1)-range which contains the predecoessor suffixes of another given suffix will be
//    referred to as the (k-1)-predecessor range.
//
// Auxiliary data structures:
// 1. The BWT is a string of characters. It gives the character before the start of each suffix.
//    Since the suffixes are sorted lexicographically and the k-mers in the SBWT index are sorted
//    colexicographically, the concatenation of the input sequences has them reversed. This means
//    that the character before a given suffix is one of the outedges in the SBWT graph of the
//    first k-mer in the SBWT index with a given (k-1) suffix. In order to support standard FM
//    index operations it is represented by 5 bitvectors - one for each of the characters in the
//    set {$, A, C, G, T}.
// 2. The LCP array gives the longest common prefix between two consecutive suffixes in the order
//    of the suffix array.
// 3. The Lengths array gives the (k+1)-lengths of all suffixes. This is in order to mark only
//    suffixes which have a true length which is exactly equal to k.
// 4. The [FullAuxiliaryData::shorter_than_k] bitvector marks suffixes with true length shorter
//    than k. Used in the step of marking which dummy suffixes and their corresponding dummy
//    k-mers are necessary.
// 5. The [FullAuxiliaryData::equal_to_k_minus_one_or_k] bitvector marks suffixes with true length
//    equal to (k-1) or k. Together with the previous bitvector it can be used to disambiguate
//    whether the true length of a suffix is (k-1) or k. If a suffix's true length is equal to k,
//    then the algorithm must check whether this suffix has a non-dummy suffix in the
//    (k-1)-predecessor range of the k-mer which the suffix represents.
// 6. The [FullAuxiliaryData::k_minus_one_ranges] and the [FullAuxiliaryData::k_ranges] bitvectors
//    mark the start of all (k-1)-ranges and k-ranges respectively. The former supports rank and
//    select operations in order to find the end of a (k-1)-range. A start of a k-range can be
//    identified if the k-length of a suffix is greater than its LCP value. This means that the
//    k-VP of this suffix is different from the k-VP of the previous suffix. For a suffix to be
//    the start of a (k-1)-range it must be the start of a k-range and its LCP value must be less
//    than (k-1). If its LCP value is greater than or equal to (k-1), then that means that it
//    shares a (k-1)-VP with the previous suffix and thus it is not the start of a (k-1)-range.
//
// ## Marking dummy k-mers
//
// If a suffix has a true length exactly equal to k, this means that it is the beginning k-mer in
// an input sequence. Such suffixes will be referred to as "start suffixes". That is, if for a
// given start suffix there doesn't exist a non-dummy predecessor suffix, we need to mark dummy
// suffixes such that when we populate the sets of the SBWT index we create a set for the dummy
// predecessors of this suffix's k-mer.
//
// For any given k-range we need to confirm whether a non-dummy predecessor suffix exists only
// once since all suffixes in this range will correspond to only one k-mer in the final SBWT. 
//
// A dummy k-mer is not needed if the k-range contains a suffix with true length greater than k.
// To show this is true, let's assume that the index of the starting character of the suffix with
// true length greater than k is I. Then the suffix starting at index I+1 has a true length equal
// to k, which makes it a non-dummy predecessor to the suffix I (no matter where it is located in
// the suffix array).
//
// A dummy k-mer is also not needed if the (k-1)-predecessor range of a given suffix contains a
// non-dummy k-mer. Here are a few prerequisites before describing the operation that can check
// this. It can be proven that an "inverse LF step" from the index in the suffix array of a given
// suffix on a k-BWT will return an index in its (k-1)-predecessor range. Given any (k-1)-range,
// if there are non-dummy suffixes in it, they will be at the end of the range due to the
// lexicographic sorting of the suffixes. Thus, given a start suffix, an inverse LF step is
// performed to find its (k-1)-predecessor range. Then, using standard bitvector operations on
// [FullAuxiliaryData::k_minus_one_ranges] the end of that range is found. Then if that range
// contains at least one non-dummy suffix, the last bit of the (k-1)-range in the
// [FullAuxiliaryData::shorter_than_k] bitvector will not be set.
//
// Otherwise, the direct dummy suffix predecessor and all of its predecessors need to be marked as
// well. This is done by performing inverse LF steps until the k-range of the root dummy node
// ($$..$). Since the inverse LF step returns the index of a predecessor suffix, this means that
// the character in the BWT at the index of the predecessor is equal to the first character in the
// successor suffix.
//
// ## Populating the sets of the SBWT index.
//
// By definition of the SBWT index, a set contains characters if the corresponding k-mer is the
// first k-mer with a given (k-1)-suffix. As such, the characters for a set are collected by
// iterating an entire (k-1)-range. For each (k-1)-range a number of sets are emitted equal to the
// number of k-ranges it contains. The first set contains the encountered characters, the rest are
// empty. That number corresponds to different k-mers which will potentially end up in the SBWT. A
// k-range that corresponds to non-dummy suffixes is always counted. A k-range that corresponds to
// dummy suffixes is counted only if there are marked suffixes. Suffixes with true length shorter
// than k-1 contribute their letter to a set only if they are marked to be kept. All other
// suffixes contribute their letter to the first emitted set.
//
// ## Generating the LCS array
//
// The LCS array can be induced from the LCP array. Each entry in it is equal to the range minimum
// from the previous kept entry. The beginning of each k-range is kept if it corresponds to dummy
// suffixes at least one of which needs to be kept or if it corresponds to non-dummy suffixes.
//
// # Construction with all dummies included
//
// This version of the algorithm is simpler. After calculating the Lengths and the LCP arrays it
// can directly populate the sets of the SBWT index. This step is the same as the other version of
// the algorithm escept that all k-ranges that correspond to dummy suffixes are kept.
//
// The Counts data structure can be created easily by counting the number of suffixes in a k-range.
//

pub mod preprocessing;
pub mod input_structures;

use input_structures::{Bwt, Lcp, CHAR_TO_INDEX};
use crate::atomic_bitmap::AtomicBitmap;
use crate::vodbg::count::{Counts, Sample};
use crate::{LcsArray, SbwtIndex, SubsetSeq};

use std::sync::{Arc, Mutex, atomic::AtomicUsize};

use bitvec::vec::BitVec;
use bitvec::field::BitField;
use simple_sds_sbwt::{bit_vector::BitVector, int_vector::IntVector};
use simple_sds_sbwt::ops::{BitVec as BitVecTrait, Rank, Select};
use simple_sds_sbwt::raw_vector::{AccessRaw, RawVector};
use simple_sds_sbwt::serialize::Serialize;

type Word = wide::u8x16;
const LANES: usize = Word::LANES as usize;
type Bitmask = u16;
const _: () = assert!(Bitmask::BITS as usize == LANES);
const FULL_SET: u8 = 0b00011110;

///
/// The result of constructing the SBWT data structure using a suffix array. The suffix array can
/// either be constructed using traditional algorithms such as SA-IS or the suffixes can be sorted
/// by a bounded prefix i.e. up to the first k-characters since that is what is strictly necessary
/// for the SBWT.
///
pub struct Output<SS: SubsetSeq + Send> {
    /// The result SBWT index.
    pub sbwt: SbwtIndex<SS>,
    /// During construction the LCS array can optionally be constructed. It is independent of
    /// whether all dummies are included or not. 
    pub lcs: Option<LcsArray>,
    /// If the SBWT index is constructed by including all dummies, then the counts of all k-mers
    /// can automatically be added during construction with the corresponding version of the
    /// algorithm.
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
            let output = par_build(threads, input, suffix_array, k, build_lcs, add_all_dummies, build_counts, false);
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
            let output = par_build::<SS>(threads, input, suffix_array, k, build_lcs, add_all_dummies, build_counts, true);
            result = Some(output);
        });
    });
    log::info!("[build_with_bounded_suffix_array] done");
    result.unwrap()
}

#[allow(clippy::too_many_arguments)]
pub fn par_build<SS: SubsetSeq + Send>(
    threads: usize,
    input: Vec<u8>,
    suffix_array: Vec<usize>,
    k: usize,
    build_lcs: bool,
    add_all_dummies: bool,
    build_counts: bool,
    is_bounded_suffix_array: bool,
) -> Output<SS> {
    log::info!(
        "[par_build] begin [length: {} | build_lcs: {} | add_all_dummies: {} | build_counts {}]",
        input.len(), build_lcs, add_all_dummies, build_counts
    );
    let output = if !add_all_dummies {
        par_build_without_redundant_dummies(threads, input, suffix_array, k, build_lcs, is_bounded_suffix_array)
    } else {
        par_build_with_all_dummies(threads, input, suffix_array, k, build_lcs, build_counts, is_bounded_suffix_array)
    };
    log::info!("[par_build] done");
    output
}

pub fn par_build_without_redundant_dummies<SS: SubsetSeq + Send>(
    threads: usize,
    input: Vec<u8>,
    bounded_context_suffix_array: Vec<usize>,
    k: usize,
    build_lcs: bool,
    is_bounded_suffix_array: bool
) -> Output<SS> {
    log::info!("[par_build_without_redundant_dummies] begin");
    let _ = threads;
    let aux = par_build_full_auxiliary_data(threads, input, bounded_context_suffix_array, k, is_bounded_suffix_array);
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

#[inline]
fn _par_build_without_redundant_dummies(
    threads: usize,
    k: usize,
    aux: FullAuxiliaryData,
    dummy_marks: DummyMarks,
    build_lcs: bool,
) -> (Vec<BitVec<u64>>, Option<IntVector>) {
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
                let result = _build_sbwt_region_without_redundant_dummies(
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

    let mut rows_parts: Vec<Vec<BitVec<u64>>> = vec![
        Vec::with_capacity(results.len()),
        Vec::with_capacity(results.len()),
        Vec::with_capacity(results.len()),
        Vec::with_capacity(results.len()),
    ];
    for row_index in 0..4 {
        rows_parts[row_index].push(bitvec::bitvec![u64, bitvec::order::Lsb0; 0; 1]);
    }

    let mut lcs_parts = {
        let capacity = if build_lcs {
            results.len()
        } else {
            0
        };
        Vec::<(Option<BitVec<u64>>, usize)>::with_capacity(capacity)
    };

    for result in results {
        let SbwtRegionResult { rows, lcs, last_lcs_value, .. } = result;
        if rows.is_empty() {
            continue;
        }
        if build_lcs {
            lcs_parts.push((lcs, last_lcs_value));
        }
        for (row_index, row) in rows.into_iter().enumerate() {
            rows_parts[row_index].push(row);
        }
    }

    log::info!("[_par_build_without_redundant_dummies] agglomerating result");
    let rows = Arc::new(Mutex::new(Vec::<(usize, BitVec<u64>)>::new()));
    let mut lcs = None;
    rayon::scope(|s| {
        let rows = &rows;
        let lcs  = &mut lcs;

        if build_lcs {
            s.spawn(move |_| {
                let bit_width = usize::BITS - (k.overflowing_sub(1).0).leading_zeros();
                let bit_width = bit_width as usize;

                let mut last_lcs_value: usize = 0;
                let mut bitvec_lcs_parts = vec![];
                bitvec_lcs_parts.push(bitvec::bitvec![u64, bitvec::order::Lsb0; 0; bit_width]);
                for (lcs_part, current_last_lcs_value) in lcs_parts {
                    let mut lcs_part = lcs_part.unwrap();
                    if lcs_part.is_empty() {
                        last_lcs_value = last_lcs_value.min(current_last_lcs_value);
                        continue;
                    }
                    let (first_value_slice, _) = lcs_part.split_at_mut(bit_width);
                    let first_value: usize = first_value_slice.load_le();
                    let new_first_value = first_value.min(last_lcs_value);
                    first_value_slice.store_le(new_first_value);
                    bitvec_lcs_parts.push(lcs_part);
                    last_lcs_value = current_last_lcs_value;
                }

                let value = par_concatenate_lcs(bit_width, bitvec_lcs_parts);
                *lcs = Some(value);
            });
        }

        for (row_index, row_parts) in rows_parts.into_iter().enumerate() {
            s.spawn(move |_| {
                let local_row = crate::util::parallel_bitvec_concat(row_parts);
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
    lcs: Option<BitVec<u64>>,
    last_lcs_value: usize,
}

fn _build_sbwt_region_without_redundant_dummies(
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

    let mut index = start;
    let mut current_lcs_value  = k - 1;

    while index < end {
        if build_lcs {
            current_lcs_value = current_lcs_value.min(aux.lcp.get(index));
        }
        if aux.k_minus_one_ranges.get(index) {
            break;
        }
        index += 1;
    }

    if index >= end {
        return SbwtRegionResult {
            start,
            rows: vec![],
            lcs: None,
            last_lcs_value: current_lcs_value
        };
    }

    let mut rows = Vec::<BitVec<u64>>::new();
    for _ in 0..4 {
        rows.push(BitVec::with_capacity(capacity));
    }

    let bit_width = usize::BITS - (k.overflowing_sub(1).0).leading_zeros();
    let bit_width = bit_width as usize;
    let mut lcs_count: usize = 0;
    let mut lcs: Option<BitVec<u64>> = if build_lcs {
        let value = bitvec::bitvec![u64, bitvec::order::Lsb0; 0; bit_width * capacity];
        Some(value)
    } else {
        None
    };

    let mut current_set = 0;
    let mut include_dummy_kmer = false;
    let mut has_dummy_kmer     = false;
    let mut k_range_count = 0;

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

        if aux.shorter_than_k.bit(index) {
            has_dummy_kmer = true;
            if dummy_marks.keep_dummy.bit(index) {
                if build_lcs && !include_dummy_kmer {
                    let (_, suffix) = lcs.as_mut().unwrap().split_at_mut(lcs_count * bit_width);
                    let (word, _)   = suffix.split_at_mut(bit_width);
                    word.store_le(current_lcs_value as u64);
                    lcs_count += 1;
                    current_lcs_value = k - 1;
                }
                include_dummy_kmer = true;
            }
            if dummy_marks.keep_dummy.bit(index) || aux.equal_to_k_minus_one_or_k.bit(index) {
                current_set = include_letter(&aux.bwt, index, current_set);
            }
        } else {
            if build_lcs && is_start_of_k_range {
                let (_, suffix) = lcs.as_mut().unwrap().split_at_mut(lcs_count * bit_width);
                let (word, _)   = suffix.split_at_mut(bit_width);
                word.store_le(current_lcs_value as u64);
                lcs_count += 1;
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

    if build_lcs {
        lcs.as_mut().unwrap().truncate(bit_width * lcs_count);
    }

    SbwtRegionResult {
        start,
        rows,
        lcs,
        last_lcs_value: current_lcs_value
    }
}

pub(crate) fn par_concatenate_lcs(bit_width: usize, bitvec_lcs_parts: Vec<BitVec<u64>>) -> IntVector {
    let mut bitvec_lcs = crate::util::parallel_bitvec_concat(bitvec_lcs_parts);

    let mut int_vector_header = [0u64, 0u64]; // len, width
    int_vector_header[0] = (bitvec_lcs.len() / bit_width) as u64;
    int_vector_header[1] = bit_width as u64;
    let mut raw_vector_header = [0u64, 0u64]; // bits, words
    raw_vector_header[0] = bitvec_lcs.len() as u64;
    raw_vector_header[1] = bitvec_lcs.len().div_ceil(64) as u64;

    let int_vector_header_bytes: &[u8] = bytemuck::cast_slice(&int_vector_header);
    let raw_vector_header_bytes: &[u8] = bytemuck::cast_slice(&raw_vector_header);

    let original_len = bitvec_lcs.len();
    bitvec_lcs.resize(original_len.next_multiple_of(64), false);
    bitvec_lcs.resize(original_len, false);

    let raw_data: &[u8] = bytemuck::cast_slice(bitvec_lcs.as_raw_slice());
    use std::io::{Cursor, Read};
    let mut data_with_headers = Cursor::new(int_vector_header_bytes)
        .chain(Cursor::new(raw_vector_header_bytes))
        .chain(Cursor::new(raw_data));

    IntVector::load(&mut data_with_headers).unwrap()
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
    let word_count = word_count(k);
    if input.len() > length {
        return;
    }
    for _ in 1..(3 * word_count * LANES) {
        input.push(b'$');
    }
}

pub fn par_bounded_context_suffix_array_bucket_sort(input: &mut Vec<u8>, k: usize, prefix_length_for_bucket_sort: usize) -> Vec<usize> {
    log::info!("[par_bounded_context_suffix_array_bucket_sort] begin");
    let length = input.len();
    let word_count = word_count(k);
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
    pub(crate) shorter_than_k: RawVector,
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

fn par_build_full_auxiliary_data(
    threads: usize,
    mut input: Vec<u8>,
    suffix_array: Vec<usize>,
    k: usize,
    is_bounded_suffix_array: bool,
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

    let length = suffix_array.len();
    pad_input(&mut input, k, length);

    let lcp = par_build_lcp(threads, &input, &suffix_array, k, is_bounded_suffix_array);
    let lengths = par_build_lengths(threads, &input, &suffix_array, k);

    log::info!("[par_build_full_auxiliary_data] done with LCP and Lengths");
    let bwt_vectors = [
        AtomicBitmap::new(length), // $
        AtomicBitmap::new(length), // A
        AtomicBitmap::new(length), // C
        AtomicBitmap::new(length), // G
        AtomicBitmap::new(length), // T
    ];
    bwt_vectors[0].set(0, true);

    let kmer_count = AtomicUsize::new(0);
    let shorter_than_k     = AtomicBitmap::new(length);
    let equal_to_k_minus_one_or_k         = AtomicBitmap::new(length);
    let k_minus_one_ranges = AtomicBitmap::new(length);
    let k_ranges           = AtomicBitmap::new(length);

    let region_size = length.div_ceil(threads);

    rayon::scope(|s| {
        let input              = &input;
        let suffix_array       = &suffix_array;
        let lcp                = &lcp;
        let lengths            = &lengths;
        let bwt_vectors        = &bwt_vectors;
        let kmer_count         = &kmer_count;
        let shorter_than_k     = &shorter_than_k;
        let equal_to_k_minus_one_or_k         = &equal_to_k_minus_one_or_k;
        let k_minus_one_ranges = &k_minus_one_ranges;
        let k_ranges           = &k_ranges;

        for thread_index in 0..threads {
            s.spawn(move |_| {
                let start = 1 + thread_index * region_size;
                let end = (start + region_size).min(length);

                for index in start..end {
                    let current_suffix_index = suffix_array[index];

                    // Set the bit in the bitvector corresponding to the previous character in the input
                    // for the (bounded) BWT.
                    let previous_character_index_in_input = (current_suffix_index + length - 1) % length;
                    let previous_character = input[previous_character_index_in_input];
                    let previous_character_index = CHAR_TO_INDEX[previous_character as usize];
                    if previous_character_index < bwt_vectors.len() {
                        bwt_vectors[previous_character_index].set(index, true);
                    }

                    let length = lengths.get(index);
                    let lcp_value = lcp.get(index);

                    if length < k {
                        shorter_than_k.set(index, true);
                        if length == k - 1 {
                            equal_to_k_minus_one_or_k.set(index, true);
                        }
                    } else if length == k {
                        equal_to_k_minus_one_or_k.set(index, true);
                    }

                    if lcp_value < length.min(k) {
                        k_ranges.set(index, true);
                        if length >= k {
                            kmer_count.fetch_add(1, std::sync::atomic::Ordering::AcqRel);
                        }
                        if length < k || lcp_value < k - 1 {
                            k_minus_one_ranges.set(index, true);
                        }
                    }
                }
            });
        }
    });

    use rayon::prelude::*;

    fn convert_atomic_bitmap(value: AtomicBitmap) -> RawVector {
        let intermediate = value.into_bitvec_u64();
        crate::util::bitvec_to_simple_sds_raw_bitvec(intermediate)
    }

    let mut bwt_op = None;
    let mut shorter_than_k_op            = None;
    let mut equal_to_k_minus_one_or_k_op = None;
    let mut k_minus_one_ranges_op        = None;
    let mut k_ranges_op                  = None;

    rayon::scope(|s| {
        s.spawn(|_| {
            log::info!("[par_build_full_auxiliary_data] bwt begin");
            let bwt_bit_vectors = bwt_vectors.into_par_iter()
                .map(|value| {
                    let intermediate = convert_atomic_bitmap(value);
                    BitVector::from(intermediate)
                })
                .collect::<Vec<_>>();
            let bwt = Bwt::new(bwt_bit_vectors);
            bwt_op = Some(bwt);
            log::info!("[par_build_full_auxiliary_data] bwt done");
        });

        s.spawn(|_| {
            let shorter_than_k = convert_atomic_bitmap(shorter_than_k);
            shorter_than_k_op = Some(shorter_than_k);
        });

        s.spawn(|_| {
            let equal_to_k_minus_one_or_k = convert_atomic_bitmap(equal_to_k_minus_one_or_k);
            equal_to_k_minus_one_or_k_op = Some(equal_to_k_minus_one_or_k);
        });

        s.spawn(|_| {
            let k_ranges = convert_atomic_bitmap(k_ranges);
            k_ranges_op = Some(k_ranges);
        });

        s.spawn(|_| {
            log::info!("[par_build_full_auxiliary_data] k_minus_one_ranges begin");
            let k_minus_one_ranges = convert_atomic_bitmap(k_minus_one_ranges);
            let mut k_minus_one_ranges = BitVector::from(k_minus_one_ranges);
            k_minus_one_ranges.enable_rank();
            k_minus_one_ranges.enable_select();
            k_minus_one_ranges_op = Some(k_minus_one_ranges);
            log::info!("[par_build_full_auxiliary_data] k_minus_one_ranges done");
        });
    });

    FullAuxiliaryData {
        kmer_count: kmer_count.load(std::sync::atomic::Ordering::Relaxed),
        bwt: bwt_op.unwrap(),
        lcp,
        shorter_than_k            : shorter_than_k_op.unwrap(),
        equal_to_k_minus_one_or_k : equal_to_k_minus_one_or_k_op.unwrap(),
        k_minus_one_ranges        : k_minus_one_ranges_op.unwrap(),
        k_ranges                  : k_ranges_op.unwrap(),
    }
}

#[allow(unused)]
fn par_build_lcp(
    threads: usize,
    input: &[u8],
    suffix_array: &[usize],
    k: usize,
    is_bounded_suffix_array: bool,
) -> Lcp {
    log::info!("[par_build_lcp] begin");
    use std::sync::atomic::Ordering;

    let length = suffix_array.len();
    let phi: Vec<AtomicUsize> = vec![0_usize; length]
        .into_iter()
        .map(AtomicUsize::new)
        .collect();

    let thread_region_length = length.div_ceil(threads);

    log::info!("[par_build_lcp] building phi array");
    rayon::scope(|s| {
        let phi = &phi;

        for thread_index in 0..threads {
            let mut start = thread_index * thread_region_length;
            let end = (start + thread_region_length).min(length);
            s.spawn(move |_| {
                if start == 0 {
                    start += 1;
                }
                for j in start..end {
                    let phi_index = suffix_array[j];
                    let value = suffix_array[j - 1];
                    phi[phi_index].store(value, Ordering::Release);
                }
            });
        }
    });

    let phi: Vec<usize> = phi
        .into_iter()
        .map(|value| value.load(Ordering::Relaxed))
        .collect();

    log::info!("[par_build_lcp] building plcp parts");
    let word_count = word_count(k);
    let lcp_width = byte_width(k);
    let plcp_parts = Arc::new(Mutex::new(Vec::<(usize, Lcp)>::with_capacity(threads)));
    rayon::scope(|s| {
        let lcp_result = &plcp_parts;
        let phi = &phi;

        for thread_index in 0..threads {
            let mut start = thread_index * thread_region_length;
            let end = (start + thread_region_length).min(length);
            s.spawn(move |_| {
                let local_region_length = end - start;
                let mut local_lcp = {
                    let mut capacity = local_region_length * lcp_width;
                    let data = Vec::<u8>::with_capacity(capacity);
                    Lcp::new_with_width(data, lcp_width)
                };

                if start == 0 {
                    local_lcp.push(0);
                    start += 1;
                }

                let mut previous_lcp_value: usize = 0;
                for i in start..end {
                    if phi[i] == 0 {
                        local_lcp.push(0);
                        previous_lcp_value = 0;
                        continue;
                    }

                    let start_lcp_value = if is_bounded_suffix_array {
                        0
                    } else {
                        previous_lcp_value.saturating_sub(1)
                    };

                    let (lcp_value, length) = find_lcp_value(
                        k,
                        input,
                        word_count,
                        start_lcp_value,
                        i,
                        phi[i]
                    );
                    local_lcp.push(lcp_value);
                    previous_lcp_value = lcp_value;
                }

                lcp_result.lock().unwrap().push((thread_index, local_lcp));
            });
        }
    });

    drop(phi);

    log::info!("[par_build_lcp] reorder:");
    let plcp_parts = {
        let mutex = Arc::try_unwrap(plcp_parts).unwrap();
        let mut unordered_plcp_parts = mutex.into_inner().unwrap();
        unordered_plcp_parts.sort_by_key(|(index, _)| *index);
        unordered_plcp_parts.into_iter().map(|(_, part)| part).collect::<Vec<_>>()
    };
    let lcp = par_reorder_parts(threads, suffix_array, lcp_width, plcp_parts);

    log::info!("[par_build_lcp] done");
    lcp
}

#[allow(unused)]
fn par_build_lengths(
    threads: usize,
    input: &[u8],
    suffix_array: &[usize],
    k: usize,
) -> Lcp {
    log::info!("[par_build_lengths] begin");
    #[derive(Debug)]
    struct LengthPart {
        start: usize,
        skipped_entries: usize,
        values: Lcp,
    }

    let lengths_width = byte_width(k + 1);
    let length = suffix_array.len();
    let thread_region_length = length.div_ceil(threads);
    let parts = Arc::new(Mutex::new(Vec::<LengthPart>::with_capacity(threads)));
    rayon::scope(|s| {
        let parts = &parts;
        for thread_index in 0..threads {
            let start = thread_index * thread_region_length;
            let end = (start + thread_region_length).min(length);

            s.spawn(move |_| {
                let mut skip_end = end;
                let mut current_length = 0;
                while skip_end > start && current_length < k + 1 {
                    skip_end -= 1;
                    if !crate::util::is_dna(input[skip_end]) {
                        current_length = 0;
                        break;
                    } else {
                        current_length += 1;
                    }
                }
                let skipped_entries = end - skip_end - 1;
                let local_buffer_length = end - start;
                let data = vec![0_u8; local_buffer_length];
                let mut local_lengths = Lcp::new_with_width(data, lengths_width);
                local_lengths.set(skip_end - start, current_length);
                for i in (start..skip_end).rev() {
                    if crate::util::is_dna(input[i]) {
                        current_length += 1;
                        current_length = current_length.min(k + 1);
                    } else {
                        current_length = 0;
                    }
                    local_lengths.set(i - start, current_length);
                }

                parts.lock().unwrap().push(LengthPart {
                    start,
                    skipped_entries,
                    values: local_lengths,
                });
            });
        }
    });

    log::info!("[par_build_lengths] fixing border values");
    let mut parts = {
        let mutex = Arc::try_unwrap(parts).unwrap();
        let mut unordered_parts = mutex.into_inner().unwrap();
        unordered_parts.sort_by_key(|part| part.start);
        unordered_parts
    };

    let mut border_value = 0;
    for part in parts.iter_mut().rev() {
        let start = part.start;
        let end = start + part.values.len();
        let skipped_entries_start = end - part.skipped_entries;
        let mut current_length = border_value;
        for i in (skipped_entries_start..end).rev() {
            if crate::util::is_dna(input[i]) {
                current_length += 1;
                current_length = current_length.min(k + 1);
            } else {
                current_length = 0;
            }
            part.values.set(i - start, current_length);
        }
        border_value = part.values.get(0);
    }

    log::info!("[par_build_lengths] reorder:");
    let parts = parts.into_iter().map(|region| region.values).collect::<Vec<_>>();
    let lengths = par_reorder_parts(threads, suffix_array, lengths_width, parts);

    log::info!("[par_build_lengths] done");
    lengths
}

fn par_reorder_parts(threads: usize, suffix_array: &[usize], width: usize, parts: Vec<Lcp>) -> Lcp {
    let length = suffix_array.len();
    let thread_region_length = length.div_ceil(threads);
    let mut buffer = vec![0_u8; length * width];
    par_concatenate_bytes(
        "par_reorder_parts: permuted to buffer",
        &mut buffer,
        parts.iter()
    );
    let permuted = Lcp::new_with_width(buffer, width);

    log::info!("[par_reorder_parts] copying permuted parts");
    let reordered_parts = Arc::new(Mutex::new(Vec::<(usize, Lcp)>::with_capacity(threads)));
    rayon::scope(|s| {
        let permuted = &permuted;
        let reordered_parts = &reordered_parts;

        for (thread_index, mut part) in parts.into_iter().enumerate() {
            let start = thread_index * thread_region_length;
            let end = (start + thread_region_length).min(length);
            s.spawn(move |_| {
                for i in start..end {
                    let value = permuted.get(suffix_array[i]);
                    part.set(i - start, value);
                }

                reordered_parts.lock().unwrap().push((thread_index, part));
            })
        }
    });

    let reordered_parts = {
        let mutex = Arc::try_unwrap(reordered_parts).unwrap();
        let mut sorted_reordered_parts = mutex.into_inner().unwrap();
        sorted_reordered_parts.sort_by_key(|(index, _)| *index);
        sorted_reordered_parts.into_iter().map(|(_, part)| part).collect::<Vec<_>>()
    };
    let mut result_buffer: Vec<u8> = permuted.into();
    par_concatenate_bytes(
        "par_reorder_parts: reordered to buffer",
        &mut result_buffer,
        reordered_parts.iter()
    );
    let result = Lcp::new_with_width(result_buffer, width);
    drop(reordered_parts);

    log::info!("[par_reorder_parts] done");
    result
}


#[inline]
fn byte_width(value: usize) -> usize {
    let bit_width = usize::BITS - value.leading_zeros();
    (bit_width.div_ceil(u8::BITS) as usize).next_power_of_two()
}

#[inline]
fn word_count(k: usize) -> usize {
    k.div_ceil(LANES)
}

fn par_concatenate_bytes<S, I>(context: &str, buffer: &mut [u8], sections: S)
where S: Iterator<Item = I> + Send, I: AsRef<[u8]> + Send,
{
    log::info!("[{}] concatenation begin", context);
    let mut slice = buffer;
    rayon::scope(|s| {
        for section in sections {
            let len = section.as_ref().len();
            let (destination_slice, the_rest) = slice.split_at_mut(len);
            slice = the_rest;
            s.spawn(move |_| {
                let section = section.as_ref();
                destination_slice.copy_from_slice(section);
            });
        }
    });
    log::info!("[{}] concatenation done", context);
}

#[inline]
fn find_lcp_value(
    k: usize,
    input: &[u8],
    mut word_count: usize,
    start_lcp_value: usize,
    mut current_index: usize,
    mut previous_index: usize
) -> (usize, usize) {
    word_count = (word_count * LANES - start_lcp_value).div_ceil(LANES);
    current_index += start_lcp_value;
    previous_index += start_lcp_value;

    let k = k as u32;
    let mut lcp_value = start_lcp_value as u32;
    let mut length    = lcp_value;
    let mut stop_accumulating_length = false;
    for _ in 0..word_count {
        let slice_for_current = &input[current_index..current_index + LANES];
        current_index += LANES;
        let simd_word_current = Word::new(slice_for_current.try_into().unwrap());

        let slice_for_previous = &input[previous_index..previous_index + LANES];
        previous_index += LANES;
        let simd_word_previous = Word::new(slice_for_previous.try_into().unwrap());

        if !stop_accumulating_length {
            let bitmask = simd_word_current.simd_eq(b'$').to_bitmask() as Bitmask;
            length += bitmask.trailing_zeros();
            stop_accumulating_length |= length >= k || bitmask != 0;
        }

        let bitmask = simd_word_current.simd_eq(simd_word_previous).to_bitmask() as Bitmask;
        lcp_value += bitmask.trailing_ones();
        if bitmask != Bitmask::MAX {
            break;
        }
    }

    length = length.min(k);
    (lcp_value.min(length) as usize, length as usize)
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
    let included_k_minus_one_range = AtomicBitmap::new(length);

    let scan_length = length - scan_start;
    let thread_region_length = scan_length.div_ceil(threads);

    rayon::scope(|s| {
        let keep_dummy  = &keep_dummy;
        let included_k_minus_one_range = &included_k_minus_one_range;

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

                    if !aux.shorter_than_k.bit(index) {
                        if aux.equal_to_k_minus_one_or_k.bit(index) {
                            // Equal to k.
                            let predecessor = bwt.inverse_lf_step(index);
                            if !predecessor_confirmed {
                                predecessor_confirmed |= has_full_kmer_predecessor(
                                    predecessor, bwt, &aux.k_minus_one_ranges, &aux.shorter_than_k
                                );
                            }

                            if !predecessor_confirmed {
                                keep_predecessors_atomic(
                                    predecessor,
                                    aux,
                                    k,
                                    keep_dummy,
                                    included_k_minus_one_range
                                );
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
    shorter_than_k: &RawVector
) -> bool {
    let range_start = predecessor;
    let one_index = k_minus_one_ranges.rank(range_start + 1);
    let range_end = if one_index == k_minus_one_ranges.count_ones() {
        bwt.len()
    } else {
        // There is at least one 1 after the current position.
        k_minus_one_ranges.select(one_index).unwrap()
    };
    !shorter_than_k.bit(range_end - 1)
}

fn keep_predecessors_atomic(
    mut predecessor: usize,
    aux: &FullAuxiliaryData,
    mut k: usize,
    keep_dummy: &AtomicBitmap,
    included_k_minus_one_range: &AtomicBitmap,
) {
    let previous_k_minus_one_range_start_index = aux.k_minus_one_ranges.rank(predecessor);
    if previous_k_minus_one_range_start_index < aux.k_minus_one_ranges.count_ones() {
        // Bit will exist.
        let start_of_k_minus_one_range = aux.k_minus_one_ranges.select(previous_k_minus_one_range_start_index).unwrap();
        if included_k_minus_one_range.get(start_of_k_minus_one_range) {
            return;
        } else {
            included_k_minus_one_range.set(start_of_k_minus_one_range, true);
        }
    }
    while k > 0 {
        keep_dummy.set(predecessor, true);
        let next = aux.bwt.inverse_lf_step(predecessor);
        predecessor = next;
        k -= 1;
    }
}

pub fn par_build_with_all_dummies<SS: SubsetSeq + Send>(
    threads: usize,
    mut input: Vec<u8>,
    suffix_array: Vec<usize>,
    k: usize,
    build_lcs: bool,
    build_counts: bool,
    is_bounded_suffix_array: bool,
) -> Output<SS> {
    log::info!("[par_build_with_all_dummies] begin");
    let length = suffix_array.len();
    pad_input(&mut input, k, length);

    let lcp = par_build_lcp(threads, &input, &suffix_array, k, is_bounded_suffix_array);
    let lengths = par_build_lengths(threads, &input, &suffix_array, k);

    let length = suffix_array.len();
    let results = Arc::new(Mutex::new(Vec::<SbwtAllDummiesRegionResult>::with_capacity(threads)));
    let threads_total_length = length - 1;
    let thread_region_length = threads_total_length.div_ceil(threads);

    log::info!("[par_build_with_all_dummies] regions begin");
    rayon::scope(|s| {
        let input         = &input;
        let suffix_array  = &suffix_array;
        let lcp           = &lcp;
        let lengths       = &lengths;
        let results       = &results;

        for thread_index in 0..threads {
            s.spawn(move |_| {
                let start = 1 + thread_index * thread_region_length;
                let end = (start + thread_region_length).min(length);
                let is_last_thread = thread_index == threads - 1;
                let result = _build_sbwt_region_with_all_dummies(
                    k,
                    start,
                    end,
                    is_last_thread,
                    input,
                    suffix_array,
                    lcp,
                    lengths,
                    build_lcs,
                    build_counts
                );
                results.lock().unwrap().push(result);
            });
        }
    });
    log::info!("[par_build_with_all_dummies] regions done");

    let mut results = {
        let mutex = Arc::try_unwrap(results).unwrap();
        mutex.into_inner().unwrap()
    };
    results.sort_by_key(|result| result.start);

    let mut total_kmer_count = 0;

    let mut rows_parts: Vec<Vec<BitVec<u64>>> = vec![
        Vec::with_capacity(results.len()),
        Vec::with_capacity(results.len()),
        Vec::with_capacity(results.len()),
        Vec::with_capacity(results.len()),
    ];
    let mut lcs_parts = {
        let capacity = if build_lcs {
            results.len()
        } else {
            0
        };
        Vec::<Option<BitVec<u64>>>::with_capacity(capacity)
    };
    let mut counts_parts = {
        let capacity = if build_counts {
            results.len()
        } else {
            0
        };
        Vec::<Option<(Vec<u8>, Vec<u64>)>>::with_capacity(capacity)
    };

    for result in results {
        let SbwtAllDummiesRegionResult { kmer_count, rows, lcs, counts, .. } = result;
        if rows.is_empty() {
            continue;
        }
        total_kmer_count += kmer_count;
        for (row_index, row) in rows.into_iter().enumerate() {
            rows_parts[row_index].push(row);
        }
        if build_lcs {
            lcs_parts.push(lcs);
        }
        if build_counts {
            counts_parts.push(counts);
        }
    }

    log::info!("[par_build_with_all_dummies] agglomerating results");
    let rows = Arc::new(Mutex::new(Vec::<(usize, BitVec<u64>)>::new()));
    let mut lcs = None;
    let mut counts = None;
    rayon::scope(|s| {
        let rows   = &rows;
        let lcs    = &mut lcs;
        let counts = &mut counts;

        for (row_index, row_parts) in rows_parts.into_iter().enumerate() {
            s.spawn(move |_| {
                let local_row = crate::util::parallel_bitvec_concat(row_parts);
                rows.lock().unwrap().push((row_index, local_row));
            });
        }

        if build_lcs {
            s.spawn(move |_| {
                let bit_width = usize::BITS - (k.overflowing_sub(1).0).leading_zeros();
                let bit_width = bit_width as usize;

                let bitvec_lcs_parts = lcs_parts.into_iter().map(Option::unwrap).collect::<Vec<_>>();
                let value = par_concatenate_lcs(bit_width, bitvec_lcs_parts);
                *lcs = Some(value);
            });
        }

        if build_counts {
            s.spawn(move |_| {
                let raw_data = counts_parts.into_iter().map(Option::unwrap).collect::<Vec<_>>();
                let mut total_individual_counts = 0;
                let mut total_large_counts = 0;
                for part in &raw_data {
                    total_individual_counts += part.0.len();
                    total_large_counts      += part.1.len();
                }
                let mut individual_counts = vec![0_u8; total_individual_counts];
                let mut large_counts = vec![0_u64; total_large_counts];
                let sample_distance = Counts::DEFAULT_SAMPLE_DISTANCE;
                let mut sample_information_op = None;

                rayon::scope(|s| {
                    let mut individual_counts_slice: &mut [u8]  = &mut individual_counts;
                    let mut large_counts_slice     : &mut [u64] = &mut large_counts;
                    for (part_individual_counts, part_large_counts) in &raw_data {
                        let (current_slice, residual) = individual_counts_slice.split_at_mut(
                            part_individual_counts.len()
                        );
                        individual_counts_slice = residual;
                        s.spawn(|_| current_slice.copy_from_slice(part_individual_counts));

                        let (current_slice, residual) = large_counts_slice.split_at_mut(
                            part_large_counts.len()
                        );
                        large_counts_slice = residual;
                        s.spawn(|_| current_slice.copy_from_slice(part_large_counts))
                    }
                    
                    s.spawn(|_| {
                        let mut index: usize = 0;
                        let sample_count = total_individual_counts.div_ceil(sample_distance) + 1;
                        let mut sample_information = Vec::<Sample>::with_capacity(sample_count);
                        let mut sample = Sample::default();
                        let mut large_count_index;

                        for (part_individual_counts, part_large_counts) in &raw_data {
                            large_count_index = 0;
                            for &individual_count in part_individual_counts {
                                if index % sample_distance == 0 {
                                    sample_information.push(sample);
                                }
                                sample.count += individual_count as u64;
                                if individual_count == u8::MAX {
                                    sample.count += part_large_counts[large_count_index];
                                    sample.large_counts_up_to_sample += 1;
                                    large_count_index += 1;
                                }
                                index += 1;
                            }
                        }
                        sample_information.push(sample);
                        sample_information_op = Some(sample_information);
                    });
                });

                let sample_information = sample_information_op.unwrap();
                *counts = Some(Counts {
                    individual_counts,
                    sample_distance,
                    sample_information,
                    large_counts,
                });
            });
        }
    });

    let mut rows_with_index = {
        let mutex = Arc::try_unwrap(rows).unwrap();
        mutex.into_inner().unwrap()
    };
    rows_with_index.sort_by_key(|(index, _)| *index);
    let rows = rows_with_index.into_iter().map(|(_, row)| row).collect::<Vec<_>>();

    let result = collect_output(k, rows, lcs, counts, total_kmer_count);
    log::info!("[par_build_with_all_dummies] done");
    result
}

#[derive(Debug)]
struct SbwtAllDummiesRegionResult {
    start: usize,
    kmer_count: usize,
    rows: Vec<BitVec<u64>>,
    lcs: Option<BitVec<u64>>,
    counts: Option<(Vec<u8>, Vec<u64>)>,
}

#[allow(clippy::too_many_arguments)]
fn _build_sbwt_region_with_all_dummies(
    k: usize,
    start: usize,
    end: usize,
    is_last_thread: bool,
    input: &[u8],
    suffix_array: &[usize],
    lcp: &Lcp,
    lengths: &Lcp,
    build_lcs: bool,
    build_counts: bool,
) -> SbwtAllDummiesRegionResult {
    let length = suffix_array.len();

    let length_lcp_outedge = |index: usize| -> (usize, usize, usize) {
        let current_suffix_index = suffix_array[index];
        let previous_character_index_in_input = (current_suffix_index + length - 1) % length;
        let previous_character = input[previous_character_index_in_input];
        let outedge = CHAR_TO_INDEX[previous_character as usize];
        (lengths.get(index).min(k), lcp.get(index), outedge)
    };

    let mut index = start;
    let mut start_of_k_range: bool;
    let mut start_of_k_minus_one_range: bool;

    if index != 1 {
        while index < end {
            let (length, lcp_value, _) = length_lcp_outedge(index);
            start_of_k_range           = lcp_value < length;
            start_of_k_minus_one_range = start_of_k_range && (length < k || lcp_value < k - 1);
            if start_of_k_minus_one_range {
                break;
            }
            index += 1;
        }
        if index >= end {
            return SbwtAllDummiesRegionResult {
                start,
                kmer_count: 0,
                rows: vec![],
                lcs: None,
                counts: None,
            }
        }
    }

    let capacity = end - start;
    let mut rows = Vec::<BitVec<u64>>::new();
    for _ in 0..4 {
        rows.push(BitVec::with_capacity(capacity));
    }

    let bit_width = usize::BITS - (k.overflowing_sub(1).0).leading_zeros();
    let bit_width = bit_width as usize;
    let mut lcs_count: usize = 0;
    let mut lcs: Option<BitVec<u64>> = if build_lcs {
        let value = bitvec::bitvec![u64, bitvec::order::Lsb0; 0; bit_width * capacity];
        Some(value)
    } else {
        None
    };

    let mut counts = if build_counts {
        let individual_counts = Vec::<u8>::with_capacity(capacity);
        let large_counts = Vec::<u64>::with_capacity(capacity);
        Some((individual_counts, large_counts))
    } else {
        None
    };

    let mut current_set: u8 = 0;

    let mut kmer_count: usize = 0;

    let mut k_range_count    : u64 = 0;
    let mut individual_count : u64 = 0;

    if start == 1 {
        lcs_count += 1;
        k_range_count = 1;
        if build_counts {
            counts.as_mut().unwrap().0.push(0);
        }
    }

    while index < length {
        let (length, lcp_value, outedge) = length_lcp_outedge(index);
        start_of_k_range           = lcp_value < length;
        start_of_k_minus_one_range = start_of_k_range && (length < k || lcp_value < k - 1);

        if start_of_k_minus_one_range {
            while k_range_count > 0 {
                push_set(&mut rows, current_set);
                current_set = 0;
                k_range_count -= 1;
            }

            current_set = 0;

            if index >= end {
                break;
            }
        }

        if start_of_k_range {
            k_range_count += 1;
            if length >= k {
                kmer_count += 1;
            }

            if build_lcs {
                let (_, suffix) = lcs.as_mut().unwrap().split_at_mut(lcs_count * bit_width);
                let (word, _)   = suffix.split_at_mut(bit_width);
                word.store_le(lcp_value as u64);
                lcs_count += 1;
            }

            if build_counts && individual_count > 0 {
                let counts = counts.as_mut().unwrap();
                if individual_count >= u8::MAX as u64 {
                    counts.0.push(u8::MAX);
                    counts.1.push(individual_count - u8::MAX as u64);
                } else {
                    counts.0.push(individual_count as u8);
                }
                individual_count = 0;
            }
        }

        if build_counts && length > 0 {
            individual_count += 1;
        }

        current_set |= (1 << outedge) as u8;

        index += 1;
    }

    if is_last_thread {
        while k_range_count > 0 {
            push_set(&mut rows, current_set);
            current_set = 0;
            k_range_count -= 1;
        }
    }

    if build_lcs {
        lcs.as_mut().unwrap().truncate(bit_width * lcs_count);
    }

    if build_counts && individual_count > 0 {
        let counts = counts.as_mut().unwrap();
        if individual_count >= u8::MAX as u64 {
            counts.0.push(u8::MAX);
            counts.1.push(individual_count - u8::MAX as u64);
        } else {
            counts.0.push(individual_count as u8);
        }
    }

    SbwtAllDummiesRegionResult {
        start,
        kmer_count,
        rows,
        lcs,
        counts
    }
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

    fn make_suffix_array_full_context(concatenation: &[u8]) -> Vec<usize> {
        let mut suffix_array = (0..concatenation.len()).collect::<Vec<_>>();
        suffix_array.sort_by_key(|index| &concatenation[*index..]);
        suffix_array
    }

    fn make_suffix_array(threads: usize, concatenation: &mut Vec<u8>, context: usize) -> Vec<usize> {
        let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(threads).build().unwrap();
        let mut result_op = None;
        thread_pool.scope(|s| {
            s.spawn(|_| {
                let result = par_bounded_context_suffix_array_bucket_sort(concatenation, context, 4);
                result_op = Some(result);
            });
        });
        result_op.unwrap()
    }

    #[inline]
    fn find_length_and_lcp_value(
        k: usize,
        input: &[u8],
        mut word_count: usize,
        mut current_suffix_index: usize,
        mut previous_suffix_index: usize
    ) -> (usize, usize) {
        if word_count * LANES == k {
            word_count += 1;
        }

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
                stop_accumulating_length |= length >= k || bitmask != 0;
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


    fn build_lcp_and_lengths(threads: usize, k: usize, seqs: Vec<Vec<u8>>, bounded: bool) {
        let mut concatenation = make_concatenation(&seqs);
        let suffix_array = if bounded {
            make_suffix_array(threads, &mut concatenation, k)
        } else {
           make_suffix_array_full_context(&concatenation)
        };
        let length = suffix_array.len();
        pad_input(&mut concatenation, k, length);

        let mut lcp_op = None;
        let mut lengths_op = None;

        let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(threads).build().unwrap();
        thread_pool.scope(|s| {
            s.spawn(|_| {
                lcp_op = Some(par_build_lcp(threads, &concatenation, &suffix_array, k, true));
            });
        });
        thread_pool.scope(|s| {
            s.spawn(|_| {
                lengths_op = Some(par_build_lengths(threads, &concatenation, &suffix_array, k));
            });
        });

        let lcp = lcp_op.unwrap();
        let lengths = lengths_op.unwrap();

        let mut ok = true;
        let word_count = word_count(k);
        for i in 1..length {
            let current_suffix = suffix_array[i];
            let previous_suffix = suffix_array[i - 1];
            let (length, lcp_value) = find_length_and_lcp_value(
                k,
                &concatenation,
                word_count,
                current_suffix,
                previous_suffix
            );

            let true_length = length.min(k + 1);
            let true_lcp_value = lcp_value.min(k);
            let built_length = lengths.get(i);
            let built_lcp_value = lcp.get(i);

            let end = (current_suffix + k).min(concatenation.len());
            let s = str::from_utf8(&concatenation[current_suffix..end]).unwrap();

            println!(
                "[{:4}] [len:{:4}|{:4}] [lcp:{:4}|{:4}] {}",
                current_suffix,
                true_length,
                built_length,
                true_lcp_value,
                built_lcp_value,
                s,
            );

            ok &= true_length == built_length;
            ok &= true_lcp_value == built_lcp_value;
        }
        assert!(ok);
    }

    #[test]
    fn build_lcp_and_lengths_01() {
        let threads = 3;
        let k = 17;

        let seqs = seqs![
            b"ACGACGACCACCGACACACAAACCCAAACGTGAACGTTAA",
            b"ACCCAAAAGTGTGTGAGAGTGTGAGCAGTGCATGATGCAA",
            b"GTGAGAGAGTGATGGACCAAAAAAAAAAAAAAAACCCGTA",
            b"GAAAAAAAAAAAAAACCAMTGAGGAGAGAGAGGGGTTTTT",
            b"ACMACMAGAMCAMMAGMAMCAMTGMAMMGATATGAGACMM",
            b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
            b"MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        ];

        build_lcp_and_lengths(threads, k, seqs.clone(), true);
        build_lcp_and_lengths(threads, k, seqs, false);
    }

    #[test]
    fn build_lcp_and_lengths_02() {
        let threads = 3;
        let k = 3;
        let seqs = seqs![
            b"ATCGTTCGTTCGTTCG",
            b"ATTTTTTTTTTTTCGTTTTTTTTTTTTTCGTTTTTTTTTTTTTCGTTTTTTTTTTTTTCG",
        ];

        build_lcp_and_lengths(threads, k, seqs.clone(), true);
        build_lcp_and_lengths(threads, k, seqs, false);
    }

    fn randomised_kmers(bounded: bool) {
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

        let mut concatenation = make_concatenation(&seqs);
        let suffix_array = if bounded {
            make_suffix_array(4, &mut concatenation, k)
        } else {
           make_suffix_array_full_context(&concatenation)
        };

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

            let correct_counts = vodbg.counts().unwrap();
            let constructed_counts = &counts.unwrap();
            assert_eq!(correct_counts, constructed_counts);

            // With full context suffix array.
            let Output {
                sbwt: mut constructed_sbwt,
                lcs: constructed_lcs,
                counts
            } = build::<SubsetMatrix>(4, concatenation, suffix_array, k, true, true, true);

            constructed_sbwt.build_select();
            let constructed_lcs = constructed_lcs.unwrap();

            assert_eq!(correct_sbwt, constructed_sbwt);
            assert_eq!(constructed_lcs, correct_lcs);
            let constructed_counts = &counts.unwrap();
            assert_eq!(correct_counts, constructed_counts);
        }
    }

    #[test]
    fn randomised_kmers_bounded_context() {
        randomised_kmers(true);
    }

    #[test]
    fn randomised_kmers_full_context() {
        randomised_kmers(false);
    }
}

// Code by Martin Kostadinov.

pub mod preprocessing;
pub mod input_structures;

use crate::vodbg::count::{Counts, Sample};
use crate::{LcsArray, SbwtIndex, SubsetSeq};
use input_structures::{Bwt, Lcp};

use std::io::Read;

use bitvec::vec::BitVec;
use simple_sds_sbwt::bit_vector::BitVector;
use simple_sds_sbwt::int_vector::IntVector;
use simple_sds_sbwt::ops::{BitVec as BitVecTrait, Push, Rank, Select};
use simple_sds_sbwt::raw_vector::{AccessRaw, RawVector};

pub struct Output<SS: SubsetSeq + Send> {
    pub sbwt: SbwtIndex<SS>,
    pub lcs: Option<LcsArray>,
    pub counts: Option<Counts>,
}

type Row = BitVec<u64>;

/// Constructs an SbwtIndex (and optionally an LCS array and the counts of the k-mers in the input
/// sequences) from the BWT and a Longest Commond Prefix array of the Suffix Array of the
/// concatenated input sequences.
///
/// Expects the BWT and the LCP to correspond to a concatenation as done by the
/// [preprocessing::concatenate_sequences] function, that is, the start of the concatenation is the
/// '#' character and then before each sequence there is a single '$' character as well as at the
/// end of the concatenated string. The input sequences should be reversed as the k-mers in the
/// SBWT are sorted colexicographically. The BWT is expected to be in ascii whereas the LCP is
/// expected to consist of 64 bit unsigned integers. The number of characters in the BWT and the
/// number of integers in the LCP should be the same.
///
/// The `build_counts` parameter is taken into account only if all dummies are included in the
/// SBWT, otherwise the Counts data structure will not answer queries about shorter k-mers
/// correctly.
pub fn build_from_input<RBwt, RLcp, SS>(
    bwt_input: &mut RBwt,
    lcp_input: &mut RLcp,
    length: usize,
    k: usize,
    build_lcs: bool,
    add_all_dummies: bool,
    build_counts: bool,
) -> std::io::Result<Output<SS>>
where
    RBwt: Read,
    RLcp: Read,
    SS: SubsetSeq + Send
{
    log::info!(
        "[build_from_input] begin [length: {} | build_lcs: {} | add_all_dummies: {} | build_counts {}]",
        length, build_lcs, add_all_dummies, build_counts
    );
    let bwt = preprocessing::ascii_to_bwt(bwt_input, length)?;
    let lcp = preprocessing::truncate_lcp::<_, true>(lcp_input, length, k)?;
    let result = if !add_all_dummies {
        build_without_redundant_dummies(&bwt, &lcp, k, build_lcs)
    } else {
        build_with_all_dummies(&bwt, &lcp, k, build_lcs, build_counts)
    };
    log::info!("[build_from_input] done");
    Ok(result)
}

pub(crate) const FULL_SET: u8 = 0b00011110;

/// The subroutine of the [build_from_input] function which handles the construction of the
/// SbwtIndex without the redundant dummies.
pub fn build_without_redundant_dummies<SS: SubsetSeq + Send>(
    bwt: &Bwt,
    lcp: &Lcp,
    k: usize,
    build_lcs: bool,
) -> Output<SS> {
    log::info!("[build_without_redundant_dummies] begin");
    let aux = build_full_auxiliary_bitvectors(bwt, lcp, k);
    let dummy_marks = build_dummy_marks(bwt, k, &aux);

    let FullAuxiliaryBitVectors {
        kmer_count,
        shorter_than_k,
        equal_to_k,
        k_minus_one_ranges,
        k_ranges
    } = aux;
    drop(equal_to_k);

    let separator_count = bwt.counts[1];
    let (rows, lcs) = _build_without_redundant_dummies(
        k,
        separator_count,
        bwt.len(),
        lcp,
        &k_minus_one_ranges,
        &k_ranges,
        &shorter_than_k,
        &dummy_marks.keep_dummy,
        build_lcs,
        |current_set: &mut u8, index| {
            if dummy_marks.keep_dummy.bit(index) || dummy_marks.keep_outedge.bit(index) {
                *current_set = include_letter(bwt, index, *current_set);
            }
        },
        |current_set: &mut u8, index| {
            *current_set = include_letter(bwt, index, *current_set);
        },
    );
    let result = collect_output(k, rows, lcs, None, kmer_count);
    log::info!("[build_without_redundant_dummies] done");
    result
}

#[allow(clippy::too_many_arguments)]
#[inline]
pub(crate) fn _build_without_redundant_dummies<FDummy, FNonDummy>(
    k: usize,
    separator_count: usize,
    length: usize,
    lcp: &Lcp,
    k_minus_one_ranges: &BitVector,
    k_ranges: &RawVector,
    shorter_than_k: &BitVector,
    keep_dummy: &RawVector,
    build_lcs: bool,
    mut add_dummy_outedge: FDummy,
    mut add_non_dummy_outedge: FNonDummy,
) -> (Vec<Row>, Option<IntVector>)
where
    FDummy: FnMut(&mut u8, usize),
    FNonDummy: FnMut(&mut u8, usize)
{
    log::info!("[_build_without_redundant_dummies] begin");
    let mut rows = Vec::<BitVec<u64>>::new();
    for _ in 0..4 {
        // note(mk): These overestimating allocations should reserve the pages in the virtual
        // memory space of the process, but it shouldn't actually use all of them. Potential
        // breaking point!
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

    let mut current_set: u8 = 0;

    for index in 1..separator_count {
        add_dummy_outedge(&mut current_set, index);
        if current_set & FULL_SET == FULL_SET {
            break;
        }
    }
    push_set(&mut rows, current_set);

    let mut current_set = 0;
    let mut current_lcs_value  = k - 1;
    let mut include_dummy_kmer = false;
    let mut has_dummy_kmer     = false;
    let mut k_range_count = 0;

    for index in separator_count..length {
        if k_minus_one_ranges.get(index) {
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
            if keep_dummy.bit(index) {
                if build_lcs && !include_dummy_kmer {
                    lcs.as_mut().unwrap().push(current_lcs_value as u64);
                    current_lcs_value = k - 1;
                }
                include_dummy_kmer = true;
            }
            add_dummy_outedge(&mut current_set, index);
        } else {
            if build_lcs && is_start_of_k_range {
                lcs.as_mut().unwrap().push(current_lcs_value as u64);
                current_lcs_value = k - 1;
            }
            add_non_dummy_outedge(&mut current_set, index);
        }
    }

    if has_dummy_kmer && !include_dummy_kmer {
        k_range_count -= 1;
    }
    while k_range_count > 0 {
        push_set(&mut rows, current_set);
        current_set = 0;
        k_range_count -= 1;
    }

    log::info!("[_build_without_redundant_dummies] begin");
    (rows, lcs)
}

pub fn build_with_all_dummies<SS: SubsetSeq + Send>(
    bwt: &Bwt,
    lcp: &Lcp,
    k: usize,
    build_lcs: bool,
    build_counts: bool,
) -> Output<SS> {
    log::info!("[build_with_all_dummies] begin");
    let aux = build_parital_auxiliary_bitvectors(bwt, lcp, k);
    let set_count = aux.set_count;

    let mut rows = Vec::<Row>::new();
    for _ in 0..4 {
        rows.push(BitVec::with_capacity(set_count));
    }

    let mut lcs: Option<IntVector> = if build_lcs {
        let bit_width = usize::BITS - (k.overflowing_sub(1).0).leading_zeros();
        let bit_width = bit_width as usize;
        let mut value = IntVector::with_capacity(set_count, bit_width).unwrap();
        value.push(0); // '$...$' dummy k-mer
        Some(value)
    } else {
        None
    };

    let mut counts: Option<Counts> = if build_counts {
        let sample_capacity = set_count / Counts::DEFAULT_SAMPLE_DISTANCE;
        let mut value = Counts {
            individual_counts: Vec::with_capacity(set_count),
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
    let separator_count = bwt.counts[1];

    for index in 1..separator_count {
        current_set = include_letter(bwt, index, current_set);
        if current_set & FULL_SET == FULL_SET {
            break;
        }
    }
    push_set(&mut rows, current_set);
    current_set = 0;

    let mut k_range_count = 0;
    let mut sbwt_index = 0;
    let mut count: u64 = 0;
    let mut individual_count: u64 = 0;
    let mut large_counts_up_to_sample: usize = 0;
    for index in separator_count..bwt.len() {
        if aux.k_minus_one_ranges.bit(index) {
            while k_range_count > 0 {
                push_set(&mut rows, current_set);
                current_set = 0;
                k_range_count -= 1;
            }

            current_set = 0;
        }

        if aux.k_ranges.bit(index) {
            k_range_count += 1;

            if build_lcs {
                lcs.as_mut().unwrap().push(lcp.get(index) as u64);
            }

            if build_counts {
                sbwt_index += 1;
                let counts = counts.as_mut().unwrap();
                if sbwt_index % Counts::DEFAULT_SAMPLE_DISTANCE == 0 {
                    counts.sample_information.push(Sample {
                        count,
                        large_counts_up_to_sample ,
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
        }

        if build_counts {
            count += 1;
            individual_count += 1;
            if individual_count == u8::MAX as u64 {
                large_counts_up_to_sample += 1;
            }
        }

        current_set = include_letter(bwt, index, current_set);
    }

    while k_range_count > 0 {
        push_set(&mut rows, current_set);
        current_set = 0;
        k_range_count -= 1;
    }

    if build_counts {
        let counts = counts.as_mut().unwrap();
        // Always create a sample at the end of the array.
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

    let result = collect_output(k, rows, lcs, counts, aux.kmer_count);
    log::info!("[build_with_all_dummies] done");
    result
}

pub(crate) fn collect_output<SS: SubsetSeq + Send>(
    k: usize,
    rows: Vec<Row>,
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

    log::info!("[collect_output] constructing done");
    Output {
        sbwt: index,
        lcs,
        counts,
    }
}


/// Auxiliary BitVectors needed for the calculation of the SbwtIndex without redundant dummies.
struct FullAuxiliaryBitVectors {
    /// It is easier to calculate the true k-mers while creating the auxiliary bitvectors rather
    /// than trying to figure out their count later.
    kmer_count: usize,
    /// Used to figure out whether a k-mer at the beginning of a sequence has a true k-mer as a
    /// predecessor in order to figure out whether the dummy k-mer is necessary. In the final pass
    /// it is used to figure out if a given (k-1)-range contains a region of dummy k-mers.
    shorter_than_k: BitVector,
    /// Used in the pass which marks the dummy k-mers that need to be kept in order to identify the
    /// k-mers which are at the beginning of an input sequence.
    equal_to_k: RawVector,
    /// Used in order to figure out the bounds at which to check the [Self::shorter_than_k]
    /// bitvector whether a non-dummy k-mer at the beginning of an input sequence has a non-dummy
    /// k-mer as a predecessor. In addition, it is used to figure out the bounds at which to
    /// "collect" outedges for the first k-mer with a given (k-1) suffix in the SBWT.
    k_minus_one_ranges: BitVector,
    /// Used to enumerate the sets of the SBWT.
    k_ranges: RawVector,
}

fn build_full_auxiliary_bitvectors(bwt: &Bwt, lcp: &Lcp, k: usize) -> FullAuxiliaryBitVectors {
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

    log::info!("[build_full_auxiliary_bitvectors] begin");
    let len = bwt.len();
    let mut shorter_than_k     = RawVector::with_len(len, false);
    let mut equal_to_k         = RawVector::with_len(len, false);
    let mut k_minus_one_ranges = RawVector::with_len(len, false);
    let mut k_ranges           = RawVector::with_len(len, false);
    let mut kmer_count = 0;
    let mut order = 1;
    let mut current_length = 0;
    for _ in 1..len {
        let (next_order, character) = bwt.lf_step(order);
        order = next_order;
        if character == b'$' {
            current_length = 0;
            shorter_than_k.set_bit(order, true);
        } else {
            current_length += 1;
            if current_length < k {
                shorter_than_k.set_bit(order, true);
            } else if current_length == k {
                equal_to_k.set_bit(order, true);
            } else {
                current_length = k;
            }
        }
        let lcp_value = lcp.get(order);
        if lcp_value < current_length {
            k_ranges.set_bit(order, true);
            if current_length >= k {
                kmer_count += 1;
            }
            if current_length < k || lcp_value < k - 1 {
                k_minus_one_ranges.set_bit(order, true);
            }
        }
    }

    // Skip the '#' region at the beginning.
    k_minus_one_ranges.set_bit(1, true);
    k_ranges.set_bit(1, true);

    log::info!("[build_full_auxiliary_bitvectors] rank for shorter than k k-mers bitvector");
    let mut shorter_than_k = BitVector::from(shorter_than_k);
    shorter_than_k.enable_rank();

    log::info!("[build_full_auxiliary_bitvectors] rank and select for (k-1)-ranges bitvector");
    let mut k_minus_one_ranges = BitVector::from(k_minus_one_ranges);
    k_minus_one_ranges.enable_rank();
    k_minus_one_ranges.enable_select();

    log::info!("[build_full_auxiliary_bitvectors] done");
    FullAuxiliaryBitVectors {
        kmer_count,
        shorter_than_k,
        equal_to_k,
        k_minus_one_ranges,
        k_ranges
    }
}

struct PartialAuxiliaryBitVectors {
    set_count: usize,
    kmer_count: usize,
    k_minus_one_ranges: RawVector,
    k_ranges: RawVector,
}

fn build_parital_auxiliary_bitvectors(bwt: &Bwt, lcp: &Lcp, k: usize) -> PartialAuxiliaryBitVectors {
    log::info!("[build_partial_auxiliary_bitvectors] begin");
    let len = bwt.len();
    let mut k_minus_one_ranges = RawVector::with_len(len, false);
    let mut k_ranges           = RawVector::with_len(len, false);
    let mut set_count = 0;
    let mut kmer_count = 0;
    let mut order = 1;
    let mut current_length = 0;
    for _ in 1..len {
        let (next_order, character) = bwt.lf_step(order);
        order = next_order;
        if character == b'$' {
            current_length = 0;
        } else {
            current_length += 1;
            if current_length > k {
                current_length = k;
            }
        }
        let lcp_value = lcp.get(order);
        if lcp_value < current_length {
            k_ranges.set_bit(order, true);
            set_count += 1;
            if current_length >= k {
                kmer_count += 1;
            }
            if current_length < k || lcp_value < k - 1 {
                k_minus_one_ranges.set_bit(order, true);
            }
        }
    }

    // Skip the '#' region at the beginning.
    k_minus_one_ranges.set_bit(1, true);
    k_ranges.set_bit(1, true);

    log::info!("[build_partial_auxiliary_bitvectors] done");
    PartialAuxiliaryBitVectors {
        set_count,
        kmer_count,
        k_minus_one_ranges,
        k_ranges,
    }
}

struct DummyMarks {
    /// We keep a dummy (k-1)-mer and all of its predecessors if a k-mer at the beginning of an
    /// input sequence would not have a non-dummy k-mer as a predecessor in the SBWT graph.
    keep_dummy: RawVector,
    /// Marks (k-1)-mers in order to keep their outedge as it is not guaranteed that the non-dummy
    /// k-mer that exists as a predecessor to a non-dummy k-mer at the beginning of a sequence has
    /// the needed label. This is needed as we don't want to include all letters of dummy k-mers
    /// whose "true" lenght is less than k-1 and there is no way to figure out which ones are
    /// those from the [FullAuxiliaryBitVectors::shorter_than_k] bitvector alone.
    keep_outedge: RawVector,
}

fn build_dummy_marks(bwt: &Bwt, k: usize, aux: &FullAuxiliaryBitVectors) -> DummyMarks {
    log::info!("[build_dummy_marks] begin");
    let len = bwt.len();
    let mut keep_dummy   = RawVector::with_len(len, false);
    let mut keep_outedge = RawVector::with_len(len, false);

    let start = bwt.counts[1];
    let mut predecessor_confirmed = false;
    for index in start..len {
        if aux.k_ranges.bit(index) {
            predecessor_confirmed = false;
        }

        if aux.equal_to_k.bit(index) {
            let (predecessor, _) = bwt.inverse_lf_step(index);
            // If we haven't found a full k-mer as a predecessor for this k-range, search for it.
            if !predecessor_confirmed {
                predecessor_confirmed |= has_full_kmer_predecessor(
                    predecessor, bwt, &aux.k_minus_one_ranges, &aux.shorter_than_k
                );
            }

            if predecessor_confirmed {
                keep_outedge.set_bit(predecessor, true);
            } else {
                keep_predecessors(predecessor, bwt, k, &mut keep_dummy);
                predecessor_confirmed = true;
            }
        } else if !aux.shorter_than_k.get(index) {
            // If the true length of the prefix of this suffix is not equal to k and it is not
            // shorter than k, this means that it is longer than k. If this is the case, this means
            // that this k-range has a predecessor.
            predecessor_confirmed = true;
        }
    }

    log::info!("[build_dummy_marks] done");
    DummyMarks {
        keep_dummy,
        keep_outedge,
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

fn keep_predecessors(mut predecessor: usize, bwt: &Bwt, mut k: usize, keep_dummy: &mut RawVector) {
    while k > 0 {
        keep_dummy.set_bit(predecessor, true);
        let (next, _) = bwt.inverse_lf_step(predecessor);
        predecessor = next;
        k -= 1;
    }
}

#[inline]
pub(crate) fn push_set(rows: &mut [Row], set: u8) {
    for i in 1..=4 {
        rows[i - 1].push(set & (1 << i) != 0);
    }
}

#[inline]
fn include_letter(bwt: &Bwt, index: usize, current_set: u8) -> u8 {
    (1 << bwt.get_char_index(index)) as u8 | current_set
}

#[cfg(test)]
pub(crate) use tests::make_concatenation;

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

    fn construct_bwt_lcp(concatenation: &[u8], k: usize) -> (Bwt, Lcp) {
        let length = concatenation.len();

        let mut suffix_array = (0..concatenation.len()).collect::<Vec<_>>();
        suffix_array.sort_by_key(|index| &concatenation[*index..]);
        let bwt = suffix_array.iter()
            .map(|index| concatenation[(index + concatenation.len() - 1) % concatenation.len()])
            .collect::<Vec<_>>();
        let bwt = preprocessing::ascii_to_bwt(&mut bwt.as_slice(), length).unwrap();

        let mut lcp = Vec::<u64>::with_capacity(length);
        lcp.push(0);
        for i in 1..bwt.len() {
            let mut it_a = suffix_array[i - 1];
            let mut it_b = suffix_array[i];
            let mut lcp_value = 0;
            while it_a < suffix_array.len() && it_b < suffix_array.len() {
                if concatenation[it_a] != concatenation[it_b] {
                    break;
                }
                lcp_value += 1;
                it_a += 1;
                it_b += 1;
            }
            lcp.push(lcp_value);
        }

        let lcp_bytes: Vec<u8> = lcp.into_iter().map(|value| value.to_le_bytes()).collect::<Vec<_>>().concat();
        let lcp = preprocessing::truncate_lcp::<_, true>(&mut lcp_bytes.as_slice(), bwt.len(), k).unwrap();

        (bwt, lcp)
    }

    #[test]
    fn bwt_lcp() {
        let seqs = seqs![b"CCACC", b"AGG", b"TACG", b"ACG", b"NOP"];
        let concatenation = make_concatenation(&seqs);
        let (bwt, mut lcp) = construct_bwt_lcp(&concatenation, 3);

        let mut order = 0;
        let mut index = bwt.len() - 1;
        for _ in 0..bwt.len() {
            let (new_order, character) = bwt.lf_step(order);
            if character != b'#' {
                // The inverse step works only for the letters in the alphabet and $, not for #.
                let (found_order, _) = bwt.inverse_lf_step(new_order);
                assert_eq!(order, found_order, "{}", index);
            }
            assert_eq!(concatenation[index], character, "{}", index);
            if index == 0 {
                index = bwt.len() - 1;
            } else {
                index -= 1;
            }
            order = new_order;
        }

        let correct = &[0, 0, 1, 2, 3, 3, 1, 1, 3, 2, 0, 2, 1, 1, 0, 1, 2, 2, 1, 2, 0, 1, 3, 1, 0];
        let values = (&mut lcp).collect::<Vec<_>>();
        assert_eq!(correct, values.as_slice());
    }

    #[test]
    fn randomised_kmers() {
        use rand_chacha::ChaCha20Rng;
        use rand_chacha::rand_core::SeedableRng;
        use rand_chacha::rand_core::RngCore;

        let k: usize = 16;
        let kmer_length: usize = 32;
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
        let (bwt, lcp) = construct_bwt_lcp(&concatenation, k);

        {
            // Without redundant dummies.
            let (correct_sbwt, correct_lcs) = BitPackedKmerSortingMem::new_from_vecs(&seqs, k)
                .build_lcs(true)
                .run();

            let Output {
                sbwt: constructed_sbwt,
                lcs: constructed_lcs,
                counts
            } = build_without_redundant_dummies::<SubsetMatrix>(&bwt, &lcp, k, true);

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
            let (mut correct_sbwt, correct_lcs) = BitPackedKmerSortingMem::new_from_vecs(&seqs, k)
                .build_lcs(true)
                .add_all_dummy_paths(true)
                .run();

            let Output {
                sbwt: constructed_sbwt,
                lcs: constructed_lcs,
                counts
            } = build_with_all_dummies::<SubsetMatrix>(&bwt, &lcp, k, true, true);

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


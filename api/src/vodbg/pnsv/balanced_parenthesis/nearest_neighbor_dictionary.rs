// Code by Martin Kostadinov adapted to Rust from the C++ library sdsl-lite.

use simple_sds_sbwt::{
    bit_vector::BitVector,
    int_vector::IntVector,
    ops::{Access, Rank},
    raw_vector::{AccessRaw, RawVector}
};

///
/// Support for rank, select, previous and next for a uniformly sparse bitvector (with respect to
/// the bits set to 1) from the paper "A simple optimal representation for balanced parentheses" by
/// Richard F. Geary, Naila Rahman, Rajeev Raman, Venkatesh Raman.
///
/// Implementation adapted from the C++ library sdsl-lite to Rust.
///
// note(mk): Should the structure be generic over the SAMPLE_OFFSET? If it is, the compiler will be
// able to optimise the division and modulo operations. If it isn't, it will be more difficult to
// dynamically create the structure.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct NearestNeighbourDictionary<const SAMPLE_OFFSET: usize = 8> {
    /// Given the index of the bits set to 1 in the bitvector - i, store the position of those bits
    /// in the bitvector whose index i is equal to 0 modulo SAMPLE_OFFSET (for 1 based indexing).
    /// Those bits will be called sampled bits. Corresponds to array A1 in the paper.
    pub sample_bit_indices: IntVector,
    
    /// Stores the distance between each bit set to 1 in the bitvector. Corresponds to array A2
    /// from the paper. In practice, it skips storing the distance between the sampled bit and its
    /// previous bit as those are not used.
    pub differences: IntVector,

    /// For each block of size SAMPLE_OFFSET record whether it contains a sampled bit. There is a
    /// minimum of (SAMPLE_OFFSET - 1) bits between each sampled bit therefore a block can contain
    /// at most 1 sampled bit. Will be guaranteed to support rank operations. Corresponds to the
    /// combined arrays A3 and A4 in the paper.
    pub block_contains_sample_bit: BitVector,

    /// The number of bits set to 1 in the bitvector.
    pub one_count: usize,
    /// The length of the original bitvector.
    pub bitvector_length: usize,
}

impl<const SAMPLE_OFFSET: usize> NearestNeighbourDictionary<SAMPLE_OFFSET> {
    pub fn empty() -> Self {
        Self {
            sample_bit_indices: IntVector::new(1).unwrap(),
            differences: IntVector::new(1).unwrap(),
            block_contains_sample_bit: RawVector::new().into(),
            one_count: 0,
            bitvector_length: 0,
        }
    }

    /// Expects an iterator over the indices of the bits set to 1 in the bitvector.
    pub fn new<B: Iterator<Item = usize> + Clone>(one_bit_indices: B, bitvector_length: usize) -> Self
    {
        let mut one_count = 0;
        let mut previous_one_index = 0;

        let mut max_distance_between_two_ones = 1;
        let mut distance;

        for i in one_bit_indices.clone() {
            one_count += 1;

            distance = if one_count == 1 {
                // Assume there is an imaginary 1 immediately before the beginning of the
                // bitvector.
                i+1
            } else {
                i - previous_one_index
            };

            if distance > max_distance_between_two_ones {
                max_distance_between_two_ones = distance;
            }

            previous_one_index = i;
        }

        // Store an additional element for the imaginary 1 immediately before the beginning of the
        // bitvector.
        let sample_bit_count = 1 + one_count / SAMPLE_OFFSET;
        let sample_bit_position_bit_width = crate::vodbg::util::bit_width(bitvector_length);
        let mut sample_bit_indices = IntVector::with_len(sample_bit_count, sample_bit_position_bit_width, 0)
            .expect("The sample_bit_position_bit_width should not be greater than 64.");

        // Skip storing the distance between the sampled bits and the previous 1 bit.
        let differences_count = one_count - one_count / SAMPLE_OFFSET;
        let differences_bit_width = crate::vodbg::util::bit_width(max_distance_between_two_ones);
        let mut differences = IntVector::with_len(differences_count, differences_bit_width, 0)
            .expect("The differences_bit_width should not be greater than 64.");

        // Ensure that the last block that may be incomplete is included in the count.
        #[allow(clippy::manual_div_ceil)]
        let block_count = (bitvector_length + SAMPLE_OFFSET - 1) / SAMPLE_OFFSET;
        let mut block_contains_sample_bit_raw = RawVector::with_len(block_count, false);

        one_count = 0;
        previous_one_index = 0;
        for i in one_bit_indices {
            one_count += 1;

            if one_count % SAMPLE_OFFSET == 0 {
                let sample_bit_index_in_array = one_count / SAMPLE_OFFSET;
                sample_bit_indices.set(sample_bit_index_in_array, i as u64);

                let block_index = i / SAMPLE_OFFSET;
                block_contains_sample_bit_raw.set_bit(block_index, true);
            } else {
                let sample_bits_before_current_bit = one_count / SAMPLE_OFFSET;
                let difference_index = (one_count - 1) - sample_bits_before_current_bit;
                let difference = i - previous_one_index;
                differences.set(difference_index, difference as u64);
            }

            previous_one_index = i;
        }

        let mut block_contains_sample_bit = BitVector::from(block_contains_sample_bit_raw);
        block_contains_sample_bit.enable_rank();

        Self {
            sample_bit_indices,
            differences,
            block_contains_sample_bit,
            one_count,
            bitvector_length,
        }
    }
    
    /// Finds the number of 1s in the bitvector before the specified index. The index supplied to
    /// this method should be 0-based.
    pub fn rank(&self, index_in_bitvector: usize) -> usize {
        // The value i' from the paper.
        let block_index = index_in_bitvector / SAMPLE_OFFSET;
        // The value r from the paper.
        let previous_sample_bit_number = self.block_contains_sample_bit.rank(block_index);

        // The number of 1 bits before the current one.
        let mut result_count = previous_sample_bit_number * SAMPLE_OFFSET;
        if result_count >= self.one_count {
            return result_count;
        }

        let mut current_bit_index = self.sample_bit_indices.get(previous_sample_bit_number) as usize;

        loop {
            // Store the next result and check whether the position of the bit corresponding to it
            // fits in the range.
            let next_result_count = result_count + 1;
            
            // The sample_bit_indices array is 0-based, however, there is one imaginary 1-bit and
            // its index in the array is 0.
            #[allow(clippy::manual_is_multiple_of)]
            if next_result_count % SAMPLE_OFFSET == 0 {
                let sample_bit_number = next_result_count / SAMPLE_OFFSET;
                current_bit_index = self.sample_bit_indices.get(sample_bit_number) as usize;
            } else {
                let sample_bits_before_next_result = next_result_count / SAMPLE_OFFSET;
                // The index of the difference corresponding to the next_result bit is 0-based. The
                // next_result variable stores the number of ones i.e. result will be the
                // corresponding index in the differences array.
                let difference_index = result_count - sample_bits_before_next_result;
                current_bit_index += self.differences.get(difference_index) as usize;
            }

            if current_bit_index >= index_in_bitvector {
                break;
            }

            result_count = next_result_count;

            if result_count >= self.one_count {
                break;
            }
        }

        result_count
    }

    /// Finds the position of the 1 bit specified by the index. The index supplied to
    /// this method should be 0-based.
    pub fn select(&self, index_of_one: usize) -> usize {
        if index_of_one >= self.one_count {
            return self.bitvector_length;
        }
        let sample_bit_number = (index_of_one + 1) / SAMPLE_OFFSET;
        let mut result_position = self.sample_bit_indices.get(sample_bit_number) as usize;

        // Subtract the numer of sample bits in order to get the beginning diference index.
        let difference_index_begin = sample_bit_number * (SAMPLE_OFFSET - 1);
        let steps_to_bit_to_select = (index_of_one + 1) % SAMPLE_OFFSET;
        let difference_index_end = difference_index_begin + steps_to_bit_to_select;

        for difference_index in difference_index_begin..difference_index_end {
            result_position += self.differences.get(difference_index) as usize;
        }

        result_position
    }

    pub fn previous(&self, index_in_bitvector: usize) -> usize {
        let rank = self.rank(index_in_bitvector);
        if rank == 0 {
            return 0;
        }
        self.select(rank - 1)
    }

    pub fn next(&self, index_in_bitvector: usize) -> usize {
        // The number of 1-bits up to and including the current index will give the 0-based index
        // of the next 1-bit.
        let rank = self.rank(index_in_bitvector + 1);
        self.select(rank)
    }

    pub const fn sample_offset(&self) -> usize {
        SAMPLE_OFFSET
    }

    pub fn len(&self) -> usize {
        self.bitvector_length
    }

    pub fn is_empty(&self) -> bool {
        // Shut up clippy.
        self.len() == 0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use simple_sds_sbwt::ops::BitVec;

    const SAMPLE_OFFSET: usize = 8;
    fn make_nnd(length: usize, step: usize) -> NearestNeighbourDictionary<SAMPLE_OFFSET> {
        let mut raw_vector = RawVector::with_len(length, false);
        let mut i = 0;
        while i < raw_vector.len() {
            raw_vector.set_bit(i, true);
            i += step;
        }

        let bitvector: BitVector = raw_vector.into();
        let ones = super::super::util::one_indices_iterator(bitvector.iter());
        NearestNeighbourDictionary::<SAMPLE_OFFSET>::new(ones, bitvector.len())
    }

    #[test]
    fn rank() {
        let step = 120;
        let nnd = make_nnd(200, step);

        let mut i = 0;
        while i < nnd.len() {
            assert_eq!(nnd.rank(i), i / step);
            i += step;
        }
    }

    #[test]
    fn select() {
        let step = 120;
        let nnd = make_nnd(2000, step);

        let mut i = 0;
        let mut position = 0;
        while position < nnd.len() {
            assert_eq!(nnd.select(i) % step, 0);
            i += 1;
            position += step;
        }
    }

    #[test]
    fn previous() {
        let step = 120;
        let nnd = make_nnd(2000, step);

        let mut position = 0;
        while position < nnd.len() {
            assert_eq!(nnd.previous(position + step / 3), position);
            position += step;
        }
    }

    #[test]
    fn next() {
        let step = 120;
        let nnd = make_nnd(2000, step);

        let mut position = 0;
        while position < nnd.len() - step {
            assert_eq!(nnd.next(position + 2 * step / 3), position + step);
            position += step;
        }
    }
}


// Code by Martin Kostadinov adapted to Rust from the C++ library sdsl-lite.

pub const LOOKUP: LookupTables = LookupTables::new();
pub const POSITION_NULL: u32 = 9;

pub struct LookupTables {
    /// Given an 8-bit word interpreted as a parenthesis sequence where 1 is an opening parenthesis
    /// and 0 is a closing parenthesis gives the excess value.
    pub excess: [i8; 256],

    /// Given an excess value x in [1,8] and a 8-bit word w interpreted as parentheses sequence
    /// max_match_pos_packed[w]:[(x-1)*4,x*4] contains the maximal position, where right excess
    /// value -x is reached and 9, if there is no such position.
    pub packed_position_for_negative_excess: [u32; 256],
}

#[inline]
pub const fn get_packed_position(word: u32, offset: u32) -> u32 {
    (word >> (offset * 4)) & 0xf
}

#[inline]
pub const fn set_packed_position(word: u32, value: u32, offset: u32) -> u32 {
    let stepped_offset = offset * 4;
    let mask = !(0xf << stepped_offset);
    (word & mask) | ((value & 0xf) << stepped_offset)
}

impl LookupTables {
    const fn new() -> Self {
        let mut excess = [0i8; 256];
        let mut packed_position_for_negative_excess = [0u32; 256];

        const PACKED_NULL: u32 = 0x99999999;

        let mut excess_value: i8;
        let mut negative_excess_value: i8;
        let mut packed_position_values: u32;

        let mut word = 0;
        let mut bit_position;

        while word < 256 {
            excess_value = 0;
            negative_excess_value = 0;
            packed_position_values = PACKED_NULL;

            bit_position = 0;
            while bit_position < 8 {
                // Increase on opening bracket, decrease on closing.
                let bit_value = if word & (1 << bit_position) == 0 { 1 } else { 0 };
                excess_value += 1 - 2 * bit_value;

                // Increase on closing bracket, decrease on opening.
                let reverse_bit_value = if word & (1 << (7 - bit_position)) != 0 { 1 } else { 0 };
                negative_excess_value += 1 - 2 * reverse_bit_value;
                if negative_excess_value < 0 {
                    // The possible values for negative_excess_value are in [-8, -1]. An additoinal
                    // -1 after the negation of that region is needed to make the index 0-based
                    // with respect to the packed positions.
                    let offset = ((-negative_excess_value) - 1) as u32;
                    if get_packed_position(packed_position_values, offset) == 9 {
                        packed_position_values = set_packed_position(
                            packed_position_values,
                            7 - bit_position as u32,
                            offset
                        );
                    }
                }

                bit_position += 1;
            }

            excess[word as usize] = excess_value;
            packed_position_for_negative_excess[word as usize] = packed_position_values;

            word += 1;
        }

        Self {
            excess,
            packed_position_for_negative_excess
        }
    }
}

/// Assumes that the bits in the bitvector are ordered such that the least significant bit in each
/// bucket has the lowest index.
fn get_bit(bit_vector: &[u64], offset: usize) -> bool {
    let bucket_index = offset >> 6; // Equivalent to division by 64.
    let bit_offset = offset & 0b111111; // Equivalent to modulo 64.
    ((bit_vector[bucket_index] >> bit_offset) & 1) != 0
}

/// Assumes that the bits in the bitvector are ordered such that the least significant bit in each
/// bucket has the lowest index.
fn get_8_bit_word(bit_vector: &[u64], offset: usize) -> u8 {
    let bucket_index = offset >> 6;
    let bit_offset = offset & 0b111111;
    let offset_word = bit_vector[bucket_index] >> bit_offset;
    (offset_word & 0xff) as u8
}

pub fn try_to_find_opening_in_block(
    balanced_parenthesis: &[u64],
    start_index: usize,
    opening_excess: usize,
    block_size: usize
) -> Option<usize> {
    let mut current_excess: i64 = 0;
    let target_excess: i64 = opening_excess as i64;

    //
    // This function searches for a matching opening parenthesis in from a given bit index to the
    // beginning of its corresponding block.
    //
    //  v index of first bit in the block
    //  v
    //  v    8 bits
    //  v  |--------|
    // [...|........|........|........|....s...]
    //     ^                          ^    ^ - start index of the search
    //     ^                          ^ - right border
    //     ^ - left border
    //

    let index_of_first_bit_in_block = ((start_index / block_size) * block_size) as isize;

    // Floor to the nearest multiple of 8.
    let right_border = (start_index & !(0b111)) as isize;

    // Ceil to the nearest multiple of 8.
    let left_border = (index_of_first_bit_in_block + 7) & !(0b111);

    // Iterator variable over the region.
    let mut it: isize;

    // This loop is for the region from the starting index of the search to the right border
    // corrected for the beginning of the block if the right border happens to be beyond it.
    it = start_index as isize;
    let bound = right_border.max(index_of_first_bit_in_block);
    while it >= bound {
        if get_bit(balanced_parenthesis, it as usize) {
            current_excess += 1;
            if current_excess == target_excess {
                return Some(it as usize);
            }
        } else {
            current_excess -= 1;
        }
        it -= 1;
    }

    // This loop is for the region from the right border of the search to the left border with a
    // step size of 8.
    it = right_border - 8;
    while it >= left_border {
        let _8_bit_window = get_8_bit_word(balanced_parenthesis, it as usize);
        if target_excess - current_excess <= 8 {
            debug_assert!(
                target_excess - current_excess > 0,
                "If the current excess is equal or greater than the target excess then \
                 we must have found the parenthesis previously. The current excess changes \
                 by only 1 (on a step of 1 bit over the bitvector) so it is \"continuous\" \
                 over the integers."
            );
            //
            // target == current iff current - target == 0
            //
            // current - target = delta
            // target - current = (-delta)
            // 
            // If there is such a position which has a negative delta equal to the difference of
            // the target excess and the current excess, then at that position the current excess
            // would be equal to the target one i.e. we have found the matching opening
            // parenthesis.
            //
            let packed_positions = LOOKUP.packed_position_for_negative_excess[_8_bit_window as usize];
            // An additional -1 to make the index 0 based.
            let index_in_packed_position = (target_excess - current_excess - 1) as u32;
            let position = get_packed_position(packed_positions, index_in_packed_position);
            if position != POSITION_NULL {
                return Some(it as usize + position as usize);
            }
        }
        current_excess += LOOKUP.excess[_8_bit_window as usize] as i64;
        it -= 8;
    }

    // This loop is for the region from the left border to the index of the first bit in the block.
    it = left_border.min(right_border) - 1;
    while it >= index_of_first_bit_in_block {
        if get_bit(balanced_parenthesis, it as usize) {
            current_excess += 1;
            if current_excess == target_excess {
                return Some(it as usize);
            }
        } else {
            current_excess -= 1;
        }
        it -= 1;
    }

    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lookup_sanity_check() {
        let packed_positions = LOOKUP.packed_position_for_negative_excess[0xff];
        let mut current_position = packed_positions;
        for i in 0..8 {
            assert_eq!((current_position & 0xf), 7-i);
            current_position >>= 4;
        }
    }

    #[test]
    fn find_opening() {
        // Parenthesis sequence:
        // 0         1         2         3         4         5         6         7         8
        // 0123456789012345678901234567890123456789012345678901234567890123456789012345678901
        // (((()()))()()()()(((())()()))(())((()()())())(((((()))())())((()())(())()()()())))
        // 1111010001010101011110010100011001110101001001111110001001001110100110010101010000
        //
        // 1111010001010101011110010100011001110101001001111110001001001110 |
        // 100110010101010000 0000000000000000000000000000000000000000000000 
        //
        
        let bit_vector = &[
            0b1111010001010101011110010100011001110101001001111110001001001110_u64.reverse_bits(),
            0b1001100101010100000000000000000000000000000000000000000000000000_u64.reverse_bits(),
        ];

        assert_eq!(try_to_find_opening_in_block(bit_vector, 28 - 1, 1, 15), Some(17));
        assert_eq!(try_to_find_opening_in_block(bit_vector, 81 - 1, 1, 15), None);
        assert_eq!(try_to_find_opening_in_block(bit_vector, 81 - 1, 1, 82), Some(0));
        assert_eq!(try_to_find_opening_in_block(bit_vector, 80 - 1, 1, 40), Some(45));
    }
}



fn get_int<const BIT_WIDTH: usize>(data: &[u64], i: usize) -> usize {
    let bit_idx = i * BIT_WIDTH; 
    let word_idx = bit_idx / 64;
    let word_offset = bit_idx % 64; // Index of the least sigfinicant bit of the bitslice that is updated
    if word_offset + BIT_WIDTH <= 64 { // Int fits in this word
        let mask = (1_u64 << BIT_WIDTH) - 1;
        let bits = (data[word_idx] & (mask << word_offset)) >> word_offset;
        bits as usize
    } else { // Combine bits from two words

        let n_bits1 = 64 - word_offset; // All of the highest-order bits in the first word
        let n_bits2 = BIT_WIDTH - n_bits1; // Rest of the bits from the start of the second word
        debug_assert!(n_bits1 + n_bits2 == BIT_WIDTH);

        let x1 = data[word_idx] >> word_offset; // Tail of the first word
        let x2 = data[word_idx + 1] & ((1_u64 << n_bits2) - 1); // Head of the second word

        (x1 | (x2 << n_bits1)) as usize // Piece together
    }
}


fn set_int<const BIT_WIDTH: usize>(data: &mut [u64], i: usize, x: usize) {
    debug_assert!((x as u64) < (1_u64 << BIT_WIDTH));
    let bit_idx = i * BIT_WIDTH; 
    let word_idx = bit_idx / 64;
    let word_offset = bit_idx % 64; // Index of the least sigfinicant bit of the bitslice that is updated
    if word_offset + BIT_WIDTH <= 64 { // Int fits in this word
        let mask = (1_u64 << BIT_WIDTH) - 1; // Hopefully computed at compile time
        data[word_idx] &= !(mask << word_offset); // Clear the bits
        data[word_idx] |= (x as u64) << word_offset ; // Set new bits
    } else { // Combine bits from two words
        let n_bits1 = 64 - word_offset; // All of the highest-order bits in the first word
        let n_bits2 = BIT_WIDTH - n_bits1; // Rest of the bits from the start of the second word

        let mask1 = (1_u64 << n_bits1) - 1;
        let clearmask1 = !(mask1 << word_offset);
        let setmask1 = (x as u64 & mask1) << word_offset;

        data[word_idx] &= clearmask1; // Clear the bits
        data[word_idx] |= setmask1; // Set the bits

        let mask2 = (1_u64 << n_bits2) - 1;
        let clearmask2 = !mask2;
        let setmask2 = x as u64 >> n_bits1;

        data[word_idx + 1] &= clearmask2; // Clear the bits
        data[word_idx + 1] |= setmask2; // Set the bits
    }
}

struct CompactIntVector<const BIT_WIDTH: usize> {
    data: Vec<u64>,
    n_elements: usize,
}

struct CompactIntVectorMutSlice<'a, const BIT_WIDTH: usize> {
    data: &'a mut [u64],
    n_elements: usize,
}

impl<const BIT_WIDTH: usize> CompactIntVectorMutSlice<'_, BIT_WIDTH> {

    fn get(&self, i: usize) -> usize {
        debug_assert!(i < self.len());
        get_int::<BIT_WIDTH>(self.data, i)
    }

    fn set(&mut self, i: usize, x: usize) {
        debug_assert!(i < self.len());
        debug_assert!(x <= self.max_allowed_value());
        set_int::<BIT_WIDTH>(self.data, i, x)
    }

    fn len(&self) -> usize {
        self.n_elements
    }

    // Hopefully this compiles to just a constant
    fn max_allowed_value(&self) -> usize {
        ((1_u64 << BIT_WIDTH) - 1) as usize
    }
}


impl<const BIT_WIDTH: usize> CompactIntVector<BIT_WIDTH> {
    fn get(&self, i: usize) -> usize {
        debug_assert!(i < self.len());
        get_int::<BIT_WIDTH>(&self.data, i)
    }

    fn set(&mut self, i: usize, x: usize) {
        debug_assert!(i < self.len());
        debug_assert!(x <= self.max_allowed_value());
        set_int::<BIT_WIDTH>(&mut self.data, i, x)
    }

    fn split_to_mut_ranges(&mut self, n_ranges: usize) {
        unimplemented!();
    }

    fn len(&self) -> usize {
        self.n_elements
    }

    // Hopefully this compiles to just a constant
    fn max_allowed_value(&self) -> usize {
        ((1_u64 << BIT_WIDTH) - 1) as usize
    }

    // Initializes all elements to zeros
    fn new(len: usize) -> CompactIntVector<BIT_WIDTH> {
        assert!(BIT_WIDTH > 0);
        CompactIntVector{data: vec![0_u64; (len * BIT_WIDTH).div_ceil(64)], n_elements: len}
    }
}

#[cfg(test)]
mod tests {
    use super::CompactIntVector;

    #[test]
    fn set_and_get() {
        let len = 300;

        // Bit width 5 does not divide 64 so we get elements
        // spanning two words.
        let mut v = CompactIntVector::<5>::new(len);
        for i in 0..v.len() {
            v.set(i, i % v.max_allowed_value()); 
        }
        for i in 0..v.len() {
            let x = i % v.max_allowed_value();
            assert_eq!(v.get(i), x) ;
        }
    }
}
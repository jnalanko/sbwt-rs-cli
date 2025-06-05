use std::mem::uninitialized;

struct CompactIntVector {
    data: Vec<u64>,
    n_elements: usize,
    bit_width: usize, // needs to be less than 64
}

impl CompactIntVector {
    fn get(&self, i: usize) -> usize {
        let bit_idx = i * self.bit_width; 
        let word_idx = bit_idx / 64;
        let word_offset = bit_idx % 64; // Index of the least sigfinicant bit of the bitslice that is updated
        if word_offset + self.bit_width <= 64 { // Int fits in this word
            let mask = (1_u64 << self.bit_width) - 1;
            let bits = (self.data[word_idx] & (mask << word_offset)) >> word_offset;
            bits as usize
        } else { // Combine bits from two words

            let n_bits1 = 64 - word_offset; // All of the highest-order bits in the first word
            let n_bits2 = self.bit_width - n_bits1; // Rest of the bits from the start of the second word

            let x1 = self.data[word_idx] >> word_offset; // Tail of the first word
            let x2 = self.data[word_idx + 1] & ((1_u64 << n_bits2) - 1); // Head of the second word

            (x1 | (x2 << n_bits1)) as usize // Piece together
        }
    }

    fn set(&mut self, i: usize, x: usize) {
        let bit_idx = i * self.bit_width; 
        let word_idx = bit_idx / 64;
        let word_offset = bit_idx % 64; // Index of the least sigfinicant bit of the bitslice that is updated
        if word_offset + self.bit_width <= 64 { // Int fits in this word
            let mask = (1_u64 << self.bit_width) - 1;
            self.data[word_idx] &= !(mask << word_offset); // Clear the bits
            self.data[word_idx] |= (x as u64) << word_offset ; // Set new bits
        } else { // Combine bits from two words
            let n_bits1 = 64 - word_offset; // All of the highest-order bits in the first word
            let n_bits2 = self.bit_width - n_bits1; // Rest of the bits from the start of the second word

            let mask1 = (1_u64 << n_bits1) - 1;
            let clearmask1 = !(mask1 << word_offset);
            let setmask1 = (x as u64 & mask1) << word_offset;

            self.data[word_idx] &= clearmask1; // Clear the bits
            self.data[word_idx] |= setmask1; // Set the bits

            let mask2 = (1_u64 << n_bits2) - 1;
            let clearmask2 = !mask2;
            let setmask2 = (x as u64 & mask2);

            self.data[word_idx + 1] &= clearmask2; // Clear the bits
            self.data[word_idx + 1] |= setmask2; // Set the bits
        }
    }

    fn split_to_mut_ranges(&mut self) {
        unimplemented!();
    }

    fn len(&self) -> usize {
        self.n_elements
    }

    fn max_allowed_value(&self) -> usize {
        ((1_u64 << self.bit_width) - 1) as usize
    }

    fn new(len: usize, bit_width: usize) -> CompactIntVector {
        assert!(bit_width > 0);
        CompactIntVector{data: vec![0_u64; (len * bit_width).div_ceil(64)], n_elements: len, bit_width}
    }
}

#[cfg(test)]
mod tests {
    use super::CompactIntVector;

    #[test]
    fn set_and_get() {
        let width = 5; // Does not divide 64 so we get elements spanning two words
        let len = 300;

        let mut v = CompactIntVector::new(len, width);
        for i in 0..v.len() {
            v.set(i, i % v.max_allowed_value()); 
        }
        for i in 0..v.len() {
            let x = i % v.max_allowed_value();
            assert_eq!(v.get(i), x) ;
        }
    }
}
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
        let word_offset = bit_idx % 64;
        if word_offset + self.bit_width <= 64 { // Int fits in this word
            let mask = (1_u64 << self.bit_width) - 1;
            let bits = (self.data[word_idx] & (mask << word_offset)) >> word_offset;
            bits as usize
        } else { // Combine bits from two words
            unimplemented!();
        }
    }

    fn set(&mut self, i: usize, x: usize) {
        let bit_idx = i * self.bit_width; 
        let word_idx = bit_idx / 64;
        let word_offset = bit_idx % 64;
        if word_offset + self.bit_width <= 64 { // Int fits in this word
            let mask = (1_u64 << self.bit_width) - 1;
            self.data[word_idx] &= !(mask << word_offset); // Clear the bits
            self.data[word_idx] |= (x as u64) << word_offset ; // Set new bits
        } else { // Combine bits from two words
            unimplemented!();
        }
    }

    fn len(&self) -> usize {
        self.n_elements
    }

    fn max_allowed_value(&self) -> usize {
        ((1_u64 << self.bit_width) - 1) as usize
    }

    fn new(len: usize, bit_width: usize) -> CompactIntVector {
        assert!(bit_width > 0);
        CompactIntVector{data: vec![0_u64; len.div_ceil(64)], n_elements: len, bit_width}
    }
}

#[cfg(test)]
mod tests {
    use super::CompactIntVector;

    #[test]
    fn set_and_get() {
        let width = 4;
        let mut v = CompactIntVector::new(300, width);
        for i in 0..v.len() {
            v.set(i, i % v.max_allowed_value()); 
        }
        for i in 0..v.len() {
            let x = i % v.max_allowed_value();
            assert_eq!(v.get(i), x) ;
        }

    }
}
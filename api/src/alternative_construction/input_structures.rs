
use simple_sds_sbwt::bit_vector::BitVector;
use simple_sds_sbwt::ops::{BitVec, Rank, Select};
use simple_sds_sbwt::serialize::Serialize;

pub const TABLE_SIZE: usize = u8::MAX as usize + 1;
pub const CHAR_TO_INDEX: [usize; TABLE_SIZE] = make_char_to_index_table();
pub const INDEX_TO_CHAR: &[u8] = b"$ACGT#";

const fn make_char_to_index_table() -> [usize; TABLE_SIZE] {
    let mut table = [u8::MAX as usize; TABLE_SIZE];
    table[b'$' as usize] = 0;
    table[b'A' as usize] = 1;
    table[b'C' as usize] = 2;
    table[b'G' as usize] = 3;
    table[b'T' as usize] = 4;
    table[b'#' as usize] = 5;
    table
}

#[inline(always)]
pub fn char_index(byte: u8) -> usize {
    CHAR_TO_INDEX[byte as usize]
}

pub struct Bwt {
    pub data: Vec<BitVector>,
    pub counts: [usize; 5],
}

impl Bwt {
    pub fn new(data: Vec<BitVector>) -> Self {
        let mut counts = [0_usize; 5];
        counts[0] = 1; // '#'
        for i in 1..5 {
            counts[i] = counts[i - 1] + data[i - 1].count_ones();
        }
        Self {
            data,
            counts
        }
    }

    #[inline]
    pub fn get_char_index(&self, index: usize) -> usize {
        assert!(index < self.data[0].len());
        for char_index in 0..self.data.len() {
            if self.data[char_index].get(index) {
                return char_index;
            }
        }
        INDEX_TO_CHAR.len() - 1
    }

    #[inline]
    pub fn character(&self, index: usize) -> u8 {
        let char_index = self.get_char_index(index);
        INDEX_TO_CHAR[char_index]
    }

    #[inline]
    pub fn lf_step(&self, index: usize) -> (usize, u8) {
        let char_index = self.get_char_index(index);
        let character = INDEX_TO_CHAR[char_index];
        if char_index < self.counts.len() {
            let order = self.counts[char_index] + self.data[char_index].rank(index);
            (order, character)
        } else {
            (0, character)
        }
    }

    #[inline]
    pub fn first_character_of_suffix(&self, index: usize) -> u8 {
        for char_index in (0..self.counts.len()).rev() {
            if index >= self.counts[char_index] {
                return INDEX_TO_CHAR[char_index];
            }
        }
        INDEX_TO_CHAR[INDEX_TO_CHAR.len() - 1]
    }

    #[inline]
    pub fn inverse_lf_step(&self, index: usize) -> usize {
        assert!(index < self.data[0].len());
        for char_index in (0..self.counts.len()).rev() {
            if index >= self.counts[char_index] {
                let rank_within_count = index - self.counts[char_index];
                let order = self.data[char_index].select(rank_within_count)
                    .expect("The given bit should exist.");
                return order;
            }
        }
        0
    }

    pub fn serialize<W: std::io::Write>(&self, output: &mut W) -> std::io::Result<usize> {
        let mut written: usize = 0;
        for vector in &self.data {
            vector.serialize(output)?;
            written += vector.size_in_bytes();
        }
        Ok(written)
    }

    pub fn load<R: std::io::Read>(input: &mut R) -> std::io::Result<Self> {
        let mut data = vec![];
        for _ in 0..5 {
            data.push(BitVector::load(input)?);
        }
        let result = Self::new(data);
        Ok(result)
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.data[0].len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.data[0].len() == 0
    }
}

pub struct Lcp {
    data: Vec<u8>,
    len: usize,
    width: usize,
    offset: usize,
}

impl Lcp {
    #[inline]
    pub fn new<E>(data: Vec<u8>) -> Self {
        assert!(size_of::<E>() <= size_of::<usize>());
        assert!(data.len() % size_of::<E>() == 0);
        let len = data.len() / size_of::<E>();
        Self {
            data,
            len,
            width: size_of::<E>(),
            offset: 0,
        }
    }

    pub fn new_with_width(data: Vec<u8>, width: usize) -> Self {
        assert!(data.len() % width == 0);
        let len = data.len() / width;
        Self {
            data,
            len,
            width,
            offset: 0
        }
    }

    #[inline]
    pub fn push(&mut self, value: usize) {
        self.len += 1;
        let bytes = &value.to_le_bytes()[..self.width];
        self.data.extend_from_slice(bytes);
    }

    #[inline]
    pub fn get(&self, index: usize) -> usize {
        let start = index * self.width;
        let end = start + self.width;
        if end > self.data.len() {
            return 0;
        }
        let mut bytes: [u8; size_of::<usize>()] = [0_u8; size_of::<usize>()];
        let src = &self.data[start..end];
        bytes[0..self.width].copy_from_slice(src);
        usize::from_le_bytes(bytes)
    }

    #[inline]
    pub fn set(&mut self, index: usize, value: usize) {
        let start = index * self.width;
        let end = start + self.width;
        if end >= self.data.len() {
            return;
        }
        let bytes = &value.to_le_bytes()[0..self.width];
        self.data[start..end].copy_from_slice(bytes);
    }

    #[inline]
    pub fn skip(&mut self, count: usize) {
        self.offset += self.width * count;
    }

    #[inline]
    pub fn reset(&mut self) {
        self.offset = 0;
    }

    #[inline]
    pub fn width(&self) -> usize {
        self.width
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.len
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }
}

impl Iterator for Lcp {
    type Item = usize;
    fn next(&mut self) -> Option<Self::Item> {
        if self.offset >= self.data.len() {
            return None;
        }
        let end = self.offset + self.width;
        let mut bytes: [u8; 8] = [0_u8; 8];
        let src = &self.data[self.offset..end];
        bytes[0..self.width].copy_from_slice(src);
        let value = usize::from_le_bytes(bytes);
        self.offset = end;
        Some(value)
    }
}

impl From<Lcp> for Vec<u8> {
    fn from(value: Lcp) -> Self {
        value.data
    }
}

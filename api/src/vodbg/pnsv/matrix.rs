// Code by Martin Kostadinov.

use super::Pnsv;

use simple_sds_sbwt::bit_vector::BitVector;
use simple_sds_sbwt::ops::{BitVec, Rank, Select};
use simple_sds_sbwt::raw_vector::{AccessRaw, RawVector};
use simple_sds_sbwt::serialize::Serialize;

pub const MAX_ROWS: usize = 12;

/// Supports previous/next smaller value queries for a range of target lengths.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Matrix {
    pub lower_bound: usize,
    pub upper_bound: usize,
    pub width: usize,
    pub rows: Vec<BitVector>,
}

impl Matrix {
    pub fn empty() -> Self {
        Self {
            lower_bound: 0,
            upper_bound: 0,
            width: 0,
            rows: vec![]
        }
    }

    pub fn from_iterator<T, I>(input: I, width: usize, lower_bound: usize, upper_bound: usize) -> Self
    where 
        T: Into<usize>,
        I: Iterator<Item = T>,
    {
        let mut rows_raw: Vec<RawVector> = Vec::with_capacity(upper_bound - lower_bound + 1);

        for _ in lower_bound..=upper_bound {
            rows_raw.push(RawVector::with_len(width, false));
        }

        let ten_percent = width / 10;
        let mut border = ten_percent;
        let mut percent_count = 0;

        for (index, item) in input.enumerate() {
            let value = item.into();
            for target_length in lower_bound..=upper_bound {
                let row_index = target_length - lower_bound;
                if value < target_length {
                    rows_raw[row_index].set_bit(index, true);
                }
            }
            if index >= border {
                border += ten_percent;
                percent_count += 1;
                log::info!("[Matrix::from_iterator] scanning... {}0%", percent_count);
            }
        }
        
        let mut rows: Vec<BitVector> = Vec::with_capacity(rows_raw.len());
        for (index, row) in rows_raw.into_iter().enumerate() {
            log::info!("[Matrix::from_iterator] building rank and select for row {}...", index);
            let mut bit_vector: BitVector = row.into();
            bit_vector.enable_rank();
            bit_vector.enable_select();
            rows.push(bit_vector);
        }

        Self {
            lower_bound,
            upper_bound,
            width,
            rows,
        }
    }

    pub fn previous(&self, index: usize, target_length: usize) -> usize {
        if index >= self.width {
            return self.width;
        }
        if target_length < self.lower_bound || target_length > self.upper_bound {
            return 0;
        }
        let row_index = target_length - self.lower_bound;
        let is_one = self.rows[row_index].get(index);
        if is_one {
            return index;
        }
        let one_rank = self.rows[row_index].rank(index);
        if one_rank == 0 {
            return 0;
        }
        // There is at least one smaller value before.
        self.rows[row_index].select(one_rank - 1).unwrap()
    }

    pub fn next(&self, index: usize, target_length: usize) -> usize {
        if index >= self.width {
            return self.width;
        }
        if target_length < self.lower_bound || target_length > self.upper_bound {
            return self.width;
        }
        let row_index = target_length - self.lower_bound;
        let is_one = self.rows[row_index].get(index);
        if is_one {
            return index;
        }
        let one_rank = self.rows[row_index].rank(index);
        if one_rank == self.rows[row_index].count_ones() {
            return self.width;
        }
        // There is at least one smaller value after.
        self.rows[row_index].select(one_rank).unwrap()
    }

    pub fn serialize<W: std::io::Write>(&self, out: &mut W) -> std::io::Result<usize> {
        let mut written: usize = 0;
        out.write_all(&(self.lower_bound as u64).to_le_bytes())?;
        out.write_all(&(self.upper_bound as u64).to_le_bytes())?;
        out.write_all(&(self.width as u64).to_le_bytes())?;
        let row_count = self.rows.len();
        out.write_all(&(row_count as u64).to_le_bytes())?;
        written += 4 * size_of::<u64>();

        for (index, row) in self.rows.iter().enumerate() {
            log::info!("[Matrix::serialize] serializing row {}...", index);
            row.serialize(out)?;
            written += row.size_in_bytes();
        }

        Ok(written)
    }

    pub fn load<R: std::io::Read>(input: &mut R) -> std::io::Result<Self> {
        use byteorder::{ReadBytesExt, LittleEndian};
        let lower_bound = input.read_u64::<LittleEndian>()? as usize;
        let upper_bound = input.read_u64::<LittleEndian>()? as usize;
        let width       = input.read_u64::<LittleEndian>()? as usize;
        let row_count   = input.read_u64::<LittleEndian>()? as usize;
        let mut rows = vec![];
        for i in 0..row_count {
            log::info!("[Matrix::load] loading level {}...", i);
            let row = BitVector::load(input)?;
            rows.push(row);
        }
        let result = Self {
            lower_bound,
            upper_bound,
            width,
            rows
        };
        Ok(result)
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.width
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.width == 0
    }
}

impl Pnsv for Matrix {
    #[inline]
    fn previous(&self, index: usize, target_length: usize) -> usize {
        self.previous(index, target_length)
    }

    #[inline]
    fn next(&self, index: usize, target_length: usize) -> usize {
        self.next(index, target_length)
    }

    fn max_target(&self) -> usize {
        self.upper_bound
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn serialize_and_load() {
        let items: &[usize] = &[
            2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9, 9,
            9, 9, 9, 9, 8, 8, 8, 8, 7, 7, 7, 7, 6, 6, 6, 6, 5, 5, 5, 5, 4, 4, 4, 4, 3, 3, 3, 3, 2, 2, 2, 2,
        ];
        let lower_bound = 4;
        let upper_bound = 6;
        let matrix = Matrix::from_iterator(items.iter().cloned(), items.len(), lower_bound, upper_bound);
        let mut buffer = Vec::<u8>::new();
        let written = matrix.serialize(&mut buffer).unwrap();
        assert_eq!(buffer.len(), written);
        let matrix_loaded = Matrix::load(&mut buffer.as_slice()).unwrap();
        assert_eq!(matrix, matrix_loaded);
    }


    fn previous(items: &[usize], index: usize, target_length: usize, lower_bound: usize, upper_bound: usize) -> usize {
        for i in (0..=index).rev() {
            let item = items[i].clamp(lower_bound - 1, upper_bound);
            if item < target_length {
                return i;
            }
        }
        0
    }

    fn next(items: &[usize], index: usize, target_length: usize, lower_bound: usize, upper_bound: usize) -> usize {
        for i in index..items.len() {
            let item = items[i].clamp(lower_bound - 1, upper_bound);
            if item < target_length {
                return i;
            }
        }
        items.len()
    }

    fn test_with_parameters(items: &[usize], lower_bound: usize, upper_bound: usize) {
        let matrix = Matrix::from_iterator(items.iter().cloned(), items.len(), lower_bound, upper_bound);
        for i in 0..items.len() {
            for target_length in lower_bound..=upper_bound {
                assert_eq!(
                    previous(items, i, target_length, lower_bound, upper_bound),
                    matrix.previous(i, target_length),
                    "previous; i: {}, target_length: {}",
                    i,
                    target_length
                );
                assert_eq!(
                    next(items, i, target_length, lower_bound, upper_bound),
                    matrix.next(i, target_length),
                    "next; i: {}, target_length: {}",
                    i,
                    target_length
                );
            }
        }
    }

    #[test]
    fn pnsv_matrix_all_01() {
        let items: &[usize] = &[
            2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9, 9,
            9, 9, 9, 9, 8, 8, 8, 8, 7, 7, 7, 7, 6, 6, 6, 6, 5, 5, 5, 5, 4, 4, 4, 4, 3, 3, 3, 3, 2, 2, 2, 2,
        ];
        test_with_parameters(items, 4, 6);
        test_with_parameters(items, 5, 7);
        test_with_parameters(items, 2, 7);
    }

    #[test]
    fn pnsv_matrix_all_02() {
        let items: &[usize] = &[
            2, 2, 2, 3, 3, 4, 3, 3, 1, 4, 6, 4, 5, 8, 8, 5, 6, 6, 7, 6, 7, 3, 2, 7, 8, 6, 6, 6, 9, 9, 9, 9,
            9, 2, 9, 9, 3, 8, 8, 8, 4, 4, 7, 7, 6, 5, 6, 6, 8, 8, 2, 5, 4, 4, 8, 8, 9, 3, 3, 3, 3, 2, 2, 2,
        ];
        test_with_parameters(items, 4, 6);
        test_with_parameters(items, 5, 7);
        test_with_parameters(items, 2, 7);
    }
}


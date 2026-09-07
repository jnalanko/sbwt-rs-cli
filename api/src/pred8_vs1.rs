// pred8_vs1.rs
//
// Optimized Rust implementation of Pred8vS1.
//
// Layout:
//
// X:
//   X[i] = rank of first element whose upper 8 bits >= i
//   X[nblocks] = n (sentinel)
//
// Y:
//   [u32 hint][32 low bytes]
//   [u32 hint][32 low bytes]
//   ...

use std::io::{self, Read, Write};

const BUCKET_SHIFT: usize = 8;

const HINT_SHIFT: usize = 5;
const HINT_INTERVAL: usize = 1 << HINT_SHIFT; // 32

const HINT_BYTES: usize = 4;
const Y_BLOCK_SIZE: usize = HINT_BYTES + HINT_INTERVAL; // 36

#[derive(Clone, Eq, PartialEq, Debug)]
pub struct Pred8vS1 {
    // Universe maximum.
    u: usize,
    // Number of elements.
    n: usize,

    // Number of 256-value buckets.
    nblocks: usize,

    // Bucket rank index.
    // Length = nblocks + 1
    // The last element is a sentinel:
    // x[nblocks] == n
    upper_level: Vec<u32>,

    // Packed hints + low bytes.
    // Every 32 values:
    // +0..3     : upper bytes hint
    // +4..35    : low bytes
    lower_level: Vec<u32>,

    // need when u>>8 >= 1<<32
    super_hint: Vec<u64>,
    // need when n >= 1<<32
    super_upper: Vec<u64>,
}

impl Pred8vS1 {
    pub fn from_sorted(data: &[u64]) -> Self {
        assert!(!data.is_empty());

        // First element in data must be 0
        // n means the universe of the lower layer is (n - 1) * 256, has n-1 buckets
        let n = data.len();
        let u = data[n - 1] as usize;

        let nblocks = ((u / 256) + 1 as usize) as usize;

        // First pass:
        // count elements per bucket.
        let mut upper_level = vec![0u32; nblocks as usize + 1];
        let mut super_hint = vec![0u64; (u >> (32 + 8)) + 1];
        let mut super_upper = vec![0u64; (n >> 32) + 1];

        for &v in data {
            let bucket = (v as usize) >> BUCKET_SHIFT;
            debug_assert!(bucket < nblocks);

            upper_level[bucket] += 1;
        }
        // Allocate Y.
        let nhints = n.div_ceil(HINT_INTERVAL);
        let num_bytes = n + nhints * HINT_BYTES;
        let y_size = num_bytes.div_ceil(4);
        let mut lower_level = vec![0u32; y_size];

        // Second pass:
        // convert bucket counts into offsets
        // and write low bits.
        let mut yi: usize = 0;
        let mut y_pos: usize = 0;
        let mut i: usize = 0;
        let mut hint_counter = 0;
        let mut next_bucket_thrs = 1;
        let mut next_byte_thrs = 1;
        while i < n {
            let v = data[i];
            let bucket = (v as usize) >> BUCKET_SHIFT;
            // super blocking hint
            if bucket >= next_bucket_thrs << 32 {
                super_hint[next_bucket_thrs] = hint_counter as u64;
                next_bucket_thrs += 1;
            }

            let bucket_count = upper_level[bucket];
            // X now stores the rank offset.
            upper_level[bucket] = yi as u32;
            if yi >= next_byte_thrs << 32 {
                super_upper[next_byte_thrs] = bucket as u64;
                next_byte_thrs += 1;
            }
            if bucket_count != 0 {
                for _ in 0..bucket_count {
                    if (yi & 31) == 0 {
                        lower_level[(hint_counter * Y_BLOCK_SIZE) / 4] = bucket as u32;
                        y_pos += HINT_BYTES;
                        hint_counter += 1;
                    }
                    let y_bytes: &mut [u8] = bytemuck::cast_slice_mut(&mut lower_level);

                    y_bytes[y_pos] = (data[i] & 255) as u8;

                    y_pos += 1;
                    yi += 1;
                    i += 1;
                }
            }
        }

        debug_assert_eq!(y_pos.div_ceil(4), lower_level.len());

        // Explicit sentinel.
        // X[nblocks] == number of elements
        upper_level[nblocks] = yi as u32;
        // Propagate empty buckets backwards.
        let mut previous = upper_level[nblocks];

        for i in (1..nblocks).rev() {
            if upper_level[i] == 0 {
                upper_level[i] = previous;
            } else {
                previous = upper_level[i];
            }
        }

        log::info!("Completed Pred8s1 construction with n {n} elements, with universe {u} and super hint level of size {} and super upper level of size {}, with {nblocks} buckets",
            super_hint.len(),
            super_upper.len(),
            );

        Self {
            u,
            n,
            nblocks,
            upper_level,
            lower_level,
            super_hint,
            super_upper,
        }
    }

    // Hot query path.
    #[inline(always)]
    pub fn select1(&self, sq: usize) -> (usize, u32) {
        if sq >= self.n {
            return ((self.u - self.n) as usize, 0);
        }
        // this means we are querying the number of elements in total in the
        if sq == self.n - 1 {
            return ((self.u - sq) as usize, 0);
        }
        // Locate hint block.
        let hint_block = sq >> HINT_SHIFT;
        let hint_pos = hint_block * Y_BLOCK_SIZE;
        let low_pos = hint_pos + HINT_BYTES + (sq & 31);
        let y_bytes: &[u8] = bytemuck::cast_slice(&self.lower_level);
        let low = y_bytes[low_pos] as u64;

        let hint_number = hint_pos / 4;
        let mut hint = self.lower_level[hint_number] as usize;
        let mut hint_offset = 0;
        let mut h = 0;
        while hint_number < self.super_hint[h] as usize {
            h += 1;
            hint_offset += 1 << 32;
        }
        hint += hint_offset;
        // Next element.
        let next_sq = sq + 1;
        let next_hint_pos = (next_sq >> HINT_SHIFT) * Y_BLOCK_SIZE;

        let next_low_pos = next_hint_pos + HINT_BYTES + (next_sq & 31);

        let next_low = y_bytes[next_low_pos] as u64;
        h = 1;
        let mut bucket_count_offset = 0;
        while sq >= (h << 32) && hint >= self.super_upper[h] as usize {
            bucket_count_offset += 1 << 32;
            h += 1;
        }
        // Scan upper_level for the correct bucket.
        // The sentinel guarantees termination.
        let mut j = 1usize;

        let mut bucket_count = self.upper_level[hint as usize + j];
        if sq > (h << 32) && hint + j >= self.super_upper[h] as usize {
            bucket_count_offset += 1 << 32;
            h += 1;
        }
        while (bucket_count as usize + bucket_count_offset) <= sq {
            j += 1;
            bucket_count = self.upper_level[hint as usize + j];
            if sq > (h << 32) && hint + j >= self.super_upper[h] as usize {
                bucket_count_offset += 1 << 32;
                h += 1;
            }
        }
        let mut j2 = j - 1;
        let mut bucket_count_2 = self.upper_level[hint as usize + j2];
        while (bucket_count_2 as usize + bucket_count_offset) <= next_sq {
            //p2 = p2.add(1);
            j2 += 1;
            bucket_count_2 = self.upper_level[hint as usize + j2];
            if next_sq > (h << 32) && hint + j2 >= self.super_upper[h] as usize {
                bucket_count_offset += 1 << 32;
                h += 1;
            }
        }

        let upper = ((hint + j - 1) << BUCKET_SHIFT) as u64;
        let next_upper = ((hint + j2 - 1) << BUCKET_SHIFT) as u64;
        let curr_pos = upper + low;
        let next_pos = next_upper + next_low;

        debug_assert!(next_pos > curr_pos);
        let bucket_count = next_pos - curr_pos - 1;
        let mut result = curr_pos - sq as u64;

        if bucket_count == 0 {
            result = result.saturating_sub(1); // it can only underflow if the first bucket is empty, which we assume it never will be due to construction of Pino, then the result is pointing correctly to itself and needs not to be subtracted by 1
        }
        (result as usize, bucket_count as u32)
    }

    // Number of bytes occupied by this structure.
    #[inline]
    pub fn size_in_bytes(&self) -> usize {
        3 * std::mem::size_of::<u64>()
            + self.upper_level.len() * std::mem::size_of::<u32>()
            + self.lower_level.len() * std::mem::size_of::<u32>()
            + self.super_upper.len() * std::mem::size_of::<u64>()
            + self.super_hint.len() * std::mem::size_of::<u64>()
    }

    pub fn num_buckets(&self) -> usize {
        self.nblocks - 1
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.n
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.n == 0
    }

    #[inline]
    pub fn universe(&self) -> usize {
        self.u
    }

    // Layout:
    // u64 u
    // u64 n
    // u64 nblocks
    // u32 upper_level[0..=nblocks]   // includes sentinel x[nblocks]
    // u8  lower_level[n + num_hints*4ß]
    #[inline]
    pub fn serialize<W: Write>(&self, mut w: W) -> io::Result<usize> {
        debug_assert_eq!(self.upper_level.len(), self.nblocks + 1);

        w.write_all(&self.u.to_le_bytes())?;
        w.write_all(&self.n.to_le_bytes())?;
        w.write_all(&self.nblocks.to_le_bytes())?;

        // X is stored as little-endian u32 values.
        // This assumes the target architecture is little endian.
        w.write_all(bytemuck::cast_slice(&self.upper_level))?;

        // Y contains both hints and payload bytes.
        w.write_all(bytemuck::cast_slice(&self.lower_level))?;

        w.write_all(bytemuck::cast_slice(&self.super_hint))?;

        w.write_all(bytemuck::cast_slice(&self.super_upper))?;

        Ok(self.size_in_bytes())
    }

    pub fn load<R: Read>(mut r: R) -> io::Result<Self> {
        fn read_u64<R: Read>(r: &mut R) -> io::Result<u64> {
            let mut buf = [0u8; 8];
            r.read_exact(&mut buf)?;
            Ok(u64::from_le_bytes(buf))
        }

        let u = read_u64(&mut r)? as usize;
        let n = read_u64(&mut r)? as usize;
        let nblocks = read_u64(&mut r)? as usize;

        let mut x = vec![0u32; nblocks + 1];

        // Read X as raw little-endian u32 storage.
        // Matches serialize().
        r.read_exact(bytemuck::cast_slice_mut(x.as_mut_slice()))?;

        let nhints = n / HINT_INTERVAL + 1;

        let num_bytes = n + nhints * HINT_BYTES;
        let y_size = num_bytes.div_ceil(4);

        let mut y = vec![0u32; y_size];
        r.read_exact(bytemuck::cast_slice_mut(&mut y))?;

        let mut super_hint = vec![0u64; (u >> 32) + 1];
        r.read_exact(bytemuck::cast_slice_mut(&mut super_hint))?;

        let mut super_upper = vec![0u64; (n >> 32) + 1];
        r.read_exact(bytemuck::cast_slice_mut(&mut super_upper))?;

        Ok(Self {
            u,
            n,
            nblocks,
            upper_level: x,
            lower_level: y,
            super_hint,
            super_upper,
        })
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn basic_select() {
        let values = vec![0, 2, 5, 255, 256, 257, 500, 696, 819, 820];

        let pred = Pred8vS1::from_sorted(&values);
        assert_eq!(pred.len(), values.len());
        assert_eq!(pred.universe(), 820);

        for (rank, pair) in values.windows(2).enumerate() {
            let pos = pair[0];
            let next = pair[1];

            let (returned_pos, count) = pred.select1(rank);

            assert!(returned_pos <= pos as usize);
            assert_eq!(count as u64, next - pos - 1);
        }
        assert_eq!(pred.upper_level[pred.nblocks], pred.n as u32);
    }

    #[test]
    fn sparse_select() {
        let values = vec![0, 2, 5, 1255, 2256, 3257, 5500, 7696, 10819, 21820];

        let pred = Pred8vS1::from_sorted(&values);
        assert_eq!(pred.len(), values.len());

        for (rank, pair) in values.windows(2).enumerate() {
            let pos = pair[0];
            let next = pair[1];

            let (returned_pos, count) = pred.select1(rank);

            assert!(returned_pos <= pos as usize);
            assert_eq!(count as u64, next - pos - 1);
        }
        assert_eq!(pred.upper_level[pred.nblocks], pred.n as u32);
    }

    #[test]
    fn dense_select() {
        let values = vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11];

        let pred = Pred8vS1::from_sorted(&values);
        assert_eq!(pred.len(), values.len());

        for (rank, pair) in values.windows(2).enumerate() {
            let pos = pair[0];
            let next = pair[1];

            let (returned_pos, count) = pred.select1(rank);

            assert_eq!(returned_pos, (pos - rank as u64) as usize);
            assert_eq!(count as u64, next - pos - 1);
        }
        assert_eq!(pred.upper_level[pred.nblocks], pred.n as u32);
    }

    #[test]
    fn pow_of_two_select() {
        let values = vec![
            0,
            256,
            2 * 256,
            3 * 256,
            4 * 256,
            5 * 256,
            6 * 256,
            7 * 256,
            8 * 256,
            9 * 256,
            10 * 256,
            11 * 256,
            12 * 256,
            13 * 256,
            14 * 256,
            15 * 256,
            16 * 256,
            17 * 256,
            18 * 256,
            19 * 256,
            20 * 256,
            21 * 256,
        ];

        let pred = Pred8vS1::from_sorted(&values);
        assert_eq!(pred.len(), values.len());

        for (rank, pair) in values.windows(2).enumerate() {
            let pos = pair[0];
            let next = pair[1];

            let (returned_pos, count) = pred.select1(rank);
            //dbg!(returned_pos, count, rank);

            assert_eq!(returned_pos, (pos - rank as u64) as usize);
            assert_eq!(count as u64, next - pos - 1);
        }
        assert_eq!(pred.upper_level[pred.nblocks], pred.n as u32);
    }

    #[test]
    fn serialize_roundtrip() {
        let values = vec![0, 10, 20, 300, 10000];
        let pred = Pred8vS1::from_sorted(&values);
        let mut bytes = Vec::new();
        pred.serialize(&mut bytes).unwrap();
        let loaded = Pred8vS1::load(bytes.as_slice()).unwrap();
        assert_eq!(pred.u, loaded.u);
        assert_eq!(pred.n, loaded.n);
        assert_eq!(pred.upper_level, loaded.upper_level);
        assert_eq!(pred.lower_level, loaded.lower_level);
    }
}

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
use std::ptr;

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
    x: Vec<u32>,

    // Packed hints + low bytes.
    // Every 32 values:
    // +0..3     : upper bits hint
    // +4..35    : low bytes
    y: Vec<u8>,
}

#[inline(always)]
unsafe fn read_hint(ptr: *const u8) -> u32 {
    ptr::read_unaligned(ptr as *const u32)
}

#[inline(always)]
unsafe fn write_hint(ptr: *mut u8, value: u32) {
    ptr::write_unaligned(ptr as *mut u32, value);
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
        let mut x = vec![0u32; nblocks as usize + 1];

        for &v in data {
            let bucket = (v as usize) >> BUCKET_SHIFT;
            debug_assert!(bucket < nblocks);

            x[bucket] += 1;
        }
        // Allocate Y.
        let nhints = n.div_ceil(HINT_INTERVAL);
        let mut y = vec![0u8; n + nhints * HINT_BYTES];

        // Second pass:
        // convert bucket counts into offsets
        // and write low bits.
        let mut yi: u32 = 0;
        let mut y_pos: usize = 0;
        let mut i: usize = 0;
        unsafe {
            let y_ptr = y.as_mut_ptr();

            while i < n {
                let v = data[i];
                let bucket = (v as usize) >> BUCKET_SHIFT;
                let bucket_count = x[bucket];
                // X now stores the rank offset.
                x[bucket] = yi;
                if bucket_count != 0 {
                    for _ in 0..bucket_count {
                        if (yi & 31) == 0 {
                            write_hint(y_ptr.add(y_pos), bucket as u32);
                            y_pos += HINT_BYTES;
                        }
                        *y_ptr.add(y_pos) = (data[i] & 255) as u8;
                        y_pos += 1;
                        yi += 1;
                        i += 1;
                    }
                }
            }
        }

        debug_assert_eq!(y_pos, y.len());

        x[nblocks] = yi as u32;
        // Propagate empty buckets backwards.
        let mut previous = x[nblocks];

        for i in (1..nblocks).rev() {
            if x[i] == 0 {
                x[i] = previous;
            } else {
                previous = x[i];
            }
        }

        // Explicit sentinel.
        // This is the important invariant:
        // X[nblocks] == number of elements
        x[nblocks] = n as u32;

        Self {
            u,
            n,
            nblocks,
            x,
            y,
        }
    }

    // Hot query path.
    #[inline(always)]
    pub fn select1(&self, sq: usize) -> (isize, u32) {
        if sq >= self.n {
            return ((self.u - sq) as isize, 0);
        }
        // this means we are querying the number of elements in total in the
        if sq == self.n - 1 {
            return ((self.u - sq) as isize, 0);
        }
        unsafe {
            // Locate hint block.
            let hint_block = sq >> HINT_SHIFT;
            let hint_pos = hint_block * Y_BLOCK_SIZE;
            let low_pos = hint_pos + HINT_BYTES + (sq & 31);
            let hint = read_hint(self.y.as_ptr().add(hint_pos)) as usize;
            let low = *self.y.get_unchecked(low_pos) as u64;

            // Next element.
            let next_sq = sq + 1;
            let next_hint_pos = (next_sq >> HINT_SHIFT) * Y_BLOCK_SIZE;

            let next_low_pos = next_hint_pos + HINT_BYTES + (next_sq & 31);

            let next_low = *self.y.get_unchecked(next_low_pos) as u64;

            // Scan X.
            // The sentinel guarantees termination.
            let mut j = 1usize;
            let mut p = self.x.as_ptr().add(hint + j);
            while (*p as usize) <= sq {
                p = p.add(1);
                j += 1;
            }
            let mut p2 = self.x.as_ptr().add(hint + j - 1);
            let mut j2 = j - 1;
            while (*p2 as usize) <= next_sq {
                p2 = p2.add(1);
                j2 += 1;
            }
            // dbg!(*p2, j2, hint + j2, self.nblocks);
            //assert!(*p2 > next_sq as u32);
            let upper = ((hint + j - 1) << BUCKET_SHIFT) as u64;
            let next_upper = ((hint + j2 - 1) << BUCKET_SHIFT) as u64;
            let curr_pos = upper + low;
            let next_pos = next_upper + next_low;
            // dbg!(curr_pos, next_pos, sq, j, j2, next_sq);
            debug_assert!(next_pos > curr_pos);
            let bucket_count = next_pos - curr_pos - 1;
            let mut result = curr_pos - sq as u64;

            //dbg!(result, bucket_count, sq);

            if bucket_count == 0 {
                result = result.saturating_sub(1); // it can only underflow if the first bucket is empty, which we assume it never will be due to construction of Pino, then the result is pointing correctly to itself and needs not to be subtracted by 1
            }
            (result as isize, bucket_count as u32)
        }
    }

    // Number of bytes occupied by this structure.
    #[inline]
    pub fn size_in_bytes(&self) -> usize {
        3 * std::mem::size_of::<u64>() + self.x.len() * std::mem::size_of::<u32>() + self.y.len()
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
    // u32 x[0..=nblocks]   // includes sentinel x[nblocks]
    // u8  y[]
    #[inline]
    pub fn serialize<W: Write>(&self, mut w: W) -> io::Result<usize> {
        debug_assert_eq!(self.x.len(), self.nblocks + 1);

        w.write_all(&self.u.to_le_bytes())?;
        w.write_all(&(self.n as u64).to_le_bytes())?;
        w.write_all(&(self.nblocks as u64).to_le_bytes())?;

        // X is stored as little-endian u32 values.
        // This assumes the target architecture is little endian.
        unsafe {
            let bytes = std::slice::from_raw_parts(
                self.x.as_ptr() as *const u8,
                self.x.len() * std::mem::size_of::<u32>(),
            );
            w.write_all(bytes)?;
        }

        // Y contains both hints and payload bytes.
        w.write_all(&self.y)?;

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
        unsafe {
            let bytes = std::slice::from_raw_parts_mut(
                x.as_mut_ptr() as *mut u8,
                x.len() * std::mem::size_of::<u32>(),
            );
            r.read_exact(bytes)?;
        }

        let nhints = n / HINT_INTERVAL + 1;
        let mut y = vec![0u8; n + nhints * HINT_BYTES];
        r.read_exact(&mut y)?;

        Ok(Self {
            u,
            n,
            nblocks,
            x,
            y,
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

            assert!(returned_pos <= pos as isize);
            assert_eq!(count as u64, next - pos - 1);
        }
        assert_eq!(pred.x[pred.nblocks], pred.n as u32);
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

            assert!(returned_pos <= pos as isize);
            assert_eq!(count as u64, next - pos - 1);
        }
        assert_eq!(pred.x[pred.nblocks], pred.n as u32);
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

            assert_eq!(returned_pos, (pos - rank as u64) as isize);
            assert_eq!(count as u64, next - pos - 1);
        }
        assert_eq!(pred.x[pred.nblocks], pred.n as u32);
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

            assert_eq!(returned_pos, (pos - rank as u64) as isize);
            assert_eq!(count as u64, next - pos - 1);
        }
        assert_eq!(pred.x[pred.nblocks], pred.n as u32);
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
        assert_eq!(pred.x, loaded.x);
        assert_eq!(pred.y, loaded.y);
    }
}

// pred8_vpino_light.rs

use crate::pred8_vs1::Pred8vS1;

#[derive(Clone, Eq, PartialEq, Debug)]
pub struct Pred8vPino {
    min: usize,
    // maximum value - minimum value
    u: usize,
    n: usize,
    nblocks: usize,
    // low 8 bits of all values, bucket concatenated
    y: Vec<u8>,
    // compressed bucket directory
    index: Pred8vS1,
}

impl Pred8vPino {
    pub fn from_sorted(data: &[u64]) -> Self {
        assert!(!data.is_empty());

        let n = data.len();
        let min = data[0] as usize;
        let u = data[n - 1] as usize - min;

        let nblocks = ((u / 256) + 1) as usize;
        let mut x = vec![0u32; nblocks + 1];

        // Count elements per bucket.
        for &value in data {
            let v = value - min as u64;
            let bucket = (v >> 8) as usize;

            debug_assert!(bucket < nblocks);
            x[bucket] += 1;
        }

        // Build payload.
        let mut y = vec![0u8; n];

        let mut yi = 0usize;
        let mut i = 0usize;

        while i < n {
            let v = data[i] - min as u64;
            let bucket = (v >> 8) as usize;

            let bucket_count = x[bucket];

            // X now becomes bucket start offset.
            x[bucket] = yi as u32;

            for _ in 0..bucket_count {
                let v = data[i] - min as u64;
                y[yi] = (v & 255) as u8;

                yi += 1;
                i += 1;
            }
        }

        debug_assert_eq!(yi, n);

        // Sentinel.
        x[nblocks] = yi as u32;
        let mut previous = x[nblocks];
        for i in (1..nblocks).rev() {
            if x[i] == 0 {
                x[i] = previous;
            } else {
                previous = x[i];
            }
        }

        let mut block_indices = Vec::with_capacity(nblocks + 2);

        let mut bi = 0u64;
        block_indices.push(bi);

        for i in 0..nblocks {
            let x_i = x[i] as u64;
            let bcount = (x[i + 1] as u64).saturating_sub(x_i);

            bi += 1 + bcount;
            block_indices.push(bi);
        }
        block_indices.push(bi + 1);

        let index = Pred8vS1::from_sorted(&block_indices);

        Self {
            min,
            u,
            n,
            nblocks,
            y,
            index,
        }
    }

    /* pub fn from_sorted(data: &[u64]) -> Self {
        assert!(!data.is_empty());

        let n = data.len();
        let min = data[0];
        let u = data[n - 1] - min;

        // let nblocks = (u as usize >> 8) + (((u & 255) != 0) as usize);
        let nblocks = (u as usize >> 8) + 1;

        // Temporary construction structure.
        // Each bucket contains the low 8 bits.
        let mut buckets = vec![Vec::<u8>::new(); nblocks];

        for &value in data {
            let v = value - min;
            let bucket = (v >> 8) as usize;
            buckets[bucket].push((v & 255) as u8);
        }
        // Flatten payload and create
        // the bucket directory.
        // The +1 is the sentinel expansion
        // used by Pred8vS1.
        let mut y = Vec::with_capacity(n);
        let mut bucket_directory = Vec::with_capacity(nblocks + 1);
        let mut pos = 0u64;
        bucket_directory.push(pos);

        for bucket in &buckets {
            y.extend_from_slice(bucket);
            pos += 1 + bucket.len() as u64;
            bucket_directory.push(pos);
        }

        debug_assert_eq!(bucket_directory.len(), nblocks + 1);

        let index = Pred8vS1::from_sorted(&bucket_directory);

        Self {
            min,
            u,
            n,
            nblocks,
            y,
            index,
        }
    } */

    // Return:
    //   (position, exact_match)
    // position:
    //   rank of predecessor
    // exact_match:
    //   true iff key exists
    pub fn get_pred(&self, key: usize) -> (usize, bool) {
        if key < self.min {
            return (usize::MAX, false);
        }
        let v = key - self.min;
        // Query beyond universe.
        // The predecessor is the last value.
        if v > self.u {
            return (self.n - 1, false);
        }
        if v == self.u {
            return (self.n - 1, true);
        }
        let bucket = (v >> 8) as usize;
        let (start, count) = self.index.select1(bucket);
        let start = start as usize;

        // Empty bucket.
        // select1 returns the predecessor
        // position in this case.
        if count == 0 {
            return (start, false);
        }
        let low = (v & 255) as u8;
        let count = count as usize;

        // Scan the bucket.
        for i in 0..count {
            let value = self.y[start + i];

            if value >= low {
                let exact = value == low;
                return (start + i - ((!exact) as usize), exact);
            }
        }

        // All elements in bucket are smaller.
        (start + count - 1, false)
    }

    // Number of elements strictly smaller than key.
    pub fn rank_contains(&self, key: usize) -> (usize, bool) {
        if self.n == 0 {
            return (0, false);
        }

        let (pos, exact) = self.get_pred(key);

        if pos == usize::MAX {
            return (0, false);
        }

        (pos + 1 - exact as usize, exact)
    }

    #[inline]
    fn fill_bucket_window(
        &self,
        bucket: usize,
        local_begin: usize,
        local_end: usize,
        window_start: usize,
        words: &mut [u64],
    ) {
        if bucket >= self.nblocks {
            return;
        }
        let (start, count) = self.index.select1(bucket);

        if count == 0 {
            return;
        }

        let start = start as usize;
        let count = count as usize;
        let bucket_base = (bucket as u64) << 8;

        for i in 0..count {
            let local = bucket_base + self.y[start + i] as u64;

            if local < local_begin as u64 || local > local_end as u64 {
                continue;
            }

            let bit = (local + self.min as u64) - window_start as u64;

            words[(bit >> 6) as usize] |= 1u64 << (bit & 63);
        }
    }

    #[inline]
    pub fn window256_starting_from_key_word_aligned(&self, key: usize) -> [u64; 4] {
        let mut words = [0u64; 4];

        let window_start = key as u64 & !63u64;
        let window_end = window_start + 255;

        // Completely outside the represented universe.
        if window_end < self.min as u64 || window_start > (self.min + self.u) as u64 {
            return words;
        }

        // Convert to Pino's local coordinates.
        let local_begin = window_start.saturating_sub(self.min as u64);
        let local_end = window_end.min((self.min + self.u) as u64) - self.min as u64;

        let first_bucket = (local_begin >> 8) as usize;
        let last_bucket = (local_end >> 8) as usize;

        self.fill_bucket_window(
            first_bucket,
            local_begin as usize,
            local_end as usize,
            window_start as usize,
            &mut words,
        );

        if last_bucket != first_bucket {
            self.fill_bucket_window(
                last_bucket,
                local_begin as usize,
                local_end as usize,
                window_start as usize,
                &mut words,
            );
        }

        words
    }

    #[inline]
    pub fn get_words_in_range(&self, range: std::ops::Range<usize>) -> Vec<u64> {
        let len = (range.start & 63) + range.end - range.start;
        let num_words = (len + 64 - 1) / 64;
        let mut words: Vec<u64> = vec![0; num_words];
        let key = range.start;
        if key >= self.u {
            return words;
        }
        let range_window_start = key as u64 & !63u64;
        let range_window_end = range_window_start + num_words as u64 * 64;
        
        // Completely outside the represented universe.
        if range_window_end < self.min as u64 || range_window_start > (self.min + self.u) as u64 {
            return words;
        }
        
        // Convert to Pino's local coordinates.
        let local_begin = range_window_start.saturating_sub(self.min as u64);
        let local_end = range_window_end.min((self.u) as u64);
        dbg!(len, num_words, range_window_start, range_window_end, local_begin, local_end);

        let first_bucket = (local_begin >> 8) as usize;
        let last_bucket = (local_end >> 8) as usize;

        for b in first_bucket..=last_bucket {
            self.fill_bucket_window(
                b,
                local_begin as usize,
                local_end as usize,
                range_window_start as usize,
                words.as_mut_slice(),
            );
        }
        words
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.n
    }

    #[inline]
    pub fn universe(&self) -> usize {
        self.u
    }

    #[inline]
    pub fn min(&self) -> usize {
        self.min
    }

    pub fn size_in_bytes(&self) -> usize {
        std::mem::size_of::<u64>() * 3
            + std::mem::size_of::<usize>() * 2
            + self.y.len()
            + self.index.size_in_bytes()
    }
}

use std::io::{self, Read, Write};

impl Pred8vPino {
    pub fn serialize<W: Write>(&self, mut w: W) -> io::Result<usize> {
        // Layout:
        // u64 u
        // u64 n
        // u64 min
        // u64 nblocks
        // u8  y[n]
        // Pred8vS1 index

        debug_assert_eq!(self.y.len(), self.n);

        let mut written = 0usize;

        w.write_all(&self.u.to_le_bytes())?;
        written += 8;

        w.write_all(&(self.n as u64).to_le_bytes())?;
        written += 8;

        w.write_all(&self.min.to_le_bytes())?;
        written += 8;

        w.write_all(&(self.nblocks as u64).to_le_bytes())?;
        written += 8;

        w.write_all(&self.y)?;
        written += self.y.len();

        written += self.index.serialize(&mut w)?;

        Ok(written)
    }

    pub fn load<R: Read>(mut r: R) -> io::Result<Self> {
        fn read_u64<R: Read>(r: &mut R) -> io::Result<u64> {
            let mut buf = [0u8; 8];
            r.read_exact(&mut buf)?;
            Ok(u64::from_le_bytes(buf))
        }

        let u = read_u64(&mut r)? as usize;
        let n = read_u64(&mut r)? as usize;
        let min = read_u64(&mut r)? as usize;
        let nblocks = read_u64(&mut r)? as usize;

        let mut y = vec![0u8; n];
        r.read_exact(&mut y)?;

        let index = Pred8vS1::load(&mut r)?;

        debug_assert_eq!(index.len() - 1, nblocks + 1);

        Ok(Self {
            min,
            u,
            n,
            nblocks,
            y,
            index,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic_pred() {
        let values = vec![10, 20, 30, 300, 500];

        let pino = Pred8vPino::from_sorted(&values);

        assert_eq!(pino.get_pred(5), (usize::MAX, false));
        assert_eq!(pino.get_pred(10), (0, true));
        assert_eq!(pino.get_pred(15), (0, false));
        assert_eq!(pino.get_pred(20), (1, true));
        assert_eq!(pino.get_pred(299), (2, false));
        assert_eq!(pino.get_pred(300), (3, true));
        assert_eq!(pino.get_pred(1000), (4, false));
    }

    #[test]
    fn test_rank() {
        let values = vec![10, 20, 30, 300, 500];

        let pino = Pred8vPino::from_sorted(&values);

        assert_eq!(pino.rank_contains(0), (0, false));
        assert_eq!(pino.rank_contains(10), (0, true));
        assert_eq!(pino.rank_contains(11), (1, false));
        assert_eq!(pino.rank_contains(20), (1, true));
        assert_eq!(pino.rank_contains(21), (2, false));
        assert_eq!(pino.rank_contains(300), (3, true));
        assert_eq!(pino.rank_contains(301), (4, false));
        assert_eq!(pino.rank_contains(1000), (5, false));
    }

    #[test]
    fn test_empty_buckets() {
        // Values deliberately separated by large gaps.
        let values = vec![0, 1, 900];

        let pino = Pred8vPino::from_sorted(&values);

        assert_eq!(pino.get_pred(0), (0, true));
        assert_eq!(pino.get_pred(1), (1, true));

        // Empty bucket between 1 and 900.
        assert_eq!(pino.get_pred(500), (1, false));

        assert_eq!(pino.get_pred(900), (2, true));
    }

    #[test]
    fn test_window256_starting_from_key_alignment() {
        let mut values = vec![1000, 1010, 1064, 1100];

        let pino = Pred8vPino::from_sorted(&values);

        let key = 1030;

        let words = pino.window256_starting_from_key_word_aligned(key);

        // key is aligned down to 1024
        let start = key & !63;

        for value in &mut values[2..4] {
            let bit = *value - start as u64;

            assert_ne!(words[(bit >> 6) as usize] & (1u64 << (bit & 63)), 0);
        }
    }
    #[test]
    fn test_window256_starting_from_key_word_aligned() {
        let values = vec![1000, 1010, 1064, 1100, 1279, 1280];

        let pino = Pred8vPino::from_sorted(&values);

        let key = 1030;
        let start = key & !63; // 1024

        let words = pino.window256_starting_from_key_word_aligned(key);

        for bit in 0..256 {
            let value = start as u64 + bit;

            let expected = values.contains(&value);

            let word = (bit / 64) as usize;
            let offset = bit % 64;

            let actual = (words[word] & (1u64 << offset)) != 0;

            assert_eq!(
                actual, expected,
                "incorrect bit for value {} (bit {})",
                value, bit
            );
        }
    }

    #[test]
    fn test_single_element() {
        let pino = Pred8vPino::from_sorted(&[123]);

        assert_eq!(pino.get_pred(122), (usize::MAX, false));
        assert_eq!(pino.get_pred(123), (0, true));
        assert_eq!(pino.get_pred(124), (0, false));

        assert_eq!(pino.rank_contains(122), (0, false));
        assert_eq!(pino.rank_contains(123), (0, true));
        assert_eq!(pino.rank_contains(124), (1, false));
    }

    #[test]
    fn test_bucket_boundaries() {
        let values = vec![255, 256, 257, 511, 512];

        let pino = Pred8vPino::from_sorted(&values);

        for &v in &values {
            assert_eq!(pino.get_pred(v as usize).1, true);
        }

        assert_eq!(pino.get_pred(254), (usize::MAX, false));
        assert_eq!(pino.get_pred(258), (2, false));
        assert_eq!(pino.get_pred(510), (2, false));
        assert_eq!(pino.get_pred(511), (3, true));
        assert_eq!(pino.get_pred(512), (4, true));
    }

    #[test]
    fn test_single_bucket() {
        let values: Vec<u64> = (1000..1100).collect();

        let pino = Pred8vPino::from_sorted(&values);

        for (i, &v) in values.iter().enumerate() {
            assert_eq!(pino.get_pred(v as usize), (i, true));
        }
    }
    #[test]
    fn test_one_value_per_bucket() {
        let values: Vec<u64> = (0..20).map(|i| i * 256).collect();

        let pino = Pred8vPino::from_sorted(&values);

        for (i, &v) in values.iter().enumerate() {
            assert_eq!(pino.get_pred(v as usize), (i, true));

            if i > 0 {
                assert_eq!(pino.get_pred((v - 1) as usize), (i - 1, false));
            }
        }
    }

    #[test]
    fn test_serialize_roundtrip() {
        let values = vec![10, 20, 300, 1000, 4000];

        let pino = Pred8vPino::from_sorted(&values);

        let mut bytes = Vec::new();
        pino.serialize(&mut bytes).unwrap();

        let loaded = Pred8vPino::load(bytes.as_slice()).unwrap();

        assert_eq!(pino, loaded);
    }

    #[test]
    fn test_random_against_reference() {
        use rand::{Rng, SeedableRng};

        let mut rng = rand::rngs::StdRng::seed_from_u64(1);

        let mut values: Vec<u64> = (0..500).map(|_| rng.gen_range(0, 10_000)).collect();
        values.sort();
        values.dedup();

        let pino = Pred8vPino::from_sorted(&values);

        for key in 0..10000usize {
            let expected_rank = values.partition_point(|&x| x < key as u64);
            let expected_exact = values.binary_search(&(key as u64)).is_ok();

            assert_eq!(
                pino.rank_contains(key),
                (expected_rank, expected_exact),
                "key={}",
                key
            );
        }
    }

    #[test]
    fn test_get_words_in_range_large_start() {
        let values = vec![1000, 1010, 1030];

        let pino = Pred8vPino::from_sorted(&values);

        let words = pino.get_words_in_range(1000..1064);

        let bit = 1000 - 960; // window starts at 960

        assert_ne!(words[bit / 64] & (1u64 << (bit % 64)), 0);
    }

    #[test]
    fn test_get_words_in_range_empty() {
        let values = vec![1000, 2000];

        let pino = Pred8vPino::from_sorted(&values);

        let words = pino.get_words_in_range(0..64);

        assert!(words.iter().all(|w| *w == 0));
    }
    #[test]
    fn test_get_words_in_range_multiple_buckets() {
        let values = vec![250, 260, 511, 512, 700];

        let pino = Pred8vPino::from_sorted(&values);

        let words = pino.get_words_in_range(240..720);

        let window_start = 192;

        for value in window_start..704 {
            let bit = value - window_start;
            let expected = values.contains(&(value as u64));

            let actual = (words[bit / 64] & (1u64 << (bit % 64))) != 0;

            assert_eq!(actual, expected, "value {}", value);
        }
    }
    #[test]
    fn test_get_words_in_range_unaligned() {
        let values = vec![70, 80, 90, 130];

        let pino = Pred8vPino::from_sorted(&values);

        let words = pino.get_words_in_range(70..131);

        assert_eq!(words.len(), 2);

        let window_start = 64;

        for value in window_start..192 {
            let bit = value - window_start;
            let expected = values.contains(&(value as u64));

            let actual = (words[bit / 64] & (1u64 << (bit % 64))) != 0;

            assert_eq!(actual, expected, "value {}", value);
        }
    }

    #[test]
    fn test_window256_random_against_reference() {
        use rand::{Rng, SeedableRng};

        let mut rng = rand::rngs::StdRng::seed_from_u64(1);

        for _case in 0..200 {
            let mut values: Vec<u64> = (0..300).map(|_| rng.gen_range(0, 5000)).collect();
            values.sort();
            values.dedup();

            if values.is_empty() {
                continue;
            }

            let pino = Pred8vPino::from_sorted(&values);

            for _ in 0..100 {
                let key = rng.gen_range(0, 5200);

                let words = pino.window256_starting_from_key_word_aligned(key);

                let window_start = key & !63;

                for bit in 0..256 {
                    let value = window_start + bit;

                    let expected = values.contains(&(value as u64));

                    let actual = (words[bit / 64] & (1u64 << (bit % 64))) != 0;

                    assert_eq!(actual, expected, "key={}, value={}", key, value);
                }
            }
        }
    }

    #[test]
    fn test_get_words_in_range_random_against_reference() {
        use rand::{Rng, SeedableRng};

        let mut rng = rand::rngs::StdRng::seed_from_u64(2);

        for _case in 0..200 {
            let mut values: Vec<u64> = (0..300).map(|_| rng.gen_range(0,5000)).collect();
            values.sort();
            values.dedup();

            if values.is_empty() {
                continue;
            }

            let pino = Pred8vPino::from_sorted(&values);

            for _ in 0..100 {
                let start = rng.gen_range(0, 5000);
                let len = rng.gen_range(1, 300);

                let range = start..start + len;

                let words = pino.get_words_in_range(range.clone());

                let window_start = start & !63;

                for value in window_start..window_start + words.len() * 64 {
                    let expected = value >= range.start
                        && value < range.end
                        && values.contains(&(value as u64));

                    let bit = value - window_start;

                    let actual = (words[bit / 64] & (1u64 << (bit % 64))) != 0;

                    assert_eq!(actual, expected, "range={:?}, value={}", range, value);
                }
            }
        }
    }

    #[test]
    fn test_fill_bucket_window_all_buckets() {
        use rand::{Rng, SeedableRng};

        let mut rng = rand::rngs::StdRng::seed_from_u64(3);

        for _case in 0..100 {
            let mut values: Vec<u64> = (0..500).map(|_| rng.gen_range(0,20000)).collect();
            values.sort();
            values.dedup();

            if values.is_empty() {
                continue;
            }

            let pino = Pred8vPino::from_sorted(&values);

            for bucket in 0..pino.nblocks {
                let mut words = [0u64; 4];

                let local_begin = bucket * 256;
                let local_end = ((bucket + 1) * 256).min(pino.u + 1);

                pino.fill_bucket_window(
                    bucket,
                    local_begin,
                    local_end,
                    pino.min + local_begin,
                    &mut words,
                );

                for bit in 0..256 {
                    let global = pino.min + local_begin + bit;

                    let expected = values.binary_search(&(global as u64)).is_ok();

                    let actual = (words[bit / 64] & (1u64 << (bit % 64))) != 0;

                    assert_eq!(actual, expected, "bucket={}, value={}", bucket, global);
                }
            }
        }
    }
}

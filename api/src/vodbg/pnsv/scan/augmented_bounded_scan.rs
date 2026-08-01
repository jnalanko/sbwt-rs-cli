use super::*;

/// Performs a bounded SIMD scan to find a previous/next smaller value. If the scan fails it falls
/// back to a NND on a bitvector which marks the first and last values which are smaller than a
/// given target length in each region determined by a SIMD word.
pub struct AugmentedBoundedScan {
    pub lcs_simd: LcsSimd8x32,
    pub target_length_lower: usize,
    pub target_length_upper: usize,
    pub scan_word_bound: usize,
    pub levels: Vec<NND<16>>,
}

impl AugmentedBoundedScan {
    pub fn from_iterator<T, I>(lcs_simd: LcsSimd8x32, input: I, scan_word_bound: usize, target_length_lower: usize, target_length_upper: usize) -> Self
    where
        T: Into<usize>,
        I: Iterator<Item = T> + Clone,
    {
        let n = lcs_simd.len();
        let mut levels = Vec::with_capacity(target_length_upper - target_length_lower + 1);

        for target_length in target_length_lower..=target_length_upper {
            let select = SelectFromSimdWord::new(input.clone(), n, target_length);
            let nnd: NND<16> = NND::new(select, n);
            levels.push(nnd);
        }

        Self {
            lcs_simd,
            target_length_lower,
            target_length_upper,
            scan_word_bound,
            levels
        }
    }
}

impl Pnsv for AugmentedBoundedScan {
    fn previous(&self, index: usize, target_length: usize) -> usize {
        if target_length > self.target_length_upper {
            return self.lcs_simd.scan_left(index, target_length);
        }
        let result = self.lcs_simd.scan_left_bounded(index, target_length, self.scan_word_bound);
        match result {
            Ok(index) => index,
            Err(continue_search_index) => {
                let nnd_index = target_length - self.target_length_lower;
                if nnd_index < self.levels.len() {
                    self.levels[nnd_index].previous(continue_search_index)
                } else {
                    index
                }
            }
        }
    }

    fn next(&self, index: usize, target_length: usize) -> usize {
        if target_length > self.target_length_upper {
            return self.lcs_simd.scan_right(index, target_length);
        }
        let result = self.lcs_simd.scan_right_bounded(index, target_length, self.scan_word_bound);
        match result {
            Ok(index) => index,
            Err(continue_search_index) => {
                let nnd_index = target_length - self.target_length_lower;
                if nnd_index < self.levels.len() {
                    self.levels[nnd_index].next(continue_search_index)
                } else {
                    index
                }
            }
        }
    }
    
    #[inline]
    fn max_target(&self) -> usize {
        self.target_length_upper
    }
}

// An iterator over the indices of the first and last value in each simd word smaller than a given
// target length.
struct SelectFromSimdWord<T, I>
where 
    T: Into<usize>,
    I: Iterator<Item = T> + Clone
{
    target_length: usize,
    index: usize,
    previous_index_of_smaller: usize,
    last_yielded: usize,
    count: usize,
    pending: Option<usize>,
    underlying_values: I,
}

// note(mk): The purpose of this iterator is only to construct a Nearest Neighbour Dictionary out
// of the first and last values which are smaller than the given target length in each SIMD word.
// Since the SIMD scan always ends on a border between two words, then the answer to a query in the
// NND will always be the first or last smaller value. Therefore, there is no need to store any
// 1-bits inbetween.
impl<T, I> SelectFromSimdWord<T, I> 
where 
    T: Into<usize>,
    I: Iterator<Item = T> + Clone
{
    fn new(underlying_values: I, count: usize, target_length: usize) -> Self {
        Self {
            target_length,
            index: 0,
            previous_index_of_smaller: count,
            last_yielded: count + 1,
            count,
            pending: None,
            underlying_values,
        }
    }
}

impl<T, I> Iterator for SelectFromSimdWord<T, I>
where 
    T: Into<usize>,
    I: Iterator<Item = T> + Clone
{
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        if self.pending.is_some() {
            self.last_yielded = self.pending.unwrap();
            return self.pending.take();
        }

        if self.index == self.count {
            return None;
        }

        #[allow(clippy::while_let_on_iterator)]
        while let Some(value) = self.underlying_values.next() {
            let value = value.into();

            let index = self.index;
            self.index += 1;

            if value >= self.target_length {
                continue;
            }

            let previous_index_of_smaller = self.previous_index_of_smaller;
            self.previous_index_of_smaller = index;

            if previous_index_of_smaller < self.count {
                let previous_word = previous_index_of_smaller / LcsSimd8x32::LANES;
                let current_word = index / LcsSimd8x32::LANES;

                if previous_word != current_word {
                    if self.last_yielded == previous_index_of_smaller {
                        self.last_yielded = index;
                        return Some(index);
                    } else {
                        self.pending = Some(index);
                        // note(mk): This next line is probably not necessary.
                        self.last_yielded = previous_index_of_smaller;
                        return Some(previous_index_of_smaller);
                    }
                }
            } else {
                return Some(index);
            }
        }

        if self.previous_index_of_smaller < self.count && self.previous_index_of_smaller != self.last_yielded {
            return Some(self.previous_index_of_smaller);
        }

        None
    }
}

// note(mk): There were some errors with the automatically derived trait. 
impl<T, I> Clone for SelectFromSimdWord<T, I>
where 
    T: Into<usize>,
    I: Iterator<Item = T> + Clone
{
    fn clone(&self) -> Self {
        Self {
            target_length: self.target_length,
            index: self.index,
            previous_index_of_smaller: self.previous_index_of_smaller,
            last_yielded: self.last_yielded,
            count: self.count,
            pending: self.pending,
            underlying_values: self.underlying_values.clone()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn select_from_simd_word() {
        let items: &[u8] = &[
            4, 3, 2, 1, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 1, 2, 3, 4,

            4, 3, 2, 1, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 1, 2, 3, 4,

            4, 3, 2, 1, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 1, 2, 3, 4,

            4, 3, 2, 1, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 1, 2, 3, 4,

            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,

            1, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,

            1, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
        ];

        let less_than_5: Vec<usize> = SelectFromSimdWord::new(items.iter().cloned(), items.len(), 5).collect();
        let less_than_4: Vec<usize> = SelectFromSimdWord::new(items.iter().cloned(), items.len(), 4).collect();
        let less_than_3: Vec<usize> = SelectFromSimdWord::new(items.iter().cloned(), items.len(), 3).collect();
        let less_than_2: Vec<usize> = SelectFromSimdWord::new(items.iter().cloned(), items.len(), 2).collect();

        assert_eq!(less_than_5, &[0, 31, 32, 63, 64, 95, 96, 127, 160, 192]);
        assert_eq!(less_than_4, &[1, 30, 33, 62, 65, 94, 97, 126, 160, 192]);
        assert_eq!(less_than_3, &[2, 29, 34, 61, 66, 93, 98, 125, 160, 192]);
        assert_eq!(less_than_2, &[3, 28, 35, 60, 67, 92, 99, 124, 160, 192]);
    }

    #[test]
    fn bounded_scan_with_fallback() {
        let items: &[u8] = &[
            4, 3, 2, 1, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 1, 2, 3, 4,

            4, 3, 2, 1, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 1, 2, 3, 4,

            4, 3, 2, 1, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 1, 2, 3, 4,

            4, 3, 2, 1, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 1, 2, 3, 4,

            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,

            1, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,

            1, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
        ];

        let target_length_lower = 2;
        let target_length_upper = 5;

        let lcs_simd = LcsSimd8x32::from_iterator(items.iter().cloned(), items.len(), 31);
        let abs = AugmentedBoundedScan::from_iterator(lcs_simd.clone(), items.iter().cloned(), 0, target_length_lower, target_length_upper);

        for target_length in target_length_lower..=target_length_upper+1 {
            for i in 0..items.len() {
                assert_eq!(
                    abs.previous(i, target_length),
                    lcs_simd.previous(i, target_length),
                    "previous; i: {}; target_length: {}",
                    i, target_length
                );
                assert_eq!(
                    abs.next(i, target_length),
                    lcs_simd.next(i, target_length),
                    "next; i: {}; target_length: {}",
                    i, target_length
                );
            }
        }
    }
}

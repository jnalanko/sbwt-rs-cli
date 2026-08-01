// Code by Martin Kostadinov adapted to Rust from the C++ library sdsl-lite.
use nearest_neighbor_dictionary::NearestNeighbourDictionary as NND;

use crate::LcsArray;
use super::Pnsv;

use simple_sds_sbwt::{
    bit_vector::BitVector,
    ops::{BitVec, Rank, Select},
    raw_vector::{AccessRaw, RawVector}
};

pub mod nearest_neighbor_dictionary;
pub mod lookup;
pub mod util;

/// A Balanced Parenthesis structure with support for the operations, rank, select, find open and
/// enclose.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Bp {
    pub levels: Vec<BpLevel>,
    pub block_size: usize,
}

/// One level of the recursive Balanced Parenthesis structure.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BpLevel {
    /// The original sequence of balanced parenthesis.
    pub parenthesis: BitVector,

    /// The extracted Nearest Neighbour Dictionary from the pioneers bit mask i.e. the bitvector
    /// with 1s in positions where a pioneer (and its match(?)) were located.
    pub pioneer_nnd: NND,
}

impl Bp {
    pub fn new(mut parenthesis_raw: RawVector, block_size: usize) -> Self {
        let mut levels: Vec<BpLevel> = vec![];
        let mut pioneers = RawVector::with_len(parenthesis_raw.len(), false);
        let mut parenthesis_stack = vec![];
        let mut previous_pioneer_count = parenthesis_raw.len();
        let mut last_level_index = 0;
        let mut len;

        loop {
            len = parenthesis_raw.len();

            let mut parenthesis = BitVector::from(parenthesis_raw);
            parenthesis.enable_rank();
            parenthesis.enable_select();

            let pioneer_count = calculate_pioneers_bitmask(&parenthesis, block_size, &mut pioneers, &mut parenthesis_stack);

            let pioneer_nnd: NND<8> = if pioneer_count > 0 {
                let bv = BitVector::from(pioneers);
                let result = NND::new(util::one_indices_iterator(bv.iter()), bv.len());
                pioneers = bv.into();
                result
            } else {
                NND::empty()
            };
            
            levels.push(BpLevel {
                parenthesis,
                pioneer_nnd,
            });

            if pioneer_count == previous_pioneer_count {
                panic!("The number of pioneers doesn't converge.");
            }

            if len <= block_size || pioneer_count == 0 {
                break;
            }

            parenthesis_raw = RawVector::with_len(pioneer_count, false);
            for i in 0..levels[last_level_index].pioneer_nnd.one_count {
                let index_in_parenthesis_vector = levels[last_level_index].pioneer_nnd.select(i);
                let parenthesis_value = levels[last_level_index].parenthesis.get(index_in_parenthesis_vector);
                parenthesis_raw.set_bit(i, parenthesis_value);
            }
            last_level_index += 1;
            previous_pioneer_count = pioneer_count;
        }

        Self {
            levels,
            block_size,
        }
    }

    #[inline]
    pub fn rank(&self, parenthesis_index: usize) -> usize {
        self.levels[0].parenthesis.rank(parenthesis_index)
    }

    #[inline]
    pub fn select(&self, item_index: usize) -> usize {
        self.levels[0].parenthesis.select(item_index).unwrap()
    }

    /// The left excess of the balanced parenthesis sequence up to and including the bracket under
    /// the given index.
    #[inline]
    fn excess(&self, index: usize, level: usize) -> usize {
        debug_assert!(!self.levels.is_empty());
        let open_count = self.levels[level].parenthesis.rank(index + 1);
        // Equivalent to:
        // open_count - close_count =
        //  = open_count - ((index + 1) - open_count)
        open_count * 2 - (index + 1)
    }

    #[inline]
    pub fn find_open(&self, index: usize) -> usize {
        self.find_open_in_level(index, 0)
    }

    fn find_open_in_level(&self, index: usize, level: usize) -> usize {
        debug_assert!(level < self.levels.len());
        debug_assert!(index < self.levels[level].parenthesis.len());
        if self.levels[level].parenthesis.get(index) {
            return index;
        }
        let balanced_parenthesis =  self.levels[level].parenthesis.as_ref().get_words();
        let potentially_matching = lookup::try_to_find_opening_in_block(balanced_parenthesis, index - 1, 1, self.block_size);
        if let Some(index) = potentially_matching {
            return index;
        }

        //
        // Gives the 0-based index of a parenthesis in the pioneer family. According to the paper
        // this is guaranteed to be a closing pioneer. Now, the match of the query parenthesis must
        // be in the same block as the match of the found closing pioneer, otherwise the query
        // parenthesis must have been itself a closing pioneer. If the query parenthesis is a
        // closing pioneer, then its match will be found either way by the latter part of this
        // procedure.
        //
        let closing_parenthesis_pioneer_family_member = self.levels[level].pioneer_nnd.rank(index);
        let matching_pioneer_from_next_level = self.find_open_in_level(closing_parenthesis_pioneer_family_member, level + 1);
        let position_of_matching_pioneer_in_current_level = self.levels[level].pioneer_nnd.select(matching_pioneer_from_next_level);
        let block_of_searched_open_parenthesis = (position_of_matching_pioneer_in_current_level / self.block_size) * self.block_size;

        // 
        // q - query parenthesis
        // m - matching parenthesis
        // s - starting parenthesis for the search in the located block
        //
        //     v block of the match to the query parenthesis,
        //     v found by locating the pioneers
        //     v 
        //     v                  v block of the query parenthesis
        //     v                  v
        // ... [..m........s] ... [.......q....] ...
        //
        //
        // The excess at the matching parenthesis is one more than the excess at the query
        // parenthesis. This is because the string of parenthesis between the query parenthesis and
        // its match must be balanced i.e. their excess must total 0. In addition, the excess at
        // the 's' parenthesis must be greater than the excess at the query parenthesis. The
        // difference between the excess at 's' and the escess at 'm' gives exactly the number of
        // opening parenthesis that have to be scanned in the found block by the helper procedure
        // in order to find the opening parenthesis of q.
        //

        let index_of_last_parenthesis_in_block = block_of_searched_open_parenthesis + self.block_size - 1;
        // let excess_of_first_parenthesis_in_next_block = self.excess(index_of_last_parenthesis_in_block + 1, level) as isize;
        let excess_of_last_parenthesis_in_block = self.excess(index_of_last_parenthesis_in_block, level) as isize;
        let excess_of_query_parenthesis = self.excess(index, level) as isize;
        // let excess_value_of_parenthesis_from_next_block: isize = 1 -
        //     2 * if self.levels[level].parenthesis.get(index_of_last_parenthesis_in_block + 1) { 1 } else { 0 };
        // let relative_excess = excess_of_first_parenthesis_in_next_block - excess_of_query_parenthesis + excess_value_of_parenthesis_from_next_block;

        //
        // Look at the relative excess calculation in the enclose procedure for an explanation why
        // there is an additional +1 in the beginning. The excess of the parenthesis that is being
        // searched is one more than the excess of the query parenthesis. And so the relative
        // excess is:
        // 
        // 1 + excess_of_last_parenthesis_in_block - (excess_of_query_parenthesis + 1)
        //
        // which is equivalent to this:
        //
        let relative_excess = excess_of_last_parenthesis_in_block - excess_of_query_parenthesis;
        let result = lookup::try_to_find_opening_in_block(
            balanced_parenthesis,
            index_of_last_parenthesis_in_block,
            relative_excess as usize,
            self.block_size
        );

        result.expect("The matching opening parenthesis should be in this block.")
    }

    #[inline]
    pub fn enclose(&self, index: usize) -> usize {
        self.enclose_in_level(index, 0)
    }

    fn enclose_in_level(&self, index: usize, level: usize) -> usize {
        if !self.levels[level].parenthesis.get(index) {
            // This is a closing parenthesis.
            return self.find_open_in_level(index, level);
        }
        let excess = self.excess(index, level) as isize;
        if excess == 1 {
            return index;
        }

        let balanced_parenthesis =  self.levels[level].parenthesis.as_ref().get_words();
        let potentially_enclosing = lookup::try_to_find_opening_in_block(balanced_parenthesis, index - 1, 1, self.block_size);
        if let Some(result) = potentially_enclosing {
            return result;
        }

        let next_pioneer_family_member_index = self.levels[level].pioneer_nnd.rank(index);
        let enclosing_for_pioneer_in_next_level = self.enclose_in_level(next_pioneer_family_member_index, level + 1);
        let position_of_pioneer_enclosing_in_this_level = self.levels[level].pioneer_nnd.select(enclosing_for_pioneer_in_next_level);
        let block_of_pioneer_enclosing = (position_of_pioneer_enclosing_in_this_level / self.block_size) * self.block_size;

        //
        // Works on the same principle as the find_open procedure. This time, however, we are
        // searching for a parenthesis whose excess is one less compared to this parenthesis'
        // excess.
        //

        let last_parenthesis_in_block = block_of_pioneer_enclosing + self.block_size - 1;
        let excess_of_last_parenthesis_in_block = self.excess(last_parenthesis_in_block, level) as isize;
        // 
        // To get the number of opening parenthesis that need to be scanned such that the last one
        // scanned is the enclosing parenthesis an additional 1 must be added. The excess of the
        // parenthesis that is being searched for is one less than the excess of the query opening
        // parenthesis. And so the relative excess is:
        //
        // 1 + excess_of_last_parenthesis_in_block - (excess - 1)
        //
        // which is equivalent to:
        //
        let relative_excess = excess_of_last_parenthesis_in_block - excess + 2;
        let result = lookup::try_to_find_opening_in_block(balanced_parenthesis, last_parenthesis_in_block, relative_excess as usize, self.block_size);

        result.expect("The enclosing parenthesis should be in this block.")
    }
}

fn calculate_pioneers_bitmask(
    balanced_parenthesis: &BitVector,
    block_size: usize, 
    pioneers: &mut RawVector,
    parenthesis_stack: &mut Vec<usize>,
) -> usize {
    pioneers.resize(balanced_parenthesis.len(), false);

    // note(mk): There must be a faster way to do this...
    for i in 0..balanced_parenthesis.len() {
        pioneers.set_bit(i, false);
    }

    // The maximum size of the parenthesis stack is k (the length of the k-mers).
    parenthesis_stack.clear();

    let mut pioneer_count = 0;

    let mut current_pioneer_block = 0;
    let mut previous_open_pioneer_index = 0;
    let mut previous_pioneer_match_index = 0;
    let mut first_index_in_block = 0;

    let mut bits_left_in_block = 0;

    let mut index = 0;
    let mut is_open;
    while index < balanced_parenthesis.len() {
        is_open = balanced_parenthesis.get(index);

        if bits_left_in_block == 0 {
            // Set the current pioneer block as each open far parenthesis in its block is a match
            // to a closing pioneer.
            current_pioneer_block = index / block_size;
            first_index_in_block = index;
            bits_left_in_block = block_size;
        }

        if is_open {
            if bits_left_in_block > 1 {
                // Check and skip for an open and closed parenthesis which are immediately next to
                // each other in the same block.
                let next_is_open = balanced_parenthesis.get(index + 1);
                if !next_is_open {
                    index += 2;
                    bits_left_in_block -= 2;
                    continue;
                }
            }
            parenthesis_stack.push(index);
        } else {
            let matching_parenthesis_index = parenthesis_stack.pop()
                .expect("The parenthesis sequence should be balanced.");
            if matching_parenthesis_index < first_index_in_block {
                let new_pioneer_block = matching_parenthesis_index / block_size;
                if new_pioneer_block == current_pioneer_block {
                    // The previous parenthesis that were marked as pioneers turned out not to be
                    // such. Reset their values in the bitmask.
                    pioneers.set_bit(previous_open_pioneer_index, false);
                    pioneers.set_bit(previous_pioneer_match_index, false);
                    pioneer_count -= 2;
                }

                pioneers.set_bit(matching_parenthesis_index, true);
                pioneers.set_bit(index, true);
                current_pioneer_block = new_pioneer_block;
                previous_open_pioneer_index = matching_parenthesis_index;
                previous_pioneer_match_index = index;

                pioneer_count += 2;
            }
        }

        index += 1;
        bits_left_in_block -= 1;
    }
    debug_assert!(parenthesis_stack.is_empty());

    // println!("============================================================================================");
    // for i in 0..balanced_parenthesis.len() {
    //     print!("{}", if pioneers.bit(i) { 1 } else { 0 });
    // }
    // println!();

    pioneer_count
}

/// Previous and Next Smaller Value using Balanced Parenthesis.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PnsvBp {
    pub previous: Bp,
    pub next: Bp,
    pub count: usize,
}

impl PnsvBp {
    /// Constructs a data structure which supports Previous and Next Smaller Value queries on the
    /// values in the original array.
    pub fn from_iterator<T, I>(values: I, count: usize, block_size: usize) -> Self
    where
        T: Ord + Copy,
        I: Iterator<Item = T> + DoubleEndedIterator + Clone,
    {
        let mut stack = vec![];
        let previous_smaller_bp_vector = Self::make_psv_bp_vector(values.clone(), count, &mut stack);
        let next_smaller_bp_vector = Self::make_psv_bp_vector(values.rev(), count, &mut stack);

        let previous = Bp::new(previous_smaller_bp_vector, block_size);
        let next = Bp::new(next_smaller_bp_vector, block_size);

        Self { previous, next, count }
    }

    /// Given an index in the original array finds the index of the previous smaller value.
    #[inline]
    pub fn previous(&self, index: usize) -> usize {
        let parenthesis_index = self.previous.select(index);
        let parent_parenthesis_index = self.previous.enclose(parenthesis_index);
        self.previous.rank(parent_parenthesis_index)
    }

    /// Given an index in the original array finds the index of the next smaller value.
    #[inline]
    pub fn next(&self, index: usize) -> usize {
        let reverse_index = self.count - 1 - index;
        let parenthesis_index = self.next.select(reverse_index);
        let parent_parenthesis_index = self.next.enclose(parenthesis_index);
        let reverse_result_index = self.next.rank(parent_parenthesis_index);
        self.count - 1 - reverse_result_index
    }

    /// Given an iterator of values constructs the balanced parenthesis representation of the
    /// Previous Smaller Value tree. A PSV tree is such that each vertex's parent is the previous
    /// smaller value in the array.
    pub fn make_psv_bp_vector<T, I>(values: I, count: usize, stack: &mut Vec<T>) -> RawVector
    where
        T: Ord + Copy,
        I: Iterator<Item = T>,
    {
        // The maximum height of the stack is the number of unique values that the array contains.
        // In the case of the Longest Common Suffix array that would be k.
        stack.clear();

        let mut result = RawVector::with_len(count * 2, false);
        let mut bit_iterator = 0;

        for value in values {
            loop {
                if stack.is_empty() {
                    break;
                }

                let previous_value = *stack.last().unwrap();
                if previous_value < value {
                    // The previous value is smaller, therefore it is this value's parent.
                    break;
                }

                // The previous value is either greater than this value in which case it is a sibling
                // or a child of a sibling of this value, or it is equal to this value in which case it
                // is a sibling of this value in the PSV tree.
                stack.pop();

                // By moving the iterator forward we leave a 0 in the resulting bitvector which is
                // equivalent to closing the parenthesis.
                bit_iterator += 1;
            }

            result.set_bit(bit_iterator, true);
            bit_iterator += 1;
            stack.push(value);
        }

        result
    }
}

pub struct LcsPnsvBp<'a> {
    pub lcs: &'a LcsArray,
    pub pnsv: PnsvBp,
}

impl<'a> LcsPnsvBp<'a> {
    pub fn new(lcs: &'a LcsArray, block_size: usize) -> Self {
        let iterator = (0..lcs.len()).map(|index| lcs.access(index));
        let pnsv = PnsvBp::from_iterator(iterator, lcs.len(), block_size);
        Self {
            lcs,
            pnsv,
        }
    }
}

impl<'a> Pnsv for LcsPnsvBp<'a> {
    fn previous(&self, mut index: usize, target_length: usize) -> usize {
        let mut jump;
        while 0 < index && self.lcs.access(index) >= target_length {
            jump = self.pnsv.previous(index);
            if index == jump {
                index = 0;
                break;
            }
            index = jump;
        }
        index
    }

    fn next(&self, mut index: usize, target_length: usize) -> usize {
        let mut jump;
        while index < self.lcs.len() && self.lcs.access(index) >= target_length {
            jump = self.pnsv.next(index);
            if index == jump {
                index = self.lcs.len();
                break;
            }
            index = jump;
        }
        index
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn find_open_simple() {
        let count = 500;
        let middle = count / 2;
        let mut balanced_parenthesis = RawVector::with_len(count, false);
        for i in 0..middle {
            balanced_parenthesis.set_bit(i, true);
        }
        let bp = Bp::new(balanced_parenthesis, 17);
        for i in middle..count {
            println!("'(': {} ')': {}", count - 1 - i, i);
            assert_eq!(count - 1 - i, bp.find_open(i));
        }
    }

    #[test]
    fn find_open_without_super_root() {
        let count = 256;
        let mut balanced_parenthesis = RawVector::with_len(count, false);
        for i in 0..count / 2 { 
            balanced_parenthesis.set_bit(2 * i, true);
        }
        let bp = Bp::new(balanced_parenthesis, 20);
        for i in 0..count / 2 { 
            println!("'(': {} ')': {}", 2 * i, 2 * i + 1);
            assert_eq!(2 * i, bp.find_open(2 * i + 1));
        }
    }

    fn bp_vector_from_string(string: &str) -> RawVector {
        let mut vector = RawVector::with_len(string.len(), false);
        for (index, c) in string.chars().enumerate() {
            if c == '(' {
                vector.set_bit(index, true);
            }
        }
        vector
    }

    #[test]
    fn specific_enclose_01() {
        let parenthesis = "()(()()(()()(((()))()((())))))";
        let balanced_parenthesis = bp_vector_from_string(parenthesis);
        let bp = Bp::new(balanced_parenthesis, 20);
        assert_eq!(7, bp.rank(bp.enclose(bp.select(12))));
    }

    #[test]
    fn specific_enclose_02() {
        let parenthesis = "(((())))";
        let balanced_parenthesis = bp_vector_from_string(parenthesis);
        let bp = Bp::new(balanced_parenthesis, 3);
        assert_eq!(2, bp.enclose(3));
    }

    #[test]
    fn enclose_simple() {
        let count = 500;
        let middle = count / 2;
        let mut balanced_parenthesis = RawVector::with_len(count, false);
        for i in 0..middle {
            balanced_parenthesis.set_bit(i, true);
        }
        let bp = Bp::new(balanced_parenthesis, 17);
        for i in 1..middle {
            println!("'(': {} '((': {}", i-1, i);
            assert_eq!(i-1, bp.enclose(i));
        }
    }

    fn make_psv_tree<T, I>(values: I, count: usize) -> Vec<usize>
    where
        T: Ord + Copy,
        I: Iterator<Item = T>,
    {
        let mut stack: Vec<(usize, T)> = vec![];

        let mut result = vec![count; count];

        for (index, value) in values.enumerate() {
            loop {
                if stack.is_empty() {
                    break;
                }
                let (previous_index, previous_value) = *stack.last().unwrap();
                if previous_value < value {
                    result[index] = previous_index;
                    break;
                }
                stack.pop();
            }

            stack.push((index, value));
        }

        result
    }

    #[test]
    fn bp_method_random_sequence() {
        use rand::distributions::uniform::Uniform;
        use rand::distributions::Distribution;
        use rand::rngs::StdRng;
        use rand::SeedableRng;

        let n = 50000;
        let k = 31u8;
        let mut rng = StdRng::from_seed([42; 32]);
        let uniform = Uniform::from(0..k);
        let numbers = uniform.sample_iter(&mut rng).take(n).collect::<Vec<_>>();

        let answers = make_psv_tree(numbers.iter(), numbers.len());

        let mut bp_value_stack = vec![];
        let bp_vector = PnsvBp::make_psv_bp_vector(numbers.iter(), numbers.len(), &mut bp_value_stack);

        for block_size in (10..=100).step_by(10) {
            let psv_tree = Bp::new(bp_vector.clone(), block_size);

            for i in 0..n {
                let parenthesis_index = psv_tree.select(i);
                let parenthesis_index_of_parent = psv_tree.enclose(parenthesis_index);
                if answers[i] == numbers.len() {
                    assert!(
                        parenthesis_index_of_parent == parenthesis_index,
                        "Failed at i: {}, parenthesis_index: {}",
                        i,
                        parenthesis_index
                    );
                    println!("{i} -> null");
                    continue;
                }
                let parent_index = psv_tree.rank(parenthesis_index_of_parent);
                assert_eq!(
                    answers[i], parent_index,
                    "Failed at i: {}, parenthesis_index: {}",
                    i, parenthesis_index
                );
                println!("{i} -> {}", parent_index);
            }
        }
    }
}


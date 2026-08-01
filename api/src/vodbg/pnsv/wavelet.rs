// Code by Martin Kostadinov.

use simple_sds_sbwt::bit_vector::BitVector;
use simple_sds_sbwt::ops::{BitVec, Rank, Select, SelectZero};
use simple_sds_sbwt::raw_vector::{PushRaw, RawVector};
use simple_sds_sbwt::serialize::Serialize;

use std::collections::vec_deque::VecDeque;

/// A wavelet tree over a continuous subset of the whole (ordered) alphabet. Supports previous
/// smaller and next smaller value operations in log(n) time where n is the size of the subset.
/// The descriptor "bounded" in the name of this data structure was not chose to avoid ambiguity
/// with the already existing BWT abbreviation.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct WindowedWaveletTree {
    pub lower_bound: usize,
    pub window_size: usize,
    pub tree: Vec<Node>,
    pub data: Vec<BitVector>,
}

/// A node in the wavelet tree. Contains indices for navigating the tree.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Node {
    pub parent: usize,
    /// The smallest character in the range this node is responsible for.
    pub lower_bound: usize,
    /// The size of the subset of the alphabet this node is responsible for.
    pub window_size: usize,
    pub left_child: usize,
    pub right_child: usize,
    pub data: usize,
}

impl Node {
    pub const NULL: usize = 0;
}

impl WindowedWaveletTree {
    pub const STACK_SIZE: usize = 32;

    pub fn empty() -> Self {
        Self {
            lower_bound: 0,
            window_size: 0,
            tree: vec![],
            data: vec![],
        }
    }

    pub fn from_iterator<T, I>(input: I, count: usize, lower_bound: usize, window_size: usize) -> Self
    where
        T: Into<usize>,
        I: Iterator<Item = T> + Clone,
    {
        let mut tree: Vec<Node> = vec![];
        Self::make_tree_nodes(&mut tree, lower_bound, window_size);

        let ten_percent = count / 10;
        let mut border = ten_percent;
        let mut percent_count = 0;

        // First pass over the items in order to get the number of bits per bitvector. When
        // allocating the bitvectors this count will be used to allocate enough space so that extra
        // allocations from the expansion of the bitvector when the bits are pushed one by one are
        // not needed.
        let upper_bound = lower_bound + window_size - 1;
        let mut bits_per_node = vec![0_usize; tree.len()];
        let mut current_node;
        for (index, item) in input.clone().enumerate() {
            let clamped_item = item.into().clamp(lower_bound, upper_bound);

            current_node = 0;

            loop {
                bits_per_node[current_node] += 1;
                let left_child_window_size =
                    (tree[current_node].window_size & 1) + (tree[current_node].window_size >> 1);

                if clamped_item < tree[current_node].lower_bound + left_child_window_size {
                    // Go to the left child of the current node.
                    current_node = tree[current_node].left_child;
                    if current_node == Node::NULL {
                        break;
                    }
                    continue;
                }

                // Go to the right child of the current node.
                current_node = tree[current_node].right_child;
                if current_node == Node::NULL {
                    break;
                }
            }

            if index >= border {
                border += ten_percent;
                percent_count += 1;
                log::info!("[WWT::from_iterator] first pass scanning... {}0%", percent_count);
            }
        }

        border = ten_percent;
        percent_count = 0;

        // Allocate the raw vectors with the needed capacity.
        let mut data_raw: Vec<RawVector> = vec![];
        for bit_count in bits_per_node {
            data_raw.push(RawVector::with_capacity(bit_count));
        }

        // Second pass to populate the bits in the bitvectors.
        for (index, item) in input.clone().enumerate() {
            let clamped_item = item.into().clamp(lower_bound, upper_bound);

            current_node = 0;

            loop {
                let left_child_window_size =
                    (tree[current_node].window_size & 1) + (tree[current_node].window_size >> 1);
                let data_index = tree[current_node].data;
                if clamped_item < tree[current_node].lower_bound + left_child_window_size {
                    data_raw[data_index].push_bit(false);
                    current_node = tree[current_node].left_child;
                    if current_node == Node::NULL {
                        break;
                    }
                    continue;
                }

                data_raw[data_index].push_bit(true);
                current_node = tree[current_node].right_child;
                if current_node == Node::NULL {
                    break;
                }
            }

            if index >= border {
                border += ten_percent;
                percent_count += 1;
                log::info!("[WWT::from_iterator] second pass scanning... {}0%", percent_count);
            }
        }

        // Transform the raw vectors into bitvectors and enable rank and select support.
        let mut data: Vec<BitVector> = Vec::with_capacity(data_raw.len());
        for raw_vector in data_raw {
            let mut bit_vector: BitVector = raw_vector.into();
            bit_vector.enable_rank();
            bit_vector.enable_select();
            bit_vector.enable_select_zero();
            data.push(bit_vector);
        }

        Self {
            lower_bound,
            window_size,
            tree,
            data,
        }
    }

    fn make_tree_nodes(tree: &mut Vec<Node>, lower_bound: usize, window_size: usize) {
        let mut queue: VecDeque<(usize, usize, usize)> = VecDeque::new();
        let mut data_index: usize = 0;

        queue.push_back((0, lower_bound, window_size));

        // Create the nodes based on the lower bound and the window size.
        while !queue.is_empty() {
            let current_node_index = tree.len();
            let (parent, lower_bound, window_size) = queue.pop_front().unwrap();

            let right_child_window_size = window_size >> 1;
            let left_child_window_size = window_size - right_child_window_size;

            let mut left_child = 0;
            let mut right_child = 0;

            // Add 1 for the node that will be pushed now.
            let mut last_node_index = tree.len() + queue.len() + 1;
            if left_child_window_size > 1 {
                left_child = last_node_index;
                last_node_index += 1;
                queue.push_back((current_node_index, lower_bound, left_child_window_size));
            }

            if right_child_window_size > 1 {
                right_child = last_node_index;
                queue.push_back((
                    current_node_index,
                    lower_bound + left_child_window_size,
                    right_child_window_size,
                ));
            }

            let data = data_index;
            data_index += 1;

            tree.push(Node {
                parent,
                lower_bound,
                window_size,
                left_child,
                right_child,
                data,
            });
        }
    }

    /// Find the previous smaller value. If the value at the index is smaller than the target
    /// length, returns that index.
    pub fn previous(&self, mut index: usize, target_length: usize) -> usize {
        let mut current_node = 0;
        let mut candidate_stack: [(usize, usize); Self::STACK_SIZE] = [(0, 0); Self::STACK_SIZE];
        let mut candidate_count = 0;

        loop {
            if self.tree[current_node].lower_bound >= target_length {
                break;
            }

            let max_value = self.tree[current_node].lower_bound + self.tree[current_node].window_size - 1;
            if max_value < target_length {
                debug_assert!(candidate_count < Self::STACK_SIZE);
                candidate_stack[candidate_count] = (current_node, index);
                candidate_count += 1;
                break;
            }

            let left_child_window_size = (self.tree[current_node].window_size & 1) + (self.tree[current_node].window_size >> 1);
            let min_value_in_right_child = self.tree[current_node].lower_bound + left_child_window_size;

            let data_index = self.tree[current_node].data;
            let is_one = self.data[data_index].get(index);
            let one_rank = self.data[data_index].rank(index);
            let zero_rank = index - one_rank;

            if min_value_in_right_child - 1 < target_length {
                // The maximum possible value in the left child is smaller than the target length.
                // Store the absolute position of the rightmost 0 up to the current index
                // (inclusive).
                if !is_one {
                    debug_assert!(candidate_count < Self::STACK_SIZE);
                    candidate_stack[candidate_count] = (current_node, index);
                    candidate_count += 1;
                    break;
                }
                if zero_rank > 0 {
                    // The number of zeroes before can be at most the number of zeroes in the given
                    // bitvector. Searching for the last one should always succeed.
                    let rightmost_zero_index_before = self.data[data_index]
                        .select_zero(zero_rank - 1)
                        .unwrap();
                    debug_assert!(candidate_count < Self::STACK_SIZE);
                    candidate_stack[candidate_count] = (current_node, rightmost_zero_index_before);
                    candidate_count += 1;
                }
            }

            if min_value_in_right_child >= target_length {
                if min_value_in_right_child - 1 < target_length {
                    // The candidates in the left child have already been checked.
                    break;
                }
                // Move to the left child as there are no possible smaller values in the right
                // child. Update the index so that it is valid in the left child.
                index = if !is_one {
                    // If this bit is a 0, move to its index in the left child.
                    zero_rank
                } else if zero_rank > 0 {
                    // If this bit is not a 0, but there are 0 bits before it, move to the
                    // corresponding index of the rightmost 0 in the left child.
                    zero_rank - 1
                } else {
                    // Otherwise, there are no more 0 bits before the current index, therefore, the
                    // left child will not contain any candidates for the result and we can
                    // exit.
                    break;
                };
                // If the maximum value in the left child is greater than or equal to the target
                // length and the minimum value in the left child which is equal to the minimum
                // value in this node is smaller than the target length, then there are more than 2
                // values in the range of the left child i.e. it is guaranteed to exist.
                current_node = self.tree[current_node].left_child;
                continue;
            }

            // min_value_in_right_child < target_length

            // The candidate indices for the previous smaller value which are in the left child
            // have been handled. Therefore, we should try moving to the right child as there might
            // still be candidates there.

            index = if is_one {
                one_rank
            } else if one_rank > 0 {
                one_rank - 1
            } else {
                break;
            };

            // If the minimum value in the right child is smaller than the target length and the
            // maximum value in the right child which is equal to the maximum value in the range of
            // this node is greater than or equal to the target length, then there must be at least
            // 2 different values in the range of the right child i.e. it is guaranteed to exist.
            current_node = self.tree[current_node].right_child;
        }

        if candidate_count == 0 {
            return 0;
        }

        candidate_count -= 1;
        let (mut result_node, mut result) = candidate_stack[candidate_count];
        while candidate_count > 0 {
            candidate_count -= 1;
            let (current_node, index_in_node) = candidate_stack[candidate_count];
            result = self.climb(result_node, result, current_node);
            result_node = current_node;
            result = result.max(index_in_node);
        }

        self.climb(result_node, result, 0)
    }

    /// Find the next smaller value. If the value at the index is smaller than the target length,
    /// returns that index.
    pub fn next(&self, mut index: usize, target_length: usize) -> usize {
        let mut current_node = 0;
        let mut candidate_stack: [(usize, usize); Self::STACK_SIZE] = [(0, 0); Self::STACK_SIZE];
        let mut candidate_count = 0;

        loop {
            if self.tree[current_node].lower_bound >= target_length {
                break;
            }

            let max_value = self.tree[current_node].lower_bound + self.tree[current_node].window_size - 1;
            if max_value < target_length {
                debug_assert!(candidate_count < Self::STACK_SIZE);
                candidate_stack[candidate_count] = (current_node, index);
                candidate_count += 1;
                break;
            }

            let left_child_window_size = (self.tree[current_node].window_size & 1) + (self.tree[current_node].window_size >> 1);
            let min_value_in_right_child = self.tree[current_node].lower_bound + left_child_window_size;

            let data_index = self.tree[current_node].data;
            let is_one = self.data[data_index].get(index);
            let one_rank = self.data[data_index].rank(index);
            let zero_rank = index - one_rank;

            if min_value_in_right_child - 1 < target_length {
                // The maximum possible value in the left child is smaller than the target length.
                // Store the absolute position of the leftmost 0 after the current index
                // (inclusive) if it is better than the previous solutions.
                if !is_one {
                    debug_assert!(candidate_count < Self::STACK_SIZE);
                    candidate_stack[candidate_count] = (current_node, index);
                    candidate_count += 1;
                    break;
                }
                if zero_rank < self.data[data_index].count_zeros() {
                    // There is at least one zero afterwards.
                    let leftmost_zero_index_after = self.data[data_index]
                        .select_zero(zero_rank)
                        .unwrap();
                    debug_assert!(candidate_count < Self::STACK_SIZE);
                    candidate_stack[candidate_count] = (current_node, leftmost_zero_index_after);
                    candidate_count += 1;
                }
            }

            if min_value_in_right_child >= target_length {
                if min_value_in_right_child - 1 < target_length {
                    // The candidates in the left child have already been checked.
                    break;
                }
                // Move to the left child as there are no possible smaller values in the right
                // child. Update the index so that it is valid in the left child.
                if is_one && zero_rank == self.data[data_index].count_zeros() {
                    // All 0s are before the current index.
                    break;
                }
                index = zero_rank;
                // If the maximum value in the left child is greater than or equal to the target
                // length and the minimum value in the left child which is equal to the minimum
                // value in this node is smaller than the target length, then there are more than 2
                // values in the range of the left child i.e. it is guaranteed to exist.
                current_node = self.tree[current_node].left_child;
                continue;
            }

            // min_value_in_right_child < target_length

            // The candidate indices for the previous smaller value which are in the left child
            // have been handled. Therefore, we should try moving to the right child as there might
            // still be candidates there.

            if !is_one && one_rank == self.data[data_index].count_ones() {
                break;
            }
            index = one_rank;
            // If the minimum value in the right child is smaller than the target length and the
            // maximum value in the right child which is equal to the maximum value in the range of
            // this node is greater than or equal to the target length, then there must be at least
            // 2 different values in the range of the right child i.e. it is guaranteed to exist.
            current_node = self.tree[current_node].right_child;
        }

        if candidate_count == 0 {
            return self.len();
        }

        candidate_count -= 1;
        let (mut result_node, mut result) = candidate_stack[candidate_count];
        while candidate_count > 0 {
            candidate_count -= 1;
            let (current_node, index_in_node) = candidate_stack[candidate_count];
            result = self.climb(result_node, result, current_node);
            result_node = current_node;
            result = result.min(index_in_node);
        }

        self.climb(result_node, result, 0)
    }

    /// Climb up the tree to find the index of a given character in the source string.
    #[inline]
    fn climb(&self, mut node: usize, mut index_in_node: usize, ancestor: usize) -> usize {
        let mut parent;
        let mut parent_data_index;
        while node > ancestor {
            parent = self.tree[node].parent;
            parent_data_index = self.tree[parent].data;
            index_in_node = if self.tree[parent].left_child == node {
                // If the node is the left child of its parent, then the bit representing the
                // current index in the parent has a value of 0.
                match self.data[parent_data_index].select_zero(index_in_node) {
                    Some(value) => value,
                    None => {
                        return self.data[0].len();
                    }
                }
            } else {
                // Otherwise it has a value of 1.
                match self.data[parent_data_index].select(index_in_node) {
                    Some(value) => value,
                    None => {
                        return self.data[0].len();
                    }
                }
            };
            node = parent;
        }
        assert!(node == ancestor);
        index_in_node
    }

    pub fn serialize<W: std::io::Write>(&self, out: &mut W) -> std::io::Result<usize> {
        let mut written: usize = 0;

        out.write_all(&(self.lower_bound as u64).to_le_bytes())?;
        out.write_all(&(self.window_size as u64).to_le_bytes())?;
        let tree_size = self.tree.len();
        assert_eq!(tree_size, self.data.len());
        out.write_all(&(tree_size as u64).to_le_bytes())?;
        written += 3 * size_of::<u64>();

        for (index, data) in self.data.iter().enumerate() {
            log::info!("[WWT::serialize] serializing node data {}...", index);
            data.serialize(out)?;
            written += data.size_in_bytes();
        }

        Ok(written)
    }

    pub fn load<R: std::io::Read>(input: &mut R) -> std::io::Result<Self> {
        use byteorder::{ReadBytesExt, LittleEndian};
        let lower_bound = input.read_u64::<LittleEndian>()? as usize;
        let window_size = input.read_u64::<LittleEndian>()? as usize;
        let tree_size   = input.read_u64::<LittleEndian>()? as usize;
        let mut tree: Vec<Node> = vec![];
        Self::make_tree_nodes(&mut tree, lower_bound, window_size);
        assert_eq!(tree_size, tree.len());
        let mut data = Vec::<BitVector>::with_capacity(tree_size);
        for index in 0..tree_size {
            log::info!("[WWT::load] loading node data {}...", index);
            let row = BitVector::load(input)?;
            data.push(row);
        }
        let result = Self {
            lower_bound,
            window_size,
            tree,
            data,
        };
        Ok(result)
    }

    #[inline]
    pub fn len(&self) -> usize {
        if self.is_empty() {
            0
        } else {
            self.data[0].len()
        }
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.window_size == 0
    }
}

impl super::Pnsv for WindowedWaveletTree {
    #[inline]
    fn previous(&self, index: usize, target_length: usize) -> usize {
        self.previous(index, target_length)
    }

    #[inline]
    fn next(&self, index: usize, target_length: usize) -> usize {
        self.next(index, target_length)
    }

    #[inline]
    fn max_target(&self) -> usize {
        (self.lower_bound + self.window_size).saturating_sub(1)
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
        let lower_bound = 3;
        let window_size = 6;
        let wavelet = WindowedWaveletTree::from_iterator(items.iter().cloned(), items.len(), lower_bound, window_size);

        let mut buffer = Vec::<u8>::new();
        let written = wavelet.serialize(&mut buffer).unwrap();
        assert_eq!(buffer.len(), written);
        let wavelet_loaded = WindowedWaveletTree::load(&mut buffer.as_slice()).unwrap();
        assert_eq!(wavelet_loaded, wavelet);
    }

    fn previous(items: &[usize], index: usize, target_length: usize, lower_bound: usize, upper_bound: usize) -> usize {
        for i in (0..=index).rev() {
            let item = items[i].clamp(lower_bound, upper_bound);
            if item < target_length {
                return i;
            }
        }
        0
    }

    fn next(items: &[usize], index: usize, target_length: usize, lower_bound: usize, upper_bound: usize) -> usize {
        for i in index..items.len() {
            let item = items[i].clamp(lower_bound, upper_bound);
            if item < target_length {
                return i;
            }
        }
        items.len()
    }

    fn test_with_parameters(items: &[usize], lower_bound: usize, window_size: usize) {
        let upper_bound = lower_bound + window_size - 1;
        let wavelet = WindowedWaveletTree::from_iterator(items.iter().cloned(), items.len(), lower_bound, window_size);
        for i in 0..items.len() {
            for target_length in lower_bound + 1..=lower_bound + window_size {
                assert_eq!(
                    previous(items, i, target_length, lower_bound, upper_bound),
                    wavelet.previous(i, target_length),
                    "ws: {window_size}, previous; i: {}, target_length: {}",
                    i,
                    target_length
                );
                assert_eq!(
                    next(items, i, target_length, lower_bound, upper_bound),
                    wavelet.next(i, target_length),
                    "ws: {window_size}, next; i: {}, target_length: {}",
                    i,
                    target_length
                );
            }
        }
    }

    #[test]
    fn windowed_wavelet_tree_all_01() {
        let items: &[usize] = &[
            2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9, 9,
            9, 9, 9, 9, 8, 8, 8, 8, 7, 7, 7, 7, 6, 6, 6, 6, 5, 5, 5, 5, 4, 4, 4, 4, 3, 3, 3, 3, 2, 2, 2, 2,
        ];
        test_with_parameters(items, 4, 3);
        test_with_parameters(items, 5, 4);
        test_with_parameters(items, 2, 6);
    }

    #[test]
    fn windowed_wavelet_tree_all_02() {
        let items: &[usize] = &[
            2, 2, 2, 3, 3, 4, 3, 3, 1, 4, 6, 4, 5, 8, 8, 5, 6, 6, 7, 6, 7, 3, 2, 7, 8, 6, 6, 6, 9, 9, 9, 9,
            9, 2, 9, 9, 3, 8, 8, 8, 4, 4, 7, 7, 6, 5, 6, 6, 8, 8, 2, 5, 4, 4, 8, 8, 9, 3, 3, 3, 3, 2, 2, 2,
        ];
        test_with_parameters(items, 4, 3);
        test_with_parameters(items, 5, 4);
        test_with_parameters(items, 2, 6);
    }
}

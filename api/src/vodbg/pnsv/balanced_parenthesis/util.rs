pub fn one_indices_iterator<B: Iterator<Item = bool> + Clone>(iterator: B) -> impl Iterator<Item = usize> + Clone {
    iterator.enumerate().filter(|(_, value)| *value).map(|(index, _)| index)
}

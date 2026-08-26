use simple_sds_sbwt::serialize::Serialize;

pub trait StreamBuilder<T: Send + Sync> {
    type Stream<'a>: Iterator<Item = T> + Send + 'a where Self: 'a;
    fn build<'a>(&'a self, offset: usize) -> Self::Stream<'a>;
}

pub struct MemoryStream<T> {
    data: Vec<T>,
}

impl<T> MemoryStream<T> {
    pub fn new(data: Vec<T>) -> Self {
        Self { 
            data
        }
    }
}

impl<T> From<Vec<T>> for MemoryStream<T> {
    fn from(value: Vec<T>) -> Self {
        Self { data: value }
    }
}

impl<T> From<MemoryStream<T>> for Vec<T> {
    fn from(value: MemoryStream<T>) -> Self {
        value.data
    }
}

use std::iter::Copied;
use std::slice::Iter;
impl<T> StreamBuilder<T> for MemoryStream<T>
where T: Send + Sync + Copy
{
    type Stream<'a> = Copied<Iter<'a, T>> where Self: 'a;
    fn build<'a>(&'a self, offset: usize) -> Self::Stream<'a> {
        self.data[offset..].iter().copied()
    }
}

pub struct DiskStream<T: CWS> {
    file: crate::tempfile::TempFile,
    _phantom: std::marker::PhantomData<T>,
}

impl<T: CWS> DiskStream<T> {
    pub fn new(data: Vec<T>, temp_file_manager: &mut crate::tempfile::TempFileManager) -> Self {
        use simple_sds_sbwt::serialize::Serialize;
        let mut file = temp_file_manager.create_new_file("suffix_array", 16, ".bin");

        // The integers from the file will be read by the same process on a machine with the same
        // architecture, so there shouldn't be worry about endianness.
        data.serialize_body(&mut file).unwrap();

        Self {
            file,
            _phantom: Default::default(),
        }
    }
}

pub struct DiskStreamIterator<'a, T: CWS> {
    reader: std::io::BufReader<std::fs::File>,
    /// Keep a reference to the stream, since if it is dropped, the temporary file will be deleted.
    _disk_stream: &'a DiskStream<T>,
}

impl<T: CWS> StreamBuilder<T> for DiskStream<T>
{
    type Stream<'a> = DiskStreamIterator<'a, T> where Self: 'a;

    fn build<'a>(&'a self, offset: usize) -> Self::Stream<'a> {
        use std::io::{Seek, SeekFrom};
        let byte_offset = (offset * T::byte_size()) as u64;
        let mut file = std::fs::File::open(&self.file.path).unwrap();
        file.seek(SeekFrom::Start(byte_offset)).unwrap();
        let reader = std::io::BufReader::new(file);
        DiskStreamIterator {
            reader,
            _disk_stream: self,
        }
    }
}

impl<'a, T: CWS> Iterator for DiskStreamIterator<'a, T> {
    type Item = T;
    fn next(&mut self) -> Option<Self::Item> {
        let result = T::load(&mut self.reader);
        result.ok()
    }
}

use constant_width_serializable as cws;
use cws::ConstantWidthSerializable as CWS;
use cws::constant_width_serializable;
constant_width_serializable!(
    CWS;
    usize,
    u64,
);

pub mod constant_width_serializable {
    pub trait ConstantWidthSerializable:
        simple_sds_sbwt::serialize::Serializable
        + Send + Sync
    {
        fn byte_size() -> usize;
    }

    #[macro_export]
    macro_rules! constant_width_serializable {
        ($cws:ident; $($t:ty),+ $(,)?) => {
            $(
                impl $cws for $t {
                    #[inline(always)]
                    fn byte_size() -> usize {
                        size_of::<$t>()
                    }
                }
            )+
        };
    }

    pub use constant_width_serializable;
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_number_vector(length: usize) -> Vec<usize> {
        let mut vector = Vec::with_capacity(length);
        vector.extend(0..length);
        vector
    }

    fn iterators_return_the_same_values(
        mut iterator_a: impl Iterator<Item = usize>,
        mut iterator_b: impl Iterator<Item = usize>,
    ) {
        loop {
            let item_a = iterator_a.next();
            let item_b = iterator_b.next();
            
            if item_a.is_none() && item_b.is_none() {
                break;
            }
            assert_eq!(item_a, item_b);
        }
    }

    #[test]
    fn memory_stream() {
        let length = 64;
        let invariant = make_number_vector(length);
        let numbers = make_number_vector(length);
        let in_memory_stream = MemoryStream::new(numbers);
        for offset in 0..length {
            iterators_return_the_same_values(
                in_memory_stream.build(offset),
                invariant.iter().copied().skip(offset)
            );
        }
    }

    #[test]
    fn disk_stream() {
        let length = 64;
        let invariant = make_number_vector(length);
        let numbers = make_number_vector(length);
        let mut temp_file_manager = crate::tempfile::TempFileManager::new(std::path::Path::new("."));
        let in_memory_stream = DiskStream::new(numbers, &mut temp_file_manager);
        for offset in 0..length {
            iterators_return_the_same_values(
                in_memory_stream.build(offset),
                invariant.iter().copied().skip(offset)
            );
        }
    }
}

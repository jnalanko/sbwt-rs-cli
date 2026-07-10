//! Disk random-access primitives for reading fixed-size `LongKmer` records (and dummy
//! kmer/length pairs) by global sorted index from the k-mer file and dummy file produced by
//! the disk-based construction pipeline. Both files are contiguous, globally sorted arrays of
//! fixed-size records, so `seek(idx * record_len)` followed by a single record read is valid
//! random access.

use std::fs::File;
use std::io::{BufReader, Read, Seek, SeekFrom};
use std::path::Path;

use crate::kmer::LongKmer;

pub(crate) fn kmer_record_len<const B: usize>() -> usize {
    LongKmer::<B>::byte_size()
}

pub(crate) fn dummy_record_len<const B: usize>() -> usize {
    LongKmer::<B>::byte_size() + 1
}

pub(crate) fn n_kmer_records<const B: usize>(path: &Path) -> usize {
    let len = std::fs::metadata(path).unwrap().len() as usize;
    let record_len = kmer_record_len::<B>();
    assert_eq!(len % record_len, 0);
    len / record_len
}

#[allow(dead_code)] // Kept for symmetry with n_kmer_records and for tests.
pub(crate) fn n_dummy_records<const B: usize>(path: &Path) -> usize {
    let len = std::fs::metadata(path).unwrap().len() as usize;
    let record_len = dummy_record_len::<B>();
    assert_eq!(len % record_len, 0);
    len / record_len
}

pub(crate) fn read_kmer_at<const B: usize>(path: &Path, idx: usize) -> LongKmer<B> {
    let record_len = kmer_record_len::<B>();
    let mut f = File::open(path).unwrap();
    f.seek(SeekFrom::Start(idx as u64 * record_len as u64)).unwrap();
    LongKmer::<B>::load(&mut f).unwrap().unwrap() // Should never be None because the caller guarantees idx is in bounds.
}

pub(crate) fn read_dummy_at<const B: usize>(path: &Path, idx: usize) -> (LongKmer<B>, u8) {
    let record_len = dummy_record_len::<B>();
    let mut f = File::open(path).unwrap();
    f.seek(SeekFrom::Start(idx as u64 * record_len as u64)).unwrap();
    let kmer = LongKmer::<B>::load(&mut f).unwrap().unwrap();
    let mut len_buf = [0_u8; 1];
    f.read_exact(&mut len_buf).unwrap();
    (kmer, len_buf[0])
}

pub(crate) fn seek_kmer_reader<const B: usize>(path: &Path, idx: usize) -> BufReader<File> {
    let record_len = kmer_record_len::<B>();
    let mut f = File::open(path).unwrap();
    f.seek(SeekFrom::Start(idx as u64 * record_len as u64)).unwrap();
    BufReader::with_capacity(1 << 23, f) // 8 MiB bufferr
}

pub(crate) fn seek_dummy_reader<const B: usize>(path: &Path, idx: usize) -> BufReader<File> {
    let record_len = dummy_record_len::<B>();
    let mut f = File::open(path).unwrap();
    f.seek(SeekFrom::Start(idx as u64 * record_len as u64)).unwrap();
    BufReader::with_capacity(1 << 23, f) // 8 MiB bufferr
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn read_kmer_and_dummy_at() {
        let mut temp_file_manager = crate::tempfile::TempFileManager::new(std::path::Path::new("/tmp"));
        let mut kmer_file = temp_file_manager.create_new_file("test-disk-access-", 10, ".kmers");
        let mut dummy_file = temp_file_manager.create_new_file("test-disk-access-", 10, ".dummies");

        let kmers = [
            LongKmer::<2>::from_ascii(b"AAAA").unwrap(),
            LongKmer::<2>::from_ascii(b"ACGT").unwrap(),
            LongKmer::<2>::from_ascii(b"TTTT").unwrap(),
        ];
        for kmer in kmers.iter() {
            kmer.serialize(&mut kmer_file).unwrap();
        }
        kmer_file.flush().unwrap();

        let dummies = [
            (LongKmer::<2>::from_ascii(b"AAAA").unwrap(), 0_u8),
            (LongKmer::<2>::from_ascii(b"ACAA").unwrap(), 2_u8),
        ];
        for (kmer, len) in dummies.iter() {
            kmer.serialize(&mut dummy_file).unwrap();
            dummy_file.write_all(&[*len]).unwrap();
        }
        dummy_file.flush().unwrap();

        assert_eq!(n_kmer_records::<2>(&kmer_file.path), 3);
        assert_eq!(n_dummy_records::<2>(&dummy_file.path), 2);

        for (i, expected) in kmers.iter().enumerate() {
            assert_eq!(read_kmer_at::<2>(&kmer_file.path, i), *expected);
        }
        for (i, expected) in dummies.iter().enumerate() {
            assert_eq!(read_dummy_at::<2>(&dummy_file.path, i), *expected);
        }

        // seek_kmer_reader/seek_dummy_reader should yield the same sequence as sequential reads from that point.
        let mut reader = seek_kmer_reader::<2>(&kmer_file.path, 1);
        assert_eq!(LongKmer::<2>::load(&mut reader).unwrap().unwrap(), kmers[1]);
        assert_eq!(LongKmer::<2>::load(&mut reader).unwrap().unwrap(), kmers[2]);
        assert_eq!(LongKmer::<2>::load(&mut reader).unwrap(), None);
    }
}

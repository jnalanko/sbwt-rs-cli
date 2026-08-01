// Code by Martin Kostadinov.

#![allow(unused)]

use crate::streaming_index::{ContractLeft, ExtendRight, StreamingIndex};
use crate::SbwtIndexVariant;
use crate::streaming_index::LcsArray;

use super::pnsv::Pnsv;

// note(mk): Temporary solution for reading input data.
pub fn read_index_and_lcs(arguments_start: usize) -> (SbwtIndexVariant, LcsArray) {
    let mut args = std::env::args().skip(arguments_start);
    let sbwt_path = args.next().expect("expected sbwt index path");
    let lcs_path = args.next().expect("expected lcs index path");

    let mut index_reader = std::io::BufReader::new(std::fs::File::open(sbwt_path).unwrap());
    let index = crate::load_sbwt_index_variant(&mut index_reader).unwrap();

    let mut lcs_reader = std::io::BufReader::new(std::fs::File::open(lcs_path).unwrap());
    let lcs = LcsArray::load(&mut lcs_reader).unwrap();

    (index, lcs)
}

pub fn read_sequences<R: std::io::BufRead + 'static + Send + Sync>(reader: R) -> Vec<Vec<u8>> {
    use crate::SeqStream;
    let mut stream = crate::JSeqIOSeqStreamWrapper {
        inner: jseqio::reader::DynamicFastXReader::new(reader).unwrap(),
    };
    let mut result: Vec<Vec<u8>> = vec![];
    while let Some(sequence) = stream.stream_next() {
        result.push(sequence.into());
    }
    result
}

pub fn benchmark_bms_separate_queries<E, C>(
    index: &StreamingIndex<'_, E, C>,
    queries: &[Vec<u8>],
    bound: usize,
) where
    E: ExtendRight,
    C: ContractLeft,
{
    let start_time = std::time::Instant::now();
    let mut checksum = 0_usize;
    let mut n_kmers_queried = 0_usize;
    for query in queries {
        for x in index.bounded_matching_statistics(query, bound).iter() {
            checksum += x.0;
        }
        n_kmers_queried += std::cmp::max(query.len() as isize - index.k as isize + 1, 0) as usize;
    }
    let end_time = std::time::Instant::now();
    print!("{:.2},", (end_time - start_time).as_secs_f64());
    let nanos_per_kmer = (end_time - start_time).as_nanos() as f64 / n_kmers_queried as f64;
    print!("{:.2},{}", nanos_per_kmer, checksum);
}

#[allow(unused)]
#[cfg(test)]
mod tests {
    // note(mk): moved to separate project. Will write something more useful here later.
}

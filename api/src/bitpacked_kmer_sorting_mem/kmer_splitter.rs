use crate::kmer::LongKmer;

use rayon::iter::IntoParallelRefIterator;
use rayon::iter::IntoParallelRefMutIterator;
use rayon::iter::ParallelIterator;
use rayon::prelude::ParallelSlice;
use rayon::prelude::ParallelSliceMut;

pub fn split_to_bins<const B: usize, IN: crate::SeqStream + Send>(
    mut seqs: IN,
    k: usize,
    dedup_batches: bool,
) -> Vec<std::io::Cursor::<Vec<LongKmer::<B>>>> {

    // Beware: changing the number of bins, or their order, messes up the
    // results and may not always break all tests.
    //
    // This function makes weird implicit assumptions that are violated if the
    // operations are performed differently?

    let mut buf = Vec::<Box<[u8]>>::new();
    while let Some(seq) = seqs.stream_next() {
        let mut seq_copy = seq.to_owned();
        seq_copy.reverse(); // Reverse to get colex sorting
        buf.push(seq_copy.into_boxed_slice());
    }

    let mut kmers = buf.par_iter().map(|seq| {
        seq.windows(k).filter_map(|bytes| {
            LongKmer::<B>::from_ascii(bytes).ok()
        }).collect::<Vec<LongKmer::<B>>>()
    }).flatten().collect::<Vec<LongKmer::<B>>>();

    kmers.par_sort_unstable_by_key(|kmer| {
        kmer.get_from_left(0) as usize * 16 + kmer.get_from_left(1) as usize * 4 + kmer.get_from_left(2) as usize
    });

    kmers.par_chunk_by(|&a, &b| {
        let x = a.get_from_left(0) as usize * 16 + a.get_from_left(1) as usize * 4 + a.get_from_left(2) as usize;
        let y = b.get_from_left(0) as usize * 16 + b.get_from_left(1) as usize * 4 + b.get_from_left(2) as usize;
        x == y}).map(|chunk| {
        if dedup_batches {
            let mut dd_chunk = chunk.to_vec();
            dd_chunk.sort_unstable();
            dd_chunk.dedup();
            std::io::Cursor::<Vec<LongKmer::<B>>>::new(dd_chunk)
        } else {
            std::io::Cursor::<Vec<LongKmer::<B>>>::new(chunk.to_vec())
        }
    }).collect::<Vec<std::io::Cursor::<Vec<LongKmer::<B>>>>>()
}

// Overwrite the bins with sorted and deduplicates files. Returns back the files after overwriting.
pub fn par_sort_and_dedup_bin_files<const B: usize>(
    bins: &mut Vec<std::io::Cursor::<Vec<LongKmer::<B>>>>,
) {
   bins.par_iter_mut().for_each(|f| {
        f.get_mut().sort_unstable();
        f.get_mut().dedup();
        f.set_position(0);
    });
}

// The original data are deleted
pub fn concat_files<const B: usize>(
    binned_kmers: &mut Vec<std::io::Cursor::<Vec<LongKmer::<B>>>>
) -> std::io::Cursor::<Vec<LongKmer::<B>>> {
    let mut concat_kmers: Vec<LongKmer::<B>> = Vec::new();
    binned_kmers.iter_mut().for_each(|file| concat_kmers.append(file.get_mut()));
    std::io::Cursor::<Vec<LongKmer::<B>>>::new(Vec::from(concat_kmers))
}

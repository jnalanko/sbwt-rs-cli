use std::io::{BufWriter, Write};

use super::kmer::LongKmer;
use crate::tempfile::TempFileManager;
use simple_sds_sbwt::raw_vector::*;
use rayon::prelude::*;

struct NullReader{}

impl std::io::Read for NullReader{
    fn read(&mut self, _buf: &mut [u8]) -> std::io::Result<usize>{
        Ok(0) // EoF
    }
}

// We take in a path and not a file object because we need multiple readers to the same file
pub fn get_sorted_dummies<const B: usize>(sorted_kmers_filepath: &std::path::Path, sigma: usize, k: usize, temp_file_manager: &mut TempFileManager) -> Vec<(LongKmer::<B>, u8)>{

    // Todo: I'm using dummy merger cursors with an empty dummy file. Should refactor things to manage without the
    // empty dummy file.

    // Number of k-mers in file
    let n = std::fs::metadata(sorted_kmers_filepath).unwrap().len() as usize / LongKmer::<B>::byte_size();

    let mut has_predecessor = simple_sds_sbwt::raw_vector::RawVector::new();
    has_predecessor.resize(n, false);

    let emptyfile = temp_file_manager.create_new_file("empty-", 10, ".bin");
    let mut char_cursors = crate::bitpacked_kmer_sorting::cursors::init_char_cursors::<B>(&emptyfile.path, sorted_kmers_filepath, k, sigma);

    let global_cursor = crate::bitpacked_kmer_sorting::cursors::DummyNodeMerger::new(
        std::io::BufReader::new(std::fs::File::open(&emptyfile.path).unwrap()),
        std::io::BufReader::new(std::fs::File::open(sorted_kmers_filepath).unwrap()),
        k,
    );

    for (x, _) in global_cursor{
        // x is reversed
        for c in 0..(sigma as u8){
            // Shiting a reversed k-mer to the right means shifting the original k-mer to the left
            let xc = x.set_from_left(k-1, 0).right_shift(1).set_from_left(0, c);
            while char_cursors[c as usize].peek().is_some(){
                match char_cursors[c as usize].peek().unwrap().0.cmp(&xc) {
                    std::cmp::Ordering::Greater => break,
                    std::cmp::Ordering::Equal => {
                        has_predecessor.set_bit(char_cursors[c as usize].nondummy_position(), true);
                        char_cursors[c as usize].next(); // Advance
                        break
                    },
                    std::cmp::Ordering::Less => {
                        char_cursors[c as usize].next(); // Advance
                        // no break
                    }
                }
            }
        }
    }

    let mut global_cursor = crate::bitpacked_kmer_sorting::cursors::DummyNodeMerger::new(
        std::io::BufReader::new(std::fs::File::open(&emptyfile.path).unwrap()),
        std::io::BufReader::new(std::fs::File::open(sorted_kmers_filepath).unwrap()),
        k,
    ); // New global cursor

    // Todo: stream to memory and sort there
    let mut required_dummies = Vec::<(LongKmer::<B>, u8)>::new(); // Pairs (data, length)

    while let Some((x, _)) = global_cursor.peek(){
        if !has_predecessor.bit(global_cursor.nondummy_position()){
            let mut prefix = x;
            for i in 0..k {
                let len = k - i - 1;
                prefix = prefix.left_shift(1);
                required_dummies.push((prefix, len as u8));
            }
        }
        global_cursor.next(); // Advance
    }

    // We always assume that the empty k-mer exists. This assumption is reflected in the C-arrya
    // later, which adds one "ghost dollar" count to all counts.
    required_dummies.push((LongKmer::<B>::from_ascii(b"").unwrap(), 0));

    required_dummies.par_sort();
    required_dummies.dedup();
    required_dummies.shrink_to_fit();
    required_dummies
    

}

pub fn write_to_disk<const B: usize>(dummies: Vec<(LongKmer::<B>, u8)>, writer: &mut std::fs::File){   
    let mut bw = BufWriter::new(writer);
    for (kmer, len) in dummies.iter(){
        kmer.serialize(&mut bw).unwrap();
        bw.write_all(&[*len]).unwrap();
    }
    bw.flush().unwrap();
}
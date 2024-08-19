use std::{io::{BufReader, Seek, Read}, fs::File, path::Path};

use simple_sds_sbwt::{ops::Access, raw_vector::AccessRaw};
use std::io::SeekFrom;
use std::cmp::min;
use super::kmer::LongKmer;
use crate::util::binary_search_leftmost_that_fulfills_pred;

pub struct DummyNodeMerger<R: std::io::Read, const B: usize> {
    dummy_reader: R, // Stream of k-mer objects
    nondummy_reader: R, // Stream of pairs (kmer, len)

    dummy_kmer: Option<(LongKmer::<B>, u8)>,
    nondummy_kmer: Option<(LongKmer::<B>, u8)>,

    k: usize,

    dummy_position: usize, // Position of the dummy cursor
    nondummy_position: usize, // Position of the nondummy cursor
}

impl <R: std::io::Read, const B: usize> DummyNodeMerger<R, B> {

    pub fn read_from_dummy_reader(dummy_reader: &mut R) -> Option<(LongKmer::<B>, u8)>{
        let kmer = match LongKmer::<B>::load(dummy_reader){
            Ok(kmer_opt) => {
                match kmer_opt{
                    Some(kmer) => kmer,
                    None => return None, // End of stream
                }   
            },
            Err(e) => panic!("IO error while streaming sorted k-mers: {}", e),
        };

        // Read length
        let mut buf = [0_u8; 1];
        dummy_reader.read_exact(&mut buf).unwrap();
        let len = u8::from_le_bytes(buf);

        Some((kmer, len))
    }

    pub fn read_from_non_dummy_reader(nondummy_reader: &mut R, k: usize) -> Option<(LongKmer::<B>, u8)>{
        match LongKmer::<B>::load(nondummy_reader){
            Ok(kmer_opt) => {
                kmer_opt.map(|kmer| (kmer, k as u8))
            },
            Err(e) => panic!("IO error while streaming sorted k-mers: {}", e),
        }
    }

    pub fn new(mut dummy_reader: R, mut nondummy_reader: R, k: usize) -> Self {
        let dummy_kmer = Self::read_from_dummy_reader(&mut dummy_reader);
        let nondummy_kmer = Self::read_from_non_dummy_reader(&mut nondummy_reader, k);

        Self {
            dummy_reader,
            nondummy_reader,
            dummy_kmer,
            nondummy_kmer,
            k,
            dummy_position: 0,
            nondummy_position: 0,
        }
    }

    // TODO: this is stupid. The given positions are just for bookkeeping that the caller might use. They don't affect the
    // cursor at all. Need to refactor.
    pub fn new_with_initial_positions(mut dummy_reader: R, mut nondummy_reader: R, k: usize, dummy_position: usize, nondummy_position: usize) -> Self {
        let dummy_kmer = Self::read_from_dummy_reader(&mut dummy_reader);
        let nondummy_kmer = Self::read_from_non_dummy_reader(&mut nondummy_reader, k);

        Self {
            dummy_reader,
            nondummy_reader,
            dummy_kmer,
            nondummy_kmer,
            k,
            dummy_position,
            nondummy_position,
        }
    }

    pub fn peek(&self) -> Option<(LongKmer::<B>, u8)>{
        match (self.dummy_kmer, self.nondummy_kmer){
            (None, None) => None,
            (Some(dummy_kmer), None) => {
                Some(dummy_kmer)
            },
            (None, Some(nondummy_kmer)) => {
                Some(nondummy_kmer)
            },
            (Some(dummy_kmer), Some(nondummy_kmer)) => {
                if dummy_kmer < nondummy_kmer {
                    Some(dummy_kmer)
                } else {
                    Some(nondummy_kmer)
                }
            }
        }
    }

    #[allow(dead_code)]
    pub fn dummy_position(&self) -> usize{
        self.dummy_position
    }

    pub fn nondummy_position(&self)  -> usize{
        self.nondummy_position
    }
}

impl<const B: usize> Iterator for DummyNodeMerger<BufReader<File>, B> {
    type Item = (LongKmer<B>, u8);

    // Produces pairs (kmer, length)
    fn next(&mut self) -> Option<(LongKmer::<B>, u8)> {
        match (self.dummy_kmer, self.nondummy_kmer){
            (None, None) => None,
            (Some(dummy_kmer), None) => {
                self.dummy_kmer = Self::read_from_dummy_reader(&mut self.dummy_reader);
                self.dummy_position += 1;
                Some(dummy_kmer)
            },
            (None, Some(nondummy_kmer)) => {
                self.nondummy_kmer = Self::read_from_non_dummy_reader(&mut self.nondummy_reader, self.k);
                self.nondummy_position += 1;
                Some(nondummy_kmer)
            },
            (Some(dummy_kmer), Some(nondummy_kmer)) => {
                if dummy_kmer < nondummy_kmer {
                    self.dummy_kmer = Self::read_from_dummy_reader(&mut self.dummy_reader);
                    self.dummy_position += 1;
                    Some(dummy_kmer)
                } else {
                    self.nondummy_kmer = Self::read_from_non_dummy_reader(&mut self.nondummy_reader, self.k);
                    self.nondummy_position += 1;
                    Some(nondummy_kmer)
                }
            }
        }
    }

}

// We take in Paths instead of a Files because we need multiple readers to the same files 
pub fn init_char_cursors<const B: usize>(dummy_filepath: &Path, nondummy_filepath: &Path, k: usize, sigma: usize)
-> Vec<DummyNodeMerger<BufReader<File>, B>>{
    let mut char_cursors = Vec::<DummyNodeMerger<BufReader<File>, B>>::new();
    for c in 0..(sigma as u8){
        log::trace!("Searching character {}", c);

        let (dummy_reader, dummy_pos) = 
        { // Seek in dummies

            let dummy_file_len = std::fs::metadata(dummy_filepath).unwrap().len() as usize;
            let dummy_record_len = LongKmer::<B>::byte_size() + 1; // Pairs (kmer, len byte)
            assert_eq!(dummy_file_len % dummy_record_len, 0);
    
            let access_fn = |pos| {
                let mut f = File::open(dummy_filepath).unwrap();
                f.seek(SeekFrom::Start(pos as u64 * dummy_record_len as u64)).unwrap();
                let kmer = LongKmer::<B>::load(&mut f).unwrap().unwrap(); // Should never be none because we know the file length

                // Read the length byte
                let mut len_buf = [0_u8; 1];
                f.read_exact(&mut len_buf).unwrap();
                let len = u8::from_le_bytes(len_buf);
                (kmer, len)
            };

            let pred_fn = |kmer: (LongKmer::<B>,u8)| {
                kmer.1 > 0 && kmer.0.get_from_left(0) >= c
            };

            let start = binary_search_leftmost_that_fulfills_pred(
                access_fn, 
                pred_fn, 
                dummy_file_len / dummy_record_len);

            let mut f = File::open(dummy_filepath).unwrap();
            f.seek(SeekFrom::Start(start as u64 * dummy_record_len as u64)).unwrap();
            (BufReader::new(f), start)
        };

        let (nondummy_reader, nondummy_pos) = 
        { // Seek in nondummies

            let nondummy_file_len = std::fs::metadata(nondummy_filepath).unwrap().len() as usize;
            let nondummy_record_len = LongKmer::<B>::byte_size();
            assert_eq!(nondummy_file_len % nondummy_record_len, 0);
    
            let access_fn = |pos| {
                let mut f = File::open(nondummy_filepath).unwrap();
                f.seek(SeekFrom::Start(pos as u64 * nondummy_record_len as u64)).unwrap();
                LongKmer::<B>::load(&mut f).unwrap().unwrap() // Should never be None because we know the file length
            };

            let pred_fn = |kmer: LongKmer::<B>| {
                kmer.get_from_left(0) >= c
            };

            let start = binary_search_leftmost_that_fulfills_pred(
                access_fn, 
                pred_fn, 
                nondummy_file_len / nondummy_record_len);
        
            let mut f = File::open(nondummy_filepath).unwrap();
            f.seek(SeekFrom::Start(start as u64 * nondummy_record_len as u64)).unwrap();
            (BufReader::new(f), start)
        };

        let cursor = DummyNodeMerger::new_with_initial_positions(dummy_reader, nondummy_reader, k, dummy_pos, nondummy_pos);
        char_cursors.push(cursor);
    }

    char_cursors

}

// Returns the SBWT bit vectors and optionally the LCS array
pub fn build_sbwt_bit_vectors<const B: usize>(
    global_cursor: DummyNodeMerger<BufReader<File>, B>, 
    mut char_cursors: Vec<DummyNodeMerger<BufReader<File>, B>>, 
    n: usize,
    k: usize, 
    sigma: usize,
    build_lcs: bool) -> (Vec<simple_sds_sbwt::raw_vector::RawVector>, Option<simple_sds_sbwt::int_vector::IntVector>)
{

    let mut rawrows = Vec::<simple_sds_sbwt::raw_vector::RawVector>::new();
    for _ in 0..sigma {
        rawrows.push(simple_sds_sbwt::raw_vector::RawVector::with_len(n, false));
    }

    let mut lcs = if build_lcs { 
        // LCS values are between 0 and k-1
        assert!(k > 0);
        let bitwidth = 64 - (k as u64 - 1).leading_zeros();
        Some(simple_sds_sbwt::int_vector::IntVector::with_len(n, bitwidth as usize, 0).unwrap()) } 
    else { 
        None 
    };

    let mut prev_kmer = LongKmer::<B>::from_ascii(b"").unwrap();
    let mut prev_len = 0_usize;
    for (kmer_idx, (kmer, len)) in global_cursor.enumerate() {
        // The k-mers enumerated are reversed

        if build_lcs && kmer_idx > 0 {
            // The longest common suffix is the longest common prefix of reversed k-mers
            let mut lcs_value = LongKmer::<B>::lcp(&prev_kmer, &kmer);
            lcs_value = min(lcs_value, min(prev_len, len as usize));
            lcs.as_mut().unwrap().set(kmer_idx, lcs_value as u64);
        }

        for c in 0..(sigma as u8) {
            let kmer_c = if len as usize == k {
                (
                    kmer.clone()
                        .set_from_left(k - 1, 0)
                        .right_shift(1)
                        .set_from_left(0, c),
                    k as u8,
                )
            } else {
                (kmer.clone().right_shift(1).set_from_left(0, c), len + 1) // Dummy
            };

            while char_cursors[c as usize].peek().is_some() && char_cursors[c as usize].peek().unwrap() < kmer_c {
                char_cursors[c as usize].next();
            }

            if char_cursors[c as usize].peek().is_some() && char_cursors[c as usize].peek().unwrap() == kmer_c {
                rawrows[c as usize].set_bit(kmer_idx, true);
                char_cursors[c as usize].next();
            }
        }
        prev_kmer = kmer;
        prev_len = len as usize;
    }

    (rawrows, lcs)

}

#[cfg(test)]
mod tests{

    use std::io::Write;

    use super::*;

    #[test]
    fn test_init_char_cursors(){
        let nondummies = [
            LongKmer::<2>::from_ascii(b"ACGT").unwrap(),
            LongKmer::<2>::from_ascii(b"AGGT").unwrap(),
            LongKmer::<2>::from_ascii(b"GGAA").unwrap(),
            LongKmer::<2>::from_ascii(b"GGGT").unwrap()
        ];
        let dummies = [
            (LongKmer::<2>::from_ascii(b"AAAA").unwrap(),0), // This is actually the empty dummy so it's not in the A-block
            (LongKmer::<2>::from_ascii(b"AAAA").unwrap(),1),
            (LongKmer::<2>::from_ascii(b"ACAA").unwrap(),2),
            (LongKmer::<2>::from_ascii(b"ACAA").unwrap(),3),
            (LongKmer::<2>::from_ascii(b"GGTT").unwrap(),3),
        ];

        let mut temp_file_manager = crate::tempfile::TempFileManager::new(std::path::Path::new("temp"));

        let mut nondummy_file = temp_file_manager.create_new_file("test-", 10, ".nondummy");
        let mut dummy_file = temp_file_manager.create_new_file("test-", 10, ".dummy");
        let nondummy_path = nondummy_file.path.clone();
        let dummy_path = dummy_file.path.clone();

        for kmer in nondummies.iter(){
            kmer.serialize(&mut nondummy_file).unwrap();
        }
        for (kmer, len) in dummies.iter(){
            kmer.serialize(&mut dummy_file).unwrap();
            let len_byte = *len as u8;
            dummy_file.write_all(&[len_byte]).unwrap();
        }

        // Flush
        dummy_file.flush().unwrap();
        nondummy_file.flush().unwrap();

        let char_cursors = init_char_cursors(&dummy_path, &nondummy_path, 4, 4);

        assert_eq!(char_cursors[0].peek(), Some((LongKmer::<2>::from_ascii(b"AAAA").unwrap(), 1))); // A
        assert_eq!(char_cursors[1].peek(), Some((LongKmer::<2>::from_ascii(b"GGAA").unwrap(), 4))); // C
        assert_eq!(char_cursors[2].peek(), Some((LongKmer::<2>::from_ascii(b"GGAA").unwrap(), 4))); // G
        assert_eq!(char_cursors[3].peek(), None); // T

    }
}
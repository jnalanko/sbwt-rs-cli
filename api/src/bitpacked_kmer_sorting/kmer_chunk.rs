use std::cmp::*;
use super::kmer::LongKmer;

// B is the number of u64 words in a k-mer
pub struct KmerChunk<const B: usize>{
    pub kmers: Vec<LongKmer::<B>>,
}

pub struct SortedKmerChunk<const B: usize>{
    pub kmers: Vec<LongKmer::<B>>,
}

impl<const B: usize> KmerChunk<B>{
    #[allow(dead_code)]
    pub fn new<I: IntoIterator<Item = LongKmer::<B>>>(input: I) -> Self{
        Self{kmers: input.into_iter().collect()}
    }

    pub fn sort(mut self) -> SortedKmerChunk<B>{
        self.kmers.sort_unstable();
        SortedKmerChunk{kmers: self.kmers}
    }

    /*pub fn par_sort_unstable(mut self) -> SortedKmerChunk<B>{
        self.kmers.par_sort_unstable();
        SortedKmerChunk{kmers: self.kmers}
    }*/

    #[allow(dead_code)]
    pub fn into_vec(self) -> Vec<LongKmer::<B>>{
        self.kmers
    }

    pub fn load<R: std::io::Read>(input: &mut R) -> std::io::Result<Self>{

        let mut kmers = Vec::<LongKmer::<B>>::new();
        while let Some(kmer) = LongKmer::<B>::load(input)?{
            kmers.push(kmer);
        }
        kmers.shrink_to_fit(); // TODO allocate up front
        Ok(Self{kmers})
    }
}

impl<const B: usize> SortedKmerChunk<B>{

    #[cfg(test)]
    pub fn into_vec(self) -> Vec<LongKmer::<B>>{
        self.kmers
    }

    pub fn dedup(&mut self){
        self.kmers.dedup();
        self.kmers.shrink_to_fit();
    }

    // Consumes 'out' to ensure flushing
    #[allow(dead_code)]
    pub fn merge_serialized<R: std::io::Read, W: std::io::Write>(mut inputs: Vec<R>, mut out: W) -> std::io::Result<()>{
        let m = inputs.len();

        let mut n_read = vec![0_usize; m]; // Positions of the "cursors" in the streams

        // Initialize a priority queue pairs (k-mer, input id)
        // Reverse is to make it a min-heap
        let mut heap = std::collections::BinaryHeap::<(Reverse<LongKmer::<B>>, usize)>::new();
        for (i, input) in inputs.iter_mut().enumerate(){
            if let Some(kmer) = LongKmer::<B>::load(input)?{
                heap.push((Reverse(kmer), i));
            }
        }

        // Merge to out
        while let Some((Reverse(kmer), i)) = heap.pop(){
            kmer.serialize(&mut out)?;
            n_read[i] += 1;
            if let Some(kmer) = LongKmer::<B>::load(&mut inputs[i])?{
                heap.push((Reverse(kmer), i));
            }
        }

        Ok(())
    }

    
    // Consumes 'out' to ensure flushing
    pub fn serialize<W: std::io::Write>(&self, mut out: W) -> std::io::Result<usize>{
        let mut n_written = 0_usize;
        
        for kmer in self.kmers.iter(){
            n_written += kmer.serialize(&mut out)?;
        }

        Ok(n_written)
    }

    #[cfg(test)]
    pub fn load<R: std::io::Read>(input: &mut R) -> std::io::Result<Self>{

        let mut kmers = Vec::<LongKmer::<B>>::new();
        while let Some(kmer) = LongKmer::<B>::load(input)?{
            kmers.push(kmer);
        }
        Ok(Self{kmers})
    }
}

#[cfg(test)]
mod tests{

    use super::*;
    use std::io::Cursor;

    #[test]
    fn chunk_sorting(){

        let seqs = [b"AGGTAAA".to_vec(), b"ACAGGTAGGAAAGGAAAGT".to_vec()];
        let k = 4;

        let mut kmer_byte_chunks: Vec<Vec<u8>> = vec![];

        let unsorted_kmers: Vec<&[u8]> = seqs.iter().flat_map(|s| s.windows(k)).collect();
        for chunk in unsorted_kmers.chunks(3){ // Testing
            let iter = chunk.iter().map(|x| LongKmer::<2>::from_ascii(x).unwrap());
            let chunk = KmerChunk::new(iter); // Shadows the loop variable
            let chunk = chunk.sort(); // Shadows again
            let buf = Cursor::new(Vec::<u8>::new());
            let mut writer = std::io::BufWriter::new(buf);
            chunk.serialize(&mut writer).unwrap();
            kmer_byte_chunks.push(writer.into_inner().unwrap().into_inner());
        }

        let readers = kmer_byte_chunks.iter().map(|chunk| Cursor::new(chunk.as_slice())).collect::<Vec<Cursor<&[u8]>>>();
        let mut out_buf = Cursor::new(Vec::<u8>::new());
        SortedKmerChunk::<2>::merge_serialized(readers, &mut out_buf).unwrap();

        let mut merged_reader = Cursor::new(out_buf.into_inner());
        let kmers = SortedKmerChunk::load(&mut merged_reader).unwrap().into_vec();

        let mut true_kmers: Vec<LongKmer::<2>> = unsorted_kmers.iter().map(|x| LongKmer::<2>::from_ascii(x).unwrap()).collect();
        true_kmers.sort();
        assert_eq!(kmers, true_kmers);
    }
}
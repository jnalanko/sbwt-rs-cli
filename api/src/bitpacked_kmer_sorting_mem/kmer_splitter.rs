use crate::kmer::LongKmer;

use rayon::iter::IntoParallelRefIterator;
use rayon::iter::IntoParallelRefMutIterator;
use rayon::iter::ParallelIterator;
use rayon::prelude::ParallelSlice;
use rayon::prelude::ParallelSliceMut;

pub fn bitpack_rev_kmers<const B: usize, IN: crate::SeqStream + Send>(
    mut seqs: IN,
    k: usize,
    n_threads: usize,
    dedup_batches: bool,
) -> Vec<LongKmer<B>> {

    let producer_buf_size = 1_000_000_usize; // TODO: respect this

    // Wrap to scope to be able to borrow seqs for the producer thread even when it's not 'static.
    let mut kmers = std::thread::scope(|scope| {
        use crossbeam::crossbeam_channel::unbounded;
        let (parser_out, encoder_in) = unbounded();
        let (encoder_out, collector_in) = unbounded();

        // Create producer
        let producer_handle = scope.spawn(move || {
            let mut buf = Vec::<Box<[u8]>>::new();
            let mut current_total_buffer_size = 0_usize;
            
            while let Some(seq) = seqs.stream_next(){
                current_total_buffer_size += seq.len();
                let mut seq_copy = seq.to_owned();
                seq_copy.reverse(); // Reverse to get colex sorting
                buf.push(seq_copy.into_boxed_slice());
                if current_total_buffer_size > producer_buf_size {
                    let mut sendbuf = Vec::<Box<[u8]>>::new();
                    std::mem::swap(&mut sendbuf, &mut buf);
                    parser_out.send(sendbuf).unwrap();
                    current_total_buffer_size = 0;
                }
            }
            
            parser_out.send(buf).unwrap();
            drop(parser_out);
        });

        // Create encoders
        let mut encoder_handles = Vec::<std::thread::JoinHandle::<()>>::new();
        for _ in 0..n_threads{
            let receiver_clone = encoder_in.clone();
            let sender_clone = encoder_out.clone();
            encoder_handles.push(std::thread::spawn(move || {
                while let Ok(batch) = receiver_clone.recv(){
                    let mut encoded_batch = Vec::<LongKmer::<B>>::new();
                    for seq in batch{
                        for kmer in seq.windows(k){
                            match LongKmer::<B>::from_ascii(kmer) {
                                Ok(kmer) => {
                                    encoded_batch.push(kmer);
                                }
                                Err(crate::kmer::KmerEncodingError::InvalidNucleotide(_)) => (), // Ignore
                                Err(crate::kmer::KmerEncodingError::TooLong(_)) => panic!("k = {} is too long", k),
                            }        
                        }
                    }
                    if dedup_batches {
                        encoded_batch.sort_unstable();
                        encoded_batch.dedup();
                    }
                    sender_clone.send(encoded_batch).unwrap();
                }
            }));
        }

        let collector_handle = std::thread::spawn( move || {
            let mut kmers = Vec::<LongKmer::<B>>::new();
            while let Ok(batch) = collector_in.recv(){
                kmers.extend(batch);
            }
            kmers
        });

        producer_handle.join().unwrap(); // Wait for the producer to finish
        drop(encoder_in); // Close the channel
        for h in encoder_handles { // Wait for the encoders to finish
            h.join().unwrap();
        }
        drop(encoder_out); // Close the channel

        collector_handle.join().unwrap() // Wait for the collector to finish and return the kmers
    });

    kmers

}

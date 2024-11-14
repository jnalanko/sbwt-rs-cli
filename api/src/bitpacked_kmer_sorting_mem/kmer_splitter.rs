use crate::kmer::LongKmer;

use rayon::iter::IntoParallelRefIterator;
use rayon::iter::IntoParallelRefMutIterator;
use rayon::iter::ParallelIterator;
use rayon::prelude::ParallelSlice;
use rayon::prelude::ParallelSliceMut;

pub fn get_bitpacked_sorted_distinct_kmers<const B: usize, IN: crate::SeqStream + Send>(
    mut seqs: IN,
    k: usize,
    n_threads: usize,
    dedup_batches: bool,
) -> Vec<LongKmer<B>> {

    let bin_prefix_len = 3_usize; // If you update this you must update all the logic below
    assert!(k >= bin_prefix_len);
    let n_bins = (4_usize).pow(bin_prefix_len as u32); // 64
    let producer_buf_size = 1_000_000_usize; // TODO: respect this
    let encoder_bin_buf_size = 1_000_000_usize; // Deduplicate after this many k-mers in bin buffer. Todo: take as parameter.

    // Wrap to scope to be able to borrow seqs for the producer thread even when it's not 'static.
    let mut bins = std::thread::scope(|scope| {
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
                    let mut bin_buffers = vec![Vec::<LongKmer::<B>>::new(); n_bins];
                    for buf in bin_buffers.iter_mut(){
                        buf.reserve_exact(encoder_bin_buf_size);
                    }
                    for seq in batch{
                        for kmer in seq.windows(k){
                            match LongKmer::<B>::from_ascii(kmer) {
                                Ok(kmer) => {
                                    let bin_id = kmer.get_from_left(0) as usize * 16 + kmer.get_from_left(1) as usize * 4 + kmer.get_from_left(2) as usize; // Interpret nucleotides in base-4
                                    bin_buffers[bin_id].push(kmer);
                                    if bin_buffers[bin_id].len() == encoder_bin_buf_size{
                                        if dedup_batches{
                                            bin_buffers[bin_id].sort_unstable();
                                            bin_buffers[bin_id].dedup();
                                        }
                                        sender_clone.send(bin_buffers[bin_id].clone()).unwrap();
                                        bin_buffers[bin_id].clear();
                                    }
                                }
                                Err(crate::kmer::KmerEncodingError::InvalidNucleotide(_)) => (), // Ignore
                                Err(crate::kmer::KmerEncodingError::TooLong(_)) => panic!("k = {} is too long", k),
                            }        
                        }
                    }

                    // Send remaining buffers
                    for mut b in bin_buffers{
                        if dedup_batches{
                            b.sort_unstable();
                            b.dedup();
                        }
                        sender_clone.send(b).unwrap();
                    }
                }
            }));
        }

        let collector_handle = std::thread::spawn( move || {
            let mut bins = vec![Vec::<LongKmer::<B>>::new(); n_bins];
            while let Ok(batch) = collector_in.recv(){
                if !batch.is_empty() {
                    let bin_id = batch[0].get_from_left(0) as usize * 16 + batch[0].get_from_left(1) as usize * 4 + batch[0].get_from_left(2) as usize; // Intepret nucleotides in base-4
                    bins[bin_id].extend(batch);
                }
            }
            bins
        });

        producer_handle.join().unwrap(); // Wait for the producer to finish
        drop(encoder_in); // Close the channel
        for h in encoder_handles { // Wait for the encoders to finish
            h.join().unwrap();
        }
        drop(encoder_out); // Close the channel

        collector_handle.join().unwrap() // Wait for the collector to finish and return the bins

    });

    // Sort bins in parallel: todo: largest first

    log::info!("Sorting k-mer bins");
    bins.par_iter_mut().for_each(|bin| {
        if !bin.is_empty() {
            let label = &bin.first().unwrap().to_string()[0..2];
            log::info!("Sorting bin {} of size {}", label, bin.len());
            bin.sort_unstable();
            bin.dedup();
        } else {
            log::info!("Empty bin -> not sorting.");
        }
    });


    let mut bin_concat = Vec::<LongKmer::<B>>::new();

    log::info!("Concatenating k-mer bins");
    // Concat bins. Each of them is sorted, and they are bucketed by the first 3-mer in order, so the concatenation is sorted.
    for bin in bins {
        bin_concat.extend(bin.iter());
    }

    bin_concat

}

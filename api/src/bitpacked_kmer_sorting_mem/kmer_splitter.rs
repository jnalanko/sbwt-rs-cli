use std::sync::Arc;
use std::sync::Mutex;

use crate::kmer::LongKmer;

use rayon::iter::IntoParallelRefMutIterator;
use rayon::iter::ParallelIterator;

fn merge_sorted_unique_in_place<const B: usize>(a: &mut Vec<LongKmer<B>>, b: Vec<LongKmer<B>>) {
    // Step 1: Extend `a` with enough space
    let original_len = a.len();
    let total_capacity = original_len + b.len();
    a.resize(total_capacity, LongKmer::from_u64_data([0; B])); // Fill with placeholders

    // Step 2: Pointers for merge
    let mut i = original_len as isize - 1; // Last element of original `a`
    let mut j = b.len() as isize - 1;      // Last element of `b`
    let mut k = total_capacity as isize - 1; // End of resized `a`

    // Step 3: Merge from back to front
    while i >= 0 && j >= 0 {
        let va = a[i as usize];
        let vb = b[j as usize];
        if va > vb {
            a[k as usize] = va;
            i -= 1;
        } else if va < vb {
            a[k as usize] = vb;
            j -= 1;
        } else {
            // Equal values, keep one
            a[k as usize] = va;
            i -= 1;
            j -= 1;
        }
        k -= 1;
    }

    // Copy remaining elements
    while i >= 0 {
        a[k as usize] = a[i as usize];
        i -= 1;
        k -= 1;
    }
    while j >= 0 {
        a[k as usize] = b[j as usize];
        j -= 1;
        k -= 1;
    }

    // Step 4: Remove duplicates and shift left. Todo: no Option 
    let mut last: Option<LongKmer<B>> = None;
    let mut write = 0;
    for read in (k + 1) as usize..total_capacity {
        if Some(a[read]) != last {
            a[write] = a[read];
            last = Some(a[read]);
            write += 1;
        }
    }

    // Step 5: Truncate the vector
    a.truncate(write);
    a.shrink_to_fit();
}

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
    let per_thread_bin_buf_size = 10_000_usize;

    let mut shared_bin_buffers_vec = Vec::<Mutex::<Vec::<LongKmer::<B>>>>::new();
    for _ in 0..n_bins {
        let buf = Vec::<LongKmer::<B>>::new();
        let b = Mutex::new(buf);
        shared_bin_buffers_vec.push(b);
    };
    let shared_bin_buf_capacity = 1_000_000;
    let shared_bin_buffers = &shared_bin_buffers_vec; // This is shared with threads

    log::info!("Bitpacking and binning k-mers");
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
        let mut encoder_handles = Vec::<std::thread::ScopedJoinHandle::<()>>::new();

        for _ in 0..n_threads {
            let receiver_clone = encoder_in.clone();
            let sender_clone = encoder_out.clone();
            encoder_handles.push(scope.spawn(move || {
                while let Ok(batch) = receiver_clone.recv(){
                    let mut this_thread_bin_buffers = vec![Vec::<LongKmer::<B>>::new(); n_bins];
                    for seq in batch{
                        for kmer in seq.windows(k){
                            match LongKmer::<B>::from_ascii(kmer) {
                                Ok(kmer) => {
                                    let bin_id = kmer.get_from_left(0) as usize * 16 + kmer.get_from_left(1) as usize * 4 + kmer.get_from_left(2) as usize; // Interpret nucleotides in base-4
                                    this_thread_bin_buffers[bin_id].push(kmer);
                                    if this_thread_bin_buffers[bin_id].len() >= per_thread_bin_buf_size {
                                        // Move this local bin buffer to a shared buffer
                                        let shared_bin = &mut shared_bin_buffers[bin_id].lock().unwrap();
                                        shared_bin.extend(&this_thread_bin_buffers[bin_id]);
                                        this_thread_bin_buffers[bin_id].clear();
                                        if shared_bin.len() >= shared_bin_buf_capacity {
                                            // Flush shared bin to the collector thread
                                            if dedup_batches {
                                                log::debug!("Sorting batch of {} kmers", shared_bin.len());
                                                shared_bin.sort_unstable();
                                                shared_bin.dedup();
                                            }
                                            sender_clone.send(shared_bin.clone()).unwrap();
                                            shared_bin.clear();
                                        }
                                    }
                                }
                                Err(crate::kmer::KmerEncodingError::InvalidNucleotide(_)) => (), // Ignore
                                Err(crate::kmer::KmerEncodingError::TooLong(_)) => panic!("k = {} is too long", k),
                            }        
                        }
                    }

                    // Send remaining internal buffers of this thread
                    for mut b in this_thread_bin_buffers.into_iter() {
                        if dedup_batches{
                            log::debug!("Sorting remaining per-thread batch of {} kmers", b.len());
                            b.sort_unstable();
                            b.dedup();
                        }
                        sender_clone.send(b).unwrap();
                    }

                    // Send remaining shared batches
                    for sb in shared_bin_buffers.iter() {
                        let mut sb = sb.lock().unwrap();
                        if dedup_batches && !sb.is_empty(){
                            log::debug!("Sorting final shared batch of {} kmers", sb.len());
                            sb.sort_unstable();
                            sb.dedup();
                        }
                        sender_clone.send(sb.clone()).unwrap();
                        sb.clear();
                    }
                }
            }));
        }

        let collector_handle = std::thread::spawn( move || {
            let mut final_bins = vec![Vec::<LongKmer::<B>>::new(); n_bins];
            while let Ok(batch) = collector_in.recv(){
                if !batch.is_empty() {
                    let bin_id = batch[0].get_from_left(0) as usize * 16 + batch[0].get_from_left(1) as usize * 4 + batch[0].get_from_left(2) as usize; // Intepret nucleotides in base-4
                    final_bins[bin_id].extend(batch);
                }
            }
            final_bins
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
            let label = &bin.first().unwrap().to_string()[0..bin_prefix_len];
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

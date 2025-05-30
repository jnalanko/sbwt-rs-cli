use std::cmp::max;
use std::sync::Mutex;

use crate::kmer::LongKmer;

use rayon::iter::IntoParallelIterator;
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

// Rayon thread pool must be initialized before calling
pub fn get_bitpacked_sorted_distinct_kmers<const B: usize, IN: crate::SeqStream + Send>(
    mut seqs: IN,
    k: usize,
    n_threads: usize,
    dedup_batches: bool,
    approx_mem_gb: usize
) -> Vec<LongKmer<B>> {

    // There is one bin for each 3-mer (64 bins). If you change the 3 to something else,
    // you must update all the logic below, and also the buffer size calculation 
    // at the caller, and the log messages about buffer sizes, and possibly more.
    const BIN_PREFIX_LEN: usize = 3_usize; 
    assert!(k >= BIN_PREFIX_LEN);
    let n_bins = (4_usize).pow(BIN_PREFIX_LEN as u32); // 64

    // Calculating suitable buffer sizes
    // * producer_buf_capacity: This determines how many k-mers are pushed to the work queue
    //   at a time. This needs to be large enough to avoid parallel contention, but small enough
    //   to distribute work quickly and evenly. A good default is 2^20.
    // * thread_local_bin_buf_capacity. There is one thread-local buffer for each of the 64
    //   distinct 3-mers. Each of these buffers up to thread_local_bin_buf_capacity k-mers.
    //   So if thread_local_bin_buf_capacity is C_t, then the total space is: 
    //      C_t * n_threads * 64 * sizeof(LongKmer<B>>
    //      = C_t * n_threads * 64 * 8B
    //      = 512 * C_t * n_threads * B bytes
    //   The bigger the C_t, the less parallel contention there will be. If for example you have
    //   1GB memory available for this over 48 threads and k = 31 (B = 1), then you'll want to set 
    //   C_t to about 40,000.
    // * shared_bin_buf_capacity: There is one shared buffer for each of the 64 3-mer bins. Each
    //   of these contains up to shared_bin_buf_capacity k-mers. If `dedup_batches` is enabled, then
    //   the shared bin buffers are sorted and deduplicated before storing to memory. The larger the
    //   buffer size, the more duplicates we will find. If the capacity is C_s, then the total space
    //   for the shared buffers will be:
    //       64 * C_s * sizeof(LongKmer<B>) = 512 * C_s * B
    //   If dedup_batches is not enabled, or your data has little to no duplicates, this buffer does 
    //   not matter so much and you can just set it to e.g. 2^20. Otherwise, it's a good idea to put
    //   all your available extra memory here to deduplicate as much as possible to decrease the peak memory. 
    //   For example, if you have 512GB available for this with k = 31 (B = 1), set C_s to 512GB / 512 = 1 GB.

    let kmer_bytes = std::mem::size_of::<crate::kmer::LongKmer<B>>();
    let producer_buf_cap = (1_usize << 23) / kmer_bytes; // 8 MB of k-mers
    let thread_local_buf_caps = 1_usize << 16;
    let shared_buf_caps = if dedup_batches { 
        let bytes_remaining = approx_mem_gb as isize * (1 << 30) - (producer_buf_cap * kmer_bytes) as isize - (thread_local_buf_caps * n_threads * 64 * kmer_bytes) as isize;
        let bytes_remaining = max(bytes_remaining, 1_isize << 30) as usize; // Use at least 1 GB

        // The total memory in the shared buffers will be:
        // 64 * buf_cap * kmer_bytes 
        // Set this equal to bytes_remaining and solve for buf_cap:
        let cap = bytes_remaining / 64 / kmer_bytes;
        assert!(cap > 0); // Should be because bytes_remaining >= 1GB
        cap
    } else {
        1 << 30 // Whatever, does not matter since we do not deduplicate
    };

    // One bin per thread per 3-mer
    let thread_local_total_bytes = thread_local_buf_caps * n_threads * 64 * kmer_bytes; 
    let shared_total_bytes = shared_buf_caps * 64 * kmer_bytes; // One bin per 3-mer
    log::info!("Producer buffer capacity: {}", human_bytes::human_bytes((producer_buf_cap * kmer_bytes) as f64));
    log::info!("Thread-local buffer capacity: {} ({} total)", human_bytes::human_bytes((thread_local_buf_caps * kmer_bytes) as f64), human_bytes::human_bytes(thread_local_total_bytes as f64)); // 64 local bins per thread
    log::info!("Shared bin buffer capacity: {} ({} total)", human_bytes::human_bytes((shared_buf_caps * kmer_bytes) as f64), human_bytes::human_bytes(shared_total_bytes as f64));
    if producer_buf_cap*kmer_bytes + thread_local_total_bytes + shared_total_bytes > approx_mem_gb * (1_usize << 30) {
        log::warn!("Exceeding memory budget");
    }

    let mut shared_bin_buffers_vec = Vec::<Mutex::<Vec::<LongKmer::<B>>>>::new();
    for _ in 0..n_bins {
        let buf = Vec::<LongKmer::<B>>::new();
        let b = Mutex::new(buf);
        shared_bin_buffers_vec.push(b);
    };
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
                if current_total_buffer_size > producer_buf_cap {
                    let mut sendbuf = Vec::<Box<[u8]>>::new();
                    std::mem::swap(&mut sendbuf, &mut buf);
                    parser_out.send(sendbuf).unwrap();
                    current_total_buffer_size = 0;
                }
            }
            
            parser_out.send(buf).unwrap();

            log::info!("Producer thread: all work pushed to work queue");
            drop(parser_out);
        });

        // Create encoders
        let mut encoder_handles = Vec::<std::thread::ScopedJoinHandle::<()>>::new();

        for _ in 0..n_threads {
            let receiver_clone = encoder_in.clone();
            let sender_clone = encoder_out.clone();
            encoder_handles.push(scope.spawn(move || {
                let mut this_thread_bin_buffers = vec![Vec::<LongKmer::<B>>::new(); n_bins];
                while let Ok(batch) = receiver_clone.recv(){
                    for seq in batch{
                        for kmer in seq.windows(k){
                            match LongKmer::<B>::from_ascii(kmer) {
                                Ok(kmer) => {
                                    let bin_id = kmer.get_from_left(0) as usize * 16 + kmer.get_from_left(1) as usize * 4 + kmer.get_from_left(2) as usize; // Interpret nucleotides in base-4
                                    this_thread_bin_buffers[bin_id].push(kmer);
                                    if this_thread_bin_buffers[bin_id].len() >= thread_local_buf_caps {
                                        // Move this local bin buffer to a shared buffer
                                        let mut shared_bin = shared_bin_buffers[bin_id].lock().unwrap();
                                        shared_bin.extend(&this_thread_bin_buffers[bin_id]);
                                        this_thread_bin_buffers[bin_id].clear();
                                        if shared_bin.len() >= shared_buf_caps {
                                            // Flush shared bin to the collector thread
                                            if dedup_batches {
                                                let mut shared_bin_copy = shared_bin.clone();
                                                shared_bin.clear();
                                                drop(shared_bin); // Release the mutex and proceed to sort
                                                let len_before = shared_bin_copy.len();
                                                shared_bin_copy.sort_unstable();
                                                shared_bin_copy.dedup();
                                                shared_bin_copy.shrink_to_fit();
                                                let len_after = shared_bin_copy.len();
                                                log::debug!("Deduplicated batch of {} kmers ({:.2}% kept)", len_before, len_after as f64 / len_before as f64 * 100.0);
                                                sender_clone.send(shared_bin_copy).unwrap();
                                            } else {
                                                sender_clone.send(shared_bin.clone()).unwrap();
                                                shared_bin.clear();
                                            }
                                        }
                                    }
                                }
                                Err(crate::kmer::KmerEncodingError::InvalidNucleotide(_)) => (), // Ignore
                                Err(crate::kmer::KmerEncodingError::TooLong(_)) => panic!("k = {} is too long", k),
                            }        
                        }
                    }
                }

                // Send remaining internal buffers of this thread
                for mut b in this_thread_bin_buffers.into_iter() {
                    if dedup_batches{
                        log::debug!("Sorting remaining thread-local batch of {} kmers", b.len());
                        b.sort_unstable();
                        b.dedup();
                    }
                    sender_clone.send(b).unwrap();
                }

            }));
        }

        // Spawn a collector that reads from the encoders and pushes to final bins
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

        // Send remaining shared batches
        shared_bin_buffers.into_par_iter().for_each(|sb| {
            let mut sb = sb.lock().unwrap();
            if dedup_batches && !sb.is_empty(){
                log::debug!("Sorting remaining shared batch of {} kmers", sb.len());
                sb.sort_unstable();
                sb.dedup();
            }
            encoder_out.send(sb.clone()).unwrap(); // TODO: no clone.
            sb.clear();
            sb.shrink_to_fit();
        });

        drop(encoder_out); // Close the channel

        collector_handle.join().unwrap() // Wait for the collector to finish and return the bins

    });

    // Sort bins in parallel: todo: largest first

    log::info!("Sorting k-mer bins");
    bins.par_iter_mut().for_each(|bin| {
        if !bin.is_empty() {
            let label = &bin.first().unwrap().to_string()[0..BIN_PREFIX_LEN];
            log::info!("Sorting bin {} of size {}", label, bin.len());
            bin.sort_unstable();
            bin.dedup();
            bin.shrink_to_fit();
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

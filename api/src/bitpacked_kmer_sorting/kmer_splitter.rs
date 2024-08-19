use super::kmer::*;
use super::kmer_chunk::KmerChunk;
use crate::tempfile::{TempFile, TempFileManager};
use crate::util::DNA_ALPHABET;
use std::io::{BufWriter, Seek, Write};
use std::thread;

fn colex_sorted_binmers(bin_prefix_len: usize) -> Vec<Vec<u8>> {
    let mut binmers = Vec::<Vec<u8>>::new();
    for i in 0..(4_usize.pow(bin_prefix_len as u32)){
        let mut binmer = Vec::<u8>::new();
        let mut j = i;
        for _ in 0..bin_prefix_len{
            binmer.push(DNA_ALPHABET[j % 4]);
            j /= 4;
        }
        binmers.push(binmer);
    }
    binmers
}

pub fn split_to_bins<const B: usize, IN: crate::SeqStream + Send>(mut seqs: IN, k: usize, mem_gb: usize, n_threads: usize, dedup_batches: bool, temp_file_manager: &mut TempFileManager) -> Vec<TempFile>{

    // Suppose we have a memory budget of m bytes and t threads.
    // Suppose each k-mer takes s bytes and there are 64 bins.
    // Let b be the number of k-mers in each splitter thread bin buffer.
    // A splitter thread uses 64bs bytes
    // In total the splitter threads use 64bst threads.
    // So, we need:
    // 64bbt = m
    // b = m / (64bt)

    // Wrap to scope to be able to borrow seqs for the producer thread even when it's not 'static.
    std::thread::scope(|scope| {

        let bin_prefix_len = 3_usize; // If you update this you must update all the logic below
        let n_bins = (4_usize).pow(bin_prefix_len as u32); // 64
        let producer_buf_size = 1_000_000_usize; // TODO: respect this
        let encoder_bin_buf_size = mem_gb * (1_usize << 30) / ((n_bins * LongKmer::<B>::byte_size()) * n_threads);

        log::info!("Splitting k-mers into {} bins", n_bins);
        log::info!("Bin buffer size: {}", encoder_bin_buf_size);

        use crossbeam::crossbeam_channel::bounded;
        let (parser_out, encoder_in) = bounded(4);
        let (encoder_out, writer_in) = bounded(4);

        // Create producer
        let producer_handle = scope.spawn(move || {
            let mut buf = Vec::<Box<[u8]>>::new();
            let mut current_total_buffer_size = 0_usize;
            
            while let Some(seq) = seqs.stream_next(){
                current_total_buffer_size += seq.len();
                let mut seq_copy = seq.to_owned();
                seq_copy.reverse(); // Reverse to get colex sorting
                buf.push(seq_copy.into_boxed_slice());
                if current_total_buffer_size > producer_buf_size{
                    let mut sendbuf = Vec::<Box<[u8]>>::new();
                    std::mem::swap(&mut sendbuf, &mut buf);
                    parser_out.send(sendbuf).unwrap();
                    current_total_buffer_size = 0;
                }
            }
            
            parser_out.send(buf).unwrap();
            drop(parser_out);
        });

        // Create encoder-splitters
        let mut encoder_handles = Vec::<thread::JoinHandle::<()>>::new();
        for _ in 0..n_threads{
            let receiver_clone = encoder_in.clone();
            let sender_clone = encoder_out.clone();
            encoder_handles.push(std::thread::spawn(move || {
                let mut bin_buffers = vec![Vec::<LongKmer::<B>>::new(); n_bins];
                for buf in bin_buffers.iter_mut(){
                    buf.reserve_exact(encoder_bin_buf_size);
                }
                while let Ok(batch) = receiver_clone.recv(){
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
                                Err(KmerEncodingError::InvalidNucleotide(_)) => (), // Ignore
                                Err(KmerEncodingError::TooLong(_)) => panic!("k = {} is too long", k),
                            }        
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
            }));
        }

        // Create writers
        let mut bin_writers = 
            Vec::<std::io::BufWriter::<TempFile>>::new();

        for binmer in colex_sorted_binmers(bin_prefix_len) {
            let name_prefix = format!("sbwt-temp-{}-", String::from_utf8(binmer).unwrap());
            let f = temp_file_manager.create_new_file(&name_prefix, 8, ".bin");
            bin_writers.push(BufWriter::new(f));
        }


        let writer_handle = thread::spawn( move || {
            while let Ok(batch) = writer_in.recv(){
                if !batch.is_empty() {
                    let bin_id = batch[0].get_from_left(0) as usize * 16 + batch[0].get_from_left(1) as usize * 4 + batch[0].get_from_left(2) as usize; // Intepret nucleotides in base-4
                    let bin_file = &mut bin_writers[bin_id];
                    for kmer in batch{
                        kmer.serialize(bin_file).unwrap(); // Todo: write all at once
                    }
                }
            }
            bin_writers
        });

        producer_handle.join().unwrap();
        drop(encoder_in); // Close the channel
        for h in encoder_handles{
            h.join().unwrap();
        }
        drop(encoder_out); // Close the channel

        // Return the TempFiles
        let writers = writer_handle.join().unwrap();
        let mut writers: Vec<TempFile> = writers.into_iter().map(|w| w.into_inner().unwrap()).collect();
        for w in writers.iter_mut(){
            w.file.seek(std::io::SeekFrom::Start(0)).unwrap();
        }
        writers

    })
}

// Overwrite the files with sorted and deduplicates files. Returns back the files after overwriting.
pub fn par_sort_and_dedup_bin_files<const B: usize>(bin_files: Vec<TempFile>, mem_gb: usize, n_threads: usize) -> Vec<TempFile> {

    let filesizes = bin_files.iter().map(|f| f.path.metadata().unwrap().len() as usize).collect::<Vec<usize>>();
    let mut files_and_sizes = bin_files.into_iter().enumerate().map(|(i, f)| (f, filesizes[i], i)).collect::<Vec<(TempFile, usize, usize)>>();
        
    files_and_sizes.sort_by_key(|(_, size, _)| *size);

    let max_mem = mem_gb * (1_usize << 30);

    log::info!("Sorting k-mer bins");

    use crossbeam::unbounded;

    // A work queue
    let (queue_sender, queue_recvr) = unbounded::<(TempFile, usize, usize)>(); // File, size, index

    // A queue to notify the producer that a bin has been processed.
    // The usize in the channel is the size of the bin that was processed.
    let (producer_notify, producer_recv_notify) = unbounded::<usize>();

    // Wrap in mutex to share between threads
    let mut total_size_in_processing = 0_usize;

    // Start the producer
    let producer_handle = thread::spawn(move || {
        while !files_and_sizes.is_empty() {
            // Push as much work to the queue as possible
            while !files_and_sizes.is_empty(){    
                let s = files_and_sizes.last().unwrap().1; // Size
                if total_size_in_processing == 0 || total_size_in_processing + s <= max_mem {
                    let (f,s,i) = files_and_sizes.pop().unwrap();
                    queue_sender.send((f, s, i)).unwrap();
                    total_size_in_processing += s;
                } else {
                    break;
                }
            }

            let s_done = producer_recv_notify.recv().unwrap(); // Wait for a notification
            total_size_in_processing -= s_done;
        }

        // All files have been pushed to the channel
        drop(queue_sender); // Close the channel
    });

    let mut consumer_handles = Vec::<thread::JoinHandle<Vec::<(TempFile, usize)>>>::new();

    // Spawn consumers
    for _ in 0..n_threads{
        let recv_clone = queue_recvr.clone();
        let producer_notify = producer_notify.clone();

        consumer_handles.push(std::thread::spawn( move || {
            let mut processed_files = Vec::<(TempFile, usize)>::new(); // File, index
            while let Ok((mut f, s, i)) = recv_clone.recv(){
                // Using debug log level as a more verbose info level
                log::debug!("Sorting bin {} of size {}", f.path.display(), s);
                let mut reader = std::io::BufReader::new(&f.file);
                let chunk = KmerChunk::<B>::load(&mut reader).unwrap();
        
                let mut chunk = chunk.sort();
                chunk.dedup();

                // Overwrite the file and seek to start
                f.file.set_len(0).unwrap();
                f.file.seek(std::io::SeekFrom::Start(0)).unwrap();
                let chunk_out = std::io::BufWriter::new(&mut f);
                chunk.serialize(chunk_out).unwrap();
                f.flush().unwrap();
                f.file.seek(std::io::SeekFrom::Start(0)).unwrap();

                // Notify the producer that s bytes have been processed and
                // new work can possibly be pushed to the queue.
                let _ = producer_notify.send(s); // This may fail if the producer has already exited. That is ok.

                processed_files.push((f,i));
            }
            processed_files // Return to owner
        }));
    }

    let mut processed_files = Vec::<(TempFile, usize)>::new();
    producer_handle.join().unwrap();
    for h in consumer_handles{
        processed_files.extend(h.join().unwrap());
    }
    processed_files.sort_by(|(_, i1), (_, i2)| i1.cmp(i2));

    processed_files.into_iter().map(|(f,_)| f).collect() // Return to owner

}

// The original files are deleted
pub fn concat_files(infiles: Vec<TempFile>, out_writer: &mut impl std::io::Write){
    let mut bw = BufWriter::new(out_writer);
    for mut fp in infiles {
        let mut reader = std::io::BufReader::new(&mut fp.file);
        std::io::copy(&mut reader, &mut bw).unwrap();
        // fp is dropped here, which deletes the file
    }
    bw.flush().unwrap();
}

mod tests {
    #[test]
    fn test_colex_sorted_binmers(){
        let binmers = super::colex_sorted_binmers(2);
        let ans = vec![b"AA", b"CA", b"GA", b"TA", b"AC", b"CC", b"GC", b"TC", b"AG", b"CG", b"GG", b"TG", b"AT", b"CT", b"GT", b"TT"];
        assert_eq!(binmers, ans);
    }
}

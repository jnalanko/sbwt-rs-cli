use std::cmp::min;
use std::fs::File;
use std::io;
use std::io::stdout;
use std::io::BufWriter;
use std::io::Write;
use std::path::PathBuf;
use sbwt::dbg::Dbg;
use sbwt::*;
use sbwt::benchmark;

use jseqio::reader::*;
use sbwt::SbwtIndex;
use bitvec::prelude::*;

struct MySeqReader {
    inner: jseqio::reader::DynamicFastXReader,
}

impl sbwt::SeqStream for MySeqReader {
    fn stream_next(&mut self) -> Option<&[u8]> {
        self.inner.read_next().unwrap().map(|rec| rec.seq)
    }
}

// sbwt is taken as mutable because we need to build select support if all_at_once is true
fn dump_kmers<SS: SubsetSeq + Sync>(sbwt: &mut SbwtIndex<SS>, outfile: Option<&std::path::Path>, all_at_once: bool, include_dummies: bool, n_threads: usize) {
    let mut fileout = outfile.map(|f| BufWriter::new(File::create(f).unwrap()));
    let mut stdout = BufWriter::new(io::stdout());
    if all_at_once {
        log::info!("Reconstructing the k-mers");
        let concat = sbwt.reconstruct_padded_spectrum(n_threads);
        assert!(concat.len() % sbwt.k() == 0);

        log::info!("Dumping to stdout");
        let n_padded_kmers = concat.len() / sbwt.k();
        for i in 0..n_padded_kmers {
            let kmer = &concat[i*sbwt.k()..(i+1)*sbwt.k()];
            if include_dummies || kmer.iter().all(|&c| c != b'$') {
                match fileout {
                    Some(ref mut fileout) => {
                        fileout.write_all(kmer).unwrap();
                        fileout.write_all(b"\n").unwrap();
                    },
                    None => {
                        stdout.write_all(kmer).unwrap();
                        stdout.write_all(b"\n").unwrap();
                    }
                }
            }
        }
    } else {
        log::info!("Building select support");
        sbwt.build_select();

        log::info!("Dumping k-mers to stdout");
        let mut buf = vec![0u8; sbwt.k() + 1];
        for i in 0..sbwt.n_sets() {
            buf.clear();
            sbwt.push_kmer_to_vec(i, &mut buf);
            buf.push(b'\n');
            if include_dummies || buf.iter().all(|&c| c != b'$') {
                match fileout {
                    Some(ref mut fileout) => {
                        fileout.write_all(&buf).unwrap();
                    },
                    None => {
                        stdout.write_all(&buf).unwrap();
                    }
                }
            }
        }
    }

}

fn dump_kmers_command(matches: &clap::ArgMatches){

    let n_threads = *matches.get_one::<usize>("threads").unwrap();
    let indexfile = matches.get_one::<std::path::PathBuf>("index").unwrap();
    let include_dummies = matches.get_flag("include-dummy-kmers");
    let all_at_once = matches.get_flag("all-at-once");
    let cpp_format = matches.get_flag("load-cpp-format");
    let outfile = matches.get_one::<std::path::PathBuf>("output");

    // Don't care if there is LCS support or not
    let mut index_reader = std::io::BufReader::new(std::fs::File::open(indexfile).unwrap());
    let index = if cpp_format {
        load_from_cpp_plain_matrix_format(&mut index_reader).unwrap()
    } else {
        load_sbwt_index_variant(&mut index_reader).unwrap()
    };

    match index {
        SbwtIndexVariant::SubsetMatrix(mut sbwt) => {
            dump_kmers(&mut sbwt, outfile.map(|f| &**f), all_at_once, include_dummies, n_threads);
        }
    };

}

#[allow(non_snake_case)]
fn build_command(matches: &clap::ArgMatches){

    let infile = matches.get_one::<std::path::PathBuf>("input").unwrap();
    let out_prefix = matches.get_one::<std::path::PathBuf>("output-prefix").unwrap();
    let build_lcs = matches.get_flag("build-lcs");
    let k = *matches.get_one::<usize>("k").unwrap();
    let mem_gb = *matches.get_one::<usize>("mem-gb").unwrap();
    let n_threads = *matches.get_one::<usize>("threads").unwrap();
    let dedup_batches = matches.get_flag("dedup-batches");
    let add_revcomp = matches.get_flag("add-revcomp");
    let precalc_length = *matches.get_one::<usize>("prefix-precalc-length").unwrap();
    let temp_dir = matches.get_one::<std::path::PathBuf>("temp-dir").unwrap();

    let reader = MySeqReader{inner: jseqio::reader::DynamicFastXReader::from_file(infile).unwrap()};

    // Need to do this to be able to append .sbwt to the filename (PathBuf can only set extension, which replaces the existing one, meaning we can't stack extensions).
    let mut sbwt_outfile = out_prefix.clone().into_os_string().into_string().unwrap(); 

    sbwt_outfile.push_str(".sbwt");
    log::info!("Sbwt output file: {}", sbwt_outfile);
    let mut sbwt_out = std::io::BufWriter::new(std::fs::File::create(sbwt_outfile).unwrap()); // Open already here to fail early if problems
 
    log::info!("Building SBWT");
    let start_time = std::time::Instant::now();
    let algorithm = BitPackedKmerSorting::new().mem_gb(mem_gb).dedup_batches(dedup_batches).temp_dir(temp_dir);
    let (sbwt, lcs) = SbwtIndexBuilder::new().k(k).n_threads(n_threads).add_rev_comp(add_revcomp).algorithm(algorithm).build_lcs(build_lcs).precalc_length(precalc_length).run(reader);
    let end_time = std::time::Instant::now();
    log::info!("Construction finished in {:.2} seconds", (end_time - start_time).as_secs_f64());

    log::info!("Serializing");
    
    let sbwt_kmers = sbwt.n_kmers();
    let sbwt_bytes = write_sbwt_index_variant(&SbwtIndexVariant::SubsetMatrix(sbwt), &mut sbwt_out).unwrap();
    log::info!("Wrote sbwt index: {} bytes ({:.2} bits / k-mer)", sbwt_bytes, sbwt_bytes as f64 * 8.0 / sbwt_kmers as f64);

    if let Some(lcs) = lcs{
        let mut lcs_outfile = out_prefix.clone().into_os_string().into_string().unwrap(); // See comment on sbwt_outfile above
        lcs_outfile.push_str(".lcs");
        let mut lcs_out = std::io::BufWriter::new(std::fs::File::create(&lcs_outfile).unwrap());
        log::info!("Lcs output file: {}", lcs_outfile);

        let lcs_bytes = lcs.serialize(&mut lcs_out).unwrap();
        log::info!("Wrote lcs array: {} bytes ({:.2} bits / k-mer)", lcs_bytes, lcs_bytes as f64 * 8.0 / sbwt_kmers as f64);
    }
}

fn build_lcs_command(matches: &clap::ArgMatches) {
    let n_threads = *matches.get_one::<usize>("threads").unwrap();
    let indexfile = matches.get_one::<std::path::PathBuf>("index").unwrap();
    let cpp_format = matches.get_flag("load-cpp-format");
    let outfile = matches.get_one::<std::path::PathBuf>("output");
    let outfile = match outfile {
        Some(out) => out.clone(),
        None => {
            let mut f = indexfile.clone();
            f.set_extension("lcs");
            f
        }
    };

    log::info!("The LCS array will be saved to {}", outfile.display());

    let mut out = BufWriter::new(File::create(&outfile).unwrap());

    // Read sbwt
    log::info!("Loading the SBWT from {}", indexfile.display());
    let mut index_reader = std::io::BufReader::new(std::fs::File::open(indexfile).unwrap());
    let index = if cpp_format {
        load_from_cpp_plain_matrix_format(&mut index_reader).unwrap()
    } else {
        load_sbwt_index_variant(&mut index_reader).unwrap()
    };

    log::info!("Computing the LCS array");
    let lcs = match index {
        SbwtIndexVariant::SubsetMatrix(sbwt) => {
            LcsArray::from_sbwt(&sbwt, n_threads)
        }
    };

    log::info!("Writing to {}", outfile.display());
    let n_bytes = lcs.serialize(&mut out).unwrap();

    log::info!("{} bytes written", n_bytes);

}

fn report_lookup_results<W: Write>(out: &mut W, colex_ranks: &[Option<usize>], membership_only: bool){
    // Reporting
    if membership_only {
        let mut out_buf = Vec::<u8>::with_capacity(colex_ranks.len() + 1); // +1 for newline at the end
        for r in colex_ranks {
            if r.is_some() {
                out_buf.push(b'1');
            } else {
                out_buf.push(b'0');
            }
        }
        out_buf.push(b'\n');
        out.write_all(&out_buf).unwrap();
    } else {
        let mut number_ascii_buffer = [0u8; 32]; // Enough space for the ascii representation of a 64-bit integer
        for (i, r) in colex_ranks.iter().enumerate() {
            if i > 0 {
                out.write_all(b" ").unwrap();
            }
            let ascii_len = if let Some(r) = r {
                usize_to_ascii(*r, &mut number_ascii_buffer)
            } else {
                number_ascii_buffer[0] = b'-';
                number_ascii_buffer[1] = b'1';
                2
            };
            out.write_all(&number_ascii_buffer[0..ascii_len]).unwrap();
        }
        out.write_all(b"\n").unwrap();
    }
}

// Pushes results to output vector
#[allow(non_snake_case)]
fn non_streaming_search<SS: SubsetSeq>(sbwt: &SbwtIndex<SS>, query: &[u8], output: &mut Vec<Option<usize>>){
    for kmer in query.windows(sbwt.k()){
        let r = sbwt.search(kmer).map(|I| {
            assert!(I.len() == 1);
            I.start
        });
        output.push(r);
    }
}

// Pushes results to output vector
#[allow(non_snake_case)]
fn streaming_search<SS: SubsetSeq>(ss: &StreamingIndex<SbwtIndex<SS>, LcsArray>, k: usize, query: &[u8], output: &mut Vec<Option<usize>>){
    // Query using matching statistics
    let ms = ss.matching_statistics(query);
    for (len, I) in ms.iter().skip(k-1) {
        if *len == k {
            assert!(I.len() == 1);
            output.push(Some(I.start));
        } else {
            output.push(None);
        }
    }
}

// Writes to stdout if the outfile is not given
fn lookup_query<SS: SubsetSeq>(sbwt: &SbwtIndex<SS>, lcs: Option<LcsArray>, queryfile: &std::path::Path, outfile: Option<&std::path::Path>, membership_only: bool) {
    log::info!("k = {}, precalc length = {}, # kmers = {}, # sbwt sets = {}", sbwt.k(), sbwt.get_lookup_table().prefix_length, sbwt.n_kmers(), sbwt.n_sets());

    let mut query_reader = DynamicFastXReader::from_file(&queryfile).unwrap();
    let mut out = outfile.map(|f| std::io::BufWriter::new(std::fs::File::create(f).unwrap()));
    let mut stdout = stdout();

    let streaming_index = lcs.as_ref().map(|lcs| StreamingIndex::new(sbwt, lcs)); 

    let start_time = std::time::Instant::now();
    let mut n_query_kmers = 0_usize;
    let mut n_found = 0_usize;
    let mut colex_ranks = Vec::<Option<usize>>::new();
    while let Some(rec) = query_reader.read_next().unwrap(){
        let seq = rec.seq;
        colex_ranks.clear();

        if let Some(streaming_index) = &streaming_index{
            streaming_search(streaming_index, sbwt.k(), seq, &mut colex_ranks);
        } else {
            non_streaming_search(sbwt, seq, &mut colex_ranks);
        };
        n_query_kmers += colex_ranks.len();
        n_found += colex_ranks.iter().fold(0_usize, |count, r| count + r.is_some() as usize);

        match out {
            Some(ref mut out) => report_lookup_results(out, &colex_ranks, membership_only),
            None => report_lookup_results(&mut stdout, &colex_ranks, membership_only)
        }

        colex_ranks.clear();

    }
    let end_time = std::time::Instant::now();
    let elapsed = end_time - start_time;
    log::info!("Queried {} k-mers", n_query_kmers);
    log::info!("{:.2}% of queried k-mers found", n_found as f64 / n_query_kmers as f64 * 100.0);
    log::info!("Elapsed time: {:.2} seconds ({:.2} ns / k-mer)", elapsed.as_secs_f64(), elapsed.as_nanos() as f64 / n_query_kmers as f64);

}

#[allow(non_snake_case)]
fn lookup_query_command(matches: &clap::ArgMatches){
    let indexfile = matches.get_one::<std::path::PathBuf>("index").unwrap();
    let outfile = matches.get_one::<std::path::PathBuf>("output");
    let queryfile = matches.get_one::<std::path::PathBuf>("query").unwrap();
    let membership_only = matches.get_flag("membership-only");
    let cpp_format = matches.get_flag("load-cpp-format");

    // Read sbwt
    let mut index_reader = std::io::BufReader::new(std::fs::File::open(indexfile).unwrap());
    let index = if cpp_format {
        load_from_cpp_plain_matrix_format(&mut index_reader).unwrap()
    } else {
        load_sbwt_index_variant(&mut index_reader).unwrap()
    };

    // Load the lcs array, if available
    let mut lcsfile = indexfile.clone();
    lcsfile.set_extension("lcs"); // Replace .sbwt with .lcs
    let lcs = match std::fs::File::open(&lcsfile) {
        Ok(f) => {
            log::info!("Loading LCS array from file {}", lcsfile.display());
            let mut lcs_reader = std::io::BufReader::new(f);
            Some(LcsArray::load(&mut lcs_reader).unwrap())
        }
        Err(_) => {
            log::info!("No LCS array found at {} -> running without LCS support", lcsfile.display());
            None
        }
    };

    match index {
        SbwtIndexVariant::SubsetMatrix(sbwt) => {
            lookup_query(&sbwt, lcs, queryfile, outfile.map(|v| &**v), membership_only)
        }
    };

}

// Return the number of bytes written
fn usize_to_ascii(mut x: usize, bytes: &mut[u8; 32]) -> usize {
    if x == 0 {
        bytes[0] = b'0';
        return 1;
    }

    let mut byte_idx = 0_usize;
    while x > 0 {
        bytes[byte_idx] = b'0' + (x % 10) as u8;
        byte_idx += 1;
        x /= 10;
    }
    bytes[0..byte_idx].reverse();
    byte_idx
}

fn matching_statistics<SS: SubsetSeq>(sbwt: &SbwtIndex<SS>, lcs: &LcsArray, queryfile: &std::path::Path, outfile: Option<&std::path::Path>) {
    let mut query_reader = DynamicFastXReader::from_file(&queryfile).unwrap();
    let mut out = outfile.map(|f| std::io::BufWriter::new(std::fs::File::create(f).unwrap()));
    let mut stdout = stdout();

    let streaming_index = StreamingIndex::new(sbwt, lcs);

    let mut total_query_length = 0_usize;
    let mut out_buffer = Vec::<u8>::new();
    let mut number_ascii_buffer = [0u8; 32]; // Enough space for the ascii representation of a 64-bit integer
    let mut elapsed_without_io = std::time::Duration::from_nanos(0);
    let start_time = std::time::Instant::now();
    while let Some(rec) = query_reader.read_next().unwrap(){
        total_query_length += rec.seq.len();
        let ms_start_time = std::time::Instant::now();

        let ms = streaming_index.matching_statistics(rec.seq);

        let ms_end_time = std::time::Instant::now();
        elapsed_without_io += ms_end_time - ms_start_time;

        // Write out in CSV format
        out_buffer.clear();
        for (i, (len, _)) in ms.iter().enumerate(){
            if i > 0 {
                out_buffer.push(b',');
            }
            let ascii_length = usize_to_ascii(*len, &mut number_ascii_buffer);
            out_buffer.extend(&number_ascii_buffer[0..ascii_length]);
        }
        out_buffer.push(b'\n');

        match out {
            Some(ref mut out) => out.write_all(out_buffer.as_slice()).unwrap(),
            None => stdout.write_all(out_buffer.as_slice()).unwrap(),
        }
    }

    let total_elapsed = std::time::Instant::now() - start_time;
    log::info!("Total query length: {} nucleotides", total_query_length);
    log::info!("Elapsed time: {:.2} seconds ({:.2} ns / nucleotide)", total_elapsed.as_secs_f64(), total_elapsed.as_nanos() as f64 / total_query_length as f64);
    log::info!("Elapsed time excluding I/O: {:.2} seconds ({:.2} ns / nucleotide)", elapsed_without_io.as_secs_f64(), elapsed_without_io.as_nanos() as f64 / total_query_length as f64);

}

fn matching_statistics_command(matches: &clap::ArgMatches){
    let indexfile = matches.get_one::<std::path::PathBuf>("index").unwrap();
    let outfile = matches.get_one::<std::path::PathBuf>("output");
    let queryfile = matches.get_one::<std::path::PathBuf>("query").unwrap();
    let cpp_format = matches.get_flag("load-cpp-format");

    // Read sbwt
    log::info!("Loading sbwt index from file {}", indexfile.display());
    let mut index_reader = std::io::BufReader::new(std::fs::File::open(indexfile).unwrap());
    let index = if cpp_format {
        load_from_cpp_plain_matrix_format(&mut index_reader).unwrap()
    } else {
        load_sbwt_index_variant(&mut index_reader).unwrap()
    };

    // Load the lcs array
    let mut lcsfile = indexfile.clone();
    lcsfile.set_extension("lcs"); // Replace .sbwt with .lcs

    log::info!("Loading LCS array from file {}", lcsfile.display());
    let lcs = match std::fs::File::open(&lcsfile) {
        Ok(f) => {
            let mut lcs_reader = std::io::BufReader::new(f);
            LcsArray::load(&mut lcs_reader).unwrap()
        }
        Err(_) => {
            log::error!("No LCS array found at {}", lcsfile.display());
            log::error!("LCS array required for matching statistics (--build-lcs in index construction)");
            return
        }
    };

    match index {
        SbwtIndexVariant::SubsetMatrix(sbwt) => {
            matching_statistics(&sbwt, &lcs, queryfile, outfile.map(|v| &**v))
        }
    };

}

fn benchmark<SS: SubsetSeq, W: Write>(sbwt: SbwtIndex<SS>, lcs: Option<LcsArray>, json_out: Option<&mut W>) {
    log::info!("benchmarking index with k = {}, precalc length = {}, # kmers = {}, # sbwt sets = {}", sbwt.k(), sbwt.get_lookup_table().prefix_length, sbwt.n_kmers(), sbwt.n_sets());
    let results = benchmark::benchmark_all(sbwt, lcs);

    // Write results in json if requested
    if let Some(out) = json_out {
        let json = serde_json::to_string(&results).unwrap();
        out.write_all(json.as_bytes()).unwrap();
        out.write_all(b"\n").unwrap();
    };
}

fn benchmark_command(matches: &clap::ArgMatches) {
    let indexfile = matches.get_one::<PathBuf>("index").unwrap();
    let precalc_length = matches.get_one::<usize>("prefix-precalc-length");
    let cpp_format = matches.get_flag("load-cpp-format");
    let json_outfile = matches.get_one::<PathBuf>("json-out");

    // Read sbwt
    let mut index_reader = std::io::BufReader::new(std::fs::File::open(indexfile).unwrap());
    let index = if cpp_format {
        load_from_cpp_plain_matrix_format(&mut index_reader).unwrap()
    } else {
        load_sbwt_index_variant(&mut index_reader).unwrap()
    };

    // Load the lcs array, if available
    let mut lcsfile = indexfile.clone();
    lcsfile.set_extension("lcs"); // Replace .sbwt with .lcs
    let lcs = match std::fs::File::open(&lcsfile) {
        Ok(f) => {
            log::info!("Loading LCS array from file {}", lcsfile.display());
            let mut lcs_reader = std::io::BufReader::new(f);
            Some(LcsArray::load(&mut lcs_reader).unwrap())
        }
        Err(_) => {
            log::info!("No LCS array found at {} -> running without LCS support", lcsfile.display());
            None
        }
    };

    let mut json_out = json_outfile.map(|path| File::create(path).unwrap());

    match index {
        SbwtIndexVariant::SubsetMatrix(mut sbwt) => {
            if let Some(p) = precalc_length {
                log::info!("Rebuilding the precalc table with prefix length {}", p);
                let new_table = PrefixLookupTable::new(&sbwt, *p);
                sbwt.set_lookup_table(new_table);
            }

            benchmark(sbwt, lcs, json_out.as_mut());
        }
    };

}

fn dump_unitigs<SS: SubsetSeq + Send + Sync>(sbwt: &mut SbwtIndex<SS>, outfile: Option<&std::path::Path>, lcs: &Option<LcsArray>, n_threads: usize) {
    let stdout = std::io::stdout();
    let stdout = std::io::BufWriter::new(stdout);
    let mut fileout = outfile.map(|f| BufWriter::new(File::create(f).unwrap()));

    log::info!("Preparing the de Bruijn graph");
    sbwt.build_select();

    let dbg = Dbg::new(sbwt, lcs.as_ref(), n_threads);

    log::info!("Dumping unitigs");
    match fileout {
        Some(ref mut fileout) => dbg.parallel_export_unitigs(fileout, n_threads),
        None => dbg.parallel_export_unitigs(stdout, n_threads),
    }
}

fn dump_unitigs_command(matches: &clap::ArgMatches) {

    let n_threads = *matches.get_one::<usize>("threads").unwrap();
    let indexfile = matches.get_one::<std::path::PathBuf>("index").unwrap();
    let cpp_format = matches.get_flag("load-cpp-format");
    let outfile = matches.get_one::<std::path::PathBuf>("output");

    // Read sbwt
    let mut index_reader = std::io::BufReader::new(std::fs::File::open(indexfile).unwrap());
    let index = if cpp_format {
        load_from_cpp_plain_matrix_format(&mut index_reader).unwrap()
    } else {
        load_sbwt_index_variant(&mut index_reader).unwrap()
    };

    // Load the lcs array, if available
    let mut lcsfile = indexfile.clone();
    lcsfile.set_extension("lcs"); // Replace .sbwt with .lcs
    let lcs = match std::fs::File::open(&lcsfile) {
        Ok(f) => {
            log::info!("Loading LCS array from file {}", lcsfile.display());
            let mut lcs_reader = std::io::BufReader::new(f);
            Some(LcsArray::load(&mut lcs_reader).unwrap())
        }
        Err(_) => {
            log::info!("No LCS array found at {} -> Some extra computation needed to initialize the DBG", lcsfile.display());
            None
        }
    };

    match index {
        SbwtIndexVariant::SubsetMatrix(mut sbwt) => {
            dump_unitigs(&mut sbwt, outfile.map(|f| &**f), &lcs, n_threads);
        }
    };

}

// Assumes there is at least one 1-bit
fn leading_zeros(s: &BitSlice) -> usize {
    s.first_one().unwrap()
}

fn refine_segmentation(s1: BitVec, s2: BitVec, chars1: &[u8], chars2: &[u8]) -> (BitVec, BitVec) {
    let mut s1_i = 0_usize; // Index in s1
    let mut s2_i = 0_usize; // Index in s2

    let mut c1_i = 0_usize; // Index in chars1
    let mut c2_i = 0_usize; // Index in chars2

    let mut out1 = bitvec![];
    let mut out2 = bitvec![];

    while s1_i < s1.len() {
        assert!(s2_i < s2.len());

        let len1 = leading_zeros(&s1[s1_i..]);
        let len2 = leading_zeros(&s2[s2_i..]);

        let c1_end = c1_i + len1; // One past the end
        let c2_end = c2_i + len2; // One past the end

        while c1_i < c1_end || c2_i < c2_end {
            let c1 = if c1_i == c1_end { u8::MAX } else { chars1[c1_i] };
            let c2 = if c2_i == c2_end { u8::MAX } else { chars2[c2_i] };
            let c = min(c1,c2);

            while c1_i < c1_end && chars1[c1_i] == c {
                c1_i += 1;
                out1.push(false); // Unary bit
            }

            while c2_i < c2_end && chars2[c2_i] == c {
                c2_i += 1;
                out2.push(false); // Unary bit
            }

            // Terminate unary representations
            out1.push(true);
            out2.push(true);

        }

        assert_eq!(c1_i, c1_end);
        assert_eq!(c2_i, c2_end);

        s1_i += len1 + 1;
        s2_i += len2 + 1;

    }

    assert_eq!(s1_i, s1.len());
    assert_eq!(s2_i, s2.len());

    (out1, out2)
}

fn jaccard_command(matches: &clap::ArgMatches) {

    let n_threads = *matches.get_one::<usize>("threads").unwrap();
    let sbwt1_path = matches.get_one::<std::path::PathBuf>("sbwt1").unwrap();
    let sbwt2_path = matches.get_one::<std::path::PathBuf>("sbwt2").unwrap();
    let cpp_format = matches.get_flag("load-cpp-format");

    // Read sbwts
    let mut index1_reader = std::io::BufReader::new(std::fs::File::open(sbwt1_path).unwrap());
    let index1 = if cpp_format {
        load_from_cpp_plain_matrix_format(&mut index1_reader).unwrap()
    } else {
        load_sbwt_index_variant(&mut index1_reader).unwrap()
    };

    let mut index2_reader = std::io::BufReader::new(std::fs::File::open(sbwt2_path).unwrap());
    let index2 = if cpp_format {
        load_from_cpp_plain_matrix_format(&mut index2_reader).unwrap()
    } else {
        load_sbwt_index_variant(&mut index2_reader).unwrap()
    };

    let SbwtIndexVariant::SubsetMatrix(index1) = index1;
    let SbwtIndexVariant::SubsetMatrix(index2) = index2;

    let k = index1.k();
    assert_eq!(k, index2.k());

    // We invert the SBWTs column by column and maintain ranges
    // in both SBWTs that are so far qual. The ranges are in increasing
    // order and the partition the SBWTs, to it's enough to just store their
    // sizes in order. The sizes are stored in concatenated unary representations.
    // Empty ranges are allowed.

    // Initialize unary concatenations with empty ranges
    let mut s1 = bitvec![0; index1.n_sets()];
    let mut s2 = bitvec![0; index2.n_sets()];
    s1.push(true);
    s2.push(true);

    let mut chars1 = index1.build_last_column();
    let mut chars2 = index2.build_last_column();

    let mut temp_char_buf_1 = vec![0u8; chars1.len()];
    let mut temp_char_buf_2 = vec![0u8; chars2.len()];

    for round in 0..k {
        log::info!("Round {}/{}", round+1, k);

        let new_arrays = refine_segmentation(s1, s2, &chars1, &chars2);
        (s1, s2) = new_arrays;

        if round != k-1 {
            index1.push_all_labels_forward(&chars1, &mut temp_char_buf_1, n_threads);
            chars1.copy_from_slice(&temp_char_buf_1);

            index2.push_all_labels_forward(&chars2, &mut temp_char_buf_2, n_threads);
            chars2.copy_from_slice(&temp_char_buf_2);
        }
    }

    let mut intersection_size = 0_usize;
    let mut union_size = 0_usize;
    let mut s1_i = 0_usize; // Index in s1
    let mut s2_i = 0_usize; // Index in s2
    let mut c1_i = 0_usize; // Index in chars1
    let mut c2_i = 0_usize; // Index in chars2

    while s1_i < s1.len() {
        assert!(s2_i < s2.len());

        let len1 = leading_zeros(&s1[s1_i..]);
        let len2 = leading_zeros(&s2[s2_i..]);
        assert!(len1 <= 1); // This is the colex range of a k-mer, so it should be empty or singleton
        assert!(len2 <= 1); // Same as above
        assert!(len1 + len2 > 0); // Should not be empty interval pairs

        let is_dummy = (len1 > 0 && chars1[c1_i] == b'$') || (len2 > 0 && chars2[c2_i] == b'$');
        if !is_dummy {
            if len1 > 0 && len2 > 0 {
                intersection_size += 1;
            }
            union_size += 1;
        }

        s1_i += len1 + 1; // Length of the unary number we just parsed
        s2_i += len2 + 1; // Same as above

        c1_i += len1;
        c2_i += len2;
    }
    assert_eq!(s1_i, s1.len());
    assert_eq!(s2_i, s2.len());

    let jaccard = intersection_size as f64 / union_size as f64;
    println!("Intersection size: {}", intersection_size);
    println!("Union size: {}", union_size);
    println!("Jaccard index: {}", jaccard);
}

fn main() {

    let cli = clap::Command::new("sbwt")
        .about("Command line tools using for the sbwt library.")
        .arg_required_else_help(true)
        .arg(clap::Arg::new("threads")
            .help("Number of threads to use for parallelized operations")
            .long("threads")
            .short('t')
            .default_value("4")
            .global(true)
            .value_parser(clap::value_parser!(usize))
        )
        .arg(clap::Arg::new("verbose")
            .help("Print more information when running.")
            .short('v')
            .long("verbose")
            .global(true)
            .action(clap::ArgAction::SetTrue)
        )
        .subcommand(clap::Command::new("build")
            .about("Build the SBWT and possibly the LCS array")
            .arg_required_else_help(true)
            .arg(clap::Arg::new("input")
                .help("Input fasta or fastq sequence file")
                .short('i')
                .long("input")
                .required(true)
                .value_parser(clap::value_parser!(std::path::PathBuf))
            )
            .arg(clap::Arg::new("output-prefix")
                .help("Prefix for the output filenames. Writes to file [prefix].sbwt, and also [prefix].lcs if --build-lcs is given.")
                .short('o')
                .long("output-prefix")
                .required(true)
                .value_parser(clap::value_parser!(std::path::PathBuf))
            )
            .arg(clap::Arg::new("temp-dir")
                .help("Directory for temporary files created during construction.")
                .long("temp-dir")
                .default_value(".")
                .value_parser(clap::value_parser!(std::path::PathBuf))
            )
            .arg(clap::Arg::new("build-lcs")
                .help("Also build the LCS array (costs about log(k) bits per SBWT node)")
                .short('l')
                .long("build-lcs")
                .action(clap::ArgAction::SetTrue)
            )
            .arg(clap::Arg::new("k")
                .help("k-mer length")
                .short('k')
                .required(true)
                .value_parser(clap::value_parser!(usize))
            )
            .arg(clap::Arg::new("mem-gb")
                .help("An approximate memory budget in GB. The actual memory usage may be higher.")
                .short('m')
                .long("mem-gb")
                .required(false)
                .value_parser(clap::value_parser!(usize))
                .default_value("8")
            )
            .arg(clap::Arg::new("dedup-batches")
                .help("Slows down the construction, but saves disk if the data is very repetitive.")
                .long("dedup-batches")
                .short('d')
                .action(clap::ArgAction::SetTrue)
            )
            .arg(clap::Arg::new("add-revcomp")
                .help("Add reverse reverse complements of all k-mers to the index.")
                .long("add-revcomp")
                .short('r')
                .action(clap::ArgAction::SetTrue))
            .arg(clap::Arg::new("prefix-precalc-length")
                .help("Prefix precalc lookup table prefix length p. Larger p speeds up queries, but the table takes 4^(p+2) bytes of space.")
                .long("prefix-precalc")
                .short('p')
                .default_value("8")
                .value_parser(clap::value_parser!(usize)))
        )
        .subcommand(clap::Command::new("build-lcs")
            .arg_required_else_help(true)
            .about("Build the LCS array given an existing SBWT. Note that it's also possible build the SBWT and the LCS at the same time with the build-command.")
            .arg(clap::Arg::new("index")
                .help("SBWT index file")
                .short('i')
                .long("index")
                .required(true)
                .value_parser(clap::value_parser!(std::path::PathBuf))
            )
            .arg(clap::Arg::new("load-cpp-format")
                .help("Load the index from the format used by the C++ API (only supports plain-matrix).")
                .long("load-cpp-format")
                .action(clap::ArgAction::SetTrue)
            )
            .arg(clap::Arg::new("output")
                .help("The output filename. By default if the input is at index.sbwt, the output will be at index.lcs.")
                .long("output")
                .short('o')
                .value_parser(clap::value_parser!(std::path::PathBuf))
            )
        )
        .subcommand(clap::Command::new("benchmark")
            .about("Benchmark queries on the index")
            .arg_required_else_help(true)
            .arg(clap::Arg::new("index")
                .help("SBWT index file")
                .short('i')
                .long("index")
                .required(true)
                .value_parser(clap::value_parser!(std::path::PathBuf))
            )
            .arg(clap::Arg::new("load-cpp-format")
                .help("Load the index from the format used by the C++ API (only supports plain-matrix).")
                .long("load-cpp-format")
                .action(clap::ArgAction::SetTrue)
            )
            .arg(clap::Arg::new("prefix-precalc-length")
                .help("Disregard the precalc table that is stored in the index, and compute a new table with this prefix length." )
                .long("prefix-precalc")
                .short('p')
                .value_parser(clap::value_parser!(usize))
            )
            .arg(clap::Arg::new("json-out")
                .help("The results are written to this file in JSON format.")
                .short('j')
                .long("json-out")
                .value_parser(clap::value_parser!(std::path::PathBuf))
            )
        )
        .subcommand(clap::Command::new("lookup")
            .about("Look up the colex ranks of all k-mers. Printed as space-separated ascii integers, one line per query sequence. If a k-mer does not exist in the index, the colex rank is -1.")
            .arg_required_else_help(true)
            .arg(clap::Arg::new("index")
                .help("SBWT index file")
                .short('i')
                .long("index")
                .required(true)
                .value_parser(clap::value_parser!(std::path::PathBuf))
            )
            .arg(clap::Arg::new("query")
                .help("Query sequences in FASTA or FASTQ format, possibly gzipped")
                .short('q')
                .long("query")
                .required(true)
                .value_parser(clap::value_parser!(std::path::PathBuf))
            )
            .arg(clap::Arg::new("output")
                .help("Output text file. If not given, prints to stdout.")
                .short('o')
                .long("output")
                .value_parser(clap::value_parser!(std::path::PathBuf))
            )
            .arg(clap::Arg::new("membership-only")
                .help("Instead of colex ranks, prints an ascii bitvector B for each query sequence, such that B[i] = 1 iff the i-th k-mer is present in the index.")
                .short('b')
                .long("membership-only")
                .action(clap::ArgAction::SetTrue)
            )
            .arg(clap::Arg::new("load-cpp-format")
                .help("Load the index from the format used by the C++ API (only supports plain-matrix).")
                .long("load-cpp-format")
                .action(clap::ArgAction::SetTrue)
            )
        )
        .subcommand(clap::Command::new("matching-statistics")
            .arg_required_else_help(true)
            .arg(clap::Arg::new("index")
                .help("SBWT index file")
                .short('i')
                .long("index")
                .required(true)
                .value_parser(clap::value_parser!(std::path::PathBuf))
            )
            .arg(clap::Arg::new("query")
                .help("Query sequences in FASTA or FASTQ format, possibly gzipped")
                .short('q')
                .long("query")
                .required(true)
                .value_parser(clap::value_parser!(std::path::PathBuf))
            )
            .arg(clap::Arg::new("output")
                .help("Output text file. Prints to stdout if not given.")
                .short('o')
                .long("output")
                .value_parser(clap::value_parser!(std::path::PathBuf))
            )
            .arg(clap::Arg::new("load-cpp-format")
                .help("Load the index from the format used by the C++ API (only supports plain-matrix).")
                .long("load-cpp-format")
                .action(clap::ArgAction::SetTrue)
            )
        )
        .subcommand(clap::Command::new("jaccard")
            .arg_required_else_help(true)
            .about("Compute the Jaccard index of the k-mer sets in two given sbwt index files.")
            .arg(clap::Arg::new("sbwt1")
                .help("The first SBWT index file")
                .required(true)
                .value_parser(clap::value_parser!(std::path::PathBuf))
            )
            .arg(clap::Arg::new("sbwt2")
                .help("The second SBWT index file")
                .required(true)
                .value_parser(clap::value_parser!(std::path::PathBuf))
            )
            .arg(clap::Arg::new("load-cpp-format")
                .help("Load the index from the format used by the C++ API (only supports plain-matrix).")
                .long("load-cpp-format")
                .action(clap::ArgAction::SetTrue)
            )
        )
        .subcommand(clap::Command::new("dump-kmers")
            .about("Prints all k-mer strings in colex order, one per line.")
            .arg_required_else_help(true)
            .arg(clap::Arg::new("index")
                .short('i')
                .long("index")
                .required(true)
                .value_parser(clap::value_parser!(std::path::PathBuf))
            )
            .arg(clap::Arg::new("output")
                .help("Output text file. Prints to stdout if not given.")
                .short('o')
                .long("output")
                .value_parser(clap::value_parser!(std::path::PathBuf))
            )
            .arg(clap::Arg::new("include-dummy-kmers")
                .help("Include the dollar-padded dummy k-mers required for the SBWT")
                .short('d')
                .long("include-dummy-kmers")
                .action(clap::ArgAction::SetTrue)
            )
            .arg(clap::Arg::new("all-at-once")
                .help("Runs a fast algorithm that that reconstructs everything into memory at once taking n*(k+2) bytes of space")
                .short('a')
                .long("all-at-once")
                .action(clap::ArgAction::SetTrue)
            )
            .arg(clap::Arg::new("load-cpp-format")
                .help("Load the index from the format used by the C++ API (only supports plain-matrix).")
                .long("load-cpp-format")
                .action(clap::ArgAction::SetTrue)
            )
        )
        .subcommand(clap::Command::new("dump-unitigs")
            .about("Print all unitigs in arbitrary order in FASTA format")
            .arg_required_else_help(true)
            .arg(clap::Arg::new("index")
                .short('i')
                .long("index")
                .required(true)
                .value_parser(clap::value_parser!(std::path::PathBuf))
            )
            .arg(clap::Arg::new("output")
                .help("Output text file. Prints to stdout if not given.")
                .short('o')
                .long("output")
                .value_parser(clap::value_parser!(std::path::PathBuf))
            )
            .arg(clap::Arg::new("load-cpp-format")
                .help("Load the index from the format used by the C++ API (only supports plain-matrix).")
                .long("load-cpp-format")
                .action(clap::ArgAction::SetTrue)
            )
        )
        ;

    let matches = cli.get_matches();

    // Initialize logging 
    let mut builder = env_logger::builder();
    if matches.get_flag("verbose") {
        builder.filter_level(log::LevelFilter::Debug);
    } else {
        builder.filter_level(log::LevelFilter::Info);
    };

    // Renaming DEBUG level to INFO. We're using the DEBUG level as
    // a more verbose INFO level. Unfortunately the log crate interface does not support
    // multiple kinds of INFO-level messages, so we have this workaround.
    builder.format(|buf, record| {
        
        let style = match record.level() {
            log::Level::Debug => buf.default_level_style(log::Level::Info), // Use INFO level styling for DEBUG messages
            _ => buf.default_level_style(record.level()),
        };

        let level_string = match record.level() {
            log::Level::Debug => style.value("INFO"),
            _ => style.value(record.level().as_str()),
        };

        let module = match record.module_path() {
            None => "unknown",
            Some(x) => x,
        };

        let time = chrono::Local::now().format("%Y-%m-%d %H:%M:%S%.3f");

        writeln!(buf, "[{}] [{}] [{}]: {}", level_string, time, module, record.args())
    });
    builder.init();

    log::info!("sbwt-rs-cli version {}, api version {}", env!("CARGO_PKG_VERSION"), sbwt::CARGO_API_VERSION);

    match matches.subcommand(){
        Some(("build", sub_matches)) => build_command(sub_matches),
        Some(("build-lcs", sub_matches)) => build_lcs_command(sub_matches),
        Some(("lookup", sub_matches)) => lookup_query_command(sub_matches),
        Some(("matching-statistics", sub_matches)) => matching_statistics_command(sub_matches),
        Some(("jaccard", sub_matches)) => jaccard_command(sub_matches),
        Some(("benchmark", sub_matches)) => benchmark_command(sub_matches),
        Some(("dump-kmers", sub_matches)) => dump_kmers_command(sub_matches),
        Some(("dump-unitigs", sub_matches)) => dump_unitigs_command(sub_matches),
        _ => unreachable!(),
    }
    
}

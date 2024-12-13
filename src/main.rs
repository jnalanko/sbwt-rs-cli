use std::fs::File;
use std::io;
use std::io::stdout;
use std::io::BufWriter;
use std::io::Write;
use sbwt::dbg::Dbg;
use sbwt::*;
use sbwt::benchmark;

use jseqio::reader::*;
use sbwt::SbwtIndex;

struct MySeqReader {
    inner: jseqio::reader::DynamicFastXReader,
}

impl sbwt::SeqStream for MySeqReader {
    fn stream_next(&mut self) -> Option<&[u8]> {
        self.inner.read_next().unwrap().map(|rec| rec.seq)
    }
}

// sbwt is taken as mutable because we need to build select support if all_at_once is true
fn dump_kmers<SS: SubsetSeq + Sync>(sbwt: &mut SbwtIndex<SS>, all_at_once: bool, include_dummies: bool, n_threads: usize) {
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
                stdout.write_all(kmer).unwrap();
                stdout.write_all(b"\n").unwrap();
            }
        }
    } else {
        log::info!("Building select support");
        sbwt.build_select();

        log::info!("Dumping k-mers to stdout");
        let mut buf = vec![0u8; sbwt.k()];
        for i in 0..sbwt.n_sets() {
            buf.clear();
            sbwt.push_kmer_to_vec(i, &mut buf);
            if include_dummies || buf.iter().all(|&c| c != b'$') {
                stdout.write_all(&buf).unwrap();
                stdout.write_all(b"\n").unwrap();
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

    // Don't care if there is LCS support or not
    let mut index_reader = std::io::BufReader::new(std::fs::File::open(indexfile).unwrap());
    let index = if cpp_format {
        load_from_cpp_plain_matrix_format(&mut index_reader).unwrap()
    } else {
        load_sbwt_index_variant(&mut index_reader).unwrap()
    };

    match index {
        SbwtIndexVariant::SubsetMatrix(mut sbwt) => {
            dump_kmers(&mut sbwt, all_at_once, include_dummies, n_threads);
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

fn matching_statistics<SS: SubsetSeq>(sbwt: &SbwtIndex<SS>, lcs: &LcsArray, queryfile: &std::path::Path, outfile: &std::path::Path) {
    let mut query_reader = DynamicFastXReader::from_file(&queryfile).unwrap();
    let mut out = std::io::BufWriter::new(std::fs::File::create(outfile).unwrap());

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
        out.write_all(out_buffer.as_slice()).unwrap();
    }

    let total_elapsed = std::time::Instant::now() - start_time;
    log::info!("Total query length: {} nucleotides", total_query_length);
    log::info!("Elapsed time: {:.2} seconds ({:.2} ns / nucleotide)", total_elapsed.as_secs_f64(), total_elapsed.as_nanos() as f64 / total_query_length as f64);
    log::info!("Elapsed time excluding I/O: {:.2} seconds ({:.2} ns / nucleotide)", elapsed_without_io.as_secs_f64(), elapsed_without_io.as_nanos() as f64 / total_query_length as f64);

}

fn matching_statistics_command(matches: &clap::ArgMatches){
    let indexfile = matches.get_one::<std::path::PathBuf>("index").unwrap();
    let outfile = matches.get_one::<std::path::PathBuf>("output").unwrap();
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
            matching_statistics(&sbwt, &lcs, queryfile, outfile)
        }
    };

}

fn benchmark<SS: SubsetSeq>(sbwt: SbwtIndex<SS>, lcs: Option<LcsArray>) {
    log::info!("benchmarking index with k = {}, precalc length = {}, # kmers = {}, # sbwt sets = {}", sbwt.k(), sbwt.get_lookup_table().prefix_length, sbwt.n_kmers(), sbwt.n_sets());
    benchmark::benchmark_all(sbwt, lcs);
}

fn benchmark_command(matches: &clap::ArgMatches) {
    let indexfile = matches.get_one::<std::path::PathBuf>("index").unwrap();
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
            benchmark(sbwt, lcs);
        }
    };

}

fn dump_unitigs<SS: SubsetSeq + Send + Sync>(sbwt: &mut SbwtIndex<SS>, lcs: &Option<LcsArray>, n_threads: usize) {
    let out = std::io::stdout();
    let mut out = std::io::BufWriter::new(out);

    log::info!("Preparing the de Bruijn graph");
    sbwt.build_select();

    let dbg = Dbg::new(sbwt, lcs.as_ref(), n_threads);

    log::info!("Dumping unitigs");
    dbg.parallel_export_unitigs(&mut out);
}

fn dump_unitigs_command(matches: &clap::ArgMatches) {

    let n_threads = *matches.get_one::<usize>("threads").unwrap();
    let indexfile = matches.get_one::<std::path::PathBuf>("index").unwrap();
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
            log::info!("No LCS array found at {} -> Some extra computation needed to initialize the DBG", lcsfile.display());
            None
        }
    };

    match index {
        SbwtIndexVariant::SubsetMatrix(mut sbwt) => {
            dump_unitigs(&mut sbwt, &lcs, n_threads);
        }
    };

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
                .help("Output text file")
                .short('o')
                .long("output")
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

    match matches.subcommand(){
        Some(("build", sub_matches)) => build_command(sub_matches),
        Some(("build-lcs", sub_matches)) => build_lcs_command(sub_matches),
        Some(("lookup", sub_matches)) => lookup_query_command(sub_matches),
        Some(("matching-statistics", sub_matches)) => matching_statistics_command(sub_matches),
        Some(("benchmark", sub_matches)) => benchmark_command(sub_matches),
        Some(("dump-kmers", sub_matches)) => dump_kmers_command(sub_matches),
        Some(("dump-unitigs", sub_matches)) => dump_unitigs_command(sub_matches),
        _ => unreachable!(),
    }
    
}

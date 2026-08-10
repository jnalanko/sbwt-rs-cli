//! `set-operations` subcommand: builds one SBWT per `--input`, then runs merge/intersect/
//! difference between them and checks the results.
//!
//! Usage:
//!   sbwt-construction-test-harness set-operations --input <FILE> --input <FILE> [--input <FILE> ...] --out-dir <DIR>
//!
//! Unlike `build`'s `--input-list` (many sequence files merged into one SBWT), each `--input`
//! occurrence here is built into its own separate SBWT, and the set operations run between
//! those SBWTs.
//!
//! Not yet implemented; `run` currently panics with a `todo!()`.

use std::path::PathBuf;

use crate::common;

/// Fully parsed and validated configuration for a `set-operations` run.
struct Cli {
    inputs: Vec<PathBuf>,
    k: usize,
    sbwt_bin: PathBuf,
    out_dir: PathBuf,
    threads: usize,
    mem_gb: usize,
    keep: bool,
}

/// Defines the `set-operations` subcommand's CLI schema.
pub fn subcommand() -> clap::Command {
    clap::Command::new("set-operations")
        .about("Builds one SBWT per --input, then runs merge/intersect/difference between \
                them and checks the results. (Not yet implemented.)")
        .arg(clap::Arg::new("input")
            .help("Input fasta or fastq sequence file. Give at least twice, once per SBWT \
                   to build; unlike `build`'s --input-list, each occurrence here produces \
                   its own separate SBWT rather than being merged into one.")
            .short('i')
            .long("input")
            .required(true)
            .action(clap::ArgAction::Append)
            .value_parser(clap::value_parser!(PathBuf)))
        .arg(clap::Arg::new("k")
            .help("k-mer length, used for every build.")
            .short('k')
            .default_value("31")
            .value_parser(clap::value_parser!(usize)))
        .arg(clap::Arg::new("sbwt-bin")
            .help("Path to the sbwt executable to test. Defaults to the `sbwt` binary \
                   built alongside this harness.")
            .long("sbwt-bin")
            .value_parser(clap::value_parser!(PathBuf)))
        .arg(clap::Arg::new("out-dir")
            .help("Directory to write built SBWTs and set-operation outputs to. Created if \
                   it doesn't exist.")
            .long("out-dir")
            .required(true)
            .value_parser(clap::value_parser!(PathBuf)))
        .arg(clap::Arg::new("threads")
            .help("Number of threads to pass to the sbwt CLI, both for building and for \
                   the set operations.")
            .long("threads")
            .short('t')
            .default_value("4")
            .value_parser(clap::value_parser!(usize)))
        .arg(clap::Arg::new("mem-gb")
            .help("Memory budget in gigabytes to pass to the sbwt CLI's --mem-gb during \
                   the builds.")
            .long("mem-gb")
            .short('m')
            .default_value("8")
            .value_parser(clap::value_parser!(usize)))
        .arg(clap::Arg::new("keep")
            .help("Do not delete the output directory when done.")
            .long("keep")
            .action(clap::ArgAction::SetTrue))
}

/// Validates parsed argv: at least two `--input` files, all of which must exist, plus a
/// resolvable `sbwt` binary and a creatable `--out-dir`.
fn parse_args(matches: &clap::ArgMatches) -> Cli {
    let inputs: Vec<PathBuf> = matches.get_many::<PathBuf>("input").unwrap().cloned().collect();
    if inputs.len() < 2 {
        eprintln!(
            "error: --input must be given at least twice \
             (need at least 2 SBWTs to run set operations between)"
        );
        std::process::exit(2);
    }
    for path in &inputs {
        if !path.is_file() {
            eprintln!("error: input file {} does not exist", path.display());
            std::process::exit(2);
        }
    }

    let k = *matches.get_one::<usize>("k").unwrap();

    let sbwt_bin = matches.get_one::<PathBuf>("sbwt-bin").cloned()
        .unwrap_or_else(common::default_sbwt_bin_path);
    if !sbwt_bin.is_file() {
        eprintln!(
            "error: sbwt executable not found at {}\n\
             Build it first (e.g. `cargo build --release --features libsais`) or pass --sbwt-bin.",
            sbwt_bin.display()
        );
        std::process::exit(2);
    }

    let out_dir = matches.get_one::<PathBuf>("out-dir").unwrap().clone();
    std::fs::create_dir_all(&out_dir).expect("failed to create --out-dir");

    let threads = *matches.get_one::<usize>("threads").unwrap();
    let mem_gb = *matches.get_one::<usize>("mem-gb").unwrap();
    let keep = matches.get_flag("keep");

    Cli { inputs, k, sbwt_bin, out_dir, threads, mem_gb, keep }
}

pub fn run(matches: &clap::ArgMatches) {
    common::require_time_binary();
    let cli = parse_args(matches);

    println!("sbwt executable: {}", cli.sbwt_bin.display());
    println!("inputs:          {}", cli.inputs.len());
    for path in &cli.inputs {
        println!("  - {}", path.display());
    }
    println!("output dir:      {}", cli.out_dir.display());
    println!("k = {}, threads = {}, mem-gb = {}, keep = {}", cli.k, cli.threads, cli.mem_gb, cli.keep);

    todo!(
        "build one SBWT per --input, then run merge/intersect/difference/jaccard between \
         them and check the results"
    )
}

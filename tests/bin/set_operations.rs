//! `set-operations` subcommand: builds one SBWT per `--input`, then runs merge (union),
//! intersect and difference between every pair of them.
//!
//! Usage:
//!   sbwt-cli-test-harness set-operations --input <FILE> --input <FILE> [--input <FILE> ...] \
//!       --out-dir <DIR> --build-args "-k 31 --add-revcomp"
//!
//! Unlike `build`'s `--input-list` (many sequence files merged into one SBWT), each `--input`
//! occurrence here is built into its own separate SBWT, and the set operations run between
//! those SBWTs.
//!
//! `sbwt build` has no default `-k`, so a k-mer length (and any other build flags) must be
//! passed via `--build-args`, a single string split on whitespace and appended as-is to every
//! `sbwt build` invocation. There's no quoting support here (unlike a real shell) - if a path
//! in --build-args needs a space, this won't handle it.
//!
//! `difference` is not symmetric, so both directions (i\j and j\i) are run for every pair;
//! `merge` and `intersect` are symmetric, so only one call per pair is needed.
//!
//! Correctness of the results is not checked yet - this just exercises the CLI and reports
//! whether every build and every operation succeeded.

use std::path::{Path, PathBuf};
use std::process::Command;

use crate::common;
use crate::common::TimedRun;

/// Fully parsed and validated configuration for a `set-operations` run.
struct Cli {
    inputs: Vec<PathBuf>,
    build_args: Vec<String>,
    sbwt_bin: PathBuf,
    out_dir: PathBuf,
    threads: usize,
    keep: bool,
}

/// Defines the `set-operations` subcommand's CLI schema.
pub fn subcommand() -> clap::Command {
    clap::Command::new("set-operations")
        .about("Builds one SBWT per --input, then runs merge/intersect/difference between \
                every pair of them. Does not yet check the results for correctness.")
        .arg(clap::Arg::new("input")
            .help("Input fasta or fastq sequence file. Give at least twice, once per SBWT \
                   to build; unlike `build`'s --input-list, each occurrence here produces \
                   its own separate SBWT rather than being merged into one.")
            .short('i')
            .long("input")
            .required(true)
            .action(clap::ArgAction::Append)
            .value_parser(clap::value_parser!(PathBuf)))
        .arg(clap::Arg::new("build-args")
            .help("Extra arguments passed to every `sbwt build` invocation, split on \
                   whitespace (no quoting support), e.g. --build-args \"-k 31 --add-revcomp\". \
                   `sbwt build` has no default -k, so it must be supplied here.")
            .long("build-args")
            .default_value("")
            .value_parser(clap::value_parser!(String)))
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
        .arg(clap::Arg::new("keep")
            .help("Do not delete the output directory when done. Currently a no-op: nothing \
                   deletes the output directory yet, since there's no correctness check to \
                   consume its contents first.")
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

    let build_args: Vec<String> = matches.get_one::<String>("build-args").unwrap()
        .split_whitespace()
        .map(String::from)
        .collect();

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
    let keep = matches.get_flag("keep");

    Cli { inputs, build_args, sbwt_bin, out_dir, threads, keep }
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
    println!("build args:      {:?}", cli.build_args);
    println!("threads = {}, keep = {}\n", cli.threads, cli.keep);

    let sbwts: Vec<PathBuf> = cli.inputs.iter().enumerate()
        .map(|(i, input)| build_one(&cli, i, input))
        .collect();

    println!();
    let mut all_ok = true;
    for i in 0..sbwts.len() {
        for j in (i + 1)..sbwts.len() {
            all_ok &= run_op(&cli, "merge", &sbwts[i], &sbwts[j],
                              &cli.out_dir.join(format!("merge_{i}_{j}.sbwt")));
            all_ok &= run_op(&cli, "intersect", &sbwts[i], &sbwts[j],
                              &cli.out_dir.join(format!("intersect_{i}_{j}.sbwt")));
            all_ok &= run_op(&cli, "difference", &sbwts[i], &sbwts[j],
                              &cli.out_dir.join(format!("difference_{i}_{j}.sbwt")));
            all_ok &= run_op(&cli, "difference", &sbwts[j], &sbwts[i],
                              &cli.out_dir.join(format!("difference_{j}_{i}.sbwt")));
        }
    }

    if all_ok {
        println!("\nAll builds and set operations completed successfully.");
    } else {
        eprintln!("\nOne or more builds or set operations failed.");
        std::process::exit(1);
    }
}

/// Builds the SBWT for `input` (the `index`-th `--input`), wrapped in `/usr/bin/time -v`.
/// Exits the process if the build fails, since every later step depends on it.
fn build_one(cli: &Cli, index: usize, input: &Path) -> PathBuf {
    print!("Building SBWT {index} from {} ... ", input.display());
    let prefix = cli.out_dir.join(format!("sbwt-{index}"));
    let mut cmd = Command::new(&cli.sbwt_bin);
    cmd.arg("--threads").arg(cli.threads.to_string())
        .arg("build")
        .arg("--input").arg(input)
        .arg("--output-prefix").arg(&prefix)
        .arg("--temp-dir").arg(&cli.out_dir)
        .args(&cli.build_args);

    let run = common::run_timed(&cmd);
    report(&run);
    if !run.success {
        eprintln!("--- stderr ---\n{}", run.stderr);
        eprintln!("\nerror: failed to build SBWT {index} from {}", input.display());
        std::process::exit(1);
    }

    common::with_ext(&prefix, "sbwt")
}

/// Runs `sbwt <op> <sbwt1> <sbwt2> --output <output>`, wrapped in `/usr/bin/time -v`.
/// Returns whether it succeeded; does not exit the process on failure, so the remaining
/// operations still get a chance to run.
fn run_op(cli: &Cli, op: &str, sbwt1: &Path, sbwt2: &Path, output: &Path) -> bool {
    print!("{op:<10} {} , {} ... ", sbwt1.display(), sbwt2.display());
    let mut cmd = Command::new(&cli.sbwt_bin);
    cmd.arg("--threads").arg(cli.threads.to_string())
        .arg(op)
        .arg(sbwt1).arg(sbwt2)
        .arg("--output").arg(output);

    let run = common::run_timed(&cmd);
    report(&run);
    if !run.success {
        eprintln!("--- stderr ---\n{}", run.stderr);
    }
    run.success
}

fn report(run: &TimedRun) {
    if run.success {
        println!("ok ({})", common::format_duration(run.elapsed));
    } else {
        println!("FAILED (exit {:?}, {})", run.status_code, common::format_duration(run.elapsed));
    }
}

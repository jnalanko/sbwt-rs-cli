use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::Duration;

use crate::common;
use crate::common::TimedRun;

/// One `--input` (a sequence file to build) or `--sbwt-input` (a prebuilt SBWT to use as-is).
enum Source {
    Seq(PathBuf),
    Sbwt(PathBuf),
}

impl Source {
    fn path(&self) -> &Path {
        match self {
            Source::Seq(p) | Source::Sbwt(p) => p,
        }
    }
}

/// Fully parsed and validated configuration for a `set-operations` run.
struct Cli {
    inputs: Vec<Source>,
    build_args: Vec<String>,
    sbwt_bin: PathBuf,
    out_dir: PathBuf,
    threads: usize,
    keep: bool,
    report_path: PathBuf,
}

/// One row of the TSV report: a build, a pairwise set operation, or a check of one, its
/// inputs (by --input index), whether it succeeded, and what it cost.
struct ReportRow {
    op: String,
    left: String,
    right: String,
    success: bool,
    elapsed: Duration,
    peak_rss_kb: Option<u64>,
}

/// Defines the `set-operations` subcommand's CLI schema.
pub fn subcommand() -> clap::Command {
    clap::Command::new("set-operations")
        .about("Builds one SBWT per --input, then runs merge/intersect/difference between \
                every pair of them, verifying each result with `sbwt check-set-operation`.")
        .arg(clap::Arg::new("input")
            .help("Input fasta/fastq sequence file to build into an SBWT. Give at least \
                   twice, once per SBWT; unlike `build`'s --input-list, each occurrence here \
                   produces its own separate SBWT rather than being merged into one. \
                   Mutually exclusive with --sbwt-input.")
            .short('i')
            .long("input")
            .required(false)
            .action(clap::ArgAction::Append)
            .conflicts_with("sbwt-input")
            .value_parser(clap::value_parser!(PathBuf)))
        .arg(clap::Arg::new("sbwt-input")
            .help("A prebuilt SBWT file to use as-is, skipping the build step. Give at least \
                   twice. Mutually exclusive with --input.")
            .long("sbwt-input")
            .required(false)
            .action(clap::ArgAction::Append)
            .conflicts_with("input")
            .value_parser(clap::value_parser!(PathBuf)))
        .arg(clap::Arg::new("build-args")
            .help("Extra arguments passed to every `sbwt build` invocation, split on \
                   whitespace (no quoting support), e.g. --build-args \"-k 31 --add-revcomp\". \
                   `sbwt build` has no default -k, so it must be supplied here.")
            .long("build-args")
            .default_value("")
            // The value is expected to start with "-" (it's a string of sbwt-build flags),
            // but clap's default heuristic refuses to consume a value that looks like a flag.
            // Without this, e.g. `--build-args "-k 31"` fails with "unexpected argument '-k'".
            .allow_hyphen_values(true)
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
        .arg(clap::Arg::new("report")
            .help("Path to write the TSV report (one row per build and per set operation: \
                   which step, the --input index/indices involved, success, time, peak \
                   memory) to. Defaults to set-operations-report.tsv inside --out-dir; a path \
                   given here is used as-is, so it doesn't have to be under --out-dir.")
            .long("report")
            .value_parser(clap::value_parser!(PathBuf)))
        .arg(clap::Arg::new("keep")
            .help("Do not delete the output directory when done. Currently a no-op: nothing \
                   deletes the output directory yet, since there's no correctness check to \
                   consume its contents first.")
            .long("keep")
            .action(clap::ArgAction::SetTrue))
}

/// Validates parsed argv: at least two `--input` XOR `--sbwt-input` sources (enforced mutually
/// exclusive by clap), all of which must exist, plus a resolvable `sbwt` binary and a
/// creatable `--out-dir`.
fn parse_args(matches: &clap::ArgMatches) -> Cli {
    let seq_inputs = matches.get_many::<PathBuf>("input").unwrap_or_default().cloned();
    let sbwt_inputs = matches.get_many::<PathBuf>("sbwt-input").unwrap_or_default().cloned();
    let inputs: Vec<Source> = seq_inputs.map(Source::Seq)
        .chain(sbwt_inputs.map(Source::Sbwt))
        .collect();
    if inputs.len() < 2 {
        eprintln!(
            "error: --input or --sbwt-input must be given at least twice \
             (need at least 2 SBWTs to run set operations between)"
        );
        std::process::exit(2);
    }
    for source in &inputs {
        if !source.path().is_file() {
            eprintln!("error: input file {} does not exist", source.path().display());
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

    let report_path = matches.get_one::<PathBuf>("report").cloned()
        .unwrap_or_else(|| out_dir.join("set-operations-report.tsv"));

    Cli { inputs, build_args, sbwt_bin, out_dir, threads, keep, report_path }
}

/// Creates (truncating any existing file) the report and writes just its header. Rows are
/// appended one at a time afterwards, by `append_report_row`, as each build or operation
/// finishes.
fn init_report(path: &Path) -> std::io::Result<()> {
    common::init_tsv_report(path, "op\tleft\tright\tsuccess\ttime_seconds\tpeak_rss_bytes")
}

/// Appends one TSV row to the report: a build's or operation's inputs (by --input index),
/// whether it succeeded, wall-clock time in seconds, and peak resident memory in bytes.
fn append_report_row(path: &Path, row: &ReportRow) -> std::io::Result<()> {
    let peak_rss_bytes = row.peak_rss_kb.map(|kb| kb * 1024);
    let line = format!(
        "{}\t{}\t{}\t{}\t{:.3}\t{}",
        row.op,
        row.left,
        row.right,
        row.success,
        row.elapsed.as_secs_f64(),
        peak_rss_bytes.map(|b| b.to_string()).unwrap_or_default(),
    );
    common::append_tsv_row(path, &line)
}

pub fn run(matches: &clap::ArgMatches) {
    common::require_time_binary();
    let cli = parse_args(matches);

    println!("sbwt executable: {}", cli.sbwt_bin.display());
    println!("inputs:          {}", cli.inputs.len());
    for source in &cli.inputs {
        let tag = match source {
            Source::Seq(_) => "seq",
            Source::Sbwt(_) => "sbwt",
        };
        println!("  - [{tag}] {}", source.path().display());
    }
    println!("output dir:      {}", cli.out_dir.display());
    println!("build args:      {:?}", cli.build_args);
    println!("threads = {}, keep = {}", cli.threads, cli.keep);

    // Written to as each build or operation finishes (see `append_report_row`), not batched
    // up until the end, so a partial report still exists if the harness is interrupted or a
    // later step hangs.
    init_report(&cli.report_path).expect("failed to create report");
    println!("report:          {} (appended to as each step finishes)\n", cli.report_path.display());

    let sbwts: Vec<PathBuf> = cli.inputs.iter().enumerate()
        .map(|(i, source)| build_one(&cli, i, source))
        .collect();

    println!();
    let mut all_ok = true;
    for i in 0..sbwts.len() {
        for j in (i + 1)..sbwts.len() {
            all_ok &= run_pair(&cli, "merge", i, &sbwts[i], j, &sbwts[j],
                                &cli.out_dir.join(format!("merge_{i}_{j}.sbwt")));
            all_ok &= run_pair(&cli, "intersect", i, &sbwts[i], j, &sbwts[j],
                                &cli.out_dir.join(format!("intersect_{i}_{j}.sbwt")));
            all_ok &= run_pair(&cli, "difference", i, &sbwts[i], j, &sbwts[j],
                                &cli.out_dir.join(format!("difference_{i}_{j}.sbwt")));
            all_ok &= run_pair(&cli, "difference", j, &sbwts[j], i, &sbwts[i],
                                &cli.out_dir.join(format!("difference_{j}_{i}.sbwt")));
        }
    }

    if all_ok {
        println!("\nAll builds, set operations, and checks completed successfully.");
    } else {
        eprintln!("\nOne or more builds, set operations, or checks failed.");
        std::process::exit(1);
    }
}

/// Runs one pairwise set operation and, if it succeeds, verifies its output with
/// `sbwt check-set-operation`. Returns whether both the operation and its check succeeded.
fn run_pair(cli: &Cli, op: &'static str, i: usize, sbwt1: &Path, j: usize, sbwt2: &Path, output: &Path) -> bool {
    if !run_op(cli, op, i, sbwt1, j, sbwt2, output) {
        return false;
    }
    check_op(cli, op, i, sbwt1, j, sbwt2, output)
}

/// Returns the SBWT for the `index`-th source: used directly if it's a `--sbwt-input`
/// (prebuilt), otherwise built from a `--input` sequence file, wrapped in `/usr/bin/time -v`.
/// Exits the process if the build fails, since every later step depends on it.
fn build_one(cli: &Cli, index: usize, source: &Source) -> PathBuf {
    let input = match source {
        Source::Sbwt(path) => {
            println!("SBWT {index}: using prebuilt {}", path.display());
            return path.clone();
        }
        Source::Seq(path) => path,
    };

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
    print_status(&run);

    let row = ReportRow {
        op: "build".to_string(), left: index.to_string(), right: String::new(),
        success: run.success, elapsed: run.elapsed, peak_rss_kb: run.peak_rss_kb,
    };
    append_report_row(&cli.report_path, &row).expect("failed to append to report");

    if !run.success {
        eprintln!("--- stderr ---\n{}", run.stderr);
        eprintln!("\nerror: failed to build SBWT {index} from {}", input.display());
        std::process::exit(1);
    }

    common::with_ext(&prefix, "sbwt")
}

/// Runs `sbwt <op> <sbwt1> <sbwt2> --output <output>`, wrapped in `/usr/bin/time -v`, where
/// `sbwt1`/`sbwt2` are the `i`-th/`j`-th built SBWTs (recorded in the report by index, not by
/// path). Returns whether it succeeded; does not exit the process on failure, so the
/// remaining operations still get a chance to run.
fn run_op(cli: &Cli, op: &'static str, i: usize, sbwt1: &Path, j: usize, sbwt2: &Path, output: &Path) -> bool {
    print!("{op:<10} {i} , {j} ... ");
    let mut cmd = Command::new(&cli.sbwt_bin);
    cmd.arg("--threads").arg(cli.threads.to_string())
        .arg(op)
        .arg(sbwt1).arg(sbwt2)
        .arg("--output").arg(output);

    let run = common::run_timed(&cmd);
    print_status(&run);

    let row = ReportRow {
        op: op.to_string(), left: i.to_string(), right: j.to_string(),
        success: run.success, elapsed: run.elapsed, peak_rss_kb: run.peak_rss_kb,
    };
    append_report_row(&cli.report_path, &row).expect("failed to append to report");

    if !run.success {
        eprintln!("--- stderr ---\n{}", run.stderr);
    }
    run.success
}

/// Runs `sbwt check-set-operation` to verify that `result` (the output of `op` on the i-th
/// and j-th built SBWTs) has exactly the expected k-mer set, computed from the original
/// `--input` sequence files where available. For a `--sbwt-input` source there's no original
/// sequence file, so `--seq1`/`--seq2` is omitted and `check-set-operation` falls back to
/// exporting that SBWT's own unitigs. Returns whether the check passed.
fn check_op(cli: &Cli, op: &'static str, i: usize, sbwt1: &Path, j: usize, sbwt2: &Path, result: &Path) -> bool {
    print!("{op:<10} {i} , {j} check ... ");
    let mut cmd = Command::new(&cli.sbwt_bin);
    cmd.arg("check-set-operation").arg("--op").arg(op);
    if let Source::Seq(path) = &cli.inputs[i] {
        cmd.arg("--seq1").arg(path);
    }
    if let Source::Seq(path) = &cli.inputs[j] {
        cmd.arg("--seq2").arg(path);
    }
    cmd.arg("--sbwt1").arg(sbwt1)
        .arg("--sbwt2").arg(sbwt2)
        .arg("--result").arg(result)
        .arg("--temp-dir").arg(&cli.out_dir);

    let run = common::run_timed(&cmd);
    print_status(&run);

    let row = ReportRow {
        op: format!("{op}-check"), left: i.to_string(), right: j.to_string(),
        success: run.success, elapsed: run.elapsed, peak_rss_kb: run.peak_rss_kb,
    };
    append_report_row(&cli.report_path, &row).expect("failed to append to report");

    if !run.success {
        eprintln!("--- stderr ---\n{}", run.stderr);
    }
    run.success
}

fn print_status(run: &TimedRun) {
    if run.success {
        println!("ok ({})", common::format_duration(run.elapsed));
    } else {
        println!("FAILED (exit {:?}, {})", run.status_code, common::format_duration(run.elapsed));
    }
}

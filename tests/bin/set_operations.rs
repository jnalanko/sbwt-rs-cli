use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::Duration;

use crate::common;
use crate::common::TimedRun;

enum Inputs {
    Sbwts(Vec<PathBuf>),
    Seqs(Vec<PathBuf>),
}

/// Fully parsed and validated configuration for a `set-operations` run.
struct Cli {
    inputs: Inputs,
    build_args: Vec<String>,
    sbwt_bin: PathBuf,
    out_dir: PathBuf,
    threads: usize,
    keep: bool,
    verbose: bool,
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

/// Builds one SBWT per --input, then runs merge/intersect/difference between every pair of
/// them, verifying each result with `sbwt check-set-operation`.
#[derive(clap::Args)]
pub struct Args {

    #[arg(long, conflicts_with = "sbwt_input", num_args = 2.., help = "Give either this or sbwt-inputs.")]
    seqfiles: Option<Vec<PathBuf>>,

    #[arg(long, num_args = 2.., help = "Give either this or --seqfiles.")]
    sbwt_inputs: Option<Vec<PathBuf>>,

    /// Extra arguments passed to every `sbwt build` invocation, split on whitespace (no
    /// quoting support), e.g. --build-args "-k 31 --add-revcomp". `sbwt build` has no default
    /// -k, so it must be supplied here.
    // The value is expected to start with "-" (it's a string of sbwt-build flags), but clap's
    // default heuristic refuses to consume a value that looks like a flag. Without
    // allow_hyphen_values, e.g. `--build-args "-k 31"` fails with "unexpected argument '-k'".
    #[arg(long, default_value = "", allow_hyphen_values = true)]
    build_args: String,

    /// Path to the sbwt executable to test. Defaults to the `sbwt` binary
    /// built alongside this harness.
    #[arg(long)]
    sbwt_bin: Option<PathBuf>,

    /// Directory to write built SBWTs and set-operation outputs to. Created if it doesn't
    /// exist.
    #[arg(long, required = true)]
    out_dir: PathBuf,

    /// Number of threads to pass to the sbwt CLI, both for building and for the set
    /// operations.
    #[arg(short, long, default_value_t = 4)]
    threads: usize,

    /// Path to write the TSV report (one row per build and per set operation: which step, the
    /// --input index/indices involved, success, time, peak memory) to. Defaults to
    /// set-operations-report.tsv inside --out-dir; a path given here is used as-is, so it
    /// doesn't have to be under --out-dir.
    #[arg(long)]
    report: Option<PathBuf>,

    /// Pass --keep-unitigs to `sbwt check-set-operation`, so that the unitig fasta files it
    /// generates in --out-dir survive the run. Only a --sbwt-input source generates any: one
    /// given as --input is checked against its own sequence file instead. Nothing else the run
    /// writes is deleted either way, so this affects nothing when there is no --sbwt-input.
    #[arg(long)]
    keep: bool,

    /// Pass -v to the sbwt CLI, so that it logs at its debug level. Its output is shown on
    /// the harness's own stdout and stderr either way.
    #[arg(short, long)]
    verbose: bool,
}

/// Validates parsed argv: at least two `--input` XOR `--sbwt-input` sources (enforced mutually
/// exclusive by clap), all of which must exist, plus a resolvable `sbwt` binary and a
/// creatable `--out-dir`.
fn parse_args(args: Args) -> Cli {

    let inputs = if let Some(seqfiles) = args.seqfiles {
        assert!(seqfiles.len() >= 2);
        Inputs::Seqs(seqfiles)
    } else if let Some(sbwts) = args.sbwt_inputs {
        assert!(sbwts.len() >= 2);
        Inputs::Sbwts(sbwts)
    } else {
        panic!("Must give either --seqfiles or --sbwts");
    };

    let build_args: Vec<String> = args.build_args
        .split_whitespace()
        .map(String::from)
        .collect();

    let sbwt_bin = args.sbwt_bin.unwrap_or_else(common::default_sbwt_bin_path);
    if !sbwt_bin.is_file() {
        log::error!(
            "sbwt executable not found at {}\n\
             Build it first (e.g. `cargo build --release --features libsais`) or pass --sbwt-bin.",
            sbwt_bin.display()
        );
        std::process::exit(2);
    }

    let out_dir = args.out_dir;
    std::fs::create_dir_all(&out_dir).expect("failed to create --out-dir");

    let threads = args.threads;
    let keep = args.keep;
    let verbose = args.verbose;

    let report_path = args.report
        .unwrap_or_else(|| out_dir.join("set-operations-report.tsv"));

    Cli { inputs, build_args, sbwt_bin, out_dir, threads, keep, verbose, report_path }
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

pub fn run(args: Args) {
    common::require_time_binary();
    let cli = parse_args(args);

    log::info!("sbwt executable: {}", cli.sbwt_bin.display());
    log::info!("output dir:      {}", cli.out_dir.display());
    log::info!("build args:      {:?}", cli.build_args);
    log::info!("threads = {}, keep = {}, verbose = {}", cli.threads, cli.keep, cli.verbose);

    // Written to as each build or operation finishes (see `append_report_row`), not batched
    // up until the end, so a partial report still exists if the harness is interrupted or a
    // later step hangs.
    init_report(&cli.report_path).expect("failed to create report");
    log::info!("report:          {} (appended to as each step finishes)", cli.report_path.display());

    let (seqfiles, sbwts) = match &cli.inputs {
        Inputs::Sbwts(sbwts) => {
            // SBWTs given, build seqfiles
            let seqfiles: Vec<PathBuf> = sbwts.iter().enumerate().map(|(i, sbwt_file)| {
                build_unitigs(&cli, i, sbwt_file)
            }).collect();
            (seqfiles, sbwts.clone())
        },
        Inputs::Seqs(seqfiles) => {
            // Seqfiles given, build SBWTs
            let sbwts: Vec<PathBuf> = seqfiles.iter().enumerate().map(|(i, seqfile)| {
                build_one(&cli, i, seqfile)
            }).collect();
            (seqfiles.clone(), sbwts)
        },
    };

    let mut all_ok = true;
    for i in 0..sbwts.len() {
        for j in (i + 1)..sbwts.len() {
            all_ok &= run_pair(&cli, "merge", i, &seqfiles[i], &sbwts[i], j, &seqfiles[j], &sbwts[j],
                                &cli.out_dir.join(format!("merge_{i}_{j}.sbwt")));
            all_ok &= run_pair(&cli, "intersect", i, &seqfiles[i], &sbwts[i], j, &seqfiles[j], &sbwts[j],
                                &cli.out_dir.join(format!("intersect_{i}_{j}.sbwt")));
            all_ok &= run_pair(&cli, "difference", i, &seqfiles[i], &sbwts[i], j, &seqfiles[j], &sbwts[j],
                                &cli.out_dir.join(format!("difference_{i}_{j}.sbwt")));
            all_ok &= run_pair(&cli, "difference", j, &seqfiles[i], &sbwts[j], i, &seqfiles[j], &sbwts[i],
                                &cli.out_dir.join(format!("difference_{j}_{i}.sbwt")));
        }
    }

    if all_ok {
        log::info!("All builds, set operations, and checks completed successfully.");
    } else {
        log::error!("One or more builds, set operations, or checks failed.");
        std::process::exit(1);
    }
}

/// Runs one pairwise set operation and, if it succeeds, verifies its output with
/// `sbwt check-set-operation`. Returns whether both the operation and its check succeeded.
fn run_pair(cli: &Cli, op: &'static str, i: usize, seqs1: &Path, sbwt1: &Path, j: usize, seqs2: &Path, sbwt2: &Path, output: &Path) -> bool {
    if !run_op(cli, op, i, sbwt1, j, sbwt2, output) {
        return false;
    }
    check_op(cli, op, i, seqs1, sbwt1, j, seqs2, sbwt2, output)
}

/// Returns a fasta file of the unitigs of the `index`-th `--sbwt-input`, dumped with
/// `sbwt dump-unitigs` and wrapped in `/usr/bin/time -v`. This is the stand-in for the original
/// sequence file, which a prebuilt SBWT doesn't come with: it spells the same k-mer set, which
/// is all the later steps need. Exits the process if the dump fails, since every later step
/// depends on it.
fn build_unitigs(cli: &Cli, index: usize, sbwtfile: &PathBuf) -> PathBuf {
    let label = format!("Dumping unitigs {index} from {}", sbwtfile.display());
    log::info!("{label} ...");
    let output = cli.out_dir.join(format!("unitigs-{index}.fna"));
    let mut cmd = Command::new(&cli.sbwt_bin);
    verbose_flag(cli, &mut cmd);
    cmd.arg("--threads").arg(cli.threads.to_string())
        .arg("dump-unitigs")
        .arg("--index").arg(sbwtfile)
        .arg("--output").arg(&output);

    let run = common::run_timed(&cmd, &cli.out_dir);
    log_status(&label, &run);

    let row = ReportRow {
        op: "dump-unitigs".to_string(), left: index.to_string(), right: String::new(),
        success: run.success, elapsed: run.elapsed, peak_rss_kb: run.peak_rss_kb,
    };
    append_report_row(&cli.report_path, &row).expect("failed to append to report");

    if !run.success {
        log::error!("failed to dump unitigs {index} from {}", sbwtfile.display());
        std::process::exit(1);
    }

    output
}

/// Returns the SBWT for the `index`-th source: used directly if it's a `--sbwt-input`
/// (prebuilt), otherwise built from a `--input` sequence file, wrapped in `/usr/bin/time -v`.
/// Exits the process if the build fails, since every later step depends on it.
fn build_one(cli: &Cli, index: usize, seqfile: &PathBuf) -> PathBuf {
    let label = format!("Building SBWT {index} from {}", seqfile.display());
    log::info!("{label} ...");
    let prefix = cli.out_dir.join(format!("sbwt-{index}"));
    let mut cmd = Command::new(&cli.sbwt_bin);
    verbose_flag(cli, &mut cmd);
    cmd.arg("--threads").arg(cli.threads.to_string())
        .arg("build")
        .arg("--input").arg(seqfile)
        .arg("--output-prefix").arg(&prefix)
        .arg("--temp-dir").arg(&cli.out_dir)
        .args(&cli.build_args);

    let run = common::run_timed(&cmd, &cli.out_dir);
    log_status(&label, &run);

    let row = ReportRow {
        op: "build".to_string(), left: index.to_string(), right: String::new(),
        success: run.success, elapsed: run.elapsed, peak_rss_kb: run.peak_rss_kb,
    };
    append_report_row(&cli.report_path, &row).expect("failed to append to report");

    if !run.success {
        log::error!("failed to build SBWT {index} from {}", seqfile.display());
        std::process::exit(1);
    }

    common::with_ext(&prefix, "sbwt")
}

/// Runs `sbwt <op> <sbwt1> <sbwt2> --output <output>`, wrapped in `/usr/bin/time -v`, where
/// `sbwt1`/`sbwt2` are the `i`-th/`j`-th built SBWTs (recorded in the report by index, not by
/// path). Returns whether it succeeded; does not exit the process on failure, so the
/// remaining operations still get a chance to run.
fn run_op(cli: &Cli, op: &'static str, i: usize, sbwt1: &Path, j: usize, sbwt2: &Path, output: &Path) -> bool {
    let label = format!("{op:<10} {i} , {j}");
    log::info!("{label} ...");
    let mut cmd = Command::new(&cli.sbwt_bin);
    verbose_flag(cli, &mut cmd);
    cmd.arg("--threads").arg(cli.threads.to_string())
        .arg(op)
        .arg(sbwt1).arg(sbwt2)
        .arg("--output").arg(output);

    let run = common::run_timed(&cmd, &cli.out_dir);
    log_status(&label, &run);

    let row = ReportRow {
        op: op.to_string(), left: i.to_string(), right: j.to_string(),
        success: run.success, elapsed: run.elapsed, peak_rss_kb: run.peak_rss_kb,
    };
    append_report_row(&cli.report_path, &row).expect("failed to append to report");

    run.success
}

/// Runs `sbwt check-set-operation` to verify that `result` (the output of `op` on the i-th
/// and j-th built SBWTs) has exactly the expected k-mer set, computed from the original
/// `--input` sequence files where available. For a `--sbwt-input` source there's no original
/// sequence file, so `--seq1`/`--seq2` is omitted and `check-set-operation` falls back to
/// exporting that SBWT's own unitigs. Returns whether the check passed.
fn check_op(cli: &Cli, op: &'static str, i: usize, seqs1: &Path, sbwt1: &Path, j: usize, seqs2: &Path, sbwt2: &Path, result: &Path) -> bool {
    let label = format!("{op:<10} {i} , {j} check");
    log::info!("{label} ...");
    let mut cmd = Command::new(&cli.sbwt_bin);
    verbose_flag(cli, &mut cmd);
    cmd.arg("check-set-operation").arg("--op").arg(op)
        .arg("--seq1").arg(seqs1)
        .arg("--seq2").arg(seqs2)
        .arg("--sbwt1").arg(sbwt1)
        .arg("--sbwt2").arg(sbwt2)
        .arg("--result").arg(result)
        .arg("--temp-dir").arg(&cli.out_dir);
    if cli.keep {
        cmd.arg("--keep-unitigs");
    }

    let run = common::run_timed(&cmd, &cli.out_dir);
    log_status(&label, &run);

    let row = ReportRow {
        op: format!("{op}-check"), left: i.to_string(), right: j.to_string(),
        success: run.success, elapsed: run.elapsed, peak_rss_kb: run.peak_rss_kb,
    };
    append_report_row(&cli.report_path, &row).expect("failed to append to report");

    run.success
}

/// Adds the sbwt CLI's own -v to `cmd` under --verbose, making the output it streams to the
/// harness's stderr the more detailed one. It has to go before the subcommand name, so this is
/// called before any other argument is pushed.
fn verbose_flag(cli: &Cli, cmd: &mut Command) {
    if cli.verbose {
        cmd.arg("-v");
    }
}

/// Logs how the step described by `label` (already logged as "<label> ..." before the step
/// started) turned out. The label is repeated so the outcome line stands on its own, since
/// other steps' output may be interleaved between the two.
fn log_status(label: &str, run: &TimedRun) {
    if run.success {
        log::info!("{label} ... ok ({})", common::format_duration(run.elapsed));
    } else {
        log::error!(
            "{label} ... FAILED (exit {:?}, {})",
            run.status_code, common::format_duration(run.elapsed)
        );
    }
}

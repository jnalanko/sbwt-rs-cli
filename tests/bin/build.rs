use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::Duration;

use crate::common;

/// Values used for the parameters the user left out. Kept here rather than in clap's
/// `default_value`, because a defaulted parameter is exactly the one --fuzz is allowed to vary,
/// so the parsed value has to stay distinguishable from an explicitly given one.
const DEFAULT_K: usize = 31;
const DEFAULT_THREADS: usize = 4;
const DEFAULT_MEM_GB: usize = 8;

struct Algorithm {
    name: &'static str,
    extra_args: &'static [&'static str],
}

// These are the only construction strategies `sbwt build` currently exposes; despite using
// unrelated internal algorithms, they're all specified to produce byte-identical output for
// the same input. "on-disk" is listed first because it becomes the baseline everything else
// is diffed against in `compare_outputs` (it's the CLI's own default, so it's the one most
// likely to already be correct if something regresses).
const ALGORITHMS: &[Algorithm] = &[
    Algorithm { name: "on-disk", extra_args: &[] },
    Algorithm { name: "in-memory", extra_args: &["--in-memory"] },
    Algorithm { name: "bounded-suffix-sort", extra_args: &["--bounded-suffix-sort"] },
    Algorithm { name: "libsais", extra_args: &["--via-libsais"] },
];

/// The build parameters that --fuzz may vary across combinations.
#[derive(Clone, Copy)]
struct Params {
    k: usize,
    threads: usize,
    add_revcomp: bool,
    add_all_dummy_paths: bool,
    dedup_batches: bool,
}

impl Params {
    fn label(&self) -> String {
        format!(
            "k={} threads={} add-revcomp={} add-all-dummy-paths={} dedup-batches={}",
            self.k, self.threads, self.add_revcomp, self.add_all_dummy_paths, self.dedup_batches
        )
    }

    fn dir_name(&self, index: usize) -> String {
        format!(
            "combo-{index}_k{}_t{}_rc{}_adp{}_db{}",
            self.k,
            self.threads,
            self.add_revcomp as u8,
            self.add_all_dummy_paths as u8,
            self.dedup_batches as u8
        )
    }
}

struct AlgoRun {
    name: &'static str,
    prefix: PathBuf,
    success: bool,
    elapsed: Duration,
    peak_rss_kb: Option<u64>,
}

/// One row of the TSV report: the parameters and algorithm that produced one build, plus
/// what it cost.
struct ReportRow {
    params: Params,
    algorithm: &'static str,
    success: bool,
    elapsed: Duration,
    peak_rss_kb: Option<u64>,
}

/// Outcome of building and comparing all four algorithms for one parameter combination.
struct ComboOutcome {
    label: String,
    failed_algorithms: Vec<&'static str>,
    mismatches: Vec<String>,
}

impl ComboOutcome {
    fn is_ok(&self) -> bool {
        self.failed_algorithms.is_empty() && self.mismatches.is_empty()
    }
}

/// Fully parsed and validated configuration for a `build` run.
struct Cli {
    input_arg: &'static str,
    input_path: PathBuf,
    mem_gb: usize,
    verbose: bool,
    sbwt_bin: PathBuf,
    out_dir: PathBuf,
    keep: bool,
    combos: Vec<Params>,
    report_path: PathBuf,
}

/// Runs all sbwt construction algorithms through the sbwt CLI on the same input and checks
/// that they produce byte-identical output.
#[derive(clap::Args)]
pub struct Args {
    /// Input fasta or fastq sequence file (same as `sbwt build --input`).
    #[arg(short, long, required_unless_present = "input_list", conflicts_with = "input_list")]
    input: Option<PathBuf>,

    /// File listing input fasta/fastq files, one path per line
    /// (same format as `sbwt build --input-list`).
    #[arg(long, required_unless_present = "input")]
    input_list: Option<PathBuf>,

    /// k-mer length [default: 31]
    #[arg(short)]
    k: Option<usize>,

    /// Path to the sbwt executable to test. Defaults to the `sbwt` binary
    /// built alongside this harness.
    #[arg(long)]
    sbwt_bin: Option<PathBuf>,

    /// Directory to write build outputs (and on-disk algorithm scratch files) to.
    /// Created if it doesn't exist.
    #[arg(long, required = true)]
    out_dir: PathBuf,

    /// Do not delete the output directory when done.
    #[arg(long)]
    keep: bool,

    /// Path to write the TSV report (one row per build: parameters, time, peak memory) to.
    /// Defaults to construction-report.tsv inside --out-dir; a path given here is used as-is,
    /// so it doesn't have to be under --out-dir.
    #[arg(long)]
    report: Option<PathBuf>,

    /// Number of threads to pass to the sbwt CLI [default: 4]
    #[arg(short, long)]
    threads: Option<usize>,

    /// Memory budget in gigabytes to pass to the sbwt CLI's --mem-gb
    /// (used by the on-disk and in-memory algorithms only). Never fuzzed [default: 8]
    #[arg(short = 'm', long)]
    mem_gb: Option<usize>,

    /// Pass --add-revcomp to every build.
    #[arg(long)]
    add_revcomp: bool,

    /// Pass --add-all-dummy-paths to every build.
    #[arg(long)]
    add_all_dummy_paths: bool,

    /// Pass --dedup-batches to every build.
    #[arg(short, long)]
    dedup_batches: bool,

    /// For each of -k, --threads, --add-revcomp, --add-all-dummy-paths and --dedup-batches
    /// that was not explicitly given, try multiple values instead of just the default, and
    /// build+compare every combination. --mem-gb is never fuzzed.
    #[arg(long)]
    fuzz: bool,

    /// Pass -v to the sbwt CLI, and print its stderr for every build (not just failed ones).
    #[arg(short, long)]
    verbose: bool,
}

/// Validates parsed argv and works out the set of parameter combinations to run (a single
/// one, unless --fuzz expands some of them).
fn parse_args(args: Args) -> Cli {
    let (input_arg, input_path) = match (args.input, args.input_list) {
        (Some(path), None) => ("--input", path),
        (None, Some(path)) => ("--input-list", path),
        _ => unreachable!("clap enforces exactly one of --input / --input-list"),
    };
    if !input_path.is_file() {
        log::error!("input file {} does not exist", input_path.display());
        std::process::exit(2);
    }

    let mem_gb = args.mem_gb.unwrap_or(DEFAULT_MEM_GB);
    let verbose = args.verbose;
    let fuzz = args.fuzz;

    // Only parameters the user left at their default get multiple candidate values; anything
    // passed explicitly is pinned to that one value. This lets a caller narrow the search
    // space to whatever they're actually unsure about instead of always paying for the full
    // combinatorial product.
    let k_candidates = candidates_usize(fuzz, args.k, DEFAULT_K, &[3, 31, 32, 100]);
    let threads_candidates = candidates_usize(fuzz, args.threads, DEFAULT_THREADS, &[1, 4]);
    let add_revcomp_candidates = candidates_bool(fuzz, args.add_revcomp);
    let add_all_dummy_paths_candidates = candidates_bool(fuzz, args.add_all_dummy_paths);
    let dedup_batches_candidates = candidates_bool(fuzz, args.dedup_batches);

    let mut combos = Vec::new();
    for &k in &k_candidates {
        for &threads in &threads_candidates {
            for &add_revcomp in &add_revcomp_candidates {
                for &add_all_dummy_paths in &add_all_dummy_paths_candidates {
                    for &dedup_batches in &dedup_batches_candidates {
                        combos.push(Params { k, threads, add_revcomp, add_all_dummy_paths, dedup_batches });
                    }
                }
            }
        }
    }

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
    let keep = args.keep;

    let report_path = args.report
        .unwrap_or_else(|| out_dir.join("construction-report.tsv"));

    Cli { input_arg, input_path, mem_gb, verbose, sbwt_bin, out_dir, keep, combos, report_path }
}

/// Candidate values for a boolean build flag: both `false` and `true` if it's being fuzzed,
/// otherwise just the value that was parsed. A flag can only ever be given explicitly as
/// `true`, so a `true` here means the user pinned it and it isn't fuzzed.
fn candidates_bool(fuzz: bool, actual: bool) -> Vec<bool> {
    if fuzz && !actual { vec![false, true] } else { vec![actual] }
}

/// Same as [`candidates_bool`], but for a numeric option, with the fuzzed value set given
/// explicitly by the caller. `actual` is `None` when the option was left out, i.e. when it is
/// eligible for fuzzing.
fn candidates_usize(fuzz: bool, actual: Option<usize>, default: usize, fuzzed: &[usize]) -> Vec<usize> {
    match actual {
        Some(value) => vec![value],
        None if fuzz => fuzzed.to_vec(),
        None => vec![default],
    }
}

pub fn run(args: Args) {
    common::require_time_binary();
    let cli = parse_args(args);

    log::info!("sbwt executable: {}", cli.sbwt_bin.display());
    log::info!("input ({}): {}", cli.input_arg, cli.input_path.display());
    log::info!("output dir:      {}", cli.out_dir.display());
    log::info!("mem-gb = {}", cli.mem_gb);
    if cli.combos.len() > 1 {
        log::info!("fuzzing {} parameter combination(s)", cli.combos.len());
    } else {
        log::info!("{}", cli.combos[0].label());
    }

    // Written to as each build finishes (see `append_report_row`), not batched up until the
    // end, so a partial report still exists if the harness is interrupted or a later build
    // hangs.
    init_report(&cli.report_path).expect("failed to create report");
    log::info!("report:          {} (appended to as each build finishes)", cli.report_path.display());

    let mut outcomes = Vec::new();
    for (i, params) in cli.combos.iter().enumerate() {
        let combo_dir = if cli.combos.len() > 1 { cli.out_dir.join(params.dir_name(i)) } else { cli.out_dir.clone() };
        let label = if cli.combos.len() > 1 { Some(params.label()) } else { None };
        outcomes.push(run_combo(&cli, params, &combo_dir, label.as_deref()));
    }

    cleanup(&cli.out_dir, cli.keep, &cli.report_path);

    if cli.combos.len() > 1 {
        log::info!("=== fuzz summary ===");
        for outcome in &outcomes {
            let status = if outcome.is_ok() { "OK" } else { "FAILED" };
            log::info!("{status:<6} {}", outcome.label);
        }
    }

    if outcomes.iter().all(ComboOutcome::is_ok) {
        log::info!("All combination(s) agree across all construction algorithms.");
    } else {
        let failed = outcomes.iter().filter(|o| !o.is_ok()).count();
        log::error!("{failed} of {} combination(s) failed.", outcomes.len());
        std::process::exit(1);
    }
}

/// Builds all four construction algorithms with the given parameters and byte-compares
/// their `.sbwt` and `.lcs` output. `label`, if given, is logged as a header (used when
/// running more than one parameter combination under --fuzz).
fn run_combo(cli: &Cli, params: &Params, combo_dir: &Path, label: Option<&str>) -> ComboOutcome {
    if let Some(label) = label {
        log::info!("--- {label} ---");
    }
    std::fs::create_dir_all(combo_dir).expect("failed to create combo output dir");

    let runs: Vec<AlgoRun> = ALGORITHMS.iter()
        .map(|algo| build_and_measure(cli, params, combo_dir, algo))
        .collect();

    log_measurements(&runs);

    let failed_algorithms: Vec<&'static str> =
        runs.iter().filter(|r| !r.success).map(|r| r.name).collect();
    let label = label.unwrap_or("(default parameters)").to_string();
    if !failed_algorithms.is_empty() {
        log::error!(
            "{} algorithm(s) failed to build: {}",
            failed_algorithms.len(),
            failed_algorithms.join(", ")
        );
        return ComboOutcome { label, failed_algorithms, mismatches: Vec::new() };
    }

    let successes: Vec<&AlgoRun> = runs.iter().filter(|r| r.success).collect();
    let mismatches = compare_outputs(&successes);
    if mismatches.is_empty() {
        log::info!("All {} construction algorithms agree.", successes.len());
    } else {
        log::error!("{} mismatch(es): {}", mismatches.len(), mismatches.join(", "));
    }

    ComboOutcome { label, failed_algorithms, mismatches }
}

/// Runs `sbwt build ...` (wrapped in `/usr/bin/time -v`) for one construction algorithm and
/// reports its success, wall-clock time and peak resident memory.
fn build_and_measure(cli: &Cli, params: &Params, combo_dir: &Path, algo: &Algorithm) -> AlgoRun {
    log::info!("Building with {} ...", algo.name);

    let prefix = combo_dir.join(algo.name);
    let mut cmd = Command::new(&cli.sbwt_bin);
    if cli.verbose {
        cmd.arg("-v");
    }
    cmd.arg("--threads").arg(params.threads.to_string())
        .arg("build")
        .arg(cli.input_arg).arg(&cli.input_path)
        .arg("-k").arg(params.k.to_string())
        .arg("--output-prefix").arg(&prefix)
        .arg("--build-lcs")
        .arg("--temp-dir").arg(combo_dir)
        .arg("--mem-gb").arg(cli.mem_gb.to_string());
    if params.add_revcomp {
        cmd.arg("--add-revcomp");
    }
    if params.add_all_dummy_paths {
        cmd.arg("--add-all-dummy-paths");
    }
    if params.dedup_batches {
        cmd.arg("--dedup-batches");
    }
    cmd.args(algo.extra_args);

    let run = common::run_timed(&cmd);

    if run.success {
        log::info!("Building with {} ... ok ({})", algo.name, common::format_duration(run.elapsed));
        if cli.verbose {
            log::info!("--- stderr for {} ---\n{}", algo.name, run.stderr);
        }
    } else {
        log::error!(
            "Building with {} ... FAILED (exit {:?}, {})",
            algo.name, run.status_code, common::format_duration(run.elapsed)
        );
        log::error!("--- stderr for {} ---\n{}", algo.name, run.stderr);
        if algo.name == "libsais" && run.stderr.contains("--features libsais") {
            log::error!(
                "Rebuild the sbwt binary with the libsais feature enabled: \
                 `cargo build --release --features \"tests libsais\"`."
            );
        }
    }

    let row = ReportRow {
        params: *params, algorithm: algo.name, success: run.success,
        elapsed: run.elapsed, peak_rss_kb: run.peak_rss_kb,
    };
    append_report_row(&cli.report_path, &row).expect("failed to append to report");

    AlgoRun { name: algo.name, prefix, success: run.success, elapsed: run.elapsed, peak_rss_kb: run.peak_rss_kb }
}

/// Byte-compares every successful run's `.sbwt` and `.lcs` output against the first
/// (baseline) run, logging a line per comparison. Returns a description of each mismatch.
fn compare_outputs(successes: &[&AlgoRun]) -> Vec<String> {
    let mut mismatches = Vec::new();
    for ext in ["sbwt", "lcs"] {
        let baseline = successes[0];
        let baseline_path = common::with_ext(&baseline.prefix, ext);
        let baseline_bytes = std::fs::read(&baseline_path)
            .unwrap_or_else(|e| panic!("failed to read {}: {e}", baseline_path.display()));

        for run in &successes[1..] {
            let path = common::with_ext(&run.prefix, ext);
            let bytes = std::fs::read(&path)
                .unwrap_or_else(|e| panic!("failed to read {}: {e}", path.display()));
            match first_difference(&baseline_bytes, &bytes) {
                None => log::info!(
                    "{:<21} .{ext} matches {} ({} bytes)",
                    run.name, baseline.name, bytes.len()
                ),
                Some(offset) => {
                    log::error!(
                        "{:<21} .{ext} MISMATCH vs {} \
                         (first differing byte at offset {offset}, sizes {} vs {})",
                        run.name, baseline.name, baseline_bytes.len(), bytes.len()
                    );
                    mismatches.push(format!("{} .{ext} vs {}", run.name, baseline.name));
                }
            }
        }
    }
    mismatches
}

fn log_measurements(runs: &[AlgoRun]) {
    log::info!("{:<21} {:>10} {:>14}", "algorithm", "time", "peak RSS");
    for r in runs {
        let time_str = common::format_duration(r.elapsed);
        let mem_str = match r.peak_rss_kb {
            Some(kb) => format!("{:.1} MiB", kb as f64 / 1024.0),
            None => "n/a".to_string(),
        };
        let name = if r.success { r.name.to_string() } else { format!("{} (failed)", r.name) };
        log::info!("{name:<21} {time_str:>10} {mem_str:>14}");
    }
}

/// Removes everything under `out_dir` except `report_path`, which is the whole point of a
/// --fuzz run and must survive even when `keep` is false.
fn cleanup(out_dir: &Path, keep: bool, report_path: &Path) {
    if keep {
        return;
    }
    let Ok(entries) = std::fs::read_dir(out_dir) else { return };
    for entry in entries.flatten() {
        let path = entry.path();
        if path == report_path {
            continue;
        }
        if path.is_dir() {
            let _ = std::fs::remove_dir_all(&path);
        } else {
            let _ = std::fs::remove_file(&path);
        }
    }
}

/// Creates (truncating any existing file) the report and writes just its header. Rows are
/// appended one at a time afterwards, by `append_report_row`, as each build finishes.
fn init_report(path: &Path) -> std::io::Result<()> {
    common::init_tsv_report(
        path,
        "k\tthreads\tadd_revcomp\tadd_all_dummy_paths\tdedup_batches\talgorithm\tsuccess\ttime_seconds\tpeak_rss_bytes"
    )
}

/// Appends one TSV row to the report: a build's parameters, whether it succeeded, wall-clock
/// time in seconds, and peak resident memory in bytes.
fn append_report_row(path: &Path, row: &ReportRow) -> std::io::Result<()> {
    let peak_rss_bytes = row.peak_rss_kb.map(|kb| kb * 1024);
    let line = format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{}",
        row.params.k,
        row.params.threads,
        row.params.add_revcomp,
        row.params.add_all_dummy_paths,
        row.params.dedup_batches,
        row.algorithm,
        row.success,
        row.elapsed.as_secs_f64(),
        peak_rss_bytes.map(|b| b.to_string()).unwrap_or_default(),
    );
    common::append_tsv_row(path, &line)
}

fn first_difference(a: &[u8], b: &[u8]) -> Option<usize> {
    if a == b {
        return None;
    }
    Some(a.iter().zip(b.iter()).position(|(x, y)| x != y).unwrap_or(a.len().min(b.len())))
}

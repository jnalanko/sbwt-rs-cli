use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::Duration;

use clap::parser::ValueSource;

use crate::common;

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

/// Defines the `build` subcommand's CLI schema.
pub fn subcommand() -> clap::Command {
    clap::Command::new("build")
        .about("Runs all sbwt construction algorithms through the sbwt CLI on the same input \
                and checks that they produce byte-identical output.")
        .arg(clap::Arg::new("input")
            .help("Input fasta or fastq sequence file (same as `sbwt build --input`).")
            .short('i')
            .long("input")
            .required_unless_present("input-list")
            .conflicts_with("input-list")
            .value_parser(clap::value_parser!(PathBuf)))
        .arg(clap::Arg::new("input-list")
            .help("File listing input fasta/fastq files, one path per line \
                   (same format as `sbwt build --input-list`).")
            .long("input-list")
            .required_unless_present("input")
            .value_parser(clap::value_parser!(PathBuf)))
        .arg(clap::Arg::new("k")
            .help("k-mer length")
            .short('k')
            .default_value("31")
            .value_parser(clap::value_parser!(usize)))
        .arg(clap::Arg::new("sbwt-bin")
            .help("Path to the sbwt executable to test. Defaults to the `sbwt` binary \
                   built alongside this harness.")
            .long("sbwt-bin")
            .value_parser(clap::value_parser!(PathBuf)))
        .arg(clap::Arg::new("out-dir")
            .help("Directory to write build outputs (and on-disk algorithm scratch files) to. \
                   Created if it doesn't exist.")
            .long("out-dir")
            .required(true)
            .value_parser(clap::value_parser!(PathBuf)))
        .arg(clap::Arg::new("keep")
            .help("Do not delete the output directory when done.")
            .long("keep")
            .action(clap::ArgAction::SetTrue))
        .arg(clap::Arg::new("report")
            .help("Path to write the TSV report (one row per build: parameters, time, peak \
                   memory) to. Defaults to construction-report.tsv inside --out-dir; a path \
                   given here is used as-is, so it doesn't have to be under --out-dir.")
            .long("report")
            .value_parser(clap::value_parser!(PathBuf)))
        .arg(clap::Arg::new("threads")
            .help("Number of threads to pass to the sbwt CLI.")
            .long("threads")
            .short('t')
            .default_value("4")
            .value_parser(clap::value_parser!(usize)))
        .arg(clap::Arg::new("mem-gb")
            .help("Memory budget in gigabytes to pass to the sbwt CLI's --mem-gb \
                   (used by the on-disk and in-memory algorithms only). Never fuzzed.")
            .long("mem-gb")
            .short('m')
            .default_value("8")
            .value_parser(clap::value_parser!(usize)))
        .arg(clap::Arg::new("add-revcomp")
            .help("Pass --add-revcomp to every build.")
            .long("add-revcomp")
            .action(clap::ArgAction::SetTrue))
        .arg(clap::Arg::new("add-all-dummy-paths")
            .help("Pass --add-all-dummy-paths to every build.")
            .long("add-all-dummy-paths")
            .action(clap::ArgAction::SetTrue))
        .arg(clap::Arg::new("dedup-batches")
            .help("Pass --dedup-batches to every build.")
            .long("dedup-batches")
            .short('d')
            .action(clap::ArgAction::SetTrue))
        .arg(clap::Arg::new("fuzz")
            .help("For each of -k, --threads, --add-revcomp, --add-all-dummy-paths and \
                   --dedup-batches that was not explicitly given, try multiple values \
                   instead of just the default, and build+compare every combination. \
                   --mem-gb is never fuzzed.")
            .long("fuzz")
            .action(clap::ArgAction::SetTrue))
        .arg(clap::Arg::new("verbose")
            .help("Pass -v to the sbwt CLI, and print its stderr for every build \
                   (not just failed ones).")
            .short('v')
            .long("verbose")
            .action(clap::ArgAction::SetTrue))
}

/// Validates parsed argv and works out the set of parameter combinations to run (a single
/// one, unless --fuzz expands some of them).
fn parse_args(matches: &clap::ArgMatches) -> Cli {
    let input = matches.get_one::<PathBuf>("input");
    let input_list = matches.get_one::<PathBuf>("input-list");
    let (input_arg, input_path) = match (input, input_list) {
        (Some(path), None) => ("--input", path.clone()),
        (None, Some(path)) => ("--input-list", path.clone()),
        _ => unreachable!("clap enforces exactly one of --input / --input-list"),
    };
    if !input_path.is_file() {
        eprintln!("error: input file {} does not exist", input_path.display());
        std::process::exit(2);
    }

    let mem_gb = *matches.get_one::<usize>("mem-gb").unwrap();
    let verbose = matches.get_flag("verbose");
    let fuzz = matches.get_flag("fuzz");

    // Only parameters the user left at their default get multiple candidate values; anything
    // passed explicitly is pinned to that one value. This lets a caller narrow the search
    // space to whatever they're actually unsure about instead of always paying for the full
    // combinatorial product.
    let explicit = |id: &str| matches.value_source(id) == Some(ValueSource::CommandLine);
    let k_candidates = candidates_usize(
        fuzz, explicit("k"), *matches.get_one::<usize>("k").unwrap(), &[3, 31, 32, 100],
    );
    let threads_candidates = candidates_usize(
        fuzz, explicit("threads"), *matches.get_one::<usize>("threads").unwrap(), &[1, 4],
    );
    let add_revcomp_candidates =
        candidates_bool(fuzz, explicit("add-revcomp"), matches.get_flag("add-revcomp"));
    let add_all_dummy_paths_candidates = candidates_bool(
        fuzz, explicit("add-all-dummy-paths"), matches.get_flag("add-all-dummy-paths"),
    );
    let dedup_batches_candidates =
        candidates_bool(fuzz, explicit("dedup-batches"), matches.get_flag("dedup-batches"));

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
    let keep = matches.get_flag("keep");

    let report_path = matches.get_one::<PathBuf>("report").cloned()
        .unwrap_or_else(|| out_dir.join("construction-report.tsv"));

    Cli { input_arg, input_path, mem_gb, verbose, sbwt_bin, out_dir, keep, combos, report_path }
}

/// Candidate values for a boolean build flag: both `false` and `true` if it's being fuzzed
/// (i.e. --fuzz was given and the user didn't pin this flag explicitly), otherwise just
/// whatever value was actually parsed.
fn candidates_bool(fuzz: bool, explicit: bool, actual: bool) -> Vec<bool> {
    if fuzz && !explicit { vec![false, true] } else { vec![actual] }
}

/// Same as [`candidates_bool`], but for a numeric option, with the fuzzed value set given
/// explicitly by the caller.
fn candidates_usize(fuzz: bool, explicit: bool, actual: usize, fuzzed: &[usize]) -> Vec<usize> {
    if fuzz && !explicit { fuzzed.to_vec() } else { vec![actual] }
}

pub fn run(matches: &clap::ArgMatches) {
    common::require_time_binary();
    let cli = parse_args(matches);

    println!("sbwt executable: {}", cli.sbwt_bin.display());
    println!("input ({}): {}", cli.input_arg, cli.input_path.display());
    println!("output dir:      {}", cli.out_dir.display());
    println!("mem-gb = {}", cli.mem_gb);
    if cli.combos.len() > 1 {
        println!("fuzzing {} parameter combination(s)\n", cli.combos.len());
    } else {
        println!("{}\n", cli.combos[0].label());
    }

    // Written to as each build finishes (see `append_report_row`), not batched up until the
    // end, so a partial report still exists if the harness is interrupted or a later build
    // hangs.
    init_report(&cli.report_path).expect("failed to create report");
    println!("report:          {} (appended to as each build finishes)\n", cli.report_path.display());

    let mut outcomes = Vec::new();
    for (i, params) in cli.combos.iter().enumerate() {
        let combo_dir = if cli.combos.len() > 1 { cli.out_dir.join(params.dir_name(i)) } else { cli.out_dir.clone() };
        let label = if cli.combos.len() > 1 { Some(params.label()) } else { None };
        outcomes.push(run_combo(&cli, params, &combo_dir, label.as_deref()));
    }

    cleanup(&cli.out_dir, cli.keep, &cli.report_path);

    if cli.combos.len() > 1 {
        println!("\n=== fuzz summary ===");
        for outcome in &outcomes {
            let status = if outcome.is_ok() { "OK" } else { "FAILED" };
            println!("{status:<6} {}", outcome.label);
        }
    }

    if outcomes.iter().all(ComboOutcome::is_ok) {
        println!("\nAll combination(s) agree across all construction algorithms.");
    } else {
        let failed = outcomes.iter().filter(|o| !o.is_ok()).count();
        eprintln!("\n{failed} of {} combination(s) failed.", outcomes.len());
        std::process::exit(1);
    }
}

/// Builds all four construction algorithms with the given parameters and byte-compares
/// their `.sbwt` and `.lcs` output. `label`, if given, is printed as a header (used when
/// running more than one parameter combination under --fuzz).
fn run_combo(cli: &Cli, params: &Params, combo_dir: &Path, label: Option<&str>) -> ComboOutcome {
    if let Some(label) = label {
        println!("--- {label} ---");
    }
    std::fs::create_dir_all(combo_dir).expect("failed to create combo output dir");

    let runs: Vec<AlgoRun> = ALGORITHMS.iter()
        .map(|algo| build_and_measure(cli, params, combo_dir, algo))
        .collect();

    print_measurements(&runs);

    let failed_algorithms: Vec<&'static str> =
        runs.iter().filter(|r| !r.success).map(|r| r.name).collect();
    let label = label.unwrap_or("(default parameters)").to_string();
    if !failed_algorithms.is_empty() {
        eprintln!(
            "\n{} algorithm(s) failed to build: {}",
            failed_algorithms.len(),
            failed_algorithms.join(", ")
        );
        return ComboOutcome { label, failed_algorithms, mismatches: Vec::new() };
    }

    println!();
    let successes: Vec<&AlgoRun> = runs.iter().filter(|r| r.success).collect();
    let mismatches = compare_outputs(&successes);
    if mismatches.is_empty() {
        println!("\nAll {} construction algorithms agree.", successes.len());
    } else {
        eprintln!("\n{} mismatch(es): {}", mismatches.len(), mismatches.join(", "));
    }

    ComboOutcome { label, failed_algorithms, mismatches }
}

/// Runs `sbwt build ...` (wrapped in `/usr/bin/time -v`) for one construction algorithm and
/// reports its success, wall-clock time and peak resident memory.
fn build_and_measure(cli: &Cli, params: &Params, combo_dir: &Path, algo: &Algorithm) -> AlgoRun {
    print!("Building with {:<21} ... ", algo.name);
    std::io::stdout().flush().ok();

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
        println!("ok ({})", common::format_duration(run.elapsed));
        if cli.verbose {
            eprintln!("--- stderr for {} ---\n{}", algo.name, run.stderr);
        }
    } else {
        println!("FAILED (exit {:?}, {})", run.status_code, common::format_duration(run.elapsed));
        eprintln!("--- stderr for {} ---\n{}", algo.name, run.stderr);
        if algo.name == "libsais" && run.stderr.contains("--features libsais") {
            eprintln!(
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
/// (baseline) run, printing a line per comparison. Returns a description of each mismatch.
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
                None => println!(
                    "{:<21} .{ext} matches {} ({} bytes)",
                    run.name, baseline.name, bytes.len()
                ),
                Some(offset) => {
                    println!(
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

fn print_measurements(runs: &[AlgoRun]) {
    println!("\n{:<21} {:>10} {:>14}", "algorithm", "time", "peak RSS");
    for r in runs {
        let time_str = common::format_duration(r.elapsed);
        let mem_str = match r.peak_rss_kb {
            Some(kb) => format!("{:.1} MiB", kb as f64 / 1024.0),
            None => "n/a".to_string(),
        };
        let name = if r.success { r.name.to_string() } else { format!("{} (failed)", r.name) };
        println!("{name:<21} {time_str:>10} {mem_str:>14}");
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

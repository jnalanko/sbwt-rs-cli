//! Standalone harness (not run by `cargo test`) that drives the `sbwt` CLI through all four
//! construction algorithms on the same input and checks that they all produce byte-identical
//! `.sbwt` and `.lcs` output. Each run is wrapped in `/usr/bin/time -v` to measure wall-clock
//! time and peak resident memory.
//!
//! Usage:
//!   cargo run --release --features libsais --bin sbwt-construction-test-harness -- --input <FILE> -k <K> --out-dir <DIR>
//!   cargo run --release --features libsais --bin sbwt-construction-test-harness -- --input-list <FILE> -k <K> --out-dir <DIR>
//!
//! With --fuzz, any of --add-revcomp, --add-all-dummy-paths, --dedup-batches and --threads
//! that were *not* explicitly given try multiple values instead of just their default, and
//! every resulting combination is built and compared separately. --mem-gb is never fuzzed.

use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::{Duration, Instant};

use clap::parser::ValueSource;

const TIME_BIN: &str = "/usr/bin/time";

struct Algorithm {
    name: &'static str,
    extra_args: &'static [&'static str],
}

const ALGORITHMS: &[Algorithm] = &[
    Algorithm { name: "on-disk", extra_args: &[] },
    Algorithm { name: "in-memory", extra_args: &["--in-memory"] },
    Algorithm { name: "bounded-suffix-sort", extra_args: &["--bounded-suffix-sort"] },
    Algorithm { name: "libsais", extra_args: &["--via-libsais"] },
];

/// The build parameters that --fuzz may vary across combinations.
#[derive(Clone, Copy)]
struct Params {
    threads: usize,
    add_revcomp: bool,
    add_all_dummy_paths: bool,
    dedup_batches: bool,
}

impl Params {
    fn label(&self) -> String {
        format!(
            "threads={} add-revcomp={} add-all-dummy-paths={} dedup-batches={}",
            self.threads, self.add_revcomp, self.add_all_dummy_paths, self.dedup_batches
        )
    }

    fn dir_name(&self, index: usize) -> String {
        format!(
            "combo-{index}_t{}_rc{}_adp{}_db{}",
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

/// Fully parsed and validated configuration for a harness run.
struct Cli {
    input_arg: &'static str,
    input_path: PathBuf,
    k: usize,
    mem_gb: usize,
    verbose: bool,
    sbwt_bin: PathBuf,
    out_dir: PathBuf,
    keep: bool,
    combos: Vec<Params>,
}

fn main() {
    let cli = parse_args();

    println!("sbwt executable: {}", cli.sbwt_bin.display());
    println!("input ({}): {}", cli.input_arg, cli.input_path.display());
    println!("output dir:      {}", cli.out_dir.display());
    println!("k = {}, mem-gb = {}", cli.k, cli.mem_gb);
    if cli.combos.len() > 1 {
        println!("fuzzing {} parameter combination(s)\n", cli.combos.len());
    } else {
        println!();
    }

    let mut outcomes = Vec::new();
    for (i, params) in cli.combos.iter().enumerate() {
        let combo_dir = if cli.combos.len() > 1 { cli.out_dir.join(params.dir_name(i)) } else { cli.out_dir.clone() };
        let label = if cli.combos.len() > 1 { Some(params.label()) } else { None };
        outcomes.push(run_combo(&cli, params, &combo_dir, label.as_deref()));
    }

    cleanup(&cli.out_dir, cli.keep);

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

/// Defines the CLI schema, parses argv, validates paths, and works out the set of parameter
/// combinations to run (a single one, unless --fuzz expands some of them).
fn parse_args() -> Cli {
    let matches = clap::Command::new("sbwt-construction-test-harness")
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
            .required(true)
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
            .help("For each of --threads, --add-revcomp, --add-all-dummy-paths and \
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
        .get_matches();

    if !Path::new(TIME_BIN).is_file() {
        eprintln!(
            "error: {TIME_BIN} not found. This harness uses GNU time to measure wall-clock \
             time and peak memory (e.g. `apt install time` on Debian/Ubuntu)."
        );
        std::process::exit(2);
    }

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

    let k = *matches.get_one::<usize>("k").unwrap();
    let mem_gb = *matches.get_one::<usize>("mem-gb").unwrap();
    let verbose = matches.get_flag("verbose");
    let fuzz = matches.get_flag("fuzz");

    let explicit = |id: &str| matches.value_source(id) == Some(ValueSource::CommandLine);
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
    for &threads in &threads_candidates {
        for &add_revcomp in &add_revcomp_candidates {
            for &add_all_dummy_paths in &add_all_dummy_paths_candidates {
                for &dedup_batches in &dedup_batches_candidates {
                    combos.push(Params { threads, add_revcomp, add_all_dummy_paths, dedup_batches });
                }
            }
        }
    }

    let sbwt_bin = matches.get_one::<PathBuf>("sbwt-bin").cloned()
        .unwrap_or_else(default_sbwt_bin_path);
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

    Cli { input_arg, input_path, k, mem_gb, verbose, sbwt_bin, out_dir, keep, combos }
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

/// Runs `/usr/bin/time -v <sbwt> build ...` for one construction algorithm and reports
/// its success, wall-clock time and peak resident memory.
fn build_and_measure(cli: &Cli, params: &Params, combo_dir: &Path, algo: &Algorithm) -> AlgoRun {
    print!("Building with {:<21} ... ", algo.name);
    std::io::stdout().flush().ok();

    let prefix = combo_dir.join(algo.name);
    let mut cmd = Command::new(TIME_BIN);
    // The first "-v" here is /usr/bin/time's own verbose flag, not sbwt's.
    cmd.arg("-v").arg(&cli.sbwt_bin);
    if cli.verbose {
        cmd.arg("-v");
    }
    cmd.arg("--threads").arg(params.threads.to_string())
        .arg("build")
        .arg(cli.input_arg).arg(&cli.input_path)
        .arg("-k").arg(cli.k.to_string())
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

    let wall_clock_start = Instant::now();
    let output = cmd.output().expect("failed to run /usr/bin/time");
    let stderr = String::from_utf8_lossy(&output.stderr);

    let success = output.status.success();
    let elapsed = parse_elapsed(&stderr).unwrap_or_else(|| wall_clock_start.elapsed());
    let peak_rss_kb = parse_max_rss_kb(&stderr);

    if success {
        println!("ok ({})", format_duration(elapsed));
        if cli.verbose {
            eprintln!("--- stderr for {} ---\n{stderr}", algo.name);
        }
    } else {
        println!("FAILED (exit {:?}, {})", output.status.code(), format_duration(elapsed));
        eprintln!("--- stderr for {} ---\n{stderr}", algo.name);
        if algo.name == "libsais" && stderr.contains("--features libsais") {
            eprintln!(
                "Rebuild the sbwt binary with the libsais feature enabled: \
                 `cargo build --release --features \"tests libsais\"`."
            );
        }
    }

    AlgoRun { name: algo.name, prefix, success, elapsed, peak_rss_kb }
}

/// Byte-compares every successful run's `.sbwt` and `.lcs` output against the first
/// (baseline) run, printing a line per comparison. Returns a description of each mismatch.
fn compare_outputs(successes: &[&AlgoRun]) -> Vec<String> {
    let mut mismatches = Vec::new();
    for ext in ["sbwt", "lcs"] {
        let baseline = successes[0];
        let baseline_path = with_ext(&baseline.prefix, ext);
        let baseline_bytes = std::fs::read(&baseline_path)
            .unwrap_or_else(|e| panic!("failed to read {}: {e}", baseline_path.display()));

        for run in &successes[1..] {
            let path = with_ext(&run.prefix, ext);
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

/// Parses the "Maximum resident set size (kbytes): N" line from `/usr/bin/time -v` output.
/// The value is whatever follows the last "): " on the line (the label's closing paren).
fn parse_max_rss_kb(time_v_output: &str) -> Option<u64> {
    let line = time_v_output.lines().find(|l| l.contains("Maximum resident set size"))?;
    line.rsplit("): ").next()?.trim().parse().ok()
}

/// Parses the "Elapsed (wall clock) time (h:mm:ss or m:ss): H:MM:SS" or "M:SS.ss" line
/// from `/usr/bin/time -v` output.
fn parse_elapsed(time_v_output: &str) -> Option<Duration> {
    let line = time_v_output.lines().find(|l| l.contains("Elapsed (wall clock) time"))?;
    let value = line.rsplit("): ").next()?.trim();

    let parts: Vec<&str> = value.split(':').collect();
    let mut seconds = 0f64;
    for (i, part) in parts.iter().rev().enumerate() {
        let unit: f64 = part.parse().ok()?;
        seconds += unit * 60f64.powi(i as i32);
    }
    Some(Duration::from_secs_f64(seconds))
}

fn print_measurements(runs: &[AlgoRun]) {
    println!("\n{:<21} {:>10} {:>14}", "algorithm", "time", "peak RSS");
    for r in runs {
        let time_str = format_duration(r.elapsed);
        let mem_str = match r.peak_rss_kb {
            Some(kb) => format!("{:.1} MiB", kb as f64 / 1024.0),
            None => "n/a".to_string(),
        };
        let name = if r.success { r.name.to_string() } else { format!("{} (failed)", r.name) };
        println!("{name:<21} {time_str:>10} {mem_str:>14}");
    }
}

fn format_duration(d: Duration) -> String {
    format!("{:.2}s", d.as_secs_f64())
}

fn cleanup(out_dir: &Path, keep: bool) {
    if !keep {
        let _ = std::fs::remove_dir_all(out_dir);
    }
}

fn with_ext(prefix: &Path, ext: &str) -> PathBuf {
    let mut s = prefix.as_os_str().to_os_string();
    s.push(".");
    s.push(ext);
    PathBuf::from(s)
}

fn first_difference(a: &[u8], b: &[u8]) -> Option<usize> {
    if a == b {
        return None;
    }
    Some(a.iter().zip(b.iter()).position(|(x, y)| x != y).unwrap_or(a.len().min(b.len())))
}

fn default_sbwt_bin_path() -> PathBuf {
    let mut path = std::env::current_exe().expect("failed to get current exe path");
    path.set_file_name(if cfg!(windows) { "sbwt.exe" } else { "sbwt" });
    path
}

//! Standalone harness (not run by `cargo test`) that drives the `sbwt` CLI through all four
//! construction algorithms on the same input and checks that they all produce byte-identical
//! `.sbwt` and `.lcs` output. Each run is wrapped in `/usr/bin/time -v` to measure wall-clock
//! time and peak resident memory.
//!
//! Usage:
//!   cargo run --release --features libsais --bin sbwt-construction-test-harness -- --input <FILE> -k <K>
//!   cargo run --release --features libsais --bin sbwt-construction-test-harness -- --input-list <FILE> -k <K>

use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::{Duration, Instant};

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

struct AlgoRun {
    name: &'static str,
    prefix: PathBuf,
    success: bool,
    elapsed: Duration,
    peak_rss_kb: Option<u64>,
}

fn main() {
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
            .help("Directory to write build outputs to. Defaults to a fresh temp directory \
                   that is deleted on success.")
            .long("out-dir")
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
        .arg(clap::Arg::new("add-revcomp")
            .help("Pass --add-revcomp to every build.")
            .long("add-revcomp")
            .action(clap::ArgAction::SetTrue))
        .arg(clap::Arg::new("add-all-dummy-paths")
            .help("Pass --add-all-dummy-paths to every build.")
            .long("add-all-dummy-paths")
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
        (Some(path), None) => ("--input", path),
        (None, Some(path)) => ("--input-list", path),
        _ => unreachable!("clap enforces exactly one of --input / --input-list"),
    };
    if !input_path.is_file() {
        eprintln!("error: input file {} does not exist", input_path.display());
        std::process::exit(2);
    }
    let k = *matches.get_one::<usize>("k").unwrap();
    let threads = *matches.get_one::<usize>("threads").unwrap();
    let add_revcomp = matches.get_flag("add-revcomp");
    let add_all_dummy_paths = matches.get_flag("add-all-dummy-paths");

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

    let out_dir_arg = matches.get_one::<PathBuf>("out-dir").cloned();
    let (out_dir, out_dir_is_temp) = match &out_dir_arg {
        Some(dir) => {
            std::fs::create_dir_all(dir).expect("failed to create --out-dir");
            (dir.clone(), false)
        }
        None => {
            let dir = std::env::temp_dir()
                .join(format!("sbwt-construction-harness-{}", std::process::id()));
            std::fs::create_dir_all(&dir).expect("failed to create temp output dir");
            (dir, true)
        }
    };
    let keep = matches.get_flag("keep") || !out_dir_is_temp;

    println!("sbwt executable: {}", sbwt_bin.display());
    println!("input ({input_arg}): {}", input_path.display());
    println!("output dir:      {}", out_dir.display());
    println!("k = {k}\n");

    let mut runs = Vec::new();

    for algo in ALGORITHMS {
        print!("Building with {:<21} ... ", algo.name);
        std::io::stdout().flush().ok();

        let prefix = out_dir.join(algo.name);
        let mut cmd = Command::new(TIME_BIN);
        cmd.arg("-v").arg(&sbwt_bin)
            .arg("--threads").arg(threads.to_string())
            .arg("build")
            .arg(input_arg).arg(input_path)
            .arg("-k").arg(k.to_string())
            .arg("--output-prefix").arg(&prefix)
            .arg("--build-lcs")
            .arg("--temp-dir").arg(&out_dir);
        if add_revcomp {
            cmd.arg("--add-revcomp");
        }
        if add_all_dummy_paths {
            cmd.arg("--add-all-dummy-paths");
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
        runs.push(AlgoRun { name: algo.name, prefix, success, elapsed, peak_rss_kb });
    }

    print_measurements(&runs);

    let failed: Vec<&str> = runs.iter().filter(|r| !r.success).map(|r| r.name).collect();
    if !failed.is_empty() {
        eprintln!("\n{} algorithm(s) failed to build: {}", failed.len(), failed.join(", "));
        cleanup(&out_dir, keep);
        std::process::exit(1);
    }

    println!();
    let successes: Vec<&AlgoRun> = runs.iter().filter(|r| r.success).collect();
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

    cleanup(&out_dir, keep);

    if mismatches.is_empty() {
        println!("\nAll {} construction algorithms agree.", successes.len());
    } else {
        eprintln!("\n{} mismatch(es): {}", mismatches.len(), mismatches.join(", "));
        std::process::exit(1);
    }
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

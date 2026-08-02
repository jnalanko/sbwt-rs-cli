//! Standalone harness (not run by `cargo test`) that drives the `sbwt` CLI through all four
//! construction algorithms on the same input and checks that they all produce byte-identical
//! `.sbwt` and `.lcs` output.
//!
//! Usage:
//!   cargo run --release --features libsais --bin sbwt-construction-test-harness -- --input <FILE> -k <K>
//!   cargo run --release --features libsais --bin sbwt-construction-test-harness -- --input-list <FILE> -k <K>

use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

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

    let mut prefixes = Vec::new();
    let mut failed = Vec::new();

    for algo in ALGORITHMS {
        print!("Building with {:<21} ... ", algo.name);
        std::io::stdout().flush().ok();

        let prefix = out_dir.join(algo.name);
        let mut cmd = Command::new(&sbwt_bin);
        cmd.arg("--threads").arg(threads.to_string())
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

        let output = cmd.output().expect("failed to run sbwt executable");
        if output.status.success() {
            println!("ok");
            prefixes.push((algo.name, prefix));
        } else {
            println!("FAILED (exit {:?})", output.status.code());
            eprintln!(
                "--- stderr for {} ---\n{}",
                algo.name,
                String::from_utf8_lossy(&output.stderr)
            );
            failed.push(algo.name);
        }
    }

    if !failed.is_empty() {
        eprintln!("\n{} algorithm(s) failed to build: {}", failed.len(), failed.join(", "));
        cleanup(&out_dir, keep);
        std::process::exit(1);
    }

    println!();
    let mut mismatches = Vec::new();
    for ext in ["sbwt", "lcs"] {
        let (baseline_name, baseline_prefix) = &prefixes[0];
        let baseline_path = with_ext(baseline_prefix, ext);
        let baseline_bytes = std::fs::read(&baseline_path)
            .unwrap_or_else(|e| panic!("failed to read {}: {e}", baseline_path.display()));

        for (name, prefix) in &prefixes[1..] {
            let path = with_ext(prefix, ext);
            let bytes = std::fs::read(&path)
                .unwrap_or_else(|e| panic!("failed to read {}: {e}", path.display()));
            match first_difference(&baseline_bytes, &bytes) {
                None => println!(
                    "{name:<21} .{ext} matches {baseline_name} ({} bytes)",
                    bytes.len()
                ),
                Some(offset) => {
                    println!(
                        "{name:<21} .{ext} MISMATCH vs {baseline_name} \
                         (first differing byte at offset {offset}, sizes {} vs {})",
                        baseline_bytes.len(),
                        bytes.len()
                    );
                    mismatches.push(format!("{name} .{ext} vs {baseline_name}"));
                }
            }
        }
    }

    cleanup(&out_dir, keep);

    if mismatches.is_empty() {
        println!("\nAll {} construction algorithms agree.", prefixes.len());
    } else {
        eprintln!("\n{} mismatch(es): {}", mismatches.len(), mismatches.join(", "));
        std::process::exit(1);
    }
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

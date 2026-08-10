use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::{Duration, Instant};

pub const TIME_BIN: &str = "/usr/bin/time";

/// Exits the process with an error if `/usr/bin/time` (used to measure wall-clock time and
/// peak memory for every build) isn't installed.
pub fn require_time_binary() {
    if !Path::new(TIME_BIN).is_file() {
        eprintln!(
            "error: {TIME_BIN} not found. This harness uses GNU time to measure wall-clock \
             time and peak memory (e.g. `apt install time` on Debian/Ubuntu)."
        );
        std::process::exit(2);
    }
}

/// Default location of the `sbwt` binary under test: alongside this harness's own executable.
pub fn default_sbwt_bin_path() -> PathBuf {
    let mut path = std::env::current_exe().expect("failed to get current exe path");
    path.set_file_name(if cfg!(windows) { "sbwt.exe" } else { "sbwt" });
    path
}

/// Appends `.{ext}` to `prefix` (e.g. an `--output-prefix` value) to get the actual output
/// file path.
pub fn with_ext(prefix: &Path, ext: &str) -> PathBuf {
    let mut s = prefix.as_os_str().to_os_string();
    s.push(".");
    s.push(ext);
    PathBuf::from(s)
}

pub fn format_duration(d: Duration) -> String {
    format!("{:.2}s", d.as_secs_f64())
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

/// Outcome of running a command wrapped in `/usr/bin/time -v`.
pub struct TimedRun {
    pub success: bool,
    pub status_code: Option<i32>,
    pub elapsed: Duration,
    pub peak_rss_kb: Option<u64>,
    pub stderr: String,
}

/// Runs `cmd` wrapped in `/usr/bin/time -v`, measuring wall-clock time and peak resident
/// memory. Falls back to a manually measured wall-clock time if `/usr/bin/time`'s own report
/// couldn't be parsed.
pub fn run_timed(cmd: &Command) -> TimedRun {
    let mut wrapped = Command::new(TIME_BIN);
    // This "-v" is /usr/bin/time's own verbose flag, not the wrapped command's.
    wrapped.arg("-v").arg(cmd.get_program()).args(cmd.get_args());
    if let Some(dir) = cmd.get_current_dir() {
        wrapped.current_dir(dir);
    }

    let wall_clock_start = Instant::now();
    let output = wrapped.output().expect("failed to run /usr/bin/time");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();

    let success = output.status.success();
    let elapsed = parse_elapsed(&stderr).unwrap_or_else(|| wall_clock_start.elapsed());
    let peak_rss_kb = parse_max_rss_kb(&stderr);

    TimedRun { success, status_code: output.status.code(), elapsed, peak_rss_kb, stderr }
}

/// Creates (truncating any existing file) `path` and writes `header` as its first line. Rows
/// are appended one at a time afterwards, by `append_tsv_row`, so a partial report survives
/// if the harness is interrupted or a later step hangs.
pub fn init_tsv_report(path: &Path, header: &str) -> std::io::Result<()> {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let mut out = std::fs::File::create(path)?;
    writeln!(out, "{header}")
}

/// Appends one already tab-separated `line` to the report at `path`. Opened and closed fresh
/// for each row (rather than keeping a file handle open across the whole run) so a row already
/// written is durable on disk even if the harness is later killed.
pub fn append_tsv_row(path: &Path, line: &str) -> std::io::Result<()> {
    let mut out = std::fs::OpenOptions::new().append(true).open(path)?;
    writeln!(out, "{line}")
}

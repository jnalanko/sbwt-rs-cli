use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};

pub const TIME_BIN: &str = "/usr/bin/time";

/// Exits the process with an error if `/usr/bin/time` (used to measure wall-clock time and
/// peak memory for every build) isn't installed.
pub fn require_time_binary() {
    if !Path::new(TIME_BIN).is_file() {
        log::error!(
            "{TIME_BIN} not found. This harness uses GNU time to measure wall-clock \
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
    /// Everything the command wrote to stderr. It has also already been echoed to the
    /// harness's own stderr as it arrived, so this is for inspecting, not for printing.
    pub stderr: String,
}

/// Runs `cmd` wrapped in `/usr/bin/time -v`, measuring wall-clock time and peak resident
/// memory. Falls back to a manually measured wall-clock time if `/usr/bin/time`'s own report
/// couldn't be parsed. `temp_dir` holds the scratch file that `/usr/bin/time` writes its
/// report to, which is removed again as soon as it has been read.
///
/// The command's output always goes to the harness's own: stdout is inherited, and stderr is
/// streamed through as it arrives, so a long-running command can be watched live rather than
/// only being reported on once it has finished. Stderr is captured into [`TimedRun::stderr`]
/// on the way past as well, for callers that want to look at what it said.
pub fn run_timed(cmd: &Command, temp_dir: &Path) -> TimedRun {
    // GNU time reports to its own stderr by default, which would both pollute the command's
    // captured stderr and bury the command's output in a screenful of resource statistics
    // after every step. "-o" sends it to a file instead, leaving the stderr stream as purely
    // the command's own.
    let report_path = temp_dir.join(format!(
        "time-report-{}-{}.txt",
        std::process::id(),
        NEXT_REPORT_ID.fetch_add(1, std::sync::atomic::Ordering::Relaxed)
    ));

    let mut wrapped = Command::new(TIME_BIN);
    // This "-v" is /usr/bin/time's own verbose flag, not the wrapped command's.
    wrapped.arg("-v").arg("-o").arg(&report_path)
        .arg(cmd.get_program()).args(cmd.get_args())
        .stdout(Stdio::inherit())
        .stderr(Stdio::piped());
    if let Some(dir) = cmd.get_current_dir() {
        wrapped.current_dir(dir);
    }

    let wall_clock_start = Instant::now();
    let mut child = wrapped.spawn().expect("failed to run /usr/bin/time");
    let stderr = tee_stderr(&mut child);
    let status = child.wait().expect("failed to wait for /usr/bin/time");

    let time_report = std::fs::read_to_string(&report_path).unwrap_or_default();
    let _ = std::fs::remove_file(&report_path);

    let elapsed = parse_elapsed(&time_report).unwrap_or_else(|| wall_clock_start.elapsed());
    let peak_rss_kb = parse_max_rss_kb(&time_report);

    TimedRun { success: status.success(), status_code: status.code(), elapsed, peak_rss_kb, stderr }
}

/// Distinguishes the scratch files of concurrent [`run_timed`] calls. Nothing runs commands in
/// parallel today, but the counter costs nothing and keeps the file name from being a shared
/// fixed path.
static NEXT_REPORT_ID: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);

/// Drains `child`'s piped stderr to EOF, returning everything it wrote. Each chunk is also
/// written straight through to the harness's stderr as it arrives — unbuffered and unprefixed,
/// so the command's output looks the same as if it had been run directly.
fn tee_stderr(child: &mut std::process::Child) -> String {
    let mut reader = child.stderr.take().expect("stderr was piped");
    let mut captured: Vec<u8> = Vec::new();
    let mut buf = [0u8; 8192];
    loop {
        let n = match reader.read(&mut buf) {
            Ok(0) | Err(_) => break,
            Ok(n) => n,
        };
        captured.extend_from_slice(&buf[..n]);
        let mut out = std::io::stderr();
        let _ = out.write_all(&buf[..n]);
        let _ = out.flush();
    }
    String::from_utf8_lossy(&captured).into_owned()
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

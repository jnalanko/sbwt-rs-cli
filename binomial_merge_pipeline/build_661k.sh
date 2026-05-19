#!/usr/bin/env bash
# =============================================================================
# ENA 661K SBWT Build Pipeline — online binomial merge
# =============================================================================
# Builds an SBWT index of all ~661K assemblies from the ENA2018-bacteria-661k
# dataset by:
#   1. Parsing sampleid_assembly_paths.txt to generate one batch list per
#      natural batch_XXX directory (already split by the dataset providers)
#   2. For each batch: download FASTAs → build SBWT → delete FASTAs →
#      immediately push into a merge stack (binomial-heap style):
#        - if merge_stack[L] is empty, store the new SBWT there
#        - if merge_stack[L] is occupied, merge the two → carry to L+1,
#          deleting both inputs
#   3. After all batches, flush the remaining stack entries into one final index
#
# Peak disk usage:
#   At most ceil(log2(N_batches)) stack slots, each 2× the size of the previous,
#   plus one batch being built and one merge in progress.
#
# Requirements:
#   - wget in PATH
#   - sbwt binary at the path set in SBWT= below
#
# Manifest (download once if missing):
#   wget -O /path/to/sampleid_assembly_paths.txt \
#        https://ftp.ebi.ac.uk/pub/databases/ENA2018-bacteria-661k/sampleid_assembly_paths.txt
# =============================================================================

set -euo pipefail
# NOTE: the FATAL trap uses log() which requires LOG_DIR to exist.  LOG_DIR is
# created by the mkdir -p below before any command that could trigger this trap.
trap '_exit_log "FATAL: command failed at line $LINENO: $BASH_COMMAND"' ERR
_exit_log() { echo "$*" >&2; [[ -n "${LOG_DIR:-}" ]] && echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >> "${LOG_DIR}/pipeline.log" || true; }

# ── Configuration ─────────────────────────────────────────────────────────────
# Set these three variables before running.
MANIFEST="/path/to/sampleid_assembly_paths.txt"   # ENA 661K manifest (see download command in header)
SBWT="${SBWT:-/path/to/sbwt-set-operations/target/release/sbwt}"
WORKDIR="/path/to/output"

K=31
MEM_GB=48
THREADS=32
DOWNLOAD_PARALLELISM=8    # simultaneous wget calls per batch; EBI FTP throttles at high concurrency

# EBI FTP base URL; assembly paths in the manifest start with /ebi/ftp/pub/...
# which maps directly to https://ftp.ebi.ac.uk/pub/...
EBI_FTP_BASE="https://ftp.ebi.ac.uk"
# The manifest paths start with /ebi/ftp — strip that prefix to get the URL path.
EBI_FTP_STRIP="/ebi/ftp"

LISTS_DIR="$WORKDIR/lists"
LOCAL_LISTS_DIR="$WORKDIR/local_lists"   # separate dir so these don't match batch_*.txt glob
FASTA_DIR="$WORKDIR/fasta"
MERGE_DIR="$WORKDIR/merge"          # all intermediate and stack SBWTs live here
LOG_DIR="$WORKDIR/logs"
FINAL_INDEX="$WORKDIR/661k_k${K}"

BENCH_TSV="$LOG_DIR/benchmarks.tsv"  # structured per-event benchmark table

# ── Setup ─────────────────────────────────────────────────────────────────────
mkdir -p "$LISTS_DIR" "$LOCAL_LISTS_DIR" "$FASTA_DIR" "$MERGE_DIR" "$LOG_DIR"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG_DIR/pipeline.log" >&2; }

# Kill any background downloads on INT/TERM.
_pipeline_cleanup() {
    log "Signal received — killing background downloads..."
    [[ -n "${prefetch_pid:-}"  ]] && kill "$prefetch_pid"  2>/dev/null || true
    [[ -n "${prefetch_pid2:-}" ]] && kill "$prefetch_pid2" 2>/dev/null || true
}
trap _pipeline_cleanup INT TERM

# Append one TSV row: timestamp, event, id, duration_s, output_bytes, n_inputs
log_bench() {
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$(date '+%Y-%m-%d %H:%M:%S')" "$1" "$2" "$3" "$4" "$5" >> "$BENCH_TSV"
}

# Write TSV header (skip if resuming an existing run)
if [[ ! -f "$BENCH_TSV" ]]; then
    printf 'timestamp\tevent\tid\tduration_s\toutput_bytes\tn_inputs\n' > "$BENCH_TSV"
fi

# ── Pre-flight checks ──────────────────────────────────────────────────────────
preflight_ok=1
[[ -x "$SBWT" ]]             || { echo "ERROR: SBWT binary not found or not executable: $SBWT"; preflight_ok=0; }
command -v wget >/dev/null   || { echo "ERROR: wget not in PATH"; preflight_ok=0; }
[[ -r "$MANIFEST" ]]         || { echo "ERROR: manifest not readable: $MANIFEST"; preflight_ok=0; }
[[ "$preflight_ok" -eq 1 ]]  || exit 1

# ── Resource limits ────────────────────────────────────────────────────────────
# Memory caps enforced via systemd cgroups. Adjust to fit your machine.
BUILD_MEM_LIMIT="80G"
MERGE_MEM_LIMIT="200G"

TIMINGS_LOG="$LOG_DIR/timings.log"

# Wrap a command with systemd memory cgroup + /usr/bin/time -v for benchmarking.
# Usage: run_timed <mem_limit> <label> <logfile> <cmd> [args...]
#   - command stdout+stderr → <logfile>
#   - /usr/bin/time -v stats → $TIMINGS_LOG  (via -o/--append, no fd juggling)
run_timed() {
    local mem_limit="$1" label="$2" logfile="$3"
    shift 3
    echo "=== $label ===" >> "$TIMINGS_LOG"

    local exit_code=0
    systemd-run --user --scope -p MemoryMax="$mem_limit" \
        -p MemorySwapMax=0B \
        -- /usr/bin/time -v -o "$TIMINGS_LOG" --append -- "$@" >> "$logfile" 2>&1 || exit_code=$?

    if [[ $exit_code -ne 0 ]]; then
        echo "[CRASH] $label failed with exit code $exit_code" | tee -a "$TIMINGS_LOG" >> "$logfile"
        echo "--- Diagnostics for $label ---" >> "$logfile"
        # dmesg may be restricted (Operation not permitted) — redirect stderr silently
        dmesg 2>/dev/null | tail -n 50 | grep -iE "oom|out of memory|killed" >> "$logfile" 2>/dev/null || true
        # Check user journal for OOM/kill events
        journalctl --user -n 100 --since "1 minute ago" >> "$logfile" 2>/dev/null || true
    fi

    return $exit_code
}

log "Starting 661K SBWT build pipeline (online binomial merge)"
log "  k=$K  mem_gb=$MEM_GB  threads=$THREADS"
log "  SBWT binary: $SBWT"

# Safely remove a directory, refusing to proceed if it's not under WORKDIR.
safe_rm_dir() {
    local dir="$1"
    if [[ "$dir" != "$WORKDIR/"* ]]; then
        log "ERROR: refusing to rm -rf '$dir' (not under WORKDIR '$WORKDIR')"
        exit 1
    fi
    rm -rf -- "$dir"
}

# ── Merge stack ────────────────────────────────────────────────────────────────
# merge_stack[L] holds the path to the SBWT sitting at level L, or "" if empty.
# Level 0 holds one batch, level 1 holds a merge of two batches, etc.
declare -a merge_stack=()
merge_counter=0  # monotonic counter for unique output filenames
MERGE_OUT=""     # global variable to capture merge result without subshells

# Merge two SBWTs, delete both inputs, set path of result to global MERGE_OUT.
do_merge() {
    local a="$1" b="$2" level="$3"
    merge_counter=$((merge_counter + 1))
    local out_prefix="$MERGE_DIR/stack_l${level}_$(printf '%06d' $merge_counter)"
    local in_bytes
    in_bytes=$(( $(du -sb "$a" | cut -f1) + $(du -sb "$b" | cut -f1) ))
    log "  [merge L${level}] $(basename "$a") + $(basename "$b") -> $(basename "${out_prefix}.sbwt") (inputs: $(numfmt --to=iec $in_bytes))"
    # Use --low-ram for level >= 5 (large indexes) and all flush merges.
    # Levels 0-4 are small enough to merge without the overhead.
    local low_ram_flag=()
    if [[ "$level" =~ ^[0-9]+$ ]] && (( level >= 5 )); then
        low_ram_flag=(--low-ram)
    elif [[ "$level" =~ ^flush ]]; then
        low_ram_flag=(--low-ram)
    fi
    local t0; t0=$(date +%s)
    run_timed "$MERGE_MEM_LIMIT" "merge:L${level}:$(basename "$out_prefix")" "$LOG_DIR/merge.log" \
        "$SBWT" merge \
            "${low_ram_flag[@]}" \
            --threads "$THREADS" \
            "$a" "$b" \
            -o "$out_prefix" || true
    # systemd-run reports exit 0 even when the process is killed by SIGKILL (OOM).
    # The output file may exist but be empty/truncated after a kill — use non-zero
    # file size as the authoritative success indicator.
    _merge_ok() { [[ -s "${out_prefix}.sbwt" ]]; }
    if ! _merge_ok; then
        rm -f "${out_prefix}.sbwt"
        if [[ ${#low_ram_flag[@]} -eq 0 ]]; then
            log "  [merge L${level}] output missing or empty (likely OOM-killed) — retrying with --low-ram"
            run_timed "$MERGE_MEM_LIMIT" "merge:L${level}:$(basename "$out_prefix").retry" "$LOG_DIR/merge.log" \
                "$SBWT" merge \
                    --low-ram \
                    --threads "$THREADS" \
                    "$a" "$b" \
                    -o "$out_prefix" || true
        fi
        # After retry (or if --low-ram was already set), check again.
        if ! _merge_ok; then
            log "WARNING: merge produced no output at ${out_prefix}.sbwt — inputs preserved (non-fatal, will keep both on stack)"
            unset -f _merge_ok
            MERGE_OUT=""
            return 1
        fi
    fi
    unset -f _merge_ok
    local t1; t1=$(date +%s)
    local out_bytes
    out_bytes=$(du -sb "${out_prefix}.sbwt" | cut -f1)
    log "  [merge L${level}] done in $((t1 - t0))s — output $(numfmt --to=iec $out_bytes)"
    log_bench "merge_L${level}" "$(basename "$out_prefix")" "$((t1 - t0))" "$out_bytes" "2"

    rm -f "$a" "$b"

    MERGE_OUT="${out_prefix}.sbwt"
}

# Push a new SBWT into the merge stack, cascading merges upward as needed.
# If a merge fails (OOM), both inputs are preserved at separate empty slots
# so the pipeline can continue; they can be merged manually later.
push_to_stack() {
    local new_sbwt="$1"
    local level=0
    while true; do
        if [[ -z "${merge_stack[$level]:-}" ]]; then
            merge_stack[$level]="$new_sbwt"
            log "  [stack] stored at level $level: $(basename "$new_sbwt")"
            break
        else
            local existing="${merge_stack[$level]}"
            merge_stack[$level]=""
            if ! do_merge "$existing" "$new_sbwt" "$level"; then
                # Merge failed (likely OOM). Place both inputs at the first two
                # empty slots at or above level+1 so processing can continue.
                # They will either be merged later by flush_stack or manually.
                log "  [stack] merge L${level} failed — parking both inputs on stack at escalated levels"
                local park=$((level + 1))
                while [[ -n "${merge_stack[$park]:-}" ]]; do park=$((park + 1)); done
                merge_stack[$park]="$existing"
                log "  [stack] parked $(basename "$existing") at level $park"
                park=$((park + 1))
                while [[ -n "${merge_stack[$park]:-}" ]]; do park=$((park + 1)); done
                merge_stack[$park]="$new_sbwt"
                log "  [stack] parked $(basename "$new_sbwt") at level $park"
                break
            fi
            new_sbwt="$MERGE_OUT"
            level=$((level + 1))
        fi
    done
}

# Flush remaining stack entries using a balanced tournament merge.
# Collects non-empty slots smallest-first, then repeatedly pairs adjacent entries
# until one remains — avoiding O(N^2) write amplification from linear accumulation.
flush_stack() {
    # Collect non-empty entries, smallest (lowest level) first
    local -a items=()
    local level
    for level in $(echo "${!merge_stack[@]}" | tr ' ' '\n' | sort -n); do
        local item="${merge_stack[$level]}"
        [[ -z "$item" ]] && continue
        items+=("$item")
        merge_stack[$level]=""
    done

    if [[ ${#items[@]} -eq 0 ]]; then
        echo ""
        return
    fi

    # Tournament: merge adjacent pairs, repeat until one item remains
    local round=0
    while [[ ${#items[@]} -gt 1 ]]; do
        round=$((round + 1))
        local -a next=()
        local i=0
        while [[ $i -lt ${#items[@]} ]]; do
            if [[ $((i + 1)) -lt ${#items[@]} ]]; then
                do_merge "${items[$i]}" "${items[$((i+1))]}" "flush${round}"
                next+=("$MERGE_OUT")
                i=$((i + 2))
            else
                next+=("${items[$i]}")
                i=$((i + 1))
            fi
        done
        items=("${next[@]}")
    done

    echo "${items[0]}"
}

# ── Step 1: Generate batch file lists from the manifest ───────────────────────
# The manifest has two tab-separated columns:
#   SAMPLE_ID   /ebi/ftp/pub/databases/ENA2018-bacteria-661k/Assemblies/batch_XXX/SAMPLE.contigs.fa.gz
# We group by the batch_XXX token and write one URL-per-line list file per batch.
BATCH_DONE_MARKER="$LISTS_DIR/.lists_generated"
if [[ ! -f "$BATCH_DONE_MARKER" ]]; then
    log "Generating batch lists from manifest..."
    # Extract batch id (e.g. batch_000) and full FTP URL from each line.
    # The /ebi/ftp prefix is stripped and replaced with $EBI_FTP_BASE.
    awk -v base="$EBI_FTP_BASE" -v strip="$EBI_FTP_STRIP" -v lists_dir="$LISTS_DIR" '
        NF >= 2 {
            path = $2
            # Extract batch_XXX from path
            match(path, /batch_[0-9]+/)
            batch = substr(path, RSTART, RLENGTH)
            # Build FTP URL: strip leading /ebi/ftp, prepend base
            sub("^" strip, "", path)
            url = base path
            print url > (lists_dir "/" batch ".txt")
        }
    ' "$MANIFEST"

    N_BATCHES=$(ls "$LISTS_DIR"/batch_*.txt | wc -l)
    log "Generated $N_BATCHES batch lists from manifest"
    touch "$BATCH_DONE_MARKER"
else
    N_BATCHES=$(ls "$LISTS_DIR"/batch_*.txt | wc -l)
    log "Batch lists already exist ($N_BATCHES batches), skipping generation"
fi

# ── Step 2: Build batches and merge on-the-fly ─────────────────────────────────
log "Building batches and merging on-the-fly..."

# Download a single FASTA from EBI FTP; called in parallel by download_batch.
# Appends the local path to local_listfile on success (O_APPEND writes are
# atomic on Linux for single lines, so no locking needed).
_download_one() {
    local url="$1" fasta_dir="$2" local_listfile="$3" logfile="$4"
    local filename local_path
    filename=$(basename "$url")
    local_path="$fasta_dir/$filename"
    # Retry up to 5 times with exponential back-off (4s, 8s, 16s, 32s).
    # wget handles its own per-attempt timeout; shell loop retries on any failure.
    local attempt
    for attempt in 1 2 3 4 5; do
        if wget -q --tries=3 --timeout=300 --waitretry=10 -O "$local_path" "$url" 2>> "$logfile"; then
            echo "$local_path" >> "$local_listfile"
            return 0
        fi
        echo "WARNING: download attempt $attempt/5 failed for $url" >> "$logfile"
        [[ $attempt -lt 5 ]] && sleep $(( 4 * 2 ** (attempt - 1) ))
        rm -f "$local_path"
    done
    echo "ERROR: all 5 download attempts failed for $url" >> "$logfile"
    rm -f "$local_path"
}
export -f _download_one

# Download all FASTAs for one batch into fasta_dir in parallel, writing local
# paths to local_listfile. Runs synchronously; call as `download_batch ... &`
# to background.
download_batch() {
    local listfile="$1" fasta_dir="$2" local_listfile="$3" logfile="$4"
    mkdir -p "$fasta_dir"
    : > "$local_listfile"
    # Each line of listfile is an EBI FTP URL; -P spawns DOWNLOAD_PARALLELISM
    # workers simultaneously. bash -c trampoline is needed so the exported
    # function is visible in the child shell.
    xargs -a "$listfile" -P "$DOWNLOAD_PARALLELISM" -I{} \
        bash -c '_download_one "$@"' _ {} "$fasta_dir" "$local_listfile" "$logfile"
}

# Fill the 2-slot prefetch queue with the next unbuilt batches at or after index $1.
# Two slots in flight hide download latency: while one batch builds, the next
# two batches download concurrently.
fill_prefetch_queue() {
    local start_i="$1"
    for (( ni = start_i; ni < N_BATCHES; ni++ )); do
        local nl="${ALL_LISTS[$ni]}"
        local nid
        nid=$(basename "$nl" .txt)

        # Skip already-built or already-being-fetched batches
        if [[ -f "$MERGE_DIR/${nid}.sbwt" ]]; then continue; fi
        if [[ "$prefetch_for"  == "$nid" || "$prefetch_for2" == "$nid" ]]; then continue; fi

        if [[ -z "$prefetch_pid" ]]; then
            log "[PREFETCH] Background download of $nid starting (slot 1)"
            download_batch "$nl" "$FASTA_DIR/$nid" \
                "$LOCAL_LISTS_DIR/${nid}_local.txt" "$LOG_DIR/${nid}.log" &
            prefetch_pid=$!
            prefetch_for="$nid"
        elif [[ -z "$prefetch_pid2" ]]; then
            log "[PREFETCH] Background download of $nid starting (slot 2)"
            download_batch "$nl" "$FASTA_DIR/$nid" \
                "$LOCAL_LISTS_DIR/${nid}_local.txt" "$LOG_DIR/${nid}.log" &
            prefetch_pid2=$!
            prefetch_for2="$nid"
            return 0  # both slots now filled
        else
            return 0  # both slots already occupied
        fi
    done
}

# Collect all batch list files into an indexed array for look-ahead prefetch.
mapfile -t ALL_LISTS < <(ls -1 "$LISTS_DIR"/batch_*.txt | sort)
N_BATCHES=${#ALL_LISTS[@]}

prefetch_pid=""    # PID of slot-1 background download (or "")
prefetch_for=""    # batch_id for slot-1 download
prefetch_pid2=""   # PID of slot-2 background download (or "")
prefetch_for2=""   # batch_id for slot-2 download

# Kick off downloads for the first two unbuilt batches immediately.
fill_prefetch_queue 0

for (( i = 0; i < N_BATCHES; i++ )); do
    listfile="${ALL_LISTS[$i]}"
    batch_id=$(basename "$listfile" .txt)
    batch_fasta_dir="$FASTA_DIR/${batch_id}"
    local_listfile="$LOCAL_LISTS_DIR/${batch_id}_local.txt"
    batch_sbwt_prefix="$MERGE_DIR/${batch_id}"
    batch_sbwt="${batch_sbwt_prefix}.sbwt"
    logfile="$LOG_DIR/${batch_id}.log"

    if [[ ! -f "$batch_sbwt" ]]; then
        log "[START] $batch_id"

        # ── Wait for this batch's download ───────────────────────────────────
        dl_t0=0 dl_t1=0
        if [[ "$prefetch_for" == "$batch_id" && -n "$prefetch_pid" ]]; then
            log "[WAIT] Waiting for download of $batch_id to finish (slot 1)..."
            dl_t0=$(date +%s)
            wait "$prefetch_pid" || true
            dl_t1=$(date +%s)
            prefetch_pid=""; prefetch_for=""
        elif [[ "$prefetch_for2" == "$batch_id" && -n "$prefetch_pid2" ]]; then
            log "[WAIT] Waiting for download of $batch_id to finish (slot 2)..."
            dl_t0=$(date +%s)
            wait "$prefetch_pid2" || true
            dl_t1=$(date +%s)
            prefetch_pid2=""; prefetch_for2=""
        else
            # No prefetch slot was aimed at this batch; download synchronously.
            log "[DOWNLOAD] $batch_id (synchronous)"
            dl_t0=$(date +%s)
            download_batch "$listfile" "$batch_fasta_dir" "$local_listfile" "$logfile"
            dl_t1=$(date +%s)
        fi

        n_downloaded=$(wc -l < "$local_listfile")
        dl_bytes=$(du -sb "$batch_fasta_dir" 2>/dev/null | cut -f1 || echo 0)
        log "[DOWNLOADED] $batch_id: $n_downloaded files, $(numfmt --to=iec ${dl_bytes:-0}), in $((dl_t1 - dl_t0))s"
        log_bench "download" "$batch_id" "$((dl_t1 - dl_t0))" "${dl_bytes:-0}" "$n_downloaded"
        if [[ "$n_downloaded" -eq 0 ]]; then
            log "[WARN] $batch_id: no files downloaded, skipping"
            echo "$batch_id" >> "$LOG_DIR/skipped_batches.log"
            safe_rm_dir "$batch_fasta_dir"
            rm -f -- "$local_listfile"
            # Start prefetch for the next unbuilt batch before continuing.
            fill_prefetch_queue $((i + 1))
            continue
        fi

        # Fill prefetch queue right after download completes, before the build
        # starts, so the next batches download concurrently with this build.
        fill_prefetch_queue $((i + 1))

        # ── Build SBWT ───────────────────────────────────────────────────────
        build_t0=$(date +%s)
        build_args=(
            "$SBWT" build
            --input-list "$local_listfile"
            --output-prefix "$batch_sbwt_prefix"
            -k "$K"
            --mem-gb "$MEM_GB"
            --threads "$THREADS"
            --add-revcomp
            --dedup-batches
            --in-memory
            -v
        )
        run_timed "$BUILD_MEM_LIMIT" "build:$batch_id" "$logfile" "${build_args[@]}"
        build_t1=$(date +%s)

        batch_bytes=$(du -sb "$batch_sbwt" | cut -f1)
        safe_rm_dir "$batch_fasta_dir"
        rm -f -- "$local_listfile"
        log "[BUILT] $batch_id in $((build_t1 - build_t0))s — $(numfmt --to=iec $batch_bytes)"
        log_bench "build" "$batch_id" "$((build_t1 - build_t0))" "$batch_bytes" "$n_downloaded"
    else
        log "[SKIP] $batch_id already built"
        # Fill prefetch queue for future batches.
        fill_prefetch_queue $((i + 1))
    fi

    # Push into online merge stack (deletes batch_sbwt when consumed)
    push_to_stack "$batch_sbwt"
done

# Clean up any lingering prefetch slots.
[[ -n "$prefetch_pid"  ]] && { wait "$prefetch_pid"  || true; }
[[ -n "$prefetch_pid2" ]] && { wait "$prefetch_pid2" || true; }

log "All batches processed. Flushing merge stack..."

# ── Step 3: Flush remaining stack entries ────────────────────────────────────
if [[ -f "$LOG_DIR/skipped_batches.log" ]]; then
    n_skipped=$(wc -l < "$LOG_DIR/skipped_batches.log")
    log "WARNING: $n_skipped batch(es) were skipped due to download failures — see $LOG_DIR/skipped_batches.log"
fi

RESULT=$(flush_stack)

# ── Final output ───────────────────────────────────────────────────────────────
log "Merge complete. Moving final index to $FINAL_INDEX.sbwt"
mv "$RESULT" "${FINAL_INDEX}.sbwt"

FINAL_BYTES=$(du -sb "${FINAL_INDEX}.sbwt" | cut -f1)
FINAL_SIZE=$(numfmt --to=iec "$FINAL_BYTES")
log "Final index: ${FINAL_INDEX}.sbwt ($FINAL_SIZE)"
log_bench "final" "661k_k${K}" "0" "$FINAL_BYTES" ""

log "Pipeline complete."
#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# FASTQ download pipeline from SRA accessions.
#
# Steps:
#   1) PREFETCH      - download .sra files if they are missing
#   2) VDB-VALIDATE  - validate downloaded .sra files
#   3) FASTQ-DUMP    - convert .sra files to FASTQ
#   4) SKIP          - skip runs whose FASTQ files already exist
#
# Expected project structure:
#
# project_directory/
# ├── samples.txt        # one accession per line; comments with "#" are allowed
# └── sra_files/
#     ├── ERR1308583/
#     │   └── ERR1308583.sra
#     └── ...
#
# Notes:
#   - Root access is not required on the host system.
#   - Docker is run with the current user UID:GID to avoid root-owned files.
#   - The current working directory is mounted as /work inside the container.
#
# Usage:
#   bash scripts/01_download_fastq_sra.sh [samples.txt] [max_jobs]
#
# Examples:
#   bash scripts/01_download_fastq_sra.sh
#   bash scripts/01_download_fastq_sra.sh samples.txt 3
# ============================================================

# ------------------------- CONFIG ---------------------------

BASEDIR="$PWD"                                  # Directory where the script is executed
RUNS_FILE="${1:-$BASEDIR/samples.txt}"          # Run accession list; one accession per line
SRA_DIR="$BASEDIR/sra_files"                    # Directory containing or receiving .sra files
OUTDIR="$BASEDIR/fastq"                         # FASTQ output directory
TMPBASE="$BASEDIR/tmp"                          # Temporary directory, organized by run
LOGDIR="$BASEDIR/logs"                          # Log directory, one log file per run
MAX_JOBS="${2:-3}"                              # Number of runs processed in parallel
IMAGE="pegi3s/sratoolkit"                       # Docker image containing SRA Toolkit

# ------------------------------------------------------------

mkdir -p "$SRA_DIR" "$OUTDIR" "$TMPBASE" "$LOGDIR"

if [[ ! -f "$RUNS_FILE" ]]; then
  echo "[ERROR] Run accession file not found: $RUNS_FILE" >&2
  exit 1
fi

# Remove blank lines and comments from the run accession file.
RUNS_CLEAN="$(mktemp)"
trap 'rm -f "$RUNS_CLEAN"' EXIT

grep -vE '^\s*($|#)' "$RUNS_FILE" > "$RUNS_CLEAN"

if [[ ! -s "$RUNS_CLEAN" ]]; then
  echo "[ERROR] No valid run accessions found in: $RUNS_FILE" >&2
  exit 1
fi

UIDGID="$(id -u):$(id -g)"

echo "BASE_DIR: $BASEDIR"
echo "RUNS:     $(wc -l < "$RUNS_CLEAN")"
echo "SRA_DIR:  $SRA_DIR"
echo "OUTDIR:   $OUTDIR"
echo "TMPBASE:  $TMPBASE"
echo "LOGDIR:   $LOGDIR"
echo "MAX_JOBS: $MAX_JOBS"
echo "IMAGE:    $IMAGE"
echo

run_one() {
  local RUN="$1"
  local run_log="$LOGDIR/${RUN}.log"
  local tmpdir="$TMPBASE/$RUN"
  local sra_path="$SRA_DIR/$RUN/$RUN.sra"

  local out1="$OUTDIR/${RUN}_1.fastq"
  local out2="$OUTDIR/${RUN}_2.fastq"
  local out1_gz="${out1}.gz"
  local out2_gz="${out2}.gz"

  mkdir -p "$tmpdir"

  {
    echo "==============================="
    echo "RUN:   $RUN"
    echo "START: $(date)"
    echo "TMP:   $tmpdir"
    echo "SRA:   $sra_path"
    echo

    # Skip this run if paired FASTQ files already exist.
    if [[ -s "$out1" && -s "$out2" ]]; then
      echo "[SKIP] FASTQ files already exist: $out1 and $out2"
      echo "END: $(date)"
      return 0
    fi

    # Also skip if compressed paired FASTQ files already exist.
    if [[ -s "$out1_gz" && -s "$out2_gz" ]]; then
      echo "[SKIP] Compressed FASTQ files already exist: $out1_gz and $out2_gz"
      echo "END: $(date)"
      return 0
    fi

    # Step 1: download the .sra file only if it is missing.
    if [[ ! -s "$sra_path" ]]; then
      echo "[1/4] PREFETCH: downloading .sra file..."
      docker run --rm -u "$UIDGID" \
        -v "$BASEDIR":/work -w /work \
        -e TMPDIR="/work/tmp/$RUN" \
        "$IMAGE" \
        prefetch "$RUN" --output-directory "/work/sra_files"
    else
      echo "[1/4] PREFETCH: local .sra file already exists"
    fi

    # Check whether prefetch generated the expected .sra file.
    if [[ ! -s "$sra_path" ]]; then
      echo "[ERROR] .sra file not found after prefetch: $sra_path" >&2
      exit 2
    fi

    # Step 2: validate the .sra file.
    echo "[2/4] VDB-VALIDATE: validating .sra file..."
    docker run --rm -u "$UIDGID" \
      -v "$BASEDIR":/work -w /work \
      "$IMAGE" \
      vdb-validate "/work/sra_files/$RUN/$RUN.sra"

    # Step 3: convert .sra to FASTQ.
    #
    # The sra_files directory is mounted inside the container at the default
    # SRA repository path to avoid VDB/path-related issues.
    echo "[3/4] FASTQ-DUMP: converting .sra to FASTQ..."
    docker run --rm -u "$UIDGID" \
      -v "$BASEDIR":/work -w /work \
      -v "$SRA_DIR":/root/ncbi/public/sra \
      -e TMPDIR="/work/tmp/$RUN" \
      "$IMAGE" \
      fastq-dump "$RUN" \
        --split-files --readids --clip --skip-technical \
        -O "/work/fastq"

    # Check whether FASTQ files were created successfully.
    if [[ -s "$out1" && -s "$out2" ]]; then
      echo "[OK] FASTQ files generated: $out1 and $out2"
    else
      echo "[ERROR] Expected FASTQ files were not found or are empty for run: $RUN" >&2
      echo "Filtered OUTDIR contents:" >&2
      ls -lah "$OUTDIR" | grep -E "${RUN}_" || true
      exit 3
    fi

    echo "[4/4] DONE"
    echo "END: $(date)"
  } >> "$run_log" 2>&1
}

# Simple parallel execution control without requiring GNU parallel.
running=0

while read -r RUN; do
  run_one "$RUN" &
  running=$((running + 1))

  # When the maximum number of parallel jobs is reached, wait for one job
  # to finish if `wait -n` is available. Otherwise, wait for all jobs.
  if [[ "$running" -ge "$MAX_JOBS" ]]; then
    if wait -n 2>/dev/null; then
      running=$((running - 1))
    else
      wait
      running=0
    fi
  fi
done < "$RUNS_CLEAN"

# Wait for remaining background jobs.
wait

# Compress downloaded FASTQ files.
echo "Compressing downloaded FASTQ files..."

shopt -s nullglob
fastq_files=("$OUTDIR"/*.fastq)

if [[ ${#fastq_files[@]} -gt 0 ]]; then
  gzip -f "${fastq_files[@]}"
  echo "[OK] FASTQ compression completed."
else
  echo "[INFO] No uncompressed FASTQ files found in: $OUTDIR"
fi

echo
echo "Pipeline finished."
echo "FASTQ directory: $OUTDIR"
echo "Log directory:   $LOGDIR"

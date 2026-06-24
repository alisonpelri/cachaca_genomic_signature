#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Run FastQC and MultiQC on FASTQ files.
#
# This script can be used for both raw and trimmed reads.
# It scans an input directory for FASTQ files, runs FastQC
# using Docker, and summarizes all FastQC reports with MultiQC.
#
# Accepted input extensions:
#   .fastq.gz
#   .fq.gz
#   .fastq
#   .fq
#
# Usage:
#   bash scripts/03_fastqc.sh <input_dir> <output_dir>
#
# Examples:
#   # Raw reads
#   bash scripts/03_fastqc.sh fastq results/fastqc_raw
#
#   # Trimmed reads
#   bash scripts/03_fastqc.sh results/trimmed results/fastqc_trimmed
#
# Requirements:
#   - Docker
#   - FASTQ files in the input directory
#
# Docker images:
#   - staphb/fastqc
#   - staphb/multiqc
# ============================================================

INPUT_DIR="${1:-}"
OUTDIR="${2:-}"

FASTQC_IMAGE="staphb/fastqc"
MULTIQC_IMAGE="staphb/multiqc"

# ------------------------------------------------------------
# Usage message
# ------------------------------------------------------------

usage() {
  cat << EOF
Usage:
  bash scripts/03_fastqc.sh <input_dir> <output_dir>

Examples:
  bash scripts/03_fastqc.sh fastq results/fastqc_raw
  bash scripts/03_fastqc.sh results/trimmed results/fastqc_trimmed
EOF
}

# ------------------------------------------------------------
# Argument and dependency checks
# ------------------------------------------------------------

if [[ -z "$INPUT_DIR" || -z "$OUTDIR" ]]; then
  usage
  exit 1
fi

if [[ ! -d "$INPUT_DIR" ]]; then
  echo "[ERROR] Input directory not found: $INPUT_DIR" >&2
  exit 1
fi

command -v docker >/dev/null 2>&1 || {
  echo "[ERROR] Docker is not installed or not available in PATH." >&2
  exit 1
}

docker info >/dev/null 2>&1 || {
  echo "[ERROR] Docker is not running or is not available for the current user." >&2
  exit 1
}

mkdir -p "$OUTDIR"

UIDGID="$(id -u):$(id -g)"

# ------------------------------------------------------------
# Find FASTQ files
# ------------------------------------------------------------

shopt -s nullglob

FASTQ_FILES=(
  "$INPUT_DIR"/*.fastq.gz
  "$INPUT_DIR"/*.fq.gz
  "$INPUT_DIR"/*.fastq
  "$INPUT_DIR"/*.fq
)

if [[ ${#FASTQ_FILES[@]} -eq 0 ]]; then
  echo "[ERROR] No FASTQ files found in: $INPUT_DIR" >&2
  echo "Accepted extensions: .fastq.gz, .fq.gz, .fastq, .fq" >&2
  exit 1
fi

# ------------------------------------------------------------
# Run FastQC
# ------------------------------------------------------------

echo "Input directory:  $INPUT_DIR"
echo "Output directory: $OUTDIR"
echo "FASTQ files:      ${#FASTQ_FILES[@]}"
echo "FastQC image:     $FASTQC_IMAGE"
echo "MultiQC image:    $MULTIQC_IMAGE"
echo

for fastq in "${FASTQ_FILES[@]}"; do
  file_name="$(basename "$fastq")"
  echo "[FastQC] Processing: $file_name"

  docker run --rm \
    --user "$UIDGID" \
    -v "$(pwd)":/data \
    -w /data \
    "$FASTQC_IMAGE" \
    fastqc \
      -o "$OUTDIR" \
      "$fastq"
done

# ------------------------------------------------------------
# Run MultiQC
# ------------------------------------------------------------

echo
echo "[MultiQC] Summarizing FastQC reports..."

docker run --rm \
  --user "$UIDGID" \
  -v "$(pwd)":/data \
  -w /data \
  "$MULTIQC_IMAGE" \
  multiqc "$OUTDIR" \
    -o "$OUTDIR"

echo
echo "FastQC/MultiQC analysis completed."
echo "Reports saved in: $OUTDIR"

#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Map trimmed reads to a reference genome using Bowtie2.
#
# This script:
#   1) Builds a Bowtie2 index for the reference genome, if needed
#   2) Maps paired-end trimmed reads
#   3) Maps single-end trimmed reads
#
# Expected paired-end input pattern:
#   SAMPLE_trimmed_R1.fq.gz
#   SAMPLE_trimmed_R2.fq.gz
#
# Expected single-end input pattern:
#   SAMPLE_trimmed.fq.gz
#   or
#   SAMPLE_trimmed.fq
#
# Usage:
#   bash scripts/04_map_reads.sh <reference_fasta> <input_dir> <output_dir> [threads] [index_prefix]
#
# Examples:
#   bash scripts/04_map_reads.sh data/reference/GCF_000146045.2_R64_genomic.fna results/trimmed results/mapping 10 R64
#
# Default:
#   threads      = 10
#   index_prefix = reference basename without extension
#
# Requirements:
#   - Docker
#   - Trimmed FASTQ files
#   - Reference genome in FASTA format
#
# Docker image:
#   - staphb/bowtie2
# ============================================================

REFERENCE="${1:-}"
INPUT_DIR="${2:-}"
OUTDIR="${3:-}"
THREADS="${4:-10}"
INDEX_PREFIX="${5:-}"

BOWTIE2_IMAGE="staphb/bowtie2"
UIDGID="$(id -u):$(id -g)"

# ------------------------------------------------------------
# Usage message
# ------------------------------------------------------------

usage() {
  cat << EOF
Usage:
  bash scripts/04_map_reads.sh <reference_fasta> <input_dir> <output_dir> [threads] [index_prefix]

Example:
  bash scripts/04_map_reads.sh data/reference/GCF_000146045.2_R64_genomic.fna results/trimmed results/mapping 10 R64
EOF
}

# ------------------------------------------------------------
# Argument and dependency checks
# ------------------------------------------------------------

if [[ -z "$REFERENCE" || -z "$INPUT_DIR" || -z "$OUTDIR" ]]; then
  usage
  exit 1
fi

if [[ ! -f "$REFERENCE" ]]; then
  echo "[ERROR] Reference FASTA not found: $REFERENCE" >&2
  exit 1
fi

if [[ ! -d "$INPUT_DIR" ]]; then
  echo "[ERROR] Input directory not found: $INPUT_DIR" >&2
  exit 1
fi

if ! [[ "$THREADS" =~ ^[0-9]+$ ]]; then
  echo "[ERROR] THREADS must be a positive integer." >&2
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

mkdir -p "$OUTDIR" "$OUTDIR/index" "$OUTDIR/logs" "$OUTDIR/sam"

REFERENCE_ABS="$(readlink -f "$REFERENCE")"
INPUT_DIR_ABS="$(readlink -f "$INPUT_DIR")"
OUTDIR_ABS="$(readlink -f "$OUTDIR")"

REFERENCE_BASENAME="$(basename "$REFERENCE_ABS")"

if [[ -z "$INDEX_PREFIX" ]]; then
  INDEX_PREFIX="${REFERENCE_BASENAME%.*}"
fi

INDEX_PATH="$OUTDIR_ABS/index/$INDEX_PREFIX"

echo "Reference FASTA: $REFERENCE_ABS"
echo "Input directory: $INPUT_DIR_ABS"
echo "Output directory: $OUTDIR_ABS"
echo "Threads:         $THREADS"
echo "Index prefix:    $INDEX_PREFIX"
echo "Docker image:    $BOWTIE2_IMAGE"
echo

# ------------------------------------------------------------
# Build Bowtie2 index, if missing
# ------------------------------------------------------------

if compgen -G "${INDEX_PATH}*.bt2" > /dev/null || compgen -G "${INDEX_PATH}*.bt2l" > /dev/null; then
  echo "[Index] Existing Bowtie2 index found. Skipping index construction."
else
  echo "[Index] Building Bowtie2 index..."

  docker run --rm \
    --user "$UIDGID" \
    -v "$(dirname "$REFERENCE_ABS")":/reference \
    -v "$OUTDIR_ABS":/output \
    "$BOWTIE2_IMAGE" \
    bowtie2-build \
      "/reference/$REFERENCE_BASENAME" \
      "/output/index/$INDEX_PREFIX"

  echo "[Index] Bowtie2 index created: $OUTDIR/index/$INDEX_PREFIX"
fi

echo

# ------------------------------------------------------------
# Map paired-end reads
# ------------------------------------------------------------

shopt -s nullglob

PE_R1_FILES=("$INPUT_DIR_ABS"/*_trimmed_R1.fq.gz)

if [[ ${#PE_R1_FILES[@]} -gt 0 ]]; then
  echo "Paired-end samples detected: ${#PE_R1_FILES[@]}"
  echo

  for r1 in "${PE_R1_FILES[@]}"; do
    sample="$(basename "$r1" _trimmed_R1.fq.gz)"
    r2="$INPUT_DIR_ABS/${sample}_trimmed_R2.fq.gz"

    if [[ ! -f "$r2" ]]; then
      echo "[WARNING] Missing R2 file for sample '$sample'. Skipping." >&2
      continue
    fi

    sam_out="$OUTDIR_ABS/sam/${sample}_PE_bt.sam"
    log_out="$OUTDIR_ABS/logs/${sample}_PE_bowtie2.log"

    if [[ -s "$sam_out" ]]; then
      echo "[SKIP] SAM file already exists for paired-end sample: $sample"
      continue
    fi

    echo "[Bowtie2 PE] Mapping sample: $sample"

    docker run --rm \
      --user "$UIDGID" \
      -v "$INPUT_DIR_ABS":/reads \
      -v "$OUTDIR_ABS":/output \
      "$BOWTIE2_IMAGE" \
      bowtie2 \
        -x "/output/index/$INDEX_PREFIX" \
        -1 "/reads/${sample}_trimmed_R1.fq.gz" \
        -2 "/reads/${sample}_trimmed_R2.fq.gz" \
        -S "/output/sam/${sample}_PE_bt.sam" \
        -p "$THREADS" \
        2> "$log_out"
  done
else
  echo "No paired-end files matching *_trimmed_R1.fq.gz were found."
fi

echo

# ------------------------------------------------------------
# Map single-end reads
# ------------------------------------------------------------

SE_FILES=(
  "$INPUT_DIR_ABS"/*_trimmed.fq.gz
  "$INPUT_DIR_ABS"/*_trimmed.fq
)

if [[ ${#SE_FILES[@]} -gt 0 ]]; then
  echo "Single-end samples detected: ${#SE_FILES[@]}"
  echo

  for se_read in "${SE_FILES[@]}"; do
    base="$(basename "$se_read")"

    if [[ "$base" == *_trimmed_R1.fq.gz || "$base" == *_trimmed_R2.fq.gz ]]; then
      continue
    fi

    sample="$base"
    sample="${sample%_trimmed.fq.gz}"
    sample="${sample%_trimmed.fq}"

    sam_out="$OUTDIR_ABS/sam/${sample}_SE_bt.sam"
    log_out="$OUTDIR_ABS/logs/${sample}_SE_bowtie2.log"

    if [[ -s "$sam_out" ]]; then
      echo "[SKIP] SAM file already exists for single-end sample: $sample"
      continue
    fi

    echo "[Bowtie2 SE] Mapping sample: $sample"

    docker run --rm \
      --user "$UIDGID" \
      -v "$INPUT_DIR_ABS":/reads \
      -v "$OUTDIR_ABS":/output \
      "$BOWTIE2_IMAGE" \
      bowtie2 \
        -x "/output/index/$INDEX_PREFIX" \
        -U "/reads/$base" \
        -S "/output/sam/${sample}_SE_bt.sam" \
        -p "$THREADS" \
        2> "$log_out"
  done
else
  echo "No single-end files matching *_trimmed.fq.gz or *_trimmed.fq were found."
fi

echo
echo "Bowtie2 mapping completed."
echo "SAM files: $OUTDIR/sam"
echo "Logs:      $OUTDIR/logs"
echo "Index:     $OUTDIR/index"

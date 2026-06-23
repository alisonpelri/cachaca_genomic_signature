```bash
#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Trim Illumina FASTQ reads using Trimmomatic.
#
# This script supports paired-end and single-end reads.
#
# Paired-end files are detected using the pattern:
#   SAMPLE_1.fastq.gz
#   SAMPLE_2.fastq.gz
#
# Single-end files are detected among remaining .fastq.gz files
# that do not match paired-end R1/R2 patterns.
#
# Usage:
#   bash scripts/03_trim_reads.sh [input_dir] [output_dir] [threads]
#
# Examples:
#   bash scripts/03_trim_reads.sh
#   bash scripts/03_trim_reads.sh fastq results/trimmed 10
#
# Default:
#   input_dir  = fastq
#   output_dir = results/trimmed
#   threads    = 10
# ============================================================

INPUT_DIR="${1:-fastq}"
OUTDIR="${2:-results/trimmed}"
THREADS="${3:-10}"

TRIMMOMATIC_IMAGE="staphb/trimmomatic"
UIDGID="$(id -u):$(id -g)"

if [[ ! -d "$INPUT_DIR" ]]; then
  echo "[ERROR] Input directory not found: $INPUT_DIR" >&2
  exit 1
fi

if ! [[ "$THREADS" =~ ^[0-9]+$ ]]; then
  echo "[ERROR] THREADS must be a positive integer." >&2
  exit 1
fi

mkdir -p "$OUTDIR"

echo "Input directory:  $INPUT_DIR"
echo "Output directory: $OUTDIR"
echo "Threads:          $THREADS"
echo

shopt -s nullglob

# ------------------------------------------------------------
# Paired-end trimming
# ------------------------------------------------------------

PE_R1_FILES=("$INPUT_DIR"/*_1.fastq.gz)
PROCESSED_PE_SAMPLES=()

if [[ ${#PE_R1_FILES[@]} -gt 0 ]]; then
  echo "Paired-end files detected: ${#PE_R1_FILES[@]}"
  echo

  for r1 in "${PE_R1_FILES[@]}"; do
    sample="$(basename "$r1" _1.fastq.gz)"
    r2="$INPUT_DIR/${sample}_2.fastq.gz"

    if [[ ! -f "$r2" ]]; then
      echo "[WARNING] Missing R2 file for sample '$sample'. Skipping paired-end trimming." >&2
      continue
    fi

    echo "[Trimmomatic PE] Processing sample: $sample"

    docker run --rm \
      --user "$UIDGID" \
      -v "$(pwd)":/data \
      -w /data \
      "$TRIMMOMATIC_IMAGE" \
      trimmomatic PE \
        -threads "$THREADS" \
        "$INPUT_DIR/${sample}_1.fastq.gz" \
        "$INPUT_DIR/${sample}_2.fastq.gz" \
        "$OUTDIR/${sample}_trimmed_R1.fq.gz" \
        "$OUTDIR/${sample}_unpaired_R1.fq.gz" \
        "$OUTDIR/${sample}_trimmed_R2.fq.gz" \
        "$OUTDIR/${sample}_unpaired_R2.fq.gz" \
        "ILLUMINACLIP:TruSeq3-PE:2:30:10" \
        SLIDINGWINDOW:4:30 \
        MINLEN:50

    PROCESSED_PE_SAMPLES+=("$sample")
  done
else
  echo "No paired-end files matching *_1.fastq.gz were found."
fi

echo

# ------------------------------------------------------------
# Single-end trimming
# ------------------------------------------------------------

ALL_FASTQ_FILES=("$INPUT_DIR"/*.fastq.gz)

if [[ ${#ALL_FASTQ_FILES[@]} -eq 0 ]]; then
  echo "[ERROR] No .fastq.gz files found in: $INPUT_DIR" >&2
  exit 1
fi

echo "Checking for single-end files..."
echo

for fastq in "${ALL_FASTQ_FILES[@]}"; do
  base="$(basename "$fastq")"

  # Skip paired-end reads already handled above.
  if [[ "$base" == *_1.fastq.gz || "$base" == *_2.fastq.gz ]]; then
    continue
  fi

  sample="${base%.fastq.gz}"

  echo "[Trimmomatic SE] Processing sample: $sample"

  docker run --rm \
    --user "$UIDGID" \
    -v "$(pwd)":/data \
    -w /data \
    "$TRIMMOMATIC_IMAGE" \
    trimmomatic SE \
      -threads "$THREADS" \
      "$INPUT_DIR/${sample}.fastq.gz" \
      "$OUTDIR/${sample}_trimmed.fq.gz" \
      SLIDINGWINDOW:4:30 \
      MINLEN:50
done

echo
echo "Read trimming completed."
echo "Trimmed reads saved in: $OUTDIR"
```


#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

###############################################################################
# Call variants from BAM files using bcftools.
#
# This script:
#   1) Indexes the reference genome with samtools faidx, if needed
#   2) Runs bcftools mpileup
#   3) Runs bcftools call
#   4) Generates one VCF file per BAM sample
#
# Expected BAM filename pattern:
#   SAMPLE_PE_mapped_sorted_dedup_q20_bt.bam
#
# Output:
#   <output_dir>/SAMPLE.vcf
#
# Usage:
#   bash scripts/07_variant_calling.sh <reference_fasta> <bam_dir> <output_dir> [threads]
#
# Example:
#   bash scripts/07_variant_calling.sh \
#     data/reference/GCF_000146045.2_R64_genomic.fna \
#     results/mapping/sam/06_q20_bam \
#     results/variants \
#     15
#
# Requirements:
#   - Docker
#   - Indexed or indexable reference FASTA
#   - BAM files generated from read mapping
#
# Docker images:
#   - staphb/samtools
#   - staphb/bcftools
###############################################################################

REFERENCE="${1:-}"
BAM_DIR="${2:-}"
OUTDIR="${3:-}"
THREADS="${4:-15}"

SAMTOOLS_IMAGE="staphb/samtools"
BCFTOOLS_IMAGE="staphb/bcftools"

show_help() {
    cat << EOF
Usage:
  bash scripts/07_variant_calling.sh <reference_fasta> <bam_dir> <output_dir> [threads]

Arguments:
  <reference_fasta>   Reference genome in FASTA format
  <bam_dir>           Directory containing BAM files
  <output_dir>        Directory where VCF files will be saved
  [threads]           Number of threads to use

Expected BAM filename pattern:
  SAMPLE_PE_mapped_sorted_dedup_q20_bt.bam

Example:
  bash scripts/07_variant_calling.sh \\
    data/reference/GCF_000146045.2_R64_genomic.fna \\
    results/mapping/sam/06_q20_bam \\
    results/variants \\
    15
EOF
}

###############################################################################
# Validate arguments and dependencies
###############################################################################

if [[ -z "$REFERENCE" || -z "$BAM_DIR" || -z "$OUTDIR" ]]; then
    show_help
    exit 1
fi

if [[ ! -f "$REFERENCE" ]]; then
    echo "[ERROR] Reference FASTA not found: $REFERENCE" >&2
    exit 1
fi

if [[ ! -d "$BAM_DIR" ]]; then
    echo "[ERROR] BAM directory not found: $BAM_DIR" >&2
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

mkdir -p "$OUTDIR"

REFERENCE_ABS="$(readlink -f "$REFERENCE")"
REFERENCE_DIR="$(dirname "$REFERENCE_ABS")"
REFERENCE_BASE="$(basename "$REFERENCE_ABS")"

BAM_DIR_ABS="$(readlink -f "$BAM_DIR")"
OUTDIR_ABS="$(readlink -f "$OUTDIR")"

UIDGID="$(id -u):$(id -g)"

###############################################################################
# Locate BAM files
###############################################################################

BAM_FILES=("${BAM_DIR_ABS}"/*_PE_mapped_sorted_dedup_q20_bt.bam)

if [[ ${#BAM_FILES[@]} -eq 0 ]]; then
    echo "[ERROR] No BAM files matching '*_PE_mapped_sorted_dedup_q20_bt.bam' found in: $BAM_DIR_ABS" >&2
    exit 1
fi

echo "Reference FASTA: $REFERENCE_ABS"
echo "BAM directory:    $BAM_DIR_ABS"
echo "Output directory: $OUTDIR_ABS"
echo "Threads:          $THREADS"
echo "BAM files:        ${#BAM_FILES[@]}"
echo

###############################################################################
# Index reference genome
###############################################################################

if [[ -f "${REFERENCE_ABS}.fai" ]]; then
    echo "[Reference] Existing FASTA index found: ${REFERENCE_ABS}.fai"
else
    echo "[Reference] Indexing reference genome with samtools faidx..."

    docker run --rm \
        --user "$UIDGID" \
        -v "$REFERENCE_DIR":/reference \
        "$SAMTOOLS_IMAGE" \
        samtools faidx "/reference/$REFERENCE_BASE"

    echo "[Reference] FASTA index created: ${REFERENCE_ABS}.fai"
fi

echo

###############################################################################
# Variant calling
###############################################################################

for bam in "${BAM_FILES[@]}"; do
    bam_base="$(basename "$bam")"
    sample="${bam_base%_PE_mapped_sorted_dedup_q20_bt.bam}"
    out_vcf="$OUTDIR_ABS/${sample}.vcf"

    if [[ -s "$out_vcf" ]]; then
        echo "[SKIP] VCF already exists for sample: $sample"
        continue
    fi

    echo "[bcftools] Calling variants for sample: $sample"

    docker run --rm \
        --user "$UIDGID" \
        -v "$REFERENCE_DIR":/reference \
        -v "$BAM_DIR_ABS":/bam \
        -v "$OUTDIR_ABS":/output \
        "$BCFTOOLS_IMAGE" \
        bash -c "
            set -euo pipefail

            bcftools mpileup \
                --threads ${THREADS} \
                -f /reference/${REFERENCE_BASE} \
                /bam/${bam_base} | \
            bcftools call \
                --threads ${THREADS} \
                -mv -Ov \
                -o /output/${sample}.vcf
        "

    if [[ -s "$out_vcf" ]]; then
        echo "[OK] VCF generated: $out_vcf"
    else
        echo "[ERROR] VCF was not generated or is empty for sample: $sample" >&2
        exit 1
    fi

    echo
done

echo "Variant calling completed."
echo "VCF files saved in: $OUTDIR_ABS"

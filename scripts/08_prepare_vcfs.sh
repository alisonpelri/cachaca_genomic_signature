```bash id="uucj63"
#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Prepare VCF files using bcftools, bgzip, and tabix with Docker.
#
# This script has two modes:
#
#   1) rename
#      - Renames VCF files based on a mapping table
#      - Replaces the sample name inside each single-sample VCF
#      - Compresses output with bgzip
#      - Indexes output with tabix
#
#   2) set-ids
#      - Creates variant IDs based on CHROM, POS, REF, and ALT
#      - Compresses output as VCF.GZ
#      - Indexes output with tabix
#
# Docker images:
#   - staphb/bcftools:latest
#   - biocontainers/tabix:v1.9-11-deb_cv1
#
# Requirements:
#   - Docker
#   - Input VCF files
#
# Usage:
#   bash scripts/09_prepare_vcfs.sh rename -m mapping.tsv [-o output_dir] [-t threads] file1.vcf [file2.vcf ...]
#   bash scripts/09_prepare_vcfs.sh set-ids -i input.vcf.gz -o output.vcf.gz
#
# Examples:
#   bash scripts/09_prepare_vcfs.sh rename \
#     -m data/metadata/rename_map.tsv \
#     -o results/variants/renamed_vcfs \
#     -t 8 \
#     results/variants/*.vcf
#
#   bash scripts/09_prepare_vcfs.sh set-ids \
#     -i results/filtered_vcf/Sace_filtered_snvs.vcf.gz \
#     -o results/filtered_vcf/Sace_pos_id.vcf.gz
###############################################################################

BCFTOOLS_IMG="staphb/bcftools:latest"
TABIX_IMG="biocontainers/tabix:v1.9-11-deb_cv1"
UID_GID="$(id -u):$(id -g)"

###############################################################################
# Help messages
###############################################################################

main_usage() {
    cat << EOF
Usage:
  bash scripts/09_prepare_vcfs.sh <mode> [options]

Modes:
  rename    Rename single-sample VCF files using a mapping table
  set-ids   Set variant IDs using CHROM, POS, REF, and ALT

Run one of the following for mode-specific help:
  bash scripts/09_prepare_vcfs.sh rename --help
  bash scripts/09_prepare_vcfs.sh set-ids --help
EOF
}

rename_usage() {
    cat << EOF
Usage:
  bash scripts/09_prepare_vcfs.sh rename -m mapping.tsv [-o output_dir] [-t threads] file1.vcf [file2.vcf file3.vcf.gz ...]

Description:
  - Renames each VCF file based on a mapping table
  - Replaces the sample name inside each VCF
  - Compresses output with bgzip
  - Indexes output with tabix using Docker

Required:
  -m   Mapping TSV file with two columns: old_id<TAB>new_id

Optional:
  -o   Output directory
       Default: renamed_vcfs

  -t   Number of threads for bgzip
       Default: 4

Examples:
  bash scripts/09_prepare_vcfs.sh rename -m rename_map.tsv ERR1308862.vcf ERR1352845.vcf
  bash scripts/09_prepare_vcfs.sh rename -m rename_map.tsv -o final_vcfs -t 8 *.vcf
EOF
}

set_ids_usage() {
    cat << EOF
Usage:
  bash scripts/09_prepare_vcfs.sh set-ids -i input.vcf.gz -o output.vcf.gz

Description:
  - Creates variant IDs based on CHROM, POS, REF, and ALT
  - Uses the following ID format:
      CHROM:POS_REF/ALT
  - Compresses output as VCF.GZ
  - Indexes output with tabix

Required:
  -i   Input VCF or VCF.GZ file
  -o   Output VCF.GZ file

Example:
  bash scripts/09_prepare_vcfs.sh set-ids \\
    -i results/filtered_vcf/Sace_filtered_snvs.vcf.gz \\
    -o results/filtered_vcf/Sace_pos_id.vcf.gz
EOF
}

###############################################################################
# Dependency checks
###############################################################################

check_docker() {
    command -v docker >/dev/null 2>&1 || {
        echo "[ERROR] Docker is not installed or not available in PATH." >&2
        exit 1
    }

    docker info >/dev/null 2>&1 || {
        echo "[ERROR] Docker is not running or is not available for the current user." >&2
        exit 1
    }
}

###############################################################################
# Shared helper functions
###############################################################################

lookup_new_id() {
    local old_id="$1"
    local mapfile_abs="$2"

    awk -F '\t' -v key="$old_id" '
        BEGIN { found=0 }
        NR==1 && ($1 ~ /^(Run|old_id|sample|sample_id)$/ || $2 ~ /^(Tree_id|new_id|new_sample|new_sample_id)$/) { next }
        $1 == key { print $2; found=1; exit }
        END { if (!found) exit 1 }
    ' "$mapfile_abs"
}

get_sample_name() {
    local file_abs="$1"
    local file_dir
    local file_base

    file_dir="$(dirname "$file_abs")"
    file_base="$(basename "$file_abs")"

    docker run --rm \
        --user "$UID_GID" \
        -v "$file_dir":/in \
        "$BCFTOOLS_IMG" \
        bcftools query -l "/in/$file_base"
}

###############################################################################
# Mode 1: rename VCF files
###############################################################################

run_rename_mode() {
    local MAPFILE=""
    local OUTDIR="renamed_vcfs"
    local THREADS=4

    while [[ $# -gt 0 ]]; do
        case "$1" in
            -m)
                MAPFILE="${2:-}"
                shift 2
                ;;
            -o)
                OUTDIR="${2:-}"
                shift 2
                ;;
            -t)
                THREADS="${2:-}"
                shift 2
                ;;
            -h|--help)
                rename_usage
                exit 0
                ;;
            --)
                shift
                break
                ;;
            -*)
                echo "[ERROR] Unknown option for rename mode: $1" >&2
                rename_usage
                exit 1
                ;;
            *)
                break
                ;;
        esac
    done

    if [[ -z "$MAPFILE" ]]; then
        echo "[ERROR] Mapping file is required." >&2
        rename_usage
        exit 1
    fi

    if [[ $# -lt 1 ]]; then
        echo "[ERROR] Provide at least one VCF file." >&2
        rename_usage
        exit 1
    fi

    if [[ ! -f "$MAPFILE" ]]; then
        echo "[ERROR] Mapping file not found: $MAPFILE" >&2
        exit 1
    fi

    if ! [[ "$THREADS" =~ ^[0-9]+$ ]]; then
        echo "[ERROR] THREADS must be a positive integer." >&2
        exit 1
    fi

    mkdir -p "$OUTDIR"

    local MAPFILE_ABS
    local OUTDIR_ABS

    MAPFILE_ABS="$(readlink -f "$MAPFILE")"
    OUTDIR_ABS="$(readlink -f "$OUTDIR")"

    echo "Mode:             rename"
    echo "Mapping file:     $MAPFILE_ABS"
    echo "Output directory: $OUTDIR_ABS"
    echo "Threads:          $THREADS"
    echo

    for input in "$@"; do
        if [[ ! -f "$input" ]]; then
            echo "[WARNING] File not found, skipping: $input" >&2
            continue
        fi

        local input_abs
        local input_dir
        local input_base
        local old_id
        local new_id
        local sample_names
        local sample_count
        local current_sample
        local tmpdir
        local sample_map
        local reheadered_vcf
        local final_vcf_gz

        input_abs="$(readlink -f "$input")"
        input_dir="$(dirname "$input_abs")"
        input_base="$(basename "$input_abs")"

        old_id="$input_base"
        old_id="${old_id%.vcf.gz}"
        old_id="${old_id%.vcf}"

        echo "Processing: $input_base"
        echo "  Old ID: $old_id"

        if ! new_id="$(lookup_new_id "$old_id" "$MAPFILE_ABS")"; then
            echo "  [WARNING] No mapping found for '$old_id'. Skipping." >&2
            echo
            continue
        fi

        echo "  New ID: $new_id"

        sample_names="$(get_sample_name "$input_abs" || true)"
        sample_count="$(printf '%s\n' "$sample_names" | sed '/^$/d' | wc -l)"

        if [[ "$sample_count" -eq 0 ]]; then
            echo "  [ERROR] No sample found inside VCF: $input_base" >&2
            echo
            continue
        elif [[ "$sample_count" -gt 1 ]]; then
            echo "  [ERROR] This script expects single-sample VCFs, but this file has more than one sample: $input_base" >&2
            echo
            continue
        fi

        current_sample="$(printf '%s\n' "$sample_names" | head -n 1)"
        echo "  Current sample in VCF: $current_sample"

        tmpdir="$(mktemp -d)"
        sample_map="$tmpdir/sample_map.txt"
        reheadered_vcf="$tmpdir/${new_id}.vcf"
        final_vcf_gz="$OUTDIR_ABS/${new_id}.vcf.gz"

        printf '%s\t%s\n' "$current_sample" "$new_id" > "$sample_map"

        echo "  Reheadering sample name..."

        docker run --rm \
            --user "$UID_GID" \
            -v "$input_dir":/in \
            -v "$tmpdir":/work \
            "$BCFTOOLS_IMG" \
            sh -c "bcftools reheader -s /work/sample_map.txt /in/$input_base -o /work/$(basename "$reheadered_vcf")"

        echo "  Compressing with bgzip..."

        docker run --rm \
            --user "$UID_GID" \
            -v "$tmpdir":/work \
            -v "$OUTDIR_ABS":/out \
            "$TABIX_IMG" \
            sh -c "bgzip -@ $THREADS -c /work/$(basename "$reheadered_vcf") > /out/$(basename "$final_vcf_gz")"

        echo "  Indexing with tabix..."

        docker run --rm \
            --user "$UID_GID" \
            -v "$OUTDIR_ABS":/out \
            "$TABIX_IMG" \
            tabix -f -p vcf "/out/$(basename "$final_vcf_gz")"

        echo "  Done:"
        echo "    $final_vcf_gz"
        echo "    ${final_vcf_gz}.tbi"

        rm -rf "$tmpdir"
        echo
    done

    echo "Rename mode completed."
}

###############################################################################
# Mode 2: set variant IDs
###############################################################################

run_set_ids_mode() {
    local INPUT_VCF=""
    local OUTPUT_VCF=""

    while [[ $# -gt 0 ]]; do
        case "$1" in
            -i)
                INPUT_VCF="${2:-}"
                shift 2
                ;;
            -o)
                OUTPUT_VCF="${2:-}"
                shift 2
                ;;
            -h|--help)
                set_ids_usage
                exit 0
                ;;
            *)
                echo "[ERROR] Unknown option for set-ids mode: $1" >&2
                set_ids_usage
                exit 1
                ;;
        esac
    done

    if [[ -z "$INPUT_VCF" || -z "$OUTPUT_VCF" ]]; then
        echo "[ERROR] Input and output VCF files are required." >&2
        set_ids_usage
        exit 1
    fi

    if [[ ! -f "$INPUT_VCF" ]]; then
        echo "[ERROR] Input VCF not found: $INPUT_VCF" >&2
        exit 1
    fi

    if [[ "$OUTPUT_VCF" != *.vcf.gz ]]; then
        echo "[ERROR] Output file must end with .vcf.gz" >&2
        exit 1
    fi

    local INPUT_ABS
    local INPUT_DIR
    local INPUT_BASE
    local OUTPUT_DIR
    local OUTPUT_BASE
    local OUTPUT_ABS

    INPUT_ABS="$(readlink -f "$INPUT_VCF")"
    INPUT_DIR="$(dirname "$INPUT_ABS")"
    INPUT_BASE="$(basename "$INPUT_ABS")"

    mkdir -p "$(dirname "$OUTPUT_VCF")"

    OUTPUT_DIR="$(readlink -f "$(dirname "$OUTPUT_VCF")")"
    OUTPUT_BASE="$(basename "$OUTPUT_VCF")"
    OUTPUT_ABS="$OUTPUT_DIR/$OUTPUT_BASE"

    echo "Mode:       set-ids"
    echo "Input VCF:  $INPUT_ABS"
    echo "Output VCF: $OUTPUT_ABS"
    echo

    echo "Creating variant IDs using CHROM:POS_REF/ALT..."

    docker run --rm \
        --user "$UID_GID" \
        -v "$INPUT_DIR":/in \
        -v "$OUTPUT_DIR":/out \
        "$BCFTOOLS_IMG" \
        bcftools annotate \
            --set-id '%CHROM:%POS\_%REF/%ALT' \
            "/in/$INPUT_BASE" \
            -Oz \
            -o "/out/$OUTPUT_BASE"

    echo "Indexing output VCF..."

    docker run --rm \
        --user "$UID_GID" \
        -v "$OUTPUT_DIR":/out \
        "$TABIX_IMG" \
        tabix -f -p vcf "/out/$OUTPUT_BASE"

    echo
    echo "Set-ids mode completed."
    echo "Output:"
    echo "  $OUTPUT_ABS"
    echo "  ${OUTPUT_ABS}.tbi"
}

###############################################################################
# Main
###############################################################################

check_docker

if [[ $# -lt 1 ]]; then
    main_usage
    exit 1
fi

MODE="$1"
shift

case "$MODE" in
    rename)
        run_rename_mode "$@"
        ;;
    set-ids)
        run_set_ids_mode "$@"
        ;;
    -h|--help)
        main_usage
        exit 0
        ;;
    *)
        echo "[ERROR] Unknown mode: $MODE" >&2
        echo
        main_usage
        exit 1
        ;;
esac
```

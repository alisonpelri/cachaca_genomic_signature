```bash
#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

###############################################################################
# Generate summary statistics for individual VCF files.
#
# This script calculates per-sample statistics for single-sample VCF files.
#
# For each VCF, it reports:
#   - total number of variants
#   - number of SNPs
#   - number of indels
#   - Ts/Tv ratio
#   - mean QUAL
#   - mean INFO/DP, if available
#   - number of homozygous reference genotypes
#   - number of heterozygous genotypes
#   - number of homozygous alternate genotypes
#   - number of missing genotypes
#
# Input:
#   Directory containing .vcf or .vcf.gz files.
#
# Output:
#   <outdir>/gz/          compressed and indexed VCF files
#   <outdir>/stats/       bcftools stats reports
#   <outdir>/query/       genotype query outputs
#   <outdir>/summary/     final summary table
#
# Usage:
#   bash scripts/09_vcf_individual_stats.sh <vcf_dir> <output_dir> [threads]
#
# Example:
#   bash scripts/09_vcf_individual_stats.sh \
#     results/variants/renamed_vcfs \
#     results/vcf_individual_stats \
#     4
#
# Requirements:
#   - Docker
#   - gawk
#
# Docker image:
#   - staphb/bcftools
###############################################################################

VCF_DIR="${1:-}"
OUTDIR="${2:-}"
THREADS="${3:-4}"

BCFTOOLS_IMAGE="staphb/bcftools"
UID_GID="$(id -u):$(id -g)"

show_help() {
    cat << EOF
Usage:
  bash scripts/09_vcf_individual_stats.sh <vcf_dir> <output_dir> [threads]

Arguments:
  <vcf_dir>      Directory containing individual .vcf or .vcf.gz files
  <output_dir>   Directory where statistics will be saved
  [threads]      Number of threads for bgzip/tabix when needed

Example:
  bash scripts/09_vcf_individual_stats.sh \\
    results/variants/renamed_vcfs \\
    results/vcf_individual_stats \\
    4
EOF
}

###############################################################################
# Validate arguments and dependencies
###############################################################################

if [[ -z "$VCF_DIR" || -z "$OUTDIR" ]]; then
    show_help
    exit 1
fi

if [[ ! -d "$VCF_DIR" ]]; then
    echo "[ERROR] VCF directory not found: $VCF_DIR" >&2
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

command -v gawk >/dev/null 2>&1 || {
    echo "[ERROR] gawk is not installed or not available in PATH." >&2
    exit 1
}

###############################################################################
# Prepare paths and output directories
###############################################################################

mkdir -p \
    "${OUTDIR}/gz" \
    "${OUTDIR}/stats" \
    "${OUTDIR}/query" \
    "${OUTDIR}/summary"

VCF_DIR_ABS="$(readlink -f "$VCF_DIR")"
OUTDIR_ABS="$(readlink -f "$OUTDIR")"

SUMMARY_FILE="${OUTDIR_ABS}/summary/vcf_summary.tsv"

vcf_plain=("${VCF_DIR_ABS}"/*.vcf)
vcf_gz=("${VCF_DIR_ABS}"/*.vcf.gz)

if [[ ${#vcf_plain[@]} -eq 0 && ${#vcf_gz[@]} -eq 0 ]]; then
    echo "[ERROR] No .vcf or .vcf.gz files found in: $VCF_DIR_ABS" >&2
    exit 1
fi

echo -e "sample\ttotal_variants\tsnps\tindels\ttstv\tmean_qual\tmean_info_dp\tn_hom_ref\tn_het\tn_hom_alt\tn_missing" \
    > "$SUMMARY_FILE"

echo "VCF directory:    $VCF_DIR_ABS"
echo "Output directory: $OUTDIR_ABS"
echo "Threads:          $THREADS"
echo "Plain VCF files:  ${#vcf_plain[@]}"
echo "Gzipped VCFs:     ${#vcf_gz[@]}"
echo

###############################################################################
# Docker helper functions
###############################################################################

run_bcftools() {
    docker run --rm \
        --user "$UID_GID" \
        -v "$OUTDIR_ABS":/out \
        "$BCFTOOLS_IMAGE" \
        bcftools "$@"
}

compress_vcf() {
    local input_vcf="$1"
    local output_gz="$2"

    local input_dir
    local input_base
    local output_dir
    local output_base

    input_dir="$(dirname "$input_vcf")"
    input_base="$(basename "$input_vcf")"
    output_dir="$(dirname "$output_gz")"
    output_base="$(basename "$output_gz")"

    docker run --rm \
        --user "$UID_GID" \
        -v "$input_dir":/in \
        -v "$output_dir":/out \
        "$BCFTOOLS_IMAGE" \
        bash -c "bgzip -@ ${THREADS} -c /in/${input_base} > /out/${output_base}"
}

index_vcf() {
    local vcf_gz="$1"

    local vcf_dir
    local vcf_base

    vcf_dir="$(dirname "$vcf_gz")"
    vcf_base="$(basename "$vcf_gz")"

    docker run --rm \
        --user "$UID_GID" \
        -v "$vcf_dir":/vcf \
        "$BCFTOOLS_IMAGE" \
        tabix -f -p vcf "/vcf/${vcf_base}"
}

bcftools_on_file() {
    local vcf_gz="$1"
    shift

    local vcf_dir
    local vcf_base

    vcf_dir="$(dirname "$vcf_gz")"
    vcf_base="$(basename "$vcf_gz")"

    docker run --rm \
        --user "$UID_GID" \
        -v "$vcf_dir":/vcf \
        "$BCFTOOLS_IMAGE" \
        bcftools "$@" "/vcf/${vcf_base}"
}

###############################################################################
# Process one VCF file
###############################################################################

process_vcf() {
    local input_vcf="$1"
    local sample="$2"

    local gzvcf="${OUTDIR_ABS}/gz/${sample}.vcf.gz"
    local stats_file="${OUTDIR_ABS}/stats/${sample}.vchk"
    local gt_file="${OUTDIR_ABS}/query/${sample}.gt.txt"

    local total
    local snps
    local indels
    local mean_qual
    local mean_info_dp
    local tstv
    local n_hom_ref
    local n_het
    local n_hom_alt
    local n_missing

    echo "[VCF stats] Processing sample: ${sample}"

    # Ensure that the VCF is compressed and indexed.
    if [[ "$input_vcf" == *.vcf ]]; then
        echo "  - Compressing VCF with bgzip"
        compress_vcf "$input_vcf" "$gzvcf"

        echo "  - Indexing VCF with tabix"
        index_vcf "$gzvcf"
    else
        echo "  - Copying compressed VCF to output directory"
        cp "$input_vcf" "$gzvcf"

        if [[ -f "${input_vcf}.tbi" ]]; then
            cp "${input_vcf}.tbi" "${gzvcf}.tbi"
        else
            echo "  - Indexing VCF with tabix"
            index_vcf "$gzvcf"
        fi
    fi

    # General VCF statistics.
    echo "  - Running bcftools stats"
    bcftools_on_file "$gzvcf" stats > "$stats_file"

    # Count variants, SNPs, and indels.
    total="$(bcftools_on_file "$gzvcf" view -H | wc -l)"
    snps="$(bcftools_on_file "$gzvcf" view -H -v snps | wc -l)"
    indels="$(bcftools_on_file "$gzvcf" view -H -v indels | wc -l)"

    # Mean QUAL.
    mean_qual="$(
        bcftools_on_file "$gzvcf" query -f '%QUAL\n' | \
        gawk 'NF>0 && $1!="." {sum+=$1; n++} END{if(n>0) printf "%.4f", sum/n; else print "NA"}'
    )"

    # Mean INFO/DP, if available.
    mean_info_dp="$(
        bcftools_on_file "$gzvcf" query -f '%INFO/DP\n' 2>/dev/null | \
        gawk 'NF>0 && $1!="." {sum+=$1; n++} END{if(n>0) printf "%.4f", sum/n; else print "NA"}'
    )"

    # Extract Ts/Tv from bcftools stats.
    tstv="$(
        gawk '
            $1=="TSTV" && $3=="0" { val=$5 }
            END {
                if (val=="") print "NA"
                else print val
            }
        ' "$stats_file"
    )"

    # Genotype summary.
    #
    # Note:
    # In yeast, the biological interpretation of heterozygosity depends on
    # ploidy and study design. Here, genotype counts are descriptive QC metrics.
    bcftools_on_file "$gzvcf" query -f '[%GT\n]' > "$gt_file" 2>/dev/null || true

    if [[ -s "$gt_file" ]]; then
        read -r n_hom_ref n_het n_hom_alt n_missing < <(
            gawk '
            {
                gt=$1

                if (gt=="0/0" || gt=="0|0") {
                    homref++
                }
                else if (gt=="0/1" || gt=="1/0" || gt=="0|1" || gt=="1|0") {
                    het++
                }
                else if (gt=="1/1" || gt=="1|1") {
                    homalt++
                }
                else if (gt=="./." || gt==".|." || gt==".") {
                    miss++
                }
            }
            END {
                print homref+0, het+0, homalt+0, miss+0
            }
            ' "$gt_file"
        )
    else
        n_hom_ref="NA"
        n_het="NA"
        n_hom_alt="NA"
        n_missing="NA"
    fi

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$sample" "$total" "$snps" "$indels" "$tstv" "$mean_qual" "$mean_info_dp" \
        "$n_hom_ref" "$n_het" "$n_hom_alt" "$n_missing" \
        >> "$SUMMARY_FILE"

    echo "  - Done"
    echo
}

###############################################################################
# Run
###############################################################################

declare -A seen_samples=()

for vcf in "${vcf_plain[@]}"; do
    sample="$(basename "$vcf" .vcf)"

    if [[ -n "${seen_samples[$sample]:-}" ]]; then
        echo "[WARNING] Duplicate sample basename detected, skipping: $sample" >&2
        continue
    fi

    seen_samples["$sample"]=1
    process_vcf "$vcf" "$sample"
done

for vcf in "${vcf_gz[@]}"; do
    sample="$(basename "$vcf" .vcf.gz)"

    if [[ -n "${seen_samples[$sample]:-}" ]]; then
        echo "[WARNING] Duplicate sample basename detected, skipping: $sample" >&2
        continue
    fi

    seen_samples["$sample"]=1
    process_vcf "$vcf" "$sample"
done

echo "Individual VCF statistics completed."
echo "Summary table: $SUMMARY_FILE"
```

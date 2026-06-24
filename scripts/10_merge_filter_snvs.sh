```bash id="d8nnd6"
#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

###############################################################################
# Merge and filter individual VCF files into a final SNV dataset.
#
# This script:
#   1) Checks and indexes input .vcf.gz files
#   2) Merges multiple individual VCF files with bcftools merge
#   3) Keeps biallelic SNPs and indels
#   4) Filters variants by QUAL >= 20
#   5) Requires minimum total ALT support >= 2 reads using FORMAT/AD
#   6) Filters variants by depth:
#        4 <= INFO/DP <= mean(INFO/DP) + 3 * sd(INFO/DP)
#   7) Applies proximity filters:
#        --SnpGap 5
#        --IndelGap 10
#   8) Keeps only SNVs in the final output VCF
#
# Important:
#   By default, this script uses bcftools merge -0, which treats missing
#   genotypes as reference genotypes. Use --no-missing-to-ref if you prefer to keep missing genotypes.
#
# Input:
#   Directory containing indexed or indexable .vcf.gz files.
#
# Output:
#   <output_dir>/<prefix>_merged_all_samples.vcf.gz
#   <output_dir>/<prefix>_step1_qual20_biallelic_allvars.vcf.gz
#   <output_dir>/<prefix>_step2_ad2_biallelic_allvars.vcf.gz
#   <output_dir>/<prefix>_step3_dp_filtered_allvars.vcf.gz
#   <output_dir>/<prefix>_step4_gap_filtered_allvars.vcf.gz
#   <output_dir>/<prefix>_filtered_snvs.vcf.gz
#   <output_dir>/<prefix>_merge_filter_summary.txt
#
# Usage:
#   bash scripts/10_merge_filter_snvs.sh -i <vcf_dir> -o <output_dir> [options]
#
# Required arguments:
#   -i, --vcf-dir       Directory containing input .vcf.gz files
#   -o, --outdir        Output directory
#
# Optional arguments:
#   -t, --threads       Number of threads
#                       Default: 8
#
#   -p, --prefix        Prefix for output files
#                       Default: Sace
#
#   --missing-to-ref    Use bcftools merge -0
#                       Default behavior
#
#   --no-missing-to-ref Do not use bcftools merge -0
#
#   -h, --help          Show this help message and exit
#
# Examples:
#   bash scripts/10_merge_filter_snvs.sh \
#     -i results/variants/renamed_vcfs \
#     -o results/filtered_vcf \
#     -t 8 \
#     -p Sace
#
#   bash scripts/10_merge_filter_snvs.sh \
#     -i results/variants/renamed_vcfs \
#     -o results/filtered_vcf \
#     --no-missing-to-ref
#
# Requirements:
#   - Docker
#   - awk
#
# Docker image:
#   - staphb/bcftools
#
# Filtering rationale:
#   The variant filtering strategy implemented here follows the approach used by
#   Rosse et al. (2017), who adapted filtering criteria from Choi et al. (2013).
#
#   In this pipeline, the implemented filters include:
#     - QUAL >= 20
#     - biallelic SNPs and indels
#     - minimum total ALT support >= 2 reads using FORMAT/AD
#     - depth filter: 4 <= INFO/DP <= mean(INFO/DP) + 3 * sd(INFO/DP)
#     - proximity filters: --SnpGap 5 and --IndelGap 10
#     - final retention of SNVs only
#
# References:
#   Rosse, I.C., Assis, J.G., Oliveira, F.S. et al. Whole genome sequencing of
#   Guzerá cattle reveals genetic variants in candidate genes for production,
#   disease resistance, and heat tolerance. Mamm Genome 28, 66–80 (2017).
#   https://doi.org/10.1007/s00335-016-9670-7
#
#   Choi, J-W., Liao, X., Park, S., Jeon, H-J., Chung, W-H., Stothard, P.,
#   Park, Y-S., Lee, J-K., Lee, K-T., Kim, S-H. (2013). Massively parallel
#   sequencing of Chikso (Korean brindle cattle) to discover genome-wide SNPs
#   and InDels. Molecular Cells 36, 203–211.
###############################################################################

VCF_DIR=""
OUTDIR=""
THREADS=8
PREFIX="Sace"
MISSING_TO_REF=1

BCFTOOLS_IMAGE="staphb/bcftools"
UID_GID="$(id -u):$(id -g)"

show_help() {
    sed -n '2,92p' "$0" | sed 's/^# \{0,1\}//'
}

###############################################################################
# Parse arguments
###############################################################################

if [[ $# -eq 0 ]]; then
    show_help
    exit 1
fi

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--vcf-dir)
            VCF_DIR="${2:-}"
            shift 2
            ;;
        -o|--outdir)
            OUTDIR="${2:-}"
            shift 2
            ;;
        -t|--threads)
            THREADS="${2:-}"
            shift 2
            ;;
        -p|--prefix)
            PREFIX="${2:-}"
            shift 2
            ;;
        --missing-to-ref)
            MISSING_TO_REF=1
            shift
            ;;
        --no-missing-to-ref)
            MISSING_TO_REF=0
            shift
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        *)
            echo "[ERROR] Unknown argument: $1" >&2
            echo
            show_help
            exit 1
            ;;
    esac
done

###############################################################################
# Validate arguments and dependencies
###############################################################################

if [[ -z "$VCF_DIR" || -z "$OUTDIR" ]]; then
    echo "[ERROR] Missing required arguments." >&2
    echo
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

command -v awk >/dev/null 2>&1 || {
    echo "[ERROR] awk is not installed or not available in PATH." >&2
    exit 1
}

mkdir -p "$OUTDIR"

VCF_DIR_ABS="$(readlink -f "$VCF_DIR")"
OUTDIR_ABS="$(readlink -f "$OUTDIR")"

TMPDIR_ABS="$(mktemp -d "${OUTDIR_ABS}/tmp.merge_filter.XXXXXX")"
trap 'rm -rf "${TMPDIR_ABS}"' EXIT

###############################################################################
# Docker helper
###############################################################################

bcftools_docker() {
    docker run --rm -i \
        --user "$UID_GID" \
        -v "$VCF_DIR_ABS":/vcfs \
        -v "$OUTDIR_ABS":/out \
        -v "$TMPDIR_ABS":/tmpwork \
        "$BCFTOOLS_IMAGE" \
        bcftools "$@"
}

###############################################################################
# Locate input VCF files
###############################################################################

VCFS=("${VCF_DIR_ABS}"/*.vcf.gz)

if [[ ${#VCFS[@]} -eq 0 ]]; then
    echo "[ERROR] No .vcf.gz files found in: $VCF_DIR_ABS" >&2
    exit 1
fi

VCF_LIST_HOST="${TMPDIR_ABS}/vcf_list.txt"
VCF_LIST_CONTAINER="/tmpwork/vcf_list.txt"

: > "$VCF_LIST_HOST"

for vcf in "${VCFS[@]}"; do
    printf "/vcfs/%s\n" "$(basename "$vcf")" >> "$VCF_LIST_HOST"
done

###############################################################################
# Output files
###############################################################################

MERGED_VCF="/out/${PREFIX}_merged_all_samples.vcf.gz"
STEP1="/out/${PREFIX}_step1_qual20_biallelic_allvars.vcf.gz"
STEP2="/out/${PREFIX}_step2_ad2_biallelic_allvars.vcf.gz"
STEP3="/out/${PREFIX}_step3_dp_filtered_allvars.vcf.gz"
STEP4="/out/${PREFIX}_step4_gap_filtered_allvars.vcf.gz"
FINAL_SNV="/out/${PREFIX}_filtered_snvs.vcf.gz"

MERGED_VCF_HOST="${OUTDIR_ABS}/${PREFIX}_merged_all_samples.vcf.gz"
STEP1_HOST="${OUTDIR_ABS}/${PREFIX}_step1_qual20_biallelic_allvars.vcf.gz"
STEP2_HOST="${OUTDIR_ABS}/${PREFIX}_step2_ad2_biallelic_allvars.vcf.gz"
STEP3_HOST="${OUTDIR_ABS}/${PREFIX}_step3_dp_filtered_allvars.vcf.gz"
STEP4_HOST="${OUTDIR_ABS}/${PREFIX}_step4_gap_filtered_allvars.vcf.gz"
FINAL_SNV_HOST="${OUTDIR_ABS}/${PREFIX}_filtered_snvs.vcf.gz"
SUMMARY_HOST="${OUTDIR_ABS}/${PREFIX}_merge_filter_summary.txt"

###############################################################################
# Report configuration
###############################################################################

echo "VCF directory:       $VCF_DIR_ABS"
echo "Output directory:    $OUTDIR_ABS"
echo "Input VCF files:     ${#VCFS[@]}"
echo "Threads:             $THREADS"
echo "Output prefix:       $PREFIX"

if [[ "$MISSING_TO_REF" -eq 1 ]]; then
    echo "Merge missing mode:  missing genotypes set to reference using bcftools merge -0"
else
    echo "Merge missing mode:  missing genotypes kept as missing"
fi

echo

###############################################################################
# Step 0: Check and index input VCF files
###############################################################################

echo "==> [0/7] Checking input VCF indexes..."

for vcf in "${VCFS[@]}"; do
    vcf_base="$(basename "$vcf")"

    if [[ ! -f "${vcf}.tbi" && ! -f "${vcf}.csi" ]]; then
        echo "    Indexing: $vcf_base"
        bcftools_docker index -f -t "/vcfs/${vcf_base}"
    fi
done

echo

###############################################################################
# Step 1: Merge VCF files
###############################################################################

echo "==> [1/7] Merging individual VCF files..."

MERGE_ARGS=(
    merge
    --threads "$THREADS"
    -l "$VCF_LIST_CONTAINER"
    -Oz
    -o "$MERGED_VCF"
)

if [[ "$MISSING_TO_REF" -eq 1 ]]; then
    MERGE_ARGS+=("-0")
fi

bcftools_docker "${MERGE_ARGS[@]}"
bcftools_docker index -f -t "$MERGED_VCF"

echo "    Merged VCF: $MERGED_VCF_HOST"
echo

###############################################################################
# Step 2: Keep biallelic SNPs/indels and filter by QUAL >= 20
###############################################################################

echo "==> [2/7] Keeping biallelic SNPs/indels and filtering QUAL >= 20..."

bcftools_docker view \
    --threads "$THREADS" \
    -m2 -M2 \
    -v snps,indels \
    "$MERGED_VCF" \
| bcftools_docker filter \
    --threads "$THREADS" \
    -i 'QUAL>=20' \
    -Oz \
    -o "$STEP1" \
    -

bcftools_docker index -f -t "$STEP1"

echo "    Step 1 VCF: $STEP1_HOST"
echo

###############################################################################
# Step 3: Filter by minimum total ALT support using FORMAT/AD
###############################################################################

echo "==> [3/7] Filtering sites with total ALT support >= 2 reads using FORMAT/AD..."

AD_QUERY_HOST="${TMPDIR_ABS}/ad_query.tsv"
KEEP_POS_HOST="${TMPDIR_ABS}/keep_by_ad.tsv"
KEEP_POS_CONTAINER="/tmpwork/keep_by_ad.tsv"

if ! bcftools_docker query -f '%CHROM\t%POS[\t%AD]\n' "$STEP1" > "$AD_QUERY_HOST"; then
    echo "[ERROR] Could not query FORMAT/AD from VCF." >&2
    echo "        Make sure the input VCF contains FORMAT/AD." >&2
    exit 1
fi

awk -F'\t' '
BEGIN { OFS="\t" }
{
    alt_sum = 0

    for (i = 3; i <= NF; i++) {
        if ($i == "." || $i == "./." || $i == "") {
            continue
        }

        n = split($i, a, ",")

        if (n >= 2 && a[2] != "." && a[2] != "") {
            alt_sum += a[2]
        }
    }

    if (alt_sum >= 2) {
        print $1, $2
    }
}
' "$AD_QUERY_HOST" > "$KEEP_POS_HOST"

if [[ ! -s "$KEEP_POS_HOST" ]]; then
    echo "[ERROR] No sites passed the FORMAT/AD filter." >&2
    echo "        Check whether FORMAT/AD exists and contains allele-depth values." >&2
    exit 1
fi

bcftools_docker view \
    --threads "$THREADS" \
    -T "$KEEP_POS_CONTAINER" \
    -Oz \
    -o "$STEP2" \
    "$STEP1"

bcftools_docker index -f -t "$STEP2"

echo "    Step 2 VCF: $STEP2_HOST"
echo

###############################################################################
# Step 4: Calculate depth threshold and filter by INFO/DP
###############################################################################

echo "==> [4/7] Calculating INFO/DP threshold..."

DP_VALUES_HOST="${TMPDIR_ABS}/dp_values.txt"
DP_STATS_HOST="${TMPDIR_ABS}/dp_stats.txt"

if ! bcftools_docker query -f '%INFO/DP\n' "$STEP2" > "$DP_VALUES_HOST"; then
    echo "[ERROR] Could not query INFO/DP from VCF." >&2
    echo "        Make sure the input VCF contains INFO/DP." >&2
    exit 1
fi

awk '
$1 != "." && $1 != "" {
    x = $1 + 0
    n++
    sum += x
    sumsq += x * x
}
END {
    if (n == 0) {
        print "ERROR"
        exit 1
    }

    mean = sum / n
    var = (sumsq / n) - (mean * mean)

    if (var < 0) {
        var = 0
    }

    sd = sqrt(var)
    maxdp = int(mean + 3 * sd)

    if (maxdp < 4) {
        maxdp = 4
    }

    printf("%.6f\t%.6f\t%d\n", mean, sd, maxdp)
}
' "$DP_VALUES_HOST" > "$DP_STATS_HOST"

if grep -q "ERROR" "$DP_STATS_HOST"; then
    echo "[ERROR] INFO/DP is missing or empty." >&2
    exit 1
fi

MEAN_DP="$(cut -f1 "$DP_STATS_HOST")"
SD_DP="$(cut -f2 "$DP_STATS_HOST")"
MAX_DP="$(cut -f3 "$DP_STATS_HOST")"
MIN_DP=4

echo "    mean(INFO/DP) = $MEAN_DP"
echo "    sd(INFO/DP)   = $SD_DP"
echo "    min DP         = $MIN_DP"
echo "    max DP         = $MAX_DP"

echo "==> [5/7] Filtering by depth..."

bcftools_docker filter \
    --threads "$THREADS" \
    -i "INFO/DP>=${MIN_DP} && INFO/DP<=${MAX_DP}" \
    -Oz \
    -o "$STEP3" \
    "$STEP2"

bcftools_docker index -f -t "$STEP3"

echo "    Step 3 VCF: $STEP3_HOST"
echo

###############################################################################
# Step 5: Apply proximity filters
###############################################################################

echo "==> [6/7] Applying proximity filters: --SnpGap 5 and --IndelGap 10..."

bcftools_docker filter \
    --threads "$THREADS" \
    --SnpGap 5 \
    --IndelGap 10 \
    -Oz \
    -o "$STEP4" \
    "$STEP3"

bcftools_docker index -f -t "$STEP4"

echo "    Step 4 VCF: $STEP4_HOST"
echo

###############################################################################
# Step 6: Keep only SNVs in the final VCF
###############################################################################

echo "==> [7/7] Keeping only SNVs in final VCF..."

bcftools_docker view \
    --threads "$THREADS" \
    -v snps \
    -Oz \
    -o "$FINAL_SNV" \
    "$STEP4"

bcftools_docker index -f -t "$FINAL_SNV"

FINAL_SNVS="$(bcftools_docker view -H "$FINAL_SNV" | wc -l | awk '{print $1}')"

echo "    Final SNVs: $FINAL_SNVS"
echo "    Final VCF:  $FINAL_SNV_HOST"
echo

###############################################################################
# Final summary
###############################################################################

cat > "$SUMMARY_HOST" << EOF
input_vcf_dir	$VCF_DIR_ABS
output_dir	$OUTDIR_ABS
prefix	$PREFIX
threads	$THREADS
n_input_vcfs	${#VCFS[@]}
missing_to_ref	$MISSING_TO_REF
mean_info_dp	$MEAN_DP
sd_info_dp	$SD_DP
min_dp	$MIN_DP
max_dp	$MAX_DP
final_snvs	$FINAL_SNVS
merged_vcf	$MERGED_VCF_HOST
step1_qual20_biallelic	$STEP1_HOST
step2_ad2	$STEP2_HOST
step3_dp_filtered	$STEP3_HOST
step4_gap_filtered	$STEP4_HOST
final_filtered_snvs	$FINAL_SNV_HOST
EOF

echo "Merge and filtering completed."
echo "Merged VCF:        $MERGED_VCF_HOST"
echo "Final filtered VCF: $FINAL_SNV_HOST"
echo "Summary:           $SUMMARY_HOST"
```

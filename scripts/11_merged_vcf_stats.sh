```bash
#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Generate summary statistics for a merged multi-sample VCF file.
#
# This script calculates dataset-level statistics for a merged VCF, especially
# useful after variant filtering.
#
# It reports:
#   - number of samples
#   - total number of variants
#   - number of SNPs
#   - number of indels
#   - number of biallelic SNPs
#   - Ts/Tv ratio
#   - mean site missingness
#   - median site missingness
#   - mean alternate allele frequency
#   - mean minor allele frequency
#   - number of sites with MAF < 0.01
#   - number of sites with MAF < 0.05
#   - number of sites with missingness > 10%
#   - number of sites with missingness > 20%
#
# Input:
#   One merged VCF file in .vcf.gz format.
#
# Output:
#   <outdir>/merged.vchk
#   <outdir>/merged.with_tags.vcf.gz
#   <outdir>/site_stats.tsv
#   <outdir>/merged_relevant_stats.tsv
#
# Usage:
#   bash scripts/11_merged_vcf_stats.sh <merged_vcf.gz> <output_dir>
#
# Example:
#   bash scripts/11_merged_vcf_stats.sh \
#     results/filtered_vcf/Sace_filtered_snvs.vcf.gz \
#     results/merged_vcf_stats
#
# Notes:
#   - AF and MAF summaries are most straightforward for biallelic sites.
#   - This script is intended for merged multi-sample VCFs, not individual VCFs.
#
# Requirements:
#   - Docker
#   - gawk
#
# Docker image:
#   - staphb/bcftools
###############################################################################

VCF="${1:-}"
OUTDIR="${2:-}"

BCFTOOLS_IMAGE="staphb/bcftools"
UID_GID="$(id -u):$(id -g)"

show_help() {
    cat << EOF
Usage:
  bash scripts/11_merged_vcf_stats.sh <merged_vcf.gz> <output_dir>

Arguments:
  <merged_vcf.gz>   Merged multi-sample VCF file
  <output_dir>      Directory where statistics will be saved

Example:
  bash scripts/11_merged_vcf_stats.sh \\
    results/filtered_vcf/Sace_filtered_snvs.vcf.gz \\
    results/merged_vcf_stats
EOF
}

###############################################################################
# Validate arguments and dependencies
###############################################################################

if [[ -z "$VCF" || -z "$OUTDIR" ]]; then
    show_help
    exit 1
fi

if [[ ! -f "$VCF" ]]; then
    echo "[ERROR] VCF file not found: $VCF" >&2
    exit 1
fi

if [[ "$VCF" != *.vcf.gz ]]; then
    echo "[ERROR] Input VCF must be compressed and end with .vcf.gz" >&2
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

mkdir -p "$OUTDIR"

VCF_ABS="$(readlink -f "$VCF")"
VCF_DIR="$(dirname "$VCF_ABS")"
VCF_BASE="$(basename "$VCF_ABS")"

OUTDIR_ABS="$(readlink -f "$OUTDIR")"

###############################################################################
# Docker helper
###############################################################################

bcftools_docker() {
    docker run --rm -i \
        --user "$UID_GID" \
        -v "$VCF_DIR":/vcf \
        -v "$OUTDIR_ABS":/out \
        "$BCFTOOLS_IMAGE" \
        bcftools "$@"
}

###############################################################################
# Output files
###############################################################################

STATS_FILE="/out/merged.vchk"
TAGS_VCF="/out/merged.with_tags.vcf.gz"
SITE_STATS_HOST="${OUTDIR_ABS}/site_stats.tsv"
SUMMARY_FILE="${OUTDIR_ABS}/merged_relevant_stats.tsv"

###############################################################################
# Start
###############################################################################

echo "Merged VCF:       $VCF_ABS"
echo "Output directory: $OUTDIR_ABS"
echo "Docker image:     $BCFTOOLS_IMAGE"
echo

###############################################################################
# Index VCF if needed
###############################################################################

if [[ ! -f "${VCF_ABS}.tbi" && ! -f "${VCF_ABS}.csi" ]]; then
    echo "[Index] Indexing merged VCF..."
    bcftools_docker index -f -t "/vcf/$VCF_BASE"
else
    echo "[Index] Existing VCF index found."
fi

echo

###############################################################################
# 1) General VCF statistics
###############################################################################

echo "[1/5] Running bcftools stats..."

bcftools_docker stats "/vcf/$VCF_BASE" > "${OUTDIR_ABS}/merged.vchk"

###############################################################################
# 2) Fill site-level tags
###############################################################################

echo "[2/5] Filling site-level tags..."

# +fill-tags calculates:
#   AN        = number of observed alleles
#   AC        = alternate allele count
#   AF        = alternate allele frequency
#   NS        = number of samples with data
#   F_MISSING = fraction of missing genotypes
bcftools_docker +fill-tags "/vcf/$VCF_BASE" \
    -Oz \
    -o "$TAGS_VCF" \
    -- -t AN,AC,AF,NS,F_MISSING

bcftools_docker index -f -t "$TAGS_VCF"

###############################################################################
# 3) Extract site-level statistics
###############################################################################

echo "[3/5] Extracting site-level statistics..."

bcftools_docker query \
    -f '%CHROM\t%POS\t%AN\t%AC\t%AF\t%NS\t%INFO/F_MISSING\n' \
    "$TAGS_VCF" \
    > "$SITE_STATS_HOST"

###############################################################################
# 4) Global VCF counts
###############################################################################

echo "[4/5] Calculating global VCF counts..."

total="$(bcftools_docker view -H "/vcf/$VCF_BASE" | wc -l | awk '{print $1}')"
snps="$(bcftools_docker view -H -v snps "/vcf/$VCF_BASE" | wc -l | awk '{print $1}')"
indels="$(bcftools_docker view -H -v indels "/vcf/$VCF_BASE" | wc -l | awk '{print $1}')"
biallelic_snps="$(bcftools_docker view -H -v snps -m2 -M2 "/vcf/$VCF_BASE" | wc -l | awk '{print $1}')"
n_samples="$(bcftools_docker query -l "/vcf/$VCF_BASE" | wc -l | awk '{print $1}')"

###############################################################################
# 5) Ts/Tv and AF/MAF/missingness summary
###############################################################################

echo "[5/5] Summarizing Ts/Tv, AF, MAF, and missingness..."

tstv="$(
    gawk '
        $1=="TSTV" && $3=="0" { val=$5 }
        END {
            if (val=="") print "NA"
            else print val
        }
    ' "${OUTDIR_ABS}/merged.vchk"
)"

# AF  = alternate allele frequency
# MAF = minor allele frequency, calculated as min(AF, 1-AF)
#
# For multiallelic sites, bcftools may return comma-separated AF values.
# This script uses the first AF value, but it is primarily intended for
# biallelic filtered VCF files.
read -r mean_missing median_missing mean_af mean_maf maf_lt_001 maf_lt_005 miss_gt_10 miss_gt_20 < <(
    gawk '
    BEGIN {
        FS = "\t"
    }

    {
        af_raw = $5
        miss = $7

        if (af_raw != "." && af_raw != "" && miss != "." && miss != "") {
            split(af_raw, af_arr, ",")
            af = af_arr[1] + 0
            miss = miss + 0

            n++

            af_sum += af
            miss_sum += miss

            maf = af
            if (maf > 0.5) {
                maf = 1 - maf
            }

            maf_sum += maf
            miss_arr[n] = miss

            if (maf < 0.01) maf001++
            if (maf < 0.05) maf005++
            if (miss > 0.10) miss10++
            if (miss > 0.20) miss20++
        }
    }

    END {
        if (n > 0) {
            asort(miss_arr)

            if (n % 2 == 1) {
                median = miss_arr[(n + 1) / 2]
            } else {
                median = (miss_arr[n / 2] + miss_arr[(n / 2) + 1]) / 2
            }

            printf "%.6f %.6f %.6f %.6f %d %d %d %d\n",
                   miss_sum / n,
                   median,
                   af_sum / n,
                   maf_sum / n,
                   maf001 + 0,
                   maf005 + 0,
                   miss10 + 0,
                   miss20 + 0
        } else {
            print "NA NA NA NA NA NA NA NA"
        }
    }
    ' "$SITE_STATS_HOST"
)

###############################################################################
# Write final summary table
###############################################################################

echo -e "dataset\tn_samples\ttotal_variants\tsnps\tindels\tbiallelic_snps\ttstv\tmean_site_missingness\tmedian_site_missingness\tmean_af\tmean_maf\tsites_maf_lt_0.01\tsites_maf_lt_0.05\tsites_missing_gt_10pct\tsites_missing_gt_20pct" \
    > "$SUMMARY_FILE"

dataset_name="$(basename "$VCF_BASE" .vcf.gz)"

printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$dataset_name" \
    "$n_samples" \
    "$total" \
    "$snps" \
    "$indels" \
    "$biallelic_snps" \
    "$tstv" \
    "$mean_missing" \
    "$median_missing" \
    "$mean_af" \
    "$mean_maf" \
    "$maf_lt_001" \
    "$maf_lt_005" \
    "$miss_gt_10" \
    "$miss_gt_20" \
    >> "$SUMMARY_FILE"

echo
echo "Merged VCF statistics completed."
echo "Summary table: $SUMMARY_FILE"
echo "Site-level table: $SITE_STATS_HOST"
echo "Tagged VCF: ${OUTDIR_ABS}/merged.with_tags.vcf.gz"
```

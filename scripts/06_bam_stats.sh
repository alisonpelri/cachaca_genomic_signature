```bash id="qsaqa4"
#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

###############################################################################
# Generate summary statistics for BAM files using Samtools.
#
# This script calculates mapping, duplication, coverage, and depth statistics
# for all BAM files in a given directory.
#
# For each BAM file, it generates:
#   - samtools flagstat output
#   - samtools coverage output
#   - samtools idxstats output
#   - samtools depth output
#   - one combined summary table
#
# Output directories:
#   <outdir>/flagstat/
#   <outdir>/coverage/
#   <outdir>/idxstats/
#   <outdir>/depth/
#   <outdir>/summary/
#
# Usage:
#   bash scripts/06_bam_stats.sh <bam_dir> <output_dir> [threads]
#
# Examples:
#   bash scripts/06_bam_stats.sh results/mapping/sam/06_q20_bam results/bam_stats 8
#   bash scripts/06_bam_stats.sh results/mapping/sam/04_dedup_bam results/bam_stats_dedup 8
#
# Default:
#   threads = 8
#
# Requirements:
#   - samtools
#   - gawk
###############################################################################

BAM_DIR="${1:-}"
OUTDIR="${2:-}"
THREADS="${3:-8}"

show_help() {
    cat << EOF
Usage:
  bash scripts/06_bam_stats.sh <bam_dir> <output_dir> [threads]

Arguments:
  <bam_dir>      Directory containing BAM files
  <output_dir>   Directory where statistics will be saved
  [threads]      Number of threads for samtools index/flagstat

Examples:
  bash scripts/06_bam_stats.sh results/mapping/sam/06_q20_bam results/bam_stats 8
  bash scripts/06_bam_stats.sh results/mapping/sam/04_dedup_bam results/bam_stats_dedup 8
EOF
}

###############################################################################
# Validate arguments and dependencies
###############################################################################

if [[ -z "$BAM_DIR" || -z "$OUTDIR" ]]; then
    show_help
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

for exe in samtools gawk; do
    command -v "$exe" >/dev/null 2>&1 || {
        echo "[ERROR] Required executable not found in PATH: $exe" >&2
        exit 1
    }
done

###############################################################################
# Prepare output directories
###############################################################################

mkdir -p \
    "${OUTDIR}/flagstat" \
    "${OUTDIR}/coverage" \
    "${OUTDIR}/idxstats" \
    "${OUTDIR}/depth" \
    "${OUTDIR}/summary"

bam_files=("${BAM_DIR}"/*.bam)

if [[ ${#bam_files[@]} -eq 0 ]]; then
    echo "[ERROR] No .bam files found in: $BAM_DIR" >&2
    exit 1
fi

SUMMARY_FILE="${OUTDIR}/summary/bam_summary.tsv"

echo -e "sample\ttotal_reads\tmapped_reads\tmapping_rate_pct\tduplicate_reads\tduplicate_rate_pct\tcontigs_in_ref\tcontigs_with_reads\tcovered_bases\tgenome_bases\tbreadth_pct\tmean_depth\tsd_depth\tmedian_depth\tp10\tp25\tp50\tp75\tp90\tp95\tp99\tpct_ge_1x\tpct_ge_5x\tpct_ge_10x\tpct_ge_20x\tpct_ge_30x" \
    > "$SUMMARY_FILE"

echo "BAM directory:   $BAM_DIR"
echo "Output directory: $OUTDIR"
echo "Threads:         $THREADS"
echo "BAM files:       ${#bam_files[@]}"
echo

###############################################################################
# Process BAM files
###############################################################################

for bam in "${bam_files[@]}"; do
    sample="$(basename "$bam" .bam)"
    echo "[BAM stats] Processing sample: ${sample}"

    # Index BAM if needed.
    if [[ ! -f "${bam}.bai" ]]; then
        echo "  - Indexing BAM"
        samtools index -@ "${THREADS}" "$bam"
    fi

    # Generate primary statistics.
    echo "  - Running samtools flagstat"
    samtools flagstat -@ "${THREADS}" "$bam" > "${OUTDIR}/flagstat/${sample}.flagstat.txt"

    echo "  - Running samtools coverage"
    samtools coverage "$bam" > "${OUTDIR}/coverage/${sample}.coverage.tsv"

    echo "  - Running samtools idxstats"
    samtools idxstats "$bam" > "${OUTDIR}/idxstats/${sample}.idxstats.tsv"

    echo "  - Running samtools depth"
    samtools depth -aa "$bam" > "${OUTDIR}/depth/${sample}.depth.tsv"

    # Summarize flagstat + coverage + depth into one line.
    gawk -v sample="$sample" -v OFS="\t" '
        BEGIN {
            total_reads="NA"; mapped_reads="NA"; mapping_rate="NA";
            duplicate_reads="NA"; duplicate_rate="NA";

            contigs_in_ref=0; contigs_with_reads=0;
            covered_bases=0; genome_bases=0; breadth_pct="NA";

            mean_depth="NA"; sd_depth="NA"; median_depth="NA";
            p10="NA"; p25="NA"; p50="NA"; p75="NA"; p90="NA"; p95="NA"; p99="NA";

            pct1="NA"; pct5="NA"; pct10="NA"; pct20="NA"; pct30="NA";
        }

        FNR == NR {
            # Parse samtools flagstat output.
            if ($0 ~ /in total/) {
                total_reads = $1
            }
            else if ($0 ~ / mapped \(/ && $0 !~ /primary mapped/) {
                mapped_reads = $1
                if (match($0, /\(([0-9.]+)%/, m)) mapping_rate = m[1]
            }
            else if ($0 ~ /duplicates/) {
                duplicate_reads = $1
                if (total_reads != "NA" && total_reads > 0) {
                    duplicate_rate = sprintf("%.4f", 100 * duplicate_reads / total_reads)
                }
            }
            next
        }

        ARGIND == 2 {
            # Parse samtools coverage output.
            # Expected columns:
            # #rname startpos endpos numreads covbases coverage meandepth meanbaseq meanmapq
            if ($1 ~ /^#rname/ || NF < 7) next

            contigs_in_ref++

            if ($4 > 0) {
                contigs_with_reads++
            }

            covered_bases += $5
            genome_bases += ($3 - $2 + 1)
            next
        }

        ARGIND == 3 {
            # Parse samtools depth output.
            # Expected columns:
            # CHROM POS DEPTH
            d[++n] = $3
            sum += $3

            if ($3 >= 1)  c1++
            if ($3 >= 5)  c5++
            if ($3 >= 10) c10++
            if ($3 >= 20) c20++
            if ($3 >= 30) c30++

            next
        }

        END {
            if (genome_bases > 0) {
                breadth_pct = sprintf("%.4f", 100 * covered_bases / genome_bases)
            }

            if (n > 0) {
                mean_depth = sum / n
                asort(d)

                idx10 = int(n * 0.10); if (idx10 < 1) idx10 = 1
                idx25 = int(n * 0.25); if (idx25 < 1) idx25 = 1
                idx50 = int(n * 0.50); if (idx50 < 1) idx50 = 1
                idx75 = int(n * 0.75); if (idx75 < 1) idx75 = 1
                idx90 = int(n * 0.90); if (idx90 < 1) idx90 = 1
                idx95 = int(n * 0.95); if (idx95 < 1) idx95 = 1
                idx99 = int(n * 0.99); if (idx99 < 1) idx99 = 1

                p10 = d[idx10]
                p25 = d[idx25]
                p50 = d[idx50]
                p75 = d[idx75]
                p90 = d[idx90]
                p95 = d[idx95]
                p99 = d[idx99]
                median_depth = p50

                for (i = 1; i <= n; i++) {
                    ss += (d[i] - mean_depth)^2
                }

                sd_depth = sqrt(ss / n)

                pct1  = sprintf("%.4f", 100 * c1 / n)
                pct5  = sprintf("%.4f", 100 * c5 / n)
                pct10 = sprintf("%.4f", 100 * c10 / n)
                pct20 = sprintf("%.4f", 100 * c20 / n)
                pct30 = sprintf("%.4f", 100 * c30 / n)

                mean_depth = sprintf("%.4f", mean_depth)
                sd_depth   = sprintf("%.4f", sd_depth)
            }

            print sample, total_reads, mapped_reads, mapping_rate,
                  duplicate_reads, duplicate_rate,
                  contigs_in_ref, contigs_with_reads,
                  covered_bases, genome_bases, breadth_pct,
                  mean_depth, sd_depth, median_depth,
                  p10, p25, p50, p75, p90, p95, p99,
                  pct1, pct5, pct10, pct20, pct30
        }
    ' \
    "${OUTDIR}/flagstat/${sample}.flagstat.txt" \
    "${OUTDIR}/coverage/${sample}.coverage.tsv" \
    "${OUTDIR}/depth/${sample}.depth.tsv" \
    >> "$SUMMARY_FILE"

    echo "  - Done"
    echo
done

echo "BAM statistics completed."
echo "Summary table: $SUMMARY_FILE"
```

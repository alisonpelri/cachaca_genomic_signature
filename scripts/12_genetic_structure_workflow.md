# Genetic Structure Analysis Workflow

This document organizes the genetic structure analyses from a filtered multi-sample VCF file.

The workflow includes:

1. Directory setup
2. PLINK conversion
3. Basic quality control
4. LD pruning
5. PCA
6. Identity-by-state distance
7. Pairwise Fst
8. Allele-frequency delta
9. Combined Fst/delta ranking
10. ADMIXTURE analysis

> **Note:** This markdown file is intended to document a workflow. For full automation, these commands can later be split into bash/R scripts.

---

## Expected input files

Place the following files in `genetic_structure_analysis/data/`:

```text
genetic_structure_analysis/
└── data/
    ├── Sace_pos_id.vcf.gz
    ├── Sace_pos_id.vcf.gz.tbi
    └── metadata.tsv
```

The metadata file should contain at least two columns:

```text
sample	group
sample_01	Cachaca
sample_02	Wine
sample_03	Beer
```

---

## 0. Create working directories

```bash
mkdir -p genetic_structure_analysis/{data,01_plink_conversion,02_qc,03_pruned,04_pca,05_ibs,06_fst,07_admixture,08_delta,09_feature_integration,10_logs}
```

Expected structure:

```text
genetic_structure_analysis/
├── data/
│   ├── Sace_pos_id.vcf.gz
│   ├── Sace_pos_id.vcf.gz.tbi
│   └── metadata.tsv
├── 01_plink_conversion/
├── 02_qc/
├── 03_pruned/
├── 04_pca/
├── 05_ibs/
├── 06_fst/
├── 07_admixture/
├── 08_delta/
├── 09_feature_integration/
└── 10_logs/
```

Change into the analysis directory before running the next commands:

```bash
cd genetic_structure_analysis
```

---

## 1. Convert VCF to PLINK format

```bash
plink \
  --vcf data/Sace_pos_id.vcf.gz \
  --double-id \
  --nonfounders \
  --allow-no-sex \
  --allow-extra-chr \
  --recode \
  --make-bed \
  --out 01_plink_conversion/Sace_dataset
```

Expected main outputs:

```text
01_plink_conversion/Sace_dataset.bed
01_plink_conversion/Sace_dataset.bim
01_plink_conversion/Sace_dataset.fam
01_plink_conversion/Sace_dataset.ped
01_plink_conversion/Sace_dataset.map
```

---

## 2. Create PLINK cluster file

The cluster file links each sample to its biological group.

```bash
awk 'BEGIN{OFS="\t"} NR>1 {print $1, $1, $2}' \
  data/metadata.tsv \
  > 02_qc/groups.clst
```

Output:

```text
02_qc/groups.clst
```

---

## 3. PLINK quality control

Filters used:

- `--geno 0.05`: remove variants with missingness greater than 5%
- `--mind 0.05`: remove samples with missingness greater than 5%
- `--maf 0.05`: remove rare variants with MAF lower than 5%

```bash
plink \
  --bfile 01_plink_conversion/Sace_dataset \
  --allow-extra-chr \
  --geno 0.05 \
  --mind 0.05 \
  --maf 0.05 \
  --make-bed \
  --out 02_qc/Sace_qc
```

Expected main outputs:

```text
02_qc/Sace_qc.bed
02_qc/Sace_qc.bim
02_qc/Sace_qc.fam
```

### 3.1. Basic QC summaries

Missingness:

```bash
plink \
  --bfile 02_qc/Sace_qc \
  --allow-extra-chr \
  --missing \
  --out 02_qc/Sace_qc_missing
```

Allele frequency:

```bash
plink \
  --bfile 02_qc/Sace_qc \
  --allow-extra-chr \
  --freq \
  --out 02_qc/Sace_qc_freq
```

---

## 4. LD pruning

This step reduces redundancy among variants in linkage disequilibrium.

Parameters:

- Window size: 50 variants
- Step size: 10 variants
- LD threshold: `r² > 0.3`

```bash
plink \
  --bfile 02_qc/Sace_qc \
  --allow-extra-chr \
  --indep-pairwise 50 10 0.3 \
  --out 03_pruned/Sace_pruned
```

Create the pruned PLINK dataset:

```bash
plink \
  --bfile 02_qc/Sace_qc \
  --allow-extra-chr \
  --extract 03_pruned/Sace_pruned.prune.in \
  --make-bed \
  --out 03_pruned/Sace_qc_pruned
```

Expected main outputs:

```text
03_pruned/Sace_qc_pruned.bed
03_pruned/Sace_qc_pruned.bim
03_pruned/Sace_qc_pruned.fam
```

---

## 5. Principal component analysis

```bash
plink \
  --bfile 03_pruned/Sace_qc_pruned \
  --allow-extra-chr \
  --pca 20 \
  --out 04_pca/Sace_pca
```

Expected outputs:

```text
04_pca/Sace_pca.eigenval
04_pca/Sace_pca.eigenvec
```

---

## 6. Identity-by-state distance

Identity-by-state, or IBS, is calculated from allele-sharing distances among samples.

```bash
plink \
  --bfile 03_pruned/Sace_qc_pruned \
  --allow-extra-chr \
  --distance square 1-ibs \
  --out 05_ibs/Sace_ibs
```

Expected outputs:

```text
05_ibs/Sace_ibs.mdist
05_ibs/Sace_ibs.mdist.id
```

---

## 7. Pairwise Fst between groups

### 7.1. Export QC-filtered dataset as VCF

```bash
plink \
  --bfile 02_qc/Sace_qc \
  --allow-extra-chr \
  --recode vcf bgz \
  --out 06_fst/Sace_qc_fst
```

Index the exported VCF:

```bash
tabix -f -p vcf 06_fst/Sace_qc_fst.vcf.gz
```

### 7.2. Create group-specific sample files

```bash
awk '$2=="Beer"{print $1}' data/metadata.tsv > 06_fst/Beer.samples
awk '$2=="Cachaca"{print $1}' data/metadata.tsv > 06_fst/Cachaca.samples
awk '$2=="Wine"{print $1}' data/metadata.tsv > 06_fst/Wine.samples
awk '$2=="Bioethanol"{print $1}' data/metadata.tsv > 06_fst/Bioethanol.samples
awk '$2=="Spirits"{print $1}' data/metadata.tsv > 06_fst/Spirits.samples
```

Remove empty group files, if any:

```bash
find 06_fst -name "*.samples" -type f -empty -delete
```

### 7.3. Optional: fix PLINK-renamed sample IDs in the VCF header

PLINK may duplicate sample IDs in VCF headers, for example `sample_sample`. If this happens, create a rename map and reheader the VCF.

```bash
bcftools query -l 06_fst/Sace_qc_fst.vcf.gz | \
awk -F'_' '{
  n = NF / 2
  out = $1
  for (i = 2; i <= n; i++) out = out "_" $i
  print $0 "\t" out
}' > 06_fst/rename_samples.txt

bcftools reheader \
  -s 06_fst/rename_samples.txt \
  -o 06_fst/Sace_qc_fst_renamed.vcf.gz \
  06_fst/Sace_qc_fst.vcf.gz

tabix -f -p vcf 06_fst/Sace_qc_fst_renamed.vcf.gz
```

Use the renamed VCF in the next step if needed:

```bash
FST_VCF="Sace_qc_fst_renamed.vcf.gz"
```

Otherwise:

```bash
FST_VCF="Sace_qc_fst.vcf.gz"
```

### 7.4. Run pairwise Fst with vcftools

```bash
mkdir -p 06_fst/fast_results

cd 06_fst

FST_VCF="Sace_qc_fst_renamed.vcf.gz"

for i in $(ls *.samples | sed 's/\.samples//'); do
  for j in $(ls *.samples | sed 's/\.samples//'); do
    if [[ "$i" < "$j" ]]; then
      echo "Running Fst: ${i} vs ${j}"
      vcftools \
        --gzvcf "$FST_VCF" \
        --weir-fst-pop "${i}.samples" \
        --weir-fst-pop "${j}.samples" \
        --out "fast_results/${i}_vs_${j}"
    fi
  done
done

cd ..
```

---

## 8. Summarize pairwise Fst results

Save the following script as `09_feature_integration/fst_summary.R`.

```r
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(purrr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop(
    "Usage: Rscript fst_summary.R <variants.bim> <fst_dir> <out.tsv>\n",
    "Example: Rscript fst_summary.R 02_qc/Sace_qc.bim 06_fst/fast_results 09_feature_integration/fst_summary.tsv"
  )
}

bim_file <- args[1]
fst_dir  <- args[2]
out_file <- args[3]

bim <- fread(bim_file, header = FALSE)

if (ncol(bim) < 4) {
  stop("Invalid .bim file: expected at least 4 columns.")
}

setnames(bim, c("CHROM", "SNP", "CM", "POS", "A1", "A2")[seq_len(ncol(bim))])

variant_map <- bim %>%
  as_tibble() %>%
  transmute(
    CHROM = as.character(CHROM),
    POS   = as.integer(POS),
    SNP   = as.character(SNP)
  )

dup_coords <- variant_map %>%
  count(CHROM, POS) %>%
  filter(n > 1)

if (nrow(dup_coords) > 0) {
  warning(
    "Duplicated CHROM+POS coordinates found in the .bim file. ",
    "Joining Fst results by coordinate may be ambiguous."
  )
}

fst_files <- list.files(
  fst_dir,
  pattern = "\\.weir\\.fst$",
  full.names = TRUE
)

if (length(fst_files) == 0) {
  stop("No .weir.fst files found in: ", fst_dir)
}

all_fst <- map_dfr(fst_files, function(f) {
  x <- fread(f)

  expected_cols <- c("CHROM", "POS", "WEIR_AND_COCKERHAM_FST")
  if (!all(expected_cols %in% names(x))) {
    stop(
      "File ", f, " does not contain expected columns: ",
      paste(expected_cols, collapse = ", ")
    )
  }

  comparison_name <- basename(f) %>%
    str_remove("\\.weir\\.fst$")

  x %>%
    as_tibble() %>%
    transmute(
      CHROM = as.character(CHROM),
      POS   = as.integer(POS),
      fst   = as.numeric(WEIR_AND_COCKERHAM_FST)
    ) %>%
    left_join(variant_map, by = c("CHROM", "POS")) %>%
    mutate(comparison = comparison_name)
})

n_unmapped <- all_fst %>%
  filter(is.na(SNP)) %>%
  nrow()

if (n_unmapped > 0) {
  warning(
    n_unmapped, " Fst rows could not be mapped to SNP IDs by CHROM+POS.\n",
    "Check whether the .bim file and .weir.fst files came from the same variant set."
  )
}

all_fst_mapped <- all_fst %>%
  filter(!is.na(SNP))

fst_summary <- all_fst_mapped %>%
  group_by(SNP) %>%
  summarise(
    CHROM = first(CHROM),
    POS = first(POS),
    mean_fst = mean(fst, na.rm = TRUE),
    median_fst = median(fst, na.rm = TRUE),
    max_fst = max(fst, na.rm = TRUE),
    min_fst = min(fst, na.rm = TRUE),
    sd_fst = sd(fst, na.rm = TRUE),
    n_comparisons = sum(!is.na(fst)),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_fst), desc(max_fst))

fwrite(fst_summary, out_file, sep = "\t")

message("File saved to: ", out_file)
message("Top 10 SNPs by mean_fst:")
print(head(fst_summary, 10))
```

Run it:

```bash
Rscript 09_feature_integration/fst_summary.R \
  02_qc/Sace_qc.bim \
  06_fst/fast_results \
  09_feature_integration/fst_summary.tsv
```

---

## 9. Allele-frequency delta

### 9.1. Allele frequency by group

```bash
plink \
  --bfile 02_qc/Sace_qc \
  --allow-extra-chr \
  --freq \
  --within 02_qc/groups.clst \
  --out 08_delta/Sace_groupfreq
```

### 9.2. Compute delta

Delta is calculated as:

```text
delta = |p_i - p_j|
```

where `p_i` and `p_j` are allele frequencies in two different groups.

Save the following script as `09_feature_integration/calc_delta.R`.

```r
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Usage: Rscript calc_delta.R <metadata.tsv> <frq.strat> <out.tsv>")
}

meta_file <- args[1]
frq_file  <- args[2]
out_file  <- args[3]

meta <- fread(meta_file)

if (!all(c("sample", "group") %in% names(meta))) {
  stop("metadata.tsv must contain columns: sample and group")
}

frq <- fread(frq_file)

required_cols <- c("CHR", "SNP", "CLST", "MAF")
missing_cols <- setdiff(required_cols, names(frq))

if (length(missing_cols) > 0) {
  stop("Missing columns in .frq.strat file: ", paste(missing_cols, collapse = ", "))
}

frq <- frq %>%
  as_tibble() %>%
  mutate(
    group = CLST,
    p = MAF
  ) %>%
  select(CHR, SNP, group, p)

groups <- sort(unique(frq$group))

if (length(groups) < 2) {
  stop("At least two groups are required to calculate delta.")
}

pair_tbl <- t(combn(groups, 2)) %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  setNames(c("group1", "group2")) %>%
  as_tibble()

delta_pairs <- pair_tbl %>%
  mutate(pair_id = paste(group1, group2, sep = "_vs_")) %>%
  pmap_dfr(function(group1, group2, pair_id) {
    x1 <- frq %>%
      filter(group == group1) %>%
      rename(p1 = p) %>%
      select(CHR, SNP, p1)

    x2 <- frq %>%
      filter(group == group2) %>%
      rename(p2 = p) %>%
      select(CHR, SNP, p2)

    full_join(x1, x2, by = c("CHR", "SNP")) %>%
      mutate(
        group1 = group1,
        group2 = group2,
        pair_id = pair_id,
        delta = abs(p1 - p2)
      )
  })

delta <- delta_pairs %>%
  group_by(CHR, SNP) %>%
  summarise(
    mean_delta = mean(delta, na.rm = TRUE),
    median_delta = median(delta, na.rm = TRUE),
    max_delta = max(delta, na.rm = TRUE),
    n_pairs = sum(!is.na(delta)),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_delta), desc(max_delta))

fwrite(delta, out_file, sep = "\t")

pair_out <- sub("\\.tsv$", "_pairwise.tsv", out_file)
fwrite(delta_pairs, pair_out, sep = "\t")

message("Main file saved to: ", out_file)
message("Pairwise comparisons saved to: ", pair_out)
message("Top 10 SNPs by mean_delta:")
print(head(delta, 10))
```

Run it:

```bash
Rscript 09_feature_integration/calc_delta.R \
  data/metadata.tsv \
  08_delta/Sace_groupfreq.frq.strat \
  09_feature_integration/delta.tsv
```

---

## 10. Integrate Fst and delta rankings

Save the following script as `09_feature_integration/integrate_fst_delta.R`.

```r
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Usage: Rscript integrate_fst_delta.R <delta.tsv> <fst.tsv> <out.tsv>")
}

delta_file <- args[1]
fst_file   <- args[2]
out_file   <- args[3]

delta <- fread(delta_file) %>% as_tibble()
fst   <- fread(fst_file) %>% as_tibble()

if (!"SNP" %in% names(delta)) {
  stop("delta.tsv must contain column: SNP")
}

if (!all(c("SNP", "mean_fst") %in% names(fst))) {
  stop("fst.tsv must contain columns: SNP and mean_fst")
}

x <- delta %>%
  left_join(fst, by = "SNP") %>%
  mutate(
    rank_delta = rank(-mean_delta, ties.method = "average"),
    rank_fst   = rank(-mean_fst, ties.method = "average"),
    rank_sum   = rank_delta + rank_fst
  ) %>%
  arrange(rank_sum, desc(mean_delta), desc(mean_fst))

fwrite(x, out_file, sep = "\t")

message("Output saved to: ", out_file)
message("Top 20 combined SNPs:")
print(head(x, 20))
```

Run it:

```bash
Rscript 09_feature_integration/integrate_fst_delta.R \
  09_feature_integration/delta.tsv \
  09_feature_integration/fst_summary.tsv \
  09_feature_integration/combined_fst_delta_ranking.tsv
```

---

## 11. ADMIXTURE analysis

ADMIXTURE requires integer chromosome codes. If the `.bim` file contains chromosome names such as `NC_001133.9`, they must be converted to integer-like labels before running ADMIXTURE.

### 11.1. Copy pruned PLINK files

```bash
cp 03_pruned/Sace_qc_pruned.bed 07_admixture/
cp 03_pruned/Sace_qc_pruned.bim 07_admixture/
cp 03_pruned/Sace_qc_pruned.fam 07_admixture/
```

### 11.2. Create chromosome map

Create `07_admixture/chr_map.tsv`:

```text
001133	NC_001133.9
001134	NC_001134.8
001135	NC_001135.5
001136	NC_001136.10
001137	NC_001137.3
001138	NC_001138.5
001139	NC_001139.9
001140	NC_001140.6
001141	NC_001141.2
001142	NC_001142.9
001143	NC_001143.9
001144	NC_001144.5
001145	NC_001145.3
001146	NC_001146.8
001147	NC_001147.6
001148	NC_001148.4
001224	NC_001224.1
```

### 11.3. Replace chromosome names in the `.bim` file

The `chr_map.tsv` file has two columns:

```text
integer_code	original_chromosome_name
```

Therefore, the replacement uses column 2 as the key and column 1 as the new value.

```bash
awk 'BEGIN{FS=OFS="\t"}
NR==FNR {map[$2]=$1; next}
{
  if ($1 in map) $1=map[$1]
  print
}' \
  07_admixture/chr_map.tsv \
  07_admixture/Sace_qc_pruned.bim \
  > 07_admixture/Sace_qc_pruned.integer_chr.bim

mv 07_admixture/Sace_qc_pruned.bim \
   07_admixture/Sace_qc_pruned.original_chr.bim

mv 07_admixture/Sace_qc_pruned.integer_chr.bim \
   07_admixture/Sace_qc_pruned.bim
```

### 11.4. Run ADMIXTURE for K = 2 to 8

```bash
cd 07_admixture

for K in 2 3 4 5 6 7 8; do
  echo "Running ADMIXTURE K=${K}"
  admixture --cv Sace_qc_pruned.bed "$K" | tee "log${K}.out"
done

cd ..
```

### 11.5. Extract cross-validation errors

```bash
grep -h "CV error" 07_admixture/log*.out | \
  sed -E 's/.*K=([0-9]+).*: ([0-9.]+)/\1\t\2/' \
  > 07_admixture/cv_errors.tsv
```

Output:

```text
07_admixture/cv_errors.tsv
```

---

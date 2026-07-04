# Three-Plot Tumor vs Healthy Methylation–Expression Pipeline

## Overview

This R script generates three PDF plots for a selected transcription factor, typically `BANP` or `NRF1`, using tumor-vs-healthy methylation and expression data.

The script is designed to compare, for each gene–motif–cancer pair:

1. DNA methylation changes between tumor and healthy samples.
2. Gene expression changes between tumor and healthy samples.
3. The joint relationship between methylation delta and expression log2 fold-change, colored by Pearson correlation.

The pipeline is intended for downstream visualization of methylation-sensitive transcription factor motif analyses, especially when integrating methylation, expression, and correlation results.

---

## Main Purpose

For a given transcription factor, the script produces:

### Plot 1: Methylation Volcano-like Plot

**x-axis:**
Mean tumor methylation − mean healthy methylation

**y-axis:**
`-log10(BH-adjusted Wilcoxon p-value)` for methylation

This plot highlights gene–motif–cancer pairs with significant and large methylation changes.

---

### Plot 2: Expression Volcano-like Plot

**x-axis:**
Expression log2 fold-change:

```r
log2((mean Tumor TPM + pseudocount) / (mean Healthy TPM + pseudocount))
```

**y-axis:**
`-log10(BH-adjusted Wilcoxon p-value)` for expression

This plot highlights gene–motif–cancer pairs with significant and large expression changes.

---

### Plot 3: Methylation Delta vs Expression log2FC

**x-axis:**
Mean tumor methylation − mean healthy methylation

**y-axis:**
Expression log2FC

Only pairs passing both methylation and expression significance/effect-size filters are shown.

Points are colored by Pearson correlation coefficient `r`.

CSV-selected gene–cancer pairs are circled in this plot.

---

## Input Argument

The script takes one command-line argument: the transcription factor name.

Example:

```bash
Rscript ./scripts/three_change_plots.R BANP
```

or:

```bash
Rscript ./scripts/three_change_plots.R NRF1
```

If no argument is provided, the script defaults to:

```r
tf_name <- "BANP"
```

---

## Required R Packages

The script requires the following R packages:

```r
data.table
ggplot2
ggrepel
grid
```
---

## Expected Directory Structure

The script expects the following project structure:

```text
results/
└── methylation/
    └── correlation_expression_methylation_2d/
        └── tumor_vs_healthy/
            └── <TF_NAME>/
                ├── <TF_NAME>_GENES-Table.csv
                ├── all_samples/
                │   ├── motif_sample_expression_methylation_all.tsv.gz
                │   └── correlation_stats/
                │       └── pearson_correlation_all_pairs_by_cancer.tsv.gz
                └── matched_patients/
                    ├── motif_sample_expression_methylation_matched_patients.tsv.gz
                    └── correlation_stats/
                        └── pearson_correlation_all_pairs_by_cancer.tsv.gz
```

For example, for `BANP`, the expected base directory is:

```text
results/methylation/correlation_expression_methylation_2d/tumor_vs_healthy/BANP/
```

---

## Required Input Files

### 1. Correlation file

For all samples:

```text
all_samples/correlation_stats/pearson_correlation_all_pairs_by_cancer.tsv.gz
```

For matched patients:

```text
matched_patients/correlation_stats/pearson_correlation_all_pairs_by_cancer.tsv.gz
```

Required columns include:

```text
gene
motif_id
cancer
```

The Pearson correlation column can be named one of:

```text
pearson_r
r
correlation
```

The Pearson FDR column is optional. If present, it can be named one of:

```text
pearson_fdr
padj_BH
fdr
FDR
adj_p
adjusted_pvalue
```

If a column named `tested` exists, only rows with `tested == TRUE` are retained.

---

### 2. Merged methylation-expression file

For all samples:

```text
all_samples/motif_sample_expression_methylation_all.tsv.gz
```

For matched patients:

```text
matched_patients/motif_sample_expression_methylation_matched_patients.tsv.gz
```

Required columns include:

```text
gene
motif_id
cancer
sample_type
```

Expression data must be provided in one of these columns:

```text
expression
log2_expr
```

Methylation data must be provided in one of these columns:

```text
meth_percent
meth_beta
```

If `meth_beta` ranges from 0 to 1, it is automatically converted to percentage points by multiplying by 100.

---

### 3. CSV-selected label file

The script also requires:

```text
<TF_NAME>_GENES-Table.csv
```

Example:

```text
BANP_GENES-Table.csv
NRF1_GENES-Table.csv
```

This file is read using semicolon separation:

```r
fread(label_file, sep = ";", fill = TRUE)
```

Required columns:

```text
gene
Cancer
Significance
```

Only rows where `Significance` is:

```text
YES
MAYBE
```

are kept.

These selected gene–cancer pairs are circled only in Plot 3.

---

## Main Parameters

The most important user-adjustable parameters are defined near the beginning of the script.

### Wilcoxon FDR cutoff

```r
wilcox_fdr_cutoff <- 0.05
```

This threshold is used for both methylation and expression significance.

---

### Pearson correlation cutoff

```r
pearson_r_cutoff <- -0.3
```

This defines anti-correlated pairs as:

```r
pearson_r < -0.3
```

By default, anti-correlation is calculated but not required for Plot 3 unless activated manually.

---

### Minimum sample number per group

```r
min_n_per_group <- 2
```

A Wilcoxon test is only performed if both the healthy and tumor groups have at least this number of valid samples.

---

### Minimum methylation delta

```r
min_abs_meth_delta <- 10
```

This means that a pair must have at least 10 percentage points difference in methylation to be considered a large methylation change.

Set to `0` to disable this filter.

---

### Minimum expression log2FC

```r
min_abs_log2fc <- 1
```

This means that a pair must have at least an absolute log2FC of 1, corresponding to at least a two-fold expression change.

Set to `0` to disable this filter.

---

### Pseudocount for expression log2FC

```r
pseudo_count <- 1
```

Expression log2FC is calculated as:

```r
log2((mean Tumor TPM + 1) / (mean Healthy TPM + 1))
```

---

### Require anti-correlation for Plot 3

```r
require_anti_correlation_for_plot3 <- FALSE
```

If set to `TRUE`, Plot 3 only keeps pairs that are:

```text
methylation Wilcoxon FDR < 0.05
expression Wilcoxon FDR < 0.05
abs(methylation delta) >= threshold
abs(expression log2FC) >= threshold
Pearson r < -0.3
```

By default, this is set to `FALSE`, so Plot 3 does not require anti-correlation.

---

## Manual Highlighting

The script allows manual highlighting of selected gene–cancer pairs.

Example:

```r
extra_label_pairs <- list(
  list(gene = "GGT6", cancer = "KIRC")
)
```

Additional pairs can be added as:

```r
extra_label_pairs <- list(
  list(gene = "GGT6", cancer = "KIRC"),
  list(gene = "TDRD1", cancer = "PRAD")
)
```

Manual pairs can be highlighted independently in each plot using:

```r
label_manual_pairs_plot1 <- TRUE
label_manual_pairs_plot2 <- TRUE
label_manual_pairs_plot3 <- TRUE
```

In Plots 1 and 2, manual pairs are circled and labeled in red.

In Plot 3, manual pairs are labeled in black.

CSV-selected pairs are circled only in Plot 3.

---

## Output Directory

All output files are saved to:

```text
results/methylation/correlation_expression_methylation_2d/tumor_vs_healthy/<TF_NAME>/delta_log2fc_wilcoxon_plots/
```

---

## Output Files

The script runs the same plotting pipeline twice:

1. On all samples.
2. On matched patients.

For each mode, it generates three PDF plots and one output table.

---

### For all samples

#### Plot 1

```text
<TF_NAME>_all_samples_01_methylation_delta_vs_methylation_wilcoxon_FDR.pdf
```

#### Plot 2

```text
<TF_NAME>_all_samples_02_expression_log2FC_vs_expression_wilcoxon_FDR.pdf
```

#### Plot 3

```text
<TF_NAME>_all_samples_03_methylation_delta_vs_expression_log2FC_ONLY_both_wilcoxon_significant_delta<MIN_METH_DELTA>_expr<MIN_LOG2FC>_colored_by_pearson_r.pdf
```

#### Plot table

```text
<TF_NAME>_all_samples_delta_log2fc_wilcoxon_pearson_plot_table.tsv.gz
```

---

### For matched patients

#### Plot 1

```text
<TF_NAME>_matched_patients_01_methylation_delta_vs_methylation_wilcoxon_FDR.pdf
```

#### Plot 2

```text
<TF_NAME>_matched_patients_02_expression_log2FC_vs_expression_wilcoxon_FDR.pdf
```

#### Plot 3

```text
<TF_NAME>_matched_patients_03_methylation_delta_vs_expression_log2FC_ONLY_both_wilcoxon_significant_delta<MIN_METH_DELTA>_expr<MIN_LOG2FC>_colored_by_pearson_r.pdf
```

#### Plot table

```text
<TF_NAME>_matched_patients_delta_log2fc_wilcoxon_pearson_plot_table.tsv.gz
```

---

## Statistical Methods

### Methylation delta

Methylation delta is calculated as:

```r
mean_meth_tumor - mean_meth_healthy
```

The value is expressed in methylation percentage points.

A positive value means that the motif is more methylated in tumor samples.

A negative value means that the motif is less methylated in tumor samples.

---

### Expression log2FC

Expression log2FC is calculated from raw TPM-like values:

```r
log2((mean Tumor TPM + pseudocount) / (mean Healthy TPM + pseudocount))
```

A positive value means that the gene is more expressed in tumor samples.

A negative value means that the gene is less expressed in tumor samples.

---

### Expression delta on log2 scale

The script also saves:

```r
expr_delta_log2 = mean log2(TPM+1) Tumor - mean log2(TPM+1) Healthy
```

This value is saved in the output table for reference but is not used for plotting.

---

### Wilcoxon tests

The script performs Wilcoxon rank-sum tests for:

1. Methylation differences between tumor and healthy samples.
2. Expression differences between tumor and healthy samples.

The expression Wilcoxon test is performed on:

```r
log2(TPM + 1)
```

The methylation Wilcoxon test is performed on methylation percentage values.

---

### Multiple testing correction

Benjamini–Hochberg FDR correction is applied independently within each cancer type:

```r
p.adjust(p, method = "BH")
```

This is done separately for expression and methylation Wilcoxon p-values.

---

### Pearson correlation

Pearson correlation values are loaded from the upstream correlation file.

The script does not recompute Pearson correlation.

The Pearson correlation coefficient is used to color points in Plot 3.

---

## Filtering Logic

### Plot 1

A point is considered methylation-relevant if:

```text
methylation Wilcoxon FDR < 0.05
abs(methylation delta) >= min_abs_meth_delta
```

These points are shown in black.

Other points are shown in grey.

---

### Plot 2

A point is considered expression-relevant if:

```text
expression Wilcoxon FDR < 0.05
abs(expression log2FC) >= min_abs_log2fc
```

These points are shown in black.

Other points are shown in grey.

---

### Plot 3

By default, a point is retained in Plot 3 if:

```text
methylation Wilcoxon FDR < 0.05
expression Wilcoxon FDR < 0.05
abs(methylation delta) >= min_abs_meth_delta
abs(expression log2FC) >= min_abs_log2fc
```

If `require_anti_correlation_for_plot3 <- TRUE`, the following condition is also required:

```text
Pearson r < pearson_r_cutoff
```

---

## Sample Type Cleaning

The script standardizes sample type labels.

The following values are converted to `Healthy`:

```text
healthy
normal
solid tissue normal
adjacent normal
11
11a
11x
```

The following values are converted to `Tumor`:

```text
tumor
tumour
primary tumor
primary tumour
cancer
01
01a
01x
```

Only cancer types with both tumor and healthy samples are retained.

---

## Console Output

During execution, the script prints useful information, including:

```text
Label file path
Number of CSV-selected gene-cancer pairs
Number of manually selected gene-cancer pairs
Cancers with both healthy and tumor samples
Cancers retained
Total plotted rows
Number of methylation significant rows
Number of expression significant rows
Number of both significant rows
Number of anti-correlated rows
Number of rows passing methylation delta threshold
Number of rows passing log2FC threshold
Number of final Plot 3 rows
Output file paths
```

This output is useful for checking whether the expected cancers and gene–cancer pairs were retained.

---

## Example Commands

Run the pipeline for `BANP`:

```bash
Rscript ./scripts/three_change_plots.R BANP
```

Run the pipeline for `NRF1`:

```bash
Rscript ./scripts/three_change_plots.R NRF1
```

If the actual script file is named differently, for example:

```text
all_corr_test2_FC.R
```

then use:

```bash
Rscript ./scripts/all_corr_test2_FC.R BANP
Rscript ./scripts/all_corr_test2_FC.R NRF1
```

---

## Important Notes

1. The script assumes that the working directory is the project root.
2. The script expects the results directory structure to already exist.
3. The Pearson correlation file is produced upstream and is not recalculated here.
4. The expression log2FC plotted in Plots 2 and 3 is the classical raw-TPM mean-based log2FC.
5. The Wilcoxon expression test is performed on `log2(TPM+1)`.
6. Methylation values are converted to percentage points if the input is beta values between 0 and 1.
7. CSV-selected pairs are circled only in Plot 3.
8. Manual pairs can be circled/labeled in Plots 1 and 2 and labeled in Plot 3.
9. Plot 3 may contain no points if no pairs pass the significance and effect-size filters.
10. If an input file is missing, the corresponding all-samples or matched-patients analysis is skipped.

---

## Troubleshooting

### Error: label file not found

Check that the following file exists:

```text
results/methylation/correlation_expression_methylation_2d/tumor_vs_healthy/<TF_NAME>/<TF_NAME>_GENES-Table.csv
```

Example:

```text
results/methylation/correlation_expression_methylation_2d/tumor_vs_healthy/BANP/BANP_GENES-Table.csv
```

---

### Error: missing required columns in label file

The label file must contain:

```text
gene
Cancer
Significance
```

The file must be semicolon-separated.

---

### Error: `motif_id` not found

The correlation and merged methylation-expression files must both contain:

```text
motif_id
```

This column is required to merge methylation-expression statistics with Pearson correlation results.

---

### No rows after merging/filtering

This usually means that the `gene`, `motif_id`, or `cancer` values do not match between the correlation file and the merged methylation-expression file.

Check that cancer names are consistent. The script removes the `TCGA-` prefix automatically.

---

### Plot 3 is empty

This means no gene–motif–cancer pair passed all Plot 3 filters.

To make the filter less strict, consider adjusting:

```r
min_abs_meth_delta <- 5
min_abs_log2fc <- 0.2
```

or disabling the anti-correlation requirement:

```r
require_anti_correlation_for_plot3 <- FALSE
```

---

## Recommended Citation in Analysis Notes

This script supports tumor-vs-healthy visualization of methylation-sensitive transcription factor motif–gene pairs by integrating methylation delta, expression log2FC, Wilcoxon FDR, and Pearson methylation-expression correlation.

It is especially useful for identifying candidate TF motif–gene pairs where cancer-associated methylation changes at predicted binding motifs are associated with altered expression of nearby or assigned target genes.

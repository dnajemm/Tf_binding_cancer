# Single Gene–Cancer Pair Methylation–Expression Plot

## Overview

This R script generates a detailed visualization for one transcription factor, one target gene, and one cancer type.

For a selected TF–gene–cancer pair, the script produces:

1. A scatterplot showing the relationship between motif methylation and gene expression.
2. A side boxplot comparing gene expression between healthy and tumor samples.
3. A bottom boxplot comparing motif methylation between healthy and tumor samples.
4. A small statistics table containing Pearson correlation and Wilcoxon test results.

The script is useful for manually inspecting specific candidate motif–gene pairs identified in a larger methylation-sensitive transcription factor analysis.

---

## Example Use Case

Example command:

```bash
Rscript ./scripts/correlation_methylation_oneplot.R BANP GGT6 KIRC
```

This generates a plot for:

```text
TF: BANP
Gene: GGT6
Cancer: KIRC
```

An optional motif can also be provided:

```bash
Rscript ./scripts/correlation_methylation_oneplot.R BANP GGT6 KIRC "chr17:4560723:4560733:+"
```

If no motif is provided, the script automatically selects the most anti-correlated motif for the selected gene and cancer.

---

## Command-Line Arguments

The script accepts up to four arguments:

```bash
Rscript ./scripts/correlation_methylation_oneplot.R <TF_NAME> <GENE> <CANCER> <MOTIF_ID>
```

### Argument 1: TF name

Example:

```text
BANP
NRF1
```

If not provided, the default is:

```r
tf_name <- "BANP"
```

The TF name is converted to uppercase.

---

### Argument 2: Target gene

Example:

```text
GGT6
TDRD1
SFRP1
```

If not provided, the default is:

```r
target_gene <- "GGT6"
```

---

### Argument 3: Target cancer

Example:

```text
KIRC
PRAD
BRCA
LUAD
```

If not provided, the default is:

```r
target_cancer <- "KIRC"
```

The cancer name is converted to uppercase.

---

### Argument 4: Target motif

This argument is optional.

Example:

```text
chr17:4560723:4560733:+
```

If no motif is provided, the script searches all motifs available for the selected gene–cancer pair and chooses the motif with the strongest negative Pearson correlation between methylation and expression.

---

## Required R Packages

The script requires the following packages:

```r
data.table
ggplot2
grid
gtable
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
                ├── all_samples/
                │   └── motif_sample_expression_methylation_all.tsv.gz
                └── matched_patients/
                    └── motif_sample_expression_methylation_matched_patients.tsv.gz
```

For example, for `BANP`, the expected all-samples input file is:

```text
results/methylation/correlation_expression_methylation_2d/tumor_vs_healthy/BANP/all_samples/motif_sample_expression_methylation_all.tsv.gz
```

---

## Data Mode

The script contains a hard-coded parameter:

```r
data_mode <- "all_samples"
```

This means that, by default, the script uses the all-samples dataset.

To use matched patients instead, change:

```r
data_mode <- "all_samples"
```

to:

```r
data_mode <- "matched_patients"
```

---

## Required Input File

Depending on `data_mode`, the script reads one of the following files.

### For all samples

```text
<TF_NAME>/all_samples/motif_sample_expression_methylation_all.tsv.gz
```

### For matched patients

```text
<TF_NAME>/matched_patients/motif_sample_expression_methylation_matched_patients.tsv.gz
```

---

## Required Columns

The merged methylation-expression file must contain the following columns:

```text
gene
motif_id
sample_barcode
sample_type
cancer
meth_beta
expression
```

The script stops with an error if any of these columns are missing.

---

## Sample Type Cleaning

The script standardizes sample type labels.

The following values are converted to `Healthy`:

```text
healthy
normal
solid tissue normal
adjacent normal
```

The following values are converted to `Tumor`:

```text
tumor
primary tumor
primary tumour
cancer
```

Only samples classified as `Healthy` or `Tumor` are retained.

---

## Methylation Processing

The input methylation column is:

```text
meth_beta
```

If methylation values appear to be beta values between 0 and 1, the script converts them to percentages:

```r
meth_percent <- meth_beta * 100
```

If methylation values are already above the 0–1 beta range, they are kept as they are.

The final methylation variable used for plotting is:

```text
meth_percent
```

---

## Expression Processing

The input expression column is:

```text
expression
```

The script transforms expression as:

```r
log2_expr <- log2(expression + 1)
```

This transformed value is used for:

1. The scatterplot y-axis.
2. The expression boxplot.
3. The expression Wilcoxon test.
4. The Pearson correlation with methylation.

---

## Motif Selection Logic

If the user provides a motif ID, the script uses that motif directly.

Example:

```bash
Rscript ./scripts/correlation_methylation_oneplot.R BANP GGT6 KIRC "chr17:4560723:4560733:+"
```

If no motif is provided, the script calculates Pearson correlation between methylation and expression for each motif available for the selected gene–cancer pair.

It then selects the motif with the most negative Pearson correlation.

In other words, it chooses the motif where higher methylation is most strongly associated with lower expression.

---

## Statistical Tests

### Pearson correlation

The script calculates Pearson correlation between:

```text
Motif methylation percentage
```

and:

```text
Gene expression log2(TPM + 1)
```

Minimum number of samples required:

```r
min_n_cor <- 3L
```

If there are fewer than 3 usable samples, Pearson correlation is returned as `NA`.

The plot displays:

```text
n
Pearson r
Pearson p-value
```

---

### Wilcoxon test for expression

The script compares gene expression between healthy and tumor samples using a Wilcoxon rank-sum test.

Tested variable:

```text
log2(TPM + 1)
```

Minimum number of samples per group:

```r
min_n_wilcox <- 2L
```

---

### Wilcoxon test for methylation

The script compares motif methylation between healthy and tumor samples using a Wilcoxon rank-sum test.

Tested variable:

```text
meth_percent
```

Minimum number of samples per group:

```r
min_n_wilcox <- 2L
```

---

## P-value Display

P-values are formatted as scientific notation.

If a p-value is smaller than R’s numerical threshold, it is displayed as:

```text
<2.2e-16
```

Significance stars are assigned as:

```text
p < 0.001   ***
p < 0.01    **
p < 0.05    *
p >= 0.05   ns
```

---

## Plot Components

The final figure contains three main panels.

---

### 1. Main scatterplot

The scatterplot shows:

```text
x-axis: motif methylation (%)
y-axis: gene expression log2(TPM + 1)
```

Tumor samples are shown as red circles.

Healthy samples are shown as black-outlined triangles.

A dashed horizontal red line is drawn at:

```r
expression_threshold <- 1
```

This corresponds to:

```text
log2(TPM + 1) = 1
```

The scatterplot also includes Pearson correlation statistics in the top-right corner.

---

### 2. Expression boxplot

The side boxplot compares expression between:

```text
Healthy
Tumor
```

The tested variable is:

```text
log2(TPM + 1)
```

The expression Wilcoxon p-value and significance stars are displayed above the boxplot.

---

### 3. Methylation boxplot

The bottom boxplot compares motif methylation between:

```text
Healthy
Tumor
```

The tested variable is:

```text
motif methylation (%)
```

The methylation Wilcoxon p-value and significance stars are displayed next to the boxplot.

---

## Axis Scaling

The methylation axis is fixed from 0 to 100:

```r
x_lim <- c(0, 100)
x_breaks <- c(0, 25, 50, 75, 100)
```

This ensures that the scatterplot and methylation boxplot use the same methylation scale.

The expression axis is calculated dynamically based on the selected data.

The expression scale is shared between the scatterplot and the expression boxplot.

---

## Output Directory

For `data_mode <- "all_samples"`, outputs are saved to:

```text
results/methylation/correlation_expression_methylation_2d/tumor_vs_healthy/<TF_NAME>/all_samples/custom_gene_cancer_pair_plots/
```

For `data_mode <- "matched_patients"`, outputs are saved to:

```text
results/methylation/correlation_expression_methylation_2d/tumor_vs_healthy/<TF_NAME>/matched_patients/custom_gene_cancer_pair_plots/
```

---

## Output Files

For each selected TF–gene–cancer–motif combination, the script saves three files:

### 1. PDF plot

```text
<TF_NAME>__<GENE>__<CANCER>__<MOTIF_ID>_scatter_axis_boxplots.pdf
```

### 2. PNG plot

```text
<TF_NAME>__<GENE>__<CANCER>__<MOTIF_ID>_scatter_axis_boxplots.png
```

### 3. Statistics table

```text
<TF_NAME>__<GENE>__<CANCER>__<MOTIF_ID>_stats.tsv
```

The script sanitizes the filename using the `safe_name()` function, replacing problematic characters with underscores.

---

## Output Statistics Table

The statistics table contains:

```text
TF
gene
cancer
motif_id
n_total
n_healthy
n_tumor
pearson_n
pearson_r
pearson_p
expr_wilcox_p
expr_wilcox_label
meth_wilcox_p
meth_wilcox_label
```

This table allows the numerical results used in the plot to be retrieved without manually reading them from the figure.

---

## Important Parameters

### Minimum number of samples for Pearson correlation

```r
min_n_cor <- 3L
```

At least 3 usable samples are required for Pearson correlation.

---

### Minimum number of samples for Wilcoxon tests

```r
min_n_wilcox <- 2L
```

At least 2 healthy and 2 tumor samples are required for each Wilcoxon test.

---

### Expression threshold line

```r
expression_threshold <- 1
```

This adds a dashed horizontal line at:

```text
log2(TPM + 1) = 1
```

---

### Side boxplot size

```r
side_box_size <- unit(6, "cm")
```

This controls both:

1. The width of the expression boxplot.
2. The height of the methylation boxplot.

---

### Top margin size

```r
top_margin_size <- unit(1.5, "cm")
```

This controls the amount of white space above the title.

---

## Example Commands

### BANP–GGT6 in KIRC

```bash
Rscript ./scripts/correlation_methylation_oneplot.R BANP GGT6 KIRC
```

### NRF1–TDRD1 in PRAD

```bash
Rscript ./scripts/correlation_methylation_oneplot.R NRF1 TDRD1 PRAD
```

### BANP–GGT6 in KIRC with a specific motif

```bash
Rscript ./scripts/correlation_methylation_oneplot.R BANP GGT6 KIRC "chr17:4560723:4560733:+"
```

---

## Console Output

When the script runs, it prints useful information such as:

```text
Selected motif
Saved PDF path
Saved PNG path
Saved stats path
Done
```

If no motif is provided, the script also prints the most anti-correlated motif selected for the pair.

---

## Troubleshooting

### Error: data file not found

Check that the input file exists in the expected directory.

For all samples:

```text
results/methylation/correlation_expression_methylation_2d/tumor_vs_healthy/<TF_NAME>/all_samples/motif_sample_expression_methylation_all.tsv.gz
```

For matched patients:

```text
results/methylation/correlation_expression_methylation_2d/tumor_vs_healthy/<TF_NAME>/matched_patients/motif_sample_expression_methylation_matched_patients.tsv.gz
```

---

### Error: missing required columns

The input file must contain:

```text
gene
motif_id
sample_barcode
sample_type
cancer
meth_beta
expression
```

---

### Error: no rows found for gene and cancer

This means the selected gene–cancer pair is not present in the merged input file.

Check:

1. Gene spelling.
2. Cancer abbreviation.
3. Whether the pair exists in the selected `data_mode`.
4. Whether the cancer names include or exclude the `TCGA-` prefix.

The script automatically removes `TCGA-` from cancer names.

---

### Error: no motif had enough usable data

This means no motif for the selected gene–cancer pair had enough finite methylation and expression values to calculate Pearson correlation.

Possible reasons:

1. Too few samples.
2. Missing methylation values.
3. Missing expression values.
4. Constant methylation or expression values.
5. Incorrect gene or cancer selection.

---

### Plot generated but Wilcoxon p-value is NA

This usually means there are too few healthy or tumor samples.

The script requires at least:

```text
2 healthy samples
2 tumor samples
```

for each Wilcoxon test.

---


## Important Notes

1. The script analyzes one TF–gene–cancer pair at a time.
2. If no motif is provided, the most anti-correlated motif is selected automatically.
3. The script currently uses `data_mode <- "all_samples"` unless manually changed.
4. Methylation is plotted as a percentage from 0 to 100.
5. Expression is plotted as `log2(TPM + 1)`.
6. Pearson correlation is computed across all retained samples for the selected motif.
7. Wilcoxon tests compare healthy vs tumor samples.
8. The script outputs both PDF and PNG versions of the plot.
9. The script also saves a TSV file with the statistics used in the plot.
10. This script is intended for visual inspection and candidate validation, not for genome-wide multiple-testing correction.

---

## Recommended Use in the Project

This script is best used after identifying candidate TF motif–gene pairs from the larger methylation-expression correlation pipeline.

For example, after detecting a candidate pair such as:

```text
BANP — GGT6 — KIRC
```

this script can be used to visually confirm whether:

1. Motif methylation differs between tumor and healthy samples.
2. Gene expression differs between tumor and healthy samples.
3. Methylation and expression are negatively correlated.
4. The selected motif appears biologically consistent with the expected regulatory model.

This makes the script useful for producing focused validation figures for presentations, reports, and candidate interpretation on one gene motif cancer pair.

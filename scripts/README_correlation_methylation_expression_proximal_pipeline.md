# TF-Centered Proximal Methylation-Expression Correlation Pipeline

## Overview

`correlation_methylation_expression_proximal_pipeline.R` performs a transcription factor-centered pan-cancer analysis linking DNA methylation at proximal motif-associated CpG sites to target-gene expression across the TCGA `2d_noOV` cohort.

Given a supported transcription factor such as `BANP` or `NRF1`, the pipeline:

1. builds a **proximal** motif-probe map from motif-associated probe-gene pairs,
2. annotates HM450 methylation samples from TCGA tumor (`01X`) and healthy (`11X`) groups,
3. computes sample-level motif methylation using the **maximum beta value** across CpGs linked to each motif instance,
4. merges motif methylation with RNA-seq expression for the same gene and sample,
5. computes pooled and cancer-specific **Pearson correlations with BH/FDR correction**,
6. extracts significantly anti-correlated gene-motif pairs,
7. builds a matched-patient subset restricted to patients with both healthy and tumor samples,
8. optionally applies by-cancer Healthy vs Tumor **Wilcoxon** filtering,
9. generates PDF reports for anti-correlated pairs in both the full and matched-patient analyses.

The pipeline writes static result tables and PDF reports under a TF-specific output directory.

## Supported Transcription Factors

- `BANP`
- `NRF1`

These are mapped internally to the motif labels:

- `BANP -> BANP_mm0to2_noCGmm`
- `NRF1 -> NRF1_mm0to2_noCGmm`

## Usage

```bash
Rscript ./scripts/correlation_methylation_expression_proximal_pipeline.R NRF1
```

If no transcription factor is provided, the pipeline defaults to:

```text
BANP
```

## Input Files

The script expects the following input files and directories:

- `./results/multi_omics/samples_2d_noOV.tsv` with this binary structure: sample  cancer  snv     methylation     rnaseq  atac
head ./results/multi_omics/samples_2d_noOV.tsv
TCGA-OR-A5J1    ACC     1       1       1       0
TCGA-OR-A5J2    ACC     1       1       1       1
TCGA-OR-A5J3    ACC     1       1       1       1
TCGA-OR-A5J5    ACC     1       1       1       0
TCGA-OR-A5J6    ACC     1       1       1       1
TCGA-OR-A5J7    ACC     1       1       1       0

- `./methylation/probe_gene_pairs_in_motifs_with_tss.tsv.gz` with this structure:
zcat ./methylation/probe_gene_pairs_in_motifs_with_tss.tsv.gz | head
TF      probe   gene    motif_instance_id       dist_to_tss     is_proximal     proximal_or_distal
NRF1_mm0to2_noCGmm      cg00000321      SFRP1   chr8:41310275:41310285:+        786     TRUE    proximal
NRF1_mm0to2_noCGmm      cg00000622      NIPA2   chr15:22838616:22838626:+       20      TRUE    proximal
BANP_mm0to2_noCGmm      cg00003345      CASZ1   chr1:10756438:10756448:+        40205   FALSE   distal
NRF1_mm0to2_noCGmm      cg00012362      S100A13 chr1:153627009:153627019:+      36      TRUE    proximal
NRF1_mm0to2_noCGmm      cg00014996      USP14   chr18:158338:158348:+   34      TRUE    proximal
NRF1_mm0to2_noCGmm      cg00015261      RP11-812E19.9   chr16:35795616:35795626:+       1950397 FALSE   distal
NRF1_mm0to2_noCGmm      cg00020474      ERRFI1  chr1:8154893:8154903:+  128595  FALSE   distal
BANP_mm0to2_noCGmm      cg00021152      DPYSL2  chr8:26509568:26509578:+        4449    FALSE   distal
NRF1_mm0to2_noCGmm      cg00032419      TP53I13 chr17:29568544:29568554:+       42      TRUE    proximal

- `./methylation/filtered_methylation/` with this structure:
ls ./methylation/filtered_methylation/ | head
HM450_TCGA-ACC-TCGA-OR-A5J1-01A_1_annotated_methylation_filtered.bed.gz
HM450_TCGA-ACC-TCGA-OR-A5J2-01A_1_annotated_methylation_filtered.bed.gz
HM450_TCGA-ACC-TCGA-OR-A5J3-01A_1_annotated_methylation_filtered.bed.gz
HM450_TCGA-ACC-TCGA-OR-A5J4-01A_1_annotated_methylation_filtered.bed.gz

- `./expression/gene_expression_matrix_2d_noOV.tsv` with this structure example for one first column (raw expression value tpm):
sample  DPM1 ...
TCGA-ACC-TCGA-OR-A5J1-01A_1     42.9586
TCGA-ACC-TCGA-OR-A5J2-01A_1     118.2528
TCGA-ACC-TCGA-OR-A5J3-01A_1     53.7845

- `./results/multi_omics/cancer_color_order_with_defined_colours.tsv` with this structure:
head ./results/multi_omics/cancer_color_order_with_defined_colours.tsv
TCGA-ACC        #107ef3
TCGA-BLCA       #25f209
TCGA-BRCA       #a95807
TCGA-CESC       #ab047c
TCGA-CHOL       #06376b

## Output Directory

Results are written under:

```text
./results/methylation/correlation_expression_methylation_2d/tumor_vs_healthy/<TF>/
```

with two major branches:

- `all_samples/`
- `matched_patients/`

## Analysis Workflow

### 1. Build proximal motif-probe map

The pipeline filters the probe-gene motif table to:

- proximal associations only (`is_proximal == TRUE`)
- the selected TF motif label

It then creates a unique long-format table with:

- `gene`
- `motif_id`
- `probe`

Output:

- `<TF>_proximal_motif_probe_map_long.tsv.gz`

### 2. Build methylation sample annotation

Methylation sample files are collected from `./methylation/filtered_methylation/` and annotated with:

- `patient_id`
- `sample_barcode`
- `sample_type`
- `cancer`
- `file`

Only patients present in the `2d_noOV` cohort file are kept.

Output:

- `methylation_sample_annotation_3d4d_noOV.tsv.gz`

Note:
The filename still contains `3d4d_noOV`, but the cohort used by the script is `2d_noOV`.

### 3. Compute motif-level methylation

For each methylation sample:

- the script reads HM450 probe beta values,
- keeps only probes linked to the selected TF proximal motifs,
- merges probe values with the motif-probe map,
- summarizes each `gene` / `motif_id` pair using the maximum beta value across its linked probes.

This pipeline uses:

```text
motif methylation = MAX(beta)
```

It also records:

- number of probes found per motif-sample pair,
- number of probes tied at the maximum beta value,
- the selected max-beta probes.

Outputs:

- `motif_sample_methylation_MAX_3d4d_noOV.tsv.gz`
- `motif_sample_methylation_MAX_selected_probes_3d4d_noOV.tsv.gz`
- `motif_sample_methylation_MAX_selected_probes_compact_3d4d_noOV.tsv.gz`
- `motif_sample_methylation_MAX_status_3d4d_noOV.tsv.gz`

### 4. Merge methylation with expression

The expression matrix is reshaped to long format, TCGA sample barcodes are extracted, and motif methylation is merged with expression by:

- `gene`
- `sample_barcode`

Output:

- `motif_sample_expression_methylation_all.tsv.gz`

### 4b. Build matched-patient subset

The matched-patient branch keeps only `gene` / `motif_id` / `patient_id` combinations where both of the following are present:

- at least one healthy sample
- at least one tumor sample

Output:

- `matched_patients/motif_sample_expression_methylation_matched_patients.tsv.gz`

### 5. Run correlation analysis

The pipeline computes Pearson correlations between:

- motif methylation (`meth_beta`)
- gene expression (`log2(expression + 1)`)

Correlations are computed:

- across all cancers combined,
- separately within each cancer type.

False discovery rate correction uses the Benjamini-Hochberg method:

- pooled correlations: corrected across all tested pairs,
- by-cancer correlations: corrected within each cancer type.

Outputs per branch:

- `pearson_correlation_all_pairs_all_cancers.tsv.gz`
- `pearson_correlation_all_pairs_by_cancer.tsv.gz`

### 6. Extract anti-correlated pairs

Significant anti-correlated pairs are defined as:

- `pearson_r < -0.3`
- `FDR < 0.05`

The pipeline extracts:

- pooled anti-correlated pairs,
- by-cancer anti-correlated pairs.

Outputs per branch:

- `anti_correlated/significant_anti_correlated_pairs_all_cancers.tsv.gz`
- `anti_correlated/significant_anti_correlated_pairs_by_cancer.tsv.gz`

### 7. Optional Wilcoxon Healthy vs Tumor filter

For by-cancer anti-correlated pairs, the pipeline can additionally test whether healthy and tumor samples differ for:

- expression
- methylation

using Wilcoxon rank-sum tests.

Current defaults in the script:

- Wilcoxon filtering enabled: `TRUE`
- Wilcoxon FDR cutoff: `0.05`
- Wilcoxon mode: `both`

Available modes:

- `both`: require both expression and methylation to be significant
- `either`: require either expression or methylation to be significant
- `expression`: require only expression to be significant
- `methylation`: require only methylation to be significant

If both healthy and tumor groups are not available in sufficient numbers for a cancer, the Wilcoxon rule is not applied to that case.

Output per branch:

- `anti_correlated/wilcox_healthy_vs_tumor_by_cancer.tsv.gz`

### 8. Generate PDF reports

The pipeline generates a PDF report for each selected anti-correlated `gene` / `motif_id` pair.

Reports include:

- all-cancer scatter plots of methylation vs expression,
- per-cancer scatter plots,
- dual-panel summaries of expression and methylation in healthy vs tumor samples,
- matched-patient versions of the same report structure for the matched branch.

Report index files are also written to summarize which plots were generated.

Outputs per branch:

- `anti_correlated_pair_plots/selected_anti_correlated_pairs_all_<TF>.tsv`
- `anti_correlated_pair_plots/generated_reports_<TF>.tsv`
- `anti_correlated_pair_plots/*_report.pdf`
- `anti_correlated_pair_plots/*_summary.tsv`

## Key Parameters

The current script defines the following main parameters:

```r
mc_cores        <- 20L
batch_size_meth <- 20L
batch_size_plot <- 5L
plot_cores      <- 5L
min_n_cor       <- 3L
min_n_wilcox    <- 2L
anti_r_cutoff   <- -0.3
anti_fdr_cutoff <- 0.05
only_significant_cancers <- FALSE
reuse_existing_intermediate_files <- TRUE
wilcox_filter_enabled <- TRUE
wilcox_fdr_cutoff     <- 0.05
wilcox_filter_mode    <- "both"
```

## Main Output Structure

Typical output layout:

```text
results/
  methylation/
    correlation_expression_methylation_2d/
      tumor_vs_healthy/
        BANP/
          all_samples/
            correlation_stats/
            anti_correlated_pair_plots/
          matched_patients/
            correlation_stats/
            anti_correlated_pair_plots/
        NRF1/
          all_samples/
          matched_patients/
```

## Important Notes

- The analysis uses `log2(expression + 1)` for correlation and plotting.
- Methylation values are plotted as percentages when beta values are in the `0-1` range.
- Existing intermediate files are reused when `reuse_existing_intermediate_files = TRUE`.
- The plotting step can include all cancers for a pair, not only cancers passing significance, because `only_significant_cancers` is currently `FALSE`.
- The matched-patient branch is not a paired statistical model; it is a subset restricted to patients represented in both healthy and tumor groups.

## Example Interpretation

For a given `gene` / `motif_id` pair, evidence in favor of a biologically interesting anti-correlation typically means:

- higher motif methylation is associated with lower gene expression,
- the association is significant after FDR correction,
- in the by-cancer setting, healthy and tumor samples may also differ significantly in expression and methylation.

## Script Entry Point

The pipeline is executed through:

```r
run_pipeline()
```

which runs the full workflow for the selected transcription factor across:

- all available samples
- the matched-patient subset

# TF-Centered Correlation Pipeline

This R pipeline performs a TF-centered pan-cancer analysis that links DNA methylation at proximal motif-associated CpG sites to target-gene expression across the TCGA `3d+4d_noOV` cohort using `01X` tumor and `11X` healthy samples.

Given a transcription factor such as `BANP` or `NRF1`, the script:

- builds a proximal motif-to-probe map
- annotates HM450 methylation samples
- computes sample-level motif methylation as the maximum beta value across motif-linked probes
- merges motif methylation with RNA-seq expression
- computes pooled and cancer-specific Pearson correlations with BH/FDR correction
- builds matched-patient subsets requiring both healthy and tumor samples for the same patient, gene, and motif
- extracts significantly anti-correlated gene-motif pairs
- generates static PDF reports
- optionally generates interactive HTML scatter plots

## Supported TFs

- `BANP`
- `NRF1`

The script maps each TF argument to the internal motif label used in the probe-gene input table:

- `BANP -> BANP_mm0to2_noCGmm`
- `NRF1 -> NRF1_mm0to2_noCGmm`

## Required R Packages

The script loads these packages:

- `data.table`
- `ggplot2`
- `plotly`
- `htmlwidgets`
- `parallel`

## Usage

Save the script as an executable R script, for example `tf_correlation_pipeline.R`, then run:

```bash
Rscript tf_correlation_pipeline.R NRF1
```

If no argument is provided, the default TF is `BANP`.

## Input Files

The pipeline expects these paths relative to the working directory:

- `./results/multi_omics/samples_3d+4d_noOV.tsv`
- `./methylation/probe_gene_pairs_in_motifs_with_tss.tsv.gz`
- `./methylation/filtered_methylation/`
- `./expression/gene_expression_matrix_3d4d_noOV.tsv`
- `./results/multi_omics/cancer_color_order_with_defined_colours.tsv`
- `./expression/expression_of_NRF1_BANP/all_samples_NRF1_BANP_expression.tsv`

## Output Layout

All results are written under:

```text
./results/methylation/correlation_expression_methylation/tumor_vs_healthy/<TF>/
```

with two main branches:

- `all_samples/`
- `matched_patients/`

## Main Parameters

These are the key user-facing options currently present in the script:

### Parallelism and batching

- `mc_cores <- 20L`
- `batch_size_meth <- 20L`
- `batch_size_plot <- 5L`
- `plot_cores <- 5L`

### Statistical thresholds

- `min_n_cor <- 3L`
- `min_n_wilcox <- 2L`
- `n_top <- 50L`
- `anti_r_cutoff <- -0.3`
- `anti_fdr_cutoff <- 0.05`

### Output behavior

- `make_interactive_html <- FALSE`
- `reuse_existing_intermediate_files <- TRUE`
- `only_significant_cancers <- FALSE`

### Wilcoxon filter for by-cancer anti-correlated pairs

- `wilcox_filter_enabled <- TRUE`
- `wilcox_fdr_cutoff <- 0.05`
- `wilcox_filter_mode <- "both"`

Allowed `wilcox_filter_mode` values:

- `"both"`: require both expression and methylation to pass Wilcoxon FDR
- `"either"`: require expression or methylation to pass
- `"expression"`: require only expression to pass
- `"methylation"`: require only methylation to pass

Important behavior:

- if a cancer lacks enough healthy or tumor samples for Wilcoxon testing, the filter is not enforced for that group and the pair is allowed through
- Wilcoxon FDR is computed within each cancer type

## Pipeline Steps

### 1. Build proximal motif-probe map

The pipeline filters the probe-gene motif table to:

- `is_proximal == TRUE`
- `TF == tf_motif_label`

and writes a long table with:

- `gene`
- `motif_id`
- `probe`

Output:

- `all_samples/<TF>_proximal_motif_probe_map_long.tsv.gz`

### 2. Build methylation sample annotation

The script scans methylation BED files, extracts:

- sample barcode
- patient ID
- sample type
- cancer code

and keeps only patients listed in `samples_3d+4d_noOV.tsv`.

Output:

- `all_samples/methylation_sample_annotation_3d4d_noOV.tsv.gz`

### 3. Compute motif methylation

For each methylation sample:

- only probes present in the motif-probe map are kept
- methylation is summarized per `gene + motif_id` using `MAX(beta)`
- the probe or probes achieving the maximum are also recorded

Outputs:

- `all_samples/motif_sample_methylation_MAX_3d4d_noOV.tsv.gz`
- `all_samples/motif_sample_methylation_MAX_selected_probes_3d4d_noOV.tsv.gz`
- `all_samples/motif_sample_methylation_MAX_selected_probes_compact_3d4d_noOV.tsv.gz`
- `all_samples/motif_sample_methylation_MAX_status_3d4d_noOV.tsv.gz`

### 4. Merge methylation with expression

The expression matrix is reshaped to long format and matched to methylation by:

- `gene`
- `sample_barcode`

Output:

- `all_samples/motif_sample_expression_methylation_all.tsv.gz`

### 4b. Build matched-patient subset

The matched branch keeps rows where the same patient has both:

- a healthy sample
- a tumor sample

for the same:

- `gene`
- `motif_id`

Output:

- `matched_patients/motif_sample_expression_methylation_matched_patients.tsv.gz`

### 5. Run correlations

The pipeline computes Pearson correlations between:

- motif methylation
- target-gene expression

Results are produced for:

- pooled all-cancer data
- each cancer separately

BH/FDR correction is applied:

- globally for pooled results
- within each cancer for by-cancer results

Outputs inside each branch:

- `correlation_stats/pearson_correlation_all_pairs_all_cancers.tsv.gz`
- `correlation_stats/pearson_correlation_all_pairs_by_cancer.tsv.gz`

### 6. Plot top pairs

The pipeline selects the top `n_top` pairs from the by-cancer statistics and creates multi-page PDF reports.

These reports include:

- a pooled scatter plot across cancers
- per-cancer scatter plots
- pooled dual-panel expression/methylation box-and-jitter plots
- per-cancer dual-panel plots

Notable plotting behavior:

- expression is plotted as `log2(expression + 1)`
- methylation is shown as percent when beta values are on the 0 to 1 scale
- a red horizontal dashed line is drawn at `log2(expression + 1) = 1`
- the line is annotated with `expression = 1`
- healthy samples are triangles and tumor samples are circles in pooled scatter plots

Output directory inside each branch:

- `top_pair_plots/`

### 7. Matched-patient plots

The same top-pair plotting logic is repeated for the matched-patient subset.

Matched reports also record patient-level counts where available.

### 8 and 9. Extract anti-correlated pairs

Anti-correlated pairs are defined using:

- `pearson_r < anti_r_cutoff`
- `pearson_fdr < anti_fdr_cutoff`

Two result sets are produced:

- pooled anti-correlated pairs across all cancers
- by-cancer anti-correlated pairs

If `wilcox_filter_enabled <- TRUE`, the by-cancer set is further filtered using healthy-vs-tumor Wilcoxon tests on:

- target expression
- motif methylation

Outputs inside each branch:

- `correlation_stats/anti_correlated/significant_anti_correlated_pairs_all_cancers.tsv.gz`
- `correlation_stats/anti_correlated/significant_anti_correlated_pairs_by_cancer.tsv.gz`
- `correlation_stats/anti_correlated/wilcox_healthy_vs_tumor_by_cancer.tsv.gz`

### 10. Plot all anti-correlated pairs

The pipeline creates PDF reports for the union of:

- pooled anti-correlated pairs
- by-cancer anti-correlated pairs

Depending on `only_significant_cancers`:

- `FALSE`: all cancers present for the pair are plotted
- `TRUE`: only cancer types significant for that pair are plotted

Output directory inside each branch:

- `anti_correlated_pair_plots/`

### 11. Optional interactive HTML

If `make_interactive_html <- TRUE`, the pipeline generates interactive Plotly scatter plots for anti-correlated pairs.

These plots include hover information for:

- TF name
- target gene
- motif ID
- sample and patient identifiers
- cancer and sample type
- motif methylation
- target expression
- TF expression
- `NRF1` and `BANP` TPM values

Output directory inside each branch:

- `interactive_target_vs_methylation_only/`

## Reuse of Intermediate Files

If:

```r
reuse_existing_intermediate_files <- TRUE
```

the pipeline reuses previously generated intermediate files when they already exist. This is especially useful when you only want to rerun downstream filtering or plotting.

## Notes on Sample Handling

The script normalizes sample-type labels to:

- `Healthy`
- `Tumor`

It also relies on TCGA barcode parsing helpers to recover:

- sample barcodes
- patient IDs
- sample type codes
- cancer identifiers

Healthy samples correspond to code `11` and tumor samples correspond to code `01`.

## Example Output Tree

```text
results/methylation/correlation_expression_methylation/
└── tumor_vs_healthy/
    └── BANP/
        ├── all_samples/
        │   ├── BANP_proximal_motif_probe_map_long.tsv.gz
        │   ├── methylation_sample_annotation_3d4d_noOV.tsv.gz
        │   ├── motif_sample_methylation_MAX_3d4d_noOV.tsv.gz
        │   ├── motif_sample_expression_methylation_all.tsv.gz
        │   ├── correlation_stats/
        │   ├── top_pair_plots/
        │   ├── anti_correlated_pair_plots/
        │   └── interactive_target_vs_methylation_only/
        └── matched_patients/
            ├── motif_sample_expression_methylation_matched_patients.tsv.gz
            ├── correlation_stats/
            ├── top_pair_plots/
            ├── anti_correlated_pair_plots/
            └── interactive_target_vs_methylation_only/
```

## Practical Notes

- The script is designed for large TCGA-scale data and uses `mclapply` for parallel execution.
- Static PDFs are always generated by the plotting steps.
- Interactive HTML is optional and disabled by default.
- The Wilcoxon-based filtering only affects the by-cancer anti-correlated output, not the pooled anti-correlated table.
- The matched-patient branch requires healthy and tumor samples from the same patient for the same gene-motif pair.

## Summary

This pipeline is a TF-specific methylation-expression correlation workflow for `BANP` and `NRF1` that:

- quantifies **proximal** motif methylation using the maximum linked-probe beta
- relates methylation to target-gene expression across pan-cancer and cancer-specific contexts
- performs matched healthy-versus-tumor patient analyses
- filters anti-correlated pairs with optional Wilcoxon rules
- generates publication-style PDF reports and optional interactive HTML outputs

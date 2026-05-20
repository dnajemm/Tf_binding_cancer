# PCA + t-SNE of DNA Methylation Profiles

This project contains an R pipeline for dimensionality reduction of TCGA DNA methylation profiles. It builds a pan-cancer methylation matrix from per-sample methylation files, filters and ranks CpG probes by variability, runs PCA, then runs t-SNE on the PCA-reduced space.

The main script is ./scripts/tSNE_methylation_pipeline.R

## What the pipeline does

The pipeline performs these major tasks:

1. Reads a cohort file containing TCGA patient IDs:
(/data/najemd/TF_binding_cancer/conda_envs/TF_binding_cancer_env) najemd@serv-bardet-dorsal:/data/najemd/TF_binding_cancer$ head ./results/multi_omics/samples_2d_noOV_noCHOL.tsv
TCGA-OR-A5J1    ACC     1       1       1       0
TCGA-OR-A5J2    ACC     1       1       1       1
2. Scans a methylation directory for files ending in `_annotated_methylation_filtered.bed.gz`.
3. Parses each filename to recover cancer type, sample barcode, patient ID, sample code, sample type, and technical replicate number.
4. Keeps only samples that:
   - belong to patients listed in the cohort file
   - have TCGA sample codes `01` to `19`
   - are not from `OV` or `CHOL`
5. Detects duplicate `sample_barcode` entries caused by technical replicates and keeps one file per sample, preferring replicate `_1` over `_2`.
6. Builds a pan-cancer methylation matrix with:
   - rows = CpG probes
   - columns = samples
   - values = methylation beta values
7. Removes CpGs with `>=50%` missing values.
8. Median-imputes the remaining missing values probe-wise.
9. Removes CpGs with zero or non-finite variance.
10. Selects the top `N` most variable CpGs.
11. Runs PCA on the retained samples using the selected CpGs.
12. Runs t-SNE on the first `initial_dims` PCA components, capped by the number of available PCs and samples.
13. Produces publication-style plots and coordinate tables for downstream analysis.

## Expected input layout

The script assumes this directory structure relative to the project root:

```text
./methylation/filtered_methylation/
./results/multi_omics/cancer_color_order_with_defined_colours.tsv
./results/multi_omics/<cohort file>
```

### Cohort file

- Passed with `--cohort_file`
- Read with `data.table::fread()`
- The first column must contain TCGA patient IDs like `TCGA-XX-YYYY`
- The second column must contain the cancer type

### Methylation files

- Located under `./methylation/filtered_methylation`
- Must end with `_annotated_methylation_filtered.bed.gz`
- The script expects filenames to contain a TCGA sample barcode such as `TCGA-BL-A0C8-01A`
- A technical replicate suffix like `_1` or `_2` is also expected when replicate files exist

### Color file

- Located at `./results/multi_omics/cancer_color_order_with_defined_colours.tsv`
- Must contain at least two columns:
  - cancer label
  - color value
(/data/najemd/TF_binding_cancer/conda_envs/TF_binding_cancer_env) najemd@serv-bardet-dorsal:/data/najemd/TF_binding_cancer$ head ./results/multi_omics/cancer_color_order_with_defined_colours.tsv
TCGA-ACC        #107ef3
TCGA-BLCA       #25f209
TCGA-BRCA       #a95807
TCGA-CESC       #ab047c
TCGA-CHOL       #06376b
TCGA-COAD       #17720b
TCGA-DLBC       #6a4118
TCGA-ESCA       #f439bf
TCGA-GBM        #4E79A7

## Command-line arguments

The script accepts these parameters:

- `--cohort_file`  
  Required. Path to the cohort file.

- `--top_cpgs`  
  Optional. Number of most variable CpGs to keep. Default: `10000`

- `--seed`  
  Optional. Random seed for t-SNE. Default: `123`

- `--initial_dims`  
  Optional. Number of PCA dimensions to pass to t-SNE. Default: `30`

- `--max_iter`  
  Optional. Maximum t-SNE iterations. Default: `1000`

## Example run

```bash
Rscript ./scripts/tSNE_methylation_pipeline.R \
  --cohort_file ./results/multi_omics/samples_2d_noOV_noCHOL.tsv \
  --top_cpgs 10000 \
  --seed 123 \
  --initial_dims 30 \
  --max_iter 1000
```

## Step-by-step outputs

For each run, outputs are written to:

```text
./results/methylation/dimred_pipeline/<cohort_basename>_top<top_cpgs>/
```

The main output files are:

- `methylation_files_0X_1X.tsv`  
  Filtered sample table used by the pipeline.

- `duplicated_methylation_sample_barcodes.tsv`  
  Written only when technical replicate duplicates are detected.

- `pan_cancer_methylation_matrix.RData`  
  The raw pan-cancer CpG-by-sample methylation matrix.

- `top<top_cpgs>_most_variable_CpGs_matrix.RData`  
  The filtered and variance-ranked CpG matrix used for PCA.

- `top<top_cpgs>_variable_CpGs.txt`  
  Plain-text list of selected CpG probe IDs.

- `CpG_filtering_summary_top<top_cpgs>.tsv`  
  Summary of filtering, imputation, variance filtering, and final counts.

- `PCA_coordinates_top<top_cpgs>.tsv`  
  PCA coordinates with sample annotations.

- `PCA_top<top_cpgs>.pdf`  
  PCA scatterplots for `PC1 vs PC2` and `PC3 vs PC4`.

- `PCA_top<top_cpgs>_prcomp.rds`  
  Saved `prcomp` object.

- `tSNE_coordinates_top<top_cpgs>.tsv`  
  t-SNE coordinates with sample annotations and t-SNE parameters.

- `tSNE_top<top_cpgs>.pdf`  
  Main t-SNE plot colored by cancer type and shaped by sample type.

- `tSNE_top<top_cpgs>_rtsne.rds`  
  Saved `Rtsne` object.

- `tSNE_top<top_cpgs>_4views_onepage.pdf`  
  One-page 4-panel t-SNE summary figure.

- `tSNE_top<top_cpgs>_4views_onepage_annotated_coordinates.tsv`  
  Annotated t-SNE coordinates including organ-group and tissue-class labels.

## Sample annotations created by the pipeline

The script derives several useful labels directly from filenames:

- `cancer`
  Parsed from the methylation filename prefix

- `sample_barcode`
  Extracted using a TCGA barcode regex

- `patient_id`
  First 12 characters of the sample barcode

- `sample_code`
  Characters 14 to 15 of the sample barcode

- `sample_type`
  Assigned as:
  - `Tumor` for sample codes `01` to `09`
  - `Healthy` for sample codes `10` to `19`

- `sample_group`
  Assigned as:
  - `0X` for sample codes `01` to `09`
  - `1X` for sample codes `10` to `19`

- `technical_replicate`
  Parsed from the trailing `_1`, `_2`, etc. in the filename

## Visualization summary

The pipeline generates:

- PCA plots:
  - `PC1 vs PC2`
  - `PC3 vs PC4`

- One main t-SNE plot:
  - colored by cancer type
  - shaped by sample type

- One 4-panel t-SNE summary figure:
  - cancer type
  - tumor versus healthy
  - affected organ group
  - tissue / cancer class

The 4-panel plot uses internal mappings for:

- organ group
- tissue / cancer class

These mappings are hard-coded in the script.

## Important assumptions and caveats

- The pipeline is designed around TCGA-style sample naming conventions.
- It assumes every methylation file has compatible probe content and ordering can be anchored from the first file.
- The methylation matrix is built in memory, so large cohorts may require substantial RAM.
- The script uses `zcat` through `fread(cmd = ...)`, so the runtime environment needs shell access to `zcat`.
- Step skipping is based on whether expected output files already exist.
- The t-SNE perplexity is chosen automatically from sample count.
- The current script structure is strong, but small-cohort edge cases may still need guarding, especially when fewer than 4 principal components are available.

## R package dependencies

The script loads these packages:

- `data.table`
- `matrixStats`
- `ggplot2`
- `Rtsne`
- `gridExtra`
- `grid`

## In plain language

This pipeline starts from many compressed methylation files, keeps the cohort-relevant TCGA samples, removes problematic probes, picks the most informative CpGs, reduces the dimensionality with PCA and t-SNE, and saves both the numeric embeddings and ready-to-use figures for biological interpretation.

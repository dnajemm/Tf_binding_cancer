# ATAC-Seq Peak Presence/Absence t-SNE Analysis

This analysis script performs t-SNE analysis on ATAC-seq peak presence/absence data across TCGA cancer types. It uses sample-level ATAC peak BED files, builds a global merged peak universe, converts each sample into a binary peak presence/absence vector, selects the most variable peaks, runs PCA followed by t-SNE, and produces PDF visualizations colored by cancer type and broader biological groupings.

The analysis has two parts:

1. Build the ATAC peak presence/absence matrix and run t-SNE.
2. Replot the same t-SNE coordinates in three views: cancer type, organ group, and tissue type.

## Overview

For all ATAC-seq peak files in:

```text
./peaks/filtered_peaks/
```

the first script:

1. Finds TCGA ATAC peak BED files.
2. Extracts the cancer type from each filename.
3. Concatenates all peaks across samples.
4. Sorts and merges overlapping peaks using `bedtools merge`.
5. Builds a binary peak presence/absence matrix for all samples.
6. Computes variability of each merged peak across samples.
7. Keeps the top 10,000 most variable peaks.
8. Scales the binary feature matrix.
9. Runs PCA using up to 30 principal components.
10. Runs t-SNE on the PCA coordinates.
11. Saves t-SNE coordinates to a TSV file.
12. Saves a cancer-colored t-SNE scatter plot as a PDF.

The second script:

1. Reads the saved t-SNE coordinates.
2. Adds organ-group and tissue-type annotations.
3. Replots the same t-SNE embedding in three panels:
   - Cancer type
   - Organ group
   - Tissue type
4. Saves a one-page PDF containing all three views.

## Requirements

### Command-Line Tools

The scripts require:

```text
bedtools
sort
cut
cat
```

### R Packages

The scripts require:

```r
data.table
Rtsne
ggplot2
Polychrome
patchwork
grid
```

## Expected Input Files

### ATAC Peak BED Files

The first script reads ATAC peak files from:

```text
./peaks/filtered_peaks/
```

It keeps files matching this pattern:

```r
^ATAC_TCGA-[A-Z]+_.*_peaks_macs\\.bed$
```

Example filenames:

```text
ATAC_TCGA-BRCA_sample1_peaks_macs.bed
ATAC_TCGA-LUAD_sample2_peaks_macs.bed
ATAC_TCGA-KIRC_sample3_peaks_macs.bed
```

Each BED file should contain at least the first three standard BED columns:

```text
chromosome
start
end
```

Only columns 1-3 are used.

### Cancer Color File

Both scripts use:

```text
./results/multi_omics/cancer_color_order_with_defined_colours.tsv
```

This file should contain at least two columns:

| Column | Description |
| --- | --- |
| Column 1 | Cancer type, for example `TCGA-BRCA` or `BRCA` |
| Column 2 | Color value, for example `#FF0000` |

The scripts remove the `TCGA-` prefix before matching colors to cancer labels.

## Part 1: Run t-SNE From ATAC Peak Files

### Main Settings

```r
peaks_dir <- "./peaks/filtered_peaks/"
file_pattern <- "^ATAC_TCGA-[A-Z]+_.*_peaks_macs\\.bed$"

out_pdf <- "./results/peaks/tSNE_ATAC_ALLCANCERS_peakPresence_10000.pdf"
out_tsv <- "./results/peaks/tSNE_ATAC_ALLCANCERS_coords_10000.tsv"

top_k_peaks <- 10000
n_pcs <- 30
seed <- 1
```

### What It Produces

| Output | Description |
| --- | --- |
| `./results/peaks/tSNE_ATAC_ALLCANCERS_coords_10000.tsv` | t-SNE coordinates and cancer labels for each ATAC sample. |
| `./results/peaks/tSNE_ATAC_ALLCANCERS_peakPresence_10000.pdf` | Cancer-colored t-SNE scatter plot. |

The coordinate TSV contains:

| Column | Description |
| --- | --- |
| `sample` | BED filename basename. |
| `cancer` | Cancer type extracted from the filename. |
| `tSNE1` | First t-SNE coordinate. |
| `tSNE2` | Second t-SNE coordinate. |

## Part 2: Create Three-View t-SNE PDF

The second script reads existing t-SNE coordinates and creates a one-page PDF with three separate panels.

### Input

The second script currently uses:

```r
in_tsv <- "./results/peaks/tSNE_ATAC_ALLCANCERS_coords.tsv"
```

If you want it to use the output from Part 1, change this to:

```r
in_tsv <- "./results/peaks/tSNE_ATAC_ALLCANCERS_coords_10000.tsv"
```

### Output

```text
./results/peaks/tSNE_ATAC_3views.pdf
```

This PDF contains three panels:

| Panel | Coloring |
| --- | --- |
| Cancer type | Individual TCGA cancer type. |
| Organ group | Manually defined organ/system group. |
| Tissue type | Broad tissue or histological class. |

## Organ Group Mapping

The second script maps cancer types to organ groups:

| Organ group | Cancer types |
| --- | --- |
| `CNS` | `GBM`, `LGG` |
| `Breast` | `BRCA` |
| `Gynecologic` | `CESC`, `OV`, `UCEC`, `UCS` |
| `Lung` | `LUAD`, `LUSC` |
| `HeadNeck` | `HNSC` |
| `ThyroidThymus` | `THCA`, `THYM` |
| `SkinEye` | `SKCM`, `UVM` |
| `GI` | `COAD`, `READ`, `ESCA`, `STAD` |
| `Hepatobiliary` | `LIHC`, `CHOL` |
| `Pancreas` | `PAAD` |
| `Kidney` | `KICH`, `KIRC`, `KIRP` |
| `Bladder` | `BLCA` |
| `Prostate` | `PRAD` |
| `Testis` | `TGCT` |
| `Adrenal` | `ACC`, `PCPG` |
| `Sarcoma` | `SARC` |
| `HemeLymph` | `LAML`, `DLBC` |
| `Mesothelioma` | `MESO` |

Cancer types not found in the mapping are assigned to:

```text
Other
```

## Tissue Type Mapping

The second script maps cancer types to broad tissue classes:

| Tissue type | Cancer types |
| --- | --- |
| `Adenocarcinoma` | `BRCA`, `COAD`, `READ`, `STAD`, `PAAD`, `PRAD`, `LUAD`, `OV`, `UCEC`, `UCS`, `CHOL`, `LIHC`, `THCA`, `KIRC`, `KIRP`, `KICH`, `ACC` |
| `Squamous` | `HNSC`, `LUSC` |
| `Mixed` | `ESCA`, `CESC`, `BLCA` |
| `NonEpithelial` | `GBM`, `LGG`, `LAML`, `DLBC`, `SKCM`, `UVM`, `SARC`, `TGCT`, `PCPG`, `MESO`, `THYM` |

Cancer types not found in the mapping are assigned to:

```text
Other
```

## How To Run From The Terminal

Replace the placeholders with the actual filenames you saved on your system. The commands are meant to be typed or pasted directly into the terminal.

## Important Assumptions

### Reproducibility
(/data/najemd/TF_binding_cancer/conda_envs/TF_binding_cancer_env) najemd@serv-bardet-dorsal:/data/najemd/TF_binding_cancer$ Rscript -e 'packages <- c("data.table", "Rtsne", "ggplot2", "Polychrome", "patchwork"); print(data.frame(package = packages, version = sapply(packages, function(pkg) as.character(packageVersion(pkg))), row.names = NULL))'

     package version
1 data.table  1.16.4
2      Rtsne    0.17
3    ggplot2   3.5.1
4 Polychrome   1.5.4
5  patchwork   1.3.0

The t-SNE seed is fixed:

```r
seed <- 1
```

This helps make the t-SNE embedding reproducible for the same input files, package versions, and environment.

## Example Outputs

After a successful run, expected outputs include:

```text
./results/peaks/tSNE_ATAC_ALLCANCERS_coords_10000.tsv
./results/peaks/tSNE_ATAC_ALLCANCERS_peakPresence_10000.pdf
./results/peaks/tSNE_ATAC_3views.pdf
```

## Interpretation

Samples that appear close together in the t-SNE plot have similar binary ATAC peak presence/absence profiles among the selected variable peaks. Clustering by cancer type, organ group, or tissue type suggests shared chromatin accessibility patterns across those biological categories.

t-SNE is primarily a visualization method. Distances between far-apart clusters and global geometry should be interpreted cautiously. For quantitative comparison between groups, use complementary analyses such as PCA, clustering, differential accessibility, or supervised classification.

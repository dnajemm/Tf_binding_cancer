 #!/usr/bin/env Rscript
"""This script runs a per-cancer PCA on TCGA HM450 methylation data using a fixed list of 10,000 CpGs, comparing tumor (01A) and normal (11A) samples and generating an interactive HTML PCA plot for each cancer type.
Input:
- Directory with filtered methylation data files for TCGA samples.
- A text file containing a list of 10,000 CpG sites (one per line).
- A text file listing cancer types that have both tumor and normal samples.
Output:
- For each cancer type with both tumor and normal samples, a TSV file with PCA coordinates and an interactive HTML PCA plot.
Output directory structure:
./results/PCA_methylation_by_cancer_interactive/TCGA_<CANCER>/
"""
suppressPackageStartupMessages({
  library(data.table)
  library(plotly)
  library(htmlwidgets)})

# ------------------ settings ------------------
meth_dir <- "./methylation/filtered_methylation"
top_cpgs_file <- "./methylation/top10000_variable_CpGs_pan_cancer.txt"   # 10k CpGs, one per line

out_root <- "./results/PCA_methylation_by_cancer_interactive"

# ------------------ load CpG list ------------------
top_cpgs <- unique(trimws(readLines(top_cpgs_file, warn = FALSE)))
top_cpgs <- top_cpgs[top_cpgs != ""]
k <- length(top_cpgs)
cat("Top CpGs:", k, "\n")

# ------------------ list all 01A/11A filtered files ------------------
files_all <- list.files(
  meth_dir,
  pattern = "_annotated_methylation_filtered\\.bed\\.gz$",
  full.names = TRUE
)
if (length(files_all) == 0) stop("No files found in: ", meth_dir)

bn_all <- basename(files_all)

# cancer code
cancer_all <- sub("^HM450_TCGA-([A-Z0-9]+)-.*$", "\\1", bn_all)

# group from filename
group_all <- ifelse(grepl("-01A_", bn_all), "Tumor",
             ifelse(grepl("-11A_", bn_all), "Normal", NA_character_))

meta_all <- data.table(
  file   = files_all,
  sample = bn_all,
  cancer = cancer_all,
  group  = group_all
)[!is.na(group)]  # keep only 01A/11A

if (nrow(meta_all) == 0) stop("No 01A/11A files matched naming pattern in: ", meth_dir)

# ------------------ cancers that have BOTH Tumor + Normal ------------------
cancers_file <- "./methylation/cancer_types_with_tumor_and_normal.txt"
cancers_with_both <- unique(trimws(readLines(cancers_file, warn = FALSE)))
cancers_with_both <- cancers_with_both[cancers_with_both != ""]
cat("Cancer types (from file):", length(cancers_with_both), "\n")

# ------------------ helper to read one sample on 10k CpGs ------------------
read_one_10k <- function(f, cpgs) {
  dt <- fread(cmd = paste("zcat", shQuote(f)),
              header = FALSE, select = c(4, 5), showProgress = FALSE)
  setnames(dt, c("CpG", "beta"))
  dt[, CpG := as.character(CpG)]
  dt[, beta := suppressWarnings(as.numeric(beta))]
  dt <- dt[CpG %chin% cpgs]
  if (nrow(dt) == 0) return(rep(NA_real_, length(cpgs)))
  setkey(dt, CpG)
  as.numeric(dt[.(cpgs), beta])  # in the same order as cpgs
}

# ------------------ loop cancers ------------------
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

for (cancer_keep in cancers_with_both) {
  cat("\n============================\n")
  cat("Processing cancer:", cancer_keep, "\n")

  meta <- meta_all[cancer == cancer_keep]
  cat("Samples:", nrow(meta), " (Tumor:", sum(meta$group == "Tumor"),
      ", Normal:", sum(meta$group == "Normal"), ")\n")

  # output folder per cancer
  out_dir <- file.path(out_root, paste0("TCGA_", cancer_keep))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  out_tsv  <- file.path(out_dir, paste0("PCA_HM450_", cancer_keep, "_01A_vs_11A.tsv"))
  out_html <- file.path(out_dir, paste0("PCA_HM450_", cancer_keep, "_01A_vs_11A.html"))

  # build CpGs x samples matrix
  X <- matrix(NA_real_, nrow = k, ncol = nrow(meta),
              dimnames = list(top_cpgs, meta$sample))

  for (i in seq_len(nrow(meta))) {
    if (i %% 25 == 0) cat("  Reading", i, "/", nrow(meta), "\n")
    X[, i] <- read_one_10k(meta$file[i], top_cpgs)
  }

  # PCA (samples x CpGs)
  M <- t(X)
  pca <- prcomp(M, center = TRUE, scale. = FALSE)

  ve <- (pca$sdev^2) / sum(pca$sdev^2)
  pc1_lab <- paste0("PC1 (", round(100 * ve[1], 1), "%)")
  pc2_lab <- paste0("PC2 (", round(100 * ve[2], 1), "%)")

  coords <- as.data.table(pca$x[, 1:2, drop = FALSE])
  setnames(coords, c("PC1", "PC2"))
  plot_dt <- cbind(meta[, .(sample, cancer, group)], coords)

  # save TSV
  fwrite(plot_dt, out_tsv, sep = "\t")

  # interactive plotly HTML
  p <- plot_ly(
    data = plot_dt,
    x = ~PC1, y = ~PC2,
    type = "scatter", mode = "markers",
    color = ~group,
    text = ~paste0("Sample: ", sample, "<br>Group: ", group, "<br>Cancer: ", cancer),
    hoverinfo = "text",
    marker = list(size = 6, opacity = 0.85)
  ) %>%
    layout(
      title = list(text = paste0("PCA — ", cancer_keep, " (01A vs 11A)")),
      xaxis = list(title = pc1_lab),
      yaxis = list(title = pc2_lab)
    )

  saveWidget(p, out_html, selfcontained = TRUE)

  cat("Wrote TSV :", out_tsv, "\n")
  cat("Wrote HTML:", out_html, "\n")
}

cat("\nDONE. Output root:\n", out_root, "\n", sep = "")

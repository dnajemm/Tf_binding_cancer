#################################################################################################################
# PCA of HM450 methylation data for tumor vs normal samples, separately by cancer type, one PDF per cancer type
#################################################################################################################
''' This script performs PCA on the methylation data for the 10,000 most variable CpGs across all samples, comparing tumor (01A) vs normal (11A) samples for each cancer type. It reads only the relevant CpGs from the filtered methylation files, fills missing values with the CpG mean, and then runs PCA to visualize the separation between tumor and normal samples in a faceted plot for each cancer type alone. The results are saved as a TSV file with PCA coordinates and a PDF with the PCA plot for all cancer type in one pdf for each cancer type.
#!/usr/bin/env Rscript

# ============================================================
# PURPOSE
# ============================================================
# Perform PCA of HM450 methylation data separately for each
# cancer type (01A tumor vs 11A normal), using the pan-cancer
# top 10,000 most variable CpGs.
#
# Cancer types are read from:
#   ./methylation/cancer_types_with_tumor_and_normal.txt
#
# For each cancer:
#   - Read only the 10k CpGs
#   - Impute missing values by CpG mean
#   - Run PCA (samples x CpGs)
#   - Save PC1/PC2 coordinates (TSV)
#   - Save PCA plot (PDF)
#
# Parallelization:
#   - Processes 3 cancer types in parallel
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(parallel)
})

# ------------------ settings ------------------
meth_dir <- "./methylation/filtered_methylation"
top_cpgs_file <- "./methylation/top10000_variable_CpGs_pan_cancer.txt"
cancers_file <- "./methylation/cancer_types_with_tumor_and_normal.txt"

# ------------------ load cancer list ------------------
cancers <- scan(cancers_file, what = "", quiet = TRUE)
cancers <- cancers[cancers != ""]
cat("Cancers to process:", paste(cancers, collapse = ", "), "\n")

# ------------------ load CpG list ------------------
top_cpgs <- unique(trimws(readLines(top_cpgs_file)))
top_cpgs <- top_cpgs[top_cpgs != ""]
k <- length(top_cpgs)
cat("Top CpGs:", k, "\n")

# ------------------ list all 01A/11A filtered files ------------------
files_all <- list.files(
  meth_dir,
  pattern = "^HM450_TCGA-[A-Z0-9]+-.*-(01A|11A)_.*_annotated_methylation_filtered\\.bed\\.gz$",
  full.names = TRUE
)

bn_all <- basename(files_all)
cancer_all <- sub("^HM450_TCGA-([A-Z0-9]+)-.*$", "\\1", bn_all)

group_all <- ifelse(grepl("-01A_", bn_all), "Tumor",
             ifelse(grepl("-11A_", bn_all), "Normal", NA_character_))
if (any(is.na(group_all))) stop("Some files are not 01A/11A")

meta_all <- data.table(
  file   = files_all,
  sample = bn_all,
  cancer = cancer_all,
  group  = group_all
)

# ============================================================
# FUNCTION: process one cancer type
# ============================================================
process_one_cancer <- function(cancer_keep) {

  cat("\n=============================\n")
  cat("Cancer:", cancer_keep, "\n")

  out_dir <- file.path("./results/methylation/PCA_methylation_one_cancer",
                       paste0("TCGA_", cancer_keep))
  out_pdf <- file.path(out_dir,
                       paste0("PCA_HM450_", cancer_keep, "_01A_vs_11A.pdf"))
  out_tsv <- file.path(out_dir,
                       paste0("PCA_HM450_", cancer_keep, "_01A_vs_11A.tsv"))

  meta <- meta_all[cancer == cancer_keep]
  if (nrow(meta) < 2) {
    cat("Skipping", cancer_keep, "(not enough samples)\n")
    return(NULL)
  }

  cat("Samples kept:", nrow(meta),
      "(Tumor:", sum(meta$group == "Tumor"),
      ", Normal:", sum(meta$group == "Normal"), ")\n")

  # ---------- read ONLY the 10k CpGs ----------
  X <- matrix(NA_real_, nrow = k, ncol = nrow(meta))
  rownames(X) <- top_cpgs
  colnames(X) <- meta$sample

  read_one_10k <- function(f, cpgs) {
    dt <- fread(cmd = paste("zcat", shQuote(f)),
                header = FALSE, select = c(4, 5), showProgress = FALSE)
    setnames(dt, c("CpG", "beta"))
    dt <- dt[CpG %chin% cpgs]
    v <- as.numeric(dt$beta)
    names(v) <- dt$CpG
    as.numeric(v[cpgs])
  }

  for (i in seq_len(nrow(meta))) {
    if (i %% 25 == 0)
      cat("[", cancer_keep, "] Reading", i, "/", nrow(meta), "\n")
    X[, i] <- read_one_10k(meta$file[i], top_cpgs)
  }

  # ---------- fill NAs by CpG mean ----------
  row_means <- rowMeans(X, na.rm = TRUE)
  na_rows <- which(rowSums(is.na(X)) > 0)
  for (r in na_rows) X[r, is.na(X[r, ])] <- row_means[r]

  # ---------- PCA ----------
  M <- t(X)
  pca <- prcomp(M, center = TRUE, scale. = FALSE)

  ve <- (pca$sdev^2) / sum(pca$sdev^2)
  pc1_lab <- paste0("PC1 (", round(100 * ve[1], 1), "%)")
  pc2_lab <- paste0("PC2 (", round(100 * ve[2], 1), "%)")

  coords <- as.data.table(pca$x[, 1:2, drop = FALSE])
  setnames(coords, c("PC1", "PC2"))

  plot_dt <- cbind(meta[, .(sample, cancer, group)], coords)

  # ---------- save outputs ----------
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(plot_dt, out_tsv, sep = "\t")

  p <- ggplot(plot_dt, aes(PC1, PC2, colour = group)) +
    geom_point(size = 1.2, alpha = 0.9) +
    coord_fixed() +
    theme_bw() +
    labs(
      title = paste0("PCA — HM450 filtered (", cancer_keep, ") (01A vs 11A)"),
      subtitle = paste0("CpGs used: ", ncol(M),
                        "; Samples: ", nrow(M)),
      x = pc1_lab, y = pc2_lab,
      colour = "Group"
    )

  ggsave(out_pdf, p, width = 8.5, height = 6.5)

  cat("Finished", cancer_keep, "\n")
  return(NULL)
}

# ============================================================
# RUN IN PARALLEL (3 cancers at a time)
# ============================================================
mclapply(
  cancers,
  process_one_cancer,
  mc.cores = 3
)

cat("\nALL DONE\n")
'''

###################################################################################################################
# PCA of HM450 methylation data for tumor vs normal samples all cancer type, one page faceted for all cancer type
###################################################################################################################
''' This script performs PCA on the methylation data for the 10,000 most variable CpGs across all samples, comparing tumor (01A) vs normal (11A) samples for each cancer type. It reads only the relevant CpGs from the filtered methylation files, fills missing values with the CpG mean, and then runs PCA to visualize the separation between tumor and normal samples in a faceted plot by cancer type. The results are saved as a TSV file with PCA coordinates and a PDF with the PCA plot for all cancer type in one fce of the pdf.
#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)

# ------------------ settings ------------------
meth_dir <- "./methylation/filtered_methylation"
top_cpgs_file <- file.path("./methylation/top10000_variable_CpGs_pan_cancer.txt")
  
out_dir <- "./results/methylation/PCA_methylation_pan_cancer"
out_pdf <- file.path(out_dir, "PCA_HM450_pan_cancer_01A_vs_11A_faceted.pdf")
out_tsv <- file.path(out_dir, "PCA_HM450_pan_cancer_01A_vs_11A_faceted.tsv")

# ------------------ load CpG list ------------------
top_cpgs <- unique(trimws(readLines(top_cpgs_file)))
top_cpgs <- top_cpgs[top_cpgs != ""]
k <- length(top_cpgs)
cat("Top CpGs:", k, "\n")

# ------------------ list all 01A/11A filtered files ------------------
files <- list.files(
  meth_dir,
  pattern = "^HM450_TCGA-[A-Z0-9]+-.*-(01A|11A)_.*_annotated_methylation_filtered\\.bed\\.gz$",
  full.names = TRUE
)

bn <- basename(files)
cancer <- sub("^HM450_TCGA-([A-Z0-9]+)-.*$", "\\1", bn)

group <- ifelse(grepl("-01A_", bn), "Tumor",
         ifelse(grepl("-11A_", bn), "Normal", NA_character_))
if (any(is.na(group))) stop("Some files are not 01A/11A:\n", paste(bn[is.na(group)], collapse = "\n"))

meta <- data.table(
  file   = files,
  sample = bn,
  cancer = cancer,
  group  = group
)

cat("Total samples:", nrow(meta), "\n")

# ------------------ read ONLY the 10k CpGs into a matrix (CpGs x samples) ------------------
X <- matrix(NA_real_, nrow = k, ncol = nrow(meta))
rownames(X) <- top_cpgs
colnames(X) <- meta$sample

read_one_10k <- function(f, cpgs) {
  dt <- fread(cmd = paste("zcat", shQuote(f)),
              header = FALSE, select = c(4, 5), showProgress = FALSE)
  setnames(dt, c("CpG", "beta"))
  dt <- dt[CpG %chin% cpgs]
  v <- as.numeric(dt$beta)
  names(v) <- dt$CpG
  as.numeric(v[cpgs])  # returns in same order as cpgs (NAs if missing)
}

for (i in seq_len(nrow(meta))) {
  if (i %% 25 == 0) cat("Reading", i, "/", nrow(meta), "\n")
  X[, i] <- read_one_10k(meta$file[i], top_cpgs)
}

# ------------------ fill remaining NAs by CpG mean (so PCA can run) ------------------
row_means <- rowMeans(X, na.rm = TRUE)
na_rows <- which(rowSums(is.na(X)) > 0)
for (r in na_rows) X[r, is.na(X[r, ])] <- row_means[r]

# ------------------ PCA (samples x CpGs) ------------------
M <- t(X)  # samples x CpGs
pca <- prcomp(M, center = TRUE, scale. = FALSE)

ve <- (pca$sdev^2) / sum(pca$sdev^2)
pc1_lab <- paste0("PC1 (", round(100 * ve[1], 1), "%)")
pc2_lab <- paste0("PC2 (", round(100 * ve[2], 1), "%)")

coords <- as.data.table(pca$x[, 1:2, drop = FALSE])
setnames(coords, c("PC1", "PC2"))

plot_dt <- cbind(meta[, .(sample, cancer, group)], coords)

# ------------------ save outputs ------------------
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(plot_dt, out_tsv, sep = "\t")

p <- ggplot(plot_dt, aes(PC1, PC2, colour = group)) +
  geom_point(size = 0.9, alpha = 0.9) +
  facet_wrap(~ cancer) +
  coord_fixed() +
  theme_bw() +
  labs(
    title = "Pan-cancer PCA — HM450 filtered (01A vs 11A)",
    subtitle = paste0("CpGs used: ", ncol(M), "; Samples: ", nrow(M)),
    x = pc1_lab, y = pc2_lab,
    colour = "Group"
  )

ggsave(out_pdf, p, width = 20, height = 30)

cat("Wrote TSV:", out_tsv, "\n")
cat("Wrote PDF:", out_pdf, "\n")
'''
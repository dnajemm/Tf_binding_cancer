#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)

# ------------------ settings ------------------
cancer_keep <- "LUAD"   # <<< change this to the cancer type 

meth_dir <- "./methylation/filtered_methylation"
top_cpgs_file <- file.path("cpg_list.txt")  # use the common CpG list

out_dir <- file.path("./results/PCA_methylation_by_cancer", paste0("TCGA_", cancer_keep))
out_pdf <- file.path(out_dir, paste0("PCA_HM450_", cancer_keep, "_01A_vs_11A.pdf"))
out_tsv <- file.path(out_dir, paste0("PCA_HM450_", cancer_keep, "_01A_vs_11A.tsv"))

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
if (any(is.na(group_all))) stop("Some files are not 01A/11A:\n", paste(bn_all[is.na(group_all)], collapse = "\n"))

meta_all <- data.table(
  file   = files_all,
  sample = bn_all,
  cancer = cancer_all,
  group  = group_all
)

# ------------------ keep ONLY one cancer type ------------------
meta <- meta_all[cancer == cancer_keep]
if (nrow(meta) < 2) stop("Not enough samples for cancer type: ", cancer_keep)

cat("Cancer:", cancer_keep, "\n")
cat("Samples kept:", nrow(meta), " (Tumor:", sum(meta$group == "Tumor"),
    ", Normal:", sum(meta$group == "Normal"), ")\n")

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
  as.numeric(v[cpgs])  # ordered like cpgs (NAs if missing)
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
  geom_point(size = 1.2, alpha = 0.9) +
  coord_fixed() +
  theme_bw() +
  labs(
    title = paste0("PCA — HM450 filtered (", cancer_keep, ") (01A vs 11A)"),
    subtitle = paste0("CpGs used: ", ncol(M), "; Samples: ", nrow(M)),
    x = pc1_lab, y = pc2_lab,
    colour = "Group"
  )

ggsave(out_pdf, p, width = 8.5, height = 6.5)

cat("Wrote TSV:", out_tsv, "\n")
cat("Wrote PDF:", out_pdf, "\n")

#!/usr/bin/env Rscript
# t-SNE of TCGA tumor methylation (HM450) — 01A only
# Uses a PROVIDED CpG list (e.g. 10,000 CpGs).
# Colors points by epithelial class. Saves PDF + coordinates TSV/RDS.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(Rtsne)
})

# ------------------ Inputs / outputs ------------------
infile      <- "./methylation/filtered_methylation_files_tumor01A_1.txt"
top10k_file <- "./methylation/top10000_variable_CpGs_pan_cancer.txt"
pan_mat_rds <- "./methylation/TCGA-all_cancer_pan_cancer_matrix.RData"

outpdf <- "./results/summarys/tSNE_methylation_01A_epithelialColors_top10k.pdf"
out_tsv <- "./results/summarys/tSNE_methylation_01A_epithelialCoords_top10k.tsv"
out_rds <- "./results/summarys/tSNE_methylation_01A_epithelialCoords_top10k.rds"

# ------------------ Parameters ------------------
n_pcs         <- 30
max_iter_tsne <- 1000
seed          <- 1

# ------------------ Load CpG list ------------------
top_probes <- unique(trimws(readLines(top10k_file, warn = FALSE)))
top_probes <- top_probes[top_probes != ""]
cat(">>> Loaded CpGs from file:", length(top_probes), "\n")
if (length(top_probes) < 100) stop("CpG list looks too small")

# ------------------ Load tumor file list (01A only) ------------------
files_all <- scan(infile, what = "", quiet = TRUE)
files_all <- files_all[file.exists(files_all)]
if (length(files_all) < 2) stop("Not enough files")

sample_names_all <- basename(files_all)
cancer_all <- sapply(strsplit(sample_names_all, "-"), `[`, 2)

# ------------------ Load pan-cancer matrix ------------------
cat(">>> Loading pan-cancer matrix...\n")
load(pan_mat_rds)   # loads pan_mat (CpGs × samples)
if (!exists("pan_mat")) stop("pan_mat not found")

# ------------------ Subset matrix: CpGs + tumor samples ------------------
common_cpgs <- intersect(top_probes, rownames(pan_mat))
cat(">>> CpGs retained:", length(common_cpgs), "\n")
if (length(common_cpgs) < 100) stop("Too few CpGs after intersection")

missing_samples <- setdiff(sample_names_all, colnames(pan_mat))
if (length(missing_samples) > 0) {
  stop("Missing samples in pan_mat:\n", paste(missing_samples, collapse = "\n"))
}

# samples × CpGs
X <- t(pan_mat[common_cpgs, sample_names_all, drop = FALSE])
rm(pan_mat); gc()

# ------------------ Impute missing + scale ------------------
cat(">>> Imputing missing values + scaling...\n")

probe_means <- colMeans(X, na.rm = TRUE)
for (j in seq_len(ncol(X))) {
  miss <- is.na(X[, j])
  if (any(miss)) X[miss, j] <- probe_means[j]
}

sds <- apply(X, 2, sd, na.rm = TRUE)
keep <- is.finite(sds) & sds > 0
cat(">>> Removing constant CpGs:", sum(!keep), "\n")
X <- X[, keep, drop = FALSE]

Xs <- scale(X)

# ------------------ PCA + t-SNE ------------------
cat(">>> Running PCA...\n")
pca <- prcomp(Xs, center = FALSE, scale. = FALSE)
X_pca <- pca$x[, seq_len(min(n_pcs, ncol(pca$x))), drop = FALSE]

n <- nrow(X_pca)
perp_max <- floor((n - 1) / 3) - 1
if (perp_max < 5) stop("Too few samples for t-SNE")
perp <- max(5, min(30, perp_max))

cat(">>> Running t-SNE (perplexity =", perp, ")...\n")
set.seed(seed)
tsne <- Rtsne(
  X_pca,
  perplexity = perp,
  check_duplicates = FALSE,
  max_iter = max_iter_tsne,
  verbose = TRUE,
  pca = FALSE
)

tsne_df <- data.frame(
  sample = rownames(X),
  cancer = factor(cancer_all),
  tSNE1  = tsne$Y[, 1],
  tSNE2  = tsne$Y[, 2]
)

# ------------------ Epithelial class mapping ------------------
epi_map <- list(
  Adenocarcinoma = c("BRCA","COAD","READ","STAD","PAAD","PRAD","LUAD","OV","UCEC","UCS",
                     "CHOL","LIHC","THCA","KIRC","KIRP","KICH","ACC"),
  Squamous = c("HNSC","LUSC"),
  Mixed = c("ESCA","CESC","BLCA"),
  NonEpithelial = c("GBM","LGG","LAML","DLBC","SKCM","UVM","SARC",
                    "TGCT","PCPG","MESO","THYM")
)

cancer_to_epi <- function(ct) {
  for (g in names(epi_map)) if (ct %in% epi_map[[g]]) return(g)
  "Other"
}

tsne_df$epithelial_class <- factor(
  vapply(as.character(tsne_df$cancer), cancer_to_epi, character(1)),
  levels = c("Adenocarcinoma","Squamous","Mixed","NonEpithelial","Other")
)

# ------------------ Save coordinates ------------------
dir.create(dirname(outpdf), recursive = TRUE, showWarnings = FALSE)
fwrite(as.data.table(tsne_df), out_tsv, sep = "\t")
saveRDS(tsne_df, out_rds)

# ------------------ Colors ------------------
epi_pal <- c(
  Adenocarcinoma = "#1b9e77",
  Squamous       = "#d95f02",
  Mixed          = "#7570b3",
  NonEpithelial  = "#e7298a",
  Other          = "grey40"
)

# ------------------ Plot ------------------
pdf(outpdf, width = 10.5, height = 8.5)

g <- ggplot(tsne_df, aes(tSNE1, tSNE2, color = epithelial_class)) +
  geom_point(size = 1.7, alpha = 0.85) +
  scale_color_manual(values = epi_pal, drop = FALSE) +
  theme_minimal() +
  labs(
    title = "t-SNE of TCGA tumor methylation (HM450) — 01A only",
    subtitle = paste0(
      "CpGs used: ", ncol(X),
      "; PCA dims: ", ncol(X_pca),
      "; perplexity=", perp
    ),
    color = "Epithelial class",
    x = "tSNE1", y = "tSNE2"
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme(legend.position = "right")

print(g)
dev.off()

cat("Wrote:", outpdf, "\n")
cat("Wrote:", out_tsv, "\n")
cat("Wrote:", out_rds, "\n")
cat(">>> DONE\n")

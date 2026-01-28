#!/usr/bin/env Rscript
# t-SNE of TCGA tumor methylation (HM450) — 01A only
# Uses a PROVIDED CpG list (e.g. 10,000 CpGs).
# Colors by cancer type; shapes by organ group; uses Polychrome if installed.
# Saves PDF + coordinates TSV/RDS.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(Rtsne)
})

has_polychrome <- requireNamespace("Polychrome", quietly = TRUE)

# ------------------ Inputs / outputs ------------------
infile      <- "./methylation/filtered_methylation_files_tumor11A_1.txt"
top10k_file <- "./methylation/top10000_variable_CpGs_pan_cancer.txt"
pan_mat_rds <- "./methylation/TCGA-all_cancer_pan_cancer_matrix.RData"

outpdf <- "./results/summarys/tSNE_methylation_11A_top10k.pdf"
out_tsv <- "./results/summarys/tSNE_methylation_11A_top10k_coords.tsv"
out_rds <- "./results/summarys/tSNE_methylation_11A_top10k_coords.rds"

# ------------------ Parameters ------------------
n_pcs         <- 30
max_iter_tsne <- 1000
seed          <- 1

# ------------------ Load CpG list ------------------
top_probes <- unique(trimws(readLines(top10k_file, warn = FALSE)))
top_probes <- top_probes[top_probes != ""]
cat(">>> Loaded CpGs from file:", length(top_probes), "\n")

# ------------------ Load tumor file list (01A only) ------------------
files_all <- scan(infile, what = "", quiet = TRUE)
files_all <- files_all[file.exists(files_all)]
cat(">>> Total existing tumor files:", length(files_all), "\n")
if (length(files_all) < 2) stop("Not enough files after filtering.")

sample_names_all <- basename(files_all)
cancer_all <- sapply(strsplit(sample_names_all, "-"), `[`, 2)

# ------------------ Load pan-cancer matrix ------------------
cat(">>> Loading pan-cancer methylation matrix...\n")
load(pan_mat_rds)   # loads object: pan_mat  (CpGs x samples)

if (!exists("pan_mat")) stop("pan_mat not found in RData")

# ------------------ Subset matrix (01A samples + top CpGs) ------------------
common_cpgs <- intersect(top_probes, rownames(pan_mat))
cat(">>> CpGs retained after intersection:", length(common_cpgs), "\n")
if (length(common_cpgs) < 100) stop("Too few CpGs after intersection")

missing_samples <- setdiff(sample_names_all, colnames(pan_mat))
if (length(missing_samples) > 0) {
  stop("These samples are missing from pan_mat:\n",
       paste(missing_samples, collapse = "\n"))
}

# samples × CpGs
X <- t(pan_mat[common_cpgs, sample_names_all, drop = FALSE])

rm(pan_mat); gc()

# ------------------ Impute missing + remove constant CpGs + scale ------------------
cat(">>> Imputing missing values (per CpG mean) + filtering constants + scaling...\n")

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
perp <- max(5, min(30, floor((n - 1) / 3) - 1))
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

# ------------------ Organ grouping ------------------
organ_map <- list(
  "CNS"           = c("GBM","LGG"),
  "Breast"        = c("BRCA"),
  "Gynecologic"   = c("CESC","OV","UCEC","UCS"),
  "Lung"          = c("LUAD","LUSC"),
  "HeadNeck"      = c("HNSC"),
  "ThyroidThymus" = c("THCA","THYM"),
  "SkinEye"       = c("SKCM","UVM"),
  "GI"            = c("COAD","READ","ESCA","STAD"),
  "Hepatobiliary" = c("LIHC","CHOL"),
  "Pancreas"      = c("PAAD"),
  "Kidney"        = c("KICH","KIRC","KIRP"),
  "Bladder"       = c("BLCA"),
  "Prostate"      = c("PRAD"),
  "Testis"        = c("TGCT"),
  "Adrenal"       = c("ACC","PCPG"),
  "Sarcoma"       = c("SARC"),
  "HemeLymph"     = c("LAML","DLBC"),
  "Mesothelioma"  = c("MESO")
)

cancer_to_organ <- function(ct) {
  for (g in names(organ_map)) if (ct %in% organ_map[[g]]) return(g)
  "Other"
}
tsne_df$organ_group <- factor(
  vapply(as.character(tsne_df$cancer), cancer_to_organ, character(1))
)

shape_values <- c(16, 17, 15, 3, 7, 8, 18, 4, 1, 2, 5, 6, 9, 10, 11, 12, 13, 14, 0)
names(shape_values) <- levels(tsne_df$organ_group)

# ------------------ Save coordinates ------------------
dir.create(dirname(outpdf), recursive = TRUE, showWarnings = FALSE)
fwrite(as.data.table(tsne_df), out_tsv, sep = "\t")
saveRDS(tsne_df, out_rds)

# ------------------ Colors ------------------
pal <- NULL
if (has_polychrome) {
  n_cancers <- nlevels(tsne_df$cancer)
  pal <- Polychrome::createPalette(
    n_cancers,
    seedcolors = c("#000000", "#E69F00", "#56B4E9", "#f93a00ff")
  )
  names(pal) <- levels(tsne_df$cancer)
}

# ------------------ Plot ------------------
pdf(outpdf, width = 11, height = 9)

g <- ggplot(tsne_df, aes(tSNE1, tSNE2, color = cancer, shape = organ_group)) +
  geom_point(size = 1.6, alpha = 0.85) +
  scale_shape_manual(values = shape_values) +
  theme_minimal() +
  labs(
    title = "t-SNE of TCGA Healthy methylation (HM450) — 11A only",
    subtitle = paste0(
      "CpGs used: ", ncol(X),
      "; PCA dims: ", ncol(X_pca),
      "; perplexity=", perp
    ),
    color = "Cancer type",
    shape = "Organ group",
    x = "tSNE1", y = "tSNE2"
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme(legend.position = "right")

if (!is.null(pal)) {
  g <- g + scale_color_manual(values = pal)
} else {
  g <- g + scale_color_hue(h = c(0, 360), l = 60, c = 120)
}

print(g)
dev.off()

cat("Wrote:", outpdf, "\n")
cat("Wrote:", out_tsv, "\n")
cat("Wrote:", out_rds, "\n")
cat(">>> DONE\n")



# tSNE for both tumor and healthy samples:
#!/usr/bin/env Rscript
# t-SNE of TCGA methylation (HM450) — paired Tumor (01A) + Normal (11A)
# Uses provided CpG list (top 10k).
# Colors by cancer; shapes by sample type (Tumor/Normal).
# Saves PDF + coordinates TSV/RDS.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(Rtsne)
})

has_polychrome <- requireNamespace("Polychrome", quietly = TRUE)

# ------------------ Inputs / outputs ------------------
pairs_tsv   <- "./methylation/sample_pairs_files/methylation_pairs_filtered.tsv"
top10k_file <- "./methylation/top10000_variable_CpGs_pan_cancer.txt"
pan_mat_rds <- "./methylation/TCGA-all_cancer_pan_cancer_matrix.RData"

outpdf  <- "./results/summarys/tSNE_methylation_01A_11A_top10k.pdf"
out_tsv <- "./results/summarys/tSNE_methylation_01A_11A_top10k_coords.tsv"
out_rds <- "./results/summarys/tSNE_methylation_01A_11A_top10k_coords.rds"

# ------------------ Parameters ------------------
n_pcs         <- 30
max_iter_tsne <- 1000
seed          <- 1

# ------------------ Load CpG list ------------------
top_probes <- unique(trimws(readLines(top10k_file, warn = FALSE)))
top_probes <- top_probes[top_probes != ""]
cat(">>> Loaded CpGs from file:", length(top_probes), "\n")

# ------------------ Load pairs ------------------
pairs <- fread(pairs_tsv)
need_cols <- c("cancer","patient","tumor_file","healthy_file")
stopifnot(all(need_cols %in% names(pairs)))

# keep only existing files
pairs <- pairs[file.exists(tumor_file) & file.exists(healthy_file)]
cat(">>> Pairs with existing tumor+healthy files:", nrow(pairs), "\n")
if (nrow(pairs) < 2) stop("Not enough pairs after filtering existing files.")

# Build sample table: one row per sample
tumor_dt <- pairs[, .(sample = basename(tumor_file),
                      file   = tumor_file,
                      cancer = sub("^TCGA-", "", cancer),  # optional cleanup
                      group  = "Tumor")]
healthy_dt <- pairs[, .(sample = basename(healthy_file),
                        file   = healthy_file,
                        cancer = sub("^TCGA-", "", cancer),
                        group  = "Normal")]

meta <- rbind(tumor_dt, healthy_dt, use.names = TRUE)
meta[, group := factor(group, levels = c("Tumor","Normal"))]
meta[, cancer := factor(cancer)]

cat(">>> Total samples (Tumor+Normal):", nrow(meta), "\n")

# ------------------ Load pan-cancer matrix ------------------
cat(">>> Loading pan-cancer methylation matrix...\n")
load(pan_mat_rds)   # expects object: pan_mat (CpGs x samples)
if (!exists("pan_mat")) stop("pan_mat not found in RData")

# ------------------ Subset matrix (top CpGs + selected samples) ------------------
common_cpgs <- intersect(top_probes, rownames(pan_mat))
cat(">>> CpGs retained after intersection:", length(common_cpgs), "\n")
if (length(common_cpgs) < 100) stop("Too few CpGs after intersection")

missing_samples <- setdiff(meta$sample, colnames(pan_mat))
if (length(missing_samples) > 0) {
  stop("These samples are missing from pan_mat:\n",
       paste(missing_samples, collapse = "\n"))
}

# samples × CpGs
X <- t(pan_mat[common_cpgs, meta$sample, drop = FALSE])

rm(pan_mat); gc()

# ------------------ Impute missing + remove constant CpGs + scale ------------------
cat(">>> Imputing missing values (per CpG mean) + filtering constants + scaling...\n")

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
perp <- max(5, min(30, floor((n - 1) / 3) - 1))
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

# ------------------ Output table ------------------
coords <- data.table(
  sample = meta$sample,
  cancer = meta$cancer,
  group  = meta$group,
  tSNE1  = tsne$Y[, 1],
  tSNE2  = tsne$Y[, 2]
)

# ------------------ Save coordinates ------------------
dir.create(dirname(outpdf), recursive = TRUE, showWarnings = FALSE)
fwrite(coords, out_tsv, sep = "\t")
saveRDS(coords, out_rds)

# ------------------ Colors ------------------
pal <- NULL
if (has_polychrome) {
  n_cancers <- nlevels(coords$cancer)
  pal <- Polychrome::createPalette(
    n_cancers,
    seedcolors = c("#000000", "#E69F00", "#56B4E9", "#f93a00ff")
  )
  names(pal) <- levels(coords$cancer)
}

# ------------------ Plot ------------------
pdf(outpdf, width = 11, height = 9)

g <- ggplot(coords, aes(tSNE1, tSNE2, color = cancer, shape = group)) +
  geom_point(size = 1.6, alpha = 0.85) +
  scale_shape_manual(values = c(Tumor = 16, Normal = 17)) +
  theme_minimal() +
  labs(
    title = "t-SNE of TCGA methylation (HM450) — Tumor (01A) + Normal (11A) pairs",
    subtitle = paste0(
      "CpGs used: ", ncol(X),
      "; PCA dims: ", ncol(X_pca),
      "; perplexity=", perp,
      "; samples=", nrow(coords)
    ),
    color = "Cancer type",
    shape = "Sample type",
    x = "tSNE1", y = "tSNE2"
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme(legend.position = "right")

if (!is.null(pal)) {
  g <- g + scale_color_manual(values = pal)
} else {
  g <- g + scale_color_hue(h = c(0, 360), l = 60, c = 120)
}

print(g)
dev.off()

cat("Wrote:", outpdf, "\n")
cat("Wrote:", out_tsv, "\n")
cat("Wrote:", out_rds, "\n")
cat(">>> DONE\n")

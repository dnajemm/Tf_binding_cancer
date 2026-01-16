#!/usr/bin/env Rscript
# t-SNE of TCGA tumor methylation (HM450) — 01A only
# Uses a PROVIDED CpG list (e.g. 10,000 CpGs from cpg_list.txt).
# Colors by cancer type; shapes by organ group; uses Polychrome if installed.
# Saves PDF + coordinates TSV/RDS.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(Rtsne)
})

has_polychrome <- requireNamespace("Polychrome", quietly = TRUE)

# ------------------ Inputs / outputs ------------------
infile      <- "./methylation/filtered_methylation_files_tumor01A_1.txt"
top10k_file <- "./cpg_list.txt"

outpdf <- "./results/summarys/tSNE_methylation_01A_top10kFromFile_polychrome.pdf"
out_tsv <- "./results/summarys/tSNE_methylation_01A_top10kFromFile_coords.tsv"
out_rds <- "./results/summarys/tSNE_methylation_01A_top10kFromFile_coords.rds"

# ------------------ Parameters ------------------
n_pcs         <- 30
max_iter_tsne <- 1000
seed          <- 1

# ------------------ Load CpG list ------------------
top_probes <- unique(trimws(readLines(top10k_file, warn = FALSE)))
top_probes <- top_probes[top_probes != ""]
cat(">>> Loaded CpGs from file:", length(top_probes), "\n")
if (length(top_probes) < 100) stop("CpG list looks too small: ", length(top_probes))

# ------------------ Load file list ------------------
files_all <- scan(infile, what = "", quiet = TRUE)
files_all <- files_all[file.exists(files_all)]
cat(">>> Total existing files:", length(files_all), "\n")
if (length(files_all) < 2) stop("Not enough files after filtering.")

sample_names_all <- basename(files_all)
cancer_all <- sapply(strsplit(sample_names_all, "-"), `[`, 2)

# ------------------ FILTER: keep ONLY 01A tumors ------------------
get_sample_type <- function(name) {
  m <- regexpr("-[0-9][0-9]A_", name)
  if (m[1] == -1) return(NA_character_)
  substr(name, m[1] + 1, m[1] + 3)  # "01A"
}
sample_type <- vapply(sample_names_all, get_sample_type, character(1))
keep_01A <- !is.na(sample_type) & sample_type == "01A"

files_all <- files_all[keep_01A]
sample_names_all <- sample_names_all[keep_01A]
cancer_all <- cancer_all[keep_01A]

cat(">>> After filtering (01A only):", length(files_all), "files\n")

# ------------------ Build matrix for all samples using provided CpGs ------------------
cat(">>> Building samples × CpGs matrix using provided list...\n")

n_samples <- length(files_all)
p <- length(top_probes)

X <- matrix(NA_real_, nrow = n_samples, ncol = p,
            dimnames = list(sample_names_all, top_probes))

t0 <- proc.time()[3]
for (i in seq_along(files_all)) {
  dt <- fread(cmd = paste("zcat", shQuote(files_all[i])),
              header = FALSE, select = c(4, 5), showProgress = FALSE)
  setnames(dt, c("CpG", "beta"))
  dt[, CpG := as.character(CpG)]
  dt[, beta := suppressWarnings(as.numeric(beta))]

  dt <- dt[CpG %chin% top_probes]
  if (nrow(dt) > 0) {
    setkey(dt, CpG)
    X[i, ] <- dt[.(top_probes), beta]   # exact CpG order
  }

  if (i %% 100 == 0) {
    cat(sprintf("   - Filled %d / %d samples (%.1f min)\n",
                i, n_samples, (proc.time()[3] - t0) / 60))
  }
}

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

# ------------------ Organ grouping for shapes ------------------
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
tsne_df$organ_group <- factor(vapply(as.character(tsne_df$cancer), cancer_to_organ, character(1)))

shape_values <- c(16, 17, 15, 3, 7, 8, 18, 4, 1, 2, 5, 6, 9, 10, 11, 12, 13, 14, 0)
names(shape_values) <- levels(tsne_df$organ_group)

# ------------------ Save coordinates ------------------
dir.create(dirname(outpdf), recursive = TRUE, showWarnings = FALSE)
fwrite(as.data.table(tsne_df), out_tsv, sep = "\t")
saveRDS(tsne_df, out_rds)

# ------------------ Colors (Polychrome if available) ------------------
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
    title = "t-SNE of TCGA tumor methylation (HM450) — 01A only",
    subtitle = paste0("CpGs used: ", ncol(X), "; PCA dims: ", min(n_pcs, ncol(X_pca)),
                      "; perplexity=", perp),
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
  message("NOTE: install Polychrome for more distinct colors: install.packages(\"Polychrome\")")
}

print(g)
dev.off()

cat("Wrote:", outpdf, "\n")
cat("Wrote:", out_tsv, "\n")
cat("Wrote:", out_rds, "\n")
cat(">>> DONE\n")

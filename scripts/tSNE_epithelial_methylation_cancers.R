#!/usr/bin/env Rscript
# t-SNE of TCGA tumor methylation (HM450) — 01A only
# Color points by epithelial class (Adenocarcinoma / Squamous / Mixed / NonEpithelial / Other)
# Shape is a single constant (no organ grouping)

library(data.table)
library(ggplot2)
library(Rtsne)

has_polychrome <- requireNamespace("Polychrome", quietly = TRUE)

# ------------------ Inputs / outputs ------------------
infile <- "./methylation/methylation_files_tumor0x.txt"
outpdf <- "./results/summarys/tSNE_methylation_01A_top2000SD_epithelialColors.pdf"
out_tsv <- "./results/summarys/tSNE_methylation_01A_top2000SD_epithelialCoords.tsv"
out_rds <- "./results/summarys/tSNE_methylation_01A_top2000SD_epithelialCoords.rds"

# ------------------ Parameters ------------------
max_per_cancer <- 200      # only for selecting variable probes
n_var_probes   <- 2000     # number of most-variable CpGs
n_pcs          <- 30       # PCA dims before tSNE
max_iter_tsne  <- 1000
seed <- 1

# ------------------ Load file list ------------------
files_all <- scan(infile, what = "", quiet = TRUE)
files_all <- files_all[file.exists(files_all)]
cat(">>> Total existing files:", length(files_all), "\n")
if (length(files_all) < 2) stop("Not enough files after filtering.")

sample_names_all <- basename(files_all)
cancer_all <- sapply(strsplit(sample_names_all, "-"), `[`, 2)

# ------------------ FILTER: keep ONLY 01A tumors ------------------
# filenames look like: ...-01A_1_annotated_methylation.bed.gz
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
if (length(files_all) < 2) stop("Not enough 01A files after filtering.")

# ------------------ Select informative samples ------------------
cat(">>> Scoring samples by within-sample beta variance...\n")

sample_var <- vapply(seq_along(files_all), function(i) {
  f <- files_all[i]
  dt <- tryCatch(
    fread(cmd = paste("zcat", shQuote(f)), header = FALSE, select = 5, showProgress = FALSE),
    error = function(e) NULL
  )
  if (is.null(dt) || ncol(dt) < 1) return(NA_real_)
  x <- suppressWarnings(as.numeric(dt[[1]]))
  stats::var(x, na.rm = TRUE)
}, numeric(1))

keep_ok <- !is.na(sample_var)
files_all <- files_all[keep_ok]
sample_names_all <- sample_names_all[keep_ok]
cancer_all <- cancer_all[keep_ok]
sample_var <- sample_var[keep_ok]
cat(">>> Files kept:", length(files_all), "\n")

idx_by_cancer <- split(seq_along(files_all), cancer_all)
subset_idx <- unlist(lapply(idx_by_cancer, function(idx) {
  idx[order(sample_var[idx], decreasing = TRUE)][seq_len(min(max_per_cancer, length(idx)))]
}), use.names = FALSE)

files_sub <- files_all[subset_idx]
cat(">>> Subset for probe selection:", length(files_sub), "files\n")

# ------------------ PASS 1: select top variable probes ------------------
cat(">>> PASS 1: Selecting top", n_var_probes, "variable probes...\n")

dt0 <- fread(cmd = paste("zcat", shQuote(files_sub[1])), header = FALSE, showProgress = FALSE)
if (ncol(dt0) < 5) stop("Bad format:", files_sub[1])

probe_ids <- as.character(dt0[[4]])
x0 <- suppressWarnings(as.numeric(dt0[[5]]))

n <- 1L
mean_vec <- x0
M2_vec <- rep(0, length(x0))

if (length(files_sub) > 1) {
  for (i in 2:length(files_sub)) {
    dt <- fread(cmd = paste("zcat", shQuote(files_sub[i])), header = FALSE, showProgress = FALSE)
    if (ncol(dt) < 5) next

    if (!identical(as.character(dt[[4]]), probe_ids)) {
      dt <- dt[, .(probe = as.character(V4), beta = suppressWarnings(as.numeric(V5)))]
      setkey(dt, probe)
      x <- dt[.(probe_ids), beta]
    } else {
      x <- suppressWarnings(as.numeric(dt[[5]]))
    }

    n_new <- n + 1L
    delta <- x - mean_vec
    mean_vec <- mean_vec + delta / n_new
    delta2 <- x - mean_vec
    M2_vec <- M2_vec + delta * delta2
    n <- n_new

    if (i %% 50 == 0) cat("   - Processed", i, "/", length(files_sub), "\n")
  }
}

sd_vec <- sqrt(M2_vec / pmax(n - 1L, 1L))
sd_vec[!is.finite(sd_vec)] <- 0
top_idx <- order(sd_vec, decreasing = TRUE)[seq_len(min(n_var_probes, length(sd_vec)))]
top_probes <- probe_ids[top_idx]

cat(">>> Selected", length(top_probes), "probes\n")

# ------------------ PASS 2: build matrix for all samples ------------------
cat(">>> PASS 2: Building samples × probes matrix...\n")

n_samples <- length(files_all)
p <- length(top_probes)

X <- matrix(NA_real_, nrow = n_samples, ncol = p,
            dimnames = list(sample_names_all, top_probes))

t0 <- proc.time()[3]
for (i in seq_along(files_all)) {
  dt <- fread(cmd = paste("zcat", shQuote(files_all[i])), header = FALSE, showProgress = FALSE)
  if (ncol(dt) < 5) next

  probes <- as.character(dt[[4]])
  betas  <- suppressWarnings(as.numeric(dt[[5]]))

  m <- match(probes, top_probes)
  ok <- !is.na(m)
  if (any(ok)) X[i, m[ok]] <- betas[ok]

  if (i %% 100 == 0)
    cat(sprintf("   - Filled %d / %d samples (%.1f min)\n",
                i, n_samples, (proc.time()[3] - t0) / 60))
}

# ------------------ Clean / impute / scale ------------------
na_rate <- colMeans(is.na(X))
X <- X[, na_rate < 0.2, drop = FALSE]
cat(">>> Probes kept after NA filter:", ncol(X), "\n")

probe_means <- colMeans(X, na.rm = TRUE)
for (j in seq_len(ncol(X))) {
  miss <- is.na(X[, j])
  if (any(miss)) X[miss, j] <- probe_means[j]
}

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
  verbose = TRUE
)

tsne_df <- data.frame(
  sample = rownames(X),
  cancer = factor(cancer_all),
  tSNE1 = tsne$Y[, 1],
  tSNE2 = tsne$Y[, 2]
)

# ------------------ Epithelial class mapping ------------------
# NOTE: ESCA + CESC are mixed cohorts in TCGA; BLCA is urothelial (transitional), kept as Mixed.
epi_map <- list(
  Adenocarcinoma = c(
    "BRCA","COAD","READ","STAD","PAAD","PRAD","LUAD","OV","UCEC","UCS",
    "CHOL","LIHC","THCA","KIRC","KIRP","KICH","ACC"
  ),
  Squamous = c(
    "HNSC","LUSC"
  ),
  Mixed = c(
    "ESCA","CESC","BLCA"
  ),
  NonEpithelial = c(
    "GBM","LGG","LAML","DLBC","SKCM","UVM","SARC",
    "TGCT","PCPG","MESO","THYM"
  )
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

# ------------------ Colors for epithelial classes ------------------
# (fixed, readable palette; independent of number of cancers)
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
    subtitle = paste0("Top ", n_var_probes, " variable probes; PCA=", min(n_pcs, ncol(X_pca)),
                      " dims; perplexity=", perp),
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


#!/usr/bin/env Rscript
# t-SNE of TCGA tumor methylation (HM450) — 01A only
# Uses a PROVIDED CpG list (e.g. 10,000 CpGs from cpg_list.txt), not re-selecting top SD probes.
# Colors points by epithelial class. Saves PDF + coordinates TSV/RDS.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(Rtsne)
})

# ------------------ Inputs / outputs ------------------
infile      <- "./methylation/filtered_methylation_files_tumor01A_1.txt"
top10k_file <- "./cpg_list.txt"

outpdf <- "./results/summarys/tSNE_methylation_01A_epithelialColors_top10kFromFile.pdf"
out_tsv <- "./results/summarys/tSNE_methylation_01A_epithelialCoords_top10kFromFile.tsv"
out_rds <- "./results/summarys/tSNE_methylation_01A_epithelialCoords_top10kFromFile.rds"

# ------------------ Parameters ------------------
n_pcs          <- 30       # PCA dims before tSNE
max_iter_tsne  <- 1000
seed           <- 1

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
if (length(files_all) < 2) stop("Not enough 01A files after filtering.")

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
    X[i, ] <- dt[.(top_probes), beta]   # fill in exact order of top_probes
  }

  if (i %% 100 == 0) {
    cat(sprintf("   - Filled %d / %d samples (%.1f min)\n",
                i, n_samples, (proc.time()[3] - t0) / 60))
  }
}

# ------------------ Impute missing + scale (safe) ------------------
cat(">>> Imputing missing values (per CpG mean) + scaling...\n")
probe_means <- colMeans(X, na.rm = TRUE)
for (j in seq_len(ncol(X))) {
  miss <- is.na(X[, j])
  if (any(miss)) X[miss, j] <- probe_means[j]
}

# remove CpGs with zero variance across samples (constant columns)
sds <- apply(X, 2, sd, na.rm = TRUE)
keep <- is.finite(sds) & sds > 0
cat(">>> Removing constant CpGs:", sum(!keep), "\n")

X <- X[, keep, drop = FALSE]

# scale after filtering
Xs <- scale(X)


# ------------------ PCA + t-SNE ------------------
cat(">>> Running PCA...\n")
pca <- prcomp(Xs, center = FALSE, scale. = FALSE)
X_pca <- pca$x[, seq_len(min(n_pcs, ncol(pca$x))), drop = FALSE]

n_embed <- nrow(X_pca)
perp_max <- floor((n_embed - 1) / 3) - 1
if (perp_max < 5) stop("Too few samples for t-SNE with perplexity>=5 (n=", n_embed, ")")
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
  Adenocarcinoma = c(
    "BRCA","COAD","READ","STAD","PAAD","PRAD","LUAD","OV","UCEC","UCS",
    "CHOL","LIHC","THCA","KIRC","KIRP","KICH","ACC"
  ),
  Squamous = c("HNSC","LUSC"),
  Mixed = c("ESCA","CESC","BLCA"),
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
    subtitle = paste0("CpGs used: ", ncol(X), "; PCA dims: ", min(n_pcs, ncol(X_pca)),
                      "; perplexity=", perp),
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

# * **Samples:** Tumor HM450 (01A), **≤ 200 samples per cancer type**
#* **Features:** **3,000 most variable CpG probes** (global SD)
#* **Scaling:** Z-score per probe
#* **UMAP:** `n_neighbors = 30`, `min_dist = 0.3`, Euclidean distance
#* **Output:** 2D embedding colored by cancer type

#!/usr/bin/env Rscript
library(data.table)
library(ggplot2)
library(uwot)

# ------------------ Inputs / outputs ------------------
infile <- "./methylation/methylation_files_tumor0x.txt"
outpdf <- "./results/summarys/UMAP_methylation_01A_top3000SD_max200perCancer.pdf"

# ------------------ Parameters ------------------
max_per_cancer <- 200
n_var_probes   <- 3000
n_neighbors    <- 30
min_dist       <- 0.3

# ------------------ Load file list ------------------
files_all <- scan(infile, what = "", quiet = TRUE)
files_all <- files_all[file.exists(files_all)]
if (length(files_all) < 2) stop("Not enough existing files in infile list.")

sample_names_all <- basename(files_all)
cancer_all <- sapply(strsplit(sample_names_all, "-"), `[`, 2)

# ------------------ Limit to max_per_cancer per cancer type ------------------
set.seed(1)
idx_by_cancer <- split(seq_along(files_all), cancer_all)

keep_idx <- unlist(lapply(idx_by_cancer, function(idx) {
  if (length(idx) <= max_per_cancer) idx else sample(idx, max_per_cancer)
}))

files_all <- files_all[keep_idx]
sample_names_all <- sample_names_all[keep_idx]
cancer_all <- cancer_all[keep_idx]

# Optional: stable ordering (nice for reproducibility/plotting)
ord <- order(cancer_all, sample_names_all)
files_all <- files_all[ord]
sample_names_all <- sample_names_all[ord]
cancer_all <- cancer_all[ord]

cat("Kept", length(files_all), "samples after max_per_cancer filter\n")
# Uncomment if you want to verify:
# print(sort(table(cancer_all), decreasing = TRUE))

# ------------------ Reader ------------------
read_one <- function(f) {
  dt <- fread(cmd = paste("zcat", shQuote(f)), header = FALSE, select = c(4, 5))
  setnames(dt, c("probe", "beta"))
  dt[, beta := as.numeric(beta)]
  dt
}

# ------------------ PASS 1: running mean/variance per probe ------------------
cat("PASS 1: computing probe SDs without building full matrix...\n")

dt0 <- read_one(files_all[1])
probes <- dt0$probe
p <- length(probes)

mean_vec <- dt0$beta
m2_vec   <- rep(0, p)
n <- 1L

for (i in 2:length(files_all)) {
  dt <- read_one(files_all[i])
  if (!identical(dt$probe, probes)) stop("Probe order mismatch: ", files_all[i])

  n <- n + 1L
  x <- dt$beta

  delta   <- x - mean_vec
  mean_vec <- mean_vec + delta / n
  delta2  <- x - mean_vec
  m2_vec  <- m2_vec + delta * delta2

  if (i %% 50 == 0) cat("  processed", i, "of", length(files_all), "files\n")
}

if (n < 2) stop("Need at least 2 samples to compute SDs.")
var_vec <- m2_vec / (n - 1L)
sd_vec  <- sqrt(var_vec)

top_idx <- order(sd_vec, decreasing = TRUE)[seq_len(min(n_var_probes, length(sd_vec)))]
top_probes <- probes[top_idx]

rm(mean_vec, m2_vec, var_vec, sd_vec); gc()
cat("Selected top probes:", length(top_probes), "\n")

# ------------------ PASS 2: build matrix only for top probes ------------------
cat("PASS 2: building matrix for top probes only...\n")

mat_small <- matrix(NA_real_, nrow = length(top_probes), ncol = length(files_all))
rownames(mat_small) <- top_probes
colnames(mat_small) <- sample_names_all

for (i in seq_along(files_all)) {
  dt <- read_one(files_all[i])
  if (!identical(dt$probe, probes)) stop("Probe order mismatch in PASS 2: ", files_all[i])

  mat_small[, i] <- dt$beta[top_idx]

  if (i %% 50 == 0) cat("  filled", i, "of", length(files_all), "files\n")
}
gc()

# samples x probes
mat_top <- t(mat_small)
rm(mat_small); gc()

# Scale features (recommended for UMAP)
mat_top <- scale(mat_top)

# ------------------ Run UMAP ------------------
cat("Running UMAP...\n")
set.seed(1)
emb <- uwot::umap(
  mat_top,
  n_neighbors = n_neighbors,
  min_dist    = min_dist,
  metric      = "euclidean"
)

plot_dt <- data.table(
  UMAP1  = emb[, 1],
  UMAP2  = emb[, 2],
  sample = sample_names_all,
  cancer = cancer_all
)

# ------------------ Plot ------------------
dir.create(dirname(outpdf), recursive = TRUE, showWarnings = FALSE)

pdf(outpdf, width = 14, height = 8)
print(
  ggplot(plot_dt, aes(UMAP1, UMAP2, color = cancer)) +
    geom_point(alpha = 0.8, size = 1.2) +
    theme_minimal(base_size = 12) +
    labs(
      title    = paste0("UMAP of TCGA tumor methylation (HM450) – top ", n_var_probes, " variable probes"),
      subtitle = paste0(
        "max_per_cancer=", max_per_cancer,
        " | n_neighbors=", n_neighbors,
        " | min_dist=", min_dist,
        " | n_samples=", nrow(plot_dt)
      ),
      color = "Cancer type"
    )
)
dev.off()

cat("Wrote:", outpdf, "\n")

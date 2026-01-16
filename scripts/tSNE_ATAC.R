''' This script performs t-SNE analysis on ATAC-seq peak presence/absence data across all TCGA cancer types.
1)It collects all ATAC-seq peak BED files in `./peaks/filtered_peaks/` (TCGA cancers), extracts the cancer type from each filename, and reports how many samples/cancers are included.
2)It builds a global merged peak set by concatenating all peaks, sorting, and running `bedtools merge` to create one unified reference peak list.
3)For each sample, it computes a binary peak presence/absence vector over the merged peaks, estimates peak variability across samples, keeps the top 3000 most variable peaks, then rebuilds a reduced matrix (peaks × samples).
4)It transposes to samples × peaks, scales the features, runs PCA (up to 30 PCs), and then runs t-SNE on the PCA scores (perplexity chosen from sample count, seed fixed).
5)Finally, it saves the t-SNE coordinates + cancer labels to a TSV and produces a colored t-SNE scatter plot PDF with one point per sample.
'''

#!/usr/bin/env Rscript
library(data.table)
library(Rtsne)
library(ggplot2)
library(Polychrome)

# =======================
# SETTINGS
# =======================
peaks_dir <- "./peaks/filtered_peaks/"
file_pattern <- "^ATAC_TCGA-[A-Z]+_.*_peaks_macs\\.bed$"

out_pdf <- "./results/summarys/tSNE_ATAC_ALLCANCERS_peakPresence.pdf"
out_tsv <- "./results/summarys/tSNE_ATAC_ALLCANCERS_coords.tsv"

top_k_peaks <- 3000
n_pcs <- 30
seed <- 1

# =======================
# HELPERS
# =======================
get_cancer_from_filename <- function(fn) {
  m <- regmatches(fn, regexpr("TCGA-[A-Z]+", fn))
  if (length(m) == 0 || is.na(m) || m == "") NA_character_ else m
}

# =======================
# LIST FILES
# =======================
files <- list.files(peaks_dir, full.names = TRUE)
files <- files[grepl(file_pattern, basename(files))]

samples <- basename(files)
cancer  <- vapply(samples, get_cancer_from_filename, character(1))

message("Total samples: ", length(files))
message("Cancer types:  ", length(unique(cancer)))

# =======================
# 1) GLOBAL MERGED PEAK SET
# =======================
dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)

tmpdir <- tempdir()
merged_bed <- file.path(tmpdir, "merged_peaks_ALLCANCERS.bed")

message("Building global merged peaks...")
cmd_merge <- paste0(
  "cat ", paste(shQuote(files), collapse = " "),
  " | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge -i stdin > ", shQuote(merged_bed)
)
if (system(cmd_merge) != 0) stop("Failed to build merged peaks.")

merged <- fread(merged_bed, header = FALSE)
setnames(merged, c("chr","start","end"))
P <- nrow(merged)

# =======================
# 2) PASS 1: variance
# =======================
mu <- numeric(P)
M2 <- numeric(P)
n_seen <- 0L

for (i in seq_along(files)) {
  cmd_int <- paste(
    "bedtools intersect -c -a", shQuote(merged_bed),
    "-b", shQuote(files[i])
  )
  dtc <- fread(cmd = cmd_int, header = FALSE)
  x <- as.numeric(dtc[[4]] > 0)

  n_seen <- n_seen + 1L
  delta <- x - mu
  mu <- mu + delta / n_seen
  M2 <- M2 + delta * (x - mu)

  if (i %% 10 == 0) message("Pass 1:", i, "/", length(files))
}

sdv <- sqrt(M2 / pmax(n_seen - 1L, 1L))
ord <- order(sdv, decreasing = TRUE)
k <- min(top_k_peaks, length(ord))
keep_idx <- ord[1:k]

# =======================
# 3) PASS 2: reduced matrix
# =======================
Xk <- matrix(0L, nrow = k, ncol = length(files))
rownames(Xk) <- paste0(merged$chr[keep_idx], ":", merged$start[keep_idx], "-", merged$end[keep_idx])
colnames(Xk) <- samples

for (i in seq_along(files)) {
  cmd_int <- paste(
    "bedtools intersect -c -a", shQuote(merged_bed),
    "-b", shQuote(files[i])
  )
  dtc <- fread(cmd = cmd_int, header = FALSE)
  Xk[, i] <- as.integer(dtc[[4]] > 0)[keep_idx]

  if (i %% 10 == 0) message("Pass 2:", i, "/", length(files))
}

# =======================
# t-SNE
# =======================
Xs <- scale(t(Xk))
pcs <- min(n_pcs, ncol(Xs) - 1, nrow(Xs) - 1)
pca <- prcomp(Xs, center = FALSE, scale. = FALSE)
Y <- pca$x[, 1:pcs, drop = FALSE]

n <- nrow(Y)
perp <- max(5, min(30, floor((n - 1) / 3)))

set.seed(seed)
ts <- Rtsne(Y, perplexity = perp, check_duplicates = FALSE, pca = FALSE)

plot_dt <- data.table(
  sample = rownames(Y),
  cancer = cancer[match(rownames(Y), samples)],
  tSNE1 = ts$Y[,1],
  tSNE2 = ts$Y[,2]
)

fwrite(plot_dt, out_tsv, sep = "\t")

# =======================
# DISTINCT COLORS (Polychrome)
# =======================
cancers <- sort(unique(plot_dt$cancer))
pal <- Polychrome::createPalette(length(cancers), seedcolors = "#000000")
names(pal) <- cancers

# =======================
# PLOT
# =======================
p <- ggplot(plot_dt, aes(tSNE1, tSNE2, color = cancer)) +
  geom_point(alpha = 0.85, size = 1.4) +
  scale_color_manual(values = pal) +
  theme_bw() +
  labs(
    title = "t-SNE on ATAC peak presence/absence (ALL TCGA cancers)",
    subtitle = paste0("Merged peaks=", P,
                      "; top variable peaks=", k,
                      "; PCA=", pcs, " dims; perplexity=", perp),
    x = "tSNE1", y = "tSNE2", color = "Cancer"
  ) +
  theme(legend.position = "right")

ggsave(out_pdf, p, width = 12, height = 8)

cat("Wrote:", out_pdf, "\n")
cat("Wrote:", out_tsv, "\n")

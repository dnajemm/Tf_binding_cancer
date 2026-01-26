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


'''
#!/usr/bin/env Rscript
# One-page PDF with 3 separate t-SNE panels (Cancer / Organ group / Tissue type),
# each with its OWN Polychrome legend placed BELOW its plot.
# Input: TSV with columns: sample, cancer, tSNE1, tSNE2

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(Polychrome)
  library(patchwork)
})

# ------------------ inputs / outputs ------------------
in_tsv  <- "./results/summarys/tSNE_ATAC_ALLCANCERS_coords.tsv"
out_pdf <- "./results/summarys/tSNE_ATAC_3views.pdf"

# ------------------ load tSNE coords ------------------
dt <- fread(in_tsv)
stopifnot(all(c("sample","cancer","tSNE1","tSNE2") %in% names(dt)))

# ------------------ mappings ------------------
organ_map <- list(
  CNS           = c("GBM","LGG"),
  Breast        = c("BRCA"),
  Gynecologic   = c("CESC","OV","UCEC","UCS"),
  Lung          = c("LUAD","LUSC"),
  HeadNeck      = c("HNSC"),
  ThyroidThymus = c("THCA","THYM"),
  SkinEye       = c("SKCM","UVM"),
  GI            = c("COAD","READ","ESCA","STAD"),
  Hepatobiliary = c("LIHC","CHOL"),
  Pancreas      = c("PAAD"),
  Kidney        = c("KICH","KIRC","KIRP"),
  Bladder       = c("BLCA"),
  Prostate      = c("PRAD"),
  Testis        = c("TGCT"),
  Adrenal       = c("ACC","PCPG"),
  Sarcoma       = c("SARC"),
  HemeLymph     = c("LAML","DLBC"),
  Mesothelioma  = c("MESO")
)

tissue_map <- list(
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

map_group <- function(ct, mapping) {
  ct2 <- gsub("^TCGA-", "", ct)
  for (g in names(mapping)) if (ct2 %in% mapping[[g]]) return(g)
  "Other"
}

# add annotations
dt[, organ_group := factor(vapply(cancer, map_group, character(1), mapping = organ_map))]
dt[, tissue_type := factor(vapply(cancer, map_group, character(1), mapping = tissue_map))]
dt[, cancer := factor(cancer)]

# ------------------ Polychrome palettes ------------------
pal_cancer <- Polychrome::createPalette(
  nlevels(dt$cancer),
  seedcolors = c("#ff0000ff", "#E69F00", "#56B4E9", "#009E73")
)
names(pal_cancer) <- levels(dt$cancer)

pal_organ <- Polychrome::createPalette(
  nlevels(dt$organ_group),
  seedcolors = c("#ff00c3ff", "#815801ff", "#24099eff", "#55ffd2ff")
)
names(pal_organ) <- levels(dt$organ_group)

pal_tissue <- Polychrome::createPalette(
  nlevels(dt$tissue_type),
  seedcolors = c("#2bff00ff", "#d9ff00ff", "#9e8dffff", "#004a36ff")
)
names(pal_tissue) <- levels(dt$tissue_type)

# ------------------ plots (each with its own legend below) ------------------
base_theme <- theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 13),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.45, "cm")
  )

p_cancer <- ggplot(dt, aes(tSNE1, tSNE2, color = cancer)) +
  geom_point(size = 1.4, alpha = 0.85) +
  scale_color_manual(values = pal_cancer) +
  labs(title = "Cancer type", x = "tSNE1", y = "tSNE2") +
  base_theme

p_organ <- ggplot(dt, aes(tSNE1, tSNE2, color = organ_group)) +
  geom_point(size = 1.4, alpha = 0.85) +
  scale_color_manual(values = pal_organ) +
  labs(title = "Organ group", x = "tSNE1", y = "tSNE2") +
  base_theme

p_tissue <- ggplot(dt, aes(tSNE1, tSNE2, color = tissue_type)) +
  geom_point(size = 1.4, alpha = 0.85) +
  scale_color_manual(values = pal_tissue) +
  labs(title = "Tissue type", x = "tSNE1", y = "tSNE2") +
  base_theme

# ------------------ export one-page PDF ------------------
dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)

pdf(out_pdf, width = 18, height = 7)
(p_cancer | p_organ | p_tissue) +
  plot_annotation(
    title = "t-SNE of ATAC-seq peak presence/absence (TCGA)",
    subtitle = "Same t-SNE coordinates, colored by: Cancer type / Organ group / Tissue type"
  )
dev.off()

cat("Wrote:", out_pdf, "\n")

''''

##################################################################################################################################
# t-SNE of TCGA tumor methylation (HM450) by beta value, using paired 01A+11A samples to cluster by methylation levels.
##################################################################################################################################
#!/usr/bin/env Rscript
# t-SNE of TCGA tumor methylation (HM450) — 11A only
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

outpdf <- "./results/methylation/tSNE_methylation/tSNE_methylation_11A_top10k.pdf"
out_tsv <- "./results/methylation/tSNE_methylation/tSNE_methylation_11A_top10k_coords.tsv"
out_rds <- "./results/methylation/tSNE_methylation/tSNE_methylation_11A_top10k_coords.rds"

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

###############################################################################################################
# The code is for paired 01A+11A samples to cluster by beta value, which is more complex but more informative.
###############################################################################################################
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

outpdf  <- "./results/methylation/tSNE_methylation/tSNE_methylation_01A_11A_top10k.pdf"
out_tsv <- "./results/methylation/tSNE_methylation/tSNE_methylation_01A_11A_top10k_coords.tsv"
out_rds <- "./results/methylation/tSNE_methylation/tSNE_methylation_01A_11A_top10k_coords.rds"

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


###############################################################################################################
# The code is for 01A samples only, colored by epithelial class (adenocarcinoma, squamous, etc.). 
###############################################################################################################
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

outpdf <- "./results/methylation/tSNE_methylation_01A_epithelialColors_top10k.pdf"
out_tsv <- "./results/methylation/tSNE_methylation_01A_epithelialCoords_top10k.tsv"
out_rds <- "./results/methylation/tSNE_methylation_01A_epithelialCoords_top10k.rds"

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


##############################################################################################
# tSNE on the 3d+4d samples only , redoing the 10k cpg most variable across those samples
##############################################################################################
mkdir -p ./results/methylation/tSNE_methylation/tsne_3d4d_noOV

##############################################################################################
# Step 1: extract the 3d+4d_noOV patient IDs from the cohort file
##############################################################################################
Rscript - <<'EOF'
suppressPackageStartupMessages(library(data.table))

infile  <- "./results/multi_omics/samples_3d+4d_noOV.tsv"
outfile <- "./results/methylation/tSNE_methylation/tsne_3d4d_noOV/patients_3d4d_noOV.tsv"

dt <- fread(infile, header = FALSE)

if (ncol(dt) < 1) stop("Input file has no columns")

setnames(dt, 1, "patient")

out <- unique(dt[, .(patient = as.character(patient))])

fwrite(out, outfile, sep = "\t")

cat("Wrote:", outfile, "\n")
cat("Number of patients:", nrow(out), "\n")
cat("First patients:\n")
print(head(out))
EOF

##############################################################################################
# Step 2: build the 3d+4d_noOV methylation sample list from filtered methylation filenames
##############################################################################################
Rscript - <<'EOF'
suppressPackageStartupMessages(library(data.table))

patient_file <- "./results/methylation/tSNE_methylation/tsne_3d4d_noOV/patients_3d4d_noOV.tsv"
meth_dir     <- "./methylation/filtered_methylation"
outfile      <- "./results/methylation/tSNE_methylation/tsne_3d4d_noOV/methylation_samples_3d4d_noOV.tsv"

patients <- fread(patient_file)
patients[, patient := as.character(patient)]

files <- list.files(meth_dir, pattern = "\\.bed\\.gz$", full.names = TRUE)
if (length(files) == 0) stop("No .bed.gz files found in: ", meth_dir)

meta <- data.table(file = files)
meta[, basename := basename(file)]

# Extract sample barcode like TCGA-OR-A5J6-01A from filenames like:
# HM450_TCGA-ACC-TCGA-OR-A5J6-01A_1_annotated_methylation_filtered.bed.gz
meta[, sample := sub("^.*(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}[A-Z]).*$", "\\1", basename)]

# Keep only rows where extraction worked
meta <- meta[grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}[A-Z]$", sample)]

# patient = first 3 TCGA fields
meta[, patient := sub("^(([^-]+-){2}[^-]+).*", "\\1", sample)]

# sample type from barcode
meta[, sample_type_code := sub("^TCGA-[^-]+-[^-]+-([0-9]{2}).*$", "\\1", sample)]
meta[, sample_type := fifelse(sample_type_code == "01", "Tumor",fifelse(sample_type_code == "11", "Healthy", "Other"))]

# keep only 3d+4d cohort patients
meta <- meta[patient %in% patients$patient]

# keep only tumor and healthy samples
meta <- meta[sample_type %in% c("Tumor", "Healthy")]

# remove duplicates if any
meta <- unique(meta[, .(file, sample, patient, sample_type)])

setorder(meta, patient, sample_type, sample)

fwrite(meta, outfile, sep = "\t")

cat("Patients in cohort:", nrow(patients), "\n")
cat("All methylation files found:", length(files), "\n")
cat("Matched methylation samples:", nrow(meta), "\n")
cat("Breakdown by sample type:\n")
print(meta[, .N, by = sample_type][order(sample_type)])
cat("Wrote:", outfile, "\n")
cat("First rows:\n")
print(head(meta))
EOF

##############################################################################################
# Step 3: keep only the 3d+4d_noOV methylation samples that are present in pan_mat
##############################################################################################
Rscript - <<'EOF'
suppressPackageStartupMessages(library(data.table))

cat("STEP 3/5 - Starting sample filtering against pan_mat\n")

sample_file  <- "./results/methylation/tSNE_methylation/tsne_3d4d_noOV/methylation_samples_3d4d_noOV.tsv"
pan_mat_rds  <- "./methylation/TCGA-all_cancer_pan_cancer_matrix.RData"
outfile      <- "./results/methylation/tSNE_methylation/tsne_3d4d_noOV/methylation_samples_3d4d_noOV_in_panmat.tsv"

cat("STEP 3/5 - Reading sample file\n")
meta <- fread(sample_file)
meta[, sample := as.character(sample)]

cat("STEP 3/5 - Number of samples from Step 2:", nrow(meta), "\n")

cat("STEP 3/5 - Loading pan_mat\n")
load(pan_mat_rds)
if (!exists("pan_mat")) stop("pan_mat not found in RData")

cat("STEP 3/5 - pan_mat loaded\n")
cat("STEP 3/5 - pan_mat dimensions:", nrow(pan_mat), "CpGs x", ncol(pan_mat), "samples\n")

cat("STEP 3/5 - Extracting TCGA sample barcodes from pan_mat column names\n")
pan_files <- colnames(pan_mat)

samples_in_panmat <- sub("^HM450_TCGA-[A-Z0-9]+-(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}[A-Z])_.*$", "\\1", pan_files)

cat("STEP 3/5 - First extracted sample barcodes from pan_mat:\n")
print(head(samples_in_panmat, 10))

cat("STEP 3/5 - Matching Step 2 samples to pan_mat columns\n")
out <- meta[sample %in% samples_in_panmat]

missing <- setdiff(meta$sample, samples_in_panmat)

cat("STEP 3/5 - Writing filtered sample table\n")
fwrite(out, outfile, sep = "\t")

cat("STEP 3/5 - DONE\n")
cat("Samples from Step 2:", nrow(meta), "\n")
cat("Samples found in pan_mat:", nrow(out), "\n")
cat("Samples missing from pan_mat:", length(missing), "\n")

if (nrow(out) > 0) {
  cat("Breakdown of kept samples by sample type:\n")
  print(out[, .N, by = sample_type][order(sample_type)])
} else {
  cat("No matched samples found\n")
}

cat("Wrote:", outfile, "\n")

if (length(missing) > 0) {
  cat("First missing samples:\n")
  print(head(missing, 20))
}
EOF


##############################################################################################
# Step 4: recompute the top 10k most variable CpGs across the 3d+4d_noOV samples
##############################################################################################
echo "Running Step 3: matching 3d+4d_noOV methylation samples to pan_mat columns"

Rscript - <<'EOF'
suppressPackageStartupMessages(library(data.table))

cat("STEP 3/5 - Starting sample filtering against pan_mat\n")

sample_file  <- "./results/methylation/tSNE_methylation/tsne_3d4d_noOV/methylation_samples_3d4d_noOV.tsv"
pan_mat_rds  <- "./methylation/TCGA-all_cancer_pan_cancer_matrix.RData"
outfile      <- "./results/methylation/tSNE_methylation/tsne_3d4d_noOV/methylation_samples_3d4d_noOV_in_panmat.tsv"

cat("STEP 3/5 - Reading sample file\n")
meta <- fread(sample_file)
meta[, sample := as.character(sample)]

cat("STEP 3/5 - Number of samples from Step 2:", nrow(meta), "\n")

cat("STEP 3/5 - Loading pan_mat\n")
load(pan_mat_rds)
if (!exists("pan_mat")) stop("pan_mat not found in RData")

cat("STEP 3/5 - pan_mat loaded\n")
cat("STEP 3/5 - pan_mat dimensions:", nrow(pan_mat), "CpGs x", ncol(pan_mat), "samples\n")

cat("STEP 3/5 - Building pan_mat annotation from column names\n")
pan_annot <- data.table(
  pan_mat_col = colnames(pan_mat)
)

pan_annot[, sample := sub(
  "^HM450_TCGA-[A-Z0-9]+-(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}[A-Z])_.*$",
  "\\1",
  pan_mat_col
)]

cat("STEP 3/5 - First extracted pan_mat sample barcodes:\n")
print(head(pan_annot, 10))

cat("STEP 3/5 - Merging Step 2 sample table with pan_mat annotation\n")
out <- merge(meta, pan_annot, by = "sample", all = FALSE)

missing <- setdiff(meta$sample, pan_annot$sample)

cat("STEP 3/5 - Writing filtered sample table\n")
fwrite(out, outfile, sep = "\t")

cat("STEP 3/5 - DONE\n")
cat("Samples from Step 2:", nrow(meta), "\n")
cat("Samples found in pan_mat:", nrow(out), "\n")
cat("Samples missing from pan_mat:", length(missing), "\n")

if (nrow(out) > 0) {
  cat("Breakdown of kept samples by sample type:\n")
  print(out[, .N, by = sample_type][order(sample_type)])
  cat("Output columns:\n")
  print(colnames(out))
}

cat("Wrote:", outfile, "\n")

if (length(missing) > 0) {
  cat("First missing samples:\n")
  print(head(missing, 20))
}
EOF

# step 4 
echo "Running Step 4: recomputing top 10k variable CpGs in 3d+4d_noOV"

Rscript - <<'EOF'
suppressPackageStartupMessages(library(data.table))

cat("STEP 4/5 - Starting top 10k CpG selection\n")

sample_file <- "./results/methylation/tSNE_methylation/tsne_3d4d_noOV/methylation_samples_3d4d_noOV_in_panmat.tsv"
pan_mat_rds <- "./methylation/TCGA-all_cancer_pan_cancer_matrix.RData"
out_file    <- "./results/methylation/tSNE_methylation/tsne_3d4d_noOV/top10000_variable_CpGs_3d4d_noOV.tsv"

cat("STEP 4/5 - Reading kept sample list\n")
meta <- fread(sample_file)

if (!"pan_mat_col" %in% colnames(meta)) {
  stop("Column 'pan_mat_col' not found in sample_file. Re-run corrected Step 3 first.")
}

samples <- unique(as.character(meta$pan_mat_col))
cat("STEP 4/5 - Number of kept pan_mat columns:", length(samples), "\n")
if (length(samples) < 2) stop("Need at least 2 samples to compute variance.")

cat("STEP 4/5 - Loading pan_mat\n")
load(pan_mat_rds)
if (!exists("pan_mat")) stop("pan_mat not found in RData")
cat("STEP 4/5 - pan_mat dimensions:", nrow(pan_mat), "CpGs x", ncol(pan_mat), "samples\n")

cat("STEP 4/5 - Checking that kept columns exist in pan_mat\n")
missing_cols <- setdiff(samples, colnames(pan_mat))
cat("STEP 4/5 - Missing columns in pan_mat:", length(missing_cols), "\n")
if (length(missing_cols) > 0) {
  cat("First missing columns:\n")
  print(head(missing_cols, 20))
  stop("Some pan_mat_col values are not present in pan_mat colnames.")
}

cat("STEP 4/5 - Subsetting pan_mat to the kept samples\n")
mat <- pan_mat[, samples, drop = FALSE]
cat("STEP 4/5 - Subset matrix dimensions:", nrow(mat), "CpGs x", ncol(mat), "samples\n")

cat("STEP 4/5 - Removing CpGs with <= 1 observed value across samples\n")
keep_obs <- rowSums(!is.na(mat)) > 1
cat("STEP 4/5 - CpGs removed for too few observed values:", sum(!keep_obs), "\n")
mat <- mat[keep_obs, , drop = FALSE]

cat("STEP 4/5 - Computing variance per CpG across kept samples\n")
vars <- apply(mat, 1, var, na.rm = TRUE)

cat("STEP 4/5 - Removing non-finite / constant CpGs\n")
keep_var <- is.finite(vars) & vars > 0
cat("STEP 4/5 - CpGs removed for zero/non-finite variance:", sum(!keep_var), "\n")
vars <- vars[keep_var]

cat("STEP 4/5 - Ranking CpGs by variance\n")
top_n <- min(10000, length(vars))
top_cpgs <- names(sort(vars, decreasing = TRUE))[seq_len(top_n)]

out <- data.table(CpG = top_cpgs)
fwrite(out, out_file, sep = "\t")

cat("STEP 4/5 - DONE\n")
cat("Top CpGs selected:", nrow(out), "\n")
cat("Wrote:", out_file, "\n")
cat("First CpGs:\n")
print(head(out))
EOF

##############################################################################################################
# step 5: run the tSNE code with the new sample list and new CpG list
##############################################################################################################
echo "Running Step 5: t-SNE on 3d+4d_noOV methylation samples"

Rscript - <<'EOF'
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(Rtsne)
})

cat("STEP 5/5 - Starting t-SNE\n")

sample_file <- "./results/methylation/tSNE_methylation/tsne_3d4d_noOV/methylation_samples_3d4d_noOV_in_panmat.tsv"
top10k_file <- "./results/methylation/tSNE_methylation/tsne_3d4d_noOV/top10000_variable_CpGs_3d4d_noOV.tsv"
pan_mat_rds <- "./methylation/TCGA-all_cancer_pan_cancer_matrix.RData"
color_file  <- "./results/multi_omics/cancer_color_order_with_defined_colours.tsv"

out_dir     <- "./results/methylation/tSNE_methylation/tsne_3d4d_noOV"
out_pdf     <- file.path(out_dir, "tSNE_methylation_3d4d_noOV_top10k.pdf")
out_tsv     <- file.path(out_dir, "tSNE_methylation_3d4d_noOV_top10k_coordinates.tsv")
out_rds     <- file.path(out_dir, "tSNE_methylation_3d4d_noOV_top10k_rtsne.rds")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("STEP 5/5 - Reading matched sample file\n")
meta <- fread(sample_file)

if (!"pan_mat_col" %in% colnames(meta)) {
  stop("Column 'pan_mat_col' not found in sample file.")
}

cat("STEP 5/5 - Reading top CpGs\n")
top_dt <- fread(top10k_file)
if (!"CpG" %in% colnames(top_dt)) {
  stop("Column 'CpG' not found in top10k file.")
}

cat("STEP 5/5 - Reading cancer color file\n")
col_dt <- fread(color_file, header = FALSE)
if (ncol(col_dt) < 2) {
  stop("color_file must have at least 2 columns: cancer and colour")
}
col_dt <- col_dt[, 1:2]
setnames(col_dt, c("cancer", "colour"))

cat("STEP 5/5 - First rows of color file:\n")
print(head(col_dt))

cat("STEP 5/5 - Building cancer annotation\n")
if (!"cancer" %in% colnames(meta)) {
  meta[, cancer := sub("^HM450_(TCGA-[A-Z0-9]+)-TCGA-.*$", "\\1", file)]
}

cat("STEP 5/5 - Unique cancers in meta:\n")
print(sort(unique(meta$cancer)))

cat("STEP 5/5 - Unique cancers in color file:\n")
print(sort(unique(col_dt$cancer)))

cancer_cols <- setNames(col_dt$colour, col_dt$cancer)

missing_palette <- setdiff(unique(meta$cancer), names(cancer_cols))
cat("STEP 5/5 - Cancers missing from palette:", length(missing_palette), "\n")
if (length(missing_palette) > 0) {
  print(missing_palette)
}

samples <- unique(as.character(meta$pan_mat_col))
cpgs    <- unique(as.character(top_dt$CpG))

cat("STEP 5/5 - Number of samples:", length(samples), "\n")
cat("STEP 5/5 - Number of top CpGs:", length(cpgs), "\n")

if (length(samples) < 2) stop("Too few samples for t-SNE.")
if (length(cpgs) < 2) stop("Too few CpGs for t-SNE.")

cat("STEP 5/5 - Loading pan_mat\n")
load(pan_mat_rds)
if (!exists("pan_mat")) stop("pan_mat not found in RData")

cat("STEP 5/5 - pan_mat dimensions:", nrow(pan_mat), "CpGs x", ncol(pan_mat), "samples\n")

missing_samples <- setdiff(samples, colnames(pan_mat))
missing_cpgs    <- setdiff(cpgs, rownames(pan_mat))

cat("STEP 5/5 - Missing samples in pan_mat:", length(missing_samples), "\n")
cat("STEP 5/5 - Missing CpGs in pan_mat:", length(missing_cpgs), "\n")

if (length(missing_samples) > 0) {
  cat("First missing samples:\n")
  print(head(missing_samples, 20))
  stop("Some samples are missing from pan_mat.")
}

cpgs <- intersect(cpgs, rownames(pan_mat))
if (length(cpgs) < 2) stop("Too few CpGs found in pan_mat.")

cat("STEP 5/5 - Subsetting pan_mat\n")
mat <- pan_mat[cpgs, samples, drop = FALSE]
cat("STEP 5/5 - Subset matrix dimensions:", nrow(mat), "CpGs x", ncol(mat), "samples\n")

cat("STEP 5/5 - Transposing to samples x CpGs\n")
x <- t(mat)

cat("STEP 5/5 - Removing CpGs with any NA across kept samples\n")
keep_complete <- colSums(is.na(x)) == 0
cat("STEP 5/5 - CpGs removed due to NA:", sum(!keep_complete), "\n")
x <- x[, keep_complete, drop = FALSE]

cat("STEP 5/5 - Final matrix dimensions for t-SNE:", nrow(x), "samples x", ncol(x), "CpGs\n")

if (nrow(x) < 2) stop("Too few samples for t-SNE.")
if (ncol(x) < 2) stop("Too few CpGs for t-SNE after NA filtering.")

set.seed(1)

perp <- min(30, floor((nrow(x) - 1) / 3))
if (perp < 2) perp <- 2
cat("STEP 5/5 - Using perplexity:", perp, "\n")

cat("STEP 5/5 - Running Rtsne\n")
ts <- Rtsne(
  x,
  dims = 2,
  perplexity = perp,
  check_duplicates = FALSE,
  pca = TRUE,
  verbose = TRUE,
  max_iter = 1000
)

coord <- data.table(
  pan_mat_col = rownames(x),
  tSNE1 = ts$Y[,1],
  tSNE2 = ts$Y[,2]
)

coord <- merge(coord, unique(meta), by = "pan_mat_col", all.x = TRUE)

cat("STEP 5/5 - Writing coordinates\n")
fwrite(coord, out_tsv, sep = "\t")

cat("STEP 5/5 - Saving Rtsne object\n")
saveRDS(ts, out_rds)

cat("STEP 5/5 - Plotting PDF\n")
p <- ggplot(coord, aes(x = tSNE1, y = tSNE2, color = cancer, shape = sample_type)) +
  geom_point(size = 2.2, alpha = 0.9) +
  scale_color_manual(values = cancer_cols, drop = FALSE, na.value = "grey70") +
  theme_bw(base_size = 12) +
  labs(
    title = "t-SNE of TCGA HM450 methylation",
    subtitle = paste0("3d+4d_noOV cohort | top ", nrow(top_dt), " variable CpGs"),
    x = "tSNE1",
    y = "tSNE2",
    color = "Cancer",
    shape = "Sample type"
  )

ggsave(out_pdf, p, width = 10, height = 7)

cat("STEP 5/5 - DONE\n")
cat("Wrote coordinates:", out_tsv, "\n")
cat("Wrote plot:", out_pdf, "\n")
cat("Wrote RDS:", out_rds, "\n")
EOF

###################################################
# adding the 3 views of the tSNE like in the ATAC 
###################################################
##############################################################################################################
# step 5: run the tSNE code with the new sample list and new CpG list
# corrected version:
# - ALWAYS rebuilds annotations from pan_mat_col / file
# - removes CpGs with zero variance across selected samples
# - imputes NA by CpG median within cancer x sample_type, only if NA fraction in that group is < 50%
# - drops CpGs still containing NA after group-wise imputation
# - 3 vertical views:
#     1) Cancer type  (keeps Tumor/Healthy separation by shape)
#     2) Organ group  (color only)
#     3) Tissue type  (color only)
# - larger PDF, slightly smaller points, larger text
##############################################################################################################
echo "Running Step 5: t-SNE on 3d+4d_noOV methylation samples"

Rscript - <<'EOF'
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(Rtsne)
  library(Polychrome)
  library(patchwork)
  library(grid)
})

cat("STEP 5/5 - Starting t-SNE\n")

sample_file <- "./results/methylation/tSNE_methylation/tsne_3d4d_noOV/methylation_samples_3d4d_noOV_in_panmat.tsv"
top10k_file <- "./results/methylation/tSNE_methylation/tsne_3d4d_noOV/top10000_variable_CpGs_3d4d_noOV.tsv"
pan_mat_rds <- "./methylation/TCGA-all_cancer_pan_cancer_matrix.RData"
color_file  <- "./results/multi_omics/cancer_color_order_with_defined_colours.tsv"

out_dir     <- "./results/methylation/tSNE_methylation/tsne_3d4d_noOV"
out_pdf     <- file.path(out_dir, "tSNE_methylation_3d4d_noOV_top10k_3views_vertical.pdf")
out_tsv     <- file.path(out_dir, "tSNE_methylation_3d4d_noOV_top10k_3viewscoordinates.tsv")
out_rds     <- file.path(out_dir, "tSNE_methylation_3d4d_noOV_top10k_3viewsrtsne.rds")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

clean_txt <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub('^"|"$', "", x)
  x <- gsub("^'|'$", "", x)
  x
}

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
  ct2 <- gsub("^TCGA-", "", as.character(ct))
  for (g in names(mapping)) {
    if (ct2 %in% mapping[[g]]) return(g)
  }
  "Other"
}

cat("STEP 5/5 - Reading matched sample file\n")
meta <- fread(sample_file)

if (!"pan_mat_col" %in% colnames(meta)) {
  stop("Column 'pan_mat_col' not found in sample file.")
}

if (!"file" %in% colnames(meta)) {
  meta[, file := pan_mat_col]
}

cat("STEP 5/5 - Reading top CpGs\n")
top_dt <- fread(top10k_file)
if (!"CpG" %in% colnames(top_dt)) {
  stop("Column 'CpG' not found in top10k file.")
}

cat("STEP 5/5 - Reading cancer color file\n")
col_dt <- fread(color_file, header = FALSE)
if (ncol(col_dt) < 2) {
  stop("color_file must have at least 2 columns: cancer and colour")
}
col_dt <- col_dt[, 1:2]
setnames(col_dt, c("cancer_raw", "colour"))
col_dt[, cancer_raw := clean_txt(cancer_raw)]
col_dt[, colour := clean_txt(colour)]
col_dt[, cancer_plot := sub("^TCGA-", "", cancer_raw)]

cat("STEP 5/5 - First rows of color file:\n")
print(head(col_dt))

cat("STEP 5/5 - Rebuilding sample annotations from pan_mat_col / file\n")

meta[, pan_mat_col := as.character(pan_mat_col)]
meta[, file := as.character(file)]

meta[, sample := sub("^HM450_[A-Z0-9-]+-(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}[A-Z])_.*$", "\\1", pan_mat_col)]
bad_sample <- is.na(meta$sample) | meta$sample == meta$pan_mat_col
if (any(bad_sample)) {
  meta[bad_sample, sample := sub("^.*?(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}[A-Z]).*$", "\\1", file)]
}

meta[, patient := sub("^(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}).*$", "\\1", sample)]

meta[, sample_code := sub("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-([0-9]{2}).*$", "\\1", sample)]
meta[, sample_type := fifelse(
  sample_code == "01", "Tumor",
  fifelse(sample_code == "11", "Healthy", paste0("Type_", sample_code))
)]

meta[, cancer := sub("^HM450_(TCGA-[A-Z0-9]+)-TCGA-.*$", "\\1", pan_mat_col)]
bad_cancer <- is.na(meta$cancer) | meta$cancer == meta$pan_mat_col
if (any(bad_cancer)) {
  meta[bad_cancer, cancer := sub("^.*HM450_(TCGA-[A-Z0-9]+)-TCGA-.*$", "\\1", file)]
}

meta[, cancer := clean_txt(cancer)]
meta[, cancer_plot := sub("^TCGA-", "", cancer)]
meta[, sample_type := factor(sample_type, levels = c("Tumor", "Healthy"))]
meta[, sample_code := NULL]

cat("STEP 5/5 - Example rebuilt annotations:\n")
print(head(meta[, .(pan_mat_col, sample, patient, sample_type, cancer, cancer_plot)]))

cat("STEP 5/5 - Unique cancers in meta:\n")
print(sort(unique(meta$cancer_plot)))

cat("STEP 5/5 - Unique cancers in color file:\n")
print(sort(unique(col_dt$cancer_plot)))

cancer_cols <- setNames(col_dt$colour, col_dt$cancer_plot)

missing_palette <- setdiff(unique(meta$cancer_plot), names(cancer_cols))
cat("STEP 5/5 - Cancers missing from palette:", length(missing_palette), "\n")
if (length(missing_palette) > 0) {
  print(missing_palette)
  for (cc in missing_palette) cancer_cols[cc] <- "#7f7f7f"
}

samples <- unique(as.character(meta$pan_mat_col))
cpgs    <- unique(as.character(top_dt$CpG))

cat("STEP 5/5 - Number of samples:", length(samples), "\n")
cat("STEP 5/5 - Number of top CpGs:", length(cpgs), "\n")

if (length(samples) < 2) stop("Too few samples for t-SNE.")
if (length(cpgs) < 2) stop("Too few CpGs for t-SNE.")

cat("STEP 5/5 - Loading pan_mat\n")
load(pan_mat_rds)
if (!exists("pan_mat")) stop("pan_mat not found in RData")

cat("STEP 5/5 - pan_mat dimensions:", nrow(pan_mat), "CpGs x", ncol(pan_mat), "samples\n")

missing_samples <- setdiff(samples, colnames(pan_mat))
missing_cpgs    <- setdiff(cpgs, rownames(pan_mat))

cat("STEP 5/5 - Missing samples in pan_mat:", length(missing_samples), "\n")
cat("STEP 5/5 - Missing CpGs in pan_mat:", length(missing_cpgs), "\n")

if (length(missing_samples) > 0) {
  cat("First missing samples:\n")
  print(head(missing_samples, 20))
  stop("Some samples are missing from pan_mat.")
}

cpgs <- intersect(cpgs, rownames(pan_mat))
if (length(cpgs) < 2) stop("Too few CpGs found in pan_mat.")

cat("STEP 5/5 - Subsetting pan_mat\n")
mat <- pan_mat[cpgs, samples, drop = FALSE]
cat("STEP 5/5 - Subset matrix dimensions:", nrow(mat), "CpGs x", ncol(mat), "samples\n")

cat("STEP 5/5 - Transposing to samples x CpGs\n")
x <- t(mat)
cat("STEP 5/5 - Initial matrix dimensions:", nrow(x), "samples x", ncol(x), "CpGs\n")

cat("STEP 5/5 - Removing zero-variance / non-finite-variance CpGs across selected samples\n")
vars0 <- apply(x, 2, function(v) var(v, na.rm = TRUE))
keep_var <- is.finite(vars0) & vars0 > 0
cat("STEP 5/5 - CpGs removed due to zero/non-finite variance:", sum(!keep_var), "\n")
x <- x[, keep_var, drop = FALSE]
cat("STEP 5/5 - Matrix after variance filter:", nrow(x), "samples x", ncol(x), "CpGs\n")

if (ncol(x) < 2) stop("Too few CpGs after zero-variance filtering.")

cat("STEP 5/5 - Group-wise NA imputation by cancer x sample_type\n")

meta_sub <- unique(meta[, .(pan_mat_col, cancer_plot, sample_type)])
setkey(meta_sub, pan_mat_col)

sample_order <- rownames(x)
meta_sub <- meta_sub[J(sample_order)]

if (any(is.na(meta_sub$cancer_plot))) stop("Missing cancer_plot after matching metadata to x rows.")
if (any(is.na(meta_sub$sample_type))) stop("Missing sample_type after matching metadata to x rows.")

groups <- unique(meta_sub[, .(cancer_plot, sample_type)])
cat("STEP 5/5 - Number of cancer x sample_type groups:", nrow(groups), "\n")

na_before_total <- sum(is.na(x))
cat("STEP 5/5 - Total NA values before imputation:", na_before_total, "\n")

n_imputed <- 0L
n_blocked_ge50 <- 0L
n_cpg_group_all_na <- 0L

for (g in seq_len(nrow(groups))) {
  this_cancer <- groups$cancer_plot[g]
  this_type   <- groups$sample_type[g]

  idx <- which(meta_sub$cancer_plot == this_cancer & meta_sub$sample_type == this_type)
  if (length(idx) == 0) next

  x_grp <- x[idx, , drop = FALSE]

  na_frac <- colMeans(is.na(x_grp))
  can_impute <- na_frac < 0.5

  if (any(!can_impute)) {
    n_blocked_ge50 <- n_blocked_ge50 + sum(!can_impute & colSums(is.na(x_grp)) > 0)
  }

  cols_to_check <- which(can_impute & colSums(is.na(x_grp)) > 0)
  if (length(cols_to_check) == 0) next

  for (j in cols_to_check) {
    med_j <- median(x_grp[, j], na.rm = TRUE)
    if (!is.finite(med_j)) {
      n_cpg_group_all_na <- n_cpg_group_all_na + 1L
      next
    }
    na_idx_local <- which(is.na(x_grp[, j]))
    if (length(na_idx_local) > 0) {
      x[idx[na_idx_local], j] <- med_j
      n_imputed <- n_imputed + length(na_idx_local)
    }
  }
}

na_after_impute <- sum(is.na(x))
cat("STEP 5/5 - NA values imputed:", n_imputed, "\n")
cat("STEP 5/5 - CpG-group combinations blocked because NA fraction >= 50%:", n_blocked_ge50, "\n")
cat("STEP 5/5 - CpG-group combinations skipped because non-missing median could not be computed:", n_cpg_group_all_na, "\n")
cat("STEP 5/5 - Total NA values after group-wise imputation:", na_after_impute, "\n")

cat("STEP 5/5 - Removing CpGs still containing NA after group-wise imputation\n")
keep_complete <- colSums(is.na(x)) == 0
cat("STEP 5/5 - CpGs removed because NA remained:", sum(!keep_complete), "\n")
x <- x[, keep_complete, drop = FALSE]

cat("STEP 5/5 - Final matrix dimensions for t-SNE:", nrow(x), "samples x", ncol(x), "CpGs\n")

if (nrow(x) < 2) stop("Too few samples for t-SNE.")
if (ncol(x) < 2) stop("Too few CpGs for t-SNE after filtering/imputation.")

set.seed(1)

perp <- min(30, floor((nrow(x) - 1) / 3))
if (perp < 2) perp <- 2
cat("STEP 5/5 - Using perplexity:", perp, "\n")

cat("STEP 5/5 - Running Rtsne\n")
ts <- Rtsne(
  x,
  dims = 2,
  perplexity = perp,
  check_duplicates = FALSE,
  pca = TRUE,
  verbose = TRUE,
  max_iter = 1000
)

coord <- data.table(
  pan_mat_col = rownames(x),
  tSNE1 = ts$Y[,1],
  tSNE2 = ts$Y[,2]
)

coord <- merge(coord, unique(meta), by = "pan_mat_col", all.x = TRUE)

coord[, cancer := clean_txt(cancer)]
coord[, cancer_plot := sub("^TCGA-", "", cancer)]
coord[, organ_group := factor(vapply(cancer_plot, map_group, character(1), mapping = organ_map))]
coord[, tissue_type := factor(vapply(cancer_plot, map_group, character(1), mapping = tissue_map))]
coord[, cancer_plot := factor(cancer_plot)]
coord[, sample_type := factor(sample_type, levels = c("Tumor", "Healthy"))]

cat("STEP 5/5 - Example coordinates after merge:\n")
print(head(coord[, .(pan_mat_col, sample, patient, sample_type, cancer, cancer_plot, organ_group, tissue_type)]))

cat("STEP 5/5 - Writing coordinates\n")
fwrite(coord, out_tsv, sep = "\t")

cat("STEP 5/5 - Saving Rtsne object\n")
saveRDS(ts, out_rds)

missing_plot_cols <- setdiff(levels(coord$cancer_plot), names(cancer_cols))
if (length(missing_plot_cols) > 0) {
  for (cc in missing_plot_cols) cancer_cols[cc] <- "#7f7f7f"
}
pal_cancer <- cancer_cols[levels(coord$cancer_plot)]

pal_organ <- Polychrome::createPalette(
  nlevels(coord$organ_group),
  seedcolors = c("#ff00c3ff", "#815801ff", "#24099eff", "#55ffd2ff")
)
names(pal_organ) <- levels(coord$organ_group)

pal_tissue <- Polychrome::createPalette(
  nlevels(coord$tissue_type),
  seedcolors = c("#2bff00ff", "#d9ff00ff", "#9e8dffff", "#004a36ff")
)
names(pal_tissue) <- levels(coord$tissue_type)

base_theme <- theme_bw(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 20),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 15),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 13),
    legend.key.size = unit(0.5, "cm")
  )

cat("STEP 5/5 - Plotting PDF\n")

p_cancer <- ggplot(coord, aes(x = tSNE1, y = tSNE2, color = cancer_plot, shape = sample_type)) +
  geom_point(size = 1.5, alpha = 0.85) +
  scale_color_manual(values = pal_cancer, drop = FALSE, na.value = "grey70") +
  labs(
    title = "Cancer type",
    x = "tSNE1",
    y = "tSNE2"
  ) +
  base_theme

p_organ <- ggplot(coord, aes(x = tSNE1, y = tSNE2, color = organ_group)) +
  geom_point(size = 1.5, alpha = 0.85) +
  scale_color_manual(values = pal_organ, drop = FALSE, na.value = "grey70") +
  labs(
    title = "Organ group",
    x = "tSNE1",
    y = "tSNE2"
  ) +
  base_theme

p_tissue <- ggplot(coord, aes(x = tSNE1, y = tSNE2, color = tissue_type)) +
  geom_point(size = 1.5, alpha = 0.85) +
  scale_color_manual(values = pal_tissue, drop = FALSE, na.value = "grey70") +
  labs(
    title = "Tissue type",
    x = "tSNE1",
    y = "tSNE2"
  ) +
  base_theme

pdf(out_pdf, width = 12, height = 24)
print(
  (p_cancer / p_organ / p_tissue) +
    plot_annotation(
      title = "t-SNE of TCGA HM450 methylation",
      subtitle = paste0(
        "samples cohort | top ", nrow(top_dt),
        " variable CpGs | n=", nrow(coord),
        " samples | p=", ncol(x), " CpGs"
      ),
      theme = theme(
        plot.title = element_text(size = 24, face = "bold"),
        plot.subtitle = element_text(size = 18)
      )
    )
)
dev.off()

cat("STEP 5/5 - DONE\n")
cat("Wrote coordinates:", out_tsv, "\n")
cat("Wrote plot:", out_pdf, "\n")
cat("Wrote RDS:", out_rds, "\n")
EOF













##########################################################################################################
# step 6 plot the interactive version of the tSNE with plotly
##############################################################################################################
echo "Running Step 6: interactive t-SNE HTML"

Rscript - <<'EOF'
suppressPackageStartupMessages({
  library(data.table)
  library(plotly)
  library(htmlwidgets)
})

cat("STEP 6/6 - Starting interactive t-SNE export\n")

coord_file <- "./results/methylation/tSNE_methylation/tsne_3d4d_noOV/tSNE_methylation_3d4d_noOV_top10k_coordinates.tsv"
color_file <- "./results/multi_omics/cancer_color_order_with_defined_colours.tsv"
out_html   <- "./results/methylation/tSNE_methylation/tsne_3d4d_noOV/tSNE_methylation_3d4d_noOV_top10k_interactive.html"

clean_txt <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub('^"|"$', "", x)
  x <- gsub("^'|'$", "", x)
  x
}

cat("STEP 6/6 - Reading coordinates\n")
dt <- fread(coord_file)

cat("STEP 6/6 - Coordinates rows:", nrow(dt), "\n")

if (!"cancer" %in% colnames(dt)) {
  if ("file" %in% colnames(dt)) {
    dt[, file_base := basename(file)]
    dt[, cancer := sub("^HM450_(TCGA-[A-Z0-9]+)-TCGA-.*$", "\\1", file_base)]
  } else {
    stop("Neither 'cancer' nor 'file' column found in coordinates file.")
  }
}

dt[, cancer := clean_txt(cancer)]
dt[, cancer_plot := sub("^TCGA-", "", cancer)]

cat("STEP 6/6 - Reading color file\n")
col_dt <- fread(color_file, header = FALSE)
col_dt <- col_dt[, 1:2]
setnames(col_dt, c("cancer_raw", "colour"))
col_dt[, cancer_raw := clean_txt(cancer_raw)]
col_dt[, colour := clean_txt(colour)]
col_dt[, cancer_plot := sub("^TCGA-", "", cancer_raw)]
col_dt <- unique(col_dt[, .(cancer_plot, colour)])

pal <- setNames(col_dt$colour, col_dt$cancer_plot)

# fallback for missing palette entries
missing_cols <- setdiff(unique(dt$cancer_plot), names(pal))
if (length(missing_cols) > 0) {
  cat("STEP 6/6 - Cancers missing from palette:", length(missing_cols), "\n")
  print(missing_cols)
  for (cc in missing_cols) pal[cc] <- "#7f7f7f"
}

# safer labels
for (nm in c("sample","patient","sample_type","file","pan_mat_col")) {
  if (!nm %in% colnames(dt)) dt[, (nm) := NA_character_]
}

dt[, hover_text := paste0(
  "Sample: ", sample,
  "<br>Patient: ", patient,
  "<br>Cancer: ", cancer_plot,
  "<br>Sample type: ", sample_type,
  "<br>pan_mat column: ", pan_mat_col,
  "<br>File: ", file
)]

cat("STEP 6/6 - Building interactive plot\n")

p <- plot_ly(
  data = dt,
  x = ~tSNE1,
  y = ~tSNE2,
  type = "scatter",
  mode = "markers",
  color = ~cancer_plot,
  colors = pal,
  symbol = ~sample_type,
  text = ~hover_text,
  hoverinfo = "text",
  marker = list(size = 7, opacity = 0.85)
) %>%
  layout(
    title = list(text = "Interactive t-SNE of TCGA HM450 methylation"),
    xaxis = list(title = "tSNE1"),
    yaxis = list(title = "tSNE2"),
    legend = list(title = list(text = "Cancer / sample type"))
  )

cat("STEP 6/6 - Writing HTML\n")
saveWidget(p, out_html, selfcontained = TRUE)

cat("STEP 6/6 - DONE\n")
cat("Wrote:", out_html, "\n")
EOF

#######################################################
# Redo the TSNE but using the Cg that are in motifs :
#######################################################
echo "Running motif-based tSNE for BANP and NRF1"

Rscript - <<'EOF'
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(Rtsne)
})

cat("STEP tSNE - Starting motif-based tSNE for BANP and NRF1\n")

# =========================
# INPUTS
# =========================
sample_file <- "./results/methylation/tSNE_methylation/tsne_3d4d_noOV/methylation_samples_3d4d_noOV_in_panmat.tsv"
pan_mat_rds <- "./methylation/TCGA-all_cancer_pan_cancer_matrix.RData"
color_file  <- "./results/multi_omics/cancer_color_order_with_defined_colours.tsv"

tf_files <- c(
  BANP = "./motifs/overlaps/intersected_motifs2mm_HM450/BANP_mm0to2_noCGmm_probeXmotif.tsv",
  NRF1 = "./motifs/overlaps/intersected_motifs2mm_HM450/NRF1_mm0to2_noCGmm_probeXmotif.tsv"
)

out_root <- "./results/methylation/tSNE_methylation_motifCpGs_allCg"

# =========================
# HELPERS
# =========================
clean_txt <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub('^"|"$', "", x)
  x <- gsub("^'|'$", "", x)
  x
}

# =========================
# LOAD SHARED INPUTS
# =========================
cat("STEP tSNE - Reading matched sample file\n")
meta <- fread(sample_file)

if (!"pan_mat_col" %in% colnames(meta)) {
  stop("Column 'pan_mat_col' not found in sample file.")
}

cat("STEP tSNE - Using all sample types:", nrow(meta), "\n")

meta[, file_base := basename(file)]

if (!"cancer" %in% colnames(meta)) {
  meta[, cancer := sub("^HM450_(TCGA-[A-Z0-9]+)-TCGA-.*$", "\\1", file_base)]
}

meta[, cancer := clean_txt(cancer)]
meta[, cancer_plot := sub("^TCGA-", "", cancer)]

cat("STEP tSNE - Reading cancer color file\n")
col_dt <- fread(color_file, header = FALSE)
if (ncol(col_dt) < 2) stop("color_file must have at least 2 columns")
col_dt <- col_dt[, 1:2]
setnames(col_dt, c("cancer_raw", "colour"))
col_dt[, cancer_raw := clean_txt(cancer_raw)]
col_dt[, colour := clean_txt(colour)]
col_dt[, cancer_plot := sub("^TCGA-", "", cancer_raw)]
col_dt <- unique(col_dt[, .(cancer_plot, colour)])
cancer_cols <- setNames(col_dt$colour, col_dt$cancer_plot)

cat("STEP tSNE - Loading pan_mat\n")
load(pan_mat_rds)
if (!exists("pan_mat")) stop("pan_mat not found in RData")
cat("STEP tSNE - pan_mat dimensions:", nrow(pan_mat), "CpGs x", ncol(pan_mat), "samples\n")

samples <- unique(as.character(meta$pan_mat_col))
missing_samples <- setdiff(samples, colnames(pan_mat))
cat("STEP tSNE - Missing samples in pan_mat:", length(missing_samples), "\n")
if (length(missing_samples) > 0) {
  print(head(missing_samples, 20))
  stop("Some samples are missing from pan_mat.")
}

# =========================
# LOOP OVER TFS
# =========================
for (tf_name in names(tf_files)) {
  cat("\n============================\n")
  cat("STEP tSNE - Processing TF:", tf_name, "\n")
  cat("============================\n")

  motif_cpg_file <- tf_files[[tf_name]]
  out_dir <- file.path(out_root, tf_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  out_pdf <- file.path(out_dir, paste0("tSNE_", tf_name, "_all_motif_CpGs_3d4d_noOV.pdf"))
  out_tsv <- file.path(out_dir, paste0("tSNE_", tf_name, "_all_motif_CpGs_3d4d_noOV_coordinates.tsv"))
  out_rds <- file.path(out_dir, paste0("tSNE_", tf_name, "_all_motif_CpGs_3d4d_noOV_rtsne.rds"))
  out_cpg <- file.path(out_dir, paste0("CpGs_used_for_tSNE_", tf_name, "_all_motif_CpGs_3d4d_noOV.tsv"))

  cat("STEP tSNE - Reading motif CpG file:", motif_cpg_file, "\n")
  mot <- fread(motif_cpg_file, header = FALSE)

  cat("STEP tSNE - motif file columns:\n")
  print(colnames(mot))

  if (!"V4" %in% colnames(mot)) {
    stop(paste("Column V4 not found in", motif_cpg_file, "- cannot extract probes"))
  }

  mot[, motif_probe := as.character(V4)]
  motif_cpgs <- unique(mot$motif_probe)

  cat("STEP tSNE - Number of unique motif CpGs in", tf_name, ":", length(motif_cpgs), "\n")

  motif_cpgs_in_mat <- intersect(motif_cpgs, rownames(pan_mat))
  cat("STEP tSNE - Motif CpGs found in pan_mat:", length(motif_cpgs_in_mat), "\n")
  if (length(motif_cpgs_in_mat) < 2) {
    cat("STEP tSNE - Skipping", tf_name, "- too few motif CpGs found in pan_mat\n")
    next
  }

  fwrite(data.table(CpG = motif_cpgs_in_mat), out_cpg, sep = "\t")

  cat("STEP tSNE - Subsetting pan_mat\n")
  mat <- pan_mat[motif_cpgs_in_mat, samples, drop = FALSE]
  cat("STEP tSNE - Subset matrix dimensions:", nrow(mat), "CpGs x", ncol(mat), "samples\n")

  cat("STEP tSNE - Transposing to samples x CpGs\n")
  x <- t(mat)

  cat("STEP tSNE - Removing CpGs with any NA across kept samples\n")
  keep_complete <- colSums(is.na(x)) == 0
  cat("STEP tSNE - CpGs removed due to NA:", sum(!keep_complete), "\n")
  x <- x[, keep_complete, drop = FALSE]

  cat("STEP tSNE - Removing zero-variance CpGs\n")
  vars <- apply(x, 2, var)
  keep_var <- is.finite(vars) & vars > 0
  cat("STEP tSNE - CpGs removed due to zero/non-finite variance:", sum(!keep_var), "\n")
  x <- x[, keep_var, drop = FALSE]

  cat("STEP tSNE - Final matrix dimensions:", nrow(x), "samples x", ncol(x), "CpGs\n")

  if (nrow(x) < 3) {
    cat("STEP tSNE - Skipping", tf_name, "- too few samples after filtering\n")
    next
  }
  if (ncol(x) < 2) {
    cat("STEP tSNE - Skipping", tf_name, "- too few CpGs after filtering\n")
    next
  }

  set.seed(1)

  perp <- min(30, floor((nrow(x) - 1) / 3))
  if (perp < 2) perp <- 2
  cat("STEP tSNE - Using perplexity:", perp, "\n")

  cat("STEP tSNE - Running Rtsne for", tf_name, "\n")
  ts <- Rtsne(
    x,
    dims = 2,
    perplexity = perp,
    check_duplicates = FALSE,
    pca = TRUE,
    verbose = TRUE,
    max_iter = 1000
  )

  coord <- data.table(
    pan_mat_col = rownames(x),
    tSNE1 = ts$Y[, 1],
    tSNE2 = ts$Y[, 2]
  )

  coord <- merge(coord, unique(meta), by = "pan_mat_col", all.x = TRUE)
  coord[, cancer_plot := clean_txt(cancer_plot)]

  cat("STEP tSNE - Writing coordinates for", tf_name, "\n")
  fwrite(coord, out_tsv, sep = "\t")

  cat("STEP tSNE - Saving Rtsne object for", tf_name, "\n")
  saveRDS(ts, out_rds)

  cat("STEP tSNE - Plotting PDF for", tf_name, "\n")
  p <- ggplot(coord, aes(x = tSNE1, y = tSNE2, color = cancer_plot, shape = sample_type)) +
    geom_point(size = 2.2, alpha = 0.9) +
    scale_color_manual(values = cancer_cols, drop = FALSE, na.value = "grey70") +
    theme_bw(base_size = 12) +
    labs(
      title = paste0("t-SNE of HM450 methylation - ", tf_name, " motif CpGs"),
      subtitle = paste0(
        "3d+4d_noOV | all sample types | all CpGs in motifs | n=", nrow(coord),
        " samples | p=", ncol(x), " CpGs"
      ),
      x = "tSNE1",
      y = "tSNE2",
      color = "Cancer",
      shape = "Sample type"
    )

  ggsave(out_pdf, p, width = 10, height = 7)

  cat("STEP tSNE - DONE for", tf_name, "\n")
  cat("  Coordinates:", out_tsv, "\n")
  cat("  Plot:", out_pdf, "\n")
  cat("  RDS:", out_rds, "\n")
  cat("  CpGs used:", out_cpg, "\n")
}

cat("\nSTEP tSNE - All TF runs completed\n")
EOF
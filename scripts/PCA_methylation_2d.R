mkdir -p ./results/methylation/pca_2d_noOV

# 1) Get the list of methylation files for the 2d noOV cohort (tumor 0X vs normal 1X): --> ./results/methylation/pca_2d_noOV/methylation_files_2d_noOV_0X_1X.tsv
Rscript -e '
suppressPackageStartupMessages({
  library(data.table)
})

meth_dir    <- "./methylation/filtered_methylation"
cohort_file <- "./results/multi_omics/samples_2d_noOV.tsv"
out_file    <- "./results/methylation/pca_2d_noOV/methylation_files_2d_noOV_0X_1X.tsv"

extract_sample_barcode <- function(x) {
  m <- regexpr("TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}[A-Z]", x)
  ifelse(m > 0, regmatches(x, m), NA_character_)
}

extract_patient_id <- function(sample_barcode) {
  sub("^(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}).*$", "\\\\1", sample_barcode)
}

extract_sample_code <- function(sample_barcode) {
  sub("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-([0-9]{2}).*$", "\\\\1", sample_barcode)
}

extract_cancer <- function(filename) {
  x <- basename(filename)
  x <- sub("^HM450_TCGA-", "", x)
  sub("-.*$", "", x)
}

# 2d noOV patients: first column, no header
cohort <- fread(cohort_file, header = FALSE)
patients_2d <- unique(as.character(cohort[[1]]))

cat(">>> 2d noOV patients:", length(patients_2d), "\\n")

files <- list.files(
  meth_dir,
  pattern = "_annotated_methylation_filtered\\\\.bed\\\\.gz$",
  full.names = TRUE
)

dt <- data.table(file = files)
dt[, filename := basename(file)]
dt[, cancer := extract_cancer(filename)]
dt[, sample_barcode := extract_sample_barcode(filename)]
dt[, patient_id := extract_patient_id(sample_barcode)]
dt[, sample_code := extract_sample_code(sample_barcode)]

dt[, sample_type := fifelse(sample_code == "01", "Tumor",
                     fifelse(sample_code == "11", "Healthy", NA_character_))]

dt_keep <- dt[
  patient_id %in% patients_2d &
  sample_code %in% sprintf("%02d", 1:19)
]

cat(">>> Total methylation files:", nrow(dt), "\\n")
cat(">>> Kept 2d noOV 01/11 files:", nrow(dt_keep), "\\n")
cat(">>> Kept cancers:", length(unique(dt_keep$cancer)), "\\n")

cat("\\n>>> Sample count by type:\\n")
print(dt_keep[, .N, by = sample_type][order(sample_type)])

cat("\\n>>> Sample count by cancer and type:\\n")
print(dt_keep[, .N, by = .(cancer, sample_type)][order(cancer, sample_type)])

fwrite(dt_keep, out_file, sep = "\\t")

cat("\\n>>> Saved:", out_file, "\\n")
'

# 2) Create the methylation matrix for the 2d noOV cohort (tumor 0X vs normal 1X) and save as RData: --> ./results/methylation/pca_2d_noOV/TCGA_2d_noOV_0X_1X_pan_cancer_matrix.RData
mkdir -p ./results/methylation/pca_2d_noOV

Rscript -e '
suppressPackageStartupMessages({
  library(data.table)
})

files_table <- "./results/methylation/pca_2d_noOV/methylation_files_2d_noOV_0X_1X.tsv"
out_rdata   <- "./results/methylation/pca_2d_noOV/TCGA_2d_noOV_0X_1X_pan_cancer_matrix.RData"

dt_files <- fread(files_table)

cat(">>> Number of methylation files:", nrow(dt_files), "\\n")

files <- dt_files$file
sample_names <- dt_files$sample_barcode

# read first file to define CpG universe
cat(">>> Reading first file to define CpG order...\\n")

dt0 <- fread(
  cmd = paste("zcat", shQuote(files[1])),
  header = FALSE,
  select = c(4, 5),
  showProgress = FALSE
)

cpgs <- as.character(dt0[[1]])

cat(">>> CpGs in first file:", length(cpgs), "\\n")

pan_mat <- matrix(
  NA_real_,
  nrow = length(cpgs),
  ncol = length(files),
  dimnames = list(cpgs, sample_names)
)

pan_mat[, 1] <- as.numeric(dt0[[2]])

cat(">>> Filling matrix...\\n")

for (i in 2:length(files)) {
  if (i %% 100 == 0) {
    cat(">>> Processed", i, "of", length(files), "files\\n")
  }

  dt <- fread(
    cmd = paste("zcat", shQuote(files[i])),
    header = FALSE,
    select = c(4, 5),
    showProgress = FALSE
  )

  vals <- setNames(as.numeric(dt[[2]]), as.character(dt[[1]]))
  pan_mat[, i] <- vals[cpgs]
}

cat(">>> Final matrix dimensions:", paste(dim(pan_mat), collapse = " x "), "\\n")

save(pan_mat, file = out_rdata)

cat(">>> Saved:", out_rdata, "\\n")
'

# 3) Select most variable CpGs across the 2d noOV cohort (tumor 0X vs normal 1X) and save as RData: --> ./results/methylation/pca_2d_noOV/TCGA_2d_noOV_0X_1X_top10k_most_variable.RData 

Rscript -e '
library(matrixStats)

# load the 2d_noOV 0X/1X pan-cancer methylation matrix
load("./results/methylation/pca_2d_noOV/TCGA_2d_noOV_0X_1X_pan_cancer_matrix.RData")
# loads object: pan_mat

# select the 10,000 most variable CpGs across the 2d_noOV 0X/1X samples
top_idx <- order(rowVars(pan_mat, na.rm = TRUE), decreasing = TRUE)[seq_len(10000)]
beta_var <- pan_mat[top_idx, , drop = FALSE]

# save results
save(beta_var, file = "./results/methylation/pca_2d_noOV/TCGA_2d_noOV_0X_1X_top10k_most_variable.RData")

writeLines(
  rownames(beta_var),
  "./results/methylation/pca_2d_noOV/top10000_variable_CpGs_2d_noOV_0X_1X.txt"
)
'

# 4) Perform PCA on the 2d noOV 0X/1X top 10k most variable CpGs 
# 4) Perform PCA on the 2d noOV 0X/1X top 10k most variable CpGs
# Color = cancer type using predefined cancer color file
# Shape = Tumor vs Healthy

Rscript -e '
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# -----------------------------
# Inputs
# -----------------------------
in_rdata <- "./results/methylation/pca_2d_noOV/TCGA_2d_noOV_0X_1X_top10k_most_variable.RData"

files_table <- "./results/methylation/pca_2d_noOV/methylation_files_2d_noOV_0X_1X.tsv"

color_file <- "./results/multi_omics/cancer_color_order_with_defined_colours.tsv"

out_coords <- "./results/methylation/pca_2d_noOV/PCA_coordinates_2d_noOV_0X_1X_top10k.tsv"
out_pdf    <- "./results/methylation/pca_2d_noOV/PCA_2d_noOV_0X_1X_top10k.pdf"
out_rds    <- "./results/methylation/pca_2d_noOV/PCA_2d_noOV_0X_1X_top10k_prcomp.rds"

# -----------------------------
# Load top variable CpG matrix
# -----------------------------
cat(">>> Loading beta_var matrix...\\n")
load(in_rdata)  # loads beta_var

cat(">>> beta_var dimensions:", paste(dim(beta_var), collapse = " x "), "\\n")
cat(">>> Rows = CpGs, columns = samples\\n")

# -----------------------------
# Remove CpGs with all NA
# -----------------------------
cat(">>> Removing CpGs with all missing values...\\n")

keep_not_all_na <- rowSums(!is.na(beta_var)) > 0
beta_var <- beta_var[keep_not_all_na, , drop = FALSE]

cat(">>> CpGs after all-NA filter:", nrow(beta_var), "\\n")

# -----------------------------
# Impute missing values by CpG mean
# -----------------------------
cat(">>> Imputing missing values by CpG mean...\\n")

row_means <- rowMeans(beta_var, na.rm = TRUE)

na_idx <- which(is.na(beta_var), arr.ind = TRUE)

if (nrow(na_idx) > 0) {
  beta_var[na_idx] <- row_means[na_idx[, 1]]
}

cat(">>> Number of imputed values:", nrow(na_idx), "\\n")

# -----------------------------
# Remove zero-variance CpGs
# -----------------------------
cat(">>> Removing zero-variance CpGs...\\n")

cpg_sd <- apply(beta_var, 1, sd, na.rm = TRUE)
valid <- is.finite(cpg_sd) & cpg_sd > 0

beta_var <- beta_var[valid, , drop = FALSE]

cat(">>> CpGs retained for PCA:", nrow(beta_var), "\\n")

# -----------------------------
# Transpose to samples x CpGs
# -----------------------------
cat(">>> Transposing matrix for PCA...\\n")

mat_pca <- t(beta_var)

cat(">>> PCA input dimensions:", paste(dim(mat_pca), collapse = " x "), "\\n")
cat(">>> Rows = samples, columns = CpGs\\n")

# -----------------------------
# Run PCA
# -----------------------------
cat(">>> Running PCA...\\n")

pca <- prcomp(
  mat_pca,
  center = TRUE,
  scale. = TRUE
)

saveRDS(pca, out_rds)

cat(">>> Saved PCA object:", out_rds, "\\n")

# -----------------------------
# Variance explained
# -----------------------------
var_explained <- pca$sdev^2 / sum(pca$sdev^2)

cat(">>> PC1 variance:", round(var_explained[1] * 100, 2), "%\\n")
cat(">>> PC2 variance:", round(var_explained[2] * 100, 2), "%\\n")

# -----------------------------
# Build PCA coordinate table
# -----------------------------
pca_df <- data.table(
  sample_barcode = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  PC3 = pca$x[, 3],
  PC4 = pca$x[, 4]
)

# -----------------------------
# Load sample metadata
# -----------------------------
cat(">>> Loading sample metadata...\\n")

meta <- fread(files_table)

# IMPORTANT:
# fread may convert 01 to 1.
# This restores leading zeros: 1 -> 01, 11 -> 11
meta[, sample_code := sprintf("%02d", as.integer(sample_code))]

# Define 0X / 1X group
meta[, sample_group := fifelse(
  sample_code %in% sprintf("%02d", 1:9), "0X",
  fifelse(sample_code %in% sprintf("%02d", 10:19), "1X", NA_character_)
)]

# Define Tumor / Healthy
meta[, sample_type_clean := fifelse(
  sample_code %in% sprintf("%02d", 1:9), "Tumor",
  fifelse(sample_code %in% sprintf("%02d", 10:19), "Healthy", NA_character_)
)]

cat(">>> Sample type table after cleaning sample_code:\\n")
print(table(meta$sample_code, meta$sample_type_clean, useNA = "ifany"))

# Keep useful metadata columns
meta_small <- unique(meta[, .(
  sample_barcode,
  patient_id,
  cancer,
  sample_code,
  sample_group,
  sample_type_clean
)])

# Merge PCA coordinates with metadata
pca_df <- merge(
  pca_df,
  meta_small,
  by = "sample_barcode",
  all.x = TRUE
)

cat(">>> PCA samples with missing cancer metadata:", sum(is.na(pca_df$cancer)), "\\n")
cat(">>> PCA samples with missing sample type:", sum(is.na(pca_df$sample_type_clean)), "\\n")

fwrite(pca_df, out_coords, sep = "\\t")

cat(">>> Saved PCA coordinates:", out_coords, "\\n")

# -----------------------------
# Load predefined cancer colors
# -----------------------------
cat(">>> Loading predefined cancer color file...\\n")

color_dt <- fread(color_file, header = FALSE)

# The file is expected to have 2 columns:
# cancer code and color
if (ncol(color_dt) < 2) {
  stop("Color file must contain at least two columns: cancer and color.")
}

color_dt <- color_dt[, .(cancer = V1, color = V2)]

# Remove possible TCGA- prefix
color_dt[, cancer := gsub("^TCGA-", "", cancer)]

# Make named color vector
cancer_colors <- setNames(color_dt$color, color_dt$cancer)

# Keep only cancers present in PCA
present_cancers <- sort(unique(pca_df$cancer))
missing_colors <- setdiff(present_cancers, names(cancer_colors))

if (length(missing_colors) > 0) {
  cat(">>> WARNING: Missing colors for cancers:\\n")
  print(missing_colors)
}

# Keep colors only for cancers present in the PCA
cancer_colors <- cancer_colors[names(cancer_colors) %in% present_cancers]

cat(">>> Cancers in PCA:", length(present_cancers), "\\n")
cat(">>> Colors available:", length(cancer_colors), "\\n")

# -----------------------------
# Plot PCA
# -----------------------------
cat(">>> Plotting PCA...\\n")

p <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    color = cancer,
    shape = sample_type_clean
  )
) +
  geom_point(size = 2, alpha = 0.85) +
  scale_color_manual(
    values = cancer_colors,
    na.translate = FALSE
  ) +
  scale_shape_manual(
    values = c(
      Tumor = 16,
      Healthy = 17
    ),
    na.translate = FALSE
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.45, "cm")
  ) +
  labs(
    title = "PCA of DNA methylation profiles in the TCGA 2d_noOV cohort",
    subtitle = "Top 10,000 most variable CpG sites across 0X and 1X samples",
    x = paste0("PC1 (", round(var_explained[1] * 100, 1), "% variance)"),
    y = paste0("PC2 (", round(var_explained[2] * 100, 1), "% variance)"),
    color = "Cancer type",
    shape = "Sample type"
  )

pdf(out_pdf, width = 13, height = 8)
print(p)
dev.off()

cat(">>> Saved PCA plot:", out_pdf, "\\n")
cat(">>> DONE\\n")
'

# 5) tSNE
# 4) Perform t-SNE on the 2d noOV 0X/1X top 10k most variable CpGs
# Color = cancer type using predefined cancer color file
# Shape = Tumor vs Healthy

Rscript -e '
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(Rtsne)
})

# -----------------------------
# Inputs
# -----------------------------
in_rdata <- "./results/methylation/pca_2d_noOV/TCGA_2d_noOV_0X_1X_top10k_most_variable.RData"

files_table <- "./results/methylation/pca_2d_noOV/methylation_files_2d_noOV_0X_1X.tsv"

color_file <- "./results/multi_omics/cancer_color_order_with_defined_colours.tsv"

out_coords <- "./results/methylation/pca_2d_noOV/tSNE_coordinates_2d_noOV_0X_1X_top10k.tsv"
out_pdf    <- "./results/methylation/pca_2d_noOV/tSNE_2d_noOV_0X_1X_top10k.pdf"
out_rds    <- "./results/methylation/pca_2d_noOV/tSNE_2d_noOV_0X_1X_top10k_rtsne.rds"

# -----------------------------
# t-SNE parameters
# -----------------------------
set.seed(123)

initial_dims <- 30
max_iter     <- 1000

# -----------------------------
# Load top variable CpG matrix
# -----------------------------
cat(">>> Loading beta_var matrix...\\n")
load(in_rdata)  # loads beta_var

cat(">>> beta_var dimensions:", paste(dim(beta_var), collapse = " x "), "\\n")
cat(">>> Rows = CpGs, columns = samples\\n")

# -----------------------------
# Remove CpGs with all NA
# -----------------------------
cat(">>> Removing CpGs with all missing values...\\n")

keep_not_all_na <- rowSums(!is.na(beta_var)) > 0
beta_var <- beta_var[keep_not_all_na, , drop = FALSE]

cat(">>> CpGs after all-NA filter:", nrow(beta_var), "\\n")

# -----------------------------
# Impute missing values by CpG mean
# -----------------------------
cat(">>> Imputing missing values by CpG mean...\\n")

row_means <- rowMeans(beta_var, na.rm = TRUE)
na_idx <- which(is.na(beta_var), arr.ind = TRUE)

if (nrow(na_idx) > 0) {
  beta_var[na_idx] <- row_means[na_idx[, 1]]
}

cat(">>> Number of imputed values:", nrow(na_idx), "\\n")

# -----------------------------
# Remove zero-variance CpGs
# -----------------------------
cat(">>> Removing zero-variance CpGs...\\n")

cpg_sd <- apply(beta_var, 1, sd, na.rm = TRUE)
valid <- is.finite(cpg_sd) & cpg_sd > 0

beta_var <- beta_var[valid, , drop = FALSE]

cat(">>> CpGs retained for t-SNE:", nrow(beta_var), "\\n")

# -----------------------------
# Transpose to samples x CpGs
# -----------------------------
cat(">>> Transposing matrix for t-SNE...\\n")

mat_tsne <- t(beta_var)

cat(">>> t-SNE input dimensions:", paste(dim(mat_tsne), collapse = " x "), "\\n")
cat(">>> Rows = samples, columns = CpGs\\n")

# -----------------------------
# Scale CpGs
# -----------------------------
cat(">>> Scaling CpGs...\\n")

mat_tsne <- scale(mat_tsne)

# Remove columns that became non-finite after scaling
finite_cols <- apply(mat_tsne, 2, function(x) all(is.finite(x)))
mat_tsne <- mat_tsne[, finite_cols, drop = FALSE]

cat(">>> CpGs retained after scaling:", ncol(mat_tsne), "\\n")

# -----------------------------
# Choose valid perplexity
# -----------------------------
n_samples <- nrow(mat_tsne)

perplexity <- min(30, floor((n_samples - 1) / 3))

if (perplexity < 2) {
  stop("Not enough samples for t-SNE. Need more samples to use perplexity >= 2.")
}

cat(">>> Number of samples:", n_samples, "\\n")
cat(">>> Using perplexity:", perplexity, "\\n")
cat(">>> Using initial_dims:", initial_dims, "\\n")
cat(">>> Using max_iter:", max_iter, "\\n")

# -----------------------------
# Run t-SNE
# -----------------------------
cat(">>> Running t-SNE...\\n")

tsne <- Rtsne(
  mat_tsne,
  dims = 2,
  perplexity = perplexity,
  initial_dims = min(initial_dims, ncol(mat_tsne), n_samples - 1),
  max_iter = max_iter,
  pca = TRUE,
  check_duplicates = FALSE,
  verbose = TRUE
)

saveRDS(tsne, out_rds)

cat(">>> Saved t-SNE object:", out_rds, "\\n")

# -----------------------------
# Build t-SNE coordinate table
# -----------------------------
tsne_df <- data.table(
  sample_barcode = rownames(mat_tsne),
  tSNE1 = tsne$Y[, 1],
  tSNE2 = tsne$Y[, 2]
)

# -----------------------------
# Load sample metadata
# -----------------------------
cat(">>> Loading sample metadata...\\n")

meta <- fread(files_table)

# fread may convert 01 to 1, so restore leading zeros
meta[, sample_code := sprintf("%02d", as.integer(sample_code))]

# Define 0X / 1X group
meta[, sample_group := fifelse(
  sample_code %in% sprintf("%02d", 1:9), "0X",
  fifelse(sample_code %in% sprintf("%02d", 10:19), "1X", NA_character_)
)]

# Define Tumor / Healthy
meta[, sample_type_clean := fifelse(
  sample_code %in% sprintf("%02d", 1:9), "Tumor",
  fifelse(sample_code %in% sprintf("%02d", 10:19), "Healthy", NA_character_)
)]

cat(">>> Sample type table after cleaning sample_code:\\n")
print(table(meta$sample_code, meta$sample_type_clean, useNA = "ifany"))

# Keep useful metadata columns
meta_small <- unique(meta[, .(
  sample_barcode,
  patient_id,
  cancer,
  sample_code,
  sample_group,
  sample_type_clean
)])

# Merge t-SNE coordinates with metadata
tsne_df <- merge(
  tsne_df,
  meta_small,
  by = "sample_barcode",
  all.x = TRUE
)

cat(">>> t-SNE samples with missing cancer metadata:", sum(is.na(tsne_df$cancer)), "\\n")
cat(">>> t-SNE samples with missing sample type:", sum(is.na(tsne_df$sample_type_clean)), "\\n")

fwrite(tsne_df, out_coords, sep = "\\t")

cat(">>> Saved t-SNE coordinates:", out_coords, "\\n")

# -----------------------------
# Load predefined cancer colors
# -----------------------------
cat(">>> Loading predefined cancer color file...\\n")

color_dt <- fread(color_file, header = FALSE)

if (ncol(color_dt) < 2) {
  stop("Color file must contain at least two columns: cancer and color.")
}

color_dt <- color_dt[, .(cancer = V1, color = V2)]

# Remove possible TCGA- prefix
color_dt[, cancer := gsub("^TCGA-", "", cancer)]

# Make named color vector
cancer_colors <- setNames(color_dt$color, color_dt$cancer)

# Keep only cancers present in t-SNE
present_cancers <- sort(unique(tsne_df$cancer))
missing_colors <- setdiff(present_cancers, names(cancer_colors))

if (length(missing_colors) > 0) {
  cat(">>> WARNING: Missing colors for cancers:\\n")
  print(missing_colors)
}

cancer_colors <- cancer_colors[names(cancer_colors) %in% present_cancers]

cat(">>> Cancers in t-SNE:", length(present_cancers), "\\n")
cat(">>> Colors available:", length(cancer_colors), "\\n")

# -----------------------------
# Plot t-SNE
# -----------------------------
cat(">>> Plotting t-SNE...\\n")

p <- ggplot(
  tsne_df,
  aes(
    x = tSNE1,
    y = tSNE2,
    color = cancer,
    shape = sample_type_clean
  )
) +
  geom_point(size = 2, alpha = 0.85) +
  scale_color_manual(
    values = cancer_colors,
    na.translate = FALSE
  ) +
  scale_shape_manual(
    values = c(
      Tumor = 16,
      Healthy = 17
    ),
    na.translate = FALSE
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.45, "cm")
  ) +
  labs(
    title = "t-SNE of DNA methylation profiles in the TCGA 2d_noOV cohort",
    subtitle = paste0(
      "Top 10,000 most variable CpG sites across 0X and 1X samples; perplexity = ",
      perplexity
    ),
    x = "t-SNE 1",
    y = "t-SNE 2",
    color = "Cancer type",
    shape = "Sample type"
  )

pdf(out_pdf, width = 13, height = 8)
print(p)
dev.off()

cat(">>> Saved t-SNE plot:", out_pdf, "\\n")
cat(">>> DONE\\n")
'
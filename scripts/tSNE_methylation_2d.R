###############################################################################
# PCA + t-SNE OF DNA METHYLATION PROFILES
# TCGA noOV_noCHOL cohort
#
# Cohort:
# - uses patients from ./results/multi_omics/samples_2d_noOV_noCHOL.tsv
# - explicitly removes OV and CHOL
# - keeps 0X tumor samples and 1X healthy/normal samples
#
# NA handling:
# - CpGs with >=50% NA across samples are removed
# - remaining NA values are imputed using CpG-wise median
#
# Dimensionality reduction:
# - top 10,000 most variable CpGs after filtering/imputation
# - PCA and t-SNE
###############################################################################
mkdir -p ./results/methylation/dimred_noOV_noCHOL

###############################################################################
# 1) Build methylation file table for TCGA noOV_noCHOL cohort
###############################################################################

mkdir -p ./results/methylation/dimred_noOV_noCHOL

Rscript -e '
suppressPackageStartupMessages({
  library(data.table)
})

meth_dir    <- "./methylation/filtered_methylation"
cohort_file <- "./results/multi_omics/samples_2d_noOV_noCHOL.tsv"
out_file    <- "./results/methylation/dimred_noOV_noCHOL/methylation_files_noOV_noCHOL_0X_1X.tsv"

extract_sample_barcode <- function(x) {
  m <- regexpr("TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}[A-Z]", x)
  out <- rep(NA_character_, length(x))
  ok <- m > 0
  out[ok] <- regmatches(x, m)[ok]
  out
}

extract_patient_id <- function(sample_barcode) {
  out <- rep(NA_character_, length(sample_barcode))
  ok <- !is.na(sample_barcode) & nchar(sample_barcode) >= 12
  out[ok] <- substr(sample_barcode[ok], 1, 12)
  out
}

extract_sample_code <- function(sample_barcode) {
  out <- rep(NA_character_, length(sample_barcode))
  ok <- !is.na(sample_barcode) & nchar(sample_barcode) >= 15
  out[ok] <- substr(sample_barcode[ok], 14, 15)
  out
}

extract_cancer <- function(filename) {
  x <- basename(filename)
  x <- sub("^HM450_TCGA-", "", x)
  sub("-.*$", "", x)
}

cat(">>> Reading cohort file...\n")

cohort <- fread(cohort_file, header = FALSE)

patients_keep <- unique(as.character(cohort[[1]]))
patients_keep <- patients_keep[grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}$", patients_keep)]

cat(">>> Patients in input noOV_noCHOL cohort:", length(patients_keep), "\n")

cat(">>> Listing methylation files...\n")

all_files <- list.files(
  meth_dir,
  full.names = TRUE
)

files <- all_files[
  endsWith(basename(all_files), "_annotated_methylation_filtered.bed.gz")
]

cat(">>> Total files matching methylation suffix:", length(files), "\n")

if (length(files) == 0) {
  stop("No methylation files found ending with _annotated_methylation_filtered.bed.gz")
}

dt <- data.table(file = files)

dt[, filename := basename(file)]
dt[, cancer := extract_cancer(filename)]
dt[, sample_barcode := extract_sample_barcode(filename)]
dt[, patient_id := extract_patient_id(sample_barcode)]
dt[, sample_code := extract_sample_code(sample_barcode)]

dt[, sample_type := fifelse(
  sample_code %in% sprintf("%02d", 1:9), "Tumor",
  fifelse(sample_code %in% sprintf("%02d", 10:19), "Healthy", NA_character_)
)]

dt[, sample_group := fifelse(
  sample_code %in% sprintf("%02d", 1:9), "0X",
  fifelse(sample_code %in% sprintf("%02d", 10:19), "1X", NA_character_)
)]

cat("\n>>> DEBUG: First few files parsed after correction:\n")
print(head(dt[, .(
  filename,
  cancer,
  sample_barcode,
  patient_id,
  sample_code,
  sample_type,
  sample_group
)], 10))

dt_keep <- dt[
  patient_id %in% patients_keep &
  sample_code %in% sprintf("%02d", 1:19) &
  !cancer %in% c("OV", "CHOL")
]

if (nrow(dt_keep) == 0) {
  cat("\n>>> DEBUG: Number of methylation patients overlapping cohort:\n")
  print(length(intersect(unique(dt$patient_id), patients_keep)))

  cat("\n>>> DEBUG: Sample codes found:\n")
  print(table(dt$sample_code, useNA = "ifany"))

  cat("\n>>> DEBUG: First few cohort patient IDs:\n")
  print(head(patients_keep, 20))

  cat("\n>>> DEBUG: First few methylation patient IDs:\n")
  print(head(unique(dt$patient_id), 20))

  stop("No methylation files retained. Check cohort file, file names, and filters.")
}

setorder(dt_keep, cancer, sample_type, sample_barcode)

cat("\n>>> Total methylation files parsed:", nrow(dt), "\n")
cat(">>> Kept noOV_noCHOL 0X/1X files:", nrow(dt_keep), "\n")
cat(">>> Kept patients:", length(unique(dt_keep$patient_id)), "\n")
cat(">>> Kept samples:", length(unique(dt_keep$sample_barcode)), "\n")
cat(">>> Kept cancers:", length(unique(dt_keep$cancer)), "\n")

cat("\n>>> Sample count by type:\n")
print(dt_keep[, .N, by = sample_type][order(sample_type)])

cat("\n>>> Sample count by cancer and type:\n")
print(dt_keep[, .N, by = .(cancer, sample_type)][order(cancer, sample_type)])

fwrite(dt_keep, out_file, sep = "\t")

cat("\n>>> Saved methylation file table:", out_file, "\n")
'

###############################################################################
# 2) Create pan-cancer methylation matrix
###############################################################################

Rscript -e '
suppressPackageStartupMessages({
  library(data.table)
})

files_table <- "./results/methylation/dimred_noOV_noCHOL/methylation_files_noOV_noCHOL_0X_1X.tsv"
out_rdata   <- "./results/methylation/dimred_noOV_noCHOL/TCGA_noOV_noCHOL_0X_1X_pan_cancer_matrix.RData"

cat(">>> Reading methylation file table...\n")
dt_files <- fread(files_table)

files <- dt_files$file
sample_names <- dt_files$sample_barcode

cat(">>> Number of methylation files:", length(files), "\n")

if (length(files) == 0) {
  stop("No methylation files found in file table.")
}

cat(">>> Reading first methylation file to define CpG order...\n")

dt0 <- fread(
  cmd = paste("zcat", shQuote(files[1])),
  header = FALSE,
  select = c(4, 5),
  showProgress = FALSE
)

cpgs <- as.character(dt0[[1]])

cat(">>> CpGs in first file:", length(cpgs), "\n")

pan_mat <- matrix(
  NA_real_,
  nrow = length(cpgs),
  ncol = length(files),
  dimnames = list(cpgs, sample_names)
)

pan_mat[, 1] <- as.numeric(dt0[[2]])

cat(">>> Filling methylation matrix...\n")

for (i in seq_along(files)) {
  if (i == 1) next

  if (i %% 100 == 0) {
    cat(">>> Processed", i, "of", length(files), "files\n")
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

cat(">>> Final methylation matrix dimensions:", paste(dim(pan_mat), collapse = " x "), "\n")

save(pan_mat, file = out_rdata)

cat(">>> Saved pan_mat:", out_rdata, "\n")
'

###############################################################################
# 3) Filter NA CpGs, median-impute, remove zero variance, select top 10k CpGs
###############################################################################

Rscript -e '
suppressPackageStartupMessages({
  library(matrixStats)
})

in_rdata  <- "./results/methylation/dimred_noOV_noCHOL/TCGA_noOV_noCHOL_0X_1X_pan_cancer_matrix.RData"
out_rdata <- "./results/methylation/dimred_noOV_noCHOL/TCGA_noOV_noCHOL_0X_1X_top10k_most_variable.RData"
out_txt   <- "./results/methylation/dimred_noOV_noCHOL/top10000_variable_CpGs_noOV_noCHOL_0X_1X.txt"
out_stats <- "./results/methylation/dimred_noOV_noCHOL/CpG_filtering_summary_noOV_noCHOL.tsv"

cat(">>> Checking input file...\n")

if (!file.exists(in_rdata)) {
  stop("Input file does not exist: ", in_rdata)
}

cat(">>> Loading pan-cancer methylation matrix...\n")
load(in_rdata)

if (!exists("pan_mat")) {
  stop("Object pan_mat was not found in: ", in_rdata)
}

cat(">>> pan_mat dimensions:", paste(dim(pan_mat), collapse = " x "), "\n")

n_cpg_initial <- nrow(pan_mat)
n_samples <- ncol(pan_mat)

###############################################################################
# 3A) Remove CpGs with >=50% missing values
###############################################################################

cat(">>> Filtering CpGs with >=50% missing values...\n")

na_fraction <- rowMeans(is.na(pan_mat))

keep_low_na <- na_fraction < 0.50

cat(">>> Initial CpGs:", n_cpg_initial, "\n")
cat(">>> CpGs removed because >=50% NA:", sum(!keep_low_na), "\n")

beta_clean <- pan_mat[keep_low_na, , drop = FALSE]

cat(">>> CpGs retained after NA filter:", nrow(beta_clean), "\n")

if (nrow(beta_clean) == 0) {
  stop("No CpGs retained after NA filtering.")
}

###############################################################################
# 3B) Median-impute remaining missing values
###############################################################################

cat(">>> Imputing remaining NA values using CpG-wise median...\n")

row_medians <- rowMedians(beta_clean, na.rm = TRUE)

names(row_medians) <- rownames(beta_clean)

na_idx <- which(is.na(beta_clean), arr.ind = TRUE)

if (nrow(na_idx) > 0) {
  beta_clean[na_idx] <- row_medians[na_idx[, 1]]
}

cat(">>> Number of imputed values:", nrow(na_idx), "\n")

if (anyNA(beta_clean)) {
  stop("NA values remain after median imputation.")
}

###############################################################################
# 3C) Remove zero-variance CpGs
###############################################################################

cat(">>> Removing zero-variance CpGs...\n")

vars <- rowVars(beta_clean)

names(vars) <- rownames(beta_clean)

valid_var <- is.finite(vars) & !is.na(vars) & vars > 0

cat(">>> CpGs removed because zero/non-finite variance:", sum(!valid_var), "\n")

beta_clean <- beta_clean[valid_var, , drop = FALSE]
vars <- vars[valid_var]

cat(">>> CpGs retained after variance filter:", nrow(beta_clean), "\n")

if (nrow(beta_clean) == 0) {
  stop("No CpGs retained after zero-variance filtering.")
}

###############################################################################
# 3D) Select top 10,000 most variable CpGs
###############################################################################

cat(">>> Selecting top 10,000 most variable CpGs...\n")

n_top <- min(10000, length(vars))

top_names <- names(sort(vars, decreasing = TRUE))[seq_len(n_top)]

beta_var <- beta_clean[top_names, , drop = FALSE]

cat(">>> beta_var dimensions:", paste(dim(beta_var), collapse = " x "), "\n")

save(beta_var, file = out_rdata)

writeLines(rownames(beta_var), out_txt)

summary_dt <- data.frame(
  step = c(
    "initial_CpGs",
    "removed_CpGs_with_ge50_percent_NA",
    "retained_after_NA_filter",
    "imputed_values",
    "removed_zero_or_nonfinite_variance_CpGs",
    "retained_after_variance_filter",
    "top_variable_CpGs_saved",
    "samples"
  ),
  value = c(
    n_cpg_initial,
    sum(!keep_low_na),
    nrow(pan_mat[keep_low_na, , drop = FALSE]),
    nrow(na_idx),
    sum(!valid_var),
    nrow(beta_clean),
    nrow(beta_var),
    n_samples
  )
)

write.table(
  summary_dt,
  file = out_stats,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat(">>> Saved top CpG matrix:", out_rdata, "\n")
cat(">>> Saved CpG list:", out_txt, "\n")
cat(">>> Saved filtering summary:", out_stats, "\n")
'

###############################################################################
# 4) PCA on cleaned top 10,000 most variable CpGs
###############################################################################

Rscript -e '
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

in_rdata <- "./results/methylation/dimred_noOV_noCHOL/TCGA_noOV_noCHOL_0X_1X_top10k_most_variable.RData"

files_table <- "./results/methylation/dimred_noOV_noCHOL/methylation_files_noOV_noCHOL_0X_1X.tsv"

color_file <- "./results/multi_omics/cancer_color_order_with_defined_colours.tsv"

out_coords <- "./results/methylation/dimred_noOV_noCHOL/PCA_coordinates_noOV_noCHOL_0X_1X_top10k.tsv"
out_pdf    <- "./results/methylation/dimred_noOV_noCHOL/PCA_noOV_noCHOL_0X_1X_top10k.pdf"
out_rds    <- "./results/methylation/dimred_noOV_noCHOL/PCA_noOV_noCHOL_0X_1X_top10k_prcomp.rds"

cat(">>> Checking PCA input file...\n")

if (!file.exists(in_rdata)) {
  stop("Input file does not exist: ", in_rdata)
}

cat(">>> Loading beta_var matrix...\n")
load(in_rdata)

if (!exists("beta_var")) {
  stop("Object beta_var was not found in: ", in_rdata)
}

cat(">>> beta_var dimensions:", paste(dim(beta_var), collapse = " x "), "\n")

if (anyNA(beta_var)) {
  stop("beta_var contains NA values. Step 3 failed.")
}

cat(">>> Transposing matrix for PCA...\n")

mat_pca <- t(beta_var)

cat(">>> PCA input dimensions:", paste(dim(mat_pca), collapse = " x "), "\n")
cat(">>> Rows = samples, columns = CpGs\n")

cat(">>> Running PCA...\n")

pca <- prcomp(
  mat_pca,
  center = TRUE,
  scale. = TRUE
)

saveRDS(pca, out_rds)

cat(">>> Saved PCA object:", out_rds, "\n")

var_explained <- pca$sdev^2 / sum(pca$sdev^2)

cat(">>> PC1 variance:", round(var_explained[1] * 100, 2), "%\n")
cat(">>> PC2 variance:", round(var_explained[2] * 100, 2), "%\n")

pca_df <- data.table(
  sample_barcode = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  PC3 = pca$x[, 3],
  PC4 = pca$x[, 4]
)

###############################################################################
# Metadata
###############################################################################

cat(">>> Loading sample metadata...\n")

meta <- fread(files_table)

meta[, sample_code := sprintf("%02d", as.integer(sample_code))]

meta[, sample_group := fifelse(
  sample_code %in% sprintf("%02d", 1:9), "0X",
  fifelse(sample_code %in% sprintf("%02d", 10:19), "1X", NA_character_)
)]

meta[, sample_type_clean := fifelse(
  sample_code %in% sprintf("%02d", 1:9), "Tumor",
  fifelse(sample_code %in% sprintf("%02d", 10:19), "Healthy", NA_character_)
)]

cat(">>> Sample type table:\n")
print(table(meta$sample_code, meta$sample_type_clean, useNA = "ifany"))

meta_small <- unique(meta[, .(
  sample_barcode,
  patient_id,
  cancer,
  sample_code,
  sample_group,
  sample_type_clean
)])

pca_df <- merge(
  pca_df,
  meta_small,
  by = "sample_barcode",
  all.x = TRUE
)

cat(">>> PCA samples with missing cancer metadata:", sum(is.na(pca_df$cancer)), "\n")
cat(">>> PCA samples with missing sample type:", sum(is.na(pca_df$sample_type_clean)), "\n")

fwrite(pca_df, out_coords, sep = "\t")

cat(">>> Saved PCA coordinates:", out_coords, "\n")

###############################################################################
# Colors
###############################################################################

cat(">>> Loading predefined cancer color file...\n")

color_dt <- fread(color_file, header = FALSE)

if (ncol(color_dt) < 2) {
  stop("Color file must contain at least two columns: cancer and color.")
}

color_dt <- color_dt[, .(cancer = V1, color = V2)]

color_dt[, cancer := gsub("^TCGA-", "", cancer)]

cancer_colors <- setNames(color_dt$color, color_dt$cancer)

present_cancers <- sort(unique(na.omit(pca_df$cancer)))
missing_colors <- setdiff(present_cancers, names(cancer_colors))

if (length(missing_colors) > 0) {
  cat(">>> WARNING: Missing colors for cancers:\n")
  print(missing_colors)
}

cancer_colors <- cancer_colors[names(cancer_colors) %in% present_cancers]

cat(">>> Cancers in PCA:", length(present_cancers), "\n")
cat(">>> Colors available:", length(cancer_colors), "\n")

###############################################################################
# Plot PCA
###############################################################################

cat(">>> Plotting PCA...\n")

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
    title = "PCA of DNA methylation profiles in the TCGA noOV_noCHOL cohort",
    subtitle = "Top 10,000 most variable CpG sites after NA filtering and median imputation",
    x = paste0("PC1 (", round(var_explained[1] * 100, 1), "% variance)"),
    y = paste0("PC2 (", round(var_explained[2] * 100, 1), "% variance)"),
    color = "Cancer type",
    shape = "Sample type"
  )

pdf(out_pdf, width = 13, height = 8)
print(p)
dev.off()

cat(">>> Saved PCA plot:", out_pdf, "\n")
cat(">>> DONE PCA\n")
'

###############################################################################
# 5) t-SNE on cleaned top 10,000 most variable CpGs
###############################################################################

Rscript -e '
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(Rtsne)
})

in_rdata <- "./results/methylation/dimred_noOV_noCHOL/TCGA_noOV_noCHOL_0X_1X_top10k_most_variable.RData"

files_table <- "./results/methylation/dimred_noOV_noCHOL/methylation_files_noOV_noCHOL_0X_1X.tsv"

color_file <- "./results/multi_omics/cancer_color_order_with_defined_colours.tsv"

out_coords <- "./results/methylation/dimred_noOV_noCHOL/tSNE_coordinates_noOV_noCHOL_0X_1X_top10k.tsv"
out_pdf    <- "./results/methylation/dimred_noOV_noCHOL/tSNE_noOV_noCHOL_0X_1X_top10k.pdf"
out_rds    <- "./results/methylation/dimred_noOV_noCHOL/tSNE_noOV_noCHOL_0X_1X_top10k_rtsne.rds"

set.seed(123)

initial_dims <- 30
max_iter     <- 1000

cat(">>> Checking t-SNE input file...\n")

if (!file.exists(in_rdata)) {
  stop("Input file does not exist: ", in_rdata)
}

cat(">>> Loading beta_var matrix...\n")
load(in_rdata)

if (!exists("beta_var")) {
  stop("Object beta_var was not found in: ", in_rdata)
}

cat(">>> beta_var dimensions:", paste(dim(beta_var), collapse = " x "), "\n")

if (anyNA(beta_var)) {
  stop("beta_var contains NA values. Step 3 failed.")
}

cat(">>> Transposing matrix for t-SNE...\n")

mat_tsne <- t(beta_var)

cat(">>> t-SNE input dimensions:", paste(dim(mat_tsne), collapse = " x "), "\n")
cat(">>> Rows = samples, columns = CpGs\n")

###############################################################################
# Scaling
###############################################################################

cat(">>> Scaling CpGs...\n")

mat_tsne <- scale(mat_tsne)

finite_cols <- apply(mat_tsne, 2, function(x) all(is.finite(x)))

mat_tsne <- mat_tsne[, finite_cols, drop = FALSE]

cat(">>> CpGs retained after scaling:", ncol(mat_tsne), "\n")

if (ncol(mat_tsne) == 0) {
  stop("No CpGs retained after scaling.")
}

###############################################################################
# t-SNE parameters
###############################################################################

n_samples <- nrow(mat_tsne)

perplexity <- min(30, floor((n_samples - 1) / 3))

if (perplexity < 2) {
  stop("Not enough samples for t-SNE. Need more samples to use perplexity >= 2.")
}

actual_initial_dims <- min(initial_dims, ncol(mat_tsne), n_samples - 1)

cat(">>> Number of samples:", n_samples, "\n")
cat(">>> Using perplexity:", perplexity, "\n")
cat(">>> Using initial_dims:", actual_initial_dims, "\n")
cat(">>> Using max_iter:", max_iter, "\n")

###############################################################################
# Run t-SNE
###############################################################################

cat(">>> Running t-SNE...\n")

tsne <- Rtsne(
  mat_tsne,
  dims = 2,
  perplexity = perplexity,
  initial_dims = actual_initial_dims,
  max_iter = max_iter,
  pca = TRUE,
  check_duplicates = FALSE,
  verbose = TRUE
)

saveRDS(tsne, out_rds)

cat(">>> Saved t-SNE object:", out_rds, "\n")

tsne_df <- data.table(
  sample_barcode = rownames(mat_tsne),
  tSNE1 = tsne$Y[, 1],
  tSNE2 = tsne$Y[, 2]
)

###############################################################################
# Metadata
###############################################################################

cat(">>> Loading sample metadata...\n")

meta <- fread(files_table)

meta[, sample_code := sprintf("%02d", as.integer(sample_code))]

meta[, sample_group := fifelse(
  sample_code %in% sprintf("%02d", 1:9), "0X",
  fifelse(sample_code %in% sprintf("%02d", 10:19), "1X", NA_character_)
)]

meta[, sample_type_clean := fifelse(
  sample_code %in% sprintf("%02d", 1:9), "Tumor",
  fifelse(sample_code %in% sprintf("%02d", 10:19), "Healthy", NA_character_)
)]

cat(">>> Sample type table:\n")
print(table(meta$sample_code, meta$sample_type_clean, useNA = "ifany"))

meta_small <- unique(meta[, .(
  sample_barcode,
  patient_id,
  cancer,
  sample_code,
  sample_group,
  sample_type_clean
)])

tsne_df <- merge(
  tsne_df,
  meta_small,
  by = "sample_barcode",
  all.x = TRUE
)

cat(">>> t-SNE samples with missing cancer metadata:", sum(is.na(tsne_df$cancer)), "\n")
cat(">>> t-SNE samples with missing sample type:", sum(is.na(tsne_df$sample_type_clean)), "\n")

fwrite(tsne_df, out_coords, sep = "\t")

cat(">>> Saved t-SNE coordinates:", out_coords, "\n")

###############################################################################
# Colors
###############################################################################

cat(">>> Loading predefined cancer color file...\n")

color_dt <- fread(color_file, header = FALSE)

if (ncol(color_dt) < 2) {
  stop("Color file must contain at least two columns: cancer and color.")
}

color_dt <- color_dt[, .(cancer = V1, color = V2)]

color_dt[, cancer := gsub("^TCGA-", "", cancer)]

cancer_colors <- setNames(color_dt$color, color_dt$cancer)

present_cancers <- sort(unique(na.omit(tsne_df$cancer)))
missing_colors <- setdiff(present_cancers, names(cancer_colors))

if (length(missing_colors) > 0) {
  cat(">>> WARNING: Missing colors for cancers:\n")
  print(missing_colors)
}

cancer_colors <- cancer_colors[names(cancer_colors) %in% present_cancers]

cat(">>> Cancers in t-SNE:", length(present_cancers), "\n")
cat(">>> Colors available:", length(cancer_colors), "\n")

###############################################################################
# Plot t-SNE
###############################################################################

cat(">>> Plotting t-SNE...\n")

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
    title = "t-SNE of DNA methylation profiles in the TCGA noOV_noCHOL cohort",
    subtitle = paste0(
      "Top 10,000 most variable CpG sites after NA filtering and median imputation; perplexity = ",
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

cat(">>> Saved t-SNE plot:", out_pdf, "\n")
cat(">>> DONE t-SNE\n")
'


########################################################################################################################
# Faceted multiple views using the files from the code above 
########################################################################################################################
Rscript - <<'EOF'
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(gridExtra)
  library(grid)
})

coords_file <- "./results/methylation/dimred_noOV_noCHOL/tSNE_coordinates_noOV_noCHOL_0X_1X_top10k.tsv"
color_file  <- "./results/multi_omics/cancer_color_order_with_defined_colours.tsv"

out_pdf <- "./results/methylation/dimred_noOV_noCHOL/tSNE_noOV_noCHOL_0X_1X_top10k_4views_onepage.pdf"

cat(">>> Reading t-SNE coordinates...\n")
tsne_df <- fread(coords_file)

###############################################################################
# Check required columns
###############################################################################

required_cols <- c("tSNE1", "tSNE2", "cancer", "sample_type_clean")
missing_cols <- setdiff(required_cols, names(tsne_df))

if (length(missing_cols) > 0) {
  stop("Missing required columns in t-SNE coordinate file: ",
       paste(missing_cols, collapse = ", "))
}

###############################################################################
# Cancer colors
###############################################################################

cat(">>> Reading cancer color file...\n")
color_dt <- fread(color_file, header = FALSE)

if (ncol(color_dt) < 2) {
  stop("Color file must contain at least two columns: cancer and color.")
}

color_dt <- color_dt[, .(
  cancer = gsub("^TCGA-", "", V1),
  color  = V2
)]

cancer_colors <- setNames(color_dt$color, color_dt$cancer)

present_cancers <- sort(unique(na.omit(tsne_df$cancer)))
missing_colors <- setdiff(present_cancers, names(cancer_colors))

if (length(missing_colors) > 0) {
  cat(">>> WARNING: Missing colors for cancers:\n")
  print(missing_colors)
}

cancer_colors <- cancer_colors[names(cancer_colors) %in% present_cancers]

###############################################################################
# Organ group mapping
###############################################################################

organ_map <- list(
  Adrenal = c("ACC", "PCPG"),
  Bladder = c("BLCA"),
  Brain_CNS = c("GBM", "LGG"),
  Breast = c("BRCA"),
  Cervix = c("CESC"),
  Colorectal = c("COAD", "READ"),
  Esophagus = c("ESCA"),
  Head_Neck = c("HNSC"),
  Hematologic = c("LAML", "DLBC"),
  Kidney = c("KICH", "KIRC", "KIRP"),
  Liver = c("LIHC"),
  Lung = c("LUAD", "LUSC"),
  Mesothelium = c("MESO"),
  Pancreas = c("PAAD"),
  Prostate = c("PRAD"),
  Sarcoma = c("SARC"),
  Skin = c("SKCM"),
  Stomach = c("STAD"),
  Testis = c("TGCT"),
  Thyroid = c("THCA"),
  Thymus = c("THYM"),
  Uterus = c("UCEC", "UCS"),
  Eye = c("UVM")
)

tsne_df[, organ_group := NA_character_]
for (nm in names(organ_map)) {
  tsne_df[cancer %in% organ_map[[nm]], organ_group := nm]
}
tsne_df[is.na(organ_group), organ_group := "Other"]

###############################################################################
# Tissue / cancer class mapping
###############################################################################

tissue_map <- list(
  Adenocarcinoma = c(
    "BRCA", "COAD", "READ", "STAD", "PAAD", "PRAD", "LUAD",
    "UCEC", "UCS", "LIHC", "THCA", "KIRC", "KIRP", "KICH", "ACC"
  ),
  Squamous = c("HNSC", "LUSC"),
  Mixed = c("ESCA", "CESC", "BLCA"),
  NonEpithelial = c(
    "GBM", "LGG", "LAML", "DLBC", "SKCM", "UVM", "SARC",
    "TGCT", "PCPG", "MESO", "THYM"
  )
)

tsne_df[, tissue_class := NA_character_]
for (nm in names(tissue_map)) {
  tsne_df[cancer %in% tissue_map[[nm]], tissue_class := nm]
}
tsne_df[is.na(tissue_class), tissue_class := "Other"]

###############################################################################
# Colors
###############################################################################

sample_type_colors <- c(
  Healthy = "darkgreen",
  Tumor = "red3"
)

organ_colors <- c(
  Adrenal = "#E69F00",
  Bladder = "#56B4E9",
  Brain_CNS = "#CC79A7",
  Breast = "#F781BF",
  Cervix = "#999999",
  Colorectal = "#009E73",
  Esophagus = "#D55E00",
  Head_Neck = "#A6761D",
  Hematologic = "#7570B3",
  Kidney = "#1B9E77",
  Liver = "#E6AB02",
  Lung = "#66A61E",
  Mesothelium = "#A6CEE3",
  Pancreas = "#FC8D62",
  Prostate = "#377EB8",
  Sarcoma = "#B2DF8A",
  Skin = "#FB9A99",
  Stomach = "#8DA0CB",
  Testis = "#CAB2D6",
  Thyroid = "#FFD92F",
  Thymus = "#BC80BD",
  Uterus = "#80B1D3",
  Eye = "#B3DE69",
  Other = "black"
)

organ_colors <- organ_colors[names(organ_colors) %in% unique(tsne_df$organ_group)]

tissue_colors <- c(
  Adenocarcinoma = "#1B9E77",
  Squamous = "#D95F02",
  Mixed = "#7570B3",
  NonEpithelial = "#E7298A",
  Other = "black"
)

tissue_colors <- tissue_colors[names(tissue_colors) %in% unique(tsne_df$tissue_class)]

###############################################################################
# Common theme
###############################################################################

base_theme <- theme_bw(base_family = "Times") +
  theme(
    text = element_text(family = "Times"),
    plot.title = element_text(hjust = 0.5, size = 19, family = "Times"),
    plot.subtitle = element_text(hjust = 0.5, size = 13, family = "Times"),
    axis.text = element_text(size = 7),                       # keep axis numbers same
    axis.title = element_text(size = 13, family = "Times"),
    legend.title = element_text(size = 14, family = "Times"),
    legend.text = element_text(size = 13, family = "Times"),
    legend.key.size = unit(0.40, "cm"),
    legend.position = "right",
    plot.margin = margin(6, 6, 6, 6)
  )

###############################################################################
# Plot 1: by cancer type, keeping tumor/healthy shape
###############################################################################

p1 <- ggplot(
  tsne_df,
  aes(x = tSNE1, y = tSNE2, color = cancer, shape = sample_type_clean)
) +
  geom_point(size = 1.4, alpha = 0.85) +
  scale_color_manual(values = cancer_colors, na.translate = FALSE) +
  scale_shape_manual(values = c(Tumor = 16, Healthy = 17), na.translate = FALSE) +
  base_theme +
  labs(
    title = "Cancer type",
    x = "t-SNE 1",
    y = "t-SNE 2",
    color = "Cancer",
    shape = "Sample"
  )

###############################################################################
# Plot 2: Healthy green, Tumor red
###############################################################################

p2 <- ggplot(
  tsne_df,
  aes(x = tSNE1, y = tSNE2, color = sample_type_clean, shape = sample_type_clean)
) +
  geom_point(size = 1.4, alpha = 0.85) +
  scale_color_manual(values = sample_type_colors, na.translate = FALSE) +
  scale_shape_manual(values = c(Tumor = 16, Healthy = 17), na.translate = FALSE) +
  base_theme +
  labs(
    title = "Tumor versus healthy",
    x = "t-SNE 1",
    y = "t-SNE 2",
    color = "Sample",
    shape = "Sample"
  )

###############################################################################
# Plot 3: by affected organ group
###############################################################################

p3 <- ggplot(
  tsne_df,
  aes(x = tSNE1, y = tSNE2, color = organ_group, shape = sample_type_clean)
) +
  geom_point(size = 1.4, alpha = 0.85) +
  scale_color_manual(values = organ_colors, na.translate = FALSE) +
  scale_shape_manual(values = c(Tumor = 16, Healthy = 17), na.translate = FALSE) +
  base_theme +
  labs(
    title = "Affected organ group",
    x = "t-SNE 1",
    y = "t-SNE 2",
    color = "Organ",
    shape = "Sample"
  )

###############################################################################
# Plot 4: by tissue / cancer class
###############################################################################

p4 <- ggplot(
  tsne_df,
  aes(x = tSNE1, y = tSNE2, color = tissue_class, shape = sample_type_clean)
) +
  geom_point(size = 1.4, alpha = 0.85) +
  scale_color_manual(values = tissue_colors, na.translate = FALSE) +
  scale_shape_manual(values = c(Tumor = 16, Healthy = 17), na.translate = FALSE) +
  base_theme +
  labs(
    title = "Tissue / cancer class",
    x = "t-SNE 1",
    y = "t-SNE 2",
    color = "Class",
    shape = "Sample"
  )

###############################################################################
# Save all 4 plots on one PDF page
###############################################################################

pdf(out_pdf, width = 18, height = 12, family = "Times")

grid.arrange(
  p1, p2,
  p3, p4,
  ncol = 2,
  nrow = 2,
  top = textGrob(
    "t-SNE of DNA methylation profiles in the 2D cohort",
    gp = gpar(fontsize = 18, fontfamily = "Times", fontface = "bold")
  )
)

dev.off()

cat(">>> Saved 4-view one-page t-SNE PDF:", out_pdf, "\n")

###############################################################################
# Optional: save annotated coordinates
###############################################################################

out_annot <- sub("\\.pdf$", "_annotated_coordinates.tsv", out_pdf)
fwrite(tsne_df, out_annot, sep = "\t")

cat(">>> Saved annotated t-SNE coordinates:", out_annot, "\n")
EOF
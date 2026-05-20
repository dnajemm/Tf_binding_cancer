###############################################################################
# PCA + t-SNE OF DNA METHYLATION PROFILES
# Parameterized R pipeline
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(matrixStats)
  library(ggplot2)
  library(Rtsne)
  library(gridExtra)
  library(grid)
})

###############################################################################
# 1) Argument parsing
###############################################################################

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL, required = FALSE) {
  idx <- match(flag, args)

  if (is.na(idx)) {
    if (required) {
      stop("Missing required argument: ", flag)
    }
    return(default)
  }

  if (idx == length(args)) {
    stop("Missing value after argument: ", flag)
  }

  args[idx + 1]
}

cohort_file <- get_arg(
  "--cohort_file",
  required = TRUE
)

top_cpgs <- as.integer(
  get_arg("--top_cpgs", default = "10000")
)

seed <- as.integer(
  get_arg("--seed", default = "123")
)

initial_dims <- as.integer(
  get_arg("--initial_dims", default = "30")
)

max_iter <- as.integer(
  get_arg("--max_iter", default = "1000")
)

###############################################################################
# 2) Validate arguments
###############################################################################

if (!file.exists(cohort_file)) {
  stop("The cohort file does not exist: ", cohort_file)
}

if (is.na(top_cpgs) || top_cpgs <= 0) {
  stop("--top_cpgs must be a positive integer.")
}

if (is.na(seed)) {
  stop("--seed must be an integer.")
}

if (is.na(initial_dims) || initial_dims <= 0) {
  stop("--initial_dims must be a positive integer.")
}

if (is.na(max_iter) || max_iter <= 0) {
  stop("--max_iter must be a positive integer.")
}

###############################################################################
# 3) Print run information
###############################################################################

cat("============================================================\n")
cat("DNA methylation dimensionality reduction pipeline\n")
cat("============================================================\n")
cat("Cohort file:        ", cohort_file, "\n", sep = "")
cat("Top variable CpGs:  ", top_cpgs, "\n", sep = "")
cat("t-SNE seed:         ", seed, "\n", sep = "")
cat("Initial dimensions: ", initial_dims, "\n", sep = "")
cat("Max iterations:     ", max_iter, "\n", sep = "")
cat("============================================================\n")

###############################################################################
# 4) Paths
###############################################################################

meth_dir <- "./methylation/filtered_methylation"

color_file <- "./results/multi_omics/cancer_color_order_with_defined_colours.tsv"

base_out <- "./results/methylation/dimred_pipeline"

run_name <- paste0(
  tools::file_path_sans_ext(basename(cohort_file)),
  "_top",
  top_cpgs
)

run_out <- file.path(base_out, run_name)

dir.create(base_out, recursive = TRUE, showWarnings = FALSE)
dir.create(run_out, recursive = TRUE, showWarnings = FALSE)

cat("Output directory:   ", run_out, "\n", sep = "")
cat("============================================================\n")

###############################################################################
# 5) Output files
###############################################################################

files_table <- file.path(
  run_out,
  "methylation_files_0X_1X.tsv"
)

pan_matrix_file <- file.path(
  run_out,
  "pan_cancer_methylation_matrix.RData"
)

top_matrix_file <- file.path(
  run_out,
  paste0("top", top_cpgs, "_most_variable_CpGs_matrix.RData")
)

top_cpg_list_file <- file.path(
  run_out,
  paste0("top", top_cpgs, "_variable_CpGs.txt")
)

cpg_summary_file <- file.path(
  run_out,
  paste0("CpG_filtering_summary_top", top_cpgs, ".tsv")
)

pca_coords_file <- file.path(
  run_out,
  paste0("PCA_coordinates_top", top_cpgs, ".tsv")
)

pca_pdf_file <- file.path(
  run_out,
  paste0("PCA_top", top_cpgs, ".pdf")
)

pca_rds_file <- file.path(
  run_out,
  paste0("PCA_top", top_cpgs, "_prcomp.rds")
)

tsne_coords_file <- file.path(
  run_out,
  paste0("tSNE_coordinates_top", top_cpgs, ".tsv")
)

tsne_pdf_file <- file.path(
  run_out,
  paste0("tSNE_top", top_cpgs, ".pdf")
)

tsne_organ_pdf_file <- file.path(
  run_out,
  paste0("tSNE_top", top_cpgs, "_affected_organ_only.pdf")
)

tsne_rds_file <- file.path(
  run_out,
  paste0("tSNE_top", top_cpgs, "_rtsne.rds")
)

tsne_4view_pdf_file <- file.path(
  run_out,
  paste0("tSNE_top", top_cpgs, "_4views_onepage.pdf")
)

###############################################################################
# 6) Step 1: Build methylation file table
###############################################################################

build_methylation_file_table <- function() {

  cat("\n============================================================\n")
  cat("STEP 1: Build methylation file table\n")
  cat("============================================================\n")

  if (file.exists(files_table)) {
    cat(">>> Output already exists, skipping Step 1:\n")
    cat(">>> ", files_table, "\n", sep = "")
    return(invisible(NULL))
  }

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

  cat(">>> Reading cohort file:\n")
  cat(">>> ", cohort_file, "\n", sep = "")

  cohort <- fread(cohort_file, header = FALSE)

  patients_keep <- unique(as.character(cohort[[1]]))

  patients_keep <- patients_keep[
    grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}$", patients_keep)
  ]

  cat(">>> Patients in input cohort:", length(patients_keep), "\n")

  if (length(patients_keep) == 0) {
    stop("No valid TCGA patient IDs found in the first column of cohort file.")
  }

  cat(">>> Listing methylation files from:\n")
  cat(">>> ", meth_dir, "\n", sep = "")

  all_files <- list.files(
    meth_dir,
    full.names = TRUE
  )

  files <- all_files[
    endsWith(
      basename(all_files),
      "_annotated_methylation_filtered.bed.gz"
    )
  ]

  cat(">>> Total methylation files matching suffix:", length(files), "\n")

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
    sample_code %in% sprintf("%02d", 1:9),
    "Tumor",
    fifelse(
      sample_code %in% sprintf("%02d", 10:19),
      "Healthy",
      NA_character_
    )
  )]

  dt[, sample_group := fifelse(
    sample_code %in% sprintf("%02d", 1:9),
    "0X",
    fifelse(
      sample_code %in% sprintf("%02d", 10:19),
      "1X",
      NA_character_
    )
  )]

###############################################################################
# Extract technical replicate suffix from filename
# Example:
# HM450_TCGA-BLCA-TCGA-BL-A0C8-01A_1_annotated_methylation_filtered.bed.gz
# HM450_TCGA-BLCA-TCGA-BL-A0C8-01A_2_annotated_methylation_filtered.bed.gz
###############################################################################

dt[, technical_replicate := sub(
  ".*_([0-9]+)_annotated_methylation_filtered\\.bed\\.gz$",
  "\\1",
  filename
)]

dt[, technical_replicate := as.integer(technical_replicate)]

  cat("\n>>> DEBUG: First parsed methylation files:\n")
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

    cat("\n>>> DEBUG: Sample codes found in methylation files:\n")
    print(table(dt$sample_code, useNA = "ifany"))

    cat("\n>>> DEBUG: First few cohort patient IDs:\n")
    print(head(patients_keep, 20))

    cat("\n>>> DEBUG: First few methylation patient IDs:\n")
    print(head(unique(dt$patient_id), 20))

    stop("No methylation files retained. Check cohort file, filenames, and filters.")
  }

    ###############################################################################
    # Remove duplicated sample barcodes caused by technical replicate files
    # Prefer replicate _1 over _2 when both are available
    ###############################################################################

    dup_dt <- dt_keep[
    sample_barcode %in% dt_keep[duplicated(sample_barcode), sample_barcode]
    ][order(sample_barcode, technical_replicate, file)]

    dup_file <- file.path(run_out, "duplicated_methylation_sample_barcodes.tsv")

    if (nrow(dup_dt) > 0) {

    cat("\n>>> WARNING: Duplicated sample_barcode values found.\n")
    cat(">>> These are likely technical replicate files distinguished by _1/_2.\n")
    cat(">>> Duplicated rows:", nrow(dup_dt), "\n")
    cat(">>> Duplicated sample_barcode values:", length(unique(dup_dt$sample_barcode)), "\n")
    cat(">>> Saving duplicated entries to:\n")
    cat(">>> ", dup_file, "\n", sep = "")

    fwrite(dup_dt, dup_file, sep = "\t")

    cat(">>> Keeping one methylation file per sample_barcode, preferring technical replicate _1.\n")

    setorder(dt_keep, sample_barcode, technical_replicate, file)

    dt_keep <- dt_keep[
        ,
        .SD[1],
        by = sample_barcode
    ]

    cat(">>> Samples retained after duplicate removal:", nrow(dt_keep), "\n")
    }

    setorder(dt_keep, cancer, sample_type, sample_barcode)

    cat("\n>>> Total methylation files parsed:", nrow(dt), "\n")
    cat(">>> Kept 0X/1X methylation files after duplicate removal:", nrow(dt_keep), "\n")
  cat(">>> Kept patients:", length(unique(dt_keep$patient_id)), "\n")
  cat(">>> Kept samples:", length(unique(dt_keep$sample_barcode)), "\n")
  cat(">>> Kept cancers:", length(unique(dt_keep$cancer)), "\n")

  cat("\n>>> Sample count by type:\n")
  print(dt_keep[, .N, by = sample_type][order(sample_type)])

  cat("\n>>> Sample count by cancer and type:\n")
  print(dt_keep[, .N, by = .(cancer, sample_type)][order(cancer, sample_type)])

  fwrite(dt_keep, files_table, sep = "\t")

  cat("\n>>> Saved methylation file table:\n")
  cat(">>> ", files_table, "\n", sep = "")

  invisible(NULL)
}

###############################################################################
# 7) Step 2: Create pan-cancer methylation matrix
###############################################################################

create_pan_cancer_matrix <- function() {

  cat("\n============================================================\n")
  cat("STEP 2: Create pan-cancer methylation matrix\n")
  cat("============================================================\n")

  if (file.exists(pan_matrix_file)) {
    cat(">>> Output already exists, skipping Step 2:\n")
    cat(">>> ", pan_matrix_file, "\n", sep = "")
    return(invisible(NULL))
  }

  if (!file.exists(files_table)) {
    stop("Step 1 output does not exist: ", files_table)
  }

  cat(">>> Reading methylation file table:\n")
  cat(">>> ", files_table, "\n", sep = "")

  dt_files <- fread(files_table)

  required_cols <- c("file", "sample_barcode")

  missing_cols <- setdiff(required_cols, names(dt_files))

  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in methylation file table: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  files <- dt_files$file
  sample_names <- dt_files$sample_barcode

  cat(">>> Number of methylation files:", length(files), "\n")

  if (length(files) == 0) {
    stop("No methylation files found in file table.")
  }

  missing_files <- files[!file.exists(files)]

  if (length(missing_files) > 0) {
    cat(">>> Missing methylation files:\n")
    print(head(missing_files, 20))
    stop("Some methylation files listed in Step 1 output do not exist.")
  }

  duplicated_samples <- sample_names[duplicated(sample_names)]

  if (length(duplicated_samples) > 0) {
    cat(">>> Duplicated sample barcodes:\n")
    print(unique(duplicated_samples))
    stop("Duplicated sample_barcode values found. Matrix column names must be unique.")
  }

  cat(">>> Reading first methylation file to define CpG order...\n")
  cat(">>> First file: ", files[1], "\n", sep = "")

  dt0 <- fread(
    cmd = paste("zcat", shQuote(files[1])),
    header = FALSE,
    select = c(4, 5),
    showProgress = FALSE
  )

  if (ncol(dt0) < 2) {
    stop("First methylation file did not return two selected columns.")
  }

  cpgs <- as.character(dt0[[1]])

  if (anyDuplicated(cpgs)) {
    dup_cpgs <- unique(cpgs[duplicated(cpgs)])
    cat(">>> Example duplicated CpGs in first file:\n")
    print(head(dup_cpgs, 20))
    stop("Duplicated CpG/probe IDs found in first methylation file.")
  }

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

    if (i %% 100 == 0 || i == length(files)) {
      cat(">>> Processed ", i, " of ", length(files), " files\n", sep = "")
    }

    dt <- fread(
      cmd = paste("zcat", shQuote(files[i])),
      header = FALSE,
      select = c(4, 5),
      showProgress = FALSE
    )

    probe_ids <- as.character(dt[[1]])
    beta_vals <- as.numeric(dt[[2]])

    if (length(probe_ids) != length(beta_vals)) {
      stop("Probe and beta lengths differ in file: ", files[i])
    }

    vals <- setNames(beta_vals, probe_ids)

    pan_mat[, i] <- vals[cpgs]
  }

  cat(">>> Final methylation matrix dimensions:", paste(dim(pan_mat), collapse = " x "), "\n")

  save(pan_mat, file = pan_matrix_file)

  cat(">>> Saved pan-cancer methylation matrix:\n")
  cat(">>> ", pan_matrix_file, "\n", sep = "")

  invisible(NULL)
}

###############################################################################
# 8) Step 3: Filter CpGs, impute NA values, select top variable CpGs
###############################################################################

select_top_variable_cpgs <- function() {

  cat("\n============================================================\n")
  cat("STEP 3: Filter CpGs and select top variable CpGs\n")
  cat("============================================================\n")

  if (file.exists(top_matrix_file)) {
    cat(">>> Output already exists, skipping Step 3:\n")
    cat(">>> ", top_matrix_file, "\n", sep = "")
    return(invisible(NULL))
  }

  if (!file.exists(pan_matrix_file)) {
    stop("Step 2 output does not exist: ", pan_matrix_file)
  }

  cat(">>> Loading pan-cancer methylation matrix:\n")
  cat(">>> ", pan_matrix_file, "\n", sep = "")

  load(pan_matrix_file)

  if (!exists("pan_mat")) {
    stop("Object pan_mat was not found in: ", pan_matrix_file)
  }

  cat(">>> pan_mat dimensions:", paste(dim(pan_mat), collapse = " x "), "\n")

  n_cpg_initial <- nrow(pan_mat)
  n_samples <- ncol(pan_mat)

  #############################################################################
  # 3A) Remove CpGs with >=50% missing values
  #############################################################################

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

  #############################################################################
  # 3B) Median-impute remaining missing values
  #############################################################################

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

  #############################################################################
  # 3C) Remove zero/non-finite variance CpGs
  #############################################################################

  cat(">>> Removing zero-variance or non-finite variance CpGs...\n")

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

  #############################################################################
  # 3D) Select top N most variable CpGs
  #############################################################################

  cat(">>> Selecting top ", top_cpgs, " most variable CpGs...\n", sep = "")

  n_top <- min(top_cpgs, length(vars))

  if (n_top < top_cpgs) {
    warning(
      "Requested ", top_cpgs,
      " CpGs, but only ", n_top,
      " valid variable CpGs are available."
    )
  }

  top_names <- names(sort(vars, decreasing = TRUE))[seq_len(n_top)]

  beta_var <- beta_clean[top_names, , drop = FALSE]

  cat(">>> beta_var dimensions:", paste(dim(beta_var), collapse = " x "), "\n")

  #############################################################################
  # Save outputs
  #############################################################################

  save(beta_var, file = top_matrix_file)

  writeLines(rownames(beta_var), top_cpg_list_file)

  summary_dt <- data.table(
    step = c(
      "initial_CpGs",
      "removed_CpGs_with_ge50_percent_NA",
      "retained_after_NA_filter",
      "imputed_values",
      "removed_zero_or_nonfinite_variance_CpGs",
      "retained_after_variance_filter",
      "top_variable_CpGs_requested",
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
      top_cpgs,
      nrow(beta_var),
      n_samples
    )
  )

  fwrite(summary_dt, cpg_summary_file, sep = "\t")

  cat(">>> Saved top CpG matrix:\n")
  cat(">>> ", top_matrix_file, "\n", sep = "")

  cat(">>> Saved CpG list:\n")
  cat(">>> ", top_cpg_list_file, "\n", sep = "")

  cat(">>> Saved filtering summary:\n")
  cat(">>> ", cpg_summary_file, "\n", sep = "")

  invisible(NULL)
}

###############################################################################
# 9) Step 4: PCA on top variable CpGs
###############################################################################

run_pca <- function() {

  cat("\n============================================================\n")
  cat("STEP 4: PCA on top variable CpGs\n")
  cat("============================================================\n")

  if (file.exists(pca_coords_file) && file.exists(pca_pdf_file) && file.exists(pca_rds_file)) {
    cat(">>> PCA outputs already exist, skipping Step 4:\n")
    cat(">>> ", pca_coords_file, "\n", sep = "")
    cat(">>> ", pca_pdf_file, "\n", sep = "")
    cat(">>> ", pca_rds_file, "\n", sep = "")
    return(invisible(NULL))
  }

  if (!file.exists(top_matrix_file)) {
    stop("Step 3 output does not exist: ", top_matrix_file)
  }

  if (!file.exists(files_table)) {
    stop("Methylation file table does not exist: ", files_table)
  }

  cat(">>> Loading top variable CpG matrix:\n")
  cat(">>> ", top_matrix_file, "\n", sep = "")

  load(top_matrix_file)

  if (!exists("beta_var")) {
    stop("Object beta_var was not found in: ", top_matrix_file)
  }

  cat(">>> beta_var dimensions:", paste(dim(beta_var), collapse = " x "), "\n")

  cat(">>> Reading methylation sample annotation:\n")
  cat(">>> ", files_table, "\n", sep = "")

  sample_annot <- fread(files_table)

  required_cols <- c("sample_barcode", "patient_id", "cancer", "sample_type")

  missing_cols <- setdiff(required_cols, names(sample_annot))

  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in sample annotation: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  sample_annot <- unique(
    sample_annot[, .(
      sample_barcode,
      patient_id,
      cancer,
      sample_type,
      sample_group
    )]
  )

  #############################################################################
  # Match annotation to beta matrix columns
  #############################################################################

  sample_names <- colnames(beta_var)

  sample_annot <- sample_annot[
    match(sample_names, sample_barcode)
  ]

  if (any(is.na(sample_annot$sample_barcode))) {
    stop("Some beta_var columns were not found in sample annotation.")
  }

  #############################################################################
  # Read cancer color file
  #############################################################################

  cat(">>> Reading cancer color file:\n")
  cat(">>> ", color_file, "\n", sep = "")

  if (!file.exists(color_file)) {
    stop("Cancer color file does not exist: ", color_file)
  }

  color_df <- fread(color_file, header = FALSE)

  if (ncol(color_df) < 2) {
    stop("Cancer color file must contain at least two columns: cancer and color.")
  }

  color_df <- color_df[, .(
    cancer_raw = as.character(V1),
    color = as.character(V2)
  )]

  color_df[, cancer := sub("^TCGA-", "", trimws(cancer_raw))]
  color_df[, color := trimws(gsub("\r", "", color))]

  cancer_colors <- setNames(color_df$color, color_df$cancer)

  missing_colors <- setdiff(unique(sample_annot$cancer), names(cancer_colors))

  if (length(missing_colors) > 0) {
    cat(">>> WARNING: Missing colors for cancers:\n")
    print(missing_colors)
  }

  #############################################################################
  # Prepare PCA input: samples x CpGs
  #############################################################################

  cat(">>> Preparing PCA input matrix: samples x CpGs...\n")

  pca_input <- t(beta_var)

  cat(">>> PCA input dimensions:", paste(dim(pca_input), collapse = " x "), "\n")

  if (anyNA(pca_input)) {
    stop("NA values found in PCA input. Step 3 should have imputed all missing values.")
  }

  #############################################################################
  # Run PCA
  #############################################################################

  cat(">>> Running PCA...\n")

  pca_res <- prcomp(
    pca_input,
    center = TRUE,
    scale. = TRUE
  )

  saveRDS(pca_res, pca_rds_file)

  cat(">>> Saved PCA object:\n")
  cat(">>> ", pca_rds_file, "\n", sep = "")

  #############################################################################
  # Variance explained
  #############################################################################

  eig <- pca_res$sdev^2
  var_explained <- eig / sum(eig) * 100

  pc1_lab <- paste0("PC1 (", round(var_explained[1], 2), "%)")
  pc2_lab <- paste0("PC2 (", round(var_explained[2], 2), "%)")
  pc3_lab <- paste0("PC3 (", round(var_explained[3], 2), "%)")
  pc4_lab <- paste0("PC4 (", round(var_explained[4], 2), "%)")

  #############################################################################
  # Save PCA coordinates
  #############################################################################

  pca_coords <- data.table(
    sample_barcode = rownames(pca_res$x),
    PC1 = pca_res$x[, 1],
    PC2 = pca_res$x[, 2],
    PC3 = pca_res$x[, 3],
    PC4 = pca_res$x[, 4]
  )

  pca_coords <- merge(
    pca_coords,
    sample_annot,
    by = "sample_barcode",
    all.x = TRUE
  )

  fwrite(pca_coords, pca_coords_file, sep = "\t")

  cat(">>> Saved PCA coordinates:\n")
  cat(">>> ", pca_coords_file, "\n", sep = "")

  #############################################################################
  # PCA plots
  #############################################################################

  p1 <- ggplot(
    pca_coords,
    aes(
      x = PC1,
      y = PC2,
      color = cancer,
      shape = sample_type
    )
  ) +
    geom_point(size = 2.2, alpha = 0.85) +
    scale_color_manual(values = cancer_colors, na.value = "grey70") +
    scale_shape_manual(values = c("Healthy" = 17, "Tumor" = 16)) +
    labs(
      title = paste0("PCA of DNA methylation profiles — top ", top_cpgs, " variable CpGs"),
      subtitle = "PC1 vs PC2",
      x = pc1_lab,
      y = pc2_lab,
      color = "Cancer type",
      shape = "Sample type"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text = element_text(color = "black"),
      legend.position = "right"
    )

  p2 <- ggplot(
    pca_coords,
    aes(
      x = PC3,
      y = PC4,
      color = cancer,
      shape = sample_type
    )
  ) +
    geom_point(size = 2.2, alpha = 0.85) +
    scale_color_manual(values = cancer_colors, na.value = "grey70") +
    scale_shape_manual(values = c("Healthy" = 17, "Tumor" = 16)) +
    labs(
      title = paste0("PCA of DNA methylation profiles — top ", top_cpgs, " variable CpGs"),
      subtitle = "PC3 vs PC4",
      x = pc3_lab,
      y = pc4_lab,
      color = "Cancer type",
      shape = "Sample type"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text = element_text(color = "black"),
      legend.position = "right"
    )

  pdf(pca_pdf_file, width = 12, height = 8)
  print(p1)
  print(p2)
  dev.off()

  cat(">>> Saved PCA PDF:\n")
  cat(">>> ", pca_pdf_file, "\n", sep = "")

  invisible(NULL)
}

###############################################################################
# 10) Step 5: t-SNE on PCA-reduced methylation data
###############################################################################

run_tsne <- function() {

  cat("\n============================================================\n")
  cat("STEP 5: t-SNE on PCA-reduced methylation data\n")
  cat("============================================================\n")

  if (
    file.exists(tsne_coords_file) &&
    file.exists(tsne_pdf_file) &&
    file.exists(tsne_rds_file)
  ) {
    cat(">>> t-SNE outputs already exist, skipping Step 5:\n")
    cat(">>> ", tsne_coords_file, "\n", sep = "")
    cat(">>> ", tsne_pdf_file, "\n", sep = "")
    cat(">>> ", tsne_rds_file, "\n", sep = "")
    return(invisible(NULL))
  }

  if (!file.exists(pca_rds_file)) {
    stop("PCA RDS file does not exist. Run Step 4 first: ", pca_rds_file)
  }

  if (!file.exists(files_table)) {
    stop("Methylation file table does not exist: ", files_table)
  }

  #############################################################################
  # Load PCA object
  #############################################################################

  cat(">>> Loading PCA object:\n")
  cat(">>> ", pca_rds_file, "\n", sep = "")

  pca_res <- readRDS(pca_rds_file)

  if (is.null(pca_res$x)) {
    stop("PCA object does not contain PCA coordinates.")
  }

  n_samples <- nrow(pca_res$x)
  n_pcs_available <- ncol(pca_res$x)

  cat(">>> PCA coordinates dimensions:", paste(dim(pca_res$x), collapse = " x "), "\n")

  #############################################################################
  # Choose number of PCs for t-SNE
  #############################################################################

  n_pcs_use <- min(initial_dims, n_pcs_available, n_samples - 1)

  if (n_pcs_use < 2) {
    stop("Not enough PCA dimensions available for t-SNE.")
  }

  cat(">>> PCs used for t-SNE:", n_pcs_use, "\n")

  tsne_input <- pca_res$x[, seq_len(n_pcs_use), drop = FALSE]

  if (anyNA(tsne_input)) {
    stop("NA values found in t-SNE input.")
  }

  if (any(!is.finite(tsne_input))) {
    stop("Non-finite values found in t-SNE input.")
  }

  #############################################################################
  # Define safe perplexity
  #############################################################################

  perplexity_value <- min(30, floor((n_samples - 1) / 3))

  if (perplexity_value < 5) {
    perplexity_value <- max(2, floor((n_samples - 1) / 3))
  }

  if (perplexity_value >= n_samples / 3) {
    perplexity_value <- floor((n_samples - 1) / 3)
  }

  cat(">>> Number of samples:", n_samples, "\n")
  cat(">>> t-SNE perplexity:", perplexity_value, "\n")
  cat(">>> t-SNE seed:", seed, "\n")
  cat(">>> t-SNE max_iter:", max_iter, "\n")

  #############################################################################
  # Run t-SNE
  #############################################################################

  set.seed(seed)

  tsne_res <- Rtsne(
    tsne_input,
    dims = 2,
    perplexity = perplexity_value,
    max_iter = max_iter,
    pca = FALSE,
    check_duplicates = FALSE,
    verbose = TRUE
  )

  saveRDS(tsne_res, tsne_rds_file)

  cat(">>> Saved t-SNE object:\n")
  cat(">>> ", tsne_rds_file, "\n", sep = "")

  #############################################################################
  # Read annotation
  #############################################################################

  sample_annot <- fread(files_table)

  required_cols <- c("sample_barcode", "patient_id", "cancer", "sample_type")

  missing_cols <- setdiff(required_cols, names(sample_annot))

  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in sample annotation: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  sample_annot <- unique(
    sample_annot[, .(
      sample_barcode,
      patient_id,
      cancer,
      sample_type,
      sample_group
    )]
  )

  sample_names <- rownames(tsne_input)

  sample_annot <- sample_annot[
    match(sample_names, sample_barcode)
  ]

  if (any(is.na(sample_annot$sample_barcode))) {
    stop("Some t-SNE samples were not found in sample annotation.")
  }

  #############################################################################
  # Read cancer colors
  #############################################################################

  cat(">>> Reading cancer color file:\n")
  cat(">>> ", color_file, "\n", sep = "")

  if (!file.exists(color_file)) {
    stop("Cancer color file does not exist: ", color_file)
  }

  color_df <- fread(color_file, header = FALSE)

  if (ncol(color_df) < 2) {
    stop("Cancer color file must contain at least two columns: cancer and color.")
  }

  color_df <- color_df[, .(
    cancer_raw = as.character(V1),
    color = as.character(V2)
  )]

  color_df[, cancer := sub("^TCGA-", "", trimws(cancer_raw))]
  color_df[, color := trimws(gsub("\r", "", color))]

  cancer_colors <- setNames(color_df$color, color_df$cancer)

  missing_colors <- setdiff(unique(sample_annot$cancer), names(cancer_colors))

  if (length(missing_colors) > 0) {
    cat(">>> WARNING: Missing colors for cancers:\n")
    print(missing_colors)
  }

  #############################################################################
  # Save t-SNE coordinates
  #############################################################################

  tsne_coords <- data.table(
    sample_barcode = sample_names,
    tSNE1 = tsne_res$Y[, 1],
    tSNE2 = tsne_res$Y[, 2],
    PCs_used = n_pcs_use,
    perplexity = perplexity_value,
    seed = seed,
    max_iter = max_iter
  )

  tsne_coords <- cbind(
    tsne_coords,
    sample_annot[, .(
      patient_id,
      cancer,
      sample_type,
      sample_group
    )]
  )

  fwrite(tsne_coords, tsne_coords_file, sep = "\t")

  cat(">>> Saved t-SNE coordinates:\n")
  cat(">>> ", tsne_coords_file, "\n", sep = "")

  #############################################################################
  # Main t-SNE plot
  #############################################################################

  p_tsne <- ggplot(
  tsne_coords,
  aes(
    x = tSNE1,
    y = tSNE2
  )
  ) +
    geom_point(
      data = tsne_coords[sample_type == "Tumor"],
      aes(color = cancer),
      shape = 16,
      size = 2.2,
      alpha = 0.85
    ) +
    geom_point(
      data = tsne_coords[sample_type == "Healthy"],
      aes(fill = cancer),
      shape = 24,
      color = "black",
      size = 2.8,
      stroke = 0.7,
      alpha = 0.95
    ) +
    scale_color_manual(values = cancer_colors, na.value = "grey70") +
    scale_fill_manual(values = cancer_colors, na.value = "grey70") +
    guides(fill = "none") +
    labs(
      title = paste0("t-SNE of DNA methylation profiles — top ", top_cpgs, " variable CpGs"),
      subtitle = paste0(
        "Input: first ", n_pcs_use,
        " PCA components | perplexity = ", perplexity_value,
        " | seed = ", seed
      ),
      x = "t-SNE 1",
      y = "t-SNE 2",
      color = "Cancer type",
      fill = "Cancer type"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text = element_text(color = "black"),
      legend.position = "right"
    )

  pdf(tsne_pdf_file, width = 12, height = 8)
  print(p_tsne)
  dev.off()

  cat(">>> Saved t-SNE PDF:\n")
  cat(">>> ", tsne_pdf_file, "\n", sep = "")

  invisible(NULL)
}

###############################################################################
# 11) Step 6: Faceted 4-view t-SNE plot
###############################################################################

run_tsne_4views <- function() {

  cat("\n============================================================\n")
  cat("STEP 6: Create 4-view t-SNE summary plot\n")
  cat("============================================================\n")

  if (file.exists(tsne_4view_pdf_file)) {
    cat(">>> 4-view t-SNE PDF already exists, skipping Step 6:\n")
    cat(">>> ", tsne_4view_pdf_file, "\n", sep = "")
    return(invisible(NULL))
  }

  if (!file.exists(tsne_coords_file)) {
    stop("t-SNE coordinates file does not exist. Run Step 5 first: ", tsne_coords_file)
  }

  if (!file.exists(color_file)) {
    stop("Cancer color file does not exist: ", color_file)
  }

  cat(">>> Reading t-SNE coordinates:\n")
  cat(">>> ", tsne_coords_file, "\n", sep = "")

  tsne_df <- fread(tsne_coords_file)

  required_cols <- c("tSNE1", "tSNE2", "cancer", "sample_type")

  missing_cols <- setdiff(required_cols, names(tsne_df))

  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in t-SNE coordinate file: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  #############################################################################
  # Clean sample type labels
  #############################################################################

  tsne_df[, sample_type := fifelse(
    sample_type %in% c("Tumor", "tumor", "Primary Tumor"),
    "Tumor",
    fifelse(
      sample_type %in% c("Healthy", "healthy", "Normal", "Solid Tissue Normal"),
      "Healthy",
      sample_type
    )
  )]

  tsne_df[, sample_type := factor(sample_type, levels = c("Tumor", "Healthy"))]

  #############################################################################
  # Cancer colors
  #############################################################################

  cat(">>> Reading cancer color file:\n")
  cat(">>> ", color_file, "\n", sep = "")

  color_dt <- fread(color_file, header = FALSE)

  if (ncol(color_dt) < 2) {
    stop("Cancer color file must contain at least two columns: cancer and color.")
  }

  color_dt <- color_dt[, .(
    cancer = sub("^TCGA-", "", trimws(as.character(V1))),
    color  = trimws(gsub("\r", "", as.character(V2)))
  )]

  cancer_colors <- setNames(color_dt$color, color_dt$cancer)

  present_cancers <- sort(unique(na.omit(tsne_df$cancer)))

  missing_colors <- setdiff(present_cancers, names(cancer_colors))

  if (length(missing_colors) > 0) {
    cat(">>> WARNING: Missing colors for cancers:\n")
    print(missing_colors)
  }

  cancer_colors <- cancer_colors[names(cancer_colors) %in% present_cancers]

  #############################################################################
  # Organ group mapping
  #############################################################################

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

  #############################################################################
  # Tissue / cancer class mapping
  #############################################################################

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

  #############################################################################
  # Colors
  #############################################################################

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

  #############################################################################
  # Common theme
  #############################################################################

  base_theme <- theme_bw(base_family = "Times") +
    theme(
      text = element_text(family = "Times"),
      plot.title = element_text(hjust = 0.5, size = 19, family = "Times"),
      plot.subtitle = element_text(hjust = 0.5, size = 13, family = "Times"),
      axis.text = element_text(size = 7, color = "black"),
      axis.title = element_text(size = 13, family = "Times"),
      legend.title = element_text(size = 14, family = "Times"),
      legend.text = element_text(size = 11, family = "Times"),
      legend.key.size = unit(0.40, "cm"),
      legend.position = "right",
      plot.margin = margin(6, 6, 6, 6)
    )

  shape_values <- c(
    Tumor = 16,
    Healthy = 17
  )

  #############################################################################
  # Plot 1: Cancer type
  #############################################################################

  p1 <- ggplot(
    tsne_df,
    aes(
      x = tSNE1,
      y = tSNE2,
      color = cancer,
      shape = sample_type
    )
  ) +
    geom_point(size = 1.4, alpha = 0.85) +
    scale_color_manual(values = cancer_colors, na.value = "grey70") +
    scale_shape_manual(values = shape_values, na.translate = FALSE) +
    base_theme +
    labs(
      title = "Cancer type",
      x = "t-SNE 1",
      y = "t-SNE 2",
      color = "Cancer",
      shape = "Sample"
    )

  #############################################################################
  # Plot 2: Tumor versus healthy
  #############################################################################

  p2 <- ggplot(
    tsne_df,
    aes(
      x = tSNE1,
      y = tSNE2,
      color = sample_type,
      shape = sample_type
    )
  ) +
    geom_point(size = 1.4, alpha = 0.85) +
    scale_color_manual(values = sample_type_colors, na.translate = FALSE) +
    scale_shape_manual(values = shape_values, na.translate = FALSE) +
    base_theme +
    labs(
      title = "Tumor versus healthy",
      x = "t-SNE 1",
      y = "t-SNE 2",
      color = "Sample",
      shape = "Sample"
    )

  #############################################################################
  # Plot 3: Organ group
  #############################################################################

  p3 <- ggplot(
    tsne_df,
    aes(
      x = tSNE1,
      y = tSNE2
    )
  ) +
    geom_point(
      data = tsne_df[sample_type == "Tumor"],
      aes(color = organ_group),
      shape = 16,
      size = 1.4,
      alpha = 0.85
    ) +
    geom_point(
      data = tsne_df[sample_type == "Healthy"],
      aes(fill = organ_group),
      shape = 24,
      color = "black",
      size = 1.8,
      stroke = 0.6,
      alpha = 0.95
    ) +
    scale_color_manual(values = organ_colors, na.translate = FALSE) +
    scale_fill_manual(values = organ_colors, na.translate = FALSE) +
    guides(fill = "none") +
    base_theme +
    labs(
      title = "Affected organ group",
      x = "t-SNE 1",
      y = "t-SNE 2",
      color = "Organ"
    )

  #############################################################################
  # Save affected organ group plot alone
  #############################################################################

  pdf(tsne_organ_pdf_file, width = 12, height = 8, family = "Times")
  print(p3)
  dev.off()

  cat(">>> Saved affected organ t-SNE PDF:\n")
  cat(">>> ", tsne_organ_pdf_file, "\n", sep = "")

  #############################################################################
  # Plot 4: Tissue / cancer class
  #############################################################################

  p4 <- ggplot(
    tsne_df,
    aes(
      x = tSNE1,
      y = tSNE2,
      color = tissue_class,
      shape = sample_type
    )
  ) +
    geom_point(size = 1.4, alpha = 0.85) +
    scale_color_manual(values = tissue_colors, na.translate = FALSE) +
    scale_shape_manual(values = shape_values, na.translate = FALSE) +
    base_theme +
    labs(
      title = "Tissue / cancer class",
      x = "t-SNE 1",
      y = "t-SNE 2",
      color = "Class",
      shape = "Sample"
    )

  #############################################################################
  # Save one-page PDF
  #############################################################################

  pdf(tsne_4view_pdf_file, width = 18, height = 12, family = "Times")

  grid.arrange(
    p1, p2,
    p3, p4,
    ncol = 2,
    nrow = 2,
    top = textGrob(
      paste0(
        "t-SNE of DNA methylation profiles — top ",
        top_cpgs,
        " variable CpGs"
      ),
      gp = gpar(
        fontsize = 18,
        fontfamily = "Times",
        fontface = "bold"
      )
    )
  )

  dev.off()

  cat(">>> Saved 4-view t-SNE PDF:\n")
  cat(">>> ", tsne_4view_pdf_file, "\n", sep = "")

  #############################################################################
  # Save annotated coordinates
  #############################################################################

  out_annot <- sub("\\.pdf$", "_annotated_coordinates.tsv", tsne_4view_pdf_file)

  fwrite(tsne_df, out_annot, sep = "\t")

  cat(">>> Saved annotated t-SNE coordinates:\n")
  cat(">>> ", out_annot, "\n", sep = "")

  invisible(NULL)
}

###############################################################################
# Run pipeline
###############################################################################

build_methylation_file_table()
create_pan_cancer_matrix()
select_top_variable_cpgs()
run_pca()
run_tsne()
run_tsne_4views()

cat("\n============================================================\n")
cat("Pipeline steps completed successfully.\n")
cat("Completed: Step 1 methylation file table\n")
cat("Completed: Step 2 pan-cancer methylation matrix\n")
cat("Completed: Step 3 top variable CpG selection\n")
cat("Completed: Step 4 PCA\n")
cat("Completed: Step 5 t-SNE\n")
cat("Completed: Step 6 4-view t-SNE summary plot\n")
cat("Output directory: ", run_out, "\n", sep = "")
cat("============================================================\n")

# exemple to run
# Rscript ./scripts/tSNE_methylation_pipeline.R --cohort_file ./results/multi_omics/samples_2d_noOV_noCHOL.tsv --top_cpgs 10000 --seed 123 --initial_dims 30 --max_iter 1000
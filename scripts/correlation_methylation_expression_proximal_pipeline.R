#!/usr/bin/env Rscript

# This pipeline performs a TF-centered pan-cancer analysis linking DNA methylation
# at proximal motif-associated CpG sites to target-gene expression across the TCGA
# 3d+4d_noOV cohort on 11X and 01X. Given a transcription factor such as BANP or NRF1,
# it builds a proximal motif–probe map, annotates HM450 methylation samples,
# computes sample-level motif methylation using the maximum beta value across
# motif-linked probes, merges these methylation values with RNA-seq expression,
# and calculates pooled as well as cancer-specific Pearson correlations with
# FDR correction. It then generates matched-patient analyses restricted to paired
# healthy and tumor samples, identifies significantly anti-correlated gene–motif
# pairs, and produces static PDF reports. Interactive HTML reports are optional.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(plotly)
  library(htmlwidgets)
  library(parallel)
})

# ============================================================
# ARGUMENTS
# ============================================================
args <- commandArgs(trailingOnly = TRUE)
tf_name <- if (length(args) >= 1) toupper(args[1]) else "BANP"

tf_map <- list(
  BANP = "BANP_mm0to2_noCGmm",
  NRF1 = "NRF1_mm0to2_noCGmm"
)

if (!tf_name %in% names(tf_map)) {
  stop("Unsupported TF: ", tf_name, ". Supported: ", paste(names(tf_map), collapse = ", "))
}

tf_motif_label <- tf_map[[tf_name]]

# ============================================================
# GLOBAL PARAMETERS
# ============================================================
cohort_file      <- "./results/multi_omics/samples_3d+4d_noOV.tsv"
probe_gene_file  <- "./methylation/probe_gene_pairs_in_motifs_with_tss.tsv.gz"
meth_dir         <- "./methylation/filtered_methylation"
expr_file        <- "./expression/gene_expression_matrix_3d4d_noOV.tsv"
color_file       <- "./results/multi_omics/cancer_color_order_with_defined_colours.tsv"
tf_expr_file     <- "./expression/expression_of_NRF1_BANP/all_samples_NRF1_BANP_expression.tsv"

base_dir <- file.path(
  ".", "results", "methylation", "correlation_expression_methylation",
  "tumor_vs_healthy", tf_name
)

all_dir     <- file.path(base_dir, "all_samples")
matched_dir <- file.path(base_dir, "matched_patients")

dir.create(all_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(matched_dir, recursive = TRUE, showWarnings = FALSE)

mc_cores        <- 20L
batch_size_meth <- 20L
batch_size_plot <- 5L
plot_cores      <- 5L
min_n_cor       <- 3L
min_n_wilcox    <- 2L
n_top           <- 50L
anti_r_cutoff   <- -0.3
anti_fdr_cutoff <- 0.05
only_significant_cancers <- FALSE

# Build interactive plotly HTML outputs?
# FALSE = skip interactive HTML generation.
# TRUE  = generate interactive HTML reports.
make_interactive_html <- FALSE

# Reuse already generated intermediate files if they exist.
# Useful when rerunning only for plots/filtering.
reuse_existing_intermediate_files <- TRUE

# Extra Wilcoxon filtering for by-cancer anti-correlated pairs
wilcox_filter_enabled <- TRUE
wilcox_fdr_cutoff     <- 0.05

# "both"        = require both expression and methylation significant
# "either"      = require either expression or methylation significant
# "expression"  = require only expression significant
# "methylation" = require only methylation significant
wilcox_filter_mode <- "both"

# ============================================================
# HELPERS
# ============================================================
msg <- function(...) cat(..., "\n")
sep_line <- function() cat(strrep("=", 70), "\n")

fmt_p <- function(x) {
  out <- rep(NA_character_, length(x))
  out[is.na(x)] <- "NA"
  out[!is.na(x) & x < 2.2e-16] <- "<2.2e-16"
  idx <- !is.na(x) & x >= 2.2e-16
  out[idx] <- formatC(x[idx], format = "e", digits = 2)
  out
}

p_to_stars <- function(p) {
  out <- rep(NA_character_, length(p))
  out[is.na(p)] <- "NA"
  out[!is.na(p) & p < 0.001] <- "***"
  out[!is.na(p) & p >= 0.001 & p < 0.01] <- "**"
  out[!is.na(p) & p >= 0.01 & p < 0.05] <- "*"
  out[!is.na(p) & p >= 0.05] <- "ns"
  out
}

fmt_sec <- function(x) {
  x <- as.numeric(x)
  if (!is.finite(x)) return("NA")
  if (x < 60) return(paste0(round(x, 1), " sec"))
  paste0(round(x / 60, 1), " min")
}

safe_name <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

clean_sample_type <- function(x) {
  x <- trimws(as.character(x))
  x[tolower(x) %in% c("healthy", "normal", "solid tissue normal", "adjacent normal")] <- "Healthy"
  x[tolower(x) %in% c("tumor", "primary tumor", "primary tumour", "cancer")] <- "Tumor"
  x
}

extract_sample_barcode_from_meth_file <- function(x) {
  x <- basename(as.character(x))
  out <- sub("^HM450_TCGA-[^-]+-(TCGA-[^-]+-[^-]+-[0-9]{2}[A-Z]).*$", "\\1", x)
  out[!grepl("^TCGA-[^-]+-[^-]+-[0-9]{2}[A-Z]$", out)] <- NA_character_
  out
}

extract_patient_id <- function(x) {
  x <- as.character(x)
  out <- sub("^(TCGA-[^-]+-[^-]+).*$", "\\1", x)
  out[!grepl("^TCGA-[^-]+-[^-]+$", out)] <- NA_character_
  out
}

extract_sample_type_tcga <- function(x) {
  x <- as.character(x)
  code <- sub("^TCGA-[^-]+-[^-]+-([0-9]{2})[A-Z]$", "\\1", x)
  fifelse(code == "01", "Tumor",
          fifelse(code == "11", "Healthy", NA_character_))
}

extract_expr_sample_barcode <- function(x) {
  x <- as.character(x)
  out <- sub("^TCGA-[^-]+-(TCGA-[^-]+-[^-]+-[0-9]{2}[A-Z])(?:_.*)?$", "\\1", x)
  out[!grepl("^TCGA-[^-]+-[^-]+-[0-9]{2}[A-Z]$", out)] <- NA_character_
  out
}

extract_expr_cancer <- function(x) {
  x <- as.character(x)
  out <- sub("^TCGA-([^-]+)-TCGA-[^-]+-[^-]+-[0-9]{2}[A-Z](?:_.*)?$", "\\1", x)
  out[!grepl("^[A-Z0-9]+$", out)] <- NA_character_
  out
}

extract_tcga_sample <- function(x) {
  x <- as.character(x)
  out <- sub(".*(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}[A-Z]).*", "\\1", x)
  out[!grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}[A-Z]$", out)] <- NA_character_
  out
}

safe_pearson <- function(x, y, min_n = 3L) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  n <- length(x)

  if (n < min_n) return(list(n = n, r = NA_real_, p = NA_real_))
  if (length(unique(x)) < 2 || length(unique(y)) < 2) {
    return(list(n = n, r = NA_real_, p = NA_real_))
  }

  ct <- tryCatch(cor.test(x, y, method = "pearson"), error = function(e) NULL)
  if (is.null(ct)) return(list(n = n, r = NA_real_, p = NA_real_))

  list(n = n, r = unname(ct$estimate), p = ct$p.value)
}

safe_wilcox_hc <- function(values, groups, min_n = 2L) {
  keep <- is.finite(values) & !is.na(groups)
  values <- values[keep]
  groups <- clean_sample_type(groups[keep])
  keep2 <- groups %in% c("Healthy", "Tumor")
  values <- values[keep2]
  groups <- groups[keep2]

  n_h <- sum(groups == "Healthy")
  n_t <- sum(groups == "Tumor")

  if (n_h < min_n || n_t < min_n) {
    return(list(p = NA_real_, n_h = n_h, n_t = n_t, status = "insufficient_n"))
  }

  wt <- tryCatch(
    suppressWarnings(wilcox.test(values ~ groups, exact = FALSE)),
    error = function(e) NULL
  )
  if (is.null(wt)) return(list(p = NA_real_, n_h = n_h, n_t = n_t, status = "wilcox_error"))

  list(p = wt$p.value, n_h = n_h, n_t = n_t, status = "ok")
}

make_dual_panel_dt <- function(d) {
  d <- copy(d)

  y1_min <- suppressWarnings(min(d$log2_expr, na.rm = TRUE))
  y1_max <- suppressWarnings(max(d$log2_expr, na.rm = TRUE))
  if (!is.finite(y1_min) || !is.finite(y1_max) || y1_min == y1_max) {
    y1_min <- 0
    y1_max <- 1
  }

  y2_min <- 0
  y2_max <- 100
  scale_factor <- (y1_max - y1_min) / (y2_max - y2_min)
  if (!is.finite(scale_factor) || scale_factor == 0) scale_factor <- 1

  expr_dt <- d[, .(
    sample_type, cancer,
    measure = "Expression",
    value_raw = log2_expr,
    value_plot = log2_expr
  )]

  meth_dt <- d[, .(
    sample_type, cancer,
    measure = "Methylation",
    value_raw = meth_percent,
    value_plot = (meth_percent - y2_min) * scale_factor + y1_min
  )]

  out <- rbindlist(list(expr_dt, meth_dt), fill = TRUE)

  out[, x_num := fifelse(
    measure == "Expression" & sample_type == "Healthy", 1,
    fifelse(
      measure == "Expression" & sample_type == "Tumor", 2,
      fifelse(
        measure == "Methylation" & sample_type == "Healthy", 4, 5
      )
    )
  )]

  list(
    dt = out,
    y1_min = y1_min,
    y1_max = y1_max,
    y2_min = y2_min,
    y2_max = y2_max,
    scale_factor = scale_factor
  )
}

load_color_map <- function(color_file) {
  col_dt <- fread(color_file, header = FALSE)
  if (ncol(col_dt) < 2) stop("color_file must have at least 2 columns")
  setnames(col_dt, c("cancer", "color"))
  col_dt <- unique(col_dt[, .(
    cancer = sub("^TCGA-", "", as.character(cancer)),
    color  = as.character(color)
  )])
  col_dt <- col_dt[!is.na(cancer) & nzchar(cancer) & !is.na(color) & nzchar(color)]
  setNames(col_dt$color, col_dt$cancer)
}

# ============================================================
# STEP 1
# ============================================================
build_motif_probe_map <- function() {
  sep_line(); msg("[Step 1] Build proximal motif-probe map")

  out_file <- file.path(all_dir, paste0(tf_name, "_proximal_motif_probe_map_long.tsv.gz"))

  if (reuse_existing_intermediate_files && file.exists(out_file)) {
    msg("Reusing existing motif-probe map:", out_file)
    return(out_file)
  }

  dt <- fread(probe_gene_file)

  dt <- dt[is_proximal == TRUE]
  dt <- dt[TF == tf_motif_label]

  out <- unique(dt[, .(
    gene = as.character(gene),
    motif_id = as.character(motif_instance_id),
    probe = as.character(probe)
  )])

  out <- out[
    !is.na(gene) & nzchar(gene) &
      !is.na(motif_id) & nzchar(motif_id) &
      !is.na(probe) & nzchar(probe)
  ]

  setorder(out, gene, motif_id, probe)
  fwrite(out, out_file, sep = "\t")

  msg("Saved:", out_file)
  msg("Rows:", nrow(out))
  msg("Unique genes:", uniqueN(out$gene))
  msg("Unique motifs:", uniqueN(out$motif_id))
  msg("Unique probes:", uniqueN(out$probe))
  out_file
}

# ============================================================
# STEP 2
# ============================================================
build_methylation_sample_annotation <- function() {
  sep_line(); msg("[Step 2] Build methylation sample annotation")

  out_file <- file.path(all_dir, "methylation_sample_annotation_3d4d_noOV.tsv.gz")

  if (reuse_existing_intermediate_files && file.exists(out_file)) {
    msg("Reusing existing methylation sample annotation:", out_file)
    return(out_file)
  }

  cohort <- fread(cohort_file, header = FALSE)
  patients_keep <- unique(as.character(cohort[[1]]))
  patients_keep <- patients_keep[grepl("^TCGA-[^-]+-[^-]+$", patients_keep)]

  files <- list.files(meth_dir, pattern = "\\.bed\\.gz$", full.names = TRUE)

  dt <- data.table(file = files, fname = basename(files))
  dt[, sample_barcode := extract_sample_barcode_from_meth_file(fname)]
  dt[, patient_id := extract_patient_id(sample_barcode)]
  dt[, sample_type := extract_sample_type_tcga(sample_barcode)]
  dt[, cancer := sub("^HM450_TCGA-([^-]+)-.*$", "\\1", fname)]

  dt <- dt[
    !is.na(sample_barcode) &
      !is.na(patient_id) &
      !is.na(sample_type)
  ]

  dt <- dt[patient_id %in% patients_keep]
  dt <- unique(dt[, .(patient_id, sample_barcode, sample_type, cancer, file)])
  setorder(dt, sample_type, cancer, patient_id, sample_barcode)

  fwrite(dt, out_file, sep = "\t")
  msg("Saved:", out_file)
  msg("Samples kept:", nrow(dt))
  msg("Tumor:", dt[sample_type == "Tumor", .N])
  msg("Healthy:", dt[sample_type == "Healthy", .N])

  out_file
}

# ============================================================
# STEP 3
# ============================================================
compute_motif_methylation <- function(motif_probe_file, meth_annot_file) {
  sep_line(); msg("[Step 3] Compute motif methylation using MAX(beta)")

  out_file_main    <- file.path(all_dir, "motif_sample_methylation_MAX_3d4d_noOV.tsv.gz")
  out_file_probes  <- file.path(all_dir, "motif_sample_methylation_MAX_selected_probes_3d4d_noOV.tsv.gz")
  out_file_compact <- file.path(all_dir, "motif_sample_methylation_MAX_selected_probes_compact_3d4d_noOV.tsv.gz")
  out_file_status  <- file.path(all_dir, "motif_sample_methylation_MAX_status_3d4d_noOV.tsv.gz")

  if (reuse_existing_intermediate_files && file.exists(out_file_main)) {
    msg("Reusing existing motif methylation file:", out_file_main)
    return(out_file_main)
  }

  mp <- fread(motif_probe_file)
  sa <- fread(meth_annot_file)

  mp <- unique(mp[!is.na(gene) & !is.na(motif_id) & !is.na(probe)])
  probe_keep <- unique(mp$probe)

  sa <- sa[file.exists(file)]
  if (nrow(sa) == 0) stop("No existing methylation files found.")

  process_one_sample <- function(i) {
    f  <- sa$file[i]
    sb <- sa$sample_barcode[i]
    st <- sa$sample_type[i]
    ca <- sa$cancer[i]

    dt <- tryCatch(
      fread(f, header = FALSE, select = c(4, 5), showProgress = FALSE),
      error = function(e) NULL
    )

    if (is.null(dt) || nrow(dt) == 0) {
      return(list(
        main = NULL,
        probes = NULL,
        status = data.table(sample_barcode = sb, ok = FALSE, reason = "file_read_failed_or_empty")
      ))
    }

    setnames(dt, c("probe", "beta"))
    dt[, probe := as.character(probe)]
    dt[, beta := suppressWarnings(as.numeric(beta))]
    dt <- dt[probe %in% probe_keep & !is.na(beta)]

    if (nrow(dt) == 0) {
      return(list(
        main = NULL,
        probes = NULL,
        status = data.table(sample_barcode = sb, ok = FALSE, reason = "no_matching_probes")
      ))
    }

    merged <- merge(mp, dt, by = "probe", all = FALSE)

    if (nrow(merged) == 0) {
      return(list(
        main = NULL,
        probes = NULL,
        status = data.table(sample_barcode = sb, ok = FALSE, reason = "empty_after_merge")
      ))
    }

    main_i <- merged[, {
      max_beta <- max(beta, na.rm = TRUE)
      .(
        meth_beta = max_beta,
        n_probes_found = .N,
        n_max_probes = sum(beta == max_beta, na.rm = TRUE)
      )
    }, by = .(gene, motif_id)]

    main_i[, `:=`(sample_barcode = sb, sample_type = st, cancer = ca)]

    probes_i <- merged[, {
      max_beta <- max(beta, na.rm = TRUE)
      .SD[beta == max_beta]
    }, by = .(gene, motif_id)]

    probes_i[, `:=`(
      sample_barcode = sb,
      sample_type = st,
      cancer = ca,
      is_max_probe = TRUE
    )]

    probes_i <- probes_i[, .(
      gene, motif_id, sample_barcode, sample_type, cancer, probe, beta, is_max_probe
    )]

    list(
      main = main_i,
      probes = probes_i,
      status = data.table(sample_barcode = sb, ok = TRUE, reason = "ok")
    )
  }

  idx_all <- seq_len(nrow(sa))
  batch_id <- ceiling(idx_all / batch_size_meth)
  batch_split <- split(idx_all, batch_id)

  main_chunks   <- vector("list", length(batch_split))
  probe_chunks  <- vector("list", length(batch_split))
  status_chunks <- vector("list", length(batch_split))

  for (b in seq_along(batch_split)) {
    idx <- batch_split[[b]]
    msg("batch", b, "/", length(batch_split), "| samples:", length(idx))

    res_batch <- mclapply(
      idx,
      process_one_sample,
      mc.cores = min(mc_cores, length(idx))
    )

    main_chunks[[b]]   <- rbindlist(lapply(res_batch, `[[`, "main"), fill = TRUE)
    probe_chunks[[b]]  <- rbindlist(lapply(res_batch, `[[`, "probes"), fill = TRUE)
    status_chunks[[b]] <- rbindlist(lapply(res_batch, `[[`, "status"), fill = TRUE)
  }

  main_res  <- rbindlist(main_chunks, fill = TRUE)
  probe_res <- rbindlist(probe_chunks, fill = TRUE)
  status_dt <- rbindlist(status_chunks, fill = TRUE)

  if (nrow(main_res) == 0) stop("No motif-level methylation values were recovered.")

  probe_summary_compact <- probe_res[, .(
    selected_probes = paste(probe, collapse = ","),
    selected_probe_betas = paste(beta, collapse = ",")
  ), by = .(gene, motif_id, sample_barcode, sample_type, cancer)]

  setcolorder(main_res, c("gene", "motif_id", "sample_barcode", "sample_type", "cancer", "meth_beta", "n_probes_found", "n_max_probes"))
  setorder(main_res, gene, motif_id, sample_type, cancer, sample_barcode)

  fwrite(main_res, out_file_main, sep = "\t")
  fwrite(probe_res, out_file_probes, sep = "\t")
  fwrite(probe_summary_compact, out_file_compact, sep = "\t")
  fwrite(status_dt, out_file_status, sep = "\t")

  msg("Saved:", out_file_main)
  msg("Main rows:", nrow(main_res))
  msg("Unique samples:", uniqueN(main_res$sample_barcode))

  out_file_main
}

# ============================================================
# STEP 4
# ============================================================
merge_with_expression <- function(meth_file) {
  sep_line(); msg("[Step 4] Merge methylation with expression")

  out_file <- file.path(all_dir, "motif_sample_expression_methylation_all.tsv.gz")

  if (reuse_existing_intermediate_files && file.exists(out_file)) {
    msg("Reusing existing merged methylation-expression file:", out_file)
    return(out_file)
  }

  meth <- fread(meth_file)
  expr <- fread(expr_file)

  if (!"sample" %in% names(expr)) stop("Expression file must contain column 'sample'")

  genes_keep <- unique(meth$gene)
  expr_cols_keep <- intersect(names(expr), c("sample", genes_keep))
  if (length(expr_cols_keep) <= 1) stop("No gene columns from methylation table found in expression matrix.")

  expr <- expr[, ..expr_cols_keep]

  expr_long <- melt(
    expr,
    id.vars = "sample",
    variable.name = "gene",
    value.name = "expression",
    variable.factor = FALSE
  )

  expr_long[, `:=`(
    sample = as.character(sample),
    gene = as.character(gene),
    expression = suppressWarnings(as.numeric(expression))
  )]

  expr_long[, sample_barcode := extract_expr_sample_barcode(sample)]
  expr_long[, cancer := extract_expr_cancer(sample)]
  expr_long[, sample_type := extract_sample_type_tcga(sample_barcode)]

  expr_long <- expr_long[
    !is.na(sample_barcode) &
      !is.na(gene) &
      !is.na(expression)
  ]

  merged <- merge(
    meth,
    expr_long[, .(gene, sample_barcode, expression)],
    by = c("gene", "sample_barcode"),
    all = FALSE
  )

  merged <- merged[!is.na(meth_beta) & !is.na(expression)]

  setcolorder(merged, c(
    "gene", "motif_id", "sample_barcode", "sample_type", "cancer",
    "meth_beta", "expression", "n_probes_found", "n_max_probes"
  ))
  setorder(merged, gene, motif_id, sample_type, cancer, sample_barcode)

  fwrite(merged, out_file, sep = "\t")
  msg("Saved:", out_file)
  msg("Merged rows:", nrow(merged))

  out_file
}

# ============================================================
# STEP 5
# ============================================================
run_correlations <- function(in_file, out_dir_corr, tumor_only = FALSE) {
  sep_line(); msg("[Step 5] Run Pearson correlations")

  dir.create(out_dir_corr, recursive = TRUE, showWarnings = FALSE)

  out_file_all       <- file.path(out_dir_corr, "pearson_correlation_all_pairs_all_cancers.tsv.gz")
  out_file_by_cancer <- file.path(out_dir_corr, "pearson_correlation_all_pairs_by_cancer.tsv.gz")

  dt <- fread(in_file)
  dt[, `:=`(
    gene = as.character(gene),
    motif_id = as.character(motif_id),
    sample_barcode = as.character(sample_barcode),
    sample_type = clean_sample_type(sample_type),
    cancer = as.character(cancer),
    meth_beta = suppressWarnings(as.numeric(meth_beta)),
    expression = suppressWarnings(as.numeric(expression))
  )]

  dt <- dt[
    !is.na(gene) & !is.na(motif_id) & !is.na(sample_barcode) &
      !is.na(cancer) & is.finite(meth_beta) & is.finite(expression)
  ]

  if (tumor_only) dt <- dt[sample_type == "Tumor"]

  res_all <- dt[, {
    z <- safe_pearson(meth_beta, expression, min_n = min_n_cor)
    .(n_all = z$n, pearson_r_all = z$r, pearson_p_all = z$p)
  }, by = .(gene, motif_id)]

  res_all[, tested := !is.na(pearson_p_all) & !is.na(n_all) & n_all >= min_n_cor]
  res_all[, pearson_fdr_all := p.adjust(pearson_p_all, method = "BH")]
  setcolorder(res_all, c("gene", "motif_id", "n_all", "pearson_r_all", "pearson_p_all", "pearson_fdr_all", "tested"))

  res_by_cancer <- dt[, {
    z <- safe_pearson(meth_beta, expression, min_n = min_n_cor)
    .(n = z$n, pearson_r = z$r, pearson_p = z$p)
  }, by = .(gene, motif_id, cancer)]

  res_by_cancer[, tested := !is.na(pearson_p) & !is.na(n) & n >= min_n_cor]

  # FDR corrected within each cancer type only.
  res_by_cancer[, pearson_fdr := p.adjust(pearson_p, method = "BH"), by = cancer]
  setcolorder(res_by_cancer, c("gene", "motif_id", "cancer", "n", "pearson_r", "pearson_p", "pearson_fdr", "tested"))
  setorder(res_by_cancer, cancer, pearson_fdr, pearson_p, gene, motif_id)

  fwrite(res_all, out_file_all, sep = "\t")
  fwrite(res_by_cancer, out_file_by_cancer, sep = "\t")

  msg("Saved:", out_file_all)
  msg("Saved:", out_file_by_cancer)

  list(all = out_file_all, by_cancer = out_file_by_cancer)
}

# ============================================================
# STEP 4b
# ============================================================
build_matched_subset <- function(in_file) {
  sep_line(); msg("[Step 4b] Build matched-patient subset")

  out_file <- file.path(matched_dir, "motif_sample_expression_methylation_matched_patients.tsv.gz")

  if (reuse_existing_intermediate_files && file.exists(out_file)) {
    msg("Reusing existing matched-patient file:", out_file)
    return(out_file)
  }

  dt <- fread(in_file)

  dt[, `:=`(
    gene = as.character(gene),
    motif_id = as.character(motif_id),
    sample_barcode = as.character(sample_barcode),
    sample_type = clean_sample_type(sample_type),
    cancer = as.character(cancer),
    meth_beta = suppressWarnings(as.numeric(meth_beta)),
    expression = suppressWarnings(as.numeric(expression))
  )]

  dt <- dt[
    !is.na(gene) & nzchar(gene) &
      !is.na(motif_id) & nzchar(motif_id) &
      !is.na(sample_barcode) & nzchar(sample_barcode) &
      sample_type %in% c("Healthy", "Tumor") &
      !is.na(cancer) & nzchar(cancer) &
      is.finite(meth_beta) & is.finite(expression)
  ]

  dt[, patient_id := sub("^(TCGA-[^-]+-[^-]+)-.*$", "\\1", sample_barcode)]

  matched_keys <- dt[, .(
    has_healthy = any(sample_type == "Healthy"),
    has_tumor   = any(sample_type == "Tumor")
  ), by = .(gene, motif_id, patient_id)][
    has_healthy == TRUE & has_tumor == TRUE,
    .(gene, motif_id, patient_id)
  ]

  dt_matched <- merge(dt, matched_keys, by = c("gene", "motif_id", "patient_id"), all = FALSE)

  setcolorder(dt_matched, c(
    "gene", "motif_id", "patient_id", "sample_barcode", "sample_type",
    "cancer", "meth_beta", "expression", "n_probes_found", "n_max_probes"
  ))
  setorder(dt_matched, gene, motif_id, cancer, patient_id, sample_type)

  fwrite(dt_matched, out_file, sep = "\t")
  msg("Saved:", out_file)
  msg("Matched rows:", nrow(dt_matched))
  msg("Unique matched patients:", uniqueN(dt_matched$patient_id))

  out_file
}

# ============================================================
# GENERIC PLOTTING CORE
# ============================================================
prepare_plot_data <- function(data_file) {
  dt <- fread(data_file)
  dt[, `:=`(
    cancer = as.character(cancer),
    gene = as.character(gene),
    motif_id = as.character(motif_id),
    sample_type = clean_sample_type(sample_type),
    meth_beta = suppressWarnings(as.numeric(meth_beta)),
    expression = suppressWarnings(as.numeric(expression))
  )]

  if ("patient_id" %in% names(dt)) dt[, patient_id := as.character(patient_id)]

  dt <- dt[
    !is.na(cancer) & nzchar(cancer) &
      !is.na(gene) & nzchar(gene) &
      !is.na(motif_id) & nzchar(motif_id) &
      sample_type %in% c("Healthy", "Tumor") &
      is.finite(meth_beta) & is.finite(expression)
  ]

  dt[, cancer := sub("^TCGA-", "", cancer)]
  dt[, log2_expr := log2(expression + 1)]

  if (max(dt$meth_beta, na.rm = TRUE) <= 1.5) {
    dt[, meth_percent := meth_beta * 100]
  } else {
    dt[, meth_percent := meth_beta]
  }

  dt
}

# ============================================================
# WILCOXON FILTER STATS
# ============================================================
compute_wilcox_filter_stats <- function(data_file) {
  msg("Computing Wilcoxon Healthy-vs-Tumor filter statistics")

  dt <- prepare_plot_data(data_file)

  wilcox_stats <- dt[, {
    expr_test <- safe_wilcox_hc(log2_expr, sample_type, min_n = min_n_wilcox)
    meth_test <- safe_wilcox_hc(meth_percent, sample_type, min_n = min_n_wilcox)

    .(
      n_healthy = sum(sample_type == "Healthy", na.rm = TRUE),
      n_tumor   = sum(sample_type == "Tumor", na.rm = TRUE),

      expr_wilcox_p = expr_test$p,
      expr_wilcox_status = expr_test$status,

      meth_wilcox_p = meth_test$p,
      meth_wilcox_status = meth_test$status
    )
  }, by = .(gene, motif_id, cancer)]

  wilcox_stats[, has_both_groups := n_healthy >= min_n_wilcox & n_tumor >= min_n_wilcox]

  # FDR within each cancer type.
  wilcox_stats[, expr_wilcox_fdr := p.adjust(expr_wilcox_p, method = "BH"), by = cancer]
  wilcox_stats[, meth_wilcox_fdr := p.adjust(meth_wilcox_p, method = "BH"), by = cancer]

  if (wilcox_filter_mode == "both") {
    wilcox_stats[, pass_wilcox_filter := fifelse(
      has_both_groups == FALSE,
      TRUE,
      !is.na(expr_wilcox_fdr) & expr_wilcox_fdr < wilcox_fdr_cutoff &
        !is.na(meth_wilcox_fdr) & meth_wilcox_fdr < wilcox_fdr_cutoff
    )]
  } else if (wilcox_filter_mode == "either") {
    wilcox_stats[, pass_wilcox_filter := fifelse(
      has_both_groups == FALSE,
      TRUE,
      (!is.na(expr_wilcox_fdr) & expr_wilcox_fdr < wilcox_fdr_cutoff) |
        (!is.na(meth_wilcox_fdr) & meth_wilcox_fdr < wilcox_fdr_cutoff)
    )]
  } else if (wilcox_filter_mode == "expression") {
    wilcox_stats[, pass_wilcox_filter := fifelse(
      has_both_groups == FALSE,
      TRUE,
      !is.na(expr_wilcox_fdr) & expr_wilcox_fdr < wilcox_fdr_cutoff
    )]
  } else if (wilcox_filter_mode == "methylation") {
    wilcox_stats[, pass_wilcox_filter := fifelse(
      has_both_groups == FALSE,
      TRUE,
      !is.na(meth_wilcox_fdr) & meth_wilcox_fdr < wilcox_fdr_cutoff
    )]
  } else {
    stop("Unsupported wilcox_filter_mode: ", wilcox_filter_mode)
  }

  wilcox_stats[, wilcox_filter_rule := fifelse(
    has_both_groups == TRUE,
    paste0("applied_", wilcox_filter_mode),
    "not_applied_missing_healthy_or_tumor"
  )]

  setorder(wilcox_stats, cancer, gene, motif_id)
  wilcox_stats[]
}

# ============================================================
# STEP 8 / 9
# ============================================================
extract_anti_correlated <- function(corr_files, out_dir_corr, data_file) {
  sep_line(); msg("[Step 8/9] Extract anti-correlated pairs")

  out_dir_anti <- file.path(out_dir_corr, "anti_correlated")
  dir.create(out_dir_anti, recursive = TRUE, showWarnings = FALSE)

  pooled_file <- file.path(out_dir_anti, "significant_anti_correlated_pairs_all_cancers.tsv.gz")
  cancer_file <- file.path(out_dir_anti, "significant_anti_correlated_pairs_by_cancer.tsv.gz")
  wilcox_file <- file.path(out_dir_anti, "wilcox_healthy_vs_tumor_by_cancer.tsv.gz")

  dt_all <- fread(corr_files$all)
  anti_all <- dt_all[
    !is.na(pearson_r_all) &
      !is.na(pearson_fdr_all) &
      pearson_r_all < anti_r_cutoff &
      pearson_fdr_all < anti_fdr_cutoff
  ][order(pearson_fdr_all, pearson_r_all)]

  fwrite(anti_all, pooled_file, sep = "\t")

  dt_bc <- fread(corr_files$by_cancer)
  dt_bc[, cancer := sub("^TCGA-", "", as.character(cancer))]

  anti_bc <- dt_bc[
    !is.na(pearson_r) &
      !is.na(pearson_fdr) &
      pearson_r < anti_r_cutoff &
      pearson_fdr < anti_fdr_cutoff
  ][order(cancer, pearson_fdr, pearson_r)]

  if (wilcox_filter_enabled && nrow(anti_bc) > 0) {
    wilcox_stats <- compute_wilcox_filter_stats(data_file)
    fwrite(wilcox_stats, wilcox_file, sep = "\t")

    anti_bc <- merge(
      anti_bc,
      wilcox_stats[, .(
        gene,
        motif_id,
        cancer,
        n_healthy,
        n_tumor,
        has_both_groups,
        expr_wilcox_p,
        expr_wilcox_fdr,
        expr_wilcox_status,
        meth_wilcox_p,
        meth_wilcox_fdr,
        meth_wilcox_status,
        wilcox_filter_rule,
        pass_wilcox_filter
      )],
      by = c("gene", "motif_id", "cancer"),
      all.x = TRUE
    )

    anti_bc[is.na(has_both_groups), has_both_groups := FALSE]
    anti_bc[is.na(pass_wilcox_filter), pass_wilcox_filter := FALSE]
    anti_bc[is.na(wilcox_filter_rule), wilcox_filter_rule := "missing_wilcox_statistics"]

    anti_bc_before <- nrow(anti_bc)
    anti_bc <- anti_bc[pass_wilcox_filter == TRUE]

    msg("Wilcoxon filter mode:", wilcox_filter_mode)
    msg("Wilcoxon FDR cutoff:", wilcox_fdr_cutoff)
    msg("By-cancer anti-correlated rows before Wilcoxon filter:", anti_bc_before)
    msg("By-cancer anti-correlated rows after Wilcoxon filter:", nrow(anti_bc))
    msg("Saved Wilcoxon table:", wilcox_file)
  } else {
    anti_bc[, `:=`(
      n_healthy = NA_integer_,
      n_tumor = NA_integer_,
      has_both_groups = NA,
      expr_wilcox_p = NA_real_,
      expr_wilcox_fdr = NA_real_,
      expr_wilcox_status = NA_character_,
      meth_wilcox_p = NA_real_,
      meth_wilcox_fdr = NA_real_,
      meth_wilcox_status = NA_character_,
      wilcox_filter_rule = "not_applied",
      pass_wilcox_filter = TRUE
    )]
  }

  setorder(anti_bc, cancer, pearson_fdr, pearson_r)

  fwrite(anti_bc, cancer_file, sep = "\t")

  msg("Saved:", pooled_file)
  msg("Saved:", cancer_file)
  msg("Anti-correlated pooled pairs:", nrow(anti_all))
  msg("Anti-correlated by-cancer pairs after Wilcoxon rule:", nrow(anti_bc))

  list(all = pooled_file, by_cancer = cancer_file)
}

# ============================================================
# PLOT ONE PAIR REPORT
# ============================================================
plot_pair_report <- function(plot_dt, pair_stats, sel_row, out_pdf, out_info, color_map_global, matched = FALSE) {
  pdf_open <- FALSE
  tryCatch({
    target_gene  <- as.character(sel_row$gene[1])
    target_motif <- as.character(sel_row$motif_id[1])

    cancers_present <- sort(unique(as.character(plot_dt$cancer)))
    col_sub <- data.table(cancer = cancers_present)
    col_sub[, color := color_map_global[cancer]]
    if (anyNA(col_sub$color)) col_sub[is.na(color), color := "grey50"]
    color_map <- setNames(col_sub$color, col_sub$cancer)

    plot_dt[, cancer := factor(as.character(cancer), levels = names(color_map))]
    plot_dt[, sample_type := factor(as.character(sample_type), levels = c("Healthy", "Tumor"))]

    count_dt <- plot_dt[, .(
      T = sum(sample_type == "Tumor", na.rm = TRUE),
      H = sum(sample_type == "Healthy", na.rm = TRUE)
    ), by = cancer]
    count_dt[, legend_lab := paste0(cancer, " (T=", T, ", H=", H, ")")]
    legend_map <- setNames(count_dt$legend_lab, count_dt$cancer)

    pair_stats <- copy(pair_stats)
    pair_stats[, signif := p_to_stars(pearson_fdr)]
    pair_stats[, cor_lab := paste0(
      "n = ", n,
      "\nr = ", ifelse(is.na(pearson_r), "NA", sprintf("%.2f", pearson_r)),
      "\np = ", fmt_p(pearson_p),
      "\nFDR = ", fmt_p(pearson_fdr),
      "\n", signif
    )]

    global_cor <- safe_pearson(plot_dt$meth_percent, plot_dt$log2_expr, min_n = min_n_cor)
    global_lab <- paste0(
      "gene = ", target_gene,
      "\nmotif = ", target_motif,
      "\nall-sample n = ", global_cor$n,
      "\nPearson r = ", ifelse(is.na(global_cor$r), "NA", sprintf("%.3f", global_cor$r)),
      "\np = ", fmt_p(global_cor$p)
    )

    per_cancer_dbg <- plot_dt[, .(
      n_total = .N,
      n_healthy = sum(sample_type == "Healthy", na.rm = TRUE),
      n_tumor = sum(sample_type == "Tumor", na.rm = TRUE),
      n_patients = if ("patient_id" %in% names(plot_dt)) uniqueN(patient_id) else NA_integer_,
      meth_unique_healthy = uniqueN(meth_percent[sample_type == "Healthy"]),
      meth_unique_tumor = uniqueN(meth_percent[sample_type == "Tumor"]),
      expr_unique_healthy = uniqueN(log2_expr[sample_type == "Healthy"]),
      expr_unique_tumor = uniqueN(log2_expr[sample_type == "Tumor"])
    ), by = cancer]

    fwrite(as.data.table(sel_row), out_info, sep = "\t")
    fwrite(per_cancer_dbg, sub("_summary.tsv$", "_per_cancer_counts.tsv", out_info), sep = "\t")

    pdf(out_pdf, width = 13, height = 10, onefile = TRUE)
    pdf_open <- TRUE

    plot_dt_p1 <- copy(plot_dt)
    plot_dt_p1[, sample_type := factor(as.character(sample_type), levels = c("Tumor", "Healthy"))]

    title1 <- if (matched) {
      paste0("Matched patients: ", target_gene, " | ", target_motif)
    } else {
      paste0("Top pair across all cancers: ", target_gene, " | ", target_motif)
    }

    p1 <- ggplot(plot_dt_p1, aes(x = meth_percent, y = log2_expr, color = cancer, shape = sample_type)) +
      geom_point(size = 2.2, alpha = 0.65) +
      geom_hline(yintercept = 1, color = "red", linewidth = 0.8, linetype = "dashed") +
      annotate("text", x = 95, y = 1, label = "expression = 1", color = "red", vjust = -0.5, hjust = 1, size = 4) +
      geom_smooth(
        data = plot_dt_p1,
        mapping = aes(x = meth_percent, y = log2_expr, group = 1),
        method = "lm", se = FALSE, inherit.aes = FALSE,
        color = "black", linewidth = 0.8
      ) +
      annotate("text", x = Inf, y = Inf, label = global_lab, hjust = 1.02, vjust = 1.05, size = 4) +
      scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100), expand = c(0.01, 0.01)) +
      scale_color_manual(values = color_map, breaks = names(color_map), labels = names(color_map), drop = FALSE, name = "Cancer type") +
      scale_shape_manual(values = c(Tumor = 16, Healthy = 17), breaks = c("Tumor", "Healthy"), labels = c("Tumor", "Healthy"), name = "Sample type") +
      guides(
        shape = guide_legend(order = 1, override.aes = list(size = 3, alpha = 1, color = "black")),
        color = guide_legend(order = 2, override.aes = list(shape = 16, size = 3, alpha = 1))
      ) +
      labs(
        title = title1,
        subtitle = "Shape shows sample type; color shows cancer type",
        x = "Motif methylation (%) [MAX probe beta]",
        y = "log2(expression + 1)"
      ) +
      theme_bw(base_size = 13) +
      theme(plot.title = element_text(face = "bold"), legend.position = "right")
    print(p1)

    for (ca in levels(plot_dt$cancer)) {
      dca <- plot_dt[cancer == ca]
      if (nrow(dca) == 0) next

      sc <- pair_stats[cancer == ca]
      sc_lab <- if (nrow(sc) == 0) {
        z <- safe_pearson(dca$meth_percent, dca$log2_expr, min_n = min_n_cor)
        paste0("n = ", z$n, "\nr = ", ifelse(is.na(z$r), "NA", sprintf("%.2f", z$r)), "\np = ", fmt_p(z$p), "\nFDR = NA\nNA")
      } else {
        sc$cor_lab[1]
      }

      dca_p <- copy(dca)
      dca_p[, draw_group := fifelse(sample_type == "Tumor", 1L, 2L)]
      setorder(dca_p, draw_group)

      subtitle_ca <- if (matched && "patient_id" %in% names(dca)) {
        paste0("Tumor n=", sum(dca$sample_type == "Tumor"),
               " ; Healthy n=", sum(dca$sample_type == "Healthy"),
               " ; Patients=", uniqueN(dca$patient_id))
      } else {
        paste0("Tumor n=", sum(dca$sample_type == "Tumor"),
               " ; Healthy n=", sum(dca$sample_type == "Healthy"))
      }

      # Red threshold line is intentionally added to every individual cancer scatter page.
      p_ca <- ggplot() +
        geom_point(data = dca_p[sample_type == "Tumor"], aes(x = meth_percent, y = log2_expr),
                   color = as.character(color_map[as.character(ca)]), shape = 16, size = 1.8, alpha = 0.40) +
        geom_point(data = dca_p[sample_type == "Healthy"], aes(x = meth_percent, y = log2_expr),
                   color = as.character(color_map[as.character(ca)]), shape = 17, size = 3.5, stroke = 1.0, alpha = 1) +
        geom_hline(yintercept = 1, color = "red", linewidth = 0.8, linetype = "dashed") +
        annotate("text", x = 95, y = 1, label = "expression = 1", color = "red", vjust = -0.5, hjust = 1, size = 4) +
        geom_smooth(data = dca, aes(x = meth_percent, y = log2_expr), method = "lm", se = FALSE,
                    color = "black", linewidth = 0.8) +
        annotate("text", x = min(dca$meth_percent, na.rm = TRUE), y = max(dca$log2_expr, na.rm = TRUE),
                 label = sc_lab, hjust = 0, vjust = 1, size = 4) +
        scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100), expand = c(0.01, 0.01)) +
        labs(
          title = paste0(if (matched) "Matched patients - " else "Scatter plot - ", ca, " | ", target_gene, " | ", target_motif),
          subtitle = subtitle_ca,
          x = "Motif methylation (%) [MAX probe beta]",
          y = "log2(expression + 1)"
        ) +
        theme_bw(base_size = 13) +
        theme(plot.title = element_text(face = "bold"))
      print(p_ca)
    }

    dual_all <- make_dual_panel_dt(plot_dt)
    d_all <- dual_all$dt
    sf_all <- dual_all$scale_factor
    y1_min_all <- dual_all$y1_min
    y2_min_all <- dual_all$y2_min

    expr_test_all <- safe_wilcox_hc(plot_dt$log2_expr, plot_dt$sample_type, min_n = min_n_wilcox)
    meth_test_all <- safe_wilcox_hc(plot_dt$meth_percent, plot_dt$sample_type, min_n = min_n_wilcox)
    expr_star_all <- p_to_stars(expr_test_all$p)
    meth_star_all <- p_to_stars(meth_test_all$p)

    expr_vals_all <- d_all[measure == "Expression"]$value_plot
    meth_vals_all <- d_all[measure == "Methylation"]$value_plot
    expr_range_all <- diff(range(expr_vals_all, na.rm = TRUE)); if (!is.finite(expr_range_all) || expr_range_all == 0) expr_range_all <- 1
    meth_range_all <- diff(range(meth_vals_all, na.rm = TRUE)); if (!is.finite(meth_range_all) || meth_range_all == 0) meth_range_all <- 1
    expr_y_all <- max(expr_vals_all, na.rm = TRUE) + 0.08 * expr_range_all
    meth_y_all <- max(meth_vals_all, na.rm = TRUE) + 0.08 * meth_range_all

    p5 <- ggplot() +
      geom_boxplot(data = d_all, aes(x = x_num, y = value_plot, group = x_num),
                   width = 0.55, outlier.shape = NA, color = "#0b3c5d", fill = NA) +
      geom_jitter(data = d_all, aes(x = x_num, y = value_plot, color = cancer),
                  width = 0.12, height = 0, size = 1.7, alpha = 0.55) +
      stat_summary(data = d_all, aes(x = x_num, y = value_plot, group = x_num),
                   fun = median, geom = "point", shape = 95, size = 10, color = "black") +
      annotate("segment", x = 1, xend = 2, y = expr_y_all, yend = expr_y_all, linewidth = 0.7) +
      annotate("text", x = 1.5, y = expr_y_all + 0.03 * expr_range_all, label = expr_star_all, size = 6) +
      annotate("segment", x = 4, xend = 5, y = meth_y_all, yend = meth_y_all, linewidth = 0.7) +
      annotate("text", x = 4.5, y = meth_y_all + 0.03 * meth_range_all, label = meth_star_all, size = 6) +
      annotate("segment", x = 3, xend = 3, y = min(d_all$value_plot, na.rm = TRUE), yend = max(d_all$value_plot, na.rm = TRUE), linewidth = 1.5) +
      annotate("text", x = 0.45, y = mean(range(d_all$value_plot, na.rm = TRUE)), label = "expression", angle = 90, size = 7) +
      annotate("text", x = 5.55, y = mean(range(d_all$value_plot, na.rm = TRUE)), label = "methylation", angle = 270, size = 7) +
      scale_x_continuous(breaks = c(1, 2, 4, 5), labels = c("H", "C", "H", "C"), limits = c(0.3, 5.7)) +
      scale_y_continuous(
        name = "log2(expression + 1)",
        sec.axis = sec_axis(trans = ~ (. - y1_min_all) / sf_all + y2_min_all,
                            name = "Motif methylation (%) [MAX probe beta]")
      ) +
      scale_color_manual(values = color_map, breaks = names(color_map), labels = legend_map[names(color_map)], drop = FALSE) +
      labs(
        title = paste0(if (matched) "Matched patients dual-panel: " else "Dual-panel summary across all cancers: ", target_gene, " | ", target_motif),
        subtitle = paste0(
          "Expression H vs C: ", expr_star_all, " (p=", fmt_p(expr_test_all$p), ", H=", expr_test_all$n_h, ", C=", expr_test_all$n_t, ") ; ",
          "Methylation H vs C: ", meth_star_all, " (p=", fmt_p(meth_test_all$p), ", H=", meth_test_all$n_h, ", C=", meth_test_all$n_t, ")"
        ),
        x = NULL, color = "Cancer"
      ) +
      theme_bw(base_size = 13) +
      theme(plot.title = element_text(face = "bold"), axis.title.x = element_blank(), legend.position = "right")
    print(p5)

    for (ca in levels(plot_dt$cancer)) {
      dca0 <- plot_dt[cancer == ca]
      if (nrow(dca0) == 0) next

      dual_ca <- make_dual_panel_dt(dca0)
      dca <- dual_ca$dt
      sf <- dual_ca$scale_factor
      y1_min <- dual_ca$y1_min
      y2_min <- dual_ca$y2_min

      expr_test <- safe_wilcox_hc(dca0$log2_expr, dca0$sample_type, min_n = min_n_wilcox)
      meth_test <- safe_wilcox_hc(dca0$meth_percent, dca0$sample_type, min_n = min_n_wilcox)
      expr_star <- p_to_stars(expr_test$p)
      meth_star <- p_to_stars(meth_test$p)

      expr_vals <- dca[measure == "Expression"]$value_plot
      meth_vals <- dca[measure == "Methylation"]$value_plot
      expr_range <- diff(range(expr_vals, na.rm = TRUE)); if (!is.finite(expr_range) || expr_range == 0) expr_range <- 1
      meth_range <- diff(range(meth_vals, na.rm = TRUE)); if (!is.finite(meth_range) || meth_range == 0) meth_range <- 1
      expr_y <- max(expr_vals, na.rm = TRUE) + 0.08 * expr_range
      meth_y <- max(meth_vals, na.rm = TRUE) + 0.08 * meth_range

      p6_ca <- ggplot() +
        geom_boxplot(data = dca, aes(x = x_num, y = value_plot, group = x_num),
                     width = 0.55, outlier.shape = NA, color = "#0b3c5d", fill = NA) +
        geom_jitter(data = dca, aes(x = x_num, y = value_plot),
                    color = as.character(color_map[as.character(ca)]), width = 0.12, height = 0, size = 1.8, alpha = 0.60) +
        stat_summary(data = dca, aes(x = x_num, y = value_plot, group = x_num),
                     fun = median, geom = "point", shape = 95, size = 10, color = "black") +
        annotate("segment", x = 1, xend = 2, y = expr_y, yend = expr_y, linewidth = 0.7) +
        annotate("text", x = 1.5, y = expr_y + 0.03 * expr_range, label = expr_star, size = 6) +
        annotate("segment", x = 4, xend = 5, y = meth_y, yend = meth_y, linewidth = 0.7) +
        annotate("text", x = 4.5, y = meth_y + 0.03 * meth_range, label = meth_star, size = 6) +
        annotate("segment", x = 3, xend = 3, y = min(dca$value_plot, na.rm = TRUE), yend = max(dca$value_plot, na.rm = TRUE), linewidth = 1.5) +
        annotate("text", x = 0.45, y = mean(range(dca$value_plot, na.rm = TRUE)), label = "expression", angle = 90, size = 7) +
        annotate("text", x = 5.55, y = mean(range(dca$value_plot, na.rm = TRUE)), label = "methylation", angle = 270, size = 7) +
        scale_x_continuous(breaks = c(1, 2, 4, 5), labels = c("H", "C", "H", "C"), limits = c(0.3, 5.7)) +
        scale_y_continuous(
          name = "log2(expression + 1)",
          sec.axis = sec_axis(trans = ~ (. - y1_min) / sf + y2_min,
                              name = "Motif methylation (%) [MAX probe beta]")
        ) +
        labs(
          title = paste0(if (matched) "Matched patients dual-panel - " else "Dual-panel summary - ", ca, " | ", target_gene, " | ", target_motif),
          subtitle = paste0(
            "Expression H vs C: ", expr_star, " (p=", fmt_p(expr_test$p), ", H=", expr_test$n_h, ", C=", expr_test$n_t, ") ; ",
            "Methylation H vs C: ", meth_star, " (p=", fmt_p(meth_test$p), ", H=", meth_test$n_h, ", C=", meth_test$n_t, ")"
          ),
          x = NULL, y = "log2(expression + 1)"
        ) +
        theme_bw(base_size = 13) +
        theme(plot.title = element_text(face = "bold"), axis.title.x = element_blank())
      print(p6_ca)
    }

    dev.off()
    TRUE
  }, error = function(e) {
    if (pdf_open) try(dev.off(), silent = TRUE)
    stop(e)
  })
}

# ============================================================
# STEP 6 / 7b
# ============================================================
plot_top_pairs <- function(data_file, stats_file, out_dir_plot, matched = FALSE) {
  sep_line(); msg(if (matched) "[Step 7b] Plot top 50 matched pairs" else "[Step 6] Plot top 50 pairs")

  dir.create(out_dir_plot, recursive = TRUE, showWarnings = FALSE)

  dt <- prepare_plot_data(data_file)
  st <- fread(stats_file)
  st[, `:=`(
    cancer = sub("^TCGA-", "", as.character(cancer)),
    gene = as.character(gene),
    motif_id = as.character(motif_id),
    n = as.integer(n),
    pearson_r = suppressWarnings(as.numeric(pearson_r)),
    pearson_p = suppressWarnings(as.numeric(pearson_p)),
    pearson_fdr = suppressWarnings(as.numeric(pearson_fdr))
  )]
  if (!"tested" %in% names(st)) st[, tested := !is.na(pearson_p) & !is.na(n) & n >= min_n_cor]

  top_pairs <- st[
    tested == TRUE & !is.na(pearson_fdr)
  ][order(gene, motif_id, pearson_fdr, -abs(pearson_r), -n)][
    , .SD[1], by = .(gene, motif_id)
  ][order(pearson_fdr, -abs(pearson_r), -n)]

  top_pairs <- top_pairs[, .(
    gene, motif_id,
    source_cancer = cancer,
    source_n = n,
    source_r = pearson_r,
    source_p = pearson_p,
    source_padj = pearson_fdr
  )]

  top_pairs <- top_pairs[seq_len(min(n_top, .N))]
  if (nrow(top_pairs) == 0) stop("No valid tested pairs found in stats_file.")
  top_pairs[, rank_idx := seq_len(.N)]

  top_index_file <- file.path(out_dir_plot, "top50_selected_pairs.tsv")
  out_index      <- file.path(out_dir_plot, "top50_generated_reports.tsv")
  fwrite(top_pairs, top_index_file, sep = "\t")

  color_map_global <- load_color_map(color_file)

  process_one_pair <- function(sel_row) {
    tryCatch({
      target_gene  <- as.character(sel_row$gene[1])
      target_motif <- as.character(sel_row$motif_id[1])
      rank_idx     <- as.integer(sel_row$rank_idx[1])

      plot_dt <- dt[gene == target_gene & motif_id == target_motif]
      if (nrow(plot_dt) == 0) {
        return(data.table(rank_idx = rank_idx, gene = target_gene, motif_id = target_motif,
                          status = "no_rows", error_message = NA_character_,
                          pdf_file = NA_character_, summary_file = NA_character_))
      }

      pair_stats <- st[gene == target_gene & motif_id == target_motif, .(
        cancer, gene, motif_id, n, pearson_r, pearson_p, pearson_fdr, tested
      )]

      missing_stat_cancers <- setdiff(as.character(unique(plot_dt$cancer)), as.character(pair_stats$cancer))
      if (length(missing_stat_cancers)) {
        extra <- plot_dt[cancer %in% missing_stat_cancers, {
          z <- safe_pearson(meth_percent, log2_expr, min_n = min_n_cor)
          .(n = z$n, pearson_r = z$r, pearson_p = z$p, pearson_fdr = NA_real_, tested = !is.na(z$p) & z$n >= min_n_cor)
        }, by = .(cancer, gene, motif_id)]
        pair_stats <- rbindlist(list(pair_stats, extra), fill = TRUE)
      }

      base_name <- paste0(sprintf("%02d", rank_idx), "_", safe_name(target_gene), "__", safe_name(target_motif))
      out_pdf  <- file.path(out_dir_plot, paste0(base_name, "_report.pdf"))
      out_info <- file.path(out_dir_plot, paste0(base_name, "_summary.tsv"))

      plot_pair_report(plot_dt, pair_stats, sel_row, out_pdf, out_info, color_map_global, matched = matched)

      data.table(rank_idx = rank_idx, gene = target_gene, motif_id = target_motif,
                 status = "ok", error_message = NA_character_,
                 pdf_file = out_pdf, summary_file = out_info)
    }, error = function(e) {
      data.table(rank_idx = as.integer(sel_row$rank_idx[1]),
                 gene = as.character(sel_row$gene[1]),
                 motif_id = as.character(sel_row$motif_id[1]),
                 status = "error", error_message = conditionMessage(e),
                 pdf_file = NA_character_, summary_file = NA_character_)
    })
  }

  idx_all <- seq_len(nrow(top_pairs))
  batch_id <- ceiling(idx_all / batch_size_plot)
  batch_split <- split(idx_all, batch_id)
  res_list <- vector("list", length(batch_split))
  t0 <- proc.time()[3]

  for (b in seq_along(batch_split)) {
    tb <- proc.time()[3]
    idx <- batch_split[[b]]
    batch_dt <- top_pairs[idx]
    msg("batch", b, "/", length(batch_split), "| pairs:", nrow(batch_dt))

    batch_rows <- lapply(seq_len(nrow(batch_dt)), function(i) batch_dt[i])

    batch_res <- mclapply(
      batch_rows, process_one_pair,
      mc.cores = min(plot_cores, length(batch_rows))
    )
    res_list[[b]] <- rbindlist(batch_res, fill = TRUE)

    elapsed_batch <- proc.time()[3] - tb
    elapsed_total <- proc.time()[3] - t0
    avg_per_batch <- elapsed_total / b
    eta <- (length(batch_split) - b) * avg_per_batch

    msg("done | batch time:", fmt_sec(elapsed_batch), "| elapsed:", fmt_sec(elapsed_total), "| ETA:", fmt_sec(eta))
  }

  res <- rbindlist(res_list, fill = TRUE)
  fwrite(res, out_index, sep = "\t")

  msg("Saved:", top_index_file)
  msg("Saved:", out_index)

  invisible(list(index = top_index_file, reports = out_index))
}

# ============================================================
# STEP 10
# ============================================================
plot_all_anti_correlated <- function(data_file, anti_files, out_dir_plot) {
  sep_line(); msg("[Step 10] Plot all anti-correlated pairs")

  dir.create(out_dir_plot, recursive = TRUE, showWarnings = FALSE)

  dt <- prepare_plot_data(data_file)
  anti_all <- fread(anti_files$all)
  anti_bc  <- fread(anti_files$by_cancer)
  anti_bc[, cancer := sub("^TCGA-", "", as.character(cancer))]

  pairs_from_all <- unique(anti_all[, .(gene, motif_id)])
  pairs_from_bc  <- unique(anti_bc[, .(gene, motif_id)])
  pair_union     <- unique(rbindlist(list(pairs_from_all, pairs_from_bc), fill = TRUE))

  all_sum <- copy(anti_all)[, .(
    pooled_n = n_all[1],
    pooled_r = pearson_r_all[1],
    pooled_p = pearson_p_all[1],
    pooled_fdr = pearson_fdr_all[1],
    in_pooled = TRUE
  ), by = .(gene, motif_id)]

  bc_sum <- anti_bc[, .(
    n_sig_cancers = uniqueN(cancer),
    sig_cancers = paste(sort(unique(cancer)), collapse = ","),
    best_cancer = cancer[order(pearson_fdr, pearson_r, -n)][1],
    best_n = n[order(pearson_fdr, pearson_r, -n)][1],
    best_r = pearson_r[order(pearson_fdr, pearson_r, -n)][1],
    best_p = pearson_p[order(pearson_fdr, pearson_r, -n)][1],
    best_fdr = pearson_fdr[order(pearson_fdr, pearson_r, -n)][1],
    in_by_cancer = TRUE
  ), by = .(gene, motif_id)]

  top_pairs <- merge(pair_union, all_sum, by = c("gene", "motif_id"), all.x = TRUE)
  top_pairs <- merge(top_pairs, bc_sum, by = c("gene", "motif_id"), all.x = TRUE)

  top_pairs[is.na(in_pooled), in_pooled := FALSE]
  top_pairs[is.na(in_by_cancer), in_by_cancer := FALSE]
  top_pairs[is.na(n_sig_cancers), n_sig_cancers := 0L]
  top_pairs[is.na(sig_cancers), sig_cancers := ""]
  top_pairs[, rank_fdr := fifelse(!is.na(pooled_fdr), pooled_fdr, best_fdr)]
  top_pairs[, rank_r   := fifelse(!is.na(pooled_r), pooled_r, best_r)]
  top_pairs <- top_pairs[order(rank_fdr, rank_r, -n_sig_cancers, gene, motif_id)]
  top_pairs[, rank_idx := seq_len(.N)]

  top_index_file <- file.path(out_dir_plot, paste0("selected_anti_correlated_pairs_all_", tf_name, ".tsv"))
  out_index      <- file.path(out_dir_plot, paste0("generated_reports_", tf_name, ".tsv"))
  fwrite(top_pairs, top_index_file, sep = "\t")

  color_map_global <- load_color_map(color_file)

  process_one_pair <- function(sel_row) {
    tryCatch({
      target_gene  <- as.character(sel_row$gene[1])
      target_motif <- as.character(sel_row$motif_id[1])
      rank_idx     <- as.integer(sel_row$rank_idx[1])

      plot_dt <- dt[gene == target_gene & motif_id == target_motif]
      if (nrow(plot_dt) == 0) {
        return(data.table(tf = tf_name, rank_idx = rank_idx, gene = target_gene, motif_id = target_motif,
                          status = "no_rows", error_message = NA_character_, pdf_file = NA_character_, summary_file = NA_character_))
      }

      sig_cancers_this_pair <- anti_bc[gene == target_gene & motif_id == target_motif, sort(unique(cancer))]
      if (only_significant_cancers && length(sig_cancers_this_pair) > 0) {
        cancers_to_plot <- sig_cancers_this_pair
      } else {
        cancers_to_plot <- sort(unique(as.character(plot_dt$cancer)))
      }
      plot_dt <- plot_dt[as.character(cancer) %in% cancers_to_plot]

      pair_stats <- anti_bc[gene == target_gene & motif_id == target_motif]
      if (nrow(pair_stats) > 0) {
        pair_stats <- unique(pair_stats[, .(cancer, gene, motif_id, n, pearson_r, pearson_p, pearson_fdr)])
      }

      missing_stat_cancers <- setdiff(cancers_to_plot, as.character(pair_stats$cancer))
      if (length(missing_stat_cancers)) {
        extra <- plot_dt[as.character(cancer) %in% missing_stat_cancers, {
          z <- safe_pearson(meth_percent, log2_expr, min_n = min_n_cor)
          .(n = z$n, pearson_r = z$r, pearson_p = z$p, pearson_fdr = NA_real_)
        }, by = .(cancer, gene, motif_id)]
        pair_stats <- rbindlist(list(pair_stats, extra), fill = TRUE)
      }

      base_name <- paste0(tf_name, "__", sprintf("%04d", rank_idx), "__", safe_name(target_gene), "__", safe_name(target_motif))
      out_pdf  <- file.path(out_dir_plot, paste0(base_name, "_report.pdf"))
      out_info <- file.path(out_dir_plot, paste0(base_name, "_summary.tsv"))

      plot_pair_report(plot_dt, pair_stats, sel_row, out_pdf, out_info, color_map_global, matched = FALSE)

      data.table(tf = tf_name, rank_idx = rank_idx, gene = target_gene, motif_id = target_motif,
                 status = "ok", error_message = NA_character_, pdf_file = out_pdf, summary_file = out_info)
    }, error = function(e) {
      data.table(tf = tf_name, rank_idx = as.integer(sel_row$rank_idx[1]),
                 gene = as.character(sel_row$gene[1]), motif_id = as.character(sel_row$motif_id[1]),
                 status = "error", error_message = conditionMessage(e),
                 pdf_file = NA_character_, summary_file = NA_character_)
    })
  }

  idx_all <- seq_len(nrow(top_pairs))
  batch_id <- ceiling(idx_all / batch_size_plot)
  batch_split <- split(idx_all, batch_id)
  res_list <- vector("list", length(batch_split))
  t0 <- proc.time()[3]

  for (b in seq_along(batch_split)) {
    tb <- proc.time()[3]
    idx <- batch_split[[b]]
    batch_dt <- top_pairs[idx]
    msg("BATCH", b, "/", length(batch_split), "| pairs:", nrow(batch_dt))

    batch_rows <- lapply(seq_len(nrow(batch_dt)), function(i) batch_dt[i])
    batch_res <- mclapply(batch_rows, process_one_pair, mc.cores = min(plot_cores, length(batch_rows)))
    res_list[[b]] <- rbindlist(batch_res, fill = TRUE)

    elapsed_batch <- proc.time()[3] - tb
    elapsed_total <- proc.time()[3] - t0
    avg_per_batch <- elapsed_total / b
    eta <- (length(batch_split) - b) * avg_per_batch
    msg("DONE | batch time:", fmt_sec(elapsed_batch), "| elapsed:", fmt_sec(elapsed_total), "| ETA:", fmt_sec(eta))
  }

  res <- rbindlist(res_list, fill = TRUE)
  fwrite(res, out_index, sep = "\t")

  msg("Saved:", top_index_file)
  msg("Saved:", out_index)
}

# ============================================================
# STEP 11
# ============================================================
build_interactive_html <- function(data_file, anti_files, out_dir_html) {
  sep_line(); msg("[Step 11] Build interactive HTML plots")

  dir.create(out_dir_html, recursive = TRUE, showWarnings = FALSE)

  dt <- fread(data_file)
  sample_col  <- intersect(names(dt), c("sample", "sample_id", "sample_barcode", "tumor_sample", "barcode"))
  patient_col <- intersect(names(dt), c("patient", "patient_id", "case_id", "submitter_id"))

  dt[, cancer := as.character(cancer)]
  dt[, gene := as.character(gene)]
  dt[, motif_id := as.character(motif_id)]
  dt[, sample_type := clean_sample_type(sample_type)]
  dt[, meth_beta := suppressWarnings(as.numeric(meth_beta))]
  dt[, expression := suppressWarnings(as.numeric(expression))]

  dt[, sample := if (length(sample_col)) as.character(get(sample_col[1])) else NA_character_]
  dt[, patient := if (length(patient_col)) as.character(get(patient_col[1])) else NA_character_]

  dt <- dt[
    !is.na(cancer) & nzchar(cancer) &
      !is.na(gene) & nzchar(gene) &
      !is.na(motif_id) & nzchar(motif_id) &
      sample_type %in% c("Healthy", "Tumor") &
      is.finite(meth_beta) & is.finite(expression)
  ]

  dt[, cancer := sub("^TCGA-", "", cancer)]
  dt[, log2_expr := log2(expression + 1)]
  if (max(dt$meth_beta, na.rm = TRUE) <= 1.5) dt[, meth_percent := meth_beta * 100] else dt[, meth_percent := meth_beta]

  if (all(is.na(dt$sample)) && !all(is.na(dt$patient))) dt[, sample := patient]
  if (all(is.na(dt$patient)) && !all(is.na(dt$sample))) dt[, patient := sample]

  dt[, sample_join := extract_tcga_sample(sample)]
  dt[, patient_join := extract_tcga_sample(patient)]

  expr_dt <- fread(tf_expr_file)
  expr_dt[, sample := as.character(sample)]
  expr_dt[, cancer := sub("^TCGA-", "", as.character(cancer))]
  expr_dt[, NRF1_tpm := suppressWarnings(as.numeric(NRF1_tpm))]
  expr_dt[, BANP_tpm := suppressWarnings(as.numeric(BANP_tpm))]
  expr_dt[, sample_join := extract_tcga_sample(sample)]
  expr_dt <- expr_dt[!is.na(sample_join)]
  expr_dt <- unique(expr_dt, by = "sample_join")

  tf_col <- paste0(tf_name, "_tpm")
  expr_sub <- expr_dt[, .(
    sample_join,
    tf_expression = get(tf_col),
    NRF1_tpm,
    BANP_tpm
  )]

  dt <- merge(dt, expr_sub, by = "sample_join", all.x = TRUE, sort = FALSE)
  dt[, tf_expr_log2 := fifelse(is.finite(tf_expression), log2(tf_expression + 1), NA_real_)]

  anti_all <- fread(anti_files$all)
  anti_bc  <- fread(anti_files$by_cancer)
  anti_bc[, cancer := sub("^TCGA-", "", as.character(cancer))]

  pairs_from_all <- unique(anti_all[, .(gene, motif_id)])
  pairs_from_bc  <- unique(anti_bc[, .(gene, motif_id)])
  pair_union     <- unique(rbindlist(list(pairs_from_all, pairs_from_bc), fill = TRUE))

  all_sum <- anti_all[, .(
    pooled_n = n_all[1],
    pooled_r = pearson_r_all[1],
    pooled_p = pearson_p_all[1],
    pooled_fdr = pearson_fdr_all[1]
  ), by = .(gene, motif_id)]

  bc_sum <- anti_bc[, .(
    n_sig_cancers = uniqueN(cancer),
    sig_cancers = paste(sort(unique(cancer)), collapse = ",")
  ), by = .(gene, motif_id)]

  top_pairs <- merge(pair_union, all_sum, by = c("gene", "motif_id"), all.x = TRUE)
  top_pairs <- merge(top_pairs, bc_sum, by = c("gene", "motif_id"), all.x = TRUE)
  top_pairs[is.na(n_sig_cancers), n_sig_cancers := 0L]
  top_pairs[is.na(sig_cancers), sig_cancers := ""]
  top_pairs[, rank_fdr := pooled_fdr]
  top_pairs[, rank_r := pooled_r]
  top_pairs <- top_pairs[order(rank_fdr, rank_r, -n_sig_cancers, gene, motif_id)]
  top_pairs[, rank_idx := seq_len(.N)]

  color_map_global <- load_color_map(color_file)

  process_one_pair <- function(sel_row) {
    tryCatch({
      target_gene  <- as.character(sel_row$gene[1])
      target_motif <- as.character(sel_row$motif_id[1])
      rank_idx     <- as.integer(sel_row$rank_idx[1])

      plot_dt <- dt[gene == target_gene & motif_id == target_motif]
      if (nrow(plot_dt) == 0) {
        return(data.table(tf = tf_name, rank_idx = rank_idx, gene = target_gene, motif_id = target_motif,
                          status = "no_rows", error_message = NA_character_, html_file = NA_character_))
      }

      sig_cancers_this_pair <- anti_bc[gene == target_gene & motif_id == target_motif, sort(unique(cancer))]
      if (only_significant_cancers && length(sig_cancers_this_pair) > 0) {
        cancers_to_plot <- sig_cancers_this_pair
      } else {
        cancers_to_plot <- sort(unique(as.character(plot_dt$cancer)))
      }

      plot_dt <- plot_dt[as.character(cancer) %in% cancers_to_plot]
      if (nrow(plot_dt) == 0) {
        return(data.table(tf = tf_name, rank_idx = rank_idx, gene = target_gene, motif_id = target_motif,
                          status = "no_rows_after_cancer_filter", error_message = NA_character_, html_file = NA_character_))
      }

      cancers_present <- sort(unique(as.character(plot_dt$cancer)))
      col_sub <- data.table(cancer = cancers_present)
      col_sub[, color := color_map_global[cancer]]
      col_sub[is.na(color), color := "grey50"]
      color_map <- setNames(col_sub$color, col_sub$cancer)

      plot_dt[, cancer := factor(as.character(cancer), levels = names(color_map))]
      plot_dt[, sample_type := factor(as.character(sample_type), levels = c("Healthy", "Tumor"))]

      global_cor <- safe_pearson(plot_dt$meth_percent, plot_dt$log2_expr, min_n = min_n_cor)

      base_name <- paste0(tf_name, "__", sprintf("%04d", rank_idx), "__", safe_name(target_gene), "__", safe_name(target_motif))
      out_html <- file.path(out_dir_html, paste0(base_name, "_interactive_target_vs_methylation.html"))

      plot_dt[, hover_text := paste0(
        "TF: ", tf_name,
        "<br>Target gene: ", gene,
        "<br>Motif: ", motif_id,
        "<br>Sample: ", fifelse(is.na(sample_join), "NA", sample_join),
        "<br>Patient: ", fifelse(is.na(patient), "NA", patient),
        "<br>Cancer: ", as.character(cancer),
        "<br>Sample type: ", as.character(sample_type),
        "<br>Methylation (%): ", ifelse(is.finite(meth_percent), sprintf('%.2f', meth_percent), "NA"),
        "<br>Target expression: ", ifelse(is.finite(expression), sprintf('%.3f', expression), "NA"),
        "<br>Target expression log2(x+1): ", ifelse(is.finite(log2_expr), sprintf('%.3f', log2_expr), "NA"),
        "<br>", tf_name, " TPM: ", ifelse(is.finite(tf_expression), sprintf('%.3f', tf_expression), "NA"),
        "<br>", tf_name, " log2(TPM+1): ", ifelse(is.finite(tf_expr_log2), sprintf('%.3f', tf_expr_log2), "NA"),
        "<br>NRF1 TPM: ", ifelse(is.finite(NRF1_tpm), sprintf('%.3f', NRF1_tpm), "NA"),
        "<br>BANP TPM: ", ifelse(is.finite(BANP_tpm), sprintf('%.3f', BANP_tpm), "NA")
      )]

      plot_df <- as.data.frame(plot_dt)
      p1 <- plotly::plot_ly(
        data = plot_df,
        x = ~meth_percent,
        y = ~log2_expr,
        type = "scatter",
        mode = "markers",
        color = ~cancer,
        colors = unname(color_map),
        symbol = ~sample_type,
        symbols = c("circle", "triangle-up"),
        text = ~hover_text,
        hoverinfo = "text",
        marker = list(size = 7, opacity = 0.70)
      )

      p1 <- plotly::layout(
        p1,
        title = list(text = paste0(
          "Interactive anti-correlated scatter: ", tf_name, " | ", target_gene, " | ", target_motif,
          "<br><sup>n=", global_cor$n,
          " | Pearson r=", ifelse(is.na(global_cor$r), "NA", sprintf("%.3f", global_cor$r)),
          " | p=", fmt_p(global_cor$p), "</sup>"
        )),
        xaxis = list(title = "Motif methylation (%)"),
        yaxis = list(title = "log2(target expression + 1)")
      )

      htmlwidgets::saveWidget(p1, out_html, selfcontained = TRUE)

      data.table(tf = tf_name, rank_idx = rank_idx, gene = target_gene, motif_id = target_motif,
                 status = "ok", error_message = NA_character_, html_file = out_html)
    }, error = function(e) {
      data.table(tf = tf_name, rank_idx = as.integer(sel_row$rank_idx[1]),
                 gene = as.character(sel_row$gene[1]),
                 motif_id = as.character(sel_row$motif_id[1]),
                 status = "error", error_message = conditionMessage(e), html_file = NA_character_)
    })
  }

  idx_all <- seq_len(nrow(top_pairs))
  batch_id <- ceiling(idx_all / 1L)
  batch_split <- split(idx_all, batch_id)
  res_list <- vector("list", length(batch_split))
  t0 <- proc.time()[3]

  for (b in seq_along(batch_split)) {
    tb <- proc.time()[3]
    idx <- batch_split[[b]]
    batch_dt <- top_pairs[idx]
    batch_rows <- lapply(seq_len(nrow(batch_dt)), function(i) batch_dt[i])
    batch_res <- lapply(batch_rows, process_one_pair)
    res_list[[b]] <- rbindlist(batch_res, fill = TRUE)

    elapsed_batch <- proc.time()[3] - tb
    elapsed_total <- proc.time()[3] - t0
    avg_per_batch <- elapsed_total / b
    eta <- (length(batch_split) - b) * avg_per_batch
    msg("BATCH", b, "/", length(batch_split), "| batch time:", fmt_sec(elapsed_batch), "| elapsed:", fmt_sec(elapsed_total), "| ETA:", fmt_sec(eta))
  }

  res <- rbindlist(res_list, fill = TRUE)
  out_index <- file.path(out_dir_html, paste0("interactive_html_index_", tf_name, ".tsv"))
  fwrite(res, out_index, sep = "\t")
  msg("Saved:", out_index)
}

# ============================================================
# MASTER PIPELINE
# ============================================================
run_pipeline <- function() {
  sep_line()
  msg("Running pipeline for TF:", tf_name)
  msg("Motif label:", tf_motif_label)
  msg("Interactive HTML:", make_interactive_html)
  msg("Reuse existing intermediate files:", reuse_existing_intermediate_files)
  sep_line()

  motif_probe_file <- build_motif_probe_map()
  meth_annot_file  <- build_methylation_sample_annotation()
  meth_file        <- compute_motif_methylation(motif_probe_file, meth_annot_file)
  merged_file      <- merge_with_expression(meth_file)

  # ==========================================================
  # ALL SAMPLES
  # ==========================================================
  corr_all <- run_correlations(
    in_file = merged_file,
    out_dir_corr = file.path(all_dir, "correlation_stats"),
    tumor_only = FALSE
  )

  plot_top_pairs(
    data_file = merged_file,
    stats_file = corr_all$by_cancer,
    out_dir_plot = file.path(all_dir, "top_pair_plots"),
    matched = FALSE
  )

  anti_files_all <- extract_anti_correlated(
    corr_files = corr_all,
    out_dir_corr = file.path(all_dir, "correlation_stats"),
    data_file = merged_file
  )

  plot_all_anti_correlated(
    data_file = merged_file,
    anti_files = anti_files_all,
    out_dir_plot = file.path(all_dir, "anti_correlated_pair_plots")
  )

  if (make_interactive_html) {
    build_interactive_html(
      data_file = merged_file,
      anti_files = anti_files_all,
      out_dir_html = file.path(all_dir, "interactive_target_vs_methylation_only")
    )
  } else {
    msg("Skipping interactive HTML plots for all samples. Set make_interactive_html <- TRUE to generate them.")
  }

  # ==========================================================
  # MATCHED PATIENTS
  # ==========================================================
  matched_file <- build_matched_subset(merged_file)

  corr_matched <- run_correlations(
    in_file = matched_file,
    out_dir_corr = file.path(matched_dir, "correlation_stats"),
    tumor_only = FALSE
  )

  plot_top_pairs(
    data_file = matched_file,
    stats_file = corr_matched$by_cancer,
    out_dir_plot = file.path(matched_dir, "top_pair_plots"),
    matched = TRUE
  )

  anti_files_matched <- extract_anti_correlated(
    corr_files = corr_matched,
    out_dir_corr = file.path(matched_dir, "correlation_stats"),
    data_file = matched_file
  )

  plot_all_anti_correlated(
    data_file = matched_file,
    anti_files = anti_files_matched,
    out_dir_plot = file.path(matched_dir, "anti_correlated_pair_plots")
  )

  if (make_interactive_html) {
    build_interactive_html(
      data_file = matched_file,
      anti_files = anti_files_matched,
      out_dir_html = file.path(matched_dir, "interactive_target_vs_methylation_only")
    )
  } else {
    msg("Skipping interactive HTML plots for matched patients. Set make_interactive_html <- TRUE to generate them.")
  }

  sep_line()
  msg("Pipeline finished for TF:", tf_name)
  msg("Main output directory:", base_dir)
  sep_line()
}

run_pipeline()

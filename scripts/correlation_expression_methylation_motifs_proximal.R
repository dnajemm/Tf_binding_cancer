
'''This script does the following steps:
Step 1: Built the proximal BANP motif–probe map from probe_gene_pairs_in_motifs_with_tss.tsv.gz by keeping only is_proximal == TRUE, restricting to TF == "BANP_mm0to2_noCGmm", and retaining unique gene, motif_id, and probe combinations.
Step 2: Built the HM450 methylation sample annotation for the 3d+4d_noOV cohort by listing methylation files, extracting sample_barcode, patient_id, sample_type (01 = Tumor, 11 = Healthy), and cancer, then keeping only samples whose patients are present in the cohort file.
Step 3: Computed sample-level motif methylation for each BANP proximal gene–motif pair by reading each sample’s HM450 beta values, restricting to mapped probes, and defining motif methylation as the MAX(beta) across linked probes; also saved the selected max probe(s), compact probe summaries, and per-sample processing status.
Step 4: Merged the motif methylation table with the gene expression matrix by reshaping expression to long format, extracting sample barcodes from expression sample names, and joining on gene + sample_barcode to create a combined methylation–expression table.
Step 5: Computed Pearson correlations between motif methylation and gene expression both across all cancers pooled and within each cancer type, then applied BH/FDR correction to identify significant associations.
Step 8: Selected the top 50 unique BANP gene–motif pairs from the cancer-specific correlation results and generated multi-page PDF reports for each pair, including pooled scatterplots, per-cancer scatterplots, and Healthy-vs-Tumor dual-panel summaries for both expression and methylation.
'''

##############################################################
# Step 1 : Build the proximal BANP motif–probe map
# input : file with probes in motifs with their TSS distance 
# ./methylation/probe_gene_pairs_in_motifs_with_tss.tsv.gz
# output : file with proximal BANP motif–probe pairs
############################################################## 
Rscript - <<'EOF'
suppressPackageStartupMessages(library(data.table))

in_file  <- "./methylation/probe_gene_pairs_in_motifs_with_tss.tsv.gz"
out_file <- "./results/methylation/correlation_expression_methylation/tumor_vs_healthy/BANP/BANP_proximal_motif_probe_map_long.tsv.gz"

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

dt <- fread(in_file)

# keep only proximal rows
dt <- dt[is_proximal == TRUE]

# keep only BANP if needed
dt <- dt[TF == "BANP_mm0to2_noCGmm"]

# keep only needed columns
out <- unique(dt[, .(
  gene = as.character(gene),
  motif_id = as.character(motif_instance_id),
  probe = as.character(probe)
)])

# remove missing / empty values
out <- out[
  !is.na(gene) & nzchar(gene) &
  !is.na(motif_id) & nzchar(motif_id) &
  !is.na(probe) & nzchar(probe)
]

setorder(out, gene, motif_id, probe)

fwrite(out, out_file, sep = "\t")

cat("Saved:", out_file, "\n")
cat("Rows:", nrow(out), "\n")
cat("Unique genes:", uniqueN(out$gene), "\n")
cat("Unique motifs:", uniqueN(out$motif_id), "\n")
cat("Unique probes:", uniqueN(out$probe), "\n")
cat("Unique gene-motif pairs:", uniqueN(out[, paste(gene, motif_id, sep='||')]), "\n")
EOF

######################################################################################
# Step 2 : Build methylation sample annotation for the 3d+4d_noOV cohort
# input : directory with methylation files, cohort file with patient IDs
# output : file with sample annotation for methylation files in the cohort
######################################################################################
Rscript - <<'EOF'
suppressPackageStartupMessages(library(data.table))

# =========================
# INPUTS
# =========================
meth_dir     <- "./methylation/filtered_methylation"
cohort_file  <- "./results/multi_omics/samples_3d+4d_noOV.tsv"
out_file     <- "./results/methylation/correlation_expression_methylation/tumor_vs_healthy/BANP/methylation_sample_annotation_3d4d_noOV.tsv.gz"

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

# =========================
# HELPERS
# =========================
extract_sample_barcode <- function(x) {
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

extract_sample_type <- function(x) {
  x <- as.character(x)
  code <- sub("^TCGA-[^-]+-[^-]+-([0-9]{2})[A-Z]$", "\\1", x)
  fifelse(code == "01", "Tumor",
    fifelse(code == "11", "Healthy", NA_character_))
}

# =========================
# LOAD COHORT
# =========================
cat("[1/4] Loading cohort file...\n")
cohort <- fread(cohort_file, header = FALSE)

# first column = patient ID
patients_keep <- unique(as.character(cohort[[1]]))
patients_keep <- patients_keep[grepl("^TCGA-[^-]+-[^-]+$", patients_keep)]

cat("Patients in cohort:", length(patients_keep), "\n")

# =========================
# LIST METHYLATION FILES
# =========================
cat("[2/4] Listing methylation files...\n")
files <- list.files(meth_dir, pattern = "\\.bed\\.gz$", full.names = TRUE)

dt <- data.table(
  file = files,
  fname = basename(files)
)

dt[, sample_barcode := extract_sample_barcode(fname)]
dt[, patient_id := extract_patient_id(sample_barcode)]
dt[, sample_type := extract_sample_type(sample_barcode)]
dt[, cancer := sub("^HM450_TCGA-([^-]+)-.*$", "\\1", fname)]

dt <- dt[
  !is.na(sample_barcode) &
  !is.na(patient_id) &
  !is.na(sample_type)
]

cat("All valid methylation samples:", nrow(dt), "\n")

# =========================
# RESTRICT TO COHORT
# =========================
cat("[3/4] Restricting to 3d+4d_noOV cohort...\n")
dt <- dt[patient_id %in% patients_keep]

dt <- unique(dt[, .(patient_id, sample_barcode, sample_type, cancer, file)])
setorder(dt, sample_type, cancer, patient_id, sample_barcode)

cat("Samples kept:", nrow(dt), "\n")
cat("Tumor:", dt[sample_type == "Tumor", .N], "\n")
cat("Healthy:", dt[sample_type == "Healthy", .N], "\n")

# =========================
# WRITE
# =========================
cat("[4/4] Writing output...\n")
fwrite(dt, out_file, sep = "\t")

cat("Done.\n")
cat("Saved:", out_file, "\n")
EOF

###############################################################################################################
# Step 3 : Compute motif methylation per sample using MAX(beta)
# input : proximal motif-probe map, methylation sample annotation, methylation files
# output : motif-sample methylation values using MAX(beta) across probes, plus probe-level details and status of each sample
# It runs with parallelization in batches of 20 samples to speed up processing 
###############################################################################################################
Rscript - <<'EOF'
suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
})

# =========================
# INPUTS
# =========================
motif_probe_file <- "./results/methylation/correlation_expression_methylation/tumor_vs_healthy/BANP/BANP_proximal_motif_probe_map_long.tsv.gz"
meth_annot_file  <- "./results/methylation/correlation_expression_methylation/tumor_vs_healthy/BANP/methylation_sample_annotation_3d4d_noOV.tsv.gz"

out_dir          <- "./results/methylation/correlation_expression_methylation/tumor_vs_healthy/BANP"
out_file_main    <- file.path(out_dir, "motif_sample_methylation_MAX_3d4d_noOV.tsv.gz")
out_file_probes  <- file.path(out_dir, "motif_sample_methylation_MAX_selected_probes_3d4d_noOV.tsv.gz")
out_file_compact <- file.path(out_dir, "motif_sample_methylation_MAX_selected_probes_compact_3d4d_noOV.tsv.gz")
out_file_status  <- file.path(out_dir, "motif_sample_methylation_MAX_status_3d4d_noOV.tsv.gz")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

mc_cores   <- 20L
batch_size <- 20L

# =========================
# LOAD MAP
# =========================
cat("[1/5] Loading proximal motif-probe map...\n")
mp <- fread(motif_probe_file)

need_mp <- c("gene", "motif_id", "probe")
miss_mp <- setdiff(need_mp, names(mp))
if (length(miss_mp)) stop("Missing columns in motif_probe_file: ", paste(miss_mp, collapse = ", "))

mp[, `:=`(
  gene = as.character(gene),
  motif_id = as.character(motif_id),
  probe = as.character(probe)
)]

mp <- unique(mp[
  !is.na(gene) & nzchar(gene) &
  !is.na(motif_id) & nzchar(motif_id) &
  !is.na(probe) & nzchar(probe)
])

probe_keep <- unique(mp$probe)

cat("Rows in map:", nrow(mp), "\n")
cat("Unique gene-motif pairs:", uniqueN(mp[, paste(gene, motif_id, sep = '||')]), "\n")
cat("Unique probes:", length(probe_keep), "\n")

# =========================
# LOAD SAMPLE ANNOTATION
# =========================
cat("[2/5] Loading methylation sample annotation...\n")
sa <- fread(meth_annot_file)

need_sa <- c("sample_barcode", "sample_type", "cancer", "file")
miss_sa <- setdiff(need_sa, names(sa))
if (length(miss_sa)) stop("Missing columns in meth_annot_file: ", paste(miss_sa, collapse = ", "))

sa[, `:=`(
  sample_barcode = as.character(sample_barcode),
  sample_type = as.character(sample_type),
  cancer = as.character(cancer),
  file = as.character(file)
)]

sa <- sa[file.exists(file)]

if (nrow(sa) == 0) stop("No existing methylation files found.")
cat("Samples to process:", nrow(sa), "\n")
cat("Parallel workers:", mc_cores, "\n")
cat("Batch size:", batch_size, "\n")

# =========================
# WORKER
# =========================
process_one_sample <- function(i) {
  f <- sa$file[i]
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

  main_i[, `:=`(
    sample_barcode = sb,
    sample_type = st,
    cancer = ca
  )]

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

# =========================
# RUN IN BATCHES OF 20
# =========================
cat("[3/5] Computing motif methylation per sample using MAX(beta)...\n")

idx_all <- seq_len(nrow(sa))
batch_id <- ceiling(idx_all / batch_size)
batch_split <- split(idx_all, batch_id)

main_chunks   <- vector("list", length(batch_split))
probe_chunks  <- vector("list", length(batch_split))
status_chunks <- vector("list", length(batch_split))

for (b in seq_along(batch_split)) {
  idx <- batch_split[[b]]
  cat("  batch", b, "/", length(batch_split), "| samples:", length(idx), "\n")

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
setorder(probe_res, gene, motif_id, sample_type, cancer, sample_barcode, probe)
setorder(probe_summary_compact, gene, motif_id, sample_type, cancer, sample_barcode)

# =========================
# WRITE
# =========================
cat("[4/5] Writing outputs...\n")
fwrite(main_res, out_file_main, sep = "\t")
fwrite(probe_res, out_file_probes, sep = "\t")
fwrite(probe_summary_compact, out_file_compact, sep = "\t")
fwrite(status_dt, out_file_status, sep = "\t")

# =========================
# SUMMARY
# =========================
cat("[5/5] Done.\n")
cat("Main output   :", out_file_main, "\n")
cat("Probe output  :", out_file_probes, "\n")
cat("Compact output:", out_file_compact, "\n")
cat("Status output :", out_file_status, "\n")
cat("Main rows:", nrow(main_res), "\n")
cat("Probe rows:", nrow(probe_res), "\n")
cat("Unique samples:", uniqueN(main_res$sample_barcode), "\n")
cat("Rows with >1 tied max probe:", main_res[n_max_probes > 1, .N], "\n")
cat("Successful samples:", status_dt[ok == TRUE, .N], "\n")
cat("Failed samples:", status_dt[ok == FALSE, .N], "\n")
EOF

################################################
# STEP 4: MERGE METHYLATION WITH EXPRESSION
################################################
Rscript - <<'EOF'
suppressPackageStartupMessages(library(data.table))

# Input from step 3
meth_file <- "./results/methylation/correlation_expression_methylation/tumor_vs_healthy/BANP/motif_sample_methylation_MAX_3d4d_noOV.tsv.gz"

# Expression matrix: rows = samples, columns = genes
expr_file <- "./expression/gene_expression_matrix_3d4d_noOV.tsv"

# Output
out_file  <- "./results/methylation/correlation_expression_methylation/tumor_vs_healthy/BANP/motif_sample_expression_methylation_all.tsv.gz"

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

# =========================
# HELPERS
# =========================

# From expression row name like:
# TCGA-ACC-TCGA-OR-A5J1-01A_1
# extract:
#   cancer = ACC
extract_expr_cancer <- function(x) {
  x <- as.character(x)
  out <- sub("^TCGA-([^-]+)-TCGA-[^-]+-[^-]+-[0-9]{2}[A-Z](?:_.*)?$", "\\1", x)
  out[!grepl("^[A-Z0-9]+$", out)] <- NA_character_
  out
}

# From expression row name like:
# TCGA-ACC-TCGA-OR-A5J1-01A_1
# extract:
#   sample_barcode = TCGA-OR-A5J1-01A
extract_expr_sample_barcode <- function(x) {
  x <- as.character(x)
  out <- sub("^TCGA-[^-]+-(TCGA-[^-]+-[^-]+-[0-9]{2}[A-Z])(?:_.*)?$", "\\1", x)
  out[!grepl("^TCGA-[^-]+-[^-]+-[0-9]{2}[A-Z]$", out)] <- NA_character_
  out
}

extract_sample_type <- function(x) {
  x <- as.character(x)
  code <- sub("^TCGA-[^-]+-[^-]+-([0-9]{2})[A-Z]$", "\\1", x)
  fifelse(code == "01", "Tumor",
    fifelse(code == "11", "Healthy", NA_character_))
}

# =========================
# LOAD METHYLATION
# =========================
cat("[1/5] Loading methylation table...\n")
meth <- fread(meth_file)

need_meth <- c("gene", "motif_id", "sample_barcode", "sample_type", "cancer", "meth_beta")
miss_meth <- setdiff(need_meth, names(meth))
if (length(miss_meth)) stop("Missing columns in meth_file: ", paste(miss_meth, collapse = ", "))

meth[, `:=`(
  gene = as.character(gene),
  motif_id = as.character(motif_id),
  sample_barcode = as.character(sample_barcode),
  sample_type = as.character(sample_type),
  cancer = as.character(cancer),
  meth_beta = as.numeric(meth_beta)
)]

cat("Methylation rows:", nrow(meth), "\n")
cat("Unique genes in methylation:", uniqueN(meth$gene), "\n")
cat("Unique samples in methylation:", uniqueN(meth$sample_barcode), "\n")

# =========================
# LOAD EXPRESSION MATRIX
# =========================
cat("[2/5] Loading expression matrix...\n")
expr <- fread(expr_file)

if (!"sample" %in% names(expr)) {
  stop("Expression file must contain a column named 'sample'")
}

# Keep only genes that appear in methylation table
genes_keep <- unique(meth$gene)
expr_cols_keep <- intersect(names(expr), c("sample", genes_keep))

if (length(expr_cols_keep) <= 1) {
  stop("No gene columns from methylation table found in expression matrix.")
}

expr <- expr[, ..expr_cols_keep]

cat("Expression rows:", nrow(expr), "\n")
cat("Genes kept from expression matrix:", length(setdiff(expr_cols_keep, "sample")), "\n")

# =========================
# RESHAPE EXPRESSION TO LONG FORMAT
# =========================
cat("[3/5] Reshaping expression matrix to long format...\n")

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
expr_long[, sample_type := extract_sample_type(sample_barcode)]

expr_long <- expr_long[
  !is.na(sample_barcode) &
  !is.na(gene) &
  !is.na(expression)
]

cat("Long expression rows:", nrow(expr_long), "\n")
cat("Unique samples in expression long:", uniqueN(expr_long$sample_barcode), "\n")

# =========================
# MERGE METHYLATION + EXPRESSION
# =========================
cat("[4/5] Merging methylation with expression...\n")

merged <- merge(
  meth,
  expr_long[, .(gene, sample_barcode, expression)],
  by = c("gene", "sample_barcode"),
  all = FALSE
)

# Optional consistency checks
merged <- merged[!is.na(meth_beta) & !is.na(expression)]

setcolorder(merged, c(
  "gene", "motif_id", "sample_barcode", "sample_type", "cancer",
  "meth_beta", "expression", "n_probes_found", "n_max_probes"
))
setorder(merged, gene, motif_id, sample_type, cancer, sample_barcode)

cat("Merged rows:", nrow(merged), "\n")
cat("Unique merged genes:", uniqueN(merged$gene), "\n")
cat("Unique merged samples:", uniqueN(merged$sample_barcode), "\n")

# =========================
# WRITE OUTPUT
# =========================
cat("[5/5] Writing output...\n")
fwrite(merged, out_file, sep = "\t")

cat("Done.\n")
cat("Saved:", out_file, "\n")
EOF

# =========================
# STEP 5: CORRELATION
#   - r_all: across all cancer types pooled
#   - r_by_cancer: within each cancer type
# =========================
Rscript - <<'EOF'
suppressPackageStartupMessages(library(data.table))

in_file <- "./results/methylation/correlation_expression_methylation/tumor_vs_healthy/BANP/motif_sample_expression_methylation_all.tsv.gz"
out_dir <- "./results/methylation/correlation_expression_methylation/tumor_vs_healthy/BANP/correlation_stats"

out_file_all       <- file.path(out_dir, "pearson_correlation_all_pairs_all_cancers.tsv.gz")
out_file_by_cancer <- file.path(out_dir, "pearson_correlation_all_pairs_by_cancer.tsv.gz")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

min_samples <- 3L

# Optional:
# set to TRUE if you want tumor-only analysis
tumor_only <- FALSE

# =========================
# LOAD
# =========================
cat("[1/5] Loading merged methylation-expression table...\n")
dt <- fread(in_file)

need <- c("gene", "motif_id", "sample_barcode", "sample_type", "cancer", "meth_beta", "expression")
miss <- setdiff(need, names(dt))
if (length(miss)) stop("Missing columns in input file: ", paste(miss, collapse = ", "))

dt[, `:=`(
  gene = as.character(gene),
  motif_id = as.character(motif_id),
  sample_barcode = as.character(sample_barcode),
  sample_type = as.character(sample_type),
  cancer = as.character(cancer),
  meth_beta = suppressWarnings(as.numeric(meth_beta)),
  expression = suppressWarnings(as.numeric(expression))
)]

dt <- dt[!is.na(gene) & !is.na(motif_id) & !is.na(sample_barcode) &
         !is.na(cancer) & !is.na(meth_beta) & !is.na(expression)]

if (tumor_only) {
  dt <- dt[sample_type == "Tumor"]
}

cat("Rows:", nrow(dt), "\n")
cat("Unique gene-motif pairs:", uniqueN(dt[, paste(gene, motif_id, sep = "||")]), "\n")
cat("Unique cancers:", uniqueN(dt$cancer), "\n")
cat("Unique samples:", uniqueN(dt$sample_barcode), "\n")

# =========================
# HELPER
# =========================
safe_pearson <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]

  n <- length(x)

  # Need enough samples and non-zero variance
  if (n < min_samples) {
    return(list(n = n, r = NA_real_, p = NA_real_))
  }
  if (length(unique(x)) < 2 || length(unique(y)) < 2) {
    return(list(n = n, r = NA_real_, p = NA_real_))
  }

  ct <- tryCatch(
    cor.test(x, y, method = "pearson"),
    error = function(e) NULL
  )

  if (is.null(ct)) {
    return(list(n = n, r = NA_real_, p = NA_real_))
  }

  list(
    n = n,
    r = unname(ct$estimate),
    p = ct$p.value
  )
}

# =========================
# r_all: ACROSS ALL CANCERS POOLED
# =========================
cat("[2/5] Computing pooled correlations (r_all)...\n")

res_all <- dt[, {
  z <- safe_pearson(meth_beta, expression)
  .(
    n_all = z$n,
    pearson_r_all = z$r,
    pearson_p_all = z$p
  )
}, by = .(gene, motif_id)]

res_all[, pearson_fdr_all := p.adjust(pearson_p_all, method = "BH")]
setorder(res_all, pearson_p_all, gene, motif_id)

cat("Rows in pooled result:", nrow(res_all), "\n")

# =========================
# r_by_cancer: WITHIN EACH CANCER
# =========================
cat("[3/5] Computing cancer-specific correlations...\n")

res_by_cancer <- dt[, {
  z <- safe_pearson(meth_beta, expression)
  .(
    n = z$n,
    pearson_r = z$r,
    pearson_p = z$p
  )
}, by = .(gene, motif_id, cancer)]

# FDR within all cancer-specific tests together
res_by_cancer[, pearson_fdr := p.adjust(pearson_p, method = "BH")]

setorder(res_by_cancer, cancer, pearson_p, gene, motif_id)

cat("Rows in cancer-specific result:", nrow(res_by_cancer), "\n")

# =========================
# WRITE
# =========================
cat("[4/5] Writing outputs...\n")
fwrite(res_all, out_file_all, sep = "\t")
fwrite(res_by_cancer, out_file_by_cancer, sep = "\t")

# =========================
# SUMMARY
# =========================
cat("[5/5] Done.\n")
cat("Saved pooled correlations   :", out_file_all, "\n")
cat("Saved cancer correlations   :", out_file_by_cancer, "\n")
cat("Significant pooled (FDR<0.05):", res_all[!is.na(pearson_fdr_all) & pearson_fdr_all < 0.05, .N], "\n")
cat("Significant by cancer (FDR<0.05):", res_by_cancer[!is.na(pearson_fdr) & pearson_fdr < 0.05, .N], "\n")
EOF

# Step 7 : plot the top 50 unique NRF1 gene–motif pairs from the cancer-specific correlation results, generating multi-page PDF reports for each pair with pooled scatterplots, per-cancer scatterplots, and Healthy-vs-Tumor dual-panel summaries for both expression and methylation.
Rscript - <<'EOF'
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
  library(parallel)
})

# =========================
# INPUTS
# =========================
data_file   <- "./results/methylation/correlation_expression_methylation/tumor_vs_healthy/NRF1/motif_sample_expression_methylation_all.tsv.gz"
stats_file  <- "./results/methylation/correlation_expression_methylation/tumor_vs_healthy/NRF1/correlation_stats/pearson_correlation_all_pairs_by_cancer.tsv.gz"
color_file  <- "./results/multi_omics/cancer_color_order_with_defined_colours.tsv"
out_dir     <- "./results/methylation/correlation_expression_methylation/tumor_vs_healthy/NRF1/top_pair_plots"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

n_top        <- 50L
mc_cores     <- 5L
batch_size   <- 5L
min_n_cor    <- 3L
min_n_wilcox <- 2L

top_index_file <- file.path(out_dir, "top50_selected_pairs.tsv")
out_index      <- file.path(out_dir, "top50_generated_reports.tsv")

# =========================
# HELPERS
# =========================
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

safe_name <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

safe_cor_test <- function(x, y, min_n = 3L) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]

  n <- length(x)
  if (n < min_n) return(list(r = NA_real_, p = NA_real_, n = n))
  if (length(unique(x)) < 2 || length(unique(y)) < 2) {
    return(list(r = NA_real_, p = NA_real_, n = n))
  }

  ct <- tryCatch(cor.test(x, y, method = "pearson"), error = function(e) NULL)
  if (is.null(ct)) return(list(r = NA_real_, p = NA_real_, n = n))

  list(r = unname(ct$estimate), p = ct$p.value, n = n)
}

safe_wilcox_hc <- function(values, groups, min_n = 2L) {
  keep <- is.finite(values) & !is.na(groups)
  values <- values[keep]
  groups <- trimws(as.character(groups[keep]))

  groups[tolower(groups) %in% c("healthy", "normal", "solid tissue normal", "adjacent normal")] <- "Healthy"
  groups[tolower(groups) %in% c("tumor", "primary tumor", "primary tumour", "cancer")] <- "Tumor"

  keep2 <- groups %in% c("Healthy", "Tumor")
  values <- values[keep2]
  groups <- groups[keep2]

  n_h <- sum(groups == "Healthy")
  n_c <- sum(groups == "Tumor")

  if (n_h < min_n || n_c < min_n) {
    return(list(p = NA_real_, n_h = n_h, n_c = n_c, status = "insufficient_n"))
  }

  wt <- tryCatch(
    suppressWarnings(wilcox.test(values ~ groups, exact = FALSE)),
    error = function(e) NULL
  )

  if (is.null(wt)) {
    return(list(p = NA_real_, n_h = n_h, n_c = n_c, status = "wilcox_error"))
  }

  list(p = wt$p.value, n_h = n_h, n_c = n_c, status = "ok")
}

fmt_sec <- function(x) {
  x <- as.numeric(x)
  if (!is.finite(x)) return("NA")
  if (x < 60) return(paste0(round(x, 1), " sec"))
  paste0(round(x / 60, 1), " min")
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
    sample_type,
    cancer,
    measure = "Expression",
    value_raw = log2_expr,
    value_plot = log2_expr
  )]

  meth_dt <- d[, .(
    sample_type,
    cancer,
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
        measure == "Methylation" & sample_type == "Healthy", 4,
        5
      )
    )
  )]

  out[, x_lab := fifelse(x_num %in% c(1, 4), "H", "C")]

  list(
    dt = out,
    y1_min = y1_min,
    y1_max = y1_max,
    y2_min = y2_min,
    y2_max = y2_max,
    scale_factor = scale_factor
  )
}

# =========================
# LOAD DATA
# =========================
cat("[1/7] Loading merged sample-level table...\n")
t0 <- proc.time()[3]
dt <- fread(data_file)

need_dt <- c("cancer", "gene", "motif_id", "sample_type", "meth_beta", "expression")
miss_dt <- setdiff(need_dt, names(dt))
if (length(miss_dt)) stop("Missing columns in data_file: ", paste(miss_dt, collapse = ", "))

dt[, `:=`(
  cancer      = as.character(cancer),
  gene        = as.character(gene),
  motif_id    = as.character(motif_id),
  sample_type = trimws(as.character(sample_type)),
  meth_beta   = suppressWarnings(as.numeric(meth_beta)),
  expression  = suppressWarnings(as.numeric(expression))
)]

dt <- dt[
  !is.na(cancer) & nzchar(cancer) &
  !is.na(gene) & nzchar(gene) &
  !is.na(motif_id) & nzchar(motif_id) &
  !is.na(sample_type) & nzchar(sample_type) &
  is.finite(meth_beta) &
  is.finite(expression)
]

dt[, cancer := sub("^TCGA-", "", cancer)]

dt[tolower(sample_type) %in% c("healthy", "normal", "solid tissue normal", "adjacent normal"), sample_type := "Healthy"]
dt[tolower(sample_type) %in% c("tumor", "primary tumor", "primary tumour", "cancer"), sample_type := "Tumor"]

cat("    sample_type counts before final filter:\n")
print(dt[, .N, by = sample_type][order(-N)])

dt <- dt[sample_type %in% c("Healthy", "Tumor")]

cat("    sample_type counts after final filter:\n")
print(dt[, .N, by = sample_type][order(-N)])

dt[, log2_expr := log2(expression + 1)]

if (max(dt$meth_beta, na.rm = TRUE) <= 1.5) {
  dt[, meth_percent := meth_beta * 100]
  cat("    methylation detected on 0-1 scale -> converted to 0-100 for plotting\n")
} else {
  dt[, meth_percent := meth_beta]
  cat("    methylation detected on 0-100 scale\n")
}

cat("    rows kept:", nrow(dt), "\n")

cat("[2/7] Loading Step 5 correlation table...\n")
st <- fread(stats_file)

need_st <- c("cancer", "gene", "motif_id", "n", "pearson_r", "pearson_p", "pearson_fdr")
miss_st <- setdiff(need_st, names(st))
if (length(miss_st)) stop("Missing columns in stats_file: ", paste(miss_st, collapse = ", "))

st[, `:=`(
  cancer      = as.character(cancer),
  gene        = as.character(gene),
  motif_id    = as.character(motif_id),
  n           = as.integer(n),
  pearson_r   = suppressWarnings(as.numeric(pearson_r)),
  pearson_p   = suppressWarnings(as.numeric(pearson_p)),
  pearson_fdr = suppressWarnings(as.numeric(pearson_fdr))
)]
st[, cancer := sub("^TCGA-", "", cancer)]

if (!("tested" %in% names(st))) {
  st[, tested := !is.na(pearson_p) & !is.na(n) & n >= min_n_cor]
}

cat("[3/7] Loading cancer colors...\n")
col_dt <- fread(color_file, header = FALSE)
if (ncol(col_dt) < 2) stop("color_file must have at least 2 columns: cancer and color")
setnames(col_dt, c("cancer", "color"))
col_dt <- unique(col_dt[, .(
  cancer = sub("^TCGA-", "", as.character(cancer)),
  color  = as.character(color)
)])
col_dt <- col_dt[!is.na(cancer) & nzchar(cancer) & !is.na(color) & nzchar(color)]
color_map_global <- setNames(col_dt$color, col_dt$cancer)

existing_keys <- character(0)

cat("[4/7] Selecting top unique pairs...\n")
top_pairs <- st[
  tested == TRUE & !is.na(pearson_fdr)
][order(gene, motif_id, pearson_fdr, -abs(pearson_r), -n)][
  , .SD[1], by = .(gene, motif_id)
][order(pearson_fdr, -abs(pearson_r), -n)]

top_pairs <- top_pairs[, .(
  gene,
  motif_id,
  source_cancer = cancer,
  source_n = n,
  source_r = pearson_r,
  source_p = pearson_p,
  source_padj = pearson_fdr
)]

top_pairs <- top_pairs[seq_len(min(n_top, .N))]
if (nrow(top_pairs) == 0) stop("No valid tested pairs found in stats_file.")
top_pairs[, rank_idx := seq_len(.N)]

fwrite(top_pairs, top_index_file, sep = "\t")

cat("    selected unique pairs:", nrow(top_pairs), "\n")
cat("    saved index:", top_index_file, "\n")

process_one_pair <- function(sel_row) {
  pdf_open <- FALSE

  tryCatch({
    target_gene  <- as.character(sel_row$gene[1])
    target_motif <- as.character(sel_row$motif_id[1])
    rank_idx     <- as.integer(sel_row$rank_idx[1])

    pair_key <- paste(target_gene, target_motif, sep = "||")
    if (pair_key %in% existing_keys) {
      return(data.table(
        rank_idx = rank_idx,
        gene = target_gene,
        motif_id = target_motif,
        status = "skipped_existing",
        error_message = NA_character_,
        pdf_file = NA_character_,
        summary_file = NA_character_
      ))
    }

    plot_dt <- dt[gene == target_gene & motif_id == target_motif]
    if (nrow(plot_dt) == 0) {
      return(data.table(
        rank_idx = rank_idx,
        gene = target_gene,
        motif_id = target_motif,
        status = "no_rows",
        error_message = NA_character_,
        pdf_file = NA_character_,
        summary_file = NA_character_
      ))
    }

    cancers_present <- sort(unique(plot_dt$cancer))
    col_sub <- data.table(cancer = cancers_present)
    col_sub[, color := color_map_global[cancer]]
    if (anyNA(col_sub$color)) col_sub[is.na(color), color := "grey50"]
    color_map <- setNames(col_sub$color, col_sub$cancer)

    plot_dt[, cancer := factor(cancer, levels = names(color_map))]
    plot_dt[, sample_type := factor(as.character(sample_type), levels = c("Healthy", "Tumor"))]

    count_dt <- plot_dt[, .(
      T = sum(sample_type == "Tumor", na.rm = TRUE),
      H = sum(sample_type == "Healthy", na.rm = TRUE)
    ), by = cancer]
    count_dt[, legend_lab := paste0(cancer, " (T=", T, ", H=", H, ")")]
    legend_map <- setNames(count_dt$legend_lab, count_dt$cancer)

    pair_stats <- st[gene == target_gene & motif_id == target_motif]
    missing_stat_cancers <- setdiff(as.character(unique(plot_dt$cancer)), as.character(pair_stats$cancer))

    if (length(missing_stat_cancers)) {
      extra <- plot_dt[cancer %in% missing_stat_cancers, {
        z <- safe_cor_test(meth_percent, log2_expr, min_n = min_n_cor)
        .(
          n = z$n,
          pearson_r = z$r,
          pearson_p = z$p,
          pearson_fdr = NA_real_,
          tested = !is.na(z$p) & z$n >= min_n_cor
        )
      }, by = .(cancer, gene, motif_id)]
      pair_stats <- rbindlist(list(pair_stats, extra), fill = TRUE)
    }

    pair_stats <- pair_stats[, .(cancer, gene, motif_id, n, pearson_r, pearson_p, pearson_fdr, tested)]
    pair_stats[, signif := p_to_stars(pearson_fdr)]
    pair_stats[, cor_lab := paste0(
      "n = ", n,
      "\nr = ", ifelse(is.na(pearson_r), "NA", sprintf("%.2f", pearson_r)),
      "\np = ", fmt_p(pearson_p),
      "\nFDR = ", fmt_p(pearson_fdr),
      "\n", signif
    )]

    global_cor <- safe_cor_test(plot_dt$meth_percent, plot_dt$log2_expr, min_n = min_n_cor)
    global_lab <- paste0(
      "gene = ", target_gene,
      "\nmotif = ", target_motif,
      "\nall-sample n = ", global_cor$n,
      "\nPearson r = ", ifelse(is.na(global_cor$r), "NA", sprintf("%.3f", global_cor$r)),
      "\np = ", fmt_p(global_cor$p)
    )

    base_name <- paste0(
      sprintf("%02d", rank_idx), "_",
      safe_name(target_gene), "__",
      safe_name(target_motif)
    )

    out_pdf  <- file.path(out_dir, paste0(base_name, "_report.pdf"))
    out_info <- file.path(out_dir, paste0(base_name, "_summary.tsv"))

    per_cancer_dbg <- plot_dt[, .(
      n_total = .N,
      n_healthy = sum(sample_type == "Healthy", na.rm = TRUE),
      n_tumor = sum(sample_type == "Tumor", na.rm = TRUE),
      meth_unique_healthy = uniqueN(meth_percent[sample_type == "Healthy"]),
      meth_unique_tumor = uniqueN(meth_percent[sample_type == "Tumor"]),
      expr_unique_healthy = uniqueN(log2_expr[sample_type == "Healthy"]),
      expr_unique_tumor = uniqueN(log2_expr[sample_type == "Tumor"])
    ), by = cancer]

    summary_main <- data.table(
      rank_idx = rank_idx,
      gene = target_gene,
      motif_id = target_motif,
      top_source_cancer = as.character(sel_row$source_cancer[1]),
      top_source_n = as.integer(sel_row$source_n[1]),
      top_source_r = as.numeric(sel_row$source_r[1]),
      top_source_p = as.numeric(sel_row$source_p[1]),
      top_source_padj = as.numeric(sel_row$source_padj[1]),
      all_sample_n = global_cor$n,
      all_sample_r = global_cor$r,
      all_sample_p = global_cor$p
    )

    fwrite(summary_main, out_info, sep = "\t")
    fwrite(per_cancer_dbg, sub("_summary.tsv$", "_per_cancer_counts.tsv", out_info), sep = "\t")

    pdf(out_pdf, width = 13, height = 10, onefile = TRUE)
    pdf_open <- TRUE

    plot_dt_p1 <- copy(plot_dt)
    plot_dt_p1[, sample_type := factor(as.character(sample_type), levels = c("Tumor", "Healthy"))]

    p1 <- ggplot(
      plot_dt_p1,
      aes(x = meth_percent, y = log2_expr, color = cancer, shape = sample_type)
    ) +
      geom_point(size = 2.2, alpha = 0.65) +
      geom_smooth(
        data = plot_dt_p1,
        mapping = aes(x = meth_percent, y = log2_expr, group = 1),
        method = "lm", se = FALSE,
        inherit.aes = FALSE,
        color = "black", linewidth = 0.8
      ) +
      annotate(
        "text",
        x = Inf, y = Inf,
        label = global_lab,
        hjust = 1.02, vjust = 1.05, size = 4
      ) +
      scale_x_continuous(
        limits = c(0, 100),
        breaks = c(0, 25, 50, 75, 100),
        expand = c(0.01, 0.01)
      ) +
      scale_color_manual(
        values = color_map,
        breaks = names(color_map),
        labels = names(color_map),
        drop = FALSE,
        name = "Cancer type"
      ) +
      scale_shape_manual(
        values = c(Tumor = 16, Healthy = 17),
        breaks = c("Tumor", "Healthy"),
        labels = c("Tumor", "Healthy"),
        name = "Sample type"
      ) +
      guides(
        shape = guide_legend(order = 1, override.aes = list(size = 3, alpha = 1, color = "black")),
        color = guide_legend(order = 2, override.aes = list(shape = 16, size = 3, alpha = 1))
      ) +
      labs(
        title = paste0("Top pair across all cancers: ", target_gene, " | ", target_motif),
        subtitle = "Shape shows sample type; color shows cancer type",
        x = "Motif methylation (%) [MAX probe beta]",
        y = "log2(expression + 1)"
      ) +
      theme_bw(base_size = 13) +
      theme(
        plot.title = element_text(face = "bold"),
        legend.position = "right"
      )
    print(p1)

    for (ca in levels(plot_dt$cancer)) {
      dca <- plot_dt[cancer == ca]
      if (nrow(dca) == 0) next

      sc <- pair_stats[cancer == ca]
      sc_lab <- if (nrow(sc) == 0) {
        z <- safe_cor_test(dca$meth_percent, dca$log2_expr, min_n = min_n_cor)
        paste0(
          "n = ", z$n,
          "\nr = ", ifelse(is.na(z$r), "NA", sprintf("%.2f", z$r)),
          "\np = ", fmt_p(z$p),
          "\nFDR = NA\nNA"
        )
      } else {
        sc$cor_lab[1]
      }

      dca_p <- copy(dca)
      dca_p[, draw_group := fifelse(sample_type == "Tumor", 1L, 2L)]
      setorder(dca_p, draw_group)

      p_ca <- ggplot() +
        geom_point(
          data = dca_p[sample_type == "Tumor"],
          aes(x = meth_percent, y = log2_expr),
          color = as.character(color_map[as.character(ca)]),
          shape = 16, size = 1.8, alpha = 0.40
        ) +
        geom_point(
          data = dca_p[sample_type == "Healthy"],
          aes(x = meth_percent, y = log2_expr),
          color = as.character(color_map[as.character(ca)]),
          shape = 17, size = 3.5, stroke = 1.0, alpha = 1
        ) +
        geom_smooth(
          data = dca,
          aes(x = meth_percent, y = log2_expr),
          method = "lm", se = FALSE,
          color = "black", linewidth = 0.8
        ) +
        annotate(
          "text",
          x = min(dca$meth_percent, na.rm = TRUE),
          y = max(dca$log2_expr, na.rm = TRUE),
          label = sc_lab,
          hjust = 0, vjust = 1, size = 4
        ) +
        scale_x_continuous(
          limits = c(0, 100),
          breaks = c(0, 25, 50, 75, 100),
          expand = c(0.01, 0.01)
        ) +
        labs(
          title = paste0("Scatter plot - ", ca, " | ", target_gene, " | ", target_motif),
          subtitle = paste0(
            "Tumor n=", sum(dca$sample_type == "Tumor"),
            " ; Healthy n=", sum(dca$sample_type == "Healthy")
          ),
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

    expr_range_all <- diff(range(expr_vals_all, na.rm = TRUE))
    meth_range_all <- diff(range(meth_vals_all, na.rm = TRUE))
    if (!is.finite(expr_range_all) || expr_range_all == 0) expr_range_all <- 1
    if (!is.finite(meth_range_all) || meth_range_all == 0) meth_range_all <- 1

    expr_y_all <- max(expr_vals_all, na.rm = TRUE) + 0.08 * expr_range_all
    meth_y_all <- max(meth_vals_all, na.rm = TRUE) + 0.08 * meth_range_all

    p5 <- ggplot() +
      geom_boxplot(
        data = d_all,
        aes(x = x_num, y = value_plot, group = x_num),
        width = 0.55, outlier.shape = NA, color = "#0b3c5d", fill = NA
      ) +
      geom_jitter(
        data = d_all,
        aes(x = x_num, y = value_plot, color = cancer),
        width = 0.12, height = 0, size = 1.7, alpha = 0.55
      ) +
      stat_summary(
        data = d_all,
        aes(x = x_num, y = value_plot, group = x_num),
        fun = median,
        geom = "point",
        shape = 95, size = 10, color = "black"
      ) +
      annotate("segment", x = 1, xend = 2, y = expr_y_all, yend = expr_y_all, linewidth = 0.7) +
      annotate("text", x = 1.5, y = expr_y_all + 0.03 * expr_range_all, label = expr_star_all, size = 6) +
      annotate("segment", x = 4, xend = 5, y = meth_y_all, yend = meth_y_all, linewidth = 0.7) +
      annotate("text", x = 4.5, y = meth_y_all + 0.03 * meth_range_all, label = meth_star_all, size = 6) +
      annotate("segment", x = 3, xend = 3, y = min(d_all$value_plot, na.rm = TRUE), yend = max(d_all$value_plot, na.rm = TRUE), linewidth = 1.5) +
      annotate("text", x = 0.45, y = mean(range(d_all$value_plot, na.rm = TRUE)), label = "expression", angle = 90, size = 7) +
      annotate("text", x = 5.55, y = mean(range(d_all$value_plot, na.rm = TRUE)), label = "methylation", angle = 270, size = 7) +
      scale_x_continuous(
        breaks = c(1, 2, 4, 5),
        labels = c("H", "C", "H", "C"),
        limits = c(0.3, 5.7)
      ) +
      scale_y_continuous(
        name = "log2(expression + 1)",
        sec.axis = sec_axis(
          trans = ~ (. - y1_min_all) / sf_all + y2_min_all,
          name = "Motif methylation (%) [MAX probe beta]"
        )
      ) +
      scale_color_manual(
        values = color_map,
        breaks = names(color_map),
        labels = legend_map[names(color_map)],
        drop = FALSE
      ) +
      labs(
        title = paste0("Dual-panel summary across all cancers: ", target_gene, " | ", target_motif),
        subtitle = paste0(
          "Expression H vs C: ", expr_star_all, " (p=", fmt_p(expr_test_all$p),
          ", H=", expr_test_all$n_h, ", C=", expr_test_all$n_c, ") ; ",
          "Methylation H vs C: ", meth_star_all, " (p=", fmt_p(meth_test_all$p),
          ", H=", meth_test_all$n_h, ", C=", meth_test_all$n_c, ")"
        ),
        x = NULL,
        color = "Cancer"
      ) +
      theme_bw(base_size = 13) +
      theme(
        plot.title = element_text(face = "bold"),
        axis.title.x = element_blank(),
        legend.position = "right"
      )
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

      expr_range <- diff(range(expr_vals, na.rm = TRUE))
      meth_range <- diff(range(meth_vals, na.rm = TRUE))
      if (!is.finite(expr_range) || expr_range == 0) expr_range <- 1
      if (!is.finite(meth_range) || meth_range == 0) meth_range <- 1

      expr_y <- max(expr_vals, na.rm = TRUE) + 0.08 * expr_range
      meth_y <- max(meth_vals, na.rm = TRUE) + 0.08 * meth_range

      p6_ca <- ggplot() +
        geom_boxplot(
          data = dca,
          aes(x = x_num, y = value_plot, group = x_num),
          width = 0.55, outlier.shape = NA,
          color = "#0b3c5d", fill = NA
        ) +
        geom_jitter(
          data = dca,
          aes(x = x_num, y = value_plot),
          color = as.character(color_map[as.character(ca)]),
          width = 0.12, height = 0, size = 1.8, alpha = 0.60
        ) +
        stat_summary(
          data = dca,
          aes(x = x_num, y = value_plot, group = x_num),
          fun = median,
          geom = "point",
          shape = 95, size = 10, color = "black"
        ) +
        annotate("segment", x = 1, xend = 2, y = expr_y, yend = expr_y, linewidth = 0.7) +
        annotate("text", x = 1.5, y = expr_y + 0.03 * expr_range, label = expr_star, size = 6) +
        annotate("segment", x = 4, xend = 5, y = meth_y, yend = meth_y, linewidth = 0.7) +
        annotate("text", x = 4.5, y = meth_y + 0.03 * meth_range, label = meth_star, size = 6) +
        annotate("segment", x = 3, xend = 3, y = min(dca$value_plot, na.rm = TRUE), yend = max(dca$value_plot, na.rm = TRUE), linewidth = 1.5) +
        annotate("text", x = 0.45, y = mean(range(dca$value_plot, na.rm = TRUE)), label = "expression", angle = 90, size = 7) +
        annotate("text", x = 5.55, y = mean(range(dca$value_plot, na.rm = TRUE)), label = "methylation", angle = 270, size = 7) +
        scale_x_continuous(
          breaks = c(1, 2, 4, 5),
          labels = c("H", "C", "H", "C"),
          limits = c(0.3, 5.7)
        ) +
        scale_y_continuous(
          name = "log2(expression + 1)",
          sec.axis = sec_axis(
            trans = ~ (. - y1_min) / sf + y2_min,
            name = "Motif methylation (%) [MAX probe beta]"
          )
        ) +
        labs(
          title = paste0("Dual-panel summary - ", ca, " | ", target_gene, " | ", target_motif),
          subtitle = paste0(
            "Expression H vs C: ", expr_star, " (p=", fmt_p(expr_test$p),
            ", H=", expr_test$n_h, ", C=", expr_test$n_c, ") ; ",
            "Methylation H vs C: ", meth_star, " (p=", fmt_p(meth_test$p),
            ", H=", meth_test$n_h, ", C=", meth_test$n_c, ")"
          ),
          x = NULL,
          y = "log2(expression + 1)"
        ) +
        theme_bw(base_size = 13) +
        theme(
          plot.title = element_text(face = "bold"),
          axis.title.x = element_blank()
        )
      print(p6_ca)
    }

    dev.off()
    pdf_open <- FALSE

    data.table(
      rank_idx = rank_idx,
      gene = target_gene,
      motif_id = target_motif,
      status = "ok",
      error_message = NA_character_,
      pdf_file = out_pdf,
      summary_file = out_info
    )

  }, error = function(e) {
    if (pdf_open) try(dev.off(), silent = TRUE)
    data.table(
      rank_idx = as.integer(sel_row$rank_idx[1]),
      gene = as.character(sel_row$gene[1]),
      motif_id = as.character(sel_row$motif_id[1]),
      status = "error",
      error_message = conditionMessage(e),
      pdf_file = NA_character_,
      summary_file = NA_character_
    )
  })
}

cat("[5/7] Running top pairs in parallel batches...\n")
idx_all <- seq_len(nrow(top_pairs))
batch_id <- ceiling(idx_all / batch_size)
batch_split <- split(idx_all, batch_id)
n_batches <- length(batch_split)

res_list <- vector("list", n_batches)

for (b in seq_along(batch_split)) {
  tb <- proc.time()[3]
  idx <- batch_split[[b]]
  batch_dt <- top_pairs[idx]

  cat("  batch", b, "/", n_batches, "| pairs:", nrow(batch_dt), "\n")

  batch_rows <- lapply(seq_len(nrow(batch_dt)), function(i) batch_dt[i])

  batch_res <- mclapply(
    batch_rows,
    process_one_pair,
    mc.cores = min(mc_cores, length(batch_rows))
  )

  res_list[[b]] <- rbindlist(batch_res, fill = TRUE)

  elapsed_batch <- proc.time()[3] - tb
  elapsed_total <- proc.time()[3] - t0
  avg_per_batch <- elapsed_total / b
  eta <- (n_batches - b) * avg_per_batch

  cat(
    "      done | batch time:", fmt_sec(elapsed_batch),
    "| elapsed:", fmt_sec(elapsed_total),
    "| ETA:", fmt_sec(eta), "\n"
  )

  if (any(res_list[[b]]$status == "error")) {
    cat("      errors in batch:\n")
    print(res_list[[b]][status == "error", .(rank_idx, gene, motif_id, error_message)])
  }
}

res <- rbindlist(res_list, fill = TRUE)

cat("[6/7] Writing output index...\n")
fwrite(res, out_index, sep = "\t")

elapsed <- proc.time()[3] - t0
cat("[7/7] Done.\n")
cat("Top-pair index:", top_index_file, "\n")
cat("Generated report index:", out_index, "\n")
cat("Elapsed:", fmt_sec(elapsed), "\n")
EOF
# create the output directory and the file with the selected motif-gene pairs for the local probe controls
mkdir -p ./results/methylation/local_probe_controls

cat > ./results/methylation/local_probe_controls/selected_motif_gene_pairs.tsv << 'EOF'
gene	cancer	TF	motif_id	chr	start	end	strand
GGT6	KIRC	BANP	chr17:4560723:4560733:+	chr17	4560723	4560733	+
TDRD1	PRAD	NRF1	chr10:114179257:114179267:+	chr10	114179257	114179267	+
EOF

# create the output directory for the nearby probe results in a selected window around the motif-gene pairs
mkdir -p ./results/methylation/local_probe_controls

Rscript - << 'EOF'
suppressPackageStartupMessages({
  library(data.table)
})

selected_file <- "./results/methylation/local_probe_controls/selected_motif_gene_pairs.tsv"
probe_file <- "./methylation/annotated_methylation_data_probes_filtered.bed"
out_file <- "./results/methylation/local_probe_controls/nearby_HM450_probes_within_1kb.tsv"

window_bp <- 1000

cat("[INFO] Reading selected motif-gene pairs...\n")
selected <- fread(selected_file)

cat("[INFO] Reading HM450 probe coordinates...\n")
probes <- fread(
  probe_file,
  header = FALSE,
  col.names = c("probe_chr", "probe_start", "probe_end", "probe")
)

selected[, start := as.integer(start)]
selected[, end := as.integer(end)]
probes[, probe_start := as.integer(probe_start)]
probes[, probe_end := as.integer(probe_end)]

cat("[INFO] Searching for probes within +/-", window_bp, "bp around each selected motif...\n")

nearby_list <- lapply(seq_len(nrow(selected)), function(i) {
  
  row <- selected[i]
  
  window_start <- max(0L, row$start - window_bp)
  window_end <- row$end + window_bp
  
  nearby <- probes[
    probe_chr == row$chr &
      probe_start >= window_start &
      probe_end <= window_end
  ]
  
  if (nrow(nearby) == 0) return(NULL)
  
  nearby[, `:=`(
    gene = row$gene,
    cancer = row$cancer,
    TF = row$TF,
    motif_id = row$motif_id,
    motif_chr = row$chr,
    motif_start = row$start,
    motif_end = row$end,
    motif_strand = row$strand,
    window_start = window_start,
    window_end = window_end,
    distance_to_motif = fifelse(
      probe_end < row$start,
      row$start - probe_end,
      fifelse(probe_start > row$end, probe_start - row$end, 0L)
    ),
    overlaps_selected_motif = probe_start < row$end & probe_end > row$start
  )]
  
  nearby[]
})

nearby_probes <- rbindlist(nearby_list, fill = TRUE)

if (nrow(nearby_probes) == 0) {
  stop("No probes found within the selected window. Try increasing window_bp.")
}

setcolorder(
  nearby_probes,
  c(
    "gene", "cancer", "TF", "motif_id",
    "motif_chr", "motif_start", "motif_end", "motif_strand",
    "window_start", "window_end",
    "probe", "probe_chr", "probe_start", "probe_end",
    "distance_to_motif", "overlaps_selected_motif"
  )
)

setorder(nearby_probes, gene, cancer, TF, motif_chr, motif_start, distance_to_motif)

fwrite(nearby_probes, out_file, sep = "\t")

cat("[INFO] Done.\n")
cat("[INFO] Output written to:", out_file, "\n")
cat("[INFO] Total nearby probes found:", nrow(nearby_probes), "\n\n")

cat("[INFO] Probe counts per selected motif:\n")
print(
  nearby_probes[
    ,
    .(
      n_nearby_probes = .N,
      n_overlapping_selected_motif = sum(overlaps_selected_motif),
      n_control_probes = sum(!overlaps_selected_motif),
      closest_control_probe_distance = suppressWarnings(min(distance_to_motif[!overlaps_selected_motif], na.rm = TRUE))
    ),
    by = .(gene, cancer, TF, motif_id)
  ]
)
EOF

# Step 3:
Rscript - << 'EOF'
suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# Paths
# ============================================================

nearby_probe_file <- "./results/methylation/local_probe_controls/nearby_HM450_probes_within_1kb.tsv"
meth_dir          <- "./methylation/filtered_methylation"
expr_file         <- "./expression/gene_expression_matrix_2d_noOV_noCHOL.tsv"
cohort_file       <- "./results/multi_omics/samples_2d_noOV_noCHOL.tsv"

outdir <- "./results/methylation/local_probe_controls/probe_level_1kb"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

out_merged <- file.path(outdir, "nearby_probe_expression_methylation_1kb.tsv.gz")
out_cor    <- file.path(outdir, "nearby_probe_correlation_results_1kb.tsv")

min_n_cor <- 3L

# ============================================================
# Helper functions
# ============================================================

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

# ============================================================
# 1. Read nearby probes
# ============================================================

cat("[INFO] Reading nearby probe table...\n")
nearby <- fread(nearby_probe_file)

nearby[, probe := as.character(probe)]
nearby[, gene := as.character(gene)]
nearby[, cancer := as.character(cancer)]
nearby[, TF := as.character(TF)]
nearby[, motif_id := as.character(motif_id)]

nearby[, probe_class := fifelse(
  overlaps_selected_motif == TRUE,
  "motif_overlapping_probe",
  "nearby_control_probe"
)]

probe_keep <- unique(nearby$probe)
genes_keep <- unique(nearby$gene)
cancers_keep <- unique(nearby$cancer)

cat("[INFO] Probes to extract: ", length(probe_keep), "\n", sep = "")
cat("[INFO] Genes: ", paste(genes_keep, collapse = ", "), "\n", sep = "")
cat("[INFO] Cancers: ", paste(cancers_keep, collapse = ", "), "\n", sep = "")

# ============================================================
# 2. Build methylation sample annotation
# ============================================================

cat("[INFO] Building methylation sample annotation...\n")

cohort <- fread(cohort_file, header = FALSE)
patients_keep <- unique(as.character(cohort[[1]]))
patients_keep <- patients_keep[grepl("^TCGA-[^-]+-[^-]+$", patients_keep)]

files <- list.files(meth_dir, pattern = "\\.bed\\.gz$", full.names = TRUE)

sa <- data.table(file = files, fname = basename(files))
sa[, sample_barcode := extract_sample_barcode_from_meth_file(fname)]
sa[, patient_id := extract_patient_id(sample_barcode)]
sa[, sample_type := extract_sample_type_tcga(sample_barcode)]
sa[, cancer := sub("^HM450_TCGA-([^-]+)-.*$", "\\1", fname)]

sa <- sa[
  !is.na(sample_barcode) &
    !is.na(patient_id) &
    !is.na(sample_type) &
    patient_id %in% patients_keep &
    cancer %in% cancers_keep &
    file.exists(file)
]

sa <- unique(sa[, .(patient_id, sample_barcode, sample_type, cancer, file)])

cat("[INFO] Methylation samples kept: ", nrow(sa), "\n", sep = "")

# ============================================================
# 3. Extract methylation beta values for nearby probes
# ============================================================

cat("[INFO] Extracting probe methylation values from HM450 files...\n")

read_one_sample <- function(i) {
  f  <- sa$file[i]
  sb <- sa$sample_barcode[i]
  st <- sa$sample_type[i]
  ca <- sa$cancer[i]

  dt <- tryCatch(
    fread(f, header = FALSE, select = c(4, 5), showProgress = FALSE),
    error = function(e) NULL
  )

  if (is.null(dt) || nrow(dt) == 0) return(NULL)

  setnames(dt, c("probe", "meth_beta"))
  dt[, probe := as.character(probe)]
  dt[, meth_beta := suppressWarnings(as.numeric(meth_beta))]

  dt <- dt[probe %in% probe_keep & is.finite(meth_beta)]

  if (nrow(dt) == 0) return(NULL)

  dt[, `:=`(
    sample_barcode = sb,
    sample_type = st,
    cancer = ca
  )]

  dt[]
}

meth_list <- lapply(seq_len(nrow(sa)), read_one_sample)
meth_long <- rbindlist(meth_list, fill = TRUE)

if (nrow(meth_long) == 0) stop("No methylation values recovered for nearby probes.")

cat("[INFO] Methylation rows recovered: ", nrow(meth_long), "\n", sep = "")

# Convert to percent when beta is 0-1
if (max(meth_long$meth_beta, na.rm = TRUE) <= 1.5) {
  meth_long[, meth_percent := meth_beta * 100]
} else {
  meth_long[, meth_percent := meth_beta]
}

# ============================================================
# 4. Read expression and keep selected genes/cancers
# ============================================================

cat("[INFO] Reading expression matrix...\n")

expr <- fread(expr_file)

if (!"sample" %in% names(expr)) {
  stop("Expression file must contain a column named 'sample'.")
}

missing_genes <- setdiff(genes_keep, names(expr))
if (length(missing_genes) > 0) {
  stop("Missing genes in expression matrix: ", paste(missing_genes, collapse = ", "))
}

expr <- expr[, c("sample", genes_keep), with = FALSE]

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
expr_long[, sample_type_expr := extract_sample_type_tcga(sample_barcode)]

expr_long <- expr_long[
  cancer %in% cancers_keep &
    gene %in% genes_keep &
    !is.na(sample_barcode) &
    is.finite(expression)
]

# ============================================================
# 5. Merge nearby probe annotation, methylation, and expression
# ============================================================

cat("[INFO] Merging methylation with nearby-probe annotation...\n")

probe_meth <- merge(
  nearby,
  meth_long,
  by = c("probe", "cancer"),
  allow.cartesian = TRUE
)

cat("[INFO] Probe methylation rows after annotation merge: ", nrow(probe_meth), "\n", sep = "")

cat("[INFO] Merging with expression...\n")

merged <- merge(
  probe_meth,
  expr_long[, .(gene, cancer, sample_barcode, expression)],
  by = c("gene", "cancer", "sample_barcode"),
  all = FALSE,
  allow.cartesian = TRUE
)

merged[, log2_expr := log2(expression + 1)]

merged <- merged[
  is.finite(meth_percent) &
    is.finite(log2_expr)
]

if (nrow(merged) == 0) {
  stop("Merged methylation-expression table is empty. Check sample barcode/cancer formats.")
}

setorder(merged, gene, cancer, TF, motif_id, probe_class, distance_to_motif, probe, sample_type, sample_barcode)

fwrite(merged, out_merged, sep = "\t")

cat("[INFO] Merged table written to: ", out_merged, "\n", sep = "")
cat("[INFO] Merged rows: ", nrow(merged), "\n", sep = "")

# ============================================================
# 6. Compute probe-level correlations
# ============================================================

cat("[INFO] Computing probe-level correlations...\n")

cor_dt <- merged[
  ,
  {
    z <- safe_pearson(meth_percent, log2_expr, min_n = min_n_cor)
    .(
      n = z$n,
      pearson_r = z$r,
      pearson_p = z$p,
      n_healthy = sum(sample_type == "Healthy", na.rm = TRUE),
      n_tumor = sum(sample_type == "Tumor", na.rm = TRUE),
      mean_meth_healthy = mean(meth_percent[sample_type == "Healthy"], na.rm = TRUE),
      mean_meth_tumor = mean(meth_percent[sample_type == "Tumor"], na.rm = TRUE),
      mean_expr_healthy = mean(log2_expr[sample_type == "Healthy"], na.rm = TRUE),
      mean_expr_tumor = mean(log2_expr[sample_type == "Tumor"], na.rm = TRUE)
    )
  },
  by = .(
    gene, cancer, TF, motif_id,
    probe, probe_class,
    distance_to_motif,
    probe_chr, probe_start, probe_end,
    overlaps_selected_motif
  )
]

cor_dt[, pearson_fdr := p.adjust(pearson_p, method = "BH"), by = .(gene, cancer, TF, motif_id)]
setorder(cor_dt, gene, cancer, TF, motif_id, probe_class, distance_to_motif, pearson_fdr)

fwrite(cor_dt, out_cor, sep = "\t")

cat("[INFO] Correlation table written to: ", out_cor, "\n\n", sep = "")

cat("[INFO] Probe-level summary:\n")
print(
  cor_dt[
    ,
    .(
      n_probes = .N,
      n_motif_overlapping = sum(probe_class == "motif_overlapping_probe"),
      n_controls = sum(probe_class == "nearby_control_probe"),
      strongest_probe = probe[which.max(abs(pearson_r))],
      strongest_probe_class = probe_class[which.max(abs(pearson_r))],
      strongest_r = pearson_r[which.max(abs(pearson_r))],
      strongest_fdr = pearson_fdr[which.max(abs(pearson_r))]
    ),
    by = .(gene, cancer, TF, motif_id)
  ]
)

cat("\n[INFO] Full correlation results:\n")
print(cor_dt)

EOF

# step 4 
Rscript - << 'EOF'
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ============================================================
# Paths
# ============================================================

merged_file <- "./results/methylation/local_probe_controls/probe_level_1kb/nearby_probe_expression_methylation_1kb.tsv.gz"
cor_file    <- "./results/methylation/local_probe_controls/probe_level_1kb/nearby_probe_correlation_results_1kb.tsv"
color_file  <- "./results/multi_omics/cancer_color_order_with_defined_colours.tsv"

outdir <- "./results/methylation/local_probe_controls/probe_level_1kb/probe_plots"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Helper functions
# ============================================================

safe_name <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

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

clean_sample_type <- function(x) {
  x <- trimws(as.character(x))
  x[tolower(x) %in% c("healthy", "normal", "solid tissue normal", "adjacent normal")] <- "Healthy"
  x[tolower(x) %in% c("tumor", "primary tumor", "primary tumour", "cancer")] <- "Tumor"
  x
}

load_color_map <- function(color_file) {
  col_dt <- fread(color_file, header = FALSE)

  if (ncol(col_dt) < 2) {
    stop("color_file must have at least 2 columns")
  }

  setnames(col_dt, c("cancer", "color"))

  col_dt <- unique(col_dt[, .(
    cancer = sub("^TCGA-", "", as.character(cancer)),
    color  = as.character(color)
  )])

  col_dt <- col_dt[
    !is.na(cancer) & nzchar(cancer) &
      !is.na(color) & nzchar(color)
  ]

  setNames(col_dt$color, col_dt$cancer)
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
    return(list(
      p = NA_real_,
      n_h = n_h,
      n_t = n_t,
      status = "insufficient_n"
    ))
  }

  wt <- tryCatch(
    suppressWarnings(wilcox.test(values ~ groups, exact = FALSE)),
    error = function(e) NULL
  )

  if (is.null(wt)) {
    return(list(
      p = NA_real_,
      n_h = n_h,
      n_t = n_t,
      status = "wilcox_error"
    ))
  }

  list(
    p = wt$p.value,
    n_h = n_h,
    n_t = n_t,
    status = "ok"
  )
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

  if (!is.finite(scale_factor) || scale_factor == 0) {
    scale_factor <- 1
  }

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

# ============================================================
# Read data
# ============================================================

cat("[INFO] Reading merged methylation-expression table...\n")
dt <- fread(merged_file)

cat("[INFO] Reading correlation table...\n")
cor_dt <- fread(cor_file)

cat("[INFO] Reading cancer color map...\n")
color_map_global <- load_color_map(color_file)

# Clean cancer names
dt[, cancer := sub("^TCGA-", "", as.character(cancer))]
cor_dt[, cancer := sub("^TCGA-", "", as.character(cancer))]

# Clean sample type
dt[, sample_type := clean_sample_type(sample_type)]

dt <- dt[
  sample_type %in% c("Healthy", "Tumor") &
    is.finite(meth_percent) &
    is.finite(log2_expr)
]

if (nrow(dt) == 0) {
  stop("No valid rows found in merged table after filtering.")
}

# Label probe class
cor_dt[, probe_class_label := fifelse(
  probe_class == "motif_overlapping_probe",
  "Motif-overlapping probe",
  "Nearby control probe"
)]

# ============================================================
# Build plot key table
# ============================================================

plot_keys <- unique(cor_dt[, .(
  gene,
  cancer,
  TF,
  motif_id,
  probe,
  probe_class,
  probe_class_label,
  distance_to_motif,
  pearson_r,
  pearson_p,
  pearson_fdr,
  n
)])

setorder(plot_keys, gene, cancer, TF, motif_id, probe_class, distance_to_motif, probe)

cat("[INFO] Number of probes to plot: ", nrow(plot_keys), "\n", sep = "")

report_index <- list()

# ============================================================
# Plot one PDF per probe
# ============================================================

for (i in seq_len(nrow(plot_keys))) {

  key <- plot_keys[i]

  d <- dt[
    gene == key$gene &
      cancer == key$cancer &
      TF == key$TF &
      motif_id == key$motif_id &
      probe == key$probe
  ]

  if (nrow(d) == 0) {
    next
  }

  d[, sample_type := factor(sample_type, levels = c("Tumor", "Healthy"))]
  d[, cancer := factor(as.character(cancer), levels = key$cancer)]

  cancer_i <- as.character(key$cancer)

  cancer_color <- color_map_global[cancer_i]

  if (is.na(cancer_color) || !nzchar(cancer_color)) {
    cancer_color <- "grey50"
  }

  probe_label <- ifelse(
    key$probe_class == "motif_overlapping_probe",
    "motif-overlapping probe",
    "nearby control probe"
  )

  stat_label <- paste0(
    "probe = ", key$probe,
    "\nclass = ", probe_label,
    "\ndistance = ", key$distance_to_motif, " bp",
    "\nn = ", key$n,
    "\nr = ", ifelse(is.na(key$pearson_r), "NA", sprintf("%.3f", key$pearson_r)),
    "\np = ", fmt_p(key$pearson_p),
    "\nFDR = ", fmt_p(key$pearson_fdr)
  )

  base_name <- paste0(
    safe_name(key$gene), "__",
    safe_name(key$cancer), "__",
    safe_name(key$TF), "__",
    safe_name(key$probe_class), "__",
    safe_name(key$probe)
  )

  out_pdf <- file.path(outdir, paste0(base_name, "_probe_report.pdf"))

  pdf(out_pdf, width = 10, height = 8)

  # ==========================================================
  # Page 1: Scatter plot
  # ==========================================================

  p1 <- ggplot(
    d,
    aes(
      x = meth_percent,
      y = log2_expr,
      shape = sample_type,
      color = cancer
    )
  ) +
    geom_point(size = 2.3, alpha = 0.65) +
    geom_hline(
      yintercept = 1,
      color = "red",
      linewidth = 0.8,
      linetype = "dashed"
    ) +
    annotate(
      "text",
      x = 95,
      y = 1,
      label = "expression = 1",
      color = "red",
      vjust = -0.5,
      hjust = 1,
      size = 4
    ) +
    geom_smooth(
      aes(group = 1),
      method = "lm",
      se = FALSE,
      color = "black",
      linewidth = 0.8
    ) +
    annotate(
      "text",
      x = Inf,
      y = Inf,
      label = stat_label,
      hjust = 1.02,
      vjust = 1.05,
      size = 4
    ) +
    scale_x_continuous(
      limits = c(0, 100),
      breaks = c(0, 25, 50, 75, 100),
      expand = c(0.01, 0.01)
    ) +
    scale_shape_manual(
      values = c(Tumor = 16, Healthy = 17),
      breaks = c("Tumor", "Healthy"),
      labels = c("Tumor", "Healthy"),
      name = "Sample type"
    ) +
    scale_color_manual(
      values = setNames(cancer_color, cancer_i),
      breaks = cancer_i,
      labels = cancer_i,
      name = "Cancer type",
      drop = FALSE
    ) +
    labs(
      title = paste0(
        key$gene,
        " in ",
        key$cancer,
        " | ",
        key$TF,
        " local CpG control"
      ),
      subtitle = paste0(
        key$probe,
        " | ",
        probe_label,
        " | distance to selected motif: ",
        key$distance_to_motif,
        " bp"
      ),
      x = "Probe methylation (%)",
      y = "log2(expression + 1)"
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )

  print(p1)

  # ==========================================================
  # Page 2: Dual-panel summary
  # ==========================================================

  dual <- make_dual_panel_dt(d)
  d_dual <- dual$dt

  sf <- dual$scale_factor
  y1_min <- dual$y1_min
  y2_min <- dual$y2_min

  expr_test <- safe_wilcox_hc(d$log2_expr, d$sample_type)
  meth_test <- safe_wilcox_hc(d$meth_percent, d$sample_type)

  expr_star <- p_to_stars(expr_test$p)
  meth_star <- p_to_stars(meth_test$p)

  expr_vals <- d_dual[measure == "Expression"]$value_plot
  meth_vals <- d_dual[measure == "Methylation"]$value_plot

  expr_range <- diff(range(expr_vals, na.rm = TRUE))
  if (!is.finite(expr_range) || expr_range == 0) {
    expr_range <- 1
  }

  meth_range <- diff(range(meth_vals, na.rm = TRUE))
  if (!is.finite(meth_range) || meth_range == 0) {
    meth_range <- 1
  }

  expr_y <- max(expr_vals, na.rm = TRUE) + 0.08 * expr_range
  meth_y <- max(meth_vals, na.rm = TRUE) + 0.08 * meth_range

  p2 <- ggplot() +
    geom_boxplot(
      data = d_dual,
      aes(
        x = x_num,
        y = value_plot,
        group = x_num
      ),
      width = 0.55,
      outlier.shape = NA,
      color = "#0b3c5d",
      fill = NA
    ) +
    geom_jitter(
      data = d_dual,
      aes(
        x = x_num,
        y = value_plot
      ),
      color = cancer_color,
      width = 0.12,
      height = 0,
      size = 1.8,
      alpha = 0.60
    ) +
    stat_summary(
      data = d_dual,
      aes(
        x = x_num,
        y = value_plot,
        group = x_num
      ),
      fun = median,
      geom = "point",
      shape = 95,
      size = 10,
      color = "black"
    ) +
    annotate(
      "segment",
      x = 1,
      xend = 2,
      y = expr_y,
      yend = expr_y,
      linewidth = 0.7
    ) +
    annotate(
      "text",
      x = 1.5,
      y = expr_y + 0.03 * expr_range,
      label = expr_star,
      size = 6
    ) +
    annotate(
      "segment",
      x = 4,
      xend = 5,
      y = meth_y,
      yend = meth_y,
      linewidth = 0.7
    ) +
    annotate(
      "text",
      x = 4.5,
      y = meth_y + 0.03 * meth_range,
      label = meth_star,
      size = 6
    ) +
    annotate(
      "segment",
      x = 3,
      xend = 3,
      y = min(d_dual$value_plot, na.rm = TRUE),
      yend = max(d_dual$value_plot, na.rm = TRUE),
      linewidth = 1.5
    ) +
    annotate(
      "text",
      x = 0.45,
      y = mean(range(d_dual$value_plot, na.rm = TRUE)),
      label = "expression",
      angle = 90,
      size = 7
    ) +
    annotate(
      "text",
      x = 5.55,
      y = mean(range(d_dual$value_plot, na.rm = TRUE)),
      label = "methylation",
      angle = 270,
      size = 7
    ) +
    scale_x_continuous(
      breaks = c(1, 2, 4, 5),
      labels = c("H", "C", "H", "C"),
      limits = c(0.3, 5.7)
    ) +
    scale_y_continuous(
      name = "log2(expression + 1)",
      sec.axis = sec_axis(
        trans = ~ (. - y1_min) / sf + y2_min,
        name = "Probe methylation (%)"
      )
    ) +
    labs(
      title = paste0(
        "Dual-panel summary: ",
        key$gene,
        " in ",
        key$cancer
      ),
      subtitle = paste0(
        key$probe,
        " | ",
        probe_label,
        " ; Expression H vs C: ",
        expr_star,
        " (p=",
        fmt_p(expr_test$p),
        ", H=",
        expr_test$n_h,
        ", C=",
        expr_test$n_t,
        ") ; Methylation H vs C: ",
        meth_star,
        " (p=",
        fmt_p(meth_test$p),
        ", H=",
        meth_test$n_h,
        ", C=",
        meth_test$n_t,
        ")"
      ),
      x = NULL
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.title.x = element_blank()
    )

  print(p2)

  dev.off()

  report_index[[length(report_index) + 1]] <- data.table(
    gene = key$gene,
    cancer = key$cancer,
    TF = key$TF,
    motif_id = key$motif_id,
    probe = key$probe,
    probe_class = key$probe_class,
    distance_to_motif = key$distance_to_motif,
    pearson_r = key$pearson_r,
    pearson_p = key$pearson_p,
    pearson_fdr = key$pearson_fdr,
    pdf_file = out_pdf
  )
}

# ============================================================
# Save report index
# ============================================================

index_dt <- rbindlist(report_index, fill = TRUE)

index_file <- file.path(outdir, "generated_probe_reports_index.tsv")
fwrite(index_dt, index_file, sep = "\t")

cat("[INFO] Done.\n")
cat("[INFO] Plot directory: ", outdir, "\n", sep = "")
cat("[INFO] Index file: ", index_file, "\n", sep = "")
EOF
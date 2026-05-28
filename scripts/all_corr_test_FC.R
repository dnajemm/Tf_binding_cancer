# =============================================================================
# Volcano-style scatter plot: Expression log2FC (or methylation delta)
# vs. Pearson correlation (motif methylation ~ gene expression)
#
# For a given transcription factor (TF), this script:
#   1. Reads pre-computed Pearson correlation statistics between TF motif
#      methylation and target gene expression, across cancer types (TCGA).
#   2. Reads a merged methylation-expression matrix and computes, per
#      gene-motif-cancer triplet:
#        - Expression log2FC  (mean Tumor - mean Healthy, in log2 TPM+1)
#        - Methylation delta  (mean Tumor - mean Healthy, in % points)
#   3. Filters to cancer types that have both Healthy and Tumor samples.
#   4. Classifies each pair as "Anti-correlation significant" when:
#        Pearson r < -0.3  AND  FDR < 0.05  AND  |expression log2FC| >= 0
#   5. Loads a CSV of curated gene-cancer pairs (YES/MAYBE significance) and
#      highlights them with red circles and labels on the plot.
#   6. Produces a PDF scatter plot (x = expression log2FC or methylation delta,
#      y = Pearson r) for two sample sets: all samples and matched patients.
#
# Usage:
#   Rscript ./scripts/all_correlation.R <TF_NAME>
#   e.g.: Rscript ./scripts/all_correlation.R BANP
#
# Inputs:
#   - ./results/methylation/correlation_expression_methylation_2d/tumor_vs_healthy/
#       <TF>/{all_samples,matched_patients}/correlation_stats/
#       pearson_correlation_all_pairs_by_cancer.tsv.gz
#   - Same subdirs: motif_sample_expression_methylation_*.tsv.gz
#   - <TF>/<TF>_GENES-Table.csv  (curated gene-cancer pairs to highlight)
#
# Output:
#   - PDF plots saved to: ./results/methylation/correlation_expression_methylation_2d/tumor_vs_healthy/<TF>/volcano_plots/
# =============================================================================

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(grid)
})

args    <- commandArgs(trailingOnly = TRUE)
tf_name <- if (length(args) >= 1) toupper(args[1]) else "BANP"

anti_r_cutoff   <- -0.3
anti_fdr_cutoff <- 0.05
log2fc_cutoff   <- 0

# Choose what you want on the x-axis:
# "expr_log2FC" = expression change Tumor - Healthy
# "meth_delta"  = methylation change Tumor - Healthy
x_axis_to_plot <- "expr_log2FC"

# Number of top significant pairs to label automatically
n_top_labels <- 0

base_dir <- file.path(
  ".",
  "results",
  "methylation",
  "correlation_expression_methylation_2d",
  "tumor_vs_healthy",
  tf_name
)

all_dir     <- file.path(base_dir, "all_samples")
matched_dir <- file.path(base_dir, "matched_patients")

corr_file_all <- file.path(
  all_dir,
  "correlation_stats",
  "pearson_correlation_all_pairs_by_cancer.tsv.gz"
)

corr_file_matched <- file.path(
  matched_dir,
  "correlation_stats",
  "pearson_correlation_all_pairs_by_cancer.tsv.gz"
)

out_dir <- file.path(base_dir, "volcano_plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Extra selected gene-cancer pairs to label manually
# These labels will be written on the plot
# ------------------------------------------------------------

extra_label_pairs <- list(
  list(gene = "TDRD1", cancer = "PRAD")
)

# ------------------------------------------------------------
# Automatically create force_label_pairs from CSV file
# ------------------------------------------------------------

label_file <- file.path(
  base_dir,
  paste0(tf_name, "_GENES-Table.csv")
)

if (!file.exists(label_file)) {
  stop("Label file not found: ", label_file)
}

labels_dt <- fread(label_file, sep = ";", fill = TRUE)

# Clean column names
setnames(labels_dt, trimws(names(labels_dt)))

# Check required columns
required_cols <- c("gene", "Cancer", "Significance")
missing_cols <- setdiff(required_cols, names(labels_dt))

if (length(missing_cols) > 0) {
  stop(
    "Missing required column(s) in label file: ",
    paste(missing_cols, collapse = ", ")
  )
}

# Keep rows with gene and cancer information
labels_dt <- labels_dt[
  !is.na(gene) &
    !is.na(Cancer) &
    trimws(as.character(gene)) != "" &
    trimws(as.character(Cancer)) != ""
]

# Keep only YES and MAYBE rows
labels_dt <- labels_dt[
  toupper(trimws(as.character(Significance))) %in% c("YES", "MAYBE")
]

# Split cancer column when several cancers are written in one cell
labels_long <- labels_dt[
  ,
  .(cancer = unlist(strsplit(as.character(Cancer), ","))),
  by = gene
]

# Clean spaces and TCGA prefixes
labels_long[, gene := trimws(as.character(gene))]
labels_long[, cancer := trimws(as.character(cancer))]
labels_long[, cancer := sub("^TCGA-", "", cancer)]

# Remove empty values and duplicates
labels_long <- unique(labels_long[
  !is.na(gene) &
    !is.na(cancer) &
    gene != "" &
    cancer != ""
])

# Create list format used later by the plot function
if (nrow(labels_long) > 0) {
  force_label_pairs <- lapply(
    seq_len(nrow(labels_long)),
    function(i) {
      list(
        gene = labels_long$gene[i],
        cancer = labels_long$cancer[i]
      )
    }
  )
} else {
  force_label_pairs <- list()
}

cat(">>> Label file:", label_file, "\n")
cat(">>> Number of force_label_pairs loaded from CSV:", length(force_label_pairs), "\n")

if (length(force_label_pairs) > 0) {
  cat(
    ">>> First force_label_pairs:",
    paste(
      head(paste(labels_long$gene, labels_long$cancer, sep = " | "), 10),
      collapse = ", "
    ),
    "\n"
  )
}

# ------------------------------------------------------------
# Plot function:
# x-axis = expression log2FC OR methylation delta
# y-axis = Pearson r
# category = anti-correlation significant or not
# ------------------------------------------------------------

plot_change_vs_r <- function(corr_file, merged_file, label) {

  if (!file.exists(corr_file)) {
    cat(">>> File not found, skipping:", corr_file, "\n")
    return(invisible(NULL))
  }

  if (!file.exists(merged_file)) {
    cat(">>> Merged file not found, skipping:", merged_file, "\n")
    return(invisible(NULL))
  }

  if (!x_axis_to_plot %in% c("expr_log2FC", "meth_delta")) {
    stop("x_axis_to_plot must be either 'expr_log2FC' or 'meth_delta'")
  }

  out_pdf <- file.path(
    out_dir,
    paste0(
      x_axis_to_plot,
      "_vs_pearson_r_anti_correlation_significant_",
      tf_name,
      "_",
      label,
      ".pdf"
    )
  )

  # ------------------------------------------------------------
  # 1. Read correlation file
  # ------------------------------------------------------------

  dt <- fread(corr_file)

  dt[, cancer := sub("^TCGA-", "", as.character(cancer))]
  dt[, gene   := as.character(gene)]

  if (!"motif_id" %in% names(dt)) {
    stop("Column motif_id not found in correlation file: ", corr_file)
  }

  # ------------------------------------------------------------
  # 2. Read merged methylation-expression file
  # ------------------------------------------------------------

  merged_dt <- fread(merged_file)

  merged_dt[, cancer      := sub("^TCGA-", "", as.character(cancer))]
  merged_dt[, gene        := as.character(gene)]
  merged_dt[, motif_id    := as.character(motif_id)]
  merged_dt[, sample_type := trimws(as.character(sample_type))]

  merged_dt[
    tolower(sample_type) %in% c(
      "healthy",
      "normal",
      "solid tissue normal",
      "adjacent normal"
    ),
    sample_type := "Healthy"
  ]

  merged_dt[
    tolower(sample_type) %in% c(
      "tumor",
      "primary tumor",
      "primary tumour",
      "cancer"
    ),
    sample_type := "Tumor"
  ]

  # ------------------------------------------------------------
  # 3. Keep only cancer types with both Healthy and Tumor samples
  # ------------------------------------------------------------

  cancers_both <- merged_dt[
    ,
    .(
      has_healthy = any(sample_type == "Healthy"),
      has_tumor   = any(sample_type == "Tumor")
    ),
    by = cancer
  ][
    has_healthy == TRUE & has_tumor == TRUE,
    cancer
  ]

  n_before <- uniqueN(dt$cancer)
  dt       <- dt[cancer %in% cancers_both]
  n_after  <- uniqueN(dt$cancer)

  cat(">>> [", label, "] Cancers with both sample types:", n_after,
      " | dropped:", n_before - n_after, "\n")

  cat(">>> [", label, "] Cancers kept:",
      paste(sort(cancers_both), collapse = ", "), "\n")

  # ------------------------------------------------------------
  # 4. Prepare expression values
  # ------------------------------------------------------------

  if ("log2_expr" %in% names(merged_dt)) {

    merged_dt[, expr_for_fc := as.numeric(log2_expr)]

  } else if ("expression" %in% names(merged_dt)) {

    merged_dt[, expr_for_fc := log2(as.numeric(expression) + 1)]

  } else {

    stop(
      "Could not find log2_expr or expression column in merged file: ",
      merged_file
    )
  }

  # ------------------------------------------------------------
  # 5. Prepare methylation values
  # ------------------------------------------------------------

  if ("meth_percent" %in% names(merged_dt)) {

    merged_dt[, meth_for_delta := as.numeric(meth_percent)]

  } else if ("meth_beta" %in% names(merged_dt)) {

    merged_dt[, meth_beta := as.numeric(meth_beta)]

    if (max(merged_dt$meth_beta, na.rm = TRUE) <= 1.5) {
      merged_dt[, meth_for_delta := meth_beta * 100]
    } else {
      merged_dt[, meth_for_delta := meth_beta]
    }

  } else {

    stop(
      "Could not find meth_percent or meth_beta column in merged file: ",
      merged_file
    )
  }

  # ------------------------------------------------------------
  # 6. Compute expression log2FC and methylation delta
  # ------------------------------------------------------------

  fc_dt <- merged_dt[
    cancer %in% cancers_both &
      sample_type %in% c("Healthy", "Tumor") &
      !is.na(expr_for_fc) &
      !is.na(meth_for_delta) &
      is.finite(expr_for_fc) &
      is.finite(meth_for_delta),
    .(
      mean_expr_healthy = mean(expr_for_fc[sample_type == "Healthy"], na.rm = TRUE),
      mean_expr_tumor   = mean(expr_for_fc[sample_type == "Tumor"], na.rm = TRUE),

      mean_meth_healthy = mean(meth_for_delta[sample_type == "Healthy"], na.rm = TRUE),
      mean_meth_tumor   = mean(meth_for_delta[sample_type == "Tumor"], na.rm = TRUE),

      n_healthy = sum(
        sample_type == "Healthy" &
          !is.na(expr_for_fc) &
          !is.na(meth_for_delta)
      ),

      n_tumor = sum(
        sample_type == "Tumor" &
          !is.na(expr_for_fc) &
          !is.na(meth_for_delta)
      )
    ),
    by = .(gene, motif_id, cancer)
  ]

  fc_dt <- fc_dt[
    n_healthy > 0 &
      n_tumor > 0 &
      is.finite(mean_expr_healthy) &
      is.finite(mean_expr_tumor) &
      is.finite(mean_meth_healthy) &
      is.finite(mean_meth_tumor)
  ]

  fc_dt[, expr_log2FC := mean_expr_tumor - mean_expr_healthy]
  fc_dt[, meth_delta  := mean_meth_tumor - mean_meth_healthy]

  # ------------------------------------------------------------
  # 7. Merge correlation values with log2FC and methylation delta
  # ------------------------------------------------------------

  dt <- merge(
    dt,
    fc_dt[, .(
      gene,
      motif_id,
      cancer,
      expr_log2FC,
      meth_delta,
      mean_expr_healthy,
      mean_expr_tumor,
      mean_meth_healthy,
      mean_meth_tumor,
      n_healthy,
      n_tumor
    )],
    by = c("gene", "motif_id", "cancer"),
    all = FALSE
  )

  dt <- dt[
    tested == TRUE &
      !is.na(pearson_r) &
      !is.na(pearson_fdr) &
      !is.na(expr_log2FC) &
      !is.na(meth_delta) &
      is.finite(pearson_r) &
      is.finite(pearson_fdr) &
      is.finite(expr_log2FC) &
      is.finite(meth_delta) &
      pearson_fdr > 0
  ]

  if (nrow(dt) == 0) {
    cat(">>> No rows after filtering, skipping:", label, "\n")
    return(invisible(NULL))
  }

  # ------------------------------------------------------------
  # 8. Define categories:
  # Anti-correlation significant or not
  #
  # Anti-correlation significant =
  # Pearson r < -0.3
  # Pearson FDR < 0.05
  # absolute expression log2FC >= 1
  # ------------------------------------------------------------

  dt[, category := fcase(
    pearson_r < anti_r_cutoff &
      pearson_fdr < anti_fdr_cutoff &
      abs(expr_log2FC) >= log2fc_cutoff,
    "Anti-correlation significant",

    default = "Not significant"
  )]

  dt[, category := factor(
    category,
    levels = c("Anti-correlation significant", "Not significant")
  )]

  n_anti <- dt[category == "Anti-correlation significant", .N]
  n_ns   <- dt[category == "Not significant", .N]
  n_cancers <- uniqueN(dt$cancer)

  # This creates a generic plotting column depending on the chosen x-axis
  dt[, plot_x := get(x_axis_to_plot)]

  x_min <- min(dt$plot_x, na.rm = TRUE)
  x_max <- max(dt$plot_x, na.rm = TRUE)

  x_pad <- 0.08 * (x_max - x_min)
  if (!is.finite(x_pad) || x_pad == 0) x_pad <- 0.5

  # ------------------------------------------------------------
  # 9. Automatically label top anti-correlation significant pairs
  # ------------------------------------------------------------

  top_label_dt <- dt[
    category == "Anti-correlation significant"
  ][
    order(pearson_fdr, pearson_r, -abs(expr_log2FC))
  ]

  if (nrow(top_label_dt) > 0 && n_top_labels > 0) {
    top_label_dt <- top_label_dt[1:min(.N, n_top_labels)]
    top_label_dt[, label_text := paste0(gene, " | ", cancer)]
  } else {
    top_label_dt <- dt[0]
    top_label_dt[, label_text := character()]
  }

  cat(">>> [", label, "] Top labelled pairs:", nrow(top_label_dt), "\n")

  # ------------------------------------------------------------
  # 10. Manually/externally selected pairs to circle in red
  # ------------------------------------------------------------

  if (length(force_label_pairs) > 0) {

    force_dt <- rbindlist(lapply(force_label_pairs, as.data.table), fill = TRUE)

    force_dt[, gene   := trimws(as.character(gene))]
    force_dt[, cancer := trimws(as.character(cancer))]
    force_dt[, cancer := sub("^TCGA-", "", cancer)]

    force_dt <- unique(force_dt[
      !is.na(gene) &
        !is.na(cancer) &
        gene != "" &
        cancer != ""
    ])

    circle_dt <- merge(
      dt,
      force_dt,
      by = c("gene", "cancer"),
      all = FALSE
    )

    # Keep only forced pairs that pass all three thresholds
    circle_dt <- circle_dt[
      pearson_r < anti_r_cutoff &
        pearson_fdr < anti_fdr_cutoff &
        abs(expr_log2FC) >= log2fc_cutoff
    ]

  } else {
    circle_dt <- dt[0]
  }

  n_circled <- nrow(circle_dt)

  if (n_circled == 0) {
    cat(">>> WARNING: none of the force_label_pairs passed all thresholds for", label, "\n")
  } else {
    cat(
      ">>> [", label, "] Red-circled pairs:",
      paste(paste(circle_dt$gene, circle_dt$cancer, sep = " | "), collapse = ", "),
      "\n"
    )
  }

  cat(">>> [", label, "] Number of red-circled pairs:", n_circled, "\n")

  # ------------------------------------------------------------
  # 11. Extra manually selected labels
  # ------------------------------------------------------------

  if (length(extra_label_pairs) > 0) {

    extra_label_dt <- rbindlist(lapply(extra_label_pairs, as.data.table), fill = TRUE)

    extra_label_dt[, gene   := trimws(as.character(gene))]
    extra_label_dt[, cancer := trimws(as.character(cancer))]
    extra_label_dt[, cancer := sub("^TCGA-", "", cancer)]

    extra_label_dt <- unique(extra_label_dt[
      !is.na(gene) &
        !is.na(cancer) &
        gene != "" &
        cancer != ""
    ])

    manual_label_dt <- merge(
      dt,
      extra_label_dt,
      by = c("gene", "cancer"),
      all = FALSE
    )

    manual_label_dt[, label_text := paste0(gene, " | ", cancer)]
    manual_label_dt[, plot_x := get(x_axis_to_plot)]

  } else {
    manual_label_dt <- dt[0]
    manual_label_dt[, label_text := character()]
  }

  n_manual_labels <- nrow(manual_label_dt)

  if (n_manual_labels == 0) {
    cat(">>> [", label, "] No extra manual labels found in data\n")
  } else {
    cat(
      ">>> [", label, "] Extra manual labels:",
      paste(manual_label_dt$label_text, collapse = ", "),
      "\n"
    )
  }

  # ------------------------------------------------------------
  # 12. Axis label depending on selected x-axis
  # ------------------------------------------------------------

  if (x_axis_to_plot == "expr_log2FC") {
    x_lab <- "Expression log2FC (Tumor - Healthy; log2(TPM + 1))"
    vline_value <- 0
  } else {
    x_lab <- "Methylation delta (Tumor - Healthy; percentage points)"
    vline_value <- 0
  }

  # ------------------------------------------------------------
  # 13. Plot
  # ------------------------------------------------------------

  pdf(out_pdf, width = 14, height = 9)

  p <- ggplot(dt, aes(x = plot_x, y = pearson_r)) +

    geom_point(
      data  = dt[category == "Not significant"],
      color = "grey75",
      size  = 0.7,
      alpha = 0.35
    ) +

    geom_point(
      data  = dt[category == "Anti-correlation significant"],
      color = "#000000",
      size  = 1.5,
      alpha = 0.90
    ) +

    geom_hline(
      yintercept = anti_r_cutoff,
      linetype   = "dashed",
      color      = "#000000",
      linewidth  = 0.55
    ) +

    geom_vline(
      xintercept = vline_value,
      linetype   = "dashed",
      color      = "grey40",
      linewidth  = 0.55
    ) +

    geom_hline(
      yintercept = 0,
      color      = "grey30",
      linewidth  = 0.35
    ) +

    geom_text_repel(
      data               = top_label_dt,
      aes(x = plot_x, y = pearson_r, label = label_text),
      inherit.aes        = FALSE,
      size               = 3.5,
      color              = "#000000",
      fontface           = "bold",
      max.overlaps       = Inf,
      segment.color      = "grey40",
      segment.size       = 0.4,
      box.padding        = 0.5,
      point.padding      = 0.25,
      min.segment.length = 0.2,
      arrow              = arrow(length = unit(0.01, "npc"))
    ) +

    geom_point(
      data        = circle_dt,
      aes(x = plot_x, y = pearson_r),
      inherit.aes = FALSE,
      color       = "red",
      size        = 5.5,
      shape       = 21,
      fill        = NA,
      stroke      = 1.5
    ) +

    geom_text_repel(
      data               = manual_label_dt,
      aes(x = plot_x, y = pearson_r, label = label_text),
      inherit.aes        = FALSE,
      size               = 3.8,
      color              = "red",
      fontface           = "bold",
      max.overlaps       = Inf,
      segment.color      = "red",
      segment.size       = 0.4,
      box.padding        = 0.6,
      point.padding      = 0.5,
      min.segment.length = 0.1
    ) +

    annotate(
      "text",
      x = x_min,
      y = -0.98,
      label = paste0("Anti-correlation significant: n = ", n_anti),
      hjust = 0,
      size = 4.0,
      color = "#000000",
      fontface = "bold"
    ) +

    annotate(
      "text",
      x = x_max,
      y = 0.95,
      label = paste0("Not significant: n = ", n_ns),
      hjust = 1,
      size = 3.8,
      color = "grey50"
    ) +

    annotate(
      "text",
      x = x_max,
      y = -0.98,
      label = paste0("Red circled: n = ", n_circled),
      hjust = 1,
      size = 4.0,
      color = "red",
      fontface = "bold"
    ) +

    scale_x_continuous(
      limits = c(x_min - x_pad, x_max + x_pad)
    ) +

    scale_y_continuous(
      limits = c(-1, 1),
      breaks = seq(-1, 1, by = 0.25)
    ) +

    labs(
      title = paste0(
        tf_name,
        " - ",
        x_axis_to_plot,
        " versus Pearson r across all gene-motif-cancer pairs"
      ),
      subtitle = paste0(
        label, "  |  ",
        n_cancers, " cancer types with both Healthy and Tumor samples  |  ",
        nrow(dt), " total pairs  |  ",
        "Anti-correlation significant: r < ", anti_r_cutoff,
        ", Pearson FDR < ", anti_fdr_cutoff,
        ", |expression log2FC| >= ", log2fc_cutoff,
        "  |  Top labels: ", n_top_labels,
        "  |  Red circled: ", n_circled
      ),
      x = x_lab,
      y = "Pearson r (motif methylation vs. gene expression)"
    ) +

    theme_bw(base_size = 15) +

    theme(
      plot.title       = element_text(face = "bold"),
      plot.subtitle    = element_text(size = 11, color = "grey40"),
      panel.grid.minor = element_blank()
    )

  print(p)
  dev.off()

  cat(">>> Saved:", out_pdf, "\n")
}

# ------------------------------------------------------------
# Run plots
# ------------------------------------------------------------

plot_change_vs_r(
  corr_file   = corr_file_all,
  merged_file = file.path(
    all_dir,
    "motif_sample_expression_methylation_all.tsv.gz"
  ),
  label = "all_samples"
)

plot_change_vs_r(
  corr_file   = corr_file_matched,
  merged_file = file.path(
    matched_dir,
    "motif_sample_expression_methylation_matched_patients.tsv.gz"
  ),
  label = "matched_patients"
)

# To run:
# Rscript ./scripts/all_correlation.R BANP
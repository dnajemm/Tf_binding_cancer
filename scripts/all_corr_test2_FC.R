# =============================================================================
# Three-plot pipeline:
#
# For a given TF, this script generates 3 PDF plots:
#
#   Plot 1:
#     x-axis = methylation delta
#              mean Tumor methylation - mean Healthy methylation
#     y-axis = -log10 adjusted Wilcoxon p-value for methylation
#
#   Plot 2:
#     x-axis = expression log2FC
#              log2((mean Tumor TPM + pseudocount) /
#                   (mean Healthy TPM + pseudocount))
#     y-axis = -log10 adjusted Wilcoxon p-value for expression
#
#   Plot 3:
#     x-axis = methylation delta
#     y-axis = expression log2FC
#     ONLY points with:
#       methylation Wilcoxon FDR < 0.05
#       expression Wilcoxon FDR < 0.05
#       absolute methylation delta >= min_abs_meth_delta percentage points
#       absolute expression log2FC  >= min_abs_log2fc
#     points colored by Pearson r
#
# Important:
#   - Expression log2FC is classical:
#       log2((mean Tumor TPM + pseudocount) /
#            (mean Healthy TPM + pseudocount))
#
#   - Wilcoxon expression test is performed on log2(TPM+1), consistent with
#     the upstream methylation-expression correlation pipeline.
#
#   - expr_delta_log2 is also saved:
#       mean log2(TPM+1) Tumor - mean log2(TPM+1) Healthy
#     but it is NOT used for plotting.
#
#   - CSV-selected gene-cancer pairs are CIRCLED ONLY in Plot 3.
#   - Pairs in extra_label_pairs can be independently highlighted in Plots 1, 2, and 3.
#   - In Plots 1 and 2, highlighted pairs are circled and labeled in red.
#   - In Plot 3, highlighted pairs are labeled in black; CSV-selected pairs remain circled in black.
#
# Usage:
#   Rscript ./scripts/three_change_plots.R BANP
#   Rscript ./scripts/three_change_plots.R NRF1
# =============================================================================

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(grid)
})

# ------------------------------------------------------------
# Arguments and parameters
# ------------------------------------------------------------

args    <- commandArgs(trailingOnly = TRUE)
tf_name <- if (length(args) >= 1) toupper(args[1]) else "BANP"

wilcox_fdr_cutoff <- 0.05
pearson_r_cutoff  <- -0.3
min_n_per_group   <- 2

# Minimum absolute methylation delta required for Plot 3.
# meth_value is in percentage points, so 5 means 5 percentage points.
# Set to 0 to disable this filter.
min_abs_meth_delta <- 10

# Minimum absolute expression log2FC required for Plot 3.
# 1.0 = at least 2-fold change.
# 0.2 = mild expression change threshold.
# Set to 0 to disable this filter.
min_abs_log2fc <- 1

# Used for classical raw-TPM log2FC:
# log2((mean Tumor TPM + pseudocount) / (mean Healthy TPM + pseudocount))
pseudo_count <- 1

# If TRUE, Plot 3 keeps only anti-correlated pairs.
# Keep FALSE if you want Plot 3 to show all pairs passing Wilcoxon + effect-size filters.
require_anti_correlation_for_plot3 <- FALSE

# Choose independently whether manually selected pairs from extra_label_pairs
# are highlighted in Plot 1, Plot 2, and Plot 3.
#
# For Plot 1 and Plot 2, TRUE adds BOTH:
#   - a red circle around the point
#   - a red text label
#
# For Plot 3, TRUE adds the text label; CSV-selected pairs are circled separately.
label_manual_pairs_plot1 <- TRUE
label_manual_pairs_plot2 <- TRUE
label_manual_pairs_plot3 <- TRUE

# ------------------------------------------------------------
# Manually chosen pairs to HIGHLIGHT
#
# CSV pairs are only circled in Plot 3.
# Pairs written here can be circled and labeled in Plots 1 and 2,
# and labeled in Plot 3, according to the three options above.
# ------------------------------------------------------------

extra_label_pairs <- list(
  list(gene = "GGT6", cancer = "KIRC")
  # Add more manually labeled pairs here:
  # list(gene = "TDRD1", cancer = "PRAD"),
  # list(gene = "GGT6", cancer = "KIRC"),
  # list(gene = "GENE2", cancer = "BRCA")
)

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

out_dir <- file.path(base_dir, "delta_log2fc_wilcoxon_plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------

clean_sample_type <- function(x) {
  x     <- trimws(as.character(x))
  x_low <- tolower(x)
  out   <- x

  out[x_low %in% c(
    "healthy", "normal", "solid tissue normal",
    "adjacent normal", "11", "11a", "11x"
  )] <- "Healthy"

  out[x_low %in% c(
    "tumor", "tumour", "primary tumor", "primary tumour",
    "cancer", "01", "01a", "01x"
  )] <- "Tumor"

  out
}

safe_wilcox <- function(x_healthy, x_tumor, min_n = 2) {
  x_healthy <- as.numeric(x_healthy)
  x_tumor   <- as.numeric(x_tumor)

  x_healthy <- x_healthy[is.finite(x_healthy)]
  x_tumor   <- x_tumor[is.finite(x_tumor)]

  if (length(x_healthy) < min_n || length(x_tumor) < min_n) return(NA_real_)
  if (length(unique(c(x_healthy, x_tumor))) < 2)             return(NA_real_)

  tryCatch(
    suppressWarnings(wilcox.test(x_tumor, x_healthy, exact = FALSE)$p.value),
    error = function(e) NA_real_
  )
}

get_first_existing_col <- function(dt, possible_cols, file_label) {
  found <- intersect(possible_cols, names(dt))
  if (length(found) == 0) {
    stop(
      "None of these columns were found in ", file_label, ": ",
      paste(possible_cols, collapse = ", ")
    )
  }
  found[1]
}

safe_neglog10 <- function(p) {
  p <- as.numeric(p)
  p[p <= 0] <- .Machine$double.xmin
  -log10(p)
}

# ------------------------------------------------------------
# Load CSV-selected gene-cancer pairs for CIRCLES ONLY, Plot 3
# ------------------------------------------------------------

label_file <- file.path(base_dir, paste0(tf_name, "_GENES-Table.csv"))

if (!file.exists(label_file)) stop("Label file not found: ", label_file)

labels_dt <- fread(label_file, sep = ";", fill = TRUE)
setnames(labels_dt, trimws(names(labels_dt)))

required_cols <- c("gene", "Cancer", "Significance")
missing_cols  <- setdiff(required_cols, names(labels_dt))

if (length(missing_cols) > 0) {
  stop(
    "Missing required column(s) in label file: ",
    paste(missing_cols, collapse = ", ")
  )
}

labels_dt[, gene         := trimws(as.character(gene))]
labels_dt[, Cancer       := trimws(as.character(Cancer))]
labels_dt[, Significance := trimws(as.character(Significance))]

labels_dt <- labels_dt[
  !is.na(gene) & !is.na(Cancer) & gene != "" & Cancer != ""
]

labels_dt <- labels_dt[toupper(Significance) %in% c("YES", "MAYBE")]

circle_pairs <- labels_dt[
  , .(cancer = unlist(strsplit(as.character(Cancer), ","))),
  by = gene
]

circle_pairs[, gene   := trimws(as.character(gene))]
circle_pairs[, cancer := trimws(as.character(cancer))]
circle_pairs[, cancer := sub("^TCGA-", "", cancer)]

circle_pairs <- unique(circle_pairs[
  !is.na(gene) & !is.na(cancer) & gene != "" & cancer != ""
])

cat(">>> Label file:", label_file, "\n")
cat(">>> Number of CSV-selected gene-cancer pairs to circle:", nrow(circle_pairs), "\n")

if (nrow(circle_pairs) > 0) {
  cat(
    ">>> First CSV circled pairs:",
    paste(head(paste(circle_pairs$gene, circle_pairs$cancer, sep = " | "), 20), collapse = ", "),
    "\n"
  )
}

# ------------------------------------------------------------
# Load manually chosen gene-cancer pairs for highlighting
# ------------------------------------------------------------

if (length(extra_label_pairs) > 0) {

  manual_label_pairs <- rbindlist(lapply(extra_label_pairs, as.data.table), fill = TRUE)

  manual_label_pairs[, gene   := trimws(as.character(gene))]
  manual_label_pairs[, cancer := trimws(as.character(cancer))]
  manual_label_pairs[, cancer := sub("^TCGA-", "", cancer)]

  manual_label_pairs <- unique(manual_label_pairs[
    !is.na(gene) & !is.na(cancer) & gene != "" & cancer != ""
  ])

} else {
  manual_label_pairs <- data.table(gene = character(), cancer = character())
}

cat(">>> Number of manually selected gene-cancer pairs to label:", nrow(manual_label_pairs), "\n")

if (nrow(manual_label_pairs) > 0) {
  cat(
    ">>> Manual labels:",
    paste(paste(manual_label_pairs$gene, manual_label_pairs$cancer, sep = " | "), collapse = ", "),
    "\n"
  )
}

# ------------------------------------------------------------
# Main function
# ------------------------------------------------------------

make_three_plots <- function(corr_file, merged_file, label) {

  cat("\n============================================================\n")
  cat(">>> Running:", label, "\n")
  cat("============================================================\n")

  if (!file.exists(corr_file)) {
    cat(">>> Correlation file not found, skipping:", corr_file, "\n")
    return(invisible(NULL))
  }

  if (!file.exists(merged_file)) {
    cat(">>> Merged file not found, skipping:", merged_file, "\n")
    return(invisible(NULL))
  }

  # ------------------------------------------------------------
  # 1. Read Pearson correlation file
  # ------------------------------------------------------------

  corr_dt <- fread(corr_file)

  corr_dt[, cancer   := sub("^TCGA-", "", as.character(cancer))]
  corr_dt[, gene     := as.character(gene)]
  corr_dt[, motif_id := as.character(motif_id)]

  if (!"motif_id" %in% names(corr_dt)) {
    stop("Column motif_id not found in correlation file: ", corr_file)
  }

  pearson_r_col <- get_first_existing_col(
    corr_dt,
    c("pearson_r", "r", "correlation"),
    corr_file
  )

  corr_dt[, pearson_r := as.numeric(get(pearson_r_col))]

  pearson_fdr_cols <- intersect(
    c("pearson_fdr", "padj_BH", "fdr", "FDR", "adj_p", "adjusted_pvalue"),
    names(corr_dt)
  )

  if (length(pearson_fdr_cols) > 0) {
    corr_dt[, pearson_fdr := as.numeric(get(pearson_fdr_cols[1]))]
  } else {
    corr_dt[, pearson_fdr := NA_real_]
  }

  if ("tested" %in% names(corr_dt)) {
    corr_dt <- corr_dt[tested == TRUE]
  }

  corr_dt <- corr_dt[!is.na(pearson_r) & is.finite(pearson_r)]

  # ------------------------------------------------------------
  # 2. Read merged methylation-expression file
  # ------------------------------------------------------------

  merged_dt <- fread(merged_file)

  merged_dt[, cancer      := sub("^TCGA-", "", as.character(cancer))]
  merged_dt[, gene        := as.character(gene)]
  merged_dt[, motif_id    := as.character(motif_id)]
  merged_dt[, sample_type := clean_sample_type(sample_type)]

  # ------------------------------------------------------------
  # 3. Prepare expression values
  #
  # expression = raw TPM-like value from upstream matrix
  # expr_value = log2(TPM+1), used for Wilcoxon
  #
  # Classical log2FC will later be computed from raw mean TPM:
  # log2((mean Tumor TPM + pseudocount) /
  #      (mean Healthy TPM + pseudocount))
  # ------------------------------------------------------------

  if ("expression" %in% names(merged_dt)) {

    merged_dt[, expr_raw := suppressWarnings(as.numeric(expression))]
    merged_dt[, expr_raw := pmax(expr_raw, 0)]
    merged_dt[, expr_value := log2(expr_raw + 1)]

  } else if ("log2_expr" %in% names(merged_dt)) {

    merged_dt[, expr_value := suppressWarnings(as.numeric(log2_expr))]

    # If only log2_expr exists, reconstruct approximate raw TPM.
    # Assumes log2_expr = log2(TPM+1).
    merged_dt[, expr_raw := pmax((2^expr_value) - 1, 0)]

  } else {
    stop(
      "Could not find expression or log2_expr column in merged file: ",
      merged_file
    )
  }

  # ------------------------------------------------------------
  # 4. Prepare methylation value
  #
  # Upstream pipeline stores meth_beta.
  # If meth_beta is 0-1, convert to percentage points.
  # ------------------------------------------------------------

  if ("meth_percent" %in% names(merged_dt)) {

    merged_dt[, meth_value := suppressWarnings(as.numeric(meth_percent))]

  } else if ("meth_beta" %in% names(merged_dt)) {

    merged_dt[, meth_beta := suppressWarnings(as.numeric(meth_beta))]

    if (max(merged_dt$meth_beta, na.rm = TRUE) <= 1.5) {
      merged_dt[, meth_value := meth_beta * 100]
    } else {
      merged_dt[, meth_value := meth_beta]
    }

  } else {
    stop(
      "Could not find meth_percent or meth_beta column in merged file: ",
      merged_file
    )
  }

  # ------------------------------------------------------------
  # 5. Keep cancer types with both Healthy and Tumor samples
  # ------------------------------------------------------------

  cancers_both <- merged_dt[
    ,
    .(
      has_healthy = any(sample_type == "Healthy"),
      has_tumor   = any(sample_type == "Tumor")
    ),
    by = cancer
  ][has_healthy == TRUE & has_tumor == TRUE, cancer]

  merged_dt <- merged_dt[
    cancer %in% cancers_both &
      sample_type %in% c("Healthy", "Tumor")
  ]

  corr_dt <- corr_dt[cancer %in% cancers_both]

  cat(">>> Cancers with both Healthy and Tumor:", length(cancers_both), "\n")
  cat(">>> Cancers kept:", paste(sort(cancers_both), collapse = ", "), "\n")

  # ------------------------------------------------------------
  # 6. Compute per gene-motif-cancer changes and Wilcoxon tests
  # ------------------------------------------------------------

  stat_dt <- merged_dt[
    !is.na(expr_value) &
      !is.na(expr_raw) &
      !is.na(meth_value) &
      is.finite(expr_value) &
      is.finite(expr_raw) &
      is.finite(meth_value),
    .(
      # Log2 expression means, used for reference and Wilcoxon consistency
      mean_log_expr_healthy = mean(expr_value[sample_type == "Healthy"], na.rm = TRUE),
      mean_log_expr_tumor   = mean(expr_value[sample_type == "Tumor"],   na.rm = TRUE),

      # Raw TPM means, used for classical log2FC
      mean_expr_raw_healthy = mean(expr_raw[sample_type == "Healthy"], na.rm = TRUE),
      mean_expr_raw_tumor   = mean(expr_raw[sample_type == "Tumor"],   na.rm = TRUE),

      # Methylation means, in percentage points
      mean_meth_healthy = mean(meth_value[sample_type == "Healthy"], na.rm = TRUE),
      mean_meth_tumor   = mean(meth_value[sample_type == "Tumor"],   na.rm = TRUE),

      n_healthy_expr = sum(sample_type == "Healthy" & is.finite(expr_value)),
      n_tumor_expr   = sum(sample_type == "Tumor"   & is.finite(expr_value)),

      n_healthy_meth = sum(sample_type == "Healthy" & is.finite(meth_value)),
      n_tumor_meth   = sum(sample_type == "Tumor"   & is.finite(meth_value)),

      # Wilcoxon test on log2(TPM+1), consistent with upstream pipeline
      expr_wilcox_p = safe_wilcox(
        expr_value[sample_type == "Healthy"],
        expr_value[sample_type == "Tumor"],
        min_n = min_n_per_group
      ),

      # Wilcoxon test on methylation percentage points
      meth_wilcox_p = safe_wilcox(
        meth_value[sample_type == "Healthy"],
        meth_value[sample_type == "Tumor"],
        min_n = min_n_per_group
      )
    ),
    by = .(gene, motif_id, cancer)
  ]

  stat_dt <- stat_dt[
    is.finite(mean_log_expr_healthy) &
      is.finite(mean_log_expr_tumor) &
      is.finite(mean_expr_raw_healthy) &
      is.finite(mean_expr_raw_tumor) &
      is.finite(mean_meth_healthy) &
      is.finite(mean_meth_tumor)
  ]

  # Classical expression log2 fold-change:
  # log2((mean Tumor TPM + pseudocount) / (mean Healthy TPM + pseudocount))
  # This is the value used for plotting.
  stat_dt[,
    expr_log2FC := log2(
      (mean_expr_raw_tumor   + pseudo_count) /
        (mean_expr_raw_healthy + pseudo_count)
    )
  ]

  # Difference on the log2(TPM+1) scale, saved only for reference.
  # This is closer to the scale used by Pearson/Wilcoxon in the upstream pipeline.
  stat_dt[,
    expr_delta_log2 := mean_log_expr_tumor - mean_log_expr_healthy
  ]

  # Methylation delta in percentage points
  stat_dt[, meth_delta := mean_meth_tumor - mean_meth_healthy]

  # BH correction applied within each cancer independently
  stat_dt[
    ,
    expr_wilcox_fdr := p.adjust(expr_wilcox_p, method = "BH"),
    by = cancer
  ]

  stat_dt[
    ,
    meth_wilcox_fdr := p.adjust(meth_wilcox_p, method = "BH"),
    by = cancer
  ]

  stat_dt[, neglog10_expr_fdr := safe_neglog10(expr_wilcox_fdr)]
  stat_dt[, neglog10_meth_fdr := safe_neglog10(meth_wilcox_fdr)]

  # ------------------------------------------------------------
  # 7. Merge with Pearson correlation results
  # ------------------------------------------------------------

  plot_dt <- merge(
    stat_dt,
    corr_dt[, .(gene, motif_id, cancer, pearson_r, pearson_fdr)],
    by  = c("gene", "motif_id", "cancer"),
    all = FALSE
  )

  plot_dt <- plot_dt[
    !is.na(expr_log2FC) &
      !is.na(expr_delta_log2) &
      !is.na(meth_delta) &
      !is.na(expr_wilcox_fdr) &
      !is.na(meth_wilcox_fdr) &
      !is.na(pearson_r) &
      is.finite(expr_log2FC) &
      is.finite(expr_delta_log2) &
      is.finite(meth_delta) &
      is.finite(expr_wilcox_fdr) &
      is.finite(meth_wilcox_fdr) &
      is.finite(pearson_r)
  ]

  if (nrow(plot_dt) == 0) {
    cat(">>> No rows after merging/filtering, skipping:", label, "\n")
    return(invisible(NULL))
  }

  # ------------------------------------------------------------
  # 8. Define significance and effect-size flags
  # ------------------------------------------------------------

  plot_dt[, meth_significant := meth_wilcox_fdr < wilcox_fdr_cutoff]
  plot_dt[, expr_significant := expr_wilcox_fdr < wilcox_fdr_cutoff]
  plot_dt[, both_significant := meth_significant & expr_significant]

  plot_dt[, anti_correlated := pearson_r < pearson_r_cutoff]

  plot_dt[, large_meth_delta := abs(meth_delta)  >= min_abs_meth_delta]
  plot_dt[, large_log2fc     := abs(expr_log2FC) >= min_abs_log2fc]

  # For Plot 1 and Plot 2 coloring:
  # relevant = statistically significant + large enough effect size.
  plot_dt[, meth_relevant := meth_significant & large_meth_delta]
  plot_dt[, expr_relevant := expr_significant & large_log2fc]

  if (require_anti_correlation_for_plot3) {
    plot_dt[, plot3_keep := both_significant &
                              anti_correlated &
                              large_meth_delta &
                              large_log2fc]
  } else {
    plot_dt[, plot3_keep := both_significant &
                              large_meth_delta &
                              large_log2fc]
  }

  # ------------------------------------------------------------
  # 9A. CSV-selected pairs to CIRCLE ONLY in Plot 3
  # ------------------------------------------------------------

  if (nrow(circle_pairs) > 0) {
    circle_dt <- merge(
      plot_dt,
      circle_pairs,
      by              = c("gene", "cancer"),
      all             = FALSE,
      allow.cartesian = TRUE
    )

    circle_dt <- unique(circle_dt, by = c("gene", "motif_id", "cancer"))

  } else {
    circle_dt <- plot_dt[0]
  }

  # ------------------------------------------------------------
  # 9B. Manually selected pairs to highlight
  # ------------------------------------------------------------

  if (nrow(manual_label_pairs) > 0) {
    label_dt <- merge(
      plot_dt,
      manual_label_pairs,
      by              = c("gene", "cancer"),
      all             = FALSE,
      allow.cartesian = TRUE
    )

    label_dt <- unique(label_dt, by = c("gene", "motif_id", "cancer"))
    label_dt[, label_text := paste0(gene, " | ", cancer)]

  } else {
    label_dt <- plot_dt[0]
    label_dt[, label_text := character()]
  }

  cat(">>> Total plotted rows:", nrow(plot_dt), "\n")
  cat(">>> Methylation significant rows:", plot_dt[meth_significant == TRUE, .N], "\n")
  cat(">>> Expression significant rows:",  plot_dt[expr_significant  == TRUE, .N], "\n")
  cat(">>> Both significant rows:",        plot_dt[both_significant  == TRUE, .N], "\n")
  cat(">>> Anti-correlated rows:",         plot_dt[anti_correlated   == TRUE, .N], "\n")

  cat(
    ">>> Rows with abs methylation delta >= ", min_abs_meth_delta, " pp: ",
    plot_dt[large_meth_delta == TRUE, .N], "\n", sep = ""
  )

  cat(
    ">>> Rows with abs log2FC >= ", min_abs_log2fc, ": ",
    plot_dt[large_log2fc == TRUE, .N], "\n", sep = ""
  )

  cat(">>> Plot 3 final rows:",              plot_dt[plot3_keep == TRUE, .N], "\n")
  cat(">>> CSV-selected rows circled:",      nrow(circle_dt), "\n")
  cat(">>> Manually selected rows labeled:", nrow(label_dt), "\n")

  if (nrow(circle_dt) > 0) {
    cat(
      ">>> CSV circled pairs found:",
      paste(unique(paste(circle_dt$gene, circle_dt$cancer, sep = " | ")), collapse = ", "),
      "\n"
    )
  }

  if (nrow(label_dt) > 0) {
    cat(
      ">>> Manual labels found:",
      paste(unique(label_dt$label_text), collapse = ", "),
      "\n"
    )
  }

  # ============================================================
  # Plot 1: methylation delta vs methylation Wilcoxon FDR
  # ============================================================

  out_pdf_1 <- file.path(
    out_dir,
    paste0(tf_name, "_", label, "_01_methylation_delta_vs_methylation_wilcoxon_FDR.pdf")
  )

  pdf(out_pdf_1, width = 12, height = 8)

  p1 <- ggplot(plot_dt, aes(x = meth_delta, y = neglog10_meth_fdr)) +

    geom_point(
      aes(color = meth_relevant),
      size  = 1.2,
      alpha = 0.65
    ) +

    {
      if (label_manual_pairs_plot1 && nrow(label_dt) > 0) {
        geom_point(
          data        = label_dt,
          aes(
            x = meth_delta,
            y = neglog10_meth_fdr
          ),
          inherit.aes = FALSE,
          shape       = 21,
          fill        = NA,
          color       = "red",
          size        = 5.5,
          stroke      = 1.4
        )
      }
    } +

    {
      if (label_manual_pairs_plot1 && nrow(label_dt) > 0) {
        geom_text_repel(
          data          = label_dt,
          aes(
            x     = meth_delta,
            y     = neglog10_meth_fdr,
            label = label_text
          ),
          inherit.aes   = FALSE,
          size          = 3.7,
          color         = "red",
          fontface      = "bold",
          max.overlaps  = Inf,
          segment.color = "red",
          segment.size  = 0.35,
          box.padding   = 0.6,
          point.padding = 0.5
        )
      }
    } +

    scale_color_manual(
      values = c("TRUE" = "black", "FALSE" = "grey75"),
      labels = c(
        "TRUE"  = paste0(
          "Methylation FDR < ", wilcox_fdr_cutoff,
          " and |delta| >= ", min_abs_meth_delta, " pp"
        ),
        "FALSE" = "Not retained"
      ),
      name = NULL
    ) +

    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey40",
      linewidth = 0.5
    ) +

    geom_vline(
      xintercept = c(-min_abs_meth_delta, min_abs_meth_delta),
      linetype = "dotted",
      color = "black",
      linewidth = 0.5
    ) +

    geom_hline(
      yintercept = -log10(wilcox_fdr_cutoff),
      linetype = "dashed",
      color = "black",
      linewidth = 0.5
    ) +

    labs(
      title = paste0(tf_name, " - methylation delta versus methylation Wilcoxon FDR"),
      subtitle = paste0(
        label,
        " | x = mean Tumor methylation - mean Healthy methylation",
        " | y = -log10(BH-adjusted Wilcoxon p-value, per cancer)",
        " | dotted black lines = +/- ", min_abs_meth_delta, " percentage points"
      ),
      x = "Methylation delta: Tumor - Healthy percentage points",
      y = expression(-log[10]("methylation Wilcoxon FDR"))
    ) +

    theme_bw(base_size = 14) +
    theme(
      plot.title       = element_text(face = "bold"),
      plot.subtitle    = element_text(size = 10, color = "grey40"),
      panel.grid.minor = element_blank(),
      legend.position  = "right"
    )

  print(p1)
  dev.off()

  cat(">>> Saved:", out_pdf_1, "\n")

  # ============================================================
  # Plot 2: expression log2FC vs expression Wilcoxon FDR
  # ============================================================

  out_pdf_2 <- file.path(
    out_dir,
    paste0(tf_name, "_", label, "_02_expression_log2FC_vs_expression_wilcoxon_FDR.pdf")
  )

  pdf(out_pdf_2, width = 12, height = 8)

  p2 <- ggplot(plot_dt, aes(x = expr_log2FC, y = neglog10_expr_fdr)) +

    geom_point(
      aes(color = expr_relevant),
      size  = 1.2,
      alpha = 0.65
    ) +

    {
      if (label_manual_pairs_plot2 && nrow(label_dt) > 0) {
        geom_point(
          data        = label_dt,
          aes(
            x = expr_log2FC,
            y = neglog10_expr_fdr
          ),
          inherit.aes = FALSE,
          shape       = 21,
          fill        = NA,
          color       = "red",
          size        = 5.5,
          stroke      = 1.4
        )
      }
    } +

    {
      if (label_manual_pairs_plot2 && nrow(label_dt) > 0) {
        geom_text_repel(
          data          = label_dt,
          aes(
            x     = expr_log2FC,
            y     = neglog10_expr_fdr,
            label = label_text
          ),
          inherit.aes   = FALSE,
          size          = 3.7,
          color         = "red",
          fontface      = "bold",
          max.overlaps  = Inf,
          segment.color = "red",
          segment.size  = 0.35,
          box.padding   = 0.6,
          point.padding = 0.5
        )
      }
    } +

    scale_color_manual(
      values = c("TRUE" = "black", "FALSE" = "grey75"),
      labels = c(
        "TRUE"  = paste0(
          "Expression FDR < ", wilcox_fdr_cutoff,
          " and |log2FC| >= ", min_abs_log2fc
        ),
        "FALSE" = "Not retained"
      ),
      name = NULL
    ) +

    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey40",
      linewidth = 0.5
    ) +

    geom_vline(
      xintercept = c(-min_abs_log2fc, min_abs_log2fc),
      linetype = "dotted",
      color = "black",
      linewidth = 0.5
    ) +

    geom_hline(
      yintercept = -log10(wilcox_fdr_cutoff),
      linetype = "dashed",
      color = "black",
      linewidth = 0.5
    ) +

    labs(
      title = paste0(tf_name, " - expression log2FC versus expression Wilcoxon FDR"),
      subtitle = paste0(
        label,
        " | x = log2((mean Tumor TPM + ", pseudo_count,
        ") / (mean Healthy TPM + ", pseudo_count, "))",
        " | y = -log10(BH-adjusted Wilcoxon p-value, per cancer)",
        " | dotted black lines = +/- ", min_abs_log2fc, " log2FC"
      ),
      x = paste0(
        "Expression log2FC: log2((mean Tumor TPM + ", pseudo_count,
        ") / (mean Healthy TPM + ", pseudo_count, "))"
      ),
      y = expression(-log[10]("expression Wilcoxon FDR"))
    ) +

    theme_bw(base_size = 14) +
    theme(
      plot.title       = element_text(face = "bold"),
      plot.subtitle    = element_text(size = 10, color = "grey40"),
      panel.grid.minor = element_blank(),
      legend.position  = "right"
    )

  print(p2)
  dev.off()

  cat(">>> Saved:", out_pdf_2, "\n")

  # ============================================================
  # Plot 3: methylation delta vs expression log2FC
  # ============================================================

  both_dt         <- plot_dt[plot3_keep == TRUE]
  circle_both_dt <- circle_dt[plot3_keep == TRUE]
  label_both_dt  <- label_dt[plot3_keep == TRUE]

  out_pdf_3 <- file.path(
    out_dir,
    paste0(
      tf_name, "_", label,
      "_03_methylation_delta_vs_expression_log2FC_ONLY_both_wilcoxon_significant_delta",
      min_abs_meth_delta,
      "_expr",
      min_abs_log2fc,
      "_colored_by_pearson_r.pdf"
    )
  )

  pdf(out_pdf_3, width = 14, height = 6)

  if (nrow(both_dt) == 0) {

    plot.new()

    extra_msg <- if (require_anti_correlation_for_plot3) {
      paste0(", Pearson r < ", pearson_r_cutoff)
    } else {
      ""
    }

    text(
      x = 0.5,
      y = 0.5,
      labels = paste0(
        "No pairs have methylation FDR < ", wilcox_fdr_cutoff,
        ", expression FDR < ", wilcox_fdr_cutoff,
        ", abs(delta methylation) >= ", min_abs_meth_delta, " pp",
        ", abs(log2FC) >= ", min_abs_log2fc,
        extra_msg
      ),
      cex = 1.2
    )

  } else {

    p3 <- ggplot(both_dt, aes(x = meth_delta, y = expr_log2FC)) +

      geom_point(
        aes(fill = pearson_r),
        shape  = 21,
        color  = "black",
        stroke = 0.2,
        size   = 2.0,
        alpha  = 0.95
      ) +

      scale_fill_gradient2(
        low      = "blue",
        mid      = "white",
        high     = "red",
        midpoint = 0,
        limits   = c(-1, 1),
        name     = "Pearson r"
      ) +

      geom_vline(
        xintercept = 0,
        linetype = "dashed",
        color = "grey40",
        linewidth = 0.5
      ) +

      geom_vline(
        xintercept = c(-min_abs_meth_delta, min_abs_meth_delta),
        linetype = "dotted",
        color = "black",
        linewidth = 0.5
      ) +

      geom_hline(
        yintercept = 0,
        linetype = "dashed",
        color = "grey40",
        linewidth = 0.5
      ) +

      geom_hline(
        yintercept = c(-min_abs_log2fc, min_abs_log2fc),
        linetype = "dotted",
        color = "black",
        linewidth = 0.5
      ) +

      geom_point(
        data        = circle_both_dt,
        aes(x = meth_delta, y = expr_log2FC),
        inherit.aes = FALSE,
        color       = "black",
        size        = 5.5,
        shape       = 21,
        fill        = NA,
        stroke      = 1.4
      ) +

      {
        if (label_manual_pairs_plot3 && nrow(label_both_dt) > 0) {
          geom_text_repel(
            data          = label_both_dt,
            aes(x = meth_delta, y = expr_log2FC, label = label_text),
            inherit.aes   = FALSE,
            size          = 3.7,
            color         = "black",
            fontface      = "bold",
            max.overlaps  = Inf,
            segment.color = "black",
            segment.size  = 0.35,
            box.padding   = 0.6,
            point.padding = 0.5
          )
        }
      } +

      labs(
        title = paste0(tf_name, " - methylation delta versus expression log2FC"),
        subtitle = paste0(
          label,
          " | ONLY pairs with methylation FDR < ", wilcox_fdr_cutoff,
          ", expression FDR < ", wilcox_fdr_cutoff,
          ", abs(delta methylation) >= ", min_abs_meth_delta, " pp",
          ", abs(log2FC) >= ", min_abs_log2fc,
          if (require_anti_correlation_for_plot3) {
            paste0(", Pearson r < ", pearson_r_cutoff)
          } else {
            ""
          },
          " | color = Pearson r",
          " | CSV circled: n = ", nrow(circle_both_dt),
          " | manually labeled: n = ", nrow(label_both_dt)
        ),
        x = "Methylation delta: Tumor - Healthy percentage points",
        y = paste0(
          "Expression log2FC: log2((mean Tumor TPM + ", pseudo_count,
          ") / (mean Healthy TPM + ", pseudo_count, "))"
        )
      ) +

      theme_bw(base_size = 14) +
      theme(
        plot.title       = element_text(face = "bold"),
        plot.subtitle    = element_text(size = 10, color = "grey40"),
        panel.grid.minor = element_blank(),
        legend.position  = "right"
      )

    print(p3)
  }

  dev.off()

  cat(">>> Saved:", out_pdf_3, "\n")

  # ------------------------------------------------------------
  # Save the table used for plotting
  # ------------------------------------------------------------

  out_table <- file.path(
    out_dir,
    paste0(tf_name, "_", label, "_delta_log2fc_wilcoxon_pearson_plot_table.tsv.gz")
  )

  fwrite(plot_dt, out_table, sep = "\t")

  cat(">>> Saved plot table:", out_table, "\n")
}

# ------------------------------------------------------------
# Run for all samples
# ------------------------------------------------------------

make_three_plots(
  corr_file   = corr_file_all,
  merged_file = file.path(all_dir, "motif_sample_expression_methylation_all.tsv.gz"),
  label       = "all_samples"
)

# ------------------------------------------------------------
# Run for matched patients
# ------------------------------------------------------------

make_three_plots(
  corr_file   = corr_file_matched,
  merged_file = file.path(matched_dir, "motif_sample_expression_methylation_matched_patients.tsv.gz"),
  label       = "matched_patients"
)

# To run:
# Rscript ./scripts/all_corr_test2_FC.R BANP
# Rscript ./scripts/all_corr_test2_FC.R NRF1
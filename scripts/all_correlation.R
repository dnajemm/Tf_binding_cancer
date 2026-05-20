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

# Number of top significant anti-correlated pairs to label automatically
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
  list(gene = "GGT6", cancer = "KIRC")
)

# To add several labels:
# extra_label_pairs <- list(
#   list(gene = "GGT6", cancer = "KIRC"),
#   list(gene = "SFRP1", cancer = "BRCA"),
#   list(gene = "ZNF582", cancer = "COAD")
# )

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

# Create list format used later by the volcano function
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
# Volcano plot function
# ------------------------------------------------------------

plot_volcano <- function(corr_file, merged_file, label) {

  if (!file.exists(corr_file)) {
    cat(">>> File not found, skipping:", corr_file, "\n")
    return(invisible(NULL))
  }

  if (!file.exists(merged_file)) {
    cat(">>> Merged file not found, skipping:", merged_file, "\n")
    return(invisible(NULL))
  }

  out_pdf <- file.path(
    out_dir,
    paste0("volcano_all_pairs_", tf_name, "_", label, ".pdf")
  )

  dt <- fread(corr_file)

  dt[, cancer := sub("^TCGA-", "", as.character(cancer))]
  dt[, gene   := as.character(gene)]

  st_dt <- fread(merged_file, select = c("cancer", "sample_type"))

  st_dt[, cancer      := sub("^TCGA-", "", as.character(cancer))]
  st_dt[, sample_type := trimws(as.character(sample_type))]

  st_dt[
    tolower(sample_type) %in% c(
      "healthy",
      "normal",
      "solid tissue normal",
      "adjacent normal"
    ),
    sample_type := "Healthy"
  ]

  st_dt[
    tolower(sample_type) %in% c(
      "tumor",
      "primary tumor",
      "primary tumour",
      "cancer"
    ),
    sample_type := "Tumor"
  ]

  cancers_both <- st_dt[
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

  dt <- dt[
    tested == TRUE &
      !is.na(pearson_r) &
      !is.na(pearson_fdr) &
      is.finite(pearson_r) &
      is.finite(pearson_fdr) &
      pearson_fdr > 0
  ]

  if (nrow(dt) == 0) {
    cat(">>> No rows after filtering, skipping:", label, "\n")
    return(invisible(NULL))
  }

  dt[, neg_log10_fdr := -log10(pearson_fdr)]

  dt[, category := fcase(
    pearson_r < anti_r_cutoff & pearson_fdr < anti_fdr_cutoff,
    "Anti-correlated (sig.)",
    default = "Not significant"
  )]

  dt[, category := factor(
    category,
    levels = c("Anti-correlated (sig.)", "Not significant")
  )]

  n_anti    <- dt[category == "Anti-correlated (sig.)", .N]
  n_ns      <- dt[category == "Not significant", .N]
  n_cancers <- uniqueN(dt$cancer)

  fdr_line <- -log10(anti_fdr_cutoff)
  y_max    <- max(dt$neg_log10_fdr, na.rm = TRUE)

  # ------------------------------------------------------------
  # 1. Automatically label top significant anti-correlated pairs
  # ------------------------------------------------------------

  top_label_dt <- dt[
    category == "Anti-correlated (sig.)"
  ][
    order(pearson_fdr, pearson_r)
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
  # 2. Manually/externally selected pairs to circle in red
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

    # Keep only forced pairs that pass the volcano thresholds
    circle_dt <- circle_dt[
      pearson_r < anti_r_cutoff &
        pearson_fdr < anti_fdr_cutoff
    ]

  } else {
    circle_dt <- dt[0]
  }

  n_circled <- nrow(circle_dt)

  if (n_circled == 0) {
    cat(">>> WARNING: none of the force_label_pairs were found in the data for", label, "\n")
  } else {
    cat(
      ">>> [", label, "] Manually circled pairs:",
      paste(paste(circle_dt$gene, circle_dt$cancer, sep = " | "), collapse = ", "),
      "\n"
    )
  }

  cat(">>> [", label, "] Number of red-circled pairs:", n_circled, "\n")
 
  # ------------------------------------------------------------
  # 3. Extra manually selected labels
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
  # 3. Plot
  # ------------------------------------------------------------

  pdf(out_pdf, width = 14, height = 9)

  p <- ggplot(dt, aes(x = pearson_r, y = neg_log10_fdr)) +

    geom_point(
      data  = dt[category == "Not significant"],
      color = "grey75",
      size  = 0.7,
      alpha = 0.35
    ) +

    geom_point(
      data  = dt[category == "Anti-correlated (sig.)"],
      color = "#000000",
      size  = 1.4,
      alpha = 0.90
    ) +

    geom_hline(
      yintercept = fdr_line,
      linetype   = "dashed",
      color      = "grey40",
      linewidth  = 0.55
    ) +

    geom_vline(
      xintercept = anti_r_cutoff,
      linetype   = "dashed",
      color      = "#000000",
      linewidth  = 0.55
    ) +

    geom_vline(
      xintercept = 0,
      color      = "grey30",
      linewidth  = 0.35
    ) +

    geom_text_repel(
      data               = top_label_dt,
      aes(label          = label_text),
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
      aes(x = pearson_r, y = neg_log10_fdr),
      inherit.aes = FALSE,
      color       = "red",
      size        = 5.5,
      shape       = 21,
      fill        = NA,
      stroke      = 1.5
    ) +

        geom_text_repel(
      data               = manual_label_dt,
      aes(x = pearson_r, y = neg_log10_fdr, label = label_text),
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
      x = -1,
      y = y_max * 0.98,
      label = paste0("Anti-corr. sig.: n = ", n_anti),
      hjust = 0,
      size = 4.2,
      color = "#000000",
      fontface = "bold"
    ) +

    annotate(
      "text",
      x = 0,
      y = y_max * 0.98,
      label = paste0("Not sig.: n = ", n_ns),
      hjust = 0.5,
      size = 3.8,
      color = "grey50"
    ) +

    annotate(
      "text",
      x = 1,
      y = y_max * 0.98,
      label = paste0("Red circled: n = ", n_circled),
      hjust = 1,
      size = 4.0,
      color = "red",
      fontface = "bold"
    ) +

    scale_x_continuous(
      limits = c(-1, 1),
      breaks = seq(-1, 1, by = 0.25)
    ) +

    scale_y_continuous(
      expand = expansion(mult = c(0.02, 0.08))
    ) +

    labs(
      title = paste0(
        tf_name,
        " - Pearson r versus -log10(FDR) across all gene-motif-cancer pairs"
      ),
      subtitle = paste0(
        label, "  |  ",
        n_cancers, " cancer types with both Healthy and Tumor samples  |  ",
        nrow(dt), " total pairs  |  ",
        "r threshold: ", anti_r_cutoff,
        "  |  FDR threshold: ", anti_fdr_cutoff,
        "  |  Top labels: ", n_top_labels,
        "  |  Red circled: ", n_circled
      ),
      x = "Pearson r  (methylation vs. log2(TPM + 1))",
      y = expression(-log[10](FDR))
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
# Run volcano plots
# ------------------------------------------------------------

plot_volcano(
  corr_file   = corr_file_all,
  merged_file = file.path(
    all_dir,
    "motif_sample_expression_methylation_all.tsv.gz"
  ),
  label = "all_samples"
)

plot_volcano(
  corr_file   = corr_file_matched,
  merged_file = file.path(
    matched_dir,
    "motif_sample_expression_methylation_matched_patients.tsv.gz"
  ),
  label = "matched_patients"
)

# to run this script: Rscript ./scripts/all_correlation.R BANP
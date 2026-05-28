#!/usr/bin/env Rscript

# ============================================================
# Plot one gene-cancer pair:
# scatter correlation + expression boxplot + methylation boxplot
#
# Example:
# Rscript ./scripts/correlation_methylation_oneplot.R BANP GGT6 KIRC
#
# Optional motif:
# Rscript ./scripts/correlation_methylation_oneplot.R BANP GGT6 KIRC "chr17:4560723:4560733:+"
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
  library(gtable)
})

# ============================================================
# ARGUMENTS
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

tf_name <- if (length(args) >= 1) toupper(args[1]) else "BANP"
target_gene <- if (length(args) >= 2) args[2] else "GGT6"
target_cancer <- if (length(args) >= 3) toupper(args[3]) else "KIRC"
target_motif <- if (length(args) >= 4) args[4] else NA_character_

# Choose "all_samples" or "matched_patients"
data_mode <- "all_samples"

# ============================================================
# PATHS
# ============================================================

base_dir <- file.path(
  ".",
  "results",
  "methylation",
  "correlation_expression_methylation_2d",
  "tumor_vs_healthy",
  tf_name
)

if (data_mode == "all_samples") {

  data_file <- file.path(
    base_dir,
    "all_samples",
    "motif_sample_expression_methylation_all.tsv.gz"
  )

  out_dir <- file.path(
    base_dir,
    "all_samples",
    "custom_gene_cancer_pair_plots"
  )

} else if (data_mode == "matched_patients") {

  data_file <- file.path(
    base_dir,
    "matched_patients",
    "motif_sample_expression_methylation_matched_patients.tsv.gz"
  )

  out_dir <- file.path(
    base_dir,
    "matched_patients",
    "custom_gene_cancer_pair_plots"
  )

} else {
  stop("Unsupported data_mode: ", data_mode)
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# PARAMETERS
# ============================================================

min_n_cor <- 3L
min_n_wilcox <- 2L

expression_threshold <- 1

# This controls BOTH:
# 1. width of the left expression boxplot
# 2. height of the bottom methylation boxplot
side_box_size <- unit(6, "cm")

# NEW: extra white margin above the title.
# Increase to unit(1.5, "cm") or unit(2, "cm") if needed.
top_margin_size <- unit(1.5, "cm")

# ============================================================
# HELPERS
# ============================================================

clean_sample_type <- function(x) {
  x <- trimws(as.character(x))

  x[tolower(x) %in% c(
    "healthy",
    "normal",
    "solid tissue normal",
    "adjacent normal"
  )] <- "Healthy"

  x[tolower(x) %in% c(
    "tumor",
    "primary tumor",
    "primary tumour",
    "cancer"
  )] <- "Tumor"

  x
}

safe_name <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

fmt_p <- function(x) {
  if (is.na(x)) return("NA")
  if (x < 2.2e-16) return("<2.2e-16")
  formatC(x, format = "e", digits = 2)
}

p_to_stars <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  "ns"
}

p_label <- function(p) {
  paste0("p = ", fmt_p(p), "\n", p_to_stars(p))
}

safe_pearson <- function(x, y, min_n = 3L) {
  ok <- is.finite(x) & is.finite(y)

  x <- x[ok]
  y <- y[ok]

  n <- length(x)

  if (n < min_n) {
    return(list(
      n = n,
      r = NA_real_,
      p = NA_real_
    ))
  }

  if (length(unique(x)) < 2 || length(unique(y)) < 2) {
    return(list(
      n = n,
      r = NA_real_,
      p = NA_real_
    ))
  }

  ct <- tryCatch(
    cor.test(x, y, method = "pearson"),
    error = function(e) NULL
  )

  if (is.null(ct)) {
    return(list(
      n = n,
      r = NA_real_,
      p = NA_real_
    ))
  }

  list(
    n = n,
    r = unname(ct$estimate),
    p = ct$p.value
  )
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
      n_t = n_t
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
      n_t = n_t
    ))
  }

  list(
    p = wt$p.value,
    n_h = n_h,
    n_t = n_t
  )
}

# ============================================================
# LOAD DATA
# ============================================================

if (!file.exists(data_file)) {
  stop("Data file not found: ", data_file)
}

dt <- fread(data_file)

required_cols <- c(
  "gene",
  "motif_id",
  "sample_barcode",
  "sample_type",
  "cancer",
  "meth_beta",
  "expression"
)

missing_cols <- setdiff(required_cols, names(dt))

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns in data file: ",
    paste(missing_cols, collapse = ", ")
  )
}

dt[, `:=`(
  gene = as.character(gene),
  motif_id = as.character(motif_id),
  sample_barcode = as.character(sample_barcode),
  sample_type = clean_sample_type(sample_type),
  cancer = sub("^TCGA-", "", as.character(cancer)),
  meth_beta = suppressWarnings(as.numeric(meth_beta)),
  expression = suppressWarnings(as.numeric(expression))
)]

dt <- dt[
  !is.na(gene) &
    !is.na(motif_id) &
    !is.na(sample_barcode) &
    sample_type %in% c("Healthy", "Tumor") &
    !is.na(cancer) &
    is.finite(meth_beta) &
    is.finite(expression)
]

dt[, log2_expr := log2(expression + 1)]

if (max(dt$meth_beta, na.rm = TRUE) <= 1.5) {
  dt[, meth_percent := meth_beta * 100]
} else {
  dt[, meth_percent := meth_beta]
}

# ============================================================
# FILTER TO GENE + CANCER
# ============================================================

d0 <- dt[
  gene == target_gene &
    cancer == target_cancer
]

if (nrow(d0) == 0) {
  stop(
    "No rows found for gene = ",
    target_gene,
    " and cancer = ",
    target_cancer,
    ". Check spelling and whether this pair exists in the merged file."
  )
}

# ============================================================
# SELECT MOTIF
# ============================================================

if (is.na(target_motif) || target_motif == "") {

  motif_stats <- d0[, {
    z <- safe_pearson(
      x = meth_percent,
      y = log2_expr,
      min_n = min_n_cor
    )

    .(
      n = z$n,
      pearson_r = z$r,
      pearson_p = z$p
    )

  }, by = motif_id]

  motif_stats <- motif_stats[
    !is.na(pearson_r) &
      !is.na(pearson_p)
  ]

  if (nrow(motif_stats) == 0) {
    stop(
      "No motif had enough usable data for Pearson correlation for ",
      target_gene,
      " in ",
      target_cancer
    )
  }

  # Select the most anti-correlated motif.
  setorder(motif_stats, pearson_r, pearson_p)

  target_motif <- motif_stats$motif_id[1]

  cat("No motif provided.\n")
  cat("Selected most anti-correlated motif:\n")
  print(motif_stats[1])

} else {
  cat("Using provided motif:", target_motif, "\n")
}

d <- d0[motif_id == target_motif]

if (nrow(d) == 0) {
  stop(
    "No rows found for gene = ",
    target_gene,
    ", cancer = ",
    target_cancer,
    ", motif = ",
    target_motif
  )
}

d[, sample_type := factor(sample_type, levels = c("Healthy", "Tumor"))]

# ============================================================
# COMMON AXIS LIMITS
# ============================================================

# Same methylation scale for scatterplot and methylation boxplot.
x_lim <- c(0, 100)
x_breaks <- c(0, 25, 50, 75, 100)

# Same expression scale for scatterplot and expression boxplot.
expr_min <- min(d$log2_expr, na.rm = TRUE)
expr_max <- max(d$log2_expr, na.rm = TRUE)

expr_range <- expr_max - expr_min
if (!is.finite(expr_range) || expr_range == 0) expr_range <- 1

expr_ylim <- c(
  max(0, expr_min - 0.08 * expr_range),
  expr_max + 0.18 * expr_range
)

expr_breaks <- pretty(expr_ylim, n = 5)

# ============================================================
# STATS
# ============================================================

cor_res <- safe_pearson(
  x = d$meth_percent,
  y = d$log2_expr,
  min_n = min_n_cor
)

expr_wilcox <- safe_wilcox_hc(
  values = d$log2_expr,
  groups = d$sample_type,
  min_n = min_n_wilcox
)

meth_wilcox <- safe_wilcox_hc(
  values = d$meth_percent,
  groups = d$sample_type,
  min_n = min_n_wilcox
)

expr_p_lab <- p_label(expr_wilcox$p)
meth_p_lab <- p_label(meth_wilcox$p)

n_healthy <- sum(d$sample_type == "Healthy", na.rm = TRUE)
n_tumor <- sum(d$sample_type == "Tumor", na.rm = TRUE)

stat_label <- paste0(
  "n = ", cor_res$n,
  "\nr = ", ifelse(is.na(cor_res$r), "NA", sprintf("%.2f", cor_res$r)),
  "\np = ", fmt_p(cor_res$p)
)

main_title <- paste0(
  "Scatter plot - ",
  target_cancer,
  " | ",
  target_gene,
  " | ",
  target_motif
)

main_subtitle <- paste0(
  "Tumor n = ",
  n_tumor,
  " ; Healthy n = ",
  n_healthy
)

# ============================================================
# OUTPUT FILES
# ============================================================

base_name <- paste0(
  tf_name,
  "__",
  safe_name(target_gene),
  "__",
  safe_name(target_cancer),
  "__",
  safe_name(target_motif)
)

out_pdf <- file.path(out_dir, paste0(base_name, "_scatter_axis_boxplots.pdf"))
out_png <- file.path(out_dir, paste0(base_name, "_scatter_axis_boxplots.png"))
out_stats <- file.path(out_dir, paste0(base_name, "_stats.tsv"))

stats_out <- data.table(
  TF = tf_name,
  gene = target_gene,
  cancer = target_cancer,
  motif_id = target_motif,
  n_total = nrow(d),
  n_healthy = n_healthy,
  n_tumor = n_tumor,
  pearson_n = cor_res$n,
  pearson_r = cor_res$r,
  pearson_p = cor_res$p,
  expr_wilcox_p = expr_wilcox$p,
  expr_wilcox_label = expr_p_lab,
  meth_wilcox_p = meth_wilcox$p,
  meth_wilcox_label = meth_p_lab
)

fwrite(stats_out, out_stats, sep = "\t")

# ============================================================
# SCATTERPLOT
# ============================================================

scatter_plot <- ggplot() +
  geom_point(
    data = d[sample_type == "Tumor"],
    aes(x = meth_percent, y = log2_expr),
    color = "red",
    shape = 16,
    size = 2.2,
    alpha = 0.65
  ) +
  geom_point(
    data = d[sample_type == "Healthy"],
    aes(x = meth_percent, y = log2_expr),
    color = "black",
    fill = "white",
    shape = 24,
    size = 3.4,
    stroke = 1.1,
    alpha = 1
  ) +
  geom_hline(
    yintercept = expression_threshold,
    color = "red",
    linewidth = 0.7,
    linetype = "dashed"
  ) +
  annotate(
    "text",
    x = 98,
    y = expression_threshold,
    label = paste0("expression = ", expression_threshold),
    color = "red",
    hjust = 1,
    vjust = -0.5,
    size = 4
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = stat_label,
    hjust = 1.05,
    vjust = 1.1,
    size = 4
  ) +
  scale_x_continuous(
    limits = x_lim,
    breaks = x_breaks,
    expand = c(0.01, 0.01)
  ) +
  scale_y_continuous(
    limits = expr_ylim,
    breaks = expr_breaks,
    expand = c(0.01, 0.01)
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 15) +
  theme(
    legend.position = "none",
    plot.margin = margin(0, 10, 0, 0),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

# ============================================================
# EXPRESSION BOXPLOT
# ============================================================

expr_boxplot <- ggplot() +
  geom_boxplot(
    data = d,
    aes(x = sample_type, y = log2_expr),
    width = 0.55,
    outlier.shape = NA,
    color = "#0b3c5d",
    fill = NA
  ) +
  geom_jitter(
    data = d[sample_type == "Tumor"],
    aes(x = sample_type, y = log2_expr),
    color = "red",
    width = 0.12,
    height = 0,
    size = 1.8,
    alpha = 0.65
  ) +
  geom_jitter(
    data = d[sample_type == "Healthy"],
    aes(x = sample_type, y = log2_expr),
    color = "black",
    fill = "white",
    shape = 24,
    width = 0.12,
    height = 0,
    size = 2.8,
    stroke = 1.0,
    alpha = 1
  ) +
  annotate(
    "segment",
    x = 1,
    xend = 2,
    y = expr_ylim[2] - 0.045 * diff(expr_ylim),
    yend = expr_ylim[2] - 0.045 * diff(expr_ylim),
    linewidth = 0.6
  ) +
  annotate(
    "text",
    x = 1.5,
    y = expr_ylim[2] - 0.005 * diff(expr_ylim),
    label = expr_p_lab,
    size = 4.3,
    lineheight = 0.9
  ) +
  scale_y_continuous(
    limits = expr_ylim,
    breaks = expr_breaks,
    expand = c(0.01, 0.01)
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 0),
    legend.position = "none",
    plot.margin = margin(0, 5, 0, 0),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

# ============================================================
# METHYLATION BOXPLOT
# ============================================================

meth_boxplot <- ggplot() +
  geom_boxplot(
    data = d,
    aes(x = meth_percent, y = sample_type),
    width = 0.55,
    outlier.shape = NA,
    color = "#0b3c5d",
    fill = NA
  ) +
  geom_jitter(
    data = d[sample_type == "Tumor"],
    aes(x = meth_percent, y = sample_type),
    color = "red",
    width = 0,
    height = 0.12,
    size = 1.8,
    alpha = 0.65
  ) +
  geom_jitter(
    data = d[sample_type == "Healthy"],
    aes(x = meth_percent, y = sample_type),
    color = "black",
    fill = "white",
    shape = 24,
    width = 0,
    height = 0.12,
    size = 2.8,
    stroke = 1.0,
    alpha = 1
  ) +
  annotate(
    "segment",
    x = 88,
    xend = 98,
    y = 1.5,
    yend = 1.5,
    linewidth = 0.6
  ) +
  annotate(
    "text",
    x = 93,
    y = 1.73,
    label = meth_p_lab,
    size = 4.3,
    lineheight = 0.9
  ) +
  scale_x_continuous(
    limits = x_lim,
    breaks = x_breaks,
    expand = c(0.01, 0.01)
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    plot.margin = margin(5, 10, 0, 0),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

# ============================================================
# COMBINE PLOTS
# ============================================================

draw_combined_plot <- function() {
  grid.newpage()

  # Convert plots to grobs.
  expr_grob <- ggplotGrob(expr_boxplot)
  scatter_grob <- ggplotGrob(scatter_plot)
  meth_grob <- ggplotGrob(meth_boxplot)

  # ------------------------------------------------------------
  # Force scatterplot and methylation boxplot to have identical
  # internal widths. This makes x = 0 and x = 100 align visually.
  # ------------------------------------------------------------

  common_x_widths <- unit.pmax(scatter_grob$widths, meth_grob$widths)
  scatter_grob$widths <- common_x_widths
  meth_grob$widths <- common_x_widths

  # ------------------------------------------------------------
  # Force expression boxplot and scatterplot to use the same
  # internal heights. This helps their y-axis scales align visually.
  # ------------------------------------------------------------

  common_y_heights <- unit.pmax(expr_grob$heights, scatter_grob$heights)
  expr_grob$heights <- common_y_heights
  scatter_grob$heights <- common_y_heights

  # ------------------------------------------------------------
  # Page layout
  #
  # col 1 = shared y-axis title
  # col 2 = expression boxplot
  # col 3 = scatterplot + methylation boxplot
  #
  # row 1 = extra top margin
  # row 2 = centered title
  # row 3 = centered subtitle
  # row 4 = expression boxplot + scatterplot
  # row 5 = methylation boxplot
  # row 6 = shared methylation x-axis label
  # ------------------------------------------------------------

  pushViewport(
    viewport(
      layout = grid.layout(
        nrow = 6,
        ncol = 3,

        widths = unit.c(
          unit(0.8, "cm"),   # shared y-axis title column
          side_box_size,     # expression boxplot width
          unit(1, "null")    # scatterplot/methylation width
        ),

        heights = unit.c(
          top_margin_size,   # NEW: extra top margin
          unit(0.65, "cm"),  # title
          unit(0.45, "cm"),  # subtitle
          unit(1, "null"),   # scatterplot/expression height
          side_box_size,     # methylation boxplot height
          unit(0.65, "cm")   # shared methylation x-axis title
        )
      )
    )
  )

  # Main title, centered across the full page.
  grid.text(
    main_title,
    x = unit(0.5, "npc"),
    y = unit(0.5, "npc"),
    just = "center",
    gp = gpar(
      fontsize = 20,
      fontface = "bold"
    ),
    vp = viewport(
      layout.pos.row = 2,
      layout.pos.col = 1:3
    )
  )

  # Subtitle, centered across the full page.
  grid.text(
    main_subtitle,
    x = unit(0.5, "npc"),
    y = unit(0.5, "npc"),
    just = "center",
    gp = gpar(fontsize = 14),
    vp = viewport(
      layout.pos.row = 3,
      layout.pos.col = 1:3
    )
  )

  # Shared expression y-axis title.
  grid.text(
    "Gene expression (log2(TPM + 1))",
    x = unit(0.5, "npc"),
    y = unit(0.5, "npc"),
    rot = 90,
    just = "center",
    gp = gpar(fontsize = 15),
    vp = viewport(
      layout.pos.row = 4,
      layout.pos.col = 1
    )
  )

  # Expression boxplot.
  grid.draw(
    editGrob(
      expr_grob,
      vp = viewport(
        layout.pos.row = 4,
        layout.pos.col = 2
      )
    )
  )

  # Scatterplot.
  grid.draw(
    editGrob(
      scatter_grob,
      vp = viewport(
        layout.pos.row = 4,
        layout.pos.col = 3
      )
    )
  )

  # Methylation boxplot.
  grid.draw(
    editGrob(
      meth_grob,
      vp = viewport(
        layout.pos.row = 5,
        layout.pos.col = 3
      )
    )
  )

  # Shared methylation x-axis title.
  grid.text(
    "Motif methylation (%)",
    x = unit(0.5, "npc"),
    y = unit(0.5, "npc"),
    just = "center",
    gp = gpar(fontsize = 15),
    vp = viewport(
      layout.pos.row = 6,
      layout.pos.col = 3
    )
  )
}

# ============================================================
# SAVE
# ============================================================

pdf(out_pdf, width = 20, height = 16)
draw_combined_plot()
dev.off()

png(out_png, width = 2600, height = 2000, res = 200)
draw_combined_plot()
dev.off()

cat("\nSaved PDF:\n", out_pdf, "\n")
cat("\nSaved PNG:\n", out_png, "\n")
cat("\nSaved stats:\n", out_stats, "\n")
cat("\nDone.\n")
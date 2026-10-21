##################################################################################################################################################################################################################
# Same integration methylation_expression pipeline as before just here this is to get one pdf page per gene-cancer pair with the scatter plot and the boxplots for methylation and expression
# to have the gene with multiple cancer type detected in the same page in the same pdf for the report file
##################################################################################################################################################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
  library(gtable)
  library(gridExtra)
  library(parallel)
})

# ── USER CONFIG ───────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop(
  "Usage: Rscript script.R <TF> <GENE> <CANCER1,CANCER2,...> [PROBE] ",
  "[FINAL_SCORE] [MEAN_PEARSON_R] [MEAN_DELTA_BETA] [MEAN_LOG2FC] ",
  "[BEST_PEARSON_P] [BEST_PADJ] [BEST_ADJ_P_VAL] [BEST_EXPR_WILCOX_P] [BEST_METH_WILCOX_P]")
arg_num <- function(i) if (length(args) >= i && nzchar(args[i])) as.numeric(args[i]) else NA_real_
arg_str <- function(i) if (length(args) >= i && nzchar(args[i])) args[i]          else NA_character_

tf_name          <- args[1]
target_gene      <- args[2]
target_cancers   <- strsplit(args[3], ",")[[1]]
target_probe     <- if (length(args) >= 4 && nzchar(args[4])) args[4] else NULL
final_score      <- arg_num(5)
arg_mean_r       <- arg_num(6)
arg_mean_db      <- arg_num(7)
arg_mean_log2fc  <- arg_num(8)
arg_best_pear_p  <- arg_str(9)
arg_best_padj    <- arg_str(10)
arg_best_adj_p   <- arg_str(11)
arg_best_expr_wp <- arg_str(12)
arg_best_meth_wp <- arg_str(13)

N_CORES  <- max(1L, detectCores() - 2L)

pairs_tsv  <- "./results/integration/meth_expr_filtered_pairs.tsv"
motif_file <- "./methylation/probe_gene_pairs_in_motifs_with_tss.tsv.gz"
beta_dir   <- paste0("./methylation/beta_matrices_", tf_name)
tpm_file   <- "./expression/gene_expression_matrix_protein_coding.tsv"
out_dir    <- file.path("./results/integration/gene_cancer_plots_NRF1/")
stats_out <- paste0("./results/integration/gene_cancer_plots_NRF1/meth_expr_wilcoxon_", tf_name,
                    "_", target_gene, "_grid_stats.tsv")

min_n_cor            <- 3L
min_n_wilcox         <- 2L
expression_threshold <- 1
side_box_size        <- unit(6, "cm")
top_margin_size      <- unit(1.5, "cm")
# ─────────────────────────────────────────────────────────────────────────────

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fmt_p <- function(x) {
  if (is.na(x)) return("NA")
  if (x == 0)   return("<2.2e-16")
  formatC(x, format = "e", digits = 2)
}
p_to_stars <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  "ns"
}
p_label <- function(p) paste0("p = ", fmt_p(p), "\n", p_to_stars(p))

safe_pearson <- function(x, y, min_n = 3L) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]; n <- length(x)
  if (n < min_n || length(unique(x)) < 2 || length(unique(y)) < 2)
    return(list(n = n, r = NA_real_, p = NA_real_))
  ct <- tryCatch(cor.test(x, y, method = "pearson"), error = function(e) NULL)
  if (is.null(ct)) return(list(n = n, r = NA_real_, p = NA_real_))
  list(n = n, r = unname(ct$estimate), p = ct$p.value)
}

safe_wilcox <- function(values, groups, min_n = 2L) {
  keep   <- is.finite(values) & groups %in% c("Healthy", "Tumor")
  values <- values[keep]; groups <- groups[keep]
  n_h <- sum(groups == "Healthy"); n_t <- sum(groups == "Tumor")
  if (n_h < min_n || n_t < min_n) return(list(p = NA_real_, n_h = n_h, n_t = n_t))
  wt <- tryCatch(suppressWarnings(wilcox.test(values ~ groups, exact = FALSE)),
                 error = function(e) NULL)
  list(p = if (is.null(wt)) NA_real_ else wt$p.value, n_h = n_h, n_t = n_t)
}

# ── Resolve probe ─────────────────────────────────────────────────────────────
if (is.null(target_probe)) {
  if (file.exists(pairs_tsv)) {
    pairs <- fread(pairs_tsv)
    pairs <- pairs[TF == tf_name & gene == target_gene & !is.na(pearson_r)]
    if (nrow(pairs) > 0) {
      target_probe <- pairs[order(-abs(pearson_r)), probe][1]
      cat("Auto-selected probe:", target_probe, "\n")
    } else {
      stop("No pairs found for gene ", target_gene, " in ", pairs_tsv,
           ". Set target_probe manually.")
    }
  } else {
    stop("pairs_tsv not found and target_probe is NULL. Set target_probe manually.")
  }
} else {
  cat("Using probe:", target_probe, "\n")
}

# ── Resolve cancers ───────────────────────────────────────────────────────────
if (is.null(target_cancers)) {
  beta_files     <- list.files(beta_dir, pattern = "\\.tsv\\.gz$")
  target_cancers <- sort(sub(paste0("_", tf_name, "_beta_matrix\\.tsv\\.gz$"), "", beta_files))
  target_cancers <- target_cancers[nchar(target_cancers) > 0]
  cat("Auto-detected cancers:", paste(target_cancers, collapse = ", "), "\n")
} else {
  cat("Using cancers:", paste(target_cancers, collapse = ", "), "\n")
}
if (length(target_cancers) == 0) stop("No cancer types found.")

# ── Load TPM matrix ───────────────────────────────────────────────────────────
cat("Loading TPM matrix…\n")
tpm_mat <- fread(tpm_file, header = TRUE, check.names = FALSE)
setnames(tpm_mat, 1L, "tpm_sample_id")
tpm_mat[, patient    := sub("^TCGA-[A-Z]+-TCGA-([A-Z0-9]+-[A-Z0-9]+)-[0-9]+[A-Z]_.*", "\\1", tpm_sample_id)]
tpm_mat[, cancer_col := sub("^TCGA-([A-Z]+)-.*", "\\1", tpm_sample_id)]
tpm_mat[, st         := fifelse(grepl("-01A_1", tpm_sample_id), "Tumor",
                        fifelse(grepl("-11A_1", tpm_sample_id), "Healthy", "other"))]
tpm_mat[, sample_key := paste(patient, st, sep = "_")]
tpm_mat <- tpm_mat[st != "other"]
cat(sprintf("TPM matrix: %d samples\n", nrow(tpm_mat)))
if (!target_gene %in% names(tpm_mat))
  stop("Gene '", target_gene, "' not found in TPM matrix.")

# ── Load data for one cancer ──────────────────────────────────────────────────
load_cancer <- function(cancer) {
  tag <- sprintf("[%s | %s]", cancer, target_gene)
  bf  <- list.files(beta_dir, pattern = paste0(cancer, ".*\\.tsv\\.gz$"), full.names = TRUE)
  if (length(bf) == 0) { cat(tag, "– no beta matrix\n"); return(NULL) }
  beta_mat <- tryCatch(fread(bf[1]), error = function(e) NULL)
  if (is.null(beta_mat)) { cat(tag, "– fread failed\n"); return(NULL) }
  if (!"probe" %in% colnames(beta_mat)) setnames(beta_mat, 1, "probe")

  n_matrix_samples <- ncol(beta_mat) - 1L

  beta_mat <- beta_mat[probe == target_probe]
  if (nrow(beta_mat) == 0) { cat(tag, "– probe not in matrix\n"); return(NULL) }
  beta_long <- melt(beta_mat, id.vars = "probe",
                    variable.name = "col_name", value.name = "meth_beta")

  n_na_probe <- sum(is.na(beta_long$meth_beta))

  beta_long[, patient     := sub("^[A-Z]+-TCGA-([A-Z0-9]+-[A-Z0-9]+)_(tumor|healthy)$",
                                 "\\1", as.character(col_name))]
  beta_long[, sample_type := sub("^.*_(tumor|healthy)$", "\\1", as.character(col_name))]
  beta_long[, sample_type := fifelse(sample_type == "tumor",   "Tumor",
                              fifelse(sample_type == "healthy", "Healthy", sample_type))]
  beta_long[, sample_key  := paste(patient, sample_type, sep = "_")]
  if (max(beta_long$meth_beta, na.rm = TRUE) <= 1.5)
    beta_long[, meth_percent := meth_beta * 100] else
    beta_long[, meth_percent := meth_beta]
  expr_dt <- tpm_mat[cancer_col == cancer, .(sample_key, count = get(target_gene))]
  if (nrow(expr_dt) == 0) { cat(tag, "– no TPM rows\n"); return(NULL) }
  d <- merge(
    beta_long[, .(sample_key, patient, sample_type, meth_percent)],
    expr_dt[,   .(sample_key, count)],
    by = "sample_key"
  )
  d <- d[sample_type %in% c("Tumor","Healthy") & is.finite(meth_percent) & is.finite(count)]
  if (nrow(d) == 0) { cat(tag, "– no matched samples\n"); return(NULL) }
  d[, log2_expr        := log2(count + 1)]
  d[, sample_type      := factor(sample_type, levels = c("Healthy","Tumor"))]
  d[, cancer           := cancer]
  d[, n_matrix_samples := n_matrix_samples]
  d[, n_na_probe       := n_na_probe]
  cat(tag, sprintf("– %d rows (T=%d H=%d) | matrix=%d NA=%d\n",
        nrow(d), sum(d$sample_type=="Tumor"), sum(d$sample_type=="Healthy"),
        n_matrix_samples, n_na_probe))
  d
}

# ── draw_cell() ───────────────────────────────────────────────────────────────
draw_cell <- function(d, cancer, expr_ylim, expr_breaks) {
  x_lim    <- c(0, 100)
  x_breaks <- c(0, 25, 50, 75, 100)
  n_tumor          <- sum(d$sample_type == "Tumor")
  n_healthy        <- sum(d$sample_type == "Healthy")
  n_matrix_samples <- d$n_matrix_samples[1]
  n_na_probe       <- d$n_na_probe[1]
  cor_res     <- safe_pearson(d$meth_percent, d$log2_expr, min_n_cor)
  expr_wilcox <- safe_wilcox(d$log2_expr,    as.character(d$sample_type), min_n_wilcox)
  meth_wilcox <- safe_wilcox(d$meth_percent, as.character(d$sample_type), min_n_wilcox)
  stat_label  <- paste0("n = ", cor_res$n,
                         "\nr = ", ifelse(is.na(cor_res$r), "NA", sprintf("%.2f", cor_res$r)),
                         "\np = ", fmt_p(cor_res$p))
  main_title    <- paste0(cancer, " | ", target_gene, " | ", target_probe)
  main_subtitle <- paste0(
    "Tumor n = ", n_tumor, "  ;  Healthy n = ", n_healthy,
    "     [matrix: ", n_matrix_samples, " samples | NA for probe: ", n_na_probe, "]"
  )

  scatter_plot <- ggplot() +
    geom_point(data = d[sample_type == "Tumor"],
               aes(x = meth_percent, y = log2_expr),
               color = "red", shape = 16, size = 2.2, alpha = 0.65) +
    geom_point(data = d[sample_type == "Healthy"],
               aes(x = meth_percent, y = log2_expr),
               color = "black", fill = "white", shape = 24, size = 3.4, stroke = 1.1, alpha = 1) +
    geom_hline(yintercept = expression_threshold, color = "red",
               linewidth = 0.7, linetype = "dashed") +
    annotate("text", x = 98, y = expression_threshold,
             label = paste0("expression = ", expression_threshold),
             color = "red", hjust = 1, vjust = -0.5, size = 14) +
    annotate("text", x = Inf, y = Inf, label = stat_label,
             hjust = 1.05, vjust = 1.1, size = 14) +
    scale_x_continuous(limits = x_lim, breaks = x_breaks, expand = c(0.01, 0.01)) +
    scale_y_continuous(limits = expr_ylim, breaks = expr_breaks, expand = c(0.01, 0.01)) +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 40) +
    theme(legend.position = "none", plot.margin = margin(0, 10, 0, 0))

  expr_boxplot <- ggplot() +
    geom_boxplot(data = d, aes(x = sample_type, y = log2_expr),
                 width = 0.55, outlier.shape = NA, color = "#0b3c5d", fill = NA) +
    geom_jitter(data = d[sample_type == "Tumor"],
                aes(x = sample_type, y = log2_expr),
                color = "red", width = 0.12, height = 0, size = 1.8, alpha = 0.65) +
    geom_jitter(data = d[sample_type == "Healthy"],
                aes(x = sample_type, y = log2_expr),
                color = "black", fill = "white", shape = 24,
                width = 0.12, height = 0, size = 2.8, stroke = 1.0, alpha = 1) +
    annotate("segment", x = 1, xend = 2,
             y = expr_ylim[2] - 0.045 * diff(expr_ylim),
             yend = expr_ylim[2] - 0.045 * diff(expr_ylim), linewidth = 0.6) +
    annotate("text", x = 1.5, y = expr_ylim[2] - 0.005 * diff(expr_ylim),
             label = p_label(expr_wilcox$p), size = 7, lineheight = 0.9) +
    scale_x_discrete(labels = c("Healthy" = "H", "Tumor" = "T")) +
    scale_y_continuous(limits = expr_ylim, breaks = expr_breaks, expand = c(0.01, 0.01)) +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 35) +
    theme(legend.position = "none", plot.margin = margin(0, 5, 0, 0))

  meth_boxplot <- ggplot() +
    geom_boxplot(data = d, aes(x = meth_percent, y = sample_type),
                 width = 0.55, outlier.shape = NA, color = "#0b3c5d", fill = NA) +
    geom_jitter(data = d[sample_type == "Tumor"],
                aes(x = meth_percent, y = sample_type),
                color = "red", width = 0, height = 0.12, size = 1.8, alpha = 0.65) +
    geom_jitter(data = d[sample_type == "Healthy"],
                aes(x = meth_percent, y = sample_type),
                color = "black", fill = "white", shape = 24,
                width = 0, height = 0.12, size = 2.8, stroke = 1.0, alpha = 1) +
    annotate("segment", x = 88, xend = 98, y = 1.35, yend = 1.35, linewidth = 0.6) +
    annotate("text", x = 93, y = 1.60,
             label = p_label(meth_wilcox$p), size = 9, lineheight = 0.9) +
    scale_x_continuous(limits = x_lim, breaks = x_breaks, expand = c(0.01, 0.01)) +
    scale_y_discrete(labels = c("Healthy" = "H", "Tumor" = "T")) +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 35) +
    theme(legend.position = "none", plot.margin = margin(5, 10, 0, 0))

  expr_grob    <- ggplotGrob(expr_boxplot)
  scatter_grob <- ggplotGrob(scatter_plot)
  meth_grob    <- ggplotGrob(meth_boxplot)
  common_x_widths      <- unit.pmax(scatter_grob$widths, meth_grob$widths)
  scatter_grob$widths  <- common_x_widths
  meth_grob$widths     <- common_x_widths
  common_y_heights     <- unit.pmax(expr_grob$heights, scatter_grob$heights)
  expr_grob$heights    <- common_y_heights
  scatter_grob$heights <- common_y_heights

  pushViewport(viewport(layout = grid.layout(
    nrow = 6, ncol = 3,
    widths  = unit.c(unit(0.8, "cm"), side_box_size, unit(1, "null")),
    heights = unit.c(unit(0.6, "cm"), unit(1.8, "cm"), unit(1.3, "cm"),
                     unit(1, "null"), side_box_size, unit(0.65, "cm"))
  )))
  grid.text(main_title, x = 0.5, y = 0.5, just = "center",
            gp = gpar(fontsize = 40, fontface = "bold"),
            vp = viewport(layout.pos.row = 2, layout.pos.col = 1:3))
  grid.text(main_subtitle, x = 0.5, y = 0.5, just = "center",
            gp = gpar(fontsize = 28),
            vp = viewport(layout.pos.row = 3, layout.pos.col = 1:3))
  grid.text("Gene expression (log2(TPM+1))", x = 0.5, y = 0.5, rot = 90, just = "center",
            gp = gpar(fontsize = 30),
            vp = viewport(layout.pos.row = 4, layout.pos.col = 1))
  grid.draw(editGrob(expr_grob,
            vp = viewport(layout.pos.row = 4, layout.pos.col = 2)))
  grid.draw(editGrob(scatter_grob,
            vp = viewport(layout.pos.row = 4, layout.pos.col = 3)))
  grid.draw(editGrob(meth_grob,
            vp = viewport(layout.pos.row = 5, layout.pos.col = 3)))
  grid.text("Probe methylation (%)", x = 0.5, y = 0.5, just = "center",
            gp = gpar(fontsize = 30),
            vp = viewport(layout.pos.row = 6, layout.pos.col = 3))
  popViewport()

  data.table(
    TF = tf_name, gene = target_gene, cancer = cancer, probe = target_probe,
    n_total = nrow(d), n_healthy, n_tumor,
    n_matrix_samples, n_na_probe,
    pearson_n = cor_res$n, pearson_r = cor_res$r, pearson_p = cor_res$p,
    expr_wilcox_p = expr_wilcox$p, meth_wilcox_p = meth_wilcox$p
  )
}

# ── Load all cancers ──────────────────────────────────────────────────────────
cat(sprintf("\nLoading %d cancers with %d cores…\n", length(target_cancers), N_CORES))
all_data <- mclapply(target_cancers, load_cancer,
                     mc.cores = N_CORES, mc.preschedule = FALSE)
names(all_data) <- target_cancers
valid   <- Filter(Negate(is.null), all_data)
n_valid <- length(valid)
cat(sprintf("%d / %d cancers have data\n", n_valid, length(target_cancers)))
if (n_valid == 0) stop("No data for any cancer type.")
valid_cancers <- names(valid)

# ── Global y-axis ─────────────────────────────────────────────────────────────
all_expr_vals   <- unlist(lapply(valid, function(d) d$log2_expr), use.names = FALSE)
global_expr_min <- min(all_expr_vals, na.rm = TRUE)
global_expr_max <- max(all_expr_vals, na.rm = TRUE)
global_expr_rng <- max(global_expr_max - global_expr_min, 1)
global_ylim_lo  <- min(max(0, global_expr_min - 0.08 * global_expr_rng),
                        expression_threshold - 0.05)
global_ylim_hi  <- max(global_expr_max + 0.18 * global_expr_rng,
                        expression_threshold + 0.05)
global_expr_ylim   <- c(global_ylim_lo, global_ylim_hi)
global_expr_breaks <- pretty(global_expr_ylim, n = 5)
cat(sprintf("Global expression y-axis: [%.3f, %.3f]\n", global_ylim_lo, global_ylim_hi))

# ── Draw PDF: page 1 = grid of plots, page 2 = summary table ─────────────────
out_pdf <- file.path(out_dir,
  paste0(tf_name, "__", target_gene, "__", target_probe, "__all_cancers_grid.pdf"))
cat("Writing PDF:", out_pdf, "\n")

N_COL         <- if (n_valid == 1) 1L else if (n_valid == 2) 2L else 3L
n_rows_outer  <- ceiling(n_valid / N_COL)
n_cols_actual <- min(n_valid, N_COL)
PDF_W <- n_cols_actual * 33
PDF_H <- n_rows_outer * 38

pdf(out_pdf, width = PDF_W, height = PDF_H, onefile = TRUE)

# ── Page 1: cancer grid ───────────────────────────────────────────────────────
grid.newpage()
pushViewport(viewport(layout = grid.layout(
  nrow    = n_rows_outer + 1L,
  ncol    = N_COL,
  heights = unit.c(unit(1.5, "cm"),
                   rep(unit(1 / n_rows_outer, "null"), n_rows_outer)),
  widths  = rep(unit(1 / N_COL, "null"), N_COL)
)))

grid.text(
  sprintf("%s — %s | %s   (%d cancer type%s)",
          tf_name, target_gene, target_probe, n_valid, ifelse(n_valid == 1, "", "s")),
  gp = gpar(fontsize = 28, fontface = "bold"),
  vp = viewport(layout.pos.row = 1, layout.pos.col = 1:N_COL)
)

stats_list <- vector("list", n_valid)
for (ci in seq_along(valid_cancers)) {
  can   <- valid_cancers[ci]
  row_i <- ceiling(ci / N_COL) + 1L
  col_i <- ((ci - 1L) %% N_COL) + 1L
  pushViewport(viewport(layout.pos.row = row_i, layout.pos.col = col_i))
  stats_list[[ci]] <- tryCatch(
    draw_cell(valid[[can]], can, global_expr_ylim, global_expr_breaks),
    error = function(e) {
      grid.text(paste("ERROR:", conditionMessage(e)), gp = gpar(col = "red", fontsize = 8))
      NULL
    }
  )
  popViewport()
}
popViewport()

# ── Page 2: summary table ─────────────────────────────────────────────────────
stats_dt <- rbindlist(Filter(Negate(is.null), stats_list), fill = TRUE)

if (nrow(stats_dt) > 0) {

  # helper: format a pre-computed value passed as string arg
  fmt_arg_p <- function(x) {
    if (is.na(x)) return("NA")
    v <- suppressWarnings(as.numeric(x))
    if (!is.na(v)) return(fmt_p(v))
    x   # already a string (e.g. "<2.2e-16")
  }

  # ── look up motif_instance_id ────────────────────────────────────────────────
  motif_id <- NA_character_
  if (file.exists(motif_file)) {
    motif_dt <- tryCatch(
      fread(cmd = paste("zcat", shQuote(motif_file))),
      error = function(e) NULL)
    if (!is.null(motif_dt)) {
      # TF column may have suffix (e.g. NRF1_mm0to2_noCGmm); match by probe + gene
      m <- motif_dt[probe == target_probe & gene == target_gene]
      if (nrow(m) == 0) m <- motif_dt[probe == target_probe]
      if (nrow(m) > 0)  motif_id <- m$motif_instance_id[1]
    }
  }

  summary_df <- data.frame(
    Metric = c("TF", "Gene", "Probe", "Motif instance", "N cancers", "Cancers",
               "Mean Pearson r", "Mean delta beta", "Mean log2FC",
               "Best Pearson p", "Best padj", "Best adj.P.Val",
               "Best expr Wilcox p", "Best meth Wilcox p",
               "Final score"),
    Value  = c(
               tf_name, target_gene, target_probe,
               ifelse(is.na(motif_id), "NA", motif_id),
               as.character(n_valid), paste(valid_cancers, collapse = ", "),
               ifelse(is.na(arg_mean_r),      as.character(round(mean(stats_dt$pearson_r, na.rm=TRUE),2)),
                                               as.character(round(arg_mean_r, 2))),
               ifelse(is.na(arg_mean_db),     "NA", as.character(round(arg_mean_db, 2))),
               ifelse(is.na(arg_mean_log2fc), "NA", as.character(round(arg_mean_log2fc, 2))),
               fmt_arg_p(arg_best_pear_p),
               fmt_arg_p(arg_best_padj),
               fmt_arg_p(arg_best_adj_p),
               fmt_arg_p(arg_best_expr_wp),
               fmt_arg_p(arg_best_meth_wp),
               as.character(round(final_score, 2)))
  )

  grid.newpage()
  grid.text(paste0(tf_name, " — ", target_gene, " | ", target_probe),
            x = 0.5, y = 0.96, just = "center",
            gp = gpar(fontsize = 22, fontface = "bold"))
  tg <- tableGrob(summary_df, rows = NULL,
                  theme = ttheme_minimal(
                    core    = list(fg_params = list(fontsize = 28),
                                   padding   = unit(c(16, 10), "mm")),
                    colhead = list(fg_params = list(col = "white", fontsize = 30,
                                                    fontface = "bold"),
                                   bg_params = list(fill = "#2c3e50", col = NA),
                                   padding   = unit(c(16, 10), "mm"))
                  ))
  # centre the table on the page
  tg_vp <- viewport(x = 0.5, y = 0.45, width = 0.85, height = 0.85,
                    just = c("centre", "centre"))
  pushViewport(tg_vp)
  grid.draw(tg)
  popViewport()

  # also write per-cancer stats TSV
  fwrite(stats_dt, stats_out, sep = "\t")
  cat("Per-cancer stats written to:", stats_out, "\n")
}

dev.off()
cat("Done. PDF written to:", out_pdf, "\n")
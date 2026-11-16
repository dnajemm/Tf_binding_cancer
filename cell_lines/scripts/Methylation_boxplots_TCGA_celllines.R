#!/usr/bin/env Rscript
# Usage: Rscript ./cell_lines/scripts/Methylation_boxplots_TCGA_celllines.R <cancer_type> [CellLine_healthy|tumor ...]
#   e.g.: Rscript s./cell_lines/scripts/Methylation_boxplots_TCGA_celllines.R COAD NCM460_healthy HCT116_tumor

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop(
  "Usage: Rscript ./cell_lines/scripts/Methylation_boxplots_TCGA_celllines.R <cancer_type> [CellLine_healthy|tumor ...]\n",
  "  e.g.: Rscript ./cell_lines/scripts/Methylation_boxplots_TCGA_celllines.R COAD NCM460_healthy HCT116_tumor\n",
  "  If no cell line arguments are given, all cell lines get distinct colours."
)
cancer_type <- args[1]
cat("Cancer type:", cancer_type, "\n")

cl_type_args <- args[-1]
cl_type_map  <- NULL
if (length(cl_type_args) > 0) {
  cl_type_map <- data.table(arg = cl_type_args)
  cl_type_map[, type      := sub(".*_([^_]+)$", "\\1", arg)]
  cl_type_map[, cell_line := sub("_[^_]+$",    "",     arg)]
  bad <- cl_type_map[!type %in% c("healthy", "tumor")]
  if (nrow(bad) > 0) stop("Cell line type must be 'healthy' or 'tumor'. Bad args: ",
                           paste(bad$arg, collapse = ", "))
  cat("Cell line types specified:\n")
  print(cl_type_map[, .(cell_line, type)])
}

project_root <- "/data/najemd/TF_binding_cancer"

# ---------------------------------------------------------------------------
# Step 1 — Get significant probes for this cancer type
# ---------------------------------------------------------------------------
integration_file <- file.path(project_root,
  "results/integration/meth_expr_consistent_pairs.tsv")

sig <- fread(integration_file, header = TRUE,
             select = c("probe", "gene", "cancer", "direction"))
sig <- sig[cancer == cancer_type]
gene_dir   <- sig[, .N, by = .(probe, gene, direction)][order(-N)][!duplicated(probe)]
sig_probes <- unique(sig$probe)
cat("Significant probes for", cancer_type, ":", length(sig_probes), "\n")

# ---------------------------------------------------------------------------
# Step 2 — Load TCGA methylation (beta matrix)
# ---------------------------------------------------------------------------
# Load both BANP and NRF1 beta matrices — significant probes may be in either
beta_files <- c(
  file.path(project_root, "methylation/beta_matrices_BANP",
            paste0(cancer_type, "_BANP_beta_matrix.tsv.gz")),
  file.path(project_root, "methylation/beta_matrices_NRF1",
            paste0(cancer_type, "_NRF1_beta_matrix.tsv.gz"))
)
beta_files <- beta_files[file.exists(beta_files)]
if (length(beta_files) == 0) stop("No beta matrix found for: ", cancer_type)
cat("Beta matrices found:", paste(basename(beta_files), collapse = ", "), "\n")

# Use first file to get sample list (both matrices have the same samples)
beta_header  <- fread(beta_files[1], nrows = 0)
probe_col    <- names(beta_header)[1]
beta_samples <- setdiff(names(beta_header), probe_col)

beta_dt <- data.table(col = beta_samples)
beta_dt[, st := fifelse(grepl("_tumor$",   col), "Tumor",
                fifelse(grepl("_healthy$", col), "Healthy", "other"))]
beta_dt <- beta_dt[st != "other"]

cat("Beta matrix samples:", nrow(beta_dt),
    "(Tumor:", sum(beta_dt$st == "Tumor"),
    "| Healthy:", sum(beta_dt$st == "Healthy"), ")\n")

# Read significant probes from both matrices and combine (drop duplicates)
beta_list <- lapply(beta_files, function(f) {
  b <- fread(f, header = TRUE)
  setnames(b, 1L, "probe")
  b[probe %in% sig_probes]
})
beta <- rbindlist(beta_list)
beta <- beta[!duplicated(probe)]
cat("Probes found in beta matrices:", nrow(beta), "/", length(sig_probes), "\n")
missing_probes <- sig_probes[!sig_probes %in% beta$probe]
if (length(missing_probes) > 0) {
  cat("Probes NOT in any beta matrix (", length(missing_probes), "):\n")
  miss_info <- gene_dir[probe %in% missing_probes, .(probe, gene, direction)]
  print(miss_info)
}

keep_cols <- intersect(c("probe", beta_dt$col), names(beta))
beta_long <- melt(beta[, ..keep_cols],
                  id.vars       = "probe",
                  variable.name = "sample_id",
                  value.name    = "beta")
beta_long <- merge(beta_long, beta_dt[, .(col, st)], by.x = "sample_id", by.y = "col")

# Auto-detect whether values are 0–1 fractions or already 0–100 percentages
max_val <- max(beta_long$beta, na.rm = TRUE)
if (max_val <= 1) {
  beta_long[, meth := beta * 100]
  cat("Beta values in 0-1 range → converted to %\n")
} else {
  beta_long[, meth := beta]
  cat("Beta values already in % range (max =", round(max_val, 2), ") → used as-is\n")
}

beta_long[, group  := st]
beta_long[, source := "TCGA"]

cat("TCGA methylation rows:", nrow(beta_long), "\n")

# ---------------------------------------------------------------------------
# Step 3 — Load RRBS methylation for cell lines via bedtools intersect
# ---------------------------------------------------------------------------
cell_line_meth_dir <- file.path(project_root, "cell_lines", cancer_type, "methylation")
probe_bed_file     <- file.path(project_root,
                       "methylation/annotated_methylation_data_probes_filtered.bed")

# Extract coordinates for significant probes
sig_bed <- fread(probe_bed_file, header = FALSE,
                 col.names = c("chr", "start", "end", "probe"))
sig_bed <- sig_bed[probe %in% sig_probes]
cat("Significant probes with RRBS coordinates:", nrow(sig_bed), "\n")

tmp_probe_bed <- tempfile(fileext = ".bed")
fwrite(sig_bed, tmp_probe_bed, sep = "\t", col.names = FALSE)

bs_files <- list.files(cell_line_meth_dir,
                       pattern    = "^(RRBS|WGBS)_.*_meth\\.bed\\.gz$",
                       full.names = TRUE)
if (length(bs_files) == 0) stop("No RRBS/WGBS meth.bed.gz files found in: ", cell_line_meth_dir)

rrbs_dt <- data.table(path = bs_files)
rrbs_dt[, basename_ := basename(path)]
rrbs_dt[, seq_type  := sub("^(RRBS|WGBS)_.*", "\\1", basename_)]
rrbs_dt[, cell_line := sub("^(?:RRBS|WGBS)_(.+)_[0-9]+_meth\\.bed\\.gz$", "\\1", basename_, perl = TRUE)]
rrbs_dt[, replicate := as.integer(sub("^(?:RRBS|WGBS)_.+_([0-9]+)_meth\\.bed\\.gz$", "\\1", basename_, perl = TRUE))]

cat("Bisulfite sequencing files found:\n")
print(rrbs_dt[, .(seq_type, cell_line, replicate, path)])

if (!is.null(cl_type_map)) {
  cl_names <- cl_type_map$cell_line
} else {
  cl_names <- sort(unique(rrbs_dt$cell_line))
}

cl_meth_list <- lapply(seq_len(nrow(rrbs_dt)), function(i) {
  tmp_out <- tempfile(fileext = ".txt")
  cmd <- sprintf(
    "bedtools intersect -a %s -b <(zcat %s) -wa -wb | awk '{print $4\"\\t\"$8\"\\t\"$9}' > %s",
    tmp_probe_bed, rrbs_dt$path[i], tmp_out
  )
  system(paste("bash -c", shQuote(cmd)))
  if (file.info(tmp_out)$size == 0) {
    cat("WARNING: No intersections for", rrbs_dt$basename_[i], "(", rrbs_dt$seq_type[i], ")\n")
    return(NULL)
  }
  dt <- fread(tmp_out, header = FALSE, col.names = c("probe", "coverage", "meth"))
  dt[, cell_line := rrbs_dt$cell_line[i]]
  dt[, replicate := rrbs_dt$replicate[i]]
  dt
})

cl_meth <- rbindlist(cl_meth_list[!sapply(cl_meth_list, is.null)])
cl_meth[, group  := cell_line]
cl_meth[, source := "CellLine"]

cat("Cell line methylation rows:", nrow(cl_meth), "\n")
for (cl in cl_names) cat(" ", cl, ":", sum(cl_meth$cell_line == cl), "probe-replicate pairs\n")

# ---------------------------------------------------------------------------
# Step 4 — Combine TCGA and cell line data
# ---------------------------------------------------------------------------
common_cols <- c("probe", "meth", "group", "source")
combined <- rbindlist(list(
  beta_long[, ..common_cols],
  cl_meth[,  ..common_cols]
))

combined[, group := factor(group, levels = c("Healthy", "Tumor", cl_names))]
combined <- merge(combined, gene_dir[, .(probe, gene, direction)], by = "probe", all.x = TRUE)

cat("Combined rows:", nrow(combined), "\n")
cat("Groups:", paste(levels(combined$group), collapse = " | "), "\n")

# ---------------------------------------------------------------------------
# Step 4b — Concordance: cell line direction vs TCGA direction
# ---------------------------------------------------------------------------
concordance_label <- NULL
gene_concordance  <- NULL

if (!is.null(cl_type_map) &&
    all(c("healthy", "tumor") %in% cl_type_map$type)) {

  cl_data_tmp <- combined[source == "CellLine"]
  cl_typed <- merge(cl_data_tmp,
                    cl_type_map[, .(group = cell_line, cl_type = type)],
                    by = "group")

  cl_med  <- cl_typed[, .(med = median(meth, na.rm = TRUE)), by = .(probe, cl_type)]
  cl_wide <- dcast(cl_med, probe ~ cl_type, value.var = "med")

  cl_wide[, cl_dir := fifelse(tumor > healthy, "Hyper",
                      fifelse(tumor < healthy, "Hypo", "Flat"))]

  cl_wide <- merge(cl_wide, gene_dir[, .(probe, direction)], by = "probe", all.x = TRUE)
  cl_wide[, tcga_dir := fifelse(grepl("Hyper", direction), "Hyper",
                        fifelse(grepl("Hypo",  direction), "Hypo", NA_character_))]

  cl_wide[, concordant := cl_dir == tcga_dir]

  n_conc  <- sum(cl_wide$concordant, na.rm = TRUE)
  n_total <- sum(!is.na(cl_wide$concordant))
  pct     <- round(100 * n_conc / n_total, 1)

  # Flag probes not covered by RRBS (against ALL significant probes, not just cl_wide)
  covered_probes  <- unique(cl_meth$probe)
  in_beta         <- unique(beta_long$probe)
  cl_wide[, covered := probe %in% covered_probes]

  n_not_cov      <- sum(!sig_probes %in% covered_probes)
  n_not_in_tcga  <- sum(!sig_probes %in% in_beta)
  concordance_label <- paste0(
    n_conc, "/", n_total, " covered probes (", pct, "%) follow the same direction as TCGA",
    "  |  ", n_not_cov, " not covered by RRBS",
    if (n_not_in_tcga > 0) paste0("  |  ", n_not_in_tcga, " absent from TCGA beta matrix") else ""
  )
  cat("Concordance:", concordance_label, "\n")
  print(cl_wide[, .(probe, cl_dir, tcga_dir, concordant, covered)])

  gene_concordance <- cl_wide[, .(probe, concordant, covered)]
}

# ---------------------------------------------------------------------------
# Step 5 — Build labels and plot
# ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(ggplot2))

if (!is.null(gene_concordance)) {
  combined <- merge(combined, gene_concordance, by = "probe", all.x = TRUE)
} else {
  combined[, concordant := NA]
  combined[, covered    := NA]
}

# Build one label per probe, uniform across all rows
dir_dt   <- unique(combined[, .(probe, gene, direction)])
conc_dt  <- unique(combined[, .(probe, concordant, covered)])[!duplicated(probe)]
label_dt <- merge(dir_dt, conc_dt, by = "probe", all.x = TRUE)
label_dt[, tag := fifelse(!is.na(covered) & !covered,              "[0]",
                  fifelse(!is.na(concordant) & concordant == TRUE,  "[+]",
                  fifelse(!is.na(concordant) & concordant == FALSE, "[-]",
                                                                     "")))]
label_dt[, probe_label := paste0(probe, "\n", gene, "\n(", direction, ") ", tag)]

combined[, probe_label := label_dt$probe_label[match(probe, label_dt$probe)]]
combined[, group := factor(group, levels = c("Healthy", "Tumor", cl_names))]

# Split AFTER labels are finalised
tcga_data <- combined[source == "TCGA"]
cl_data   <- combined[source == "CellLine"]

# Colours
fallback_palette <- c("#762A83", "#01665E", "#8C510A", "#4393C3",
                      "#F46D43", "#74ADD1", "#A50026", "#313695")
fallback_i <- 1L
cl_colours <- sapply(cl_names, function(cl) {
  if (!is.null(cl_type_map) && cl %in% cl_type_map$cell_line) {
    type <- cl_type_map[cell_line == cl, type]
    if (type == "healthy") "#1B7837" else "#D6604D"
  } else {
    col <- fallback_palette[fallback_i]
    fallback_i <<- fallback_i + 1L
    col
  }
})
group_colours <- c(Healthy = "#1B7837", Tumor = "#D6604D", cl_colours)

p <- ggplot(combined, aes(x = group, y = meth, fill = group, colour = group)) +

  geom_boxplot(data = tcga_data,
               outlier.shape = NA,
               width = 0.6, alpha = 0.4, linewidth = 0.4) +

  geom_point(data = cl_data, size = 2.5, alpha = 0.9) +

  stat_summary(data = cl_data, fun = median, geom = "crossbar",
               width = 0.4, linewidth = 0.5) +

  geom_jitter(data = tcga_data, width = 0.2, size = 0.3, alpha = 0.15, shape = 16) +

  scale_fill_manual(values   = group_colours) +
  scale_colour_manual(values = group_colours) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  coord_cartesian(ylim = c(0, 100)) +

  facet_wrap(~ probe_label, ncol = 4) +

  labs(x = NULL,
       y = "Methylation (%)",
       title    = paste0(cancer_type, " — significant probes: methylation by group"),
       subtitle = if (!is.null(concordance_label)) concordance_label else NULL) +

  theme_bw(base_size = 10) +
  theme(
    legend.position  = "none",
    axis.text.x      = element_text(angle = 35, hjust = 1, size = 8),
    strip.background = element_rect(fill = "grey92"),
    strip.text       = element_text(face = "bold", size = 8),
    panel.grid.minor = element_blank()
  )

out_dir <- file.path(project_root, "cell_lines", cancer_type, "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_pdf <- file.path(out_dir, paste0(cancer_type, "_methylation_boxplots.pdf"))
ggsave(out_pdf, plot = p,
       width  = max(12, 4 * 2.5),
       height = ceiling(length(sig_probes) / 4) * 3.5,
       units  = "in", device = "pdf")
cat("Plot saved to:", out_pdf, "\n")
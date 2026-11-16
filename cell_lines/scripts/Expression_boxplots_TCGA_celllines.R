#!/usr/bin/env Rscript
# Usage: Rscript Expression_boxplots_TCGA_celllines.R <cancer_type>

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop(
  "Usage: Rscript Expression_boxplots_TCGA_celllines.R <cancer_type> [CellLine_healthy|tumor ...]\n",
  "  e.g.: Rscript script.R COAD NCM460_healthy HCT116_tumor\n",
  "  e.g.: Rscript script.R PRAD LNCaP_tumor PrEC_healthy\n",
  "  If no cell line arguments are given, all cell lines get distinct colours."
)
cancer_type <- args[1]
cat("Cancer type:", cancer_type, "\n")

# Parse optional cell line type arguments (format: CellLineName_healthy or CellLineName_tumor)
# Cell line name may contain underscores; type is always the LAST underscore-delimited token
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
# Step 1 — Load TCGA expression data and filter to matched samples
# ---------------------------------------------------------------------------
tpm_file <- file.path(project_root, "expression/gene_expression_matrix_protein_coding.tsv")

beta_file <- file.path(project_root, "methylation/beta_matrices_BANP",
                       paste0(cancer_type, "_BANP_beta_matrix.tsv.gz"))
if (!file.exists(beta_file)) {
  beta_file <- file.path(project_root, "methylation/beta_matrices_NRF1",
                         paste0(cancer_type, "_NRF1_beta_matrix.tsv.gz"))
}
cat("Reading beta matrix header from:", beta_file, "\n")

beta_header <- fread(beta_file, nrows = 0)
beta_samples <- setdiff(names(beta_header), names(beta_header)[1])

beta_dt <- data.table(col = beta_samples)
beta_dt[, patient  := sub("_[a-z]+$", "", col)]
beta_dt[, st       := fifelse(grepl("_tumor$",   col), "Tumor",
                     fifelse(grepl("_healthy$", col), "Healthy", "other"))]
beta_dt <- beta_dt[st != "other"]
beta_dt[, pair_key := paste(patient, st, sep = "__")]

cat("Beta matrix samples:", nrow(beta_dt),
    "(Tumor:", sum(beta_dt$st == "Tumor"),
    "| Healthy:", sum(beta_dt$st == "Healthy"), ")\n")

cat("Loading TPM matrix...\n")
tpm <- fread(tpm_file, header = TRUE, check.names = FALSE)
setnames(tpm, 1L, "sample_id")

tpm[, cancer   := sub("^TCGA-([A-Z]+)-.*", "\\1", sample_id)]
tpm[, st       := fifelse(grepl("-01A_1", sample_id), "Tumor",
                  fifelse(grepl("-11A_1", sample_id), "Healthy", "other"))]
tpm[, patient  := sub("^TCGA-([A-Z]+-TCGA-[A-Z0-9]+-[A-Z0-9]+)-[0-9]+[A-Z]_.*",
                       "\\1", sample_id)]
tpm[, pair_key := paste(patient, st, sep = "__")]

tpm <- tpm[cancer == cancer_type & st != "other" & pair_key %in% beta_dt$pair_key]

cat("TCGA samples kept:", nrow(tpm),
    "(Tumor:", sum(tpm$st == "Tumor"),
    "| Healthy:", sum(tpm$st == "Healthy"), ")\n")

# ---------------------------------------------------------------------------
# Step 2 — Load cell line expression data
# ---------------------------------------------------------------------------
cell_line_dir <- file.path(project_root, "cell_lines", cancer_type, "expression")

all_files <- list.files(cell_line_dir, pattern = "^RNA_.*_genes_count\\.txt$",
                        full.names = TRUE)
if (length(all_files) == 0) stop("No expression files found in: ", cell_line_dir)

file_dt <- data.table(path = all_files)
file_dt[, basename_ := basename(path)]
file_dt[, cell_line := sub("^RNA_(.+)_[0-9]+_genes_count\\.txt$", "\\1", basename_)]
file_dt[, replicate := as.integer(sub("^RNA_.+_([0-9]+)_genes_count\\.txt$", "\\1", basename_))]

cat("Cell line files found:\n")
print(file_dt[, .(cell_line, replicate, path)])

cl_long <- rbindlist(lapply(seq_len(nrow(file_dt)), function(i) {
  dt <- fread(file_dt$path[i], header = TRUE, select = c("geneName", "TPM"))
  dt[, cell_line := file_dt$cell_line[i]]
  dt[, replicate := file_dt$replicate[i]]
  dt
}))

if (!is.null(cl_type_map)) {
  cl_names <- cl_type_map$cell_line
} else {
  cl_names <- sort(unique(cl_long$cell_line))
}
cat("Cell line rows loaded:", nrow(cl_long), "\n")
for (cl in cl_names) cat(" ", cl, ":", sum(cl_long$cell_line == cl), "\n")
cat("Unique genes (cell lines):", uniqueN(cl_long$geneName), "\n")

# ---------------------------------------------------------------------------
# Step 3 — Select significant genes and build a unified long-format table
# ---------------------------------------------------------------------------
integration_file <- file.path(project_root,
  "results/integration/meth_expr_consistent_pairs_cancermine_labeled.tsv")

sig <- fread(integration_file, header = TRUE, select = c("gene", "cancer", "direction"))
sig <- sig[cancer == cancer_type]
gene_dir  <- sig[, .N, by = .(gene, direction)][order(-N)][!duplicated(gene)]
sig_genes <- unique(sig$gene)
cat("Significant genes for", cancer_type, ":", length(sig_genes), "\n")

meta_cols <- c("sample_id", "cancer", "st", "patient", "pair_key")
gene_cols  <- intersect(sig_genes, setdiff(names(tpm), meta_cols))
cat("Significant genes found in TCGA matrix:", length(gene_cols), "\n")

tpm_long <- melt(tpm[, c(meta_cols, gene_cols), with = FALSE],
                 id.vars       = meta_cols,
                 variable.name = "gene",
                 value.name    = "TPM")
tpm_long[, group := st]
tpm_long[, source := "TCGA"]

cl_sig <- cl_long[geneName %in% sig_genes]
cl_sig[, group  := cell_line]
cl_sig[, source := "CellLine"]
setnames(cl_sig, "geneName", "gene")

common_cols <- c("gene", "TPM", "group", "source")
combined <- rbindlist(list(
  tpm_long[, ..common_cols],
  cl_sig[,  ..common_cols]
))

combined[, group := factor(group, levels = c("Healthy", "Tumor", cl_names))]
combined <- merge(combined, gene_dir[, .(gene, direction)], by = "gene", all.x = TRUE)

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

  cl_med  <- cl_typed[, .(med = median(TPM, na.rm = TRUE)), by = .(gene, cl_type)]
  cl_wide <- dcast(cl_med, gene ~ cl_type, value.var = "med")

  cl_wide[, cl_dir := fifelse(tumor > healthy, "Up",
                      fifelse(tumor < healthy, "Down", "Flat"))]

  cl_wide <- merge(cl_wide, gene_dir[, .(gene, direction)], by = "gene", all.x = TRUE)
  cl_wide[, tcga_dir := fifelse(grepl("Down", direction), "Down",
                        fifelse(grepl("Up",   direction), "Up", NA_character_))]

  cl_wide[, concordant := cl_dir == tcga_dir]

  n_conc  <- sum(cl_wide$concordant, na.rm = TRUE)
  n_total <- sum(!is.na(cl_wide$concordant))
  pct     <- round(100 * n_conc / n_total, 1)

  # Check expression in cell lines
  cl_expr <- cl_long[geneName %in% sig_genes,
                     .(median_cl = median(TPM, na.rm = TRUE)), by = geneName]
  setnames(cl_expr, "geneName", "gene")
  cl_wide <- merge(cl_wide, cl_expr, by = "gene", all.x = TRUE)
  cl_wide[, expressed := median_cl >= 1]

  n_not_expr <- sum(!cl_wide$expressed, na.rm = TRUE)
  concordance_label <- paste0(
    n_conc, "/", n_total, " genes (", pct, "%) follow the same direction as TCGA",
    "  |  ", n_not_expr, " gene(s) not expressed in cell lines (median TPM < 1)"
  )
  cat("Concordance:", concordance_label, "\n")
  print(cl_wide[, .(gene, cl_dir, tcga_dir, concordant, expressed, median_cl)])

  gene_concordance <- cl_wide[, .(gene, concordant, expressed)]
}

# ---------------------------------------------------------------------------
# Step 4 — Build final gene labels, then split data for plotting
# ---------------------------------------------------------------------------

# Merge concordance into combined
if (!is.null(gene_concordance)) {
  combined <- merge(combined, gene_concordance, by = "gene", all.x = TRUE)
} else {
  combined[, concordant := NA]
  combined[, expressed  := NA]
}

# Build one label per gene (uniform across all rows)
dir_dt     <- unique(combined[, .(gene, direction)])
conc_dt    <- unique(combined[, .(gene, concordant, expressed)])
conc_dt    <- conc_dt[!duplicated(gene)]
emoji_dt   <- merge(dir_dt, conc_dt, by = "gene", all.x = TRUE)
emoji_dt[, emoji := fifelse(!is.na(expressed) & !expressed,           "[0]",
                    fifelse(!is.na(concordant) & concordant == TRUE,  "[+]",
                    fifelse(!is.na(concordant) & concordant == FALSE, "[-]",
                                                                       "")))]
emoji_dt[, gene_label := paste0(gene, "\n(", direction, ") ", emoji)]

# Apply uniform label to every row in combined
combined[, gene_label := emoji_dt$gene_label[match(gene, emoji_dt$gene)]]

# NOW split — after gene_label is finalised
tcga_data <- combined[source == "TCGA"]
cl_data   <- combined[source == "CellLine"]

# ---------------------------------------------------------------------------
# Step 5 — Plot
# ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(ggplot2))

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

p <- ggplot(combined, aes(x = group, y = log2(TPM + 1), fill = group, colour = group)) +

  geom_boxplot(data      = tcga_data,
               outlier.shape = NA,
               width     = 0.6,
               alpha     = 0.4,
               linewidth = 0.4) +

  geom_point(data  = cl_data,
             size  = 2.5,
             alpha = 0.9) +

  stat_summary(data  = cl_data,
               fun   = median,
               geom  = "crossbar",
               width = 0.4,
               linewidth = 0.5) +

  geom_jitter(data   = tcga_data,
              width  = 0.2,
              size   = 0.3,
              alpha  = 0.15,
              shape  = 16) +

  scale_fill_manual(values   = group_colours) +
  scale_colour_manual(values = group_colours) +

  facet_wrap(~ gene_label, scales = "free_y", ncol = 4) +

  labs(x = NULL,
       y = expression(log[2](TPM + 1)),
       title    = paste0(cancer_type, " — significant genes: expression by group"),
       subtitle = if (!is.null(concordance_label)) concordance_label else NULL) +

  theme_bw(base_size = 10) +
  theme(
    legend.position  = "none",
    axis.text.x      = element_text(angle = 35, hjust = 1, size = 8),
    strip.background = element_rect(fill = "grey92"),
    strip.text       = element_text(face = "bold", size = 9),
    panel.grid.minor = element_blank()
  )

out_dir <- file.path(project_root, "cell_lines", cancer_type, "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_pdf <- file.path(out_dir, paste0(cancer_type, "_expression_boxplots.pdf"))
ggsave(out_pdf, plot = p,
       width  = max(12, 4 * 2.5),
       height = ceiling(length(sig_genes) / 4) * 3.5,
       units  = "in", device = "pdf")
cat("Plot saved to:", out_pdf, "\n")
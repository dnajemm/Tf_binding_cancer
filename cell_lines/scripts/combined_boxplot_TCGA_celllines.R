#!/usr/bin/env Rscript
# Combined methylation + expression boxplots (TCGA + cell lines)
# One facet panel per probe-gene pair:
#   left half  -> methylation (%)          -- left y-axis
#   right half -> expression (log2 TPM+1) -- right y-axis (sec_axis)
#   dashed vertical line separates the two halves
# Layout: 2 columns per row
#
# Usage:
#   Rscript /data/najemd/TF_binding_cancer/cell_lines/scripts/combined_boxplot_TCGA_celllines.R <cancer_type> [meth_prefix:expr_prefix:healthy|tumor ...]
#   e.g.: Rscript /data/najemd/TF_binding_cancer/cell_lines/scripts/combined_boxplot_TCGA_celllines.R COAD RRBS_NCM460:RNA_NCM460:healthy RRBS_HCT116:RNA_HCT116:tumor
#   e.g.: Rscript /data/najemd/TF_binding_cancer/cell_lines/scripts/combined_boxplot_TCGA_celllines.R PRAD WGBS_PreC:RNA_PreC:healthy WGBS_LNCaP:RNA_LNCaP:tumor
#   e.g.: Rscript /data/najemd/TF_binding_cancer/cell_lines/scripts/combined_boxplot_TCGA_celllines.R BRCA WGBS_MCF10A:RNA_MCF10A:healthy RRBS_HCC1954:RNA_HCC1954:tumor
#   e.g.: Rscript /data/najemd/TF_binding_cancer/cell_lines/scripts/combined_boxplot_TCGA_celllines.R BRCA WGBS_HMEC:RNA_HMEC:healthy WGBS_HCC1954:RNA_HCC1954:tumor

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  if (!requireNamespace("ggh4x", quietly = TRUE))
    install.packages("ggh4x", repos = "https://cloud.r-project.org", quiet = TRUE)
  library(ggh4x)
})

# ===========================================================================
# 0 - Arguments
# ===========================================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop(
  "Usage: Rscript ... <cancer_type> [meth_prefix:expr_prefix:healthy|tumor ...]\n",
  "  e.g.: Rscript ... BRCA WGBS_MCF10A:RNA_MCF10A:healthy RRBS_HCC1954:RNA_HCC1954:tumor\n"
)
cancer_type <- args[1]
sample_args <- if (length(args) >= 2) args[-1] else character(0)
cat("Cancer type:", cancer_type, "\n")

cl_type_map <- NULL
if (length(sample_args) > 0) {
  parts <- strsplit(sample_args, ":")
  cl_type_map <- rbindlist(lapply(parts, function(p) {
    if (length(p) != 3) stop(
      "Each sample arg must be meth_prefix:expr_prefix:type  (got: ",
      paste(p, collapse = ":"), ")")
    list(meth_prefix = p[1], expr_prefix = p[2], type = p[3])
  }))
  cl_type_map[, cell_line := sub("^RNA_", "", expr_prefix)]
  bad <- cl_type_map[!type %in% c("healthy", "tumor")]
  if (nrow(bad) > 0) stop(
    "type must be 'healthy' or 'tumor'. Bad args: ",
    paste(bad$meth_prefix, collapse = ", "))
  cat("Cell-line map:\n"); print(cl_type_map[, .(meth_prefix, expr_prefix, cell_line, type)])
}

project_root <- "/data/najemd/TF_binding_cancer"

# ===========================================================================
# Step 1 - Load significant probe-gene pairs
# ===========================================================================
integration_file <- file.path(project_root,
  "results/integration/meth_expr_consistent_pairs.tsv")

sig_all  <- fread(integration_file, header = TRUE,
                  select = c("probe", "gene", "cancer", "direction"))
sig      <- sig_all[cancer == cancer_type]
pair_dir <- sig[, .N, by = .(probe, gene, direction)][order(-N)][!duplicated(probe)]
sig_probes <- unique(pair_dir$probe)
sig_genes  <- unique(pair_dir$gene)
pair_dir[, pair_rank := .I]
cat("Significant probe-gene pairs:", nrow(pair_dir),
    "| Unique probes:", length(sig_probes),
    "| Unique genes:", length(sig_genes), "\n")

# ===========================================================================
# Step 2 - TCGA methylation (beta matrices)
# ===========================================================================
beta_files <- c(
  file.path(project_root, "methylation/beta_matrices_BANP",
            paste0(cancer_type, "_BANP_beta_matrix.tsv.gz")),
  file.path(project_root, "methylation/beta_matrices_NRF1",
            paste0(cancer_type, "_NRF1_beta_matrix.tsv.gz"))
)
beta_files <- beta_files[file.exists(beta_files)]
if (length(beta_files) == 0) stop("No beta matrix found for: ", cancer_type)
cat("Beta matrix files:", paste(basename(beta_files), collapse = ", "), "\n")

beta_header  <- fread(beta_files[1], nrows = 0)
beta_samp_dt <- data.table(col = setdiff(names(beta_header), names(beta_header)[1]))
beta_samp_dt[, st := fifelse(grepl("_tumor$",   col), "Tumor",
                     fifelse(grepl("_healthy$", col), "Healthy", "other"))]
beta_samp_dt <- beta_samp_dt[st != "other"]
beta_samp_dt[, patient  := sub("_[a-z]+$", "", col)]
beta_samp_dt[, pair_key := paste(patient, st, sep = "__")]
cat("Beta samples:", nrow(beta_samp_dt),
    "(Tumor:", sum(beta_samp_dt$st == "Tumor"),
    "| Healthy:", sum(beta_samp_dt$st == "Healthy"), ")\n")

beta_list <- lapply(beta_files, function(f) {
  b <- fread(f, header = TRUE); setnames(b, 1L, "probe"); b[probe %in% sig_probes]
})
beta_raw <- rbindlist(beta_list)[!duplicated(probe)]
cat("Probes found in beta matrices:", nrow(beta_raw), "/", length(sig_probes), "\n")
missing_probes <- sig_probes[!sig_probes %in% beta_raw$probe]
if (length(missing_probes) > 0) {
  cat("Probes absent from beta matrices (", length(missing_probes), "):\n")
  print(pair_dir[probe %in% missing_probes, .(probe, gene, direction)])
}

keep_beta <- intersect(c("probe", beta_samp_dt$col), names(beta_raw))
beta_long <- melt(beta_raw[, ..keep_beta],
                  id.vars       = "probe",
                  variable.name = "sample_id",
                  value.name    = "raw_beta")
beta_long <- merge(beta_long, beta_samp_dt[, .(col, st)],
                   by.x = "sample_id", by.y = "col")
max_beta  <- max(beta_long$raw_beta, na.rm = TRUE)
beta_long[, meth := if (max_beta <= 1) raw_beta * 100 else raw_beta]
beta_long[, c("raw_beta", "sample_id") := NULL]
beta_long[, group  := st]
beta_long[, source := "TCGA"]
cat("TCGA methylation rows:", nrow(beta_long), "\n")

# ===========================================================================
# Step 3 - Cell-line methylation (bedtools intersect against RRBS/WGBS)
# ===========================================================================
cell_line_meth_dir <- file.path(project_root, "cell_lines", cancer_type, "methylation")
probe_bed_file     <- file.path(project_root,
  "methylation/annotated_methylation_data_probes_filtered.bed")

sig_bed <- fread(probe_bed_file, header = FALSE,
                 col.names = c("chr", "start", "end", "probe"))[probe %in% sig_probes]
cat("Significant probes with RRBS coordinates:", nrow(sig_bed), "\n")
tmp_bed <- tempfile(fileext = ".bed")
fwrite(sig_bed, tmp_bed, sep = "\t", col.names = FALSE)

bs_files <- list.files(cell_line_meth_dir,
                       pattern    = "^(RRBS|WGBS)_.*_meth\\.bed\\.gz$",
                       full.names = TRUE)
if (length(bs_files) == 0) stop("No RRBS/WGBS files in: ", cell_line_meth_dir)
if (!is.null(cl_type_map)) {
  pfx_pattern <- paste0("^(", paste(cl_type_map$meth_prefix, collapse = "|"), ")_")
  bs_files    <- bs_files[grepl(pfx_pattern, basename(bs_files))]
}
cat("Bisulfite files matched:", length(bs_files), "\n")

rrbs_dt <- data.table(path = bs_files)
rrbs_dt[, basename_ := basename(path)]
if (!is.null(cl_type_map)) {
  rrbs_dt[, cell_line := sapply(basename_, function(f)
    cl_type_map$cell_line[which(startsWith(f, paste0(cl_type_map$meth_prefix, "_")))[1]])]
} else {
  rrbs_dt[, cell_line := sub("^(?:RRBS|WGBS)_(.+)_[0-9]+_meth\\.bed\\.gz$",
                              "\\1", basename_, perl = TRUE)]
}
rrbs_dt[, replicate := as.integer(sub("^(?:RRBS|WGBS)_.+_([0-9]+)_meth\\.bed\\.gz$",
                                       "\\1", basename_, perl = TRUE))]
cat("Bisulfite files found:\n"); print(rrbs_dt[, .(cell_line, replicate, basename_)])

cl_names <- if (!is.null(cl_type_map)) cl_type_map$cell_line else sort(unique(rrbs_dt$cell_line))

cl_meth <- rbindlist(Filter(Negate(is.null), lapply(seq_len(nrow(rrbs_dt)), function(i) {
  tmp_out <- tempfile(fileext = ".txt")
  system(paste("bash -c", shQuote(sprintf(
    "bedtools intersect -a %s -b <(zcat %s) -wa -wb | awk '{print $4\"\\t\"$8\"\\t\"$9}' > %s",
    tmp_bed, rrbs_dt$path[i], tmp_out))))
  if (file.info(tmp_out)$size == 0) {
    cat("WARNING: no intersections for", rrbs_dt$basename_[i], "\n"); return(NULL)
  }
  dt <- fread(tmp_out, header = FALSE, col.names = c("probe", "coverage", "meth"))
  dt[, cell_line := rrbs_dt$cell_line[i]]
  dt[, replicate := rrbs_dt$replicate[i]]
  dt
})))
cl_meth[, group  := cell_line]
cl_meth[, source := "CellLine"]
cat("Cell-line methylation rows:", nrow(cl_meth), "\n")
for (cl in cl_names) cat(" ", cl, ":", sum(cl_meth$cell_line == cl),
                          "probe-replicate pairs\n")

# ===========================================================================
# Step 4 - TCGA expression (filtered to samples matched in beta matrix)
# ===========================================================================
tpm_file <- file.path(project_root, "expression/gene_expression_matrix_protein_coding.tsv")
cat("Loading TPM matrix...\n")
tpm <- fread(tpm_file, header = TRUE, check.names = FALSE)
setnames(tpm, 1L, "sample_id")
tpm[, cancer   := sub("^TCGA-([A-Z]+)-.*", "\\1", sample_id)]
tpm[, st       := fifelse(grepl("-01A_1", sample_id), "Tumor",
                  fifelse(grepl("-11A_1", sample_id), "Healthy", "other"))]
tpm[, patient  := sub("^TCGA-([A-Z]+-TCGA-[A-Z0-9]+-[A-Z0-9]+)-[0-9]+[A-Z]_.*",
                       "\\1", sample_id)]
tpm[, pair_key := paste(patient, st, sep = "__")]
tpm <- tpm[cancer == cancer_type & st != "other" & pair_key %in% beta_samp_dt$pair_key]
cat("TCGA expression samples kept:", nrow(tpm),
    "(Tumor:", sum(tpm$st == "Tumor"),
    "| Healthy:", sum(tpm$st == "Healthy"), ")\n")

meta_cols <- c("sample_id", "cancer", "st", "patient", "pair_key")
gene_cols <- intersect(sig_genes, setdiff(names(tpm), meta_cols))
cat("Significant genes found in TCGA matrix:", length(gene_cols), "\n")

tpm_long <- melt(tpm[, c(meta_cols, gene_cols), with = FALSE],
                 id.vars       = meta_cols,
                 variable.name = "gene",
                 value.name    = "TPM")
tpm_long[, group  := st]
tpm_long[, source := "TCGA"]

# ===========================================================================
# Step 5 - Cell-line expression
# ===========================================================================
cl_expr_dir <- file.path(project_root, "cell_lines", cancer_type, "expression")
expr_files  <- list.files(cl_expr_dir, pattern = "^RNA_.*_genes_count\\.txt$",
                          full.names = TRUE)
if (length(expr_files) == 0) stop("No expression files in: ", cl_expr_dir)
if (!is.null(cl_type_map)) {
  pfx_pattern <- paste0("^(", paste(cl_type_map$expr_prefix, collapse = "|"), ")_")
  expr_files  <- expr_files[grepl(pfx_pattern, basename(expr_files))]
}
cat("Expression files matched:", length(expr_files), "\n")

expr_file_dt <- data.table(path = expr_files)
expr_file_dt[, basename_ := basename(path)]
if (!is.null(cl_type_map)) {
  expr_file_dt[, cell_line := sapply(basename_, function(f)
    cl_type_map$cell_line[which(startsWith(f, paste0(cl_type_map$expr_prefix, "_")))[1]])]
} else {
  expr_file_dt[, cell_line := sub("^RNA_(.+)_[0-9]+_genes_count\\.txt$", "\\1", basename_)]
}
expr_file_dt[, replicate := as.integer(sub("^RNA_.+_([0-9]+)_genes_count\\.txt$",
                                            "\\1", basename_))]
cat("Cell-line expression files:\n"); print(expr_file_dt[, .(cell_line, replicate, basename_)])

cl_expr <- rbindlist(lapply(seq_len(nrow(expr_file_dt)), function(i) {
  dt <- fread(expr_file_dt$path[i], header = TRUE, select = c("geneName", "TPM"))
  dt <- dt[geneName %in% sig_genes]
  dt[, cell_line := expr_file_dt$cell_line[i]]
  dt[, replicate := expr_file_dt$replicate[i]]
  dt
}))
setnames(cl_expr, "geneName", "gene")
cl_expr[, group  := cell_line]
cl_expr[, source := "CellLine"]
cat("Cell-line expression rows:", nrow(cl_expr), "\n")

# ===========================================================================
# Step 6 - Concordance (only when both healthy AND tumor cell lines given)
# ===========================================================================
meth_conc_dt    <- NULL
expr_conc_dt    <- NULL
meth_conc_label <- NULL
expr_conc_label <- NULL
both_label      <- NULL

if (!is.null(cl_type_map) && all(c("healthy", "tumor") %in% cl_type_map$type)) {

  # ---- Methylation concordance ----
  meth_typed <- merge(
    cl_meth[probe %in% sig_probes, .(probe, meth, group, source)],
    cl_type_map[, .(group = cell_line, cl_type = type)],
    by = "group"
  )
  cl_meth_med  <- meth_typed[, .(med = median(meth, na.rm = TRUE)), by = .(probe, cl_type)]
  cl_meth_wide <- dcast(cl_meth_med, probe ~ cl_type, value.var = "med")
  if (!all(c("healthy", "tumor") %in% names(cl_meth_wide))) {
    cat("Skipping meth concordance: no data for",
        paste(setdiff(c("healthy", "tumor"), names(cl_meth_wide)), collapse = ", "),
        "cell line type\n")
  } else {
    cl_meth_wide[, cl_dir := fifelse(tumor > healthy, "Hyper",
                              fifelse(tumor < healthy, "Hypo", "Flat"))]
    cl_meth_wide <- merge(cl_meth_wide, pair_dir[, .(probe, direction)],
                          by = "probe", all.x = TRUE)
    cl_meth_wide[, tcga_dir := fifelse(grepl("Hyper", direction), "Hyper",
                               fifelse(grepl("Hypo",  direction), "Hypo", NA_character_))]
    cl_meth_wide[, concordant := cl_dir == tcga_dir]
    covered_probes <- unique(cl_meth$probe)
    cl_meth_wide[, covered := probe %in% covered_probes]

    n_conc_m <- sum(cl_meth_wide$concordant, na.rm = TRUE)
    n_tot_m  <- sum(!is.na(cl_meth_wide$concordant))
    n_uncov  <- sum(!sig_probes %in% covered_probes)
    meth_conc_label <- paste0(
      "Meth: ", n_conc_m, "/", n_tot_m,
      " (", round(100 * n_conc_m / n_tot_m, 1), "%) concordant with TCGA",
      if (n_uncov > 0) paste0("  |  ", n_uncov, " probe(s) not covered by WGBS") else ""
    )
    meth_conc_dt <- cl_meth_wide[, .(probe, meth_concordant = concordant, covered)]
    cat("Meth concordance:", meth_conc_label, "\n")
  } # end meth concordance else

  # ---- Expression concordance ----
  expr_typed <- merge(
    cl_expr[gene %in% sig_genes, .(gene, TPM, group, source)],
    cl_type_map[, .(group = cell_line, cl_type = type)],
    by = "group"
  )
  cl_expr_med  <- expr_typed[, .(med = median(TPM, na.rm = TRUE)), by = .(gene, cl_type)]
  cl_expr_wide <- dcast(cl_expr_med, gene ~ cl_type, value.var = "med")
  if (!all(c("healthy", "tumor") %in% names(cl_expr_wide))) {
    cat("Skipping expr concordance: no data for",
        paste(setdiff(c("healthy", "tumor"), names(cl_expr_wide)), collapse = ", "),
        "cell line type\n")
  } else {
    cl_expr_wide[, cl_dir := fifelse(tumor > healthy, "Up",
                              fifelse(tumor < healthy, "Down", "Flat"))]
    gene_dir_uniq <- pair_dir[!duplicated(gene), .(gene, direction)]
    cl_expr_wide  <- merge(cl_expr_wide, gene_dir_uniq, by = "gene", all.x = TRUE)
    cl_expr_wide[, tcga_dir := fifelse(grepl("Down", direction), "Down",
                               fifelse(grepl("Up",   direction), "Up", NA_character_))]
    cl_expr_wide[, expr_concordant := cl_dir == tcga_dir]

    n_conc_e <- sum(cl_expr_wide$expr_concordant, na.rm = TRUE)
    n_tot_e  <- sum(!is.na(cl_expr_wide$expr_concordant))
    expr_conc_label <- paste0(
      "Expr: ", n_conc_e, "/", n_tot_e,
      " (", round(100 * n_conc_e / n_tot_e, 1), "%) concordant with TCGA"
    )
    expr_conc_dt <- cl_expr_wide[, .(gene, expr_concordant)]
    cat("Expr concordance:", expr_conc_label, "\n")
  } # end expr concordance else

  # ---- Both concordant ----
  if (!is.null(meth_conc_dt) && !is.null(expr_conc_dt)) {
    both_conc <- merge(
      pair_dir[, .(probe, gene)],
      meth_conc_dt[, .(probe, meth_concordant, covered)],
      by = "probe", all.x = TRUE
    )
    both_conc <- merge(
      both_conc,
      expr_conc_dt[, .(gene, expr_concordant)],
      by = "gene", all.x = TRUE
    )
    n_both <- sum(
      both_conc$meth_concordant == TRUE &
      both_conc$expr_concordant == TRUE &
      both_conc$covered         == TRUE,
      na.rm = TRUE
    )
    n_tot_both <- sum(
      both_conc$covered == TRUE,
      na.rm = TRUE
    )
    n_both_genes <- uniqueN(both_conc[
      meth_concordant == TRUE & expr_concordant == TRUE &
      covered == TRUE, gene])
    n_tot_genes <- uniqueN(both_conc[
      covered == TRUE, gene])
    both_label <- paste0(
      "Both concordant: ", n_both, "/", n_tot_both,
      " (", round(100 * n_both / n_tot_both, 1), "%) probe-gene pairs",
      " [", n_both_genes, "/", n_tot_genes, " unique genes]"
    )
        cat("Both concordant:", both_label, "\n")

    # Save concordant genes list
    concordant_genes <- both_conc[
      meth_concordant == TRUE & expr_concordant == TRUE & covered == TRUE,
      .(gene, probe)
    ]
    conc_out <- file.path(project_root, "cell_lines", cancer_type, "results",
                          paste0(cancer_type, "_concordant_genes_HCC1954_HMEC.tsv"))
    dir.create(dirname(conc_out), recursive = TRUE, showWarnings = FALSE)
    fwrite(concordant_genes, conc_out, sep = "\t")
    cat("Concordant genes saved to:", conc_out, "\n")

  } # end both concordant
} # end concordance block


# ===========================================================================
# Step 7 - Build unified panel data (one panel per probe-gene, dual y-axis)
# ===========================================================================
group_levels <- c("Healthy", "Tumor", cl_names)
n_groups     <- length(group_levels)

meth_xpos <- setNames(seq_len(n_groups),                 group_levels)
expr_xpos <- setNames(seq_len(n_groups) + n_groups + 1L, group_levels)
vline_at  <- n_groups + 1.5

# ---- Meth data ----
meth_long_all <- rbindlist(list(
  beta_long[probe %in% sig_probes, .(probe, y_raw = meth, group, source)],
  cl_meth[probe %in% sig_probes,   .(probe, y_raw = meth, group, source)]
))
meth_long_all <- merge(meth_long_all, pair_dir[, .(probe, gene, direction, pair_rank)],
                       by = "probe", all.x = TRUE)
if (!is.null(meth_conc_dt)) {
  meth_long_all <- merge(meth_long_all, meth_conc_dt, by = "probe", all.x = TRUE)
} else {
  meth_long_all[, meth_concordant := NA]
  meth_long_all[, covered         := NA]
}
meth_long_all[, meth_tag := fifelse(!is.na(covered)         & !covered,        "[0]",
                            fifelse(!is.na(meth_concordant) & meth_concordant,  "[+]",
                            fifelse(!is.na(meth_concordant) & !meth_concordant, "[-]", "")))]
meth_long_all[, block    := "Methylation"]
meth_long_all[, x_pos    := meth_xpos[as.character(group)]]
meth_long_all[, y_scaled := y_raw]  # meth is already 0-100 %

# ---- Expr data ----
expr_long_all <- rbindlist(list(
  tpm_long[gene %in% sig_genes, .(gene, y_raw = log2(TPM + 1), group, source)],
  cl_expr[gene %in% sig_genes,  .(gene, y_raw = log2(TPM + 1), group, source)]
))
expr_long_all <- merge(expr_long_all,
                       pair_dir[, .(gene, probe, direction, pair_rank)],
                       by = "gene", allow.cartesian = TRUE)
if (!is.null(expr_conc_dt)) {
  expr_long_all <- merge(expr_long_all, expr_conc_dt, by = "gene", all.x = TRUE)
} else {
  expr_long_all[, expr_concordant := NA]
}
expr_long_all[, expr_tag := fifelse(!is.na(expr_concordant) & expr_concordant,  "[+]",
                            fifelse(!is.na(expr_concordant) & !expr_concordant, "[-]", ""))]
expr_long_all[, block := "Expression"]
expr_long_all[, x_pos := expr_xpos[as.character(group)]]

expr_global_max <- max(expr_long_all$y_raw, na.rm = TRUE)
expr_long_all[, y_scaled := y_raw / expr_global_max * 100]
cat("Expression scaling: max log2(TPM+1) =", round(expr_global_max, 2),
    "-> mapped to 100 on the shared y-axis\n")

# ---- Combine ----
meth_cols <- c("probe", "gene", "pair_rank", "direction", "block",
               "y_raw", "y_scaled", "x_pos", "group", "source", "meth_tag")
expr_cols <- c("probe", "gene", "pair_rank", "direction", "block",
               "y_raw", "y_scaled", "x_pos", "group", "source", "expr_tag")

meth_sub <- meth_long_all[, ..meth_cols]
meth_sub[, expr_tag := ""]
expr_sub <- expr_long_all[, ..expr_cols]
expr_sub[, meth_tag := ""]

all_panels <- rbindlist(list(meth_sub, expr_sub), use.names = TRUE, fill = TRUE)

# Build per-probe-gene concordance tag line
tag_lookup <- unique(pair_dir[, .(probe, gene, direction, pair_rank)])
if (!is.null(meth_conc_dt)) {
  tag_lookup <- merge(tag_lookup, meth_conc_dt, by = "probe", all.x = TRUE)
  tag_lookup[, m_tag := fifelse(!is.na(covered) & !covered, "meth[0]",
                        fifelse(!is.na(meth_concordant) & meth_concordant,  "meth[+]",
                        fifelse(!is.na(meth_concordant) & !meth_concordant, "meth[-]", "")))]
} else {
  tag_lookup[, m_tag := ""]
}
if (!is.null(expr_conc_dt)) {
  tag_lookup <- merge(tag_lookup, expr_conc_dt, by = "gene", all.x = TRUE)
  tag_lookup[, e_tag := fifelse(!is.na(expr_concordant) & expr_concordant,  "expr[+]",
                        fifelse(!is.na(expr_concordant) & !expr_concordant, "expr[-]", ""))]
} else {
  tag_lookup[, e_tag := ""]
}
tag_lookup[, combined_tag := trimws(paste(m_tag, e_tag))]
tag_lookup[, panel_label  := paste0(probe, " / ", gene,
                                    "\n(", direction, ")",
                                    fifelse(nchar(combined_tag) > 0,
                                            paste0("  ", combined_tag), ""))]

all_panels <- merge(all_panels,
                    tag_lookup[, .(probe, gene, pair_rank, panel_label)],
                    by = c("probe", "gene", "pair_rank"), all.x = TRUE)

panel_lvls <- unique(all_panels[order(pair_rank), panel_label])
all_panels[, panel_label := factor(panel_label, levels = panel_lvls)]
all_panels[, group       := factor(group, levels = group_levels)]

tcga_data <- all_panels[source == "TCGA"]
cl_data   <- all_panels[source == "CellLine"]

cat("Final panels:", nlevels(all_panels$panel_label), "(one per probe-gene pair)\n")
cat("Groups:", paste(group_levels, collapse = " | "), "\n")

# ===========================================================================
# Step 8 - Colours
# ===========================================================================
fallback_palette <- c("#762A83", "#01665E", "#8C510A", "#4393C3",
                      "#F46D43", "#74ADD1", "#A50026", "#313695")
fallback_i <- 1L
cl_colours <- sapply(cl_names, function(cl) {
  if (!is.null(cl_type_map) && cl %in% cl_type_map$cell_line) {
    type <- cl_type_map[cell_line == cl, type]
    if (type == "healthy") "#1B7837" else "#D6604D"
  } else {
    col <- fallback_palette[fallback_i]; fallback_i <<- fallback_i + 1L; col
  }
})
group_colours <- c(Healthy = "#1B7837", Tumor = "#D6604D", cl_colours)

# ===========================================================================
# Step 9 - Subtitle
# ===========================================================================
subtitle <- NULL
if (!is.null(meth_conc_label) || !is.null(expr_conc_label)) {
  subtitle <- paste(c(meth_conc_label, expr_conc_label, both_label), collapse = "\n")
}

# ===========================================================================
# Step 10 - Plot
# ===========================================================================
expr_breaks_actual <- pretty(c(0, expr_global_max), n = 5)
expr_breaks_actual <- expr_breaks_actual[expr_breaks_actual <= expr_global_max]
expr_breaks_scaled <- expr_breaks_actual / expr_global_max * 100

x_breaks <- c(unname(meth_xpos), unname(expr_xpos))
x_labels <- rep(group_levels, 2)

meth_centre <- mean(meth_xpos)
expr_centre <- mean(expr_xpos)

p <- ggplot(all_panels, aes(x = x_pos, y = y_scaled, fill = group, colour = group)) +

  annotate("rect",
           xmin = min(meth_xpos) - 0.5, xmax = max(meth_xpos) + 0.5,
           ymin = -Inf, ymax = Inf, fill = "#F2EBF5", alpha = 0.5) +
  annotate("rect",
           xmin = min(expr_xpos) - 0.5, xmax = max(expr_xpos) + 0.5,
           ymin = -Inf, ymax = Inf, fill = "#EBF2F5", alpha = 0.5) +

  annotate("text", x = meth_centre, y = 97,
           label = "Methylation (%)", size = 2.6,
           colour = "grey25", fontface = "italic", hjust = 0.5) +
  annotate("text", x = expr_centre, y = 97,
           label = "Expression (log2 TPM+1)", size = 2.6,
           colour = "grey25", fontface = "italic", hjust = 0.5) +

  geom_jitter(data    = tcga_data,
              aes(group = x_pos),
              width   = 0.18, size = 0.25, alpha = 0.12, shape = 16) +
  geom_boxplot(data          = tcga_data,
               aes(group     = x_pos),
               outlier.shape = NA, fill = NA,
               width         = 0.55, linewidth = 0.55) +

  geom_point(data = cl_data, size = 2.5, alpha = 0.9) +
  stat_summary(data     = cl_data, fun = median,
               geom     = "crossbar", width = 0.4, linewidth = 0.5) +

  geom_vline(xintercept = vline_at,
             linetype = "dashed", colour = "grey40", linewidth = 0.5) +

  scale_fill_manual(values   = group_colours) +
  scale_colour_manual(values = group_colours) +

  scale_x_continuous(breaks = x_breaks, labels = x_labels) +
  scale_y_continuous(
    name     = "Methylation (%)",
    limits   = c(0, 100),
    breaks   = seq(0, 100, 25),
    sec.axis = sec_axis(
      transform = ~ . * expr_global_max / 100,
      name      = expression(log[2](TPM + 1)),
      breaks    = expr_breaks_actual,
      labels    = round(expr_breaks_actual, 1)
    )
  ) +

  facet_wrap2(~ panel_label, ncol = 2, axes = "all") +

  labs(x        = NULL,
       title    = paste0(cancer_type,
                         " -- methylation (left) & expression (right): paired view"),
       subtitle = subtitle) +

  theme_bw(base_size = 10) +
  theme(
    legend.position    = "none",
    axis.text.x        = element_text(angle = 35, hjust = 1, size = 7.5),
    axis.title.y.left  = element_text(colour = "#6B3D7A"),
    axis.title.y.right = element_text(colour = "#2E6B8A"),
    strip.background   = element_rect(fill = "grey92"),
    strip.text         = element_text(face = "bold", size = 9),
    panel.grid.minor   = element_blank(),
    plot.subtitle      = element_text(size = 8, colour = "grey30")
  )

# ===========================================================================
# Step 11 - Save PDF
# ===========================================================================
n_pairs <- nrow(pair_dir)
out_dir <- file.path(project_root, "cell_lines", cancer_type, "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cl_tag  <- if (!is.null(cl_type_map)) paste(cl_type_map$cell_line, collapse = "_") else "all"
out_pdf <- file.path(out_dir,
  paste0(cancer_type, "_", cl_tag, "_meth_expr_combined_boxplots.pdf"))
panel_width <- (n_groups * 2 + 1) * 0.55 + 1.5
ggsave(out_pdf, plot = p,
       width  = 2 * panel_width,
       height = ceiling(n_pairs / 2) * 3.8,
       units  = "in", device = "pdf",
       limitsize = FALSE)
cat("Plot saved to:", out_pdf, "\n")
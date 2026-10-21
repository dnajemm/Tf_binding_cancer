suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
})

# ── Arguments ─────────────────────────────────────────────────────────────────
args    <- commandArgs(trailingOnly = TRUE)
tf_name <- args[1]
if (is.na(tf_name) || !tf_name %in% c("NRF1", "BANP"))
  stop("Usage: Rscript correlation_controle_probe_test.R <TF>   (NRF1 or BANP)")
cat(sprintf("[INFO] TF: %s\n", tf_name))

# ── Paths ─────────────────────────────────────────────────────────────────────
ranking_file <- paste0("./results/integration/gene_ranking_", tf_name, ".tsv")
pairs_file   <- "./results/integration/meth_expr_filtered_pairs.tsv"
motif_file   <- "./methylation/probe_gene_pairs_in_motifs_with_tss.tsv.gz"
beta_dir     <- paste0("./methylation/beta_matrices_", tf_name)
meth_dir     <- "./methylation/filtered_methylation"
counts_dir   <- "./expression/deseq2/counts_all"
shared_file  <- "./results/multi_omics/expr_meth_shared_samples.txt"
out_dir      <- "./results/integration/controle_probe"
ckpt_dir     <- file.path(out_dir, paste0("checkpoints_", tf_name))
dir.create(out_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(ckpt_dir, recursive = TRUE, showWarnings = FALSE)

ckpt <- function(name) file.path(ckpt_dir, paste0(name, ".rds"))
load_or_run <- function(name, expr) {
  f <- ckpt(name)
  if (file.exists(f)) {
    cat(sprintf("[SKIP] %s already done — loading from checkpoint\n", name))
    return(readRDS(f))
  }
  result <- expr
  saveRDS(result, f)
  cat(sprintf("[SAVE] %s checkpoint saved\n", name))
  result
}

# ══════════════════════════════════════════════════════════════════════════════
# STEP 1 — Load and join
# ══════════════════════════════════════════════════════════════════════════════
todo <- load_or_run("step1_todo", {
  cat("[STEP 1] Loading gene ranking...\n")
  ranking    <- fread(ranking_file)
  genes_keep <- unique(ranking$gene)
  cat(sprintf("  %d genes in ranking\n", length(genes_keep)))

  cat("[STEP 1] Loading filtered pairs...\n")
  pairs <- fread(pairs_file)
  pairs <- pairs[TF == tf_name & gene %in% genes_keep]
  cat(sprintf("  %d gene-cancer-probe rows after filter\n", nrow(pairs)))

  cat("[STEP 1] Loading motif instance IDs...\n")
  motifs <- fread(cmd = paste("zcat", motif_file))
  motifs <- motifs[grepl(paste0("^", tf_name), TF) & gene %in% genes_keep]
  motifs <- unique(motifs[, .(TF, probe, gene, motif_instance_id)])
  cat(sprintf("  %d probe-gene-motif rows\n", nrow(motifs)))

  cat("[STEP 1] Joining...\n")
  td <- merge(pairs, motifs, by = c("probe", "gene"), all.x = FALSE)
  td[, motif_chr    := sub(":.*", "", motif_instance_id)]
  td[, motif_start  := as.integer(sub("^[^:]+:([0-9]+):.*", "\\1", motif_instance_id))]
  td[, motif_end    := as.integer(sub("^[^:]+:[0-9]+:([0-9]+):.*", "\\1", motif_instance_id))]
  td[, motif_strand := sub(".*:", "", motif_instance_id)]
  td
})

cat(sprintf("[STEP 1] %d rows\n", nrow(todo)))
print(todo[1:min(5,.N), .(gene, cancer, probe, motif_chr, motif_start, motif_end, motif_strand)])

# ══════════════════════════════════════════════════════════════════════════════
# STEP 2 — Find nearby control probes
# ══════════════════════════════════════════════════════════════════════════════
controls <- load_or_run("step2_controls", {
  cat("[STEP 2] Loading probe BED file...\n")
  all_probes <- fread(
    "./methylation/annotated_methylation_data_probes_filtered.bed",
    header = FALSE, col.names = c("probe_chr","probe_start","probe_end","probe"))
  all_probes[, probe_start := as.integer(probe_start)]
  all_probes[, probe_end   := as.integer(probe_end)]
  cat(sprintf("  %d probes loaded\n", nrow(all_probes)))

  todo[, win_start := fifelse(motif_strand == "+", motif_start - 500L,  motif_start - 1500L)]
  todo[, win_end   := fifelse(motif_strand == "+", motif_end   + 1500L, motif_end   + 500L)]
  todo[, win_start := pmax(0L, win_start)]

  cat("[STEP 2] Finding nearby control probes...\n")
  ctrl_list <- lapply(seq_len(nrow(todo)), function(i) {
    row <- todo[i]
    nearby <- all_probes[
      probe_chr   == row$motif_chr &
      probe_start >= row$win_start &
      probe_end   <= row$win_end &
      probe       != row$probe &
      !(probe_start < row$motif_end &
        probe_end > row$motif_start)
    ]
    if (nrow(nearby) == 0) return(NULL)
    nearby[, `:=`(gene         = row$gene,
                  cancer       = row$cancer,
                  target_probe = row$probe,
                  motif_id     = row$motif_instance_id,
                  win_start    = row$win_start,
                  win_end      = row$win_end)]
    nearby[]
  })
  rbindlist(ctrl_list, fill = TRUE)
})

cat(sprintf("[STEP 2] %d control probe rows | %d unique control probes\n",
            nrow(controls), uniqueN(controls$probe)))
print(controls[, .N, by = .(gene, cancer, target_probe)][order(-N)][1:min(10,.N)])

# ══════════════════════════════════════════════════════════════════════════════
# STEP 3 — Load beta values
# Option B: sample_key uses FULL patient ID (CANCER-TCGA-PART1-PART2)
#           to match original pipeline format exactly
# ══════════════════════════════════════════════════════════════════════════════
beta_long <- load_or_run("step3_beta_long", {
  shared         <- fread(shared_file, header = FALSE)$V1  # full IDs: CANCER-TCGA-PART1-PART2
  cancers_needed <- unique(todo$cancer)
  control_probes <- unique(controls$probe)
  target_probes  <- unique(todo$probe)

  # 3a — target probes from beta matrices
  # Column names in beta matrices: CANCER-TCGA-PART1-PART2_tumor / _healthy (lowercase)
  cat("[STEP 3a] Loading target probe betas from beta matrices...\n")
  target_beta <- rbindlist(lapply(cancers_needed, function(cx) {
    bf <- list.files(beta_dir, pattern = paste0(cx, ".*\\.tsv\\.gz$"), full.names = TRUE)
    if (length(bf) == 0) { cat(sprintf("  [WARN] No beta matrix for %s\n", cx)); return(NULL) }
    mat <- fread(bf[1])
    if (!"probe" %in% colnames(mat)) setnames(mat, 1, "probe")
    mat <- mat[probe %in% target_probes]
    if (nrow(mat) == 0) return(NULL)
    long <- melt(mat, id.vars = "probe", variable.name = "col_name", value.name = "meth_beta")
    long[, col_name := as.character(col_name)]
    # Extract full patient (CANCER-TCGA-PART1-PART2) — strip trailing _tumor / _healthy
    long[, full_patient := sub("_(tumor|healthy)$", "", col_name)]
    long[, sample_type  := sub(".*_(tumor|healthy)$", "\\1", col_name)]   # lowercase
    long[, sample_key   := paste(full_patient, sample_type, sep = "_")]   # CANCER-TCGA-XX-YY_tumor
    long <- long[sample_type %in% c("tumor","healthy") & full_patient %in% shared]
    if (nrow(long) == 0) return(NULL)
    long[, cancer := cx]
    if (max(long$meth_beta, na.rm = TRUE) <= 1.5)
      long[, meth_pct := meth_beta * 100] else
      long[, meth_pct := meth_beta]
    long[, .(cancer, probe, sample_key, meth_pct)]
  }), fill = TRUE)
  cat(sprintf("  Target beta rows: %d | Unique probes: %d\n",
              nrow(target_beta), uniqueN(target_beta$probe)))

  # 3b — control probes from raw methylation files (10 workers)
  # sample_key format: CANCER-TCGA-PART1-PART2_tumor  (full, lowercase)
  raw_files <- list.files(meth_dir,
    pattern = "(01A_1|11A_1)_annotated_methylation_filtered\\.bed\\.gz$",
    full.names = TRUE)
  raw_meta <- data.table(file = raw_files, fname = basename(raw_files))
  raw_meta[, cancer       := sub("^HM450_TCGA-([A-Z]+)-.*", "\\1", fname)]
  raw_meta[, full_patient := sub("^HM450_TCGA-([A-Z]+-TCGA-[A-Z0-9]+-[A-Z0-9]+)-.*", "\\1", fname)]
  raw_meta[, sample_type  := fifelse(grepl("-01A_1_", fname), "tumor",
                             fifelse(grepl("-11A_1_", fname), "healthy", NA_character_))]
  raw_meta[, sample_key   := paste(full_patient, sample_type, sep = "_")]
  raw_meta <- raw_meta[cancer %in% cancers_needed & !is.na(sample_type) & full_patient %in% shared]
  cat(sprintf("[STEP 3b] %d raw files to read (10 workers)...\n", nrow(raw_meta)))

  read_one <- function(i) {
    dt <- tryCatch(
      fread(cmd = paste("zcat", shQuote(raw_meta$file[i])),
            header = FALSE, select = c(4L, 5L),
            col.names = c("probe", "meth_beta"), showProgress = FALSE),
      error = function(e) NULL)
    if (is.null(dt) || nrow(dt) == 0) return(NULL)
    dt[, meth_beta := suppressWarnings(as.numeric(meth_beta))]
    dt <- dt[probe %chin% control_probes & is.finite(meth_beta)]
    if (nrow(dt) == 0) return(NULL)
    dt[, `:=`(cancer     = raw_meta$cancer[i],
               sample_key = raw_meta$sample_key[i])]
    if (max(dt$meth_beta, na.rm = TRUE) <= 1.5)
      dt[, meth_pct := meth_beta * 100] else
      dt[, meth_pct := meth_beta]
    dt[, .(cancer, probe, sample_key, meth_pct)]
  }

  n_files    <- nrow(raw_meta)
  chunk_size <- 10L
  chunks     <- split(seq_len(n_files), ceiling(seq_len(n_files) / chunk_size))
  n_chunks   <- length(chunks)

  cl <- makeCluster(10L)
  clusterExport(cl, varlist = c("raw_meta", "control_probes"), envir = environment())
  clusterEvalQ(cl, library(data.table))

  res_list <- vector("list", n_chunks)
  for (k in seq_along(chunks)) {
    res_list[[k]] <- parLapply(cl, chunks[[k]], read_one)
    cat(sprintf("\r  Progress: %d / %d files (%.0f%%)",
                min(k * chunk_size, n_files), n_files, 100 * k / n_chunks))
    flush.console()
  }
  stopCluster(cl)
  cat("\n")

  control_beta <- rbindlist(unlist(res_list, recursive = FALSE), fill = TRUE)
  cat(sprintf("  Control beta rows: %d | Unique probes: %d\n",
              nrow(control_beta), uniqueN(control_beta$probe)))

  # 3c — combine (no probe_type here — added in Step 5)
  rbindlist(list(target_beta, control_beta), fill = TRUE)
})

cat(sprintf("[STEP 3] Total beta rows: %d\n", nrow(beta_long)))

# ══════════════════════════════════════════════════════════════════════════════
# STEP 4 — Load expression from raw STAR count files (Option B)
# Exactly matches the original pipeline: log2(count+1), col 3 of count files
# sample_key format: CANCER-TCGA-PART1-PART2_tumor  (full, lowercase)
# ══════════════════════════════════════════════════════════════════════════════
count_long <- load_or_run("step4_count_long", {
  genes_needed   <- unique(todo$gene)
  cancers_needed <- unique(todo$cancer)
  shared         <- fread(shared_file, header = FALSE)$V1   # full IDs

  cat(sprintf("[STEP 4] Reading count files — %d cancers, %d genes...\n",
              length(cancers_needed), length(genes_needed)))

  expr_all <- rbindlist(lapply(cancers_needed, function(cx) {
    cnt_files <- list.files(counts_dir,
      pattern = paste0("^RNA_TCGA-", cx, "-.*_genes_count\\.txt$"),
      full.names = TRUE)
    if (length(cnt_files) == 0) {
      cat(sprintf("  [WARN] No count files for %s\n", cx))
      return(NULL)
    }
    bn  <- basename(cnt_files)
    # Full patient: CANCER-TCGA-PART1-PART2
    pt  <- sub("^RNA_TCGA-([A-Z]+)-TCGA-([A-Z0-9]+-[A-Z0-9]+)-[0-9]+[A-Z]_.*",
               "\\1-TCGA-\\2", bn)
    st  <- fifelse(grepl("-01[A-Z]_", bn), "tumor",
           fifelse(grepl("-11[A-Z]_", bn), "healthy", NA_character_))
    fkey <- paste(pt, st, sep = "_")   # CANCER-TCGA-PART1-PART2_tumor

    keep      <- !is.na(st) & pt %in% shared
    cnt_files <- cnt_files[keep]
    fkey      <- fkey[keep]
    if (length(cnt_files) == 0) return(NULL)

    cat(sprintf("  %s: %d files\n", cx, length(cnt_files)))
    rbindlist(lapply(seq_along(cnt_files), function(i) {
      dt <- tryCatch(
        fread(cnt_files[i], header = FALSE, select = c(2L, 3L),
              col.names = c("gene", "count"), showProgress = FALSE),
        error = function(e) NULL)
      if (is.null(dt) || nrow(dt) == 0) return(NULL)
      dt[, count := suppressWarnings(as.numeric(count))]
      dt <- dt[gene %chin% genes_needed & is.finite(count)]
      if (nrow(dt) == 0) return(NULL)
      dt[, `:=`(cancer = cx, sample_key = fkey[i])]
      dt[]
    }), fill = TRUE)
  }), fill = TRUE)

  expr_all[, log2_expr := log2(count + 1)]
  expr_all[is.finite(log2_expr)]
})

cat(sprintf("[STEP 4] %d expression rows | %d samples | %d genes\n",
            nrow(count_long), uniqueN(count_long$sample_key), uniqueN(count_long$gene)))

# ══════════════════════════════════════════════════════════════════════════════
# STEP 5 — Correlations
# ══════════════════════════════════════════════════════════════════════════════
cor_results <- load_or_run("step5_cor_results", {
  cat("[STEP 5] Building probe x sample table...\n")

  safe_pearson <- function(x, y, min_n = 5L) {
    ok <- is.finite(x) & is.finite(y)
    x <- x[ok]; y <- y[ok]; n <- length(x)
    if (n < min_n || length(unique(x)) < 2 || length(unique(y)) < 2)
      return(list(n = n, r = NA_real_, p = NA_real_))
    ct <- tryCatch(cor.test(x, y, method = "pearson"), error = function(e) NULL)
    if (is.null(ct)) return(list(n = n, r = NA_real_, p = NA_real_))
    list(n = n, r = unname(ct$estimate), p = ct$p.value)
  }

  all_probe_rows <- rbind(
    todo[,     .(gene, cancer, probe, probe_type = "target")],
    controls[, .(gene, cancer, probe, probe_type = "control")]
  )

  # Join methylation
  all_probe_rows <- merge(all_probe_rows, beta_long,
                          by = c("cancer", "probe"), allow.cartesian = TRUE)
  cat(sprintf("  After beta join: %d rows\n", nrow(all_probe_rows)))

  # Join expression (gene + cancer + sample_key — all full-format, match exactly)
  all_probe_rows <- merge(all_probe_rows,
                          count_long[, .(gene, cancer, sample_key, log2_expr)],
                          by = c("gene", "cancer", "sample_key"),
                          all = FALSE)
  cat(sprintf("  After expression join: %d rows\n", nrow(all_probe_rows)))

  all_probe_rows[
    is.finite(meth_pct) & is.finite(log2_expr),
    { cr <- safe_pearson(meth_pct, log2_expr)
      .(n = cr$n, pearson_r = cr$r, pearson_p = cr$p) },
    by = .(gene, cancer, probe, probe_type)
  ]
})

cat(sprintf("[STEP 5] Correlations — target: %d | control: %d\n",
    cor_results[probe_type == "target"  & !is.na(pearson_r), .N],
    cor_results[probe_type == "control" & !is.na(pearson_r), .N]))
cat("[INFO] Preview (target):\n")
print(cor_results[probe_type == "target"][1:min(10,.N),
      .(gene, cancer, probe, probe_type, n, pearson_r, pearson_p)])
cat("[INFO] Preview (control):\n")
print(cor_results[probe_type == "control"][1:min(10,.N),
      .(gene, cancer, probe, probe_type, n, pearson_r, pearson_p)])

# ══════════════════════════════════════════════════════════════════════════════
# STEP 6 — Save
# ══════════════════════════════════════════════════════════════════════════════
out_file <- file.path(out_dir, paste0("control_probe_correlations_", tf_name, ".tsv"))
fwrite(cor_results, out_file, sep = "\t")
cat(sprintf("[DONE] Results saved to %s\n", out_file))

# ══════════════════════════════════════════════════════════════════════════════
# STEP 7 — Compare target r vs mean control r
# For each target probe: mean(control pearson_r) vs target pearson_r
# delta_r = target_r - mean_ctrl_r  (negative = target more correlated than controls)
# ══════════════════════════════════════════════════════════════════════════════
comparison <- load_or_run("step7_comparison", {
  cat("[STEP 7] Comparing target r vs mean control r...\n")
 
  # Link each control probe back to its target probe
  ctrl_map <- unique(controls[, .(gene, cancer, probe, target_probe)])
 
  ctrl_r <- merge(
    cor_results[probe_type == "control" & !is.na(pearson_r),
                .(gene, cancer, probe, ctrl_r = pearson_r)],
    ctrl_map,
    by = c("gene", "cancer", "probe")
  )
 
  # Mean control r per (gene, cancer, target_probe)
  ctrl_mean <- ctrl_r[, .(mean_ctrl_r = mean(ctrl_r), n_ctrl = .N),
                       by = .(gene, cancer, target_probe)]
 
  # Target r
  target_r <- cor_results[probe_type == "target" & !is.na(pearson_r),
                           .(gene, cancer, probe, target_r = pearson_r)]
 
  comp <- merge(target_r, ctrl_mean,
                by.x = c("gene", "cancer", "probe"),
                by.y = c("gene", "cancer", "target_probe"))
 
  comp[, delta_r := target_r - mean_ctrl_r]
  comp
})
 
cat("\n=== Target r vs mean control r ===\n")
print(comparison[order(delta_r), .(gene, cancer, probe, target_r, mean_ctrl_r, n_ctrl, delta_r)])
 
cat("\n=== Summary of delta_r ===\n")
print(summary(comparison$delta_r))
 
cat(sprintf("\nTarget probes more negative than controls (delta_r < 0): %d / %d\n",
    comparison[delta_r < 0, .N], nrow(comparison)))
 
comp_file <- file.path(out_dir, paste0("target_vs_control_mean_", tf_name, ".tsv"))
fwrite(comparison, comp_file, sep = "\t")
cat(sprintf("[STEP 7] Comparison saved to %s\n", comp_file))
 
# ══════════════════════════════════════════════════════════════════════════════
# STEP 8 — Scatter plot: target r vs mean control r
# ══════════════════════════════════════════════════════════════════════════════
cat("[STEP 8] Generating scatter plot...\n")
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
 
thresh <- 0.1   # 10% r change threshold
 
lim <- range(c(comparison$target_r, comparison$mean_ctrl_r), na.rm = TRUE)
pad <- diff(lim) * 0.07
lim <- c(lim[1] - pad, lim[2] + pad)
 
n_below   <- comparison[delta_r <  0,      .N]
n_above   <- comparison[delta_r >= 0,      .N]
n_sig_neg <- comparison[delta_r < -thresh, .N]
n_sig_pos <- comparison[delta_r >  thresh, .N]
 
# Label only top 10 by absolute delta_r
top10_idx <- order(abs(comparison$delta_r), decreasing = TRUE)[1:min(10L, nrow(comparison))]
comparison[, label := fifelse(.I %in% top10_idx,
                               paste(gene, cancer, sep = "\n"), NA_character_)]
 
p <- ggplot(comparison, aes(x = mean_ctrl_r, y = target_r, label = label)) +
  # zero lines
  geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.4) +
  geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.4) +
  # threshold bands (±10% r change)
  geom_abline(slope = 1, intercept =  thresh, linetype = "dotted",
              colour = "#4575b4", linewidth = 0.7) +
  geom_abline(slope = 1, intercept = -thresh, linetype = "dotted",
              colour = "#d73027", linewidth = 0.7) +
  # identity diagonal
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey30", linewidth = 0.6) +
  # points
  geom_point(aes(colour = delta_r), size = 1, alpha = 0.85) +
  # labels — only top 10 by |delta_r|, plain text no box
  geom_text_repel(na.rm = TRUE, size = 2.5, max.overlaps = 20,
                  box.padding = 0.35, point.padding = 0.2,
                  segment.colour = "grey40", segment.size = 0.3) +
  scale_colour_gradient2(
    low = "#d73027", mid = "grey85", high = "#4575b4",
    midpoint = 0, limits = c(-0.5, 0.5), oob = scales::squish,
    name = "delta_r\n(target − ctrl)") +
  coord_equal(xlim = lim, ylim = lim) +
  annotate("text", x = lim[1] + 0.02, y = lim[2] - 0.02,
           label = sprintf(">10%% above: %d", n_sig_pos),
           hjust = 0, vjust = 1, size = 3.5, colour = "#4575b4") +
  annotate("text", x = lim[2] - 0.02, y = lim[1] + 0.02,
           label = sprintf(">10%% below: %d", n_sig_neg),
           hjust = 1, vjust = 0, size = 3.5, colour = "#d73027") +
  labs(
    title    = sprintf("%s — target probe r vs mean control probe r", tf_name),
    subtitle = sprintf("%d / %d probes below diagonal  |  ±10%% threshold: %d negative, %d positive",
                       n_below, nrow(comparison), n_sig_neg, n_sig_pos),
    x = "Mean control probe r",
    y = "Target probe r") +
  theme_bw(base_size = 13) +
  theme(legend.position = "right", plot.subtitle = element_text(size = 10))
 
plot_file <- file.path(out_dir, paste0("plot_target_vs_ctrl_mean_", tf_name, ".png"))
ggsave(plot_file, p, width = 6, height = 6, dpi = 150)
cat(sprintf("[DONE] Plot saved to %s\n", plot_file))

# ══════════════════════════════════════════════════════════════════════════════
# STEP 9 — Fisher z-test: is target r significantly more negative than controls?
#
# For each target probe:
#   z_target    = atanh(target_r),  SE = 1/sqrt(n_target - 3)
#   z_ctrl_i    = atanh(ctrl_r_i)  for each control probe
#   z_mean_ctrl = mean(z_ctrl_i),   SE = sd(z_ctrl_i) / sqrt(n_ctrl_probes)
#   z_stat      = (z_target - z_mean_ctrl) / sqrt(SE_target^2 + SE_mean_ctrl^2)
#   p_fisher    = pnorm(z_stat)  [one-tailed: target MORE negative than controls]
# ══════════════════════════════════════════════════════════════════════════════
fisher_results <- load_or_run("step9_fisher", {
  cat("[STEP 9] Fisher z-test: target r vs control r distribution...\n")
 
  # Target probe info (r + n)
  target_info <- cor_results[probe_type == "target" & !is.na(pearson_r),
                              .(gene, cancer, probe,
                                target_r = pearson_r, n_target = n)]
  target_info[, z_target  := atanh(target_r)]
  target_info[, se_target := 1 / sqrt(pmax(n_target - 3L, 1L))]
 
  # Control probe r values linked to each target probe
  ctrl_map  <- unique(controls[, .(gene, cancer, probe, target_probe)])
  ctrl_info <- merge(
    cor_results[probe_type == "control" & !is.na(pearson_r),
                .(gene, cancer, probe, ctrl_r = pearson_r)],
    ctrl_map,
    by = c("gene", "cancer", "probe")
  )
  ctrl_info[, z_ctrl := atanh(ctrl_r)]
 
  # Summarise control z distribution per target probe
  ctrl_summary <- ctrl_info[, {
    n_p   <- .N
    z_mu  <- mean(z_ctrl)
    z_sd  <- if (n_p > 1L) sd(z_ctrl) else NA_real_
    se_mu <- if (n_p > 1L) z_sd / sqrt(n_p) else NA_real_
    .(n_ctrl_probes = n_p,
      z_mean_ctrl   = z_mu,
      z_sd_ctrl     = z_sd,
      se_mean_ctrl  = se_mu)
  }, by = .(gene, cancer, target_probe)]
 
  # Merge and compute test statistic
  res <- merge(target_info, ctrl_summary,
               by.x = c("gene", "cancer", "probe"),
               by.y = c("gene", "cancer", "target_probe"))
 
  # z_stat: NA when only 1 control probe (no SD estimate)
  res[, z_stat   := (z_target - z_mean_ctrl) /
                     sqrt(se_target^2 + se_mean_ctrl^2)]
  # One-tailed: P(target more negative than controls)
  res[, p_fisher := pnorm(z_stat)]
  res[, p_adj    := p.adjust(p_fisher, method = "BH")]
 
  res[order(p_fisher)]
})
 
cat("\n=== Fisher z-test results ===\n")
print(fisher_results[, .(gene, cancer, probe, target_r, n_ctrl_probes,
                          z_stat, p_fisher, p_adj)])
 
cat(sprintf("\nSignificant (p_fisher < 0.05): %d / %d\n",
    fisher_results[p_fisher < 0.05, .N], nrow(fisher_results)))
cat(sprintf("Significant after BH correction (p_adj < 0.05): %d / %d\n",
    fisher_results[p_adj < 0.05, .N], nrow(fisher_results)))
 
fisher_file <- file.path(out_dir, paste0("fisher_test_", tf_name, ".tsv"))
fwrite(fisher_results, fisher_file, sep = "\t")
cat(sprintf("[DONE] Fisher results saved to %s\n", fisher_file))

# ══════════════════════════════════════════════════════════════════════════════
# STEP 10 — Filter gene ranking by Fisher p_adj
# Keep only gene-cancer pairs where the target probe correlation is
# significantly more negative than its control probes (p_adj < 0.05)
# ══════════════════════════════════════════════════════════════════════════════
cat("[STEP 10] Filtering gene ranking by Fisher p_adj < 0.05...\n")
 
ranking <- fread(ranking_file)
cat(sprintf("  Original ranking: %d genes\n", nrow(ranking)))
 
# Aggregate ALL Fisher results to gene level (for annotated file)
best_p_all <- fisher_results[, .(best_p_adj       = min(p_adj, na.rm = TRUE),
                                  best_cancer      = cancer[which.min(p_adj)],
                                  best_probe       = probe[which.min(p_adj)],
                                  n_cancers_tested = uniqueN(cancer)),
                              by = gene]
 
# Preserve original ranking order after merge
ranking[, .rank_order := .I]
ranking_annotated <- merge(ranking, best_p_all, by = "gene", all.x = TRUE)
setorder(ranking_annotated, .rank_order)
ranking_annotated[, .rank_order := NULL]
 
cat(sprintf("  Genes with Fisher result:    %d\n", sum(!is.na(ranking_annotated$best_p_adj))))
cat(sprintf("  After filter (p_adj < 0.05): %d genes\n",
            ranking_annotated[best_p_adj < 0.05, .N]))
 
# Load source files to recalculate summary columns for significant cancers only
pairs_all  <- fread(pairs_file)   # has mean_delta_beta, log2foldchange, adj.P.Val, padj, pearson_r, pearson_p
stats_file <- paste0("./results/integration/meth_expr_wilcoxon_", tf_name, "_stats.tsv")
stats_all  <- fread(stats_file)   # has expr_wilcox_p, meth_wilcox_p
 
# Significant gene-cancer pairs from Fisher test
sig_pairs <- fisher_results[p_adj < 0.05, .(gene, cancer)]
 
# Filter both source files to significant cancers only
pairs_sig <- merge(pairs_all[TF == tf_name], sig_pairs, by = c("gene", "cancer"))
stats_sig <- merge(stats_all[TF == tf_name], sig_pairs, by = c("gene", "cancer"))
 
# Recalculate summary columns per gene from significant cancers
pairs_summary <- pairs_sig[, .(
  n_cancers      = uniqueN(cancer),
  cancers        = paste(sort(unique(cancer)), collapse = ", "),
  mean_pearson_r = mean(pearson_r,       na.rm = TRUE),
  mean_delta_beta= mean(mean_delta_beta, na.rm = TRUE),
  mean_log2FC    = mean(log2foldchange,  na.rm = TRUE),
  best_pearson_p = min(pearson_p,        na.rm = TRUE),
  best_padj      = min(padj,             na.rm = TRUE),
  best_adj_P_Val = min(adj.P.Val,        na.rm = TRUE)
), by = gene]
 
stats_summary <- stats_sig[, .(
  best_expr_wilcox_p = min(expr_wilcox_p, na.rm = TRUE),
  best_meth_wilcox_p = min(meth_wilcox_p, na.rm = TRUE)
), by = gene]
 
# Fisher-level summary (best cancer, probe, p_adj from Fisher test)
best_p_sig <- fisher_results[p_adj < 0.05,
                              .(best_p_adj  = min(p_adj, na.rm = TRUE),
                                best_cancer = cancer[which.min(p_adj)],
                                best_probe  = probe[which.min(p_adj)]),
                              by = gene]
 
# Drop columns being replaced, then merge all summaries
cols_replace <- c("n_cancers", "cancers", "mean_pearson_r", "mean_delta_beta",
                  "mean_log2FC", "best_pearson_p", "best_padj", "best_adj_P_Val",
                  "best_expr_wilcox_p", "best_meth_wilcox_p")
ranking[, .rank_order := .I]
ranking_filtered <- merge(ranking[, !cols_replace, with = FALSE], pairs_summary,  by = "gene")
ranking_filtered <- merge(ranking_filtered, stats_summary, by = "gene")
ranking_filtered <- merge(ranking_filtered, best_p_sig,   by = "gene")
setorder(ranking_filtered, .rank_order)
ranking_filtered[, .rank_order := NULL]
 
cat(sprintf("  Genes passing filter (p_adj < 0.05): %d\n", nrow(ranking_filtered)))
 
# Save outputs
filtered_file  <- file.path(out_dir, paste0("gene_ranking_", tf_name, "_ctrl_filtered.tsv"))
annotated_file <- file.path(out_dir, paste0("gene_ranking_", tf_name, "_ctrl_annotated.tsv"))
fwrite(ranking_filtered,  filtered_file,  sep = "\t")
fwrite(ranking_annotated, annotated_file, sep = "\t")
 
cat(sprintf("[STEP 10] Filtered ranking saved to  %s\n", filtered_file))
cat(sprintf("[STEP 10] Annotated ranking saved to %s\n", annotated_file))
 
cat("\n=== Top filtered hits ===\n")
print(ranking_filtered[order(best_p_adj)][1:min(20, .N),
      .(gene, best_cancer, best_probe, n_cancers, best_p_adj)])
##################################
# Sample inventory — run before anything else
##################################

##################################
# Create lists of expression and methylation files, and a list of shared samples to use the correlation analysis on them 
# only 11A and 01A and the first replicates _1 were chosen for the analysis
##################################

EXPR_LIST="./results/multi_omics/expr_files.txt"
METH_LIST="./results/multi_omics/meth_files.txt"
SHARED_LIST="./results/multi_omics/expr_meth_shared_samples.txt"

mkdir -p ./results/multi_omics

# 1. Expression files
echo "=== Expression files ==="
> "$EXPR_LIST"
for f in ./expression/RNA_TCGA-*-01A_1.augmented_star_gene_counts_filtered.tsv \
          ./expression/RNA_TCGA-*-11A_1.augmented_star_gene_counts_filtered.tsv; do
    [ -e "$f" ] || continue
    echo "  $f"
    echo "$f" >> "$EXPR_LIST"
done
sort -o "$EXPR_LIST" "$EXPR_LIST"

# 2. Methylation files
echo ""
echo "=== Methylation files ==="
> "$METH_LIST"
for f in ./methylation/filtered_methylation/*-01A_1_annotated_methylation_filtered.bed.gz \
          ./methylation/filtered_methylation/*-11A_1_annotated_methylation_filtered.bed.gz; do
    [ -e "$f" ] || continue
    echo "  $f"
    echo "$f" >> "$METH_LIST"
done
sort -o "$METH_LIST" "$METH_LIST"

# 3. Shared samples
awk -F'/' '{split($NF,a,"-"); print a[2]"-"a[3]"-"a[4]"-"a[5]}' "$EXPR_LIST" | sort -u > /tmp/expr_ids.txt
awk -F'/' '{split($NF,a,"-"); print a[2]"-"a[3]"-"a[4]"-"a[5]}' "$METH_LIST" | sort -u > /tmp/meth_ids.txt
comm -12 /tmp/expr_ids.txt /tmp/meth_ids.txt > "$SHARED_LIST"

echo ""
echo "=== Summary ==="
echo "Expression files : $(wc -l < $EXPR_LIST)"
echo "Methylation files: $(wc -l < $METH_LIST)"
echo "Shared samples   : $(wc -l < $SHARED_LIST)"

##################################
# Prepare Deseq2 files
##################################
############################################
# Create the design file for DESeq2 analysis only takes the 01A and 11A samples (tumor and healthy) and the first replicate
############################################
SHARED_LIST="./results/multi_omics/expr_meth_shared_samples.txt"

mkdir -p ./expression/deseq2
echo -e "sample\tcondition\treplicate" > ./expression/deseq2/design_file.txt

total=$(ls ./expression/RNA_TCGA-*-01A_1.augmented_star_gene_counts_filtered.tsv \
           ./expression/RNA_TCGA-*-11A_1.augmented_star_gene_counts_filtered.tsv \
           2>/dev/null | wc -l)
count=0
included=0

for f in ./expression/RNA_TCGA-*-01A_1.augmented_star_gene_counts_filtered.tsv \
          ./expression/RNA_TCGA-*-11A_1.augmented_star_gene_counts_filtered.tsv; do
    [ -e "$f" ] || continue
    count=$((count + 1))

    patient_id=$(basename "$f" | grep -oP 'TCGA-\K[A-Z]+-[A-Z0-9]+-[A-Z0-9]+-[A-Z0-9]+')
    if ! grep -q "^${patient_id}$" "$SHARED_LIST"; then
        echo "  Skipping  $count/$total: $(basename $f)  [no matching methylation]"
        continue
    fi

    included=$((included + 1))
    echo "  Including $count/$total: $(basename $f)"

    sample=$(basename "$f" | cut -d'.' -f1)
    condition=$(echo "$sample" | cut -d'-' -f6 | cut -d'_' -f1)
    replicate=$(echo "$sample" | cut -d'-' -f6 | cut -d'_' -f2)
    echo -e "$sample\t$condition\t$replicate" >> ./expression/deseq2/design_file.txt
done

echo "Done! $included/$total samples included."

awk 'BEGIN{OFS="\t"}
NR==1{print; next}
$2=="01A"{$2="tumor";   print}
$2=="11A"{$2="healthy"; print}' ./expression/deseq2/design_file.txt \
    > ./expression/deseq2/design_file_annotated.txt

#############################################
# Create one design file per cancer type
#############################################
for cancer in $(tail -n +2 ./expression/deseq2/design_file_annotated.txt | cut -d'-' -f2 | sort -u); do
    head -1 ./expression/deseq2/design_file_annotated.txt > ./expression/deseq2/design_${cancer}.txt
    grep "\-${cancer}\-" ./expression/deseq2/design_file_annotated.txt >> ./expression/deseq2/design_${cancer}.txt
    echo "Created: design_${cancer}.txt"
done

# Barplot the number of samples per cancer type

Rscript - <<'EOF'
library(ggplot2)
library(data.table)

beta_dir <- "./methylation/beta_matrices_BANP"
files <- list.files(beta_dir, pattern = "\\.tsv\\.gz$", full.names = TRUE)

dt <- rbindlist(lapply(files, function(f) {
  cols <- colnames(fread(f, nrows = 0))
  cols <- cols[cols != "probe"]
  cancer <- sub("_BANP_beta_matrix\\.tsv\\.gz$", "", basename(f))
  data.table(sample = cols, cancer = cancer)
}))

dt[, condition := ifelse(grepl("_tumor$", sample), "tumor", "healthy")]

counts <- dt[, .N, by = .(cancer, condition)]
totals <- counts[, .(total = sum(N)), by = cancer]
counts[, cancer := factor(cancer, levels = totals[order(-total), cancer])]

ggplot(counts, aes(x = cancer, y = N, fill = condition)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.9, preserve = "single")) +
  geom_text(aes(label = N),
          position = position_dodge(width = 0.9, preserve = "single"),
          hjust = 0, vjust = -0.5, size = 3)+
  scale_fill_manual(values = c(tumor = "#C0392B", healthy = "#2471A3")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Samples per cancer type (expression + methylation shared)",
       subtitle = paste0("Total: ", nrow(dt), " samples"),
       x = NULL, y = "Samples", fill = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA))
ggsave("./results/multi_omics/samples_2d_per_cancer.png",
       width = 10, height = 5, dpi = 150)
EOF


####################################
# Deseq2 analysis for all genes (not just NRF1 or BANP) for each cancer type USED A PVALUE OF 1e-10
####################################

# build set of valid patient IDs from column 1 only for samples that have methylation expression and mutation data 
mkdir -p ./expression/deseq2/counts_all
for f in ./expression/RNA_TCGA-*.augmented_star_gene_counts_filtered.tsv; do
    sample=$(basename "$f" | cut -d'.' -f1)
    if grep -q "^$sample" ./expression/deseq2/design_file_annotated.txt; then
        awk -F'\t' 'BEGIN{OFS="\t"} {print $1,$2,$4,1000}' "$f" \
            > ./expression/deseq2/counts_all/${sample}_genes_count.txt
        echo "Processed: $sample"
    fi
done

echo "Done. $(ls ./expression/deseq2/counts_all | wc -l) files created."

# Run the Deseq2 pipeline for each cancer type that has at least 10 tumor and 10 healthy sample (01A and 11A) chose the cancer type from the /TF_binding_cancer/results/multi_omics/sample_meth_expr_counts.txt
mkdir -p ./results/expression/deseq2/all_genes

# Get cancer types with >= 10 tumor AND >= 10 healthy samples in the design file
cancer_list=$(
    # count tumor per cancer
    tail -n +2 ./expression/deseq2/design_file_annotated.txt | \
        awk '$2=="tumor"{split($1,a,"-"); print a[2]}' | sort | uniq -c | \
        awk '$1>=10{print $2}' > /tmp/tumor_ok.txt

    # count healthy per cancer
    tail -n +2 ./expression/deseq2/design_file_annotated.txt | \
        awk '$2=="healthy"{split($1,a,"-"); print a[2]}' | sort | uniq -c | \
        awk '$1>=10{print $2}' > /tmp/healthy_ok.txt

    # keep only cancers in both
    comm -12 <(sort /tmp/tumor_ok.txt) <(sort /tmp/healthy_ok.txt)
)

echo "Cancer types passing filter:"
echo "$cancer_list"
batch=0
for cancer in $cancer_list; do
    Rscript --vanilla ./scripts/Deseq2.R \
        -n ${cancer} \
        -d ./expression/deseq2/design_${cancer}.txt \
        -c ./expression/deseq2/counts_all \
        -o ./results/expression/deseq2/all_genes \
        > ./results/expression/deseq2/all_genes/${cancer}.log 2>&1 &

    batch=$((batch+1))
    if [ $batch -eq 2 ]; then
        wait
        batch=0
    fi
done
wait
echo "All done!"

##################
# Combine methylation beta values for all samples into one file per cancer type and put them in a beta_matrices directory
##################

Rscript -e '
library(data.table)
library(parallel)

tf_files <- list(
  BANP = "./methylation/overlaps/intersected_motifs2mm_HM450/BANP_mm0to2_noCGmm_probes.txt",
  NRF1 = "./methylation/overlaps/intersected_motifs2mm_HM450/NRF1_mm0to2_noCGmm_probes.txt"
)
tf_probes <- lapply(tf_files, function(f) fread(f, header = FALSE)$V1)
all_probes <- unique(c(tf_probes$BANP, tf_probes$NRF1))  # union — read once

shared <- fread("./results/multi_omics/expr_meth_shared_samples.txt", header = FALSE)$V1

files <- list.files("./methylation/filtered_methylation",
                    pattern = "(01A_1|11A_1)_annotated_methylation_filtered[.]bed[.]gz$",
                    full.names = TRUE)

parse_filename <- function(f) {
  bn    <- basename(f)
  parts <- strsplit(bn, "-")[[1]]
  condition_raw <- strsplit(parts[6], "_")[[1]][1]
  list(cancer  = parts[2],
       patient = paste(parts[2], parts[3], parts[4], parts[5], sep = "-"),
       type    = ifelse(condition_raw == "01A", "tumor", "healthy"))
}

meta <- rbindlist(lapply(files, function(f) {
  p <- parse_filename(f)
  data.table(file = f, cancer = p$cancer, patient = p$patient, type = p$type)
}))

meta <- meta[patient %in% shared]

design <- fread("./expression/deseq2/design_file_annotated.txt")
design[, expr_patient := sub("^RNA_TCGA-", "", sample)]
design[, expr_patient := sub("-01A_1$|-11A_1$", "", expr_patient)]
allowed_pairs <- paste(design$expr_patient, design$condition, sep = "__")
meta[, pair := paste(patient, type, sep = "__")]
meta <- meta[pair %in% allowed_pairs][, pair := NULL]
message("Samples after filtering: ", nrow(meta))

# ── Read ALL files ONCE in parallel ──────────────────────────────────────────
n_cores <- min(detectCores() - 1, 8)
message("Reading ", nrow(meta), " files on ", n_cores, " cores...")

read_one <- function(i) {
  d <- fread(cmd = paste("zcat", shQuote(meta$file[i])),
             select = c(4, 5), col.names = c("probe", "beta"), header = FALSE)
  d[, beta := as.numeric(beta)]
  d <- d[!is.na(beta) & probe %chin% all_probes]   # %chin% = fast character %in%
  if (nrow(d) == 0L) return(NULL)
  d <- d[, .(beta = mean(beta)), by = probe]
  d[, `:=`(cancer = meta$cancer[i],
           sample = paste0(meta$patient[i], "_", meta$type[i]))]
  d
}

dt_all <- rbindlist(mclapply(seq_len(nrow(meta)), read_one, mc.cores = n_cores))
message("Done reading. Rows in memory: ", nrow(dt_all))

# ── Write matrices per TF × cancer (no more file I/O) ────────────────────────
for (tf in names(tf_probes)) {
  probes_keep <- tf_probes[[tf]]
  outdir <- paste0("./methylation/beta_matrices_", tf)
  dir.create(outdir, showWarnings = FALSE)
  message("\n=== TF: ", tf, " ===")

  dt_tf <- dt_all[probe %chin% probes_keep]

  for (cx in sort(unique(meta$cancer))) {
    outfile <- file.path(outdir, paste0(cx, "_", tf, "_beta_matrix.tsv.gz"))
    dt_cx <- dt_tf[cancer == cx]
    if (nrow(dt_cx) == 0L) { message("  Skipping ", cx); next }
    mat  <- dcast(dt_cx, probe ~ sample, value.var = "beta", fun.aggregate = mean)
    vals <- unlist(mat[, -1]); vals <- vals[!is.na(vals)]
    message("  ", cx, ": ", nrow(mat), " probes x ", ncol(mat)-1, " samples",
            " | beta ", round(min(vals),2), "–", round(max(vals),2))
    fwrite(mat, outfile, sep = "\t", na = "NA")
    rm(dt_cx, mat)
  }
  rm(dt_tf); gc()
}
message("All done.")
'

####################################
# Do a differential methylation analysis using limma
####################################
mkdir -p ./results/methylation/limma
Rscript - <<'EOF'
suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(ggplot2)
  library(ggrepel)
})

out_dir <- "./results/methylation/limma"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

probe_gene <- fread("./methylation/probe_gene_pairs_in_motifs.tsv")
setnames(probe_gene, tolower(colnames(probe_gene)))
probe_gene <- unique(probe_gene[, .(tf, probe, gene)])

tfs <- c("BANP", "NRF1")

run_limma <- function(mat, probes, tumor_cols, healthy_cols) {
  beta_mat <- as.matrix(mat[, c(tumor_cols, healthy_cols), with = FALSE])
  rownames(beta_mat) <- probes
  beta_clipped <- beta_mat
  beta_clipped[!is.na(beta_clipped) & beta_clipped < 1]  <- 1
  beta_clipped[!is.na(beta_clipped) & beta_clipped > 99] <- 99
  M_mat  <- log2(beta_clipped / (100 - beta_clipped))
  group  <- factor(c(rep("tumor", length(tumor_cols)), rep("healthy", length(healthy_cols))))
  design <- model.matrix(~0 + group)
  colnames(design) <- c("healthy", "tumor")
  fit  <- lmFit(M_mat, design)
  cont <- makeContrasts(tumor - healthy, levels = design)
  fit2 <- eBayes(contrasts.fit(fit, cont))
  res  <- as.data.table(topTable(fit2, number = Inf, sort.by = "none"))
  res[, probe := probes]
  res[, mean_delta_beta := rowMeans(beta_mat[, tumor_cols, drop = FALSE], na.rm = TRUE) -
                           rowMeans(beta_mat[, healthy_cols, drop = FALSE], na.rm = TRUE)]
  res
}

for (tf in tfs) {
  mat_dir   <- paste0("./methylation/beta_matrices_", tf)
  mat_files <- list.files(mat_dir, pattern = "\\.tsv\\.gz$", full.names = TRUE)
  cancers   <- sub(paste0("_", tf, "_beta_matrix\\.tsv\\.gz$"), "", basename(mat_files))

  all_res <- list()

  for (cancer in cancers) {
    mat_file <- file.path(mat_dir, paste0(cancer, "_", tf, "_beta_matrix.tsv.gz"))
    if (!file.exists(mat_file)) next

    cat("Processing:", tf, cancer, "\n")
    mat <- fread(mat_file)
    probes <- mat$probe
    mat[, probe := NULL]

    all_cols     <- colnames(mat)
    tumor_cols   <- grep("_tumor$",   all_cols, value = TRUE)
    healthy_cols <- grep("_healthy$", all_cols, value = TRUE)

    if (length(tumor_cols) == 0 || length(healthy_cols) == 0) {
      cat("  Skipping", cancer, "— missing tumor or healthy samples\n")
      next
    }

    res <- tryCatch(
      run_limma(mat, probes, tumor_cols, healthy_cols),
      error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL }
    )
    if (is.null(res)) next

    res[, cancer    := cancer]
    res[, TF        := tf]
    res[, n_tumor   := length(tumor_cols)]
    res[, n_healthy := length(healthy_cols)]

    tsv_out <- file.path(out_dir, paste0(tf, "_", cancer, "_limma_results.tsv"))
    fwrite(res, tsv_out, sep = "\t")
    cat("  Saved:", tsv_out, "\n")

    all_res[[cancer]] <- res
  }

  if (length(all_res) == 0) next

  # ── per-cancer PDF (one page each) ────────────────────────────────────────
  pdf_out <- file.path(out_dir, paste0(tf, "_limma_volcanos.pdf"))
  pdf(pdf_out, width = 8, height = 6)

  for (cancer in names(all_res)) {
    res <- all_res[[cancer]]
    res[, sig := "NS"]
    res[adj.P.Val < 0.05 & mean_delta_beta >  10, sig := "Hyper"]
    res[adj.P.Val < 0.05 & mean_delta_beta < -10, sig := "Hypo"]

    top_label <- res[sig != "NS"][order(adj.P.Val)][seq_len(min(15, .N))]

    p <- ggplot(res, aes(mean_delta_beta, -log10(adj.P.Val), colour = sig)) +
      geom_point(size = 0.8, alpha = 0.6) +
      scale_colour_manual(values = c(NS = "grey70", Hyper = "#d73027", Hypo = "#4575b4")) +
      geom_text_repel(data = top_label, aes(label = probe), size = 2.5,
                      max.overlaps = 20, show.legend = FALSE) +
      geom_vline(xintercept = c(-10, 10), linetype = "dashed", colour = "grey40") +
      geom_hline(yintercept = -log10(0.05),  linetype = "dashed", colour = "grey40") +
      labs(title = paste0(tf, " — ", cancer,
                          "  (tumor: ", unique(res$n_tumor),
                          ", healthy: ", unique(res$n_healthy), ")"),
           x = "Mean delta-beta (tumor - healthy)",
           y = "-log10(adj. p-value)",
           colour = NULL) +
      theme_bw(base_size = 11)

    print(p)
    cat("  Plotted", cancer, "| Hyper:", sum(res$sig == "Hyper"),
        "| Hypo:", sum(res$sig == "Hypo"), "\n")
  }

  dev.off()
  cat("PDF saved:", pdf_out, "\n")

  # ── combined PDF: all cancers on one page ──────────────────────────────────
  combined <- rbindlist(all_res)
  combined[, sig := "NS"]
  combined[adj.P.Val < 0.05 & mean_delta_beta >  10, sig := "Hyper"]
  combined[adj.P.Val < 0.05 & mean_delta_beta < -10, sig := "Hypo"]

  top_combined <- combined[sig != "NS"][order(adj.P.Val)][
    , .SD[seq_len(min(3, .N))], by = cancer]

  counts <- combined[, .(
    Hyper     = sum(sig == "Hyper"),
    Hypo      = sum(sig == "Hypo"),
    n_tumor   = unique(n_tumor),
    n_healthy = unique(n_healthy)
  ), by = cancer]
  counts[, label := paste0("Hyper: ", Hyper, "\nHypo: ", Hypo,
                            "\nn=", n_tumor, "T / ", n_healthy, "H")]

  x_pos <- combined[, min(mean_delta_beta, na.rm = TRUE)]
  y_pos <- combined[, max(-log10(adj.P.Val), na.rm = TRUE)]
  counts[, x := x_pos]
  counts[, y := y_pos]

  n_cancers <- length(all_res)
  nc <- ceiling(sqrt(n_cancers))
  nr <- ceiling(n_cancers / nc)

  p_all <- ggplot(combined, aes(mean_delta_beta, -log10(adj.P.Val), colour = sig)) +
    geom_point(size = 0.5, alpha = 0.5) +
    scale_colour_manual(values = c(NS = "grey70", Hyper = "#d73027", Hypo = "#4575b4")) +
    geom_text_repel(data = top_combined, aes(label = probe), size = 2,
                    max.overlaps = 10, show.legend = FALSE) +
    geom_text(data = counts, aes(x = x, y = y, label = label),
              colour = "black", size = 2.2, hjust = 0, vjust = 1, inherit.aes = FALSE) +
    geom_vline(xintercept = c(-10, 10), linetype = "dashed", colour = "grey40", linewidth = 0.3) +
    geom_hline(yintercept = -log10(0.05),  linetype = "dashed", colour = "grey40", linewidth = 0.3) +
    facet_wrap(~cancer, ncol = nc) +
    labs(title = paste(tf, "— all cancers"),
         x = "Mean delta-beta (tumor - healthy)",
         y = "-log10(adj. p-value)",
         colour = NULL) +
    theme_bw(base_size = 8) +
    theme(legend.position = "bottom")

  combined_pdf <- file.path(out_dir, paste0(tf, "_limma_all_cancers.pdf"))
  ggsave(combined_pdf, p_all, width = nc * 3.5 + 1, height = nr * 3 + 1, limitsize = FALSE)
  cat("Combined PDF saved:", combined_pdf, "\n")
}

cat("\nDone.\n")
EOF

###########################
# Combine the limma results with the DESeq2 results for each TF and cancer type and plot the pearson correlation between delta-beta and log2FC 
###########################
Rscript - <<'EOF'
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

# ── Thresholds ────────────────────────────────────────────────────────────────
DELTA_BETA_THRESH <- 10
METH_PADJ         <- 0.05
LOG2FC_THRESH     <- 1
EXPR_PADJ         <- 1e-10
PROX_THRESH       <- 2000
MIN_SAMPLES       <- 10

# ── Paths ─────────────────────────────────────────────────────────────────────
limma_dir  <- "./results/methylation/limma"
deseq_dir  <- "./results/expression/deseq2/all_genes"
pg_file    <- "./methylation/probe_gene_pairs_in_motifs_with_tss.tsv.gz"
counts_dir <- "./expression/deseq2/counts_all"
tpm_file   <- "./expression/gene_expression_matrix_protein_coding.tsv"
out_dir    <- "./results/integration"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cat("Loading TPM matrix...\n")
tpm_mat <- fread(tpm_file, header = TRUE, check.names = FALSE)
setnames(tpm_mat, 1L, "tpm_sample_id")
tpm_mat[, cancer_col := sub("^TCGA-([A-Z]+)-.*", "\\1", tpm_sample_id)]
tpm_mat[, st := fifelse(grepl("-01A_1", tpm_sample_id), "Tumor",
                fifelse(grepl("-11A_1", tpm_sample_id), "Healthy", "other"))]
tpm_mat <- tpm_mat[st != "other"]
tpm_mat[, patient    := sub("^TCGA-[A-Z]+-TCGA-([A-Z0-9]+-[A-Z0-9]+)-[0-9]+[A-Z]_.*", "\\1", tpm_sample_id)]
tpm_mat[, sample_key := paste(patient, st, sep = "_")]
tsv_out <- file.path(out_dir, "meth_expr_filtered_pairs.tsv")

if (file.exists(tsv_out)) {
  cat("TSV already exists, loading and going straight to plots ...\n")
  full_dt <- fread(tsv_out)
} else {

  cat("Loading probe-gene pairs ...\n")
  pg <- fread(pg_file)
  setnames(pg, tolower(colnames(pg)))
  pg <- pg[abs(dist_to_tss) <= PROX_THRESH]
  pg <- unique(pg[, .(tf, probe, gene, dist_to_tss)])
  cat(sprintf("  %d unique proximal probe-gene pairs\n", nrow(pg)))

  tfs         <- c("BANP", "NRF1")
  all_results <- list()

  for (tf in tfs) {
    cat("\n========", tf, "========\n")
    TF    <- tf
    tf_pattern <- paste0("^", tf)
    pg_tf <- unique(pg[grepl(tf_pattern, pg$tf), .(probe, gene, dist_to_tss)])
    cat("  pg_tf rows:", nrow(pg_tf), "\n")
    cat("  Example tf values in pg:", head(unique(pg$tf)), "\n")
    cat(sprintf("  %d probe-gene pairs for %s\n", nrow(pg_tf), tf))

    beta_dir <- paste0("./methylation/beta_matrices_", tf)

    limma_files <- list.files(limma_dir,
                               pattern    = paste0("^", tf, "_.*_limma_results\\.tsv$"),
                               full.names = TRUE)
    cancers <- sub(paste0("^", tf, "_"), "",
                   sub("_limma_results\\.tsv$", "", basename(limma_files)))
    cat(sprintf("  %d cancers found\n", length(cancers)))

    for (cancer in cancers) {
      cat(" ", cancer, "...")

      bf <- list.files(beta_dir,
                  pattern    = paste0(cancer, ".*\\.tsv\\.gz$"),
                  full.names = TRUE)
      if (length(bf) == 0) { cat(" no beta matrix\n"); next }

      beta_cols     <- colnames(fread(bf[1], nrows = 0))
      beta_patients <- sub("_(tumor|healthy)$", "", beta_cols[beta_cols != "probe"])
      n_meth_tum    <- sum(grepl("_tumor$",   beta_cols))
      n_meth_hlt    <- sum(grepl("_healthy$", beta_cols))

      cf_all     <- list.files(counts_dir,
                                pattern    = paste0("^RNA_TCGA-", cancer, "-.*_genes_count\\.txt$"),
                                full.names = FALSE)
      pt_expr    <- sub("^RNA_TCGA-([A-Z]+)-TCGA-([A-Z0-9]+-[A-Z0-9]+)-[0-9]+[A-Z]_.*",
                        "\\1-TCGA-\\2", cf_all)
      cf_matched <- cf_all[pt_expr %in% beta_patients]
      n_expr_tum <- sum(grepl("-01[A-Z]_", cf_matched))
      n_expr_hlt <- sum(grepl("-11[A-Z]_", cf_matched))

      if (n_meth_tum < MIN_SAMPLES || n_meth_hlt < MIN_SAMPLES ||
          n_expr_tum < MIN_SAMPLES || n_expr_hlt < MIN_SAMPLES) {
        cat(sprintf(" skipped (meth T=%d H=%d | expr T=%d H=%d)\n",
                    n_meth_tum, n_meth_hlt, n_expr_tum, n_expr_hlt))
        next
      }
      cat(sprintf(" meth T=%d H=%d | expr T=%d H=%d",
                  n_meth_tum, n_meth_hlt, n_expr_tum, n_expr_hlt))
      
      # ── 1. Limma ──────────────────────────────────────────────────────────
      lf        <- file.path(limma_dir, paste0(tf, "_", cancer, "_limma_results.tsv"))
      limma_res <- fread(lf)
      limma_res <- limma_res[adj.P.Val < METH_PADJ & abs(mean_delta_beta) >= DELTA_BETA_THRESH]
      if (nrow(limma_res) == 0) { cat(" no sig methylation\n"); next }

      # ── 2. DESeq2 ─────────────────────────────────────────────────────────
      deseq_files <- list.files(deseq_dir,
                                 pattern    = paste0("^", cancer, "_.*_deseq2_diff_table\\.txt$"),
                                 full.names = TRUE)
      if (length(deseq_files) == 0) { cat(" no DESeq2 file\n"); next }
      deseq_res <- fread(deseq_files[1])
      setnames(deseq_res, tolower(colnames(deseq_res)))
      if (!"gene" %in% colnames(deseq_res)) setnames(deseq_res, 2, "gene")
      deseq_res <- deseq_res[!is.na(padj) & padj < EXPR_PADJ & abs(log2foldchange) >= LOG2FC_THRESH]
      if (nrow(deseq_res) == 0) { cat(" no sig expression\n"); next }

      # ── 3. Merge probe-gene pairs ─────────────────────────────────────────
      merged <- merge(
        limma_res[, .(probe, mean_delta_beta, adj.P.Val)],
        pg_tf, by = "probe"
      )
      merged <- merge(
        merged,
        deseq_res[, .(gene, log2foldchange, padj)],
        by = "gene"
      )
      if (nrow(merged) == 0) { cat(" no overlapping pairs\n"); next }

      # ── 4. Load beta matrix and expression ────────────────────────────────
      beta_mat <- fread(bf[1])
      if (!"probe" %in% colnames(beta_mat)) setnames(beta_mat, 1, "probe")
      beta_mat  <- beta_mat[probe %in% merged$probe]
      beta_long <- melt(beta_mat, id.vars = "probe",
                        variable.name = "col_name", value.name = "beta")
      beta_long[, sample_key := as.character(col_name)]

      cnt_files <- list.files(counts_dir,
                               pattern    = paste0("^RNA_TCGA-", cancer, "-.*_genes_count\\.txt$"),
                               full.names = TRUE)
      bn   <- basename(cnt_files)
      pt <- sub("^RNA_TCGA-([A-Z]+)-TCGA-([A-Z0-9]+-[A-Z0-9]+)-[0-9]+[A-Z]_.*","\\1-TCGA-\\2", bn)
      st   <- ifelse(grepl("-01[A-Z]_", bn), "tumor",
             ifelse(grepl("-11[A-Z]_", bn), "healthy", "other"))
      fkey <- paste(pt, st, sep = "_")

      target_genes <- unique(merged$gene)
      expr_list <- lapply(seq_along(cnt_files), function(i) {
        dt <- tryCatch(
          fread(cnt_files[i], header = FALSE,
                select = c(2L, 3L), col.names = c("gene", "count"),
                showProgress = FALSE),
          error = function(e) NULL
        )
        if (is.null(dt)) return(NULL)
        dt <- dt[gene %in% target_genes]
        if (nrow(dt) == 0) return(NULL)
        dt[, sample_key := fkey[i]]
        dt
      })
      expr_dt <- rbindlist(Filter(Negate(is.null), expr_list))
      expr_dt[, log2_tpm := log2(count + 1)]
      
      # ── 5. Check matched patient overlap ──────────────────────────────────
      test_merge <- merge(
        unique(beta_long[, .(sample_key)]),
        unique(expr_dt[,   .(sample_key)]),
        by = "sample_key"
      )
      n_matched_tum <- sum(grepl("_tumor$",   test_merge$sample_key))
      n_matched_hlt <- sum(grepl("_healthy$", test_merge$sample_key))
      cat(sprintf(" | matched T=%d H=%d", n_matched_tum, n_matched_hlt))
      if (n_matched_tum < MIN_SAMPLES || n_matched_hlt < MIN_SAMPLES) {
        cat(sprintf(" skipped: too few matched patients (T=%d H=%d)\n",
                    n_matched_tum, n_matched_hlt))
        next
      }

      # ── 6+7. Pearson r + means (vectorized) ───────────────────────────────
      beta_long[, group := ifelse(grepl("_tumor$",   sample_key), "tumor",
                          ifelse(grepl("_healthy$", sample_key), "healthy", "other"))]
      expr_dt[,  group := ifelse(grepl("_tumor$",   sample_key), "tumor",
                          ifelse(grepl("_healthy$", sample_key), "healthy", "other"))]

      # means
      beta_means <- beta_long[group != "other",
        .(mean_beta_tumor     = mean(beta[group == "tumor"],   na.rm = TRUE),
          mean_beta_healthy   = mean(beta[group == "healthy"], na.rm = TRUE),
          median_beta_tumor   = median(beta[group == 'tumor'],   na.rm = TRUE),
          median_beta_healthy = median(beta[group == 'healthy'], na.rm = TRUE)),
        by = probe]

      expr_means <- expr_dt[group != "other",
        .(mean_expr_tumor     = mean(log2_tpm[group == "tumor"],   na.rm = TRUE),
          mean_expr_healthy   = mean(log2_tpm[group == "healthy"], na.rm = TRUE),
          median_expr_tumor   = median(log2_tpm[group == 'tumor'],   na.rm = TRUE),
          median_expr_healthy = median(log2_tpm[group == 'healthy'], na.rm = TRUE)),
        by = gene]

      merged <- merge(merged, beta_means, by = "probe",  all.x = TRUE)
      merged <- merge(merged, expr_means, by = "gene",   all.x = TRUE)

      # TPM means — matched patients only
      tpm_genes <- intersect(merged$gene, names(tpm_mat))
      if (length(tpm_genes) > 0) {
        # extract patient IDs from beta_long sample keys e.g. "BLCA-TCGA-2F-A9KO_tumor" -> "2F-A9KO"
        beta_patients_matched <- unique(sub(paste0("^", cancer, "-TCGA-([A-Z0-9]+-[A-Z0-9]+)_(tumor|healthy)$"),
                                            "\\1", unique(beta_long$sample_key)))
        expr_patients_matched <- unique(sub(paste0("^", cancer, "-TCGA-([A-Z0-9]+-[A-Z0-9]+)_(tumor|healthy)$"),
                                    "\\1", unique(expr_dt$sample_key)))
        matched_patients      <- intersect(beta_patients_matched, expr_patients_matched)

        tpm_matched <- tpm_mat[patient %in% matched_patients & cancer_col == cancer]
        tpm_long <- melt(tpm_matched[, c("patient", "st", tpm_genes), with = FALSE],
                        id.vars = c("patient", "st"),
                        variable.name = "gene", value.name = "tpm")
        tpm_means <- tpm_long[st != "other",
          .(mean_tpm_tumor   = mean(tpm[st == "Tumor"],   na.rm = TRUE),
            mean_tpm_healthy = mean(tpm[st == "Healthy"], na.rm = TRUE)),
          by = gene]
        merged <- merge(merged, tpm_means, by = "gene", all.x = TRUE)
      } else {
        merged[, mean_tpm_tumor := NA_real_]
        merged[, mean_tpm_healthy := NA_real_]
      }

      merged <- merged[mean_tpm_tumor > 1 | mean_tpm_healthy > 1]

      # Pearson r (vectorized via merge)
      both <- merge(
        beta_long[, .(probe, sample_key, beta)],
        expr_dt[,   .(gene,  sample_key, log2_tpm)],
        by = "sample_key",
        allow.cartesian = TRUE
      )
      both <- merge(both, merged[, .(probe, gene)], by = c("probe", "gene"))

      cor_res <- both[, {
        ok <- is.finite(beta) & is.finite(log2_tpm)
        x <- beta[ok]; y <- log2_tpm[ok]
        if (length(x) >= 5 && length(unique(x)) > 1 && length(unique(y)) > 1) {
          ct <- tryCatch(cor.test(x, y, method = "pearson"), error = function(e) NULL)
          if (!is.null(ct)) .(pearson_r = unname(ct$estimate), pearson_p = ct$p.value)
          else .(pearson_r = NA_real_, pearson_p = NA_real_)
        } else .(pearson_r = NA_real_, pearson_p = NA_real_)
      }, by = .(probe, gene)]

      merged[, c("pearson_r", "pearson_p") := NULL]
      merged <- merge(merged, cor_res, by = c("probe", "gene"), all.x = TRUE)
      # ── Wilcoxon test per probe-gene pair ─────────────────────────────────────
      wilcox_res <- both[, {
        grp      <- ifelse(grepl("_tumor$",   sample_key), "tumor", "healthy")
        expr_ok  <- is.finite(log2_tpm)
        meth_ok  <- is.finite(beta)

        expr_p <- if (sum(grp == "tumor" & expr_ok) >= 2 &&
                      sum(grp == "healthy" & expr_ok) >= 2) {
          wt <- tryCatch(suppressWarnings(
            wilcox.test(log2_tpm[grp == "tumor"   & expr_ok],
                        log2_tpm[grp == "healthy" & expr_ok], exact = FALSE)),
            error = function(e) NULL)
          if (!is.null(wt)) wt$p.value else NA_real_
        } else NA_real_

        meth_p <- if (sum(grp == "tumor" & meth_ok) >= 2 &&
                      sum(grp == "healthy" & meth_ok) >= 2) {
          wt <- tryCatch(suppressWarnings(
            wilcox.test(beta[grp == "tumor"   & meth_ok],
                        beta[grp == "healthy" & meth_ok], exact = FALSE)),
            error = function(e) NULL)
          if (!is.null(wt)) wt$p.value else NA_real_
        } else NA_real_

        .(expr_wilcox_p = expr_p, meth_wilcox_p = meth_p)
      }, by = .(probe, gene)]

      merged <- merge(merged, wilcox_res, by = c("probe", "gene"), all.x = TRUE)
      merged <- merged[!is.na(expr_wilcox_p) & !is.na(meth_wilcox_p) &
                        expr_wilcox_p < 0.05 & meth_wilcox_p < 0.05]

      merged[, cancer := cancer]
      merged[, TF     := tf]
      n <- nrow(merged)
      cat(sprintf(" | %d pairs | r computed for %d\n", n, sum(!is.na(merged$pearson_r))))
      all_results[[paste(tf, cancer, sep = "_")]] <- merged
    }
  }

  if (length(all_results) == 0) stop("No results found.")

  full_dt <- rbindlist(all_results, fill = TRUE)
  fwrite(full_dt, tsv_out, sep = "\t")
  cat("\nTable saved:", tsv_out, "\n")
}

# ── Plot ──────────────────────────────────────────────────────────────────────
r_scale <- scale_colour_gradient2(
  low = "#2166ac", mid = "white", high = "#d6604d",
  midpoint = 0, limits = c(-1, 1), na.value = "grey70", name = "Pearson r"
)

tfs <- unique(full_dt$TF)

for (tf in tfs) {
  dt <- full_dt[TF == tf]
  if (nrow(dt) == 0) next

  top_labels <- dt[!is.na(pearson_r),
    ][, score := abs(log2foldchange) * abs(pearson_r)
    ][order(-score), head(.SD, 20), by = cancer]
  top_labels[, label := paste0(gene, " | ", cancer)]

  n_cancers <- length(unique(dt$cancer))
  nc <- ceiling(sqrt(n_cancers))
  nr <- ceiling(n_cancers / nc)

  p <- ggplot(dt, aes(mean_delta_beta, log2foldchange, colour = pearson_r)) +
    geom_point(alpha = 0.8, size = 1.8) +
    geom_text_repel(
      data          = top_labels,
      aes(label     = label),
      colour        = "black",
      size          = 2.2,
      fontface      = "italic",
      max.overlaps  = 30,
      segment.size  = 0.3,
      segment.color = "grey50",
      box.padding   = 0.3,
      show.legend   = FALSE
    ) +
    r_scale +
    facet_wrap(~cancer, scales = "free") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
    labs(
      title    = paste0(tf, " methylation delta vs expression log2FC"),
      subtitle = paste0("Proximal probes =< ", PROX_THRESH, " bp | |deltabeta| >= ", DELTA_BETA_THRESH,
                        "%, padj < ", METH_PADJ, " | |log2FC| >= ", LOG2FC_THRESH,
                        ", padj < ", EXPR_PADJ,
                        " | >=", MIN_SAMPLES, " matched tumor+healthy patients",
                        " | top 20 labelled by |log2FC|x|r|"),
      x = "Methylation delta (tumor - healthy, %)",
      y = "Expression log2FC (tumor vs healthy)"
    ) +
    theme_bw(base_size = 9) +
    theme(strip.text      = element_text(face = "bold"),
          legend.position = "right",
          plot.subtitle   = element_text(size = 7, colour = "grey40"))

  out_pdf <- file.path(out_dir, paste0(tf, "_meth_expr_scatter.pdf"))
  ggsave(out_pdf, p,
         width = nc * 4 + 1.5, height = nr * 3.5 + 1, limitsize = FALSE)
  cat("Saved:", out_pdf, "\n")
}

cat("All done.\n")
EOF


#####################################################################################################################################
# Script to generate scatter plots for each gene in tsv file it takes the list of genes that are in the first column of the tsv file: 
#####################################################################################################################################
Rscript - << 'EOF'
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
  library(gtable)
  library(parallel)
})

# ── Config ────────────────────────────────────────────────────────────────────
N_CORES    <- max(1L, detectCores() - 2L)   # or set e.g. 6L
tf_name    <- "NRF1"
pairs_tsv  <- "./results/integration/meth_expr_filtered_pairs.tsv"
beta_dir   <- paste0("./methylation/beta_matrices_", tf_name)
tpm_file <- "./expression/gene_expression_matrix_protein_coding.tsv"
out_dir    <- file.path("./results/integration/gene_cancer_plots", tf_name)
stats_out  <- paste0("./results/integration/meth_expr_wilcoxon_", tf_name, "_stats.tsv")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

min_n_cor        <- 3L
min_n_wilcox     <- 2L
expression_threshold <- 1
side_box_size    <- unit(6, "cm")
top_margin_size  <- unit(1.5, "cm")

cat(sprintf("Running with %d cores\n", N_CORES))

# ── Helpers ───────────────────────────────────────────────────────────────────
fmt_p <- function(x) {
  if (is.na(x)) return("NA")
  if (x < 2.2e-16) return("<2.2e-16")
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

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
  gsub("(^_|_$)", "", gsub("_+", "_", x))
}

# ── Load pairs, pick best probe per TF-gene-cancer ───────────────────────────
if (!file.exists(pairs_tsv)) stop("pairs TSV not found: ", pairs_tsv)
pairs <- fread(pairs_tsv)
pairs <- pairs[TF == tf_name]

todo <- pairs[!is.na(pearson_r)][order(-abs(pearson_r))][, .SD[1], by = .(TF, gene, cancer)]
cat(sprintf("Processing %d gene-cancer pairs\n", nrow(todo)))

# ── Load TPM matrix once (shared across all workers via fork) ─────────────────
cat("Loading TPM matrix...\n")
tpm_mat <- fread(tpm_file, header = TRUE, check.names = FALSE)
setnames(tpm_mat, 1L, "tpm_sample_id")
tpm_mat[, patient    := sub("^TCGA-[A-Z]+-TCGA-([A-Z0-9]+-[A-Z0-9]+)-[0-9]+[A-Z]_.*", "\\1", tpm_sample_id)]
tpm_mat[, cancer_col := sub("^TCGA-([A-Z]+)-.*", "\\1", tpm_sample_id)]
tpm_mat[, st := fifelse(grepl("-01A_1", tpm_sample_id), "Tumor",
                fifelse(grepl("-11A_1", tpm_sample_id), "Healthy", "other"))]
tpm_mat[, sample_key := paste(patient, st, sep = "_")]
tpm_mat <- tpm_mat[st != "other"]
cat(sprintf("TPM matrix loaded: %d samples, %d genes\n",
            nrow(tpm_mat), ncol(tpm_mat) - 5L))

# ── Worker: one pair → stats data.table (or NULL) ────────────────────────────
process_pair <- function(i) {
  row           <- todo[i]
  target_gene   <- row$gene
  target_cancer <- row$cancer
  target_probe  <- row$probe
  tag <- sprintf("[%d/%d] %s | %s", i, nrow(todo), target_cancer, target_gene)

  # Load beta (one probe row)
  bf <- list.files(beta_dir,
    pattern = paste0(target_cancer, ".*\\.tsv\\.gz$"), full.names = TRUE)
  if (length(bf) == 0) { cat(tag, "– no beta matrix\n"); return(NULL) }

  beta_mat <- tryCatch(fread(bf[1]), error = function(e) NULL)
  if (is.null(beta_mat)) { cat(tag, "– fread failed\n"); return(NULL) }
  if (!"probe" %in% colnames(beta_mat)) setnames(beta_mat, 1, "probe")
  .probe   <- target_probe
  beta_mat <- beta_mat[probe == .probe]
  if (nrow(beta_mat) == 0) { cat(tag, "– probe missing\n"); return(NULL) }

  beta_long <- melt(beta_mat, id.vars = "probe",
                    variable.name = "col_name", value.name = "meth_beta")
  beta_long[, patient     := sub("^[A-Z]+-TCGA-([A-Z0-9]+-[A-Z0-9]+)_(tumor|healthy)$",
                                "\\1", as.character(col_name))]
  beta_long[, sample_type := sub("^.*_(tumor|healthy)$", "\\1", as.character(col_name))]
  beta_long[, sample_type := fifelse(sample_type == "tumor", "Tumor", "Healthy")]
  beta_long[, sample_key  := paste(patient, sample_type, sep = "_")]
  if (max(beta_long$meth_beta, na.rm = TRUE) <= 1.5)
    beta_long[, meth_percent := meth_beta * 100] else
    beta_long[, meth_percent := meth_beta]

  # Load expression from TPM matrix
  .gene <- target_gene
  if (!(.gene %in% names(tpm_mat))) { cat(tag, "– gene not in TPM matrix\n"); return(NULL) }

  expr_dt <- tpm_mat[cancer_col == target_cancer, .(sample_key, count = get(.gene))]
  if (nrow(expr_dt) == 0) { cat(tag, "– no TPM samples for this cancer\n"); return(NULL) }

  # Merge
  d <- merge(
    beta_long[, .(sample_key, patient, sample_type, meth_percent)],
    expr_dt[,   .(sample_key, count)],
    by = "sample_key"
  )
  d <- d[sample_type %in% c("Tumor", "Healthy") & is.finite(meth_percent) & is.finite(count)]
  d[, log2_expr   := log2(count + 1)]
  d[, sample_type := factor(sample_type, levels = c("Healthy", "Tumor"))]
  if (nrow(d) == 0) { cat(tag, "– no matched samples\n"); return(NULL) }
  med_tumor_expr   <- median(d[sample_type == "Tumor",   log2_expr],    na.rm = TRUE)
  med_healthy_expr <- median(d[sample_type == "Healthy", log2_expr],    na.rm = TRUE)
  med_tumor_meth   <- median(d[sample_type == "Tumor",   meth_percent], na.rm = TRUE)
  med_healthy_meth <- median(d[sample_type == "Healthy", meth_percent], na.rm = TRUE)

  # Stats
  cor_res     <- safe_pearson(d$meth_percent, d$log2_expr, min_n_cor)
  expr_wilcox <- safe_wilcox(d$log2_expr,    as.character(d$sample_type), min_n_wilcox)
  meth_wilcox <- safe_wilcox(d$meth_percent, as.character(d$sample_type), min_n_wilcox)
  n_healthy   <- sum(d$sample_type == "Healthy")
  n_tumor     <- sum(d$sample_type == "Tumor")

  stats_row <- data.table(
    TF = tf_name, gene = target_gene, cancer = target_cancer, probe = target_probe,
    n_total = nrow(d), n_healthy, n_tumor,
    pearson_n = cor_res$n, pearson_r = cor_res$r, pearson_p = cor_res$p,
    expr_wilcox_p = expr_wilcox$p, meth_wilcox_p = meth_wilcox$p
  )

  # Wilcoxon filter
  if (is.na(expr_wilcox$p) || is.na(meth_wilcox$p) ||
      expr_wilcox$p >= 0.05 || meth_wilcox$p >= 0.05) {
    cat(tag, "– skipped (Wilcoxon ns)\n")
    return(stats_row)
  }

  # Build plots
  x_lim    <- c(0, 100)
  x_breaks <- c(0, 25, 50, 75, 100)
  expr_min <- min(d$log2_expr, na.rm = TRUE)
  expr_max <- max(d$log2_expr, na.rm = TRUE)
  expr_rng <- max(expr_max - expr_min, 1)
  expr_ylim   <- c(max(0, expr_min - 0.08 * expr_rng), expr_max + 0.18 * expr_rng)
  expr_breaks <- pretty(expr_ylim, n = 5)

  stat_label    <- paste0("n = ", cor_res$n,
                          "\nr = ", ifelse(is.na(cor_res$r), "NA", sprintf("%.2f", cor_res$r)),
                          "\np = ", fmt_p(cor_res$p))
  main_title    <- paste0("Scatter plot - ", target_cancer, " | ", target_gene, " | ", target_probe)
  main_subtitle <- paste0("Tumor n = ", n_tumor, " ; Healthy n = ", n_healthy)

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
             color = "red", hjust = 1, vjust = -0.5, size = 4) +
    annotate("text", x = Inf, y = Inf, label = stat_label,
             hjust = 1.05, vjust = 1.1, size = 4) +
    scale_x_continuous(limits = x_lim, breaks = x_breaks, expand = c(0.01, 0.01)) +
    scale_y_continuous(limits = expr_ylim, breaks = expr_breaks, expand = c(0.01, 0.01)) +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 15) +
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
             label = p_label(expr_wilcox$p), size = 4.3, lineheight = 0.9) +
        annotate("text", x = 1, y = expr_ylim[1] + 0.04 * diff(expr_ylim),
             label = sprintf("med=%.2f", med_healthy_expr),
             size = 4.5, colour = "black", hjust = 0.5) +
    annotate("text", x = 2, y = expr_ylim[1] + 0.04 * diff(expr_ylim),
             label = sprintf("med=%.2f", med_tumor_expr),
             size = 4.5, colour = "black", hjust = 0.5) +
    scale_y_continuous(limits = expr_ylim, breaks = expr_breaks, expand = c(0.01, 0.01)) +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 12) +
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
    annotate("segment", x = 88, xend = 98, y = 1.5, yend = 1.5, linewidth = 0.6) +
    annotate("text", x = 93, y = 1.73,
             label = p_label(meth_wilcox$p), size = 4.3, lineheight = 0.9) +
    annotate("text", x = med_healthy_meth, y = 0.6,
             label = sprintf("med=%.1f%%", med_healthy_meth),
             size = 4.5, colour = "black", hjust = 0.5) +
    annotate("text", x = med_tumor_meth, y = 1.6,
             label = sprintf("med=%.1f%%", med_tumor_meth),
             size =4.5, colour = "black", hjust = 0.5) +
    scale_x_continuous(limits = x_lim, breaks = x_breaks, expand = c(0.01, 0.01)) +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none", plot.margin = margin(5, 10, 0, 0))

  draw_combined_plot <- function() {
    grid.newpage()
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
      heights = unit.c(top_margin_size, unit(0.65, "cm"), unit(0.45, "cm"),
                       unit(1, "null"), side_box_size, unit(0.65, "cm"))
    )))
    grid.text(main_title, x = 0.5, y = 0.5, just = "center",
              gp = gpar(fontsize = 20, fontface = "bold"),
              vp = viewport(layout.pos.row = 2, layout.pos.col = 1:3))
    grid.text(main_subtitle, x = 0.5, y = 0.5, just = "center",
              gp = gpar(fontsize = 14),
              vp = viewport(layout.pos.row = 3, layout.pos.col = 1:3))
    grid.text("Gene expression (log2(TPM + 1))", x = 0.5, y = 0.5, rot = 90, just = "center",
              gp = gpar(fontsize = 15),
              vp = viewport(layout.pos.row = 4, layout.pos.col = 1))
    grid.draw(editGrob(expr_grob,
              vp = viewport(layout.pos.row = 4, layout.pos.col = 2)))
    grid.draw(editGrob(scatter_grob,
              vp = viewport(layout.pos.row = 4, layout.pos.col = 3)))
    grid.draw(editGrob(meth_grob,
              vp = viewport(layout.pos.row = 5, layout.pos.col = 3)))
    grid.text("Probe methylation (%)", x = 0.5, y = 0.5, just = "center",
              gp = gpar(fontsize = 15),
              vp = viewport(layout.pos.row = 6, layout.pos.col = 3))
  }

  base_name <- paste0(tf_name, "__", safe_name(target_gene), "__",
                      safe_name(target_cancer), "__", safe_name(target_probe))
  out_pdf   <- file.path(out_dir, paste0(base_name, "_scatter_axis_boxplots.pdf"))

  tryCatch({
    pdf(out_pdf, width = 20, height = 16)
    draw_combined_plot()
    dev.off()
    cat(tag, "– PDF saved\n")
  }, error = function(e) {
    try(dev.off(), silent = TRUE)
    cat(tag, "– PDF ERROR:", conditionMessage(e), "\n")
  })

  stats_row
}

# ── Run in parallel ───────────────────────────────────────────────────────────
t0 <- proc.time()

results <- mclapply(
  seq_len(nrow(todo)),
  process_pair,
  mc.cores       = N_CORES,
  mc.preschedule = FALSE
)

elapsed <- proc.time() - t0
cat(sprintf("\nDone in %.1f s (%.1f min)\n", elapsed["elapsed"], elapsed["elapsed"] / 60))

# ── Write stats TSV once ──────────────────────────────────────────────────────
ok        <- sapply(results, is.data.table)
all_stats <- rbindlist(results[ok])
cat(sprintf("Pairs processed: %d / %d | stats rows: %d\n",
            sum(ok), nrow(todo), nrow(all_stats)))

if (nrow(all_stats) > 0) {
  fwrite(all_stats, stats_out, sep = "\t")
  cat("Stats written to:", stats_out, "\n")
}

errs <- results[!ok]
if (length(errs) > 0) {
  cat(sprintf("\n%d worker errors:\n", length(errs)))
  for (e in errs) cat(" ", conditionMessage(e), "\n")
}
EOF

###############################
# hyper + down & hypo + up
###############################
Rscript - <<'EOF'
suppressPackageStartupMessages(library(data.table))
dt <- fread("./results/integration/meth_expr_filtered_pairs.tsv")
dt[, direction := fifelse(mean_delta_beta > 0 & log2foldchange < 0, "Hyper_Down",
                 fifelse(mean_delta_beta < 0 & log2foldchange > 0, "Hypo_Up", "Other"))]
consistent <- dt[direction != "Other"]
cat("=== Per direction ===\n")
print(consistent[, .(n_genes = uniqueN(gene)), by = .(TF, direction)][order(TF, direction)])
cat("\n=== Unique genes total ===\n")
print(consistent[, .(n_genes = uniqueN(gene)), by = TF])
EOF

# === Per direction ===
#        TF  direction n_genes
#    <char>     <char>   <int>
# 1:   BANP Hyper_Down      40
# 2:   BANP    Hypo_Up      11
# 3:   NRF1 Hyper_Down     132
# 4:   NRF1    Hypo_Up      25
# unique genes is lower of just the sum because some genes are in both hyper_down and hypo_up. but in different cancer types so the final sum is lower 
# === Unique genes total ===
#        TF n_genes
#    <char>   <int>
# 1:   BANP      49
# 2:   NRF1     153

######################################################################################################
# Create consistent pairs TSV (Hyper+Down and Hypo+Up only) + filter on median 
######################################################################################################
Rscript - <<'EOF'
suppressPackageStartupMessages(library(data.table))

full_dt <- fread("./results/integration/meth_expr_filtered_pairs.tsv")

full_dt[, direction := fifelse(mean_delta_beta > 0 & log2foldchange < 0, "Hyper_Down",
                      fifelse(mean_delta_beta < 0 & log2foldchange > 0, "Hypo_Up",
                              "Other"))]

filtered_dt <- full_dt[direction %in% c("Hyper_Down", "Hypo_Up")]

cat("=== Before median delta beta filter ===\n")
print(filtered_dt[, .(n_genes = uniqueN(gene)), by = .(TF, direction)][order(TF, direction)])
cat("Unique genes total:\n")
print(filtered_dt[, .(n_genes = uniqueN(gene)), by = TF])

# Filter on median delta beta >= 10%
filtered_dt[, median_delta_beta := median_beta_tumor - median_beta_healthy]
filtered_dt <- filtered_dt[abs(median_delta_beta) >= 10]

cat("\n=== After median delta beta filter (>= 10%) ===\n")
print(filtered_dt[, .(n_genes = uniqueN(gene)), by = .(TF, direction)][order(TF, direction)])
cat("Unique genes total:\n")
print(filtered_dt[, .(n_genes = uniqueN(gene)), by = TF])

# Filter on |delta median log2(TPM+1)| >= 1
cat("[EXPR] Loading TPM matrix...\n")
tpm_mat <- fread("./expression/gene_expression_matrix_protein_coding.tsv",
                 header = TRUE, check.names = FALSE)
setnames(tpm_mat, 1L, "tpm_sample_id")
tpm_mat[, cancer_col := sub("^TCGA-([A-Z]+)-.*", "\\1", tpm_sample_id)]
tpm_mat[, st         := fifelse(grepl("-01A_1", tpm_sample_id), "Tumor",
                        fifelse(grepl("-11A_1", tpm_sample_id), "Healthy", "other"))]
tpm_mat <- tpm_mat[st != "other"]

genes_needed <- intersect(unique(filtered_dt$gene), names(tpm_mat))
cat("[EXPR] Genes in pairs:", uniqueN(filtered_dt$gene),
    "| found in TPM matrix:", length(genes_needed), "\n")

tpm_long <- melt(
  tpm_mat[, c("cancer_col", "st", genes_needed), with = FALSE],
  id.vars       = c("cancer_col", "st"),
  variable.name = "gene",
  value.name    = "tpm"
)[, gene := as.character(gene)]

med_expr <- tpm_long[, .(
  median_log2tpm_healthy = median(log2(tpm[st == "Healthy"] + 1), na.rm = TRUE),
  median_log2tpm_tumor   = median(log2(tpm[st == "Tumor"]   + 1), na.rm = TRUE)
), by = .(gene, cancer = cancer_col)]
med_expr[, delta_median_log2tpm := median_log2tpm_tumor - median_log2tpm_healthy]

filtered_dt <- merge(filtered_dt, med_expr, by = c("gene", "cancer"), all.x = TRUE)
filtered_dt <- filtered_dt[abs(delta_median_log2tpm) >= 1]

cat("\n=== After |delta median log2(TPM+1)| >= 1 ===\n")
print(filtered_dt[, .(n_genes = uniqueN(gene)), by = .(TF, direction)][order(TF, direction)])
cat("Unique genes total:\n")
print(filtered_dt[, .(n_genes = uniqueN(gene)), by = TF])

# Add motif_instance_id from probe-gene pairs file
pg <- fread("./methylation/probe_gene_pairs_in_motifs_with_tss.tsv.gz")
setnames(pg, tolower(colnames(pg)))
pg <- unique(pg[, .(probe, motif_instance_id)])

filtered_dt <- merge(filtered_dt, pg, by = "probe", all.x = TRUE)

# Reorder so motif_instance_id is after probe
cols       <- colnames(filtered_dt)
other_cols <- setdiff(cols, "motif_instance_id")
probe_pos  <- which(other_cols == "probe")
new_order  <- c(other_cols[1:probe_pos], "motif_instance_id", other_cols[(probe_pos + 1):length(other_cols)])
setcolorder(filtered_dt, new_order)

fwrite(filtered_dt, "./results/integration/meth_expr_consistent_pairs.tsv", sep = "\t")
cat("\nDone. Kept", nrow(filtered_dt), "pairs (", uniqueN(filtered_dt$gene), "unique genes)\n")
EOF

# === Before median delta beta filter ===
#        TF  direction n_genes
#    <char>     <char>   <int>
# 1:   BANP Hyper_Down      40
# 2:   BANP    Hypo_Up      11
# 3:   NRF1 Hyper_Down     132
# 4:   NRF1    Hypo_Up      25
# Unique genes total:
#        TF n_genes
#    <char>   <int>
# 1:   BANP      49
# 2:   NRF1     153

# === After median delta beta filter (>= 10%) ===
#        TF  direction n_genes
#    <char>     <char>   <int>
# 1:   BANP Hyper_Down      30
# 2:   BANP    Hypo_Up      10
# 3:   NRF1 Hyper_Down     106
# 4:   NRF1    Hypo_Up      24
# Unique genes total:
#        TF n_genes
#    <char>   <int>
# 1:   BANP      38
# 2:   NRF1     126
# [EXPR] Loading TPM matrix...
# [EXPR] Genes in pairs: 161 | found in TPM matrix: 161

# === After |delta median log2(TPM+1)| >= 1 ===
#        TF  direction n_genes
#    <char>     <char>   <int>
# 1:   BANP Hyper_Down      23
# 2:   BANP    Hypo_Up       3
# 3:   NRF1 Hyper_Down      81
# 4:   NRF1    Hypo_Up      17
# Unique genes total:
#        TF n_genes
#    <char>   <int>
# 1:   NRF1      95
# 2:   BANP      25


##################################
# Copy consistent gene-cancer plots to a new directory 
##################################
Rscript - <<'EOF'
suppressPackageStartupMessages(library(data.table))

consistent <- fread("./results/integration/meth_expr_consistent_pairs.tsv")

tf_name  <- "NRF1"
in_dir   <- file.path("./results/integration/gene_cancer_plots", tf_name)
out_dir  <- file.path("./results/integration/gene_cancer_plots_consistent", tf_name)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pairs <- unique(consistent[TF == tf_name, .(gene, cancer)])

copied <- 0
missing <- 0
for (i in seq_len(nrow(pairs))) {
  gene   <- pairs$gene[i]
  cancer <- pairs$cancer[i]
  pattern <- paste0(tf_name, "__", gene, "__", cancer, "__.*\\.pdf$")
  files <- list.files(in_dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) { missing <- missing + 1; next }
  file.copy(files, out_dir)
  copied <- copied + length(files)
}
if (missing > 0) {
  cat("Missing pairs:\n")
  for (i in seq_len(nrow(pairs))) {
    gene   <- pairs$gene[i]
    cancer <- pairs$cancer[i]
    pattern <- paste0(tf_name, "__", gene, "__", cancer, "__.*\\.pdf$")
    files <- list.files(in_dir, pattern = pattern, full.names = TRUE)
    if (length(files) == 0) cat(sprintf("  %s  %s\n", gene, cancer))
  }
}
cat(sprintf("Copied %d PDFs | %d gene-cancer pairs had no matching file\n", copied, missing))
EOF

##############################################################################
# CpG Probe Distance to TSS Distribution 
# Description : Plot the density of CpG probe distances to TSS for
#               methylation-expression consistent pairs (Hyper+Down / Hypo+Up),
#               stratified by transcription factor.
# Input       : results/integration/meth_expr_consistent_pairs.tsv
# Output      : results/integration/dist_to_tss_distribution.pdf
##############################################################################
Rscript - <<'EOF'
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

dt <- fread("./results/integration/meth_expr_consistent_pairs.tsv")

dt_uniq <- unique(dt[, .(TF, direction, dist_to_tss, gene)])[
  , .SD[1], by = .(TF, direction, gene)]
dt_uniq[, TF_label := paste0(TF, " (n=", uniqueN(gene), " genes)"), by = TF]

p <- ggplot(dt_uniq, aes(x = dist_to_tss, fill = direction, colour = direction)) +
  geom_density(alpha = 0.35, linewidth = 0.7) +
  geom_rug(alpha = 0.8, linewidth = 0.8, length = unit(0.06, "npc"), sides = "b") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey30", linewidth = 0.5) +
  scale_fill_manual(values = c(Hyper_Down = "#d73027", Hypo_Up = "#4575b4"),
                    labels = c(Hyper_Down = "Hyper + Down", Hypo_Up = "Hypo + Up")) +
  scale_colour_manual(values = c(Hyper_Down = "#d73027", Hypo_Up = "#4575b4"),
                      labels = c(Hyper_Down = "Hyper + Down", Hypo_Up = "Hypo + Up")) +
  scale_x_continuous(limits = c(-200, 1000),
                     breaks = seq(-200, 1000, 200)) +
  facet_wrap(~TF_label, ncol = 1) +
  labs(
    x      = "Distance to TSS (bp)",
    y      = "Density",
    fill   = "Direction",
    colour = "Direction",
    title  = "CpG probe distance to TSS — consistent pairs"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position  = "top",
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 13)
  )

ggsave("./results/integration/dist_to_tss_distribution.pdf",
       p, width = 7, height = 5)
cat("Plot saved.\n")
EOF

##################################
# CancerMine Database analysis
##################################
#####################################################
# Look into the conistent genes to see if they are labeled as tumor suppressors or oncogenes in the CancerMine database
#####################################################
wget -O ./results/integration/cancermine.tsv https://zenodo.org/records/23597453/files/cancermine_collated.tsv?download=1

# label the consistent genes as tumor suppressors or oncogenes based on the CancerMine database
# seeing if its in the same cancer type as the ones i identified
# Redo this 
awk 'BEGIN {FS=OFS="\t"}
NR==FNR {
    if (FNR>1) genes[$3]=1
    next
}
FNR>1 && ($7 in genes) &&
($2=="Oncogene" || $2=="Tumor_Suppressor" || $2=="Driver") {
    print tolower($4)
}' \
/data/najemd/TF_binding_cancer/results/integration/meth_expr_consistent_pairs.tsv \
/data/najemd/TF_binding_cancer/results/integration/cancermine.tsv \
| sort -u \
> /data/najemd/TF_binding_cancer/results/integration/cancermine_cancer_names_to_review.txt

# added the TCGA cancer codes to the CancerMine cancer names to review file manually

awk 'BEGIN{FS=OFS="\t"
  # load the cancer name → TCGA code mapping before processing any file
  map_file = "/data/najemd/TF_binding_cancer/results/integration/cancermine_cancer_names_to_review.txt"
  getline < map_file   # skip header
  while ((getline line < map_file) > 0) {
    split(line, f, "\t")
    cm_codes[tolower(f[1])] = f[2]   # e.g. "kidney cancer" → "KIRC/KIRP"
  }
  close(map_file)
}
NR==FNR{
  if(NR>1 && ($2=="Oncogene"||$2=="Tumor_Suppressor"||$2=="Driver")){
    g=$7; r=$2; c=tolower($4); n=$8
    if(!seen_role[g,r]++)  roles[g]   = roles[g]   ? roles[g]"|"r              : r
    if(!seen_det[g,r,c]++) details[g] = details[g] ? details[g]";"r":"$4"("n")" : r":"$4"("n")"
    if(c in cm_codes){
      ncodes = split(cm_codes[c], code_arr, "/")
      for(k=1; k<=ncodes; k++){
        code = code_arr[k]
        if(!seen_match[g,code,r]++)
          match_role[g,code] = match_role[g,code] ? match_role[g,code]"|"r : r
      }
    }
  }
  next
}
FNR==1{ print $0,"cancermine_role","cancermine_details","cancermine_role_same_cancer"; next }
{
  g=$3; cancer=$4
  print $0,
        (g in roles              ? roles[g]               : "NA"),
        (g in details            ? details[g]             : "NA"),
        ((g SUBSEP cancer) in match_role ? match_role[g,cancer] : "NA")
}
' \
/data/najemd/TF_binding_cancer/results/integration/cancermine.tsv \
/data/najemd/TF_binding_cancer/results/integration/meth_expr_consistent_pairs.tsv \
> /data/najemd/TF_binding_cancer/results/integration/meth_expr_consistent_pairs_cancermine_labeled.tsv

# summarize the number of consistent genes that are labeled as tumor suppressors or oncogenes in the CancerMine database
awk 'BEGIN{FS=OFS="\t"}
FNR==1{
  for(i=1;i<=NF;i++){
    if($i=="TF")                          tc=i
    if($i=="gene")                        gc=i
    if($i=="cancermine_role")             rc=i
    if($i=="cancermine_role_same_cancer") sc=i
  }
  print "TF","n_Oncogene","n_Tumor_Suppressor","n_Driver","n_any_labeled","n_matched_cancer","n_NA"
  next
}
{
  tf=$tc; g=$gc; role=$rc; same=$sc
  tfs[tf]=1

  if(role=="NA" || role==""){
    if(!seen_na[tf,g]++) na[tf]++
    next
  }
  if(!seen_any[tf,g]++) total[tf]++
  n=split(role, roles, "|")
  for(i=1;i<=n;i++){
    r=roles[i]
    if(!seen[tf,r,g]++) cnt[tf,r]++
  }
  if(same!="NA" && same!="")
    if(!seen_same[tf,g]++) matched[tf]++
}
END{
  for(tf in tfs)
    print tf,
          (tf,"Oncogene")         in cnt ? cnt[tf,"Oncogene"]         : 0,
          (tf,"Tumor_Suppressor") in cnt ? cnt[tf,"Tumor_Suppressor"] : 0,
          (tf,"Driver")           in cnt ? cnt[tf,"Driver"]           : 0,
          total[tf]+0,
          tf in matched ? matched[tf] : 0,
          tf in na      ? na[tf]      : 0
}
' /data/najemd/TF_binding_cancer/results/integration/meth_expr_consistent_pairs_cancermine_labeled.tsv \
| { IFS= read -r header; echo "$header"; sort -t$'\t' -k5 -rn; } \
| column -t -s $'\t'

# TF    n_Oncogene  n_Tumor_Suppressor  n_Driver  n_any_labeled  n_matched_cancer  n_NA
# NRF1  52          49                  30        67             26                28
# BANP  14          14                  10        18             6                 7
# STILL NOT FINISHED IM HERE!!

# plot the distribution of CancerMine roles for the consistent genes in a stacked bar plot 
Rscript - << 'EOF'
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

dat_file <- "./results/integration/meth_expr_consistent_pairs_cancermine_labeled.tsv"
out_dir  <- "./results/integration"

dat <- fread(dat_file)

# ── Plot 1: stacked bar with counts and percentages ────────────────────────
role_summary <- dat[, .SD[1], by = .(TF, gene)][
  , .(TF, gene, role = cancermine_role, same = cancermine_role_same_cancer)]

role_summary[, category := fcase(
  role == "NA" | role == "" | is.na(role), "Unlabeled",
  grepl("Oncogene", role) & !is.na(same) & same != "NA", "Oncogene (same cancer)",
  grepl("Tumor_Suppressor", role) & !is.na(same) & same != "NA", "TSG (same cancer)",
  grepl("Oncogene", role), "Oncogene (other cancer)",
  grepl("Tumor_Suppressor", role), "TSG (other cancer)",
  default = "Driver only"
)]
role_summary[, category := factor(category, levels = c(
  "Oncogene (same cancer)", "TSG (same cancer)",
  "Oncogene (other cancer)", "TSG (other cancer)",
  "Driver only", "Unlabeled"))]

bar_dat <- role_summary[, .N, by = .(TF, category)]
bar_dat[, total := sum(N), by = TF]
bar_dat[, pct   := round(100 * N / total, 1)]

cols <- c(
  "Oncogene (same cancer)"  = "#d73027",
  "TSG (same cancer)"       = "#4575b4",
  "Oncogene (other cancer)" = "#f4a582",
  "TSG (other cancer)"      = "#abd9e9",
  "Driver only"             = "#fee090",
  "Unlabeled"               = "grey85")

p1 <- ggplot(bar_dat, aes(x = TF, y = N, fill = category)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = N),
            position = position_stack(vjust = 0.5),
            size = 3.5, colour = "grey20", fontface = "bold") +
  scale_fill_manual(values = cols, name = "CancerMine role") +
  labs(title = "CancerMine roles of consistent genes",
       x = NULL, y = "Number of unique genes") +
  theme_bw(base_size = 13)

ggsave(file.path(out_dir, "cancermine_stacked_bar.png"), p1,
       width = 6, height = 5, dpi = 150)
cat("[DONE] cancermine_stacked_bar.png\n")
EOF
# plot the lollipop plot of log2 fold change for the consistent genes, colored by whether they are labeled as oncogenes or tumor suppressors in the CancerMine database
Rscript - << 'EOF'
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

dat_file <- "./results/integration/meth_expr_consistent_pairs_cancermine_labeled.tsv"
out_dir  <- "./results/integration"
dat <- fread(dat_file)

# ── Build one row per gene-cancer pair ─────────────────────────────────────
anno <- unique(dat[, .(TF, gene, cancer,
                        log2FC = log2foldchange,
                        role   = cancermine_role,
                        same   = cancermine_role_same_cancer)])

anno[, og_same   := grepl("Oncogene",         same) & !is.na(same) & same != "NA"]
anno[, tsg_same  := grepl("Tumor_Suppressor", same) & !is.na(same) & same != "NA"]
anno[, og_other  := grepl("Oncogene",         role) & !(og_same)]
anno[, tsg_other := grepl("Tumor_Suppressor", role) & !(tsg_same)]

anno <- anno[order(TF, log2FC)]
anno[, row_id := paste0(gene, " (", cancer, ")")]
anno[, row_id := factor(row_id, levels = unique(row_id))]

# ── Compute percentages per TF ─────────────────────────────────────────────
anno[, concordant_same  := (og_same  & log2FC > 0) | (tsg_same  & log2FC < 0)]
anno[, concordant_other := !concordant_same &
                           ((og_other & log2FC > 0) | (tsg_other & log2FC < 0))]

pct_summary <- anno[, .(
  pct_same  = round(100 * sum(concordant_same)  / .N, 1),
  pct_other = round(100 * sum(concordant_other) / .N, 1),
  pct_rest  = round(100 * sum(!concordant_same & !concordant_other) / .N, 1),
  x_pos     = min(log2FC) * 0.95
), by = TF]

pct_summary[, label := paste0(
  "Same cancer concordant: ",  pct_same,  "%\n",
  "Other cancer concordant: ", pct_other, "%\n",
  "Other: ",                   pct_rest,  "%")]

# ── Plot A: lollipop ───────────────────────────────────────────────────────
pA <- ggplot(anno, aes(y = row_id)) +
  geom_segment(aes(x = 0, xend = log2FC,
                   y = row_id, yend = row_id,
                   colour = log2FC > 0), linewidth = 0.6) +
  geom_point(aes(x = log2FC, colour = log2FC > 0), size = 2.5) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_text(data = pct_summary,
            aes(x = x_pos, y = Inf, label = label),
            inherit.aes = FALSE,
            hjust = 0, vjust = 1.3, size = 2.8,
            colour = "grey20", fontface = "italic") +
  scale_colour_manual(values = c("TRUE" = "#d73027", "FALSE" = "#4575b4"),
                      labels = c("TRUE" = "Up (OG-like)", "FALSE" = "Down (TSG-like)"),
                      name = "Expression") +
  facet_wrap(~TF, scales = "free_y", ncol = 1) +
  labs(
    title    = "Expression direction of consistent genes vs CancerMine role",
    subtitle = paste0(
      "Each row is one gene-cancer pair, sorted by log2 fold change (tumor vs healthy).\n",
      "Lollipop direction: red = upregulated (oncogene-like), blue = downregulated (TSG-like).\n",
      "Tiles: dark red = OG in same cancer | dark blue = TSG in same cancer | ",
      "light red = OG in other cancer | light blue = TSG in other cancer | grey = not in CancerMine.\n",
      "Concordant = expression direction matches the CancerMine label. ",
      "Same-cancer genes are excluded from the other-cancer percentage."),
    x = "log2 fold change (tumor / healthy)",
    y = NULL) +
  theme_bw(base_size = 10) +
  theme(strip.text         = element_text(face = "bold"),
        legend.position    = "bottom",
        panel.grid.major.y = element_line(colour = "grey93"),
        plot.subtitle      = element_text(size = 7.5, colour = "grey30"))

# ── Plot B: annotation tiles ───────────────────────────────────────────────
tile_vars   <- c("og_same", "tsg_same", "og_other", "tsg_other")
tile_labels <- c(og_same   = "OG\n(same)",  tsg_same  = "TSG\n(same)",
                 og_other  = "OG\n(other)", tsg_other = "TSG\n(other)")
tile_cols   <- c(og_same   = "#d73027", tsg_same  = "#4575b4",
                 og_other  = "#f4a582", tsg_other = "#abd9e9")

tile_long <- melt(anno[, c("TF", "row_id", tile_vars), with = FALSE],
                  id.vars = c("TF", "row_id"),
                  variable.name = "anno_type", value.name = "present")
tile_long[, anno_type := factor(anno_type, levels = tile_vars)]

pB <- ggplot(tile_long, aes(x = anno_type, y = row_id,
                             fill = ifelse(present, as.character(anno_type), NA))) +
  geom_tile(colour = "white", linewidth = 0.4) +
  scale_fill_manual(values = tile_cols, na.value = "grey93", guide = "none") +
  scale_x_discrete(labels = tile_labels) +
  facet_wrap(~TF, scales = "free_y", ncol = 1) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 10) +
  theme(axis.text.y      = element_blank(),
        axis.ticks.y     = element_blank(),
        strip.text       = element_text(face = "bold", colour = "white"),
        strip.background = element_rect(fill = "white", colour = "white"),
        panel.grid       = element_blank())

# ── Combine and save ───────────────────────────────────────────────────────
p_combined <- pA + pB + plot_layout(widths = c(4, 1))

n_genes <- uniqueN(anno$gene)
h <- max(6, n_genes * 0.28)

ggsave(file.path(out_dir, "cancermine_lollipop_annotated.pdf"),
       p_combined, width = 10, height = h, limitsize = FALSE)
cat("[DONE] cancermine_lollipop_annotated.png\n")
EOF

# =============================================================================
# CancerMine Enrichment Control Test
# =============================================================================
# PURPOSE:
#   Tests whether the genes linked to significant NRF1/BANP probe-gene pairs are more enriched for known cancer genes (from CancerMine) than what would be expected by chance, using the full set of motif-linked genes as background.
#
# INPUT FILES:
#   1. probe_gene_pairs_in_motifs_with_tss.tsv.gz
#      → All probe-gene pairs where the probe overlaps an NRF1 or BANP motif.
#        This is the UNIVERSE of testable pairs before any significance filter.
#        Columns used: TF (stripped to NRF1/BANP), gene.
#
#   2. meth_expr_consistent_pairs_cancermine_labeled.tsv
#      → The SIGNIFICANT pairs: probe-gene pairs that passed all filters
#        (concordant methylation + expression change in TCGA + cell lines).
#        Columns used: TF, gene.
#
#   3. cancermine.tsv
#      → CancerMine database of literature-mined cancer gene roles.
#        Col 2 = role (Oncogene / Tumor_Suppressor / Driver)
#        Col 7 = gene symbol
#        Only Oncogene, Tumor_Suppressor and Driver roles are kept.
#
# ANALYSIS STEPS:
#
#   STEP 1 — Build the background set
#     Take all unique TF-gene pairs from the motif file (bg_all).
#     Filter to PROXIMAL probes only (is_proximal == TRUE) to match the
#     same probe scope used for the significant pairs.
#     Remove any pair that also appears in the significant list (anti-join).
#     Result: bg_only = genes linked to proximal NRF1/BANP motif probes
#     but NOT significant. This is the null expectation: what % of proximal
#     motif-linked genes are cancer genes?
#
#   STEP 2 — CancerMine labeling
#     For both the significant set and the background set, merge with CancerMine.
#     A gene is called "labeled" if it appears in CancerMine with ANY role
#     (Oncogene, Tumor Suppressor, or Driver) in ANY cancer type.
#     Multiple roles per gene are collapsed: one row per TF-gene pair with
#     boolean flags any_labeled, is_OG, is_TSG.
#
#   STEP 3 — Per-TF summary table
#     For each TF (NRF1, BANP) and each group (Significant, Background):
#     count total genes, count labeled genes, compute %.
#
#   STEP 4 — Fisher's exact test (per TF)
#     Builds a 2x2 contingency table:
#
#                  | CancerMine-labeled | Not labeled |
#     Significant  |        k           |    N - k    |
#     Background   |        j           |    M - j    |
#
#     Fisher's exact test asks: is the proportion k/N significantly higher
#     than j/M? (one-sided enrichment in the significant set)
#     Output: odds ratio (OR), 95% CI, p-value.
#     OR > 1 means significant genes are more cancer-enriched than background.
#     Saved to: cancermine_enrichment_fisher.tsv
#
#   STEP 5 — Plot
#     Bar chart per TF showing % labeled in Background (grey) vs Significant
#     (blue), with n shown on each bar.
#     Above each panel: significance stars, exact p-value, OR [95% CI].
#     Subtitle and caption clarify that background is restricted to
#     proximal probes (is_proximal == TRUE) only.
#     Saved to: cancermine_enrichment_vs_background.pdf
#
# OUTPUT FILES:
#   - cancermine_enrichment_fisher.tsv  → table with OR, CI, p-value per TF
#   - cancermine_enrichment_vs_background.pdf → bar plot
# =============================================================================

Rscript - << 'EOF'
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

project_root <- "/data/najemd/TF_binding_cancer"
out_dir      <- file.path(project_root, "results/integration")

# ── Load files ───────────────────────────────────────────────────────────────
bg_raw <- fread(cmd = paste("zcat",
  file.path(project_root, "methylation/probe_gene_pairs_in_motifs_with_tss.tsv.gz")))
bg_raw[, TF := sub("_.*", "", TF)]
cat("[INFO] Background pairs:", nrow(bg_raw), "\n")

sig <- fread(file.path(out_dir, "meth_expr_consistent_pairs_cancermine_labeled.tsv"))
cat("[INFO] Significant pairs:", nrow(sig), "\n")

# CancerMine: rename by position (col 2 = role, col 7 = gene symbol)
cm <- fread(file.path(out_dir, "cancermine.tsv"))
cat("[INFO] CancerMine columns:", paste(names(cm), collapse = ", "), "\n")
setnames(cm, c(2, 7), c("role", "gene_symbol"))
cm <- cm[role %in% c("Oncogene", "Tumor_Suppressor", "Driver")]
cm_labeled <- unique(cm[, .(gene = gene_symbol, role)])

# ── Build background-only gene sets ──────────────────────────────────────────
sig_pairs <- unique(sig[, .(TF, gene)])
bg_raw <- bg_raw[is_proximal == TRUE]
bg_all    <- unique(bg_raw[, .(TF, gene)])
bg_only   <- bg_all[!sig_pairs, on = .(TF, gene)]
cat("[INFO] Sig pairs:", nrow(sig_pairs), " | BG all:", nrow(bg_all),
    " | BG only:", nrow(bg_only), "\n")

# ── CancerMine labeling function ──────────────────────────────────────────────
label_genes <- function(dt) {
  dt <- unique(dt[, .(TF, gene)])
  dt <- merge(dt, cm_labeled, by = "gene", all.x = TRUE, allow.cartesian = TRUE)
  dt[, .(any_labeled = any(!is.na(role)),
         is_OG       = any(role == "Oncogene",        na.rm = TRUE),
         is_TSG      = any(role == "Tumor_Suppressor", na.rm = TRUE)),
     by = .(TF, gene)]
}

sig_lab <- label_genes(sig_pairs)
bg_lab  <- label_genes(bg_only)

# ── Per-TF summary ────────────────────────────────────────────────────────────
make_summary <- function(dt, group_name) {
  dt[, .(group       = group_name,
         n_genes     = .N,
         n_labeled   = sum(any_labeled),
         pct_labeled = round(100 * sum(any_labeled) / .N, 1),
         pct_OG      = round(100 * sum(is_OG)  / .N, 1),
         pct_TSG     = round(100 * sum(is_TSG) / .N, 1)),
     by = TF]
}

summ <- rbindlist(list(make_summary(sig_lab, "Significant"),
                       make_summary(bg_lab,  "Background")))
cat("\n── Summary ─────────────────────────────────────────────────────────\n")
print(summ[order(TF, group)])

# ── Fisher's exact test ───────────────────────────────────────────────────────
fisher_res <- rbindlist(lapply(unique(sig_lab$TF), function(tf) {
  s   <- sig_lab[TF == tf]
  b   <- bg_lab[TF == tf]
  mat <- matrix(c(sum(s$any_labeled),  sum(!s$any_labeled),
                  sum(b$any_labeled),  sum(!b$any_labeled)),
                nrow = 2, byrow = TRUE,
                dimnames = list(c("Significant","Background"),
                                c("Labeled","Not_labeled")))
  cat("\n[Fisher] TF =", tf, "\n"); print(mat)
  ft <- fisher.test(mat)
  data.table(TF        = tf,
             OR        = round(ft$estimate, 2),
             CI_low    = round(ft$conf.int[1], 2),
             CI_high   = round(ft$conf.int[2], 2),
             p_value   = signif(ft$p.value, 3),
             sig_label = ifelse(ft$p.value < 0.001, "***",
                         ifelse(ft$p.value < 0.01,  "**",
                         ifelse(ft$p.value < 0.05,  "*", "ns"))))
}))
cat("\n── Fisher results ──────────────────────────────────────────────────\n")
print(fisher_res)
fwrite(fisher_res, file.path(out_dir, "cancermine_enrichment_fisher.tsv"), sep = "\t")

# ── Plot ──────────────────────────────────────────────────────────────────────
plot_dat <- merge(summ,
                  fisher_res[, .(TF, sig_label, OR, CI_low, CI_high, p_value)],
                  by = "TF")
plot_dat[, group := factor(group, levels = c("Background", "Significant"))]

# Always show actual p-value (scientific for very small) + stars
plot_dat[, p_label := paste0(
  sig_label, "\n",                                             # stars first line
  ifelse(p_value < 0.001,
         paste0("p = ", formatC(p_value, format = "e", digits = 2)),
         paste0("p = ", formatC(p_value, format = "f", digits = 3))),
  "\nOR = ", OR, " [", CI_low, "–", CI_high, "]"
)]

# Fix annotation at a consistent height well above bars, only on Significant bar
star_dat <- plot_dat[group == "Significant", .(TF, y = 95, p_label)]

p <- ggplot(plot_dat, aes(x = group, y = pct_labeled, fill = group)) +
  geom_col(width = 0.55) +
  geom_text(aes(label = paste0(pct_labeled, "%\n(n=", n_genes, ")")),
            vjust = -0.4, size = 3.4) +
  geom_text(data = star_dat,
            aes(x = 1.5, y = y, label = p_label),   # centered between both bars
            inherit.aes = FALSE, size = 3.6, fontface = "bold", hjust = 0.5) +
  scale_fill_manual(values = c(Background = "grey70", Significant = "#4575b4"),
                    guide = "none") +
  facet_wrap(~TF) +
  scale_y_continuous(limits = c(0, 105),
                     breaks = seq(0, 100, 25),
                     labels = function(x) paste0(x, "%")) +
  labs(
    title    = "CancerMine enrichment: significant vs background genes",
    subtitle = paste0(
      "Fisher's exact test: 2\u00d72 table comparing the proportion of CancerMine-labeled genes\n",
      "in the SIGNIFICANT set (probe-gene pairs passing all filters) vs the BACKGROUND\n",
      "(all proximal motif-linked probe-gene pairs for NRF1/BANP that did NOT pass the filters)."
    ),
    x       = NULL,
    y       = "% genes with any CancerMine label",
    caption = paste0(
      "Background = every unique gene linked to a PROXIMAL motif probe in the probe_gene_pairs_in_motifs_with_tss file, EXCLUDING genes already in the significant list.\n",
      "A gene is 'labeled' if CancerMine lists it as an Oncogene, Tumor Suppressor, or Driver in any cancer type.\n",
      "Odds ratio > 1 means significant genes are more enriched for cancer genes than expected by chance."
    )
  ) +
  theme_bw(base_size = 13) +
  theme(strip.text            = element_text(face = "bold"),
        plot.subtitle         = element_text(size = 8.5, colour = "grey20", lineheight = 1.3),
        plot.caption          = element_text(size = 7.5, colour = "grey40", hjust = 0,
                                             lineheight = 1.3),
        plot.caption.position = "plot")

ggsave(file.path(out_dir, "cancermine_enrichment_vs_background.pdf"), p,
       width = 9, height = 7)
cat("[DONE] cancermine_enrichment_vs_background.pdf\n")
cat("[DONE] cancermine_enrichment_fisher.tsv\n")
EOF


###################################################################
# Methylation landscape plots — target probe vs control probes
###################################################################
# Pipeline:
#   Step 1  — Load significant pairs and parse motif coordinates
#   Step 2  — Find nearby control probes (TSS ± 2000 bp, outside motif)
#   Step 3  — Load beta values (target probes from beta matrices,
#              control probes from raw methylation files)
#   Step 4  — Compute delta beta (mean tumor − mean healthy) per probe
#   Step 5  — Wilcoxon test: target delta beta vs control probe delta betas
#              (one-sample, mu = target delta beta; BH correction)
#   Step 6  — Generate landscape plots (methylation profile + expression boxplot)
#
# Input:
#   ./results/integration/meth_expr_consistent_pairs.tsv
#   ./methylation/annotated_methylation_data_probes_filtered.bed
#   ./results/multi_omics/expr_meth_shared_samples.txt
#   ./expression/gene_expression_matrix_protein_coding.tsv
#   ./methylation/beta_matrices_<TF>/
#   ./methylation/filtered_methylation/
#   /data/genome/annotations/hg38_tss.bed
#
# Output:
#   ./results/integration/landscape_plots_<TF>/landscape_<TF>_<gene>_<cancer>.png
#   ./results/integration/wilcox_delta_beta_<TF>.tsv
#
# ══════════════════════════════════════════════════════════════════════════════
Rscript - <<'EOF'
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(parallel)
})

# ── Helpers ───────────────────────────────────────────────────────────────────
fmt_p <- function(x) {
  if (is.na(x)) return("NA")
  if (x < 2.2e-16) return("<2.2e-16")
  formatC(x, format = "e", digits = 2)
}
p_label <- function(p) {
  stars <- if (is.na(p)) "NA" else if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"
  paste0("p = ", fmt_p(p), "\n", stars)
}

# ── Paths ─────────────────────────────────────────────────────────────────────
pairs_file  <- "./results/integration/meth_expr_consistent_pairs.tsv"
bed_file    <- "./methylation/annotated_methylation_data_probes_filtered.bed"
shared_file <- "./results/multi_omics/expr_meth_shared_samples.txt"
tpm_file    <- "./expression/gene_expression_matrix_protein_coding.tsv"
meth_dir    <- "./methylation/filtered_methylation"
tss_file    <- "/data/genome/annotations/hg38_tss.bed"
out_dir     <- "./results/integration/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (tf_name in c("NRF1", "BANP")) {
  cat(sprintf("\n========== TF: %s ==========\n", tf_name))

  beta_dir      <- paste0("./methylation/beta_matrices_", tf_name)
  landscape_dir <- file.path(out_dir, paste0("landscape_plots_", tf_name))
  dir.create(landscape_dir, recursive = TRUE, showWarnings = FALSE)

  # ── STEP 1 — Load pairs and parse motif coordinates ────────────────────────
  cat("[STEP 1] Loading pairs...\n")
  pairs <- fread(pairs_file)
  pairs <- pairs[TF == tf_name]
  pairs[, motif_chr    := sub(":.*", "", motif_instance_id)]
  pairs[, motif_start  := as.integer(sub("^[^:]+:([0-9]+):.*",         "\\1", motif_instance_id))]
  pairs[, motif_end    := as.integer(sub("^[^:]+:[0-9]+:([0-9]+):.*",  "\\1", motif_instance_id))]
  pairs[, motif_strand := sub(".*:", "", motif_instance_id)]
  pairs[, motif_mid    := as.integer((motif_start + motif_end) / 2)]
  cat(sprintf("  %d pairs loaded\n", nrow(pairs)))

  # ── STEP 2 — Find control probes ───────────────────────────────────────────
  cat("[STEP 2] Finding control probes...\n")
  probe_bed <- fread(bed_file, header = FALSE,
                     col.names = c("probe_chr", "probe_start", "probe_end", "probe"))
  probe_bed[, probe_start := as.integer(probe_start)]
  probe_bed[, probe_end   := as.integer(probe_end)]

  # Load TSS annotation
  tss_bed <- fread(tss_file, header = FALSE,
                   col.names = c("tss_chr", "tss_start", "tss_end",
                                 "gene_name", "ensembl_gene", "ensembl_tx", "strand"))

  # For each pair, find the TSS closest to the motif midpoint also saves gene_strand
  pairs[, c("tss_pos", "gene_strand") := {
    gn  <- gene[1]
    chr <- motif_chr[1]
    mid <- motif_mid[1]
    candidates <- tss_bed[gene_name == gn & tss_chr == chr]
    if (nrow(candidates) == 0) list(NA_integer_, NA_character_)
    else {
      best <- candidates[which.min(abs(tss_start - mid))]
      list(best$tss_start, best$strand)
    }
  }, by = .(gene, motif_chr, motif_mid)]

  n_missing <- pairs[is.na(tss_pos), .N]
  if (n_missing > 0)
    cat(sprintf("  WARNING: %d pairs have no TSS match\n", n_missing))

  # Window: TSS ± 2000 bp
  pairs[, win_start := pmax(0L, tss_pos - 2000L)]
  pairs[, win_end   :=          tss_pos + 2000L  ]

  ctrl_list <- lapply(seq_len(nrow(pairs)), function(i) {
    row <- pairs[i]
    if (is.na(row$tss_pos)) return(NULL)
    nearby <- probe_bed[
      probe_chr   == row$motif_chr   &
      probe_start >= row$win_start   &
      probe_end   <= row$win_end     &
      probe       != row$probe       &
      !(probe_start < row$motif_end  &
        probe_end   > row$motif_start)
    ]
    if (nrow(nearby) == 0) return(NULL)
    nearby[, `:=`(gene         = row$gene,
                  cancer       = row$cancer,
                  target_probe = row$probe)]
    nearby[]
  })
  controls <- rbindlist(ctrl_list, fill = TRUE)
  cat(sprintf("  %d control probe rows\n", nrow(controls)))

  # ── STEP 3 — Load beta values ───────────────────────────────────────────────
  cat("[STEP 3] Loading beta values...\n")
  shared         <- fread(shared_file, header = FALSE)$V1
  cancers_needed <- unique(pairs$cancer)
  target_probes  <- unique(pairs$probe)
  control_probes <- unique(controls$probe)

  # Target probes — from beta matrices
  target_beta <- rbindlist(lapply(cancers_needed, function(cx) {
    bf <- list.files(beta_dir, pattern = paste0(cx, ".*\\.tsv\\.gz$"), full.names = TRUE)
    if (length(bf) == 0) return(NULL)
    mat <- fread(bf[1])
    if (!"probe" %in% colnames(mat)) setnames(mat, 1, "probe")
    mat <- mat[probe %in% target_probes]
    if (nrow(mat) == 0) return(NULL)
    long <- melt(mat, id.vars = "probe", variable.name = "col_name", value.name = "meth_beta")
    long[, col_name     := as.character(col_name)]
    long[, full_patient := sub("_(tumor|healthy)$", "", col_name)]
    long[, sample_type  := sub(".*_(tumor|healthy)$", "\\1", col_name)]
    long <- long[sample_type %in% c("tumor","healthy") & full_patient %in% shared]
    if (nrow(long) == 0) return(NULL)
    long[, cancer   := cx]
    long[, meth_pct := if (max(meth_beta, na.rm = TRUE) <= 1.5) meth_beta * 100 else meth_beta]
    long[, .(cancer, probe, sample_type, meth_pct)]
  }), fill = TRUE)

  # Control probes — from raw methylation files
  raw_files <- list.files(meth_dir,
    pattern = "(01A_1|11A_1)_annotated_methylation_filtered\\.bed\\.gz$",
    full.names = TRUE)
  raw_meta <- data.table(file = raw_files, fname = basename(raw_files))
  raw_meta[, cancer       := sub("^HM450_TCGA-([A-Z]+)-.*", "\\1", fname)]
  raw_meta[, full_patient := sub("^HM450_TCGA-([A-Z]+-TCGA-[A-Z0-9]+-[A-Z0-9]+)-.*", "\\1", fname)]
  raw_meta[, sample_type  := fifelse(grepl("-01A_1_", fname), "tumor",
                             fifelse(grepl("-11A_1_", fname), "healthy", NA_character_))]
  raw_meta <- raw_meta[cancer %in% cancers_needed & !is.na(sample_type) & full_patient %in% shared]

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
    dt[, `:=`(cancer      = raw_meta$cancer[i],
              sample_type = raw_meta$sample_type[i])]
    dt[, meth_pct := if (max(meth_beta, na.rm = TRUE) <= 1.5) meth_beta * 100 else meth_beta]
    dt[, .(cancer, probe, sample_type, meth_pct)]
  }

  cl <- makeCluster(10L)
  clusterExport(cl, c("raw_meta", "control_probes"), envir = environment())
  clusterEvalQ(cl, library(data.table))
  control_beta <- rbindlist(parLapply(cl, seq_len(nrow(raw_meta)), read_one), fill = TRUE)
  stopCluster(cl)

  beta_long <- rbindlist(list(target_beta, control_beta), fill = TRUE)
  cat(sprintf("  Total beta rows: %d\n", nrow(beta_long)))

  # ── STEP 4 — Delta beta per probe ──────────────────────────────────────────
  cat("[STEP 4] Computing delta beta...\n")
  all_probe_rows <- rbind(
    pairs[,    .(gene, cancer, probe, probe_type = "target")],
    controls[, .(gene, cancer, probe, probe_type = "control")]
  )
  all_probe_rows <- merge(all_probe_rows, beta_long,
                          by = c("cancer", "probe"), allow.cartesian = TRUE)
  all_probe_rows <- all_probe_rows[sample_type %in% c("tumor","healthy") & is.finite(meth_pct)]

  delta_db <- all_probe_rows[,
    .(mean_beta_tumor   = mean(meth_pct[sample_type == "tumor"],   na.rm = TRUE),
      mean_beta_healthy = mean(meth_pct[sample_type == "healthy"], na.rm = TRUE)),
    by = .(gene, cancer, probe, probe_type)
  ][, delta_beta := mean_beta_tumor - mean_beta_healthy]

  # ── STEP 5 — Wilcoxon: target delta beta vs control delta betas ────────────
  cat("[STEP 5] Wilcoxon test (target probe vs control probes)...\n")
  ctrl_map <- unique(controls[, .(gene, cancer, probe, target_probe)])

  ctrl_db <- merge(
    delta_db[probe_type == "control" & is.finite(delta_beta),
             .(gene, cancer, probe, ctrl_delta_beta = delta_beta)],
    ctrl_map, by = c("gene", "cancer", "probe")
  )
  target_db <- delta_db[probe_type == "target" & is.finite(delta_beta),
                        .(gene, cancer, probe, target_delta_beta = delta_beta)]

  wilcox_results <- merge(target_db, ctrl_db,
                          by.x = c("gene", "cancer", "probe"),
                          by.y = c("gene", "cancer", "target_probe"),
                          allow.cartesian = TRUE)[, {
    tdb   <- target_delta_beta[1]
    cdbs  <- ctrl_delta_beta
    wt <- tryCatch(wilcox.test(cdbs, mu = tdb, exact = FALSE),
               error = function(e) NULL)
    .(target_delta_beta    = tdb,
      mean_ctrl_delta_beta = mean(cdbs, na.rm = TRUE),
      n_ctrl_probes        = length(cdbs),
      delta_diff           = tdb - mean(cdbs, na.rm = TRUE),
      p_wilcox             = if (!is.null(wt)) wt$p.value else NA_real_)
  }, by = .(gene, cancer, probe)]

  wilcox_results[, p_adj := p.adjust(p_wilcox, method = "BH")]

  wilcox_results <- merge(wilcox_results, pairs, by = c("gene", "cancer", "probe"))

  fwrite(wilcox_results,
         file.path(out_dir, paste0("wilcox_delta_beta_", tf_name, ".tsv")),
         sep = "\t")
  cat(sprintf("  Significant (p_adj < 0.05): %d / %d pairs\n",
              wilcox_results[p_adj < 0.05, .N], nrow(wilcox_results)))

  # ── STEP 6 — Landscape plots ───────────────────────────────────────────────
  cat("[STEP 6] Loading TPM matrix...\n")
  tpm_mat <- fread(tpm_file, header = TRUE, check.names = FALSE)
  setnames(tpm_mat, 1L, "tpm_sample_id")
  tpm_mat[, cancer_col := sub("^TCGA-([A-Z]+)-.*", "\\1", tpm_sample_id)]
  tpm_mat[, st         := fifelse(grepl("-01A_1", tpm_sample_id), "Tumor",
                          fifelse(grepl("-11A_1", tpm_sample_id), "Healthy", "other"))]
  tpm_mat <- tpm_mat[st != "other"]

  pairs_to_plot <- unique(wilcox_results[, .(gene, cancer, probe)])
  cat(sprintf("[STEP 6] %d landscape plots to generate\n", nrow(pairs_to_plot)))

  for (i in seq_len(nrow(pairs_to_plot))) {
    gn  <- pairs_to_plot$gene[i]
    cx  <- pairs_to_plot$cancer[i]
    tpr <- pairs_to_plot$probe[i]

    cpr    <- ctrl_map[gene == gn & cancer == cx & target_probe == tpr, probe]
    all_pr <- c(tpr, cpr)

    motif_row <- pairs[gene == gn & cancer == cx & probe == tpr][1]
    if (nrow(motif_row) == 0) next

    pos   <- probe_bed[probe %in% all_pr, .(probe, probe_start, probe_chr)]
    betas <- beta_long[probe %in% all_pr & cancer == cx & is.finite(meth_pct)]
    if (nrow(betas) == 0) next

    med <- betas[, .(median_meth = median(meth_pct, na.rm = TRUE),
                 sd_meth     = sd(meth_pct,     na.rm = TRUE)),
             by = .(probe, sample_type)]
    med <- merge(med, pos[, .(probe, probe_start)], by = "probe")
    med[, probe_type := fifelse(probe == tpr, "target", "control")]
    setorder(med, probe_start)

    wx <- wilcox_results[gene == gn & cancer == cx & probe == tpr]

    ctrl_subtitle <- if (nrow(wx) > 0)
      sprintf(
        "Target probe ΔBeta = %.1f%%  |  Mean control probes ΔBeta = %.1f%%  |  diff = %.1f%%  |  Wilcoxon p_adj = %.3g  (%d control probes)",
        wx$target_delta_beta[1],
        wx$mean_ctrl_delta_beta[1],
        wx$delta_diff[1],
        wx$p_adj[1],
        wx$n_ctrl_probes[1]
      ) else ""

    p <- ggplot(med, aes(x = probe_start, y = median_meth,
                         colour = sample_type, group = sample_type)) +
      # Motif region shading
      annotate("rect",
               xmin = motif_row$motif_start, xmax = motif_row$motif_end,
               ymin = -Inf, ymax = Inf, fill = "#f4a261", alpha = 0.25) +
      annotate("text",
               x     = (motif_row$motif_start + motif_row$motif_end) / 2,
               y     = Inf,
               label = sprintf("%s motif", tf_name),
               vjust = 1.4, size = 2.8, colour = "#e76f51", fontface = "italic") +
      # TSS vertical line
      geom_vline(xintercept = motif_row$tss_pos,
                 linetype = "solid", colour = "#2d6a4f", linewidth = 0.7) +
      annotate("text",
               x     = motif_row$tss_pos,
               y     = Inf,
               label = "TSS",
               vjust = 1.4, hjust = -0.15,
               size  = 2.8, colour = "#2d6a4f", fontface = "bold") +
            # Gene orientation arrow — 300 bp arm from TSS in transcription direction
      annotate("segment",
               x    = motif_row$tss_pos,
               xend = motif_row$tss_pos +
                        fifelse(motif_row$gene_strand == "+", 300L, -300L),
               y = 4, yend = 4,
               colour = "#2d6a4f", linewidth = 1.2,
               arrow = arrow(length = unit(0.25, "cm"),
                             type = "closed", ends = "last")) +
      annotate("text",
               x     = motif_row$tss_pos +
                          fifelse(motif_row$gene_strand == "+", 320L, -320L),
               y     = 4,
               label = fifelse(motif_row$gene_strand == "+", "3′ →", "← 3′"),
               hjust = fifelse(motif_row$gene_strand == "+", 0, 1),
               size  = 3, colour = "#2d6a4f", fontface = "bold") +
      geom_ribbon(aes(ymin = median_meth - sd_meth,
                ymax = median_meth + sd_meth,
                fill = sample_type),
            alpha = 0.15, colour = NA) +
      geom_line(linewidth = 0.55, alpha = 0.8) +
      geom_point(data = med[probe_type == "control"],
                 aes(fill = sample_type),
                 shape = 21, size = 2.5, colour = "white", stroke = 0.4) +
      geom_point(data = med[probe_type == "target"],
                 aes(fill = sample_type),
                 shape = 23, size = 4.5, colour = "black", stroke = 0.8) +
      geom_vline(xintercept = med[probe_type == "target", probe_start[1]],
                 linetype = "dashed", colour = "grey30", linewidth = 0.5) +
            # ── 100 bp scale bar — bottom-right corner ──────────────────────────
      annotate("segment",                          # horizontal bar
               x = motif_row$tss_pos + 1700L, xend = motif_row$tss_pos + 1800L,
               y = 8, yend = 8,
               colour = "grey20", linewidth = 1) +
      annotate("segment",                          # left tick
               x = motif_row$tss_pos + 1700L, xend = motif_row$tss_pos + 1700L,
               y = 6, yend = 10,
               colour = "grey20", linewidth = 1) +
      annotate("segment",                          # right tick
               x = motif_row$tss_pos + 1800L, xend = motif_row$tss_pos + 1800L,
               y = 6, yend = 10,
               colour = "grey20", linewidth = 1) +
      annotate("text",
               x     = motif_row$tss_pos + 1750L,
               y     = 12,
               label = "100 bp",
               size = 3, colour = "grey20", hjust = 0.5) +
      scale_colour_manual(values = c(tumor = "#d73027", healthy = "#4575b4"), name = NULL) +
      scale_fill_manual(values   = c(tumor = "#d73027", healthy = "#4575b4"), name = NULL) +
      coord_cartesian(xlim = c(motif_row$tss_pos - 2000L, motif_row$tss_pos + 2000L)) +
      labs(
        title    = sprintf("%s  |  %s  |  %s  (◆ = target probe  ●  = control probes)",
                           gn, cx, tf_name),
        subtitle = ctrl_subtitle,
        x        = sprintf("Genomic position (%s)", motif_row$motif_chr),
        y        = "Median methylation (%)") +
      theme_bw(base_size = 12) +
      theme(legend.position  = "top",
            plot.subtitle    = element_text(size = 8, colour = "grey20", face = "bold"),
            panel.grid.minor = element_blank())

    # Expression boxplot
    expr_panel <- NULL
    if (gn %in% names(tpm_mat)) {
      expr_dt <- tpm_mat[cancer_col == cx,
                         .(sample_type = st, log2_expr = log2(get(gn) + 1))]
      expr_dt <- expr_dt[sample_type %in% c("Tumor", "Healthy")]
      expr_dt[, sample_type := factor(sample_type, levels = c("Healthy", "Tumor"))]

      if (nrow(expr_dt) >= 4) {
        expr_subtitle <- sprintf("log2FC = %.2f  |  padj = %.3g",
                                 wx$log2foldchange[1], wx$padj[1])
        expr_min  <- min(expr_dt$log2_expr, na.rm = TRUE)
        expr_max  <- max(expr_dt$log2_expr, na.rm = TRUE)
        expr_rng  <- max(expr_max - expr_min, 1)
        expr_ylim <- c(max(0, expr_min - 0.08 * expr_rng), expr_max + 0.28 * expr_rng)
        med_h     <- wx$median_log2tpm_healthy[1]
        med_t     <- wx$median_log2tpm_tumor[1]

        expr_panel <- ggplot(expr_dt, aes(x = sample_type, y = log2_expr)) +
          geom_boxplot(width = 0.5, outlier.shape = NA, color = "#0b3c5d", fill = NA) +
          geom_jitter(data = expr_dt[sample_type == "Tumor"],
                      color = "#d73027", width = 0.12, size = 1.8, alpha = 0.65) +
          geom_jitter(data = expr_dt[sample_type == "Healthy"],
                      color = "#4575b4", width = 0.12, size = 1.8, alpha = 0.65) +
          annotate("text", x = 1.35, y = med_h,
                   label = sprintf("%.2f", med_h), size = 3.2, colour = "#4575b4", hjust = 0) +
          annotate("text", x = 2.35, y = med_t,
                   label = sprintf("%.2f", med_t), size = 3.2, colour = "#d73027", hjust = 0) +
          scale_y_continuous(limits = expr_ylim) +
          scale_x_discrete(expand = expansion(add = c(0.5, 1.2))) +
          labs(title    = "Expression",
               subtitle = expr_subtitle,
               x = NULL, y = "log2(TPM+1)") +
          theme_bw(base_size = 11) +
          theme(legend.position  = "none",
                plot.subtitle    = element_text(size = 8, colour = "grey20", face = "bold"),
                plot.title       = element_text(size = 10))
      }
    }

    combined   <- if (!is.null(expr_panel)) (p | expr_panel) + plot_layout(widths = c(3, 1)) else p
    plot_width <- if (!is.null(expr_panel)) 13 else 9
    ggsave(file.path(landscape_dir,
                     sprintf("landscape_%s_%s_%s.png", tf_name, gn, cx)),
           combined, width = plot_width, height = 5, dpi = 150)

    if (i %% 20 == 0 || i == nrow(pairs_to_plot))
      cat(sprintf("  %d / %d plots done\n", i, nrow(pairs_to_plot)))
  }

  cat(sprintf("[DONE] %s — %d plots saved to %s/\n",
              tf_name, nrow(pairs_to_plot), landscape_dir))
}
EOF

#######################################################################################
# Control probe analysis — target vs nearby control probe DELTA BETA
#
# Pipeline:
#   Step 1  — Load filtered pairs and motif instance IDs
#   Step 2  — Find nearby control probes (TSS ± 2000 bp, outside motif)
#   Step 3  — Load beta values (target probes from beta matrices,
#              control probes from raw methylation files)
#   Step 4  — Compute delta beta (mean tumor − mean healthy) per probe
#   Step 5  — Wilcoxon test: target vs control delta betas (BH correction)
#   Step 6  — Scatter plot: target delta beta vs mean control delta beta
#   Step 7  — Filter significant pairs (p_adj < 0.05) and save
#
# Input:
#   ./results/integration/meth_expr_filtered_pairs.tsv
#   ./methylation/annotated_methylation_data_probes_filtered.bed
#   ./methylation/beta_matrices_<TF>/
#   ./methylation/filtered_methylation/
#   ./results/multi_omics/expr_meth_shared_samples.txt
#   /data/genome/annotations/hg38_tss.bed
#
# Output:
#   ./results/integration/controle_meth_probe/wilcox_delta_beta_<TF>.tsv
#   ./results/integration/controle_meth_probe/plot_target_vs_ctrl_delta_beta_<TF>.png
#   ./results/integration/controle_meth_probe/filtered_pairs_delta_beta_<TF>.tsv
#
#######################################################################################

Rscript - <<'EOF'
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(parallel)
})

# ── Helpers ───────────────────────────────────────────────────────────────────
fmt_p <- function(x) {
  if (is.na(x)) return("NA")
  if (x < 2.2e-16) return("<2.2e-16")
  formatC(x, format = "e", digits = 2)
}
p_label <- function(p) {
  stars <- if (is.na(p)) "NA" else if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"
  paste0("p = ", fmt_p(p), "\n", stars)
}

# ──  paths & parameters ─────────────────────────────────────────────────
pairs_file  <- "./results/integration/meth_expr_consistent_pairs.tsv"
motif_file  <- "./methylation/probe_gene_pairs_in_motifs_with_tss.tsv.gz"
meth_dir    <- "./methylation/filtered_methylation"
shared_file <- "./results/multi_omics/expr_meth_shared_samples.txt"
out_dir     <- "./results/integration/controle_meth_probe"
thresh      <- 10   # delta beta difference threshold (%)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (tf_name in c("NRF1", "BANP")) {
  cat(sprintf("\n========== TF: %s ==========\n", tf_name))

  beta_dir <- paste0("./methylation/beta_matrices_", tf_name)
  ckpt_dir <- file.path(out_dir, paste0("checkpoints_", tf_name))
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
    # ── STEP 1 — Load consistent pairs and parse motif coordinates ─────────────
  todo <- load_or_run("step1_todo", {
    cat("[STEP 1] Loading consistent pairs...\n")
    pairs <- fread(pairs_file)
    pairs <- pairs[TF == tf_name]
    pairs[, motif_chr    := sub(":.*", "", motif_instance_id)]
    pairs[, motif_start  := as.integer(sub("^[^:]+:([0-9]+):.*",        "\\1", motif_instance_id))]
    pairs[, motif_end    := as.integer(sub("^[^:]+:[0-9]+:([0-9]+):.*", "\\1", motif_instance_id))]
    pairs[, motif_strand := sub(".*:", "", motif_instance_id)]
    cat(sprintf("  %d pairs loaded\n", nrow(pairs)))
    pairs
  })
  todo[, motif_mid := as.integer((motif_start + motif_end) / 2L)]
  cat(sprintf("[STEP 1] %d rows\n", nrow(todo)))
  print(todo[1:min(5, .N), .(gene, cancer, probe, motif_chr, motif_start, motif_end, motif_strand)])
  # ── STEP 2 — Find nearby control probes ────────────────────────────────────
  controls <- load_or_run("step2_controls", {
    cat("[STEP 2] Loading probe BED file...\n")
    all_probes <- fread("./methylation/annotated_methylation_data_probes_filtered.bed",
                        header = FALSE,
                        col.names = c("probe_chr", "probe_start", "probe_end", "probe"))
    all_probes[, probe_start := as.integer(probe_start)]
    all_probes[, probe_end   := as.integer(probe_end)]
    cat(sprintf("  %d probes loaded\n", nrow(all_probes)))

    tss_bed <- fread("/data/genome/annotations/hg38_tss.bed", header = FALSE,
                 col.names = c("tss_chr", "tss_start", "tss_end",
                               "gene_name", "ensembl_gene", "ensembl_tx", "strand"))

    todo[, tss_pos := {
      gn  <- gene[1]; chr <- motif_chr[1]; mid <- motif_mid[1]
      candidates <- tss_bed[gene_name == gn & tss_chr == chr]
      if (nrow(candidates) == 0L) NA_integer_
      else candidates[which.min(abs(tss_start - mid)), tss_start]
    }, by = .(gene, motif_chr, motif_mid)]

    todo <- todo[!is.na(tss_pos)]   # drop any pair with no TSS match
    todo[, win_start := pmax(0L, tss_pos - 2000L)]
    todo[, win_end   :=           tss_pos + 2000L ]

    cat("[STEP 2] Finding nearby control probes...\n")
    ctrl_list <- lapply(seq_len(nrow(todo)), function(i) {
      row <- todo[i]
      nearby <- all_probes[
        probe_chr   == row$motif_chr &
        probe_start >= row$win_start &
        probe_end   <= row$win_end   &
        probe       != row$probe     &
        !(probe_start < row$motif_end & probe_end > row$motif_start)
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
  print(controls[, .N, by = .(gene, cancer, target_probe)][order(-N)][1:min(10, .N)])

  # ── STEP 3 — Load beta values ───────────────────────────────────────────────
  beta_long <- load_or_run("step3_beta_long", {
    shared         <- fread(shared_file, header = FALSE)$V1
    cancers_needed <- unique(todo$cancer)
    control_probes <- unique(controls$probe)
    target_probes  <- unique(todo$probe)

    cat("[STEP 3a] Loading target probe betas...\n")
    target_beta <- rbindlist(lapply(cancers_needed, function(cx) {
      bf <- list.files(beta_dir, pattern = paste0(cx, ".*\\.tsv\\.gz$"), full.names = TRUE)
      if (length(bf) == 0) { cat(sprintf("  [WARN] No beta matrix for %s\n", cx)); return(NULL) }
      mat <- fread(bf[1])
      if (!"probe" %in% colnames(mat)) setnames(mat, 1, "probe")
      mat <- mat[probe %in% target_probes]
      if (nrow(mat) == 0) return(NULL)
      long <- melt(mat, id.vars = "probe", variable.name = "col_name", value.name = "meth_beta")
      long[, col_name     := as.character(col_name)]
      long[, full_patient := sub("_(tumor|healthy)$", "", col_name)]
      long[, sample_type  := sub(".*_(tumor|healthy)$", "\\1", col_name)]
      long[, sample_key   := paste(full_patient, sample_type, sep = "_")]
      long <- long[sample_type %in% c("tumor", "healthy") & full_patient %in% shared]
      if (nrow(long) == 0) return(NULL)
      long[, cancer  := cx]
      long[, meth_pct := if (max(meth_beta, na.rm = TRUE) <= 1.5) meth_beta * 100 else meth_beta]
      long[, .(cancer, probe, sample_key, sample_type, meth_pct)]
    }), fill = TRUE)
    cat(sprintf("  Target beta rows: %d\n", nrow(target_beta)))

    cat("[STEP 3b] Loading control probe betas...\n")
    raw_files <- list.files(meth_dir,
      pattern = "(01A_1|11A_1)_annotated_methylation_filtered\\.bed\\.gz$",
      full.names = TRUE)
    raw_meta <- data.table(file = raw_files, fname = basename(raw_files))
    raw_meta[, cancer       := sub("^HM450_TCGA-([A-Z]+)-.*",                          "\\1", fname)]
    raw_meta[, full_patient := sub("^HM450_TCGA-([A-Z]+-TCGA-[A-Z0-9]+-[A-Z0-9]+)-.*","\\1", fname)]
    raw_meta[, sample_type  := fifelse(grepl("-01A_1_", fname), "tumor",
                               fifelse(grepl("-11A_1_", fname), "healthy", NA_character_))]
    raw_meta[, sample_key   := paste(full_patient, sample_type, sep = "_")]
    raw_meta <- raw_meta[cancer %in% cancers_needed & !is.na(sample_type) & full_patient %in% shared]
    cat(sprintf("  %d raw files to read (10 workers)...\n", nrow(raw_meta)))

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
      dt[, `:=`(cancer      = raw_meta$cancer[i],
                sample_key  = raw_meta$sample_key[i],
                sample_type = raw_meta$sample_type[i])]
      dt[, meth_pct := if (max(meth_beta, na.rm = TRUE) <= 1.5) meth_beta * 100 else meth_beta]
      dt[, .(cancer, probe, sample_key, sample_type, meth_pct)]
    }

    n_files    <- nrow(raw_meta)
    chunk_size <- 10L
    chunks     <- split(seq_len(n_files), ceiling(seq_len(n_files) / chunk_size))
    cl <- makeCluster(10L)
    clusterExport(cl, c("raw_meta", "control_probes"), envir = environment())
    clusterEvalQ(cl, library(data.table))
    res_list <- vector("list", length(chunks))
    for (k in seq_along(chunks)) {
      res_list[[k]] <- parLapply(cl, chunks[[k]], read_one)
      cat(sprintf("\r  Progress: %d / %d files (%.0f%%)",
                  min(k * chunk_size, n_files), n_files, 100 * k / length(chunks)))
      flush.console()
    }
    stopCluster(cl)
    cat("\n")

    control_beta <- rbindlist(unlist(res_list, recursive = FALSE), fill = TRUE)
    cat(sprintf("  Control beta rows: %d\n", nrow(control_beta)))
    rbindlist(list(target_beta, control_beta), fill = TRUE)
  })
  cat(sprintf("[STEP 3] Total beta rows: %d\n", nrow(beta_long)))

  # ── STEP 4 — Delta beta per probe ──────────────────────────────────────────
  delta_beta_results <- load_or_run("step4_delta_beta", {
    cat("[STEP 4] Computing delta beta per probe...\n")
    all_probe_rows <- rbind(
      todo[,     .(gene, cancer, probe, probe_type = "target")],
      controls[, .(gene, cancer, probe, probe_type = "control")]
    )
    all_probe_rows <- merge(all_probe_rows, beta_long,
                            by = c("cancer", "probe"), allow.cartesian = TRUE)
    all_probe_rows <- all_probe_rows[sample_type %in% c("tumor", "healthy") & is.finite(meth_pct)]
    all_probe_rows[,
      .(mean_beta_tumor   = mean(meth_pct[sample_type == "tumor"],   na.rm = TRUE),
        mean_beta_healthy = mean(meth_pct[sample_type == "healthy"], na.rm = TRUE),
        n_tumor           = sum(sample_type == "tumor"),
        n_healthy         = sum(sample_type == "healthy")),
      by = .(gene, cancer, probe, probe_type)
    ][, delta_beta := mean_beta_tumor - mean_beta_healthy]
  })
  cat(sprintf("[STEP 4] Delta beta — target: %d | control: %d\n",
      delta_beta_results[probe_type == "target",  .N],
      delta_beta_results[probe_type == "control", .N]))
  print(delta_beta_results[probe_type == "target"][1:min(10, .N),
        .(gene, cancer, probe, mean_beta_tumor, mean_beta_healthy, delta_beta)])

  # ── STEP 5 — Wilcoxon test ─────────────────────────────────────────────────
  wilcox_results <- load_or_run("step5_wilcox", {
    cat("[STEP 5] Wilcoxon: target delta_beta vs control delta_betas...\n")
    ctrl_map  <- unique(controls[, .(gene, cancer, probe, target_probe)])
    ctrl_db   <- merge(
      delta_beta_results[probe_type == "control" & is.finite(delta_beta),
                         .(gene, cancer, probe, ctrl_delta_beta = delta_beta)],
      ctrl_map, by = c("gene", "cancer", "probe"))
    target_db <- delta_beta_results[probe_type == "target" & is.finite(delta_beta),
                                     .(gene, cancer, probe, target_delta_beta = delta_beta)]
    res <- merge(target_db, ctrl_db,
                 by.x = c("gene", "cancer", "probe"),
                 by.y = c("gene", "cancer", "target_probe"),
                 allow.cartesian = TRUE)
    res <- res[, {
      tdb  <- target_delta_beta[1]
      cdbs <- ctrl_delta_beta
      wt <- tryCatch(wilcox.test(cdbs, mu = tdb, exact = FALSE),
               error = function(e) NULL)
      .(target_delta_beta    = tdb,
        mean_ctrl_delta_beta = mean(cdbs, na.rm = TRUE),
        n_ctrl_probes        = length(cdbs),
        delta_diff           = tdb - mean(cdbs, na.rm = TRUE),
        p_wilcox             = if (!is.null(wt)) wt$p.value else NA_real_)
    }, by = .(gene, cancer, probe)]
    res[, p_adj := p.adjust(p_wilcox, method = "BH")]
    res[order(p_wilcox)]
  })

  cat("\n=== Wilcoxon results ===\n")
  print(wilcox_results[, .(gene, cancer, probe, target_delta_beta,
                            mean_ctrl_delta_beta, n_ctrl_probes, delta_diff, p_wilcox, p_adj)])
  cat(sprintf("\nSignificant (p_wilcox < 0.05): %d / %d\n",
      wilcox_results[p_wilcox < 0.05, .N], nrow(wilcox_results)))
  cat(sprintf("Significant after BH (p_adj < 0.05): %d / %d\n",
      wilcox_results[p_adj   < 0.05, .N], nrow(wilcox_results)))

  wilcox_file <- file.path(out_dir, paste0("wilcox_delta_beta_", tf_name, "_thresh", thresh, ".tsv"))
  fwrite(wilcox_results, wilcox_file, sep = "\t")
  cat(sprintf("[DONE] Wilcoxon results saved to %s\n", wilcox_file))

  # ── STEP 6 — Scatter plot ──────────────────────────────────────────────────
  cat("[STEP 6] Generating scatter plot...\n")
  lim <- range(c(wilcox_results$target_delta_beta,
                 wilcox_results$mean_ctrl_delta_beta), na.rm = TRUE)
  pad <- diff(lim) * 0.07
  lim <- c(lim[1] - pad, lim[2] + pad)

  n_below   <- wilcox_results[delta_diff <  0,      .N]
  n_sig_neg <- wilcox_results[delta_diff < -thresh, .N]
  n_sig_pos <- wilcox_results[delta_diff >  thresh, .N]

  top10_idx <- order(abs(wilcox_results$delta_diff), decreasing = TRUE)[1:min(10L, nrow(wilcox_results))]
  wilcox_results[, label := fifelse(.I %in% top10_idx,
                                     paste(gene, cancer, sep = "\n"), NA_character_)]
  n_sig_genes <- uniqueN(wilcox_results[p_adj < 0.05 & abs(delta_diff) >= thresh, gene])

  p <- ggplot(wilcox_results, aes(x = mean_ctrl_delta_beta, y = target_delta_beta, label = label)) +
    geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.4) +
    geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.4) +
    geom_abline(slope = 1, intercept =  thresh, linetype = "dotted", colour = "#4575b4", linewidth = 0.7) +
    geom_abline(slope = 1, intercept = -thresh, linetype = "dotted", colour = "#d73027", linewidth = 0.7) +
    geom_abline(slope = 1, intercept = 0,       linetype = "dashed",  colour = "grey30",  linewidth = 0.6) +
    geom_point(aes(colour = -log10(p_adj)), size = 1.8, alpha = 0.85) +
    geom_text_repel(na.rm = TRUE, size = 2.5, max.overlaps = 20,
                    box.padding = 0.35, point.padding = 0.2,
                    segment.colour = "grey40", segment.size = 0.3) +
    scale_colour_gradient(low = "grey85", high = "#d73027", name = "-log10(p_adj)") +
    coord_equal(xlim = lim, ylim = lim) +
    annotate("text", x = lim[1] + 1, y = lim[2] - 1,
             label = sprintf(">%d%% above: %d", thresh, n_sig_pos),
             hjust = 0, vjust = 1, size = 3.5, colour = "#4575b4") +
    annotate("text", x = lim[2] - 1, y = lim[1] + 1,
             label = sprintf(">%d%% below: %d", thresh, n_sig_neg),
             hjust = 1, vjust = 0, size = 3.5, colour = "#d73027") +
    labs(
      title    = sprintf("%s — target probe delta beta vs mean control probe delta beta", tf_name),
      subtitle = sprintf("%d / %d probes below diagonal  \n ±%d%% threshold: %d negative, %d positive  \n  %d unique genes (p_adj < 0.05)",
                         n_below, nrow(wilcox_results), thresh, n_sig_neg, n_sig_pos, n_sig_genes),
      x = "Mean control probe delta beta (%)",
      y = "Target probe delta beta (%)") +
    theme_bw(base_size = 13) +
    theme(legend.position = "right", plot.subtitle = element_text(size = 10))

  plot_file <- file.path(out_dir, paste0("plot_target_vs_ctrl_delta_beta_", tf_name, "_thresh", thresh, ".png"))
  ggsave(plot_file, p, width = 6, height = 6, dpi = 150)
  cat(sprintf("[DONE] Plot saved to %s\n", plot_file))

  # ── STEP 7 — Filter significant pairs and save ─────────────────────────────
  cat("[STEP 7] Filtering pairs by p_adj < 0.05 and |delta_diff| >= thresh...\n")
  sig_pairs <- wilcox_results[p_adj < 0.05 & abs(delta_diff) >= thresh,
                               .(gene, cancer, probe, target_delta_beta,
                                 mean_ctrl_delta_beta, n_ctrl_probes,
                                 delta_diff, p_wilcox, p_adj)]
  cat(sprintf("  Significant pairs: %d / %d\n", nrow(sig_pairs), nrow(wilcox_results)))

  pairs_all <- fread(pairs_file)
  pairs_sig <- merge(pairs_all[TF == tf_name], sig_pairs, by = c("gene", "cancer", "probe"))

  filtered_file <- file.path(out_dir, paste0("filtered_pairs_delta_beta_", tf_name, "_thresh", thresh, ".tsv"))
  fwrite(pairs_sig, filtered_file, sep = "\t")
  cat(sprintf("[DONE] Filtered pairs saved to %s\n", filtered_file))

  cat("\n=== Top hits by delta_diff ===\n")
  print(pairs_sig[order(-abs(delta_diff))][1:min(20, .N),
        .(gene, cancer, probe, target_delta_beta, mean_ctrl_delta_beta, delta_diff, p_adj)])

  cat(sprintf("\n========== DONE: %s ==========\n", tf_name))
}
EOF


#############################################
# DepMap :
#############################################
mkdir -p ./results/integration/depmap

# step 1 : get the list of candidate genes from the filtered pairs
cut -f3,4 ./results/integration/meth_expr_consistent_pairs.tsv | tail -n +2 | sort -u > ./results/integration/depmap/candidate_genes.txt
wc -l ./results/integration/depmap/candidate_genes.txt

# step 2 : download the dependency and expression data from depmap portal
# Download CRISPR dependency scores from this link https://depmap.org/portal/data_page/?tab=allData&releasename=DepMap%20Public%2026Q3&filename=CRISPRGeneDependency.csv
# Download cell line metadata (needed to match cancer types) from this link https://depmap.org/portal/data_page/?tab=allData&releasename=DepMap%20Public%2026Q3&filename=Model.csv 
# get the expression data from this link : https://depmap.org/portal/data_page/?tab=allData&releasename=DepMap%20Public%2026Q3&filename=OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv 


# Step 3: Integrate DepMap expression and dependency
# using cancer-specific, matched cell lines
Rscript - << 'EOF'
base <- "/data/najemd/TF_binding_cancer/results/integration/depmap"

candidates <- unique(read.table(
  file.path(base, "candidate_genes.txt"),
  header = FALSE, sep = "\t",
  col.names = c("gene", "cancer"),
  stringsAsFactors = FALSE,
  quote = "", comment.char = ""
))

candidates <- candidates[
  !candidates$gene %in% c("PCDHGA5", "PCDHGA6", "SPON2"), ,
  drop = FALSE
]

model <- read.csv(
  file.path(base, "Model.csv"),
  stringsAsFactors = FALSE
)

tcga_codes <- list(
  BLCA = "BLCA",
  BRCA = c("BRCANOS", "BRCNOS", "IDC", "ILC", "MDLC", "BRCA"),
  COAD = "COAD",
  HNSC = c("HNSC", "HPHSC", "LXSC", "OCSC", "OPHSC"),
  KIRC = "CCRCC",
  KIRP = "PRCC",
  LIHC = "HCC",
  LUAD = "LUAD",
  PRAD = "PRAD",
  THCA = "THPA",
  UCEC = c("UCEC", "UEC", "USC", "UMEC")
)

stopifnot(all(candidates$cancer %in% names(tcga_codes)))

# Load only metadata and candidate-gene columns
read_candidates <- function(filename, expression = FALSE) {
  path <- file.path(base, filename)
  headers <- names(read.csv(path, nrows = 0, check.names = FALSE))

  symbols <- sub(" \\([0-9]+\\)$", "", headers)
  gene_cols <- which(
    grepl(" \\([0-9]+\\)$", headers) &
    symbols %in% candidates$gene
  )

  metadata_cols <- if (expression) {
    match(c("ModelID", "IsDefaultEntryForModel"), headers)
  } else {
    1L
  }

  stopifnot(!anyNA(metadata_cols))

  classes <- rep("NULL", length(headers))
  classes[metadata_cols] <- "character"
  classes[gene_cols] <- "numeric"

  dat <- read.csv(
    path, colClasses = classes,
    check.names = FALSE, stringsAsFactors = FALSE
  )

  if (!expression) names(dat)[1] <- "ModelID"
  dat
}

cat("Loading dependency data...\n")
dep <- read_candidates("CRISPRGeneDependency.csv")

cat("Loading expression data...\n")
expr <- read_candidates(
  "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv",
  expression = TRUE
)

expr <- expr[
  !is.na(expr$IsDefaultEntryForModel) &
  expr$IsDefaultEntryForModel == "Yes", ,
  drop = FALSE
]

valid_ids <- function(x) !is.na(x) & nzchar(x)
dep <- dep[valid_ids(dep$ModelID), , drop = FALSE]
expr <- expr[valid_ids(expr$ModelID), , drop = FALSE]

stopifnot(
  !anyDuplicated(dep$ModelID),
  !anyDuplicated(expr$ModelID),
  !anyDuplicated(model$ModelID)
)

# Match symbols exactly, rather than taking the first regex hit
gene_column <- function(dat, gene) {
  headers <- names(dat)
  symbols <- sub(" \\([0-9]+\\)$", "", headers)
  hits <- which(
    grepl(" \\([0-9]+\\)$", headers) & symbols == gene
  )
  if (length(hits) != 1) {
    stop("Expected exactly one column for gene: ", gene)
  }
  headers[hits]
}

shared_ids <- intersect(dep$ModelID, expr$ModelID)
summary_rows <- list()
detail_rows <- list()

cat("Calculating matched summaries...\n")

for (i in seq_len(nrow(candidates))) {
  g <- candidates$gene[i]
  cancer <- candidates$cancer[i]

  ids <- model$ModelID[
    model$OncotreeCode %in% tcga_codes[[cancer]]
  ]
  ids <- intersect(ids, shared_ids)

  dcol <- gene_column(dep, g)
  ecol <- gene_column(expr, g)

  # Confirm that both columns refer to the same numeric gene ID
  get_id <- function(x) sub(".* \\(([0-9]+)\\)$", "\\1", x)
  stopifnot(get_id(dcol) == get_id(ecol))

  d <- dep[[dcol]][match(ids, dep$ModelID)]
  e <- expr[[ecol]][match(ids, expr$ModelID)]

  # Keep identical models with valid values for both assays
  keep <- is.finite(d) & is.finite(e)
  ids_valid <- ids[keep]
  d <- d[keep]
  e <- e[keep]

  if (any(d < 0 | d > 1)) {
    stop("Dependency probabilities outside 0–1 for ", g)
  }

  n <- length(d)

  summary_rows[[i]] <- data.frame(
    gene = g,
    cancer = cancer,
    n_shared_models = length(ids),
    n_models_used = n,
    n_missing_either = length(ids) - n,
    median_expr = if (n > 0) median(e) else NA_real_,
    median_dependency = if (n > 0) median(d) else NA_real_,
    max_dependency = if (n > 0) max(d) else NA_real_,
    n_dependent = if (n > 0) sum(d > 0.5) else NA_integer_,
    fraction_dependent = if (n > 0) mean(d > 0.5) else NA_real_,
    status = if (n == 0) "no_paired_values" else
             if (n == 1) "single_model" else "multiple_models"
  )

  if (n > 0) {
    idx <- match(ids_valid, model$ModelID)
    cell_names <- if ("CellLineName" %in% names(model)) {
      model$CellLineName[idx]
    } else {
      rep(NA_character_, n)
    }

    detail_rows[[length(detail_rows) + 1L]] <- data.frame(
      gene = g,
      cancer = cancer,
      ModelID = ids_valid,
      CellLineName = cell_names,
      OncotreeCode = model$OncotreeCode[idx],
      expression_log2TPM1 = e,
      dependency_probability = d,
      dependent = d > 0.5
    )
  }
}

results <- do.call(rbind, summary_rows)
results <- results[
  order(-results$median_dependency, na.last = TRUE), ,
  drop = FALSE
]

write.table(
  results,
  file.path(base, "depmap_matched_summary.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE, na = "NA"
)

if (length(detail_rows) > 0) {
  details <- do.call(rbind, detail_rows)
  write.table(
    details,
    file.path(base, "depmap_matched_cell_lines.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE, na = "NA"
  )
}

cat("\nSummary rows:", nrow(results), "\n")
cat("\nStatus counts:\n")
print(table(results$status))

cat("\nTop 15 pairs by median dependency probability:\n")
print(head(results, 15), row.names = FALSE)

cat("\nSaved corrected results in:", base, "\n")
EOF


# Step 4: Plot matched expression vs dependency
Rscript - << 'EOF'
library(ggplot2)

base <- "/data/najemd/TF_binding_cancer/results/integration/depmap"

results <- read.table(
  file.path(base, "depmap_matched_summary.tsv"),
  header = TRUE, sep = "\t",
  stringsAsFactors = FALSE
)

# Plot only pairs with valid values for both measurements
plot_data <- results[
  results$n_models_used > 0 &
  is.finite(results$median_expr) &
  is.finite(results$median_dependency), ,
  drop = FALSE
]

if (nrow(plot_data) == 0) {
  stop("No valid gene-cancer pairs available for plotting.")
}

cat("Pairs plotted:", nrow(plot_data), "\n")
cat("Pairs omitted:", nrow(results) - nrow(plot_data), "\n")

# Load your existing cancer colours
colors_df <- read.table(
  "/data/najemd/TF_binding_cancer/results/multi_omics/cancer_color_order_with_defined_colours.tsv",
  header = FALSE, sep = "\t",
  col.names = c("cancer_full", "color"),
  stringsAsFactors = FALSE,
  comment.char = "", quote = ""
)

colors_df$cancer <- sub("^TCGA-", "", colors_df$cancer_full)
color_map <- setNames(colors_df$color, colors_df$cancer)

missing_colors <- setdiff(unique(plot_data$cancer), names(color_map))
if (length(missing_colors) > 0) {
  stop("Missing colours for: ", paste(missing_colors, collapse = ", "))
}

# Distinguish single-model results
plot_data$model_group <- factor(
  ifelse(
    plot_data$n_models_used == 1,
    "Single model", "Multiple models"
  ),
  levels = c("Multiple models", "Single model")
)

# Label the top 10 pairs; this is a ranking, not an essentiality call
plot_data <- plot_data[
  order(-plot_data$median_dependency, plot_data$gene, plot_data$cancer), ,
  drop = FALSE
]

labels <- head(plot_data, 10)
labels$label <- paste0(labels$gene, " (", labels$cancer, ")")

p <- ggplot(
  plot_data,
  aes(
    x = median_expr,
    y = median_dependency,
    color = cancer,
    shape = model_group
  )
) +
  geom_point(size = 2.8, alpha = 0.85) +
  scale_color_manual(values = color_map, name = "Cancer type") +
  scale_shape_manual(
    values = c("Multiple models" = 16, "Single model" = 17),
    name = "Model availability"
  ) +
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0.03, 0.15))
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.05, 0.12))
  ) +
  labs(
    title = "Expression vs dependency in DepMap cell lines",
    subtitle = "Candidate genes from TCGA methylation–expression analysis",
    x = "Median expression (log2(TPM + 1))",
    y = "Median probability of dependency",
    caption = paste(
      "Each point represents one gene–cancer pair;",
      "both medians use the same models with valid paired measurements.",
      "\nTriangles represent a single model.",
      "Labels identify the 10 highest medians, not confirmed essential genes."
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.caption = element_text(hjust = 0, size = 9),
    legend.position = "right"
  )

# Use non-overlapping labels if ggrepel is available
if (requireNamespace("ggrepel", quietly = TRUE)) {
  p <- p + ggrepel::geom_text_repel(
    data = labels,
    aes(label = label),
    size = 3,
    seed = 42,
    max.overlaps = Inf,
    box.padding = 0.5,
    show.legend = FALSE
  )
} else {
  p <- p + geom_text(
    data = labels,
    aes(label = label),
    size = 3,
    nudge_y = 0.008,
    check_overlap = TRUE,
    show.legend = FALSE
  )
  cat("ggrepel is unavailable: some overlapping labels may be omitted.\n")
}

ggsave(
  file.path(base, "depmap_matched_scatter.pdf"),
  plot = p, width = 11, height = 7.5
)

ggsave(
  file.path(base, "depmap_matched_scatter.png"),
  plot = p, width = 11, height = 7.5,
  dpi = 300, bg = "white"
)

cat("Saved depmap_matched_scatter.pdf and depmap_matched_scatter.png\n")
EOF




#######################################
# Adding Mutation information
#######################################
##################################
# Sample inventory — run before anything else
##################################

##################################
# Create lists of snv files, and a list of shared samples methylation with mutation
##################################

SHARED_LIST="./results/multi_omics/expr_meth_shared_samples.txt"
SNV_LIST="./results/multi_omics/SNV_files.txt"
SHARED_SNV_LIST="./results/multi_omics/expr_meth_snv_shared_samples.txt"

# 1. List SNV files (vcf.gz only, not .tbi)
echo "=== SNV files ==="
find ./snv/ -name "*.vcf.gz" ! -name "*.tbi" | sort > "$SNV_LIST"
echo "  $(wc -l < $SNV_LIST) files found"

# 2. Extract sample IDs to match format: ACC-TCGA-OR-A5J1
sed 's|.*/SNV_TCGA-||' "$SNV_LIST" \
  | grep -oP '^[A-Z]+-TCGA-[A-Z0-9]+-[A-Z0-9]+-[A-Z0-9]+' \
  | sed 's/-[0-9][0-9][A-Z]$//' \
  | sort -u > /tmp/snv_ids.txt
echo "  $(wc -l < /tmp/snv_ids.txt) unique patient IDs"
echo "  Sample: $(head -3 /tmp/snv_ids.txt | tr '\n' ' ')"

# 3. Intersect with existing shared samples
sort "$SHARED_LIST" > /tmp/shared_sorted.txt
comm -12 /tmp/shared_sorted.txt /tmp/snv_ids.txt > "$SHARED_SNV_LIST"

echo ""
echo "=== Summary ==="
echo "Expr+Meth shared : $(wc -l < $SHARED_LIST)"
echo "SNV files        : $(wc -l < $SNV_LIST)"
echo "All three shared : $(wc -l < $SHARED_SNV_LIST)"

# barplot the 3d 
Rscript - <<'EOF'
library(ggplot2)
library(data.table)

# ── Expr + Meth samples (from beta matrices) ──────────────────────────────
beta_dir <- "./methylation/beta_matrices_BANP"
files <- list.files(beta_dir, pattern = "\\.tsv\\.gz$", full.names = TRUE)

dt_meth <- rbindlist(lapply(files, function(f) {
  cols <- colnames(fread(f, nrows = 0))
  cols <- cols[cols != "probe"]
  cancer <- sub("_BANP_beta_matrix\\.tsv\\.gz$", "", basename(f))
  data.table(sample = cols, cancer = cancer)
}))
dt_meth[, condition := ifelse(grepl("_tumor$", sample), "Tumor (expr+meth)", "Healthy (expr+meth)")]

# ── Expr + Meth + SNV samples ─────────────────────────────────────────────
snv_shared <- fread("./results/multi_omics/expr_meth_snv_shared_samples.txt",
                    header = FALSE)$V1
# Format: ACC-TCGA-OR-A5J1 — extract cancer type
dt_snv <- data.table(sample = snv_shared)
dt_snv[, cancer    := sub("-TCGA-.*", "", sample)]
dt_snv[, condition := "SNV (expr+meth+snv)"]

# ── Combine ───────────────────────────────────────────────────────────────
dt <- rbind(dt_meth[, .(cancer, condition)],
            dt_snv[,  .(cancer, condition)])

counts <- dt[, .N, by = .(cancer, condition)]
totals <- counts[, .(total = sum(N)), by = cancer]
counts[, cancer := factor(cancer, levels = totals[order(-total), cancer])]
counts[, condition := factor(condition,
                             levels = c("Tumor (expr+meth)",
                                        "Healthy (expr+meth)",
                                        "SNV (expr+meth+snv)"))]

ggplot(counts, aes(x = cancer, y = N, fill = condition)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.9, preserve = "single")) +
  geom_text(aes(label = N),
            position = position_dodge(width = 0.9, preserve = "single"),
            hjust = 0.5, vjust = -0.5, size = 2.5) +
  scale_fill_manual(values = c(
    "Tumor (expr+meth)"     = "#C0392B",
    "Healthy (expr+meth)"   = "#2471A3",
    "SNV (expr+meth+snv)"   = "#27AE60")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title    = "Samples per cancer type across data types",
       subtitle = paste0("Expr+Meth: ", nrow(dt_meth),
                         "  |  Expr+Meth+SNV: ", nrow(dt_snv)),
       x = NULL, y = "Samples", fill = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA))

ggsave("./results/multi_omics/samples_3d_per_cancer.png",
       width = 12, height = 5, dpi = 150)
cat("Saved: ./results/multi_omics/samples_3d_per_cancer.png\n")
EOF


#!/bin/bash
# ══════════════════════════════════════════════════════════════════════════════
# Mutation enrichment analysis at TF motif level — Fisher's exact test
#
# Goal: Test whether motif instances containing a consistent probe
#       (significant methylation-expression dysregulation) are more frequently
#       mutated than motif instances without a consistent probe.
#       Uses Fisher's exact test on the proportion of motifs carrying
#       at least one PASS somatic mutation across all TCGA samples,
#       both pooled across all cancer types and per cancer type separately.
#
# Steps:
#   Step 1 — Build motif-level BED from all motif instances
#   Step 2 — Intersect PASS somatic mutations with motif windows
#   Step 3 — Fisher's exact test (consistent vs pure background, per TF,
#             pooled across all cancers and per cancer type)
#   Step 4 — Plot results (one PDF: pooled page + one page per cancer type)
#
# Input:
#   ./methylation/probe_gene_pairs_in_motifs.tsv
#   ./results/integration/meth_expr_consistent_pairs.tsv
#   ./results/multi_omics/expr_meth_snv_shared_samples.txt
#   ./snv/*.vcf.gz
#
# Output:
#   ./results/integration/mutation_analysis/
#     motif_windows.bed
#     mutation_motif_hits_all.tsv
#     motif_mutation_fisher.tsv
#     motif_mutation_enrichment.pdf
# ══════════════════════════════════════════════════════════════════════════════

set -euo pipefail

OUT_DIR="./results/integration/mutation_analysis"
MOTIF_PAIRS="./methylation/probe_gene_pairs_in_motifs.tsv"
CONSISTENT="./results/integration/meth_expr_consistent_pairs.tsv"
SHARED_SNV="./results/multi_omics/expr_meth_snv_shared_samples.txt"
SNV_DIR="./snv"
THREADS=40

mkdir -p "$OUT_DIR"

# ══════════════════════════════════════════════════════════════════════════════
# STEP 1 — Build motif-level BED (one row per unique motif instance)
#
# Purpose: Extract motif coordinates from motif_instance_id column.
#          Format: chr:start:end:strand → BED (0-based).
#          Each motif is ~10 bp — used as-is, no padding needed.
#          Skipped if already done.
# ══════════════════════════════════════════════════════════════════════════════
if [ -f "$OUT_DIR/motif_windows.bed" ] && [ -s "$OUT_DIR/motif_windows.bed" ]; then
  echo "[STEP 1] motif_windows.bed already exists, skipping."
  echo "  $(wc -l < "$OUT_DIR/motif_windows.bed") motif windows"
else
  echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
  echo "[STEP 1] Building motif-level BED..."
  echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
  tail -n +2 "$MOTIF_PAIRS" \
    | awk 'BEGIN{OFS="\t"} {
        split($4, a, ":");
        print a[1], a[2], a[3], $4, $1
      }' \
    | sort -k1,1 -k2,2n \
    | sort -u -k4,4 > "$OUT_DIR/motif_windows.bed"
  echo "  $(wc -l < "$OUT_DIR/motif_windows.bed") unique motif instances"
fi

# ══════════════════════════════════════════════════════════════════════════════
# STEP 2 — Intersect PASS somatic mutations with motif windows
#
# Purpose: For each sample with SNV data, extract PASS-filtered mutations,
#          convert VCF 1-based coordinates to 0-based BED, and intersect
#          with motif windows. Output: motif_instance_id, sample_id, cancer.
#          Counter messages are sent to stderr to keep the hits file clean.
#          Skipped if already done.
# ══════════════════════════════════════════════════════════════════════════════
if [ -f "$OUT_DIR/mutation_motif_hits_all.tsv" ] && [ -s "$OUT_DIR/mutation_motif_hits_all.tsv" ]; then
  echo "[STEP 2] mutation_motif_hits_all.tsv already exists, skipping."
  echo "  $(wc -l < "$OUT_DIR/mutation_motif_hits_all.tsv") hits"
else
  echo ""
  echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
  echo "[STEP 2] Intersecting mutations with motif windows ($THREADS cores)..."
  echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

  HITS_FILE="$OUT_DIR/mutation_motif_hits_all.tsv"
  MOTIF_WINDOWS="$OUT_DIR/motif_windows.bed"
  > "$HITS_FILE"

  process_sample() {
    sample_id="$1"
    SNV_DIR="$2"
    MOTIF_WINDOWS="$3"
    cancer=$(echo "$sample_id" | cut -d'-' -f1)
    patient=$(echo "$sample_id" | cut -d'-' -f2-)
    vcf=$(find "$SNV_DIR" -maxdepth 1 -name "*${patient}-01A*.vcf.gz" ! -name "*.tbi" 2>/dev/null | head -1)
    [ -z "$vcf" ] && return
    bcftools view -f PASS "$vcf" 2>/dev/null \
      | awk '!/^#/ {print $1"\t"($2-1)"\t"$2}' \
      | bedtools intersect -a "$MOTIF_WINDOWS" -b - \
      | awk -v sid="$sample_id" -v cx="$cancer" \
          'BEGIN{OFS="\t"} {print $4, sid, cx}'
  }
  export -f process_sample

  n_samples=$(wc -l < "$SHARED_SNV")
  echo "  Processing $n_samples samples..."

  COUNTER_DIR=$(mktemp -d)
  cat "$SHARED_SNV" \
    | xargs -P "$THREADS" -I{} bash -c '
        process_sample "$@"
        touch "'"$COUNTER_DIR"'"/done_$$
        n=$(ls "'"$COUNTER_DIR"'" | wc -l)
        if [ $((n % 500)) -eq 0 ]; then
          echo "  ... $n / '"$n_samples"' samples processed" >&2
        fi
      ' _ {} "$SNV_DIR" "$MOTIF_WINDOWS" \
    >> "$HITS_FILE"

  rm -rf "$COUNTER_DIR"
  echo "  $(wc -l < "$HITS_FILE") mutation-motif hits → $HITS_FILE"
fi

# ══════════════════════════════════════════════════════════════════════════════
# STEP 3 — Fisher's exact test (R)
#
# Purpose: For each motif, determine whether it carries at least one PASS
#          somatic mutation (binary: mutated / not mutated).
#          Consistent motifs = contain a significant probe in a given cancer.
#          Pure background = motifs never consistent in any cancer type.
#          Runs two tests per TF:
#            (a) Pooled — mutations across all cancer types combined
#            (b) Per cancer — mutations restricted to that cancer type,
#                consistent motifs specific to that cancer vs pure background
#          BH correction applied across all tests.
# ══════════════════════════════════════════════════════════════════════════════
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "[STEP 3] Fisher's exact test (R)..."
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

Rscript - <<'EOF'
library(data.table)

out_dir    <- "./results/integration/mutation_analysis"
hits_file  <- file.path(out_dir, "mutation_motif_hits_all.tsv")
motif_file <- "./methylation/probe_gene_pairs_in_motifs.tsv"
cons_file  <- "./results/integration/meth_expr_consistent_pairs.tsv"

cat("[R] Loading data...\n")
hits  <- fread(hits_file, header = FALSE,
               col.names = c("motif_instance_id", "sample_id", "cancer"))
motif <- fread(motif_file)
cons  <- fread(cons_file)

motif[, TF_clean := sub("_.*", "", TF)]

# All unique motifs with TF label
all_motifs <- unique(motif[, .(motif_instance_id, TF_clean)])

# Motifs consistent in ANY cancer (used to define pure background)
cons_any <- unique(cons[, .(motif_instance_id, TF_clean = TF)])
cons_any[, ever_consistent := TRUE]
all_motifs <- merge(all_motifs, cons_any,
                    by = c("motif_instance_id", "TF_clean"), all.x = TRUE)
all_motifs[is.na(ever_consistent), ever_consistent := FALSE]

# Pure background = never consistent in any cancer
pure_background <- all_motifs[ever_consistent == FALSE]

# ── POOLED Fisher test (all cancers) ─────────────────────────────────────────
mutated_all <- unique(hits[, .(motif_instance_id)])
mutated_all[, any_mutation := TRUE]

all_motifs_pool <- merge(all_motifs, mutated_all,
                         by = "motif_instance_id", all.x = TRUE)
all_motifs_pool[is.na(any_mutation), any_mutation := FALSE]

cons_motifs_pool <- unique(cons[, .(motif_instance_id, TF_clean = TF)])
cons_motifs_pool[, group := "consistent"]
all_motifs_pool <- merge(all_motifs_pool, cons_motifs_pool,
                         by = c("motif_instance_id", "TF_clean"), all.x = TRUE)
all_motifs_pool[is.na(group), group := "background"]

pooled_results <- rbindlist(lapply(c("BANP", "NRF1"), function(tf) {
  sub       <- all_motifs_pool[TF_clean == tf]
  con_mut   <- sum(sub[group == "consistent",  any_mutation])
  con_nomut <- sum(!sub[group == "consistent", any_mutation])
  bg_mut    <- sum(sub[group == "background",  any_mutation])
  bg_nomut  <- sum(!sub[group == "background", any_mutation])
  mat <- matrix(c(con_mut, con_nomut, bg_mut, bg_nomut), nrow = 2)
  ft  <- fisher.test(mat, alternative = "greater")
  cat(sprintf("\n=== POOLED %s === Fisher p = %.4f | OR = %.3f\n",
              tf, ft$p.value, ft$estimate))
  data.table(
    cancer      = "ALL",
    TF          = tf,
    con_mutated = con_mut,
    con_total   = con_mut + con_nomut,
    bg_mutated  = bg_mut,
    bg_total    = bg_mut + bg_nomut,
    pct_mut_con = round(100 * con_mut / (con_mut + con_nomut), 2),
    pct_mut_bg  = round(100 * bg_mut  / (bg_mut  + bg_nomut),  2),
    odds_ratio  = ft$estimate,
    fisher_p    = ft$p.value
  )
}))

# ── PER-CANCER Fisher test ────────────────────────────────────────────────────
# Consistent motifs are cancer-specific; background is pure (never consistent
# in any cancer). Mutations are restricted to the cancer type being tested.
cancers <- unique(cons$cancer)

per_cancer_results <- rbindlist(lapply(cancers, function(cx) {
  rbindlist(lapply(c("BANP", "NRF1"), function(tf) {

    cons_cx   <- unique(cons[cancer == cx & TF == tf, .(motif_instance_id)])
    if (nrow(cons_cx) == 0) return(NULL)

    bg_motifs <- pure_background[TF_clean == tf, .(motif_instance_id)]
    if (nrow(bg_motifs) == 0) return(NULL)

    hits_cx    <- hits[cancer == cx]
    mutated_cx <- unique(hits_cx[, .(motif_instance_id)])
    mutated_cx[, any_mutation := TRUE]

    con_dt <- merge(cons_cx,   mutated_cx, by = "motif_instance_id", all.x = TRUE)
    bg_dt  <- merge(bg_motifs, mutated_cx, by = "motif_instance_id", all.x = TRUE)
    con_dt[is.na(any_mutation), any_mutation := FALSE]
    bg_dt[is.na(any_mutation),  any_mutation := FALSE]

    con_mut   <- sum(con_dt$any_mutation)
    con_nomut <- sum(!con_dt$any_mutation)
    bg_mut    <- sum(bg_dt$any_mutation)
    bg_nomut  <- sum(!bg_dt$any_mutation)

    if ((con_mut + con_nomut) == 0 || (bg_mut + bg_nomut) == 0) return(NULL)

    mat <- matrix(c(con_mut, con_nomut, bg_mut, bg_nomut), nrow = 2)
    ft  <- tryCatch(fisher.test(mat, alternative = "greater"), error = function(e) NULL)
    if (is.null(ft)) return(NULL)

    cat(sprintf("  %s | %s: Fisher p = %.4f | OR = %.3f | con=%d/%d | bg=%d/%d\n",
                cx, tf, ft$p.value, ft$estimate,
                con_mut, con_mut + con_nomut,
                bg_mut,  bg_mut  + bg_nomut))

    data.table(
      cancer      = cx,
      TF          = tf,
      con_mutated = con_mut,
      con_total   = con_mut + con_nomut,
      bg_mutated  = bg_mut,
      bg_total    = bg_mut + bg_nomut,
      pct_mut_con = round(100 * con_mut / (con_mut + con_nomut), 2),
      pct_mut_bg  = round(100 * bg_mut  / (bg_mut  + bg_nomut),  2),
      odds_ratio  = ft$estimate,
      fisher_p    = ft$p.value
    )
  }), fill = TRUE)
}), fill = TRUE)

# Combine, apply BH correction, save
all_results <- rbind(pooled_results, per_cancer_results, fill = TRUE)
all_results[, p_adj := p.adjust(fisher_p, method = "BH")]
fwrite(all_results, file.path(out_dir, "motif_mutation_fisher.tsv"), sep = "\t")
cat("\n[R] Results saved.\n")
EOF

# ══════════════════════════════════════════════════════════════════════════════
# STEP 4 — Plot (R)
#
# Purpose: Bar charts showing percentage of motifs mutated in consistent vs
#          pure background groups for each TF, y-axis fixed 0–100%.
#          Output: one PDF with:
#            - Page 1: pooled results (BANP + NRF1 side by side)
#            - One page per cancer type (BANP + NRF1 side by side),
#              annotated with Fisher p, BH-adjusted p, and odds ratio.
# ══════════════════════════════════════════════════════════════════════════════
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "[STEP 4] Plotting (R)..."
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

Rscript - <<'EOF'
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

out_dir <- "./results/integration/mutation_analysis"
res     <- fread(file.path(out_dir, "motif_mutation_fisher.tsv"))

pdf(file.path(out_dir, "motif_mutation_enrichment.pdf"), width = 8, height = 5)

# ── Page 1: pooled (all cancers) ─────────────────────────────────────────────
plots_pooled <- lapply(c("BANP", "NRF1"), function(tf) {
  r <- res[TF == tf & cancer == "ALL"]
  if (nrow(r) == 0) return(NULL)
  lab <- sprintf("Fisher p = %.4f\nOR = %.2f", r$fisher_p, r$odds_ratio)
  dat <- data.table(
    group   = factor(c("consistent", "background"),
                     levels = c("background", "consistent")),
    pct_mut = c(r$pct_mut_con, r$pct_mut_bg)
  )
  ggplot(dat, aes(x = group, y = pct_mut, fill = group)) +
    geom_col(width = 0.5) +
    scale_fill_manual(values = c(background = "#AAAAAA", consistent = "#E05A2B")) +
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
    annotate("text", x = 1.5, y = 85, label = lab, size = 3.5) +
    labs(title = sprintf("%s — All cancers pooled", tf),
         x = NULL, y = "% motifs with ≥1 PASS mutation") +
    theme_classic(base_size = 12) +
    theme(legend.position = "none")
})
print(wrap_plots(plots_pooled, ncol = 2) +
        plot_annotation(title = "Pooled — all cancer types"))

# ── One page per cancer type ──────────────────────────────────────────────────
cancers <- sort(unique(res[cancer != "ALL", cancer]))

for (cx in cancers) {
  plots_cx <- lapply(c("BANP", "NRF1"), function(tf) {
    r <- res[TF == tf & cancer == cx]
    if (nrow(r) == 0) return(ggplot() + theme_void() +
                               labs(title = sprintf("%s — %s: no data", tf, cx)))
    lab <- sprintf("Fisher p = %.4f\np_adj = %.4f\nOR = %.2f",
                   r$fisher_p, r$p_adj, r$odds_ratio)
    dat <- data.table(
      group   = factor(c("consistent", "background"),
                       levels = c("background", "consistent")),
      pct_mut = c(r$pct_mut_con, r$pct_mut_bg)
    )
    ggplot(dat, aes(x = group, y = pct_mut, fill = group)) +
      geom_col(width = 0.5) +
      scale_fill_manual(values = c(background = "#AAAAAA", consistent = "#E05A2B")) +
      scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
      annotate("text", x = 1.5, y = 85, label = lab, size = 3.2) +
      labs(title = sprintf("%s — %s", tf, cx),
           x = NULL, y = "% motifs with ≥1 PASS mutation") +
      theme_classic(base_size = 12) +
      theme(legend.position = "none")
  })
  print(wrap_plots(plots_cx, ncol = 2) +
          plot_annotation(title = sprintf("Cancer: %s", cx)))
}

dev.off()
cat("[DONE] PDF saved.\n")
EOF

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "[DONE] All outputs in $OUT_DIR"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"


###############################################################################
# SNV IN KNOWN GENE NFX1 
###############################################################################

###############################################################################
# GETTING NFX1 GENOMIC COORDINATES
###############################################################################
# genes annotations found in  : /data/genome/annotations/hg38_genes.bed
mkdir ./snv/NFX1_gene/
cat /data/genome/annotations/hg38_genes.bed | grep 'ENSG00000086102' > ./snv/NFX1_gene/NFX1_genomic_coordinates.txt
# chr9    33290510        33371157        NFX1    +       ENSG00000086102

###############################################################################
# FINDING SNV IN BANP MOTIF THAT OVERLAP WITH THE NFX1 GENE 
###############################################################################
# looking into the file ./snv/overlaps/intersected_motifs_snv_filtered/BANP_intersected_SNVs.bed
cat ./snv/overlaps/intersected_motifs_snv_filtered/BANP_intersected_SNVs.bed | grep 'chr9' | awk '$2 > 33290510 && $3 < 33371157 {print $4}'
# chr9:33329363:C:T
# chr9:33329368:G:A
# chr9:33339406:T:G
# chr9:33339410:A:C
# chr9:33339411:C:T 

cat ./motifs/BANP_filtered_6mer.MA2503.1.bed | grep 'chr9' | grep 'TATCGCGAGA'
# chr9    20670569        20670579        TATCGCGAGA      5.01e-06        +
# chr9    33290683        33290693        TATCGCGAGA      5.01e-06        + <-- THIS SHOULD BE WHERE I HAVE THE SNV
# chr9    38354804        38354814        TATCGCGAGA      5.01e-06        +
# chr9    127451451       127451461       TATCGCGAGA      5.01e-06        +

cat ./snv/overlaps/intersected_motifs_snv_filtered/BANP_intersected_SNVs.bed |awk '$1=="chr9" && $2 <= 33290690 && $3 >= 33290690 {print $0}'
# pos_dbSNP = 33290690
cat ./snv/all_unique_variants_across_cancers.bed | grep 'chr9' | awk ' $2 >=33290688 && $3 <=33290695 {print $0}'
# chr9    33290694        33290695        chr9:33290695:T:C

# LOOKED INTO THE PAPER WHERE IT WAS FOUND , THEY MAPPED ON THE HG19 GENOME WHILE MINE ARE ON HG38. PLUS THE FOUND SNP IS A GERMLINE VARIATION WHILE MUTECT2 IS FOR SOMATIC MUTATIONS.

find ./snv/ -name "SNV_TCGA*.vcf.gz" | \
  parallel -j 10 --bar 'zgrep -P "^chr9\t33290690\t" {} 2>/dev/null | sed "s|^|{}\t|"' \
  > ./results/integration/mutation_analysis/results_rs4879682.txt




##############################################################
# Overlapping my results with cg island and gene body position
##############################################################
zcat /data/genome/annotations/hg38_cpgIsland.bed.gz

# Extract motif coordinates + TF name
awk 'NR>1 {
    split($2, a, ":");
    print a[1] "\t" a[2] "\t" a[3] "\t" $2 "\t" $1 "\t" $24
}' /data/najemd/TF_binding_cancer/results/integration/meth_expr_consistent_pairs_cancermine_labeled.tsv \
> /data/najemd/TF_binding_cancer/results/integration/consistent_motifs_TF.bed

# Intersect with CpG islands
bedtools intersect \
  -a /data/najemd/TF_binding_cancer/results/integration/consistent_motifs_TF.bed \
  -b /data/genome/annotations/hg38_cpgIsland.bed.gz -u \
> /data/najemd/TF_binding_cancer/results/integration/motifs_in_cpgi_TF.bed

# Per-TF summary table: total | in_CpGI | percentage
echo -e "TF\ttotal\tin_CpGI\tpct_in_CpGI" \
> /data/najemd/TF_binding_cancer/results/integration/cpgi_per_TF.tsv

awk 'NR==FNR {total[$6]++; next} {cpgi[$6]++}
END {
    for (tf in total)
        printf "%s\t%d\t%d\t%.1f\n", tf, total[tf], cpgi[tf]+0, (cpgi[tf]+0)/total[tf]*100
}' \
  /data/najemd/TF_binding_cancer/results/integration/consistent_motifs_TF.bed \
  /data/najemd/TF_binding_cancer/results/integration/motifs_in_cpgi_TF.bed \
| sort -t$'\t' -k4 -rn \
>> /data/najemd/TF_binding_cancer/results/integration/cpgi_per_TF.tsv

# Preview
cat /data/najemd/TF_binding_cancer/results/integration/cpgi_per_TF.tsv 

# TF      total   in_CpGI pct_in_CpGI
# NRF1    176     124     70.5
# BANP    44      23      52.3

# add column for cpg island overlap to the main results table
# Extract just the motif_instance_ids that overlap CpG islands (column 4)
bedtools intersect \
  -a /data/najemd/TF_binding_cancer/results/integration/consistent_motifs_TF.bed \
  -b /data/genome/annotations/hg38_cpgIsland.bed.gz -u \
| awk '{print $4}' | sort -u \
> /data/najemd/TF_binding_cancer/results/integration/cpgi_motif_ids.txt

wc -l /data/najemd/TF_binding_cancer/results/integration/cpgi_motif_ids.txt
# 80 /data/najemd/TF_binding_cancer/results/integration/cpgi_motif_ids.txt

python3 - << 'EOF'
with open("/data/najemd/TF_binding_cancer/results/integration/cpgi_motif_ids.txt") as f:
    cpgi_set = set(line.strip() for line in f)

print(f"CpGI motifs loaded: {len(cpgi_set)}")

with open("/data/najemd/TF_binding_cancer/results/integration/meth_expr_consistent_pairs_cancermine_labeled.tsv") as fin, \
     open("/data/najemd/TF_binding_cancer/results/integration/meth_expr_consistent_pairs_cancermine_cpgi.tsv", "w") as fout:
    
    for i, line in enumerate(fin):
        line = line.rstrip("\n")
        if i == 0:
            fout.write(line + "\tin_CpGI\n")
            continue
        fields = line.split("\t")
        motif_id = fields[1]
        label = "yes" if motif_id in cpgi_set else "no"
        fout.write(line + "\t" + label + "\n")

print("Done!")
EOF

# Barplot the results 
Rscript -e '
library(ggplot2)
library(data.table)

cpgi <- fread("/data/najemd/TF_binding_cancer/results/integration/cpgi_per_TF.tsv")
cpgi[, not_in_CpGI := total - in_CpGI]

cpgi_long <- melt(cpgi[, .(TF, in_CpGI, not_in_CpGI)],
                  id.vars       = "TF",
                  variable.name = "category",
                  value.name    = "count")
cpgi_long[, category := factor(category,
                                levels = c("not_in_CpGI", "in_CpGI"),
                                labels = c("Not in CpG island", "In CpG island"))]

label_dt <- cpgi[, .(TF, total, pct_in_CpGI)]

p <- ggplot(cpgi_long, aes(x = TF, y = count, fill = category)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(data = label_dt,
            aes(x = TF, y = total, label = paste0(pct_in_CpGI, "%")),
            inherit.aes = FALSE,
            vjust = -0.8, fontface = "bold", size = 5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  scale_fill_manual(values = c("Not in CpG island" = "#B2DFDB",
                               "In CpG island"     = "#00796B")) +
  labs(x = NULL, y = "Number of motifs",
       title = "TF motifs overlapping CpG islands",
       fill  = NULL) +
  theme_bw(base_size = 14) +
  theme(legend.position  = "top",
        panel.grid.major.x = element_blank(),
        plot.title       = element_text(size = 15, hjust = 0.5),
        plot.margin      = margin(t = 20, r = 20, b = 10, l = 10))

ggsave("/data/najemd/TF_binding_cancer/results/integration/motifs_in_cpgi_barplot.pdf",
       plot = p, width = 5, height = 6)
'

# Enrichment test to see if i take all motifs of NRF1 and BANP that overlap with a probe how much falls in a CpG islandIM HERE LOOK AT TH EOUTPUT 
Rscript -e '
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

project_root <- "/data/najemd/TF_binding_cancer"
out_dir      <- file.path(project_root, "results/integration")

bg_raw <- fread(cmd = paste("zcat", file.path(project_root, "methylation/probe_gene_pairs_in_motifs_with_tss.tsv.gz")))
bg_raw[, TF := sub("_.*", "", TF)]
bg_raw <- bg_raw[is_proximal == TRUE]

sig       <- fread(file.path(out_dir, "meth_expr_consistent_pairs_cancermine_cpgi.tsv"))
sig_pairs <- unique(sig[, .(TF, gene)])
bg_data   <- bg_raw[!sig_pairs, on = .(TF, gene)]
cat("[INFO] Sig pairs:", nrow(sig_pairs), "| BG only pairs:", nrow(unique(bg_data[, .(TF, gene)])), "\n")

sig_cpgi <- unique(sig[, .(TF, gene, in_CpGI)])
if (is.character(sig_cpgi$in_CpGI)) sig_cpgi[, in_CpGI := in_CpGI == "yes"]
sig_cpgi <- sig_cpgi[, .(in_CpGI = any(in_CpGI)), by = .(TF, gene)]
cat("[INFO] Sig in CpGI:", sum(sig_cpgi$in_CpGI), "/ total:", nrow(sig_cpgi), "\n")

bg_bed <- unique(bg_data[, .(TF, gene, motif_id = motif_instance_id)])
cat("[INFO] Sample motif_id:", bg_bed$motif_id[1], "\n")
bg_bed[, c("chr","start","end") := tstrsplit(motif_id, ":", keep=1:3)]
tmp_a <- tempfile(fileext=".bed"); tmp_b <- tempfile(fileext=".bed")
fwrite(bg_bed[, .(chr,start,end,motif_id,TF,gene)], tmp_a, sep="\t", col.names=FALSE)
system(paste("bedtools intersect -a", tmp_a, "-b /data/genome/annotations/hg38_cpgIsland.bed.gz -u >", tmp_b))
cpgi_res <- fread(tmp_b, header=FALSE)
cat("[INFO] BG motifs in CpGI:", nrow(cpgi_res), "\n")
bg_bed[, in_CpGI := motif_id %in% (if(nrow(cpgi_res)>0) cpgi_res[[4]] else character(0))]
bg_cpgi <- bg_bed[, .(in_CpGI = any(in_CpGI)), by = .(TF, gene)]
cat("[INFO] BG in CpGI:", sum(bg_cpgi$in_CpGI), "/ total:", nrow(bg_cpgi), "\n")

make_summ <- function(dt, grp) dt[, .(group=grp, n_pairs=.N, n_in_cpgi=sum(in_CpGI), pct_in_cpgi=round(100*sum(in_CpGI)/.N,1)), by=TF]
summ <- rbindlist(list(make_summ(sig_cpgi,"Significant"), make_summ(bg_cpgi,"Background")))
print(summ[order(TF,group)])

fisher_res <- rbindlist(lapply(unique(sig_cpgi$TF), function(tf) {
  s <- sig_cpgi[TF==tf]; b <- bg_cpgi[TF==tf]
  if (nrow(b)==0) return(NULL)
  mat <- matrix(c(sum(s$in_CpGI),sum(!s$in_CpGI),sum(b$in_CpGI),sum(!b$in_CpGI)), nrow=2, byrow=TRUE,
                dimnames=list(c("Significant","Background"),c("In_CpGI","Not_in_CpGI")))
  cat("\n[Fisher] TF =", tf, "\n"); print(mat)
  ft <- fisher.test(mat)
  data.table(TF=tf, n_sig=nrow(s), n_bg=nrow(b), OR=round(ft$estimate,2),
             CI_low=round(ft$conf.int[1],2), CI_high=round(ft$conf.int[2],2),
             p_value=signif(ft$p.value,3),
             sig_label=ifelse(ft$p.value<0.001,"***",ifelse(ft$p.value<0.01,"**",ifelse(ft$p.value<0.05,"*","ns"))))
}))
print(fisher_res)
fwrite(fisher_res, file.path(out_dir,"cpgi_enrichment_fisher.tsv"), sep="\t")

plot_dat <- merge(summ, fisher_res[,.(TF,sig_label,OR,CI_low,CI_high,p_value)], by="TF")
plot_dat[, group := factor(group, levels=c("Background","Significant"))]
star_dat <- plot_dat[, .(y=max(pct_in_cpgi)+8, p_label=paste0(sig_label[1],"\n",
  ifelse(p_value[1]<0.001, paste0("p = ",formatC(p_value[1],format="e",digits=2)), paste0("p = ",formatC(p_value[1],format="f",digits=3))),
  "\nOR = ",OR[1]," [",CI_low[1],"–",CI_high[1],"]")), by=TF]

p <- ggplot(plot_dat, aes(x=group, y=pct_in_cpgi, fill=group)) +
  geom_col(width=0.55) +
  geom_text(aes(label=paste0(pct_in_cpgi,"%\n(n=",n_pairs,")")), vjust=-0.4, size=3.4) +
  geom_text(data=star_dat, aes(x=1.5, y=y, label=p_label), inherit.aes=FALSE, size=3.6, fontface="bold", hjust=0.5) +
  scale_fill_manual(values=c(Background="grey70",Significant="#00796B"), guide="none") +
  facet_wrap(~TF) +
  scale_y_continuous(expand=expansion(mult=c(0,0.25)), breaks=seq(0,100,25), labels=function(x) paste0(x,"%")) +
  labs(title="CpG island enrichment: significant vs background motif-linked pairs",
       subtitle="Per TF: proportion of proximal probe-gene pairs whose motif overlaps a CpG island.\nBackground = all proximal motif-linked pairs NOT in the significant set.",
       x=NULL, y="% pairs with motif in CpG island",
       caption="OR > 1: significant pairs more likely to fall in CpG islands than background.\nFisher'\''s exact test; pairs aggregated at TF–gene level (any motif in CpGI = TRUE).") +
  theme_bw(base_size=13) +
  theme(strip.text=element_text(face="bold"),
        plot.subtitle=element_text(size=8.5,colour="grey20",lineheight=1.3),
        plot.caption=element_text(size=7.5,colour="grey40",hjust=0,lineheight=1.3),
        plot.caption.position="plot", plot.margin=margin(t=15,r=15,b=10,l=10))

ggsave(file.path(out_dir,"cpgi_enrichment_vs_background.pdf"), p, width=9, height=7)
cat("[DONE] cpgi_enrichment_vs_background.pdf\n")
cat("[DONE] cpgi_enrichment_fisher.tsv\n")
'
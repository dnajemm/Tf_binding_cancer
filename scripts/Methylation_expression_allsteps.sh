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
out_dir    <- "./results/integration"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tsv_out <- file.path(out_dir, "meth_expr_filtered_pairs.tsv")

if (file.exists(tsv_out)) {
  cat("TSV already exists, loading and going straight to plots ...\n")
  full_dt <- fread(tsv_out)
} else {

  cat("Loading probe-gene pairs ...\n")
  pg <- fread(pg_file)
  setnames(pg, tolower(colnames(pg)))
  pg <- pg[abs(dist_to_tss) <= PROX_THRESH]
  pg <- unique(pg[, .(tf, probe, gene)])
  cat(sprintf("  %d unique proximal probe-gene pairs\n", nrow(pg)))

  tfs         <- c("BANP", "NRF1")
  all_results <- list()

  for (tf in tfs) {
    cat("\n========", tf, "========\n")

    TF    <- tf
    pg_tf <- unique(pg[grepl(paste0("^", TF, "(_|$)"), pg$tf), .(probe, gene)])
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

      # ── 6. Pearson r per probe-gene pair ──────────────────────────────────
      merged[, pearson_r := NA_real_]
      merged[, pearson_p := NA_real_]

      for (i in seq_len(nrow(merged))) {
        pr   <- merged$probe[i]
        gn   <- merged$gene[i]
        bv   <- beta_long[probe == pr, .(sample_key, beta)]
        ev   <- expr_dt[gene == gn,    .(sample_key, count)]
        both <- merge(bv, ev, by = "sample_key")
        if (nrow(both) < 5) next
        ct <- tryCatch(
          cor.test(both$beta, log2(both$count + 1), method = "pearson"),
          error = function(e) NULL
        )
        if (!is.null(ct)) {
          merged[i, pearson_r := as.numeric(ct$estimate)]
          merged[i, pearson_p := ct$p.value]
        }
      }

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


# Script to generate scatter plots for each gene in tsv file it takes the list of genes that are in the first column of the tsv file: 
# ============================================================
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
stats_out  <- "./results/integration/meth_expr_wilcoxon_NRF1_stats.tsv"

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

########################
# Add a Score to each gene-cancer pair based on the correlation and the Wilcoxon test results
########################
Rscript - << 'EOF'
suppressPackageStartupMessages(library(data.table))

# ── Config ─────────────────────────────────────────────────────────────────────
tf_name    <- "NRF1"   # change to NRF1 or loop
pairs_tsv  <- "./results/integration/meth_expr_filtered_pairs.tsv"
wilcox_tsv <- paste0("./results/integration/meth_expr_wilcoxon_", tf_name, "_stats.tsv")
out_tsv    <- paste0("./results/integration/gene_ranking_", tf_name, ".tsv")

# ── Load ───────────────────────────────────────────────────────────────────────
pairs  <- fread(pairs_tsv)[TF == tf_name]
wilcox <- fread(wilcox_tsv)

# ── Merge pairs + wilcoxon stats ───────────────────────────────────────────────
dt <- merge(pairs, wilcox[, .(TF, gene, cancer, probe,
                               n_tumor, n_healthy,
                               expr_wilcox_p, meth_wilcox_p)],
            by = c("TF", "gene", "cancer", "probe"), all.x = TRUE)

# ── Filter: both wilcoxon significant ─────────────────────────────────────────
dt <- dt[!is.na(expr_wilcox_p) & !is.na(meth_wilcox_p) &
          expr_wilcox_p < 0.05  & meth_wilcox_p < 0.05]

cat(sprintf("Pairs after wilcoxon filter: %d | unique genes: %d\n",
            nrow(dt), uniqueN(dt$gene)))

# ── Helper: normalize vector to 0-1 ───────────────────────────────────────────
norm01 <- function(x) {
  r <- range(x, na.rm = TRUE)
  if (diff(r) == 0) return(rep(0, length(x)))
  (x - r[1]) / diff(r)
}

# ── Score per pair ─────────────────────────────────────────────────────────────
# Each component normalized 0-1, then weighted
dt[, score :=
  2 * norm01(abs(pearson_r))           +
  norm01(-log10(pearson_p))        +
  2 * norm01(abs(log2foldchange))      +
  norm01(-log10(padj))             +
  2 * norm01(abs(mean_delta_beta))     +
  norm01(-log10(adj.P.Val))        +
  norm01(-log10(expr_wilcox_p))    +
  norm01(-log10(meth_wilcox_p))
]

# ── Collapse per gene across cancers ───────────────────────────────────────────
gene_scores <- dt[, .(
  n_cancers       = .N,
  best_score      = max(score),
  mean_score      = round(mean(score), 4),
  cancers         = paste(sort(cancer), collapse = ", "),
  mean_pearson_r  = round(mean(pearson_r), 3),
  mean_delta_beta = round(mean(mean_delta_beta), 2),
  mean_log2FC     = round(mean(log2foldchange), 3),
  best_pearson_p  = min(pearson_p),
  best_padj       = min(padj),
  best_adj_P_Val  = min(adj.P.Val),
  best_expr_wilcox_p = min(expr_wilcox_p),
  best_meth_wilcox_p = min(meth_wilcox_p)
), by = gene][order(-best_score)]

# bonus: genes recurring in multiple cancers get a bump
gene_scores[, final_score := best_score * (1 + 0.2 * (n_cancers - 1))]
setorder(gene_scores, -final_score)

cat(sprintf("\nTop 20 genes for %s:\n", tf_name))
print(gene_scores[1:min(20, .N), .(gene, n_cancers, final_score, cancers,
                                    mean_pearson_r, mean_delta_beta, mean_log2FC)])

fwrite(gene_scores, out_tsv, sep = "\t")
cat(sprintf("\nFull ranking written to: %s\n", out_tsv))
EOF

###################################
# Gene enrichment GO analysis for the genes
###################################
# 1) Extract gene lists
awk -F'\t' '{print $1}' ./results/integration/gene_ranking_NRF1.tsv | tail -n +2 | sort -u > ./results/integration/NRF1_genes_cancer.txt
awk -F'\t' '{print $1}' ./results/integration/gene_ranking_BANP.tsv | tail -n +2 | sort -u > ./results/integration/BANP_genes_cancer.txt

NRF1_file="./results/integration/NRF1_genes_cancer.txt" 
BANP_file="./results/integration/BANP_genes_cancer.txt" 


mkdir -p ./results/integration/GO_plots/
mkdir -p ./results/integration/GO_plots/NRF1/plots
mkdir -p ./results/integration/GO_plots/BANP/plots
# 2) Run enrichment for both TFs
# 2) run the enrichment pipeline for each of the gene lists:

file="./results/integration/BANP_genes_cancer.txt"
output_dir="./results/integration/GO_plots/BANP"
column=1
species="hsa"
format="name"
pvalue=0.001
database="go_biological_process"

while getopts "hf:s:n:c:p:o:d:" opt; do
    case $opt in
	h)
	    echo -e $help >&2
	    exit
	    ;;
	f)
	    file=$OPTARG
	    ;;
	s)
	    species=$OPTARG
	    ;;
	n)
	    format=$OPTARG
	    ;;
	c)
	    column=$OPTARG
	    ;;
	p)
	    pvalue=$OPTARG
	    ;;
	o)
	    output_dir=$OPTARG
	    ;;
	d)
	    database=$OPTARG
	    ;;
	\?)
            echo -e $help >&2
            exit
	    ;;
    esac
done

# Check file name and species are provided
if [[ -z $file || -z $species ]]; then echo -e $help >&2; exit; fi

# Check file exist and is not empty
if [[ ! -s $file ]]; then echo ${file}" file does not exist or is empty" >&2; exit 1; fi

# Check species is in the list
if [[ $species != "hsa" && $species != "mmu" && $species != "bta" && $species != "cfa" && $species != "dre" && $species != "gga" && $species != "ssc" ]]; then echo ${species}" species is not supported" >&2; echo -e $help >&2; exit 1; fi

# Assign species to DB file
if [[ $database = "go_biological_process" ]]; then
    DB=/data/genome/annotations/${species}_go_biological_process.txt
elif [[ $database = "wikipathways" ]]; then
    DB=/data/genome/annotations/${species}_wikipathways.txt
elif [[ $database = "tf_dbd" ]]; then
    DB=/data/genome/annotations/${species}_tf_dbd.txt
else echo $database" database not available"; exit 1
fi
if [[ ! -s $DB ]]; then echo "Database file not available: "$DB; exit 1; fi
						   

# Check column exists in file
if [[ ! $column =~ ^[0-9]+$ ]]; then echo "COLUMN "$column" is not a number"; exit 1; fi
status=$(awk -vN=$column '(N>NF){print "KO"}' $file | uniq)
if [[ $status = "KO" ]]; then echo $file" file does not contain column "$column" on all lines"; exit 1; fi

# Check p-value is between 0 and 1
status=$(awk -vP=$pvalue 'BEGIN{if(P<0||P>1){print "KO"}}')
if [[ $pvalue = "KO" ]]; then echo "P-value threshold should be between 0 and 1"; exit 1; fi

# TMP file containing gene names
if [[ $format != "name" && $format != "id" ]]; then echo -e $help >&2; exit 1; fi
TMP=$(mktemp -u -p /scratch/$USER/tmp/)
if [[ $format = "name" ]]; then
    cat $file | awk -vN=$column '{print $N}' | sort -u > $TMP
# Convert ensembl IDs to gene names
elif [[ $format = "id" ]]; then
    if [[ $species = "hsa" ]]; then gene_file=/data/genome/annotations/hg38_genes.bed
    elif [[ $species = "mmu" ]]; then gene_file=/data/genome/annotations/mm10_genes.bed
    elif [[ $species = "bta" ]]; then gene_file=/data/genome/annotations/bosTau8_genes.bed
    elif [[ $species = "cfa" ]]; then gene_file=/data/genome/annotations/canFam3_genes.bed
    elif [[ $species = "dre" ]]; then gene_file=/data/genome/annotations/danRer10_genes.bed
    elif [[ $species = "gga" ]]; then gene_file=/data/genome/annotations/galGal5_genes.bed
    elif [[ $species = "ssc" ]]; then gene_file=/data/genome/annotations/susScr11_genes.bed
    fi
    cat $file | awk -vN=$column '{print $N}' | sort -u | awk -vF=$gene_file -vTMP=$TMP -vFS='\t' 'BEGIN{while(getline<F){G[$6]=$4}}{if(G[$1]){print G[$1]>TMP;c++}}END{print c" IDs could be converted to gene names ("int(c*100/NR)"%)"}'
fi

echo $(date)": gene_enrichment start processing "$file

# Check gene names exist
total=$(cat $TMP | sort -u | wc -l)
count=$(cat $TMP | awk -vDB=$DB 'BEGIN{while(getline<DB){G[$1]=1}}{if($1 in G){c++}}END{print c+0}')
if [[ $count -eq 0 ]]; then  echo "None of your genes are present in the current annotation database "$DB" (check gene name/id format)"; exit 1;
else echo $count" genes are present in the current annotation database ("$(awk -vC=$count -vT=$total 'BEGIN{print int(C*100/T)}')"%)";fi

# Output directory
if [[ -n $output_dir ]]; then mkdir -p $output_dir
else output_dir=$(echo $file | awk '{if($1!~"/"){print "."}else{split($1,X,"/");L=X[1];for(i=2;i<=length(X)-1;i++){L=L"/"X[i]};print L}}')
fi

# Output file
output_file=$(echo $file | awk '{split($1,X,"/");split(X[length(X)],Y,".");print Y[1]}')

# Run annotation enrichment
(echo -e "#db_category\tgenes_in_category\tgenes\tdb_category\tdb\tp-value\tcategory_name\tgene_list"
 cat $TMP | awk -vF=$DB -vOFS='\t' '{G[$1]=1}END{while(getline<F){DBt[$2]++;FUNC[$2]=$3;Gt[$1]=1;if($1 in G){Gg[$1]=1;DBg[$2]++;if(!Lg[$2]){Lg[$2]=$1}else{Lg[$2]=Lg[$2]","$1}}};for(db in DBt){if(!Lg[db]){Lg[db]="."};print db,DBg[db]+0,length(Gg)+0,DBt[db]+0,length(Gt)+0,FUNC[db],Lg[db]}}' | Rscript /data/software/dorsal/scripts/hyper.R -i - | awk -vOFS='\t' -vP=$pvalue '($8<=P){print $1,$2,$3,$4,$5,$8,$6,$7}' | sort -k6,6g -k2,2nr
)> ${output_dir}/${output_file}"_"$database"_table.txt"

# Number of categories enriched
NB=$(awk '(NR>1)' ${output_dir}/${output_file}"_"$database"_table.txt" | wc -l)
echo $NB" categories significantly enriched (p-value <= "$pvalue") (top ones below)"
if [[ $NB -gt 0 ]]; then 
    cat ${output_dir}/${output_file}"_"$database"_table.txt" | head | column -t
else
    rm ${output_dir}/${output_file}"_"$database"_table.txt"
fi

# Plot results
if [[ $NB -ge 20 ]]; then NB=20; fi
if [[ $NB -gt 0 ]]; then 
    Rscript /data/bardeta/scripts/gene_enrichment_plot.R -i ${output_dir}/${output_file}"_"$database"_table.txt" -n $NB -o ${output_dir}/plots/${output_file}"_"$database"_barplot.pdf"
    echo "Plot generated: display "${output_dir}/plots/${output_file}"_"$database"_barplot.pdf"
fi
rm $TMP

echo $(date)": gene_enrichment processed "$file" successfully"

# NOT DONE FOR ALL OF THEM YET
##################################################################################################################################################################################################################
# Same integration methylation_expression pipeline as before just here this is to get one pdf page per gene-cancer pair with the scatter plot and the boxplots for methylation and expression
# to have the gene with multiple cancer type detected in the same page in the same pdf for the report file
##################################################################################################################################################################################################################
Rscript - << 'EOF'
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
  library(gtable)
  library(parallel)
})

# ── USER CONFIG ───────────────────────────────────────────────────────────────
target_gene    <- "CAMK2A"
target_probe   <- NULL
target_cancers <- c("BLCA","BRCA","KIRP","UCEC")

tf_name  <- "BANP"
N_COL    <- 3L
PDF_W    <- 100
PDF_H    <- 75
N_CORES  <- max(1L, detectCores() - 2L)

pairs_tsv <- "./results/integration/meth_expr_filtered_pairs.tsv"
beta_dir  <- paste0("./methylation/beta_matrices_", tf_name)
tpm_file  <- "./expression/gene_expression_matrix_protein_coding.tsv"
out_dir   <- file.path("./results/integration/gene_cancer_plots_BANP/")
stats_out <- paste0("./results/integration//gene_cancer_plots_BANP/meth_expr_wilcoxon_", tf_name,
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
# FIX 2: only 01A_1 (tumor) and 11A_1 (healthy)
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

  # count total samples in matrix before any filtering
  n_matrix_samples <- ncol(beta_mat) - 1L

  beta_mat <- beta_mat[probe == target_probe]
  if (nrow(beta_mat) == 0) { cat(tag, "– probe not in matrix\n"); return(NULL) }
  beta_long <- melt(beta_mat, id.vars = "probe",
                    variable.name = "col_name", value.name = "meth_beta")

  # count NAs for this probe across all samples
  n_na_probe <- sum(is.na(beta_long$meth_beta))

  # FIX 1: columns are CANCER-TCGA-XX-XXXX_tumor, not TCGA-XX-XXXX_tumor
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

# ── Draw grid PDF ─────────────────────────────────────────────────────────────
out_pdf <- file.path(out_dir,
  paste0(tf_name, "__", target_gene, "__", target_probe, "__all_cancers_grid.pdf"))
cat("Writing PDF:", out_pdf, "\n")

pdf(out_pdf, width = PDF_W, height = PDF_H, onefile = FALSE)
grid.newpage()
n_rows_outer <- ceiling(n_valid / N_COL)

pushViewport(viewport(layout = grid.layout(
  nrow    = n_rows_outer + 1L,
  ncol    = N_COL,
  heights = unit.c(unit(1.5, "cm"),
                   rep(unit(1 / n_rows_outer, "null"), n_rows_outer)),
  widths  = rep(unit(1 / N_COL, "null"), N_COL)
)))

grid.text(
  sprintf("%s — %s | %s   (%d cancer type%s)",
          tf_name, target_gene, target_probe, n_valid, ifelse(n_valid==1,"","s")),
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
dev.off()
cat("Done. PDF written to:", out_pdf, "\n")

# ── Write stats ───────────────────────────────────────────────────────────────
stats_dt <- rbindlist(Filter(Negate(is.null), stats_list), fill = TRUE)
if (nrow(stats_dt) > 0) {
  fwrite(stats_dt, stats_out, sep = "\t")
  cat("Stats written to:", stats_out, "\n")
}
EOF
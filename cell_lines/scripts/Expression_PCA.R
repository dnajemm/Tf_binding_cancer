#!/usr/bin/Rscript --vanilla

## PCA — TCGA samples (matched to beta matrix) + cell lines
##
## Usage:
##   Rscript Expression_PCA.R -t TCGA_TPM -d CELLLINE_DIR -c CELLLINES -k CANCER_TYPE [-o OUTDIR] [-n NGENES] [--tumor-only]
##
## -t           : TCGA TPM matrix — .tsv/.tsv.gz (auto-transposed) OR pre-cached .rds (faster)
## -d           : Directory containing cell line *_genes_count.txt files
## -c           : Cell line names to include, comma-separated (e.g. "HCC1954,HMEC")
## -k           : Cancer type (e.g. BRCA) — used to find beta matrix and filter TCGA samples
## -o           : Output directory (default: .)
## -n           : Number of top variable genes for PCA (default: all)
## -f           : Output filename (default: CANCERTYPE_pca_tpm.pdf)
## -s           : Optional subtypes file (e.g. BRCA_subtypes.txt with Sample + PAM50 columns)
##                If provided, TCGA tumor samples are colored by subtype instead of "tumor"
## --tumor-only : Keep only tumor TCGA samples (drop healthy)

args <- commandArgs(TRUE)

help <- function(){
    cat("\nTPM PCA — TCGA (matched pairs) + cell lines\n")
    cat("Usage: Rscript Expression_PCA.R -t TCGA_TPM -d CELLLINE_DIR -c CELLLINES -k CANCER_TYPE [-o OUTDIR] [-n NGENES] [--tumor-only]\n\n")
    q()
}

if(length(args)==0 || !is.na(charmatch("-help", args))) help()

tumor_only <- "--tumor-only" %in% args
args       <- args[args != "--tumor-only"]

for(ii in seq_along(args)){
    if(grepl("^-", args[ii]) && args[ii] != "-"){
        if(ii+1 <= length(args) && (!grepl("^-", args[ii+1]) || args[ii+1] == "-"))
            assign(gsub("-","", args[ii]), args[ii+1])
        else
            assign(gsub("-","", args[ii]), 1)
    }
}

if(!exists("t")) { cat("ERROR: -t (TCGA TPM matrix) is required\n"); q() }
if(!exists("d")) { cat("ERROR: -d (cell line directory) is required\n"); q() }
if(!exists("c")) { cat("ERROR: -c (cell line names, e.g. 'HCC1954,HMEC') is required\n"); q() }
if(!exists("k")) { cat("ERROR: -k (cancer type, e.g. BRCA) is required\n"); q() }
if(!exists("o")) o <- "."
if(!exists("n")) n <- 99999 else n <- as.integer(n)
dir.create(o, recursive=TRUE, showWarnings=FALSE)

cancer_type  <- k
project_root <- "/data/najemd/TF_binding_cancer"
cl_wanted    <- trimws(strsplit(c, ",")[[1]])

cat("Cancer type :", cancer_type, "\n")
cat("Cell lines  :", paste(cl_wanted, collapse=", "), "\n")
cat("Tumor only  :", tumor_only, "\n\n")

suppressMessages({
    library(ggplot2)
    library(ggrepel)
    library(data.table)
})

## ── 1. Find beta matrix & get matched sample list ─────────────────────────────
cat("Finding beta matrix for", cancer_type, "...\n")
beta_file <- file.path(project_root, "methylation/beta_matrices_BANP",
                       paste0(cancer_type, "_BANP_beta_matrix.tsv.gz"))
if(!file.exists(beta_file)){
    beta_file <- file.path(project_root, "methylation/beta_matrices_NRF1",
                           paste0(cancer_type, "_NRF1_beta_matrix.tsv.gz"))
}
if(!file.exists(beta_file)){
    cat("WARNING: no beta matrix found — will use ALL TCGA samples for this cancer type\n")
    beta_dt <- NULL
} else {
    cat("  Using:", beta_file, "\n")
    beta_header  <- fread(beta_file, nrows=0)
    beta_samples <- setdiff(names(beta_header), names(beta_header)[1])
    beta_dt <- data.table(col=beta_samples)
    beta_dt[, patient  := sub("_[a-z]+$", "", col)]
    beta_dt[, st       := fifelse(grepl("_tumor$",   col), "Tumor",
                         fifelse(grepl("_healthy$", col), "Healthy", "other"))]
    beta_dt <- beta_dt[st != "other"]
    beta_dt[, pair_key := paste(patient, st, sep="__")]
    cat(paste0("  Beta matrix: ", sum(beta_dt$st=="Tumor"), " tumor | ",
               sum(beta_dt$st=="Healthy"), " healthy matched samples\n"))
}

## ── 2. Load TCGA TPM matrix ───────────────────────────────────────────────────
cat("\nLoading TCGA TPM matrix...\n")

rds_cache <- sub("\\.(tsv|tsv\\.gz)$", ".rds", t)

if(grepl("\\.rds$", t)){
    # Fast path: pre-cached RDS (genes x samples)
    cat("  Loading from RDS cache...\n")
    tcga <- readRDS(t)
} else {
    # Slow path: read TSV, transpose, then save RDS for next time
    tcga_raw <- read.table(t, header=TRUE, sep="\t", row.names=1, check.names=FALSE)
    if(ncol(tcga_raw) > nrow(tcga_raw)){
        cat("  Detected samples-as-rows format — transposing...\n")
        tcga <- as.data.frame(t(tcga_raw))
    } else {
        tcga <- tcga_raw
    }
    rm(tcga_raw)
    # Save RDS cache for future runs
    cat(paste0("  Saving RDS cache to: ", rds_cache, " (use -t ", rds_cache, " next time)\n"))
    saveRDS(tcga, rds_cache)
}
cat(paste0("  ", nrow(tcga), " genes x ", ncol(tcga), " TCGA samples (before filtering)\n"))

## ── 3. Filter TCGA samples ────────────────────────────────────────────────────
tcga_samples_all <- colnames(tcga)

tcga_cancer  <- sub("^TCGA-([A-Z]+)-.*", "\\1", tcga_samples_all)
tcga_st      <- ifelse(grepl("-01A_1", tcga_samples_all), "Tumor",
               ifelse(grepl("-11A_1", tcga_samples_all), "Healthy", "other"))
tcga_patient <- sub("^TCGA-([A-Z]+-TCGA-[A-Z0-9]+-[A-Z0-9]+)-[0-9]+[A-Z]_.*",
                    "\\1", tcga_samples_all)
tcga_pk      <- paste(tcga_patient, tcga_st, sep="__")

if(!is.null(beta_dt)){
    keep <- tcga_cancer == cancer_type &
            tcga_st     != "other"     &
            tcga_pk     %in% beta_dt$pair_key
} else {
    keep <- tcga_cancer == cancer_type & tcga_st != "other"
}

# Apply --tumor-only flag
if(tumor_only){
    keep <- keep & tcga_st == "Tumor"
    cat("  --tumor-only: keeping tumor samples only\n")
}

tcga         <- tcga[, keep, drop=FALSE]
tcga_st_kept <- tcga_st[keep]
cat(paste0("  After filtering: ",
           sum(tcga_st_kept=="Tumor"),   " tumor | ",
           sum(tcga_st_kept=="Healthy"), " healthy TCGA samples kept\n"))

if(ncol(tcga) == 0){ cat("ERROR: no TCGA samples passed filtering\n"); q() }

## ── 4. Load & merge cell line files ──────────────────────────────────────────
cat("\nLoading cell line files from:", d, "\n")
cat("  Requested:", paste(cl_wanted, collapse=", "), "\n")

pattern  <- paste0("(", paste(cl_wanted, collapse="|"), ").*_genes_count\\.txt$")
cl_files <- list.files(d, pattern=pattern, full.names=TRUE)
if(length(cl_files) == 0){
    cat("ERROR: no files found for:", paste(cl_wanted, collapse=", "), "in", d, "\n"); q()
}

cl_names <- sub("^RNA_", "", sub("_genes_count\\.txt$", "", basename(cl_files)))
names(cl_files) <- cl_names

cat("  Files matched:\n")
for(nm in cl_names) cat("    -", nm, "\n")

cl_list <- lapply(cl_names, function(nm){
    f  <- cl_files[nm]
    d2 <- read.table(f, header=TRUE, sep="\t", stringsAsFactors=FALSE, comment.char="")
    d2 <- d2[order(d2$TPM, decreasing=TRUE), ]
    d2 <- d2[!duplicated(d2$geneName), ]
    tpm_vec <- d2$TPM
    names(tpm_vec) <- d2$geneName
    tpm_vec
})
names(cl_list) <- cl_names

cl_common <- Reduce(intersect, lapply(cl_list, names))
cl_mat    <- do.call(cbind, lapply(cl_list, function(x) x[cl_common]))
rownames(cl_mat) <- cl_common
cat(paste0("  Cell line matrix: ", nrow(cl_mat), " genes x ", ncol(cl_mat), " samples\n"))

## ── 5. Find common genes & combine ───────────────────────────────────────────
cat("\nFinding common genes...\n")
common_genes <- intersect(rownames(tcga), rownames(cl_mat))
cat(paste0("  ", length(common_genes), " genes in common\n"))
if(length(common_genes) < 100)
    cat("WARNING: very few common genes — check both use gene symbols\n")

mat <- cbind(tcga[common_genes, ], cl_mat[common_genes, ])
cat(paste0("  Combined matrix: ", nrow(mat), " genes x ", ncol(mat), " samples\n"))

## ── 6. Annotate samples ───────────────────────────────────────────────────────
cl_samples   <- cl_names
tcga_samples <- colnames(tcga)

# Base TCGA group: tumor / healthy
tcga_group <- ifelse(grepl("-01[A-Z0-9]", tcga_samples), "tumor",
              ifelse(grepl("-11[A-Z0-9]", tcga_samples), "healthy", "tumor"))
names(tcga_group) <- tcga_samples

# Optional PAM50 subtypes — replace "tumor" with subtype label
subtype_levels <- NULL
if(exists("s") && file.exists(s)){
    cat("\nLoading PAM50 subtypes from:", s, "\n")
    raw_lines <- readLines(s)
    raw_lines <- gsub("\r", "", raw_lines)
    raw_lines <- raw_lines[nchar(raw_lines) > 0]
    subtypes  <- read.table(text=paste(raw_lines, collapse="\n"),
                            header=TRUE, sep="\t", stringsAsFactors=FALSE)
    # Extract standard TCGA barcode from our sample names (e.g. TCGA-A1-A0SK-01A)
    extract_std_barcode <- function(x){
        m <- regmatches(x, gregexpr("TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}[A-Z]", x))
        sapply(m, function(v) if(length(v) >= 1) v[length(v)] else NA_character_)
    }
    std_barcodes         <- extract_std_barcode(tcga_samples)
    pam50_lookup         <- setNames(subtypes$PAM50, subtypes$Sample)
    pam50_for_samples    <- pam50_lookup[std_barcodes]
    # Only replace tumor samples that have a subtype match
    has_subtype          <- !is.na(pam50_for_samples) & tcga_group == "tumor"
    tcga_group[has_subtype] <- pam50_for_samples[has_subtype]
    subtype_levels <- c("LumA","LumB","Basal","Her2","Normal")
    cat(paste0("  PAM50 subtypes matched: ", sum(has_subtype), " / ",
               sum(tcga_group == "tumor" | has_subtype), " tumor samples\n"))
    cat("  Subtype distribution:\n"); print(table(tcga_group[has_subtype]))
}

cl_group        <- sub("_[0-9]+$", "", cl_samples)
names(cl_group) <- cl_samples

group       <- c(tcga_group, cl_group)
source_type <- ifelse(colnames(mat) %in% cl_samples, "cell_line", "TCGA")

cat("\nSample group counts:\n")
print(table(group))

## ── 7. Log2 transform & top variable genes ───────────────────────────────────
cat(paste0("\nLog2(TPM+1) + top ", n, " variable genes...\n"))
log_mat   <- log2(mat + 1)
gene_sd   <- apply(log_mat, 1, sd)
gene_sd   <- gene_sd[gene_sd > 0]          # remove zero-variance genes
top_genes <- names(sort(gene_sd, decreasing=TRUE))[1:min(n, length(gene_sd))]
log_mat   <- log_mat[top_genes, ]
cat(paste0("  Using ", length(top_genes), " genes for PCA\n"))

## ── 8. PCA ────────────────────────────────────────────────────────────────────
cat("Running PCA...\n")
pca     <- prcomp(t(log_mat), scale.=TRUE)
pct_var <- round(100 * summary(pca)$importance[2, 1:2], 1)
cat(paste0("  PC1: ", pct_var[1], "%  PC2: ", pct_var[2], "%\n"))

pca_df <- data.frame(
    PC1    = pca$x[, 1],
    PC2    = pca$x[, 2],
    group  = factor(group[colnames(mat)]),
    source = factor(source_type, levels=c("TCGA","cell_line")),
    sample = colnames(mat)
)

## ── 9. Colors & shapes ────────────────────────────────────────────────────────
cl_line_names <- unique(sub("_[0-9]+$", "", cl_samples))
cl_palette    <- c("#2CA86A","#F5A623","#9B59B6","#1ABC9C","#E67E22")
cl_colors     <- setNames(cl_palette[seq_along(cl_line_names)], cl_line_names)

if(!is.null(subtype_levels)){
    # PAM50 color palette (CVD-safe)
    pam50_colors <- c(
        LumA    = "#2a78d6",
        LumB    = "#eb6834",
        Basal   = "#1baf7a",
        Her2    = "#e87ba4",
        Normal  = "#008300",
        tumor   = "#E05C5C",   # unmatched tumors
        healthy = "#4A90D9"
    )
    colors <- c(pam50_colors, cl_colors)
} else {
    colors <- c(tumor="#E05C5C", healthy="#4A90D9", cl_colors)
}
shapes <- c(TCGA=16, cell_line=17)

## ── 10. Subtitle with sample counts ──────────────────────────────────────────
grp_counts   <- table(group[colnames(mat)])
n_tumor      <- grp_counts["tumor"]
n_healthy    <- grp_counts["healthy"]
cl_counts    <- sapply(cl_line_names, function(cl)
                   sum(grepl(paste0("^", cl), cl_samples)))
cl_count_str <- paste(paste0(cl_line_names, " n=", cl_counts), collapse=" | ")
tcga_str <- if(tumor_only){
    paste0("TCGA ", cancer_type, ": tumor n=", n_tumor)
} else {
    paste0("TCGA ", cancer_type, ": tumor n=", n_tumor, " | healthy n=", n_healthy)
}
subtitle_str <- paste0(tcga_str, "    Cell lines: ", cl_count_str)

## ── 11. Plot ──────────────────────────────────────────────────────────────────
p <- ggplot(pca_df, aes(x=PC1, y=PC2, color=group, shape=source)) +
    geom_point(data=subset(pca_df, source=="TCGA"),      size=1.8, alpha=0.65) +
    geom_point(data=subset(pca_df, source=="cell_line"), size=4,   alpha=1   ) +
    geom_text_repel(
        data        = subset(pca_df, source=="cell_line"),
        aes(label   = sample),
        size        = 3,
        show.legend = FALSE,
        max.overlaps= 20
    ) +
    scale_color_manual(values=colors, name="Group") +
    scale_shape_manual(values=shapes, name="Source") +
    xlab(paste0("PC1: ", pct_var[1], "% variance")) +
    ylab(paste0("PC2: ", pct_var[2], "% variance")) +
    ggtitle(paste0("Expression PCA — TCGA ", cancer_type, " + cell lines"),
            subtitle=subtitle_str) +
    theme_bw(base_size=13) +
    theme(
        legend.position  = "right",
        panel.grid.minor = element_blank(),
        plot.subtitle    = element_text(size=9, color="grey40")
    ) +
    guides(
        color = guide_legend(override.aes=list(size=3.5)),
        shape = guide_legend(override.aes=list(size=3.5))
    )

if(!exists("f")) f <- paste0(cancer_type, "_pca_tpm.pdf")
out_pdf <- file.path(o, f)
pdf(out_pdf, width=8, height=6.5)
print(p)
dev.off()
cat(paste0("\n=> ", out_pdf, "\n"))
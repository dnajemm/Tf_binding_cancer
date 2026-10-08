#!/usr/bin/Rscript --vanilla
######
# BRCA PCA colored by PAM50 subtype
# Plain PCA from raw count files — no DESeq2
######

## ── Arguments ──────────────────────────────────────────────────────────────
args = commandArgs(TRUE)

help <- function(){
    cat("\nBRCA PCA with PAM50 subtype colors\n")
    cat("Usage: BRCA_PCA_PAM50.R -s SUBTYPES_FILE -c COUNTS_DIR [-o OUTPUT_FILE] [-t TOP_GENES]\n")
    cat("-s : BRCA_subtypes.txt (col4 = PAM50)\n")
    cat("-c : Directory containing *_genes_count.txt files\n")
    cat("-o : Output PDF path (default: BRCA_PCA_PAM50.pdf)\n")
    cat("-g : Number of most-variable genes to use for PCA (default: 5000)\n")
    cat("\n")
    q()
}

if(length(args) == 0 || !is.na(charmatch("-help", args))){
    help()
} else {
    for(ii in 1:length(args)){
        if(grepl("^-", args[ii]) && args[ii] != "-"){
            if(ii+1 <= length(args) && (!grepl("^-", args[ii+1]) || args[ii+1] == "-")){
                assign(gsub("-","",args[ii]), args[ii+1])
            }
        }
    }
}

if(!exists("s")){ cat("No subtypes file specified (-s)\n"); q() }
if(!exists("c")){ cat("No counts directory specified (-c)\n"); q() }
if(!exists("o")) o <- "BRCA_PCA_PAM50.pdf"
top_genes <- if(exists("g")) as.integer(g) else 5000

# Column index for TPM in augmented STAR files:
# gene_id | gene_name | gene_type | counts_unstranded | counts_fwd | counts_rev | TPM_unstranded | TPM_fwd | TPM_rev
TPM_COL <- 7   # unstranded TPM — change to 8 (fwd) or 9 (rev) for stranded libraries

suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))

## ── Load PAM50 subtypes ────────────────────────────────────────────────────
cat("Loading PAM50 subtypes...\n")

raw_lines <- readLines(s)
raw_lines <- gsub("\r", "", raw_lines)          # fix CR-only line endings
raw_lines <- raw_lines[nchar(raw_lines) > 0]

subtypes <- read.table(
    text       = paste(raw_lines, collapse = "\n"),
    header     = TRUE,
    sep        = "\t",
    stringsAsFactors = FALSE
)
# Keep Sample (full barcode) and PAM50
subtypes <- subtypes[, c("Sample", "bcr_patient_barcode", "PAM50")]
cat(sprintf("  %d samples in subtypes file\n", nrow(subtypes)))

## ── Find BRCA count files ──────────────────────────────────────────────────
cat("Scanning counts directory...\n")
all_files <- list.files(c, pattern = "\\.augmented_star_gene_counts_filtered\\.tsv$", full.names = TRUE)
cat(sprintf("  %d count files found\n", length(all_files)))

# Extract TCGA barcode from filename
# Expected pattern: RNA_TCGA-BRCA-TCGA-A1-A0SK-01A_1_genes_count.txt
#                or TCGA-A1-A0SK-01A_genes_count.txt
extract_barcode <- function(f){
    base <- basename(f)
    base <- sub("\\.augmented_star_gene_counts_filtered\\.tsv$", "", base)
    # Match the last TCGA barcode pattern in the name (TCGA-XX-XXXX-01A)
    m <- regmatches(base, gregexpr("TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}[A-Z]", base))[[1]]
    if(length(m) >= 1) return(m[length(m)])  # take last match
    return(base)
}

file_barcodes <- sapply(all_files, extract_barcode)

# Keep BRCA tumor (01A) matched to subtypes + all BRCA normal adjacent (11A)
in_subtypes  <- file_barcodes %in% subtypes$Sample
is_normal_adj <- grepl("-11[A-Z]$", file_barcodes)
keep_files   <- in_subtypes | is_normal_adj

brca_files    <- all_files[keep_files]
brca_barcodes <- file_barcodes[keep_files]

cat(sprintf("  %d files kept: %d PAM50-annotated tumors + %d normal adjacent (11A)\n",
            length(brca_files), sum(in_subtypes), sum(is_normal_adj)))
if(length(brca_files) == 0){
    cat("ERROR: No files found. Check filenames and subtypes file.\n")
    cat("Example file barcode extracted:", head(file_barcodes, 3), "\n")
    q()
}

## ── Build TPM matrix (genes x samples) ────────────────────────────────────
# File format: gene_id | gene_name | gene_type | counts_unstranded | counts_fwd | counts_rev | TPM_unstranded | TPM_fwd | TPM_rev
cat("Loading TPM values from augmented STAR files...\n")

tpm_list <- vector("list", length(brca_files))
for(i in seq_along(brca_files)){
    x <- read.table(brca_files[i], header = FALSE, sep = "\t",
                    stringsAsFactors = FALSE)
    tpm_list[[i]] <- x[, TPM_COL]
}

# Gene IDs from the first file
gene_ids <- read.table(brca_files[1], header = FALSE, sep = "\t")[, 1]

tpm_mat <- do.call(cbind, tpm_list)
rownames(tpm_mat) <- gene_ids
colnames(tpm_mat) <- brca_barcodes

cat(sprintf("  Matrix: %d genes x %d samples\n", nrow(tpm_mat), ncol(tpm_mat)))

## ── Normalise: log2(TPM + 1) ──────────────────────────────────────────────
cat("Applying log2(TPM + 1)...\n")

cat(sprintf("  Using all %d genes (no filtering)\n", nrow(tpm_mat)))

log_mat <- log2(tpm_mat + 1)

## ── Select top variable genes ──────────────────────────────────────────────
gene_vars   <- apply(log_mat, 1, var)
top_idx     <- order(gene_vars, decreasing = TRUE)[1:min(top_genes, nrow(log_mat))]
log_mat_top <- log_mat[top_idx, ]
cat(sprintf("  Using top %d most-variable genes for PCA\n", nrow(log_mat_top)))

## ── PCA ───────────────────────────────────────────────────────────────────
cat("Running PCA...\n")
pca    <- prcomp(t(log_mat_top), center = TRUE, scale. = TRUE)
pct    <- round(100 * summary(pca)$importance[2, 1:2], 1)  # % variance PC1, PC2

pca_df <- data.frame(
    Sample = colnames(log_mat_top),
    PC1    = pca$x[, 1],
    PC2    = pca$x[, 2],
    stringsAsFactors = FALSE
)

## ── Attach PAM50 labels ────────────────────────────────────────────────────
pca_df$PAM50 <- subtypes$PAM50[match(pca_df$Sample, subtypes$Sample)]

# Label 11A samples as "Normal adjacent" instead of Unknown
is_11a <- grepl("-11[A-Z]$", pca_df$Sample)
pca_df$PAM50[is_11a] <- "Normal adjacent"
pca_df$PAM50[is.na(pca_df$PAM50)] <- "Unknown"

pca_df$PAM50 <- factor(pca_df$PAM50,
                        levels = c("LumA","LumB","Basal","Her2","Normal","Normal adjacent","Unknown"))

cat("PAM50 distribution in PCA:\n")
print(table(pca_df$PAM50))

## ── Colour palette (validated, CVD-safe) ──────────────────────────────────
pam50_colors <- c(
    "LumA"            = "#2a78d6",   # blue
    "LumB"            = "#eb6834",   # orange
    "Basal"           = "#1baf7a",   # aqua/teal
    "Her2"            = "#e87ba4",   # magenta
    "Normal"          = "#008300",   # green
    "Normal adjacent" = "#4a3aa7",   # violet
    "Unknown"         = "#b0b0b0"    # grey
)
pam50_shapes <- c(
    "LumA" = 16, "LumB" = 17, "Basal" = 15,
    "Her2" = 18, "Normal" = 8, "Normal adjacent" = 3, "Unknown" = 4
)

## ── Plot ───────────────────────────────────────────────────────────────────
cat("Generating plot...\n")

p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = PAM50, shape = PAM50)) +
    geom_point(size = 2.5, alpha = 0.80, stroke = 0.3) +
    scale_color_manual(values = pam50_colors, drop = FALSE) +
    scale_shape_manual(values = pam50_shapes, drop = FALSE) +
    labs(
        title    = "BRCA gene expression PCA",
        subtitle = sprintf("log2(TPM+1), top %d variable genes", nrow(log_mat_top)),
        x        = sprintf("PC1: %s%% variance", pct[1]),
        y        = sprintf("PC2: %s%% variance", pct[2]),
        color    = "PAM50",
        shape    = "PAM50"
    ) +
    coord_fixed() +
    theme_classic(base_size = 13) +
    theme(
        plot.title    = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(color = "#52514e", size = 10),
        legend.title  = element_text(face = "bold"),
        legend.key    = element_blank(),
        panel.grid.major = element_line(color = "#e1e0d9", linewidth = 0.3),
        axis.line     = element_line(color = "#c3c2b7"),
        axis.ticks    = element_line(color = "#c3c2b7")
    )

pdf(o, width = 7, height = 6, useDingbats = FALSE)
suppressWarnings(print(p))
dev.off()

cat(sprintf("Done! => %s\n", o))
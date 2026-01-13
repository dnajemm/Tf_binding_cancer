#!/usr/bin/env Rscript

library(data.table)
library(Rtsne)
library(ggplot2)
library(plotly)
library(htmlwidgets)

# ------------------ arguments ------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript tsne_methylation.R <TCGA-CANCER>\nExample: Rscript tsne_methylation.R TCGA-COAD")
}
cancer <- args[1]

# ------------------ settings ------------------
k_var  <- 3000
n_pcs  <- 30
seed   <- 1

out_dir  <- file.path("./results/tSNE_interactive_plots_methylation", cancer)
out_pdf  <- file.path(out_dir, paste0("tSNE_HM450_", cancer, "_01A_vs_11A.pdf"))
out_tsv  <- file.path(out_dir, paste0("tSNE_HM450_", cancer, "_01A_vs_11A.tsv"))
out_html <- file.path(out_dir, paste0("tsne_", cancer, "_plot.html"))

# ------------------ files (01A + 11A only) ------------------
files <- Sys.glob(paste0("/data/papers/tcga/", cancer, "/*/*_annotated_methylation.bed.gz"))
files <- files[grepl("-01A_|-11A_", basename(files))]

if (length(files) < 3) stop("Not enough 01A/11A files for ", cancer)

group <- ifelse(grepl("-01A_", basename(files)), "Tumor", "Normal")

# ------------------ read methylation ------------------
read_one <- function(f) {
  fread(cmd = paste("zcat", shQuote(f)), header = FALSE, select = c(4, 5))
}

d1 <- read_one(files[1])
probes <- d1[[1]]

X <- sapply(files, function(f) read_one(f)[[2]])
X <- as.matrix(X)
rownames(X) <- probes

# ------------------ select variable CpGs ------------------
sds <- apply(X, 1, sd, na.rm = TRUE)
top <- order(sds, decreasing = TRUE)[1:min(k_var, length(sds))]

Xs <- scale(t(X[top, ]))   # samples x CpGs

# ------------------ PCA + tSNE ------------------
pcs <- min(n_pcs, ncol(Xs) - 1)
if (pcs < 2) stop("Too few samples for PCA/tSNE after filtering for ", cancer)

pca <- prcomp(Xs)
Y <- pca$x[, 1:pcs, drop = FALSE]

set.seed(seed)
perp <- max(5, min(30, floor((nrow(Y) - 1) / 3)))

ts <- Rtsne(Y, perplexity = perp, pca = FALSE, check_duplicates = FALSE)

# ------------------ results table ------------------
plot_dt <- data.table(
  sample = basename(files),
  group  = group,
  tSNE1  = ts$Y[, 1],
  tSNE2  = ts$Y[, 2]
)

# ------------------ save outputs ------------------
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(plot_dt, out_tsv, sep = "\t")

# ------------------ PDF plot ------------------
p_static <- ggplot(plot_dt, aes(tSNE1, tSNE2, color = group)) +
  geom_point(size = 2) +
  theme_bw() +
  labs(
    title = paste0("t-SNE ", cancer, " HM450 (01A vs 11A)"),
    subtitle = paste0("Top ", k_var, " variable CpGs; PCA=", pcs, " dims; perplexity=", perp),
    x = "tSNE1",
    y = "tSNE2",
    color = "Group"
  )

ggsave(out_pdf, p_static, width = 8, height = 6)

# ------------------ interactive plot ------------------
p_interactive <- ggplotly(
  p_static + aes(text = sample),
  tooltip = c("text", "color")
)

saveWidget(p_interactive, out_html, selfcontained = TRUE)

cat("Wrote TSV:  ", out_tsv, "\n")
cat("Wrote PDF:  ", out_pdf, "\n")
cat("Wrote HTML: ", out_html, "\n")

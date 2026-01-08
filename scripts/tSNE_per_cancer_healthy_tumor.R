library(data.table)
library(Rtsne)
library(ggplot2)

cancer  <- "TCGA-COAD"
out_pdf <- "./results/summarys/tSNE_HM450_TCGA-COAD_tumor_vs_normal.pdf"
out_tsv <- "./results/summarys/tSNE_HM450_TCGA-COAD_coords.tsv"

n_var_probes <- 3000
n_pcs <- 30
seed <- 1

# ------------------ collect files ------------------
files <- Sys.glob(file.path("/data/papers/tcga", cancer, "*", "*_annotated_methylation.bed.gz"))
if (length(files) < 3) stop("Not enough annotated methylation files for ", cancer)

get_sample_type <- function(f) {
  m <- regexpr("-[0-9][0-9]A_", basename(f))
  if (m[1] == -1) return(NA_character_)
  substr(basename(f), m[1] + 1, m[1] + 3)  # e.g. 01A, 11A, 06A, 10A
}

sample_type <- vapply(files, get_sample_type, character(1))
keep <- !is.na(sample_type) & (grepl("^[01][0-9]A$", sample_type))  # keep 0*A and 1*A only
files <- files[keep]
sample_type <- sample_type[keep]

if (length(files) < 3) stop("Not enough 0*A/1*A samples after filtering for ", cancer)

# ------------------ define group + shapes ------------------
group <- ifelse(substr(sample_type, 1, 1) == "0", "Tumor", "Normal")

# Shapes: keep tumor sample types distinct (01A/02A/06A/etc), collapse normals into one shape
shape_class <- ifelse(group == "Tumor", sample_type, "Normal")

# ------------------ read methylation (probe+meth only) ------------------
read_probe_meth <- function(f) {
  dt <- fread(cmd = paste("zcat", shQuote(f)), header = FALSE, select = c(4,5))
  setnames(dt, c("probe","meth"))
  dt[, meth := as.numeric(meth)]
  dt
}

message("Reading ", length(files), " files for ", cancer, " ...")

dt1 <- read_probe_meth(files[1])
probes <- dt1$probe
X <- matrix(NA_real_, nrow = length(probes), ncol = length(files))
X[, 1] <- dt1$meth

if (length(files) > 1) {
  for (i in 2:length(files)) {
    dti <- read_probe_meth(files[i])
    if (nrow(dti) != length(probes) || !identical(dti$probe, probes)) {
      stop("Probe order mismatch in: ", files[i], "\n(Your files must have identical probe order.)")
    }
    X[, i] <- dti$meth
    if (i %% 20 == 0) message("... loaded ", i, "/", length(files))
  }
}

# remove probes with any NA
ok <- rowSums(is.na(X)) == 0
X <- X[ok, , drop = FALSE]

# ------------------ select variable probes ------------------
sds <- apply(X, 1, sd)
ord <- order(sds, decreasing = TRUE)
k <- min(n_var_probes, length(ord))
Xv <- X[ord[1:k], , drop = FALSE]

# samples x probes, scale
Xs <- scale(t(Xv))

# PCA pre-reduction
pcs <- min(n_pcs, ncol(Xs) - 1, nrow(Xs) - 1)
if (pcs < 2) stop("Too few samples for PCA/tSNE after filtering.")
pca <- prcomp(Xs, center = FALSE, scale. = FALSE)
Y <- pca$x[, 1:pcs, drop = FALSE]

# tSNE perplexity must be < (n-1)/3
n <- nrow(Y)
perp <- min(30, floor((n - 1) / 3))
perp <- max(perp, 5)

set.seed(seed)
ts <- Rtsne(Y, perplexity = perp, check_duplicates = FALSE, pca = FALSE, max_iter = 1000)

meta <- data.table(
  sample = basename(files),
  sample_type = sample_type,
  group = group,
  shape_class = shape_class,
  file = files
)

plot_dt <- cbind(meta, data.table(tSNE1 = ts$Y[,1], tSNE2 = ts$Y[,2]))
dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
fwrite(plot_dt, out_tsv, sep = "\t")

# ------------------ plot ------------------
p <- ggplot(plot_dt, aes(tSNE1, tSNE2, color = group, shape = shape_class)) +
  geom_point(alpha = 0.75, size = 1.8) +
  theme_bw() +
  labs(
    title = "t-SNE of TCGA-BRCA HM450 methylation (tumor + normals)",
    subtitle = paste0("Top ", k, " variable CpGs; PCA=", pcs, " dims; perplexity=", perp),
    x = "tSNE1", y = "tSNE2", color = "Group", shape = "Sample type"
  ) +
  theme(legend.position = "right")

ggsave(out_pdf, p, width = 10, height = 7)

cat("Wrote:", out_pdf, "\n")
cat("Wrote:", out_tsv, "\n")

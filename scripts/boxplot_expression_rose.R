# plot boxplot of each selected TF separately
# color = cancer type
# shape = Tumor vs Healthy
# TPM is read directly from RNA-seq files, column 7
# no restriction to 2D cohort

Rscript - <<'EOF'

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(grid)
})

# ============================================================
# USER PARAMETERS
# ============================================================

tf_list <- c("TERT")

rnaseq_dir <- "./expression/"

rnaseq_pattern <- "augmented_star_gene_counts_filtered\\.tsv$"

gene_col <- 2
tpm_col  <- 7

out_dir <- "./results/"

# ============================================================
# Read mapping file
# ============================================================

cat("[INFO] Reading mapping file...\n")

map <- read_tsv(
  "./all_tcga_samples_with_cancer.tsv",
  col_names = c("sample_id", "cancer"),
  show_col_types = FALSE
)

map <- map %>%
  mutate(
    case_id = str_replace(sample_id, "-[0-9]{2}[A-Za-z]$", ""),
    cancer = str_replace(cancer, "^TCGA-", "")
  )

# ============================================================
# Read cancer color file
# ============================================================

cat("[INFO] Reading cancer color file...\n")

color_df <- read_tsv(
  "./results/multi_omics/cancer_color_order_with_defined_colours.tsv",
  col_names = c("cancer_raw", "color"),
  show_col_types = FALSE
) %>%
  mutate(
    cancer = str_replace(cancer_raw, "^TCGA-", ""),
    color = str_trim(color)
  )

cancer_palette <- setNames(color_df$color, color_df$cancer)

# ============================================================
# List RNA-seq files
# ============================================================

cat("[INFO] Listing RNA-seq files...\n")

files <- list.files(
  rnaseq_dir,
  pattern = rnaseq_pattern,
  full.names = TRUE
)

cat("[INFO] Found ", length(files), " RNA-seq files\n", sep = "")

if (length(files) == 0) {
  stop("No RNA-seq files found in: ", rnaseq_dir)
}

# ============================================================
# Helper function to extract sample ID from file name
# ============================================================

extract_sample_id <- function(filename) {
  str_extract(
    basename(filename),
    "TCGA-[A-Za-z0-9]{2}-[A-Za-z0-9]{4}-[0-9]{2}[A-Za-z]"
  )
}

# ============================================================
# Read RNA-seq files directly
# ============================================================

cat("[INFO] Reading RNA-seq expression files directly...\n")

df <- bind_rows(lapply(files, function(f) {

  cat("[READ] ", basename(f), "\n", sep = "")

  sample_id <- extract_sample_id(f)

  if (is.na(sample_id)) {
    warning("Could not extract TCGA sample ID from file name: ", basename(f))
    return(NULL)
  }

  x <- read_tsv(
    f,
    col_names = FALSE,
    show_col_types = FALSE
  )

  if (ncol(x) < max(gene_col, tpm_col)) {
    warning(
      "Skipping file with fewer columns than expected: ",
      basename(f),
      " | ncol = ",
      ncol(x)
    )
    return(NULL)
  }

  tibble(
    sample_id = sample_id,
    gene_name = as.character(x[[gene_col]]),
    tpm_unstranded = suppressWarnings(as.numeric(x[[tpm_col]]))
  ) %>%
    filter(gene_name %in% tf_list)
}))

if (nrow(df) == 0) {
  stop("No selected TFs were found in the RNA-seq files.")
}

# ============================================================
# Prepare table and join cancer annotation
# ============================================================

cat("[INFO] Preparing expression table and joining cancer annotation...\n")

df <- df %>%
  mutate(
    case_id = str_replace(sample_id, "-[0-9]{2}[A-Za-z]$", ""),
    code = str_match(sample_id, "-([0-9]{2})[A-Za-z]$")[, 2],
    sample_type = case_when(
      code %in% c("01", "02", "03", "05", "06", "07", "08", "09") ~ "Tumor",
      code %in% c("10", "11", "12", "13", "14") ~ "Healthy",
      TRUE ~ "Other"
    )
  ) %>%
  left_join(
    map %>% select(case_id, cancer) %>% distinct(),
    by = "case_id"
  ) %>%
  filter(!is.na(cancer)) %>%
  filter(!is.na(tpm_unstranded)) %>%
  filter(sample_type %in% c("Healthy", "Tumor"))

cat(
  "[INFO] Rows after keeping all available annotated RNA-seq samples: ",
  nrow(df),
  "\n",
  sep = ""
)

if (nrow(df) == 0) {
  stop("No annotated RNA-seq samples left after filtering.")
}

cat("[INFO] Sample code breakdown:\n")
print(table(df$code, df$sample_type, useNA = "ifany"))

# ============================================================
# Print summary
# ============================================================

cat(
  "[INFO] Tumor samples: ",
  length(unique(df$sample_id[df$sample_type == "Tumor"])),
  "\n",
  sep = ""
)

cat(
  "[INFO] Healthy samples: ",
  length(unique(df$sample_id[df$sample_type == "Healthy"])),
  "\n",
  sep = ""
)

cat(
  "[INFO] Cancers: ",
  length(unique(df$cancer)),
  "\n",
  sep = ""
)

cat("[INFO] TFs found after filtering:\n")
print(table(df$gene_name))

# ============================================================
# Prepare colors
# ============================================================

present_cancers <- sort(unique(df$cancer))
missing_colors <- setdiff(present_cancers, names(cancer_palette))

if (length(missing_colors) > 0) {
  warning(
    "Missing colors for cancers: ",
    paste(missing_colors, collapse = ", "),
    ". They will appear in grey."
  )

  cancer_palette <- c(
    cancer_palette,
    setNames(rep("grey70", length(missing_colors)), missing_colors)
  )
}

cancer_palette <- cancer_palette[present_cancers]

# ============================================================
# Plot function
# ============================================================

make_plot <- function(d, gene) {

  if (nrow(d) == 0) {
    stop("No data to plot for gene: ", gene)
  }

  ymax <- max(d$tpm_unstranded, na.rm = TRUE)

  ggplot(d, aes(x = cancer, y = tpm_unstranded)) +

    geom_boxplot(
      aes(color = cancer),
      outlier.shape = NA,
      linewidth = 0.5,
      show.legend = FALSE
    ) +

    geom_jitter(
      data = d %>% filter(sample_type == "Tumor"),
      aes(color = cancer, shape = sample_type),
      width = 0.2,
      alpha = 0.65,
      size = 2
    ) +

    geom_jitter(
      data = d %>% filter(sample_type == "Healthy"),
      aes(fill = cancer, shape = sample_type),
      width = 0.2,
      alpha = 0.65,
      size = 2,
      color = "black",
      stroke = 0.6
    ) +

    geom_hline(
      yintercept = 1,
      color = "blue",
      linewidth = 0.8
    ) +

    scale_color_manual(
      values = cancer_palette,
      drop = FALSE
    ) +

    scale_fill_manual(
      values = cancer_palette,
      drop = FALSE
    ) +

    scale_shape_manual(
      values = c(
        Healthy = 24,
        Tumor = 16
      )
    ) +

    guides(
      fill = "none",
      color = guide_legend(
        override.aes = list(size = 4, linewidth = 1)
      ),
      shape = guide_legend(
        override.aes = list(size = 4)
      )
    ) +

    coord_cartesian(ylim = c(0, ymax)) +

    theme_bw() +

    theme(
      axis.text.x = element_text(
        angle = 60,
        hjust = 1
      ),

      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      ),

      plot.subtitle = element_text(
        hjust = 0.5
      ),

      legend.position = "right",

      legend.title = element_text(
        size = 14,
        face = "bold"
      ),

      legend.text = element_text(
        size = 12
      ),

      legend.key.size = unit(0.7, "cm")
    ) +
    labs(
      title = paste0(gene, " TPM across all available RNA-seq samples"),
      subtitle = "Healthy vs Tumor | all available RNA-seq samples",
      x = "Cancer type",
      y = "TPM",
      color = "Cancer type",
      shape = "Sample type"
    )
}

# ============================================================
# Build and save one PDF per TF
# ============================================================

cat("[INFO] Building plots for selected TFs...\n")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (tf in tf_list) {

  d_tf <- df %>%
    filter(gene_name == tf)

  if (nrow(d_tf) == 0) {
    warning("No data found for TF after filtering: ", tf)
    next
  }

  p_tf <- make_plot(
    d = d_tf,
    gene = tf
  )

  out_tf <- file.path(
    out_dir,
    paste0(
      "boxplot_TPM_",
      tf,
      "_by_cancer_Healthy_vs_Tumor_all_available_RNAseq.pdf"
    )
  )

  cat("[INFO] Writing PDF for ", tf, "...\n", sep = "")

  pdf(out_tf, width = 14, height = 6, onefile = TRUE)
  print(p_tf)
  dev.off()

  cat("[DONE] Wrote ", out_tf, "\n", sep = "")
}

cat("[DONE] All selected TF plots completed.\n")

EOF
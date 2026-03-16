# Plots for expected results poster M2: 

Rscript -e ' 
# Expected Results: schematic methylation–expression relationship
library(ggplot2)

# Fake conceptual data for a clean, illustrative plot
set.seed(42)

normal <- data.frame(
  methylation = seq(0.05, 0.45, length.out = 40),
  expression  = 20 - 2 * seq(0.05, 0.45, length.out = 40) + rnorm(40, 0, 0.3),
  group = "Normal")

tumor <- data.frame(
  methylation = seq(0.40, 0.95, length.out = 40),
  expression  = 17 - 5 * (seq(0.40, 0.95, length.out = 40) - 0.40) + rnorm(40, 0, 0.3),
  group = "Tumor")

df <- rbind(normal, tumor)

# Plot
ggplot(df, aes(x = methylation, y = expression, color = group)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, size = 1.2) +
  scale_color_manual(values = c("Normal" = "#1B9E77", "Tumor" = "#D95F02")) +
  labs(
    x = "Methylation",
    y = "Gene Expression",
    title = "Expected Relationship: Hypermethylation → Reduced Expression",
    subtitle = "NRF1/BANP target genes are expected to show decreased expression\nwhen promoter CpGs become hypermethylated"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 11),
    legend.position = "top"
  )
ggsave("./results/summarys/Expected_methylation_expression_relationship_NRF1_BANP_targets.pdf", width = 8, height = 6, dpi = 300, bg = "white")
'

#expected Variant effect on TF binding plot
Rscript -e '
library(ggplot2)
set.seed(123)
# 1. Simulate expression
n <- 80

# Motif with variant: lower expression
expr_variant <- rlnorm(n, meanlog = 0.6, sdlog = 0.6)

# WT motif: higher expression
expr_wt      <- rlnorm(n, meanlog = 2.5, sdlog = 0.4)

df <- rbind(
  data.frame(group = "Motif with variant", expression = expr_variant),
  data.frame(group = "WT motif",           expression = expr_wt))

# 2. Plot
p <- ggplot(df, aes(x = group, y = expression)) +
  # violins (pastel)
  geom_violin(aes(fill = group), trim = FALSE, alpha = 0.6, width = 0.9, colour = NA) +
  # boxplots (darker, on top)
  geom_boxplot(aes(fill = group),
               width = 0.18,
               outlier.size = 0.7,
               colour = "black",
               alpha = 0.9) +
  scale_fill_manual(values = c(
    "Motif with variant" = "#E67E22",  # peach
    "WT motif"           = "#2ECC71"   # mint
  )) +
  scale_y_log10(
    breaks = c(1, 3, 10, 30),
    labels = c("1", "3", "10", "30")
  ) +
  labs(
    title = "Effect of Variants in NRF1/BANP Motifs on\nGene Expression for Each Cancer Type",
    x     = "",
    y     = "Gene expression"
  ) +
  theme_bw(base_size = 16) +
  theme(
    plot.title   = element_text(hjust = 0.5, face = "bold", size = 20),
    axis.title.y = element_text(size = 18),
    axis.text.x  = element_text(size = 18),
    axis.text.y  = element_text(size = 14),
    legend.position = "none"
  )

# 3. Save PDF
ggsave(
  filename = "./results/summarys/Effect_variants_NRF1_BANP_gene_expression.pdf",
  plot     = p,
  width    = 8,
  height   = 5
)
'

# expected gene network plot :
Rscript -e '
library(igraph)
# DEFINE NODES (MULTIPLE TFs + MANY GENES)
tf_names   <- paste0("TF_", 1:4)
gene_names <- paste0("GENE_", LETTERS[1:12])

nodes <- data.frame(
  name   = c(tf_names, gene_names),
  type   = c(rep("TF", length(tf_names)), rep("Gene", length(gene_names))),
  status = c(
    rep("TF", length(tf_names)),  # TFs
    "Altered",  "Altered",  "Unchanged", "Unchanged",
    "Altered",  "Unchanged", "Altered",   "Unchanged",
    "Altered",  "Unchanged", "Unchanged", "Altered"
  ),
  stringsAsFactors = FALSE)

# DEFINE EDGES (MULTIPLE TF -> MANY GENES, SOME SHARED TARGETS)
edge_list <- rbind(
  data.frame(from = "TF_1", to = c("GENE_A", "GENE_B", "GENE_C", "GENE_D", "GENE_E")),
  data.frame(from = "TF_2", to = c("GENE_C", "GENE_F", "GENE_G", "GENE_H")),
  data.frame(from = "TF_3", to = c("GENE_A", "GENE_D", "GENE_I", "GENE_J", "GENE_K")),
  data.frame(from = "TF_4", to = c("GENE_E", "GENE_G", "GENE_H", "GENE_L")))

# add a few gene–gene edges
gene_gene_edges <- data.frame(
  from = c("GENE_B", "GENE_F", "GENE_I"),
  to   = c("GENE_C", "GENE_G", "GENE_J"))

edges <- rbind(edge_list, gene_gene_edges)

# CREATE GRAPH
g <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)

# COLOR / SHAPE / SIZE
vertex_colors <- ifelse(
  V(g)$type == "TF", "#1F77B4",
  ifelse(V(g)$status == "Altered", "#D62728", "lightgrey"))

vertex_shapes <- ifelse(V(g)$type == "TF", "square", "circle")
vertex_sizes  <- ifelse(V(g)$type == "TF", 28, 23)

# LAYOUT: TFs INNER CIRCLE, GENES OUTER CIRCLE
n_tf   <- length(tf_names)
n_gene <- length(gene_names)

theta_tf   <- seq(0, 2*pi, length.out = n_tf + 1)[-(n_tf + 1)]
theta_gene <- seq(0, 2*pi, length.out = n_gene + 1)[-(n_gene + 1)]

layout_tf <- cbind(0.6 * cos(theta_tf), 0.6 * sin(theta_tf))
layout_gene <- cbind(1.4 * cos(theta_gene), 1.4 * sin(theta_gene))

layout_final <- rbind(layout_tf, layout_gene)

# PLOT TO PDF
pdf("./results/summarys/TF_multicenter_gene_network_expected.pdf", width = 6, height = 6)

plot(
  g,
  layout             = layout_final,
  vertex.color       = vertex_colors,
  vertex.shape       = vertex_shapes,
  vertex.size        = vertex_sizes,
  vertex.label       = V(g)$name,
  vertex.label.cex   = 0.8,
  vertex.label.color = "black",
  edge.width         = 2,
  edge.color         = "grey50",
  main = "Expected Gene Network: Multi-TF model\nCancer-altered targets highlighted in red"
)

dev.off()
'
# expected methylation gene expression plot 
Rscript -e 'library(ggplot2);

set.seed(123);

# 1. Simulate delta methylation (%)
n <- 120;   # fewer points
delta_meth <- runif(n, -100, 100);

# 2. Simulate log2FC with coupling
log2fc <- -0.03 * delta_meth + rnorm(n, sd = 0.5);

df <- data.frame(
  delta_meth = delta_meth,
  log2FC     = log2fc
);

# 3. Fit linear model
fit <- lm(log2FC ~ delta_meth, data = df);

line_df <- data.frame(
  intercept = coef(fit)[1],
  slope     = coef(fit)[2]
);

# 4. Plot (large colored points)
p <- ggplot(df, aes(x = delta_meth, y = log2FC, color = delta_meth)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_point(alpha = 0.8, size = 2) +       
  scale_color_gradient2(
    low = "#1f77b4", mid = "white", high = "#d62728",
    midpoint = 0,
    name = "Methylation Change"
  ) +
  geom_abline(
    data = line_df,
    aes(intercept = intercept, slope = slope),
    colour = "black",
    linewidth = 1
  ) +
  labs(
    title    = "Expected Cancer-Specific Impact of Promoter Motif Alterations",
    subtitle = "Promoter hypermethylation at TF motifs reduces gene expression",
    x        = "Change in promoter methylation (Cancer - Healthy, %)",
    y        = "Change in expression (log2 fold change)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 11));

# 5. Save PDF
ggsave(
  filename = "./results/summarys/Expected_methylation_expression_coupling.pdf",
  plot     = p,
  width    = 8,
  height   = 5
);'



############################################
# IMCbio presentation expected results
############################################
# NRF1 and BANP gene expression for imcbio presentation
Rscript -e '
library(data.table)
library(ggplot2)

in_file  <- "./expression/all_samples_NRF1_BANP_expression.tsv"
col_file <- "./results/multi_omics/cancer_color_order_with_defined_colours.tsv"
out_file <- "./results/expression/NRF1_BANP_boxplot_expression_sameScale_TPMge1_cap50.pdf"

dt <- fread(in_file)

# --- Load cancer colors (2 columns: cancer, color) robustly ---
cc <- fread(col_file)
if (ncol(cc) < 2) stop("Color file must have at least 2 columns: cancer and color")

# Keep only first two cols if extra exist
cc <- cc[, .(cancer = get(names(cc)[1]), color = get(names(cc)[2]))]

# Normalize cancer naming: remove possible "TCGA-" prefix
cc[, cancer := sub("^TCGA-","", cancer)]
dt[, cancer := sub("^TCGA-","", cancer)]

# --- Filter expression table ---
dt <- dt[gene %in% c("NRF1","BANP")]
dt <- dt[tpm_unstranded >= 1]

# Intersect cancers present in both
common <- intersect(unique(dt$cancer), unique(cc$cancer))
dt <- dt[cancer %in% common]
cc <- cc[cancer %in% common]

# Safety check
if (nrow(dt) == 0) {
  stop("After filtering, dt is empty. Check that cancer names match between expression file and color file.")
}

# Named vector: cancer -> color
cols <- setNames(cc$color, cc$cancer)

# Keep a consistent cancer order (from color file)
dt[, cancer := factor(cancer, levels = cc$cancer)]

p <- ggplot(dt, aes(x=cancer, y=tpm_unstranded, fill=cancer)) +
  geom_boxplot(outlier.shape=NA, width=0.7) +
  geom_hline(yintercept=1, color="red", linewidth=1) +
  facet_wrap(~gene, scales="fixed") +
  scale_fill_manual(values=cols, drop=FALSE) +
  coord_cartesian(ylim=c(1,50)) +
  theme_bw(base_size=14) +
  theme(
    axis.text.x = element_text(angle=60, hjust=1),
    legend.position="none",
    strip.text = element_text(size=14, face="bold")
  ) +
  labs(
    x="Cancer type",
    y="TPM expression (TPM ≥ 1)",
    title="NRF1 and BANP Expression Across Cancers"
  )

ggsave(out_file, p, width=14, height=6)
cat("Saved:", out_file, "\n")
'

# expected healthy vs cancer boxplot for imcbio presentation
Rscript -e '
library(data.table)
library(ggplot2)

set.seed(1)

# ---------------------------
# Synthetic "expected" data
# ---------------------------
n_healthy <- 60
n_cancer  <- 120

effect_size <- 1.2   # cancer shift (negative = downregulation)
sd_common   <- 0.8

expr_healthy <- rnorm(n_healthy, mean=4.0, sd=sd_common)
expr_cancer  <- rnorm(n_cancer,  mean=4.0 + effect_size, sd=sd_common)

dt <- data.table(
  group = factor(c(rep("Healthy", n_healthy), rep("Cancer", n_cancer)),
                 levels=c("Cancer","Healthy")),
  expr  = c(expr_healthy, expr_cancer)
)

# Colors you asked for
cols <- c("Healthy" = "#2EAD64",  # green
          "Cancer"  = "#E53935")  # red

out_file <- "./expected_healthy_vs_cancer_boxplot_ggplot_style.pdf"
dir.create(dirname(out_file), recursive=TRUE, showWarnings=FALSE)

p <- ggplot(dt, aes(x=group, y=expr, fill=group)) +
  geom_boxplot(outlier.shape=NA, width=0.65, color="black") +
  scale_fill_manual(values=cols) +
  coord_flip() +
  theme_bw(base_size=14) +
  theme(
    legend.position="none",
    strip.text = element_text(size=14, face="bold")
  ) +
  labs(
    x="",
    y="Synthetic gene expression (log2(TPM+1)-like)",
    title="Expected expression difference (synthetic)"
  )

ggsave(out_file, p, width=5.6, height=4.2)

cat("Saved:", out_file, "\n")
'

# cancer type rewiring score plot 
Rscript -e '
library(data.table)
library(ggplot2)

col_file <- "./results/multi_omics/cancer_color_order_with_defined_colours.tsv"

cc <- fread(col_file)
if (ncol(cc) < 2) stop("Color file must have at least 2 columns: cancer and color")

cc <- cc[, .(cancer = get(names(cc)[1]), color = get(names(cc)[2]))]
cc[, cancer := sub("^TCGA-","", cancer)]

cols <- setNames(cc$color, cc$cancer)

df <- data.frame(
  cancer = c("BRCA","LUAD","COAD","KIRC","ACC","LAML"),
  rewiring_score = c(0.92, 0.68, 0.48, 0.31, 0.18, 0.12)
)

df$cancer <- sub("^TCGA-","", df$cancer)
df <- df[df$cancer %in% names(cols), ]
if (nrow(df) == 0) stop("After matching with the color file, df is empty. Check cancer names.")

top_n <- 6
df <- df[order(df$rewiring_score, decreasing = TRUE), ]
df <- head(df, top_n)

# Order cancers left-to-right by decreasing score
df$cancer <- factor(df$cancer, levels = df$cancer)

p <- ggplot(df, aes(x = cancer, y = rewiring_score, fill = cancer)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = cols, drop = FALSE) +
  scale_y_continuous(
    limits = c(0, max(df$rewiring_score) * 1.1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "Cancer-specific promoter rewiring",
    x = NULL,
    y = "Rewiring score"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(face = "bold", angle = 45, hjust = 1),
    legend.position = "none"
  )

ggsave("./right_panel_rewiring.png", p, width = 7.5, height = 4.2, dpi = 300, bg = "white")
cat("Saved: ./right_panel_rewiring.png\n")
'
#############################################################################################s


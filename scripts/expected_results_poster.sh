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

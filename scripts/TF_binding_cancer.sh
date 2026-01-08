#!/bin/bash
# To connect to the server dorsal : ssh najemd@serv-bardet-dorsal.igbmc.u-strasbg.fr
cd TF_binding_cancer

# 1) Filter motifs with p-value threshold and sort them.
mkdir -p ./motifs

zcat /data/genome/motifs/jaspar_2024/fimo/vertebrata/hg38/NRF1.MA0506.3.bed.gz | awk '($5<=1/4^6)' > ./motifs/NRF1_filtered_6mer.MA0506.3.bed
zcat /data/genome/motifs/jaspar_2024/fimo/vertebrata/hg38/ZBTB33.MA0527.2.bed.gz | awk '($5<=1/4^6)' > ./motifs/ZBTB33_filtered_6mer.MA0527.2.bed

# barplot of the number of motifs before and after filtering 
Rscript -e '
library(ggplot2)
motifs <- c("NRF1", "BANP")
before_counts <- c( nrow(read.table(gzfile("/data/genome/motifs/jaspar_2024/fimo/vertebrata/hg38/NRF1.MA0506.3.bed.gz"))),
                    nrow(read.table(gzfile("/data/genome/motifs/jaspar_2024/fimo/vertebrata/hg38/ZBTB33.MA0527.2.bed.gz"))) )
after_counts  <- c( nrow(read.table("./motifs/NRF1_filtered_6mer.MA0506.3.bed")),
                    nrow(read.table("./motifs/ZBTB33_filtered_6mer.MA0527.2.bed")) )

statuses <- c("Before 6mer Filtering", "After 6mer Filtering")

df <- data.frame(
  Motif  = rep(motifs, times = 2),                 # NRF1, BANP, NRF1, BANP
  Status = rep(statuses, each = length(motifs)),   # Before, Before, After, After
  Count  = c(before_counts, after_counts)          # b_NRF1, b_BANP, a_NRF1, a_BANP
)

df$Status <- factor(df$Status,levels = c("Before 6mer Filtering", "After 6mer Filtering"))

# Plot
ggplot(df, aes(x = Motif, y = Count, fill = Status)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = Count),
          position = position_dodge(width = 0.9),
          vjust = -0.3,
          size = 3) +
  theme_minimal() +
  labs(title = "Number of motifs before and after 6mer filtering",
       subtitle = "Filtering threshold: p-value <= 1/4^6",
       y = "Count",
       x = "Motif") +
  theme(plot.title = element_text(hjust = 0.7 , size = 12),
      plot.subtitle = element_text(hjust = 0.7 , size = 9))
# Save plot
ggsave("./results/summarys/motif_counts_before_after_6mer_filtering.pdf", width = 6, height = 4, dpi = 300, bg = "white")
'
# barplot to see only the number of motifs after filtering
Rscript -e '
library(ggplot2)
motifs <- c("NRF1", "BANP")
after_counts  <- c( nrow(read.table("./motifs/NRF1_filtered_6mer.MA0506.3.bed")),
                    nrow(read.table("./motifs/ZBTB33_filtered_6mer.MA0527.2.bed")) )
df <- data.frame(
  Motif  = motifs,
  Count  = after_counts
)
# Plot
ggplot(df, aes(x = Motif, y = Count, fill = Motif)) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_text(aes(label = Count),
          vjust = -0.3,
          size = 4) +
  theme_minimal() +
  labs(title = "Number of motifs after 6mer filtering",
       subtitle = "Filtering threshold: p-value <= 1/4^6",
       y = "Count",
       x = "Motif") +
  theme(plot.title = element_text(hjust = 0.7 , size = 12),
      plot.subtitle = element_text(hjust = 0.7 , size = 9),
      legend.position = "none")
# Save plot
ggsave("./results/summarys/motif_counts_after_6mer_filtering.pdf", width = 5, height = 4, dpi = 300, bg = "white")
'


# Distance of motifs to TSS : generate a file per sample with the distances of each motif to the nearest TSS
mkdir -p ./motifs/dist_to_tss/

for f in ./motifs/*_filtered_6mer.*.bed; do
    base=$(basename "$f")                     
    sample=${base%_filtered_6mer.*.bed} 
    echo "Processing $sample ..."

    closestBed -t first -d -a "$f" -b /data/genome/annotations/hg38_tss.bed | awk '{print $NF}' > "./motifs/dist_to_tss/${sample}_dist_tss.txt"
done

# Plot histogram of distances to TSS for all samples combined

Rscript -e '
library(ggplot2);

motif_files <- list.files("motifs", pattern = "_filtered_6mer.*.bed$", full.names = TRUE)
samples <- sub("_filtered_6mer.*.bed$", "", basename(motif_files))

pdf("./results/summarys/dist_tss_motifs_hist.pdf") 
par(bg = "white")

for (sample in samples) {
  dist_file <- paste0("./motifs/dist_to_tss/", sample, "_dist_tss.txt")
  distances <- read.table(dist_file)[,1]

  # compute proximal/distal percentages
  proximal <- sum(distances < 2000)
  distal   <- sum(distances >= 2000)
  total    <- length(distances)
  pct_prox <- round((proximal / total) * 100, 1)
  pct_dist <- round((distal   / total) * 100, 1)

# rename sample for display to rename ZBTB33 to BANP
  display_name <- sample
  if (sample == "ZBTB33") {
    display_name <- "BANP"} 

  # plot
  hist(log10(distances + 1),
       xlim = c(0, 8),
       xlab = "Distance to TSS (log10)",
       main = paste0(display_name, "  |  Proximal = ", pct_prox, "%   |   Distal = ", pct_dist, "%"),
       breaks = 50,
       col = "steelblue",
       border = "black")
  # red vertical line at log10(2000) --> promotor 
  abline(v = log10(2000), col = "red", lwd = 2)}
dev.off()
' 


mkdir -p ./methylation
# 2) Download the HM450 probes bed file :
#wget -O ./methylation/HM450.hg38.manifest.tsv.gz \
#https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/HM450/HM450.hg38.manifest.tsv.gz

#3) Setting up R environment with required packages
mkdir -p ./conda_envs/conda_pkgs
export CONDA_PKGS_DIRS=/data/najemd/TF_binding_cancer/conda_envs/conda_pkgs
#conda create --prefix ./conda_envs/TF_binding_cancer_env r-base -c conda-forge -y
#conda env export --prefix ./conda_envs/TF_binding_cancer_env > ./conda_envs/TF_binding_cancer_env.yml
conda activate ./conda_envs/TF_binding_cancer_env
#conda install -c conda-forge r-venndiagram r-eulerr r-ggplot2 r-dplyr r-tidyr r-hexbin r-mass r-kernsmooth r-venndir --yes

# 4) Overlap 6mer motifs with annotated methylation data to find intersecting regions and 
# see if the motifs fall within methylated probes regions.

mkdir -p ./motifs/overlaps/intersected_motifs_HM450

# transform the HM450 manifest annotated methylation data to bed format. chr starrt end probe_id
zcat ./methylation/HM450.hg38.manifest.tsv.gz | awk -vOFS="\t" ' NR>1 && $1 != "NA" && $2 != "NA" && $3 != "NA" {print $1, $2, $3, $9}' | sort -k1,1 -k2,2n  > ./methylation/annotated_methylation_data_probes.bed

# intersect NRF1 6mer motifs with annotated methylation data
bedtools intersect -u -a ./motifs/NRF1_filtered_6mer.MA0506.3.bed -b ./methylation/annotated_methylation_data_probes.bed > ./motifs/overlaps/intersected_motifs_HM450/NRF1_intersected_methylation.bed
bedtools intersect -u -a ./motifs/ZBTB33_filtered_6mer.MA0527.2.bed -b ./methylation/annotated_methylation_data_probes.bed > ./motifs/overlaps/intersected_motifs_HM450/ZBTB33_intersected_methylation.bed

# 5) Summarize the intersected results : 1   NRF1 30234 ZBTB33 6119
mkdir -p ./results/summarys

# Count lines in bash and save to variables
n_NRF1=$(cat ./motifs/overlaps/intersected_motifs_HM450/NRF1_intersected_methylation.bed | wc -l)
n_ZBTB33=$(cat ./motifs/overlaps/intersected_motifs_HM450/ZBTB33_intersected_methylation.bed | wc -l)
# Create a summary table using R
Rscript -e '
table <- data.frame(
  Motif = c("NRF1","ZBTB33"),
  Overlapping_regions = c('"$n_NRF1"','"$n_ZBTB33"'))
print(table)
'
# Create a Venn diagram to visualize the overlaps between motifs and methylation for NRF1 

# Count total number of 6mer motifs
n_motifs_NRF1=$(wc -l < ./motifs/NRF1_filtered_6mer.MA0506.3.bed)

# Count total number of peaks 
n_probes=$(wc -l < ./methylation/annotated_methylation_data_probes.bed)

# Count number of 6mer NRF1 motifs that fall within peaks
n_overlap_NRF1=$(wc -l < ./motifs/overlaps/intersected_motifs_HM450/NRF1_intersected_methylation.bed)

Rscript -e '
library(VennDiagram)
library(grid)

motifs  <- '"$n_motifs_NRF1"'
probes   <- '"$n_probes"'
overlap <- '"$n_overlap_NRF1"'

venn.obj <- draw.pairwise.venn(
  area1      = motifs,
  area2      = probes,
  cross.area = overlap,
  category   = c("NRF1 6mer motifs", "Methylation probes"),
  fill       = c("salmon", "skyblue"),
  alpha      = 0.7,
  cex        = 2,
  cat.cex    = 2.2,
  fontfamily = "Helvetica",
  cat.pos    = c(-15, 20),
  cat.dist   = c(0.03, 0.02),
  cat.fontfamily = "Helvetica",
  print.mode = c("raw","percent"),
  sigdig     = 2,
  ind        = FALSE)

png("./results/summarys/NRF1_motifs_vs_methylation_venn.png", width = 1300, height = 1300)
grid.newpage()
# Add title
grid.text(
  "Overlap between NRF1 6mer motifs and Methylation probes",
  x = 0.5, y = 0.96,
  gp = gpar(fontsize = 36,fontfamily = "Helvetica"))

pushViewport(viewport(y = 0.45))
grid.draw(venn.obj)
dev.off()
'

# Create a Venn diagram to visualize the overlaps between motifs and methylation for BANP 

# Count total number of 6mer motifs
n_motifs_BANP=$(wc -l < ./motifs/ZBTB33_filtered_6mer.MA0527.2.bed)

# Count total number of probes 
n_probes=$(wc -l < ./methylation/annotated_methylation_data_probes.bed)

# Count number of 6mer BANP motifs that fall within methylation probes
n_overlap_BANP=$(wc -l < ./motifs/overlaps/intersected_motifs_HM450/ZBTB33_intersected_methylation.bed)

Rscript -e '
library(VennDiagram)
library(grid)

motifs  <- '"$n_motifs_BANP"'
probes   <- '"$n_probes"'
overlap <- '"$n_overlap_BANP"'

# Generate venn WITHOUT category labels
venn.obj <- draw.pairwise.venn(
  area1      = motifs,
  area2      = probes,
  cross.area = overlap,
  category   = c("", ""),               
  fill       = c("salmon", "skyblue"),
  alpha      = 0.7,
  cex        = 1.9,       # size of the numbers inside the circles
  cat.cex    = 2.2,
  fontfamily = "Helvetica",   # font family for numbers
  print.mode = c("raw","percent"),
  sigdig     = 2,
  ind        = FALSE)

png("./results/summarys/BANP_motifs_vs_methylation_venn.png", width = 1300, height = 1300)
grid.newpage()

#Title
grid.text("Overlap between BANP 6mer motifs and methylation probes",
          x = 0.5, y = 0.95,
          gp = gpar(fontsize = 36,fontfamily = "Helvetica"))
# place the category labels above the circles
grid.text("Methylation probes",
          x = 0.37, y = 0.73,           
          gp = gpar(fontsize = 28,fontfamily = "Helvetica"))
grid.text("BANP 6mer motifs",
          x = 0.75, y = 0.65,         
          gp = gpar(fontsize = 28,fontfamily = "Helvetica"))
# Venn diagram
pushViewport(viewport(y = 0.45))
grid.draw(venn.obj)

dev.off()
'

# 6) link the peak files from the previous analysis of /data/hichamif/pred_tf_cancer/peak/ to this folder :
#mkdir -p ./peaks
#for file in /data/hichamif/pred_tf_cancer/peak/ATAC_TCGA*_peaks_macs.bed; do
#    ln -s "$file" ./peaks/
#done

#combine all peak files into one and sort them
#cat ./peaks/ATAC_TCGA*_peaks_macs.bed | sort -k1,1 -k2,2n | bedtools merge -i stdin > ./peaks/merged_peaks.bed

# 7) Overlap the 6mer motifs with the peaks to see which motifs fall within the peaks regions.
mkdir -p ./motifs/overlaps/motif_peak_overlaps

# Loop through each 6mer motif file and intersect with all peak files from peaks
for file in ./motifs/*_filtered_6mer*.bed; do
   motif_name=$(basename "$file" .bed)
   bedtools intersect -u -a "$file" -b ./peaks/merged_peaks.bed > ./motifs/overlaps/motif_peak_overlaps/"${motif_name}_peak_overlaps.bed"
done

# Summarize the intersected results : NRF1 384997 ZBTB33 69795

# Count lines in bash and save to variables
n_NRF1=$(cat ./motifs/overlaps/motif_peak_overlaps/NRF1_filtered_6mer.MA0506.3_peak_overlaps.bed | wc -l)
n_ZBTB33=$(cat ./motifs/overlaps/motif_peak_overlaps/ZBTB33_filtered_6mer.MA0527.2_peak_overlaps.bed | wc -l)
# Create a summary table using R
Rscript -e '
table <- data.frame(
  Motif = c("NRF1","ZBTB33"),
  Overlapping_regions = c('"$n_NRF1"','"$n_ZBTB33"'))
print(table)
'
# Create a Venn diagram to visualize the overlaps between motifs and peaks for NRF1 

# Count total number of 6mer motifs
n_motifs_NRF1=$(wc -l < ./motifs/NRF1_filtered_6mer.MA0506.3.bed)

# Count total number of peaks 
n_peaks=$(wc -l < ./peaks/merged_peaks.bed)

# Count number of 6mer NRF1 motifs that fall within peaks
n_overlap_NRF1=$(wc -l < ./motifs/overlaps/motif_peak_overlaps/NRF1_filtered_6mer.MA0506.3_peak_overlaps.bed)

Rscript -e '
library(VennDiagram)
library(grid)

motifs  <- '"$n_motifs_NRF1"'
peaks   <- '"$n_peaks"'
overlap <- '"$n_overlap_NRF1"'

venn.obj <- draw.pairwise.venn(
  area1      = motifs,
  area2      = peaks,
  cross.area = overlap,
  category   = c("NRF1 6mer motifs", "ATAC peaks"),
  fill       = c("salmon", "skyblue"),
  alpha      = 0.7,
  cex        = 2,
  cat.cex    = 2.2,
  fontfamily = "Helvetica",
  cat.pos    = c(-15, 20),
  cat.dist   = c(0.03, 0.02),
  cat.fontfamily = "Helvetica",
  print.mode = c("raw","percent"),
  sigdig     = 2,
  ind        = FALSE)

png("./results/summarys/NRF1_motifs_vs_peaks_venn.png", width = 1300, height = 1300)
grid.newpage()
# Add title
grid.text(
  "Overlap between NRF1 6mer motifs and ATAC peaks",
  x = 0.5, y = 0.96,
  gp = gpar(fontsize = 36,fontfamily = "Helvetica"))

pushViewport(viewport(y = 0.45))
grid.draw(venn.obj)
dev.off()
'

# Create a Venn diagram to visualize the overlaps between motifs and peaks for BANP 

# Count total number of 6mer motifs
n_motifs_BANP=$(wc -l < ./motifs/ZBTB33_filtered_6mer.MA0527.2.bed)

# Count total number of peaks (e.g., merged peaks)
n_peaks=$(wc -l < ./peaks/merged_peaks.bed)

# Count number of 6mer BANP motifs that fall within peaks
n_overlap_BANP=$(wc -l < ./motifs/overlaps/motif_peak_overlaps/ZBTB33_filtered_6mer.MA0527.2_peak_overlaps.bed)

Rscript -e '
library(VennDiagram)
library(grid)

motifs  <- '"$n_motifs_BANP"'
peaks   <- '"$n_peaks"'
overlap <- '"$n_overlap_BANP"'

# Generate venn 
venn.obj <- draw.pairwise.venn(
  area1      = motifs,
  area2      = peaks,
  cross.area = overlap,
  category   = c("", ""),               
  fill       = c("salmon", "skyblue"),
  alpha      = 0.7,
  cex        = 1.9,       # size of the numbers inside the circles
  cat.cex    = 2.2,
  fontfamily = "Helvetica",   # font family for numbers
  print.mode = c("raw","percent"),
  sigdig     = 2,
  ind        = FALSE)

png("./results/summarys/BANP_motifs_vs_peaks_venn.png", width = 1300, height = 1300)
grid.newpage()

#Title
grid.text("Overlap between BANP 6mer motifs and ATAC peaks",
          x = 0.5, y = 0.95,
          gp = gpar(fontsize = 36,fontfamily = "Helvetica"))
# place the category labels above the circles
grid.text("ATAC peaks",
          x = 0.37, y = 0.82,           
          gp = gpar(fontsize = 28,fontfamily = "Helvetica"))
grid.text("BANP 6mer motifs",
          x = 0.75, y = 0.70,         
          gp = gpar(fontsize = 28,fontfamily = "Helvetica"))
# Venn diagram
pushViewport(viewport(y = 0.40))
grid.draw(venn.obj)

dev.off()
'

# 8) Overlap the intersected methylation motifs with the peak overlapping motifs to find common regions
mkdir -p ./motifs/overlaps/intersected_overlaps

bedtools intersect -u -a ./motifs/overlaps/intersected_motifs_HM450/NRF1_intersected_methylation.bed -b ./motifs/overlaps/motif_peak_overlaps/NRF1_filtered_6mer.MA0506.3_peak_overlaps.bed > ./motifs/overlaps/intersected_overlaps/NRF1_methylation_peak_overlap.bed  
bedtools intersect -u -a ./motifs/overlaps/intersected_motifs_HM450/ZBTB33_intersected_methylation.bed -b ./motifs/overlaps/motif_peak_overlaps/ZBTB33_filtered_6mer.MA0527.2_peak_overlaps.bed > ./motifs/overlaps/intersected_overlaps/ZBTB33_methylation_peak_overlap.bed

# Summarize the intersected results : NRF1 24325 ZBTB33 4863

# Count lines in bash and save to variables
n_NRF1=$(cat ./motifs/overlaps/intersected_overlaps/NRF1_methylation_peak_overlap.bed | wc -l)
n_ZBTB33=$(cat ./motifs/overlaps/intersected_overlaps/ZBTB33_methylation_peak_overlap.bed | wc -l)
# Create a summary table using R
Rscript -e '
table <- data.frame(
  Motif = c("NRF1","ZBTB33"),
  Overlapping_regions = c('"$n_NRF1"','"$n_ZBTB33"'))
print(table)
'
# create a Venn diagram to visualize the overlaps between motifs in peaks and motifs on HM450 for NRF1 and ZBTB33

Rscript -e '
library(VennDiagram)
library(grid)

# Motifs NRF1 in the peaks ATAC
motif_peak <- read.table("./motifs/overlaps/motif_peak_overlaps/NRF1_filtered_6mer.MA0506.3_peak_overlaps.bed")
peak_ids <- with(motif_peak, paste(V1, V2, V3, sep=":"))

# Motifs NRF1 in methylation probes HM450
motif_hm450 <- read.table("./motifs/overlaps/intersected_motifs_HM450/NRF1_intersected_methylation.bed")
hm450_ids <- with(motif_hm450, paste(V1, V2, V3, sep=":"))

# Venn : 
nrf1_lists <- list(
  "NRF1 motifs in ATAC peaks" = peak_ids,
  "NRF1 motifs in HM450"      = hm450_ids)

venn.plot <- venn.diagram(
  x = nrf1_lists,
  category.names = c("In ATAC peaks", "In HM450"),
  filename = NULL,
  fill = c("salmon", "skyblue"),
  alpha = 0.7,
  main = "NRF1 6mer motifs: ATAC peaks vs HM450",
  cex = 2,
  cat.cex = 2,
  main.cex = 3,
  cat.pos  = c(-20, -10), 
  cat.dist = c(0.02, 0.015),
  print.mode = c("raw", "percent"))

png("./results/summarys/NRF1_peaks_vs_HM450_venn.png", width = 1000, height = 1000)
grid.draw(venn.plot)
dev.off()
'
rm ./VennDiagram*.log  # remove temporary files created by VennDiagram package

Rscript -e '
library(VennDiagram)
library(grid)

## 1) Motifs ZBTB33 dans les pics ATAC
motif_peak <- read.table("./motifs/overlaps/motif_peak_overlaps/ZBTB33_filtered_6mer.MA0527.2_peak_overlaps.bed")
peak_ids <- with(motif_peak, paste(V1, V2, V3, sep=":"))

## 2) Motifs ZBTB33 qui overlappent HM450
motif_hm450 <- read.table("./motifs/overlaps/intersected_motifs_HM450/ZBTB33_intersected_methylation.bed")
hm450_ids <- with(motif_hm450, paste(V1, V2, V3, sep=":"))

## 3) Venn à 2 ensembles : 
##    A = motifs in ATAC peaks
##    B = motifs on HM450
zbtb33_lists <- list(
  "ZBTB33 motifs in ATAC peaks" = peak_ids,
  "ZBTB33 motifs on HM450"      = hm450_ids)

venn.plot <- venn.diagram(
  x = zbtb33_lists,
  category.names = c("In ATAC peaks", "In HM450"),
  filename = NULL,
  fill = c("salmon", "skyblue"),
  alpha = 0.7,
  main = "BANP 6mer motifs: ATAC peaks vs HM450",
  cex = 2,
  cat.cex = 2,
  main.cex = 3,
  cat.pos  = c(-20, -10), 
  cat.dist = c(0.02, 0.015),
  print.mode = c("raw", "percent"))

png("./results/summarys/BANP_peaks_vs_HM450_venn.png", width = 1000, height = 1000)
grid.draw(venn.plot)
dev.off()
'
rm ./VennDiagram*.log  # remove temporary files created by VennDiagram package

# create a barplot : motifs , peaks and HM450 for NRF1 
Rscript -e '
library(ggplot2)

## 1) Read motif sets and build unique IDs (chr:start:end)
motif_all   <- read.table("./motifs/NRF1_filtered_6mer.MA0506.3.bed")
motif_ids   <- unique(with(motif_all, paste(V1, V2, V3, sep=":")))

motif_peak  <- read.table("./motifs/overlaps/motif_peak_overlaps/NRF1_filtered_6mer.MA0506.3_peak_overlaps.bed")
peak_ids    <- unique(with(motif_peak, paste(V1, V2, V3, sep=":")))

motif_hm450 <- read.table("./motifs/overlaps/intersected_motifs_HM450/NRF1_intersected_methylation.bed")
hm450_ids   <- unique(with(motif_hm450, paste(V1, V2, V3, sep=":")))

## 2) Counts
total   <- length(motif_ids)
in_atac <- length(peak_ids)
in_hm450<- length(hm450_ids)
in_both <- length(intersect(peak_ids, hm450_ids))

## 3) Build data frame with percentages + pretty labels
df <- data.frame(
  category = c("All NRF1 motifs",
               "Motifs in ATAC peaks",
               "Motifs with HM450 probes",
               "Motifs in both ATAC + HM450"),
  count = c(total, in_atac, in_hm450, in_both))

df$category <- factor(
  df$category,
  levels = c("All NRF1 motifs",
             "Motifs in ATAC peaks",
             "Motifs with HM450 probes",
             "Motifs in both ATAC + HM450"))

df$percent <- df$count / total * 100
df$label   <- sprintf("%d (%.2f%%)", df$count, df$percent)

## 4) Plot
p <- ggplot(df, aes(x = reorder(category, -count), y = count, fill = category)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = label), vjust = -0.4, size = 4) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 20, hjust = 1),
    plot.title = element_text(hjust = 0.5),
  plot.subtitle = element_text(hjust = 0.5)
  ) +
  labs(
    title = "NRF1 6-mer motif distribution across datasets",
    subtitle = sprintf("all NRF1 6-mer motifs (n = %d)", total),
    y = "Number of motifs",
    x = NULL)

print(p)

ggsave("./results/summarys/NRF1_3D_vs_peaks_vs_HM450_barplot.pdf",plot = p, width = 10, height = 8, dpi = 300, bg = "white")
'

# create a barplot : motifs , peaks and HM450 for BANP 
Rscript -e '
library(ggplot2)

## 1) Read motif sets and build unique IDs (chr:start:end)
motif_all   <- read.table("./motifs/ZBTB33_filtered_6mer.MA0527.2.bed")
motif_ids   <- unique(with(motif_all, paste(V1, V2, V3, sep=":")))

motif_peak  <- read.table("./motifs/overlaps/motif_peak_overlaps/ZBTB33_filtered_6mer.MA0527.2_peak_overlaps.bed")
peak_ids    <- unique(with(motif_peak, paste(V1, V2, V3, sep=":")))

motif_hm450 <- read.table("./motifs/overlaps/intersected_motifs_HM450/ZBTB33_intersected_methylation.bed")
hm450_ids   <- unique(with(motif_hm450, paste(V1, V2, V3, sep=":")))

## 2) Counts
total   <- length(motif_ids)
in_atac <- length(peak_ids)
in_hm450<- length(hm450_ids)
in_both <- length(intersect(peak_ids, hm450_ids))

## 3) Build data frame with percentages + pretty labels
df <- data.frame(
  category = c("All BANP motifs",
               "Motifs in ATAC peaks",
               "Motifs with HM450 probes",
               "Motifs in both ATAC + HM450"),
  count = c(total, in_atac, in_hm450, in_both))

df$category <- factor(
  df$category,
  levels = c("All BANP motifs",
             "Motifs in ATAC peaks",
             "Motifs with HM450 probes",
             "Motifs in both ATAC + HM450"))

df$percent <- df$count / total * 100
df$label   <- sprintf("%d (%.2f%%)", df$count, df$percent)

## 4) Plot
p <- ggplot(df, aes(x = reorder(category, -count), y = count, fill = category)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = label), vjust = -0.4, size = 4) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 20, hjust = 1),
    plot.title = element_text(hjust = 0.5),
  plot.subtitle = element_text(hjust = 0.5)
  ) +
  labs(
    title = "BANP 6-mer motif distribution across datasets",
    subtitle = sprintf("all BANP 6-mer motifs (n = %d)", total),
    y = "Number of motifs",
    x = NULL)
print(p)

ggsave("./results/summarys/BANP_3D_vs_peaks_vs_HM450_barplot.pdf",plot = p, width = 10, height = 8, dpi = 300, bg = "white")
'



conda deactivate

# 9) Linking the methylation data from /data/papers/tcga to the methylation folder 
mkdir -p ./methylation/methylation_data/

for file in /data/papers/tcga/TCGA*/*/HM450*.txt; do
    ln -s "$file" ./methylation/methylation_data/
done

# 10) Linking the methylation data with the HM450 manifest file for annotation : chr start end probe_id(gc) beta_value

mkdir -p ./methylation/annotated_methylation_data/
# for each methylation file, annotate it and save it as a gzipped bed file in the same input directory
total=$(ls /data/papers/tcga/TCGA*/*/HM450*.txt | wc -l)
count=0

for file in /data/papers/tcga/TCGA*/*/HM450*.txt; do
    count=$((count+1))
    sample_name=$(basename "$file" .txt)
    echo "[$count / $total] Processing file: $file"
    if [ -f "$(dirname "$file")/${sample_name}_annotated_methylation.bed.gz" ] \
    && [ "$(zcat "$(dirname "$file")/${sample_name}_annotated_methylation.bed.gz" | wc -l)" -ne 0 ]; then # dirname to get the directory of the input file 
        echo "Skipping $(dirname "$file")/${sample_name}_annotated_methylation.bed.gz (already exists)"
        continue
    fi
    awk 'BEGIN{while(getline<"./methylation/annotated_methylation_data_probes.bed"){a[$4]=$0}}($1 in a){print a[$1],$2*100}' $file  | sort -k1,1 -k2,2n | gzip > "$(dirname "$file")/${sample_name}_annotated_methylation.bed.gz"
done

#for one file : 
#file=/data/papers/tcga/TCGA-STAD/TCGA-VQ-A91Z/HM450_TCGA-STAD-TCGA-VQ-A91Z-01A_1.txt
#sample_name=$(basename "$file" .txt)
#awk 'BEGIN{while(getline<"./methylation/annotated_methylation_data_probes.bed"){a[$4]=$0}}($1 in a){print a[$1],$2*100}' $file  | sort -k1,1 -k2,2n | gzip > "$(dirname "$file")/${sample_name}_annotated_methylation.bed.gz"

#to see how many were generated : ls /data/papers/tcga/TCGA*/*/*_annotated_methylation.bed.gz | wc -l. # should be 9812 files

# 11) Verify that all generated files have the same length: 485569 lines
for f in /data/papers/tcga/TCGA*/*/*_annotated_methylation.bed.gz; do
    echo -n "$f : "
    zcat "$f" | wc -l
done > ./methylation/annotated_methylation_file_line_counts.txt

cat ./methylation/annotated_methylation_file_line_counts.txt | grep 485569 | wc -l # should be 9812 

# 12) Statistics for Peaks : nbr of cancer types with peaks, nbr of peaks per cancer type, total nbr of peaks

# nbr of cancer types with peaks --> 23 types 
ls ./peaks/ATAC_TCGA*bed | cut -d'-' -f2 | cut -d'_' -f1 | sort | uniq -c | awk '{print $2 "\t" $1}' > ./peaks/cancer_counts.tsv

conda activate ./conda_envs/TF_binding_cancer_env

# Pie chart of the distribution of peaks per cancer type
Rscript -e '
library(ggplot2)

# Load TSV file and assign column names
df <- read.table("./peaks/cancer_counts.tsv",col.names = c("cancer_type", "count"))

# Total count
total_n <- sum(df$count)

# Compute percentages
df$perc <- df$count / total_n * 100

# Legend labels: "BRCA (n=75, 18.3%)"
df$legend_label <- paste0(df$cancer_type," (n=", df$count, ", ",round(df$perc, 1), "%)")

# Pie chart 
p <- ggplot(df, aes(x = "", y = count, fill = cancer_type)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(hjust = 0.5, size = 14),  # centered title
    plot.title.position = "plot",
    legend.title = element_text(size = 9),
    legend.text  = element_text(size = 6)) +
  ggtitle(paste0(
    "Distribution of ATAC-seq Peaks per Cancer Type\n",
    "Total cases: ", total_n)) +
  scale_fill_discrete(
    name = "Cancer type (n : nbr of peaks, %)",
    labels = df$legend_label,
    breaks = df$cancer_type)

# Save PDF
ggsave("./results/summarys/pie_chart_cancer_peaks_count.pdf",plot = p,width = 8,height = 5,dpi = 300,bg = "white")
'
# Plot histogram of nbr of Peaks per cancer type
# count peaks number per cancer type we have to merge all peak files per cancer type first and take the count of merged peaks


# Get list of cancer types from peak filenames
# Filenames look like: ATAC_TCGA-BRCA_TCGA-AR-A0TP_1_peaks_macs.bed
mkdir -p ./peaks/peaks_counts_per_cancer_type/

# 0) empty the output file (or create it)
> ./peaks/peaks_counts_per_cancer_type/peaks_unique_per_cancer.tsv

# 1) list cancer codes from filenames
#    ATAC_TCGA-BRCA_TCGA-AR-A0TP_1_peaks_macs.bed  →  TCGA-BRCA
cancers=$(ls ./peaks/ATAC_TCGA-*_*_peaks_macs.bed \
          | sed 's#.*/ATAC_##; s/_.*//' \
          | sort -u)

# 2) for each cancer, merge all peaks and count merged intervals
for cancer in $cancers; do
    echo "→ Processing $cancer ..."

    count=$(cat ./peaks/ATAC_${cancer}_*_peaks_macs.bed \
            | cut -f1-3 \
            | sort -k1,1 -k2,2n \
            | bedtools merge -i stdin \
            | wc -l)

    printf "%s\t%d\n" "$cancer" "$count" >> ./peaks/peaks_counts_per_cancer_type/peaks_unique_per_cancer.tsv
done

cut -d "-" -f2 ./peaks/peaks_counts_per_cancer_type/peaks_unique_per_cancer.tsv > ./peaks/peaks_counts_per_cancer_type/tmp.tsv
mv ./peaks/peaks_counts_per_cancer_type/tmp.tsv ./peaks/peaks_counts_per_cancer_type/peaks_unique_per_cancer.tsv

# create a cancer color order file based on the order of cancers in the peaks and snv unique per cancer files
cat ./snv/snv_counts_per_cancer_type/snv_unique_per_cancer.tsv \
    ./peaks/peaks_counts_per_cancer_type/peaks_unique_per_cancer.tsv \
    | cut -f1 \
    | sort -u \
    > ./results/summarys/cancer_color_order.txt


Rscript -e '

data <- read.table("./peaks/peaks_counts_per_cancer_type/peaks_unique_per_cancer.tsv")
                   
cancer_names <- data[,1]
peak_counts  <- data[,2]
# Sort from biggest to smallest
ord <- order(peak_counts, decreasing = TRUE)
cancer_names <- cancer_names[ord]
peak_counts  <- peak_counts[ord]
# reverse order for horizontal barplot (largest at top)
cancer_names <- rev(cancer_names)
peak_counts  <- rev(peak_counts)
# Convert to thousands
peak_thousands <- peak_counts / 1e3
# Build consistent colors
all_cancers <- scan("./results/summarys/cancer_color_order.txt", what = "")
palette <- rainbow(length(all_cancers))
names(palette) <- all_cancers
bar_colors <- palette[cancer_names]

pdf("./results/summarys/peaks_per_cancer_type_barplot.pdf", width=12, height=9 , bg = "white")

barplot(
  peak_thousands,
  names.arg=cancer_names,
  horiz=TRUE,                          # horizontal bars
  las=1,                               # labels readable
  col = bar_colors,
  border="black",
  main="Peaks per Cancer Type",
  cex.main = 3,
  cex.lab = 1.5,
  xlab="Total Peaks (thousands)")

dev.off()
'

# Distance of peaks to TSS : generate a file per sample with the distances of each peak to the nearest TSS
mkdir -p ./peaks/dist_to_tss/

for f in ./peaks/ATAC_TCGA*peaks_macs.bed; do
    sample=$(basename "$f" _peaks_macs.bed)

    echo "Processing $sample ..."

    closestBed -t first -d -a "$f" -b /data/genome/annotations/hg38_tss.bed | awk '{print $NF}' > "./peaks/dist_to_tss/${sample}_dist_tss.txt"
done

# Plot histogram of distances to TSS for all samples combined

Rscript -e 'library(ggplot2);
peak_files <- list.files("peaks", pattern = "_peaks_macs.bed$", full.names = TRUE);
samples <- sub("_peaks_macs.bed$", "", basename(peak_files));
pdf("./results/summarys/dist_tss_peaks_hist.pdf", width = 11.7, height = 8.3);
par(bg = "white");
for (sample in samples) {
  dist_file <- paste0("./peaks/dist_to_tss/", sample, "_dist_tss.txt");
  distances <- read.table(dist_file)[,1];

# compute proximal/distal percentages
  proximal <- sum(distances < 2000)
  distal   <- sum(distances >= 2000)
  total    <- length(distances)
  pct_prox <- round((proximal / total) * 100, 1)
  pct_dist <- round((distal   / total) * 100, 1)

  hist(log10(distances + 1),
       xlim = c(0, 8),
       xlab = "Distance to TSS (log10)",
       main = paste0(sample," | Proximal = ", pct_prox, "%| Distal = ", pct_dist, "%"),
       breaks = 50,
       col = "steelblue",
       border = "black");
  # red vertical line at log10(2000) --> promotor 
  abline(v = log10(2000), col = "red", lwd = 2)};
dev.off();'

# plot a boxplot of the percentage of proximal vs distal peaks across all samples
Rscript -e '
library(ggplot2)

# Load all samples 
peak_files <- list.files("peaks", pattern = "_peaks_macs.bed$", full.names = TRUE)
samples <- sub("_peaks_macs.bed$", "", basename(peak_files))

# Data frame to store percentages
df <- data.frame(
  sample = character(),
  region = character(),
  percent = numeric(),
  stringsAsFactors = FALSE)

# Compute percentages for each sample 
for (sample in samples) {
  dist_file <- paste0("./peaks/dist_to_tss/", sample, "_dist_tss.txt")
  distances <- read.table(dist_file)[,1] 
  proximal <- sum(distances < 2000)
  distal   <- sum(distances >= 2000)
  total    <- length(distances)
  pct_prox <- (proximal / total) * 100
  pct_dist <- (distal   / total) * 100

  # Add to dataframe
  df <- rbind(df,data.frame(sample = sample, region = "Proximal", percent = pct_prox),data.frame(sample = sample, region = "Distal",percent = pct_dist))}

# Boxplot of percentages 
pdf("./results/summarys/boxplot_ATAC_peaks_proximal_distal.pdf", width = 7, height = 5)
ggplot(df, aes(x = region, y = percent, fill = region)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.5) +
  theme_bw(base_size = 14) +
  labs(title = "Distribution of ATAC Peak Distances Across Samples",x = "", y = "Percentage of peaks")
dev.off()
'


# 13) statistics methylation probes per sample : 

# Plot histogram of methylation percentage for one sample and highlight UMR, LMR, FMR regions with different colors and add legend with percentages of each region.

# /data/papers/tcga/TCGA-BRCA/TCGA-AR-A0TP/HM450_TCGA-BRCA-TCGA-AR-A0TP-01A_1_annotated_methylation.bed.gz
# /data/papers/tcga/TCGA-ACC/TCGA-OR-A5J1/HM450_TCGA-ACC-TCGA-OR-A5J1-01A_1_annotated_methylation.bed.gz
# /data/papers/tcga/TCGA-CHOL/TCGA-3X-AAV9/HM450_TCGA-CHOL-TCGA-3X-AAV9-01A_1_annotated_methylation.bed.gz
# /data/papers/tcga/TCGA-GBM/TCGA-OX-A56R/HM450_TCGA-GBM-TCGA-OX-A56R-01A_1_annotated_methylation.bed.gz

Rscript -e '
data <- read.table("/data/papers/tcga/TCGA-BRCA/TCGA-AR-A0TP/HM450_TCGA-BRCA-TCGA-AR-A0TP-01A_1_annotated_methylation.bed.gz")
meth <- data[, ncol(data)]   # last column = methylation %

# Compute frequencies for each region
total <- length(meth)
n_UMR <- sum(meth >= 0  & meth < 10)
n_LMR <- sum(meth >= 10 & meth < 50)
n_FMR <- sum(meth >= 50 & meth <= 100)
pct_UMR <- round(100 * n_UMR / total, 1)
pct_LMR <- round(100 * n_LMR / total, 1)
pct_FMR <- round(100 * n_FMR / total, 1)

png("./results/summarys/Methylation_distribution_TCGA-BRCA-TCGA-AR-A0TP-01A_1.png", width=1000, height=800)

# Base histogram with no color 
hist(meth,
     breaks=50,
     col="white",
     border="black",
     main="Methylation distribution of TCGA-BRCA-TCGA-AR-A0TP-01A_1 sample",
     xlab="Methylation (%)")

# Add colored background rectangles for each region
usr <- par("usr")  # plot boundaries: (x1,x2,y1,y2)

# UMR: 0–10%
rect(0, usr[3], 10, usr[4], col=rgb(0,0,1,0.35), border=NA)

# LMR: 10–50%
rect(10, usr[3], 50, usr[4], col=rgb(1,0,0,0.35), border=NA)

# FMR: 50–100%
rect(50, usr[3], 100, usr[4], col=rgb(0,1,0,0.35), border=NA)

# Redraw histogram lines on top
hist(meth,
     breaks=50,
     col=NA,
     border="black",
     add=TRUE)

# Add vertical separation lines
abline(v = 10, col="black", lwd=2)
abline(v = 50, col="black", lwd=2)

# Legend with percentages
legend("topright",
       legend=c(
         paste0("UMR: 0–10% (", pct_UMR, "%)"),
         paste0("LMR: 10–50% (", pct_LMR, "%)"),
         paste0("FMR: 50–100% (", pct_FMR, "%)")),
       fill=c(rgb(0,0,1,0.25), rgb(1,0,0,0.25), rgb(0,1,0,0.25)),
       border=NA,
       cex=1.2)
dev.off()
'
# Compare UMR, LMR, FMR percentages across all samples per cancer type using a barplot
# 1) Create a summary table with percentages of UMR, LMR, FMR per sample and cancer type
mkdir -p ./methylation/methylation_counts/
# This file will contain: sample, cancer_type, pct_UMR, pct_LMR, pct_FMR only on cancer samples (not healthy)
files=(/data/papers/tcga/TCGA*/*/*_annotated_methylation.bed.gz)
total_files=${#files[@]}
count=0

for file in "${files[@]}"; do
    count=$((count + 1))
    echo "→ Processing sample $count / $total_files : $(basename "$file")"
    sample=$(basename "$file" _annotated_methylation.bed.gz)
    cancer_type=$(echo "$sample" | cut -d'-' -f2)
    # Extract TCGA sample type (01=tumor, 11=normal, 06=metastasis…)
    type=$(echo "$sample" | cut -d'-' -f4 | cut -c1-2)
    # Skip healthy 
    [[ "$type" == "11" ]] && continue
    # Proceed only with TUMOR samples
    meth_values=$(zcat "$file" | awk '{print $NF}')
    total=$(echo "$meth_values" | wc -l)
    n_UMR=$(echo "$meth_values" | awk '$1 >= 0 && $1 < 10' | wc -l)
    n_LMR=$(echo "$meth_values" | awk '$1 >= 10 && $1 < 50' | wc -l)
    n_FMR=$(echo "$meth_values" | awk '$1 >= 50 && $1 <= 100' | wc -l)
    pct_UMR=$(awk -v n="$n_UMR" -v t="$total" 'BEGIN {printf "%.2f", (n / t) * 100}')
    pct_LMR=$(awk -v n="$n_LMR" -v t="$total" 'BEGIN {printf "%.2f", (n / t) * 100}')
    pct_FMR=$(awk -v n="$n_FMR" -v t="$total" 'BEGIN {printf "%.2f", (n / t) * 100}')
    echo -e "$sample\t$cancer_type\t$pct_UMR\t$pct_LMR\t$pct_FMR" >> ./methylation/methylation_counts/methylation_region_percentages_per_cancer_samples.tsv
done

# 2)add a header to the file
awk 'BEGIN{print "sample\tcancer_type\tpct_UMR\tpct_LMR\tpct_FMR"}{print}' ./methylation/methylation_counts/methylation_region_percentages_per_cancer_samples.tsv > ./methylation/methylation_counts/methylation_region_percentages_per_sample.tsv

# 3) combine the samples of the same cancer type and taking the median percentage for each region into a new tsv file 
Rscript -e '
df <- read.table("./methylation/methylation_counts/methylation_region_percentages_per_sample.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

medians_by_cancer <- aggregate(
  df[ , c("pct_UMR", "pct_LMR", "pct_FMR")],
  by = list(Cancer = df$cancer_type),
  FUN = median,
  na.rm = TRUE)

write.table(medians_by_cancer,file = "./methylation/methylation_counts/methylation_region_percentages_medians_by_cancer.tsv",sep = "\t", quote = FALSE, row.names = FALSE)
'

# 4) Plot barplot using ggplot2
Rscript -e '
library(ggplot2)
library(reshape2)
library(RColorBrewer)

df <- read.table("./methylation/methylation_counts/methylation_region_percentages_medians_by_cancer.tsv",header = TRUE)

df_long <- melt(df, id.vars = "Cancer",variable.name = "region", value.name = "percentage")

# reorder cancers alphabetically or by median total percentage
df_long$Cancer <- factor(df_long$Cancer, levels = sort(unique(df_long$Cancer)))

# color palette
palette <- brewer.pal(n = length(unique(df_long$region)), "Set2")

pdf("./results/summarys/methylation_region_percentages_barplot.pdf",width = 14, height = 8)

ggplot(df_long, aes(x = Cancer, y = percentage, fill = region)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = palette) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, size = 10),
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(hjust = 0.5)
) +
  labs(
    title = "Methylation Region Percentages per Cancer Type",
    x = "",
    y = "Percentage (%)")
dev.off()
'

# Compare the methylation distribution between cancer types using boxplot for UMR, LMR, FMR to visualize each point to see if we have outliers 
Rscript -e '
library(ggplot2)
library(reshape2)

# 1) Load per-sample percentages
df <- read.table("./methylation/methylation_counts/methylation_region_percentages_per_sample.tsv",header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 2) Go to long format: one line = (sample, cancer_type, region, percentage)
df_long <- melt(df,
                id.vars = c("sample", "cancer_type"),
                variable.name = "region",
                value.name   = "percentage")

# 3) Order cancer types alphabetically (for x-axis)
df_long$cancer_type <- factor(df_long$cancer_type,levels = sort(unique(df_long$cancer_type)))

# 4) Boxplot + jitter per region to see outliers ; jitter soo the points do not overlap and are visible
pdf("./results/summarys/methylation_region_percentages_boxplot_per_sample.pdf",width = 16, height = 9)

ggplot(df_long, aes(x = cancer_type, y = percentage, color = region)) +
  geom_boxplot(outlier.colour = NA) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 0.7) +
  facet_wrap(~ region, nrow = 3, scales = "fixed") +
  scale_color_brewer(palette = "Set2") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 60, hjust = 1, size = 7),
    strip.text = element_text(size = 12),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    title = "UMR / LMR / FMR methylation percentages per sample",
    x = "Cancer type",
    y = "Percentage (%)",
    color = "Region")

dev.off()
'

# Distance of methylation sites to TSS : generate a file per sample with the distances of each methylation to the nearest TSS

mkdir -p ./methylation/dist_to_tss/
file='/data/papers/tcga/TCGA-GBM/TCGA-OX-A56R/HM450_TCGA-GBM-TCGA-OX-A56R-01A_1_annotated_methylation.bed.gz'
closestBed -t first -d -a "$file" -b /data/genome/annotations/hg38_tss.bed | awk '{print $NF}' > "./methylation/dist_to_tss/TCGA-GBM-TCGA-OX-A56R-01A_1_dist_tss.txt"

# Plot histogram of distances to TSS for one sample 

Rscript -e '
dist_file <- "./methylation/dist_to_tss/TCGA-BRCA-TCGA-AR-A0TP-01A_1_dist_tss.txt"
distances <- read.table(dist_file)[,1]

png("./results/summarys/dist_tss_methylation_TCGA-BRCA-TCGA-AR-A0TP-01A_1_hist.png", width=1000, height=800)
par(bg = "white")

# compute proximal/distal percentages
  proximal <- sum(distances < 2000)
  distal   <- sum(distances >= 2000)
  total    <- length(distances)
  pct_prox <- round((proximal / total) * 100, 1)
  pct_dist <- round((distal   / total) * 100, 1)

hist(log10(distances + 1),
     xlim = c(0, 8),
     xlab = "Distance to TSS (log10)",
     main = paste0("TCGA-BRCA-TCGA-AR-A0TP-01A_1", "  | Proximal = ", pct_prox, "% Distal = ", pct_dist, "%"),
     breaks = 50,
     col = "steelblue",
     border = "black")
  abline(v = log10(2000), col = "red", lwd = 2)
dev.off()
'

# number of methylation probes per cancer type
mkdir -p ./methylation/methylation_counts/
for file in /data/papers/tcga/TCGA*/*/*_annotated_methylation.bed.gz; do
    sample=$(basename "$file" _annotated_methylation.bed.gz)
    count=$(zcat "$file" | wc -l)
    echo -e "$sample\t$count" >> ./methylation/methylation_counts/methylation_probes_per_sample.tsv
done

# Comparing the methylation distribution for the same sample between healthy and cancer samples

# 2>/dev/null = hide all error messages from this command.

#to find samples with healthy and cancer data : 
for dir in /data/papers/tcga/TCGA-BRCA/*/; do
    # Any tumor sample (01,02,03,05,06,07,08,09)
    tumor_files=( "$dir"/HM450*-0[1-9]*.bed.gz )
    # Any normal sample (10,11,12,13)
    normal_files=( "$dir"/HM450*-1[0-3]*.bed.gz )
    # Check that both tumor and normal exist
    if [ -e "${tumor_files[0]}" ] && [ -e "${normal_files[0]}" ]; then
        echo ">>> $dir has tumor + normal HM450 methylation:"
        echo "Tumor:"
        ls "${tumor_files[@]}"
        echo "Normal:"
        ls "${normal_files[@]}"
        echo ""
    fi
done

# to find samples with no healthy data :
mkdir -p ./methylation/methylation_counts/
mkdir -p ./methylation/methylation_counts/

# (optional) add a header once
echo -e "cancer\ttumor_count\thealthy_count" > ./methylation/methylation_counts/methylation_presence.tsv

for cancer_dir in /data/papers/tcga/TCGA-*; do
    cancer=$(basename "$cancer_dir" | sed 's/^TCGA-//')

    # Any tumor sample: 0[1-9] = 01,02,03,05,06,07,08,09...
    tumor_count=$(find "$cancer_dir" -type f \
        -name "HM450*-0[1-9]*_1_annotated_methylation.bed.gz" | wc -l)

    # Any normal sample: 1[0-3] = 10,11,12,13
    healthy_count=$(find "$cancer_dir" -type f \
        -name "HM450*-1[0-3]*_1_annotated_methylation.bed.gz" | wc -l)

    if [ "$tumor_count" -gt 0 ] && [ "$healthy_count" -gt 0 ]; then
        echo "$cancer : has BOTH tumor (0[1-9]) and healthy (1[0-3]) methylation"
    else
        echo "$cancer : tumor = $tumor_count, healthy = $healthy_count"
    fi

    printf "%s\t%d\t%d\n" "$cancer" "$tumor_count" "$healthy_count" \
        >> ./methylation/methylation_counts/methylation_presence.tsv
done


# Visualize the results in a barplot using R
Rscript -e '
library(ggplot2)
library(reshape2)

# Load the table with counts
df <- read.table("./methylation/methylation_counts/methylation_presence.tsv", header = TRUE)
colnames(df) <- c("cancer", "tumor", "healthy")

# Add flag for "no healthy"
df$no_healthy <- df$healthy == 0

# Compute total samples per cancer type
df$total <- df$tumor + df$healthy

# Reorder cancers by total samples (descending)
df$cancer <- factor(df$cancer, levels = df$cancer[order(-df$total)])

# Build label dataframe AFTER reordering + after no_healthy exists
label_df <- df[, c("cancer", "no_healthy")]

# Convert to long format (only tumor + healthy)
df_long <- melt(df,
                id.vars = "cancer",
                measure.vars = c("tumor", "healthy"),
                variable.name = "type",
                value.name = "count")

# Create the barplot
pdf("./results/summarys/methylation_counts_barplot.pdf", width=12, height=8)

ggplot(df_long, aes(x = cancer, y = count, fill = type)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("tumor" = "salmon", "healthy" = "steelblue"), name = "") +
  labs(
    title   = "Tumor vs Healthy Methylation Sample Counts per Cancer Type",
    x       = "Cancer Type",
    y       = "Number of Samples",
    caption = "Red cancer names = cancer types with no healthy samples") +
  theme_minimal(base_size = 13) +
  theme(
    plot.title   = element_text(hjust = 0.5),
    axis.text.x  = element_blank(),     # hide default x labels
    axis.ticks.x = element_blank(),
    plot.caption = element_text(hjust = 0.5, color = "red", size = 12),
    plot.margin  = margin(t = 10, r = 20, b = 60, l = 20)  # more space for labels) +
  coord_cartesian(clip = "off") +   # prevent clipping
  geom_text(
    data = label_df,
    aes(
      x     = cancer,
      y     = 0,
      label = cancer,
      color = no_healthy),
    angle       = 45,
    hjust       = 1,
    vjust       = 1.2,
    inherit.aes = FALSE,
    size        = 4) +
  scale_color_manual(
    values = c(`TRUE` = "red", `FALSE` = "black"),
    guide  = "none")

dev.off()
'

#example of outputs from BRCA : 
#/data/papers/tcga/TCGA-BRCA/TCGA-A7-A0D9/ has 2 bed.gz files : 
# /data/papers/tcga/TCGA-BRCA/TCGA-A7-A0D9/HM450_TCGA-BRCA-TCGA-A7-A0D9-11A_1_annotated_methylation.bed.gz --> healthy
# /data/papers/tcga/TCGA-BRCA/TCGA-A7-A0D9/HM450_TCGA-BRCA-TCGA-A7-A0D9-01A_1_annotated_methylation.bed.gz --> cancer
# /data/papers/tcga/TCGA-BRCA/TCGA-A7-A13F//HM450_TCGA-BRCA-TCGA-A7-A13F-01A_1_annotated_methylation.bed.gz
# /data/papers/tcga/TCGA-BRCA/TCGA-A7-A13F//HM450_TCGA-BRCA-TCGA-A7-A13F-11A_1_annotated_methylation.bed.gz

# between two replicates of cancer :
# /data/papers/tcga/TCGA-BRCA/TCGA-A7-A26J//HM450_TCGA-BRCA-TCGA-A7-A26J-01A_1_annotated_methylation.bed.gz
# /data/papers/tcga/TCGA-BRCA/TCGA-A7-A26J//HM450_TCGA-BRCA-TCGA-A7-A26J-01A_2_annotated_methylation.bed.gz

# /data/papers/tcga/TCGA-ESCA/TCGA-IC-A6RE/ has 2 HM450 methylation files
# /data/papers/tcga/TCGA-ESCA/TCGA-IC-A6RE//HM450_TCGA-ESCA-TCGA-IC-A6RE-01A_1_annotated_methylation.bed.gz
# /data/papers/tcga/TCGA-ESCA/TCGA-IC-A6RE//HM450_TCGA-ESCA-TCGA-IC-A6RE-11A_1_annotated_methylation.bed.gz
# /data/papers/tcga/TCGA-ESCA/TCGA-L5-A4OQ/ has 2 HM450 methylation files
# /data/papers/tcga/TCGA-ESCA/TCGA-L5-A4OQ//HM450_TCGA-ESCA-TCGA-L5-A4OQ-01A_1_annotated_methylation.bed.gz
# /data/papers/tcga/TCGA-ESCA/TCGA-L5-A4OQ//HM450_TCGA-ESCA-TCGA-L5-A4OQ-11A_1_annotated_methylation.bed.gz

# plot the methylation distribution for both healthy and cancer samples of TCGA-BRCA-TCGA-A7-A0D9 using a delta methylation histogram :
Rscript -e '
library(ggplot2)
# Load files (same format: chr start end probe meth)
cancer  <- read.table("/data/papers/tcga/TCGA-BRCA/TCGA-A7-A13F//HM450_TCGA-BRCA-TCGA-A7-A13F-01A_1_annotated_methylation.bed.gz",  header=FALSE)
healthy <- read.table("/data/papers/tcga/TCGA-BRCA/TCGA-A7-A13F//HM450_TCGA-BRCA-TCGA-A7-A13F-11A_1_annotated_methylation.bed.gz", header=FALSE)

# Extract probe ID and methylation (col 4 and 5)
cancer_df  <- data.frame(probe=cancer[,4],  meth_cancer=cancer[,5])
healthy_df <- data.frame(probe=healthy[,4], meth_healthy=healthy[,5])

# Merge on probe ID
merged <- merge(cancer_df, healthy_df, by="probe")
cat("Cancer rows:", nrow(cancer_df), "\n")
cat("Healthy rows:", nrow(healthy_df), "\n")
cat("Merged rows:", nrow(merged), "\n")

# Compute delta methylation
delta <- merged$meth_cancer - merged$meth_healthy

# Compute percentages
total <- length(delta)
pct_negative <- round(sum(delta < 0) / total * 100, 1)
pct_positive <- round(sum(delta > 0) / total * 100, 1)

png("./results/summarys/delta_methylation_hist_TCGA-BRCA-TCGA-A7-A13F.png", width=1000, height=800)

hist(delta,
     breaks = 100,
     col = "steelblue",
     border = "black",
     main = paste0("Δ Methylation (Cancer – Healthy TCGA-BRCA-TCGA-A7-A13F )\n",
                   pct_negative, "% CpGs hypomethylated (Δ<0) | ",
                   pct_positive, "% hypermethylated (Δ>0)"),
     xlab = "Δ methylation (%)")

# Add reference line at 0
abline(v = 0, col = "red", lwd = 2)

# Add text labels inside the plot
usr <- par("usr")
text(x = usr[1] + 5,
     y = usr[4] * 0.9,
     labels = paste0("Hypomethylated (Δ < 0): ", pct_negative, "%"),
     col = "blue", adj = 0, cex = 1.2)

text(x = usr[1] + 5,
     y = usr[4] * 0.82,
     labels = paste0("Hypermethylated (Δ > 0): ", pct_positive, "%"),
     col = "darkred", adj = 0, cex = 1.2)

dev.off()
'
# using a smooth scatter plot to visualize the correlation between healthy and cancer methylation levels for the same sample : 
Rscript -e '
library(KernSmooth)

# Load files
cancer  <- read.table("/data/papers/tcga/TCGA-BRCA/TCGA-A7-A0D9/HM450_TCGA-BRCA-TCGA-A7-A0D9-01A_1_annotated_methylation.bed.gz",  header=FALSE)
healthy <- read.table("/data/papers/tcga/TCGA-BRCA/TCGA-A7-A0D9/HM450_TCGA-BRCA-TCGA-A7-A0D9-11A_1_annotated_methylation.bed.gz", header=FALSE)

# Extract probe + methylation
cancer_df  <- data.frame(probe=cancer[,4],  meth_cancer=cancer[,5])
healthy_df <- data.frame(probe=healthy[,4], meth_healthy=healthy[,5])

# Merge
merged <- merge(cancer_df, healthy_df, by="probe")

# Delta methylation
delta <- merged$meth_cancer - merged$meth_healthy
# delta thresholds
thr <-10  # 10%

# classification
hypo  <- merged$meth_cancer < (merged$meth_healthy - thr)
hyper <- merged$meth_cancer > (merged$meth_healthy + thr)
equal <- !hypo & !hyper   # everything in between

# Hypo/hyper definitions 
n_hypo  <- sum(hypo,  na.rm=TRUE)
n_hyper <- sum(hyper, na.rm=TRUE)
n_equal <- sum(equal, na.rm=TRUE)

tot <- n_hypo + n_hyper + n_equal

pct_hypo  <- round(100 * n_hypo  / tot, 1)
pct_hyper <- round(100 * n_hyper / tot, 1)
pct_equal <- round(100 * n_equal / tot, 1)

pdf("./results/summarys/smoothScatter_delta_TCGA-BRCA-TCGA-A7-A0D9.pdf", width=8, height=8)

smoothScatter(
  merged$meth_healthy,
  merged$meth_cancer,
  xlab="Healthy methylation (%)",
  ylab="Cancer methylation (%)",
  main="Cancer vs Healthy Methylation\n TCGA-BRCA-TCGA-A7-A0D9. \n delta = cancer - healthy, hypo: delta < 10%, hyper: delta > 10%")


# Identity line
abline(0, 1, col="red", lwd=2)

# ±10% shift lines
abline(a = -10, b = 1, col = "blue", lwd = 2)   # 10% lower
abline(a =  10, b = 1, col = "blue", lwd = 2) # 10% higher

# Add annotation
legend("left", bty="n",cex = 0.8, legend=c(
  paste0("Hypomethylated: ", n_hypo, " (", pct_hypo, "%)"),
  paste0("Hypermethylated: ", n_hyper, " (", pct_hyper, "%)"),
  paste0("Unchanged: ", n_equal, " (", pct_equal, "%)")))
dev.off()
'

# create a list of all methylation files for cancer samples includes all tumor types and samples while excluding healthy samples and including only primary vial of tumor samples (A)
find /data/papers/tcga/TCGA-* -type f -name "HM450*_annotated_methylation.bed.gz" | grep -E -- '-0[0-9]A_' > ./methylation/methylation_files_tumor0x.txt


# SNV statistics : number of SNVs per sample
mkdir -p ./snv/snv_counts/
for file in ./snv/SNV_TCGA*vcf.gz; do
    counter=$((counter+1))   # increment counter
    sample=$(basename "$file" .vcf.gz | sed 's/SNV_//') # remove SNV_ prefix
    count=$(zcat "$file" | grep -v "^#" | wc -l) # count non-header lines --> number of SNVs
    echo -e "$counter\t$sample\t$count" >> ./snv/snv_counts/snv_per_sample.tsv
    echo "Processed $counter : $sample with $count SNVs"
done

# Plot histogram of SNVs per sample 
Rscript -e '
data <- read.table("./snv/snv_counts/snv_per_sample.tsv", header=FALSE)
snv_counts <- data[,3]   # 3rd column = SNV counts
png("./results/summarys/snv_per_sample_hist.png", width=1000, height=800)
par(bg = "white")

hist(log10(snv_counts + 1),
     xlim = c(0, 8),
     xlab = "SNVs per sample",
     main = "SNVs per sample",
     breaks = 50,
     col = "steelblue",
     border = "black")
dev.off()
'

# number of SNVs per cancer type with unique SNVs only (removing duplicates across samples of the same cancer type)
mkdir -p ./snv/snv_counts_per_cancer_type/
# unique SNVs per cancer type
# sed to extract cancer type from filename , first remove prefix and then removes everything after the first '-'
cancers=$(find ./snv -maxdepth 1 -name 'SNV_TCGA-*.vcf.gz' \
          | sed 's#.*/SNV_TCGA-##; s/-.*//' \
          | sort -u)

for cancer in $cancers; do
    echo "→ Processing $cancer..."
    count=$(zgrep -hv "^#" ./snv/SNV_TCGA-"$cancer"-*.vcf.gz \
            | awk '{print $1":"$2":"$4":"$5}' \
            | sort -u \
            | wc -l)
    echo "Done $cancer ($count unique SNVs)"
    printf "%s\t%s\n" "$cancer" "$count" >> ./snv/snv_counts_per_cancer_type/snv_unique_per_cancer.tsv
done

# Plot histogram of SNVs per cancer type
Rscript -e '

data <- read.table("./snv/snv_counts_per_cancer_type/snv_unique_per_cancer.tsv")

cancer_names <- data[,1]
snv_counts   <- data[,2]

# Sort from biggest to smallest
ord <- order(snv_counts, decreasing = TRUE)
cancer_names <- cancer_names[ord]
snv_counts  <- snv_counts[ord]
# reverse order for horizontal barplot (largest at top)
cancer_names <- rev(cancer_names)
snv_counts  <- rev(snv_counts)

snv_millions <- snv_counts / 1e6
# Build consistent colors
all_cancers <- scan("./results/summarys/cancer_color_order.txt", what = "")
palette <- rainbow(length(all_cancers))
names(palette) <- all_cancers
bar_colors <- palette[cancer_names]

pdf("./results/summarys/snv_per_cancer_type_barplot2.pdf",width=12, height=9, bg = "white")

nice_max <- max(pretty(snv_millions))
barplot(
  snv_millions,
  names.arg = cancer_names,
  horiz = TRUE,
  las = 1,
  col = bar_colors,
  border = "black",
  main = "SNVs per Cancer Type",
  cex.main = 3,
  cex.lab  = 1.5,
  xlab = "Total SNVs (millions)",
  xlim = c(0, nice_max))
dev.off()
'





# Distance of SNVs sites to TSS : generate a file per sample with the distances of each SNV to the nearest TSS
# ./snv/SNV_TCGA-BRCA-TCGA-AR-A0TP-01A_vs_TCGA-AR-A0TP-10A_1.vcf.gz
# ./snv/SNV_TCGA-ACC-TCGA-OR-A5J1-10A_vs_TCGA-OR-A5J1-01A_1.vcf.gz
# ./snv/SNV_TCGA-CHOL-TCGA-3X-AAV9-10A_vs_TCGA-3X-AAV9-01A_1.vcf.gz
# ./snv/SNV_TCGA-GBM-TCGA-OX-A56R-10A_vs_TCGA-OX-A56R-01A_1.vcf.gz 

# Convert VCF → BED-like (0-based start, 1-based end) and sort
mkdir -p ./snv/dist_to_tss/
zcat ./snv/SNV_TCGA-GBM-TCGA-OX-A56R-10A_vs_TCGA-OX-A56R-01A_1.vcf.gz  \
    | grep -v '^#' \
    | awk '{print $1"\t"$2-1"\t"$2}' \
    | sort -k1,1 -k2,2n \
    | closestBed -t first -d -a - -b /data/genome/annotations/hg38_tss.bed | awk '{print $NF}' > "./snv/dist_to_tss/TCGA-GBM-TCGA-OX-A56R-10A_1_dist_tss.txt"
# Plot histogram of distances to TSS for one sample 

Rscript -e '
dist_file <- "./snv/dist_to_tss/TCGA-BRCA-TCGA-AR-A0TP-01A_1_dist_tss.txt"
distances <- read.table(dist_file)[,1]

png("./results/summarys/dist_tss_snv_TCGA-BRCA-TCGA-AR-A0TP-01A_1_hist.png", width=1000, height=800)
par(bg = "white")

# compute proximal/distal percentages
  proximal <- sum(distances < 2000)
  distal   <- sum(distances >= 2000)
  total    <- length(distances)
  pct_prox <- round((proximal / total) * 100, 1)
  pct_dist <- round((distal   / total) * 100, 1)

hist(log10(distances + 1),
     xlim = c(0, 8),
     xlab = "Distance to TSS (log10)",
     main = paste0("TCGA-BRCA-TCGA-AR-A0TP-01A_1", "  |  Proximal = ", pct_prox, "%  | Distal = ", pct_dist, "%"),
     breaks = 50,
     col = "steelblue",
     border = "black")
     # red vertical line at log10(2000) --> promotor 
  abline(v = log10(2000), col = "red", lwd = 2)
dev.off()
'

# 15) Counting cancer types with both peaks and SNV data
mkdir -p ./results/summarys

# List of cancer types with methylation data (from the methylation file names)
ls ./methylation/methylation_data/ | cut -d'-' -f2 | sort | uniq > ./methylation/cancer_types_with_methylation.txt

# List of cancer types with ATAC peaks (cancer_counts.tsv)
cut -f1 ./peaks/cancer_counts.tsv | sort > ./peaks/cancer_types_with_peaks.txt

# List of cancer types with SNV data (snv_total_per_cancer.tsv)
cut -f1 ./snv/snv_counts_per_cancer_type/snv_total_per_cancer.tsv | sort > ./snv/cancer_types_with_snv.txt

# heatmap of cancer types with ATAC vs SNV vs Methylation data
Rscript -e 'library(ggplot2)

# Read lists (one cancer type per line)
atac <- scan("./peaks/cancer_types_with_peaks.txt", what="")
snv  <- scan("./snv/cancer_types_with_snv.txt", what="")
methylation <- scan("./methylation/cancer_types_with_methylation.txt", what="") 
# Combine
all_types <- sort(unique(c(atac, snv , methylation)))

# Build matrix manually (super simple)
df <- data.frame(
  cancer_type = all_types,
  ATAC = ifelse(all_types %in% atac, 1, 0),
  SNV  = ifelse(all_types %in% snv,  1, 0) ,
  Methylation = ifelse(all_types %in% methylation, 1, 0))

# Convert to long format manually (still simple)
df_long <- data.frame(
  cancer_type = rep(df$cancer_type, times = 3),
  dataset     = rep(c("ATAC", "Methylation", "SNV"), each = nrow(df)),
  available   = c(df$ATAC, df$Methylation, df$SNV))

# Heatmap with salmon = no, steelblue = yes
p <- ggplot(df_long, aes(x = dataset, y = cancer_type, fill = factor(available))) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("0" = "salmon", "1" = "steelblue"),labels = c("No", "Yes")) +
  theme_minimal() +
  theme(axis.title = element_blank(),
    axis.text.x = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5),
    axis.text.y = element_text(size = 8),
    legend.title = element_blank(),
    panel.grid = element_blank()) +
  ggtitle("Availability of ATAC, SNV, and Methylation Data per Cancer Type")

ggsave("./results/summarys/cancer_types_ATAC_SNV_Methylation_heatmap.png",p, width = 10, height = 10, dpi = 300, bg = "white")
'
# clean up intermediate files
rm ./snv/cancer_types_with_snv.txt
rm ./methylation/cancer_types_with_methylation.txt
rm ./peaks/cancer_types_with_peaks.txt

# 16) Overlap SNV positions with ATAC peak positions per sample and count the number of SNVs that overlap peaks per sample
mkdir -p ./snv/overlaps/snv_peak_overlap/

for snv_file in ./snv/SNV_TCGA-*.vcf.gz; do
    snv_base=$(basename "$snv_file" .vcf.gz | sed 's/^SNV_//')
  
    cancer=$(echo "$snv_base" | cut -d'-' -f1-2)          # TCGA-ACC
    rest=${snv_base#${cancer}-}                           # TCGA-OR-A5J2-01A_vs_TCGA-OR-A5J2-10A_1
    patient=$(echo "$rest" | cut -d'-' -f1-3)             # TCGA-OR-A5J2

    peak_file="./peaks/ATAC_${cancer}_${patient}_1_peaks_macs.bed"

    if [ -f "$peak_file" ]; then
        echo "Processing $snv_base with peaks $peak_file"
        # Convert VCF → BED
        zcat "$snv_file" | grep -v '^#' | awk '{print $1"\t"$2-1"\t"$2}' \
            | sort -k1,1 -k2,2n > "./snv/overlaps/snv_peak_overlap/${snv_base}_snv.bed"

        # Overlap SNVs with peaks
        bedtools intersect -u -a "./snv/overlaps/snv_peak_overlap/${snv_base}_snv.bed" -b "$peak_file" \
            | wc -l > "./snv/overlaps/snv_peak_overlap/${snv_base}_snv_in_peaks.txt"
        rm "./snv/overlaps/snv_peak_overlap/${snv_base}_snv.bed"
        echo "Done $snv_base"
    else
        echo "No peak file for $snv_base (looked for $peak_file), skipping."
    fi
done

# plot barplot of number of SNV in ATAC peak data
Rscript -e '
library(ggplot2)

# 1) Files with SNVs-in-peaks counts
overlap_files <- list.files("./snv/overlaps/snv_peak_overlap",pattern = "_snv_in_peaks.txt$",full.names = TRUE)

# extract sample names that actually have ATAC since only 410 ATAC peaks files exist 
samples_with_peaks <- sub("_snv_in_peaks.txt$", "", basename(overlap_files))
# basename() removes all folder paths → keeps only the filename
# sub() removes the suffix _snv_in_peaks.txt

# 2) Total SNVs for those samples only 
snv_counts <- read.table("./snv/snv_counts/snv_per_sample.tsv",header = FALSE, stringsAsFactors = FALSE)
colnames(snv_counts) <- c("idx","sample","count")
# subset to samples with ATAC data
snv_sub <- snv_counts[snv_counts$sample %in% samples_with_peaks, ]
total_snv <- sum(snv_sub$count)
cat("Total SNVs in samples with ATAC:", total_snv, "\n")

# 3) Sum SNVs in ATAC peaks : for each overlap file, read the count and sum them
counts <- vapply(
  overlap_files,
  function(f) {as.numeric(readLines(f, n = 1))}, # read a single line and convert to numeric 
  numeric(1)) # output type numeric of length 1

snv_in_peaks <- sum(counts)
cat("Total SNVs in ATAC peaks:", snv_in_peaks, "\n")

# 4) Barplot in vs out
df <- data.frame(category = c("SNVs in ATAC peaks", "SNVs outside ATAC peaks"),count = c(snv_in_peaks, total_snv - snv_in_peaks))
df$percent <- round(100 * df$count / sum(df$count), 2)

p <- ggplot(df, aes(x = category, y = count, fill = category)) +
  geom_col() +
  geom_text(aes(label = paste0(count, " (", percent, "%)")),vjust = -0.3) +
  theme_minimal() +
  labs(title = "Proportion of SNVs overlapping ATAC peaks",x = "", y = "Number of SNVs") +
  theme(axis.text.x = element_text(angle = 15, hjust = 1),legend.position = "none")

ggsave("./results/summarys/SNV_in_vs_out_ATAC_barplot.pdf",p, width = 7, height = 5, bg = "white")
'

# 17) Overlap SNVs with methylation sites per sample and count the number of SNVs that overlap methylation sites per sample
mkdir -p ./snv/overlaps/snv_methylation_overlap/

for snv_file in ./snv/SNV_TCGA-*.vcf.gz; do
    snv_base=$(basename "$snv_file" .vcf.gz)   # SNV_TCGA-ACC-TCGA-OR-A5J2-01A_vs_TCGA-OR-A5J2-10A_1
    
    # Remove SNV_ prefix for easier parsing
    core=${snv_base#SNV_}                      # TCGA-ACC-TCGA-OR-A5J2-01A_vs_TCGA-OR-A5J2-10A_1
    # cancer = TCGA-ACC
    cancer=$(echo "$core" | cut -d'-' -f1-2 )  # ACC
    # patient = TCGA-OR-A5J2
    patient=$(echo "$core" | cut -d'-' -f3-5)
    # tumor code = 01A (from "01A_vsTCGA-..." part)
    tumor_block=$(echo "$core" | cut -d'-' -f6)   # 01A_vs_TCGA-OR-A5J2-10A_1
    tumor_code=${tumor_block%%_*}                 # 01A

    # Path to annotated methylation file of the tumor sample
    methylation_file="/data/papers/tcga/${cancer}/${patient}/HM450_${cancer}-${patient}-${tumor_code}_1_annotated_methylation.bed.gz"
    
    if [ -f "$methylation_file" ]; then
        echo "Processing $snv_base with methylation $methylation_file"

        # Convert SNVs VCF -> BED (0-based start, 1-based end)
        zcat "$snv_file" | grep -v '^#' | awk '{print $1"\t"$2-1"\t"$2}'| sort -k1,1 -k2,2n > "./snv/overlaps/snv_methylation_overlap/${snv_base}_snv.bed"

        # Overlap SNVs with methylation probes
        bedtools intersect -u -a "./snv/overlaps/snv_methylation_overlap/${snv_base}_snv.bed" -b "$methylation_file" | wc -l > "./snv/overlaps/snv_methylation_overlap/${snv_base}_snv_in_methylation.txt"
        rm "./snv/overlaps/snv_methylation_overlap/${snv_base}_snv.bed"
        echo "Done $snv_base"
    else
        echo "No methylation file for $snv_base (looked for $methylation_file), skipping."
    fi
done

# plot barplot of number of SNV in methylation data
Rscript -e '
library(ggplot2)
# 1) Files with SNVs-in-methylation counts
overlap_files <- list.files("./snv/overlaps/snv_methylation_overlap",pattern = "_snv_in_methylation.txt$",full.names = TRUE)

samples_with_methylation <- basename(overlap_files)
samples_with_methylation <- sub("_snv_in_methylation.txt$", "", samples_with_methylation)
samples_with_methylation <- sub("^SNV_", "", samples_with_methylation)  # <-- IMPORTANT

# 2) Total SNVs for those samples only
snv_counts <- read.table("./snv/snv_counts/snv_per_sample.tsv",header = FALSE, stringsAsFactors = FALSE)
colnames(snv_counts) <- c("idx","sample","count")
# subset to samples with methylation data
snv_sub <- snv_counts[snv_counts$sample %in% samples_with_methylation, ]
total_snv <- sum(snv_sub$count)
cat("Total SNVs in samples with methylation:", total_snv, "\n")

# 3) Sum SNVs in methylation probes : for each overlap file, read the count and sum them
counts <- vapply(
  overlap_files,
  function(f) {as.numeric(readLines(f, n = 1))}, # read a single line and convert to numeric
  numeric(1)) # output type numeric of length 1
snv_in_methylation <- sum(counts) 
cat("Total SNVs in methylation probes:", snv_in_methylation, "\n")

# 4) Barplot in vs out
df <- data.frame(category = c("SNVs in methylation probes", "SNVs outside methylation probes"),count = c(snv_in_methylation, total_snv - snv_in_methylation))
df$percent <- round(100 * df$count / sum(df$count), 2)
p <- ggplot(df, aes(x = category, y = count, fill = category)) +
  geom_col() +
  geom_text(aes(label = paste0(count, " (", percent, "%)")),  vjust = -0.3) +
  theme_minimal() +
  labs(title = "Proportion of SNVs overlapping Methylation probes",x = "", y = "Number of SNVs") +
  theme(axis.text.x = element_text(angle = 15, hjust = 1),legend.position = "none")
ggsave("./results/summarys/SNV_in_vs_out_Methylation_barplot.pdf",p, width = 7, height = 5, bg = "white")
'
# 18) Overlap SNV in TF motifs and ATAC peaks 
mkdir -p ./snv/overlaps/intersected_motifs_snv/
motifs_in_peaks_NRF1=./motifs/overlaps/motif_peak_overlaps/NRF1_filtered_6mer.MA0506.3_peak_overlaps.bed
motifs_in_peaks_BANP=./motifs/overlaps/motif_peak_overlaps/ZBTB33_filtered_6mer.MA0527.2_peak_overlaps.bed
bedtools intersect -u -a ./snv/all_unique_SNVs_across_cancers.bed -b $motifs_in_peaks_NRF1 > ./snv/overlaps/intersected_motifs_snv/NRF1_SNVs_in_motifs_in_peaks.bed
bedtools intersect -u -a ./snv/all_unique_SNVs_across_cancers.bed -b $motifs_in_peaks_BANP > ./snv/overlaps/intersected_motifs_snv/BANP_SNVs_in_motifs_in_peaks.bed

# create a barplot : motifs , peaks and snv for TF : the % is calculated based on the total number of motifs of each TF
Rscript -e '
library(ggplot2)

## 1) Read motif sets and build unique IDs (chr:start:end)
motif_all   <- read.table("./motifs/ZBTB33_filtered_6mer.MA0527.2.bed")
motif_ids   <- unique(with(motif_all, paste(V1, V2, V3, sep=":")))

motif_peak  <- read.table("./motifs/overlaps/motif_peak_overlaps/ZBTB33_filtered_6mer.MA0527.2_peak_overlaps.bed")
peak_ids    <- unique(with(motif_peak, paste(V1, V2, V3, sep=":")))

motif_snv <- read.table("./snv/overlaps/intersected_motifs_snv/BANP_intersected_SNVs_motifskept.bed")
snv_ids   <- unique(with(motif_snv, paste(V1, V2, V3, sep=":")))

## 2) Counts
total      <- length(motif_ids)
in_atac    <- length(peak_ids)
in_variant <- length(snv_ids)
in_both    <- length(intersect(peak_ids, snv_ids))

## 3) Build data frame with percentages + pretty labels
df <- data.frame(
  category = c("All BANP motifs",
               "Motifs with variants",
               "Motifs in ATAC peaks",
               "Motifs in ATAC peaks with variants"),
  count = c(total, in_variant, in_atac, in_both))

# make sure factor levels MATCH the labels above and are in a nice order
df$category <- factor(
  df$category,
  levels = c("All BANP motifs",
             "Motifs with variants",
             "Motifs in ATAC peaks",
             "Motifs in ATAC peaks with variants"))

df$percent <- df$count / total * 100
df$label   <- sprintf("%d (%.2f%%)", df$count, df$percent)

## 4) Plot
p <- ggplot(df, aes(x = category, y = count, fill = category)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = label), vjust = -0.4, size = 4) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 20, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  ) +
  labs(
    title = "BANP 6-mer motif distribution across datasets",
    subtitle = sprintf("All BANP 6-mer motifs (n = %d)", total),
    y = "Number of motifs",
    x = NULL)

print(p)

ggsave("./results/summarys/BANP_3D_vs_peaks_vs_SNV_barplot.pdf",plot = p, width = 10, height = 8, dpi = 300, bg = "white")
'

# One Barplot too have motifs , peaks and snv for TF : the % is calculated based on the total number of motifs of each TF
Rscript -e '
library(ggplot2)

## 1) Read motif sets and build unique IDs (chr:start:end)
motif_all   <- read.table("./motifs/ZBTB33_filtered_6mer.MA0527.2.bed")
motif_ids   <- unique(with(motif_all, paste(V1, V2, V3, sep=":")))

motif_peak  <- read.table("./motifs/overlaps/motif_peak_overlaps/ZBTB33_filtered_6mer.MA0527.2_peak_overlaps.bed")
peak_ids    <- unique(with(motif_peak, paste(V1, V2, V3, sep=":")))

motif_hm450 <- read.table("./motifs/overlaps/intersected_motifs_HM450/ZBTB33_intersected_methylation.bed")
hm450_ids   <- unique(with(motif_hm450, paste(V1, V2, V3, sep=":")))

motif_snv <- read.table("./snv/overlaps/intersected_motifs_snv/BANP_intersected_SNVs_motifskept.bed")
snv_ids   <- unique(with(motif_snv, paste(V1, V2, V3, sep=":")))

## 2) Counts
total          <- length(motif_ids)
in_atac        <- length(peak_ids)
in_hm450       <- length(hm450_ids)
in_snv         <- length(snv_ids)
in_peak_hm450  <- length(intersect(peak_ids, hm450_ids))
in_peak_snv    <- length(intersect(peak_ids, snv_ids))

## 3) Build data frame with percentages + pretty labels
df <- data.frame(
  category = c("All BANP motifs",
               "Motifs with variants",
               "Motifs in ATAC peaks",
               "Motifs in ATAC peaks with variants",
               "Motifs with HM450 probes",
               "Motifs in ATAC peaks and HM450"),
  count = c(total,
            in_snv,
            in_atac,
            in_peak_snv,
            in_hm450,
            in_peak_hm450))

df$category <- factor(
  df$category,
  levels = c("All BANP motifs",
             "Motifs with variants",
             "Motifs in ATAC peaks",
             "Motifs in ATAC peaks with variants",
             "Motifs with HM450 probes",
             "Motifs in ATAC peaks and HM450"))

df$percent <- df$count / total * 100
df$label   <- sprintf("%d (%.2f%%)", df$count, df$percent)

## 4) Plot
p <- ggplot(df, aes(x = category, y = count, fill = category)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = label), vjust = -0.4, size = 4) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 20, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  ) +
  labs(
    title = "BANP 6-mer motif distribution across datasets",
    subtitle = sprintf("All BANP 6-mer motifs (n = %d)", total),
    y = "Number of motifs",
    x = NULL)

print(p)

ggsave("./results/summarys/BANP_3D_vs_peaks_vs_HM450_vs_SNV_barplot.pdf",plot = p, width = 10, height = 8, dpi = 300, bg = "white")
'
# one barplot for both TFs NRF1 and BANP together : this script creates a small barplot showing the percentage of motifs in each category for both TFs 
Rscript -e '
library(ggplot2)

# Function to compute stats
compute_tf_stats <- function(tf_name, motif_file, peak_file, hm450_file, snv_file) {

  motif_all <- read.table(motif_file)
  motif_ids <- unique(with(motif_all, paste(V1, V2, V3, sep=":")))
  total <- length(motif_ids)

  motif_peak <- read.table(peak_file)
  peak_ids   <- unique(with(motif_peak, paste(V1, V2, V3, sep=":")))
  in_atac    <- length(peak_ids)

  motif_hm450 <- read.table(hm450_file)
  hm450_ids   <- unique(with(motif_hm450, paste(V1, V2, V3, sep=":")))
  in_hm450    <- length(hm450_ids)

  motif_snv <- read.table(snv_file)
  snv_ids   <- unique(with(motif_snv, paste(V1, V2, V3, sep=":")))
  in_snv    <- length(snv_ids)

  in_peak_snv   <- length(intersect(peak_ids, snv_ids))
  in_peak_hm450 <- length(intersect(peak_ids, hm450_ids))

  data.frame(
    TF = tf_name,
    category = c(
      "All motifs",
      "Motifs with variants",
      "Motifs in ATAC peaks",
      "Motifs in ATAC peaks + variants",
      "Motifs with HM450 probes",
      "Motifs in ATAC peaks + HM450"
    ),
    count = c(total, in_snv, in_atac, in_peak_snv, in_hm450, in_peak_hm450),
    total = total
  )
}

# Load NRF1
df_nrf1 <- compute_tf_stats(
  "NRF1",
  "./motifs/NRF1_filtered_6mer.MA0506.3.bed",
  "./motifs/overlaps/motif_peak_overlaps/NRF1_filtered_6mer.MA0506.3_peak_overlaps.bed",
  "./motifs/overlaps/intersected_motifs_HM450/NRF1_intersected_methylation.bed",
  "./snv/overlaps/intersected_motifs_snv/NRF1_intersected_SNVs_motifskept.bed")

# Load BANP
df_banp <- compute_tf_stats(
  "BANP",
  "./motifs/ZBTB33_filtered_6mer.MA0527.2.bed",
  "./motifs/overlaps/motif_peak_overlaps/ZBTB33_filtered_6mer.MA0527.2_peak_overlaps.bed",
  "./motifs/overlaps/intersected_motifs_HM450/ZBTB33_intersected_methylation.bed",
  "./snv/overlaps/intersected_motifs_snv/BANP_intersected_SNVs_motifskept.bed")

# Merge datasets
df <- rbind(df_nrf1, df_banp)

df$category <- factor(df$category,
  levels = c("All motifs",
             "Motifs with variants",
             "Motifs in ATAC peaks",
             "Motifs in ATAC peaks + variants",
             "Motifs with HM450 probes",
             "Motifs in ATAC peaks + HM450"))

# Percent values
df$percent <- df$count / df$total * 100


# Short X-axis labels
df$category_short <- factor(
  df$category,
  labels = c("All",
             "Variants",
             "ATAC",
             "ATAC+Var",
             "HM450",
             "ATAC+HM450"),
  levels = levels(df$category))

# Compact, small plot

p <- ggplot(df, aes(x = category_short, y = percent, fill = TF)) +
  geom_col(position = position_dodge(width = 0.3), width = 0.25) +   # thin, close bars
  theme_minimal(base_size = 6) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 5),
    axis.text.y = element_text(size = 5),
    axis.title.x = element_text(size = 6),
    axis.title.y = element_text(size = 6),
    plot.title  = element_text(size = 7, hjust = 0.5),
    plot.subtitle = element_text(size = 6, hjust = 0.5),
    legend.position = "top",
    legend.text = element_text(size = 5),
    legend.title = element_text(size = 6)) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(
    title = "NRF1 and BANP motif distribution across datasets",
    subtitle = "Percentage of total motifs per TF",
    x = "",
    y = "Percentage of motifs (%)",
    fill = "TF")

# Save figure
ggsave("./results/summarys/NRF1_BANP_combined_barplot_PERCENT_ULTRASMALL.pdf",p,width = 4,height = 2.5,dpi = 300,bg = "white")
'

# 19) Heatmap to show SNVs across differents cancer types 

# 1)Build a list of unique SNVs per cancer type
mkdir -p ./snv/snv_unique_per_cancer/
# extract cancer types from filenames
cancers=$(find ./snv -maxdepth 1 -name 'SNV_TCGA-*.vcf.gz' | sed 's#.*/SNV_TCGA-##; s/-.*//' | sort -u)

for cancer in $cancers; do
    echo "→ Processing $cancer"
    zgrep -hv "^#" ./snv/SNV_TCGA-"$cancer"-*.vcf.gz \
        | awk '{print $1":"$2":"$4":"$5}' \
        | sort -u \
        | gzip > ./snv/snv_unique_per_cancer/${cancer}_unique_SNVs.txt.gz
    echo "Done $cancer"
done

# 2)Create a unified list of all unique SNVs across all cancer types
zcat ./snv/snv_unique_per_cancer/*_unique_SNVs.txt.gz | sort -u | gzip > ./snv/all_unique_SNVs_across_cancers.txt.gz

# 3) Build a presence/absence matrix of SNVs (rows = SNVs, columns = cancer types)
outfile="./snv/snv_presence_matrix.tsv"

# header : SNV    ACC  BLCA  BRCA  COAD ...
echo -ne "SNV" > $outfile # first column name of the header is SNV
for f in ./snv/snv_unique_per_cancer/*_unique_SNVs.txt.gz; do
    cancer=$(basename "$f" _unique_SNVs.txt.gz)
    echo -ne "\t${cancer}" >> $outfile # add a tab + the cancer name for each type of cancer
done
echo "" >> $outfile # add a newline at the end of the header
  
# Filter only SNV that overlap with motifs of NRF1 and BANP: 
# overlap SNV with motifs of NRF1 and BANP 
mkdir -p ./snv/overlaps/intersected_motifs_snv/ 
zcat ./snv/all_unique_SNVs_across_cancers.txt.gz | awk -F':' -v OFS='\t' '{print $1, $2-1, $2, $0}' | sort -k1,1 -k2,2n > ./snv/all_unique_SNVs_across_cancers.bed

# transform the SNV data file to a bed file : chrom pos-1 pos nom:chrom:pos:reference:variation intersecting with the motif bed files and keep only the SNVs that overlap with the motifs
bedtools intersect -u -a ./snv/all_unique_SNVs_across_cancers.bed -b ./motifs/NRF1_filtered_6mer.MA0506.3.bed > ./snv/overlaps/intersected_motifs_snv/NRF1_intersected_SNVs.bed 
bedtools intersect -u -a ./snv/all_unique_SNVs_across_cancers.bed -b ./motifs/ZBTB33_filtered_6mer.MA0527.2.bed > ./snv/overlaps/intersected_motifs_snv/BANP_intersected_SNVs.bed

# 4) Keep only the motifs that have SNV inside (instead of keeping the SNV that overlap with motifs) and plot the number of motifs that have SNV inside for NRF1 and BANP
# transform the SNV data file to a bed file but intersecting with the motif bed keeping the motif that have SNV inside
# NRF1_intersected_SNVs_motifskept.bed contains the motifs that have SNV inside 
bedtools intersect -u -a ./motifs/NRF1_filtered_6mer.MA0506.3.bed -b ./snv/all_unique_SNVs_across_cancers.bed > ./snv/overlaps/intersected_motifs_snv/NRF1_intersected_SNVs_motifskept.bed 
bedtools intersect -u -a ./motifs/ZBTB33_filtered_6mer.MA0527.2.bed -b ./snv/all_unique_SNVs_across_cancers.bed > ./snv/overlaps/intersected_motifs_snv/BANP_intersected_SNVs_motifskept.bed

# count number of SNVs that overlap with NRF1 and BANP motifs
cat ./snv/overlaps/intersected_motifs_snv/NRF1_intersected_SNVs_motifskept.bed  | wc -l > ./snv/overlaps/intersected_motifs_snv/NRF1_intersected_SNVs_count_kept.tsv
cat ./snv/overlaps/intersected_motifs_snv/BANP_intersected_SNVs_motifskept.bed  | wc -l > cat ./snv/overlaps/intersected_motifs_snv/BANP_intersected_SNVs_motifskept.bed  | wc -l > ./snv/overlaps/intersected_motifs_snv/BANP_intersected_SNVs_count_kept.tsv
# Plot barplot of SNV counts in NRF1 and BANP motifs
Rscript -e '
library(ggplot2)

# Motif names
motifs <- c("NRF1", "BANP")

# Read the numeric counts from the .tsv files (SNVs in motifs)
snv_counts <- c(
  as.numeric(readLines("./snv/overlaps/intersected_motifs_snv/NRF1_intersected_SNVs_count_kept.tsv")),
  as.numeric(readLines("./snv/overlaps/intersected_motifs_snv/BANP_intersected_SNVs_count_kept.tsv")))

# TOTAL motif counts
total_motifs <- c(
  843539,  
  168557)

# Build data frame
df <- data.frame(
  Motif       = motifs,
  SNV_count   = snv_counts,
  TotalMotifs = total_motifs)

# Percentage of SNVs relative to all motifs
df$Percent <- df$SNV_count / df$TotalMotifs * 100

# Label to show **inside** the bar: value + percentage
df$BarLabel <- sprintf("%d\n(%.2f%%)", df$SNV_count, df$Percent)

# X-axis label with total motifs **under** motif name
df$MotifLabel <- sprintf("%s\n(n = %d motifs)", df$Motif, df$TotalMotifs)

# Plot
p <- ggplot(df, aes(x = MotifLabel, y = SNV_count, fill = Motif)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = BarLabel),
            vjust = -0.3,
            size = 3) +
  theme_minimal() +
  labs(
    title = "NRF1 and BANP motifs that contain variants",
    y = "Number of variants in motifs",
    x = "Motif (TF) and total motif count") +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 12),
    axis.text    = element_text(size = 10),
    axis.title   = element_text(size = 11),
    legend.position = "none")

print(p)

ggsave("./results/summarys/SNV_counts_in_motifs_NRF1_BANP_kept.pdf",plot = p, width = 10, height = 12, dpi = 300, bg = "white")
'

# 20) Build presence/absence matrix of SNVs in NRF1 and BANP motifs per cancer type using these files ./snv/overlaps/intersected_motifs_snv/NRF1_intersected_SNVs.bed and ./snv/overlaps/intersected_motifs_snv/BANP_intersected_SNVs.bed
# transform the intersection of SNV in motifs bed file to a list of SNV ids only
# note that the SNV id is in the 4th column of the bed file and that it will generate a unique list of SNV ids only so less lines than the bed file
# presence and absence matrix separtly for each type of cancer type for NRF1 and BANP motifs

# 1) make a list of SNV ids that overlap with NRF1 and BANP motifs
cut -f4 ./snv/overlaps/intersected_motifs_snv/NRF1_intersected_SNVs.bed | sort -u > ./snv/overlaps/intersected_motifs_snv/NRF1_intersected_SNVs_ids.txt
snv_ids_NRF1="./snv/overlaps/intersected_motifs_snv/NRF1_intersected_SNVs_ids.txt"
# folder for intermediate columns
mkdir -p ./snv/snv_presence_columns_NRF1/

build_column_NRF1() {
    local snv_list="$1"      # ./snv/overlaps/intersected_motifs_snv/NRF1_intersected_SNVs_ids.txt
    local cancer_file="$2"   # ./snv/snv_unique_per_cancer/BRCA_unique_SNVs.txt.gz
    local out_file="$3"      # ./snv/snv_presence_columns_NRF1/BRCA.col

    local cancer
    cancer=$(basename "$cancer_file" _unique_SNVs.txt.gz)
    echo "  → building column for $cancer"

    awk -v OFS='\t' '
        NR==FNR { has[$1]=1; next }    # first file = cancer SNVs → store in hash
        {
            val = ($1 in has ? 1 : 0); # if this SNV appears in cancer → 1 else 0
            print $1, val;
        }
    ' <(zcat "$cancer_file") "$snv_list" > "$out_file"
}


# 2) loop over cancers in parallel
threads=5    # threads minimum to keep free 
max_jobs=2   # max parallel jobs
job_count=0

for f in ./snv/snv_unique_per_cancer/*_unique_SNVs.txt.gz; do
    cancer=$(basename "$f" _unique_SNVs.txt.gz)
    out="./snv/snv_presence_columns_NRF1/${cancer}.col"

    # check available CPU before launching a new job 
    CPU=$(top -b -n 2 | awk '($1~"%Cpu"){cpu=$2}END{print int(48*(100-cpu)/100)}')
    while [ "$CPU" -le "$threads" ]; do
        echo "[$cancer] attente CPU : seulement $CPU CPU libres (min = $threads)"
        sleep 60
        CPU=$(top -b -n 2 | awk '($1~"%Cpu"){cpu=$2}END{print int(48*(100-cpu)/100)}')
    done

    echo "→ Lancement job pour $cancer (CPU libres = $CPU)"
    build_column_NRF1 "$snv_ids_NRF1" "$f" "$out" &

    job_count=$((job_count + 1))
    if [ "$job_count" -ge "$max_jobs" ]; then
        # wait for current jobs to finish before launching new ones
        wait
        job_count=0
    fi
done

# final 
wait
echo "All jobs for NRF1 columns are done."

# 3) merge toutes les colonnes en une seule matrice
first_file=./snv/snv_presence_columns_NRF1/ACC.col

# 1) Extract SNV IDs from first file → this must become tmp.matrix
awk '{print $1}' "$first_file" > ./snv/tmp.matrix

# 2) Loop through all columns and add only the 0/1 column
for f in ./snv/snv_presence_columns_NRF1/*.col; do
    base=$(basename "$f" .col)

    echo "→ Adding $base"

    # Extract column 2 robustly (cut may fail if spacing inconsistent)
    awk '{print $2}' "$f" > "./snv/${base}.values"

    # paste SNV column + new cancer column
    paste ./snv/tmp.matrix "./snv/${base}.values" > ./snv/tmp2.matrix

    mv ./snv/tmp2.matrix ./snv/tmp.matrix
    rm "./snv/${base}.values"
done

# After loop → rename tmp.matrix to final matrix
mv ./snv/tmp.matrix ./snv/snv_presence_matrix_NRF1.tsv

# 3) add header
header=./snv/header.tmp
echo -ne "SNV" > "$header"

for f in ./snv/snv_presence_columns_NRF1/*.col; do
    base=$(basename "$f" .col)
    echo -ne "\t$base" >> "$header"
done
echo >> "$header"

# Combine header + matrix into final file
cat "$header" ./snv/snv_presence_matrix_NRF1.tsv > ./snv/snv_presence_matrix_NRF1.with_header.tsv

# Replace old file
mv ./snv/snv_presence_matrix_NRF1.with_header.tsv ./snv/snv_presence_matrix_NRF1.tsv
rm "$header"


# 11) heatmap of raw shared SNVs between cancer types IN A LOG10 SCALE
# it also adds the clustering groups cancers that have similar patterns of shared NRF1-motif SNVs across all cancer types
Rscript -e '
library(pheatmap)

data <- read.table("./snv/snv_presence_matrix_NRF1.tsv",header = TRUE, stringsAsFactors = FALSE)
mat <- as.matrix(data[ , -1])
rownames(mat) <- data$SNV

# shared[i,j] = number of SNVs present in both cancer i and j
shared <- t(mat) %*% mat

# ---- log10 transform to enhance contrast ----
shared_log <- log10(shared + 1)   # +1 to avoid log10(0)

# remove diagonal from the matrix we plot
diag(shared_log) <- NA

# compute color range only on off-diagonal
min_val <- min(shared_log, na.rm = TRUE)
max_val <- max(shared_log, na.rm = TRUE)
cat("Color range (log10 scale):", min_val, "to", max_val, "\n")

pdf("./results/summarys/SNV_shared_NRF1_variants_between_cancer_types_heatmap_rescaled.pdf",width = 10, height = 8, bg = "white")

pheatmap(shared_log,
         main = "Number of shared NRF1-motif variants\nbetween cancer types (log10 scale)",
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = colorRampPalette(c("white","steelblue","navy"))(200),
         breaks = seq(min_val, max_val, length.out = 201),
         na_col = "grey90")   # diagonal = grey

dev.off()
'

# 21) Heatmap with the percentage oof shared variants between the different cancer types : 
Rscript -e '
library(pheatmap)

# 1) Read presence/absence matrix
#    Rows = SNVs, columns = SNV + cancer types
data <- read.table("./snv/snv_presence_matrix_BANP.tsv",header = TRUE,stringsAsFactors = FALSE)

# 2) Sort cancer types alphabetically (exclude first column "SNV")
cancer_cols  <- colnames(data)[-1]
cancer_order <- sort(cancer_cols)

# Rebuild matrix with columns in alphabetical order
mat <- as.matrix(data[, cancer_order])
rownames(mat) <- data$SNV

# 3) Shared SNV counts between cancer types
#    shared[i,j] = number of SNVs present in BOTH cancer i and j
shared <- t(mat) %*% mat

# 4) Convert counts to PERCENT of SNVs in the row cancer type
#    pct_shared[row = A, col = B] = (shared(A,B) / total SNVs in A) * 100
totals <- diag(shared)                          # total SNVs per cancer
pct_shared <- sweep(shared, 1, totals, FUN = "/") * 100

# Remove diagonal (100%) so scale focuses on between-cancer similarity
diag(pct_shared) <- NA

# 5) Color range based on off-diagonal values only
min_val <- min(pct_shared, na.rm = TRUE)
max_val <- max(pct_shared, na.rm = TRUE)
cat("NRF1: percentage range:", min_val, "to", max_val, "\n")

pdf("./results/summarys/SNV_shared_BANP_percentage_between_cancer_types_heatmap.pdf",width = 10, height = 8, bg = "white")

pheatmap(pct_shared,
         main = "Percentage of shared BANP-motif variants\nbetween cancer types",
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = colorRampPalette(c("white","steelblue","navy"))(200),
         breaks = seq(min_val, max_val, length.out = 201),
         na_col = "grey90")

dev.off()
'
# 22) Venn diagram of motifs with peaks and methylation probes or SNVs
Rscript -e '
library(VennDiagram)
library(grid)

# 1) Load NRF1 motif sets
motif_all <- unique(with(read.table("./motifs/NRF1_filtered_6mer.MA0506.3.bed"),paste(V1, V2, V3, sep=":")))

motif_peak <- unique(with(read.table("./motifs/overlaps/motif_peak_overlaps/NRF1_filtered_6mer.MA0506.3_peak_overlaps.bed"),paste(V1, V2, V3, sep=":")))

motif_variants <- unique(with(read.table("./snv/overlaps/intersected_motifs_snv/NRF1_intersected_SNVs_motifskept.bed"),paste(V1, V2, V3, sep=":")))

# 2) 3-way Venn
venn.plot <- venn.diagram(
  x = list(
    "All NRF1 motifs"          = motif_all,
    "Motifs in ATAC peaks"     = motif_peak,
    "Motifs with variants" = motif_variants),
  filename = NULL,
  fill = c("steelblue", "lightgreen", "salmon"),
  alpha = 0.5,
  lwd = 2,
  cex = 1.3,
  cat.cex = 1.3,
  main = "NRF1 motif overlap: genomic, accessible, and variant-covered")

# 3) Save to PDF
pdf("./results/summarys/NRF1_motifs_ATAC_variant_3D_venn.pdf", width = 6, height = 6)
grid.draw(venn.plot)
dev.off()
'
# 23) Proportional Euler diagram of motifs with peaks and methylation probes

Rscript -e '
library(eulerr)

# Charger les sets de motif (6-mer)
motif_all <- unique(with(read.table("./motifs/ZBTB33_filtered_6mer.MA0527.2.bed"),paste(V1, V2, V3, sep = ":")))

motif_peak <- unique(with(read.table("./motifs/overlaps/motif_peak_overlaps/ZBTB33_filtered_6mer.MA0527.2_peak_overlaps.bed"),paste(V1, V2, V3, sep = ":")))

motif_variants <- unique(with(read.table("./snv/overlaps/intersected_motifs_snv/BANP_intersected_SNVs_motifskept.bed"),paste(V1, V2, V3, sep = ":")))
# Liste nommée de sets 
sets <- list(
  "All NRF1 motifs"    = motif_all,
  "Motifs in ATAC"     = motif_peak,
  "Motifs with variants"  = motif_variants)

# Calcul du diagramme proportionnel avec eulerr
fit <- euler(sets)

pdf("./results/summarys/BANP_motifs_ATAC_variants_3D_venn.pdf", width = 7, height = 7)
plot(
  fit,
  fills = c("steelblue", "lightgreen", "red"),  # ta palette BANP
  edges = "black",
  quantities = list(type = "counts", cex = 0.9),
  labels = list(cex = 1.1),
  main = "BANP motifs: genome-wide, accessible,\n and covered by HM450 probes")
dev.off()
'

# 24) Plots for expected results poster: 

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

# 25) Range of mutations for each cancer type : cancer_type min_mutation_length max_mutation_length
# The loop goes through each cancer type, reads all its compressed VCF files, calculates the size of each mutation in absolute value, keeps track of the smallest and largest values, and writes one summary line per cancer to the output file.

#create a table of min and max mutation ranges for each cancer type
echo -e "cancer\tmin\tmax" > ./snv/snv_min_max_cancer_type.tsv
#loop for each file and extract the min and max length of mutation found
for cancer in $(ls ./snv/SNV_TCGA-*.vcf.gz | cut -d'-' -f2 | sort -u)
do
  zcat ./snv/SNV_TCGA-${cancer}-*.vcf.gz | \
  awk 'BEGIN{FS="\t"} !/^#/ {
    v=length($5)-length($4);
    if(v<0) v=-v;
    print $1, $2, v}' | \
  # For each line of input, it looks at the third column ($3, the mutation length), and keeps track of the smallest and largest value seen so far by looking at each line and comparing the value with the already seen values.
  awk '{
    v=$3        
    if(NR==1 || v<min) min=v
    if(NR==1 || v>max) max=v
}
  END{print "'$cancer'\t" min "\t" max}' >> ./snv/snv_min_max_cancer_type.tsv
done

# dumbbell plot /snv_min_max_cancer_type.tsv to have the minimum and maximum length of variation for each cancer type
Rscript -e '
library(ggplot2)
# Read data
df <- read.table("./snv/snv_min_max_cancer_type.tsv",header = TRUE, sep = "\t")
# Plot
p <- ggplot(df, aes(y = reorder(cancer, max))) +
  geom_segment(aes(x = min, xend = max, y = cancer, yend = cancer),color = "grey70", linewidth = 0.6) +
  geom_point(aes(x = min, color = "min"), size = 2) +
  geom_point(aes(x = max, color = "max"), size = 2) +
  scale_color_manual(
    name = "Variation type",
    values = c(min = "#264653", max = "#E76F51")) +
  theme_minimal(base_size = 12) +
  labs(
    x = "Variation length (bp)",
    y = "Cancer type")
# Save plot
ggsave("./results/summarys/snv_min_max_dumbbell.pdf",plot = p, width = 7, height = 9, dpi = 300, bg = "white")'

# 26) frequency of mutations found in cancer types in SNV small indels and structural variants 
# This script loops over each TCGA cancer type, scans all corresponding VCF files, and counts SNVs (1bp), small indels (<50 bp), and large (≥50 bp) indels.
# It runs the analysis in parallel (up to two cancers at a time) and writes one summary line per cancer to a TSV file.

echo -e "cancer\tSNV\tsmall_indel\tstruct_indel" > ./snv/snv_type_freq.tsv
max_jobs=2
job_count=0
for cancer in $(ls ./snv/SNV_TCGA-*.vcf.gz | cut -d'-' -f2 | sort -u)
do
  echo "→ Launching $cancer"
  (
    zcat ./snv/SNV_TCGA-${cancer}-*.vcf.gz | \
    awk 'BEGIN{FS="\t"} !/^#/ {
      d = length($5) - length($4)
      if (d < 0) d = -d

      if (length($4)==1 && length($5)==1) snv++
      else if (d >= 50) struct++
      else small++
    }
    END {
      print "'$cancer'\t" snv "\t" small "\t" struct
    }'
  ) >> ./snv/snv_type_freq.tsv &

  job_count=$((job_count + 1))
  if [ "$job_count" -ge "$max_jobs" ]; then
    wait
    job_count=0
  fi
done
wait
echo "All jobs done."

# plot stacked barplot of mutation types per cancer type in log10 scale
Rscript -e '
library(ggplot2)
library(tidyr)

# Read data
df <- read.table("./snv/snv_type_freq.tsv", header = TRUE, sep = "\t")

# Long format
df_long <- pivot_longer(
  df,
  cols = c(SNV, small_indel, struct_indel),
  names_to = "type",
  values_to = "count"
)

# Plot
p <- ggplot(df_long, aes(x = cancer, y = count, fill = type)) +
  geom_col() +
  scale_y_log10() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    legend.position = "top"
  ) +
  labs(
    x = "Cancer type",
    y = "Number of variants (log10)",
    fill = "Variant type"
  )

# Save
ggsave(
  "./results/summarys/stacked_mutation_types_log10.pdf",
  p,
  width = 12,
  height = 7,
  dpi = 300,
  bg = "white")'

# Compare the variation distribution between samples of cancer types using boxplot for SNV, small indels, structural indels to visualize each point to see if we have outliers 
mkdir -p ./snv/snv_counts_per_sample

echo -e "sample\tcancer_type\tSNV\tsmall_indel\tstruct_indel" \
  > ./snv/snv_counts_per_sample/snv_type_per_sample.tsv

for vcf in ./snv/SNV_TCGA-*.vcf.gz; do
    sample=$(basename "$vcf" .vcf.gz | sed 's/^SNV_//')
    cancer=$(echo "$sample" | cut -d'-' -f2)

    zcat "$vcf" | awk -v sample="$sample" -v cancer="$cancer" '
      BEGIN{FS="\t"; snv=0; small=0; struct=0}
      !/^#/ {
        d = length($5) - length($4)
        if (d < 0) d = -d

        if (length($4)==1 && length($5)==1) snv++
        else if (d >= 50) struct++
        else small++
      }
      END {
        printf "%s\t%s\t%d\t%d\t%d\n", sample, cancer, snv, small, struct
      }'
done >> ./snv/snv_counts_per_sample/snv_type_per_sample.tsv

# compute the boxplot (variants in %)

Rscript -e '
library(ggplot2)
library(reshape2)

# 1) Load per-sample variant COUNTS 
df <- read.table("./snv/snv_counts_per_sample/snv_type_per_sample.tsv",header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 2) Convert counts -> percentages per sample
df$total <- df$SNV + df$small_indel + df$struct_indel
df <- df[df$total > 0, ]  # safety check to avoid division by zero (if any samples have zero variants)

df$pct_SNV          <- df$SNV / df$total * 100
df$pct_small_indel  <- df$small_indel / df$total * 100
df$pct_struct_indel <- df$struct_indel / df$total * 100

# Keep only % columns + sample + cancer_type
df_pct <- df[, c("sample", "cancer_type","pct_SNV", "pct_small_indel", "pct_struct_indel")]

# 3) Go to long format
df_long <- melt(df_pct,
                id.vars = c("sample", "cancer_type"),
                variable.name = "region",
                value.name   = "percentage")

# nicer facet labels
df_long$region <- factor(df_long$region,
                         levels = c("pct_SNV","pct_small_indel","pct_struct_indel"),
                         labels = c("SNV", "Small indel (<50bp)", "Structural indel (>=50bp)"))

# 4) Order cancer types alphabetically (for x-axis)
df_long$cancer_type <- factor(df_long$cancer_type,levels = sort(unique(df_long$cancer_type)))

# 5) Boxplot + jitter 
pdf("./results/summarys/variant_type_percentages_boxplot_per_sample.pdf",width = 16, height = 9)

ggplot(df_long, aes(x = cancer_type, y = percentage, color = region)) +
  geom_boxplot(outlier.colour = NA) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 0.7) +
  facet_wrap(~ region, nrow = 3, scales = "fixed") +
  scale_color_brewer(palette = "Set2") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 60, hjust = 1, size = 7),
    strip.text = element_text(size = 12),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    title = "SNV / small indel / structural indel percentages per sample",
    x = "Cancer type",
    y = "Percentage (%)",
    color = "Variant type")

dev.off()
'

# 27) frequency of variants found within each cancer type : number of samples that share each variant within a cancer type
mkdir -p ./snv/within_cancer_freq

for cancer in $(ls ./snv/SNV_TCGA-*.vcf.gz | cut -d'-' -f2 | sort -u); do
  nsamples=$(ls ./snv/SNV_TCGA-${cancer}-*.vcf.gz | wc -l)
  echo "→ $cancer ($nsamples samples)"

  for vcf in ./snv/SNV_TCGA-${cancer}-*.vcf.gz; do
    echo "   processing $(basename "$vcf")"
    zcat "$vcf"
  done \
  | awk 'BEGIN{FS="\t"} !/^#/ {print $1":"$2":"$4":"$5}' \
  | sort \
  | uniq -c \
  | awk -v N="$nsamples" '{print $2"\t"$1"\t"$1/N}' \
  | sort -k2,2nr \
  > ./snv/within_cancer_freq/${cancer}_variant_samplefreq.tsv

done

# 28) plot the PCA or tSNE of the samples based on the presence/absence matrix of SNVs in NRF1 and BANP motifs per cancer type
# first filter SNV that are in motifs and fall in peaks
mkdir -p ./snv/overlaps/intersected_motifs_in_peaks_snv

# NRF1: SNVs that overlap NRF1 motifs that are in peaks
bedtools intersect -u \
  -a ./snv/all_unique_SNVs_across_cancers.bed \
  -b ./motifs/overlaps/motif_peak_overlaps/NRF1_filtered_6mer.MA0506.3_peak_overlaps.bed \
  > ./snv/overlaps/intersected_motifs_in_peaks_snv/NRF1_SNVs_in_motifs_in_peaks.bed

cut -f4 ./snv/overlaps/intersected_motifs_in_peaks_snv/NRF1_SNVs_in_motifs_in_peaks.bed | sort -u \
  > ./snv/overlaps/intersected_motifs_in_peaks_snv/NRF1_SNVs_in_motifs_in_peaks_ids.txt

wc -l ./snv/overlaps/intersected_motifs_in_peaks_snv/NRF1_SNVs_in_motifs_in_peaks_ids.txt

# BANP: SNVs that overlap BANP motifs that are in peaks
bedtools intersect -u \
  -a ./snv/all_unique_SNVs_across_cancers.bed \
  -b ./motifs/overlaps/motif_peak_overlaps/ZBTB33_filtered_6mer.MA0527.2_peak_overlaps.bed \
  > ./snv/overlaps/intersected_motifs_in_peaks_snv/BANP_SNVs_in_motifs_in_peaks.bed

cut -f4 ./snv/overlaps/intersected_motifs_in_peaks_snv/BANP_SNVs_in_motifs_in_peaks.bed | sort -u \
  > ./snv/overlaps/intersected_motifs_in_peaks_snv/BANP_SNVs_in_motifs_in_peaks_ids.txt

wc -l ./snv/overlaps/intersected_motifs_in_peaks_snv/BANP_SNVs_in_motifs_in_peaks_ids.txt

# filter the ones that are in mootiifs and in peaks and are in at least 3 samples of a cancer type
ids=./snv/overlaps/intersected_motifs_in_peaks_snv/NRF1_SNVs_in_motifs_in_peaks_ids.txt
outdir=./snv/within_cancer_freq/filtered_NRF1_inMotifsInPeaks_min3
mkdir -p "$outdir"

for f in ./snv/within_cancer_freq/*_variant_samplefreq.tsv; do
  cancer=$(basename "$f" _variant_samplefreq.tsv)

  awk -v OFS="\t" '
    NR==FNR {
      keep[$1]=1
      next
    }
    {
      # split variant ID: chr:pos:REF:ALT
      n = split($1, a, ":")
      ref = a[3]
      alt = a[4]
      len = length(ref) - length(alt)
      if (len < 0) len = -len
    }
    ($1 in keep) && ($2 >= 3) && (len <= 50) {
      print
    }
  ' "$ids" "$f" \
  > "${outdir}/${cancer}_NRF1_inMotifsInPeaks_min3_lenLE50.tsv"
done


# redo Ven diagram overlaps betweenn moottifs atac peaks and variants whith filtering structural variants 

# filtering VCF files to remove structural variants (indels >= 50 bp)
mkdir ./snv/snv_filtered_without_structural_variants/
for vcf in ./snv/SNV_TCGA-*.vcf.gz; do
  out=./snv/snv_filtered_without_structural_variants/$(basename "$vcf")

  zcat "$vcf" | awk '
    BEGIN{FS="\t"; OFS="\t"}
    /^#/ {print; next}
    {
      d = length($5) - length($4)
      if (d < 0) d = -d
      if (d <= 50) print
    }
  ' | gzip > "$out"

  echo "Filtered $(basename "$vcf")"
done

# overlap filtered VCFs with motifs in peaks for NRF1

mkdir -p ./snv/overlaps_filtered/
mkdir -p ./snv/snv_filtered_unique_per_cancer/

# make one file with all unique SNVs across cancers from filtered VCFs

zgrep -hv "^#" ./snv/snv_filtered_without_structural_variants/SNV_TCGA-*.vcf.gz \
| awk '{print $1":"$2":"$4":"$5}' \
| sort -u \
| gzip > ./snv/snv_filtered_unique_per_cancer/ALL_unique_SNVs.txt.gz

#convert vcf to bed for overlap
zcat ./snv/snv_filtered_unique_per_cancer/ALL_unique_SNVs.txt.gz \
| awk -F':' '{
    chr=$1;
    pos=$2;
    start=pos-1;
    end=pos;
    id=$0;
    print chr "\t" start "\t" end "\t" id
}' \
| sort -k1,1 -k2,2n \
> ./snv/snv_filtered_unique_per_cancer/ALL_unique_SNVs.bed

# overlap with NRF1 motifs in peaks 
bedtools intersect -u \
  -a ./motifs/overlaps/motif_peak_overlaps/NRF1_filtered_6mer.MA0506.3_peak_overlaps.bed \
  -b ./snv/snv_filtered_unique_per_cancer/ALL_unique_SNVs.bed \
  > ./snv/overlaps_filtered/NRF1_filtered_SNVs_in_motifs_in_peaks.bed

# overlap with BANP motifs in peaks
bedtools intersect -u \
  -a ./motifs/overlaps/motif_peak_overlaps/ZBTB33_filtered_6mer.MA0527.2_peak_overlaps.bed \
  -b ./snv/snv_filtered_unique_per_cancer/ALL_unique_SNVs.bed \
  > ./snv/overlaps_filtered/BANP_filtered_SNVs_in_motifs_in_peaks.bed

# overlap with NRF1 motifs 
bedtools intersect -u \
  -a ./motifs/NRF1_filtered_6mer.MA0506.3.bed \
  -b ./snv/snv_filtered_unique_per_cancer/ALL_unique_SNVs.bed \
  > ./snv/overlaps_filtered/NRF1_filtered_SNVs_in_motifs.bed

# overlap with BANP motifs 
bedtools intersect -u \
  -a ./motifs/ZBTB33_filtered_6mer.MA0527.2.bed \
  -b ./snv/snv_filtered_unique_per_cancer/ALL_unique_SNVs.bed \
  > ./snv/overlaps_filtered/BANP_filtered_SNVs_in_motifs.bed

# Venn diagram 3D with filtered SNVs in motifs and in peak 
# Proportional Euler diagram of motifs with peaks and methylation probes

Rscript -e '
library(eulerr)

# Charger les sets de motif (6-mer)
motif_all <- unique(with(read.table("./motifs/ZBTB33_filtered_6mer.MA0527.2.bed"),paste(V1, V2, V3, sep = ":")))

motif_peak <- unique(with(read.table("./motifs/overlaps/motif_peak_overlaps/ZBTB33_filtered_6mer.MA0527.2_peak_overlaps.bed"),paste(V1, V2, V3, sep = ":")))

motif_variants <- unique(with(read.table("./snv/overlaps_filtered/BANP_filtered_SNVs_in_motifs.bed"),paste(V1, V2, V3, sep = ":")))
# Liste nommée de sets 
sets <- list(
  "All ZBTB33 motifs"    = motif_all,
  "Motifs in ATAC"     = motif_peak,
  "Motifs with variants"  = motif_variants)

# Calcul du diagramme proportionnel avec eulerr
fit <- euler(sets)

pdf("./results/summarys/BANP_motifs_ATAC_variants_filtered_3D_venn.pdf", width = 7, height = 7)
plot(
  fit,
  fills = c("steelblue", "lightgreen", "red"),  # ta palette BANP
  edges = "black",
  quantities = list(type = "counts", cex = 0.9),
  labels = list(cex = 1.1),
  main = "BANP motifs: genome-wide, accessible,\n and covered by HM450 probes")
dev.off()
'

# 29) Methylation distribution between cancer and healthy samples to see change in CpG methylation in NRF1 and BANP motifs
# To look into hypermethylation hypomethylation and unchanged methylation patterns in cancer compared to healthy samples at NRF1 and BANP motifs

# 1) build pairs files with sample healthy cancer pairs for each cancer type
mkdir -p ./methylation/sample_pairs_files/

echo -e "folder\ttumor_file\thealthy_file" > ./methylation/sample_pairs_files/methylation_pairs.tsv

for d in /data/papers/tcga/TCGA-*/*; do
  tumors=($d/*-0*A*_annotated_methylation.bed.gz)
  normals=($d/*-1*A*_annotated_methylation.bed.gz)

  for t in "${tumors[@]}"; do
    for n in "${normals[@]}"; do
      [ -f "$t" ] && [ -f "$n" ] && echo -e "$d\t$t\t$n"
    done
  done
done >> ./methylation/sample_pairs_files/methylation_pairs.tsv
# this file contains the list of all tumor-healthy pairs for methylation analysis for 24 cancer types that had at least one tumor-healthy pair.

# 2) compute delta methylation per pair and classify into hypo hyper unchanged and make boxplot per cancer type
Rscript -e '
library(data.table)
library(ggplot2)

pairs_file <- "./methylation/sample_pairs_files/methylation_pairs.tsv"
out_tsv    <- "./methylation/sample_pairs_files/delta_methylation_ALL_CpG_per_pair_by_cancer.tsv"
out_pdf    <- "./results/summarys/boxplot_delta_methylation_ALL_CpG_by_cancer_oneplot.pdf"

thr <- 10  # 10% (beta in [0,1])

read_probe_meth <- function(f) {
  dt <- fread(cmd = paste("zcat", shQuote(f)), header = FALSE, select = c(4,5))
  setnames(dt, c("probe","meth"))
  dt[, meth := as.numeric(meth)]
  dt
}

get_cancer <- function(folder) {
  parts <- strsplit(folder, "/")[[1]]
  parts[which(grepl("^TCGA-[A-Z]+$", parts))[1]]
}

pairs <- fread(pairs_file)[!is.na(tumor_file) & !is.na(healthy_file)]
pairs[, cancer := vapply(folder, get_cancer, character(1))]

res <- vector("list", nrow(pairs))

for (i in seq_len(nrow(pairs))) {

  tumor  <- read_probe_meth(pairs$tumor_file[i])
  normal <- read_probe_meth(pairs$healthy_file[i])

  # safety check for probe order
  if (nrow(tumor) != nrow(normal) || !identical(tumor$probe, normal$probe)) {
    stop("Probe order mismatch at pair ", i,
         "\nTumor: ", pairs$tumor_file[i],
         "\nNormal: ", pairs$healthy_file[i])
  }

  delta <- tumor$meth - normal$meth

  hypo  <- delta < -thr
  hyper <- delta >  thr
  same  <- !hypo & !hyper

  n <- sum(!is.na(delta))

  res[[i]] <- data.table(
    cancer = pairs$cancer[i],
    folder = pairs$folder[i],
    tumor_file = pairs$tumor_file[i],
    healthy_file = pairs$healthy_file[i],
    n_cpg = n,
    pct_hypo  = round(100 * sum(hypo,  na.rm=TRUE) / n, 2),
    pct_hyper = round(100 * sum(hyper, na.rm=TRUE) / n, 2),
    pct_same  = round(100 * sum(same,  na.rm=TRUE) / n, 2)
  )

  if (i %% 50 == 0) cat("Processed", i, "of", nrow(pairs), "\n")
}

res_dt <- rbindlist(res)
dir.create(dirname(out_tsv), recursive=TRUE, showWarnings=FALSE)
fwrite(res_dt, out_tsv, sep="\t")

# --------- Make plot in the SAME STYLE as your UMR/LMR/FMR code ---------
# add counts in x labels
counts <- res_dt[, .N, by = cancer]
res_dt <- merge(res_dt, counts, by="cancer", all.x=TRUE)
res_dt[, cancer_label := paste0(cancer, " (n=", N, ")")]

# alphabetical ordering like your other plot
res_dt[, cancer_label := factor(cancer_label, levels = sort(unique(cancer_label)))]

# long format like reshape2::melt
long <- melt(
  res_dt,
  id.vars = c("cancer","cancer_label","folder","tumor_file","healthy_file","n_cpg","N"),
  measure.vars = c("pct_hypo","pct_hyper","pct_same"),
  variable.name = "class",
  value.name = "percentage"
)

long[, class := fifelse(class=="pct_hypo","Hypomethylated",fifelse(class=="pct_hyper","Hypermethylated","Unchanged"))]

dir.create(dirname(out_pdf), recursive=TRUE, showWarnings=FALSE)

pdf(out_pdf, width=16, height=9)

ggplot(long, aes(x = cancer_label, y = percentage, color = class)) +
  geom_boxplot(outlier.colour = NA) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 0.7) +
  facet_wrap(~ class, nrow = 3, scales = "fixed") +
  scale_color_brewer(palette = "Set2") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 60, hjust = 1, size = 7),
    strip.text = element_text(size = 12),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    title = "Tumor vs Healthy methylation change categories per pair",
    x = "Cancer type",
    y = "Percentage (%)",
    color = "Class"
  )

dev.off()

cat("Wrote:", out_tsv, "\n")
cat("Wrote:", out_pdf, "\n")
'

# Do the boxplot but on only CpGs that are in NRF1 motifs and BANP motifs 
mkdir -p ./methylation/overlaps/intersected_motifs_HM450/
bedtools intersect -u -a ./methylation/annotated_methylation_data_probes.bed -b ./motifs/NRF1_filtered_6mer.MA0506.3.bed > ./methylation/overlaps/intersected_motifs_HM450/NRF1_intersected_methylation.bed
bedtools intersect -u -a ./methylation/annotated_methylation_data_probes.bed -b ./motifs/ZBTB33_filtered_6mer.MA0527.2.bed  > ./methylation/overlaps/intersected_motifs_HM450/BANP_intersected_methylation.bed

Rscript -e '
library(data.table)
library(ggplot2)

pairs_file <- "./methylation/sample_pairs_files/methylation_pairs.tsv"
out_tsv    <- "./methylation/sample_pairs_files/delta_methylation_MOTIF_CpG_per_pair_by_cancer.tsv"
out_pdf    <- "./results/summarys/boxplot_delta_methylation_MOTIF_CpG_by_cancer_oneplot.pdf"

thr <- 10  # 10% (your values are in 0–100)

# =========================
# ADD THIS: choose motif CpGs file (NRF1 or BANP)
# =========================
motif_bed <- "./methylation/overlaps/intersected_motifs_HM450/NRF1_intersected_methylation.bed"
# motif_bed <- "./methylation/overlaps/intersected_motifs_HM450/BANP_intersected_methylation.bed"

motif_probes <- fread(motif_bed, header = FALSE)[[4]]
motif_probes <- unique(as.character(motif_probes))
if (length(motif_probes) == 0) stop("Motif probe list is empty: ", motif_bed)
cat("Keeping", length(motif_probes), "CpGs from:", motif_bed, "\n")

read_probe_meth <- function(f) {
  dt <- fread(cmd = paste("zcat", shQuote(f)), header = FALSE, select = c(4,5))
  setnames(dt, c("probe","meth"))
  dt[, meth := as.numeric(meth)]
  dt
}

get_cancer <- function(folder) {
  parts <- strsplit(folder, "/")[[1]]
  parts[which(grepl("^TCGA-[A-Z]+$", parts))[1]]
}

pairs <- fread(pairs_file)[!is.na(tumor_file) & !is.na(healthy_file)]
pairs[, cancer := vapply(folder, get_cancer, character(1))]

res <- vector("list", nrow(pairs))

for (i in seq_len(nrow(pairs))) {

  tumor  <- read_probe_meth(pairs$tumor_file[i])
  normal <- read_probe_meth(pairs$healthy_file[i])

  # safety check for probe order
  if (nrow(tumor) != nrow(normal) || !identical(tumor$probe, normal$probe)) {
    stop("Probe order mismatch at pair ", i,
         "\nTumor: ", pairs$tumor_file[i],
         "\nNormal: ", pairs$healthy_file[i])
  }

  # =========================
  # ADD THIS: filter to motif CpGs
  # =========================
  idx <- tumor$probe %chin% motif_probes
  delta <- tumor$meth[idx] - normal$meth[idx]

  hypo  <- delta < -thr
  hyper <- delta >  thr
  same  <- !hypo & !hyper

  n <- sum(!is.na(delta))

  res[[i]] <- data.table(
    cancer = pairs$cancer[i],
    folder = pairs$folder[i],
    tumor_file = pairs$tumor_file[i],
    healthy_file = pairs$healthy_file[i],
    n_cpg = n,
    pct_hypo  = round(100 * sum(hypo,  na.rm=TRUE) / n, 2),
    pct_hyper = round(100 * sum(hyper, na.rm=TRUE) / n, 2),
    pct_same  = round(100 * sum(same,  na.rm=TRUE) / n, 2)
  )

  if (i %% 50 == 0) cat("Processed", i, "of", nrow(pairs), "\n")
}

res_dt <- rbindlist(res)
dir.create(dirname(out_tsv), recursive=TRUE, showWarnings=FALSE)
fwrite(res_dt, out_tsv, sep="\t")

# --------- Make plot in the SAME STYLE as your UMR/LMR/FMR code ---------
counts <- res_dt[, .N, by = cancer]
res_dt <- merge(res_dt, counts, by="cancer", all.x=TRUE)
res_dt[, cancer_label := paste0(cancer, " (n=", N, ")")]

res_dt[, cancer_label := factor(cancer_label, levels = sort(unique(cancer_label)))]

long <- melt(
  res_dt,
  id.vars = c("cancer","cancer_label","folder","tumor_file","healthy_file","n_cpg","N"),
  measure.vars = c("pct_hypo","pct_hyper","pct_same"),
  variable.name = "class",
  value.name = "percentage"
)

long[, class := fifelse(class=="pct_hypo","Hypomethylated",
                 fifelse(class=="pct_hyper","Hypermethylated","Unchanged"))]

dir.create(dirname(out_pdf), recursive=TRUE, showWarnings=FALSE)

pdf(out_pdf, width=16, height=9)

ggplot(long, aes(x = cancer_label, y = percentage, color = class)) +
  geom_boxplot(outlier.colour = NA) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 0.7) +
  facet_wrap(~ class, nrow = 3, scales = "fixed") +
  scale_color_brewer(palette = "Set2") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 60, hjust = 1, size = 7),
    strip.text = element_text(size = 12),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    title = paste0("Tumor vs Healthy methylation change categories per pair (motif CpGs)\n", basename(motif_bed)),
    x = "Cancer type",
    y = "Percentage (%)",
    color = "Class"
  )

dev.off()

cat("Wrote:", out_tsv, "\n")
cat("Wrote:", out_pdf, "\n")
'

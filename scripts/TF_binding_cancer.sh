#!/bin/bash
cd TF_binding_cancer

# 1) Filter motifs with p-value threshold and sort them.
mkdir -p ./motifs

zcat /data/genome/motifs/jaspar_2024/fimo/vertebrata/hg38/NRF1.MA0506.3.bed.gz | awk '($5<=1/4^6)' > ./motifs/NRF1_filtered_6mer.MA0506.3.bed
zcat /data/genome/motifs/jaspar_2024/fimo/vertebrata/hg38/ZBTB33.MA0527.2.bed.gz | awk '($5<=1/4^6)' > ./motifs/ZBTB33_filtered_6mer.MA0527.2.bed

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
#conda install -c conda-forge r-venndiagram r-eulerr r-ggplot2 r-dplyr r-tidyr --yes

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

# Save PNG
ggsave("./results/summarys/pie_chart_cancer_peaks_count.png",plot = p,width = 8,height = 5,dpi = 300,bg = "white")
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

#to find samples with healthy and cancer data : 
for dir in /data/papers/tcga/TCGA-ESCA/*/; do
    count=$(ls "$dir"/HM450*.bed.gz 2>/dev/null | wc -l)
    if [ "$count" -gt 1 ]; then
        echo "$dir has $count HM450 methylation files"
        ls "$dir"/HM450*.bed.gz
        echo ""
    fi
done

#example of outputs from BRCA : 
#/data/papers/tcga/TCGA-BRCA/TCGA-A7-A0D9/ has 2 bed.gz files : 
# /data/papers/tcga/TCGA-BRCA/TCGA-A7-A0D9/HM450_TCGA-BRCA-TCGA-A7-A0D9-11A_1_annotated_methylation.bed.gz --> healthy
# /data/papers/tcga/TCGA-BRCA/TCGA-A7-A0D9/HM450_TCGA-BRCA-TCGA-A7-A0D9-01A_1_annotated_methylation.bed.gz --> cancer
# /data/papers/tcga/TCGA-ESCA/TCGA-IC-A6RE/ has 2 HM450 methylation files
# /data/papers/tcga/TCGA-ESCA/TCGA-IC-A6RE//HM450_TCGA-ESCA-TCGA-IC-A6RE-01A_1_annotated_methylation.bed.gz
# /data/papers/tcga/TCGA-ESCA/TCGA-IC-A6RE//HM450_TCGA-ESCA-TCGA-IC-A6RE-11A_1_annotated_methylation.bed.gz

# plot the methylation distribution for both healthy and cancer samples of TCGA-BRCA-TCGA-A7-A0D9 using a delta methylation histogram :
Rscript -e '
library(ggplot2)
# Load files (same format: chr start end probe meth)
cancer  <- read.table("/data/papers/tcga/TCGA-ESCA/TCGA-IC-A6RE//HM450_TCGA-ESCA-TCGA-IC-A6RE-01A_1_annotated_methylation.bed.gz",  header=FALSE)
healthy <- read.table("/data/papers/tcga/TCGA-ESCA/TCGA-IC-A6RE//HM450_TCGA-ESCA-TCGA-IC-A6RE-11A_1_annotated_methylation.bed.gz", header=FALSE)

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

png("./results/summarys/delta_methylation_hist_TCGA-ESCA-TCGA-IC-A6RE.png", width=1000, height=800)

hist(delta,
     breaks = 100,
     col = "steelblue",
     border = "black",
     main = paste0("Δ Methylation (Cancer – Healthy TCGA-ESCA-TCGA-IC-A6RE )\n",
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

# 14) statistics of SNV data :
mkdir -p ./snv

# Linking the SNV data from /data/papers/tcga to the snv folder  
for file in /data/papers/tcga/TCGA-*/*/SNV_TCGA*vcf.gz; do
    ln -s "$file" ./snv/
done

# number of samples with SNV data --> 11204 samples
ls ./snv/SNV_TCGA*vcf.gz | wc -l 

# number of SNVs per sample
mkdir -p ./snv/snv_counts/
counter=0
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

###### IMPORTANT !!!!!
# number of SNVs per cancer type
mkdir -p ./snv/snv_counts_per_cancer_type/
# unique SNVs per cancer type
cancers=$(ls ./snv/SNV_TCGA-*.vcf.gz \
          | sed 's#./snv/SNV_TCGA-##; s/-.*//' \
          | sort -u)
for cancer in $cancers; do
    count=$(zgrep -hv "^#" ./snv/SNV_TCGA-"$cancer"-*.vcf.gz \
            | awk '{print $1":"$2":"$4":"$5}' \
            | sort -u \
            | wc -l)
    echo -e "${cancer}\t${count}" >> ./snv/snv_counts_per_cancer_type/snv_unique_per_cancer.tsv
done

cut -f2,3 ./snv/snv_counts/snv_per_sample.tsv | awk -F'\t' '{
    sample=$1;           # column 1 = sample name
    snv=$2;              # column 2 = SNV count
    split(sample,a,"-"); # split sample name using '-'
    cancer=a[2];         # the 2nd element is the cancer type
    print cancer "\t" snv;}' \
    | awk '{sum[$1] += $2} END {for (c in sum) print c "\t" sum[c]}' | sort -r -k1,1  > ./snv/snv_counts_per_cancer_type/snv_total_per_cancer.tsv

# Plot histogram of SNVs per cancer type
Rscript -e '

data <- read.table("./snv/snv_counts_per_cancer_type/snv_total_per_cancer.tsv")
                   
cancer_names <- data[,1]
snv_counts   <- data[,2]
snv_millions <- snv_counts / 1e6
png("./results/summarys/snv_per_cancer_type_barplot.png", width=1200, height=900 , bg = "white")

barplot(
  snv_millions,
  names.arg=cancer_names,
  horiz=TRUE,                          # horizontal bars
  las=1,                               # labels readable
  col = rainbow(length(snv_counts)),
  border="black",
  main="SNVs per Cancer Type",
  cex.main = 3,
  cex.lab = 1.5,
  xlab="Total SNVs (millions)")

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

# heatmap of cancer types with ATAC vs SNV
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
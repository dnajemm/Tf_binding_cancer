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

Rscript -e 'library(ggplot2);
motif_files <- list.files("motifs", pattern = "_filtered_6mer.*.bed$", full.names = TRUE);
samples <- sub("_filtered_6mer.*.bed$", "", basename(motif_files));
pdf("./results/summarys/dist_tss_motifs_hist.pdf");
par(bg = "white");
for (sample in samples) {
  dist_file <- paste0("./motifs/dist_to_tss/", sample, "_dist_tss.txt");
  distances <- read.table(dist_file)[,1];
  hist(log10(distances + 1),
       xlim = c(0, 8),
       xlab = "Distance to TSS (log10)",
       main = sample,
       breaks = 50,
       col = "steelblue",
       border = "black");};
dev.off();'

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

# 4) Overlap filtered motifs with annotated methylation data to find intersecting regions and 
# see if the motifs fall within methylated probes regions.

mkdir -p ./motifs/overlaps/intersected_motifs_HM450

# transform the HM450 manifest annotated methylation data to bed format. chr starrt end probe_id
zcat ./methylation/HM450.hg38.manifest.tsv.gz | awk -vOFS="\t" ' NR>1 && $1 != "NA" && $2 != "NA" && $3 != "NA" {print $1, $2, $3, $9}' | sort -k1,1 -k2,2n  > ./methylation/annotated_methylation_data_probes.bed

# intersect NRF1 filtered motifs with annotated methylation data
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

# Count total number of filtered motifs
n_motifs_NRF1=$(wc -l < ./motifs/NRF1_filtered_6mer.MA0506.3.bed)

# Count total number of peaks 
n_probes=$(wc -l < ./methylation/annotated_methylation_data_probes.bed)

# Count number of filtered NRF1 motifs that fall within peaks
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
  category   = c("Filtered NRF1 motifs", "Methylation probes"),
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
  "Overlap between NRF1 filtered motifs and Methylation probes",
  x = 0.5, y = 0.96,
  gp = gpar(fontsize = 36,fontfamily = "Helvetica"))

pushViewport(viewport(y = 0.45))
grid.draw(venn.obj)
dev.off()
'

# Create a Venn diagram to visualize the overlaps between motifs and methylation for BANP 

# Count total number of filtered motifs
n_motifs_BANP=$(wc -l < ./motifs/ZBTB33_filtered_6mer.MA0527.2.bed)

# Count total number of probes 
n_probes=$(wc -l < ./methylation/annotated_methylation_data_probes.bed)

# Count number of filtered BANP motifs that fall within methylation probes
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
grid.text("Overlap between BANP filtered motifs and methylation probes",
          x = 0.5, y = 0.95,
          gp = gpar(fontsize = 36,fontfamily = "Helvetica"))
# place the category labels above the circles
grid.text("Methylation probes",
          x = 0.37, y = 0.73,           
          gp = gpar(fontsize = 28,fontfamily = "Helvetica"))
grid.text("Filtered BANP motifs",
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

# 7) Overlap the filtered motifs with the peaks to see which motifs fall within the peaks regions.
mkdir -p ./motifs/overlaps/motif_peak_overlaps

# Loop through each filtered motif file and intersect with all peak files from peaks
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

# Count total number of filtered motifs
n_motifs_NRF1=$(wc -l < ./motifs/NRF1_filtered_6mer.MA0506.3.bed)

# Count total number of peaks 
n_peaks=$(wc -l < ./peaks/merged_peaks.bed)

# Count number of filtered NRF1 motifs that fall within peaks
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
  category   = c("Filtered NRF1 motifs", "ATAC peaks"),
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
  "Overlap between NRF1 filtered motifs and ATAC peaks",
  x = 0.5, y = 0.96,
  gp = gpar(fontsize = 36,fontfamily = "Helvetica"))

pushViewport(viewport(y = 0.45))
grid.draw(venn.obj)
dev.off()
'

# Create a Venn diagram to visualize the overlaps between motifs and peaks for BANP 

# Count total number of filtered motifs
n_motifs_BANP=$(wc -l < ./motifs/ZBTB33_filtered_6mer.MA0527.2.bed)

# Count total number of peaks (e.g., merged peaks)
n_peaks=$(wc -l < ./peaks/merged_peaks.bed)

# Count number of filtered BANP motifs that fall within peaks
n_overlap_BANP=$(wc -l < ./motifs/overlaps/motif_peak_overlaps/ZBTB33_filtered_6mer.MA0527.2_peak_overlaps.bed)

Rscript -e '
library(VennDiagram)
library(grid)

motifs  <- '"$n_motifs_BANP"'
peaks   <- '"$n_peaks"'
overlap <- '"$n_overlap_BANP"'

# Generate venn WITHOUT category labels
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
grid.text("Overlap between BANP filtered motifs and ATAC peaks",
          x = 0.5, y = 0.95,
          gp = gpar(fontsize = 36,fontfamily = "Helvetica"))
# place the category labels above the circles
grid.text("ATAC peaks",
          x = 0.37, y = 0.82,           
          gp = gpar(fontsize = 28,fontfamily = "Helvetica"))
grid.text("Filtered BANP motifs",
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
  main = "NRF1 motifs: ATAC peaks vs HM450",
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
  main = "BANP motifs: ATAC peaks vs HM450",
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
pdf("./results/summarys/dist_tss_peaks_hist.pdf");
par(bg = "white");
for (sample in samples) {
  dist_file <- paste0("./peaks/dist_to_tss/", sample, "_dist_tss.txt");
  distances <- read.table(dist_file)[,1];
  hist(log10(distances + 1),
       xlim = c(0, 8),
       xlab = "Distance to TSS (log10)",
       main = sample,
       breaks = 50,
       col = "steelblue",
       border = "black");};
dev.off();'

# plot the median distance to TSS across all samples
Rscript -e '
library(ggplot2)
peak_files <- list.files("peaks", pattern = "_peaks_macs.bed$", full.names = TRUE);
samples <- sub("_peaks_macs.bed$", "", basename(peak_files));
median_distances <- data.frame(Sample = character(), Median_Distance = numeric())
for (sample in samples) {
  dist_file <- paste0("./peaks/dist_to_tss/", sample, "_dist_tss.txt");
  distances <- read.table(dist_file)[,1];
  median_distance <- median(distances);
  median_distances <- rbind(median_distances, data.frame(Sample = sample, Median_Distance = median_distance));
}
pdf("./results/summarys/median_dist_tss_peaks.pdf", width=12, height=18)
ggplot(median_distances, aes(x = Sample, y = Median_Distance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1 , size =12), axis.text.y = element_text(size = 4, margin = margin(r = 12))) +
  labs(title = "Median Distance to TSS for ATAC Peaks per Sample", x = "Sample", y = "Median Distance to TSS (bp)")
dev.off()
' 


# 13) statistics methylation probes per sample : 

# Plot histogram of methylation percentage for one sample

# /data/papers/tcga/TCGA-BRCA/TCGA-AR-A0TP/HM450_TCGA-BRCA-TCGA-AR-A0TP-01A_1_annotated_methylation.bed.gz
# /data/papers/tcga/TCGA-ACC/TCGA-OR-A5J1/HM450_TCGA-ACC-TCGA-OR-A5J1-01A_1_annotated_methylation.bed.gz
# /data/papers/tcga/TCGA-CHOL/TCGA-3X-AAV9/HM450_TCGA-CHOL-TCGA-3X-AAV9-01A_1_annotated_methylation.bed.gz
# /data/papers/tcga/TCGA-GBM/TCGA-OX-A56R/HM450_TCGA-GBM-TCGA-OX-A56R-01A_1_annotated_methylation.bed.gz
Rscript -e '
data <- read.table("/data/papers/tcga/TCGA-GBM/TCGA-OX-A56R/HM450_TCGA-GBM-TCGA-OX-A56R-01A_1_annotated_methylation.bed.gz")
meth <- data[, ncol(data)]   # last column = methylation %

png("./results/summarys/Methylation_distribution_TCGA-GBM-TCGA-OX-A56R-01A_1.png", width=1000, height=800)
hist(meth,
     breaks=50,
     col="steelblue",
     border="black",
     main="Methylation distribution of TCGA-GBM-TCGA-OX-A56R-01A_1 sample",
     xlab="Methylation (%)")
dev.off()
'

# Distance of methylation sites to TSS : generate a file per sample with the distances of each methylation to the nearest TSS

mkdir -p ./methylation/dist_to_tss/
file='/data/papers/tcga/TCGA-GBM/TCGA-OX-A56R/HM450_TCGA-GBM-TCGA-OX-A56R-01A_1_annotated_methylation.bed.gz'
closestBed -t first -d -a "$file" -b /data/genome/annotations/hg38_tss.bed | awk '{print $NF}' > "./methylation/dist_to_tss/TCGA-GBM-TCGA-OX-A56R-01A_1_dist_tss.txt"

# Plot histogram of distances to TSS for one sample 

Rscript -e '
dist_file <- "./methylation/dist_to_tss/TCGA-GBM-TCGA-OX-A56R-01A_1_dist_tss.txt"
distances <- read.table(dist_file)[,1]

png("./results/summarys/dist_tss_methylation_TCGA-GBM-TCGA-OX-A56R-01A_1_hist.png", width=1000, height=800)
par(bg = "white")

hist(log10(distances + 1),
     xlim = c(0, 8),
     xlab = "Distance to TSS (log10)",
     main = "TCGA-GBM-TCGA-OX-A56R-01A_1",
     breaks = 50,
     col = "steelblue",
     border = "black")
dev.off()
'

# number oof methylation probes per cancer type
mkdir -p ./methylation/methylation_counts/
for file in /data/papers/tcga/TCGA*/*/*_annotated_methylation.bed.gz; do
    sample=$(basename "$file" _annotated_methylation.bed.gz)
    count=$(zcat "$file" | wc -l)
    echo -e "$sample\t$count" >> ./methylation/methylation_counts/methylation_probes_per_sample.tsv
done

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

# number of SNVs per cancer type
mkdir -p ./snv/snv_counts_per_cancer_type/
cut -f2,3 ./snv/snv_counts/snv_per_sample.tsv | awk -F'\t' '{
    sample=$1;           # column 1 = sample name
    snv=$2;              # column 2 = SNV count
    split(sample,a,"-"); # split sample name using '-'
    cancer=a[2];         # the 2nd element is the cancer type
    print cancer "\t" snv;}' \
    | awk '{sum[$1] += $2} END {for (c in sum) print c "\t" sum[c]}'  > ./snv/snv_counts_per_cancer_type/snv_total_per_cancer.tsv
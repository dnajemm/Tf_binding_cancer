#!/bin/bash
cd TF_binding_cancer

# 1) Filter motifs with p-value threshold and sort them.
mkdir -p ./motifs

zcat /data/genome/motifs/jaspar_2024/fimo/vertebrata/hg38/NRF1.MA0506.3.bed.gz | awk '($5<=1/4^6)' > ./motifs/NRF1_filtered_6mer.MA0506.3.bed
zcat /data/genome/motifs/jaspar_2024/fimo/vertebrata/hg38/ZBTB33.MA0527.2.bed.gz | awk '($5<=1/4^6)' > ./motifs/ZBTB33_filtered_6mer.MA0527.2.bed

mkdir -p ./methylation
# 2) Download the HM450 probes bed file :
wget -O ./methylation/HM450.hg38.manifest.tsv.gz \
https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/HM450/HM450.hg38.manifest.tsv.gz

# 3) Overlap filtered motifs with annotated methylation data to find intersecting regions and 
# see if the motifs fall within methylated probes regions.

mkdir -p ./motifs/overlaps/intersected_motifs_HM450

# transform the HM450 manifest annotated methylation data to bed format
zcat ./methylation/HM450.hg38.manifest.tsv.gz | awk -vOFS="\t" ' NR>1 && $1 != "NA" && $2 != "NA" && $3 != "NA" {print $1, $2, $3}' | sort -k1,1 -k2,2n  > ./methylation/annotated_methylation_data_probes.bed

# intersect NRF1 filtered motifs with annotated methylation data
bedtools intersect -u -a ./motifs/NRF1_filtered_6mer.MA0506.3.bed -b ./methylation/annotated_methylation_data_probes.bed > ./motifs/overlaps/intersected_motifs_HM450/NRF1_intersected_methylation.bed
bedtools intersect -u -a ./motifs/ZBTB33_filtered_6mer.MA0527.2.bed -b ./methylation/annotated_methylation_data_probes.bed > ./motifs/overlaps/intersected_motifs_HM450/ZBTB33_intersected_methylation.bed

# 4) Summarize the intersected results : 1   NRF1 30234 ZBTB33 6119
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

# 5) link the peak files from the previous analysis of /data/hichamif/pred_tf_cancer/peak/ to this folder :
mkdir -p ./peaks
for file in /data/hichamif/pred_tf_cancer/peak/ATAC_TCGA*_peaks_macs.bed; do
    ln -s "$file" ./peaks/
done

#combine all peak files into one and sort them
cat ./peaks/ATAC_TCGA*_peaks_macs.bed | sort -k1,1 -k2,2n | bedtools merge -i stdin > ./peaks/merged_peaks.bed

# 6) Overlap the filtered motifs with the peaks to see which motifs fall within the peaks regions.
mkdir -p ./motifs/overlaps/motif_peak_overlaps

# Loop through each filtered motif file and intersect with all peak files from peaks
for file in ./motifs/*_filtered_6mer*.bed; do
   motif_name=$(basename "$file" .bed)
   bedtools intersect -u -a "$file" -b ./peaks/merged_peaks.bed > ./motifs/overlaps/motif_peak_overlaps/"${motif_name}_peak_overlaps.bed"
done

# 7) Summarize the intersected results : NRF1 384997 ZBTB33 69795

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

#8) Overlap the intersected methylation motifs with the peak overlapping motifs to find common regions
mkdir -p ./motifs/overlaps/intersected_overlaps

bedtools intersect -u -a ./motifs/overlaps/intersected_motifs_HM450/NRF1_intersected_methylation.bed -b ./motifs/overlaps/motif_peak_overlaps/NRF1_filtered_6mer.MA0506.3_peak_overlaps.bed > ./motifs/overlaps/intersected_overlaps/NRF1_methylation_peak_overlap.bed  
bedtools intersect -u -a ./motifs/overlaps/intersected_motifs_HM450/ZBTB33_intersected_methylation.bed -b ./motifs/overlaps/motif_peak_overlaps/ZBTB33_filtered_6mer.MA0527.2_peak_overlaps.bed > ./motifs/overlaps/intersected_overlaps/ZBTB33_methylation_peak_overlap.bed

# 9) Summarize the intersected results : NRF1 24325 ZBTB33 4863

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

#10) Setting up R environment with required packages
mkdir -p ./conda_envs/conda_pkgs
export CONDA_PKGS_DIRS=/data/najemd/TF_binding_cancer/conda_envs/conda_pkgs
conda create --prefix ./conda_envs/TF_binding_cancer_env r-base -c conda-forge -y
conda env export --prefix ./conda_envs/TF_binding_cancer_env > ./conda_envs/TF_binding_cancer_env.yml
conda activate ./conda_envs/TF_binding_cancer_env
conda install -c conda-forge r-venndiagram r-eulerr r-ggplot2 r-dplyr r-tidyr --yes

# 11) create a Venn diagram to visualize the overlaps between the three datasets for each motif

Rscript -e '
library(VennDiagram)
library(grid)

## 1) Motifs NRF1 dans les pics ATAC
motif_peak <- read.table("./motifs/overlaps/motif_peak_overlaps/NRF1_filtered_6mer.MA0506.3_peak_overlaps.bed")
peak_ids <- with(motif_peak, paste(V1, V2, V3, sep=":"))

## 2) Motifs NRF1 qui overlappent HM450
motif_hm450 <- read.table("./motifs/overlaps/intersected_motifs_HM450/NRF1_intersected_methylation.bed")
hm450_ids <- with(motif_hm450, paste(V1, V2, V3, sep=":"))

## 3) Venn à 2 ensembles : 
##    A = motifs in ATAC peaks
##    B = motifs on HM450
nrf1_lists <- list(
  "NRF1 motifs in ATAC peaks" = peak_ids,
  "NRF1 motifs on HM450"      = hm450_ids)

venn.plot <- venn.diagram(
  x = nrf1_lists,
  category.names = c("In ATAC peaks", "On HM450"),
  filename = NULL,
  fill = c("salmon", "skyblue"),
  alpha = 0.7,
  main = "NRF1 motifs: ATAC peaks vs HM450",
  cex = 2,
  cat.cex = 2,
  main.cex = 3,
  cat.pos  = c(-20, -10), 
  cat.dist = c(0.02, 0.015))

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

## 2) Motifs ZBTB33 qui overlappent HM450 (méthylés)
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
  cat.dist = c(0.02, 0.015))

png("./results/summarys/BANP_peaks_vs_HM450_venn.png", width = 1000, height = 1000)
grid.draw(venn.plot)
dev.off()
'
rm ./VennDiagram*.log  # remove temporary files created by VennDiagram package
conda deactivate
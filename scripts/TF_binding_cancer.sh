###############################################################################
# CONNECT & ENVIRONMENT SETUP
###############################################################################
#!/bin/bash
# To connect to the server dorsal : ssh najemd@serv-bardet-dorsal.igbmc.u-strasbg.fr
cd TF_binding_cancer

# Setting up R environment with required packages
mkdir -p ./conda_envs/conda_pkgs
export CONDA_PKGS_DIRS=/data/najemd/TF_binding_cancer/conda_envs/conda_pkgs
#conda create --prefix ./conda_envs/TF_binding_cancer_env r-base -c conda-forge -y
#conda env export --prefix ./conda_envs/TF_binding_cancer_env > ./conda_envs/TF_binding_cancer_env.yml
conda activate ./conda_envs/TF_binding_cancer_env
#conda install -c conda-forge r-venndiagram r-eulerr r-ggplot2 r-dplyr r-tidyr r-hexbin r-mass r-kernsmooth r-venndir r-polychrome r-plotly r-FNN bioconda bioconductor-illuminaHumanMethylation450kanno.ilmn12.hg19 bioconductor-minfi bioconductor-sesame bioconductor-sesamedata --yes 

###############################################################################
# MOTIFS SCRIPT
###############################################################################


# Install MEM BANP of JASPAR 2025 and put it in the folder of MEME:

#wget "https://jaspar.elixir.no/api/v1/matrix/MA2503.1.meme" -O /data/genome/motifs/jaspar_2024/meme/vertebrata/MA2503.1.meme.  <--DIDNT DO IT BECAUSE NO PERMISSION

wget "https://jaspar.elixir.no/api/v1/matrix/MA2503.1.meme" -O ./motifs/MA2503.1.meme

###############################################################################
# BANP FIMO SCANNING
###############################################################################
#!/usr/bin/env bash
set -euo pipefail

#  module load meme/5.5.7

# ===============================
# VARIABLES (ONLY THINGS TO SET)
# ===============================
motif=BANP
genome=hg38
pval=0.001
# Location of the MEME file
group=motifs
id=MA2503.1         
GENOME_FA=/data/genome/genomes/${genome}.fa

# Output directories
mkdir -p ./motifs/

# ===============================
# STEP 1 — RUN FIMO
# ===============================
fimo --text --thresh ${pval} ./${group}/${id}.meme ${GENOME_FA} 2> /dev/null | awk '$1=="motif_id"{next}{ print $3, $4-1, $5, toupper($9), $8, $6}' \
| awk 'BEGIN{f=0}{if(s!="" && ($2!=s || $3!=e) && f==0){print l} else if($2==s && $3==e){ if($5<p){print $0}else{print l}f=1} else{ f=0 } s=$2; e=$3; p=$5; l=$0} \
  END{if(f==0){print l}}' \
| gzip > ./motifs/${motif}.bed.gz

echo "[OK] FIMO finished for ${motif}"

# ===============================
# STEP 2 — ADD MATCHED SEQUENCE
# ===============================
TMP=$(mktemp)

paste \
  <(zcat ./motifs/${motif}.bed.gz) \
  <(
    zcat ./motifs/${motif}.bed.gz \
    | fastaFromBed -s -fi ${GENOME_FA} -bed stdin \
    | grep -v ">"
  ) \
| awk '{print $1, $2, $3, toupper($7), $5, $6}' \
| gzip > "${TMP}"

mv "${TMP}" ./motifs/${motif}.bed.gz

echo "[DONE] BANP motif BED with sequence:"
echo "       ./motifs/${motif}.bed.gz"

###############################################################################
# 1) MOTIF FILTERING (6-MER STRINGENCY)
###############################################################################
# Filter motifs with p-value threshold and sort them.
mkdir -p ./motifs

zcat /data/genome/motifs/jaspar_2024/fimo/vertebrata/hg38/NRF1.MA0506.3.bed.gz | awk '($5<=1/4^6)' > ./motifs/NRF1_filtered_6mer.MA0506.3.bed
zcat ./motifs/BANP.bed.gz | awk '($5<=1/4^6)' > ./motifs/BANP_filtered_6mer.MA2503.1.bed

###############################################################################
# MOTIF COUNT BEFORE / AFTER FILTERING
###############################################################################
# barplot of the number of motifs before and after filtering 
Rscript -e '
library(ggplot2)
motifs <- c("NRF1", "BANP")
before_counts <- c( nrow(read.table(gzfile("/data/genome/motifs/jaspar_2024/fimo/vertebrata/hg38/NRF1.MA0506.3.bed.gz"))),
                    nrow(read.table(gzfile("./motifs/BANP.bed.gz"))) )
after_counts  <- c( nrow(read.table("./motifs/NRF1_filtered_6mer.MA0506.3.bed")),
                    nrow(read.table("./motifs/BANP_filtered_6mer.MA2503.1.bed")) )

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

###############################################################################
# MOTIF COUNT AFTER FILTERING
###############################################################################
# barplot to see only the number of motifs after filtering
Rscript -e '
library(ggplot2)
motifs <- c("NRF1", "BANP")
after_counts  <- c( nrow(read.table("./motifs/NRF1_filtered_6mer.MA0506.3.bed")),
                    nrow(read.table("./motifs/BANP_filtered_6mer.MA2503.1.bed")) )
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

###############################################################################
# MOTIF LIST FOR EACH TF 
###############################################################################
mkdir -p ./motifs/motifs_list

for file in ./motifs/*_filtered_6mer.*.bed; do
  awk '{print $4}' "$file" | sort -u > "./motifs/motifs_list/$(basename "$file" .bed)_motifs.txt"
done

###############################################################################
# MOTIF MISMATCH POSITION LIST FOR EACH TF 
###############################################################################
Rscript -e '
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(tidyr)
})

# ===============================
# INPUTS
# ===============================
files <- list(
  NRF1 = "./motifs/motifs_list/NRF1_filtered_6mer.MA0506.3_motifs.txt",
  BANP = "./motifs/motifs_list/BANP_filtered_6mer.MA2503.1_motifs.txt"
)

NRF1motif <- "CTGCGCATGCGC"   # length 12
BANPmotif <- "TCTCGCGAGA"     # length 10

out_pdf <- "./results/summarys/motif_variant_mismatch_visual_summary.pdf"

# ===============================
# HELPERS
# ===============================
split_seq <- function(x) do.call(rbind, strsplit(x, ""))

# ===============================
# PDF OUTPUT
# ===============================
pdf(out_pdf, width = 10, height = 7)

for (mot in names(files)) {

  # -----------------------------
  # LOAD UNIQUE MOTIF VARIANTS
  # -----------------------------
  seqs <- unique(fread(files[[mot]], header = FALSE)$V1)

  ref <- if (mot == "NRF1") NRF1motif else BANPmotif
  L   <- nchar(ref)

  stopifnot(all(nchar(seqs) == L))

  hit_mat <- split_seq(seqs)
  ref_mat <- split_seq(rep(ref, length(seqs)))

  mismatch_mat <- (hit_mat != ref_mat)
  n_variants <- nrow(mismatch_mat)

  # -----------------------------
  # CORE POSITIONS
  # -----------------------------
  if (mot == "NRF1") {
    core_pos <- c(4, 5, 10, 11)   # CG + CG
  } else {
    core_pos <- c(4, 5, 6, 7)     # CGCG
  }

  all_pos    <- seq_len(L)
  before_pos <- all_pos[all_pos < min(core_pos)]
  after_pos  <- all_pos[all_pos > max(core_pos)]

  before <- rowSums(mismatch_mat[, before_pos, drop = FALSE])
  core   <- rowSums(mismatch_mat[, core_pos, drop = FALSE])
  after  <- rowSums(mismatch_mat[, after_pos, drop = FALSE])

  dt <- data.table(
    variant = seqs,
    before  = before,
    core    = core,
    after   = after
  )
  dt[, total := before + core + after]

  # -----------------------------
  # REMOVE PERFECT MATCHES
  # -----------------------------
  dt <- dt[total > 0]

  # =====================================================
  # PLOT 1 — POSITIONAL MISMATCH TOLERANCE
  # =====================================================
  pos_df <- data.frame(
    position   = seq_len(L),
    mismatches = colSums(mismatch_mat),
    total      = n_variants,
    base       = strsplit(ref, "")[[1]]
  )

  print(
    ggplot(pos_df, aes(x = position)) +
      geom_col(aes(y = total), fill = "grey85") +
      geom_col(aes(y = mismatches), fill = "firebrick") +
      geom_text(aes(y = -0.5, label = base), size = 5) +
      theme_minimal() +
      labs(
        title = paste(mot, "- positional mismatch tolerance"),
        subtitle = "Red = variants with mismatch, Grey = total variants",
        x = "Motif position",
        y = "Number of variants"
      ) +
      coord_cartesian(ylim = c(-1, n_variants)) +
      theme(
        panel.grid.major.x = element_blank(),
        axis.text.x = element_blank()
      )
  )

  # =====================================================
  # CONFIGURATION PER VARIANT (ONLY CHANGE IS HERE)
  # =====================================================
  dt[, config := fifelse(
    core == total, "all_core",
    fifelse(core == 0, "no_core",
            fifelse(core == 1, "mixed_1CG", "mixed_2CG"))
  )]

  cfg_pct <- dt[, .N, by = .(total, config)]
  cfg_pct[, pct := round(100 * N / sum(N)), by = total]
  cfg_pct <- cfg_pct[total > 0]

  labels_df <- cfg_pct[, .(
    label = paste0(
      "Core: ",
      ifelse(length(pct[config == "all_core"]) == 0, 0, pct[config == "all_core"]),
      "%\n",
      "Mixed 2CG: ",
      ifelse(length(pct[config == "mixed_2CG"]) == 0, 0, pct[config == "mixed_2CG"]),
      "%\n",
      "Mixed 1CG: ",
      ifelse(length(pct[config == "mixed_1CG"]) == 0, 0, pct[config == "mixed_1CG"]),
      "%\n",
      "NoCore: ",
      ifelse(length(pct[config == "no_core"]) == 0, 0, pct[config == "no_core"]),
      "%"
    )
  ), by = total]

  # =====================================================
  # PLOT 2 — TOTAL MISMATCHES PER VARIANT
  # =====================================================
  print(
    ggplot(dt, aes(x = total)) +
      geom_bar(fill = "grey40") +
      geom_text(
        data = labels_df,
        aes(x = total, y = Inf, label = label),
        vjust = 1.15,
        size = 3.1,
        lineheight = 0.95
      ) +
      theme_minimal() +
      labs(
        title = paste(mot, "- mismatch configuration per variant"),
        subtitle = "Core = CG only | Mixed 2CG = ≥2 CG + flanks | Mixed 1CG = 1 CG + flanks | NoCore = flanks only",
        x = "Number of mismatches",
        y = "Number of variants"
      )
  )
}

dev.off()

cat("[DONE] PDF written:", out_pdf, "\n")
'

###############################################################################
# 2) MOTIF DISTANCE TO TSS
###############################################################################
# Distance of motifs to TSS : generate a file per sample with the distances of each motif to the nearest TSS
mkdir -p ./motifs/dist_to_tss/

for f in ./motifs/*_filtered_6mer.*.bed; do
    base=$(basename "$f")                     
    sample=${base%_filtered_6mer.*.bed} 
    echo "Processing $sample ..."

    closestBed -t first -d -a "$f" -b /data/genome/annotations/hg38_tss.bed | awk '{print $NF}' > "./motifs/dist_to_tss/${sample}_dist_tss.txt"
done

###############################################################################
# HISTOGRAMS OF DISTANCE TO TSS (ONE PAGE PER TF)
###############################################################################
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

  # plot
  hist(log10(distances + 1),
       xlim = c(0, 8),
       xlab = "Distance to TSS (log10)",
       main = paste0(sample, "  |  Proximal = ", pct_prox, "%   |   Distal = ", pct_dist, "%"),
       breaks = 50,
       col = "steelblue",
       border = "black")
  # red vertical line at log10(2000) --> promotor 
  abline(v = log10(2000), col = "red", lwd = 2)}
dev.off()
' 

###############################################################################
# 3) MOTIF ∩ HM450 METHYLATION PROBES
###############################################################################
# Overlap 6mer motifs with annotated methylation data to find intersecting regions and 
# see if the motifs fall within methylated probes regions.

mkdir -p ./motifs/overlaps/intersected_motifs_HM450

# intersect NRF1 6mer motifs with annotated methylation data
bedtools intersect -u -a ./motifs/NRF1_filtered_6mer.MA0506.3.bed -b ./methylation/annotated_methylation_data_probes_filtered.bed > ./motifs/overlaps/intersected_motifs_HM450/NRF1_intersected_methylation.bed
bedtools intersect -u -a ./motifs/BANP_filtered_6mer.MA2503.1.bed -b ./methylation/annotated_methylation_data_probes_filtered.bed > ./motifs/overlaps/intersected_motifs_HM450/BANP_intersected_methylation.bed

# Summarize the intersected results : 1   NRF1 29401 BANP 2919
mkdir -p ./results/summarys

# Count lines in bash and save to variables
n_NRF1=$(cat ./motifs/overlaps/intersected_motifs_HM450/NRF1_intersected_methylation.bed | wc -l)
n_BANP=$(cat ./motifs/overlaps/intersected_motifs_HM450/BANP_intersected_methylation.bed | wc -l)
# Create a summary table using R
Rscript -e '
table <- data.frame(
  Motif = c("NRF1","BANP"),
  Overlapping_regions = c('"$n_NRF1"','"$n_BANP"'))
print(table)
'
###############################################################################
# VEN DIAGRAM OF MOTIF ∩ HM450 METHYLATION PROBES
###############################################################################
# Create a Venn diagram to visualize the overlaps between motifs and methylation for NRF1 and BANP : output pdf file with 2 pages : page 1 NRF1 , page 2 BANP.
Rscript -e '
library(VennDiagram)
library(grid)

motifs_NRF1 <- as.numeric(
  system("wc -l < ./motifs/NRF1_filtered_6mer.MA0506.3.bed", intern = TRUE)
)

overlap_NRF1 <- as.numeric(
  system("wc -l < ./motifs/overlaps/intersected_motifs_HM450/NRF1_intersected_methylation.bed", intern = TRUE)
)

motifs_BANP <- as.numeric(
  system("wc -l < ./motifs/BANP_filtered_6mer.MA2503.1.bed", intern = TRUE)
)

overlap_BANP <- as.numeric(
  system("wc -l < ./motifs/overlaps/intersected_motifs_HM450/BANP_intersected_methylation.bed", intern = TRUE)
)

probes <- as.numeric(
  system("wc -l < ./methylation/annotated_methylation_data_probes_filtered.bed", intern = TRUE)
)

# -------------------- OUTPUT PDF --------------------
out_pdf <- "./results/summarys/NRF1_BANP_motifs_vs_methylation_venn.pdf"
pdf(out_pdf, width = 13, height = 13)

# ==================== PAGE 1: NRF1 ====================
grid.newpage()

venn_NRF1 <- draw.pairwise.venn(
  area1      = motifs_NRF1,
  area2      = probes,
  cross.area = overlap_NRF1,
  category   = c("NRF1 6mer motifs", "Methylation probes"),
  fill       = c("salmon", "skyblue"),
  alpha      = 0.7,
  cex        = 2,
  cat.cex    = 2.2,
  fontfamily = "Helvetica",
  cat.pos    = c(-15, 20),
  cat.dist   = c(0.03, 0.02),
  print.mode = c("raw","percent"),
  sigdig     = 2,
  ind        = FALSE
)

grid.text(
  "Overlap between NRF1 6mer motifs and methylation probes",
  x = 0.5, y = 0.96,
  gp = gpar(fontsize = 30, fontfamily = "Helvetica")
)

pushViewport(viewport(y = 0.45))
grid.draw(venn_NRF1)
popViewport()

# ==================== PAGE 2: BANP ====================
grid.newpage()

venn_BANP <- draw.pairwise.venn(
  area1      = motifs_BANP,
  area2      = probes,
  cross.area = overlap_BANP,
  category   = c("", ""),
  fill       = c("salmon", "skyblue"),
  alpha      = 0.7,
  cex        = 1.9,
  fontfamily = "Helvetica",
  print.mode = c("raw","percent"),
  sigdig     = 2,
  ind        = FALSE
)

grid.text(
  "Overlap between BANP 6mer motifs and methylation probes",
  x = 0.5, y = 0.95,
  gp = gpar(fontsize = 30, fontfamily = "Helvetica")
)

grid.text(
  "Methylation probes",
  x = 0.37, y = 0.73,
  gp = gpar(fontsize = 28, fontfamily = "Helvetica")
)

grid.text(
  "BANP 6mer motifs",
  x = 0.75, y = 0.65,
  gp = gpar(fontsize = 28, fontfamily = "Helvetica")
)

pushViewport(viewport(y = 0.45))
grid.draw(venn_BANP)
popViewport()

# -------------------- CLOSE PDF --------------------
dev.off()
cat("Wrote:", out_pdf, "\n")
'

###############################################################################
# 4) MOTIF ∩ ATAC PEAKS
###############################################################################
# Overlap the 6mer motifs with the peaks to see which motifs fall within the peaks regions.
mkdir -p ./motifs/overlaps/motif_peak_overlaps

# Loop through each 6mer motif file and intersect with all peak files from peaks
for file in ./motifs/*_filtered_6mer*.bed; do
   motif_name=$(basename "$file" .bed)
   bedtools intersect -u -a "$file" -b ./peaks/filtered_peaks/merged_peaks.bed > ./motifs/overlaps/motif_peak_overlaps/"${motif_name}_peak_overlaps.bed"
done

# Summarize the intersected results : NRF1 376062 BANP 33226

# Count lines in bash and save to variables
n_NRF1=$(cat ./motifs/overlaps/motif_peak_overlaps/NRF1_filtered_6mer.MA0506.3_peak_overlaps.bed | wc -l)
n_BANP=$(cat ./motifs/overlaps/motif_peak_overlaps/BANP_filtered_6mer.MA2503.1_peak_overlaps.bed | wc -l)
# Create a summary table using R
Rscript -e '
table <- data.frame(
  Motif = c("NRF1","BANP"),
  Overlapping_regions = c('"$n_NRF1"','"$n_BANP"'))
print(table)
'

###############################################################################
# 5) MOTIF ∩ PEAK ∩ HM450
###############################################################################
# 6) Overlap the intersected methylation motifs with the peak overlapping motifs to find common regions
mkdir -p ./motifs/overlaps/intersected_overlaps

bedtools intersect -u -a ./motifs/overlaps/intersected_motifs_HM450/NRF1_intersected_methylation.bed -b ./motifs/overlaps/motif_peak_overlaps/NRF1_filtered_6mer.MA0506.3_peak_overlaps.bed > ./motifs/overlaps/intersected_overlaps/NRF1_methylation_peak_overlap.bed  
bedtools intersect -u -a ./motifs/overlaps/intersected_motifs_HM450/BANP_intersected_methylation.bed -b ./motifs/overlaps/motif_peak_overlaps/BANP_filtered_6mer.MA2503.1_peak_overlaps.bed > ./motifs/overlaps/intersected_overlaps/BANP_methylation_peak_overlap.bed

# Summarize the intersected results : NRF1 23416 BANP 2217

# Count lines in bash and save to variables
n_NRF1=$(cat ./motifs/overlaps/intersected_overlaps/NRF1_methylation_peak_overlap.bed | wc -l)
n_BANP=$(cat ./motifs/overlaps/intersected_overlaps/BANP_methylation_peak_overlap.bed | wc -l)
# Create a summary table using R
Rscript -e '
table <- data.frame(
  Motif = c("NRF1","BANP"),
  Overlapping_regions = c('"$n_NRF1"','"$n_BANP"'))
print(table)
'

###############################################################################
# 6) BARPLOTS – MOTIFS ACROSS DATASETS
###############################################################################

###############################################################################
# ONE BARPLOT PER TF : MOTIFS, ATAC, HM450, SNVs 
###############################################################################
# create a barplot : One page per TF: motifs, ATAC, HM450, SNVs
Rscript -e '
library(ggplot2)

out_pdf <- "./results/summarys/NRF1_BANP_3D_vs_peaks_vs_HM450_barplots.pdf"
dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)

pdf(out_pdf, width = 10, height = 8)

# ===================== PAGE 1: NRF1 =====================
## 1) Read motif sets and build unique IDs
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

## 3) Data frame
df <- data.frame(
  category = c("All NRF1 motifs",
               "Motifs in ATAC peaks",
               "Motifs with HM450 probes",
               "Motifs in both ATAC + HM450"),
  count = c(total, in_atac, in_hm450, in_both)
)

df$category <- factor(df$category, levels = df$category)
df$percent  <- df$count / total * 100
df$label    <- sprintf("%d (%.2f%%)", df$count, df$percent)

## 4) Plot
p1 <- ggplot(df, aes(x = reorder(category, -count), y = count, fill = category)) +
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
    subtitle = sprintf("All NRF1 6-mer motifs (n = %d)", total),
    y = "Number of motifs",
    x = NULL
  )

print(p1)

# ===================== PAGE 2: BANP =====================
## 1) Read motif sets and build unique IDs
motif_all   <- read.table("./motifs/BANP_filtered_6mer.MA2503.1.bed")
motif_ids   <- unique(with(motif_all, paste(V1, V2, V3, sep=":")))

motif_peak  <- read.table("./motifs/overlaps/motif_peak_overlaps/BANP_filtered_6mer.MA2503.1_peak_overlaps.bed")
peak_ids    <- unique(with(motif_peak, paste(V1, V2, V3, sep=":")))

motif_hm450 <- read.table("./motifs/overlaps/intersected_motifs_HM450/BANP_intersected_methylation.bed")
hm450_ids   <- unique(with(motif_hm450, paste(V1, V2, V3, sep=":")))

## 2) Counts
total   <- length(motif_ids)
in_atac <- length(peak_ids)
in_hm450<- length(hm450_ids)
in_both <- length(intersect(peak_ids, hm450_ids))

## 3) Data frame
df <- data.frame(
  category = c("All BANP motifs",
               "Motifs in ATAC peaks",
               "Motifs with HM450 probes",
               "Motifs in both ATAC + HM450"),
  count = c(total, in_atac, in_hm450, in_both)
)

df$category <- factor(df$category, levels = df$category)
df$percent  <- df$count / total * 100
df$label    <- sprintf("%d (%.2f%%)", df$count, df$percent)

## 4) Plot
p2 <- ggplot(df, aes(x = reorder(category, -count), y = count, fill = category)) +
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
    x = NULL
  )

print(p2)

# ===================== CLOSE =====================
dev.off()
cat("Wrote:", out_pdf, "\n")
'
# !!! DONE !!!
###############################################################################
# ONE BARPLOT FOR BOTH TF : MOTIFS, ATAC, HM450, SNVs 
###############################################################################
# one barplot for both TFs NRF1 and BANP together : this script creates a small barplot showing the percentage of motifs in each category for both TFs one tf per page
Rscript -e '
library(ggplot2)

# ===============================
# TF-specific file definitions
# ===============================
tf_files <- list(
  BANP = list(
    motif = "./motifs/BANP_filtered_6mer.MA2503.1.bed",
    peak  = "./motifs/overlaps/motif_peak_overlaps/BANP_filtered_6mer.MA2503.1_peak_overlaps.bed",
    hm450 = "./motifs/overlaps/intersected_motifs_HM450/BANP_intersected_methylation.bed",
    snv   = "./snv/overlaps/intersected_motifs_snv_filtered/BANP_filtered_SNVs_in_motifskept.bed"
  ),
  NRF1 = list(
    motif = "./motifs/NRF1_filtered_6mer.MA0506.3.bed",
    peak  = "./motifs/overlaps/motif_peak_overlaps/NRF1_filtered_6mer.MA0506.3_peak_overlaps.bed",
    hm450 = "./motifs/overlaps/intersected_motifs_HM450/NRF1_intersected_methylation.bed",
    snv   = "./snv/overlaps/intersected_motifs_snv_filtered/NRF1_filtered_SNVs_in_motifskept.bed"
  )
)

out_pdf <- "./results/summarys/BANP_NRF1_motif_peaks_probes_snv_barplots.pdf"
pdf(out_pdf, width = 10, height = 8)

for (tf in names(tf_files)) {

  ## 1) Motifs
  motif_all <- read.table(tf_files[[tf]]$motif)
  motif_ids <- unique(with(motif_all, paste(V1, V2, V3, sep=":")))

  ## 2) ATAC peak overlaps
  motif_peak <- read.table(tf_files[[tf]]$peak)
  peak_ids <- unique(with(motif_peak, paste(V1, V2, V3, sep=":")))

  ## 3) HM450 overlaps
  motif_hm450 <- read.table(tf_files[[tf]]$hm450)
  hm450_ids <- unique(with(motif_hm450, paste(V1, V2, V3, sep=":")))

  ## 4) SNV overlaps
  motif_snv <- read.table(tf_files[[tf]]$snv)
  snv_ids <- unique(with(motif_snv, paste(V1, V2, V3, sep=":")))

  ## 5) Counts
  total         <- length(motif_ids)
  in_atac       <- length(peak_ids)
  in_hm450      <- length(hm450_ids)
  in_snv        <- length(snv_ids)
  in_peak_hm450 <- length(intersect(peak_ids, hm450_ids))
  in_peak_snv   <- length(intersect(peak_ids, snv_ids))

  ## 6) Data frame (UNCHANGED logic)
  df <- data.frame(
    category = c("All motifs",
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
              in_peak_hm450)
  )

  df$category <- factor(df$category, levels = df$category)
  df$percent  <- df$count / total * 100
  df$label    <- sprintf("%d (%.2f%%)", df$count, df$percent)

  ## 7) Plot (exact same style)
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
      title = paste0(tf, " 6-mer motif distribution across datasets"),
      subtitle = sprintf("All %s 6-mer motifs (n = %d)", tf, total),
      y = "Number of motifs",
      x = NULL
    )

  print(p)
}

dev.off()
cat("DONE →", out_pdf, "\n")
'
#!!! DONE  !!!
###############################################################################
# PROPORTIONAL EULER DIAGRAM PER TF : MOTIFS, ATAC, HM450, SNVs 
###############################################################################
# Proportional Euler diagram of motifs with peaks and methylation probes 
# Venn diagram 3D with filtered SNVs in motifs and in peak 
# Proportional Euler diagram of motifs with peaks and methylation probes

Rscript -e '
library(eulerr)

out_pdf <- "./results/summarys/Euler_VENN_NRF1_BANP_Motifs_ATAC_HM450_SNV.pdf"
dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)

pdf(out_pdf, width = 10, height = 10)

# margins: bottom, left, top, right
par(mar = c(4, 4, 6, 2))

################################
## PAGE 1 — NRF1 : ATAC + SNVs ##
################################
plot.new()   

motif_all <- unique(with(
  read.table("./motifs/NRF1_filtered_6mer.MA0506.3.bed"),
  paste(V1, V2, V3, sep=":")
))
motif_peak <- unique(with(
  read.table("./motifs/overlaps/motif_peak_overlaps/NRF1_filtered_6mer.MA0506.3_peak_overlaps.bed"),
  paste(V1, V2, V3, sep=":")
))
motif_variants <- unique(with(
  read.table("./snv/overlaps/intersected_motifs_snv_filtered/NRF1_filtered_SNVs_in_motifskept.bed"),
  paste(V1, V2, V3, sep=":")
))

fit <- euler(list(
  "All NRF1 motifs"       = motif_all,
  "Motifs in ATAC"        = motif_peak,
  "Motifs with variants" = motif_variants
))

plot(fit,
     fills = c("steelblue","lightgreen","red"),
     edges = "black",
     quantities = list(type="counts", cex=0.9),
     labels = list(cex=1.1),
     main = NULL)

mtext("NRF1 motifs: genome-wide, accessible,\nand overlapping SNVs",
      side = 3, line = 2, cex = 1.2)

##################################
## PAGE 2 — NRF1 : ATAC + HM450 ##
##################################
plot.new()

motif_hm450 <- unique(with(
  read.table("./motifs/overlaps/intersected_motifs_HM450/NRF1_intersected_methylation.bed"),
  paste(V1, V2, V3, sep=":")
))

fit <- euler(list(
  "All NRF1 motifs"          = motif_all,
  "Motifs in ATAC"           = motif_peak,
  "Motifs with HM450 probes" = motif_hm450
))

plot(fit,
     fills = c("steelblue","lightgreen","orange"),
     edges = "black",
     quantities = list(type="counts", cex=0.9),
     labels = list(cex=1.1),
     main = NULL)

mtext("NRF1 motifs: genome-wide, accessible,\nand covered by HM450 probes",
      side = 3, line = 2, cex = 1.2)

################################
## PAGE 3 — BANP : ATAC + SNVs ##
################################
plot.new()

motif_all <- unique(with(
  read.table("./motifs/BANP_filtered_6mer.MA2503.1.bed"),
  paste(V1, V2, V3, sep=":")
))
motif_peak <- unique(with(
  read.table("./motifs/overlaps/motif_peak_overlaps/BANP_filtered_6mer.MA2503.1_peak_overlaps.bed"),
  paste(V1, V2, V3, sep=":")
))
motif_variants <- unique(with(
  read.table("./snv/overlaps/intersected_motifs_snv_filtered/BANP_filtered_SNVs_in_motifskept.bed"),
  paste(V1, V2, V3, sep=":")
))

fit <- euler(list(
  "All BANP motifs"       = motif_all,
  "Motifs in ATAC"        = motif_peak,
  "Motifs with variants" = motif_variants
))

plot(fit,
     fills = c("steelblue","lightgreen","red"),
     edges = "black",
     quantities = list(type="counts", cex=0.9),
     labels = list(cex=1.1),
     main = NULL)

mtext("BANP motifs: genome-wide, accessible,\nand overlapping SNVs",
      side = 3, line = 2, cex = 1.2)

##################################
## PAGE 4 — BANP : ATAC + HM450 ##
##################################
plot.new()

motif_hm450 <- unique(with(
  read.table("./motifs/overlaps/intersected_motifs_HM450/BANP_intersected_methylation.bed"),
  paste(V1, V2, V3, sep=":")
))

fit <- euler(list(
  "All BANP motifs"          = motif_all,
  "Motifs in ATAC"           = motif_peak,
  "Motifs with HM450 probes" = motif_hm450
))

plot(fit,
     fills = c("steelblue","lightgreen","orange"),
     edges = "black",
     quantities = list(type="counts", cex=0.9),
     labels = list(cex=1.1),
     main = NULL)

mtext("BANP motifs: genome-wide, accessible,\nand covered by HM450 probes",
      side = 3, line = 2, cex = 1.2)

dev.off()
cat("DONE →", out_pdf, "\n")
'

###############################################################################

###############################################################################
# ATAC PEAKS SCRIPT
###############################################################################

###############################################################################
# 1) LINK PEAK FILES FROM /data/hichamif/pred_tf_cancer/peak/
###############################################################################
# link the peak files from the previous analysis of /data/hichamif/pred_tf_cancer/peak/ to this folder :
# mkdir -p ./peaks
# for file in /data/hichamif/pred_tf_cancer/peak/ATAC_TCGA*_peaks_macs.bed; do
#    ln -s "$file" ./peaks/
# done

###############################################################################
# 2) FILTERING OUT PEAKS FILES --> 217 FILES LEFT
###############################################################################
# Filter out the peak files corresponding to failing samples from the previous analysis of /data/hichamif/pred_tf_cancer/QC_results/failing_samples_by_step/all_failing_samples.txt
# Filters : NFR 15% Mapping rate > 80% Peakcount >20,000 FRIP > 0.2
mkdir -p ./peaks/filtered_peaks/

FAIL_TSV="./all_failing_samples.txt"

for file in ./peaks/ATAC_TCGA-*_peaks_macs.bed; do
    base=$(basename "$file")
    sample=${base%_peaks_macs.bed}

    echo "Processing $sample ..."

    if awk -F'\t' -v s="$sample" 'NR>1 && $1==s {found=1} END{exit !found}' "$FAIL_TSV"; then
        echo "Skipping $sample (in failing list)."
        continue
    fi

    cp "$file" ./peaks/filtered_peaks/
done

###############################################################################
# 3) COMBINING ALL PEAK FILES INTO ONE 
###############################################################################
# combine all peak files into one and sort them
cat ./peaks/filtered_peaks/ATAC_TCGA*_peaks_macs.bed | sort -k1,1 -k2,2n | bedtools merge -i stdin > ./peaks/filtered_peaks/merged_peaks.bed

###############################################################################
# 4) STATISTICS ON FILTERED PEAKS DATA
###############################################################################

###############################################################################
# NBR OF CANCER TYPES WITH PEAKS DATA
###############################################################################
# nbr of cancer types with peaks --> 22 types 
ls ./peaks/filtered_peaks/ATAC_TCGA*bed | cut -d'-' -f2 | cut -d'_' -f1 | sort | uniq -c | awk '{print $2 "\t" $1}' > ./peaks/cancer_counts.tsv

###############################################################################
# CANCER COLOR ORDER FILE 
###############################################################################
# create a cancer color order file based on the order of cancers in the peaks and snv unique per cancer files
cat ./snv/snv_counts_per_cancer_type/snv_unique_per_cancer.tsv \
    ./peaks/peaks_counts_per_cancer_type/peaks_unique_per_cancer.tsv \
    | cut -f1 \
    | sort -u \
    > ./results/summarys/cancer_color_order.txt

###############################################################################
# HISTOGRAM OF NBR OF PEAKS PER CANCER TYPE
###############################################################################
# Plot histogram of nbr of Peaks per cancer type
# count peaks number per cancer type we have to merge all peak files per cancer type first and take the count of merged peaks

# Get list of cancer types from peak filenames
# Filenames look like: ATAC_TCGA-BRCA_TCGA-AR-A0TP_1_peaks_macs.bed
mkdir -p ./peaks/peaks_counts_per_cancer_type/

# 0) empty the output file (or create it)
> ./peaks/peaks_counts_per_cancer_type/peaks_unique_per_cancer.tsv

# 1) list cancer codes from filenames
#    ATAC_TCGA-BRCA_TCGA-AR-A0TP_1_peaks_macs.bed  →  TCGA-BRCA
cancers=$(ls ./peaks/filtered_peaks/ATAC_TCGA-*_*_peaks_macs.bed \
          | sed 's#.*/ATAC_##; s/_.*//' \
          | sort -u)

# 2) for each cancer, merge all peaks and count merged intervals
for cancer in $cancers; do
    echo "→ Processing $cancer ..."

    count=$(cat ./peaks/filtered_peaks/ATAC_${cancer}_*_peaks_macs.bed \
            | cut -f1-3 \
            | sort -k1,1 -k2,2n \
            | bedtools merge -i stdin \
            | wc -l)

    printf "%s\t%d\n" "$cancer" "$count" >> ./peaks/peaks_counts_per_cancer_type/peaks_unique_per_cancer.tsv
done

cut -d "-" -f2 ./peaks/peaks_counts_per_cancer_type/peaks_unique_per_cancer.tsv > ./peaks/peaks_counts_per_cancer_type/tmp.tsv
mv ./peaks/peaks_counts_per_cancer_type/tmp.tsv ./peaks/peaks_counts_per_cancer_type/peaks_unique_per_cancer.tsv

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

###############################################################################
# NBR OF PEAKS PER SAMPLE PER CANCER TYPE
###############################################################################
# compute number of peaks per sample per cancer type :
mkdir -p ./peaks/peak_counts_per_sample/
echo -e "sample\tcancer_type\tn_peaks" > ./peaks/peak_counts_per_sample/peak_counts_per_sample.tsv
for peakfile in ./peaks/filtered_peaks/ATAC_TCGA-*-*.bed; do
  sample=$(basename "$peakfile" .bed | sed 's/^Peaks_//')
  cancer=$(echo "$sample" | cut -d'-' -f2 | cut -d'_' -f1)
  n_peaks=$(cat "$peakfile" | grep -v "^#" | wc -l)

  echo -e "${sample}\t${cancer}\t${n_peaks}" >> ./peaks/peak_counts_per_sample/peak_counts_per_sample.tsv
done

# barplot of number of peaks per sample per cancer type
Rscript -e '
library(ggplot2)

in_tsv  <- "./peaks/peak_counts_per_sample/peak_counts_per_sample.tsv"
out_pdf <- "./results/summarys/ATAC_peak_counts_per_sample_by_cancer.pdf"

df <- read.table(in_tsv, header=FALSE, sep="\t", stringsAsFactors=FALSE)
colnames(df) <- c("sample","cancer","n_peaks")
df$n_peaks <- as.numeric(df$n_peaks)

dir.create(dirname(out_pdf), recursive=TRUE, showWarnings=FALSE)

pdf(out_pdf, width=12, height=7)

for (c in sort(unique(df$cancer))) {
  sub <- df[df$cancer == c, ]
  sub <- sub[order(sub$n_peaks, decreasing=TRUE), ]
  sub$sample <- factor(sub$sample, levels=sub$sample)

  p <- ggplot(sub, aes(x=sample, y=n_peaks)) +
    geom_bar(stat="identity") +
    theme_minimal(base_size=12) +
    theme(
      axis.text.x = element_text(angle=60, hjust=1, size=6),
      plot.title  = element_text(hjust=0.5)
    ) +
    labs(
      title = paste0(c, " — number of ATAC peaks per sample"),
      x = "Sample",
      y = "Number of peaks"
    )

  print(p)
}

dev.off()

cat("Wrote:", out_pdf, "\n")
'

###############################################################################
# 5) DISTANCE OF FILTERED PEAKS TO TSS 
###############################################################################
# Distance of filtered peaks to TSS : generate a file per sample with the distances of each peak to the nearest TSS
mkdir -p ./peaks/dist_to_tss/

for f in ./peaks/filtered_peaks/ATAC_TCGA*peaks_macs.bed; do
    sample=$(basename "$f" _peaks_macs.bed)

    echo "Processing $sample ..."

    closestBed -t first -d -a "$f" -b /data/genome/annotations/hg38_tss.bed | awk '{print $NF}' > "./peaks/dist_to_tss/${sample}_dist_tss.txt"
done

###############################################################################
# HISTOGRAM OF DISTANCE OF FILTERED PEAKS TO TSS 
###############################################################################
# Plot histogram of distances to TSS for all samples combined

Rscript -e 'library(ggplot2);
peak_files <- list.files("peaks/filtered_peaks", pattern = "_peaks_macs.bed$", full.names = TRUE);
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

###############################################################################
#  BOXPLOT OF % PROXIMAL AND DISTAL PEAKS ACROSS ALL SAMPLES 
###############################################################################
# plot a boxplot of the percentage of proximal vs distal peaks across all samples
Rscript -e '
library(ggplot2)

# Load all samples 
peak_files <- list.files("peaks/filtered_peaks", pattern = "_peaks_macs.bed$", full.names = TRUE)
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

###############################################################################
# 6) ATAC PEAKS ∩ MOTIFS OVERLAPPING PEAKS KEPT
###############################################################################
# Look into the ATAC peaks overlapping with motifs
mkdir -p ./peaks/overlaps/intersected_motifs/NRF1

motif_sites="./motifs/NRF1_filtered_6mer.MA0506.3.bed"

for f in ./peaks/filtered_peaks/ATAC_TCGA*peaks_macs.bed; do
  sample=$(basename "$f" _peaks_macs.bed)
  out="./peaks/overlaps/intersected_motifs/NRF1/${sample}_NRF1_intersected_peaks.bed"

  bedtools intersect -u -a "$f" -b "$motif_sites" > "$out"
  echo "Wrote: $out"
done

mkdir -p ./peaks/overlaps/intersected_motifs/BANP

motif_sites="./motifs/BANP_filtered_6mer.MA2503.1.bed"

for f in ./peaks/filtered_peaks/ATAC_TCGA*peaks_macs.bed; do
  sample=$(basename "$f" _peaks_macs.bed)
  out="./peaks/overlaps/intersected_motifs/BANP/${sample}_BANP_intersected_peaks.bed"
  bedtools intersect -u -a "$f" -b "$motif_sites" > "$out"
  echo "Wrote: $out"
done


###############################################################################
# 7) MOTIF THAT OVERLAP WITH PEAKS DISTANCE TO TSS
###############################################################################
# Distance of motifs in peaks to TSS : generate a file per sample with the distances of each motif to the nearest TSS

# NRF1 
mkdir -p ./peaks/dist_to_tss_motifs/NRF1

for f in ./peaks/overlaps/intersected_motifs/NRF1/*_NRF1_intersected_peaks.bed; do
  base=$(basename "$f" _NRF1_intersected_peaks.bed)
  out="./peaks/dist_to_tss_motifs/NRF1/${base}_NRF1_dist_tss.txt"

  closestBed -t first -d \
    -a "$f" \
    -b /data/genome/annotations/hg38_tss.bed \
  | awk '{print $NF}' > "$out"

  echo "Wrote: $out"
done

# BANP
mkdir -p ./peaks/dist_to_tss_motifs/BANP

for f in ./peaks/overlaps/intersected_motifs/BANP/*_BANP_intersected_peaks.bed; do
  base=$(basename "$f" _BANP_intersected_peaks.bed)
  out="./peaks/dist_to_tss_motifs/BANP/${base}_BANP_dist_tss.txt"

  closestBed -t first -d \
    -a "$f" \
    -b /data/genome/annotations/hg38_tss.bed \
  | awk '{print $NF}' > "$out"

  echo "Wrote: $out"
done

# Plot histogram of distances to TSS for NRF1 and BANP motifs per cancer type
Rscript -e '
library(data.table)

out_pdf <- "./results/summarys/ATAC_motifs_TSS/dist_tss_ATAC_NRF1_BANP_byCancer.pdf"
dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)

pdf(out_pdf, width = 11.7, height = 8.3)
par(bg = "white")

# =========================
# RUN ONCE FOR EACH TF
# =========================
for (TF in c("NRF1", "BANP")) {

  dist_dir <- paste0("./peaks/dist_to_tss_motifs/", TF)
  out_tsv  <- paste0("./results/summarys/ATAC_motifs_TSS/", TF,"/summary_", TF, "_byCancer.tsv")

  files <- list.files(dist_dir, pattern = paste0("_", TF, "_dist_tss.txt$"),full.names = TRUE)
  if (!length(files)) {
    plot.new()
    title(main = paste0(TF, ": no distance files found"))
    next
  }

  # extract cancer code from filename: TCGA-XXX-
  get_cancer <- function(x) {
    m <- regexpr("TCGA-[A-Z0-9]+", x)
    if (m[1] == -1) return(NA_character_)
    sub("^TCGA-", "", regmatches(x, m))
  }


  cancers <- vapply(basename(files), get_cancer, character(1))
  ok <- !is.na(cancers) & cancers != ""
  files <- files[ok]
  cancers <- cancers[ok]

  res <- list()
  uniq_cancers <- sort(unique(cancers))

  for (c in uniq_cancers) {
    fc <- files[cancers == c]

    d_list <- lapply(fc, function(f)
      suppressWarnings(as.numeric(readLines(f))))
    d <- unlist(d_list, use.names = FALSE)
    d <- d[is.finite(d)]

    n_samples <- length(fc)
    n_vals <- length(d)

    if (n_vals == 0) {
      plot.new()
      title(main = paste0("TCGA-", c, " | ", TF, " motif CpGs: none found"))
      res[[c]] <- data.frame(cancer=c, n_samples=n_samples, n_values=0,pct_prox=NA, pct_dist=NA, median=NA)
      next
    }

    proximal <- sum(d < 2000)
    distal   <- sum(d >= 2000)
    pct_prox <- round(100 * proximal / n_vals, 1)
    pct_dist <- round(100 * distal   / n_vals, 1)

    hist(log10(d + 1),
         xlim = c(0, 8),
         xlab = "Distance to TSS (log10)",
         main = paste0("TCGA-", c, " | ", TF,
                       " | samples=", n_samples,
                       " | n=", n_vals,
                       " | Proximal=", pct_prox,
                       "% | Distal=", pct_dist, "%"),
         breaks = 50,
         col = "steelblue",
         border = "black")
    abline(v = log10(2000), col = "red", lwd = 2)

    res[[c]] <- data.frame(cancer=c, n_samples=n_samples, n_values=n_vals,pct_prox=pct_prox, pct_dist=pct_dist, median=median(d))
  }

  dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)
  write.table(do.call(rbind, res), out_tsv,sep = "\t", quote = FALSE, row.names = FALSE)

}

dev.off()

cat("Wrote:", out_pdf, "\n")
'

###############################################################################

###############################################################################
# METHYLATION SCRIPT
###############################################################################

###############################################################################
# SECTION 1 — SETUP & DATA LINKING
###############################################################################

# 1) Create methylation folder
mkdir -p ./methylation
# 2) Download the HM450 probes bed file :
#wget -O ./methylation/HM450.hg38.manifest.tsv.gz \
#https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/HM450/HM450.hg38.manifest.tsv.gz

###############################################################################
# SECTION 2 — LINK RAW TCGA HM450 FILES
###############################################################################

# Linking the methylation data from /data/papers/tcga to the methylation folder 
mkdir -p ./methylation/methylation_data/

for file in /data/papers/tcga/TCGA*/*/HM450*.txt; do
    ln -s "$file" ./methylation/methylation_data/
done
##!!!Done!!!
###############################################################################
# SECTION 3 — ANNOTATE HM450 FILES WITH GENOMIC COORDINATES FROM MANIFEST FILE 
###############################################################################
# Linking the methylation data with the HM450 manifest file for annotation : chr start end probe_id(gc) beta_value
# for each methylation file, annotate it and save it as a gzipped bed file in the same input directory
ls /data/papers/tcga/TCGA*/*/HM450*.txt \
| xargs -n 1 -P 5 bash -c '
file="$1"
sample_name=$(basename "$file" .txt)
out="$(dirname "$file")/${sample_name}_annotated_methylation.bed.gz"

echo "[PID $$] Processing $file"

awk "
BEGIN {
  FS = OFS = \"\t\"
  while (getline < \"./methylation/annotated_methylation_data_probes.bed\") {
    a[\$4] = \$0
  }
}
(\$1 in a) {
  if (\$2 == \"NA\" || \$2 == \"\") {
    print a[\$1], \"NA\"
  } else {
    print a[\$1], \$2 * 100
  }
}
" "$file" \
| sort -k1,1 -k2,2n \
| gzip > "$out"
' _


# for one file : 
# file=/data/papers/tcga/TCGA-STAD/TCGA-VQ-A91Z/HM450_TCGA-STAD-TCGA-VQ-A91Z-01A_1.txt
# sample_name=$(basename "$file" .txt)
# awk 'BEGIN{while(getline<"./methylation/annotated_methylation_data_probes.bed"){a[$4]=$0}}($1 in a){print a[$1],$2*100}' $file  | sort -k1,1 -k2,2n | gzip > "$(dirname "$file")/${sample_name}_annotated_methylation.bed.gz"

# build pairs files with sample healthy cancer pairs for each cancer type
set -euo pipefail
shopt -s nullglob

mkdir -p ./methylation/sample_pairs_files
out=./methylation/sample_pairs_files/methylation_pairs.tsv
echo -e "cancer\tpatient\ttumor_file\thealthy_file" > "$out"

for cancer_dir in /data/papers/tcga/TCGA-*; do
  cancer=$(basename "$cancer_dir" | sed 's/^TCGA-//')

  for patient_dir in "$cancer_dir"/*; do
    [ -d "$patient_dir" ] || continue
    patient=$(basename "$patient_dir")

    tumors=(  "$patient_dir"/*-01A_*_annotated_methylation.bed.gz )
    normals=( "$patient_dir"/*-11A_*_annotated_methylation.bed.gz )

    # Skip unless exactly one of each (01A vs 11A only, one pair per patient)
    [ "${#tumors[@]}" -eq 1 ] || continue
    [ "${#normals[@]}" -eq 1 ] || continue

    echo -e "$cancer\t$patient\t${tumors[0]}\t${normals[0]}" >> "$out"
  done
done

echo "Wrote: $out"

# transform the HM450 manifest annotated methylation data to bed format. chr start end probe_id --> 485569 probes
zcat ./methylation/HM450.hg38.manifest.tsv.gz | awk -vOFS="\t" ' NR>1 && $1 != "NA" && $2 != "NA" && $3 != "NA" {print $1, $2, $3, $9}' | sort -k1,1 -k2,2n  > ./methylation/annotated_methylation_data_probes.bed

###############################################################################
# SECTION 4 — FILTER METHYLATION FILES AND CREATE SAMPLE PAIRS FILES
###############################################################################

# Filter methylation files: filter CpG probes to remove chrX chrY chrM and keep only unique probes based on probe ID (4th column) and keep only probes starting with "cg" (to remove control probes)
# after filtering we have 470851 probes (from original 485569 probes).
# filter methylation files
# removes chrX chrY chrM, keeps only unique probes based on probe ID (4th column), keeps only probes starting with "cg" (to remove control probes)
mkdir -p ./methylation/filtered_methylation/
ls /data/papers/tcga/TCGA*/*/*_annotated_methylation.bed.gz | xargs -n 1 -P 5 bash -c '
file="$0"
out=./methylation/filtered_methylation/$(basename "$file" .bed.gz)_filtered.bed.gz

zcat "$file" | awk "BEGIN{FS=OFS=\"\t\"}
       \$1!=\"chrX\" && \$1!=\"chrY\" && \$1!=\"chrM\" &&
       !seen[\$4]++ &&
       \$4 ~ /^cg/ {print}" \
| gzip > "$out"

echo "[PID $$] Filtered $(basename "$file") -> $(basename "$out")"
'

# filter out probes that are in chrM and chrY and chrX and duplicated probes and keep only cg probes --> 470851 probes
cat ./methylation/annotated_methylation_data_probes.bed | awk 'BEGIN{FS=OFS="\t"} $1!="chrX" && $1!="chrY" && $1!="chrM" && !seen[$4]++ && $4 ~ /^cg/ {print}' > ./methylation/annotated_methylation_data_probes_filtered.bed

# create a list of all methylation files for cancer samples includes all tumor types and samples while excluding healthy samples and including only primary vial of tumor samples (A)
find ./methylation/filtered_methylation/* -type f -name "HM450*_annotated_methylation_filtered.bed.gz" | grep -E -- '-0[0-9]A_' > ./methylation/methylation_files_tumor0x.txt

# create a list of all filtered methylation files for cancer samples only 01A_1 tumor type first replicate and samples while excluding healthy samples 
find ./methylation/filtered_methylation/ -type f -name "HM450*_annotated_methylation_filtered.bed.gz" | grep -E -- '-01A_1' > ./methylation/filtered_methylation_files_tumor01A_1.txt

# build pairs files with sample healthy cancer pairs for each cancer type from filtered methylation files only 01A tumor and 11A healthy first replicates _1
# this file contains the list of all tumor-healthy pairs for methylation analysis for 22 cancer types that had at least one tumor-healthy pair.
set -euo pipefail
shopt -s nullglob

in_dir="./methylation/filtered_methylation"
mkdir -p ./methylation/sample_pairs_files
out="./methylation/sample_pairs_files/methylation_pairs_filtered.tsv"

echo -e "cancer\tpatient\ttumor_file\thealthy_file" > "$out"

# Loop ONLY over tumor files ending with -01A_1_...
for tumor in "$in_dir"/HM450_TCGA-*-TCGA-*-01A_1_annotated_methylation_filtered.bed.gz; do
  bn=$(basename "$tumor")

  # Extract cancer (e.g., TCGA-BRCA)
  cancer=$(echo "$bn" | sed -n 's/^HM450_\(TCGA-[A-Z0-9]\+\)-.*$/\1/p')

  # Extract patient barcode (e.g., TCGA-A2-A25C)
  patient=$(echo "$bn" | sed -n 's/^HM450_TCGA-[A-Z0-9]\+-\(TCGA-[A-Z0-9]\+-[A-Z0-9]\+\)-01A_1_.*$/\1/p')

  [[ -n "$cancer" && -n "$patient" ]] || continue

  # Require normal ALSO to be exactly -11A_1_...
  normal="$in_dir/HM450_${cancer}-${patient}-11A_1_annotated_methylation_filtered.bed.gz"

  # Skip if that exact normal _1 file doesn't exist
  [[ -f "$normal" ]] || continue

  echo -e "${cancer}\t${patient}\t${tumor}\t${normal}" >> "$out"
done

echo "Wrote: $out"

# build pairs files of replicates cancer type from filtered methylation files --> ./methylation/sample_pairs_files/methylation_replicates_pairs_filtered.tsv
# Build replicate tumor–tumor (01A_1 vs 01A_2) methylation pairs

set -euo pipefail
shopt -s nullglob

in_dir="./methylation/filtered_methylation"
out_dir="./methylation/sample_pairs_files"
out="$out_dir/methylation_replicates_pairs_filtered.tsv"

mkdir -p "$out_dir"

# Write header ONCE
echo -e "cancer\tpatient\treplicate_1\treplicate_2" > "$out"

# Loop over tumor replicate 1 files ONLY
for rep1 in "$in_dir"/HM450_TCGA-*-TCGA-*-01A_1_annotated_methylation_filtered.bed.gz; do
  bn=$(basename "$rep1")

  # Extract cancer (e.g. TCGA-BRCA)
  cancer=$(echo "$bn" | sed -E 's/^HM450_(TCGA-[A-Z0-9]+)-.*/\1/')

  # Extract patient (e.g. TCGA-A2-A25C)
  patient=$(echo "$bn" | sed -E 's/^HM450_TCGA-[A-Z0-9]+-(TCGA-[A-Z0-9]+-[A-Z0-9]+)-01A_1_.*/\1/')

  # Safety check
  [[ -n "$cancer" && -n "$patient" ]] || continue

  # Expected replicate 2 file
  rep2="$in_dir/HM450_${cancer}-${patient}-01A_2_annotated_methylation_filtered.bed.gz"

  # Skip if replicate 2 does not exist
  [[ -f "$rep2" ]] || continue

  echo -e "${cancer}\t${patient}\t${rep1}\t${rep2}" >> "$out"
done

echo "Wrote: $out"

###############################################################################
# SECTION 5 — QC: CONSISTENCY OF PROBE COUNTS IN THE /PAPERS/TCGA/ FOLDER
###############################################################################

# to see how many were generated : ls /data/papers/tcga/TCGA*/*/*_annotated_methylation.bed.gz | wc -l. # should be 9812 files
# Verify that all generated files have the same length: 485569 lines
for f in /data/papers/tcga/TCGA*/*/*_annotated_methylation.bed.gz; do
    echo -n "$f : "
    zcat "$f" | wc -l
done > ./methylation/annotated_methylation_file_line_counts.txt

cat ./methylation/annotated_methylation_file_line_counts.txt | grep 485569 | wc -l # should be 9812
rm ./methylation/annotated_methylation_file_line_counts.txt

###############################################################################
# SECTION 6 — STATISTICS AND GLOBAL METHYLATION DISTRIBUTION (UMR / LMR / FMR) 
###############################################################################
# ON FILTERED METHYLATION FILES
# to look into the range oof the beta values to make sure it between 0 and 100 for one file
zcat ./methylation/filtered_methylation/HM450_TCGA-ACC-TCGA-OR-A5KX-01A_1_annotated_methylation_filtered.bed.gz | awk '{print $5}' | sort -n | awk 'NR==1{min=$1} END{print "min:", min, "max:", $1}'
# min: 0 max: 99.3833

###############################################################################
# NUMBER OF PROBES PER CANCER TYPE AND NUMBER OF HEALTHY VS CANCER SAMPLES
###############################################################################
# number of methylation probes per cancer type --> should be the same for all samples : 470851 probes
mkdir -p ./methylation/methylation_counts/
for file in ./methylation/filtered_methylation/HM450*_annotated_methylation_filtered.bed.gz; do
    sample=$(basename "$file" _annotated_methylation_filtered.bed.gz)
    count=$(zcat "$file" | wc -l)
    echo -e "$sample\t$count" >> ./methylation/methylation_counts/methylation_probes_per_sample.tsv
done

# number of healthy vs cancer samples per cancer type 
# Comparing the methylation distribution for the same sample between healthy and cancer samples
# 2>/dev/null = hide all error messages from this command.

# to find samples with no healthy data any type : 0XA --> tumor ; 1XA --> healthy
mkdir -p ./methylation/methylation_counts
out=./methylation/methylation_counts/methylation_presence.tsv

echo -e "cancer\ttumor_count(01-09A)\thealthy_count(10-19A)" > "$out"

for cancer in $(find ./methylation/filtered_methylation \
  -name "HM450_TCGA-*_annotated_methylation_filtered.bed.gz" \
  | sed -E 's|.*/HM450_TCGA-([^-]+)-.*|\1|' \
  | sort -u); do

  tumor_count=$(find ./methylation/filtered_methylation \
    -name "HM450_TCGA-${cancer}-*-0[1-9]A*_annotated_methylation_filtered.bed.gz" \
    | wc -l)

  healthy_count=$(find ./methylation/filtered_methylation \
    -name "HM450_TCGA-${cancer}-*-1[0-9]A*_annotated_methylation_filtered.bed.gz" \
    | wc -l)

  printf "%s\t%d\t%d\n" "$cancer" "$tumor_count" "$healthy_count" >> "$out"
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
    title   = "Tumor (0XA) vs Healthy (1XA) Methylation Sample Counts per Cancer Type",
    x       = "Cancer Type",
    y       = "Number of Samples",
    caption = "Red cancer names = cancer types with no healthy samples") +
  theme_minimal(base_size = 13) +
  theme(
    plot.title   = element_text(hjust = 0.5),
    axis.text.x  = element_blank(),     # hide default x labels
    axis.ticks.x = element_blank(),
    plot.caption = element_text(hjust = 0.5, color = "red", size = 12),
    plot.margin  = margin(t = 10, r = 20, b = 60, l = 20) ) +  # more space for labels
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

# Statistics methylation probes per sample : 
#!!! DONE !!!
###############################################################################
# HISTOGRAM OF METHYLATION PERCENTAGE FOR ALL SAMPLES
###############################################################################
# Plot histogram of methylation percentage for all samples and highlight UMR, LMR, FMR regions with different colors and add legend with percentages of each region.
Rscript -e '
# One PDF with all samples, one per page — same plot logic

out_pdf <- "./results/summarys/Methylation_distribution_UMR_LMR_FMR_ALL_samples_histogram.pdf"
dir.create(dirname(out_pdf), recursive=TRUE, showWarnings=FALSE)

files <- list.files("./methylation/filtered_methylation/",
                    pattern="HM450_.*_annotated_methylation_filtered.bed.gz$",
                    recursive=TRUE,
                    full.names=TRUE)

pdf(out_pdf, width=11.7, height=8.3)
par(bg="white")

for (f in files) {

  sample <- sub("_annotated_methylation_filtered.bed.gz$", "", basename(f))
  cat("→ Processing", f, "/", length(files), ":", sample, "\n")

  data <- read.table(f)

  # -------- MINIMAL NA HANDLING (ADDED) --------
  meth_all <- as.numeric(data[, ncol(data)])   # last column = methylation %
  n_NA <- sum(is.na(meth_all))
  pct_NA <- round(100 * n_NA / length(meth_all), 1)

  meth <- meth_all[!is.na(meth_all)]           # remove NA only for calculations
  # --------------------------------------------

  # Compute frequencies for each region (NA-free)
  total <- length(meth)
  n_UMR <- sum(meth >= 0  & meth < 10)
  n_LMR <- sum(meth >= 10 & meth < 50)
  n_FMR <- sum(meth >= 50 & meth <= 100)

  pct_UMR <- round(100 * n_UMR / total, 1)
  pct_LMR <- round(100 * n_LMR / total, 1)
  pct_FMR <- round(100 * n_FMR / total, 1)

  # Base histogram (unchanged)
  hist(meth,
       breaks=50,
       col="white",
       border="black",
       main=paste0("Methylation distribution of ", sample, " sample"),
       xlab="Methylation (%)")

  # Colored background regions
  usr <- par("usr")

  rect(0,  usr[3], 10,  usr[4], col=rgb(0,0,1,0.35), border=NA)   # UMR
  rect(10, usr[3], 50,  usr[4], col=rgb(1,0,0,0.35), border=NA)   # LMR
  rect(50, usr[3], 100, usr[4], col=rgb(0,1,0,0.35), border=NA)   # FMR

  # Redraw histogram lines on top
  hist(meth,
       breaks=50,
       col=NA,
       border="black",
       add=TRUE)

  # Vertical separation lines
  abline(v = 10, col="black", lwd=2)
  abline(v = 50, col="black", lwd=2)

  # -------- LEGEND (NA % ADDED) --------
  legend("topright",
         legend=c(
           paste0("UMR: 0–10% (", pct_UMR, "%)"),
           paste0("LMR: 10–50% (", pct_LMR, "%)"),
           paste0("FMR: 50–100% (", pct_FMR, "%)"),
           paste0("NA: ", pct_NA, "%")
         ),
         fill=c(
           rgb(0,0,1,0.25),
           rgb(1,0,0,0.25),
           rgb(0,1,0,0.25),
           "white"
         ),
         border=c(NA, NA, NA, "black"),
         cex=1.2)
}

dev.off()
cat("Wrote:", out_pdf, "\n")
'

###############################################################################
# BARPLOT OF METHYLATION PERCENTAGE FOR ALL SAMPLES
###############################################################################
# Compare UMR, LMR, FMR percentages across all samples per cancer type using a barplot ONLY FOR CANCER SAMPLES (NOT HEALTHY)

# 1) Create a summary table with percentages of UMR, LMR, FMR per sample and cancer type
mkdir -p ./methylation/methylation_counts/
# This file will contain: sample, cancer_type, pct_UMR, pct_LMR, pct_FMR only on cancer samples (not healthy)
shopt -s nullglob
files=(./methylation/filtered_methylation/HM450_TCGA-*_annotated_methylation_filtered.bed.gz)
total_files=${#files[@]}
count=0

for file in "${files[@]}"; do
    count=$((count + 1))
    echo "→ Processing sample $count / $total_files : $(basename "$file")"
    sample=$(basename "$file" _annotated_methylation_filtered.bed.gz)
    cancer_type=$(echo "$sample" | cut -d'-' -f2)
    # Extract TCGA sample type (01=tumor, 11=normal, 06=metastasis…)
    type=$(echo "$sample" | cut -d'-' -f4 | cut -c1-2)
    # Skip healthy 
    [[ "$type" == "11" ]] && continue
    # Proceed only with TUMOR samples
    total=$(zcat "$file" | wc -l)
    n_UMR=$(zcat "$file" | awk '{print $NF}' | awk '$1 >= 0  && $1 < 10'  | wc -l)
    n_LMR=$(zcat "$file" | awk '{print $NF}' | awk '$1 >= 10 && $1 < 50'  | wc -l)
    n_FMR=$(zcat "$file" | awk '{print $NF}' | awk '$1 >= 50 && $1 <= 100' | wc -l)
    pct_UMR=$(awk -v n="$n_UMR" -v t="$total" 'BEGIN {printf "%.2f", (n / t) * 100}')
    pct_LMR=$(awk -v n="$n_LMR" -v t="$total" 'BEGIN {printf "%.2f", (n / t) * 100}')
    pct_FMR=$(awk -v n="$n_FMR" -v t="$total" 'BEGIN {printf "%.2f", (n / t) * 100}')
    echo -e "$sample\t$cancer_type\t$pct_UMR\t$pct_LMR\t$pct_FMR" >> ./methylation/methylation_counts/methylation_region_percentages_per_cancer_samples.tsv
done

# 2)add a header to the file
awk 'BEGIN{print "sample\tcancer_type\tpct_UMR\tpct_LMR\tpct_FMR"}{print}' ./methylation/methylation_counts/methylation_region_percentages_per_cancer_samples.tsv > ./methylation/methylation_counts/methylation_region_percentages_per_sample.tsv
rm ./methylation/methylation_counts/methylation_region_percentages_per_cancer_samples.tsv

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
# THE TOTAL FOR THE CPGS HAVE THE NUMBER OF CPG WITH NA VALUES INCLUDED SO IF WE SUM UMR LMR AND FMR IT DOESNT COMBINE TO 100%.
# 4) Plot barplot using ggplot2 of methylation_region_percentages
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

###############################################################################
# BOXPLOT OF METHYLATION PERCENTAGE FOR ALL SAMPLES
###############################################################################
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

###############################################################################
# SECTION 7 — DISTANCE OF CpG PROBES TO TSS
###############################################################################

# Distance of methylation sites to TSS : generate a file per sample with the distances of each methylation to the nearest TSS
# Since its the same probes for all samples we can do it for one sample only and reuse the file for all samples
mkdir -p ./methylation/dist_to_tss/
file='./methylation/filtered_methylation/HM450_TCGA-GBM-TCGA-OX-A56R-01A_1_annotated_methylation_filtered.bed.gz'
closestBed -t first -d -a "$file" -b /data/genome/annotations/hg38_tss.bed | awk '{print $NF}' > "./methylation/dist_to_tss/probes_dist_tss.txt"

# Plot histogram of distances to TSS for one sample 

Rscript -e '
dist_file <- "./methylation/dist_to_tss/probes_dist_tss.txt"
distances <- read.table(dist_file)[,1]

pdf("./results/summarys/dist_tss_methylation_histogram.pdf", width=10, height=10)
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
     main = paste0("All comon CpGs", "  | Proximal = ", pct_prox, "% Distal = ", pct_dist, "%"),
     breaks = 50,
     col = "steelblue",
     border = "black")
  abline(v = log10(2000), col = "red", lwd = 2)
dev.off()
'


###############################################################################
# SECTION 8 — DELTA METHYLATION (ALL CpGs after filtering))
###############################################################################

###########################################################################################
# SMOOTH SCATTER PLOT FOR ALL CANCER-HEALTHY SAMPLE PAIRS FILTERED
###########################################################################################
# plot the methylation distribution for both healthy and cancer samples of all samples using a delta methylation histogram :
#This script generates a single PDF where each page contains three plots for one tumor–normal pair: global CpG methylation, CpGs in NRF1 motifs, and CpGs in BANP motifs.
#In the second and third plots, CpGs overlapping NRF1 or BANP binding sites with strong methylation changes (|Δβ| ≥ 10%) are highlighted in red to assess motif-specific epigenetic disruption.
Rscript -e '
suppressPackageStartupMessages({
  library(data.table)
  library(KernSmooth)
})

# =========================
# INPUTS
# =========================
pairs_file <- "./methylation/sample_pairs_files/methylation_pairs_filtered.tsv"

motif_files <- list(
  NRF1 = "./methylation/overlaps/intersected_motifs_HM450/NRF1_intersected_methylation.bed",
  BANP = "./methylation/overlaps/intersected_motifs_HM450/BANP_intersected_methylation.bed"
)

out_pdf <- "./results/summarys/smoothScatter_ALL_vs_NRF1_vs_BANP_per_pair.pdf"
thr <- 10

dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)

# =========================
# LOAD DATA
# =========================
pairs <- fread(pairs_file)
stopifnot(all(c("cancer","patient","tumor_file","healthy_file") %in% colnames(pairs)))

motif_cpgs <- lapply(motif_files, function(f) {
  unique(fread(f, header = FALSE)[[4]])
})

# =========================
# START PDF
# =========================
pdf(out_pdf, width = 16, height = 6)
par(mfrow = c(1, 3), mar = c(4, 4, 4, 1))

# =========================
# LOOP OVER PAIRS
# =========================
for (i in seq_len(nrow(pairs))) {

  cat("Processing:", pairs$cancer[i], pairs$patient[i], "\n")

  tumor  <- fread(cmd = paste("zcat", pairs$tumor_file[i]),   header = FALSE)
  normal <- fread(cmd = paste("zcat", pairs$healthy_file[i]), header = FALSE)

  tumor_df  <- data.frame(probe = tumor[[4]],  meth_tumor  = tumor[[5]])
  normal_df <- data.frame(probe = normal[[4]], meth_normal = normal[[5]])

  merged <- merge(tumor_df, normal_df, by = "probe")
  if (nrow(merged) == 0) next

  delta <- merged$meth_tumor - merged$meth_normal

  # =========================
  # PLOT 1 — ALL CpGs
  # =========================
  smoothScatter(
    merged$meth_normal,
    merged$meth_tumor,
    xlab = "Healthy methylation (%)",
    ylab = "Tumor methylation (%)",
    main = paste0(
      pairs$cancer[i], " | ", pairs$patient[i], "\n",
      "All CpGs"
    )
  )
  abline(0, 1, col = "black", lwd = 2)
  abline(a = -thr, b = 1, col = "blue", lwd = 2)
  abline(a =  thr, b = 1, col = "blue", lwd = 2)

  # =========================
  # PLOT 2 — NRF1 motif CpGs (|Δ| ≥ 10%)
  # =========================
  smoothScatter(
    merged$meth_normal,
    merged$meth_tumor,
    xlab = "Healthy methylation (%)",
    ylab = "Tumor methylation (%)",
    main = "NRF1 motif CpGs in red (|DeltaBeta| ≥ 10%)"
  )

  is_NRF1 <- merged$probe %in% motif_cpgs$NRF1 &
             abs(delta) >= thr

  points(
    merged$meth_normal[is_NRF1],
    merged$meth_tumor[is_NRF1],
    pch = 16,
    cex = 0.5,
    col = rgb(1, 0, 0, 0.6)
  )

  abline(0, 1, col = "black", lwd = 2)
  abline(a = -thr, b = 1, col = "blue", lwd = 2)
  abline(a =  thr, b = 1, col = "blue", lwd = 2)

  # =========================
  # PLOT 3 — BANP motif CpGs (|Δ| ≥ 10%)
  # =========================
  smoothScatter(
    merged$meth_normal,
    merged$meth_tumor,
    xlab = "Healthy methylation (%)",
    ylab = "Tumor methylation (%)",
    main = "BANP motif CpGs in red (|DeltaBeta| ≥ 10%)"
  )

  is_BANP <- merged$probe %in% motif_cpgs$BANP &
             abs(delta) >= thr

  points(
    merged$meth_normal[is_BANP],
    merged$meth_tumor[is_BANP],
    pch = 16,
    cex = 0.5,
    col = rgb(1, 0, 0, 0.6)
  )

  abline(0, 1, col = "black", lwd = 2)
  abline(a = -thr, b = 1, col = "blue", lwd = 2)
  abline(a =  thr, b = 1, col = "blue", lwd = 2)
}

dev.off()
cat("DONE →", out_pdf, "\n")
'
#!!! Done !!!
###########################################################################################
# SMOOTH SCATTER PLOT FOR REPLICATES
###########################################################################################
# This code does the same as the one above just does it for pairs of tumor replicates not between healthy and cancer. 
Rscript -e '
suppressPackageStartupMessages({
  library(data.table)
  library(KernSmooth)
})

# =========================
# INPUTS
# =========================
pairs_file <- "./methylation/sample_pairs_files/methylation_replicates_pairs_filtered.tsv"

motif_files <- list(
  NRF1 = "./methylation/overlaps/intersected_motifs_HM450/NRF1_intersected_methylation.bed",
  BANP = "./methylation/overlaps/intersected_motifs_HM450/BANP_intersected_methylation.bed"
)

out_pdf <- "./results/summarys/smoothScatter_ALL_vs_NRF1_vs_BANP_per_pair_replicates_cancer.pdf"
thr <- 10

dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)

# =========================
# LOAD PAIRS
# =========================
pairs <- fread(pairs_file)
stopifnot(all(c("cancer","patient","replicate_1","replicate_2") %in% colnames(pairs)))

# =========================
# LOAD MOTIF CpGs
# =========================
motif_cpgs <- lapply(motif_files, function(f) {
  unique(fread(f, header = FALSE)[[4]])
})

# =========================
# START PDF
# =========================
pdf(out_pdf, width = 16, height = 6)
par(mfrow = c(1, 3), mar = c(4, 4, 4, 1))

# =========================
# LOOP OVER REPLICATE PAIRS
# =========================
for (i in seq_len(nrow(pairs))) {

  cat("Processing:", pairs$cancer[i], pairs$patient[i], "\n")

  # -------------------------
  # READ REPLICATES
  # -------------------------
  rep1 <- fread(cmd = paste("zcat", pairs$replicate_1[i]), header = FALSE)
  rep2 <- fread(cmd = paste("zcat", pairs$replicate_2[i]), header = FALSE)

  # -------------------------
  # EXTRACT PROBE + BETA
  # -------------------------
  rep1_df <- data.frame(
    probe      = rep1[[4]],
    meth_rep1  = rep1[[5]]
  )

  rep2_df <- data.frame(
    probe      = rep2[[4]],
    meth_rep2  = rep2[[5]]
  )

  merged <- merge(rep1_df, rep2_df, by = "probe")
  if (nrow(merged) == 0) next

  # -------------------------
  # DROP CpGs WITH ANY NA
  # -------------------------
  merged <- merged[
    !is.na(merged$meth_rep1) &
    !is.na(merged$meth_rep2),
  ]

  if (nrow(merged) == 0) {
    plot.new()
    title(main = paste0(
      pairs$cancer[i], " | ", pairs$patient[i],
      "\nNo CpGs with beta in both replicates"
    ))
    next
  }

  # -------------------------
  # DELTA
  # -------------------------
  delta <- merged$meth_rep1 - merged$meth_rep2

  # =========================
  # PLOT 1 — ALL CpGs
  # =========================
  smoothScatter(
    merged$meth_rep2,
    merged$meth_rep1,
    xlab = "Replicate 2 methylation (%)",
    ylab = "Replicate 1 methylation (%)",
    main = paste0(
      pairs$cancer[i], " | ", pairs$patient[i], "\n",
      "All CpGs"
    )
  )

  abline(0, 1, col = "black", lwd = 2)
  abline(a = -thr, b = 1, col = "blue", lwd = 2)
  abline(a =  thr, b = 1, col = "blue", lwd = 2)

  # =========================
  # PLOT 2 — NRF1 motif CpGs
  # =========================
  smoothScatter(
    merged$meth_rep2,
    merged$meth_rep1,
    xlab = "Replicate 2 methylation (%)",
    ylab = "Replicate 1 methylation (%)",
    main = "NRF1 motif CpGs in red (|Deltabeta| ≥ 10%)"
  )

  is_NRF1 <- merged$probe %in% motif_cpgs$NRF1 &
             abs(delta) >= thr

  points(
    merged$meth_rep2[is_NRF1],
    merged$meth_rep1[is_NRF1],
    pch = 16,
    cex = 0.5,
    col = rgb(1, 0, 0, 0.6)
  )

  abline(0, 1, col = "black", lwd = 2)
  abline(a = -thr, b = 1, col = "blue", lwd = 2)
  abline(a =  thr, b = 1, col = "blue", lwd = 2)

  # =========================
  # PLOT 3 — BANP motif CpGs
  # =========================
  smoothScatter(
    merged$meth_rep2,
    merged$meth_rep1,
    xlab = "Replicate 2 methylation (%)",
    ylab = "Replicate 1 methylation (%)",
    main = "BANP motif CpGs in red (|Deltabeta| ≥ 10%)"
  )

  is_BANP <- merged$probe %in% motif_cpgs$BANP &
             abs(delta) >= thr

  points(
    merged$meth_rep2[is_BANP],
    merged$meth_rep1[is_BANP],
    pch = 16,
    cex = 0.5,
    col = rgb(1, 0, 0, 0.6)
  )

  abline(0, 1, col = "black", lwd = 2)
  abline(a = -thr, b = 1, col = "blue", lwd = 2)
  abline(a =  thr, b = 1, col = "blue", lwd = 2)
}

dev.off()
cat("DONE →", out_pdf, "\n")
'
#!!!Done!!
###########################################################################################
# DELTA BOXPLOT PLOT FOR ALL CANCER-HEALTHY TYPES PAIRS FILTERED WITH ALL CPG AND CPG IN TF --> PER CANCER TYPE
###########################################################################################
mkdir -p ./methylation/overlaps/intersected_motifs_HM450/
# get the CpGs that are in NRF1 motifs and BANP motifs
bedtools intersect -u -a ./methylation/annotated_methylation_data_probes_filtered.bed -b ./motifs/NRF1_filtered_6mer.MA0506.3.bed > ./methylation/overlaps/intersected_motifs_HM450/NRF1_intersected_methylation.bed
bedtools intersect -u -a ./methylation/annotated_methylation_data_probes_filtered.bed -b ./motifs/BANP_filtered_6mer.MA2503.1.bed  > ./methylation/overlaps/intersected_motifs_HM450/BANP_intersected_methylation.bed
# PLOT THE BOXPLOT FOR MOTIFS TAKIG THE MEDIAN FOR EACH CANCER TYPE
# compute delta methylation per pair and classify into hypo hyper unchanged and make boxplot per cancer type for NRF1 and BANP motifs
Rscript -e '
library(data.table)
library(ggplot2)

# =========================
# INPUTS
# =========================
pairs_file <- "./methylation/sample_pairs_files/methylation_pairs_filtered.tsv"
out_pdf    <- "./results/summarys/stacked_median_composition_delta_methylation_per_cancer.pdf"
thr <- 10

# =========================
# HELPERS
# =========================
read_probe_meth <- function(f) {
  dt <- fread(cmd = paste("zcat", shQuote(f)), header = FALSE, select = c(4,5))
  setnames(dt, c("probe","meth"))
  dt[, meth := as.numeric(meth)]
  dt
}

# =========================
# LOAD PAIRS
# =========================
pairs <- fread(pairs_file)[!is.na(tumor_file) & !is.na(healthy_file)]

# =========================
# LOAD MOTIF CpGs
# =========================
motif_cpgs <- list(
  NRF1 = unique(fread(
    "./methylation/overlaps/intersected_motifs_HM450/NRF1_intersected_methylation.bed",
    header = FALSE
  )[[4]]),
  BANP = unique(fread(
    "./methylation/overlaps/intersected_motifs_HM450/BANP_intersected_methylation.bed",
    header = FALSE
  )[[4]])
)

# =========================
# COMPUTE PER-PATIENT COMPOSITIONS
# =========================
res <- list()

for (i in seq_len(nrow(pairs))) {
    cat(
      sprintf(
        "Processing %d / %d — cancer: %s\n",
        i, nrow(pairs), pairs$cancer[i]
      )
    )

  tumor  <- read_probe_meth(pairs$tumor_file[i])
  normal <- read_probe_meth(pairs$healthy_file[i])

  if (!identical(tumor$probe, normal$probe))
    stop("Probe order mismatch at pair ", i)

  delta <- tumor$meth - normal$meth

  class_idx <- list(
    Unchanged = delta >= -thr & delta <= thr,
    Hypomethylated = delta < -thr,
    Hypermethylated = delta > thr
  )

  for (cl in names(class_idx)) {

    idx <- class_idx[[cl]] & !is.na(delta)
    n_all <- sum(idx)
    if (n_all == 0) next

    n_nrf1 <- sum(idx & tumor$probe %chin% motif_cpgs$NRF1)
    n_banp <- sum(idx & tumor$probe %chin% motif_cpgs$BANP)
    n_other <- n_all - n_nrf1 - n_banp

    res[[length(res)+1]] <- data.table(
      cancer = pairs$cancer[i],
      class = cl,
      context = "NRF1",
      percentage = 100 * n_nrf1 / n_all
    )
    res[[length(res)+1]] <- data.table(
      cancer = pairs$cancer[i],
      class = cl,
      context = "BANP",
      percentage = 100 * n_banp / n_all
    )
    res[[length(res)+1]] <- data.table(
      cancer = pairs$cancer[i],
      class = cl,
      context = "Other",
      percentage = 100 * n_other / n_all
    )
  }
}

dt <- rbindlist(res)

# =========================
# MEDIAN PER CANCER
# =========================
dt_med <- dt[
  , .(percentage = median(percentage, na.rm = TRUE)),
  by = .(cancer, class, context)
]

dt_med[, class := factor(class,
  levels = c("Unchanged","Hypomethylated","Hypermethylated")
)]
dt_med[, context := factor(context,
  levels = c("NRF1","BANP","Other")
)]

# =========================
# COLORS
# =========================
cols <- c(
  NRF1  = "#1F78B4",
  BANP  = "#E31A1C",
  Other = "#BDBDBD"
)

# =========================
# PLOTTING
# =========================
pdf(out_pdf, width = 14, height = 5)

for (c in sort(unique(dt_med$cancer))) {

  sub <- dt_med[cancer == c]

  p <- ggplot(sub, aes(x = class, y = percentage, fill = context)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = cols) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(size = 11),
      plot.title = element_text(hjust = 0.5)
    ) +
    labs(
      title = paste0(
        "Median CpG context composition within Delta-methylation classes — ", c
      ),
      x = "",
      y = "Median percentage of CpGs (%)",
      fill = "CpG context"
    ) +
    geom_text(
      data = sub[context %in% c("NRF1","BANP")],
      aes(
        label = paste0(context, ": ", round(percentage, 1), "%"),
        y = cumsum(percentage) - percentage/2 +
            ifelse(context == "NRF1", 5, -3)
      ),
      color = "black",
      size = 3.5
    )
  print(p)
}
dev.off()
cat("Wrote:", out_pdf, "\n")
'

## !!! Done !!!
###############################################################################
# SECTION 9 — Distance OF CpG PROBES IN MOTIFS TO TSS
###############################################################################
# Look into the distance of CpG probes in motifs to othe TSS to see if they are proximal or distal
# produce a tsv file per sample with two columns: CpG probe ID and distance to nearest TSS in the folder ./methylation/dist_to_tss/
mkdir -p ./methylation/dist_to_tss/
for file in ./methylation/filtered_methylation/HM450*_annotated_methylation_filtered.bed.gz; do
  sample=$(basename "$file" .bed.gz)
  out="./methylation/dist_to_tss/${sample}_cpg_dist_tss.tsv"

  closestBed -t first -d \
    -a "$file" \
    -b /data/genome/annotations/hg38_tss.bed \
  | awk 'BEGIN{FS=OFS="\t"} {print $4, $NF}' \
  > "$out"

  echo "Wrote: $out"
done
##!! DONE !!
# for each tsv file in ./methylation/dist_to_tss/ filter only those in motif of NRF1 and make a histogram of distances to TSS for each sample
# NRF1 
mkdir -p ./results/summarys/methylation_motifs_TSS/NRF1
mkdir -p ./methylation/dist_to_tss_motif/NRF1

awk 'BEGIN{FS=OFS="\t"} {print $4}' \
  ./methylation/overlaps/intersected_motifs_HM450/NRF1_intersected_methylation.bed \
| grep -E '^cg' | sort -u \
> ./methylation/overlaps/intersected_motifs_HM450/NRF1_probes.txt

for f in ./methylation/dist_to_tss/*_cpg_dist_tss.tsv; do
  base=$(basename "$f")
  out="./methylation/dist_to_tss_motif/NRF1/${base/_cpg_dist_tss.tsv/_NRF1_dist_tss.txt}"

  awk 'NR==FNR{a[$1]=1; next} ($1 in a){print $2}' \
    ./methylation/overlaps/intersected_motifs_HM450/NRF1_probes.txt \
    "$f" > "$out"
done
# BANP 
mkdir -p ./results/summarys/methylation_motifs_TSS/BANP
mkdir -p ./methylation/dist_to_tss_motif/BANP

awk 'BEGIN{FS=OFS="\t"} {print $4}' \
  ./methylation/overlaps/intersected_motifs_HM450/BANP_intersected_methylation.bed \
| grep -E '^cg' | sort -u \
> ./methylation/overlaps/intersected_motifs_HM450/BANP_probes.txt

for f in ./methylation/dist_to_tss/*_cpg_dist_tss.tsv; do
  base=$(basename "$f")
  out="./methylation/dist_to_tss_motif/BANP/${base/_cpg_dist_tss.tsv/_BANP_dist_tss.txt}"
  awk 'NR==FNR{a[$1]=1; next} ($1 in a){print $2}' \
    ./methylation/overlaps/intersected_motifs_HM450/BANP_probes.txt \
    "$f" > "$out"
done
#!! DONE!!!
# Plot histogram of distances to TSS for NRF1 and BANP motifs per cancer type
Rscript -e '
library(data.table)

out_pdf <- "./results/summarys/methylation_motifs_TSS/dist_tss_methylation_NRF1_BANP_byCancer.pdf"
dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)

pdf(out_pdf, width = 11.7, height = 8.3)
par(bg = "white")

# =========================
# RUN ONCE FOR EACH TF
# =========================
for (TF in c("NRF1", "BANP")) {

  dist_dir <- paste0("./methylation/dist_to_tss_motif/", TF)
  out_tsv  <- paste0("./results/summarys/methylation_motifs_TSS/", TF,"/summary_", TF, "_byCancer.tsv")

  files <- list.files(dist_dir, pattern = paste0("_", TF, "_dist_tss.txt$"),full.names = TRUE)
  if (!length(files)) {
    plot.new()
    title(main = paste0(TF, ": no distance files found"))
    next
  }

  # extract cancer code from filename: TCGA-XXX-
  get_cancer <- function(x) {
    m <- regexpr("TCGA-[A-Z0-9]+-", x)
    if (m[1] == -1) return(NA_character_)
    sub("^TCGA-|-$", "", regmatches(x, m))
  }

  cancers <- vapply(basename(files), get_cancer, character(1))
  ok <- !is.na(cancers) & cancers != ""
  files <- files[ok]
  cancers <- cancers[ok]

  res <- list()
  uniq_cancers <- sort(unique(cancers))

  for (c in uniq_cancers) {
    fc <- files[cancers == c]

    d_list <- lapply(fc, function(f)
      suppressWarnings(as.numeric(readLines(f))))
    d <- unlist(d_list, use.names = FALSE)
    d <- d[is.finite(d)]

    n_samples <- length(fc)
    n_vals <- length(d)

    if (n_vals == 0) {
      plot.new()
      title(main = paste0("TCGA-", c, " | ", TF, " motif CpGs: none found"))
      res[[c]] <- data.frame(cancer=c, n_samples=n_samples, n_values=0,pct_prox=NA, pct_dist=NA, median=NA)
      next
    }

    proximal <- sum(d < 2000)
    distal   <- sum(d >= 2000)
    pct_prox <- round(100 * proximal / n_vals, 1)
    pct_dist <- round(100 * distal   / n_vals, 1)

    hist(log10(d + 1),
         xlim = c(0, 8),
         xlab = "Distance to TSS (log10)",
         main = paste0("TCGA-", c, " | ", TF,
                       " | samples=", n_samples,
                       " | n=", n_vals,
                       " | Proximal=", pct_prox,
                       "% | Distal=", pct_dist, "%"),
         breaks = 50,
         col = "steelblue",
         border = "black")
    abline(v = log10(2000), col = "red", lwd = 2)

    res[[c]] <- data.frame(cancer=c, n_samples=n_samples, n_values=n_vals,pct_prox=pct_prox, pct_dist=pct_dist, median=median(d))
  }

  dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)
  write.table(do.call(rbind, res), out_tsv,sep = "\t", quote = FALSE, row.names = FALSE)

}

dev.off()

cat("Wrote:", out_pdf, "\n")
'

# !!LAUNCHED RSCRIPT PCA SCRIPT !!!
###############################################################################
# SECTION 10 — PCA USING TOP 10000 VARIABLE CPGs PAN-CANCER
###############################################################################

###############################################################################
# GET TOP 10000 VARIABLE CPGs PAN-CANCER
###############################################################################
#  get list of top 10000 variable CpGs pan-cancer from yspill data 
# What he did to get top 10000 variable CpGs pan-cancer: /shared/misc/tinman/yspill/2018-08-10_pca_of_normals/pca_of_normals.R
# 1) For each TCGA cancer type, a pre-normalized DNA methylation beta-value matrix (CpGs × samples), generated using NOOB normalization, is loaded; each matrix contains methylation values for all samples of a given cancer type.
# 2) CpG probe identifiers (row names) are extracted from each cancer-specific matrix and intersected across all cancer types, retaining only CpGs that are present in every dataset, thereby enforcing a common and consistent CpG feature space across cancers.
# 3) Each cancer-specific matrix is subset to this shared CpG set and all matrices are concatenated column-wise to produce a single pan-cancer methylation matrix containing all samples from all cancer types with identical CpG ordering.
# 4) The variance of each CpG is computed across all samples in the merged pan-cancer matrix, CpGs are ranked by decreasing variance, and the 10,000 most variable CpGs are selected; the merged matrix is then restricted to these CpGs and saved for downstream dimensionality-reduction analyses such as PCA or SVD.

Rscript -e '
library(data.table)
load("/data/yspill/2018-08-10_pca_of_normals/TCGA-all_cancer_top10k_most_variable.RData")
writeLines(rownames(beta_var), "./cpg_list.txt")
'
# create a list of cancer types with methylation data available --> 33 
find ./methylation/filtered_methylation -maxdepth 1 -type f -name "*_annotated_methylation_filtered.bed.gz" -printf "%f\n" \
| awk '
{
  cancer=$0
  sub(/^HM450_TCGA-/, "", cancer)   # remove prefix
  sub(/-.*/, "", cancer)            # keep up to first dash -> cancer code
  print cancer
}
' | sort -u > ./methylation/cancer_types_with_methylation_data.txt

# create a list of cancer types with both normal and tumor samples for methylation analysis --> 23
find ./methylation/filtered_methylation -maxdepth 1 -type f -name "*_annotated_methylation_filtered.bed.gz" -printf "%f\n" \
| awk '
{
  cancer=$0
  sub(/^HM450_TCGA-/, "", cancer)   # remove prefix
  sub(/-.*/, "", cancer)            # keep up to first dash -> cancer code

  if ($0 ~ /-01A_/) t[cancer]=1
  if ($0 ~ /-11A_/) n[cancer]=1
}
END{
  for (c in t) if (n[c]) print c
}
' | sort > ./methylation/cancer_types_with_tumor_and_normal.txt


#!!! Done !!
# CREATE NOOB NORMALIZED MATRICES FOR METHYLATION DATA PER CANCER TYPE :
#This script takes **filtered HM450 methylation files** from the input directory `./methylation/filtered_methylation` and generates **one normalized methylation matrix per cancer type**, saving the results to the output directory `./methylation/noob_normalized_matrices`. For each cancer listed, it collects all corresponding methylation files, reads the CpG identifiers and beta values, and constructs a matrix where **rows correspond to CpG sites and columns to samples**, using a consistent CpG order across samples. Cancer types with fewer than two samples are skipped. Each per-cancer matrix (`myNorm`) is saved as an `.RData` file in the output directory for downstream pan-cancer analyses.
#!/usr/bin/env Rscript
mkdir -p ./methylation/noob_normalized_matrices
Rscript -e '
library(data.table)

meth_dir <- "./methylation/filtered_methylation"
out_dir  <- "./methylation/noob_normalized_matrices"

# list of cancer types
cancers <- scan("./methylation/cancer_types_with_methylation_data.txt",what = "", quiet = TRUE)
cancers <- cancers[cancers != ""]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (cancer in cancers) {
  message(">>> Processing ", cancer)

  files <- list.files(
    meth_dir,
    pattern = paste0("^HM450_TCGA-", cancer, "-.*_annotated_methylation_filtered.bed.gz$"),
    full.names = TRUE
  )

  if (length(files) < 2) {
    message("  Skipping (not enough samples)")
    next
  }

  # read first file to define CpG order
  dt0  <- fread(cmd = paste("zcat", shQuote(files[1])),header = FALSE, select = c(4, 5), showProgress = FALSE)
  cpgs <- as.character(dt0[[1]])

  myNorm <- matrix(
    NA_real_,
    nrow = length(cpgs),
    ncol = length(files),
    dimnames = list(cpgs, basename(files))
)

  myNorm[, 1] <- as.numeric(dt0[[2]])

  for (i in 2:length(files)) {
    dt <- fread(cmd = paste("zcat", shQuote(files[i])),header = FALSE, select = c(4, 5), showProgress = FALSE)
    myNorm[, i] <- as.numeric(dt[[2]])
  }

  save(myNorm,file = file.path(out_dir,paste0("TCGA-", cancer,"_cancer_normalization_noob.RData")))
}
'
#!! DONE !!
# Merges all per-cancer filtered myNorm_filt matrices into a single pan-cancer matrix with only common CpGs across all cancers.
# Saves the pan-cancer matrix to an RData file.
Rscript -e '
suppressPackageStartupMessages({
  library(purrr)
  library(abind)
})

# ------------------ inputs ------------------
cancers_file <- "./methylation/cancer_types_with_methylation_data.txt"
norm_dir     <- "./methylation/noob_normalized_matrices"
out_rdata    <- "./methylation/TCGA-all_cancer_pan_cancer_matrix.RData"

# ------------------ load cancer list ------------------
cancers <- scan(cancers_file, what = "", quiet = TRUE)
cancers <- cancers[cancers != ""]
cat(">>> Cancers:", length(cancers), "\n")

# ------------------ load all matrices (SEQUENTIAL) ------------------
cat(">>> Loading normalized matrices (sequential)...\n")

matrices <- lapply(cancers, function(cancer) {
  f <- file.path(norm_dir, paste0("TCGA-", cancer, "_cancer_normalization_noob.RData"))
  if (!file.exists(f)) return(NULL)
  load(f)  # loads myNorm
  myNorm
})
names(matrices) <- cancers

# ------------------ safety check ------------------
bad <- vapply(matrices, function(x) is.null(x) || !is.matrix(x) || is.null(rownames(x)), logical(1))
if (any(bad)) {
  cat(">>> ERROR: missing/bad matrices for cancers:\n")
  cat("    ", paste(names(matrices)[bad], collapse = ", "), "\n")
  stop("Fix missing/bad .RData files before continuing.")
}
cat(">>> Matrices loaded:", length(matrices), "\n")

# ------------------ common CpGs ------------------
cat(">>> Computing common CpGs...\n")
common_cpgs <- reduce(lapply(matrices, rownames), intersect)
cat(">>> Common CpGs:", length(common_cpgs), "\n")
if (length(common_cpgs) == 0) stop("No common CpGs found.")

# ------------------ subset & merge ------------------
cat(">>> Subsetting matrices to common CpGs...\n")
reduced_matrices <- lapply(matrices, function(x) x[common_cpgs, , drop = FALSE])
rm(matrices); gc()

cat(">>> Merging matrices (CpGs x samples)...\n")
pan_mat <- abind(reduced_matrices, along = 2)
rm(reduced_matrices); gc()

cat(">>> Pan-cancer matrix dims:", paste(dim(pan_mat), collapse = " x "), "\n")

# ------------------ save ------------------
cat(">>> Saving pan-cancer matrix to:", out_rdata, "\n")
save(pan_mat, file = out_rdata)

cat(">>> DONE\n")
'

#!! DONE !!
# Load the pan-cancer matrix 
library(matrixStats)

# load the full pan-cancer matrix
load("./methylation/TCGA-all_cancer_pan_cancer_matrix.RData")  # loads object: pan_mat

# select the 10,000 most variable CpGs
top_idx <- order(rowVars(pan_mat), decreasing = TRUE)[seq_len(10000)]
beta_var <- pan_mat[top_idx, , drop = FALSE]

# save results
save(beta_var, file = "./methylation/TCGA-all_cancer_top10k_most_variable.RData")
writeLines(rownames(beta_var), "./methylation/top10000_variable_CpGs_pan_cancer.txt")

# verify how many of the CpGs in cpg_list.txt are in the top 10000 variable CpGs pan-cancer filtered list
awk 'NR==FNR {a[$0]; next} ($0 in a) {c++} END{print c}' ./methylation/top10000_variable_CpGs_pan_cancer.txt ./10k_yspill_cpg_list.txt


###################################################################################

###############################################################################
# SNV SCRIPT
###############################################################################

###############################################################################
# 1) RAW SNV INVENTORY (UNFILTERED)
###############################################################################

# SNV statistics : number of SNVs per sample
mkdir -p ./snv/snv_counts/
counter=0
for file in ./snv/SNV_TCGA*vcf.gz; do
    counter=$((counter+1))
    sample=$(basename "$file" .vcf.gz | sed 's/SNV_//')
    count=$(zcat "$file" | grep -v "^#" | wc -l)
    echo -e "$counter\t$sample\t$count" >> ./snv/snv_counts/snv_per_sample.tsv
    echo "Processed $counter : $sample with $count SNVs"
done

# Plot histogram of SNVs per sample 
Rscript -e '
data <- read.table("./snv/snv_counts/snv_per_sample.tsv", header=FALSE)
snv_counts <- data[,3]   # 3rd column = SNV counts
pdf("./results/summarys/snv_per_sample_hist.pdf", width=10, height=8)
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

# Range of mutations for each cancer type : cancer_type min_mutation_length max_mutation_length
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


##############################################################################
# UNIQUE SNVs PER CANCER TYPE (RAW)
###############################################################################

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

pdf("./results/summarys/snv_per_cancer_type_barplot.pdf",width=12, height=9, bg = "white")

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

###############################################################################
# VARIANT SIZE CLASSES (RAW)
###############################################################################

# frequency of mutations found in cancer types in SNV small indels and structural variants 
# This script loops over each TCGA cancer type, scans all corresponding VCF files, and counts SNVs (1bp), small indels (<50 bp), and large (≥50 bp) indels.
# It runs the analysis in parallel (up to two cancers at a time) and writes one summary line per cancer to a TSV file.
# output file: snv_type_freq.tsv with columns: cancer, SNV, small_indel, struct_indel

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
ggsave("./results/summarys/stacked_mutation_types_log10.pdf",p,width = 12,height = 7,dpi = 300,bg = "white")'

# Compare the variation distribution between samples of cancer types using boxplot for SNV, small indels, structural indels to visualize each point to see if we have outliers 
# output file: snv_type_per_sample.tsv with columns: sample, cancer_type, SNV, small_indel, struct_indel
mkdir -p ./snv/snv_counts_per_sample

echo -e "sample\tcancer_type\tSNV\tsmall_indel\tstruct_indel" > ./snv/snv_counts_per_sample/snv_type_per_sample.tsv
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

# compute the boxplot (variants in %) of SNV, small indels, structural indels per sample across cancer types
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

###############################################################################
# 2) FILTER STRUCTURAL VARIANTS
###############################################################################

# Since structural variants (indels >= 50 bp) are less relevant for TF binding site analysis, we filter them out.
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
      if (d < 50) print
    }
  ' | gzip > "$out"

  echo "Filtered $(basename "$vcf")"
done

###############################################################################
# SNVs PER SAMPLE AFTER FILTERING
###############################################################################
# number of SNVs  per sample after
counter=0
for file in ./snv/snv_filtered_without_structural_variants/SNV_TCGA*vcf.gz; do
    counter=$((counter+1))   # increment counter
    sample=$(basename "$file" .vcf.gz | sed 's/SNV_//') # remove SNV_ prefix
    count=$(zcat "$file" | grep -v "^#" | wc -l) # count non-header lines --> number of SNVs
    echo -e "$counter\t$sample\t$count" >> ./snv/snv_counts/snv_filtered_per_sample.tsv
    echo "Processed $counter : $sample with $count SNVs"
done

###############################################################################
# Frequency of variants found within each cancer type after filtering structural variants
###############################################################################
#!!!Done!!
# frequency of variants found within each cancer type : number of samples that share each variant within a cancer type
mkdir -p ./snv/within_cancer_freq

for cancer in $(ls ./snv/snv_filtered_without_structural_variants/SNV_TCGA-*.vcf.gz | cut -d'-' -f2 | sort -u); do
  nsamples=$(ls ./snv/snv_filtered_without_structural_variants/SNV_TCGA-${cancer}-*.vcf.gz | wc -l)
  echo "→ $cancer ($nsamples samples)"

  for vcf in ./snv/snv_filtered_without_structural_variants/SNV_TCGA-${cancer}-*.vcf.gz; do
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
#!!! Done !!!
###############################################################################
# 3) DISTANCE TO TSS (FILTERED)
###############################################################################
# Plot histogram of SNVs per sample after filtering structural variants
# Distance of SNVs sites to TSS : generate a file per sample with the distances of each SNV to the nearest TSS
# Convert VCF --> BED-like (0-based start, 1-based end)--> $2-1 and sort before using closestBed to find nearest TSS
mkdir -p ./snv/dist_to_tss

for vcf in ./snv/snv_filtered_without_structural_variants/SNV_TCGA-*_vs_*.vcf.gz; do
  [ -e "$vcf" ] || continue

  base=$(basename "$vcf" .vcf.gz)
  out="./snv/dist_to_tss/${base}_dist_tss.txt"

  # Skip if output already exists and is non-empty
  if [ -s "$out" ]; then
    echo "Skipping (already done): $base"
    continue
  fi

  echo "Processing: $base"

  zcat "$vcf" \
    | grep -v '^#' \
    | awk '{print $1"\t"$2-1"\t"$2}' \
    | sort -k1,1 -k2,2n \
    | closestBed -t first -d -a - -b /data/genome/annotations/hg38_tss.bed \
    | awk '{print $NF}' \
    > "$out"
done
#!!!Done!!!
# Plot histogram of distances to TSS for all samples in one pdf 
Rscript -e '
dist_files <- list.files(
  "./snv/dist_to_tss/",
  pattern = "_dist_tss.txt$",
  full.names = TRUE
)

out_pdf <- "./results/summarys/dist_tss_snv_histograms.pdf"

pdf(out_pdf, width = 10, height = 8)
par(bg = "white")

for (dist_file in dist_files) {

  distances <- read.table(dist_file)[,1]

  if (length(distances) == 0) next

  proximal <- sum(distances < 2000)
  distal   <- sum(distances >= 2000)
  total    <- length(distances)

  pct_prox <- round((proximal / total) * 100, 1)
  pct_dist <- round((distal   / total) * 100, 1)

  sample_name <- gsub("_dist_tss.txt$", "", basename(dist_file))

  hist(
    log10(distances + 1),
    xlim   = c(0, 8),
    breaks = 50,
    col    = "steelblue",
    border = "black",
    xlab   = "Distance to TSS (log10)",
    main   = paste0(
      sample_name,
      "  |  Proximal = ", pct_prox,
      "%  |  Distal = ", pct_dist, "%"
    )
  )

  abline(v = log10(2000), col = "red", lwd = 2)
}

dev.off()
cat("DONE →", out_pdf, "\n")
'

###############################################################################
# 4) UNIQUE FILTERED SNVs ACROSS ALL CANCERS
###############################################################################
#!!!DONE!!!
# Create a unique list of SNVs across all cancer types (after filtering structural variants)
zcat ./snv/snv_filtered_without_structural_variants/*.vcf.gz \
  | awk 'BEGIN{OFS="\t"} !/^#/ {
      snv_id = $1 ":" $2 ":" $4 ":" $5
      print $1, $2-1, $2, snv_id
  }' \
  | sort -k1,1 -k2,2n -u \
  > ./snv/all_unique_SNVs_across_cancers.bed
#DONE
###############################################################################
# 5) SNV ∩ ATAC PEAKS
###############################################################################
# Overlap SNV positions with ATAC peak positions per sample and count the number of SNVs that overlap peaks per sample
mkdir -p ./snv/overlaps/snv_peak_overlap/

for snv_file in ./snv/snv_filtered_without_structural_variants/SNV_TCGA-*.vcf.gz; do
    snv_base=$(basename "$snv_file" .vcf.gz | sed 's/^SNV_//')
  
    cancer=$(echo "$snv_base" | cut -d'-' -f1-2)          # TCGA-ACC
    rest=${snv_base#${cancer}-}                           # TCGA-OR-A5J2-01A_vs_TCGA-OR-A5J2-10A_1
    patient=$(echo "$rest" | cut -d'-' -f1-3)             # TCGA-OR-A5J2

    peak_file="./peaks/filtered_peaks/ATAC_${cancer}_${patient}_1_peaks_macs.bed"

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

# extract sample names that actually have ATAC since only 217 ATAC peaks files exist 
samples_with_peaks <- sub("_snv_in_peaks.txt$", "", basename(overlap_files))
# basename() removes all folder paths → keeps only the filename
# sub() removes the suffix _snv_in_peaks.txt

# 2) Total SNVs for those samples only 
snv_counts <- read.table("./snv/snv_counts/snv_filtered_per_sample.tsv",header = FALSE, stringsAsFactors = FALSE)
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
###############################################################################
# 6) SNV ∩ METHYLATION
###############################################################################
#!!DONE!!!
# Overlap filtered SNVs with methylation sites per sample and count the number of SNVs that overlap methylation sites per sample
mkdir -p ./snv/overlaps/snv_methylation_overlap/

for snv_file in ./snv/snv_filtered_without_structural_variants/SNV_TCGA-*_vs_*_1.vcf.gz; do
    snv_base=$(basename "$snv_file" .vcf.gz)

    core=${snv_base#SNV_}

    # cancer = ACC
    cancer=$(echo "$core" | cut -d'-' -f2)

    # patient = TCGA-OR-A5J2
    patient=$(echo "$core" | cut -d'-' -f3-5)

    # tumor code = 01A
    tumor_code=$(echo "$core" | sed -E 's/.*-([0-9]{2}[A-Z]).*/\1/')

    methylation_file="./methylation/filtered_methylation/HM450_TCGA-${cancer}-${patient}-${tumor_code}_1_annotated_methylation_filtered.bed.gz"

    if [ -f "$methylation_file" ]; then
        echo "Processing $snv_base with $methylation_file"

        zcat "$snv_file" \
          | grep -v '^#' \
          | awk '{print $1"\t"$2-1"\t"$2}' \
          | sort -k1,1 -k2,2n \
          > "./snv/overlaps/snv_methylation_overlap/${snv_base}_snv.bed"

        bedtools intersect -u \
          -a "./snv/overlaps/snv_methylation_overlap/${snv_base}_snv.bed" \
          -b "$methylation_file" \
          | wc -l \
          > "./snv/overlaps/snv_methylation_overlap/${snv_base}_snv_in_methylation.txt"

        rm "./snv/overlaps/snv_methylation_overlap/${snv_base}_snv.bed"
        echo "Done $snv_base"
    else
        echo "No methylation file for $snv_base (expected $methylation_file)"
    fi
done
#!!Done!!!
# plot barplot of number of filtered SNV in methylation data
Rscript -e '
library(ggplot2)
# 1) Files with SNVs-in-methylation counts
overlap_files <- list.files("./snv/overlaps/snv_methylation_overlap",pattern = "_snv_in_methylation.txt$",full.names = TRUE)

samples_with_methylation <- basename(overlap_files)
samples_with_methylation <- sub("_snv_in_methylation.txt$", "", samples_with_methylation)
samples_with_methylation <- sub("^SNV_", "", samples_with_methylation)  # <-- IMPORTANT

# 2) Total SNVs for those samples only
snv_counts <- read.table("./snv/snv_counts/snv_filtered_per_sample.tsv",header = FALSE, stringsAsFactors = FALSE)
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
###############################################################################
# 7) MOTIFS (NRF1 / BANP)
###############################################################################
# kept means kept the motifs after the intersection.
# !!!Done!!!
# Overlap filtered SNV in TF motifs and ATAC peaks 
mkdir -p ./snv/overlaps/snv_motifs_overlap/
motifs_in_peaks_NRF1=./motifs/overlaps/motif_peak_overlaps/NRF1_filtered_6mer.MA0506.3_peak_overlaps.bed
motifs_in_peaks_BANP=./motifs/overlaps/motif_peak_overlaps/BANP_filtered_6mer.MA2503.1_peak_overlaps.bed
bedtools intersect -u -a ./snv/all_unique_SNVs_across_cancers.bed -b $motifs_in_peaks_NRF1 > ./snv/overlaps/snv_motifs_overlap/NRF1_SNVs_in_motifs_in_peaks.bed
bedtools intersect -u -a ./snv/all_unique_SNVs_across_cancers.bed -b $motifs_in_peaks_BANP > ./snv/overlaps/snv_motifs_overlap/BANP_SNVs_in_motifs_in_peaks.bed

# Overlap filtered SNV in TF motifs only
bedtools intersect -u -a ./snv/all_unique_SNVs_across_cancers.bed -b ./motifs/NRF1_filtered_6mer.MA0506.3.bed > ./snv/overlaps/intersected_motifs_snv_filtered/NRF1_intersected_SNVs.bed 
bedtools intersect -u -a ./snv/all_unique_SNVs_across_cancers.bed -b ./motifs/BANP_filtered_6mer.MA2503.1.bed > ./snv/overlaps/intersected_motifs_snv_filtered/BANP_intersected_SNVs.bed
# !! DONE !!

# overlap with NRF1 motifs 
bedtools intersect -u \
  -a ./motifs/NRF1_filtered_6mer.MA0506.3.bed \
  -b ./snv/all_unique_SNVs_across_cancers.bed \
  > ./snv/overlaps/intersected_motifs_snv_filtered/NRF1_filtered_SNVs_in_motifskept.bed

# overlap with BANP motifs 
bedtools intersect -u \
  -a ./motifs/BANP_filtered_6mer.MA2503.1.bed \
  -b ./snv/all_unique_SNVs_across_cancers.bed \
  > ./snv/overlaps/intersected_motifs_snv_filtered/BANP_filtered_SNVs_in_motifskept.bed

# overlap with NRF1 motifs in peaks 
bedtools intersect -u \
  -a ./motifs/overlaps/motif_peak_overlaps/NRF1_filtered_6mer.MA0506.3_peak_overlaps.bed \
  -b ./snv/all_unique_SNVs_across_cancers.bed \
  > ./snv/overlaps/intersected_motifs_snv_filtered/NRF1_filtered_SNVs_in_motifs_in_peaks_kept.bed

# overlap with BANP motifs in peaks
bedtools intersect -u \
  -a ./motifs/overlaps/motif_peak_overlaps/BANP_filtered_6mer.MA2503.1_peak_overlaps.bed \
  -b ./snv/all_unique_SNVs_across_cancers.bed \
  > ./snv/overlaps/intersected_motifs_snv_filtered/BANP_filtered_SNVs_in_motifs_in_peaks_kept.bed

###############################################################################
# 8) SNV ∩ PEAKS ∩ METHYLATION
###############################################################################

# Counting cancer types with both peaks and SNV data
mkdir -p ./results/summarys

# List of cancer types with methylation data (from the methylation file names) --> 33 cancer types
ls ./methylation/filtered_methylation/ | grep -- 'annotated_methylation_filtered.bed.gz' | cut -d'-' -f2 | sort | uniq > ./methylation/cancer_types_with_methylation.txt

# List of cancer types with ATAC peaks (cancer_counts.tsv) --> 22 cancer types
cut -f1 ./peaks/cancer_counts.tsv | sort > ./peaks/cancer_types_with_peaks.txt

# List of cancer types with SNV data (snv_filtered_total_per_cancer.tsv) --> 32 cancer types
cut -f2,3 ./snv/snv_counts/snv_filtered_per_sample.tsv | awk -F'\t' '{
   sample=$1;           # column 1 = sample name
   snv=$2;              # column 2 = SNV count
   split(sample,a,"-"); # split sample name using '-'
   cancer=a[2];         # the 2nd element is the cancer type
   print cancer "\t" snv;}' \
   | awk '{sum[$1] += $2} END {for (c in sum) print c "\t" sum[c]}' | sort -r -k1,1  > ./snv/snv_counts_per_cancer_type/snv_filtered_total_per_cancer.tsv
cut -f1 ./snv/snv_counts_per_cancer_type/snv_filtered_total_per_cancer.tsv | sort > ./snv/cancer_types_with_snv_filtered.txt

# ============================================================
# Heatmap of cancer types SNV ∩ PEAKS ∩ METHYLATION
# ============================================================

# heatmap of cancer types with ATAC vs SNV vs Methylation data
Rscript -e 'library(ggplot2)

# Read lists (one cancer type per line)
atac <- scan("./peaks/cancer_types_with_peaks.txt", what="")
snv  <- scan("./snv/cancer_types_with_snv_filtered.txt", what="")
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

ggsave("./results/summarys/cancer_types_ATAC_SNV_Methylation_heatmap.pdf",p, width = 10, height = 10, dpi = 300, bg = "white")
'
# clean up intermediate files
rm ./snv/cancer_types_with_snv.txt
rm ./methylation/cancer_types_with_methylation.txt
rm ./peaks/cancer_types_with_peaks.txt

###############################################################################
# 9) SNV HEATMAP ACROSS CANCER TYPES
###############################################################################
##!!! Done !!!
# Heatmap to show SNVs across differents cancer types 
# 1)Build a list of unique SNVs per cancer type
mkdir -p ./snv/snv_filtered_unique_per_cancer/

# extract cancer types from FILTERED filenames
cancers=$(find ./snv/snv_filtered_without_structural_variants \
  -name 'SNV_TCGA-*.vcf.gz' \
  | sed 's#.*/SNV_TCGA-##; s/-.*//' \
  | sort -u)

for cancer in $cancers; do
    echo "→ Processing $cancer"

    zgrep -hv "^#" ./snv/snv_filtered_without_structural_variants/SNV_TCGA-"$cancer"-*.vcf.gz \
        | awk '{print $1":"$2":"$4":"$5}' \
        | sort -u \
        | gzip > ./snv/snv_filtered_unique_per_cancer/${cancer}_unique_SNVs.txt.gz

    echo "Done $cancer"
done
# !! Done!!  
# 2) Build a presence/absence matrix of SNVs (rows = SNVs, columns = cancer types)
outfile="./snv/snv_presence_matrix.tsv"

# header : SNV    ACC  BLCA  BRCA  COAD ...
echo -ne "SNV" > $outfile # first column name of the header is SNV
for f in ./snv/snv_filtered_unique_per_cancer/*_unique_SNVs.txt.gz; do
    cancer=$(basename "$f" _unique_SNVs.txt.gz)
    echo -ne "\t${cancer}" >> $outfile # add a tab + the cancer name for each type of cancer
done
echo "" >> $outfile # add a newline at the end of the header

# 3) Keep only the motifs that have SNV inside (instead of keeping the SNV that overlap with motifs) and plot the number of motifs that have SNV inside for NRF1 and BANP
# transform the SNV data file to a bed file but intersecting with the motif bed keeping the motif that have SNV inside
#!!DONE!!
# count number of SNVs that overlap with NRF1 and BANP motifs
cat ./snv/overlaps/intersected_motifs_snv_filtered/NRF1_filtered_SNVs_in_motifskept.bed  | wc -l > ./snv/overlaps/intersected_motifs_snv_filtered/NRF1_intersected_SNVs_count_kept.tsv 
cat ./snv/overlaps/intersected_motifs_snv_filtered/BANP_filtered_SNVs_in_motifskept.bed  | wc -l > ./snv/overlaps/intersected_motifs_snv_filtered/BANP_intersected_SNVs_count_kept.tsv
# Plot barplot of SNV counts in NRF1 and BANP motifs
Rscript -e '
library(ggplot2)

# Motif names
motifs <- c("NRF1", "BANP")

# Read the numeric counts from the .tsv files (SNVs in motifs)
snv_counts <- c(
  as.numeric(readLines("./snv/overlaps/intersected_motifs_snv_filtered/NRF1_intersected_SNVs_count_kept.tsv")),
  as.numeric(readLines("./snv/overlaps/intersected_motifs_snv_filtered/BANP_intersected_SNVs_count_kept.tsv")))

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

# 4) Build presence/absence matrix of SNVs in NRF1 and BANP motifs per cancer type using these files ./snv/overlaps/intersected_motifs_snv_filtered/NRF1_intersected_SNVs.bed and ./snv/overlaps/intersected_motifs_snv_filtered/BANP_intersected_SNVs.bed
# transform the intersection of SNV in motifs bed file to a list of SNV ids only
# note that the SNV id is in the 4th column of the bed file and that it will generate a unique list of SNV ids only so less lines than the bed file
# presence and absence matrix separtly for each type of cancer type for NRF1 and BANP motifs
#!! DONE!!
# 1) make a list of SNV ids that overlap with NRF1 and BANP motifs
cut -f4 ./snv/overlaps/intersected_motifs_snv_filtered/NRF1_intersected_SNVs.bed | sort -u > ./snv/overlaps/intersected_motifs_snv_filtered/NRF1_intersected_SNVs_ids.txt
#!! DONE!!
snv_ids_NRF1="./snv/overlaps/intersected_motifs_snv_filtered/NRF1_intersected_SNVs_ids.txt"
# folder for intermediate columns
mkdir -p ./snv/snv_presence_columns_NRF1/

build_column_NRF1() {
    local snv_list="$1"      # BANP_intersected_SNVs_ids.txt
    local cancer_file="$2"   # ACC_unique_SNVs.txt.gz
    local out_file="$3"

    local cancer
    cancer=$(basename "$cancer_file" _unique_SNVs.txt.gz)
    echo "  → building column for $cancer"

    awk -v OFS='\t' '
        NR==FNR { has[$1]=1; next }     # FIRST = cancer SNVs
        {
            val = ($1 in has ? 1 : 0);  # check if NRF1 SNV exists in cancer
            print $1, val;
        }
    ' <(zcat "$cancer_file") "$snv_list" > "$out_file"
}

# 2) loop over cancers in parallel
threads=5    # threads minimum to keep free 
max_jobs=2   # max parallel jobs
job_count=0

for f in ./snv/snv_filtered_unique_per_cancer/*_unique_SNVs.txt.gz; do
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
echo "All jobs for BANP columns are done."
#!!DONE!!
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

#!!! Done !!
# Plot heatmap of raw shared SNVs between cancer types for TF motifs

# heatmap of raw shared SNVs between cancer types IN A LOG10 SCALE
# it also adds the clustering groups cancers that have similar patterns of shared NRF1-motif SNVs across all cancer types
Rscript -e '
library(pheatmap)

pdf("./results/summarys/SNV_shared_NRF1_BANP_variants_between_cancer_types_heatmap_rescaled.pdf",
    width = 10, height = 8, bg = "white")

# ===================== NRF1 =====================
data <- read.table("./snv/snv_presence_matrix_NRF1.tsv",
                   header = TRUE, stringsAsFactors = FALSE)
mat <- as.matrix(data[ , -1])
rownames(mat) <- data$SNV

shared <- t(mat) %*% mat
shared_log <- log10(shared + 1)
diag(shared_log) <- NA

min_val <- min(shared_log, na.rm = TRUE)
max_val <- max(shared_log, na.rm = TRUE)

pheatmap(shared_log,
         main = "Number of shared NRF1-motif variants\nbetween cancer types (log10 scale)",
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = colorRampPalette(c("white","steelblue","navy"))(200),
         breaks = seq(min_val, max_val, length.out = 201),
         na_col = "grey90")

# ===================== BANP =====================
data <- read.table("./snv/snv_presence_matrix_BANP.tsv",
                   header = TRUE, stringsAsFactors = FALSE)
mat <- as.matrix(data[ , -1])
rownames(mat) <- data$SNV

shared <- t(mat) %*% mat
shared_log <- log10(shared + 1)
diag(shared_log) <- NA

min_val <- min(shared_log, na.rm = TRUE)
max_val <- max(shared_log, na.rm = TRUE)

pheatmap(shared_log,
         main = "Number of shared BANP-motif variants\nbetween cancer types (log10 scale)",
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = colorRampPalette(c("white","steelblue","navy"))(200),
         breaks = seq(min_val, max_val, length.out = 201),
         na_col = "grey90")

dev.off()
'
#!!! Done !!
# Plot heatmap of percentage of shared SNVs between cancer types for BANP motifs
# Heatmap with the percentage oof shared variants between the different cancer types : 
# Heatmap with the percentage of shared variants between cancer types
Rscript -e '
library(pheatmap)

pdf("./results/summarys/SNV_shared_NRF1_BANP_percentage_between_cancer_types_heatmap.pdf",
    width = 10, height = 8, bg = "white")

# ===================== NRF1 =====================
data <- read.table("./snv/snv_presence_matrix_NRF1.tsv",
                   header = TRUE, stringsAsFactors = FALSE)

cancer_cols  <- colnames(data)[-1]
cancer_order <- sort(cancer_cols)

mat <- as.matrix(data[, cancer_order])
rownames(mat) <- data$SNV

shared <- t(mat) %*% mat

totals <- diag(shared)
pct_shared <- sweep(shared, 1, totals, FUN = "/") * 100
diag(pct_shared) <- NA

min_val <- min(pct_shared, na.rm = TRUE)
max_val <- max(pct_shared, na.rm = TRUE)
cat("NRF1: percentage range:", min_val, "to", max_val, "\n")

pheatmap(pct_shared,
         main = "Percentage of shared NRF1-motif variants\nbetween cancer types",
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = colorRampPalette(c("white","steelblue","navy"))(200),
         breaks = seq(min_val, max_val, length.out = 201),
         na_col = "grey90")

# ===================== BANP =====================
data <- read.table("./snv/snv_presence_matrix_BANP.tsv",
                   header = TRUE, stringsAsFactors = FALSE)

cancer_cols  <- colnames(data)[-1]
cancer_order <- sort(cancer_cols)

mat <- as.matrix(data[, cancer_order])
rownames(mat) <- data$SNV

shared <- t(mat) %*% mat

totals <- diag(shared)
pct_shared <- sweep(shared, 1, totals, FUN = "/") * 100
diag(pct_shared) <- NA

min_val <- min(pct_shared, na.rm = TRUE)
max_val <- max(pct_shared, na.rm = TRUE)
cat("BANP: percentage range:", min_val, "to", max_val, "\n")

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
#!!! NOT DONE !!!
###############################################################################
# 10) SNV IN KNOWN GENE NFX1 
###############################################################################

###############################################################################
# GETTING NFX1 GENOMIC COORDINATES
###############################################################################
# genes annotations found in  : /data/genome/annotations/hg38_genes.bed
mkdir ./snv/NFX1_gene/
cat /data/genome/annotations/hg38_genes.bed | grep 'ENSG00000086102' > ./snv/NFX1_gene/NFX1_genomic_coordinates.txt
# chr9    33290510        33371157        NFX1    +       ENSG00000086102

###############################################################################
# FINDING SNV IN BANP MOTIF THAT OVERLAP WITH THE NFX1 GENE 
###############################################################################
# looking into the file ./snv/overlaps/intersected_motifs_snv_filtered/BANP_intersected_SNVs.bed
cat ./snv/overlaps/intersected_motifs_snv_filtered/BANP_intersected_SNVs.bed | grep 'chr9' | awk '$2 > 33290510 && $3 < 33371157 {print $4}'
# chr9:33329363:C:T
# chr9:33329368:G:A
# chr9:33339406:T:G
# chr9:33339410:A:C
# chr9:33339411:C:T 

cat ./motifs/BANP_filtered_6mer.MA2503.1.bed | grep 'chr9' | grep 'TATCGCGAGA'
# chr9    20670569        20670579        TATCGCGAGA      5.01e-06        +
# chr9    33290683        33290693        TATCGCGAGA      5.01e-06        + <-- THIS SHOULD BE WHERE I HAVE THE SNV
# chr9    38354804        38354814        TATCGCGAGA      5.01e-06        +
# chr9    127451451       127451461       TATCGCGAGA      5.01e-06        +



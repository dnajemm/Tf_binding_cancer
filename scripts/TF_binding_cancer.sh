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

###############################################################################
# BANP 
###############################################################################

# Install MEM BANP of JASPAR 2026 and put it in the folder of MEME:

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
# NRF1 MOTIF CHANGE FORM JASPAR GCGCATGCGC
###############################################################################

###############################################################################
# NRF1 MEME GENERATION
###############################################################################
# Copy MEME to our folder from /data/bardeta/ralph_cancer/GCGCATGCGC.meme
cp /data/bardeta/ralph_cancer/GCGCATGCGC.meme /data/najemd/TF_binding_cancer/motifs/

###############################################################################
# NRF1 FIMO SCANNING
###############################################################################

# ===============================
# VARIABLES (ONLY THINGS TO SET)
# ===============================
motif=NRF1
genome=hg38
pval=0.001
GENOME_FA=/data/genome/genomes/${genome}.fa
MOTIF_MEME="/data/najemd/TF_binding_cancer/motifs/GCGCATGCGC.meme"

# ===============================
# STEP 1 — RUN FIMO
# ===============================
fimo --text --thresh ${pval} ${MOTIF_MEME} ${GENOME_FA} 2> /dev/null | awk '$1=="motif_id"{next}{ print $3, $4-1, $5, toupper($9), $8, $6}' \
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

echo "[DONE] NRF1 motif BED with sequence:"
echo "       ./motifs/${motif}.bed.gz"


###############################################################################
# 1) MOTIF FILTERING SYSTEMATICALLY
###############################################################################
motif="NRF1"
ref="GCGCATGCGC"   # length 10

zcat ./motifs/${motif}.bed.gz \
| awk -v REF="$ref" -v MOTIF="$motif" '
BEGIN{L=length(REF)}
{
  s=toupper($4)
  if(length(s)!=L) next

  # protect CG sites: positions 2-3 and 8-9 must match exactly
  if(substr(s,2,2) != substr(REF,2,2)) next
  if(substr(s,8,2) != substr(REF,8,2)) next

  # count mismatches only in flanks (positions 1,4-7,10)
  mm=0
  if(substr(s,1,1)  != substr(REF,1,1))  mm++
  for(i=4;i<=7;i++) if(substr(s,i,1) != substr(REF,i,1)) mm++
  if(substr(s,10,1) != substr(REF,10,1)) mm++

  if(mm>=0 && mm<=3) c[mm]++
}
END{
  for(k=0;k<=3;k++) printf("%s\tmm%d\t%d\n", MOTIF, k, c[k]+0)
  printf("%s\ttotal_mm0to3\t%d\n", MOTIF, (c[0]+c[1]+c[2]+c[3])+0)
}'


# NRF1    mm0     969
# NRF1    mm1     6084
# NRF1    mm2     27232
# NRF1    mm3     40804
# NRF1    total_mm0to3   75089

motif="NRF1"
ref="GCGCATGCGC"

zcat ./motifs/${motif}.bed.gz \
| awk -v REF="$ref" '
BEGIN{L=length(REF)}
{
  s=toupper($4)
  if(length(s)!=L) next
  mm=0
  for(i=1;i<=L;i++) if(substr(s,i,1)!=substr(REF,i,1)) mm++
  if(mm>=0 && mm<=3) c[mm]++
}
END{
  for(k=0;k<=3;k++) printf("%s\tmm%d\t%d\n","'"$motif"'",k,c[k]+0)
  printf("%s\ttotal_mm0to3\t%d\n","'"$motif"'",(c[0]+c[1]+c[2]+c[3])+0)
}'

# NRF1    mm0     969
# NRF1    mm1     12144
# NRF1    mm2     247229
# NRF1    mm3     100018
# NRF1    total_mm0to3    360360

motif="BANP"
ref="TCTCGCGAGA"   # length 10

zcat ./motifs/${motif}.bed.gz \
| awk -v REF="$ref" -v MOTIF="$motif" '
BEGIN{L=length(REF)}
{
  s=toupper($4)
  if(length(s)!=L) next

  # protect CGCG block (positions 4-7): must match exactly
  if(substr(s,4,4) != substr(REF,4,4)) next

  # count mismatches only in flanks (1-3 and 8-10)
  mm=0
  for(i=1;i<=3;i++)  if(substr(s,i,1)!=substr(REF,i,1)) mm++
  for(i=8;i<=10;i++) if(substr(s,i,1)!=substr(REF,i,1)) mm++

  if(mm>=0 && mm<=3) c[mm]++
}
END{
  for(k=0;k<=3;k++) printf("%s\tmm%d\t%d\n", MOTIF, k, c[k]+0)
  printf("%s\ttotal_mm0to3\t%d\n", MOTIF, (c[0]+c[1]+c[2]+c[3])+0)
}'

# BANP    mm0     309
# BANP    mm1     1218
# BANP    mm2     9441
# BANP    mm3     44619
# BANP    total_mm0to3   55587

motif="BANP"
ref="TCTCGCGAGA"

zcat ./motifs/${motif}.bed.gz \
| awk -v REF="$ref" '
BEGIN{L=length(REF)}
{
  s=toupper($4)
  if(length(s)!=L) next
  mm=0
  for(i=1;i<=L;i++) if(substr(s,i,1)!=substr(REF,i,1)) mm++
  if(mm>=0 && mm<=3) c[mm]++
}
END{
  for(k=0;k<=3;k++) printf("%s\tmm%d\t%d\n","'"$motif"'",k,c[k]+0)
  printf("%s\ttotal_mm0to3\t%d\n","'"$motif"'",(c[0]+c[1]+c[2]+c[3])+0)
}'

# BANP    mm0     309
# BANP    mm1     12694
# BANP    mm2     231331
# BANP    mm3     290316
# BANP    total_mm0to3    534650

##################################################
# FOUND THAT WE SHOULD ALLOW MAXIMUM 2 MISMATCHES 
##################################################
##################################################
# NOT IN CG SITES
##################################################
# Filter for NRF1 --> 34,285 occurrences left
motif="NRF1"
ref="GCGCATGCGC"   # length 10
in="./motifs/${motif}.bed.gz"
out="./motifs/${motif}_mm0to2_noCGmm.bed.gz"

zcat "$in" | awk -v OFS="\t" -v REF="$ref" -v MOTIF="$motif" '
BEGIN{L=length(REF)}
{
  seq=toupper($4)
  if(length(seq)!=L) next

  # protect CG sites: positions 2-3 and 8-9 must match exactly
  if(substr(seq,2,2) != substr(REF,2,2)) next
  if(substr(seq,8,2) != substr(REF,8,2)) next

  # count mismatches only in flanks (positions 1,4-7,10)
  mm = 0
  if(substr(seq,1,1)  != substr(REF,1,1))  mm++
  for(i=4;i<=7;i++) if(substr(seq,i,1) != substr(REF,i,1)) mm++
  if(substr(seq,10,1) != substr(REF,10,1)) mm++

  # allow max 2 mismatches in flanks
  if(mm <= 2){
    c[mm]++
    print
  }
}
END{
  for(k=0;k<=2;k++) printf("%s\tmm%d\t%d\n", MOTIF, k, c[k]+0) > "/dev/stderr"
  printf("%s\ttotal_allowed_mm0to2\t%d\n", MOTIF, (c[0]+c[1]+c[2])+0) > "/dev/stderr"
}' | gzip > "$out"

echo "[DONE] Wrote: $out"

# Filter for BANP --> 10,968 occurrences left
motif="BANP"
ref="TCTCGCGAGA"   # length 10
in="./motifs/${motif}.bed.gz"
out="./motifs/${motif}_mm0to2_noCGmm.bed.gz"

zcat "$in" | awk -v OFS="\t" -v REF="$ref" -v MOTIF="$motif" '
BEGIN{L=length(REF)}
{
  seq=toupper($4)
  if(length(seq)!=L) next

  # protect CGCG block (positions 4-7): must match exactly
  if(substr(seq,4,4) != substr(REF,4,4)) next

  # count mismatches only in flanks (1-3 and 8-10)
  mm=0
  for(i=1;i<=3;i++)  if(substr(seq,i,1)!=substr(REF,i,1)) mm++
  for(i=8;i<=10;i++) if(substr(seq,i,1)!=substr(REF,i,1)) mm++

  # allow max 2 mismatches in flanks
  if(mm <= 2){
    c[mm]++
    print
  }
}
END{
  for(k=0;k<=2;k++) printf("%s\tmm%d\t%d\n", MOTIF, k, c[k]+0) > "/dev/stderr"
  printf("%s\ttotal_allowed_mm0to2\t%d\n", MOTIF, (c[0]+c[1]+c[2])+0) > "/dev/stderr"
}' | gzip > "$out"

echo "[DONE] Wrote: $out"

##################################################
# IN CG SITES
##################################################
# NRF1 --> 260,342 
motif="NRF1"
ref="GCGCATGCGC"
in="./motifs/${motif}.bed.gz"
out="./motifs/${motif}_mm0to2.bed.gz"

zcat "$in" \
| awk -v OFS="\t" -v REF="$ref" '
BEGIN{L=length(REF)}
{
  s=toupper($4)
  if(length(s)!=L) next

  mm=0
  for(i=1;i<=L;i++) if(substr(s,i,1)!=substr(REF,i,1)) mm++

  if(mm <= 2) print
}' \
| gzip > "$out"

echo "[DONE] Wrote: $out"
# BANP --> 244,334
motif="BANP"
ref="TCTCGCGAGA"
in="./motifs/${motif}.bed.gz"
out="./motifs/${motif}_mm0to2.bed.gz"

zcat "$in" \
| awk -v OFS="\t" -v REF="$ref" '
BEGIN{L=length(REF)}
{
  s=toupper($4)
  if(length(s)!=L) next

  mm=0
  for(i=1;i<=L;i++) if(substr(s,i,1)!=substr(REF,i,1)) mm++

  if(mm <= 2) print
}' \
| gzip > "$out"

echo "[DONE] Wrote: $out"

# Allowing 2 mismatches not in a CG site :
# NRF1 : 34,285
# BANP : 10,968
# Allowing 2 mismatches including in a CG site :
# NRF1 : 260,342
# BANP : 244,334


# TO GENERATE A MEME LOGO USE THIS :
iupac2meme -logodds MOTIF > ./motifs/cg.meme

###############################################################################
# MOTIF COUNT BEFORE / AFTER FILTERING Allowing 2 mismatches
###############################################################################
# barplot of the number of motifs before and after filtering 
Rscript -e '
library(ggplot2)
library(scales)

motifs   <- c("NRF1", "BANP")
statuses <- c("Original", "mm0to2_noCGmm", "mm0to2")

files_map <- list(
  NRF1 = c("./motifs/NRF1.bed.gz",
           "./motifs/NRF1_mm0to2_noCGmm.bed.gz",
           "./motifs/NRF1_mm0to2.bed.gz"),
  BANP = c("./motifs/BANP.bed.gz",
           "./motifs/BANP_mm0to2_noCGmm.bed.gz",
           "./motifs/BANP_mm0to2.bed.gz")
)

nrf1_counts <- sapply(files_map[["NRF1"]], function(f) nrow(read.table(gzfile(f))))
banp_counts <- sapply(files_map[["BANP"]], function(f) nrow(read.table(gzfile(f))))

df <- data.frame(
  Motif  = rep(motifs, each = length(statuses)),
  Status = rep(statuses, times = length(motifs)),
  Count  = c(nrf1_counts, banp_counts)
)

df$Status <- factor(df$Status, levels = statuses)

# order the 3 bars within each Motif from bigger to smaller 
df$StatusWithin <- paste(df$Status, df$Motif, sep="__")
lev <- with(df, StatusWithin[order(Motif, -Count)])
df$StatusWithin <- factor(df$StatusWithin, levels = unique(lev))
out_pdf <- "./results/motifs/motif_counts_NRF1_BANP_barplot_mm0to2mm.pdf"
dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)

p <- ggplot(df, aes(x = Motif, y = Count, fill = Status, group = StatusWithin)) +
  geom_bar(stat = "identity", position = position_dodge2(width = 0.9, preserve = "single")) +
  geom_text(aes(label = comma(Count)),
            position = position_dodge2(width = 0.9, preserve = "single"),
            vjust = -0.3, size = 3) +
  scale_y_continuous(labels = comma) +
  theme_minimal() +
  labs(title = "Number of motifs across filtering mismatches",
       subtitle = "Original vs mm<=2 (with/without CG-protection)",
       y = "Count",
       x = "Motif") +
  theme(plot.title = element_text(hjust = 0.7, size = 12),
        plot.subtitle = element_text(hjust = 0.7, size = 9))

ggsave(out_pdf, plot = p, width = 7, height = 4, dpi = 300, bg = "white")
cat("[DONE] Wrote:", out_pdf, "\n")
'
###################################################################
# Filter out chrX and chrY motifs
###################################################################
for f in ./motifs/*_mm0to2.bed.gz; do
    base=$(basename "$f")                     
    sample=${base%_mm0to2.bed.gz} 
    echo "Filtering $sample ..."
    zcat "$f" | awk '$1 != "chrX" && $1 != "chrY"' | gzip > "./motifs/${sample}_mm0to2_noXY.bed.gz"
done

for f in ./motifs/*_noCGmm.bed.gz; do
    base=$(basename "$f")                     
    sample=${base%_mm0to2_noCGmm.bed.gz} 
    echo "Filtering $sample ..."
    zcat "$f" | awk '$1 != "chrX" && $1 != "chrY"' | gzip > "./motifs/${sample}_mm0to2_noCGmm_noXY.bed.gz"
done


zcat ./motifs/BANP_mm0to2_noCGmm.bed.gz | wc -l
# 10968
zcat ./motifs/BANP_mm0to2_noCGmm_noXY.bed.gz | wc -l
# 10418 (no chrX/Y motifs )
zcat ./motifs/NRF1_mm0to2_noCGmm.bed.gz | wc -l
# 34285
zcat ./motifs/NRF1_mm0to2_noCGmm_noXY.bed.gz | wc -l
# 32800 (no chrX/Y motifs )

#plot the number of motifs before and after filtering out chrX/Y
Rscript -e '
library(ggplot2)
library(scales)

motifs <- c("BANP", "NRF1")
before_counts <- c(10968, 34285)
after_counts  <- c(10418, 32800)
statuses <- c("Before chrX/Y filtering", "After chrX/Y filtering")

df <- data.frame(
  Motif  = rep(motifs, times = 2),
  Status = rep(statuses, each = length(motifs)),
  Count  = c(before_counts, after_counts)
)

df$Status <- factor(
  df$Status,
  levels = c("Before chrX/Y filtering", "After chrX/Y filtering")
)

p <- ggplot(df, aes(x = Motif, y = Count, fill = Status)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(
    aes(label = comma(Count)),
    position = position_dodge(width = 0.9),
    vjust = -0.3,
    size = 3
  ) +
  theme_minimal() +
  labs(
    title = "Number of motifs before and after chrX/Y filtering",
    subtitle = "Filtering out motifs on chrX and chrY",
    y = "Count",
    x = "Motif"
  ) +
  theme(
    plot.title = element_text(hjust = 0.7, size = 12),
    plot.subtitle = element_text(hjust = 0.7, size = 9)
  )

ggsave(
  "./results/motifs/motif_counts_before_after_chrX_Y_filtering.pdf",
  plot = p,
  width = 6,
  height = 4,
  dpi = 300,
  bg = "white"
)
'


###############################################################################
# 2) MOTIF DISTANCE TO TSS
###############################################################################

# Distance of 2mm motifs to TSS : generate a file per sample with the distances of each motif to the nearest TSS
mkdir -p ./motifs/dist_to_tss/

for f in ./motifs/*_mm0to2_noCGmm.bed.gz ./motifs/*_mm0to2.bed.gz; do
  [ -e "$f" ] || continue

  base=$(basename "$f")
  sample=${base%.bed.gz}

  echo "Processing $sample ..."

  closestBed -t first -d -a "$f" -b /data/genome/annotations/hg38_tss.bed \
    | awk '{print $NF}' > "./motifs/dist_to_tss/${sample}_dist_tss.txt"
done

###############################################################################
# HISTOGRAMS OF DISTANCE TO TSS (ONE PAGE PER TF)
###############################################################################
# Plot histogram of distances to TSS for all samples combined
Rscript -e '
dist_files <- list.files("./motifs/dist_to_tss", full.names = TRUE)
dist_files <- dist_files[endsWith(dist_files, "_dist_tss.txt")]

if(length(dist_files) == 0) stop("No *_dist_tss.txt files found in ./motifs/dist_to_tss/")

pdf("./results/motifs/dist_tss_motifs_hist.pdf")
par(bg = "white")

for (dist_file in dist_files) {

  sample <- sub("_dist_tss.txt", "", basename(dist_file), fixed = TRUE)
  distances <- read.table(dist_file)[,1]

  proximal <- sum(distances < 2000)
  distal   <- sum(distances >= 2000)
  total    <- length(distances)

  pct_prox <- round((proximal / total) * 100, 1)
  pct_dist <- round((distal   / total) * 100, 1)

  hist(log10(distances + 1),
       xlim = c(0, 8),
       xlab = "Distance to TSS (log10)",
       main = paste0(sample, "  |  Proximal = ", pct_prox, "%   |   Distal = ", pct_dist, "%"),
       breaks = 50,
       col = "steelblue",
       border = "black")

  abline(v = log10(2000), col = "red", lwd = 2)
}

dev.off()
cat("[DONE] Wrote: ./results/motifs/dist_tss_motifs_hist.pdf\n")
'

###############################################################################
# 3) MOTIF ∩ HM450 METHYLATION PROBES
###############################################################################

mkdir -p ./motifs/overlaps/intersected_motifs2mm_HM450

# intersect 2mm motifs with annotated methylation data
bedtools intersect -u -a <(zcat ./motifs/NRF1_mm0to2.bed.gz) -b ./methylation/annotated_methylation_data_probes_filtered.bed > ./motifs/overlaps/intersected_motifs2mm_HM450/NRF1_intersected_methylation.bed
bedtools intersect -u -a <(zcat ./motifs/BANP_mm0to2.bed.gz) -b ./methylation/annotated_methylation_data_probes_filtered.bed > ./motifs/overlaps/intersected_motifs2mm_HM450/BANP_intersected_methylation.bed
bedtools intersect -u -a <(zcat ./motifs/NRF1_mm0to2_noCGmm.bed.gz) -b ./methylation/annotated_methylation_data_probes_filtered.bed > ./motifs/overlaps/intersected_motifs2mm_HM450/NRF1_noCGmm_intersected_methylation.bed
bedtools intersect -u -a <(zcat ./motifs/BANP_mm0to2_noCGmm.bed.gz) -b ./methylation/annotated_methylation_data_probes_filtered.bed > ./motifs/overlaps/intersected_motifs2mm_HM450/BANP_noCGmm_intersected_methylation.bed

# Summarize the intersected results :  NRF1mm0to2 7099 BANPmm0to2 2531 NRF1mm0to2_noCGmm 3245 BANPmm0to2_noCGmm 919
mkdir -p ./results/motifs

# Count lines in bash and save to variables
n_NRF1_2mm=$(wc -l < ./motifs/overlaps/intersected_motifs2mm_HM450/NRF1_intersected_methylation.bed)
n_BANP_2mm=$(wc -l < ./motifs/overlaps/intersected_motifs2mm_HM450/BANP_intersected_methylation.bed)
n_NRF1_2mm_noCGmm=$(wc -l < ./motifs/overlaps/intersected_motifs2mm_HM450/NRF1_noCGmm_intersected_methylation.bed)
n_BANP_2mm_noCGmm=$(wc -l < ./motifs/overlaps/intersected_motifs2mm_HM450/BANP_noCGmm_intersected_methylation.bed)

Rscript -e "
tab <- data.frame(
  Motif = c('NRF1','BANP'),
  mm0to2 = c($n_NRF1_2mm, $n_BANP_2mm),
  mm0to2_noCGmm = c($n_NRF1_2mm_noCGmm, $n_BANP_2mm_noCGmm)
)
print(tab)
"

###############################################################################
# VEN DIAGRAM OF MOTIF ∩ HM450 METHYLATION PROBES
###############################################################################
# Create a Venn diagram to visualize the overlaps between 2mm motifs and methylation for NRF1 and BANP : output pdf file with 4 pages : 1 page NRF1 mm0to2, 1 page NRF1 mm0to_noCGmm,  1 page BANP mm0to2, 1 page BANP mm0to_noCGmm.
Rscript -e '
library(VennDiagram)
library(grid)

probes <- as.numeric(
  system("wc -l < ./methylation/annotated_methylation_data_probes_filtered.bed", intern = TRUE)
)

variants <- c("mm0to2", "mm0to2_noCGmm")

# Output
out_pdf <- "./results/motifs/NRF1_BANP_motifs_2mm_vs_methylation_venn.pdf"
pdf(out_pdf, width = 13, height = 13)

for (motif in c("NRF1","BANP")) {
  for (var in variants) {

    # pick motif file + overlap file
    motif_file <- if (var == "mm0to2") {
      paste0("./motifs/", motif, "_mm0to2_noXY.bed.gz")
    } else {
      paste0("./motifs/", motif, "_mm0to2_noCGmm_noXY.bed.gz")
    }

    overlap_file <- if (var == "mm0to2") {
      paste0("./motifs/overlaps/intersected_motifs2mm_HM450/", motif, "_intersected_methylation.bed")
    } else {
      paste0("./motifs/overlaps/intersected_motifs2mm_HM450/", motif, "_noCGmm_intersected_methylation.bed")
    }

    motifs_n <- as.numeric(system(paste0("zcat ", motif_file, " | wc -l"), intern = TRUE))
    overlap_n <- as.numeric(system(paste0("wc -l < ", overlap_file), intern = TRUE))

    grid.newpage()

    venn_obj <- draw.pairwise.venn(
      area1      = motifs_n,
      area2      = probes,
      cross.area = overlap_n,
      category   = c("", ""),
      fill       = c("salmon", "skyblue"),
      alpha      = 0.7,
      cex        = 2,
      cat.cex    = 2.2,
      fontfamily = "Helvetica",
      print.mode = c("raw","percent"),
      sigdig     = 2,
      ind        = FALSE
    )

    grid.text(
      paste0("Overlap between ", motif, " ", var, " motifs and methylation probes"),
      x = 0.5, y = 0.96,
      gp = gpar(fontsize = 30, fontfamily = "Helvetica")
    )

    # keep your manual labels for BANP-like layout (works fine for all pages)
    grid.text("Methylation probes", x = 0.37, y = 0.73,
              gp = gpar(fontsize = 28, fontfamily = "Helvetica"))
    grid.text(paste0(motif, " ", var, " motifs"), x = 0.75, y = 0.65,
              gp = gpar(fontsize = 28, fontfamily = "Helvetica"))

    pushViewport(viewport(y = 0.45))
    grid.draw(venn_obj)
    popViewport()
  }
}

dev.off()
cat("Wrote:", out_pdf, "\n")
'

###############################################################################
# 4) MOTIF ∩ ATAC PEAKS --> merged peaks from all samples has 534,646 peaks
###############################################################################
# Overlap the 2mm motifs with the peaks to see which motifs fall within the peaks regions.
mkdir -p ./motifs/overlaps/motif2mm_peak_overlaps

# Loop through each 2mm motif file and intersect with all peak files from peaks

for file in ./motifs/*_mm0to2_noCGmm_noXY.bed.gz ./motifs/*_mm0to2_noXY.bed.gz; do
  [ -e "$file" ] || continue

  motif_name=$(basename "$file" .bed.gz)

  bedtools intersect -u \
    -a <(zcat "$file") \
    -b ./peaks/filtered_peaks/merged_peaks.bed \
    > ./motifs/overlaps/motif2mm_peak_overlaps/"${motif_name}_peak_overlaps.bed"
done

# Summarize the intersected results : NRF1mm0to2 99073 BANPmm0to2 64685 NRF1mm0to2_noCGmm 23962 BANPmm0to2_noCGmm 5845
# Count lines in bash and save to variables
n_NRF1_2mm=$(wc -l < ./motifs/overlaps/motif2mm_peak_overlaps/NRF1_mm0to2_peak_overlaps.bed)
n_BANP_2mm=$(wc -l < ./motifs/overlaps/motif2mm_peak_overlaps/BANP_mm0to2_peak_overlaps.bed)
n_NRF1_2mm_noCGmm=$(wc -l < ./motifs/overlaps/motif2mm_peak_overlaps/NRF1_mm0to2_noCGmm_peak_overlaps.bed)
n_BANP_2mm_noCGmm=$(wc -l < ./motifs/overlaps/motif2mm_peak_overlaps/BANP_mm0to2_noCGmm_peak_overlaps.bed)

Rscript -e "
tab <- data.frame(
  Motif = c('NRF1','BANP'),
  mm0to2 = c($n_NRF1_2mm, $n_BANP_2mm),
  mm0to2_noCGmm = c($n_NRF1_2mm_noCGmm, $n_BANP_2mm_noCGmm)
)
print(tab)
"

# plot the euler diagram of the 2mm motifs and peaks overlaps for NRF1 and BANP : output pdf file with 4 pages : 1 page NRF1 mm0to2, 1 page NRF1 mm0to_noCGmm,  1 page BANP mm0to2, 1 page BANP mm0to_noCGmm.
Rscript -e '
library(data.table)
library(eulerr)

id3 <- function(f) {
  if (endsWith(f, ".gz")) {
    unique(fread(cmd=paste("zcat", f), header=FALSE)[, paste(V1,V2,V3,sep=":")])
  } else {
    unique(fread(f, header=FALSE)[, paste(V1,V2,V3,sep=":")])
  }
}

dir.create("./results/multi_omics", recursive=TRUE, showWarnings=FALSE)
pdf("./results/multi_omics/euler_2mmmotif_peak_overlaps.pdf", width=7, height=7)

# NRF1
m <- length(id3("./motifs/NRF1_mm0to2_noCGmm_noXY.bed.gz"))
p <- length(id3("./peaks/filtered_peaks/merged_peaks.bed"))
o <- length(id3("./motifs/overlaps/motif2mm_peak_overlaps/NRF1_mm0to2_noCGmm_peak_overlaps.bed"))
plot(euler(c(Motifs=m-o, Peaks=p-o, "Motifs&Peaks"=o)), quantities=TRUE, main="NRF1 mm0to2 noCGmm")

# BANP
m <- length(id3("./motifs/BANP_mm0to2_noCGmm_noXY.bed.gz"))
p <- length(id3("./peaks/filtered_peaks/merged_peaks.bed"))
o <- length(id3("./motifs/overlaps/motif2mm_peak_overlaps/BANP_mm0to2_noCGmm_peak_overlaps.bed"))
plot(euler(c(Motifs=m-o, Peaks=p-o, "Motifs&Peaks"=o)), quantities=TRUE, main="BANP mm0to2 noCGmm")

dev.off()
'


###############################################################################
# 5) MOTIF ∩ PEAK ∩ HM450
###############################################################################
# Overlap methylation-intersected 2mm motifs with peak-overlapping 2mm motifs 
mkdir -p ./motifs/overlaps/intersected_overlaps2mm_peaks_HM450

# ---- mm0to2 ----
bedtools intersect -u \
  -a ./motifs/overlaps/intersected_motifs2mm_HM450/NRF1_intersected_methylation.bed \
  -b ./motifs/overlaps/motif2mm_peak_overlaps/NRF1_mm0to2_peak_overlaps.bed \
  > ./motifs/overlaps/intersected_overlaps2mm_peaks_HM450/NRF1_mm0to2_methylation_peak_overlap.bed

bedtools intersect -u \
  -a ./motifs/overlaps/intersected_motifs2mm_HM450/BANP_intersected_methylation.bed \
  -b ./motifs/overlaps/motif2mm_peak_overlaps/BANP_mm0to2_peak_overlaps.bed \
  > ./motifs/overlaps/intersected_overlaps2mm_peaks_HM450/BANP_mm0to2_methylation_peak_overlap.bed

# ---- mm0to2_noCGmm ----
bedtools intersect -u \
  -a ./motifs/overlaps/intersected_motifs2mm_HM450/NRF1_noCGmm_intersected_methylation.bed \
  -b ./motifs/overlaps/motif2mm_peak_overlaps/NRF1_mm0to2_noCGmm_peak_overlaps.bed \
  > ./motifs/overlaps/intersected_overlaps2mm_peaks_HM450/NRF1_mm0to2_noCGmm_methylation_peak_overlap.bed

bedtools intersect -u \
  -a ./motifs/overlaps/intersected_motifs2mm_HM450/BANP_noCGmm_intersected_methylation.bed \
  -b ./motifs/overlaps/motif2mm_peak_overlaps/BANP_mm0to2_noCGmm_peak_overlaps.bed \
  > ./motifs/overlaps/intersected_overlaps2mm_peaks_HM450/BANP_mm0to2_noCGmm_methylation_peak_overlap.bed

# ---- counts ----
n_NRF1_2mm=$(wc -l < ./motifs/overlaps/intersected_overlaps2mm_peaks_HM450/NRF1_mm0to2_methylation_peak_overlap.bed)
n_BANP_2mm=$(wc -l < ./motifs/overlaps/intersected_overlaps2mm_peaks_HM450/BANP_mm0to2_methylation_peak_overlap.bed)
n_NRF1_2mm_noCGmm=$(wc -l < ./motifs/overlaps/intersected_overlaps2mm_peaks_HM450/NRF1_mm0to2_noCGmm_methylation_peak_overlap.bed)
n_BANP_2mm_noCGmm=$(wc -l < ./motifs/overlaps/intersected_overlaps2mm_peaks_HM450/BANP_mm0to2_noCGmm_methylation_peak_overlap.bed)

# Summarize the intersected results : NRF1mm0to2 5556 BANPmm0to2 1781 NRF1mm0to2_noCGmm 2837 BANPmm0to2_noCGmm 752
# ---- summary table ----
Rscript -e "
tab <- data.frame(
  Motif = c('NRF1','BANP'),
  mm0to2 = c($n_NRF1_2mm, $n_BANP_2mm),
  mm0to2_noCGmm = c($n_NRF1_2mm_noCGmm, $n_BANP_2mm_noCGmm)
)
print(tab)
"

###############################################################################
# 6) BARPLOTS – MOTIFS ACROSS DATASETS
###############################################################################

###############################################################################
# ONE BARPLOT PER TF : MOTIFS, ATAC, HM450, SNVs 
###############################################################################
# create a barplot : One page per TF: 2mm motifs, ATAC, HM450, SNVs
Rscript -e '
library(ggplot2)
library(scales)

out_pdf <- "./results/motifs/NRF1_BANP_mm0to2_variants_vs_peaks_vs_HM450_SNV_barplots.pdf"
dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)

pdf(out_pdf, width = 10, height = 8)

plot_one <- function(motif, variant){

  motif_all_file <- if(variant == "mm0to2"){
    paste0("./motifs/", motif, "_mm0to2_noXY.bed.gz")
  } else {
    paste0("./motifs/", motif, "_mm0to2_noCGmm_noXY.bed.gz")
  }

  motif_peak_file <- if(variant == "mm0to2"){
    paste0("./motifs/overlaps/motif2mm_peak_overlaps/", motif, "_mm0to2_peak_overlaps.bed")
  } else {
    paste0("./motifs/overlaps/motif2mm_peak_overlaps/", motif, "_mm0to2_noCGmm_peak_overlaps.bed")
  }

  motif_hm450_file <- if(variant == "mm0to2"){
    paste0("./motifs/overlaps/intersected_motifs2mm_HM450/", motif, "_intersected_methylation.bed")
  } else {
    paste0("./motifs/overlaps/intersected_motifs2mm_HM450/", motif, "_noCGmm_intersected_methylation.bed")
  }

  motif_snv_file <- if(variant == "mm0to2"){
    paste0("./motifs/overlaps/motif2mm_snv_overlaps/", motif, "_mm0to2_motifs_with_SNVs.bed")
  } else {
    paste0("./motifs/overlaps/motif2mm_snv_overlaps/", motif, "_mm0to2_noCGmm_motifs_with_SNVs.bed")
  }

  motif_peak_snv_file <- if(variant == "mm0to2"){
    paste0("./motifs/overlaps/motifs2mm_in_peaks_snv/", motif, "_mm0to2_motifs_in_peaks_with_SNVs.bed")
  } else {
    paste0("./motifs/overlaps/motifs2mm_in_peaks_snv/", motif, "_mm0to2_noCGmm_motifs_in_peaks_with_SNVs.bed")
  }

  motif_all <- read.table(gzfile(motif_all_file))
  motif_ids <- unique(with(motif_all, paste(V1, V2, V3, sep=":")))

  motif_peak <- read.table(motif_peak_file)
  peak_ids   <- unique(with(motif_peak, paste(V1, V2, V3, sep=":")))

  motif_hm450 <- read.table(motif_hm450_file)
  hm450_ids   <- unique(with(motif_hm450, paste(V1, V2, V3, sep=":")))

  motif_snv <- read.table(motif_snv_file)
  snv_motif_ids <- unique(with(motif_snv, paste(V1, V2, V3, sep=":")))

  motif_peak_snv <- read.table(motif_peak_snv_file)
  snv_peak_motif_ids <- unique(with(motif_peak_snv, paste(V1, V2, V3, sep=":")))

  total        <- length(motif_ids)
  in_atac      <- length(peak_ids)
  in_hm450     <- length(hm450_ids)
  in_both      <- length(intersect(peak_ids, hm450_ids))
  in_snv       <- length(snv_motif_ids)
  in_atac_snv  <- length(snv_peak_motif_ids)

  df <- data.frame(
    category = c(paste0("All ", motif, " motifs (", variant, ")"),
                 "Motifs in ATAC peaks",
                 "Motifs with HM450 probes",
                 "Motifs in both ATAC + HM450",
                 "Motifs containing SNVs",
                 "Motifs in ATAC peaks containing SNVs"),
    count = c(total, in_atac, in_hm450, in_both, in_snv, in_atac_snv)
  )

  df$category <- factor(df$category, levels = df$category)
  df$percent  <- df$count / total * 100
  df$label    <- sprintf("%s (%.2f%%)", comma(df$count), df$percent)

  p <- ggplot(df, aes(x = reorder(category, -count), y = count, fill = category)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = label), vjust = -0.4, size = 4) +
    scale_y_continuous(labels = comma) +   # commas on y-axis
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 20, hjust = 1),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    ) +
    labs(
      title = paste0(motif, " motif distribution across datasets"),
      subtitle = sprintf("All %s motifs (%s) (n = %s)", motif, variant, comma(total)),
      y = "Number of motifs",
      x = NULL
    )

  print(p)
}

plot_one("NRF1", "mm0to2")
plot_one("NRF1", "mm0to2_noCGmm")
plot_one("BANP", "mm0to2")
plot_one("BANP", "mm0to2_noCGmm")

dev.off()
cat("Wrote:", out_pdf, "\n")
'


###############################################################################
# ONE BARPLOT FOR BOTH TF : MOTIFS, ATAC, HM450, SNVs NOT DONE
###############################################################################
# one barplot for both TFs NRF1 and BANP together : this script creates a small barplot showing the percentage of motifs in each category for both TFs one tf per page
Rscript -e '
library(ggplot2)

# ===============================
# TF-specific file definitions
# ===============================
tf_files <- list(
  BANP = list(
    motif = "./motifs/BANP_mm0to2_noCGmm_noXY.bed.gz",
    hm450 = "./motifs/overlaps/intersected_motifs2mm_HM450/BANP_noCGmm_intersected_methylation.bed",
    snv   = "./motifs/overlaps/motif2mm_snv_overlaps//BANP_mm0to2_noCGmm_motifs_with_SNVs.bed"
  ),
  NRF1 = list(
    motif = "./motifs/NRF1_mm0to2_noCGmm_noXY.bed.gz",
    hm450 = "./motifs/overlaps/intersected_motifs2mm_HM450/NRF1_noCGmm_intersected_methylation.bed",
    snv   = "./motifs/overlaps/motif2mm_snv_overlaps//NRF1_mm0to2_noCGmm_motifs_with_SNVs.bed"
  )
)

out_pdf <- "./results/motifs/BANP_NRF1_motif_probes_snv_barplots_oneaxis.pdf"
dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)

all_df <- list()

for (tf in names(tf_files)) {

  ## 1) All motifs
  motif_all <- read.table(tf_files[[tf]]$motif)
  motif_ids <- unique(with(motif_all, paste(V1, V2, V3, sep = ":")))

  ## 2) HM450 overlaps
  motif_hm450 <- read.table(tf_files[[tf]]$hm450)
  hm450_ids <- unique(with(motif_hm450, paste(V1, V2, V3, sep = ":")))

  ## 3) SNV overlaps
  motif_snv <- read.table(tf_files[[tf]]$snv)
  snv_ids <- unique(with(motif_snv, paste(V1, V2, V3, sep = ":")))

  ## 4) Counts
  total    <- length(motif_ids)
  in_snv   <- length(snv_ids)
  in_hm450 <- length(hm450_ids)

  ## 5) Per-TF table
  df <- data.frame(
    tf = tf,
    category = c("All motifs", "Motifs with variants", "Motifs with HM450 probes"),
    count = c(total, in_snv, in_hm450),
    stringsAsFactors = FALSE
  )

  df$percent <- df$count / total * 100
  df$label <- sprintf(
    "%s (%.2f%%)",
    format(df$count, big.mark = ",", scientific = FALSE, trim = TRUE),
    df$percent
  )

  all_df[[tf]] <- df
}

plot_df <- do.call(rbind, all_df)
rownames(plot_df) <- NULL

plot_df$tf <- factor(plot_df$tf, levels = c("BANP", "NRF1"))
plot_df$category <- factor(
  plot_df$category,
  levels = c("All motifs", "Motifs with variants", "Motifs with HM450 probes")
)

pdf(out_pdf, width = 11, height = 7)

p <- ggplot(plot_df, aes(x = tf, y = count, fill = category)) +
  geom_col(position = position_dodge(width = 0.6), width = 0.45) +
  geom_text(
    aes(label = label),
    position = position_dodge(width = 0.6),
    vjust = -0.25,
    size = 3
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  ) +
  labs(
    title = "BANP and NRF1 motif distribution across datasets",
    subtitle = "Counts and percentages relative to all motifs for each TF",
    x = "Transcription factor",
    y = "Number of motifs",
    fill = NULL
  ) +
  scale_y_continuous(
    labels = function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
  )

print(p)
dev.off()

cat("DONE ->", out_pdf, "\n")
'

###############################################################################
# PROPORTIONAL EULER DIAGRAM PER TF : 2mm MOTIFS, ATAC, HM450, SNVs 
###############################################################################
Rscript -e '
library(eulerr)

out_pdf <- "./results/motifs/Euler_VENN_NRF1_BANP_mm0to2_Motifs_ATAC_HM450_SNV.pdf"
dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)

pdf(out_pdf, width = 14, height = 14)
par(mar = c(4, 4, 6, 2))

motifs   <- c("NRF1","BANP")
variants <- c("mm0to2","mm0to2_noCGmm")

read_ids <- function(path, gz=FALSE){
  if(!file.exists(path)){
    warning("Missing file: ", path)
    return(character(0))
  }
  x <- if(gz) read.table(gzfile(path)) else read.table(path)
  unique(with(x, paste(V1, V2, V3, sep=":")))
}

for(motif in motifs){
  for(variant in variants){

    # ---------- inputs for this motif/variant ----------
    motif_all_file  <- paste0("./motifs/", motif, "_", variant, "noCGmm_noXY.bed.gz")
    motif_peak_file <- paste0("./motifs/overlaps/motif2mm_peak_overlaps/", motif, "_", variant, "_peak_overlaps.bed")
    motif_snv_file  <- paste0("./motifs/overlaps/motif2mm_snv_overlaps/", motif, "_", variant, "_motifs_with_SNVs.bed")

    motif_hm450_file <- if(variant == "mm0to2"){
      paste0("./motifs/overlaps/intersected_motifs2mm_HM450/", motif, "_intersected_methylation.bed")
    } else {
      paste0("./motifs/overlaps/intersected_motifs2mm_HM450/", motif, "_noCGmm_intersected_methylation.bed")
    }

    motif_all      <- read_ids(motif_all_file, gz=TRUE)
    motif_peak     <- read_ids(motif_peak_file)
    motif_variants <- read_ids(motif_snv_file)
    motif_hm450    <- read_ids(motif_hm450_file)

    all_name <- paste0("All ", motif, " motifs (", variant, ")")

    ############################################
    ## PAGE A — (same as PAGE 1/3): ATAC + SNVs ##
    ############################################
    plot.new()

    fit <- euler(c(
      setNames(list(motif_all), all_name),
      list(
        "Motifs in ATAC"         = motif_peak,
        "Motifs with variants"   = motif_variants
      )
    ))

    p <- plot(fit,
              fills = c("steelblue","lightgreen","red"),
              edges = "black",
              quantities = list(type="counts", cex=0.9),
              labels = list(cex=1.1),
              main = NULL)
    print(p)

    mtext(paste0(motif, " motifs (", variant, "): genome-wide, accessible,\nand overlapping SNVs"),
          side = 3, line = 2, cex = 1.2)

    #############################################
    ## PAGE B — (same as PAGE 2/4): ATAC + HM450 ##
    #############################################
    plot.new()

    fit <- euler(c(
      setNames(list(motif_all), all_name),
      list(
        "Motifs in ATAC"            = motif_peak,
        "Motifs with HM450 probes"  = motif_hm450
      )
    ))

    p <- plot(fit,
              fills = c("steelblue","lightgreen","orange"),
              edges = "black",
              quantities = list(type="counts", cex=0.9),
              labels = list(cex=1.1),
              main = NULL)
    print(p)

    text(x = 0.8, y = 0.95,
        labels = "NRF1 motifs: genome-wide, accessible,\nand overlapping SNVs",
        cex = 1.2)
  }
}

dev.off()
cat("DONE ->", out_pdf, "\n")
'
###############################################################################
# 7) Linking motifs with the closest gene using bedtools closest 
###############################################################################
# Gene annotation file : /data/genome/annotations/hg38_genes.bed
for files in ./motifs/*_mm0to2_noCGmm_noXY.bed.gz; do 
    echo "Processing $files ..."
    zcat "$files" | sort -k1,1 -k2,2n | closestBed -t "first" -d -a stdin -b /data/genome/annotations/hg38_tss.bed > "${files%.bed.gz}_closest_genes.bed"
done

# add header to the closest genes files
for file in ./motifs/*_mm0to2_noCGmm_noXY_closest_genes.bed; do
  tmp="${file}.tmp"
  { echo -e "chrom_motif\tstart_motif\tend_motif\tmotif\tpvalue\tstrand_motif\tchrom_gene\tstart_gene\tend_gene\tgene\tstrand_gene\tgene_id\tdistance"; cat "$file"; } > "$tmp" && mv "$tmp" "$file"
done


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

########################################################################################
# 2) FILTERING OUT PEAKS FILES --> 217 FILES LEFT and filter out peaks on chrX and chrY
########################################################################################
# Filter out the peak files corresponding to failing samples from the previous analysis of /data/hichamif/pred_tf_cancer/QC_results/failing_samples_by_step/all_failing_samples.txt
# Filters : NFR > 15% Mapping rate > 80% Peakcount >20,000 FRIP > 0.2
mkdir -p ./peaks/filtered_peaks/

FAIL_TSV="./all_failing_samples.txt"

for file in ./peaks/ATAC_TCGA-*_peaks_macs.bed; do
    base=$(basename "$file")
    sample=${base%_peaks_macs.bed}
    out="./peaks/filtered_peaks/$base"

    echo "Processing $sample ..."

    if awk -F'\t' -v s="$sample" 'NR>1 && $1==s {found=1} END{exit !found}' "$FAIL_TSV"; then
        echo "Skipping $sample (in failing list)."
        continue
    fi

    awk 'BEGIN{FS=OFS="\t"} $1 != "chrX" && $1 != "chrY"' "$file" > "$out"
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

# plot the barplot of number of unique peaks per cancer type
Rscript -e '
# 1. Read the peak counts
data <- read.table("./peaks/peaks_counts_per_cancer_type/peaks_unique_per_cancer.tsv", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

cancer_names <- data[,1]
peak_counts  <- data[,2]

# Sort and reverse for horizontal barplot (largest at top)
ord <- order(peak_counts, decreasing = TRUE)
cancer_names <- rev(cancer_names[ord])
peak_counts  <- rev(peak_counts[ord])
peak_thousands <- peak_counts / 1e3

# 2. Read the color file - CRITICAL: added comment.char = ""
color_df <- read.table("./results/multi_omics/cancer_color_order_with_defined_colours.tsv", 
                       header = FALSE, 
                       sep = "\t", 
                       stringsAsFactors = FALSE, 
                       comment.char = "") # This prevents # hex codes from being ignored

colnames(color_df) <- c("Cancer_full", "Color")

# 3. Clean names and Build Palette
# Strip TCGA-, trim whitespace, and force uppercase for a perfect match
color_df$Clean_Name <- toupper(trimws(sub("^TCGA-", "", color_df$Cancer_full)))
palette <- color_df$Color
names(palette) <- color_df$Clean_Name

# 4. Map colors to the peaks data
clean_peak_names <- toupper(trimws(cancer_names))
bar_colors <- palette[clean_peak_names]

# Fallback for any unmatched names
bar_colors[is.na(bar_colors)] <- "gray"

# 5. Plotting
pdf("./results/peaks/peaks_per_cancer_type_barplot.pdf", width = 12, height = 9, bg = "white")

# Give the left margin (2nd value) more room for cancer names
par(mar = c(5, 10, 4, 2)) 

barplot(
  peak_thousands,
  names.arg = cancer_names,
  horiz = TRUE,
  las = 1,
  col = bar_colors,
  border = "black",
  main = "Peaks per Cancer Type",
  cex.main = 2,
  xlab = "Total Peaks (thousands)"
)

dev.off()

cat("Success! Plot generated with hex colors.\n")
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
out_pdf <- "./results/peaks/ATAC_peak_counts_per_sample_by_cancer.pdf"

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
# NBR OF samples with peaks per cancer type
###############################################################################
mkdir -p ./peaks/peak_counts_per_sample/

{
  echo -e "cancer\tsample_count"
  for peakfile in ./peaks/filtered_peaks/ATAC_TCGA-*-*.bed; do
    sample=$(basename "$peakfile" .bed)
    cancer=$(echo "$sample" | cut -d'-' -f2 | cut -d'_' -f1)
    echo "$cancer"
  done | sort | uniq -c | awk '{print "TCGA-"$2"\t"$1}'
} > ./peaks/peak_counts_per_sample/samples_per_cancer.tsv

# plot the barplot of number of samples with peaks per cancer type
Rscript -e '
data <- read.table("./peaks/peak_counts_per_sample/samples_per_cancer.tsv",header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)
data$sample_count <- as.numeric(data$sample_count)
data <- data[order(data$sample_count, decreasing=FALSE), ]
cols <- read.table("./results/multi_omics/cancer_color_order_with_defined_colours.tsv",header=FALSE, sep="\t", stringsAsFactors=FALSE,quote="", comment.char="")
colnames(cols) <- c("cancer","color")
palette <- setNames(cols$color, cols$cancer)
bar_colors <- unname(palette[data$cancer])

# validate colors
bad <- sapply(bar_colors, function(x) {
  inherits(try(col2rgb(x), silent=TRUE), "try-error")
})
if (any(bad)) {
  cat("[ERROR] These colors are invalid:\n")
  print(unique(bar_colors[bad]))
  quit(status=1)
}

pdf("./results/peaks/samples_with_peaks_per_cancer_barplot.pdf", width=9, height=9, bg="white")

nice_max <- max(pretty(c(0, data$sample_count)))
labels <- sub("^TCGA-", "", data$cancer)

par(mar=c(5,8,4,2))
barplot(
  data$sample_count,
  names.arg=labels,
  horiz=TRUE,
  las=1,
  col=bar_colors,
  border="black",
  main="Samples with peaks per cancer type",
  cex.main=2,
  cex.lab=1.1,
  cex.names=0.8,
  xlab="Number of samples",
  xlim=c(0, nice_max)
)

dev.off()
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
pdf("./results/peaks/dist_tss_peaks_hist.pdf", width = 11.7, height = 8.3);
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
pdf("./results/peaks/boxplot_ATAC_peaks_proximal_distal.pdf", width = 7, height = 5)
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
#########################################################
# Look into the ATAC peaks overlapping with 2mm motifs
#########################################################
mkdir -p ./peaks/overlaps/intersected_motifs/NRF1_mm0to2

motif_sites="./motifs/NRF1_mm0to2_noXY.bed.gz"

for f in ./peaks/filtered_peaks/ATAC_TCGA*peaks_macs.bed; do
  sample=$(basename "$f" _peaks_macs.bed)
  out="./peaks/overlaps/intersected_motifs/NRF1_mm0to2/${sample}_NRF1_intersected_peaks.bed"

  bedtools intersect -u -a "$f" -b <(zcat "$motif_sites") > "$out"
  echo "Wrote: $out"
done

mkdir -p ./peaks/overlaps/intersected_motifs/BANP_mm0to2

motif_sites="./motifs/BANP_mm0to2_noXY.bed.gz"

for f in ./peaks/filtered_peaks/ATAC_TCGA*peaks_macs.bed; do
  sample=$(basename "$f" _peaks_macs.bed)
  out="./peaks/overlaps/intersected_motifs/BANP_mm0to2/${sample}_BANP_intersected_peaks.bed"

  bedtools intersect -u -a "$f" -b <(zcat "$motif_sites") > "$out"
  echo "Wrote: $out"
done

mkdir -p ./peaks/overlaps/intersected_motifs/NRF1_mm0to2_noCGmm

motif_sites="./motifs/NRF1_mm0to2_noCGmm_noXY.bed.gz"

for f in ./peaks/filtered_peaks/ATAC_TCGA*peaks_macs.bed; do
  sample=$(basename "$f" _peaks_macs.bed)
  out="./peaks/overlaps/intersected_motifs/NRF1_mm0to2_noCGmm/${sample}_NRF1_intersected_peaks.bed"

  bedtools intersect -u -a "$f" -b <(zcat "$motif_sites") > "$out"
  echo "Wrote: $out"
done

mkdir -p ./peaks/overlaps/intersected_motifs/BANP_mm0to2_noCGmm

motif_sites="./motifs/BANP_mm0to2_noCGmm_noXY.bed.gz"

for f in ./peaks/filtered_peaks/ATAC_TCGA*peaks_macs.bed; do
  sample=$(basename "$f" _peaks_macs.bed)
  out="./peaks/overlaps/intersected_motifs/BANP_mm0to2_noCGmm/${sample}_BANP_intersected_peaks.bed"

  bedtools intersect -u -a "$f" -b <(zcat "$motif_sites") > "$out"
  echo "Wrote: $out"
done

###############################################################################
# 7) 2mm MOTIF THAT OVERLAP WITH PEAKS DISTANCE TO TSS
###############################################################################
# Distance of motifs in peaks to TSS : generate a file per sample with the distances of each motif to the nearest TSS

# NRF1_mm0to2
mkdir -p ./peaks/dist_to_tss_motifs/NRF1_mm0to2

for f in ./peaks/overlaps/intersected_motifs/NRF1_mm0to2/*_NRF1_intersected_peaks.bed; do
  base=$(basename "$f" _NRF1_intersected_peaks.bed)
  out="./peaks/dist_to_tss_motifs/NRF1_mm0to2/${base}_NRF1_dist_tss.txt"

  closestBed -t first -d \
    -a "$f" \
    -b /data/genome/annotations/hg38_tss.bed \
  | awk '{print $NF}' > "$out"

  echo "Wrote: $out"
done

# NRF1_mm0to2_noCGmm
mkdir -p ./peaks/dist_to_tss_motifs/NRF1_mm0to2_noCGmm

for f in ./peaks/overlaps/intersected_motifs/NRF1_mm0to2_noCGmm/*_NRF1_intersected_peaks.bed; do
  base=$(basename "$f" _NRF1_intersected_peaks.bed)
  out="./peaks/dist_to_tss_motifs/NRF1_mm0to2_noCGmm/${base}_NRF1_dist_tss.txt"

  closestBed -t first -d \
    -a "$f" \
    -b /data/genome/annotations/hg38_tss.bed \
  | awk '{print $NF}' > "$out"

  echo "Wrote: $out"
done

# BANP_mm0to2
mkdir -p ./peaks/dist_to_tss_motifs/BANP_mm0to2

for f in ./peaks/overlaps/intersected_motifs/BANP_mm0to2/*_BANP_intersected_peaks.bed; do
  base=$(basename "$f" _BANP_intersected_peaks.bed)
  out="./peaks/dist_to_tss_motifs/BANP_mm0to2/${base}_BANP_dist_tss.txt"

  closestBed -t first -d \
    -a "$f" \
    -b /data/genome/annotations/hg38_tss.bed \
  | awk '{print $NF}' > "$out"

  echo "Wrote: $out"
done

# BANP_mm0to2_noCGmm
mkdir -p ./peaks/dist_to_tss_motifs/BANP_mm0to2_noCGmm

for f in ./peaks/overlaps/intersected_motifs/BANP_mm0to2_noCGmm/*_BANP_intersected_peaks.bed; do
  base=$(basename "$f" _BANP_intersected_peaks.bed)
  out="./peaks/dist_to_tss_motifs/BANP_mm0to2_noCGmm/${base}_BANP_dist_tss.txt"

  closestBed -t first -d \
    -a "$f" \
    -b /data/genome/annotations/hg38_tss.bed \
  | awk '{print $NF}' > "$out"

  echo "Wrote: $out"
done

# Plot histogram of distances to TSS for NRF1 and BANP motifs per cancer type --> takes the median of all the samples per cancer type.
Rscript -e '
library(data.table)

out_pdf <- "./results/peaks/ATAC_motifs_TSS/dist_tss_ATAC_NRF1_BANP_mm0to2_byCancer.pdf"
dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)

pdf(out_pdf, width = 11.7, height = 8.3)
par(bg = "white")

# Folders encode the variant
targets <- c("NRF1_mm0to2","NRF1_mm0to2_noCGmm","BANP_mm0to2","BANP_mm0to2_noCGmm")

get_cancer <- function(x) {
  m <- regexpr("TCGA-[A-Z0-9]+", x)
  if (m[1] == -1) return(NA_character_)
  sub("^TCGA-", "", regmatches(x, m))
}

for (target in targets) {

  TF <- sub("_mm0to2(_noCGmm)?$", "", target)   # "NRF1" or "BANP"
  dist_dir <- file.path("./peaks/dist_to_tss_motifs", target)

  out_tsv <- file.path("./results/peaks/ATAC_motifs_TSS", target,
                       paste0("summary_", target, "_byCancer.tsv"))
  dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)

  # Your dist files end with _NRF1_dist_tss.txt or _BANP_dist_tss.txt
  files <- list.files(dist_dir, pattern = paste0("_", TF, "_dist_tss.txt$"), full.names = TRUE)

  if (!length(files)) {
    plot.new()
    title(main = paste0(target, ": no distance files found in ", dist_dir))
    next
  }

  cancers <- vapply(basename(files), get_cancer, character(1))
  ok <- !is.na(cancers) & cancers != ""
  files <- files[ok]
  cancers <- cancers[ok]

  uniq_cancers <- sort(unique(cancers))
  res <- list()

  for (c in uniq_cancers) {
    fc <- files[cancers == c]

    d_list <- lapply(fc, function(f) suppressWarnings(as.numeric(readLines(f))))
    d <- unlist(d_list, use.names = FALSE)
    d <- d[is.finite(d)]

    n_samples <- length(fc)
    n_vals <- length(d)

    if (n_vals == 0) {
      plot.new()
      title(main = paste0("TCGA-", c, " | ", target, ": none found"))
      res[[c]] <- data.frame(cancer=c, n_samples=n_samples, n_values=0,
                             pct_prox=NA, pct_dist=NA, median=NA)
      next
    }

    proximal <- sum(d < 2000)
    distal   <- sum(d >= 2000)
    pct_prox <- round(100 * proximal / n_vals, 1)
    pct_dist <- round(100 * distal   / n_vals, 1)

    hist(log10(d + 1),
         xlim = c(0, 8),
         xlab = "Distance to TSS (log10)",
         main = paste0("TCGA-", c, " | ", target,
                       " | samples=", n_samples,
                       " | n=", n_vals,
                       " | Proximal=", pct_prox,
                       "% | Distal=", pct_dist, "%"),
         breaks = 50,
         col = "steelblue",
         border = "black")
    abline(v = log10(2000), col = "red", lwd = 2)

    res[[c]] <- data.frame(cancer=c, n_samples=n_samples, n_values=n_vals,
                           pct_prox=pct_prox, pct_dist=pct_dist, median=median(d))
  }

  write.table(do.call(rbind, res), out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
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
mkdir -p ./methylation/

for file in /data/papers/tcga/TCGA*/*/HM450*.txt; do
    ln -s "$file" ./methylation/
done

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

# create a list of all filtered methylation files for healthy samples only 11A_1 tumor type first replicate and samples while excluding cancer samples 
find ./methylation/filtered_methylation/ -type f -name "HM450*_annotated_methylation_filtered.bed.gz" | grep -E -- '-11A_1' > ./methylation/filtered_methylation_files_tumor11A_1.txt

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
pdf("./results/methylation/methylation_counts_barplot.pdf", width=12, height=8)

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
###############################################################################
# HISTOGRAM OF METHYLATION PERCENTAGE FOR ALL SAMPLES
###############################################################################
# Plot histogram of methylation percentage for all samples and highlight UMR, LMR, FMR regions with different colors and add legend with percentages of each region.
Rscript -e '
# One PDF with all samples, one per page — same plot logic

out_pdf <- "./results/methylation/Methylation_distribution_UMR_LMR_FMR_ALL_samples_histogram.pdf"
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

pdf("./results/methylation/methylation_region_percentages_barplot.pdf",width = 14, height = 8)

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
pdf("./results/methylation/methylation_region_percentages_boxplot_per_sample.pdf",width = 16, height = 9)

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
distances <- read.table(dist_file, header = FALSE)[,1]

# Remove missing values if any
distances <- distances[!is.na(distances)]

pdf("./results/methylation/dist_tss_methylation_histogram.pdf", width = 10, height = 10)
par(bg = "white")

# Compute proximal/distal percentages using absolute distance to TSS
proximal <- sum(abs(distances) < 2000)
distal   <- sum(abs(distances) >= 2000)
total    <- length(distances)

pct_prox <- round((proximal / total) * 100, 1)
pct_dist <- round((distal   / total) * 100, 1)

hist(
  log10(abs(distances) + 1),
  xlim = c(0, 8),
  xlab = "Absolute distance to TSS (log10)",
  main = paste0(
    "All common CpGs",
    " | Proximal = ", pct_prox, "%",
    " Distal = ", pct_dist, "%"
  ),
  breaks = 50,
  col = "steelblue",
  border = "black"
)

abline(v = log10(2000 + 1), col = "red", lwd = 2)

dev.off()
'


###############################################################################
# SECTION 8 — DELTA METHYLATION (ALL CpGs after filtering))
###############################################################################
###########################################################################################
# Overlaps Methylation probes with 2mm motifs
########################################################################################### 
mkdir -p ./methylation/overlaps/intersected_motifs2mm_HM450

PROBES=./methylation/annotated_methylation_data_probes_filtered.bed

# NRF1 mm0to2
bedtools intersect -u \
  -a "$PROBES" \
  -b <(zcat ./motifs/NRF1_mm0to2_noXY.bed.gz) \
  > ./methylation/overlaps/intersected_motifs2mm_HM450/NRF1_mm0to2_intersected_methylation.bed

# NRF1 mm0to2 noCGmm
bedtools intersect -u \
  -a "$PROBES" \
  -b <(zcat ./motifs/NRF1_mm0to2_noCGmm_noXY.bed.gz) \
  > ./methylation/overlaps/intersected_motifs2mm_HM450/NRF1_mm0to2_noCGmm_intersected_methylation.bed

# BANP mm0to2
bedtools intersect -u \
  -a "$PROBES" \
  -b <(zcat ./motifs/BANP_mm0to2_noXY.bed.gz) \
  > ./methylation/overlaps/intersected_motifs2mm_HM450/BANP_mm0to2_intersected_methylation.bed

# BANP mm0to2 noCGmm
bedtools intersect -u \
  -a "$PROBES" \
  -b <(zcat ./motifs/BANP_mm0to2_noCGmm_noXY.bed.gz) \
  > ./methylation/overlaps/intersected_motifs2mm_HM450/BANP_mm0to2_noCGmm_intersected_methylation.bed


###########################################################################################
# Overlaps Methylaiton probes with 2mm motifs in peaks 
########################################################################################### 
mkdir -p ./methylation/overlaps/intersected_motifs2mm_HM450_peaks

PROBES=./methylation/annotated_methylation_data_probes_filtered.bed

# NRF1 mm0to2
bedtools intersect -u \
  -a "$PROBES" \
  -b <(cat ./motifs/overlaps/motif2mm_peak_overlaps/NRF1_mm0to2_peak_overlaps.bed) \
  > ./methylation/overlaps/intersected_motifs2mm_HM450_peaks/NRF1_mm0to2_intersected_methylation.bed

# NRF1 mm0to2 noCGmm
bedtools intersect -u \
  -a "$PROBES" \
  -b <(cat ./motifs/overlaps/motif2mm_peak_overlaps/NRF1_mm0to2_noCGmm_peak_overlaps.bed) \
  > ./methylation/overlaps/intersected_motifs2mm_HM450_peaks/NRF1_mm0to2_noCGmm_intersected_methylation.bed

# BANP mm0to2
bedtools intersect -u \
  -a "$PROBES" \
  -b <(cat ./motifs/overlaps/motif2mm_peak_overlaps/BANP_mm0to2_peak_overlaps.bed) \
  > ./methylation/overlaps/intersected_motifs2mm_HM450_peaks/BANP_mm0to2_intersected_methylation.bed

# BANP mm0to2 noCGmm
bedtools intersect -u \
  -a "$PROBES" \
  -b <(cat ./motifs/overlaps/motif2mm_peak_overlaps/BANP_mm0to2_noCGmm_peak_overlaps.bed) \
  > ./methylation/overlaps/intersected_motifs2mm_HM450_peaks/BANP_mm0to2_noCGmm_intersected_methylation.bed  


#############################################################################################################################
# SMOOTH SCATTER PLOT FOR ALL CANCER-HEALTHY SAMPLE PAIRS FILTERED with 2mm motif for one cancer type , healthy vs tumor
#############################################################################################################################
# plot the methylation distribution for both healthy and cancer samples of all samples using a delta methylation histogram :
#This script generates a single PDF where each page contains three plots for one tumor–normal pair: global CpG methylation, CpGs in NRF1 motifs, and CpGs in BANP motifs.
#In the second and third plots, CpGs overlapping NRF1 or BANP binding sites with strong methylation changes (|Δβ| ≥ 20%) are highlighted in red to assess motif-specific epigenetic disruption.

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
  NRF1_mm0to2        = "./methylation/overlaps/intersected_motifs2mm_HM450/NRF1_mm0to2_intersected_methylation.bed",
  BANP_mm0to2        = "./methylation/overlaps/intersected_motifs2mm_HM450/BANP_mm0to2_intersected_methylation.bed",
  NRF1_noCGmm_mm0to2 = "./methylation/overlaps/intersected_motifs2mm_HM450/NRF1_mm0to2_noCGmm_intersected_methylation.bed",
  BANP_noCGmm_mm0to2 = "./methylation/overlaps/intersected_motifs2mm_HM450/BANP_mm0to2_noCGmm_intersected_methylation.bed"
)

out_pdf <- "./results/methylation/smoothscatter/smoothScatter_mm0to2_vs_noCGmm_mm0to2_per_pair.pdf"
thr <- 20
dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)

# =========================
# LOAD DATA
# =========================
pairs <- fread(pairs_file)
stopifnot(all(c("cancer","patient","tumor_file","healthy_file") %in% colnames(pairs)))

motif_cpgs <- lapply(motif_files, function(f){
  if(!file.exists(f)) stop(paste("Missing motif file:", f))
  unique(fread(f, header = FALSE)[[4]])
})

# =========================
# HELPERS
# =========================
pct_fmt <- function(x, n) {
  if (n == 0) return("0%")
  sprintf("%.1f%%", 100 * x / n)
}

add_class_legend <- function(delta_vec, thr, title = "All CpGs") {
  n <- length(delta_vec)
  hypo <- sum(delta_vec <= -thr, na.rm = TRUE)
  hyper <- sum(delta_vec >=  thr, na.rm = TRUE)
  unch <- sum(abs(delta_vec) < thr, na.rm = TRUE)

  lab <- paste0(
    title, "\n",
    "Hypo (Delta≤-", thr, "): ", pct_fmt(hypo, n), "\n",
    "Unchanged (|Delta|<", thr, "): ", pct_fmt(unch, n), "\n",
    "Hyper (Delta≥+", thr, "): ", pct_fmt(hyper, n)
  )

  usr <- par("usr")
  x <- usr[1] + 0.02 * (usr[2] - usr[1])
  y <- usr[4] - 0.02 * (usr[4] - usr[3])
  text(x, y, labels = lab, adj = c(0, 1), cex = 0.85)
}

read_bed_gz <- function(path) {
  cmd <- paste("zcat", shQuote(path))
  fread(cmd = cmd, header = FALSE, showProgress = FALSE)
}

plot_all <- function(merged, delta, thr, main_title){
  smoothScatter(
    merged$meth_normal, merged$meth_tumor,
    xlab = "Healthy methylation (%)",
    ylab = "Tumor methylation (%)",
    main = main_title
  )
  abline(0, 1, col = "black", lwd = 2)
  abline(a = -thr, b = 1, col = "blue", lwd = 2)
  abline(a =  thr, b = 1, col = "blue", lwd = 2)
  add_class_legend(delta, thr, title = "All CpGs")
}

plot_motif_red <- function(merged, delta, thr, in_motif, label_title){
  smoothScatter(
    merged$meth_normal, merged$meth_tumor,
    xlab = "Healthy methylation (%)",
    ylab = "Tumor methylation (%)",
    main = paste0(label_title, " in red (|DeltaBeta| ≥ ", thr, "%)")
  )
  is_red <- in_motif & (abs(delta) >= thr)
  points(
    merged$meth_normal[is_red],
    merged$meth_tumor[is_red],
    pch = 16, cex = 0.5, col = rgb(1, 0, 0, 0.6)
  )
  abline(0, 1, col = "black", lwd = 2)
  abline(a = -thr, b = 1, col = "blue", lwd = 2)
  abline(a =  thr, b = 1, col = "blue", lwd = 2)
  add_class_legend(delta[in_motif], thr, title = label_title)
}

# =========================
# START PDF
# =========================
pdf(out_pdf, width = 16, height = 6, colormodel = "rgb", useDingbats = FALSE)
par(mfrow = c(1, 3), mar = c(4, 4, 4, 1))
pages_written <- 0L

if (nrow(pairs) == 0) {
  plot.new(); title(main = paste0("No rows in pairs file:\n", pairs_file))
  pages_written <- 1L
} else {

  for (i in seq_len(nrow(pairs))) {

    cat("Processing:", pairs$cancer[i], pairs$patient[i], "\n")

    tryCatch({

      tf <- pairs$tumor_file[i]
      nf <- pairs$healthy_file[i]

      if (!file.exists(tf) || !file.exists(nf)) {
        plot.new()
        title(main = paste0(
          pairs$cancer[i], " | ", pairs$patient[i], "\nMissing file(s)\n",
          if (!file.exists(tf)) paste0("tumor_file: ", tf, "\n") else "",
          if (!file.exists(nf)) paste0("healthy_file: ", nf) else ""
        ))
        pages_written <- pages_written + 1L
        next
      }

      tumor  <- read_bed_gz(tf)
      normal <- read_bed_gz(nf)

      tumor_df  <- data.frame(probe = tumor[[4]],  meth_tumor  = tumor[[5]])
      normal_df <- data.frame(probe = normal[[4]], meth_normal = normal[[5]])

      merged <- merge(tumor_df, normal_df, by = "probe")
      merged <- merged[!is.na(merged$meth_tumor) & !is.na(merged$meth_normal), ]
      if (nrow(merged) == 0) {
        plot.new()
        title(main = paste0(pairs$cancer[i], " | ", pairs$patient[i], "\nMerged = 0 usable CpGs"))
        pages_written <- pages_written + 1L
        next
      }

      delta <- merged$meth_tumor - merged$meth_normal

      # =========================
      # PAGE 1 — mm0to2 (NRF1 + BANP)
      # =========================
      plot_all(merged, delta, thr,
               paste0(pairs$cancer[i], " | ", pairs$patient[i], "\nAll CpGs (mm0to2 page)"))

      in_NRF1 <- merged$probe %in% motif_cpgs$NRF1_mm0to2
      plot_motif_red(merged, delta, thr, in_NRF1, "NRF1_mm0to2 CpGs")

      in_BANP <- merged$probe %in% motif_cpgs$BANP_mm0to2
      plot_motif_red(merged, delta, thr, in_BANP, "BANP_mm0to2 CpGs")

      pages_written <- pages_written + 1L

      # =========================
      # PAGE 2 — mm0to2_noCGmm (NRF1 + BANP)
      # =========================
      plot_all(merged, delta, thr,
               paste0(pairs$cancer[i], " | ", pairs$patient[i], "\nAll CpGs (noCGmm page)"))

      in_NRF1_nc <- merged$probe %in% motif_cpgs$NRF1_noCGmm_mm0to2
      plot_motif_red(merged, delta, thr, in_NRF1_nc, "NRF1_mm0to2_noCGmm CpGs")

      in_BANP_nc <- merged$probe %in% motif_cpgs$BANP_noCGmm_mm0to2
      plot_motif_red(merged, delta, thr, in_BANP_nc, "BANP_mm0to2_noCGmm CpGs")

      pages_written <- pages_written + 1L

    }, error = function(e) {
      plot.new()
      title(main = paste0(pairs$cancer[i], " | ", pairs$patient[i], "\nERROR:\n", conditionMessage(e)))
      pages_written <- pages_written + 1L
    })
  }
}

dev.off()
cat("DONE -> ", out_pdf, " (pages:", pages_written, ")\n", sep="")
'



#############################################################################
# create an annotated file with probe_id distance and label proximal distal
#############################################################################
awk 'BEGIN{OFS="\t"} {label=(($2 < 0 ? -$2 : $2) <= 2000 ? "proximal" : "distal"); print $1, $2, label}' \
./methylation/dist_to_tss/HM450_TCGA-ACC-TCGA-OR-A5J1-01A_1_annotated_methylation_filtered_cpg_dist_tss.tsv \
> ./methylation/probes_tss_annotation.tsv


#smoothscatterplot separating proxial vs distal probes to see if we have more methylation changes in one group vs the other
Rscript -e '
suppressPackageStartupMessages({
  library(data.table)
  library(KernSmooth)
})

# =========================
# INPUTS
# =========================
pairs_file <- "./methylation/sample_pairs_files/methylation_pairs_filtered.tsv"
annot_file <- "./methylation/probes_tss_annotation.tsv"
out_pdf    <- "./results/methylation/smoothscatter/smoothScatter_mm0to2_proxdist_thr20.pdf"
thr <- 20

motif_files <- list(
  NRF1 = "./methylation/overlaps/intersected_motifs2mm_HM450/NRF1_mm0to2_intersected_methylation.bed",
  BANP = "./methylation/overlaps/intersected_motifs2mm_HM450/BANP_mm0to2_intersected_methylation.bed"
)

dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)

# =========================
# HELPERS
# =========================
pct_fmt <- function(x, n) {
  if (n == 0) return("0%")
  sprintf("%.1f%%", 100 * x / n)
}

read_bed_gz <- function(f) {
  fread(cmd = paste("zcat", shQuote(f)), header = FALSE, showProgress = FALSE)
}

sum_txt <- function(d, nm, thr) {
  n <- length(d)
  hypo  <- sum(d <= -thr, na.rm = TRUE)
  hyper <- sum(d >=  thr, na.rm = TRUE)
  unch  <- sum(abs(d) < thr, na.rm = TRUE)

  paste0(
    nm, " (n=", n, ")\n",
    "Hypo: ", hypo, " (", pct_fmt(hypo, n), ")\n",
    "Unchanged: ", unch, " (", pct_fmt(unch, n), ")\n",
    "Hyper: ", hyper, " (", pct_fmt(hyper, n), ")"
  )
}

plot_all <- function(m, d, title, thr) {
  smoothScatter(
    m$meth_normal, m$meth_tumor,
    xlab = "Healthy methylation (%)",
    ylab = "Tumor methylation (%)",
    main = title
  )

  abline(a = 0, b = 1, lwd = 2, col = "black")
  abline(a = -thr, b = 1, lwd = 2, col = "blue")
  abline(a =  thr, b = 1, lwd = 2, col = "blue")

  u <- par("usr")
  text(
    u[1] + 0.02 * (u[2] - u[1]),
    u[4] - 0.02 * (u[4] - u[3]),
    labels = sum_txt(d, "All CpGs", thr),
    adj = c(0, 1),
    cex = 0.75
  )
}

plot_motif <- function(m, d, motif_probes, title, thr) {
  is_prox <- m$probe %in% motif_probes & m$cg_proximity == "proximal"
  is_dist <- m$probe %in% motif_probes & m$cg_proximity == "distal"

  smoothScatter(
    m$meth_normal, m$meth_tumor,
    xlab = "Healthy methylation (%)",
    ylab = "Tumor methylation (%)",
    main = title
  )

  points(
    m$meth_normal[is_prox], m$meth_tumor[is_prox],
    pch = 16, cex = 0.5, col = rgb(0, 0, 1, 0.6)
  )

  points(
    m$meth_normal[is_dist], m$meth_tumor[is_dist],
    pch = 16, cex = 0.5, col = rgb(1, 0.4, 0.7, 0.6)
  )

  abline(a = 0, b = 1, lwd = 2, col = "black")
  abline(a = -thr, b = 1, lwd = 2, col = "blue")
  abline(a =  thr, b = 1, lwd = 2, col = "blue")

  u <- par("usr")

  # résumé proximal (bleu)
  text(
    u[1] + 0.02 * (u[2] - u[1]),
    u[4] - 0.02 * (u[4] - u[3]),
    labels = sum_txt(d[is_prox], "Proximal", thr),
    adj = c(0, 1),
    cex = 0.70,
    col = "blue"
  )

  # résumé distal (rose), placé à côté du bleu
  text(
    u[1] + 0.30 * (u[2] - u[1]),
    u[4] - 0.02 * (u[4] - u[3]),
    labels = sum_txt(d[is_dist], "Distal", thr),
    adj = c(0, 1),
    cex = 0.70,
    col = rgb(1, 0.4, 0.7, 1)
  )

  legend(
    "bottomright",
    legend = c(
      paste0("Proximal (n=", sum(is_prox, na.rm = TRUE), ")"),
      paste0("Distal (n=", sum(is_dist, na.rm = TRUE), ")")
    ),
    col = c(rgb(0, 0, 1, 0.8), rgb(1, 0.4, 0.7, 0.8)),
    pch = 16,
    bty = "n",
    cex = 0.8
  )
}

# =========================
# LOAD DATA
# =========================
pairs <- fread(pairs_file)
stopifnot(all(c("cancer", "patient", "tumor_file", "healthy_file") %in% names(pairs)))

annot <- fread(annot_file, header = FALSE)
if (ncol(annot) < 3) stop("Annotation file must have at least 3 columns.")
setnames(annot, c("probe", "dist_to_tss", "cg_proximity"))
annot <- unique(annot[, .(probe, cg_proximity)])

motifs <- lapply(motif_files, function(f) {
  if (!file.exists(f)) stop(paste("Missing motif file:", f))
  unique(fread(f, header = FALSE)[[4]])
})

# =========================
# PLOT
# =========================
pdf(out_pdf, width = 16, height = 6, useDingbats = FALSE)
pages_written <- 0L

for (i in seq_len(nrow(pairs))) {
  cat("Processing:", pairs$cancer[i], pairs$patient[i], "\n")
  flush.console()

  par(mfrow = c(1, 3), mar = c(4, 4, 4, 1))

  tryCatch({
    tf <- pairs$tumor_file[i]
    nf <- pairs$healthy_file[i]

    if (!file.exists(tf) || !file.exists(nf)) {
      plot.new()
      title(main = paste0(
        pairs$cancer[i], " | ", pairs$patient[i], "\nMissing file(s)"
      ))
      plot.new()
      plot.new()
      pages_written <- pages_written + 1L
      next
    }

    tumor  <- read_bed_gz(tf)
    normal <- read_bed_gz(nf)

    m <- merge(
      data.table(probe = tumor[[4]],  meth_tumor  = as.numeric(tumor[[5]])),
      data.table(probe = normal[[4]], meth_normal = as.numeric(normal[[5]])),
      by = "probe"
    )

    m <- m[!is.na(meth_tumor) & !is.na(meth_normal)]
    m <- merge(m, annot, by = "probe", all.x = TRUE)

    if (nrow(m) == 0) {
      plot.new()
      title(main = paste0(pairs$cancer[i], " | ", pairs$patient[i], "\nMerged = 0 usable CpGs"))
      plot.new()
      plot.new()
      pages_written <- pages_written + 1L
      next
    }

    d <- m$meth_tumor - m$meth_normal
    lab <- paste0(pairs$cancer[i], " | ", pairs$patient[i])

    plot_all(
      m, d,
      paste0(lab, "\nAll CpGs"),
      thr
    )

    plot_motif(
      m, d, motifs$NRF1,
      paste0(lab, "\nNRF1 motifs"),
      thr
    )

    plot_motif(
      m, d, motifs$BANP,
      paste0(lab, "\nBANP motifs"),
      thr
    )

    pages_written <- pages_written + 1L

  }, error = function(e) {
    plot.new()
    title(main = paste0(
      pairs$cancer[i], " | ", pairs$patient[i], "\nERROR:\n", conditionMessage(e)
    ))
    plot.new()
    plot.new()
    pages_written <- pages_written + 1L
  })
}

dev.off()
cat("DONE -> ", out_pdf, " (pages:", pages_written, ")\n", sep = "")
'

#################################################################################################################
# Resume methylation changes in a probes matrix for all samples in filtered pairs methylation and in 3d+4d tsv :
#################################################################################################################
Rscript -e '
library(data.table)

pairs <- fread("./methylation/sample_pairs_files/methylation_pairs_filtered.tsv")
annot <- fread("./methylation/probes_tss_annotation.tsv")
setnames(annot, names(annot)[1:2], c("probe","dist_to_tss"))
annot <- annot[, .(probe, dist_to_tss)]
annot[, cg_proximity := ifelse(abs(dist_to_tss) <= 2000, "proximal", "distal")]

read_bed <- function(f) fread(f)[, .(probe=V4, beta=as.numeric(V5))]

all_dt <- rbindlist(lapply(1:nrow(pairs), function(i){
  t <- read_bed(pairs$tumor_file[i]);   setnames(t, "beta", "beta_tumor")
  h <- read_bed(pairs$healthy_file[i]); setnames(h, "beta", "beta_healthy")
  merge(t, h, by="probe")[, .(
    probe,
    cancer = pairs$cancer[i],
    patient = pairs$patient[i],
    delta_beta = beta_tumor - beta_healthy
  )]
}))

res <- merge(all_dt, annot, by="probe", all.x=TRUE)[
  , .(n_pairs=.N, median_delta_beta=median(delta_beta, na.rm=TRUE)),
  by=.(probe, cancer, dist_to_tss, cg_proximity)
]

fwrite(res, "./methylation/median_delta_beta_per_probe_per_cancer.tsv", sep="\t")
'

# To get per cancer type all cpg proximal summary:
Rscript -e '
library(data.table)

res <- fread("./methylation/median_delta_beta_per_probe_per_cancer.tsv")

final_prox <- res[cg_proximity == "distal",
  .(
    n_cpg = .N,
    final_median_delta_beta = median(median_delta_beta, na.rm = TRUE)
  ),
  by = cancer
][order(final_median_delta_beta)]

final_prox
'

# plot only one with no threshhold to see the overall distribution of methylation changes in motifs vs non motifs for one pair of samples
Rscript -e '
suppressPackageStartupMessages({
  library(data.table)
})

# =========================
# INPUTS: set one pair here
# =========================
tumor_file   <- "./methylation/filtered_methylation/HM450_TCGA-BRCA-TCGA-E9-A1RH-01A_1_annotated_methylation_filtered.bed.gz"
healthy_file <- "./methylation/filtered_methylation/HM450_TCGA-BRCA-TCGA-E9-A1RH-11A_1_annotated_methylation_filtered.bed.gz"
pair_label   <- "BRCA | TCGA-E9-A1RH"

motif_files <- list(
  NRF1_mm0to2        = "./methylation/overlaps/intersected_motifs2mm_HM450/NRF1_mm0to2_intersected_methylation.bed",
  BANP_mm0to2        = "./methylation/overlaps/intersected_motifs2mm_HM450/BANP_mm0to2_intersected_methylation.bed",
  NRF1_noCGmm_mm0to2 = "./methylation/overlaps/intersected_motifs2mm_HM450/NRF1_mm0to2_noCGmm_intersected_methylation.bed",
  BANP_noCGmm_mm0to2 = "./methylation/overlaps/intersected_motifs2mm_HM450/BANP_mm0to2_noCGmm_intersected_methylation.bed"
)

out_pdf <- "./results/methylation/smoothscatter/smoothScatter_one_pair_mm0to2_vs_noCGmm_nothreshold.pdf"
dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)

# =========================
# LOAD MOTIF PROBES
# =========================
motif_cpgs <- lapply(motif_files, function(f){
  if(!file.exists(f)) stop(paste("Missing motif file:", f))
  unique(fread(f, header = FALSE)[[4]])
})

# =========================
# HELPERS
# =========================
read_bed_gz <- function(path) {
  cmd <- paste("zcat", shQuote(path))
  fread(cmd = cmd, header = FALSE, showProgress = FALSE)
}

plot_all <- function(merged, main_title){
  smoothScatter(
    merged$meth_normal, merged$meth_tumor,
    xlab = "Healthy methylation (%)",
    ylab = "Tumor methylation (%)",
    main = main_title
  )
  abline(0, 1, col = "black", lwd = 2)
}

plot_motif_red <- function(merged, in_motif, label_title){
  smoothScatter(
    merged$meth_normal, merged$meth_tumor,
    xlab = "Healthy methylation (%)",
    ylab = "Tumor methylation (%)",
    main = paste0(label_title, " in red")
  )
  points(
    merged$meth_normal[in_motif],
    merged$meth_tumor[in_motif],
    pch = 16, cex = 0.5, col = rgb(1, 0, 0, 0.6)
  )
  abline(0, 1, col = "black", lwd = 2)
}

# =========================
# LOAD PAIR
# =========================
if (!file.exists(tumor_file)) stop(paste("Missing tumor file:", tumor_file))
if (!file.exists(healthy_file)) stop(paste("Missing healthy file:", healthy_file))

cat("Reading tumor file...\n")
tumor <- read_bed_gz(tumor_file)

cat("Reading healthy file...\n")
normal <- read_bed_gz(healthy_file)

tumor_df  <- data.frame(probe = tumor[[4]],  meth_tumor  = tumor[[5]])
normal_df <- data.frame(probe = normal[[4]], meth_normal = normal[[5]])

merged <- merge(tumor_df, normal_df, by = "probe")
merged <- merged[!is.na(merged$meth_tumor) & !is.na(merged$meth_normal), ]

cat("Merged CpGs:", nrow(merged), "\n")
if (nrow(merged) == 0) stop("Merged table has 0 usable CpGs.")

# =========================
# PDF
# =========================
pdf(out_pdf, width = 16, height = 6, colormodel = "rgb", useDingbats = FALSE)
par(mfrow = c(1, 3), mar = c(4, 4, 4, 1))

# PAGE 1 — mm0to2
plot_all(merged, paste0(pair_label, "\nAll CpGs (mm0to2 page)"))

in_NRF1 <- merged$probe %in% motif_cpgs$NRF1_mm0to2
plot_motif_red(merged, in_NRF1, "NRF1_mm0to2 CpGs")

in_BANP <- merged$probe %in% motif_cpgs$BANP_mm0to2
plot_motif_red(merged, in_BANP, "BANP_mm0to2 CpGs")

# PAGE 2 — noCGmm
plot_all(merged, paste0(pair_label, "\nAll CpGs (noCGmm page)"))

in_NRF1_nc <- merged$probe %in% motif_cpgs$NRF1_noCGmm_mm0to2
plot_motif_red(merged, in_NRF1_nc, "NRF1_mm0to2_noCGmm CpGs")

in_BANP_nc <- merged$probe %in% motif_cpgs$BANP_noCGmm_mm0to2
plot_motif_red(merged, in_BANP_nc, "BANP_mm0to2_noCGmm CpGs")

dev.off()

cat("DONE -> ", out_pdf, "\n", sep="")
'

# HM450_TCGA-BRCA-TCGA-E9-A1RH-01A_1_annotated_methylation_filtered.bed.gz
# HM450_TCGA-BRCA-TCGA-E9-A1RH-11A_1_annotated_methylation_filtered.bed.gz
# HM450_TCGA-BLCA-TCGA-BL-A13J-01A_1_annotated_methylation_filtered.bed.gz
# HM450_TCGA-BLCA-TCGA-BL-A13J-11A_1_annotated_methylation_filtered.bed.gz
#############################################################################################################################
# SMOOTH SCATTER PLOT FOR one CANCER-CANCER SAMPLE between two cancer type 
#############################################################################################################################
Rscript -e '
suppressPackageStartupMessages({
  library(data.table)
  library(KernSmooth)
})

# =========================
# INPUTS (edit these 2!)
# =========================
tumorA_file <- "./methylation/filtered_methylation/HM450_TCGA-BLCA-TCGA-BL-A13J-11A_1_annotated_methylation_filtered.bed.gz" 
tumorB_file <- "./methylation/filtered_methylation/HM450_TCGA-BLCA-TCGA-BT-A20U-11A_1_annotated_methylation_filtered.bed.gz"  
labelA <- "BLCA | BL-A13J | 11A_1"
labelB <- "BLCA | BT-A20U | 11A_1"

thr <- 20
out_pdf <- "./results/methylation/smoothscatter/smoothScatter_HEALTHY_vs_HEALTHY_one_cancer_type.pdf"
dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)

# motif CpG sets
motif_files <- list(
  NRF1_mm0to2        = "./methylation/overlaps/intersected_motifs2mm_HM450/NRF1_mm0to2_intersected_methylation.bed",
  BANP_mm0to2        = "./methylation/overlaps/intersected_motifs2mm_HM450/BANP_mm0to2_intersected_methylation.bed",
  NRF1_noCGmm_mm0to2 = "./methylation/overlaps/intersected_motifs2mm_HM450/NRF1_mm0to2_noCGmm_intersected_methylation.bed",
  BANP_noCGmm_mm0to2 = "./methylation/overlaps/intersected_motifs2mm_HM450/BANP_mm0to2_noCGmm_intersected_methylation.bed"
)

# =========================
# HELPERS
# =========================
pct_fmt <- function(x, n) {
  if (n == 0) return("0%")
  sprintf("%.1f%%", 100 * x / n)
}

add_class_legend <- function(delta_vec, thr, title = "All CpGs") {
  n <- length(delta_vec)
  hypo <- sum(delta_vec <= -thr, na.rm = TRUE)
  hyper <- sum(delta_vec >=  thr, na.rm = TRUE)
  unch <- sum(abs(delta_vec) < thr, na.rm = TRUE)

  lab <- paste0(
    title, "\n",
    "HYPO (Delta≤-", thr, "): ", pct_fmt(hypo, n), "\n",
    "Unchanged (|Delta|<", thr, "): ", pct_fmt(unch, n), "\n",
    "HYPER (Delta≥+", thr, "): ", pct_fmt(hyper, n)
  )

  usr <- par("usr")
  x <- usr[1] + 0.02 * (usr[2] - usr[1])
  y <- usr[4] - 0.02 * (usr[4] - usr[3])
  text(x, y, labels = lab, adj = c(0, 1), cex = 0.85)
}

read_bed_gz <- function(path) {
  if(!file.exists(path)) stop(paste("Missing:", path))
  fread(cmd = paste("zcat", shQuote(path)), header = FALSE, showProgress = FALSE)
}

plot_all <- function(merged, delta, thr, main_title){
  smoothScatter(
    merged$meth_A, merged$meth_B,
    xlab = paste0(labelA, " methylation (%)"),
    ylab = paste0(labelB, " methylation (%)"),
    main = main_title
  )
  abline(0, 1, col = "black", lwd = 2)
  abline(a = -thr, b = 1, col = "blue", lwd = 2)
  abline(a =  thr, b = 1, col = "blue", lwd = 2)
  add_class_legend(delta, thr, title = "All CpGs")
}

plot_motif_red <- function(merged, delta, thr, in_motif, label_title){
  smoothScatter(
    merged$meth_A, merged$meth_B,
    xlab = paste0(labelA, " methylation (%)"),
    ylab = paste0(labelB, " methylation (%)"),
    main = paste0(label_title, " in red (|DeltaBeta| ≥ ", thr, "%)")
  )
  is_red <- in_motif & (abs(delta) >= thr)
  points(
    merged$meth_A[is_red],
    merged$meth_B[is_red],
    pch = 16, cex = 0.5, col = rgb(1, 0, 0, 0.6)
  )
  abline(0, 1, col = "black", lwd = 2)
  abline(a = -thr, b = 1, col = "blue", lwd = 2)
  abline(a =  thr, b = 1, col = "blue", lwd = 2)
  add_class_legend(delta[in_motif], thr, title = label_title)
}

# =========================
# LOAD MOTIF SETS
# =========================
motif_cpgs <- lapply(motif_files, function(f){
  if(!file.exists(f)) stop(paste("Missing motif file:", f))
  unique(fread(f, header = FALSE, showProgress = FALSE)[[4]])
})

# =========================
# LOAD TWO TUMORS + MERGE
# =========================
cat("[INFO] Reading tumor A:", tumorA_file, "\n")
A <- read_bed_gz(tumorA_file)
cat("[INFO] Reading tumor B:", tumorB_file, "\n")
B <- read_bed_gz(tumorB_file)

A_dt <- data.table(probe=A[[4]], meth_A=as.numeric(A[[5]]))
B_dt <- data.table(probe=B[[4]], meth_B=as.numeric(B[[5]]))

merged <- merge(A_dt, B_dt, by="probe")
merged <- merged[!is.na(meth_A) & !is.na(meth_B)]
if(nrow(merged) == 0) stop("Merged = 0 usable CpGs (no overlap or all NA)")

delta <- merged$meth_B - merged$meth_A

# =========================
# PLOT PDF (2 pages, each has 3 panels)
# =========================
pdf(out_pdf, width = 16, height = 6, colormodel = "rgb", useDingbats = FALSE)
par(mfrow = c(1, 3), mar = c(4, 4, 4, 1))

# Page 1: mm0to2
plot_all(merged, delta, thr, paste0(labelA, "  vs  ", labelB, "\nAll CpGs (mm0to2 page)"))
plot_motif_red(merged, delta, thr, merged$probe %in% motif_cpgs$NRF1_mm0to2, "NRF1_mm0to2 CpGs")
plot_motif_red(merged, delta, thr, merged$probe %in% motif_cpgs$BANP_mm0to2, "BANP_mm0to2 CpGs")

# Page 2: noCGmm_mm0to2
plot_all(merged, delta, thr, paste0(labelA, "  vs  ", labelB, "\nAll CpGs (noCGmm page)"))
plot_motif_red(merged, delta, thr, merged$probe %in% motif_cpgs$NRF1_noCGmm_mm0to2, "NRF1_mm0to2_noCGmm CpGs")
plot_motif_red(merged, delta, thr, merged$probe %in% motif_cpgs$BANP_noCGmm_mm0to2, "BANP_mm0to2_noCGmm CpGs")

dev.off()
cat("DONE -> ", out_pdf, "\n", sep="")
'

###########################################################################################################################
# plot heatmap of % of changed CpGs in each motif∩peaks set comparing samples from same cancer type or other or same organ..
###########################################################################################################################
# create a metadata file with the pairs and their categories (same cancer type, same organ, different organ)
Rscript -e '
library(data.table)
pairs <- fread("./methylation/sample_pairs_files/methylation_pairs_filtered.tsv")

# small cancer -> organ mapping
cancer_to_organ <- c(
  "TCGA-BLCA" = "Bladder",
  "TCGA-BRCA" = "Breast",
  "TCGA-COAD" = "Colon",
  "TCGA-READ" = "Colon",
  "TCGA-GBM"  = "Brain",
  "TCGA-LGG"  = "Brain",
  "TCGA-HNSC" = "HeadNeck",
  "TCGA-KICH" = "Kidney",
  "TCGA-KIRC" = "Kidney",
  "TCGA-KIRP" = "Kidney",
  "TCGA-LIHC" = "Liver",
  "TCGA-LUAD" = "Lung",
  "TCGA-LUSC" = "Lung",
  "TCGA-OV"   = "Ovary",
  "TCGA-PAAD" = "Pancreas",
  "TCGA-PRAD" = "Prostate",
  "TCGA-SKCM" = "Skin",
  "TCGA-STAD" = "Stomach",
  "TCGA-THCA" = "Thyroid",
  "TCGA-UCEC" = "Uterus"
)

tumor <- pairs[, .(
  cancer = cancer,
  patient = patient,
  type = "Tumor",
  file = tumor_file
)]

healthy <- pairs[, .(
  cancer = cancer,
  patient = patient,
  type = "Healthy",
  file = healthy_file
)]

meta <- rbindlist(list(tumor, healthy), use.names = TRUE)

# extract sample_id from file name
meta[, sample_id := sub(".*HM450_[^-]+-[^-]+-(TCGA-[A-Z0-9-]+)_[^/]*$", "\\1", file)]

# if the above does not work with your file names, use this instead:
# meta[, sample_id := sub(".*HM450_[^-]+-(TCGA-[A-Z0-9-]+-[0-9]{2}[A-Z])_1_.*", "\\1", file)]

meta[, organ := cancer_to_organ[cancer]]

# reorder columns
meta <- meta[, .(sample_id, cancer, organ, type, file, patient)]

# remove missing files / duplicates
meta <- unique(meta[!is.na(file) & file != ""])

fwrite(meta, "./methylation/sample_metadata_for_comparisons.tsv", sep = "\t")
'

# plot the heatmap using the output from the previous code 
Rscript -e '
suppressPackageStartupMessages({
  library(data.table)
  library(pheatmap)
})

meta_file <- "./methylation/sample_pairs_files/methylation_pairs_filtered.tsv"
delta_thr <- 20

motifs <- list(
  NRF1 = "./motifs/overlaps/intersected_motifs2mm_HM450/NRF1_mm0to2_noCGmm_probeXmotif.tsv",
  BANP = "./motifs/overlaps/intersected_motifs2mm_HM450/BANP_mm0to2_noCGmm_probeXmotif.tsv"
)

outdir <- "./results/methylation/pairwise_cancer_heatmaps_FULL"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

find_probe_col <- function(dt) {
  cand <- c("probe","Probe","probe_id","Probe_ID","cg","cg_id","ID_REF","IlmnID","V4")
  hit <- cand[cand %in% names(dt)]
  if (length(hit) == 0) stop("No probe column found: ", paste(names(dt), collapse=", "))
  hit[1]
}

find_beta_col <- function(dt) {
  cand <- c("beta","Beta","beta_value","Beta_value","b","V5")
  hit <- cand[cand %in% names(dt)]
  if (length(hit) == 0) stop("No beta column found: ", paste(names(dt), collapse=", "))
  hit[1]
}

read_meth <- function(f) {
  x <- fread(cmd = paste("zcat", shQuote(f)))
  probe_col <- find_probe_col(x)
  beta_col  <- find_beta_col(x)
  x <- x[, .(probe = get(probe_col), beta = suppressWarnings(as.numeric(get(beta_col))))]
  x <- x[!is.na(probe) & !is.na(beta)]
  x
}

meta <- fread(meta_file)

tumor_by_cancer   <- split(meta$tumor_file, meta$cancer)
healthy_by_cancer <- split(meta$healthy_file, meta$cancer)

cancers <- sort(intersect(names(tumor_by_cancer), names(healthy_by_cancer)))
tumor_by_cancer   <- lapply(tumor_by_cancer[cancers], unique)
healthy_by_cancer <- lapply(healthy_by_cancer[cancers], unique)

cat("Cancers included:\n")
print(cancers)

all_files <- unique(c(unlist(tumor_by_cancer), unlist(healthy_by_cancer)))
cat("Loading", length(all_files), "methylation files...\n")

methylation_data <- setNames(vector("list", length(all_files)), all_files)
for (k in seq_along(all_files)) {
  f <- all_files[k]
  cat("  ", k, "/", length(all_files), " ", basename(f), "\n", sep = "")
  methylation_data[[f]] <- read_meth(f)
}

for (nm in names(motifs)) {
  cat("\n============================\n")
  cat("Processing ", nm, "\n", sep = "")
  cat("============================\n")

  motif_dt <- fread(motifs[[nm]])
  motif_probe_col <- find_probe_col(motif_dt)
  probes_keep <- unique(motif_dt[[motif_probe_col]])
  probes_keep <- probes_keep[!is.na(probes_keep)]

  cat("Motif probes: ", length(probes_keep), "\n", sep = "")

  cat("Building tumor vs healthy matrix...\n")
  mat_TH <- matrix(NA_real_, length(cancers), length(cancers),dimnames = list(cancers, cancers))

  for (i in seq_along(cancers)) {
    ci <- cancers[i]
    cat("\n  Row cancer (tumor): ", ci, "\n", sep = "")
    files_i <- tumor_by_cancer[[ci]]

    for (j in seq_along(cancers)) {
      cj <- cancers[j]
      files_j <- healthy_by_cancer[[cj]]

      cat("    vs column cancer (healthy): ", cj,
          " | comparing ", length(files_i), " tumor samples vs ",
          length(files_j), " healthy samples\n", sep = "")

      vals <- c()

      for (f1 in files_i) {
        x1 <- methylation_data[[f1]][probe %in% probes_keep]

        for (f2 in files_j) {
          x2 <- methylation_data[[f2]][probe %in% probes_keep]

          m <- merge(x1, x2, by = "probe", suffixes = c(".1", ".2"))
          m <- m[!is.na(beta.1) & !is.na(beta.2)]
          if (nrow(m) == 0) next

          vals <- c(vals, mean(abs(m$beta.1 - m$beta.2) >= delta_thr) * 100)
        }
      }

      mat_TH[i, j] <- if (length(vals)) median(vals, na.rm = TRUE) else NA_real_

      cat("      -> cell value: ", ifelse(is.na(mat_TH[i, j]), "NA", sprintf("%.2f", mat_TH[i, j])), "\n", sep = "")
    }
  }

  fwrite(as.data.table(mat_TH, keep.rownames = "Cancer"),file.path(outdir, paste0(nm, "_tumor_vs_healthy_matrix.tsv")),sep = "\t")

  pdf(file.path(outdir, paste0(nm, "_tumor_vs_healthy_heatmap.pdf")), width = 12, height = 10)
  pheatmap(mat_TH,
           main = paste0(nm, " tumor vs healthy\nMedian % CpGs with absolute beta difference >= 20"),
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           display_numbers = TRUE,
           number_format = "%.1f",
           na_col = "grey90")
  dev.off()

  cat("\nSaved tumor vs healthy outputs for ", nm, "\n", sep = "")

  cat("\nBuilding tumor vs tumor matrix...\n")
  mat_TT <- matrix(NA_real_, length(cancers), length(cancers),dimnames = list(cancers, cancers))

  for (i in seq_along(cancers)) {
    ci <- cancers[i]
    cat("\n  Row cancer (tumor): ", ci, "\n", sep = "")
    files_i <- tumor_by_cancer[[ci]]

    for (j in seq_along(cancers)) {
      cj <- cancers[j]
      files_j <- tumor_by_cancer[[cj]]

      cat("    vs column cancer (tumor): ", cj,
          " | comparing ", length(files_i), " tumor samples vs ",
          length(files_j), " tumor samples\n", sep = "")

      vals <- c()

      for (f1 in files_i) {
        x1 <- methylation_data[[f1]][probe %in% probes_keep]

        for (f2 in files_j) {
          x2 <- methylation_data[[f2]][probe %in% probes_keep]

          m <- merge(x1, x2, by = "probe", suffixes = c(".1", ".2"))
          m <- m[!is.na(beta.1) & !is.na(beta.2)]
          if (nrow(m) == 0) next

          vals <- c(vals, mean(abs(m$beta.1 - m$beta.2) >= delta_thr) * 100)
        }
      }

      mat_TT[i, j] <- if (length(vals)) median(vals, na.rm = TRUE) else NA_real_

      cat("      -> cell value: ", ifelse(is.na(mat_TT[i, j]), "NA", sprintf("%.2f", mat_TT[i, j])), "\n", sep = "")
    }
  }

  fwrite(as.data.table(mat_TT, keep.rownames = "Cancer"),file.path(outdir, paste0(nm, "_tumor_vs_tumor_matrix.tsv")),sep = "\t")

  pdf(file.path(outdir, paste0(nm, "_tumor_vs_tumor_heatmap.pdf")), width = 12, height = 10)
  pheatmap(mat_TT,
           main = paste0(nm, " tumor vs tumor\nMedian % CpGs with absolute beta difference >= 20"),
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           display_numbers = TRUE,
           number_format = "%.1f",
           na_col = "grey90")
  dev.off()

  cat("\nSaved tumor vs tumor outputs for ", nm, "\n", sep = "")
}

cat("\nFULL run done.\n")
'



###########################################################################################################################
# QUANTIFICATION OF % OF CpGs WITH |DELTA| >= 20% IN NRF1 AND BANP MOTIFS (mm0to2 vs noCGmm mm0to2) AND IN ATAC PEAKS for only pairs with healthy vs tumor samples (filtered pairs file)
###########################################################################################################################
Rscript -e '
suppressPackageStartupMessages({library(data.table)})

pairs_file <- "./methylation/sample_pairs_files/methylation_pairs_filtered.tsv"
thr <- 20

# motif ∩ peaks CpG sets (from  bedtools intersect outputs)
motif_peak_files <- list(
  NRF1_mm0to2        = "./methylation/overlaps/intersected_motifs2mm_HM450_peaks/NRF1_mm0to2_intersected_methylation.bed",
  NRF1_mm0to2_noCGmm = "./methylation/overlaps/intersected_motifs2mm_HM450_peaks/NRF1_mm0to2_noCGmm_intersected_methylation.bed",
  BANP_mm0to2        = "./methylation/overlaps/intersected_motifs2mm_HM450_peaks/BANP_mm0to2_intersected_methylation.bed",
  BANP_mm0to2_noCGmm = "./methylation/overlaps/intersected_motifs2mm_HM450_peaks/BANP_mm0to2_noCGmm_intersected_methylation.bed"
)

out_tsv <- "./results/methylation/count_CpGs_motifPeaks_delta20_per_pair.tsv"
dir.create(dirname(out_tsv), recursive=TRUE, showWarnings=FALSE)

read_bed_gz <- function(path) fread(cmd=paste("zcat", shQuote(path)), header=FALSE, showProgress=FALSE)

pairs <- fread(pairs_file)

# load probe IDs for each motif∩peaks set (column 4 = probe id in PROBES bed)
motif_peak_cpgs <- lapply(motif_peak_files, function(f){
  if(!file.exists(f)) stop(paste("Missing:", f))
  unique(fread(f, header=FALSE)[[4]])
})

res <- vector("list", nrow(pairs))
for(i in seq_len(nrow(pairs))){

  tf <- pairs$tumor_file[i]
  nf <- pairs$healthy_file[i]
  if(!file.exists(tf) || !file.exists(nf)){
    res[[i]] <- data.table(cancer=pairs$cancer[i], patient=pairs$patient[i])
    next
  }
  tumor  <- read_bed_gz(tf)
  normal <- read_bed_gz(nf)
  tumor_dt  <- data.table(probe=tumor[[4]],  meth_tumor=as.numeric(tumor[[5]]))
  normal_dt <- data.table(probe=normal[[4]], meth_normal=as.numeric(normal[[5]]))
  merged <- merge(tumor_dt, normal_dt, by="probe")
  merged <- merged[!is.na(meth_tumor) & !is.na(meth_normal)]
  merged[, delta := meth_tumor - meth_normal]
  merged[, big := abs(delta) >= thr]

  row <- data.table(cancer=pairs$cancer[i], patient=pairs$patient[i])

  # counts for each motif∩peaks set
  for(nm in names(motif_peak_cpgs)){
    row[[paste0(nm, "_nCpGs")]] <- sum(merged$probe %in% motif_peak_cpgs[[nm]])
    row[[paste0(nm, "_nDelta20")]] <- sum((merged$probe %in% motif_peak_cpgs[[nm]]) & merged$big)}

  # optional: overall context
  row[["AllCpGs_nDelta20"]] <- sum(merged$big)

  res[[i]] <- row
}

res_dt <- rbindlist(res, fill=TRUE)
fwrite(res_dt, out_tsv, sep="\\t")

# totals across all pairs
num_cols <- setdiff(names(res_dt), c("cancer","patient"))
tot <- res_dt[, lapply(.SD, function(x) sum(as.numeric(x), na.rm=TRUE)), .SDcols=num_cols]
fwrite(tot, sub("\\\\.tsv$", "_TOTALS.tsv", out_tsv), sep="\\t")

cat("DONE ->", out_tsv, "\\n")
'
# boxplot for the % of CpGs with |delta| >= 20% in each motif∩peaks set across all pairs.
Rscript -e '
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(Polychrome)
})

in_tsv  <- "./results/methylation/count_CpGs_motifPeaks_delta20_per_pair.tsv"
out_pdf <- "./results/methylation/boxplot_per_cancer_motifPeaks_delta20.pdf"
dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)

df <- read_tsv(in_tsv, show_col_types = FALSE)

# Keep only the columns we need
keep_cols <- c(
  "cancer","patient",
  "NRF1_mm0to2_nCpGs","NRF1_mm0to2_nDelta20",
  "NRF1_mm0to2_noCGmm_nCpGs","NRF1_mm0to2_noCGmm_nDelta20",
  "BANP_mm0to2_nCpGs","BANP_mm0to2_nDelta20",
  "BANP_mm0to2_noCGmm_nCpGs","BANP_mm0to2_noCGmm_nDelta20"
)
df <- df %>% select(any_of(keep_cols))

# Long format with nCpGs + nDelta20 directly 
long <- df %>%
  pivot_longer(
    cols = -c(cancer, patient),
    names_to = c("TF", "nocg", ".value"),
    names_pattern = "(NRF1|BANP)_mm0to2(_noCGmm)?_(nCpGs|nDelta20)",
    values_drop_na = TRUE
  ) %>%
  mutate(
    set = ifelse(is.na(nocg) | nocg == "", "mm0to2", "mm0to2_noCGmm"),
    nCpGs = as.numeric(nCpGs),
    nDelta20 = as.numeric(nDelta20),
    pct = 100 * nDelta20 / nCpGs
  ) %>%
  mutate(
    panel = paste(TF, set, sep = " | "),
    cancer = factor(cancer)
  )

# ---- Polychrome palette  ----
n_c <- nlevels(long$cancer)
poly_cols <- createPalette(n_c, seedcolors = c("#3B4992", "#EE0000", "#008B45"))
cancer_cols <- setNames(poly_cols, levels(long$cancer))

make_panel <- function(dat, title_txt){
  ggplot(dat, aes(x = cancer, y = pct, fill = cancer, color = cancer)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.25) +
    geom_jitter(width = 0.2, size = 1.8, alpha = 0.8) +
    coord_flip() +
    theme_bw() +
    scale_fill_manual(values = cancer_cols) +
    scale_color_manual(values = cancer_cols) +
    labs(
      x = NULL,
      y = "% CpGs with |DeltaBeta| ≥ 20",
      title = title_txt,
      caption = paste0(
        "Each dot represents one patient with matched tumor and normal samples.\n",
        "Y-axis: percentage of CpGs within the TF motif ∩ ATAC-peak set\n",
        "showing an absolute methylation difference |DeltaBeta| ≥ 20% (tumor vs normal).\n",
        "Boxplots summarize the distribution across patients for each cancer type."
      )
    ) +
    theme(
      plot.title   = element_text(size = 11, face = "bold"),
      plot.caption = element_text(size = 9, hjust = 0),
      legend.position = "none"
    )
}

p1 <- make_panel(filter(long, TF=="NRF1", set=="mm0to2"),        "NRF1 | mm0to2 (motif ∩ peaks)")
p2 <- make_panel(filter(long, TF=="NRF1", set=="mm0to2_noCGmm"), "NRF1 | mm0to2_noCGmm (motif ∩ peaks)")
p3 <- make_panel(filter(long, TF=="BANP", set=="mm0to2"),        "BANP | mm0to2 (motif ∩ peaks)")
p4 <- make_panel(filter(long, TF=="BANP", set=="mm0to2_noCGmm"), "BANP | mm0to2_noCGmm (motif ∩ peaks)")

# One page PDF with 4 plots (2x2)
pdf(out_pdf, width = 16, height = 10, useDingbats = FALSE)
par(mfrow = c(2,2), mar = c(4, 8, 3, 1))
print(p1)
print(p2)
print(p3)
print(p4)
dev.off()

cat("DONE -> ", out_pdf, "\n", sep="")
'

# write out the unique CpG probe IDs that have |delta| >= 20% in each motif∩peaks set across all pairs (one file per set). Cg in the output file are in atleast one pair with |delta| >= 20% and are in the motif∩peaks set. 
Rscript -e 'q
suppressPackageStartupMessages({library(data.table)})

pairs_file <- "./methylation/sample_pairs_files/methylation_pairs_filtered.tsv"
thr <- 20

motif_peak_files <- list(
  NRF1_mm0to2        = "./methylation/overlaps/intersected_motifs2mm_HM450_peaks/NRF1_mm0to2_intersected_methylation.bed",
  NRF1_mm0to2_noCGmm = "./methylation/overlaps/intersected_motifs2mm_HM450_peaks/NRF1_mm0to2_noCGmm_intersected_methylation.bed",
  BANP_mm0to2        = "./methylation/overlaps/intersected_motifs2mm_HM450_peaks/BANP_mm0to2_intersected_methylation.bed",
  BANP_mm0to2_noCGmm = "./methylation/overlaps/intersected_motifs2mm_HM450_peaks/BANP_mm0to2_noCGmm_intersected_methylation.bed"
)

out_dir <- "./methylation/delta20_cg_lists"
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

read_bed_gz <- function(path){
  fread(cmd=paste("zcat", shQuote(path)),
        header=FALSE, showProgress=FALSE)
}

pairs <- fread(pairs_file)

# load CpG probe IDs for each motif∩peaks set (col 4 = probe id)
motif_peak_cpgs <- lapply(motif_peak_files, function(f){
  if(!file.exists(f)) stop(paste("Missing:", f))
  unique(fread(f, header=FALSE, showProgress=FALSE)[[4]])
})

# store CpGs passing threshold per set (collect in lists for speed)
hit_sets <- setNames(vector("list", length(motif_peak_cpgs)), names(motif_peak_cpgs))
for(nm in names(hit_sets)) hit_sets[[nm]] <- character(0)

n_pairs <- nrow(pairs)
for(i in seq_len(n_pairs)){

  if(i %% 20 == 0) cat("[", i, "/", n_pairs, "]\n")

  tf <- pairs$tumor_file[i]
  nf <- pairs$healthy_file[i]
  if(!file.exists(tf) || !file.exists(nf)) next

  tumor  <- read_bed_gz(tf)
  normal <- read_bed_gz(nf)

  # BED columns: chr, start, end, probe, beta
  tumor_dt  <- data.table(probe=tumor[[4]],  meth_tumor = as.numeric(tumor[[5]]))
  normal_dt <- data.table(probe=normal[[4]], meth_normal= as.numeric(normal[[5]]))

  merged <- merge(tumor_dt, normal_dt, by="probe")
  merged <- merged[!is.na(meth_tumor) & !is.na(meth_normal)]
  merged[, delta := meth_tumor - meth_normal]

  # CpGs with |Δβ| >= thr
  big_probes <- merged[abs(delta) >= thr, probe]
  if(length(big_probes) == 0) next

  # add to each set (use %chin% for speed)
  for(nm in names(motif_peak_cpgs)){
    hit_sets[[nm]] <- c(hit_sets[[nm]], big_probes[big_probes %chin% motif_peak_cpgs[[nm]]])
  }
}

# write unique IDs per set
for(nm in names(hit_sets)){
  ids <- sort(unique(hit_sets[[nm]]))
  out <- file.path(out_dir, paste0(nm, "_delta20_cg_ids.txt"))
  fwrite(data.table(probe=ids), out, sep="\t", col.names=FALSE)
  cat("WROTE", length(ids), "CpGs ->", out, "\n")
}
'
# add the motif associated with each cg in the above output files (one file per set). This is for downstream use when we want to check the motif associated with each cg in the list of cg with |delta| >= 20% and in the motif∩peaks set. The output file will have two columns: probe and motif and chr start end.
PROBES="./methylation/annotated_methylation_data_probes_filtered.bed"
MOTIF_DIR="./motifs"
OUTDIR="./motifs/overlaps/intersected_motifs2mm_HM450"
mkdir -p "$OUTDIR"

# ensure probes are sorted once
sort -k1,1 -k2,2n "$PROBES" > "$OUTDIR/probes.sorted.bed"

for MF in "$MOTIF_DIR"/*mm0to2*.bed.gz; do
  echo "[INFO] Processing $MF ..."
  base=$(basename "$MF" .bed.gz)
  # Decompress + sort motif bed (bedtools likes sorted input)
  zcat "$MF" \
    | sort -k1,1 -k2,2n \
    | bedtools intersect -wa -wb \
        -a "$OUTDIR/probes.sorted.bed" \
        -b stdin \
    > "$OUTDIR/${base}_probeXmotif.tsv"
done

# Build a matrix  of presence/absence (0/1) of each cg (rows) with |delta| >= 20% in each cancer-patient pair (columns) for the CpGs in the list of cg with |delta| >= 20% and in the motif∩peaks set. We will have one matrix per cancer type (with all patients of that cancer as columns).
Rscript -e '
suppressPackageStartupMessages({library(data.table)})

pairs_file <- "./methylation/sample_pairs_files/methylation_pairs_filtered.tsv"
thr <- 20

cg_list_file <- "./methylation/delta20_cg_lists/BANP_mm0to2_noCGmm_delta20_cg_ids.txt"
set_name <- sub("_delta20_cg_ids[.]txt$", "", basename(cg_list_file))

out_dir <- file.path("./methylation/delta20_cg_lists/cg_presence_matrices", set_name)
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

read_bed_gz <- function(path){
  fread(cmd=paste("zcat", shQuote(path)), header=FALSE, showProgress=FALSE)
}

pairs <- fread(pairs_file)
stopifnot(all(c("cancer","patient","tumor_file","healthy_file") %in% names(pairs)))

cgs <- fread(cg_list_file, header=FALSE, col.names=c("probe"))
cgs <- unique(cgs)
cat("[INFO] CpGs in list:", nrow(cgs), "\n")

for (C in sort(unique(pairs$cancer))) {

  subp <- pairs[cancer == C]
  if (nrow(subp) == 0) next
  cat("[CANCER]", C, "pairs:", nrow(subp), "\n")

  long_list <- vector("list", nrow(subp))
  k <- 0L

  for (i in seq_len(nrow(subp))) {

    pat <- subp$patient[i]
    tf <- subp$tumor_file[i]
    nf <- subp$healthy_file[i]
    if (!file.exists(tf) || !file.exists(nf)) next

    tumor  <- read_bed_gz(tf)
    normal <- read_bed_gz(nf)

    tumor_dt  <- data.table(probe=as.character(tumor[[4]]),  meth_tumor = as.numeric(tumor[[5]]))
    normal_dt <- data.table(probe=as.character(normal[[4]]), meth_normal= as.numeric(normal[[5]]))

    merged <- merge(tumor_dt, normal_dt, by="probe")
    merged <- merged[!is.na(meth_tumor) & !is.na(meth_normal)]
    if (nrow(merged) == 0) next

    merged[, delta := meth_tumor - meth_normal]

    merged <- merged[probe %chin% cgs$probe]
    if (nrow(merged) == 0) next

    hits <- merged[abs(delta) >= thr, .(probe)]
    if (nrow(hits) == 0) next

    hits[, `:=`(patient = pat, value = 1L)]
    k <- k + 1L
    long_list[[k]] <- hits
  }

  long <- rbindlist(long_list[seq_len(k)], fill=TRUE)

  if (nrow(long) == 0) {
    mat <- copy(cgs)
    out_file <- file.path(out_dir, paste0("matrix_", C, ".tsv"))
    fwrite(mat, out_file, sep="\t")
    cat("[WARN] No hits found for", C, "-> wrote probe-only file\n")
    next
  }

  # build CpG x patient matrix (fill missing with 0 directly)
  mat <- dcast(long, probe ~ patient, value.var="value", fill=0L)

  # ensure all CpGs are present as rows
  mat <- merge(cgs, mat, by="probe", all.x=TRUE)
  for (j in 2:ncol(mat)) set(mat, which(is.na(mat[[j]])), j, 0L)

  out_file <- file.path(out_dir, paste0("matrix_", C, ".tsv"))
  fwrite(mat, out_file, sep="\t")
  cat("[DONE] Wrote", out_file, " (rows:", nrow(mat), " cols:", ncol(mat), ")\n")
}
'
# plot heatmap of the above matrices (one per cancer type) to visualize the presence/absence of each cg with |delta| >= 20% in each patient for the CpGs in the list of cg with |delta| >= 20% and in the motif∩peaks set. 
# This script does : For each CpG-set folder, it reads each cancer’s CpG×patient 0/1 matrix, keeps CpGs present at least once, optionally keeps the top 521 most frequent CpGs, and writes heatmaps for all cancers into a single PDF.
Rscript -e '
library(data.table)
library(pheatmap)

in_root <- "./methylation/delta20_cg_lists/cg_presence_matrices"
out_dir <- "./results/methylation/heatmap_delta20_cg_motifnpeak"
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

wrap_label <- function(x, width=35) paste(strwrap(x, width=width), collapse="\n")
TOP_N <- 521

for (setdir in list.dirs(in_root, recursive=FALSE, full.names=TRUE)) {

  set_name <- basename(setdir)

  # NO backslash regex
  mats <- list.files(setdir, pattern="^matrix_.*[.]tsv$", full.names=TRUE)
  if (length(mats) == 0) next

  pdf(file.path(out_dir, paste0("HEATMAP_", set_name, "_per_cancer.pdf")),
      width=24, height=16)

  for (f in mats) {

    # NO backslash regex
    cancer <- sub("^matrix_", "", sub("[.]tsv$", "", basename(f)))

    dt <- fread(f)

    cg_ids <- dt[[1]]
    m <- as.matrix(dt[, -1, with=FALSE])
    colnames(m) <- names(dt)[-1]
    rownames(m) <- cg_ids
    storage.mode(m) <- "numeric"
    m[is.na(m)] <- 0

    if (ncol(m) == 0 || nrow(m) == 0) next

    keep <- rowSums(m) > 0
    m <- m[keep, , drop=FALSE]
    cg_ids <- cg_ids[keep]
    if (nrow(m) == 0) next

    rs <- rowSums(m)
    if (nrow(m) > TOP_N) {
      o <- order(rs, decreasing=TRUE)[1:TOP_N]
      m <- m[o, , drop=FALSE]
      cg_ids <- cg_ids[o]
    }

    n_samples <- ncol(m)
    n_cgs <- nrow(m)
    lab <- vapply(cg_ids, wrap_label, character(1), width=35)

    pheatmap(
      m,
      color = c("skyblue", "yellow"),
      breaks = c(-0.5, 0.5, 1.5),
      labels_row = lab,
      show_rownames = TRUE,
      show_colnames = TRUE,
      fontsize_row = 2,
      fontsize_col = 3.7,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      border_color = "grey60",
      legend_breaks = c(0, 1),
      legend_labels = c("absent (0)", "present (1)"),
      main = paste0(set_name, " — ", cancer,
                    " | patients=", n_samples,
                    " | CpGs in motifnpeaks=", n_cgs)
    )
  }

  dev.off()
}
'

# Combine all matrixes into one big matrix (rows = unique CpGs across all cancers, cols = all patients across all cancers) for the CpGs in the list of cg with |delta| >= 20% and in the motif∩peaks set. Also write a separate metadata file with patient and cancer info for each column in the big matrix.
Rscript -e '
suppressPackageStartupMessages({library(data.table)})

pairs_file <- "./methylation/sample_pairs_files/methylation_pairs_filtered.tsv"
out_meta   <- "./methylation/delta20_cg_lists/cg_presence_matrices/NRF1_mm0to2_noCGmm/big_matrix_colmeta.tsv"

pairs <- fread(pairs_file)
meta <- unique(pairs[, .(patient, cancer)])
setorder(meta, cancer, patient)

fwrite(meta, out_meta, sep="\t")
cat("[DONE] wrote -> ", out_meta, "\n", sep="")
'

Rscript -e '
suppressPackageStartupMessages({library(data.table)})

set_dir <- "./methylation/delta20_cg_lists/cg_presence_matrices/NRF1_mm0to2_noCGmm"
set_name <- basename(set_dir)

out_gz <- file.path("./methylation/delta20_cg_lists/cg_presence_matrices/NRF1_mm0to2_noCGmm", paste0("big_matrix_", set_name, ".tsv.gz"))
tmp_tsv <- sub("[.]gz$", "", out_gz)
dir.create("./methylation/delta20_cg_lists/cg_presence_matrices/NRF1_mm0to2_noCGmm", recursive=TRUE, showWarnings=FALSE)

files <- list.files(set_dir, pattern="^matrix_.*[.]tsv$", full.names=TRUE)
if(length(files) == 0) stop("No matrix files found")

long_list <- vector("list", length(files))

for(i in seq_along(files)){
  f <- files[[i]]
  dt <- fread(f)
  if(ncol(dt) < 2) next

  mdt <- melt(dt,
              id.vars=names(dt)[1],
              variable.name="patient",
              value.name="value")
  setnames(mdt, names(mdt)[1], "probe")
  mdt[, patient := as.character(patient)]
  long_list[[i]] <- mdt[, .(probe, patient, value)]
}

long <- rbindlist(long_list, fill=TRUE)
mat <- dcast(long, probe ~ patient, value.var="value", fill=0L)

# write then gzip (compatible with older data.table)
fwrite(mat, tmp_tsv, sep="\\t")
system(paste("gzip -f", shQuote(tmp_tsv)))

cat("[DONE] big matrix -> ", out_gz, "\\n", sep="")
'
# plot heatmap of the big matrix (rows = unique CpGs across all cancers, cols = all patients across all cancers) for the CpGs in the list of cg with |delta| >= 20% and in the motif∩peaks set. The columns (patients) should be renamed to include cancer type info from the metadata file.
Rscript -e '
library(data.table)
library(pheatmap)

in_root <- "./methylation/delta20_cg_lists/cg_presence_matrices"
out_dir <- "./results/methylation/heatmap_delta20_cg_motifnpeak"
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

wrap_label <- function(x, width=35) paste(strwrap(x, width=width), collapse="\n")
TOP_N <- 521

for (setdir in list.dirs(in_root, recursive=FALSE, full.names=TRUE)) {
  set_name <- basename(setdir)
  # read the big matrix
  big_file <- file.path(setdir, paste0("big_matrix_", set_name, ".tsv.gz"))
  if (!file.exists(big_file)) next
  dt <- fread(cmd=paste("zcat", shQuote(big_file)))
  cg_ids <- dt[[1]]
  m <- as.matrix(dt[, -1, with=FALSE])
  colnames(m) <- names(dt)[-1]
  rownames(m) <- cg_ids
  storage.mode(m) <- "numeric"
  m[is.na(m)] <- 0
  if (ncol(m) == 0 || nrow(m) == 0) next
  keep <- rowSums(m) > 0
  m <- m[keep, , drop=FALSE]
  cg_ids <- cg_ids[keep]
  if (nrow(m) == 0) next
  rs <- rowSums(m)
  if (nrow(m) > TOP_N) {
    o <- order(rs, decreasing=TRUE)[1:TOP_N]
    m <- m[o, , drop=FALSE]
    cg_ids <- cg_ids[o]
  }
  n_samples <- ncol(m)
  n_cgs <- nrow(m)
  lab <- vapply(cg_ids, wrap_label, character(1), width=35)
  cairo_pdf(file.path(out_dir, paste0("HEATMAP_big_", set_name, ".pdf")),width=35, height=38)

  # read cancer mapping
  pairs <- fread("./methylation/sample_pairs_files/methylation_pairs_filtered.tsv")
  meta <- unique(pairs[, .(patient, cancer)])

  # keep only patients present in matrix
  meta <- meta[patient %in% colnames(m)]

  # create named vector: old_name -> new_name
  new_names <- setNames(
    paste0(meta$cancer, "-", meta$patient),
    meta$patient
  )

  # rename columns
  colnames(m) <- new_names[colnames(m)]

  pheatmap(
    m,
    color = c("skyblue", "yellow"),
    breaks = c(-0.5, 0.5, 1.5),
    labels_row = lab,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 3,
    fontsize_col = 3,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    border_color = "grey60",
    legend_breaks = c(0, 1),
    legend_labels = c("absent (0)", "present (1)"),
    main = paste0(set_name,
                  " | patients=", n_samples,
                  " | CpGs in motifnpeaks=", n_cgs)
  )

  dev.off()
}
'

# merge the big matrix so that probes that affect the same motif are merged together and the motif associated with each probe is indicated in the rownames. 
Rscript -e '
library(data.table)

in_root  <- "./methylation/delta20_cg_lists/cg_presence_matrices"
map_root <- "./motifs/overlaps/intersected_motifs2mm_HM450"
out_root <- "./methylation/delta20_cg_lists/motif_presence_matrices"
dir.create(out_root, recursive=TRUE, showWarnings=FALSE)

setdirs <- list.dirs(in_root, recursive=FALSE, full.names=TRUE)

for (k in seq_along(setdirs)) {
  setdir <- setdirs[k]
  tf <- basename(setdir)

  big_file <- file.path(setdir, paste0("big_matrix_", tf, ".tsv.gz"))
  map_file <- file.path(map_root, paste0(tf, "_probeXmotif.tsv"))
  out_dir  <- file.path(out_root, tf)
  out_file <- file.path(out_dir, paste0("big_matrix_", tf, "_motif_count.tsv.gz"))

  cat("\n[", k, "/", length(setdirs), "] ", tf, "\n", sep="")

  if (!file.exists(big_file) || !file.exists(map_file)) {
    cat("  missing file, skipping\n")
    next
  }

  dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

  cat("  reading big matrix...\n")
  dt <- fread(cmd=paste("zcat", shQuote(big_file)))
  m <- as.matrix(dt[, -1, with=FALSE])
  rownames(m) <- dt[[1]]
  colnames(m) <- names(dt)[-1]
  storage.mode(m) <- "numeric"
  m[is.na(m)] <- 0
  cat("  probes in matrix:", nrow(m), "| patients:", ncol(m), "\n")

  cat("  reading probeXmotif...\n")
  mp <- fread(map_file, header=FALSE)
  mp <- unique(mp[, .(
    probe = V4,
    motif = paste(V5, V6, V7, V8, sep=":")
  )])
  cat("  probe-motif pairs:", nrow(mp), "| motifs:", uniqueN(mp$motif), "\n")

  mp <- mp[probe %in% rownames(m)]
  cat("  mapped probes kept:", uniqueN(mp$probe), "\n")
  if (nrow(mp) == 0) {
    cat("  no overlap, skipping\n")
    next
  }

  cat("  merging...\n")
  idx <- split(match(mp$probe, rownames(m)), mp$motif)

  merged <- t(sapply(idx, function(i) {
    if (length(i) == 1) m[i, ] else colSums(m[i, , drop=FALSE])
  }))

  merged <- as.data.table(merged, keep.rownames="motif")
  cat("  merged motifs:", nrow(merged), "\n")

  cat("  writing output...\n")
  fwrite(merged, out_file, sep="\t")
  cat("  saved:", out_file, "\n")
}
'
# plot the merged altered motifs
Rscript -e '
library(data.table)
library(pheatmap)

in_root  <- "./methylation/delta20_cg_lists/motif_presence_matrices"
map_root <- "./motifs/overlaps/intersected_motifs2mm_HM450"
out_dir  <- "./results/methylation/heatmap_motif_binary_gene_filtered"
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

pairs <- fread("./methylation/sample_pairs_files/methylation_pairs_filtered.tsv")
meta  <- unique(pairs[, .(patient, cancer)])

pg <- fread("./methylation/probe_gene_pairs_in_motifs.tsv")[, .(probe, gene)]
pg <- unique(pg)

for (d in list.dirs(in_root, recursive=FALSE, full.names=TRUE)) {

  tf <- basename(d)

  f_count <- file.path(d, paste0("big_matrix_", tf, "_motif_count.tsv.gz"))
  f_map   <- file.path(map_root, paste0(tf, "_probeXmotif.tsv"))

  if (!file.exists(f_count) || !file.exists(f_map)) next

  cat("Processing:", tf, "\n")

  dt <- fread(cmd=paste("zcat", shQuote(f_count)))
  m <- as.matrix(dt[, -1, with=FALSE])
  rownames(m) <- dt[[1]]
  colnames(m) <- names(dt)[-1]
  storage.mode(m) <- "numeric"
  m[is.na(m)] <- 0

  # binary
  m <- (m > 0) * 1

  # fix motif IDs
  rownames(m) <- sub(":[ACGT]+$", ":+", rownames(m))

  # remove empty motifs
  m <- m[rowSums(m) > 0, , drop=FALSE]
  if (nrow(m) == 0) next

  # motif -> gene
  mp <- fread(f_map, header=FALSE)[, .(
    probe = V4,
    motif = paste(V5, V6, V7, V10, sep=":")
  )]

  mg <- unique(merge(mp, pg, by="probe")[, .(motif, gene)])

  g <- mg$gene[match(rownames(m), mg$motif)]
  g[is.na(g)] <- "no_gene"

  row_lab <- paste0(g, " | ", rownames(m))

  # FILTER: keep motifs shared in >=10 patients
  keep <- rowSums(m) >= 10
  m <- m[keep, , drop=FALSE]
  row_lab <- row_lab[keep]
  if (nrow(m) == 0) next

  # rename columns
  meta_sub <- meta[match(colnames(m), patient)]
  colnames(m) <- paste0(meta_sub$cancer, "-", meta_sub$patient)

  # order rows
  o <- order(rowSums(m), decreasing=TRUE, rownames(m))
  m <- m[o, , drop=FALSE]
  row_lab <- row_lab[o]

  # group columns by cancer
  ordc <- order(meta_sub$cancer, meta_sub$patient)
  m <- m[, ordc, drop=FALSE]

  pdf(file.path(out_dir, paste0(tf, "_min10.pdf")), width=14, height=20)
  pheatmap(
    m,
    color = c("white", "yellow"),
    breaks = c(-0.5, 0.5, 1.5),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    labels_row = row_lab,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 4,
    fontsize_col = 4,
    border_color = "grey70",
    main = paste0(tf, " | motifs≥10 samples=", nrow(m))
  )
  dev.off()
}
'

# to build cancer specific matrix : an altered motif is considered present in a cancer if it is present in at least 50% of the patients of that cancer. 
Rscript -e '
library(data.table)
library(pheatmap)

# paths
in_root  <- "./methylation/delta20_cg_lists/motif_presence_matrices"
map_root <- "./motifs/overlaps/intersected_motifs2mm_HM450"
out_dir  <- "./results/methylation/heatmap_motif_by_cancer_50pct_gene"
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

# sample -> cancer
pairs <- fread("./methylation/sample_pairs_files/methylation_pairs_filtered.tsv")
meta  <- unique(pairs[, .(patient, cancer)])

# probe -> gene
pg <- fread("./methylation/probe_gene_pairs_in_motifs.tsv")[, .(probe, gene)]
pg <- unique(pg)

for (d in list.dirs(in_root, recursive=FALSE, full.names=TRUE)) {

  tf <- basename(d)

  f_count <- file.path(d, paste0("big_matrix_", tf, "_motif_count.tsv.gz"))
  f_map   <- file.path(map_root, paste0(tf, "_probeXmotif.tsv"))

  if (!file.exists(f_count) || !file.exists(f_map)) next

  cat("Processing:", tf, "\n")

  # read matrix
  dt <- fread(cmd=paste("zcat", shQuote(f_count)))
  m <- as.matrix(dt[, -1, with=FALSE])
  rownames(m) <- dt[[1]]
  colnames(m) <- names(dt)[-1]
  storage.mode(m) <- "numeric"
  m[is.na(m)] <- 0

  # binary
  m <- (m > 0) * 1

  # fix motif IDs
  rownames(m) <- sub(":[ACGT]+$", ":+", rownames(m))

  # remove empty motifs
  m <- m[rowSums(m) > 0, , drop=FALSE]
  if (nrow(m) == 0) next

  # motif -> gene
  mp <- fread(f_map, header=FALSE)[, .(
    probe = V4,
    motif = paste(V5, V6, V7, V10, sep=":")
  )]

  mg <- unique(merge(mp, pg, by="probe")[, .(motif, gene)])

  g <- mg$gene[match(rownames(m), mg$motif)]
  g[is.na(g)] <- "no_gene"

  row_lab <- paste0(g, " | ", rownames(m))

  # build cancer-level matrix
  meta_sub <- meta[match(colnames(m), patient)]
  cancers <- unique(meta_sub$cancer)

  mc <- sapply(cancers, function(ca) {
    cols <- which(meta_sub$cancer == ca)
    rowSums(m[, cols, drop=FALSE])
  })
  mc <- as.matrix(mc)
  colnames(mc) <- cancers

  # proportions
  n_patients <- table(meta_sub$cancer)
  prop <- sweep(mc, 2, n_patients, "/")

  # ≥50% rule
  present <- (prop >= 0.5) * 1

  # remove empty motifs
  present <- present[rowSums(present) > 0, , drop=FALSE]
  if (nrow(present) == 0) next

  # order rows
  o <- order(rowSums(present), decreasing=TRUE, rownames(present))
  present <- present[o, , drop=FALSE]
  row_lab <- row_lab[o]

  # order columns
  present <- present[, sort(colnames(present)), drop=FALSE]

  # ADD sample size directly to column names
  colnames(present) <- paste0(
    colnames(present),
    " (N=", n_patients[colnames(present)], ")"
  )

  # plot
  pdf(file.path(out_dir, paste0(tf, "_by_cancer_gene.pdf")), width=14, height=25)
  pheatmap(
    present,
    color = c("white", "red"),
    breaks = c(-0.5, 0.5, 1.5),
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    labels_row = row_lab,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 5,
    fontsize_col = 8,
    border_color = "grey70",
    main = paste0(tf, " | ≥50% patients per cancer")
  )
  dev.off()
}
'

#################################################################################################
# Redo altered motif with cg change more than 20% across all samples in 3d+4d without OV 
#################################################################################################









##################################################################################################################################################
# Barplot of the number of CpGs with |delta| >= 20% in each motif∩peaks set across all pairs, grouped by cancer type.
##################################################################################################################################################
Rscript -e '
library(data.table)

# ===== INPUT =====
big_mat_file <- "./methylation/delta20_cg_lists/cg_presence_matrices/NRF1_mm0to2_noCGmm/big_matrix_NRF1_mm0to2_noCGmm.tsv.gz"
pairs_file   <- "./methylation/sample_pairs_files/methylation_pairs_filtered.tsv"

out_file <- "./methylation/delta20_cg_lists/cg_presence_matrices/NRF1_mm0to2_noCGmm/CpG_counts_per_cancer_NRF1_mm0to2_noCGmm.tsv"

# ===== LOAD MATRIX =====
dt <- fread(cmd=paste("zcat", shQuote(big_mat_file)))
setnames(dt, 1, "probe")

# ===== LOAD META =====
pairs <- fread(pairs_file)
meta <- unique(pairs[, .(patient, cancer)])

# ===== MELT MATRIX =====
long <- melt(dt,
             id.vars="probe",
             variable.name="patient",
             value.name="value")

# keep only CpGs present (value == 1)
long <- long[value == 1]

# attach cancer type
long <- merge(long, meta, by="patient")

# ===== COUNT PER CpG PER CANCER =====
counts <- long[, .(n_samples = .N), by=.(probe, cancer)]

# optional: reshape to CpG x Cancer matrix
counts_mat <- dcast(counts, probe ~ cancer, value.var="n_samples", fill=0L)

fwrite(counts_mat, out_file, sep="\t")

cat("[DONE] -> ", out_file, "\n", sep="")
'

# plot barplot of the above counts (one per cancer type) to visualize the number of samples with |delta| >= 20% for each CpG in the list of cg with |delta| >= 20% and in the motif∩peaks set. Bars should be ordered by the most frequent CpGs across all cancers.

Rscript -e '
library(data.table)

# ===== INPUT =====
TF <- "NRF1"

in_file    <- "./methylation/delta20_cg_lists/cg_presence_matrices/NRF1_mm0to2_noCGmm/CpG_counts_per_cancer_NRF1_mm0to2_noCGmm.tsv"
pairs_file <- "./methylation/sample_pairs_files/methylation_pairs_filtered.tsv"
dist_file  <- "./methylation/dist_to_tss/HM450_TCGA-ACC-TCGA-OR-A5J2-01A_1_annotated_methylation_filtered_cpg_dist_tss.tsv"
out_pdf    <- "./results/methylation/heatmap_delta20_cg_motifnpeak/BARPLOT_per_CpG_NRF1_mm0to2_noCGmm.pdf"

PROX_BP <- 2000

dt <- fread(in_file)
cancers <- names(dt)[-1]

pairs <- fread(pairs_file)
tot_dt <- pairs[, .(total_samples = uniqueN(patient)), by=cancer]
tot_vec <- setNames(tot_dt$total_samples, tot_dt$cancer)

dist_dt <- fread(dist_file)
setnames(dist_dt, 1:2, c("probe","dist_to_tss"))
dist_dt[, dist_to_tss := as.numeric(dist_to_tss)]
dist_map <- setNames(dist_dt$dist_to_tss, dist_dt$probe)

dt[, total_presence := rowSums(dt[, -1, with=FALSE])]
dt <- dt[order(-total_presence)]

pdf(out_pdf, width=10, height=10)
par(mar = c(10, 4, 4, 2))

for(i in 1:nrow(dt)){
  cg <- dt$probe[i]
  counts <- as.numeric(dt[i, ..cancers])

  o <- order(counts, decreasing=TRUE)
  counts_sorted  <- counts[o]
  cancers_sorted <- cancers[o]
  ymax <- max(counts_sorted, na.rm=TRUE)

  bp <- barplot(
    counts_sorted,
    names.arg = cancers_sorted,
    las = 2,
    col = "steelblue",
    ylab = "Number of samples with DeltaBeta > 20%",
    ylim = c(0, ymax * 1.20 + 1)
  )

  totals_sorted <- tot_vec[cancers_sorted]
  pct <- ifelse(totals_sorted > 0, 100 * counts_sorted / totals_sorted, NA_real_)

  labs <- ifelse(is.na(pct),
                 paste0(counts_sorted, "/", totals_sorted),
                 paste0(counts_sorted, "/", totals_sorted, "\n", sprintf("%.1f%%", pct)))

  text(x = bp, y = counts_sorted + ymax*0.03 + 0.2, labels = labs, cex = 0.65)

  d <- dist_map[[cg]]
  pd <- ifelse(is.na(d), "NA", ifelse(abs(d) < PROX_BP, "proximal", "distal"))

  title(main = paste0(TF, " | ", cg,
                      " (", pd, ", distTSS=", ifelse(is.na(d),"NA",d), "bp) : samples per cancer"))

  grid()
}

dev.off()
cat("[DONE] -> ", out_pdf, "\n", sep="")
'

# count how many CpGs with |delta| >= 20% in each motif∩peaks set are proximal (< 2000bp) vs distal (>= 2000bp) to TSS. This is for the CpGs in the list of cg with |delta| >= 20% and in the motif∩peaks set. We will use the distance to TSS file we generated before and check the distance for each CpG in the list of cg with |delta| >= 20% and in the motif∩peaks set. We will then count how many of these CpGs are proximal vs distal to TSS.

Rscript -e '
library(data.table)

# files
counts_file <- "./methylation/delta20_cg_lists/cg_presence_matrices/BANP_mm0to2_noCGmm/CpG_counts_per_cancer_BANP_mm0to2_noCGmm.tsv"
dist_file   <- "./methylation/dist_to_tss/HM450_TCGA-ACC-TCGA-OR-A5J2-01A_1_annotated_methylation_filtered_cpg_dist_tss.tsv"

PROX_BP <- 2000

# load CpGs in your set
dt <- fread(counts_file)
cg_set <- dt$probe

# load distances
dist_dt <- fread(dist_file)
setnames(dist_dt, 1:2, c("probe","dist_to_tss"))
dist_dt[, dist_to_tss := as.numeric(dist_to_tss)]

# keep only CpGs in your set
dist_sub <- dist_dt[probe %in% cg_set]

# counts
n_prox <- sum(abs(dist_sub$dist_to_tss) < PROX_BP, na.rm=TRUE)
n_dist <- sum(abs(dist_sub$dist_to_tss) >= PROX_BP, na.rm=TRUE)

cat("Proximal CpGs (<", PROX_BP, "bp): ", n_prox, "\n", sep="")
cat("Distal CpGs (≥", PROX_BP, "bp): ", n_dist, "\n", sep="")
'
# For NRF1
# Proximal CpGs (<2000bp): 1264
# Distal CpGs (≥2000bp): 529
# For BANP
# Proximal CpGs (<2000bp): 378
# Distal CpGs (≥2000bp): 143

# to merge cg that affect the same motif instance : 
Rscript -e '
library(data.table)

tf <- "BANP_mm0to2_noCGmm"
prox_only <- TRUE
PROX_BP <- 2000

map <- fread(paste0("./motifs/overlaps/intersected_motifs2mm_HM450/", tf, "_probeXmotif.tsv"), header=FALSE)
cnt <- fread(paste0("./methylation/delta20_cg_lists/cg_presence_matrices/", tf, "/CpG_counts_per_cancer_", tf, ".tsv"))

# probe = V4 ; motif coords = V6,V7,V8
map <- map[, .(probe = V4, motif_id = paste(V6, V7, V8, sep=":"))]

if (prox_only) {
  dist <- fread("./methylation/dist_to_tss/HM450_TCGA-ACC-TCGA-OR-A5J2-01A_1_annotated_methylation_filtered_cpg_dist_tss.tsv")
  setnames(dist, 1:2, c("probe","dist_to_tss"))
  dist[, dist_to_tss := as.numeric(dist_to_tss)]
  cnt <- merge(cnt, dist[, .(probe, dist_to_tss)], by="probe")
  cnt <- cnt[abs(dist_to_tss) < PROX_BP]
  cnt[, dist_to_tss := NULL]
}

dt <- merge(map, cnt, by="probe")
res <- dt[, c(list(n_cpg = uniqueN(probe)),
              lapply(.SD, sum, na.rm=TRUE)),
          by = motif_id,
          .SDcols = setdiff(names(dt), c("probe","motif_id"))]

outdir <- paste0("./methylation/delta20_cg_lists/motif_presence_matrices/", tf)
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

suffix <- if (prox_only) "_proximal" else ""
fwrite(res, paste0(outdir, "/Motif_counts_per_cancer_", tf, suffix, ".tsv"), sep="\t")
'
# BANP motif instances with at least one CpG with |delta| >= 20% in proximal region :  334
# NRF1 motif instances with at least one CpG with |delta| >= 20% in proximal region :  1112

# To add the gene that is linked to that motif instance (based on the TSS that is closest to the motif instance) :
Rscript -e '
library(data.table)

tf <- "NRF1_mm0to2_noCGmm"

mot <- fread(paste0(
  "./methylation/delta20_cg_lists/motif_presence_matrices/", tf,
  "/Motif_counts_per_cancer_", tf, "_proximal.tsv"
))

g <- fread(paste0("./motifs/", tf, "_noXY_closest_genes.bed"), header=FALSE, skip=1)

g <- g[, .(
  motif_id = paste(V2, V3, V4, sep=":"),
  gene = V10,
  gene_id = V11,
  transcript_id = V12,
  gene_strand = V13,
  distance = V14
)]

g <- unique(g)

res <- merge(mot, g, by="motif_id", all.x=TRUE)

fwrite(
  res,
  paste0(
    "./methylation/delta20_cg_lists/motif_presence_matrices/", tf,
    "/Motif_counts_per_cancer_", tf, "_proximal_withClosestGene.tsv"
  ),
  sep="\t"
)
'
# merge the ones thayt affect the same gene :
Rscript -e '
library(data.table)

tf <- "NRF1_mm0to2_noCGmm"

dt <- fread(paste0(
  "./methylation/delta20_cg_lists/motif_presence_matrices/", tf,
  "/Motif_counts_per_cancer_", tf, "_proximal_withClosestGene.tsv"
))

# keep only motifs with a gene
dt <- dt[!is.na(gene) & gene != ""]

meta_cols <- c("motif_id","gene","gene_id","transcript_id","gene_strand","distance","n_cpg")
cancer_cols <- setdiff(names(dt), meta_cols)

res <- dt[, c(
  list(
    n_motifs = .N,
    total_cpg = sum(n_cpg, na.rm=TRUE),
    min_distance = min(distance, na.rm=TRUE)
  ),
  lapply(.SD, sum, na.rm=TRUE)
), by = .(gene, gene_id), .SDcols = cancer_cols]

fwrite(
  res,
  paste0(
    "./methylation/delta20_cg_lists/motif_presence_matrices/", tf,
    "/Gene_counts_per_cancer_", tf, "_proximal.tsv"
  ),
  sep = "\t"
)
'

# BANP has 324 genes that have at least one motif instance with a CpG with |delta| >= 20% in proximal region
# NRF1 has 1022 genes that have at least one motif instance with a CpG with |delta| >= 20% in proximal region





###########################################################################################
# SMOOTH SCATTER PLOT FOR REPLICATES for 6 mer motifs
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

out_pdf <- "./results/methylation/smoothScatter_ALL_vs_NRF1_vs_BANP_per_pair_replicates_cancer.pdf"
thr <- 20

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
# HELPERS
# =========================
pct_fmt <- function(x, n) {
  if (n == 0) return("NA")
  sprintf("%.1f%%", 100 * x / n)
}

add_class_legend <- function(delta_vec, thr, title = "All CpGs") {
  n <- length(delta_vec)
  hypo <- sum(delta_vec <= -thr, na.rm = TRUE)
  hyper <- sum(delta_vec >=  thr, na.rm = TRUE)
  unch <- sum(abs(delta_vec) < thr, na.rm = TRUE)

  lab <- paste0(
    title, "\n",
    "Hypo (Delta≤-", thr, "): ", pct_fmt(hypo, n), "\n",
    "Unchanged (|Delta|<", thr, "): ", pct_fmt(unch, n), "\n",
    "Hyper (Delta≥+", thr, "): ", pct_fmt(hyper, n)
  )

  usr <- par("usr")
  x <- usr[1] + 0.02 * (usr[2] - usr[1])
  y <- usr[4] - 0.02 * (usr[4] - usr[3])

  text(x, y, labels = lab, adj = c(0, 1), cex = 0.85)
}

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
  rep1_df <- data.frame(probe = rep1[[4]], meth_rep1 = rep1[[5]])
  rep2_df <- data.frame(probe = rep2[[4]], meth_rep2 = rep2[[5]])

  merged <- merge(rep1_df, rep2_df, by = "probe")
  if (nrow(merged) == 0) next

  # -------------------------
  # DROP CpGs WITH ANY NA
  # -------------------------
  merged <- merged[!is.na(merged$meth_rep1) & !is.na(merged$meth_rep2), ]

  if (nrow(merged) == 0) {
    plot.new()
    title(main = paste0(
      pairs$cancer[i], " | ", pairs$patient[i],
      "\nNo CpGs with beta in both replicates"
    ))
    next
  }

  # -------------------------
  # DELTA (rep1 - rep2)
  # -------------------------
  delta <- merged$meth_rep1 - merged$meth_rep2

  # =========================
  # PLOT 1 — ALL CpGs + %
  # =========================
  smoothScatter(
    merged$meth_rep2,
    merged$meth_rep1,
    xlab = "Replicate 2 methylation (%)",
    ylab = "Replicate 1 methylation (%)",
    main = paste0(pairs$cancer[i], " | ", pairs$patient[i], "\nAll CpGs")
  )
  abline(0, 1, col = "black", lwd = 2)
  abline(a = -thr, b = 1, col = "blue", lwd = 2)
  abline(a =  thr, b = 1, col = "blue", lwd = 2)
  add_class_legend(delta, thr, title = "All CpGs")

  # =========================
  # PLOT 2 — NRF1 motif CpGs + % within NRF1 CpGs
  # =========================
  smoothScatter(
    merged$meth_rep2,
    merged$meth_rep1,
    xlab = "Replicate 2 methylation (%)",
    ylab = "Replicate 1 methylation (%)",
    main = "NRF1 motif CpGs in red (|Deltabeta| ≥ 20%)"
  )
  in_NRF1 <- merged$probe %in% motif_cpgs$NRF1
  is_NRF1 <- in_NRF1 & (abs(delta) >= thr)
  points(
    merged$meth_rep2[is_NRF1],
    merged$meth_rep1[is_NRF1],
    pch = 16, cex = 0.5, col = rgb(1, 0, 0, 0.6)
  )
  abline(0, 1, col = "black", lwd = 2)
  abline(a = -thr, b = 1, col = "blue", lwd = 2)
  abline(a =  thr, b = 1, col = "blue", lwd = 2)
  add_class_legend(delta[in_NRF1], thr, title = paste0("NRF1 CpGs"))

  # =========================
  # PLOT 3 — BANP motif CpGs + % within BANP CpGs
  # =========================
  smoothScatter(
    merged$meth_rep2,
    merged$meth_rep1,
    xlab = "Replicate 2 methylation (%)",
    ylab = "Replicate 1 methylation (%)",
    main = "BANP motif CpGs in red (|Deltabeta| ≥ 20%)"
  )
  in_BANP <- merged$probe %in% motif_cpgs$BANP
  is_BANP <- in_BANP & (abs(delta) >= thr)
  points(
    merged$meth_rep2[is_BANP],
    merged$meth_rep1[is_BANP],
    pch = 16, cex = 0.5, col = rgb(1, 0, 0, 0.6)
  )
  abline(0, 1, col = "black", lwd = 2)
  abline(a = -thr, b = 1, col = "blue", lwd = 2)
  abline(a =  thr, b = 1, col = "blue", lwd = 2)
  add_class_legend(delta[in_BANP], thr, title = paste0("BANP CpGs"))
}

dev.off()
cat("DONE →", out_pdf, "\n")
'

###########################################################################################
# SMOOTH SCATTER PLOT FOR REPLICATES for 2mm motifs
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

# These must be PROBE lines (cg IDs in column 4), i.e. probes as -a in bedtools intersect
motif_files <- list(
  NRF1_mm0to2        = "./methylation/overlaps/intersected_motifs2mm_HM450/NRF1_mm0to2_intersected_methylation.bed",
  BANP_mm0to2        = "./methylation/overlaps/intersected_motifs2mm_HM450/BANP_mm0to2_intersected_methylation.bed",
  NRF1_mm0to2_noCGmm = "./methylation/overlaps/intersected_motifs2mm_HM450/NRF1_mm0to2_noCGmm_intersected_methylation.bed",
  BANP_mm0to2_noCGmm = "./methylation/overlaps/intersected_motifs2mm_HM450/BANP_mm0to2_noCGmm_intersected_methylation.bed"
)

out_pdf <- "./results/methylation/smoothscatter/smoothScatter_replicates_mm0to2_vs_noCGmm_mm0to2.pdf"
thr <- 20

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
  if(!file.exists(f)) stop(paste("Missing motif probe file:", f))
  unique(trimws(as.character(fread(f, header = FALSE)[[4]])))
})

# =========================
# HELPERS
# =========================
pct_fmt <- function(x, n) {
  if (n == 0) return("0.0%")
  sprintf("%.1f%%", 100 * x / n)
}

add_class_legend <- function(delta_vec, thr, title = "All CpGs") {
  delta_vec <- delta_vec[is.finite(delta_vec)]
  n <- length(delta_vec)

  hypo  <- sum(delta_vec <= -thr)
  hyper <- sum(delta_vec >=  thr)
  unch  <- sum(abs(delta_vec) < thr)

  lab <- paste0(
    title, " (n=", n, ")\n",
    "Hypo (Delta≤-", thr, "): ", pct_fmt(hypo, n), "\n",
    "Unchanged (|Delta|<", thr, "): ", pct_fmt(unch, n), "\n",
    "Hyper (Delta≥+", thr, "): ", pct_fmt(hyper, n)
  )

  usr <- par("usr")
  x <- usr[1] + 0.02 * (usr[2] - usr[1])
  y <- usr[4] - 0.02 * (usr[4] - usr[3])
  text(x, y, labels = lab, adj = c(0, 1), cex = 0.85)
}

plot_all <- function(merged, delta, thr, main_title){
  smoothScatter(
    merged$meth_rep2, merged$meth_rep1,
    xlab = "Replicate 2 methylation (%)",
    ylab = "Replicate 1 methylation (%)",
    main = main_title
  )
  abline(0, 1, col = "black", lwd = 2)
  abline(a = -thr, b = 1, col = "blue", lwd = 2)
  abline(a =  thr, b = 1, col = "blue", lwd = 2)
  add_class_legend(delta, thr, title = "All CpGs")
}

plot_motif_red <- function(merged, delta, thr, in_motif, label_title){
  smoothScatter(
    merged$meth_rep2, merged$meth_rep1,
    xlab = "Replicate 2 methylation (%)",
    ylab = "Replicate 1 methylation (%)",
    main = paste0(label_title, " in red (|DeltaBeta| ≥ ", thr, "%)")
  )
  is_red <- in_motif & (abs(delta) >= thr)
  points(
    merged$meth_rep2[is_red],
    merged$meth_rep1[is_red],
    pch = 16, cex = 0.5, col = rgb(1, 0, 0, 0.6)
  )
  abline(0, 1, col = "black", lwd = 2)
  abline(a = -thr, b = 1, col = "blue", lwd = 2)
  abline(a =  thr, b = 1, col = "blue", lwd = 2)
  add_class_legend(delta[in_motif], thr, title = label_title)
}

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

  rep1_path <- pairs$replicate_1[i]
  rep2_path <- pairs$replicate_2[i]

  if(!file.exists(rep1_path) || !file.exists(rep2_path)){
    plot.new()
    title(main = paste0(pairs$cancer[i], " | ", pairs$patient[i], "\nMissing replicate file(s)"))
    next
  }

  rep1 <- fread(cmd = paste("zcat", shQuote(rep1_path)), header = FALSE)
  rep2 <- fread(cmd = paste("zcat", shQuote(rep2_path)), header = FALSE)

  rep1_df <- data.frame(probe = trimws(as.character(rep1[[4]])), meth_rep1 = rep1[[5]])
  rep2_df <- data.frame(probe = trimws(as.character(rep2[[4]])), meth_rep2 = rep2[[5]])

  merged <- merge(rep1_df, rep2_df, by = "probe")
  if (nrow(merged) == 0) next

  merged <- merged[!is.na(merged$meth_rep1) & !is.na(merged$meth_rep2), ]
  if (nrow(merged) == 0) {
    plot.new()
    title(main = paste0(pairs$cancer[i], " | ", pairs$patient[i],
                        "\nNo CpGs with beta in both replicates"))
    next
  }

  delta <- merged$meth_rep1 - merged$meth_rep2

  # =========================
  # PAGE 1 — mm0to2
  # =========================
  plot_all(merged, delta, thr,
           paste0(pairs$cancer[i], " | ", pairs$patient[i], "\nAll CpGs (mm0to2 page)"))

  in_NRF1 <- merged$probe %in% motif_cpgs$NRF1_mm0to2
  plot_motif_red(merged, delta, thr, in_NRF1, "NRF1_mm0to2 CpGs")

  in_BANP <- merged$probe %in% motif_cpgs$BANP_mm0to2
  plot_motif_red(merged, delta, thr, in_BANP, "BANP_mm0to2 CpGs")

  # =========================
  # PAGE 2 — mm0to2_noCGmm
  # =========================
  plot_all(merged, delta, thr,
           paste0(pairs$cancer[i], " | ", pairs$patient[i], "\nAll CpGs (noCGmm page)"))

  in_NRF1_nc <- merged$probe %in% motif_cpgs$NRF1_mm0to2_noCGmm
  plot_motif_red(merged, delta, thr, in_NRF1_nc, "NRF1_mm0to2_noCGmm CpGs")

  in_BANP_nc <- merged$probe %in% motif_cpgs$BANP_mm0to2_noCGmm
  plot_motif_red(merged, delta, thr, in_BANP_nc, "BANP_mm0to2_noCGmm CpGs")
}

dev.off()
cat("DONE -> ", out_pdf, "\n", sep="")
'

###############################################################################
# SECTION 9 — Distance OF CpG PROBES IN 2mm MOTIFS TO TSS
###############################################################################
# Look into the distance of CpG probes in motifs to othe TSS to see if they are proximal or distal
# for each tsv file in ./methylation/dist_to_tss/ filter only those in motif of NRF1 and make a histogram of distances to TSS for each sample
# NRF1_mm0to2 
mkdir -p ./results/methylation/methylation_motifs_TSS/NRF1_mm0to2
mkdir -p ./methylation/dist_to_tss_motif/NRF1_mm0to2

awk 'BEGIN{FS=OFS="\t"} {print $4}' \
  ./methylation/overlaps/intersected_motifs2mm_HM450/NRF1_mm0to2_intersected_methylation.bed \
| grep -E '^cg' | sort -u \
> ./methylation/overlaps/intersected_motifs2mm_HM450/NRF1_mm0to2_probes.txt

for f in ./methylation/dist_to_tss/*_cpg_dist_tss.tsv; do
  base=$(basename "$f")
  out="./methylation/dist_to_tss_motif/NRF1_mm0to2/${base/_cpg_dist_tss.tsv/_NRF1_mm0to2_dist_tss.txt}"

  awk 'NR==FNR{a[$1]=1; next} ($1 in a){print $2}' \
    ./methylation/overlaps/intersected_motifs2mm_HM450/NRF1_mm0to2_probes.txt \
    "$f" > "$out"
done
# BANP_mm0to2 
mkdir -p ./results/methylation/methylation_motifs_TSS/BANP_mm0to2
mkdir -p ./methylation/dist_to_tss_motif/BANP_mm0to2

awk 'BEGIN{FS=OFS="\t"} {print $4}' \
  ./methylation/overlaps/intersected_motifs2mm_HM450/BANP_mm0to2_intersected_methylation.bed \
| grep -E '^cg' | sort -u \
> ./methylation/overlaps/intersected_motifs2mm_HM450/BANP_mm0to2_probes.txt

for f in ./methylation/dist_to_tss/*_cpg_dist_tss.tsv; do
  base=$(basename "$f")
  out="./methylation/dist_to_tss_motif/BANP_mm0to2/${base/_cpg_dist_tss.tsv/_BANP_mm0to2_dist_tss.txt}"

  awk 'NR==FNR{a[$1]=1; next} ($1 in a){print $2}' \
    ./methylation/overlaps/intersected_motifs2mm_HM450/BANP_mm0to2_probes.txt \
    "$f" > "$out"
done

# NRF1_mm0to2_noCGmm 
mkdir -p ./results/methylation/methylation_motifs_TSS/NRF1_mm0to2_noCGmm
mkdir -p ./methylation/dist_to_tss_motif/NRF1_mm0to2_noCGmm

awk 'BEGIN{FS=OFS="\t"} {print $4}' \
  ./methylation/overlaps/intersected_motifs2mm_HM450/NRF1_mm0to2_noCGmm_intersected_methylation.bed \
| grep -E '^cg' | sort -u \
> ./methylation/overlaps/intersected_motifs2mm_HM450/NRF1_mm0to2_noCGmm_probes.txt

for f in ./methylation/dist_to_tss/*_cpg_dist_tss.tsv; do
  base=$(basename "$f")
  out="./methylation/dist_to_tss_motif/NRF1_mm0to2_noCGmm/${base/_cpg_dist_tss.tsv/_NRF1_mm0to2_noCGmm_dist_tss.txt}"

  awk 'NR==FNR{a[$1]=1; next} ($1 in a){print $2}' \
    ./methylation/overlaps/intersected_motifs2mm_HM450/NRF1_mm0to2_noCGmm_probes.txt \
    "$f" > "$out"
done
# BANP_mm0to2_noCGmm 
mkdir -p ./results/methylation/methylation_motifs_TSS/BANP_mm0to2_noCGmm
mkdir -p ./methylation/dist_to_tss_motif/BANP_mm0to2_noCGmm

awk 'BEGIN{FS=OFS="\t"} {print $4}' \
  ./methylation/overlaps/intersected_motifs2mm_HM450/BANP_mm0to2_noCGmm_intersected_methylation.bed \
| grep -E '^cg' | sort -u \
> ./methylation/overlaps/intersected_motifs2mm_HM450/BANP_mm0to2_noCGmm_probes.txt

for f in ./methylation/dist_to_tss/*_cpg_dist_tss.tsv; do
  base=$(basename "$f")
  out="./methylation/dist_to_tss_motif/BANP_mm0to2_noCGmm/${base/_cpg_dist_tss.tsv/_BANP_mm0to2_noCGmm_dist_tss.txt}"

  awk 'NR==FNR{a[$1]=1; next} ($1 in a){print $2}' \
    ./methylation/overlaps/intersected_motifs2mm_HM450/BANP_mm0to2_noCGmm_probes.txt \
    "$f" > "$out"
done

# Plot histogram of distances to TSS for NRF1 and BANP motifs per cancer type
Rscript -e '
library(data.table)

out_pdf <- "./results/methylation/methylation_motifs_TSS/dist_tss_methylation_NRF1_BANP_mm0to2_byCancer.pdf"
dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)

pdf(out_pdf, width = 11.7, height = 8.3)
par(bg = "white")

cases <- c(
  "NRF1_mm0to2",
  "BANP_mm0to2",
  "NRF1_mm0to2_noCGmm",
  "BANP_mm0to2_noCGmm"
)

for (CASE in cases) {

  dist_dir <- file.path("./methylation/dist_to_tss_motif", CASE)
  out_tsv  <- file.path("./results/methylation/methylation_motifs_TSS", CASE,
                        paste0("summary_", CASE, "_byCancer.tsv"))

  files <- list.files(dist_dir, pattern = paste0("_", CASE, "_dist_tss.txt$"),
                      full.names = TRUE)

  if (!length(files)) {
    plot.new()
    title(main = paste0(CASE, ": no distance files found"))
    next
  }

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

    d_list <- lapply(fc, function(f) suppressWarnings(as.numeric(readLines(f))))
    d <- unlist(d_list, use.names = FALSE)
    d <- d[is.finite(d)]

    n_samples <- length(fc)
    n_vals <- length(d)

    if (n_vals == 0) {
      plot.new()
      title(main = paste0("TCGA-", c, " | ", CASE, ": none found"))
      res[[c]] <- data.frame(cancer=c, n_samples=n_samples, n_values=0,
                             pct_prox=NA, pct_dist=NA, median=NA)
      next
    }

    proximal <- sum(d < 2000)
    distal   <- sum(d >= 2000)
    pct_prox <- round(100 * proximal / n_vals, 1)
    pct_dist <- round(100 * distal   / n_vals, 1)

    hist(log10(d + 1),
         xlim = c(0, 8),
         xlab = "Distance to TSS (log10)",
         main = paste0("TCGA-", c, " | ", CASE,
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


###########################################
# Find altered motifs in the 3d+4d samples this outputs a matrix with motif by sample table for all altered motifs (altered cpg change more than 20% in methylation) in the tumor samples from the 3d+4d samples. 
###########################################
Rscript - <<'EOF'
suppressPackageStartupMessages(library(data.table))

samples_file <- "./results/multi_omics/samples_3d+4d_noOV.tsv"
motif_file   <- "./motifs/overlaps/intersected_motifs2mm_HM450/BANP_mm0to2_noCGmm_probeXmotif.tsv"
out_file     <- "./results/methylation/Cpg_altered/BANP_motif_sample_matrix_noOV.tsv.gz"

delta_thr <- 20
min_normals_same <- 3
min_normals_related <- 3

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

fallback_ref <- c(
  ACC=NA, BLCA=NA, BRCA=NA, CESC="UCEC", CHOL="LIHC", COAD="READ", DLBC=NA,
  ESCA="STAD", GBM="LGG", HNSC="ESCA", KICH="KIRP", KIRC="KIRP", KIRP="KIRC",
  LAML=NA, LGG="GBM", LIHC="CHOL", LUAD="LUSC", LUSC="LUAD", MESO=NA, PAAD=NA,
  PCPG="ACC", PRAD=NA, READ="COAD", SARC=NA, SKCM=NA, STAD="ESCA", TGCT=NA,
  THCA=NA, THYM=NA, UCEC="CESC", UCS="UCEC", UVM="SKCM"
)

read_meth <- function(f) {
  x <- fread(cmd = paste("zcat", shQuote(f)), header = FALSE,
             na.strings = c("NA","","NaN"), showProgress = FALSE)
  x[, .(probe = as.character(V4), beta = as.numeric(V5))]
}

extract_sample_type <- function(x) {
  x <- basename(x)
  if (grepl("-01A([_.-]|$)", x)) return("tumor")
  if (grepl("-11A([_.-]|$)", x)) return("normal")
  NA_character_
}
extract_patient_id <- function(x) sub("^HM450_TCGA-[^-]+-(TCGA-[^-]+-[^-]+)-.*$", "\\1", basename(x))
extract_cancer     <- function(x) sub("^HM450_TCGA-([^-]+)-.*$", "\\1", basename(x))

samples <- fread(samples_file, header = FALSE)
setnames(samples, c("sample","cancer","SNV","METH","RNA","ATAC"))
samples <- samples[METH == 1, .(sample, cancer)]

motif_raw <- fread(motif_file, header = FALSE, showProgress = FALSE)
motif_map <- unique(motif_raw[, .(
  probe    = as.character(V4),
  motif_id = paste(V5, V6, V7, V10, sep=":")
)])

all_files <- list.files(
  "./methylation/filtered_methylation",
  pattern = "annotated_methylation_filtered[.]bed[.]gz$",
  full.names = TRUE
)

file_index <- data.table(
  meth_file   = all_files,
  patient_id  = vapply(all_files, extract_patient_id, character(1)),
  cancer      = vapply(all_files, extract_cancer, character(1)),
  sample_type = vapply(all_files, extract_sample_type, character(1))
)[!is.na(sample_type)]

tumor_index  <- file_index[sample_type == "tumor"]
normal_index <- file_index[sample_type == "normal"]

tumor_samples <- merge(
  samples, tumor_index,
  by.x = c("sample","cancer"),
  by.y = c("patient_id","cancer")
)

if (nrow(tumor_samples) == 0) stop("No tumor methylation files matched to cohort.")
if (nrow(normal_index) == 0) stop("No normal methylation files found.")

cat("Tumor files matched :", nrow(tumor_samples), "\n")
cat("Normal files found  :", nrow(normal_index), "\n")
cat("Unique motif probes :", uniqueN(motif_map$probe), "\n")
cat("Unique motifs       :", uniqueN(motif_map$motif_id), "\n")

# =========================
# READ NORMALS ONCE
# =========================
t_ref1 <- Sys.time()
cat("Reading normal methylation files once...\n")

norm_list <- vector("list", nrow(normal_index))

for (i in seq_len(nrow(normal_index))) {
  if (i %% 10 == 0 || i == 1 || i == nrow(normal_index)) {
    cat("[", i, "/", nrow(normal_index), "] normals\n", sep = "")
  }
  f  <- normal_index$meth_file[i]
  cc <- normal_index$cancer[i]

  dt <- read_meth(f)
  dt <- motif_map[dt, on = "probe", nomatch = 0]
  dt <- dt[!is.na(beta)]
  if (nrow(dt) == 0) next
  dt[, cancer := cc]
  norm_list[[i]] <- dt[, .(cancer, probe, beta)]
}

norm_dt <- rbindlist(norm_list, fill = TRUE)
if (nrow(norm_dt) == 0) stop("No motif-overlapping normal beta values found.")

cat("Building cancer-specific and pan-cancer references...\n")

ref_by_cancer <- norm_dt[, .(ref_beta = median(beta, na.rm = TRUE)), by = .(cancer, probe)]
ref_pan       <- norm_dt[, .(ref_beta = median(beta, na.rm = TRUE)), by = probe]

t_ref2 <- Sys.time()
cat("Reference build time:", round(as.numeric(difftime(t_ref2, t_ref1, units = "mins")), 2), "minutes\n")

# counts of available normal files by cancer
normal_counts <- normal_index[, .N, by = cancer]
setnames(normal_counts, "N", "n_normals")

# =========================
# PROCESS TUMORS
# =========================
cat("Processing tumor samples...\n")
t_tum1 <- Sys.time()

res <- vector("list", nrow(tumor_samples))

for (i in seq_len(nrow(tumor_samples))) {
  s <- tumor_samples[i]
  sample_id <- s$sample
  ccur <- s$cancer

  if (i %% 25 == 0 || i == 1 || i == nrow(tumor_samples)) {
    cat("[", i, "/", nrow(tumor_samples), "] ", sample_id, " ", ccur, "\n", sep = "")
  }

  dt <- read_meth(s$meth_file)
  dt <- motif_map[dt, on = "probe", nomatch = 0]
  dt <- dt[!is.na(beta)]
  if (nrow(dt) == 0) next

  n_same <- normal_counts[cancer == ccur, n_normals]
  if (length(n_same) == 0) n_same <- 0L

  related <- if (ccur %in% names(fallback_ref)) fallback_ref[[ccur]] else NA_character_
  n_related <- if (!is.na(related)) normal_counts[cancer == related, n_normals] else 0L
  if (length(n_related) == 0) n_related <- 0L

  if (n_same >= min_normals_same && any(ref_by_cancer$cancer == ccur)) {
    ref_dt <- ref_by_cancer[cancer == ccur, .(probe, ref_beta)]
  } else if (!is.na(related) && n_related >= min_normals_related &&
             any(ref_by_cancer$cancer == related)) {
    ref_dt <- ref_by_cancer[cancer == related, .(probe, ref_beta)]
  } else {
    ref_dt <- ref_pan
  }

  dt <- ref_dt[dt, on = "probe", nomatch = 0]
  if (nrow(dt) == 0) next

  dt <- dt[abs(beta - ref_beta) >= delta_thr]
  if (nrow(dt) == 0) next

  res[[i]] <- unique(dt[, .(motif_id, sample_id)])[, motif_altered := 1L]
}

motif_sample_dt <- rbindlist(res, fill = TRUE)
if (nrow(motif_sample_dt) == 0) stop("No altered motifs detected.")

motif_matrix <- dcast(
  motif_sample_dt,
  motif_id ~ sample_id,
  value.var = "motif_altered",
  fill = 0
)

fwrite(motif_matrix, out_file, sep = "\t")
cat("[DONE] Wrote:", out_file, "\n")

t_tum2 <- Sys.time()
cat("Tumor processing time:", round(as.numeric(difftime(t_tum2, t_tum1, units = "mins")), 2), "minutes\n")
cat("Total time:", round(as.numeric(difftime(t_tum2, t_ref1, units = "mins")), 2), "minutes\n")
EOF

# plot all motifs 
Rscript - <<'EOF'
suppressPackageStartupMessages({
  library(data.table)
  library(pheatmap)
})

mat_file     <- "./results/methylation/Cpg_altered/NRF1_motif_sample_matrix_noOV.tsv.gz"
samples_file <- "./results/multi_omics/samples_3d+4d_noOV.tsv"
col_file     <- "./results/multi_omics/cancer_color_order_with_defined_colours.tsv"
out_pdf      <- "./results/methylation/Cpg_altered/NRF1_motif_sample_matrix_noOV_heatmap_allMotifs.pdf"

dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)

# =========================
# LOAD MATRIX
# =========================
mat_dt <- fread(mat_file)
if (ncol(mat_dt) < 2) stop("Matrix file has no sample columns: ", mat_file)

motif_ids <- mat_dt$motif_id
mat <- as.matrix(mat_dt[, -1, with = FALSE])
rownames(mat) <- motif_ids
mode(mat) <- "numeric"

cat("Number of motifs plotted :", nrow(mat), "\n")
cat("Number of samples plotted:", ncol(mat), "\n")

# =========================
# LOAD SAMPLE ANNOTATION
# =========================
samples <- fread(samples_file, header = FALSE)
setnames(samples, c("sample", "cancer", "SNV", "METH", "RNA", "ATAC"))
samples <- unique(samples[METH == 1, .(sample, cancer)])

ann <- samples[sample %in% colnames(mat)]
if (nrow(ann) == 0) stop("No matrix columns matched sample annotation.")

ann[, Cancer := as.character(cancer)]

# keep only annotated samples and order matrix accordingly
sample_order <- ann$sample[ann$sample %in% colnames(mat)]
mat <- mat[, sample_order, drop = FALSE]

# =========================
# LOAD CANCER COLORS / ORDER
# =========================
cols <- fread(col_file, header = FALSE)[, 1:2]
setnames(cols, c("Cancer", "Color"))
cols[, Cancer := gsub("\r", "", trimws(as.character(Cancer)))]
cols[, Color  := gsub("\r", "", trimws(as.character(Color)))]
cols[, Cancer := sub("^TCGA-", "", Cancer)]

cancer_levels <- unique(cols$Cancer)
cancer_cols_all <- setNames(cols$Color, cols$Cancer)

missing_cancers <- setdiff(unique(ann$Cancer), names(cancer_cols_all))
if (length(missing_cancers) > 0) {
  extra_cols <- grDevices::rainbow(length(missing_cancers))
  names(extra_cols) <- missing_cancers
  cancer_cols_all <- c(cancer_cols_all, extra_cols)
  cancer_levels <- c(cancer_levels, missing_cancers)
}

ann[, Cancer := factor(Cancer, levels = cancer_levels)]
setorder(ann, Cancer, sample)

sample_order <- ann$sample
mat <- mat[, sample_order, drop = FALSE]

ann_df <- data.frame(
  Cancer = ann$Cancer,
  row.names = ann$sample,
  stringsAsFactors = FALSE
)

# =========================
# COLUMN GROUPS / LABELS
# =========================
group_sizes <- ann[, .N, by = Cancer][order(match(Cancer, cancer_levels))]$N
gaps_col <- cumsum(group_sizes)
gaps_col <- gaps_col[gaps_col < ncol(mat)]

labels_col <- rep("", ncol(mat))
group_dt <- ann[, .(
  start = .I[1],
  end   = .I[.N]
), by = Cancer]
group_dt[, mid := floor((start + end) / 2)]
labels_col[group_dt$mid] <- as.character(group_dt$Cancer)

ann_colors <- list(Cancer = cancer_cols_all)

# =========================
# PLOT
# =========================
pdf(out_pdf, width = 20, height = 30)

pheatmap(
  mat,
  color = c("white", "red"),
  breaks = c(-0.001, 0.5, 1.001),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  gaps_col = gaps_col,
  show_rownames = FALSE,
  show_colnames = TRUE,
  labels_col = labels_col,
  angle_col = 90,
  fontsize_col = 8,
  annotation_col = ann_df,
  annotation_colors = ann_colors,
  main = "NRF1 altered motifs across 3d+4d noOV tumor samples\nBinary matrix: 1 = altered (|Delta beta| >= 20), 0 = not altered",
  border_color = NA
)

dev.off()

cat("[DONE] Wrote:", out_pdf, "\n")
EOF

# the same heatmap but only the first 50 altered motifs to have a more readable heatmap for presentations 
Rscript - <<'EOF'
suppressPackageStartupMessages({
  library(data.table)
  library(pheatmap)
})

mat_file     <- "./results/methylation/Cpg_altered/BANP_motif_sample_matrix_noOV.tsv.gz"
samples_file <- "./results/multi_omics/samples_3d+4d_noOV.tsv"
col_file     <- "./results/multi_omics/cancer_color_order_with_defined_colours.tsv"
out_pdf      <- "./results/methylation/Cpg_altered/BANP_motif_sample_matrix_noOV_heatmap_first50.pdf"

dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)

# =========================
# LOAD MATRIX
# =========================
mat_dt <- fread(mat_file)
if (ncol(mat_dt) < 2) stop("Matrix file has no sample columns: ", mat_file)

mat_dt <- mat_dt[1:min(50, .N)]   # keep only first 50 motifs

motif_ids <- mat_dt$motif_id
mat <- as.matrix(mat_dt[, -1, with = FALSE])
rownames(mat) <- motif_ids
mode(mat) <- "numeric"

cat("Number of motifs plotted :", nrow(mat), "\n")
cat("Number of samples plotted:", ncol(mat), "\n")

# =========================
# LOAD SAMPLE ANNOTATION
# =========================
samples <- fread(samples_file, header = FALSE)
setnames(samples, c("sample", "cancer", "SNV", "METH", "RNA", "ATAC"))
samples <- unique(samples[METH == 1, .(sample, cancer)])

ann <- samples[sample %in% colnames(mat)]
if (nrow(ann) == 0) stop("No matrix columns matched sample annotation.")

ann[, Cancer := as.character(cancer)]

sample_order <- ann$sample[ann$sample %in% colnames(mat)]
mat <- mat[, sample_order, drop = FALSE]

# =========================
# LOAD CANCER COLORS / ORDER
# =========================
cols <- fread(col_file, header = FALSE)[, 1:2]
setnames(cols, c("Cancer", "Color"))
cols[, Cancer := gsub("\r", "", trimws(as.character(Cancer)))]
cols[, Color  := gsub("\r", "", trimws(as.character(Color)))]
cols[, Cancer := sub("^TCGA-", "", Cancer)]

cancer_levels <- unique(cols$Cancer)
cancer_cols_all <- setNames(cols$Color, cols$Cancer)

missing_cancers <- setdiff(unique(ann$Cancer), names(cancer_cols_all))
if (length(missing_cancers) > 0) {
  extra_cols <- grDevices::rainbow(length(missing_cancers))
  names(extra_cols) <- missing_cancers
  cancer_cols_all <- c(cancer_cols_all, extra_cols)
  cancer_levels <- c(cancer_levels, missing_cancers)
}

ann[, Cancer := factor(Cancer, levels = cancer_levels)]
setorder(ann, Cancer, sample)

sample_order <- ann$sample
mat <- mat[, sample_order, drop = FALSE]

ann_df <- data.frame(
  Cancer = ann$Cancer,
  row.names = ann$sample,
  stringsAsFactors = FALSE
)

group_sizes <- ann[, .N, by = Cancer][order(match(Cancer, cancer_levels))]$N
gaps_col <- cumsum(group_sizes)
gaps_col <- gaps_col[gaps_col < ncol(mat)]

labels_col <- rep("", ncol(mat))
group_dt <- ann[, .(start = .I[1], end = .I[.N]), by = Cancer]
group_dt[, mid := floor((start + end) / 2)]
labels_col[group_dt$mid] <- as.character(group_dt$Cancer)

ann_colors <- list(Cancer = cancer_cols_all)

# =========================
# PLOT
# =========================
pdf(out_pdf, width = 20, height = 12)

pheatmap(
  mat,
  color = c("white", "red"),
  breaks = c(-0.001, 0.5, 1.001),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  gaps_col = gaps_col,
  show_rownames = TRUE,
  show_colnames = TRUE,
  labels_col = labels_col,
  angle_col = 90,
  fontsize_row = 6,
  fontsize_col = 8,
  annotation_col = ann_df,
  annotation_colors = ann_colors,
  main = "BANP altered motifs across 3d+4d noOV tumor samples\nFirst 50 motifs | 1 = altered (|Delta beta| >= 20), 0 = not altered",
  border_color = NA
)

dev.off()

cat("[DONE] Wrote:", out_pdf, "\n")
EOF

# look at the altered motifs to see if the expression of their linked genes change :
# 1) First generate the long format of the motif-sample matrix to have one line per motif-sample combination for the altered motifs (motif_altered == 1) and then we can link this to the expression data of the linked genes to see if there is a correlation between altered motifs and gene expression changes.
Rscript - <<'EOF'
suppressPackageStartupMessages(library(data.table))

mat_file  <- "./results/methylation/Cpg_altered/NRF1_motif_sample_matrix_noOV.tsv.gz"
out_long  <- "./results/methylation/Cpg_altered/NRF1_motif_sample_matrix_noOV_long.tsv.gz"

dt <- fread(mat_file)

long_dt <- melt(
  dt,
  id.vars = "motif_id",
  variable.name = "sample_id",
  value.name = "motif_altered"
)

long_dt <- long_dt[motif_altered == 1]

fwrite(long_dt, out_long, sep = "\t")
cat("[DONE] Wrote:", out_long, "\n")
cat("Rows:", nrow(long_dt), "\n")
EOF
# 2) Adding gene information to the long format motif-sample table by merging with a motif-to-gene mapping file.
Rscript - <<'EOF'
suppressPackageStartupMessages(library(data.table))

mat_file <- "./results/methylation/Cpg_altered/NRF1_motif_sample_matrix_noOV.tsv.gz"
map_file <- "./methylation/probe_gene_pairs_in_motifs.tsv"
out_file <- "./results/methylation/Cpg_altered/NRF1_motif_sample_gene_pairs_noOV.tsv.gz"

# 1) Load matrix
mat <- fread(mat_file)
mat[, motif_id := trimws(as.character(motif_id))]

# 2) Convert matrix to long altered motif/sample pairs
alt <- melt(
  mat,
  id.vars = "motif_id",
  variable.name = "sample_id",
  value.name = "motif_altered"
)[motif_altered == 1, .(motif_id, sample_id)]

# 3) Load motif -> gene map
map <- fread(map_file)
map[, TF := trimws(as.character(TF))]
map[, motif_instance_id := trimws(as.character(motif_instance_id))]
map[, gene := trimws(as.character(gene))]

map <- map[TF == "NRF1_mm0to2_noCGmm", .(
  motif_id = motif_instance_id,
  gene
)]

# keep one gene per motif
map <- unique(map)
map1 <- map[, .(gene = unique(gene)[1]), by = motif_id]

# 4) Merge = add linked gene column
x <- merge(alt, map1, by = "motif_id", all.x = TRUE)

# 5) Quick checks
cat("Rows in altered long table:", nrow(alt), "\n")
cat("Rows after merge:", nrow(x), "\n")
cat("Rows with gene found:", sum(!is.na(x$gene) & x$gene != ""), "\n")
cat("Rows missing gene:", sum(is.na(x$gene) | x$gene == ""), "\n")

cat("\nExample rows:\n")
print(head(x, 10))

# 6) Save
fwrite(x, out_file, sep = "\t")
cat("\n[DONE] Wrote:", out_file, "\n")
EOF

# Number of altered motifs for each TF 
zcat ./results/methylation/Cpg_altered/NRF1_motif_sample_matrix_noOV.tsv.gz | awk -F'\t' 'NR>1 {print $1}' | sort | uniq -c | sort -k1,1nr | wc -l 
# 2416
zcat ./results/methylation/Cpg_altered/BANP_motif_sample_matrix_noOV.tsv.gz | awk -F'\t' 'NR>1 {print $1}' | sort | uniq -c | sort -k1,1nr | wc -l 
# 718

# to choose the gene based on the highest number of altered samples linked to the gene (e.g. BANP) we can do :
zcat ./results/methylation/Cpg_altered/NRF1_motif_sample_gene_pairs_noOV.tsv.gz \
| awk -F'\t' 'NR>1 {print $3 "\t" $2}' \
| sort -u \
| cut -f1 \
| sort \
| uniq -c \
| sort -k1,1nr \
| head
zcat ./results/methylation/Cpg_altered/BANP_motif_sample_gene_pairs_noOV.tsv.gz \
| awk -F'\t' 'NR>1 {print $3 "\t" $2}' \
| sort -u \
| cut -f1 \
| sort \
| uniq -c \
| sort -k1,1nr \
| head

# 3) plot the boxplot of the expression of one gene of interest (e.g. TBX15) in the samples with altered motifs vs the  cancer samples without altered motifs across all cancers to see if there is a correlation between altered motifs and gene expression changes. this also draws the p-value of the wilcoxon test comparing the expression of the gene in the altered vs non-altered samples for each cancer type. This will help us understand if the altered motifs are associated with changes in gene expression in the linked gene across different cancers and also draws the density.
Rscript - <<'EOF'
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# =========================
# INPUTS
# =========================
gene_of_interest <- "PCDHGA5"

alt_file  <- "./results/methylation/Cpg_altered/BANP_motif_sample_gene_pairs_noOV.tsv.gz"
expr_file <- "./expression/gene_expression_matrix_3d4d_noOV.tsv"
out_pdf   <- paste0("./results/methylation/Cpg_altered/BANP_", gene_of_interest, "_expression_violin_allCancers_tumorOnly.pdf")

dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)

# =========================
# HELPERS
# =========================
extract_patient_id_expr <- function(x) {
  sub("^TCGA-[^-]+-(TCGA-[^-]+-[^-]+)-.*$", "\\1", x)
}

is_tumor_expr <- function(x) {
  grepl("-01A([_.-]|$)", x)
}

fmt_p <- function(x) {
  ifelse(is.na(x), "NA", formatC(x, format = "e", digits = 2))
}

p_to_star <- function(p) {
  if (is.na(p)) "ns"
  else if (p < 0.001) "***"
  else if (p < 0.01) "**"
  else if (p < 0.05) "*"
  else "ns"
}

# =========================
# LOAD ALTERED SAMPLE-GENE TABLE
# =========================
alt <- fread(alt_file)
alt_gene <- unique(alt[gene == gene_of_interest, .(sample_id)])
alt_gene[, altered := 1L]

cat("Samples with altered motifs linked to", gene_of_interest, ":", nrow(alt_gene), "\n")

# =========================
# LOAD EXPRESSION
# =========================
expr <- fread(expr_file)

if (!("sample" %in% names(expr))) stop("Expression file must contain a 'sample' column.")
if (!(gene_of_interest %in% names(expr))) stop(paste("Gene not found in expression matrix:", gene_of_interest))

plot_dt <- expr[, .(
  expr_sample = sample,
  expression  = suppressWarnings(as.numeric(get(gene_of_interest)))
)]

# keep tumor samples only
plot_dt <- plot_dt[is_tumor_expr(expr_sample)]

# derive IDs
plot_dt[, sample_id := extract_patient_id_expr(expr_sample)]
plot_dt[, cancer := sub("^TCGA-([^-]+)-.*$", "\\1", expr_sample)]

# altered / non-altered
plot_dt <- merge(
  plot_dt,
  alt_gene,
  by = "sample_id",
  all.x = TRUE
)

plot_dt[is.na(altered), altered := 0L]
plot_dt[, altered_label := ifelse(altered == 1L, "Altered motif", "No altered motif")]
plot_dt <- plot_dt[!is.na(expression) & !is.na(cancer)]

cat("Total tumor samples with expression:", nrow(plot_dt), "\n")
cat("Altered tumor samples:", sum(plot_dt$altered == 1), "\n")
cat("Non-altered tumor samples:", sum(plot_dt$altered == 0), "\n")

# keep only cancer types where both groups exist
keep_cancers <- plot_dt[, .(
  n_alt = sum(altered == 1),
  n_non = sum(altered == 0)
), by = cancer][n_alt > 0 & n_non > 0, cancer]

plot_dt <- plot_dt[cancer %in% keep_cancers]

if (nrow(plot_dt) == 0) {
  stop("No cancer type has both altered and non-altered tumor samples for this gene.")
}

# =========================
# WILCOXON TEST PER CANCER
# =========================
wilcox_dt <- rbindlist(lapply(sort(unique(plot_dt$cancer)), function(cc) {
  subdt <- plot_dt[cancer == cc]

  pval <- tryCatch(
    wilcox.test(expression ~ altered_label, data = subdt)$p.value,
    error = function(e) NA_real_
  )

  data.table(
    cancer = cc,
    n_alt = sum(subdt$altered == 1),
    n_non = sum(subdt$altered == 0),
    p_value = pval
  )
}), fill = TRUE)

wilcox_dt[, p_adj := p.adjust(p_value, method = "BH")]
wilcox_dt[, signif_label := vapply(p_adj, p_to_star, character(1))]
wilcox_dt[, facet_label := paste0(
  cancer,
  "\nAlt=", n_alt, " | Non=", n_non,
  "\nBH p=", fmt_p(p_adj)
)]

print(wilcox_dt[order(p_adj)])

# add facet labels
plot_dt <- merge(
  plot_dt,
  wilcox_dt[, .(cancer, facet_label, p_adj, signif_label)],
  by = "cancer",
  all.x = TRUE
)

facet_levels <- wilcox_dt[order(cancer)]$facet_label
plot_dt[, facet_label := factor(facet_label, levels = facet_levels)]

# =========================
# STAR POSITIONS
# =========================
star_dt <- plot_dt[, .(
  y_min = min(expression, na.rm = TRUE),
  y_max = max(expression, na.rm = TRUE)
), by = .(cancer, facet_label)]

star_dt[, y := y_max - 0.08 * (y_max - y_min)]

star_dt <- merge(
  star_dt,
  unique(wilcox_dt[, .(cancer, signif_label)]),
  by = "cancer",
  all.x = TRUE
)

star_dt[, x := 1.5]

# =========================
# PLOT
# =========================
p <- ggplot(plot_dt, aes(x = altered_label, y = expression)) +
  geom_violin(aes(fill = altered_label), trim = FALSE, scale = "width", alpha = 0.5) +
  geom_boxplot(aes(fill = altered_label), width = 0.18, outlier.shape = NA, alpha = 0.9) +
  geom_jitter(aes(color = altered_label), width = 0.12, size = 0.7, alpha = 0.7) +
  geom_text(
    data = star_dt,
    aes(x = x, y = y, label = signif_label),
    inherit.aes = FALSE,
    size = 5
  ) +
  facet_wrap(~ facet_label, scales = "free_y") +
  scale_fill_manual(values = c("No altered motif" = "black", "Altered motif" = "red")) +
  scale_color_manual(values = c("No altered motif" = "black", "Altered motif" = "red")) +
  labs(
    title = paste0(gene_of_interest, " expression by cancer type"),
    subtitle = "Tumor samples only | Red = altered linked motif | Black = no altered linked motif",
    x = "",
    y = "Expression"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

ggsave(out_pdf, p, width = 16, height = 12)
cat("[DONE] Wrote:", out_pdf, "\n")

out_stats <- sub("[.]pdf$", "_wilcoxon_stats.tsv", out_pdf)
fwrite(wilcox_dt, out_stats, sep = "\t")
cat("[DONE] Wrote:", out_stats, "\n")
EOF

# to do the pdf for all genes :
Rscript - <<'EOF'
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# =========================
# INPUTS
# =========================
alt_file   <- "./results/methylation/Cpg_altered/NRF1_motif_sample_gene_pairs_noOV.tsv.gz"
expr_file  <- "./expression/gene_expression_matrix_3d4d_noOV.tsv"
out_pdf    <- "./results/methylation/Cpg_altered/NRF1_sigGenes_expression_violin_allCancers_tumorOnly.pdf"
out_tsv    <- "./results/methylation/Cpg_altered/NRF1_sigGenes_expression_stats_allCancers_tumorOnly.tsv.gz"
out_sig_tsv <- "./results/methylation/Cpg_altered/NRF1_sigOnly_expression_stats_allCancers_tumorOnly.tsv.gz"

dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)

# =========================
# HELPERS
# =========================
extract_patient_id_expr <- function(x) {
  sub("^TCGA-[^-]+-(TCGA-[^-]+-[^-]+)-.*$", "\\1", x)
}

is_tumor_expr <- function(x) {
  grepl("-01A([_.-]|$)", x)
}

fmt_p <- function(x) {
  ifelse(is.na(x), "NA", formatC(x, format = "e", digits = 2))
}

p_to_star <- function(p) {
  if (is.na(p)) "ns"
  else if (p < 0.001) "***"
  else if (p < 0.01) "**"
  else if (p < 0.05) "*"
  else "ns"
}

collapse_unique <- function(x) {
  x <- unique(na.omit(as.character(x)))
  x <- x[nzchar(x)]
  if (length(x) == 0) NA_character_ else paste(sort(x), collapse = ", ")
}

short_motif_label <- function(x, max_show = 3) {
  x <- unique(na.omit(as.character(x)))
  x <- x[nzchar(x)]
  if (length(x) == 0) return("NA")
  if (length(x) <= max_show) return(paste(sort(x), collapse = ", "))
  paste0(paste(sort(x)[1:max_show], collapse = ", "), " ... (n=", length(x), ")")
}

# =========================
# LOAD DATA
# =========================
cat("[1/4] Loading altered table...\n")
alt <- fread(alt_file)
req_alt <- c("motif_id", "sample_id", "gene")
miss_alt <- setdiff(req_alt, names(alt))
if (length(miss_alt)) stop("Missing columns in altered table: ", paste(miss_alt, collapse=", "))
cat("  Rows:", nrow(alt), "\n")
cat("  Unique genes:", uniqueN(alt$gene), "\n")

cat("[2/4] Loading expression matrix...\n")
expr <- fread(expr_file)
if (!("sample" %in% names(expr))) stop("Expression file must contain a 'sample' column.")
cat("  Rows:", nrow(expr), "\n")
cat("  Columns:", ncol(expr), "\n")

cat("[3/4] Building gene list...\n")
genes_all <- sort(intersect(unique(alt$gene), setdiff(names(expr), "sample")))
if (length(genes_all) == 0) stop("No overlapping genes between altered table and expression matrix.")
cat("  Genes to test:", length(genes_all), "\n")

cat("[4/4] Starting analysis...\n")

# =========================
# MAIN LOOP
# =========================
all_stats <- list()
pdf_open <- FALSE
n_pages <- 0L
n_sig_genes <- 0L

for (i in seq_along(genes_all)) {
  g <- genes_all[i]
  cat("\n[Gene", i, "/", length(genes_all), "]", g, "\n")

  alt_gene_full <- alt[gene == g]
  alt_gene <- unique(alt_gene_full[, .(sample_id)])
  alt_gene[, altered := 1L]

  motifs_full  <- collapse_unique(alt_gene_full$motif_id)
  motifs_short <- short_motif_label(alt_gene_full$motif_id, max_show = 3)
  n_motifs     <- uniqueN(alt_gene_full$motif_id)

  cat("  Altered samples:", nrow(alt_gene), "\n")
  cat("  Motif instances:", n_motifs, "\n")

  # IMPORTANT: build directly from expr to avoid length mismatch
  plot_dt <- expr[, .(
    expr_sample = sample,
    expression = suppressWarnings(as.numeric(get(g)))
  )]

  plot_dt <- plot_dt[is_tumor_expr(expr_sample)]
  plot_dt <- plot_dt[!is.na(expression)]
  plot_dt[, sample_id := extract_patient_id_expr(expr_sample)]
  plot_dt[, cancer := sub("^TCGA-([^-]+)-.*$", "\\1", expr_sample)]

  cat("  Tumor samples with expression:", nrow(plot_dt), "\n")

  # merge altered status
  plot_dt <- merge(plot_dt, alt_gene, by = "sample_id", all.x = TRUE)
  plot_dt[is.na(altered), altered := 0L]
  plot_dt[, altered_label := ifelse(altered == 1L, "Altered motif", "No altered motif")]

  # only cancers with both groups
  keep_stats <- plot_dt[, .(
    n_alt = sum(altered == 1L),
    n_non = sum(altered == 0L)
  ), by = cancer]

  keep_cancers <- keep_stats[n_alt > 0 & n_non > 0, cancer]
  plot_dt <- plot_dt[cancer %in% keep_cancers]

  cat("  Cancer types with both groups:", length(keep_cancers), "\n")

  if (nrow(plot_dt) == 0) {
    cat("  -> Skipped: no cancer with both groups\n")
    next
  }

  wilcox_dt <- rbindlist(lapply(sort(unique(plot_dt$cancer)), function(cc) {
    subdt <- plot_dt[cancer == cc]

    pval <- tryCatch(
      wilcox.test(expression ~ altered_label, data = subdt)$p.value,
      error = function(e) NA_real_
    )

    data.table(
      gene = g,
      motifs = motifs_full,
      n_motif_instances = n_motifs,
      cancer = cc,
      n_alt = sum(subdt$altered == 1L),
      n_non = sum(subdt$altered == 0L),
      median_alt = median(subdt[altered == 1L]$expression, na.rm = TRUE),
      median_non = median(subdt[altered == 0L]$expression, na.rm = TRUE),
      mean_alt = mean(subdt[altered == 1L]$expression, na.rm = TRUE),
      mean_non = mean(subdt[altered == 0L]$expression, na.rm = TRUE),
      p_value = pval
    )
  }), fill = TRUE)

  wilcox_dt[, p_adj := p.adjust(p_value, method = "BH")]
  wilcox_dt[, signif_label := vapply(p_adj, p_to_star, character(1))]
  wilcox_dt[, is_significant := !is.na(p_adj) & p_adj < 0.05]
  wilcox_dt[, direction := fifelse(median_alt > median_non, "higher_in_altered",
                            fifelse(median_alt < median_non, "lower_in_altered", "equal_median"))]

  all_stats[[length(all_stats) + 1L]] <- wilcox_dt

  n_sig_this_gene <- sum(wilcox_dt$is_significant, na.rm = TRUE)
  cat("  Significant cancers (BH<0.05):", n_sig_this_gene, "\n")

  if (n_sig_this_gene == 0L) {
    cat("  -> Not added to PDF\n")
    next
  }

  n_sig_genes <- n_sig_genes + 1L

  wilcox_dt[, facet_label := paste0(
    cancer,
    "\nAlt=", n_alt, " | Non=", n_non,
    "\nBH p=", fmt_p(p_adj)
  )]

  plot_dt <- merge(
    plot_dt,
    wilcox_dt[, .(cancer, facet_label, signif_label)],
    by = "cancer",
    all.x = TRUE
  )

  facet_levels <- wilcox_dt[order(cancer)]$facet_label
  plot_dt[, facet_label := factor(facet_label, levels = facet_levels)]

  star_dt <- plot_dt[, .(
    y_min = min(expression, na.rm = TRUE),
    y_max = max(expression, na.rm = TRUE)
  ), by = .(cancer, facet_label)]

  star_dt[, y := y_max - 0.12 * (y_max - y_min)]
  star_dt <- merge(
    star_dt,
    unique(wilcox_dt[, .(cancer, signif_label)]),
    by = "cancer",
    all.x = TRUE
  )
  star_dt[, x := 1.5]

  if (!pdf_open) {
    pdf(out_pdf, width = 16, height = 12)
    pdf_open <- TRUE
  }

  p <- ggplot(plot_dt, aes(x = altered_label, y = expression)) +
    geom_violin(aes(fill = altered_label), trim = FALSE, scale = "width", alpha = 0.5) +
    geom_boxplot(aes(fill = altered_label), width = 0.18, outlier.shape = NA, alpha = 0.9) +
    geom_jitter(aes(color = altered_label), width = 0.12, size = 0.7, alpha = 0.7) +
    geom_text(
      data = star_dt,
      aes(x = x, y = y, label = signif_label),
      inherit.aes = FALSE,
      size = 5,
      fontface = "bold"
    ) +
    facet_wrap(~ facet_label, scales = "free_y") +
    scale_fill_manual(values = c("No altered motif" = "black", "Altered motif" = "red")) +
    scale_color_manual(values = c("No altered motif" = "black", "Altered motif" = "red")) +
    labs(
      title = paste0(g, " | motifs: ", motifs_short),
      subtitle = "Tumor samples only | Gene kept only if significant in at least one cancer type (BH < 0.05)",
      x = "",
      y = "Expression"
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      strip.text = element_text(size = 8),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    )

  print(p)
  n_pages <- n_pages + 1L
  cat("  -> Added to PDF\n")
}

# =========================
# SAVE OUTPUTS
# =========================
stats_dt <- if (length(all_stats)) rbindlist(all_stats, fill = TRUE) else data.table()

if (nrow(stats_dt) > 0) {
  fwrite(stats_dt, out_tsv, sep = "\t")
  fwrite(stats_dt[is_significant == TRUE], out_sig_tsv, sep = "\t")
  cat("\n[STATS] Wrote:", out_tsv, "\n")
  cat("[SIG STATS] Wrote:", out_sig_tsv, "\n")
  cat("[STATS] Rows:", nrow(stats_dt), "\n")
  cat("[SIG STATS] Rows:", nrow(stats_dt[is_significant == TRUE]), "\n")
} else {
  cat("\n[STATS] No statistics to write.\n")
}

if (pdf_open) {
  dev.off()
  cat("[PDF] Wrote:", out_pdf, "\n")
} else {
  cat("[PDF] No gene significant in any cancer. No PDF created.\n")
}

cat("\nSummary:\n")
cat("  Tested genes:", length(genes_all), "\n")
cat("  Genes added to PDF:", n_sig_genes, "\n")
cat("  PDF pages:", n_pages, "\n")
EOF



##########################################################################################
# Redo the altered motifs of above but also separating up and down methylation changes 
##########################################################################################
mkdir -p ./results/methylation/Cpg_altered_direction/
# This script generates a motif-sample matrix for TF where the values indicate the direction of methylation change (1 for higher methylation in altered samples, -1 for lower methylation in altered samples, and 0 for no alteration) based on a specified delta beta threshold. It uses the same cohort and motif mapping as before but now distinguishes between hypermethylation and hypomethylation alterations in the motifs. The output is saved as a gzipped TSV file for downstream analysis and visualization and put in the ./results/methylation/Cpg_altered_direction/ directory.
# it keeps only cancer types with both tumor and normal samples in the cohort to ensure a valid comparison for determining the direction of methylation changes in the motifs.
Rscript - <<'EOF'
suppressPackageStartupMessages(library(data.table))

# =========================
# INPUTS
# =========================
samples_file <- "./results/multi_omics/samples_3d+4d_noOV.tsv"
motif_file   <- "./motifs/overlaps/intersected_motifs2mm_HM450/NRF1_mm0to2_noCGmm_probeXmotif.tsv"
out_file     <- "./results/methylation/Cpg_altered_direction/NRF1_motif_sample_matrix_sameCancer_signed_noOV.tsv.gz"

delta_thr <- 20

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

# =========================
# HELPERS
# =========================
read_meth <- function(f) {
  x <- fread(cmd = paste("zcat", shQuote(f)),
             header = FALSE,
             na.strings = c("NA","","NaN"),
             showProgress = FALSE)
  x[, .(probe = as.character(V4), beta = as.numeric(V5))]
}

extract_sample_type <- function(x) {
  x <- basename(x)
  if (grepl("-01A([_.-]|$)", x)) return("tumor")
  if (grepl("-11A([_.-]|$)", x)) return("normal")
  NA_character_
}

extract_patient_id <- function(x) {
  sub("^HM450_TCGA-[^-]+-(TCGA-[^-]+-[^-]+)-.*$", "\\1", basename(x))
}

extract_cancer <- function(x) {
  sub("^HM450_TCGA-([^-]+)-.*$", "\\1", basename(x))
}

# =========================
# COHORT
# =========================
samples <- fread(samples_file, header = FALSE)
setnames(samples, c("sample","cancer","SNV","METH","RNA","ATAC"))
samples <- samples[METH == 1, .(sample, cancer)]

# =========================
# MOTIF MAP
# =========================
motif_raw <- fread(motif_file, header = FALSE, showProgress = FALSE)
motif_map <- unique(motif_raw[, .(
  probe    = as.character(V4),
  motif_id = paste(V5, V6, V7, V10, sep=":")
)])

cat("Unique motif probes :", uniqueN(motif_map$probe), "\n")
cat("Unique motifs       :", uniqueN(motif_map$motif_id), "\n")

# =========================
# INDEX ALL METH FILES
# =========================
all_files <- list.files(
  "./methylation/filtered_methylation",
  pattern = "annotated_methylation_filtered[.]bed[.]gz$",
  full.names = TRUE
)

file_index <- data.table(
  meth_file   = all_files,
  patient_id  = vapply(all_files, extract_patient_id, character(1)),
  cancer      = vapply(all_files, extract_cancer, character(1)),
  sample_type = vapply(all_files, extract_sample_type, character(1))
)[!is.na(sample_type)]

tumor_index  <- file_index[sample_type == "tumor"]
normal_index <- file_index[sample_type == "normal"]

# =========================
# MATCH TO COHORT
# =========================
tumor_samples <- merge(
  samples, tumor_index,
  by.x = c("sample","cancer"),
  by.y = c("patient_id","cancer")
)

normal_samples <- merge(
  samples, normal_index,
  by.x = c("sample","cancer"),
  by.y = c("patient_id","cancer")
)

if (nrow(tumor_samples) == 0) stop("No tumor methylation files matched to cohort.")
if (nrow(normal_samples) == 0) stop("No normal methylation files matched to cohort.")

# =========================
# KEEP ONLY CANCERS WITH BOTH TUMOR + NORMAL
# =========================
eligible_cancers <- intersect(unique(tumor_samples$cancer), unique(normal_samples$cancer))

tumor_samples  <- tumor_samples[cancer %in% eligible_cancers]
normal_samples <- normal_samples[cancer %in% eligible_cancers]

if (length(eligible_cancers) == 0) {
  stop("No cancer types with both tumor and normal samples in the cohort.")
}

cat("Eligible cancers:", length(eligible_cancers), "\n")
cat("Tumor files matched :", nrow(tumor_samples), "\n")
cat("Normal files matched:", nrow(normal_samples), "\n")

# =========================
# BUILD SAME-CANCER NORMAL REFERENCE
# =========================
cat("Reading normal methylation files...\n")
t_ref1 <- Sys.time()

norm_list <- vector("list", nrow(normal_samples))

for (i in seq_len(nrow(normal_samples))) {
  if (i %% 10 == 0 || i == 1 || i == nrow(normal_samples)) {
    cat("[", i, "/", nrow(normal_samples), "] normals\n", sep = "")
  }

  f  <- normal_samples$meth_file[i]
  cc <- normal_samples$cancer[i]

  dt <- read_meth(f)
  dt <- motif_map[dt, on = "probe", nomatch = 0]
  dt <- dt[!is.na(beta)]
  if (nrow(dt) == 0) next

  norm_list[[i]] <- dt[, .(cancer = cc, probe, beta)]
}

norm_dt <- rbindlist(norm_list, fill = TRUE)
if (nrow(norm_dt) == 0) stop("No motif-overlapping normal beta values found.")

ref_by_cancer <- norm_dt[, .(ref_beta = median(beta, na.rm = TRUE)), by = .(cancer, probe)]

t_ref2 <- Sys.time()
cat("Reference build time:",
    round(as.numeric(difftime(t_ref2, t_ref1, units = "mins")), 2),
    "minutes\n")

# =========================
# PROCESS TUMOR SAMPLES
# RULE:
# for each motif_id x sample_id, keep probe with max abs(delta)
# state = +1 if delta >= threshold
# state = -1 if delta <= -threshold
# state = 0 otherwise
# =========================
cat("Processing tumor samples...\n")
t_tum1 <- Sys.time()

res <- vector("list", nrow(tumor_samples))

for (i in seq_len(nrow(tumor_samples))) {
  s <- tumor_samples[i]
  sample_id <- s$sample
  ccur <- s$cancer

  if (i %% 25 == 0 || i == 1 || i == nrow(tumor_samples)) {
    cat("[", i, "/", nrow(tumor_samples), "] ", sample_id, " ", ccur, "\n", sep = "")
  }

  dt <- read_meth(s$meth_file)
  dt <- motif_map[dt, on = "probe", nomatch = 0]
  dt <- dt[!is.na(beta)]
  if (nrow(dt) == 0) next

  ref_dt <- ref_by_cancer[cancer == ccur, .(probe, ref_beta)]
  if (nrow(ref_dt) == 0) next

  dt <- ref_dt[dt, on = "probe", nomatch = 0]
  if (nrow(dt) == 0) next

  dt[, delta := beta - ref_beta]

  # keep strongest probe-level signal per motif in this sample
  strongest <- dt[, .SD[which.max(abs(delta))], by = motif_id]

  strongest[, state := fifelse(delta >=  delta_thr,  1L,
                        fifelse(delta <= -delta_thr, -1L, 0L))]

  strongest <- strongest[state != 0L, .(motif_id, sample_id, state)]

  if (nrow(strongest) > 0) {
    res[[i]] <- strongest
  }
}

motif_sample_dt <- rbindlist(res, fill = TRUE)

if (nrow(motif_sample_dt) == 0) {
  stop("No altered motifs detected.")
}

# safety: remove duplicate motif/sample if any
motif_sample_dt <- unique(motif_sample_dt, by = c("motif_id", "sample_id"))

motif_matrix <- dcast(
  motif_sample_dt,
  motif_id ~ sample_id,
  value.var = "state",
  fill = 0
)

fwrite(motif_matrix, out_file, sep = "\t")
cat("[DONE] Wrote:", out_file, "\n")

t_tum2 <- Sys.time()
cat("Tumor processing time:",
    round(as.numeric(difftime(t_tum2, t_tum1, units = "mins")), 2),
    "minutes\n")
cat("Total time:",
    round(as.numeric(difftime(t_tum2, t_ref1, units = "mins")), 2),
    "minutes\n")

EOF
####################################
# plot the heatmap of all motifs 
####################################
Rscript - <<'EOF'
suppressPackageStartupMessages({
  library(data.table)
  library(pheatmap)
})

mat_file     <- "./results/methylation/Cpg_altered_direction/BANP_motif_sample_matrix_sameCancer_signed_noOV.tsv.gz"
samples_file <- "./results/multi_omics/samples_3d+4d_noOV.tsv"
col_file     <- "./results/multi_omics/cancer_color_order_with_defined_colours.tsv"
out_pdf      <- "./results/methylation/Cpg_altered_direction/Heatmap_BANP_motif_sample_matrix_allMotifs_directional.pdf"

dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)

# =========================
# LOAD MATRIX
# =========================
mat_dt <- fread(mat_file)
if (ncol(mat_dt) < 2) stop("Matrix file has no sample columns: ", mat_file)

motif_ids <- mat_dt$motif_id
mat <- as.matrix(mat_dt[, -1, with = FALSE])
rownames(mat) <- motif_ids
mode(mat) <- "numeric"

# optional safety check
vals <- sort(unique(as.numeric(mat)))
cat("Unique matrix values:", paste(vals, collapse = ", "), "\n")

cat("Number of motifs plotted :", nrow(mat), "\n")
cat("Number of samples plotted:", ncol(mat), "\n")

# =========================
# LOAD SAMPLE ANNOTATION
# =========================
samples <- fread(samples_file, header = FALSE)
setnames(samples, c("sample", "cancer", "SNV", "METH", "RNA", "ATAC"))
samples <- unique(samples[METH == 1, .(sample, cancer)])

ann <- samples[sample %in% colnames(mat)]
if (nrow(ann) == 0) stop("No matrix columns matched sample annotation.")

ann[, Cancer := as.character(cancer)]

# keep only annotated samples and order matrix accordingly
sample_order <- ann$sample[ann$sample %in% colnames(mat)]
mat <- mat[, sample_order, drop = FALSE]

# =========================
# LOAD CANCER COLORS / ORDER
# =========================
cols <- fread(col_file, header = FALSE)[, 1:2]
setnames(cols, c("Cancer", "Color"))
cols[, Cancer := gsub("\r", "", trimws(as.character(Cancer)))]
cols[, Color  := gsub("\r", "", trimws(as.character(Color)))]
cols[, Cancer := sub("^TCGA-", "", Cancer)]

cancer_levels <- unique(cols$Cancer)
cancer_cols_all <- setNames(cols$Color, cols$Cancer)

missing_cancers <- setdiff(unique(ann$Cancer), names(cancer_cols_all))
if (length(missing_cancers) > 0) {
  extra_cols <- grDevices::rainbow(length(missing_cancers))
  names(extra_cols) <- missing_cancers
  cancer_cols_all <- c(cancer_cols_all, extra_cols)
  cancer_levels <- c(cancer_levels, missing_cancers)
}

ann[, Cancer := factor(Cancer, levels = cancer_levels)]
setorder(ann, Cancer, sample)

sample_order <- ann$sample
mat <- mat[, sample_order, drop = FALSE]

ann_df <- data.frame(
  Cancer = ann$Cancer,
  row.names = ann$sample,
  stringsAsFactors = FALSE
)

# =========================
# COLUMN GROUPS / LABELS
# =========================
group_sizes <- ann[, .N, by = Cancer][order(match(Cancer, cancer_levels))]$N
gaps_col <- cumsum(group_sizes)
gaps_col <- gaps_col[gaps_col < ncol(mat)]

labels_col <- rep("", ncol(mat))
group_dt <- ann[, .(
  start = .I[1],
  end   = .I[.N]
), by = Cancer]
group_dt[, mid := floor((start + end) / 2)]
labels_col[group_dt$mid] <- as.character(group_dt$Cancer)

ann_colors <- list(Cancer = cancer_cols_all)

# =========================
# PLOT
# =========================
pdf(out_pdf, width = 20, height = 30)

pheatmap(
  mat,
  color = c("blue", "white", "red"),
  breaks = c(-1.5, -0.5, 0.5, 1.5),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  gaps_col = gaps_col,
  show_rownames = FALSE,
  show_colnames = TRUE,
  labels_col = labels_col,
  angle_col = 90,
  fontsize_col = 8,
  annotation_col = ann_df,
  annotation_colors = ann_colors,
  main = "BANP altered motifs across 3d+4d noOV tumor samples\n-1 = hypomethylated, 0 = unchanged, +1 = hypermethylated",
  border_color = NA,
  legend_breaks = c(-1, 0, 1),
  legend_labels = c("Hypo", "Unchanged", "Hyper")
)

dev.off()

cat("[DONE] Wrote:", out_pdf, "\n")
EOF

# plot the heatmap of motifs that are altered in at least % of samples of at least one cancer type
Rscript - <<'EOF'
suppressPackageStartupMessages({
  library(data.table)
  library(pheatmap)
})

mat_file     <- "./results/methylation/Cpg_altered_direction/BANP_motif_sample_matrix_sameCancer_signed_noOV.tsv.gz"
samples_file <- "./results/multi_omics/samples_3d+4d_noOV.tsv"
col_file     <- "./results/multi_omics/cancer_color_order_with_defined_colours.tsv"
out_pdf      <- "./results/methylation/Cpg_altered_direction/Heatmap_BANP_motif_sample_matrix_min50pctOneCancer.pdf"

dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)

min_frac <- 0.50   # keep motif if altered in >=50% of samples in at least one cancer

cat("[1/6] Loading matrix...\n")

# =========================
# LOAD MATRIX
# =========================
mat_dt <- fread(mat_file)
if (ncol(mat_dt) < 2) stop("Matrix file has no sample columns: ", mat_file)

motif_ids <- mat_dt$motif_id
mat <- as.matrix(mat_dt[, -1, with = FALSE])
rownames(mat) <- motif_ids
mode(mat) <- "numeric"

vals <- sort(unique(as.numeric(mat)))
cat("Unique matrix values:", paste(vals, collapse = ", "), "\n")
cat("Number of motifs loaded :", nrow(mat), "\n")
cat("Number of samples loaded:", ncol(mat), "\n")
cat("Total NA values in matrix:", sum(is.na(mat)), "\n")

cat("[2/6] Loading sample annotation...\n")

# =========================
# LOAD SAMPLE ANNOTATION
# =========================
samples <- fread(samples_file, header = FALSE)
setnames(samples, c("sample", "cancer", "SNV", "METH", "RNA", "ATAC"))
samples <- unique(samples[METH == 1, .(sample, cancer)])

ann <- samples[sample %in% colnames(mat)]
if (nrow(ann) == 0) stop("No matrix columns matched sample annotation.")

ann[, Cancer := as.character(cancer)]

# keep only annotated samples
sample_order <- ann$sample[ann$sample %in% colnames(mat)]
mat <- mat[, sample_order, drop = FALSE]

cat("[3/6] Loading cancer colors and ordering samples...\n")

# =========================
# LOAD CANCER COLORS / ORDER
# =========================
cols <- fread(col_file, header = FALSE)[, 1:2]
setnames(cols, c("Cancer", "Color"))
cols[, Cancer := gsub("\r", "", trimws(as.character(Cancer)))]
cols[, Color  := gsub("\r", "", trimws(as.character(Color)))]
cols[, Cancer := sub("^TCGA-", "", Cancer)]

cancer_levels <- unique(cols$Cancer)
cancer_cols_all <- setNames(cols$Color, cols$Cancer)

missing_cancers <- setdiff(unique(ann$Cancer), names(cancer_cols_all))
if (length(missing_cancers) > 0) {
  extra_cols <- grDevices::rainbow(length(missing_cancers))
  names(extra_cols) <- missing_cancers
  cancer_cols_all <- c(cancer_cols_all, extra_cols)
  cancer_levels <- c(cancer_levels, missing_cancers)
}

ann[, Cancer := factor(Cancer, levels = cancer_levels)]
setorder(ann, Cancer, sample)

sample_order <- ann$sample
mat <- mat[, sample_order, drop = FALSE]

if (!identical(colnames(mat), ann$sample)) {
  stop("Matrix columns do not match annotation sample order.")
}

cat("[4/6] Filtering motifs by cancer-specific frequency...\n")

# =========================
# FILTER MOTIFS:
# keep motif if altered (!=0) in >=50% of samples of at least one cancer type
# =========================
cancer_sample_idx <- split(seq_len(nrow(ann)), ann$Cancer)

n_motifs <- nrow(mat)
keep_motif <- logical(n_motifs)

for (i in seq_len(n_motifs)) {
  cat(sprintf("Processing motif %d / %d\r", i, n_motifs))
  x <- mat[i, ]

  keep_motif[i] <- any(vapply(cancer_sample_idx, function(idx) {
    vals <- x[idx]
    n_valid <- sum(!is.na(vals))
    if (n_valid == 0) return(FALSE)
    frac_alt <- sum(vals != 0, na.rm = TRUE) / n_valid
    frac_alt >= min_frac
  }, logical(1)), na.rm = TRUE)
}
cat("\n")

cat("Motifs before filter:", nrow(mat), "\n")
cat("Motifs kept         :", sum(keep_motif, na.rm = TRUE), "\n")
cat("Motifs removed      :", sum(!keep_motif, na.rm = TRUE), "\n")

mat <- mat[keep_motif, , drop = FALSE]

if (nrow(mat) == 0) {
  stop("No motifs passed the >=50% in at least one cancer type filter.")
}

n_kept_motifs <- nrow(mat)

cat("[5/6] Preparing annotations and labels...\n")

# =========================
# COLUMN ANNOTATION
# =========================
ann_df <- data.frame(
  Cancer = ann$Cancer,
  row.names = ann$sample,
  stringsAsFactors = FALSE
)

group_sizes <- ann[, .N, by = Cancer][order(match(Cancer, cancer_levels))]$N
gaps_col <- cumsum(group_sizes)
gaps_col <- gaps_col[gaps_col < ncol(mat)]

labels_col <- rep("", ncol(mat))
group_dt <- ann[, .(
  start = .I[1],
  end   = .I[.N]
), by = Cancer]
group_dt[, mid := floor((start + end) / 2)]
labels_col[group_dt$mid] <- as.character(group_dt$Cancer)

ann_colors <- list(Cancer = cancer_cols_all)

cat("[6/6] Plotting heatmap...\n")

# =========================
# PLOT
# =========================
pdf(out_pdf, width = 22, height = max(12, nrow(mat) * 0.18))

pheatmap(
  mat,
  color = c("blue", "white", "red"),
  breaks = c(-1.5, -0.5, 0.5, 1.5),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  gaps_col = gaps_col,
  show_rownames = TRUE,
  show_colnames = TRUE,
  labels_row = rownames(mat),
  labels_col = labels_col,
  angle_col = 90,
  fontsize_row = max(4, min(8, 100 / max(1, nrow(mat)))),
  fontsize_col = 8,
  annotation_col = ann_df,
  annotation_colors = ann_colors,
  main = paste0(
    "BANP altered motifs across 3d+4d noOV tumor samples\n",
    "Retained motifs: ", n_kept_motifs, " / ", length(motif_ids), "\n",
    "Filtered: altered in >= ", min_frac * 100, "% of samples in at least one cancer type\n",
    "-1 = hypomethylated, 0 = unchanged, +1 = hypermethylated"
  ),
  border_color = NA,
  legend_breaks = c(-1, 0, 1),
  legend_labels = c("Hypo", "Unchanged", "Hyper")
)

dev.off()

cat("[DONE] Wrote:", out_pdf, "\n")
EOF

# reshape the matrix to long format for easier plotting of alteration frequencies in each cancer type and direction (hypo vs hyper)
Rscript - <<'EOF'
suppressPackageStartupMessages(library(data.table))

mat_file <- "./results/methylation/Cpg_altered_direction/BANP_motif_sample_matrix_sameCancer_signed_noOV.tsv.gz"
out_long <- "./results/methylation/Cpg_altered_direction/BANP_motif_sample_matrix_sameCancer_signed_noOV_long.tsv.gz"

dt <- fread(mat_file)
dt[, motif_id := trimws(as.character(motif_id))]

long_dt <- melt(
  dt,
  id.vars = "motif_id",
  variable.name = "sample_id",
  value.name = "motif_altered"
)

# keep all altered motifs: -1 or +1
long_dt <- long_dt[!is.na(motif_altered) & motif_altered != 0]

fwrite(long_dt, out_long, sep = "\t")
cat("[DONE] Wrote:", out_long, "\n")
cat("Rows:", nrow(long_dt), "\n")
cat("Hypo rows:", sum(long_dt$motif_altered == -1, na.rm = TRUE), "\n")
cat("Hyper rows:", sum(long_dt$motif_altered == 1, na.rm = TRUE), "\n")
EOF

##############################
# add linked gene information
############################## 
Rscript - <<'EOF'
suppressPackageStartupMessages(library(data.table))

mat_file  <- "./results/methylation/Cpg_altered_direction/NRF1_motif_sample_matrix_sameCancer_signed_noOV.tsv.gz"
map_file  <- "./methylation/probe_gene_pairs_in_motifs.tsv"
dist_file <- "./methylation/dist_to_tss/HM450_TCGA-ACC-TCGA-OR-A5J2-01A_1_annotated_methylation_filtered_cpg_dist_tss.tsv"
out_file  <- "./results/methylation/Cpg_altered_direction/NRF1_motif_sample_gene_pairs_sameCancer_signed_noOV.tsv.gz"

prox_thr <- 2000

collapse_unique <- function(x) {
  x <- unique(trimws(as.character(x)))
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) return(NA_character_)
  paste(sort(x), collapse = ",")
}

# =========================
# 1) Load matrix
# =========================
mat <- fread(mat_file)
mat[, motif_id := trimws(as.character(motif_id))]

# =========================
# 2) Convert matrix to long altered motif/sample pairs
# =========================
alt <- melt(
  mat,
  id.vars = "motif_id",
  variable.name = "sample_id",
  value.name = "motif_altered"
)[!is.na(motif_altered) & motif_altered != 0, .(motif_id, sample_id, motif_altered)]

alt[, motif_id := trimws(as.character(motif_id))]
alt[, sample_id := trimws(as.character(sample_id))]
alt[, motif_altered := as.integer(motif_altered)]

# =========================
# 3) Load motif -> probe -> gene map
# =========================
map <- fread(map_file)
names(map) <- trimws(names(map))

req_map <- c("TF", "probe", "gene", "motif_instance_id")
miss_map <- setdiff(req_map, names(map))
if (length(miss_map)) stop("Missing columns in map_file: ", paste(miss_map, collapse = ", "))

map[, TF := trimws(as.character(TF))]
map[, probe := trimws(as.character(probe))]
map[, gene := trimws(as.character(gene))]
map[, motif_instance_id := trimws(as.character(motif_instance_id))]

map <- map[TF == "NRF1_mm0to2_noCGmm", .(
  motif_id = motif_instance_id,
  probe,
  gene
)]

map <- unique(map)

# =========================
# 4) Load probe -> TSS distance
# =========================
dist_dt <- fread(dist_file, header = FALSE, col.names = c("probe", "dist_to_tss"))
dist_dt[, probe := trimws(as.character(probe))]
dist_dt[, dist_to_tss := as.numeric(dist_to_tss)]
dist_dt[, is_proximal_probe := !is.na(dist_to_tss) & abs(dist_to_tss) <= prox_thr]
dist_dt <- unique(dist_dt, by = "probe")

# =========================
# 5) Annotate map with TSS distance
# =========================
map <- merge(map, dist_dt, by = "probe", all.x = TRUE)

# =========================
# 6) Collapse to one row per motif
# =========================
map_motif <- map[, .(
  gene = unique(gene)[1],
  probes = collapse_unique(probe),
  n_probes = uniqueN(probe),
  probe_distances = collapse_unique(
    ifelse(is.na(dist_to_tss), NA_character_, as.character(dist_to_tss))
  ),
  min_abs_dist_to_tss = if (all(is.na(dist_to_tss))) NA_real_ else min(abs(dist_to_tss), na.rm = TRUE),
  is_proximal = any(is_proximal_probe %in% TRUE, na.rm = TRUE)
), by = motif_id]

# optional: check motifs with multiple genes
motif_gene_check <- map[, .(n_genes = uniqueN(gene)), by = motif_id]
cat("Motifs mapping to >1 gene:", sum(motif_gene_check$n_genes > 1), "\n")

# =========================
# 7) Merge with altered motif/sample pairs
# =========================
x <- merge(alt, map_motif, by = "motif_id", all.x = TRUE)

# optional readable label
x[, motif_direction := fifelse(
  motif_altered == -1L, "hypo",
  fifelse(motif_altered == 1L, "hyper", NA_character_)
)]

setcolorder(x, c(
  "motif_id", "gene", "sample_id", "motif_altered", "motif_direction",
  "probes", "n_probes", "probe_distances", "min_abs_dist_to_tss", "is_proximal"
))

# =========================
# 8) Checks
# =========================
cat("Rows in altered long table:", nrow(alt), "\n")
cat("Rows after collapsed merge:", nrow(x), "\n")
cat("Rows with gene found:", sum(!is.na(x$gene) & x$gene != ""), "\n")
cat("Rows proximal motifs:", sum(x$is_proximal %in% TRUE, na.rm = TRUE), "\n")
cat("Rows distal motifs:", sum(x$is_proximal %in% FALSE, na.rm = TRUE), "\n")
cat("Hypo rows:", sum(x$motif_altered == -1L, na.rm = TRUE), "\n")
cat("Hyper rows:", sum(x$motif_altered == 1L, na.rm = TRUE), "\n")

cat("\nExample rows:\n")
print(head(x, 10))

# =========================
# 9) Save
# =========================
fwrite(x, out_file, sep = "\t")
cat("\n[DONE] Wrote:", out_file, "\n")
EOF

########################################################################################################################################
# plot the expression of genes 
# This script tests whether expression of genes linked to altered proximal NRF1 motifs differs in tumors according to motif methylation status. For each gene, all altered proximal motif instances are collapsed at the sample level into hypo, hyper, unchanged, or mixed status. Mixed cases are excluded from testing and saved separately. Expression is measured as tumor-only log2 fold change relative to the median of unchanged tumors within the same cancer. Pairwise Wilcoxon tests compare hypo vs unchanged and hyper vs unchanged for each cancer, with BH correction. Only genes with at least 30% altered samples in at least one cancer are retained, and only genes significant in at least one cancer are plotted.
########################################################################################################################################
Rscript - <<'EOF'
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# =========================
# INPUTS
# =========================
alt_file      <- "./results/methylation/Cpg_altered_direction/NRF1_motif_sample_gene_pairs_sameCancer_signed_noOV.tsv.gz"
expr_file     <- "./expression/gene_expression_matrix_3d4d_noOV.tsv"
out_pdf       <- "./results/methylation/Cpg_altered_direction/NRF1_sigGenes_log2FC_violin_allCancers_tumorOnly_pairwiseWilcox_geneMin30pctAnyAltered_PROXIMALonly.pdf"
out_tsv       <- "./results/methylation/Cpg_altered_direction/NRF1_sigGenes_log2FC_stats_allCancers_tumorOnly_pairwiseWilcox_geneMin30pctAnyAltered_PROXIMALonly.tsv.gz"
out_sig_tsv   <- "./results/methylation/Cpg_altered_direction/NRF1_sigOnly_log2FC_stats_allCancers_tumorOnly_pairwiseWilcox_geneMin30pctAnyAltered_PROXIMALonly.tsv.gz"
out_mixed_tsv <- "./results/methylation/Cpg_altered_direction/NRF1_mixed_gene_sample_cases_with_cancer_PROXIMALonly.tsv.gz"

dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)

min_frac_gene     <- 0.30
pseudocount       <- 1
use_proximal_only <- TRUE

# =========================
# HELPERS
# =========================
extract_patient_id_expr <- function(x) {
  sub("^TCGA-[^-]+-(TCGA-[^-]+-[^-]+)-.*$", "\\1", x)
}

is_tumor_expr <- function(x) {
  grepl("-01A([_.-]|$)", x)
}

fmt_p <- function(x) {
  ifelse(is.na(x), "NA", formatC(x, format = "e", digits = 2))
}

p_to_star <- function(p) {
  if (is.na(p)) "ns"
  else if (p < 0.001) "***"
  else if (p < 0.01) "**"
  else if (p < 0.05) "*"
  else "ns"
}

collapse_unique <- function(x) {
  x <- unique(na.omit(as.character(x)))
  x <- trimws(x)
  x <- x[nzchar(x)]
  if (length(x) == 0) NA_character_ else paste(sort(x), collapse = ", ")
}

short_label <- function(x, max_show = 3) {
  x <- unique(na.omit(as.character(x)))
  x <- trimws(x)
  x <- x[nzchar(x)]
  if (length(x) == 0) return("NA")
  if (length(x) <= max_show) return(paste(sort(x), collapse = ", "))
  paste0(paste(sort(x)[1:max_show], collapse = ", "), " ... (n=", length(x), ")")
}

collapse_status <- function(x) {
  x <- unique(na.omit(as.integer(x)))
  x <- x[x != 0]
  if (length(x) == 0) return(0L)
  if (all(x == 1L)) return(1L)
  if (all(x == -1L)) return(-1L)
  return(2L)  # mixed hypo+hyper
}

status_to_label <- function(x) {
  fifelse(
    x == -1L, "Hypo altered motif",
    fifelse(
      x == 0L, "Unchanged motif",
      fifelse(
        x == 1L, "Hyper altered motif",
        "Mixed altered motif"
      )
    )
  )
}

safe_wilcox <- function(x, y) {
  if (length(x) == 0 || length(y) == 0) return(NA_real_)
  if (all(is.na(x)) || all(is.na(y))) return(NA_real_)
  tryCatch(
    wilcox.test(x, y)$p.value,
    error = function(e) NA_real_
  )
}

# =========================
# LOAD DATA
# =========================
cat("[1/5] Loading altered table...\n")
alt <- fread(alt_file)

req_alt <- c("motif_id", "sample_id", "gene", "motif_altered")
miss_alt <- setdiff(req_alt, names(alt))
if (length(miss_alt)) stop("Missing columns in altered table: ", paste(miss_alt, collapse = ", "))

alt[, motif_altered := as.integer(motif_altered)]
alt[, motif_id := as.character(motif_id)]
alt[, gene := as.character(gene)]
alt[, sample_id := as.character(sample_id)]

if ("is_proximal" %in% names(alt)) {
  alt[, is_proximal := as.logical(is_proximal)]
} else {
  stop("Column 'is_proximal' not found in altered table. Rebuild upstream file first.")
}

cat("  Rows before proximal filter:", nrow(alt), "\n")
cat("  Unique genes before proximal filter:", uniqueN(alt$gene), "\n")

if (use_proximal_only) {
  alt <- alt[!is.na(is_proximal) & is_proximal == TRUE]
  cat("  Rows after proximal filter:", nrow(alt), "\n")
  cat("  Unique genes after proximal filter:", uniqueN(alt$gene), "\n")
}

if (nrow(alt) == 0) stop("No rows left after proximal filter.")

cat("[2/5] Loading expression matrix...\n")
expr <- fread(expr_file)
if (!("sample" %in% names(expr))) stop("Expression file must contain a 'sample' column.")
cat("  Rows:", nrow(expr), "\n")
cat("  Columns:", ncol(expr), "\n")

cat("[3/5] Preparing tumor sample -> cancer map from expression...\n")
expr_samples <- unique(expr[, .(expr_sample = sample)])
expr_samples <- expr_samples[is_tumor_expr(expr_sample)]
expr_samples[, sample_id := extract_patient_id_expr(expr_sample)]
expr_samples[, cancer := sub("^TCGA-([^-]+)-.*$", "\\1", expr_sample)]
sample_cancer_map <- unique(expr_samples[, .(sample_id, cancer)])
cat("  Tumor sample-cancer pairs:", nrow(sample_cancer_map), "\n")

cat("[4/5] Building mixed motif cases table...\n")

# =========================
# MIXED GENE-SAMPLE CASES
# =========================
extra_cols <- intersect(c("probes", "n_probes", "min_abs_dist_to_tss"), names(alt))

mixed_cases <- alt[, {
  base_list <- list(
    motif_status   = collapse_status(motif_altered),
    n_hypo_motifs  = sum(motif_altered == -1L, na.rm = TRUE),
    n_hyper_motifs = sum(motif_altered == 1L, na.rm = TRUE),
    motif_ids_hypo = collapse_unique(motif_id[motif_altered == -1L]),
    motif_ids_hyper = collapse_unique(motif_id[motif_altered == 1L]),
    motif_ids_all  = collapse_unique(motif_id)
  )

  if ("probes" %in% extra_cols) {
    base_list$probes_all <- collapse_unique(probes)
  }
  if ("n_probes" %in% extra_cols) {
    base_list$total_n_probes <- sum(as.numeric(n_probes), na.rm = TRUE)
  }
  if ("min_abs_dist_to_tss" %in% extra_cols) {
    mad <- suppressWarnings(as.numeric(min_abs_dist_to_tss))
    base_list$min_abs_dist_to_tss <- if (all(is.na(mad))) NA_real_ else min(mad, na.rm = TRUE)
  }

  base_list
}, by = .(gene, sample_id)]

mixed_cases <- mixed_cases[motif_status == 2L]
mixed_cases <- merge(mixed_cases, sample_cancer_map, by = "sample_id", all.x = TRUE)

core_order <- c(
  "gene", "sample_id", "cancer", "motif_status",
  "n_hypo_motifs", "n_hyper_motifs",
  "motif_ids_hypo", "motif_ids_hyper", "motif_ids_all"
)
extra_order <- intersect(c("probes_all", "total_n_probes", "min_abs_dist_to_tss"), names(mixed_cases))
setcolorder(mixed_cases, c(core_order, extra_order))

fwrite(mixed_cases, out_mixed_tsv, sep = "\t")
cat("  Mixed gene-sample cases:", nrow(mixed_cases), "\n")
cat("  [DONE] Wrote mixed cases:", out_mixed_tsv, "\n")

cat("[5/5] Starting expression analysis...\n")

# =========================
# GENE LIST
# =========================
genes_all <- sort(intersect(unique(alt$gene), setdiff(names(expr), "sample")))
if (length(genes_all) == 0) stop("No overlapping genes between altered table and expression matrix.")
cat("  Genes to test:", length(genes_all), "\n")

# =========================
# MAIN LOOP
# =========================
all_stats <- list()
pdf_open <- FALSE
n_pages <- 0L
n_sig_genes <- 0L

for (i in seq_along(genes_all)) {
  g <- genes_all[i]
  cat("\n[Gene", i, "/", length(genes_all), "]", g, "\n")

  alt_gene_full <- alt[gene == g]
  if (nrow(alt_gene_full) == 0) next

  motifs_full  <- collapse_unique(alt_gene_full$motif_id)
  motifs_short <- short_label(alt_gene_full$motif_id, max_show = 3)
  n_motifs     <- uniqueN(alt_gene_full$motif_id)

  probes_full  <- if ("probes" %in% names(alt_gene_full)) collapse_unique(alt_gene_full$probes) else NA_character_
  probes_short <- if ("probes" %in% names(alt_gene_full)) short_label(unlist(strsplit(collapse_unique(alt_gene_full$probes), ",\\s*")), max_show = 3) else NA_character_

  alt_gene <- alt_gene_full[, .(
    motif_status = collapse_status(motif_altered)
  ), by = sample_id]

  cat("  Samples linked to altered proximal motifs:", nrow(alt_gene), "\n")
  cat("  Motif instances:", n_motifs, "\n")
  cat("  Hypo samples:", sum(alt_gene$motif_status == -1L), "\n")
  cat("  Hyper samples:", sum(alt_gene$motif_status == 1L), "\n")
  cat("  Mixed samples:", sum(alt_gene$motif_status == 2L), "\n")

  plot_dt <- expr[, .(
    expr_sample = sample,
    raw_expression = suppressWarnings(as.numeric(get(g)))
  )]

  plot_dt <- plot_dt[is_tumor_expr(expr_sample)]
  plot_dt <- plot_dt[!is.na(raw_expression)]
  plot_dt[, sample_id := extract_patient_id_expr(expr_sample)]
  plot_dt[, cancer := sub("^TCGA-([^-]+)-.*$", "\\1", expr_sample)]

  cat("  Tumor samples with expression:", nrow(plot_dt), "\n")

  plot_dt <- merge(plot_dt, alt_gene, by = "sample_id", all.x = TRUE)
  plot_dt[is.na(motif_status), motif_status := 0L]
  plot_dt[, altered_label := status_to_label(motif_status)]

  # exclude mixed from testing
  plot_dt <- plot_dt[motif_status != 2L]

  if (nrow(plot_dt) == 0) {
    cat("  -> Skipped: no usable rows after excluding mixed\n")
    next
  }

  # =========================
  # GENE-LEVEL PRE-FILTER
  # =========================
  gene_prev <- plot_dt[, .(
    n_total = .N,
    n_altered = sum(motif_status != 0L, na.rm = TRUE),
    frac_altered = sum(motif_status != 0L, na.rm = TRUE) / .N
  ), by = cancer]

  max_frac <- if (nrow(gene_prev) > 0) max(gene_prev$frac_altered, na.rm = TRUE) else NA_real_
  cat("  Max altered fraction across cancers:", sprintf("%.3f", max_frac), "\n")

  if (!any(gene_prev$frac_altered >= min_frac_gene, na.rm = TRUE)) {
    cat("  -> Skipped: gene altered in <30% of samples in every cancer type\n")
    next
  }

  # =========================
  # KEEP CANCERS WITH VALID GROUPS
  # =========================
  keep_stats <- plot_dt[, .(
    n_total = .N,
    n_hypo = sum(motif_status == -1L, na.rm = TRUE),
    n_unch = sum(motif_status == 0L, na.rm = TRUE),
    n_hyper = sum(motif_status == 1L, na.rm = TRUE),
    n_altered = sum(motif_status != 0L, na.rm = TRUE),
    frac_altered = sum(motif_status != 0L, na.rm = TRUE) / .N
  ), by = cancer]

  keep_cancers <- keep_stats[
    n_unch > 0 & (n_hypo > 0 | n_hyper > 0),
    cancer
  ]

  plot_dt <- plot_dt[cancer %in% keep_cancers]

  cat("  Cancer types kept:", length(keep_cancers), "\n")

  if (nrow(plot_dt) == 0) {
    cat("  -> Skipped: no cancer has valid comparison groups\n")
    next
  }

  # =========================
  # COMPUTE per-sample log2FC relative to unchanged median in same cancer
  # =========================
  baseline_dt <- plot_dt[motif_status == 0L, .(
    unchanged_median = median(raw_expression, na.rm = TRUE)
  ), by = cancer]

  plot_dt <- merge(plot_dt, baseline_dt, by = "cancer", all.x = TRUE)
  plot_dt <- plot_dt[!is.na(unchanged_median)]

  if (nrow(plot_dt) == 0) {
    cat("  -> Skipped: no unchanged baseline available\n")
    next
  }

  plot_dt[, expression := log2((raw_expression + pseudocount) / (unchanged_median + pseudocount))]

  stats_gene <- rbindlist(lapply(sort(unique(plot_dt$cancer)), function(cc) {
    subdt <- copy(plot_dt[cancer == cc])

    expr_hypo  <- subdt[motif_status == -1L]$expression
    expr_unch  <- subdt[motif_status == 0L]$expression
    expr_hyper <- subdt[motif_status == 1L]$expression

    p_hypo_vs_unch  <- safe_wilcox(expr_hypo, expr_unch)
    p_hyper_vs_unch <- safe_wilcox(expr_hyper, expr_unch)

    data.table(
      gene = g,
      motifs = motifs_full,
      motifs_short = motifs_short,
      n_motif_instances = n_motifs,
      probes = probes_full,
      probes_short = probes_short,
      cancer = cc,
      n_hypo = length(expr_hypo),
      n_unchanged = length(expr_unch),
      n_hyper = length(expr_hyper),
      n_altered = length(expr_hypo) + length(expr_hyper),
      frac_altered = (length(expr_hypo) + length(expr_hyper)) / nrow(subdt),
      unchanged_median_raw = unique(subdt$unchanged_median)[1],
      median_log2fc_hypo = if (length(expr_hypo) > 0) median(expr_hypo, na.rm = TRUE) else NA_real_,
      median_log2fc_unchanged = if (length(expr_unch) > 0) median(expr_unch, na.rm = TRUE) else NA_real_,
      median_log2fc_hyper = if (length(expr_hyper) > 0) median(expr_hyper, na.rm = TRUE) else NA_real_,
      mean_log2fc_hypo = if (length(expr_hypo) > 0) mean(expr_hypo, na.rm = TRUE) else NA_real_,
      mean_log2fc_unchanged = if (length(expr_unch) > 0) mean(expr_unch, na.rm = TRUE) else NA_real_,
      mean_log2fc_hyper = if (length(expr_hyper) > 0) mean(expr_hyper, na.rm = TRUE) else NA_real_,
      p_hypo_vs_unch = p_hypo_vs_unch,
      p_hyper_vs_unch = p_hyper_vs_unch
    )
  }), fill = TRUE)

  stats_gene[, p_adj_hypo_vs_unch := p.adjust(p_hypo_vs_unch, method = "BH")]
  stats_gene[, p_adj_hyper_vs_unch := p.adjust(p_hyper_vs_unch, method = "BH")]

  stats_gene[, signif_hypo_vs_unch := vapply(p_adj_hypo_vs_unch, p_to_star, character(1))]
  stats_gene[, signif_hyper_vs_unch := vapply(p_adj_hyper_vs_unch, p_to_star, character(1))]

  stats_gene[, sig_hypo_vs_unch := !is.na(p_adj_hypo_vs_unch) & p_adj_hypo_vs_unch < 0.05]
  stats_gene[, sig_hyper_vs_unch := !is.na(p_adj_hyper_vs_unch) & p_adj_hyper_vs_unch < 0.05]
  stats_gene[, is_significant := sig_hypo_vs_unch | sig_hyper_vs_unch]

  all_stats[[length(all_stats) + 1L]] <- stats_gene

  n_sig_this_gene <- sum(stats_gene$is_significant, na.rm = TRUE)
  cat("  Significant cancers (BH<0.05 in at least one pairwise test):", n_sig_this_gene, "\n")

  if (n_sig_this_gene == 0L) {
    cat("  -> Not added to PDF\n")
    next
  }

  n_sig_genes <- n_sig_genes + 1L

  stats_gene[, facet_label := paste0(
    cancer,
    "\nAlt=", n_altered, " (", sprintf("%.1f", 100 * frac_altered), "%)",
    " | Hypo=", n_hypo, " | Unch=", n_unchanged, " | Hyper=", n_hyper,
    "\nHypo vs Unch BH=", fmt_p(p_adj_hypo_vs_unch),
    " | Hyper vs Unch BH=", fmt_p(p_adj_hyper_vs_unch)
  )]

  plot_dt <- merge(
    plot_dt,
    stats_gene[, .(
      cancer, facet_label,
      signif_hypo_vs_unch, signif_hyper_vs_unch
    )],
    by = "cancer",
    all.x = TRUE
  )

  plot_dt[, altered_label := factor(
    altered_label,
    levels = c("Hypo altered motif", "Unchanged motif", "Hyper altered motif")
  )]

  facet_levels <- stats_gene[order(cancer)]$facet_label
  plot_dt[, facet_label := factor(facet_label, levels = facet_levels)]

  star_dt <- plot_dt[, .(
    y_min = min(expression, na.rm = TRUE),
    y_max = max(expression, na.rm = TRUE)
  ), by = .(cancer, facet_label)]

  star_dt[, spread := y_max - y_min]
  star_dt[is.na(spread) | spread == 0, spread := 1]
  star_dt[, y1 := y_max - 0.08 * spread]
  star_dt[, y2 := y_max - 0.18 * spread]

  star_dt <- merge(
    star_dt,
    unique(stats_gene[, .(
      cancer, signif_hypo_vs_unch, signif_hyper_vs_unch
    )]),
    by = "cancer",
    all.x = TRUE
  )

  if (!pdf_open) {
    pdf(out_pdf, width = 16, height = 12)
    pdf_open <- TRUE
  }

  title_txt <- paste0(g, " | proximal motifs: ", motifs_short)
  if (!is.na(probes_short)) {
    title_txt <- paste0(title_txt, " | probes: ", probes_short)
  }

  p <- ggplot(plot_dt, aes(x = altered_label, y = expression)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_violin(aes(fill = altered_label), trim = FALSE, scale = "width", alpha = 0.5) +
    geom_boxplot(aes(fill = altered_label), width = 0.18, outlier.shape = NA, alpha = 0.9) +
    geom_jitter(aes(color = altered_label), width = 0.12, size = 0.7, alpha = 0.7) +
    geom_text(
      data = star_dt,
      aes(x = 1.5, y = y1, label = paste0("Hypo vs Unch: ", signif_hypo_vs_unch)),
      inherit.aes = FALSE,
      size = 4,
      fontface = "bold"
    ) +
    geom_text(
      data = star_dt,
      aes(x = 2.5, y = y2, label = paste0("Hyper vs Unch: ", signif_hyper_vs_unch)),
      inherit.aes = FALSE,
      size = 4,
      fontface = "bold"
    ) +
    facet_wrap(~ facet_label, scales = "free_y") +
    scale_fill_manual(values = c(
      "Hypo altered motif" = "blue",
      "Unchanged motif" = "black",
      "Hyper altered motif" = "red"
    )) +
    scale_color_manual(values = c(
      "Hypo altered motif" = "blue",
      "Unchanged motif" = "black",
      "Hyper altered motif" = "red"
    )) +
    labs(
      title = title_txt,
      subtitle = paste0(
        "Tumor samples only | proximal motifs only | y = log2((sample expression + ", pseudocount,
        ") / (median unchanged expression in same cancer + ", pseudocount, ")) | ",
        "Gene kept if any cancer has >= ", 100 * min_frac_gene, "% altered | ",
        "All cancers with valid comparison groups are shown | Mixed cases excluded"
      ),
      x = "",
      y = "log2 fold change vs unchanged median"
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      strip.text = element_text(size = 8),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    )

  print(p)
  n_pages <- n_pages + 1L
  cat("  -> Added to PDF\n")
}

# =========================
# SAVE OUTPUTS
# =========================
stats_dt <- if (length(all_stats)) rbindlist(all_stats, fill = TRUE) else data.table()

if (nrow(stats_dt) > 0) {
  fwrite(stats_dt, out_tsv, sep = "\t")
  fwrite(stats_dt[is_significant == TRUE], out_sig_tsv, sep = "\t")
  cat("\n[STATS] Wrote:", out_tsv, "\n")
  cat("[SIG STATS] Wrote:", out_sig_tsv, "\n")
  cat("[STATS] Rows:", nrow(stats_dt), "\n")
  cat("[SIG STATS] Rows:", nrow(stats_dt[is_significant == TRUE]), "\n")
} else {
  cat("\n[STATS] No statistics to write.\n")
}

if (pdf_open) {
  dev.off()
  cat("[PDF] Wrote:", out_pdf, "\n")
} else {
  cat("[PDF] No gene significant in any cancer. No PDF created.\n")
}

cat("\nSummary:\n")
cat("  Tested genes:", length(genes_all), "\n")
cat("  Genes added to PDF:", n_sig_genes, "\n")
cat("  PDF pages:", n_pages, "\n")
cat("  Mixed cases file:", out_mixed_tsv, "\n")
EOF

#################################################################################################################################################
# plot the correlation across all samples across all cancer between expression of gene and methylation of linked probes --> this code is for one gene at a time, and will output a pdf with 3 pages: 1) all samples, 2) tumor samples only, 3) normal samples only. Each page will have a scatter plot of expression vs methylation, colored by cancer type, and will show the correlation coefficient and p-value in the title. This is to check if the correlation is consistent across all samples, or if it is driven by tumor or normal samples only.
#################################################################################################################################################
# LAUNCHED FOR NRF1 and BANP
###################################################################
# Plot for multiple genes
###################################################################
Rscript - <<'EOF'
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(parallel)
})

# =========================
# INPUTS
# =========================
alt_file    <- "./results/methylation/Cpg_altered_direction/BANP_motif_sample_gene_pairs_sameCancer_signed_noOV.tsv.gz"
expr_file   <- "./expression/gene_expression_matrix_3d4d_noOV.tsv"
meth_dir    <- "./methylation/filtered_methylation"
color_file  <- "./results/multi_omics/cancer_color_order_with_defined_colours.tsv"
out_dir     <- "./results/methylation/Cpg_altered_direction/correlation_gene_expression_methylation/BANP_all_genes"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

use_proximal_only <- TRUE
batch_size <- 20L

# =========================
# HELPERS
# =========================
extract_patient_id_expr <- function(x) {
  out <- sub("^.*(TCGA-[^-]+-[^-]+)-01[A-Z].*$", "\\1", as.character(x))
  out[!grepl("^TCGA-[^-]+-[^-]+$", out)] <- NA_character_
  out
}

is_tumor_expr <- function(x) {
  grepl("-01[A-Z]([_.-]|$)", as.character(x))
}

collapse_status <- function(x) {
  x <- unique(na.omit(as.integer(x)))
  if (length(x) == 0) return(0L)

  x_no0 <- x[x != 0L]
  if (length(x_no0) == 0) return(0L)
  if (all(x_no0 == 1L)) return(1L)
  if (all(x_no0 == -1L)) return(-1L)
  return(2L)
}

status_to_label <- function(x) {
  fifelse(
    x == -1L, "Hypo altered",
    fifelse(
      x == 0L, "Unchanged",
      fifelse(
        x == 1L, "Hyper altered",
        "Mixed"
      )
    )
  )
}

safe_cor <- function(x, y, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(list(r = NA_real_, p = NA_real_))
  z <- tryCatch(cor.test(x[ok], y[ok], method = method), error = function(e) NULL)
  if (is.null(z)) return(list(r = NA_real_, p = NA_real_))
  list(r = unname(z$estimate), p = z$p.value)
}

extract_patient_id_file <- function(x) {
  sub("^HM450_TCGA-[^-]+-(TCGA-[^-]+-[^-]+)-.*$", "\\1", basename(x))
}

fmt_sec <- function(x) {
  x <- as.numeric(x)
  if (!is.finite(x)) return("NA")
  if (x < 60) return(paste0(round(x, 1), " sec"))
  paste0(round(x / 60, 1), " min")
}

read_sample_meth_fast <- function(f, probes_keep) {
  if (!file.exists(f) || length(probes_keep) == 0) return(NULL)

  probes_keep <- unique(trimws(as.character(probes_keep)))
  probes_keep <- probes_keep[nzchar(probes_keep)]
  if (length(probes_keep) == 0) return(NULL)

  probes_esc <- gsub("([][{}()+*^$.|?\\\\-])", "\\\\\\1", probes_keep)
  patt <- paste(probes_esc, collapse = "|")

  cmd <- sprintf(
    "zcat %s | awk -F '\\t' '$4 ~ /^(%s)$/ {print}'",
    shQuote(f),
    patt
  )

  dt <- tryCatch(
    fread(cmd = cmd, sep = "\t", header = FALSE, showProgress = FALSE),
    error = function(e) NULL
  )

  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  if (ncol(dt) < 5) return(NULL)

  out <- data.table(
    probe = as.character(dt[[4]]),
    beta  = suppressWarnings(as.numeric(dt[[5]]))
  )

  out <- out[probe %in% probes_keep & !is.na(beta)]
  if (nrow(out) == 0) return(NULL)
  out
}

# =========================
# LOAD STATIC DATA ONCE
# =========================
cat("[1/6] Loading cancer colors...\n")
cols_dt <- fread(color_file, header = FALSE)
if (ncol(cols_dt) < 2) stop("Color file must have at least 2 columns.")
cols_dt <- cols_dt[, 1:2]
setnames(cols_dt, c("cancer_full", "color"))
cols_dt[, `:=`(
  cancer_full = as.character(cancer_full),
  color = as.character(color)
)]
cols_dt[, cancer := sub("^TCGA-", "", cancer_full)]
cols_dt <- unique(cols_dt[!is.na(cancer) & !is.na(color) & nzchar(cancer) & nzchar(color)])
cancer_cols <- setNames(cols_dt$color, cols_dt$cancer)

cat("[2/6] Loading altered table...\n")
alt <- fread(alt_file)
need_cols <- c("gene", "sample_id", "motif_altered", "probes")
if (!all(need_cols %in% names(alt))) {
  stop("alt_file must contain columns: ", paste(need_cols, collapse = ", "))
}

alt[, `:=`(
  gene = as.character(gene),
  sample_id = as.character(sample_id),
  motif_altered = as.integer(motif_altered),
  probes = as.character(probes)
)]

if (use_proximal_only) {
  if (!("is_proximal" %in% names(alt))) stop("Column 'is_proximal' not found in alt_file.")
  alt[, is_proximal := as.logical(is_proximal)]
  alt <- alt[!is.na(is_proximal) & is_proximal == TRUE]
}

cat("[3/6] Loading expression...\n")
expr <- fread(expr_file)
if (!("sample" %in% names(expr))) stop("Expression file must contain column 'sample'.")

expr_samples <- data.table(
  expr_sample = as.character(expr$sample)
)
expr_samples <- expr_samples[is_tumor_expr(expr_sample)]
expr_samples[, sample_id := extract_patient_id_expr(expr_sample)]
expr_samples[, cancer := sub("^TCGA-([^-]+)-.*$", "\\1", expr_sample)]
expr_samples <- expr_samples[!is.na(sample_id)]
expr_samples <- expr_samples[, .(
  expr_sample = expr_sample[1],
  cancer = cancer[1]
), by = sample_id]

cat("[4/6] Indexing methylation filenames...\n")
meth_files <- list.files(meth_dir, pattern = "\\.bed\\.gz$", full.names = TRUE)
if (length(meth_files) == 0) stop("No methylation .bed.gz files found in: ", meth_dir)

meth_index <- data.table(file = meth_files, fname = basename(meth_files))
meth_index[, sample_id := extract_patient_id_file(fname)]
meth_index <- meth_index[!is.na(sample_id) & grepl("^TCGA-[^-]+-[^-]+$", sample_id)]
setorder(meth_index, fname)
meth_index <- meth_index[, .SD[1], by = sample_id]

all_genes <- sort(unique(alt$gene))
all_genes <- all_genes[!is.na(all_genes) & nzchar(all_genes)]
all_genes <- all_genes[all_genes %in% names(expr)]

cat("[5/6] Genes to process:", length(all_genes), "\n")
cat("[6/6] Starting loop...\n")

summary_list <- vector("list", length(all_genes))

for (g in seq_along(all_genes)) {
  gene_name <- all_genes[g]

  cat("\n====================================================\n")
  cat("[", g, "/", length(all_genes), "] Gene:", gene_name, "\n")
  cat("====================================================\n")

  summary_list[[g]] <- tryCatch({

    alt_gene <- alt[gene == gene_name]
    if (nrow(alt_gene) == 0) stop("No rows found for gene")

    gene_status <- alt_gene[, .(
      motif_status = collapse_status(motif_altered)
    ), by = sample_id]
    gene_status[, alteration := status_to_label(motif_status)]

    probes_keep <- unique(unlist(strsplit(paste(unique(alt_gene$probes), collapse = ","), ",\\s*")))
    probes_keep <- trimws(probes_keep)
    probes_keep <- probes_keep[nzchar(probes_keep)]
    if (length(probes_keep) == 0) stop("No probes found for gene")

    expr_dt <- expr[, .(
      expr_sample = as.character(sample),
      expression  = suppressWarnings(as.numeric(get(gene_name)))
    )]
    expr_dt <- expr_dt[is_tumor_expr(expr_sample) & !is.na(expression)]
    expr_dt[, sample_id := extract_patient_id_expr(expr_sample)]
    expr_dt[, cancer := sub("^TCGA-([^-]+)-.*$", "\\1", expr_sample)]
    expr_dt <- expr_dt[!is.na(sample_id)]

    expr_dt <- expr_dt[, .(
      expr_sample = expr_sample[1],
      expression  = expression[1],
      cancer      = cancer[1]
    ), by = sample_id]

    samples_needed <- unique(expr_dt$sample_id)
    match_dt <- merge(
      data.table(sample_id = samples_needed),
      meth_index[, .(sample_id, file)],
      by = "sample_id",
      all.x = TRUE
    )
    match_dt <- match_dt[!is.na(file)]
    if (nrow(match_dt) == 0) stop("No matched methylation files for kept samples")

    sid_vec  <- match_dt$sample_id
    file_vec <- match_dt$file
    n_files  <- length(file_vec)

    batch_id    <- ceiling(seq_len(n_files) / batch_size)
    batch_split <- split(seq_len(n_files), batch_id)
    n_batches   <- length(batch_split)

    meth_chunks <- vector("list", n_batches)
    t0 <- proc.time()[3]

    for (b in seq_along(batch_split)) {
      idx <- batch_split[[b]]
      batch_start <- proc.time()[3]

      batch_res <- mclapply(
        idx,
        function(i) {
          dtm <- read_sample_meth_fast(file_vec[i], probes_keep)
          if (!is.null(dtm) && nrow(dtm) > 0) {
            data.table(
              sample_id = sid_vec[i],
              meth_beta = median(dtm$beta, na.rm = TRUE),
              n_probes_found = nrow(dtm)
            )
          } else {
            NULL
          }
        },
        mc.cores = min(batch_size, length(idx))
      )

      meth_chunks[[b]] <- rbindlist(batch_res, fill = TRUE)

      elapsed_total    <- as.numeric(proc.time()[3] - t0)
      elapsed_batch    <- as.numeric(proc.time()[3] - batch_start)
      done_files       <- max(idx)
      avg_per_file     <- elapsed_total / done_files
      remaining        <- (n_files - done_files) * avg_per_file
      recovered_so_far <- sum(vapply(meth_chunks[seq_len(b)], nrow, integer(1)))

      cat(
        "  batch", b, "/", n_batches,
        "| files:", length(idx),
        "| done:", done_files, "/", n_files,
        "| recovered:", recovered_so_far,
        "| batch time:", fmt_sec(elapsed_batch),
        "| elapsed:", fmt_sec(elapsed_total),
        "| ETA:", fmt_sec(remaining), "\n"
      )
    }

    meth_dt <- rbindlist(meth_chunks, fill = TRUE)
    if (nrow(meth_dt) == 0) stop("No methylation values recovered")

    meth_dt <- meth_dt[, .(
      meth_beta      = median(meth_beta, na.rm = TRUE),
      n_probes_found = max(n_probes_found, na.rm = TRUE)
    ), by = sample_id]

    plot_dt <- merge(expr_dt, meth_dt, by = "sample_id", all = FALSE)
    plot_dt <- merge(
      plot_dt,
      gene_status[, .(sample_id, motif_status, alteration)],
      by = "sample_id",
      all.x = TRUE
    )

    plot_dt[is.na(motif_status), motif_status := 0L]
    plot_dt[is.na(alteration), alteration := "Unchanged"]
    plot_dt <- plot_dt[alteration %in% c("Unchanged", "Hypo altered", "Hyper altered")]

    if (nrow(plot_dt) == 0) stop("No samples left after merging")

    plot_dt[, methylation_pct := meth_beta]
    plot_dt[, log_expr := log2(expression + 1)]
    plot_dt[, alteration := factor(alteration, levels = c("Unchanged", "Hypo altered", "Hyper altered"))]
    plot_dt[, meth_bin := fifelse(methylation_pct < 50, "Unmethylated", "Methylated")]
    plot_dt[, meth_bin := factor(meth_bin, levels = c("Unmethylated", "Methylated"))]

    missing_cols <- setdiff(unique(plot_dt$cancer), names(cancer_cols))
    if (length(missing_cols)) stop("Missing colors for: ", paste(missing_cols, collapse = ", "))

    plot_dt[, cancer := factor(cancer, levels = intersect(unique(cols_dt$cancer), unique(plot_dt$cancer)))]

    cor_res <- safe_cor(plot_dt$methylation_pct, plot_dt$log_expr, method = "pearson")
    p_bin <- tryCatch(
      wilcox.test(log_expr ~ meth_bin, data = plot_dt)$p.value,
      error = function(e) NA_real_
    )

    shape_values <- c(
      "Unchanged"     = 16,
      "Hypo altered"  = 17,
      "Hyper altered" = 15
    )

    out_pdf <- paste0(out_dir, "/", gene_name, "_expr_vs_methylation_3pages.pdf")
    out_tsv <- paste0(out_dir, "/", gene_name, "_expr_vs_methylation_3pages.tsv.gz")

    fwrite(plot_dt, out_tsv, sep = "\t")

    p1 <- ggplot() +
      geom_point(
        data = plot_dt[alteration == "Unchanged"],
        aes(x = methylation_pct, y = log_expr, color = cancer, shape = alteration),
        size = 1.2, alpha = 0.10
      ) +
      geom_point(
        data = plot_dt[alteration != "Unchanged"],
        aes(x = methylation_pct, y = log_expr, color = cancer, shape = alteration),
        size = 2.2, alpha = 0.90
      ) +
      geom_smooth(
        data = plot_dt,
        mapping = aes(x = methylation_pct, y = log_expr),
        method = "lm",
        se = FALSE,
        inherit.aes = FALSE,
        color = "black",
        linewidth = 0.8
      ) +
      scale_color_manual(values = cancer_cols, drop = FALSE) +
      scale_shape_manual(values = shape_values, drop = FALSE) +
      labs(
        title = paste0(gene_name, ": pooled across all cancer types"),
        subtitle = paste0(
          "log2(expression + 1) | faded unchanged, highlighted altered | Pearson r = ",
          ifelse(is.finite(cor_res$r), round(cor_res$r, 3), "NA"),
          " ; p = ",
          ifelse(is.finite(cor_res$p), formatC(cor_res$p, format = "e", digits = 2), "NA")
        ),
        x = "Median methylation across proximal gene-linked probes (%)",
        y = "log2(expression + 1)"
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank()
      )

    p2 <- ggplot() +
      geom_point(
        data = plot_dt[alteration == "Unchanged"],
        aes(x = methylation_pct, y = log_expr, color = cancer, shape = alteration),
        size = 0.9, alpha = 0.08
      ) +
      geom_point(
        data = plot_dt[alteration != "Unchanged"],
        aes(x = methylation_pct, y = log_expr, color = cancer, shape = alteration),
        size = 1.6, alpha = 0.85
      ) +
      geom_smooth(
        data = plot_dt,
        mapping = aes(x = methylation_pct, y = log_expr, group = 1),
        method = "lm",
        se = FALSE,
        color = "black",
        linewidth = 0.5,
        inherit.aes = FALSE
      ) +
      scale_color_manual(values = cancer_cols, drop = FALSE) +
      scale_shape_manual(values = shape_values, drop = FALSE) +
      facet_wrap(~ cancer, scales = "free_y") +
      labs(
        title = paste0(gene_name, ": faceted by cancer type"),
        subtitle = "Full data with unchanged faded and altered highlighted",
        x = "Median methylation across proximal gene-linked probes (%)",
        y = "log2(expression + 1)"
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        legend.position = "bottom"
      )

    p3 <- ggplot(
      plot_dt,
      aes(x = meth_bin, y = log_expr, fill = meth_bin)
    ) +
      geom_violin(trim = TRUE, scale = "width", alpha = 0.8) +
      geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white") +
      geom_jitter(
        data = plot_dt[alteration == "Unchanged"],
        aes(shape = alteration),
        width = 0.15, alpha = 0.08, size = 0.8, color = "black"
      ) +
      geom_jitter(
        data = plot_dt[alteration != "Unchanged"],
        aes(shape = alteration),
        width = 0.15, alpha = 0.35, size = 1.2, color = "black"
      ) +
      scale_shape_manual(values = shape_values, drop = FALSE) +
      labs(
        title = paste0(gene_name, ": expression by binary methylation bin"),
        subtitle = paste0(
          "< 50% = Unmethylated ; >= 50% = Methylated",
          if (is.finite(p_bin)) paste0(" | Wilcoxon p = ", formatC(p_bin, format = "e", digits = 2)) else ""
        ),
        x = "Methylation bin",
        y = "log2(expression + 1)"
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank()
      ) +
      guides(fill = "none")

    pdf(out_pdf, width = 12, height = 8)
    print(p1)
    print(p2)
    print(p3)
    dev.off()

    data.table(
      gene = gene_name,
      status = "OK",
      n_samples = nrow(plot_dt),
      n_unchanged = plot_dt[alteration == "Unchanged", .N],
      n_hypo = plot_dt[alteration == "Hypo altered", .N],
      n_hyper = plot_dt[alteration == "Hyper altered", .N],
      cor_r = cor_res$r,
      cor_p = cor_res$p
    )

  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    data.table(
      gene = gene_name,
      status = paste("ERROR:", conditionMessage(e)),
      n_samples = NA_integer_,
      n_unchanged = NA_integer_,
      n_hypo = NA_integer_,
      n_hyper = NA_integer_,
      cor_r = NA_real_,
      cor_p = NA_real_
    )
  })
}
cat("\nDone.\n")
EOF

# Summary of above codes:
# These scripts form a small pipeline to study the link between **motif methylation alterations** and **gene expression**. The first script builds the alteration calls by comparing, for each tumor sample, the methylation level of probes overlapping NRF1/BANP motifs to the **median methylation of normal samples from the same cancer type**, computes a delta for each probe, keeps the strongest probe signal per motif, and classifies the motif in each sample as **hyper altered** (`+1`), **hypo altered** (`-1`), or **unchanged** (`0`) using the chosen threshold. The second script then uses this alteration information together with the tumor expression matrix and per-sample methylation files to perform a **per-gene analysis restricted to proximal motif-linked events**: for each gene linked to altered proximal motifs, it aggregates the sample alteration status, extracts the relevant proximal methylation probes, computes a sample-level median methylation value, merges this with gene expression, and generates plots and summary tables showing how methylation and motif alteration status relate to expression across tumors and cancer types.


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
pdf("./results/snv/snv_per_sample_hist.pdf", width=10, height=8)
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
ggsave("./results/snv/snv_min_max_dumbbell.pdf",plot = p, width = 7, height = 9, dpi = 300, bg = "white")'


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
data <- read.table("./snv/snv_counts_per_cancer_type/snv_unique_per_cancer.tsv",header = FALSE, sep = "\t", stringsAsFactors = FALSE)
cancer_names <- trimws(data[,1])
snv_counts   <- as.numeric(data[,2])

# Sort from biggest to smallest
ord <- order(snv_counts, decreasing = TRUE)
cancer_names <- cancer_names[ord]
snv_counts   <- snv_counts[ord]

# Reverse order for horizontal barplot (largest at top)
cancer_names <- rev(cancer_names)
snv_counts   <- rev(snv_counts)

snv_millions <- snv_counts / 1e6

# Read defined colors
color_df <- read.table("./results/multi_omics/cancer_color_order_with_defined_colours.tsv",header = FALSE, sep = "\t", stringsAsFactors = FALSE,comment.char = "", quote = "")
colnames(color_df) <- c("Cancer_full", "Color")
color_df$Cancer <- sub("^TCGA-", "", trimws(color_df$Cancer_full))
color_df$Color  <- trimws(gsub("\r", "", color_df$Color))

# Build named palette
palette <- color_df$Color
names(palette) <- color_df$Cancer

# Match colors to cancers in plot
bar_colors <- palette[cancer_names]
pdf("./results/snv/snv_per_cancer_type_barplot.pdf", width = 12, height = 9, bg = "white")
nice_max <- max(pretty(snv_millions))
bp <- barplot(
  snv_millions,
  names.arg = cancer_names,
  horiz = TRUE,
  las = 1,
  col = unname(bar_colors),
  border = "black",
  main = "SNVs per Cancer Type",
  cex.main = 3,
  cex.lab  = 1.5,
  xlab = "Total SNVs (millions)",
  xlim = c(0, nice_max * 1.08)
)
text(
  x = snv_millions,
  y = bp,
  labels = format(snv_counts, big.mark = ",", scientific = FALSE),
  pos = 4,
  cex = 0.8,
  xpd = TRUE
)

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
ggsave("./results/snv/stacked_mutation_types_log10.pdf",p,width = 12,height = 7,dpi = 300,bg = "white")'

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
pdf("./results/snv/variant_type_percentages_boxplot_per_sample.pdf",width = 16, height = 9)

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
# 2) FILTER STRUCTURAL VARIANTS and chrX and Y and M 
###############################################################################

# Since structural variants (indels >= 50 bp) are less relevant for TF binding site analysis, we filter them out.
# filtering VCF files to remove structural variants (indels >= 50 bp) and ones at chrX, Y and M 
mkdir -p ./snv/snv_filtered_without_structural_variants/

for vcf in ./snv/SNV_TCGA-*.vcf.gz; do
  out=./snv/snv_filtered_without_structural_variants/$(basename "$vcf")

  zcat "$vcf" | awk '
    BEGIN{FS="\t"; OFS="\t"}
    /^#/ {print; next}
    {
      chr=$1
      gsub(/^chr/,"",chr)                 # allow chrX or X
      if (chr=="X" || chr=="Y" || chr=="M" || chr=="MT") next

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

#######################################################################################################
# SAMPLES WITH SNV PER CANCER TYPE --> ./snv/snv_sample_counts_per_cancer/samples_per_cancer.tsv
#######################################################################################################

VCF_DIR="./snv/snv_filtered_without_structural_variants"
OUT_DIR="./snv/snv_sample_counts_per_cancer"
OUT_FILE="${OUT_DIR}/samples_per_cancer.tsv"

mkdir -p "$OUT_DIR"

# Extract cancer type from filename SNV_TCGA-<CANCER>-..., count occurrences
ls "${VCF_DIR}"/SNV_TCGA-*.vcf.gz 2>/dev/null \
  | sed -E 's|.*/SNV_(TCGA-[A-Z0-9]+)-.*|\1|' \
  | sort \
  | uniq -c \
  | awk -v OFS="\t" '{print $2, $1}' \
  > "$OUT_FILE"

echo "[DONE] Wrote: $OUT_FILE"

# barplot 
Rscript -e '
data <- read.table("./snv/snv_sample_counts_per_cancer/samples_per_cancer.tsv",header=FALSE, sep="\t", stringsAsFactors=FALSE)
colnames(data) <- c("cancer","sample_count")
data <- data[order(data$sample_count, decreasing=FALSE), ]

cols <- read.table("./results/multi_omics/cancer_color_order_with_defined_colours.tsv",header=FALSE, sep="\t", stringsAsFactors=FALSE, quote="", comment.char="")
colnames(cols) <- c("cancer","color")

palette <- setNames(cols$color, cols$cancer)
bar_colors <- unname(palette[data$cancer])

# validate colors (will error if any are not valid)
bad <- which(!is.na(bar_colors) & tryCatch(is.na(col2rgb(bar_colors)), error=function(e) TRUE))
if (length(bad) > 0) {
  cat("[ERROR] These colors are invalid:\n")
  print(unique(bar_colors[bad]))
  quit(status=1)
}

pdf("./results/snv/samples_with_snvs_per_cancer_barplot.pdf", width=9, height=9, bg="white")
nice_max <- max(pretty(data$sample_count))
labels <- sub("^TCGA-", "", data$cancer)

barplot(data$sample_count, names.arg=labels, horiz=TRUE, las=1,
        col=bar_colors, border="black",
        main="Samples per Cancer Type",
        cex.main=3, cex.lab=1.5, xlab="Number of samples",
        xlim=c(0, nice_max))

dev.off()
'


###############################################################################
# Frequency of variants found within each cancer type after filtering structural variants
###############################################################################

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

# Plot histogram of distances to TSS for all samples in one pdf 
Rscript -e '
dist_files <- list.files(
  "./snv/dist_to_tss/",
  pattern = "_dist_tss.txt$",
  full.names = TRUE
)

out_pdf <- "./results/snv/dist_tss_snv_histograms.pdf"

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

# Create a unique list of SNVs across all cancer types (after filtering structural variants)
zcat ./snv/snv_filtered_without_structural_variants/*.vcf.gz \
| awk 'BEGIN{OFS="\t"} !/^#/ {
    snv_id = $1 ":" $2 ":" $4 ":" $5
    print $1, $2-1, $2, snv_id
}' \
| sort -k1,1 -k2,2n -k3,3n -k4,4 \
| uniq \
> ./snv/all_unique_variants_across_cancers.bed

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

ggsave("./results/snv/SNV_in_vs_out_ATAC_barplot.pdf",p, width = 7, height = 5, bg = "white")
'
###############################################################################
# 6) SNV ∩ METHYLATION
###############################################################################
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
ggsave("./results/snv/SNV_in_vs_out_Methylation_barplot.pdf",p, width = 7, height = 5, bg = "white")
'

######################################################
# Overlap filtered SNV in TF 2mm motifs and ATAC peaks
###################################################### 
# -------------------------
# FOLDERS
# -------------------------
mkdir -p ./snv/overlaps/snv_in_motifs2mm
mkdir -p ./snv/overlaps/snv_in_motifs_in_peaks2mm

mkdir -p ./motifs/overlaps/motif2mm_snv_overlaps
mkdir -p ./motifs/overlaps/motifs2mm_in_peaks_snv

# -------------------------
# 1) SNVs that are in 2mm motifs  (SNV output)
# -------------------------
# NRF1 
bedtools intersect -u \
  -a ./snv/all_unique_variants_across_cancers.bed \
  -b <(zcat ./motifs/NRF1_mm0to2_noXY.bed.gz) \
  > ./snv/overlaps/snv_in_motifs2mm/NRF1_mm0to2_SNVs_in_motifs.bed

bedtools intersect -u \
  -a ./snv/all_unique_variants_across_cancers.bed \
  -b <(zcat ./motifs/NRF1_mm0to2_noCGmm_noXY.bed.gz) \
  > ./snv/overlaps/snv_in_motifs2mm/NRF1_mm0to2_noCGmm_SNVs_in_motifs.bed

# BANP
bedtools intersect -u \
  -a ./snv/all_unique_variants_across_cancers.bed \
  -b <(zcat ./motifs/BANP_mm0to2_noXY.bed.gz) \
  > ./snv/overlaps/snv_in_motifs2mm/BANP_mm0to2_SNVs_in_motifs.bed

bedtools intersect -u \
  -a ./snv/all_unique_variants_across_cancers.bed \
  -b <(zcat ./motifs/BANP_mm0to2_noCGmm_noXY.bed.gz) \
  > ./snv/overlaps/snv_in_motifs2mm/BANP_mm0to2_noCGmm_SNVs_in_motifs.bed


# -------------------------
# 2) Motifs that contain SNVs (motif output)
# -------------------------
# NRF1 
bedtools intersect -u \
  -a <(zcat ./motifs/NRF1_mm0to2_noXY.bed.gz) \
  -b ./snv/all_unique_variants_across_cancers.bed \
  > ./motifs/overlaps/motif2mm_snv_overlaps/NRF1_mm0to2_motifs_with_SNVs.bed
bedtools intersect -u \
  -a <(zcat ./motifs/NRF1_mm0to2_noCGmm_noXY.bed.gz) \
  -b ./snv/all_unique_variants_across_cancers.bed \
  > ./motifs/overlaps/motif2mm_snv_overlaps/NRF1_mm0to2_noCGmm_motifs_with_SNVs.bed

# BANP
bedtools intersect -u \
  -a <(zcat ./motifs/BANP_mm0to2_noXY.bed.gz) \
  -b ./snv/all_unique_variants_across_cancers.bed \
  > ./motifs/overlaps/motif2mm_snv_overlaps/BANP_mm0to2_motifs_with_SNVs.bed

bedtools intersect -u \
  -a <(zcat ./motifs/BANP_mm0to2_noCGmm_noXY.bed.gz) \
  -b ./snv/all_unique_variants_across_cancers.bed \
  > ./motifs/overlaps/motif2mm_snv_overlaps/BANP_mm0to2_noCGmm_motifs_with_SNVs.bed


# -------------------------
# Motifs that contain SNVs snvidxmotif output
# -------------------------
# NRF1 
bedtools intersect -wa -wb \
  -a <(zcat ./motifs/NRF1_mm0to2_noXY.bed.gz) \
  -b ./snv/all_unique_variants_across_cancers.bed \
  > ./snv/overlaps/snv_in_motifs2mm/NRF1_mm0to2_SNVs_in_motifs_SNVidxmotif.bed
bedtools intersect -wa -wb\
  -a <(zcat ./motifs/NRF1_mm0to2_noCGmm_noXY.bed.gz) \
  -b ./snv/all_unique_variants_across_cancers.bed \
  > ./snv/overlaps/snv_in_motifs2mm/NRF1_mm0to2_noCGmm_SNVs_in_motifs_SNVidxmotif.bed
# BANP
bedtools intersect -wa -wb \
  -a <(zcat ./motifs/BANP_mm0to2_noXY.bed.gz) \
  -b ./snv/all_unique_variants_across_cancers.bed \
  > ./snv/overlaps/snv_in_motifs2mm/BANP_mm0to2_SNVs_in_motifs_SNVidxmotif.bed

bedtools intersect -wa -wb \
  -a <(zcat ./motifs/BANP_mm0to2_noCGmm_noXY.bed.gz) \
  -b ./snv/all_unique_variants_across_cancers.bed \
  > ./snv/overlaps/snv_in_motifs2mm/BANP_mm0to2_noCGmm_SNVs_in_motifs_SNVidxmotif.bed

# -------------------------
# 3) SNVs in motifs-in-peaks (SNV output)
# -------------------------
# NRF1 
bedtools intersect -u \
  -a ./snv/all_unique_variants_across_cancers.bed \
  -b ./motifs/overlaps/motif2mm_peak_overlaps/NRF1_mm0to2_peak_overlaps.bed \
  > ./snv/overlaps/snv_in_motifs_in_peaks2mm/NRF1_mm0to2_SNVs_in_motifs_in_peaks.bed
bedtools intersect -u \
  -a ./snv/all_unique_variants_across_cancers.bed \
  -b ./motifs/overlaps/motif2mm_peak_overlaps/NRF1_mm0to2_noCGmm_peak_overlaps.bed \
  > ./snv/overlaps/snv_in_motifs_in_peaks2mm/NRF1_mm0to2_noCGmm_SNVs_in_motifs_in_peaks.bed

# BANP 
bedtools intersect -u \
  -a ./snv/all_unique_variants_across_cancers.bed \
  -b ./motifs/overlaps/motif2mm_peak_overlaps/BANP_mm0to2_peak_overlaps.bed \
  > ./snv/overlaps/snv_in_motifs_in_peaks2mm/BANP_mm0to2_SNVs_in_motifs_in_peaks.bed
bedtools intersect -u \
  -a ./snv/all_unique_variants_across_cancers.bed \
  -b ./motifs/overlaps/motif2mm_peak_overlaps/BANP_mm0to2_noCGmm_peak_overlaps.bed \
  > ./snv/overlaps/snv_in_motifs_in_peaks2mm/BANP_mm0to2_noCGmm_SNVs_in_motifs_in_peaks.bed


# -------------------------
# 4) Motifs-in-peaks that contain SNVs (motif output)
# -------------------------
# NRF1 
bedtools intersect -u \
  -a ./motifs/overlaps/motif2mm_peak_overlaps/NRF1_mm0to2_peak_overlaps.bed \
  -b ./snv/all_unique_variants_across_cancers.bed \
  > ./motifs/overlaps/motifs2mm_in_peaks_snv/NRF1_mm0to2_motifs_in_peaks_with_SNVs.bed

bedtools intersect -u \
  -a ./motifs/overlaps/motif2mm_peak_overlaps/NRF1_mm0to2_noCGmm_peak_overlaps.bed \
  -b ./snv/all_unique_variants_across_cancers.bed \
  > ./motifs/overlaps/motifs2mm_in_peaks_snv/NRF1_mm0to2_noCGmm_motifs_in_peaks_with_SNVs.bed

# BANP
bedtools intersect -u \
  -a ./motifs/overlaps/motif2mm_peak_overlaps/BANP_mm0to2_peak_overlaps.bed \
  -b ./snv/all_unique_variants_across_cancers.bed \
  > ./motifs/overlaps/motifs2mm_in_peaks_snv/BANP_mm0to2_motifs_in_peaks_with_SNVs.bed

bedtools intersect -u \
  -a ./motifs/overlaps/motif2mm_peak_overlaps/BANP_mm0to2_noCGmm_peak_overlaps.bed \
  -b ./snv/all_unique_variants_across_cancers.bed \
  > ./motifs/overlaps/motifs2mm_in_peaks_snv/BANP_mm0to2_noCGmm_motifs_in_peaks_with_SNVs.bed


###############################################################################
# 8) SNV ∩ PEAKS ∩ METHYLATION
###############################################################################

# Counting cancer types with both peaks and SNV data
mkdir -p ./results/multi_omics

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

ggsave("./results/multi_omics/cancer_types_ATAC_SNV_Methylation_heatmap.pdf",p, width = 10, height = 10, dpi = 300, bg = "white")
'
# clean up intermediate files
rm ./methylation/cancer_types_with_methylation.txt
rm ./peaks/cancer_types_with_peaks.txt

###############################################################################
# 9) SNV HEATMAP ACROSS CANCER TYPES for 6 mer motifs
###############################################################################
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
# count number of SNVs that overlap with NRF1 and BANP motifs
cat ./motifs/overlaps/motifs_snv_overlap/NRF1_filtered_SNVs_in_motifskept.bed  | wc -l > ./snv/overlaps/snv_motifs_overlap/NRF1_intersected_SNVs_count_kept.tsv 
cat ./motifs/overlaps/motifs_snv_overlap/BANP_filtered_SNVs_in_motifskept.bed  | wc -l > ./snv/overlaps/snv_motifs_overlap/BANP_intersected_SNVs_count_kept.tsv
# Plot barplot of SNV counts in NRF1 and BANP motifs
Rscript -e '
library(ggplot2)

# Motif names
motifs <- c("NRF1", "BANP")

# Read the numeric counts from the .tsv files (SNVs in motifs)
snv_counts <- c(
  as.numeric(readLines("./snv/overlaps/snv_motifs_overlap/NRF1_intersected_SNVs_count_kept.tsv")),
  as.numeric(readLines("./snv/overlaps/snv_motifs_overlap/BANP_intersected_SNVs_count_kept.tsv")))

# TOTAL motif counts
total_motifs <- c(
  843539,  
  86505)

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

ggsave("./results/snv/SNV_counts_in_motifs_NRF1_BANP_kept.pdf",plot = p, width = 10, height = 12, dpi = 300, bg = "white")
'

###############################################################################
# SNV HEATMAP ACROSS CANCER TYPES for 2mm motifs
###############################################################################

###############################################################################
# CONFIG
###############################################################################
vcf_dir="./snv/snv_filtered_without_structural_variants"
per_cancer_dir="./snv/snv_filtered_unique_per_cancer"
out_dir="./snv"
results_dir="./results/snv"

mkdir -p "$per_cancer_dir" "$results_dir"

# Cancer list / order (consistent across all matrices)
cancers=$(ls "${per_cancer_dir}"/*_unique_SNVs.txt.gz | sed 's#.*/##; s/_unique_SNVs\.txt\.gz##' | sort)
echo "[INFO] Cancers: $(echo "$cancers" | wc -l)"

###############################################################################
# 2) Build 4 SNV-ID lists (from your 2mm overlap beds)
###############################################################################
echo "[RUN] Building SNV ID lists from 2mm motif-overlap BEDs..."

# Takes as input these BED files:
# ./snv/overlaps/snv_in_motifs2mm/NRF1_mm0to2_SNVs_in_motifs.bed
# ./snv/overlaps/snv_in_motifs2mm/BANP_mm0to2_SNVs_in_motifs.bed
# ./snv/overlaps/snv_in_motifs2mm/NRF1_mm0to2_noCGmm_SNVs_in_motifs.bed
# ./snv/overlaps/snv_in_motifs2mm/BANP_mm0to2_noCGmm_SNVs_in_motifs.bed

cut -f4 ./snv/overlaps/snv_in_motifs2mm/NRF1_mm0to2_SNVs_in_motifs.bed | sort -u \
  > "./snv/overlaps/snv_in_motifs2mm/NRF1_mm0to2_SNV_ids.txt"

cut -f4 ./snv/overlaps/snv_in_motifs2mm/BANP_mm0to2_SNVs_in_motifs.bed | sort -u \
  > "./snv/overlaps/snv_in_motifs2mm/BANP_mm0to2_SNV_ids.txt"

cut -f4 ./snv/overlaps/snv_in_motifs2mm/NRF1_mm0to2_noCGmm_SNVs_in_motifs.bed | sort -u \
  > "./snv/overlaps/snv_in_motifs2mm/NRF1_mm0to2_noCGmm_SNV_ids.txt"

cut -f4 ./snv/overlaps/snv_in_motifs2mm/BANP_mm0to2_noCGmm_SNVs_in_motifs.bed | sort -u \
  > "./snv/overlaps/snv_in_motifs2mm/BANP_mm0to2_noCGmm_SNV_ids.txt"

#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# INPUTS 
###############################################################################
# These SNV-id lists must be one SNV per line, like: chr1:123:A:T

snv_ids_NRF1_mm0to2="./snv/overlaps/snv_in_motifs2mm/NRF1_mm0to2_SNV_ids.txt"
snv_ids_BANP_mm0to2="./snv/overlaps/snv_in_motifs2mm/BANP_mm0to2_SNV_ids.txt"
snv_ids_NRF1_mm0to2_noCGmm="./snv/overlaps/snv_in_motifs2mm/NRF1_mm0to2_noCGmm_SNV_ids.txt"
snv_ids_BANP_mm0to2_noCGmm="./snv/overlaps/snv_in_motifs2mm/BANP_mm0to2_noCGmm_SNV_ids.txt"

per_cancer_dir="./snv/snv_filtered_unique_per_cancer"   # contains *_unique_SNVs.txt.gz
mkdir -p "$per_cancer_dir"

###############################################################################
# SETTINGS (keep your original approach)
###############################################################################
threads=5
max_jobs=2

###############################################################################
# 1) Generic function to build columns (same as yours)
###############################################################################
build_column() {
    local snv_list="$1"      # list of SNV ids (one per line)
    local cancer_file="$2"   # e.g. ACC_unique_SNVs.txt.gz
    local out_file="$3"

    local cancer
    cancer=$(basename "$cancer_file" _unique_SNVs.txt.gz)
    echo "  → building column for $cancer"

    awk -v OFS='\t' '
        NR==FNR { has[$1]=1; next }     # FIRST = cancer SNVs
        {
            val = ($1 in has ? 1 : 0);
            print $1, val;
        }
    ' <(zcat "$cancer_file") "$snv_list" > "$out_file"
}

###############################################################################
# 2) Build all columns in parallel for a given tag + SNV list
###############################################################################
build_all_columns() {
    local snv_list="$1"
    local tag="$2"
    local col_dir="./snv/snv_presence_columns_${tag}"

    mkdir -p "$col_dir"

    echo "[RUN] Building presence columns for: $tag"
    job_count=0

    for f in "${per_cancer_dir}"/*_unique_SNVs.txt.gz; do
        cancer=$(basename "$f" _unique_SNVs.txt.gz)
        out="${col_dir}/${cancer}.col"

        # Optional CPU gate (same idea as yours, keep if you want)
        CPU=$(top -b -n 2 | awk '($1~"%Cpu"){cpu=$2}END{print int(48*(100-cpu)/100)}')
        while [ "$CPU" -le "$threads" ]; do
            echo "[$cancer] attente CPU : seulement $CPU CPU libres (min = $threads)"
            sleep 60
            CPU=$(top -b -n 2 | awk '($1~"%Cpu"){cpu=$2}END{print int(48*(100-cpu)/100)}')
        done

        echo "→ Lancement job pour $cancer (CPU libres = $CPU)"
        build_column "$snv_list" "$f" "$out" &

        job_count=$((job_count + 1))
        if [ "$job_count" -ge "$max_jobs" ]; then
            wait
            job_count=0
        fi
    done

    wait
    echo "[DONE] All jobs for ${tag} columns are done."
}

###############################################################################
# 3) Merge columns into a matrix (exactly your working merge pattern)
###############################################################################
merge_columns() {
    local tag="$1"
    local col_dir="./snv/snv_presence_columns_${tag}"
    local out="./snv/snv_presence_matrix_${tag}.tsv"

    local first_file="${col_dir}/ACC.col"
    if [ ! -f "$first_file" ]; then
        # fallback to first .col found if ACC is not present
        first_file=$(ls "${col_dir}"/*.col | head -n 1)
    fi

    echo "[RUN] Merging columns for: $tag"
    echo "      first_file = $first_file"

    tmp=./snv/tmp.matrix
    snv_ids=./snv/snv_ids.tmp
    header=./snv/header.tmp

    # 1) SNV IDs reference (from first file)
    awk '{print $1}' "$first_file" > "$snv_ids"

    # 2) Initialize matrix with SNV + first column values
    paste "$snv_ids" <(awk '{print $2}' "$first_file") > "$tmp"

    # 3) Append other cancers (skip first file)
    for f in "${col_dir}"/*.col; do
      [ "$f" = "$first_file" ] && continue
      base=$(basename "$f" .col)
      echo "→ Adding $base"

      # Check SNV IDs match exactly (line-by-line)
      if ! diff -q "$snv_ids" <(awk '{print $1}' "$f") >/dev/null; then
        echo "ERROR: SNV list/order mismatch in $f (must match $first_file)" >&2
        exit 1
      fi

      paste "$tmp" <(awk '{print $2}' "$f") > "${tmp}.2"
      mv "${tmp}.2" "$tmp"
    done

    # 4) Header
    echo -ne "SNV" > "$header"
    for f in "${col_dir}"/*.col; do
      base=$(basename "$f" .col)
      echo -ne "\t$base" >> "$header"
    done
    echo >> "$header"

    cat "$header" "$tmp" > "$out"

    rm -f "$tmp" "$snv_ids" "$header"
    echo "[DONE] Wrote: $out"
}

###############################################################################
# 4) RUN for the 4 mm0to2 cases
###############################################################################
# NRF1 mm0to2
build_all_columns "$snv_ids_NRF1_mm0to2" "NRF1_mm0to2"
merge_columns "NRF1_mm0to2"

# BANP mm0to2
build_all_columns "$snv_ids_BANP_mm0to2" "BANP_mm0to2"
merge_columns "BANP_mm0to2"

# NRF1 mm0to2_noCGmm
build_all_columns "$snv_ids_NRF1_mm0to2_noCGmm" "NRF1_mm0to2_noCGmm"
merge_columns "NRF1_mm0to2_noCGmm"

# BANP mm0to2_noCGmm
build_all_columns "$snv_ids_BANP_mm0to2_noCGmm" "BANP_mm0to2_noCGmm"
merge_columns "BANP_mm0to2_noCGmm"

echo "[ALL DONE]"

###############################################################################
# 4) Heatmaps for the two pairs
###############################################################################
Rscript -e '
library(pheatmap)

plot_pair <- function(nrf1_mat, banp_mat, out_pdf_log, out_pdf_pct, tag){

  # ---------------- LOG10 shared counts ----------------
  pdf(out_pdf_log, width=10, height=8, bg="white")

  # NRF1
  data <- read.table(nrf1_mat, header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
  mat <- as.matrix(data[,-1]); rownames(mat) <- data$SNV
  shared <- t(mat) %*% mat
  shared_log <- log10(shared + 1)
  diag(shared_log) <- NA

  min_val <- min(shared_log, na.rm=TRUE)
  max_val <- max(shared_log, na.rm=TRUE)

  pheatmap(shared_log,
           main = paste0("Number of shared NRF1-motif variants between cancer types (log10)\\n", tag),
           cluster_cols=TRUE, cluster_rows=TRUE,
           show_rownames=TRUE, show_colnames=TRUE,
           color=colorRampPalette(c("white","steelblue","navy"))(200),
           breaks=seq(min_val, max_val, length.out=201),
           na_col="grey90")

  # BANP
  data <- read.table(banp_mat, header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
  mat <- as.matrix(data[,-1]); rownames(mat) <- data$SNV
  shared <- t(mat) %*% mat
  shared_log <- log10(shared + 1)
  diag(shared_log) <- NA

  min_val <- min(shared_log, na.rm=TRUE)
  max_val <- max(shared_log, na.rm=TRUE)

  pheatmap(shared_log,
           main = paste0("Number of shared BANP-motif variants between cancer types (log10)\\n", tag),
           cluster_cols=TRUE, cluster_rows=TRUE,
           show_rownames=TRUE, show_colnames=TRUE,
           color=colorRampPalette(c("white","steelblue","navy"))(200),
           breaks=seq(min_val, max_val, length.out=201),
           na_col="grey90")

  dev.off()

  # ---------------- % shared (row-normalized by diagonal) ----------------
  pdf(out_pdf_pct, width=10, height=8, bg="white")

  # NRF1
  data <- read.table(nrf1_mat, header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
  cancer_cols  <- colnames(data)[-1]
  cancer_order <- sort(cancer_cols)

  mat <- as.matrix(data[, cancer_order])
  rownames(mat) <- data$SNV

  shared <- t(mat) %*% mat
  totals <- diag(shared)
  pct_shared <- sweep(shared, 1, totals, FUN="/") * 100
  diag(pct_shared) <- NA

  min_val <- min(pct_shared, na.rm=TRUE)
  max_val <- max(pct_shared, na.rm=TRUE)

  pheatmap(pct_shared,
           main = paste0("Percentage of shared NRF1-motif variants between cancer types\\n", tag),
           cluster_cols=TRUE, cluster_rows=TRUE,
           show_rownames=TRUE, show_colnames=TRUE,
           color=colorRampPalette(c("white","steelblue","navy"))(200),
           breaks=seq(min_val, max_val, length.out=201),
           na_col="grey90")

  # BANP
  data <- read.table(banp_mat, header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
  cancer_cols  <- colnames(data)[-1]
  cancer_order <- sort(cancer_cols)

  mat <- as.matrix(data[, cancer_order])
  rownames(mat) <- data$SNV

  shared <- t(mat) %*% mat
  totals <- diag(shared)
  pct_shared <- sweep(shared, 1, totals, FUN="/") * 100
  diag(pct_shared) <- NA

  min_val <- min(pct_shared, na.rm=TRUE)
  max_val <- max(pct_shared, na.rm=TRUE)

  pheatmap(pct_shared,
           main = paste0("Percentage of shared BANP-motif variants between cancer types\\n", tag),
           cluster_cols=TRUE, cluster_rows=TRUE,
           show_rownames=TRUE, show_colnames=TRUE,
           color=colorRampPalette(c("white","steelblue","navy"))(200),
           breaks=seq(min_val, max_val, length.out=201),
           na_col="grey90")

  dev.off()
}

# Pair A: mm0to2
plot_pair(
  "./snv/snv_presence_matrix_NRF1_mm0to2.tsv",
  "./snv/snv_presence_matrix_BANP_mm0to2.tsv",
  "./results/snv/SNV_shared_NRF1_BANP_mm0to2_log10_heatmap.pdf",
  "./results/snv/SNV_shared_NRF1_BANP_mm0to2_percentage_heatmap.pdf",
  "mm0to2"
)

# Pair B: mm0to2_noCGmm
plot_pair(
  "./snv/snv_presence_matrix_NRF1_mm0to2_noCGmm.tsv",
  "./snv/snv_presence_matrix_BANP_mm0to2_noCGmm.tsv",
  "./results/snv/SNV_shared_NRF1_BANP_mm0to2_noCGmm_log10_heatmap.pdf",
  "./results/snv/SNV_shared_NRF1_BANP_mm0to2_noCGmm_percentage_heatmap.pdf",
  "mm0to2_noCGmm"
)
'
echo "[DONE] Heatmaps written to ${results_dir}"

###############################################################################
# Number of SNVs in NRF1 and BANP motifs across all cancer types 
###############################################################################
#!/usr/bin/env bash
set -euo pipefail

IN_DIR="./snv/overlaps/snv_in_motifs2mm"
OUT="./results/snv/SNV_counts_in_motifs_NRF1_BANP_mm0to2_noCGmm.tsv"
mkdir -p "$(dirname "$OUT")"

NRF1_BED="$IN_DIR/NRF1_mm0to2_noCGmm_SNVs_in_motifs.bed"
BANP_BED="$IN_DIR/BANP_mm0to2_noCGmm_SNVs_in_motifs.bed"

n_nrf1=$(wc -l < "$NRF1_BED" | tr -d " ")
n_banp=$(wc -l < "$BANP_BED" | tr -d " ")

echo -e "NRF1_mm0to2_noCGmm\tBANP_mm0to2_noCGmm" > "$OUoT"
echo -e "${n_nrf1}\t${n_banp}" >> "$OUT"

echo "[DONE] Wrote $OUT"

# barplot snv counts in NRF1 and BANP motifs across all cancer types
Rscript -e '
library(ggplot2)
# Read counts
data <- read.table("./results/snv/SNV_counts_in_motifs_NRF1_BANP_mm0to2_noCGmm.tsv", header=TRUE, stringsAsFactors=FALSE)
# Transform to long format
df_long <- data.frame(
  Motif = c("NRF1_mm0to2_noCGmm", "BANP_mm0to2_noCGmm"),
  SNV_Count = c(data$NRF1_mm0to2_noCGmm, data$BANP_mm0to2_noCGmm)
)
# Plot
p <- ggplot(df_long, aes(x=Motif, y=SNV_Count, fill=Motif)) +
  geom_bar(stat="identity", width=0.6) +
  theme_minimal() +
  labs(title="Number of SNVs in NRF1 and BANP motifs (mm0to2_noCGmm)", y="Number of SNVs", x="Motif") +
  theme(plot.title = element_text(hjust=0.5, size=12), axis.text=element_text(size=10), axis.title=element_text(size=11), legend.position="none")
print(p)
ggsave("./results/snv/SNV_counts_in_motifs_NRF1_BANP_mm0to2_noCGmm_barplot.pdf", plot=p, width=8, height=6, dpi=300, bg="white")
' 
cat ./results/snv/SNV_counts_in_motifs_NRF1_BANP_mm0to2_noCGmm.tsv
# NRF1_mm0to2_noCGmm      BANP_mm0to2_noCGmm
# 62436   20329


###############################################################################
# Number of SNVs in NRF1 and BANP motifs and in peaks across all cancer types 
###############################################################################
#!/usr/bin/env bash
set -euo pipefail

IN_DIR="./snv/overlaps/snv_in_motifs_in_peaks2mm"
OUT="./results/snv/SNV_counts_in_motifsnpeaks_NRF1_BANP_mm0to2_noCGmm.tsv"
mkdir -p "$(dirname "$OUT")"

NRF1_BED="$IN_DIR/NRF1_mm0to2_noCGmm_SNVs_in_motifs_in_peaks.bed"
BANP_BED="$IN_DIR/BANP_mm0to2_noCGmm_SNVs_in_motifs_in_peaks.bed"

n_nrf1=$(wc -l < "$NRF1_BED" | tr -d " ")
n_banp=$(wc -l < "$BANP_BED" | tr -d " ")

echo -e "NRF1_mm0to2_noCGmm\tBANP_mm0to2_noCGmm" > "$OUT"
echo -e "${n_nrf1}\t${n_banp}" >> "$OUT"

echo "[DONE] Wrote $OUT"

# barplot snv counts in NRF1 and BANP motifs across all cancer types
Rscript -e '
library(ggplot2)
# Read counts
data <- read.table("./results/snv/SNV_counts_in_motifsnpeaks_NRF1_BANP_mm0to2_noCGmm.tsv", header=TRUE, stringsAsFactors=FALSE)
# Transform to long format
df_long <- data.frame(
  Motif = c("NRF1_mm0to2_noCGmm", "BANP_mm0to2_noCGmm"),
  SNV_Count = c(data$NRF1_mm0to2_noCGmm, data$BANP_mm0to2_noCGmm)
)
# Plot
p <- ggplot(df_long, aes(x=Motif, y=SNV_Count, fill=Motif)) +
  geom_bar(stat="identity", width=0.6) +
  theme_minimal() +
  labs(title="Number of SNVs in NRF1 and BANP motifs and peaks (mm0to2_noCGmm)", y="Number of SNVs", x="Motifnpeaks") +
  theme(plot.title = element_text(hjust=0.5, size=12), axis.text=element_text(size=10), axis.title=element_text(size=11), legend.position="none")
print(p)
ggsave("./results/snv/SNV_counts_in_motifsnpeaks_NRF1_BANP_mm0to2_noCGmm_barplot.pdf", plot=p, width=8, height=6, dpi=300, bg="white")
' 
cat ./results/snv/SNV_counts_in_motifsnpeaks_NRF1_BANP_mm0to2_noCGmm.tsv
# NRF1_mm0to2_noCGmm      BANP_mm0to2_noCGmm
# 34843   8787



#!!! DONE SNV NOT FOUND !!
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

cat ./snv/overlaps/intersected_motifs_snv_filtered/BANP_intersected_SNVs.bed |awk '$1=="chr9" && $2 <= 33290690 && $3 >= 33290690 {print $0}'
# pos_dbSNP = 33290690
cat ./snv/all_unique_variants_across_cancers.bed | grep 'chr9' | awk ' $2 >=33290688 && $3 <=33290695 {print $0}'
# chr9    33290694        33290695        chr9:33290695:T:C

# LOOKED INTO THE PAPER WHERE IT WAS FOUND , THEY MAPPED ON THE HG19 GENOME WHILE MINE ARE ON HG38. PLUS THE FOUND SNP IS A GERMLINE VARIATION WHILE MUTECT2 IS FOR SOMATIC MUTATIONS.



####################################################################################################
# Build a presence mutation only SNV matrix across the 34+4d samples across cancer type
####################################################################################################
# 1) filter the vcf files to keep only snv --> ./snv/filtered_SNV/
mkdir -p ./snv/filtered_SNV

max_jobs=10

for file in ./snv/snv_filtered_without_structural_variants/*.vcf.gz; do
    echo "Processing $file"
    bcftools view -v snps "$file" -Oz -o "./snv/filtered_SNV/$(basename "$file")" &

    while [ "$(jobs -r | wc -l)" -ge "$max_jobs" ]; do
        sleep 2
    done
done

wait

# 2) build a unique list of SNVs in motifs per cancer type 
awk -F'\t' '{
    # Split the last column by colon
    split($NF, a, ":"); 
    ref = a[3]; 
    alt = a[4]; 
    # Keep only if both Ref and Alt are 1 character long
    if (length(ref) == 1 && length(alt) == 1) print $0 
}' ./snv/overlaps/snv_in_motifs2mm/BANP_mm0to2_noCGmm_SNVs_in_motifs_SNVidxmotif.bed > ./snv/overlaps/snv_in_motifs2mm/BANP_mm0to2_noCGmm_SNVs_in_motifs_SNVidxmotif_filtered_snvs.bed

awk -F'\t' '{
    # Split the last column by colon
    split($NF, a, ":"); 
    ref = a[3]; 
    alt = a[4]; 
    # Keep only if both Ref and Alt are 1 character long
    if (length(ref) == 1 && length(alt) == 1) print $0 
}' ./snv/overlaps/snv_in_motifs2mm/NRF1_mm0to2_noCGmm_SNVs_in_motifs_SNVidxmotif.bed > ./snv/overlaps/snv_in_motifs2mm/NRF1_mm0to2_noCGmm_SNVs_in_motifs_SNVidxmotif_filtered_snvs.bed

# filter out OV samples from the metadata file to keep only 3d and 4d samples for the presence/absence matrix
awk '$2!="OV"' ./results/multi_omics/samples_3d+4d.tsv > ./results/multi_omics/samples_3d+4d_noOV.tsv

# 3)build a presence/absence matrix of SNVs in motifs across samples and cancer types (rows = SNVs, columns = samples, values = 1 if SNV is present in sample, 0 otherwise)
Rscript -e '
library(data.table)
library(parallel)

# ===== choose motif =====
# motif_name <- "NRF1_mm0to2_noCGmm"
motif_name <- "BANP_mm0to2_noCGmm"

overlap_file <- paste0("./snv/overlaps/snv_in_motifs2mm/", motif_name, "_SNVs_in_motifs_SNVidxmotif_filtered_snvs.bed")
meta_file    <- "./results/multi_omics/samples_3d+4d_noOV.tsv"
vcf_dir      <- "./snv/filtered_SNV"
outdir       <- file.path("./snv/SNVs_presence_matrix_motifs2mm", motif_name)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ===== helper: extract TCGA patient/sample ID like TCGA-OR-A5J1 =====
get_tcga_id <- function(x) {
  m <- regexpr("TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}", x)
  if (m[1] == -1) return(NA_character_)
  regmatches(x, m)
}

# ===== unique SNVs overlapping motifs =====
ov <- fread(overlap_file, sep = "\t", header = FALSE)
snvs <- sort(unique(ov[[ncol(ov)]]))   # last column = chr:pos:ref:alt

# ===== sample / cancer metadata =====
meta <- fread(meta_file, header = FALSE)
setnames(meta, c("sample","cancer","SNV","METH","RNA","ATAC"))

# ===== keep only VCFs whose extracted TCGA ID is in metadata =====
vcfs <- list.files(vcf_dir, pattern = "vcf.gz", full.names = TRUE)
vcf_samples <- vapply(basename(vcfs), get_tcga_id, character(1))
vcfs <- vcfs[vcf_samples %in% meta$sample]

cat("Number of motif SNVs:", length(snvs), "\n")
cat("Number of matched VCFs:", length(vcfs), "\n")
flush.console()

# ===== for one sample VCF, keep motif-SNVs present in that sample =====
get_hits <- function(f) {
  sample <- get_tcga_id(basename(f))
  cat("Processing sample:", sample, "\n")
  flush.console()

  dt <- fread(cmd = paste("zcat", shQuote(f)), sep = "\t", header = FALSE, skip = "#CHROM")
  if (nrow(dt) == 0) return(NULL)

  dt[, SNV_ID := paste(V1, V2, V4, V5, sep = ":")]
  hit <- unique(dt[SNV_ID %in% snvs, .(SNV = SNV_ID)])
  if (nrow(hit) == 0) return(NULL)

  hit[, sample := sample]
  hit
}

# ===== scan sample VCFs in parallel =====
res <- rbindlist(mclapply(vcfs, get_hits, mc.cores = 10), fill = TRUE)

cat("Number of sample-SNV hits found:", nrow(res), "\n")
flush.console()

if (nrow(res) > 0) {
  res <- merge(res, meta[, .(sample, cancer)], by = "sample")
}

# ===== build one matrix per cancer type =====
for (cc in sort(unique(meta$cancer))) {
  cat("Writing matrix for cancer type:", cc, "\n")
  flush.console()

  samples_cc <- meta[cancer == cc, sample]
  sub <- if (nrow(res) > 0) res[cancer == cc] else data.table()

  if (nrow(sub) > 0) {
    sub[, value := 1L]
    mat <- dcast(sub, SNV ~ sample, value.var = "value", fill = 0)
  } else {
    mat <- data.table(SNV = snvs)
  }

  missing_samples <- setdiff(samples_cc, names(mat))
  for (s in missing_samples) mat[, (s) := 0L]

  mat <- merge(data.table(SNV = snvs), mat, by = "SNV", all.x = TRUE)
  for (j in setdiff(names(mat), "SNV")) {
    set(mat, which(is.na(mat[[j]])), j, 0L)
  }

  fwrite(mat, file.path(outdir, paste0(cc, "_SNV_presence_matrix.tsv")), sep = "\t")
}
'

# 4) build one big matrix with all cancers together (rows = SNVs, columns = TCGA-cancer-samples, values = 1 if SNV is present in sample, 0 otherwise)
Rscript -e '
library(data.table)

indir <- "./snv/SNVs_presence_matrix_motifs2mm/BANP_mm0to2_noCGmm"
outfile <- "./snv/SNVs_presence_matrix_motifs2mm/BANP_mm0to2_noCGmm/all_cancers_SNV_presence_matrix.tsv"

files <- list.files(indir, pattern = "_SNV_presence_matrix.tsv$", full.names = TRUE)

lst <- lapply(files, function(f) {
  cancer <- sub("_SNV_presence_matrix.tsv", "", basename(f), fixed = TRUE)
  cat("Adding cancer type:", cancer, "\n")
  flush.console()

  dt <- fread(f)
  sample_cols <- setdiff(names(dt), "SNV")
  setnames(dt, sample_cols, paste(cancer, sample_cols, sep = "_"))
  dt
})

mat <- Reduce(function(x, y) merge(x, y, by = "SNV", all = TRUE), lst)

for (j in setdiff(names(mat), "SNV")) {
  set(mat, which(is.na(mat[[j]])), j, 0L)
}

fwrite(mat, outfile, sep = "\t")
'

# add motif information to the SNV presence matrix (chr:pos:ref:alt) by merging with the original overlap file that contains the motif information for each SNV
Rscript -e '
library(data.table)

motif_name   <- "NRF1_mm0to2_noCGmm"
overlap_file <- paste0("./snv/overlaps/snv_in_motifs2mm/", motif_name, "_SNVs_in_motifs_SNVidxmotif_filtered_snvs.bed")
matrix_file  <- paste0("./snv/SNVs_presence_matrix_motifs2mm/", motif_name, "/all_cancers_SNV_presence_matrix.tsv")
out_file     <- paste0("./snv/SNVs_presence_matrix_motifs2mm/", motif_name, "/all_cancers_SNV_presence_matrix_annotated.tsv")

cat("Step 1: Loading overlap file...\n"); flush.console()
ov <- fread(overlap_file, header = FALSE)

setnames(ov, c("motif_chr","motif_start","motif_end","motif_seq","motif_score","motif_strand","snv_chr","snv_start","snv_end","SNV"))

cat("Number of overlap rows:", nrow(ov), "\n"); flush.console()

cat("Step 2: Building SNV -> motif annotation...\n"); flush.console()
annot <- ov[, .(
  motif_name   = motif_name,
  motif_chr    = paste(unique(motif_chr), collapse=";"),
  motif_start  = paste(unique(motif_start), collapse=";"),
  motif_end    = paste(unique(motif_end), collapse=";"),
  motif_seq    = paste(unique(motif_seq), collapse=";"),
  motif_score  = paste(unique(motif_score), collapse=";"),
  motif_strand = paste(unique(motif_strand), collapse=";")
), by = SNV]

cat("Number of unique SNVs annotated:", nrow(annot), "\n"); flush.console()

cat("Step 3: Loading SNV matrix...\n"); flush.console()
mat <- fread(matrix_file)
cat("Matrix dimensions:", nrow(mat), "rows,", ncol(mat), "columns\n"); flush.console()

cat("Step 4: Merging annotation with matrix...\n"); flush.console()
mat_annot <- merge(annot, mat, by = "SNV", all.y = TRUE)

cat("Annotated matrix dimensions:", nrow(mat_annot), "rows,", ncol(mat_annot), "columns\n"); flush.console()

cat("Step 5: Writing output file...\n"); flush.console()
fwrite(mat_annot, out_file, sep = "\t")

cat("Done. Output written to:", out_file, "\n"); flush.console()
'

# add the linked gene information to the SNV presence matrix by merging with the original overlap file that contains the linked gene information for each SNV (if available)
Rscript -e '
library(data.table)

# motif_name   <- "NRF1_mm0to2_noCGmm"
motif_name <- "BANP_mm0to2_noCGmm"

overlap_file <- paste0("./snv/overlaps/snv_in_motifs2mm/", motif_name, "_SNVs_in_motifs_SNVidxmotif_filtered_snvs.bed")
gene_file    <- paste0("./motifs/", motif_name, "_noXY_closest_genes.bed")
matrix_file  <- paste0("./snv/SNVs_presence_matrix_motifs2mm/", motif_name, "/all_cancers_SNV_presence_matrix.tsv")
out_file     <- paste0("./snv/SNVs_presence_matrix_motifs2mm/", motif_name, "/all_cancers_SNV_presence_matrix_annotated_with_gene.tsv")

cat("Step 1: Loading overlap file...\n"); flush.console()
ov <- fread(overlap_file, header = FALSE)
setnames(ov, c("motif_chr","motif_start","motif_end","motif_seq","motif_score","motif_strand",
               "snv_chr","snv_start","snv_end","SNV"))
cat("Overlap rows:", nrow(ov), "\n"); flush.console()

cat("Step 2: Loading closest-gene file...\n"); flush.console()
g <- fread(gene_file, header = FALSE, skip = 1)
setnames(g, c("chrom_motif","start_motif","end_motif","motif","pvalue","strand_motif",
              "chrom_gene","start_gene","end_gene","gene","gene_id","transcript_id","strand_gene","distance"))
cat("Closest-gene rows:", nrow(g), "\n"); flush.console()

cat("Step 3: Merging motif overlaps with closest-gene annotations...\n"); flush.console()
ann0 <- merge(
  ov,
  g,
  by.x = c("motif_chr","motif_start","motif_end","motif_seq","motif_strand"),
  by.y = c("chrom_motif","start_motif","end_motif","motif","strand_motif"),
  all.x = TRUE
)
cat("Merged annotation rows:", nrow(ann0), "\n"); flush.console()

cat("Step 4: Collapsing to one annotation row per SNV...\n"); flush.console()
annot <- ann0[, .(
  motif_name    = motif_name,
  motif_chr     = paste(unique(na.omit(motif_chr)), collapse=";"),
  motif_start   = paste(unique(na.omit(motif_start)), collapse=";"),
  motif_end     = paste(unique(na.omit(motif_end)), collapse=";"),
  motif_seq     = paste(unique(na.omit(motif_seq)), collapse=";"),
  motif_score   = paste(unique(na.omit(motif_score)), collapse=";"),
  motif_strand  = paste(unique(na.omit(motif_strand)), collapse=";"),
  gene          = paste(unique(na.omit(gene)), collapse=";"),
  gene_id       = paste(unique(na.omit(gene_id)), collapse=";"),
  transcript_id = paste(unique(na.omit(transcript_id)), collapse=";"),
  gene_chr      = paste(unique(na.omit(chrom_gene)), collapse=";"),
  gene_start    = paste(unique(na.omit(start_gene)), collapse=";"),
  gene_end      = paste(unique(na.omit(end_gene)), collapse=";"),
  gene_strand   = paste(unique(na.omit(strand_gene)), collapse=";"),
  distance      = paste(unique(na.omit(distance)), collapse=";")
), by = SNV]
cat("Unique annotated SNVs:", nrow(annot), "\n"); flush.console()

cat("Step 5: Loading all-cancers matrix...\n"); flush.console()
mat <- fread(matrix_file)
cat("Matrix dimensions:", nrow(mat), "rows,", ncol(mat), "columns\n"); flush.console()

cat("Step 6: Merging annotation into matrix...\n"); flush.console()
mat_annot <- merge(annot, mat, by = "SNV", all.y = TRUE)
cat("Annotated matrix dimensions:", nrow(mat_annot), "rows,", ncol(mat_annot), "columns\n"); flush.console()

cat("Step 7: Writing output...\n"); flush.console()
fwrite(mat_annot, out_file, sep = "\t")

cat("Done. Output written to:", out_file, "\n"); flush.console()
'

#filter out SNV that are not present in any samples to reduce the size of the matrix.
Rscript -e '
library(data.table)

infile  <- "./snv/SNVs_presence_matrix_motifs2mm/BANP_mm0to2_noCGmm/all_cancers_SNV_presence_matrix_annotated_with_gene.tsv"
outfile <- "./snv/SNVs_presence_matrix_motifs2mm/BANP_mm0to2_noCGmm/all_cancers_SNV_presence_matrix_annotated_with_gene_filtered.tsv"

cat("Step 1: Loading matrix...\n"); flush.console()
dt <- fread(infile)

annot_cols <- c("SNV","motif_name","motif_chr","motif_start","motif_end","motif_seq","motif_score",
                "motif_strand","gene","gene_id","transcript_id","gene_chr","gene_start",
                "gene_end","gene_strand","distance")

sample_cols <- setdiff(names(dt), annot_cols)

cat("Rows before filtering:", nrow(dt), "\n"); flush.console()

dt[, total_presence := rowSums(.SD), .SDcols = sample_cols]
dt_filt <- dt[total_presence > 0]
dt_filt[, total_presence := NULL]

cat("Rows after filtering:", nrow(dt_filt), "\n"); flush.console()
cat("Rows removed:", nrow(dt) - nrow(dt_filt), "\n"); flush.console()

cat("Step 2: Writing filtered matrix...\n"); flush.console()
fwrite(dt_filt, outfile, sep = "\t")

cat("Done. Output written to:", outfile, "\n"); flush.console()
'

# Build the motif alteration matrix , where rows are motifs (instead of SNVs) and values indicate whether any SNV in that motif is present in each sample.
Rscript -e '
library(data.table)

infile  <- "./snv/SNVs_presence_matrix_motifs2mm/NRF1_mm0to2_noCGmm/all_cancers_SNV_presence_matrix_annotated_with_gene_filtered.tsv"
outfile <- "./snv/SNVs_presence_matrix_motifs2mm/NRF1_mm0to2_noCGmm/all_cancers_MOTIF_presence_matrix.tsv"

cat("Step 1: Loading SNV matrix...\n"); flush.console()
dt <- fread(infile)

annot_cols <- c(
  "SNV","motif_name","motif_chr","motif_start","motif_end","motif_seq","motif_score",
  "motif_strand","gene","gene_id","transcript_id","gene_chr","gene_start",
  "gene_end","gene_strand","distance"
)

sample_cols <- setdiff(names(dt), annot_cols)

cat("SNV rows before collapse:", nrow(dt), "\n"); flush.console()

cat("Step 2: Building motif IDs...\n"); flush.console()
dt[, motif_id := paste(motif_name, motif_chr, motif_start, motif_end, motif_strand, sep="|")]

cat("Unique motifs before collapse:", uniqueN(dt$motif_id), "\n"); flush.console()

cat("Step 3: Collapsing SNVs into motif-level presence...\n"); flush.console()
motif_dt <- dt[, c(
  .(
    motif_name    = first(motif_name),
    motif_chr     = first(motif_chr),
    motif_start   = first(motif_start),
    motif_end     = first(motif_end),
    motif_seq     = first(motif_seq),
    motif_score   = first(motif_score),
    motif_strand  = first(motif_strand),
    gene          = first(gene),
    gene_id       = first(gene_id),
    transcript_id = first(transcript_id),
    gene_chr      = first(gene_chr),
    gene_start    = first(gene_start),
    gene_end      = first(gene_end),
    gene_strand   = first(gene_strand),
    distance      = first(distance)
  ),
  lapply(.SD, max, na.rm = TRUE)
), by = motif_id, .SDcols = sample_cols]

cat("Rows after motif collapse:", nrow(motif_dt), "\n"); flush.console()

cat("Step 4: Writing motif-level matrix...\n"); flush.console()
fwrite(motif_dt, outfile, sep = "\t")

cat("Done. Output written to:", outfile, "\n"); flush.console()
'
# motif alteration frequency per cancer type by averaging the presence/absence values of all motifs in each cancer type.
Rscript -e '
library(data.table)

infile  <- "./snv/SNVs_presence_matrix_motifs2mm/NRF1_mm0to2_noCGmm/all_cancers_MOTIF_presence_matrix.tsv"
outfile <- "./snv/SNVs_presence_matrix_motifs2mm/NRF1_mm0to2_noCGmm/motif_frequency_per_cancer.tsv"

cat("Step 1: Loading motif-level matrix...\n"); flush.console()
dt <- fread(infile)

annot_cols <- c(
  "motif_id","motif_name","motif_chr","motif_start","motif_end","motif_seq","motif_score",
  "motif_strand","gene","gene_id","transcript_id","gene_chr","gene_start",
  "gene_end","gene_strand","distance"
)

sample_cols <- setdiff(names(dt), annot_cols)

cat("Number of motifs:", nrow(dt), "\n"); flush.console()
cat("Number of sample columns:", length(sample_cols), "\n"); flush.console()

m <- as.matrix(dt[, ..sample_cols])
rownames(m) <- dt$motif_id

cancer <- sub("_.*", "", colnames(m))

cat("Cancer types found:", paste(unique(cancer), collapse=", "), "\n"); flush.console()

cat("Step 2: Computing motif frequency per cancer...\n"); flush.console()
res <- rbindlist(lapply(unique(cancer), function(cc){
  cols <- which(cancer == cc)
  data.table(
    motif_id = rownames(m),
    cancer   = cc,
    freq     = rowMeans(m[, cols, drop=FALSE], na.rm=TRUE)
  )
}))

annot <- unique(dt[, .(
  motif_id, motif_name, motif_chr, motif_start, motif_end, motif_strand,
  motif_seq, motif_score, gene
)])

res <- merge(res, annot, by="motif_id", all.x=TRUE)

cat("Step 3: Writing frequency table...\n"); flush.console()
fwrite(res, outfile, sep="\t")

cat("Done. Output written to:", outfile, "\n"); flush.console()
'

# gene alteration 
Rscript -e '
library(data.table)

infile  <- "./snv/SNVs_presence_matrix_motifs2mm/BANP_mm0to2_noCGmm/all_cancers_MOTIF_presence_matrix.tsv"
outfile <- "./snv/SNVs_presence_matrix_motifs2mm/BANP_mm0to2_noCGmm/all_cancers_GENE_presence_matrix.tsv"

cat("Step 1: Loading motif matrix...\n"); flush.console()
dt <- fread(infile)

annot_cols <- c(
  "motif_id","motif_name","motif_chr","motif_start","motif_end","motif_seq",
  "motif_score","motif_strand","gene","gene_id","transcript_id","gene_chr",
  "gene_start","gene_end","gene_strand","distance"
)

sample_cols <- setdiff(names(dt), annot_cols)

cat("Motifs:", nrow(dt), "\n"); flush.console()
cat("True sample columns:", length(sample_cols), "\n"); flush.console()

# force sample columns to numeric
dt[, (sample_cols) := lapply(.SD, as.numeric), .SDcols = sample_cols]

cat("Step 2: Collapsing motifs into gene-level...\n"); flush.console()
gene_dt <- dt[, lapply(.SD, max, na.rm = TRUE), by = gene, .SDcols = sample_cols]

cat("Genes:", nrow(gene_dt), "\n"); flush.console()

cat("Step 3: Writing output...\n"); flush.console()
fwrite(gene_dt, outfile, sep = "\t")

cat("Done. Output written to:", outfile, "\n"); flush.console()
'

# frequency of gene alteration per cancer type by averaging the presence/absence values of all genes in each cancer type.
Rscript -e '
library(data.table)

infile  <- "./snv/SNVs_presence_matrix_motifs2mm/BANP_mm0to2_noCGmm/all_cancers_GENE_presence_matrix.tsv"
outfile <- "./snv/SNVs_presence_matrix_motifs2mm/BANP_mm0to2_noCGmm/gene_frequency_per_cancer.tsv"

cat("Step 1: Loading gene matrix...\n"); flush.console()
dt <- fread(infile)

sample_cols <- setdiff(names(dt), "gene")

cat("Genes:", nrow(dt), "\n"); flush.console()
cat("Samples:", length(sample_cols), "\n"); flush.console()

dt[, (sample_cols) := lapply(.SD, as.numeric), .SDcols = sample_cols]

m <- as.matrix(dt[, ..sample_cols])
rownames(m) <- dt$gene

cancer <- sub("_.*", "", colnames(m))

cat("Cancer types found:", paste(unique(cancer), collapse=", "), "\n"); flush.console()

cat("Step 2: Computing frequency per cancer...\n"); flush.console()
res <- rbindlist(lapply(unique(cancer), function(cc) {
  cols <- which(cancer == cc)
  data.table(
    gene   = rownames(m),
    cancer = cc,
    freq   = rowMeans(m[, cols, drop = FALSE], na.rm = TRUE)
  )
}))

cat("Step 3: Writing output...\n"); flush.console()
fwrite(res, outfile, sep = "\t")

cat("Done. Output written to:", outfile, "\n"); flush.console()
'

# plot the frequency of the top 100 genes with the highest alteration frequency across cancer types (bubble plot where x-axis = cancer type, y-axis = gene, size/color of points = alteration frequency)
Rscript -e '
library(data.table)
library(ggplot2)

infile  <- "./snv/SNVs_presence_matrix_motifs2mm/NRF1_mm0to2_noCGmm/gene_frequency_per_cancer.tsv"
outfile <- "./results/snv/gene_NRF1_frequency_bubbleplot_top100.pdf"

dt <- fread(infile)

top_genes <- dt[, .(max_freq = max(freq, na.rm=TRUE)), by=gene][order(-max_freq)][1:100, gene]
plot_dt <- dt[gene %in% top_genes & freq > 0]

gene_order <- plot_dt[, .(max_freq = max(freq, na.rm=TRUE)), by=gene][order(max_freq), gene]
plot_dt[, gene := factor(gene, levels = gene_order)]

pdf(outfile, width=15, height=18)

ggplot(plot_dt, aes(x=cancer, y=gene, size=freq)) +
  geom_point() +
  theme_bw() +
  labs(
    title = "Gene alteration frequency per cancer",
    x = "Cancer type",
    y = "Gene",
    size = "Frequency"
  ) +
  theme(
    axis.text.x = element_text(angle=45, hjust=1),
    panel.grid.major = element_line(linewidth=0.2),
    panel.grid.minor = element_blank()
  )

dev.off()
'
# Number of genes with snv in at least one sample
# NRF1
cat ./snv/SNVs_presence_matrix_motifs2mm/NRF1_mm0to2_noCGmm/gene_frequency_per_cancer.tsv | awk '{print $1}' | sort | uniq -cd | wc -l
# 9244
cat ./snv/SNVs_presence_matrix_motifs2mm/BANP_mm0to2_noCGmm/gene_frequency_per_cancer.tsv | awk '{print $1}' | sort | uniq -cd | wc -l
# 5039

################################################
# Multiple Mutations needed analysis 
################################################
# look into if there are mutations that coexist 








###############################################################################
# RNAseq script 
###############################################################################

###############################################################################
# 1) Link RNAseq files from papers to a expression directory. --> 11,505 files
###############################################################################
# link the expression files from the /data/papers/tcga/ to this folder :
mkdir -p ./expression/
for file in /data/papers/tcga/TCGA*/*/RNA_TCGA*.augmented_star_gene_counts.tsv; do
    ln -s "$file" ./expression/
done

###############################################################################
# 2) Count Number of RNAseq data per cancer type
###############################################################################
echo -e "Cancer_type\tNumber_of_RNAseq_files" > ./expression/number_of_rnaseq_per_cancer.txt
ls ./expression/ | grep augmented_star_gene_counts.tsv | cut -d'_' -f2 | cut -d'-' -f2 | sort | uniq -c | awk '{print $2 "\t" $1}' >> ./expression/number_of_rnaseq_per_cancer.txt

# Plot the number of RNAseq data per cancer type (barplot)
Rscript -e '
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
})

in_file <- "./expression/number_of_rnaseq_per_cancer.txt"
out_pdf <- "./results/expression/number_of_rnaseq_per_cancer.pdf"
dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)

# Read tab-delimited file
df <- read.delim(in_file, header = TRUE, stringsAsFactors = FALSE) %>%
  rename(
    Count  = Number_of_RNAseq_files,
    Cancer = Cancer_type
  )

# ---- total number of RNA-seq files
total_files <- sum(df$Count, na.rm = TRUE)

# Read defined color file
color_df <- read.table("./results/multi_omics/cancer_color_order_with_defined_colours.tsv",
                       header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                       comment.char = "", quote = "")

colnames(color_df) <- c("Cancer_full", "Color")
color_df$Cancer <- sub("^TCGA-", "", trimws(color_df$Cancer_full))
color_df$Color  <- trimws(gsub("\r", "", color_df$Color))

palette_main <- color_df$Color
names(palette_main) <- color_df$Cancer

# Keep all cancers (unknown ones will be grey)
palette_plot <- c(
  palette_main,
  setNames(
    rep("grey70", length(setdiff(unique(df$Cancer), names(palette_main)))),
    setdiff(unique(df$Cancer), names(palette_main))
  )
)

# Order bars by count (biggest → smallest)
df <- df %>% arrange(Count)
df$Cancer <- factor(df$Cancer, levels = df$Cancer)

p <- ggplot(df, aes(x = Cancer, y = Count, fill = Cancer)) +
  geom_col(color = "black") +
  coord_flip() +
  theme_bw() +
  scale_fill_manual(values = palette_plot) +
  theme(legend.position = "none") +
  labs(
    x = NULL,
    y = "Number of RNA-seq files",
    title = paste0("RNA-seq files per cancer type (total = ", total_files, ")")
  )

ggsave(out_pdf, plot = p, width = 10, height = 8)
'

###############################################################################
# 3) Count Number of RNAseq data per cancer type per healthy and cancer samples
###############################################################################

mkdir -p ./expression/expression_counts
out=./expression/expression_counts/expression_presence.tsv

echo -e "cancer\ttumor_count(01-09A)\thealthy_count(10-19A)" > "$out"

# get cancer list 
cancers=$(ls -1 ./expression/ \
  | grep augmented_star_gene_counts.tsv \
  | cut -d'_' -f2 \
  | cut -d'-' -f2 \
  | sort -u)

for cancer in $cancers; do
  # tumor: 01A..09A
  tumor_count=$(ls -1 ./expression/ \
    | grep augmented_star_gene_counts.tsv \
    | grep -E "_TCGA-${cancer}-.*-0[1-9]A" \
    | wc -l)

  # healthy: 10A..19A
  healthy_count=$(ls -1 ./expression/ \
    | grep augmented_star_gene_counts.tsv \
    | grep -E "_TCGA-${cancer}-.*-1[0-9]A" \
    | wc -l)

  printf "%s\t%d\t%d\n" "$cancer" "$tumor_count" "$healthy_count" >> "$out"
done

echo "[DONE] Wrote $out"

# Plot the number of RNAseq data per cancer type for tumor vs healthy samples (grouped barplot)
Rscript -e '
library(ggplot2)
library(reshape2)

# Load RNA-seq counts
df <- read.table("./expression/expression_counts/expression_presence.tsv",header = TRUE, sep = "\t")
colnames(df) <- c("cancer", "tumor", "healthy")

# Flag cancers with no healthy samples
df$no_healthy <- df$healthy == 0

# Total samples per cancer
df$total <- df$tumor + df$healthy

# Reorder cancers by total (descending)
df$cancer <- factor(df$cancer, levels = df$cancer[order(-df$total)])

# Label dataframe (after reordering)
label_df <- df[, c("cancer", "no_healthy")]

# Long format (tumor + healthy only)
df_long <- melt(df,
                id.vars = "cancer",
                measure.vars = c("tumor", "healthy"),
                variable.name = "type",
                value.name = "count")

pdf("./results/expression/expression_counts_healthy_cancer_per_cancer_type_barplot.pdf",width = 12, height = 8)

ggplot(df_long, aes(x = cancer, y = count, fill = type)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(
    values = c("tumor" = "salmon", "healthy" = "steelblue"),
    name = ""
  ) +
  labs(
    title   = "Tumor (0XA) vs Healthy (1XA) RNA-seq Sample Counts per Cancer Type",
    x       = "Cancer Type",
    y       = "Number of Samples",
    caption = "Red cancer names = cancer types with no healthy samples"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title   = element_text(hjust = 0.5),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    plot.caption = element_text(hjust = 0.5, color = "red", size = 12),
    plot.margin  = margin(t = 10, r = 20, b = 60, l = 20)
  ) +
  coord_cartesian(clip = "off") +
  geom_text(
    data = label_df,
    aes(
      x     = cancer,
      y     = 0,
      label = cancer,
      color = no_healthy
    ),
    angle       = 45,
    hjust       = 1,
    vjust       = 1.2,
    inherit.aes = FALSE,
    size        = 4
  ) +
  scale_color_manual(
    values = c(`TRUE` = "red", `FALSE` = "black"),
    guide  = "none"
  )

dev.off()
'

###############################################################################
# 4) Count Number of patients with RNA seq data 
###############################################################################
ls -1 ./expression/*augmented_star_gene_counts.tsv 2>/dev/null \
| sed 's#.*/##' \
| grep -oE 'TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}' \
| sort -u \
| wc -l

# 10091

###############################################################################
# 5) Count Number of GENE TYPE IN EACH FILE RNAseq data 
###############################################################################
awk -F'\t' 'NR>1 && $1 ~ /^ENS/ {print $3}' ./expression/RNA_TCGA-ACC-TCGA-OR-A5J1-01A_1.augmented_star_gene_counts.tsv | sort | uniq -c 
# plot the distribution of gene types in the RNA-seq files (barplot)
Rscript -e '
library(data.table)

infile <- "./expression/RNA_TCGA-ACC-TCGA-OR-A5J1-01A_1.augmented_star_gene_counts.tsv"
outpdf <- "./results/expression/barplot_gene_types_RNA_TCGA-ACC-TCGA-OR-A5J1-01A_1.pdf"
dir.create(dirname(outpdf), recursive=TRUE, showWarnings=FALSE)

dt <- fread(infile)
# keep gene rows only (those that start with ENS)
dt <- dt[grepl("^ENS", dt[[1]])]

# gene_type is column 3
cnt <- dt[, .N, by=.(gene_type = dt[[3]])]
setorder(cnt, -N)

# -----  colors + y-axis scaling -----
n <- nrow(cnt)
cols <- grDevices::hcl.colors(n, palette="Dark2", alpha=0.9)

ymax <- max(cnt$N, na.rm=TRUE)
ylim <- c(0, ymax * 1.12)  # headroom for text labels
# ---------------------------------------

pdf(outpdf, width=10, height=8)
par(mar=c(10,4,4,2))
bp <- barplot(
  cnt$N,
  names.arg=cnt$gene_type,
  las=2,
  ylab="Number of genes",
  main="Gene biotypes in file",
  col=cols,
  border=NA,
  ylim=ylim,
  cex.names = 0.6
)

text(bp, cnt$N, labels=cnt$N, pos=3, cex=0.5)
grid(nx=NA, ny=NULL)
dev.off()

cat("[DONE] -> ", outpdf, "\n", sep="")
'

###############################################################################
# 6) Range of tpm_unstranded for protein_coding
###############################################################################

echo -e "tpm_unstranded_values" > ./expression/tpm_unstranded_values.txt

for f in ./expression/*augmented_star_gene_counts.tsv; do
  awk '$1 ~ /^ENS/ && $3 ~ /^protein_coding$/ {print $7}' "$f" >> ./expression/tpm_unstranded_values.txt
done
cat ./expression/tpm_unstranded_values.txt | grep -v 'tpm' | sort -un > ./expression/tpm_unstranded_values_unique.txt 
head -1 ./expression/tpm_unstranded_values_unique.txt
# 0.0004
tail -1 ./expression/tpm_unstranded_values_unique.txt
# 337184.2864

# LOOK IF THERE ARE ANY PROTEIN CODING GENES WITH NA TPM_UNSTRANDED VALUES : NOTE: THERE ARE NO PROTEIN CODING GENES WITH NA TPM_UNSTRANDED VALUES.
for f in ./expression/*augmented_star_gene_counts.tsv; do
  echo "processing $f"
  awk '$3=="protein_coding" && $7=="NA" {print $0}' "$f" >> ./expression/protein_coding_genes_with_NA_tpm_unstranded.txt
done
cat ./expression/protein_coding_genes_with_NA_tpm_unstranded.txt | wc -l
# 0
rm ./expression/protein_coding_genes_with_NA_tpm_unstranded.txt


###############################################################################
# 7) Filter out genes on chrX,chrY
###############################################################################
# 1) Create a list of genes on chrX,chrY 
cat /data/genome/annotations/hg38_genes.bed | awk '$1=="chrY" || $1=="chrX" {print $4}' > ./expression/genes_chrX_chrY_to_filter_out.txt

# 2) Filter out genes that are on chrX,chrY from the expression files that I have:
for files in ./expression/*.augmented_star_gene_counts.tsv; do
  echo "Processing $files"
  grep '^ENS' "$files" \
    | awk 'NR==FNR {genes[$1]=1; next} !($2 in genes)' \
        ./expression/genes_chrX_chrY_to_filter_out.txt - \
    > "${files%.tsv}_filtered.tsv"
done

# 3) Check how many genes are left after filtering out chrX,chrY genes 
cat ./expression/RNA_TCGA-ACC-TCGA-OR-A5J2-01A_1.augmented_star_gene_counts_filtered.tsv | grep '^ENS' | grep protein_coding | wc -l
# 19117



##########################################################################################################################################################################################
# Create a summary table of presence/absence of omics data types (methylation, SNV, RNA-seq) for each sample, to know which samples have which data types available for later analyses.
##########################################################################################################################################################################################

#!/usr/bin/env bash
set -euo pipefail
# list of all samples from tcga --> samples_id cancer_type
ls -d /data/papers/tcga/TCGA-*/*/ 2>/dev/null | awk -F'/' '{c=$(NF-2); sub(/^TCGA-/,"",c); s=$(NF-1); print s"\t"c}' > ./all_tcga_samples_with_cancer.tsv

# DATA SUMMARY FOR EACH SAMPLE
# For each sample, check if there is a SNV file, a methylation file, and an RNA-seq file, and compile into a summary TSV.
set -euo pipefail
OUT="./results/multi_omics/sample_data_summary.tsv"
mkdir -p "$(dirname "$OUT")"
echo -e "sample\tcancer\tsnv\tmethylation\trnaseq\tatac" > "$OUT"

while IFS=$'\t' read -r sample_id cancer; do
  [[ -n "$sample_id" ]] || continue
  snv_file=$(ls -1 ./snv/snv_filtered_without_structural_variants/SNV_*"${sample_id}"*.vcf.gz 2>/dev/null || true)
  meth_file=$(ls -1 ./methylation/filtered_methylation/HM450_*"${sample_id}"*_annotated_methylation_filtered.bed.gz 2>/dev/null || true)
  rnaseq_file=$(ls -1 ./expression/RNA_*"${sample_id}"*.augmented_star_gene_counts.tsv 2>/dev/null || true)
  atac_file=$(ls -1 ./peaks/ATAC_TCGA*"${sample_id}"* 2>/dev/null || true)
  echo -e "${sample_id}\t${cancer}\t$( [[ -n "$snv_file" ]] && echo 1 || echo 0 )\t$( [[ -n "$meth_file" ]] && echo 1 || echo 0 )\t$( [[ -n "$rnaseq_file" ]] && echo 1 || echo 0 )\t$( [[ -n "$atac_file" ]] && echo 1 || echo 0 )" >> "$OUT"
done < ./all_tcga_samples_with_cancer.tsv

# heatmap of data availability across samples (samples on x-axis, data types on y-axis, colored by presence/absence)
Rscript -e '
library(data.table)
library(pheatmap)
library(ggplot2)
library(grid)
library(gridExtra)
in_tsv   <- "./results/multi_omics/sample_data_summary.tsv"
col_file <- "./results/multi_omics/cancer_color_order_with_defined_colours.tsv"
out_pdf  <- "./results/multi_omics/sample_data_summary_heatmap.pdf"
dt <- fread(in_tsv)
mat <- t(as.matrix(dt[, .(snv, methylation, rnaseq, atac)]))
colnames(mat) <- dt$sample
rownames(mat) <- c("Variant data","Methylation data","RNAseq","ATACseq")
ord <- order(dt$cancer, dt$sample)
mat <- mat[, ord, drop=FALSE]
ann <- data.frame(
  Cancer = paste0("TCGA-", as.character(dt$cancer[ord])),
  stringsAsFactors=FALSE
)
rownames(ann) <- colnames(mat)
cols <- fread(col_file, header=FALSE)[, 1:2]
setnames(cols, c("Cancer","Color"))
cols[, Cancer := gsub("\r","", trimws(as.character(Cancer)))]
cols[, Color  := gsub("\r","", trimws(as.character(Color)))]
cancer_cols_all <- setNames(cols$Color, cols$Cancer)
present <- sort(unique(ann$Cancer))
cancer_cols <- cancer_cols_all[present]
ann_colors <- list(Cancer = cancer_cols)

# ===== summaries =====
counts <- colSums(dt[, .(snv, methylation, rnaseq, atac)], na.rm=TRUE)
counts_dt <- data.table(dataset=names(counts), n=as.integer(counts))
counts_dt[, dataset := factor(dataset, levels=c("snv","methylation","rnaseq","atac"),labels=c("SNV","Methylation","RNA-seq","ATAC-seq"))]

p_counts <- ggplot(counts_dt, aes(x=dataset, y=n, fill=dataset)) +
  geom_col() +
  theme_bw() +
  labs(title="Number of samples with each dataset", x=NULL, y="Samples") +
  theme(axis.text.x = element_text(angle=30, hjust=1),
        legend.title = element_blank()) +
  scale_fill_brewer(palette="Set2")

# combos (SNV Meth RNA ATAC)
dt[, combo := paste0(snv, methylation, rnaseq, atac)]
combo_counts <- dt[, .N, by=combo][order(-N)]
setnames(combo_counts, "N", "n_samples")

# ===== many distinct, vivid colors =====
combo_levels <- combo_counts$combo
combo_cols <- setNames(
  hcl(h=seq(15, 375, length.out=length(combo_levels)+1)[1:length(combo_levels)],c=100, l=55),
  combo_levels
)

p_combo <- ggplot(combo_counts,aes(x=reorder(combo, -n_samples), y=n_samples, fill=combo)) +
  geom_col() +
  theme_bw() +
  labs(title="Sample counts per dataset combination (SNV Meth RNA ATAC)",
       x="Combination code (SNV Meth RNA ATAC)", y="Samples") +
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.title = element_blank()) +
  scale_fill_manual(values=combo_cols)

combo_counts[, combo_label := paste0("SNV=", substr(combo,1,1),
                                    " Meth=", substr(combo,2,2),
                                    " RNA=", substr(combo,3,3),
                                    " ATAC=", substr(combo,4,4))]
combo_tbl <- combo_counts[, .(combo, combo_label, n_samples)]
tg <- tableGrob(combo_tbl, rows=NULL, theme=ttheme_minimal(base_size=10))

# ===== one PDF, multiple pages =====
pdf(out_pdf, width=24, height=10)

# PAGE 1: heatmap
ph <- pheatmap(mat,
         cluster_rows=FALSE,
         cluster_cols=FALSE,
         annotation_col=ann,
         annotation_colors=ann_colors,
         show_colnames=FALSE,
         fontsize_row=12,
         border_color="grey90",
         color=c("skyblue","yellow"),
         breaks=c(-0.5, 0.5, 1.5),
         silent=TRUE)
grid.draw(ph$gtable)

# PAGE 2: dataset counts
grid.newpage()
print(p_counts)

# PAGE 3: combo barplot + table
grid.newpage()
grid.arrange(p_combo, tg, ncol=2, widths=c(1.1, 1.4))
dev.off()
cat("[DONE] ", out_pdf, "\n", sep="")
'

########################################################
# Expression of NRF1 and BANP in the different samples 
########################################################
# grep -w --> to find matches of the while word (gene name) and not partial matches
for files in ./expression/RNA_TCGA-*-*.augmented_star_gene_counts.tsv; do
  sample=$(basename "$files" .augmented_star_gene_counts.tsv)
  out="./expression/expression_of_NRF1_BANP/${sample}_NRF1_BANP_expression.tsv"
  mkdir -p "$(dirname "$out")"
  echo -e "gene_id\tgene_name\tgene_type\tunstranded\tstranded_first\tstranded_second\ttpm_unstranded\tfpkm_unstranded\tfpkm_uq_unstranded" > "$out"
  grep -w 'NRF1' "$files"  >> "$out"
  grep -w 'BANP' "$files"  >> "$out"
done


#############################################################################################################################################################################################################################################################
# NRF1 and BANP are ubiquitously expressed across TCGA cancer types (definition: median TPM > 1, all sample types included) --> create a badge for each gene with the percentage of cancer types in which it is expressed, and a subtitle with the number of cancer types in which it is expressed out of the total number of cancer types. 
#############################################################################################################################################################################################################################################################
# This is for the presentation to show that these factors are ubiquitously expressed across all samples.
Rscript -e '
library(data.table)
library(ggplot2)

expr <- fread("./expression/gene_expression_matrix_protein_coding.tsv.gz")
setnames(expr, 1, "sample")

genes <- c("NRF1","BANP")
THR_TPM <- 1

dt <- expr[, c("sample", genes), with=FALSE]
long <- melt(dt, id.vars="sample", variable.name="gene", value.name="TPM")
long[, TPM := as.numeric(TPM)]

# cancer type (ACC, BLCA, ...)
long[, cancer := tstrsplit(sample, "-", fixed=TRUE)[[2]]]

# expressed per cancer if median TPM > threshold
med <- long[, .(median_TPM = median(TPM, na.rm=TRUE)), by=.(cancer, gene)]
med[, expressed := median_TPM > THR_TPM]

# summarize per gene
sumdt <- med[, .(
  n_cancers = .N,
  n_expr    = sum(expressed, na.rm=TRUE),
  pct       = round(100 * mean(expressed, na.rm=TRUE))
), by=gene]

sumdt[, subtitle := sprintf("Expressed in %d / %d cancer types", n_expr, n_cancers)]
sumdt[, big := paste0(pct, "%")]
sumdt[, gene := factor(gene, levels=genes)]

# small helper to place two tiles nicely
sumdt[, x := as.integer(gene)]
sumdt[, y := 1]

p <- ggplot(sumdt, aes(x=x, y=y)) +
  # colorful cards
  geom_tile(aes(fill=gene), width=0.9, height=0.55, color="white", linewidth=1.2) +
  # big TF name
  geom_text(aes(label=as.character(gene)), y=1.12, fontface="bold", size=10, color="white") +
  # huge percent
  geom_text(aes(label=big), y=1.00, fontface="bold", size=14, color="white") +
  # subtitle
  geom_text(aes(label=subtitle), y=0.86, size=4.2, color="white") +
  # global headline
  annotate("text", x=1.5, y=1.55,
           label="Ubiquitously expressed across TCGA cancer types",
           fontface="bold", size=6) +
  annotate("text", x=1.5, y=1.40,
           label=paste0("Definition: median TPM > ", THR_TPM, " (all sample types included)"),
           size=4) +
  scale_fill_manual(values=c("NRF1"="#7B2FF7", "BANP"="#FFB000")) +
  coord_cartesian(xlim=c(0.5, 2.5), ylim=c(0.6, 1.7), expand=FALSE) +
  theme_void() +
  theme(legend.position="none")

dir.create("./results/expression", recursive=TRUE, showWarnings=FALSE)
ggsave("./results/expression/NRF1_BANP_ubiquity_badges_slide.pdf",
       p, width=9, height=3.2, device=cairo_pdf)

cat("WROTE -> ./results/expression/NRF1_BANP_ubiquity_badges_slide.pdf\n")
'
##################################################################################################################################################
# plot  the expression of NRF1 and BANP in the different samples, grouped by cancer type one page for healthy and one for tumor samples (boxplot)
##################################################################################################################################################
Rscript -e '
library(readr); library(dplyr); library(stringr); library(ggplot2)

cat("[INFO] Reading mapping file...\n")
map <- read_tsv("./all_tcga_samples_with_cancer.tsv", col_names=c("sample_id","cancer"), show_col_types=FALSE)
cat("[INFO] Listing NRF1/BANP per-sample files...\n")
files <- list.files("./expression/expression_of_NRF1_BANP",pattern = "_NRF1_BANP_expression.tsv$",full.names = TRUE)
cat("[INFO] Found ", length(files), " files\n", sep="")
cat("[INFO] Reading expression files...\n")
df <- bind_rows(lapply(files, function(f){
  cat("[READ] ", basename(f), "\n", sep="")
  read_tsv(f, show_col_types=FALSE) %>%
    mutate(sample_id = str_extract(basename(f), "TCGA-[A-Za-z0-9]{2}-[A-Za-z0-9]{4}-[0-9]{2}[A-Za-z]"))
}))

cat("[INFO] Preparing table + joining cancer...\n")
map <- map %>% mutate(case_id = str_replace(sample_id, "-[0-9]{2}[A-Za-z]$", ""))

df <- df %>%
  select(sample_id, gene_name, tpm_unstranded) %>%
  filter(gene_name %in% c("NRF1","BANP")) %>%
  mutate(
    tpm_unstranded = as.numeric(tpm_unstranded),   # raw TPM
    case_id = str_replace(sample_id, "-[0-9]{2}[A-Za-z]$", ""),
    code = str_match(sample_id, "-([0-9]{2})[A-Za-z]$")[,2],
    sample_type = case_when(
      str_detect(code, "^0") ~ "Tumor",
      str_detect(code, "^1") ~ "Normal",
      TRUE ~ "Other"
    )
  ) %>%
  left_join(map %>% select(case_id, cancer), by="case_id") %>%
  filter(!is.na(cancer), !is.na(tpm_unstranded)) %>%
  filter(sample_type %in% c("Normal","Tumor"))

# =====  gene-specific cutoffs + removed counts + filter =====
df <- df %>%
  mutate(cutoff = case_when(
    gene_name == "NRF1" ~ 50,
    gene_name == "BANP" ~ 30,
    TRUE ~ Inf
  ))
removed <- df %>%
  group_by(gene_name) %>%
  summarise(n_removed = sum(tpm_unstranded > cutoff, na.rm=TRUE), .groups="drop")
get_removed <- function(g) {
  x <- removed$n_removed[removed$gene_name == g]
  ifelse(length(x) == 0, 0, x)
}

df <- df %>% filter(tpm_unstranded <= cutoff)
# ======================================================================

cat("[INFO] Samples Tumor (0X):  ", length(unique(df$sample_id[df$sample_type=="Tumor"])), "\n", sep="")
cat("[INFO] Samples Normal (1X): ", length(unique(df$sample_id[df$sample_type=="Normal"])), "\n", sep="")
cat("[INFO] Cancers: ", length(unique(df$cancer)), "\n", sep="")

# ===== keep same y-scale across both plots  =====
ymax <- max(df$tpm_unstranded, na.rm=TRUE)

make_plot <- function(d, title_txt, removed_n) {
  ggplot(d, aes(cancer, tpm_unstranded)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color = sample_type), width = 0.2, alpha = 0.6, size = 0.7) +
    geom_hline(yintercept = 1, color = "blue", linewidth = 0.8) +
    scale_color_manual(values = c(Normal="black", Tumor="red")) +
    coord_cartesian(ylim = c(0, ymax)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle=60, hjust=1)) +
    labs(
      title = title_txt,
      x="Cancer type", y="TPM",
      color = paste0("Sample type (removed above cutoff: ", removed_n, ")")
    )
}

cat("[INFO] Building plots...\n")
p_nrf1 <- make_plot(filter(df, gene_name=="NRF1"),
                    "NRF1 TPM — Normal (black) vs Tumor (red) | cutoff=50 removed",
                    get_removed("NRF1"))
p_banp <- make_plot(filter(df, gene_name=="BANP"),
                    "BANP TPM — Normal (black) vs Tumor (red) | cutoff=30 removed",
                    get_removed("BANP"))

out <- "./results/expression/boxplot_TPM_NRF1_BANP_by_cancer_Normal_vs_Tumor.pdf"
dir.create(dirname(out), recursive=TRUE, showWarnings=FALSE)

cat("[INFO] Writing 2-page PDF...\n")
pdf(out, width=14, height=6)
if (nrow(filter(df, gene_name=="NRF1")) > 0) print(p_nrf1) else {plot.new(); title("No NRF1 samples found")}
if (nrow(filter(df, gene_name=="BANP")) > 0) print(p_banp) else {plot.new(); title("No BANP samples found")}
dev.off()

cat("[DONE] Wrote ", out, "\n", sep="")
'

####################################################################################################################################################################################################
# Build a list of the expression of NRF1 and BANP in all samples and then filter to keep only the samples for which we have both tumor and healthy samples
####################################################################################################################################################################################################
# 1) Build one tsv with all samples expression of NRF1 and BANP with columns : samples cancer gene tpm_unstranded
echo -e 'samples\tcancer\tsample_type\tgene\ttpm_unstranded' > ./expression/all_samples_NRF1_BANP_expression.tsv
OUT="./expression/all_samples_NRF1_BANP_expression.tsv"
for f in ./expression/expression_of_NRF1_BANP/RNA_TCGA-*.tsv; do
  echo "[READ] $f"
  sample_id=$(basename "$f" _NRF1_BANP_expression.tsv)     
  sample_id=${sample_id#RNA_TCGA-}          
  cancer=$(echo "$sample_id" | cut -d'-' -f1) 
  sample_type=$(echo "$sample_id" | cut -d'-' -f5)
  awk -v sid="$sample_id" -v c="$cancer" -v type="$sample_type" 'BEGIN{OFS="\t"} NR>1 {print sid, c,type, $2, $7}' "$f" >> "$OUT"
done
# 2) List the samples for which we have both tumor and healthy samples, and their expression of NRF1 and BANP (for later plotting)
OUT="./expression/patients_with_both_tumor_and_healthy_samples_with_tpm.tsv"
echo -e "case_id\tcancer\tgene\ttumor_sample\ttumor_tpm\tnormal_sample\tnormal_tpm" > "$OUT"
tail -n +2 ./expression/all_samples_NRF1_BANP_expression.tsv \
| awk -F'\t' 'BEGIN{OFS="\t"}
{
  cancer = $2
  st = substr($3,1,2); if(st!="01" && st!="11") next;
  split($1,a,"-");
  case_id = a[1] "-" a[2] "-" a[3] "-" a[4];   
  g = $4
  tpm = $5
  if(!(case_id in c_of)) c_of[case_id] = cancer
  if(st=="01"){ tum[case_id]=$1; tt[case_id,g]=tpm }
  if(st=="11"){ nor[case_id]=$1; nt[case_id,g]=tpm }
}
END{
  for (k in tt) {
    split(k,b,SUBSEP); cid=b[1]; g=b[2];
    if ((cid in nor) && ((cid SUBSEP g) in nt)) {
      print cid, c_of[cid], g, tum[cid], tt[cid,g], nor[cid], nt[cid,g]
    }
  }
}' | sort >> "$OUT"

# Samples with NRF1 and BANP expression < 1 :
cat ./expression/all_samples_NRF1_BANP_expression.tsv | awk '$5<1 {print $0}' > ./expression/samples_with_verylow_expression.tsv
cat ./expression/samples_with_verylow_expression.tsv | wc -l 
# 46

# number of each cancer type in those samples:
cat ./expression/samples_with_verylow_expression.tsv | awk '{print $2}' | uniq -c

# 9 ACC <-- there a 2 for the same sample so there is 8 samples to remove here (1 for NRF1 and 1 for BANP)
# 1 BRCA
# 1 CESC
# 1 COAD
# 1 HNSC
# 3 KICH
# 8 KIRC
# 4 KIRP
# 5 LIHC
# 1 PAAD
# 3 PCPG
# 2 STAD
# 4 THCA
# 1 UCEC
# 2 UVM

# number of samples in total:
cat ./expression/all_samples_NRF1_BANP_expression.tsv| awk '{print $2}' | uniq -c

#158 ACC
#    862 BLCA
#   2462 BRCA
#    618 CESC
#    88 CHOL
#   1048 COAD
#    96 DLBC
#    396 ESCA
#    782 GBM
#   1132 HNSC
#    182 KICH
#   1228 KIRC
#    646 KIRP
#    302 LAML
#   1068 LGG
#    848 LIHC
#   1202 LUAD
#   1124 LUSC
#    174 MESO
#    868 OV
#    366 PAAD
#    374 PCPG
#   1108 PRAD
#    354 READ
#    530 SARC
#    946 SKCM
#   896 STAD
#   312 TGCT
#  1144 THCA
#   244 THYM
#  1178 UCEC
#   114 UCS
#   160 UVM

##################################################################################################################################################
# PCA of the expression of genes across all samples, to see if they cluster by cancer type or by healthy vs tumor status
# it removes genes that are completely flat (variance = 0) across all samples, PCA = 50, seed 1 , perp uses a dynamic calculation min 30.
##################################################################################################################################################
# create a file with the color order for the cancers (to use the same colors as before in the PCA plot)
Rscript -e '
library(data.table)

in_map  <- "./all_tcga_samples_with_cancer.tsv"
out_col <- "./results/multi_omics/cancer_color_order_with_defined_colours.tsv"
dir.create(dirname(out_col), recursive=TRUE, showWarnings=FALSE)

m <- fread(in_map, header=FALSE)
setnames(m, c("sample_id","cancer"))

cancers <- sort(unique(paste0("TCGA-", m$cancer)))

# nice 4-color set repeated (simple + readable)
base_cols <- c("#4E79A7", "#59A14F", "#F28E2B", "#B07AA1")
cols <- rep(base_cols, length.out=length(cancers))

dt <- data.table(cancer=cancers, color=cols)
fwrite(dt, out_col, sep="\t", col.names=FALSE)

cat("[DONE] Wrote ", out_col, " with ", nrow(dt), " cancers\n", sep="")
'
# create matrix of TPM values for protein_coding genes across all samples (genes on rows, samples on columns)
Rscript -e '
library(data.table)
TPM_COL <- 7
files <- list.files("./expression/",pattern="augmented_star_gene_counts_filtered.tsv$",full.names=TRUE)
cat("[STEP] Reading first file to define protein_coding gene set...\n")
dt0 <- fread(files[1], header=FALSE, showProgress=FALSE)
# V2 = gene_name, V3 = gene_type
pc_idx <- which(dt0$V3 == "protein_coding")
genes  <- dt0$V2[pc_idx]
G <- length(genes)
S <- length(files)
cat("[INFO] Protein_coding genes:", G, " | Samples:", S, "\n")
cat("[STEP] Allocating genes x samples matrix...\n")
m <- matrix(NA_real_, nrow=G, ncol=S)
rownames(m) <- genes
cat("[STEP] Filling matrix with TPM column ", TPM_COL, "...\n", sep="")
for (i in seq_along(files)) {
  if (i %% 50 == 0) cat("[READ]", i, "/", S, "\n")
  dt <- fread(files[i], header=FALSE, showProgress=FALSE)
  # safety: ensure row counts match reference
  if (nrow(dt) != nrow(dt0)) {
    cat("[SKIP] Row mismatch:", basename(files[i]), "\n")
    next
  }
  m[, i] <- suppressWarnings(as.numeric(dt[[TPM_COL]][pc_idx]))
}
cat("[STEP] Converting to samples x genes matrix...\n")
X <- t(m)  # samples x genes
colnames(X) <- genes
sample_ids <- sub("^RNA_", "",sub("[.]augmented_star_gene_counts_filtered[.]tsv$","", basename(files)))
rownames(X) <- sample_ids
cat("[INFO] Final matrix dims (samples x genes):", nrow(X), "x", ncol(X), "\n")
cat("[STEP] Saving to gene_expression_matrix_protein_coding.tsv.gz ...\n")
gene_expression_matrix <- X
fwrite(as.data.table(gene_expression_matrix, keep.rownames="sample"),file="./expression/gene_expression_matrix_protein_coding.tsv.gz",sep="\t")
cat("[DONE] Saved ./expression/gene_expression_matrix_protein_coding.tsv.gz\n")
'

# plot tSNE of the samples using the matrix of TPM values for protein_coding genes, colored by cancer type (using the same colors as before) as input.
Rscript -e '
library(data.table)
library(ggplot2)

cat("[INFO] Loading required packages...\n")
if (!requireNamespace("irlba", quietly=TRUE)) stop("Install irlba: install.packages(\"irlba\")")
if (!requireNamespace("Rtsne", quietly=TRUE)) stop("Install Rtsne: install.packages(\"Rtsne\")")
if (!requireNamespace("plotly", quietly=TRUE)) stop("Install plotly: install.packages(\"plotly\")")
if (!requireNamespace("htmlwidgets", quietly=TRUE)) stop("Install htmlwidgets: install.packages(\"htmlwidgets\")")

# ===== INPUT MATRIX (samples x genes) =====
mat_file <- "./expression/gene_expression_matrix_protein_coding.tsv.gz"
col_file <- "./results/multi_omics/cancer_color_order_with_defined_colours.tsv"
out_html <- "./results/expression/TSNE_from_matrix_TPM_interactive.html"
cat("[INFO] Loading cancer colors from:", col_file, "\n")
col_dt <- fread(col_file, header=FALSE)
if (ncol(col_dt) < 2) stop("color file must have >=2 columns")
col_dt <- col_dt[, 1:2]
setnames(col_dt, c("cancer","color"))
cancer_levels <- as.character(col_dt$cancer)
cancer_cols   <- setNames(as.character(col_dt$color), as.character(col_dt$cancer))
cat("[INFO] Reading matrix:", mat_file, "\n")
dtm <- fread(mat_file)          # expects first col named "sample" (as you saved it)
samples <- as.character(dtm$sample)
dtm[, sample := NULL]
xs <- as.matrix(dtm)
mode(xs) <- "numeric"
rownames(xs) <- samples
cat("[INFO] Matrix dims:", nrow(xs), "samples x", ncol(xs), "genes\n")

cat("[INFO] Removing zero-variance genes...\n")
v <- apply(xs, 2, var)
xs <- xs[, v > 0, drop=FALSE]
cat("[INFO] Removed", sum(v == 0), "zero-variance genes; remaining:", ncol(xs), "\n")

cat("[INFO] Running PCA (50 dims) with irlba...\n")
set.seed(1)
pca50 <- irlba::prcomp_irlba(xs, n=50, center=TRUE, scale.=TRUE)

cat("[INFO] Running t-SNE...\n")
n <- nrow(xs)
perp <- min(30, floor((n - 1) / 3))
if (perp < 5) perp <- 5
cat("[INFO] t-SNE perplexity =", perp, " (n_samples =", n, ")\n")

ts <- Rtsne::Rtsne(pca50$x, dims=2, perplexity=perp, pca=FALSE,check_duplicates=FALSE, verbose=TRUE)

cat("[INFO] Building plot...\n")
df <- data.frame(sample=rownames(xs), TSNE1=ts$Y[,1], TSNE2=ts$Y[,2], stringsAsFactors=FALSE)

tok <- strsplit(df$sample, "-", fixed=TRUE)
df$cancer <- vapply(tok, function(v) paste0("TCGA-", v[2]), character(1))

last_tok <- vapply(tok, function(v) v[length(v)], character(1))
code2 <- substr(last_tok, 1, 2)
df$type <- ifelse(substr(code2,1,1)=="0","Tumor",
           ifelse(substr(code2,1,1)=="1","Normal","Other"))
df$type <- factor(df$type, levels=c("Tumor","Normal","Other"))
df$cancer <- factor(df$cancer, levels=cancer_levels)

df$hover <- paste0("Sample: ", df$sample,
                   "<br>Cancer: ", as.character(df$cancer),
                   "<br>Type: ", as.character(df$type))

p <- ggplot(df, aes(TSNE1, TSNE2, color=cancer, shape=type)) +
  geom_point(aes(tooltip = hover), size=2, alpha=0.85) +
  scale_shape_manual(values=c(Tumor=16, Normal=17, Other=4), drop=FALSE) +
  scale_color_manual(values=cancer_cols, breaks=cancer_levels, drop=FALSE, na.value="grey70") +
  theme_bw() +
  labs(
    title="Interactive t-SNE — from expression matrix",
    subtitle=paste0("Genes used=", ncol(xs), " | Samples=", nrow(df),
                    " | PCA dims=50 | perplexity=", perp),
    color="Cancer", shape="Type"
  )

cat("[INFO] Saving HTML to:", out_html, "\n")
g <- plotly::ggplotly(p, tooltip = "tooltip")
dir.create(dirname(out_html), recursive=TRUE, showWarnings=FALSE)
htmlwidgets::saveWidget(g, out_html, selfcontained=TRUE)

cat("[DONE] Wrote ", out_html, "\n", sep="")
'
# add this above if i want the log2(TPM+1) transform for the tSNE  before this cat("[INFO] Removing zero-variance genes...\n")
# cat("[INFO] log2(TPM+1) transform...\n")
# xs <- log2(xs + 1)


####################################################################################################
# TF binding in motifs : link the TF motif instances to the genes assigned to them and the probes that overlap them, to get a table of (TF, probe, gene, motif_instance_id) for all probes in motifs and their assigned genes.
#####################################################################################################
# Create table cpgxgene : TF probe gene motif_instance_id(chr:start:end:strand) --> ./methylation/probe_gene_pairs_in_motifs.tsv
# This code links each CpG probe that overlaps a TF motif instance (from the overlap files) to the gene assigned to that motif instance (from the closest_genes.bed files). It produces a table of (TF, probe, gene, motif_instance_id) for all probes in motifs and their assigned genes.
# columns in the output file in the header:
# TF      probe   gene    motif_instance_id(chr:start:end:strand)
Rscript -e '
library(data.table)
TF_LIST <- c("BANP_mm0to2_noCGmm","NRF1_mm0to2_noCGmm")
overlap_dir <- "./motifs/overlaps/intersected_motifs2mm_HM450"
gene_dir    <- "./motifs"
out_file    <- "./methylation/probe_gene_pairs_in_motifs.tsv"
dir.create(dirname(out_file), recursive=TRUE, showWarnings=FALSE)

# ---- motif_instance -> gene from closest_genes.bed
gene_lookup <- rbindlist(lapply(TF_LIST, function(tf){
  f <- file.path(gene_dir, paste0(tf, "_noXY_closest_genes.bed"))
  if (!file.exists(f)) stop("Missing: ", f)
  g <- fread(f, header=FALSE, select=c(1,2,3,6,10))
  setnames(g, c("V1","V2","V3","V6","V10"), c("chr","s","e","st","gene"))
  g[, motif_instance_id := paste(chr, as.integer(s), as.integer(e), st, sep=":")]
  unique(g[, .(TF=tf, motif_instance_id, gene)])
}), fill=TRUE)

pg_list <- list()
k <- 0L

for (tf in TF_LIST) {
  f <- file.path(overlap_dir, paste0(tf, "_probeXmotif.tsv"))
  if (!file.exists(f)) {
    cat("[SKIP] missing overlap:", f, "\n")
    next
  }

  # overlap file has no header; columns:
  # 1-3 probe chr/start/end, 4 probeID
  # 5-7 motif chr/start/end, 8 motif sequence, 9 pval, 10 strand
  x <- fread(f, header=FALSE)
  if (ncol(x) < 10) stop("Unexpected columns in: ", f, " (ncol=", ncol(x), ")")

  dt <- x[, .(
    TF = tf,
    probe = as.character(V4),
    motif_instance_id = paste(as.character(V5), as.integer(V6), as.integer(V7), as.character(V10), sep=":")
  )]
  dt <- unique(dt)
  dt <- merge(dt, gene_lookup[TF==tf], by=c("TF","motif_instance_id"), all.x=TRUE)

  # keep only rows with a gene assigned
  dt <- dt[!is.na(gene) & gene != ""]
  if (nrow(dt) == 0) {
    cat("[WARN] no probe-gene pairs after join for", tf, "\n")
    next
  }

  k <- k + 1L
  pg_list[[k]] <- unique(dt[, .(TF, probe, gene, motif_instance_id)])
  cat("[OK]", tf, "pairs:", nrow(pg_list[[k]]), "\n")
}

pg <- rbindlist(pg_list, fill=TRUE)
if (nrow(pg) == 0) stop("No probe-gene pairs produced.")

# Save
fwrite(pg, out_file, sep="\t")
cat("WROTE -> ", out_file, " (rows=", nrow(pg), ", unique probes=", uniqueN(pg$probe), ")\n", sep="")
'
# number of pairs gene probe:
# BANP_mm0to2_noCGmm pairs: 1053 
# NRF1_mm0to2_noCGmm pairs: 3875

# number of gene with multiple probes in them:
cat ./methylation/probe_gene_pairs_in_motifs.tsv | grep BANP_mm0to2_noCGmm | awk '{print $3}' | uniq -c | sort -k1,1nr | awk '$1>1 {print $0}' | wc -l
# 169 genes for BANP that have multiple probes in them 
cat ./methylation/probe_gene_pairs_in_motifs.tsv | grep NRF1_mm0to2_noCGmm | awk '{print $3}' | uniq -c | sort -k1,1nr | awk '$1>1 {print $0}' | wc -l
# 801 genes for NRF1 that have multiple probes in them 

# Add the distance information in the ./methylation/probe_gene_pairs_in_motifs.tsv file with columns : TF probe gene motif_instance_id dist_to_tss is_proximal proximal_or_distal
Rscript - <<'EOF'
suppressPackageStartupMessages(library(data.table))

# input files
pairs_file <- "./methylation/probe_gene_pairs_in_motifs.tsv"
dist_file  <- "./methylation/dist_to_tss/HM450_TCGA-ACC-TCGA-OR-A5J1-01A_1_annotated_methylation_filtered_cpg_dist_tss.tsv"

# output
out_file <- "./methylation/probe_gene_pairs_in_motifs_with_tss.tsv.gz"

# load motif-probe-gene table
pairs <- fread(pairs_file)

# load distance table (no header in your example)
dist <- fread(dist_file, header = FALSE)
setnames(dist, c("probe", "dist_to_tss"))

# clean types
pairs[, probe := as.character(probe)]
dist[, probe := as.character(probe)]
dist[, dist_to_tss := as.numeric(dist_to_tss)]

# if duplicates exist in the distance file, keep one per probe
dist <- unique(dist, by = "probe")

# merge
out <- merge(pairs, dist, by = "probe", all.x = TRUE)

# define proximal/distal
out[, is_proximal := !is.na(dist_to_tss) & abs(dist_to_tss) < 2000]
out[, proximal_or_distal := fifelse(
  is.na(dist_to_tss), NA_character_,
  fifelse(abs(dist_to_tss) < 2000, "proximal", "distal")
)]

# reorder columns
setcolorder(out, c("TF", "probe", "gene", "motif_instance_id", "dist_to_tss", "is_proximal", "proximal_or_distal"))

# write
fwrite(out, out_file, sep = "\t")

cat("Saved:", out_file, "\n")
cat("Rows:", nrow(out), "\n")
cat("Matched probes:", out[!is.na(dist_to_tss), .N], "\n")
cat("Unmatched probes:", out[is.na(dist_to_tss), .N], "\n")
cat("Proximal:", out[is_proximal == TRUE, .N], "\n")
cat("Distal:", out[proximal_or_distal == 'distal', .N], "\n")
EOF




# Create a table of the samples for which we have 3d+4d data (SNV, METH, RNA) and whether we have ATAC data for them or not, to see how many samples we have with 3d+4d data and how many of those have ATAC data.
# create a tsv of the samples with 3d+4d:
cat ./results/multi_omics/sample_data_summary.tsv | awk '$3==1 && $4==1 && $5==1 && ($6==0 ||$6==1) {print $0}' > ./results/multi_omics/samples_3d+4d.tsv

# barplot the number of samples witht 3d+4d per cancer type:
Rscript -e '
library(data.table)
library(ggplot2)

dt <- fread("./results/multi_omics/samples_3d+4d.tsv", header = FALSE)
setnames(dt, c("sample","cancer","SNV","METH","RNA","ATAC"))

counts <- dt[, .N, by = cancer]
total_samples <- nrow(dt)

# Read defined color file
color_df <- fread(
  "./results/multi_omics/cancer_color_order_with_defined_colours.tsv",
  header = FALSE
)
setnames(color_df, c("Cancer_full", "Color"))
color_df[, Cancer := sub("^TCGA-", "", trimws(Cancer_full))]
color_df[, Color  := trimws(gsub("\r", "", Color))]
palette <- setNames(color_df$Color, color_df$Cancer)

# Match colors to cancers in the plot
counts[, Color := palette[cancer]]
counts[is.na(Color), Color := "grey70"]

p <- ggplot(counts, aes(x = reorder(cancer, N), y = N, fill = cancer)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = N), hjust = -0.15, size = 3) +
  coord_flip() +
  scale_fill_manual(values = setNames(counts$Color, counts$cancer)) +
  labs(
    x = "Cancer type",
    y = "Number of samples",
    title = "Samples per Cancer Type",
    subtitle = paste0("Multi-omics cohort (n = ", total_samples, ")")
  ) +
  expand_limits(y = max(counts$N) * 1.08) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave("./results/multi_omics/samples_3+4datasets_barplot.pdf",plot = p, width = 5, height = 5, dpi = 300, bg = "white")
'

#################################################################################################
# create pearson correlation plot between delta_beta and delta_expr for the probe-gene pairs
#################################################################################################
# a) create a gene expressionmatrix of the 3d+4d samples with no OV cancer type:
f="./results/multi_omics/samples_3d+4d.tsv"

zcat ./expression/gene_expression_matrix_protein_coding.tsv.gz | \
awk '
  NR==FNR {
    # cohort file: keep patient IDs except OV
    if ($2 != "OV") keep[$1]=1
    next
  }
  FNR==1 { print; next }  # keep header

  {
    if (FNR % 500 == 0)
      printf("[awk] processed %d expression rows | kept %d\n", FNR, kept) > "/dev/stderr"

    # extract patient barcode from expression sample id (no {n} intervals)
    if (match($1, /TCGA-[A-Z0-9][A-Z0-9]-[A-Z0-9][A-Z0-9][A-Z0-9][A-Z0-9]/)) {
      pat = substr($1, RSTART, RLENGTH)
      if (pat in keep) { print; kept++ }
    }
  }
  END {
    printf("[awk] DONE | kept %d rows (plus header)\n", kept) > "/dev/stderr"
  }
' "$f" - \
> ./expression/gene_expression_matrix_3d4d_noOV.tsv
echo "Done! Saved: ./expression/gene_expression_matrix_3d4d_noOV.tsv"


# to look at overlapping motifs 
awk 'BEGIN{FS=OFS="\t"}
NR==1 {next}
{
  if ($1==prev_chr && $2 < prev_end) print prev_line "\n" $0 "\n"
  prev_chr=$1
  prev_end=$3
  prev_line=$0
}' ./motifs/BANP_mm0to2_noCGmm_noXY_closest_genes.bed | awk '$1!="chrY" && $1!="chrM" && $1!="chrX" {print $0}' 


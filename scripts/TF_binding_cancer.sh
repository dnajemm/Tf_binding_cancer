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
# 1) MOTIF FILTERING (6-MER STRINGENCY) , I DONT USE IT.
###############################################################################
# Filter motifs with p-value threshold and sort them.
mkdir -p ./motifs

# zcat /data/genome/motifs/jaspar_2024/fimo/vertebrata/hg38/NRF1.MA0506.3.bed.gz | awk '($5<=1/4^6)' > ./motifs/NRF1_filtered_6mer.MA0506.3.bed
# zcat ./motifs/BANP.bed.gz | awk '($5<=1/4^6)' > ./motifs/BANP_filtered_6mer.MA2503.1.bed

###############################################################################
# MOTIF COUNT BEFORE / AFTER FILTERING THE 6mer
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
ggsave("./results/motifs/motif_counts_before_after_6mer_filtering.pdf", width = 6, height = 4, dpi = 300, bg = "white")
'
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

###############################################################################
# MOTIF COUNT AFTER FILTERING 6mer
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
ggsave("./results/motifs/motif_counts_after_6mer_filtering.pdf", width = 5, height = 4, dpi = 300, bg = "white")
'

###############################################################################
# MOTIF LISTS (UNIQUE + NOT UNIQUE) for the 6mer
###############################################################################
mkdir -p ./motifs/motifs_list

for file in ./motifs/*_filtered_6mer.*.bed; do
  base=$(basename "$file" .bed)

  # not-unique (keeps all occurrences), sorted for easier counting
  awk '{print $4}' "$file" | sort > "./motifs/motifs_list/${base}_motifs_not_unique.txt"

  # unique (one per distinct sequence)
  awk '{print $4}' "$file" | sort -u > "./motifs/motifs_list/${base}_motifs_unique.txt"
done

###############################################################################
# MOTIF MISMATCH SUMMARY for the 6 mer (3 PLOTS / TF)
#  - Plot 1: positional mismatch tolerance (unique variants)
#  - Plot 2: number of unique variants by total mismatches (includes bin 0)
#  - Plot 3: number of motif occurrences (hits) by total mismatches (includes bin 0)
###############################################################################
Rscript -e '
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ===============================
# INPUTS
# ===============================
files <- list(
  NRF1 = "./motifs/motifs_list/NRF1_filtered_6mer.GCGCATGCGC_motifs_not_unique.txt",
  BANP = "./motifs/motifs_list/BANP_filtered_6mer.MA2503.1_motifs_not_unique.txt"
)

NRF1motif <- "GCGCATGCGC"   # length 10
BANPmotif <- "TCTCGCGAGA"     # length 10

out_pdf <- "./results/motifs/motif_variant_mismatch_visual_summary_TEST.pdf"
dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)

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
  # LOAD NOT-UNIQUE HITS + COUNT OCCURRENCES PER VARIANT
  # -----------------------------
  all_hits <- fread(files[[mot]], header = FALSE)$V1

  hit_counts <- as.data.table(table(all_hits))
  setnames(hit_counts, c("variant", "hits"))
  hit_counts[, hits := as.integer(hits)]

  # Unique variants used for mismatch matrix
  seqs <- hit_counts$variant

  ref <- if (mot == "NRF1") NRF1motif else BANPmotif
  L   <- nchar(ref)

  # ensure same length as reference
  if (!all(nchar(seqs) == L)) {
    bad <- unique(nchar(seqs)[nchar(seqs) != L])
    stop(paste0("[ERROR] ", mot, ": some sequences have length != ", L,
                " (found lengths: ", paste(bad, collapse = ","), "). Check your input column 4 content."))
  }

  hit_mat <- split_seq(seqs)
  ref_mat <- split_seq(rep(ref, length(seqs)))

  mismatch_mat <- (hit_mat != ref_mat)
  n_variants <- nrow(mismatch_mat)

  # -----------------------------
  # CORE POSITIONS
  # -----------------------------
  if (mot == "NRF1") {
    core_pos <- c(2, 3, 8, 9)   # CG + CG (for 12-mer NRF1 ref)
  } else {
    core_pos <- c(4, 5, 6, 7)     # CGCG (for 10-mer BANP ref)
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

  # attach occurrence counts
  dt <- merge(dt, hit_counts, by = "variant", all.x = TRUE)
  dt[is.na(hits), hits := 0L]

  # -----------------------------
  # PLOT 1 — POSITIONAL MISMATCH TOLERANCE (UNIQUE VARIANTS)
  # -----------------------------
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
        y = "Number of unique variants"
      ) +
      coord_cartesian(ylim = c(-1, n_variants)) +
      theme(
        panel.grid.major.x = element_blank(),
        axis.text.x = element_blank()
      )
  )

  # -----------------------------
  # CONFIGURATION PER VARIANT
  # -----------------------------
  dt[, config := fifelse(
    total == 0, "perfect_match",
    fifelse(core == total, "all_core",
            fifelse(core == 0, "no_core",
                    fifelse(core == 1, "mixed_1CG", "mixed_2CG")))
  )]

  cfg_pct <- dt[, .N, by = .(total, config)]
  cfg_pct[, pct := round(100 * N / sum(N)), by = total]

  labels_df <- cfg_pct[, .(
    label = paste0(
      "Perfect: ", ifelse(length(pct[config == "perfect_match"]) == 0, 0, pct[config == "perfect_match"]), "%\n",
      "Core: ",    ifelse(length(pct[config == "all_core"]) == 0, 0, pct[config == "all_core"]), "%\n",
      "Mixed 2CG: ", ifelse(length(pct[config == "mixed_2CG"]) == 0, 0, pct[config == "mixed_2CG"]), "%\n",
      "Mixed 1CG: ", ifelse(length(pct[config == "mixed_1CG"]) == 0, 0, pct[config == "mixed_1CG"]), "%\n",
      "NoCore: ",  ifelse(length(pct[config == "no_core"]) == 0, 0, pct[config == "no_core"]), "%"
    )
  ), by = total]

  # -----------------------------
  # PLOT 2 — UNIQUE VARIANTS BY TOTAL MISMATCHES
  # -----------------------------
  breaks_vec <- 0:max(dt$total)

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
      scale_x_continuous(breaks = breaks_vec, limits = c(-0.5, max(breaks_vec) + 0.5)) +
      theme_minimal() +
      labs(
        title = paste(mot, "- unique variants by mismatch count"),
        subtitle = "Counts unique motif sequences (variants). Labels show config % within each mismatch bin.",
        x = "Number of mismatches",
        y = "Number of unique variants"
      )
  )

  # -----------------------------
  # PLOT 3 — OCCURRENCES (HITS) BY TOTAL MISMATCHES
  # -----------------------------
  occ_df <- dt[, .(occurrences = sum(hits)), by = total][order(total)]

  print(
    ggplot(occ_df, aes(x = total, y = occurrences)) +
      geom_col(fill = "grey40") +
      geom_text(
        data = occ_df[total == 0],
        aes(label = occurrences),
        vjust = -0.4
      ) +
      scale_x_continuous(breaks = breaks_vec, limits = c(-0.5, max(breaks_vec) + 0.5)) +
      theme_minimal() +
      labs(
        title = paste(mot, "- motif occurrences by mismatch count"),
        subtitle = "Counts total motif hits (not unique sequences).",
        x = "Number of mismatches",
        y = "Number of motif occurrences"
      )
  )
}

dev.off()

cat("[DONE] PDF written:", out_pdf, "\\n")
'

###############################################################################
# 2) MOTIF DISTANCE TO TSS
###############################################################################
# Distance of 6mer motifs to TSS : generate a file per sample with the distances of each motif to the nearest TSS
mkdir -p ./motifs/dist_to_tss/

for f in ./motifs/*_filtered_6mer.*.bed; do
    base=$(basename "$f")                     
    sample=${base%_filtered_6mer.*.bed} 
    echo "Processing $sample ..."

    closestBed -t first -d -a "$f" -b /data/genome/annotations/hg38_tss.bed | awk '{print $NF}' > "./motifs/dist_to_tss/${sample}_dist_tss.txt"
done

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
# Overlap 6mer motifs with annotated methylation data to find intersecting regions and 
# see if the motifs fall within methylated probes regions.

mkdir -p ./motifs/overlaps/intersected_motifs_HM450

# intersect NRF1 6mer motifs with annotated methylation data
bedtools intersect -u -a ./motifs/NRF1_filtered_6mer.MA0506.3.bed -b ./methylation/annotated_methylation_data_probes_filtered.bed > ./motifs/overlaps/intersected_motifs_HM450/NRF1_intersected_methylation.bed
bedtools intersect -u -a ./motifs/BANP_filtered_6mer.MA2503.1.bed -b ./methylation/annotated_methylation_data_probes_filtered.bed > ./motifs/overlaps/intersected_motifs_HM450/BANP_intersected_methylation.bed

# Summarize the intersected results : 1   NRF1 29401 BANP 2919
mkdir -p ./results/motifs

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

mkdir -p ./motifs/overlaps/intersected_motifs2mm_HM450

# intersect NRF1 6mer motifs with annotated methylation data
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
out_pdf <- "./results/motifs/NRF1_BANP_motifs_vs_methylation_venn.pdf"
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
      paste0("./motifs/", motif, "_mm0to2.bed.gz")
    } else {
      paste0("./motifs/", motif, "_mm0to2_noCGmm.bed.gz")
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

# Overlap the 2mm motifs with the peaks to see which motifs fall within the peaks regions.
mkdir -p ./motifs/overlaps/motif2mm_peak_overlaps

# Loop through each 6mer motif file and intersect with all peak files from peaks

for file in ./motifs/*_mm0to2_noCGmm.bed.gz ./motifs/*_mm0to2.bed.gz; do
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

###############################################################################
# 5) MOTIF ∩ PEAK ∩ HM450
###############################################################################
# Overlap the intersected methylation 6mer motifs with the peak overlapping motifs to find common regions
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
# create a barplot : One page per TF: 6 mer motifs, ATAC, HM450, SNVs
Rscript -e '
library(ggplot2)

out_pdf <- "./results/motifs/NRF1_BANP_3D_vs_peaks_vs_HM450_barplots.pdf"
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

# create a barplot : One page per TF: 2mm motifs, ATAC, HM450, SNVs
Rscript -e '
library(ggplot2)
library(scales)

out_pdf <- "./results/motifs/NRF1_BANP_mm0to2_variants_vs_peaks_vs_HM450_SNV_barplots.pdf"
dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)

pdf(out_pdf, width = 10, height = 8)

plot_one <- function(motif, variant){

  motif_all_file <- if(variant == "mm0to2"){
    paste0("./motifs/", motif, "_mm0to2.bed.gz")
  } else {
    paste0("./motifs/", motif, "_mm0to2_noCGmm.bed.gz")
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
    snv   = "./motifs/overlaps/motifs_snv_overlap/BANP_filtered_SNVs_in_motifskept.bed"
  ),
  NRF1 = list(
    motif = "./motifs/NRF1_filtered_6mer.MA0506.3.bed",
    peak  = "./motifs/overlaps/motif_peak_overlaps/NRF1_filtered_6mer.MA0506.3_peak_overlaps.bed",
    hm450 = "./motifs/overlaps/intersected_motifs_HM450/NRF1_intersected_methylation.bed",
    snv   = "./motifs/overlaps/motifs_snv_overlap/NRF1_filtered_SNVs_in_motifskept.bed"
  )
)

out_pdf <- "./results/motifs/BANP_NRF1_motif_peaks_probes_snv_barplots.pdf"
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

###############################################################################
# PROPORTIONAL EULER DIAGRAM PER TF : 6 mer MOTIFS, ATAC, HM450, SNVs 
###############################################################################
# Proportional Euler diagram of motifs with peaks and methylation probes 
# Venn diagram 3D with filtered SNVs in motifs and in peak 
# Proportional Euler diagram of motifs with peaks and methylation probes

Rscript -e '
library(eulerr)

out_pdf <- "./results/motifs/Euler_VENN_NRF1_BANP_Motifs_ATAC_HM450_SNV.pdf"
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
  read.table("./motifs/overlaps/motifs_snv_overlap/NRF1_filtered_SNVs_in_motifskept.bed"),
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
  read.table("./motifs/overlaps/motifs_snv_overlap/BANP_filtered_SNVs_in_motifskept.bedd"),
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
    motif_all_file  <- paste0("./motifs/", motif, "_", variant, ".bed.gz")
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
for files in ./motifs/*_mm0to2*.bed.gz; do 
    echo "Processing $files ..."
    zcat "$files" | sort -k1,1 -k2,2n | closestBed -t "first" -d -a stdin -b /data/genome/annotations/hg38_tss.bed > "${files%.bed.gz}_closest_genes.bed"
done

# add header to the closest genes files
for file in ./motifs/*_mm0to2*_closest_genes.bed; do
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
    > ./results/multi_omics/cancer_color_order.txt

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
all_cancers <- scan("./results/multi_omics/cancer_color_order.txt", what = "")
palette <- rainbow(length(all_cancers))
names(palette) <- all_cancers
bar_colors <- palette[cancer_names]

pdf("./results/peaks/peaks_per_cancer_type_barplot.pdf", width=12, height=9 , bg = "white")

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
# Look into the ATAC peaks overlapping with 6mer motifs
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

# Look into the ATAC peaks overlapping with 2mm motifs
mkdir -p ./peaks/overlaps/intersected_motifs/NRF1_mm0to2

motif_sites="./motifs/NRF1_mm0to2.bed.gz"

for f in ./peaks/filtered_peaks/ATAC_TCGA*peaks_macs.bed; do
  sample=$(basename "$f" _peaks_macs.bed)
  out="./peaks/overlaps/intersected_motifs/NRF1_mm0to2/${sample}_NRF1_intersected_peaks.bed"

  bedtools intersect -u -a "$f" -b <(zcat "$motif_sites") > "$out"
  echo "Wrote: $out"
done

mkdir -p ./peaks/overlaps/intersected_motifs/BANP_mm0to2

motif_sites="./motifs/BANP_mm0to2.bed.gz"

for f in ./peaks/filtered_peaks/ATAC_TCGA*peaks_macs.bed; do
  sample=$(basename "$f" _peaks_macs.bed)
  out="./peaks/overlaps/intersected_motifs/BANP_mm0to2/${sample}_BANP_intersected_peaks.bed"

  bedtools intersect -u -a "$f" -b <(zcat "$motif_sites") > "$out"
  echo "Wrote: $out"
done

mkdir -p ./peaks/overlaps/intersected_motifs/NRF1_mm0to2_noCGmm

motif_sites="./motifs/NRF1_mm0to2_noCGmm.bed.gz"

for f in ./peaks/filtered_peaks/ATAC_TCGA*peaks_macs.bed; do
  sample=$(basename "$f" _peaks_macs.bed)
  out="./peaks/overlaps/intersected_motifs/NRF1_mm0to2_noCGmm/${sample}_NRF1_intersected_peaks.bed"

  bedtools intersect -u -a "$f" -b <(zcat "$motif_sites") > "$out"
  echo "Wrote: $out"
done

mkdir -p ./peaks/overlaps/intersected_motifs/BANP_mm0to2_noCGmm

motif_sites="./motifs/BANP_mm0to2_noCGmm.bed.gz"

for f in ./peaks/filtered_peaks/ATAC_TCGA*peaks_macs.bed; do
  sample=$(basename "$f" _peaks_macs.bed)
  out="./peaks/overlaps/intersected_motifs/BANP_mm0to2_noCGmm/${sample}_BANP_intersected_peaks.bed"

  bedtools intersect -u -a "$f" -b <(zcat "$motif_sites") > "$out"
  echo "Wrote: $out"
done


###############################################################################
# 7) 6 mer MOTIF THAT OVERLAP WITH PEAKS DISTANCE TO TSS
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

out_pdf <- "./results/peaks/ATAC_motifs_TSS/dist_tss_ATAC_NRF1_BANP_byCancer.pdf"
dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)

pdf(out_pdf, width = 11.7, height = 8.3)
par(bg = "white")

# =========================
# RUN ONCE FOR EACH TF
# =========================
for (TF in c("NRF1", "BANP")) {

  dist_dir <- paste0("./peaks/dist_to_tss_motifs/", TF)
  out_tsv  <- paste0("./results/peaks/ATAC_motifs_TSS/", TF,"/summary_", TF, "_byCancer.tsv")

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
distances <- read.table(dist_file)[,1]

pdf("./results/methylation/dist_tss_methylation_histogram.pdf", width=10, height=10)
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
# SMOOTH SCATTER PLOT FOR ALL CANCER-HEALTHY SAMPLE PAIRS FILTERED with 6mer motif
###########################################################################################
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
  NRF1 = "./methylation/overlaps/intersected_motifs_HM450/NRF1_intersected_methylation.bed",
  BANP = "./methylation/overlaps/intersected_motifs_HM450/BANP_intersected_methylation.bed"
)

out_pdf <- "./results/methylation/smoothScatter_vs_NRF1_vs_BANP_per_pair.pdf"
thr <- 20

dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)

# =========================
# LOAD DATA
# =========================
pairs <- fread(pairs_file)
stopifnot(all(c("cancer","patient","tumor_file","healthy_file") %in% colnames(pairs)))

motif_cpgs <- lapply(motif_files, function(f) unique(fread(f, header = FALSE)[[4]]))

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

read_bed_gz <- function(path) {
  cmd <- paste("zcat", shQuote(path))
  fread(cmd = cmd, header = FALSE, showProgress = FALSE)
}

# =========================
# START PDF
# =========================
pdf(out_pdf, width = 16, height = 6, colormodel = "rgb", useDingbats = FALSE)
par(mfrow = c(1, 3), mar = c(4, 4, 4, 1))

pages_written <- 0L

# =========================
# LOOP OVER PAIRS
# =========================
if (nrow(pairs) == 0) {
  plot.new()
  title(main = paste0("No rows in pairs file:\n", pairs_file))
  pages_written <- 1L
} else {

  for (i in seq_len(nrow(pairs))) {

    cat("Processing:", pairs$cancer[i], pairs$patient[i], "\n")

    tryCatch({

      # ---- read files (robust quoting) ----
      tf <- pairs$tumor_file[i]
      nf <- pairs$healthy_file[i]

      if (!file.exists(tf) || !file.exists(nf)) {
        plot.new()
        title(main = paste0(
          pairs$cancer[i], " | ", pairs$patient[i], "\n",
          "Missing file(s):\n",
          if (!file.exists(tf)) paste0("tumor_file: ", tf, "\n") else "",
          if (!file.exists(nf)) paste0("healthy_file: ", nf) else ""
        ))
        pages_written <- pages_written + 1L
        next
      }

      tumor  <- read_bed_gz(tf)
      normal <- read_bed_gz(nf)

      # ---- extract probe + beta ----
      tumor_df  <- data.frame(probe = tumor[[4]],  meth_tumor  = tumor[[5]])
      normal_df <- data.frame(probe = normal[[4]], meth_normal = normal[[5]])

      merged <- merge(tumor_df, normal_df, by = "probe")
      if (nrow(merged) == 0) {
        plot.new()
        title(main = paste0(
          pairs$cancer[i], " | ", pairs$patient[i], "\n",
          "Merged = 0 CpGs (no shared probes)"
        ))
        pages_written <- pages_written + 1L
        next
      }

      # (optional but recommended) drop NA beta
      merged <- merged[!is.na(merged$meth_tumor) & !is.na(merged$meth_normal), ]
      if (nrow(merged) == 0) {
        plot.new()
        title(main = paste0(
          pairs$cancer[i], " | ", pairs$patient[i], "\n",
          "All shared probes had NA beta"
        ))
        pages_written <- pages_written + 1L
        next
      }

      delta <- merged$meth_tumor - merged$meth_normal

      # =========================
      # PLOT 1 — ALL CpGs + %
      # =========================
      smoothScatter(
        merged$meth_normal,
        merged$meth_tumor,
        xlab = "Healthy methylation (%)",
        ylab = "Tumor methylation (%)",
        main = paste0(pairs$cancer[i], " | ", pairs$patient[i], "\nAll CpGs")
      )
      abline(0, 1, col = "black", lwd = 2)
      abline(a = -thr, b = 1, col = "blue", lwd = 2)
      abline(a =  thr, b = 1, col = "blue", lwd = 2)
      add_class_legend(delta, thr, title = "All CpGs")

      # =========================
      # PLOT 2 — NRF1
      # =========================
      smoothScatter(
        merged$meth_normal,
        merged$meth_tumor,
        xlab = "Healthy methylation (%)",
        ylab = "Tumor methylation (%)",
        main = "NRF1 motif CpGs in red (|DeltaBeta| ≥ 20%)"
      )
      in_NRF1 <- merged$probe %in% motif_cpgs$NRF1
      is_NRF1_red <- in_NRF1 & (abs(delta) >= thr)
      points(
        merged$meth_normal[is_NRF1_red],
        merged$meth_tumor[is_NRF1_red],
        pch = 16, cex = 0.5, col = rgb(1, 0, 0, 0.6)
      )
      abline(0, 1, col = "black", lwd = 2)
      abline(a = -thr, b = 1, col = "blue", lwd = 2)
      abline(a =  thr, b = 1, col = "blue", lwd = 2)
      add_class_legend(delta[in_NRF1], thr, title = "NRF1 CpGs")

      # =========================
      # PLOT 3 — BANP
      # =========================
      smoothScatter(
        merged$meth_normal,
        merged$meth_tumor,
        xlab = "Healthy methylation (%)",
        ylab = "Tumor methylation (%)",
        main = "BANP motif CpGs in red (|DeltaBeta| ≥ 20%)"
      )
      in_BANP <- merged$probe %in% motif_cpgs$BANP
      is_BANP_red <- in_BANP & (abs(delta) >= thr)
      points(
        merged$meth_normal[is_BANP_red],
        merged$meth_tumor[is_BANP_red],
        pch = 16, cex = 0.5, col = rgb(1, 0, 0, 0.6)
      )
      abline(0, 1, col = "black", lwd = 2)
      abline(a = -thr, b = 1, col = "blue", lwd = 2)
      abline(a =  thr, b = 1, col = "blue", lwd = 2)
      add_class_legend(delta[in_BANP], thr, title = "BANP CpGs")

      pages_written <- pages_written + 1L

    }, error = function(e) {
      plot.new()
      title(main = paste0(
        pairs$cancer[i], " | ", pairs$patient[i], "\nERROR:\n", conditionMessage(e)
      ))
      pages_written <- pages_written + 1L
    })
  }
}

# If literally nothing got plotted , add a diagnostic page
if (pages_written == 0L) {
  plot.new()
  title(main = paste0(
    "No pages produced.\nCheck pairs file and file paths:\n", pairs_file
  ))
}

dev.off()
cat("DONE →", out_pdf, " (pages:", pages_written, ")\n")
'

###########################################################################################
# Overlaps Mehtylaiton probes with 2mm motifs
########################################################################################### 
mkdir -p ./methylation/overlaps/intersected_motifs2mm_HM450

PROBES=./methylation/annotated_methylation_data_probes_filtered.bed

# NRF1 mm0to2
bedtools intersect -u \
  -a "$PROBES" \
  -b <(zcat ./motifs/NRF1_mm0to2.bed.gz) \
  > ./methylation/overlaps/intersected_motifs2mm_HM450/NRF1_mm0to2_intersected_methylation.bed

# NRF1 mm0to2 noCGmm
bedtools intersect -u \
  -a "$PROBES" \
  -b <(zcat ./motifs/NRF1_mm0to2_noCGmm.bed.gz) \
  > ./methylation/overlaps/intersected_motifs2mm_HM450/NRF1_mm0to2_noCGmm_intersected_methylation.bed

# BANP mm0to2
bedtools intersect -u \
  -a "$PROBES" \
  -b <(zcat ./motifs/BANP_mm0to2.bed.gz) \
  > ./methylation/overlaps/intersected_motifs2mm_HM450/BANP_mm0to2_intersected_methylation.bed

# BANP mm0to2 noCGmm
bedtools intersect -u \
  -a "$PROBES" \
  -b <(zcat ./motifs/BANP_mm0to2_noCGmm.bed.gz) \
  > ./methylation/overlaps/intersected_motifs2mm_HM450/BANP_mm0to2_noCGmm_intersected_methylation.bed


###########################################################################################
# Overlaps Mehtylaiton probes with 2mm motifs in peaks 
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

out_pdf <- "./results/methylation/smoothScatter_mm0to2_vs_noCGmm_mm0to2_per_pair.pdf"
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
out_pdf <- "./results/methylation/smoothScatter_HEALTHY_vs_HEALTHY_one_cancer_type.pdf"
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
# QUANTIFICATION OF % OF CpGs WITH |DELTA| >= 20% IN NRF1 AND BANP MOTIFS (mm0to2 vs noCGmm mm0to2) AND IN ATAC PEAKS
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
counts_file <- "./methylation/delta20_cg_lists/cg_presence_matrices/NRF1_mm0to2_noCGmm/CpG_counts_per_cancer_NRF1_mm0to2_noCGmm.tsv"
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

# Proximal CpGs (<2000bp): 1264
# Distal CpGs (≥2000bp): 529


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

out_pdf <- "./results/methylation/smoothScatter_replicates_mm0to2_vs_noCGmm_mm0to2.pdf"
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
out_pdf    <- "./results/methylation/stacked_median_composition_delta_methylation_per_cancer.pdf"
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

###############################################################################
# SECTION 9 — Distance OF CpG PROBES IN 6 mer MOTIFS TO TSS
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

# for each tsv file in ./methylation/dist_to_tss/ filter only those in motif of NRF1 and make a histogram of distances to TSS for each sample
# NRF1 
mkdir -p ./results/methylation/methylation_motifs_TSS/NRF1
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
mkdir -p ./results/methylation/methylation_motifs_TSS/BANP
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

# Plot histogram of distances to TSS for NRF1 and BANP motifs per cancer type
Rscript -e '
library(data.table)

out_pdf <- "./results/methylation/methylation_motifs_TSS/dist_tss_methylation_NRF1_BANP_byCancer.pdf"
dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)

pdf(out_pdf, width = 11.7, height = 8.3)
par(bg = "white")

# =========================
# RUN ONCE FOR EACH TF
# =========================
for (TF in c("NRF1", "BANP")) {

  dist_dir <- paste0("./methylation/dist_to_tss_motif/", TF)
  out_tsv  <- paste0("./results/methylation/methylation_motifs_TSS/", TF,"/summary_", TF, "_byCancer.tsv")

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
all_cancers <- scan("./results/multi_omics/cancer_color_order.txt", what = "")
palette <- rainbow(length(all_cancers))
names(palette) <- all_cancers
bar_colors <- palette[cancer_names]

pdf("./results/snv/snv_per_cancer_type_barplot.pdf",width=12, height=9, bg = "white")

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
###############################################################################
# 7) MOTIFS (NRF1 / BANP)
###############################################################################
# kept means kept the motifs after the intersection.
# !!!NOT Done AGAIN , DONE FOR THE MM0TO2M!!!
# Overlap filtered SNV in TF 6mer motifs and ATAC peaks 
mkdir -p ./snv/overlaps/snv_motifs_peaks_overlap/
bedtools intersect -u -a ./snv/all_unique_variants_across_cancers.bed -b ./motifs/overlaps/motif_peak_overlaps/NRF1_filtered_6mer.MA0506.3_peak_overlaps.bed > ./snv/overlaps/snv_motifs_peaks_overlap/NRF1_SNVs_in_motifs_in_peaks.bed
bedtools intersect -u -a ./snv/all_unique_variants_across_cancers.bed -b ./motifs/overlaps/motif_peak_overlaps/BANP_filtered_6mer.MA2503.1_peak_overlaps.bed > ./snv/overlaps/snv_motifs_peaks_overlap/BANP_SNVs_in_motifs_in_peaks.bed
# Overlap filtered SNV in TF motifs only
mkdir -p ./snv/overlaps/snv_motifs_overlap/
bedtools intersect -u -a ./snv/all_unique_variants_across_cancers.bed -b ./motifs/NRF1_filtered_6mer.MA0506.3.bed > ./snv/overlaps/snv_motifs_overlap/NRF1_intersected_SNVs.bed 
bedtools intersect -u -a ./snv/all_unique_variants_across_cancers.bed -b ./motifs/BANP_filtered_6mer.MA2503.1.bed > ./snv/overlaps/snv_motifs_overlap/BANP_intersected_SNVs.bed

mkdir -p ./motifs/overlaps/motifs_snv_overlap/
# overlap with NRF1 motifs 
bedtools intersect -u \
  -a ./motifs/NRF1_filtered_6mer.MA0506.3.bed \
  -b ./snv/all_unique_variants_across_cancers.bed \
  > ./motifs/overlaps/motifs_snv_overlap/NRF1_filtered_SNVs_in_motifskept.bed

# overlap with BANP motifs 
bedtools intersect -u \
  -a ./motifs/BANP_filtered_6mer.MA2503.1.bed \
  -b ./snv/all_unique_variants_across_cancers.bed \
  > ./motifs/overlaps/motifs_snv_overlap/BANP_filtered_SNVs_in_motifskept.bed

# overlap with NRF1 motifs in peaks 
bedtools intersect -u \
  -a ./motifs/overlaps/motif_peak_overlaps/NRF1_filtered_6mer.MA0506.3_peak_overlaps.bed \
  -b ./snv/all_unique_variants_across_cancers.bed \
  > ./snv/overlaps/snv_motifs_peaks_overlap/NRF1_filtered_SNVs_in_motifs_in_peaks_kept.bed

# overlap with BANP motifs in peaks
bedtools intersect -u \
  -a ./motifs/overlaps/motif_peak_overlaps/BANP_filtered_6mer.MA2503.1_peak_overlaps.bed \
  -b ./snv/all_unique_variants_across_cancers.bed \
  > ./snv/overlaps/snv_motifs_peaks_overlap/BANP_filtered_SNVs_in_motifs_in_peaks_kept.bed

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
  -b <(zcat ./motifs/NRF1_mm0to2.bed.gz) \
  > ./snv/overlaps/snv_in_motifs2mm/NRF1_mm0to2_SNVs_in_motifs.bed

bedtools intersect -u \
  -a ./snv/all_unique_variants_across_cancers.bed \
  -b <(zcat ./motifs/NRF1_mm0to2_noCGmm.bed.gz) \
  > ./snv/overlaps/snv_in_motifs2mm/NRF1_mm0to2_noCGmm_SNVs_in_motifs.bed

# BANP
bedtools intersect -u \
  -a ./snv/all_unique_variants_across_cancers.bed \
  -b <(zcat ./motifs/BANP_mm0to2.bed.gz) \
  > ./snv/overlaps/snv_in_motifs2mm/BANP_mm0to2_SNVs_in_motifs.bed

bedtools intersect -u \
  -a ./snv/all_unique_variants_across_cancers.bed \
  -b <(zcat ./motifs/BANP_mm0to2_noCGmm.bed.gz) \
  > ./snv/overlaps/snv_in_motifs2mm/BANP_mm0to2_noCGmm_SNVs_in_motifs.bed


# -------------------------
# 2) Motifs that contain SNVs (motif output)
# -------------------------
# NRF1 
bedtools intersect -u \
  -a <(zcat ./motifs/NRF1_mm0to2.bed.gz) \
  -b ./snv/all_unique_variants_across_cancers.bed \
  > ./motifs/overlaps/motif2mm_snv_overlaps/NRF1_mm0to2_motifs_with_SNVs.bed
bedtools intersect -u \
  -a <(zcat ./motifs/NRF1_mm0to2_noCGmm.bed.gz) \
  -b ./snv/all_unique_variants_across_cancers.bed \
  > ./motifs/overlaps/motif2mm_snv_overlaps/NRF1_mm0to2_noCGmm_motifs_with_SNVs.bed

# BANP
bedtools intersect -u \
  -a <(zcat ./motifs/BANP_mm0to2.bed.gz) \
  -b ./snv/all_unique_variants_across_cancers.bed \
  > ./motifs/overlaps/motif2mm_snv_overlaps/BANP_mm0to2_motifs_with_SNVs.bed

bedtools intersect -u \
  -a <(zcat ./motifs/BANP_mm0to2_noCGmm.bed.gz) \
  -b ./snv/all_unique_variants_across_cancers.bed \
  > ./motifs/overlaps/motif2mm_snv_overlaps/BANP_mm0to2_noCGmm_motifs_with_SNVs.bed


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

############################################################################
# 11) LOOKING INTO IF MUTATION AFFECT METHYLATION LEVELS AT A CG SITE 
############################################################################
# 1) Get the list of HM450 probes that overlap with NRF1 motif (2 CG in motif)
set -euo pipefail

PROBES_TXT="./methylation/overlaps/intersected_motifs2mm_HM450/NRF1_mm0to2_noCGmm_probes.txt"
ALL_CPG_BED="./methylation/annotated_methylation_data_probes_filtered.bed"
MOTIF_BED="./motifs/NRF1_mm0to2_noCGmm.bed.gz"
OUTDIR="./results/methylation/methylation_motif_mutation_overlap/NRF1_mm0to2_noCGmm_probes"
mkdir -p "$OUTDIR"

# Output: chr  start  end  cg_id  motif_chr  motif_start  motif_end  
bedtools intersect -wa -wb \
  -a <(
        awk 'BEGIN{FS=OFS="\t"}
             NR==FNR {keep[$1]=1; next}
             ($4 in keep) {print $1,$2,$3,$4}' \
          <(grep -v "^[[:space:]]*$" "$PROBES_TXT" | awk '{print $1}' | sort -u) \
          "$ALL_CPG_BED" \
        | sort -k1,1 -k2,2n
      ) \
  -b <(zcat "$MOTIF_BED") \
> "$OUTDIR/NRF1_mm0to2_noCGmm_probes_with_coords_and_motif.tsv"

set -euo pipefail

PROBES_TXT="./methylation/overlaps/intersected_motifs2mm_HM450/BANP_mm0to2_noCGmm_probes.txt"
ALL_CPG_BED="./methylation/annotated_methylation_data_probes_filtered.bed"
MOTIF_BED="./motifs/BANP_mm0to2_noCGmm.bed.gz"
OUTDIR="./results/methylation/methylation_motif_mutation_overlap/BANP_mm0to2_noCGmm_probes"
mkdir -p "$OUTDIR"

# Output: chr  start  end  cg_id  motif_chr  motif_start  motif_end  
bedtools intersect -wa -wb \
  -a <(
        awk 'BEGIN{FS=OFS="\t"}
             NR==FNR {keep[$1]=1; next}
             ($4 in keep) {print $1,$2,$3,$4}' \
          <(grep -v "^[[:space:]]*$" "$PROBES_TXT" | awk '{print $1}' | sort -u) \
          "$ALL_CPG_BED" \
        | sort -k1,1 -k2,2n
      ) \
  -b <(zcat "$MOTIF_BED") \
> "$OUTDIR/BANP_mm0to2_noCGmm_probes_with_coords_and_motif.tsv"

# 2) Get the list of SNVs that overlap with these probes in motifs (using the TSV above) and the VCF files per cancer type (filtered SNVs without structural variants)
#!/usr/bin/env bash
set -euo pipefail
VCF_GLOB="./snv/snv_filtered_without_structural_variants/*.vcf.gz"
TSV="./results/methylation/methylation_motif_mutation_overlap/NRF1_mm0to2_noCGmm_probes/NRF1_mm0to2_noCGmm_probes_with_coords_and_motif.tsv"
OUT_DIR="./results/methylation/methylation_motif_mutation_overlap/vcf_hits_in_NRF1_mm0to2_noCGmm_motifs"
mkdir -p "$OUT_DIR"
TMPDIR="$(mktemp -d)"
trap 'rm -rf "$TMPDIR"' EXIT
# Build a BED with motif intervals FROM TSV (these are the motifs that overlap cg probes) and attach cg as 4th column for bedtools -wb output.
# TSV columns: cg=$4, motif_chr=$5, motif_start=$6, motif_end=$7
TSV_BED_CG="$TMPDIR/tsv_motifs_with_cg.bed"
awk 'BEGIN{FS=OFS="\t"} NR>1 && $4!="" && $5!="" {print $5,$6,$7,$4}' "$TSV" | sort -k1,1 -k2,2n -k3,3n -k4,4 -u > "$TSV_BED_CG"
for vcf in $VCF_GLOB; do
  sample=$(basename "$vcf" .vcf.gz)
  # cancer code from "SNV_TCGA-ACC-...." -> ACC
  cancer=$(echo "$sample" | sed -n 's/^SNV_TCGA-\([A-Za-z0-9]\+\)-.*/\1/p')
  [[ -n "${cancer:-}" ]]
  out_subdir="$OUT_DIR/$cancer"
  mkdir -p "$out_subdir"
  out="$out_subdir/${sample}.tsv"
  echo "[RUN] $cancer / $sample"
  # Variant stream: chrom, start0, end1, chrom, pos, ref, alt
  bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\n' "$vcf" \
    | awk 'BEGIN{OFS="\t"}{print $1,$2-1,$2,$1,$2,$3,$4}' \
    | bedtools intersect -wa -wb \
        -a "$TSV_BED_CG" \
        -b - \
    | awk 'BEGIN{OFS="\t"}{
        # A (from TSV_BED_CG): $1 $2 $3 $4 = motif_chr motif_start motif_end cg
        # B (from variants): starts at $5
        #   $8=var_chr  $9=var_pos  $10=ref  $11=alt
        print $1,$2,$3,$4,$8,$9,$10">"$11
      }' > "$out"
  echo -e "[OK] Wrote $out\t(hits: $(wc -l < "$out"))"
done
echo "[DONE] Outputs in $OUT_DIR/<CANCER>/"

#!/usr/bin/env bash
set -euo pipefail
VCF_GLOB="./snv/snv_filtered_without_structural_variants/*.vcf.gz"
TSV="./results/methylation/methylation_motif_mutation_overlap/BANP_mm0to2_noCGmm_probes/BANP_mm0to2_noCGmm_probes_with_coords_and_motif.tsv"
OUT_DIR="./results/methylation/methylation_motif_mutation_overlap/vcf_hits_in_BANP_mm0to2_noCGmm_motifs"
mkdir -p "$OUT_DIR"
TMPDIR="$(mktemp -d)"
trap 'rm -rf "$TMPDIR"' EXIT
# Build a BED with motif intervals FROM TSV (these are the motifs that overlap cg probes) and attach cg as 4th column for bedtools -wb output.
# TSV columns: cg=$4, motif_chr=$5, motif_start=$6, motif_end=$7
TSV_BED_CG="$TMPDIR/tsv_motifs_with_cg.bed"
awk 'BEGIN{FS=OFS="\t"} NR>1 && $4!="" && $5!="" {print $5,$6,$7,$4}' "$TSV" | sort -k1,1 -k2,2n -k3,3n -k4,4 -u > "$TSV_BED_CG"
for vcf in $VCF_GLOB; do
  sample=$(basename "$vcf" .vcf.gz)
  # cancer code from "SNV_TCGA-ACC-...." -> ACC
  cancer=$(echo "$sample" | sed -n 's/^SNV_TCGA-\([A-Za-z0-9]\+\)-.*/\1/p')
  [[ -n "${cancer:-}" ]]
  out_subdir="$OUT_DIR/$cancer"
  mkdir -p "$out_subdir"
  out="$out_subdir/${sample}.tsv"
  echo "[RUN] $cancer / $sample"
  # Variant stream: chrom, start0, end1, chrom, pos, ref, alt
  bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\n' "$vcf" \
    | awk 'BEGIN{OFS="\t"}{print $1,$2-1,$2,$1,$2,$3,$4}' \
    | bedtools intersect -wa -wb \
        -a "$TSV_BED_CG" \
        -b - \
    | awk 'BEGIN{OFS="\t"}{
        # A (from TSV_BED_CG): $1 $2 $3 $4 = motif_chr motif_start motif_end cg
        # B (from variants): starts at $5
        #   $8=var_chr  $9=var_pos  $10=ref  $11=alt
        print $1,$2,$3,$4,$8,$9,$10">"$11
      }' > "$out"
  echo -e "[OK] Wrote $out\t(hits: $(wc -l < "$out"))"
done
echo "[DONE] Outputs in $OUT_DIR/<CANCER>/"

# 3) For each hit, get the beta values for the mutated sample and a non-mutated sample (if available) from the methylation files, and compile into a TSV for statistical testing.
#!/usr/bin/env bash
set -euo pipefail
HITS_ROOT="./results/methylation/methylation_motif_mutation_overlap/vcf_hits_in_NRF1_mm0to2_noCGmm_motifs"
METH_DIR="./methylation/filtered_methylation"
OUT="./results/methylation/methylation_motif_mutation_overlap/NRF1_one_pair_per_cancer_mut_vs_nomut_with_beta.tsv"
mkdir -p "$(dirname "$OUT")"
TMPDIR="$(mktemp -d)"
trap 'rm -rf "$TMPDIR"' EXIT
echo -e "cancer\tsample_mut\tsample_nomut\tmotif_chr\tmotif_start\tmotif_end\tcg\tvar_chr\tvar_pos\tref_alt\tbeta_nomut\tbeta_mut" > "$OUT"
extract_cancer_from_sample() {
  echo "$1" | sed -n 's/^SNV_TCGA-\([A-Za-z0-9]\+\)-.*/\1/p'
}

extract_tumor_id_01A() {
  local base
  base=$(echo "$1" | grep -oE 'TCGA-[A-Za-z0-9]{2}-[A-Za-z0-9]{4}-01A' | head -n 1 || true)
  [[ -n "${base:-}" ]] && echo "${base}_1" || echo ""
}

find_meth_file() {
  local cancer="$1"
  local tumor_id="$2"
  ls -1 "$METH_DIR"/HM450_TCGA-"$cancer"-"$tumor_id"_annotated_methylation_filtered.bed.gz 2>/dev/null | head -n 1 || true
}

# Build cached cg->beta map for a methylation file (no sorting randomness: keep first seen)
build_beta_map() {
  local meth="$1"
  local key out
  key=$(basename "$meth" .bed.gz)
  out="$TMPDIR/${key}.cg_beta.tsv"
  [[ -s "$out" ]] && { echo "$out"; return; }

  zcat "$meth" | awk 'BEGIN{FS=OFS="\t"}
    {
      # Find cg anywhere in line, then pick the numeric beta immediately after it
      cg=""; beta=""
      for(i=1;i<=NF;i++){
        if($i ~ /^cg[0-9]+$/){ cg=$i; cg_i=i; break }
      }
      if(cg!=""){
        for(j=cg_i+1;j<=NF;j++){
          if($j ~ /^[0-9]*\.?[0-9]+$/){ beta=$j; break }
        }
      }
      if(cg!="" && beta!="" && !seen[cg]++){ print cg,beta }
    }' > "$out"
  echo "$out"
}

lookup_beta() {
  local map="$1"
  local cg="$2"
  awk -v cg="$cg" 'BEGIN{FS="\t"} $1==cg {print $2; exit}' "$map"
}

for cancer_dir in "$HITS_ROOT"/*; do
  [[ -d "$cancer_dir" ]] || continue
  cancer=$(basename "$cancer_dir")
  echo "[CANCER] $cancer"
  mapfile -t files < <(ls -1 "$cancer_dir"/*.tsv 2>/dev/null | sort)
  [[ ${#files[@]} -ge 2 ]] || { echo "  [SKIP] <2 samples"; continue; }
  found=0

  # Try mutated sample files first (non-empty)
  for mut_file in "${files[@]}"; do
    [[ -s "$mut_file" ]] || continue
    sample_mut=$(basename "$mut_file" .tsv)
    # Find methylation file for the mutated sample (01A)
    c_mut=$(extract_cancer_from_sample "$sample_mut")
    id_mut=$(extract_tumor_id_01A "$sample_mut")
    mut_meth=$(find_meth_file "$c_mut" "$id_mut")
    [[ -n "${mut_meth:-}" ]] || continue
    mut_map=$(build_beta_map "$mut_meth")

    # Candidate non-mutated files: simplest definition = empty hit file
    for nomut_file in "${files[@]}"; do
      [[ "$nomut_file" == "$mut_file" ]] && continue
      [[ ! -s "$nomut_file" ]] || continue
      sample_nomut=$(basename "$nomut_file" .tsv)
      c_nomut=$(extract_cancer_from_sample "$sample_nomut")
      id_nomut=$(extract_tumor_id_01A "$sample_nomut")
      nomut_meth=$(find_meth_file "$c_nomut" "$id_nomut")
      [[ -n "${nomut_meth:-}" ]] || continue
      nomut_map=$(build_beta_map "$nomut_meth")

      # Now iterate over hits in mutated file until we find a cg with non-NA beta in both
      while IFS=$'\t' read -r motif_chr motif_start motif_end cg var_chr var_pos ref_alt; do
        [[ -n "${cg:-}" ]] || continue
        beta_mut=$(lookup_beta "$mut_map" "$cg" || true)
        beta_nomut=$(lookup_beta "$nomut_map" "$cg" || true)
        [[ -n "${beta_mut:-}" && -n "${beta_nomut:-}" ]] || continue
        echo -e "${cancer}\t${sample_mut}\t${sample_nomut}\t${motif_chr}\t${motif_start}\t${motif_end}\t${cg}\t${var_chr}\t${var_pos}\t${ref_alt}\t${beta_nomut}\t${beta_mut}" >> "$OUT"
        echo "  [OK] pair chosen: $sample_mut vs $sample_nomut (cg=$cg)"
        found=1
        break
      done < "$mut_file"
      [[ "$found" -eq 1 ]] && break
    done
    [[ "$found" -eq 1 ]] && break
  done
  [[ "$found" -eq 0 ]] && echo "  [WARN] no valid pair with beta found for $cancer"
done

echo "[DONE] Wrote $OUT"

# plot the difference in beta values between mutated and non-mutated samples for the cg sites in the motifs, per cancer type (boxplot or similar)
Rscript -e '
library(readr)
library(dplyr)
library(ggplot2)

df <- read_tsv("./results/methylation/methylation_motif_mutation_overlap/NRF1_one_pair_per_cancer_mut_vs_nomut_with_beta.tsv") %>%
  mutate(
    beta_mut   = as.numeric(beta_mut),
    beta_nomut = as.numeric(beta_nomut),
    delta = beta_mut - beta_nomut
  )

p <- ggplot(df, aes(x = reorder(cancer, delta), y = delta)) +
  geom_hline(yintercept = 0) +
  geom_point() +
  coord_flip() +
  theme_bw() +
  labs(x = NULL, y = "deltabeta = beta_mut - beta_nomut",title = "Difference in DNA methylation at TF binding motifs when a mutation is present")

ggsave("./results/methylation/methylation_motif_mutation_overlap/NRF1_delta_beta_per_cancer.pdf", plot = p, width = 12, height = 12)
'


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

# Load canonical cancer list/order
cancer_order <- scan("./results/multi_omics/cancer_color_order.txt", what = "", quiet = TRUE)

# Same palette as peaks plot
palette_main <- rainbow(length(cancer_order))
names(palette_main) <- cancer_order

# Keep all cancers (unknown ones will be grey)
palette_plot <- c(
  palette_main,
  setNames(
    rep("grey70", length(setdiff(unique(df$Cancer), cancer_order))),
    setdiff(unique(df$Cancer), cancer_order)
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

##################################################################################################
# Boxplot the distribution of gene expression of NRF1 and BANP between healthy and tumor samples -per pairs
##################################################################################################
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

# 3) Plot the distribution of TPM values for NRF1 and BANP in tumor vs healthy samples (boxplot)
Rscript -e '
library(data.table)
library(ggplot2)
IN  <- "./expression/patients_with_both_tumor_and_healthy_samples_with_tpm.tsv"
OUT <- "./results/expression/boxplot_distribution_NRF1_BANP_tumor_vs_healthy_per_cancer_wilcox.pdf"
df <- fread(IN)

# Long format
d <- rbind(
  df[, .(case_id, cancer, gene, group="tumor",   tpm=as.numeric(tumor_tpm))],
  df[, .(case_id, cancer, gene, group="healthy", tpm=as.numeric(normal_tpm))])

d[, group := factor(group, levels=c("healthy","tumor"))]
d[, expr := log2(tpm + 1)]

# Paired table
tum <- d[group=="tumor",   .(case_id, cancer, gene, tumor=expr)]
hea <- d[group=="healthy", .(case_id, cancer, gene, healthy=expr)]
pair <- merge(tum, hea, by=c("case_id","cancer","gene"))

# Wilcoxon paired p-values per (cancer, gene)
pvals <- pair[, .(p = wilcox.test(tumor, healthy, paired=TRUE)$p.value), by=.(cancer, gene)]

# p-value -> stars (blank if not significant)
p_to_stars <- function(p){
  if (is.na(p)) return("NA")
  if (p < 1e-4) return("****")
  if (p < 1e-3) return("***")
  if (p < 1e-2) return("**")
  if (p < 5e-2) return("*")
  return("")}
pvals[, stars := vapply(p, p_to_stars, character(1))]

# Plotting data restricted to paired set only
dplot <- merge(d, pair[, .(case_id, cancer, gene)], by=c("case_id","cancer","gene"))

# Add n per cancer (paired cases) and facet label
n_df <- pair[, .(n=.N), by=.(cancer, gene)]
dplot <- merge(dplot, n_df, by=c("cancer","gene"))
dplot[, cancer_lab := paste0(cancer, " (n=", n, ")")]

# Annotation positions (top of each facet) + SAME facet label column
ymax <- dplot[, .(y=max(expr, na.rm=TRUE)), by=.(cancer, gene, cancer_lab)]
ann  <- merge(pvals, ymax, by=c("cancer","gene"))
ann  <- merge(ann, n_df, by=c("cancer","gene"))
ann[, cancer_lab := paste0(cancer, " (n=", n, ")")]
ann[, `:=`(x=1.5, y=y+0.15)]

# Colors
fill_cols <- c(healthy="#4C78A8", tumor="#E45756")
pdf(OUT, width=14, height=8)

for (g in sort(unique(dplot$gene))) {
  p <- ggplot(dplot[gene==g], aes(x=group, y=expr, fill=group)) +
    geom_violin(trim=FALSE, alpha=0.35, color=NA) +
    geom_boxplot(width=0.22, outlier.shape=NA, alpha=0.85) +
    geom_jitter(width=0.12, size=0.9, alpha=0.45) +
    facet_wrap(~cancer_lab, scales="free_y") +
    geom_text(data=ann[gene==g], aes(x=x, y=y, label=stars), inherit.aes=FALSE, size=4) +
    scale_fill_manual(values=fill_cols) +
    theme_classic(base_size=13) +
    theme(
      legend.position="top",
      strip.background = element_rect(fill="grey95", color=NA),
      strip.text = element_text(face="bold")
    ) +
    labs(
      title = paste0(g, ": Tumor vs Healthy per cancer type (paired Wilcoxon)"),
      subtitle = "Each facet is a cancer type. Each dot is one patient (paired tumor and normal). Boxes show median and IQR; violins show distribution.",
      x = NULL,
      y = "log2(TPM + 1) expression",
      fill = "Sample type",
      caption = "Stars indicate paired Wilcoxon significance (* p<0.05, ** p<0.01, *** p<0.001, **** p<1e-4)."
    )
  print(p)
}
dev.off()
'

#!!IM HERE!!
##################################################################################################################################################
# PCA of the expression of top 10000 variable genes across all samples, to see if they cluster by cancer type or by healthy vs tumor status
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
#!!LAUNCHED!!
#!!LAUNCHED!!
Rscript -e '
library(data.table)
library(ggplot2)

cat("[INFO] Loading required packages...\n")
if (!requireNamespace("irlba", quietly=TRUE)) stop("Install irlba: install.packages(\"irlba\")")
if (!requireNamespace("Rtsne", quietly=TRUE)) stop("Install Rtsne: install.packages(\"Rtsne\")")
if (!requireNamespace("plotly", quietly=TRUE)) stop("Install plotly: install.packages(\"plotly\")")
if (!requireNamespace("htmlwidgets", quietly=TRUE)) stop("Install htmlwidgets: install.packages(\"htmlwidgets\")")

TPM_COL <- 7  # <-- set to the column index of tpm_unstranded in your filtered files

col_file <- "./results/multi_omics/cancer_color_order_with_defined_colours.tsv"
out_html <- "./results/expression/TSNE_TPM_protein_coding_ALL_rawTPM_interactive.html"

cat("[INFO] Loading cancer colors from:", col_file, "\n")
col_dt <- fread(col_file, header=FALSE)
if (ncol(col_dt) < 2) stop("color file must have >=2 columns")
col_dt <- col_dt[, 1:2]
setnames(col_dt, c("cancer","color"))
cancer_levels <- as.character(col_dt$cancer)
cancer_cols   <- setNames(as.character(col_dt$color), as.character(col_dt$cancer))

cat("[INFO] Listing filtered expression files...\n")
files <- list.files("./expression/expression_files",pattern="[.]augmented_star_gene_counts_filtered[.]tsv$",full.names=TRUE)
cat("[INFO] Found", length(files), "files\n")

cat("[INFO] Reading first file (schema reference):\n  ", files[1], "\n", sep="")
dt0 <- fread(files[1], header=FALSE, showProgress=FALSE)
cat("[INFO] Columns detected:", ncol(dt0), " | Rows:", nrow(dt0), "\n")
cat("[INFO] Assuming: V2=gene_name, V3=gene_type, TPM_COL=", TPM_COL, "\n", sep="")

pc_idx <- which(dt0$V3 == "protein_coding")
g0 <- dt0$V2[pc_idx]
G <- length(g0)
S <- length(files)
cat("[INFO] Protein_coding genes:", G, " | Samples(files):", S, "\n")

cat("[INFO] Allocating matrix...\n")
m <- matrix(0, nrow=G, ncol=S)
rownames(m) <- g0

skipped <- 0L
used    <- 0L

cat("[INFO] Reading TPM column across files...\n")
for (i in seq_along(files)) {
  if (i %% 25 == 0) cat("[READ]", i, "/", S, "\n")
  dt <- fread(files[i], header=FALSE, showProgress=FALSE)

  if (nrow(dt) != nrow(dt0)) {
    cat("[SKIP] Row mismatch:", basename(files[i]), " rows=", nrow(dt), " ref=", nrow(dt0), "\n", sep="")
    skipped <- skipped + 1L
    next
  }

  x <- suppressWarnings(as.numeric(dt[[TPM_COL]][pc_idx]))
  if (length(x) != G) {
    cat("[SKIP] protein_coding length mismatch:", basename(files[i]), " len=", length(x), " ref=", G, "\n", sep="")
    skipped <- skipped + 1L
    next
  }

  m[, i] <- x
  used <- used + 1L
}

cat("[INFO] Done reading. Used:", used, " | Skipped:", skipped, "\n")

cat("[INFO] Setting sample names...\n")
sids <- sub("^RNA_", "", sub("[.]augmented_star_gene_counts_filtered[.]tsv$", "", basename(files)))
colnames(m) <- sids

m[is.na(m)] <- 0
xs <- t(m)
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
df$type <- ifelse(substr(code2,1,1)=="0","Tumor",ifelse(substr(code2,1,1)=="1","Normal","Other"))
df$type <- factor(df$type, levels=c("Tumor","Normal","Other"))
df$cancer <- factor(df$cancer, levels=cancer_levels)

# IMPORTANT: do NOT name this column "text" (collides with text() function)
df$hover <- paste0("Sample: ", df$sample,
                   "<br>Cancer: ", as.character(df$cancer),
                   "<br>Type: ", as.character(df$type))

p <- ggplot(df, aes(TSNE1, TSNE2, color=cancer, shape=type)) +
  geom_point(aes(tooltip = hover), size=2, alpha=0.85) +
  scale_shape_manual(values=c(Tumor=16, Normal=17, Other=4)) +
  scale_color_manual(values=cancer_cols, breaks=cancer_levels, drop=FALSE, na.value="grey70") +
  theme_bw() +
  labs(
    title="Interactive t-SNE — protein_coding (no header files)",
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

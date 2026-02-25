# This script processes filtered SNV VCF files to generate presence/absence matrices of SNVs in NRF1 motif regions across different cancer types, and then creates heatmaps to visualize the results.
#!/usr/bin/env bash
set -euo pipefail
# 1. For each filtered SNV VCF file, extract unique variant IDs in the format CHROM:POS:REF:ALT
VCF_DIR="./snv/snv_filtered_without_structural_variants"
OUT_ROOT="./snv/heatmap"
mkdir -p "$OUT_ROOT"

for vcf in "${VCF_DIR}"/SNV_TCGA-*.vcf.gz; do
  bn=$(basename "$vcf" .vcf.gz)

  # extract cancer type from filename: SNV_TCGA-<CANCER>-...
  CANCER=$(echo "$bn" | sed -E 's/^SNV_TCGA-([A-Z0-9]+)-.*$/\1/')

  OUT_DIR="${OUT_ROOT}/${CANCER}"
  mkdir -p "$OUT_DIR"

  echo "[RUNNING] ${bn}  -> ${OUT_DIR}/${bn}.ids"

  bcftools query -f "%CHROM:%POS:%REF:%ALT\n" "$vcf" \
    | awk -F: 'BEGIN{OFS=":"} {c=$1; if(c !~ /^chr/) c="chr"c; $1=c; print}' \
    | sort -u \
    > "${OUT_DIR}/${bn}.ids"

  echo "[DONE] ${bn}  (n=$(wc -l < "${OUT_DIR}/${bn}.ids"))"
done

#2.Filter those not in the NRF1 motif regions
#!/usr/bin/env bash
set -euo pipefail

OUT_ROOT="./snv/heatmap"
MOTIF_DIR="./snv/overlaps/snv_in_motifs2mm"
TF_LIST=("BANP_mm0to2_noCGmm" "BANP_mm0to2" "NRF1_mm0to2_noCGmm" "NRF1_mm0to2")

MAX_JOBS=2
job_count=0

run_one_tf() {
  local TF="$1"
  local MOTIF_IDS="${MOTIF_DIR}/${TF}_SNV_ids.txt"
  [[ -f "$MOTIF_IDS" ]] || { echo "[SKIP] missing motif list: $MOTIF_IDS"; return 0; }
  echo "[CASE] $TF"
  for dir in "${OUT_ROOT}"/*; do
    [[ -d "$dir" ]] || continue

    for id_file in "${dir}"/SNV_TCGA-*.ids; do
      [[ -f "$id_file" ]] || continue
      [[ "$id_file" == *"_motif.ids" ]] && continue
      local bn
      bn=$(basename "$id_file" .ids)
      local OUT_FILE="${dir}/${bn}_in_${TF}_motif.ids"
      grep -F -x -f "$MOTIF_IDS" "$id_file" > "$OUT_FILE" || true
    done
  done
  echo "[DONE] $TF"
}

for TF in "${TF_LIST[@]}"; do
  run_one_tf "$TF" &   # run in background
  job_count=$((job_count+1))

  # if we have MAX_JOBS running, wait for them to finish before launching more
  if [[ "$job_count" -ge "$MAX_JOBS" ]]; then
    wait
    job_count=0
  fi
done

# wait for any remaining jobs
wait
echo "ALL DONE"

# 3. Build presence/absence matrices for each cancer type and TF
#!/usr/bin/env bash
set -euo pipefail
OUT_ROOT="./snv/heatmap"
TF_LIST=("BANP_mm0to2_noCGmm" "BANP_mm0to2" "NRF1_mm0to2_noCGmm" "NRF1_mm0to2")

for CANCER_DIR in "${OUT_ROOT}"/*; do
  [[ -d "$CANCER_DIR" ]] || continue
  CANCER=$(basename "$CANCER_DIR")
  echo "=============================="
  echo "[CANCER] $CANCER"
  echo "=============================="
  for TF in "${TF_LIST[@]}"; do
    echo "[TF] $CANCER / $TF"
    shopt -s nullglob
    files=("${CANCER_DIR}"/SNV_TCGA-*_in_"${TF}"_motif.ids)
    shopt -u nullglob
    if [[ ${#files[@]} -eq 0 ]]; then
      echo "  [SKIP] no *_in_${TF}_motif.ids files"
      continue
    fi

    # 1) Build ROWS
    ROWS="${CANCER_DIR}/ROWS_${CANCER}_in_${TF}_motif.txt"
    cat "${files[@]}" | sort -u > "$ROWS"
    nrows=$(wc -l < "$ROWS" | tr -d ' ')
    if [[ "$nrows" -eq 0 ]]; then
      echo "  [SKIP] empty rows"
      continue
    fi
    echo "  [INFO] rows=$nrows samples=${#files[@]}"

    # 2) Header
    OUT_MATRIX="${CANCER_DIR}/MATRIX_${CANCER}_in_${TF}_motif.tsv"
    {
      printf "SNV_ID"
      for f in "${files[@]}"; do
        bn=$(basename "$f" .ids)
        printf "\t%s" "$bn"
      done
      printf "\n"
    } > "$OUT_MATRIX"

    # 3) Build columns
    tmpdir=$(mktemp -d)
    # ensure tmpdir removed even if a command fails
    cleanup() { rm -rf "$tmpdir"; }
    trap cleanup EXIT
    cp "$ROWS" "${tmpdir}/col0"
    total=${#files[@]}
    coln=1
    for f in "${files[@]}"; do
      awk '
        NR==FNR { has[$1]=1; next }
        { print ($1 in has ? 1 : 0) }
      ' "$f" "$ROWS" > "${tmpdir}/col${coln}"
      coln=$((coln+1))
    done
    # Paste EXACT columns: col0 + col1..colN 
    paste_args=("${tmpdir}/col0")
    for ((i=1; i<=total; i++)); do
      paste_args+=("${tmpdir}/col${i}")
    done
    paste "${paste_args[@]}" >> "$OUT_MATRIX"

    # remove trap for next iteration and cleanup now
    trap - EXIT
    cleanup
    echo "  [DONE] $(basename "$OUT_MATRIX")"
  done
done
echo "ALL MATRICES BUILT"

# 4. Generate heatmaps from the presence/absence matrices
Rscript -e '
library(data.table)
library(pheatmap)

in_root <- "./snv/heatmap"
out_dir <- "./results/snv/heatmap_snv"
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

wrap_label <- function(x, width=35) paste(strwrap(x, width=width), collapse="\n")

TOP_N <- 30 # top N most frequent SNVs to show

for (cdir in list.dirs(in_root, recursive=FALSE, full.names=TRUE)) {
  cancer <- basename(cdir)
  mats <- list.files(cdir, pattern="motif.tsv$", full.names=TRUE)
  if (length(mats) == 0) next

  pdf(file.path(out_dir, paste0("HEATMAP_", cancer, "_01A_vs_10A.pdf")), width=24, height=16)

  for (f in mats) {
    dt <- fread(f)
    snv_ids <- dt[[1]]
    m <- as.matrix(dt[, -1, with=FALSE])
    colnames(m) <- names(dt)[-1]
    rownames(m) <- snv_ids
    storage.mode(m) <- "numeric"
    m[is.na(m)] <- 0

    # keep ONLY columns that are tumor 01A paired with normal 10A and only the first replicate
    keep_cols <- grepl("-01A_vs_", colnames(m)) & grepl("_vs_.*-10A_1_in_", colnames(m)) &
             !grepl("-(01B|10B|10C|11B|11C)", colnames(m)) &
             grepl("_1_in_", colnames(m))
    m <- m[, keep_cols, drop=FALSE]
    if (ncol(m) == 0) next

    # shorten sample column names to: TCGA-OR-A5LN-10A_vs_TCGA-OR-A5LN-01A_1
    cn <- colnames(m)
    cn <- sub("^SNV_TCGA-[A-Z0-9]+-", "", cn)  # remove "SNV_TCGA-ACC-" part
    cn <- sub("_in_.*$", "", cn)              # remove "_in_<TF>_motif" part
    colnames(m) <- cn

    # keep SNVs present at least once within the selected columns
    keep <- rowSums(m) > 0
    m <- m[keep, , drop=FALSE]
    snv_ids <- snv_ids[keep]

    # keep top N most frequent SNVs
    rs <- rowSums(m)
    if (nrow(m) > TOP_N) {
      o <- order(rs, decreasing=TRUE)[1:TOP_N]
      m <- m[o, , drop=FALSE]
      snv_ids <- snv_ids[o]
    }

    # counts that match what is plotted
    n_samples <- ncol(m)
    n_snvs <- nrow(m)

    lab <- vapply(snv_ids, wrap_label, character(1), width=35)

    pheatmap(
      m,
      color = c("skyblue", "yellow"),
      breaks = c(-0.5, 0.5, 1.5),
      labels_row = lab,
      show_rownames = TRUE,
      show_colnames = TRUE,
      fontsize_row = 6,
      fontsize_col = 3.7, 
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      border_color = "grey60",
      legend_breaks = c(0, 1),
      legend_labels = c("absent (0)", "present (1)"),
      main = paste0(cancer, " — ", basename(f)," | 01A_vs_10A pairs=", n_samples," | SNVs(top)=", n_snvs)
    )
  }

  dev.off()
}
'

# nb_cancers : the number of cancer types where that SNV appears in the top 30 most recurrent SNVs (within the filtered columns) for that cancer and that TF.
# This script generates a TSV file with columns: TF, snv_id, nb_cancers, where each row corresponds to a unique SNV that appears in the top 30 most recurrent SNVs for at least one cancer type and one TF. The nb_cancers column indicates how many different cancer types that SNV appears in within the top 30 list for that TF.
Rscript -e '
library(data.table)
in_root <- "./snv/heatmap"
out_dir <- "./snv/top30_panCancer"
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
TF_LIST <- c("BANP_mm0to2_noCGmm","BANP_mm0to2","NRF1_mm0to2_noCGmm","NRF1_mm0to2")
TOP_N <- 30
filter_cols <- function(cols) {
  grepl("-01A_vs_", cols) &
    grepl("_vs_.*-10A_1_in_", cols) &
    !grepl("-(01B|10B|10C|11B|11C)", cols) &
    grepl("_1_in_", cols)
}

all_top <- data.table()
cdirs <- list.dirs(in_root, recursive=FALSE, full.names=TRUE)
for (cdir in cdirs) {
  cancer <- basename(cdir)

  for (TF in TF_LIST) {
    f <- file.path(cdir, sprintf("MATRIX_%s_in_%s_motif.tsv", cancer, TF))
    if (!file.exists(f)) next
    dt <- fread(f)
    if (ncol(dt) < 2) next
    snv_ids <- dt[[1]]
    cols <- names(dt)[-1]
    keep <- filter_cols(cols)
    if (!any(keep)) next
    m <- as.matrix(dt[, cols[keep], with=FALSE])
    storage.mode(m) <- "numeric"
    m[is.na(m)] <- 0
    rs <- rowSums(m)
    ok <- rs > 0
    if (!any(ok)) next
    snv_ids <- snv_ids[ok]
    rs <- rs[ok]
    o <- order(rs, decreasing=TRUE, snv_ids)
    if (length(o) > TOP_N) o <- o[1:TOP_N]
    # store ONLY TF, cancer, snv_id 
    all_top <- rbind(all_top, data.table(TF=TF, cancer=cancer, snv_id=snv_ids[o]))
  }
}

# Count how many cancers each SNV appears in (per TF)
shared <- all_top[, .(nb_cancers = uniqueN(cancer)), by=.(TF, snv_id)]
setorder(shared, TF, -nb_cancers, snv_id)

out_file <- file.path(out_dir, paste0("shared_across_cancers_top", TOP_N, "_nbCancersOnly.tsv"))
fwrite(shared, out_file, sep="\t")
cat("DONE: wrote ", out_file, "\n", sep="")
'
# Plot barplots for the most shared SNVs across cancers, faceted by TF. Each bar is one SNV, and the height is the number of cancer types in which that SNV appears among the top 30 most recurrent SNVs (per cancer) within the motif-filtered set for that TF.
Rscript -e '
library(data.table)
library(ggplot2)

in_file <- "./snv/top30_panCancer/shared_across_cancers_top30_nbCancersOnly.tsv"
out_pdf <- "./results/snv/SNV_shared_top30_cancer_barplots.pdf"

dt <- fread(in_file)

topK <- 30
top <- dt[, head(.SD[order(-nb_cancers, snv_id)], topK), by=TF]

pdf(out_pdf, width=20, height=14)
for(tf in unique(top$TF)){
  d <- top[TF==tf]
  p <- ggplot(d, aes(x=reorder(snv_id, nb_cancers), y=nb_cancers)) +
    geom_col() +
    coord_flip() +
    theme_bw() +
    labs(
      title = paste0(tf, ": most shared Variations across cancers"),
      subtitle = paste0(
        "Each bar is one Variation. The height is the number of cancer types in which this Variant appears ",
        "among the top ", topK, " most frequent Variants (per cancer) within the ", tf, " motif-filtered set."
      ),
      x = "Variant id",
      y = "# cancer types"
    ) +
    theme(
      plot.title = element_text(face="bold",hjust=0.5),
      plot.subtitle = element_text(size=10,hjust=0)
    )
  print(p)
}
dev.off()

cat("DONE -> ", out_pdf, "\n", sep="")
'
# Create a TSV file with columns: cancer, TF, snv_id, n_samples, n_total_samples, where each row corresponds to a unique SNV that appears in the top 30 most recurrent SNVs for that cancer type and that TF. The n_samples column indicates in how many samples (columns) that SNV is present within the filtered columns (tumor 01A paired with normal 10A), and n_total_samples is the total number of such filtered columns for that cancer and TF.
# NA is used for n_samples if the SNV is not present in any sample (after filtering), and NA is used for n_total_samples if there are no columns that pass the filtering criteria for that cancer and TF. 
IN_ROOT=./snv/heatmap
OUT=./snv/top30_panCancer/top30_SNV_nSamples_perCancer_perTF.tsv
SHARED=./snv/top30_panCancer/shared_across_cancers_top30_nbCancersOnly.tsv
TF_LIST=("BANP_mm0to2_noCGmm" "BANP_mm0to2" "NRF1_mm0to2_noCGmm" "NRF1_mm0to2")
TOPK=30
mkdir -p "$(dirname "$OUT")"
echo -e "cancer\tTF\tsnv_id\tn_samples\tn_total_samples" > "$OUT"
for cdir in "$IN_ROOT"/*; do
  cancer=$(basename "$cdir")
  for TF in "${TF_LIST[@]}"; do
    f="$cdir/MATRIX_${cancer}_in_${TF}_motif.tsv"
    [[ -f "$f" ]] || continue
    awk -F'\t' -v OFS='\t' -v cancer="$cancer" -v TF="$TF" -v K="$TOPK" -v S="$SHARED" '
      BEGIN{
        n=0
        while((getline < S)>0){
          if($1==TF){
            n++
            if(n<=K){ top[$2]=1; top_list[n]=$2 }
            else break
          }
        }
        close(S)
      }
      NR==1{
        nk=0
        for(i=2;i<=NF;i++){
          if($i ~ /-01A_vs_/ && $i ~ /_vs_.*-10A_1_in_/ && $i !~ /-(01B|10B|10C|11B|11C)/ && $i ~ /_1_in_/){
            keep[i]=1; nk++
          }
        }
        next
      }
      ($1 in top){
        c=0
        for(i=2;i<=NF;i++) if(keep[i] && ($i+0)>0) c++
        cnt[$1]=c
      }
      END{
        for(i=1;i<=n && i<=K;i++){
          snv=top_list[i]
          val = (snv in cnt) ? cnt[snv] : "NA"
          print cancer, TF, snv, val, nk
        }
      }
    ' "$f" >> "$OUT"
  done
done

# Barplot showing the percentage of samples with each SNV present (after filtering columns to keep only 01A vs 10A pairs), per cancer type, faceted by TF. Each bar is one SNV, and the height is the percentage of samples with that SNV present within the filtered columns for that cancer and TF. The count of samples with the SNV present (n_samples) and the total number of filtered samples (n_total_samples) are also shown as labels on the bars. Cancers types with no samples will not be plotted for that SNV. Only SNVs that are present in at least one sample (after filtering) will be plotted. The pages are ordered by the total number of samples with that SNV present across all cancers across all TF (most shared SNVs first).
Rscript -e '
library(data.table)
library(ggplot2)
in_tsv  <- "./snv/top30_panCancer/top30_SNV_nSamples_perCancer_perTF.tsv"
out_pdf <- "./results/snv/barplot_top30_SNV_nSamples_perCancer_perTF.pdf"
dt <- fread(in_tsv)
dt[, n_samples := suppressWarnings(as.numeric(n_samples))]
dt[, n_total_samples := suppressWarnings(as.numeric(n_total_samples))]
keep_vars <- dt[!is.na(n_samples) & n_samples > 0, unique(snv_id)]
dt <- dt[snv_id %in% keep_vars]
tmp <- dt[, .(tot = sum(n_samples, na.rm=TRUE)), by=snv_id]
setorder(tmp, -tot, snv_id)
var_order <- tmp[["snv_id"]]

pdf(out_pdf, width=16, height=9, useDingbats=FALSE)

for (v in var_order) {
  d <- dt[snv_id == v]
  d <- d[!is.na(n_samples) & n_samples > 0 & !is.na(n_total_samples) & n_total_samples > 0]
  if (nrow(d) == 0) next
  d[, pct := 100 * n_samples / n_total_samples]
  ord <- d[, .(tot = sum(n_samples, na.rm=TRUE)), by=cancer]
  setorder(ord, -tot, cancer)  
  d[, cancer := factor(cancer, levels = ord[["cancer"]])]

  p <- ggplot(d, aes(x=reorder(cancer, pct, FUN=mean), y=pct)) +
    geom_col() +
    coord_flip() +
    facet_wrap(~TF, ncol=2, scales="free_x") +
    theme_bw() +
    labs(
      title = paste0("Variant: ", v),
      subtitle = "Bar height = % of samples with this variant present (after filtering columns to keep only 01A vs 10A pairs), per cancer type.",
      x = "Cancer type",
      y = "% samples with variant"
    ) +
    theme(
      plot.title = element_text(face="bold", hjust=0),
      plot.subtitle = element_text(hjust=0),
      strip.text = element_text(face="bold")
    ) +
    geom_text(aes(label = paste0(n_samples, "/", n_total_samples)),hjust = -0.05, size = 3) +
    expand_limits(y = max(d$pct, na.rm=TRUE) * 1.12)
  print(p)
}
dev.off()
cat("DONE -> ", out_pdf, "\n", sep="")
'

# build  a tsv file with the full list of snv in motifs across all cancer_types and TFs, with the count of samples with that SNV present (after filtering columns to keep only 01A vs 10A pairs) and the total number of filtered samples for that cancer and TF. This file can be used to explore any SNV of interest, not just the most shared ones.
# columns: TF, snv_id, n_samples, cancer, n_total_samples pct_samples_with_snv
# the n_total_samples is the total number of filtered samples in the subset of cancers where this SNV occurs (≥1 time) after filtering columns to keep only 01A vs 10A pairs, for that TF.
Rscript -e '
library(data.table)

in_root <- "./snv/heatmap"
out_dir <- "./snv/top30_panCancer"
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

TF_LIST <- c("BANP_mm0to2_noCGmm","BANP_mm0to2","NRF1_mm0to2_noCGmm","NRF1_mm0to2")

filter_cols <- function(cols) {
  grepl("-01A_vs_", cols) &
    grepl("_vs_.*-10A_1_in_", cols) &
    !grepl("-(01B|10B|10C|11B|11C)", cols) &
    grepl("_1_in_", cols)
}

global <- data.table()
cdirs <- list.dirs(in_root, recursive=FALSE, full.names=TRUE)
for (cdir in cdirs) {
  cancer <- basename(cdir)
  for (TF in TF_LIST) {
    f <- file.path(cdir, sprintf("MATRIX_%s_in_%s_motif.tsv", cancer, TF))
    if (!file.exists(f)) next
    dt <- fread(f)
    if (ncol(dt) < 2) next
    snv_ids <- dt[[1]]
    cols <- names(dt)[-1]
    keep <- filter_cols(cols)
    if (!any(keep)) next
    m <- as.matrix(dt[, cols[keep], with=FALSE])
    storage.mode(m) <- "numeric"
    m[is.na(m)] <- 0
    rs <- rowSums(m)              # number of samples (filtered) with SNV present, within this cancer
    ok <- rs > 0                  # keep SNVs present at least once
    if (!any(ok)) next
    global <- rbind(
      global,
      data.table(TF=TF, snv_id=snv_ids[ok], n_samples=rs[ok], cancer=cancer, n_total_samples=ncol(m))
    )
  }
}

# 1) Global frequency across ALL cancers (per TF, per snv_id)
global_sum <- global[, .(
  total_samples_with_snv = sum(n_samples, na.rm=TRUE),
  cancers_with_snv       = uniqueN(cancer),
  total_samples_scanned  = sum(n_total_samples, na.rm=TRUE)
), by=.(TF, snv_id)]

global_sum[, prevalence_pct := 100 * total_samples_with_snv / total_samples_scanned]
setorder(global_sum, TF, -total_samples_with_snv, -cancers_with_snv, snv_id)

out1 <- file.path(out_dir, "global_motifSNV_frequency_byTF.tsv")
fwrite(global_sum, out1, sep="\t")
cat("DONE\\n- ", out1, "\\n- ", sep="")
'
# NOTE I DIDNT GET THE SAME NUMBER AS THE COUNT OF SNVS IN MOTIFS BECAUSE THE COUNT IN MOTIFS TAKEcat S ALL SAMPLES REGARDLESS IF THEY ARE 10A/01A PAIRED OR NOT, WHEREAS THIS COUNT IS ONLY FOR THE FILTERED SET OF SAMPLES (10A/01A PAIRED). SO SOME SNVS THAT ARE PRESENT IN THE FULL SET OF SAMPLES MIGHT NOT BE PRESENT IN THE FILTERED SET, HENCE THE DISCREPANCY.

# Cluster SNV based on the motif they affect:
#!/usr/bin/env bash
set -euo pipefail

IN="./snv/top30_panCancer/global_motifSNV_frequency_byTF.tsv"
OUT="./snv/top30_panCancer/global_motifSNV_frequency_byTF_clustered.tsv"
DIR="./snv/overlaps/snv_in_motifs2mm"

awk -F'\t' -v OFS='\t' -v DIR="$DIR" '
BEGIN{
  # build mapping TF+snv_id -> motif_sequence (first one seen)
  split("BANP_mm0to2_noCGmm BANP_mm0to2 NRF1_mm0to2_noCGmm NRF1_mm0to2", TFs, " ")
  for (i in TFs) {
    tf = TFs[i]
    f = DIR "/" tf "_SNVs_in_motifs_SNVidxmotif.bed"
    while ((getline line < f) > 0) {
      n = split(line, a, "\t")
      snv = a[10]      # snv_id
      motif = a[4]     # motif sequence
      key = tf "\t" snv
      if (!(key in map)) map[key] = motif
    }
    close(f)
  }
}
NR==1 { print $0, "motif"; next }
{
  key = $1 "\t" $2
  print $0, (key in map ? map[key] : "NA")
}
' "$IN" > "$OUT"

echo "[DONE] $OUT"

########################################################
# To combine snv based on the motif they affect.
########################################################
Rscript -e '
library(data.table)

# --- 1. Configuration ---
in_root  <- "./snv/heatmap"
map_dir  <- "./snv/overlaps/snv_in_motifs2mm"
gene_dir <- "./motifs"
out_dir  <- "./snv/top30_panCancer"
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

TF_LIST <- c("BANP_mm0to2_noCGmm", "BANP_mm0to2", "NRF1_mm0to2_noCGmm", "NRF1_mm0to2")

filter_cols <- function(cols) {
  grepl("-01A_vs_", cols) &
    grepl("_vs_.*-10A_1_in_", cols) &
    !grepl("-(01B|10B|10C|11B|11C)", cols) &
    grepl("_1_in_", cols)
}

# --- 2. Loading SNV->motif-instance map ---
cat("[1/5] Loading SNV-to-motif maps...\n")
snv2motif <- rbindlist(lapply(TF_LIST, function(tf_name){
  f <- file.path(map_dir, paste0(tf_name, "_SNVs_in_motifs_SNVidxmotif.bed"))
  if (!file.exists(f)) return(NULL)
  x <- fread(f, header=FALSE)
  m <- unique(x[, .(
    TF=tf_name,
    snv_id=V10,
    motif_chr=V1,
    motif_start_bed=as.integer(V2),
    motif_end_bed=as.integer(V3),
    motif_sequence=V4,
    motif_strand=V6
  )])
  m[, motif_instance_id := paste(motif_chr, motif_start_bed, motif_end_bed, motif_strand, sep=":")]
  m
}), fill=TRUE)

# --- 3. Loading Gene Map Lookup ---
cat("[2/5] Loading Gene name lookup...\n")
gene_lookup <- rbindlist(lapply(TF_LIST, function(tf_name){
  f <- file.path(gene_dir, paste0(tf_name, "_closest_genes.bed"))
  if (!file.exists(f)) return(NULL)
  g <- fread(f, select=c(1,2,3,6,10), header=FALSE)
  setnames(g, c("V1","V2","V3","V6","V10"), c("chr","s","e","st","gene"))
  g[, .(TF=tf_name, motif_instance_id=paste(chr,s,e,st,sep=":"), gene)]
}), fill=TRUE)

# --- 4. Main Processing Loop ---
cdirs <- list.dirs(in_root, recursive=FALSE, full.names=TRUE)
results <- list()

# Store denominators: for each cancer and TF, how many kept sample columns exist
cancer_counts_list <- list()

for (TF_val in TF_LIST) {
  cat("[TF]", TF_val, "\n")
  map_tf <- snv2motif[TF == TF_val]
  if (nrow(map_tf) == 0) next

  tf_hits_list <- list()
  k <- 0L

  for (cdir in cdirs) {
    cancer <- basename(cdir)
    f <- file.path(cdir, sprintf("MATRIX_%s_in_%s_motif.tsv", cancer, TF_val))
    if (!file.exists(f)) next

    dt <- fread(f)
    if (ncol(dt) < 2) next

    all_cols <- names(dt)
    snv_col <- all_cols[1]
    sample_cols <- all_cols[-1]
    keep_cols <- sample_cols[filter_cols(sample_cols)]
    if (length(keep_cols) == 0) next

    # Denominator (per cancer, per TF): total kept samples
    if (!(cancer %in% names(cancer_counts_list))) cancer_counts_list[[cancer]] <- list()
    cancer_counts_list[[cancer]][[TF_val]] <- length(keep_cols)

    # Long format + keep present
    dt_long <- melt(dt,
                    id.vars=snv_col,
                    measure.vars=keep_cols,
                    variable.name="sample",
                    value.name="score")

    present <- dt_long[score > 0 & get(snv_col) %in% map_tf$snv_id]
    if (nrow(present) == 0) next

    setnames(present, snv_col, "snv_id")
    present[, sample_uid := paste(cancer, sample, sep="|")]

    # Join with motif metadata (snv_id -> motif_instance_id etc.)
    present <- merge(present, map_tf, by="snv_id", all.x=TRUE, allow.cartesian=TRUE)

    k <- k + 1L
    tf_hits_list[[k]] <- unique(present[, .(TF=TF_val, motif_instance_id, sample_uid, cancer, snv_id)])
  }

  if (k == 0) next
  tf_hits_all <- rbindlist(tf_hits_list, use.names=TRUE, fill=TRUE)

  # ---- Per-motif overall metrics (pan-cancer) ----
  summary_motif <- tf_hits_all[, .(
    samples_with_motif_snv        = uniqueN(sample_uid),
    cancers_with_motif_snv        = uniqueN(cancer),
    unique_snvs_mapped_to_motif   = uniqueN(snv_id)
  ), by=.(TF, motif_instance_id)]

  # Attach coordinates and sequence
  meta <- unique(snv2motif[TF == TF_val,
                           .(TF, motif_instance_id, motif_chr, motif_start_bed, motif_end_bed, motif_strand, motif_sequence)])
  summary_motif <- merge(summary_motif, meta, by=c("TF", "motif_instance_id"), all.x=TRUE)

  # ---- FIXED: per-cancer numerator = # unique samples in that cancer with motif-SNV ----
  cancer_num <- tf_hits_all[, .(
    samples_with_motif_snv_in_cancer = uniqueN(sample_uid)
  ), by=.(TF, motif_instance_id, cancer)]

  # Wide: one column per cancer with sample counts
  wide_counts <- dcast(
    cancer_num,
    TF + motif_instance_id ~ cancer,
    value.var = "samples_with_motif_snv_in_cancer",
    fill = 0
  )

  summary_motif <- merge(summary_motif, wide_counts, by=c("TF", "motif_instance_id"), all.x=TRUE)

  # Convert raw per-cancer counts into per-cancer percentages using per-cancer denominators
  cancer_names <- setdiff(names(summary_motif),
                          c("TF","motif_instance_id","samples_with_motif_snv","cancers_with_motif_snv",
                            "unique_snvs_mapped_to_motif","motif_chr","motif_start_bed","motif_end_bed",
                            "motif_strand","motif_sequence"))

  for (cn in cancer_names) {
    denom <- NULL
    if (!is.null(cancer_counts_list[[cn]]) && !is.null(cancer_counts_list[[cn]][[TF_val]])) {
      denom <- cancer_counts_list[[cn]][[TF_val]]
    }
    summary_motif[[paste0(cn, "_pct")]] <- if (!is.null(denom) && denom > 0) {
      100 * (summary_motif[[cn]] / denom)
    } else {
      NA_real_
    }
    summary_motif[[cn]] <- NULL
  }

  results[[TF_val]] <- summary_motif
}

# --- 5. Aggregating and Gene-level info ---
cat("[3/5] Final aggregation...\n")
final_motif_table <- rbindlist(results, fill=TRUE, use.names=TRUE)

# Join gene names
final_motif_table <- merge(final_motif_table, gene_lookup, by=c("TF", "motif_instance_id"), all.x=TRUE)

# --- 6. Formatting & Writing ---
# Put key columns first (keep everything else after)
first_cols <- c("TF", "gene", "motif_instance_id", "motif_chr", "motif_start_bed", "motif_end_bed", "motif_strand", "motif_sequence",
                "samples_with_motif_snv", "cancers_with_motif_snv", "unique_snvs_mapped_to_motif")
first_cols <- intersect(first_cols, names(final_motif_table))
setcolorder(final_motif_table, c(first_cols, setdiff(names(final_motif_table), first_cols)))

setorder(final_motif_table, -samples_with_motif_snv)

out_file <- file.path(out_dir, "motif_instance_level_panCancer_wide.tsv")
fwrite(final_motif_table, out_file, sep="\\t")
cat("[4/5] DONE -> ", out_file, "\\n", sep="")
'
##!!NOT CORRRECT YET IM HERE!!
##############################
# plot heatmap
##############################
Rscript -e '
library(data.table)
library(ggplot2)

in_file  <- "./snv/top30_panCancer/motif_instance_level_panCancer_wide.tsv"
out_pdf  <- "./results/snv/heatmap_gene_by_cancer_perTF.pdf"
dir.create(dirname(out_pdf), recursive=TRUE, showWarnings=FALSE)

dt <- fread(in_file)
if (!("gene" %in% names(dt))) stop("Your wide TSV has no gene column.")

# --- wide -> long for *_pct columns
pct_cols <- grep("_pct$", names(dt), value=TRUE)
if (length(pct_cols) == 0) stop("No *_pct columns found (expected columns like BRCA_pct).")

# --- CLEAN gene names (critical for true combining)
dt[, gene := trimws(gene)]
dt[gene %in% c("", ".", "NA", "NaN"), gene := NA_character_]
dt <- dt[!is.na(gene)]

# --- long format
long <- melt(
  dt,
  id.vars = c("TF","gene"),
  measure.vars = pct_cols,
  variable.name = "cancer",
  value.name = "pct_samples"
)

long[, cancer := sub("_pct$", "", cancer)]
long[, pct_samples := suppressWarnings(as.numeric(pct_samples))]

# --- COMBINE motif instances per same gene (per TF, cancer)
# using MAX: "best motif instance hit for that gene in that cancer"
agg <- long[, .(pct_samples = suppressWarnings(max(pct_samples, na.rm=TRUE))),
            by=.(TF, cancer, gene)]
agg[!is.finite(pct_samples), pct_samples := NA_real_]

# --- Sanity checks to prove combining happened
cat("\n[Check] Duplicates per (TF,cancer,gene) AFTER aggregation (should be none):\n")
dup <- agg[, .N, by=.(TF,cancer,gene)][N > 1]
print(dup)

cat("\n#genes per TF (after combining):\n")
print(agg[, .(n_genes = uniqueN(gene)), by=TF][order(TF)])

# --- Gene order per TF (no filtering)
# Order genes by their overall max across cancers (within each TF)
gene_order <- agg[, .(mx = max(pct_samples, na.rm=TRUE)), by=.(TF, gene)]
gene_order[!is.finite(mx), mx := NA_real_]
setorder(gene_order, TF, -mx, gene)

pdf(out_pdf, width=100, height=10, useDingbats=FALSE)

for (tf in unique(agg$TF)) {
  d <- agg[TF == tf]

  # apply ordering (all genes, no top selection)
  go <- gene_order[TF == tf & !is.na(mx), gene]
  if (length(go) > 0) d[, gene := factor(gene, levels=go)]

  # cancer order alphabetical
  d[, cancer := factor(cancer, levels=sort(unique(cancer)))]

  p <- ggplot(d, aes(x=gene, y=cancer, fill=pct_samples)) +
    geom_tile() +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=2),
      axis.text.y = element_text(size=8),
      plot.title  = element_text(face="bold")
    ) +
    labs(
      title = paste0(tf, " — % samples with ≥1 motif-SNV mapped to gene "),
      x = "Gene",
      y = "Cancer type",
      fill = "% samples"
    )

  print(p)
}

dev.off()
cat("DONE -> ", out_pdf, "\n", sep="")
'

############################################################################################################################################
# get gene list for each TF with at least one SNV in motif present in at least one sample (after filtering)
#############################################################################################################################################
cat ./snv/top30_panCancer/motif_instance_level_panCancer_wide.tsv | grep -w "NRF1_mm0to2_noCGmm" | awk '{print $2}' | sort -u > ./snv/top30_panCancer/genes_with_SNV_in_motif_NRF1_mm0to2_noCGmm.txt
cat ./snv/top30_panCancer/motif_instance_level_panCancer_wide.tsv | grep -w "NRF1_mm0to2" | awk '{print $2}' | sort -u > ./snv/top30_panCancer/genes_with_SNV_in_motif_NRF1_mm0to2.txt
cat ./snv/top30_panCancer/motif_instance_level_panCancer_wide.tsv | grep -w "BANP_mm0to2_noCGmm" | awk '{print $2}' | sort -u > ./snv/top30_panCancer/genes_with_SNV_in_motif_BANP_mm0to2_noCGmm.txt
cat ./snv/top30_panCancer/motif_instance_level_panCancer_wide.tsv | grep -w "BANP_mm0to2" | awk '{print $2}' | sort -u > ./snv/top30_panCancer/genes_with_SNV_in_motif_BANP_mm0to2.txt

############################################################################################################################################
# Look into expression of the genes with SNV in motif
#############################################################################################################################################
# generate matrix of expression for the genes with SNV in motif, for the same samples (10A/01A pairs) used in the SNV heatmap, and plot heatmap of expression for those genes across those samples, faceted by TF. by filtering the matrix generated by gene expression
Rscript -e '
library(data.table)

mat_file <- "./expression/gene_expression_matrix_protein_coding.tsv.gz"
gene_dir <- "./snv/top30_panCancer"
out_dir  <- "./results/expression_SNV_motif_genes_01A11A_from_BIGMATRIX"
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

TF_LIST <- c("BANP_mm0to2_noCGmm","BANP_mm0to2","NRF1_mm0to2_noCGmm","NRF1_mm0to2")

cat("[1/4] Loading big expression matrix...\n")
X <- fread(mat_file)  # columns: sample + many genes
if (!("sample" %in% names(X))) stop("Expected first column named sample in big matrix (keep.rownames=\"sample\").")

cat("[2/4] Filtering samples to -01A or -11A...\n")
X <- X[grepl("-(01A|11A)$", sample)]
cat("[INFO] Samples kept:", nrow(X), "\n")

# load gene lists
cat("[3/4] Loading gene lists...\n")
gene_sets <- lapply(TF_LIST, function(tf){
  f <- file.path(gene_dir, paste0("genes_with_SNV_in_motif_", tf, ".txt"))
  if (!file.exists(f)) stop("Missing gene list: ", f)
  g <- fread(f, header=FALSE)[[1]]
  g <- trimws(g)
  g <- g[g != "" & g != "gene" & g != "." & !is.na(g)]
  unique(g)
})
names(gene_sets) <- TF_LIST

all_genes <- sort(unique(unlist(gene_sets)))

# which genes exist in the big matrix?
available <- intersect(all_genes, setdiff(names(X), "sample"))
missing   <- setdiff(all_genes, available)
cat("[INFO] Genes requested:", length(all_genes),
    " | available in matrix:", length(available),
    " | missing:", length(missing), "\n")

if (length(missing) > 0) {
  fwrite(data.table(missing_gene=missing),
         file.path(out_dir, "genes_missing_from_big_matrix.tsv"),
         sep="\t")
  cat("[WARN] Missing genes written to genes_missing_from_big_matrix.tsv\n")
}

# write union matrix
cat("[4/4] Writing matrices...\n")
X_union <- X[, c("sample", available), with=FALSE]
fwrite(X_union, file.path(out_dir, "TPM_matrix_unionGenes_01A11A.tsv"), sep="\t")

# per TF matrices
for (tf in TF_LIST) {
  gtf <- intersect(gene_sets[[tf]], setdiff(names(X), "sample"))
  out <- file.path(out_dir, paste0("TPM_matrix_", tf, "_01A11A.tsv"))
  fwrite(X[, c("sample", gtf), with=FALSE], out, sep="\t")
  cat("WROTE -> ", out, " | genes=", length(gtf), " | samples=", nrow(X), "\n", sep="")
}

cat("DONE -> ", out_dir, "\n", sep="")
'

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
out_dir <- "./results/heatmap_snv"
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
out_pdf <- "./results/summarys/SNV_shared_top30_cancer_barplots.pdf"

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
#!/usr/bin/env bash
set -euo pipefail
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
out_pdf <- "./results/summarys/barplot_top30_SNV_nSamples_perCancer_perTF.pdf"
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

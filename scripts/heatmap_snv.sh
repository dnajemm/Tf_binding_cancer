#!!IM HERE LAUNCHED CHECK OUTPUT!!
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

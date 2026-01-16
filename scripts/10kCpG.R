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
    pattern = paste0("^HM450_TCGA-", cancer,"-.*_annotated_methylation_filtered\\.bed\\.gz$"),
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

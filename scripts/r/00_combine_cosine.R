#!/usr/bin/env Rscript
## Combine per-tool cosine TSVs and run smoke QC checks.

suppressPackageStartupMessages({
  library(jsonlite)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
source(file.path(script_dir, "helpers.R"))

root <- .mywork_root()
cfg <- load_pipeline_config(root)
cos_dir <- file.path(root, "results", "cosine")
dir.create(cos_dir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(cos_dir, pattern = "^CosineResults_.*_ALL\\.tsv$", full.names = TRUE)
files <- files[!grepl("CosineResults_AllTools_AllCOSMIC\\.tsv$", files)]

if (!length(files)) {
  stop("No CosineResults_*_ALL.tsv files found in ", cos_dir)
}

dfs <- lapply(files, function(f) {
  read.table(f, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
})
all_df <- do.call(rbind, dfs)
rownames(all_df) <- NULL

required <- c("Sample", "Cosine_Similarity", "Dataset", "Data_Type", "Tool", "COSMIC_Version", "Mutation_Count")
missing_cols <- setdiff(required, names(all_df))
if (length(missing_cols)) stop("Missing columns: ", paste(missing_cols, collapse = ", "))

out_path <- file.path(cos_dir, "CosineResults_AllTools_AllCOSMIC.tsv")
write.table(all_df, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
message("Wrote ", out_path, " (", nrow(all_df), " rows)")

# ---- QC ----
qc <- list()
qc$n_rows <- nrow(all_df)
qc$tools <- sort(unique(all_df$Tool))
qc$cosmic <- sort(unique(all_df$COSMIC_Version))
qc$datasets <- sort(unique(all_df$Dataset))

sample_counts <- aggregate(Sample ~ Dataset + Tool + COSMIC_Version, data = all_df, FUN = length)
names(sample_counts)[4] <- "n"
qc$sample_count_table <- sample_counts

# Expected n per dataset
expected <- sapply(names(cfg$datasets), function(d) cfg$datasets[[d]]$expected_n_samples)
mismatches <- sample_counts[sample_counts$n != expected[as.character(sample_counts$Dataset)], , drop = FALSE]
qc$sample_count_mismatches <- mismatches

# Contribution column-sum check (spot-check existing normalized files)
contrib_root <- file.path(root, "results", "contributions")
norm_files <- list.files(contrib_root, pattern = "_contribution_normalized\\.tsv$", recursive = TRUE, full.names = TRUE)
colsum_issues <- list()
for (nf in norm_files) {
  m <- as.matrix(read.table(nf, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE))
  cs <- colSums(m)
  bad <- which(abs(cs - 1) > 1e-3 & cs > 0)
  if (length(bad)) {
    colsum_issues[[nf]] <- unname(cs[bad][seq_len(min(5, length(bad)))])
  }
}
qc$normalized_colsum_issues <- length(colsum_issues)
qc$example_colsum_issues <- head(colsum_issues, 5)

qc_path <- file.path(cos_dir, "cosine_qc_summary.json")
write_json(qc, qc_path, pretty = TRUE, auto_unbox = TRUE, dataframe = "rows")
message("QC summary: ", qc_path)
if (nrow(mismatches)) {
  message("WARNING: sample-count mismatches vs expected:")
  print(mismatches)
} else {
  message("Sample counts match expected sizes for all Tool×COSMIC×Dataset blocks present.")
}
if (length(colsum_issues)) {
  message("WARNING: ", length(colsum_issues), " normalized contribution files have column sums != 1")
}

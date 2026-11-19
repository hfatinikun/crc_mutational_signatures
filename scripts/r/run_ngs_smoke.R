#!/usr/bin/env Rscript
## Phase 2 smoke: MutationalPatterns × COSMIC v3.2 on WGS + TSO-500.
## Checks cosine median >= 0.9 and positive ΔMMR on both data types.

suppressPackageStartupMessages({
  library(jsonlite)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
source(file.path(script_dir, "helpers.R"))

root <- .mywork_root()
cfg <- load_pipeline_config(root)
groups <- load_signature_groups(root)

meta_path <- file.path(root, "metadata", "sample_phenotypes.tsv")
if (!file.exists(meta_path)) {
  system2("Rscript", file.path(script_dir, "00_build_metadata.R"))
}

message("Running Phase-2 NGS smoke fit (MP × v3.2 × WGS+TSO)...")
# Use relative script path from mywork root so spaces in absolute paths cannot break system2
old_wd <- getwd()
setwd(root)
on.exit(setwd(old_wd), add = TRUE)
fit_status <- system2(
  "Rscript",
  args = c(
    "scripts/r/run_r_fitting.R",
    "--tools", "MutationalPatterns",
    "--cosmic", "v3.2",
    "--phase", "2"
  ),
  stdout = "results/logs/ngs_smoke_stdout.txt",
  stderr = "results/logs/ngs_smoke_stderr.txt"
)
if (!identical(as.integer(fit_status), 0L)) {
  message("Fitting failed with status ", fit_status)
  quit(status = 1)
}

pheno <- read.table(meta_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
msi <- pheno$Sample[pheno$Dataset == "Preclinical" & pheno$MSI_status == "MSI"]
mss <- pheno$Sample[pheno$Dataset == "Preclinical" & pheno$MSI_status == "MSS"]

checks <- list()
ok_all <- TRUE
for (ds in c("Preclinical_WGS", "Preclinical_TSO-500")) {
  dt <- cfg$datasets[[ds]]$data_type
  cos_file <- file.path(root, "results", "cosine", "CosineResults_MutationalPatterns_ALL.tsv")
  cos_df <- read.table(cos_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  sub <- cos_df[cos_df$Dataset == ds & cos_df$COSMIC_Version == "v3.2", ]
  med <- median(sub$Cosine_Similarity, na.rm = TRUE)
  frac <- mean(sub$Cosine_Similarity >= 0.9, na.rm = TRUE)

  norm_path <- file.path(root, "results", "contributions", "MutationalPatterns", "v3.2",
                         paste0(ds, "_contribution_normalized.tsv"))
  mat <- as.matrix(read.table(norm_path, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE))
  mmr_sigs <- resolve_signature_group("MMR", "v3.2", rownames(mat), groups)
  mmr <- colSums(mat[mmr_sigs, , drop = FALSE])
  delta <- median(mmr[intersect(msi, names(mmr))], na.rm = TRUE) -
    median(mmr[intersect(mss, names(mmr))], na.rm = TRUE)

  pass_cos <- is.finite(med) && med >= 0.9
  pass_bio <- is.finite(delta) && delta > 0
  checks[[ds]] <- list(
    data_type = dt,
    n = nrow(sub),
    median_cosine = med,
    frac_assigned = frac,
    delta_MMR = delta,
    pass_cosine = pass_cos,
    pass_delta_MMR = pass_bio
  )
  ok_all <- ok_all && pass_cos && pass_bio
  message(sprintf("[%s/%s] n=%d med=%.4f frac=%.3f dMMR=%.4f cos_ok=%s bio_ok=%s",
                  ds, dt, nrow(sub), med, frac, delta, pass_cos, pass_bio))
}

dir.create(file.path(root, "results", "downstream", "ngs_workflow"), recursive = TRUE, showWarnings = FALSE)
write_json(list(passed = ok_all, checks = checks),
           file.path(root, "results", "downstream", "ngs_workflow", "ngs_smoke_validation.json"),
           pretty = TRUE, auto_unbox = TRUE)

if (!ok_all) {
  message("NGS SMOKE VALIDATION FAILED")
  quit(status = 1)
}
message("NGS SMOKE VALIDATION PASSED")

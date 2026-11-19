#!/usr/bin/env Rscript
## Smoke validation: MutationalPatterns × v2/v3.2 on Preclinical WES.
## Checks cosine median >= 0.9 and positive ΔMMR (MSI > MSS).

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

# Ensure metadata exists
meta_path <- file.path(root, "metadata", "sample_phenotypes.tsv")
if (!file.exists(meta_path)) {
  message("Building metadata...")
  system2("Rscript", file.path(script_dir, "00_build_metadata.R"))
}

# Run fitting smoke (inherit working directory)
message("Running MutationalPatterns smoke fit...")
fit_status <- system2(
  "Rscript",
  c(file.path(script_dir, "run_r_fitting.R"),
    "--tools", "MutationalPatterns",
    "--cosmic", "v2,v3.2",
    "--datasets", "Preclinical"),
  stdout = file.path(root, "results", "logs", "mp_smoke_fit_stdout.txt"),
  stderr = file.path(root, "results", "logs", "mp_smoke_fit_stderr.txt")
)
if (!identical(as.integer(fit_status), 0L)) {
  message("Fitting subprocess failed with status ", fit_status)
  message("See results/logs/mp_smoke_fit_stdout.txt and mp_smoke_fit_stderr.txt")
  quit(status = 1)
}

pheno <- read.table(meta_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
msi <- pheno$Sample[pheno$Dataset == "Preclinical" & pheno$MSI_status == "MSI"]
mss <- pheno$Sample[pheno$Dataset == "Preclinical" & pheno$MSI_status == "MSS"]

checks <- list()
ok_all <- TRUE

for (ck in c("v2", "v3.2")) {
  cos_file <- file.path(root, "results", "cosine", "CosineResults_MutationalPatterns_ALL.tsv")
  cos_df <- read.table(cos_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  sub <- cos_df[cos_df$Dataset == "Preclinical" & cos_df$COSMIC_Version == ck, ]
  med <- median(sub$Cosine_Similarity, na.rm = TRUE)
  frac <- mean(sub$Cosine_Similarity >= 0.9, na.rm = TRUE)

  norm_path <- file.path(root, "results", "contributions", "MutationalPatterns", ck,
                         "Preclinical_contribution_normalized.tsv")
  mat <- as.matrix(read.table(norm_path, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE))
  mmr_sigs <- resolve_signature_group("MMR", ck, rownames(mat), groups)
  mmr <- if (length(mmr_sigs) == 1) as.numeric(mat[mmr_sigs, ]) else colSums(mat[mmr_sigs, , drop = FALSE])
  names(mmr) <- colnames(mat)
  delta <- median(mmr[intersect(msi, names(mmr))], na.rm = TRUE) -
    median(mmr[intersect(mss, names(mmr))], na.rm = TRUE)

  pass_cos <- is.finite(med) && med >= 0.9
  pass_bio <- is.finite(delta) && delta > 0
  checks[[ck]] <- list(
    median_cosine = med,
    frac_assigned = frac,
    delta_MMR = delta,
    pass_cosine = pass_cos,
    pass_delta_MMR = pass_bio
  )
  ok_all <- ok_all && pass_cos && pass_bio
  message(sprintf("[%s] median_cosine=%.4f frac>=0.9=%.3f delta_MMR=%.4f  cosine_ok=%s bio_ok=%s",
                  ck, med, frac, delta, pass_cos, pass_bio))
}

out <- list(passed = ok_all, checks = checks)
dir.create(file.path(root, "results", "downstream"), recursive = TRUE, showWarnings = FALSE)
write_json(out, file.path(root, "results", "downstream", "mp_smoke_validation.json"),
           pretty = TRUE, auto_unbox = TRUE)

if (!ok_all) {
  message("SMOKE VALIDATION FAILED")
  quit(status = 1)
}
message("SMOKE VALIDATION PASSED")

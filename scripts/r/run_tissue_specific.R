#!/usr/bin/env Rscript
## Phase-3 expansion: CRC tissue-specific reference (paper Reference analysis).
## Fits MutationalPatterns on Preclinical + Clinical WES using COSMIC key "TS".
##
## Usage (from mywork/):
##   Rscript scripts/r/run_tissue_specific.R

# Paths with spaces break system2(absolute); run relative from mywork cwd.
status <- system2(
  "Rscript",
  c(
    "scripts/r/run_r_fitting.R",
    "--tools", "MutationalPatterns",
    "--cosmic", "TS",
    "--datasets", "Preclinical,Clinical"
  ),
  stdout = "", stderr = ""
)
if (!identical(as.integer(status), 0L)) stop("Tissue-specific fitting failed with status ", status)

message("Tissue-specific (TS) fitting done. Next:")
message("  Rscript scripts/r/00_combine_cosine.R")
message("  Rscript scripts/r/downstream/run_downstream.R")
message("  Rscript scripts/r/downstream/run_reference_expansion.R")

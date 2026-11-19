#!/usr/bin/env Rscript
## Build MSI/MSS/POLE sample lists from Table S1 into metadata/

suppressPackageStartupMessages({
  library(readxl)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
source(file.path(script_dir, "helpers.R"))

root <- .mywork_root()
cfg <- load_pipeline_config(root)
aliases <- load_sample_aliases(root)
# drop comment keys
aliases <- aliases[names(aliases) != "_comment"]

xlsx <- file.path(root, cfg$paths$input_dir, "tables1_bbae249.xlsx")
out_dir <- file.path(root, cfg$paths$metadata_dir)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Preclinical ----
pre_cat <- load_catalogue("Preclinical", cfg, root)
df1 <- as.data.frame(read_excel(xlsx, sheet = 1, skip = 1))
names(df1)[1:3] <- c("Sample_table", "MSI_status", "POLE_status")
df1$Sample <- vapply(df1$Sample_table, function(s) {
  if (!is.null(aliases[[s]])) aliases[[s]] else s
}, character(1))

missing_pre <- setdiff(df1$Sample, colnames(pre_cat))
if (length(missing_pre)) {
  warning("Preclinical samples in Table S1 not in matrix: ", paste(missing_pre, collapse = ", "))
}
df1 <- df1[df1$Sample %in% colnames(pre_cat), , drop = FALSE]

msi_pre <- df1$Sample[df1$MSI_status == "MSI"]
mss_pre <- df1$Sample[df1$MSI_status == "MSS"]
pole_mut_pre <- df1$Sample[df1$POLE_status == "MUT" & df1$MSI_status == "MSS"]
pole_wt_mss_pre <- setdiff(mss_pre, pole_mut_pre)

writeLines(msi_pre, file.path(out_dir, "preclinical_MSI.txt"))
writeLines(mss_pre, file.path(out_dir, "preclinical_MSS.txt"))
writeLines(pole_mut_pre, file.path(out_dir, "preclinical_POLE_MUT.txt"))
writeLines(pole_wt_mss_pre, file.path(out_dir, "preclinical_POLE_WT_MSS.txt"))

# phenotype table
pheno_pre <- data.frame(
  Sample = df1$Sample,
  Dataset = "Preclinical",
  MSI_status = df1$MSI_status,
  POLE_status = df1$POLE_status,
  stringsAsFactors = FALSE
)

# ---- Clinical ----
clin_cat <- load_catalogue("Clinical", cfg, root)
df2 <- as.data.frame(read_excel(xlsx, sheet = 2, skip = 1))
names(df2)[1:3] <- c("Sample_table", "MSI_status", "POLE_status")

map_clinical <- function(id, colnames_clin) {
  hits <- colnames_clin[colnames_clin == id | startsWith(colnames_clin, paste0(id, "-"))]
  if (length(hits) == 1L) return(hits)
  if (length(hits) > 1L) return(hits[1])
  NA_character_
}

df2$Sample <- vapply(df2$Sample_table, map_clinical, character(1), colnames_clin = colnames(clin_cat))
if (anyNA(df2$Sample)) {
  stop("Unmapped clinical samples: ", paste(df2$Sample_table[is.na(df2$Sample)], collapse = ", "))
}

# Paper treats MSI-L with MSS for some analyses; keep raw label, also binary MMRd
df2$MMRd <- ifelse(df2$MSI_status == "MSI", "MSI", "MSS")

msi_clin <- df2$Sample[df2$MMRd == "MSI"]
mss_clin <- df2$Sample[df2$MMRd == "MSS"]
pole_mut_clin <- df2$Sample[df2$POLE_status == "MUT" & df2$MMRd == "MSS"]
pole_wt_mss_clin <- setdiff(mss_clin, pole_mut_clin)

writeLines(msi_clin, file.path(out_dir, "clinical_MSI.txt"))
writeLines(mss_clin, file.path(out_dir, "clinical_MSS.txt"))
writeLines(pole_mut_clin, file.path(out_dir, "clinical_POLE_MUT.txt"))
writeLines(pole_wt_mss_clin, file.path(out_dir, "clinical_POLE_WT_MSS.txt"))

pheno_clin <- data.frame(
  Sample = df2$Sample,
  Dataset = "Clinical",
  MSI_status = df2$MSI_status,
  POLE_status = df2$POLE_status,
  stringsAsFactors = FALSE
)

pheno <- rbind(pheno_pre, pheno_clin)
write.table(pheno, file.path(out_dir, "sample_phenotypes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

summary <- list(
  preclinical = list(
    n = nrow(pheno_pre),
    MSI = length(msi_pre),
    MSS = length(mss_pre),
    POLE_MUT = length(pole_mut_pre)
  ),
  clinical = list(
    n = nrow(pheno_clin),
    MSI = length(msi_clin),
    MSS = length(mss_clin),
    POLE_MUT = length(pole_mut_clin)
  )
)
write_json(summary, file.path(out_dir, "phenotype_summary.json"), pretty = TRUE, auto_unbox = TRUE)

message("Wrote metadata to ", out_dir)
message(toJSON(summary, auto_unbox = TRUE, pretty = TRUE))

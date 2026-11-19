#!/usr/bin/env Rscript
## Build thesis figure pack + paper-vs-replicate comparison table.
## Run from mywork/: Rscript scripts/r/downstream/build_thesis_pack.R

suppressPackageStartupMessages(library(jsonlite))

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
source(file.path(dirname(script_dir), "helpers.R"))

root <- .mywork_root()
out_dir <- file.path(root, "results", "thesis")
fig_dir <- file.path(out_dir, "figures")
if (dir.exists(fig_dir)) unlink(fig_dir, recursive = TRUE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

delta <- read.table(file.path(root, "results", "downstream", "tables", "delta_MMR_POLE.tsv"),
                    header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
asn <- read.table(file.path(root, "results", "downstream", "tables", "assigned_fraction.tsv"),
                  header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
ngs <- read.table(file.path(root, "results", "downstream", "ngs_workflow", "tables", "ngs_delta_MMR_POLE.tsv"),
                  header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)

pick_delta <- function(tool, cosmic, dataset = "Preclinical", data_type = "WES") {
  row <- delta[delta$Tool == tool & delta$COSMIC_Version == cosmic &
                 delta$Dataset == dataset & delta$Data_Type == data_type, , drop = FALSE]
  if (!nrow(row)) return(NA_real_)
  row$delta_MMR[[1]]
}

pick_assigned <- function(tool, cosmic, dataset = "Preclinical") {
  row <- asn[asn$Tool == tool & asn$COSMIC_Version == cosmic & asn$Dataset == dataset, , drop = FALSE]
  if (!nrow(row)) return(c(n_assigned = NA_real_, n = NA_real_, frac = NA_real_))
  c(n_assigned = row$n_assigned[[1]], n = row$n[[1]], frac = row$frac_assigned[[1]])
}

pick_ngs <- function(tool, cosmic, data_type) {
  row <- ngs[ngs$Tool == tool & ngs$COSMIC_Version == cosmic & ngs$Data_Type == data_type, , drop = FALSE]
  if (!nrow(row)) return(c(delta_MMR = NA_real_, median_flat_MSI = NA_real_))
  c(delta_MMR = row$delta_MMR[[1]], median_flat_MSI = row$median_flat_MSI[[1]])
}

rows <- list(
  data.frame(
    Analysis = "Algorithm", Cohort = "Preclinical", Data_Type = "WES",
    Tool = "MutationalPatterns", COSMIC_Version = "v2", Metric = "delta_MMR",
    Paper = 0.34, Paper_note = "approx. from paper text (~0.34)",
    Replicate = pick_delta("MutationalPatterns", "v2"),
    stringsAsFactors = FALSE
  ),
  data.frame(
    Analysis = "Algorithm", Cohort = "Preclinical", Data_Type = "WES",
    Tool = "MutationalPatterns", COSMIC_Version = "v3.2", Metric = "delta_MMR",
    Paper = 0.26, Paper_note = "approx. from paper text (~0.26)",
    Replicate = pick_delta("MutationalPatterns", "v3.2"),
    stringsAsFactors = FALSE
  ),
  data.frame(
    Analysis = "Algorithm", Cohort = "Preclinical", Data_Type = "WES",
    Tool = "SigProfilerAssignment", COSMIC_Version = "v2", Metric = "delta_MMR",
    Paper = NA_real_, Paper_note = "highest ΔMMR among paper tool set",
    Replicate = pick_delta("SigProfilerAssignment", "v2"),
    stringsAsFactors = FALSE
  ),
  {
    spa <- pick_assigned("SigProfilerAssignment", "v2")
    data.frame(
      Analysis = "Algorithm", Cohort = "Preclinical", Data_Type = "WES",
      Tool = "SigProfilerAssignment", COSMIC_Version = "v2", Metric = "frac_assigned_ge_0.9",
      Paper = 0.79, Paper_note = "approx. ~79% assigned (cosine ≥ 0.9)",
      Replicate = unname(spa[["frac"]]),
      stringsAsFactors = FALSE
    )
  },
  {
    w <- pick_ngs("MutationalPatterns", "v3.2", "WES")
    data.frame(
      Analysis = "NGS_workflow", Cohort = "Preclinical", Data_Type = "WES",
      Tool = "MutationalPatterns", COSMIC_Version = "v3.2", Metric = "delta_MMR",
      Paper = NA_real_, Paper_note = "Fig. 3 pattern: intermediate vs TSO/WGS",
      Replicate = unname(w[["delta_MMR"]]),
      stringsAsFactors = FALSE
    )
  },
  {
    w <- pick_ngs("MutationalPatterns", "v3.2", "WGS")
    data.frame(
      Analysis = "NGS_workflow", Cohort = "Preclinical", Data_Type = "WGS",
      Tool = "MutationalPatterns", COSMIC_Version = "v3.2", Metric = "delta_MMR",
      Paper = NA_real_, Paper_note = "Fig. 3 pattern: lowest ΔMMR",
      Replicate = unname(w[["delta_MMR"]]),
      stringsAsFactors = FALSE
    )
  },
  {
    w <- pick_ngs("MutationalPatterns", "v3.2", "TSO-500")
    data.frame(
      Analysis = "NGS_workflow", Cohort = "Preclinical", Data_Type = "TSO-500",
      Tool = "MutationalPatterns", COSMIC_Version = "v3.2", Metric = "delta_MMR",
      Paper = NA_real_, Paper_note = "Fig. 3 pattern: highest ΔMMR",
      Replicate = unname(w[["delta_MMR"]]),
      stringsAsFactors = FALSE
    )
  },
  {
    w <- pick_ngs("MutationalPatterns", "v3.2", "WGS")
    data.frame(
      Analysis = "NGS_workflow", Cohort = "Preclinical", Data_Type = "WGS",
      Tool = "MutationalPatterns", COSMIC_Version = "v3.2", Metric = "median_flat_MSI",
      Paper = NA_real_, Paper_note = "Fig. 3D: flat signatures elevated on WGS MSI",
      Replicate = unname(w[["median_flat_MSI"]]),
      stringsAsFactors = FALSE
    )
  }
)

cmp <- do.call(rbind, rows)
cmp$Abs_diff <- ifelse(is.finite(cmp$Paper) & is.finite(cmp$Replicate),
                       abs(cmp$Replicate - cmp$Paper), NA_real_)
cmp$Match <- NA

num_idx <- is.finite(cmp$Paper) & is.finite(cmp$Replicate) &
  cmp$Metric %in% c("delta_MMR", "frac_assigned_ge_0.9")
cmp$Match[num_idx] <- abs(cmp$Replicate[num_idx] - cmp$Paper[num_idx]) < 0.02

paper_tools <- c("MutationalPatterns", "deconstructSigs", "signature.tools.lib",
                 "SigProfilerAssignment", "MuSiCal")
wes_v2 <- delta[delta$Dataset == "Preclinical" & delta$Data_Type == "WES" &
                  delta$COSMIC_Version == "v2" & delta$Tool %in% paper_tools, ]
spa_highest <- nrow(wes_v2) > 0 &&
  wes_v2$Tool[which.max(wes_v2$delta_MMR)] == "SigProfilerAssignment"
spa_row <- cmp$Tool == "SigProfilerAssignment" & cmp$Metric == "delta_MMR" &
  cmp$Analysis == "Algorithm"
cmp$Match[spa_row] <- spa_highest

ngs_ord <- ngs[ngs$Tool == "MutationalPatterns" & ngs$COSMIC_Version == "v3.2", ]
ngs_ord <- ngs_ord[match(c("TSO-500", "WES", "WGS"), ngs_ord$Data_Type), ]
order_ok <- isTRUE(all(diff(ngs_ord$delta_MMR) < 0))
cmp$Match[cmp$Analysis == "NGS_workflow" & cmp$Metric == "delta_MMR"] <- order_ok
cmp$Match[cmp$Metric == "median_flat_MSI"] <-
  isTRUE(cmp$Replicate[cmp$Metric == "median_flat_MSI"][[1]] > 0.1)

write.table(cmp, file.path(out_dir, "paper_vs_replicate.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

## Prefer paper_style PNG/SVG panels; fall back to overview grids.
paper_ds <- file.path(root, "results", "downstream", "plots", "paper_style")
paper_ngs <- file.path(root, "results", "downstream", "ngs_workflow", "plots", "paper_style")
legacy_ds <- file.path(root, "results", "downstream", "plots")
legacy_ngs <- file.path(root, "results", "downstream", "ngs_workflow", "plots")

pick_src <- function(...) {
  cands <- c(...)
  for (p in cands) if (file.exists(p)) return(p)
  NA_character_
}

map <- list(
  "Fig01_alg_cosine_by_tool_v2" = pick_src(
    file.path(paper_ds, "alg_cosine_by_tool_v2.png"),
    file.path(legacy_ds, "cosine_by_tool_cosmic.png")
  ),
  "Fig02_alg_MMR_boxplot_by_tool_v2" = pick_src(
    file.path(paper_ds, "alg_MMR_boxplot_by_tool_v2.png")
  ),
  "Fig03_alg_delta_MMR_violin_v2" = pick_src(
    file.path(paper_ds, "alg_delta_MMR_violin_v2.png"),
    file.path(legacy_ds, "delta_MMR_by_tool_cosmic.png")
  ),
  "Fig04_alg_POLE_boxplot_by_tool_v2" = pick_src(
    file.path(paper_ds, "alg_POLE_boxplot_by_tool_v2.png"),
    file.path(legacy_ds, "delta_POLE_by_tool_cosmic.png")
  ),
  "Fig05_alg_delta_POLE_violin_v2" = pick_src(
    file.path(paper_ds, "alg_delta_POLE_violin_v2.png")
  ),
  "Fig06_ref_MMR_composition_MP" = pick_src(
    file.path(paper_ds, "ref_MMR_composition_MP.png"),
    file.path(legacy_ds, "MMR_composition_MSI.png")
  ),
  "Fig07_ngs_cosine_MP_v3.2" = pick_src(
    file.path(paper_ngs, "ngs_cosine_MP_v3.2.png"),
    file.path(legacy_ngs, "ngs_cosine_MP_v3.2.png")
  ),
  "Fig08_ngs_MMR_boxplot_MP_v3.2" = pick_src(
    file.path(paper_ngs, "ngs_MMR_boxplot_MP_v3.2.png"),
    file.path(legacy_ngs, "ngs_delta_MMR_MP_v3.2.png")
  ),
  "Fig09_ngs_POLE_boxplot_MP_v3.2" = pick_src(
    file.path(paper_ngs, "ngs_POLE_boxplot_MP_v3.2.png")
  ),
  "Fig10_ngs_signature_classes_MP_v3.2" = pick_src(
    file.path(paper_ngs, "ngs_signature_classes_MP_v3.2.png"),
    file.path(legacy_ngs, "ngs_flat_MSI_MP_v3.2.png")
  ),
  "FigS1_cosine_by_tool_cosmic" = pick_src(file.path(legacy_ds, "cosine_by_tool_cosmic.png")),
  "FigS2_delta_MMR_by_tool_cosmic" = pick_src(file.path(legacy_ds, "delta_MMR_by_tool_cosmic.png")),
  "FigS3_ngs_cosine_by_datatype" = pick_src(file.path(legacy_ngs, "ngs_cosine_by_datatype.png")),
  "FigS4_ngs_delta_MMR" = pick_src(file.path(legacy_ngs, "ngs_delta_MMR.png")),
  "FigS5_ngs_delta_POLE" = pick_src(file.path(legacy_ngs, "ngs_delta_POLE.png")),
  "Fig11_ref_cosine_v2_v32_TS" = pick_src(
    file.path(root, "results", "downstream", "reference_expansion", "plots", "ref_cosine_v2_v32_TS.png")
  ),
  "Fig12_ref_MMR_boxplot_preclinical" = pick_src(
    file.path(root, "results", "downstream", "reference_expansion", "plots", "ref_MMR_boxplot_preclinical.png")
  ),
  "Fig13_ref_MMR_composition_v2_v32_TS" = pick_src(
    file.path(root, "results", "downstream", "reference_expansion", "plots", "ref_MMR_composition_v2_v32_TS.png")
  ),
  "Fig14_downsample_Preclinical" = pick_src(
    file.path(root, "results", "downstream", "mutation_threshold", "plots", "Fig6A_preclinical_cosine_vs_mutations.png"),
    file.path(root, "results", "downstream", "mutation_threshold", "plots", "downsample_Preclinical.png")
  ),
  "Fig15_downsample_Clinical" = pick_src(
    file.path(root, "results", "downstream", "mutation_threshold", "plots", "Fig6B_clinical_cosine_vs_mutations.png"),
    file.path(root, "results", "downstream", "mutation_threshold", "plots", "downsample_Clinical.png")
  ),
  "Fig6A_preclinical_cosine_vs_mutations" = pick_src(
    file.path(root, "results", "downstream", "mutation_threshold", "plots", "Fig6A_preclinical_cosine_vs_mutations.png")
  ),
  "Fig6B_clinical_cosine_vs_mutations" = pick_src(
    file.path(root, "results", "downstream", "mutation_threshold", "plots", "Fig6B_clinical_cosine_vs_mutations.png")
  ),
  "Fig6C_preclinical_threshold_by_cosmic" = pick_src(
    file.path(root, "results", "downstream", "mutation_threshold", "plots", "Fig6C_preclinical_threshold_by_cosmic.png")
  ),
  "Fig6D_clinical_threshold_by_cosmic" = pick_src(
    file.path(root, "results", "downstream", "mutation_threshold", "plots", "Fig6D_clinical_threshold_by_cosmic.png")
  ),
  "Fig6E_preclinical_assigned_vs_mutations" = pick_src(
    file.path(root, "results", "downstream", "mutation_threshold", "plots", "Fig6E_preclinical_assigned_vs_mutations.png")
  ),
  "Fig6F_clinical_assigned_vs_mutations" = pick_src(
    file.path(root, "results", "downstream", "mutation_threshold", "plots", "Fig6F_clinical_assigned_vs_mutations.png")
  )
)

copied <- character()
missing <- character()
for (nm in names(map)) {
  src <- map[[nm]]
  if (is.na(src) || !nzchar(src) || !file.exists(src)) {
    missing <- c(missing, nm)
    next
  }
  ext <- tools::file_ext(src)
  dst <- file.path(fig_dir, paste0(nm, ".", ext))
  file.copy(src, dst, overwrite = TRUE)
  # also copy sibling SVG when present
  svg_src <- sub("\\.png$", ".svg", src)
  if (file.exists(svg_src)) {
    file.copy(svg_src, file.path(fig_dir, paste0(nm, ".svg")), overwrite = TRUE)
  }
  copied <- c(copied, nm)
}

manifest <- paste0(
  "# Thesis figure pack (paper style)\n\n",
  "Generated by `scripts/r/downstream/build_thesis_pack.R`.\n",
  "Style follows Battuello et al. *Brief Bioinform* 2024 ([bbae249](https://doi.org/10.1093/bib/bbae249)): ",
  "`theme_bw`, MSI/MSS (`lightblue1`/`salmon2`), POLE (`darkgoldenrod1`/`olivedrab4`), ",
  "dashed cosine = 0.9, SVG+PNG export.\n\n",
  "## Main figures\n\n",
  "| File | Paper analogue |\n|------|----------------|\n",
  "| `Fig01_alg_cosine_by_tool_v2` | Fig. 4B — cosine by algorithm |\n",
  "| `Fig02_alg_MMR_boxplot_by_tool_v2` | Fig. 4C — MMR contribution MSI vs MSS |\n",
  "| `Fig03_alg_delta_MMR_violin_v2` | Fig. 4D — ΔMMR violin |\n",
  "| `Fig04_alg_POLE_boxplot_by_tool_v2` | Fig. 4E — POLE contribution |\n",
  "| `Fig05_alg_delta_POLE_violin_v2` | Fig. 4F — ΔPOLE violin |\n",
  "| `Fig06_ref_MMR_composition_MP` | Fig. 5D — MMR SBS mix by reference |\n",
  "| `Fig07_ngs_cosine_MP_v3.2` | Fig. 3A — cosine by NGS (horizontal grey) |\n",
  "| `Fig08_ngs_MMR_boxplot_MP_v3.2` | Fig. 3B — MMR by NGS |\n",
  "| `Fig09_ngs_POLE_boxplot_MP_v3.2` | Fig. 3C — POLE by NGS |\n",
  "| `Fig10_ngs_signature_classes_MP_v3.2` | Fig. 3D — signature classes in MSI |\n",
  "| `Fig11_ref_cosine_v2_v32_TS` | Fig. 5A — cosine by reference (incl. TS) |\n",
  "| `Fig12_ref_MMR_boxplot_preclinical` | Fig. 5B — MMR by reference |\n",
  "| `Fig13_ref_MMR_composition_v2_v32_TS` | Fig. 5D — MMR SBS mix by reference |\n",
  "| `Fig14_downsample_Preclinical` | Fig. 6A alias |\n",
  "| `Fig15_downsample_Clinical` | Fig. 6B alias |\n",
  "| `Fig6A_preclinical_cosine_vs_mutations` | Fig. 6A — preclinical cosine vs mutations (all COSMIC) |\n",
  "| `Fig6B_clinical_cosine_vs_mutations` | Fig. 6B — clinical cosine vs mutations |\n",
  "| `Fig6C_preclinical_threshold_by_cosmic` | Fig. 6C — preclinical threshold by COSMIC |\n",
  "| `Fig6D_clinical_threshold_by_cosmic` | Fig. 6D — clinical threshold by COSMIC |\n",
  "| `Fig6E_preclinical_assigned_vs_mutations` | Fig. 6E — preclinical assigned fraction |\n",
  "| `Fig6F_clinical_assigned_vs_mutations` | Fig. 6F — clinical assigned fraction |\n\n",
  "## Supplementary (full grids)\n\n",
  "| File | Role |\n|------|------|\n",
  "| `FigS1_cosine_by_tool_cosmic` | All tools × COSMIC cosine |\n",
  "| `FigS2_delta_MMR_by_tool_cosmic` | ΔMMR bar overview |\n",
  "| `FigS3_ngs_cosine_by_datatype` | Full NGS cosine grid |\n",
  "| `FigS4_ngs_delta_MMR` | Full NGS ΔMMR grid |\n",
  "| `FigS5_ngs_delta_POLE` | Full NGS ΔPOLE grid |\n\n",
  "## Comparison table\n\n",
  "`paper_vs_replicate.tsv`\n\n",
  sprintf("Copied: %d / %d figures.\n", length(copied), length(map)),
  if (length(missing)) paste0("Missing: ", paste(missing, collapse = ", "), "\n") else ""
)
writeLines(manifest, file.path(out_dir, "FIGURES.md"))

message("Wrote ", file.path(out_dir, "paper_vs_replicate.tsv"))
message("Figures: ", length(copied), "/", length(map), " -> ", fig_dir)
if (length(missing)) message("Missing: ", paste(missing, collapse = ", "))
message("Manifest: ", file.path(out_dir, "FIGURES.md"))

#!/usr/bin/env Rscript
## Run R fitting grid: datasets × COSMIC × enabled R tools
## Usage (from mywork/):
##   Rscript scripts/r/run_r_fitting.R
##   Rscript scripts/r/run_r_fitting.R --tools MutationalPatterns --cosmic v2,v3.2 --datasets Preclinical

suppressPackageStartupMessages({
  library(jsonlite)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
source(file.path(script_dir, "helpers.R"))
source(file.path(script_dir, "fitters_r.R"))

parse_cli <- function(argv) {
  out <- list(tools = NULL, cosmic = NULL, datasets = NULL, phase = NULL)
  i <- 1L
  while (i <= length(argv)) {
    if (argv[i] == "--tools" && i < length(argv)) {
      out$tools <- strsplit(argv[i + 1], ",", fixed = TRUE)[[1]]
      i <- i + 2L
    } else if (argv[i] == "--cosmic" && i < length(argv)) {
      out$cosmic <- strsplit(argv[i + 1], ",", fixed = TRUE)[[1]]
      i <- i + 2L
    } else if (argv[i] == "--datasets" && i < length(argv)) {
      out$datasets <- strsplit(argv[i + 1], ",", fixed = TRUE)[[1]]
      i <- i + 2L
    } else if (argv[i] == "--phase" && i < length(argv)) {
      out$phase <- as.integer(argv[i + 1])
      i <- i + 2L
    } else {
      i <- i + 1L
    }
  }
  out
}

cli <- parse_cli(commandArgs(trailingOnly = TRUE))
root <- .mywork_root()
cfg <- load_pipeline_config(root)
log_dir <- file.path(root, "results", "logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(log_dir, paste0("r_fitting_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")
on.exit({
  sink(type = "message")
  sink(type = "output")
  close(log_con)
}, add = TRUE)

r_tools <- names(Filter(function(t) isTRUE(t$enabled) && identical(t$language, "R"), cfg$tools))
if (!is.null(cli$tools)) r_tools <- intersect(r_tools, cli$tools)
cosmic_keys <- names(cfg$cosmic)
if (!is.null(cli$cosmic)) cosmic_keys <- intersect(cosmic_keys, cli$cosmic)
datasets <- names(cfg$datasets)
if (!is.null(cli$phase)) {
  datasets <- datasets[vapply(datasets, function(d) {
    ph <- cfg$datasets[[d]]$phase
    !is.null(ph) && as.integer(ph) == cli$phase
  }, logical(1))]
}
if (!is.null(cli$datasets)) datasets <- intersect(datasets, cli$datasets)

message("Tools: ", paste(r_tools, collapse = ", "))
message("COSMIC: ", paste(cosmic_keys, collapse = ", "))
message("Datasets: ", paste(datasets, collapse = ", "))
if (!length(datasets)) stop("No datasets selected (check --phase / --datasets).")

# Preload catalogues
catalogues <- lapply(datasets, function(d) load_catalogue(d, cfg, root))
names(catalogues) <- datasets

for (tool in r_tools) {
  if (!tool %in% names(FITTERS_R)) {
    message("SKIP unknown fitter: ", tool)
    next
  }
  fit_fun <- FITTERS_R[[tool]]
  message("\n==== Tool: ", tool, " ====")

  for (ck in cosmic_keys) {
    message("-- COSMIC ", ck)
    cosmic <- tryCatch(load_cosmic(ck, cfg, root), error = function(e) {
      message("  ERROR loading COSMIC: ", conditionMessage(e))
      NULL
    })
    if (is.null(cosmic)) next

    for (ds in datasets) {
      message("  Dataset: ", ds)
      data_type <- cfg$datasets[[ds]]$data_type
      cat_mat <- catalogues[[ds]]
      fit_raw <- tryCatch(fit_fun(cat_mat, cosmic), error = function(e) {
        message("    FIT ERROR: ", conditionMessage(e))
        NULL
      })
      if (is.null(fit_raw)) next

      fit_u <- tryCatch(ensure_uniform_fit(fit_raw, cat_mat, cosmic), error = function(e) {
        message("    UNIFORMIZE ERROR: ", conditionMessage(e))
        NULL
      })
      if (is.null(fit_u)) next

      cos_df <- tryCatch(
        save_fit_outputs(fit_u, tool, ck, ds, data_type, root),
        error = function(e) {
          message("    SAVE ERROR: ", conditionMessage(e))
          NULL
        }
      )
      if (is.null(cos_df)) next
      append_cosine_tsv(cos_df, tool, root)
      message("    OK: median cosine=", signif(median(cos_df$Cosine_Similarity, na.rm = TRUE), 4),
              " assigned(>=0.9)=", sum(cos_df$Cosine_Similarity >= 0.9, na.rm = TRUE), "/", nrow(cos_df))
    }
  }
}

message("\nDone. Log: ", log_file)

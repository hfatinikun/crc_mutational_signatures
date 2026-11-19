#!/usr/bin/env Rscript
## Shared loaders / path helpers for the CRC signature pipeline.
## Assumes getwd() == mywork/

suppressPackageStartupMessages({
  library(jsonlite)
})

.mywork_root <- function() {
  # Prefer cwd if it looks like mywork; else walk up looking for config/pipeline.json
  if (file.exists(file.path(getwd(), "config", "pipeline.json"))) {
    return(normalizePath(getwd()))
  }
  cur <- normalizePath(getwd())
  for (i in seq_len(6)) {
    cand <- file.path(cur, "mywork")
    if (file.exists(file.path(cand, "config", "pipeline.json"))) {
      return(normalizePath(cand))
    }
    parent <- dirname(cur)
    if (identical(parent, cur)) break
    cur <- parent
  }
  stop("Cannot locate mywork/ (config/pipeline.json). Run from mywork/.")
}

load_pipeline_config <- function(root = .mywork_root()) {
  cfg_path <- file.path(root, "config", "pipeline.json")
  fromJSON(cfg_path, simplifyVector = TRUE)
}

load_signature_groups <- function(root = .mywork_root()) {
  fromJSON(file.path(root, "config", "signature_groups.json"), simplifyVector = TRUE)
}

load_sample_aliases <- function(root = .mywork_root()) {
  fromJSON(file.path(root, "config", "sample_aliases.json"), simplifyVector = TRUE)
}

cosmic_version_label <- function(cosmic_key) {
  # keys are already "v2", "v3.2", ...
  as.character(cosmic_key)
}

load_catalogue <- function(dataset_name, cfg = load_pipeline_config(), root = .mywork_root()) {
  stopifnot(dataset_name %in% names(cfg$datasets))
  rel <- cfg$datasets[[dataset_name]]$file
  path <- file.path(root, cfg$paths$input_dir, rel)
  mat <- as.matrix(read.table(path, header = TRUE, sep = "\t", row.names = 1,
                              check.names = FALSE, quote = ""))
  storage.mode(mat) <- "numeric"
  mat
}

load_cosmic <- function(cosmic_key, cfg = load_pipeline_config(), root = .mywork_root()) {
  stopifnot(cosmic_key %in% names(cfg$cosmic))
  path <- file.path(root, cfg$paths$cosmic_dir, cfg$cosmic[[cosmic_key]])
  mat <- as.matrix(read.table(path, header = TRUE, sep = "\t", row.names = 1,
                              check.names = FALSE, quote = ""))
  storage.mode(mat) <- "numeric"
  mat
}

align_contexts <- function(catalogue, cosmic) {
  catalogue <- as.matrix(catalogue)
  cosmic <- as.matrix(cosmic)
  common <- intersect(rownames(catalogue), rownames(cosmic))
  if (length(common) == 0L) stop("No overlapping mutation contexts.")
  # Prefer COSMIC context order for stability
  common <- rownames(cosmic)[rownames(cosmic) %in% common]
  list(
    catalogue = catalogue[common, , drop = FALSE],
    cosmic = cosmic[common, , drop = FALSE]
  )
}

normalize_contributions <- function(contribution) {
  contribution <- as.matrix(contribution)
  cs <- colSums(contribution)
  cs[cs == 0] <- 1
  sweep(contribution, 2, cs, "/")
}

cosine_per_sample <- function(observed, reconstructed) {
  observed <- as.matrix(observed)
  reconstructed <- as.matrix(reconstructed)
  common_ctx <- intersect(rownames(observed), rownames(reconstructed))
  common_smp <- intersect(colnames(observed), colnames(reconstructed))
  if (length(common_ctx) == 0L || length(common_smp) == 0L) {
    stop("Cannot compute cosine: no overlapping contexts/samples.")
  }
  X <- observed[common_ctx, common_smp, drop = FALSE]
  Y <- reconstructed[common_ctx, common_smp, drop = FALSE]
  # Prefer MutationalPatterns if available
  if (requireNamespace("MutationalPatterns", quietly = TRUE)) {
    return(diag(MutationalPatterns::cos_sim_matrix(X, Y)))
  }
  # Fallback: column-wise cosine
  vapply(seq_along(common_smp), function(j) {
    x <- X[, j]; y <- Y[, j]
    nx <- sqrt(sum(x * x)); ny <- sqrt(sum(y * y))
    if (nx == 0 || ny == 0) return(NA_real_)
    sum(x * y) / (nx * ny)
  }, numeric(1))
}

mutation_counts <- function(catalogue) {
  colSums(as.matrix(catalogue))
}

write_matrix_tsv <- function(mat, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(mat, file = path, sep = "\t", quote = FALSE, col.names = NA)
}

ensure_uniform_fit <- function(fit, catalogue, cosmic) {
  # Enforce contribution = sigs x samples, reconstructed = contexts x samples
  stopifnot(!is.null(fit$contribution), !is.null(fit$reconstructed))
  contrib <- as.matrix(fit$contribution)
  recon <- as.matrix(fit$reconstructed)

  # Fix orientation if samples x sigs
  if (nrow(contrib) == ncol(catalogue) && ncol(contrib) == ncol(cosmic)) {
    contrib <- t(contrib)
  }
  if (is.null(rownames(contrib))) rownames(contrib) <- colnames(cosmic)
  if (is.null(colnames(contrib))) colnames(contrib) <- colnames(catalogue)

  common_sigs <- intersect(colnames(cosmic), rownames(contrib))
  if (length(common_sigs) == 0L) {
    # nameless fallback: assume order matches
    if (nrow(contrib) == ncol(cosmic)) {
      rownames(contrib) <- colnames(cosmic)
      common_sigs <- colnames(cosmic)
    } else {
      stop("Contribution matrix has no overlapping signature names with COSMIC.")
    }
  }
  cosmic_use <- cosmic[, common_sigs, drop = FALSE]
  contrib <- contrib[common_sigs, , drop = FALSE]

  # Align sample columns to catalogue
  common_smp <- intersect(colnames(catalogue), colnames(contrib))
  if (length(common_smp) == 0L) {
    if (ncol(contrib) == ncol(catalogue)) {
      colnames(contrib) <- colnames(catalogue)
      common_smp <- colnames(catalogue)
    } else {
      stop("Contribution samples do not match catalogue.")
    }
  }
  contrib <- contrib[, common_smp, drop = FALSE]
  catalogue_use <- catalogue[, common_smp, drop = FALSE]

  # Rebuild reconstruction for consistency
  recon <- cosmic_use %*% contrib
  rownames(recon) <- rownames(cosmic_use)
  colnames(recon) <- colnames(contrib)

  list(
    contribution = contrib,
    reconstructed = recon,
    catalogue_aligned = catalogue_use,
    cosmic_aligned = cosmic_use
  )
}

save_fit_outputs <- function(fit_uniform, tool, cosmic_key, dataset_name, data_type,
                             root = .mywork_root()) {
  ver <- cosmic_version_label(cosmic_key)
  contrib_dir <- file.path(root, "results", "contributions", tool, ver)
  recon_dir <- file.path(root, "results", "reconstructed", tool, ver)
  dir.create(contrib_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(recon_dir, recursive = TRUE, showWarnings = FALSE)

  prefix <- dataset_name
  write_matrix_tsv(fit_uniform$contribution,
                   file.path(contrib_dir, paste0(prefix, "_contribution.tsv")))
  write_matrix_tsv(normalize_contributions(fit_uniform$contribution),
                   file.path(contrib_dir, paste0(prefix, "_contribution_normalized.tsv")))
  write_matrix_tsv(fit_uniform$reconstructed,
                   file.path(recon_dir, paste0(prefix, "_reconstructed.tsv")))

  cos_vals <- cosine_per_sample(fit_uniform$catalogue_aligned, fit_uniform$reconstructed)
  mut_n <- mutation_counts(fit_uniform$catalogue_aligned)
  data.frame(
    Sample = names(cos_vals),
    Cosine_Similarity = as.numeric(cos_vals),
    Dataset = dataset_name,
    Data_Type = data_type,
    Tool = tool,
    COSMIC_Version = ver,
    Mutation_Count = as.numeric(mut_n[names(cos_vals)]),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

append_cosine_tsv <- function(df, tool, root = .mywork_root()) {
  out_dir <- file.path(root, "results", "cosine")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(out_dir, paste0("CosineResults_", tool, "_ALL.tsv"))
  if (file.exists(path)) {
    old <- read.table(path, header = TRUE, sep = "\t", check.names = FALSE,
                      stringsAsFactors = FALSE)
    # Include Data_Type so WES/WGS/TSO blocks do not overwrite each other
    key_new <- paste(df$Dataset, df$Data_Type, df$COSMIC_Version, sep = "||")
    key_old <- paste(old$Dataset, old$Data_Type, old$COSMIC_Version, sep = "||")
    old <- old[!key_old %in% unique(key_new), , drop = FALSE]
    df <- rbind(old, df)
  }
  write.table(df, file = path, sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(path)
}

dataset_cohort <- function(dataset_name, cfg = load_pipeline_config()) {
  if (!dataset_name %in% names(cfg$datasets)) return(dataset_name)
  cohort <- cfg$datasets[[dataset_name]]$cohort
  if (is.null(cohort) || !nzchar(cohort)) dataset_name else as.character(cohort)
}

dataset_data_type <- function(dataset_name, cfg = load_pipeline_config()) {
  if (!dataset_name %in% names(cfg$datasets)) return(NA_character_)
  as.character(cfg$datasets[[dataset_name]]$data_type)
}

resolve_signature_group <- function(group_name, cosmic_key, present_sigs,
                                    groups = load_signature_groups()) {
  key <- if (cosmic_key %in% names(groups)) cosmic_key else {
    # fallback to nearest known
    if (startsWith(cosmic_key, "v3")) "v3.2" else "v2"
  }
  wanted <- groups[[key]][[group_name]]
  if (is.null(wanted) || length(wanted) == 0L) return(character(0))
  wanted <- as.character(unlist(wanted, use.names = FALSE))
  intersect(wanted, present_sigs)
}

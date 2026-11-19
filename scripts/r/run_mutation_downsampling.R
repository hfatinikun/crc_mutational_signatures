#!/usr/bin/env Rscript
## Mutation downsampling across ALL configured COSMIC references (user analysis).
## Produces Fig. 6A–F under results/downstream/mutation_threshold/plots/.
##
## Panels:
##   A  Preclinical: median cosine vs median mutations (lines = COSMIC versions)
##   B  Clinical: same
##   C  Preclinical: inferred mutation threshold by COSMIC
##   D  Clinical: inferred mutation threshold by COSMIC
##   E  Preclinical: fraction assigned (cosine≥0.9) vs median mutations
##   F  Clinical: same
##
## Paper Fig. 6C–D (metanormal) remain blocked without PBMC data — see BLOCKED_AND_NOTES.md.
##
## Usage (from mywork/):
##   Rscript scripts/r/run_mutation_downsampling.R
##   Rscript scripts/r/run_mutation_downsampling.R --cosmic v2,v3.2,TS --reps 3

suppressPackageStartupMessages({
  library(jsonlite)
  library(ggplot2)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
source(file.path(script_dir, "helpers.R"))
source(file.path(script_dir, "fitters_r.R"))
source(file.path(script_dir, "downstream", "plot_style.R"))

parse_cli <- function(argv) {
  out <- list(datasets = c("Preclinical", "Clinical"), cosmic = NULL, reps = 5L, seed = 1L)
  i <- 1L
  while (i <= length(argv)) {
    if (argv[i] == "--datasets" && i < length(argv)) {
      out$datasets <- strsplit(argv[i + 1], ",", fixed = TRUE)[[1]]; i <- i + 2L
    } else if (argv[i] == "--cosmic" && i < length(argv)) {
      out$cosmic <- strsplit(argv[i + 1], ",", fixed = TRUE)[[1]]; i <- i + 2L
    } else if (argv[i] == "--reps" && i < length(argv)) {
      out$reps <- as.integer(argv[i + 1]); i <- i + 2L
    } else if (argv[i] == "--seed" && i < length(argv)) {
      out$seed <- as.integer(argv[i + 1]); i <- i + 2L
    } else i <- i + 1L
  }
  out
}

cli <- parse_cli(commandArgs(trailingOnly = TRUE))
root <- .mywork_root()
cfg <- load_pipeline_config(root)
set.seed(cli$seed)

cosmic_keys <- names(cfg$cosmic)
if (!is.null(cli$cosmic)) cosmic_keys <- intersect(cosmic_keys, cli$cosmic)
if (!length(cosmic_keys)) stop("No COSMIC versions selected.")

out_dir <- file.path(root, "results", "downstream", "mutation_threshold")
plot_dir <- file.path(out_dir, "plots")
tab_dir <- file.path(out_dir, "tables")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

fracs <- seq(0.05, 0.95, by = 0.05)
thr <- cfg$cosine_threshold

subsample_catalogue <- function(mat, frac) {
  mat <- as.matrix(mat)
  out <- mat
  for (j in seq_len(ncol(mat))) {
    n <- sum(mat[, j])
    keep <- max(1L, as.integer(round(n * frac)))
    if (n <= 0 || keep >= n) next
    probs <- mat[, j] / n
    out[, j] <- as.integer(stats::rmultinom(1, size = keep, prob = probs))
  }
  out
}

infer_threshold <- function(med_mut, med_cos, threshold = 0.9) {
  ord <- order(med_mut)
  med_mut <- med_mut[ord]; med_cos <- med_cos[ord]
  hit <- which(med_cos >= threshold)
  if (!length(hit)) return(NA_real_)
  i <- hit[[1]]
  if (i == 1L) return(med_mut[[i]])
  x0 <- med_mut[[i - 1]]; y0 <- med_cos[[i - 1]]
  x1 <- med_mut[[i]]; y1 <- med_cos[[i]]
  if (isTRUE(all.equal(y1, y0))) x1 else x0 + (threshold - y0) * (x1 - x0) / (y1 - y0)
}

# Preload catalogues + cosmics
catalogues <- lapply(cli$datasets, function(d) load_catalogue(d, cfg, root))
names(catalogues) <- cli$datasets
cosmics <- lapply(cosmic_keys, function(ck) load_cosmic(ck, cfg, root))
names(cosmics) <- cosmic_keys

message("Datasets: ", paste(cli$datasets, collapse = ", "))
message("COSMIC: ", paste(cosmic_keys, collapse = ", "))
message("Replicates: ", cli$reps)

all_rows <- list()

for (ds in cli$datasets) {
  cat0 <- catalogues[[ds]]
  message("==== ", ds, " (n=", ncol(cat0), ", median mut=", median(colSums(cat0)), ") ====")
  for (ck in cosmic_keys) {
    cosmic <- cosmics[[ck]]
    message("  -- ", ck)
    for (frac in fracs) {
      for (rep in seq_len(cli$reps)) {
        sub <- subsample_catalogue(cat0, frac)
        fit <- tryCatch(FITTERS_R$MutationalPatterns(sub, cosmic), error = function(e) NULL)
        if (is.null(fit)) next
        cos <- cosine_per_sample(sub, fit$reconstructed)
        mut_sub <- colSums(sub)
        for (i in seq_along(cos)) {
          all_rows[[length(all_rows) + 1L]] <- data.frame(
            Dataset = ds, COSMIC_Version = ck, Fraction = frac, Replicate = rep,
            Sample = names(cos)[i],
            Mutation_Count = unname(mut_sub[i]),
            Cosine_Similarity = unname(cos[i]),
            Assigned = is.finite(cos[i]) && cos[i] >= thr,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
}

long_df <- do.call(rbind, all_rows)
write.table(long_df, file.path(tab_dir, "downsampling_cosine_long.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

agg <- aggregate(
  cbind(Cosine_Similarity, Mutation_Count, Assigned) ~ Dataset + COSMIC_Version + Fraction,
  data = long_df,
  FUN = median, na.rm = TRUE
)
names(agg)[names(agg) == "Cosine_Similarity"] <- "median_cosine"
names(agg)[names(agg) == "Mutation_Count"] <- "median_mutations"
names(agg)[names(agg) == "Assigned"] <- "frac_assigned"
agg <- agg[order(agg$Dataset, agg$COSMIC_Version, agg$Fraction), ]
write.table(agg, file.path(tab_dir, "downsampling_median_by_fraction.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Thresholds per Dataset × COSMIC
thresh_rows <- list()
paper_targets <- c(Preclinical = 323, Clinical = 64)
for (ds in unique(agg$Dataset)) {
  for (ck in unique(agg$COSMIC_Version)) {
    sub <- agg[agg$Dataset == ds & agg$COSMIC_Version == ck, ]
    tmut <- infer_threshold(sub$median_mutations, sub$median_cosine, thr)
    thresh_rows[[length(thresh_rows) + 1L]] <- data.frame(
      Dataset = ds,
      COSMIC_Version = ck,
      threshold_mutations = tmut,
      paper_reported = unname(paper_targets[ds]),
      abs_diff_vs_paper = if (is.finite(tmut) && is.finite(paper_targets[ds])) abs(tmut - paper_targets[ds]) else NA_real_,
      stringsAsFactors = FALSE
    )
  }
}
thresh_df <- do.call(rbind, thresh_rows)
write.table(thresh_df, file.path(tab_dir, "mutation_count_threshold.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Colour map for COSMIC versions
cosmic_levels <- cosmic_keys
agg$COSMIC_Version <- factor(agg$COSMIC_Version, levels = cosmic_levels)
thresh_df$COSMIC_Version <- factor(thresh_df$COSMIC_Version, levels = cosmic_levels)
n_cv <- length(cosmic_levels)
cv_cols <- setNames(grDevices::hcl.colors(max(n_cv, 3L), "Dark 2")[seq_len(n_cv)], cosmic_levels)

plot_curve <- function(ds, y_col, ylab, fname, title) {
  sub <- agg[agg$Dataset == ds, ]
  p <- ggplot(sub, aes(x = median_mutations, y = .data[[y_col]], colour = COSMIC_Version)) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.6) +
    { if (identical(y_col, "median_cosine")) geom_hline(yintercept = thr, linetype = "dashed") else NULL } +
    { if (identical(y_col, "frac_assigned")) geom_hline(yintercept = 0.5, linetype = "dotted", colour = "grey40") else NULL } +
    scale_colour_manual(values = cv_cols) +
    paper_theme(10) +
    labs(title = title, x = "Median mutation count (subsample)", y = ylab, colour = "COSMIC")
  save_paper_plot(p, file.path(plot_dir, fname), width_cm = 12, height_cm = 8)
}

plot_threshold_bars <- function(ds, fname, title) {
  sub <- thresh_df[thresh_df$Dataset == ds, ]
  p <- ggplot(sub, aes(x = COSMIC_Version, y = threshold_mutations, fill = COSMIC_Version)) +
    geom_col(width = 0.7) +
    { if (is.finite(paper_targets[ds])) {
      geom_hline(yintercept = paper_targets[ds], linetype = "dashed", colour = "firebrick3")
    } else NULL } +
    scale_fill_manual(values = cv_cols) +
    paper_theme(10) +
    theme(legend.position = "none") +
    labs(title = title, x = "COSMIC version", y = "Inferred mutation threshold (cosine ≥ 0.9)")
  save_paper_plot(p, file.path(plot_dir, fname), width_cm = 12, height_cm = 8)
}

# Fig. 6A–F
plot_curve("Preclinical", "median_cosine", "Median cosine similarity",
           "Fig6A_preclinical_cosine_vs_mutations",
           "Fig. 6A — Preclinical: cosine vs mutation count (all COSMIC)")
plot_curve("Clinical", "median_cosine", "Median cosine similarity",
           "Fig6B_clinical_cosine_vs_mutations",
           "Fig. 6B — Clinical: cosine vs mutation count (all COSMIC)")
plot_threshold_bars("Preclinical", "Fig6C_preclinical_threshold_by_cosmic",
                    "Fig. 6C — Preclinical: mutation threshold by COSMIC")
plot_threshold_bars("Clinical", "Fig6D_clinical_threshold_by_cosmic",
                    "Fig. 6D — Clinical: mutation threshold by COSMIC")
plot_curve("Preclinical", "frac_assigned", "Fraction assigned (cosine ≥ 0.9)",
           "Fig6E_preclinical_assigned_vs_mutations",
           "Fig. 6E — Preclinical: assigned fraction vs mutation count")
plot_curve("Clinical", "frac_assigned", "Fraction assigned (cosine ≥ 0.9)",
           "Fig6F_clinical_assigned_vs_mutations",
           "Fig. 6F — Clinical: assigned fraction vs mutation count")

# Also keep legacy names pointing at A/B for older scripts
file.copy(file.path(plot_dir, "Fig6A_preclinical_cosine_vs_mutations.png"),
          file.path(plot_dir, "downsample_Preclinical.png"), overwrite = TRUE)
file.copy(file.path(plot_dir, "Fig6B_clinical_cosine_vs_mutations.png"),
          file.path(plot_dir, "downsample_Clinical.png"), overwrite = TRUE)

note <- c(
  "# Figure 6 expansion (multi-COSMIC)",
  "",
  "This analysis fits **MutationalPatterns** on mutation-subsampled catalogues for **every COSMIC key** in `config/pipeline.json`.",
  "",
  "## Panels produced",
  "- **Fig6A** Preclinical cosine vs mutations (lines = COSMIC versions)",
  "- **Fig6B** Clinical cosine vs mutations",
  "- **Fig6C** Preclinical inferred threshold by COSMIC (dashed = paper 323)",
  "- **Fig6D** Clinical inferred threshold by COSMIC (dashed = paper 64)",
  "- **Fig6E** Preclinical assigned fraction vs mutations",
  "- **Fig6F** Clinical assigned fraction vs mutations",
  "",
  "## Paper metanormal panels",
  "Paper Fig. 6C–D used a 21-PBMC metanormal. Those exact panels remain **blocked** (no PBMC/metanormal data locally).",
  "Our Fig6C–F instead report the multi-COSMIC expansion of the downsampling experiment.",
  "",
  paste0("COSMIC versions: ", paste(cosmic_keys, collapse = ", ")),
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
)
writeLines(note, file.path(out_dir, "BLOCKED_AND_NOTES.md"))

message("Wrote Fig6A–F under ", plot_dir)
print(thresh_df)

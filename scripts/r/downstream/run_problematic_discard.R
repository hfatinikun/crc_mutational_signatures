#!/usr/bin/env Rscript
## Problematic-sample QC + discard experiment (Maura-style caveats × Battuello readouts).
##
## Flags samples that are problematic under one or more rules, then recomputes
## ΔMMR / ΔPOLE with vs without those samples for every Tool × COSMIC × Dataset.
##
## Rules (configurable):
##   1. cosine < 0.9                         (poor reconstruction / unassigned)
##   2. Mutation_Count < cohort×COSMIC threshold from Fig.6 downsampling
##   3. flat-signature contribution ≥ flat_cut (default 0.25; Maura flat ambiguity)
##
## Optional: if metadata/clinical_survival.tsv exists (or is fetched), also compare
## Kaplan–Meier / Cox association of MMR contribution with OS before vs after discard.
##
## Usage (from mywork/):
##   Rscript scripts/r/downstream/run_problematic_discard.R
##   Rscript scripts/r/downstream/run_problematic_discard.R --flat-cut 0.25

suppressPackageStartupMessages({
  library(jsonlite)
  library(ggplot2)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
source(file.path(dirname(script_dir), "helpers.R"))
source(file.path(script_dir, "plot_style.R"))

parse_cli <- function(argv) {
  out <- list(flat_cut = 0.25, cosine_thr = NULL)
  i <- 1L
  while (i <= length(argv)) {
    if (argv[i] == "--flat-cut" && i < length(argv)) {
      out$flat_cut <- as.numeric(argv[i + 1]); i <- i + 2L
    } else if (argv[i] == "--cosine-thr" && i < length(argv)) {
      out$cosine_thr <- as.numeric(argv[i + 1]); i <- i + 2L
    } else i <- i + 1L
  }
  out
}

cli <- parse_cli(commandArgs(trailingOnly = TRUE))
root <- .mywork_root()
cfg <- load_pipeline_config(root)
groups <- load_signature_groups(root)
cos_thr <- if (is.null(cli$cosine_thr)) cfg$cosine_threshold else cli$cosine_thr
flat_cut <- cli$flat_cut

out_dir <- file.path(root, "results", "downstream", "problematic_discard")
dir.create(file.path(out_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

# ---------- Inputs ----------
cos_path <- file.path(root, "results", "cosine", "CosineResults_AllTools_AllCOSMIC.tsv")
if (!file.exists(cos_path)) stop("Missing ", cos_path)
cos_df <- read.table(cos_path, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)

pheno <- read.table(file.path(root, "metadata", "sample_phenotypes.tsv"),
                    header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

thr_path <- file.path(root, "results", "downstream", "mutation_threshold",
                      "tables", "mutation_count_threshold.tsv")
thr_df <- if (file.exists(thr_path)) {
  read.table(thr_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
} else NULL

sum_group <- function(norm_mat, sigs) {
  sigs <- intersect(sigs, rownames(norm_mat))
  if (!length(sigs)) return(setNames(rep(0, ncol(norm_mat)), colnames(norm_mat)))
  if (length(sigs) == 1) setNames(as.numeric(norm_mat[sigs, ]), colnames(norm_mat)) else
    colSums(norm_mat[sigs, , drop = FALSE])
}

# Restrict to WES algorithm/reference axis (Preclinical + Clinical) for biological readout
focus_ds <- c("Preclinical", "Clinical")
cos_focus <- cos_df[cos_df$Dataset %in% focus_ds &
                      (is.na(cos_df$Data_Type) | cos_df$Data_Type == "WES"), , drop = FALSE]

contrib_root <- file.path(root, "results", "contributions")
flag_rows <- list()
delta_rows <- list()

tools <- sort(unique(cos_focus$Tool))
cosmics <- sort(unique(cos_focus$COSMIC_Version))

for (tool in tools) {
  for (ck in cosmics) {
    for (ds in focus_ds) {
      nf <- file.path(contrib_root, tool, ck, paste0(ds, "_contribution_normalized.tsv"))
      if (!file.exists(nf)) next

      mat <- as.matrix(read.table(nf, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE))
      cohort <- dataset_cohort(ds, cfg)
      ph <- pheno[pheno$Dataset == cohort & pheno$Sample %in% colnames(mat), , drop = FALSE]
      if (!nrow(ph)) next

      cos_sub <- cos_focus[cos_focus$Tool == tool & cos_focus$COSMIC_Version == ck &
                             cos_focus$Dataset == ds, , drop = FALSE]
      cos_map <- setNames(cos_sub$Cosine_Similarity, cos_sub$Sample)
      mut_map <- setNames(cos_sub$Mutation_Count, cos_sub$Sample)

      mmr_sigs <- as.character(resolve_signature_group("MMR", ck, rownames(mat), groups))
      pole_sigs <- as.character(resolve_signature_group("POLE", ck, rownames(mat), groups))
      flat_sigs <- as.character(resolve_signature_group("flat", ck, rownames(mat), groups))

      mmr_score <- sum_group(mat, mmr_sigs)
      pole_score <- sum_group(mat, pole_sigs)
      flat_score <- sum_group(mat, flat_sigs)

      # Mutation threshold for this cohort×COSMIC (fallback: paper-ish defaults)
      mut_thr <- NA_real_
      if (!is.null(thr_df)) {
        hit <- thr_df$threshold_mutations[thr_df$Dataset == ds & thr_df$COSMIC_Version == ck]
        if (length(hit)) mut_thr <- hit[[1]]
      }
      if (!is.finite(mut_thr)) {
        mut_thr <- if (identical(ds, "Preclinical")) 323 else 64
      }

      samples <- ph$Sample
      for (s in samples) {
        cos_v <- unname(cos_map[s])
        mut_v <- unname(mut_map[s])
        flat_v <- unname(flat_score[s])
        if (!is.finite(cos_v)) cos_v <- NA_real_
        if (!is.finite(mut_v)) mut_v <- NA_real_
        if (!is.finite(flat_v)) flat_v <- 0

        flag_cos <- is.finite(cos_v) && cos_v < cos_thr
        flag_mut <- is.finite(mut_v) && mut_v < mut_thr
        flag_flat <- flat_v >= flat_cut
        problematic <- flag_cos || flag_mut || flag_flat

        ms <- ph$MSI_status[ph$Sample == s][[1]]
        ps <- ph$POLE_status[ph$Sample == s][[1]]
        flag_rows[[length(flag_rows) + 1L]] <- data.frame(
          Tool = tool, COSMIC_Version = ck, Dataset = ds, Sample = s,
          MSI_status = ms, POLE_status = ps,
          Cosine_Similarity = cos_v, Mutation_Count = mut_v,
          Flat_Contribution = flat_v,
          mut_threshold = mut_thr, flat_cut = flat_cut, cosine_threshold = cos_thr,
          flag_low_cosine = flag_cos,
          flag_low_mutations = flag_mut,
          flag_high_flat = flag_flat,
          problematic = problematic,
          MMR_Contribution = unname(mmr_score[s]),
          POLE_Contribution = unname(pole_score[s]),
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

flags <- do.call(rbind, flag_rows)
write.table(flags, file.path(out_dir, "tables", "sample_problematic_flags.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Recompute deltas from flags table (cleaner)
delta_rows <- list()
for (tool in unique(flags$Tool)) {
  for (ck in unique(flags$COSMIC_Version)) {
    for (ds in unique(flags$Dataset)) {
      sub <- flags[flags$Tool == tool & flags$COSMIC_Version == ck & flags$Dataset == ds, ]
      if (!nrow(sub)) next

      mmr <- setNames(sub$MMR_Contribution, sub$Sample)
      pole <- setNames(sub$POLE_Contribution, sub$Sample)

      make_row <- function(keep, label) {
        sk <- sub[sub$Sample %in% keep, ]
        msi <- sk$Sample[sk$MSI_status == "MSI"]
        mss <- sk$Sample[sk$MSI_status %in% c("MSS", "MSI-L")]
        if (!length(mss)) mss <- sk$Sample[sk$MSI_status != "MSI"]
        pole_mut <- sk$Sample[sk$POLE_status == "MUT" & sk$Sample %in% mss]
        pole_wt <- setdiff(mss, pole_mut)
        data.frame(
          Tool = tool, COSMIC_Version = ck, Dataset = ds, Filter = label,
          n_total = nrow(sub), n_kept = length(keep),
          n_discarded = nrow(sub) - length(keep),
          frac_discarded = (nrow(sub) - length(keep)) / nrow(sub),
          n_MSI = length(msi), n_MSS = length(mss), n_POLE_MUT = length(pole_mut),
          delta_MMR = if (length(msi) && length(mss))
            median(mmr[msi], na.rm = TRUE) - median(mmr[mss], na.rm = TRUE) else NA_real_,
          p_MMR = if (length(msi) >= 2 && length(mss) >= 2)
            suppressWarnings(wilcox.test(mmr[msi], mmr[mss])$p.value) else NA_real_,
          delta_POLE = if (length(pole_mut) && length(pole_wt))
            median(pole[pole_mut], na.rm = TRUE) - median(pole[pole_wt], na.rm = TRUE) else NA_real_,
          p_POLE = if (length(pole_mut) >= 2 && length(pole_wt) >= 2)
            suppressWarnings(wilcox.test(pole[pole_mut], pole[pole_wt])$p.value) else NA_real_,
          stringsAsFactors = FALSE
        )
      }

      all_s <- sub$Sample
      keep_any <- sub$Sample[!sub$problematic]
      keep_cos <- sub$Sample[!sub$flag_low_cosine]
      keep_mut <- sub$Sample[!sub$flag_low_mutations]
      keep_flat <- sub$Sample[!sub$flag_high_flat]

      delta_rows[[length(delta_rows) + 1L]] <- make_row(all_s, "all_samples")
      delta_rows[[length(delta_rows) + 1L]] <- make_row(keep_cos, "drop_low_cosine")
      delta_rows[[length(delta_rows) + 1L]] <- make_row(keep_mut, "drop_low_mutations")
      delta_rows[[length(delta_rows) + 1L]] <- make_row(keep_flat, "drop_high_flat")
      delta_rows[[length(delta_rows) + 1L]] <- make_row(keep_any, "drop_any_problematic")
    }
  }
}

deltas <- do.call(rbind, delta_rows)
write.table(deltas, file.path(out_dir, "tables", "delta_MMR_POLE_by_filter.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Improvement summary: drop_any vs all
wide <- merge(
  deltas[deltas$Filter == "all_samples",
         c("Tool", "COSMIC_Version", "Dataset", "delta_MMR", "p_MMR", "delta_POLE", "n_kept", "frac_discarded")],
  deltas[deltas$Filter == "drop_any_problematic",
         c("Tool", "COSMIC_Version", "Dataset", "delta_MMR", "p_MMR", "delta_POLE", "n_kept", "frac_discarded")],
  by = c("Tool", "COSMIC_Version", "Dataset"),
  suffixes = c("_all", "_filtered")
)
wide$delta_MMR_change <- wide$delta_MMR_filtered - wide$delta_MMR_all
wide$delta_POLE_change <- wide$delta_POLE_filtered - wide$delta_POLE_all
wide$improved_MMR <- wide$delta_MMR_change > 0
wide$improved_POLE <- wide$delta_POLE_change > 0
write.table(wide, file.path(out_dir, "tables", "filter_effect_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Flag prevalence
prev <- aggregate(cbind(flag_low_cosine, flag_low_mutations, flag_high_flat, problematic) ~
                    Tool + COSMIC_Version + Dataset,
                  data = flags, FUN = mean, na.rm = TRUE)
write.table(prev, file.path(out_dir, "tables", "flag_prevalence.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---------- Plots ----------
# 1) ΔMMR all vs filtered (focus preclinical, facet COSMIC)
pre <- wide[wide$Dataset == "Preclinical", ]
if (nrow(pre)) {
  long <- rbind(
    data.frame(Tool = pre$Tool, COSMIC_Version = pre$COSMIC_Version,
               Filter = "All samples", delta_MMR = pre$delta_MMR_all, stringsAsFactors = FALSE),
    data.frame(Tool = pre$Tool, COSMIC_Version = pre$COSMIC_Version,
               Filter = "Drop problematic", delta_MMR = pre$delta_MMR_filtered, stringsAsFactors = FALSE)
  )
  p1 <- ggplot(long, aes(x = Tool, y = delta_MMR, fill = Filter)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    facet_wrap(~COSMIC_Version, nrow = 2) +
    paper_theme(9, rotate_x = TRUE) +
    labs(title = "ΔMMR before vs after discarding problematic samples (Preclinical WES)",
         y = "ΔMMR (MSI − MSS)", x = NULL)
  save_paper_plot(p1, file.path(out_dir, "plots", "delta_MMR_all_vs_filtered_preclinical"),
                  width_cm = 16, height_cm = 10)
}

# 2) Change in ΔMMR after filtering
if (nrow(pre)) {
  p2 <- ggplot(pre, aes(x = Tool, y = delta_MMR_change, fill = COSMIC_Version)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    paper_theme(9, rotate_x = TRUE) +
    labs(title = "Change in ΔMMR after discarding problematic samples (Preclinical)",
         y = "ΔMMR_filtered − ΔMMR_all", x = NULL, fill = "COSMIC")
  save_paper_plot(p2, file.path(out_dir, "plots", "delta_MMR_change_preclinical"),
                  width_cm = 14, height_cm = 8)
}

# 3) Flag prevalence for MP × v3.2 as example
mp_flags <- flags[flags$Tool == "MutationalPatterns" & flags$COSMIC_Version == "v3.2", ]
if (nrow(mp_flags)) {
  tall <- data.frame(
    Dataset = rep(mp_flags$Dataset, 3),
    Rule = rep(c("Low cosine", "Low mutations", "High flat"), each = nrow(mp_flags)),
    Flagged = c(mp_flags$flag_low_cosine, mp_flags$flag_low_mutations, mp_flags$flag_high_flat)
  )
  agg <- aggregate(Flagged ~ Dataset + Rule, data = tall, FUN = mean)
  p3 <- ggplot(agg, aes(x = Rule, y = Flagged, fill = Dataset)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0),
                       labels = function(x) paste0(round(100 * x), "%")) +
    paper_theme(10) +
    labs(title = "Problematic-sample rule prevalence (MutationalPatterns · v3.2)",
         y = "Fraction flagged", x = NULL)
  save_paper_plot(p3, file.path(out_dir, "plots", "flag_prevalence_MP_v3.2"),
                  width_cm = 12, height_cm = 7)
}

# ---------- Optional survival (TCGA clinical) ----------
surv_path <- file.path(root, "metadata", "clinical_survival.tsv")
surv_note <- "No local survival table; MSI/POLE used as clinical endpoints."
if (file.exists(surv_path)) {
  surv <- read.table(surv_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  clin <- flags[flags$Dataset == "Clinical" & flags$Tool == "MutationalPatterns" &
                  flags$COSMIC_Version == "v3.2", ]
  spa <- flags[flags$Dataset == "Clinical" & flags$Tool == "SigProfilerAssignment" &
                 flags$COSMIC_Version == "v2", ]
  km_ok <- requireNamespace("survival", quietly = TRUE)

  run_surv <- function(flag_sub, label) {
    if (!nrow(flag_sub) || !nrow(surv)) return(NULL)
    m <- merge(flag_sub, surv, by = "Sample", all = FALSE)
    if (nrow(m) < 20) return(NULL)
    analyse <- function(dat, tag) {
      out <- data.frame(Setting = label, Cohort = tag, n = nrow(dat),
                        n_events = sum(dat$OS_event, na.rm = TRUE),
                        HR = NA_real_, p_cox = NA_real_, p_logrank = NA_real_,
                        stringsAsFactors = FALSE)
      if (nrow(dat) < 20 || sum(dat$OS_event, na.rm = TRUE) < 5) return(out)
      med <- median(dat$MMR_Contribution, na.rm = TRUE)
      dat$MMR_high <- as.integer(dat$MMR_Contribution >= med)
      if (km_ok) {
        sf <- survival::Surv(dat$OS_time, dat$OS_event)
        cx <- tryCatch(survival::coxph(sf ~ MMR_high, data = dat), error = function(e) NULL)
        if (!is.null(cx)) {
          s <- summary(cx)
          out$HR <- unname(exp(stats::coef(cx))[1])
          out$p_cox <- s$coefficients[1, "Pr(>|z|)"]
        }
        sd <- tryCatch(survival::survdiff(sf ~ MMR_high, data = dat), error = function(e) NULL)
        if (!is.null(sd)) out$p_logrank <- 1 - stats::pchisq(sd$chisq, length(sd$n) - 1)
      }
      out
    }
    rbind(analyse(m, "all_samples"),
          analyse(m[!m$problematic, , drop = FALSE], "drop_any_problematic"))
  }

  surv_tab <- rbind(run_surv(clin, "MutationalPatterns_v3.2"),
                    run_surv(spa, "SigProfilerAssignment_v2"))
  if (!is.null(surv_tab) && nrow(surv_tab)) {
    write.table(surv_tab, file.path(out_dir, "tables", "survival_MMR_high_vs_low.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    surv_note <- sprintf("TCGA OS n=%d (%d events); see survival_MMR_high_vs_low.tsv%s",
                         nrow(surv), sum(surv$OS_event),
                         if (!km_ok) " [install R 'survival' for Cox/log-rank]" else "")
  }
  writeLines(c("# Survival analysis", surv_note),
             file.path(out_dir, "tables", "survival_status.txt"))
} else {
  writeLines(c(
    "# Survival status",
    "",
    "No `metadata/clinical_survival.tsv` found.",
    "Clinical endpoints used: MSI/MSS and POLE (ΔMMR / ΔPOLE)."
  ), file.path(out_dir, "tables", "survival_status.txt"))
}

# ---------- Summary notes ----------
n_imp_mmr <- sum(wide$improved_MMR, na.rm = TRUE)
n_tot <- sum(is.finite(wide$delta_MMR_change))
med_change <- median(wide$delta_MMR_change[wide$Dataset == "Preclinical"], na.rm = TRUE)
frac_prob <- mean(flags$problematic, na.rm = TRUE)

notes <- c(
  "# Problematic-sample discard experiment",
  "",
  "## Rules",
  paste0("- Low cosine: Cosine_Similarity < ", cos_thr),
  "- Low mutations: Mutation_Count < Fig.6 threshold for Dataset×COSMIC (fallback 323/64)",
  paste0("- High flat: flat-signature contribution ≥ ", flat_cut),
  "- Combined filter `drop_any_problematic`: discard if any rule fires",
  "",
  "## Headline",
  sprintf("- Overall fraction flagged (any rule): %.1f%% of sample×tool×COSMIC rows", 100 * frac_prob),
  sprintf("- Preclinical median ΔMMR change after discard: %+.4f", med_change),
  sprintf("- Tool×COSMIC×Dataset blocks with ΔMMR increase: %d / %d", n_imp_mmr, n_tot),
  paste0("- Survival: ", surv_note),
  "",
  "## Outputs",
  "- tables/sample_problematic_flags.tsv",
  "- tables/delta_MMR_POLE_by_filter.tsv",
  "- tables/filter_effect_summary.tsv",
  "- tables/flag_prevalence.tsv",
  "- plots/delta_MMR_all_vs_filtered_preclinical.(png|svg)",
  "- plots/delta_MMR_change_preclinical.(png|svg)",
  "- plots/flag_prevalence_MP_v3.2.(png|svg)",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
)
writeLines(notes, file.path(out_dir, "README.md"))

message("Wrote problematic-discard results to ", out_dir)
message(sprintf("Flagged rows: %.1f%% | Preclinical median ΔMMR change: %+.4f | improved blocks: %d/%d",
                100 * frac_prob, med_change, n_imp_mmr, n_tot))
print(head(wide[order(-wide$delta_MMR_change),
                c("Tool", "COSMIC_Version", "Dataset", "frac_discarded_filtered",
                  "delta_MMR_all", "delta_MMR_filtered", "delta_MMR_change")], 12))

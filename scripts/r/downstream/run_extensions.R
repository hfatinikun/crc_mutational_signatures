#!/usr/bin/env Rscript
## Extensions the Battuello paper does not report.
## All four blocks read existing fits. Nothing is refit.
##
##   1. Paired assays on the cell lines present in WES, WGS and TSO-500.
##   2. ROC / AUC of MMR (and, descriptively, POLE) contribution.
##   3. Do tools rank the same tumours high for MMR?
##   4. SPA flat-cutoff sweep, v3.2–v3.6.
##   5. Does the MSI label itself associate with TCGA overall survival?
##
## Usage (from mywork/):
##   Rscript scripts/r/downstream/run_extensions.R

suppressPackageStartupMessages({
  library(jsonlite)
  library(ggplot2)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
source(file.path(dirname(script_dir), "helpers.R"))
source(file.path(script_dir, "plot_style.R"))

root <- .mywork_root()
cfg <- load_pipeline_config(root)
groups <- load_signature_groups(root)

out_dir <- file.path(root, "results", "downstream", "extensions")
plot_dir <- file.path(out_dir, "plots")
tab_dir <- file.path(out_dir, "tables")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

cosmic_levels <- c("v2", "v3.2", "v3.4", "v3.5", "v3.6", "TS")
assay_levels <- c("TSO-500", "WES", "WGS")

fmt4 <- function(x) ifelse(is.finite(x), sprintf("%.4f", x), "NA")

spearman_pair <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  out <- c(rho = NA_real_, p = NA_real_, n = sum(ok))
  if (sum(ok) < 5L) return(out)
  if (stats::sd(x[ok]) == 0 || stats::sd(y[ok]) == 0) return(out)
  ct <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  out["rho"] <- unname(ct$estimate)
  out["p"] <- ct$p.value
  out
}

paired_diff <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  d <- a[ok] - b[ok]
  p <- NA_real_
  if (sum(ok) >= 5L && any(d != 0)) {
    p <- suppressWarnings(stats::wilcox.test(a[ok], b[ok], paired = TRUE, exact = FALSE)$p.value)
  }
  c(median_diff = stats::median(d), p = p, n = sum(ok))
}

delta_mmr <- function(status, score) {
  msi <- score[status == "MSI"]
  mss <- score[status == "MSS"]
  msi <- msi[is.finite(msi)]
  mss <- mss[is.finite(mss)]
  if (!length(msi) || !length(mss)) return(NA_real_)
  stats::median(msi) - stats::median(mss)
}

median_in <- function(status, score, label) {
  v <- score[status == label]
  v <- v[is.finite(v)]
  if (!length(v)) return(NA_real_)
  stats::median(v)
}

mann_auc <- function(score, is_pos) {
  ok <- is.finite(score) & !is.na(is_pos)
  score <- score[ok]
  is_pos <- is_pos[ok]
  n1 <- sum(is_pos)
  n0 <- sum(!is_pos)
  if (n1 < 1L || n0 < 1L) return(NA_real_)
  if (length(unique(score)) < 2L) return(0.5)
  r <- rank(score, ties.method = "average")
  u <- sum(r[is_pos]) - n1 * (n1 + 1) / 2
  as.numeric(u / (n1 * n0))
}

roc_curve <- function(score, is_pos, tool) {
  ok <- is.finite(score) & !is.na(is_pos)
  score <- score[ok]
  is_pos <- is_pos[ok]
  ord <- order(score, decreasing = TRUE)
  score <- score[ord]
  is_pos <- is_pos[ord]
  n1 <- sum(is_pos)
  n0 <- sum(!is_pos)
  data.frame(
    Tool = tool,
    FPR = c(0, cumsum(!is_pos) / n0),
    TPR = c(0, cumsum(is_pos) / n1),
    stringsAsFactors = FALSE
  )
}

order_label <- function(named_vals) {
  named_vals <- named_vals[is.finite(named_vals)]
  if (length(named_vals) < 2L) return(NA_character_)
  paste(names(sort(named_vals, decreasing = TRUE)), collapse = " > ")
}

# =============================================================================
# 1. Paired assay comparison
# =============================================================================
ngs_path <- file.path(root, "results", "downstream", "ngs_workflow", "tables",
                      "ngs_sample_MMR_POLE_scores.tsv")
if (!file.exists(ngs_path)) stop("Missing ", ngs_path, " — run run_ngs_workflow.R first.")
ngs <- read.table(ngs_path, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)

metrics <- c("MMR_Contribution", "POLE_Contribution", "Flat_Contribution")
pair_rows <- list()
line_rows <- list()

for (tool in sort(unique(ngs$Tool))) {
  for (ck in intersect(cosmic_levels, unique(ngs$COSMIC_Version))) {
    sub <- ngs[ngs$Tool == tool & ngs$COSMIC_Version == ck, , drop = FALSE]
    by_assay <- lapply(assay_levels, function(a) sub[sub$Data_Type == a, , drop = FALSE])
    names(by_assay) <- assay_levels
    if (any(vapply(by_assay, nrow, 1L) == 0L)) next
    ids <- Reduce(intersect, lapply(by_assay, function(d) d$Sample))
    if (length(ids) < 10L) next

    wide <- by_assay[["WES"]][match(ids, by_assay[["WES"]]$Sample), c("Sample", "MS_Status", "POLE_Status")]
    for (a in assay_levels) {
      hit <- by_assay[[a]][match(ids, by_assay[[a]]$Sample), metrics, drop = FALSE]
      names(hit) <- paste(c("MMR", "POLE", "Flat"), a, sep = "_")
      wide <- cbind(wide, hit)
    }
    wide$Tool <- tool
    wide$COSMIC_Version <- ck
    line_rows[[length(line_rows) + 1L]] <- wide

    full_delta <- setNames(vapply(assay_levels, function(a) {
      delta_mmr(by_assay[[a]]$MS_Status, by_assay[[a]]$MMR_Contribution)
    }, numeric(1)), assay_levels)
    shared_delta <- setNames(vapply(assay_levels, function(a) {
      delta_mmr(wide$MS_Status, wide[[paste0("MMR_", a)]])
    }, numeric(1)), assay_levels)
    flat_msi <- setNames(vapply(assay_levels, function(a) {
      median_in(wide$MS_Status, wide[[paste0("Flat_", a)]], "MSI")
    }, numeric(1)), assay_levels)
    mmr_msi <- setNames(vapply(assay_levels, function(a) {
      median_in(wide$MS_Status, wide[[paste0("MMR_", a)]], "MSI")
    }, numeric(1)), assay_levels)

    msi <- wide$MS_Status == "MSI"
    sp_mmr_wg <- spearman_pair(wide$MMR_WES, wide$MMR_WGS)
    sp_mmr_tw <- spearman_pair(wide[["MMR_TSO-500"]], wide$MMR_WES)
    sp_flat_wg <- spearman_pair(wide$Flat_WES, wide$Flat_WGS)
    sp_pole_wg <- spearman_pair(wide$POLE_WES, wide$POLE_WGS)
    d_mmr <- paired_diff(wide$MMR_WGS, wide$MMR_WES)
    d_mmr_msi <- paired_diff(wide$MMR_WGS[msi], wide$MMR_WES[msi])
    d_flat <- paired_diff(wide$Flat_WGS, wide$Flat_WES)
    d_flat_msi <- paired_diff(wide$Flat_WGS[msi], wide$Flat_WES[msi])
    d_mmr_tso <- paired_diff(wide$`MMR_TSO-500`, wide$MMR_WES)
    d_mmr_tso_msi <- paired_diff(wide$`MMR_TSO-500`[msi], wide$MMR_WES[msi])

    pair_rows[[length(pair_rows) + 1L]] <- data.frame(
      Tool = tool,
      COSMIC_Version = ck,
      n_shared = length(ids),
      n_MSI = sum(msi),
      n_MSS = sum(wide$MS_Status == "MSS"),
      n_POLE_MUT = sum(wide$POLE_Status == "POLE MUT", na.rm = TRUE),
      delta_MMR_TSO_full = unname(full_delta["TSO-500"]),
      delta_MMR_WES_full = unname(full_delta["WES"]),
      delta_MMR_WGS_full = unname(full_delta["WGS"]),
      order_full = order_label(full_delta),
      delta_MMR_TSO_shared = unname(shared_delta["TSO-500"]),
      delta_MMR_WES_shared = unname(shared_delta["WES"]),
      delta_MMR_WGS_shared = unname(shared_delta["WGS"]),
      order_shared = order_label(shared_delta),
      same_order = identical(order_label(full_delta), order_label(shared_delta)),
      median_MMR_MSI_TSO = unname(mmr_msi["TSO-500"]),
      median_MMR_MSI_WES = unname(mmr_msi["WES"]),
      median_MMR_MSI_WGS = unname(mmr_msi["WGS"]),
      median_flat_MSI_TSO = unname(flat_msi["TSO-500"]),
      median_flat_MSI_WES = unname(flat_msi["WES"]),
      median_flat_MSI_WGS = unname(flat_msi["WGS"]),
      rho_MMR_WES_WGS = unname(sp_mmr_wg["rho"]),
      p_MMR_WES_WGS = unname(sp_mmr_wg["p"]),
      rho_MMR_TSO_WES = unname(sp_mmr_tw["rho"]),
      p_MMR_TSO_WES = unname(sp_mmr_tw["p"]),
      rho_POLE_WES_WGS = unname(sp_pole_wg["rho"]),
      rho_flat_WES_WGS = unname(sp_flat_wg["rho"]),
      p_flat_WES_WGS = unname(sp_flat_wg["p"]),
      median_MMR_WGS_minus_WES = unname(d_mmr["median_diff"]),
      p_MMR_WGS_vs_WES = unname(d_mmr["p"]),
      median_MMR_WGS_minus_WES_MSI = unname(d_mmr_msi["median_diff"]),
      p_MMR_WGS_vs_WES_MSI = unname(d_mmr_msi["p"]),
      median_flat_WGS_minus_WES = unname(d_flat["median_diff"]),
      p_flat_WGS_vs_WES = unname(d_flat["p"]),
      median_flat_WGS_minus_WES_MSI = unname(d_flat_msi["median_diff"]),
      p_flat_WGS_vs_WES_MSI = unname(d_flat_msi["p"]),
      median_MMR_TSO_minus_WES = unname(d_mmr_tso["median_diff"]),
      median_MMR_TSO_minus_WES_MSI = unname(d_mmr_tso_msi["median_diff"]),
      stringsAsFactors = FALSE
    )
  }
}

pair_df <- do.call(rbind, pair_rows)
line_df <- do.call(rbind, line_rows)
write.table(pair_df, file.path(tab_dir, "paired_assay_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(line_df, file.path(tab_dir, "paired_assay_samples.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

mp_pairs <- pair_df[pair_df$Tool == "MutationalPatterns", , drop = FALSE]
mp_long <- do.call(rbind, lapply(seq_len(nrow(mp_pairs)), function(i) {
  data.frame(
    COSMIC_Version = mp_pairs$COSMIC_Version[i],
    Assay = assay_levels,
    Cohort = rep(c("All samples of that assay", "Same 63 lines"), each = 3),
    delta_MMR = c(
      mp_pairs$delta_MMR_TSO_full[i], mp_pairs$delta_MMR_WES_full[i], mp_pairs$delta_MMR_WGS_full[i],
      mp_pairs$delta_MMR_TSO_shared[i], mp_pairs$delta_MMR_WES_shared[i], mp_pairs$delta_MMR_WGS_shared[i]
    ),
    stringsAsFactors = FALSE
  )
}))
mp_long$Assay <- factor(mp_long$Assay, levels = assay_levels)
mp_long$COSMIC_Version <- factor(mp_long$COSMIC_Version, levels = cosmic_levels)
p_delta <- ggplot(mp_long, aes(x = Assay, y = delta_MMR, fill = Cohort)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  facet_wrap(~COSMIC_Version, nrow = 1) +
  paper_theme(9, rotate_x = TRUE) +
  labs(title = "ΔMMR on every sample versus the 63 lines assayed three ways",
       subtitle = "MutationalPatterns. Same-line bars ask whether the ordering is the assay or the subset.",
       y = "ΔMMR", x = NULL, fill = NULL)
save_paper_plot(p_delta, file.path(plot_dir, "paired_delta_MMR_full_vs_shared_MP"),
                width_cm = 28, height_cm = 9)

focus_lines <- line_df[line_df$Tool == "MutationalPatterns" & line_df$COSMIC_Version == "v3.2" &
                         line_df$MS_Status == "MSI", , drop = FALSE]
spaghetti <- rbind(
  data.frame(Sample = focus_lines$Sample, Metric = "MMR contribution",
             Assay = "TSO-500", Value = focus_lines$`MMR_TSO-500`, stringsAsFactors = FALSE),
  data.frame(Sample = focus_lines$Sample, Metric = "MMR contribution",
             Assay = "WES", Value = focus_lines$MMR_WES, stringsAsFactors = FALSE),
  data.frame(Sample = focus_lines$Sample, Metric = "MMR contribution",
             Assay = "WGS", Value = focus_lines$MMR_WGS, stringsAsFactors = FALSE),
  data.frame(Sample = focus_lines$Sample, Metric = "Flat contribution",
             Assay = "TSO-500", Value = focus_lines$`Flat_TSO-500`, stringsAsFactors = FALSE),
  data.frame(Sample = focus_lines$Sample, Metric = "Flat contribution",
             Assay = "WES", Value = focus_lines$Flat_WES, stringsAsFactors = FALSE),
  data.frame(Sample = focus_lines$Sample, Metric = "Flat contribution",
             Assay = "WGS", Value = focus_lines$Flat_WGS, stringsAsFactors = FALSE)
)
spaghetti$Assay <- factor(spaghetti$Assay, levels = assay_levels)
p_sp <- ggplot(spaghetti, aes(x = Assay, y = Value, group = Sample)) +
  geom_line(colour = "salmon2", alpha = 0.7) +
  geom_point(colour = "salmon2", size = 1.4) +
  facet_wrap(~Metric) +
  scale_y_unit() +
  paper_theme(10) +
  labs(title = "Same MSI lines across assays (MutationalPatterns · v3.2)",
       subtitle = paste0(nrow(focus_lines), " MSI lines that have WES, WGS and TSO-500"),
       y = "Contribution", x = NULL)
save_paper_plot(p_sp, file.path(plot_dir, "paired_MSI_lines_MP_v32"),
                width_cm = 18, height_cm = 9)

# =============================================================================
# 2. Classification: AUC of contribution versus the label
# =============================================================================
score_path <- file.path(root, "results", "downstream", "tables", "sample_MMR_POLE_scores.tsv")
if (!file.exists(score_path)) stop("Missing ", score_path)
scores <- read.table(score_path, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
delta_path <- file.path(root, "results", "downstream", "tables", "delta_MMR_POLE.tsv")
deltas <- read.table(delta_path, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)

wes <- scores[scores$Dataset %in% c("Preclinical", "Clinical"), , drop = FALSE]
auc_rows <- list()
roc_focus <- list()
for (ds in c("Preclinical", "Clinical")) {
  for (tool in sort(unique(wes$Tool))) {
    for (ck in intersect(cosmic_levels, unique(wes$COSMIC_Version))) {
      sub <- wes[wes$Dataset == ds & wes$Tool == tool & wes$COSMIC_Version == ck, , drop = FALSE]
      if (!nrow(sub)) next
      sub <- sub[!is.na(sub$MS_Status) & is.finite(sub$MMR_Contribution), , drop = FALSE]
      is_msi <- sub$MS_Status == "MSI"
      auc_mmr <- mann_auc(sub$MMR_Contribution, is_msi)
      pole <- sub[sub$POLE_Status %in% c("POLE MUT", "POLE WT") & is.finite(sub$POLE_Contribution), , drop = FALSE]
      auc_pole <- mann_auc(pole$POLE_Contribution, pole$POLE_Status == "POLE MUT")
      dhit <- deltas[deltas$Dataset == ds & deltas$Tool == tool & deltas$COSMIC_Version == ck, , drop = FALSE]
      auc_rows[[length(auc_rows) + 1L]] <- data.frame(
        Dataset = ds, Tool = tool, COSMIC_Version = ck,
        n = nrow(sub), n_MSI = sum(is_msi), n_MSS = sum(!is_msi),
        delta_MMR = if (nrow(dhit)) dhit$delta_MMR[1] else NA_real_,
        AUC_MMR = auc_mmr,
        n_POLE_MUT = sum(pole$POLE_Status == "POLE MUT"),
        n_POLE_WT_MSS = sum(pole$POLE_Status == "POLE WT"),
        delta_POLE = if (nrow(dhit)) dhit$delta_POLE[1] else NA_real_,
        AUC_POLE = auc_pole,
        stringsAsFactors = FALSE
      )
      if (ds == "Preclinical" && ck == "v2" &&
          tool %in% c("MutationalPatterns", "SigProfilerAssignment", "sigLASSO")) {
        roc_focus[[length(roc_focus) + 1L]] <- roc_curve(sub$MMR_Contribution, is_msi, tool)
      }
    }
  }
}
auc_df <- do.call(rbind, auc_rows)
auc_df$rank_delta <- ave(-auc_df$delta_MMR, auc_df$Dataset, auc_df$COSMIC_Version,
                         FUN = function(z) rank(z, ties.method = "min", na.last = "keep"))
auc_df$rank_AUC <- ave(-auc_df$AUC_MMR, auc_df$Dataset, auc_df$COSMIC_Version,
                       FUN = function(z) rank(z, ties.method = "min", na.last = "keep"))
auc_df$rank_gap <- auc_df$rank_AUC - auc_df$rank_delta
write.table(auc_df, file.path(tab_dir, "auc_vs_delta.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

pre_v2 <- auc_df[auc_df$Dataset == "Preclinical" & auc_df$COSMIC_Version == "v2", , drop = FALSE]
p_auc <- ggplot(pre_v2, aes(x = delta_MMR, y = AUC_MMR, label = Tool)) +
  geom_point(size = 2.4, colour = "firebrick3") +
  geom_text(vjust = -0.6, size = 2.6) +
  paper_theme(10) +
  labs(title = "Cohort gap versus classification (preclinical WES · COSMIC v2)",
       subtitle = "AUC is the chance a random MSI line outranks a random MSS line on MMR weight",
       x = "ΔMMR", y = "AUC")
save_paper_plot(p_auc, file.path(plot_dir, "auc_vs_delta_preclinical_v2"),
                width_cm = 16, height_cm = 12)

if (length(roc_focus)) {
  roc_df <- do.call(rbind, roc_focus)
  p_roc <- ggplot(roc_df, aes(x = FPR, y = TPR, colour = Tool)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    geom_path(linewidth = 0.7) +
    paper_theme(10) +
    labs(title = "ROC of MMR contribution for the MSI label",
         subtitle = "Preclinical WES · COSMIC v2",
         x = "False positive rate", y = "True positive rate", colour = NULL)
  save_paper_plot(p_roc, file.path(plot_dir, "roc_MMR_preclinical_v2"),
                  width_cm = 14, height_cm = 12)
}

# =============================================================================
# 3. Do tools agree on which samples are high-MMR?
# =============================================================================
agree_rows <- list()
discord_rows <- list()
discord_samples <- list()
pre <- wes[wes$Dataset == "Preclinical", , drop = FALSE]
for (ck in intersect(cosmic_levels, unique(pre$COSMIC_Version))) {
  sub <- pre[pre$COSMIC_Version == ck & is.finite(pre$MMR_Contribution), , drop = FALSE]
  tools_here <- sort(unique(sub$Tool))
  wide <- reshape(
    sub[, c("Sample", "MS_Status", "Tool", "MMR_Contribution")],
    idvar = c("Sample", "MS_Status"),
    timevar = "Tool",
    direction = "wide"
  )
  names(wide) <- sub("^MMR_Contribution\\.", "", names(wide))
  for (i in seq_along(tools_here)) {
    for (j in seq_along(tools_here)) {
      if (j <= i) next
      a <- tools_here[i]
      b <- tools_here[j]
      if (!all(c(a, b) %in% names(wide))) next
      sp <- spearman_pair(wide[[a]], wide[[b]])
      agree_rows[[length(agree_rows) + 1L]] <- data.frame(
        Dataset = "Preclinical", COSMIC_Version = ck,
        Tool_A = a, Tool_B = b,
        n = unname(sp["n"]), rho = unname(sp["rho"]), p = unname(sp["p"]),
        stringsAsFactors = FALSE
      )
    }
  }
  need <- c("MutationalPatterns", "SigProfilerAssignment")
  if (all(need %in% names(wide))) {
    mp <- wide$MutationalPatterns
    spa <- wide$SigProfilerAssignment
    ok <- is.finite(mp) & is.finite(spa)
    top_spa <- spa >= stats::quantile(spa[ok], 0.75)
    bot_mp <- mp <= stats::quantile(mp[ok], 0.25)
    top_mp <- mp >= stats::quantile(mp[ok], 0.75)
    bot_spa <- spa <= stats::quantile(spa[ok], 0.25)
    both_top <- top_spa & top_mp & ok
    discord_rows[[length(discord_rows) + 1L]] <- data.frame(
      COSMIC_Version = ck,
      n = sum(ok),
      n_top_SPA_and_bottom_MP = sum(top_spa & bot_mp & ok, na.rm = TRUE),
      n_top_MP_and_bottom_SPA = sum(top_mp & bot_spa & ok, na.rm = TRUE),
      n_both_top_quartile = sum(both_top, na.rm = TRUE),
      n_top_quartile = sum(top_mp & ok, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    if (ck == "v2") {
      hit <- wide[which(top_spa & bot_mp & ok), c("Sample", "MS_Status")]
      if (nrow(hit)) {
        hit$MMR_MutationalPatterns <- mp[top_spa & bot_mp & ok]
        hit$MMR_SigProfilerAssignment <- spa[top_spa & bot_mp & ok]
        hit$COSMIC_Version <- ck
        discord_samples[[length(discord_samples) + 1L]] <- hit
      }
    }
  }
}
agree_df <- do.call(rbind, agree_rows)
discord_df <- do.call(rbind, discord_rows)
write.table(agree_df, file.path(tab_dir, "tool_agreement_spearman.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(discord_df, file.path(tab_dir, "tool_quartile_discordance.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
if (length(discord_samples)) {
  write.table(do.call(rbind, discord_samples),
              file.path(tab_dir, "discordant_samples_SPA_top_MP_bottom_v2.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

ctrl_tools <- c("MutationalPatterns", "sigminer", "YAPSA", "SigProfilerAssignment", "sigLASSO")
heat <- agree_df[agree_df$COSMIC_Version == "v2" &
                   agree_df$Tool_A %in% ctrl_tools & agree_df$Tool_B %in% ctrl_tools, , drop = FALSE]
if (nrow(heat)) {
  heat_sym <- rbind(
    heat,
    data.frame(Dataset = heat$Dataset, COSMIC_Version = heat$COSMIC_Version,
               Tool_A = heat$Tool_B, Tool_B = heat$Tool_A,
               n = heat$n, rho = heat$rho, p = heat$p, stringsAsFactors = FALSE),
    data.frame(Dataset = "Preclinical", COSMIC_Version = "v2",
               Tool_A = ctrl_tools, Tool_B = ctrl_tools,
               n = NA_integer_, rho = 1, p = NA_real_, stringsAsFactors = FALSE)
  )
  heat_sym$Tool_A <- factor(heat_sym$Tool_A, levels = ctrl_tools)
  heat_sym$Tool_B <- factor(heat_sym$Tool_B, levels = ctrl_tools)
  p_heat <- ggplot(heat_sym, aes(x = Tool_A, y = Tool_B, fill = rho)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%.2f", rho)), size = 3) +
    scale_fill_gradient(limits = c(0, 1), low = "grey90", high = "firebrick3", na.value = "grey90") +
    paper_theme(9, rotate_x = TRUE) +
    labs(title = "Spearman correlation of per-sample MMR weight",
         subtitle = "Preclinical WES · COSMIC v2. sigminer and YAPSA are the negative control.",
         x = NULL, y = NULL, fill = "rho")
  save_paper_plot(p_heat, file.path(plot_dir, "tool_agreement_spearman_v2"),
                  width_cm = 16, height_cm = 14)
}

# =============================================================================
# 4. Flat-cutoff sweep for SPA, COSMIC v3.2–v3.6
# =============================================================================
sum_group <- function(norm_mat, sigs) {
  sigs <- intersect(sigs, rownames(norm_mat))
  if (!length(sigs)) return(setNames(rep(0, ncol(norm_mat)), colnames(norm_mat)))
  if (length(sigs) == 1L) setNames(as.numeric(norm_mat[sigs, ]), colnames(norm_mat)) else
    colSums(norm_mat[sigs, , drop = FALSE])
}

pheno <- read.table(file.path(root, "metadata", "sample_phenotypes.tsv"),
                    header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
cuts <- seq(0.05, 0.50, by = 0.05)
sweep_rows <- list()
for (ck in c("v3.2", "v3.4", "v3.5", "v3.6")) {
  for (ds in c("Preclinical", "Clinical")) {
    nf <- file.path(root, "results", "contributions", "SigProfilerAssignment", ck,
                    paste0(ds, "_contribution_normalized.tsv"))
    if (!file.exists(nf)) next
    mat <- as.matrix(read.table(nf, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE))
    ph <- pheno[pheno$Dataset == ds & pheno$Sample %in% colnames(mat), , drop = FALSE]
    if (!nrow(ph)) next
    mmr <- sum_group(mat, resolve_signature_group("MMR", ck, rownames(mat), groups))
    flat <- sum_group(mat, resolve_signature_group("flat", ck, rownames(mat), groups))
    ph$MMR <- unname(mmr[ph$Sample])
    ph$Flat <- unname(flat[ph$Sample])
    ph$MS <- ifelse(ph$MSI_status == "MSI", "MSI", "MSS")
    n_all <- nrow(ph)
    base_delta <- delta_mmr(ph$MS, ph$MMR)
    sweep_rows[[length(sweep_rows) + 1L]] <- data.frame(
      Tool = "SigProfilerAssignment", COSMIC_Version = ck, Dataset = ds,
      flat_cut = NA_real_, rule = "no flat filter",
      n_kept = n_all, frac_kept = 1,
      n_MSI_kept = sum(ph$MS == "MSI"), n_MSS_kept = sum(ph$MS == "MSS"),
      delta_MMR = base_delta, stringsAsFactors = FALSE
    )
    for (cut in cuts) {
      keep <- is.finite(ph$Flat) & ph$Flat < cut
      kept <- ph[keep, , drop = FALSE]
      sweep_rows[[length(sweep_rows) + 1L]] <- data.frame(
        Tool = "SigProfilerAssignment", COSMIC_Version = ck, Dataset = ds,
        flat_cut = cut, rule = "drop flat >= cut",
        n_kept = nrow(kept), frac_kept = nrow(kept) / n_all,
        n_MSI_kept = sum(kept$MS == "MSI"), n_MSS_kept = sum(kept$MS == "MSS"),
        delta_MMR = delta_mmr(kept$MS, kept$MMR),
        stringsAsFactors = FALSE
      )
    }
  }
}
sweep_df <- do.call(rbind, sweep_rows)
write.table(sweep_df, file.path(tab_dir, "spa_flat_cutoff_sweep.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

sweep_plot <- sweep_df[sweep_df$rule == "drop flat >= cut" & is.finite(sweep_df$delta_MMR), , drop = FALSE]
sweep_plot$COSMIC_Version <- factor(sweep_plot$COSMIC_Version, levels = c("v3.2", "v3.4", "v3.5", "v3.6"))
mark <- sweep_plot[abs(sweep_plot$flat_cut - 0.25) < 1e-9, , drop = FALSE]
p_sw <- ggplot(sweep_plot, aes(x = frac_kept, y = delta_MMR, colour = COSMIC_Version)) +
  geom_path(linewidth = 0.7) +
  geom_point(size = 1.6) +
  geom_point(data = mark, shape = 1, size = 3.5, stroke = 0.7) +
  facet_wrap(~Dataset) +
  scale_x_continuous(limits = c(0, 1), labels = function(x) paste0(round(100 * x), "%")) +
  scale_colour_manual(values = paper_cols$reference[c("v3.2", "v3.4", "v3.5", "v3.6")]) +
  paper_theme(10) +
  labs(title = "SPA: ΔMMR after dropping high flat load",
       subtitle = "Cuts 0.05 to 0.50. Open circles are the thesis cut of 0.25.",
       x = "Samples kept", y = "ΔMMR", colour = "COSMIC")
save_paper_plot(p_sw, file.path(plot_dir, "spa_flat_cutoff_sweep"),
                width_cm = 18, height_cm = 9)

# =============================================================================
# 5. MSI label versus TCGA overall survival
# =============================================================================
surv_path <- file.path(root, "metadata", "clinical_survival.tsv")
surv_df <- NULL
if (file.exists(surv_path) && requireNamespace("survival", quietly = TRUE)) {
  surv <- read.table(surv_path, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  clin <- pheno[pheno$Dataset == "Clinical", c("Sample", "MSI_status", "POLE_status")]
  dat <- merge(clin, surv, by = "Sample", all = FALSE)
  dat <- dat[is.finite(dat$OS_time) & !is.na(dat$OS_event), , drop = FALSE]
  dat$MSI_group <- ifelse(dat$MSI_status == "MSI", "MSI", "MSS")
  dat$MSI_bin <- as.integer(dat$MSI_group == "MSI")
  sf <- survival::Surv(dat$OS_time, dat$OS_event)
  cx <- survival::coxph(sf ~ MSI_bin, data = dat)
  cx_s <- summary(cx)
  sd <- survival::survdiff(sf ~ MSI_group, data = dat)
  p_lr <- 1 - stats::pchisq(sd$chisq, length(sd$n) - 1)
  arm <- function(g) {
    x <- dat[dat$MSI_group == g, , drop = FALSE]
    data.frame(n = nrow(x), events = sum(x$OS_event))
  }
  a_msi <- arm("MSI")
  a_mss <- arm("MSS")
  surv_df <- data.frame(
    contrast = "MSI vs MSS (label, not signature)",
    n = nrow(dat),
    n_MSI = a_msi$n, events_MSI = a_msi$events,
    n_MSS = a_mss$n, events_MSS = a_mss$events,
    HR_MSI = unname(exp(stats::coef(cx))[1]),
    p_cox = cx_s$coefficients[1, "Pr(>|z|)"],
    p_logrank = p_lr,
    stringsAsFactors = FALSE
  )
  write.table(surv_df, file.path(tab_dir, "survival_MSI_label.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  fit <- survival::survfit(survival::Surv(OS_time, OS_event) ~ MSI_group, data = dat)
  strata_n <- fit$strata
  km <- data.frame(
    time = fit$time,
    surv = fit$surv,
    group = rep(sub("^MSI_group=", "", names(strata_n)), times = strata_n),
    stringsAsFactors = FALSE
  )
  km <- rbind(
    data.frame(time = 0, surv = 1, group = unique(km$group), stringsAsFactors = FALSE),
    km
  )
  p_km <- ggplot(km, aes(x = time, y = surv, colour = group)) +
    geom_step(linewidth = 0.7) +
    scale_colour_manual(values = c(MSS = "lightblue1", MSI = "salmon2")) +
    scale_y_continuous(limits = c(0, 1)) +
    paper_theme(10) +
    labs(title = "TCGA overall survival by the MSI label",
         subtitle = sprintf("n = %d, events = %d. Cox HR (MSI vs MSS) = %.2f, p = %.3f",
                            nrow(dat), sum(dat$OS_event), surv_df$HR_MSI, surv_df$p_cox),
         x = "Time (days)", y = "Overall survival", colour = NULL)
  save_paper_plot(p_km, file.path(plot_dir, "km_MSI_label_OS"),
                  width_cm = 14, height_cm = 10)
} else if (!file.exists(surv_path)) {
  message("No clinical_survival.tsv; skipped the MSI survival check.")
} else {
  message("R package 'survival' is not installed; skipped the MSI survival check.")
}

# =============================================================================
# Notes
# =============================================================================
mp32 <- pair_df[pair_df$Tool == "MutationalPatterns" & pair_df$COSMIC_Version == "v3.2", , drop = FALSE]
pre_v2 <- pre_v2[order(-pre_v2$AUC_MMR), ]
ctrl <- agree_df[agree_df$COSMIC_Version == "v2" & agree_df$Tool_A == "MutationalPatterns" &
                   agree_df$Tool_B %in% c("sigminer", "YAPSA", "SigProfilerAssignment", "sigLASSO"), , drop = FALSE]
disc_v2 <- discord_df[discord_df$COSMIC_Version == "v2", , drop = FALSE]
sw25 <- sweep_df[sweep_df$Dataset == "Preclinical" & is.finite(sweep_df$flat_cut) &
                   abs(sweep_df$flat_cut - 0.25) < 1e-9, , drop = FALSE]
n_same <- sum(pair_df$same_order, na.rm = TRUE)
line_fmt <- function(df, fun) {
  if (!nrow(df)) return(character())
  vapply(seq_len(nrow(df)), function(i) fun(df[i, , drop = FALSE]), character(1))
}

notes <- c(
  "# Extensions beyond Battuello et al.",
  "",
  "## 1. Paired assays",
  sprintf("- Shared lines: %d (%d MSI, %d MSS, %d POLE-mutant).",
          mp32$n_shared, mp32$n_MSI, mp32$n_MSS, mp32$n_POLE_MUT),
  sprintf("- MutationalPatterns · v3.2, all samples: %s (Δ %.3f / %.3f / %.3f).",
          mp32$order_full, mp32$delta_MMR_TSO_full, mp32$delta_MMR_WES_full, mp32$delta_MMR_WGS_full),
  sprintf("- Same tool and catalogue, shared 63: %s (Δ %.3f / %.3f / %.3f).",
          mp32$order_shared, mp32$delta_MMR_TSO_shared, mp32$delta_MMR_WES_shared, mp32$delta_MMR_WGS_shared),
  sprintf("- Within MSI lines, median MMR(WGS) − MMR(WES) = %.3f (paired p = %.3g); median flat difference = %.3f (p = %.3g).",
          mp32$median_MMR_WGS_minus_WES_MSI, mp32$p_MMR_WGS_vs_WES_MSI,
          mp32$median_flat_WGS_minus_WES_MSI, mp32$p_flat_WGS_vs_WES_MSI),
  sprintf("- Spearman MMR, WES vs WGS, all 63 lines: rho = %.3f.", mp32$rho_MMR_WES_WGS),
  sprintf("- Tool × catalogue blocks whose full-cohort order matches the shared-63 order: %d / %d.",
          n_same, nrow(pair_df)),
  "",
  "## 2. AUC versus ΔMMR (preclinical WES · v2, highest AUC first)",
  line_fmt(pre_v2, function(r) {
    sprintf("- %s: ΔMMR %s (rank %d), AUC %s (rank %d).",
            r$Tool, fmt4(r$delta_MMR), r$rank_delta, fmt4(r$AUC_MMR), r$rank_AUC)
  }),
  "",
  "## 3. Tool agreement on preclinical WES · v2",
  line_fmt(ctrl, function(r) {
    sprintf("- Spearman MMR, MutationalPatterns vs %s: rho = %s.", r$Tool_B, fmt4(r$rho))
  }),
  sprintf("- Samples in the top SPA quartile and the bottom MutationalPatterns quartile: %s.",
          if (nrow(disc_v2)) disc_v2$n_top_SPA_and_bottom_MP else "NA"),
  sprintf("- Samples in both top quartiles: %s of %s.",
          if (nrow(disc_v2)) disc_v2$n_both_top_quartile else "NA",
          if (nrow(disc_v2)) disc_v2$n_top_quartile else "NA"),
  "",
  "## 4. SPA flat cut 0.25, preclinical",
  line_fmt(sw25, function(r) {
    sprintf("- %s: keep %.0f%%, ΔMMR %s.", r$COSMIC_Version, 100 * r$frac_kept, fmt4(r$delta_MMR))
  }),
  "",
  "## 5. MSI label versus overall survival",
  if (is.null(surv_df)) {
    "- Not run."
  } else {
    sprintf("- n = %d (%d MSI / %d MSS), events = %d (%d / %d). HR = %.2f, Cox p = %.3f, log-rank p = %.3f.",
            surv_df$n, surv_df$n_MSI, surv_df$n_MSS,
            surv_df$events_MSI + surv_df$events_MSS, surv_df$events_MSI, surv_df$events_MSS,
            surv_df$HR_MSI, surv_df$p_cox, surv_df$p_logrank)
  },
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
)
writeLines(notes, file.path(out_dir, "README.md"))
message(paste(notes, collapse = "\n"))
message("Wrote ", out_dir)

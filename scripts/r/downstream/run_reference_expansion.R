#!/usr/bin/env Rscript
## Paper Fig. 5-style reference comparison: COSMIC v2 / v3.2 / tissue-specific (TS).
## Requires TS fitting outputs already present (run_tissue_specific.R).

suppressPackageStartupMessages({
  library(ggplot2)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
source(file.path(dirname(script_dir), "helpers.R"))
source(file.path(script_dir, "plot_style.R"))

root <- .mywork_root()
out_dir <- file.path(root, "results", "downstream", "reference_expansion")
dir.create(file.path(out_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

refs <- c("v2", "v3.2", "TS")
ref_labels <- c(v2 = "COSMIC v2", `v3.2` = "COSMIC v3.2", TS = "Tissue-specific")
ref_cols <- c(`COSMIC v2` = "#1B8148", `COSMIC v3.2` = "Orange", `Tissue-specific` = "#3852BD")

cos_path <- file.path(root, "results", "cosine", "CosineResults_AllTools_AllCOSMIC.tsv")
if (!file.exists(cos_path)) stop("Missing combined cosine; run 00_combine_cosine.R")
cos_df <- read.table(cos_path, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)

delta_path <- file.path(root, "results", "downstream", "tables", "delta_MMR_POLE.tsv")
if (!file.exists(delta_path)) stop("Missing delta_MMR_POLE.tsv; run run_downstream.R")
delta <- read.table(delta_path, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)

comp_path <- file.path(root, "results", "downstream", "tables", "MMR_composition_MSI.tsv")
comp <- if (file.exists(comp_path)) {
  read.table(comp_path, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
} else NULL

sample_path <- file.path(root, "results", "downstream", "tables", "sample_MMR_POLE_scores.tsv")
samp <- if (file.exists(sample_path)) {
  read.table(sample_path, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
} else NULL

# Cosine: MutationalPatterns × WES × refs
cos_sub <- cos_df[cos_df$Tool == "MutationalPatterns" &
                    cos_df$Dataset %in% c("Preclinical", "Clinical") &
                    cos_df$COSMIC_Version %in% refs &
                    (is.na(cos_df$Data_Type) | cos_df$Data_Type == "WES"), ]
if (!nrow(cos_sub)) stop("No MP×TS/v2/v3.2 WES cosine rows — run tissue-specific fitting + combine.")
cos_sub$Reference <- factor(ref_labels[cos_sub$COSMIC_Version], levels = unname(ref_labels))

p1 <- ggplot(cos_sub, aes(x = Dataset, y = Cosine_Similarity, colour = Reference)) +
  geom_boxplot(linewidth = 0.9) +
  scale_colour_manual(values = ref_cols) +
  geom_hline(yintercept = 0.9, linetype = "dashed") +
  coord_flip() +
  scale_y_continuous(breaks = seq(0.825, 1, 0.05), limits = c(0.825, 1), expand = c(0, 0)) +
  paper_theme(10) +
  labs(title = "Cosine by reference (MutationalPatterns · WES)",
       y = "Cosine Similarity", x = NULL, colour = "Reference Dataset")
save_paper_plot(p1, file.path(out_dir, "plots", "ref_cosine_v2_v32_TS"), width_cm = 12, height_cm = 7)

# ΔMMR / ΔPOLE summary table for refs
dsub <- delta[delta$Tool == "MutationalPatterns" &
                delta$Dataset %in% c("Preclinical", "Clinical") &
                delta$COSMIC_Version %in% refs &
                (is.na(delta$Data_Type) | delta$Data_Type == "WES"), ]
dsub$Reference <- ref_labels[dsub$COSMIC_Version]
write.table(dsub, file.path(out_dir, "tables", "ref_delta_MMR_POLE.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Paper-reported preclinical ΔMMR: C2=0.34, C3=0.26, TS=0.27
paper_cmp <- data.frame(
  Dataset = "Preclinical",
  COSMIC_Version = c("v2", "v3.2", "TS"),
  Metric = "delta_MMR",
  Paper = c(0.34, 0.26, 0.27),
  stringsAsFactors = FALSE
)
paper_cmp$Replicate <- vapply(paper_cmp$COSMIC_Version, function(cv) {
  v <- dsub$delta_MMR[dsub$Dataset == "Preclinical" & dsub$COSMIC_Version == cv]
  if (length(v)) v[[1]] else NA_real_
}, numeric(1))
paper_cmp$Abs_diff <- abs(paper_cmp$Replicate - paper_cmp$Paper)
paper_cmp$Match <- paper_cmp$Abs_diff < 0.03
write.table(paper_cmp, file.path(out_dir, "tables", "ref_paper_vs_replicate.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# MMR contribution boxplots by reference (preclinical)
if (!is.null(samp)) {
  s <- samp[samp$Tool == "MutationalPatterns" & samp$Dataset == "Preclinical" &
              samp$COSMIC_Version %in% refs & !is.na(samp$MS_Status), ]
  if (nrow(s)) {
    s$Reference <- factor(ref_labels[s$COSMIC_Version], levels = unname(ref_labels))
    s$MS_Status <- factor(s$MS_Status, levels = c("MSS", "MSI"))
    p2 <- ggplot(s, aes(x = Reference, y = MMR_Contribution, fill = MS_Status)) +
      geom_boxplot(outlier.size = 0.4) +
      scale_fill_manual(values = paper_cols$msi_mss) +
      scale_y_unit() +
      paper_theme(10) +
      labs(title = "MMR contribution by reference (MP · Preclinical WES)",
           x = "Reference Dataset", y = "Overall MMR Signatures Contribution",
           fill = "Microsatellite Status")
    cmp <- maybe_compare_means(method = "wilcox.test", label = "p.signif")
    if (!is.null(cmp)) p2 <- p2 + cmp
    save_paper_plot(p2, file.path(out_dir, "plots", "ref_MMR_boxplot_preclinical"),
                    width_cm = 12, height_cm = 8)
  }
}

# MMR composition stacked bars
if (!is.null(comp)) {
  csub <- comp[comp$Tool == "MutationalPatterns" & comp$Dataset == "Preclinical" &
                 comp$COSMIC_Version %in% refs, ]
  if (nrow(csub)) {
    csub$Reference <- factor(ref_labels[csub$COSMIC_Version], levels = unname(ref_labels))
    fill_vals <- paper_cols$mmr_sbs
    miss <- setdiff(unique(csub$Signature), names(fill_vals))
    if (length(miss)) {
      fill_vals <- c(fill_vals, setNames(grDevices::hcl.colors(length(miss), "Dark 3"), miss))
    }
    p3 <- ggplot(csub, aes(x = Reference, y = Mean_Contribution, fill = Signature)) +
      geom_col(width = 0.7) +
      scale_fill_manual(values = fill_vals) +
      scale_y_unit() +
      paper_theme(10) +
      labs(title = "MMR signature mix in MSI by reference (MP · Preclinical)",
           x = "Reference Dataset", y = "Relative contribution within MMR set",
           fill = "Signature") +
      guides(fill = guide_legend(ncol = 2, byrow = TRUE, title.position = "top"))
    save_paper_plot(p3, file.path(out_dir, "plots", "ref_MMR_composition_v2_v32_TS"),
                    width_cm = 12, height_cm = 9)
  }
}

message("Reference expansion written to ", out_dir)
print(paper_cmp)

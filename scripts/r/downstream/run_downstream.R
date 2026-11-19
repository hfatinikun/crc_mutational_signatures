#!/usr/bin/env Rscript
## Algorithm + Reference downstream — paper-style figures (Battuello et al. bbae249).

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

out_dir <- file.path(root, "results", "downstream")
dir.create(file.path(out_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
paper_dir <- file.path(out_dir, "plots", "paper_style")
dir.create(paper_dir, recursive = TRUE, showWarnings = FALSE)

cos_path <- file.path(root, "results", "cosine", "CosineResults_AllTools_AllCOSMIC.tsv")
if (!file.exists(cos_path)) {
  message("Combined cosine file missing; running 00_combine_cosine.R")
  system2("Rscript", file.path(dirname(script_dir), "00_combine_cosine.R"), stdout = TRUE, stderr = TRUE)
}
cos_df <- read.table(cos_path, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
cos_df$Assigned <- cos_df$Cosine_Similarity >= cfg$cosine_threshold

pheno <- read.table(file.path(root, "metadata", "sample_phenotypes.tsv"),
                    header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

# ---------- Overview grids (supplement) ----------
p1 <- ggplot(cos_df, aes(x = Tool, y = Cosine_Similarity, fill = COSMIC_Version)) +
  geom_boxplot(outlier.size = 0.4, position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = cfg$cosine_threshold, linetype = "dashed") +
  facet_wrap(~Dataset) +
  paper_theme(11, rotate_x = TRUE) +
  labs(title = "Cosine similarity by tool and COSMIC version",
       y = "Cosine similarity", x = NULL)
ggsave(file.path(out_dir, "plots", "cosine_by_tool_cosmic.png"), p1, width = 12, height = 6, dpi = 150)

assigned_tab <- aggregate(Assigned ~ Tool + COSMIC_Version + Dataset, data = cos_df,
                          FUN = function(x) c(n = length(x), assigned = sum(x), frac = mean(x)))
assigned_df <- do.call(data.frame, assigned_tab)
names(assigned_df) <- c("Tool", "COSMIC_Version", "Dataset", "n", "n_assigned", "frac_assigned")
write.table(assigned_df, file.path(out_dir, "tables", "assigned_fraction.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

wilcox_rows <- list()
for (ds in unique(cos_df$Dataset)) {
  for (cv in unique(cos_df$COSMIC_Version)) {
    sub <- cos_df[cos_df$Dataset == ds & cos_df$COSMIC_Version == cv, ]
    tools <- sort(unique(sub$Tool))
    if (length(tools) < 2) next
    for (i in seq_len(length(tools) - 1)) {
      for (j in (i + 1):length(tools)) {
        a <- sub$Cosine_Similarity[sub$Tool == tools[i]]
        b <- sub$Cosine_Similarity[sub$Tool == tools[j]]
        if (length(a) < 2 || length(b) < 2) next
        wt <- suppressWarnings(wilcox.test(a, b))
        wilcox_rows[[length(wilcox_rows) + 1]] <- data.frame(
          Dataset = ds, COSMIC_Version = cv,
          Tool_A = tools[i], Tool_B = tools[j],
          p_value = wt$p.value,
          median_A = median(a, na.rm = TRUE),
          median_B = median(b, na.rm = TRUE),
          stringsAsFactors = FALSE
        )
      }
    }
  }
}
if (length(wilcox_rows)) {
  write.table(do.call(rbind, wilcox_rows),
              file.path(out_dir, "tables", "cosine_wilcoxon_tool_pairs.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

# ---------- Contributions: aggregates + sample-level ----------
sum_group <- function(norm_mat, sigs) {
  sigs <- intersect(sigs, rownames(norm_mat))
  if (!length(sigs)) return(setNames(rep(0, ncol(norm_mat)), colnames(norm_mat)))
  if (length(sigs) == 1) setNames(as.numeric(norm_mat[sigs, ]), colnames(norm_mat)) else colSums(norm_mat[sigs, , drop = FALSE])
}

contrib_root <- file.path(root, "results", "contributions")
norm_files <- list.files(contrib_root, pattern = "_contribution_normalized\\.tsv$",
                         recursive = TRUE, full.names = TRUE)

bio_rows <- list()
mmr_comp_rows <- list()
sample_rows <- list()

for (nf in norm_files) {
  parts <- strsplit(nf, .Platform$file.sep, fixed = TRUE)[[1]]
  fname <- parts[length(parts)]
  cosmic_key <- parts[length(parts) - 1]
  tool <- parts[length(parts) - 2]
  dataset <- sub("_contribution_normalized\\.tsv$", "", fname)

  mat <- as.matrix(read.table(nf, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE))
  cohort <- dataset_cohort(dataset, cfg)
  data_type <- dataset_data_type(dataset, cfg)
  ph <- pheno[pheno$Dataset == cohort, , drop = FALSE]
  if (!nrow(ph)) next

  mmr_sigs <- resolve_signature_group("MMR", cosmic_key, rownames(mat), groups)
  pole_sigs <- resolve_signature_group("POLE", cosmic_key, rownames(mat), groups)

  mmr_score <- sum_group(mat, mmr_sigs)
  pole_score <- sum_group(mat, pole_sigs)

  ph <- ph[ph$Sample %in% colnames(mat), , drop = FALSE]
  if (!nrow(ph)) next

  msi <- ph$Sample[ph$MSI_status == "MSI"]
  mss <- ph$Sample[ph$MSI_status %in% c("MSS", "MSI-L")]
  if (!length(mss)) mss <- ph$Sample[ph$MSI_status != "MSI"]
  pole_mut <- ph$Sample[ph$POLE_status == "MUT" & ph$Sample %in% mss]
  pole_wt <- setdiff(mss, pole_mut)

  for (s in ph$Sample) {
    ms <- if (s %in% msi) "MSI" else if (s %in% mss) "MSS" else NA_character_
    ps <- if (s %in% pole_mut) "POLE MUT" else if (s %in% pole_wt) "POLE WT" else NA_character_
    sample_rows[[length(sample_rows) + 1L]] <- data.frame(
      Tool = tool, COSMIC_Version = cosmic_key, Dataset = dataset,
      Cohort = cohort, Data_Type = if (is.na(data_type)) NA_character_ else data_type,
      Sample = s, MS_Status = ms, POLE_Status = ps,
      MMR_Contribution = unname(mmr_score[s]),
      POLE_Contribution = unname(pole_score[s]),
      stringsAsFactors = FALSE
    )
  }

  bio_rows[[length(bio_rows) + 1]] <- data.frame(
    Tool = tool,
    COSMIC_Version = cosmic_key,
    Dataset = dataset,
    Cohort = cohort,
    Data_Type = if (is.na(data_type)) NA_character_ else data_type,
    n_MSI = length(msi),
    n_MSS = length(mss),
    n_POLE_MUT = length(pole_mut),
    median_MMR_MSI = if (length(msi)) median(mmr_score[msi], na.rm = TRUE) else NA_real_,
    median_MMR_MSS = if (length(mss)) median(mmr_score[mss], na.rm = TRUE) else NA_real_,
    delta_MMR = if (length(msi) && length(mss)) median(mmr_score[msi], na.rm = TRUE) - median(mmr_score[mss], na.rm = TRUE) else NA_real_,
    p_MMR = if (length(msi) >= 2 && length(mss) >= 2) suppressWarnings(wilcox.test(mmr_score[msi], mmr_score[mss])$p.value) else NA_real_,
    median_POLE_MUT = if (length(pole_mut)) median(pole_score[pole_mut], na.rm = TRUE) else NA_real_,
    median_POLE_WT = if (length(pole_wt)) median(pole_score[pole_wt], na.rm = TRUE) else NA_real_,
    delta_POLE = if (length(pole_mut) && length(pole_wt)) median(pole_score[pole_mut], na.rm = TRUE) - median(pole_score[pole_wt], na.rm = TRUE) else NA_real_,
    p_POLE = if (length(pole_mut) >= 2 && length(pole_wt) >= 2) suppressWarnings(wilcox.test(pole_score[pole_mut], pole_score[pole_wt])$p.value) else NA_real_,
    stringsAsFactors = FALSE
  )

  if (length(msi) && length(mmr_sigs)) {
    sub <- mat[mmr_sigs, msi, drop = FALSE]
    means <- rowMeans(sub)
    if (sum(means) > 0) means <- means / sum(means)
    for (sg in names(means)) {
      mmr_comp_rows[[length(mmr_comp_rows) + 1]] <- data.frame(
        Tool = tool, COSMIC_Version = cosmic_key, Dataset = dataset,
        Cohort = cohort, Data_Type = if (is.na(data_type)) NA_character_ else data_type,
        Signature = sg, Mean_Contribution = as.numeric(means[sg]),
        stringsAsFactors = FALSE
      )
    }
  }
}

sample_df <- if (length(sample_rows)) do.call(rbind, sample_rows) else NULL
if (!is.null(sample_df)) {
  write.table(sample_df, file.path(out_dir, "tables", "sample_MMR_POLE_scores.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

if (length(bio_rows)) {
  bio_df <- do.call(rbind, bio_rows)
  write.table(bio_df, file.path(out_dir, "tables", "delta_MMR_POLE.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  p2 <- ggplot(bio_df, aes(x = Tool, y = delta_MMR, fill = COSMIC_Version)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    facet_wrap(~Dataset) +
    paper_theme(11, rotate_x = TRUE) +
    labs(title = "ΔMMR (median MSI − median MSS)", y = "ΔMMR", x = NULL)
  ggsave(file.path(out_dir, "plots", "delta_MMR_by_tool_cosmic.png"), p2, width = 12, height = 5, dpi = 150)

  p3 <- ggplot(bio_df, aes(x = Tool, y = delta_POLE, fill = COSMIC_Version)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    facet_wrap(~Dataset) +
    paper_theme(11, rotate_x = TRUE) +
    labs(title = "ΔPOLE (median MUT − median WT MSS)", y = "ΔPOLE", x = NULL)
  ggsave(file.path(out_dir, "plots", "delta_POLE_by_tool_cosmic.png"), p3, width = 12, height = 5, dpi = 150)
}

if (length(mmr_comp_rows)) {
  comp_df <- do.call(rbind, mmr_comp_rows)
  write.table(comp_df, file.path(out_dir, "tables", "MMR_composition_MSI.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  focus_tool <- if ("MutationalPatterns" %in% unique(comp_df$Tool)) "MutationalPatterns" else unique(comp_df$Tool)[1]
  subc <- comp_df[comp_df$Tool == focus_tool, , drop = FALSE]
  p4 <- ggplot(subc, aes(x = COSMIC_Version, y = Mean_Contribution, fill = Signature)) +
    geom_col(position = "stack", width = 0.7) +
    facet_wrap(~Dataset) +
    paper_theme(11) +
    labs(title = paste0("MMR signature composition in MSI (", focus_tool, ")"),
         y = "Relative contribution within MMR set", x = "COSMIC version")
  ggsave(file.path(out_dir, "plots", "MMR_composition_MSI.png"), p4, width = 10, height = 5, dpi = 150)
}

# ---------- Paper-style focus panels (Algorithm: Preclinical WES · v2) ----------
if (!is.null(sample_df)) {
  alg <- sample_df[sample_df$Dataset == "Preclinical" & sample_df$COSMIC_Version == "v2" &
                     (is.na(sample_df$Data_Type) | sample_df$Data_Type == "WES"), , drop = FALSE]
  cos_alg <- cos_df[cos_df$Dataset == "Preclinical" & cos_df$COSMIC_Version == "v2" &
                      (is.na(cos_df$Data_Type) | cos_df$Data_Type == "WES"), , drop = FALSE]

  if (nrow(cos_alg)) {
    p_cos <- ggplot(cos_alg, aes(x = Tool, y = Cosine_Similarity, fill = Tool)) +
      geom_boxplot(outlier.size = 0.5) +
      geom_hline(yintercept = cfg$cosine_threshold, linetype = "dashed") +
      scale_y_continuous(breaks = seq(0.75, 1, 0.05), limits = c(0.75, 1), expand = c(0, 0)) +
      paper_theme(10, rotate_x = TRUE) +
      theme(legend.position = "none") +
      labs(title = "Cosine similarity by tool (Preclinical WES · COSMIC v2)",
           y = "Cosine Similarity", x = NULL)
    save_paper_plot(p_cos, file.path(paper_dir, "alg_cosine_by_tool_v2"), width_cm = 14, height_cm = 8)
  }

  if (nrow(alg) && any(!is.na(alg$MS_Status))) {
    alg_ms <- alg[!is.na(alg$MS_Status), ]
    alg_ms$MS_Status <- factor(alg_ms$MS_Status, levels = c("MSS", "MSI"))
    p_mmr <- ggplot(alg_ms, aes(x = Tool, y = MMR_Contribution, fill = MS_Status)) +
      geom_boxplot(outlier.size = 0.4) +
      scale_fill_manual(values = paper_cols$msi_mss) +
      scale_y_unit() +
      paper_theme(10, rotate_x = TRUE) +
      labs(title = "MMRd signature contribution by tool (Preclinical · v2)",
           y = "Overall MMR Signatures Contribution", x = NULL, fill = "Microsatellite Status")
    cmp <- maybe_compare_means(method = "wilcox.test", label = "p.signif")
    if (!is.null(cmp)) p_mmr <- p_mmr + cmp
    save_paper_plot(p_mmr, file.path(paper_dir, "alg_MMR_boxplot_by_tool_v2"), width_cm = 14, height_cm = 8)

    # ΔMMR violin (paper construction)
    viol_rows <- list()
    for (tool in unique(alg_ms$Tool)) {
      msi_v <- alg_ms$MMR_Contribution[alg_ms$Tool == tool & alg_ms$MS_Status == "MSI"]
      mss_v <- alg_ms$MMR_Contribution[alg_ms$Tool == tool & alg_ms$MS_Status == "MSS"]
      if (!length(msi_v) || !length(mss_v)) next
      d <- paper_delta_distribution(msi_v, mss_v)
      viol_rows[[length(viol_rows) + 1L]] <- data.frame(Tool = tool, Delta = d, stringsAsFactors = FALSE)
    }
    if (length(viol_rows)) {
      viol <- do.call(rbind, viol_rows)
      p_v <- ggplot(viol, aes(x = Tool, y = Delta)) +
        geom_violin(fill = paper_cols$delta_mmr, colour = "brown3", adjust = 4) +
        geom_boxplot(width = 0.1, colour = "black", fill = "white", outlier.colour = NA) +
        scale_y_unit() +
        paper_theme(10, rotate_x = TRUE) +
        theme(legend.position = "none") +
        labs(title = "ΔMMR distribution by tool (Preclinical · v2)",
             y = "ΔMMR", x = NULL)
      save_paper_plot(p_v, file.path(paper_dir, "alg_delta_MMR_violin_v2"), width_cm = 14, height_cm = 8)
    }
  }

  if (nrow(alg) && any(!is.na(alg$POLE_Status))) {
    alg_p <- alg[!is.na(alg$POLE_Status), ]
    alg_p$POLE_Status <- factor(alg_p$POLE_Status, levels = c("POLE WT", "POLE MUT"))
    p_pole <- ggplot(alg_p, aes(x = Tool, y = POLE_Contribution, fill = POLE_Status)) +
      geom_boxplot(outlier.size = 0.4) +
      scale_fill_manual(values = paper_cols$pole) +
      scale_y_unit() +
      paper_theme(10, rotate_x = TRUE) +
      labs(title = "POLE signature contribution by tool (Preclinical · v2)",
           y = "Overall POLE Signatures Contribution", x = NULL, fill = "POLE status")
    cmp <- maybe_compare_means(method = "wilcox.test", label = "p.signif")
    if (!is.null(cmp)) p_pole <- p_pole + cmp
    save_paper_plot(p_pole, file.path(paper_dir, "alg_POLE_boxplot_by_tool_v2"), width_cm = 14, height_cm = 8)

    viol_rows <- list()
    for (tool in unique(alg_p$Tool)) {
      mut_v <- alg_p$POLE_Contribution[alg_p$Tool == tool & alg_p$POLE_Status == "POLE MUT"]
      wt_v <- alg_p$POLE_Contribution[alg_p$Tool == tool & alg_p$POLE_Status == "POLE WT"]
      if (!length(mut_v) || !length(wt_v)) next
      d <- paper_delta_distribution(mut_v, wt_v)
      viol_rows[[length(viol_rows) + 1L]] <- data.frame(Tool = tool, Delta = d, stringsAsFactors = FALSE)
    }
    if (length(viol_rows)) {
      viol <- do.call(rbind, viol_rows)
      p_v <- ggplot(viol, aes(x = Tool, y = Delta)) +
        geom_violin(fill = paper_cols$delta_pole, colour = "Dark Blue", adjust = 8) +
        geom_boxplot(width = 0.1, colour = "black", fill = "white", outlier.colour = NA) +
        scale_y_unit() +
        paper_theme(10, rotate_x = TRUE) +
        theme(legend.position = "none") +
        labs(title = "ΔPOLE distribution by tool (Preclinical · v2)",
             y = "ΔPOLE", x = NULL)
      save_paper_plot(p_v, file.path(paper_dir, "alg_delta_POLE_violin_v2"), width_cm = 14, height_cm = 8)
    }
  }

  # Reference focus: MutationalPatterns · Preclinical WES · all COSMIC
  if (length(mmr_comp_rows)) {
    ref_comp <- comp_df[comp_df$Tool == "MutationalPatterns" & comp_df$Dataset == "Preclinical", ]
    if (nrow(ref_comp)) {
      fill_vals <- paper_cols$mmr_sbs
      miss <- setdiff(unique(ref_comp$Signature), names(fill_vals))
      if (length(miss)) {
        extra <- setNames(grDevices::hcl.colors(length(miss), "Dark 3"), miss)
        fill_vals <- c(fill_vals, extra)
      }
      p_comp <- ggplot(ref_comp, aes(x = COSMIC_Version, y = Mean_Contribution, fill = Signature)) +
        geom_col(stat = "identity", width = 0.7) +
        scale_fill_manual(values = fill_vals) +
        scale_y_unit() +
        paper_theme(10) +
        labs(title = "MMR signature mix in MSI (MutationalPatterns · Preclinical)",
             y = "MMR Mutational Signatures Contribution", x = "Reference",
             fill = "MMR Mutational Signatures") +
        guides(fill = guide_legend(ncol = 2, byrow = TRUE, title.position = "top"))
      save_paper_plot(p_comp, file.path(paper_dir, "ref_MMR_composition_MP"), width_cm = 12, height_cm = 9)
    }
  }

  cos_ref <- cos_df[cos_df$Tool == "MutationalPatterns" & cos_df$Dataset %in% c("Preclinical", "Clinical") &
                      (is.na(cos_df$Data_Type) | cos_df$Data_Type == "WES"), ]
  if (nrow(cos_ref)) {
    p_ref <- ggplot(cos_ref, aes(x = Dataset, y = Cosine_Similarity, colour = COSMIC_Version)) +
      geom_boxplot(linewidth = 0.9) +
      scale_colour_manual(values = paper_cols$reference) +
      geom_hline(yintercept = cfg$cosine_threshold, linetype = "dashed") +
      coord_flip() +
      scale_y_continuous(breaks = seq(0.825, 1, 0.05), limits = c(0.825, 1), expand = c(0, 0)) +
      paper_theme(10) +
      labs(title = "Cosine similarity by reference (MutationalPatterns · WES)",
           y = "Cosine Similarity", x = NULL, colour = "Reference Dataset")
    save_paper_plot(p_ref, file.path(paper_dir, "ref_cosine_MP_by_cosmic"), width_cm = 12, height_cm = 7)
  }
}

message("Downstream outputs written to ", out_dir)
message("Paper-style panels: ", paper_dir)

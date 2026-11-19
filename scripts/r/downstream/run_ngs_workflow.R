#!/usr/bin/env Rscript
## Phase 2 — NGS workflow comparison in Battuello et al. Fig. 3 grammar.

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

out_dir <- file.path(root, "results", "downstream", "ngs_workflow")
paper_dir <- file.path(out_dir, "plots", "paper_style")
dir.create(file.path(out_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(paper_dir, recursive = TRUE, showWarnings = FALSE)

cos_path <- file.path(root, "results", "cosine", "CosineResults_AllTools_AllCOSMIC.tsv")
if (!file.exists(cos_path)) {
  system2("Rscript", file.path(dirname(script_dir), "00_combine_cosine.R"))
}
cos_df <- read.table(cos_path, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)

ngs_types <- c("WGS", "WES", "TSO-500")  # paper order for horizontal cosine
pre_cos <- cos_df[cos_df$Data_Type %in% ngs_types & grepl("^Preclinical", cos_df$Dataset), , drop = FALSE]
if (!nrow(pre_cos)) {
  stop("No preclinical NGS cosine rows found. Run fitting for Preclinical_WGS / Preclinical_TSO-500 first.")
}
pre_cos$Assigned <- pre_cos$Cosine_Similarity >= cfg$cosine_threshold
pre_cos$Data_Type <- factor(pre_cos$Data_Type, levels = ngs_types)

# ---------- Supplement grids ----------
p1 <- ggplot(pre_cos, aes(x = Data_Type, y = Cosine_Similarity, fill = Data_Type)) +
  geom_boxplot(outlier.size = 0.4) +
  geom_hline(yintercept = cfg$cosine_threshold, linetype = "dashed") +
  facet_grid(Tool ~ COSMIC_Version) +
  paper_theme(10, rotate_x = TRUE) +
  theme(legend.position = "none") +
  labs(title = "Cosine similarity by NGS workflow (preclinical)",
       y = "Cosine similarity", x = NULL)
ggsave(file.path(out_dir, "plots", "ngs_cosine_by_datatype.png"), p1, width = 14, height = 10, dpi = 150)

asn <- aggregate(Assigned ~ Tool + COSMIC_Version + Data_Type, data = pre_cos,
                 FUN = function(x) c(n = length(x), assigned = sum(x, na.rm = TRUE), frac = mean(x, na.rm = TRUE)))
asn_df <- do.call(data.frame, asn)
names(asn_df) <- c("Tool", "COSMIC_Version", "Data_Type", "n", "n_assigned", "frac_assigned")
write.table(asn_df, file.path(out_dir, "tables", "ngs_assigned_fraction.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

wilcox_rows <- list()
for (tool in unique(pre_cos$Tool)) {
  for (cv in unique(pre_cos$COSMIC_Version)) {
    sub <- pre_cos[pre_cos$Tool == tool & pre_cos$COSMIC_Version == cv, ]
    dts <- levels(droplevels(sub$Data_Type))
    if (length(dts) < 2) next
    for (i in seq_len(length(dts) - 1L)) {
      for (j in (i + 1L):length(dts)) {
        a <- sub$Cosine_Similarity[sub$Data_Type == dts[i]]
        b <- sub$Cosine_Similarity[sub$Data_Type == dts[j]]
        if (sum(is.finite(a)) < 2 || sum(is.finite(b)) < 2) next
        wt <- suppressWarnings(wilcox.test(a, b))
        wilcox_rows[[length(wilcox_rows) + 1L]] <- data.frame(
          Tool = tool, COSMIC_Version = cv,
          Data_Type_A = dts[i], Data_Type_B = dts[j],
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
              file.path(out_dir, "tables", "ngs_cosine_wilcoxon.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

# ---------- Sample-level contributions ----------
pheno <- read.table(file.path(root, "metadata", "sample_phenotypes.tsv"),
                    header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
pheno_pre <- pheno[pheno$Dataset == "Preclinical", ]

sum_group <- function(norm_mat, sigs) {
  sigs <- intersect(sigs, rownames(norm_mat))
  if (!length(sigs)) return(setNames(rep(0, ncol(norm_mat)), colnames(norm_mat)))
  if (length(sigs) == 1) setNames(as.numeric(norm_mat[sigs, ]), colnames(norm_mat)) else colSums(norm_mat[sigs, , drop = FALSE])
}

dt_to_dataset <- list()
for (dn in names(cfg$datasets)) {
  if (!identical(cfg$datasets[[dn]]$cohort, "Preclinical")) next
  dt_to_dataset[[cfg$datasets[[dn]]$data_type]] <- dn
}

bio_rows <- list()
flat_rows <- list()
sample_rows <- list()
class_rows <- list()

for (tool_dir in list.dirs(file.path(root, "results", "contributions"), recursive = FALSE, full.names = TRUE)) {
  tool <- basename(tool_dir)
  for (cv_dir in list.dirs(tool_dir, recursive = FALSE, full.names = TRUE)) {
    cosmic_key <- basename(cv_dir)
    for (dt in ngs_types) {
      ds_name <- dt_to_dataset[[dt]]
      if (is.null(ds_name)) next
      nf <- file.path(cv_dir, paste0(ds_name, "_contribution_normalized.tsv"))
      if (!file.exists(nf) && identical(dt, "WES")) {
        nf <- file.path(cv_dir, "Preclinical_contribution_normalized.tsv")
        ds_name <- "Preclinical"
      }
      if (!file.exists(nf)) next

      mat <- as.matrix(read.table(nf, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE))
      ph <- pheno_pre[pheno_pre$Sample %in% colnames(mat), , drop = FALSE]
      if (!nrow(ph)) next

      mmr_sigs <- as.character(resolve_signature_group("MMR", cosmic_key, rownames(mat), groups))
      pole_sigs <- as.character(resolve_signature_group("POLE", cosmic_key, rownames(mat), groups))
      flat_sigs <- as.character(resolve_signature_group("flat", cosmic_key, rownames(mat), groups))
      art_sigs <- as.character(resolve_signature_group("artefact", cosmic_key, rownames(mat), groups))

      mmr_score <- sum_group(mat, mmr_sigs)
      pole_score <- sum_group(mat, pole_sigs)
      flat_score <- sum_group(mat, flat_sigs)
      art_score <- sum_group(mat, art_sigs)

      msi <- ph$Sample[ph$MSI_status == "MSI"]
      mss <- ph$Sample[ph$MSI_status %in% c("MSS", "MSI-L")]
      if (!length(mss)) mss <- ph$Sample[ph$MSI_status != "MSI"]
      pole_mut <- ph$Sample[ph$POLE_status == "MUT" & ph$Sample %in% mss]
      pole_wt <- setdiff(mss, pole_mut)

      for (s in ph$Sample) {
        ms <- if (s %in% msi) "MSI" else if (s %in% mss) "MSS" else NA_character_
        ps <- if (s %in% pole_mut) "POLE MUT" else if (s %in% pole_wt) "POLE WT" else NA_character_
        sample_rows[[length(sample_rows) + 1L]] <- data.frame(
          Tool = tool, COSMIC_Version = cosmic_key, Data_Type = dt, Sample = s,
          MS_Status = ms, POLE_Status = ps,
          MMR_Contribution = unname(mmr_score[s]),
          POLE_Contribution = unname(pole_score[s]),
          Flat_Contribution = unname(flat_score[s]),
          stringsAsFactors = FALSE
        )
      }

      bio_rows[[length(bio_rows) + 1L]] <- data.frame(
        Tool = tool, COSMIC_Version = cosmic_key, Dataset = ds_name, Data_Type = dt,
        n_samples = ncol(mat), n_MSI = length(msi), n_MSS = length(mss), n_POLE_MUT = length(pole_mut),
        delta_MMR = if (length(msi) && length(mss)) median(mmr_score[msi], na.rm = TRUE) - median(mmr_score[mss], na.rm = TRUE) else NA_real_,
        p_MMR = if (length(msi) >= 2 && length(mss) >= 2) suppressWarnings(wilcox.test(mmr_score[msi], mmr_score[mss])$p.value) else NA_real_,
        delta_POLE = if (length(pole_mut) && length(pole_wt)) median(pole_score[pole_mut], na.rm = TRUE) - median(pole_score[pole_wt], na.rm = TRUE) else NA_real_,
        p_POLE = if (length(pole_mut) >= 2 && length(pole_wt) >= 2) suppressWarnings(wilcox.test(pole_score[pole_mut], pole_score[pole_wt])$p.value) else NA_real_,
        median_flat_MSI = if (length(msi)) median(flat_score[msi], na.rm = TRUE) else NA_real_,
        median_flat_all = median(flat_score, na.rm = TRUE),
        stringsAsFactors = FALSE
      )

      if (length(msi)) {
        flat_rows[[length(flat_rows) + 1L]] <- data.frame(
          Tool = tool, COSMIC_Version = cosmic_key, Data_Type = dt,
          median_flat_MSI = median(flat_score[msi], na.rm = TRUE),
          mean_flat_MSI = mean(flat_score[msi], na.rm = TRUE),
          stringsAsFactors = FALSE
        )

        # Paper Fig. 3D: median per-signature then renormalize class sums among MSI
        med_sig <- apply(mat[, msi, drop = FALSE], 1, median)
        if (sum(med_sig) > 0) med_sig <- med_sig / sum(med_sig)
        mmr_c <- sum(med_sig[intersect(names(med_sig), mmr_sigs)])
        art_c <- sum(med_sig[intersect(names(med_sig), art_sigs)])
        flat_c <- sum(med_sig[intersect(names(med_sig), flat_sigs)])
        other_c <- max(0, 1 - (mmr_c + art_c + flat_c))
        for (pair in list(
          c(mmr_c, "MMR signatures"),
          c(art_c, "Artefacts signatures"),
          c(flat_c, "Flat signatures"),
          c(other_c, "Other signatures")
        )) {
          class_rows[[length(class_rows) + 1L]] <- data.frame(
            Tool = tool, COSMIC_Version = cosmic_key, Data_Type = dt,
            Signatures_type = pair[[2]], Signature_Contributions = as.numeric(pair[[1]]),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
}

sample_df <- if (length(sample_rows)) do.call(rbind, sample_rows) else NULL
if (!is.null(sample_df)) {
  write.table(sample_df, file.path(out_dir, "tables", "ngs_sample_MMR_POLE_scores.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

if (length(bio_rows)) {
  bio_df <- do.call(rbind, bio_rows)
  bio_df$Data_Type <- factor(bio_df$Data_Type, levels = ngs_types)
  write.table(bio_df, file.path(out_dir, "tables", "ngs_delta_MMR_POLE.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  p2 <- ggplot(bio_df, aes(x = Data_Type, y = delta_MMR, fill = Data_Type)) +
    geom_col(width = 0.7) +
    facet_grid(Tool ~ COSMIC_Version) +
    paper_theme(10, rotate_x = TRUE) +
    theme(legend.position = "none") +
    labs(title = "ΔMMR by NGS workflow (preclinical)", y = "ΔMMR (MSI − MSS)", x = NULL)
  ggsave(file.path(out_dir, "plots", "ngs_delta_MMR.png"), p2, width = 14, height = 10, dpi = 150)

  p3 <- ggplot(bio_df, aes(x = Data_Type, y = delta_POLE, fill = Data_Type)) +
    geom_col(width = 0.7) +
    facet_grid(Tool ~ COSMIC_Version) +
    paper_theme(10, rotate_x = TRUE) +
    theme(legend.position = "none") +
    labs(title = "ΔPOLE by NGS workflow (preclinical)", y = "ΔPOLE (MUT − WT)", x = NULL)
  ggsave(file.path(out_dir, "plots", "ngs_delta_POLE.png"), p3, width = 14, height = 10, dpi = 150)
}

if (length(flat_rows)) {
  flat_df <- do.call(rbind, flat_rows)
  flat_df$Data_Type <- factor(flat_df$Data_Type, levels = ngs_types)
  write.table(flat_df, file.path(out_dir, "tables", "ngs_flat_signatures_MSI.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

if (length(class_rows)) {
  class_df <- do.call(rbind, class_rows)
  write.table(class_df, file.path(out_dir, "tables", "ngs_signature_classes_MSI.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

# ---------- Paper Fig. 3-style focus: MutationalPatterns · v3.2 ----------
focus_tool <- "MutationalPatterns"
focus_cv <- "v3.2"

mp_cos <- pre_cos[pre_cos$Tool == focus_tool & pre_cos$COSMIC_Version == focus_cv, ]
if (nrow(mp_cos)) {
  # Fig 3A — horizontal grey boxplot
  p_a <- ggplot(mp_cos, aes(x = Data_Type, y = Cosine_Similarity)) +
    geom_boxplot(colour = "black", fill = paper_cols$cosine_fill, outlier.size = 0.6) +
    coord_flip() +
    scale_y_continuous(breaks = seq(0.9, 1, 0.025), limits = c(0.89, 1), expand = c(0, 0)) +
    geom_hline(yintercept = cfg$cosine_threshold, linetype = "dashed") +
    paper_theme(10) +
    labs(title = "Cosine similarity by NGS workflow (MP · v3.2)",
         x = "NGS Workflow", y = "Cosine Similarity")
  save_paper_plot(p_a, file.path(paper_dir, "ngs_cosine_MP_v3.2"), width_cm = 10, height_cm = 7)
  # also keep legacy focus path
  ggsave(file.path(out_dir, "plots", "ngs_cosine_MP_v3.2.png"), p_a, width = 5, height = 3.5, dpi = 150)
}

if (!is.null(sample_df)) {
  mp_s <- sample_df[sample_df$Tool == focus_tool & sample_df$COSMIC_Version == focus_cv, ]
  mp_s$Data_Type <- factor(mp_s$Data_Type, levels = ngs_types)

  if (nrow(mp_s) && any(!is.na(mp_s$MS_Status))) {
    ms <- mp_s[!is.na(mp_s$MS_Status), ]
    ms$MS_Status <- factor(ms$MS_Status, levels = c("MSS", "MSI"))
    p_b <- ggplot(ms, aes(x = Data_Type, y = MMR_Contribution, fill = MS_Status)) +
      geom_boxplot(outlier.size = 0.5) +
      scale_fill_manual(values = paper_cols$msi_mss) +
      scale_y_unit() +
      paper_theme(10) +
      labs(title = "MMR signatures by NGS workflow (MP · v3.2)",
           x = "NGS Workflow", y = "Overall MMR Signatures Contribution",
           fill = "Microsatellite Status")
    cmp <- maybe_compare_means(method = "wilcox.test", label = "p.signif")
    if (!is.null(cmp)) p_b <- p_b + cmp
    save_paper_plot(p_b, file.path(paper_dir, "ngs_MMR_boxplot_MP_v3.2"), width_cm = 10, height_cm = 7.5)
    ggsave(file.path(out_dir, "plots", "ngs_delta_MMR_MP_v3.2.png"), p_b, width = 5, height = 4, dpi = 150)
  }

  if (nrow(mp_s) && any(!is.na(mp_s$POLE_Status))) {
    ps <- mp_s[!is.na(mp_s$POLE_Status), ]
    ps$POLE_Status <- factor(ps$POLE_Status, levels = c("POLE WT", "POLE MUT"))
    p_c <- ggplot(ps, aes(x = Data_Type, y = POLE_Contribution, fill = POLE_Status)) +
      geom_boxplot(outlier.size = 0.5) +
      scale_fill_manual(values = paper_cols$pole) +
      scale_y_unit() +
      paper_theme(10) +
      labs(title = "POLE signatures by NGS workflow (MP · v3.2)",
           x = "NGS Workflow", y = "Overall POLE Signatures Contribution",
           fill = "POLE status")
    cmp <- maybe_compare_means(method = "wilcox.test", label = "p.signif")
    if (!is.null(cmp)) p_c <- p_c + cmp
    save_paper_plot(p_c, file.path(paper_dir, "ngs_POLE_boxplot_MP_v3.2"), width_cm = 10, height_cm = 7.5)
  }
}

if (length(class_rows)) {
  class_df <- do.call(rbind, class_rows)
  mp_cl <- class_df[class_df$Tool == focus_tool & class_df$COSMIC_Version == focus_cv, ]
  if (nrow(mp_cl)) {
    mp_cl$Signatures_type <- factor(mp_cl$Signatures_type,
                                    levels = c("Other signatures", "Artefacts signatures",
                                               "Flat signatures", "MMR signatures"))
    mp_cl$Data_Type <- factor(mp_cl$Data_Type, levels = ngs_types)
    p_d <- ggplot(mp_cl, aes(x = Data_Type, y = Signature_Contributions, fill = Signatures_type)) +
      geom_col(stat = "identity", width = 0.7) +
      scale_fill_manual(values = paper_cols$sig_class) +
      scale_y_unit() +
      paper_theme(10) +
      labs(title = "Signature classes in MSI (MP · v3.2)",
           x = "NGS Workflow", y = "Mutational Signatures Contribution",
           fill = "Mutational Signatures Types") +
      guides(fill = guide_legend(ncol = 2, nrow = 3, byrow = TRUE, title.position = "top"))
    save_paper_plot(p_d, file.path(paper_dir, "ngs_signature_classes_MP_v3.2"), width_cm = 10, height_cm = 7.5)
    ggsave(file.path(out_dir, "plots", "ngs_flat_MSI_MP_v3.2.png"), p_d, width = 5, height = 4, dpi = 150)
  }
}

message("NGS workflow downstream written to ", out_dir)
message("Paper-style panels: ", paper_dir)

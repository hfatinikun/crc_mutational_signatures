## Shared Battuello et al. (bbae249) ggplot grammar for thesis figures.

paper_cols <- list(
  msi_mss = c(MSS = "lightblue1", MSI = "salmon2"),
  pole = c(`POLE WT` = "olivedrab4", `POLE MUT` = "darkgoldenrod1",
           POLE_WT = "olivedrab4", POLE_MUT = "darkgoldenrod1"),
  sig_class = c(
    `Other signatures` = "gray80",
    `Artefacts signatures` = "deepskyblue3",
    `Flat signatures` = "goldenrod1",
    `MMR signatures` = "firebrick3"
  ),
  mmr_sbs = c(
    Other_SBS = "gray80",
    SBS6 = "seagreen1", SBS14 = "royalblue3", SBS15 = "orchid3",
    SBS20 = "turquoise2", SBS21 = "goldenrod1", SBS26 = "springgreen4",
    SBS44 = "orangered2",
    Signature_6 = "seagreen1", Signature_14 = "royalblue3",
    Signature_15 = "orchid3", Signature_20 = "turquoise2",
    Signature_21 = "goldenrod1", Signature_26 = "springgreen4"
  ),
  reference = c(v2 = "#1B8148", v3.2 = "Orange", v3.4 = "#3852BD",
                v3.5 = "#8B008B", v3.6 = "#A0522D"),
  delta_mmr = "coral",
  delta_pole = "aliceblue",
  cosine_fill = "grey"
)

paper_theme <- function(base_size = 10, rotate_x = FALSE) {
  th <- ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = base_size),
      axis.text = ggplot2::element_text(size = base_size),
      legend.text = ggplot2::element_text(size = base_size),
      legend.title = ggplot2::element_text(size = base_size),
      legend.position = "top",
      panel.border = ggplot2::element_rect(colour = "black", linewidth = 0.4),
      plot.title = ggplot2::element_text(size = base_size + 1, face = "bold", hjust = 0.5)
    )
  if (isTRUE(rotate_x)) {
    th <- th + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = -45, vjust = 0.5, hjust = 0.1))
  }
  th
}

save_paper_plot <- function(plot, path_no_ext, width_cm = 10, height_cm = 7.5, dpi = 150) {
  dir.create(dirname(path_no_ext), recursive = TRUE, showWarnings = FALSE)
  w_in <- width_cm / 2.54
  h_in <- height_cm / 2.54
  ggplot2::ggsave(paste0(path_no_ext, ".png"), plot, width = w_in, height = h_in, dpi = dpi)
  tryCatch(
    ggplot2::ggsave(paste0(path_no_ext, ".svg"), plot, width = w_in, height = h_in),
    error = function(e) message("SVG save skipped: ", conditionMessage(e))
  )
  invisible(path_no_ext)
}

## Paper Δ distribution used for violin panels (Fig. 4D/F style).
paper_delta_distribution <- function(group_high, group_low) {
  med_high <- stats::median(group_high, na.rm = TRUE)
  med_low <- stats::median(group_low, na.rm = TRUE)
  dist <- c(group_high - med_low, abs(group_low - med_high))
  abs(dist)
}

scale_y_unit <- function(limits = c(0, 1), breaks = seq(0, 1, 0.25)) {
  ggplot2::scale_y_continuous(breaks = breaks, limits = limits, expand = c(0, 0))
}

maybe_compare_means <- function(mapping = NULL, ...) {
  if (requireNamespace("ggpubr", quietly = TRUE)) {
    return(ggpubr::stat_compare_means(mapping = mapping, ...))
  }
  NULL
}

maybe_force_panelsizes <- function(rows_cm = 6, cols_cm = 8) {
  if (requireNamespace("ggh4x", quietly = TRUE)) {
    return(ggh4x::force_panelsizes(rows = grid::unit(rows_cm, "cm"),
                                   cols = grid::unit(cols_cm, "cm")))
  }
  NULL
}

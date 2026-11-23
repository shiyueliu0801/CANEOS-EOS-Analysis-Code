# ==============================================================================
# Module 07A: Pure-R Figure 2 Generation (TIFF + Vector PDF + Editable Vector PPTX)
# ==============================================================================
# Consumes outputs/03_regimen_coverage/Figure2_Plot_Data.csv produced dynamically
# by Module 03. Renders:
#   - Panel A: Regimen family grouping, evaluable denominators (n/N), tri-state
#     isolate composition (Covered / Uncovered / Unknown), and design-weighted
#     coverage (95% CI).
#   - Panel B: Same-isolate paired coverage differences (pp, 95% CI, paired n),
#     full-cohort sharp logical bounds, Joint 95% lower confidence limit, and
#     Simultaneous (4-contrast family) 95% lower confidence limit.
#
# Optional PPTX export uses `officer` and `rvg` to emit native DrawingML
# vector shapes and text (never degrades to an embedded PNG raster).
# ==============================================================================

render_figure2 <- function(plot_data_csv, out_dir, export_pptx = FALSE) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dat <- readr::read_csv(plot_data_csv, show_col_types = FALSE)

  a <- dat[dat$Panel == "A", , drop = FALSE]
  a <- a[order(a$Order), , drop = FALSE]
  a$Family <- ifelse(grepl("GEN", a$Label, fixed = TRUE),
                     "WHO-recommended (Gentamicin)", "China-recommended (3GC)")
  a$YLabel <- sprintf("%s\n(n = %s/%s)",
                      a$Label,
                      format(a$Evaluable_n, big.mark = ",", trim = TRUE),
                      format(a$Candidate_N, big.mark = ",", trim = TRUE))
  a$YLabel_f <- factor(a$YLabel, levels = rev(a$YLabel))

  b <- dat[dat$Panel == "B", , drop = FALSE]
  b <- b[order(b$Order), , drop = FALSE]
  b$YLabel <- sprintf("%s\n(paired n = %s)",
                      b$Label,
                      format(b$Paired_n, big.mark = ",", trim = TRUE))
  b$YLabel_f <- factor(b$YLabel, levels = rev(b$YLabel))

  # --- Panel A Left: Tri-state composition (Covered / Uncovered / Unknown) ---
  a_long <- tidyr::pivot_longer(
    a[, c("YLabel_f", "Family", "Covered_share", "Uncovered_share", "Unknown_share")],
    cols = c("Covered_share", "Uncovered_share", "Unknown_share"),
    names_to = "Status", values_to = "Share"
  )
  a_long$Status <- factor(
    a_long$Status,
    levels = c("Unknown_share", "Uncovered_share", "Covered_share"),
    labels = c("Unknown", "Uncovered", "Covered")
  )
  a_long$Pct <- 100 * a_long$Share
  a_long$PctLabel <- ifelse(a_long$Pct >= 7, sprintf("%.0f%%", a_long$Pct),
                            ifelse(a_long$Pct >= 2.0, sprintf("%.1f%%", a_long$Pct), ""))

  pal_a <- c("Covered" = "#2B6CB0", "Uncovered" = "#CBD5E0", "Unknown" = "#ECC94B")

  p_a1 <- ggplot2::ggplot(a_long, ggplot2::aes(x = Pct, y = YLabel_f, fill = Status)) +
    ggplot2::geom_hline(yintercept = 4.5, linetype = "dotted", colour = "#4A5568", linewidth = 0.5) +
    ggplot2::geom_col(width = 0.64, colour = "white", linewidth = 0.3) +
    ggplot2::geom_text(
      ggplot2::aes(label = PctLabel),
      position = ggplot2::position_stack(vjust = 0.5),
      size = 2.8, colour = "#1A202C", fontface = "bold"
    ) +
    ggplot2::scale_fill_manual(values = pal_a, breaks = c("Covered", "Uncovered", "Unknown")) +
    ggplot2::scale_x_continuous(limits = c(0, 102), expand = c(0, 0),
                                breaks = seq(0, 100, 25), labels = function(x) paste0(x, "%")) +
    ggplot2::labs(
      title = "A. Isolate-level regimen status (N = 1,196) & weighted coverage",
      subtitle = "Top two rows: WHO gentamicin regimens; bottom four rows: China 3GC regimens",
      x = "Share of all isolates (%)", y = NULL, fill = NULL
    ) +
    ggplot2::theme_minimal(base_size = 9.5) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 10.5, hjust = 0),
      plot.subtitle = ggplot2::element_text(size = 8.2, colour = "#4A5568"),
      legend.position = "bottom",
      legend.margin = ggplot2::margin(t = 0, b = 0),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(face = "bold", size = 8.4, colour = "#1A202C", lineheight = 0.95)
    )

  # --- Panel A Right: Design-weighted coverage among evaluable isolates (95% CI) ---
  a$Cov_pct <- 100 * a$Weighted_coverage
  a$Low_pct <- 100 * a$Coverage_CI_low
  a$Hi_pct  <- 100 * a$Coverage_CI_high
  a$TextLab <- sprintf("%.1f%% (%.1f\u2013%.1f%%)", a$Cov_pct, a$Low_pct, a$Hi_pct)

  p_a2 <- ggplot2::ggplot(a, ggplot2::aes(x = Cov_pct, y = YLabel_f, colour = Family)) +
    ggplot2::geom_hline(yintercept = 4.5, linetype = "dotted", colour = "#4A5568", linewidth = 0.5) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = Low_pct, xmax = Hi_pct),
                            height = 0.22, linewidth = 0.75) +
    ggplot2::geom_point(size = 2.8) +
    ggplot2::geom_text(ggplot2::aes(x = Hi_pct + 1.8, label = TextLab),
                       hjust = 0, size = 2.8, colour = "#1A202C", show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = c(
      "WHO-recommended (Gentamicin)" = "#1B4F72",
      "China-recommended (3GC)"      = "#2E86AB"
    )) +
    ggplot2::scale_x_continuous(limits = c(55, 122), breaks = seq(60, 100, 10),
                                labels = function(x) paste0(x, "%")) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(
      title = " ",
      subtitle = "Design-weighted coverage among evaluable isolates (95% CI)",
      x = "Weighted coverage (%)", y = NULL, colour = NULL
    ) +
    ggplot2::theme_minimal(base_size = 9.5) +
    ggplot2::theme(
      plot.subtitle = ggplot2::element_text(size = 8.2, colour = "#4A5568"),
      legend.position = "bottom",
      legend.margin = ggplot2::margin(t = 0, b = 0),
      axis.text.y = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(t = 5, r = 14, b = 5, l = 4)
    )

  # --- Panel B: Paired differences, sharp logical bounds, Joint & Simultaneous 95% LCL ---
  fmt_signed <- function(x) sub("^-", "\u2212", sprintf("%.1f", x))
  b$DiffLab  <- sprintf("+%.1f pp (%.1f\u2013%.1f)",
                        b$Difference_pp, b$Difference_CI_low_pp, b$Difference_CI_high_pp)
  b$BoundLab <- sprintf("Bounds: %s\u2013%s | Joint LCL: %s | Sim. LCL: %s",
                        fmt_signed(b$Sample_lower_pp), fmt_signed(b$Sample_upper_pp),
                        fmt_signed(b$Joint_lower_pp),  fmt_signed(b$Simultaneous_lower_pp))

  lcl_pts <- rbind(
    data.frame(YLabel_f = b$YLabel_f, x = b$Joint_lower_pp,        Metric = "Joint 95% lower limit",        stringsAsFactors = FALSE),
    data.frame(YLabel_f = b$YLabel_f, x = b$Simultaneous_lower_pp, Metric = "Simultaneous 95% lower limit", stringsAsFactors = FALSE)
  )
  lcl_pts$Metric <- factor(lcl_pts$Metric, levels = c("Joint 95% lower limit", "Simultaneous 95% lower limit"))

  p_b <- ggplot2::ggplot(b, ggplot2::aes(y = YLabel_f)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "#718096", linewidth = 0.5) +
    # Full-cohort sharp logical identification interval
    ggplot2::geom_segment(
      ggplot2::aes(x = Sample_lower_pp, xend = Sample_upper_pp, yend = YLabel_f),
      linewidth = 5.2, colour = "#CBD5E0", alpha = 0.9
    ) +
    # Joint 95% lower limit (amber diamond) and Simultaneous 95% lower limit (red triangle)
    ggplot2::geom_point(
      data = lcl_pts,
      ggplot2::aes(x = x, y = YLabel_f, shape = Metric, colour = Metric),
      size = 2.9
    ) +
    # Paired complete-case 95% CI & point estimate
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = Difference_CI_low_pp, xmax = Difference_CI_high_pp),
      height = 0.22, linewidth = 0.85, colour = "#1B4F72"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = Difference_pp),
      size = 3.0, colour = "#1B4F72"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = 27.5, label = DiffLab),
      hjust = 0, vjust = -0.35, size = 2.85, colour = "#1A202C", fontface = "bold"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = 27.5, label = BoundLab),
      hjust = 0, vjust = 1.25, size = 2.55, colour = "#4A5568"
    ) +
    ggplot2::scale_shape_manual(values = c(
      "Joint 95% lower limit"        = 18,
      "Simultaneous 95% lower limit" = 17
    )) +
    ggplot2::scale_colour_manual(values = c(
      "Joint 95% lower limit"        = "#D97706",
      "Simultaneous 95% lower limit" = "#C53030"
    )) +
    ggplot2::scale_x_continuous(limits = c(-7, 68), breaks = seq(-5, 25, 5)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(
      title = "B. Same-isolate paired coverage differences (percentage points) & sharp logical bounds",
      subtitle = "Blue points/bars: paired difference (95% CI); grey band: full-cohort sharp logical bounds (N = 1,196)",
      x = "Coverage difference (percentage points, Gentamicin vs 3GC regimen)",
      y = NULL, shape = NULL, colour = NULL
    ) +
    ggplot2::theme_minimal(base_size = 9.5) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 10.5, hjust = 0),
      plot.subtitle = ggplot2::element_text(size = 8.2, colour = "#4A5568"),
      legend.position = "bottom",
      legend.margin = ggplot2::margin(t = 0, b = 0),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(face = "bold", size = 8.4, colour = "#1A202C", lineheight = 0.95),
      plot.margin = ggplot2::margin(t = 6, r = 14, b = 6, l = 5)
    )

  fig2 <- (p_a1 | p_a2) / p_b +
    patchwork::plot_layout(heights = c(1.22, 1.0), widths = c(1.1, 1.0)) &
    ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", colour = NA))

  tiff_path <- file.path(out_dir, "Figure2_300dpi.tiff")
  pdf_path  <- file.path(out_dir, "Figure2.pdf")
  pptx_path <- file.path(out_dir, "Figure2_Editable.pptx")

  ggplot2::ggsave(tiff_path, plot = fig2, width = 11.2, height = 7.8, units = "in",
                  dpi = 300, compression = "lzw", bg = "white", device = "tiff")
  ggplot2::ggsave(pdf_path,  plot = fig2, width = 11.2, height = 7.8, units = "in",
                  bg = "white", device = grDevices::cairo_pdf)

  if (isTRUE(export_pptx)) {
    tryCatch({
      if (!requireNamespace("officer", quietly = TRUE) || !requireNamespace("rvg", quietly = TRUE))
        stop("Optional PPTX requires officer and rvg.")
  # Export strictly native editable DrawingML vector PPTX via rvg::dml()
  ppt <- officer::read_pptx()
  ppt <- officer::add_slide(ppt, layout = "Blank", master = "Office Theme")
  ppt <- officer::ph_with(
    ppt,
    value = rvg::dml(ggobj = fig2, bg = "white"),
    location = officer::ph_location(left = 0.20, top = 0.25, width = 9.60, height = 6.70)
  )
  print(ppt, target = pptx_path)

    }, error = function(e) warning("Optional Figure 2 PPTX was not generated: ", conditionMessage(e), call. = FALSE))
  }
  if (!file.exists(tiff_path) || file.info(tiff_path)$size <= 0) stop("Required Figure 2 TIFF missing or empty.", call. = FALSE)
  message("Rendered mandatory Figure 2 TIFF; optional PPTX is checked separately: ", out_dir)
  invisible(list(tiff = tiff_path, pdf = pdf_path, pptx = if (file.exists(pptx_path) && !dir.exists(pptx_path) && file.info(pptx_path)$size > 0) pptx_path else NULL))
}

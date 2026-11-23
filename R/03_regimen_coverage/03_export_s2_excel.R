# ==============================================================================
# Module 03: Pure-R Export of S2 Fig Heatmap Workbook (S2_Heatmap.xlsx)
# ==============================================================================
# Replaces the legacy Python exporter with native R (openxlsx), reproducing all
# 23 x 23 coverage cells, evaluable counts, and the 100-step white-to-blue
# colour gradient (#FFFFFF -> #0066CC).
# ==============================================================================

export_s2_heatmap_excel <- function(rates_df, counts_df, drug_key_df, palette, out_xlsx) {
  wb <- openxlsx::createWorkbook()

  # Sheet 1: Coverage heatmap (23 x 23 formatted cells)
  openxlsx::addWorksheet(wb, "Coverage heatmap", gridLines = FALSE)
  openxlsx::freezePane(wb, "Coverage heatmap", firstActiveRow = 5, firstActiveCol = 2)

  title_txt <- "S2 Fig. In vitro coverage of individual antibiotics and two-drug combinations in culture-confirmed early-onset sepsis, 2015\u20132024."
  note_txt  <- "The diagonal and upper triangle are unweighted. The lower triangle uses macro-region\u2013year calibration weights. Each n is the evaluable isolate count."

  openxlsx::writeData(wb, "Coverage heatmap", title_txt, startRow = 1, startCol = 1, colNames = FALSE)
  openxlsx::mergeCells(wb, "Coverage heatmap", cols = 1:24, rows = 1)
  openxlsx::addStyle(wb, "Coverage heatmap",
                     openxlsx::createStyle(fontName = "Times New Roman", fontSize = 13, textDecoration = "bold"),
                     rows = 1, cols = 1)

  openxlsx::writeData(wb, "Coverage heatmap", note_txt, startRow = 2, startCol = 1, colNames = FALSE)
  openxlsx::mergeCells(wb, "Coverage heatmap", cols = 1:24, rows = 2)
  openxlsx::addStyle(wb, "Coverage heatmap",
                     openxlsx::createStyle(fontName = "Times New Roman", fontSize = 10, wrapText = TRUE),
                     rows = 2, cols = 1)
  openxlsx::setRowHeights(wb, "Coverage heatmap", rows = 2, heights = 28)

  drugs <- names(rates_df)[-1L]
  stopifnot(length(drugs) == 23L, nrow(rates_df) == 23L, nrow(counts_df) == 23L)
  rate_mat  <- as.matrix(rates_df[, -1L, drop = FALSE])
  count_mat <- as.matrix(counts_df[, -1L, drop = FALSE])
  storage.mode(rate_mat)  <- "numeric"
  storage.mode(count_mat) <- "numeric"

  low  <- min(rate_mat, na.rm = TRUE)
  high <- max(rate_mat, na.rm = TRUE)

  hdr_style <- openxlsx::createStyle(fontName = "Times New Roman", fontSize = 10,
                                     textDecoration = "bold", halign = "center", valign = "center")
  for (j in seq_along(drugs)) {
    openxlsx::writeData(wb, "Coverage heatmap", drugs[j], startRow = 4, startCol = j + 1L, colNames = FALSE)
    openxlsx::addStyle(wb, "Coverage heatmap", hdr_style, rows = 4, cols = j + 1L)
  }
  openxlsx::setColWidths(wb, "Coverage heatmap", cols = 1, widths = 10)
  openxlsx::setColWidths(wb, "Coverage heatmap", cols = 2:24, widths = 9.3)

  # Pre-build 100 fill styles for fast styling
  cell_styles <- lapply(palette, function(col_hex) {
    openxlsx::createStyle(fontName = "Times New Roman", fontSize = 10,
                          halign = "center", valign = "center", wrapText = TRUE,
                          fgFill = col_hex)
  })
  row_hdr_style <- openxlsx::createStyle(fontName = "Times New Roman", fontSize = 10,
                                         textDecoration = "bold", valign = "center")

  for (i in seq_len(23L)) {
    r_idx <- i + 4L
    openxlsx::writeData(wb, "Coverage heatmap", rates_df[[1L]][i], startRow = r_idx, startCol = 1, colNames = FALSE)
    openxlsx::addStyle(wb, "Coverage heatmap", row_hdr_style, rows = r_idx, cols = 1)
    openxlsx::setRowHeights(wb, "Coverage heatmap", rows = r_idx, heights = 31)
    for (j in seq_len(23L)) {
      val <- rate_mat[i, j]
      cnt <- count_mat[i, j]
      txt <- sprintf("%.1f\nn=%d", 100 * val, as.integer(round(cnt)))
      p_idx <- max(1L, min(length(palette), as.integer(round(1 + (val - low) / (high - low) * (length(palette) - 1L)))))
      openxlsx::writeData(wb, "Coverage heatmap", txt, startRow = r_idx, startCol = j + 1L, colNames = FALSE)
      openxlsx::addStyle(wb, "Coverage heatmap", cell_styles[[p_idx]], rows = r_idx, cols = j + 1L)
    }
  }
  openxlsx::writeData(wb, "Coverage heatmap",
                      "Darker blue denotes higher coverage. The Drug key sheet gives all generic names.",
                      startRow = 29, startCol = 1, colNames = FALSE)

  # Sheet 2: Drug key
  openxlsx::addWorksheet(wb, "Drug key", gridLines = FALSE)
  openxlsx::writeData(wb, "Drug key", drug_key_df, startRow = 1, startCol = 1)

  # Sheet 3: Evaluable n
  openxlsx::addWorksheet(wb, "Evaluable n", gridLines = FALSE)
  openxlsx::writeData(wb, "Evaluable n", counts_df, startRow = 1, startCol = 1)

  # Sheet 4: Colour legend
  legend_df <- data.frame(
    Palette_index = seq_along(palette),
    Coverage_pct  = 100 * (low + (high - low) * (seq_along(palette) - 1L) / (length(palette) - 1L)),
    Fill_hex      = palette,
    stringsAsFactors = FALSE
  )
  openxlsx::addWorksheet(wb, "Colour legend", gridLines = FALSE)
  openxlsx::writeData(wb, "Colour legend", legend_df, startRow = 1, startCol = 1)
  for (k in seq_along(palette)) {
    openxlsx::addStyle(wb, "Colour legend", cell_styles[[k]], rows = k + 1L, cols = 3)
  }

  # This workbook has no drawings or printer settings. Remove openxlsx default
  # orphan relationships so standard OOXML readers can open it without repair.
  for (i in seq_along(wb$worksheets)) {
    wb$worksheets[[i]]$drawing <- character()
    wb$worksheets_rels[[i]] <- character()
  }
  openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)
  message("[PASS] Exported pure-R S2_Heatmap.xlsx (529 cells, 4 sheets): ", out_xlsx)
  invisible(out_xlsx)
}

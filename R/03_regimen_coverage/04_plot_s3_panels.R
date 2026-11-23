if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
# ==============================================================================
# Module 03: Render S3 Fig Four Separate 600-dpi TIFF Panels
# Usage: Rscript --vanilla 04_plot_s3_panels.R S3_Fig_Data.csv OUTPUT_DIR panels
# ==============================================================================
# Display-only adaptation of 11.9/S3Fig_重绘.R.
# No susceptibility classification, fits, tests, or reweighting.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L, args[3] %in% c("combined", "panels"))
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
})
data <- read.csv(args[1], check.names = FALSE)
regimens <- c("PEN + CAZ" = "Penicillin + Ceftazidime",
              "PEN + CTX" = "Penicillin + Cefotaxime",
              "AMP + CAZ" = "Ampicillin + Ceftazidime",
              "AMP + CTX" = "Ampicillin + Cefotaxime")
stopifnot(nrow(data) == 40L, !anyDuplicated(data[c("Regimen", "Year")]),
          setequal(data$Regimen, names(regimens)),
          setequal(data$Year, 2015:2024),
          all(is.finite(data$Unweighted_coverage)),
          all(is.finite(data$Weighted_coverage)),
          all(data$Unweighted_coverage >= 0 & data$Unweighted_coverage <= 1),
          all(data$Weighted_coverage >= 0 & data$Weighted_coverage <= 1))
yearly_all <- data %>% transmute(
  year = Year, regimen_en = unname(regimens[Regimen]),
  coverage_unw = Unweighted_coverage, coverage_w = Weighted_coverage
)

plot_coverage_trend <- function(yearly_df, regimen_name, y_limits = c(0, 1)) {
  dat <- yearly_df %>% filter(regimen_en == regimen_name) %>% arrange(year)
  dat_long <- dat %>%
    select(year, coverage_unw, coverage_w) %>%
    pivot_longer(starts_with("coverage_"), names_to = "weighting", values_to = "coverage") %>%
    mutate(weighting = recode(weighting, coverage_unw = "Unweighted", coverage_w = "Weighted"))
  yrs <- seq(min(dat_long$year), max(dat_long$year))
  dat_long <- complete(dat_long, year = yrs, weighting)
  ggplot(dat_long, aes(x = year, y = coverage, group = weighting,
                       linetype = weighting, shape = weighting)) +
    geom_line(linewidth = 0.9, na.rm = TRUE) +
    geom_point(size = 2.2, na.rm = TRUE) +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       limits = y_limits, breaks = pretty(y_limits, n = 6)) +
    scale_x_continuous(breaks = yrs, minor_breaks = NULL) +
    labs(title = regimen_name, x = "Year", y = "Coverage",
         linetype = NULL, shape = NULL) +
    theme_bw(base_size = 13) +
    theme(legend.position = "top", panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"))
}
plots_s3 <- lapply(unname(regimens), function(reg) plot_coverage_trend(yearly_all, reg))
draw_s3 <- function() {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 2)))
  for (i in seq_along(plots_s3)) {
    print(plots_s3[[i]], newpage = FALSE,
          vp = grid::viewport(layout.pos.row = (i-1) %/% 2+1, layout.pos.col = (i-1) %% 2+1))
  }
  grid::popViewport()
}
if (args[3] == "combined") {
  grDevices::png(file.path(args[2], "S3_Fig.png"), width = 16, height = 10,
                 units = "in", res = 300, bg = "white")
  draw_s3()
  invisible(grDevices::dev.off())
  grDevices::tiff(file.path(args[2], "S3_Fig.tiff"), width = 16, height = 10,
                  units = "in", res = 300, compression = "lzw", bg = "white")
  draw_s3()
  invisible(grDevices::dev.off())
  capture.output(sessionInfo(), file = file.path(args[2], "S3_graphics_session.txt"))
}

# Export the existing ggplot objects at their original 8 x 5 inch panel size.
# PNG matches the 300 dpi combined figure; TIFF is publication resolution.
# Start in a separate R process to reproduce the accepted Windows device state.
if (args[3] == "panels") {
  for (i in seq_along(plots_s3)) {
    stem <- file.path(args[2], paste0("S3_panel_", i))
    grDevices::png(paste0(stem, ".png"), width = 8, height = 5,
                   units = "in", res = 300, bg = "white")
    print(plots_s3[[i]])
    invisible(grDevices::dev.off())

    grDevices::tiff(paste0(stem, ".tiff"), width = 8, height = 5,
                    units = "in", res = 600, compression = "lzw", bg = "white")
    print(plots_s3[[i]])
    invisible(grDevices::dev.off())
  }
}

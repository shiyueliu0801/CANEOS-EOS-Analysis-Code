# ==============================================================================
# CANEOS R Package Dependency Installer
# ==============================================================================
# For exact version reproduction, use `renv::restore()` with `renv.lock`.
# Otherwise, run `Rscript --vanilla install_packages.R` to install all required
# CRAN packages.
# ==============================================================================

required_packages <- c(
  # Core statistical analysis & complex survey design
  "survey", "readr", "digest", "dplyr", "tidyr", "stringr", "scales",
  # Excel & native editable vector PowerPoint export (S2 Fig Excel, Figure 2 PPTX)
  "openxlsx", "writexl",
  # Publication graphics & cartography (Fig 2, Fig 3, S1 Fig, S3 Fig)
  "ggplot2", "patchwork", "cowplot", "ragg", "ggtext", "systemfonts", "sf",
  # Reproducible environment management
  "renv"
)

missing_req <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_req) > 0L) {
  message("Installing missing required packages: ", paste(missing_req, collapse = ", "))
  install.packages(missing_req, repos = "https://cloud.r-project.org")
} else {
  message("All core R packages are already installed.")
}

# Optional Figure 2 PPTX: install.packages(c("officer", "rvg")) only when requested.

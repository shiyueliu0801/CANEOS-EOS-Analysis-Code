module_started <- Sys.time()
if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
# ==============================================================================
# Module 03 Entrypoint: Empirical Regimen Coverage, S2 Heatmap & S3 Panels
# Final manuscript outputs: S9, S10, S14 Tables; Figure 2 data; S2 Fig (Excel);
#                           S3 Fig (4 x 600-dpi TIFF panels)
# Usage: Rscript --vanilla R/03_regimen_coverage/run.R [--config PATH]
# ==============================================================================

local({
  cli <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  mod_dir <- if (length(cli) >= 1L) dirname(normalizePath(sub("^--file=", "", cli[1L]), winslash = "/", mustWork = TRUE)) else file.path(getwd(), "R", "03_regimen_coverage")
  repo_root <- normalizePath(file.path(mod_dir, "..", ".."), winslash = "/", mustWork = TRUE)

  source(file.path(repo_root, "R", "common", "io_helpers.R"), local = TRUE, encoding = "UTF-8")
  source(file.path(repo_root, "R", "common", "validation_helpers.R"), local = TRUE, encoding = "UTF-8")
  source(file.path(repo_root, "R", "common", "survey_helpers.R"), local = TRUE, encoding = "UTF-8")

  ensure_utf8_locale()
  require_packages(c("survey", "readr", "digest", "dplyr", "writexl", "openxlsx", "ggplot2", "tidyr", "scales"))
  set_caneos_survey_options()

  cfg <- load_caneos_config(repo_root)
  static_cfg <- load_static_configs(repo_root)
  output_dir <- ensure_module_out_dir(cfg, "03_regimen_coverage")

  # Step 1 & 2: Statistical coverage analysis (S9, S10, S14, Fig 2 backing data, S2 & S3 data)
  env <- new.env(parent = globalenv())
  env$cfg <- cfg
  env$static_cfg <- static_cfg
  env$output_dir <- output_dir
  env$mod_dir <- mod_dir
  env$read_csv_safe <- read_csv_safe
  env$write_csv_utf8 <- write_csv_utf8
  source(file.path(mod_dir, "02_coverage_analysis.R"), local = env, encoding = "UTF-8")

  # Step 3: Pure-R S2 Fig Heatmap Excel workbook
  source(file.path(mod_dir, "03_export_s2_excel.R"), local = TRUE, encoding = "UTF-8")
  export_s2_heatmap_excel(
    rates_df    = env$s2_heatmap_rates,
    counts_df   = env$s2_heatmap_n,
    drug_key_df = env$s2_drug_key,
    palette     = static_cfg$get_s2_heatmap_palette(),
    out_xlsx    = file.path(output_dir, "S2_Heatmap.xlsx")
  )

  # Step 4: Render S3 Fig four separate 600-dpi TIFF panels
  s3_data <- file.path(output_dir, "S3_Fig_Data.csv")
  s3_out  <- file.path(output_dir, "S3_Fig_Panels")
  dir.create(s3_out, recursive = TRUE, showWarnings = FALSE)
  rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  status <- system2(rscript, shQuote(c("--vanilla", file.path(mod_dir, "04_plot_s3_panels.R"), s3_data, s3_out, "panels")))
  if (status != 0L) stop("S3 panel rendering failed with status ", status, call. = FALSE)

  cat("Module 03 completed. Outputs written to:", output_dir, "\n")
  assert_module_outputs(cfg,"03",module_started)
})

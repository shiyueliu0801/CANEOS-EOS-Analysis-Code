module_started <- Sys.time()
if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
# ==============================================================================
# Module 01 Entrypoint: Baseline Characteristics, Incidence & Trends (BH-FDR)
# Final manuscript outputs: Table 1, S3 Table, S11 Table
# Usage: Rscript --vanilla R/01_baseline_incidence/run.R [--config PATH]
# ==============================================================================

local({
  cli <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  mod_dir <- if (length(cli) >= 1L) dirname(normalizePath(sub("^--file=", "", cli[1L]), winslash = "/", mustWork = TRUE)) else file.path(getwd(), "R", "01_baseline_incidence")
  repo_root <- normalizePath(file.path(mod_dir, "..", ".."), winslash = "/", mustWork = TRUE)

  source(file.path(repo_root, "R", "common", "io_helpers.R"), local = TRUE, encoding = "UTF-8")
  source(file.path(repo_root, "R", "common", "validation_helpers.R"), local = TRUE, encoding = "UTF-8")
  source(file.path(repo_root, "R", "common", "survey_helpers.R"), local = TRUE, encoding = "UTF-8")

  ensure_utf8_locale()
  require_packages(c("readr", "survey"))
  set_caneos_survey_options()
  options(digits = 15)

  cfg <- load_caneos_config(repo_root)
  static_cfg <- load_static_configs(repo_root)
  out_dir <- ensure_module_out_dir(cfg, "01_baseline_incidence")

  analysis <- new.env(parent = globalenv())
  analysis$cfg <- cfg
  analysis$static_cfg <- static_cfg
  analysis$output_dir <- out_dir
  analysis$project <- mod_dir
  analysis$read_csv_safe <- read_csv_safe
  analysis$write_csv_utf8 <- write_csv_utf8
  analysis$bh_adjust_planned <- bh_adjust_planned
  analysis$format_p_display <- format_p_display

  for (script in c("00_shared.R", "01_prepare.R", "02_incidence.R",
                   "03_baseline.R", "04_trends.R", "05_bh_fdr.R", "06_supporting_tables.R")) {
    source(file.path(mod_dir, script), local = analysis, encoding = "UTF-8")
  }
  # Write verified final-numbered aliases alongside original stems (fail fast if any source is missing)
  aliases <- list(
    c("S4_Incidence.csv",           "S3_Table_Hospital_Incidence.csv"),
    c("S4_Incidence_Estimates.csv", "S3_Table_Hospital_Incidence_Estimates.csv"),
    c("Baseline_Trends.csv",        "S11_Table_Baseline_Trends.csv")
  )
  for (pair in aliases) {
    src_f <- file.path(out_dir, pair[1L])
    dst_f <- file.path(out_dir, pair[2L])
    if (!file.exists(src_f)) stop("Expected Module 01 output missing: ", src_f, call. = FALSE)
    ok <- file.copy(src_f, dst_f, overwrite = TRUE)
    if (!isTRUE(ok) || !file.exists(dst_f)) stop("Failed to create output alias: ", dst_f, call. = FALSE)
  }
  cat("Module 01 completed. Outputs written to:", out_dir, "\n")
  assert_module_outputs(cfg,"01",module_started)
})

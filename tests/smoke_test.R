if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
# ==============================================================================
# CANEOS Level-1 Code-Layer Smoke Test (Zero Real Study Data Required)
# ==============================================================================
# Verifies:
#   1. Every .R script in the repository parses with zero syntax errors.
#   2. All declared required R packages are available.
#   3. All static methodological configs load with expected exact dimensions.
#   4. BH-FDR helper reproduces planned-M conservative adjustment exactly.
#   5. Common UTF-8 / GB18030 I/O helpers and config loader work on temporary
#      synthetic files created strictly inside tempdir().
#
# Note: Passing this smoke test confirms software integrity only; reproducing
# submitted CANEOS numerical results requires authorised study inputs.
# ==============================================================================

local({
  cli <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  repo_root <- if (length(cli) >= 1L) {
    normalizePath(file.path(dirname(sub("^--file=", "", cli[1L])), ".."), winslash = "/", mustWork = TRUE)
  } else {
    normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  }

  source(file.path(repo_root, "R", "common", "io_helpers.R"), local = TRUE, encoding = "UTF-8")
  source(file.path(repo_root, "R", "common", "validation_helpers.R"), local = TRUE, encoding = "UTF-8")
  source(file.path(repo_root, "R", "common", "survey_helpers.R"), local = TRUE, encoding = "UTF-8")
  ensure_utf8_locale()

  cat("[1/6] Checking syntax of all R scripts...\n")
  r_files <- list.files(repo_root, pattern = "\\.[Rr]$", recursive = TRUE, full.names = TRUE)
  r_files <- r_files[!grepl("paths\\.local\\.R$", r_files)]
  for (f in r_files) {
    tryCatch(parse(file = f, encoding = "UTF-8"),
             error = function(e) stop("Syntax error in ", f, ": ", conditionMessage(e), call. = FALSE))
  }
  cat(sprintf("      PASS: %d R scripts parsed cleanly.\n", length(r_files)))

  cat("[2/6] Checking required R package availability...\n")
  require_packages(c("survey", "readr", "digest", "dplyr", "tidyr",
                     "stringr", "scales", "openxlsx", "writexl",
                     "ggplot2", "patchwork", "cowplot", "ragg", "ggtext",
                     "systemfonts", "sf"))
  cat("      PASS: All 16 required CRAN packages are installed.\n")

  cat("[3/6] Checking static methodological configurations...\n")
  sc <- load_static_configs(repo_root)
  stopifnot(
    nrow(sc$get_calibration_weights()) == 40L,
    nrow(sc$get_ast_intrinsic_rules()) == 77L,
    nrow(sc$get_coverage_intrinsic_rules()) == 83L,
    nrow(sc$get_s17_ast_applicability()) == 8L,
    nrow(sc$get_raw_field_dictionary()) == 95L,
    nrow(sc$get_baseline_case_columns()) == 44L,
    nrow(sc$get_ast_fields()) == 29L,
    nrow(sc$get_organism_groups()) == 16L,
    nrow(sc$get_ast_reporting_scope()) == 348L,
    nrow(sc$get_amr_trend_scope()) == 222L,
    nrow(sc$get_figure3_palette()) == 10L,
    length(sc$get_s2_heatmap_palette()) == 100L,
    nrow(sc$get_s1_province_anchors()) == 20L,
    nrow(sc$get_s19_missingness_dictionary()) == 62L
  )
  cat("      PASS: All 14 static configuration tables verified.\n")

  cat("[4/6] Checking BH-FDR planned-family adjustment math...\n")
  p_test <- c(0.001, 0.02, NA_real_, 0.04)
  q_got  <- bh_adjust_planned(p_test, planned_m = 5L)
  q_exp  <- c(stats::p.adjust(c(0.001, 0.02, 0.04), method = "BH", n = 5L)[1:2],
              NA_real_,
              stats::p.adjust(c(0.001, 0.02, 0.04), method = "BH", n = 5L)[3])
  stopifnot(isTRUE(all.equal(q_got, q_exp, tolerance = 1e-15)))
  cat("      PASS: bh_adjust_planned() matches stats::p.adjust(..., n = M) exactly.\n")

  cat("[5/6] Checking config resolution & UTF-8 Chinese CSV reading on synthetic tempdir...\n")
  td <- tempfile("caneos_smoke_")
  dir.create(td, recursive = TRUE)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  syn_csv <- file.path(td, "synthetic_cohort.csv")
  syn_df  <- data.frame("ID" = "SYN-001", "医院名称" = "Synthetic Site A",
                        "入院时间" = "2020", "分层" = "一层", check.names = FALSE)
  readr::write_csv(syn_df, syn_csv)
  syn_cfg <- file.path(td, "paths.smoke.R")
  writeLines(c(
    "get_caneos_paths <- function(repo_root) {",
    sprintf('  list(input_files = list(cohort_data = "%s"), output_dir = "%s")',
            normalizePath(syn_csv, winslash = "/"), normalizePath(file.path(td, "out"), winslash = "/", mustWork = FALSE)),
    "}"
  ), syn_cfg)
  cfg <- load_caneos_config(repo_root, config_path = syn_cfg)
  rd  <- read_csv_safe(cfg$input_files$cohort_data, encoding = "auto", expected_cols = c("ID", "医院名称", "入院时间", "分层"))
  stopifnot(identical(rd[["分层"]], "一层"), identical(rd[["医院名称"]], "Synthetic Site A"))
  cat("      PASS: Config loader and Chinese-header CSV reader verified on temporary synthetic input.\n")

  cat("[6/6] Running standalone Module 00 (S1 Appendix site priority ordering)...\n")
  rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  st00 <- system2(rscript, shQuote(c("--vanilla", file.path(repo_root, "R", "00_site_selection", "run.R"), "--config", syn_cfg)))
  stopifnot(st00 == 0L)
  s1_app <- readr::read_csv(file.path(td, "out", "00_site_selection", "S1_Appendix_Site_Priority_Order.csv"), show_col_types = FALSE)
  stopifnot(nrow(s1_app) == 94L, !anyDuplicated(s1_app$SiteCode),
            identical(as.integer(table(s1_app$Region)[c("Municipalities", "Eastern", "Central", "Western")]),
                      c(16L, 28L, 26L, 24L)))
  cat("      PASS: Module 00 reproduced all 94 anonymised candidate site rankings.\n")

  cat("\n=== ALL 6 SMOKE TESTS PASSED ===\n")
})

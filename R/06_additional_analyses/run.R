module_started <- Sys.time()
if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
# ==============================================================================
# Module 06 Entrypoint: Additional & Sensitivity Analyses
#   - Submodule 6A (regimen_ast_outcomes)    -> S17 Table
#   - Submodule 6B (missingness_evaluability)-> S19 Table
#   - Submodule 6C (deff_site_influence)     -> S20 Table, S21 Table
# Usage: Rscript --vanilla R/06_additional_analyses/run.R [--config PATH]
# ==============================================================================

local({
  cli <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  mod_dir <- if (length(cli) >= 1L) dirname(normalizePath(sub("^--file=", "", cli[1L]), winslash = "/", mustWork = TRUE)) else file.path(getwd(), "R", "06_additional_analyses")
  repo_root <- normalizePath(file.path(mod_dir, "..", ".."), winslash = "/", mustWork = TRUE)

  source(file.path(repo_root, "R", "common", "io_helpers.R"), local = TRUE, encoding = "UTF-8")
  source(file.path(repo_root, "R", "common", "validation_helpers.R"), local = TRUE, encoding = "UTF-8")
  source(file.path(repo_root, "R", "common", "survey_helpers.R"), local = TRUE, encoding = "UTF-8")

  ensure_utf8_locale()
  require_packages(c("survey", "readr", "digest"))
  set_caneos_survey_options()

  cfg <- load_caneos_config(repo_root)
  static_cfg <- load_static_configs(repo_root)
  base_out <- ensure_module_out_dir(cfg, "06_additional_analyses")
  data_file <- normalizePath(cfg$input_files$cohort_data, winslash = "/", mustWork = TRUE)

  # --- Submodule 6A: S17 Table (Initial Regimens & Paired CAZ/CTX AST vs Outcomes) ---
  cat("\n--- Running Module 06A: S17 Table (Regimen & Paired CAZ/CTX AST vs Outcomes) ---\n")
  out_s17 <- file.path(base_out, "s17_regimen_ast_outcomes")
  dir.create(out_s17, recursive = TRUE, showWarnings = FALSE)
  e6a <- new.env(parent = globalenv())
  source(file.path(mod_dir, "regimen_ast_outcomes", "01_data_preparation.R"), local = e6a, encoding = "UTF-8")
  source(file.path(mod_dir, "regimen_ast_outcomes", "02_analysis.R"),         local = e6a, encoding = "UTF-8")
  dat17 <- e6a$prepare_data(data_file, static_cfg$get_s17_ast_applicability())
  design17 <- survey::svydesign(ids = ~HospitalID, strata = ~Region, weights = ~w_unit,
                                data = dat17$infants, nest = TRUE)
  res17 <- e6a$analyse_s18(dat17, design17)
  for (nm in names(res17)) {
    out_nm <- sub("^S18_", "S17_", nm)
    readr::write_csv(res17[[nm]], file.path(out_s17, paste0(out_nm, ".csv")), na = "")
  }
  cat("Module 06A (S17 Table) completed:", nrow(res17$S18_overall_tests), "Rao-Scott tests.\n")

  # --- Submodule 6B: S19 Table (Missingness & Regimen Evaluability) ---
  cat("\n--- Running Module 06B: S19 Table (Missingness & Evaluability) ---\n")
  out_s19 <- file.path(base_out, "s19_missingness_evaluability")
  dir.create(out_s19, recursive = TRUE, showWarnings = FALSE)
  e6b <- new.env(parent = globalenv())
  e6b$static_cfg <- static_cfg
  for (f in c("00_shared_survey_functions.R", "01_data_preparation.R",
              "02_clinical_missingness.R", "03_coverage_evaluability.R",
              "04_combine_final_results.R")) {
    source(file.path(mod_dir, "missingness_evaluability", f), local = e6b, encoding = "UTF-8")
  }
  dat19 <- e6b$prepare_data(data_file, mod_dir)
  a19   <- e6b$analyse_clinical(dat19)
  b19   <- e6b$analyse_coverage(dat19, mod_dir)
  res19 <- e6b$combine_results(dat19, a19, b19)
  for (nm in names(res19)) {
    readr::write_csv(res19[[nm]], file.path(out_s19, paste0(nm, ".csv")), na = "")
  }
  static_cfg$copy_required(file.path(out_s19, "final_missingness_table.csv"),
            file.path(out_s19, "S19_Table_Missingness_Evaluability.csv"), overwrite = TRUE)
  cat("Module 06B (S19 Table) completed:", nrow(a19$overall), "clinical indicators;",
      nrow(b19$overall), "regimens.\n")

  # --- Submodule 6C: S20 & S21 Tables (Design Effects & Site 23 Sensitivity) ---
  cat("\n--- Running Module 06C: S20 & S21 Tables (DEFF & Site 23 Exclusion) ---\n")
  out_s20 <- file.path(base_out, "s20_s21_deff_site_influence")
  dir.create(out_s20, recursive = TRUE, showWarnings = FALSE)
  e6c <- new.env(parent = globalenv())
  for (f in c("01_data_preparation.R", "02_s20_deff.R",
              "03_s21_site23_sensitivity.R", "04_write_results.R")) {
    source(file.path(mod_dir, "deff_site_influence", f), local = e6c, encoding = "UTF-8")
  }
  inp20 <- e6c$load_inputs(cfg, static_cfg)
  d20   <- e6c$prepare_data(inp20)
  s20   <- e6c$compute_s21(d20$infant)
  s21   <- e6c$compute_s22(d20)
  e6c$write_results(out_s20, s20, s21, d20)
  # Final-manuscript numbered aliases (S20 Table = DEFF, S21 Table = Site 23 exclusion)
  static_cfg$copy_required(file.path(out_s20, "S21_final_results.csv"), file.path(out_s20, "S20_Table_DEFF_Final.csv"), overwrite = TRUE)
  static_cfg$copy_required(file.path(out_s20, "S22_final_results.csv"), file.path(out_s20, "S21_Table_Site23_Sensitivity_Final.csv"), overwrite = TRUE)
  cat("Module 06C (S20 & S21 Tables) completed.\n")

  cat("\nModule 06 completed. Outputs written to:", base_out, "\n")
  assert_module_outputs(cfg,"06",module_started)
})

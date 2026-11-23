module_started <- Sys.time()
if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
# ==============================================================================
# Module 02 Entrypoint: Pathogen Distribution, ESBL, AST & Temporal Trends
# Final manuscript outputs: Table 2, Table 3, S4, S5, S7, S12, S13 Tables
# Usage: Rscript --vanilla R/02_pathogens_ast/run.R [--config PATH]
# ==============================================================================

local({
  cli <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  mod_dir <- if (length(cli) >= 1L) dirname(normalizePath(sub("^--file=", "", cli[1L]), winslash = "/", mustWork = TRUE)) else file.path(getwd(), "R", "02_pathogens_ast")
  repo_root <- normalizePath(file.path(mod_dir, "..", ".."), winslash = "/", mustWork = TRUE)

  source(file.path(repo_root, "R", "common", "io_helpers.R"), local = TRUE, encoding = "UTF-8")
  source(file.path(repo_root, "R", "common", "validation_helpers.R"), local = TRUE, encoding = "UTF-8")
  source(file.path(repo_root, "R", "common", "survey_helpers.R"), local = TRUE, encoding = "UTF-8")

  ensure_utf8_locale()
  require_packages(c("survey", "readr"))
  set_caneos_survey_options()

  cfg <- load_caneos_config(repo_root)
  static_cfg <- load_static_configs(repo_root)
  output_dir <- ensure_module_out_dir(cfg, "02_pathogens_ast")

  inputs <- c(
    Cases = normalizePath(cfg$input_files$cohort_data, winslash = "/", mustWork = TRUE),
    Sites = normalizePath(cfg$input_files$site_map, winslash = "/", mustWork = TRUE)
  )

  for (script in c("01_inputs_and_design.R", "02_pathogens.R", "03_ast_and_esbl.R", "04_write_results.R")) {
    source(file.path(mod_dir, script), local = environment(), encoding = "UTF-8")
  }
  cat("Module 02 completed. Outputs written to:", output_dir, "\n")
  assert_module_outputs(cfg,"02",module_started)
})

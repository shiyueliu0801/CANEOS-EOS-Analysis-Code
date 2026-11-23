module_started <- Sys.time()
if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
# ==============================================================================
# Module 05 Entrypoint: CoNS-Related Sensitivity Analyses
# Final manuscript output: S18 Table
# Usage: Rscript --vanilla R/05_cons_sensitivity/run.R [--config PATH]
# ==============================================================================

cli <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
mod_dir <- if (length(cli) >= 1L) dirname(normalizePath(sub("^--file=", "", cli[1L]), winslash = "/", mustWork = TRUE)) else file.path(getwd(), "R", "05_cons_sensitivity")
repo_root <- normalizePath(file.path(mod_dir, "..", ".."), winslash = "/", mustWork = TRUE)

source(file.path(repo_root, "R", "common", "io_helpers.R"), encoding = "UTF-8")
source(file.path(repo_root, "R", "common", "validation_helpers.R"), encoding = "UTF-8")
source(file.path(repo_root, "R", "common", "survey_helpers.R"), encoding = "UTF-8")

ensure_utf8_locale()
require_packages(c("survey", "readr", "digest"))
set_caneos_survey_options()

cfg <- load_caneos_config(repo_root)
static_cfg <- load_static_configs(repo_root)
output <- ensure_module_out_dir(cfg, "05_cons_sensitivity")

for (script in c("01_inputs.R", "02_survey_design.R", "03_isolate_results.R",
                 "04_infant_results.R", "05_export_results.R")) {
  source(file.path(mod_dir, script), encoding = "UTF-8")
}
cat("Module 05 completed. Outputs written to:", output, "\n")
  assert_module_outputs(cfg,"05",module_started)

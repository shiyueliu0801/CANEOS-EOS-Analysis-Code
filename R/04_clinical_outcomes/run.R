module_started <- Sys.time()
if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
# ==============================================================================
# Module 04 Entrypoint: Clinical Management, Outcomes, Trends & BH-FDR
# Final manuscript outputs: Table 4, S6, S8, S15, S16 Tables; Fig 3 trajectory
# Usage: Rscript --vanilla R/04_clinical_outcomes/run.R [--config PATH]
# ==============================================================================

cli <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
mod_dir <- if (length(cli) >= 1L) dirname(normalizePath(sub("^--file=", "", cli[1L]), winslash = "/", mustWork = TRUE)) else file.path(getwd(), "R", "04_clinical_outcomes")
repo_root <- normalizePath(file.path(mod_dir, "..", ".."), winslash = "/", mustWork = TRUE)

source(file.path(repo_root, "R", "common", "io_helpers.R"), encoding = "UTF-8")
source(file.path(repo_root, "R", "common", "validation_helpers.R"), encoding = "UTF-8")
source(file.path(repo_root, "R", "common", "survey_helpers.R"), encoding = "UTF-8")

ensure_utf8_locale()
require_packages(c("readr", "survey", "digest"))
set_caneos_survey_options()

cfg <- load_caneos_config(repo_root)
output_path <- ensure_module_out_dir(cfg, "04_clinical_outcomes")
case_path <- normalizePath(cfg$input_files$cohort_data, winslash = "/", mustWork = TRUE)
site_path <- normalizePath(cfg$input_files$site_map,    winslash = "/", mustWork = TRUE)

for (script in c("01_cohort.R", "02_models.R", "03_clinical_comparisons.R",
                 "04_antimicrobial_use.R", "05_clinical_trends.R",
                 "06_antimicrobial_trends.R", "07_bh_fdr_and_trajectories.R")) {
  source(file.path(mod_dir, script), encoding = "UTF-8")
}

# Write final-manuscript numbered tables (and preserve audit tables)
release_tables <- list(
  Table4                         = Table4,
  S6_Table_Antimicrobial_Use     = S7,
  S8_Table_Risk_Ratios           = S9,
  S15_Table_Antimicrobial_Trends = S16,
  S16_Table_Clinical_Trends      = S17
)
audit_tables <- c(clinical_audit, antimicrobial_audit, clinical_trend_audit,
                  s16_results[setdiff(names(s16_results), "S16_Antimicrobial_Trends")])
audit_path <- file.path(output_path, "audit")
dir.create(audit_path, recursive = TRUE, showWarnings = FALSE)
for (nm in names(release_tables)) readr::write_csv(release_tables[[nm]], file.path(output_path, paste0(nm, ".csv")), na = "")
for (nm in names(audit_tables))   readr::write_csv(audit_tables[[nm]],   file.path(audit_path,  paste0(nm, ".csv")), na = "")

cat("Module 04 completed:", nrow(inf), "infants in", psu_count(inf), "hospitals. Outputs:", output_path, "\n")
  assert_module_outputs(cfg,"04",module_started)

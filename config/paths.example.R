# ==============================================================================
# CANEOS Public Code Repository — Local Input Path Configuration Template
# ==============================================================================
# Copy this file to `config/paths.local.R` (ignored by Git) OR save it outside
# the repository and set the environment variable `CANEOS_CONFIG=/path/to/file.R`
# (or pass `--config /path/to/file.R` to `run_all.R` or any module `run.R`).
#
# Do NOT commit real file paths or restricted study datasets to Git.
# Original raw Chinese filenames do NOT need to be renamed; map them via the
# logical aliases below.
# ==============================================================================

get_caneos_paths <- function(repo_root) {
  # Base directory holding authorised, restricted study inputs outside the repo
  data_dir <- Sys.getenv("CANEOS_DATA_DIR", unset = "/path/to/restricted_caneos_data")

  # Output directory (created automatically; recommended outside the repository)
  output_dir <- Sys.getenv("CANEOS_OUTPUT_DIR", unset = file.path(dirname(repo_root), "CANEOS_restricted_outputs"))

  list(
    # --------------------------------------------------------------------------
    # Controlled study inputs (obtained via institutional ethics approval;
    # see docs/RESTRICTED_DATA_ACCESS.md and docs/INPUT_SCHEMA.md)
    # --------------------------------------------------------------------------
    input_files = list(
      # Controlled 1196-row sidecar; values are never committed.
      coverage_scope = file.path(data_dir, "isolate_analysis_scope.csv"),

      # Core 95-column isolate/infant-level dataset (UTF-8, Chinese headers)
      cohort_data = file.path(data_dir, "controlled_cohort_data.csv"),

      # Hospital x year neonatal admissions denominator table (2015-2024;
      # GB18030 or UTF-8; used in Module 01 for Table 1, S3, S11)
      hospital_admissions = file.path(data_dir, "hospital_admissions.csv"),

      # Hospital name <-> anonymised Site ID (Site 1..Site 30) <-> Macro-region
      # crosswalk (UTF-8; used in Modules 02, 03, 04, 06)
      site_map = file.path(data_dir, "site_map.csv"),

      # Required for Module 07: Hospital <-> Province crosswalk (columns: Hospital, province_cn,
      # province_en, ne_adm1_code, ne_iso_3166_2) used ONLY by Module 07 (S1 Fig)
      # to dynamically aggregate provincial counts from cohort_data.
      hospital_province_map = file.path(data_dir, "hospital_province_map.csv")
    ),

    # Required: External pre-built S1 geography GeoPackage. Prepare it explicitly
    # with R/07_figures/s1_figure/prepare_geography.R as described in README.
    # Full execution fails if it is absent; it does not silently build or skip.
    geography_gpkg = file.path(output_dir, "07_figures", "s1_figure", "S1_geography.gpkg"),

    # Root directory where all analysis outputs and figures will be written
    # Extra workbooks contain restricted row-level data; enable only internally.
    export_internal_audit = FALSE,
    export_pptx = FALSE,
    output_dir = output_dir
  )
}

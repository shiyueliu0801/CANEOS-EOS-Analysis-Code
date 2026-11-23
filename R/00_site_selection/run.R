module_started <- Sys.time()
if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
# ==============================================================================
# Module 00 Entrypoint: Stratified Site-Selection Random Priority Ordering
# Final manuscript output: S1 Appendix (anonymised candidate site priority list)
# Requires ZERO restricted clinical or hospital data; runs completely standalone.
# Usage: Rscript --vanilla R/00_site_selection/run.R [--config PATH]
# ==============================================================================

local({
  cli <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  mod_dir <- if (length(cli) >= 1L) {
    dirname(normalizePath(sub("^--file=", "", cli[1L]), winslash = "/", mustWork = TRUE))
  } else {
    file.path(getwd(), "R", "00_site_selection")
  }
  repo_root <- normalizePath(file.path(mod_dir, "..", ".."), winslash = "/", mustWork = TRUE)

  source(file.path(repo_root, "R", "common", "io_helpers.R"), local = TRUE, encoding = "UTF-8")
  source(file.path(repo_root, "R", "common", "validation_helpers.R"), local = TRUE, encoding = "UTF-8")
  ensure_utf8_locale()
  require_packages(c("readr"))

  # Allow standalone execution without a local restricted-data config
  cfg <- tryCatch(
    load_caneos_config(repo_root),
    error = function(e) {
      out_root <- Sys.getenv("CANEOS_OUTPUT_DIR", unset = file.path(dirname(repo_root), "CANEOS_restricted_outputs"))
      dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
      list(repo_root = repo_root, output_dir = normalizePath(out_root, winslash = "/", mustWork = TRUE))
    }
  )
  out_dir <- ensure_module_out_dir(cfg, "00_site_selection")

  # Anonymised sampling-stratum classification for the 94 candidate SiteCodes
  county_codes <- c(
    # Eastern Region County-Level Unit (12 sites)
    "SITE-EA-2CF172", "SITE-EA-DBD24B", "SITE-EA-E1570C", "SITE-EA-BD64F2",
    "SITE-EA-54B3CF", "SITE-EA-33DA12", "SITE-EA-AB9AD4", "SITE-EA-73B697",
    "SITE-EA-EB2845", "SITE-EA-A6ED4E", "SITE-EA-13F8F8", "SITE-EA-798185",
    # Central Region County-Level Unit (3 sites)
    "SITE-CE-AD8E68", "SITE-CE-14B960", "SITE-CE-7B82A4",
    # Western Region County-Level Unit (5 sites)
    "SITE-WE-14FA00", "SITE-WE-B5D72A", "SITE-WE-2F0F70", "SITE-WE-90291F",
    "SITE-WE-2EEAA7"
  )

  run_stratum <- function(script_name, region_label) {
    env <- new.env(parent = globalenv())
    invisible(capture.output(source(file.path(mod_dir, script_name), local = env, encoding = "UTF-8")))
    res <- env$result
    code_col <- intersect(c("SiteCode", "Site_Code", "Site"), names(res))[1L]
    if (is.na(code_col)) stop("Missing site code column in ", script_name, call. = FALSE)
    sc <- as.character(res[[code_col]])
    strata <- if (region_label == "Municipalities") {
      rep("Municipalities", length(sc))
    } else {
      ifelse(sc %in% county_codes,
             paste(region_label, "Region County-Level Unit"),
             paste(region_label, "Region Provincial Capital District"))
    }
    data.frame(
      Region          = region_label,
      Sampling_strata = strata,
      SiteCode        = sc,
      Rank            = as.integer(res$Rank),
      stringsAsFactors = FALSE
    )
  }

  combined <- rbind(
    run_stratum("01_municipality_priority_order.R", "Municipalities"),
    run_stratum("02_eastern_priority_order.R",      "Eastern"),
    run_stratum("03_central_priority_order.R",      "Central"),
    run_stratum("04_western_priority_order.R",      "Western")
  )
  stopifnot(nrow(combined) == 94L, !anyDuplicated(combined$SiteCode))
  write_csv_utf8(combined, file.path(out_dir, "S1_Appendix_Site_Priority_Order.csv"))
  cat("Module 00 completed:", nrow(combined), "anonymised candidate sites ranked. Output:", out_dir, "\n")
  if (exists("assert_module_outputs", mode = "function")) {
    assert_module_outputs(cfg, "00", module_started)
  }
})

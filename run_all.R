if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
# ==============================================================================
# CANEOS Master Replication Entrypoint (Pure R: site ordering plus Modules 01-07)
# ==============================================================================
# Usage:
#   Rscript --vanilla run_all.R [--config /path/to/paths.local.R]
#                               [--modules 01,02,03,04,05,06,07]
#
# Configuration resolution order:
#   1. --config <path>
#   2. Environment variable CANEOS_CONFIG
#   3. Repository-local ignored file config/paths.local.R
# ==============================================================================

local({
  cli <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  repo_root <- if (length(cli) >= 1L) {
    dirname(normalizePath(sub("^--file=", "", cli[1L]), winslash = "/", mustWork = TRUE))
  } else {
    normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  }

  source(file.path(repo_root, "R", "common", "io_helpers.R"), local = TRUE, encoding = "UTF-8")
  ensure_utf8_locale()

  args <- commandArgs(trailingOnly = TRUE)
  cfg_path <- parse_config_arg(args)
  cfg <- load_caneos_config(repo_root, config_path = cfg_path)

  # Optional subset selection via --modules 01,04,07
  mod_filter <- NULL
  m_idx <- which(args == "--modules")
  if (length(m_idx) >= 1L && length(args) >= m_idx[1L] + 1L) {
    mod_filter <- strsplit(args[m_idx[1L] + 1L], ",", fixed = TRUE)[[1L]]
  }
  m_eq <- grep("^--modules=", args, value = TRUE)
  if (length(m_eq) >= 1L) {
    mod_filter <- strsplit(sub("^--modules=", "", m_eq[1L]), ",", fixed = TRUE)[[1L]]
  }

  modules <- list(
    "00" = "R/00_site_selection/run.R",
    "01" = "R/01_baseline_incidence/run.R",
    "02" = "R/02_pathogens_ast/run.R",
    "03" = "R/03_regimen_coverage/run.R",
    "04" = "R/04_clinical_outcomes/run.R",
    "05" = "R/05_cons_sensitivity/run.R",
    "06" = "R/06_additional_analyses/run.R",
    "07" = "R/07_figures/run.R"
  )

  selected <- if (is.null(mod_filter)) names(modules) else intersect(names(modules), trimws(mod_filter))
  if (length(selected) == 0L) stop("No valid modules selected. Choose from: 00,01,02,03,04,05,06,07.", call. = FALSE)

  rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  t0 <- Sys.time()
  cat("========================================================================\n")
  cat("CANEOS Master Analysis Pipeline\n")
  cat("Config :", cfg$config_file, "\n")
  cat("Outputs:", cfg$output_dir, "\n")
  cat("Modules:", paste(selected, collapse = ", "), "\n")
  cat("========================================================================\n\n")

  for (id in selected) {
    script_path <- file.path(repo_root, modules[[id]])
    cat(sprintf(">>> [%s] Running %s ...\n", format(Sys.time(), "%H:%M:%S"), modules[[id]]))
    module_started <- Sys.time()
    status <- system2(rscript, shQuote(c("--vanilla", script_path, "--config", cfg$config_file)))
    if (status != 0L) {
      stop(sprintf("Module %s (%s) exited with status %d.", id, modules[[id]], status), call. = FALSE)
    }
    assert_module_outputs(cfg, id, module_started)
    cat(sprintf("<<< [%s] Module %s executed; required artifacts verified.\n\n", format(Sys.time(), "%H:%M:%S"), id))
  }

  elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  cat(sprintf("All selected modules (%s) completed successfully in %.1f s.\n",
              paste(selected, collapse = ", "), elapsed))
})

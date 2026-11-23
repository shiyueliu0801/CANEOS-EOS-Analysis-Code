module_started <- Sys.time()
if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
# ==============================================================================
# Module 07 Entrypoint: Publication Figures (Pure R)
#   - Figure 2 (figure2/plot_figure2.R)     -> required 300-dpi TIFF + PDF; optional native PPTX (rvg)
#   - Figure 3 (figure3/plot_figure3.R)     -> 300-dpi TIFF + PDF (+ preview PNG)
#   - S1 Fig   (s1_figure/plot_s1_figure.R) -> 300-dpi TIFF + PDF
# Note: S2 Fig (TIFF + S2_Heatmap.xlsx) & S3 Fig (4 x 600-dpi TIFF) are generated in Module 03.
# Usage: Rscript --vanilla R/07_figures/run.R [--config PATH] [--skip-s1]
# ==============================================================================

local({
  cli <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  mod_dir <- if (length(cli) >= 1L) dirname(normalizePath(sub("^--file=", "", cli[1L]), winslash = "/", mustWork = TRUE)) else file.path(getwd(), "R", "07_figures")
  repo_root <- normalizePath(file.path(mod_dir, "..", ".."), winslash = "/", mustWork = TRUE)

  source(file.path(repo_root, "R", "common", "io_helpers.R"), local = TRUE, encoding = "UTF-8")
  source(file.path(repo_root, "R", "common", "validation_helpers.R"), local = TRUE, encoding = "UTF-8")

  ensure_utf8_locale()
  require_packages(c("readr", "dplyr", "tidyr", "ggplot2", "patchwork", "stringr",
                     "scales", "sf", "cowplot", "ragg", "ggtext", "systemfonts"))

  args <- commandArgs(trailingOnly = TRUE)
  skip_s1 <- "--skip-s1" %in% args

  cfg <- load_caneos_config(repo_root)
  static_cfg <- load_static_configs(repo_root)
  fig_out <- ensure_module_out_dir(cfg, "07_figures")

  # --- 1. Figure 2 (consumes Module 03 Figure2_Plot_Data.csv) ---
  f2_csv <- file.path(cfg$output_dir, "03_regimen_coverage", "Figure2_Plot_Data.csv")
  if (!file.exists(f2_csv)) stop("Run Module 03 first to generate ", f2_csv, call. = FALSE)
  source(file.path(mod_dir, "figure2", "plot_figure2.R"), local = TRUE, encoding = "UTF-8")
  f2_res <- render_figure2(f2_csv, file.path(fig_out, "figure2"), export_pptx = isTRUE(cfg$export_pptx))
  stopifnot(file.exists(f2_res$tiff), file.info(f2_res$tiff)$size > 0)

  # --- 2. Figure 3 (consumes Module 04 Figure3_Three_Phase_Trajectories.csv) ---
  f3_dir <- file.path(cfg$output_dir, "04_clinical_outcomes")
  if (!file.exists(file.path(f3_dir, "Figure3_Three_Phase_Trajectories.csv"))) {
    stop("Run Module 04 first to generate Figure3_Three_Phase_Trajectories.csv.", call. = FALSE)
  }
  e3 <- new.env(parent = globalenv())
  e3$repo_root <- repo_root
  e3$static_cfg <- static_cfg
  e3$fig3_data_dir <- f3_dir
  e3$fig3_out_dir  <- file.path(fig_out, "figure3")
  source(file.path(mod_dir, "figure3", "plot_figure3.R"), local = e3, encoding = "UTF-8")
  stopifnot(file.exists(file.path(fig_out, "figure3", "Figure3_Sankey_Final_v3_300dpi.tiff")))

  # --- 3. S1 Fig (required by default; fails fast unless --skip-s1 is explicitly passed) ---
  hpm  <- cfg$input_files$hospital_province_map
  gpkg <- cfg$geography_gpkg
  if (!skip_s1) {
    if (is.null(hpm) || !nzchar(hpm) || !file.exists(hpm)) {
      stop("Missing required S1 Fig input `hospital_province_map`: ", hpm,
           "\nProvide it in your local config or pass --skip-s1 to omit S1 Fig explicitly.", call. = FALSE)
    }
    if (is.null(gpkg) || !nzchar(gpkg) || !file.exists(gpkg)) {
      stop("Missing required S1 Fig geography file `geography_gpkg`: ", gpkg,
           "\nBuild it via R/07_figures/s1_figure/prepare_geography.R or pass --skip-s1.", call. = FALSE)
    }
    e_s1 <- new.env(parent = globalenv())
    source(file.path(mod_dir, "s1_figure", "plot_s1_figure.R"), local = e_s1, encoding = "UTF-8")
    s1_res <- e_s1$render_s1_figure(
      cohort_csv    = cfg$input_files$cohort_data,
      hosp_prov_csv = hpm,
      anchors_df    = static_cfg$get_s1_province_anchors(),
      gpkg_path     = gpkg,
      out_dir       = file.path(fig_out, "s1_figure")
    )
    stopifnot(file.exists(s1_res$tiff))
  } else {
    message("[NOTICE] S1 Fig skipped explicitly via --skip-s1 flag.")
  }

  cat("Module 07 completed. Outputs written to:", fig_out, "\n")
  assert_module_outputs(cfg,"07",module_started, skip_s1 = skip_s1)
})

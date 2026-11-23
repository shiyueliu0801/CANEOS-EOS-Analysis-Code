# ==============================================================================
# CANEOS Common Helpers: Configuration Resolution & Encoding-Safe I/O
# ==============================================================================

#' Locate repository root from a running Rscript or working directory
resolve_repo_root <- function(fallback = getwd()) {
  cli <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(cli) >= 1L) {
    script_dir <- dirname(normalizePath(sub("^--file=", "", cli[1L]), winslash = "/", mustWork = TRUE))
    cur <- script_dir
    for (i in seq_len(6L)) {
      if (file.exists(file.path(cur, "run_all.R")) && dir.exists(file.path(cur, "config"))) {
        return(normalizePath(cur, winslash = "/", mustWork = TRUE))
      }
      parent <- dirname(cur)
      if (identical(parent, cur)) break
      cur <- parent
    }
  }
  cur <- normalizePath(fallback, winslash = "/", mustWork = TRUE)
  if (file.exists(file.path(cur, "run_all.R")) && dir.exists(file.path(cur, "config"))) {
    return(cur)
  }
  stop("Unable to locate repository root containing run_all.R and config/.", call. = FALSE)
}

#' Parse optional --config <path> argument from commandArgs(trailingOnly = TRUE)
parse_config_arg <- function(args = commandArgs(trailingOnly = TRUE)) {
  idx <- which(args == "--config")
  if (length(idx) >= 1L) {
    if (length(args) < idx[1L] + 1L) stop("--config requires a file path.", call. = FALSE)
    return(args[idx[1L] + 1L])
  }
  eq <- grep("^--config=", args, value = TRUE)
  if (length(eq) >= 1L) return(sub("^--config=", "", eq[1L]))
  NULL
}

#' Gently request UTF-8 locale on Windows without hard-failing on unsupported systems
ensure_utf8_locale <- function() {
  if (.Platform$OS.type == "windows" && !isTRUE(l10n_info()[["UTF-8"]])) {
    for (loc in c("English_United States.utf8", "C.UTF-8", "Chinese (Simplified)_China.utf8")) {
      res <- suppressWarnings(tryCatch(Sys.setlocale("LC_CTYPE", loc), error = function(e) ""))
      if (nzchar(res) && isTRUE(l10n_info()[["UTF-8"]])) break
    }
  }
  invisible(isTRUE(l10n_info()[["UTF-8"]]))
}

#' Load local paths configuration using standard precedence:
#'   1. Explicit config_path argument (or --config CLI flag)
#'   2. Environment variable CANEOS_CONFIG
#'   3. Repository-local ignored file config/paths.local.R
load_caneos_config <- function(repo_root = resolve_repo_root(), config_path = NULL) {
  if (is.null(config_path)) config_path <- parse_config_arg()
  if (is.null(config_path) || !nzchar(config_path)) {
    env_cfg <- Sys.getenv("CANEOS_CONFIG", unset = "")
    if (nzchar(env_cfg)) config_path <- env_cfg
  }
  if (is.null(config_path) || !nzchar(config_path)) {
    local_cfg <- file.path(repo_root, "config", "paths.local.R")
    if (file.exists(local_cfg)) config_path <- local_cfg
  }
  if (is.null(config_path) || !file.exists(config_path)) {
    stop(
      paste0(
        "No local CANEOS configuration file found.\n",
        "Copy config/paths.example.R to config/paths.local.R (ignored by Git),\n",
        "or set CANEOS_CONFIG=/path/to/your_paths.R, or pass --config <path>.\n",
        "See docs/RESTRICTED_DATA_ACCESS.md and docs/INPUT_SCHEMA.md."
      ),
      call. = FALSE
    )
  }
  env <- new.env(parent = baseenv())
  source(config_path, local = env, encoding = "UTF-8")
  if (!exists("get_caneos_paths", envir = env, mode = "function")) {
    stop("Configuration file must define get_caneos_paths(repo_root): ", config_path, call. = FALSE)
  }
  cfg <- env$get_caneos_paths(repo_root)
  cfg$repo_root <- normalizePath(repo_root, winslash = "/", mustWork = TRUE)
  cfg$config_file <- normalizePath(config_path, winslash = "/", mustWork = TRUE)
  dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
  cfg$output_dir <- normalizePath(cfg$output_dir, winslash = "/", mustWork = TRUE)
  cfg
}

#' Load all static methodological configurations (zero patient/hospital records)
load_static_configs <- function(repo_root = resolve_repo_root()) {
  env <- new.env(parent = globalenv())
  for (f in c("calibration_weights.R", "intrinsic_resistance_rules.R", "field_mappings.R")) {
    source(file.path(repo_root, "config", f), local = env, encoding = "UTF-8")
  }
  env$copy_required <- copy_required
  env
}

#' Read CSV with explicit UTF-8 / GB18030 detection and Chinese header verification
read_csv_safe <- function(path, encoding = c("auto", "UTF-8", "GB18030"),
                          expected_cols = NULL) {
  encoding <- match.arg(encoding)
  if (!file.exists(path)) {
    stop("Required input file not found: ", path,
         "\nSee docs/RESTRICTED_DATA_ACCESS.md and config/paths.example.R.", call. = FALSE)
  }
  try_read <- function(enc) {
    df <- suppressWarnings(
      readr::read_csv(
        path,
        col_types = readr::cols(.default = readr::col_character()),
        locale = readr::locale(encoding = enc),
        na = character(),
        trim_ws = FALSE,
        show_col_types = FALSE,
        progress = FALSE,
        name_repair = "minimal"
      )
    )
    list(data = df, problems = nrow(readr::problems(df)))
  }
  if (encoding != "auto") {
    res <- try_read(encoding)
    if (res$problems > 0L) stop("CSV parsing failed (", encoding, "): ", basename(path), call. = FALSE)
    return(res$data)
  }
  for (enc in c("UTF-8", "GB18030")) {
    res <- tryCatch(try_read(enc), error = function(e) NULL)
    if (!is.null(res) && res$problems == 0L) {
      hdr <- names(res$data)
      if (all(validUTF8(hdr))) {
        if (is.null(expected_cols) || all(expected_cols %in% hdr)) {
          return(res$data)
        }
      }
    }
  }
  stop("Unable to decode CSV cleanly with UTF-8 or GB18030: ", basename(path), call. = FALSE)
}

#' Ensure module output directory exists and return normalized path
ensure_module_out_dir <- function(cfg, module_name) {
  d <- file.path(cfg$output_dir, module_name)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  normalizePath(d, winslash = "/", mustWork = TRUE)
}

#' Write UTF-8 CSV cleanly with full double precision where requested
write_csv_utf8 <- function(df, path, full_precision = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  out <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  if (isTRUE(full_precision)) {
    for (j in seq_along(out)) {
      if (is.numeric(out[[j]]) && !is.integer(out[[j]])) {
        out[[j]] <- ifelse(is.na(out[[j]]), "", sprintf("%.17g", out[[j]]))
      }
    }
  }
  readr::write_csv(out, path, na = "")
  invisible(path)
}

copy_required <- function(from, to, overwrite = TRUE) {
  if (!file.exists(from) || is.na(file.info(from)$size) || file.info(from)$size <= 0)
    stop("Required copy source missing or empty: ", from, call. = FALSE)
  if (!isTRUE(file.copy(from, to, overwrite = overwrite)) || !file.exists(to) ||
      !identical(unname(tools::md5sum(from)), unname(tools::md5sum(to))))
    stop("Required copy failed: ", basename(to), call. = FALSE)
  invisible(TRUE)
}

assert_module_outputs <- function(cfg, id, started, skip_s1 = FALSE) {
  e <- new.env(parent = baseenv())
  source(file.path(cfg$repo_root,"config","output_contract.R"),local=e)
  rel <- e$get_output_contract()[[id]]
  if (isTRUE(skip_s1)) rel <- rel[!grepl("/s1_figure/",rel,fixed=TRUE)]
  paths <- file.path(cfg$output_dir,rel);info <- file.info(paths)
  bad <- !file.exists(paths) | is.na(info$size) | info$size <= 0 | is.na(info$mtime) | info$mtime < started - 2
  if (any(bad)) stop("Required artifacts missing, empty or stale: ", paste(basename(paths[bad]),collapse=", "),call.=FALSE)
  invisible(TRUE)
}

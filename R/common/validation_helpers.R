# ==============================================================================
# CANEOS Common Helpers: Schema Validation, Continuation Rows, and BH-FDR
# ==============================================================================

#' Require a logical condition with an informative error message
require_check <- function(cond, msg) {
  if (!isTRUE(cond)) stop(msg, call. = FALSE)
  invisible(TRUE)
}

#' Check that required R packages are installed before running a module
require_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop("Missing required R package(s): ", paste(missing, collapse = ", "),
         "\nRun `Rscript --vanilla install_packages.R` or restore via renv.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Compute Benjamini-Hochberg q-values retaining the planned family size M.
#' When a planned hypothesis test is withheld or non-estimable (NA), its p-value
#' is conservatively treated as 1 for adjusting the estimable tests while the
#' returned q-value for the withheld test remains NA_real_.
bh_adjust_planned <- function(p_values, planned_m = length(p_values)) {
  p_num <- as.numeric(p_values)
  n_obs <- length(p_num)
  require_check(planned_m >= n_obs,
                sprintf("planned_m (%d) cannot be smaller than length(p_values) (%d).", planned_m, n_obs))
  valid <- is.finite(p_num)
  q_out <- rep(NA_real_, n_obs)
  if (!any(valid)) return(q_out)
  # Pad withheld/non-estimable and omitted planned tests with p = 1.0
  p_padded <- c(p_num[valid], rep(1.0, planned_m - sum(valid)))
  q_padded <- stats::p.adjust(p_padded, method = "BH")
  q_out[valid] <- q_padded[seq_len(sum(valid))]
  q_out
}

#' Format P / BH-q values to publication display convention (<0.001 or 3 decimals)
format_p_display <- function(x) {
  ifelse(is.finite(x), ifelse(x < 0.001, "<0.001", sprintf("%.3f", x)), "\u2014")
}

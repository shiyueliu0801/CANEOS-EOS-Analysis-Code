# ==============================================================================
# CANEOS Common Helpers: Complex Survey Design Defaults
# ==============================================================================

#' Apply locked global survey options used across all CANEOS modules
set_caneos_survey_options <- function() {
  options(
    survey.lonely.psu = "fail",
    survey.adjust.domain.lonely = FALSE,
    stringsAsFactors = FALSE
  )
  invisible(TRUE)
}



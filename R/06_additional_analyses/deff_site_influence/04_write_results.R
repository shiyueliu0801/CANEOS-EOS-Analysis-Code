# Export the manuscript displays and their full precision backing separately.
# Standard display precision follows the release specification. Text variants
# in the authoritative Word (e.g., 7.2 versus 7.20) compare numerically.
nfmt <- function(x) formatC(x, digits = 0, format = "f", big.mark = ",")
pm <- function(x, k) ifelse(is.na(x), "", sprintf(paste0("%.", k, "f"), x))
ci_display <- function(p, lo, hi) {
  ifelse(is.na(lo) | is.na(hi), paste0(pm(100*p, 1), " (CI unavailable)"),
    sprintf("%.1f (%.1f\u2013%.1f)", 100*p, 100*lo, 100*hi))
}

display_s21 <- function(x) data.frame(Panel = x$Panel, Row = x$Row,
  Hospitals = as.character(x$Hospitals), Positive_n_N = paste(nfmt(x$Positive), nfmt(x$N), sep = "/"),
  Proportion_CI = ci_display(x$Proportion, x$Lower95, x$Upper95), SE_pp = pm(100*x$SE, 2),
  CI_half_width_pp = pm(100*x$HalfWidth95, 2), DEFF = pm(x$DEFF, 2),
  Effective_n = pm(x$EffectiveN, 2), Precision_target_met = x$Precision_target_met,
  Change_in_proportion_pp = pm(x$Delta_proportion_pp, 2))

display_s22 <- function(x) {
  a <- x[x$Scenario == "A", ]
  ast_names <- c("Third-generation cephalosporins" = "3GC resistance",
    Ampicillin = "Ampicillin resistance", Gentamicin = "Gentamicin resistance",
    Carbapenems = "Carbapenem resistance", "Reported ESBL positivity" = "Reported ESBL positivity")
  labels <- a$Endpoint
  k <- a$Module == "AST_ESBL"
  labels[k] <- unname(ast_names[labels[k]])
  out <- data.frame(Row_ID = a$Row_ID, Module = a$Module, Region = a$Region, Endpoint = labels)
  for (s in c("A","B","C")) {
    z <- x[x$Scenario == s, ]
    z <- z[match(a$Row_ID,z$Row_ID), ]
    n <- nfmt(z$N)
    nonpair <- z$Module != "Paired"
    n[nonpair] <- paste(nfmt(z$Events[nonpair]), n[nonpair], sep = "/")
    out[[s]] <- paste(n, ci_display(z$Estimate,z$Low,z$High), sep = " / ")
  }
  out$C_minus_A_pp <- pm(a$C_minus_A_pp, 2)
  out
}

write_results <- function(out, a, b, d) {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  check(dir.exists(out), "Could not create output directory.")
  files <- list(S21_final_results = display_s21(a$results), S22_final_results = display_s22(b$results),
    S21_full_precision = a$results, S22_full_precision = b$results,
    S21_domain_validation = a$domains, S22_weight_calibration = b$weights,
    S22_scenarios = b$scenarios, Cohort_validation = d$counts,
    Input_hash_verification = d$hashes)
  checks <- data.frame(Check = c("Input SHA-256", "Infant-isolate linkage", "Hospital and time roster",
    "S21 DEFF Taylor and SRSWR checks", "S21 regional domain checks", "S21 influence designs",
    "S22 scenario design and calibration", "S22 AST and coverage denominators",
    "S22 paired covariance and Taylor check", "Output row counts"),
    Status = "PASS",
    Detail = c("Three external files and two bundled rule tables match input lock",
      sprintf("%d infants; %d isolates; %d infants with two isolates", nrow(d$infant), nrow(d$iso), d$double_infections),
      "30 PSUs; four strata; Central Site 17–23; calendar years 2015–2024",
      "Analytical hospital-score variance, SRSWR and effective n agree with survey for all 15 rows",
      "All four domains agree with standalone designs, including effective n and df",
      "Seven Central exclusions and one overall exclusion each rebuild the design",
      "A/B/C; 40 cells with preserved within-year birth shares in C; 30/29/29 PSUs",
      "E. coli valid AST and ESBL denominators checked; Strict-S and three-state counts checked",
      "Common denominator; paired covariance, direct weighting, survey GLM and manual PSU scores agree",
      "S21: 5 + 10 rows; S22: 8 + 5 + 6 + 4 rows, each in A/B/C"))
  files$verification_results <- checks
  for (k in names(files)) readr::write_csv(files[[k]], file.path(out, paste0(k, ".csv")), na = "")
  writeLines(capture.output(sessionInfo()), file.path(out, "R_sessionInfo.txt"))
}

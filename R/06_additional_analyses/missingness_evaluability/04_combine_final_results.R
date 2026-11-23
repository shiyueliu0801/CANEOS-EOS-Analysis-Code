fmt_p <- function(x) ifelse(is.na(x), "—", ifelse(x < 0.001, "<0.001", sprintf("%.3f", x)))
fmt_range <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return("—")
  paste(sprintf("%.1f", range(x)), collapse = "–")
}
combine_results <- function(dat, a, b) {
  clinical <- a$overall
  clinical$Module <- "Clinical"
  clinical$Display <- sprintf("%d/%d (%.1f)", clinical$Missing, clinical$Eligible, clinical$Percent)
  coverage <- b$overall
  coverage$Module <- "Coverage"
  coverage$Display <- sprintf("%d/%d (%.1f)", coverage$Unknown, coverage$Records, coverage$WeightedPercent)
  tests <- rbind(a$tests, b$tests)
  table_rows <- range_rows <- list()
  for (module in c("Clinical", "Coverage")) {
    res <- if (module == "Clinical") clinical else coverage
    counts <- if (module == "Clinical") a$counts else b$counts
    for (i in seq_len(nrow(res))) {
      k <- res$Key[i]
      row <- data.frame(Module = module, Key = k,
        Variable = if (module == "Clinical") res$Variable[i] else k, MissingDisplay = res$Display[i])
      for (axis in c("Year", "Region", "Outcome")) {
        z <- counts[counts$Key == k & counts$Dimension == axis, ]
        t <- tests[tests$Module == module & tests$Key == k & tests$Dimension == axis &
          tests$Weighting == if (module == "Clinical") "Unit" else "S3", ]
        rates <- if (module == "Clinical") z$Percent else z$WeightedPercent
        row[[paste0(axis, "Range")]] <- fmt_range(rates)
        rates <- rates[is.finite(rates)]
        range_rows[[length(range_rows) + 1L]] <- data.frame(Module = module, Key = k,
          Dimension = axis, AvailableGroups = length(rates),
          MinPercent = if (length(rates)) min(rates) else NA_real_,
          MaxPercent = if (length(rates)) max(rates) else NA_real_)
        row[[paste0(axis, "P")]] <- fmt_p(t$P)
        row[[paste0(axis, "PAnnotated")]] <- paste0(fmt_p(t$P),
          if (t$SparseMissing && is.finite(t$P)) "†" else "")
        row[[paste0(axis, "q")]] <- fmt_p(t$q_BH)
        row[[paste0(axis, "Sparse")]] <- t$SparseMissing
      }
      table_rows[[length(table_rows) + 1L]] <- row
    }
  }
  # The clinical q columns retain the original internal sensitivity analysis;
  # the accepted 8-column clinical table reports nominal P values only.
  validation <- data.frame(Check = c("Input SHA-256", "Infants", "Isolates", "Additional isolates",
    "Hospitals", "Survey strata", "Known infant outcomes", "Known isolate outcomes",
    "Early death blank late records", "Early withdrawal blank late records"),
    Value = c(dat$hash, nrow(dat$infants), nrow(dat$context), length(dat$extra), nrow(dat$sites),
      nlevels(dat$context$region), sum(!is.na(dat$infants$outcome)), sum(!is.na(dat$context$outcome)),
      sum(a$states$abx_late == "Not evaluable"), sum(a$states$abx_late == "Uncertain")))
  denominators <- bind(lapply(c("Infant", "Isolate"), function(unit) {
    d <- if (unit == "Infant") dat$infants else dat$context
    bind(lapply(names(groups(d)), function(axis) {
      g <- groups(d)[[axis]]
      data.frame(Unit = unit, Dimension = axis, Level = levels(g), N = as.integer(table(g)))
    }))
  }))
  list(clinical_missingness_results = clinical, clinical_group_counts = a$counts,
    coverage_missingness_results = coverage, coverage_group_counts = b$counts,
    final_missingness_table = bind(table_rows), statistical_test_details = tests,
    data_validation = validation, cohort_denominators = denominators,
    missingness_ranges = bind(range_rows))
}

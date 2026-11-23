bind <- function(x) do.call(rbind, x)
capture_fit <- function(expr) {
  messages <- character()
  value <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleWarning")
  }), error = function(e) {messages <<- c(messages, conditionMessage(e)); NULL})
  list(value = value, diagnostic = paste(unique(messages), collapse = " | "))
}
make_design <- function(d, w) {
  d$analysis_weight <- w
  survey::svydesign(ids = ~hospital, strata = ~region, weights = ~analysis_weight,
    data = d, nest = TRUE)
}
groups <- function(d) list(Year = d$year_factor, Region = d$region, Outcome = d$outcome)
# The full design precedes every eligibility or outcome-domain subset.
rs_test <- function(d, eligible, y, group, weight, module, key, dimension,
                    circular = FALSE) {
  keep <- eligible & !is.na(group)
  z <- d[keep, , drop = FALSE]
  f <- droplevels(group[keep]); yy <- y[keep]
  psus <- length(unique(z$hospital)); nstrata <- length(unique(z$region))
  diagnostic <- ""
  df <- psus - nstrata
  F <- df1 <- df2 <- p <- NA_real_
  cells <- table(f, factor(yy, levels = 0:1))
  missing_psus <- length(unique(z$hospital[yy == 1L]))
  observed_psus <- length(unique(z$hospital[yy == 0L]))
  if (circular) diagnostic <- "Outcome group defined by this variable."
  else if (!length(yy) || length(unique(yy)) < 2L || nlevels(f) < 2L) {
    diagnostic <- "No variation in missingness or fewer than two eligible groups."
  } else if (module == "Clinical" && (missing_psus < 2L || observed_psus < 2L)) {
    diagnostic <- "Fewer than two contributing PSUs on one side of the missingness indicator."
  } else if (module == "Clinical" && df < nlevels(f) - 1L) {
    diagnostic <- "Insufficient design degrees of freedom for the requested contrasts."
  } else if (module == "Coverage" && (any(rowSums(cells) == 0) || any(colSums(cells) == 0) || df <= 0)) {
    diagnostic <- "Empty margin or insufficient design df"
  } else {
    full <- make_design(d, weight)
    full$variables$Group <- group
    full$variables$y <- if (module == "Clinical") y else factor(y, levels = 0:1)
    dom <- subset(full, keep)
    dom$variables$Group <- f
    formula <- if (module == "Clinical") ~y + Group else ~Group + y
    fit <- capture_fit(survey::svychisq(formula, dom, statistic = "F"))
    diagnostic <- fit$diagnostic
    if (!is.null(fit$value)) {
      v <- unname(c(fit$value$statistic, fit$value$parameter, fit$value$p.value))
      if (length(v) == 4L && all(is.finite(v)) && all(v[2:3] > 0)) {
        F <- v[1]; df1 <- v[2]; df2 <- v[3]; p <- v[4]
      } else diagnostic <- paste(diagnostic, "Non-finite Rao-Scott statistic or degrees of freedom.")
    }
    stopifnot(survey::degf(dom) == df)
  }
  hospitals <- table(factor(unique(z[c("hospital", "region")])$region, levels = regions))
  data.frame(Module = module, Key = key, Dimension = dimension,
    Weighting = if (all(weight == 1)) "Unit" else "S3", N = sum(keep), Missing = sum(yy),
    PSUs = psus, DesignDF = df, MissingPSUs = missing_psus,
    MunicipalityPSUs = hospitals[1], EasternPSUs = hospitals[2],
    CentralPSUs = hospitals[3], WesternPSUs = hospitals[4],
    F = F, DF1 = df1, DF2 = df2, P = p, q_BH = NA_real_,
    MinimumObservedCell = if (length(cells)) min(cells) else NA_real_,
    ZeroCells = sum(cells == 0), SparseMissing = sum(y[eligible]) < 10L,
    Status = if (circular) "Not applicable" else if (is.finite(p)) "Estimated" else "Not estimable",
    Diagnostic = diagnostic, row.names = NULL)
}
adjust_families <- function(tests) {
  family <- interaction(tests$Module, tests$Dimension, tests$Weighting, drop = TRUE)
  for (g in levels(family)) {
    j <- which(family == g & is.finite(tests$P))
    tests$q_BH[j] <- p.adjust(tests$P[j], method = "BH",
      n = if (tests$Module[j[1]] %in% "Coverage") 6L else length(j))
  }
  tests
}

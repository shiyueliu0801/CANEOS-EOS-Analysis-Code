# Fit the published continuous-year models using the 2015–2024 calendar year.
binary_problem <- function(d) {
  if (sum(d$x == 1) < 3L || sum(d$x == 0) < 3L) return("Fewer than three events or non-events.")
  if (length(unique(d$HospitalID[d$x == 1])) < 2L || length(unique(d$HospitalID[d$x == 0])) < 2L) {
    return("Event or non-event support in fewer than two hospitals.")
  }
  r0 <- range(d$t[d$x == 0])
  r1 <- range(d$t[d$x == 1])
  if (r0[2] <= r1[1] || r1[2] <= r0[1]) return("Complete or quasi-complete separation by calendar year.")
  ""
}
fit_trend <- function(d, type) {
  d <- d[!is.na(d$x), ]
  method <- switch(type, continuous = "Survey linear regression", binary = "Survey logistic regression",
                   ordinal = "Survey proportional-odds regression")
  effect <- switch(type, continuous = "Mean change", binary = "OR", ordinal = "Common OR")
  a <- trend_skeleton(d, method, effect)
  a$Iterations <- NA_integer_
  a$FittedMin <- a$FittedMax <- NA_real_
  if (type == "binary") {
    a$Events <- sum(d$x == 1)
    a$NonEvents <- sum(d$x == 0)
    a$EventHospitals <- length(unique(d$HospitalID[d$x == 1]))
    a$NonEventHospitals <- length(unique(d$HospitalID[d$x == 0]))
  }
  if (a$Years < 2L || a$DF <= 0 || a$N < 10L || length(unique(d$x)) < 2L) {
    a$Diagnostic <- "Insufficient variation, eligible observations or design degrees of freedom."
    return(a)
  }
  if (type == "binary") {
    problem <- binary_problem(d)
    if (nzchar(problem)) {
      a$Diagnostic <- problem
      return(a)
    }
  }
  z <- capture(switch(type,
    continuous = survey::svyglm(x ~ t, design(d), family = gaussian()),
    binary = survey::svyglm(x ~ t, design(d), family = quasibinomial(),
                            control = stats::glm.control(epsilon = 1e-11, maxit = 50)),
    ordinal = survey::svyolr(x ~ t, design(d), method = "logistic",
                            control = list(reltol = 1e-12, maxit = 1000))))
  a$Diagnostic <- z$Message
  if (is.null(z$Fit)) return(a)
  if (type == "binary") {
    a$Iterations <- z$Fit$iter
    a$FittedMin <- min(fitted(z$Fit))
    a$FittedMax <- max(fitted(z$Fit))
  }
  a$Converged <- if (type == "ordinal") z$Fit$convergence == 0L else isTRUE(z$Fit$converged)
  if (!a$Converged || grepl("0 or 1|rank-deficient", z$Message, fixed = FALSE)) {
    a$Diagnostic <- combine_text(a$Diagnostic, "Convergence, separation or rank check failed.")
    return(a)
  }
  beta <- unname(coef(z$Fit)["t"])
  se <- sqrt(vcov(z$Fit)["t", "t"])
  if (type != "continuous" && (!is.finite(beta) || !is.finite(se) || abs(beta) >= 5 || se >= 5)) {
    a$Diagnostic <- combine_text(a$Diagnostic, "Unstable log odds or standard error.")
    return(a)
  }
  fill_trend(a, beta, se)
}

ordinal_check <- function(d) {
  d <- d[!is.na(d$x), ]
  d$x <- droplevels(d$x)
  levels_x <- levels(d$x)
  k <- length(levels_x) - 1L
  if (k < 2L) return(list(P = NA_real_, Omnibus = NA_real_, Note = "Fewer than three observed ordinal categories.", Thresholds = list()))
  parts <- lapply(seq_len(k), function(i) {
    v <- d
    v$x <- as.integer(d$x) > i
    v$Cut <- factor(i, levels = seq_len(k))
    v
  })
  thresholds <- setNames(lapply(parts, fit_trend, type = "binary"), paste0("Above ", levels_x[seq_len(k)]))
  empty <- function(s) list(P = NA_real_, Omnibus = NA_real_, Note = s, Thresholds = thresholds)
  if (!all(vapply(thresholds, function(r) r$Status == "Estimated", logical(1)))) {
    return(empty("Proportional-odds assumption could not be checked at one or more sparse thresholds."))
  }
  df <- design_df(d) - k
  if (df <= 0) return(empty("Too few design degrees of freedom for the threshold slope comparison."))
  stack <- do.call(rbind, parts)
  z <- capture(survey::svyglm(x ~ Cut * t, design(stack), family = quasibinomial(),
                             control = stats::glm.control(epsilon = 1e-11, maxit = 50)))
  if (is.null(z$Fit) || !isTRUE(z$Fit$converged) || nzchar(z$Message)) return(empty(combine_text("Stacked model did not pass diagnostics.", z$Message)))
  beta <- coef(z$Fit)
  V <- vcov(z$Fit)
  i <- which(grepl(":t$", names(beta)))
  wald <- function(j) {
    A <- V[j, j, drop = FALSE]
    if (length(j) == 0L || any(!is.finite(A)) || qr(A)$rank < length(j)) return(NA_real_)
    stat <- as.numeric(crossprod(beta[j], solve(A, beta[j]))) / length(j)
    pf(stat, length(j), df, lower.tail = FALSE)
  }
  p <- wald(i)
  list(P = p, Omnibus = wald(c(which(names(beta) == "t"), i)),
       Note = if (is.na(p)) "Threshold covariance is not full rank." else "",
       Thresholds = thresholds)
}

# The order and labels follow the approved baseline-trend table body.
columns <- c(5L, 8:16, 18:35, 37:38, 40:43)
labels <- c("Sex (male; ref: female)", "Gestational age, weeks", "Gestational age category (extremely preterm to term)",
  "Postnatal age at EOS onset (<24, 24 to <48, 48–72 h)", "Birth weight, g",
  "Mode of delivery (caesarean; ref: vaginal)", "Pregnancy plurality (singleton; ref: multiple gestation)",
  "Apgar score at 1 min", "Apgar score at 5 min", "Apgar score at 10 min", "PROM", "PROM ≥18 h", "MSAF",
  "MSAF grade among MSAF-positive infants (Grade I to III)", "Maternal age, years",
  "Gestational diabetes", "Gestational hypertension", "Intrahepatic cholestasis of pregnancy",
  "Chorioamnionitis", "Antenatal corticosteroid exposure", "None or ≤2 h before delivery",
  "GBS-specific >2 h before delivery", "Broad-spectrum >2 h before delivery",
  "Conventional mechanical ventilation", "High-frequency ventilation", "CPAP", "NIPPV", "HFNC",
  "Glucocorticoid therapy", "Surfactant administration", "Seizures within 72 h of birth",
  "Meningitis", "RDS", "Vasoactive therapy within 72 h of birth")
require_data(length(columns) == length(labels), "Trend variable specification is inconsistent.")
trend_value <- function(j) {
  switch(as.character(j), `5` = binary(j, "男性", "女性"),
    `9` = ordered(coded(j, group_age), levels = group_age),
    `10` = ordered(coded(j, onset), levels = onset),
    `12` = binary(j, "剖宫产", "自然分娩"), `13` = binary(j, "单胎", "双胎"),
    `20` = msaf_binary, `21` = ordered(positive_grade, 1:3, labels = grade_labels),
    `41` = meningitis_binary,
    if (j %in% c(8L, 11L, 14:16, 22L)) num(j)
    else if (j %in% 23:26) binary(j, "有", "无")
    else binary(j))
}
format_trend <- function(label, region, fit, reported_p) {
  show <- fit$Status == "Estimated"
  data.frame(Characteristic = label, Region = region,
    `β (per year)` = if (show) fmt_dec(fit$Beta) else "—",
    SE = if (show) fmt_dec(fit$SE) else "—",
    `Estimate per year (95% CI)` = if (show) fmt_interval(fit$Annual, fit$AnnualLower, fit$AnnualUpper) else "—",
    `Estimate over 5 years (95% CI)` = if (show) fmt_interval(fit$Five, fit$FiveLower, fit$FiveUpper) else "—",
    `P value` = fmt_p(reported_p), check.names = FALSE)
}
trend_display <- lapply(seq_along(regions), function(i) {
  a <- incidence_trends[i, ]
  format_trend(if (i == 1L) "Hospital-based incidence of culture-confirmed EOS" else "",
               a$Region, a, a$P)
})
trend_fits <- list()
po_fits <- list()
for (j in seq_along(columns)) {
  col <- columns[j]
  d <- inf
  d$x <- trend_value(col)
  type <- if (is.ordered(d$x)) "ordinal" else if (col %in% c(8L, 11L, 14:16, 22L)) "continuous" else "binary"
  for (i in seq_along(regions)) {
    g <- regions[i]
    domain <- d[region_rows(d, g), ]
    fit <- fit_trend(domain, type)
    p <- fit$P
    if (type == "ordinal") {
      po <- ordinal_check(domain)
      fit$POTestP <- po$P
      fit$Diagnostic <- combine_text(fit$Diagnostic, po$Note)
      if (is.na(po$P) || po$P < .05 || fit$Status != "Estimated") {
        fit$Status <- "Common OR withheld"
        p <- po$Omnibus
      }
      for (threshold in names(po$Thresholds)) {
        po_fits[[length(po_fits) + 1L]] <- cbind(Variable = labels[j], Region = g,
           Threshold = threshold, CommonSlopeP = po$P, JointSlopeP = po$Omnibus, po$Thresholds[[threshold]])
      }
    }
    trend_fits[[length(trend_fits) + 1L]] <- cbind(Variable = labels[j], SourceField = case_columns[col],
      Region = g, ReportedP = p, Statistic = fit$Beta/fit$SE, fit)
    trend_display[[length(trend_display) + 1L]] <- format_trend(if (i == 1L) labels[j] else "", g, fit, p)
  }
}
write_result(do.call(rbind, trend_display), "Baseline_Trends")
write_result(do.call(rbind, trend_fits), "Baseline_Trend_Models")
write_result(do.call(rbind, po_fits), "Ordinal_Diagnostics")

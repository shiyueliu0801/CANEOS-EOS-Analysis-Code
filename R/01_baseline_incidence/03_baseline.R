# Compare available-case characteristics between the four macro-regions.
test_regions <- function(x) {
  d <- inf
  d$x <- x
  d <- d[!is.na(d$x), ]
  continuous <- is.numeric(x)
  base <- data.frame(Method = if (continuous) "Design-based Kruskal–Wallis rank test" else "Rao–Scott second-order F",
    Statistic = NA_real_, DF1 = NA_real_, DF2 = NA_real_, P = NA_real_,
    Hospitals = psu_count(d), N = nrow(d), Status = "Not estimable", Diagnostic = "")
  if (length(unique(d$Region)) != 4L || design_df(d) < 3L || length(unique(d$x)) < 2L) {
    base$Diagnostic <- "Insufficient represented strata, degrees of freedom or outcome variation."
    return(base)
  }
  if (!continuous) {
    d$x <- droplevels(d$x)
    h <- vapply(split(d$HospitalID, d$x), function(v) length(unique(v)), integer(1))
    if (min(h) < 2L) {
      base$Diagnostic <- "At least one category represented in fewer than two hospitals."
      return(base)
    }
    n <- table(d$x, d$Region)
    if (any(n < 5L)) base$Diagnostic <- "Small observed cell(s); category support verified across hospitals."
    t <- capture(survey::svychisq(~x + Region, design(d), statistic = "F"))
  } else {
    t <- capture(survey::svyranktest(x ~ Region, design(d), test = "KruskalWallis"))
  }
  base$Diagnostic <- combine_text(base$Diagnostic, t$Message)
  if (!is.null(t$Fit)) {
    stat <- t$Fit$statistic
    parameter <- t$Fit$parameter
      if (continuous) {
      base$DF1 <- as.numeric(stat[1])
      base$DF2 <- t$Fit$ddf
      base$Statistic <- as.numeric(parameter[1]) / base$DF1
    } else {
      base$Statistic <- as.numeric(stat[1])
      base$DF1 <- as.numeric(parameter[1])
      base$DF2 <- as.numeric(parameter[2])
    }
    p <- t$Fit$p.value
    if (is.finite(p) && p >= 0 && p <= 1 && all(is.finite(c(base$Statistic, base$DF1, base$DF2))) &&
        abs(p - pf(base$Statistic, base$DF1, base$DF2, lower.tail = FALSE)) < 1e-9) {
      base$P <- p
      base$Status <- "Estimated"
    }
  }
  if (base$Status != "Estimated") base$Diagnostic <- combine_text(base$Diagnostic, "Design test unavailable.")
  base
}
desc <- function(x, display) {
  valid <- x[!is.na(x)]
  n <- length(valid)
  if (!n) return("—")
  if (is.numeric(x)) {
    q <- quantile(valid, c(.25, .50, .75), type = 2)
    f <- function(z) format(round(z, 2), trim = TRUE, scientific = FALSE)
    return(paste0(f(q[2]), " (", f(q[1]), "–", f(q[3]), ")"))
  }
  k <- sum(valid == display)
  sprintf("%s (%.1f)", format(k, big.mark = ",", trim = TRUE), 100 * k/n)
}

coded <- function(j, levels) {
  x <- num(j)
  require_data(all(na.omit(x) %in% seq_along(levels)), paste("Invalid category code in", case_columns[j]))
  factor(x, levels = seq_along(levels), labels = levels)
}
binary_f <- function(j, yes = "是", no = "否") {
  factor(binary(j, yes, no), levels = 0:1, labels = c("No", "Yes"))
}
apgar <- function(j) {
  x <- num(j)
  require_data(all(na.omit(x) %in% 0:10), "Invalid Apgar score.")
  factor(ifelse(x < 7, "<7", "7–10"), levels = c("<7", "7–10"))
}
grade_labels <- c("Grade I", "Grade II", "Grade III")
group_age <- c("Extremely preterm (<28 wk)", "Very preterm (28 to <32 wk)",
               "Moderate to late preterm (32 to <37 wk)", "Term (≥37 wk)")
onset <- c("<24 h", "24 to <48 h", "48–72 h")
sex <- factor(v(5L), c("男性", "女性"), c("Male", "Female"))
delivery <- factor(v(12L), c("自然分娩", "剖宫产"), c("Vaginal delivery", "Caesarean delivery"))
plurality <- factor(v(13L), c("单胎", "双胎"), c("Singleton", "Multiple gestation"))
meningitis_class <- factor(ifelse(meningitis_code %in% c("疑诊", "可疑"), "Suspected",
  ifelse(meningitis_code == "是", "Confirmed", ifelse(meningitis_code == "否", "No", NA_character_))),
  levels = c("Confirmed", "Suspected", "No"))

region_tests <- list()
baseline_long <- list()
summary_values <- function(g, x, label) {
  v <- x[region_rows(inf, g)]
  v <- v[!is.na(v)]
  n <- length(v)
  if (is.numeric(v)) {
    q <- if (n) unname(quantile(v, c(.25, .50, .75), type = 2)) else rep(NA_real_, 3)
    shown <- if (n) sprintf("%s (%s–%s)", count_text(q[2L]), count_text(q[1L]), count_text(q[3L])) else "—"
    return(data.frame(Region = g, Category = "", ValidN = n, Events = NA_integer_, Percentage = NA_real_,
                      LowerQuartile = q[1L], Median = q[2L], UpperQuartile = q[3L], Display = shown))
  }
  events <- sum(v == label)
  pct <- if (n) 100 * events/n else NA_real_
  data.frame(Region = g, Category = label, ValidN = n, Events = events, Percentage = pct,
             LowerQuartile = NA_real_, Median = NA_real_, UpperQuartile = NA_real_,
             Display = if (n) sprintf("%s (%.1f)", count_text(events), pct) else "—")
}
characteristic <- function(label, x, table = "Table 1", show = if (is.factor(x)) levels(x) else "") {
  fit <- test_regions(x)
  region_tests[[length(region_tests) + 1L]] <<- cbind(Table = table, Variable = label, fit)
  labels <- if (is.factor(x)) levels(x) else show
  summaries <- do.call(rbind, lapply(labels, function(level) {
    do.call(rbind, lapply(regions, summary_values, x = x, label = level))
  }))
  baseline_long[[length(baseline_long) + 1L]] <<- cbind(Table = table, Variable = label, summaries)
  title <- sprintf("%s (n = %s)", label, count_text(fit$N))
  p <- fmt_p(fit$P)
  if (length(show) == 1L) return(pub_row(title, summaries$Display[summaries$Category == show], p))
  do.call(rbind, c(list(pub_row(title, p = p)), lapply(show, function(level) {
    pub_row(level, summaries$Display[summaries$Category == level])
  })))
}
indicator <- function(label, j, positive = "是", negative = "否", table = "Table 1") {
  characteristic(label, binary_f(j, positive, negative), table, show = "Yes")
}
t1 <- rbind(
  pub_row("Culture-confirmed EOS cases, n", count_text(pooled_eos$Cases)),
  pub_row("Incidence per 1,000 neonatal admissions (95% CI)",
          fmt_interval(pooled_eos$Rate, pooled_eos$Lower, pooled_eos$Upper), fmt_p(inc_region$P)),
  pub_row("Neonatal characteristics"),
  characteristic("Sex", sex),
  characteristic("Gestational age", coded(9L, group_age)),
  characteristic("Birth weight, g", num(11L)),
  characteristic("Postnatal age at EOS onset", coded(10L, onset)),
  characteristic("Mode of delivery", delivery),
  characteristic("Pregnancy plurality", plurality),
  characteristic("Apgar score at 1 min", apgar(14L)),
  characteristic("Apgar score at 5 min", apgar(15L)),
  characteristic("Apgar score at 10 min", apgar(16L)),
  pub_row("Maternal and obstetric characteristics"),
  characteristic("Maternal age, years", num(22L)),
  indicator("PROM", 18L),
  indicator("PROM ≥18 h", 19L),
  indicator("MSAF", 20L, "有", "无"),
  characteristic("MSAF grade", factor(positive_grade, levels = 1:3, labels = grade_labels)),
  indicator("Gestational diabetes", 23L, "有", "无"),
  indicator("Gestational hypertension", 24L, "有", "无"),
  indicator("Intrahepatic cholestasis of pregnancy", 25L, "有", "无"),
  indicator("Chorioamnionitis", 26L, "有", "无"),
  indicator("Antenatal corticosteroid exposure", 27L),
  pub_row("Intrapartum antibiotic exposure"),
  indicator("None or ≤2 h before delivery", 28L),
  indicator("GBS-specific >2 h before delivery", 29L),
  indicator("Broad-spectrum >2 h before delivery", 30L))
write_result(t1, "Table1_Baseline")

# These eleven early clinical and management characteristics were moved from Table 1 to Table 4.
t4_indicator <- function(label, j) indicator(label, j, table = "Table 4")
t4 <- rbind(
  pub_row("Early clinical features and management"),
  pub_row("Respiratory support within 72 h of birth"),
  t4_indicator("Conventional mechanical ventilation", 31L),
  t4_indicator("High-frequency ventilation", 32L),
  t4_indicator("CPAP", 33L),
  t4_indicator("NIPPV", 34L),
  t4_indicator("HFNC", 35L),
  t4_indicator("Glucocorticoid therapy", 37L),
  t4_indicator("Surfactant administration", 38L),
  t4_indicator("Seizures within 72 h of birth", 40L),
  characteristic("Meningitis", meningitis_class, "Table 4", show = c("Confirmed", "Suspected")),
  t4_indicator("RDS", 42L),
  t4_indicator("Vasoactive therapy within 72 h of birth", 43L))
write_result(t4, "Table4_Baseline_Section")
write_result(do.call(rbind, baseline_long), "Baseline_Descriptives")
write_result(do.call(rbind, region_tests), "Baseline_Regional_Tests")


# Final descriptive rows, full-design domains, Rao-Scott tests and one contrast.
format_count <- function(events, valid_n) sprintf("%d/%d (%.1f)", events, valid_n, 100 * events / valid_n)
format_p <- function(p) ifelse(p < .001, "<0.001", sprintf("%.3f", p))
eligible_field <- c(Enterobacterales = "Eligible_Enterobacterales", "E. coli" = "Eligible_Ecoli")
group_definitions <- function(z) {
  defs <- list(
    list("Initial regimens", "All infants", "3GC-containing", z$Initial %in% c("3GC", "3GC+PEN")),
    list("Initial regimens", "All infants", "3GC+PEN", z$Initial == "3GC+PEN"),
    list("Initial regimens", "All infants", "BSP", z$Initial == "BSP"),
    list("Initial regimens", "All infants", "Carbapenem-containing", z$Initial %in% c("CARB", "CARB+GLY/LZD")),
    list("Initial regimens", "All infants", "CARB", z$Initial == "CARB"))
  for (population in names(eligible_field)) {
    for (g in patterns) {
      defs[[length(defs) + 1L]] <- list("Susceptibility", population, g,
        z[[eligible_field[[population]]]] & z$Complete_AST & z$TestGroup %in% g)
    }
  }
  defs
}
descriptive_results <- function(z) {
  defs <- group_definitions(z)
  do.call(rbind, lapply(seq_along(defs), function(i) {
    d <- defs[[i]]
    do.call(rbind, lapply(names(endpoint_labels), function(k) {
      v <- z[[k]][d[[4L]]]
      n <- sum(!is.na(v))
      e <- sum(v, na.rm = TRUE)
      data.frame(Row = i, Section = d[[1L]], Population = d[[2L]], Group = d[[3L]],
        Infants = sum(d[[4L]]), Code = k, Indicator = endpoint_labels[[k]],
        Events = e, Valid_N = n, Percent = 100 * e/n, Display = format_count(e, n))
    }))
  }))
}
capture_fit <- function(expr) {
  warnings <- character()
  error <- ""
  fit <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) { error <<- conditionMessage(e); NULL })
  list(fit = fit, warnings = paste(unique(warnings), collapse = " | "), error = error)
}
design_support <- function(d) {
  z <- d$variables
  vapply(regions, function(r) length(unique(z$HospitalID[z$Region == r])), 0L)
}
design_information <- function(d) {
  h <- design_support(d)
  as.data.frame(c(list(N = nrow(d$variables), Hospitals = sum(h), Strata = sum(h > 0L),
    Design_df = survey::degf(d)), as.list(setNames(unname(h), paste0(regions, "_hospitals")))))
}
fit_association <- function(d, code) {
  capture_fit(survey::svychisq(reformulate(c("TestGroup", paste0("Outcome_", code))), d, statistic = "F"))
}
extract_test <- function(model) {
  if (is.null(model$fit)) return(c(F = NA_real_, Numerator_df = NA_real_, Denominator_df = NA_real_, P = NA_real_))
  c(F = as.numeric(model$fit$statistic), Numerator_df = as.numeric(model$fit$parameter[1L]),
    Denominator_df = as.numeric(model$fit$parameter[2L]), P = as.numeric(model$fit$p.value))
}
test_domain <- function(full_design, population, code) {
  z <- full_design$variables
  include <- z[[eligible_field[[population]]]] & z$Complete_AST & z[[paste0("OK_", code)]]
  d <- subset(full_design, include)
  v <- d$variables
  xt <- table(v$TestGroup, v[[paste0("Outcome_", code)]])
  require_true(all(dim(xt) == c(3L, 2L)) && sum(xt) == nrow(v), "Invalid 3-by-2 domain table.")
  expected <- outer(rowSums(xt), colSums(xt)) / sum(xt)
  model <- fit_association(d, code)
  statistics <- extract_test(model)
  require_true(all(is.finite(statistics)), paste("Non-estimable Rao-Scott test", population, code, model$error))
  # No singleton correction is applied. Record support before interpreting a fit.
  require_true(all(design_support(d) >= 2L), "A domain has a singleton or empty stratum. Review before changing survey options.")
  cells <- do.call(rbind, lapply(seq_along(patterns), function(i) {
    group <- patterns[i]
    sub <- subset(d, TestGroup == group)
    y <- sub$variables[[code]]
    h <- sub$variables$HospitalID
    mean_fit <- survey::svymean(reformulate(code), sub)
    data.frame(Population = population, Code = code, Group = group, Events = sum(y),
      Non_events = sum(y == 0L), Valid_N = length(y), Expected_events = expected[i, "Yes"],
      Expected_non_events = expected[i, "No"], Event_hospitals = length(unique(h[y == 1L])),
      Non_event_hospitals = length(unique(h[y == 0L])), Design_SE = as.numeric(survey::SE(mean_fit)))
  }))
  hospitals <- sort(unique(v$HospitalID))
  loo <- do.call(rbind, lapply(hospitals, function(h) {
    mask <- v$HospitalID != h
    dropped <- subset(d, mask)
    vv <- dropped$variables
    tab <- table(vv$TestGroup, vv[[paste0("Outcome_", code)]])
    ex <- outer(rowSums(tab), colSums(tab))/sum(tab)
    omitted <- fit_association(dropped, code)
    s <- extract_test(omitted)
    data.frame(Population = population, Code = code,
      Omission = sprintf("Cluster %02d", match(h, levels(z$HospitalID))),
      design_information(dropped), Events = sum(vv[[code]]), Min_cell = min(tab),
      Min_expected = min(ex), as.list(s), Warnings = omitted$warnings, Error = omitted$error)
  }))
  sparse <- min(expected) < 5
  tests <- data.frame(Population = population, Code = code, Indicator = endpoint_labels[[code]],
    Family = if (code %in% c("Death", "Composite")) "Clinical outcome" else "Treatment process",
    design_information(d), Total_events = sum(v[[code]]), Min_cell = min(xt),
    Min_expected = min(expected), as.list(statistics), Sparse = sparse,
    Warnings = model$warnings, Error = model$error,
    Min_event_hospitals = min(cells$Event_hospitals), Min_non_event_hospitals = min(cells$Non_event_hospitals),
    Leave_one_min_P = min(loo$P, na.rm = TRUE), Leave_one_max_P = max(loo$P, na.rm = TRUE),
    Leave_one_warned = sum(nzchar(loo$Warnings)), Leave_one_failed = sum(!is.finite(loo$P)),
    Singleton_method = "Full-design lonely.psu=fail; domain adjustment off",
    Display_P = paste0(format_p(statistics[["P"]]), if (sparse) "†" else ""))
  list(tests = tests, cells = cells, loo = loo)
}
effect_size <- function(full_design) {
  z <- full_design$variables
  include <- z$Eligible_Enterobacterales & z$Complete_AST & z$OK_Carb7 & z$TestGroup %in% patterns[c(1L, 3L)]
  d <- subset(full_design, include)
  estimates <- survey::svyby(~Carb7, ~TestGroup, d, survey::svymean, covmat = TRUE)
  contrast <- survey::svycontrast(estimates, setNames(c(-1, 1), patterns[c(1L, 3L)]))
  point <- as.numeric(coef(contrast))
  se <- as.numeric(survey::SE(contrast))
  df <- survey::degf(d)
  ci <- point + c(-1, 1) * qt(.975, df) * se
  v <- d$variables
  e <- tapply(v$Carb7, v$TestGroup, sum)
  n <- table(v$TestGroup)
  data.frame(Population = "Enterobacterales", Code = "Carb7", Contrast = "Both R minus Both S",
    Both_S_events = unname(e[patterns[1L]]), Both_S_N = unname(n[patterns[1L]]),
    Both_R_events = unname(e[patterns[3L]]), Both_R_N = unname(n[patterns[3L]]),
    Both_S_percent = 100 * coef(estimates)[patterns[1L]], Both_R_percent = 100 * coef(estimates)[patterns[3L]],
    Difference_pp = 100 * point, Design_SE_pp = 100 * se,
    CI_95_lower_pp = 100 * ci[1L], CI_95_upper_pp = 100 * ci[2L], design_information(d),
    Interval = "Nominal unadjusted Student t interval; two-group domain design df")
}
cohort_summary <- function(dat, design, tests) {
  z <- dat$infants
  ent <- z$Eligible_Enterobacterales
  ec <- z$Eligible_Ecoli
  value <- c(Infants = nrow(z), Isolates = nrow(dat$isolates), Hospitals = nlevels(z$HospitalID),
    Macroregions = nlevels(z$Region), First_year = min(z$Year), Last_year = max(z$Year),
    Dual_pathogen_infants = sum(z$Isolates == 2L), Full_design_df = survey::degf(design),
    Initial_3GC_containing = sum(z$Initial %in% c("3GC", "3GC+PEN")),
    Initial_3GC_PEN = sum(z$Initial == "3GC+PEN"), Guideline_specific = sum(z$Guideline),
    Initial_BSP = sum(z$Initial == "BSP"),
    Initial_Carbapenem_containing = sum(z$Initial %in% c("CARB", "CARB+GLY/LZD")), Initial_CARB = sum(z$Initial == "CARB"),
    Other_known_regimens = sum(z$Initial %in% c("2GC", "PEN", "1GC", "GLY/LZD", "Other")),
    Unknown_initial_regimen = sum(z$Initial == "UNK"),
    Enterobacterales_eligible = sum(ent), Enterobacterales_complete = sum(ent & z$Complete_AST),
    Enterobacterales_excluded_I = sum(ent & z$AST_exclusion %in% "Contains I"),
    Enterobacterales_excluded_unavailable = sum(ent & z$AST_exclusion %in% "Unavailable or uninterpretable paired AST"),
    Ecoli_eligible = sum(ec), Ecoli_complete = sum(ec & z$Complete_AST),
    Ecoli_excluded_I = sum(ec & z$AST_exclusion %in% "Contains I"),
    Ecoli_excluded_unavailable = sum(ec & z$AST_exclusion %in% "Unavailable or uninterpretable paired AST"),
    Polymicrobial_excluded_from_AST = sum(z$Initial %in% c("3GC", "3GC+PEN") & z$Isolates > 1L &
      z$ID %in% dat$isolates$ID[dat$isolates$Taxon %in% dat$rules$Taxon]),
    Enterobacterales_death_min_expected = tests$Min_expected[tests$Population == "Enterobacterales" & tests$Code == "Death"],
    Ecoli_death_min_expected = tests$Min_expected[tests$Population == "E. coli" & tests$Code == "Death"])
  data.frame(Item = names(value), Value = unname(value))
}
submission_values <- function(descriptive, tests) {
  rows <- split(descriptive, factor(descriptive$Row, levels = unique(descriptive$Row)))
  result <- do.call(rbind, lapply(rows, function(x) {
    row <- x[1L, c("Section", "Population", "Group", "Infants")]
    for (k in names(endpoint_labels)) row[[k]] <- x$Display[x$Code == k]
    row
  }))
  for (population in names(eligible_field)) {
    p <- tests[tests$Population == population, ]
    row <- data.frame(Section = "Susceptibility", Population = population, Group = "Overall P", Infants = "—")
    for (k in names(endpoint_labels)) row[[k]] <- p$Display_P[p$Code == k]
    last <- max(which(result$Population == population))
    result <- rbind(result[seq_len(last), ], row, result[seq_len(nrow(result)) > last, ])
  }
  rownames(result) <- NULL
  result
}
analyse_s18 <- function(dat, design) {
  z <- dat$infants
  descriptive <- descriptive_results(z)
  fits <- lapply(names(eligible_field), function(p) lapply(names(endpoint_labels), function(k) test_domain(design, p, k)))
  fits <- unlist(fits, recursive = FALSE)
  tests <- do.call(rbind, lapply(fits, `[[`, "tests"))
  tests$BH_q <- NA_real_
  tests$Family_size <- NA_integer_
  for (family in unique(tests$Family)) {
    use <- tests$Family == family
    tests$BH_q[use] <- p.adjust(tests$P[use], method = "BH")
    tests$Family_size[use] <- sum(use)
  }
  bh <- do.call(rbind, lapply(unique(tests$Family), function(family) {
    x <- tests[tests$Family == family, ]
    data.frame(Family = family, Tests = nrow(x), Minimum_q = min(x$BH_q), Below_0_05 = sum(x$BH_q < .05))
  }))
  directional <- do.call(rbind, lapply(names(eligible_field), function(p) {
    ids <- z$ID[z[[eligible_field[[p]]]]]
    i <- match(ids, dat$isolates$ID)
    a <- dat$isolates$CAZ[i]
    b <- dat$isolates$CTX[i]
    data.frame(Population = p, CAZ_S_CTX_R = sum(a %in% "S" & b %in% "R"),
      CAZ_R_CTX_S = sum(a %in% "R" & b %in% "S"))
  }))
  states <- do.call(rbind, lapply(names(endpoint_labels), function(k) {
    x <- table(z[[paste0("State_", k)]])
    data.frame(Code = k, State = names(x), Infants = as.integer(x))
  }))
  display <- submission_values(descriptive, tests)
  list(S18_descriptive_results = descriptive,
    S18_descriptive_rows = display[display$Group != "Overall P", ],
    S18_overall_tests = tests, S18_BH_corrections = tests[c("Population", "Code", "Family", "Family_size", "P", "BH_q")],
    S18_BH_summary = bh, S18_effect_size = effect_size(design),
    S18_3x2_cells = do.call(rbind, lapply(fits, `[[`, "cells")),
    S18_hospital_deletion = do.call(rbind, lapply(fits, `[[`, "loo")),
    S18_cohort_and_notes = cohort_summary(dat, design, tests),
    S18_directional_AST = directional, S18_evaluability = states,
    S18_submission_values = display)
}

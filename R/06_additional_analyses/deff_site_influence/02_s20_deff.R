# Sources: original CANEOS_Achieved_DEFF_S3.R estimate(); accepted 11.16
# extension's precision, regional domains and hospital influence analyses.
infant_design <- function(data) survey::svydesign(ids = ~hospital_id, strata = ~macroregion,
  weights = ~w_unit, data = data, nest = TRUE)

deff_estimate <- function(data, group) {
  y <- data$target
  wgt <- data$w_unit
  n <- length(y)
  design <- infant_design(data)
  psus <- nrow(unique(data[c("hospital_id", "macroregion")]))
  df <- survey::degf(design)
  check(df == psus - length(unique(data$macroregion)), paste("Unexpected S21 design df:", group))
  fit <- survey::svymean(~target, design, deff = "replace")
  p_hat <- as.numeric(coef(fit))
  v <- as.numeric(vcov(fit))
  deff <- as.numeric(survey::deff(fit))
  p_manual <- sum(wgt * y) / sum(wgt)
  v_srs <- sum(wgt * (y - p_manual)^2) / sum(wgt) / (n - 1)
  hospital <- aggregate(wgt * (y - p_manual) / sum(wgt),
    data[c("macroregion", "hospital_id")], sum)
  v_manual <- sum(vapply(split(hospital$x, hospital$macroregion, drop = TRUE),
    function(u) length(u) * var(u), numeric(1)))
  deff_manual <- v_manual / v_srs
  check(is.finite(deff) && deff > 0 && df > 0 && near(p_hat, p_manual) &&
    near(v, v_manual) && near(deff, deff_manual) && near(n / deff, n * v_srs / v),
    paste("S21 independent variance check failed:", group))
  se <- sqrt(v)
  half <- qt(0.975, df) * se
  data.frame(Region = group, Hospitals = psus, N = n, Positive = sum(y), Proportion = p_hat,
    DesignVariance = v, SRSVariance = v_srs, DEFF = deff, DEFF_manual = deff_manual,
    DesignVariance_manual = v_manual, EffectiveN = n / deff, DF = df, SE = se,
    Lower95 = p_hat - half, Upper95 = p_hat + half, HalfWidth95 = half)
}

hospital_influence <- function(data, group, sites) {
  baseline <- deff_estimate(data, group)
  one <- function(h) {
    removed <- data$hospital_id == h
    est <- if (h == "None") baseline else deff_estimate(data[!removed, ], group)
    if (h != "None") check(sum(removed) > 0 && est$Hospitals == baseline$Hospitals - 1L &&
      est$DF == baseline$DF - 1L, paste("S21 deletion did not remove one PSU:", h))
    data.frame(Omitted_PSU = h, Removed_n = sum(removed), Removed_positive = sum(data$target[removed]),
      est, Delta_proportion_pp = 100 * (est$Proportion - baseline$Proportion))
  }
  bind(lapply(c("None", sites), one))
}

compute_s21 <- function(infant) {
  central <- infant[infant$macroregion == "Central", ]
  group_rows <- function(g) if (g == "Overall") infant else infant[infant$macroregion == g, ]
  a <- bind(lapply(c("Overall", regions), function(g) deff_estimate(group_rows(g), g)))
  margin_pp <- 7.32
  a$Precision_target_met <- ifelse(a$Region == "Overall", "N/A",
    ifelse(100 * a$HalfWidth95 <= margin_pp, "Yes", "No"))
  check(identical(a$Precision_target_met, c("N/A","Yes","Yes","No","Yes")), "Unexpected S21 margin decisions.")
  full <- infant_design(infant)
  domain_qc <- bind(lapply(regions, function(g) {
    s <- subset(full, macroregion == g)
    f <- survey::svymean(~target, s, deff = "replace")
    p <- a[a$Region == g, ]
    n <- nrow(s$variables)
    v <- as.numeric(vcov(f))
    d <- as.numeric(survey::deff(f))
    srs <- as.numeric(survey::svyvar(~target, s)) / n
    observed <- c(coef(f), v, srs, d, survey::degf(s), n, n / d)
    expected <- c(p$Proportion, p$DesignVariance, p$SRSVariance, p$DEFF, p$DF, p$N, p$EffectiveN)
    check(near(observed, expected), paste("S21 regional domain verification failed:", g))
    data.frame(Region = g, N = n, DF = survey::degf(s), Domain_proportion = as.numeric(coef(f)),
      Domain_variance = v, Domain_SRS_WR = srs, Domain_DEFF = d, Domain_effective_n = n / d)
  }))
  site_num <- sort(unique(as.integer(sub("^Site ", "", central$hospital_id))))
  b <- rbind(hospital_influence(central, "Central", paste("Site", site_num)),
    hospital_influence(infant, "Overall", "Site 23"))
  check(nrow(b) == 10L && all(b$N[b$Omitted_PSU == "Site 23"] == c(215L, 1103L)) &&
    all(b$Positive[b$Omitted_PSU == "Site 23"] == c(45L, 234L)) &&
    100 * b$HalfWidth95[b$Omitted_PSU == "Site 23" & b$Region == "Central"] > margin_pp,
    "Unexpected Site 23 exclusion results.")
  a$Panel <- "A"
  b$Panel <- "B"
  a$Row <- a$Region
  b$Row <- paste(b$Region, ifelse(b$Omitted_PSU == "None", "all hospitals",
    paste("excluding", b$Omitted_PSU)), sep = ", ")
  for (f in setdiff(names(b), names(a))) a[[f]] <- NA
  for (f in setdiff(names(a), names(b))) b[[f]] <- NA
  results <- rbind(a, b)[c("Panel", "Row", setdiff(names(a), c("Panel","Row")))]
  results$Method <- "Taylor linearization, Student t; with-replacement SRS reference"
  list(results = results, domains = domain_qc)
}


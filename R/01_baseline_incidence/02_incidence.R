# Aggregate the same hospitals across all years and retain zero-case hospital-years.
hospital_year <- do.call(rbind, lapply(seq_len(nrow(site_map)), function(i) {
  s <- site_map[i, ]
  a <- inf[inf$HospitalID == s$HospitalID, ]
  counts <- function(years) as.integer(table(factor(years, levels = year_set)))
  data.frame(HospitalID = s$HospitalID, Region = s$Region, Year = year_set, t = year_set - min(year_set), w = 1,
    Admissions = as.integer(den[i, ]), EOS = counts(a$Year),
    Ecoli = counts(a$Year[a$Ecoli == 1L]), GBS = counts(a$Year[a$GBS == 1L]),
    Isolates = counts(iso$Year[iso$HospitalID == s$HospitalID]))
}))
require_data(sum(hospital_year$EOS) == nrow(inf) && sum(hospital_year$Isolates) == nrow(iso) &&
               sum(hospital_year$Admissions) == sum(den),
             "Hospital-year counts failed aggregation checks.")
hospital_total <- aggregate(hospital_year[c("Admissions", "EOS", "Ecoli", "GBS")],
                            hospital_year[c("HospitalID", "Region")], sum)
hospital_total <- hospital_total[match(site_map$HospitalID, hospital_total$HospitalID), ]
hospital_total$w <- 1
hospital_year$RatePer1000 <- 1000 * hospital_year$EOS / hospital_year$Admissions
write_result(hospital_year[c("HospitalID", "Region", "Year", "Admissions", "EOS", "Ecoli", "GBS",
                            "Isolates", "RatePer1000")], "Hospital_Year")
pooled_rate <- function(d, endpoint) {
  f <- reformulate(endpoint)
  fit <- survey::svyratio(f, ~Admissions, design(d))
  point <- as.numeric(coef(fit))
  se <- as.numeric(survey::SE(fit))
  scores <- (d[[endpoint]] - point * d$Admissions) / sum(d$Admissions)
  independent <- sqrt(sum(vapply(split(scores, d$Region, drop = TRUE),
                                function(x) length(x) * var(x), numeric(1))))
  if (!isTRUE(all.equal(se, independent, tolerance = 1e-12))) stop("Hospital ratio variance failed independent check.")
  df <- design_df(d)
  if (df <= 0 || point <= 0 || !is.finite(se)) stop("Ratio interval not estimable: ", endpoint)
  margin <- qt(0.975, df) * se / point
  data.frame(Cases = sum(d[[endpoint]]), Admissions = sum(d$Admissions),
    Rate = 1000 * point, Lower = 1000 * point * exp(-margin), Upper = 1000 * point * exp(margin),
    SE = 1000 * se, Hospitals = psu_count(d), DF = df,
    Method = "Hospital-stratified ratio, Taylor log-delta t CI")
}
incidence_comparison <- function() {
  d <- hospital_total
  fit <- survey::svyglm(EOS ~ Region + offset(log(Admissions)), design(d), family = quasipoisson())
  tst <- survey::regTermTest(fit, ~Region, df = design_df(d), method = "Wald")
  data.frame(Method = "Hospital-level quasi-Poisson, survey Wald F", Statistic = as.numeric(tst$Ftest),
    DF1 = as.numeric(tst$df), DF2 = as.numeric(tst$ddf), P = as.numeric(tst$p),
    Hospitals = nrow(d), N = nrow(d), Status = "Estimated",
    Diagnostic = "Hospital totals, log admissions offset, unit weights, region strata")
}
site_rate <- function(d, endpoint) {
  n <- d[[endpoint]]
  N <- d$Admissions
  lim <- poisson.test(n, T = N)$conf.int * 1000
  data.frame(Cases = n, Admissions = N, Rate = 1000 * n / N, Lower = lim[1L], Upper = lim[2L],
             SE = NA_real_, Hospitals = 1L, DF = NA_real_, Method = "Exact Poisson CI")
}
s4_rate <- function(d, region, site = "") {
  endpoints <- c("EOS", "Ecoli", "GBS")
  do.call(rbind, lapply(endpoints, function(endpoint) {
    rate <- if (nzchar(site)) site_rate(d, endpoint) else pooled_rate(d, endpoint)
    cbind(Region = region, Site = site, Endpoint = endpoint, rate)
  }))
}
s4_display <- function(x, region_label, site_label) {
  cells <- rbind(paste(x$Cases, x$Admissions, sep = "/"), fmt_interval(x$Rate, x$Lower, x$Upper))
  names_table <- c("Stratum", "Centre", "EOS n/N", "EOS per 1000 (95% CI)",
                  "E. coli n/N", "E. coli per 1000 (95% CI)", "GBS n/N", "GBS per 1000 (95% CI)")
  as.data.frame(setNames(as.list(c(region_label, site_label, as.vector(cells))), names_table),
                check.names = FALSE, stringsAsFactors = FALSE)
}
overall <- s4_rate(hospital_total, "Overall")
s4_models <- list(overall)
s4_rows <- list(s4_display(overall, "Overall (all centres)", ""))
for (g in region_names) {
  d <- hospital_total[region_rows(hospital_total, g), ]
  for (i in seq_len(nrow(d))) {
    single <- s4_rate(d[i, ], g, d$HospitalID[i])
    s4_models[[length(s4_models) + 1L]] <- single
    s4_rows[[length(s4_rows) + 1L]] <- s4_display(single, if (i == 1L) g else "", d$HospitalID[i])
  }
  subtotal <- s4_rate(d, g)
  s4_models[[length(s4_models) + 1L]] <- subtotal
  s4_rows[[length(s4_rows) + 1L]] <- s4_display(subtotal, "", "Subtotal")
}
s4_numeric <- do.call(rbind, s4_models)
pooled_eos <- s4_numeric[s4_numeric$Site == "" & s4_numeric$Endpoint == "EOS", ]
inc_region <- incidence_comparison()
write_result(do.call(rbind, s4_rows), "S4_Incidence")
write_result(s4_numeric, "S4_Incidence_Estimates")
write_result(inc_region, "Incidence_Regional_Comparison")

# The pooled trend includes macro-region main effects, as in the revised paper.
run_inc_trend <- function(g) {
  d <- hospital_year[region_rows(hospital_year, g), ]
  overall <- g == "Overall"
  fm <- if (overall) EOS ~ t + Region + offset(log(Admissions)) else EOS ~ t + offset(log(Admissions))
  method <- if (overall) "Quasi-Poisson, year + macro-region + offset" else "Quasi-Poisson, year + offset"
  a <- trend_skeleton(d, method, "IRR", if (overall) 4L else 1L)
  a$Events <- sum(d$EOS)
  a$EventHospitals <- length(unique(d$HospitalID[d$EOS > 0]))
  fit <- capture(survey::svyglm(fm, design(d), family = quasipoisson(),
                               control = stats::glm.control(epsilon = 1e-11, maxit = 50)))
  a$Converged <- !is.null(fit$Fit) && isTRUE(fit$Fit$converged)
  a$Diagnostic <- fit$Message
  if (a$Converged) {
    a <- fill_trend(a, coef(fit$Fit)[["t"]], sqrt(vcov(fit$Fit)["t", "t"]))
  }
  cbind(Region = g, Admissions = sum(d$Admissions), ZeroCaseHospitalYears = sum(d$EOS == 0), Statistic = a$Beta/a$SE, a)
}
incidence_trends <- do.call(rbind, lapply(regions, run_inc_trend))
write_result(incidence_trends, "Incidence_Trends")
region_cohort <- do.call(rbind, lapply(regions, function(g) data.frame(Region = g,
  Infants = sum(region_rows(inf, g)), Isolates = sum(region_rows(iso, g)),
  Hospitals = psu_count(hospital_year[region_rows(hospital_year, g), ]))))
write_result(region_cohort, "Cohort")


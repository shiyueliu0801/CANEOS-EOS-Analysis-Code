# Taylor estimators and Rao–Scott tests from the 11.16 coverage, pathogen and S19 releases.
near <- function(a, b, tolerance = 1e-10) {
  isTRUE(all.equal(unname(a), unname(b), tolerance = tolerance, check.attributes = FALSE))
}
bind <- function(x) do.call(rbind, unname(x))

make_design <- function(data) {
  survey::svydesign(ids = ~Site, strata = ~Region, weights = ~Weight, data = data, nest = TRUE)
}
isolate_design <- make_design(iso)
infant_design <- make_design(infants)
stopifnot(survey::degf(isolate_design) == cohort$Hospitals[1L] - cohort$Strata[1L],
  survey::degf(isolate_design) == survey::degf(infant_design))

domain <- function(design, keep, ...) {
  stopifnot(length(keep) == nrow(design$variables), !anyNA(keep))
  added <- list(...)
  for (label in names(added)) design$variables[[label]] <- added[[label]]
  subset(design, keep)
}

psu_variance <- function(design, y, keep, estimate) {
  # Include zero scores for hospitals without evaluable domain observations.
  original <- design$variables
  w <- weights(design)
  influence <- numeric(length(keep))
  influence[keep] <- w[keep]*(y[keep] - estimate)/sum(w[keep])
  score <- rowsum(influence, original$Site, reorder = FALSE)
  strata <- original$Region[match(as.integer(rownames(score)), original$Site)]
  sum(vapply(split(score[,1L], strata), function(x) length(x)*var(x), 0))
}

proportion <- function(design, y, keep, clinical = FALSE) {
  domain_rows <- keep & !is.na(y)
  d <- domain(design, domain_rows, Response = y)
  n <- sum(domain_rows)
  e <- sum(y[domain_rows])
  fit <- survey::svymean(~Response, d)
  estimate <- if (clinical) e/n else weighted.mean(y[domain_rows], weights(d))
  se <- as.numeric(survey::SE(fit))
  df <- survey::degf(d)
  stopifnot(n > 0L, near(coef(fit), estimate), near(se^2, psu_variance(design, y, domain_rows, estimate)))
  interval <- c(NA_real_, NA_real_)
  if (estimate > 0 && estimate < 1 && df > 0 && se > 0) {
    interval <- plogis(qlogis(estimate) + c(-1, 1)*qt(.975, df)*se/(estimate*(1-estimate)))
    if (clinical) {
      original <- confint(survey::svyciprop(~Response, d, method = "xlogit", df = df))
      stopifnot(near(as.numeric(original), interval))
      interval <- as.numeric(original)
    }
  }
  data.frame(Events = e, N = n, Estimate = estimate, Lower = interval[1L], Upper = interval[2L],
    SE = se, DF = df, Candidate = sum(keep), Hospitals = length(unique(d$variables$Site)))
}

wide <- function(groups) {
  do.call(cbind, lapply(names(groups), function(g) {
    d <- groups[[g]]
    names(d) <- paste(g, names(d), sep = "_")
    d
  }))
}

rao_scott <- function(y) {
  d <- domain(infant_design, !is.na(y), State = factor(y),
    Group = factor(infants$Group, levels = c("CoNS", "Non-CoNS")))
  d$variables$State <- droplevels(d$variables$State)
  tab <- table(d$variables$Group, d$variables$State)
  stopifnot(all(dim(tab) == c(2L, 2L)), all(rowSums(tab) > 0L))
  fit <- survey::svychisq(~State + Group, d, statistic = "F")
  ans <- data.frame(Test = "Rao–Scott second-order F", Test_N = nrow(d$variables),
    Design_DF = survey::degf(d), F = unname(fit$statistic), DF1 = unname(fit$parameter[1L]),
    DF2 = unname(fit$parameter[2L]), P = fit$p.value, Min_cell = min(tab),
    Min_expected_unweighted = min(outer(rowSums(tab), colSums(tab))/sum(tab)))
  stopifnot(all(is.finite(unlist(ans[c("F","DF1","DF2","P")]))),
    near(ans$P, pf(ans$F, ans$DF1, ans$DF2, lower.tail = FALSE)))
  ans
}

adjust_bh <- function(x, family) {
  x$Family <- family
  x$Planned_tests <- ave(rep(1L, nrow(x)), family, FUN = sum)
  x$Estimated_tests <- ave(as.integer(is.finite(x$P)), family, FUN = sum)
  x$q <- ave(x$P, family, FUN = function(p) p.adjust(p, "BH", n = length(p)))
  x
}

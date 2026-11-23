# Survey design, clustered covariance, marginal risks and trend diagnostics.

capture <- function(expr) {
  messages <- character()
  value <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    messages <<- c(messages, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) { messages <<- c(messages, conditionMessage(e)); NULL })
  list(Fit = value, Message = paste(unique(messages), collapse = " | "))
}
combine_text <- function(...) paste(Filter(nzchar, list(...)), collapse = " | ")
design <- function(d) survey::svydesign(ids = ~HospitalID, strata = ~Region, weights = ~w, nest = TRUE, data = d)
psu_count <- function(d) length(unique(d$HospitalID))
design_df <- function(d) psu_count(d) - length(unique(d$Region))
region_rows <- function(d, g) if (g == "Overall") rep(TRUE, nrow(d)) else d$Region == g
fmt_p <- function(p) ifelse(is.na(p), "—", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
fmt_dec <- function(x, n = 2) ifelse(is.na(x), "—", sprintf(paste0("%.", n, "f"), ifelse(abs(x) < 0.5 * 10^-n, 0, x)))
fmt_interval <- function(x, l, u, n = 2) ifelse(is.na(x), "—",
  ifelse(is.na(l) | is.na(u), paste0(fmt_dec(x, n), " (—)"),
         paste0(fmt_dec(x, n), " (", fmt_dec(l, n), "–", fmt_dec(u, n), ")")))
controls <- glm.control(epsilon = 1e-14, maxit = 50)

manual_vcov <- function(fit, d) {
  x <- model.matrix(fit)
  pr <- fitted(fit)
  score <- x * as.numeric(d$y - pr)
  cluster <- rowsum(score, d$HospitalID)
  stratum <- d$Region[match(rownames(cluster), d$HospitalID)]
  meat <- matrix(0, ncol(x), ncol(x))
  for (m in split(seq_len(nrow(cluster)), stratum, drop = TRUE)) {
    stopifnot(length(m) > 1L)
    centered <- scale(cluster[m, , drop = FALSE], center = TRUE, scale = FALSE)
    meat <- meat + length(m) / (length(m) - 1L) * crossprod(centered)
  }
  bread <- solve(crossprod(x, x * as.numeric(pr * (1 - pr))))
  bread %*% meat %*% bread
}

# Binary logistic model with the working-independence score and unit weights.
fit_binary <- function(d, by, levels, id, adjust = FALSE) {
  d$Group <- factor(d[[by]], levels = levels)
  d$y <- d[[id]]
  d <- d[!is.na(d$Group) & !is.na(d$y) & (!adjust | !is.na(d$GAWeeks)), ]
  a <- data.frame(Outcome = labels[[id]], Code = id, Comparison = by, N = nrow(d), Events = sum(d$y),
    Hospitals = psu_count(d), EventHospitals = length(unique(d$HospitalID[d$y == 1L])),
    NonEventHospitals = length(unique(d$HospitalID[d$y == 0L])), DesignDF = design_df(d),
    ResidualDF = design_df(d) - length(levels) + 1L - adjust, NumDF = length(levels) - 1L,
    Method = "Working-independence marginal logistic; Taylor hospital sandwich",
    CovarianceMaxError = NA_real_, CoefficientMaxError = NA_real_, F = NA_real_,
    P = NA_real_, Converged = FALSE, Diagnostic = "")
  support <- table(d$Group, factor(d$y, levels = 0:1))
  by_stratum <- table(d$Region, d$HospitalID)
  if (any(rowSums(by_stratum > 0) < 2) || a$ResidualDF <= 0L) {
    a$Diagnostic <- "Insufficient hospital design degrees of freedom."
    return(list(audit = a, data = d))
  }
  if (any(support == 0)) {
    a$Diagnostic <- "At least one comparison group has no events or no non-events."
    return(list(audit = a, data = d))
  }
  f <- if (adjust) y ~ Group + GAWeeks else y ~ Group
  dd <- design(d)
  stopifnot(survey::degf(dd) == a$DesignDF)
  z <- capture(survey::svyglm(f, dd, family = quasibinomial(), control = controls))
  fit <- z$Fit
  if (is.null(fit) || !isTRUE(fit$converged) || nzchar(z$Message) || any(!is.finite(coef(fit)))) {
    a$Diagnostic <- combine_text("Model did not pass convergence or rank checks.", z$Message)
    return(list(audit = a, data = d))
  }
  old_fit <- glm(f, data = d, family = binomial(), control = controls)
  a$CoefficientMaxError <- max(abs(coef(fit) - coef(old_fit)))
  a$CovarianceMaxError <- max(abs(vcov(fit) - manual_vcov(fit, d)))
  if (a$CoefficientMaxError > 1e-6 || a$CovarianceMaxError > 1e-6) {
    stop("Hospital sandwich QC failed for ", id, " / ", by, " | beta error ",
         a$CoefficientMaxError, " | covariance error ", a$CovarianceMaxError)
  }
  test <- survey::regTermTest(fit, ~Group, df = a$ResidualDF, method = "Wald")
  ix <- grep("^Group", names(coef(fit)))
  b <- coef(fit)[ix]
  v <- vcov(fit)[ix, ix, drop = FALSE]
  f_manual <- drop(crossprod(b, solve(v, b))) / length(ix)
  a$F <- as.numeric(test$Ftest)
  a$P <- as.numeric(test$p)
  a$Converged <- TRUE
  stopifnot(length(ix) == a$NumDF, isTRUE(all.equal(a$F, f_manual)),
            isTRUE(all.equal(a$P, pf(f_manual, length(ix), a$ResidualDF, lower.tail = FALSE))))
  list(audit = a, fit = fit, data = d)
}
estimate_ci <- function(value, gradient, covariance, df, link = "identity") {
  se <- sqrt(drop(crossprod(gradient, covariance %*% gradient)))
  if (!is.finite(se) || df <= 0L) stop("Non-finite delta SE or insufficient degrees of freedom.")
  span <- qt(0.975, df) * se
  transform <- switch(link, identity = identity, log = exp, logit = plogis)
  setNames(transform(c(value, value - span, value + span)), c("Estimate", "Lower", "Upper"))
}
delta_validate <- function(f, b, g) {
  h <- 1e-5
  numeric <- vapply(seq_along(b), function(j) {
    a <- b; z <- b
    a[j] <- a[j] + h; z[j] <- z[j] - h
    (f(a) - f(z)) / (2 * h)
  }, 0)
  if (max(abs(numeric - g)) > 1e-5 * max(1, max(abs(g)))) stop("Gradient verification failed.")
}

# Standardise to the same outcome-specific model population for every group.
# Risk differences and ratios use the full coefficient covariance.
marginal <- function(model) {
  if (!isTRUE(model$audit$Converged)) stop("Required marginal model did not converge: ", model$audit$Code)
  f <- model$fit; d <- model$data; df <- model$audit$ResidualDF
  b <- coef(f); v <- vcov(f); lv <- levels(d$Group)
  mx <- lapply(lv, function(g) {
    new <- transform(d, Group = factor(g, levels = lv))
    x <- model.matrix(delete.response(terms(f)), new, contrasts.arg = f$contrasts)
    stopifnot(identical(colnames(x), names(b)))
    x
  })
  mean_risk <- function(x, z) mean(plogis(drop(x %*% z)))
  risk <- vapply(mx, mean_risk, 0, z = b)
  grad <- lapply(mx, function(x) {
    p <- plogis(drop(x %*% b))
    colMeans(x * (p * (1 - p)))
  })
  out_rows <- lapply(seq_along(lv), function(i) {
    prob <- risk[i]; ref <- risk[1]; g <- grad[[i]]; g0 <- grad[[1]]
    delta_validate(function(z) qlogis(mean_risk(mx[[i]], z)), b, g / (prob * (1 - prob)))
    ci <- estimate_ci(qlogis(prob), g / (prob * (1 - prob)), v, df, "logit")
    rd_g <- 100 * (g - g0)
    rr_g <- g / prob - g0 / ref
    delta_validate(function(z) 100 * (mean_risk(mx[[i]], z) - mean_risk(mx[[1]], z)), b, rd_g)
    delta_validate(function(z) log(mean_risk(mx[[i]], z) / mean_risk(mx[[1]], z)), b, rr_g)
    rd <- estimate_ci(100 * (prob - ref), rd_g, v, df)
    rr <- estimate_ci(log(prob / ref), rr_g, v, df, "log")
    e <- numeric(length(b))
    if (i > 1L) e[grep("^Group", names(b))[i - 1L]] <- 1
    or <- estimate_ci(sum(b * e), e, v, df, "log")
    r <- d$Group == lv[i]
    data.frame(Group = lv[i], Reference = lv[1], GroupN = sum(r), Events = sum(d$y[r]),
      ObservedRisk = mean(d$y[r]), StandardN = nrow(d), Risk = ci[1], RiskLower = ci[2], RiskUpper = ci[3],
      RD_pp = rd[1], RDLower_pp = rd[2], RDUpper_pp = rd[3], RR = rr[1], RRLower = rr[2], RRUpper = rr[3],
      OR = or[1], ORLower = or[2], ORUpper = or[3], OmnibusP = model$audit$P, DF = df, row.names = NULL)
  })
  do.call(rbind, out_rows)
}

base_model <- function(k, by, groups) {
  d <- inf
  d$Group <- factor(d[[by]], levels=groups)
  d$y <- d[[k]]
  d[!is.na(d$Group) & !is.na(d$y), ]
}

fit_ga_adjusted <- function(k, by, groups) {
  d <- base_model(k,by,groups)
  a <- data.frame(Outcome=labels[[k]],Code=k,Comparison=by,N=nrow(d),Events=sum(d$y),
    Hospitals=psu_count(d),EventHospitals=length(unique(d$HospitalID[d$y==1L])),
    NonEventHospitals=length(unique(d$HospitalID[d$y==0L])),DesignDF=design_df(d),
    ResidualDF=design_df(d)-(length(groups)-1L)-(length(region_names)-1L),NumDF=length(groups)-1L,
    Method="Working-independence marginal logistic; Taylor hospital sandwich",
    CovarianceMaxError=NA_real_,CoefficientMaxError=NA_real_,F=NA_real_,P=NA_real_,Converged=FALSE,Diagnostic="")
  if (any(table(d$Group, factor(d$y,levels=0:1))==0L) || a$ResidualDF<=0L) {
    a$Diagnostic <- "Zero comparison cell or insufficient design degrees of freedom."
    return(list(audit=a,data=d))
  }
  formula <- y ~ Group + Region
  dd <- design(d)
  z <- capture(survey::svyglm(formula,dd,family=quasibinomial(),control=controls))
  f <- z$Fit
  if(is.null(f) || !isTRUE(f$converged) || f$rank!=length(coef(f)) || nzchar(z$Message)) {
    a$Diagnostic <- combine_text("Convergence, rank or separation problem.",z$Message)
    return(list(audit=a,data=d))
  }
  simple <- glm(formula,data=d,family=binomial(),control=controls)
  a$CoefficientMaxError <- max(abs(coef(simple)-coef(f)))
  a$CovarianceMaxError <- max(abs(manual_vcov(f,d)-vcov(f)))
  stopifnot(a$CoefficientMaxError<1e-6,a$CovarianceMaxError<1e-6,
            survey::degf(dd)==a$DesignDF,f$df.residual==a$ResidualDF)
  j <- grep("^Group",names(coef(f)))
  q <- coef(f)[j]
  cov <- vcov(f)[j,j,drop=FALSE]
  direct <- drop(crossprod(q,solve(cov,q)))/length(j)
  test <- survey::regTermTest(f,~Group,df=a$ResidualDF,method="Wald")
  a$F <- unname(test$Ftest[1L])
  a$P <- unname(test$p[1L])
  a$Converged <- TRUE
  stopifnot(abs(a$F-direct)<1e-8,abs(a$P-pf(direct,a$NumDF,a$ResidualDF,lower.tail=FALSE))<1e-10)
  list(audit=a,fit=f,data=d)
}

support_by_group <- function(d, groups) do.call(rbind,lapply(groups,function(g) {
  z <- d[d$Group==g, ]
  data.frame(Group=g,N=nrow(z),Events=sum(z$y),NonEvents=sum(z$y==0L),
    EventHospitals=length(unique(z$HospitalID[z$y==1L])),
    NonEventHospitals=length(unique(z$HospitalID[z$y==0L])))
}))

sparse_rr <- function(model, effects, support) {
  g0 <- support[1L, ]
  j <- match(effects$Group,support$Group)
  g1 <- support[j, ]
  least_events <- pmin(g0$Events,g1$Events,g0$NonEvents,g1$NonEvents)
  least_hosp <- pmin(g0$EventHospitals,g1$EventHospitals,g0$NonEventHospitals,g1$NonEventHospitals)
  se <- (log(effects$RRUpper)-log(effects$RRLower))/(2*qt(.975,effects$DF))
  i <- seq_len(nrow(effects))
  relative <- i>1L
  bad_support <- relative & (least_events<5L | least_hosp<3L)
  very_precise <- relative & least_events<10L & is.finite(se) & se<0.10
  check <- which(relative & least_events<10L & !bad_support & is.finite(se) & se>0)
  max_shift <- rep(NA_real_, length(i))
  failures <- rep(NA_integer_, length(i))
  bad_delete <- rep(FALSE,length(i))
  if (length(check) && !is.null(model$fit)) {
    full <- model$data
    formula <- formula(model$fit)
    ref <- effects$Reference[1L]
    replace_group <- function(g) transform(full,Group=factor(g,levels=levels(full$Group)))
    ref_data <- replace_group(ref)
    alt <- lapply(effects$Group[check],replace_group)
    deleted <- vapply(unique(full$HospitalID),function(h) {
      z <- capture(glm(formula,data=full[full$HospitalID!=h, ],family=binomial(),control=controls))
      f <- z$Fit
      if (is.null(f) || !f$converged || f$boundary || nzchar(z$Message) || anyNA(coef(f))) {
        return(rep(NA_real_,length(check)))
      }
      p0 <- mean(predict(f,ref_data,type="response"))
      vapply(alt,function(new) log(mean(predict(f,new,type="response"))/p0),0)
    },numeric(length(check)))
    deleted <- matrix(deleted,nrow=length(check))
    for (k in seq_along(check)) {
      a <- check[k]
      r <- deleted[k, ]
      failures[a] <- sum(!is.finite(r))
      if (any(is.finite(r))) max_shift[a] <- max(abs(r[is.finite(r)]-log(effects$RR[a])))
      bad_delete[a] <- failures[a]>0L || isTRUE(max_shift[a]>2*se[a])
    }
  }
  issue <- vapply(i,function(a) {
    if (!relative[a]) return("")
    combine_text(
      if(!isTRUE(model$audit$Converged)) model$audit$Diagnostic else "",
      if(bad_support[a]) "Fewer than five events/non-events or three contributing hospitals in a comparison group." else "",
      if(very_precise[a]) "Sparse comparison with log-RR SE below 0.10." else "",
      if(bad_delete[a]) "Sparse comparison unstable after deleting one hospital." else "",
      if(!is.finite(se[a]) || se[a]<=0) "Invalid log-RR variance." else "")
  },"")
  data.frame(MinEventOrNonEvent=least_events,MinHospitalSupport=least_hosp,LogRR_SE=se,
    RRMaxDeletionShift=max_shift,RRDeletionFailures=failures,RRDisplayed=relative & !nzchar(issue),RRReason=issue)
}

same_beta <- function(d) {
  d <- d[!is.na(d$x), ]
  z <- capture(glm(x ~ t, data = d, family = binomial(), control = controls))
  if (is.null(z$Fit) || !isTRUE(z$Fit$converged) || grepl("0 or 1", z$Message, fixed = TRUE)) return(NA_real_)
  unname(coef(z$Fit)["t"])
}
trend_skeleton <- function(d) {
  data.frame(Method = "Working-independence marginal logistic; Taylor hospital sandwich",
    Effect = "OR", N = nrow(d), Hospitals = psu_count(d),
    DesignDF = design_df(d), DF = design_df(d) - 1L,
    FirstYear = if (nrow(d)) min(d$Year) else NA_real_, LastYear = if (nrow(d)) max(d$Year) else NA_real_,
    Years = length(unique(d$Year)), Events = NA_real_, NonEvents = NA_real_,
    EventHospitals = NA_real_, NonEventHospitals = NA_real_,
    Beta = NA_real_, SE = NA_real_, P = NA_real_, Annual = NA_real_, AnnualLower = NA_real_, AnnualUpper = NA_real_,
    Five = NA_real_, FiveLower = NA_real_, FiveUpper = NA_real_,
    SameDataOrdinaryBeta = NA_real_, CovarianceMaxError = NA_real_, Converged = NA, Status = "Not estimable",
    Diagnostic = "")
}
fill_trend <- function(a, b, se, ordinary = NA_real_) {
  if (length(b) != 1L || length(se) != 1L || !is.finite(b) || !is.finite(se) || se <= 0 || a$DF <= 0) {
    a$Diagnostic <- combine_text(a$Diagnostic, "Non-finite estimate, variance or denominator degrees of freedom.")
    return(a)
  }
  a$Beta <- b
  a$SE <- se
  a$P <- 2 * pt(-abs(b / se), a$DF)
  delta <- qt(0.975, a$DF) * se
  z <- exp(c(b, b - delta, b + delta, 5 * b, 5 * (b - delta), 5 * (b + delta)))
  if (any(!is.finite(z))) stop("Non-finite transformed result.")
  a[c("Annual", "AnnualLower", "AnnualUpper", "Five", "FiveLower", "FiveUpper")] <- as.list(z)
  a$SameDataOrdinaryBeta <- ordinary
  a$Status <- "Estimated"
  a
}
validate_trends <- function(t, label) {
  e <- t[t$Status == "Estimated", ]
  ok <- isTRUE(all.equal(e$Five, e$Annual^5)) &&
    isTRUE(all.equal(e$P, 2 * pt(-abs(e$Beta / e$SE), e$DF))) &&
    all(e$AnnualLower <= e$Annual & e$Annual <= e$AnnualUpper) &&
    all(e$FiveLower <= e$Five & e$Five <= e$FiveUpper) &&
    all(e$Hospitals >= 2 & e$DF > 0 & e$Converged)
  if (!ok) stop(label, ": inconsistent interval, trend scaling or variance.")
}
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
fit_trend <- function(d) {
  d <- d[!is.na(d$x), ]
  a <- trend_skeleton(d)
  a$Events <- sum(d$x == 1)
  a$NonEvents <- sum(d$x == 0)
  a$EventHospitals <- length(unique(d$HospitalID[d$x == 1]))
  a$NonEventHospitals <- length(unique(d$HospitalID[d$x == 0]))
  if (a$Years < 2L || a$DF <= 0 || a$N < 10L || length(unique(d$x)) < 2L) {
    a$Diagnostic <- "Insufficient variation, eligible observations or design degrees of freedom."
    return(a)
  }
  problem <- binary_problem(d)
  if (nzchar(problem)) {
    a$Diagnostic <- problem
    return(a)
  }
  d$y <- d$x
  z <- capture(survey::svyglm(y ~ t, design(d), family = quasibinomial(), control = controls))
  a$Diagnostic <- z$Message
  if (is.null(z$Fit)) return(a)
  a$Converged <- isTRUE(z$Fit$converged)
  if (!a$Converged || grepl("0 or 1|rank-deficient", z$Message, fixed = FALSE)) {
    a$Diagnostic <- combine_text(a$Diagnostic, "Convergence, separation or rank check failed.")
    return(a)
  }
  beta <- unname(coef(z$Fit)["t"])
  se <- sqrt(vcov(z$Fit)["t", "t"])
  if (!is.finite(beta) || !is.finite(se) || abs(beta) >= 5 || se >= 5) {
    a$Diagnostic <- combine_text(a$Diagnostic, "Unstable log odds or standard error.")
    return(a)
  }
  a$CovarianceMaxError <- max(abs(vcov(z$Fit) - manual_vcov(z$Fit, d)))
  fill_trend(a, beta, se, same_beta(d))
}

trend_cells <- function(a, omnibus = NA_real_) {
  if (a$Status != "Estimated") return(c(rep("—", 4), fmt_p(omnibus)))
  c(fmt_dec(a$Beta), fmt_dec(a$SE), fmt_interval(a$Annual, a$AnnualLower, a$AnnualUpper),
    fmt_interval(a$Five, a$FiveLower, a$FiveUpper), fmt_p(a$P))
}

trend_support <- function(d, a) {
  d <- d[!is.na(d$x), ]
  warnings <- character()
  get_glm <- function(z) {
    fit <- suppressWarnings(glm(x ~ t, family = binomial(), data = z, control = controls))
    if (!isTRUE(fit$converged) || fit$boundary ||
        any(!is.finite(coef(fit))) || length(unique(z$x)) < 2) return(NULL)
    fit
  }
  naive <- if (a$Status == "Estimated") get_glm(d) else NULL
  ordinary_se <- if (is.null(naive)) NA_real_ else sqrt(vcov(naive)["t", "t"])
  if (!is.null(naive)) {
    verify_diff <- abs(unname(coef(naive)["t"]) - a$Beta)
    stopifnot(verify_diff < 1e-7, a$CovarianceMaxError < 1e-6)
  }
  ratio <- a$SE / ordinary_se
  delete <- if (is.null(naive)) numeric() else vapply(unique(d$HospitalID), function(h) {
    z <- get_glm(d[d$HospitalID != h, ])
    if (is.null(z)) NA_real_ else unname(coef(z)["t"])
  }, numeric(1))
  max_shift <- if (length(delete) && any(is.finite(delete))) max(abs(delete - a$Beta), na.rm = TRUE) else NA_real_
  sep <- if (length(unique(d$x)) == 2L) {
    pos <- range(d$t[d$x == 1L]); neg <- range(d$t[d$x == 0L])
    pos[1] >= neg[2] || neg[1] >= pos[2]
  } else NA
  if (a$Status != "Estimated") warnings <- c(warnings, a$Diagnostic)
  if (min(a$Events, a$NonEvents) < 10L) warnings <- c(warnings, "Fewer than 10 events or non-events")
  if (min(a$EventHospitals, a$NonEventHospitals) < 4L)
    warnings <- c(warnings, "Events or non-events in fewer than four hospitals")
  # Check markedly small sandwich SEs against hospital-deletion influence on the same slope.
  if (is.finite(ratio) && ratio < 0.5 && (anyNA(delete) || isTRUE(max_shift > a$SE)))
    warnings <- c(warnings, "Sandwich SE below half the ordinary SE with an unstable hospital-deletion check")
  by_year <- table(factor(d$Year, levels = year_set), factor(d$x, levels = 0:1))
  data.frame(OrdinarySE = ordinary_se, SandwichToOrdinarySE = ratio,
    SeparationByYear = sep, MaxHospitalDeletionShift = max_shift,
    HospitalDeletionFailures = sum(is.na(delete)),
    EligibleByYear = paste(rowSums(by_year), collapse = "/"),
    EventsByYear = paste(by_year[, "1"], collapse = "/"),
    Reporting = if (length(warnings)) "Withheld" else "Report",
    ReportingReason = paste(unique(warnings[nzchar(warnings)]), collapse = ". "))
}

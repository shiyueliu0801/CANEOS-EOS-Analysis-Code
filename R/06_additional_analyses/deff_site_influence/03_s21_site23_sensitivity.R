# Sources: accepted 11.16 pathogen inputs_and_design.R and ast_and_esbl.R;
# CANEOS_Coverage_11.15.5.R (activity, combine, prop_fit, pair_fit);
# CANEOS_DEFF_S3_Extended_11.16.R (S22 scenarios and recalibration).

# AST denominators are specific to the reportable E. coli endpoint. The
# accepted Table 3 scope and intrinsic matrix have no exclusion for these four
# endpoints or their component agents. Both agents in each class are eligible.
derive_ecoli_ast <- function(raw, iso, fields) {
  sir <- c("敏感" = 1L, "中介" = 2L, "耐药" = 3L, S = 1L, I = 2L, R = 3L)
  mapping <- setNames(fields$Field, fields$Drug)
  observed <- function(agent) unname(sir[canon(raw[[mapping[[agent]]]])])
  # Rebuild class from observed applicable components: R > I > S; all missing stays missing.
  endpoints <- list(
    "Third-generation cephalosporins" = pmax(observed("Ceftazidime"), observed("Cefotaxime"), na.rm = TRUE),
    "Ampicillin" = observed("Ampicillin"), "Gentamicin" = observed("Gentamicin"),
    "Carbapenems" = pmax(observed("Imipenem"), observed("Meropenem"), na.rm = TRUE))
  endpoints <- lapply(endpoints, function(v) replace(v, iso$Organism != "E. coli", NA_integer_))
  # Study-specific reporting convention: a recorded positive flag defines ESBL
  # positivity among all E. coli tested routinely; no flag denotes nonpositive.
  list(values = endpoints, esbl_positive = as.integer(canon(raw[["ESBL"]]) %in% "是"))
}

activity <- function(x) {
  z <- matrix(NA_integer_, nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x))
  z[x %in% c("敏感", "S")] <- 1L
  z[x %in% c("耐药", "中介", "I", "R")] <- 0L
  z
}

combine <- function(a, b) {
  z <- rep(NA_integer_, length(a))
  z[a %in% 0L & b %in% 0L] <- 0L
  z[a %in% 1L | b %in% 1L] <- 1L
  z
}

regimens <- data.frame(Regimen = c("PEN + CAZ", "PEN + CTX", "AMP + CAZ", "AMP + CTX", "PEN + GEN", "AMP + GEN"),
  A = c("Penicillin", "Penicillin", "Ampicillin", "Ampicillin", "Penicillin", "Ampicillin"),
  B = c("Ceftazidime", "Cefotaxime", "Ceftazidime", "Cefotaxime", "Gentamicin", "Gentamicin"))
pairs <- data.frame(ID = paste0("C", 1:4), A = regimens$Regimen[c(6,6,5,5)], B = regimens$Regimen[c(4,3,2,1)])
pairs$Comparison <- paste(pairs$A, "\u2212", pairs$B)

regimen_status <- function(raw, iso, fields, rules) {
  roster <- unique(c(regimens$A, regimens$B))
  mapping <- setNames(fields$Field, fields$Drug)
  check(!anyDuplicated(fields$Drug) && all(mapping %in% names(raw)) &&
    all(c(roster,"Imipenem","Meropenem") %in% names(mapping)) && "ESBL" %in% names(raw),
    "Required AST field is absent.")
  # Coverage uses the exact raw tokens without normalizing whitespace.
  raw_ast <- as.matrix(as.data.frame(raw[unname(mapping[roster])], check.names = FALSE))
  colnames(raw_ast) <- roster
  ast <- raw_ast
  check(nrow(rules) == 24L && !anyDuplicated(rules[c("Organism","Drug")]) &&
    all(rules$Drug %in% roster & rules$Organism %in% org_names & rules$Value == "R"),
    "Unexpected coverage protocol rule set.")
  for (j in seq_len(nrow(rules))) {
    ast[iso$Organism == rules$Organism[j], rules$Drug[j]] <- "R"
  }
  grid <- expand.grid(A = c(0L,1L,NA_integer_), B = c(0L,1L,NA_integer_))
  check(identical(combine(grid$A, grid$B), c(0L,1L,NA,1L,1L,1L,NA,1L,NA)), "Strict-S truth table failed.")
  binary <- activity(ast)
  out <- sapply(seq_len(nrow(regimens)), function(i) combine(binary[,regimens$A[i]], binary[,regimens$B[i]]))
  colnames(out) <- regimens$Regimen
  # Accepted S10: fungi are structurally uncovered by these antibacterial regimens.
  out[iso$Organism == "Fungi", ] <- 0L
  out
}

weighted_design <- function(z) survey::svydesign(ids = ~Site, strata = ~Region, weights = ~wg, data = z, nest = TRUE)

# Center and aggregate endpoint scores in the full scenario. PSUs outside an
# endpoint domain contribute zero scores, preserving domain Taylor variance.
score_variance <- function(des, keep, response, estimate) {
  z <- des$variables
  u <- numeric(nrow(z))
  u[keep] <- z$wg[keep] * (response[keep] - estimate) / sum(z$wg[keep])
  h <- aggregate(u, z[c("Region","Site")], sum)
  sum(vapply(split(h$x, h$Region, drop = TRUE), function(v) length(v) * var(v), 0))
}

weighted_binary <- function(des, y, eligible, label, coverage = FALSE) {
  check(length(y) == nrow(des$variables) && length(y) == length(eligible) && !anyNA(eligible),
    paste("Invalid weighted endpoint:", label))
  keep <- !is.na(y) & eligible
  n <- sum(keep)
  events <- sum(y[keep])
  des$variables$Response <- as.integer(y)
  sub <- subset(des, keep)
  z <- sub$variables
  df <- survey::degf(sub)
  p <- if (n) {
    if (coverage) sum(z$wg[z$Response == 1L]) / sum(z$wg) else weighted.mean(z$Response,z$wg)
  } else NA_real_
  ans <- data.frame(N = n, Events = events, Unweighted = if (n) events/n else NA_real_,
    Estimate = p, SE = NA_real_, Low = NA_real_, High = NA_real_, DF = df,
    PSU = length(unique(z$Site)), Method = "Taylor linearization, logit t",
    CI_status = "No evaluable results")
  if (!n) return(ans)
  fit <- survey::svymean(~Response, sub)
  se <- sqrt(as.numeric(vcov(fit)))
  ans$SE <- se
  check(near(coef(fit), p), paste("Weighted mean verification failed:", label))
  if (p <= 0 || p >= 1) {
    ans$CI_status <- "Exact boundary: logit CI undefined"
    return(ans)
  }
  if (df <= 0 || !is.finite(se) || se <= 0) {
    ans$CI_status <- "Insufficient design df or no positive estimable variance"
    return(ans)
  }
  ci <- plogis(qlogis(p) + c(-1,1) * qt(.975,df) * se / (p*(1-p)))
  ans$Low <- ci[1]
  ans$High <- ci[2]
  ans$CI_status <- if (min(events, n-events) < 5L) "Estimated; fewer than five events or non-events" else "Estimated"
  check(near(score_variance(des, keep, y, p), se^2, 1e-10),
    paste("Independent weighted Taylor variance check failed:", label))
  if (coverage) {
    interval <- survey::svyciprop(~Response, sub, method = "logit", df = df)
    check(near(as.numeric(confint(interval)), ci, 1e-7), paste("Coverage logit CI check failed:", label))
  }
  ans
}

pair_fit <- function(i, mat, des) {
  ab <- cbind(A = mat[,pairs$A[i]], B = mat[,pairs$B[i]])
  keep <- complete.cases(ab)
  diff <- ab[,1] - ab[,2]
  des$variables$Diff <- diff
  des$variables$A_cov <- ab[,1]
  des$variables$B_cov <- ab[,2]
  sub <- subset(des, keep)
  df <- survey::degf(sub)
  check(df > 0L, paste("No design degrees of freedom in paired domain:", pairs$ID[i]))
  m <- survey::svymean(~Diff + A_cov + B_cov, sub)
  est <- unname(coef(m))
  v <- vcov(m)
  se <- sqrt(v[1,1])
  ci <- est[1] + c(-1,1)*qt(.975,df)*se
  glm <- survey::svyglm(Diff ~ 1, sub)
  check(is.finite(se) && se > 0 && near(est[1], est[2] - est[3]) &&
    near(est[1], weighted.mean(diff[keep], sub$variables$wg)) &&
    near(est[1], coef(glm)[1]) && near(se^2, vcov(glm)[1,1]) &&
    near(se^2, v[2,2] + v[3,3] - 2*v[2,3]) &&
    near(se^2, score_variance(des, keep, diff, est[1])),
    paste("Same-isolate paired variance check failed:", pairs$ID[i]))
  data.frame(N = sum(keep), Events = NA_integer_, Unweighted = mean(diff[keep]), Estimate = est[1],
    SE = se, Low = ci[1], High = ci[2], DF = df, PSU = length(unique(sub$variables$Site)),
    Method = "Taylor linearization, absolute difference t", CI_status = "Estimated",
    Candidate = nrow(mat), A_coverage = est[2], B_coverage = est[3], Cov_AB = v[2,3], Contrast = pairs$ID[i])
}

compute_s22 <- function(d) {
  iso <- d$iso
  cov <- regimen_status(d$raw, iso, d$ast_fields, d$coverage_rules)
  ast <- derive_ecoli_ast(d$raw, iso, d$ast_fields)
  keep <- iso$Site != "Site 23"
  cells <- d$weight_cells
  frac <- function(x) x / ave(x, cells$Year, FUN = sum)
  cells$Retained_n <- tabulate(d$iw[keep], nbins = nrow(cells))
  check(all(cells$Retained_n > 0), "Site 23 exclusion leaves a calibration cell without support.")
  cells$Recalibrated_weight <- cells$Target_birth_share / frac(cells$Retained_n)
  cells$Full_weight_share <- frac(cells$Original_weight * cells$Full_n)
  cells$Fixed_weight_share <- frac(cells$Original_weight * cells$Retained_n)
  cells$Recalibrated_share <- frac(cells$Recalibrated_weight * cells$Retained_n)
  check(near(cells$Full_weight_share, cells$Target_birth_share, 1e-8) &&
    near(cells$Recalibrated_share, cells$Target_birth_share), "Calibration target check failed.")

  scenario <- function(s) {
    rows <- if (s == "A") rep(TRUE,nrow(iso)) else keep
    z <- iso[rows, ]
    if (s == "C") z$wg <- cells$Recalibrated_weight[d$iw[rows]]
    des <- weighted_design(z)
    metadata <- data.frame(Scenario = s, Hospitals = length(unique(z$Site)),
      Isolates = nrow(z), Infants = length(unique(z$InfantID)), Design_df = survey::degf(des),
      Weighting = c(A = "Original S3 weights", B = "Original S3 weights retained",
                    C = "Recalibrated to original S3 targets")[[s]])
    full <- rep(TRUE,nrow(z))
    pathogen <- bind(lapply(c("Overall","Central"), function(g) {
      included <- if (g == "Overall") full else z$Region == g
      bind(lapply(org_names[1:4], function(o) data.frame(Module = "Pathogen", Region = g, Endpoint = o,
        weighted_binary(des, z$Organism == o, included, paste(s,g,o)))))
    }))
    ec <- z$Organism == "E. coli"
    ast_results <- bind(lapply(names(ast$values), function(name) {
      v <- ast$values[[name]][rows]
      fit <- weighted_binary(des, v == 3L, !is.na(v), paste(s,name))
      check(sum(v %in% 1:3) == fit$N, paste("Invalid AST denominator:", name))
      data.frame(Module = "AST_ESBL", Region = "Overall", Endpoint = name,
        fit, S = sum(v %in% 1L), I = sum(v %in% 2L), Organism_n = sum(ec))
    }))
    ast_results <- rbind(ast_results, data.frame(Module = "AST_ESBL", Region = "Overall",
      Endpoint = "Reported ESBL positivity", weighted_binary(des, ast$esbl_positive[rows], ec, paste(s,"ESBL")),
      S = NA_integer_, I = NA_integer_, Organism_n = sum(ec)))
    mat <- cov[rows, , drop = FALSE]
    coverage <- bind(lapply(colnames(mat), function(name) {
      y <- mat[,name]
      f <- weighted_binary(des, y, full, paste(s,name), coverage = TRUE)
      n0 <- sum(y %in% 0L)
      nu <- sum(is.na(y))
      check(f$Events + n0 + nu == length(y) && f$Events + n0 == f$N, paste("Coverage count check:", name))
      data.frame(Module = "Coverage", Region = "Overall", Endpoint = name, f, Candidate = length(y),
        Covered = f$Events, Uncovered = n0, Unknown = nu, Evaluable = f$N)
    }))
    paired <- bind(lapply(seq_len(nrow(pairs)), function(i) data.frame(Module = "Paired", Region = "Overall",
      Endpoint = pairs$Comparison[i], pair_fit(i, mat, des))))
    list(metadata = metadata, results = lapply(list(pathogen, ast_results, coverage, paired),
      function(x) cbind(Scenario = s, x)))
  }

  estimates <- lapply(c("A","B","C"), scenario)
  metadata <- bind(lapply(estimates, "[[", "metadata"))
  check(identical(metadata$Hospitals, c(30L,29L,29L)) &&
    identical(metadata$Isolates, c(1196L,1110L,1110L)) && identical(metadata$Infants, c(1189L,1103L,1103L)) &&
    all(metadata$Design_df == c(26,25,25)), "Unexpected S22 scenario cohort or design counts.")
  pieces <- lapply(seq_len(4L), function(k) bind(lapply(estimates, function(x) x$results[[k]])))
  cols <- unique(unlist(lapply(pieces, names)))
  add_cols <- function(x) { for (f in setdiff(cols,names(x))) x[[f]] <- NA; x[cols] }
  results <- bind(lapply(pieces, add_cols))
  ids <- paste(results$Module, results$Region, results$Endpoint, sep = "/")
  a <- results$Scenario == "A"
  c <- results$Scenario == "C"
  a_est <- setNames(results$Estimate[a], ids[a])
  c_est <- setNames(results$Estimate[c], ids[c])
  # Estimate and SE are stored on the probability/difference scale throughout.
  results$Change_vs_A_pp <- 100 * (results$Estimate - unname(a_est[ids]))
  results$C_minus_A_pp <- 100 * unname(c_est[ids] - a_est[ids])
  results$Display_unit <- ifelse(results$Module == "Paired", "percentage points", "percent")
  results$Row_ID <- ids
  check(sum(a) == 23L && nrow(results) == 69L && !anyDuplicated(paste(ids, results$Scenario)),
    "Incorrect S22 output dimensions or row keys.")
  es <- results$Endpoint == "Reported ESBL positivity"
  check(all(results$Events[es] == c(65L,43L,43L)) && all(results$N[es] == c(275L,234L,234L)),
    "ESBL flag or E. coli denominator differs from the accepted release.")
  list(results = results, weights = cells, scenarios = metadata)
}


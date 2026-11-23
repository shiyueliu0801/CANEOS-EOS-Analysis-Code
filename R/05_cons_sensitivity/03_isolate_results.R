# Taxon estimates use the accepted Table 2 mapping and S19 exclusion domains.
pathogen_ids <- c("Gram-positive", "Gram-negative", "E. coli", "GBS", "CoNS",
  "Klebsiella spp.", "Enterococcus spp.", "L. monocytogenes")
pathogen_labels <- c("Gram-positive bacteria", "Gram-negative bacteria", "E. coli", "GBS", "CoNS",
  "Klebsiella spp.", "Enterococcus spp.", "Listeria monocytogenes")
pathogens <- bind(lapply(seq_along(pathogen_ids), function(i) {
  id <- pathogen_ids[i]
  y <- as.integer((if (i <= 2L) iso$Gram else iso$Organism) == id)
  props <- lapply(isolate_domains, function(keep) proportion(isolate_design, y, keep))
  excluded <- id == "CoNS"
  data.frame(Code = id, Indicator = pathogen_labels[i], Excluded = excluded, wide(props),
    Change_pp = if (excluded) NA_real_ else 100*(props$NonCoNS$Estimate-props$Overall$Estimate))
}))

# Nine named source columns feed the 12 accepted S10 regimens. The complete
# approved rule matrix is verified so this subset cannot silently use a different matrix.
ast_fields <- c(Penicillin = "青霉素（窄谱青霉素）", Ampicillin = "氨苄西林（广谱青霉素）",
  Ceftazidime = "头孢他啶（三代）", Cefotaxime = "头孢噻肟（三代）", Gentamicin = "庆大霉素（氨基糖苷）",
  Amikacin = "阿米卡星（氨基糖苷）", Vancomycin = "万古霉素", Imipenem = "亚胺培南", Meropenem = "美罗培南")
stopifnot(all(ast_fields %in% names(raw)[46:74]))
ast <- as.matrix(raw[unname(ast_fields)])
ast[is.na(ast)] <- ""
colnames(ast) <- names(ast_fields)
stopifnot(all(ast %in% c("", "敏感", "中介", "耐药", "S", "I", "R")))

rules_raw <- static_cfg$get_coverage_intrinsic_rules()
rules <- data.frame(Organism = rules_raw$Organism, Drug = rules_raw$Antibiotic, stringsAsFactors = FALSE)
stopifnot(nrow(rules) == 83L, !anyDuplicated(rules), all(rules$Organism %in% iso$Organism))
for (i in which(rules$Drug %in% colnames(ast))) {
  ast[iso$Organism == rules$Organism[i], rules$Drug[i]] <- "R"
}

# Strict-S tri-state activity and coverage are extracted from the accepted
# CANEOS_Coverage_11.15.5.R. Fungal isolates are structural non-coverage.
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
truth <- expand.grid(A = c(0L, 1L, NA), B = c(0L, 1L, NA))
stopifnot(identical(combine(truth$A, truth$B), c(0L, 1L, NA, 1L, 1L, 1L, NA, 1L, NA)))
regimens <- data.frame(Code = c("PEN + CAZ", "PEN + CTX", "AMP + CAZ", "AMP + CTX", "PEN + GEN", "PEN + AMK",
  "AMP + GEN", "AMP + AMK", "VAN + IPM", "VAN + MEM", "VAN + GEN", "VAN + AMK"),
  A = c("Penicillin", "Penicillin", "Ampicillin", "Ampicillin", "Penicillin", "Penicillin",
    "Ampicillin", "Ampicillin", rep("Vancomycin", 4L)),
  B = c("Ceftazidime", "Cefotaxime", "Ceftazidime", "Cefotaxime", "Gentamicin", "Amikacin",
    "Gentamicin", "Amikacin", "Imipenem", "Meropenem", "Gentamicin", "Amikacin"))
binary <- activity(ast)
fungus <- iso$Organism == "Fungi"
covered <- sapply(seq_len(nrow(regimens)), function(i) combine(binary[,regimens$A[i]], binary[,regimens$B[i]]))
covered[fungus,] <- 0L
colnames(covered) <- regimens$Code

coverage <- bind(lapply(seq_len(nrow(regimens)), function(i) {
  y <- covered[,i]
  props <- lapply(isolate_domains, function(keep) {
    a <- proportion(isolate_design, y, keep)
    a$Unknown <- a$Candidate - a$N
    a$Unknown_rate <- a$Unknown/a$Candidate
    a$Weighted_unknown_rate <- sum(iso$Weight[keep & is.na(y)])/sum(iso$Weight[keep])
    a$Evaluable_rate <- a$N/a$Candidate
    a
  })
  cons_n <- sum(isolate_cons)
  cons_unknown <- sum(isolate_cons & is.na(y))
  data.frame(Code = regimens$Code[i], Indicator = paste(regimens$A[i], regimens$B[i], sep = " + "),
    wide(props), CoNS_Candidate = cons_n, CoNS_Evaluable = cons_n-cons_unknown,
    CoNS_Unknown = cons_unknown, CoNS_Evaluable_rate = (cons_n-cons_unknown)/cons_n,
    Change_pp = 100*(props$NonCoNS$Estimate - props$Overall$Estimate))
}))
stopifnot(all(coverage$Overall_Unknown == coverage$NonCoNS_Unknown + coverage$CoNS_Unknown),
  all(coverage$Overall_N == coverage$NonCoNS_N + coverage$CoNS_Evaluable))

# Exact same-isolate differences and joint missingness bounds from S11.
pairs <- data.frame(Code = paste0("C", 1:4),
  A = c("AMP + GEN", "AMP + GEN", "PEN + GEN", "PEN + GEN"),
  B = c("AMP + CTX", "AMP + CAZ", "PEN + CTX", "PEN + CAZ"),
  Shared = rep(c("Ampicillin", "Penicillin"), each = 2L),
  Right = rep(c("Cefotaxime", "Ceftazidime"), 2L))
pairs$Indicator <- paste(pairs$A, "−", pairs$B)

paired_fit <- function(a, b, eligible) {
  keep <- eligible & !is.na(a) & !is.na(b)
  y <- a - b
  d <- domain(isolate_design, keep, Difference = y, A = a, B = b)
  fit <- survey::svymean(~Difference + A + B, d)
  est <- unname(coef(fit))
  se <- as.numeric(survey::SE(fit)[1L])
  df <- survey::degf(d)
  stopifnot(is.finite(se), se > 0, df > 0, near(est[1L], est[2L]-est[3L]),
    near(est[1L], weighted.mean(y[keep], iso$Weight[keep])),
    near(se^2, psu_variance(isolate_design, y, keep, est[1L])))
  lim <- 100*(est[1L] + c(-1, 1)*qt(.975, df)*se)
  data.frame(Candidate = sum(eligible), N = sum(keep), A_coverage = est[2L], B_coverage = est[3L],
    Difference_pp = 100*est[1L], Lower = lim[1L], Upper = lim[2L], SE_pp = 100*se,
    t_statistic = est[1L]/se, DF = df, P = 2*pt(-abs(est[1L]/se), df))
}

paired_domains <- lapply(isolate_domains, function(eligible) {
  r <- bind(lapply(seq_len(nrow(pairs)), function(i) paired_fit(covered[,pairs$A[i]], covered[,pairs$B[i]], eligible)))
  adjust_bh(r, rep("Four overall paired comparisons", nrow(r)))
})
paired <- data.frame(pairs[c("Code","Indicator")], wide(paired_domains),
  Change_pp = paired_domains$NonCoNS$Difference_pp - paired_domains$Overall$Difference_pp)

pair_limits <- function(x) {
  grid <- expand.grid(Shared = 0:1, Left = 0:1, Right = 0:1)
  d <- matrix(NA_integer_, nrow(x), nrow(grid))
  for (j in seq_len(nrow(grid))) {
    vals <- unlist(grid[j,])
    may <- rowSums(sweep(x, 2, vals, "!="), na.rm = TRUE) == 0L
    d[may,j] <- as.integer(max(vals[c(1,2)]) - max(vals[c(1,3)]))
  }
  stopifnot(all(rowSums(!is.na(d)) > 0L))
  lo <- apply(d,1,min,na.rm = TRUE)
  hi <- apply(d,1,max,na.rm = TRUE)
  lo[fungus] <- hi[fungus] <- 0L
  list(Lower = lo, Upper = hi)
}
bound_mean <- function(y, keep) {
  d <- domain(isolate_design, keep, Result = y)
  m <- survey::svymean(~Result, d)
  stopifnot(near(coef(m), weighted.mean(y[keep], iso$Weight[keep])),
    near(as.numeric(vcov(m)), psu_variance(isolate_design, y, keep, as.numeric(coef(m)))))
  c(Estimate = as.numeric(coef(m))*100, SE = as.numeric(survey::SE(m))*100, DF = survey::degf(d))
}
paired_bounds <- bind(lapply(names(isolate_domains), function(group) {
  keep <- isolate_domains[[group]]
  bind(lapply(seq_len(nrow(pairs)), function(i) {
    x <- pair_limits(binary[,c(pairs$Shared[i], "Gentamicin", pairs$Right[i]),drop = FALSE])
    both <- keep & complete.cases(covered[,c(pairs$A[i],pairs$B[i])])
    observed <- covered[both,pairs$A[i]] - covered[both,pairs$B[i]]
    stopifnot(all(x$Lower[both] == observed), all(x$Upper[both] == observed))
    l <- bound_mean(x$Lower, keep)
    u <- bound_mean(x$Upper, keep)
    df <- unname(l["DF"])
    stopifnot(df == u["DF"], l["Estimate"] <= u["Estimate"])
    t <- qt(.975, df)
    tf <- qt(1-.05/(2*nrow(pairs)), df)
    data.frame(Group = group, Code = pairs$Code[i], Indicator = pairs$Indicator[i], Target_N = sum(keep),
      Common_N = sum(both), Uncertain_isolates = sum(keep & x$Lower < x$Upper), DF = df,
      Sample_lower_pp = unname(l["Estimate"]), Sample_upper_pp = unname(u["Estimate"]),
      Lower_SE = unname(l["SE"]), Upper_SE = unname(u["SE"]),
      Joint95_low = unname(l["Estimate"]-t*l["SE"]), Joint95_high = unname(u["Estimate"]+t*u["SE"]),
      Family95_low = unname(l["Estimate"]-tf*l["SE"]), Family95_high = unname(u["Estimate"]+tf*u["SE"]))
  }))
}))

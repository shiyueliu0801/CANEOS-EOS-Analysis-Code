# ==============================================================================
# Module 03: Empirical Antibiotic Coverage Statistical Core
# ==============================================================================
# Preserves all verified statistical functions, survey designs, strict-S rules,
# paired differences, sharp logical bounds, S9/S10/S14 tables, and S2/S3/Fig 2
# backing data from the locked revision package.
# ==============================================================================

source(file.path(mod_dir, "01_prepare_inputs.R"), local = TRUE, encoding = "UTF-8")
check(nrow(rules) == 83L && all(rules$Matrix_value == "R"),
      "Coverage intrinsic rules", "All 83 locked coverage-only intrinsic resistance rules loaded from static config.")
cell_key <- function(d) paste(d$Isolate_ID, d$Antibiotic, sep = "|")
cell_pos <- function(d) cbind(match(d$Isolate_ID, iso$Isolate_ID), match(d$Antibiotic, drug_fields))
positions <- cell_pos(deferred)
rule_id <- matrix(NA_integer_,n,length(drug_fields),dimnames=dimnames(raw_ast))
for (j in seq_len(nrow(rules))) {
  k <- iso$Organism_mapped == rules$Organism[j]
  rule_id[k,match(rules$Antibiotic[j],drug_fields)] <- j
}
rule_positions <- !is.na(rule_id)
mask24 <- mask6 <- matrix(FALSE,n,length(drug_fields),dimnames=dimnames(raw_ast))
mask24[positions] <- TRUE
restored_groups <- c("GBS","Enterococcus spp.","Other Streptococcus spp.")
six_rules <- which(rules$Organism %in% restored_groups & rules$Antibiotic %in% c("Gentamicin","Amikacin"))
if (length(six_rules)!=6L) {
  gbs_group <- unique(rules$Organism[grepl("GBS|agalactiae|Group B",rules$Organism,ignore.case=TRUE)])
  restored_groups[1] <- gbs_group
  six_rules <- which(rules$Organism %in% restored_groups & rules$Antibiotic %in% c("Gentamicin","Amikacin"))
}
mask6[] <- rule_id %in% six_rules
check(nrow(deferred)==24L && !anyDuplicated(cell_key(deferred)) &&
  all(rule_positions[positions]) && !any(mask6 & mask24) && length(six_rules)==6L &&
  identical(rules$Matrix_cell[rule_id[positions]],deferred$Matrix_cell),
  "Six matrix cells and 24 isolate cells", "The six low-level aminoglycoside rules and the 24 approved species-scope isolate cells are distinct. Both rollbacks use the current raw cells.")
apply_rules <- function(ast, mask) {
  ast[mask] <- "R"
  ast
}
final_ast <- apply_rules(raw_ast,rule_positions)
no24_ast <- apply_rules(raw_ast,rule_positions & !mask24)
no6_ast <- apply_rules(raw_ast,rule_positions & !mask6)
ast_code <- function(x) {
  unname(c("敏感"="S","中介"="I","耐药"="R",S="S",I="I",R="R")[x])
}
change_type <- function(x) paste0(ifelse(x=="","Missing",ast_code(x))," to R")
rule_coords <- which(rule_positions,arr.ind=TRUE)
ri <- rule_id[rule_coords]
assignment <- data.frame(iso[rule_coords[,1],c("Isolate_ID","Source_row","Site","Hospital","Year","Region",
    "Organism_mapped","Species_mapped")],
  Antibiotic=drug_fields[rule_coords[,2]], Matrix_cell=rules$Matrix_cell[ri],
  Raw_AST=raw_ast[rule_coords], Analysis_AST=final_ast[rule_coords],
  Disposition=change_type(raw_ast[rule_coords]),
  Approved_24_cell=mask24[rule_coords], Six_low_level_AG_rule=mask6[rule_coords],
  Rationale=rules$Rule[ri], row.names=NULL)
counts_rule <- function(x) {
  v <- ast_code(x)
  c(Applications=length(v), Missing_to_R=sum(x==""), S_to_R=sum(v %in% "S"),
    I_to_R=sum(v %in% "I"), Retained_R=sum(v %in% "R"))
}
rule_summary <- data.frame(rules[c("Matrix_cell","Organism","Antibiotic","Field_type","Rule","Scope","Reference")],
  t(vapply(seq_len(nrow(rules)),function(j) counts_rule(raw_ast[which(rule_id==j)]),numeric(5))),
  In_six_rule_rollback=seq_len(nrow(rules)) %in% six_rules, row.names=NULL)
rule_totals <- data.frame(Release="11.15.5 Phase 2.2",as.list(counts_rule(raw_ast[rule_positions])))
check(identical(final_ast[!rule_positions],raw_ast[!rule_positions]) &&
  identical(final_ast[fungus | unidentified,],raw_ast[fungus | unidentified,]) &&
  identical(no24_ast[mask24],raw_ast[mask24]) && identical(no6_ast[mask6],raw_ast[mask6]) &&
  all(rowSums(rule_summary[c("Missing_to_R","S_to_R","I_to_R","Retained_R")]) == rule_summary$Applications), "Matrix scope and AST invariance",
  sprintf("%d eligible cells. %d original S/I. All other cells and raw AST preserved.", nrow(assignment),
    sum(assignment$Disposition %in% c("S to R","I to R"))))
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
encode <- function(x) c("Uncovered", "Covered", "Unknown")[match(x, c(0L, 1L, NA_integer_))]
st <- c(0L,1L,NA_integer_)
grid <- expand.grid(A = st, B = st)
truth <- c(0L,1L,NA,1L,1L,1L,NA,1L,NA)
check(identical(combine(grid$A,grid$B), truth), "Strict-S truth table",
  "All nine possible combinations verified independently of the isolate data.")
regimens <- data.frame(Regimen = c("PEN + CAZ", "PEN + CTX", "AMP + CAZ", "AMP + CTX", "PEN + GEN", "AMP + GEN"),
  A = c("Penicillin", "Penicillin", "Ampicillin", "Ampicillin", "Penicillin", "Ampicillin"),
  B = c("Ceftazidime", "Cefotaxime", "Ceftazidime", "Cefotaxime", "Gentamicin", "Gentamicin"))
regimens$Family <- rep(c("China", "WHO"), c(4, 2))
regs <- regimens$Regimen
pairs <- data.frame(ID = 1:4, Regimen_A = regs[c(6,6,5,5)], Regimen_B = regs[c(4,3,2,1)],
  Shared = rep(c("Ampicillin", "Penicillin"), each = 2), Left = "Gentamicin",
  Right = rep(c("Cefotaxime", "Ceftazidime"), 2))
pairs$Comparison <- paste0("C", pairs$ID)
binary <- activity(final_ast)
raw_binary <- activity(raw_ast)
no24_binary <- activity(no24_ast)
no6_binary <- activity(no6_ast)
coverage <- function(ast, roster, fungi_zero) {
  out <- sapply(seq_len(nrow(roster)), function(i) combine(ast[, roster$A[i]], ast[, roster$B[i]]))
  colnames(out) <- roster$Regimen
  if (fungi_zero) out[fungus, ] <- 0L
  out
}
main <- coverage(binary,regimens,TRUE)
no24 <- coverage(no24_binary,regimens,TRUE)
no6 <- coverage(no6_binary,regimens,TRUE)
table_counts <- function(m) t(vapply(seq_len(ncol(m)), function(i) {
  x <- m[,i]
  c(Covered = sum(x %in% 1L), Uncovered = sum(x %in% 0L), Unknown = sum(is.na(x)))
}, numeric(3)))
check(all(rowSums(table_counts(main))==n) && all(main[fungus,]==0L) &&
  identical(main[unidentified,],coverage(raw_binary,regimens,TRUE)[unidentified,]),
  "Six exhaustive statuses and special organisms", "All six tri-state totals equal 1,196. 20 fungi are structural zeros only for combinations. Four unidentified isolates retain their original AST.")



# Complex-survey analysis
w_region <- wg$Region
w_year <- as.integer(wg$Year)
wk <- paste(w_region,w_year,sep="/")
wi <- match(paste(iso$Region,iso$Year,sep="/"),wk)
site_map <- unique(iso[,c("Site","Region")])
site_map <- site_map[order(as.integer(sub("^Site ","",site_map$Site))),]
site_n <- table(factor(site_map$Region,levels=regions))
weight_n <- tabulate(wi,nbins=length(wk))
births <- wg$Births
fyear <- function(x) x/ave(x,w_year,FUN=sum)
check(nrow(site_map)==30L && !anyDuplicated(site_map$Site) &&
  identical(as.integer(site_n),c(9L,7L,7L,7L)) && nrow(wg)==40L &&
  !anyNA(wi) && !anyDuplicated(wk) && all(is.finite(w)&w>0) &&
  near(w,wg$Weight[wi],1e-12) && all(weight_n==wg$Isolates) &&
  near(wg$Weight,fyear(births)/fyear(weight_n),1e-8),
  "PSUs, strata and S3 weights","Thirty hospitals, four strata and all forty locked calibration cells")
weight_cells <- data.frame(Region=w_region,Year=w_year,Isolates=weight_n,
  S3_sample_share=fyear(weight_n),S3_birth_share=fyear(births),S3_weight=wg$Weight)
weight_distribution <- data.frame(Statistic=c("Min","P01","P05","Median","P95","P99","Max","Sum","Mean"),
  S3_weight=unname(c(quantile(w,c(0,.01,.05,.5,.95,.99,1)),sum(w),mean(w))))
iso$Unit_weight <- 1
iso$Year_centered <- iso$Year - 2015L
options(survey.lonely.psu = "fail", survey.adjust.domain.lonely = FALSE)
design_full <- survey::svydesign(ids = ~Site, strata = ~Region, weights = ~S3_weight, data = iso, nest = TRUE)
design_unit <- survey::svydesign(ids = ~Site, strata = ~Region, weights = ~Unit_weight, data = iso, nest = TRUE)
check(survey::degf(design_full) == 26L && survey::degf(design_unit) == 26L, "Full design degrees of freedom",
  "30 hospitals minus four strata gives 26 design df for both weighting methods.")
psu_qc <- list()
domain <- function(des, keep, label, detail = TRUE) {
  stopifnot(nrow(des$variables) == n, length(keep) == n, !anyNA(keep))
  z <- subset(des, keep)
  if (detail) {
    ds <- z$variables
    counts <- table(factor(unique(ds[c("Site","Region")])$Region, levels = regions))
    psu_qc[[length(psu_qc) + 1L]] <<- data.frame(Domain = label, Isolates = nrow(ds), Design_df = survey::degf(z),
      Contributing_PSUs = sum(counts), Contributing_strata = sum(counts > 0), Singleton_domains = sum(counts == 1),
      Municipality = counts[1], Eastern = counts[2], Central = counts[3], Western = counts[4], row.names = NULL)
  }
  z
}
injected <- function(des, ...) {
  extra <- list(...)
  for (k in names(extra)) des$variables[[k]] <- extra[[k]]
  des
}
wm <- function(x, wt = w) sum(x * wt) / sum(wt)
nz <- function(num, den) num / ifelse(den > 0, den, NA_real_)
prop_fit <- function(y, eligible, label, des = design_full, infer = TRUE, detail = TRUE) {
  den <- !is.na(y) & eligible
  wt <- weights(des)
  p <- nz(sum(wt[den & y %in% 1L]), sum(wt[den]))
  result <- c(Coverage = p, SE = NA_real_, Low = NA_real_, High = NA_real_, DF = NA_real_)
  if (!any(den) || !infer) return(list(est = result, method = "Descriptive"))
  d <- domain(injected(des, Response = y), den, label, detail)
  df <- survey::degf(d)
  m <- survey::svymean(~Response, d)
  s <- sqrt(as.numeric(vcov(m)))
  stopifnot(near(as.numeric(coef(m)), p))
  result[c("SE", "DF")] <- c(s, df)
  if (p <= 0 || p >= 1 || df <= 0 || !is.finite(s) || s <= 0) {
    return(list(est = result, method = "CI unavailable (boundary or no estimable variance)"))
  }
  direct <- plogis(qlogis(p) + c(-1,1) * qt(.975, df) * s / (p * (1 - p)))
  if (detail) {
    interval <- survey::svyciprop(~Response, d, method = "logit", df = df)
    ci <- as.numeric(confint(interval))
    stopifnot(near(ci, direct, 1e-7))
  }
  stopifnot(all(is.finite(direct)), direct[1] <= p, p <= direct[2])
  result[c("Low", "High")] <- direct
  list(est = result, method = "Taylor linearisation, t-logit")
}
cov_row <- function(y, k, elig, subgroup, label, infer = TRUE, detail = TRUE) {
  a <- elig
  c1 <- sum(a & y %in% 1L)
  c0 <- sum(a & y %in% 0L)
  cu <- sum(a & is.na(y))
  f <- prop_fit(y, a, label, infer = infer, detail = detail)
  u <- prop_fit(y, a, paste(label, "Unit"), design_unit, infer = infer, detail = FALSE)
  data.frame(Subgroup = subgroup, Regimen = k, Candidate = sum(a), Covered = c1, Uncovered = c0,
    Unknown = cu, Evaluable = c1+c0, Unweighted_coverage = nz(c1,c1+c0), Unit_CI_low = u$est[3], Unit_CI_high = u$est[4],
    Unit_SE = u$est[2], Unit_DF = u$est[5], Unit_CI_method = u$method, Weighted_coverage = f$est[1],
    Weighted_CI_low = f$est[3], Weighted_CI_high = f$est[4], Weighted_SE = f$est[2],
    Unknown_rate = nz(cu,sum(a)), Weighted_unknown_rate = nz(sum(w[a & is.na(y)]),sum(w[a])),
    Covered_share_all = nz(c1,sum(a)), Weighted_covered_share_all = nz(sum(w[a & y %in% 1L]),sum(w[a])),
    DF = f$est[5], CI_method = f$method, row.names = NULL)
}
all_isolates <- rep(TRUE,n)
coverage_rows <- function(m, elig = all_isolates, name = "Overall", prefix = "A", infer = TRUE, detail = TRUE) {
  bind(lapply(seq_len(ncol(m)), function(i) cov_row(m[,i], colnames(m)[i], elig, name,
    paste(prefix, name, colnames(m)[i], sep = " / "), infer, detail)))
}
table_a <- coverage_rows(main)[,-1]
check(all(is.finite(table_a$Weighted_CI_low)) && all(table_a$DF == 26L),
  "Primary coverage intervals", "Six t-logit intervals match explicit delta-method intervals and all use 26 design df.")
manual_variance <- function(d, keep, wt, estimate) {
  u <- numeric(n)
  u[keep] <- wt[keep] * (d[keep] - estimate) / sum(wt[keep])
  z <- as.numeric(rowsum(u, iso$Site, reorder = FALSE)[match(site_map$Site, unique(iso$Site)),1])
  sum(vapply(split(z, site_map$Region), function(v) length(v) * stats::var(v), 0))
}
p_index <- function(i, m = main) cbind(A = m[,pairs$Regimen_A[i]], B = m[,pairs$Regimen_B[i]])
pair_label <- function(i) paste(pairs$Regimen_A[i], "\u2212", pairs$Regimen_B[i])
crosscheck <- list()
pair_fit <- function(i, mat, elig = all_isolates, des = design_full, label = "Primary") {
  ab <- p_index(i, mat)
  keep <- elig & complete.cases(ab)
  d <- ab[,1] - ab[,2]
  wt <- weights(des)
  sub <- domain(injected(des, Diff = d, A_cov = ab[,1], B_cov = ab[,2]), keep, paste(label, pairs$Comparison[i], sep = " / "))
  df <- survey::degf(sub)
  if (df <= 0L) stop("No degrees of freedom in paired domain: ", label, " ", pairs$Comparison[i])
  m <- survey::svymean(~Diff + A_cov + B_cov, sub)
  est <- unname(coef(m))
  se <- sqrt(vcov(m)[1,1])
  glm <- survey::svyglm(Diff ~ 1, sub)
  t <- est[1] / se
  ci <- est[1] + c(-1,1)*qt(.975,df)*se
  p <- 2 * pt(-abs(t),df)
  vm <- manual_variance(d, keep, wt, est[1])
  direct <- wm(d[keep], wt[keep])
  stopifnot(is.finite(se), se > 0, near(est[1], est[2] - est[3]), near(est[1], direct),
    near(est[1], coef(glm)[1]), near(se^2, vm), near(se^2, vcov(glm)[1,1]),
    near(p, 2 * pt(-abs(coef(glm)[1]/sqrt(vcov(glm)[1,1])),df)))
  crosscheck[[length(crosscheck)+1L]] <<- data.frame(Analysis = label, Contrast = pairs$Comparison[i],
    Mean_difference = est[1], Mean_A_minus_B = est[2]-est[3], Weighted_sum_ratio = direct,
    GLM_intercept = unname(coef(glm)[1]), Taylor_SE_pp = 100*se, Manual_SE_pp = 100*sqrt(vm),
    GLM_SE_pp = 100*sqrt(vcov(glm)[1,1]))
  data.frame(Analysis = label, Contrast = pairs$Comparison[i], Comparison = pair_label(i), Candidate = sum(elig),
    Both_evaluable = sum(keep), Retention_rate = mean(keep), A_coverage = est[2], B_coverage = est[3],
    Difference_pp = 100*est[1], SE_pp = 100*se, CI_low_pp = 100*ci[1],
    CI_high_pp = 100*ci[2], t_statistic = t, DF = df, P_value = p)
}
fit_pairs <- function(mat, des, label, full_panel = FALSE) {
  bind(lapply(seq_len(nrow(pairs)), function(i) {
    eligible <- if (full_panel) complete.cases(binary[,unlist(pairs[i,c("Shared","Left","Right")])]) else all_isolates
    pair_fit(i, mat, eligible, des, label)
  }))
}
table_b <- fit_pairs(main, design_full, "Primary S3")
table_b$BH_q <- p.adjust(table_b$P_value, "BH")
pair_unit <- fit_pairs(main,design_unit,"Unit")
pair_panel <- fit_pairs(main,design_full,"Three determined components",full_panel = TRUE)
check(all(table_b$DF==26L) && identical(pair_unit$Both_evaluable,table_b$Both_evaluable),
  "Current paired domains","Four same-isolate domains retain the accepted full-design construction")
case_order <- c("Both evaluable", "A only evaluable", "B only evaluable", "Neither evaluable")
availability <- function(a, b) {
  factor(ifelse(is.na(a), ifelse(is.na(b), case_order[4], case_order[3]),
    ifelse(is.na(b), case_order[2], case_order[1])), levels = case_order)
}
mix_rows <- function(i, keep, label) {
  a <- table(factor(iso$Organism_mapped, levels = study_groups))
  b <- table(factor(iso$Organism_mapped[keep], levels = study_groups))
  wn <- vapply(study_groups, function(k) sum(w[iso$Organism_mapped == k]), 0)
  we <- vapply(study_groups, function(k) sum(w[keep & iso$Organism_mapped == k]), 0)
  data.frame(Contrast = pairs$Comparison[i], Domain = label, Organism = study_groups,
    Candidate = as.integer(a), Evaluable = as.integer(b), Cohort_share = as.numeric(a)/n,
    Domain_share = as.numeric(b)/sum(keep), Weighted_cohort_share = wn/sum(w),
    Weighted_domain_share = we/sum(w[keep]), Retention = as.numeric(b)/as.numeric(a), row.names = NULL)
}
pair_struct <- lapply(seq_len(nrow(pairs)), function(i) {
  ab <- p_index(i)
  ag <- availability(ab[,1], ab[,2])
  keep <- complete.cases(ab)
  three_known <- complete.cases(binary[, unlist(pairs[i,c("Shared","Left","Right")])])
  stopifnot(!any(three_known & fungus), all(!three_known | keep))
  et <- table(ag)
  ct <- table(factor(ab[keep,1], c(1,0)), factor(ab[keep,2], c(1,0)))
  ev <- data.frame(Contrast = pairs$Comparison[i], Comparison = pair_label(i), Relation = names(et),
    N = as.integer(et), Share_all = as.integer(et)/n,
    Weighted_share_all = sapply(names(et), function(x) sum(w[ag == x])/sum(w)), row.names = NULL)
  ct_rows <- data.frame(Contrast = pairs$Comparison[i], Comparison = pair_label(i), A_status = rep(c("Covered","Uncovered"), each = 2),
    B_status = rep(c("Covered","Uncovered"), 2), N = as.numeric(t(unclass(ct))), Common_N = sum(keep))
  list(ev = ev, ct = ct_rows, mix = bind(list(mix_rows(i, keep, "Common evaluable"),
    mix_rows(i, three_known, "Three determined components"))), keep = keep, panel = three_known)
})
extract <- function(x, name) bind(lapply(x, function(g) g[[name]]))
pair_evaluability <- extract(pair_struct,"ev")
pair_coverage <- extract(pair_struct,"ct")
pair_mix <- extract(pair_struct,"mix")
pair_mix_metrics <- bind(lapply(seq_len(nrow(pairs)), function(i) {
  g <- pair_mix[pair_mix$Contrast == pairs$Comparison[i], ]
  a <- g$Domain == "Common evaluable"
  d <- g$Domain_share[a] - g$Cohort_share[a]
  p <- g$Domain_share[!a] - g$Domain_share[a]
  data.frame(Contrast = pairs$Comparison[i], Comparison = pair_label(i),
    Main_vs_cohort_TV = sum(abs(d))/2,
    Weighted_main_vs_cohort_TV = sum(abs(g$Weighted_domain_share[a] - g$Weighted_cohort_share[a]))/2,
    Largest_main_shift = g$Organism[a][which.max(abs(d))], Main_shift_pp = 100*d[which.max(abs(d))],
    Panel_vs_main_TV = sum(abs(p))/2, Largest_panel_shift = g$Organism[a][which.max(abs(p))],
    Panel_shift_pp = 100*p[which.max(abs(p))])
}))
check(all(tapply(pair_evaluability$N, pair_evaluability$Contrast, sum) == n) &&
  all(as.numeric(tapply(pair_coverage$N, pair_coverage$Contrast, sum)) == table_b$Both_evaluable) &&
  all(vapply(pair_struct,function(z) sum(z$panel),0L) == pair_panel$Both_evaluable),
  "Paired tables and panel nesting", paste("Evaluability totals", paste(tapply(pair_evaluability$N, pair_evaluability$Contrast, sum), collapse = "/"),
    "| coverage totals", paste(tapply(pair_coverage$N, pair_coverage$Contrast, sum), collapse = "/"),
    "| three-component counts", paste(vapply(pair_struct,function(z) sum(z$panel),0L), collapse = "/"),
    "| All panel records lie within the main domains."))
pattern_names <- c("A non-S / B missing", "A missing / B non-S", "Both missing")
pat <- function(a, b, cov) {
  x <- rep(NA_character_, length(a))
  x[is.na(cov) & a %in% 0L & is.na(b)] <- pattern_names[1]
  x[is.na(cov) & is.na(a) & b %in% 0L] <- pattern_names[2]
  x[is.na(cov) & is.na(a) & is.na(b)] <- pattern_names[3]
  stopifnot(identical(is.na(x), !is.na(cov)))
  x
}
missing_patterns <- sapply(seq_len(nrow(regimens)), function(i) pat(binary[,regimens$A[i]], binary[,regimens$B[i]], main[,i]))
colnames(missing_patterns) <- regs
pat_records <- bind(lapply(seq_len(nrow(regimens)), function(i) {
  k <- which(!is.na(missing_patterns[,i]))
  data.frame(iso[k,c("Isolate_ID","Source_row","Site","Year","Region","Organism_mapped","Source_species_label","S3_weight")],
    Regimen = regs[i], Component_A = regimens$A[i], Component_B = regimens$B[i],
    A_AST = final_ast[k,regimens$A[i]], B_AST = final_ast[k,regimens$B[i]], Pattern = missing_patterns[k,i])
}))
tab_pattern <- function(lev, labels = unique(lev), variable = "Group") {
  bind(lapply(seq_len(nrow(regimens)), function(i) {
    a <- missing_patterns[,i]
    bind(lapply(labels, function(g) {
      j <- lev == g
      t <- table(factor(a[j], levels = pattern_names))
      num <- as.integer(t)
      wz <- vapply(pattern_names, function(p) sum(w[j & a %in% p]), 0)
      data.frame(Regimen = regs[i], Group = g, Pattern = pattern_names, Candidate = sum(j),
        Unknown = sum(is.na(main[j,i])), N = num, Share_of_unknown_in_group = nz(num, sum(t)),
        Share_of_group = nz(num, sum(j)), Weighted_share_of_group = nz(wz,sum(w[j])),
        Weighted_share_of_group_unknown = nz(wz,sum(w[j & is.na(main[,i])])),
        Share_of_pattern_across_groups = nz(num, as.integer(table(factor(a, levels = pattern_names)))),
        row.names = NULL)
    }))
  })) |>
    (\(d) { names(d)[2] <- variable; d })()
}
pattern_overall <- tab_pattern(rep("Overall",n))[,-2]
pattern_pathogen <- tab_pattern(iso$Organism_mapped, study_groups, "Organism")
pattern_sites <- tab_pattern(iso$Site, site_map$Site, "Site")
pattern_years <- tab_pattern(iso$Year, 2015:2024, "Year")
pattern_regions <- tab_pattern(iso$Region, regions, "Region")
pattern_overall$Top_organisms <- vapply(seq_len(nrow(pattern_overall)), function(i) {
  d <- pattern_pathogen[pattern_pathogen$Regimen == pattern_overall$Regimen[i] &
    pattern_pathogen$Pattern == pattern_overall$Pattern[i] & pattern_pathogen$N > 0, ]
  d <- utils::head(d[order(-d$N, d$Organism), ], 3)
  if (!nrow(d)) return("None")
  paste(paste0(d$Organism, " (", d$N, ")"), collapse = ", ")
}, "")
check(nrow(pat_records) == sum(table_a$Unknown) &&
  all(vapply(list(pattern_overall,pattern_pathogen,pattern_sites,pattern_years,pattern_regions),
    function(g) all(tapply(g$N,factor(g$Regimen, levels=regs),sum) == table_a$Unknown), TRUE)),
  "Three mutually exclusive missing patterns", "All Unknowns are exhausted by the three component states in the overall, pathogen, hospital, calendar-year and macro-region tables.")
check(all(near(pattern_overall$Share_of_unknown_in_group, pattern_overall$N / pattern_overall$Unknown),
  near(as.numeric(tapply(pattern_overall$Share_of_unknown_in_group, pattern_overall$Regimen, sum)), rep(1,6)),
  near(pattern_overall$Share_of_group, pattern_overall$N / n)),
  "Missing pattern shares", "Pattern percentages use their stated denominators and each regimen's conditional shares sum to 100%.")
staph <- c("CoNS", "S. aureus")
staph_rows <- iso$Organism_mapped %in% staph
target <- staph_rows & is.na(binary[,"Cefotaxime"])
ctx_regimens <- c("PEN + CTX","AMP + CTX")
ctx_contrib <- pattern_pathogen[pattern_pathogen$Regimen %in% ctx_regimens & pattern_pathogen$Organism %in% staph, ]
ctx_contrib$Share_of_all_regimen_unknown <- ctx_contrib$N/table_a$Unknown[match(ctx_contrib$Regimen, regs)]
ctx_contrib$Species_CTX_missing <- vapply(ctx_contrib$Organism, function(k) sum(target & iso$Organism_mapped == k), 0L)
asym <- bind(lapply(seq_len(nrow(pairs)), function(i) {
  a <- pairs$Regimen_A[i]
  b <- pairs$Regimen_B[i]
  ab <- p_index(i)
  include <- !complete.cases(ab)
  k <- data.frame(Contrast = pairs$Comparison[i], Relation = availability(ab[,1], ab[,2]),
    Organism = iso$Organism_mapped, A_status = encode(ab[,1]), B_status = encode(ab[,2]),
    A_unknown_pattern = ifelse(is.na(missing_patterns[,a]), "Evaluable", missing_patterns[,a]),
    B_unknown_pattern = ifelse(is.na(missing_patterns[,b]), "Evaluable", missing_patterns[,b]), W = w)
  k <- k[include, ]
  groups <- setdiff(names(k),"W")
  as.data.frame(dplyr::summarise(dplyr::group_by(k,dplyr::across(dplyr::all_of(groups))),
    N = dplyr::n(), Sum_weight = sum(W), .groups = "drop"))
}))
check(all(as.numeric(tapply(asym$N,asym$Contrast,sum)) == n - table_b$Both_evaluable),
  "Asymmetric evaluability audit", "All non-overlap records are reconciled by pathogen, the other regimen status and each applicable missing pattern.")

# Positivity and selection diagnostics. Condition on the observed status of
# the other component; preserve original vs filled testing provenance.
comp_status <- function(x) ifelse(is.na(x), "Missing", ifelse(x == 1L, "S", "I/R"))
observed_support <- function(other, focal, group, w0, raw_other, raw_focal) {
  a <- other[group]
  b <- focal[group]
  q <- w0[group]
  rr <- raw_other[group]
  ff <- raw_focal[group]
  known <- !is.na(b)
  c(N = sum(group), Focal_S = sum(b %in% 1L), Focal_nonS = sum(b %in% 0L), Focal_missing = sum(!known),
    Focal_missing_rate = nz(sum(!known), sum(group)),
    Weighted_focal_missing_rate = nz(sum(q[!known]),sum(q)),
    S_among_known = nz(sum(b %in% 1L),sum(known)),
    Weighted_S_among_known = nz(sum(q[b %in% 1L]),sum(q[known])),
    Focal_raw_observed = sum(!is.na(ff)), Both_raw_observed = sum(!is.na(rr) & !is.na(ff)),
    Missing_filled_by_rule = sum(!is.na(b) & is.na(ff)), Recorded_S_overridden = sum(ff %in% 1L & b %in% 0L))
}
support_drug <- function(i, direction, lev, labels) {
  col <- unlist(regimens[i,c("A","B")], use.names = FALSE)
  if (direction == "A given B") col <- rev(col)
  oc <- comp_status(binary[,col[1]])
  bind(lapply(labels, function(g) {
    keep <- lev == g & !fungus
    bind(lapply(c("S","I/R","Missing"), function(st) {
      idx <- keep & oc == st
      j <- observed_support(binary[,col[1]], binary[,col[2]], idx, w, raw_binary[,col[1]], raw_binary[,col[2]])
      sites_known <- length(unique(iso$Site[idx & !is.na(binary[,col[2]])]))
      years_known <- length(unique(iso$Year[idx & !is.na(binary[,col[2]])]))
      sites_missing <- unique(iso$Site[idx & is.na(binary[,col[2]])])
      site_overlap <- sum(sites_missing %in% iso$Site[idx & !is.na(binary[,col[2]])])
      data.frame(Regimen = regs[i], Group = g, Direction = direction, Focal = col[2], Other = col[1],
        Other_status = st, as.list(j), Sites_with_known = sites_known, Years_with_known = years_known,
        Sites_with_missing = length(sites_missing), Sites_with_both = site_overlap)
    }))
  }))
}
conditional <- bind(lapply(seq_len(nrow(regimens)), function(i) {
  bind(lapply(c("B given A", "A given B"), function(d) support_drug(i,d,iso$Organism_mapped,setdiff(study_groups,"Fungi"))))
}))
names(conditional)[names(conditional) == "Group"] <- "Organism"
ctx_species <- bind(lapply(which(regs %in% ctx_regimens), function(i) support_drug(i,"B given A",iso$Source_species_label,
  unique(iso$Source_species_label[staph_rows]))))
joint <- bind(lapply(seq_len(nrow(regimens)), function(i) {
  a <- regimens$A[i]
  b <- regimens$B[i]
  states <- expand.grid(A_status = c("S","I/R","Missing"), B_status = c("S","I/R","Missing"), stringsAsFactors = FALSE)
  m <- comp_status(binary[,a])
  k <- comp_status(binary[,b])
  rm <- comp_status(raw_binary[,a])
  rk <- comp_status(raw_binary[,b])
  bind(lapply(setdiff(study_groups,"Fungi"), function(g) {
    keep <- iso$Organism_mapped == g
    tally <- function(x,y) vapply(seq_len(nrow(states)),
      function(j) sum(keep & x == states$A_status[j] & y == states$B_status[j]),0L)
    data.frame(Regimen = regs[i], Organism = g, Component_A = a, Component_B = b, states,
      Recorded_n = tally(rm,rk), Analysis_n = tally(m,k))
  }))
}))
check(sum(conditional$N) == 2L*6L*sum(!fungus) &&
  all(rowSums(conditional[c("Focal_S","Focal_nonS","Focal_missing")]) == conditional$N) &&
  sum(joint$Recorded_n) == sum(joint$Analysis_n) && sum(joint$Analysis_n) == 6L*sum(!fungus),
  "Conditional drug support", "Observed and analysis AST availability reconciles in all non-fungal pathogen-by-companion-state cells. Unknown cells use no fitted probabilities.")
focus_support <- conditional[conditional$Organism %in% c("CoNS","S. aureus","E. coli") &
  conditional$Focal %in% c("Cefotaxime","Gentamicin"),]
completion_grid <- expand.grid(Shared = 0:1, Left = 0:1, Right = 0:1)
marginal_min <- function(y) ifelse(is.na(y), 0L, y)
marginal_max <- function(y) ifelse(is.na(y), 1L, y)
three <- function(i, mat = binary) {
  a <- mat[, unlist(pairs[i,c("Shared","Left","Right")]), drop = FALSE]
  colnames(a) <- c("Shared","Left","Right")
  a
}
limits <- function(x) {
  d <- matrix(NA_integer_, n, nrow(completion_grid))
  for (j in seq_len(nrow(completion_grid))) {
    vals <- unlist(completion_grid[j,])
    may <- rowSums(sweep(x, 2, vals, "!="), na.rm = TRUE) == 0L
    d[may,j] <- as.integer(max(vals[c(1,2)]) - max(vals[c(1,3)]))
  }
  plausible <- rowSums(!is.na(d))
  stopifnot(all(plausible > 0L))
  # Closed form from d=(1-shared)*(left-right) cross-checks enumeration.
  lower <- marginal_min(x[,2]) - marginal_max(x[,3])
  upper <- marginal_max(x[,2]) - marginal_min(x[,3])
  lo_formula <- ifelse(x[,1] %in% 1L, 0L, ifelse(is.na(x[,1]), pmin(lower,0L), lower))
  hi_formula <- ifelse(x[,1] %in% 1L, 0L, ifelse(is.na(x[,1]), pmax(upper,0L), upper))
  lo <- apply(d,1,min, na.rm = TRUE)
  hi <- apply(d,1,max, na.rm = TRUE)
  stopifnot(all(lo == lo_formula), all(hi == hi_formula))
  lo[fungus] <- hi[fungus] <- 0L
  plausible[fungus] <- 1L
  data.frame(Min = lo, Max = hi, Compatible_completions = plausible)
}
unit_bounds <- bind(lapply(seq_len(ncol(main)), function(i) {
  a <- main[,i]
  data.frame(Regimen = regs[i], N = n, Known_covered = sum(a %in% 1L), Unknown = sum(is.na(a)),
    Unit_lower = mean(marginal_min(a)), Unit_upper = mean(marginal_max(a)),
    S3_lower = wm(marginal_min(a)), S3_upper = wm(marginal_max(a)))
}))
b_check <- lapply(seq_len(nrow(pairs)), function(i) {
  a <- three(i)
  obs <- p_index(i)
  exact <- limits(a)
  crude_lo <- marginal_min(obs[,1]) - marginal_max(obs[,2])
  crude_hi <- marginal_max(obs[,1]) - marginal_min(obs[,2])
  stopifnot(all(exact$Min <= exact$Max), all(exact$Min >= crude_lo), all(exact$Max <= crude_hi),
    all(exact$Min[complete.cases(obs)] == obs[complete.cases(obs),1]-obs[complete.cases(obs),2]),
    all(exact$Min[complete.cases(obs)] == exact$Max[complete.cases(obs)]))
  row <- data.frame(Contrast = pairs$Comparison[i], Comparison = pair_label(i), Target_N = n,
    Main_denominator = table_b$Both_evaluable[i], Exact_uncertain_isolates = sum(exact$Min < exact$Max),
    Exact_unit_min_pp = 100*mean(exact$Min), Exact_unit_max_pp = 100*mean(exact$Max),
    Exact_S3_min_pp = 100*wm(exact$Min), Exact_S3_max_pp = 100*wm(exact$Max),
    Crude_unit_min_pp = 100*mean(crude_lo), Crude_unit_max_pp = 100*mean(crude_hi),
    Crude_S3_min_pp = 100*wm(crude_lo), Crude_S3_max_pp = 100*wm(crude_hi))
  row$S3_sign <- if (row$Exact_S3_min_pp > 0) "Positive throughout" else if (row$Exact_S3_max_pp < 0) "Negative throughout" else "Includes zero"
  list(summary = row, data = data.frame(Contrast = pairs$Comparison[i],
    iso[c("Isolate_ID","Source_row","Organism_mapped","S3_weight")],
    Shared_drug = pairs$Shared[i], Left_drug = pairs$Left[i], Right_drug = pairs$Right[i],
    Shared_AST = final_ast[,pairs$Shared[i]], Left_AST = final_ast[,pairs$Left[i]], Right_AST = final_ast[,pairs$Right[i]],
    A_status = encode(obs[,1]), B_status = encode(obs[,2]), Fungal_zero = fungus,
    Exact_min = exact$Min, Exact_max = exact$Max, Compatible_completions = exact$Compatible_completions,
    Crude_min = crude_lo, Crude_max = crude_hi))
})
table_d <- extract(b_check,"summary")
bound_records <- extract(b_check,"data")
# PSU-centered endpoint inference. Intervals are Taylor approximations with
# finite design degrees of freedom. Two Bonferroni one-sided limits cover the
# full identified set for a contrast. Eight limits cover all four sets jointly.
mean_design <- function(y, des, keep = all_isolates) {
  z <- subset(injected(des, Result = y), keep)
  fit <- survey::svymean(~Result,z)
  est <- as.numeric(coef(fit))
  var <- as.numeric(vcov(fit))
  wt <- weights(des)
  stopifnot(near(est,wm(y[keep],wt[keep])),near(var,manual_variance(y,keep,wt,est)))
  c(Estimate=est,SE=sqrt(var),DF=survey::degf(z))
}
endpoint_cols <- function(lo,hi,des,scale=100) {
  l <- mean_design(lo,des)
  u <- mean_design(hi,des)
  stopifnot(l[3]==u[3],l[1]<=u[1]+1e-12)
  t <- qt(.975,l[3])
  tf <- qt(1-.05/8,l[3])
  data.frame(DF=l[3], Lower_SE=scale*l[2], Upper_SE=scale*u[2],
    Lower_pointwise_low=scale*(l[1]-t*l[2]), Lower_pointwise_high=scale*(l[1]+t*l[2]),
    Upper_pointwise_low=scale*(u[1]-t*u[2]), Upper_pointwise_high=scale*(u[1]+t*u[2]),
    Joint95_low=scale*(l[1]-t*l[2]), Joint95_high=scale*(u[1]+t*u[2]),
    Family95_low=scale*(l[1]-tf*l[2]), Family95_high=scale*(u[1]+tf*u[2]),
    Joint95_strictly_positive=l[1]-t*l[2]>0, Family95_strictly_positive=l[1]-tf*l[2]>0, row.names=NULL)
}
bound_inference <- function(mat, label, count_ref) {
  bind(lapply(seq_len(nrow(pairs)),function(i) {
    lim <- limits(three(i,mat))
    bind(lapply(c("S3","Unit"),function(method) {
      des <- if (method=="S3") design_full else design_unit
      wt <- weights(des)
      data.frame(Analysis=label, Contrast=pairs$Comparison[i], Comparison=pair_label(i),
        Weighting=method, Target_N=n, Jointly_evaluable=count_ref$Both_evaluable[i],
        Exact_uncertain_isolates=sum(lim$Min < lim$Max), Sample_lower_pp=100*wm(lim$Min,wt),
        Sample_upper_pp=100*wm(lim$Max,wt),endpoint_cols(lim$Min,lim$Max,des),row.names=NULL)
    }))
  }))
}
infer_bounds <- bound_inference(binary,"11.15 main",table_b)
bound_range <- function(m,label) {
  bind(lapply(seq_along(regs), function(i) {
    z <- m[,i]
    data.frame(Analysis=label,Regimen=regs[i],N=n,Unit_lower=mean(marginal_min(z)),Unit_upper=mean(marginal_max(z)),
      S3_lower=wm(marginal_min(z)),S3_upper=wm(marginal_max(z)))
  }))
}
check(near(infer_bounds$Sample_lower_pp[infer_bounds$Weighting=="S3"],table_d$Exact_S3_min_pp) &&
  near(infer_bounds$Sample_upper_pp[infer_bounds$Weighting=="S3"],table_d$Exact_S3_max_pp) &&
  all(infer_bounds$Family95_low <= infer_bounds$Joint95_low & infer_bounds$Joint95_high <= infer_bounds$Family95_high),
  "Sharp endpoints and simultaneous limits", "All endpoints agree with closed-form enumeration and independent stratum-centred PSU variance. Contrast-wise 95% and four-contrast family limits are explicit.")
sensitivity <- lapply(list(
  list(label="A  Restore six low-level GEN/AMK rules",m=no6,mat=no6_binary,mask=mask6,ast=no6_ast),
  list(label="B  Roll back 24 species-scope isolate cells",m=no24,mat=no24_binary,mask=mask24,ast=no24_ast)),function(v) {
    a <- coverage_rows(v$m,prefix=v$label)[,-1]
    a$Analysis <- v$label
    a$Delta_weighted_pp <- 100*(a$Weighted_coverage-table_a$Weighted_coverage)
    a$Delta_evaluable <- a$Evaluable-table_a$Evaluable
    b <- fit_pairs(v$m,design_full,v$label)
    b$BH_q <- p.adjust(b$P_value,"BH")
    b$Delta_paired_pp <- b$Difference_pp-table_b$Difference_pp
    b$Delta_common_n <- b$Both_evaluable-table_b$Both_evaluable
    affected <- assignment[v$mask[rule_coords], ]
    affected$Restored_AST <- affected$Raw_AST
    affected$Analysis <- v$label
    stopifnot(identical(v$ast[!v$mask],final_ast[!v$mask]))
    list(cov=a,pair=b,regbounds=bound_range(v$m,v$label),bound=bound_inference(v$mat,v$label,b),cells=affected)
})
sens_cov <- extract(sensitivity,"cov")
sens_pairs <- extract(sensitivity,"pair")
sens_bounds <- extract(sensitivity,"bound")
sens_margins <- extract(sensitivity,"regbounds")
check(nrow(sens_cov)==12L && nrow(sens_pairs)==8L && nrow(sens_bounds)==16L &&
  nrow(sensitivity[[2]]$cells)==24L && nrow(sensitivity[[1]]$cells)==sum(mask6),
  "Two separate rule rollbacks", "Each alternative recalculates six coverages, the four paired domains and all logical bounds using only its own rollback.")

# Retain the published S10/S14 scopes and the Phase 2 trend model.
supplement_regimens <- rbind(regimens, data.frame(Regimen = c("PEN + AMK","AMP + AMK","VAN + IPM",
  "VAN + MEM","VAN + GEN","VAN + AMK"), A = c("Penicillin","Ampicillin",rep("Vancomycin",4)),
  B = c("Amikacin","Amikacin","Imipenem","Meropenem","Gentamicin","Amikacin"), Family = "Descriptive S10"))
supplement_regimens <- supplement_regimens[c(1:5,7,6,8:12), ]
s10_mat <- coverage(binary, supplement_regimens, TRUE)
region_rows <- setNames(lapply(regions, function(r) iso$Region == r), regions)
ga_rows <- setNames(lapply(ga, function(g) iso$GA %in% g), ga)
s10 <- bind(list(
  bind(lapply(names(region_rows), function(k) coverage_rows(s10_mat, region_rows[[k]], k, "S10 region"))),
  bind(lapply(names(ga_rows), function(k) coverage_rows(s10_mat, ga_rows[[k]], k, "S10 GA")))
))
s10$Domain_type <- rep(c("Macro-region", "Gestational age"), each=4L*ncol(s10_mat))
s10 <- s10[c("Domain_type",setdiff(names(s10),"Domain_type"))]
years <- bind(lapply(2015:2024, function(k) {
  d <- coverage_rows(main, iso$Year == k, k, "Annual")
  names(d)[1] <- "Year"
  d
}))
by_site <- bind(lapply(site_map$Site, function(k) {
  d <- coverage_rows(main, iso$Site == k, k, "Hospital", infer = FALSE)
  names(d)[1] <- "Site"
  d$Region <- unique(iso$Region[iso$Site == k])
  d[c("Site","Region","Regimen","Candidate","Covered","Uncovered","Unknown","Evaluable",
    "Unweighted_coverage","Weighted_coverage","Unknown_rate","Weighted_unknown_rate")]
}))
pathogen_results <- bind(lapply(study_groups, function(k) {
  d <- coverage_rows(main, iso$Organism_mapped == k,k, "Pathogen", infer = FALSE)
  names(d)[1] <- "Organism"
  d$Cohort_share <- d$Candidate/n
  d$Evaluable_share <- d$Evaluable/table_a$Evaluable[match(d$Regimen,regs)]
  d$Evaluable_minus_cohort_pp <- 100*(d$Evaluable_share-d$Cohort_share)
  d[c("Organism","Regimen","Candidate","Covered","Uncovered","Unknown","Evaluable","Unweighted_coverage",
    "Weighted_coverage","Unknown_rate","Weighted_unknown_rate","Cohort_share","Evaluable_share","Evaluable_minus_cohort_pp")]
}))
trend_domains <- c(list(Overall = all_isolates), region_rows, rev(ga_rows))
trend_one <- function(i,k,wt,des,y=main[,i]) {
  cand <- trend_domains[[k]]
  usable <- cand & !is.na(y)
  ny <- length(unique(iso$Year[usable]))
  label <- paste("S14",wt,k,regs[i],sep=" / ")
  ds <- domain(injected(des, Response=y),usable,label)
  df <- survey::degf(ds)
  dp <- df-1L
  header <- data.frame(Subgroup = k, Regimen = regs[i], Family = regimens$Family[i], Weighting = wt,
    Candidate = sum(cand), Evaluable = sum(usable), Unknown = sum(cand & is.na(y)), Years = ny,
    Design_df = df, Model_df = dp)
  result <- c(Beta_year=NA_real_, SE_year=NA_real_, PR_year=NA_real_, Year_low=NA_real_, Year_high=NA_real_,
    PR_5yr=NA_real_, Fiveyr_low=NA_real_, Fiveyr_high=NA_real_, t_statistic=NA_real_, P_value=NA_real_)
  problem <- if (sum(usable) < 10L || ny < 3L || length(unique(y[usable])) < 2L || dp <= 0L) "Sparse or constant response" else ""
  warnings <- character()
  if (!nzchar(problem)) {
    fit <- tryCatch(withCallingHandlers(survey::svyglm(Response ~ Year_centered, ds, family = quasipoisson(link = "log")),
      warning = function(e) {
        warnings <<- c(warnings, conditionMessage(e))
        invokeRestart("muffleWarning")
      }), error = function(e) {
        problem <<- conditionMessage(e)
        NULL
      })
    if (!is.null(fit) && !isTRUE(fit$converged)) problem <- "Model did not converge"
    if (!is.null(fit) && isTRUE(fit$converged)) {
      b <- coef(fit)[2]
      s <- sqrt(vcov(fit)[2,2])
      if (!is.finite(b) || !is.finite(s) || s <= 0) {
        problem <- "Non-estimable coefficient or variance"
      } else {
        tt <- b/s
        lim <- b+c(-1,1)*qt(.975,dp)*s
        pv <- 2*pt(-abs(tt),dp)
        stopifnot(near(summary(fit,df.resid=dp)$coefficients[2,], c(b,s,tt,pv)),
          (pv < .05) == (lim[1] > 0 || lim[2] < 0))
        result[] <- c(b,s,exp(b),exp(lim),exp(5*b),exp(5*lim),tt,pv)
      }
    }
  }
  data.frame(header, as.list(result), Status = if (nzchar(problem)) problem else "Estimated",
    Warnings = paste(unique(warnings),collapse = " | "))
}
trend <- bind(lapply(c("S3","Unit"), function(wt) {
  des <- if (wt == "S3") design_full else design_unit
  bind(lapply(names(trend_domains), function(k) bind(lapply(seq_along(regs), function(i) trend_one(i,k,wt,des)))))
}))
trend$BH_q <- ave(trend$P_value, trend$Weighting,trend$Family,
  FUN = function(ps) p.adjust(ps, method="BH", n=length(ps)))
trend$Family_tests <- ave(rep(1L,nrow(trend)), trend$Weighting,trend$Family,FUN=sum)
trend_families <- unique(trend[c("Weighting","Family","Family_tests")])
trend_families$Estimated <- vapply(seq_len(nrow(trend_families)), function(i) {
  g <- trend$Family == trend_families$Family[i] & trend$Weighting == trend_families$Weighting[i]
  sum(is.finite(trend$P_value[g]))
},0L)
trend_families$BH_below_005 <- vapply(seq_len(nrow(trend_families)), function(i) {
  g <- trend$Family == trend_families$Family[i] & trend$Weighting == trend_families$Weighting[i]
  sum(trend$BH_q[g] < .05,na.rm=TRUE)
},0L)
check(nrow(s10) == 96L && nrow(years) == 60L && nrow(trend) == 108L &&
  all(tapply(years$Candidate,years$Regimen,sum) == n) &&
  identical(as.integer(trend_families$Family_tests),c(36L,18L,36L,18L)) &&
  all(trend$Evaluable + trend$Unknown == trend$Candidate),
  "S10 S14 and annual scope", "All 12 original S10 regimens, four region and four GA rows per regimen. S14 has the original 54 models per weighting and 36/18 BH families.")
annual_domains <- bind(lapply(names(trend_domains)[-1],function(k) {
  bind(lapply(2015:2024,function(yr) {
    a <- coverage_rows(main,trend_domains[[k]] & iso$Year==yr,k,"Annual domain",detail=FALSE)
    a$Year <- yr
    a
  }))
}))
annual_all <- bind(list(data.frame(Subgroup="Overall",years),annual_domains))
idx <- seq_len(23L)
drug_pairs <- which(upper.tri(matrix(0,23,23),diag=TRUE), arr.ind=TRUE)
s2 <- bind(lapply(seq_len(nrow(drug_pairs)), function(i) {
  a <- drug_pairs[i,1]
  b <- drug_pairs[i,2]
  z <- if (a == b) binary[,a] else combine(binary[,a],binary[,b])
  if (a != b) z[fungus] <- 0L
  data.frame(Drug_A = drug_fields[a], Drug_B = drug_fields[b], Candidate = n, Covered = sum(z %in% 1L),
    Uncovered = sum(z %in% 0L), Unknown = sum(is.na(z)), Evaluable = sum(!is.na(z)),
    Unit_coverage = mean(z,na.rm=TRUE), S3_coverage = wm(z[!is.na(z)],w[!is.na(z)]), A_index = a, B_index = b)
}))
s2_main <- vapply(seq_along(regs), function(i) {
  v <- s2$Drug_A == regimens$A[i] & s2$Drug_B == regimens$B[i] |
    s2$Drug_B == regimens$A[i] & s2$Drug_A == regimens$B[i]
  stopifnot(sum(v)==1L)
  s2$S3_coverage[v]
},0)
check(nrow(s2) == choose(24,2) && near(s2_main,table_a$Weighted_coverage),
  "S2 heatmap", "All 23 singles and 253 distinct pairs recalculated; the six core cells match Table A.")
s2_ab <- c("PEN","AMP","PIP","SAM","PIS","TZP","VAN","TEC","LZD","CFZ","CAZ","CTX",
  "FEP","CSL","IPM","MEM","GEN","AMK","CIP","LVX","AZM","ERY","ATM")
color_coverage <- function(x) {
  ramp <- grDevices::colorRamp(c("#ECF5F4","#174E67"),space="Lab")
  result <- rep("#F1F2F2",length(x))
  valid <- !is.na(x)
  z <- ramp(x[valid])
  result[valid] <- grDevices::rgb(z[,1],z[,2],z[,3],maxColorValue=255)
  result
}
s2$Unit_fill <- color_coverage(s2$Unit_coverage)
s2$S3_fill <- color_coverage(s2$S3_coverage)
s2_legend <- data.frame(Coverage=c(0,.2,.4,.6,.8,1,NA), Fill=color_coverage(c(0,.2,.4,.6,.8,1,NA)))
# Exploratory paired domains use the original full design and identical pair definitions.
subgroup_pairs <- bind(lapply(names(c(region_rows,ga_rows)), function(k) {
  eligible <- c(region_rows,ga_rows)[[k]]
  bind(lapply(seq_len(nrow(pairs)), function(i) {
    ab <- p_index(i)
    keep <- eligible & complete.cases(ab)
    d <- ab[,1]-ab[,2]
    psus <- length(unique(iso$Site[keep]))
    strata <- length(unique(iso$Region[keep]))
    df <- psus-strata
    cross <- table(factor(ab[keep,1],levels=0:1), factor(ab[keep,2],levels=0:1))
    estimate <- nz(sum(w[keep]*d[keep]),sum(w[keep]))
    se <- lo <- hi <- tt <- pv <- NA_real_
    reason <- if (!any(keep)) "No common observations" else if (sum(keep)<10L) "Fewer than 10 common isolates" else
      if (df<=0L) "Too few PSUs" else if (sum(cross[2,1],cross[1,2])==0L) "No discordant pairs" else ""
    if (!nzchar(reason)) {
      des <- domain(injected(design_full,Diff=d),keep,paste("Exploratory",k,i))
      stopifnot(survey::degf(des)==df)
      fit <- survey::svymean(~Diff,des)
      se <- sqrt(as.numeric(vcov(fit)))
      stopifnot(near(coef(fit)[1],estimate),near(se^2,manual_variance(d,keep,w,estimate)))
      if (is.finite(se) && se>0) {
        interval <- estimate+c(-1,1)*qt(.975,df)*se
        lo <- interval[1]; hi <- interval[2]
        tt <- estimate/se; pv <- 2*pt(-abs(tt),df)
      } else reason <- "Variance not estimable"
    }
    data.frame(Subgroup=k, Domain_type=if(k %in% regions) "Macro-region" else "Gestational age",
      Contrast=pairs$Comparison[i],Comparison=pair_label(i),Candidate=sum(eligible),Both_evaluable=sum(keep),
      A0_B0=cross[1,1],A0_B1=cross[1,2],A1_B0=cross[2,1],A1_B1=cross[2,2],
      A_coverage=nz(sum(w[keep]*ab[keep,1]),sum(w[keep])),
      B_coverage=nz(sum(w[keep]*ab[keep,2]),sum(w[keep])),
      Difference_pp=100*estimate,SE_pp=100*se,CI_low_pp=100*lo,CI_high_pp=100*hi,
      Contributing_PSUs=psus,DF=df,t_statistic=tt,P_value=pv,
      Status=ifelse(nzchar(reason),"NE","Estimated"),Reason=reason)
  }))
}))
subgroup_pairs$BH_q <- p.adjust(subgroup_pairs$P_value,method="BH",n=32)
check(nrow(subgroup_pairs)==32 && all(rowSums(subgroup_pairs[c("A0_B0","A0_B1","A1_B0","A1_B1")])==subgroup_pairs$Both_evaluable),
      "32 paired subgroup contrasts","All four factorial counts reconcile; BH family comprises 32 planned exploratory comparisons.")

# Refit the same S14 model for each currently supported trend and both full-domain completions.
selected_trends <- trend[!is.na(trend$BH_q) & trend$BH_q<.05,]
scenarios <- bind(lapply(seq_len(nrow(selected_trends)),function(j) {
  base <- selected_trends[j,]
  des <- if (base$Weighting=="S3") design_full else design_unit
  i <- match(base$Regimen,regs)
  bind(lapply(c("Evaluable","Unknown all Covered","Unknown all Uncovered"),function(scenario) {
    y <- main[,i]
    if (scenario!="Evaluable") y[is.na(y)] <- as.integer(scenario=="Unknown all Covered")
    fit <- trend_one(i,base$Subgroup,base$Weighting,des,y)
    fit$Scenario <- scenario
    fit$Main_BH_q <- if (scenario=="Evaluable") base$BH_q else NA_real_
    fit$Shift_PR_year <- fit$PR_year-base$PR_year
    fit$Shift_PR_5yr <- fit$PR_5yr-base$PR_5yr
    fit$Direction_retained <- sign(fit$Beta_year)==sign(base$Beta_year)
    fit$Nominal_P_below_05 <- fit$P_value<.05
    fit
  }))
}))
if (!nrow(scenarios)) scenarios <- data.frame(Status="No S14 trend has BH q < 0.05")
check(nrow(selected_trends)==0L || nrow(scenarios)==3L*nrow(selected_trends),
      "S14 selected trend scenarios","Selection is based on the current release BH results; completion fits report nominal P and CI.")

# S19 disposition definitions are taken from the original clinical missingness method.
outcome_levels <- c("Survived to discharge","In-hospital death","Withdrawal of care")
outcome_rule <- setNames(c(outcome_levels,NA_character_),c("否/否","是/是","否/是","不详/不详"))
outcome_pair <- paste(trimws(linked[,87]),trimws(linked[,88]),sep="/")
check(all(outcome_pair %in% names(outcome_rule)),
      "S19 outcome coding","Joint coded death and death-or-withdrawal status; unknown stays out of outcome test domains.")
outcome <- factor(unname(outcome_rule[outcome_pair]),levels=outcome_levels)
s19_link <- data.frame(iso[c("Isolate_ID","Patient_ID","Source_row","Continuation","Site","Region","Year","S3_weight")],
  Death_record=linked[,87],Death_or_withdrawal_record=linked[,88],Discharge_outcome=outcome)
s19_groups <- list(Year=factor(iso$Year,levels=2015:2024),Region=factor(iso$Region,levels=regions),Outcome=outcome)
s19_units <- bind(lapply(seq_along(regs),function(i) {
  unknown <- is.na(main[,i])
  bind(lapply(names(s19_groups),function(dimn) {
    grp <- s19_groups[[dimn]]
    bind(lapply(levels(grp),function(g) {
      keep <- !is.na(grp) & grp==g
      data.frame(Regimen=regs[i],Dimension=dimn,Level=g,Unknown=sum(unknown&keep),Denominator=sum(keep),
        Unit_unknown_rate=nz(sum(unknown&keep),sum(keep)),
        Weighted_unknown_rate=nz(sum(w[unknown&keep]),sum(w[keep])),PSUs=length(unique(iso$Site[keep])))
    }))
  }))
}))
s19_test <- function(i,dimn,wt) {
  grp <- s19_groups[[dimn]]
  keep <- !is.na(grp)
  indicator <- factor(is.na(main[,i]),levels=c(FALSE,TRUE))
  des <- if(wt=="S3") design_full else design_unit
  dom <- domain(injected(des,Unknown_status=indicator,Group=grp),keep,paste("S19",wt,dimn,regs[i]))
  cells <- table(grp[keep],indicator[keep])
  est <- s19_units[s19_units$Regimen==regs[i] & s19_units$Dimension==dimn,]
  prop <- if(wt=="S3") est$Weighted_unknown_rate else est$Unit_unknown_rate
  warnings <- character()
  reason <- ""
  f <- d1 <- d2 <- p <- NA_real_
  if (any(rowSums(cells)==0) || any(colSums(cells)==0) || survey::degf(dom)<=0)
    reason <- "Empty margin or insufficient design df"
  if (!nzchar(reason)) {
    fit <- tryCatch(withCallingHandlers(survey::svychisq(~Group+Unknown_status,dom,statistic="F"),
      warning=function(e) { warnings <<- c(warnings,conditionMessage(e)); invokeRestart("muffleWarning") }),
      error=function(e) { reason <<- conditionMessage(e); NULL })
    if(!is.null(fit)) {
      values <- unname(c(fit$statistic,fit$parameter,fit$p.value))
      if(length(values)==4L && all(is.finite(values)) && all(values[2:3]>0)) {
        f <- values[1]; d1 <- values[2]; d2 <- values[3]; p <- values[4]
        stopifnot(near(p,pf(f,d1,d2,lower.tail=FALSE)))
      } else reason <- "Non-finite Rao-Scott result"
    }
  }
  data.frame(Regimen=regs[i],Dimension=dimn,Weighting=wt,Denominator=sum(keep),Excluded=sum(!keep),
    Unknown_n=sum(keep & is.na(main[,i])),Design_df=survey::degf(dom),Min_rate=min(prop),Max_rate=max(prop),
    F=f,Numerator_df=d1,Denominator_df=d2,P_value=p,Min_observed_cell=min(cells),
    Small_cells_below_5=sum(cells<5),Zero_cells=sum(cells==0),
    Status=ifelse(nzchar(reason),"NE","Estimated"),Reason=reason,Warnings=paste(warnings,collapse=" | "))
}
s19_tests <- bind(lapply(seq_along(regs),function(i) {
  bind(lapply(names(s19_groups),function(d) bind(lapply(c("S3","Unit"),function(wt)s19_test(i,d,wt)))))
}))
s19_tests$BH_q <- ave(s19_tests$P_value,s19_tests$Dimension,s19_tests$Weighting,
    FUN=function(p) p.adjust(p,"BH",n=6L))
check(nrow(s19_tests)==36L && nrow(s19_units)==102L &&
      all(s19_tests$Denominator+s19_tests$Excluded==n),
      "S19 dimensions and denominators","18 S3 and 18 unit-weight tests. Each dimension has a separate six-regimen BH family.")
s19_docs <- data.frame(Regimen=regs,Unknown=table_a$Unknown,Denominator=n,
    Unit_unknown_rate=table_a$Unknown_rate,Weighted_unknown_rate=table_a$Weighted_unknown_rate)



s3_data <- years[years$Regimen %in% regs[1:4],]
options(OutDec=".")
write.csv(s3_data,file.path(out,"S3_Fig_Data.csv"),row.names=FALSE,na="",fileEncoding="UTF-8")

# Outputs
iso_core <- iso[setdiff(names(iso),c("Unit_weight","Year_centered"))]
status_frame <- function(m) as.data.frame(lapply(as.data.frame(m,check.names=FALSE),encode),check.names=FALSE)
status_labeled <- data.frame(iso_core,status_frame(main),check.names=FALSE)
raw_columns <- c("CSV ID","CSV Hospital","CSV Year","CSV Stratum","CSV GA weeks","CSV GA category","CSV organism group","CSV Gram")
raw_source_meta <- as.data.frame(csv[,c(1:4,8L,6L,44:45)],check.names=FALSE)
names(raw_source_meta) <- raw_columns
original_data <- data.frame(Isolate_ID=iso$Isolate_ID,Source_row=iso$Source_row,raw_source_meta,raw_ast,check.names=FALSE)
final_data <- data.frame(iso_core,final_ast,check.names=FALSE)
linked_data <- data.frame(Isolate_ID=iso$Isolate_ID,Source_row=iso$Source_row,Site=iso$Site,Region=iso$Region,
  S3_weight=w,linked[,patient_cols],check.names=FALSE)
all_sir <- c("S","I","R","Missing")
truth16 <- expand.grid(A=all_sir,B=all_sir,stringsAsFactors=FALSE)
binary_state <- setNames(c(1L,0L,0L,NA_integer_),all_sir)
truth16$Expected <- ifelse(truth16$A=="S"|truth16$B=="S","Covered",
  ifelse(truth16$A %in% c("I","R") & truth16$B %in% c("I","R"),"Uncovered","Unknown"))
truth16$Observed <- encode(combine(binary_state[truth16$A],binary_state[truth16$B]))
check(identical(truth16$Expected,truth16$Observed),"All 16 recorded state combinations",
  "Distinguish S, I, R and missing in both drugs, including R + missing and I + missing.")
special <- status_labeled[fungus|unidentified,c("Isolate_ID","Source_row","Organism_mapped","Site","Year",regs)]
special$Decision <- ifelse(fungus[fungus|unidentified],"Structural regimen zero. Raw single-agent AST retained.",
  "Use observed AST. Any unresolved regimen status is Unknown.")
design_summary <- data.frame(Region=c(regions,"Overall"),Hospitals=c(as.integer(site_n),30),
  Isolates=c(as.integer(table(factor(iso$Region,levels=regions))),n),Design_df=c(as.integer(site_n)-1L,26))
psu_domains <- bind(psu_qc)
boundary_rows <- bind(list(s10,data.frame(Subgroup=as.character(years$Year),years)))
boundary_rows <- boundary_rows[grepl("unavailable",boundary_rows$CI_method) | grepl("unavailable",boundary_rows$Unit_CI_method),]
check(identical(vapply(paths,sha,""),hashes),
  "Input SHA-256 on completion","All declared authorised and public input files remain unchanged")
record("Coverage method scope","Class summaries and phenotype screens do not supply missing single-agent AST. Table 3 / S6 and noncoverage CoNS analyses are outside this run.")
def_table <- function(...) {
  a <- list(...)
  data.frame(Item=names(a),Definition=unlist(a,use.names=FALSE),row.names=NULL)
}
defs <- def_table(
  "Source rows"="Source_row is the 1-based data record excluding header. CSV_record includes the header as record 1.",
  "Strict S"="After the approved intrinsic resistance determination, at least one S is Covered, both I/R is Uncovered, and no S with any missing component is Unknown.",
  "Fungi"="Structural Uncovered for combinations only. Single-agent and class AST remain recorded values.",
  "Design"="Construct full design from 1,196 isolates. Hospital PSU, macro-region stratum and region-year S3 weights; subset that design for domains.",
  "Finite degrees of freedom"="Domain PSUs minus represented strata. S14 residual degrees of freedom subtract one for the continuous year slope.",
  "CIs"="Proportion intervals use t-logit Taylor linearization. Boundary or no-variance cells have explicitly unavailable intervals.",
  "Multiplicities"="Four paired contrasts share one BH family in each analysis. S14 has four separate test families: China 36 and WHO 18 for each of S3 and unit weights.",
  "Coverage limits"="Sample full-cohort lower covered/all; upper (covered + unknown)/all, with the same weights as the target.",
  "S10 scope"="Four macro-regions, four gestational age groups. Six core regimens plus the six pre-existing descriptive rows.",
  "S14 scope"="Continuous year quasi-Poisson log link; annual and 5-year PR with t tests. Year-specific coverage uses the evaluable domain for that year.")
defs_bounds <- def_table(
  "Sharp bound"="Enumerate all compatible binary states of the three distinct drugs per isolate. Shared drug takes the same state in both regimens. Full cohort n=1196.",
  "Pointwise t intervals"="The variance of each endpoint uses the full design; subtract or add t(0.975, df) times endpoint SE.",
  "Joint95"="A per-contrast conservative 95% joint interval for the entire population identification set. Each endpoint has a one-sided 2.5% error allowance.",
  "Family95"="For each weighting and scenario, eight one-sided tests across all four contrasts use 0.05/8 each.",
  "Tipping root uncertainty"="Invert the simultaneous band: values of p for which zero remains compatible. A sample crossing can be absent even when zero is sampling-compatible.",
  "Local domain"="Fixed set of jointly determined records at both p=0 and p=1. Included as a diagnostic with no full-cohort interpretation.",
  "Conditional joint support"="Counts use the companion state before imputation. Genuine recorded AST and rule-determined values are tabulated separately.")
heatmap_notes <- def_table(
  "Figure title"="S2 Fig. In vitro coverage of individual agents and two-agent regimens.",
  "Matrix"="Same Phase 2 order of all 23 agents. Upper triangle and diagonal use unweighted proportions. Lower triangle uses S3-weighted proportions.",
  "Denominator"="Percentage among evaluable isolates. Companion worksheet shows corresponding unweighted evaluable n in both triangles.",
  "Diagonals"="Single-drug susceptibility among known recorded or approved assigned results. Weighted single-drug estimates appear in Backing_data.",
  "Colour scale"="Fixed fill reproduces the Phase 2 Lab gradient from #ECF5F4 at 0% to #174E67 at 100%. Grey is not estimable.",
  "Strict S"="Any S for a combination is Covered. Two I/R are Uncovered. Other incomplete results are Unknown. Fungi are structural Uncovered for combinations.")
books <- list(
  Data_Audit=list(Input_hashes=inputs,Continuation_7=link_summary,Inherited_field_log=link_detail,
    Definitions=def_table("Scope"="Current-input and restoration checks only. Historical release comparisons are archived separately.")),
  FinalData=list(Isolate_coverage=status_labeled,Analysis_AST=final_data,
    Original_AST=original_data,Linked_patient_fields=linked_data,Matrix_assignments=assignment,
    Approved_24_assignments=assignment[assignment$Approved_24_cell,],
    S_I_conflicts=assignment[assignment$Disposition %in% c("S to R","I to R"),],
    Matrix_83_rules=rule_summary,Rule_release_counts=rule_totals,Special_organisms=special,
    Continuation_7=link_summary,Inputs=inputs),
  Results=list(Six_regimens=table_a,Four_paired_contrasts=table_b,
    Paired_coverage_2x2=pair_coverage,Paired_evaluability_2x2=pair_evaluability,
    Paired_pathogen_mix=pair_mix,Composition_shift=pair_mix_metrics,
    Overall_regimen_bounds=unit_bounds,S10_region_GA=s10,S14_trends=trend,
    S14_selected_trend_sensitivity=scenarios,By_hospital=by_site,By_pathogen=pathogen_results,
    Definitions=defs),
  Missing_Bounds=list(Exact_paired_endpoints=infer_bounds,Exact_and_crude_comparison=table_d,
    Regimen_bounds=unit_bounds,Unknown_overall=pattern_overall,
    Unknown_pathogen=pattern_pathogen,Unknown_hospital=pattern_sites,
    Unknown_region=pattern_regions,Unknown_year=pattern_years,Unknown_Staph_CTX=ctx_contrib,
    Unknown_isolate_log=pat_records,Asymmetric_evaluability=asym,
    Conditional_AST_support=conditional,CTX_species_support=ctx_species,Focused_AST_support=focus_support,
    Joint_AST_3x3=joint,Bounds_by_isolate=bound_records,
    Definitions=defs_bounds),
  Sensitivity=list(Coverage_rollbacks=sens_cov,Paired_rollbacks=sens_pairs,Paired_bounds_rollbacks=sens_bounds,
    Regimen_bounds_rollbacks=sens_margins,Restore_6_matrix_cells=sensitivity[[1]]$cells,
    Restore_24_isolate_cells=sensitivity[[2]]$cells,Complete_three_components=pair_panel,
    Unit_weight_pairs=pair_unit,Main_bounds=infer_bounds,Definitions=def_table(
      "A"="Restore GEN/AMK raw values only for GBS, Other Streptococcus spp. and Enterococcus spp.",
      "B"="Roll back only the 24 specified species-scope cells to their 11.15 raw values.",
      "Independent alternatives"="A and B are each applied separately to the 11.15 main analysis. Report-only phenotype screens are not changes.",
      "Pair q-values"="BH adjustment across the four tests within each sensitivity analysis.")),
  Subgroup_Paired=list(Exploratory_32=subgroup_pairs,Definitions=def_table(
    "Family"="One exploratory BH family of 32 comparisons, including NE slots.",
    "Inference"="S3-weighted within-isolate difference. Complete design domains with contributing PSU minus represented stratum df.",
    "NE"="No common isolates, fewer than ten common isolates, too few PSUs, no discordance or variance unavailable. No interaction test.")),
  S19_Coverage_Missingness_Detail=list(Panel_E_overall=s19_docs,Group_rates=s19_units,
    Rao_Scott_tests=s19_tests,Outcome_linkage=s19_link,Annual_coverage_overall=years,Annual_coverage_all_domains=annual_all,
    Definitions=def_table("Unknown"="Strict-S combination is indeterminate. A known S makes a regimen evaluable.",
    "S19 design"="Second-order Rao-Scott F. Hospital PSU, macro-region strata. Separate S3 and unit weight analyses.",
    "Missing outcome"="Exclude only from outcome groups and tests. Join to the infant of the preceding row for each of seven continuations.",
    "Multiplicity"="Six p values per dimension and weighting; BH. Small observed cells flagged. Failed design tests retained as NE.",
    "S14"="Coverage among evaluable isolates, unlike S19 Unknown proportion among all eligible isolates.")),
  S2_Heatmap=list(Backing_data=s2,Colour_legend=s2_legend,
    Drug_key=data.frame(Order=seq_len(23),Abbreviation=s2_ab,Agent=drug_fields[1:23]),
    Title_and_legend=heatmap_notes),
  QC=list(Assertions=bind(qc),Paired_independent_Taylor=bind(crosscheck),
    Bound_endpoint_crosscheck=infer_bounds,
    PSU_by_stratum=design_summary,S3_weight_summary=weight_distribution,S3_calibration_40_cells=weight_cells,
    Site_Region_lookup=site_lookup,Design_domains=psu_domains,Inference_boundaries=boundary_rows,
    Strict_S_16_cases=truth16,Trend_BH_families=trend_families,Input_hashes=inputs))
names(books) <- ifelse(names(books)=="S19_Coverage_Missingness_Detail",names(books),paste0("CANEOS_Coverage_11.15.5_",names(books)))
make_heatmap <- function(n_values=FALSE) {
  m <- matrix(NA_real_,23,23,dimnames=list(s2_ab,s2_ab))
  m[cbind(s2$A_index,s2$B_index)] <- if(n_values) s2$Evaluable else s2$Unit_coverage
  off <- s2$A_index!=s2$B_index
  m[cbind(s2$B_index[off],s2$A_index[off])] <- if(n_values) s2$Evaluable[off] else s2$S3_coverage[off]
  data.frame(Agent=s2_ab,m,check.names=FALSE,row.names=NULL)
}
heat <- "CANEOS_Coverage_11.15.5_S2_Heatmap"
books[[heat]] <- c(list(Coverage_heatmap=make_heatmap(),Evaluable_n=make_heatmap(TRUE)),books[[heat]])
stopifnot(all(vapply(books,function(b) all(nchar(names(b))<=31L),TRUE)))
capture.output(sessionInfo(),file=file.path(out,"R_sessionInfo.txt"))
# Only the subgroup workbook is required by the formal output map. All other
# coverage workbooks may contain restricted row-level audit information.
export_books <- if (isTRUE(cfg$export_internal_audit)) books else
  books["CANEOS_Coverage_11.15.5_Subgroup_Paired"]
for (k in names(export_books)) writexl::write_xlsx(
  export_books[[k]], file.path(out,paste0(k,".xlsx")), format_headers=TRUE)
readr::write_csv(bind(qc),file.path(out,"Current_Input_Checks.csv"))
print(table_a[c("Regimen","Covered","Uncovered","Unknown","Weighted_coverage","Weighted_CI_low","Weighted_CI_high")],row.names=FALSE)
print(table_b[c("Contrast","Both_evaluable","Difference_pp","CI_low_pp","CI_high_pp","BH_q")],row.names=FALSE)
print(infer_bounds[infer_bounds$Weighting=="S3",c("Contrast","Sample_lower_pp","Sample_upper_pp","Joint95_low","Family95_low")],row.names=FALSE)


# ------------------------------------------------------------------------------
# Write clean final-manuscript CSVs & downstream plotting inputs
# ------------------------------------------------------------------------------
write_csv_utf8(s10, file.path(out, "S9_Table_Empirical_Coverage.csv"), full_precision = TRUE)
write_csv_utf8(table_b, file.path(out, "S10_Table_Paired_Coverage.csv"), full_precision = TRUE)
write_csv_utf8(infer_bounds, file.path(out, "S10_Table_Paired_Logical_Bounds.csv"), full_precision = TRUE)
write_csv_utf8(trend, file.path(out, "S14_Table_Coverage_Trends.csv"), full_precision = TRUE)
write_csv_utf8(table_a, file.path(out, "Figure2_Six_Regimens.csv"), full_precision = TRUE)

# Export Figure 2 Panel A / Panel B plotting table consumed by R/07_figures/figure2
fig2_order <- c("AMP + GEN", "PEN + GEN", "AMP + CAZ", "AMP + CTX", "PEN + CTX", "PEN + CAZ")
fa <- table_a[match(fig2_order, table_a$Regimen), , drop = FALSE]
fb_order <- c("C1", "C2", "C3", "C4")
fb <- table_b[match(fb_order, table_b$Contrast), , drop = FALSE]
fb_b <- infer_bounds[infer_bounds$Weighting == "S3", , drop = FALSE]
fb_b <- fb_b[match(fb_order, fb_b$Contrast), , drop = FALSE]
stopifnot(!anyNA(fa$Regimen), !anyNA(fb$Contrast), !anyNA(fb_b$Contrast))

plot_a <- data.frame(
  Panel = "A", Order = 1:6, Label = gsub(" ", "", fa$Regimen, fixed = TRUE),
  Candidate_N = fa$Candidate, Covered_n = fa$Covered, Uncovered_n = fa$Uncovered,
  Unknown_n = fa$Unknown, Evaluable_n = fa$Evaluable,
  Covered_share = fa$Covered / fa$Candidate, Uncovered_share = fa$Uncovered / fa$Candidate,
  Unknown_share = fa$Unknown / fa$Candidate,
  Weighted_coverage = fa$Weighted_coverage, Coverage_CI_low = fa$Weighted_CI_low,
  Coverage_CI_high = fa$Weighted_CI_high, stringsAsFactors = FALSE
)
plot_b <- data.frame(
  Panel = "B", Order = 1:4, Label = gsub(" + ", "+", fb$Comparison, fixed = TRUE),
  Candidate_N = fb_b$Target_N, Paired_n = fb$Both_evaluable,
  Paired_A_coverage = fb$A_coverage, Paired_B_coverage = fb$B_coverage,
  Difference_pp = fb$Difference_pp, Difference_CI_low_pp = fb$CI_low_pp,
  Difference_CI_high_pp = fb$CI_high_pp,
  Sample_lower_pp = fb_b$Sample_lower_pp, Sample_upper_pp = fb_b$Sample_upper_pp,
  Joint_lower_pp = fb_b$Joint95_low, Simultaneous_lower_pp = fb_b$Family95_low,
  stringsAsFactors = FALSE
)
all_cols <- union(names(plot_a), names(plot_b))
for (col in setdiff(all_cols, names(plot_a))) plot_a[[col]] <- NA_real_
for (col in setdiff(all_cols, names(plot_b))) plot_b[[col]] <- NA_real_
fig2_plot_data <- rbind(plot_a[all_cols], plot_b[all_cols])
write_csv_utf8(fig2_plot_data, file.path(out, "Figure2_Plot_Data.csv"), full_precision = TRUE)

# Export S2 Heatmap backing data (consumed by 03_export_s2_excel.R)
s2_heatmap_rates <- make_heatmap(FALSE)
s2_heatmap_n     <- make_heatmap(TRUE)
s2_drug_key      <- data.frame(Order = seq_len(23), Abbreviation = s2_ab, Antibiotic = drug_fields[1:23], stringsAsFactors = FALSE)
write_csv_utf8(s2_heatmap_rates, file.path(out, "S2_Heatmap_Rates.csv"), full_precision = TRUE)
write_csv_utf8(s2_heatmap_n,     file.path(out, "S2_Heatmap_Evaluable_N.csv"), full_precision = FALSE)
write_csv_utf8(s2_drug_key,      file.path(out, "S2_Heatmap_Drug_Key.csv"), full_precision = FALSE)

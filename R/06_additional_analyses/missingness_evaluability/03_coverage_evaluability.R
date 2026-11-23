# Extracted strict-S functions and the 24 applicable cells of the approved
# 83-cell coverage matrix. No AST-reportability matrix is substituted here.
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
classify_regimen_evaluability <- function(dat, root) {
  mapping <- static_cfg$get_organism_groups()
  all_rules <- static_cfg$get_coverage_intrinsic_rules()
  five_drugs <- c("Penicillin", "Ampicillin", "Ceftazidime", "Cefotaxime", "Gentamicin")
  rules_sub <- all_rules[all_rules$Antibiotic %in% five_drugs, , drop = FALSE]
  drug_columns <- c(Penicillin = 46L, Ampicillin = 47L, Ceftazidime = 56L,
    Cefotaxime = 57L, Gentamicin = 62L)
  rules <- data.frame(Organism = rules_sub$Organism, Drug = rules_sub$Antibiotic,
    Column = unname(drug_columns[rules_sub$Antibiotic]), Assignment = rules_sub$Matrix_value,
    stringsAsFactors = FALSE)
  organism <- mapping$Organism[match(dat$raw[[44]], mapping$Raw_group)]
  stopifnot(!anyNA(organism), nrow(rules) == 24L,
    !anyDuplicated(rules[c("Organism", "Drug")]), all(rules$Assignment == "R"))
  drug_columns <- c(Penicillin = 46L, Ampicillin = 47L, Ceftazidime = 56L,
    Cefotaxime = 57L, Gentamicin = 62L)
  ast <- as.matrix(dat$raw[, drug_columns])
  colnames(ast) <- names(drug_columns)
  stopifnot(all(ast %in% c("", "敏感", "中介", "耐药", "S", "I", "R")),
    all(unname(drug_columns[rules$Drug]) == rules$Column))
  for (j in seq_len(nrow(rules))) ast[organism == rules$Organism[j], rules$Drug[j]] <- "R"
  regimens <- data.frame(Key = c("PEN+CAZ", "PEN+CTX", "AMP+CAZ", "AMP+CTX", "PEN+GEN", "AMP+GEN"),
    A = rep(c("Penicillin", "Ampicillin", "Penicillin", "Ampicillin"), c(2, 2, 1, 1)),
    B = c("Ceftazidime", "Cefotaxime", "Ceftazidime", "Cefotaxime", "Gentamicin", "Gentamicin"))
  binary <- activity(ast)
  status <- sapply(seq_len(nrow(regimens)), function(i) combine(binary[, regimens$A[i]], binary[, regimens$B[i]]))
  colnames(status) <- regimens$Key
  # Fungi remain structural Uncovered, as in the approved coverage module.
  status[organism == "Fungi", ] <- 0L
  weights <- static_cfg$get_calibration_weights()
  names(weights)[names(weights) == "Isolates"] <- "SampleN"
  names(weights)[names(weights) == "Births"]   <- "PopulationN"
  key <- paste(weights$Region, weights$Year)
  j <- match(paste(dat$context$region, dat$context$year), key)
  sample_n <- tabulate(j, nbins = nrow(weights))
  fraction <- function(x) x / ave(x, weights$Year, FUN = sum)
  stopifnot(nrow(weights) == 40L, !anyDuplicated(key), !anyNA(j),
    identical(sample_n, as.integer(weights$SampleN)), all(weights$Weight > 0),
    max(abs(weights$Weight - fraction(weights$PopulationN) / fraction(sample_n))) < 1e-8)
  list(status = status, weight = weights$Weight[j], regimens = regimens,
    organism = organism, rules = rules, mapping = mapping, calibration = weights)
}
analyse_coverage <- function(dat, root) {
  classified <- classify_regimen_evaluability(dat, root)
  status <- classified$status; weight <- classified$weight; d <- dat$context
  counts <- tests <- list()
  for (k in colnames(status)) {
    y <- as.integer(is.na(status[, k]))
    axes <- c(list(Overall = factor(rep("Overall", nrow(d)))), groups(d))
    for (axis in names(axes)) {
      g <- axes[[axis]]
      counts[[length(counts) + 1L]] <- bind(lapply(levels(g), function(level) {
        keep <- !is.na(g) & g == level; a <- status[keep, k]; w <- weight[keep]
        C <- sum(a %in% 1L); U <- sum(a %in% 0L); M <- sum(is.na(a)); N <- sum(keep)
        data.frame(Key = k, Dimension = axis, Level = level, Records = N,
          Covered = C, Uncovered = U, Unknown = M, Evaluable = C + U,
          UnitPercent = 100 * M / N, WeightedPercent = 100 * sum(w[is.na(a)]) / sum(w),
          PSUs = length(unique(d$hospital[keep])))
      }))
    }
    for (axis in names(groups(d))) for (wt in c("S3", "Unit")) {
      tests[[length(tests) + 1L]] <- rs_test(d, rep(TRUE, nrow(d)), y,
        groups(d)[[axis]], if (wt == "S3") weight else rep(1, nrow(d)),
        "Coverage", k, axis)
    }
  }
  counts <- bind(counts)
  list(overall = counts[counts$Dimension == "Overall", ], counts = counts,
    tests = adjust_families(bind(tests)), classified = classified)
}

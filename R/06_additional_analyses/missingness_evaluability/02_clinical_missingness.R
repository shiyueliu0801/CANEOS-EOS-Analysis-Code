classify_clinical_missingness <- function(dat) {
  raw <- dat$raw; primary <- dat$primary; dictionary <- dat$dictionary
  v <- function(j) raw[[j]][primary]
  late_blank <- blank(v(80))
  early_death <- v(83) == "是"
  early_withdrawal <- v(83) == "否" & v(84) == "是"
  stopifnot(identical(late_blank, blank(v(82))),
    all(late_blank == (early_death | early_withdrawal)),
    all(v(20) != "无" | v(21) == "0"),
    all(v(20) != "有" | v(21) %in% c("", "1", "2", "3")),
    all(!is_missing(v(20)) | blank(v(21))),
    all(v(18) != "否" | (v(17) == "0" & v(19) == "否")))
  states <- lapply(seq_len(nrow(dictionary)), function(i) {
    info <- dictionary[i, ]
    z <- if (info$Rule == "isolate") raw[[info$Column]] else v(info$Column)
    miss <- is_missing(z) | (info$Rule == "isolate" & z == "菌不详")
    a <- ifelse(miss, "Missing", "Observed")
    a[tolower(z) %in% c("not applicable", "不适用")] <- "Not applicable"
    a[tolower(z) %in% c("not evaluable", "不可评价")] <- "Not evaluable"
    observed <- z[a == "Observed"]
    if (info$Type == "numeric" && any(!is.finite(suppressWarnings(as.numeric(observed))))) {
      stop("Unrecognised numeric code: ", info$Key)
    }
    if (info$Type == "regimen" && !all(tolower(observed) %in%
      c("3gc+pen", "bsp", "carb", "2gc", "pen", "3gc", "1gc", "carb+gly/lzd", "other", "gly/lzd"))) {
      stop("Unrecognised regimen code: ", info$Key)
    }
    if (info$Type == "change" && !all(tolower(observed) %in% c("change", "unchange", "unchanged"))) {
      stop("Unrecognised modification code: ", info$Key)
    }
    if (info$Rule == "msaf_positive") {
      a[v(20) == "无"] <- "Not applicable"
      a[is_missing(v(20))] <- "Uncertain"
    }
    if (info$Rule == "late") {
      a[late_blank & early_death] <- "Not evaluable"
      a[late_blank & early_withdrawal] <- "Uncertain"
    }
    a
  })
  names(states) <- dictionary$Key
  states
}
analyse_clinical <- function(dat) {
  states <- classify_clinical_missingness(dat)
  counts <- tests <- list()
  for (i in seq_len(nrow(dat$dictionary))) {
    info <- dat$dictionary[i, ]; st <- states[[i]]
    d <- if (info$Unit == "Isolate") dat$context else dat$infants
    axes <- c(list(Overall = factor(rep("Overall", nrow(d)))), groups(d))
    for (axis in names(axes)) {
      if (axis == "Outcome" && info$OutcomeCircular) next
      g <- axes[[axis]]
      counts[[length(counts) + 1L]] <- bind(lapply(levels(g), function(level) {
        keep <- which(g == level); a <- st[keep]
        n <- sum(a == "Missing"); N <- sum(a %in% c("Missing", "Observed"))
        data.frame(Key = info$Key, Variable = info$Variable, Unit = info$Unit,
          Dimension = axis, Level = level, Records = length(keep), Missing = n, Eligible = N,
          Percent = if (N > 0) 100 * n / N else NA_real_, Observed = sum(a == "Observed"),
          NotApplicable = sum(a == "Not applicable"), NotEvaluable = sum(a == "Not evaluable"),
          Uncertain = sum(a == "Uncertain"))
      }))
    }
    for (axis in names(groups(d))) {
      tests[[length(tests) + 1L]] <- rs_test(d, st %in% c("Missing", "Observed"),
        as.integer(st == "Missing"), groups(d)[[axis]], rep(1, nrow(d)),
        "Clinical", info$Key, axis, axis == "Outcome" && info$OutcomeCircular)
    }
  }
  counts <- bind(counts)
  list(overall = counts[counts$Dimension == "Overall", ], counts = counts,
    tests = adjust_families(bind(tests)), states = states)
}

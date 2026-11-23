# The binary coding and denominators follow the accepted clinical cohort,
# antimicrobial_use and manuscript_results scripts and the S19 family plan.
values <- function(field) raw[[field]][primary]
categorise <- function(x, code) {
  x <- canon(x)
  suspected <- if (code %in% c("IVH", "NEC")) c("疑诊", "可疑") else character()
  if (any(!is.na(x) & !x %in% c("否", "是", suspected))) stop("Unexpected outcome code: ", code)
  ifelse(is.na(x), NA_integer_, ifelse(x == "否", 0L, ifelse(x == "是", 1L, 2L)))
}
states <- setNames(lapply(names(source_fields), function(k) categorise(values(source_fields[[k]]), k)), names(source_fields))
outcome <- lapply(states, function(x) as.integer(x == 1L))
for (window in list(c("D72","C72"), c("D7","C7"), c("Death","Composite"))) {
  death <- outcome[[window[1L]]]
  composite <- outcome[[window[2L]]]
  stopifnot(identical(is.na(death), is.na(composite)), all(death <= composite, na.rm = TRUE))
}
outcome$Withdrawal <- as.integer(outcome$Death == 0L & outcome$Composite == 1L)
outcome$IVH_all <- as.integer(states$IVH > 0L)
outcome$NEC_all <- as.integer(states$NEC > 0L)
outcome_labels <- c(D72 = "Death within 72 h of birth", D7 = "Death within 7 days of birth",
  Death = "All-cause in-hospital death", C72 = "Death or withdrawal of treatment within 72 h of birth",
  C7 = "Death or withdrawal of treatment within 7 days of birth", Withdrawal = "Withdrawal of treatment overall",
  Composite = "Death or withdrawal of treatment overall", BPD = "BPD", Shock = "Shock",
  LOS = "In-hospital late-onset sepsis", ROP = "ROP", IVH = "Confirmed IVH", NEC = "Confirmed NEC",
  IVH_all = "Confirmed or suspected IVH", NEC_all = "Confirmed or suspected NEC")
disposition <- names(outcome_labels)[seq_len(7L)]

regimen_codes <- c("PEN","BSP","1GC","2GC","3GC","3GC+PEN","CARB","GLY/LZD","CARB+GLY/LZD","Other")
regimen_labels <- c("Narrow-spectrum penicillin", "Broad-spectrum penicillin", "First-generation cephalosporin",
  "Second-generation cephalosporin", "Third-generation cephalosporin",
  "Third-generation cephalosporin plus a penicillin-class agent", "Carbapenem", "Glycopeptide or linezolid",
  "Carbapenem plus glycopeptide or linezolid", "Other regimens")
normalize_regimen <- function(x) {
  x <- clean(x)
  x[is.na(x) | x == "" | toupper(x) %in% c("UNK","UNKNOWN") | x == "不详"] <- "UNK"
  if (any(!x %in% c(regimen_codes, "UNK"))) stop("Unexpected antimicrobial regimen code.")
  x
}
normalize_change <- function(x) {
  x <- tolower(clean(x))
  x[x %in% c("", "unk", "unknown", "不详")] <- NA_character_
  if (any(!is.na(x) & !x %in% c("change","unchange"))) stop("Unexpected modification code.")
  x
}
plans <- setNames(lapply(regimen_fields, function(f) normalize_regimen(values(f))), c("Initial","Early","Late"))
change <- lapply(modification_fields, function(f) normalize_change(values(f)))
not_assessable <- !is.na(outcome$C72) & outcome$C72 == 1L
late_raw <- clean(values(regimen_fields[3L]))
stopifnot(identical(is.na(change[[2L]]), not_assessable),
  identical(is.na(late_raw) | late_raw == "", not_assessable))
plans$Late[not_assessable] <- NA_character_
classic <- clean(values(classic_field))
stopifnot(!anyNA(classic), all(classic %in% c("", "经典方案")),
  all(classic != "经典方案" | plans$Initial == "3GC+PEN"))
modified <- change[[1L]] %in% "change" | change[[2L]] %in% "change"
determinate <- modified | (change[[1L]] %in% "unchange" & change[[2L]] %in% "unchange")
any_regimen <- function(codes) Reduce(`|`, lapply(plans, function(v) v %in% codes))
third <- c("3GC","3GC+PEN")
carb <- c("CARB","CARB+GLY/LZD")
gly <- c("GLY/LZD","CARB+GLY/LZD")

# All 40 tests are needed for the accepted BH results even when a row
# contributes only to the multiplicity calculation.
treatment <- list(ClassicFlag = classic == "经典方案",
  Initial_3GC_any = plans$Initial %in% third, Any_3GC = any_regimen(third),
  Change_early = change[[1L]] == "change", Change_late = change[[2L]] == "change",
  Any_change = ifelse(determinate, modified, NA), Any_CARB = any_regimen(carb),
  Any_GLY = any_regimen(gly), Any_broad = any_regimen(union(carb,gly)),
  Prophylaxis = categorise(values(prophylaxis_field), "Prophylaxis"))
treatment_labels <- setNames(c("Classic empirical regimen", "Initial third-generation cephalosporin regimen",
  "Third-generation cephalosporin category in any recorded window", "Modification within 72 h",
  "Modification during 72 h–7 days", "At least one modification among infants with determinable status",
  "Carbapenem category in any recorded window", "Glycopeptide or linezolid category in any recorded window",
  "Carbapenem, glycopeptide or linezolid category in any recorded window", "Antifungal prophylaxis"), names(treatment))
for (stage in names(plans)) {
  for (i in seq_along(regimen_codes)) {
    key <- paste(stage, regimen_codes[i], sep = "_")
    treatment[[key]] <- plans[[stage]] == regimen_codes[i]
    treatment_labels[key] <- regimen_labels[i]
  }
}
stopifnot(length(treatment) == 40L, identical(names(treatment_labels), names(treatment)))
selected <- c("3GC", "3GC+PEN", "CARB", "GLY/LZD", "CARB+GLY/LZD")
display_groups <- list(
  Summary = c("ClassicFlag","Any_3GC","Change_early","Change_late","Any_change","Any_CARB","Any_GLY"),
  Initial = paste0("Initial_", regimen_codes),
  Early = paste0("Early_", selected), Late = paste0("Late_", selected),
  Additional = "Prophylaxis")
display_codes <- unlist(display_groups, use.names = FALSE)
treatment_tests <- adjust_bh(bind(lapply(names(treatment), function(k) {
  data.frame(Code = k, Indicator = treatment_labels[[k]], rao_scott(as.integer(treatment[[k]])))
})), rep("Treatment (40)", length(treatment)))
treatment_tests$Display <- ifelse(treatment_tests$Code %in% setdiff(display_codes, "Prophylaxis"), "Word",
  ifelse(treatment_tests$Code == "Prophylaxis", "Requested additional indicator", "BH family only"))

infant_results <- function(y, labels, codes, tests) {
  bind(lapply(codes, function(k) {
    vals <- as.integer(y[[k]])
    props <- lapply(infant_domains, function(keep) proportion(infant_design, vals, keep, clinical = TRUE))
    stopifnot(props$Overall$Events == props$CoNS$Events + props$NonCoNS$Events,
      props$Overall$N == props$CoNS$N + props$NonCoNS$N)
    test <- tests[match(k, tests$Code), ]
    data.frame(Code = k, Indicator = labels[[k]], wide(props),
      test[setdiff(names(test), c("Code","Indicator"))],
      Change_pp = 100*(props$NonCoNS$Estimate - props$Overall$Estimate), row.names = NULL)
  }))
}
treatments <- infant_results(treatment, treatment_labels, display_codes, treatment_tests)
treatments$Block <- rep(names(display_groups), lengths(display_groups))
outcome_tests <- adjust_bh(bind(lapply(names(outcome_labels), function(k) {
  data.frame(Code = k, rao_scott(outcome[[k]]))
})), ifelse(names(outcome_labels) %in% c("IVH_all","NEC_all"), "Suspected diagnosis (2)", "Primary outcomes (13)"))
outcomes <- infant_results(outcome, outcome_labels, names(outcome_labels), outcome_tests)

# Reconstruct the infant-isolate relationship from the locked source CSV.
# Only anonymous keys and necessary fields are retained in memory.
trim_value <- function(x) trimws(as.character(x), whitespace = "[\\h\\v]")
is_unknown <- function(x) is.na(x) | tolower(trim_value(x)) %in%
  c("", "na", "n/a", "unknown", "unk", "不详", "未知", "缺失")
canonical <- function(x) {
  x <- trim_value(x)
  x[is_unknown(x)] <- NA_character_
  x
}
require_true <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}
read_source <- function(path) {
  x <- readr::read_csv(path, col_types = readr::cols(.default = readr::col_character()),
    locale = readr::locale(encoding = "UTF-8"), na = character(),
    name_repair = "minimal", show_col_types = FALSE)
  require_true(nrow(readr::problems(x)) == 0L, "Malformed source CSV.")
  x
}
regions <- c("Municipality", "Eastern", "Central", "Western")
patterns <- c("Both S", "One S, one R", "Both R")
endpoint_labels <- c(Mod0 = "Regimen modification within 72 h",
  Carb0 = "Carbapenem-containing regimen at 0-72 h",
  Mod7 = "Regimen modification at 72 h-7 d",
  Carb7 = "Carbapenem-containing regimen at 72 h-7 d",
  Death = "In-hospital death", Composite = "In-hospital death or withdrawal of care")
endpoint_fields <- c(Mod0 = "72h内是否更改", Carb0 = "Use within 72 h",
  Mod7 = "72h-7d是否更改", Carb7 = "Use within 72 h–7 days",
  Death = "死亡", Composite = "死亡或放弃治疗")
regimen_codes <- c("3GC+PEN", "BSP", "CARB", "CARB+GLY/LZD", "3GC", "2GC", "PEN", "1GC", "GLY/LZD", "Other", "UNK")
regimen_value <- function(x) {
  x <- trim_value(x)
  x[is_unknown(x)] <- "UNK"
  require_true(all(x %in% regimen_codes), "Unexpected treatment category.")
  x
}
change_value <- function(x) {
  x <- tolower(canonical(x))
  require_true(all(na.omit(x) %in% c("change", "unchange")), "Unexpected treatment-modification code.")
  as.integer(x == "change")
}
outcome_value <- function(x) {
  x <- canonical(x)
  require_true(all(na.omit(x) %in% c("是", "否")), "Unexpected clinical-outcome code.")
  as.integer(x == "是")
}
carb_value <- function(x) ifelse(x == "UNK", NA_integer_, as.integer(x %in% c("CARB", "CARB+GLY/LZD")))

prepare_data <- function(data_file, rules_df) {
  raw <- read_source(data_file)
  keys <- c("ID", "医院名称", "入院时间", "分层")
  taxon_field <- "菌群整理归类（归类到属）"
  required <- c(keys, "Gram", taxon_field, "是否经典方案", "First use after Birth",
    "生后72小时内是否死亡", "生后72小时内是否死亡或放弃", unname(endpoint_fields),
    "头孢他啶（三代）", "头孢噻肟（三代）")
  require_true(all(required %in% names(raw)), "Required source fields are missing.")
  taxon_columns <- which(names(raw) == taxon_field)
  require_true(length(taxon_columns) == 2L && identical(trim_value(raw[[taxon_columns[1L]]]),
    trim_value(raw[[taxon_columns[2L]]])), "Duplicated taxon columns disagree.")
  primary <- which(!is_unknown(raw$ID))
  continuation <- which(is_unknown(raw$ID))
  require_true(length(primary) > 0L && primary[1L] == 1L && !anyDuplicated(trim_value(raw$ID[primary])),
    "Infant primary keys are missing or duplicated.")
  require_true(all(vapply(raw[continuation, keys], function(x) all(is_unknown(x)), TRUE)),
    "Continuation rows contain partially populated linkage keys.")
  anchor <- cummax(ifelse(!is_unknown(raw$ID), seq_len(nrow(raw)), 0L))
  # The source stores extra isolates immediately after the infant's primary row.
  # Validate every populated shared clinical field before restoring its keys.
  shared <- c(seq.int(5L, taxon_columns[1L] - 1L),
              seq.int(match("是否经典方案", names(raw)), ncol(raw)))
  require_true(all(vapply(shared, function(j) {
    current <- canonical(raw[[j]][continuation])
    parent <- canonical(raw[[j]][anchor[continuation]])
    all(is.na(current) | (!is.na(parent) & current == parent))
  }, TRUE)), "A continuation row has conflicting or additional infant-level information.")
  linked <- as.data.frame(lapply(raw[keys], trim_value), check.names = FALSE)
  linked[continuation, keys] <- linked[anchor[continuation], keys]
  require_true(!anyNA(canonical(unlist(linked))) &&
    nrow(unique(linked)) == length(primary), "Recovered infant keys are inconsistent.")
  region_map <- c("一层" = "Municipality", "二层" = "Eastern", "三层" = "Central", "四层" = "Western")
  region <- unname(region_map[linked$分层])
  hospital <- linked$医院名称
  year <- suppressWarnings(as.integer(linked$入院时间))
  require_true(!anyNA(region) && all(!is.na(year) & year %in% 2015:2024) &&
    all(vapply(split(region, hospital), function(x) length(unique(x)) == 1L, TRUE)),
    "Hospital, macro-region or year coding is invalid.")
  # Hospital names themselves define the same clusters as the validated site map.
  # Factor labels are never exported. Public omission labels use input-order indices.
  infants <- data.frame(ID = linked$ID[primary],
    HospitalID = factor(hospital[primary], levels = unique(hospital)),
    Region = factor(region[primary], levels = regions), Year = year[primary], w_unit = 1)
  counts <- table(linked$ID)
  infants$Isolates <- as.integer(counts[infants$ID])
  infants$Initial <- regimen_value(raw[["First use after Birth"]][primary])
  infants$Reg0 <- regimen_value(raw[[endpoint_fields[["Carb0"]]]][primary])
  infants$Reg7 <- regimen_value(raw[[endpoint_fields[["Carb7"]]]][primary])
  infants$Mod0 <- change_value(raw[[endpoint_fields[["Mod0"]]]][primary])
  infants$Mod7 <- change_value(raw[[endpoint_fields[["Mod7"]]]][primary])
  infants$Carb0 <- carb_value(infants$Reg0)
  infants$Carb7 <- carb_value(infants$Reg7)
  infants$Death <- outcome_value(raw[[endpoint_fields[["Death"]]]][primary])
  infants$Composite <- outcome_value(raw[[endpoint_fields[["Composite"]]]][primary])
  infants$Death72 <- outcome_value(raw[["生后72小时内是否死亡"]][primary])
  infants$Composite72 <- outcome_value(raw[["生后72小时内是否死亡或放弃"]][primary])
  guideline <- canonical(raw[["是否经典方案"]][primary])
  require_true(all(na.omit(guideline) == "经典方案"), "Unexpected guideline-regimen flag.")
  infants$Guideline <- guideline %in% "经典方案"
  require_true(all(infants$Initial[infants$Guideline] == "3GC+PEN") &&
    identical(is.na(infants$Death), is.na(infants$Composite)) &&
    !any(infants$Death > infants$Composite, na.rm = TRUE) &&
    !any(infants$Death72 > infants$Composite72, na.rm = TRUE), "Outcome or nested-regimen consistency failed.")

  rules <- rules_df; rules$Applicable <- ifelse(as.logical(rules$Applicable), "TRUE", "FALSE")
  rule_key <- paste(rules$Taxon, rules$Drug)
  require_true(!anyDuplicated(rule_key) && all(rules$Applicable %in% c("TRUE", "FALSE")) &&
    all(c("CAZ", "CTX") %in% rules$Drug), "Invalid AST applicability metadata.")
  orgs <- unique(rules$Taxon[rules$Applicable == "TRUE"])
  require_true(all(vapply(orgs, function(x) all(paste(x, c("CAZ", "CTX")) %in% rule_key), TRUE)),
    "The applicability metadata do not cover both agents for every included taxon.")
  isolates <- data.frame(ID = linked$ID, HospitalID = hospital, Region = region, Year = year,
    Taxon = trim_value(raw[[taxon_columns[1L]]]), CAZraw = trim_value(raw[["头孢他啶（三代）"]]),
    CTXraw = trim_value(raw[["头孢噻肟（三代）"]]))
  ast_codes <- c("敏感" = "S", "中介" = "I", "耐药" = "R", S = "S", I = "I", R = "R")
  for (agent in c("CAZ", "CTX")) {
    index <- match(paste(isolates$Taxon, agent), rule_key)
    applicable <- !is.na(index) & rules$Applicable[index] %in% "TRUE"
    isolates[[paste0(agent, "Applicable")]] <- applicable
    isolates[[agent]] <- ifelse(applicable, unname(ast_codes[canonical(isolates[[paste0(agent, "raw")]])]), NA_character_)
  }
  started_3gc <- infants$Initial %in% c("3GC", "3GC+PEN")
  any_ent <- infants$ID %in% isolates$ID[isolates$Taxon %in% orgs]
  infants$Eligible_Enterobacterales <- started_3gc & any_ent & infants$Isolates == 1L
  infants$Eligible_Ecoli <- infants$Eligible_Enterobacterales &
    infants$ID %in% isolates$ID[isolates$Taxon == "大肠埃希菌"]
  # Every isolate is retained. Only monomicrobial infants receive a paired phenotype.
  single <- match(infants$ID, isolates$ID)
  paired <- paste0(isolates$CAZ[single], isolates$CTX[single])
  ast <- unname(c(SS = patterns[1L], SR = patterns[2L], RS = patterns[2L], RR = patterns[3L])[paired])
  infants$Complete_AST <- infants$Eligible_Enterobacterales & !is.na(ast)
  infants$TestGroup <- factor(ifelse(infants$Complete_AST, ast, NA_character_), levels = patterns)
  infants$AST_exclusion <- ifelse(infants$Eligible_Enterobacterales & !infants$Complete_AST,
    ifelse(isolates$CAZ[single] %in% "I" | isolates$CTX[single] %in% "I", "Contains I", "Unavailable or uninterpretable paired AST"), NA_character_)
  for (k in names(endpoint_labels)) {
    infants[[paste0("OK_", k)]] <- !is.na(infants[[k]])
    infants[[paste0("Outcome_", k)]] <- factor(infants[[k]], levels = 0:1, labels = c("No", "Yes"))
  }
  early <- infants$Composite72 %in% 1L
  late_status <- function(raw_value, value, late = FALSE, regimen = FALSE) {
    explicit <- toupper(trim_value(raw_value)) %in% c("UNK", "UNKNOWN")
    status <- ifelse(is.na(value), "Missing/uninterpretable record", "Evaluable")
    status[is.na(value) & explicit] <- if (regimen) "Regimen recorded as UNK" else "Modification flag unknown"
    blank_late <- late & early & is.na(value) & !explicit
    status[blank_late & infants$Death72 %in% 1L] <- "Early death; late window blank"
    status[blank_late & infants$Death72 %in% 0L] <- "Early withdrawal; late window blank"
    status[blank_late & is.na(infants$Death72)] <- "Early composite; component unknown; late window blank"
    status
  }
  for (k in names(endpoint_labels)) {
    infants[[paste0("State_", k)]] <- late_status(raw[[endpoint_fields[[k]]]][primary], infants[[k]],
      late = k %in% c("Mod7", "Carb7"), regimen = k %in% c("Carb0", "Carb7"))
  }
  list(infants = infants, isolates = isolates, raw = raw, primary = primary,
    continuation = continuation, rules = rules, rules_sha256 = "static_config")
}

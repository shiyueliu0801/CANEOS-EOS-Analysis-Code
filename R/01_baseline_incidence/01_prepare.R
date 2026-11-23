# ==============================================================================
# Module 01: Data Preparation (Baseline Characteristics & Hospital-Based Incidence)
# ==============================================================================
# Preserves all verified cohort construction, continuation-row validation, and
# denominator checks from the locked revision package while using unified
# logical input paths and static methodological column definitions.
# ==============================================================================

case_columns <- static_cfg$get_baseline_case_columns()$Source_header

# Build one infant cohort and one isolate linkage view from the cohort CSV
raw <- read_csv_safe(cfg$input_files$cohort_data, encoding = "UTF-8")
require_data(ncol(raw) >= length(case_columns) &&
               identical(names(raw)[seq_along(case_columns)], case_columns),
             "Invalid cohort_data schema. See docs/INPUT_SCHEMA.md.")

organism <- case_columns[44L]
clinical_columns <- 5:43
keys <- case_columns[1:4]
current <- extract(raw, c(keys, organism))
primary <- which(!is_missing(current[[1L]]))
continuation <- which(is_missing(current[[1L]]))

require_data(length(primary) > 0 && !anyDuplicated(current[[1L]][primary]) &&
               all(continuation - 1L %in% primary) &&
               all(vapply(current[continuation, keys, drop = FALSE], function(x) all(is_missing(x)), logical(1))) &&
               !any(is_missing(current[[organism]])),
             "Unlinked case row, non-unique infant key, or missing isolate identity.")

for (j in clinical_columns) {
  a <- canon(raw[[j]][continuation])
  b <- canon(raw[[j]][continuation - 1L])
  require_data(all(is.na(a) | (!is.na(b) & a == b)),
               paste("Continuation row differs from its preceding infant in", names(raw)[j]))
}
current[continuation, keys] <- current[continuation - 1L, keys]

if (ncol(raw) >= 76L && names(raw)[76L] == organism) {
  require_data(identical(clean(raw[[44L]]), clean(raw[[76L]])),
               "Duplicated isolate classification fields disagree.")
}
require_data(!anyDuplicated(current[c(keys[1L], organism)]),
             "Identical organism repeated for one infant.")

year <- suppressWarnings(as.integer(current[[3L]]))
reg <- unname(region_map[current[[4L]]])
require_data(all(year %in% year_set) && !anyNA(reg),
             "Case years or macro-region codes are invalid.")

# Join hospital annual admissions with the controlled anonymised site map
read_admissions <- function() {
  adm <- read_csv_safe(cfg$input_files$hospital_admissions, encoding = "auto",
                       expected_cols = c(keys[2L], as.character(year_set)))
  sm <- read_csv_safe(cfg$input_files$site_map, encoding = "auto")

  # Support both site_map formats: (S4 center, S4 stratum, Hospital) or (HospitalID, Region, CaseKey)
  if (all(c("S4 center", "S4 stratum", "Hospital") %in% names(sm))) {
    sm_rows <- sm[nzchar(clean(sm$Hospital)) & grepl("^Site [0-9]+$", clean(sm[["S4 center"]])), , drop = FALSE]
    site_df <- data.frame(
      HospitalID = clean(sm_rows[["S4 center"]]),
      CaseKey    = clean(sm_rows[["Hospital"]]),
      Region     = clean(sm_rows[["S4 stratum"]]),
      stringsAsFactors = FALSE
    )
  } else if (all(c("HospitalID", "CaseKey", "Region") %in% names(sm))) {
    site_df <- data.frame(
      HospitalID = clean(sm$HospitalID),
      CaseKey    = clean(sm$CaseKey),
      Region     = clean(sm$Region),
      stringsAsFactors = FALSE
    )
  } else {
    stop("Unsupported site_map schema. See docs/INPUT_SCHEMA.md.", call. = FALSE)
  }

  idx_num <- suppressWarnings(as.integer(sub("^Site ", "", site_df$HospitalID)))
  require_data(!anyNA(idx_num) && identical(sort(idx_num), seq_len(nrow(site_df))),
               "HospitalID must be unique sequential Site identifiers.")
  site_df <- site_df[order(idx_num), , drop = FALSE]

  adm_keys <- clean(adm[[keys[2L]]])
  mi <- match(site_df$CaseKey, adm_keys)
  require_data(nrow(adm) == nrow(site_df) && !anyNA(mi) && !anyDuplicated(mi),
               "hospital_admissions and site_map do not form a 1-to-1 hospital crosswalk.")

  names_year <- as.character(year_set)
  denom <- vapply(names_year, function(y) as.numeric(clean(adm[[y]][mi])), numeric(nrow(site_df)))
  dimnames(denom) <- list(NULL, names_year)
  list(site = site_df, denominators = denom)
}

ad <- read_admissions()
site_map <- ad$site
den <- ad$denominators
site_map$Region <- factor(site_map$Region, levels = region_names)
require_data(!anyNA(site_map) && !any(is_missing(site_map$CaseKey)) && !anyDuplicated(site_map$CaseKey) &&
               all(is.finite(den)) && all(den > 0) && all(den == trunc(den)) &&
               all(table(site_map$Region) >= 2L),
             "Admissions require positive integer denominators, unique keys and >=2 hospitals per stratum.")

site_map$Admissions <- rowSums(den)
site_idx <- match(current[[2L]], site_map$CaseKey)
require_data(!anyNA(site_idx) && all(as.character(site_map$Region[site_idx]) == reg),
             "A case hospital is absent from the admissions input or has an inconsistent macro-region.")

inf <- data.frame(ID = current[[1L]][primary], HospitalID = site_map$HospitalID[site_idx[primary]],
                  Region = site_map$Region[site_idx[primary]], Year = year[primary],
                  t = year[primary] - min(year_set), w = 1)
iso <- data.frame(ID = current[[1L]], HospitalID = site_map$HospitalID[site_idx],
                  Region = site_map$Region[site_idx], Year = year, Organism = current[[organism]])
inf$Ecoli <- as.integer(inf$ID %in% iso$ID[iso$Organism == "大肠埃希菌"])
inf$GBS <- as.integer(inf$ID %in% iso$ID[iso$Organism == "GBS"])
inf$EOS <- 1L

v <- function(j) canon(raw[[j]][primary])
num <- function(j) {
  x <- v(j)
  z <- suppressWarnings(as.numeric(x))
  if (any(!is.na(x) & !is.finite(z))) stop("Invalid numeric field: ", case_columns[j])
  z
}
binary <- function(j, positive = "是", negative = "否") {
  x <- v(j)
  if (!all(na.omit(x) %in% c(positive, negative))) stop("Invalid binary field: ", case_columns[j])
  ifelse(is.na(x), NA_integer_, as.integer(x == positive))
}

msaf_binary <- binary(20L, "有", "无")
msaf_grade <- num(21L)
graded <- !is.na(msaf_grade)
require_data(all(na.omit(msaf_grade) %in% 0:3) &&
               all(msaf_binary[graded] == as.integer(msaf_grade[graded] > 0)),
             "Known MSAF grades conflict with the recorded presence indicator.")
positive_grade <- ifelse(msaf_binary == 1L & graded, msaf_grade, NA_real_)
meningitis_code <- v(41L)
require_data(all(na.omit(meningitis_code) %in% c("是", "否", "疑诊", "可疑")),
             "Invalid meningitis code.")
meningitis_binary <- ifelse(is.na(meningitis_code), NA_integer_, as.integer(meningitis_code != "否"))

cat("Module 01 validated ", nrow(inf), " infants and ", nrow(iso), " isolates in ",
    nrow(site_map), " hospitals.\n", sep = "")

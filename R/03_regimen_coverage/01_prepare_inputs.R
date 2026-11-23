# ==============================================================================
# Module 03: Prepare Inputs for Empirical Antibiotic Regimen Coverage Analysis
# ==============================================================================
# Uses unified logical inputs (cohort_data, site_map) and static methodological
# configurations (config/field_mappings.R, config/intrinsic_resistance_rules.R,
# config/calibration_weights.R). Reads the controlled isolate-scope sidecar; no record memberships are embedded.
# ==============================================================================

raw_path  <- normalizePath(cfg$input_files$cohort_data, winslash = "/", mustWork = TRUE)
site_path <- normalizePath(cfg$input_files$site_map,    winslash = "/", mustWork = TRUE)
out       <- output_dir

sha <- function(path) digest::digest(file = path, algo = "sha256")
paths  <- c(Raw = raw_path, SiteMap = site_path)
hashes <- vapply(paths, sha, "")

qc <- list()
record <- function(check, detail, status = "PASS") {
  qc[[length(qc) + 1L]] <<- data.frame(Check = check, Status = status, Detail = detail)
  message("[", status, "] ", check, ": ", detail)
}
check <- function(ok, label, detail) {
  if (!isTRUE(ok)) stop(label, " | ", detail, call. = FALSE)
  record(label, detail)
}
bind <- function(x) as.data.frame(dplyr::bind_rows(x))
near <- function(a, b, tol = 1e-10) isTRUE(all.equal(unname(a), unname(b), tolerance = tol, check.attributes = FALSE))
as_token <- function(x) {
  x <- as.matrix(x); storage.mode(x) <- "character"; x[is.na(x)] <- ""; x
}

raw_tbl <- read_csv_safe(raw_path, encoding = "UTF-8")
stopifnot(identical(dim(raw_tbl), c(1196L, 95L)))
csv <- as_token(raw_tbl)
n <- nrow(csv)

dictionary <- static_cfg$get_raw_field_dictionary()
check(identical(colnames(csv), dictionary$Source_field),
      "Input schema", "All 95 headers and their order match the field dictionary")

continued <- which(csv[, 1] == "")
parents <- continued - 1L
patient_cols <- c(1:43, 77:95)
check(length(continued) == 7L && all(continued > 1L) &&
        all(csv[continued, 1:4] == "") && all(csv[parents, 1] != ""),
      "Continuation roster", "Seven approved second-isolate records retain separate AST")

linked <- csv
linked[continued, patient_cols] <- csv[parents, patient_cols]
check(identical(linked[, -patient_cols], csv[, -patient_cols]) &&
        identical(linked[-continued, ], csv[-continued, ]),
      "Restoration scope", "Exactly 62 patient fields inherited on continuation rows")

# Map hospitals to anonymised Site 1..Site 30 from controlled site_map
sm_raw <- read_csv_safe(site_path, encoding = "auto")
sm <- sm_raw[nzchar(trimws(sm_raw[["Hospital"]])) & grepl("^Site [0-9]+$", trimws(sm_raw[["S4 center"]])), , drop = FALSE]
site_idx <- match(trimws(linked[, 2]), trimws(sm[["Hospital"]]))
check(!anyNA(site_idx) && nrow(sm) == 30L,
      "Site crosswalk", "All 30 hospitals mapped to anonymised Site IDs")

groups     <- static_cfg$get_organism_groups()
ast_fields <- static_cfg$get_ast_fields()
drug_fields <- ast_fields$Antibiotic
rules      <- static_cfg$get_coverage_intrinsic_rules()
wg         <- static_cfg$get_calibration_weights()

regions    <- c("Municipality", "Eastern", "Central", "Western")
region_map <- setNames(regions, c("一层", "二层", "三层", "四层"))
reg_vec    <- unname(region_map[linked[, 4]])
wi         <- match(paste(reg_vec, as.numeric(linked[, 3]), sep = "/"),
                    paste(wg$Region, wg$Year, sep = "/"))

scope_path <- cfg$input_files$coverage_scope
if (is.null(scope_path) || !file.exists(scope_path)) stop("Required controlled coverage_scope input missing.", call. = FALSE)
scope <- read.csv(scope_path, check.names = FALSE, na.strings = "", stringsAsFactors = FALSE, fileEncoding = "UTF-8")
stopifnot(identical(names(scope), c("Source_row", "Site", "Organism_original", "Species_mapped", "Approved_24_drugs")),
          identical(as.integer(scope$Source_row), seq_len(n)), !anyNA(scope$Site), !anyNA(scope$Organism_original),
          identical(scope$Site, trimws(sm[["S4 center"]][site_idx])))
org_mapped <- groups$Organism[match(csv[, 44], groups$Raw_group)]
iso <- data.frame(
  Isolate_ID          = sprintf("ISO-%04d", seq_len(n)),
  Patient_ID          = linked[, 1],
  Source_row          = seq_len(n) + 0,
  Site                = trimws(sm[["S4 center"]][site_idx]),
  Hospital            = linked[, 2],
  Year                = as.numeric(linked[, 3]),
  Region              = reg_vec,
  Raw_pathogen_group  = csv[, 44],
  Organism_original   = scope$Organism_original,
  Organism_mapped     = org_mapped,
  Species_mapped      = scope$Species_mapped,
  Gram                = csv[, 45],
  GA_weeks            = as.numeric(linked[, 8]),
  GA_category         = ifelse(linked[, 6] == "", NA_character_, linked[, 6]),
  Continuation        = seq_len(n) %in% continued,
  S3_weight           = wg$Weight[wi],
  stringsAsFactors    = FALSE
)

raw_ast <- csv[, ast_fields$Position]
colnames(raw_ast) <- drug_fields
w <- iso$S3_weight

site_lookup <- unique(iso[c("Hospital", "Site", "Region")])
site_number <- as.integer(sub("^Site ", "", iso$Site))
patient_rows <- split(seq_len(n), iso$Patient_ID)
check(length(patient_rows) == 1189L && identical(as.integer(table(lengths(patient_rows))), c(1182L, 7L)) &&
        nrow(site_lookup) == 30L && !anyDuplicated(site_lookup$Hospital) && !anyDuplicated(site_lookup$Site) &&
        identical(sort(unique(site_number)), 1:30) && setequal(iso$Year, 2015:2024) &&
        !anyNA(iso$Organism_mapped) && !anyNA(iso$Region) && all(is.finite(w) & w > 0),
      "Cohort and hospital nesting", "1196 isolates, 1189 infants, 30 hospital PSUs and four macro-regions")

gach <- c("超早产儿", "极早产儿", "中晚期早产儿", "足月儿")
ga   <- c("Extremely preterm infants", "Very preterm infants", "Moderate to late preterm infants", "Term infants")
iso$GA <- unname(setNames(ga, gach)[iso$GA_category])
species_label <- iso$Species_mapped
species_label[is.na(species_label)] <- trimws(iso$Organism_original[is.na(species_label)])
iso$Source_species_label <- species_label

study_groups <- c(unique(rules$Organism), setdiff(unique(iso$Organism_mapped), rules$Organism))
fungus       <- iso$Organism_mapped == "Fungi"
unidentified <- iso$Organism_mapped == "Unidentified organism(s)"
stopifnot(length(study_groups) == 16L, length(drug_fields) == 29L,
          all(raw_ast %in% c("", "敏感", "中介", "耐药", "S", "I", "R")),
          sum(fungus) == 20L, sum(unidentified) == 4L, sum(is.na(iso$GA)) == 3L,
          identical(as.integer(table(factor(iso$GA, levels = ga))), c(69L, 200L, 306L, 618L)))

# Load the previously verified applicability memberships from the controlled input.
deferred <- bind(lapply(which(!is.na(scope$Approved_24_drugs)), function(i) {
  drugs <- strsplit(scope$Approved_24_drugs[i], "|", fixed = TRUE)[[1]]
  rule_keys <- paste(rules$Organism, rules$Antibiotic, sep = "/")
  indices <- match(paste(iso$Organism_mapped[i], drugs, sep = "/"), rule_keys)
  stopifnot(!anyNA(indices))
  data.frame(Isolate_ID = iso$Isolate_ID[i], Antibiotic = drugs, Matrix_cell = rules$Matrix_cell[indices])
}))
stopifnot(nrow(deferred) == 24L)

inputs <- data.frame(
  Role   = names(paths),
  File   = basename(paths),
  SHA256 = unname(hashes),
  Bytes  = unname(file.info(paths)$size),
  Source = c("Controlled 11.15.5 cohort CSV", "Controlled hospital-SiteID map"),
  stringsAsFactors = FALSE
)

link_detail <- bind(lapply(seq_along(continued), function(j) {
  i <- continued[j]; a <- csv[i, patient_cols]; b <- linked[i, patient_cols]
  data.frame(Isolate_ID = iso$Isolate_ID[i], Data_row = i, CSV_record = i + 1L,
             Inherited_from_data_row = parents[j], Inherited_from_CSV_record = parents[j] + 1L,
             Parent_isolate_ID = iso$Isolate_ID[parents[j]], Field = colnames(csv)[patient_cols],
             Continuation_raw = a, Parent_value = b, Linked_value = b,
             Outcome = ifelse(a == b, ifelse(a == "", "Both blank", "Confirmed same"),
                              ifelse(a == "", "Restored blank", "Replaced continuation value")),
             row.names = NULL)
}))

link_summary <- bind(lapply(seq_along(continued), function(j) {
  i <- continued[j]; f <- link_detail[link_detail$Data_row == i, ]
  data.frame(Isolate_ID = iso$Isolate_ID[i], Data_row = i, CSV_record = i + 1L,
             Parent_data_row = parents[j], Parent_CSV_record = parents[j] + 1L,
             Parent_isolate_ID = iso$Isolate_ID[parents[j]],
             Patient_ID = linked[i, 1], Hospital = iso$Hospital[i], Site = iso$Site[i],
             Year = iso$Year[i], Region = iso$Region[i], GA_weeks = iso$GA_weeks[i],
             S3_weight = iso$S3_weight[i], Patient_fields_inherited = nrow(f),
             Filled_fields = sum(f$Outcome == "Restored blank"),
             Nonblank_disagreements = sum(f$Outcome == "Replaced continuation value"),
             Preserved_fields = "Pathogen group, Gram, AST, ESBL, duplicate pathogen field, isolate ID")
}))

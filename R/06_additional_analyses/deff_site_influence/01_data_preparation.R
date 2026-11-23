# Sources: accepted 11.16 pathogen inputs_and_design.R; DEFF 11.16 extension.
# Only the original release inputs and the two bundled protocol extracts are read.
trim <- function(x) trimws(as.character(x), whitespace = "[\\h\\v]")
blank <- function(x) is.na(x) | tolower(trim(x)) %in% c("", "na", "n/a", "unknown", "unk", "不详", "未知", "缺失")
canon <- function(x) {
  x <- trim(x)
  x[blank(x)] <- NA_character_
  x
}
near <- function(a, b, tol = 1e-11) isTRUE(all.equal(unname(a), unname(b),
  tolerance = tol, check.attributes = FALSE))
bind <- function(x) do.call(rbind, x)
check <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)

read_input_csv <- function(path, encoding = "UTF-8") {
  x <- readr::read_csv(path, locale = readr::locale(encoding = encoding),
    show_col_types = FALSE, col_types = readr::cols(.default = readr::col_character()),
    na = character(), trim_ws = FALSE, name_repair = "minimal", progress = FALSE)
  if (nrow(readr::problems(x))) stop("Malformed CSV: ", basename(path), call. = FALSE)
  x
}

load_inputs <- function(cfg, static_cfg) {
  cases_path <- normalizePath(cfg$input_files$cohort_data, winslash = "/", mustWork = TRUE)
  sites_path <- normalizePath(cfg$input_files$site_map,    winslash = "/", mustWork = TRUE)
  w_pub <- static_cfg$get_calibration_weights()
  rev_reg <- setNames(c("一层", "二层", "三层", "四层"), c("Municipality", "Eastern", "Central", "Western"))
  w_df <- data.frame(
    "入院时间" = as.character(w_pub$Year),
    "分层"     = unname(rev_reg[w_pub$Region]),
    n_s        = as.character(w_pub$Isolates),
    n_p        = as.character(w_pub$Births),
    wg         = as.character(w_pub$Weight),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  seven_drugs <- c("Penicillin", "Ampicillin", "Ceftazidime", "Cefotaxime", "Gentamicin", "Imipenem", "Meropenem")
  af <- static_cfg$get_ast_fields()
  af7 <- af[match(seven_drugs, af$Drug), c("Drug", "Field"), drop = FALSE]
  all_cov <- static_cfg$get_coverage_intrinsic_rules()
  five_drugs <- c("Penicillin", "Ampicillin", "Ceftazidime", "Cefotaxime", "Gentamicin")
  cr <- all_cov[all_cov$Antibiotic %in% five_drugs, , drop = FALSE]
  cr24 <- data.frame(Matrix_cell = cr$Matrix_cell, Organism = cr$Organism,
    Drug = cr$Antibiotic, Value = cr$Matrix_value, stringsAsFactors = FALSE)
  hashes <- data.frame(
    Role = c("Cases", "Sites"),
    Filename = basename(c(cases_path, sites_path)),
    Observed_SHA256 = vapply(c(cases_path, sites_path), function(p) digest::digest(file = p, algo = "sha256"), ""),
    Status = "PASS", stringsAsFactors = FALSE
  )
  list(
    Cases          = read_input_csv(cases_path, "UTF-8"),
    Sites          = read_input_csv(sites_path, "UTF-8"),
    Weights        = w_df,
    AST_fields     = af7,
    Coverage_rules = cr24,
    hashes         = hashes
  )
}

pick <- function(x, cols) {
  i <- match(cols, names(x))
  if (anyNA(i)) stop("Input fields missing: ", paste(cols[is.na(i)], collapse = ", "), call. = FALSE)
  as.data.frame(setNames(lapply(i, function(j) trim(x[[j]])), cols), optional = TRUE)
}

regions <- c("Municipality", "Eastern", "Central", "Western")
region_map <- setNames(regions, c("一层", "二层", "三层", "四层"))
orgs <- c("大肠埃希菌", "GBS", "CoNS", "克雷伯菌属", "肠球菌属", "产单核细胞李斯特菌", "链球菌属",
  "金黄色葡萄球菌", "肠杆菌属", "不动杆菌属", "假单胞菌属", "沙雷菌属", "其他革兰氏阳性菌",
  "其他革兰氏阴性菌", "真菌", "菌不详")
org_names <- c("E. coli", "GBS", "CoNS", "Klebsiella spp.", "Enterococcus spp.", "L. monocytogenes",
  "Other Streptococcus spp.", "S. aureus", "Enterobacter spp.", "Acinetobacter spp.", "Pseudomonas spp.",
  "Serratia spp.", "Other Gram-positive bacteria", "Other Gram-negative bacteria", "Fungi", "Unidentified organism(s)")

prepare_data <- function(inputs) {
  raw <- inputs$Cases
  taxon <- "菌群整理归类（归类到属）"
  key <- c("ID", "医院名称", "入院时间", "分层")
  linked <- pick(raw, c(key, taxon, "Gram"))
  prim <- which(!blank(linked$ID))
  extra <- which(blank(linked$ID))
  check(identical(dim(raw), c(1196L, 95L)) && length(prim) == 1189L && length(extra) == 7L &&
    !anyDuplicated(linked$ID[prim]) && identical(trim(raw[[44]]), trim(raw[[76]])) &&
    all(vapply(linked[extra, key], function(v) all(blank(v)), TRUE)) &&
    all(extra > 1L) && all((extra - 1L) %in% prim), "Locked cohort or continuation-row structure changed.")
  # The validated source stores the additional isolate immediately after its primary row.
  # Only the infant-level linking fields are inherited; organism and AST stay isolate-specific.
  linked[extra, key] <- linked[extra - 1L, key]
  check(identical(as.integer(table(table(linked$ID))), c(1182L, 7L)), "Invalid isolate-to-infant linkage.")

  sites <- inputs$Sites
  sites <- sites[grepl("^Site [0-9]+$", sites[["S4 center"]]), ]
  index <- match(linked[["医院名称"]], sites[["Hospital"]])
  check(nrow(sites) == 30L && !anyNA(index) && !anyDuplicated(sites[["Hospital"]]) &&
    setequal(sites[["S4 center"]], paste("Site", 1:30)) &&
    all(sites[["S4 stratum"]][index] == region_map[linked[["分层"]]]), "Hospital or macro-region mapping failed.")
  iso <- data.frame(InfantID = linked$ID, Site = sites[["S4 center"]][index],
    Region = factor(sites[["S4 stratum"]][index], levels = regions),
    Year = as.integer(linked[["入院时间"]]),
    Organism = factor(unname(setNames(org_names, orgs)[linked[[taxon]]]), levels = org_names))
  ordered_sites <- function(x) sort(as.integer(sub("^Site ", "", unique(x))))
  check(!anyNA(iso) && setequal(iso$Year, 2015:2024) &&
    identical(ordered_sites(iso$Site[iso$Region == "Central"]), 17:23),
    "Unclassified organism, incorrect year range or unexpected Central codes.")
  roster <- unique(iso[c("Site", "Region")])
  counts <- data.frame(Region = regions, Infants = as.integer(table(iso$Region[prim])),
    Isolates = as.integer(table(iso$Region)), Hospitals = as.integer(table(roster$Region)))
  check(identical(counts$Infants, c(299L,279L,301L,310L)) &&
    identical(counts$Isolates, c(299L,281L,304L,312L)) &&
    identical(counts$Hospitals, c(9L,7L,7L,7L)), "Unexpected macro-region cohort counts.")
  infant <- data.frame(hospital_id = iso$Site[prim], macroregion = iso$Region[prim],
    w_unit = 1, target = as.integer(iso$InfantID[prim] %in% iso$InfantID[iso$Organism == "E. coli"]))
  check(identical(as.integer(rowsum(infant$target, infant$macroregion)), c(56L,62L,86L,71L)),
    "Unexpected E. coli-positive infant counts.")

  # Original S3 weights are ratios of live-birth and sample shares within each year.
  weights <- inputs$Weights
  wk <- function(y, s) paste(y, s, sep = "/")
  year <- as.integer(weights[["入院时间"]])
  key_w <- wk(year, weights[["分层"]])
  iw <- match(wk(iso$Year, linked[["分层"]]), key_w)
  n_cell <- tabulate(iw, nbins = length(key_w))
  pop <- as.numeric(weights[["n_p"]])
  wg <- as.numeric(weights[["wg"]])
  frac <- function(x) x / ave(x, year, FUN = sum)
  check(!anyNA(iw) && length(key_w) == 40L && !anyDuplicated(key_w) &&
    all(is.finite(wg) & wg > 0 & is.finite(pop) & pop > 0) && all(n_cell > 0) &&
    all(n_cell == as.integer(weights[["n_s"]])) &&
    max(abs(wg - frac(pop) / frac(n_cell))) < 1e-8, "Calibration input, support or weight mismatch.")
  iso$wg <- wg[iw]
  weight_cells <- data.frame(Region = unname(region_map[weights[["分层"]]]), Year = year,
    Full_n = n_cell, Target_birth_share = frac(pop), Original_weight = wg)

  list(iso = iso, infant = infant, counts = counts, raw = raw, iw = iw,
    weight_cells = weight_cells, hashes = inputs$hashes,
    ast_fields = inputs$AST_fields, coverage_rules = inputs$Coverage_rules, double_infections = length(extra))
}


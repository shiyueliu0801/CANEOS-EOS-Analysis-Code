# Cohort and calibration adapted from the accepted pathogen inputs_and_design.R
# and clinical cohort.R. Hospital labels in the extract identify distinct PSUs.
input_paths <- c(Cases = normalizePath(cfg$input_files$cohort_data, winslash = "/", mustWork = TRUE))

read_text <- function(path, encoding = "UTF-8") {
  x <- readr::read_csv(path, locale = readr::locale(encoding = encoding),
    col_types = readr::cols(.default = readr::col_character()), show_col_types = FALSE,
    name_repair = "minimal", trim_ws = FALSE, na = character(), progress = FALSE)
  stopifnot(nrow(readr::problems(x)) == 0L)
  x
}
clean <- function(x) trimws(as.character(x), whitespace = "[\\h\\v]")
canon <- function(x) {
  x <- clean(x)
  x[tolower(x) %in% c("", "na", "n/a", "unk", "unknown", "不详", "未知", "缺失")] <- NA_character_
  x
}

source_fields <- c(D72 = "生后72小时内是否死亡", C72 = "生后72小时内是否死亡或放弃",
  D7 = "生后7天内是否死亡", C7 = "生后7天内是否死亡或放弃", Death = "死亡", Composite = "死亡或放弃治疗",
  BPD = "BPD（支气管肺发育不良）", Shock = "休克", LOS = "院内晚发败血症",
  ROP = "ROP（早产儿视网膜病变）", IVH = "IVH（早产儿脑室内出血）", NEC = "NEC（新生儿坏死性小肠结肠炎）")
regimen_fields <- c("First use after Birth", "Use within 72 h", "Use within 72 h–7 days")
modification_fields <- c("72h内是否更改", "72h-7d是否更改")
prophylaxis_field <- "是否预防性使用抗真菌药物"
classic_field <- "是否经典方案"
key <- c("ID", "医院名称", "入院时间", "分层")
taxon <- "菌群整理归类（归类到属）"
raw <- read_text(input_paths[["Cases"]])
clinical_fields <- c(source_fields, regimen_fields, modification_fields, prophylaxis_field, classic_field)
required <- c(key, taxon, "Gram", clinical_fields)
if (!all(required %in% names(raw))) stop("Missing case fields: ", paste(setdiff(required, names(raw)), collapse = ", "))
stopifnot(ncol(raw) == 95L, identical(which(duplicated(names(raw))), 76L),
  identical(names(raw)[c(44L, 76L)], rep(taxon, 2L)), identical(clean(raw[[44L]]), clean(raw[[76L]])))

ids <- canon(raw$ID)
primary <- which(!is.na(ids))
extra <- which(is.na(ids))
stopifnot(!anyDuplicated(ids[primary]), all(extra > 1L), all((extra - 1L) %in% primary))
linked <- lapply(raw[key], canon)
for (field in key) {
  x <- linked[[field]]
  stopifnot(!anyNA(x[primary]), all(is.na(x[extra])))
  x[extra] <- x[extra - 1L]
  linked[[field]] <- x
}
for (field in setdiff(clinical_fields, classic_field)) {
  x <- canon(raw[[field]])
  stopifnot(isTRUE(all(is.na(x[extra]) | x[extra] == x[extra - 1L])))
}

regions <- c("Municipality", "Eastern", "Central", "Western")
region <- factor(match(linked[["分层"]], c("一层","二层","三层","四层")), levels = 1:4, labels = regions)
hospital <- linked[["医院名称"]]
# Match in record order to give portable integer labels for unchanged hospital partitions.
site <- match(hospital, unique(hospital))
year <- as.integer(linked[["入院时间"]])
stopifnot(!anyNA(region), all(year %in% 2015:2024), !anyDuplicated(unique(data.frame(site, region))$site))

taxa <- c("大肠埃希菌" = "E. coli", "GBS" = "GBS", "CoNS" = "CoNS", "克雷伯菌属" = "Klebsiella spp.",
  "肠球菌属" = "Enterococcus spp.", "产单核细胞李斯特菌" = "L. monocytogenes", "链球菌属" = "Other Streptococcus spp.",
  "金黄色葡萄球菌" = "S. aureus", "肠杆菌属" = "Enterobacter spp.", "不动杆菌属" = "Acinetobacter spp.",
  "假单胞菌属" = "Pseudomonas spp.", "沙雷菌属" = "Serratia spp.", "其他革兰氏阳性菌" = "Other Gram-positive bacteria",
  "其他革兰氏阴性菌" = "Other Gram-negative bacteria", "真菌" = "Fungi", "菌不详" = "Unidentified organism(s)")
organism <- unname(taxa[clean(raw[[taxon]])])
gram <- unname(c("G+" = "Gram-positive", "G-" = "Gram-negative", "真菌" = "Fungi", "不详" = "Unidentified")[clean(raw$Gram)])
stopifnot(!anyNA(organism), !anyNA(gram))

calibration <- static_cfg$get_calibration_weights()
weight_key <- paste(calibration$Year, calibration$Region, sep = "/")
weight_index <- match(paste(year, as.character(region), sep = "/"), weight_key)
cell_n <- tabulate(weight_index, nbins = length(weight_key))
cell_weight <- as.numeric(calibration$Weight)
cell_births <- as.numeric(calibration$Births)
year_fraction <- function(x) x / ave(x, calibration$Year, FUN = sum)
stopifnot(length(weight_key) == 40L, !anyDuplicated(weight_key), !anyNA(weight_index),
  all(is.finite(cell_weight) & cell_weight > 0), all(cell_n == as.integer(calibration$Isolates)),
  max(abs(cell_weight - year_fraction(cell_births)/year_fraction(cell_n))) < 1e-8)
iso <- data.frame(Site = site, Region = region, Weight = cell_weight[weight_index],
  Organism = organism, Gram = gram)
isolate_cons <- organism == "CoNS"
infant_cons <- tapply(isolate_cons, linked$ID, any)[ids[primary]]
mixed <- tapply(isolate_cons, linked$ID, function(x) any(x) && !all(x))[ids[primary]]
stopifnot(!anyNA(infant_cons), !anyNA(mixed), !any(mixed))
infants <- data.frame(Site = site[primary], Region = region[primary], Weight = 1,
  Group = ifelse(infant_cons, "CoNS", "Non-CoNS"))
cohort <- data.frame(Unit = c("Isolate", "Infant"), Total = c(nrow(iso), nrow(infants)),
  CoNS = c(sum(isolate_cons), sum(infant_cons)), Non_CoNS = c(sum(!isolate_cons), sum(!infant_cons)),
  Polymicrobial_infants = sum(table(linked$ID) > 1L), Mixed_CoNS_infants = sum(mixed),
  Hospitals = length(unique(site)), Strata = nlevels(region))
stopifnot(all(cohort$Total == cohort$CoNS + cohort$Non_CoNS),
  identical(as.integer(table(table(linked$ID))), c(1182L, 7L)),
  all(!isolate_cons[extra] & !isolate_cons[extra - 1L]),
  identical(cohort$Total, c(1196L, 1189L)), identical(cohort$CoNS, c(193L, 193L)),
  cohort$Hospitals[1L] == 30L, all(tabulate(as.integer(region)) > 0L))
isolate_domains <- list(Overall = rep(TRUE, nrow(iso)), NonCoNS = !isolate_cons)
infant_domains <- list(Overall = rep(TRUE, nrow(infants)), CoNS = infant_cons, NonCoNS = !infant_cons)

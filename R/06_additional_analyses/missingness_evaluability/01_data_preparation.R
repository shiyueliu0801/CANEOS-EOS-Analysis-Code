clean <- function(x) trimws(as.character(x), whitespace = "[\\h\\v]")
blank <- function(x) is.na(x) | clean(x) == ""
is_missing <- function(x) blank(x) | tolower(clean(x)) %in%
  c("na", "n/a", "n.a.", "null", "nan", "unk", "unknown", "不详", "未知", "缺失")
regions <- c("Municipality", "Eastern", "Central", "Western")
outcomes <- c("Survived to discharge", "In-hospital death", "Withdrawal of care")
read_text <- function(path) {
  x <- readr::read_csv(path, col_types = readr::cols(.default = readr::col_character()),
    locale = readr::locale(encoding = "UTF-8"), na = character(), trim_ws = FALSE,
    name_repair = "minimal", show_col_types = FALSE)
  if (nrow(readr::problems(x))) stop("Malformed CSV: ", basename(path))
  as.data.frame(x, check.names = FALSE)
}
prepare_data <- function(path, root) {
  hash <- digest::digest(file = path, algo = "sha256")
  if (hash != "19aedb6d6668a11d4caceb4a9c57b04355c3a13bfb0df02df5e5107ef44cbe48") {
    stop("The case file does not match the locked 11.15.5 SHA-256. No analysis performed.")
  }
  raw <- read_text(path)
  stopifnot(identical(dim(raw), c(1196L, 95L)))
  raw[] <- lapply(raw, clean)
  extra <- which(blank(raw[[1]])); primary <- which(!blank(raw[[1]]))
  stopifnot(length(primary) == 1189L, length(extra) == 7L,
    !anyDuplicated(raw[[1]][primary]), all(extra > 1L), all((extra - 1L) %in% primary),
    all(as.matrix(raw[extra, 1:4]) == ""), identical(raw[[44]], raw[[76]]))
  linked <- raw
  # The accepted adjacent continuation rule restores infant metadata only.
  patient_columns <- c(1:43, 77:95)
  linked[extra, patient_columns] <- raw[extra - 1L, patient_columns]
  stopifnot(identical(as.integer(table(table(linked[[1]]))), c(1182L, 7L)))
  region <- unname(setNames(regions, c("一层", "二层", "三层", "四层"))[linked[[4]]])
  year <- suppressWarnings(as.integer(linked[[3]]))
  stopifnot(!anyNA(region), !anyNA(year), setequal(year, 2015:2024), !any(blank(linked[[2]])))
  sites <- unique(data.frame(hospital = linked[[2]], region = region))
  stopifnot(nrow(sites) == 30L, !anyDuplicated(sites$hospital),
    identical(as.integer(table(factor(sites$region, levels = regions))), c(9L, 7L, 7L, 7L)))
  pair <- paste(linked[[87]], linked[[88]], sep = "/")
  group_map <- setNames(c(outcomes, NA_character_), c("否/否", "是/是", "否/是", "不详/不详"))
  stopifnot(all(pair %in% names(group_map)))
  context <- data.frame(hospital = factor(linked[[2]]), year = year,
    region = factor(region, levels = regions), unit_weight = 1,
    outcome = factor(unname(group_map[pair]), levels = outcomes),
    year_factor = factor(year, levels = 2015:2024))
  infants <- context[primary, ]
  stopifnot(identical(as.integer(table(infants$region)), c(299L, 279L, 301L, 310L)),
    identical(as.integer(table(context$region)), c(299L, 281L, 304L, 312L)),
    identical(as.integer(table(infants$outcome, useNA = "always")), c(954L, 98L, 128L, 9L)),
    identical(as.integer(table(context$outcome, useNA = "always")), c(960L, 99L, 128L, 9L)))
  dictionary <- static_cfg$get_s19_missingness_dictionary()
  dictionary <- dictionary[dictionary$Module == "Clinical", ]
  dictionary$Column <- as.integer(dictionary$Column)
  dictionary$OutcomeCircular <- dictionary$OutcomeCircular == "TRUE"
  stopifnot(nrow(dictionary) == 56L, !anyDuplicated(dictionary$Key),
    identical(names(raw)[dictionary$Column], dictionary$SourceField))
  list(raw = raw, linked = linked, primary = primary, extra = extra, context = context,
    infants = infants, dictionary = dictionary, sites = sites, hash = hash, path = path)
}

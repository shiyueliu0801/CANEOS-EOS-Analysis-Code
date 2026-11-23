# Shared analysis functions retained from the accepted revision.
require_data <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}
year_set <- 2015:2024
region_names <- c("Municipality", "Eastern", "Central", "Western")
regions <- c("Overall", region_names)
region_map <- setNames(region_names, c("一层", "二层", "三层", "四层"))
clean <- function(v) trimws(as.character(v), whitespace = "[\\h\\v]")
missing_tokens <- c("", "na", "n/a", "unk", "unknown", "不详", "未知", "缺失")
is_missing <- function(v) is.na(v) | tolower(clean(v)) %in% missing_tokens
canon <- function(v) {
  v <- clean(v)
  v[is_missing(v)] <- NA_character_
  v
}
read_source <- function(path, encoding = "UTF-8") {
  tab <- readr::read_csv(path, col_types = readr::cols(.default = readr::col_character()),
    locale = readr::locale(encoding = encoding), name_repair = "minimal",
    na = character(), trim_ws = FALSE, show_col_types = FALSE)
  if (nrow(readr::problems(tab))) stop("CSV 解析错误：", path)
  tab
}
extract <- function(tab, fields) {
  j <- match(fields, names(tab))
  if (anyNA(j)) stop("字段不存在：", paste(fields[is.na(j)], collapse = "、"))
  as.data.frame(setNames(lapply(j, function(i) clean(tab[[i]])), fields), optional = TRUE)
}
capture <- function(expr) {
  messages <- character()
  value <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    messages <<- c(messages, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) { messages <<- c(messages, conditionMessage(e)); NULL })
  list(Fit = value, Message = paste(unique(messages), collapse = " | "))
}
combine_text <- function(...) paste(Filter(nzchar, list(...)), collapse = " | ")
design <- function(d) survey::svydesign(ids = ~HospitalID, strata = ~Region, weights = ~w, nest = TRUE, data = d)
psu_count <- function(d) length(unique(d$HospitalID))
design_df <- function(d) psu_count(d) - length(unique(d$Region))
region_rows <- function(d, g) if (g == "Overall") rep(TRUE, nrow(d)) else d$Region == g
fmt_p <- function(p) ifelse(is.na(p), "—", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
fmt_dec <- function(x, n = 2) ifelse(is.na(x), "—", sprintf(paste0("%.", n, "f"), ifelse(abs(x) < 0.5 * 10^-n, 0, x)))
fmt_interval <- function(x, l, u, n = 2) ifelse(is.na(x), "—",
  ifelse(is.na(l) | is.na(u), paste0(fmt_dec(x, n), " (—)"),
         paste0(fmt_dec(x, n), " (", fmt_dec(l, n), "–", fmt_dec(u, n), ")")))
count_text <- function(x) format(x, trim = TRUE, big.mark = ",", scientific = FALSE)
write_result <- function(x, stem) {
  readr::write_csv(x, file.path(output_dir, paste0(stem, ".csv")), na = "")
}
pub_row <- function(label, values = rep("", 5L), p = "") {
  as.data.frame(setNames(as.list(c(label, unname(values), p)), c("Characteristic", regions, "P value")),
                check.names = FALSE, stringsAsFactors = FALSE)
}
trend_skeleton <- function(d, method, effect, predictors = 1L) {
  data.frame(Method = method, Effect = effect, N = nrow(d), Hospitals = psu_count(d),
    DesignDF = design_df(d), DF = design_df(d) - predictors,
    FirstYear = if (nrow(d)) min(d$Year) else NA_real_, LastYear = if (nrow(d)) max(d$Year) else NA_real_,
    Years = length(unique(d$Year)), Events = NA_real_, NonEvents = NA_real_,
    EventHospitals = NA_real_, NonEventHospitals = NA_real_,
    Beta = NA_real_, SE = NA_real_, P = NA_real_, Annual = NA_real_, AnnualLower = NA_real_, AnnualUpper = NA_real_,
    Five = NA_real_, FiveLower = NA_real_, FiveUpper = NA_real_,
    POTestP = NA_real_, Converged = NA, Status = "Not estimable",
    Diagnostic = "")
}
fill_trend <- function(a, b, se) {
  if (length(b) != 1L || length(se) != 1L || !is.finite(b) || !is.finite(se) || se <= 0 || a$DF <= 0) {
    a$Diagnostic <- combine_text(a$Diagnostic, "Non-finite estimate, variance or denominator degrees of freedom.")
    return(a)
  }
  a$Beta <- b
  a$SE <- se
  a$P <- 2 * pt(-abs(b / se), a$DF)
  delta <- qt(0.975, a$DF) * se
  z <- c(b, b - delta, b + delta, 5 * b, 5 * (b - delta), 5 * (b + delta))
  if (a$Effect != "Mean change") z <- exp(z)
  if (any(!is.finite(z))) stop("Non-finite transformed result.")
  a[c("Annual", "AnnualLower", "AnnualUpper", "Five", "FiveLower", "FiveUpper")] <- as.list(z)
  a$Status <- "Estimated"
  a
}

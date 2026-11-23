# All display values are derived from unrounded R results.
format_number <- function(x, digits = 1L) {
  if (is.na(x)) return("—")
  if (abs(x) < .5*10^-digits) x <- 0
  sub("-", "−", sprintf(paste0("%.", digits, "f"), x), fixed = TRUE)
}
format_count <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
format_nN <- function(n, N) paste0(format_count(n), "/", format_count(N))
format_p <- function(p) if (is.na(p)) "—" else if (p < .001) "<0.001" else format_number(p, 3L)
format_change <- function(x) {
  s <- format_number(x)
  if (!is.na(x) && x > 0 && s != "0.0") paste0("+", s) else s
}
format_interval <- function(x, lo, hi, scale = 1) {
  p <- format_number(scale*x)
  if (is.na(lo) || is.na(hi)) return(p)
  paste0(p, " (", format_number(scale*lo), " to ", format_number(scale*hi), ")")
}

sections <- c("Pathogen distribution", "Empirical antibiotic regimen coverage", "Paired coverage differences",
  "Postnatal antimicrobial treatment and modification", "Clinical outcomes and complications")
display_row <- function(r, section, on_word = TRUE) {
  pair <- section == sections[3L]
  clinical <- section %in% sections[4:5]
  group <- function(prefix) {
    f <- function(key) r[[paste0(prefix, "_", key)]]
    n <- if (pair) format_count(f("N")) else format_nN(f("Events"), f("N"))
    x <- if (pair) format_interval(f("Difference_pp"), f("Lower"), f("Upper")) else {
      show_ci <- !clinical || r$Code %in% disposition
      if (show_ci) format_interval(f("Estimate"), f("Lower"), f("Upper"), 100) else format_number(100*f("Estimate"))
    }
    c(n = n, result = x)
  }
  all <- group("Overall")
  non <- group("NonCoNS")
  cons <- if (clinical) group("CoNS") else c(n = "", result = "")
  if (section == sections[1L] && isTRUE(r$Excluded)) non <- c(n = "Excluded", result = "—")
  data.frame(Section = section, Code = r$Code, Indicator = r$Indicator, Displayed_in_Word = on_word,
    Overall_n = all[1L], Overall_statistic = all[2L], CoNS_n = cons[1L], CoNS_statistic = cons[2L],
    NonCoNS_n = non[1L], NonCoNS_statistic = non[2L], P = if (clinical) format_p(r$P) else "",
    q = if (clinical) format_p(r$q) else "", Overall_pair_P = if (pair) format_p(r$Overall_P) else "",
    Overall_pair_q = if (pair) format_p(r$Overall_q) else "", NonCoNS_pair_P = if (pair) format_p(r$NonCoNS_P) else "",
    NonCoNS_pair_q = if (pair) format_p(r$NonCoNS_q) else "", Change_pp = format_change(r$Change_pp), row.names = NULL)
}

numeric_tables <- list(Pathogen = pathogens, Coverage = coverage, Paired = paired, Treatment = treatments, Outcomes = outcomes)
display <- bind(lapply(seq_along(numeric_tables), function(i) {
  t <- numeric_tables[[i]]
  bind(lapply(seq_len(nrow(t)), function(j) display_row(t[j,], sections[i], t$Code[j] != "Prophylaxis")))
}))
display$Word_data_row <- NA_integer_
display$Word_data_row[display$Displayed_in_Word] <- seq_len(sum(display$Displayed_in_Word))
stopifnot(nrow(display) == 67L, sum(display$Displayed_in_Word) == 66L,
  identical(as.integer(table(factor(display$Section, levels = sections))), c(8L,12L,4L,28L,15L)),
  all(treatment_tests$Planned_tests == 40L & treatment_tests$Estimated_tests == 40L),
  identical(as.integer(table(outcomes$Family)), c(13L,2L)),
  all(paired_bounds$Family95_low <= paired_bounds$Joint95_low),
  all(paired_bounds$Joint95_high <= paired_bounds$Family95_high))

dir.create(output, recursive = TRUE, showWarnings = FALSE)
result <- c(numeric_tables, list(
  Paired_Bounds       = paired_bounds,
  Treatment_BH_Family = treatment_tests,
  Cohort              = cohort,
  S18_Table_Display   = display
))
for (name in names(result)) readr::write_csv(result[[name]], file.path(output, paste0(name, ".csv")), na = "")
message("Module 05 (S18 Table CoNS sensitivity) completed: ", nrow(display), " indicators across 5 sections.")

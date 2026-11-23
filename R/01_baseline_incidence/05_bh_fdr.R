# ==============================================================================
# Module 01: BH-FDR Adjustment for Table 1 (F1) and S11 Table (F3-F6)
# ==============================================================================
# Consumes full-precision unadjusted P values directly produced by 02_incidence,
# 03_baseline, and 04_trends. Retains the locked planned family sizes:
#   F1 (Table 1) : All macro-region omnibus tests including incidence (M = 23)
#   F3 (S11)     : Overall EOS incidence trend (M = 1)
#   F4 (S11)     : Four regional EOS incidence trends (M = 4)
#   F5 (S11)     : Overall baseline & early clinical trends (M = 34)
#   F6 (S11)     : Four regions jointly: baseline & early clinical trends (M = 136)
# ==============================================================================

inc_reg <- readr::read_csv(file.path(output_dir, "Incidence_Regional_Comparison.csv"), show_col_types = FALSE)
base_reg <- readr::read_csv(file.path(output_dir, "Baseline_Regional_Tests.csv"), show_col_types = FALSE)
inc_trd <- readr::read_csv(file.path(output_dir, "Incidence_Trends.csv"), show_col_types = FALSE)
base_trd <- readr::read_csv(file.path(output_dir, "Baseline_Trend_Models.csv"), show_col_types = FALSE)

# --- Family F1: Table 1 (M = 23) ---
t1_base <- base_reg[base_reg$Table == "Table 1", , drop = FALSE]
t1_pq <- data.frame(
  Table    = "Table 1",
  FamilyID = "F1",
  Variable = c("Hospital-based EOS incidence", t1_base$Variable),
  Subgroup = "Regional comparison",
  FullP    = c(as.numeric(inc_reg$P[1L]), as.numeric(t1_base$P)),
  PlannedM = 23L,
  stringsAsFactors = FALSE
)
require_data(nrow(t1_pq) == 23L, paste("Expected 23 planned tests in F1 (Table 1), got", nrow(t1_pq)))
t1_pq$BH_q     <- bh_adjust_planned(t1_pq$FullP, planned_m = 23L)
t1_pq$DisplayP <- format_p_display(t1_pq$FullP)
t1_pq$DisplayQ <- format_p_display(t1_pq$BH_q)
write_csv_utf8(t1_pq, file.path(output_dir, "Table1_BH_FDR.csv"), full_precision = TRUE)

# --- Families F3, F4, F5, F6: S11 Table (Temporal Trends) ---
s11_inc_all <- inc_trd[inc_trd$Region == "Overall", , drop = FALSE]
s11_inc_reg <- inc_trd[inc_trd$Region != "Overall", , drop = FALSE]
s11_bt_all  <- base_trd[base_trd$Region == "Overall", , drop = FALSE]
s11_bt_reg  <- base_trd[base_trd$Region != "Overall", , drop = FALSE]

require_data(nrow(s11_inc_all) == 1L && nrow(s11_inc_reg) == 4L &&
               nrow(s11_bt_all) == 34L && nrow(s11_bt_reg) == 136L,
             "Unexpected row counts for S11 Table BH-FDR families F3-F6.")

make_fam_df <- function(fam_id, vars, subs, pvals, m) {
  q <- bh_adjust_planned(pvals, planned_m = m)
  data.frame(
    Table    = "S11 Table",
    FamilyID = fam_id,
    Variable = vars,
    Subgroup = subs,
    FullP    = as.numeric(pvals),
    PlannedM = as.integer(m),
    BH_q     = q,
    DisplayP = format_p_display(pvals),
    DisplayQ = format_p_display(q),
    stringsAsFactors = FALSE
  )
}

s11_pq <- rbind(
  make_fam_df("F3", "Hospital-based EOS incidence", s11_inc_all$Region, s11_inc_all$P, 1L),
  make_fam_df("F4", rep("Hospital-based EOS incidence", 4L), s11_inc_reg$Region, s11_inc_reg$P, 4L),
  make_fam_df("F5", s11_bt_all$Variable, s11_bt_all$Region, s11_bt_all$ReportedP, 34L),
  make_fam_df("F6", s11_bt_reg$Variable, s11_bt_reg$Region, s11_bt_reg$ReportedP, 136L)
)
write_csv_utf8(s11_pq, file.path(output_dir, "S11_Table_BH_FDR.csv"), full_precision = TRUE)

cat("Module 01 BH-FDR completed: F1 (23 tests), F3 (1), F4 (4), F5 (34), F6 (136).\n")

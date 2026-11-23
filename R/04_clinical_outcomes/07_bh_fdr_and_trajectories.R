# ==============================================================================
# Module 04: BH-FDR Adjustment (F2, F7-F13) & Dynamic Figure 3 Trajectory Build
# ==============================================================================
# Computes exact BH-FDR q-values retaining locked planned family sizes:
#   F2  (Table 4)  : All macro-region omnibus tests (M = 24)
#   F7  (S15)      : Main overall antimicrobial trends (M = 34)
#   F8  (S15)      : Main four-region antimicrobial trends (M = 136)
#   F9  (S15)      : Supplementary overall antimicrobial trends (M = 3)
#   F10 (S15)      : Supplementary four-region antimicrobial trends (M = 12)
#   F11 (S16)      : Overall clinical outcome trends (M = 13)
#   F12 (S16)      : Four regions jointly: clinical outcome trends (M = 52)
#   F13 (S16)      : Four gestational-age groups jointly: outcome trends (M = 52)
#
# Also builds Figure3_Three_Phase_Trajectories.csv directly from infant-level
# treatment records so R/07_figures/figure3 requires zero static real-data CSVs.
# ==============================================================================

# --- Family F2: Table 4 (M = 24) ---
c4 <- clinical_audit$Table4_Tests
stopifnot(nrow(c4) == 24L)
t4_pq <- data.frame(
  Table      = "Table 4",
  FamilyID   = "F2",
  VariableID = c4$Code,
  Variable   = c4$Code,
  Subgroup   = "Regional comparison",
  FullP      = as.numeric(c4$P),
  PlannedM   = 24L,
  stringsAsFactors = FALSE
)
t4_pq$BH_q     <- bh_adjust_planned(t4_pq$FullP, planned_m = 24L)
t4_pq$DisplayP <- format_p_display(t4_pq$FullP)
t4_pq$DisplayQ <- format_p_display(t4_pq$BH_q)
write_csv_utf8(t4_pq, file.path(output_path, "Table4_BH_FDR.csv"), full_precision = TRUE)

# --- Families F7-F10: S15 Table (Antimicrobial Trends; audit name S16_Model_Audit) ---
at <- s16_results$S16_Model_Audit
antifamily <- function(type, group) {
  ifelse(type == "Main",
         ifelse(group == "Overall", "F7", "F8"),
         ifelse(group == "Overall", "F9", "F10"))
}
at$FamilyID <- antifamily(at$FamilyType, at$Subgroup)
at$FullP    <- ifelse(at$Reporting == "Report" & is.finite(at$P), as.numeric(at$P), NA_real_)
at$RecalcQ  <- NA_real_
s15_m <- c(F7 = 34L, F8 = 136L, F9 = 3L, F10 = 12L)
for (fid in names(s15_m)) {
  idx_f <- which(at$FamilyID == fid)
  stopifnot(length(idx_f) == s15_m[[fid]])
  at$RecalcQ[idx_f] <- bh_adjust_planned(at$FullP[idx_f], planned_m = s15_m[[fid]])
  stopifnot(isTRUE(all.equal(at$RecalcQ[idx_f], at$BH_q[idx_f], tolerance = 1e-14)))
}
s15_pq <- data.frame(
  Table      = "S15 Table",
  FamilyID   = at$FamilyID,
  Stage      = at$Stage,
  VariableID = at$Code,
  Variable   = at$Indicator,
  Subgroup   = at$Subgroup,
  Reporting  = at$Reporting,
  FullP      = at$FullP,
  PlannedM   = at$FamilyN,
  BH_q       = at$RecalcQ,
  DisplayP   = format_p_display(at$FullP),
  DisplayQ   = format_p_display(at$RecalcQ),
  stringsAsFactors = FALSE
)
write_csv_utf8(s15_pq, file.path(output_path, "S15_Table_BH_FDR.csv"), full_precision = TRUE)

# --- Families F11-F13: S16 Table (Clinical Outcome Trends; audit name S17_Model_Audit) ---
ct_raw <- clinical_trend_audit$S17_Model_Audit
ct <- ct_raw[ct_raw$Role == "Primary", , drop = FALSE]
reg_list <- c("Municipality", "Eastern", "Central", "Western")
ct$FamilyID <- ifelse(ct$Subgroup == "Overall", "F11",
                      ifelse(ct$Subgroup %in% reg_list, "F12", "F13"))
ct$FullP <- ifelse(ct$Reporting == "Report" & is.finite(ct$P), as.numeric(ct$P), NA_real_)
ct$BH_q  <- NA_real_
s16_m <- c(F11 = 13L, F12 = 52L, F13 = 52L)
for (fid in names(s16_m)) {
  idx_f <- which(ct$FamilyID == fid)
  stopifnot(length(idx_f) == s16_m[[fid]])
  ct$BH_q[idx_f] <- bh_adjust_planned(ct$FullP[idx_f], planned_m = s16_m[[fid]])
}
s16_pq <- data.frame(
  Table      = "S16 Table",
  FamilyID   = ct$FamilyID,
  VariableID = ct$Code,
  Variable   = ct$Outcome,
  Subgroup   = ct$Subgroup,
  Reporting  = ct$Reporting,
  FullP      = ct$FullP,
  PlannedM   = unname(s16_m[ct$FamilyID]),
  BH_q       = ct$BH_q,
  DisplayP   = format_p_display(ct$FullP),
  DisplayQ   = format_p_display(ct$BH_q),
  stringsAsFactors = FALSE
)
write_csv_utf8(s16_pq, file.path(output_path, "S16_Table_BH_FDR.csv"), full_precision = TRUE)

# --- Dynamic Figure 3 Three-Phase Trajectory Aggregation ---
collapse_fig3 <- function(x) {
  x <- as.character(x)
  x[x %in% c("1GC", "2GC")] <- "1GC/2GC"
  x[x %in% c("Other", "UNK")] <- "Other/UNK"
  x
}
t1 <- collapse_fig3(treatment[[regimen_fields[1L]]])
# Sparse T1-only glycopeptide/linezolid regimens (<2% of cohort) are grouped into Other/UNK at T1
t1[t1 %in% c("CARB+GLY/LZD", "GLY/LZD")] <- "Other/UNK"
t2 <- collapse_fig3(treatment[[regimen_fields[2L]]])
t3 <- ifelse(not_assessable, "Death or treatment withdrawal \u226472 h",
             collapse_fig3(treatment$LateRegimen))

traj_tbl <- as.data.frame(table(T1 = t1, T2 = t2, T3 = t3), stringsAsFactors = FALSE)
traj_tbl <- traj_tbl[traj_tbl$Freq > 0L, , drop = FALSE]
names(traj_tbl)[4L] <- "n"
traj_tbl <- traj_tbl[order(traj_tbl$T1, traj_tbl$T2, traj_tbl$T3), , drop = FALSE]
rownames(traj_tbl) <- NULL
stopifnot(sum(traj_tbl$n) == nrow(inf))
write_csv_utf8(traj_tbl, file.path(output_path, "Figure3_Three_Phase_Trajectories.csv"))

cat("Module 04 BH-FDR completed: F2 (24), F7-F10 (185), F11-F13 (117).\n")
cat("Module 04 built Figure 3 trajectory table:", nrow(traj_tbl), "unique trajectories,", sum(traj_tbl$n), "infants.\n")

# Formal results, code, and actual outputs

All output paths are relative to the configured output directory. Outputs contain restricted analyses and must remain outside the public archive. Script paths are relative to this repository.

| Item | Actual generating script(s) | Actual required output(s) | Notes |
|---|---|---|---|
| Table 1 | R/01_baseline_incidence/02_incidence.R; R/01_baseline_incidence/03_baseline.R; R/01_baseline_incidence/05_bh_fdr.R | 01_baseline_incidence/Table1_Baseline.csv; 01_baseline_incidence/Incidence_Regional_Comparison.csv; 01_baseline_incidence/Table1_BH_FDR.csv |  |
| Table 2 | R/02_pathogens_ast/02_pathogens.R; R/02_pathogens_ast/04_write_results.R | 02_pathogens_ast/Table2_Pathogen.csv; 02_pathogens_ast/Table2_Display.csv; 02_pathogens_ast/Table2_Display_P.csv |  |
| Table 3 | R/02_pathogens_ast/03_ast_and_esbl.R; R/02_pathogens_ast/04_write_results.R | 02_pathogens_ast/Table3_S6_AST.csv; 02_pathogens_ast/Table3_Display.csv |  |
| Table 4 | R/04_clinical_outcomes/03_clinical_comparisons.R; R/04_clinical_outcomes/07_bh_fdr_and_trajectories.R | 04_clinical_outcomes/Table4.csv; 04_clinical_outcomes/Table4_BH_FDR.csv |  |
| S1 Appendix | R/00_site_selection/01_municipality_priority_order.R – 04_western_priority_order.R; R/00_site_selection/run.R | 00_site_selection/S1_Appendix_Site_Priority_Order.csv | Standalone anonymised random priority ordering for all 94 candidate SiteCodes (`set.seed(20230701)`). |
| S1 Table | R/01_baseline_incidence/06_supporting_tables.R | 01_baseline_incidence/S1_Table_Site_Year_Isolates.csv | Dynamic site-year isolate tabulation. |
| S2 Table | R/01_baseline_incidence/06_supporting_tables.R; config/calibration_weights.R | 01_baseline_incidence/S2_Table_Calibration.csv | Authorised 40-cell aggregate configuration, full precision. |
| S3 Table | R/01_baseline_incidence/02_incidence.R; R/01_baseline_incidence/run.R | 01_baseline_incidence/S3_Table_Hospital_Incidence.csv; 01_baseline_incidence/S3_Table_Hospital_Incidence_Estimates.csv |  |
| S4 Table | R/02_pathogens_ast/03_ast_and_esbl.R; R/02_pathogens_ast/04_write_results.R | 02_pathogens_ast/S4_Table_ESBL_Display.csv; 02_pathogens_ast/S5_ESBL.csv |  |
| S5 Table | R/02_pathogens_ast/03_ast_and_esbl.R; R/02_pathogens_ast/04_write_results.R | 02_pathogens_ast/S5_Table_AST_Full_Display.csv; 02_pathogens_ast/Table3_S6_AST.csv |  |
| S6 Table | R/04_clinical_outcomes/04_antimicrobial_use.R | 04_clinical_outcomes/S6_Table_Antimicrobial_Use.csv |  |
| S7 Table | R/02_pathogens_ast/02_pathogens.R; R/02_pathogens_ast/04_write_results.R | 02_pathogens_ast/S7_Table_PanelA_Display.csv; 02_pathogens_ast/S7_Table_PanelB_Display.csv; 02_pathogens_ast/S8_Tests.csv |  |
| S8 Table | R/04_clinical_outcomes/03_clinical_comparisons.R | 04_clinical_outcomes/S8_Table_Risk_Ratios.csv |  |
| S9 Table | R/03_regimen_coverage/02_coverage_analysis.R | 03_regimen_coverage/S9_Table_Empirical_Coverage.csv |  |
| S10 Table | R/03_regimen_coverage/02_coverage_analysis.R | 03_regimen_coverage/S10_Table_Paired_Coverage.csv; 03_regimen_coverage/S10_Table_Paired_Logical_Bounds.csv; 03_regimen_coverage/CANEOS_Coverage_11.15.5_Subgroup_Paired.xlsx | Four primary CSV contrasts; workbook sheet Exploratory_32 contains 32 subgroup contrasts. |
| S11 Table | R/01_baseline_incidence/02_incidence.R; R/01_baseline_incidence/04_trends.R; R/01_baseline_incidence/05_bh_fdr.R | 01_baseline_incidence/S11_Table_Baseline_Trends.csv; 01_baseline_incidence/Incidence_Trends.csv; 01_baseline_incidence/S11_Table_BH_FDR.csv |  |
| S12 Table | R/02_pathogens_ast/02_pathogens.R; R/02_pathogens_ast/04_write_results.R | 02_pathogens_ast/S12_Table_Pathogen_Period_Display.csv; 02_pathogens_ast/Pathogen_Period.csv |  |
| S13 Table | R/02_pathogens_ast/03_ast_and_esbl.R; R/02_pathogens_ast/04_write_results.R | 02_pathogens_ast/S13_Table_AMR_ESBL_Trends_Display.csv; 02_pathogens_ast/AMR_Trend.csv; 02_pathogens_ast/ESBL_Trend.csv |  |
| S14 Table | R/03_regimen_coverage/02_coverage_analysis.R | 03_regimen_coverage/S14_Table_Coverage_Trends.csv |  |
| S15 Table | R/04_clinical_outcomes/06_antimicrobial_trends.R; R/04_clinical_outcomes/07_bh_fdr_and_trajectories.R | 04_clinical_outcomes/S15_Table_Antimicrobial_Trends.csv; 04_clinical_outcomes/S15_Table_BH_FDR.csv |  |
| S16 Table | R/04_clinical_outcomes/05_clinical_trends.R; R/04_clinical_outcomes/07_bh_fdr_and_trajectories.R | 04_clinical_outcomes/S16_Table_Clinical_Trends.csv; 04_clinical_outcomes/S16_Table_BH_FDR.csv |  |
| S17 Table | R/06_additional_analyses/regimen_ast_outcomes/02_analysis.R | 06_additional_analyses/s17_regimen_ast_outcomes/S17_submission_values.csv; 06_additional_analyses/s17_regimen_ast_outcomes/S17_overall_tests.csv |  |
| S18 Table | R/05_cons_sensitivity/05_export_results.R | 05_cons_sensitivity/S18_Table_Display.csv; 05_cons_sensitivity/Pathogen.csv; 05_cons_sensitivity/Coverage.csv; 05_cons_sensitivity/Paired.csv; 05_cons_sensitivity/Treatment_BH_Family.csv; 05_cons_sensitivity/Outcomes.csv |  |
| S19 Table | R/06_additional_analyses/missingness_evaluability/00_shared_survey_functions.R; R/06_additional_analyses/missingness_evaluability/04_combine_final_results.R | 06_additional_analyses/s19_missingness_evaluability/S19_Table_Missingness_Evaluability.csv; 06_additional_analyses/s19_missingness_evaluability/statistical_test_details.csv |  |
| S20 Table | R/06_additional_analyses/deff_site_influence/02_s20_deff.R | 06_additional_analyses/s20_s21_deff_site_influence/S20_Table_DEFF_Final.csv |  |
| S21 Table | R/06_additional_analyses/deff_site_influence/03_s21_site23_sensitivity.R | 06_additional_analyses/s20_s21_deff_site_influence/S21_Table_Site23_Sensitivity_Final.csv |  |
| Figure 1 |  |  | Manual recruitment/cohort flowchart; final analytic cohort counts checked by R. Screening/exclusion source records are outside the locked analytic inputs. |
| Figure 2 | R/07_figures/figure2/plot_figure2.R | 07_figures/figure2/Figure2_300dpi.tiff | Data: 03_regimen_coverage/Figure2_Plot_Data.csv; PDF also emitted; native PPTX optional (export_pptx=TRUE). |
| Figure 3 | R/07_figures/figure3/plot_figure3.R | 07_figures/figure3/Figure3_Sankey_Final_v3_300dpi.tiff | Data: 04_clinical_outcomes/Figure3_Three_Phase_Trajectories.csv; TIFF required, no PPTX required. |
| S1 Fig | R/07_figures/s1_figure/prepare_geography.R; R/07_figures/s1_figure/plot_s1_figure.R | 07_figures/s1_figure/S1_Fig_300dpi.tiff; 07_figures/s1_figure/S1_Provincial_Counts_Dynamic.csv | Explicit external geography and restricted hospital-province map required; default full run fails if missing. |
| S2 Fig | R/03_regimen_coverage/03_export_s2_excel.R | 03_regimen_coverage/S2_Heatmap.xlsx; 03_regimen_coverage/S2_Heatmap_Rates.csv; 03_regimen_coverage/S2_Heatmap_Evaluable_N.csv | Pure-R formatted 23×23 coverage heatmap Excel workbook (4 sheets; upper/diagonal unit weights, lower S3 weights; evaluable n). |
| S3 Fig | R/03_regimen_coverage/04_plot_s3_panels.R | 03_regimen_coverage/S3_Fig_Panels/S3_panel_1.tiff; 03_regimen_coverage/S3_Fig_Panels/S3_panel_2.tiff; 03_regimen_coverage/S3_Fig_Panels/S3_panel_3.tiff; 03_regimen_coverage/S3_Fig_Panels/S3_panel_4.tiff | Four 600 dpi TIFF panels. |

## Audit-only optional exports

Module03 defaults to the formal subgroup workbook, formatted S2 workbook and required CSV/TIFF outputs. Other coverage workbooks are disabled unless local config sets `export_internal_audit = TRUE`; they can contain isolate records and must remain internal. The historical all-book JSON export is removed. S2 Fig remains a workbook-based workflow: this package computes its full matrix/colours and workbook, with final figure export performed manually; it does not claim an R-generated S2 TIFF.

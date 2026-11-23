# Statistical Methods, Methodological Rules, and Study-Wide BH-FDR Registry

## 1. Study Design and Complex Survey Inference

- **Study design and site selection**: Multicenter, geographically stratified retrospective cohort of 30 participating tertiary hospitals across 4 geographic macro-region strata in mainland China (Municipality [9 hospitals], Eastern [7], Central [7], Western [7]; 2015–2024). Analyses follow the study protocol and revision specifications documented in S1 Text; no separate prospective standalone statistical analysis plan (SAP) was pre-registered prior to analysis, and detailed inferential models were finalized during analysis and peer review (consistent with study revision responses). No post-protocol model is described here as prospectively prespecified. Final hospital inclusion probabilities are unknown; recruitment is non-probability, with consecutive eligible cases within participating hospitals.
- **Cohort size**: 1,189 unique infants contributing 1,196 culture-confirmed early-onset sepsis (EOS) isolates (1,182 monomicrobial infants + 7 polymicrobial infants with 2 distinct isolates each).
- **Continuation-row rule**: Seven second-isolate records appear immediately below their primary infant row with blank linkage keys (`ID`, `医院名称`, `入院时间`, `分层`); scripts deterministically inherit the 62 infant-level fields from the preceding row while preserving separate isolate-level pathogen, Gram, AST, and ESBL fields.
- **Survey variance estimation**: All design-based estimates use `survey::svydesign(ids = ~Site, strata = ~Region, nest = TRUE)` with `options(survey.lonely.psu = "fail", survey.adjust.domain.lonely = FALSE)`, Taylor-series linearization, and finite design degrees of freedom (`df = contributing PSUs - represented strata`; `df = 26` for the full 30-hospital cohort).
- **Weighting**:
  - **Infant-level analyses** (Table 1, Table 4, S3, S6, S8, S11, S15, S16, S17, S20): Unweighted (`w = 1`) hospital-clustered inference.
  - **Isolate-level pathogen distribution, antimicrobial resistance, ESBL positivity, and estimated in vitro empirical coverage analyses** (Table 2, Table 3, S4, S5, S7, S9, S10, S12, S13, S14, S18, S21): Macro-region × year poststratification calibration weights (`config/calibration_weights.R`, matching submitted S2 Table; public aggregate disclosure explicitly authorised), calibrating pooled estimates from participating hospitals to national macro-region live-birth shares within each year. This does not establish national probability representativeness or eliminate hospital-selection bias. In vitro coverage does not estimate clinical treatment efficacy.

S18 combines isolate-level microbiology/coverage (calibrated weights) with infant-level treatment/outcomes (unit weights); S19 and S21 retain the analysis unit and weighting of each component rather than applying one table-wide weight.

---

## 2. Two Distinct Intrinsic Resistance Rule Sets (`config/intrinsic_resistance_rules.R`)

Two separate, locked rule sets are intentionally maintained and must never be conflated:

1. **AST Reporting Intrinsic Resistance (`get_ast_intrinsic_rules()`, 77 pathogen-drug pairs)**:
   - Used in Module `02_pathogens_ast` (Table 3, S5 Table, S13 Table).
   - Suppresses intrinsically resistant organism-agent cells from antimicrobial resistance reporting and excludes intrinsically resistant components before constructing class-level resistance summaries.
2. **In Vitro Coverage Intrinsic Resistance (`get_coverage_intrinsic_rules()`, 83 pathogen-drug pairs)**:
   - Used in Modules `03_regimen_coverage`, `05_cons_sensitivity`, and `06_additional_analyses` (S9, S10, S14, S18, S19, S21 Tables; Fig 2; S2–S3 Figs).
   - Assigns non-susceptibility (`R`) to intrinsically resistant organism-drug pairs under the **Strict-S** two-drug coverage rule: at least one `S` = Covered; both `I`/`R` = Uncovered; no `S` with any missing component = Unknown; fungi (`n = 20`) are treated as structural Uncovered for antibacterial combinations while retaining raw single-agent AST.
   - Rules follow CLSI M100 (34th ed.) and EUCAST Expected Resistant Phenotypes (v1.2).

---

## 3. Complete Study-Wide Benjamini–Hochberg (BH) FDR Family Registry

All BH-FDR adjustments are computed from **full double-precision unadjusted P values** (never from rounded display P values) directly inside their respective analysis modules. For F1–F13 and other explicitly planned families, when hypotheses are withheld or non-estimable (`NA`), the planned family size `M` is retained via `stats::p.adjust(p_valid, method = "BH", n = M)` (equivalently padding absent tests with `1.0` without displaying P/q for withheld cells). Existing Module 02 families use the locked estimable-endpoint policy: default `p.adjust` excludes NA P values. S19 clinical missingness likewise adjusts estimable tests; coverage missingness retains six slots. These policies were audited and were not redefined.

| Family ID | Final Table | Module | Planned Scope | Planned `M` |
|---|---|---|---|---|
| **F1** | Table 1 | `01_baseline_incidence` | All macro-region omnibus comparisons including EOS incidence | 23 |
| **F3** | S11 Table | `01_baseline_incidence` | Overall EOS incidence temporal trend | 1 |
| **F4** | S11 Table | `01_baseline_incidence` | Four regional EOS incidence temporal trends | 4 |
| **F5** | S11 Table | `01_baseline_incidence` | Overall baseline and early clinical characteristic trends | 34 |
| **F6** | S11 Table | `01_baseline_incidence` | Four regions jointly: baseline and early clinical characteristic trends | 136 |
| **F2** | Table 4 | `04_clinical_outcomes` | All macro-region omnibus comparisons (management & outcomes) | 24 |
| **F7** | S15 Table | `04_clinical_outcomes` | Main overall postnatal antimicrobial regimen trends | 34 |
| **F8** | S15 Table | `04_clinical_outcomes` | Main four-region postnatal antimicrobial regimen trends | 136 |
| **F9** | S15 Table | `04_clinical_outcomes` | Supplementary overall antimicrobial trends | 3 |
| **F10** | S15 Table | `04_clinical_outcomes` | Supplementary four-region antimicrobial trends | 12 |
| **F11** | S16 Table | `04_clinical_outcomes` | Overall clinical outcome and complication trends | 13 |
| **F12** | S16 Table | `04_clinical_outcomes` | Four regions jointly: clinical outcome and complication trends | 52 |
| **F13** | S16 Table | `04_clinical_outcomes` | Four gestational-age groups jointly: clinical outcome trends | 52 |
| **M2-S7** | S7 Table | `02_pathogens_ast` | Pathogen comparisons across Region, GA, Delivery, and <35w vs ≥35w (per dimension, 12 taxa) | 12 taxa per block; adjustment over estimable P values |
| **M2-S12** | S12 Table | `02_pathogens_ast` | Pathogen distribution shift (2020–2024 vs 2015–2019), weighted & unweighted | Estimable taxa per block × level × weighting |
| **M2-S13** | S13 Table | `02_pathogens_ast` | Temporal trends in AMR (within reporting scope) and exploratory ESBL positivity | Estimable agents per pathogen × weighting; exploratory ESBL uses 2 taxa per weighting |
| **M3-S10** | S10 Table | `03_regimen_coverage` | Four primary same-isolate paired coverage contrasts (`C1–C4`) per weighting/sensitivity; plus 32 exploratory subgroup paired contrasts | 4 (primary) / 32 (subgroup) |
| **M3-S14** | S14 Table | `03_regimen_coverage` | Temporal trends in empirical coverage: China regimens (36) and WHO regimens (18), separately for S3 and unit weights | 36 / 18 |
| **M5-S18** | S18 Table | `05_cons_sensitivity` | CoNS sensitivity: paired contrasts (4), treatment comparisons (40), clinical outcomes (13 primary / 2 sensitivity) | 4 / 40 / 13 primary / 2 sensitivity |
| **M6A-S17** | S17 Table | `06_additional_analyses` | Paired CAZ/CTX AST: treatment-process endpoints pooled across both populations, separately from clinical endpoints | 8 treatment-process / 4 clinical |
| **M6B-S19** | S19 Table | `06_additional_analyses` | Missingness Rao–Scott tests across calendar year, macro-region, and outcome (clinical & coverage families) | Per dimension & weighting |

---

## 4. Third-Party Data & Geographic Sources (S1 Fig)

The S1 Fig basemap (`R/07_figures/s1_figure/prepare_geography.R`) is constructed from:
- **Natural Earth** (1:10m cultural/physical country boundaries, provincial boundaries, land, and coastlines): Public Domain ([terms of use](https://www.naturalearthdata.com/about/terms-of-use/)).
- **Supeset/China-GeoData** (supplementary maritime-line and island features, commit `5822c4c0a0bdfd73327f9454976c8661bfd6ad9f`): MIT License.
- Provincial infant/isolate/hospital counts are **not** bundled in this repository; `plot_s1_figure.R` aggregates them dynamically at runtime from authorised `cohort_data` and `hospital_province_map`.

## Historical and redistribution limits

The original dated protocol and any historical SAP require separate author verification. This code records the locked final analysis, not a claim that every inferential method was specified before data inspection. Method rules reproduce two distinct author-supplied locked workbooks; the audit established numeric/mapping equivalence, while redistribution provenance for standards-derived material remains to be confirmed (see THIRD_PARTY_NOTICES.md).

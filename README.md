# CANEOS Statistical Analysis and Figure Generation Code (v1.0.0)

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

Reproducible R analysis, site-selection ordering, and figure generation code for the **China Antimicrobial Resistance and Neonatal Early-Onset Sepsis (CANEOS)** multicenter cohort study (2015–2024), prepared for submission to *PLOS Medicine*.

> **Data Privacy Notice**: This repository contains pure-R statistical code, methodological rule configurations, the 40 published macro-region × year aggregate calibration cells matching `S2 Table`, and documentation only. No patient-level records, real hospital names, or hospital-level denominators are included.

---

## Repository Architecture

- `R/00_site_selection/`: Standalone stratified random priority ordering for the 94 anonymised candidate SiteCodes (`S1 Appendix`; requires zero restricted data).
- `R/01_baseline_incidence/`: Hospital-based EOS incidence, baseline characteristics, temporal trends, and BH-FDR (`Table 1`, `S1–S3 Tables`, `S11 Table`).
- `R/02_pathogens_ast/`: Pathogen distribution, ESBL positivity, antimicrobial susceptibility, subgroup comparisons, and temporal trends (`Table 2`, `Table 3`, `S4`, `S5`, `S7`, `S12`, `S13 Tables`).
- `R/03_regimen_coverage/`: Empirical two-drug coverage, same-isolate paired contrasts, sharp logical bounds, temporal trends, formatted `S2_Heatmap.xlsx`, and four 600-dpi `S3 Fig` TIFF panels (`S9`, `S10`, `S14 Tables`, `Fig 2` data, `S2 Fig`, `S3 Fig`).
- `R/04_clinical_outcomes/`: Postnatal antimicrobial therapy, clinical outcomes, risk ratios, temporal trends with BH-FDR, and dynamic `Figure 3` three-phase treatment trajectories (`Table 4`, `S6`, `S8`, `S15`, `S16 Tables`).
- `R/05_cons_sensitivity/`: CoNS-related sensitivity analyses (`S18 Table`).
- `R/06_additional_analyses/`: Paired CAZ–CTX AST vs outcomes (`S17 Table`), clinical & coverage missingness (`S19 Table`), achieved design effects (`S20 Table`), and Site 23 exclusion sensitivity (`S21 Table`).
- `R/07_figures/`: Publication figures (`Figure 2`, `Figure 3`, `S1 Fig`).
- `config/`: Path template (`paths.example.R`), required output checks (`output_contract.R`), and pure-R static methodological tables (`calibration_weights.R`, `intrinsic_resistance_rules.R`, `field_mappings.R`).
- `docs/`: Complete English documentation (`INPUT_SCHEMA.md`, `OUTPUT_MAP.md`, `ANALYSIS_METHODS.md`, `RESTRICTED_DATA_ACCESS.md`, `THIRD_PARTY_NOTICES.md`).
- `tests/smoke_test.R`: Zero-restricted-data software integrity and standalone `Module 00` verification script.

---

## Setup and Execution

Use **R 4.4.1** and restore exact package versions with `renv::restore()` (using `renv.lock`), or install core CRAN packages via `Rscript --vanilla install_packages.R`.

### 1. Public Software Check & Standalone Site-Selection Run (No Restricted Data Required)

```bash
Rscript --vanilla tests/smoke_test.R
Rscript --vanilla R/00_site_selection/run.R
```

### 2. Full Pipeline Reproduction (With Authorised Restricted Inputs)

1. Obtain the restricted inputs documented in [docs/INPUT_SCHEMA.md](docs/INPUT_SCHEMA.md) through the controlled data access process described in [docs/RESTRICTED_DATA_ACCESS.md](docs/RESTRICTED_DATA_ACCESS.md).
2. Copy `config/paths.example.R` to a private file outside this repository and map the logical aliases to your local input files (original filenames do **not** need to be renamed). Supply it via `--config` or the `CANEOS_CONFIG` environment variable.
3. Build the external `S1 Fig` geography GeoPackage from the pinned public geographic sources specified in `R/07_figures/s1_figure/prepare_geography.R`:

```bash
Rscript --vanilla R/07_figures/s1_figure/prepare_geography.R "/path/to/raw_geography" "/path/to/S1_geography.gpkg"
```

4. Run the full pipeline (`00`–`07`) into an external output directory (defaults to a sibling `CANEOS_restricted_outputs` folder outside this repository, or set `CANEOS_OUTPUT_DIR`):

```bash
Rscript --vanilla run_all.R --config "/path/to/paths.local.R"
```

Use `--modules "01,03"` to run a subset, or invoke any module directly via `R/<module>/run.R`.

---

## Outputs, Privacy & Interpretation

- See [docs/OUTPUT_MAP.md](docs/OUTPUT_MAP.md) for the complete table/figure-to-script mapping.
- Runtime outputs contain study derivatives and must always remain outside any public archive.
- Geographic stratification and macro-region × year live-birth calibration align pooled estimates across the 30 participating tertiary hospitals and do not establish a national probability sample or eliminate hospital-selection bias. In vitro coverage estimates do not equate to clinical efficacy. See [docs/ANALYSIS_METHODS.md](docs/ANALYSIS_METHODS.md) and [docs/THIRD_PARTY_NOTICES.md](docs/THIRD_PARTY_NOTICES.md) for methodological details and third-party attributions.

---

## License

Author-generated analysis code is released under the [MIT License](LICENSE) (`Copyright (c) 2026 CANEOS Study Group`).

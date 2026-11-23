# Third-Party Notices and Attribution

## 1. Author-Generated Analysis Code

All author-generated R scripts and documentation in this repository are distributed under the [MIT License](../LICENSE) (`Copyright (c) 2026 CANEOS Study Group`).

## 2. Geographic Basemap Sources (`S1 Fig`)

No third-party geographic binary datasets are bundled in this repository. `R/07_figures/s1_figure/prepare_geography.R` reconstructs the `S1 Fig` GeoPackage from the following external public sources (verifying pinned input hashes):

- **Natural Earth** (1:10m cultural and physical country boundaries, provincial boundaries, land, and coastlines): Public Domain ([Terms of Use](https://www.naturalearthdata.com/about/terms-of-use/)).
- **Supeset/China-GeoData** (supplementary maritime-line and island features, commit `5822c4c0a0bdfd73327f9454976c8661bfd6ad9f`): [MIT License](https://raw.githubusercontent.com/Supeset/China-GeoData/5822c4c0a0bdfd73327f9454976c8661bfd6ad9f/LICENSE) (`Copyright (c) 2025 圈集`; [repository](https://github.com/Supeset/China-GeoData/tree/5822c4c0a0bdfd73327f9454976c8661bfd6ad9f)). Retain upstream license notices if redistributing any derived geographic files.

## 3. Antimicrobial Susceptibility & Intrinsic Resistance References

`config/intrinsic_resistance_rules.R` encodes the study's analytical pathogen-by-agent classification pairs (77 AST reporting exclusion rules, 83 coverage intrinsic-resistance rules, and 8 Enterobacterales CAZ/CTX applicability entries). No copyrighted CLSI or EUCAST PDFs or full guideline tables are bundled. Analytical classifications follow:

- **CLSI**: Clinical and Laboratory Standards Institute. *Performance Standards for Antimicrobial Susceptibility Testing*, 34th ed. CLSI supplement M100. Wayne, PA: CLSI; 2024 ([CLSI Copyright Policy](https://support.clsi.org/what-is-clsis-copyright-policy)).
- **EUCAST**: European Committee on Antimicrobial Susceptibility Testing. *Expected Resistant Phenotypes*, Version 1.2 (January 2023) ([EUCAST Expected Phenotypes](https://www.eucast.org/bacteria/important-additional-information/expected-phenotypes/)).

# Controlled Data Access

## Restricted Study Inputs

Individual-level infant and isolate records, hospital-specific annual neonatal admission denominators, hospital-to-anonymised-Site-ID crosswalks, hospital-to-province crosswalks, and the isolate-level coverage scope sidecar contain sensitive institutional or clinical records and cannot be redistributed publicly under multicenter ethics and data-governance agreements.

Only the 40 macro-region × year aggregate calibration cells (published in `S2 Table`) and non-record methodological lookup configurations are included in this repository.

## How Qualified Researchers Can Request Access

Qualified researchers may submit data-access requests through the institutional ethics committee and data-governance process at the:

> **National Clinical Research Center for Children and Adolescents' Health and Diseases**  
> **Children's Hospital of Chongqing Medical University**

Please refer to the **Data Availability Statement** in the accompanying *PLOS Medicine* manuscript for the official institutional contact address and submission channel. Requests should outline the research objectives, ethics review status at the requesting institution, and data-security arrangements, and are subject to institutional ethics approval and execution of a Data Use Agreement.

## Testing Without Restricted Data

To verify script syntax, R package availability, static methodological configurations, BH-FDR helper logic, UTF-8 Chinese CSV reading, and standalone site-selection ordering (`Module 00`) without restricted inputs, run:

```bash
Rscript --vanilla tests/smoke_test.R
```

Note that `tests/smoke_test.R` uses temporary synthetic inputs for software verification only; reproducing the clinical and microbiological tables and figures requires the authorised inputs documented in [INPUT_SCHEMA.md](INPUT_SCHEMA.md).

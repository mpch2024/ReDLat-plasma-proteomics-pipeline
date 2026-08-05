# DEP refactor changelog

## Structural changes

- Renamed the 12 active scripts into a single canonical execution order.
- Removed the exact duplicate sensitivity script.
- Moved the superseded fragmented Supplementary Tables generator to `archive/DEP_legacy/`.
- Added clean-process runners for analysis, reporting and the complete DEP workflow.

## Reproducibility changes

- Replaced local path discovery and directory inference with one configuration object.
- Removed all automatic package installation calls; the pipeline now expects `renv::restore()`.
- Removed `setwd()` and Windows/user-specific absolute paths.
- Fixed the analysis-root inference used by the sensitivity and robustness scripts.
- Disabled duplicated country/site robustness calculations in Script 01; these are performed by Script 04.
- Enforced the expected 9,638-row fixed gene–SOMAmer map in the fixed-map rebuild.

## Privacy changes

- Raw inputs are supplied only through local configuration or environment variables.
- Participant-set exports are routed to `result/private/participant_sets/`.
- Identifier-bearing PCA outputs are routed to `result/private/pca/`.
- Reporting uses a privacy-safe PCA copy without stable identifiers.
- Source Data workbooks stop if direct or stable identifier columns are detected.
- `data_private/`, `result/` and `publication_candidate/` are ignored by Git.
- Added repository and publication-candidate audits for identifiers and local absolute paths.

## Scientific scope

No statistical formula, model coefficient, FDR threshold, p-tau217 quantile, matching specification, fixed-map selection rule or figure artboard was intentionally changed in this refactor.

## Version 2
- Replaced inherited development headers with concise professional script headers.
- Removed review-stage and version-history language from introductory titles.
- Removed development-only attribution and retained only scientific and technical documentation.

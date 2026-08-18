# ReDLat plasma proteomics pipeline

Reproducible analysis code for the study **“Plasma proteomics maps the molecular architecture of Alzheimer’s disease in Latin America”**, including differential protein abundance analysis (DEP), weighted gene co-expression network analysis (WGCNA), sensitivity analyses, and version-locked manuscript-reproduction machine-learning workflows. A separate strict-CV reference implementation is provided in `python/redlat_ml/strict_cv.py`.

This repository contains analysis code, configuration helpers, environment specifications, validation utilities, and privacy audits. **Participant-level ReDLat data are not included in the public repository.** Access to governed data remains subject to ReDLat ethics approvals, informed-consent conditions, and data-use agreements.

## Processed reviewer-data interface

The codebase has been validated against a minimal processed and pseudonymized reviewer dataset consisting of exactly three files:

- `ReDLat_proteomics_somamer_log2.csv` — 653 participants × 10,751 human SOMAmer features, represented as log2-transformed normalized RFU values;
- `ReDLat_metadata_deidentified.csv` — pseudonymized demographic, clinical, biomarker, and analysis-inclusion variables required by the reported workflows;
- `ReDLat_feature_annotation.csv` — SOMAmer/protein annotation and the feature-selection mappings required to reconstruct the 9,638-gene analysis universes.

These files are **not committed to GitHub**. For an authorized local rerun, place them under:

```text
data_private/reviewer_inputs/
```

or point `REDLAT_REVIEWER_DATA_DIR` to an approved local directory containing the same three files.

The processed-data compatibility layer does not require or recreate the original SomaScan ADAT. For legacy R steps that historically performed a single log2 transformation, the adapter supplies `2^(released log2 RFU)` in memory so the original numerical workflow is preserved without double transformation. The historical machine-learning master is reconstructed locally from the released feature map and normalized-RFU-compatible values; no raw participant linkage file is required.

## Repository structure

| Path | Description |
|---|---|
| `scripts/DEP/` | Differential abundance, sensitivity, robustness, reporting, and figure-generation scripts |
| `scripts/WGCNA/` | Network construction, module association, preservation, reporting, and audit scripts |
| `scripts/ML/` | Nested-CV, biomarker, APOE, LOCO, matching, reporting, and publication-audit scripts |
| `scripts/validation/` | Reviewer-data preflight, locked-count checks, static privacy audit, and dependency map |
| `R/` | Shared R bootstrap utilities and processed-data compatibility layer |
| `python/redlat_ml/` | Shared Python configuration and privacy utilities |
| `config/` | Public project-relative configuration |
| `environment/` | Package requirements used for rerun environments |
| `data_private/reviewer_inputs/` | Untracked local location for governed reviewer inputs |
| `outputs/` | Generated local analysis products; excluded from version control |

Existing repository-level workflows, tests, CI configuration, license files, and historical utilities remain part of the repository unless explicitly superseded by a later validated update.

## Environment

The manuscript analyses were developed and executed using **R 4.5.2** and **Python 3.13.11**. Required Python packages are listed in `environment/requirements.txt`; the R packages used by the reviewer-data rerun are listed in `environment/R_PACKAGES_REQUIRED.txt`. The repository-level `renv.lock` and any existing environment files should be retained as the historical dependency record.

Scripts do not install or update packages automatically.

## First validation step

After placing the three authorized reviewer files in `data_private/reviewer_inputs/`, run from the repository root:

**Windows**

```bat
RUN_PRECHECK.bat
```

**macOS/Linux**

```bash
bash RUN_PRECHECK.sh
```

The precheck:

1. validates the three-file data contract and reconstructs local ML-compatible inputs;
2. checks the repository for personal absolute paths, `setwd()` calls, active ADAT reads, runtime package installation, and transient patch artifacts;
3. verifies locked manuscript counts and feature dimensions;
4. generates a script dependency map.

Expected locked quantities include 653 total participants, 639 primary complete cases, 579 APOE-eligible participants, 343 AT(N)-eligible participants, the selected 191 CN + 191 AD matched cohort, 10,751 SOMAmers, and 9,638 gene-level features.

## Running the analyses

The validated direct execution order is recorded in `RUN_ANALYSES_IN_ORDER.txt`. In brief:

1. run the DEP scripts in numerical order;
2. run WGCNA scripts in numerical order after DEP inputs are available;
3. run the ML R/Python scripts in numerical order after `RUN_PRECHECK` has generated the local ML compatibility files.

The ML workflow preserves the version-locked implementation used to reproduce the reported results. A separate strict-CV reference implementation is provided in `python/redlat_ml/strict_cv.py`, where preprocessing is encapsulated within training CV pipelines. Country-held-out analyses are interpreted as internal geographic robustness rather than external validation.

## Validation status

The processed reviewer-data interface was audited against the manuscript-level analytical quantities. The validation included:

- exact primary cohort and feature counts;
- recovery of the submitted seven-protein ML panel (`SPC25`, `CPLX2`, `TCP11L1`, `ACHE`, `ODC1`, `SPON1`, `RTN4RL1`);
- APOE, p-tau217, LOCO, matched-cohort, and LOCO meta-analysis reruns;
- WGCNA construction, association, stability, and structural-preservation workflows;
- static privacy checks and publication-output audits.

See `docs/FINAL_CODE_VALIDATION.md` for the code-update audit and `docs/REVIEWER_DATA_INTERFACE.md` for the local data contract.

## Data governance

Do not commit participant-level data, pseudonymous IDs, linkage keys, raw ADAT files, local absolute paths, credentials, or private analysis outputs. The public repository contains code only. Authorized reviewer datasets should be distributed through the controlled/private mechanism approved for the peer-review process.

## License

The analysis code is released under the MIT License. The existing repository-level `LICENSE` file should be retained.

## Contact

For questions regarding the study, data access, or analytical workflow:

**Claudia Durán-Aniotz**
Universidad Adolfo Ibáñez
Email: claudia.duran@uai.cl

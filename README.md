# ReDLat plasma proteomics pipeline

Reproducible analysis code for the ReDLat plasma proteomics study, including differential-abundance analysis (DEP), weighted gene co-expression network analysis (WGCNA), and leakage-controlled machine-learning workflows.

This repository contains analysis code, workflow runners, environment specifications, documentation, tests, and static audits. Participant-level clinical, biomarker, genotype, and SOMAscan data are not included and remain subject to ReDLat governance, informed consent, local ethics approvals, and data-use agreements.

## Repository structure

| Path | Description |
|---|---|
| `scripts/DEP/` | Primary differential-abundance, sensitivity, robustness, reporting, and figure-generation scripts |
| `scripts/WGCNA/` | WGCNA network construction, association, preservation, reporting, and audit scripts |
| `scripts/ML/` | Strict nested-cross-validation, biomarker, APOE, LOCO, matching, reporting, and audit scripts |
| `workflow/` | Module-level analysis and reporting runners, environment capture, and release audits |
| `R/` | Shared R bootstrap utilities |
| `python/redlat_ml/` | Shared Python utilities for strict machine-learning workflows |
| `config/` | Public configuration templates |
| `environment/` | Python environment specifications |
| `docs/` | Pipeline maps, input contracts, reproducibility records, manifests, and audit documentation |
| `tests/` | DEP infrastructure and machine-learning smoke tests |
| `data_private/` | Placeholder for governed local inputs; participant-level data must never be committed |
| `result/` | Local private analysis outputs; excluded from version control |
| `publication_candidate/` | Privacy-audited release candidates; excluded from version control |

## System requirements

The workflows use R, Python, and standard command-line tools. No non-standard hardware is required by the code, although the full WGCNA, nested-cross-validation, and simulation workflows can be computationally intensive.

### R environment

The authoritative R dependency record is the repository-level `renv.lock` file. From the repository root, restore the environment with:

```r
install.packages("renv")
renv::restore()
```

Scripts do not install or update packages automatically.

### Python environment

A conservative cross-platform rerun environment is specified in `environment/environment.yml`:

```bash
conda env create -f environment/environment.yml
conda activate redlat-ml-strict
```

Alternatively:

```bash
python -m pip install -r environment/requirements.txt
```

The dependency ranges in `environment/` are installation constraints for the strict rerun environment and should not be interpreted as an exact record of the historical submission environment. Local software versions can be captured with:

```bash
python workflow/00_capture_software_environment.py
```

## Local configuration

Copy `.Renviron.example` to an untracked local `.Renviron` file and enter only governed local paths:

```bash
cp .Renviron.example .Renviron
```

On Windows, the file may also be copied manually in File Explorer.

The canonical private inputs are:

- harmonized clinical metadata;
- SOMAscan ADAT data;
- analysis-specific biomarker, APOE, exclusion, or matching files where required.

The required variables and accepted fields are documented in `config/README.md` and `docs/DEP/INPUT_DATA_CONTRACT.md`.

Never place participant identifiers, credentials, governed exclusions, populated local paths, or raw data in tracked source files.

## Running the workflows

Run commands from the repository root. The environment variable `REDLAT_PROJECT_ROOT` may be used to define the repository location; otherwise, the workflow runners resolve the current project root.

### Differential-abundance analysis

Run the primary DEP analyses:

```bash
Rscript --vanilla workflow/run_dep_analysis.R
```

Run the reporting workflow:

```bash
Rscript --vanilla workflow/run_dep_reporting.R
```

Run the complete DEP workflow:

```bash
Rscript --vanilla workflow/run_dep_all.R
```

The retrospective simulation-based detectable-effect analysis is optional because of its computational cost. It is not an a priori sample-size calculation and is not run by default. To include it explicitly:

```bash
Rscript --vanilla workflow/run_dep_analysis.R --with-simulation
```

or:

```bash
Rscript --vanilla workflow/run_dep_all.R --with-simulation
```

The canonical DEP order and scientific role of each script are documented in `docs/DEP/DEP_PIPELINE_MAP.md`.

### WGCNA

Run the analysis workflow:

```bash
Rscript --vanilla workflow/run_wgcna_analysis.R
```

Run reporting and submission-package generation:

```bash
Rscript --vanilla workflow/run_wgcna_reporting.R
```

Run the complete WGCNA workflow:

```bash
Rscript --vanilla workflow/run_wgcna_all.R
```

The network is constructed before clinical association testing. Country and site analyses assess internal stability and are not interpreted as external validation.

### Machine learning

Run the strict analysis workflow:

```bash
python workflow/run_ml_analysis.py
```

Run reporting and privacy audits:

```bash
python workflow/run_ml_reporting.py
```

Run the complete machine-learning workflow:

```bash
python workflow/run_ml_all.py
```

The strict workflow performs fold-specific differential-abundance filtering, preprocessing, feature selection, and hyperparameter optimization within training data to prevent information leakage. The analytical sequence is documented in `docs/ML/ML_PIPELINE_MAP.md`.

## Tests and static audits

Run the DEP infrastructure test:

```bash
Rscript --vanilla tests/test_dep_infrastructure.R
```

Run the synthetic machine-learning smoke test:

```bash
python tests/test_strict_cv_smoke.py
```

Run repository and release-candidate audits:

```bash
Rscript --vanilla workflow/audit_dep_repository.R
Rscript --vanilla workflow/audit_publication_candidate.R
Rscript --vanilla workflow/audit_wgcna_publication.R
python scripts/ML/25_ML_audit_strict_pipeline.py
```

GitHub Actions also provide static checks for DEP, WGCNA, and machine-learning source files. Static checks do not access or execute governed participant-level data.

## Expected outputs

- `result/`: private intermediate objects, model outputs, workflow logs, and participant-level analytical products;
- `publication_candidate/DEP/`: DEP figures, Source Data, Supplementary Tables, and Supplementary Data awaiting privacy and consistency review;
- `publication_candidate/WGCNA/`: WGCNA release candidates and submission-package outputs;
- `publication_candidate/ML/`: machine-learning figures, Source Data, and audited aggregate outputs.

Generated outputs are excluded from version control. Nothing is promoted from `publication_candidate/` to a public repository or archive automatically. Every release candidate requires independent scientific, consistency, and privacy review.

## Runtime

Runtime depends on the governed dataset, hardware, operating system, and selected workflow. Full nested-cross-validation, WGCNA preservation, and the optional DEP simulation can require substantially longer than static audits or smoke tests. No universal runtime benchmark is claimed; users should record elapsed time and software information for each validated local run.

## Data availability and governance

This repository does not provide access to restricted ReDLat data and cannot reconstruct participant identities. Participant-level data are stored only in governed local environments. Access requests are evaluated under the ReDLat consortium's ethics approvals, informed-consent conditions, and data-use agreements.

The privacy and release model is documented in:

- `docs/DEP/PRIVACY_AND_RELEASE.md`;
- `docs/DEP/OUTPUT_CLASSIFICATION.md`;
- the module-specific audit reports under `docs/`.

## Reproducibility notes

- Canonical scripts are executed in clean processes by the workflow runners where specified.
- Prespecified random seeds are retained within the originating scripts.
- Local inputs are configured explicitly; scripts do not search personal desktops or historical project folders.
- Country- and site-based analyses are internal robustness assessments rather than independent external validation.
- The optional DEP simulation is a retrospective operating-characteristic analysis and must not be described as prospective power determination or observed post hoc power.

## License

The analysis code is released under the MIT License. See `LICENSE`.

## Code repository

https://github.com/mpch2024/ReDLat-plasma-proteomics-pipeline

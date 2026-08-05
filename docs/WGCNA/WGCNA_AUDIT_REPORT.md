# WGCNA workflow audit

## Canonical decisions

- The workbook-corrected structural-preservation script is retained.
- The earlier structural-preservation script is archived.
- The no-plotmath Extended Data generator is retained.
- The compact-classic v6 generator is archived.
- The historical delivery runner is archived and replaced by analysis and reporting runners.

## Reproducibility corrections

- Removed personal Windows paths and external folder names.
- Added a local, git-ignored configuration file.
- Removed automatic CRAN and Bioconductor installation.
- Removed working-directory changes from ZIP creation.
- Centralized analysis outputs under `result/WGCNA/`.
- Centralized publication outputs under `publication_candidate/WGCNA/`.
- Added package and publication-tree privacy audits.
- Added explicit analysis and reporting runners.

## Privacy classification

The following outputs are private intermediates and must not be committed or deposited publicly:

- gene-collapsed participant-level expression matrices;
- WGCNA participant metadata;
- per-participant module eigengenes;
- participant-removal and matching tables;
- analysis workspaces containing participant-level objects.

Only audited aggregate tables, figures and Source Data without direct identifiers may enter the submission package.

## Scientific scope

The transformation does not intentionally change network parameters, statistical formulas, model families, FDR definitions, module assignments or figure dimensions. Full numerical validation requires execution in the locked R environment with the governed ReDLat inputs.

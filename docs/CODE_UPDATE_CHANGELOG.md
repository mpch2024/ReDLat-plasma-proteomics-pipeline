# Processed reviewer-data code update changelog

## Shared data interface

- Added `R/reviewer_data_adapter.R` to load the three processed/pseudonymized reviewer inputs without requiring the original ADAT.
- Added project-relative public configuration for DEP, WGCNA and ML.
- Added validation/precheck scripts that reconstruct local ML-compatible derivatives and verify locked manuscript counts.

## DEP

- Consolidated the validated DEP scripts 01–12.
- Retained explicit namespace use where required to avoid masked `rename()` behavior.
- Preserved the fixed feature-selection mapping used in sensitivity and matching workflows.

## WGCNA

- Consolidated validated input-selection handling, workspace validation, structural-preservation checkpoint reuse, optional historical-comparison handling, and final figure/package references.

## Machine learning

- Preserved the historical 9,638-gene ML feature mapping through the controlled feature-annotation flag used by the local precheck.
- Standardized p-tau217 compatibility and pseudonymous `Study_ID`/`SampleId` handling.
- LOCO country membership is read from the explicit `Country` metadata field rather than inferred from original identifier prefixes.
- The matched sensitivity audits/restores the exact prespecified pseudonymized 191 CN + 191 AD cohort instead of silently rematching a tie-sensitive solution.
- Consolidated fixed-panel, LOCO meta-analysis, output-directory, plotting, and publication-output privacy fixes.
- Added public Source Data guards against participant identifiers, ReDLat pseudonyms, and personal filesystem paths.

## Public-repository boundary

- Reviewer CSVs, pseudonymous IDs, raw ADAT files, linkage keys, private outputs, and reviewer sharing links remain outside Git.
- Existing public workflow runners, tests, CI, license, renv records, DEP script 13, and ML script 26 remain preserved.
- ML script 25 was superseded after CI review so that it blocks public-repository safety/portability violations while transparently reporting version-locked manuscript-reproduction implementation signatures without changing the model-fitting scripts.

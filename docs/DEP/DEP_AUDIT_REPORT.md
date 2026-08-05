# DEP pipeline refactor audit report

## Status

The DEP module has been reorganized into 12 canonical scripts, three clean-process runners, shared configuration, privacy controls, static CI and release-candidate audits.

## Passed static checks

- 12 canonical analysis/reporting scripts present in one execution order.
- No `install.packages()` or `BiocManager::install()` calls in canonical code.
- No `setwd()` calls.
- No user-specific Windows, macOS or Linux home-directory paths.
- No original private metadata or ADAT filenames embedded in canonical code.
- Direct participant-set exports are routed to git-ignored private directories.
- All canonical R files passed a delimiter/string lexical balance check.
- Input/output folder references were reconciled across the 12 stages.
- The exact duplicate Script 02 was removed; the competing legacy Script 08 was archived.

## Scientific changes intentionally made

- Script 01 no longer reruns country/site robustness analyses already implemented canonically in Script 04.
- The fixed-map rebuild now requires the expected 9,638 gene–SOMAmer pairs rather than treating this check as optional.

No model formula, FDR threshold, p-tau217 quantile, matching specification, effect definition or figure dimension was intentionally changed.

## Privacy controls

- Source Data generation is blocked by default and requires an explicit, untracked governance flag.
- Generated data and publication candidates are excluded from Git.
- Source Data workbooks are checked for direct identifier columns before writing.
- A separate workbook/CSV audit checks publication candidates for identifier tokens and local absolute paths.

## Validation limitation

This environment does not contain R, the locked package library or restricted ReDLat inputs. Therefore, the full pipeline could not be executed end to end here. Static syntax structure, cross-script paths, privacy patterns and file organization were audited. Before replacing the GitHub branch, run the pipeline on a clean local clone using the repository `renv.lock` and compare all canonical numerical checkpoints with the frozen Supporting Information package.

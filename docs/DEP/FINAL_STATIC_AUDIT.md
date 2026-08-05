# DEP final static audit

## Status

**Conditionally ready for local execution.** The canonical scripts are structurally organized, use relative or privately configured paths, and separate participant-level outputs from publication candidates.

## Passed checks

- Twelve canonical scripts are present in the documented order.
- No legacy scripts are included in the release package.
- No user-specific absolute path, personal identifier or development attribution is embedded in canonical code.
- Runtime package installation and `setwd()` are absent.
- Direct participant identifiers are written only under the private result tree.
- Participant-level Source Data export is disabled by default and requires an untracked governance setting.
- The environment-capture utility records runtime versions and source hashes without committing local paths.
- Static delimiter and quote balance passed for all R files.

## Required before release

- Parse and execute all R scripts in the locked local R environment.
- Run `renv::restore()` from the repository-level lock file.
- Run `workflow/00_capture_software_environment.py` immediately before the final production run.
- Compare numerical checkpoints and generated workbooks with the frozen Supporting Information.
- Run both repository and publication-candidate privacy audits.

No scientific result was recomputed during this static audit.

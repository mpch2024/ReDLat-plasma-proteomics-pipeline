# WGCNA final static audit

## Status

**Conditionally ready for local execution.** The active workflow is organized into analysis and reporting stages and no longer contains legacy files with personal paths or development annotations.

## Passed checks

- Fifteen canonical scripts are present in the documented order.
- Legacy preservation, figure and delivery scripts were removed from the release package.
- No user-specific absolute path, participant identifier or development attribution is embedded in canonical code.
- Runtime package installation and `setwd()` are absent.
- Analysis outputs remain under the local WGCNA result root; publication artifacts are built separately.
- The environment-capture utility records R, system-tool and source-file versions for the internal archive.
- Static delimiter and quote balance passed for all R files.

## Required before release

- Parse and execute all R scripts in the locked local R environment.
- Restore the repository-level `renv.lock` and record the final R/Bioconductor environment.
- Rebuild the submission package from a clean result directory.
- Run Scripts 15 and `workflow/audit_wgcna_publication.R` and inspect every workbook for indirect re-identification risk.
- Compare network, module, preservation and figure checkpoints with the frozen Supporting Information.

No WGCNA model was recomputed during this static audit.

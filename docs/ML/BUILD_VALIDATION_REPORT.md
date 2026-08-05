# Build validation report

## Completed in the build environment

- All Python files compile successfully.
- The strict implementation audit passes.
- A synthetic missing-data smoke test passes for training-only feature filtering, nested elastic-net selection, pipeline-based imputation and scaling, SVM tuning, probability calibration and held-out permutation importance.
- No active script contains a personal absolute path or automatic package installation, `setwd()`, recursive historical-input discovery or fold-wise score standardization.
- Historical `Result_*` directory overrides were removed from the four fold-specific DEP scripts.
- The matching script now keeps its configuration in scope and writes only to the configured private result root.

## Not executable in the build environment

R is not installed, and governed ReDLat inputs are unavailable. R syntax/runtime validation, full end-to-end execution and numerical checkpoint comparison must therefore be completed locally.

## Release gate

Do not replace the historical manuscript values until all five strict branches have been rerun, aggregate outputs have been reviewed, and the figures, Source Data and Supplementary Data have been regenerated from the strict results.

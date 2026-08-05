# Strict ML audit report

## Corrections implemented

1. Replaced pre-scaled `LogisticRegressionCV` with training-only elastic-net searches whose imputation and scaling are part of the estimator pipeline.
2. Replaced pre-scaled SVM searches with pipelines fitted inside inner cross-validation.
3. Added training-only protein missingness and variance filters.
4. Removed global protein complete-case exclusion from p-tau217, APOE and fixed-panel analyses.
5. Replaced fold-wise z-normalized decision scores with calibrated OOF probabilities.
6. Replaced full-cohort permutation importance with outer-test-fold permutation importance.
7. Added fold-level leakage audit tables and effective sample-size reports.
8. Corrected DEP-fold scripts that previously overwrote configured fold directories with historical paths.
9. Rebuilt the matching script so configuration remains available throughout execution.
10. Removed recursive discovery of private historical inputs from matching.

## Required validation

This package has passed Python compilation, static method checks and a synthetic smoke test. Numerical validation still requires the governed ReDLat data, R, the locked R environment and the configured private inputs. The strict rerun is expected to be compared against the historical Supporting Information; numerical equality is not assumed.

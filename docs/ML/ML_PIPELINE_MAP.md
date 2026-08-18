# Machine-learning pipeline map

Scripts 1–24 form the version-locked manuscript-reproduction path. They preserve
the analytical implementation used to generate the reported results. The
separate `python/redlat_ml/strict_cv.py` module provides reusable stricter-CV
components with preprocessing fitted inside training CV splits; it is not
silently substituted into the manuscript-reproduction path because doing so
would define a different analysis.

| Step | Role | Evaluation / reproducibility scope | Main output |
|---:|---|---|---|
| 1–3 | Primary protein classifier | Outer held-out folds with fold-specific DEP; historical selection, scaling and tuning implementation retained for manuscript reproduction | AUC, recurrent panel, OOF scores |
| 4–6 | Protein, p-tau217 and combined models | Common biomarker-availability cohort defined before CV; outer folds held out; historical model-comparison implementation retained | Paired model comparison and DeLong tests |
| 7–9 | Protein and protein + APOE models | APOE-availability cohort with outer held-out folds; historical model-comparison implementation retained | Paired comparison and DeLong test |
| 10–12 | Country-held-out robustness | Entire country excluded before DEP and model fitting | Country-specific AUCs |
| 13 | Prespecified seven-protein panel | Fixed-panel biomarker comparison; no feature reselection | Fixed-panel biomarker comparison |
| 14 | Clinical alignment | Descriptive analysis separated from classifier evaluation | Regression tables |
| 15 | Matching | Explicit configured inputs; exact prespecified matched cohort can be restored from the controlled inclusion flag | 191 CN + 191 AD matched cohort |
| 16 | Matched sensitivity | Exact matched cohort; historical nested and fixed-panel branches retained for manuscript reproduction | Matched AUCs and OOF scores |
| 17 | LOCO synthesis | Aggregate country-level estimates only | Random-effects synthesis |
| 18–24 | Reporting and privacy audit | Reads completed analysis outputs; no refitting in reporting scripts | Publication-candidate package |
| 25 | Public-repository/static audit | Blocks personal absolute paths and active runtime installers; reports known historical implementation signatures without rewriting analytical scripts | Static audit report |
| 26 | Extended Data Figure 10 generation | Reporting/figure generation from validated matched-analysis outputs | Extended Data Figure 10 |

## Strict-CV reference utility

`python/redlat_ml/strict_cv.py` contains a separate reference implementation
using scikit-learn pipelines so imputation/scaling are fitted within training CV
splits, together with calibrated prediction and fold-level utilities. It is
tested separately and should not be interpreted as a silent replacement for
the version-locked manuscript-reproduction scripts.

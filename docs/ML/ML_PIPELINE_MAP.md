# Strict machine-learning pipeline map

| Step | Role | Leakage control | Main output |
|---:|---|---|---|
| 1–3 | Primary protein classifier | Fold-specific DEP; preprocessing and tuning inside training CV | AUC, recurrent panel, OOF probabilities |
| 4–6 | Protein, p-tau217 and combined models | Biomarker-defined cohort; protein missingness imputed inside CV | Paired model comparison and DeLong tests |
| 7–9 | Protein and protein + APOE models | APOE-defined cohort; preprocessing inside CV | Paired comparison and DeLong test |
| 10–12 | Country-held-out robustness | Entire country excluded before DEP and model fitting | Country-specific AUCs |
| 13 | Prespecified seven-protein panel | Nested tuning; no feature selection | Fixed-panel biomarker comparison |
| 14 | Clinical alignment | Descriptive analysis separated from classifier evaluation | Regression tables |
| 15 | Matching | Explicit configured inputs; no recursive file discovery | 191 CN + 191 AD matched cohort |
| 16 | Matched sensitivity | Training-only DEP and nested preprocessing; fixed-panel branch prespecified | Matched AUCs and OOF probabilities |
| 17 | LOCO synthesis | Aggregate country-level estimates only | Random-effects synthesis |
| 18–24 | Reporting and privacy audit | Reads completed analysis outputs; no refitting | Publication candidate package |
| 25 | Strict implementation audit | Blocks historical leakage patterns and personal paths | Static audit report |

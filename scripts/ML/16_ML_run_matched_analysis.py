###############################################################################
# ReDLat plasma proteomics — strict machine-learning workflow
# 16. Run strict matched-sample sensitivity
# Requires: Script 15 selected matching output and private master matrix
# Produces: reselected-panel, fixed-panel and p-tau217 matched comparisons
# Data policy: participant-level files remain under the private result root.
###############################################################################

from pathlib import Path
import sys
_REPO_HINT = Path(__file__).resolve()
for _candidate in [_REPO_HINT.parent, *_REPO_HINT.parents]:
    if (_candidate / ".redlat-root").exists() or (_candidate / ".here").exists():
        PROJECT_ROOT = _candidate
        break
else:
    raise FileNotFoundError("Repository root not found. Set REDLAT_PROJECT_ROOT.")
sys.path.insert(0, str(PROJECT_ROOT / "python"))


import numpy as np
import pandas as pd
from scipy import stats
from redlat_ml.config import load_config, require_files
from redlat_ml.strict_cv import (
    META_COLUMNS, StrictCVSettings, aggregate_importance, aggregate_panel_frequency,
    delong_pairwise, feature_audit, fit_tuned_calibrated_svm, fold_permutation_importance,
    metric_row, predict_probability, resolve_column, roc_table, safe_n_splits,
    select_stable_panel, settings_table,
)

CONFIG = load_config(__file__); SETTINGS = StrictCVSettings.from_env()
ROOT = CONFIG.private_root / "matched"; ROOT.mkdir(parents=True, exist_ok=True)
RESELECT = ROOT / "02_nested_cv"; FIXED = ROOT / "03_primary_fixed_panel_cv"; PTAU = ROOT / "04_primary_panel_ptau"
for path in (RESELECT, FIXED, PTAU): path.mkdir(parents=True, exist_ok=True)
require_files([(CONFIG.master_file, "ML master matrix"), (CONFIG.matched_ids_file, "selected matched IDs")])
PANEL = ["SPC25", "CPLX2", "TCP11L1", "ACHE", "ODC1", "SPON1", "RTN4RL1"]
master = pd.read_csv(CONFIG.master_file, low_memory=False); matched_ids = pd.read_csv(CONFIG.matched_ids_file)
master["SampleId"] = master["SampleId"].astype(str).str.strip(); matched_ids["SampleId"] = matched_ids["SampleId"].astype(str).str.strip()
matched = matched_ids[["SampleId", "SampleGroup"]].merge(master, on="SampleId", how="left", suffixes=("_match", ""))
if "SampleGroup" not in matched.columns: matched["SampleGroup"] = matched["SampleGroup_match"]
matched = matched.drop_duplicates("SampleId").reset_index(drop=True)
counts = matched.SampleGroup.value_counts()
if int(counts.get("CN", 0)) != 191 or int(counts.get("AD", 0)) != 191: raise RuntimeError(f"Matched cohort must contain 191 CN and 191 AD; observed {counts.to_dict()}")
y = matched.SampleGroup.map({"CN": 0, "AD": 1}).astype(int)
protein_candidates = [c for c in master.columns if c not in META_COLUMNS and c not in {"SampleId", "SampleGroup"}]
protein_candidates, protein_audit = feature_audit(matched, protein_candidates, SETTINGS); protein_audit.to_csv(RESELECT / "candidate_protein_audit.csv", index=False)
outer = list(StratifiedKFold(n_splits=5, shuffle=True, random_state=SETTINGS.random_state).split(matched, y))

def bh(values):
    values = np.asarray(values, dtype=float); order = np.argsort(values); ranked = values[order]; adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1); adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]; out = np.empty_like(adjusted); out[order] = np.clip(adjusted, 0, 1); return out

def training_dep(train: pd.DataFrame, proteins: list[str], fold: int):
    metadata = train[["SampleGroup", "Age", "Sex", "Education", "Country"]].copy()
    if metadata.isna().any(axis=1).any(): raise RuntimeError(f"Fold {fold}: missing DEP covariates in matched training data.")
    design = pd.DataFrame({"Intercept": 1.0, "AD": (metadata.SampleGroup == "AD").astype(float), "Age": pd.to_numeric(metadata.Age), "Education": pd.to_numeric(metadata.Education)}, index=train.index)
    design = pd.concat([design, pd.get_dummies(metadata.Sex.astype(str), prefix="Sex", drop_first=True, dtype=float), pd.get_dummies(metadata.Country.astype(str), prefix="Country", drop_first=True, dtype=float)], axis=1)
    design = design.loc[:, [c for c in design if c == "Intercept" or design[c].nunique(dropna=False) > 1]]
    eligible, _ = feature_audit(train, proteins, SETTINGS); raw = train[eligible].apply(pd.to_numeric, errors="coerce").to_numpy(float); raw[raw <= 0] = np.nan; expression = np.log2(raw)
    medians = np.nanmedian(expression, axis=0); medians[~np.isfinite(medians)] = 0; missing = np.where(~np.isfinite(expression)); expression[missing] = np.take(medians, missing[1])
    X = design.to_numpy(float); inverse = np.linalg.pinv(X.T @ X); beta = inverse @ X.T @ expression; residual = expression - X @ beta; rank = np.linalg.matrix_rank(X); df_resid = max(len(train) - rank, 1); sigma2 = np.sum(residual ** 2, axis=0) / df_resid; ad_idx = list(design.columns).index("AD"); se = np.sqrt(np.maximum(sigma2 * inverse[ad_idx, ad_idx], 1e-300)); t = beta[ad_idx] / se; p = 2 * stats.t.sf(np.abs(t), df=df_resid); q = bh(p)
    table = pd.DataFrame({"Protein": eligible, "logFC": beta[ad_idx], "t": t, "P.Value": p, "adj.P.Val": q}).sort_values(["adj.P.Val", "P.Value"])
    selected = table.loc[table["adj.P.Val"] < 0.05, "Protein"].tolist(); fallback = False
    if len(selected) < 20: selected = table.head(min(100, len(table))).Protein.tolist(); fallback = True
    return table, selected, fallback

# Reslected matched protein model.
oof_reselect = pd.Series(np.nan, index=matched.index); metrics, rocs, panels, importance, dep_summaries = [], [], [], [], []
for fold, (train_idx, test_idx) in enumerate(outer, start=1):
    train, test = matched.iloc[train_idx], matched.iloc[test_idx]; dep, candidates, fallback_dep = training_dep(train, protein_candidates, fold); dep.to_csv(RESELECT / f"DEP_Fold{fold}.csv", index=False)
    selection = select_stable_panel(train, y.iloc[train_idx], candidates, SETTINGS, SETTINGS.random_state + fold * 100); panel = selection["panel"]; panels.append(panel)
    selection["frequency"].to_csv(RESELECT / f"frequency_fold_{fold}.csv", index=False); pd.DataFrame({"Protein": panel}).to_csv(RESELECT / f"panel_fold_{fold}.csv", index=False)
    estimator, params = fit_tuned_calibrated_svm(train[panel], y.iloc[train_idx], SETTINGS, SETTINGS.random_state + fold); probability = predict_probability(estimator, test[panel]); oof_reselect.iloc[test_idx] = probability
    metrics.append(metric_row(y.iloc[test_idx], probability, Fold=fold, Model="Matched reselected panel", N_train=len(train_idx), DEP_proteins=len(candidates), Consensus_panel=len(panel)))
    rocs.append(roc_table(y.iloc[test_idx], probability, Fold=fold, Model="Matched reselected panel")); importance.append(fold_permutation_importance(estimator, test[panel], y.iloc[test_idx], panel, SETTINGS, str(fold))); dep_summaries.append({"Fold": fold, "N_DEP_candidates": len(candidates), "DEP_fallback_top100": fallback_dep, "N_panel": len(panel), **params})
metrics_df = pd.DataFrame(metrics); metrics_df.to_csv(RESELECT / "metrics_nestedCV.csv", index=False); pd.concat(rocs, ignore_index=True).to_csv(RESELECT / "roc_curves.csv", index=False); pd.DataFrame(dep_summaries).to_csv(RESELECT / "best_params_outer.csv", index=False)
panel_frequency = aggregate_panel_frequency(panels); panel_frequency.rename("Selection_frequency").to_csv(RESELECT / "outer_panel_frequency.csv"); final_dir = RESELECT / "Final_Panels"; final_dir.mkdir(exist_ok=True); final_panel = panel_frequency[panel_frequency >= .8].index.tolist() or panel_frequency.head(SETTINGS.fallback_top_n).index.tolist(); pd.DataFrame({"Protein": final_panel}).to_csv(final_dir / "panel_final_freq08.csv", index=False); aggregate_importance(importance, panel_frequency).to_csv(final_dir / "permutation_importance_freq08.csv", index=False); pd.DataFrame({"SampleId": matched.SampleId, "y_true": y, "oof_probability": oof_reselect}).to_csv(RESELECT / "oof_predictions.csv", index=False)

# Prespecified primary panel on the full matched cohort.
for protein in PANEL:
    if protein not in matched.columns: raise KeyError(f"Fixed-panel protein not found: {protein}")
oof_fixed = pd.Series(np.nan, index=matched.index); fixed_metrics, fixed_rocs, fixed_importance, fixed_params = [], [], [], []
for fold, (train_idx, test_idx) in enumerate(outer, start=1):
    estimator, params = fit_tuned_calibrated_svm(matched.iloc[train_idx][PANEL], y.iloc[train_idx], SETTINGS, SETTINGS.random_state + 1000 + fold); probability = predict_probability(estimator, matched.iloc[test_idx][PANEL]); oof_fixed.iloc[test_idx] = probability
    fixed_metrics.append(metric_row(y.iloc[test_idx], probability, Fold=fold, Model="Primary fixed panel", N_train=len(train_idx), Consensus_panel=len(PANEL))); fixed_rocs.append(roc_table(y.iloc[test_idx], probability, Fold=fold, Model="Primary fixed panel")); fixed_importance.append(fold_permutation_importance(estimator, matched.iloc[test_idx][PANEL], y.iloc[test_idx], PANEL, SETTINGS, str(fold))); fixed_params.append({"Fold": fold, **params})
pd.DataFrame(fixed_metrics).to_csv(FIXED / "metrics_nestedCV.csv", index=False); pd.concat(fixed_rocs, ignore_index=True).to_csv(FIXED / "roc_curves.csv", index=False); pd.DataFrame(fixed_params).to_csv(FIXED / "best_params_outer_fixed_panel.csv", index=False); aggregate_importance(fixed_importance).to_csv(FIXED / "permutation_importance.csv", index=False); pd.DataFrame({"SampleId": matched.SampleId, "y_true": y, "oof_probability": oof_fixed}).to_csv(FIXED / "oof_predictions_fixed_panel.csv", index=False)

# Matched p-tau217 complete-case comparison; protein missingness is imputed inside CV.
ptau_col = resolve_column(matched, ("p.tau217", "p-tau217", "pTau217", "ptau217")); matched[ptau_col] = pd.to_numeric(matched[ptau_col], errors="coerce"); subset = matched[matched[ptau_col].notna()].reset_index(drop=True); y_sub = subset.SampleGroup.map({"CN": 0, "AD": 1}).astype(int); outer_sub = StratifiedKFold(n_splits=safe_n_splits(y_sub, 5), shuffle=True, random_state=SETTINGS.random_state)
model_features = {"Primary fixed panel": PANEL, "p-tau217": [ptau_col], "Combined": PANEL + [ptau_col]}; oof_models = {name: pd.Series(np.nan, index=subset.index) for name in model_features}; all_metrics, all_rocs, all_params = [], [], []
for fold, (train_idx, test_idx) in enumerate(outer_sub.split(subset, y_sub), start=1):
    for offset, (name, features) in enumerate(model_features.items(), start=1):
        estimator, params = fit_tuned_calibrated_svm(subset.iloc[train_idx][features], y_sub.iloc[train_idx], SETTINGS, SETTINGS.random_state + 2000 + fold * 10 + offset); probability = predict_probability(estimator, subset.iloc[test_idx][features]); oof_models[name].iloc[test_idx] = probability; all_metrics.append(metric_row(y_sub.iloc[test_idx], probability, Fold=fold, Model=name, N_train=len(train_idx))); all_rocs.append(roc_table(y_sub.iloc[test_idx], probability, Fold=fold, Model=name)); all_params.append({"Fold": fold, "Model": name, **params})
pd.DataFrame(all_metrics).to_csv(PTAU / "metrics_all_models.csv", index=False); pd.concat(all_rocs, ignore_index=True).to_csv(PTAU / "roc_curves.csv", index=False); pd.DataFrame(all_params).to_csv(PTAU / "best_params_outer.csv", index=False)
oof_ptau = pd.DataFrame({"SampleId": subset.SampleId, "y_true": y_sub, **{name: values for name, values in oof_models.items()}}); oof_ptau.to_csv(PTAU / "oof_predictions.csv", index=False)
comparisons = []
for label, a, b in (("Panel_vs_pTau", "Primary fixed panel", "p-tau217"), ("Panel_vs_Combined", "Primary fixed panel", "Combined"), ("pTau_vs_Combined", "p-tau217", "Combined")): comparisons.append({"Comparison": label, **delong_pairwise(y_sub, oof_models[a], oof_models[b])})
pd.DataFrame(comparisons).to_csv(PTAU / "delong_results.csv", index=False)
summary = pd.DataFrame(all_metrics).groupby("Model").agg(Mean_AUC=("AUC", "mean"), SD_AUC=("AUC", "std"), N_effective=("N_test", "sum")).reset_index(); summary.to_csv(PTAU / "model_comparison_summary.csv", index=False)
pd.DataFrame({"matched_total": [len(matched)], "matched_CN": [int((y == 0).sum())], "matched_AD": [int((y == 1).sum())], "ptau_complete_total": [len(subset)], "ptau_complete_CN": [int((y_sub == 0).sum())], "ptau_complete_AD": [int((y_sub == 1).sum())], "protein_complete_case_required": [False]}).to_csv(PTAU / "cohort_summary.csv", index=False)
pd.DataFrame([
    {"Analysis": "matched_reselected", "Matching_fixed_before_CV": True, "DEP_training_only": True, "Feature_filter_training_only": True, "Preprocessing_inside_inner_CV": True, "OOF_foldwise_rescaling": False, "Importance_outer_test_only": True},
    {"Analysis": "matched_fixed_panel", "Matching_fixed_before_CV": True, "Panel_prespecified": True, "Preprocessing_inside_inner_CV": True, "OOF_foldwise_rescaling": False, "Importance_outer_test_only": True},
    {"Analysis": "matched_ptau", "Matching_fixed_before_CV": True, "Cohort_defined_by_ptau_availability": True, "Protein_complete_case_required": False, "Preprocessing_inside_inner_CV": True, "OOF_foldwise_rescaling": False},
]).to_csv(ROOT / "leakage_audit.csv", index=False)
settings_table(SETTINGS).to_csv(ROOT / "strict_pipeline_parameters.csv", index=False)
pd.DataFrame({"Analysis": ["Matched reselected", "Matched fixed panel"], "Mean_AUC": [metrics_df.AUC.mean(), pd.DataFrame(fixed_metrics).AUC.mean()], "SD_AUC": [metrics_df.AUC.std(), pd.DataFrame(fixed_metrics).AUC.std()], "N": [len(matched), len(matched)]}).to_csv(ROOT / "matched_model_summary.csv", index=False)
print(pd.read_csv(ROOT / "matched_model_summary.csv")); print(summary)

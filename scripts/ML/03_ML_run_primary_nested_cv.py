###############################################################################
# ReDLat plasma proteomics — strict machine-learning workflow
# 03. Run primary nested classification
# Requires: Scripts 01–02 outputs
# Produces: strict nested-CV metrics, recurrent panel and private OOF probabilities
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
from redlat_ml.config import load_config, require_files
from redlat_ml.strict_cv import (
    StrictCVSettings, aggregate_importance, aggregate_panel_frequency,
    assert_fold_integrity, fit_tuned_calibrated_svm, fold_permutation_importance,
    load_candidate_genes, metric_row, predict_probability, read_ids, roc_table,
    select_stable_panel, settings_table,
)

CONFIG = load_config(__file__)
SETTINGS = StrictCVSettings.from_env()
RESULTS_DIR = CONFIG.private_root / "primary" / "02_nested_cv"
FOLDS_DIR = CONFIG.private_root / "primary" / "00_folds"
DEP_DIR = CONFIG.private_root / "primary" / "01_dep_folds"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
require_files([(CONFIG.master_file, "ML master matrix")])

df = pd.read_csv(CONFIG.master_file, low_memory=False)
if not {"SampleId", "SampleGroup"}.issubset(df.columns):
    raise KeyError("The master matrix must contain SampleId and SampleGroup.")
df["SampleId"] = df["SampleId"].astype(str).str.strip()
df = df[df["SampleGroup"].isin(["CN", "AD"])].drop_duplicates("SampleId").set_index("SampleId")
y = df["SampleGroup"].map({"CN": 0, "AD": 1}).astype(int)

oof = pd.Series(np.nan, index=df.index, dtype=float)
metrics, best_params, selection_summary = [], [], []
roc_tables, outer_panels, importance_tables, leakage_rows = [], [], [], []

for fold in range(1, 6):
    train_path = FOLDS_DIR / f"fold_{fold}_train_ids.csv"
    test_path = FOLDS_DIR / f"fold_{fold}_test_ids.csv"
    candidate_path = DEP_DIR / f"fold_{fold}" / "candidate_gene_symbols.csv"
    require_files([(train_path, "training IDs"), (test_path, "test IDs"), (candidate_path, "fold-specific DEP candidates")])
    train_ids, test_ids = read_ids(train_path), read_ids(test_path)
    assert_fold_integrity(train_ids, test_ids, df.index, f"Fold {fold}")
    candidates = load_candidate_genes(candidate_path)
    selection = select_stable_panel(df.loc[train_ids], y.loc[train_ids], candidates, SETTINGS, SETTINGS.random_state + fold * 100)
    panel = selection["panel"]
    selection["frequency"].to_csv(RESULTS_DIR / f"frequency_fold_{fold}.csv", index=False)
    selection["selector_audit"].assign(Fold=fold).to_csv(RESULTS_DIR / f"selector_audit_fold_{fold}.csv", index=False)
    pd.DataFrame({"Protein": panel}).to_csv(RESULTS_DIR / f"panel_fold_{fold}.csv", index=False)
    estimator, params = fit_tuned_calibrated_svm(df.loc[train_ids, panel], y.loc[train_ids], SETTINGS, SETTINGS.random_state + fold)
    probability = predict_probability(estimator, df.loc[test_ids, panel])
    oof.loc[test_ids] = probability
    metrics.append(metric_row(y.loc[test_ids], probability, Fold=fold, Model="Proteins", N_train=len(train_ids), DEP_proteins=len(candidates), Consensus_panel=len(panel)))
    best_params.append({"Fold": fold, **params, "DEP_proteins": len(candidates), "Consensus_panel": len(panel), "AUC": metrics[-1]["AUC"]})
    selection_summary.append({"Fold": fold, "N_features": len(panel), "Fallback_used": selection["fallback_used"], "N_training_eligible": selection["n_outer_eligible"]})
    roc_tables.append(roc_table(y.loc[test_ids], probability, Fold=fold, Model="Proteins"))
    importance_tables.append(fold_permutation_importance(estimator, df.loc[test_ids, panel], y.loc[test_ids], panel, SETTINGS, str(fold)))
    outer_panels.append(panel)
    leakage_rows.append({
        "Fold": fold, "N_train": len(train_ids), "N_test": len(test_ids), "Train_test_overlap": 0,
        "DEP_candidates_training_only": True, "Feature_filter_training_only": True,
        "Imputation_inside_inner_CV": True, "Scaling_inside_inner_CV": True,
        "Elastic_net_tuning_inside_training": True, "SVM_tuning_inside_training": True,
        "OOF_foldwise_rescaling": False, "Importance_outer_test_only": True,
    })

if oof.isna().any():
    raise RuntimeError(f"Missing OOF probabilities for {int(oof.isna().sum())} participants.")
metrics_df = pd.DataFrame(metrics)
metrics_df.to_csv(RESULTS_DIR / "metrics_nestedCV.csv", index=False)
pd.DataFrame(best_params).to_csv(RESULTS_DIR / "best_params_outer.csv", index=False)
pd.DataFrame(selection_summary).to_csv(RESULTS_DIR / "feature_counts.csv", index=False)
pd.DataFrame(leakage_rows).to_csv(RESULTS_DIR / "leakage_audit.csv", index=False)
settings_table(SETTINGS).to_csv(RESULTS_DIR / "pipeline_parameters.csv", index=False)
metrics_df.mean(numeric_only=True).to_csv(RESULTS_DIR / "metrics_mean.csv")

panel_frequency = aggregate_panel_frequency(outer_panels)
panel_frequency.rename("Selection_frequency").to_csv(RESULTS_DIR / "outer_panel_frequency.csv")
final_dir = RESULTS_DIR / "Final_Panels"
final_dir.mkdir(exist_ok=True)
for tag, threshold in (("freq06", 0.6), ("freq08", 0.8), ("freq10", 1.0)):
    panel = panel_frequency[panel_frequency >= threshold].index.tolist()
    if not panel:
        panel = panel_frequency.head(SETTINGS.fallback_top_n).index.tolist()
    pd.DataFrame({"Protein": panel}).to_csv(final_dir / f"panel_final_{tag}.csv", index=False)

importance = aggregate_importance(importance_tables, panel_frequency)
importance.to_csv(final_dir / "permutation_importance_freq08.csv", index=False)
pd.concat(roc_tables, ignore_index=True).to_csv(RESULTS_DIR / "roc_curves.csv", index=False)
pd.DataFrame({"SampleId": df.index, "y_true": y, "oof": oof}).to_csv(RESULTS_DIR / "oof_predictions.csv", index=False)
pd.DataFrame({
    "Model": ["SomaScan"], "Mean_CN": [oof[y == 0].mean()], "Mean_AD": [oof[y == 1].mean()],
    "Score_type": ["calibrated OOF probability"],
}).to_csv(RESULTS_DIR / "score_summary.csv", index=False)
pd.DataFrame({
    "Stage": ["Eligible clinical cohort", "OOF predictions"],
    "N": [len(df), int(oof.notna().sum())],
    "CN": [int((y == 0).sum()), int((y[oof.notna()] == 0).sum())],
    "AD": [int((y == 1).sum()), int((y[oof.notna()] == 1).sum())],
}).to_csv(RESULTS_DIR / "cohort_flow.csv", index=False)
print(metrics_df.groupby("Model")["AUC"].agg(["mean", "std"]))

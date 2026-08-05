###############################################################################
# ReDLat plasma proteomics — strict machine-learning workflow
# 06. Run p-tau217 nested comparison
# Requires: Scripts 04–05 outputs
# Produces: strict protein, p-tau217 and combined OOF comparisons
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
    assert_fold_integrity, delong_pairwise, fit_tuned_calibrated_svm,
    fold_permutation_importance, load_candidate_genes, metric_row,
    predict_probability, read_ids, resolve_column, roc_table,
    select_stable_panel, settings_table,
)

CONFIG = load_config(__file__)
SETTINGS = StrictCVSettings.from_env()
RESULTS_DIR = CONFIG.private_root / "ptau" / "02_nested_cv"
FOLDS_DIR = CONFIG.private_root / "ptau" / "00_folds"
DEP_DIR = CONFIG.private_root / "ptau" / "01_dep_folds"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
require_files([(CONFIG.master_file, "ML master matrix")])

df = pd.read_csv(CONFIG.master_file, low_memory=False)
df["SampleId"] = df["SampleId"].astype(str).str.strip()
df = df[df["SampleGroup"].isin(["CN", "AD"])].drop_duplicates("SampleId").set_index("SampleId")
ptau_col = resolve_column(df, ("p.tau217", "p-tau217", "pTau217", "ptau217"))
df[ptau_col] = pd.to_numeric(df[ptau_col], errors="coerce")
y = df["SampleGroup"].map({"CN": 0, "AD": 1}).astype(int)
models = ["Proteins", "pTau217", "Proteins + p-tau217"]
oof = {model: pd.Series(np.nan, index=df.index, dtype=float) for model in models}
metrics = {model: [] for model in models}
best = {model: [] for model in models}
roc_tables, outer_panels, importance_tables, leakage_rows = [], [], [], []
all_fold_ids = set()

for fold in range(1, 6):
    train_path = FOLDS_DIR / f"fold_{fold}_train_ids.csv"
    test_path = FOLDS_DIR / f"fold_{fold}_test_ids.csv"
    candidate_path = DEP_DIR / f"fold_{fold}" / "candidate_gene_symbols.csv"
    require_files([(train_path, "training IDs"), (test_path, "test IDs"), (candidate_path, "fold-specific DEP candidates")])
    train_ids, test_ids = read_ids(train_path), read_ids(test_path)
    assert_fold_integrity(train_ids, test_ids, df.index, f"Fold {fold}")
    all_fold_ids.update(train_ids); all_fold_ids.update(test_ids)
    missing_ptau = df.loc[list(set(train_ids) | set(test_ids)), ptau_col].isna()
    if missing_ptau.any():
        raise RuntimeError(f"Fold {fold}: p-tau217 is missing for {int(missing_ptau.sum())} prespecified fold participants.")
    candidates = load_candidate_genes(candidate_path)
    selection = select_stable_panel(df.loc[train_ids], y.loc[train_ids], candidates, SETTINGS, SETTINGS.random_state + fold * 100)
    panel = selection["panel"]
    selection["frequency"].to_csv(RESULTS_DIR / f"frequency_fold_{fold}.csv", index=False)
    pd.DataFrame({"Protein": panel}).to_csv(RESULTS_DIR / f"panel_fold_{fold}.csv", index=False)
    outer_panels.append(panel)

    feature_sets = {
        "Proteins": panel,
        "pTau217": [ptau_col],
        "Proteins + p-tau217": panel + [ptau_col],
    }
    for offset, (model, features) in enumerate(feature_sets.items(), start=1):
        estimator, params = fit_tuned_calibrated_svm(df.loc[train_ids, features], y.loc[train_ids], SETTINGS, SETTINGS.random_state + fold * 10 + offset)
        probability = predict_probability(estimator, df.loc[test_ids, features])
        oof[model].loc[test_ids] = probability
        row = metric_row(y.loc[test_ids], probability, Fold=fold, Model=model, N_train=len(train_ids), DEP_proteins=len(candidates), Consensus_panel=len(panel))
        metrics[model].append(row)
        best[model].append({"Fold": fold, "Model": model, **params, "AUC": row["AUC"], "Consensus_panel": len(panel)})
        roc_tables.append(roc_table(y.loc[test_ids], probability, Fold=fold, Model=model))
        if model == "Proteins + p-tau217":
            importance_tables.append(fold_permutation_importance(estimator, df.loc[test_ids, features], y.loc[test_ids], features, SETTINGS, str(fold)))
    leakage_rows.append({
        "Fold": fold, "N_train": len(train_ids), "N_test": len(test_ids), "Train_test_overlap": 0,
        "Cohort_defined_before_protein_missingness": True, "Protein_missingness_imputed_inside_CV": True,
        "Biomarker_complete_by_prespecified_fold": True, "DEP_candidates_training_only": True,
        "Preprocessing_inside_inner_CV": True, "OOF_foldwise_rescaling": False,
    })

fold_ids = sorted(all_fold_ids)
for model in models:
    if oof[model].loc[fold_ids].isna().any():
        raise RuntimeError(f"{model}: incomplete OOF probability coverage.")

name_map = {"Proteins": "soma", "pTau217": "ptau", "Proteins + p-tau217": "combo"}
summary_rows = []
for model in models:
    table = pd.DataFrame(metrics[model])
    suffix = name_map[model]
    table.to_csv(RESULTS_DIR / f"metrics_{suffix}.csv", index=False)
    pd.DataFrame(best[model]).to_csv(RESULTS_DIR / f"best_params_outer_{suffix}.csv", index=False)
    table.mean(numeric_only=True).to_csv(RESULTS_DIR / f"metrics_mean_{suffix}.csv")
    summary_rows.append({
        "Model": "SomaScan" if model == "Proteins" else ("pTau217" if model == "pTau217" else "Combo"),
        "Mean_AUC": table["AUC"].mean(), "SD_AUC": table["AUC"].std(),
        "Mean_Accuracy": table["Accuracy"].mean(), "Mean_Sensitivity": table["Sensitivity"].mean(),
        "Mean_Specificity": table["Specificity"].mean(), "N_effective": len(fold_ids),
    })
pd.DataFrame(summary_rows).to_csv(RESULTS_DIR / "model_comparison_summary.csv", index=False)

panel_frequency = aggregate_panel_frequency(outer_panels)
panel_frequency.rename("Selection_frequency").to_csv(RESULTS_DIR / "outer_panel_frequency.csv")
final_dir = RESULTS_DIR / "Final_Panels"; final_dir.mkdir(exist_ok=True)
for tag, threshold in (("freq06", 0.6), ("freq08", 0.8), ("freq10", 1.0)):
    panel = panel_frequency[panel_frequency >= threshold].index.tolist() or panel_frequency.head(SETTINGS.fallback_top_n).index.tolist()
    pd.DataFrame({"Protein": panel}).to_csv(final_dir / f"panel_final_{tag}.csv", index=False)
aggregate_importance(importance_tables, panel_frequency).to_csv(final_dir / "permutation_importance_freq08.csv", index=False)
pd.concat(roc_tables, ignore_index=True).to_csv(RESULTS_DIR / "roc_curves.csv", index=False)

out = pd.DataFrame({"SampleId": fold_ids, "y_true": y.loc[fold_ids].values})
for model in models:
    out[name_map[model] if model != "Proteins + p-tau217" else "oof_combo"] = oof[model].loc[fold_ids].values
out = out.rename(columns={"soma": "oof_soma", "ptau": "oof_ptau"})
out.to_csv(RESULTS_DIR / "oof_predictions.csv", index=False)

comparisons = []
for label, a, b in (("Soma_vs_pTau", "Proteins", "pTau217"), ("Soma_vs_Combo", "Proteins", "Proteins + p-tau217"), ("pTau_vs_Combo", "pTau217", "Proteins + p-tau217")):
    comparisons.append({"Comparison": label, **delong_pairwise(y.loc[fold_ids], oof[a].loc[fold_ids], oof[b].loc[fold_ids])})
pd.DataFrame(comparisons).to_csv(RESULTS_DIR / "delong_results.csv", index=False)
pd.DataFrame(leakage_rows).to_csv(RESULTS_DIR / "leakage_audit.csv", index=False)
settings_table(SETTINGS).to_csv(RESULTS_DIR / "pipeline_parameters.csv", index=False)
pd.DataFrame({
    "Stage": ["p-tau217 available and fold-eligible", "OOF predictions"],
    "N": [len(fold_ids), len(fold_ids)],
    "CN": [int((y.loc[fold_ids] == 0).sum())] * 2,
    "AD": [int((y.loc[fold_ids] == 1).sum())] * 2,
    "Protein_complete_case_required": [False, False],
}).to_csv(RESULTS_DIR / "selected_subjects.csv", index=False)
print(pd.DataFrame(summary_rows)[["Model", "Mean_AUC", "SD_AUC", "N_effective"]])

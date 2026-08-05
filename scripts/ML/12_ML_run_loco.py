###############################################################################
# ReDLat plasma proteomics — strict machine-learning workflow
# 12. Run country-held-out classification
# Requires: Scripts 10–11 outputs
# Produces: strict country-held-out metrics and private probabilities
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

CONFIG = load_config(__file__); SETTINGS = StrictCVSettings.from_env()
RESULTS_DIR = CONFIG.private_root / "loco" / "02_nested_cv"; RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FOLDS_DIR = CONFIG.private_root / "loco" / "00_folds"; DEP_DIR = CONFIG.private_root / "loco" / "01_dep_folds"
require_files([(CONFIG.master_file, "ML master matrix")])
df = pd.read_csv(CONFIG.master_file, low_memory=False); df["SampleId"] = df["SampleId"].astype(str).str.strip()
df = df[df["SampleGroup"].isin(["CN", "AD"])].drop_duplicates("SampleId").set_index("SampleId")
y = df["SampleGroup"].map({"CN": 0, "AD": 1}).astype(int)
fold_files = sorted(FOLDS_DIR.glob("*_test_ids.csv"))
if not fold_files: raise FileNotFoundError(f"No LOCO fold files found under {FOLDS_DIR}")
metrics, best_params, selection_summary, roc_tables, outer_panels, importance_tables, leakage_rows = [], [], [], [], [], [], []
oof = pd.Series(np.nan, index=df.index, dtype=float)
for test_file in fold_files:
    country = test_file.name.removesuffix("_test_ids.csv")
    train_file = FOLDS_DIR / f"{country}_train_ids.csv"; candidate_file = DEP_DIR / country / "candidate_gene_symbols.csv"
    require_files([(train_file, "LOCO training IDs"), (candidate_file, "LOCO DEP candidates")])
    train_ids, test_ids = read_ids(train_file), read_ids(test_file); assert_fold_integrity(train_ids, test_ids, df.index, country)
    candidates = load_candidate_genes(candidate_file)
    selection = select_stable_panel(df.loc[train_ids], y.loc[train_ids], candidates, SETTINGS, SETTINGS.random_state + len(metrics) * 100)
    panel = selection["panel"]; outer_panels.append(panel)
    selection["frequency"].to_csv(RESULTS_DIR / f"frequency_{country}.csv", index=False); pd.DataFrame({"Protein": panel}).to_csv(RESULTS_DIR / f"panel_{country}.csv", index=False)
    estimator, params = fit_tuned_calibrated_svm(df.loc[train_ids, panel], y.loc[train_ids], SETTINGS, SETTINGS.random_state + len(metrics) + 1)
    probability = predict_probability(estimator, df.loc[test_ids, panel]); oof.loc[test_ids] = probability
    row = metric_row(y.loc[test_ids], probability, Country=country, Fold=country, Model="Proteins", N_train=len(train_ids), DEP_proteins=len(candidates), Consensus_panel=len(panel))
    metrics.append(row); best_params.append({"Country": country, **params, "AUC": row["AUC"], "Consensus_panel": len(panel)})
    selection_summary.append({"Country": country, "N_features": len(panel), "Fallback_used": selection["fallback_used"]})
    roc_tables.append(roc_table(y.loc[test_ids], probability, Country=country, Fold=country, Model="Proteins"))
    importance_tables.append(fold_permutation_importance(estimator, df.loc[test_ids, panel], y.loc[test_ids], panel, SETTINGS, country))
    leakage_rows.append({"Country": country, "N_train": len(train_ids), "N_test": len(test_ids), "Train_test_overlap": 0, "Held_out_country_absent_from_training": True, "DEP_candidates_training_only": True, "Preprocessing_inside_inner_CV": True, "OOF_foldwise_rescaling": False})
metrics_df = pd.DataFrame(metrics); metrics_df.to_csv(RESULTS_DIR / "metrics_nestedCV.csv", index=False); metrics_df.mean(numeric_only=True).to_csv(RESULTS_DIR / "metrics_mean.csv")
pd.DataFrame(best_params).to_csv(RESULTS_DIR / "best_params_outer.csv", index=False); pd.DataFrame(selection_summary).to_csv(RESULTS_DIR / "feature_counts.csv", index=False); pd.DataFrame(leakage_rows).to_csv(RESULTS_DIR / "leakage_audit.csv", index=False)
panel_frequency = aggregate_panel_frequency(outer_panels); panel_frequency.rename("Selection_frequency").to_csv(RESULTS_DIR / "outer_panel_frequency.csv")
final_dir = RESULTS_DIR / "Final_Panels"; final_dir.mkdir(exist_ok=True)
for tag, threshold in (("freq06", .6), ("freq08", .8), ("freq10", 1.0)):
    panel = panel_frequency[panel_frequency >= threshold].index.tolist() or panel_frequency.head(SETTINGS.fallback_top_n).index.tolist(); pd.DataFrame({"Protein": panel}).to_csv(final_dir / f"panel_final_{tag}.csv", index=False)
aggregate_importance(importance_tables, panel_frequency).to_csv(final_dir / "permutation_importance_freq08.csv", index=False)
pd.concat(roc_tables, ignore_index=True).to_csv(RESULTS_DIR / "roc_curves.csv", index=False)
covered = oof.notna(); pd.DataFrame({"SampleId": df.index[covered], "y_true": y[covered], "oof": oof[covered]}).to_csv(RESULTS_DIR / "oof_predictions.csv", index=False)
pd.DataFrame({"Model": ["SomaScan"], "Mean_CN": [oof[(y == 0) & covered].mean()], "Mean_AD": [oof[(y == 1) & covered].mean()], "Score_type": ["calibrated OOF probability"]}).to_csv(RESULTS_DIR / "score_summary.csv", index=False)
settings_table(SETTINGS).to_csv(RESULTS_DIR / "pipeline_parameters.csv", index=False)
print(metrics_df[["Country", "AUC", "N_test"]])

###############################################################################
# ReDLat plasma proteomics — machine-learning workflow
# 23. Generate ML Source Data
# Requires: validated aggregate outputs from Scripts 03, 06, 09, 13, 14 and 17
# Produces: Source Data for Figure 4 and Supplementary Data 7
# Data policy: participant-level inputs and predictions remain local.
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
from redlat_ml.config import load_config, require_files
CONFIG = load_config(__file__)

from pathlib import Path
import pandas as pd

from redlat_ml.config import assert_public_table


def read_csv(path: Path, label: str, required: bool = True) -> pd.DataFrame:
    if not path.exists():
        if required:
            raise FileNotFoundError(f"{label} not found: {path}")
        return pd.DataFrame()
    return pd.read_csv(path)


def first_existing(paths: list[Path]) -> Path | None:
    return next((path for path in paths if path.exists()), None)


def safe_sheet(writer, name: str, title: str, data: pd.DataFrame, note: str | None = None) -> None:
    if data is None or data.empty:
        return
    assert_public_table(data, name)
    start = 2 if note is None else 3
    pd.DataFrame([[title]]).to_excel(writer, sheet_name=name, index=False, header=False)
    if note is not None:
        pd.DataFrame([[note]]).to_excel(writer, sheet_name=name, index=False, header=False, startrow=1)
    data.to_excel(writer, sheet_name=name, index=False, startrow=start)


private = CONFIG.private_root
publication = CONFIG.publication_root
source_dir = publication / "source_data"
supp_dir = publication / "supplementary_data"
source_dir.mkdir(parents=True, exist_ok=True)
supp_dir.mkdir(parents=True, exist_ok=True)

primary = private / "primary" / "02_nested_cv"
ptau = private / "ptau" / "02_nested_cv"
apoe = private / "apoe" / "02_nested_cv"
fixed = private / "fixed_panel_ptau"
reg = private / "clinical_regression"
loco = private / "loco" / "02_nested_cv"
matched = private / "matched"
loco_meta = publication / "loco_meta_analysis" / "source_data"

panel_file = first_existing([
    primary / "Final_Panels" / "panel_final_freq08.csv",
    primary / "Final_Panels" / "panel_final_08.csv",
])
importance_file = first_existing([
    primary / "Final_Panels" / "permutation_importance_freq08.csv",
    primary / "Final_Panels" / "permutation_importance_08.csv",
])
loco_importance_file = first_existing([
    loco / "Final_Panels" / "permutation_importance_freq08.csv",
    loco / "Final_Panels" / "permutation_importance_08.csv",
])

primary_metrics = read_csv(primary / "metrics_nestedCV.csv", "Primary metrics")
primary_roc = read_csv(primary / "roc_curves.csv", "Primary ROC")
primary_frequency = read_csv(primary / "outer_panel_frequency.csv", "Primary panel frequency")
primary_panel = read_csv(panel_file, "Primary final panel") if panel_file else pd.DataFrame()
primary_importance = read_csv(importance_file, "Primary permutation importance") if importance_file else pd.DataFrame()
ptau_summary = read_csv(ptau / "model_comparison_summary.csv", "p-tau217 summary")
ptau_soma = read_csv(ptau / "metrics_soma.csv", "p-tau217 protein metrics")
ptau_ptau = read_csv(ptau / "metrics_ptau.csv", "p-tau217 metrics")
ptau_combo = read_csv(ptau / "metrics_combo.csv", "combined metrics")
ptau_fold = pd.concat([
    ptau_soma.assign(Model="Proteins"),
    ptau_ptau.assign(Model="p-tau217"),
    ptau_combo.assign(Model="Combined"),
], ignore_index=True)
ptau_roc = read_csv(ptau / "roc_curves.csv", "p-tau217 ROC")
ptau_delong = read_csv(ptau / "delong_results.csv", "p-tau217 DeLong")
apoe_summary = read_csv(apoe / "model_comparison_summary.csv", "APOE summary")
apoe_soma = read_csv(apoe / "metrics_soma.csv", "APOE protein metrics")
apoe_combo = read_csv(apoe / "metrics_combo.csv", "APOE combined metrics")
apoe_fold = pd.concat([
    apoe_soma.assign(Model="Proteins"),
    apoe_combo.assign(Model="Proteins + APOE ε4"),
], ignore_index=True)
apoe_roc = read_csv(apoe / "roc_curves.csv", "APOE ROC")
apoe_delong = read_csv(apoe / "delong_results.csv", "APOE DeLong")
apoe_importance_path = first_existing([
    apoe / "Final_Panels" / "permutation_importance_freq08.csv",
    apoe / "Final_Panels" / "permutation_importance_08.csv",
])
apoe_importance = read_csv(apoe_importance_path, "APOE importance") if apoe_importance_path else pd.DataFrame()
fixed_auc = read_csv(fixed / "auc_summary.csv", "Fixed-panel AUC", required=False)
fixed_stats = read_csv(fixed / "model_statistics.csv", "Fixed-panel statistics", required=False)
fixed_lrt = read_csv(fixed / "likelihood_ratio_tests.csv", "Fixed-panel LRT", required=False)
fixed_delong = read_csv(fixed / "delong_results.csv", "Fixed-panel DeLong", required=False)
fixed_coefs = pd.concat([
    read_csv(path, path.stem, required=False).assign(Source_file=path.name)
    for path in sorted(fixed.glob("coefficients_*.csv"))
], ignore_index=True) if list(fixed.glob("coefficients_*.csv")) else pd.DataFrame()
clinical_files = sorted(reg.glob("Regression_CN_AD_*.csv"))
clinical_coefs = pd.concat([
    read_csv(path, path.stem, required=False).assign(Outcome=path.stem.replace("Regression_CN_AD_", ""))
    for path in clinical_files
], ignore_index=True) if clinical_files else pd.DataFrame()
clinical_summary = read_csv(reg / "model_summary.csv", "Clinical summary", required=False)
loco_metrics = read_csv(loco / "metrics_nestedCV.csv", "LOCO metrics")
loco_roc = read_csv(loco / "roc_curves.csv", "LOCO ROC")
loco_importance = read_csv(loco_importance_file, "LOCO importance") if loco_importance_file else pd.DataFrame()
country_auc = read_csv(loco_meta / "LOCO_country_AUC_estimates.csv", "Country AUC", required=False)
meta_summary = read_csv(loco_meta / "LOCO_AUC_meta_analysis_summary.csv", "LOCO meta-analysis", required=False)

fig4_path = source_dir / "Source_Data_Fig_4.xlsx"
with pd.ExcelWriter(fig4_path, engine="openpyxl") as writer:
    safe_sheet(writer, "Fig4a_Primary_ROC", "Fig. 4a primary protein-only ROC curves", primary_roc)
    safe_sheet(writer, "Fig4a_Primary_Metrics", "Fig. 4a primary protein-only fold metrics", primary_metrics)
    safe_sheet(writer, "Fig4a_Final_Panel", "Fig. 4a final recurrent protein panel", primary_panel)
    safe_sheet(writer, "Fig4a_Importance", "Fig. 4a final-panel permutation importance", primary_importance)
    safe_sheet(writer, "Fig4b_pTau_ROC", "Fig. 4b p-tau217 comparison ROC curves", ptau_roc)
    safe_sheet(writer, "Fig4b_pTau_Summary", "Fig. 4b p-tau217 model comparison summary", ptau_summary)
    safe_sheet(writer, "Fig4b_DeLong", "Fig. 4b DeLong comparisons", ptau_delong)
    safe_sheet(writer, "Fig4c_Model_Summary", "Fig. 4c clinical-alignment model statistics", clinical_summary)
    safe_sheet(writer, "Fig4c_Coefficients", "Fig. 4c clinical-alignment protein coefficients", clinical_coefs)
    safe_sheet(writer, "Fig4d_LOCO_ROC", "Fig. 4d leave-one-country-out ROC curves", loco_roc)
    safe_sheet(writer, "Fig4d_LOCO_Importance", "Fig. 4d LOCO permutation importance", loco_importance)
    safe_sheet(writer, "Fig4e_Country_AUC", "Fig. 4e country-specific LOCO AUC estimates", country_auc)
    safe_sheet(writer, "Fig4e_Meta_Analysis", "Fig. 4e LOCO AUC meta-analysis", meta_summary)

# Matched files are generated by Script 16. Only aggregate tables are included.
matched_candidates = {
    "Match_Balance_191": matched / "00_matching" / "balance_summary.csv",
    "Match_Reselect_Summary": matched / "02_nested_cv" / "model_summary.csv",
    "Match_Reselect_Metrics": matched / "02_nested_cv" / "metrics_nestedCV.csv",
    "Match_Reselect_ROC": matched / "02_nested_cv" / "roc_curves.csv",
    "Match_Reselect_Frequency": matched / "02_nested_cv" / "outer_panel_frequency.csv",
    "Match_Reselect_Panel": matched / "02_nested_cv" / "Final_Panels" / "panel_final_freq08.csv",
    "Match_Reselect_Import": matched / "02_nested_cv" / "Final_Panels" / "permutation_importance_freq08.csv",
    "Match_Fixed_Summary": matched / "03_primary_fixed_panel_cv" / "model_summary.csv",
    "Match_Fixed_Metrics": matched / "03_primary_fixed_panel_cv" / "metrics_nestedCV.csv",
    "Match_Fixed_ROC": matched / "03_primary_fixed_panel_cv" / "roc_curves.csv",
    "Match_Fixed_Import": matched / "03_primary_fixed_panel_cv" / "permutation_importance.csv",
    "Match_Panel_Overlap": matched / "matched_panel_overlap.csv",
    "Match_pTau_Cohort": matched / "04_primary_panel_ptau" / "cohort_summary.csv",
    "Match_pTau_Summary": matched / "04_primary_panel_ptau" / "model_comparison_summary.csv",
    "Match_pTau_Metrics": matched / "04_primary_panel_ptau" / "metrics_all_models.csv",
    "Match_pTau_ROC": matched / "04_primary_panel_ptau" / "roc_curves.csv",
    "Match_pTau_DeLong": matched / "04_primary_panel_ptau" / "delong_results.csv",
}

supp_path = supp_dir / "Supplementary_Data_7.xlsx"
with pd.ExcelWriter(supp_path, engine="openpyxl") as writer:
    tables = {
        "Primary_Fold_Metrics": primary_metrics,
        "Primary_ROC": primary_roc,
        "Primary_Panel_Frequency": primary_frequency,
        "Primary_Final_Panel": primary_panel,
        "Primary_Permutation": primary_importance,
        "pTau_Model_Summary": ptau_summary,
        "pTau_Fold_Metrics": ptau_fold,
        "pTau_ROC": ptau_roc,
        "pTau_DeLong": ptau_delong,
        "APOE_Model_Summary": apoe_summary,
        "APOE_Fold_Metrics": apoe_fold,
        "APOE_ROC": apoe_roc,
        "APOE_DeLong": apoe_delong,
        "APOE_Permutation": apoe_importance,
        "Final_Logistic_AUC": fixed_auc,
        "Final_Logistic_Stats": fixed_stats,
        "Final_Logistic_LRT": fixed_lrt,
        "Final_Logistic_DeLong": fixed_delong,
        "Final_Logistic_Coefs": fixed_coefs,
        "Clinical_Model_Summary": clinical_summary,
        "Clinical_Coefficients": clinical_coefs,
        "LOCO_Fold_Metrics": loco_metrics,
        "LOCO_ROC": loco_roc,
        "LOCO_Permutation": loco_importance,
        "LOCO_Country_AUC": country_auc,
        "LOCO_Meta_Analysis": meta_summary,
    }
    for sheet, path in matched_candidates.items():
        tables[sheet] = read_csv(path, sheet, required=False)
    for sheet, data in tables.items():
        safe_sheet(writer, sheet[:31], sheet.replace("_", " "), data)

print(f"Generated: {fig4_path}")
print(f"Generated: {supp_path}")

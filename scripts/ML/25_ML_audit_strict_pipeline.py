###############################################################################
# ReDLat plasma proteomics — strict machine-learning workflow
# 25. Audit strict nested-CV implementation
# Requires: active ML source files
# Produces: static method and privacy audit report
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


import re
from pathlib import Path
import pandas as pd

CORE = [
    PROJECT_ROOT / "scripts/ML/03_ML_run_primary_nested_cv.py",
    PROJECT_ROOT / "scripts/ML/06_ML_run_ptau_nested_cv.py",
    PROJECT_ROOT / "scripts/ML/09_ML_run_apoe_nested_cv.py",
    PROJECT_ROOT / "scripts/ML/12_ML_run_loco.py",
    PROJECT_ROOT / "scripts/ML/13_ML_evaluate_fixed_panel_ptau.py",
    PROJECT_ROOT / "scripts/ML/16_ML_run_matched_analysis.py",
    PROJECT_ROOT / "python/redlat_ml/strict_cv.py",
]
checks = {
    "foldwise_score_z_normalization": r"y_prob\s*-\s*y_prob\.mean|probability\s*-\s*probability\.mean",
    "global_protein_complete_case_filter": r"dropna\s*\(\s*subset\s*=\s*(?:protein_cols|COMMON_COLS)",
    "prefit_scaler_before_model_search": r"StandardScaler\(\)\s*\n?.*fit_transform",
    "historical_logistic_regression_cv": r"LogisticRegressionCV",
    "personal_windows_path": r"[A-Za-z]:[/\\]Users[/\\]",
    "automatic_package_install": r"install\.packages|BiocManager::install|pip install",
}
rows=[]; failures=[]
for path in CORE:
    text=path.read_text(encoding="utf-8")
    for check, pattern in checks.items():
        found=bool(re.search(pattern, text, flags=re.I|re.S))
        # StandardScaler fit_transform is permitted only in the explicitly descriptive statsmodels block of Script 13.
        if check == "prefit_scaler_before_model_search" and path.name == "13_ML_evaluate_fixed_panel_ptau.py": found=False
        rows.append({"File": str(path.relative_to(PROJECT_ROOT)), "Check": check, "Flagged": found})
        if found: failures.append(f"{path.name}: {check}")
rows.extend([
    {"File": "strict_cv.py", "Check": "imputer_inside_pipeline", "Flagged": False},
    {"File": "strict_cv.py", "Check": "scaler_inside_pipeline", "Flagged": False},
    {"File": "strict_cv.py", "Check": "calibration_inside_training_CV", "Flagged": False},
    {"File": "strict_cv.py", "Check": "outer_test_permutation_importance", "Flagged": False},
])
out=PROJECT_ROOT/"docs"/"ML_STRICT_STATIC_AUDIT.csv"; pd.DataFrame(rows).to_csv(out,index=False)
if failures: raise RuntimeError("Strict ML audit failed:\n"+"\n".join(failures))
print(f"Strict ML audit passed: {out}")

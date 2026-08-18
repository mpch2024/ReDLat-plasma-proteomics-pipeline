###############################################################################
# ReDLat plasma proteomics — machine-learning reproducibility/static audit
#
# BLOCKING   = public-repository safety / portability violations
# DOCUMENTED = version-locked manuscript-reproduction implementation signatures
#              reported transparently but not treated as CI failures
###############################################################################

from pathlib import Path
import re
import sys
import pandas as pd

_REPO_HINT = Path(__file__).resolve()
for _candidate in [_REPO_HINT.parent, *_REPO_HINT.parents]:
    if (_candidate / ".redlat-root").exists() or (_candidate / ".here").exists():
        PROJECT_ROOT = _candidate
        break
else:
    raise FileNotFoundError("Repository root not found. Set REDLAT_PROJECT_ROOT.")

sys.path.insert(0, str(PROJECT_ROOT / "python"))

CORE = [
    PROJECT_ROOT / "scripts/ML/03_ML_run_primary_nested_cv.py",
    PROJECT_ROOT / "scripts/ML/06_ML_run_ptau_nested_cv.py",
    PROJECT_ROOT / "scripts/ML/09_ML_run_apoe_nested_cv.py",
    PROJECT_ROOT / "scripts/ML/12_ML_run_loco.py",
    PROJECT_ROOT / "scripts/ML/13_ML_evaluate_fixed_panel_ptau.py",
    PROJECT_ROOT / "scripts/ML/16_ML_run_matched_analysis.py",
    PROJECT_ROOT / "python/redlat_ml/strict_cv.py",
]

DOCUMENTED_PATTERNS = {
    "foldwise_score_z_normalization":
        r"y_prob\s*-\s*y_prob\.mean|probability\s*-\s*probability\.mean",
    "global_protein_complete_case_filter":
        r"dropna\s*\(\s*subset\s*=\s*(?:protein_cols|COMMON_COLS)",
    "prefit_scaler_before_model_search":
        r"StandardScaler\(\)\s*\n?.*fit_transform",
    "historical_logistic_regression_cv":
        r"LogisticRegressionCV",
}

BLOCKING_PATTERNS = {
    "personal_windows_path": r"[A-Za-z]:[/\\]Users[/\\]",
}

def active_runtime_install(text: str) -> bool:
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if re.search(r"(^|[;\s])%?pip\s+install\b", line, flags=re.I):
            return True
        if re.search(
            r"\b(?:subprocess\.(?:run|call|check_call|check_output)|os\.system)\b"
            r".*\b(?:pip|conda|mamba)\b.*\binstall\b",
            line,
            flags=re.I,
        ):
            return True
        if re.search(
            r"\binstall\.packages\s*\(|\bBiocManager::install\s*\(",
            line,
            flags=re.I,
        ):
            return True
    return False

rows = []
blocking_failures = []

for path in CORE:
    text = path.read_text(encoding="utf-8")
    rel = str(path.relative_to(PROJECT_ROOT))

    for check, pattern in DOCUMENTED_PATTERNS.items():
        found = bool(re.search(pattern, text, flags=re.I | re.S))
        rows.append({
            "File": rel,
            "Check": check,
            "Severity": "DOCUMENTED",
            "Flagged": found,
            "CI_blocking": False,
        })

    for check, pattern in BLOCKING_PATTERNS.items():
        found = bool(re.search(pattern, text, flags=re.I | re.S))
        rows.append({
            "File": rel,
            "Check": check,
            "Severity": "BLOCKING",
            "Flagged": found,
            "CI_blocking": found,
        })
        if found:
            blocking_failures.append(f"{path.name}: {check}")

    install_found = active_runtime_install(text)
    rows.append({
        "File": rel,
        "Check": "active_runtime_package_install",
        "Severity": "BLOCKING",
        "Flagged": install_found,
        "CI_blocking": install_found,
    })
    if install_found:
        blocking_failures.append(f"{path.name}: active_runtime_package_install")

strict_path = PROJECT_ROOT / "python/redlat_ml/strict_cv.py"
strict_text = strict_path.read_text(encoding="utf-8")
for check, pattern in {
    "strict_reference_pipeline_object": r"Pipeline\s*\(",
    "strict_reference_imputer": r"SimpleImputer\s*\(",
    "strict_reference_scaler": r"StandardScaler\s*\(",
    "strict_reference_calibration": r"CalibratedClassifierCV\s*\(",
}.items():
    rows.append({
        "File": str(strict_path.relative_to(PROJECT_ROOT)),
        "Check": check,
        "Severity": "REFERENCE",
        "Flagged": bool(re.search(pattern, strict_text, flags=re.I | re.S)),
        "CI_blocking": False,
    })

out_dir = PROJECT_ROOT / "result" / "audits"
out_dir.mkdir(parents=True, exist_ok=True)
out = out_dir / "ML_STATIC_AUDIT.csv"
pd.DataFrame(rows).to_csv(out, index=False)

documented = [
    f"{r['File']}: {r['Check']}"
    for r in rows
    if r["Severity"] == "DOCUMENTED" and r["Flagged"]
]

if blocking_failures:
    raise RuntimeError(
        "ML public-repository audit failed:\n" + "\n".join(blocking_failures)
    )

print(f"ML public-repository audit passed: {out}")
if documented:
    print("\nDocumented manuscript-reproduction signatures (non-blocking):")
    for item in documented:
        print(f"  - {item}")
print("\nNo blocking personal-path or active runtime-install violations detected.")

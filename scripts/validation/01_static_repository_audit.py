#!/usr/bin/env python3
from pathlib import Path
import os
import json
import re
import pandas as pd

HERE = Path(__file__).resolve()
ROOT = next(p for p in [HERE.parent, *HERE.parents] if (p / ".redlat-root").exists() or (p / ".here").exists())
OUT = ROOT / "outputs" / "validation"
OUT.mkdir(parents=True, exist_ok=True)
DATA = Path(os.environ.get("REDLAT_REVIEWER_DATA_DIR", ROOT / "data_private" / "reviewer_inputs")).expanduser().resolve()

analysis_roots = [ROOT / "scripts", ROOT / "R", ROOT / "config", ROOT / "python"]
paths = []
for base in analysis_roots:
    if not base.exists():
        continue
    paths.extend(base.rglob("*.R"))
    paths.extend(base.rglob("*.py"))

rows = []
for path in sorted(set(paths)):
    text = path.read_text(encoding="utf-8", errors="ignore")
    rows.append({
        "file": str(path.relative_to(ROOT)),
        "personal_absolute_path": bool(re.search(r"[A-Za-z]:[/\\]Users[/\\][A-Za-z0-9._-]+|/Users/[A-Za-z0-9._-]+/|/home/[A-Za-z0-9._-]+/", text)),
        "setwd_call": bool(re.search(r"^\s*setwd\s*\(", text, flags=re.M)),
        "active_read_adat_call": bool(re.search(r"SomaDataIO::read_adat\s*\(", text)),
        "runtime_install_command": bool(re.search(r"(^|\s)(install\.packages|BiocManager::install)\s*\(", text, flags=re.M)),
    })
audit = pd.DataFrame(rows)
audit.to_csv(OUT / "static_script_audit.csv", index=False)

forbidden_artifacts = []
for path in ROOT.rglob("*"):
    if not path.is_file():
        continue
    rel = str(path.relative_to(ROOT))
    name = path.name.lower()
    if '.bak' in name or name.endswith('.pyc') or name.startswith('patch_') or name.startswith('recover_historical_'):
        forbidden_artifacts.append(rel)

csvs = sorted(p.name for p in DATA.glob("*.csv"))
expected_csvs = sorted([
    "ReDLat_feature_annotation.csv",
    "ReDLat_metadata_deidentified.csv",
    "ReDLat_proteomics_somamer_log2.csv",
])


failures = {
    "personal_absolute_paths": int(audit["personal_absolute_path"].sum()),
    "setwd_calls": int(audit["setwd_call"].sum()),
    "active_read_adat_calls": int(audit["active_read_adat_call"].sum()),
    "runtime_install_commands": int(audit["runtime_install_command"].sum()),
    "forbidden_transient_artifacts": len(forbidden_artifacts),
    "data_csv_set_exact": csvs == expected_csvs,
}
summary = {
    "scripts_and_config_files_checked": int(len(audit)),
    **failures,
    "forbidden_transient_files": forbidden_artifacts,
    "data_csv_files": csvs,
}
(OUT / "static_script_audit.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
print(json.dumps(summary, indent=2))

hard_fail = (
    failures["personal_absolute_paths"] > 0
    or failures["setwd_calls"] > 0
    or failures["active_read_adat_calls"] > 0
    or failures["runtime_install_commands"] > 0
    or failures["forbidden_transient_artifacts"] > 0
    or not failures["data_csv_set_exact"]
)
if hard_fail:
    raise SystemExit("Static reviewer-package audit failed.")

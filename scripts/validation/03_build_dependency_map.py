#!/usr/bin/env python3
from pathlib import Path
import re
import pandas as pd
import json

HERE = Path(__file__).resolve()
ROOT = next(p for p in [HERE.parent, *HERE.parents] if (p/".redlat-root").exists() or (p/".here").exists())
OUT = ROOT/"outputs"/"validation"
OUT.mkdir(parents=True, exist_ok=True)

patterns = {
    "reads_csv": r"(read_csv|read\.csv|fread)\s*\(",
    "reads_xlsx": r"(read_excel|read\.xlsx|load_workbook)\s*\(",
    "reads_rds_workspace": r"(readRDS|load)\s*\(",
    "writes_files": r"(write_csv|write\.csv|fwrite|write_xlsx|write\.xlsx|saveRDS|save)\s*\(",
    "raw_adat_call": r"SomaDataIO::read_adat\s*\(",
    "setwd": r"^\s*setwd\s*\(",
}
rows=[]
analysis_paths = []
for analysis_root in [ROOT/"scripts"/"DEP", ROOT/"scripts"/"WGCNA", ROOT/"scripts"/"ML"]:
    analysis_paths.extend(analysis_root.rglob("*.R"))
    analysis_paths.extend(analysis_root.rglob("*.py"))
for path in sorted(analysis_paths):
    txt=path.read_text(encoding="utf-8", errors="ignore")
    row={"script":str(path.relative_to(ROOT))}
    for k,pat in patterns.items():
        row[k]=bool(re.search(pat,txt,flags=re.M))
    row["mentions_reviewer_proteomics"]="ReDLat_proteomics_somamer_log2.csv" in txt or "proteomics_file" in txt
    row["mentions_reviewer_metadata"]="ReDLat_metadata_deidentified.csv" in txt or "metadata_file" in txt
    row["mentions_reviewer_annotation"]="ReDLat_feature_annotation.csv" in txt or "annotation_file" in txt
    rows.append(row)

df=pd.DataFrame(rows)
df.to_csv(OUT/"script_dependency_map.csv",index=False)
summary={
    "scripts":len(df),
    "active_raw_adat_calls":int(df.raw_adat_call.sum()),
    "setwd_calls":int(df.setwd.sum()),
}
(OUT/"script_dependency_map.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
print(json.dumps(summary,indent=2))

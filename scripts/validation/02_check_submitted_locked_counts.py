#!/usr/bin/env python3
"""Checks locked processed-data quantities used by the reproducibility workflow."""
from pathlib import Path
import os
import pandas as pd

HERE = Path(__file__).resolve()
ROOT = next(p for p in [HERE.parent, *HERE.parents] if (p / ".redlat-root").exists() or (p / ".here").exists())
OUT = ROOT / "outputs" / "validation"
OUT.mkdir(parents=True, exist_ok=True)
DATA = Path(os.environ.get("REDLAT_REVIEWER_DATA_DIR", ROOT / "data_private" / "reviewer_inputs")).expanduser().resolve()
meta = pd.read_csv(DATA / "ReDLat_metadata_deidentified.csv")
ann = pd.read_csv(DATA / "ReDLat_feature_annotation.csv")
prot = pd.read_csv(DATA / "ReDLat_proteomics_somamer_log2.csv", nrows=2)

truthy = {"true", "1", "t", "yes", "y"}
def flag(x):
    if pd.api.types.is_bool_dtype(x): return x.fillna(False)
    if pd.api.types.is_numeric_dtype(x): return pd.to_numeric(x, errors="coerce").fillna(0).ne(0)
    return x.astype(str).str.strip().str.lower().isin(truthy)

primary = flag(meta["include_primary"])
matched = flag(meta["include_matched_selected"])
mlmap = flag(ann["Selected_for_ML_master"])
checks = [
    ("All participants", len(meta), 653),
    ("CN all", int((meta.SampleGroup == "CN").sum()), 327),
    ("AD all", int((meta.SampleGroup == "AD").sum()), 326),
    ("Primary complete cases", int(primary.sum()), 639),
    ("Primary CN", int(((meta.SampleGroup == "CN") & primary).sum()), 313),
    ("Primary AD", int(((meta.SampleGroup == "AD") & primary).sum()), 326),
    ("APOE analytic n", int(flag(meta["include_APOE"]).sum()), 579),
    ("AT(N) analytic n", int(flag(meta["include_ATN"]).sum()), 343),
    ("Selected matched n", int(matched.sum()), 382),
    ("Selected matched CN", int(((meta.SampleGroup == "CN") & matched).sum()), 191),
    ("Selected matched AD", int(((meta.SampleGroup == "AD") & matched).sum()), 191),
    ("SOMAmer features", len([c for c in prot.columns if c.startswith("seq.") or c.startswith("seq_")]), 10751),
    ("Gene-level universe", int(ann["EntrezGeneSymbol"].nunique()), 9638),
    ("DEP/WGCNA selected map", int(flag(ann["Selected_for_gene_matrix"]).sum()), 9638),
    ("Historical ML selected map", int(mlmap.sum()), 9638),
    ("Historical ML selected genes unique", int(ann.loc[mlmap, "EntrezGeneSymbol"].nunique()), 9638),
]
df = pd.DataFrame(checks, columns=["check", "observed", "expected"])
df["status"] = df.apply(lambda r: "PASS" if int(r.observed) == int(r.expected) else "FAIL", axis=1)
df.to_csv(OUT / "submitted_manuscript_locked_counts.csv", index=False)
print(df.to_string(index=False))
if (df.status == "FAIL").any():
    raise SystemExit("One or more locked-count checks failed.")

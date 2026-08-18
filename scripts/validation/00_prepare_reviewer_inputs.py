#!/usr/bin/env python3
"""
Processed-data reproducibility preflight and ML input construction.

Inputs: exactly the three processed/pseudonymized reviewer CSVs in data/.
Outputs under outputs/ML/private/derived/:
  - ReDLat_metadata_ML_compat.csv
  - ReDLat_ML_gene_master_RAW.csv

The ML master reproduces the historical 9,638-gene representative-SOMAmer map
stored as Selected_for_ML_master in ReDLat_feature_annotation.csv and converts
released log2 RFU back to normalized-RFU-compatible values using 2**x. This is
required because the submitted historical ML pipeline operated on that RAW-scale
gene master. DEP/WGCNA continue to use the released log2 data through the R
reviewer adapter.
"""
from pathlib import Path
import os
import json
import hashlib
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
ROOT = next(
    p for p in [HERE.parent, *HERE.parents]
    if (p / ".redlat-root").exists() or (p / ".here").exists()
)
DATA = Path(os.environ.get("REDLAT_REVIEWER_DATA_DIR", ROOT / "data_private" / "reviewer_inputs")).expanduser().resolve()
DERIVED = ROOT / "outputs" / "ML" / "private" / "derived"
VALID = ROOT / "outputs" / "validation"
DERIVED.mkdir(parents=True, exist_ok=True)
VALID.mkdir(parents=True, exist_ok=True)

META_FILE = DATA / "ReDLat_metadata_deidentified.csv"
PROT_FILE = DATA / "ReDLat_proteomics_somamer_log2.csv"
ANN_FILE = DATA / "ReDLat_feature_annotation.csv"

for path in (META_FILE, PROT_FILE, ANN_FILE):
    if not path.exists():
        raise FileNotFoundError(path)

meta = pd.read_csv(META_FILE, low_memory=False)
ann = pd.read_csv(ANN_FILE, low_memory=False)
prot = pd.read_csv(PROT_FILE, low_memory=False)

required_meta = [
    "Study_ID", "SampleGroup", "Site", "Country", "Sex", "Age", "Education",
    "APOE4_carrier", "cdr_global", "cdr_boxscore", "mmse_total", "udsfaq_total",
    "cog_benson", "cog_tmt_a", "cog_tmt_b", "cog_craft_verb_delayed",
    "NPI", "Mini.SEA", "T.ADLQ", "p.tau217", "p.tau181", "NfL",
    "ratio.AB42.40", "include_primary", "include_matched_selected",
]
missing_meta = [c for c in required_meta if c not in meta.columns]
if missing_meta:
    raise ValueError(f"Metadata missing: {missing_meta}")

if len(meta) != 653 or len(prot) != 653:
    raise ValueError(f"Expected 653 rows; metadata={len(meta)}, proteomics={len(prot)}")
meta["Study_ID"] = meta["Study_ID"].astype(str)
prot["Study_ID"] = prot["Study_ID"].astype(str)
if meta["Study_ID"].duplicated().any() or prot["Study_ID"].duplicated().any():
    raise ValueError("Duplicated Study_ID detected.")
if set(meta["Study_ID"]) != set(prot["Study_ID"]):
    raise ValueError("Study_ID sets differ between metadata and proteomics.")

seq_cols = [c for c in prot.columns if c.startswith("seq.") or c.startswith("seq_")]
if len(seq_cols) != 10751:
    raise ValueError(f"Expected 10,751 SOMAmer columns; observed {len(seq_cols)}")
if len(ann) != 10751:
    raise ValueError(f"Expected 10,751 annotation rows; observed {len(ann)}")

truthy = {"true", "1", "t", "yes", "y"}
def as_bool(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    if pd.api.types.is_numeric_dtype(series):
        return pd.to_numeric(series, errors="coerce").fillna(0).ne(0)
    return series.astype(str).str.strip().str.lower().isin(truthy)

# Historical ML representative-SOMAmer map, encoded directly in annotation.
if "Selected_for_ML_master" not in ann.columns:
    raise ValueError("Annotation is missing Selected_for_ML_master.")
ml_selected = ann.loc[
    as_bool(ann["Selected_for_ML_master"]),
    ["EntrezGeneSymbol", "AptName"],
].dropna().copy()
ml_selected["EntrezGeneSymbol"] = ml_selected["EntrezGeneSymbol"].astype(str)
ml_selected["AptName"] = ml_selected["AptName"].astype(str)
if len(ml_selected) != 9638:
    raise ValueError(f"Expected 9,638 historical ML representatives; observed {len(ml_selected)}")
if ml_selected["EntrezGeneSymbol"].duplicated().any():
    raise ValueError("Historical ML map contains duplicated gene symbols.")
if ml_selected["AptName"].duplicated().any():
    raise ValueError("Historical ML map contains duplicated AptName values.")
missing_apt = sorted(set(ml_selected["AptName"]) - set(prot.columns))
if missing_apt:
    raise ValueError(f"Historical ML AptNames missing from proteomics: {len(missing_apt)}")

# Primary ML cohort is the submitted 639-participant complete-case cohort.
primary_mask = as_bool(meta["include_primary"])
m = meta.loc[primary_mask].copy()
if len(m) != 639:
    raise ValueError(f"Expected 639 primary participants; observed {len(m)}")
cn = int((m["SampleGroup"] == "CN").sum())
ad = int((m["SampleGroup"] == "AD").sum())
if (cn, ad) != (313, 326):
    raise ValueError(f"Expected primary CN/AD=313/326; observed {cn}/{ad}")

# Exact pseudonymized matched subset used by the submitted matched sensitivity.
matched_mask = as_bool(meta["include_matched_selected"])
matched_all = meta.loc[matched_mask, ["Study_ID", "SampleGroup"]].copy()
if len(matched_all) != 382:
    raise ValueError(f"Expected 382 selected matched participants; observed {len(matched_all)}")
matched_counts = matched_all["SampleGroup"].value_counts()
if int(matched_counts.get("CN", 0)) != 191 or int(matched_counts.get("AD", 0)) != 191:
    raise ValueError("Selected matched cohort must be exactly 191 CN + 191 AD.")
if not set(matched_all["Study_ID"]).issubset(set(m["Study_ID"])):
    raise ValueError("Selected matched cohort is not a subset of the 639-person ML cohort.")

# ML-compatible metadata. SampleId is a pseudonymous compatibility alias only.
m["SampleId"] = m["Study_ID"]
m["ApoE"] = pd.Series(pd.NA, index=m.index, dtype="object")
m.loc[m["APOE4_carrier"].eq(1), "ApoE"] = "e3/e4"
m.loc[m["APOE4_carrier"].eq(0), "ApoE"] = "e3/e3"
m["APOE_group"] = "Unknown"
m.loc[m["APOE4_carrier"].eq(1), "APOE_group"] = "E4 carrier"
m.loc[m["APOE4_carrier"].eq(0), "APOE_group"] = "Non-E4"
m["GFAP_1"] = (
    pd.to_numeric(m["GFAP_biomarker"], errors="coerce")
    if "GFAP_biomarker" in m.columns else np.nan
)

ml_meta_cols = [
    "SampleId", "SampleGroup", "Site", "Country", "Sex", "Age", "Education",
    "ApoE", "APOE_group", "APOE4_carrier",
    "cdr_global", "cdr_boxscore", "mmse_total", "udsfaq_total",
    "cog_benson", "cog_tmt_a", "cog_tmt_b", "cog_craft_verb_delayed",
    "NPI", "Mini.SEA", "T.ADLQ",
    "p.tau217", "p.tau181", "NfL", "ratio.AB42.40", "GFAP_1",
]
meta_compat = m[ml_meta_cols + ["include_matched_selected"]].reset_index(drop=True)
meta_compat.to_csv(DERIVED / "ReDLat_metadata_ML_compat.csv", index=False)

# Preserve primary participant order, choose historical ML SOMAmer per gene,
# then invert the released log2 transform to normalized RFU.
prot_primary = (
    m[["Study_ID"]]
    .merge(prot, on="Study_ID", how="left", validate="one_to_one")
    .reset_index(drop=True)
)
apt_order = ml_selected["AptName"].tolist()
gene_order = ml_selected["EntrezGeneSymbol"].tolist()
log2_values = prot_primary[apt_order].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
raw_values = np.exp2(log2_values)
gene_values = pd.DataFrame(raw_values, columns=gene_order)
master = pd.concat([m[ml_meta_cols].reset_index(drop=True), gene_values], axis=1)
if master.shape != (639, 9664):
    raise ValueError(f"Expected ML master shape (639, 9664); observed {master.shape}")
if master.columns.duplicated().any():
    dup = master.columns[master.columns.duplicated()].tolist()
    raise ValueError(f"Duplicated columns in ML master: {dup[:20]}")
master_path = DERIVED / "ReDLat_ML_gene_master_RAW.csv"
master.to_csv(master_path, index=False, float_format="%.12g")

panel = ["SPC25", "CPLX2", "TCP11L1", "ACHE", "ODC1", "SPON1", "RTN4RL1"]
missing_panel = [x for x in panel if x not in master.columns]

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

report = {
    "participants": int(len(meta)),
    "CN_all": int((meta["SampleGroup"] == "CN").sum()),
    "AD_all": int((meta["SampleGroup"] == "AD").sum()),
    "primary_n": int(primary_mask.sum()),
    "primary_CN": cn,
    "primary_AD": ad,
    "APOE_n": int(as_bool(meta["include_APOE"]).sum()) if "include_APOE" in meta else int(meta["APOE4_carrier"].notna().sum()),
    "ATN_n": int(as_bool(meta["include_ATN"]).sum()) if "include_ATN" in meta else None,
    "matched_selected_n": int(matched_mask.sum()),
    "matched_selected_CN": int((matched_all["SampleGroup"] == "CN").sum()),
    "matched_selected_AD": int((matched_all["SampleGroup"] == "AD").sum()),
    "somamers": len(seq_cols),
    "historical_ML_gene_features": int(len(ml_selected)),
    "ML_master_shape": list(master.shape),
    "ML_master_scale": "normalized RFU reconstructed as 2**released_log2_RFU",
    "seven_protein_panel_present": len(missing_panel) == 0,
    "missing_panel": missing_panel,
    "sha256": {
        META_FILE.name: sha256(META_FILE),
        PROT_FILE.name: sha256(PROT_FILE),
        ANN_FILE.name: sha256(ANN_FILE),
    },
}
(VALID / "reviewer_data_preflight.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
pd.DataFrame([{k: v for k, v in report.items() if k != "sha256"}]).to_csv(
    VALID / "reviewer_data_preflight.csv", index=False
)

print(json.dumps(report, indent=2))
print("\nGenerated:")
print(DERIVED / "ReDLat_metadata_ML_compat.csv")
print(master_path)

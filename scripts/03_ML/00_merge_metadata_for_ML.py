import pandas as pd

PROJECT_DIR = r"C:\Clasificacion"

DATA_DIR = rf"{PROJECT_DIR}\Data"

BASE_FILE = rf"{DATA_DIR}\somascan_filter_DEP_final_FDR.csv"

META_FILE = rf"{DATA_DIR}\metadata.csv"

OUT_FILE = rf"{DATA_DIR}\somascan_filter_DEP_final_FDR_metadata.csv"

# ------------------------
# LOAD
# ------------------------

base = pd.read_csv(BASE_FILE)

meta = pd.read_csv(META_FILE)

# ------------------------
# ID
# ------------------------

base["SampleId"] = (base["SampleId"].astype(str))

meta["SampleId"] = (meta["SampleId"].astype(str))

# ------------------------
# KEEP ONLY NEW COLUMNS
# ------------------------

meta_cols = [
    c
    for c in meta.columns
    if c != "SampleId"
]

# ------------------------
# MERGE
# ------------------------

merged = base.merge(meta[["SampleId"] + meta_cols],on="SampleId",how="left")

# ------------------------
# REMOVE SUBJECT
# ------------------------

merged = merged[merged["SampleId"] != "MA205"]

# ------------------------
# CHECK
# ------------------------

print("\nBase:",base.shape)
print("Meta:",meta.shape)
print("Merged:", merged.shape)
print("\nSubjects after drop:")
print( merged.shape[0])
print("\nMissing metadata")
print(merged[ meta_cols].isna().sum())


# ------------------------
# SAVE
# ------------------------

merged.to_csv(OUT_FILE,index=False)

print("\nSaved")
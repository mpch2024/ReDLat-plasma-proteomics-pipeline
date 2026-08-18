# Processed reviewer-data interface

The reproducibility rerun was validated using a minimal three-file, processed and pseudonymized data interface. The files are not distributed through GitHub.

## Required local files

| File | Role |
|---|---|
| `ReDLat_proteomics_somamer_log2.csv` | 653 × 10,751 SOMAmer log2-transformed normalized RFU matrix plus pseudonymous study ID |
| `ReDLat_metadata_deidentified.csv` | Pseudonymized metadata and analysis-inclusion flags required by DEP/WGCNA/ML |
| `ReDLat_feature_annotation.csv` | SOMAmer/protein annotation and 9,638-gene feature mappings |

The default local location is `data_private/reviewer_inputs/`. An approved alternative directory can be supplied with `REDLAT_REVIEWER_DATA_DIR`.

## Compatibility logic

Legacy DEP/WGCNA code historically operated on normalized RFU and applied one log2 transformation internally. The R compatibility layer therefore supplies `2^(released log2 RFU)` in memory, preserving the original single-transform workflow without requiring the original ADAT.

The historical ML workflow used a 9,638-gene normalized-RFU master. `scripts/validation/00_prepare_reviewer_inputs.py` reconstructs that local master from the feature marked `Selected_for_ML_master` for each gene and `2**released_log2_RFU` values. The local derivative is written under `outputs/ML/private/derived/` and must remain untracked.

The exact pseudonymized 191 CN + 191 AD matched subset used by the submitted matched sensitivity is encoded by the controlled metadata flag `include_matched_selected`; no original linkage key is required by the public code.

## Privacy boundary

The public repository must not contain participant-level data or pseudonymous IDs. Publication/source-data scripts include privacy checks that reject participant identifiers, ReDLat pseudonyms, and personal filesystem paths from public outputs.

# Governed reviewer inputs

This directory is intentionally empty in the public repository.

For an authorized local rerun, place exactly these three processed/pseudonymized files here:

1. `ReDLat_proteomics_somamer_log2.csv`
2. `ReDLat_metadata_deidentified.csv`
3. `ReDLat_feature_annotation.csv`

Alternatively set `REDLAT_REVIEWER_DATA_DIR` to the approved local directory containing those files.

Do **not** commit the files themselves, raw ADAT data, participant linkage keys, original identifiers, credentials, or governed local paths.

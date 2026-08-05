# Local configuration

R and Python share paths and strict-analysis settings through environment variables. Copy `.Renviron.example` to a private `.Renviron`, enter only governed local paths, and keep that file untracked.

Required inputs:

- `REDLAT_ML_METADATA_FILE`
- `REDLAT_ML_PTAU_METADATA_FILE`
- `REDLAT_ML_ADAT_FILE`
- `REDLAT_ML_MASTER_FILE`
- `REDLAT_ML_EXCLUDED_IDS_FILE`, when exclusions are governed separately

The matching outputs default to `result/ML/private/matching/`. Override `REDLAT_ML_MATCHED_IDS_FILE` and `REDLAT_ML_MATCHED_FULL_FILE` only when the selected outputs are stored in another governed private location.

Strict nested-CV parameters are also listed in `.Renviron.example`. Record any deviations in the internal software-environment archive and in the analysis manifest.

Never place participant identifiers, credentials, exclusions or local absolute paths in tracked source files.

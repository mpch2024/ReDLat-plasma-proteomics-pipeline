# Copy this file to config/dep_config.R and edit only the local paths.
# config/dep_config.R is ignored by Git and must never contain credentials.

list(
  data_dir = file.path(project_root, "data_private"),
  metadata_file = file.path(project_root, "data_private", "clinical_metadata.csv"),
  adat_file = file.path(project_root, "data_private", "proteomics.adat"),
  result_root = file.path(project_root, "result"),
  publication_root = file.path(project_root, "publication_candidate", "DEP"),
  allow_participant_level_exports = FALSE
)

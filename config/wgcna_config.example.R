# Copy to config/wgcna_config.R and edit only local paths.
# The local file is ignored by Git and must not contain credentials.

WGCNA_LOCAL_CONFIG <- list(
  dep_project_root = project_root,
  dep_result_root = file.path(project_root, "result"),
  result_root = file.path(project_root, "result", "WGCNA"),
  publication_root = file.path(project_root, "publication_candidate", "WGCNA"),
  allow_participant_level_public_exports = FALSE
)

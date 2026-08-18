ml_project_config <- function(project_root) {
  data_dir <- Sys.getenv(
    "REDLAT_REVIEWER_DATA_DIR",
    unset = file.path(project_root, "data_private", "reviewer_inputs")
  )
  private_root <- file.path(project_root, "outputs", "ML", "private")
  list(
    project_root = project_root,
    data_dir = data_dir,
    metadata_file = file.path(private_root, "derived", "ReDLat_metadata_ML_compat.csv"),
    master_file = file.path(private_root, "derived", "ReDLat_ML_gene_master_RAW.csv"),
    proteomics_file = file.path(data_dir, "ReDLat_proteomics_somamer_log2.csv"),
    annotation_file = file.path(data_dir, "ReDLat_feature_annotation.csv"),
    # Legacy compatibility alias only; this is NOT an ADAT.
    adat_file = file.path(data_dir, "ReDLat_proteomics_somamer_log2.csv"),
    no_exclusions_file = "",
    matching_focus_id = "",
    private_root = private_root,
    publication_root = file.path(project_root, "outputs", "ML", "publication"),
    public_root = file.path(project_root, "outputs", "ML", "publication")
  )
}

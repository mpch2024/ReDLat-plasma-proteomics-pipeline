wgcna_project_config <- function(project_root) {
  list(
    project_root = project_root,
    dep_project_root = project_root,
    dep_result_root = file.path(project_root, "outputs", "DEP"),
    result_root = file.path(project_root, "outputs", "WGCNA"),
    publication_root = file.path(project_root, "outputs", "publication", "WGCNA")
  )
}

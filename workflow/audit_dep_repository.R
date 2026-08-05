# ReDLat DEP workflow — repository privacy audit
project_root <- Sys.getenv("REDLAT_PROJECT_ROOT", unset = normalizePath(getwd(), winslash = "/", mustWork = TRUE))
source(file.path(project_root, "R", "dep_bootstrap.R"))
hits <- dep_scan_text_for_private_paths(project_root)
if (length(hits) > 0L) {
  details <- unlist(Map(function(file, lines) paste0(file, ":", lines), names(hits), hits))
  stop("Private/local path patterns were found:\n", paste(details, collapse = "\n"), call. = FALSE)
}
tracked_candidates <- c(file.path(project_root, "data_private"), file.path(project_root, "result"), file.path(project_root, "publication_candidate"))
message("Repository path audit passed. Generated data directories remain git-ignored.")

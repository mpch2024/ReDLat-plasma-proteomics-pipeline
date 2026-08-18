wgcna_load_config <- function(project_root) {
  cfg <- file.path(project_root, "config", "wgcna_config.R")
  if (!file.exists(cfg)) stop("WGCNA config not found: ", cfg, call. = FALSE)
  source(cfg, local = TRUE)
  out <- wgcna_project_config(project_root)
  dir.create(out$result_root, recursive = TRUE, showWarnings = FALSE)
  dir.create(out$publication_root, recursive = TRUE, showWarnings = FALSE)
  out
}

wgcna_assert_public_tree <- function(root) {
  if (!dir.exists(root)) stop("Public tree not found: ", root, call. = FALSE)
  files <- list.files(root, recursive = TRUE, full.names = TRUE)
  bad_names <- grep(
    "SampleId|SubjectID|linkage|PRIVATE_DO_NOT_SHARE|\\.adat$",
    basename(files), ignore.case = TRUE, value = TRUE
  )
  if (length(bad_names)) {
    stop("Public WGCNA tree contains privacy-sensitive filenames: ",
         paste(bad_names, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

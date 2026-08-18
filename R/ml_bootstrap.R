ml_load_config <- function(project_root) {
  cfg <- file.path(project_root, "config", "ml_config.R")
  if (!file.exists(cfg)) stop("ML config not found: ", cfg, call. = FALSE)
  source(cfg, local = TRUE)
  out <- ml_project_config(project_root)
  dir.create(out$private_root, recursive = TRUE, showWarnings = FALSE)
  dir.create(out$publication_root, recursive = TRUE, showWarnings = FALSE)
  out
}

ml_read_excluded_ids <- function(config) {
  if (is.null(config$no_exclusions_file) || !nzchar(config$no_exclusions_file) ||
      !file.exists(config$no_exclusions_file)) {
    return(character(0))
  }
  x <- readr::read_csv(config$no_exclusions_file, show_col_types = FALSE)
  id_col <- intersect(c("SampleId", "Study_ID"), names(x))[1]
  if (is.na(id_col) || length(id_col) == 0L) return(character(0))
  unique(as.character(x[[id_col]]))
}

ml_assert_private_path <- function(path, config) {
  root <- normalizePath(config$private_root, winslash = "/", mustWork = FALSE)
  target <- normalizePath(path, winslash = "/", mustWork = FALSE)
  if (!startsWith(target, root)) {
    stop("ML participant-level output must remain under private_root: ", target,
         call. = FALSE)
  }
  invisible(TRUE)
}

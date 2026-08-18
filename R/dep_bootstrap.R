dep_load_config <- function(project_root) {
  cfg <- file.path(project_root, "config", "dep_config.R")
  if (!file.exists(cfg)) stop("DEP config not found: ", cfg, call. = FALSE)
  source(cfg, local = TRUE)
  out <- dep_project_config(project_root)
  dir.create(out$result_root, recursive = TRUE, showWarnings = FALSE)
  dir.create(out$publication_root, recursive = TRUE, showWarnings = FALSE)
  out
}

dep_assert_no_direct_identifiers <- function(x, label = "object") {
  forbidden <- tolower(c(
    "SubjectID", "Subject_ID", "original_sampleid", "barcode",
    "scannerid", "plateposition", "run_date", "date_of_birth"
  ))
  hit <- intersect(tolower(names(x)), forbidden)
  if (length(hit)) stop(label, " contains forbidden direct/private columns: ",
                        paste(hit, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

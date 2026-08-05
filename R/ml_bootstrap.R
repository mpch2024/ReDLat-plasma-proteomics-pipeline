# Shared infrastructure for the machine-learning workflow.

ml_env_flag <- function(name, default = FALSE) {
  value <- tolower(trimws(Sys.getenv(name, unset = as.character(default))))
  value %in% c("1", "true", "t", "yes", "y")
}

ml_project_root <- function() {
  explicit <- trimws(Sys.getenv("REDLAT_PROJECT_ROOT", unset = ""))
  if (nzchar(explicit)) return(normalizePath(explicit, winslash = "/", mustWork = TRUE))
  if (!requireNamespace("here", quietly = TRUE)) {
    stop("Package 'here' is required. Restore the project environment before running the workflow.", call. = FALSE)
  }
  normalizePath(here::here(), winslash = "/", mustWork = TRUE)
}

ml_load_config <- function(project_root = ml_project_root()) {
  data_dir <- Sys.getenv("REDLAT_ML_DATA_DIR", unset = file.path(project_root, "data_private"))
  result_root <- Sys.getenv("REDLAT_ML_RESULT_ROOT", unset = file.path(project_root, "result", "ML"))
  publication_root <- Sys.getenv(
    "REDLAT_ML_PUBLICATION_ROOT",
    unset = file.path(project_root, "publication_candidate", "ML")
  )
  config <- list(
    project_root = project_root,
    data_dir = normalizePath(data_dir, winslash = "/", mustWork = FALSE),
    metadata_file = normalizePath(
      Sys.getenv("REDLAT_ML_METADATA_FILE", unset = file.path(data_dir, "clinical_metadata.csv")),
      winslash = "/", mustWork = FALSE
    ),
    adat_file = normalizePath(
      Sys.getenv("REDLAT_ML_ADAT_FILE", unset = file.path(data_dir, "proteomics.adat")),
      winslash = "/", mustWork = FALSE
    ),
    master_file = normalizePath(
      Sys.getenv("REDLAT_ML_MASTER_FILE", unset = file.path(data_dir, "ml_master_matrix.csv")),
      winslash = "/", mustWork = FALSE
    ),
    ptau_metadata_file = normalizePath(
      Sys.getenv("REDLAT_ML_PTAU_METADATA_FILE", unset = Sys.getenv(
        "REDLAT_ML_METADATA_FILE", unset = file.path(data_dir, "clinical_metadata.csv")
      )),
      winslash = "/", mustWork = FALSE
    ),
    matched_ids_file = normalizePath(
      Sys.getenv("REDLAT_ML_MATCHED_IDS_FILE", unset = file.path(
        result_root, "private", "matching", "matched_ids_SELECTED.csv"
      )),
      winslash = "/", mustWork = FALSE
    ),
    matched_full_file = normalizePath(
      Sys.getenv("REDLAT_ML_MATCHED_FULL_FILE", unset = file.path(
        result_root, "private", "matching", "Matched_Output_SELECTED.csv"
      )),
      winslash = "/", mustWork = FALSE
    ),
    no_exclusions_file = trimws(Sys.getenv("REDLAT_ML_NO_EXCLUSIONS_FILE", unset = "")),
    excluded_ids_file = trimws(Sys.getenv("REDLAT_ML_EXCLUDED_IDS_FILE", unset = "")),
    result_root = normalizePath(result_root, winslash = "/", mustWork = FALSE),
    publication_root = normalizePath(publication_root, winslash = "/", mustWork = FALSE),
    allow_participant_exports = ml_env_flag("REDLAT_ALLOW_PARTICIPANT_LEVEL_EXPORTS", FALSE),
    matching_focus_id = trimws(Sys.getenv("REDLAT_ML_MATCHING_FOCUS_ID", unset = ""))
  )
  config$private_root <- file.path(config$result_root, "private")
  config$public_root <- file.path(config$result_root, "public")
  invisible(lapply(
    c(config$data_dir, config$result_root, config$private_root, config$public_root, config$publication_root),
    dir.create, recursive = TRUE, showWarnings = FALSE
  ))
  config
}

ml_require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop(
      "Missing required R packages: ", paste(missing, collapse = ", "),
      ". Restore the locked R environment before running this script.",
      call. = FALSE
    )
  }
  invisible(packages)
}

ml_read_excluded_ids <- function(config) {
  path <- config$excluded_ids_file
  if (!nzchar(path) || !file.exists(path)) return(character(0))
  ids <- trimws(readLines(path, warn = FALSE))
  ids <- sub(",.*$", "", ids)
  ids[nzchar(ids) & !tolower(ids) %in% c("sampleid", "sample_id", "id")]
}

ml_assert_private_path <- function(path, config) {
  private_root <- normalizePath(config$private_root, winslash = "/", mustWork = FALSE)
  target <- normalizePath(path, winslash = "/", mustWork = FALSE)
  if (!startsWith(target, paste0(private_root, "/")) && target != private_root) {
    stop("Participant-level output must remain under the private result root: ", target, call. = FALSE)
  }
  invisible(path)
}

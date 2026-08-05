# ReDLat DEP workflow — shared configuration and privacy functions

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || !nzchar(as.character(x)[1])) y else x

dep_env_flag <- function(name, default = FALSE) {
  value <- tolower(trimws(Sys.getenv(name, unset = as.character(default))))
  value %in% c("1", "true", "t", "yes", "y")
}

dep_load_config <- function(project_root) {
  project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
  local_config <- file.path(project_root, "config", "dep_config.R")
  config <- list(
    data_dir = Sys.getenv("REDLAT_DEP_DATA_DIR", unset = file.path(project_root, "data_private")),
    metadata_file = Sys.getenv("REDLAT_DEP_METADATA_FILE", unset = ""),
    adat_file = Sys.getenv("REDLAT_DEP_ADAT_FILE", unset = ""),
    result_root = Sys.getenv("REDLAT_DEP_RESULT_ROOT", unset = file.path(project_root, "result")),
    publication_root = Sys.getenv("REDLAT_DEP_PUBLICATION_ROOT", unset = file.path(project_root, "publication_candidate", "DEP")),
    allow_participant_level_exports = dep_env_flag("REDLAT_ALLOW_PARTICIPANT_LEVEL_EXPORTS", FALSE)
  )
  if (file.exists(local_config)) {
    config_env <- new.env(parent = baseenv())
    config_env$project_root <- project_root
    local_values <- source(local_config, local = config_env)$value
    if (!is.list(local_values)) stop("config/dep_config.R must return a named list.", call. = FALSE)
    config[names(local_values)] <- local_values
  }
  if (!nzchar(config$metadata_file)) config$metadata_file <- file.path(config$data_dir, "clinical_metadata.csv")
  if (!nzchar(config$adat_file)) config$adat_file <- file.path(config$data_dir, "proteomics.adat")
  for (nm in c("data_dir", "result_root", "publication_root")) {
    config[[nm]] <- normalizePath(config[[nm]], winslash = "/", mustWork = FALSE)
  }
  config$metadata_file <- normalizePath(config$metadata_file, winslash = "/", mustWork = FALSE)
  config$adat_file <- normalizePath(config$adat_file, winslash = "/", mustWork = FALSE)
  dir.create(config$result_root, recursive = TRUE, showWarnings = FALSE)
  dir.create(config$publication_root, recursive = TRUE, showWarnings = FALSE)
  config
}

dep_require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop("Missing required packages: ", paste(missing, collapse = ", "),
         ". Restore the locked environment with renv::restore().", call. = FALSE)
  }
  invisible(packages)
}

dep_sensitive_column_pattern <- paste(
  c("^sample_?id$", "^participant_?id$", "^subject_?id$", "^record_?id$",
    "^patient_?id$", "^individual_?id$", "^family_?id$", "^subclass$",
    "^match_pair_id$", "^source_row$"), collapse = "|")

dep_assert_no_direct_identifiers <- function(x, label = "table") {
  if (is.null(x) || !is.data.frame(x)) return(invisible(TRUE))
  bad <- names(x)[grepl(dep_sensitive_column_pattern, names(x), ignore.case = TRUE)]
  # Synthetic Record columns created by the public Source Data generator are allowed.
  bad <- setdiff(bad, "Record")
  if (length(bad) > 0L) {
    stop(label, " contains direct/stable identifier columns: ", paste(bad, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

dep_scan_text_for_private_paths <- function(root) {
  files <- list.files(
    root,
    recursive = TRUE,
    full.names = TRUE,
    pattern = "\\.(R|Rmd|py|md|txt|ya?ml|json)$",
    ignore.case = TRUE
  )
  
  # Ignore package libraries restored locally or by GitHub Actions.
  files <- files[
    !grepl(
      "[/\\\\]renv[/\\\\](library|staging)[/\\\\]",
      files,
      perl = TRUE
    )
  ]
  
  # These files contain the privacy-detection patterns themselves.
  scanner_files <- c(
    "dep_bootstrap.R",
    "wgcna_bootstrap.R",
    "ml_bootstrap.R",
    "audit_dep_repository.R",
    "audit_publication_candidate.R",
    "audit_wgcna_publication.R"
  )
  
  patterns <- c(
    "[A-Za-z]:[/\\\\]Users[/\\\\]",
    "/Users/",
    "/home/[^/]+/",
    "OneDrive",
    "Desktop"
  )
  
  hits <- list()
  
  for (f in files) {
    if (basename(f) %in% scanner_files) next
    
    lines <- tryCatch(
      readLines(f, warn = FALSE),
      error = function(e) character(0)
    )
    
    idx <- which(
      vapply(
        lines,
        function(line) {
          any(
            vapply(
              patterns,
              function(pattern) grepl(pattern, line, perl = TRUE),
              logical(1)
            )
          )
        },
        logical(1)
      )
    )
    
    if (length(idx) > 0L) {
      hits[[f]] <- idx
    }
  }
  
  hits
}


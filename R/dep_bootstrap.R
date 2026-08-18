# ReDLat DEP workflow — shared configuration and privacy functions

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || !nzchar(as.character(x)[1])) y else x
}

dep_env_flag <- function(name, default = FALSE) {
  value <- tolower(trimws(Sys.getenv(name, unset = as.character(default))))
  value %in% c("1", "true", "t", "yes", "y")
}

dep_load_config <- function(project_root) {
  project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
  cfg <- file.path(project_root, "config", "dep_config.R")
  if (!file.exists(cfg)) stop("DEP config not found: ", cfg, call. = FALSE)

  source(cfg, local = TRUE)
  out <- dep_project_config(project_root)

  if (!is.list(out)) {
    stop("dep_project_config() must return a named list.", call. = FALSE)
  }

  required <- c(
    "data_dir", "metadata_file", "adat_file",
    "result_root", "publication_root"
  )
  missing <- setdiff(required, names(out))
  if (length(missing)) {
    stop(
      "DEP config is missing required fields: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  dir.create(out$result_root, recursive = TRUE, showWarnings = FALSE)
  dir.create(out$publication_root, recursive = TRUE, showWarnings = FALSE)
  out
}

dep_require_packages <- function(packages) {
  missing <- packages[
    !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing) > 0L) {
    stop(
      "Missing required packages: ",
      paste(missing, collapse = ", "),
      ". Restore the locked environment with renv::restore().",
      call. = FALSE
    )
  }
  invisible(packages)
}

dep_sensitive_column_pattern <- paste(
  c(
    "^sample_?id$",
    "^participant_?id$",
    "^subject_?id$",
    "^record_?id$",
    "^patient_?id$",
    "^individual_?id$",
    "^family_?id$",
    "^original_sampleid$",
    "^barcode$",
    "^scannerid$",
    "^plateposition$",
    "^run_date$",
    "^date_of_birth$",
    "^subclass$",
    "^match_pair_id$",
    "^source_row$"
  ),
  collapse = "|"
)

dep_assert_no_direct_identifiers <- function(x, label = "table") {
  if (is.null(x) || !is.data.frame(x)) return(invisible(TRUE))

  bad <- names(x)[
    grepl(
      dep_sensitive_column_pattern,
      names(x),
      ignore.case = TRUE
    )
  ]

  bad <- setdiff(bad, "Record")

  if (length(bad) > 0L) {
    stop(
      label,
      " contains direct/stable identifier columns: ",
      paste(bad, collapse = ", "),
      call. = FALSE
    )
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

  files <- files[
    !grepl(
      "[/\\\\]renv[/\\\\](library|staging)[/\\\\]",
      files,
      perl = TRUE
    )
  ]

  scanner_files <- c(
    "dep_bootstrap.R",
    "wgcna_bootstrap.R",
    "ml_bootstrap.R",
    "audit_dep_repository.R",
    "audit_publication_candidate.R",
    "audit_wgcna_publication.R",
    "24_ML_audit_publication_outputs.py",
    "25_ML_audit_strict_pipeline.py",
    "00_capture_software_environment.py",
    "01_static_repository_audit.py",
    "config.py"
  )

  patterns <- c(
    "[A-Za-z]:[/\\\\]Users[/\\\\]",
    "/Users/",
    "/home/[^/]+/",
    "OneDrive"
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

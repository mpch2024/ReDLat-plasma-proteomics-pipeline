# Shared infrastructure for the WGCNA workflow.

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || !nzchar(as.character(x)[1])) y else x
}

wgcna_env_flag <- function(name, default = FALSE) {
  value <- tolower(trimws(Sys.getenv(name, unset = as.character(default))))
  value %in% c("1", "true", "t", "yes", "y")
}

wgcna_load_config <- function(project_root) {
  project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
  config <- list(
    project_root = project_root,
    dep_project_root = Sys.getenv("REDLAT_DEP_PROJECT_ROOT", unset = project_root),
    dep_result_root = Sys.getenv("REDLAT_DEP_RESULT_ROOT", unset = file.path(project_root, "result")),
    result_root = Sys.getenv("REDLAT_WGCNA_RESULT_ROOT", unset = file.path(project_root, "result", "WGCNA")),
    publication_root = Sys.getenv("REDLAT_WGCNA_PUBLICATION_ROOT", unset = file.path(project_root, "publication_candidate", "WGCNA")),
    allow_participant_level_public_exports = wgcna_env_flag("REDLAT_ALLOW_PARTICIPANT_LEVEL_EXPORTS", FALSE)
  )

  local_config <- file.path(project_root, "config", "wgcna_config.R")
  if (file.exists(local_config)) {
    env <- new.env(parent = baseenv())
    sys.source(local_config, envir = env)
    values <- if (exists("WGCNA_LOCAL_CONFIG", envir = env, inherits = FALSE)) {
      get("WGCNA_LOCAL_CONFIG", envir = env)
    } else if (exists("config", envir = env, inherits = FALSE)) {
      get("config", envir = env)
    } else {
      stop("config/wgcna_config.R must define WGCNA_LOCAL_CONFIG or config as a named list.", call. = FALSE)
    }
    if (!is.list(values)) stop("The local WGCNA configuration must be a named list.", call. = FALSE)
    config[names(values)] <- values
  }

  for (nm in c("project_root", "dep_project_root", "dep_result_root", "result_root", "publication_root")) {
    config[[nm]] <- normalizePath(config[[nm]], winslash = "/", mustWork = FALSE)
  }
  dir.create(config$result_root, recursive = TRUE, showWarnings = FALSE)
  dir.create(config$publication_root, recursive = TRUE, showWarnings = FALSE)
  config
}

wgcna_require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop("Missing required packages: ", paste(missing, collapse = ", "),
         ". Restore the locked environment with renv::restore().", call. = FALSE)
  }
  invisible(packages)
}

wgcna_direct_identifier_pattern <- paste(
  c("^sample_?id$", "^participant_?id$", "^subject_?id$", "^patient_?id$",
    "^individual_?id$", "^record_?id$", "^family_?id$", "^source_row$",
    "^subclass$", "^match_pair_id$"),
  collapse = "|"
)

wgcna_assert_no_direct_identifiers <- function(x, label = "table") {
  if (is.null(x) || !is.data.frame(x)) return(invisible(TRUE))
  bad <- names(x)[grepl(wgcna_direct_identifier_pattern, names(x), ignore.case = TRUE)]
  bad <- setdiff(bad, c("Record", "record"))
  if (length(bad) > 0L) {
    stop(label, " contains direct or stable identifier columns: ",
         paste(bad, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

wgcna_scan_workbook_headers <- function(path) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) return(character(0))
  sheets <- tryCatch(openxlsx::getSheetNames(path), error = function(e) character(0))
  hits <- character(0)
  for (sheet in sheets) {
    raw <- tryCatch(
      openxlsx::read.xlsx(path, sheet = sheet, rows = 1:15, colNames = FALSE,
                          skipEmptyRows = FALSE, skipEmptyCols = FALSE),
      error = function(e) NULL
    )
    if (is.null(raw)) next
    cells <- trimws(tolower(as.character(unlist(raw, use.names = FALSE))))
    exact <- c("sampleid", "sample_id", "participantid", "participant_id",
               "subjectid", "subject_id", "patientid", "patient_id",
               "subclass", "match_pair_id")
    if (any(cells %in% exact, na.rm = TRUE)) hits <- c(hits, paste0(basename(path), "::", sheet))
  }
  unique(hits)
}

wgcna_assert_public_tree <- function(root) {
  if (!dir.exists(root)) stop("Public output directory does not exist: ", root, call. = FALSE)
  files <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = TRUE)
  files <- files[file.info(files)$isdir %in% FALSE]

  path_patterns <- c("[A-Za-z]:[/\\\\]Users[/\\\\]", "/Users/", "/home/[^/]+/", "OneDrive", "Desktop")
  text_files <- files[grepl("\\.(txt|md|csv|tsv|json|ya?ml|R|py)$", files, ignore.case = TRUE)]
  path_hits <- character(0)
  identifier_hits <- character(0)

  for (file in text_files) {
    lines <- tryCatch(readLines(file, warn = FALSE), error = function(e) character(0))
    if (length(lines) && any(vapply(lines, function(line) any(vapply(path_patterns, grepl, logical(1), x = line, perl = TRUE)), logical(1)))) {
      path_hits <- c(path_hits, file)
    }
    if (grepl("\\.(csv|tsv)$", file, ignore.case = TRUE)) {
      sep <- if (grepl("\\.tsv$", file, ignore.case = TRUE)) "\\t" else ","
      header <- tryCatch(strsplit(readLines(file, n = 1L, warn = FALSE), sep)[[1]], error = function(e) character(0))
      bad <- header[grepl(wgcna_direct_identifier_pattern, header, ignore.case = TRUE)]
      bad <- setdiff(bad, c("Record", "record"))
      if (length(bad)) identifier_hits <- c(identifier_hits, paste0(file, ": ", paste(bad, collapse = ", ")))
    }
  }

  xlsx_files <- files[grepl("\\.xlsx$", files, ignore.case = TRUE)]
  workbook_hits <- unlist(lapply(xlsx_files, wgcna_scan_workbook_headers), use.names = FALSE)

  if (length(path_hits) || length(identifier_hits) || length(workbook_hits)) {
    stop(
      "Public-output audit failed.",
      if (length(path_hits)) paste0("\nPrivate path references: ", paste(path_hits, collapse = "; ")) else "",
      if (length(identifier_hits)) paste0("\nIdentifier columns: ", paste(identifier_hits, collapse = "; ")) else "",
      if (length(workbook_hits)) paste0("\nWorkbook identifier headers: ", paste(workbook_hits, collapse = "; ")) else "",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

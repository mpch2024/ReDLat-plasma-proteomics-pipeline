# ReDLat DEP workflow — publication-file audit
project_root <- Sys.getenv(
  "REDLAT_PROJECT_ROOT",
  unset = normalizePath(getwd(), winslash = "/", mustWork = TRUE)
)
source(file.path(project_root, "R", "dep_bootstrap.R"))
config <- dep_load_config(project_root)
root <- config$publication_root
if (!dir.exists(root)) stop("Publication candidate directory does not exist: ", root, call. = FALSE)

sensitive_tokens <- c(
  "SampleId", "SampleID", "ParticipantId", "ParticipantID", "SubjectId",
  "SubjectID", "PatientId", "PatientID", "match_pair_id", "subclass"
)

csv_files <- list.files(root, recursive = TRUE, full.names = TRUE, pattern = "\\.csv$", ignore.case = TRUE)
for (file in csv_files) {
  header <- names(utils::read.csv(file, nrows = 1, check.names = FALSE))
  bad <- intersect(header, sensitive_tokens)
  if (length(bad)) {
    stop("Direct identifier columns found in ", file, ": ", paste(bad, collapse = ", "), call. = FALSE)
  }
}

xlsx_files <- list.files(root, recursive = TRUE, full.names = TRUE, pattern = "\\.xlsx$", ignore.case = TRUE)
if (length(xlsx_files) > 0L) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package 'openxlsx' is required to audit XLSX publication candidates.", call. = FALSE)
  }
  for (file in xlsx_files) {
    for (sheet in openxlsx::getSheetNames(file)) {
      preview <- tryCatch(
        openxlsx::read.xlsx(file, sheet = sheet, rows = 1:12, colNames = FALSE, skipEmptyRows = FALSE),
        error = function(e) NULL
      )
      if (is.null(preview)) next
      values <- trimws(as.character(unlist(preview, use.names = FALSE)))
      bad <- intersect(values, sensitive_tokens)
      if (length(bad)) {
        stop("Direct identifier token found in ", basename(file), " / ", sheet,
             ": ", paste(unique(bad), collapse = ", "), call. = FALSE)
      }
      private_path_hits <- grep("([A-Za-z]:[/\\\\]Users[/\\\\]|/Users/|/home/[^/]+/)", values, value = TRUE, perl = TRUE)
      if (length(private_path_hits)) {
        stop("Local absolute path found in ", basename(file), " / ", sheet, call. = FALSE)
      }
    }
  }
}
message("Publication-candidate identifier and local-path audit passed.")

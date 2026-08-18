###############################################################################
# ReDLat plasma proteomics — machine-learning workflow
# 15. Restore and audit the submitted matched sensitivity cohort
# Requires: reviewer-derived ML metadata/master from RUN_PRECHECK
# Produces: exact pseudonymized 191 CN + 191 AD matched-ID set and balance audit
#
# Reproducibility note:
#   Nearest-neighbour propensity matching can choose different equally valid
#   controls when ties/order/software state differ. To reproduce the submitted
#   downstream matched analysis exactly, the deidentified metadata carries the
#   prespecified pseudonymous membership flag `include_matched_selected`.
#   This script audits and restores that exact submitted cohort; it does not
#   silently reselect a different 191/191 solution.
###############################################################################

.project_root <- if (nzchar(Sys.getenv("REDLAT_PROJECT_ROOT", unset = ""))) {
  normalizePath(Sys.getenv("REDLAT_PROJECT_ROOT"), winslash = "/", mustWork = TRUE)
} else if (requireNamespace("here", quietly = TRUE)) {
  normalizePath(here::here(), winslash = "/", mustWork = TRUE)
} else {
  stop("Package 'here' is required. Restore the R environment before running the workflow.", call. = FALSE)
}
source(file.path(.project_root, "R", "ml_bootstrap.R"), local = FALSE)
ML_CONFIG <- ml_load_config(.project_root)

required_pkgs <- c("dplyr", "readr", "tibble")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) {
  stop("Missing required R packages: ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
}
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

OUT_DIR <- file.path(ML_CONFIG$private_root, "Result_matching_rebuild_v5")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
ml_assert_private_path(OUT_DIR, ML_CONFIG)

as_flag <- function(x) {
  if (is.logical(x)) return(replace(x, is.na(x), FALSE))
  if (is.numeric(x)) return(!is.na(x) & x != 0)
  tolower(trimws(as.character(x))) %in% c("true", "1", "t", "yes", "y")
}

meta <- readr::read_csv(ML_CONFIG$metadata_file, show_col_types = FALSE, progress = FALSE)
required <- c(
  "SampleId", "SampleGroup", "Sex", "Age", "Education", "Country",
  "include_matched_selected"
)
missing <- setdiff(required, names(meta))
if (length(missing)) stop("ML metadata missing: ", paste(missing, collapse = ", "), call. = FALSE)
if (anyDuplicated(meta$SampleId)) stop("Duplicated SampleId in ML metadata.", call. = FALSE)

selected <- meta %>%
  dplyr::filter(as_flag(include_matched_selected)) %>%
  dplyr::mutate(
    SampleId = as.character(SampleId),
    SampleGroup = as.character(SampleGroup),
    Sex = as.character(Sex),
    Country = as.character(Country),
    Age = suppressWarnings(as.numeric(Age)),
    Education = suppressWarnings(as.numeric(Education))
  )

counts <- table(selected$SampleGroup)
if (nrow(selected) != 382L || as.integer(counts["CN"]) != 191L || as.integer(counts["AD"]) != 191L) {
  stop(
    "Submitted matched cohort must be exactly 382 participants (191 CN + 191 AD). Observed N=",
    nrow(selected), ", CN=", as.integer(counts["CN"]), ", AD=", as.integer(counts["AD"]),
    call. = FALSE
  )
}
if (anyDuplicated(selected$SampleId)) stop("Duplicated selected matched IDs.", call. = FALSE)
if (anyNA(selected[, c("SampleGroup", "Sex", "Age", "Education", "Country")])) {
  stop("Selected matched cohort contains missing required matching/model metadata.", call. = FALSE)
}

# Lightweight check that every selected pseudonym is present in the canonical ML master.
master_ids <- readr::read_csv(
  ML_CONFIG$master_file,
  col_select = c(SampleId, SampleGroup),
  show_col_types = FALSE,
  progress = FALSE
)
if (anyDuplicated(master_ids$SampleId)) stop("Duplicated SampleId in canonical ML master.", call. = FALSE)
missing_ids <- setdiff(selected$SampleId, master_ids$SampleId)
if (length(missing_ids)) stop("Selected matched IDs missing from canonical ML master.", call. = FALSE)
check_groups <- selected %>%
  dplyr::select(SampleId, SampleGroup_selected = SampleGroup) %>%
  dplyr::left_join(master_ids, by = "SampleId")
if (any(check_groups$SampleGroup_selected != check_groups$SampleGroup)) {
  stop("Diagnosis mismatch between selected matching flag and canonical ML master.", call. = FALSE)
}

selected_ids <- selected %>% dplyr::select(SampleId, SampleGroup)
readr::write_csv(selected_ids, file.path(OUT_DIR, "matched_ids_SELECTED.csv"))

balance <- selected %>%
  dplyr::group_by(SampleGroup) %>%
  dplyr::summarise(
    N = dplyr::n(),
    Age_mean = mean(Age),
    Age_sd = stats::sd(Age),
    Education_mean = mean(Education),
    Education_sd = stats::sd(Education),
    Female_n = sum(toupper(Sex) %in% c("F", "FEMALE")),
    Female_percent = 100 * mean(toupper(Sex) %in% c("F", "FEMALE")),
    .groups = "drop"
  )
readr::write_csv(balance, file.path(OUT_DIR, "matching_balance_SELECTED.csv"))

country_counts <- selected %>%
  dplyr::count(Country, SampleGroup, name = "N")
readr::write_csv(country_counts, file.path(OUT_DIR, "matching_country_counts_SELECTED.csv"))

provenance <- tibble::tibble(
  Item = c(
    "Selection source",
    "Selected total",
    "Selected CN",
    "Selected AD",
    "Participant identifier",
    "Selection policy"
  ),
  Value = c(
    "ReDLat_metadata_deidentified.csv::include_matched_selected",
    "382", "191", "191",
    "pseudonymous Study_ID/SampleId compatibility alias",
    "Restore exact submitted cohort; do not reselect an alternative nearest-neighbour solution"
  )
)
readr::write_csv(provenance, file.path(OUT_DIR, "matching_selection_provenance.csv"))

cat("\nMATCHED COHORT RESTORE/AUDIT PASSED\n")
cat("===================================\n")
cat("Selected matched output: 191 CN + 191 AD (N=382)\n")
cat("Source: deidentified metadata flag include_matched_selected\n")
cat("Outputs written to:\n", OUT_DIR, "\n", sep = "")

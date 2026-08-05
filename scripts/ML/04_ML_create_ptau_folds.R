###############################################################################
# ReDLat plasma proteomics — machine-learning workflow
# 04. Create p-tau217 nested folds
# Requires: private master matrix
# Produces: five stratified outer folds for the biomarker subset
# Data policy: participant-level inputs and predictions remain local.
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
EXCLUDED_SAMPLE_IDS <- ml_read_excluded_ids(ML_CONFIG)

project_root <- ML_CONFIG$project_root
csv_file <- ML_CONFIG$metadata_file
OUT_DIR <- file.path(ML_CONFIG$private_root, "ptau", "00_folds")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
ml_assert_private_path(OUT_DIR, ML_CONFIG)

library(here)
library(dplyr)
library(readr)
library(caret)

set.seed(42)

###############################################################################
# PATHS
###############################################################################

###############################################################################
# LOAD METADATA
###############################################################################

meta <- read.csv(
  csv_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

###############################################################################
# KEEP CN / AD ONLY
###############################################################################

ptau_candidates <- c("p-tau217", "p.tau217", "pTau217", "ptau217")
ptau_col <- ptau_candidates[ptau_candidates %in% names(meta)][1]
if (is.na(ptau_col)) stop("No p-tau217 column found in metadata.", call. = FALSE)
meta <- meta %>%
  filter(SampleGroup %in% c("CN", "AD")) %>%
  filter(!SampleId %in% EXCLUDED_SAMPLE_IDS) %>%
  filter(
    !is.na(Age),
    !is.na(Sex),
    !is.na(Country),
    !is.na(Education),
    !is.na(.data[[ptau_col]])
  )

cat("\n=====================\n")
cat("Subjects:", nrow(meta), "\n")
print(table(meta$SampleGroup))
cat("=====================\n")

readr::write_csv(tibble::tibble(Stage = "Fold-eligible cohort", N = nrow(meta), CN = sum(meta$SampleGroup == "CN"), AD = sum(meta$SampleGroup == "AD"), Excluded_IDs_configured = length(EXCLUDED_SAMPLE_IDS)), file.path(OUT_DIR, "cohort_flow.csv"))

###############################################################################
# CREATE 5 STRATIFIED FOLDS
###############################################################################

folds <- createFolds(
  y = meta$SampleGroup,
  k = 5,
  list = TRUE,
  returnTrain = FALSE
)

###############################################################################
# EXPORT
###############################################################################
all_test_ids <- c()

for (fold in seq_along(folds)) {
  test_idx <- folds[[fold]]
  test_ids <- meta$SampleId[test_idx]
  all_test_ids <- c(all_test_ids, test_ids)
  train_ids <- meta$SampleId[
    !(meta$SampleId %in% test_ids)
  ]
  overlap <- intersect(train_ids, test_ids)
  cat("\nFold", fold, "- overlap train/test:", length(overlap), "\n")
  write_csv(
    tibble(
      SampleId = train_ids
    ),
    file.path(
      OUT_DIR,
      paste0("fold_", fold, "_train_ids.csv")
    )
  )
  write_csv(
    tibble(
      SampleId = test_ids
    ),
    file.path(
      OUT_DIR,
      paste0("fold_", fold, "_test_ids.csv")
    )
  )
  cat("\nFold", fold, "\nTrain:", length(train_ids), "\nTest:", length(test_ids), "\n")
}

cat("\n====================================\n")

cat("Subjects in metadata:", nrow(meta), "\n")

cat("Unique test IDs:", length(unique(all_test_ids)), "\n")

cat("Total test IDs:", length(all_test_ids), "\n")

###############################################################################
# AUDIT
###############################################################################

audit <- bind_rows(
  lapply(seq_along(folds), function(fold){

    test_idx <- folds[[fold]]

    tmp <- meta[test_idx, ]

    tibble(
      Fold = fold,
      N = nrow(tmp),

      CN = sum(tmp$SampleGroup == "CN"),
      AD = sum(tmp$SampleGroup == "AD"),

      MeanAge = round(mean(tmp$Age),1),
      SDAge = round(sd(tmp$Age),1),

      MeanEducation = round(mean(tmp$Education),1),
      SDEducation = round(sd(tmp$Education),1),

      Male = sum(tmp$Sex == "2"),
      Female = sum(tmp$Sex == "1")
    )

  })
)

write_csv(
  audit,
  file.path(
    OUT_DIR,
    "fold_balance_audit.csv"
  )
)

###############################################################################
# COUNTRY AUDIT
###############################################################################

for (fold in seq_along(folds)){

  cat("\n====================================\n")
  cat("Fold", fold, "\n")
  print(table(meta$Country[folds[[fold]]]))
}

cat("\nDone.\n")
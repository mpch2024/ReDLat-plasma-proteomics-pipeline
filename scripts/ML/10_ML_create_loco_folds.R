###############################################################################
# ReDLat plasma proteomics — machine-learning workflow
# 10. Create country-held-out folds
# Requires: private master matrix
# Produces: one held-out fold per country
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
OUT_DIR <- file.path(ML_CONFIG$private_root, "loco", "00_folds")
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

table(meta$Country, meta$SampleGroup)

###############################################################################
# KEEP CN / AD ONLY
###############################################################################

meta <- meta %>%
  filter(
    SampleGroup %in% c("CN", "AD")
  ) %>%
  filter(
    !is.na(Age),
    !is.na(Sex),
    !is.na(Country),
    !is.na(Education)
  )

###############################################################################
# CREATE 5 FOLDS COUNTRIES
###############################################################################

countries <- sort(unique(meta$Country))

folds <- vector(
  "list",
  length(countries)
)

names(folds) <- countries

for(i in seq_along(countries)){

  test_country <- countries[i]

  test_idx <- which(
    meta$Country == test_country
  )

  train_idx <- which(
    meta$Country != test_country
  )

  folds[[i]] <- list(
    train = train_idx,
    test = test_idx,
    test_country = test_country
  )
}

all_test_ids <- c()

for (country in names(folds)) {

  test_idx  <- folds[[country]]$test
  train_idx <- folds[[country]]$train

  test_ids  <- meta$SampleId[test_idx]
  train_ids <- meta$SampleId[train_idx]
  all_test_ids <- c(all_test_ids, test_ids)

  overlap <- intersect(train_ids, test_ids)

  cat("\n", country, "- overlap:", length(overlap), "\n")

  write_csv(
    tibble(SampleId = train_ids),
    file.path(
      OUT_DIR,
      paste0(country, "_train_ids.csv")
    )
  )

  write_csv(
    tibble(SampleId = test_ids),
    file.path(
      OUT_DIR,
      paste0(country, "_test_ids.csv")
    )
  )

  cat(
    "\nCountry:",
    country,
    "\nTrain:",
    length(train_ids),
    "\nTest:",
    length(test_ids),
    "\n"
  )
}

audit <- bind_rows(
  lapply(folds, function(x) {

    train_meta <- meta[x$train,]
    test_meta  <- meta[x$test,]

    data.frame(
      TestCountry = x$test_country,

      Train_N = nrow(train_meta),
      Test_N  = nrow(test_meta),

      Train_AD = sum(train_meta$SampleGroup=="AD"),
      Train_CN = sum(train_meta$SampleGroup=="CN"),

      Test_AD = sum(test_meta$SampleGroup=="AD"),
      Test_CN = sum(test_meta$SampleGroup=="CN")
    )
  })
)

write_csv(audit, file.path(OUT_DIR, "LOCO_fold_audit.csv"))
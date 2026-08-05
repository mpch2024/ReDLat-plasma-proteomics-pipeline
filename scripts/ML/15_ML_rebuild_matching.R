###############################################################################
# ReDLat plasma proteomics — strict machine-learning workflow
# 15. Rebuild and audit propensity-score matching
# Requires: explicitly configured private master and metadata files
# Produces: selected matched IDs and aggregate balance diagnostics
# Data policy: participant-level matching files remain private.
###############################################################################

.project_root <- if (nzchar(Sys.getenv("REDLAT_PROJECT_ROOT", unset = ""))) {
  normalizePath(Sys.getenv("REDLAT_PROJECT_ROOT"), winslash = "/", mustWork = TRUE)
} else if (requireNamespace("here", quietly = TRUE)) {
  normalizePath(here::here(), winslash = "/", mustWork = TRUE)
} else stop("Package 'here' is required.", call. = FALSE)
source(file.path(.project_root, "R", "ml_bootstrap.R"), local = FALSE)
ML_CONFIG <- ml_load_config(.project_root)
ml_require_packages(c("MatchIt", "dplyr", "readr", "tibble", "stringr"))
suppressPackageStartupMessages({library(MatchIt); library(dplyr); library(readr); library(tibble); library(stringr)})
HAS_COBALT <- requireNamespace("cobalt", quietly = TRUE)

OUT_DIR <- file.path(ML_CONFIG$private_root, "matching")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
ml_assert_private_path(OUT_DIR, ML_CONFIG)
SEED <- 1111L; CALIPER <- 0.20

read_source <- function(path, label) {
  if (!file.exists(path)) stop(label, " not found: ", path, call. = FALSE)
  x <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE, name_repair = "minimal")
  candidates <- list(
    SampleId = c("SampleId", "sample_id", "ParticipantId", "ID"),
    SampleGroup = c("SampleGroup", "Diagnosis", "Group"),
    Sex = c("Sex", "sex", "Gender"), Age = c("Age", "age"),
    Education = c("Education", "education", "YearsEducation"), Country = c("Country", "country")
  )
  pick <- function(names) { hit <- names[names %in% colnames(x)]; if (length(hit)) hit[[1]] else NA_character_ }
  cols <- vapply(candidates, pick, character(1))
  if (any(is.na(cols[c("SampleId", "SampleGroup", "Sex", "Age", "Education")]))) stop(label, " lacks required matching columns.", call. = FALSE)
  clean_sex <- function(z) dplyr::case_when(toupper(trimws(as.character(z))) %in% c("1", "F", "FEMALE", "MUJER") ~ "F", toupper(trimws(as.character(z))) %in% c("2", "M", "MALE", "HOMBRE") ~ "M", TRUE ~ NA_character_)
  tibble(
    SampleId = trimws(as.character(x[[cols[["SampleId"]]]])),
    SampleGroup = trimws(as.character(x[[cols[["SampleGroup"]]]])),
    Sex = clean_sex(x[[cols[["Sex"]]]]),
    Age = suppressWarnings(as.numeric(x[[cols[["Age"]]]])),
    Education = suppressWarnings(as.numeric(x[[cols[["Education"]]]])),
    Country = if (is.na(cols[["Country"]])) NA_character_ else trimws(as.character(x[[cols[["Country"]]]]))
  ) %>% distinct(SampleId, .keep_all = TRUE)
}

current <- read_source(ML_CONFIG$master_file, "Configured master matrix") %>%
  filter(SampleGroup %in% c("CN", "AD"), !SampleId %in% ml_read_excluded_ids(ML_CONFIG))
complete <- current %>% filter(!is.na(Sex), !is.na(Age), !is.na(Education)) %>% mutate(SampleGroup_bin = as.integer(SampleGroup == "AD"), Sex = factor(Sex))
if (anyDuplicated(complete$SampleId)) stop("Duplicated SampleId values in matching input.", call. = FALSE)
set.seed(SEED)
fit <- MatchIt::matchit(SampleGroup_bin ~ Sex + Age + Education, data = complete, method = "nearest", distance = "logit", replace = FALSE, caliper = CALIPER)
matched <- MatchIt::match.data(fit) %>% as_tibble() %>% mutate(match_pair_id = if ("subclass" %in% names(.)) as.character(subclass) else NA_character_)
counts <- table(matched$SampleGroup)
summary <- tibble(
  input_total = nrow(complete), input_CN = sum(complete$SampleGroup == "CN"), input_AD = sum(complete$SampleGroup == "AD"),
  matched_total = nrow(matched), matched_CN = as.integer(counts[["CN"]]), matched_AD = as.integer(counts[["AD"]]),
  matched_pairs = n_distinct(na.omit(matched$match_pair_id)), exact_191_191 = matched_CN == 191L & matched_AD == 191L,
  seed = SEED, caliper = CALIPER, method = "nearest", distance = "logit", replace = FALSE,
  formula = "SampleGroup_bin ~ Sex + Age + Education", MatchIt_version = as.character(packageVersion("MatchIt")), R_version = R.version.string
)
readr::write_csv(summary, file.path(OUT_DIR, "matching_summary.csv"))
if (!isTRUE(summary$exact_191_191) && !ml_env_flag("REDLAT_ML_ALLOW_NONCANONICAL_MATCH_COUNTS", FALSE)) {
  stop("Matching did not reproduce 191 CN + 191 AD. Review matching_summary.csv before proceeding.", call. = FALSE)
}
readr::write_csv(matched, ML_CONFIG$matched_full_file)
readr::write_csv(matched %>% select(SampleId, SampleGroup, SampleGroup_bin, any_of(c("subclass", "match_pair_id", "distance", "weights", "Age", "Sex", "Education", "Country"))), ML_CONFIG$matched_ids_file)

balance <- bind_rows(
  complete %>% group_by(SampleGroup) %>% summarise(stage = "before", n = n(), age_mean = mean(Age), age_sd = sd(Age), education_mean = mean(Education), education_sd = sd(Education), female_percent = 100 * mean(Sex == "F"), .groups = "drop"),
  matched %>% group_by(SampleGroup) %>% summarise(stage = "after", n = n(), age_mean = mean(Age), age_sd = sd(Age), education_mean = mean(Education), education_sd = sd(Education), female_percent = 100 * mean(Sex == "F"), .groups = "drop")
) %>% select(stage, everything())
readr::write_csv(balance, file.path(OUT_DIR, "balance_summary.csv"))
if (HAS_COBALT) {
  cobalt_balance <- as.data.frame(cobalt::bal.tab(fit, un = TRUE, m.threshold = 0.1)$Balance) %>% tibble::rownames_to_column("Variable")
  readr::write_csv(cobalt_balance, file.path(OUT_DIR, "cobalt_balance.csv"))
}
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "matching_session_info.txt"))
message("Selected matching output: ", summary$matched_CN, " CN + ", summary$matched_AD, " AD.")

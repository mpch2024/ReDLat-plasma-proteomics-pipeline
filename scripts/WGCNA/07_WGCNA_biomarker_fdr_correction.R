###############################################################################
# ReDLat plasma proteomics — WGCNA workflow
# 07. Correct biomarker model multiplicity
# Requires: outputs from Scripts 05–06
# Produces: family-wide FDR and stability summaries
# Data policy: participant-level inputs and intermediate outputs remain local.
###############################################################################

rm(list = ls())


# -----------------------------------------------------------------------------
# Repository configuration
# -----------------------------------------------------------------------------
.project_root_env <- Sys.getenv("REDLAT_PROJECT_ROOT", unset = "")
if (nzchar(.project_root_env)) {
  project_root <- normalizePath(.project_root_env, winslash = "/", mustWork = TRUE)
} else if (requireNamespace("here", quietly = TRUE)) {
  project_root <- normalizePath(here::here(), winslash = "/", mustWork = TRUE)
} else {
  stop("Package 'here' is required. Restore the project environment with renv::restore().", call. = FALSE)
}
source(file.path(project_root, "R", "wgcna_bootstrap.R"), local = FALSE)
WGCNA_CONFIG <- wgcna_load_config(project_root)

###############################################################################
# 1) PACKAGES
###############################################################################

cran_pkgs <- c(
  "dplyr",
  "tidyr",
  "readr",
  "tibble",
  "stringr",
  "purrr",
  "ggplot2",
  "scales",
  "sandwich",
  "openxlsx",
  "forcats"
)

cran_missing <- cran_pkgs[
  !vapply(
    cran_pkgs,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(cran_missing) > 0L) {
  stop("Missing required packages: ", paste(cran_missing, collapse = ", "),
       ". Run renv::restore() before this script.", call. = FALSE)
}

invisible(lapply(
  cran_pkgs,
  library,
  character.only = TRUE
))

options(stringsAsFactors = FALSE)
options(error = traceback)

###############################################################################
# 2) PATHS
###############################################################################

BASE_DIR <- WGCNA_CONFIG$project_root

SCRIPT13B_DIR <- file.path(WGCNA_CONFIG$result_root,
  "05_sensitivity"
)

SCRIPT14_DIR <- file.path(WGCNA_CONFIG$result_root,
  "06_stability"
)

INPUT_FILE <- file.path(
  SCRIPT14_DIR,
  "tables",
  "sensitivity_input_clean.csv"
)

MODULE_REFERENCE_FILE <- file.path(WGCNA_CONFIG$result_root,
  "04_module_traits",
  "tables",
  "module_biological_label_reference.csv"
)

FULL_REFERENCE_32_FILE <- file.path(
  SCRIPT13B_DIR,
  "tables",
  "biomarker_robustness",
  "biomarker_models_log_HC3_primary_sensitivity.csv"
)

DOWNSAMPLING_MANIFEST_FILE <- file.path(
  SCRIPT14_DIR,
  "tables",
  "downsampling",
  "balanced_downsampling_sample_manifest.csv"
)

SITE_DELETION_MAP_FILE <- file.path(
  SCRIPT14_DIR,
  "tables",
  "within_country_loso",
  "country_site_deletion_map.csv"
)

OLD_FULL_12_FILE <- file.path(
  SCRIPT14_DIR,
  "tables",
  "full",
  "full_focal_log_biomarker_HC3_models_recomputed.csv"
)

OLD_LOCO_12_FILE <- file.path(
  SCRIPT14_DIR,
  "tables",
  "model_stability",
  "biomarker",
  "loco_focal_log_biomarker_HC3_results.csv"
)

OLD_WC_LOSO_12_FILE <- file.path(
  SCRIPT14_DIR,
  "tables",
  "model_stability",
  "biomarker",
  "within_country_loso_focal_log_biomarker_HC3_results.csv"
)

OLD_DOWNSAMPLING_12_FILE <- file.path(
  SCRIPT14_DIR,
  "tables",
  "model_stability",
  "biomarker",
  "balanced_downsampling_focal_log_biomarker_HC3_results.csv"
)

OUTDIR <- file.path(WGCNA_CONFIG$result_root,
  "07_biomarker_fdr"
)

SUBDIRS <- c(
  "tables",
  "tables/full",
  "tables/loco",
  "tables/within_country_loso",
  "tables/downsampling",
  "tables/focal_modules",
  "tables/fdr_family_audit",
  "figures",
  "figures/full",
  "figures/focal_modules",
  "figures/fdr_family_audit",
  "workspace"
)

invisible(lapply(
  file.path(OUTDIR, SUBDIRS),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

required_files <- c(
  INPUT_FILE,
  MODULE_REFERENCE_FILE,
  FULL_REFERENCE_32_FILE,
  DOWNSAMPLING_MANIFEST_FILE,
  SITE_DELETION_MAP_FILE,
  OLD_FULL_12_FILE,
  OLD_LOCO_12_FILE,
  OLD_WC_LOSO_12_FILE,
  OLD_DOWNSAMPLING_12_FILE
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  stop(
    "Missing required Script 13b/14 files:\n",
    paste(missing_files, collapse = "\n"),
    "\n\nRun Scripts 13b and 14 first.",
    call. = FALSE
  )
}

###############################################################################
# 3) PARAMETERS
###############################################################################

EXPECTED_N_SAMPLES <- 639L
EXPECTED_N_MODULES <- 8L
EXPECTED_N_BIOMARKERS <- 4L
EXPECTED_FULL_MODELS <- 32L
EXPECTED_N_COUNTRIES <- 5L
EXPECTED_WC_LOSO_UNITS <- 4L
EXPECTED_DOWNSAMPLING_ITERATIONS <- 500L
EXPECTED_DOWNSAMPLING_N_PER_ITERATION <- 170L
EXPECTED_DOWNSAMPLING_MANIFEST_ROWS <- 85000L

REQUIRE_EXPECTED_DIMENSIONS <- TRUE
REQUIRE_FULL_REFERENCE_MATCH <- TRUE

REFERENCE_TOLERANCE <- 1e-8
OLD_MODEL_ESTIMATE_TOLERANCE <- 1e-8

MIN_N_FOR_MODEL <- 30L
ROBUST_VCOV_TYPE <- "HC3"
PROGRESS_EVERY <- 25L

MODULES_OF_INTEREST <- NULL

FOCAL_MODULES <- c(
  "blue",
  "green",
  "brown"
)

BIOMARKERS <- c(
  "p_tau181",
  "p_tau217",
  "NfL",
  "ratio_AB42_40"
)

BIOMARKER_LABELS <- c(
  p_tau181 = "p-tau181",
  p_tau217 = "p-tau217",
  NfL = "NfL",
  ratio_AB42_40 = "Aβ42/40"
)

SITE_CANDIDATES <- c(
  "site",
  "Site",
  "Center",
  "center",
  "Cohort",
  "cohort",
  "RecruitmentSite",
  "recruitment_site",
  "site_id",
  "Site_ID"
)

MODULE_COLORS <- c(
  black       = "#2B2B2B",
  brown       = "#9C6B00",
  yellow      = "#B7D500",
  blue        = "#1535E8",
  green       = "#3E8F4E",
  red         = "#B4334A",
  purple      = "#8A2BE2",
  magenta     = "#C75ACD",
  pink        = "#F3A6D6",
  greenyellow = "#ADFF2F",
  turquoise   = "#40E0D0",
  cyan        = "#00BFC4",
  tan         = "#D2B48C",
  salmon      = "#FA8072",
  midnightblue = "#191970",
  lightcyan   = "#E0FFFF",
  grey        = "#9E9E9E"
)

DPI <- 300L

###############################################################################
# 4) HELPERS
###############################################################################

safe_read_csv <- function(file) {
  if (!file.exists(file)) {
    stop(
      "Missing required file:\n",
      file,
      call. = FALSE
    )
  }

  readr::read_csv(
    file,
    show_col_types = FALSE,
    guess_max = 100000
  )
}

safe_write_csv <- function(x, file) {
  dir.create(
    dirname(file),
    recursive = TRUE,
    showWarnings = FALSE
  )

  readr::write_csv(
    tibble::as_tibble(x),
    file
  )
}

clean_text_na <- function(x) {
  x <- as.character(x)
  x <- trimws(x)

  x[x %in% c(
    "",
    "NA",
    "NaN",
    "NULL",
    "null",
    "N/A",
    "nan"
  )] <- NA_character_

  x
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(x))
}

safe_mean <- function(x) {
  x <- safe_numeric(x)
  x <- x[is.finite(x)]

  if (length(x) == 0) {
    return(NA_real_)
  }

  mean(x)
}

safe_median <- function(x) {
  x <- safe_numeric(x)
  x <- x[is.finite(x)]

  if (length(x) == 0) {
    return(NA_real_)
  }

  stats::median(x)
}

safe_min <- function(x) {
  x <- safe_numeric(x)
  x <- x[is.finite(x)]

  if (length(x) == 0) {
    return(NA_real_)
  }

  min(x)
}

safe_max <- function(x) {
  x <- safe_numeric(x)
  x <- x[is.finite(x)]

  if (length(x) == 0) {
    return(NA_real_)
  }

  max(x)
}

safe_quantile <- function(x, prob) {
  x <- safe_numeric(x)
  x <- x[is.finite(x)]

  if (length(x) == 0) {
    return(NA_real_)
  }

  unname(
    stats::quantile(
      x,
      probs = prob,
      na.rm = TRUE
    )
  )
}

safe_proportion <- function(x) {
  x <- x[!is.na(x)]

  if (length(x) == 0) {
    return(NA_real_)
  }

  mean(
    as.logical(x)
  )
}

get_module_colors <- function(modules) {
  modules <- as.character(modules)
  cols <- MODULE_COLORS[modules]

  missing_modules <- modules[
    is.na(cols)
  ]

  if (length(missing_modules) > 0) {
    extra_cols <- grDevices::rainbow(
      length(unique(missing_modules)),
      s = 0.45,
      v = 0.75
    )

    names(extra_cols) <- unique(
      missing_modules
    )

    cols[is.na(cols)] <- extra_cols[
      missing_modules
    ]
  }

  names(cols) <- modules
  cols
}

detect_site_variable <- function(df) {
  hit <- SITE_CANDIDATES[
    SITE_CANDIDATES %in%
      names(df)
  ][1]

  if (
    length(hit) == 0 ||
    is.na(hit)
  ) {
    return(NA_character_)
  }

  hit
}

add_sig_stars <- function(q) {
  dplyr::case_when(
    is.na(q) ~ "",
    q < 0.001 ~ "***",
    q < 0.01 ~ "**",
    q < 0.05 ~ "*",
    TRUE ~ ""
  )
}

###############################################################################
# 5) LOG-HC3 MODEL FUNCTIONS
###############################################################################

extract_hc3_term <- function(
    fit,
    term,
    vcov_type = ROBUST_VCOV_TYPE
) {
  if (is.null(fit)) {
    return(tibble::tibble(
      estimate = NA_real_,
      robust_se = NA_real_,
      robust_statistic = NA_real_,
      robust_p_value = NA_real_,
      robust_conf_low = NA_real_,
      robust_conf_high = NA_real_,
      model_status = "model_failed"
    ))
  }

  coef_names <- names(
    stats::coef(fit)
  )

  if (!term %in% coef_names) {
    return(tibble::tibble(
      estimate = NA_real_,
      robust_se = NA_real_,
      robust_statistic = NA_real_,
      robust_p_value = NA_real_,
      robust_conf_low = NA_real_,
      robust_conf_high = NA_real_,
      model_status = "term_not_found"
    ))
  }

  robust_vcov <- tryCatch(
    sandwich::vcovHC(
      fit,
      type = vcov_type
    ),
    error = function(e) NULL
  )

  if (is.null(robust_vcov)) {
    return(tibble::tibble(
      estimate = stats::coef(
        fit
      )[[term]],
      robust_se = NA_real_,
      robust_statistic = NA_real_,
      robust_p_value = NA_real_,
      robust_conf_low = NA_real_,
      robust_conf_high = NA_real_,
      model_status = "robust_vcov_failed"
    ))
  }

  estimate <- stats::coef(
    fit
  )[[term]]

  robust_se <- sqrt(
    diag(
      robust_vcov
    )
  )[[term]]

  residual_df <- stats::df.residual(
    fit
  )

  robust_statistic <- estimate /
    robust_se

  robust_p_value <- 2 *
    stats::pt(
      -abs(
        robust_statistic
      ),
      df = residual_df
    )

  critical_value <- stats::qt(
    0.975,
    df = residual_df
  )

  tibble::tibble(
    estimate = estimate,
    robust_se = robust_se,
    robust_statistic =
      robust_statistic,
    robust_p_value =
      robust_p_value,
    robust_conf_low =
      estimate -
      critical_value *
      robust_se,
    robust_conf_high =
      estimate +
      critical_value *
      robust_se,
    model_status = "ok"
  )
}

fit_log_biomarker_hc3_single <- function(
    df,
    module,
    biomarker
) {
  required_vars <- c(
    module,
    biomarker,
    "Age",
    "Sex_bin",
    "Education",
    "Country"
  )

  if (!all(
    required_vars %in%
      names(df)
  )) {
    return(tibble::tibble(
      Module = module,
      Biomarker = biomarker,
      N = 0L,
      estimate_log_scale = NA_real_,
      robust_se = NA_real_,
      robust_statistic = NA_real_,
      robust_p_value = NA_real_,
      robust_conf_low = NA_real_,
      robust_conf_high = NA_real_,
      module_SD = NA_real_,
      log_outcome_SD = NA_real_,
      standardized_beta = NA_real_,
      percent_change_per_1SD_module =
        NA_real_,
      model_status =
        "required_variables_missing"
    ))
  }

  d <- df %>%
    dplyr::select(
      dplyr::all_of(
        required_vars
      )
    )

  d[[module]] <- safe_numeric(
    d[[module]]
  )

  d[[biomarker]] <- safe_numeric(
    d[[biomarker]]
  )

  d$.log_outcome <- ifelse(
    is.finite(
      d[[biomarker]]
    ) &
      d[[biomarker]] > 0,
    log(
      d[[biomarker]]
    ),
    NA_real_
  )

  d <- d %>%
    tidyr::drop_na(
      dplyr::all_of(
        c(
          module,
          ".log_outcome",
          "Age",
          "Sex_bin",
          "Education",
          "Country"
        )
      )
    )

  if (
    nrow(d) < MIN_N_FOR_MODEL ||
    dplyr::n_distinct(
      d[[module]]
    ) <= 1 ||
    dplyr::n_distinct(
      d$.log_outcome
    ) <= 1
  ) {
    return(tibble::tibble(
      Module = module,
      Biomarker = biomarker,
      N = nrow(d),
      estimate_log_scale = NA_real_,
      robust_se = NA_real_,
      robust_statistic = NA_real_,
      robust_p_value = NA_real_,
      robust_conf_low = NA_real_,
      robust_conf_high = NA_real_,
      module_SD = stats::sd(
        d[[module]],
        na.rm = TRUE
      ),
      log_outcome_SD = stats::sd(
        d$.log_outcome,
        na.rm = TRUE
      ),
      standardized_beta = NA_real_,
      percent_change_per_1SD_module =
        NA_real_,
      model_status =
        "insufficient_or_constant"
    ))
  }

  d$Country <- droplevels(
    factor(d$Country)
  )

  formula_text <- paste0(
    ".log_outcome ~ ",
    module,
    " + Age + Sex_bin + Education + Country"
  )

  fit <- tryCatch(
    stats::lm(
      stats::as.formula(
        formula_text
      ),
      data = d
    ),
    error = function(e) NULL
  )

  robust_term <- extract_hc3_term(
    fit = fit,
    term = module
  )

  module_sd <- stats::sd(
    d[[module]],
    na.rm = TRUE
  )

  log_outcome_sd <- stats::sd(
    d$.log_outcome,
    na.rm = TRUE
  )

  effect_per_1sd <- robust_term$estimate *
    module_sd

  standardized_beta <- if (
    is.finite(log_outcome_sd) &&
    log_outcome_sd > 0
  ) {
    effect_per_1sd /
      log_outcome_sd
  } else {
    NA_real_
  }

  percent_change_per_1sd <- if (
    is.finite(effect_per_1sd)
  ) {
    100 *
      (
        exp(
          effect_per_1sd
        ) - 1
      )
  } else {
    NA_real_
  }

  tibble::tibble(
    Module = module,
    Biomarker = biomarker,
    N = nrow(d),
    estimate_log_scale =
      robust_term$estimate,
    robust_se =
      robust_term$robust_se,
    robust_statistic =
      robust_term$robust_statistic,
    robust_p_value =
      robust_term$robust_p_value,
    robust_conf_low =
      robust_term$robust_conf_low,
    robust_conf_high =
      robust_term$robust_conf_high,
    module_SD = module_sd,
    log_outcome_SD =
      log_outcome_sd,
    standardized_beta =
      standardized_beta,
    percent_change_per_1SD_module =
      percent_change_per_1sd,
    model_status =
      robust_term$model_status
  )
}

run_log_biomarker_family32 <- function(
    df,
    modules,
    biomarkers,
    analysis,
    run_id,
    excluded_unit = NA_character_,
    iteration = NA_integer_
) {
  result <- purrr::map_dfr(
    modules,
    function(mod) {
      purrr::map_dfr(
        biomarkers,
        function(bio) {
          fit_log_biomarker_hc3_single(
            df = df,
            module = mod,
            biomarker = bio
          )
        }
      )
    }
  )

  expected_rows <- length(modules) *
    length(biomarkers)

  if (nrow(result) != expected_rows) {
    stop(
      "Expected ",
      expected_rows,
      " biomarker models in run ",
      run_id,
      ", but obtained ",
      nrow(result),
      ".",
      call. = FALSE
    )
  }

  result %>%
    dplyr::mutate(
      FDR_family32 = p.adjust(
        robust_p_value,
        method = "BH"
      ),
      FDR_family_size = expected_rows,
      Analysis = analysis,
      Run_ID = run_id,
      Excluded_unit =
        excluded_unit,
      Iteration = iteration,
      N_total_run = nrow(df)
    )
}

###############################################################################
# 6) STABILITY HELPERS
###############################################################################

attach_full_reference <- function(
    sensitivity_tbl,
    full_tbl
) {
  sensitivity_tbl %>%
    dplyr::left_join(
      full_tbl %>%
        dplyr::select(
          Module,
          Biomarker,
          standardized_beta_full =
            standardized_beta,
          robust_p_full =
            robust_p_value,
          FDR_full_family32 =
            FDR_family32
        ),
      by = c(
        "Module",
        "Biomarker"
      )
    ) %>%
    dplyr::mutate(
      delta_standardized_beta =
        standardized_beta -
        standardized_beta_full,
      abs_delta_standardized_beta =
        abs(
          delta_standardized_beta
        ),
      same_direction =
        dplyr::case_when(
          is.na(
            standardized_beta
          ) |
            is.na(
              standardized_beta_full
            ) ~ NA,
          standardized_beta == 0 |
            standardized_beta_full == 0 ~
            NA,
          TRUE ~
            sign(
              standardized_beta
            ) ==
              sign(
                standardized_beta_full
              )
        ),
      nominal_significant =
        !is.na(
          robust_p_value
        ) &
        robust_p_value < 0.05,
      FDR_significant_family32 =
        !is.na(
          FDR_family32
        ) &
        FDR_family32 < 0.05,
      full_FDR_significant_family32 =
        !is.na(
          FDR_full_family32
        ) &
        FDR_full_family32 < 0.05
    )
}

summarize_stability_by_model <- function(
    delta_tbl
) {
  summary_tbl <- delta_tbl %>%
    dplyr::group_by(
      Module,
      Biomarker
    ) %>%
    dplyr::summarise(
      standardized_beta_full =
        dplyr::first(
          standardized_beta_full
        ),
      robust_p_full =
        dplyr::first(
          robust_p_full
        ),
      FDR_full_family32 =
        dplyr::first(
          FDR_full_family32
        ),
      n_runs =
        dplyr::n_distinct(
          Run_ID
        ),
      mean_standardized_beta =
        safe_mean(
          standardized_beta
        ),
      median_standardized_beta =
        safe_median(
          standardized_beta
        ),
      empirical_beta_2_5 =
        safe_quantile(
          standardized_beta,
          0.025
        ),
      empirical_beta_97_5 =
        safe_quantile(
          standardized_beta,
          0.975
        ),
      mean_abs_delta_beta =
        safe_mean(
          abs_delta_standardized_beta
        ),
      median_abs_delta_beta =
        safe_median(
          abs_delta_standardized_beta
        ),
      max_abs_delta_beta =
        safe_max(
          abs_delta_standardized_beta
        ),
      direction_consistency =
        safe_proportion(
          same_direction
        ),
      nominal_significance_rate =
        safe_proportion(
          nominal_significant
        ),
      FDR32_significance_rate =
        safe_proportion(
          FDR_significant_family32
        ),
      full_significant_direction_preservation =
        dplyr::case_when(
          !dplyr::first(
            full_FDR_significant_family32
          ) ~ NA_real_,
          TRUE ~
            safe_proportion(
              same_direction
            )
        ),
      full_significant_direction_and_nominal_preservation =
        dplyr::case_when(
          !dplyr::first(
            full_FDR_significant_family32
          ) ~ NA_real_,
          TRUE ~
            safe_proportion(
              same_direction &
                nominal_significant
            )
        ),
      full_beta_inside_empirical_95_interval =
        dplyr::case_when(
          is.na(
            standardized_beta_full
          ) |
            is.na(
              empirical_beta_2_5
            ) |
            is.na(
              empirical_beta_97_5
            ) ~ NA,
          TRUE ~
            standardized_beta_full >=
              empirical_beta_2_5 &
            standardized_beta_full <=
              empirical_beta_97_5
        ),
      .groups = "drop"
    )

  worst_run <- delta_tbl %>%
    dplyr::group_by(
      Module,
      Biomarker
    ) %>%
    dplyr::slice_max(
      order_by =
        abs_delta_standardized_beta,
      n = 1,
      with_ties = FALSE
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(
      Module,
      Biomarker,
      worst_run_ID = Run_ID,
      worst_excluded_unit =
        Excluded_unit,
      worst_abs_delta_beta =
        abs_delta_standardized_beta
    )

  summary_tbl %>%
    dplyr::left_join(
      worst_run,
      by = c(
        "Module",
        "Biomarker"
      )
    )
}

summarize_stability_by_module <- function(
    model_summary
) {
  model_summary %>%
    dplyr::group_by(Module) %>%
    dplyr::summarise(
      n_biomarkers =
        dplyr::n(),
      n_full_FDR32_lt_0_05 =
        sum(
          FDR_full_family32 <
            0.05,
          na.rm = TRUE
        ),
      mean_direction_consistency =
        safe_mean(
          direction_consistency
        ),
      minimum_direction_consistency =
        safe_min(
          direction_consistency
        ),
      mean_FDR32_significance_rate =
        safe_mean(
          FDR32_significance_rate
        ),
      mean_abs_delta_beta =
        safe_mean(
          mean_abs_delta_beta
        ),
      maximum_abs_delta_beta =
        safe_max(
          max_abs_delta_beta
        ),
      .groups = "drop"
    )
}

###############################################################################
# 7) LOAD DATA
###############################################################################

analysis_df <- safe_read_csv(
  INPUT_FILE
)

module_reference <- safe_read_csv(
  MODULE_REFERENCE_FILE
)

full_reference_32 <- safe_read_csv(
  FULL_REFERENCE_32_FILE
)

downsampling_manifest <- safe_read_csv(
  DOWNSAMPLING_MANIFEST_FILE
)

site_deletion_map <- safe_read_csv(
  SITE_DELETION_MAP_FILE
)

old_full_12 <- safe_read_csv(
  OLD_FULL_12_FILE
)

old_loco_12 <- safe_read_csv(
  OLD_LOCO_12_FILE
)

old_wc_loso_12 <- safe_read_csv(
  OLD_WC_LOSO_12_FILE
)

old_downsampling_12 <- safe_read_csv(
  OLD_DOWNSAMPLING_12_FILE
)

required_input_cols <- c(
  "SampleId",
  "Country",
  "Age",
  "Sex_bin",
  "Education"
)

missing_input_cols <- setdiff(
  required_input_cols,
  names(analysis_df)
)

if (length(missing_input_cols) > 0) {
  stop(
    "The Script 14 clean input is missing required columns: ",
    paste(
      missing_input_cols,
      collapse = ", "
    ),
    call. = FALSE
  )
}

analysis_df <- analysis_df %>%
  dplyr::mutate(
    SampleId = as.character(
      SampleId
    ),
    Country = droplevels(
      factor(
        clean_text_na(
          Country
        )
      )
    ),
    Age = safe_numeric(Age),
    Sex_bin = safe_numeric(
      Sex_bin
    ),
    Education = safe_numeric(
      Education
    )
  ) %>%
  dplyr::distinct(
    SampleId,
    .keep_all = TRUE
  )

if ("Country_numeric" %in%
    names(analysis_df)) {
  analysis_df$Country_numeric <- NULL
}

site_var <- detect_site_variable(
  analysis_df
)

if (
  is.na(site_var) ||
  !site_var %in%
    names(analysis_df)
) {
  stop(
    "A recruitment-site variable was not detected.",
    call. = FALSE
  )
}

analysis_df[[site_var]] <- droplevels(
  factor(
    clean_text_na(
      analysis_df[[site_var]]
    )
  )
)

module_reference <- module_reference %>%
  dplyr::mutate(
    Module = as.character(Module)
  )

module_candidates <- module_reference$Module[
  module_reference$Module %in%
    names(analysis_df)
]

module_candidates <- module_candidates[
  vapply(
    analysis_df[module_candidates],
    is.numeric,
    logical(1)
  )
]

if (is.null(MODULES_OF_INTEREST)) {
  modules_use <- module_candidates
} else {
  modules_use <- intersect(
    MODULES_OF_INTEREST,
    module_candidates
  )
}

focal_modules_use <- intersect(
  FOCAL_MODULES,
  modules_use
)

if (length(modules_use) == 0) {
  stop(
    "No module eigengene columns were found.",
    call. = FALSE
  )
}

if (length(focal_modules_use) !=
    length(FOCAL_MODULES)) {
  stop(
    "Expected focal modules were not all detected: ",
    paste(
      setdiff(
        FOCAL_MODULES,
        focal_modules_use
      ),
      collapse = ", "
    ),
    call. = FALSE
  )
}

missing_biomarkers <- setdiff(
  BIOMARKERS,
  names(analysis_df)
)

if (length(missing_biomarkers) > 0) {
  stop(
    "Missing biomarker columns: ",
    paste(
      missing_biomarkers,
      collapse = ", "
    ),
    call. = FALSE
  )
}

for (mod in modules_use) {
  analysis_df[[mod]] <- safe_numeric(
    analysis_df[[mod]]
  )
}

for (bio in BIOMARKERS) {
  analysis_df[[bio]] <- safe_numeric(
    analysis_df[[bio]]
  )
}

countries <- sort(
  unique(
    as.character(
      analysis_df$Country
    )
  )
)

countries <- countries[
  !is.na(countries)
]

primary_wc_loso_units <- site_deletion_map %>%
  dplyr::filter(
    country_remains_after_site_exclusion
  ) %>%
  dplyr::mutate(
    Country = as.character(
      Country
    ),
    Site = as.character(Site),
    Country_Site = as.character(
      Country_Site
    )
  )

downsampling_manifest <- downsampling_manifest %>%
  dplyr::mutate(
    Iteration = as.integer(
      Iteration
    ),
    Run_ID = as.character(
      Run_ID
    ),
    SampleId = as.character(
      SampleId
    )
  )

###############################################################################
# 8) STRICT INPUT AUDITS
###############################################################################

manifest_iteration_summary <- downsampling_manifest %>%
  dplyr::group_by(
    Iteration,
    Run_ID
  ) %>%
  dplyr::summarise(
    N = dplyr::n(),
    n_unique_samples =
      dplyr::n_distinct(
        SampleId
      ),
    .groups = "drop"
  )

missing_manifest_samples <- setdiff(
  unique(
    downsampling_manifest$SampleId
  ),
  analysis_df$SampleId
)

if (length(missing_manifest_samples) > 0) {
  stop(
    "The downsampling manifest contains SampleId values absent from ",
    "the Script 14 clean input.",
    call. = FALSE
  )
}

if (REQUIRE_EXPECTED_DIMENSIONS) {
  if (nrow(analysis_df) !=
      EXPECTED_N_SAMPLES) {
    stop(
      "Expected ",
      EXPECTED_N_SAMPLES,
      " samples, but found ",
      nrow(analysis_df),
      ".",
      call. = FALSE
    )
  }

  if (length(modules_use) !=
      EXPECTED_N_MODULES) {
    stop(
      "Expected ",
      EXPECTED_N_MODULES,
      " modules, but found ",
      length(modules_use),
      ".",
      call. = FALSE
    )
  }

  if (length(BIOMARKERS) !=
      EXPECTED_N_BIOMARKERS) {
    stop(
      "Expected ",
      EXPECTED_N_BIOMARKERS,
      " biomarkers.",
      call. = FALSE
    )
  }

  if (length(countries) !=
      EXPECTED_N_COUNTRIES) {
    stop(
      "Expected ",
      EXPECTED_N_COUNTRIES,
      " countries, but found ",
      length(countries),
      ".",
      call. = FALSE
    )
  }

  if (nrow(primary_wc_loso_units) !=
      EXPECTED_WC_LOSO_UNITS) {
    stop(
      "Expected ",
      EXPECTED_WC_LOSO_UNITS,
      " valid WC-LOSO units, but found ",
      nrow(primary_wc_loso_units),
      ".",
      call. = FALSE
    )
  }

  if (
    dplyr::n_distinct(
      downsampling_manifest$Iteration
    ) !=
      EXPECTED_DOWNSAMPLING_ITERATIONS
  ) {
    stop(
      "Expected ",
      EXPECTED_DOWNSAMPLING_ITERATIONS,
      " downsampling iterations.",
      call. = FALSE
    )
  }

  if (nrow(downsampling_manifest) !=
      EXPECTED_DOWNSAMPLING_MANIFEST_ROWS) {
    stop(
      "Expected ",
      EXPECTED_DOWNSAMPLING_MANIFEST_ROWS,
      " downsampling-manifest rows, but found ",
      nrow(downsampling_manifest),
      ".",
      call. = FALSE
    )
  }

  if (any(
    manifest_iteration_summary$N !=
      EXPECTED_DOWNSAMPLING_N_PER_ITERATION
  )) {
    stop(
      "At least one downsampling iteration does not contain exactly ",
      EXPECTED_DOWNSAMPLING_N_PER_ITERATION,
      " participant rows.",
      call. = FALSE
    )
  }

  if (any(
    manifest_iteration_summary$
      n_unique_samples !=
      EXPECTED_DOWNSAMPLING_N_PER_ITERATION
  )) {
    stop(
      "At least one downsampling iteration contains duplicated participants.",
      call. = FALSE
    )
  }
}

input_audit <- tibble::tibble(
  metric = c(
    "base_dir",
    "script13b_dir",
    "script14_dir",
    "output_dir",
    "n_samples",
    "n_modules",
    "modules",
    "n_biomarkers",
    "biomarkers",
    "full_FDR_family_size",
    "n_countries_LOCO",
    "countries",
    "n_primary_WC_LOSO_units",
    "primary_WC_LOSO_units",
    "n_downsampling_iterations",
    "downsampling_n_per_iteration",
    "downsampling_manifest_rows",
    "focal_modules_extracted_after_FDR",
    "Country_numeric_present",
    "robust_vcov_type"
  ),
  value = c(
    BASE_DIR,
    SCRIPT13B_DIR,
    SCRIPT14_DIR,
    OUTDIR,
    as.character(
      nrow(analysis_df)
    ),
    as.character(
      length(modules_use)
    ),
    paste(
      modules_use,
      collapse = ", "
    ),
    as.character(
      length(BIOMARKERS)
    ),
    paste(
      BIOMARKERS,
      collapse = ", "
    ),
    as.character(
      length(modules_use) *
        length(BIOMARKERS)
    ),
    as.character(
      length(countries)
    ),
    paste(
      countries,
      collapse = ", "
    ),
    as.character(
      nrow(
        primary_wc_loso_units
      )
    ),
    paste(
      primary_wc_loso_units$
        Country_Site,
      collapse = ", "
    ),
    as.character(
      dplyr::n_distinct(
        downsampling_manifest$
          Iteration
      )
    ),
    as.character(
      unique(
        manifest_iteration_summary$N
      )
    ),
    as.character(
      nrow(
        downsampling_manifest
      )
    ),
    paste(
      focal_modules_use,
      collapse = ", "
    ),
    as.character(
      "Country_numeric" %in%
        names(analysis_df)
    ),
    ROBUST_VCOV_TYPE
  )
)

safe_write_csv(
  input_audit,
  file.path(
    OUTDIR,
    "tables",
    "script14b_input_audit.csv"
  )
)

safe_write_csv(
  manifest_iteration_summary,
  file.path(
    OUTDIR,
    "tables",
    "downsampling",
    "reused_downsampling_manifest_iteration_audit.csv"
  )
)

cat("Samples:", nrow(analysis_df), "\n")
cat("Modules:", paste(modules_use, collapse = ", "), "\n")
cat("Biomarkers:", paste(BIOMARKERS, collapse = ", "), "\n")
cat("FDR family per run:", length(modules_use) * length(BIOMARKERS), "\n")
cat("LOCO countries:", paste(countries, collapse = ", "), "\n")
cat(
  "WC-LOSO units:",
  paste(
    primary_wc_loso_units$Country_Site,
    collapse = ", "
  ),
  "\n"
)
cat(
  "Reused downsampling iterations:",
  dplyr::n_distinct(
    downsampling_manifest$Iteration
  ),
  "\n\n"
)

###############################################################################
# 9) FULL SAMPLE — 32-MODEL FAMILY AND SCRIPT 13b AUDIT
###############################################################################

full_family32 <- run_log_biomarker_family32(
  df = analysis_df,
  modules = modules_use,
  biomarkers = BIOMARKERS,
  analysis = "full",
  run_id = "full"
) %>%
  dplyr::mutate(
    Biomarker_label =
      dplyr::recode(
        Biomarker,
        !!!BIOMARKER_LABELS,
        .default = Biomarker
      )
  ) %>%
  dplyr::left_join(
    module_reference,
    by = "Module"
  )

if (nrow(full_family32) !=
    EXPECTED_FULL_MODELS) {
  stop(
    "The full 32-model family does not contain exactly ",
    EXPECTED_FULL_MODELS,
    " rows.",
    call. = FALSE
  )
}

full_reference_clean <- full_reference_32 %>%
  dplyr::transmute(
    Module = as.character(Module),
    Biomarker =
      as.character(Biomarker),
    N_reference =
      safe_numeric(N),
    estimate_reference =
      safe_numeric(estimate),
    robust_se_reference =
      safe_numeric(robust_se),
    robust_p_reference =
      safe_numeric(
        robust_p_value
      ),
    standardized_beta_reference =
      safe_numeric(
        standardized_beta
      ),
    FDR_reference_family32 =
      safe_numeric(
        FDR_within_transform_family
      )
  )

full_reference_audit <- full_family32 %>%
  dplyr::left_join(
    full_reference_clean,
    by = c(
      "Module",
      "Biomarker"
    )
  ) %>%
  dplyr::mutate(
    N_matches =
      N == N_reference,
    abs_delta_estimate =
      abs(
        estimate_log_scale -
          estimate_reference
      ),
    abs_delta_robust_se =
      abs(
        robust_se -
          robust_se_reference
      ),
    abs_delta_robust_p =
      abs(
        robust_p_value -
          robust_p_reference
      ),
    abs_delta_standardized_beta =
      abs(
        standardized_beta -
          standardized_beta_reference
      ),
    abs_delta_FDR =
      abs(
        FDR_family32 -
          FDR_reference_family32
      )
  )

full_reference_audit_summary <- tibble::tibble(
  metric = c(
    "n_models_recomputed",
    "n_reference_models",
    "n_models_matched",
    "all_N_match",
    "max_abs_delta_estimate",
    "max_abs_delta_robust_se",
    "max_abs_delta_robust_p",
    "max_abs_delta_standardized_beta",
    "max_abs_delta_FDR",
    "tolerance",
    "reference_match_passed"
  ),
  value = c(
    as.character(
      nrow(full_family32)
    ),
    as.character(
      nrow(full_reference_clean)
    ),
    as.character(
      sum(
        !is.na(
          full_reference_audit$
            estimate_reference
        )
      )
    ),
    as.character(
      all(
        full_reference_audit$
          N_matches,
        na.rm = TRUE
      )
    ),
    as.character(
      safe_max(
        full_reference_audit$
          abs_delta_estimate
      )
    ),
    as.character(
      safe_max(
        full_reference_audit$
          abs_delta_robust_se
      )
    ),
    as.character(
      safe_max(
        full_reference_audit$
          abs_delta_robust_p
      )
    ),
    as.character(
      safe_max(
        full_reference_audit$
          abs_delta_standardized_beta
      )
    ),
    as.character(
      safe_max(
        full_reference_audit$
          abs_delta_FDR
      )
    ),
    as.character(
      REFERENCE_TOLERANCE
    ),
    as.character(
      all(
        full_reference_audit$
          N_matches,
        na.rm = TRUE
      ) &&
        safe_max(
          full_reference_audit$
            abs_delta_estimate
        ) <=
          REFERENCE_TOLERANCE &&
        safe_max(
          full_reference_audit$
            abs_delta_robust_se
        ) <=
          REFERENCE_TOLERANCE &&
        safe_max(
          full_reference_audit$
            abs_delta_robust_p
        ) <=
          REFERENCE_TOLERANCE &&
        safe_max(
          full_reference_audit$
            abs_delta_standardized_beta
        ) <=
          REFERENCE_TOLERANCE &&
        safe_max(
          full_reference_audit$
            abs_delta_FDR
        ) <=
          REFERENCE_TOLERANCE
    )
  )
)

full_reference_match_passed <-
  all(
    full_reference_audit$N_matches,
    na.rm = TRUE
  ) &&
  safe_max(
    full_reference_audit$
      abs_delta_estimate
  ) <=
    REFERENCE_TOLERANCE &&
  safe_max(
    full_reference_audit$
      abs_delta_robust_se
  ) <=
    REFERENCE_TOLERANCE &&
  safe_max(
    full_reference_audit$
      abs_delta_robust_p
  ) <=
    REFERENCE_TOLERANCE &&
  safe_max(
    full_reference_audit$
      abs_delta_standardized_beta
  ) <=
    REFERENCE_TOLERANCE &&
  safe_max(
    full_reference_audit$
      abs_delta_FDR
  ) <=
    REFERENCE_TOLERANCE

safe_write_csv(
  full_family32,
  file.path(
    OUTDIR,
    "tables",
    "full",
    "full_all_modules_log_HC3_family32.csv"
  )
)

safe_write_csv(
  full_reference_audit,
  file.path(
    OUTDIR,
    "tables",
    "full",
    "full_family32_vs_script13b_reference_audit.csv"
  )
)

safe_write_csv(
  full_reference_audit_summary,
  file.path(
    OUTDIR,
    "tables",
    "full",
    "full_family32_vs_script13b_reference_audit_summary.csv"
  )
)

if (
  REQUIRE_FULL_REFERENCE_MATCH &&
  !full_reference_match_passed
) {
  stop(
    "The recomputed full 32-model family does not match Script 13b ",
    "within tolerance. Review the full reference-audit table.",
    call. = FALSE
  )
}

###############################################################################
# 10) LOCO — 32 MODELS AND FDR WITHIN EACH COUNTRY-DELETION RUN
###############################################################################

loco_results <- list()

for (cc in countries) {
  cat(
    "14b LOCO: excluding ",
    cc,
    "\n",
    sep = ""
  )

  df_run <- analysis_df %>%
    dplyr::filter(
      as.character(Country) !=
        cc
    )

  df_run$Country <- droplevels(
    factor(df_run$Country)
  )

  run_id <- paste0(
    "LOCO_",
    cc
  )

  loco_results[[cc]] <-
    run_log_biomarker_family32(
      df = df_run,
      modules = modules_use,
      biomarkers = BIOMARKERS,
      analysis = "LOCO_country",
      run_id = run_id,
      excluded_unit = cc
    )
}

loco_family32 <- dplyr::bind_rows(
  loco_results
) %>%
  dplyr::mutate(
    Biomarker_label =
      dplyr::recode(
        Biomarker,
        !!!BIOMARKER_LABELS,
        .default = Biomarker
      )
  ) %>%
  dplyr::left_join(
    module_reference,
    by = "Module"
  )

expected_loco_rows <- length(countries) *
  EXPECTED_FULL_MODELS

if (nrow(loco_family32) !=
    expected_loco_rows) {
  stop(
    "LOCO row-count mismatch: expected ",
    expected_loco_rows,
    ", observed ",
    nrow(loco_family32),
    ".",
    call. = FALSE
  )
}

loco_delta32 <- attach_full_reference(
  sensitivity_tbl =
    loco_family32,
  full_tbl =
    full_family32
)

loco_summary_by_model32 <-
  summarize_stability_by_model(
    loco_delta32
  ) %>%
  dplyr::mutate(
    Analysis =
      "LOCO country"
  ) %>%
  dplyr::left_join(
    module_reference,
    by = "Module"
  )

loco_summary_by_module32 <-
  summarize_stability_by_module(
    loco_summary_by_model32
  ) %>%
  dplyr::mutate(
    Analysis =
      "LOCO country"
  ) %>%
  dplyr::left_join(
    module_reference,
    by = "Module"
  )

safe_write_csv(
  loco_family32,
  file.path(
    OUTDIR,
    "tables",
    "loco",
    "loco_all_modules_log_HC3_family32_results.csv"
  )
)

safe_write_csv(
  loco_delta32,
  file.path(
    OUTDIR,
    "tables",
    "loco",
    "loco_all_modules_log_HC3_family32_delta_results.csv"
  )
)

safe_write_csv(
  loco_summary_by_model32,
  file.path(
    OUTDIR,
    "tables",
    "loco",
    "loco_family32_stability_by_module_biomarker.csv"
  )
)

safe_write_csv(
  loco_summary_by_module32,
  file.path(
    OUTDIR,
    "tables",
    "loco",
    "loco_family32_stability_by_module.csv"
  )
)

###############################################################################
# 11) CORRECTED WITHIN-COUNTRY LOSO — 32 MODELS PER RUN
###############################################################################

wc_loso_results <- list()

for (ii in seq_len(
  nrow(primary_wc_loso_units)
)) {
  country_i <-
    primary_wc_loso_units$
      Country[[ii]]

  site_i <-
    primary_wc_loso_units$
      Site[[ii]]

  unit_i <-
    primary_wc_loso_units$
      Country_Site[[ii]]

  cat(
    "14b WC-LOSO: excluding ",
    unit_i,
    "\n",
    sep = ""
  )

  df_run <- analysis_df %>%
    dplyr::filter(
      !(
        as.character(Country) ==
          country_i &
          as.character(
            .data[[site_var]]
          ) ==
            site_i
      )
    )

  df_run$Country <- droplevels(
    factor(df_run$Country)
  )

  df_run[[site_var]] <- droplevels(
    factor(
      df_run[[site_var]]
    )
  )

  run_id <- paste0(
    "WC_LOSO_",
    gsub(
      "::",
      "_",
      unit_i
    )
  )

  wc_loso_results[[unit_i]] <-
    run_log_biomarker_family32(
      df = df_run,
      modules = modules_use,
      biomarkers = BIOMARKERS,
      analysis =
        "WC_LOSO_primary",
      run_id = run_id,
      excluded_unit = unit_i
    )
}

wc_loso_family32 <- dplyr::bind_rows(
  wc_loso_results
) %>%
  dplyr::mutate(
    Biomarker_label =
      dplyr::recode(
        Biomarker,
        !!!BIOMARKER_LABELS,
        .default = Biomarker
      )
  ) %>%
  dplyr::left_join(
    module_reference,
    by = "Module"
  )

expected_wc_loso_rows <-
  nrow(primary_wc_loso_units) *
  EXPECTED_FULL_MODELS

if (nrow(wc_loso_family32) !=
    expected_wc_loso_rows) {
  stop(
    "WC-LOSO row-count mismatch: expected ",
    expected_wc_loso_rows,
    ", observed ",
    nrow(wc_loso_family32),
    ".",
    call. = FALSE
  )
}

wc_loso_delta32 <- attach_full_reference(
  sensitivity_tbl =
    wc_loso_family32,
  full_tbl =
    full_family32
)

wc_loso_summary_by_model32 <-
  summarize_stability_by_model(
    wc_loso_delta32
  ) %>%
  dplyr::mutate(
    Analysis =
      "WC-LOSO primary"
  ) %>%
  dplyr::left_join(
    module_reference,
    by = "Module"
  )

wc_loso_summary_by_module32 <-
  summarize_stability_by_module(
    wc_loso_summary_by_model32
  ) %>%
  dplyr::mutate(
    Analysis =
      "WC-LOSO primary"
  ) %>%
  dplyr::left_join(
    module_reference,
    by = "Module"
  )

safe_write_csv(
  wc_loso_family32,
  file.path(
    OUTDIR,
    "tables",
    "within_country_loso",
    "within_country_loso_all_modules_log_HC3_family32_results.csv"
  )
)

safe_write_csv(
  wc_loso_delta32,
  file.path(
    OUTDIR,
    "tables",
    "within_country_loso",
    "within_country_loso_all_modules_log_HC3_family32_delta_results.csv"
  )
)

safe_write_csv(
  wc_loso_summary_by_model32,
  file.path(
    OUTDIR,
    "tables",
    "within_country_loso",
    "within_country_loso_family32_stability_by_module_biomarker.csv"
  )
)

safe_write_csv(
  wc_loso_summary_by_module32,
  file.path(
    OUTDIR,
    "tables",
    "within_country_loso",
    "within_country_loso_family32_stability_by_module.csv"
  )
)

###############################################################################
# 12) EXACT REUSE OF THE 500 SCRIPT 14 DOWNSAMPLING SAMPLES
###############################################################################

downsampling_results <- vector(
  "list",
  EXPECTED_DOWNSAMPLING_ITERATIONS
)

iteration_values <- sort(
  unique(
    downsampling_manifest$Iteration
  )
)

for (jj in seq_along(
  iteration_values
)) {
  iter <- iteration_values[[jj]]

  manifest_iter <- downsampling_manifest %>%
    dplyr::filter(
      Iteration == iter
    )

  run_ids <- unique(
    manifest_iter$Run_ID
  )

  if (length(run_ids) != 1) {
    stop(
      "Iteration ",
      iter,
      " contains more than one Run_ID.",
      call. = FALSE
    )
  }

  run_id <- run_ids[[1]]

  sample_ids <- manifest_iter$SampleId

  if (anyDuplicated(sample_ids) > 0) {
    stop(
      "Iteration ",
      iter,
      " contains duplicated SampleId values.",
      call. = FALSE
    )
  }

  df_run <- analysis_df %>%
    dplyr::filter(
      SampleId %in%
        sample_ids
    )

  if (nrow(df_run) !=
      length(sample_ids)) {
    stop(
      "Iteration ",
      iter,
      " could not be reconstructed exactly from the manifest.",
      call. = FALSE
    )
  }

  df_run <- df_run[
    match(
      sample_ids,
      df_run$SampleId
    ),
    ,
    drop = FALSE
  ]

  if (!all(
    df_run$SampleId ==
      sample_ids
  )) {
    stop(
      "Iteration ",
      iter,
      " sample order could not be reconstructed.",
      call. = FALSE
    )
  }

  df_run$Country <- droplevels(
    factor(df_run$Country)
  )

  downsampling_results[[jj]] <-
    run_log_biomarker_family32(
      df = df_run,
      modules = modules_use,
      biomarkers = BIOMARKERS,
      analysis =
        "balanced_downsampling",
      run_id = run_id,
      iteration = iter
    )

  if (
    jj %% PROGRESS_EVERY == 0 ||
    jj ==
      length(iteration_values)
  ) {
    cat(
      "14b family32 downsampling completed: ",
      jj,
      "/",
      length(iteration_values),
      "\n",
      sep = ""
    )
  }
}

downsampling_family32 <- dplyr::bind_rows(
  downsampling_results
) %>%
  dplyr::mutate(
    Biomarker_label =
      dplyr::recode(
        Biomarker,
        !!!BIOMARKER_LABELS,
        .default = Biomarker
      )
  ) %>%
  dplyr::left_join(
    module_reference,
    by = "Module"
  )

expected_downsampling_rows <-
  EXPECTED_DOWNSAMPLING_ITERATIONS *
  EXPECTED_FULL_MODELS

if (nrow(downsampling_family32) !=
    expected_downsampling_rows) {
  stop(
    "Downsampling row-count mismatch: expected ",
    expected_downsampling_rows,
    ", observed ",
    nrow(downsampling_family32),
    ".",
    call. = FALSE
  )
}

downsampling_run_audit <- downsampling_family32 %>%
  dplyr::group_by(
    Iteration,
    Run_ID
  ) %>%
  dplyr::summarise(
    n_models = dplyr::n(),
    n_modules =
      dplyr::n_distinct(
        Module
      ),
    n_biomarkers =
      dplyr::n_distinct(
        Biomarker
      ),
    FDR_family_sizes =
      paste(
        sort(
          unique(
            FDR_family_size
          )
        ),
        collapse = ", "
      ),
    N_total_run =
      dplyr::first(
        N_total_run
      ),
    n_model_failures =
      sum(
        model_status != "ok",
        na.rm = TRUE
      ),
    .groups = "drop"
  )

if (any(
  downsampling_run_audit$
    n_models !=
    EXPECTED_FULL_MODELS
)) {
  stop(
    "At least one downsampling run does not contain exactly 32 models.",
    call. = FALSE
  )
}

if (any(
  downsampling_run_audit$
    N_total_run !=
    EXPECTED_DOWNSAMPLING_N_PER_ITERATION
)) {
  stop(
    "At least one reconstructed downsampling run does not contain 170 participants.",
    call. = FALSE
  )
}

downsampling_delta32 <-
  attach_full_reference(
    sensitivity_tbl =
      downsampling_family32,
    full_tbl =
      full_family32
  )

downsampling_summary_by_model32 <-
  summarize_stability_by_model(
    downsampling_delta32
  ) %>%
  dplyr::mutate(
    Analysis =
      "Balanced downsampling"
  ) %>%
  dplyr::left_join(
    module_reference,
    by = "Module"
  )

downsampling_summary_by_module32 <-
  summarize_stability_by_module(
    downsampling_summary_by_model32
  ) %>%
  dplyr::mutate(
    Analysis =
      "Balanced downsampling"
  ) %>%
  dplyr::left_join(
    module_reference,
    by = "Module"
  )

safe_write_csv(
  downsampling_family32,
  file.path(
    OUTDIR,
    "tables",
    "downsampling",
    "balanced_downsampling_all_modules_log_HC3_family32_results.csv"
  )
)

safe_write_csv(
  downsampling_delta32,
  file.path(
    OUTDIR,
    "tables",
    "downsampling",
    "balanced_downsampling_all_modules_log_HC3_family32_delta_results.csv"
  )
)

safe_write_csv(
  downsampling_summary_by_model32,
  file.path(
    OUTDIR,
    "tables",
    "downsampling",
    "balanced_downsampling_family32_stability_by_module_biomarker.csv"
  )
)

safe_write_csv(
  downsampling_summary_by_module32,
  file.path(
    OUTDIR,
    "tables",
    "downsampling",
    "balanced_downsampling_family32_stability_by_module.csv"
  )
)

safe_write_csv(
  downsampling_run_audit,
  file.path(
    OUTDIR,
    "tables",
    "downsampling",
    "balanced_downsampling_family32_run_audit.csv"
  )
)

###############################################################################
# 13) FOCAL MODULES EXTRACTED AFTER THE 32-MODEL CORRECTION
###############################################################################

full_focal32 <- full_family32 %>%
  dplyr::filter(
    Module %in%
      focal_modules_use
  )

loco_focal32 <- loco_delta32 %>%
  dplyr::filter(
    Module %in%
      focal_modules_use
  )

wc_loso_focal32 <- wc_loso_delta32 %>%
  dplyr::filter(
    Module %in%
      focal_modules_use
  )

downsampling_focal32 <-
  downsampling_delta32 %>%
  dplyr::filter(
    Module %in%
      focal_modules_use
  )

loco_focal_summary32 <-
  loco_summary_by_model32 %>%
  dplyr::filter(
    Module %in%
      focal_modules_use
  )

wc_loso_focal_summary32 <-
  wc_loso_summary_by_model32 %>%
  dplyr::filter(
    Module %in%
      focal_modules_use
  )

downsampling_focal_summary32 <-
  downsampling_summary_by_model32 %>%
  dplyr::filter(
    Module %in%
      focal_modules_use
  )

integrated_focal_stability32 <-
  full_focal32 %>%
  dplyr::select(
    Module,
    Biomarker,
    Biomarker_label,
    Module_color,
    Biological_label,
    full_standardized_beta =
      standardized_beta,
    full_percent_change_per_1SD =
      percent_change_per_1SD_module,
    full_robust_p =
      robust_p_value,
    full_FDR_family32 =
      FDR_family32
  ) %>%
  dplyr::left_join(
    loco_focal_summary32 %>%
      dplyr::select(
        Module,
        Biomarker,
        loco_direction_consistency =
          direction_consistency,
        loco_nominal_rate =
          nominal_significance_rate,
        loco_FDR32_rate =
          FDR32_significance_rate,
        loco_mean_abs_delta_beta =
          mean_abs_delta_beta,
        loco_max_abs_delta_beta =
          max_abs_delta_beta,
        loco_beta_2_5 =
          empirical_beta_2_5,
        loco_beta_97_5 =
          empirical_beta_97_5,
        loco_worst_excluded_unit =
          worst_excluded_unit
      ),
    by = c(
      "Module",
      "Biomarker"
    )
  ) %>%
  dplyr::left_join(
    wc_loso_focal_summary32 %>%
      dplyr::select(
        Module,
        Biomarker,
        wc_loso_direction_consistency =
          direction_consistency,
        wc_loso_nominal_rate =
          nominal_significance_rate,
        wc_loso_FDR32_rate =
          FDR32_significance_rate,
        wc_loso_mean_abs_delta_beta =
          mean_abs_delta_beta,
        wc_loso_max_abs_delta_beta =
          max_abs_delta_beta,
        wc_loso_beta_2_5 =
          empirical_beta_2_5,
        wc_loso_beta_97_5 =
          empirical_beta_97_5,
        wc_loso_worst_excluded_unit =
          worst_excluded_unit
      ),
    by = c(
      "Module",
      "Biomarker"
    )
  ) %>%
  dplyr::left_join(
    downsampling_focal_summary32 %>%
      dplyr::select(
        Module,
        Biomarker,
        downsampling_direction_consistency =
          direction_consistency,
        downsampling_nominal_rate =
          nominal_significance_rate,
        downsampling_FDR32_rate =
          FDR32_significance_rate,
        downsampling_mean_abs_delta_beta =
          mean_abs_delta_beta,
        downsampling_max_abs_delta_beta =
          max_abs_delta_beta,
        downsampling_beta_2_5 =
          empirical_beta_2_5,
        downsampling_beta_97_5 =
          empirical_beta_97_5
      ),
    by = c(
      "Module",
      "Biomarker"
    )
  ) %>%
  dplyr::arrange(
    Module,
    full_FDR_family32
  )

safe_write_csv(
  full_focal32,
  file.path(
    OUTDIR,
    "tables",
    "focal_modules",
    "full_focal_modules_after_family32_correction.csv"
  )
)

safe_write_csv(
  loco_focal32,
  file.path(
    OUTDIR,
    "tables",
    "focal_modules",
    "loco_focal_modules_after_family32_correction.csv"
  )
)

safe_write_csv(
  wc_loso_focal32,
  file.path(
    OUTDIR,
    "tables",
    "focal_modules",
    "within_country_loso_focal_modules_after_family32_correction.csv"
  )
)

safe_write_csv(
  downsampling_focal32,
  file.path(
    OUTDIR,
    "tables",
    "focal_modules",
    "balanced_downsampling_focal_modules_after_family32_correction.csv"
  )
)

safe_write_csv(
  integrated_focal_stability32,
  file.path(
    OUTDIR,
    "tables",
    "focal_modules",
    "integrated_focal_log_HC3_stability_family32.csv"
  )
)

###############################################################################
# 14) OLD 12-MODEL VS CORRECTED 32-MODEL FDR AUDIT
###############################################################################

prepare_old_family12 <- function(
    old_tbl,
    analysis_label
) {
  old_tbl %>%
    dplyr::transmute(
      Module = as.character(Module),
      Biomarker =
        as.character(Biomarker),
      Run_ID = as.character(
        Run_ID
      ),
      Iteration = suppressWarnings(
        as.integer(Iteration)
      ),
      Excluded_unit =
        as.character(
          Excluded_unit
        ),
      old_estimate_log_scale =
        safe_numeric(
          estimate_log_scale
        ),
      old_robust_p =
        safe_numeric(
          robust_p_value
        ),
      old_FDR_family12 =
        safe_numeric(
          FDR_HC3
        ),
      old_analysis =
        analysis_label
    )
}

prepare_new_focal32 <- function(
    new_tbl,
    analysis_label
) {
  new_tbl %>%
    dplyr::filter(
      Module %in%
        focal_modules_use
    ) %>%
    dplyr::transmute(
      Module = as.character(Module),
      Biomarker =
        as.character(Biomarker),
      Run_ID = as.character(
        Run_ID
      ),
      Iteration = suppressWarnings(
        as.integer(Iteration)
      ),
      Excluded_unit =
        as.character(
          Excluded_unit
        ),
      corrected_estimate_log_scale =
        safe_numeric(
          estimate_log_scale
        ),
      corrected_robust_p =
        safe_numeric(
          robust_p_value
        ),
      corrected_FDR_family32 =
        safe_numeric(
          FDR_family32
        ),
      corrected_analysis =
        analysis_label
    )
}

compare_family12_vs_family32 <- function(
    old_tbl,
    new_tbl,
    analysis_label
) {
  old_clean <- prepare_old_family12(
    old_tbl,
    analysis_label
  )

  new_clean <- prepare_new_focal32(
    new_tbl,
    analysis_label
  )

  join_keys <- c(
    "Module",
    "Biomarker",
    "Run_ID"
  )

  if (
    analysis_label ==
      "balanced_downsampling"
  ) {
    join_keys <- c(
      join_keys,
      "Iteration"
    )
  }

  old_clean %>%
    dplyr::full_join(
      new_clean,
      by = join_keys,
      suffix = c(
        "_old",
        "_new"
      )
    ) %>%
    dplyr::mutate(
      Analysis =
        analysis_label,
      abs_delta_estimate =
        abs(
          old_estimate_log_scale -
            corrected_estimate_log_scale
        ),
      abs_delta_robust_p =
        abs(
          old_robust_p -
            corrected_robust_p
        ),
      FDR_inflation =
        corrected_FDR_family32 -
        old_FDR_family12,
      FDR_ratio =
        corrected_FDR_family32 /
        old_FDR_family12,
      old_FDR12_lt_0_05 =
        !is.na(
          old_FDR_family12
        ) &
        old_FDR_family12 <
          0.05,
      corrected_FDR32_lt_0_05 =
        !is.na(
          corrected_FDR_family32
        ) &
        corrected_FDR_family32 <
          0.05,
      significance_class_change =
        dplyr::case_when(
          old_FDR12_lt_0_05 &
            corrected_FDR32_lt_0_05 ~
            "Significant in both",
          old_FDR12_lt_0_05 &
            !corrected_FDR32_lt_0_05 ~
            "Lost after family32 correction",
          !old_FDR12_lt_0_05 &
            corrected_FDR32_lt_0_05 ~
            "Gained after family32 correction",
          TRUE ~
            "Not significant in either"
        )
    )
}

full_family_audit <-
  compare_family12_vs_family32(
    old_tbl = old_full_12,
    new_tbl = full_family32,
    analysis_label = "full"
  )

loco_family_audit <-
  compare_family12_vs_family32(
    old_tbl = old_loco_12,
    new_tbl = loco_family32,
    analysis_label = "LOCO_country"
  )

wc_loso_family_audit <-
  compare_family12_vs_family32(
    old_tbl = old_wc_loso_12,
    new_tbl = wc_loso_family32,
    analysis_label =
      "WC_LOSO_primary"
  )

downsampling_family_audit <-
  compare_family12_vs_family32(
    old_tbl =
      old_downsampling_12,
    new_tbl =
      downsampling_family32,
    analysis_label =
      "balanced_downsampling"
  )

family_audit_all <- dplyr::bind_rows(
  full_family_audit,
  loco_family_audit,
  wc_loso_family_audit,
  downsampling_family_audit
)

if (
  safe_max(
    family_audit_all$
      abs_delta_estimate
  ) >
    OLD_MODEL_ESTIMATE_TOLERANCE ||
  safe_max(
    family_audit_all$
      abs_delta_robust_p
  ) >
    OLD_MODEL_ESTIMATE_TOLERANCE
) {
  stop(
    "Old focal and corrected family32 models differ in estimates or robust ",
    "P values. This correction should alter only the FDR family.",
    call. = FALSE
  )
}

family_audit_summary <- family_audit_all %>%
  dplyr::group_by(Analysis) %>%
  dplyr::summarise(
    n_models =
      dplyr::n(),
    max_abs_delta_estimate =
      safe_max(
        abs_delta_estimate
      ),
    max_abs_delta_robust_p =
      safe_max(
        abs_delta_robust_p
      ),
    mean_FDR_inflation =
      safe_mean(
        FDR_inflation
      ),
    median_FDR_inflation =
      safe_median(
        FDR_inflation
      ),
    maximum_FDR_inflation =
      safe_max(
        FDR_inflation
      ),
    n_significant_family12 =
      sum(
        old_FDR12_lt_0_05,
        na.rm = TRUE
      ),
    n_significant_family32 =
      sum(
        corrected_FDR32_lt_0_05,
        na.rm = TRUE
      ),
    n_lost_after_family32 =
      sum(
        significance_class_change ==
          "Lost after family32 correction",
        na.rm = TRUE
      ),
    n_gained_after_family32 =
      sum(
        significance_class_change ==
          "Gained after family32 correction",
        na.rm = TRUE
      ),
    .groups = "drop"
  )

family_audit_by_model <- family_audit_all %>%
  dplyr::group_by(
    Analysis,
    Module,
    Biomarker
  ) %>%
  dplyr::summarise(
    n_runs =
      dplyr::n(),
    mean_old_FDR12 =
      safe_mean(
        old_FDR_family12
      ),
    mean_corrected_FDR32 =
      safe_mean(
        corrected_FDR_family32
      ),
    mean_FDR_inflation =
      safe_mean(
        FDR_inflation
      ),
    family12_significance_rate =
      safe_proportion(
        old_FDR12_lt_0_05
      ),
    family32_significance_rate =
      safe_proportion(
        corrected_FDR32_lt_0_05
      ),
    significance_rate_change =
      family32_significance_rate -
      family12_significance_rate,
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Biomarker_label =
      dplyr::recode(
        Biomarker,
        !!!BIOMARKER_LABELS,
        .default = Biomarker
      )
  )

safe_write_csv(
  family_audit_all,
  file.path(
    OUTDIR,
    "tables",
    "fdr_family_audit",
    "old_family12_vs_corrected_family32_all_runs.csv"
  )
)

safe_write_csv(
  family_audit_summary,
  file.path(
    OUTDIR,
    "tables",
    "fdr_family_audit",
    "old_family12_vs_corrected_family32_summary.csv"
  )
)

safe_write_csv(
  family_audit_by_model,
  file.path(
    OUTDIR,
    "tables",
    "fdr_family_audit",
    "old_family12_vs_corrected_family32_by_model.csv"
  )
)

###############################################################################
# 15) FIGURES — FULL 32-MODEL FAMILY
###############################################################################

full_plot_df <- full_family32 %>%
  dplyr::mutate(
    Module = factor(
      Module,
      levels = rev(
        modules_use
      )
    ),
    Biomarker_label = factor(
      Biomarker_label,
      levels = unname(
        BIOMARKER_LABELS[
          BIOMARKERS
        ]
      )
    ),
    significance =
      add_sig_stars(
        FDR_family32
      ),
    cell_label =
      dplyr::case_when(
        is.na(
          standardized_beta
        ) ~ "",
        significance == "" ~
          sprintf(
            "%.2f",
            standardized_beta
          ),
        TRUE ~ paste0(
          sprintf(
            "%.2f",
            standardized_beta
          ),
          "\n",
          significance
        )
      )
  )

p_full_heatmap <- ggplot(
  full_plot_df,
  aes(
    x = Biomarker_label,
    y = Module,
    fill = standardized_beta
  )
) +
  geom_tile(
    colour = "white",
    linewidth = 0.5
  ) +
  geom_text(
    aes(
      label = cell_label
    ),
    size = 3.4
  ) +
  scale_fill_gradient2(
    low = "#4682B4",
    mid = "white",
    high = "#F46D43",
    midpoint = 0,
    name = "Standardized\nbeta",
    oob = scales::squish
  ) +
  labs(
    title =
      "Full log-HC3 module-biomarker associations",
    subtitle = paste0(
      "Stars indicate BH-FDR across all 32 module-biomarker models"
    ),
    x = NULL,
    y = NULL
  ) +
  theme_bw(
    base_size = 12
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 25,
      hjust = 1
    ),
    axis.text.y = element_text(
      face = "bold"
    ),
    plot.title = element_text(
      face = "bold"
    )
  )

ggsave(
  file.path(
    OUTDIR,
    "figures",
    "full",
    "full_log_HC3_family32_heatmap.pdf"
  ),
  p_full_heatmap,
  width = 7.5,
  height = 5.8
)

ggsave(
  file.path(
    OUTDIR,
    "figures",
    "full",
    "full_log_HC3_family32_heatmap.png"
  ),
  p_full_heatmap,
  width = 7.5,
  height = 5.8,
  dpi = DPI
)

###############################################################################
# 16) FIGURES — FOCAL STABILITY AFTER FAMILY32 CORRECTION
###############################################################################

focal_direction_plot <- integrated_focal_stability32 %>%
  dplyr::select(
    Module,
    Biomarker_label,
    LOCO =
      loco_direction_consistency,
    WC_LOSO =
      wc_loso_direction_consistency,
    Downsampling =
      downsampling_direction_consistency
  ) %>%
  tidyr::pivot_longer(
    cols = c(
      LOCO,
      WC_LOSO,
      Downsampling
    ),
    names_to = "Analysis",
    values_to =
      "Direction_consistency"
  ) %>%
  dplyr::mutate(
    Analysis = dplyr::recode(
      Analysis,
      LOCO = "LOCO country",
      WC_LOSO =
        "WC-LOSO primary",
      Downsampling =
        "Balanced downsampling"
    ),
    Analysis = factor(
      Analysis,
      levels = c(
        "LOCO country",
        "WC-LOSO primary",
        "Balanced downsampling"
      )
    ),
    Module_Biomarker = paste(
      Module,
      Biomarker_label,
      sep = " - "
    ),
    Module_Biomarker =
      factor(
        Module_Biomarker,
        levels = rev(
          unique(
            Module_Biomarker
          )
        )
      )
  )

p_direction <- ggplot(
  focal_direction_plot,
  aes(
    x = Analysis,
    y = Module_Biomarker,
    fill = Direction_consistency
  )
) +
  geom_tile(
    colour = "white",
    linewidth = 0.5
  ) +
  geom_text(
    aes(
      label = scales::percent(
        Direction_consistency,
        accuracy = 1
      )
    ),
    size = 3
  ) +
  scale_fill_gradient(
    low = "#FEE2E2",
    high = "#166534",
    limits = c(0, 1),
    oob = scales::squish,
    labels = scales::percent,
    name = "Direction\nconsistency"
  ) +
  labs(
    title =
      "Focal log-HC3 biomarker direction stability",
    subtitle = paste0(
      "Blue, green and brown extracted after the 32-model FDR correction"
    ),
    x = NULL,
    y = NULL
  ) +
  theme_bw(
    base_size = 11
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(
      size = 9
    ),
    plot.title = element_text(
      face = "bold"
    )
  )

ggsave(
  file.path(
    OUTDIR,
    "figures",
    "focal_modules",
    "focal_family32_direction_stability.pdf"
  ),
  p_direction,
  width = 8,
  height = 7
)

ggsave(
  file.path(
    OUTDIR,
    "figures",
    "focal_modules",
    "focal_family32_direction_stability.png"
  ),
  p_direction,
  width = 8,
  height = 7,
  dpi = DPI
)

focal_fdr_plot <- integrated_focal_stability32 %>%
  dplyr::select(
    Module,
    Biomarker_label,
    LOCO =
      loco_FDR32_rate,
    WC_LOSO =
      wc_loso_FDR32_rate,
    Downsampling =
      downsampling_FDR32_rate
  ) %>%
  tidyr::pivot_longer(
    cols = c(
      LOCO,
      WC_LOSO,
      Downsampling
    ),
    names_to = "Analysis",
    values_to =
      "FDR32_significance_rate"
  ) %>%
  dplyr::mutate(
    Analysis = dplyr::recode(
      Analysis,
      LOCO = "LOCO country",
      WC_LOSO =
        "WC-LOSO primary",
      Downsampling =
        "Balanced downsampling"
    ),
    Analysis = factor(
      Analysis,
      levels = c(
        "LOCO country",
        "WC-LOSO primary",
        "Balanced downsampling"
      )
    ),
    Module_Biomarker = paste(
      Module,
      Biomarker_label,
      sep = " - "
    ),
    Module_Biomarker =
      factor(
        Module_Biomarker,
        levels = rev(
          unique(
            Module_Biomarker
          )
        )
      )
  )

p_fdr_rate <- ggplot(
  focal_fdr_plot,
  aes(
    x = Analysis,
    y = Module_Biomarker,
    fill = FDR32_significance_rate
  )
) +
  geom_tile(
    colour = "white",
    linewidth = 0.5
  ) +
  geom_text(
    aes(
      label = scales::percent(
        FDR32_significance_rate,
        accuracy = 1
      )
    ),
    size = 3
  ) +
  scale_fill_gradient(
    low = "white",
    high = "#6A1B9A",
    limits = c(0, 1),
    oob = scales::squish,
    labels = scales::percent,
    name = "FDR < 0.05\nrun rate"
  ) +
  labs(
    title =
      "Focal biomarker FDR stability after family32 correction",
    subtitle = paste0(
      "Each run applies BH-FDR across all 8 modules x 4 biomarkers"
    ),
    x = NULL,
    y = NULL
  ) +
  theme_bw(
    base_size = 11
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(
      size = 9
    ),
    plot.title = element_text(
      face = "bold"
    )
  )

ggsave(
  file.path(
    OUTDIR,
    "figures",
    "focal_modules",
    "focal_family32_FDR_significance_stability.pdf"
  ),
  p_fdr_rate,
  width = 8,
  height = 7
)

ggsave(
  file.path(
    OUTDIR,
    "figures",
    "focal_modules",
    "focal_family32_FDR_significance_stability.png"
  ),
  p_fdr_rate,
  width = 8,
  height = 7,
  dpi = DPI
)

###############################################################################
# 17) FIGURE — FDR FAMILY CORRECTION AUDIT
###############################################################################

family_change_plot_df <- family_audit_by_model %>%
  dplyr::filter(
    Analysis %in% c(
      "LOCO_country",
      "WC_LOSO_primary",
      "balanced_downsampling"
    )
  ) %>%
  dplyr::mutate(
    Analysis = dplyr::recode(
      Analysis,
      LOCO_country =
        "LOCO country",
      WC_LOSO_primary =
        "WC-LOSO primary",
      balanced_downsampling =
        "Balanced downsampling"
    ),
    Module_Biomarker = paste(
      Module,
      Biomarker_label,
      sep = " - "
    ),
    Module_Biomarker =
      factor(
        Module_Biomarker,
        levels = rev(
          unique(
            Module_Biomarker
          )
        )
      )
  )

p_family_change <- ggplot(
  family_change_plot_df,
  aes(
    x = Analysis,
    y = Module_Biomarker,
    fill = significance_rate_change
  )
) +
  geom_tile(
    colour = "white",
    linewidth = 0.5
  ) +
  geom_text(
    aes(
      label = scales::percent(
        significance_rate_change,
        accuracy = 1
      )
    ),
    size = 3
  ) +
  scale_fill_gradient2(
    low = "#B2182B",
    mid = "white",
    high = "#2166AC",
    midpoint = 0,
    name = "Family32 minus\nfamily12 rate",
    labels = scales::percent
  ) +
  labs(
    title =
      "Impact of correcting the biomarker FDR family",
    subtitle = paste0(
      "Negative values indicate fewer significant runs after correction ",
      "from 12 to 32 models"
    ),
    x = NULL,
    y = NULL
  ) +
  theme_bw(
    base_size = 11
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(
      size = 9
    ),
    plot.title = element_text(
      face = "bold"
    )
  )

ggsave(
  file.path(
    OUTDIR,
    "figures",
    "fdr_family_audit",
    "family12_vs_family32_significance_rate_change.pdf"
  ),
  p_family_change,
  width = 8,
  height = 7
)

ggsave(
  file.path(
    OUTDIR,
    "figures",
    "fdr_family_audit",
    "family12_vs_family32_significance_rate_change.png"
  ),
  p_family_change,
  width = 8,
  height = 7,
  dpi = DPI
)

###############################################################################
# 18) METHODS-READY WORDING
###############################################################################

methods_wording <- tibble::tribble(
  ~section,
  ~text,

  "Methods - biomarker model family",
  paste(
    "For biomarker-model stability analyses, log-transformed biomarker",
    "outcomes were modeled against each of the eight module eigengenes,",
    "adjusting for age, sex, education and country as a categorical factor",
    "and using HC3 robust standard errors. In every full-sample,",
    "leave-one-country-out, within-country leave-one-site-out and balanced",
    "downsampling run, false-discovery-rate correction was applied across",
    "the complete family of 32 module-biomarker models",
    "(eight modules by four biomarkers)."
  ),

  "Methods - focal module display",
  paste(
    "The blue, green and brown modules were extracted for focused stability",
    "summaries only after model fitting and multiple-testing correction across",
    "the complete 32-model family."
  ),

  "Methods - downsampling reuse",
  paste(
    "The biomarker correction reused the exact participant manifests from",
    "the 500 country-by-diagnosis balanced samples generated in the primary",
    "association-stability workflow; no new random samples were drawn."
  ),

  "Interpretation",
  paste(
    "Direction consistency and effect-size preservation were interpreted",
    "jointly with the proportion of sensitivity runs surviving the complete",
    "32-model FDR correction. Loss of FDR significance after balanced",
    "downsampling was not treated as direction instability because the",
    "balanced samples intentionally reduced effective sample size."
  ),

  "Audit statement",
  paste(
    "The previous focal 12-model FDR values were retained only for audit.",
    "All manuscript-facing biomarker significance results should use the",
    "corrected 32-model family."
  )
)

safe_write_csv(
  methods_wording,
  file.path(
    OUTDIR,
    "tables",
    "script14b_methods_and_interpretation_wording.csv"
  )
)

###############################################################################
# 19) EXCEL WORKBOOK
###############################################################################

workbook_tables <- list(
  Input_audit =
    input_audit,
  Full_family32 =
    full_family32,
  Full_reference_audit =
    full_reference_audit_summary,
  LOCO_model_summary =
    loco_summary_by_model32,
  LOCO_module_summary =
    loco_summary_by_module32,
  WC_LOSO_model_summary =
    wc_loso_summary_by_model32,
  WC_LOSO_module_summary =
    wc_loso_summary_by_module32,
  Downsampling_model_summary =
    downsampling_summary_by_model32,
  Downsampling_module_summary =
    downsampling_summary_by_module32,
  Focal_integrated =
    integrated_focal_stability32,
  Family12_vs_32_summary =
    family_audit_summary,
  Family12_vs_32_by_model =
    family_audit_by_model,
  Methods_wording =
    methods_wording
)

openxlsx::write.xlsx(
  workbook_tables,
  file = file.path(
    OUTDIR,
    "WGCNA_14b_Biomarker_FDR32_Stability_Correction.xlsx"
  ),
  overwrite = TRUE
)

###############################################################################
# 20) FINAL SUMMARY AND OUTPUT MANIFEST
###############################################################################

script14b_summary <- tibble::tibble(
  metric = c(
    "base_dir",
    "output_dir",
    "n_samples_full",
    "n_modules",
    "modules",
    "n_biomarkers",
    "full_FDR_family_size",
    "full_reference_match_passed",
    "n_full_models",
    "n_LOCO_runs",
    "n_LOCO_models",
    "n_WC_LOSO_runs",
    "n_WC_LOSO_models",
    "n_downsampling_iterations",
    "n_downsampling_models",
    "n_downsampling_model_failures",
    "focal_modules_extracted_after_correction",
    "maximum_old_vs_new_estimate_difference",
    "maximum_old_vs_new_P_difference",
    "Country_numeric_included"
  ),
  value = c(
    BASE_DIR,
    OUTDIR,
    as.character(
      nrow(analysis_df)
    ),
    as.character(
      length(modules_use)
    ),
    paste(
      modules_use,
      collapse = ", "
    ),
    as.character(
      length(BIOMARKERS)
    ),
    as.character(
      EXPECTED_FULL_MODELS
    ),
    as.character(
      full_reference_match_passed
    ),
    as.character(
      nrow(full_family32)
    ),
    as.character(
      length(countries)
    ),
    as.character(
      nrow(loco_family32)
    ),
    as.character(
      nrow(
        primary_wc_loso_units
      )
    ),
    as.character(
      nrow(wc_loso_family32)
    ),
    as.character(
      length(iteration_values)
    ),
    as.character(
      nrow(
        downsampling_family32
      )
    ),
    as.character(
      sum(
        downsampling_family32$
          model_status != "ok",
        na.rm = TRUE
      )
    ),
    paste(
      focal_modules_use,
      collapse = ", "
    ),
    as.character(
      safe_max(
        family_audit_all$
          abs_delta_estimate
      )
    ),
    as.character(
      safe_max(
        family_audit_all$
          abs_delta_robust_p
      )
    ),
    as.character(
      "Country_numeric" %in%
        names(analysis_df)
    )
  )
)

safe_write_csv(
  script14b_summary,
  file.path(
    OUTDIR,
    "tables",
    "script14b_final_summary.csv"
  )
)

output_manifest <- tibble::tibble(
  output_file = c(
    "tables/script14b_input_audit.csv",
    "tables/full/full_all_modules_log_HC3_family32.csv",
    "tables/full/full_family32_vs_script13b_reference_audit.csv",
    "tables/full/full_family32_vs_script13b_reference_audit_summary.csv",
    "tables/loco/loco_all_modules_log_HC3_family32_results.csv",
    "tables/loco/loco_all_modules_log_HC3_family32_delta_results.csv",
    "tables/loco/loco_family32_stability_by_module_biomarker.csv",
    "tables/loco/loco_family32_stability_by_module.csv",
    "tables/within_country_loso/within_country_loso_all_modules_log_HC3_family32_results.csv",
    "tables/within_country_loso/within_country_loso_all_modules_log_HC3_family32_delta_results.csv",
    "tables/within_country_loso/within_country_loso_family32_stability_by_module_biomarker.csv",
    "tables/within_country_loso/within_country_loso_family32_stability_by_module.csv",
    "tables/downsampling/balanced_downsampling_all_modules_log_HC3_family32_results.csv",
    "tables/downsampling/balanced_downsampling_all_modules_log_HC3_family32_delta_results.csv",
    "tables/downsampling/balanced_downsampling_family32_stability_by_module_biomarker.csv",
    "tables/downsampling/balanced_downsampling_family32_stability_by_module.csv",
    "tables/downsampling/balanced_downsampling_family32_run_audit.csv",
    "tables/focal_modules/integrated_focal_log_HC3_stability_family32.csv",
    "tables/fdr_family_audit/old_family12_vs_corrected_family32_all_runs.csv",
    "tables/fdr_family_audit/old_family12_vs_corrected_family32_summary.csv",
    "tables/fdr_family_audit/old_family12_vs_corrected_family32_by_model.csv",
    "figures/full/full_log_HC3_family32_heatmap.pdf/png",
    "figures/focal_modules/focal_family32_direction_stability.pdf/png",
    "figures/focal_modules/focal_family32_FDR_significance_stability.pdf/png",
    "figures/fdr_family_audit/family12_vs_family32_significance_rate_change.pdf/png",
    "tables/script14b_methods_and_interpretation_wording.csv",
    "WGCNA_14b_Biomarker_FDR32_Stability_Correction.xlsx",
    "workspace/script14b_biomarker_FDR32_stability_workspace.RData",
    "sessionInfo.txt"
  ),
  description = c(
    "Strict sample, module, manifest and FDR-family audit.",
    "Full-sample 8-module by 4-biomarker log-HC3 family.",
    "Row-level comparison with the Script 13b 32-model reference.",
    "Summary of exact full-sample reference agreement.",
    "All 32 models for each leave-one-country-out run.",
    "LOCO effects joined to the full 32-model reference.",
    "LOCO stability by module-biomarker pair.",
    "LOCO stability summarized by module.",
    "All 32 models for each valid within-country site deletion.",
    "WC-LOSO effects joined to the full reference.",
    "WC-LOSO stability by module-biomarker pair.",
    "WC-LOSO stability summarized by module.",
    "All 16,000 models from the exact 500 balanced samples.",
    "Downsampling effects joined to the full reference.",
    "Empirical intervals and stability by module-biomarker.",
    "Downsampling stability summarized by module.",
    "Per-iteration verification of 32 models and 170 participants.",
    "Blue, green and brown extracted after family32 correction.",
    "Complete audit of old 12-model versus corrected 32-model FDR values.",
    "Family-level count of significance changes.",
    "Model-level change in significance rate.",
    "Full 32-model standardized-beta heatmap.",
    "Focal direction stability after correct FDR handling.",
    "Focal FDR-significance run rates using family32.",
    "Impact of changing the FDR family from 12 to 32 models.",
    "Methods-ready wording and interpretation boundaries.",
    "Integrated reviewer-facing workbook.",
    "Workspace for final figure and table generation.",
    "R session information."
  )
)

safe_write_csv(
  output_manifest,
  file.path(
    OUTDIR,
    "script14b_output_manifest.csv"
  )
)

###############################################################################
# 21) SAVE WORKSPACE
###############################################################################

save(
  analysis_df,
  module_reference,
  modules_use,
  focal_modules_use,
  BIOMARKERS,
  countries,
  site_var,
  primary_wc_loso_units,
  downsampling_manifest,
  manifest_iteration_summary,
  full_family32,
  full_reference_audit,
  full_reference_audit_summary,
  loco_family32,
  loco_delta32,
  loco_summary_by_model32,
  loco_summary_by_module32,
  wc_loso_family32,
  wc_loso_delta32,
  wc_loso_summary_by_model32,
  wc_loso_summary_by_module32,
  downsampling_family32,
  downsampling_delta32,
  downsampling_summary_by_model32,
  downsampling_summary_by_module32,
  downsampling_run_audit,
  full_focal32,
  loco_focal32,
  wc_loso_focal32,
  downsampling_focal32,
  integrated_focal_stability32,
  family_audit_all,
  family_audit_summary,
  family_audit_by_model,
  methods_wording,
  script14b_summary,
  output_manifest,
  file = file.path(
    OUTDIR,
    "workspace",
    "script14b_biomarker_FDR32_stability_workspace.RData"
  )
)

writeLines(
  capture.output(
    utils::sessionInfo()
  ),
  con = file.path(
    OUTDIR,
    "sessionInfo.txt"
  )
)

###############################################################################
# 22) FINAL MESSAGE
###############################################################################

cat("\nScript 14b finished successfully.\n")
cat("Main output directory:\n", OUTDIR, "\n")
cat("Samples in full analysis: ", nrow(analysis_df), "\n", sep = "")
cat("Modules modeled per run: ", length(modules_use), "\n", sep = "")
cat("Biomarkers modeled per run: ", length(BIOMARKERS), "\n", sep = "")
cat("Corrected FDR family size: ", EXPECTED_FULL_MODELS, "\n", sep = "")
cat("Full Script 13b reference match: ", full_reference_match_passed, "\n", sep = "")
cat("LOCO models completed: ", nrow(loco_family32), "\n", sep = "")
cat("WC-LOSO models completed: ", nrow(wc_loso_family32), "\n", sep = "")
cat("Downsampling iterations reused: ", length(iteration_values), "\n", sep = "")
cat("Downsampling models completed: ", nrow(downsampling_family32), "\n", sep = "")
cat(
  "Downsampling model failures: ",
  sum(
    downsampling_family32$model_status != "ok",
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)
cat(
  "Maximum old-vs-new estimate difference: ",
  safe_max(
    family_audit_all$abs_delta_estimate
  ),
  "\n",
  sep = ""
)
cat(
  "Maximum old-vs-new robust P difference: ",
  safe_max(
    family_audit_all$abs_delta_robust_p
  ),
  "\n",
  sep = ""
)
cat("Country_numeric included: FALSE\n")
cat("\nKey outputs:\n")
cat("- tables/full/full_all_modules_log_HC3_family32.csv\n")
cat("- tables/loco/loco_family32_stability_by_module_biomarker.csv\n")
cat("- tables/within_country_loso/within_country_loso_family32_stability_by_module_biomarker.csv\n")
cat("- tables/downsampling/balanced_downsampling_family32_stability_by_module_biomarker.csv\n")
cat("- tables/focal_modules/integrated_focal_log_HC3_stability_family32.csv\n")
cat("- tables/fdr_family_audit/old_family12_vs_corrected_family32_summary.csv\n")
cat("- WGCNA_14b_Biomarker_FDR32_Stability_Correction.xlsx\n")

###############################################################################
# END
###############################################################################


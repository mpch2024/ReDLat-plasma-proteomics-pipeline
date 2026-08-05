###############################################################################
# ReDLat plasma proteomics — WGCNA workflow
# 06. Assess association stability
# Requires: outputs from Scripts 04–05
# Produces: LOCO, LOSO and balanced-downsampling results
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
  "lmtest",
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

SCRIPT13_DIR <- file.path(WGCNA_CONFIG$result_root,
  "04_module_traits"
)

SCRIPT13B_DIR <- file.path(WGCNA_CONFIG$result_root,
  "05_sensitivity"
)

INPUT_FILE <- file.path(
  SCRIPT13_DIR,
  "tables",
  "module_trait_input_clean.csv"
)

MODULE_REFERENCE_FILE <- file.path(
  SCRIPT13_DIR,
  "tables",
  "module_biological_label_reference.csv"
)

FULL_CORRELATION_REFERENCE_FILE <- file.path(
  SCRIPT13_DIR,
  "tables",
  "correlations",
  "module_trait_results_long.csv"
)

PRIORITIZATION_FILE <- file.path(
  SCRIPT13_DIR,
  "tables",
  "prioritization",
  "final_module_prioritization_table.csv"
)

DIAGNOSIS_REFERENCE_FILE <- file.path(
  SCRIPT13B_DIR,
  "tables",
  "diagnosis_robustness",
  "adjusted_diagnosis_module_models_HC3.csv"
)

BIOMARKER_REFERENCE_FILE <- file.path(
  SCRIPT13B_DIR,
  "tables",
  "biomarker_robustness",
  "biomarker_models_log_HC3_primary_sensitivity.csv"
)

SITE_STRUCTURE_REFERENCE_FILE <- file.path(
  SCRIPT13B_DIR,
  "tables",
  "context",
  "site_nested_structure_audit.csv"
)

OUTDIR <- file.path(WGCNA_CONFIG$result_root,
  "06_stability"
)

SUBDIRS <- c(
  "tables",
  "tables/full",
  "tables/loco",
  "tables/within_country_loso",
  "tables/all_site_deletion_audit",
  "tables/downsampling",
  "tables/model_stability",
  "tables/model_stability/diagnosis",
  "tables/model_stability/biomarker",
  "tables/core_modules",
  "figures",
  "figures/loco",
  "figures/within_country_loso",
  "figures/downsampling",
  "figures/model_stability",
  "figures/core_modules",
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
  FULL_CORRELATION_REFERENCE_FILE,
  PRIORITIZATION_FILE,
  DIAGNOSIS_REFERENCE_FILE,
  BIOMARKER_REFERENCE_FILE,
  SITE_STRUCTURE_REFERENCE_FILE
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  stop(
    "Missing required Script 13/13b files:\n",
    paste(missing_files, collapse = "\n"),
    "\n\nRun Scripts 13 and 13b first.",
    call. = FALSE
  )
}

###############################################################################
# 3) PARAMETERS
###############################################################################

MIN_N_FOR_CORR <- 6L
MIN_N_FOR_MODEL <- 30L

EXPECTED_N_SAMPLES <- 639L
EXPECTED_N_MODULES <- 8L
EXPECTED_N_TRAITS <- 16L
EXPECTED_N_COUNTRIES <- 5L
EXPECTED_N_SITE_LEVELS <- 7L
EXPECTED_N_COUNTRY_DIAGNOSIS_STRATA <- 10L

REQUIRE_EXPECTED_DIMENSIONS <- TRUE
REQUIRE_FULL_CORRELATION_MATCH <- TRUE
FULL_CORRELATION_TOLERANCE <- 1e-9
MODEL_REFERENCE_TOLERANCE <- 1e-6

N_DOWNSAMPLING_ITER <- 500L
SEED <- 123L
PROGRESS_EVERY <- 50L

ROBUST_VCOV_TYPE <- "HC3"

EXCLUDE_MODULES <- c(
  "grey",
  "MEgrey"
)

CORE_MODULES <- c(
  "blue",
  "green",
  "brown"
)

FOCAL_BIOMARKER_MODULES <- c(
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

TRAITS_USE <- c(
  "SampleGroup_bin",
  "cdr_global",
  "cdr_boxscore",
  "mmse_total",
  "udsfaq_total",
  "NPI",
  "Mini_SEA",
  "T_ADLQ",
  "p_tau181",
  "p_tau217",
  "NfL",
  "ratio_AB42_40",
  "Age",
  "Sex_bin",
  "Education",
  "APOE4_carrier"
)

TRAIT_LABELS <- c(
  SampleGroup_bin = "Diagnosis",
  cdr_global = "CDR global",
  cdr_boxscore = "CDR-SB",
  mmse_total = "MMSE",
  udsfaq_total = "PFAQ",
  NPI = "NPI-Q",
  Mini_SEA = "Mini-SEA",
  T_ADLQ = "T-ADLQ",
  p_tau181 = "p-tau181",
  p_tau217 = "p-tau217",
  NfL = "NfL",
  ratio_AB42_40 = "Aβ42/40",
  Age = "Age",
  Sex_bin = "Sex",
  Education = "Education",
  APOE4_carrier = "APOE ε4"
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
# 4) BASIC HELPERS
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

safe_proportion <- function(x) {
  x <- x[!is.na(x)]

  if (length(x) == 0) {
    return(NA_real_)
  }

  mean(
    as.logical(x)
  )
}

###############################################################################
# 5) CORRELATION HELPERS
###############################################################################

safe_spearman <- function(
    df,
    x,
    y,
    min_n = MIN_N_FOR_CORR
) {
  d <- df %>%
    dplyr::select(
      dplyr::all_of(c(
        x,
        y
      ))
    ) %>%
    tidyr::drop_na()

  n <- nrow(d)

  if (n < min_n) {
    return(tibble::tibble(
      rho = NA_real_,
      p_value = NA_real_,
      N = n
    ))
  }

  if (
    dplyr::n_distinct(
      d[[x]]
    ) <= 1 ||
    dplyr::n_distinct(
      d[[y]]
    ) <= 1
  ) {
    return(tibble::tibble(
      rho = NA_real_,
      p_value = NA_real_,
      N = n
    ))
  }

  test <- suppressWarnings(
    stats::cor.test(
      d[[x]],
      d[[y]],
      method = "spearman",
      exact = FALSE
    )
  )

  tibble::tibble(
    rho = unname(
      test$estimate
    ),
    p_value = test$p.value,
    N = n
  )
}

compute_module_trait <- function(
    df,
    modules,
    traits,
    analysis,
    run_id,
    excluded_unit = NA_character_,
    iteration = NA_integer_
) {
  out <- expand.grid(
    Module = modules,
    Trait = traits,
    stringsAsFactors = FALSE
  ) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
      result = purrr::map2(
        Module,
        Trait,
        ~ safe_spearman(
          df,
          .x,
          .y,
          min_n = MIN_N_FOR_CORR
        )
      )
    ) %>%
    tidyr::unnest(result) %>%
    dplyr::mutate(
      FDR = p.adjust(
        p_value,
        method = "BH"
      ),
      Analysis = analysis,
      Run_ID = run_id,
      Excluded_unit =
        excluded_unit,
      Iteration = iteration,
      N_total_run = nrow(df)
    )

  out
}

attach_full_correlation_reference <- function(
    sensitivity_tbl,
    full_tbl
) {
  sensitivity_tbl %>%
    dplyr::left_join(
      full_tbl %>%
        dplyr::select(
          Module,
          Trait,
          rho_full = rho,
          p_full = p_value,
          FDR_full = FDR,
          N_full = N
        ),
      by = c(
        "Module",
        "Trait"
      )
    ) %>%
    dplyr::mutate(
      delta_rho = rho - rho_full,
      abs_delta_rho = abs(
        delta_rho
      ),
      abs_rho_full = abs(
        rho_full
      ),
      abs_rho_sensitivity = abs(
        rho
      ),
      same_direction = dplyr::case_when(
        is.na(rho) |
          is.na(rho_full) ~ NA,
        rho == 0 |
          rho_full == 0 ~ NA,
        TRUE ~
          sign(rho) ==
            sign(rho_full)
      ),
      full_FDR_lt_0_05 =
        !is.na(FDR_full) &
        FDR_full < 0.05,
      sensitivity_nominal_lt_0_05 =
        !is.na(p_value) &
        p_value < 0.05,
      sensitivity_FDR_lt_0_05 =
        !is.na(FDR) &
        FDR < 0.05,
      preserved_direction_and_nominal_if_full_sig =
        dplyr::case_when(
          !full_FDR_lt_0_05 ~ NA,
          is.na(same_direction) ~ NA,
          TRUE ~
            same_direction &
            sensitivity_nominal_lt_0_05
        ),
      preserved_direction_if_full_sig =
        dplyr::case_when(
          !full_FDR_lt_0_05 ~ NA,
          TRUE ~
            same_direction
        )
    )
}

summarize_correlation_delta_by_module <- function(
    delta_tbl
) {
  delta_tbl %>%
    dplyr::group_by(Module) %>%
    dplyr::summarise(
      n_runs =
        dplyr::n_distinct(
          Run_ID
        ),
      n_comparisons =
        dplyr::n(),
      mean_abs_delta_rho =
        safe_mean(
          abs_delta_rho
        ),
      median_abs_delta_rho =
        safe_median(
          abs_delta_rho
        ),
      q90_abs_delta_rho =
        safe_quantile(
          abs_delta_rho,
          0.90
        ),
      max_abs_delta_rho =
        safe_max(
          abs_delta_rho
        ),
      direction_consistency_all =
        safe_proportion(
          same_direction
        ),
      direction_consistency_full_abs_rho_ge_0_10 =
        safe_proportion(
          same_direction[
            abs_rho_full >= 0.10
          ]
        ),
      direction_consistency_full_FDR_lt_0_05 =
        safe_proportion(
          same_direction[
            full_FDR_lt_0_05
          ]
        ),
      full_significant_direction_and_nominal_preservation =
        safe_proportion(
          preserved_direction_and_nominal_if_full_sig
        ),
      full_significant_direction_preservation =
        safe_proportion(
          preserved_direction_if_full_sig
        ),
      sensitivity_nominal_significance_rate =
        safe_proportion(
          sensitivity_nominal_lt_0_05
        ),
      sensitivity_FDR_significance_rate =
        safe_proportion(
          sensitivity_FDR_lt_0_05
        ),
      .groups = "drop"
    ) %>%
    dplyr::arrange(
      mean_abs_delta_rho
    )
}

summarize_correlation_delta_by_module_trait <- function(
    delta_tbl
) {
  summary_tbl <- delta_tbl %>%
    dplyr::group_by(
      Module,
      Trait
    ) %>%
    dplyr::summarise(
      rho_full = dplyr::first(
        rho_full
      ),
      FDR_full = dplyr::first(
        FDR_full
      ),
      n_runs =
        dplyr::n_distinct(
          Run_ID
        ),
      mean_rho_sensitivity =
        safe_mean(rho),
      median_rho_sensitivity =
        safe_median(rho),
      minimum_rho_sensitivity =
        safe_min(rho),
      maximum_rho_sensitivity =
        safe_max(rho),
      empirical_rho_2_5 =
        safe_quantile(
          rho,
          0.025
        ),
      empirical_rho_97_5 =
        safe_quantile(
          rho,
          0.975
        ),
      mean_abs_delta_rho =
        safe_mean(
          abs_delta_rho
        ),
      median_abs_delta_rho =
        safe_median(
          abs_delta_rho
        ),
      q90_abs_delta_rho =
        safe_quantile(
          abs_delta_rho,
          0.90
        ),
      max_abs_delta_rho =
        safe_max(
          abs_delta_rho
        ),
      direction_consistency =
        safe_proportion(
          same_direction
        ),
      nominal_significance_rate =
        safe_proportion(
          sensitivity_nominal_lt_0_05
        ),
      FDR_significance_rate =
        safe_proportion(
          sensitivity_FDR_lt_0_05
        ),
      full_significant_direction_and_nominal_preservation =
        safe_proportion(
          preserved_direction_and_nominal_if_full_sig
        ),
      full_rho_inside_empirical_95_interval =
        dplyr::case_when(
          is.na(rho_full) |
            is.na(
              empirical_rho_2_5
            ) |
            is.na(
              empirical_rho_97_5
            ) ~ NA,
          TRUE ~
            rho_full >=
              empirical_rho_2_5 &
            rho_full <=
              empirical_rho_97_5
        ),
      .groups = "drop"
    )

  worst_runs <- delta_tbl %>%
    dplyr::group_by(
      Module,
      Trait
    ) %>%
    dplyr::slice_max(
      order_by = abs_delta_rho,
      n = 1,
      with_ties = FALSE
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(
      Module,
      Trait,
      worst_run_ID = Run_ID,
      worst_excluded_unit =
        Excluded_unit,
      worst_abs_delta_rho =
        abs_delta_rho
    )

  summary_tbl %>%
    dplyr::left_join(
      worst_runs,
      by = c(
        "Module",
        "Trait"
      )
    )
}

summarize_correlation_delta_by_run <- function(
    delta_tbl
) {
  delta_tbl %>%
    dplyr::group_by(
      Analysis,
      Run_ID,
      Excluded_unit,
      Iteration
    ) %>%
    dplyr::summarise(
      N_total_run =
        dplyr::first(
          N_total_run
        ),
      n_tests =
        dplyr::n(),
      mean_abs_delta_rho =
        safe_mean(
          abs_delta_rho
        ),
      median_abs_delta_rho =
        safe_median(
          abs_delta_rho
        ),
      q90_abs_delta_rho =
        safe_quantile(
          abs_delta_rho,
          0.90
        ),
      max_abs_delta_rho =
        safe_max(
          abs_delta_rho
        ),
      direction_consistency =
        safe_proportion(
          same_direction
        ),
      direction_consistency_full_FDR_lt_0_05 =
        safe_proportion(
          same_direction[
            full_FDR_lt_0_05
          ]
        ),
      full_significant_direction_and_nominal_preservation =
        safe_proportion(
          preserved_direction_and_nominal_if_full_sig
        ),
      .groups = "drop"
    )
}

###############################################################################
# 6) HC3 MODEL HELPERS
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

fit_diagnosis_hc3_single <- function(
    df,
    module
) {
  required_vars <- c(
    module,
    "SampleGroup_bin",
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
      N = 0L,
      estimate_diagnosis = NA_real_,
      robust_se = NA_real_,
      robust_p_value = NA_real_,
      robust_conf_low = NA_real_,
      robust_conf_high = NA_real_,
      module_SD = NA_real_,
      standardized_difference =
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
    ) %>%
    tidyr::drop_na()

  if (
    nrow(d) < MIN_N_FOR_MODEL ||
    dplyr::n_distinct(
      d$SampleGroup_bin
    ) < 2 ||
    dplyr::n_distinct(
      d[[module]]
    ) <= 1
  ) {
    return(tibble::tibble(
      Module = module,
      N = nrow(d),
      estimate_diagnosis = NA_real_,
      robust_se = NA_real_,
      robust_p_value = NA_real_,
      robust_conf_low = NA_real_,
      robust_conf_high = NA_real_,
      module_SD = stats::sd(
        d[[module]],
        na.rm = TRUE
      ),
      standardized_difference =
        NA_real_,
      model_status =
        "insufficient_or_constant"
    ))
  }

  d$Country <- droplevels(
    factor(d$Country)
  )

  fit <- tryCatch(
    stats::lm(
      stats::as.formula(
        paste0(
          module,
          " ~ SampleGroup_bin + Age + Sex_bin + Education + Country"
        )
      ),
      data = d
    ),
    error = function(e) NULL
  )

  robust_term <- extract_hc3_term(
    fit,
    term = "SampleGroup_bin"
  )

  module_sd <- stats::sd(
    d[[module]],
    na.rm = TRUE
  )

  standardized_difference <- if (
    is.finite(module_sd) &&
    module_sd > 0
  ) {
    robust_term$estimate /
      module_sd
  } else {
    NA_real_
  }

  tibble::tibble(
    Module = module,
    N = nrow(d),
    estimate_diagnosis =
      robust_term$estimate,
    robust_se =
      robust_term$robust_se,
    robust_p_value =
      robust_term$robust_p_value,
    robust_conf_low =
      robust_term$robust_conf_low,
    robust_conf_high =
      robust_term$robust_conf_high,
    module_SD = module_sd,
    standardized_difference =
      standardized_difference,
    model_status =
      robust_term$model_status
  )
}

run_diagnosis_hc3_models <- function(
    df,
    modules,
    analysis,
    run_id,
    excluded_unit = NA_character_,
    iteration = NA_integer_
) {
  purrr::map_dfr(
    modules,
    ~ fit_diagnosis_hc3_single(
      df,
      .x
    )
  ) %>%
    dplyr::mutate(
      FDR_HC3 = p.adjust(
        robust_p_value,
        method = "BH"
      ),
      Analysis = analysis,
      Run_ID = run_id,
      Excluded_unit =
        excluded_unit,
      Iteration = iteration,
      N_total_run = nrow(df)
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

  fit <- tryCatch(
    stats::lm(
      stats::as.formula(
        paste0(
          ".log_outcome ~ ",
          module,
          " + Age + Sex_bin + Education + Country"
        )
      ),
      data = d
    ),
    error = function(e) NULL
  )

  robust_term <- extract_hc3_term(
    fit,
    term = module
  )

  module_sd <- stats::sd(
    d[[module]],
    na.rm = TRUE
  )

  outcome_sd <- stats::sd(
    d$.log_outcome,
    na.rm = TRUE
  )

  effect_per_sd <- robust_term$estimate *
    module_sd

  standardized_beta <- if (
    is.finite(outcome_sd) &&
    outcome_sd > 0
  ) {
    effect_per_sd /
      outcome_sd
  } else {
    NA_real_
  }

  percent_change_per_sd <- if (
    is.finite(effect_per_sd)
  ) {
    100 *
      (
        exp(
          effect_per_sd
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
    robust_p_value =
      robust_term$robust_p_value,
    robust_conf_low =
      robust_term$robust_conf_low,
    robust_conf_high =
      robust_term$robust_conf_high,
    module_SD = module_sd,
    log_outcome_SD = outcome_sd,
    standardized_beta =
      standardized_beta,
    percent_change_per_1SD_module =
      percent_change_per_sd,
    model_status =
      robust_term$model_status
  )
}

run_log_biomarker_hc3_models <- function(
    df,
    modules,
    biomarkers,
    analysis,
    run_id,
    excluded_unit = NA_character_,
    iteration = NA_integer_
) {
  purrr::map_dfr(
    modules,
    function(mod) {
      purrr::map_dfr(
        biomarkers,
        function(bio) {
          fit_log_biomarker_hc3_single(
            df,
            module = mod,
            biomarker = bio
          )
        }
      )
    }
  ) %>%
    dplyr::mutate(
      FDR_HC3 = p.adjust(
        robust_p_value,
        method = "BH"
      ),
      Analysis = analysis,
      Run_ID = run_id,
      Excluded_unit =
        excluded_unit,
      Iteration = iteration,
      N_total_run = nrow(df)
    )
}

attach_full_model_reference <- function(
    sensitivity_tbl,
    full_tbl,
    keys,
    effect_col,
    full_effect_col = "effect_full"
) {
  full_reference <- full_tbl %>%
    dplyr::select(
      dplyr::all_of(
        c(
          keys,
          effect_col
        )
      )
    )

  names(full_reference)[
    names(full_reference) ==
      effect_col
  ] <- full_effect_col

  sensitivity_tbl %>%
    dplyr::left_join(
      full_reference,
      by = keys
    ) %>%
    dplyr::mutate(
      delta_effect =
        .data[[effect_col]] -
        .data[[full_effect_col]],
      abs_delta_effect =
        abs(delta_effect),
      same_direction = dplyr::case_when(
        is.na(
          .data[[effect_col]]
        ) |
          is.na(
            .data[[full_effect_col]]
          ) ~ NA,
        .data[[effect_col]] == 0 |
          .data[[full_effect_col]] == 0 ~ NA,
        TRUE ~
          sign(
            .data[[effect_col]]
          ) ==
            sign(
              .data[[full_effect_col]]
            )
      ),
      nominal_significant =
        !is.na(robust_p_value) &
        robust_p_value < 0.05,
      FDR_significant =
        !is.na(FDR_HC3) &
        FDR_HC3 < 0.05
    )
}

summarize_model_stability <- function(
    delta_tbl,
    group_vars,
    effect_col,
    full_effect_col = "effect_full"
) {
  delta_tbl %>%
    dplyr::group_by(
      dplyr::across(
        dplyr::all_of(
          group_vars
        )
      )
    ) %>%
    dplyr::summarise(
      full_effect =
        dplyr::first(
          .data[[full_effect_col]]
        ),
      n_runs =
        dplyr::n_distinct(
          Run_ID
        ),
      mean_effect =
        safe_mean(
          .data[[effect_col]]
        ),
      median_effect =
        safe_median(
          .data[[effect_col]]
        ),
      empirical_effect_2_5 =
        safe_quantile(
          .data[[effect_col]],
          0.025
        ),
      empirical_effect_97_5 =
        safe_quantile(
          .data[[effect_col]],
          0.975
        ),
      mean_abs_delta_effect =
        safe_mean(
          abs_delta_effect
        ),
      max_abs_delta_effect =
        safe_max(
          abs_delta_effect
        ),
      direction_consistency =
        safe_proportion(
          same_direction
        ),
      nominal_significance_rate =
        safe_proportion(
          nominal_significant
        ),
      FDR_significance_rate =
        safe_proportion(
          FDR_significant
        ),
      full_effect_inside_empirical_95_interval =
        dplyr::case_when(
          is.na(full_effect) |
            is.na(
              empirical_effect_2_5
            ) |
            is.na(
              empirical_effect_97_5
            ) ~ NA,
          TRUE ~
            full_effect >=
              empirical_effect_2_5 &
            full_effect <=
              empirical_effect_97_5
        ),
      .groups = "drop"
    )
}

###############################################################################
# 7) LOAD AND PREPARE DATA
###############################################################################

analysis_df <- safe_read_csv(
  INPUT_FILE
)

module_reference <- safe_read_csv(
  MODULE_REFERENCE_FILE
)

full_correlation_reference <- safe_read_csv(
  FULL_CORRELATION_REFERENCE_FILE
)

prioritization_reference <- safe_read_csv(
  PRIORITIZATION_FILE
)

diagnosis_reference <- safe_read_csv(
  DIAGNOSIS_REFERENCE_FILE
)

biomarker_reference <- safe_read_csv(
  BIOMARKER_REFERENCE_FILE
)

site_structure_reference <- safe_read_csv(
  SITE_STRUCTURE_REFERENCE_FILE
)

required_input_cols <- c(
  "SampleId",
  "SampleGroup",
  "SampleGroup_bin",
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
    "module_trait_input_clean.csv is missing required columns: ",
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
    SampleGroup = clean_text_na(
      SampleGroup
    ),
    SampleGroup_bin = safe_numeric(
      SampleGroup_bin
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

modules_use <- module_reference$Module[
  module_reference$Module %in%
    names(analysis_df)
]

modules_use <- modules_use[
  !modules_use %in%
    EXCLUDE_MODULES
]

modules_use <- modules_use[
  vapply(
    analysis_df[modules_use],
    is.numeric,
    logical(1)
  )
]

traits_available <- intersect(
  TRAITS_USE,
  names(analysis_df)
)

core_modules_use <- intersect(
  CORE_MODULES,
  modules_use
)

focal_biomarker_modules_use <- intersect(
  FOCAL_BIOMARKER_MODULES,
  modules_use
)

if (length(modules_use) == 0) {
  stop(
    "No module eigengene columns were detected.",
    call. = FALSE
  )
}

if (length(traits_available) == 0) {
  stop(
    "No eligible traits were detected.",
    call. = FALSE
  )
}

if (length(core_modules_use) == 0) {
  stop(
    "None of the prespecified core modules were detected.",
    call. = FALSE
  )
}

if (
  length(focal_biomarker_modules_use) ==
    0
) {
  stop(
    "None of the focal biomarker modules were detected.",
    call. = FALSE
  )
}

for (mod in modules_use) {
  analysis_df[[mod]] <- safe_numeric(
    analysis_df[[mod]]
  )
}

for (trait in traits_available) {
  analysis_df[[trait]] <- safe_numeric(
    analysis_df[[trait]]
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

country_site_map <- analysis_df %>%
  dplyr::filter(
    !is.na(Country),
    !is.na(.data[[site_var]])
  ) %>%
  dplyr::distinct(
    Country,
    Site = .data[[site_var]]
  ) %>%
  dplyr::group_by(Country) %>%
  dplyr::mutate(
    n_sites_in_country =
      dplyr::n()
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    Country = as.character(
      Country
    ),
    Site = as.character(Site),
    Country_Site = paste(
      Country,
      Site,
      sep = "::"
    ),
    country_remains_after_site_exclusion =
      n_sites_in_country >= 2,
    equivalent_to_country_removal =
      n_sites_in_country == 1
  ) %>%
  dplyr::arrange(
    Country,
    Site
  )

primary_wc_loso_units <- country_site_map %>%
  dplyr::filter(
    country_remains_after_site_exclusion
  )

all_site_units <- country_site_map

strata_counts <- analysis_df %>%
  dplyr::count(
    Country,
    SampleGroup_bin,
    name = "N"
  ) %>%
  tidyr::drop_na() %>%
  dplyr::arrange(
    Country,
    SampleGroup_bin
  )

min_stratum_n <- min(
  strata_counts$N,
  na.rm = TRUE
)

downsampling_total_n <- min_stratum_n *
  nrow(strata_counts)

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

  if (length(traits_available) !=
      EXPECTED_N_TRAITS) {
    stop(
      "Expected ",
      EXPECTED_N_TRAITS,
      " traits, but found ",
      length(traits_available),
      ".",
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

  if (
    dplyr::n_distinct(
      analysis_df[[site_var]]
    ) !=
      EXPECTED_N_SITE_LEVELS
  ) {
    stop(
      "Expected ",
      EXPECTED_N_SITE_LEVELS,
      " site levels, but found ",
      dplyr::n_distinct(
        analysis_df[[site_var]]
      ),
      ".",
      call. = FALSE
    )
  }

  if (nrow(strata_counts) !=
      EXPECTED_N_COUNTRY_DIAGNOSIS_STRATA) {
    stop(
      "Expected ",
      EXPECTED_N_COUNTRY_DIAGNOSIS_STRATA,
      " Country × diagnosis strata, but found ",
      nrow(strata_counts),
      ".",
      call. = FALSE
    )
  }
}

if (
  min_stratum_n <
    MIN_N_FOR_CORR
) {
  stop(
    "The smallest Country × diagnosis stratum contains only ",
    min_stratum_n,
    " participants, which is insufficient for balanced downsampling.",
    call. = FALSE
  )
}

safe_write_csv(
  analysis_df,
  file.path(
    OUTDIR,
    "tables",
    "sensitivity_input_clean.csv"
  )
)

safe_write_csv(
  country_site_map,
  file.path(
    OUTDIR,
    "tables",
    "within_country_loso",
    "country_site_deletion_map.csv"
  )
)

safe_write_csv(
  strata_counts,
  file.path(
    OUTDIR,
    "tables",
    "downsampling",
    "country_diagnosis_strata_counts.csv"
  )
)

input_audit <- tibble::tibble(
  metric = c(
    "base_dir",
    "n_samples",
    "n_modules",
    "modules",
    "n_traits",
    "traits",
    "n_countries",
    "countries",
    "site_variable",
    "n_site_levels",
    "n_primary_within_country_site_deletions",
    "primary_within_country_site_deletions",
    "n_all_site_deletion_audit_units",
    "n_country_diagnosis_strata",
    "minimum_country_diagnosis_stratum_n",
    "balanced_downsampling_total_n",
    "downsampling_iterations",
    "random_seed",
    "Country_numeric_present",
    "wgcna_input_level"
  ),
  value = c(
    BASE_DIR,
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
      length(traits_available)
    ),
    paste(
      traits_available,
      collapse = ", "
    ),
    as.character(
      length(countries)
    ),
    paste(
      countries,
      collapse = ", "
    ),
    site_var,
    as.character(
      dplyr::n_distinct(
        analysis_df[[site_var]]
      )
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
      nrow(all_site_units)
    ),
    as.character(
      nrow(strata_counts)
    ),
    as.character(
      min_stratum_n
    ),
    as.character(
      downsampling_total_n
    ),
    as.character(
      N_DOWNSAMPLING_ITER
    ),
    as.character(SEED),
    as.character(
      "Country_numeric" %in%
        names(analysis_df)
    ),
    "GENE-COLLAPSED, outcome-independent SOMAmer selection"
  )
)

safe_write_csv(
  input_audit,
  file.path(
    OUTDIR,
    "tables",
    "script14_input_alignment_audit.csv"
  )
)

cat("Samples:", nrow(analysis_df), "\n")
cat("Modules:", paste(modules_use, collapse = ", "), "\n")
cat("Core modules:", paste(core_modules_use, collapse = ", "), "\n")
cat("Traits:", paste(traits_available, collapse = ", "), "\n")
cat("Countries:", paste(countries, collapse = ", "), "\n")
cat(
  "Primary WC-LOSO units:",
  paste(
    primary_wc_loso_units$
      Country_Site,
    collapse = ", "
  ),
  "\n"
)
cat(
  "Balanced downsampling:",
  min_stratum_n,
  "per Country × diagnosis stratum;",
  downsampling_total_n,
  "participants per iteration.\n\n"
)

###############################################################################
# 8) FULL MODULE-TRAIT CORRELATION RECOMPUTATION AND AUDIT
###############################################################################

full_module_trait <- compute_module_trait(
  df = analysis_df,
  modules = modules_use,
  traits = traits_available,
  analysis = "full",
  run_id = "full"
) %>%
  dplyr::mutate(
    Trait_label = dplyr::recode(
      Trait,
      !!!TRAIT_LABELS,
      .default = Trait
    )
  )

full_reference_for_audit <- full_correlation_reference %>%
  dplyr::mutate(
    Module = as.character(Module),
    Trait = as.character(Trait)
  ) %>%
  dplyr::select(
    Module,
    Trait,
    rho_reference = rho,
    p_reference = p_value,
    FDR_reference = FDR,
    N_reference = N
  )

full_correlation_audit <- full_module_trait %>%
  dplyr::left_join(
    full_reference_for_audit,
    by = c(
      "Module",
      "Trait"
    )
  ) %>%
  dplyr::mutate(
    abs_delta_rho_reference =
      abs(
        rho -
          rho_reference
      ),
    abs_delta_p_reference =
      abs(
        p_value -
          p_reference
      ),
    abs_delta_FDR_reference =
      abs(
        FDR -
          FDR_reference
      ),
    N_matches_reference =
      N ==
        N_reference
  )

correlation_reference_summary <- tibble::tibble(
  metric = c(
    "n_recomputed_tests",
    "n_reference_tests",
    "n_matched_tests",
    "max_abs_delta_rho",
    "max_abs_delta_p",
    "max_abs_delta_FDR",
    "all_N_match",
    "tolerance",
    "reference_match_passed"
  ),
  value = c(
    as.character(
      nrow(full_module_trait)
    ),
    as.character(
      nrow(
        full_correlation_reference
      )
    ),
    as.character(
      sum(
        !is.na(
          full_correlation_audit$
            rho_reference
        )
      )
    ),
    as.character(
      safe_max(
        full_correlation_audit$
          abs_delta_rho_reference
      )
    ),
    as.character(
      safe_max(
        full_correlation_audit$
          abs_delta_p_reference
      )
    ),
    as.character(
      safe_max(
        full_correlation_audit$
          abs_delta_FDR_reference
      )
    ),
    as.character(
      all(
        full_correlation_audit$
          N_matches_reference,
        na.rm = TRUE
      )
    ),
    as.character(
      FULL_CORRELATION_TOLERANCE
    ),
    as.character(
      safe_max(
        full_correlation_audit$
          abs_delta_rho_reference
      ) <=
        FULL_CORRELATION_TOLERANCE &&
        safe_max(
          full_correlation_audit$
            abs_delta_p_reference
        ) <=
          FULL_CORRELATION_TOLERANCE &&
        safe_max(
          full_correlation_audit$
            abs_delta_FDR_reference
        ) <=
          FULL_CORRELATION_TOLERANCE &&
        all(
          full_correlation_audit$
            N_matches_reference,
          na.rm = TRUE
        )
    )
  )
)

safe_write_csv(
  full_module_trait,
  file.path(
    OUTDIR,
    "tables",
    "full",
    "full_module_trait_correlations_recomputed.csv"
  )
)

safe_write_csv(
  full_correlation_audit,
  file.path(
    OUTDIR,
    "tables",
    "full",
    "full_module_trait_reference_audit.csv"
  )
)

safe_write_csv(
  correlation_reference_summary,
  file.path(
    OUTDIR,
    "tables",
    "full",
    "full_module_trait_reference_audit_summary.csv"
  )
)

reference_match_passed <-
  safe_max(
    full_correlation_audit$
      abs_delta_rho_reference
  ) <=
    FULL_CORRELATION_TOLERANCE &&
  safe_max(
    full_correlation_audit$
      abs_delta_p_reference
  ) <=
    FULL_CORRELATION_TOLERANCE &&
  safe_max(
    full_correlation_audit$
      abs_delta_FDR_reference
  ) <=
    FULL_CORRELATION_TOLERANCE &&
  all(
    full_correlation_audit$
      N_matches_reference,
    na.rm = TRUE
  )

if (
  REQUIRE_FULL_CORRELATION_MATCH &&
  !reference_match_passed
) {
  stop(
    "The independently recomputed full module-trait matrix does not match ",
    "Script 13 within tolerance. Review tables/full/",
    "full_module_trait_reference_audit.csv.",
    call. = FALSE
  )
}

###############################################################################
# 9) FULL HC3 MODEL REFERENCES AND AUDIT
###############################################################################

full_diagnosis_models <- run_diagnosis_hc3_models(
  df = analysis_df,
  modules = modules_use,
  analysis = "full",
  run_id = "full"
)

full_biomarker_models <- run_log_biomarker_hc3_models(
  df = analysis_df,
  modules =
    focal_biomarker_modules_use,
  biomarkers = BIOMARKERS,
  analysis = "full",
  run_id = "full"
)

diagnosis_reference_audit <- full_diagnosis_models %>%
  dplyr::left_join(
    diagnosis_reference %>%
      dplyr::transmute(
        Module = as.character(Module),
        estimate_reference =
          safe_numeric(
            estimate_diagnosis
          ),
        standardized_reference =
          safe_numeric(
            standardized_difference
          ),
        p_reference =
          safe_numeric(
            robust_p_value
          ),
        FDR_reference =
          safe_numeric(
            FDR_HC3
          )
      ),
    by = "Module"
  ) %>%
  dplyr::mutate(
    abs_delta_estimate =
      abs(
        estimate_diagnosis -
          estimate_reference
      ),
    abs_delta_standardized =
      abs(
        standardized_difference -
          standardized_reference
      ),
    abs_delta_p =
      abs(
        robust_p_value -
          p_reference
      ),
    abs_delta_FDR =
      abs(
        FDR_HC3 -
          FDR_reference
      )
  )

biomarker_reference_audit <- full_biomarker_models %>%
  dplyr::left_join(
    biomarker_reference %>%
      dplyr::filter(
        Module %in%
          focal_biomarker_modules_use,
        Biomarker %in%
          BIOMARKERS
      ) %>%
      dplyr::transmute(
        Module = as.character(Module),
        Biomarker =
          as.character(Biomarker),
        estimate_reference =
          safe_numeric(estimate),
        standardized_reference =
          safe_numeric(
            standardized_beta
          ),
        p_reference =
          safe_numeric(
            robust_p_value
          ),
        FDR_reference =
          safe_numeric(
            FDR_within_transform_family
          )
      ),
    by = c(
      "Module",
      "Biomarker"
    )
  ) %>%
  dplyr::mutate(
    abs_delta_estimate =
      abs(
        estimate_log_scale -
          estimate_reference
      ),
    abs_delta_standardized =
      abs(
        standardized_beta -
          standardized_reference
      ),
    abs_delta_p =
      abs(
        robust_p_value -
          p_reference
      ),
    abs_delta_FDR =
      abs(
        FDR_HC3 -
          FDR_reference
      )
  )

model_reference_summary <- tibble::tibble(
  model_family = c(
    "diagnosis_HC3",
    "focal_log_biomarker_HC3"
  ),
  n_models = c(
    nrow(
      diagnosis_reference_audit
    ),
    nrow(
      biomarker_reference_audit
    )
  ),
  max_abs_delta_estimate = c(
    safe_max(
      diagnosis_reference_audit$
        abs_delta_estimate
    ),
    safe_max(
      biomarker_reference_audit$
        abs_delta_estimate
    )
  ),
  max_abs_delta_standardized = c(
    safe_max(
      diagnosis_reference_audit$
        abs_delta_standardized
    ),
    safe_max(
      biomarker_reference_audit$
        abs_delta_standardized
    )
  ),
  max_abs_delta_p = c(
    safe_max(
      diagnosis_reference_audit$
        abs_delta_p
    ),
    safe_max(
      biomarker_reference_audit$
        abs_delta_p
    )
  ),
  max_abs_delta_FDR = c(
    safe_max(
      diagnosis_reference_audit$
        abs_delta_FDR
    ),
    safe_max(
      biomarker_reference_audit$
        abs_delta_FDR
    )
  ),
  tolerance = MODEL_REFERENCE_TOLERANCE,
  match_passed = c(
    safe_max(
      diagnosis_reference_audit$
        abs_delta_estimate
    ) <=
      MODEL_REFERENCE_TOLERANCE,
    safe_max(
      biomarker_reference_audit$
        abs_delta_estimate
    ) <=
      MODEL_REFERENCE_TOLERANCE
  )
)

safe_write_csv(
  full_diagnosis_models,
  file.path(
    OUTDIR,
    "tables",
    "full",
    "full_diagnosis_HC3_models_recomputed.csv"
  )
)

safe_write_csv(
  full_biomarker_models,
  file.path(
    OUTDIR,
    "tables",
    "full",
    "full_focal_log_biomarker_HC3_models_recomputed.csv"
  )
)

safe_write_csv(
  diagnosis_reference_audit,
  file.path(
    OUTDIR,
    "tables",
    "full",
    "full_diagnosis_HC3_reference_audit.csv"
  )
)

safe_write_csv(
  biomarker_reference_audit,
  file.path(
    OUTDIR,
    "tables",
    "full",
    "full_focal_log_biomarker_HC3_reference_audit.csv"
  )
)

safe_write_csv(
  model_reference_summary,
  file.path(
    OUTDIR,
    "tables",
    "full",
    "full_HC3_model_reference_audit_summary.csv"
  )
)

if (
  any(
    !model_reference_summary$
      match_passed
  )
) {
  warning(
    "At least one recomputed HC3 model family differs from Script 13b ",
    "beyond the reference tolerance. The script continues, but the audit ",
    "tables must be reviewed."
  )
}

###############################################################################
# 10) LEAVE-ONE-COUNTRY-OUT
###############################################################################

loco_correlation_results <- list()
loco_diagnosis_results <- list()
loco_biomarker_results <- list()

for (cc in countries) {
  cat(
    "LOCO: excluding country ",
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

  loco_correlation_results[[cc]] <-
    compute_module_trait(
      df = df_run,
      modules = modules_use,
      traits = traits_available,
      analysis = "LOCO_country",
      run_id = run_id,
      excluded_unit = cc
    )

  loco_diagnosis_results[[cc]] <-
    run_diagnosis_hc3_models(
      df = df_run,
      modules = modules_use,
      analysis = "LOCO_country",
      run_id = run_id,
      excluded_unit = cc
    )

  loco_biomarker_results[[cc]] <-
    run_log_biomarker_hc3_models(
      df = df_run,
      modules =
        focal_biomarker_modules_use,
      biomarkers = BIOMARKERS,
      analysis = "LOCO_country",
      run_id = run_id,
      excluded_unit = cc
    )
}

loco_correlation_long <- dplyr::bind_rows(
  loco_correlation_results
)

loco_correlation_delta <- attach_full_correlation_reference(
  sensitivity_tbl =
    loco_correlation_long,
  full_tbl =
    full_module_trait
)

loco_summary_by_module <-
  summarize_correlation_delta_by_module(
    loco_correlation_delta
  )

loco_summary_by_module_trait <-
  summarize_correlation_delta_by_module_trait(
    loco_correlation_delta
  )

loco_summary_by_excluded_country <-
  summarize_correlation_delta_by_run(
    loco_correlation_delta
  )

loco_diagnosis_long <- dplyr::bind_rows(
  loco_diagnosis_results
)

loco_diagnosis_delta <- attach_full_model_reference(
  sensitivity_tbl =
    loco_diagnosis_long,
  full_tbl =
    full_diagnosis_models,
  keys = "Module",
  effect_col =
    "standardized_difference",
  full_effect_col =
    "standardized_difference_full"
)

loco_diagnosis_summary <-
  summarize_model_stability(
    delta_tbl =
      loco_diagnosis_delta,
    group_vars = "Module",
    effect_col =
      "standardized_difference",
    full_effect_col =
      "standardized_difference_full"
  )

loco_biomarker_long <- dplyr::bind_rows(
  loco_biomarker_results
)

loco_biomarker_delta <- attach_full_model_reference(
  sensitivity_tbl =
    loco_biomarker_long,
  full_tbl =
    full_biomarker_models,
  keys = c(
    "Module",
    "Biomarker"
  ),
  effect_col =
    "standardized_beta",
  full_effect_col =
    "standardized_beta_full"
)

loco_biomarker_summary <-
  summarize_model_stability(
    delta_tbl =
      loco_biomarker_delta,
    group_vars = c(
      "Module",
      "Biomarker"
    ),
    effect_col =
      "standardized_beta",
    full_effect_col =
      "standardized_beta_full"
  )

safe_write_csv(
  loco_correlation_long,
  file.path(
    OUTDIR,
    "tables",
    "loco",
    "loco_country_module_trait_results.csv"
  )
)

safe_write_csv(
  loco_correlation_delta,
  file.path(
    OUTDIR,
    "tables",
    "loco",
    "loco_country_delta_results.csv"
  )
)

safe_write_csv(
  loco_summary_by_module,
  file.path(
    OUTDIR,
    "tables",
    "loco",
    "loco_country_summary_by_module.csv"
  )
)

safe_write_csv(
  loco_summary_by_module_trait,
  file.path(
    OUTDIR,
    "tables",
    "loco",
    "loco_country_summary_by_module_trait.csv"
  )
)

safe_write_csv(
  loco_summary_by_excluded_country,
  file.path(
    OUTDIR,
    "tables",
    "loco",
    "loco_country_summary_by_excluded_country.csv"
  )
)

safe_write_csv(
  loco_diagnosis_long,
  file.path(
    OUTDIR,
    "tables",
    "model_stability",
    "diagnosis",
    "loco_diagnosis_HC3_results.csv"
  )
)

safe_write_csv(
  loco_diagnosis_delta,
  file.path(
    OUTDIR,
    "tables",
    "model_stability",
    "diagnosis",
    "loco_diagnosis_HC3_delta_results.csv"
  )
)

safe_write_csv(
  loco_diagnosis_summary,
  file.path(
    OUTDIR,
    "tables",
    "model_stability",
    "diagnosis",
    "loco_diagnosis_HC3_summary_by_module.csv"
  )
)

safe_write_csv(
  loco_biomarker_long,
  file.path(
    OUTDIR,
    "tables",
    "model_stability",
    "biomarker",
    "loco_focal_log_biomarker_HC3_results.csv"
  )
)

safe_write_csv(
  loco_biomarker_delta,
  file.path(
    OUTDIR,
    "tables",
    "model_stability",
    "biomarker",
    "loco_focal_log_biomarker_HC3_delta_results.csv"
  )
)

safe_write_csv(
  loco_biomarker_summary,
  file.path(
    OUTDIR,
    "tables",
    "model_stability",
    "biomarker",
    "loco_focal_log_biomarker_HC3_summary.csv"
  )
)

###############################################################################
# 11) PRIMARY CORRECTED WITHIN-COUNTRY LOSO
###############################################################################

wc_loso_correlation_results <- list()
wc_loso_diagnosis_results <- list()
wc_loso_biomarker_results <- list()

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
    "WC-LOSO: excluding ",
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

  wc_loso_correlation_results[[
    unit_i
  ]] <- compute_module_trait(
    df = df_run,
    modules = modules_use,
    traits = traits_available,
    analysis =
      "WC_LOSO_primary",
    run_id = run_id,
    excluded_unit = unit_i
  )

  wc_loso_diagnosis_results[[
    unit_i
  ]] <- run_diagnosis_hc3_models(
    df = df_run,
    modules = modules_use,
    analysis =
      "WC_LOSO_primary",
    run_id = run_id,
    excluded_unit = unit_i
  )

  wc_loso_biomarker_results[[
    unit_i
  ]] <-
    run_log_biomarker_hc3_models(
      df = df_run,
      modules =
        focal_biomarker_modules_use,
      biomarkers = BIOMARKERS,
      analysis =
        "WC_LOSO_primary",
      run_id = run_id,
      excluded_unit = unit_i
    )
}

wc_loso_correlation_long <-
  dplyr::bind_rows(
    wc_loso_correlation_results
  )

wc_loso_correlation_delta <-
  attach_full_correlation_reference(
    sensitivity_tbl =
      wc_loso_correlation_long,
    full_tbl =
      full_module_trait
  )

wc_loso_summary_by_module <-
  summarize_correlation_delta_by_module(
    wc_loso_correlation_delta
  )

wc_loso_summary_by_module_trait <-
  summarize_correlation_delta_by_module_trait(
    wc_loso_correlation_delta
  )

wc_loso_summary_by_excluded_site <-
  summarize_correlation_delta_by_run(
    wc_loso_correlation_delta
  )

wc_loso_diagnosis_long <-
  dplyr::bind_rows(
    wc_loso_diagnosis_results
  )

wc_loso_diagnosis_delta <-
  attach_full_model_reference(
    sensitivity_tbl =
      wc_loso_diagnosis_long,
    full_tbl =
      full_diagnosis_models,
    keys = "Module",
    effect_col =
      "standardized_difference",
    full_effect_col =
      "standardized_difference_full"
  )

wc_loso_diagnosis_summary <-
  summarize_model_stability(
    delta_tbl =
      wc_loso_diagnosis_delta,
    group_vars = "Module",
    effect_col =
      "standardized_difference",
    full_effect_col =
      "standardized_difference_full"
  )

wc_loso_biomarker_long <-
  dplyr::bind_rows(
    wc_loso_biomarker_results
  )

wc_loso_biomarker_delta <-
  attach_full_model_reference(
    sensitivity_tbl =
      wc_loso_biomarker_long,
    full_tbl =
      full_biomarker_models,
    keys = c(
      "Module",
      "Biomarker"
    ),
    effect_col =
      "standardized_beta",
    full_effect_col =
      "standardized_beta_full"
  )

wc_loso_biomarker_summary <-
  summarize_model_stability(
    delta_tbl =
      wc_loso_biomarker_delta,
    group_vars = c(
      "Module",
      "Biomarker"
    ),
    effect_col =
      "standardized_beta",
    full_effect_col =
      "standardized_beta_full"
  )

safe_write_csv(
  wc_loso_correlation_long,
  file.path(
    OUTDIR,
    "tables",
    "within_country_loso",
    "within_country_loso_module_trait_results.csv"
  )
)

safe_write_csv(
  wc_loso_correlation_delta,
  file.path(
    OUTDIR,
    "tables",
    "within_country_loso",
    "within_country_loso_delta_results.csv"
  )
)

safe_write_csv(
  wc_loso_summary_by_module,
  file.path(
    OUTDIR,
    "tables",
    "within_country_loso",
    "within_country_loso_summary_by_module.csv"
  )
)

safe_write_csv(
  wc_loso_summary_by_module_trait,
  file.path(
    OUTDIR,
    "tables",
    "within_country_loso",
    "within_country_loso_summary_by_module_trait.csv"
  )
)

safe_write_csv(
  wc_loso_summary_by_excluded_site,
  file.path(
    OUTDIR,
    "tables",
    "within_country_loso",
    "within_country_loso_summary_by_excluded_site.csv"
  )
)

safe_write_csv(
  wc_loso_diagnosis_long,
  file.path(
    OUTDIR,
    "tables",
    "model_stability",
    "diagnosis",
    "within_country_loso_diagnosis_HC3_results.csv"
  )
)

safe_write_csv(
  wc_loso_diagnosis_delta,
  file.path(
    OUTDIR,
    "tables",
    "model_stability",
    "diagnosis",
    "within_country_loso_diagnosis_HC3_delta_results.csv"
  )
)

safe_write_csv(
  wc_loso_diagnosis_summary,
  file.path(
    OUTDIR,
    "tables",
    "model_stability",
    "diagnosis",
    "within_country_loso_diagnosis_HC3_summary_by_module.csv"
  )
)

safe_write_csv(
  wc_loso_biomarker_long,
  file.path(
    OUTDIR,
    "tables",
    "model_stability",
    "biomarker",
    "within_country_loso_focal_log_biomarker_HC3_results.csv"
  )
)

safe_write_csv(
  wc_loso_biomarker_delta,
  file.path(
    OUTDIR,
    "tables",
    "model_stability",
    "biomarker",
    "within_country_loso_focal_log_biomarker_HC3_delta_results.csv"
  )
)

safe_write_csv(
  wc_loso_biomarker_summary,
  file.path(
    OUTDIR,
    "tables",
    "model_stability",
    "biomarker",
    "within_country_loso_focal_log_biomarker_HC3_summary.csv"
  )
)

###############################################################################
# 12) ALL-SITE DELETION AUDIT
###############################################################################

all_site_correlation_results <- list()

for (ii in seq_len(
  nrow(all_site_units)
)) {
  country_i <-
    all_site_units$
      Country[[ii]]

  site_i <-
    all_site_units$
      Site[[ii]]

  unit_i <-
    all_site_units$
      Country_Site[[ii]]

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

  run_id <- paste0(
    "ALL_SITE_DELETE_",
    gsub(
      "::",
      "_",
      unit_i
    )
  )

  all_site_correlation_results[[
    unit_i
  ]] <- compute_module_trait(
    df = df_run,
    modules = modules_use,
    traits = traits_available,
    analysis =
      "all_site_deletion_audit",
    run_id = run_id,
    excluded_unit = unit_i
  ) %>%
    dplyr::mutate(
      Excluded_country =
        country_i,
      Excluded_site = site_i,
      Country_remains_after_exclusion =
        all_site_units$
          country_remains_after_site_exclusion[
            ii
          ],
      Equivalent_to_country_removal =
        all_site_units$
          equivalent_to_country_removal[
            ii
          ]
    )
}

all_site_correlation_long <-
  dplyr::bind_rows(
    all_site_correlation_results
  )

all_site_correlation_delta <-
  attach_full_correlation_reference(
    sensitivity_tbl =
      all_site_correlation_long,
    full_tbl =
      full_module_trait
  )

all_site_summary_by_module <-
  summarize_correlation_delta_by_module(
    all_site_correlation_delta
  )

all_site_summary_by_module_trait <-
  summarize_correlation_delta_by_module_trait(
    all_site_correlation_delta
  )

all_site_summary_by_unit <-
  summarize_correlation_delta_by_run(
    all_site_correlation_delta
  ) %>%
  dplyr::left_join(
    all_site_units %>%
      dplyr::select(
        Excluded_unit =
          Country_Site,
        Country,
        Site,
        country_remains_after_site_exclusion,
        equivalent_to_country_removal
      ),
    by = "Excluded_unit"
  )

safe_write_csv(
  all_site_correlation_long,
  file.path(
    OUTDIR,
    "tables",
    "all_site_deletion_audit",
    "all_site_deletion_module_trait_results.csv"
  )
)

safe_write_csv(
  all_site_correlation_delta,
  file.path(
    OUTDIR,
    "tables",
    "all_site_deletion_audit",
    "all_site_deletion_delta_results.csv"
  )
)

safe_write_csv(
  all_site_summary_by_module,
  file.path(
    OUTDIR,
    "tables",
    "all_site_deletion_audit",
    "all_site_deletion_summary_by_module.csv"
  )
)

safe_write_csv(
  all_site_summary_by_module_trait,
  file.path(
    OUTDIR,
    "tables",
    "all_site_deletion_audit",
    "all_site_deletion_summary_by_module_trait.csv"
  )
)

safe_write_csv(
  all_site_summary_by_unit,
  file.path(
    OUTDIR,
    "tables",
    "all_site_deletion_audit",
    "all_site_deletion_summary_by_unit.csv"
  )
)

all_site_interpretation <- tibble::tibble(
  analysis = "all_site_deletion_audit",
  valid_primary_site_deletions =
    paste(
      primary_wc_loso_units$
        Country_Site,
      collapse = ", "
    ),
  site_deletions_equivalent_to_country_removal =
    paste(
      all_site_units$
        Country_Site[
          all_site_units$
            equivalent_to_country_removal
        ],
      collapse = ", "
    ),
  interpretation = paste(
    "All-site deletion is retained only as an audit.",
    "Deletions from countries with a single site are equivalent to",
    "leave-one-country-out and must not be interpreted as independent",
    "site robustness. Primary site sensitivity is the corrected",
    "within-country LOSO analysis."
  )
)

safe_write_csv(
  all_site_interpretation,
  file.path(
    OUTDIR,
    "tables",
    "all_site_deletion_audit",
    "all_site_deletion_interpretation.csv"
  )
)

###############################################################################
# 13) BALANCED DOWNSAMPLING
###############################################################################

set.seed(SEED)

downsampling_correlation_results <-
  vector(
    "list",
    N_DOWNSAMPLING_ITER
  )

downsampling_diagnosis_results <-
  vector(
    "list",
    N_DOWNSAMPLING_ITER
  )

downsampling_biomarker_results <-
  vector(
    "list",
    N_DOWNSAMPLING_ITER
  )

downsampling_sample_manifest <-
  vector(
    "list",
    N_DOWNSAMPLING_ITER
  )

for (iter in seq_len(
  N_DOWNSAMPLING_ITER
)) {
  sampled_df <- analysis_df %>%
    dplyr::group_by(
      Country,
      SampleGroup_bin
    ) %>%
    dplyr::slice_sample(
      n = min_stratum_n,
      replace = FALSE
    ) %>%
    dplyr::ungroup()

  sampled_df$Country <- droplevels(
    factor(
      sampled_df$Country
    )
  )

  run_id <- paste0(
    "DOWNSAMPLE_",
    sprintf(
      "%03d",
      iter
    )
  )

  downsampling_sample_manifest[[
    iter
  ]] <- sampled_df %>%
    dplyr::transmute(
      Iteration = iter,
      Run_ID = run_id,
      SampleId,
      Country =
        as.character(Country),
      SampleGroup,
      SampleGroup_bin
    )

  downsampling_correlation_results[[
    iter
  ]] <- compute_module_trait(
    df = sampled_df,
    modules = modules_use,
    traits = traits_available,
    analysis =
      "balanced_downsampling",
    run_id = run_id,
    iteration = iter
  ) %>%
    dplyr::mutate(
      n_per_country_diagnosis =
        min_stratum_n
    )

  downsampling_diagnosis_results[[
    iter
  ]] <- run_diagnosis_hc3_models(
    df = sampled_df,
    modules = modules_use,
    analysis =
      "balanced_downsampling",
    run_id = run_id,
    iteration = iter
  ) %>%
    dplyr::mutate(
      n_per_country_diagnosis =
        min_stratum_n
    )

  downsampling_biomarker_results[[
    iter
  ]] <-
    run_log_biomarker_hc3_models(
      df = sampled_df,
      modules =
        focal_biomarker_modules_use,
      biomarkers = BIOMARKERS,
      analysis =
        "balanced_downsampling",
      run_id = run_id,
      iteration = iter
    ) %>%
    dplyr::mutate(
      n_per_country_diagnosis =
        min_stratum_n
    )

  if (
    iter %% PROGRESS_EVERY == 0 ||
    iter ==
      N_DOWNSAMPLING_ITER
  ) {
    cat(
      "Balanced downsampling completed: ",
      iter,
      "/",
      N_DOWNSAMPLING_ITER,
      "\n",
      sep = ""
    )
  }
}

downsampling_sample_manifest_tbl <-
  dplyr::bind_rows(
    downsampling_sample_manifest
  )

downsampling_correlation_long <-
  dplyr::bind_rows(
    downsampling_correlation_results
  )

downsampling_correlation_delta <-
  attach_full_correlation_reference(
    sensitivity_tbl =
      downsampling_correlation_long,
    full_tbl =
      full_module_trait
  )

downsampling_summary_by_module <-
  summarize_correlation_delta_by_module(
    downsampling_correlation_delta
  )

downsampling_summary_by_module_trait <-
  summarize_correlation_delta_by_module_trait(
    downsampling_correlation_delta
  )

downsampling_summary_by_iteration <-
  summarize_correlation_delta_by_run(
    downsampling_correlation_delta
  )

downsampling_diagnosis_long <-
  dplyr::bind_rows(
    downsampling_diagnosis_results
  )

downsampling_diagnosis_delta <-
  attach_full_model_reference(
    sensitivity_tbl =
      downsampling_diagnosis_long,
    full_tbl =
      full_diagnosis_models,
    keys = "Module",
    effect_col =
      "standardized_difference",
    full_effect_col =
      "standardized_difference_full"
  )

downsampling_diagnosis_summary <-
  summarize_model_stability(
    delta_tbl =
      downsampling_diagnosis_delta,
    group_vars = "Module",
    effect_col =
      "standardized_difference",
    full_effect_col =
      "standardized_difference_full"
  )

downsampling_biomarker_long <-
  dplyr::bind_rows(
    downsampling_biomarker_results
  )

downsampling_biomarker_delta <-
  attach_full_model_reference(
    sensitivity_tbl =
      downsampling_biomarker_long,
    full_tbl =
      full_biomarker_models,
    keys = c(
      "Module",
      "Biomarker"
    ),
    effect_col =
      "standardized_beta",
    full_effect_col =
      "standardized_beta_full"
  )

downsampling_biomarker_summary <-
  summarize_model_stability(
    delta_tbl =
      downsampling_biomarker_delta,
    group_vars = c(
      "Module",
      "Biomarker"
    ),
    effect_col =
      "standardized_beta",
    full_effect_col =
      "standardized_beta_full"
  )

expected_downsampling_corr_rows <-
  N_DOWNSAMPLING_ITER *
  length(modules_use) *
  length(traits_available)

expected_downsampling_diagnosis_rows <-
  N_DOWNSAMPLING_ITER *
  length(modules_use)

expected_downsampling_biomarker_rows <-
  N_DOWNSAMPLING_ITER *
  length(
    focal_biomarker_modules_use
  ) *
  length(BIOMARKERS)

if (
  nrow(
    downsampling_correlation_long
  ) !=
    expected_downsampling_corr_rows
) {
  stop(
    "Balanced-downsampling correlation row count mismatch.",
    call. = FALSE
  )
}

if (
  nrow(
    downsampling_diagnosis_long
  ) !=
    expected_downsampling_diagnosis_rows
) {
  stop(
    "Balanced-downsampling diagnosis-model row count mismatch.",
    call. = FALSE
  )
}

if (
  nrow(
    downsampling_biomarker_long
  ) !=
    expected_downsampling_biomarker_rows
) {
  stop(
    "Balanced-downsampling biomarker-model row count mismatch.",
    call. = FALSE
  )
}

safe_write_csv(
  downsampling_sample_manifest_tbl,
  file.path(
    OUTDIR,
    "tables",
    "downsampling",
    "balanced_downsampling_sample_manifest.csv"
  )
)

safe_write_csv(
  downsampling_correlation_long,
  file.path(
    OUTDIR,
    "tables",
    "downsampling",
    "balanced_downsampling_module_trait_results.csv"
  )
)

safe_write_csv(
  downsampling_correlation_delta,
  file.path(
    OUTDIR,
    "tables",
    "downsampling",
    "balanced_downsampling_delta_results.csv"
  )
)

safe_write_csv(
  downsampling_summary_by_module,
  file.path(
    OUTDIR,
    "tables",
    "downsampling",
    "balanced_downsampling_summary_by_module.csv"
  )
)

safe_write_csv(
  downsampling_summary_by_module_trait,
  file.path(
    OUTDIR,
    "tables",
    "downsampling",
    "balanced_downsampling_summary_by_module_trait.csv"
  )
)

safe_write_csv(
  downsampling_summary_by_iteration,
  file.path(
    OUTDIR,
    "tables",
    "downsampling",
    "balanced_downsampling_summary_by_iteration.csv"
  )
)

safe_write_csv(
  downsampling_diagnosis_long,
  file.path(
    OUTDIR,
    "tables",
    "model_stability",
    "diagnosis",
    "balanced_downsampling_diagnosis_HC3_results.csv"
  )
)

safe_write_csv(
  downsampling_diagnosis_delta,
  file.path(
    OUTDIR,
    "tables",
    "model_stability",
    "diagnosis",
    "balanced_downsampling_diagnosis_HC3_delta_results.csv"
  )
)

safe_write_csv(
  downsampling_diagnosis_summary,
  file.path(
    OUTDIR,
    "tables",
    "model_stability",
    "diagnosis",
    "balanced_downsampling_diagnosis_HC3_summary_by_module.csv"
  )
)

safe_write_csv(
  downsampling_biomarker_long,
  file.path(
    OUTDIR,
    "tables",
    "model_stability",
    "biomarker",
    "balanced_downsampling_focal_log_biomarker_HC3_results.csv"
  )
)

safe_write_csv(
  downsampling_biomarker_delta,
  file.path(
    OUTDIR,
    "tables",
    "model_stability",
    "biomarker",
    "balanced_downsampling_focal_log_biomarker_HC3_delta_results.csv"
  )
)

safe_write_csv(
  downsampling_biomarker_summary,
  file.path(
    OUTDIR,
    "tables",
    "model_stability",
    "biomarker",
    "balanced_downsampling_focal_log_biomarker_HC3_summary.csv"
  )
)

###############################################################################
# 14) CORE-MODULE CORRELATION STABILITY SUMMARY
###############################################################################

core_loco <- loco_summary_by_module %>%
  dplyr::filter(
    Module %in%
      core_modules_use
  ) %>%
  dplyr::mutate(
    Analysis =
      "LOCO country"
  )

core_wc_loso <-
  wc_loso_summary_by_module %>%
  dplyr::filter(
    Module %in%
      core_modules_use
  ) %>%
  dplyr::mutate(
    Analysis =
      "WC-LOSO primary"
  )

core_downsampling <-
  downsampling_summary_by_module %>%
  dplyr::filter(
    Module %in%
      core_modules_use
  ) %>%
  dplyr::mutate(
    Analysis =
      "Balanced downsampling"
  )

core_correlation_stability <-
  dplyr::bind_rows(
    core_loco,
    core_wc_loso,
    core_downsampling
  ) %>%
  dplyr::left_join(
    module_reference,
    by = "Module"
  ) %>%
  dplyr::arrange(
    Analysis,
    mean_abs_delta_rho
  )

safe_write_csv(
  core_correlation_stability,
  file.path(
    OUTDIR,
    "tables",
    "core_modules",
    "core_modules_correlation_stability_summary.csv"
  )
)

core_module_trait_stability <-
  dplyr::bind_rows(
    loco_summary_by_module_trait %>%
      dplyr::mutate(
        Analysis =
          "LOCO country"
      ),
    wc_loso_summary_by_module_trait %>%
      dplyr::mutate(
        Analysis =
          "WC-LOSO primary"
      ),
    downsampling_summary_by_module_trait %>%
      dplyr::mutate(
        Analysis =
          "Balanced downsampling"
      )
  ) %>%
  dplyr::filter(
    Module %in%
      core_modules_use
  ) %>%
  dplyr::mutate(
    Trait_label =
      dplyr::recode(
        Trait,
        !!!TRAIT_LABELS,
        .default = Trait
      )
  ) %>%
  dplyr::left_join(
    module_reference,
    by = "Module"
  )

safe_write_csv(
  core_module_trait_stability,
  file.path(
    OUTDIR,
    "tables",
    "core_modules",
    "core_modules_trait_level_stability_summary.csv"
  )
)

###############################################################################
# 15) INTEGRATED DIAGNOSIS AND BIOMARKER MODEL STABILITY
###############################################################################

diagnosis_stability_integrated <-
  module_reference %>%
  dplyr::left_join(
    full_diagnosis_models %>%
      dplyr::select(
        Module,
        full_standardized_difference =
          standardized_difference,
        full_robust_p =
          robust_p_value,
        full_FDR =
          FDR_HC3
      ),
    by = "Module"
  ) %>%
  dplyr::left_join(
    loco_diagnosis_summary %>%
      dplyr::select(
        Module,
        loco_direction_consistency =
          direction_consistency,
        loco_mean_abs_delta =
          mean_abs_delta_effect,
        loco_max_abs_delta =
          max_abs_delta_effect,
        loco_nominal_rate =
          nominal_significance_rate,
        loco_FDR_rate =
          FDR_significance_rate,
        loco_effect_2_5 =
          empirical_effect_2_5,
        loco_effect_97_5 =
          empirical_effect_97_5
      ),
    by = "Module"
  ) %>%
  dplyr::left_join(
    wc_loso_diagnosis_summary %>%
      dplyr::select(
        Module,
        wc_loso_direction_consistency =
          direction_consistency,
        wc_loso_mean_abs_delta =
          mean_abs_delta_effect,
        wc_loso_max_abs_delta =
          max_abs_delta_effect,
        wc_loso_nominal_rate =
          nominal_significance_rate,
        wc_loso_FDR_rate =
          FDR_significance_rate,
        wc_loso_effect_2_5 =
          empirical_effect_2_5,
        wc_loso_effect_97_5 =
          empirical_effect_97_5
      ),
    by = "Module"
  ) %>%
  dplyr::left_join(
    downsampling_diagnosis_summary %>%
      dplyr::select(
        Module,
        downsampling_direction_consistency =
          direction_consistency,
        downsampling_mean_abs_delta =
          mean_abs_delta_effect,
        downsampling_max_abs_delta =
          max_abs_delta_effect,
        downsampling_nominal_rate =
          nominal_significance_rate,
        downsampling_FDR_rate =
          FDR_significance_rate,
        downsampling_effect_2_5 =
          empirical_effect_2_5,
        downsampling_effect_97_5 =
          empirical_effect_97_5
      ),
    by = "Module"
  ) %>%
  dplyr::arrange(
    full_FDR,
    dplyr::desc(
      abs(
        full_standardized_difference
      )
    )
  )

safe_write_csv(
  diagnosis_stability_integrated,
  file.path(
    OUTDIR,
    "tables",
    "model_stability",
    "diagnosis",
    "integrated_diagnosis_HC3_stability_summary.csv"
  )
)

biomarker_stability_integrated <-
  full_biomarker_models %>%
  dplyr::select(
    Module,
    Biomarker,
    full_standardized_beta =
      standardized_beta,
    full_percent_change_per_1SD =
      percent_change_per_1SD_module,
    full_robust_p =
      robust_p_value,
    full_FDR =
      FDR_HC3
  ) %>%
  dplyr::left_join(
    loco_biomarker_summary %>%
      dplyr::select(
        Module,
        Biomarker,
        loco_direction_consistency =
          direction_consistency,
        loco_mean_abs_delta =
          mean_abs_delta_effect,
        loco_max_abs_delta =
          max_abs_delta_effect,
        loco_nominal_rate =
          nominal_significance_rate,
        loco_FDR_rate =
          FDR_significance_rate,
        loco_beta_2_5 =
          empirical_effect_2_5,
        loco_beta_97_5 =
          empirical_effect_97_5
      ),
    by = c(
      "Module",
      "Biomarker"
    )
  ) %>%
  dplyr::left_join(
    wc_loso_biomarker_summary %>%
      dplyr::select(
        Module,
        Biomarker,
        wc_loso_direction_consistency =
          direction_consistency,
        wc_loso_mean_abs_delta =
          mean_abs_delta_effect,
        wc_loso_max_abs_delta =
          max_abs_delta_effect,
        wc_loso_nominal_rate =
          nominal_significance_rate,
        wc_loso_FDR_rate =
          FDR_significance_rate,
        wc_loso_beta_2_5 =
          empirical_effect_2_5,
        wc_loso_beta_97_5 =
          empirical_effect_97_5
      ),
    by = c(
      "Module",
      "Biomarker"
    )
  ) %>%
  dplyr::left_join(
    downsampling_biomarker_summary %>%
      dplyr::select(
        Module,
        Biomarker,
        downsampling_direction_consistency =
          direction_consistency,
        downsampling_mean_abs_delta =
          mean_abs_delta_effect,
        downsampling_max_abs_delta =
          max_abs_delta_effect,
        downsampling_nominal_rate =
          nominal_significance_rate,
        downsampling_FDR_rate =
          FDR_significance_rate,
        downsampling_beta_2_5 =
          empirical_effect_2_5,
        downsampling_beta_97_5 =
          empirical_effect_97_5
      ),
    by = c(
      "Module",
      "Biomarker"
    )
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
  ) %>%
  dplyr::arrange(
    Module,
    full_FDR
  )

safe_write_csv(
  biomarker_stability_integrated,
  file.path(
    OUTDIR,
    "tables",
    "model_stability",
    "biomarker",
    "integrated_focal_log_biomarker_HC3_stability_summary.csv"
  )
)

###############################################################################
# 16) FIGURES — MODULE-LEVEL CORRELATION STABILITY
###############################################################################

plot_module_delta <- function(
    summary_tbl,
    title,
    subtitle,
    file_stub,
    output_subdir
) {
  p <- summary_tbl %>%
    dplyr::mutate(
      Module = factor(
        Module,
        levels = Module[
          order(
            mean_abs_delta_rho,
            decreasing = TRUE
          )
        ]
      )
    ) %>%
    ggplot(
      aes(
        x = Module,
        y = mean_abs_delta_rho,
        fill = Module
      )
    ) +
    geom_col(
      show.legend = FALSE
    ) +
    coord_flip() +
    scale_fill_manual(
      values = get_module_colors(
        summary_tbl$Module
      )
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = NULL,
      y = "Mean absolute change in Spearman rho"
    ) +
    theme_bw(
      base_size = 12
    ) +
    theme(
      panel.grid.major.y =
        element_blank(),
      plot.title = element_text(
        face = "bold"
      )
    )

  ggsave(
    file.path(
      OUTDIR,
      "figures",
      output_subdir,
      paste0(
        file_stub,
        ".pdf"
      )
    ),
    p,
    width = 7,
    height = 5
  )

  ggsave(
    file.path(
      OUTDIR,
      "figures",
      output_subdir,
      paste0(
        file_stub,
        ".png"
      )
    ),
    p,
    width = 7,
    height = 5,
    dpi = DPI
  )

  p
}

plot_module_delta(
  summary_tbl =
    loco_summary_by_module,
  title =
    "Leave-one-country-out association stability",
  subtitle =
    "All 8 modules and 16 traits; module definitions held fixed",
  file_stub =
    "loco_country_mean_abs_delta_rho",
  output_subdir = "loco"
)

plot_module_delta(
  summary_tbl =
    wc_loso_summary_by_module,
  title =
    "Corrected within-country leave-one-site-out stability",
  subtitle = paste0(
    "Primary site-deletion analysis: ",
    paste(
      primary_wc_loso_units$
        Country_Site,
      collapse = ", "
    )
  ),
  file_stub =
    "within_country_loso_mean_abs_delta_rho",
  output_subdir =
    "within_country_loso"
)

plot_module_delta(
  summary_tbl =
    downsampling_summary_by_module,
  title =
    "Country-diagnosis balanced downsampling stability",
  subtitle = paste0(
    N_DOWNSAMPLING_ITER,
    " iterations; ",
    min_stratum_n,
    " participants per Country × diagnosis stratum; total n = ",
    downsampling_total_n
  ),
  file_stub =
    "balanced_downsampling_mean_abs_delta_rho",
  output_subdir =
    "downsampling"
)

###############################################################################
# 17) FIGURES — CORE MODULES
###############################################################################

p_core <- core_correlation_stability %>%
  dplyr::mutate(
    Analysis = factor(
      Analysis,
      levels = c(
        "LOCO country",
        "WC-LOSO primary",
        "Balanced downsampling"
      )
    )
  ) %>%
  ggplot(
    aes(
      x = Module,
      y = mean_abs_delta_rho,
      fill = Module
    )
  ) +
  geom_col(
    show.legend = FALSE
  ) +
  facet_wrap(
    ~ Analysis,
    scales = "free_y"
  ) +
  scale_fill_manual(
    values = get_module_colors(
      core_modules_use
    )
  ) +
  labs(
    title =
      "Core-module association stability",
    subtitle = paste0(
      "Blue, green and brown modules; lower absolute change indicates ",
      "greater effect-size stability"
    ),
    x = NULL,
    y = "Mean absolute change in Spearman rho"
  ) +
  theme_bw(
    base_size = 12
  ) +
  theme(
    panel.grid.major.x =
      element_blank(),
    strip.text = element_text(
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
    "core_modules",
    "core_modules_correlation_stability.pdf"
  ),
  p_core,
  width = 10,
  height = 4.5
)

ggsave(
  file.path(
    OUTDIR,
    "figures",
    "core_modules",
    "core_modules_correlation_stability.png"
  ),
  p_core,
  width = 10,
  height = 4.5,
  dpi = DPI
)

core_heatmap_df <- core_module_trait_stability %>%
  dplyr::mutate(
    Module = factor(
      Module,
      levels = rev(
        core_modules_use
      )
    ),
    Trait_label = factor(
      Trait_label,
      levels = unname(
        dplyr::recode(
          traits_available,
          !!!TRAIT_LABELS,
          .default =
            traits_available
        )
      )
    )
  )

for (analysis_name in unique(
  core_heatmap_df$Analysis
)) {
  plot_df <- core_heatmap_df %>%
    dplyr::filter(
      Analysis ==
        analysis_name
    )

  p_heat <- ggplot(
    plot_df,
    aes(
      x = Trait_label,
      y = Module,
      fill = mean_abs_delta_rho
    )
  ) +
    geom_tile(
      colour = "white",
      linewidth = 0.4
    ) +
    geom_text(
      aes(
        label = sprintf(
          "%.02f",
          direction_consistency
        )
      ),
      size = 2.8
    ) +
    scale_fill_gradient(
      low = "white",
      high = "#B2182B",
      name = "Mean |Δrho|"
    ) +
    labs(
      title = paste0(
        "Core module-trait stability — ",
        analysis_name
      ),
      subtitle =
        "Cell text indicates direction-consistency proportion",
      x = NULL,
      y = NULL
    ) +
    theme_bw(
      base_size = 10
    ) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
      axis.text.y = element_text(
        face = "bold"
      ),
      plot.title = element_text(
        face = "bold"
      )
    )

  file_tag <- gsub(
    "[^A-Za-z0-9]+",
    "_",
    tolower(
      analysis_name
    )
  )

  ggsave(
    file.path(
      OUTDIR,
      "figures",
      "core_modules",
      paste0(
        "core_module_trait_stability_",
        file_tag,
        ".pdf"
      )
    ),
    p_heat,
    width = 12,
    height = 4.5
  )

  ggsave(
    file.path(
      OUTDIR,
      "figures",
      "core_modules",
      paste0(
        "core_module_trait_stability_",
        file_tag,
        ".png"
      )
    ),
    p_heat,
    width = 12,
    height = 4.5,
    dpi = DPI
  )
}

###############################################################################
# 18) FIGURES — DIAGNOSIS MODEL STABILITY
###############################################################################

diagnosis_plot_df <- diagnosis_stability_integrated %>%
  dplyr::filter(
    Module %in%
      modules_use
  ) %>%
  dplyr::mutate(
    Module = factor(
      Module,
      levels = rev(
        Module[
          order(
            full_standardized_difference,
            na.last = TRUE
          )
        ]
      )
    )
  )

p_diag_stability <- ggplot(
  diagnosis_plot_df,
  aes(
    y = Module,
    colour = Module
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = 2,
    colour = "grey55"
  ) +
  geom_errorbarh(
    aes(
      xmin =
        downsampling_effect_2_5,
      xmax =
        downsampling_effect_97_5
    ),
    height = 0.18,
    linewidth = 0.8,
    alpha = 0.75
  ) +
  geom_point(
    aes(
      x =
        full_standardized_difference
    ),
    size = 3
  ) +
  scale_colour_manual(
    values = get_module_colors(
      modules_use
    )
  ) +
  labs(
    title =
      "Adjusted diagnosis-effect stability",
    subtitle = paste0(
      "Points: full HC3-adjusted effect; horizontal intervals: ",
      "2.5th–97.5th percentiles across balanced downsampling"
    ),
    x =
      "Diagnosis effect on module eigengene (SD units)",
    y = NULL
  ) +
  theme_bw(
    base_size = 12
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(
      face = "bold"
    )
  )

ggsave(
  file.path(
    OUTDIR,
    "figures",
    "model_stability",
    "diagnosis_HC3_downsampling_stability.pdf"
  ),
  p_diag_stability,
  width = 8,
  height = 5.5
)

ggsave(
  file.path(
    OUTDIR,
    "figures",
    "model_stability",
    "diagnosis_HC3_downsampling_stability.png"
  ),
  p_diag_stability,
  width = 8,
  height = 5.5,
  dpi = DPI
)

###############################################################################
# 19) FIGURES — FOCAL BIOMARKER MODEL STABILITY
###############################################################################

biomarker_direction_plot <- biomarker_stability_integrated %>%
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
    Module_Biomarker = paste(
      Module,
      Biomarker_label,
      sep = " — "
    ),
    Module_Biomarker =
      factor(
        Module_Biomarker,
        levels = rev(
          unique(
            Module_Biomarker
          )
        )
      ),
    Analysis = factor(
      Analysis,
      levels = c(
        "LOCO country",
        "WC-LOSO primary",
        "Balanced downsampling"
      )
    )
  )

p_biomarker_stability <- ggplot(
  biomarker_direction_plot,
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
      "Focal log-HC3 biomarker association stability",
    subtitle =
      "Blue, green and brown modules across recruitment-context and balancing perturbations",
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
    "model_stability",
    "focal_log_HC3_biomarker_direction_stability.pdf"
  ),
  p_biomarker_stability,
  width = 8,
  height = 7
)

ggsave(
  file.path(
    OUTDIR,
    "figures",
    "model_stability",
    "focal_log_HC3_biomarker_direction_stability.png"
  ),
  p_biomarker_stability,
  width = 8,
  height = 7,
  dpi = DPI
)

###############################################################################
# 20) MANUSCRIPT-FACING METHODS AND INTERPRETATION
###############################################################################

methods_wording <- tibble::tribble(
  ~section,
  ~text,

  "Methods - association stability",
  paste(
    "Module definitions and eigengenes were held fixed while module-trait",
    "associations were re-estimated under leave-one-country-out,",
    "within-country leave-one-site-out and balanced-downsampling perturbations.",
    "For each perturbation, Spearman correlations were compared with the",
    "full-sample estimates using absolute changes in correlation magnitude,",
    "direction consistency and preservation of nominal significance."
  ),

  "Methods - LOCO",
  paste(
    "For leave-one-country-out analyses, each recruitment country was excluded",
    "in turn and the complete module-trait matrix was recomputed. False-discovery",
    "rates were recalculated within each leave-one-country-out run."
  ),

  "Methods - corrected WC-LOSO",
  paste(
    "Because recruitment site was nested within country, the primary",
    "leave-one-site-out analysis was restricted to site deletions within",
    "countries containing at least two sites. Thus, the corresponding country",
    "remained represented after each site deletion. Deletions of sites from",
    "single-site countries were retained only as an audit because they are",
    "equivalent to leave-one-country-out analyses."
  ),

  "Methods - balanced downsampling",
  paste(
    "Country-diagnosis imbalance was evaluated using 500 reproducible",
    "downsampling iterations. In each iteration, the same number of participants",
    "was sampled without replacement from every country-by-diagnosis stratum,",
    "using the smallest observed stratum size. Module-trait associations were",
    "recomputed in each balanced sample."
  ),

  "Methods - model-based stability",
  paste(
    "Adjusted diagnosis-module models and prespecified log-transformed plasma",
    "biomarker models were additionally repeated using HC3 robust standard",
    "errors. Diagnosis models were evaluated across all eight modules, whereas",
    "biomarker-model stability focused on the prespecified blue, green and brown",
    "modules for p-tau181, p-tau217, NfL and Aβ42/40."
  ),

  "Interpretation",
  paste(
    "These analyses quantify stability of module associations with fixed module",
    "definitions. They do not evaluate structural network preservation and do",
    "not constitute external validation. Loss of nominal significance in",
    "balanced subsamples was interpreted jointly with effect direction and",
    "magnitude because downsampling intentionally reduced statistical power."
  )
)

safe_write_csv(
  methods_wording,
  file.path(
    OUTDIR,
    "tables",
    "script14_methods_and_interpretation_wording.csv"
  )
)

###############################################################################
# 21) EXCEL WORKBOOK
###############################################################################

workbook_tables <- list(
  Input_audit =
    input_audit,
  Full_reference_audit =
    correlation_reference_summary,
  Full_correlations =
    full_module_trait,
  LOCO_by_module =
    loco_summary_by_module,
  LOCO_by_module_trait =
    loco_summary_by_module_trait,
  LOCO_by_country =
    loco_summary_by_excluded_country,
  WC_LOSO_by_module =
    wc_loso_summary_by_module,
  WC_LOSO_by_module_trait =
    wc_loso_summary_by_module_trait,
  WC_LOSO_by_site =
    wc_loso_summary_by_excluded_site,
  All_site_audit =
    all_site_summary_by_unit,
  Downsampling_by_module =
    downsampling_summary_by_module,
  Downsampling_by_trait =
    downsampling_summary_by_module_trait,
  Downsampling_by_iteration =
    downsampling_summary_by_iteration,
  Core_modules =
    core_correlation_stability,
  Diagnosis_stability =
    diagnosis_stability_integrated,
  Biomarker_stability =
    biomarker_stability_integrated,
  Methods_wording =
    methods_wording
)

openxlsx::write.xlsx(
  workbook_tables,
  file = file.path(
    OUTDIR,
    "WGCNA_14_Association_Stability_LOCO_LOSO_Downsampling.xlsx"
  ),
  overwrite = TRUE
)

###############################################################################
# 22) FINAL SUMMARY
###############################################################################

script14_summary <- tibble::tibble(
  metric = c(
    "base_dir",
    "output_dir",
    "wgcna_input_level",
    "n_samples_full",
    "n_modules",
    "modules",
    "n_traits",
    "n_full_module_trait_tests",
    "full_reference_match_passed",
    "n_countries_LOCO",
    "n_primary_WC_LOSO_units",
    "primary_WC_LOSO_units",
    "n_all_site_audit_units",
    "n_downsampling_iterations_requested",
    "n_downsampling_iterations_completed",
    "minimum_country_diagnosis_stratum_n",
    "downsampling_total_n_per_iteration",
    "downsampling_correlation_rows",
    "downsampling_diagnosis_model_rows",
    "downsampling_biomarker_model_rows",
    "core_modules",
    "focal_biomarker_modules",
    "random_seed",
    "Country_numeric_included"
  ),
  value = c(
    BASE_DIR,
    OUTDIR,
    "GENE-COLLAPSED, outcome-independent SOMAmer selection",
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
      length(traits_available)
    ),
    as.character(
      nrow(full_module_trait)
    ),
    as.character(
      reference_match_passed
    ),
    as.character(
      length(countries)
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
      nrow(all_site_units)
    ),
    as.character(
      N_DOWNSAMPLING_ITER
    ),
    as.character(
      dplyr::n_distinct(
        downsampling_correlation_long$
          Iteration
      )
    ),
    as.character(
      min_stratum_n
    ),
    as.character(
      downsampling_total_n
    ),
    as.character(
      nrow(
        downsampling_correlation_long
      )
    ),
    as.character(
      nrow(
        downsampling_diagnosis_long
      )
    ),
    as.character(
      nrow(
        downsampling_biomarker_long
      )
    ),
    paste(
      core_modules_use,
      collapse = ", "
    ),
    paste(
      focal_biomarker_modules_use,
      collapse = ", "
    ),
    as.character(SEED),
    as.character(
      "Country_numeric" %in%
        names(analysis_df)
    )
  )
)

safe_write_csv(
  script14_summary,
  file.path(
    OUTDIR,
    "tables",
    "script14_final_summary.csv"
  )
)

output_manifest <- tibble::tibble(
  output_file = c(
    "tables/script14_input_alignment_audit.csv",
    "tables/full/full_module_trait_correlations_recomputed.csv",
    "tables/full/full_module_trait_reference_audit.csv",
    "tables/loco/loco_country_module_trait_results.csv",
    "tables/loco/loco_country_delta_results.csv",
    "tables/loco/loco_country_summary_by_module.csv",
    "tables/loco/loco_country_summary_by_module_trait.csv",
    "tables/loco/loco_country_summary_by_excluded_country.csv",
    "tables/within_country_loso/within_country_loso_module_trait_results.csv",
    "tables/within_country_loso/within_country_loso_delta_results.csv",
    "tables/within_country_loso/within_country_loso_summary_by_module.csv",
    "tables/within_country_loso/within_country_loso_summary_by_module_trait.csv",
    "tables/within_country_loso/within_country_loso_summary_by_excluded_site.csv",
    "tables/all_site_deletion_audit/all_site_deletion_*.csv",
    "tables/downsampling/balanced_downsampling_sample_manifest.csv",
    "tables/downsampling/balanced_downsampling_module_trait_results.csv",
    "tables/downsampling/balanced_downsampling_delta_results.csv",
    "tables/downsampling/balanced_downsampling_summary_by_module.csv",
    "tables/downsampling/balanced_downsampling_summary_by_module_trait.csv",
    "tables/downsampling/balanced_downsampling_summary_by_iteration.csv",
    "tables/model_stability/diagnosis/*.csv",
    "tables/model_stability/biomarker/*.csv",
    "tables/core_modules/core_modules_correlation_stability_summary.csv",
    "tables/core_modules/core_modules_trait_level_stability_summary.csv",
    "tables/script14_methods_and_interpretation_wording.csv",
    "figures/loco/loco_country_mean_abs_delta_rho.pdf/png",
    "figures/within_country_loso/within_country_loso_mean_abs_delta_rho.pdf/png",
    "figures/downsampling/balanced_downsampling_mean_abs_delta_rho.pdf/png",
    "figures/core_modules/core_modules_correlation_stability.pdf/png",
    "figures/model_stability/diagnosis_HC3_downsampling_stability.pdf/png",
    "figures/model_stability/focal_log_HC3_biomarker_direction_stability.pdf/png",
    "WGCNA_14_Association_Stability_LOCO_LOSO_Downsampling.xlsx",
    "workspace/script14_association_stability_workspace.RData",
    "sessionInfo.txt"
  ),
  description = c(
    "Strict input, country, site and downsampling-structure audit.",
    "Independently recomputed full module-trait matrix.",
    "Exact comparison with the Script 13 full matrix.",
    "All LOCO module-trait estimates.",
    "LOCO estimates joined to full effects and stability metrics.",
    "LOCO stability summarized by module.",
    "LOCO stability summarized by module-trait pair.",
    "Overall impact of excluding each country.",
    "Primary corrected within-country site-deletion estimates.",
    "WC-LOSO estimates joined to full effects.",
    "WC-LOSO stability summarized by module.",
    "WC-LOSO stability summarized by module-trait pair.",
    "Overall impact of each valid within-country site deletion.",
    "Secondary all-site audit with single-site-country warnings.",
    "Exact participant membership in all 500 balanced samples.",
    "All balanced-downsampling module-trait estimates.",
    "Downsampling estimates joined to full effects.",
    "Downsampling stability summarized by module.",
    "Downsampling empirical intervals and stability by module-trait.",
    "Overall perturbation magnitude by iteration.",
    "Full, LOCO, WC-LOSO and downsampling HC3 diagnosis stability.",
    "Focal blue/green/brown log-HC3 biomarker stability.",
    "Integrated correlation stability for blue, green and brown.",
    "Trait-level core-module stability across all perturbations.",
    "Methods-ready wording and interpretation boundaries.",
    "Module-level LOCO stability figure.",
    "Corrected within-country LOSO stability figure.",
    "Balanced-downsampling stability figure.",
    "Core-module comparison across sensitivity families.",
    "Full versus downsampling-adjusted diagnosis effects.",
    "Direction consistency of focal robust biomarker models.",
    "Integrated reviewer-facing workbook.",
    "Workspace for subsequent structural preservation and figures.",
    "R session information."
  )
)

safe_write_csv(
  output_manifest,
  file.path(
    OUTDIR,
    "script14_output_manifest.csv"
  )
)

###############################################################################
# 23) SAVE WORKSPACE
###############################################################################

save(
  analysis_df,
  module_reference,
  prioritization_reference,
  modules_use,
  traits_available,
  core_modules_use,
  focal_biomarker_modules_use,
  site_var,
  country_site_map,
  primary_wc_loso_units,
  all_site_units,
  strata_counts,
  min_stratum_n,
  downsampling_total_n,
  full_module_trait,
  full_correlation_audit,
  full_diagnosis_models,
  full_biomarker_models,
  loco_correlation_long,
  loco_correlation_delta,
  loco_summary_by_module,
  loco_summary_by_module_trait,
  loco_summary_by_excluded_country,
  loco_diagnosis_long,
  loco_diagnosis_delta,
  loco_diagnosis_summary,
  loco_biomarker_long,
  loco_biomarker_delta,
  loco_biomarker_summary,
  wc_loso_correlation_long,
  wc_loso_correlation_delta,
  wc_loso_summary_by_module,
  wc_loso_summary_by_module_trait,
  wc_loso_summary_by_excluded_site,
  wc_loso_diagnosis_long,
  wc_loso_diagnosis_delta,
  wc_loso_diagnosis_summary,
  wc_loso_biomarker_long,
  wc_loso_biomarker_delta,
  wc_loso_biomarker_summary,
  all_site_correlation_long,
  all_site_correlation_delta,
  all_site_summary_by_module,
  all_site_summary_by_module_trait,
  all_site_summary_by_unit,
  downsampling_sample_manifest_tbl,
  downsampling_correlation_long,
  downsampling_correlation_delta,
  downsampling_summary_by_module,
  downsampling_summary_by_module_trait,
  downsampling_summary_by_iteration,
  downsampling_diagnosis_long,
  downsampling_diagnosis_delta,
  downsampling_diagnosis_summary,
  downsampling_biomarker_long,
  downsampling_biomarker_delta,
  downsampling_biomarker_summary,
  core_correlation_stability,
  core_module_trait_stability,
  diagnosis_stability_integrated,
  biomarker_stability_integrated,
  methods_wording,
  script14_summary,
  output_manifest,
  file = file.path(
    OUTDIR,
    "workspace",
    "script14_association_stability_workspace.RData"
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
# 24) FINAL MESSAGE
###############################################################################

cat("\nScript 14 finished successfully.\n")
cat("Main output directory:\n", OUTDIR, "\n")
cat("Input level: GENE-COLLAPSED, outcome-independent SOMAmer selection.\n")
cat("Samples in full analysis: ", nrow(analysis_df), "\n", sep = "")
cat("Modules: ", length(modules_use), "\n", sep = "")
cat("Traits: ", length(traits_available), "\n", sep = "")
cat("Full module-trait tests: ", nrow(full_module_trait), "\n", sep = "")
cat("Full Script 13 reference match: ", reference_match_passed, "\n", sep = "")
cat("LOCO countries: ", length(countries), "\n", sep = "")
cat(
  "Primary WC-LOSO units: ",
  nrow(primary_wc_loso_units),
  " (",
  paste(
    primary_wc_loso_units$
      Country_Site,
    collapse = ", "
  ),
  ")\n",
  sep = ""
)
cat("All-site deletion audit units: ", nrow(all_site_units), "\n", sep = "")
cat("Balanced-downsampling iterations: ", N_DOWNSAMPLING_ITER, "\n", sep = "")
cat("Participants per balanced iteration: ", downsampling_total_n, "\n", sep = "")
cat("Country_numeric included: FALSE\n")
cat("\nKey outputs:\n")
cat("- tables/loco/loco_country_summary_by_module.csv\n")
cat("- tables/loco/loco_country_summary_by_module_trait.csv\n")
cat("- tables/within_country_loso/within_country_loso_summary_by_module.csv\n")
cat("- tables/within_country_loso/within_country_loso_summary_by_module_trait.csv\n")
cat("- tables/downsampling/balanced_downsampling_summary_by_module.csv\n")
cat("- tables/downsampling/balanced_downsampling_summary_by_module_trait.csv\n")
cat("- tables/model_stability/diagnosis/integrated_diagnosis_HC3_stability_summary.csv\n")
cat("- tables/model_stability/biomarker/integrated_focal_log_biomarker_HC3_stability_summary.csv\n")
cat("- tables/core_modules/core_modules_correlation_stability_summary.csv\n")
cat("- WGCNA_14_Association_Stability_LOCO_LOSO_Downsampling.xlsx\n")

###############################################################################
# END
###############################################################################


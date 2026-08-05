###############################################################################
# ReDLat plasma proteomics — WGCNA workflow
# 05. Evaluate site and biomarker sensitivity
# Requires: outputs from Script 04
# Produces: HC3 models, nested-site analyses and sensitivity tables
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
  "broom",
  "sandwich",
  "lmtest",
  "openxlsx",
  "scales",
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

ORIGINAL_ADJUSTED_BIOMARKER_FILE <- file.path(
  SCRIPT13_DIR,
  "tables",
  "regression",
  "adjusted_module_models.csv"
)

ORIGINAL_DIAGNOSIS_FILE <- file.path(
  SCRIPT13_DIR,
  "tables",
  "regression",
  "adjusted_diagnosis_module_models.csv"
)

OLD_SITE_MODEL_FILE <- file.path(
  SCRIPT13_DIR,
  "tables",
  "context",
  "site_categorical_adjusted_lm_by_module.csv"
)

COUNTRY_MODEL_FILE <- file.path(
  SCRIPT13_DIR,
  "tables",
  "context",
  "country_categorical_adjusted_lm_by_module.csv"
)

OUTDIR <- file.path(WGCNA_CONFIG$result_root,
  "05_sensitivity"
)

SUBDIRS <- c(
  "tables",
  "tables/context",
  "tables/biomarker_robustness",
  "tables/diagnosis_robustness",
  "figures",
  "figures/context",
  "figures/biomarker_robustness",
  "figures/diagnosis_robustness",
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
  ORIGINAL_ADJUSTED_BIOMARKER_FILE,
  ORIGINAL_DIAGNOSIS_FILE,
  OLD_SITE_MODEL_FILE,
  COUNTRY_MODEL_FILE
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  stop(
    "Missing required Script 13 files:\n",
    paste(missing_files, collapse = "\n"),
    "\n\nRun Script 13 first.",
    call. = FALSE
  )
}

###############################################################################
# 3) PARAMETERS
###############################################################################

EXPECTED_N_SAMPLES <- 639L
EXPECTED_N_MODULES <- 8L
REQUIRE_EXPECTED_DIMENSIONS <- TRUE

MODULES_OF_INTEREST <- NULL

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

CONTEXT_COVARS <- c(
  "SampleGroup_bin",
  "Age",
  "Sex_bin",
  "Education"
)

BIOMARKER_MODEL_COVARS <- c(
  "Age",
  "Sex_bin",
  "Education",
  "Country"
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

MIN_N_FOR_MODEL <- 30L
MIN_N_PER_SITE <- 10L

ROBUST_VCOV_TYPE <- "HC3"
ALPHA <- 0.05

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

safe_mean <- function(x) {
  x <- safe_numeric(x)
  x <- x[is.finite(x)]

  if (length(x) == 0) {
    return(NA_real_)
  }

  mean(x)
}

safe_sd <- function(x) {
  x <- safe_numeric(x)
  x <- x[is.finite(x)]

  if (length(x) < 2) {
    return(NA_real_)
  }

  stats::sd(x)
}

safe_skewness <- function(x) {
  x <- safe_numeric(x)
  x <- x[is.finite(x)]

  if (length(x) < 3) {
    return(NA_real_)
  }

  sx <- stats::sd(x)

  if (!is.finite(sx) || sx == 0) {
    return(NA_real_)
  }

  mean(
    (
      (x - mean(x)) / sx
    )^3
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

add_sig_stars <- function(q) {
  dplyr::case_when(
    is.na(q) ~ "",
    q < 0.001 ~ "***",
    q < 0.01 ~ "**",
    q < 0.05 ~ "*",
    TRUE ~ ""
  )
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

formula_text <- function(outcome, predictors) {
  paste0(
    outcome,
    " ~ ",
    paste(
      predictors,
      collapse = " + "
    )
  )
}

model_rank <- function(fit) {
  if (is.null(fit)) {
    return(NA_integer_)
  }

  fit$rank
}

model_rss <- function(fit) {
  if (is.null(fit)) {
    return(NA_real_)
  }

  sum(
    stats::residuals(fit)^2
  )
}

nested_model_comparison <- function(
    fit_reduced,
    fit_full
) {
  if (
    is.null(fit_reduced) ||
    is.null(fit_full)
  ) {
    return(tibble::tibble(
      rank_reduced = NA_integer_,
      rank_full = NA_integer_,
      added_df = NA_integer_,
      residual_df_full = NA_integer_,
      rss_reduced = NA_real_,
      rss_full = NA_real_,
      F_statistic = NA_real_,
      p_value = NA_real_,
      partial_R2 = NA_real_,
      model_status = "model_failed"
    ))
  }

  rank_reduced <- model_rank(
    fit_reduced
  )

  rank_full <- model_rank(
    fit_full
  )

  added_df <- rank_full -
    rank_reduced

  residual_df_full <- stats::df.residual(
    fit_full
  )

  rss_reduced <- model_rss(
    fit_reduced
  )

  rss_full <- model_rss(
    fit_full
  )

  if (
    !is.finite(rss_reduced) ||
    !is.finite(rss_full) ||
    added_df <= 0 ||
    residual_df_full <= 0 ||
    rss_reduced <= 0 ||
    rss_full <= 0
  ) {
    return(tibble::tibble(
      rank_reduced = rank_reduced,
      rank_full = rank_full,
      added_df = added_df,
      residual_df_full = residual_df_full,
      rss_reduced = rss_reduced,
      rss_full = rss_full,
      F_statistic = NA_real_,
      p_value = NA_real_,
      partial_R2 = NA_real_,
      model_status = "nonpositive_or_non_nested_df"
    ))
  }

  F_statistic <- (
    (rss_reduced - rss_full) /
      added_df
  ) / (
    rss_full /
      residual_df_full
  )

  if (
    !is.finite(F_statistic) ||
    F_statistic < 0
  ) {
    p_value <- NA_real_
  } else {
    p_value <- stats::pf(
      F_statistic,
      df1 = added_df,
      df2 = residual_df_full,
      lower.tail = FALSE
    )
  }

  partial_R2 <- (
    rss_reduced - rss_full
  ) / rss_reduced

  tibble::tibble(
    rank_reduced = rank_reduced,
    rank_full = rank_full,
    added_df = added_df,
    residual_df_full = residual_df_full,
    rss_reduced = rss_reduced,
    rss_full = rss_full,
    F_statistic = F_statistic,
    p_value = p_value,
    partial_R2 = partial_R2,
    model_status = "ok"
  )
}

fit_site_nested_all_sample <- function(
    df,
    module,
    site_var,
    covars = CONTEXT_COVARS
) {
  required_vars <- unique(
    c(
      module,
      "Country",
      site_var,
      covars
    )
  )

  available_vars <- required_vars[
    required_vars %in%
      names(df)
  ]

  if (!all(
    c(
      module,
      "Country",
      site_var
    ) %in%
      available_vars
  )) {
    return(tibble::tibble(
      Module = module,
      analysis_scope = "all_samples_nested_site_within_country",
      N = 0L,
      n_countries = NA_integer_,
      n_country_site_levels = NA_integer_,
      n_multisite_countries = NA_integer_,
      expected_added_site_df = NA_integer_,
      rank_reduced = NA_integer_,
      rank_full = NA_integer_,
      added_df = NA_integer_,
      residual_df_full = NA_integer_,
      F_statistic = NA_real_,
      p_value = NA_real_,
      partial_R2_site_within_country = NA_real_,
      full_model_R2 = NA_real_,
      full_model_adj_R2 = NA_real_,
      reduced_formula = NA_character_,
      full_formula = NA_character_,
      model_status = "required_variables_missing"
    ))
  }

  covars_use <- covars[
    covars %in%
      names(df)
  ]

  d <- df %>%
    dplyr::select(
      dplyr::all_of(
        unique(
          c(
            module,
            "Country",
            site_var,
            covars_use
          )
        )
      )
    ) %>%
    tidyr::drop_na()

  if (nrow(d) < MIN_N_FOR_MODEL) {
    return(tibble::tibble(
      Module = module,
      analysis_scope = "all_samples_nested_site_within_country",
      N = nrow(d),
      n_countries = NA_integer_,
      n_country_site_levels = NA_integer_,
      n_multisite_countries = NA_integer_,
      expected_added_site_df = NA_integer_,
      rank_reduced = NA_integer_,
      rank_full = NA_integer_,
      added_df = NA_integer_,
      residual_df_full = NA_integer_,
      F_statistic = NA_real_,
      p_value = NA_real_,
      partial_R2_site_within_country = NA_real_,
      full_model_R2 = NA_real_,
      full_model_adj_R2 = NA_real_,
      reduced_formula = NA_character_,
      full_formula = NA_character_,
      model_status = "insufficient_n"
    ))
  }

  d$Country <- droplevels(
    factor(d$Country)
  )

  d[[site_var]] <- droplevels(
    factor(d[[site_var]])
  )

  d$Country_Site <- droplevels(
    interaction(
      d$Country,
      d[[site_var]],
      sep = "::",
      drop = TRUE
    )
  )

  sites_per_country <- d %>%
    dplyr::distinct(
      Country,
      .data[[site_var]]
    ) %>%
    dplyr::count(
      Country,
      name = "n_sites"
    )

  expected_added_df <- sum(
    pmax(
      sites_per_country$n_sites - 1L,
      0L
    )
  )

  predictors_reduced <- c(
    covars_use,
    "Country"
  )

  predictors_full <- c(
    covars_use,
    "Country",
    "Country_Site"
  )

  reduced_formula <- formula_text(
    module,
    predictors_reduced
  )

  full_formula <- formula_text(
    module,
    predictors_full
  )

  fit_reduced <- tryCatch(
    stats::lm(
      stats::as.formula(
        reduced_formula
      ),
      data = d
    ),
    error = function(e) NULL
  )

  fit_full <- tryCatch(
    stats::lm(
      stats::as.formula(
        full_formula
      ),
      data = d
    ),
    error = function(e) NULL
  )

  cmp <- nested_model_comparison(
    fit_reduced,
    fit_full
  )

  fit_summary <- if (
    is.null(fit_full)
  ) {
    NULL
  } else {
    summary(fit_full)
  }

  status <- cmp$model_status[[1]]

  if (
    status == "ok" &&
    cmp$added_df[[1]] !=
      expected_added_df
  ) {
    status <- paste0(
      "ok_but_added_df_",
      cmp$added_df[[1]],
      "_differs_from_expected_",
      expected_added_df
    )
  }

  tibble::tibble(
    Module = module,
    analysis_scope = "all_samples_nested_site_within_country",
    N = nrow(d),
    n_countries = nlevels(
      d$Country
    ),
    n_country_site_levels = nlevels(
      d$Country_Site
    ),
    n_multisite_countries = sum(
      sites_per_country$n_sites >= 2
    ),
    expected_added_site_df = expected_added_df,
    rank_reduced = cmp$rank_reduced,
    rank_full = cmp$rank_full,
    added_df = cmp$added_df,
    residual_df_full = cmp$residual_df_full,
    F_statistic = cmp$F_statistic,
    p_value = cmp$p_value,
    partial_R2_site_within_country =
      cmp$partial_R2,
    full_model_R2 = if (
      is.null(fit_summary)
    ) {
      NA_real_
    } else {
      fit_summary$r.squared
    },
    full_model_adj_R2 = if (
      is.null(fit_summary)
    ) {
      NA_real_
    } else {
      fit_summary$adj.r.squared
    },
    reduced_formula = reduced_formula,
    full_formula = full_formula,
    model_status = status
  )
}

fit_site_within_country <- function(
    df,
    module,
    country_value,
    site_var,
    covars = CONTEXT_COVARS
) {
  covars_use <- covars[
    covars %in%
      names(df)
  ]

  required_vars <- unique(
    c(
      module,
      "Country",
      site_var,
      covars_use
    )
  )

  d <- df %>%
    dplyr::filter(
      as.character(Country) ==
        country_value
    ) %>%
    dplyr::select(
      dplyr::all_of(
        required_vars
      )
    ) %>%
    tidyr::drop_na()

  if (
    nrow(d) < MIN_N_FOR_MODEL
  ) {
    return(tibble::tibble(
      Country = country_value,
      Module = module,
      N = nrow(d),
      n_sites = NA_integer_,
      minimum_site_n = NA_integer_,
      rank_reduced = NA_integer_,
      rank_full = NA_integer_,
      added_df = NA_integer_,
      residual_df_full = NA_integer_,
      F_statistic = NA_real_,
      p_value = NA_real_,
      partial_R2_site = NA_real_,
      full_model_R2 = NA_real_,
      full_model_adj_R2 = NA_real_,
      reduced_formula = NA_character_,
      full_formula = NA_character_,
      model_status = "insufficient_n"
    ))
  }

  d[[site_var]] <- droplevels(
    factor(d[[site_var]])
  )

  site_counts <- table(
    d[[site_var]]
  )

  n_sites <- nlevels(
    d[[site_var]]
  )

  minimum_site_n <- min(
    as.integer(site_counts)
  )

  if (n_sites < 2) {
    return(tibble::tibble(
      Country = country_value,
      Module = module,
      N = nrow(d),
      n_sites = n_sites,
      minimum_site_n = minimum_site_n,
      rank_reduced = NA_integer_,
      rank_full = NA_integer_,
      added_df = NA_integer_,
      residual_df_full = NA_integer_,
      F_statistic = NA_real_,
      p_value = NA_real_,
      partial_R2_site = NA_real_,
      full_model_R2 = NA_real_,
      full_model_adj_R2 = NA_real_,
      reduced_formula = NA_character_,
      full_formula = NA_character_,
      model_status = "fewer_than_two_sites"
    ))
  }

  predictors_reduced <- covars_use
  predictors_full <- c(
    site_var,
    covars_use
  )

  reduced_formula <- formula_text(
    module,
    predictors_reduced
  )

  full_formula <- formula_text(
    module,
    predictors_full
  )

  fit_reduced <- tryCatch(
    stats::lm(
      stats::as.formula(
        reduced_formula
      ),
      data = d
    ),
    error = function(e) NULL
  )

  fit_full <- tryCatch(
    stats::lm(
      stats::as.formula(
        full_formula
      ),
      data = d
    ),
    error = function(e) NULL
  )

  cmp <- nested_model_comparison(
    fit_reduced,
    fit_full
  )

  fit_summary <- if (
    is.null(fit_full)
  ) {
    NULL
  } else {
    summary(fit_full)
  }

  status <- cmp$model_status[[1]]

  if (
    minimum_site_n <
      MIN_N_PER_SITE
  ) {
    status <- paste0(
      status,
      "_small_site_cell"
    )
  }

  tibble::tibble(
    Country = country_value,
    Module = module,
    N = nrow(d),
    n_sites = n_sites,
    minimum_site_n = minimum_site_n,
    rank_reduced = cmp$rank_reduced,
    rank_full = cmp$rank_full,
    added_df = cmp$added_df,
    residual_df_full = cmp$residual_df_full,
    F_statistic = cmp$F_statistic,
    p_value = cmp$p_value,
    partial_R2_site = cmp$partial_R2,
    full_model_R2 = if (
      is.null(fit_summary)
    ) {
      NA_real_
    } else {
      fit_summary$r.squared
    },
    full_model_adj_R2 = if (
      is.null(fit_summary)
    ) {
      NA_real_
    } else {
      fit_summary$adj.r.squared
    },
    reduced_formula = reduced_formula,
    full_formula = full_formula,
    model_status = status
  )
}

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
      residual_df = NA_real_,
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
      residual_df = stats::df.residual(
        fit
      ),
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
      residual_df = stats::df.residual(
        fit
      ),
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
    1 - ALPHA / 2,
    df = residual_df
  )

  tibble::tibble(
    estimate = estimate,
    robust_se = robust_se,
    robust_statistic = robust_statistic,
    robust_p_value = robust_p_value,
    robust_conf_low = estimate -
      critical_value *
      robust_se,
    robust_conf_high = estimate +
      critical_value *
      robust_se,
    residual_df = residual_df,
    model_status = "ok"
  )
}

fit_biomarker_hc3 <- function(
    df,
    module,
    biomarker,
    transform = c(
      "raw",
      "log"
    ),
    covars = BIOMARKER_MODEL_COVARS
) {
  transform <- match.arg(
    transform
  )

  covars_use <- covars[
    covars %in%
      names(df)
  ]

  required_vars <- unique(
    c(
      biomarker,
      module,
      covars_use
    )
  )

  if (!all(
    c(
      biomarker,
      module
    ) %in%
      names(df)
  )) {
    return(tibble::tibble(
      Module = module,
      Biomarker = biomarker,
      Transform = transform,
      N = 0L,
      n_nonpositive_excluded = NA_integer_,
      estimate = NA_real_,
      robust_se = NA_real_,
      robust_statistic = NA_real_,
      robust_p_value = NA_real_,
      robust_conf_low = NA_real_,
      robust_conf_high = NA_real_,
      module_SD = NA_real_,
      outcome_SD = NA_real_,
      effect_per_1SD_module = NA_real_,
      effect_per_1SD_module_conf_low = NA_real_,
      effect_per_1SD_module_conf_high = NA_real_,
      standardized_beta = NA_real_,
      percent_change_per_1SD_module = NA_real_,
      percent_change_conf_low = NA_real_,
      percent_change_conf_high = NA_real_,
      raw_change_per_1SD_module = NA_real_,
      raw_change_conf_low = NA_real_,
      raw_change_conf_high = NA_real_,
      Breusch_Pagan_p = NA_real_,
      model_R2 = NA_real_,
      model_adj_R2 = NA_real_,
      Formula = NA_character_,
      model_status = "required_variables_missing"
    ))
  }

  d <- df %>%
    dplyr::select(
      dplyr::all_of(
        required_vars
      )
    )

  d[[biomarker]] <- safe_numeric(
    d[[biomarker]]
  )

  d[[module]] <- safe_numeric(
    d[[module]]
  )

  n_nonpositive_excluded <- sum(
    !is.na(d[[biomarker]]) &
      d[[biomarker]] <= 0
  )

  if (transform == "log") {
    d$.outcome_model <- ifelse(
      is.finite(d[[biomarker]]) &
        d[[biomarker]] > 0,
      log(
        d[[biomarker]]
      ),
      NA_real_
    )
  } else {
    d$.outcome_model <- d[[biomarker]]
  }

  d <- d %>%
    tidyr::drop_na(
      dplyr::all_of(
        c(
          ".outcome_model",
          module,
          covars_use
        )
      )
    )

  if (
    nrow(d) < MIN_N_FOR_MODEL ||
    dplyr::n_distinct(
      d$.outcome_model
    ) <= 1 ||
    dplyr::n_distinct(
      d[[module]]
    ) <= 1
  ) {
    return(tibble::tibble(
      Module = module,
      Biomarker = biomarker,
      Transform = transform,
      N = nrow(d),
      n_nonpositive_excluded =
        n_nonpositive_excluded,
      estimate = NA_real_,
      robust_se = NA_real_,
      robust_statistic = NA_real_,
      robust_p_value = NA_real_,
      robust_conf_low = NA_real_,
      robust_conf_high = NA_real_,
      module_SD = safe_sd(
        d[[module]]
      ),
      outcome_SD = safe_sd(
        d$.outcome_model
      ),
      effect_per_1SD_module = NA_real_,
      effect_per_1SD_module_conf_low = NA_real_,
      effect_per_1SD_module_conf_high = NA_real_,
      standardized_beta = NA_real_,
      percent_change_per_1SD_module = NA_real_,
      percent_change_conf_low = NA_real_,
      percent_change_conf_high = NA_real_,
      raw_change_per_1SD_module = NA_real_,
      raw_change_conf_low = NA_real_,
      raw_change_conf_high = NA_real_,
      Breusch_Pagan_p = NA_real_,
      model_R2 = NA_real_,
      model_adj_R2 = NA_real_,
      Formula = NA_character_,
      model_status = "insufficient_or_constant"
    ))
  }

  if ("Country" %in% names(d)) {
    d$Country <- droplevels(
      factor(d$Country)
    )
  }

  form_txt <- formula_text(
    ".outcome_model",
    c(
      module,
      covars_use
    )
  )

  fit <- tryCatch(
    stats::lm(
      stats::as.formula(
        form_txt
      ),
      data = d
    ),
    error = function(e) NULL
  )

  robust_term <- extract_hc3_term(
    fit,
    term = module,
    vcov_type = ROBUST_VCOV_TYPE
  )

  module_sd <- safe_sd(
    d[[module]]
  )

  outcome_sd <- safe_sd(
    d$.outcome_model
  )

  effect_per_sd <- robust_term$estimate *
    module_sd

  effect_low_per_sd <- robust_term$robust_conf_low *
    module_sd

  effect_high_per_sd <- robust_term$robust_conf_high *
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

  if (transform == "log") {
    percent_change_per_sd <- 100 *
      (
        exp(
          effect_per_sd
        ) - 1
      )

    percent_change_low <- 100 *
      (
        exp(
          effect_low_per_sd
        ) - 1
      )

    percent_change_high <- 100 *
      (
        exp(
          effect_high_per_sd
        ) - 1
      )

    raw_change_per_sd <- NA_real_
    raw_change_low <- NA_real_
    raw_change_high <- NA_real_
  } else {
    percent_change_per_sd <- NA_real_
    percent_change_low <- NA_real_
    percent_change_high <- NA_real_

    raw_change_per_sd <- effect_per_sd
    raw_change_low <- effect_low_per_sd
    raw_change_high <- effect_high_per_sd
  }

  bp_p <- if (is.null(fit)) {
    NA_real_
  } else {
    tryCatch(
      lmtest::bptest(
        fit
      )$p.value,
      error = function(e) NA_real_
    )
  }

  fit_summary <- if (is.null(fit)) {
    NULL
  } else {
    summary(fit)
  }

  tibble::tibble(
    Module = module,
    Biomarker = biomarker,
    Transform = transform,
    N = nrow(d),
    n_nonpositive_excluded =
      n_nonpositive_excluded,
    estimate = robust_term$estimate,
    robust_se = robust_term$robust_se,
    robust_statistic =
      robust_term$robust_statistic,
    robust_p_value =
      robust_term$robust_p_value,
    robust_conf_low =
      robust_term$robust_conf_low,
    robust_conf_high =
      robust_term$robust_conf_high,
    module_SD = module_sd,
    outcome_SD = outcome_sd,
    effect_per_1SD_module =
      effect_per_sd,
    effect_per_1SD_module_conf_low =
      effect_low_per_sd,
    effect_per_1SD_module_conf_high =
      effect_high_per_sd,
    standardized_beta =
      standardized_beta,
    percent_change_per_1SD_module =
      percent_change_per_sd,
    percent_change_conf_low =
      percent_change_low,
    percent_change_conf_high =
      percent_change_high,
    raw_change_per_1SD_module =
      raw_change_per_sd,
    raw_change_conf_low =
      raw_change_low,
    raw_change_conf_high =
      raw_change_high,
    Breusch_Pagan_p = bp_p,
    model_R2 = if (
      is.null(fit_summary)
    ) {
      NA_real_
    } else {
      fit_summary$r.squared
    },
    model_adj_R2 = if (
      is.null(fit_summary)
    ) {
      NA_real_
    } else {
      fit_summary$adj.r.squared
    },
    Formula = form_txt,
    model_status =
      robust_term$model_status
  )
}

fit_diagnosis_hc3 <- function(
    df,
    module
) {
  covars <- c(
    "Age",
    "Sex_bin",
    "Education",
    "Country"
  )

  covars_use <- covars[
    covars %in%
      names(df)
  ]

  required_vars <- unique(
    c(
      module,
      "SampleGroup_bin",
      covars_use
    )
  )

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
    ) < 2
  ) {
    return(tibble::tibble(
      Module = module,
      N = nrow(d),
      estimate_diagnosis = NA_real_,
      robust_se = NA_real_,
      robust_statistic = NA_real_,
      robust_p_value = NA_real_,
      robust_conf_low = NA_real_,
      robust_conf_high = NA_real_,
      module_SD = safe_sd(
        d[[module]]
      ),
      standardized_difference =
        NA_real_,
      model_R2 = NA_real_,
      model_adj_R2 = NA_real_,
      Formula = NA_character_,
      model_status =
        "insufficient_or_constant"
    ))
  }

  d$Country <- droplevels(
    factor(d$Country)
  )

  form_txt <- formula_text(
    module,
    c(
      "SampleGroup_bin",
      covars_use
    )
  )

  fit <- tryCatch(
    stats::lm(
      stats::as.formula(
        form_txt
      ),
      data = d
    ),
    error = function(e) NULL
  )

  robust_term <- extract_hc3_term(
    fit,
    term = "SampleGroup_bin",
    vcov_type = ROBUST_VCOV_TYPE
  )

  module_sd <- safe_sd(
    d[[module]]
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

  fit_summary <- if (is.null(fit)) {
    NULL
  } else {
    summary(fit)
  }

  tibble::tibble(
    Module = module,
    N = nrow(d),
    estimate_diagnosis =
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
    standardized_difference =
      standardized_difference,
    model_R2 = if (
      is.null(fit_summary)
    ) {
      NA_real_
    } else {
      fit_summary$r.squared
    },
    model_adj_R2 = if (
      is.null(fit_summary)
    ) {
      NA_real_
    } else {
      fit_summary$adj.r.squared
    },
    Formula = form_txt,
    model_status =
      robust_term$model_status
  )
}

classify_diagnosis_evidence <- function(
    p,
    fdr
) {
  dplyr::case_when(
    is.na(p) |
      is.na(fdr) ~
      "Not estimable",
    fdr < 0.05 ~
      "FDR-significant",
    p < 0.05 &
      fdr < 0.10 ~
      "Nominal and FDR-borderline",
    p < 0.05 ~
      "Nominal only",
    TRUE ~
      "Not nominally significant"
  )
}

###############################################################################
# 5) LOAD AND PREPARE DATA
###############################################################################

analysis_df <- safe_read_csv(
  INPUT_FILE
)

module_reference <- safe_read_csv(
  MODULE_REFERENCE_FILE
)

original_adjusted <- safe_read_csv(
  ORIGINAL_ADJUSTED_BIOMARKER_FILE
)

original_diagnosis <- safe_read_csv(
  ORIGINAL_DIAGNOSIS_FILE
)

old_site_models <- safe_read_csv(
  OLD_SITE_MODEL_FILE
)

country_models <- safe_read_csv(
  COUNTRY_MODEL_FILE
)

required_input_cols <- c(
  "SampleId",
  "Country",
  "SampleGroup_bin",
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
    Country = droplevels(
      factor(
        clean_text_na(
          Country
        )
      )
    ),
    SampleGroup_bin = safe_numeric(
      SampleGroup_bin
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

if (length(modules_use) == 0) {
  stop(
    "No module eigengene columns were found.",
    call. = FALSE
  )
}

for (mod in modules_use) {
  analysis_df[[mod]] <- safe_numeric(
    analysis_df[[mod]]
  )
}

for (bio in intersect(
  BIOMARKERS,
  names(analysis_df)
)) {
  analysis_df[[bio]] <- safe_numeric(
    analysis_df[[bio]]
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
}

missing_biomarkers <- setdiff(
  BIOMARKERS,
  names(analysis_df)
)

if (length(missing_biomarkers) > 0) {
  stop(
    "Missing required biomarker columns: ",
    paste(
      missing_biomarkers,
      collapse = ", "
    ),
    call. = FALSE
  )
}

###############################################################################
# 6) COUNTRY-SITE STRUCTURE AUDIT
###############################################################################

country_site_counts <- analysis_df %>%
  dplyr::filter(
    !is.na(Country),
    !is.na(.data[[site_var]])
  ) %>%
  dplyr::count(
    Country,
    Site = .data[[site_var]],
    name = "N"
  ) %>%
  dplyr::arrange(
    Country,
    Site
  )

sites_per_country <- country_site_counts %>%
  dplyr::count(
    Country,
    name = "n_sites"
  ) %>%
  dplyr::left_join(
    country_site_counts %>%
      dplyr::group_by(Country) %>%
      dplyr::summarise(
        total_N = sum(N),
        minimum_site_N = min(N),
        maximum_site_N = max(N),
        sites = paste(
          Site,
          collapse = ", "
        ),
        .groups = "drop"
      ),
    by = "Country"
  ) %>%
  dplyr::mutate(
    contributes_within_country_site_df =
      pmax(
        n_sites - 1L,
        0L
      ),
    contributes_to_nested_site_test =
      n_sites >= 2
  )

multisite_countries <- sites_per_country %>%
  dplyr::filter(
    n_sites >= 2
  ) %>%
  dplyr::pull(Country) %>%
  as.character()

expected_nested_site_df <- sum(
  sites_per_country$
    contributes_within_country_site_df
)

if (length(multisite_countries) == 0) {
  stop(
    "No country contains more than one recruitment site; ",
    "site-within-country effects are not estimable.",
    call. = FALSE
  )
}

safe_write_csv(
  country_site_counts,
  file.path(
    OUTDIR,
    "tables",
    "context",
    "country_site_sample_counts.csv"
  )
)

safe_write_csv(
  sites_per_country,
  file.path(
    OUTDIR,
    "tables",
    "context",
    "site_nested_structure_audit.csv"
  )
)

context_structure_summary <- tibble::tibble(
  metric = c(
    "site_variable",
    "n_samples",
    "n_countries",
    "n_site_levels",
    "n_observed_country_site_combinations",
    "n_multisite_countries",
    "multisite_countries",
    "expected_site_within_country_df",
    "single_site_countries",
    "interpretation"
  ),
  value = c(
    site_var,
    as.character(nrow(analysis_df)),
    as.character(
      dplyr::n_distinct(
        analysis_df$Country,
        na.rm = TRUE
      )
    ),
    as.character(
      dplyr::n_distinct(
        analysis_df[[site_var]],
        na.rm = TRUE
      )
    ),
    as.character(
      nrow(country_site_counts)
    ),
    as.character(
      length(multisite_countries)
    ),
    paste(
      multisite_countries,
      collapse = ", "
    ),
    as.character(
      expected_nested_site_df
    ),
    paste(
      sites_per_country$Country[
        sites_per_country$n_sites == 1
      ],
      collapse = ", "
    ),
    paste(
      "Only within-country site contrasts are interpretable.",
      "Countries containing one site contribute to country adjustment",
      "but not to additional site-within-country degrees of freedom."
    )
  )
)

safe_write_csv(
  context_structure_summary,
  file.path(
    OUTDIR,
    "tables",
    "context",
    "site_nested_analysis_structure_summary.csv"
  )
)

###############################################################################
# 7) CORRECTED ALL-SAMPLE NESTED SITE MODELS
###############################################################################

nested_site_all <- purrr::map_dfr(
  modules_use,
  function(mod) {
    fit_site_nested_all_sample(
      df = analysis_df,
      module = mod,
      site_var = site_var,
      covars = CONTEXT_COVARS
    )
  }
) %>%
  dplyr::mutate(
    FDR_site_within_country = p.adjust(
      p_value,
      method = "BH"
    )
  ) %>%
  dplyr::left_join(
    module_reference,
    by = "Module"
  ) %>%
  dplyr::arrange(
    FDR_site_within_country,
    p_value
  )

if (
  any(
    nested_site_all$added_df !=
      expected_nested_site_df,
    na.rm = TRUE
  )
) {
  warning(
    "At least one module has an added nested-site model rank ",
    "different from the expected site-within-country degrees of freedom. ",
    "Review tables/context/corrected_nested_site_all_samples.csv."
  )
}

safe_write_csv(
  nested_site_all,
  file.path(
    OUTDIR,
    "tables",
    "context",
    "corrected_nested_site_all_samples.csv"
  )
)

###############################################################################
# 8) WITHIN-MULTISITE-COUNTRY SENSITIVITY
###############################################################################

site_within_country <- purrr::map_dfr(
  multisite_countries,
  function(country_value) {
    purrr::map_dfr(
      modules_use,
      function(mod) {
        fit_site_within_country(
          df = analysis_df,
          module = mod,
          country_value =
            country_value,
          site_var = site_var,
          covars = CONTEXT_COVARS
        )
      }
    )
  }
) %>%
  dplyr::group_by(Country) %>%
  dplyr::mutate(
    FDR_within_country = p.adjust(
      p_value,
      method = "BH"
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(
    module_reference,
    by = "Module"
  ) %>%
  dplyr::arrange(
    Country,
    FDR_within_country,
    p_value
  )

safe_write_csv(
  site_within_country,
  file.path(
    OUTDIR,
    "tables",
    "context",
    "corrected_site_models_within_multisite_countries.csv"
  )
)

site_within_country_summary <- site_within_country %>%
  dplyr::group_by(Module) %>%
  dplyr::summarise(
    n_multisite_countries_tested =
      dplyr::n(),
    n_nominal_p_lt_0_05 = sum(
      p_value < 0.05,
      na.rm = TRUE
    ),
    n_FDR_lt_0_05 = sum(
      FDR_within_country < 0.05,
      na.rm = TRUE
    ),
    mean_partial_R2 =
      safe_mean(
        partial_R2_site
      ),
    max_partial_R2 =
      safe_max(
        partial_R2_site
      ),
    same_direction_not_applicable =
      "Site is a multi-level factor; direction is not summarized as a single sign.",
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    module_reference,
    by = "Module"
  )

safe_write_csv(
  site_within_country_summary,
  file.path(
    OUTDIR,
    "tables",
    "context",
    "corrected_site_within_country_summary_by_module.csv"
  )
)

###############################################################################
# 9) COMPARISON WITH THE OLD NON-NESTED SITE MODEL
###############################################################################

old_site_comparison <- old_site_models %>%
  dplyr::transmute(
    Module = as.character(Module),
    old_non_nested_site_p =
      safe_numeric(
        p_value_factor
      ),
    old_non_nested_site_FDR =
      safe_numeric(
        FDR_factor
      ),
    old_non_nested_site_partial_R2 =
      safe_numeric(
        partial_R2_factor
      ),
    old_non_nested_model_status =
      as.character(
        model_status
      )
  ) %>%
  dplyr::right_join(
    nested_site_all %>%
      dplyr::select(
        Module,
        corrected_nested_site_p =
          p_value,
        corrected_nested_site_FDR =
          FDR_site_within_country,
        corrected_nested_site_partial_R2 =
          partial_R2_site_within_country,
        corrected_added_df =
          added_df,
        corrected_model_status =
          model_status
      ),
    by = "Module"
  ) %>%
  dplyr::mutate(
    partial_R2_attenuation =
      old_non_nested_site_partial_R2 -
      corrected_nested_site_partial_R2,
    corrected_to_old_partial_R2_ratio =
      corrected_nested_site_partial_R2 /
      old_non_nested_site_partial_R2,
    interpretation = paste(
      "The corrected estimate isolates site-within-country variation.",
      "The old estimate combined country and site variation and is retained",
      "only for comparison."
    )
  ) %>%
  dplyr::left_join(
    module_reference,
    by = "Module"
  ) %>%
  dplyr::arrange(
    corrected_nested_site_FDR,
    corrected_nested_site_p
  )

safe_write_csv(
  old_site_comparison,
  file.path(
    OUTDIR,
    "tables",
    "context",
    "old_non_nested_vs_corrected_nested_site_comparison.csv"
  )
)

###############################################################################
# 10) BIOMARKER DISTRIBUTION AND LOG-TRANSFORMATION AUDIT
###############################################################################

biomarker_transform_audit <- purrr::map_dfr(
  BIOMARKERS,
  function(bio) {
    x <- safe_numeric(
      analysis_df[[bio]]
    )

    x_log <- ifelse(
      is.finite(x) &
        x > 0,
      log(x),
      NA_real_
    )

    tibble::tibble(
      Biomarker = bio,
      Biomarker_label = dplyr::recode(
        bio,
        !!!BIOMARKER_LABELS,
        .default = bio
      ),
      N_nonmissing_raw = sum(
        is.finite(x)
      ),
      N_positive = sum(
        is.finite(x) &
          x > 0
      ),
      N_zero = sum(
        is.finite(x) &
          x == 0
      ),
      N_negative = sum(
        is.finite(x) &
          x < 0
      ),
      minimum_raw = safe_min(x),
      median_raw = stats::median(
        x,
        na.rm = TRUE
      ),
      maximum_raw = safe_max(x),
      skewness_raw =
        safe_skewness(x),
      minimum_log =
        safe_min(x_log),
      median_log = stats::median(
        x_log,
        na.rm = TRUE
      ),
      maximum_log =
        safe_max(x_log),
      skewness_log =
        safe_skewness(x_log),
      absolute_skewness_reduction =
        abs(
          safe_skewness(x)
        ) -
        abs(
          safe_skewness(x_log)
        ),
      log_transform_valid =
        sum(
          is.finite(x) &
            x <= 0
        ) == 0
    )
  }
)

safe_write_csv(
  biomarker_transform_audit,
  file.path(
    OUTDIR,
    "tables",
    "biomarker_robustness",
    "biomarker_log_transformation_audit.csv"
  )
)

if (
  any(
    biomarker_transform_audit$
      N_zero > 0 |
      biomarker_transform_audit$
        N_negative > 0
  )
) {
  warning(
    "At least one biomarker contains nonpositive values. ",
    "Those values are excluded from natural-log sensitivity models."
  )
}

###############################################################################
# 11) RAW-HC3 AND LOG-HC3 BIOMARKER MODELS
###############################################################################

biomarker_robust_models <- purrr::map_dfr(
  modules_use,
  function(mod) {
    purrr::map_dfr(
      BIOMARKERS,
      function(bio) {
        dplyr::bind_rows(
          fit_biomarker_hc3(
            df = analysis_df,
            module = mod,
            biomarker = bio,
            transform = "raw",
            covars =
              BIOMARKER_MODEL_COVARS
          ),
          fit_biomarker_hc3(
            df = analysis_df,
            module = mod,
            biomarker = bio,
            transform = "log",
            covars =
              BIOMARKER_MODEL_COVARS
          )
        )
      }
    )
  }
) %>%
  dplyr::group_by(Transform) %>%
  dplyr::mutate(
    FDR_within_transform_family =
      p.adjust(
        robust_p_value,
        method = "BH"
      )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    FDR_across_all_robust_models =
      p.adjust(
        robust_p_value,
        method = "BH"
      ),
    Biomarker_label = dplyr::recode(
      Biomarker,
      !!!BIOMARKER_LABELS,
      .default = Biomarker
    ),
    significance_stars =
      add_sig_stars(
        FDR_within_transform_family
      ),
    heteroskedasticity_flag =
      !is.na(Breusch_Pagan_p) &
      Breusch_Pagan_p < 0.05
  ) %>%
  dplyr::left_join(
    module_reference,
    by = "Module"
  ) %>%
  dplyr::arrange(
    Transform,
    FDR_within_transform_family,
    robust_p_value
  )

expected_robust_models <- length(
  modules_use
) *
  length(
    BIOMARKERS
  ) *
  2L

if (nrow(biomarker_robust_models) !=
    expected_robust_models) {
  stop(
    "Expected ",
    expected_robust_models,
    " robust biomarker models, but obtained ",
    nrow(biomarker_robust_models),
    ".",
    call. = FALSE
  )
}

safe_write_csv(
  biomarker_robust_models,
  file.path(
    OUTDIR,
    "tables",
    "biomarker_robustness",
    "biomarker_models_raw_and_log_HC3.csv"
  )
)

safe_write_csv(
  biomarker_robust_models %>%
    dplyr::filter(
      Transform == "log"
    ),
  file.path(
    OUTDIR,
    "tables",
    "biomarker_robustness",
    "biomarker_models_log_HC3_primary_sensitivity.csv"
  )
)

safe_write_csv(
  biomarker_robust_models %>%
    dplyr::filter(
      Transform == "raw"
    ),
  file.path(
    OUTDIR,
    "tables",
    "biomarker_robustness",
    "biomarker_models_raw_HC3_sensitivity.csv"
  )
)

###############################################################################
# 12) ORIGINAL OLS VS ROBUST MODEL COMPARISON
###############################################################################

original_biomarker_ols <- original_adjusted %>%
  dplyr::filter(
    Outcome %in%
      BIOMARKERS
  ) %>%
  dplyr::transmute(
    Module = as.character(Module),
    Biomarker = as.character(Outcome),
    original_raw_OLS_N = safe_numeric(N),
    original_raw_OLS_estimate =
      safe_numeric(estimate),
    original_raw_OLS_p =
      safe_numeric(p.value),
    original_raw_OLS_FDR =
      safe_numeric(FDR)
  )

biomarker_model_comparison <- biomarker_robust_models %>%
  dplyr::left_join(
    original_biomarker_ols,
    by = c(
      "Module",
      "Biomarker"
    )
  ) %>%
  dplyr::mutate(
    sign_matches_original_raw_OLS =
      dplyr::case_when(
        is.na(estimate) |
          is.na(
            original_raw_OLS_estimate
          ) ~ NA,
        estimate == 0 |
          original_raw_OLS_estimate == 0 ~
          estimate ==
            original_raw_OLS_estimate,
        TRUE ~
          sign(estimate) ==
            sign(
              original_raw_OLS_estimate
            )
      ),
    robust_FDR_lt_0_05 =
      FDR_within_transform_family <
      0.05,
    original_OLS_FDR_lt_0_05 =
      original_raw_OLS_FDR <
      0.05,
    significance_preservation =
      dplyr::case_when(
        is.na(
          robust_FDR_lt_0_05
        ) |
          is.na(
            original_OLS_FDR_lt_0_05
          ) ~ "Not estimable",
        original_OLS_FDR_lt_0_05 &
          robust_FDR_lt_0_05 ~
          "Significant in both",
        original_OLS_FDR_lt_0_05 &
          !robust_FDR_lt_0_05 ~
          "Original OLS only",
        !original_OLS_FDR_lt_0_05 &
          robust_FDR_lt_0_05 ~
          "Robust sensitivity only",
        TRUE ~
          "Not significant in either"
      )
  ) %>%
  dplyr::arrange(
    Transform,
    FDR_within_transform_family
  )

safe_write_csv(
  biomarker_model_comparison,
  file.path(
    OUTDIR,
    "tables",
    "biomarker_robustness",
    "original_OLS_vs_HC3_model_comparison.csv"
  )
)

log_hc3_summary_by_module <- biomarker_robust_models %>%
  dplyr::filter(
    Transform == "log"
  ) %>%
  dplyr::group_by(Module) %>%
  dplyr::summarise(
    n_log_HC3_models =
      dplyr::n(),
    n_log_HC3_FDR_lt_0_05 =
      sum(
        FDR_within_transform_family <
          0.05,
        na.rm = TRUE
      ),
    minimum_log_HC3_FDR =
      safe_min(
        FDR_within_transform_family
      ),
    maximum_abs_standardized_beta =
      safe_max(
        abs(
          standardized_beta
        )
      ),
    n_heteroskedasticity_flags =
      sum(
        heteroskedasticity_flag,
        na.rm = TRUE
      ),
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    module_reference,
    by = "Module"
  )

safe_write_csv(
  log_hc3_summary_by_module,
  file.path(
    OUTDIR,
    "tables",
    "biomarker_robustness",
    "log_HC3_summary_by_module.csv"
  )
)

###############################################################################
# 13) ROBUST DIAGNOSIS MODELS
###############################################################################

diagnosis_hc3 <- purrr::map_dfr(
  modules_use,
  function(mod) {
    fit_diagnosis_hc3(
      df = analysis_df,
      module = mod
    )
  }
) %>%
  dplyr::mutate(
    FDR_HC3 = p.adjust(
      robust_p_value,
      method = "BH"
    ),
    evidence_class_HC3 =
      classify_diagnosis_evidence(
        robust_p_value,
        FDR_HC3
      )
  ) %>%
  dplyr::left_join(
    module_reference,
    by = "Module"
  ) %>%
  dplyr::arrange(
    FDR_HC3,
    robust_p_value
  )

original_diagnosis_clean <- original_diagnosis %>%
  dplyr::transmute(
    Module = as.character(Module),
    original_OLS_N = safe_numeric(N),
    original_OLS_estimate =
      safe_numeric(
        estimate_diagnosis
      ),
    original_OLS_p =
      safe_numeric(p.value),
    original_OLS_FDR =
      safe_numeric(FDR),
    original_OLS_evidence_class =
      classify_diagnosis_evidence(
        original_OLS_p,
        original_OLS_FDR
      )
  )

diagnosis_comparison <- diagnosis_hc3 %>%
  dplyr::left_join(
    original_diagnosis_clean,
    by = "Module"
  ) %>%
  dplyr::mutate(
    sign_matches_original_OLS =
      dplyr::case_when(
        is.na(
          estimate_diagnosis
        ) |
          is.na(
            original_OLS_estimate
          ) ~ NA,
        TRUE ~
          sign(
            estimate_diagnosis
          ) ==
            sign(
              original_OLS_estimate
            )
      ),
    interpretation = dplyr::case_when(
      evidence_class_HC3 ==
        "FDR-significant" ~
        "Adjusted diagnosis association remains significant with HC3 robust inference.",
      evidence_class_HC3 ==
        "Nominal and FDR-borderline" ~
        paste(
          "Adjusted diagnosis association is nominal and close to the",
          "multiple-testing threshold; report as borderline rather than",
          "as a definitive independent association."
        ),
      evidence_class_HC3 ==
        "Nominal only" ~
        paste(
          "Adjusted diagnosis association is nominal only and does not",
          "survive multiple-testing correction."
        ),
      TRUE ~
        paste(
          "No robust adjusted diagnosis evidence after covariate and",
          "multiple-testing control."
        )
    )
  ) %>%
  dplyr::arrange(
    FDR_HC3,
    robust_p_value
  )

safe_write_csv(
  diagnosis_hc3,
  file.path(
    OUTDIR,
    "tables",
    "diagnosis_robustness",
    "adjusted_diagnosis_module_models_HC3.csv"
  )
)

safe_write_csv(
  diagnosis_comparison,
  file.path(
    OUTDIR,
    "tables",
    "diagnosis_robustness",
    "original_OLS_vs_HC3_diagnosis_comparison.csv"
  )
)

###############################################################################
# 14) FIGURES — CORRECTED NESTED SITE RESULTS
###############################################################################

p_nested_site <- nested_site_all %>%
  dplyr::mutate(
    Module = factor(
      Module,
      levels = rev(
        Module[
          order(
            partial_R2_site_within_country,
            na.last = TRUE
          )
        ]
      ),
    ),
    significant = dplyr::case_when(
      FDR_site_within_country <
        0.05 ~ "FDR < 0.05",
      p_value < 0.05 ~
        "Nominal P < 0.05",
      TRUE ~ "Not significant"
    )
  ) %>%
  ggplot(
    aes(
      x = Module,
      y = partial_R2_site_within_country,
      fill = Module
    )
  ) +
  geom_col(
    width = 0.78
  ) +
  geom_text(
    aes(
      label = add_sig_stars(
        FDR_site_within_country
      )
    ),
    hjust = -0.15,
    size = 4
  ) +
  coord_flip(
    clip = "off"
  ) +
  scale_fill_manual(
    values = get_module_colors(
      modules_use
    )
  ) +
  scale_y_continuous(
    labels = scales::percent_format(
      accuracy = 1
    ),
    expand = expansion(
      mult = c(0, 0.12)
    )
  ) +
  labs(
    title = "Corrected site-within-country effects on module eigengenes",
    subtitle = paste0(
      "Incremental partial R² after diagnosis, age, sex, education and country; ",
      "site df supplied by ",
      paste(
        multisite_countries,
        collapse = " and "
      )
    ),
    x = NULL,
    y = "Partial R² for site within country"
  ) +
  theme_bw(
    base_size = 12
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(
      face = "bold"
    ),
    plot.margin = margin(
      5.5,
      25,
      5.5,
      5.5
    )
  )

ggsave(
  file.path(
    OUTDIR,
    "figures",
    "context",
    "corrected_nested_site_partial_R2.pdf"
  ),
  p_nested_site,
  width = 8,
  height = 5.5
)

ggsave(
  file.path(
    OUTDIR,
    "figures",
    "context",
    "corrected_nested_site_partial_R2.png"
  ),
  p_nested_site,
  width = 8,
  height = 5.5,
  dpi = DPI
)

###############################################################################
# 15) FIGURES — LOG-HC3 BIOMARKER RESULTS
###############################################################################

log_hc3_plot_df <- biomarker_robust_models %>%
  dplyr::filter(
    Transform == "log"
  ) %>%
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
    cell_label = dplyr::case_when(
      is.na(
        standardized_beta
      ) ~ "",
      significance_stars == "" ~
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
        significance_stars
      )
    )
  )

p_log_hc3_heatmap <- ggplot(
  log_hc3_plot_df,
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
    size = 3.5
  ) +
  scale_fill_gradient2(
    low = "#4682B4",
    mid = "white",
    high = "#F46D43",
    midpoint = 0,
    name = "Standardized\nβ",
    oob = scales::squish
  ) +
  labs(
    title = "Log-transformed plasma biomarker sensitivity models",
    subtitle = paste0(
      "Outcome ~ module + age + sex + education + categorical country; ",
      "HC3 robust inference and BH-FDR across 32 log-scale models"
    ),
    x = NULL,
    y = NULL
  ) +
  theme_bw(
    base_size = 12
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(
      face = "bold"
    ),
    axis.text.x = element_text(
      angle = 25,
      hjust = 1
    ),
    plot.title = element_text(
      face = "bold"
    )
  )

ggsave(
  file.path(
    OUTDIR,
    "figures",
    "biomarker_robustness",
    "log_HC3_biomarker_standardized_beta_heatmap.pdf"
  ),
  p_log_hc3_heatmap,
  width = 7.5,
  height = 5.5
)

ggsave(
  file.path(
    OUTDIR,
    "figures",
    "biomarker_robustness",
    "log_HC3_biomarker_standardized_beta_heatmap.png"
  ),
  p_log_hc3_heatmap,
  width = 7.5,
  height = 5.5,
  dpi = DPI
)

###############################################################################
# 16) FIGURES — ROBUST DIAGNOSIS RESULTS
###############################################################################

p_diagnosis <- diagnosis_comparison %>%
  dplyr::mutate(
    Module = factor(
      Module,
      levels = rev(
        Module[
          order(
            standardized_difference,
            na.last = TRUE
          )
        ]
      )
    )
  ) %>%
  ggplot(
    aes(
      x = Module,
      y = standardized_difference,
      colour = Module
    )
  ) +
  geom_hline(
    yintercept = 0,
    linetype = 2,
    colour = "grey50"
  ) +
  geom_point(
    size = 3
  ) +
  geom_text(
    aes(
      label = add_sig_stars(
        FDR_HC3
      )
    ),
    hjust = -0.7,
    size = 4
  ) +
  coord_flip(
    clip = "off"
  ) +
  scale_colour_manual(
    values = get_module_colors(
      modules_use
    )
  ) +
  labs(
    title = "Adjusted diagnosis-module associations with HC3 inference",
    subtitle = paste0(
      "Positive values indicate higher eigengene values in clinically diagnosed AD; ",
      "stars indicate BH-FDR"
    ),
    x = NULL,
    y = "Adjusted diagnosis difference (module SD units)"
  ) +
  theme_bw(
    base_size = 12
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(
      face = "bold"
    ),
    plot.margin = margin(
      5.5,
      25,
      5.5,
      5.5
    )
  )

ggsave(
  file.path(
    OUTDIR,
    "figures",
    "diagnosis_robustness",
    "adjusted_diagnosis_HC3_standardized_effects.pdf"
  ),
  p_diagnosis,
  width = 8,
  height = 5.5
)

ggsave(
  file.path(
    OUTDIR,
    "figures",
    "diagnosis_robustness",
    "adjusted_diagnosis_HC3_standardized_effects.png"
  ),
  p_diagnosis,
  width = 8,
  height = 5.5,
  dpi = DPI
)

###############################################################################
# 17) INTEGRATED CORRECTION SUMMARY
###############################################################################

correction_summary_by_module <- module_reference %>%
  dplyr::left_join(
    nested_site_all %>%
      dplyr::select(
        Module,
        corrected_nested_site_p =
          p_value,
        corrected_nested_site_FDR =
          FDR_site_within_country,
        corrected_site_partial_R2 =
          partial_R2_site_within_country,
        corrected_site_added_df =
          added_df,
        corrected_site_model_status =
          model_status
      ),
    by = "Module"
  ) %>%
  dplyr::left_join(
    old_site_comparison %>%
      dplyr::select(
        Module,
        old_non_nested_site_FDR,
        old_non_nested_site_partial_R2,
        partial_R2_attenuation,
        corrected_to_old_partial_R2_ratio
      ),
    by = "Module"
  ) %>%
  dplyr::left_join(
    log_hc3_summary_by_module %>%
      dplyr::select(
        Module,
        n_log_HC3_FDR_lt_0_05,
        minimum_log_HC3_FDR,
        maximum_abs_standardized_beta,
        n_heteroskedasticity_flags
      ),
    by = "Module"
  ) %>%
  dplyr::left_join(
    diagnosis_comparison %>%
      dplyr::select(
        Module,
        diagnosis_HC3_beta =
          estimate_diagnosis,
        diagnosis_HC3_standardized =
          standardized_difference,
        diagnosis_HC3_p =
          robust_p_value,
        diagnosis_HC3_FDR =
          FDR_HC3,
        diagnosis_HC3_evidence =
          evidence_class_HC3
      ),
    by = "Module"
  ) %>%
  dplyr::arrange(
    corrected_nested_site_FDR,
    minimum_log_HC3_FDR
  )

safe_write_csv(
  correction_summary_by_module,
  file.path(
    OUTDIR,
    "tables",
    "integrated_13b_correction_summary_by_module.csv"
  )
)

###############################################################################
# 18) METHODS-READY WORDING
###############################################################################

methods_wording <- tibble::tribble(
  ~section,
  ~text,

  "Methods - site nested within country",
  paste(
    "Because recruitment sites were nested within countries, site-related",
    "variation was evaluated using nested model comparisons rather than by",
    "treating site as an independent cohort-wide factor. For each module",
    "eigengene, a reduced model including diagnosis, age, sex, education and",
    "country was compared with a full model additionally containing the",
    "observed country-by-site factor. The incremental partial R-squared therefore",
    "quantified estimable site-within-country variation. Country-specific site",
    "models were additionally fitted within countries containing at least two",
    "recruitment sites."
  ),

  "Methods - robust biomarker sensitivity",
  paste(
    "Adjusted module-biomarker associations were repeated as sensitivity",
    "analyses using heteroskedasticity-consistent HC3 standard errors. Models",
    "were fitted on both the original biomarker scale and after natural-log",
    "transformation of p-tau181, p-tau217, NfL and the Aβ42/40 ratio. Each",
    "model adjusted for age, sex, education and country as a categorical factor.",
    "False-discovery-rate correction was applied separately across the 32",
    "raw-scale HC3 models and the 32 log-scale HC3 models."
  ),

  "Methods - robust diagnosis sensitivity",
  paste(
    "Diagnosis-related module differences were repeated with HC3 robust",
    "standard errors using module eigengene as the outcome and clinical",
    "diagnosis as the predictor of interest, adjusting for age, sex, education",
    "and country. P values were corrected across the eight modules."
  ),

  "Interpretation - site",
  paste(
    "Corrected site estimates represent internal site-within-country",
    "recruitment-context sensitivity. They do not constitute external validation",
    "and should not be interpreted as evidence of site-specific biology."
  ),

  "Interpretation - diagnosis",
  paste(
    "Adjusted diagnosis effects with nominal P < 0.05 but FDR close to 0.05",
    "should be described as nominal and borderline after multiple-testing",
    "correction rather than as either definitive independent associations or",
    "complete absence of signal."
  )
)

safe_write_csv(
  methods_wording,
  file.path(
    OUTDIR,
    "tables",
    "script13b_methods_and_interpretation_wording.csv"
  )
)

###############################################################################
# 19) EXCEL WORKBOOK
###############################################################################

workbook_tables <- list(
  Analysis_structure =
    context_structure_summary,
  Country_site_counts =
    country_site_counts,
  Sites_per_country =
    sites_per_country,
  Nested_site_all =
    nested_site_all,
  Site_within_country =
    site_within_country,
  Site_by_module =
    site_within_country_summary,
  Old_vs_corrected_site =
    old_site_comparison,
  Biomarker_transform =
    biomarker_transform_audit,
  Biomarker_raw_log_HC3 =
    biomarker_robust_models,
  OLS_vs_HC3 =
    biomarker_model_comparison,
  Log_HC3_by_module =
    log_hc3_summary_by_module,
  Diagnosis_HC3 =
    diagnosis_hc3,
  Diagnosis_comparison =
    diagnosis_comparison,
  Integrated_summary =
    correction_summary_by_module,
  Methods_wording =
    methods_wording
)

openxlsx::write.xlsx(
  workbook_tables,
  file = file.path(
    OUTDIR,
    "WGCNA_13b_Nested_Site_Biomarker_Robustness.xlsx"
  ),
  overwrite = TRUE
)

###############################################################################
# 20) FINAL SUMMARY AND MANIFEST
###############################################################################

script13b_summary <- tibble::tibble(
  metric = c(
    "base_dir",
    "input_script13_dir",
    "output_dir",
    "n_samples",
    "n_modules",
    "modules",
    "site_variable",
    "n_countries",
    "n_site_levels",
    "multisite_countries",
    "expected_site_within_country_df",
    "observed_added_df_min",
    "observed_added_df_max",
    "n_corrected_nested_site_FDR_lt_0_05",
    "n_old_non_nested_site_FDR_lt_0_05",
    "n_raw_HC3_biomarker_FDR_lt_0_05",
    "n_log_HC3_biomarker_FDR_lt_0_05",
    "n_HC3_diagnosis_FDR_lt_0_05",
    "n_HC3_diagnosis_nominal_borderline",
    "robust_vcov_type",
    "primary_biomarker_sensitivity"
  ),
  value = c(
    BASE_DIR,
    SCRIPT13_DIR,
    OUTDIR,
    as.character(nrow(analysis_df)),
    as.character(length(modules_use)),
    paste(
      modules_use,
      collapse = ", "
    ),
    site_var,
    as.character(
      dplyr::n_distinct(
        analysis_df$Country
      )
    ),
    as.character(
      dplyr::n_distinct(
        analysis_df[[site_var]]
      )
    ),
    paste(
      multisite_countries,
      collapse = ", "
    ),
    as.character(
      expected_nested_site_df
    ),
    as.character(
      safe_min(
        nested_site_all$added_df
      )
    ),
    as.character(
      safe_max(
        nested_site_all$added_df
      )
    ),
    as.character(
      sum(
        nested_site_all$
          FDR_site_within_country <
          0.05,
        na.rm = TRUE
      )
    ),
    as.character(
      sum(
        old_site_models$FDR_factor <
          0.05,
        na.rm = TRUE
      )
    ),
    as.character(
      sum(
        biomarker_robust_models$
          Transform == "raw" &
          biomarker_robust_models$
            FDR_within_transform_family <
            0.05,
        na.rm = TRUE
      )
    ),
    as.character(
      sum(
        biomarker_robust_models$
          Transform == "log" &
          biomarker_robust_models$
            FDR_within_transform_family <
            0.05,
        na.rm = TRUE
      )
    ),
    as.character(
      sum(
        diagnosis_hc3$FDR_HC3 <
          0.05,
        na.rm = TRUE
      )
    ),
    as.character(
      sum(
        diagnosis_hc3$
          evidence_class_HC3 ==
          "Nominal and FDR-borderline",
        na.rm = TRUE
      )
    ),
    ROBUST_VCOV_TYPE,
    paste(
      "Natural-log biomarker outcome with HC3 robust standard errors;",
      "age, sex, education and categorical country adjustment."
    )
  )
)

safe_write_csv(
  script13b_summary,
  file.path(
    OUTDIR,
    "tables",
    "script13b_final_summary.csv"
  )
)

output_manifest <- tibble::tibble(
  output_file = c(
    "tables/context/country_site_sample_counts.csv",
    "tables/context/site_nested_structure_audit.csv",
    "tables/context/site_nested_analysis_structure_summary.csv",
    "tables/context/corrected_nested_site_all_samples.csv",
    "tables/context/corrected_site_models_within_multisite_countries.csv",
    "tables/context/corrected_site_within_country_summary_by_module.csv",
    "tables/context/old_non_nested_vs_corrected_nested_site_comparison.csv",
    "tables/biomarker_robustness/biomarker_log_transformation_audit.csv",
    "tables/biomarker_robustness/biomarker_models_raw_and_log_HC3.csv",
    "tables/biomarker_robustness/biomarker_models_log_HC3_primary_sensitivity.csv",
    "tables/biomarker_robustness/biomarker_models_raw_HC3_sensitivity.csv",
    "tables/biomarker_robustness/original_OLS_vs_HC3_model_comparison.csv",
    "tables/biomarker_robustness/log_HC3_summary_by_module.csv",
    "tables/diagnosis_robustness/adjusted_diagnosis_module_models_HC3.csv",
    "tables/diagnosis_robustness/original_OLS_vs_HC3_diagnosis_comparison.csv",
    "tables/integrated_13b_correction_summary_by_module.csv",
    "tables/script13b_methods_and_interpretation_wording.csv",
    "figures/context/corrected_nested_site_partial_R2.pdf/png",
    "figures/biomarker_robustness/log_HC3_biomarker_standardized_beta_heatmap.pdf/png",
    "figures/diagnosis_robustness/adjusted_diagnosis_HC3_standardized_effects.pdf/png",
    "WGCNA_13b_Nested_Site_Biomarker_Robustness.xlsx",
    "workspace/script13b_nested_site_biomarker_robustness_workspace.RData",
    "sessionInfo.txt"
  ),
  description = c(
    "Observed participant counts for every country-site combination.",
    "Number of sites and estimable within-country site degrees of freedom.",
    "Global audit of the nested recruitment structure.",
    "Primary corrected all-sample site-within-country nested models.",
    "Adjusted site models within each country containing at least two sites.",
    "Country-specific site robustness summarized by module.",
    "Direct comparison of invalid non-nested and corrected nested site estimates.",
    "Positivity and skewness audit before and after natural-log transformation.",
    "Complete 64-model raw/log HC3 sensitivity table.",
    "Primary 32-model natural-log HC3 biomarker sensitivity.",
    "Raw-scale HC3 biomarker sensitivity.",
    "Comparison with Script 13 conventional raw-scale OLS inference.",
    "Log-HC3 biomarker robustness summarized by module.",
    "Adjusted diagnosis models with HC3 robust standard errors.",
    "Comparison of conventional OLS and HC3 diagnosis inference.",
    "Integrated nested-site, biomarker and diagnosis correction summary.",
    "Methods-ready wording and interpretation boundaries.",
    "Corrected site-within-country partial R2 figure.",
    "Standardized log-HC3 biomarker association heatmap.",
    "Adjusted robust diagnosis-effect figure.",
    "Integrated reviewer-facing workbook.",
    "Workspace for Script 14 robustness analyses.",
    "R session information."
  )
)

safe_write_csv(
  output_manifest,
  file.path(
    OUTDIR,
    "script13b_output_manifest.csv"
  )
)

###############################################################################
# 21) SAVE WORKSPACE
###############################################################################

save(
  analysis_df,
  module_reference,
  modules_use,
  site_var,
  country_site_counts,
  sites_per_country,
  multisite_countries,
  expected_nested_site_df,
  context_structure_summary,
  nested_site_all,
  site_within_country,
  site_within_country_summary,
  old_site_comparison,
  biomarker_transform_audit,
  biomarker_robust_models,
  biomarker_model_comparison,
  log_hc3_summary_by_module,
  diagnosis_hc3,
  diagnosis_comparison,
  correction_summary_by_module,
  methods_wording,
  script13b_summary,
  output_manifest,
  file = file.path(
    OUTDIR,
    "workspace",
    "script13b_nested_site_biomarker_robustness_workspace.RData"
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

cat("\nScript 13b finished successfully.\n")
cat("Main output directory:\n", OUTDIR, "\n")
cat("Samples analyzed: ", nrow(analysis_df), "\n", sep = "")
cat("Modules analyzed: ", length(modules_use), "\n", sep = "")
cat("Site variable: ", site_var, "\n", sep = "")
cat(
  "Multisite countries: ",
  paste(
    multisite_countries,
    collapse = ", "
  ),
  "\n",
  sep = ""
)
cat(
  "Expected site-within-country added df: ",
  expected_nested_site_df,
  "\n",
  sep = ""
)
cat(
  "Observed added df range: ",
  safe_min(nested_site_all$added_df),
  " to ",
  safe_max(nested_site_all$added_df),
  "\n",
  sep = ""
)
cat(
  "Corrected nested-site models at FDR < 0.05: ",
  sum(
    nested_site_all$
      FDR_site_within_country <
      0.05,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)
cat(
  "Raw-HC3 biomarker models at FDR < 0.05: ",
  sum(
    biomarker_robust_models$
      Transform == "raw" &
      biomarker_robust_models$
        FDR_within_transform_family <
        0.05,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)
cat(
  "Log-HC3 biomarker models at FDR < 0.05: ",
  sum(
    biomarker_robust_models$
      Transform == "log" &
      biomarker_robust_models$
        FDR_within_transform_family <
        0.05,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)
cat(
  "HC3 diagnosis models at FDR < 0.05: ",
  sum(
    diagnosis_hc3$FDR_HC3 <
      0.05,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)
cat(
  "HC3 diagnosis models nominal and FDR-borderline: ",
  sum(
    diagnosis_hc3$
      evidence_class_HC3 ==
      "Nominal and FDR-borderline",
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)
cat("\nKey outputs:\n")
cat("- tables/context/corrected_nested_site_all_samples.csv\n")
cat("- tables/context/corrected_site_models_within_multisite_countries.csv\n")
cat("- tables/context/old_non_nested_vs_corrected_nested_site_comparison.csv\n")
cat("- tables/biomarker_robustness/biomarker_models_log_HC3_primary_sensitivity.csv\n")
cat("- tables/biomarker_robustness/original_OLS_vs_HC3_model_comparison.csv\n")
cat("- tables/diagnosis_robustness/original_OLS_vs_HC3_diagnosis_comparison.csv\n")
cat("- tables/integrated_13b_correction_summary_by_module.csv\n")
cat("- WGCNA_13b_Nested_Site_Biomarker_Robustness.xlsx\n")

###############################################################################
# END
###############################################################################


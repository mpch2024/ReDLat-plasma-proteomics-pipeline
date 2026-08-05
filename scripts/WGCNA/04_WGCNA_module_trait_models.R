###############################################################################
# ReDLat plasma proteomics — WGCNA workflow
# 04. Test module–trait associations
# Requires: outputs from Scripts 01–03
# Produces: correlations, adjusted models and context summaries
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
  "pheatmap",
  "ggrepel",
  "scales",
  "broom",
  "forcats",
  "openxlsx"
)

cran_missing <- cran_pkgs[
  !vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)
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

SCRIPT10_DIR <- file.path(WGCNA_CONFIG$result_root,
  "01_input"
)

SCRIPT11_DIR <- file.path(WGCNA_CONFIG$result_root,
  "02_network"
)

SCRIPT12_DIR <- file.path(WGCNA_CONFIG$result_root,
  "03_modules"
)

META_FILE <- file.path(
  SCRIPT10_DIR,
  "wgcna_sample_metadata.csv"
)

EIGENGENE_FILE <- file.path(
  SCRIPT11_DIR,
  "eigengenes",
  "module_eigengenes_per_sample.csv"
)

INTEGRATED_BIOLOGY_FILE <- file.path(
  SCRIPT12_DIR,
  "tables",
  "integrated_module_biology_summary.csv"
)

HUB_SUMMARY_FILE <- file.path(
  SCRIPT12_DIR,
  "tables",
  "module_hub_summary.csv"
)

DEP_BURDEN_FILE <- file.path(
  SCRIPT12_DIR,
  "tables",
  "module_DEP_burden_dual_definition.csv"
)

DEP_OVERREP_FILE <- file.path(
  SCRIPT12_DIR,
  "tables",
  "module_DEP_overrepresentation_dual_definition.csv"
)

ENRICHMENT_SUMMARY_FILE <- file.path(
  SCRIPT12_DIR,
  "tables",
  "enrichment_summary_by_module.csv"
)

TOP_HUBS_FILE <- file.path(
  SCRIPT12_DIR,
  "tables",
  "top_hubs_all_modules_combined.csv"
)

OUTDIR <- file.path(WGCNA_CONFIG$result_root,
  "04_module_traits"
)

SUBDIRS <- c(
  "tables",
  "tables/correlations",
  "tables/regression",
  "tables/context",
  "tables/prioritization",
  "figures",
  "figures/module_trait",
  "figures/scatter",
  "figures/context",
  "figures/prioritization",
  "workspace"
)

invisible(lapply(
  file.path(OUTDIR, SUBDIRS),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

required_files <- c(
  META_FILE,
  EIGENGENE_FILE,
  INTEGRATED_BIOLOGY_FILE,
  HUB_SUMMARY_FILE,
  DEP_BURDEN_FILE,
  DEP_OVERREP_FILE,
  ENRICHMENT_SUMMARY_FILE
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  stop(
    "Missing required input files:\n",
    paste(missing_files, collapse = "\n"),
    "\n\nRun Scripts 10, 11 and 12 first.",
    call. = FALSE
  )
}

###############################################################################
# 3) PARAMETERS
###############################################################################

MIN_N_FOR_CORR <- 6L
MIN_N_FOR_REG <- 15L
MIN_N_PER_CONTEXT_LEVEL_DESCRIPTIVE <- 1L

EXPECTED_N_SAMPLES <- 639L
EXPECTED_N_MODULES <- 8L
REQUIRE_EXPECTED_DIMENSIONS <- TRUE

MODULES_OF_INTEREST <- NULL
EXCLUDE_MODULES <- c("grey", "MEgrey")

# Country_numeric is deliberately absent.
TRAIT_ORDER <- c(
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

MAIN_TRAITS_FOR_INTEGRATION <- c(
  "SampleGroup_bin",
  "cdr_boxscore",
  "mmse_total",
  "udsfaq_total",
  "p_tau181",
  "p_tau217",
  "NfL",
  "ratio_AB42_40",
  "Age"
)

SCATTER_PRIORITY_TRAITS <- c(
  "SampleGroup_bin",
  "cdr_boxscore",
  "mmse_total",
  "udsfaq_total",
  "p_tau181",
  "p_tau217",
  "NfL",
  "ratio_AB42_40",
  "Age"
)

TOP_SCATTERS_PER_MODULE <- 4L

# Original continuous clinical/biomarker outcomes retained.
REGRESSION_OUTCOMES <- c(
  "p_tau181",
  "p_tau217",
  "NfL",
  "ratio_AB42_40",
  "cdr_boxscore",
  "mmse_total"
)

# Adjusted model covariates.
# Country is included as a factor in clinical models.
CLINICAL_MODEL_COVARS <- c(
  "Age",
  "Sex_bin",
  "Education",
  "Country"
)

# Country/site nested models adjust for diagnosis and demographics.
CONTEXT_MODEL_COVARS <- c(
  "SampleGroup_bin",
  "Age",
  "Sex_bin",
  "Education"
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

# Provisional biological names based on the audited Script 12 outputs.
# These labels are descriptive and can be updated without rerunning WGCNA.
MODULE_BIOLOGICAL_LABELS <- c(
  blue = "AD-elevated high-differential-burden",
  green = "Neuronal connectivity and extracellular matrix",
  brown = "RNA processing and ribonucleoprotein",
  magenta = "Broad intracellular trafficking and proteostasis",
  pink = "Limited enrichment / compact module",
  purple = "Limited enrichment / compact module",
  black = "Limited enrichment / JAK-STAT-related signal",
  red = "Limited enrichment / compact module"
)

DPI <- 300L
HEATMAP_WIDTH <- 12
HEATMAP_HEIGHT <- 7

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

safe_read_optional_csv <- function(file) {
  if (!file.exists(file)) {
    return(NULL)
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

standardize_module_name <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  gsub("^ME", "", x)
}

get_module_colors <- function(modules) {
  modules <- as.character(modules)
  cols <- MODULE_COLORS[modules]

  missing_mods <- modules[is.na(cols)]

  if (length(missing_mods) > 0) {
    extra_cols <- grDevices::rainbow(
      length(unique(missing_mods)),
      s = 0.45,
      v = 0.75
    )

    names(extra_cols) <- unique(missing_mods)
    cols[is.na(cols)] <- extra_cols[missing_mods]
  }

  names(cols) <- modules
  cols
}

get_biological_label <- function(module) {
  module <- as.character(module)
  out <- MODULE_BIOLOGICAL_LABELS[module]
  out[is.na(out)] <- "Biological label pending"
  unname(out)
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

coalesce_numeric_aliases <- function(df, target, candidates) {
  present <- unique(
    c(target, candidates)[
      c(target, candidates) %in% names(df)
    ]
  )

  if (length(present) == 0) {
    df[[target]] <- NA_real_
    return(df)
  }

  values <- lapply(
    present,
    function(cc) safe_numeric(df[[cc]])
  )

  df[[target]] <- Reduce(
    dplyr::coalesce,
    values
  )

  df
}

standardize_metadata_columns <- function(df) {
  df <- df %>%
    dplyr::mutate(
      SampleId = as.character(SampleId)
    )

  df <- coalesce_numeric_aliases(
    df,
    "Mini_SEA",
    c("Mini-SEA", "Mini.SEA")
  )

  df <- coalesce_numeric_aliases(
    df,
    "T_ADLQ",
    c("T-ADLQ", "T.ADLQ")
  )

  df <- coalesce_numeric_aliases(
    df,
    "p_tau181",
    c("p-tau181", "p.tau181")
  )

  df <- coalesce_numeric_aliases(
    df,
    "p_tau217",
    c("p-tau217", "p.tau217")
  )

  df <- coalesce_numeric_aliases(
    df,
    "ratio_AB42_40",
    c("ratio AB42/40", "ratio.AB42.40")
  )

  df
}

encode_basic_variables <- function(df) {
  if ("SampleGroup" %in% names(df)) {
    df <- df %>%
      dplyr::mutate(
        SampleGroup = as.character(SampleGroup),
        SampleGroup_bin = dplyr::case_when(
          SampleGroup == "CN" ~ 0,
          SampleGroup == "AD" ~ 1,
          TRUE ~ NA_real_
        )
      )
  }

  if ("Sex" %in% names(df)) {
    df <- df %>%
      dplyr::mutate(
        Sex_chr = as.character(Sex),
        Sex_bin = dplyr::case_when(
          Sex_chr %in% c(
            "1",
            "M",
            "Male",
            "male",
            "Hombre"
          ) ~ 0,
          Sex_chr %in% c(
            "2",
            "F",
            "Female",
            "female",
            "Mujer"
          ) ~ 1,
          TRUE ~ safe_numeric(Sex_chr)
        )
      )
  }

  if (
    "ApoE" %in% names(df) &&
    !"APOE_group" %in% names(df)
  ) {
    df <- df %>%
      dplyr::mutate(
        ApoE = as.character(ApoE),
        APOE_group = dplyr::case_when(
          ApoE %in% c(
            "e2/e4",
            "e3/e4",
            "e4/e4"
          ) ~ "E4 carrier",
          ApoE %in% c(
            "e2/e2",
            "e2/e3",
            "e3/e3"
          ) ~ "Non-E4",
          TRUE ~ NA_character_
        )
      )
  }

  if (
    "APOE_group" %in% names(df) &&
    !"APOE4_carrier" %in% names(df)
  ) {
    df <- df %>%
      dplyr::mutate(
        APOE4_carrier = dplyr::case_when(
          APOE_group == "E4 carrier" ~ 1,
          APOE_group == "Non-E4" ~ 0,
          TRUE ~ NA_real_
        )
      )
  }

  numeric_candidates <- c(
    "Age",
    "Education",
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
    "APOE4_carrier",
    "SampleGroup_bin",
    "Sex_bin"
  )

  for (cc in intersect(
    numeric_candidates,
    names(df)
  )) {
    df[[cc]] <- safe_numeric(df[[cc]])
  }

  if ("Country" %in% names(df)) {
    df <- df %>%
      dplyr::mutate(
        Country = factor(
          clean_text_na(Country)
        )
      )
  }

  # Explicitly remove any inherited numeric encoding.
  if ("Country_numeric" %in% names(df)) {
    df$Country_numeric <- NULL
  }

  df
}

detect_site_variable <- function(df) {
  candidates <- c(
    "Site",
    "site",
    "Center",
    "center",
    "Cohort",
    "cohort",
    "RecruitmentSite",
    "recruitment_site",
    "site_id",
    "Site_ID"
  )

  hit <- candidates[
    candidates %in% names(df)
  ][1]

  if (
    length(hit) == 0 ||
    is.na(hit)
  ) {
    return(NA_character_)
  }

  hit
}

safe_spearman <- function(
    df,
    xcol,
    ycol,
    min_n = MIN_N_FOR_CORR
) {
  d <- df %>%
    dplyr::select(
      dplyr::all_of(c(
        xcol,
        ycol
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
    dplyr::n_distinct(d[[xcol]]) <= 1 ||
    dplyr::n_distinct(d[[ycol]]) <= 1
  ) {
    return(tibble::tibble(
      rho = NA_real_,
      p_value = NA_real_,
      N = n
    ))
  }

  test <- suppressWarnings(
    stats::cor.test(
      d[[xcol]],
      d[[ycol]],
      method = "spearman",
      exact = FALSE
    )
  )

  tibble::tibble(
    rho = unname(test$estimate),
    p_value = test$p.value,
    N = n
  )
}

compute_module_trait <- function(
    df,
    modules,
    traits,
    label = "full"
) {
  expand.grid(
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
      Analysis = label
    )
}

safe_kruskal <- function(
    df,
    outcome,
    factor_var
) {
  d <- df %>%
    dplyr::select(
      dplyr::all_of(c(
        outcome,
        factor_var
      ))
    ) %>%
    tidyr::drop_na()

  if (nrow(d) < MIN_N_FOR_CORR) {
    return(tibble::tibble(
      statistic = NA_real_,
      p_value = NA_real_,
      n = nrow(d),
      n_levels = NA_integer_
    ))
  }

  d[[factor_var]] <- droplevels(
    factor(d[[factor_var]])
  )

  n_levels <- nlevels(
    d[[factor_var]]
  )

  if (n_levels < 2) {
    return(tibble::tibble(
      statistic = NA_real_,
      p_value = NA_real_,
      n = nrow(d),
      n_levels = n_levels
    ))
  }

  kt <- tryCatch(
    stats::kruskal.test(
      stats::as.formula(
        paste0(
          outcome,
          " ~ ",
          factor_var
        )
      ),
      data = d
    ),
    error = function(e) NULL
  )

  if (is.null(kt)) {
    return(tibble::tibble(
      statistic = NA_real_,
      p_value = NA_real_,
      n = nrow(d),
      n_levels = n_levels
    ))
  }

  tibble::tibble(
    statistic = unname(kt$statistic),
    p_value = kt$p.value,
    n = nrow(d),
    n_levels = n_levels
  )
}

safe_adjusted_factor_lm <- function(
    df,
    outcome,
    factor_var,
    covars = CONTEXT_MODEL_COVARS
) {
  available_covars <- covars[
    covars %in% names(df)
  ]

  model_vars <- unique(
    c(
      outcome,
      factor_var,
      available_covars
    )
  )

  d <- df %>%
    dplyr::select(
      dplyr::all_of(model_vars)
    ) %>%
    tidyr::drop_na()

  if (nrow(d) < MIN_N_FOR_REG) {
    return(tibble::tibble(
      n = nrow(d),
      n_levels = NA_integer_,
      residual_df_full = NA_real_,
      covariates = paste(
        available_covars,
        collapse = ", "
      ),
      p_value_factor = NA_real_,
      partial_R2_factor = NA_real_,
      full_model_R2 = NA_real_,
      full_model_adj_R2 = NA_real_,
      model_status = "insufficient_n"
    ))
  }

  d[[factor_var]] <- droplevels(
    factor(d[[factor_var]])
  )

  n_levels <- nlevels(
    d[[factor_var]]
  )

  if (n_levels < 2) {
    return(tibble::tibble(
      n = nrow(d),
      n_levels = n_levels,
      residual_df_full = NA_real_,
      covariates = paste(
        available_covars,
        collapse = ", "
      ),
      p_value_factor = NA_real_,
      partial_R2_factor = NA_real_,
      full_model_R2 = NA_real_,
      full_model_adj_R2 = NA_real_,
      model_status = "fewer_than_two_levels"
    ))
  }

  rhs_full <- paste(
    c(
      factor_var,
      available_covars
    ),
    collapse = " + "
  )

  rhs_reduced <- if (
    length(available_covars) > 0
  ) {
    paste(
      available_covars,
      collapse = " + "
    )
  } else {
    "1"
  }

  form_full <- stats::as.formula(
    paste0(
      outcome,
      " ~ ",
      rhs_full
    )
  )

  form_reduced <- stats::as.formula(
    paste0(
      outcome,
      " ~ ",
      rhs_reduced
    )
  )

  fit_full <- tryCatch(
    stats::lm(
      form_full,
      data = d
    ),
    error = function(e) NULL
  )

  fit_reduced <- tryCatch(
    stats::lm(
      form_reduced,
      data = d
    ),
    error = function(e) NULL
  )

  if (
    is.null(fit_full) ||
    is.null(fit_reduced)
  ) {
    return(tibble::tibble(
      n = nrow(d),
      n_levels = n_levels,
      residual_df_full = NA_real_,
      covariates = paste(
        available_covars,
        collapse = ", "
      ),
      p_value_factor = NA_real_,
      partial_R2_factor = NA_real_,
      full_model_R2 = NA_real_,
      full_model_adj_R2 = NA_real_,
      model_status = "model_failed"
    ))
  }

  cmp <- tryCatch(
    stats::anova(
      fit_reduced,
      fit_full
    ),
    error = function(e) NULL
  )

  if (
    is.null(cmp) ||
    nrow(cmp) < 2
  ) {
    p_val <- NA_real_
    partial_r2 <- NA_real_
  } else {
    p_val <- cmp$`Pr(>F)`[2]

    rss_reduced <- cmp$RSS[1]
    rss_full <- cmp$RSS[2]

    partial_r2 <- (
      rss_reduced - rss_full
    ) / rss_reduced
  }

  fit_summary <- summary(fit_full)

  tibble::tibble(
    n = nrow(d),
    n_levels = n_levels,
    residual_df_full = stats::df.residual(
      fit_full
    ),
    covariates = paste(
      available_covars,
      collapse = ", "
    ),
    p_value_factor = p_val,
    partial_R2_factor = partial_r2,
    full_model_R2 = fit_summary$r.squared,
    full_model_adj_R2 = fit_summary$adj.r.squared,
    model_status = "ok"
  )
}

make_factor_descriptives <- function(
    df,
    modules,
    factor_var
) {
  purrr::map_dfr(
    modules,
    function(m) {
      df %>%
        dplyr::select(
          dplyr::all_of(c(
            m,
            factor_var
          ))
        ) %>%
        tidyr::drop_na() %>%
        dplyr::mutate(
          Level = as.character(
            .data[[factor_var]]
          )
        ) %>%
        dplyr::group_by(Level) %>%
        dplyr::summarise(
          Module = m,
          n = dplyr::n(),
          mean_eigengene = mean(
            .data[[m]],
            na.rm = TRUE
          ),
          sd_eigengene = stats::sd(
            .data[[m]],
            na.rm = TRUE
          ),
          median_eigengene = stats::median(
            .data[[m]],
            na.rm = TRUE
          ),
          q25_eigengene = stats::quantile(
            .data[[m]],
            0.25,
            na.rm = TRUE
          ),
          q75_eigengene = stats::quantile(
            .data[[m]],
            0.75,
            na.rm = TRUE
          ),
          .groups = "drop"
        ) %>%
        dplyr::select(
          Module,
          Level,
          dplyr::everything()
        )
    }
  )
}

safe_continuous_module_model <- function(
    df,
    outcome,
    module,
    covars = CLINICAL_MODEL_COVARS
) {
  available_covars <- covars[
    covars %in% names(df)
  ]

  model_vars <- unique(
    c(
      outcome,
      module,
      available_covars
    )
  )

  d <- df %>%
    dplyr::select(
      dplyr::all_of(model_vars)
    ) %>%
    tidyr::drop_na()

  if (
    nrow(d) < MIN_N_FOR_REG ||
    dplyr::n_distinct(d[[outcome]]) <= 1 ||
    dplyr::n_distinct(d[[module]]) <= 1
  ) {
    return(tibble::tibble(
      Module = module,
      Outcome = outcome,
      N = nrow(d),
      estimate = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_,
      adj_r_squared = NA_real_,
      Formula = NA_character_,
      model_status = "insufficient_or_constant"
    ))
  }

  if ("Country" %in% names(d)) {
    d$Country <- droplevels(
      factor(d$Country)
    )
  }

  form_txt <- paste0(
    outcome,
    " ~ ",
    module,
    if (length(available_covars) > 0) {
      paste0(
        " + ",
        paste(
          available_covars,
          collapse = " + "
        )
      )
    } else {
      ""
    }
  )

  fit <- tryCatch(
    stats::lm(
      stats::as.formula(form_txt),
      data = d
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(tibble::tibble(
      Module = module,
      Outcome = outcome,
      N = nrow(d),
      estimate = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_,
      adj_r_squared = NA_real_,
      Formula = form_txt,
      model_status = "model_failed"
    ))
  }

  term_tbl <- broom::tidy(
    fit,
    conf.int = TRUE
  ) %>%
    dplyr::filter(
      term == module
    )

  if (nrow(term_tbl) != 1) {
    return(tibble::tibble(
      Module = module,
      Outcome = outcome,
      N = nrow(d),
      estimate = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_,
      adj_r_squared = summary(fit)$adj.r.squared,
      Formula = form_txt,
      model_status = "module_term_not_found"
    ))
  }

  term_tbl %>%
    dplyr::transmute(
      Module = module,
      Outcome = outcome,
      N = nrow(d),
      estimate,
      std.error,
      statistic,
      p.value,
      conf.low,
      conf.high,
      adj_r_squared = summary(fit)$adj.r.squared,
      Formula = form_txt,
      model_status = "ok"
    )
}

safe_diagnosis_module_model <- function(
    df,
    module
) {
  required <- c(
    module,
    "SampleGroup_bin"
  )

  if (!all(required %in% names(df))) {
    return(tibble::tibble(
      Module = module,
      N = 0L,
      estimate_diagnosis = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_,
      adj_r_squared = NA_real_,
      Formula = NA_character_,
      model_status = "required_variables_missing"
    ))
  }

  covars <- intersect(
    CLINICAL_MODEL_COVARS,
    names(df)
  )

  model_vars <- unique(
    c(
      module,
      "SampleGroup_bin",
      covars
    )
  )

  d <- df %>%
    dplyr::select(
      dplyr::all_of(model_vars)
    ) %>%
    tidyr::drop_na()

  if (
    nrow(d) < MIN_N_FOR_REG ||
    dplyr::n_distinct(d$SampleGroup_bin) < 2
  ) {
    return(tibble::tibble(
      Module = module,
      N = nrow(d),
      estimate_diagnosis = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_,
      adj_r_squared = NA_real_,
      Formula = NA_character_,
      model_status = "insufficient_or_constant"
    ))
  }

  if ("Country" %in% names(d)) {
    d$Country <- droplevels(
      factor(d$Country)
    )
  }

  form_txt <- paste0(
    module,
    " ~ SampleGroup_bin",
    if (length(covars) > 0) {
      paste0(
        " + ",
        paste(
          covars,
          collapse = " + "
        )
      )
    } else {
      ""
    }
  )

  fit <- tryCatch(
    stats::lm(
      stats::as.formula(form_txt),
      data = d
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(tibble::tibble(
      Module = module,
      N = nrow(d),
      estimate_diagnosis = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_,
      adj_r_squared = NA_real_,
      Formula = form_txt,
      model_status = "model_failed"
    ))
  }

  term_tbl <- broom::tidy(
    fit,
    conf.int = TRUE
  ) %>%
    dplyr::filter(
      term == "SampleGroup_bin"
    )

  if (nrow(term_tbl) != 1) {
    return(tibble::tibble(
      Module = module,
      N = nrow(d),
      estimate_diagnosis = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_,
      adj_r_squared = summary(fit)$adj.r.squared,
      Formula = form_txt,
      model_status = "diagnosis_term_not_found"
    ))
  }

  term_tbl %>%
    dplyr::transmute(
      Module = module,
      N = nrow(d),
      estimate_diagnosis = estimate,
      std.error,
      statistic,
      p.value,
      conf.low,
      conf.high,
      adj_r_squared = summary(fit)$adj.r.squared,
      Formula = form_txt,
      model_status = "ok"
    )
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

percentile_rank_safe <- function(x) {
  x_num <- safe_numeric(x)
  out <- rep(NA_real_, length(x_num))
  ok <- is.finite(x_num)

  if (!any(ok)) {
    return(out)
  }

  if (dplyr::n_distinct(x_num[ok]) == 1) {
    out[ok] <- 0.5
    return(out)
  }

  out[ok] <- dplyr::percent_rank(
    x_num[ok]
  )

  out
}

###############################################################################
# 5) LOAD DATA
###############################################################################

eig <- safe_read_csv(
  EIGENGENE_FILE
)

meta <- safe_read_csv(
  META_FILE
)

integrated_biology <- safe_read_csv(
  INTEGRATED_BIOLOGY_FILE
)

hub_summary <- safe_read_csv(
  HUB_SUMMARY_FILE
)

dep_burden <- safe_read_csv(
  DEP_BURDEN_FILE
)

dep_overrep <- safe_read_csv(
  DEP_OVERREP_FILE
)

enrichment_summary <- safe_read_csv(
  ENRICHMENT_SUMMARY_FILE
)

top_hubs <- safe_read_optional_csv(
  TOP_HUBS_FILE
)

cat("INPUT DIMENSIONS\n")
cat("eigengenes             :", dim(eig), "\n")
cat("metadata               :", dim(meta), "\n")
cat("integrated biology     :", dim(integrated_biology), "\n")
cat("hub summary            :", dim(hub_summary), "\n")
cat("DEP burden             :", dim(dep_burden), "\n")
cat("DEP overrepresentation :", dim(dep_overrep), "\n")
cat("enrichment summary     :", dim(enrichment_summary), "\n\n")

###############################################################################
# 6) PREPARE METADATA AND EIGENGENES
###############################################################################

if (!"SampleId" %in% names(eig)) {
  stop(
    "Eigengene table does not contain SampleId.",
    call. = FALSE
  )
}

if (!"SampleId" %in% names(meta)) {
  stop(
    "Metadata does not contain SampleId.",
    call. = FALSE
  )
}

eig <- eig %>%
  dplyr::mutate(
    SampleId = as.character(SampleId)
  ) %>%
  dplyr::distinct(
    SampleId,
    .keep_all = TRUE
  )

meta <- meta %>%
  standardize_metadata_columns() %>%
  encode_basic_variables() %>%
  dplyr::distinct(
    SampleId,
    .keep_all = TRUE
  )

if (
  anyDuplicated(eig$SampleId) > 0 ||
  anyDuplicated(meta$SampleId) > 0
) {
  stop(
    "Duplicated SampleId values remained after cleaning.",
    call. = FALSE
  )
}

if (!setequal(
  eig$SampleId,
  meta$SampleId
)) {
  stop(
    "Eigengene and metadata files do not contain the same SampleId set.",
    call. = FALSE
  )
}

meta <- meta[
  match(
    eig$SampleId,
    meta$SampleId
  ),
  ,
  drop = FALSE
]

if (!all(
  eig$SampleId ==
    meta$SampleId
)) {
  stop(
    "Eigengene and metadata sample order could not be aligned.",
    call. = FALSE
  )
}

eig_cols_original <- setdiff(
  names(eig),
  "SampleId"
)

eig_cols_original <- eig_cols_original[
  !eig_cols_original %in%
    EXCLUDE_MODULES
]

eig_clean <- eig

names(eig_clean)[
  names(eig_clean) %in%
    eig_cols_original
] <- standardize_module_name(
  eig_cols_original
)

module_cols <- setdiff(
  names(eig_clean),
  "SampleId"
)

module_cols <- module_cols[
  vapply(
    eig_clean[module_cols],
    is.numeric,
    logical(1)
  )
]

module_cols <- setdiff(
  module_cols,
  standardize_module_name(
    EXCLUDE_MODULES
  )
)

if (is.null(MODULES_OF_INTEREST)) {
  modules_use <- module_cols
} else {
  modules_use <- intersect(
    MODULES_OF_INTEREST,
    module_cols
  )
}

if (length(modules_use) == 0) {
  stop(
    "No valid module eigengene columns were detected.",
    call. = FALSE
  )
}

analysis_df <- eig_clean %>%
  dplyr::inner_join(
    meta,
    by = "SampleId"
  )

if (nrow(analysis_df) == 0) {
  stop(
    "No shared samples remained between eigengenes and metadata.",
    call. = FALSE
  )
}

for (mc in modules_use) {
  analysis_df[[mc]] <- safe_numeric(
    analysis_df[[mc]]
  )
}

if ("Country_numeric" %in% names(analysis_df)) {
  analysis_df$Country_numeric <- NULL
}

if (!"Country" %in% names(analysis_df)) {
  stop(
    "Country was not found in metadata.",
    call. = FALSE
  )
}

analysis_df$Country <- droplevels(
  factor(
    clean_text_na(
      analysis_df$Country
    )
  )
)

site_var <- detect_site_variable(
  analysis_df
)

if (
  !is.na(site_var) &&
  site_var %in% names(analysis_df)
) {
  analysis_df[[site_var]] <- droplevels(
    factor(
      clean_text_na(
        analysis_df[[site_var]]
      )
    )
  )
}

traits_available <- intersect(
  TRAIT_ORDER,
  names(analysis_df)
)

if (length(traits_available) == 0) {
  stop(
    "No eligible traits were available for module-trait correlations.",
    call. = FALSE
  )
}

if ("Country_numeric" %in% traits_available) {
  stop(
    "Country_numeric was incorrectly included among correlation traits.",
    call. = FALSE
  )
}

if (REQUIRE_EXPECTED_DIMENSIONS) {
  if (nrow(analysis_df) != EXPECTED_N_SAMPLES) {
    stop(
      "Expected ",
      EXPECTED_N_SAMPLES,
      " aligned samples, but found ",
      nrow(analysis_df),
      ".",
      call. = FALSE
    )
  }

  if (length(modules_use) != EXPECTED_N_MODULES) {
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

module_reference <- tibble::tibble(
  Module = modules_use,
  Module_color = unname(
    get_module_colors(
      modules_use
    )
  ),
  Biological_label = get_biological_label(
    modules_use
  )
)

safe_write_csv(
  analysis_df,
  file.path(
    OUTDIR,
    "tables",
    "module_trait_input_clean.csv"
  )
)

safe_write_csv(
  module_reference,
  file.path(
    OUTDIR,
    "tables",
    "module_biological_label_reference.csv"
  )
)

input_audit <- tibble::tibble(
  metric = c(
    "base_dir",
    "n_samples_eigengenes",
    "n_samples_metadata",
    "n_samples_analysis",
    "n_modules",
    "modules",
    "n_traits_available",
    "traits_available",
    "country_levels",
    "site_variable_detected",
    "site_levels",
    "Country_numeric_present_in_analysis_df",
    "Country_numeric_present_in_traits",
    "wgcna_input_level"
  ),
  value = c(
    BASE_DIR,
    as.character(nrow(eig)),
    as.character(nrow(meta)),
    as.character(nrow(analysis_df)),
    as.character(length(modules_use)),
    paste(
      modules_use,
      collapse = ", "
    ),
    as.character(length(traits_available)),
    paste(
      traits_available,
      collapse = ", "
    ),
    paste(
      levels(
        droplevels(
          analysis_df$Country
        )
      ),
      collapse = ", "
    ),
    ifelse(
      is.na(site_var),
      "No",
      site_var
    ),
    ifelse(
      is.na(site_var),
      "Not available",
      paste(
        levels(
          droplevels(
            analysis_df[[site_var]]
          )
        ),
        collapse = ", "
      )
    ),
    as.character(
      "Country_numeric" %in%
        names(analysis_df)
    ),
    as.character(
      "Country_numeric" %in%
        traits_available
    ),
    "GENE-COLLAPSED, outcome-independent SOMAmer selection"
  )
)

safe_write_csv(
  input_audit,
  file.path(
    OUTDIR,
    "tables",
    "script13_input_alignment_audit.csv"
  )
)

cat("Samples in analysis:", nrow(analysis_df), "\n")
cat("Modules analyzed:", paste(modules_use, collapse = ", "), "\n")
cat("Traits available:", paste(traits_available, collapse = ", "), "\n")
cat(
  "Site variable:",
  ifelse(is.na(site_var), "not detected", site_var),
  "\n\n"
)

###############################################################################
# 7) MODULE-TRAIT CORRELATIONS
###############################################################################

module_trait_long <- compute_module_trait(
  df = analysis_df,
  modules = modules_use,
  traits = traits_available,
  label = "full"
) %>%
  dplyr::mutate(
    FDR = p.adjust(
      p_value,
      method = "BH"
    ),
    stars = add_sig_stars(
      FDR
    ),
    Trait_label = dplyr::recode(
      Trait,
      !!!TRAIT_LABELS,
      .default = Trait
    ),
    Biological_label = get_biological_label(
      Module
    )
  ) %>%
  dplyr::arrange(
    FDR,
    p_value,
    dplyr::desc(abs(rho))
  )

expected_test_count <- length(
  modules_use
) * length(
  traits_available
)

if (nrow(module_trait_long) != expected_test_count) {
  stop(
    "Expected ",
    expected_test_count,
    " module-trait tests, but obtained ",
    nrow(module_trait_long),
    ".",
    call. = FALSE
  )
}

safe_write_csv(
  module_trait_long,
  file.path(
    OUTDIR,
    "module_trait_results_long.csv"
  )
)

safe_write_csv(
  module_trait_long,
  file.path(
    OUTDIR,
    "tables",
    "correlations",
    "module_trait_results_long.csv"
  )
)

rho_mat <- module_trait_long %>%
  dplyr::select(
    Module,
    Trait,
    rho
  ) %>%
  tidyr::pivot_wider(
    names_from = Trait,
    values_from = rho
  ) %>%
  tibble::column_to_rownames(
    "Module"
  ) %>%
  as.matrix()

p_mat <- module_trait_long %>%
  dplyr::select(
    Module,
    Trait,
    p_value
  ) %>%
  tidyr::pivot_wider(
    names_from = Trait,
    values_from = p_value
  ) %>%
  tibble::column_to_rownames(
    "Module"
  ) %>%
  as.matrix()

fdr_mat <- module_trait_long %>%
  dplyr::select(
    Module,
    Trait,
    FDR
  ) %>%
  tidyr::pivot_wider(
    names_from = Trait,
    values_from = FDR
  ) %>%
  tibble::column_to_rownames(
    "Module"
  ) %>%
  as.matrix()

annot_mat <- module_trait_long %>%
  dplyr::select(
    Module,
    Trait,
    stars
  ) %>%
  tidyr::pivot_wider(
    names_from = Trait,
    values_from = stars
  ) %>%
  tibble::column_to_rownames(
    "Module"
  ) %>%
  as.matrix()

safe_write_csv(
  as.data.frame(rho_mat) %>%
    tibble::rownames_to_column(
      "Module"
    ),
  file.path(
    OUTDIR,
    "tables",
    "correlations",
    "module_trait_correlations_final.csv"
  )
)

safe_write_csv(
  as.data.frame(p_mat) %>%
    tibble::rownames_to_column(
      "Module"
    ),
  file.path(
    OUTDIR,
    "tables",
    "correlations",
    "module_trait_pvalues_final.csv"
  )
)

safe_write_csv(
  as.data.frame(fdr_mat) %>%
    tibble::rownames_to_column(
      "Module"
    ),
  file.path(
    OUTDIR,
    "tables",
    "correlations",
    "module_trait_fdr_final.csv"
  )
)

safe_write_csv(
  as.data.frame(annot_mat) %>%
    tibble::rownames_to_column(
      "Module"
    ),
  file.path(
    OUTDIR,
    "tables",
    "correlations",
    "module_trait_annotations_final.csv"
  )
)

###############################################################################
# 8) MODULE-TRAIT HEATMAPS
###############################################################################

col_order <- intersect(
  TRAIT_ORDER,
  colnames(rho_mat)
)

rho_mat_plot <- rho_mat[
  modules_use,
  col_order,
  drop = FALSE
]

annot_mat_plot <- annot_mat[
  modules_use,
  col_order,
  drop = FALSE
]

colnames(rho_mat_plot) <- dplyr::recode(
  colnames(rho_mat_plot),
  !!!TRAIT_LABELS,
  .default = colnames(rho_mat_plot)
)

colnames(annot_mat_plot) <- colnames(
  rho_mat_plot
)

annotation_row <- data.frame(
  Module = factor(
    rownames(rho_mat_plot),
    levels = rownames(rho_mat_plot)
  )
)

rownames(annotation_row) <- rownames(
  rho_mat_plot
)

annotation_colors <- list(
  Module = get_module_colors(
    rownames(rho_mat_plot)
  )
)

pdf(
  file.path(
    OUTDIR,
    "figures",
    "module_trait",
    "module_trait_heatmap_fdr.pdf"
  ),
  width = HEATMAP_WIDTH,
  height = HEATMAP_HEIGHT
)

pheatmap::pheatmap(
  rho_mat_plot,
  color = grDevices::colorRampPalette(
    c(
      "#4682B4",
      "white",
      "#F46D43"
    )
  )(101),
  breaks = seq(
    -1,
    1,
    length.out = 102
  ),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  display_numbers = annot_mat_plot,
  number_color = "black",
  fontsize = 10,
  fontsize_number = 12,
  border_color = "grey90",
  annotation_row = annotation_row,
  annotation_colors = annotation_colors,
  main = paste0(
    "WGCNA module-trait associations\n",
    "Spearman rho; stars indicate global BH-FDR significance"
  )
)

dev.off()

png(
  file.path(
    OUTDIR,
    "figures",
    "module_trait",
    "module_trait_heatmap_fdr.png"
  ),
  width = HEATMAP_WIDTH,
  height = HEATMAP_HEIGHT,
  units = "in",
  res = DPI
)

pheatmap::pheatmap(
  rho_mat_plot,
  color = grDevices::colorRampPalette(
    c(
      "#4682B4",
      "white",
      "#F46D43"
    )
  )(101),
  breaks = seq(
    -1,
    1,
    length.out = 102
  ),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  display_numbers = annot_mat_plot,
  number_color = "black",
  fontsize = 10,
  fontsize_number = 12,
  border_color = "grey90",
  annotation_row = annotation_row,
  annotation_colors = annotation_colors,
  main = paste0(
    "WGCNA module-trait associations\n",
    "Spearman rho; stars indicate global BH-FDR significance"
  )
)

dev.off()

heatmap_long <- module_trait_long %>%
  dplyr::mutate(
    Module = factor(
      Module,
      levels = rev(modules_use)
    ),
    Trait = factor(
      Trait,
      levels = col_order
    ),
    Trait_label = factor(
      dplyr::recode(
        as.character(Trait),
        !!!TRAIT_LABELS,
        .default = as.character(Trait)
      ),
      levels = dplyr::recode(
        col_order,
        !!!TRAIT_LABELS,
        .default = col_order
      )
    ),
    cell_label = dplyr::case_when(
      is.na(rho) ~ "",
      stars == "" ~ sprintf("%.2f", rho),
      TRUE ~ paste0(
        sprintf("%.2f", rho),
        "\n",
        stars
      )
    )
  )

p_heatmap_manuscript <- ggplot(
  heatmap_long,
  aes(
    x = Trait_label,
    y = Module,
    fill = rho
  )
) +
  geom_tile(
    colour = "white",
    linewidth = 0.4
  ) +
  geom_text(
    aes(label = cell_label),
    size = 3
  ) +
  scale_fill_gradient2(
    low = "#4682B4",
    mid = "white",
    high = "#F46D43",
    midpoint = 0,
    limits = c(-1, 1),
    oob = scales::squish,
    name = "Spearman\nrho"
  ) +
  labs(
    title = "WGCNA module-trait associations",
    subtitle = paste0(
      "Outcome-independent gene-collapsed network; ",
      "stars indicate global BH-FDR"
    ),
    x = NULL,
    y = NULL
  ) +
  theme_bw(
    base_size = 11
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
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
    "module_trait",
    "module_trait_heatmap_manuscript_style.pdf"
  ),
  p_heatmap_manuscript,
  width = HEATMAP_WIDTH,
  height = HEATMAP_HEIGHT
)

ggsave(
  file.path(
    OUTDIR,
    "figures",
    "module_trait",
    "module_trait_heatmap_manuscript_style.png"
  ),
  p_heatmap_manuscript,
  width = HEATMAP_WIDTH,
  height = HEATMAP_HEIGHT,
  dpi = DPI
)

###############################################################################
# 9) FOCALIZED MODULE-TRAIT PLOTS
###############################################################################

scatter_candidates <- module_trait_long %>%
  dplyr::filter(
    Trait %in%
      SCATTER_PRIORITY_TRAITS,
    !is.na(rho)
  ) %>%
  dplyr::group_by(Module) %>%
  dplyr::arrange(
    FDR,
    dplyr::desc(abs(rho)),
    .by_group = TRUE
  ) %>%
  dplyr::slice_head(
    n = TOP_SCATTERS_PER_MODULE
  ) %>%
  dplyr::ungroup()

safe_write_csv(
  scatter_candidates,
  file.path(
    OUTDIR,
    "tables",
    "scatter_candidates_top_module_trait.csv"
  )
)

for (ii in seq_len(
  nrow(scatter_candidates)
)) {
  row_i <- scatter_candidates[
    ii,
    ,
    drop = FALSE
  ]

  mod <- row_i$Module[[1]]
  trait <- row_i$Trait[[1]]

  trait_lab <- dplyr::recode(
    trait,
    !!!TRAIT_LABELS,
    .default = trait
  )

  if (
    !mod %in% names(analysis_df) ||
    !trait %in% names(analysis_df)
  ) {
    next
  }

  df_plot <- analysis_df %>%
    dplyr::select(
      SampleId,
      dplyr::all_of(c(
        mod,
        trait
      ))
    ) %>%
    tidyr::drop_na()

  if (nrow(df_plot) < MIN_N_FOR_CORR) {
    next
  }

  module_col <- unname(
    get_module_colors(mod)
  )

  if (trait == "SampleGroup_bin") {
    df_plot <- df_plot %>%
      dplyr::mutate(
        Diagnosis = factor(
          .data[[trait]],
          levels = c(0, 1),
          labels = c("CN", "AD")
        )
      )

    p <- ggplot(
      df_plot,
      aes(
        x = Diagnosis,
        y = .data[[mod]]
      )
    ) +
      geom_boxplot(
        outlier.shape = NA,
        fill = module_col,
        alpha = 0.45
      ) +
      geom_jitter(
        width = 0.15,
        alpha = 0.35,
        size = 1
      ) +
      labs(
        title = paste0(
          mod,
          " module vs ",
          trait_lab
        ),
        subtitle = paste0(
          "Spearman rho = ",
          sprintf("%.2f", row_i$rho[[1]]),
          "; FDR = ",
          format(
            row_i$FDR[[1]],
            scientific = TRUE,
            digits = 2
          )
        ),
        x = trait_lab,
        y = paste0(
          mod,
          " module eigengene"
        )
      ) +
      theme_bw(
        base_size = 12
      )
  } else {
    p <- ggplot(
      df_plot,
      aes(
        x = .data[[trait]],
        y = .data[[mod]]
      )
    ) +
      geom_point(
        alpha = 0.55,
        colour = module_col
      ) +
      geom_smooth(
        method = "lm",
        se = TRUE,
        colour = "black",
        linewidth = 0.7
      ) +
      labs(
        title = paste0(
          mod,
          " module vs ",
          trait_lab
        ),
        subtitle = paste0(
          "Spearman rho = ",
          sprintf("%.2f", row_i$rho[[1]]),
          "; FDR = ",
          format(
            row_i$FDR[[1]],
            scientific = TRUE,
            digits = 2
          )
        ),
        x = trait_lab,
        y = paste0(
          mod,
          " module eigengene"
        )
      ) +
      theme_bw(
        base_size = 12
      )
  }

  file_tag <- paste0(
    "scatter_",
    mod,
    "_",
    trait
  )

  ggsave(
    file.path(
      OUTDIR,
      "figures",
      "scatter",
      paste0(
        file_tag,
        ".pdf"
      )
    ),
    p,
    width = 5.5,
    height = 4.5
  )

  ggsave(
    file.path(
      OUTDIR,
      "figures",
      "scatter",
      paste0(
        file_tag,
        ".png"
      )
    ),
    p,
    width = 5.5,
    height = 4.5,
    dpi = DPI
  )
}

###############################################################################
# 10) ADJUSTED CLINICAL AND BIOMARKER MODELS
###############################################################################

adjusted_module_models <- purrr::map_dfr(
  modules_use,
  function(mod) {
    purrr::map_dfr(
      REGRESSION_OUTCOMES,
      function(outcome) {
        if (!outcome %in% names(analysis_df)) {
          return(NULL)
        }

        safe_continuous_module_model(
          df = analysis_df,
          outcome = outcome,
          module = mod,
          covars = CLINICAL_MODEL_COVARS
        )
      }
    )
  }
)

if (nrow(adjusted_module_models) > 0) {
  adjusted_module_models <- adjusted_module_models %>%
    dplyr::mutate(
      FDR = p.adjust(
        p.value,
        method = "BH"
      ),
      Outcome_label = dplyr::recode(
        Outcome,
        !!!TRAIT_LABELS,
        .default = Outcome
      ),
      Biological_label = get_biological_label(
        Module
      )
    ) %>%
    dplyr::arrange(
      FDR,
      p.value
    )
}

safe_write_csv(
  adjusted_module_models,
  file.path(
    OUTDIR,
    "tables",
    "regression",
    "adjusted_module_models.csv"
  )
)

diagnosis_module_models <- purrr::map_dfr(
  modules_use,
  function(mod) {
    safe_diagnosis_module_model(
      df = analysis_df,
      module = mod
    )
  }
) %>%
  dplyr::mutate(
    FDR = p.adjust(
      p.value,
      method = "BH"
    ),
    Biological_label = get_biological_label(
      Module
    )
  ) %>%
  dplyr::arrange(
    FDR,
    p.value
  )

safe_write_csv(
  diagnosis_module_models,
  file.path(
    OUTDIR,
    "tables",
    "regression",
    "adjusted_diagnosis_module_models.csv"
  )
)

###############################################################################
# 11) COUNTRY AS A CATEGORICAL RECRUITMENT-CONTEXT FACTOR
###############################################################################

country_kruskal <- purrr::map_dfr(
  modules_use,
  function(m) {
    safe_kruskal(
      analysis_df,
      outcome = m,
      factor_var = "Country"
    ) %>%
      dplyr::mutate(
        Module = m,
        factor = "Country",
        .before = 1
      )
  }
) %>%
  dplyr::mutate(
    FDR = p.adjust(
      p_value,
      method = "BH"
    ),
    Biological_label = get_biological_label(
      Module
    )
  ) %>%
  dplyr::arrange(
    FDR,
    p_value
  )

country_lm <- purrr::map_dfr(
  modules_use,
  function(m) {
    safe_adjusted_factor_lm(
      analysis_df,
      outcome = m,
      factor_var = "Country",
      covars = CONTEXT_MODEL_COVARS
    ) %>%
      dplyr::mutate(
        Module = m,
        factor = "Country",
        .before = 1
      )
  }
) %>%
  dplyr::mutate(
    FDR_factor = p.adjust(
      p_value_factor,
      method = "BH"
    ),
    Biological_label = get_biological_label(
      Module
    )
  ) %>%
  dplyr::arrange(
    FDR_factor,
    p_value_factor
  )

country_descriptives <- make_factor_descriptives(
  analysis_df,
  modules = modules_use,
  factor_var = "Country"
) %>%
  dplyr::rename(
    Country = Level
  ) %>%
  dplyr::mutate(
    Biological_label = get_biological_label(
      Module
    )
  )

safe_write_csv(
  country_kruskal,
  file.path(
    OUTDIR,
    "tables",
    "context",
    "country_categorical_kruskal_by_module.csv"
  )
)

safe_write_csv(
  country_lm,
  file.path(
    OUTDIR,
    "tables",
    "context",
    "country_categorical_adjusted_lm_by_module.csv"
  )
)

safe_write_csv(
  country_descriptives,
  file.path(
    OUTDIR,
    "tables",
    "context",
    "country_module_descriptives.csv"
  )
)

p_country <- ggplot(
  country_descriptives,
  aes(
    x = Country,
    y = median_eigengene,
    group = Module,
    colour = Module
  )
) +
  geom_point(
    size = 2
  ) +
  geom_line(
    alpha = 0.65
  ) +
  facet_wrap(
    ~ Module,
    scales = "free_y"
  ) +
  scale_colour_manual(
    values = get_module_colors(
      modules_use
    )
  ) +
  labs(
    title = "Module eigengene distributions across countries",
    subtitle = paste0(
      "Descriptive recruitment-context comparison; ",
      "not evidence of country-specific biology"
    ),
    x = "Country",
    y = "Median module eigengene"
  ) +
  theme_bw(
    base_size = 11
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

ggsave(
  file.path(
    OUTDIR,
    "figures",
    "context",
    "country_module_median_eigengenes.pdf"
  ),
  p_country,
  width = 10,
  height = 7
)

ggsave(
  file.path(
    OUTDIR,
    "figures",
    "context",
    "country_module_median_eigengenes.png"
  ),
  p_country,
  width = 10,
  height = 7,
  dpi = DPI
)

###############################################################################
# 12) OPTIONAL SITE AS A CATEGORICAL RECRUITMENT-CONTEXT FACTOR
###############################################################################

site_kruskal <- tibble::tibble()
site_lm <- tibble::tibble()
site_descriptives <- tibble::tibble()

if (
  !is.na(site_var) &&
  site_var %in% names(analysis_df)
) {
  site_kruskal <- purrr::map_dfr(
    modules_use,
    function(m) {
      safe_kruskal(
        analysis_df,
        outcome = m,
        factor_var = site_var
      ) %>%
        dplyr::mutate(
          Module = m,
          factor = site_var,
          .before = 1
        )
    }
  ) %>%
    dplyr::mutate(
      FDR = p.adjust(
        p_value,
        method = "BH"
      ),
      Biological_label = get_biological_label(
        Module
      )
    ) %>%
    dplyr::arrange(
      FDR,
      p_value
    )

  site_lm <- purrr::map_dfr(
    modules_use,
    function(m) {
      safe_adjusted_factor_lm(
        analysis_df,
        outcome = m,
        factor_var = site_var,
        covars = CONTEXT_MODEL_COVARS
      ) %>%
        dplyr::mutate(
          Module = m,
          factor = site_var,
          .before = 1
        )
    }
  ) %>%
    dplyr::mutate(
      FDR_factor = p.adjust(
        p_value_factor,
        method = "BH"
      ),
      Biological_label = get_biological_label(
        Module
      )
    ) %>%
    dplyr::arrange(
      FDR_factor,
      p_value_factor
    )

  site_descriptives <- make_factor_descriptives(
    analysis_df,
    modules = modules_use,
    factor_var = site_var
  ) %>%
    dplyr::rename(
      Site = Level
    ) %>%
    dplyr::mutate(
      Biological_label = get_biological_label(
        Module
      )
    )

  safe_write_csv(
    site_kruskal,
    file.path(
      OUTDIR,
      "tables",
      "context",
      "site_categorical_kruskal_by_module.csv"
    )
  )

  safe_write_csv(
    site_lm,
    file.path(
      OUTDIR,
      "tables",
      "context",
      "site_categorical_adjusted_lm_by_module.csv"
    )
  )

  safe_write_csv(
    site_descriptives,
    file.path(
      OUTDIR,
      "tables",
      "context",
      "site_module_descriptives.csv"
    )
  )
} else {
  safe_write_csv(
    tibble::tibble(
      status = "SITE_VARIABLE_NOT_DETECTED"
    ),
    file.path(
      OUTDIR,
      "tables",
      "context",
      "site_analysis_status.csv"
    )
  )
}

###############################################################################
# 13) CONTEXT EFFECT SUMMARY
###############################################################################

context_effect_summary <- module_reference %>%
  dplyr::left_join(
    country_kruskal %>%
      dplyr::select(
        Module,
        country_kruskal_p = p_value,
        country_kruskal_FDR = FDR
      ),
    by = "Module"
  ) %>%
  dplyr::left_join(
    country_lm %>%
      dplyr::select(
        Module,
        country_adjusted_p = p_value_factor,
        country_adjusted_FDR = FDR_factor,
        country_partial_R2 = partial_R2_factor,
        country_full_model_adj_R2 =
          full_model_adj_R2,
        country_model_status = model_status
      ),
    by = "Module"
  )

if (
  nrow(site_lm) > 0 &&
  nrow(site_kruskal) > 0
) {
  context_effect_summary <- context_effect_summary %>%
    dplyr::left_join(
      site_kruskal %>%
        dplyr::select(
          Module,
          site_kruskal_p = p_value,
          site_kruskal_FDR = FDR
        ),
      by = "Module"
    ) %>%
    dplyr::left_join(
      site_lm %>%
        dplyr::select(
          Module,
          site_adjusted_p = p_value_factor,
          site_adjusted_FDR = FDR_factor,
          site_partial_R2 = partial_R2_factor,
          site_full_model_adj_R2 =
            full_model_adj_R2,
          site_model_status = model_status
        ),
      by = "Module"
    )
}

safe_write_csv(
  context_effect_summary,
  file.path(
    OUTDIR,
    "tables",
    "context",
    "module_recruitment_context_effect_summary.csv"
  )
)

###############################################################################
# 14) INTEGRATED MODULE PRIORITIZATION
###############################################################################

main_trait_summary <- module_trait_long %>%
  dplyr::filter(
    Trait %in%
      MAIN_TRAITS_FOR_INTEGRATION
  ) %>%
  dplyr::group_by(Module) %>%
  dplyr::summarise(
    n_main_traits_tested = dplyr::n(),
    n_trait_fdr_0_05 = sum(
      !is.na(FDR) &
        FDR < 0.05
    ),
    min_trait_FDR = safe_min(
      FDR
    ),
    max_abs_trait_rho = safe_max(
      abs(rho)
    ),
    mean_abs_trait_rho = safe_mean(
      abs(rho)
    ),
    .groups = "drop"
  )

adjusted_clinical_summary <- adjusted_module_models %>%
  dplyr::group_by(Module) %>%
  dplyr::summarise(
    n_adjusted_outcomes_tested = sum(
      !is.na(p.value)
    ),
    n_adjusted_outcomes_fdr_0_05 = sum(
      !is.na(FDR) &
        FDR < 0.05
    ),
    min_adjusted_outcome_FDR = safe_min(
      FDR
    ),
    max_abs_adjusted_beta = safe_max(
      abs(estimate)
    ),
    .groups = "drop"
  )

diagnosis_summary <- diagnosis_module_models %>%
  dplyr::select(
    Module,
    adjusted_diagnosis_beta =
      estimate_diagnosis,
    adjusted_diagnosis_p =
      p.value,
    adjusted_diagnosis_FDR =
      FDR,
    adjusted_diagnosis_N = N
  )

# Begin with Script 12 integrated biology table and add clinical/context domains.
final_prioritization <- integrated_biology %>%
  dplyr::mutate(
    Module = as.character(Module)
  ) %>%
  dplyr::right_join(
    module_reference,
    by = "Module"
  ) %>%
  dplyr::left_join(
    main_trait_summary,
    by = "Module"
  ) %>%
  dplyr::left_join(
    adjusted_clinical_summary,
    by = "Module"
  ) %>%
  dplyr::left_join(
    diagnosis_summary,
    by = "Module"
  ) %>%
  dplyr::left_join(
    context_effect_summary %>%
      dplyr::select(
        -Module_color,
        -Biological_label
      ),
    by = "Module"
  )

# Transparent descriptive prioritization:
# each domain is converted to a percentile rank and all available domain ranks
# receive equal weight. This score is descriptive, not inferential.
final_prioritization <- final_prioritization %>%
  dplyr::mutate(
    dep_domain_rank = percentile_rank_safe(
      primary_DEP_overlap_prop
    ),
    clinical_domain_rank = percentile_rank_safe(
      max_abs_trait_rho
    ),
    adjusted_clinical_domain_rank =
      percentile_rank_safe(
        n_adjusted_outcomes_fdr_0_05
      ),
    enrichment_domain_rank =
      percentile_rank_safe(
        log1p(
          total_significant_enrichment_terms
        )
      ),
    cohesion_domain_rank =
      percentile_rank_safe(
        median_abs_kME
      ),
    priority_score_descriptive = rowMeans(
      cbind(
        dep_domain_rank,
        clinical_domain_rank,
        adjusted_clinical_domain_rank,
        enrichment_domain_rank,
        cohesion_domain_rank
      ),
      na.rm = TRUE
    )
  ) %>%
  dplyr::arrange(
    dplyr::desc(
      priority_score_descriptive
    ),
    dplyr::desc(
      primary_DEP_overlap_prop
    ),
    dplyr::desc(
      max_abs_trait_rho
    )
  ) %>%
  dplyr::mutate(
    priority_rank = dplyr::row_number(),
    prioritization_note = paste(
      "Descriptive integration of equal-weight percentile ranks for",
      "primary DEP burden, unadjusted clinical association strength,",
      "adjusted clinical significance count, enrichment and module cohesion.",
      "Not a hypothesis test."
    )
  )

safe_write_csv(
  main_trait_summary,
  file.path(
    OUTDIR,
    "tables",
    "prioritization",
    "module_main_trait_summary.csv"
  )
)

safe_write_csv(
  adjusted_clinical_summary,
  file.path(
    OUTDIR,
    "tables",
    "prioritization",
    "module_adjusted_clinical_summary.csv"
  )
)

safe_write_csv(
  final_prioritization,
  file.path(
    OUTDIR,
    "tables",
    "prioritization",
    "final_module_prioritization_table.csv"
  )
)

p_prio <- ggplot(
  final_prioritization,
  aes(
    x = max_abs_trait_rho,
    y = primary_DEP_overlap_prop,
    size = total_significant_enrichment_terms,
    colour = Module,
    label = Module
  )
) +
  geom_point(
    alpha = 0.85
  ) +
  ggrepel::geom_text_repel(
    show.legend = FALSE,
    max.overlaps = Inf
  ) +
  scale_colour_manual(
    values = get_module_colors(
      final_prioritization$Module
    )
  ) +
  scale_size_continuous(
    range = c(3, 12)
  ) +
  labs(
    title = "Integrated WGCNA module prioritization",
    subtitle = paste0(
      "DEP burden, clinical association strength and ",
      "functional enrichment"
    ),
    x = "Maximum absolute module-trait Spearman rho",
    y = "Primary DEP overlap proportion",
    size = "Significant\nenrichment terms",
    colour = "Module"
  ) +
  theme_bw(
    base_size = 12
  ) +
  theme(
    plot.title = element_text(
      face = "bold"
    )
  )

ggsave(
  file.path(
    OUTDIR,
    "figures",
    "prioritization",
    "module_prioritization_scatter.pdf"
  ),
  p_prio,
  width = 7,
  height = 5.5
)

ggsave(
  file.path(
    OUTDIR,
    "figures",
    "prioritization",
    "module_prioritization_scatter.png"
  ),
  p_prio,
  width = 7,
  height = 5.5,
  dpi = DPI
)

###############################################################################
# 15) MANUSCRIPT-FACING METHODS AND INTERPRETATION
###############################################################################

methods_wording <- tibble::tribble(
  ~section, ~text,

  "Methods - module-trait correlations",
  paste(
    "Module eigengenes were correlated with clinical, cognitive, demographic,",
    "genetic and plasma biomarker variables using pairwise-complete Spearman",
    "correlations. False-discovery-rate correction was applied globally across",
    "the complete module-trait matrix. Country was not included as a numerically",
    "encoded trait because it is a categorical recruitment-context variable."
  ),

  "Methods - adjusted clinical models",
  paste(
    "For selected continuous clinical and plasma biomarker outcomes, adjusted",
    "linear models were fitted with the outcome as the dependent variable and",
    "the module eigengene as the exposure of interest, accounting for age, sex,",
    "education and country as a categorical factor."
  ),

  "Methods - adjusted diagnosis models",
  paste(
    "Diagnosis-related module differences were additionally evaluated using",
    "linear models with each module eigengene as the dependent variable and",
    "clinical diagnosis as the predictor of interest, adjusting for age, sex,",
    "education and country."
  ),

  "Methods - recruitment context",
  paste(
    "Recruitment-context effects were evaluated separately by treating country",
    "and recruitment site, when available, as categorical factors. For each",
    "module eigengene, Kruskal-Wallis tests compared distributions across",
    "contexts, and nested linear models quantified the additional variance",
    "associated with the categorical context after accounting for diagnosis,",
    "age, sex and education."
  ),

  "Interpretation",
  paste(
    "Country- and site-level analyses were interpreted as internal",
    "recruitment-context sensitivity analyses rather than as evidence of",
    "country-specific biology or external validation."
  ),

  "Module naming",
  paste(
    "Biological module labels were assigned provisionally from the combined",
    "pattern of module enrichment, hub composition, differential-protein burden",
    "and clinical association. WGCNA color labels were retained as stable",
    "internal identifiers."
  )
)

safe_write_csv(
  methods_wording,
  file.path(
    OUTDIR,
    "tables",
    "manuscript_methods_and_interpretation_wording.csv"
  )
)

###############################################################################
# 16) INTEGRATED EXCEL WORKBOOK
###############################################################################

workbook_tables <- list(
  Analysis_summary = input_audit,
  Module_reference = module_reference,
  Module_trait_long = module_trait_long,
  Module_trait_rho = as.data.frame(rho_mat) %>%
    tibble::rownames_to_column("Module"),
  Module_trait_FDR = as.data.frame(fdr_mat) %>%
    tibble::rownames_to_column("Module"),
  Adjusted_clinical = adjusted_module_models,
  Adjusted_diagnosis = diagnosis_module_models,
  Country_Kruskal = country_kruskal,
  Country_Adjusted_LM = country_lm,
  Country_Descriptives = country_descriptives,
  Context_summary = context_effect_summary,
  Main_trait_summary = main_trait_summary,
  Final_prioritization = final_prioritization,
  Methods_wording = methods_wording
)

if (nrow(site_kruskal) > 0) {
  workbook_tables$Site_Kruskal <- site_kruskal
}

if (nrow(site_lm) > 0) {
  workbook_tables$Site_Adjusted_LM <- site_lm
}

if (nrow(site_descriptives) > 0) {
  workbook_tables$Site_Descriptives <- site_descriptives
}

openxlsx::write.xlsx(
  workbook_tables,
  file = file.path(
    OUTDIR,
    "WGCNA_Module_Trait_Context_Integration.xlsx"
  ),
  overwrite = TRUE
)

###############################################################################
# 17) FINAL SUMMARY, MANIFEST AND WORKSPACE
###############################################################################

script13_summary <- tibble::tibble(
  metric = c(
    "base_dir",
    "outdir",
    "wgcna_input_level",
    "n_samples",
    "n_modules",
    "modules",
    "n_traits",
    "traits",
    "n_module_trait_tests",
    "n_module_trait_fdr_0_05",
    "n_adjusted_clinical_models",
    "n_adjusted_clinical_fdr_0_05",
    "n_adjusted_diagnosis_models",
    "n_adjusted_diagnosis_fdr_0_05",
    "n_country_kruskal_fdr_0_05",
    "n_country_adjusted_fdr_0_05",
    "site_variable_detected",
    "n_site_kruskal_fdr_0_05",
    "n_site_adjusted_fdr_0_05",
    "Country_numeric_removed",
    "top_priority_module"
  ),
  value = c(
    BASE_DIR,
    OUTDIR,
    "GENE-COLLAPSED, outcome-independent SOMAmer selection",
    as.character(nrow(analysis_df)),
    as.character(length(modules_use)),
    paste(
      modules_use,
      collapse = ", "
    ),
    as.character(length(traits_available)),
    paste(
      traits_available,
      collapse = ", "
    ),
    as.character(nrow(module_trait_long)),
    as.character(sum(
      module_trait_long$FDR < 0.05,
      na.rm = TRUE
    )),
    as.character(nrow(adjusted_module_models)),
    as.character(sum(
      adjusted_module_models$FDR < 0.05,
      na.rm = TRUE
    )),
    as.character(nrow(diagnosis_module_models)),
    as.character(sum(
      diagnosis_module_models$FDR < 0.05,
      na.rm = TRUE
    )),
    as.character(sum(
      country_kruskal$FDR < 0.05,
      na.rm = TRUE
    )),
    as.character(sum(
      country_lm$FDR_factor < 0.05,
      na.rm = TRUE
    )),
    ifelse(
      is.na(site_var),
      "No",
      site_var
    ),
    as.character(if (
      nrow(site_kruskal) > 0
    ) {
      sum(
        site_kruskal$FDR < 0.05,
        na.rm = TRUE
      )
    } else {
      0
    }),
    as.character(if (
      nrow(site_lm) > 0
    ) {
      sum(
        site_lm$FDR_factor < 0.05,
        na.rm = TRUE
      )
    } else {
      0
    }),
    as.character(
      !"Country_numeric" %in%
        names(analysis_df)
    ),
    as.character(
      final_prioritization$Module[[1]]
    )
  )
)

safe_write_csv(
  script13_summary,
  file.path(
    OUTDIR,
    "tables",
    "script13_final_summary.csv"
  )
)

output_manifest <- tibble::tibble(
  output_file = c(
    "tables/script13_input_alignment_audit.csv",
    "tables/module_biological_label_reference.csv",
    "tables/module_trait_input_clean.csv",
    "tables/correlations/module_trait_results_long.csv",
    "tables/correlations/module_trait_correlations_final.csv",
    "tables/correlations/module_trait_pvalues_final.csv",
    "tables/correlations/module_trait_fdr_final.csv",
    "tables/correlations/module_trait_annotations_final.csv",
    "figures/module_trait/module_trait_heatmap_fdr.pdf/png",
    "figures/module_trait/module_trait_heatmap_manuscript_style.pdf/png",
    "tables/regression/adjusted_module_models.csv",
    "tables/regression/adjusted_diagnosis_module_models.csv",
    "tables/context/country_categorical_kruskal_by_module.csv",
    "tables/context/country_categorical_adjusted_lm_by_module.csv",
    "tables/context/country_module_descriptives.csv",
    "tables/context/site_categorical_*.csv",
    "tables/context/module_recruitment_context_effect_summary.csv",
    "tables/prioritization/module_main_trait_summary.csv",
    "tables/prioritization/module_adjusted_clinical_summary.csv",
    "tables/prioritization/final_module_prioritization_table.csv",
    "tables/manuscript_methods_and_interpretation_wording.csv",
    "WGCNA_Module_Trait_Context_Integration.xlsx",
    "workspace/script13_module_trait_context_workspace.RData",
    "sessionInfo.txt"
  ),
  description = c(
    "Strict sample, module, trait and Country_numeric audit.",
    "Stable module colors and provisional biological labels.",
    "Aligned module eigengenes and standardized metadata.",
    "Global module-trait Spearman results with BH-FDR.",
    "Wide Spearman correlation matrix.",
    "Wide nominal P-value matrix.",
    "Wide global BH-FDR matrix.",
    "Wide significance-star matrix.",
    "Classic module-trait heatmap.",
    "Compact manuscript-style module-trait heatmap.",
    "Continuous outcome ~ module + age + sex + education + categorical country.",
    "Module ~ diagnosis + age + sex + education + categorical country.",
    "Unadjusted categorical country comparisons.",
    "Adjusted nested categorical country models and partial R2.",
    "Module eigengene descriptives by country.",
    "Optional site-level categorical analyses.",
    "Integrated country/site context-effect summary.",
    "Main module-trait signal summary.",
    "Adjusted clinical signal summary.",
    "Biology, DEP, clinical and context integration.",
    "Methods-ready wording and interpretation boundaries.",
    "Integrated reviewer-facing workbook.",
    "Workspace for Script 14 robustness and figures.",
    "R session information."
  )
)

safe_write_csv(
  output_manifest,
  file.path(
    OUTDIR,
    "script13_output_manifest.csv"
  )
)

save(
  analysis_df,
  modules_use,
  traits_available,
  site_var,
  module_reference,
  module_trait_long,
  rho_mat,
  p_mat,
  fdr_mat,
  annot_mat,
  scatter_candidates,
  adjusted_module_models,
  diagnosis_module_models,
  country_kruskal,
  country_lm,
  country_descriptives,
  site_kruskal,
  site_lm,
  site_descriptives,
  context_effect_summary,
  integrated_biology,
  main_trait_summary,
  adjusted_clinical_summary,
  final_prioritization,
  methods_wording,
  script13_summary,
  output_manifest,
  file = file.path(
    OUTDIR,
    "workspace",
    "script13_module_trait_context_workspace.RData"
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
# 18) FINAL MESSAGE
###############################################################################

cat("\nScript 13 finished successfully.\n")
cat("Main output directory:\n", OUTDIR, "\n")
cat("Input level: GENE-COLLAPSED, outcome-independent SOMAmer selection.\n")
cat("Samples analyzed: ", nrow(analysis_df), "\n", sep = "")
cat("Modules analyzed: ", length(modules_use), "\n", sep = "")
cat("Module labels: ", paste(modules_use, collapse = ", "), "\n", sep = "")
cat("Traits analyzed: ", length(traits_available), "\n", sep = "")
cat("Module-trait tests: ", nrow(module_trait_long), "\n", sep = "")
cat(
  "Module-trait associations at FDR < 0.05: ",
  sum(module_trait_long$FDR < 0.05, na.rm = TRUE),
  "\n",
  sep = ""
)
cat(
  "Adjusted clinical models at FDR < 0.05: ",
  sum(adjusted_module_models$FDR < 0.05, na.rm = TRUE),
  "\n",
  sep = ""
)
cat(
  "Adjusted diagnosis models at FDR < 0.05: ",
  sum(diagnosis_module_models$FDR < 0.05, na.rm = TRUE),
  "\n",
  sep = ""
)
cat(
  "Country categorical adjusted models at FDR < 0.05: ",
  sum(country_lm$FDR_factor < 0.05, na.rm = TRUE),
  "\n",
  sep = ""
)
cat(
  "Site variable detected: ",
  ifelse(is.na(site_var), "No", site_var),
  "\n",
  sep = ""
)
cat("Country_numeric included: FALSE\n")
cat(
  "Top descriptive priority module: ",
  final_prioritization$Module[[1]],
  "\n",
  sep = ""
)
cat("\nKey outputs:\n")
cat("- tables/correlations/module_trait_results_long.csv\n")
cat("- figures/module_trait/module_trait_heatmap_fdr.pdf/png\n")
cat("- tables/regression/adjusted_module_models.csv\n")
cat("- tables/regression/adjusted_diagnosis_module_models.csv\n")
cat("- tables/context/country_categorical_adjusted_lm_by_module.csv\n")
cat("- tables/context/module_recruitment_context_effect_summary.csv\n")
cat("- tables/prioritization/final_module_prioritization_table.csv\n")
cat("- WGCNA_Module_Trait_Context_Integration.xlsx\n")

###############################################################################
# END
###############################################################################


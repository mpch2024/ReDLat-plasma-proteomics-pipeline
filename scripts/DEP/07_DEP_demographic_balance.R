###############################################################################
# ReDLat plasma proteomics — DEP workflow
# 07. Demographic balance by recruitment context
# Requires: Analysis workspace from Script 01
# Produces: Country- and site-level balance summaries
# Data policy: participant-level data and intermediate outputs remain local.
###############################################################################

rm(list = ls())

# -----------------------------------------------------------------------------
# Project setup
# -----------------------------------------------------------------------------
.project_root_env <- Sys.getenv("REDLAT_PROJECT_ROOT", unset = "")
if (nzchar(.project_root_env)) {
  project_root <- normalizePath(.project_root_env, winslash = "/", mustWork = TRUE)
} else if (requireNamespace("here", quietly = TRUE)) {
  project_root <- normalizePath(here::here(), winslash = "/", mustWork = TRUE)
} else {
  stop("Package 'here' is required. Restore the project environment with renv::restore().", call. = FALSE)
}
source(file.path(project_root, "R", "dep_bootstrap.R"), local = FALSE)
DEP_CONFIG <- dep_load_config(project_root)
analysis_root <- DEP_CONFIG$result_root
publication_root <- DEP_CONFIG$publication_root


PROJECT_ROOT <- project_root
OUTDIR <- project_root

required_pkgs <- c(
  "dplyr", "tidyr", "purrr", "tibble", "readr", "stringr",
  "ggplot2", "openxlsx"
)
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
}
invisible(lapply(required_pkgs, library, character.only = TRUE))
options(stringsAsFactors = FALSE)

# ----------------------------- User parameters -----------------------------
# Used only when Sex is numerically coded as 1/2 and there are no text labels.
# ReDLat coding: 1 = female and 2 = male. The script also writes
# logs/sex_coding_audit.csv so the resulting classification can be verified.
SEX_NUMERIC_12_FEMALE_CODE <- 1

# Primary complete-case sample should match the main DEP model covariates.
PRIMARY_COVARIATES_REQUIRED <- c("SampleGroup", "Age", "Sex", "Country", "Education")

# SMD interpretation cutoffs used only for descriptive reviewer flags.
SMD_SMALL <- 0.10
SMD_MODERATE <- 0.25
SMD_LARGE <- 0.50

# Cell-size flags for country/site descriptive interpretation.
VERY_SMALL_CELL_N <- 5
SMALL_CELL_N <- 10

# ------------------------------- Helpers -----------------------------------
first_existing_file <- function(paths) {
  hit <- paths[file.exists(paths)][1]
  if (length(hit) == 0 || is.na(hit)) return(NA_character_)
  hit
}

safe_write_csv <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(x, file)
}

save_xlsx_safe <- function(named_list, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  wb <- openxlsx::createWorkbook()
  used <- character(0)
  for (nm in names(named_list)) {
    sheet <- substr(gsub("[^A-Za-z0-9_]+", "_", nm), 1, 31)
    sheet0 <- sheet
    i <- 1
    while (sheet %in% used) {
      suffix <- paste0("_", i)
      sheet <- paste0(substr(sheet0, 1, 31 - nchar(suffix)), suffix)
      i <- i + 1
    }
    used <- c(used, sheet)
    openxlsx::addWorksheet(wb, sheet)
    openxlsx::writeData(wb, sheet, named_list[[nm]])
  }
  openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
}

pick_col <- function(df, candidates) {
  if (is.null(df) || length(names(df)) == 0) return(NA_character_)
  exact <- candidates[candidates %in% names(df)][1]
  if (!is.na(exact)) return(exact)
  clean <- function(x) tolower(gsub("[^a-z0-9]+", "", x))
  nms_clean <- clean(names(df))
  cand_clean <- clean(candidates)
  idx <- match(cand_clean, nms_clean, nomatch = 0)
  idx <- idx[idx > 0][1]
  if (is.na(idx) || length(idx) == 0) return(NA_character_)
  names(df)[idx]
}

safe_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))

fmt_num <- function(x, digits = 3) {
  x <- suppressWarnings(as.numeric(x))
  if (!is.finite(x)) return("NA")
  format(round(x, digits), nsmall = digits, trim = TRUE)
}

fmt_num1 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (!is.finite(x)) return("NA")
  format(round(x, 1), nsmall = 1, trim = TRUE)
}

finite_max_abs <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  max(abs(x), na.rm = TRUE)
}

# Standardizes sex to a binary variable where 1 = female and 0 = male.
# This is used only for descriptive SMD reporting. The original Sex variable is
# still used for complete-case filtering because that is what enters the model.
standardize_sex_binary <- function(x, female_code_12 = SEX_NUMERIC_12_FEMALE_CODE) {
  x_chr <- trimws(as.character(x))
  x_low <- tolower(x_chr)
  num <- suppressWarnings(as.numeric(x_chr))

  out <- rep(NA_real_, length(x_chr))

  female_text <- x_low %in% c("female", "f", "mujer", "woman", "women", "femenino", "feminino")
  male_text <- x_low %in% c("male", "m", "hombre", "man", "men", "masculino")
  out[female_text] <- 1
  out[male_text] <- 0

  # Infer numeric coding only from still-unassigned entries. This avoids the
  # common 1/2 coding being misread as 0/1 coding.
  numeric_unassigned <- is.na(out) & !is.na(num)
  numeric_values <- sort(unique(num[numeric_unassigned]))

  if (length(numeric_values) > 0) {
    uses_01 <- all(numeric_values %in% c(0, 1)) && any(numeric_values == 0)
    uses_12 <- all(numeric_values %in% c(1, 2)) && any(numeric_values == 2) && !any(numeric_values == 0)

    if (uses_01) {
      out[numeric_unassigned & num %in% c(0, 1)] <- num[numeric_unassigned & num %in% c(0, 1)]
    } else if (uses_12) {
      if (female_code_12 == 2) {
        out[numeric_unassigned & num %in% c(1, 2)] <- ifelse(num[numeric_unassigned & num %in% c(1, 2)] == 2, 1, 0)
      } else if (female_code_12 == 1) {
        out[numeric_unassigned & num %in% c(1, 2)] <- ifelse(num[numeric_unassigned & num %in% c(1, 2)] == 1, 1, 0)
      }
    }
  }

  out
}

smd_continuous <- function(x_cn, x_ad) {
  x_cn <- safe_numeric(x_cn)
  x_ad <- safe_numeric(x_ad)
  n_cn <- sum(!is.na(x_cn))
  n_ad <- sum(!is.na(x_ad))
  if (n_cn < 2 || n_ad < 2) return(NA_real_)

  m_cn <- mean(x_cn, na.rm = TRUE)
  m_ad <- mean(x_ad, na.rm = TRUE)
  sd_cn <- stats::sd(x_cn, na.rm = TRUE)
  sd_ad <- stats::sd(x_ad, na.rm = TRUE)
  pooled <- sqrt(((n_cn - 1) * sd_cn^2 + (n_ad - 1) * sd_ad^2) / (n_cn + n_ad - 2))
  if (!is.finite(pooled) || pooled == 0) return(NA_real_)
  (m_ad - m_cn) / pooled
}

# Standardized mean difference for a binary variable using the average of the
# group-specific Bernoulli variances. This is the usual balance-diagnostics form.
smd_binary <- function(x_cn, x_ad) {
  x_cn <- suppressWarnings(as.numeric(x_cn))
  x_ad <- suppressWarnings(as.numeric(x_ad))
  x_cn <- x_cn[x_cn %in% c(0, 1)]
  x_ad <- x_ad[x_ad %in% c(0, 1)]
  if (length(x_cn) == 0 || length(x_ad) == 0) return(NA_real_)

  p_cn <- mean(x_cn, na.rm = TRUE)
  p_ad <- mean(x_ad, na.rm = TRUE)
  denom <- sqrt((p_cn * (1 - p_cn) + p_ad * (1 - p_ad)) / 2)
  if (!is.finite(denom) || denom == 0) return(NA_real_)
  (p_ad - p_cn) / denom
}

add_balance_flags <- function(tbl) {
  tbl %>%
    dplyr::mutate(
      min_cell_n = pmin(n_CN, n_AD),
      cell_size_flag = dplyr::case_when(
        n_CN == 0 | n_AD == 0 ~ "single_group_only",
        min_cell_n < VERY_SMALL_CELL_N ~ "very_small_cell",
        min_cell_n < SMALL_CELL_N ~ "small_cell",
        TRUE ~ "adequate_cell_size"
      ),
      max_abs_SMD = purrr::pmap_dbl(
        list(age_SMD, education_SMD, sex_SMD),
        ~ finite_max_abs(c(...))
      ),
      # Context-level classification: reflects the largest absolute SMD across
      # age, education and sex within the country/site.
      context_balance_flag = dplyr::case_when(
        n_CN == 0 | n_AD == 0 ~ "single_group_only",
        !is.finite(max_abs_SMD) ~ "not_estimable",
        max_abs_SMD >= SMD_LARGE ~ "large_imbalance",
        max_abs_SMD >= SMD_MODERATE ~ "moderate_imbalance",
        max_abs_SMD >= SMD_SMALL ~ "small_imbalance",
        TRUE ~ "well_balanced"
      ),
      # Backward-compatible alias retained for downstream scripts/workbooks that
      # may still expect a column named balance_flag.
      balance_flag = context_balance_flag
    )
}

# Variable-level classification: applies the thresholds to one SMD only.
# This is intentionally separate from context_balance_flag.
classify_individual_smd <- function(x) {
  x <- abs(suppressWarnings(as.numeric(x)))

  dplyr::case_when(
    !is.finite(x) ~ "not_estimable",
    x >= SMD_LARGE ~ "large_imbalance",
    x >= SMD_MODERATE ~ "moderate_imbalance",
    x >= SMD_SMALL ~ "small_imbalance",
    TRUE ~ "well_balanced"
  )
}

summarize_one_context <- function(df, context_vars) {
  df <- df %>%
    dplyr::filter(SampleGroup %in% c("CN", "AD")) %>%
    dplyr::mutate(SampleGroup = factor(SampleGroup, levels = c("CN", "AD")))

  df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(context_vars))) %>%
    dplyr::group_modify(function(.x, .y) {
      cn <- .x %>% dplyr::filter(SampleGroup == "CN")
      ad <- .x %>% dplyr::filter(SampleGroup == "AD")

      tibble::tibble(
        n_total = nrow(.x),
        n_CN = nrow(cn),
        n_AD = nrow(ad),

        age_CN_mean = mean(cn$Age_num, na.rm = TRUE),
        age_CN_sd = stats::sd(cn$Age_num, na.rm = TRUE),
        age_AD_mean = mean(ad$Age_num, na.rm = TRUE),
        age_AD_sd = stats::sd(ad$Age_num, na.rm = TRUE),
        age_diff_AD_minus_CN = age_AD_mean - age_CN_mean,
        age_SMD = smd_continuous(cn$Age_num, ad$Age_num),

        education_CN_mean = mean(cn$Education_num, na.rm = TRUE),
        education_CN_sd = stats::sd(cn$Education_num, na.rm = TRUE),
        education_AD_mean = mean(ad$Education_num, na.rm = TRUE),
        education_AD_sd = stats::sd(ad$Education_num, na.rm = TRUE),
        education_diff_AD_minus_CN = education_AD_mean - education_CN_mean,
        education_SMD = smd_continuous(cn$Education_num, ad$Education_num),

        female_CN_n = sum(cn$Sex_binary == 1, na.rm = TRUE),
        female_CN_pct = mean(cn$Sex_binary == 1, na.rm = TRUE) * 100,
        female_AD_n = sum(ad$Sex_binary == 1, na.rm = TRUE),
        female_AD_pct = mean(ad$Sex_binary == 1, na.rm = TRUE) * 100,
        sex_SMD = smd_binary(cn$Sex_binary, ad$Sex_binary)
      )
    }) %>%
    dplyr::ungroup() %>%
    add_balance_flags()
}

make_key_imbalance_table <- function(tbl, context_level, context_expr) {
  if (!all(c("age_SMD", "education_SMD", "sex_SMD") %in% names(tbl))) {
    return(tibble::tibble(note = paste("Cannot build key imbalance table for", context_level)))
  }

  dplyr::bind_rows(
    tbl %>%
      dplyr::transmute(
        context_level = context_level,
        context = {{ context_expr }},
        variable = "age",
        n_total, n_CN, n_AD, min_cell_n, cell_size_flag,
        mean_CN = age_CN_mean,
        mean_AD = age_AD_mean,
        difference_AD_minus_CN = age_diff_AD_minus_CN,
        SMD = age_SMD,
        context_balance_flag = balance_flag
      ),
    tbl %>%
      dplyr::transmute(
        context_level = context_level,
        context = {{ context_expr }},
        variable = "education",
        n_total, n_CN, n_AD, min_cell_n, cell_size_flag,
        mean_CN = education_CN_mean,
        mean_AD = education_AD_mean,
        difference_AD_minus_CN = education_diff_AD_minus_CN,
        SMD = education_SMD,
        context_balance_flag = balance_flag
      ),
    tbl %>%
      dplyr::transmute(
        context_level = context_level,
        context = {{ context_expr }},
        variable = "female_sex",
        n_total, n_CN, n_AD, min_cell_n, cell_size_flag,
        mean_CN = female_CN_pct,
        mean_AD = female_AD_pct,
        difference_AD_minus_CN = female_AD_pct - female_CN_pct,
        SMD = sex_SMD,
        context_balance_flag = balance_flag
      )
  ) %>%
    dplyr::mutate(
      abs_SMD = abs(SMD),
      variable_balance_flag = classify_individual_smd(SMD)
    ) %>%
    dplyr::arrange(dplyr::desc(abs_SMD))
}

# ------------------------------- Load data ---------------------------------
workspace_file <- first_existing_file(c(
  file.path(analysis_root, "workspace", "analysis_workspace.RData")
))
if (is.na(workspace_file)) {
  stop("Could not find analysis_workspace.RData. Run the DEP differential analysis first.", call. = FALSE)
}
load(workspace_file)
if (!exists("dep_df")) stop("Object dep_df not found in workspace.", call. = FALSE)

OUT_ROOT <- file.path(analysis_root, "06_robustness", "demographic_balance_reviewer_v3")
TABLE_DIR <- file.path(OUT_ROOT, "tables")
FIG_DIR <- file.path(OUT_ROOT, "figures")
LOG_DIR <- file.path(OUT_ROOT, "logs")
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

# --------------------------- Column detection -------------------------------
group_col <- pick_col(dep_df, c("SampleGroup", "sample_group", "Group", "Diagnosis", "Dx", "diagnosis"))
age_col <- pick_col(dep_df, c("Age", "age", "Edad", "edad"))
sex_col <- pick_col(dep_df, c("Sex", "sex", "Sexo", "sexo", "Gender", "gender"))
country_col <- pick_col(dep_df, c("Country", "country", "Pais", "País", "pais", "Country_site"))
education_col <- pick_col(dep_df, c("Education", "education", "Educ", "educ", "YearsEducation", "years_education", "Escolaridad", "escolaridad"))
site_var <- pick_col(dep_df, c("Site", "site", "Center", "center", "Cohort", "cohort", "RecruitmentSite", "Recruitment_Site", "recruitment_site", "Hospital", "hospital"))

col_audit <- tibble::tibble(
  standard_column = c("SampleGroup", "Age", "Sex", "Country", "Education", "Site_optional"),
  detected_column = c(group_col, age_col, sex_col, country_col, education_col, site_var),
  required_for_primary_balance = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE)
)
safe_write_csv(col_audit, file.path(LOG_DIR, "column_detection_audit.csv"))

required_detected <- c(group_col, age_col, sex_col, country_col, education_col)
if (any(is.na(required_detected))) {
  stop(
    "Missing one or more required columns for demographic balance. See logs/column_detection_audit.csv.",
    call. = FALSE
  )
}

# Standardized analytic dataframe.
dep_std <- dep_df %>%
  dplyr::transmute(
    SampleGroup_raw = .data[[group_col]],
    SampleGroup = stringr::str_trim(as.character(.data[[group_col]])),
    Age_raw = .data[[age_col]],
    Age_num = safe_numeric(.data[[age_col]]),
    Sex_raw = .data[[sex_col]],
    Sex_binary = standardize_sex_binary(.data[[sex_col]]),
    Country_raw = .data[[country_col]],
    Country = stringr::str_trim(as.character(.data[[country_col]])),
    Education_raw = .data[[education_col]],
    Education_num = safe_numeric(.data[[education_col]]),
    Site_detected = if (!is.na(site_var)) stringr::str_trim(as.character(.data[[site_var]])) else NA_character_
  ) %>%
  dplyr::mutate(
    SampleGroup = dplyr::case_when(
      SampleGroup %in% c("CN", "Control", "Controls", "Cognitively normal", "Cognitively Normal") ~ "CN",
      SampleGroup %in% c("AD", "Alzheimer", "Alzheimer disease", "Alzheimer's disease", "clinically diagnosed AD") ~ "AD",
      TRUE ~ SampleGroup
    )
  )

sex_coding_audit <- dep_std %>%
  dplyr::count(Sex_raw, Sex_binary, name = "n") %>%
  dplyr::arrange(dplyr::desc(n))
safe_write_csv(sex_coding_audit, file.path(LOG_DIR, "sex_coding_audit.csv"))

missingness_audit <- dep_std %>%
  dplyr::filter(SampleGroup %in% c("CN", "AD")) %>%
  dplyr::summarise(
    n_all_CN_AD = dplyr::n(),
    missing_Age = sum(is.na(Age_num)),
    missing_Sex_original = sum(is.na(Sex_raw) | trimws(as.character(Sex_raw)) == ""),
    missing_Sex_binary_after_standardization = sum(is.na(Sex_binary)),
    missing_Country = sum(is.na(Country) | Country == ""),
    missing_Education = sum(is.na(Education_num)),
    complete_primary_covariates = sum(
      !is.na(Age_num) &
        !(is.na(Sex_raw) | trimws(as.character(Sex_raw)) == "") &
        !(is.na(Country) | Country == "") &
        !is.na(Education_num)
    )
  )
safe_write_csv(missingness_audit, file.path(LOG_DIR, "primary_covariate_missingness_audit.csv"))

# Main sample used for reviewer-facing balance metrics.
analysis_df_all_CN_AD <- dep_std %>%
  dplyr::filter(SampleGroup %in% c("CN", "AD"))

analysis_df <- analysis_df_all_CN_AD %>%
  dplyr::filter(
    !is.na(Age_num),
    !(is.na(Sex_raw) | trimws(as.character(Sex_raw)) == ""),
    !(is.na(Country) | Country == ""),
    !is.na(Education_num)
  )

sample_tracking <- tibble::tibble(
  sample_set = c("all_CN_AD", "primary_complete_case_Age_Sex_Country_Education"),
  n_total = c(nrow(analysis_df_all_CN_AD), nrow(analysis_df)),
  n_CN = c(sum(analysis_df_all_CN_AD$SampleGroup == "CN"), sum(analysis_df$SampleGroup == "CN")),
  n_AD = c(sum(analysis_df_all_CN_AD$SampleGroup == "AD"), sum(analysis_df$SampleGroup == "AD")),
  note = c(
    "All CN/AD participants available in dep_df before primary covariate complete-case filtering.",
    "Main demographic balance sample, matching the primary DEP covariate structure."
  )
)
safe_write_csv(sample_tracking, file.path(LOG_DIR, "sample_tracking.csv"))

# ----------------------------- Main analyses --------------------------------
overall_tbl <- summarize_one_context(analysis_df %>% dplyr::mutate(All = "Overall"), "All")
country_tbl <- summarize_one_context(analysis_df, "Country") %>%
  dplyr::arrange(dplyr::desc(abs(age_SMD)), dplyr::desc(abs(education_SMD)), dplyr::desc(abs(sex_SMD)))

if (!is.na(site_var)) {
  site_tbl <- summarize_one_context(analysis_df, c("Country", "Site_detected")) %>%
    dplyr::arrange(Country, dplyr::desc(max_abs_SMD), dplyr::desc(abs(age_SMD)), dplyr::desc(abs(education_SMD)))
} else {
  site_tbl <- tibble::tibble(note = "No site/center/cohort variable detected in dep_df.")
}

safe_write_csv(overall_tbl, file.path(TABLE_DIR, "demographic_balance_overall_primary_complete_case.csv"))
safe_write_csv(country_tbl, file.path(TABLE_DIR, "demographic_balance_by_country_primary_complete_case.csv"))
safe_write_csv(site_tbl, file.path(TABLE_DIR, "demographic_balance_by_country_site_primary_complete_case.csv"))

# Long SMD tables.
long_smd_country <- country_tbl %>%
  dplyr::select(
    Country,
    n_total, n_CN, n_AD,
    min_cell_n,
    cell_size_flag,
    age_SMD,
    education_SMD,
    sex_SMD,
    context_balance_flag = balance_flag
  ) %>%
  tidyr::pivot_longer(
    cols = c(age_SMD, education_SMD, sex_SMD),
    names_to = "variable",
    values_to = "SMD"
  ) %>%
  dplyr::mutate(
    abs_SMD = abs(SMD),
    variable_balance_flag = classify_individual_smd(SMD)
  )
safe_write_csv(long_smd_country, file.path(TABLE_DIR, "demographic_balance_country_SMD_long_primary_complete_case.csv"))

country_key_imbalances <- make_key_imbalance_table(country_tbl, "country", Country)
safe_write_csv(country_key_imbalances, file.path(TABLE_DIR, "demographic_balance_country_key_imbalances_ranked.csv"))

if (!is.na(site_var) && all(c("Country", "Site_detected", "age_SMD", "education_SMD", "sex_SMD") %in% names(site_tbl))) {
  site_key_imbalances <- make_key_imbalance_table(site_tbl, "site", paste(Country, Site_detected, sep = " / "))
  safe_write_csv(site_key_imbalances, file.path(TABLE_DIR, "demographic_balance_site_key_imbalances_ranked.csv"))

  long_smd_site <- site_tbl %>%
    dplyr::mutate(context = paste(Country, Site_detected, sep = " / ")) %>%
    dplyr::select(
      context,
      Country,
      Site_detected,
      n_total, n_CN, n_AD,
      min_cell_n,
      cell_size_flag,
      age_SMD,
      education_SMD,
      sex_SMD,
      context_balance_flag = balance_flag
    ) %>%
    tidyr::pivot_longer(
      cols = c(age_SMD, education_SMD, sex_SMD),
      names_to = "variable",
      values_to = "SMD"
    ) %>%
    dplyr::mutate(
      abs_SMD = abs(SMD),
      variable_balance_flag = classify_individual_smd(SMD)
    )
  safe_write_csv(long_smd_site, file.path(TABLE_DIR, "demographic_balance_site_SMD_long_primary_complete_case.csv"))
} else {
  site_key_imbalances <- tibble::tibble(note = "No site-level imbalance table generated because no site variable was detected.")
  long_smd_site <- tibble::tibble(note = "No site-level long SMD table generated because no site variable was detected.")
}

# ------------------------------- Figures -----------------------------------
p_country <- ggplot(long_smd_country, aes(x = reorder(Country, abs_SMD, FUN = max), y = SMD)) +
  geom_hline(yintercept = c(-SMD_SMALL, SMD_SMALL), linetype = "dotted") +
  geom_hline(yintercept = c(-SMD_MODERATE, SMD_MODERATE), linetype = "dashed") +
  geom_hline(yintercept = c(-SMD_LARGE, SMD_LARGE), linetype = "solid") +
  geom_col() +
  coord_flip() +
  facet_wrap(~ variable, scales = "free_x") +
  labs(
    title = "Demographic balance by country",
    subtitle = "Primary complete-case analytic sample; SMD > 0 indicates higher values in AD",
    x = NULL,
    y = "Standardized mean difference"
  ) +
  theme_classic(base_size = 9)

ggsave(file.path(FIG_DIR, "demographic_balance_by_country_SMD_v3.pdf"), p_country, width = 7.5, height = 4.8)
ggsave(file.path(FIG_DIR, "demographic_balance_by_country_SMD_v3.png"), p_country, width = 7.5, height = 4.8, dpi = 300)

if (nrow(long_smd_site) > 0 && "SMD" %in% names(long_smd_site)) {
  p_site <- long_smd_site %>%
    dplyr::group_by(context) %>%
    dplyr::mutate(max_abs_SMD_context = max(abs_SMD, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    ggplot(aes(x = reorder(context, max_abs_SMD_context), y = SMD)) +
    geom_hline(yintercept = c(-SMD_SMALL, SMD_SMALL), linetype = "dotted") +
    geom_hline(yintercept = c(-SMD_MODERATE, SMD_MODERATE), linetype = "dashed") +
    geom_hline(yintercept = c(-SMD_LARGE, SMD_LARGE), linetype = "solid") +
    geom_col() +
    coord_flip() +
    facet_wrap(~ variable, scales = "free_x") +
    labs(
      title = "Demographic balance by recruitment site",
      subtitle = "Primary complete-case analytic sample; SMD > 0 indicates higher values in AD",
      x = NULL,
      y = "Standardized mean difference"
    ) +
    theme_classic(base_size = 8)

  ggsave(file.path(FIG_DIR, "demographic_balance_by_site_SMD_v3.pdf"), p_site, width = 8.5, height = 6.0)
  ggsave(file.path(FIG_DIR, "demographic_balance_by_site_SMD_v3.png"), p_site, width = 8.5, height = 6.0, dpi = 300)
}

# ----------------------- Reviewer/manuscript summaries ----------------------
top_age_country <- country_key_imbalances %>%
  dplyr::filter(variable == "age", is.finite(SMD), n_CN > 0, n_AD > 0) %>%
  dplyr::slice_max(abs_SMD, n = 1, with_ties = FALSE)

top_edu_country <- country_key_imbalances %>%
  dplyr::filter(variable == "education", is.finite(SMD), n_CN > 0, n_AD > 0) %>%
  dplyr::slice_max(abs_SMD, n = 1, with_ties = FALSE)

top_sex_country <- country_key_imbalances %>%
  dplyr::filter(variable == "female_sex", is.finite(SMD), n_CN > 0, n_AD > 0) %>%
  dplyr::slice_max(abs_SMD, n = 1, with_ties = FALSE)

n_small_site_contexts <- if ("cell_size_flag" %in% names(site_tbl)) {
  sum(site_tbl$cell_size_flag %in% c("very_small_cell", "small_cell"), na.rm = TRUE)
} else NA_integer_

site_note <- if (!is.na(site_var)) {
  paste0(
    "Site-level analyses used the detected site variable '", site_var,
    "' and should be interpreted descriptively, especially for contexts with small cells."
  )
} else {
  "No site variable was detected; site-level balance could not be summarized."
}

demographic_reviewer_takehome <- tibble::tibble(
  finding = c(
    "analytic_sample",
    "overall_balance",
    "largest_age_imbalance_country",
    "largest_education_imbalance_country",
    "largest_sex_imbalance_country",
    "site_level_cell_size_note",
    "recommended_interpretation"
  ),
  value = c(
    paste0(
      "Primary complete-case sample n=", nrow(analysis_df),
      " (CN=", sum(analysis_df$SampleGroup == "CN"),
      "; AD=", sum(analysis_df$SampleGroup == "AD"), ")."
    ),
    paste0(
      "Overall age SMD=", fmt_num(overall_tbl$age_SMD[1], 3),
      "; education SMD=", fmt_num(overall_tbl$education_SMD[1], 3),
      "; sex SMD=", fmt_num(overall_tbl$sex_SMD[1], 3), "."
    ),
    if (nrow(top_age_country) > 0) paste0(top_age_country$context[1], " (SMD=", fmt_num(top_age_country$SMD[1], 3), ")") else "not estimable",
    if (nrow(top_edu_country) > 0) paste0(top_edu_country$context[1], " (SMD=", fmt_num(top_edu_country$SMD[1], 3), ")") else "not estimable",
    if (nrow(top_sex_country) > 0) paste0(top_sex_country$context[1], " (SMD=", fmt_num(top_sex_country$SMD[1], 3), ")") else "not estimable",
    paste0(site_note, " Small/very small site cells: ", ifelse(is.na(n_small_site_contexts), "NA", n_small_site_contexts), "."),
    "Age and education imbalance varied across recruitment contexts, supporting covariate adjustment, propensity-score matching and country/site robustness analyses. Site-level SMDs are descriptive and should be interpreted in relation to cell size."
  )
)
safe_write_csv(demographic_reviewer_takehome, file.path(TABLE_DIR, "demographic_balance_reviewer_takehome.csv"))

# Main manuscript-ready Results paragraph.
results_text <- paste0(
  "To characterize demographic balance across recruitment contexts, we quantified standardized mean differences for age, education and sex in the primary complete-case analytic sample used for the differential-abundance models (n = ",
  nrow(analysis_df), "; CN = ", sum(analysis_df$SampleGroup == "CN"), "; clinically diagnosed AD = ", sum(analysis_df$SampleGroup == "AD"), "). ",
  "At the overall cohort level, diagnostic groups showed marked age imbalance (SMD = ", fmt_num(overall_tbl$age_SMD[1], 2),
  ") and smaller education imbalance (SMD = ", fmt_num(overall_tbl$education_SMD[1], 2),
  "), whereas sex distribution was comparatively balanced (SMD = ", fmt_num(overall_tbl$sex_SMD[1], 2), "). ",
  "Country- and site-stratified balance diagnostics showed that these differences were not uniform across recruitment contexts, with the largest country-level age imbalance observed in ",
  if (nrow(top_age_country) > 0) as.character(top_age_country$context[1]) else "not estimable",
  " and the largest country-level education imbalance observed in ",
  if (nrow(top_edu_country) > 0) as.character(top_edu_country$context[1]) else "not estimable",
  ". These findings support the use of covariate adjustment, propensity-score matching and country/site robustness analyses to evaluate the stability of the clinical AD-associated proteomic profile across recruitment contexts."
)

methods_text <- paste0(
  "Demographic balance between cognitively normal and clinically diagnosed AD participants was evaluated using standardized mean differences for age, education and sex. ",
  "Balance diagnostics were computed in the primary complete-case analytic sample defined by availability of diagnosis group, age, sex, country and education. ",
  "Continuous-variable SMDs were calculated using pooled standard deviations, whereas binary sex SMDs were calculated using the average of group-specific Bernoulli variances. ",
  "Diagnostics were summarized overall and after stratification by country and recruitment site when site information was available. Site-level estimates were interpreted descriptively and flagged when one diagnostic group was absent or when the minimum group cell size was small."
)

supplementary_legend_text <- paste0(
  "Supplementary Table XX. Demographic balance by country and recruitment site in the primary complete-case analytic sample. ",
  "Standardized mean differences are shown for age, education and sex. Positive SMD values indicate higher values or proportions in clinically diagnosed AD relative to cognitively normal participants. ",
  "Cell-size flags identify recruitment contexts with single-group-only, very small or small diagnostic cells."
)

writeLines(results_text, file.path(OUT_ROOT, "manuscript_ready_demographic_balance_Results_text.txt"))
writeLines(methods_text, file.path(OUT_ROOT, "manuscript_ready_demographic_balance_Methods_text.txt"))
writeLines(supplementary_legend_text, file.path(OUT_ROOT, "supplementary_table_legend_demographic_balance.txt"))

# ------------------------------- Workbook ----------------------------------
summary_tbl <- tibble::tibble(
  item = c(
    "workspace_file",
    "site_variable_detected",
    "n_all_CN_AD",
    "n_primary_complete_case",
    "n_CN_primary_complete_case",
    "n_AD_primary_complete_case",
    "n_countries_primary_complete_case",
    "max_abs_age_SMD_country",
    "max_abs_education_SMD_country",
    "max_abs_sex_SMD_country",
    "sex_numeric_12_female_code_parameter",
    "output_directory"
  ),
  value = c(
    workspace_file,
    ifelse(is.na(site_var), "none", site_var),
    as.character(nrow(analysis_df_all_CN_AD)),
    as.character(nrow(analysis_df)),
    as.character(sum(analysis_df$SampleGroup == "CN")),
    as.character(sum(analysis_df$SampleGroup == "AD")),
    as.character(dplyr::n_distinct(analysis_df$Country)),
    as.character(max(abs(country_tbl$age_SMD), na.rm = TRUE)),
    as.character(max(abs(country_tbl$education_SMD), na.rm = TRUE)),
    as.character(max(abs(country_tbl$sex_SMD), na.rm = TRUE)),
    as.character(SEX_NUMERIC_12_FEMALE_CODE),
    OUT_ROOT
  )
)
safe_write_csv(summary_tbl, file.path(OUT_ROOT, "demographic_balance_reviewer_summary.csv"))

save_xlsx_safe(
  list(
    summary = summary_tbl,
    sample_tracking = sample_tracking,
    column_detection = col_audit,
    covariate_missingness = missingness_audit,
    sex_coding_audit = sex_coding_audit,
    overall_primary = overall_tbl,
    by_country_primary = country_tbl,
    by_country_site_primary = site_tbl,
    country_SMD_long = long_smd_country,
    site_SMD_long = long_smd_site,
    country_key_imbalances = country_key_imbalances,
    site_key_imbalances = site_key_imbalances,
    reviewer_takehome = demographic_reviewer_takehome
  ),
  file.path(OUT_ROOT, "Demographic_Balance_Country_Site_Reviewer_Audit_v3.xlsx")
)

script_metadata <- tibble::tibble(
  item = c(
    "script",
    "project_root",
    "workspace_file",
    "output_directory",
    "primary_sample_definition",
    "date_time"
  ),
  value = c(
    "06B_demographic_balance_by_country_site_REVIEWER_v4_FIXED_BALANCE_FLAGS.R",
    PROJECT_ROOT,
    workspace_file,
    OUT_ROOT,
    "CN/AD participants complete for SampleGroup, Age, Sex, Country and Education.",
    as.character(Sys.time())
  )
)
safe_write_csv(script_metadata, file.path(LOG_DIR, "script_metadata.csv"))
writeLines(capture.output(utils::sessionInfo()), file.path(LOG_DIR, "sessionInfo.txt"))

message("Demographic balance audit complete.")
message("Primary complete-case sample n = ", nrow(analysis_df),
        " (CN = ", sum(analysis_df$SampleGroup == "CN"),
        "; AD = ", sum(analysis_df$SampleGroup == "AD"), ").")
message("Results text saved to: ", file.path(OUT_ROOT, "manuscript_ready_demographic_balance_Results_text.txt"))
message("Methods text saved to: ", file.path(OUT_ROOT, "manuscript_ready_demographic_balance_Methods_text.txt"))
message("Outputs saved to: ", OUT_ROOT)


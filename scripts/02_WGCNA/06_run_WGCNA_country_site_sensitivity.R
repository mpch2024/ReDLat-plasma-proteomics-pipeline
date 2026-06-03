###############################################################################
# 05_WGCNA_country_site_downsampling_sensitivity_FINAL_v2.R
#
# TITLE
# Country, site and balanced downsampling sensitivity analyses for
# GENE-COLLAPSED WGCNA module-trait associations
#
# PURPOSE
# This script evaluates robustness of WGCNA module-trait associations using:
#   1) full module-trait Spearman correlations
#   2) leave-one-country-out sensitivity
#   3) leave-one-site-out sensitivity if a site/center variable exists
#   4) balanced downsampling by Country and diagnosis
#
# IMPORTANT
# This script DOES NOT rebuild WGCNA.
# It uses module eigengenes from Script 02 and metadata from Script 01.
#
# WGCNA input level:
#   GENE-COLLAPSED, not aptamer-level.
#
# MAIN MODULES FOR MANUSCRIPT VISUALIZATION
#   yellow = clinically integrated / neuronal-ECM module
#   brown  = highest DEP-burden module
#   black  = narrow biomarker-linked module
#
# INPUTS
#   Script 01:
#     results/01_define_wgcna_input_from_DEP/wgcna_sample_metadata.csv
#
#   Script 02:
#     results/02_WGCNA_core_collapsed_genes/eigengenes/module_eigengenes_per_sample.csv
#
#   Script 04:
#     results/04_WGCNA_module_trait_clinical_integration/tables/module_trait_results_long.csv
#     results/04_WGCNA_module_trait_clinical_integration/tables/final_module_prioritization_table.csv
#
# OUTPUTS
#   results/05_WGCNA_country_site_downsampling_sensitivity/
#
# AUTHOR
# Matías Pizarro + ChatGPT support
###############################################################################

rm(list = ls())

###############################################################################
# 1) PACKAGES
###############################################################################

cran_pkgs <- c(
  "dplyr",
  "tidyr",
  "readr",
  "stringr",
  "purrr",
  "tibble",
  "ggplot2",
  "scales"
)

cran_missing <- cran_pkgs[
  !sapply(cran_pkgs, requireNamespace, quietly = TRUE)
]

if (length(cran_missing) > 0) {
  install.packages(cran_missing)
}

invisible(lapply(cran_pkgs, library, character.only = TRUE))

options(stringsAsFactors = FALSE)
options(error = traceback)

###############################################################################
# 2) PATHS
###############################################################################

BASE_DIR <- "C:/Users/mnpiz/Desktop/WGCNA_Workflow_april_V2"

SCRIPT1_DIR <- file.path(
  BASE_DIR,
  "results",
  "01_define_wgcna_input_from_DEP"
)

SCRIPT2_DIR <- file.path(
  BASE_DIR,
  "results",
  "02_WGCNA_core_collapsed_genes"
)

SCRIPT4_DIR <- file.path(
  BASE_DIR,
  "results",
  "04_WGCNA_module_trait_clinical_integration"
)

EIGENGENE_FILE <- file.path(
  SCRIPT2_DIR,
  "eigengenes",
  "module_eigengenes_per_sample.csv"
)

META_FILE <- file.path(
  SCRIPT1_DIR,
  "wgcna_sample_metadata.csv"
)

MODULE_TRAIT_FILE <- file.path(
  SCRIPT4_DIR,
  "tables",
  "module_trait_results_long.csv"
)

PRIORITIZATION_FILE <- file.path(
  SCRIPT4_DIR,
  "tables",
  "final_module_prioritization_table.csv"
)

OUTDIR <- file.path(
  BASE_DIR,
  "results",
  "05_WGCNA_country_site_downsampling_sensitivity"
)

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

SUBDIRS <- c(
  "tables",
  "tables/loco",
  "tables/loso",
  "tables/downsampling",
  "figures",
  "figures/loco",
  "figures/loso",
  "figures/downsampling",
  "figures/core_modules",
  "workspace"
)

invisible(lapply(
  file.path(OUTDIR, SUBDIRS),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

REQUIRED_FILES <- c(
  EIGENGENE_FILE,
  META_FILE
)

missing_files <- REQUIRED_FILES[!file.exists(REQUIRED_FILES)]

if (length(missing_files) > 0) {
  stop(
    "Faltan archivos requeridos:\n",
    paste(missing_files, collapse = "\n"),
    call. = FALSE
  )
}

###############################################################################
# 3) PARAMETERS
###############################################################################

MIN_N_FOR_CORR <- 6
N_DOWNSAMPLING_ITER <- 500
SEED <- 123

EXCLUDE_MODULES <- c("grey", "MEgrey")

CORE_MODULES <- c(
  "yellow",
  "brown",
  "black"
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
  cdr_global      = "CDR global",
  cdr_boxscore    = "CDR-SB",
  mmse_total      = "MMSE",
  udsfaq_total    = "PFAQ",
  NPI             = "NPI-Q",
  Mini_SEA        = "Mini-SEA",
  T_ADLQ          = "T-ADLQ",
  p_tau181        = "p-tau181",
  p_tau217        = "p-tau217",
  NfL             = "NfL",
  ratio_AB42_40   = "Aβ42/40",
  Age             = "Age",
  Sex_bin         = "Sex",
  Education       = "Education",
  APOE4_carrier   = "APOE ε4"
)

SITE_CANDIDATES <- c(
  "Site",
  "site",
  "Center",
  "center",
  "RecruitmentSite",
  "recruitment_site",
  "StudySite",
  "study_site",
  "Cohort",
  "cohort"
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
  grey        = "#9E9E9E"
)

###############################################################################
# 4) HELPERS
###############################################################################

safe_write_csv <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(x, file)
}

clean_text_na <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA", "NaN", "NULL", "null", "N/A", "nan")] <- NA_character_
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
    extra_cols <- grDevices::rainbow(length(unique(missing_mods)), s = 0.45, v = 0.75)
    names(extra_cols) <- unique(missing_mods)
    cols[is.na(cols)] <- extra_cols[missing_mods]
  }

  names(cols) <- modules
  cols
}

standardize_metadata_columns <- function(df) {
  rename_map <- c(
    "Mini-SEA" = "Mini_SEA",
    "Mini.SEA" = "Mini_SEA",
    "T-ADLQ" = "T_ADLQ",
    "T.ADLQ" = "T_ADLQ",
    "p-tau181" = "p_tau181",
    "p.tau181" = "p_tau181",
    "p-tau217" = "p_tau217",
    "p.tau217" = "p_tau217",
    "ratio AB42/40" = "ratio_AB42_40",
    "ratio.AB42.40" = "ratio_AB42_40"
  )

  for (old in names(rename_map)) {
    new <- rename_map[[old]]
    if (old %in% names(df) && !new %in% names(df)) {
      names(df)[names(df) == old] <- new
    }
  }

  df
}

encode_basic_variables <- function(df) {
  df <- df %>%
    dplyr::mutate(SampleId = as.character(SampleId))

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
          Sex_chr %in% c("1", "M", "Male", "male") ~ 0,
          Sex_chr %in% c("2", "F", "Female", "female") ~ 1,
          TRUE ~ safe_numeric(Sex_chr)
        )
      )
  }

  if ("ApoE" %in% names(df) && !"APOE_group" %in% names(df)) {
    df <- df %>%
      dplyr::mutate(
        ApoE = as.character(ApoE),
        APOE_group = dplyr::case_when(
          ApoE %in% c("e2/e4", "e3/e4", "e4/e4") ~ "E4 carrier",
          ApoE %in% c("e2/e2", "e2/e3", "e3/e3") ~ "Non-E4",
          TRUE ~ NA_character_
        )
      )
  }

  if ("APOE_group" %in% names(df) && !"APOE4_carrier" %in% names(df)) {
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
    "Age", "Education", "cdr_global", "cdr_boxscore", "mmse_total",
    "udsfaq_total", "NPI", "Mini_SEA", "T_ADLQ", "p_tau181",
    "p_tau217", "NfL", "ratio_AB42_40", "APOE4_carrier"
  )

  for (cc in intersect(numeric_candidates, names(df))) {
    df[[cc]] <- safe_numeric(df[[cc]])
  }

  df
}

safe_spearman <- function(df, x, y, min_n = 6) {
  d <- df %>%
    dplyr::select(dplyr::all_of(c(x, y))) %>%
    tidyr::drop_na()

  n <- nrow(d)

  if (n < min_n) {
    return(tibble::tibble(rho = NA_real_, p_value = NA_real_, N = n))
  }

  if (length(unique(d[[x]])) <= 1 || length(unique(d[[y]])) <= 1) {
    return(tibble::tibble(rho = NA_real_, p_value = NA_real_, N = n))
  }

  test <- suppressWarnings(
    stats::cor.test(d[[x]], d[[y]], method = "spearman", exact = FALSE)
  )

  tibble::tibble(
    rho = unname(test$estimate),
    p_value = test$p.value,
    N = n
  )
}

compute_module_trait <- function(df, modules, traits, label = "full") {
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
        ~ safe_spearman(df, .x, .y, min_n = MIN_N_FOR_CORR)
      )
    ) %>%
    tidyr::unnest(result) %>%
    dplyr::mutate(Analysis = label)
}

make_delta_summary <- function(full_tbl, sens_tbl) {
  sens_tbl %>%
    dplyr::left_join(
      full_tbl %>%
        dplyr::select(Module, Trait, rho_full = rho, p_full = p_value, N_full = N),
      by = c("Module", "Trait")
    ) %>%
    dplyr::mutate(
      delta_rho = rho - rho_full,
      abs_delta_rho = abs(delta_rho),
      direction_full = sign(rho_full),
      direction_sensitivity = sign(rho),
      same_direction = direction_full == direction_sensitivity
    )
}

summarize_delta_by_module <- function(delta_tbl) {
  delta_tbl %>%
    dplyr::group_by(Module) %>%
    dplyr::summarise(
      mean_abs_delta_rho = mean(abs_delta_rho, na.rm = TRUE),
      median_abs_delta_rho = median(abs_delta_rho, na.rm = TRUE),
      max_abs_delta_rho = max(abs_delta_rho, na.rm = TRUE),
      direction_consistency = mean(same_direction, na.rm = TRUE),
      n_tests = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::arrange(mean_abs_delta_rho)
}

summarize_delta_by_module_trait <- function(delta_tbl) {
  delta_tbl %>%
    dplyr::group_by(Module, Trait) %>%
    dplyr::summarise(
      mean_abs_delta_rho = mean(abs_delta_rho, na.rm = TRUE),
      median_abs_delta_rho = median(abs_delta_rho, na.rm = TRUE),
      max_abs_delta_rho = max(abs_delta_rho, na.rm = TRUE),
      direction_consistency = mean(same_direction, na.rm = TRUE),
      n_sensitivity_runs = dplyr::n(),
      .groups = "drop"
    )
}

###############################################################################
# 5) LOAD AND PREPARE DATA
###############################################################################

eig <- readr::read_csv(EIGENGENE_FILE, show_col_types = FALSE)
meta <- readr::read_csv(META_FILE, show_col_types = FALSE)

if (file.exists(MODULE_TRAIT_FILE)) {
  module_trait_reference <- readr::read_csv(MODULE_TRAIT_FILE, show_col_types = FALSE)
} else {
  module_trait_reference <- NULL
}

if (file.exists(PRIORITIZATION_FILE)) {
  prioritization_reference <- readr::read_csv(PRIORITIZATION_FILE, show_col_types = FALSE)
} else {
  prioritization_reference <- NULL
}

eig <- eig %>%
  dplyr::mutate(SampleId = as.character(SampleId)) %>%
  dplyr::distinct(SampleId, .keep_all = TRUE)

meta <- meta %>%
  standardize_metadata_columns() %>%
  encode_basic_variables() %>%
  dplyr::distinct(SampleId, .keep_all = TRUE)

eig_cols <- setdiff(names(eig), "SampleId")
eig_cols <- eig_cols[!eig_cols %in% EXCLUDE_MODULES]

eig_renamed <- eig
names(eig_renamed)[names(eig_renamed) %in% eig_cols] <- paste0(
  "ME",
  standardize_module_name(eig_cols)
)

module_cols <- setdiff(names(eig_renamed), "SampleId")
module_cols <- module_cols[!module_cols %in% c("MEgrey", "grey")]

analysis_df <- eig_renamed %>%
  dplyr::inner_join(meta, by = "SampleId")

if (nrow(analysis_df) == 0) {
  stop("No hay muestras compartidas entre eigengenes y metadata.", call. = FALSE)
}

for (mc in module_cols) {
  analysis_df[[mc]] <- safe_numeric(analysis_df[[mc]])
}

module_names_clean <- standardize_module_name(module_cols)
names(analysis_df)[match(module_cols, names(analysis_df))] <- module_names_clean
modules_use <- module_names_clean

traits_available <- intersect(TRAITS_USE, names(analysis_df))

if (length(traits_available) == 0) {
  stop("No hay traits disponibles para sensibilidad.", call. = FALSE)
}

core_modules_use <- intersect(CORE_MODULES, modules_use)

safe_write_csv(
  analysis_df,
  file.path(OUTDIR, "tables", "sensitivity_input_clean.csv")
)

cat("Samples:", nrow(analysis_df), "\n")
cat("Modules:", paste(modules_use, collapse = ", "), "\n")
cat("Core modules:", paste(core_modules_use, collapse = ", "), "\n")
cat("Traits:", paste(traits_available, collapse = ", "), "\n\n")

###############################################################################
# 6) FULL MODULE-TRAIT CORRELATIONS
###############################################################################

full_module_trait <- compute_module_trait(
  df = analysis_df,
  modules = modules_use,
  traits = traits_available,
  label = "full"
) %>%
  dplyr::mutate(
    FDR = p.adjust(p_value, method = "BH")
  )

safe_write_csv(
  full_module_trait,
  file.path(OUTDIR, "tables", "full_module_trait_correlations.csv")
)

###############################################################################
# 7) LEAVE-ONE-COUNTRY-OUT SENSITIVITY
###############################################################################

loco_country_results <- list()

if ("Country" %in% names(analysis_df)) {

  countries <- sort(unique(as.character(analysis_df$Country)))
  countries <- countries[!is.na(countries)]

  for (cc in countries) {

    df_cc <- analysis_df %>%
      dplyr::filter(as.character(Country) != cc)

    sens_tbl <- compute_module_trait(
      df = df_cc,
      modules = modules_use,
      traits = traits_available,
      label = paste0("leave_out_country_", cc)
    ) %>%
      dplyr::mutate(Excluded_country = cc)

    loco_country_results[[cc]] <- sens_tbl
  }

  loco_country_long <- dplyr::bind_rows(loco_country_results)

  loco_country_delta <- make_delta_summary(
    full_tbl = full_module_trait,
    sens_tbl = loco_country_long
  )

  loco_country_summary_by_module <- summarize_delta_by_module(loco_country_delta)
  loco_country_summary_by_module_trait <- summarize_delta_by_module_trait(loco_country_delta)

  safe_write_csv(
    loco_country_long,
    file.path(OUTDIR, "tables", "loco", "loco_country_module_trait_results.csv")
  )

  safe_write_csv(
    loco_country_delta,
    file.path(OUTDIR, "tables", "loco", "loco_country_delta_results.csv")
  )

  safe_write_csv(
    loco_country_summary_by_module,
    file.path(OUTDIR, "tables", "loco", "loco_country_summary_by_module.csv")
  )

  safe_write_csv(
    loco_country_summary_by_module_trait,
    file.path(OUTDIR, "tables", "loco", "loco_country_summary_by_module_trait.csv")
  )

  p_loco_all <- ggplot(
    loco_country_summary_by_module,
    aes(
      x = reorder(Module, mean_abs_delta_rho),
      y = mean_abs_delta_rho,
      fill = Module
    )
  ) +
    geom_col(show.legend = FALSE) +
    scale_fill_manual(values = get_module_colors(loco_country_summary_by_module$Module)) +
    coord_flip() +
    labs(
      title = "Leave-one-country-out robustness",
      x = "Module",
      y = "Mean absolute delta rho"
    ) +
    theme_bw(base_size = 12) +
    theme(panel.grid.major.y = element_blank())

  ggsave(
    file.path(OUTDIR, "figures", "loco", "loco_country_mean_abs_delta_rho_all_modules.pdf"),
    p_loco_all,
    width = 7,
    height = 5
  )

  ggsave(
    file.path(OUTDIR, "figures", "loco", "loco_country_mean_abs_delta_rho_all_modules.png"),
    p_loco_all,
    width = 7,
    height = 5,
    dpi = 300
  )

  if (length(core_modules_use) > 0) {
    p_loco_core <- loco_country_summary_by_module %>%
      dplyr::filter(Module %in% core_modules_use) %>%
      ggplot(aes(
        x = reorder(Module, mean_abs_delta_rho),
        y = mean_abs_delta_rho,
        fill = Module
      )) +
      geom_col(show.legend = FALSE) +
      scale_fill_manual(values = get_module_colors(core_modules_use)) +
      coord_flip() +
      labs(
        title = "Leave-one-country-out robustness — core modules",
        x = "Module",
        y = "Mean absolute delta rho"
      ) +
      theme_bw(base_size = 12) +
      theme(panel.grid.major.y = element_blank())

    ggsave(
      file.path(OUTDIR, "figures", "core_modules", "loco_country_core_modules.pdf"),
      p_loco_core,
      width = 6,
      height = 4
    )

    ggsave(
      file.path(OUTDIR, "figures", "core_modules", "loco_country_core_modules.png"),
      p_loco_core,
      width = 6,
      height = 4,
      dpi = 300
    )
  }

} else {

  loco_country_long <- tibble::tibble()
  loco_country_delta <- tibble::tibble()
  loco_country_summary_by_module <- tibble::tibble(
    Module = modules_use,
    mean_abs_delta_rho = NA_real_,
    median_abs_delta_rho = NA_real_,
    max_abs_delta_rho = NA_real_,
    direction_consistency = NA_real_,
    n_tests = NA_integer_
  )
  loco_country_summary_by_module_trait <- tibble::tibble()
}

###############################################################################
# 8) LEAVE-ONE-SITE-OUT SENSITIVITY
###############################################################################

site_var <- SITE_CANDIDATES[SITE_CANDIDATES %in% names(analysis_df)]
site_var <- ifelse(length(site_var) > 0, site_var[[1]], NA_character_)

loso_site_results <- list()

if (!is.na(site_var)) {

  sites <- sort(unique(as.character(analysis_df[[site_var]])))
  sites <- sites[!is.na(sites)]

  if (length(sites) >= 2) {

    for (ss in sites) {

      df_ss <- analysis_df %>%
        dplyr::filter(as.character(.data[[site_var]]) != ss)

      sens_tbl <- compute_module_trait(
        df = df_ss,
        modules = modules_use,
        traits = traits_available,
        label = paste0("leave_out_site_", ss)
      ) %>%
        dplyr::mutate(
          Excluded_site = ss,
          Site_variable = site_var
        )

      loso_site_results[[ss]] <- sens_tbl
    }

    loso_site_long <- dplyr::bind_rows(loso_site_results)

    loso_site_delta <- make_delta_summary(
      full_tbl = full_module_trait,
      sens_tbl = loso_site_long
    )

    loso_site_summary_by_module <- summarize_delta_by_module(loso_site_delta)
    loso_site_summary_by_module_trait <- summarize_delta_by_module_trait(loso_site_delta)

    safe_write_csv(
      loso_site_long,
      file.path(OUTDIR, "tables", "loso", "loso_site_module_trait_results.csv")
    )

    safe_write_csv(
      loso_site_delta,
      file.path(OUTDIR, "tables", "loso", "loso_site_delta_results.csv")
    )

    safe_write_csv(
      loso_site_summary_by_module,
      file.path(OUTDIR, "tables", "loso", "loso_site_summary_by_module.csv")
    )

    safe_write_csv(
      loso_site_summary_by_module_trait,
      file.path(OUTDIR, "tables", "loso", "loso_site_summary_by_module_trait.csv")
    )

    p_loso <- ggplot(
      loso_site_summary_by_module,
      aes(
        x = reorder(Module, mean_abs_delta_rho),
        y = mean_abs_delta_rho,
        fill = Module
      )
    ) +
      geom_col(show.legend = FALSE) +
      scale_fill_manual(values = get_module_colors(loso_site_summary_by_module$Module)) +
      coord_flip() +
      labs(
        title = paste0("Leave-one-site-out robustness: ", site_var),
        x = "Module",
        y = "Mean absolute delta rho"
      ) +
      theme_bw(base_size = 12) +
      theme(panel.grid.major.y = element_blank())

    ggsave(
      file.path(OUTDIR, "figures", "loso", "loso_site_mean_abs_delta_rho.pdf"),
      p_loso,
      width = 7,
      height = 5
    )

    ggsave(
      file.path(OUTDIR, "figures", "loso", "loso_site_mean_abs_delta_rho.png"),
      p_loso,
      width = 7,
      height = 5,
      dpi = 300
    )

  } else {

    loso_site_long <- tibble::tibble()
    loso_site_delta <- tibble::tibble()
    loso_site_summary_by_module <- tibble::tibble()
    loso_site_summary_by_module_trait <- tibble::tibble()
  }

} else {

  loso_site_long <- tibble::tibble()
  loso_site_delta <- tibble::tibble()
  loso_site_summary_by_module <- tibble::tibble()
  loso_site_summary_by_module_trait <- tibble::tibble()
}

###############################################################################
# 9) BALANCED DOWNSAMPLING BY COUNTRY AND DIAGNOSIS
###############################################################################

set.seed(SEED)

downsampling_results <- list()

can_downsample <- all(c("Country", "SampleGroup") %in% names(analysis_df))

if (can_downsample) {

  strata_counts <- analysis_df %>%
    dplyr::count(Country, SampleGroup, name = "n") %>%
    tidyr::drop_na()

  min_stratum_n <- min(strata_counts$n, na.rm = TRUE)

  if (is.finite(min_stratum_n) && min_stratum_n >= MIN_N_FOR_CORR) {

    for (iter in seq_len(N_DOWNSAMPLING_ITER)) {

      sampled_ids <- analysis_df %>%
        dplyr::group_by(Country, SampleGroup) %>%
        dplyr::slice_sample(n = min_stratum_n, replace = FALSE) %>%
        dplyr::ungroup() %>%
        dplyr::pull(SampleId)

      df_iter <- analysis_df %>%
        dplyr::filter(SampleId %in% sampled_ids)

      iter_tbl <- compute_module_trait(
        df = df_iter,
        modules = modules_use,
        traits = traits_available,
        label = paste0("balanced_downsampling_", iter)
      ) %>%
        dplyr::mutate(
          Iteration = iter,
          n_per_country_group = min_stratum_n,
          n_sampled = nrow(df_iter)
        )

      downsampling_results[[iter]] <- iter_tbl
    }

    downsampling_long <- dplyr::bind_rows(downsampling_results)

    downsampling_delta <- make_delta_summary(
      full_tbl = full_module_trait,
      sens_tbl = downsampling_long
    )

    downsampling_summary_by_module <- summarize_delta_by_module(downsampling_delta)
    downsampling_summary_by_module_trait <- summarize_delta_by_module_trait(downsampling_delta)

    safe_write_csv(
      strata_counts,
      file.path(OUTDIR, "tables", "downsampling", "country_group_strata_counts.csv")
    )

    safe_write_csv(
      downsampling_long,
      file.path(OUTDIR, "tables", "downsampling", "balanced_downsampling_module_trait_results.csv")
    )

    safe_write_csv(
      downsampling_delta,
      file.path(OUTDIR, "tables", "downsampling", "balanced_downsampling_delta_results.csv")
    )

    safe_write_csv(
      downsampling_summary_by_module,
      file.path(OUTDIR, "tables", "downsampling", "balanced_downsampling_summary_by_module.csv")
    )

    safe_write_csv(
      downsampling_summary_by_module_trait,
      file.path(OUTDIR, "tables", "downsampling", "balanced_downsampling_summary_by_module_trait.csv")
    )

    p_down <- ggplot(
      downsampling_summary_by_module,
      aes(
        x = reorder(Module, mean_abs_delta_rho),
        y = mean_abs_delta_rho,
        fill = Module
      )
    ) +
      geom_col(show.legend = FALSE) +
      scale_fill_manual(values = get_module_colors(downsampling_summary_by_module$Module)) +
      coord_flip() +
      labs(
        title = "Balanced downsampling robustness",
        subtitle = paste0(
          N_DOWNSAMPLING_ITER,
          " iterations; n per country-diagnosis stratum = ",
          min_stratum_n
        ),
        x = "Module",
        y = "Mean absolute delta rho"
      ) +
      theme_bw(base_size = 12) +
      theme(panel.grid.major.y = element_blank())

    ggsave(
      file.path(OUTDIR, "figures", "downsampling", "balanced_downsampling_mean_abs_delta_rho.pdf"),
      p_down,
      width = 7,
      height = 5
    )

    ggsave(
      file.path(OUTDIR, "figures", "downsampling", "balanced_downsampling_mean_abs_delta_rho.png"),
      p_down,
      width = 7,
      height = 5,
      dpi = 300
    )

  } else {

    downsampling_long <- tibble::tibble()
    downsampling_delta <- tibble::tibble()
    downsampling_summary_by_module <- tibble::tibble()
    downsampling_summary_by_module_trait <- tibble::tibble()

    warning(
      "Balanced downsampling was skipped because the minimum country-diagnosis stratum is too small: ",
      min_stratum_n
    )
  }

} else {

  downsampling_long <- tibble::tibble()
  downsampling_delta <- tibble::tibble()
  downsampling_summary_by_module <- tibble::tibble()
  downsampling_summary_by_module_trait <- tibble::tibble()
}

###############################################################################
# 10) CORE MODULE STABILITY SUMMARY
###############################################################################

core_loco <- loco_country_summary_by_module %>%
  dplyr::filter(Module %in% core_modules_use) %>%
  dplyr::mutate(Analysis = "LOCO_country")

core_loso <- loso_site_summary_by_module %>%
  dplyr::filter(Module %in% core_modules_use) %>%
  dplyr::mutate(Analysis = "LOSO_site")

core_down <- downsampling_summary_by_module %>%
  dplyr::filter(Module %in% core_modules_use) %>%
  dplyr::mutate(Analysis = "Balanced_downsampling")

core_stability_summary <- dplyr::bind_rows(
  core_loco,
  core_loso,
  core_down
) %>%
  dplyr::select(
    Analysis,
    Module,
    dplyr::everything()
  )

safe_write_csv(
  core_stability_summary,
  file.path(OUTDIR, "tables", "core_modules_stability_summary.csv")
)

if (nrow(core_stability_summary) > 0) {

  p_core <- ggplot(
    core_stability_summary,
    aes(
      x = Module,
      y = mean_abs_delta_rho,
      fill = Module
    )
  ) +
    geom_col(show.legend = FALSE) +
    scale_fill_manual(values = get_module_colors(core_stability_summary$Module)) +
    facet_wrap(~ Analysis, scales = "free_y") +
    labs(
      title = "Core WGCNA module stability",
      x = "Module",
      y = "Mean absolute delta rho"
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major.x = element_blank(),
      strip.text = element_text(face = "bold")
    )

  ggsave(
    file.path(OUTDIR, "figures", "core_modules", "core_modules_stability_summary.pdf"),
    p_core,
    width = 9,
    height = 4
  )

  ggsave(
    file.path(OUTDIR, "figures", "core_modules", "core_modules_stability_summary.png"),
    p_core,
    width = 9,
    height = 4,
    dpi = 300
  )
}

###############################################################################
# 11) FINAL SUMMARY
###############################################################################

script5_summary <- tibble::tibble(
  metric = c(
    "base_dir",
    "outdir",
    "wgcna_input_level",
    "n_samples",
    "n_modules",
    "n_traits",
    "core_modules",
    "country_variable_available",
    "site_variable_used",
    "n_downsampling_iterations_requested",
    "random_seed"
  ),
  value = c(
    BASE_DIR,
    OUTDIR,
    "GENE-COLLAPSED, not aptamer-level",
    as.character(nrow(analysis_df)),
    as.character(length(modules_use)),
    as.character(length(traits_available)),
    paste(core_modules_use, collapse = ", "),
    as.character("Country" %in% names(analysis_df)),
    as.character(site_var),
    as.character(N_DOWNSAMPLING_ITER),
    as.character(SEED)
  )
)

safe_write_csv(
  script5_summary,
  file.path(OUTDIR, "tables", "script5_final_summary.csv")
)

save(
  analysis_df,
  full_module_trait,
  loco_country_long,
  loco_country_delta,
  loco_country_summary_by_module,
  loco_country_summary_by_module_trait,
  loso_site_long,
  loso_site_delta,
  loso_site_summary_by_module,
  loso_site_summary_by_module_trait,
  downsampling_long,
  downsampling_delta,
  downsampling_summary_by_module,
  downsampling_summary_by_module_trait,
  core_stability_summary,
  script5_summary,
  file = file.path(OUTDIR, "workspace", "script5_sensitivity_workspace.RData")
)

writeLines(
  capture.output(utils::sessionInfo()),
  con = file.path(OUTDIR, "sessionInfo.txt")
)

###############################################################################
# 12) FINAL MESSAGE
###############################################################################

cat("\nScript 05 terminado correctamente.\n")
cat("Directorio principal de salida:\n", OUTDIR, "\n")
cat("Input level: GENE-COLLAPSED, not aptamer-level.\n")
cat("Core modules:", paste(core_modules_use, collapse = ", "), "\n")
cat("\nArchivos clave:\n")
cat("- tables/full_module_trait_correlations.csv\n")
cat("- tables/loco/loco_country_summary_by_module.csv\n")
cat("- tables/loso/loso_site_summary_by_module.csv\n")
cat("- tables/downsampling/balanced_downsampling_summary_by_module.csv\n")
cat("- tables/core_modules_stability_summary.csv\n")
cat("- figures/core_modules/core_modules_stability_summary.pdf/png\n")

###############################################################################
# END
###############################################################################

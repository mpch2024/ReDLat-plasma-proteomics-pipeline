###############################################################################
# ReDLat plasma proteomics — WGCNA workflow
# 08. Assess structural preservation
# Requires: outputs from Scripts 02–04
# Produces: fixed-gene country and site preservation results
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

suppressPackageStartupMessages(
  library(WGCNA)
)

WGCNA::allowWGCNAThreads()

options(stringsAsFactors = FALSE)
options(error = traceback)

###############################################################################
# 2) PATHS
###############################################################################

BASE_DIR <- WGCNA_CONFIG$project_root

SCRIPT11_DIR <- file.path(WGCNA_CONFIG$result_root,
  "02_network"
)

SCRIPT13_DIR <- file.path(WGCNA_CONFIG$result_root,
  "04_module_traits"
)

SCRIPT15_DIR <- file.path(WGCNA_CONFIG$result_root,
  "15_WGCNA_structural_module_preservation_outcome_independent"
)

WORKSPACE_FILE <- file.path(
  SCRIPT11_DIR,
  "workspace",
  "wgcna_core_light_workspace.RData"
)

MODULE_REFERENCE_FILE <- file.path(
  SCRIPT13_DIR,
  "tables",
  "module_biological_label_reference.csv"
)

OLD_COUNTRY_FILE <- file.path(
  SCRIPT15_DIR,
  "tables",
  "country",
  "country_module_preservation_primary_statistics.csv"
)

OLD_SITE_FILE <- file.path(
  SCRIPT15_DIR,
  "tables",
  "site",
  "site_module_preservation_primary_statistics.csv"
)

OUTDIR <- file.path(WGCNA_CONFIG$result_root,
  "08_preservation"
)

SUBDIRS <- c(
  "tables",
  "tables/input",
  "tables/country",
  "tables/site",
  "tables/comparison",
  "figures",
  "figures/country",
  "figures/site",
  "figures/comparison",
  "workspace",
  "workspace/raw_module_preservation_objects",
  "workspace/permutation_statistics",
  "logs"
)

invisible(lapply(
  file.path(OUTDIR, SUBDIRS),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

required_files <- c(
  WORKSPACE_FILE,
  MODULE_REFERENCE_FILE
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  stop(
    "Missing required files:\n",
    paste(missing_files, collapse = "\n"),
    "\n\nRun Scripts 11 and 13 first.",
    call. = FALSE
  )
}

###############################################################################
# 3) PARAMETERS
###############################################################################

EXPECTED_N_SAMPLES <- 639L
EXPECTED_N_GENES_FULL <- 9638L
EXPECTED_N_MODULES <- 8L
EXPECTED_N_COUNTRIES <- 5L
EXPECTED_N_SITE_LEVELS <- 7L
EXPECTED_N_MULTISITE_COUNTRIES <- 2L
EXPECTED_N_RECIPROCAL_SITE_RUNS <- 4L
EXPECTED_N_FIXED_GENES <- 4239L

REQUIRE_EXPECTED_DIMENSIONS <- TRUE
REQUIRE_EXPECTED_FIXED_GENE_COUNT <- TRUE

NETWORK_TYPE <- "signed"
COR_FUNCTION <- "cor"
COR_OPTIONS <- "use = 'p'"

N_PERMUTATIONS <- 100L

# One seed selects the fixed module gene sets.
FIXED_GENESET_SEED <- 15501L

# One common preservation seed is deliberately reused in every context.
# Because the expression matrices differ, the analyses remain context-specific,
# while the random permutation mechanism is held constant.
MODULE_PRESERVATION_SEED <- 15511L

MAX_FIXED_GENES_PER_MODULE <- 1000L
MAX_MODULE_SIZE <- 1000L
MAX_GOLD_MODULE_SIZE <- 1000L

QUICK_COR <- 1
CALCULATE_CLUSTER_COEFFICIENT <- FALSE
CALCULATE_COR_KIM_ALL <- FALSE
INCLUDE_KME_ALL_IN_SUMMARY <- FALSE
PARALLEL_MODULE_PRESERVATION <- FALSE

RESUME_EXISTING_RUNS <- TRUE
OVERWRITE_EXISTING_RUNS <- FALSE

MIN_REFERENCE_SAMPLES <- 30L
MIN_TEST_SAMPLES <- 30L

GREY_NAME <- "grey"
GOLD_NAME <- "gold"

CHECKPOINT_VERSION_TAG <- "FIXEDGENESET_V1"

MODULE_COLORS <- c(
  black = "#2B2B2B",
  brown = "#9C6B00",
  yellow = "#B7D500",
  blue = "#1535E8",
  green = "#3E8F4E",
  red = "#B4334A",
  purple = "#8A2BE2",
  magenta = "#C75ACD",
  pink = "#F3A6D6",
  greenyellow = "#ADFF2F",
  turquoise = "#40E0D0",
  cyan = "#00BFC4",
  tan = "#D2B48C",
  salmon = "#FA8072",
  midnightblue = "#191970",
  lightcyan = "#E0FFFF",
  grey = "#9E9E9E"
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

safe_proportion <- function(x) {
  x <- x[!is.na(x)]

  if (length(x) == 0) {
    return(NA_real_)
  }

  mean(as.logical(x))
}

sanitize_file_tag <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
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
  candidates <- c(
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

classify_zsummary <- function(z) {
  dplyr::case_when(
    is.na(z) ~ "Not estimable",
    z >= 10 ~ "Strong evidence",
    z >= 2 ~ "Weak-to-moderate evidence",
    TRUE ~ "No evidence"
  )
}

###############################################################################
# 5) FIXED GENE-SET CONSTRUCTION
###############################################################################

build_fixed_gene_manifest <- function(
    module_colors,
    modules,
    max_genes_per_module,
    seed
) {
  set.seed(seed)

  purrr::map_dfr(
    modules,
    function(mod) {
      genes_mod <- names(
        module_colors
      )[
        module_colors == mod
      ]

      genes_mod <- sort(
        genes_mod
      )

      original_size <- length(
        genes_mod
      )

      if (
        original_size >
          max_genes_per_module
      ) {
        selected <- sample(
          genes_mod,
          size =
            max_genes_per_module,
          replace = FALSE
        )

        selected <- sort(selected)
        selection_method <-
          "Fixed random subset"
      } else {
        selected <- genes_mod
        selection_method <-
          "All module genes"
      }

      tibble::tibble(
        Module = mod,
        Gene = selected,
        Original_module_size =
          original_size,
        Fixed_module_size =
          length(selected),
        Module_was_restricted =
          original_size >
          max_genes_per_module,
        Selection_method =
          selection_method,
        Fixed_geneset_seed = seed,
        Max_fixed_genes_per_module =
          max_genes_per_module
      )
    }
  )
}

###############################################################################
# 6) MODULE-PRESERVATION OUTPUT EXTRACTION
###############################################################################

extract_nested_pair <- function(
    x,
    reference_index = 1L,
    test_index = 2L
) {
  if (is.null(x)) {
    return(NULL)
  }

  if (length(x) < reference_index) {
    return(NULL)
  }

  ref_component <- x[[reference_index]]

  if (
    is.null(ref_component) ||
    length(ref_component) < test_index
  ) {
    return(NULL)
  }

  ref_component[[test_index]]
}

matrix_to_module_table <- function(
    x,
    prefix = NULL
) {
  if (is.null(x)) {
    return(tibble::tibble())
  }

  x <- as.data.frame(
    x,
    check.names = FALSE
  )

  if (nrow(x) == 0) {
    return(tibble::tibble())
  }

  out <- x %>%
    tibble::rownames_to_column(
      "Module"
    ) %>%
    tibble::as_tibble()

  if (!is.null(prefix)) {
    cols_rename <- setdiff(
      names(out),
      "Module"
    )

    names(out)[
      match(
        cols_rename,
        names(out)
      )
    ] <- paste0(
      prefix,
      cols_rename
    )
  }

  out
}

extract_module_preservation_tables <- function(mp_object) {
  preservation_z <- extract_nested_pair(
    mp_object$preservation$Z
  )

  preservation_observed <- extract_nested_pair(
    mp_object$preservation$observed
  )

  preservation_log_p <- extract_nested_pair(
    mp_object$preservation$log.p
  )

  preservation_log_p_bonf <- extract_nested_pair(
    mp_object$preservation$log.pBonf
  )

  quality_z <- extract_nested_pair(
    mp_object$quality$Z
  )

  quality_observed <- extract_nested_pair(
    mp_object$quality$observed
  )

  reference_separability_z <- extract_nested_pair(
    mp_object$referenceSeparability$Z
  )

  test_separability_z <- extract_nested_pair(
    mp_object$testSeparability$Z
  )

  tbl <- matrix_to_module_table(
    preservation_z,
    prefix = NULL
  )

  if (nrow(tbl) == 0) {
    stop(
      "WGCNA modulePreservation did not return a preservation Z table.",
      call. = FALSE
    )
  }

  other_tables <- list(
    matrix_to_module_table(
      preservation_observed,
      prefix = "observed."
    ),
    matrix_to_module_table(
      preservation_log_p,
      prefix = "log10p."
    ),
    matrix_to_module_table(
      preservation_log_p_bonf,
      prefix = "log10pBonf."
    ),
    matrix_to_module_table(
      quality_z,
      prefix = "quality."
    ),
    matrix_to_module_table(
      quality_observed,
      prefix = "qualityObserved."
    ),
    matrix_to_module_table(
      reference_separability_z,
      prefix = "referenceSeparability."
    ),
    matrix_to_module_table(
      test_separability_z,
      prefix = "testSeparability."
    )
  )

  for (other_tbl in other_tables) {
    if (nrow(other_tbl) > 0) {
      tbl <- tbl %>%
        dplyr::left_join(
          other_tbl,
          by = "Module"
        )
    }
  }

  tbl
}

###############################################################################
# 7) FIXED-GENE-SET MODULE-PRESERVATION RUNNER
###############################################################################

run_fixed_module_preservation_pair <- function(
    reference_expr,
    test_expr,
    fixed_module_colors,
    run_id,
    analysis_type,
    reference_label,
    test_label
) {
  if (nrow(reference_expr) < MIN_REFERENCE_SAMPLES) {
    stop(
      "Reference set ",
      reference_label,
      " contains only ",
      nrow(reference_expr),
      " samples.",
      call. = FALSE
    )
  }

  if (nrow(test_expr) < MIN_TEST_SAMPLES) {
    stop(
      "Test set ",
      test_label,
      " contains only ",
      nrow(test_expr),
      " samples.",
      call. = FALSE
    )
  }

  common_genes <- Reduce(
    intersect,
    list(
      colnames(reference_expr),
      colnames(test_expr),
      names(fixed_module_colors)
    )
  )

  if (length(common_genes) == 0) {
    stop(
      "No common fixed genes were found for run ",
      run_id,
      ".",
      call. = FALSE
    )
  }

  # Preserve the exact preselected fixed-gene order.
  fixed_order <- names(
    fixed_module_colors
  )

  common_genes <- fixed_order[
    fixed_order %in% common_genes
  ]

  reference_expr <- as.matrix(
    reference_expr[
      ,
      common_genes,
      drop = FALSE
    ]
  )

  test_expr <- as.matrix(
    test_expr[
      ,
      common_genes,
      drop = FALSE
    ]
  )

  colors_use <- fixed_module_colors[
    common_genes
  ]

  if (
    !identical(
      colnames(reference_expr),
      colnames(test_expr)
    ) ||
    !identical(
      colnames(reference_expr),
      names(colors_use)
    )
  ) {
    stop(
      "Fixed-gene alignment failed for run ",
      run_id,
      ".",
      call. = FALSE
    )
  }

  if (
    anyNA(reference_expr) ||
    anyNA(test_expr)
  ) {
    stop(
      "Missing expression values detected in run ",
      run_id,
      ".",
      call. = FALSE
    )
  }

  run_module_sizes <- table(
    colors_use
  )

  if (any(
    run_module_sizes >
      MAX_MODULE_SIZE
  )) {
    stop(
      "At least one fixed module still exceeds MAX_MODULE_SIZE in run ",
      run_id,
      ".",
      call. = FALSE
    )
  }

  run_tag <- paste(
    CHECKPOINT_VERSION_TAG,
    sanitize_file_tag(run_id),
    sep = "_"
  )

  raw_object_file <- file.path(
    OUTDIR,
    "workspace",
    "raw_module_preservation_objects",
    paste0(
      run_tag,
      "_modulePreservation.rds"
    )
  )

  permutation_file <- file.path(
    OUTDIR,
    "workspace",
    "permutation_statistics",
    paste0(
      run_tag,
      "_permuted_statistics.RData"
    )
  )

  result_csv <- file.path(
    OUTDIR,
    "workspace",
    "raw_module_preservation_objects",
    paste0(
      run_tag,
      "_extracted_results.csv"
    )
  )

  can_resume <- (
    RESUME_EXISTING_RUNS &&
      !OVERWRITE_EXISTING_RUNS &&
      file.exists(raw_object_file) &&
      file.exists(result_csv)
  )

  if (can_resume) {
    cat(
      "Resuming completed fixed-gene preservation run: ",
      run_id,
      "\n",
      sep = ""
    )

    mp_object <- readRDS(
      raw_object_file
    )

    extracted <- safe_read_csv(
      result_csv
    )

    return(list(
      object = mp_object,
      table = extracted,
      resumed = TRUE,
      raw_object_file = raw_object_file,
      permutation_file = permutation_file
    ))
  }

  if (
    OVERWRITE_EXISTING_RUNS &&
    file.exists(raw_object_file)
  ) {
    unlink(
      raw_object_file,
      force = TRUE
    )
  }

  if (file.exists(permutation_file)) {
    unlink(
      permutation_file,
      force = TRUE
    )
  }

  cat("\n")
  cat("============================================================\n")
  cat("Fixed-gene structural preservation run: ", run_id, "\n", sep = "")
  cat("Analysis type: ", analysis_type, "\n", sep = "")
  cat("Reference: ", reference_label, " (n=", nrow(reference_expr), ")\n", sep = "")
  cat("Test: ", test_label, " (n=", nrow(test_expr), ")\n", sep = "")
  cat("Fixed genes: ", ncol(reference_expr), "\n", sep = "")
  cat("Permutations: ", N_PERMUTATIONS, "\n", sep = "")
  cat("Common preservation seed: ", MODULE_PRESERVATION_SEED, "\n", sep = "")
  cat("============================================================\n")

  multi_data <- list(
    reference = list(
      data = reference_expr
    ),
    test = list(
      data = test_expr
    )
  )

  multi_color <- list(
    reference = colors_use
  )

  start_time <- Sys.time()

  mp_object <- WGCNA::modulePreservation(
    multiData = multi_data,
    multiColor = multi_color,
    dataIsExpr = TRUE,
    networkType = NETWORK_TYPE,
    corFnc = COR_FUNCTION,
    corOptions = COR_OPTIONS,
    referenceNetworks = 1,
    testNetworks = 2,
    nPermutations = N_PERMUTATIONS,
    includekMEallInSummary =
      INCLUDE_KME_ALL_IN_SUMMARY,
    restrictSummaryForGeneralNetworks = TRUE,
    calculateQvalue = FALSE,
    randomSeed =
      MODULE_PRESERVATION_SEED,
    maxGoldModuleSize =
      MAX_GOLD_MODULE_SIZE,
    maxModuleSize =
      MAX_MODULE_SIZE,
    quickCor = QUICK_COR,
    calculateCor.kIMall =
      CALCULATE_COR_KIM_ALL,
    calculateClusterCoeff =
      CALCULATE_CLUSTER_COEFFICIENT,
    useInterpolation = FALSE,
    checkData = TRUE,
    greyName = GREY_NAME,
    goldName = GOLD_NAME,
    savePermutedStatistics = TRUE,
    loadPermutedStatistics = FALSE,
    permutedStatisticsFile =
      permutation_file,
    discardInvalidOutput = TRUE,
    parallelCalculation =
      PARALLEL_MODULE_PRESERVATION,
    verbose = 3
  )

  end_time <- Sys.time()

  saveRDS(
    mp_object,
    raw_object_file
  )

  extracted <- extract_module_preservation_tables(
    mp_object
  ) %>%
    dplyr::mutate(
      Analysis_type = analysis_type,
      Run_ID = run_id,
      Reference_label =
        reference_label,
      Test_label = test_label,
      N_reference =
        nrow(reference_expr),
      N_test = nrow(test_expr),
      N_fixed_genes =
        ncol(reference_expr),
      N_permutations =
        N_PERMUTATIONS,
      Fixed_geneset_seed =
        FIXED_GENESET_SEED,
      Module_preservation_seed =
        MODULE_PRESERVATION_SEED,
      Max_module_size =
        MAX_MODULE_SIZE,
      Max_gold_module_size =
        MAX_GOLD_MODULE_SIZE,
      Checkpoint_version =
        CHECKPOINT_VERSION_TAG,
      Runtime_minutes =
        as.numeric(
          difftime(
            end_time,
            start_time,
            units = "mins"
          )
        ),
      Resumed_existing_run =
        FALSE
    )

  safe_write_csv(
    extracted,
    result_csv
  )

  list(
    object = mp_object,
    table = extracted,
    resumed = FALSE,
    raw_object_file = raw_object_file,
    permutation_file = permutation_file
  )
}

###############################################################################
# 8) LOAD SCRIPT 11 WORKSPACE AND MODULE REFERENCE
###############################################################################

load(
  WORKSPACE_FILE
)

required_objects <- c(
  "datExpr_clean",
  "meta_final",
  "softPower",
  "mergedColors",
  "module_count_tbl",
  "gene_module_assignment"
)

missing_objects <- required_objects[
  !vapply(
    required_objects,
    exists,
    logical(1),
    inherits = FALSE
  )
]

if (length(missing_objects) > 0) {
  stop(
    "The Script 11 light workspace is missing required objects: ",
    paste(
      missing_objects,
      collapse = ", "
    ),
    call. = FALSE
  )
}

module_reference <- safe_read_csv(
  MODULE_REFERENCE_FILE
)

###############################################################################
# 9) ALIGNMENT AUDIT
###############################################################################

datExpr_clean <- as.data.frame(
  datExpr_clean,
  check.names = FALSE
)

if (!"SampleId" %in% names(meta_final)) {
  stop(
    "meta_final does not contain SampleId.",
    call. = FALSE
  )
}

meta_final <- meta_final %>%
  dplyr::mutate(
    SampleId = as.character(
      SampleId
    )
  )

if (!identical(
  rownames(datExpr_clean),
  meta_final$SampleId
)) {
  stop(
    "datExpr_clean rows are not aligned with meta_final$SampleId.",
    call. = FALSE
  )
}

if (
  is.null(names(mergedColors)) ||
  !identical(
    names(mergedColors),
    colnames(datExpr_clean)
  )
) {
  stop(
    "mergedColors is not aligned with datExpr_clean columns.",
    call. = FALSE
  )
}

if (anyNA(datExpr_clean)) {
  stop(
    "datExpr_clean contains missing values.",
    call. = FALSE
  )
}

if (!"Country" %in% names(meta_final)) {
  stop(
    "Country was not found in meta_final.",
    call. = FALSE
  )
}

meta_final$Country <- droplevels(
  factor(
    clean_text_na(
      meta_final$Country
    )
  )
)

site_var <- detect_site_variable(
  meta_final
)

if (
  is.na(site_var) ||
  !site_var %in% names(meta_final)
) {
  stop(
    "A recruitment-site variable was not detected in meta_final.",
    call. = FALSE
  )
}

meta_final[[site_var]] <- droplevels(
  factor(
    clean_text_na(
      meta_final[[site_var]]
    )
  )
)

if ("Country_numeric" %in% names(meta_final)) {
  meta_final$Country_numeric <- NULL
}

full_module_colors <- as.character(
  mergedColors
)

names(full_module_colors) <- names(
  mergedColors
)

modules_use <- sort(
  unique(
    full_module_colors[
      !full_module_colors %in%
        c(
          GREY_NAME,
          GOLD_NAME
        )
    ]
  )
)

module_reference <- module_reference %>%
  dplyr::mutate(
    Module = as.character(Module)
  )

countries <- sort(
  unique(
    as.character(
      meta_final$Country
    )
  )
)

countries <- countries[
  !is.na(countries)
]

country_site_counts <- meta_final %>%
  dplyr::filter(
    !is.na(Country),
    !is.na(
      .data[[site_var]]
    )
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
    name = "N_sites"
  )

multisite_countries <- sites_per_country %>%
  dplyr::filter(
    N_sites >= 2
  ) %>%
  dplyr::pull(Country) %>%
  as.character()

site_pair_plan <- purrr::map_dfr(
  multisite_countries,
  function(country_value) {
    sites_country <- country_site_counts %>%
      dplyr::filter(
        as.character(Country) ==
          country_value
      ) %>%
      dplyr::pull(Site) %>%
      as.character()

    expand.grid(
      Reference_site =
        sites_country,
      Test_site =
        sites_country,
      stringsAsFactors = FALSE
    ) %>%
      dplyr::filter(
        Reference_site !=
          Test_site
      ) %>%
      dplyr::mutate(
        Country = country_value,
        Run_ID = paste0(
          "SITE_",
          sanitize_file_tag(
            country_value
          ),
          "_",
          sanitize_file_tag(
            Reference_site
          ),
          "_to_",
          sanitize_file_tag(
            Test_site
          )
        ),
        Reference_label = paste(
          country_value,
          Reference_site,
          sep = "::"
        ),
        Test_label = paste(
          country_value,
          Test_site,
          sep = "::"
        )
      ) %>%
      dplyr::select(
        Country,
        Reference_site,
        Test_site,
        Reference_label,
        Test_label,
        Run_ID
      )
  }
)

###############################################################################
# 10) CREATE THE SINGLE FIXED GENE SET
###############################################################################

fixed_gene_manifest <- build_fixed_gene_manifest(
  module_colors =
    full_module_colors,
  modules = modules_use,
  max_genes_per_module =
    MAX_FIXED_GENES_PER_MODULE,
  seed =
    FIXED_GENESET_SEED
)

fixed_module_summary <- fixed_gene_manifest %>%
  dplyr::group_by(Module) %>%
  dplyr::summarise(
    Original_module_size =
      dplyr::first(
        Original_module_size
      ),
    Fixed_module_size =
      dplyr::n(),
    Module_was_restricted =
      dplyr::first(
        Module_was_restricted
      ),
    Selection_method =
      dplyr::first(
        Selection_method
      ),
    Fixed_geneset_seed =
      dplyr::first(
        Fixed_geneset_seed
      ),
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    module_reference,
    by = "Module"
  ) %>%
  dplyr::arrange(
    dplyr::desc(
      Original_module_size
    )
  )

fixed_genes <- fixed_gene_manifest$Gene

if (anyDuplicated(fixed_genes) > 0) {
  stop(
    "The fixed gene manifest contains duplicated genes.",
    call. = FALSE
  )
}

missing_fixed_genes <- setdiff(
  fixed_genes,
  colnames(datExpr_clean)
)

if (length(missing_fixed_genes) > 0) {
  stop(
    "The fixed manifest contains genes absent from datExpr_clean.",
    call. = FALSE
  )
}

datExpr_fixed <- datExpr_clean[
  ,
  fixed_genes,
  drop = FALSE
]

fixed_module_colors <- fixed_gene_manifest$Module
names(fixed_module_colors) <-
  fixed_gene_manifest$Gene

if (!identical(
  colnames(datExpr_fixed),
  names(fixed_module_colors)
)) {
  stop(
    "The fixed expression matrix and fixed module colors are not aligned.",
    call. = FALSE
  )
}

fixed_module_sizes <- table(
  fixed_module_colors
)

if (any(
  fixed_module_sizes >
    MAX_MODULE_SIZE
)) {
  stop(
    "At least one fixed module exceeds MAX_MODULE_SIZE.",
    call. = FALSE
  )
}

if (REQUIRE_EXPECTED_DIMENSIONS) {
  if (nrow(datExpr_clean) !=
      EXPECTED_N_SAMPLES) {
    stop(
      "Expected ",
      EXPECTED_N_SAMPLES,
      " samples, but found ",
      nrow(datExpr_clean),
      ".",
      call. = FALSE
    )
  }

  if (ncol(datExpr_clean) !=
      EXPECTED_N_GENES_FULL) {
    stop(
      "Expected ",
      EXPECTED_N_GENES_FULL,
      " full-network genes, but found ",
      ncol(datExpr_clean),
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
      meta_final[[site_var]]
    ) !=
      EXPECTED_N_SITE_LEVELS
  ) {
    stop(
      "Expected ",
      EXPECTED_N_SITE_LEVELS,
      " site levels, but found ",
      dplyr::n_distinct(
        meta_final[[site_var]]
      ),
      ".",
      call. = FALSE
    )
  }

  if (length(multisite_countries) !=
      EXPECTED_N_MULTISITE_COUNTRIES) {
    stop(
      "Expected ",
      EXPECTED_N_MULTISITE_COUNTRIES,
      " multisite countries, but found ",
      length(multisite_countries),
      ".",
      call. = FALSE
    )
  }

  if (nrow(site_pair_plan) !=
      EXPECTED_N_RECIPROCAL_SITE_RUNS) {
    stop(
      "Expected ",
      EXPECTED_N_RECIPROCAL_SITE_RUNS,
      " reciprocal site runs, but found ",
      nrow(site_pair_plan),
      ".",
      call. = FALSE
    )
  }
}

if (
  REQUIRE_EXPECTED_FIXED_GENE_COUNT &&
  ncol(datExpr_fixed) !=
    EXPECTED_N_FIXED_GENES
) {
  stop(
    "Expected ",
    EXPECTED_N_FIXED_GENES,
    " fixed genes, but found ",
    ncol(datExpr_fixed),
    ". Review module sizes and the fixed-gene manifest.",
    call. = FALSE
  )
}

safe_write_csv(
  fixed_gene_manifest,
  file.path(
    OUTDIR,
    "tables",
    "input",
    "fixed_gene_manifest_all_modules.csv"
  )
)

safe_write_csv(
  fixed_module_summary,
  file.path(
    OUTDIR,
    "tables",
    "input",
    "fixed_gene_manifest_summary_by_module.csv"
  )
)

safe_write_csv(
  country_site_counts,
  file.path(
    OUTDIR,
    "tables",
    "input",
    "country_site_sample_counts.csv"
  )
)

safe_write_csv(
  site_pair_plan,
  file.path(
    OUTDIR,
    "tables",
    "input",
    "reciprocal_within_country_site_run_plan.csv"
  )
)

input_audit <- tibble::tibble(
  metric = c(
    "base_dir",
    "workspace_file",
    "output_dir",
    "n_samples",
    "n_full_network_genes",
    "n_fixed_network_genes",
    "n_modules",
    "modules",
    "fixed_module_sizes",
    "restricted_modules",
    "fixed_geneset_seed",
    "common_module_preservation_seed",
    "checkpoint_version",
    "n_permutations",
    "max_module_size",
    "any_fixed_module_exceeds_max",
    "n_countries",
    "countries",
    "site_variable",
    "n_site_levels",
    "multisite_countries",
    "n_reciprocal_site_runs",
    "Country_numeric_present",
    "wgcna_input_level"
  ),
  value = c(
    BASE_DIR,
    WORKSPACE_FILE,
    OUTDIR,
    as.character(
      nrow(datExpr_fixed)
    ),
    as.character(
      ncol(datExpr_clean)
    ),
    as.character(
      ncol(datExpr_fixed)
    ),
    as.character(
      length(modules_use)
    ),
    paste(
      modules_use,
      collapse = ", "
    ),
    paste0(
      names(fixed_module_sizes),
      "=",
      as.integer(
        fixed_module_sizes
      ),
      collapse = "; "
    ),
    paste(
      fixed_module_summary$Module[
        fixed_module_summary$
          Module_was_restricted
      ],
      collapse = ", "
    ),
    as.character(
      FIXED_GENESET_SEED
    ),
    as.character(
      MODULE_PRESERVATION_SEED
    ),
    CHECKPOINT_VERSION_TAG,
    as.character(
      N_PERMUTATIONS
    ),
    as.character(
      MAX_MODULE_SIZE
    ),
    as.character(
      any(
        fixed_module_sizes >
          MAX_MODULE_SIZE
      )
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
        meta_final[[site_var]]
      )
    ),
    paste(
      multisite_countries,
      collapse = ", "
    ),
    as.character(
      nrow(site_pair_plan)
    ),
    as.character(
      "Country_numeric" %in%
        names(meta_final)
    ),
    "GENE-COLLAPSED, outcome-independent SOMAmer selection"
  )
)

safe_write_csv(
  input_audit,
  file.path(
    OUTDIR,
    "tables",
    "input",
    "script15b_input_alignment_audit.csv"
  )
)

cat("Samples:", nrow(datExpr_fixed), "\n")
cat("Full-network genes:", ncol(datExpr_clean), "\n")
cat("Fixed-network genes:", ncol(datExpr_fixed), "\n")
cat(
  "Fixed module sizes:",
  paste0(
    names(fixed_module_sizes),
    "=",
    as.integer(fixed_module_sizes),
    collapse = ", "
  ),
  "\n"
)
cat(
  "Restricted modules:",
  paste(
    fixed_module_summary$Module[
      fixed_module_summary$
        Module_was_restricted
    ],
    collapse = ", "
  ),
  "\n"
)
cat("Countries:", paste(countries, collapse = ", "), "\n")
cat("Reciprocal site runs:", nrow(site_pair_plan), "\n")
cat("Permutations per run:", N_PERMUTATIONS, "\n")
cat("Common preservation seed:", MODULE_PRESERVATION_SEED, "\n\n")

###############################################################################
# 11) COUNTRY-HELD-OUT FIXED-GENE PRESERVATION
###############################################################################

country_mp_objects <- list()
country_results <- list()
country_run_log <- list()

for (ii in seq_along(countries)) {
  held_out_country <- countries[[ii]]

  reference_index <- as.character(
    meta_final$Country
  ) != held_out_country

  test_index <- as.character(
    meta_final$Country
  ) == held_out_country

  reference_expr <- datExpr_fixed[
    reference_index,
    ,
    drop = FALSE
  ]

  test_expr <- datExpr_fixed[
    test_index,
    ,
    drop = FALSE
  ]

  run_id <- paste0(
    "COUNTRY_REST_TO_",
    sanitize_file_tag(
      held_out_country
    )
  )

  reference_label <- paste0(
    "All countries except ",
    held_out_country
  )

  test_label <- held_out_country

  run_result <- run_fixed_module_preservation_pair(
    reference_expr = reference_expr,
    test_expr = test_expr,
    fixed_module_colors =
      fixed_module_colors,
    run_id = run_id,
    analysis_type =
      "country_held_out_fixed_geneset",
    reference_label =
      reference_label,
    test_label = test_label
  )

  country_mp_objects[[
    held_out_country
  ]] <- run_result$object

  country_results[[
    held_out_country
  ]] <- run_result$table %>%
    dplyr::mutate(
      Held_out_country =
        held_out_country
    )

  country_run_log[[
    held_out_country
  ]] <- tibble::tibble(
    Analysis_type =
      "country_held_out_fixed_geneset",
    Run_ID = run_id,
    Reference_label =
      reference_label,
    Test_label = test_label,
    Held_out_country =
      held_out_country,
    N_reference =
      nrow(reference_expr),
    N_test =
      nrow(test_expr),
    N_fixed_genes =
      ncol(reference_expr),
    Fixed_geneset_seed =
      FIXED_GENESET_SEED,
    Module_preservation_seed =
      MODULE_PRESERVATION_SEED,
    Checkpoint_version =
      CHECKPOINT_VERSION_TAG,
    Resumed_existing_run =
      run_result$resumed,
    Raw_object_file =
      run_result$raw_object_file,
    Permutation_file =
      run_result$permutation_file
  )

  rm(
    reference_expr,
    test_expr,
    run_result
  )

  gc(verbose = FALSE)
}

country_preservation_raw <- dplyr::bind_rows(
  country_results
)

country_run_log_tbl <- dplyr::bind_rows(
  country_run_log
)

###############################################################################
# 12) CLEAN AND SUMMARIZE COUNTRY RESULTS
###############################################################################

country_preservation <- country_preservation_raw %>%
  dplyr::filter(
    !Module %in%
      c(
        GREY_NAME,
        GOLD_NAME
      )
  ) %>%
  dplyr::mutate(
    Module = as.character(Module)
  ) %>%
  dplyr::left_join(
    fixed_module_summary,
    by = "Module"
  )

required_preservation_columns <- c(
  "Zsummary.pres",
  "Zdensity.pres",
  "Zconnectivity.pres",
  "observed.medianRank.pres"
)

missing_preservation_columns <- setdiff(
  required_preservation_columns,
  names(country_preservation)
)

if (length(missing_preservation_columns) > 0) {
  stop(
    "Expected WGCNA preservation columns were not found: ",
    paste(
      missing_preservation_columns,
      collapse = ", "
    ),
    call. = FALSE
  )
}

country_preservation <- country_preservation %>%
  dplyr::mutate(
    Preservation_class =
      classify_zsummary(
        .data[["Zsummary.pres"]]
      ),
    Strong_preservation =
      !is.na(
        .data[["Zsummary.pres"]]
      ) &
      .data[["Zsummary.pres"]] >=
        10,
    At_least_moderate_preservation =
      !is.na(
        .data[["Zsummary.pres"]]
      ) &
      .data[["Zsummary.pres"]] >=
        2,
    Preservation_rank_within_run =
      rank(
        .data[[
          "observed.medianRank.pres"
        ]],
        ties.method = "average",
        na.last = "keep"
      ),
    Fixed_gene_set_identical_across_runs =
      TRUE
  ) %>%
  dplyr::arrange(
    Held_out_country,
    dplyr::desc(
      .data[["Zsummary.pres"]]
    )
  )

country_summary_by_module <- country_preservation %>%
  dplyr::group_by(Module) %>%
  dplyr::summarise(
    Original_module_size =
      dplyr::first(
        Original_module_size
      ),
    Fixed_module_size =
      dplyr::first(
        Fixed_module_size
      ),
    Module_was_restricted =
      dplyr::first(
        Module_was_restricted
      ),
    N_country_tests =
      dplyr::n(),
    Mean_Zsummary =
      safe_mean(
        .data[["Zsummary.pres"]]
      ),
    Median_Zsummary =
      safe_median(
        .data[["Zsummary.pres"]]
      ),
    Minimum_Zsummary =
      safe_min(
        .data[["Zsummary.pres"]]
      ),
    Maximum_Zsummary =
      safe_max(
        .data[["Zsummary.pres"]]
      ),
    Mean_Zdensity =
      safe_mean(
        .data[["Zdensity.pres"]]
      ),
    Minimum_Zdensity =
      safe_min(
        .data[["Zdensity.pres"]]
      ),
    Mean_Zconnectivity =
      safe_mean(
        .data[["Zconnectivity.pres"]]
      ),
    Minimum_Zconnectivity =
      safe_min(
        .data[["Zconnectivity.pres"]]
      ),
    Mean_medianRank =
      safe_mean(
        .data[[
          "observed.medianRank.pres"
        ]]
      ),
    Worst_medianRank =
      safe_max(
        .data[[
          "observed.medianRank.pres"
        ]]
      ),
    N_strong =
      sum(
        Strong_preservation,
        na.rm = TRUE
      ),
    N_at_least_moderate =
      sum(
        At_least_moderate_preservation,
        na.rm = TRUE
      ),
    Strong_preservation_rate =
      safe_proportion(
        Strong_preservation
      ),
    At_least_moderate_rate =
      safe_proportion(
        At_least_moderate_preservation
      ),
    .groups = "drop"
  )

country_worst_context <- country_preservation %>%
  dplyr::group_by(Module) %>%
  dplyr::slice_min(
    order_by =
      .data[["Zsummary.pres"]],
    n = 1,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(
    Module,
    Worst_country =
      Held_out_country,
    Worst_country_Zsummary =
      .data[["Zsummary.pres"]],
    Worst_country_Zdensity =
      .data[["Zdensity.pres"]],
    Worst_country_Zconnectivity =
      .data[[
        "Zconnectivity.pres"
      ]],
    Worst_country_medianRank =
      .data[[
        "observed.medianRank.pres"
      ]],
    Worst_country_class =
      Preservation_class
  )

country_summary_by_module <- country_summary_by_module %>%
  dplyr::left_join(
    country_worst_context,
    by = "Module"
  ) %>%
  dplyr::left_join(
    module_reference %>%
      dplyr::select(
        Module,
        Module_color,
        Biological_label
      ),
    by = "Module"
  ) %>%
  dplyr::arrange(
    dplyr::desc(
      Minimum_Zsummary
    ),
    Mean_medianRank
  )

country_summary_by_context <- country_preservation %>%
  dplyr::group_by(
    Held_out_country
  ) %>%
  dplyr::summarise(
    N_modules =
      dplyr::n(),
    Mean_Zsummary =
      safe_mean(
        .data[["Zsummary.pres"]]
      ),
    Median_Zsummary =
      safe_median(
        .data[["Zsummary.pres"]]
      ),
    Minimum_Zsummary =
      safe_min(
        .data[["Zsummary.pres"]]
      ),
    Mean_medianRank =
      safe_mean(
        .data[[
          "observed.medianRank.pres"
        ]]
      ),
    N_strong =
      sum(
        Strong_preservation,
        na.rm = TRUE
      ),
    N_at_least_moderate =
      sum(
        At_least_moderate_preservation,
        na.rm = TRUE
      ),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    Mean_Zsummary
  )

safe_write_csv(
  country_preservation_raw,
  file.path(
    OUTDIR,
    "tables",
    "country",
    "country_fixed_geneset_all_WGCNA_statistics_raw.csv"
  )
)

safe_write_csv(
  country_preservation,
  file.path(
    OUTDIR,
    "tables",
    "country",
    "country_fixed_geneset_primary_statistics.csv"
  )
)

safe_write_csv(
  country_summary_by_module,
  file.path(
    OUTDIR,
    "tables",
    "country",
    "country_fixed_geneset_summary_by_module.csv"
  )
)

safe_write_csv(
  country_summary_by_context,
  file.path(
    OUTDIR,
    "tables",
    "country",
    "country_fixed_geneset_summary_by_held_out_country.csv"
  )
)

safe_write_csv(
  country_run_log_tbl,
  file.path(
    OUTDIR,
    "tables",
    "country",
    "country_fixed_geneset_run_log.csv"
  )
)

###############################################################################
# 13) RECIPROCAL WITHIN-COUNTRY SITE PRESERVATION
###############################################################################

site_mp_objects <- list()
site_results <- list()
site_run_log <- list()

for (ii in seq_len(
  nrow(site_pair_plan)
)) {
  country_i <-
    site_pair_plan$Country[[ii]]

  reference_site_i <-
    site_pair_plan$Reference_site[[ii]]

  test_site_i <-
    site_pair_plan$Test_site[[ii]]

  run_id_i <-
    site_pair_plan$Run_ID[[ii]]

  reference_label_i <-
    site_pair_plan$Reference_label[[ii]]

  test_label_i <-
    site_pair_plan$Test_label[[ii]]

  reference_index <- (
    as.character(
      meta_final$Country
    ) == country_i &
      as.character(
        meta_final[[site_var]]
      ) == reference_site_i
  )

  test_index <- (
    as.character(
      meta_final$Country
    ) == country_i &
      as.character(
        meta_final[[site_var]]
      ) == test_site_i
  )

  reference_expr <- datExpr_fixed[
    reference_index,
    ,
    drop = FALSE
  ]

  test_expr <- datExpr_fixed[
    test_index,
    ,
    drop = FALSE
  ]

  run_result <- run_fixed_module_preservation_pair(
    reference_expr =
      reference_expr,
    test_expr =
      test_expr,
    fixed_module_colors =
      fixed_module_colors,
    run_id = run_id_i,
    analysis_type =
      "within_country_site_reciprocal_fixed_geneset",
    reference_label =
      reference_label_i,
    test_label =
      test_label_i
  )

  site_mp_objects[[
    run_id_i
  ]] <- run_result$object

  site_results[[
    run_id_i
  ]] <- run_result$table %>%
    dplyr::mutate(
      Country = country_i,
      Reference_site =
        reference_site_i,
      Test_site =
        test_site_i
    )

  site_run_log[[
    run_id_i
  ]] <- tibble::tibble(
    Analysis_type =
      "within_country_site_reciprocal_fixed_geneset",
    Run_ID = run_id_i,
    Country = country_i,
    Reference_site =
      reference_site_i,
    Test_site = test_site_i,
    Reference_label =
      reference_label_i,
    Test_label =
      test_label_i,
    N_reference =
      nrow(reference_expr),
    N_test =
      nrow(test_expr),
    N_fixed_genes =
      ncol(reference_expr),
    Fixed_geneset_seed =
      FIXED_GENESET_SEED,
    Module_preservation_seed =
      MODULE_PRESERVATION_SEED,
    Checkpoint_version =
      CHECKPOINT_VERSION_TAG,
    Resumed_existing_run =
      run_result$resumed,
    Raw_object_file =
      run_result$raw_object_file,
    Permutation_file =
      run_result$permutation_file
  )

  rm(
    reference_expr,
    test_expr,
    run_result
  )

  gc(verbose = FALSE)
}

site_preservation_raw <- dplyr::bind_rows(
  site_results
)

site_run_log_tbl <- dplyr::bind_rows(
  site_run_log
)

###############################################################################
# 14) CLEAN AND SUMMARIZE SITE RESULTS
###############################################################################

site_preservation <- site_preservation_raw %>%
  dplyr::filter(
    !Module %in%
      c(
        GREY_NAME,
        GOLD_NAME
      )
  ) %>%
  dplyr::mutate(
    Module = as.character(Module)
  ) %>%
  dplyr::left_join(
    fixed_module_summary,
    by = "Module"
  )

missing_site_columns <- setdiff(
  required_preservation_columns,
  names(site_preservation)
)

if (length(missing_site_columns) > 0) {
  stop(
    "Expected site-preservation columns were not found: ",
    paste(
      missing_site_columns,
      collapse = ", "
    ),
    call. = FALSE
  )
}

site_preservation <- site_preservation %>%
  dplyr::mutate(
    Preservation_class =
      classify_zsummary(
        .data[["Zsummary.pres"]]
      ),
    Strong_preservation =
      !is.na(
        .data[["Zsummary.pres"]]
      ) &
      .data[["Zsummary.pres"]] >=
        10,
    At_least_moderate_preservation =
      !is.na(
        .data[["Zsummary.pres"]]
      ) &
      .data[["Zsummary.pres"]] >=
        2,
    Preservation_rank_within_run =
      rank(
        .data[[
          "observed.medianRank.pres"
        ]],
        ties.method = "average",
        na.last = "keep"
      ),
    Fixed_gene_set_identical_across_runs =
      TRUE,
    Site_direction = paste(
      Reference_site,
      "->",
      Test_site
    )
  ) %>%
  dplyr::arrange(
    Country,
    Reference_site,
    Test_site,
    dplyr::desc(
      .data[["Zsummary.pres"]]
    )
  )

site_summary_by_module <- site_preservation %>%
  dplyr::group_by(Module) %>%
  dplyr::summarise(
    Original_module_size =
      dplyr::first(
        Original_module_size
      ),
    Fixed_module_size =
      dplyr::first(
        Fixed_module_size
      ),
    Module_was_restricted =
      dplyr::first(
        Module_was_restricted
      ),
    N_directional_site_tests =
      dplyr::n(),
    Mean_Zsummary =
      safe_mean(
        .data[["Zsummary.pres"]]
      ),
    Median_Zsummary =
      safe_median(
        .data[["Zsummary.pres"]]
      ),
    Minimum_Zsummary =
      safe_min(
        .data[["Zsummary.pres"]]
      ),
    Maximum_Zsummary =
      safe_max(
        .data[["Zsummary.pres"]]
      ),
    Mean_Zdensity =
      safe_mean(
        .data[["Zdensity.pres"]]
      ),
    Minimum_Zdensity =
      safe_min(
        .data[["Zdensity.pres"]]
      ),
    Mean_Zconnectivity =
      safe_mean(
        .data[["Zconnectivity.pres"]]
      ),
    Minimum_Zconnectivity =
      safe_min(
        .data[["Zconnectivity.pres"]]
      ),
    Mean_medianRank =
      safe_mean(
        .data[[
          "observed.medianRank.pres"
        ]]
      ),
    Worst_medianRank =
      safe_max(
        .data[[
          "observed.medianRank.pres"
        ]]
      ),
    N_strong =
      sum(
        Strong_preservation,
        na.rm = TRUE
      ),
    N_at_least_moderate =
      sum(
        At_least_moderate_preservation,
        na.rm = TRUE
      ),
    Strong_preservation_rate =
      safe_proportion(
        Strong_preservation
      ),
    At_least_moderate_rate =
      safe_proportion(
        At_least_moderate_preservation
      ),
    .groups = "drop"
  )

site_worst_context <- site_preservation %>%
  dplyr::group_by(Module) %>%
  dplyr::slice_min(
    order_by =
      .data[["Zsummary.pres"]],
    n = 1,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(
    Module,
    Worst_site_country =
      Country,
    Worst_site_reference =
      Reference_site,
    Worst_site_test =
      Test_site,
    Worst_site_Zsummary =
      .data[["Zsummary.pres"]],
    Worst_site_Zdensity =
      .data[["Zdensity.pres"]],
    Worst_site_Zconnectivity =
      .data[[
        "Zconnectivity.pres"
      ]],
    Worst_site_medianRank =
      .data[[
        "observed.medianRank.pres"
      ]],
    Worst_site_class =
      Preservation_class
  )

site_summary_by_module <- site_summary_by_module %>%
  dplyr::left_join(
    site_worst_context,
    by = "Module"
  ) %>%
  dplyr::left_join(
    module_reference %>%
      dplyr::select(
        Module,
        Module_color,
        Biological_label
      ),
    by = "Module"
  ) %>%
  dplyr::arrange(
    dplyr::desc(
      Minimum_Zsummary
    ),
    Mean_medianRank
  )

site_reciprocal_summary <- site_preservation %>%
  dplyr::group_by(
    Country,
    Module
  ) %>%
  dplyr::summarise(
    N_directions =
      dplyr::n(),
    Mean_Zsummary =
      safe_mean(
        .data[["Zsummary.pres"]]
      ),
    Minimum_Zsummary =
      safe_min(
        .data[["Zsummary.pres"]]
      ),
    Maximum_Zsummary =
      safe_max(
        .data[["Zsummary.pres"]]
      ),
    Mean_Zdensity =
      safe_mean(
        .data[["Zdensity.pres"]]
      ),
    Mean_Zconnectivity =
      safe_mean(
        .data[["Zconnectivity.pres"]]
      ),
    Mean_medianRank =
      safe_mean(
        .data[[
          "observed.medianRank.pres"
        ]]
      ),
    Both_directions_strong =
      all(
        Strong_preservation,
        na.rm = TRUE
      ),
    Both_directions_at_least_moderate =
      all(
        At_least_moderate_preservation,
        na.rm = TRUE
      ),
    Directional_classes =
      paste(
        Reference_site,
        "->",
        Test_site,
        ": ",
        Preservation_class,
        collapse = "; "
      ),
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    module_reference,
    by = "Module"
  ) %>%
  dplyr::arrange(
    Country,
    dplyr::desc(
      Minimum_Zsummary
    )
  )

safe_write_csv(
  site_preservation_raw,
  file.path(
    OUTDIR,
    "tables",
    "site",
    "site_fixed_geneset_all_WGCNA_statistics_raw.csv"
  )
)

safe_write_csv(
  site_preservation,
  file.path(
    OUTDIR,
    "tables",
    "site",
    "site_fixed_geneset_primary_statistics.csv"
  )
)

safe_write_csv(
  site_summary_by_module,
  file.path(
    OUTDIR,
    "tables",
    "site",
    "site_fixed_geneset_summary_by_module.csv"
  )
)

safe_write_csv(
  site_reciprocal_summary,
  file.path(
    OUTDIR,
    "tables",
    "site",
    "site_fixed_geneset_reciprocal_summary_by_country_module.csv"
  )
)

safe_write_csv(
  site_run_log_tbl,
  file.path(
    OUTDIR,
    "tables",
    "site",
    "site_fixed_geneset_run_log.csv"
  )
)

###############################################################################
# 15) OPTIONAL COMPARISON WITH SCRIPT 15
###############################################################################

old_country_available <- file.exists(
  OLD_COUNTRY_FILE
)

old_site_available <- file.exists(
  OLD_SITE_FILE
)

country_old_vs_fixed <- tibble::tibble()
site_old_vs_fixed <- tibble::tibble()
comparison_summary <- tibble::tibble()

if (old_country_available) {
  old_country <- safe_read_csv(
    OLD_COUNTRY_FILE
  ) %>%
    dplyr::transmute(
      Module = as.character(Module),
      Held_out_country =
        as.character(
          Held_out_country
        ),
      Old_Zsummary =
        safe_numeric(
          .data[["Zsummary.pres"]]
        ),
      Old_Zdensity =
        safe_numeric(
          .data[["Zdensity.pres"]]
        ),
      Old_Zconnectivity =
        safe_numeric(
          .data[[
            "Zconnectivity.pres"
          ]]
        ),
      Old_medianRank =
        safe_numeric(
          .data[[
            "observed.medianRank.pres"
          ]]
        ),
      Old_class =
        classify_zsummary(
          Old_Zsummary
        )
    )

  country_old_vs_fixed <- country_preservation %>%
    dplyr::transmute(
      Module,
      Held_out_country =
        as.character(
          Held_out_country
        ),
      Fixed_Zsummary =
        safe_numeric(
          .data[["Zsummary.pres"]]
        ),
      Fixed_Zdensity =
        safe_numeric(
          .data[["Zdensity.pres"]]
        ),
      Fixed_Zconnectivity =
        safe_numeric(
          .data[[
            "Zconnectivity.pres"
          ]]
        ),
      Fixed_medianRank =
        safe_numeric(
          .data[[
            "observed.medianRank.pres"
          ]]
        ),
      Fixed_class =
        Preservation_class,
      Original_module_size,
      Fixed_module_size,
      Module_was_restricted
    ) %>%
    dplyr::left_join(
      old_country,
      by = c(
        "Module",
        "Held_out_country"
      )
    ) %>%
    dplyr::mutate(
      Delta_Zsummary =
        Fixed_Zsummary -
        Old_Zsummary,
      Abs_delta_Zsummary =
        abs(
          Delta_Zsummary
        ),
      Delta_Zdensity =
        Fixed_Zdensity -
        Old_Zdensity,
      Delta_Zconnectivity =
        Fixed_Zconnectivity -
        Old_Zconnectivity,
      Delta_medianRank =
        Fixed_medianRank -
        Old_medianRank,
      Preservation_class_unchanged =
        Fixed_class ==
        Old_class,
      Crossed_Zsummary_2_threshold =
        (
          Old_Zsummary < 2 &
            Fixed_Zsummary >= 2
        ) |
        (
          Old_Zsummary >= 2 &
            Fixed_Zsummary < 2
        ),
      Crossed_Zsummary_10_threshold =
        (
          Old_Zsummary < 10 &
            Fixed_Zsummary >= 10
        ) |
        (
          Old_Zsummary >= 10 &
            Fixed_Zsummary < 10
        )
    ) %>%
    dplyr::left_join(
      module_reference %>%
        dplyr::select(
          Module,
          Module_color,
          Biological_label
        ),
      by = "Module"
    )

  safe_write_csv(
    country_old_vs_fixed,
    file.path(
      OUTDIR,
      "tables",
      "comparison",
      "script15_vs_script15b_country_comparison.csv"
    )
  )
}

if (old_site_available) {
  old_site <- safe_read_csv(
    OLD_SITE_FILE
  ) %>%
    dplyr::transmute(
      Module = as.character(Module),
      Country = as.character(Country),
      Reference_site =
        as.character(
          Reference_site
        ),
      Test_site =
        as.character(
          Test_site
        ),
      Old_Zsummary =
        safe_numeric(
          .data[["Zsummary.pres"]]
        ),
      Old_Zdensity =
        safe_numeric(
          .data[["Zdensity.pres"]]
        ),
      Old_Zconnectivity =
        safe_numeric(
          .data[[
            "Zconnectivity.pres"
          ]]
        ),
      Old_medianRank =
        safe_numeric(
          .data[[
            "observed.medianRank.pres"
          ]]
        ),
      Old_class =
        classify_zsummary(
          Old_Zsummary
        )
    )

  site_old_vs_fixed <- site_preservation %>%
    dplyr::transmute(
      Module,
      Country =
        as.character(Country),
      Reference_site =
        as.character(
          Reference_site
        ),
      Test_site =
        as.character(
          Test_site
        ),
      Fixed_Zsummary =
        safe_numeric(
          .data[["Zsummary.pres"]]
        ),
      Fixed_Zdensity =
        safe_numeric(
          .data[["Zdensity.pres"]]
        ),
      Fixed_Zconnectivity =
        safe_numeric(
          .data[[
            "Zconnectivity.pres"
          ]]
        ),
      Fixed_medianRank =
        safe_numeric(
          .data[[
            "observed.medianRank.pres"
          ]]
        ),
      Fixed_class =
        Preservation_class,
      Original_module_size,
      Fixed_module_size,
      Module_was_restricted
    ) %>%
    dplyr::left_join(
      old_site,
      by = c(
        "Module",
        "Country",
        "Reference_site",
        "Test_site"
      )
    ) %>%
    dplyr::mutate(
      Delta_Zsummary =
        Fixed_Zsummary -
        Old_Zsummary,
      Abs_delta_Zsummary =
        abs(
          Delta_Zsummary
        ),
      Delta_Zdensity =
        Fixed_Zdensity -
        Old_Zdensity,
      Delta_Zconnectivity =
        Fixed_Zconnectivity -
        Old_Zconnectivity,
      Delta_medianRank =
        Fixed_medianRank -
        Old_medianRank,
      Preservation_class_unchanged =
        Fixed_class ==
        Old_class,
      Crossed_Zsummary_2_threshold =
        (
          Old_Zsummary < 2 &
            Fixed_Zsummary >= 2
        ) |
        (
          Old_Zsummary >= 2 &
            Fixed_Zsummary < 2
        ),
      Crossed_Zsummary_10_threshold =
        (
          Old_Zsummary < 10 &
            Fixed_Zsummary >= 10
        ) |
        (
          Old_Zsummary >= 10 &
            Fixed_Zsummary < 10
        )
    ) %>%
    dplyr::left_join(
      module_reference %>%
        dplyr::select(
          Module,
          Module_color,
          Biological_label
        ),
      by = "Module"
    )

  safe_write_csv(
    site_old_vs_fixed,
    file.path(
      OUTDIR,
      "tables",
      "comparison",
      "script15_vs_script15b_site_comparison.csv"
    )
  )
}

comparison_summary <- dplyr::bind_rows(
  if (nrow(country_old_vs_fixed) > 0) {
    country_old_vs_fixed %>%
      dplyr::mutate(
        Analysis = "Country"
      )
  } else {
    tibble::tibble()
  },
  if (nrow(site_old_vs_fixed) > 0) {
    site_old_vs_fixed %>%
      dplyr::mutate(
        Analysis = "Site"
      )
  } else {
    tibble::tibble()
  }
) %>%
  dplyr::group_by(
    Analysis,
    Module
  ) %>%
  dplyr::summarise(
    N_comparisons =
      dplyr::n(),
    Original_module_size =
      dplyr::first(
        Original_module_size
      ),
    Fixed_module_size =
      dplyr::first(
        Fixed_module_size
      ),
    Module_was_restricted =
      dplyr::first(
        Module_was_restricted
      ),
    Mean_delta_Zsummary =
      safe_mean(
        Delta_Zsummary
      ),
    Mean_abs_delta_Zsummary =
      safe_mean(
        Abs_delta_Zsummary
      ),
    Maximum_abs_delta_Zsummary =
      safe_max(
        Abs_delta_Zsummary
      ),
    Preservation_class_agreement_rate =
      safe_proportion(
        Preservation_class_unchanged
      ),
    N_crossed_Zsummary_2 =
      sum(
        Crossed_Zsummary_2_threshold,
        na.rm = TRUE
      ),
    N_crossed_Zsummary_10 =
      sum(
        Crossed_Zsummary_10_threshold,
        na.rm = TRUE
      ),
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    module_reference,
    by = "Module"
  )

if (nrow(comparison_summary) > 0) {
  safe_write_csv(
    comparison_summary,
    file.path(
      OUTDIR,
      "tables",
      "comparison",
      "script15_vs_script15b_comparison_summary_by_module.csv"
    )
  )
}

###############################################################################
# 16) FIGURES — CORRECTED COUNTRY PRESERVATION
###############################################################################

country_plot_df <- country_preservation %>%
  dplyr::mutate(
    Held_out_country = factor(
      Held_out_country,
      levels = countries
    ),
    Module = factor(
      Module,
      levels = rev(
        modules_use
      )
    )
  )

country_fill_limits <- range(
  country_plot_df[[
    "Zsummary.pres"
  ]],
  na.rm = TRUE
)

country_fill_breaks <- unique(
  sort(
    c(
      country_fill_limits[[1]],
      2,
      5,
      10,
      country_fill_limits[[2]]
    )
  )
)

country_fill_values <- scales::rescale(
  country_fill_breaks,
  from = country_fill_limits
)

country_fill_colors <- grDevices::colorRampPalette(
  c(
    "#B2182B",
    "#FDDDBC",
    "#FFF7EC",
    "#74ADD1",
    "#2166AC"
  )
)(
  length(country_fill_breaks)
)

p_country <- ggplot(
  country_plot_df,
  aes(
    x = Held_out_country,
    y = Module,
    fill =
      .data[["Zsummary.pres"]]
  )
) +
  geom_tile(
    colour = "white",
    linewidth = 0.5
  ) +
  geom_text(
    aes(
      label = sprintf(
        "%.1f",
        .data[["Zsummary.pres"]]
      )
    ),
    size = 3.2
  ) +
  scale_fill_gradientn(
    colours =
      country_fill_colors,
    values =
      country_fill_values,
    name = "Zsummary",
    oob = scales::squish
  ) +
  labs(
    title =
      "Country-held-out structural preservation using one fixed gene set",
    subtitle = paste0(
      "The same ",
      ncol(datExpr_fixed),
      " genes and common permutation seed were used in all country runs"
    ),
    x = "Held-out test country",
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
    plot.title = element_text(
      face = "bold"
    )
  )

ggsave(
  file.path(
    OUTDIR,
    "figures",
    "country",
    "country_fixed_geneset_Zsummary_heatmap.pdf"
  ),
  p_country,
  width = 8,
  height = 5.8
)

ggsave(
  file.path(
    OUTDIR,
    "figures",
    "country",
    "country_fixed_geneset_Zsummary_heatmap.png"
  ),
  p_country,
  width = 8,
  height = 5.8,
  dpi = DPI
)

###############################################################################
# 17) FIGURES — CORRECTED SITE PRESERVATION
###############################################################################

site_plot_df <- site_preservation %>%
  dplyr::mutate(
    Module = factor(
      Module,
      levels = rev(
        modules_use
      )
    )
  )

site_fill_limits <- range(
  site_plot_df[[
    "Zsummary.pres"
  ]],
  na.rm = TRUE
)

site_fill_breaks <- unique(
  sort(
    c(
      site_fill_limits[[1]],
      2,
      5,
      10,
      site_fill_limits[[2]]
    )
  )
)

site_fill_values <- scales::rescale(
  site_fill_breaks,
  from = site_fill_limits
)

site_fill_colors <- grDevices::colorRampPalette(
  c(
    "#B2182B",
    "#FDDDBC",
    "#FFF7EC",
    "#74ADD1",
    "#2166AC"
  )
)(
  length(site_fill_breaks)
)

p_site <- ggplot(
  site_plot_df,
  aes(
    x = Site_direction,
    y = Module,
    fill =
      .data[["Zsummary.pres"]]
  )
) +
  geom_tile(
    colour = "white",
    linewidth = 0.5
  ) +
  geom_text(
    aes(
      label = sprintf(
        "%.1f",
        .data[["Zsummary.pres"]]
      )
    ),
    size = 3.1
  ) +
  facet_grid(
    . ~ Country,
    scales = "free_x",
    space = "free_x"
  ) +
  scale_fill_gradientn(
    colours =
      site_fill_colors,
    values =
      site_fill_values,
    name = "Zsummary",
    oob = scales::squish
  ) +
  labs(
    title =
      "Reciprocal site preservation using one fixed gene set",
    subtitle = paste0(
      "Each arrow denotes reference site -> test site; ",
      ncol(datExpr_fixed),
      " identical genes were used in all runs"
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
      angle = 35,
      hjust = 1
    ),
    axis.text.y = element_text(
      face = "bold"
    ),
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
    "site",
    "site_fixed_geneset_Zsummary_heatmap.pdf"
  ),
  p_site,
  width = 9,
  height = 6
)

ggsave(
  file.path(
    OUTDIR,
    "figures",
    "site",
    "site_fixed_geneset_Zsummary_heatmap.png"
  ),
  p_site,
  width = 9,
  height = 6,
  dpi = DPI
)

###############################################################################
# 18) FIGURES — SCRIPT 15 VERSUS SCRIPT 15b
###############################################################################

if (nrow(country_old_vs_fixed) > 0) {
  p_country_compare <- ggplot(
    country_old_vs_fixed,
    aes(
      x = Old_Zsummary,
      y = Fixed_Zsummary,
      colour = Module,
      shape =
        Module_was_restricted
    )
  ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = 2,
      colour = "grey50"
    ) +
    geom_hline(
      yintercept = 2,
      linetype = 3,
      colour = "grey65"
    ) +
    geom_hline(
      yintercept = 10,
      linetype = 3,
      colour = "grey65"
    ) +
    geom_vline(
      xintercept = 2,
      linetype = 3,
      colour = "grey65"
    ) +
    geom_vline(
      xintercept = 10,
      linetype = 3,
      colour = "grey65"
    ) +
    geom_point(
      size = 2.8,
      alpha = 0.85
    ) +
    facet_wrap(
      ~ Held_out_country
    ) +
    scale_colour_manual(
      values = get_module_colors(
        modules_use
      )
    ) +
    scale_shape_manual(
      values = c(
        `FALSE` = 16,
        `TRUE` = 17
      ),
      labels = c(
        `FALSE` =
          "All genes retained",
        `TRUE` =
          "Fixed 1,000-gene subset"
      ),
      name = "Module handling"
    ) +
    labs(
      title =
        "Country Zsummary before and after fixed-gene correction",
      subtitle =
        "The diagonal indicates identical Zsummary values",
      x =
        "Script 15 Zsummary",
      y =
        "Script 15b fixed-gene Zsummary",
      colour = "Module"
    ) +
    theme_bw(
      base_size = 11
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
      "comparison",
      "script15_vs_script15b_country_Zsummary.pdf"
    ),
    p_country_compare,
    width = 10,
    height = 6.5
  )

  ggsave(
    file.path(
      OUTDIR,
      "figures",
      "comparison",
      "script15_vs_script15b_country_Zsummary.png"
    ),
    p_country_compare,
    width = 10,
    height = 6.5,
    dpi = DPI
  )
}

if (nrow(site_old_vs_fixed) > 0) {
  site_old_vs_fixed <- site_old_vs_fixed %>%
    dplyr::mutate(
      Site_direction = paste(
        Reference_site,
        "->",
        Test_site
      )
    )

  p_site_compare <- ggplot(
    site_old_vs_fixed,
    aes(
      x = Old_Zsummary,
      y = Fixed_Zsummary,
      colour = Module,
      shape =
        Module_was_restricted
    )
  ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = 2,
      colour = "grey50"
    ) +
    geom_hline(
      yintercept = 2,
      linetype = 3,
      colour = "grey65"
    ) +
    geom_hline(
      yintercept = 10,
      linetype = 3,
      colour = "grey65"
    ) +
    geom_vline(
      xintercept = 2,
      linetype = 3,
      colour = "grey65"
    ) +
    geom_vline(
      xintercept = 10,
      linetype = 3,
      colour = "grey65"
    ) +
    geom_point(
      size = 2.8,
      alpha = 0.85
    ) +
    facet_grid(
      . ~ Country,
      scales = "free"
    ) +
    scale_colour_manual(
      values = get_module_colors(
        modules_use
      )
    ) +
    scale_shape_manual(
      values = c(
        `FALSE` = 16,
        `TRUE` = 17
      ),
      labels = c(
        `FALSE` =
          "All genes retained",
        `TRUE` =
          "Fixed 1,000-gene subset"
      ),
      name = "Module handling"
    ) +
    labs(
      title =
        "Site Zsummary before and after fixed-gene correction",
      subtitle =
        "The diagonal indicates identical Zsummary values",
      x =
        "Script 15 Zsummary",
      y =
        "Script 15b fixed-gene Zsummary",
      colour = "Module"
    ) +
    theme_bw(
      base_size = 11
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
      "comparison",
      "script15_vs_script15b_site_Zsummary.pdf"
    ),
    p_site_compare,
    width = 9,
    height = 5.8
  )

  ggsave(
    file.path(
      OUTDIR,
      "figures",
      "comparison",
      "script15_vs_script15b_site_Zsummary.png"
    ),
    p_site_compare,
    width = 9,
    height = 5.8,
    dpi = DPI
  )
}

###############################################################################
# 19) DEFINITIVE INTEGRATED SUMMARY
###############################################################################

definitive_summary <- fixed_module_summary %>%
  dplyr::select(
    Module,
    Original_module_size,
    Fixed_module_size,
    Module_was_restricted,
    Selection_method,
    Module_color,
    Biological_label
  ) %>%
  dplyr::left_join(
    country_summary_by_module %>%
      dplyr::select(
        Module,
        Country_mean_Zsummary =
          Mean_Zsummary,
        Country_median_Zsummary =
          Median_Zsummary,
        Country_minimum_Zsummary =
          Minimum_Zsummary,
        Country_maximum_Zsummary =
          Maximum_Zsummary,
        Country_mean_Zdensity =
          Mean_Zdensity,
        Country_minimum_Zdensity =
          Minimum_Zdensity,
        Country_mean_Zconnectivity =
          Mean_Zconnectivity,
        Country_minimum_Zconnectivity =
          Minimum_Zconnectivity,
        Country_mean_medianRank =
          Mean_medianRank,
        Country_strong_rate =
          Strong_preservation_rate,
        Country_at_least_moderate_rate =
          At_least_moderate_rate,
        Worst_country,
        Worst_country_class
      ),
    by = "Module"
  ) %>%
  dplyr::left_join(
    site_summary_by_module %>%
      dplyr::select(
        Module,
        Site_mean_Zsummary =
          Mean_Zsummary,
        Site_median_Zsummary =
          Median_Zsummary,
        Site_minimum_Zsummary =
          Minimum_Zsummary,
        Site_maximum_Zsummary =
          Maximum_Zsummary,
        Site_mean_Zdensity =
          Mean_Zdensity,
        Site_minimum_Zdensity =
          Minimum_Zdensity,
        Site_mean_Zconnectivity =
          Mean_Zconnectivity,
        Site_minimum_Zconnectivity =
          Minimum_Zconnectivity,
        Site_mean_medianRank =
          Mean_medianRank,
        Site_strong_rate =
          Strong_preservation_rate,
        Site_at_least_moderate_rate =
          At_least_moderate_rate,
        Worst_site_country,
        Worst_site_reference,
        Worst_site_test,
        Worst_site_class
      ),
    by = "Module"
  ) %>%
  dplyr::mutate(
    Country_minimum_preservation_class =
      classify_zsummary(
        Country_minimum_Zsummary
      ),
    Site_minimum_preservation_class =
      classify_zsummary(
        Site_minimum_Zsummary
      ),
    Fixed_gene_comparability =
      "Identical genes used in every country and site run",
    Interpretation = paste(
      "Canonical WGCNA preservation using a single fixed module gene set;",
      "report with medianRank and component Z statistics."
    )
  ) %>%
  dplyr::arrange(
    dplyr::desc(
      Country_minimum_Zsummary
    ),
    dplyr::desc(
      Site_minimum_Zsummary
    )
  )

safe_write_csv(
  definitive_summary,
  file.path(
    OUTDIR,
    "tables",
    "fixed_geneset_structural_preservation_definitive_summary.csv"
  )
)

###############################################################################
# 20) METHODS-READY WORDING
###############################################################################

methods_wording <- tibble::tribble(
  ~section,
  ~text,

  "Methods - fixed module gene sets",
  paste(
    "To ensure direct comparability of preservation statistics across",
    "recruitment contexts, a single gene set was selected before any",
    "preservation analysis. All genes were retained for modules containing",
    "1,000 or fewer genes, whereas one reproducible random subset of 1,000",
    "genes was selected from each larger module. The same fixed genes were",
    "used in every country and site comparison."
  ),

  "Methods - country preservation",
  paste(
    "For each country, expression data from all remaining countries served",
    "as the reference network and the held-out country served as the test",
    "network. Reference and test participants were therefore non-overlapping.",
    "Module assignments were fixed to those obtained from the",
    "outcome-independent full-cohort WGCNA network."
  ),

  "Methods - site preservation",
  paste(
    "Recruitment-site preservation was evaluated reciprocally within Chile",
    "and Colombia, the two countries containing more than one recruitment",
    "site. Each site was used once as the reference and once as the test set."
  ),

  "Methods - preservation inference",
  paste(
    "WGCNA modulePreservation statistics were calculated using signed",
    "Pearson-correlation networks and 100 permutations. The same permutation",
    "seed was used in every context. Primary statistics included Zsummary,",
    "Zdensity, Zconnectivity and medianRank. Zsummary values of at least 10",
    "were described as strong evidence of preservation, values from 2 to",
    "below 10 as weak-to-moderate evidence and values below 2 as no evidence."
  ),

  "Methods - large-module correction",
  paste(
    "Because all module gene sets were restricted to at most 1,000 genes",
    "before calling modulePreservation, the function did not perform",
    "context-specific subsampling of large modules. Thus, blue, brown and",
    "magenta were evaluated using identical genes across all runs."
  ),

  "Interpretation",
  paste(
    "The analysis assesses internal structural reproducibility of fixed",
    "full-cohort modules. It does not constitute external validation because",
    "module assignments were learned using the complete cohort. Structural",
    "preservation was interpreted separately from module-trait association",
    "stability."
  )
)

safe_write_csv(
  methods_wording,
  file.path(
    OUTDIR,
    "tables",
    "script15b_methods_and_interpretation_wording.csv"
  )
)

###############################################################################
# 21) EXCEL WORKBOOK
###############################################################################

workbook_tables <- list(
  Input_audit =
    input_audit,
  Fixed_module_summary =
    fixed_module_summary,
  Country_primary =
    country_preservation,
  Country_by_module =
    country_summary_by_module,
  Country_by_context =
    country_summary_by_context,
  Site_primary =
    site_preservation,
  Site_by_module =
    site_summary_by_module,
  Site_reciprocal =
    site_reciprocal_summary,
  Definitive_summary =
    definitive_summary,
  Country_run_log =
    country_run_log_tbl,
  Site_run_log =
    site_run_log_tbl,
  Methods_wording =
    methods_wording
)

if (nrow(country_old_vs_fixed) > 0) {
  workbook_tables[["Old_vs_fixed_country"]] <-
    country_old_vs_fixed
}

if (nrow(site_old_vs_fixed) > 0) {
  workbook_tables[["Old_vs_fixed_site"]] <-
    site_old_vs_fixed
}

if (nrow(comparison_summary) > 0) {
  workbook_tables[["Comparison_summary"]] <-
    comparison_summary
}

openxlsx::write.xlsx(
  workbook_tables,
  file = file.path(
    OUTDIR,
    "WGCNA_15b_Fixed_Gene_Set_Structural_Preservation.xlsx"
  ),
  overwrite = TRUE
)

###############################################################################
# 22) FINAL SUMMARY AND MANIFEST
###############################################################################

script15b_summary <- tibble::tibble(
  metric = c(
    "base_dir",
    "output_dir",
    "wgcna_input_level",
    "n_samples",
    "n_full_network_genes",
    "n_fixed_network_genes",
    "n_modules",
    "modules",
    "fixed_module_sizes",
    "restricted_modules",
    "fixed_geneset_seed",
    "common_module_preservation_seed",
    "checkpoint_version",
    "n_permutations",
    "n_country_runs",
    "n_country_module_results",
    "n_reciprocal_site_runs",
    "n_site_module_results",
    "country_runs_resumed",
    "site_runs_resumed",
    "n_country_strong_results",
    "n_country_at_least_moderate_results",
    "n_site_strong_results",
    "n_site_at_least_moderate_results",
    "old_country_comparison_available",
    "old_site_comparison_available",
    "Country_numeric_included"
  ),
  value = c(
    BASE_DIR,
    OUTDIR,
    "GENE-COLLAPSED, outcome-independent SOMAmer selection",
    as.character(
      nrow(datExpr_fixed)
    ),
    as.character(
      ncol(datExpr_clean)
    ),
    as.character(
      ncol(datExpr_fixed)
    ),
    as.character(
      length(modules_use)
    ),
    paste(
      modules_use,
      collapse = ", "
    ),
    paste0(
      names(fixed_module_sizes),
      "=",
      as.integer(
        fixed_module_sizes
      ),
      collapse = "; "
    ),
    paste(
      fixed_module_summary$Module[
        fixed_module_summary$
          Module_was_restricted
      ],
      collapse = ", "
    ),
    as.character(
      FIXED_GENESET_SEED
    ),
    as.character(
      MODULE_PRESERVATION_SEED
    ),
    CHECKPOINT_VERSION_TAG,
    as.character(
      N_PERMUTATIONS
    ),
    as.character(
      length(countries)
    ),
    as.character(
      nrow(country_preservation)
    ),
    as.character(
      nrow(site_pair_plan)
    ),
    as.character(
      nrow(site_preservation)
    ),
    as.character(
      sum(
        country_run_log_tbl$
          Resumed_existing_run,
        na.rm = TRUE
      )
    ),
    as.character(
      sum(
        site_run_log_tbl$
          Resumed_existing_run,
        na.rm = TRUE
      )
    ),
    as.character(
      sum(
        country_preservation$
          Strong_preservation,
        na.rm = TRUE
      )
    ),
    as.character(
      sum(
        country_preservation$
          At_least_moderate_preservation,
        na.rm = TRUE
      )
    ),
    as.character(
      sum(
        site_preservation$
          Strong_preservation,
        na.rm = TRUE
      )
    ),
    as.character(
      sum(
        site_preservation$
          At_least_moderate_preservation,
        na.rm = TRUE
      )
    ),
    as.character(
      old_country_available
    ),
    as.character(
      old_site_available
    ),
    as.character(
      "Country_numeric" %in%
        names(meta_final)
    )
  )
)

safe_write_csv(
  script15b_summary,
  file.path(
    OUTDIR,
    "tables",
    "script15b_final_summary.csv"
  )
)

output_manifest <- tibble::tibble(
  output_file = c(
    "tables/input/script15b_input_alignment_audit.csv",
    "tables/input/fixed_gene_manifest_all_modules.csv",
    "tables/input/fixed_gene_manifest_summary_by_module.csv",
    "tables/input/country_site_sample_counts.csv",
    "tables/input/reciprocal_within_country_site_run_plan.csv",
    "tables/country/country_fixed_geneset_all_WGCNA_statistics_raw.csv",
    "tables/country/country_fixed_geneset_primary_statistics.csv",
    "tables/country/country_fixed_geneset_summary_by_module.csv",
    "tables/country/country_fixed_geneset_summary_by_held_out_country.csv",
    "tables/country/country_fixed_geneset_run_log.csv",
    "tables/site/site_fixed_geneset_all_WGCNA_statistics_raw.csv",
    "tables/site/site_fixed_geneset_primary_statistics.csv",
    "tables/site/site_fixed_geneset_summary_by_module.csv",
    "tables/site/site_fixed_geneset_reciprocal_summary_by_country_module.csv",
    "tables/site/site_fixed_geneset_run_log.csv",
    "tables/comparison/script15_vs_script15b_country_comparison.csv",
    "tables/comparison/script15_vs_script15b_site_comparison.csv",
    "tables/comparison/script15_vs_script15b_comparison_summary_by_module.csv",
    "tables/fixed_geneset_structural_preservation_definitive_summary.csv",
    "figures/country/country_fixed_geneset_Zsummary_heatmap.pdf/png",
    "figures/site/site_fixed_geneset_Zsummary_heatmap.pdf/png",
    "figures/comparison/script15_vs_script15b_country_Zsummary.pdf/png",
    "figures/comparison/script15_vs_script15b_site_Zsummary.pdf/png",
    "tables/script15b_methods_and_interpretation_wording.csv",
    "WGCNA_15b_Fixed_Gene_Set_Structural_Preservation.xlsx",
    "workspace/raw_module_preservation_objects/*.rds",
    "workspace/permutation_statistics/*.RData",
    "workspace/script15b_fixed_geneset_structural_preservation_workspace.RData",
    "sessionInfo.txt"
  ),
  description = c(
    "Strict sample, gene, module, country, site and seed audit.",
    "Exact gene-level manifest used identically in every preservation run.",
    "Original and fixed module sizes and restriction status.",
    "Participant counts for every country-site combination.",
    "Reciprocal within-country site run plan.",
    "Complete WGCNA output for all fixed-gene country runs.",
    "Primary fixed-gene country preservation statistics.",
    "Country preservation summarized by module.",
    "Country preservation summarized by held-out context.",
    "Country run checkpoints, sample sizes and common seeds.",
    "Complete WGCNA output for reciprocal fixed-gene site runs.",
    "Primary fixed-gene site preservation statistics.",
    "Site preservation summarized by module.",
    "Reciprocal preservation summarized within Chile and Colombia.",
    "Site run checkpoints, sample sizes and common seeds.",
    "Direct country comparison between Scripts 15 and 15b when available.",
    "Direct site comparison between Scripts 15 and 15b when available.",
    "Module-level summary of changes caused by the fixed-gene correction.",
    "Definitive structural-preservation summary using identical genes.",
    "Corrected country Zsummary heatmap.",
    "Corrected reciprocal-site Zsummary heatmap.",
    "Old versus corrected country Zsummary comparison.",
    "Old versus corrected site Zsummary comparison.",
    "Methods-ready wording and interpretation boundaries.",
    "Integrated reviewer-facing workbook.",
    "New fixed-gene checkpoint objects; Script 15 objects are not reused.",
    "New permutation-statistic files for reproducible resumption.",
    "Workspace for final figures and supplementary tables.",
    "R session information."
  )
)

safe_write_csv(
  output_manifest,
  file.path(
    OUTDIR,
    "script15b_output_manifest.csv"
  )
)

###############################################################################
# 23) SAVE LIGHT WORKSPACE
###############################################################################

save(
  meta_final,
  softPower,
  modules_use,
  countries,
  site_var,
  country_site_counts,
  sites_per_country,
  multisite_countries,
  site_pair_plan,
  fixed_gene_manifest,
  fixed_module_summary,
  fixed_module_colors,
  country_preservation_raw,
  country_preservation,
  country_summary_by_module,
  country_summary_by_context,
  country_run_log_tbl,
  site_preservation_raw,
  site_preservation,
  site_summary_by_module,
  site_reciprocal_summary,
  site_run_log_tbl,
  country_old_vs_fixed,
  site_old_vs_fixed,
  comparison_summary,
  definitive_summary,
  methods_wording,
  script15b_summary,
  output_manifest,
  file = file.path(
    OUTDIR,
    "workspace",
    "script15b_fixed_geneset_structural_preservation_workspace.RData"
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

cat("\nScript 15b finished successfully.\n")
cat("Main output directory:\n", OUTDIR, "\n")
cat("Samples: ", nrow(datExpr_fixed), "\n", sep = "")
cat("Full-network genes: ", ncol(datExpr_clean), "\n", sep = "")
cat("Fixed-network genes: ", ncol(datExpr_fixed), "\n", sep = "")
cat(
  "Fixed module sizes: ",
  paste0(
    names(fixed_module_sizes),
    "=",
    as.integer(
      fixed_module_sizes
    ),
    collapse = ", "
  ),
  "\n",
  sep = ""
)
cat(
  "Restricted modules: ",
  paste(
    fixed_module_summary$Module[
      fixed_module_summary$
        Module_was_restricted
    ],
    collapse = ", "
  ),
  "\n",
  sep = ""
)
cat("Fixed gene-set seed: ", FIXED_GENESET_SEED, "\n", sep = "")
cat("Common preservation seed: ", MODULE_PRESERVATION_SEED, "\n", sep = "")
cat("Country structural runs: ", length(countries), "\n", sep = "")
cat("Reciprocal site structural runs: ", nrow(site_pair_plan), "\n", sep = "")
cat("Permutations per run: ", N_PERMUTATIONS, "\n", sep = "")
cat(
  "Country module-results with strong preservation: ",
  sum(
    country_preservation$
      Strong_preservation,
    na.rm = TRUE
  ),
  "/",
  nrow(country_preservation),
  "\n",
  sep = ""
)
cat(
  "Country module-results with at least moderate preservation: ",
  sum(
    country_preservation$
      At_least_moderate_preservation,
    na.rm = TRUE
  ),
  "/",
  nrow(country_preservation),
  "\n",
  sep = ""
)
cat(
  "Site module-results with strong preservation: ",
  sum(
    site_preservation$
      Strong_preservation,
    na.rm = TRUE
  ),
  "/",
  nrow(site_preservation),
  "\n",
  sep = ""
)
cat(
  "Site module-results with at least moderate preservation: ",
  sum(
    site_preservation$
      At_least_moderate_preservation,
    na.rm = TRUE
  ),
  "/",
  nrow(site_preservation),
  "\n",
  sep = ""
)
cat("Country_numeric included: FALSE\n")
cat("\nKey outputs:\n")
cat("- tables/input/fixed_gene_manifest_all_modules.csv\n")
cat("- tables/country/country_fixed_geneset_primary_statistics.csv\n")
cat("- tables/country/country_fixed_geneset_summary_by_module.csv\n")
cat("- tables/site/site_fixed_geneset_primary_statistics.csv\n")
cat("- tables/site/site_fixed_geneset_reciprocal_summary_by_country_module.csv\n")
cat("- tables/comparison/script15_vs_script15b_comparison_summary_by_module.csv\n")
cat("- tables/fixed_geneset_structural_preservation_definitive_summary.csv\n")
cat("- WGCNA_15b_Fixed_Gene_Set_Structural_Preservation.xlsx\n")
cat("\nNew checkpoint files are used. Script 15 checkpoints are never loaded.\n")
cat("Interrupted Script 15b runs can be resumed by sourcing this script again.\n")

###############################################################################
# END
###############################################################################


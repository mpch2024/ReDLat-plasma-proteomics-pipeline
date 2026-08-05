###############################################################################
# ReDLat plasma proteomics — WGCNA workflow
# 09. Assess network quality
# Requires: outputs from Scripts 02 and 04
# Produces: modularity, topology and network-quality diagnostics
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

SCRIPT15B_DIR <- file.path(WGCNA_CONFIG$result_root,
  "08_preservation"
)

FULL_WORKSPACE_FILE <- file.path(
  SCRIPT11_DIR,
  "workspace",
  "wgcna_core_collapsed_workspace.RData"
)

SOFT_SCAN_FILE <- file.path(
  SCRIPT11_DIR,
  "soft_threshold",
  "soft_threshold_scan.csv"
)

SOFT_SELECTED_FILE <- file.path(
  SCRIPT11_DIR,
  "soft_threshold",
  "soft_threshold_selected.csv"
)

MODULE_COUNTS_FILE <- file.path(
  SCRIPT11_DIR,
  "tables",
  "module_counts.csv"
)

CORE_SUMMARY_FILE <- file.path(
  SCRIPT11_DIR,
  "tables",
  "wgcna_core_summary.csv"
)

NETWORK_AUDIT_FILE <- file.path(
  SCRIPT11_DIR,
  "qc",
  "network_matrix_alignment_audit.csv"
)

MODULE_REFERENCE_FILE <- file.path(
  SCRIPT13_DIR,
  "tables",
  "module_biological_label_reference.csv"
)

PRESERVATION_SUMMARY_FILE <- file.path(
  SCRIPT15B_DIR,
  "tables",
  "fixed_geneset_structural_preservation_definitive_summary.csv"
)

OUTDIR <- file.path(WGCNA_CONFIG$result_root,
  "09_network_quality"
)

SUBDIRS <- c(
  "tables",
  "tables/input",
  "tables/soft_threshold",
  "tables/modularity",
  "tables/matrix_mixing",
  "tables/connectivity",
  "tables/kME",
  "tables/eigengenes",
  "tables/integration",
  "figures",
  "figures/soft_threshold",
  "figures/modules",
  "figures/modularity",
  "figures/connectivity",
  "figures/eigengenes",
  "figures/integration",
  "workspace",
  "logs"
)

invisible(lapply(
  file.path(OUTDIR, SUBDIRS),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

required_files <- c(
  FULL_WORKSPACE_FILE,
  SOFT_SCAN_FILE,
  SOFT_SELECTED_FILE,
  MODULE_COUNTS_FILE,
  CORE_SUMMARY_FILE,
  NETWORK_AUDIT_FILE,
  MODULE_REFERENCE_FILE
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  stop(
    "Missing required Script 11/13 files:\n",
    paste(missing_files, collapse = "\n"),
    "\n\nRun Scripts 11 and 13 first.",
    call. = FALSE
  )
}

###############################################################################
# 3) PARAMETERS
###############################################################################

EXPECTED_N_SAMPLES <- 639L
EXPECTED_N_GENES <- 9638L
EXPECTED_N_MODULES <- 8L
EXPECTED_SOFT_POWER <- 4

REQUIRE_EXPECTED_DIMENSIONS <- TRUE
REQUIRE_NO_GREY_GENES <- TRUE
REQUIRE_EXACT_MATRIX_ALIGNMENT <- TRUE

# Reproducible edge sample used only for:
# - off-diagonal distributional quantiles;
# - sampled symmetry checks;
# - adjacency-TOM edge concordance.
EDGE_SAMPLE_SIZE <- 1000000L
EDGE_SAMPLE_SEED <- 16001L

TOP_HUBS_PER_MODULE <- 30L

# Module merging in Script 11 used mergeCutHeight = 0.25, corresponding to
# eigengene correlation > 0.75 under 1 - correlation dissimilarity.
EIGENGENE_MERGE_COR_THRESHOLD <- 0.75

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

safe_numeric <- function(x) {
  suppressWarnings(
    as.numeric(
      as.character(x)
    )
  )
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

safe_sd <- function(x) {
  x <- safe_numeric(x)
  x <- x[is.finite(x)]

  if (length(x) < 2) {
    return(NA_real_)
  }

  stats::sd(x)
}

safe_quantile <- function(
    x,
    prob
) {
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

safe_cor <- function(
    x,
    y,
    method = "pearson"
) {
  ok <- is.finite(x) &
    is.finite(y)

  if (
    sum(ok) < 3 ||
    dplyr::n_distinct(
      x[ok]
    ) <= 1 ||
    dplyr::n_distinct(
      y[ok]
    ) <= 1
  ) {
    return(NA_real_)
  }

  suppressWarnings(
    stats::cor(
      x[ok],
      y[ok],
      method = method
    )
  )
}

pick_col <- function(
    df,
    candidates
) {
  if (
    is.null(df) ||
    nrow(
      as.data.frame(df)
    ) == 0
  ) {
    return(NA_character_)
  }

  exact <- candidates[
    candidates %in%
      names(df)
  ][1]

  if (
    length(exact) > 0 &&
    !is.na(exact)
  ) {
    return(exact)
  }

  clean_name <- function(x) {
    tolower(
      gsub(
        "[^a-z0-9]+",
        "",
        x
      )
    )
  }

  names_clean <- clean_name(
    names(df)
  )

  candidates_clean <- clean_name(
    candidates
  )

  idx <- match(
    candidates_clean,
    names_clean,
    nomatch = 0
  )

  idx <- idx[
    idx > 0
  ][1]

  if (
    length(idx) == 0 ||
    is.na(idx)
  ) {
    return(NA_character_)
  }

  names(df)[idx]
}

get_module_colors <- function(modules) {
  modules <- as.character(modules)
  cols <- MODULE_COLORS[modules]

  missing_modules <- modules[
    is.na(cols)
  ]

  if (length(missing_modules) > 0) {
    extra_cols <- grDevices::rainbow(
      length(
        unique(
          missing_modules
        )
      ),
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

safe_sheet_name <- function(
    x,
    used
) {
  base <- substr(
    gsub(
      "[^A-Za-z0-9_]+",
      "_",
      x
    ),
    1,
    31
  )

  candidate <- base
  index <- 1L

  while (candidate %in% used) {
    suffix <- paste0(
      "_",
      index
    )

    candidate <- paste0(
      substr(
        base,
        1,
        31 -
          nchar(suffix)
      ),
      suffix
    )

    index <- index + 1L
  }

  candidate
}

write_workbook_safe <- function(
    named_tables,
    file
) {
  wb <- openxlsx::createWorkbook()
  used <- character()

  for (nm in names(named_tables)) {
    sheet <- safe_sheet_name(
      nm,
      used
    )

    used <- c(
      used,
      sheet
    )

    openxlsx::addWorksheet(
      wb,
      sheet
    )

    openxlsx::writeData(
      wb,
      sheet,
      named_tables[[nm]]
    )

    openxlsx::freezePane(
      wb,
      sheet,
      firstRow = TRUE
    )

    if (
      ncol(
        named_tables[[nm]]
      ) > 0
    ) {
      openxlsx::setColWidths(
        wb,
        sheet,
        cols = seq_len(
          ncol(
            named_tables[[nm]]
          )
        ),
        widths = "auto"
      )
    }
  }

  openxlsx::saveWorkbook(
    wb,
    file,
    overwrite = TRUE
  )
}

matrix_sample_quantiles <- function(
    x,
    matrix_name
) {
  probs <- c(
    0,
    0.01,
    0.05,
    0.25,
    0.50,
    0.75,
    0.95,
    0.99,
    1
  )

  values <- x[
    is.finite(x)
  ]

  quantiles <- stats::quantile(
    values,
    probs = probs,
    na.rm = TRUE,
    names = FALSE
  )

  tibble::tibble(
    Matrix = matrix_name,
    Probability = probs,
    Quantile = as.numeric(
      quantiles
    ),
    N_sampled_edges =
      length(values)
  )
}

###############################################################################
# 5) LOAD TABULAR SCRIPT 11/13/15b OUTPUTS
###############################################################################

soft_scan <- safe_read_csv(
  SOFT_SCAN_FILE
)

soft_selected <- safe_read_csv(
  SOFT_SELECTED_FILE
)

module_counts_file_tbl <- safe_read_csv(
  MODULE_COUNTS_FILE
)

core_summary <- safe_read_csv(
  CORE_SUMMARY_FILE
)

network_audit_script11 <- safe_read_csv(
  NETWORK_AUDIT_FILE
)

module_reference <- safe_read_csv(
  MODULE_REFERENCE_FILE
) %>%
  dplyr::mutate(
    Module = as.character(Module)
  )

preservation_available <- file.exists(
  PRESERVATION_SUMMARY_FILE
)

preservation_summary <- if (
  preservation_available
) {
  safe_read_csv(
    PRESERVATION_SUMMARY_FILE
  ) %>%
    dplyr::mutate(
      Module = as.character(Module)
    )
} else {
  tibble::tibble()
}

###############################################################################
# 6) LOAD FULL SCRIPT 11 WORKSPACE WITH MEMORY CONTROL
###############################################################################

workspace_size_gb <- file.info(
  FULL_WORKSPACE_FILE
)$size /
  1024^3

cat(
  "Loading full Script 11 workspace (",
  round(
    workspace_size_gb,
    2
  ),
  " GB on disk)...\n",
  sep = ""
)

core_env <- new.env(
  parent = emptyenv()
)

load(
  FULL_WORKSPACE_FILE,
  envir = core_env
)

required_objects <- c(
  "datExpr_clean",
  "meta_final",
  "softPower",
  "mergedColors",
  "mergedMEs",
  "adjacency_mat",
  "TOM"
)

missing_objects <- required_objects[
  !vapply(
    required_objects,
    exists,
    logical(1),
    envir = core_env,
    inherits = FALSE
  )
]

if (length(missing_objects) > 0) {
  stop(
    "The full Script 11 workspace is missing required objects: ",
    paste(
      missing_objects,
      collapse = ", "
    ),
    call. = FALSE
  )
}

keep_objects <- unique(
  c(
    required_objects,
    "dynamicColors",
    "module_count_tbl",
    "gene_module_assignment",
    "summary_tbl",
    "network_matrix_audit",
    "module_assignment_audit",
    "eigengene_alignment_audit"
  )
)

remove_objects <- setdiff(
  ls(core_env),
  keep_objects
)

if (length(remove_objects) > 0) {
  rm(
    list = remove_objects,
    envir = core_env
  )
}

gc(verbose = FALSE)

datExpr_clean <- core_env$datExpr_clean
meta_final <- core_env$meta_final
softPower <- safe_numeric(
  core_env$softPower
)[1]
mergedColors <- as.character(
  core_env$mergedColors
)
names(mergedColors) <- names(
  core_env$mergedColors
)
mergedMEs <- as.data.frame(
  core_env$mergedMEs,
  check.names = FALSE
)
adjacency_mat <- core_env$adjacency_mat
TOM <- core_env$TOM

dynamicColors <- if (
  exists(
    "dynamicColors",
    envir = core_env,
    inherits = FALSE
  )
) {
  as.character(
    core_env$dynamicColors
  )
} else {
  NULL
}

###############################################################################
# 7) STRICT ALIGNMENT AUDIT
###############################################################################

datExpr_clean <- as.data.frame(
  datExpr_clean,
  check.names = FALSE
)

gene_ids <- colnames(
  datExpr_clean
)

sample_ids <- rownames(
  datExpr_clean
)

if (
  is.null(gene_ids) ||
  anyDuplicated(gene_ids) > 0
) {
  stop(
    "datExpr_clean gene names are absent or duplicated.",
    call. = FALSE
  )
}

if (
  is.null(sample_ids) ||
  anyDuplicated(sample_ids) > 0
) {
  stop(
    "datExpr_clean sample names are absent or duplicated.",
    call. = FALSE
  )
}

if (
  is.null(names(mergedColors)) ||
  !identical(
    names(mergedColors),
    gene_ids
  )
) {
  stop(
    "mergedColors is not exactly aligned with datExpr_clean columns.",
    call. = FALSE
  )
}

if (
  nrow(adjacency_mat) !=
    length(gene_ids) ||
  ncol(adjacency_mat) !=
    length(gene_ids)
) {
  stop(
    "adjacency_mat dimensions do not match datExpr_clean.",
    call. = FALSE
  )
}

if (
  nrow(TOM) !=
    length(gene_ids) ||
  ncol(TOM) !=
    length(gene_ids)
) {
  stop(
    "TOM dimensions do not match datExpr_clean.",
    call. = FALSE
  )
}

adjacency_names_exact <- (
  !is.null(
    rownames(adjacency_mat)
  ) &&
  !is.null(
    colnames(adjacency_mat)
  ) &&
  identical(
    rownames(adjacency_mat),
    gene_ids
  ) &&
  identical(
    colnames(adjacency_mat),
    gene_ids
  )
)

tom_names_exact <- (
  !is.null(rownames(TOM)) &&
  !is.null(colnames(TOM)) &&
  identical(
    rownames(TOM),
    gene_ids
  ) &&
  identical(
    colnames(TOM),
    gene_ids
  )
)

if (
  REQUIRE_EXACT_MATRIX_ALIGNMENT &&
  (
    !adjacency_names_exact ||
    !tom_names_exact
  )
) {
  stop(
    "Dense matrices are not exactly name-aligned with datExpr_clean. ",
    "The script will not reorder a 9,638 x 9,638 matrix silently.",
    call. = FALSE
  )
}

if (
  nrow(mergedMEs) !=
    nrow(datExpr_clean)
) {
  stop(
    "mergedMEs rows do not match datExpr_clean samples.",
    call. = FALSE
  )
}

if (
  !is.null(
    rownames(mergedMEs)
  ) &&
  !identical(
    rownames(mergedMEs),
    sample_ids
  )
) {
  stop(
    "mergedMEs is not exactly sample-aligned with datExpr_clean.",
    call. = FALSE
  )
}

modules_use <- sort(
  unique(
    mergedColors[
      !is.na(mergedColors) &
        mergedColors != "grey"
    ]
  )
)

grey_gene_count <- sum(
  mergedColors == "grey",
  na.rm = TRUE
)

if (
  REQUIRE_NO_GREY_GENES &&
  grey_gene_count != 0
) {
  stop(
    "Expected zero grey genes, but found ",
    grey_gene_count,
    ".",
    call. = FALSE
  )
}

if (REQUIRE_EXPECTED_DIMENSIONS) {
  if (
    nrow(datExpr_clean) !=
      EXPECTED_N_SAMPLES
  ) {
    stop(
      "Expected ",
      EXPECTED_N_SAMPLES,
      " samples, but found ",
      nrow(datExpr_clean),
      ".",
      call. = FALSE
    )
  }

  if (
    ncol(datExpr_clean) !=
      EXPECTED_N_GENES
  ) {
    stop(
      "Expected ",
      EXPECTED_N_GENES,
      " genes, but found ",
      ncol(datExpr_clean),
      ".",
      call. = FALSE
    )
  }

  if (
    length(modules_use) !=
      EXPECTED_N_MODULES
  ) {
    stop(
      "Expected ",
      EXPECTED_N_MODULES,
      " modules, but found ",
      length(modules_use),
      ".",
      call. = FALSE
    )
  }

  if (
    !isTRUE(
      all.equal(
        softPower,
        EXPECTED_SOFT_POWER
      )
    )
  ) {
    stop(
      "Expected soft power ",
      EXPECTED_SOFT_POWER,
      ", but found ",
      softPower,
      ".",
      call. = FALSE
    )
  }
}

input_alignment_audit <- tibble::tibble(
  Metric = c(
    "Workspace",
    "Workspace_size_GB",
    "N_samples",
    "N_genes",
    "N_modules_excluding_grey",
    "Modules",
    "Grey_gene_count",
    "Soft_power",
    "Expression_missing_values",
    "MergedColors_exact_alignment",
    "Adjacency_dimensions",
    "Adjacency_names_exact_alignment",
    "TOM_dimensions",
    "TOM_names_exact_alignment",
    "Eigengene_sample_alignment",
    "Country_numeric_used"
  ),
  Value = c(
    FULL_WORKSPACE_FILE,
    as.character(
      workspace_size_gb
    ),
    as.character(
      nrow(datExpr_clean)
    ),
    as.character(
      ncol(datExpr_clean)
    ),
    as.character(
      length(modules_use)
    ),
    paste(
      modules_use,
      collapse = ", "
    ),
    as.character(
      grey_gene_count
    ),
    as.character(
      softPower
    ),
    as.character(
      sum(
        is.na(datExpr_clean)
      )
    ),
    as.character(
      identical(
        names(mergedColors),
        gene_ids
      )
    ),
    paste(
      dim(adjacency_mat),
      collapse = " x "
    ),
    as.character(
      adjacency_names_exact
    ),
    paste(
      dim(TOM),
      collapse = " x "
    ),
    as.character(
      tom_names_exact
    ),
    as.character(
      nrow(mergedMEs) ==
        nrow(datExpr_clean)
    ),
    "FALSE"
  )
)

safe_write_csv(
  input_alignment_audit,
  file.path(
    OUTDIR,
    "tables",
    "input",
    "script16_input_alignment_audit.csv"
  )
)

safe_write_csv(
  network_audit_script11,
  file.path(
    OUTDIR,
    "tables",
    "input",
    "script11_network_matrix_alignment_audit.csv"
  )
)

safe_write_csv(
  core_summary,
  file.path(
    OUTDIR,
    "tables",
    "input",
    "script11_wgcna_core_summary.csv"
  )
)

###############################################################################
# 8) SOFT-THRESHOLD DIAGNOSTICS
###############################################################################

power_col <- pick_col(
  soft_scan,
  c(
    "Power",
    "power",
    "softPower",
    "beta"
  )
)

signed_r2_col <- pick_col(
  soft_scan,
  c(
    "signed_R2",
    "signed.R2",
    "Signed_R2",
    "signed_scale_free_R2",
    "SFT.R.sq"
  )
)

raw_r2_col <- pick_col(
  soft_scan,
  c(
    "SFT.R.sq",
    "Rsquared.SFT",
    "scale_free_R2",
    "ScaleFreeFit",
    "R2"
  )
)

slope_col <- pick_col(
  soft_scan,
  c(
    "slope",
    "Slope"
  )
)

mean_k_col <- pick_col(
  soft_scan,
  c(
    "mean.k.",
    "mean.k",
    "mean_connectivity",
    "MeanConnectivity",
    "mean.k.."
  )
)

median_k_col <- pick_col(
  soft_scan,
  c(
    "median.k.",
    "median.k",
    "median_connectivity"
  )
)

max_k_col <- pick_col(
  soft_scan,
  c(
    "max.k.",
    "max.k",
    "max_connectivity"
  )
)

if (
  is.na(power_col) ||
  is.na(mean_k_col) ||
  (
    is.na(signed_r2_col) &&
    is.na(raw_r2_col)
  )
) {
  stop(
    "Required soft-threshold columns were not detected.",
    call. = FALSE
  )
}

soft_scan_clean <- soft_scan %>%
  dplyr::transmute(
    Power =
      safe_numeric(
        .data[[power_col]]
      ),
    Signed_scale_free_R2 =
      if (!is.na(
        signed_r2_col
      )) {
        safe_numeric(
          .data[[signed_r2_col]]
        )
      } else {
        safe_numeric(
          .data[[raw_r2_col]]
        )
      },
    Raw_scale_free_R2 =
      if (!is.na(raw_r2_col)) {
        safe_numeric(
          .data[[raw_r2_col]]
        )
      } else {
        NA_real_
      },
    Slope =
      if (!is.na(slope_col)) {
        safe_numeric(
          .data[[slope_col]]
        )
      } else {
        NA_real_
      },
    Mean_connectivity =
      safe_numeric(
        .data[[mean_k_col]]
      ),
    Median_connectivity =
      if (!is.na(
        median_k_col
      )) {
        safe_numeric(
          .data[[median_k_col]]
        )
      } else {
        NA_real_
      },
    Maximum_connectivity =
      if (!is.na(max_k_col)) {
        safe_numeric(
          .data[[max_k_col]]
        )
      } else {
        NA_real_
      }
  ) %>%
  dplyr::filter(
    is.finite(Power)
  ) %>%
  dplyr::arrange(Power)

selected_soft_row <- soft_scan_clean %>%
  dplyr::filter(
    abs(
      Power -
        softPower
    ) <
      1e-9
  )

if (nrow(selected_soft_row) != 1) {
  stop(
    "The selected soft-power row was not uniquely identified.",
    call. = FALSE
  )
}

first_power_R2_ge_0_80 <- soft_scan_clean %>%
  dplyr::filter(
    Signed_scale_free_R2 >=
      0.80
  ) %>%
  dplyr::slice_min(
    Power,
    n = 1,
    with_ties = FALSE
  ) %>%
  dplyr::pull(Power)

if (
  length(
    first_power_R2_ge_0_80
  ) == 0
) {
  first_power_R2_ge_0_80 <-
    NA_real_
}

first_power_R2_ge_0_90 <- soft_scan_clean %>%
  dplyr::filter(
    Signed_scale_free_R2 >=
      0.90
  ) %>%
  dplyr::slice_min(
    Power,
    n = 1,
    with_ties = FALSE
  ) %>%
  dplyr::pull(Power)

if (
  length(
    first_power_R2_ge_0_90
  ) == 0
) {
  first_power_R2_ge_0_90 <-
    NA_real_
}

soft_threshold_summary <- tibble::tibble(
  Metric = c(
    "Selected_power",
    "Signed_scale_free_R2_at_selected_power",
    "Raw_scale_free_R2_at_selected_power",
    "Slope_at_selected_power",
    "Mean_connectivity_at_selected_power",
    "Median_connectivity_at_selected_power",
    "Maximum_connectivity_at_selected_power",
    "First_power_with_signed_R2_ge_0_80",
    "First_power_with_signed_R2_ge_0_90"
  ),
  Value = c(
    softPower,
    selected_soft_row$
      Signed_scale_free_R2,
    selected_soft_row$
      Raw_scale_free_R2,
    selected_soft_row$Slope,
    selected_soft_row$
      Mean_connectivity,
    selected_soft_row$
      Median_connectivity,
    selected_soft_row$
      Maximum_connectivity,
    first_power_R2_ge_0_80,
    first_power_R2_ge_0_90
  )
)

safe_write_csv(
  soft_scan_clean,
  file.path(
    OUTDIR,
    "tables",
    "soft_threshold",
    "soft_threshold_scan_clean.csv"
  )
)

safe_write_csv(
  selected_soft_row,
  file.path(
    OUTDIR,
    "tables",
    "soft_threshold",
    "soft_threshold_selected_recomputed.csv"
  )
)

safe_write_csv(
  soft_threshold_summary,
  file.path(
    OUTDIR,
    "tables",
    "soft_threshold",
    "soft_threshold_quality_summary.csv"
  )
)

###############################################################################
# 9) MODULE ARCHITECTURE AND CONCENTRATION
###############################################################################

module_size_tbl <- tibble::tibble(
  Module = modules_use,
  N_genes = as.integer(
    table(
      factor(
        mergedColors,
        levels = modules_use
      )
    )
  )
) %>%
  dplyr::mutate(
    Proportion_of_network =
      N_genes /
      sum(N_genes),
    Module_color =
      unname(
        get_module_colors(Module)
      )
  ) %>%
  dplyr::left_join(
    module_reference %>%
      dplyr::select(
        Module,
        Biological_label
      ),
    by = "Module"
  ) %>%
  dplyr::arrange(
    dplyr::desc(N_genes)
  )

module_proportions <- module_size_tbl$
  Proportion_of_network

module_entropy <- -sum(
  module_proportions *
    log(module_proportions)
)

normalized_module_entropy <- module_entropy /
  log(
    length(module_proportions)
  )

effective_module_number <- exp(
  module_entropy
)

module_HHI <- sum(
  module_proportions^2
)

module_architecture_summary <- tibble::tibble(
  Metric = c(
    "N_modules",
    "N_grey_genes",
    "Largest_module",
    "Largest_module_size",
    "Largest_module_fraction",
    "Smallest_module",
    "Smallest_module_size",
    "Module_size_median",
    "Module_size_IQR",
    "Shannon_entropy",
    "Normalized_Shannon_entropy",
    "Effective_number_of_modules",
    "Herfindahl_Hirschman_index",
    "Initial_dynamic_module_count",
    "Final_merged_module_count"
  ),
  Value = c(
    length(modules_use),
    grey_gene_count,
    module_size_tbl$Module[[1]],
    module_size_tbl$N_genes[[1]],
    module_size_tbl$
      Proportion_of_network[[1]],
    module_size_tbl$Module[[
      nrow(module_size_tbl)
    ]],
    module_size_tbl$N_genes[[
      nrow(module_size_tbl)
    ]],
    stats::median(
      module_size_tbl$N_genes
    ),
    stats::IQR(
      module_size_tbl$N_genes
    ),
    module_entropy,
    normalized_module_entropy,
    effective_module_number,
    module_HHI,
    if (!is.null(dynamicColors)) {
      dplyr::n_distinct(
        dynamicColors[
          !is.na(dynamicColors) &
            dynamicColors != "grey"
        ]
      )
    } else {
      NA_integer_
    },
    length(modules_use)
  )
)

safe_write_csv(
  module_size_tbl,
  file.path(
    OUTDIR,
    "tables",
    "modularity",
    "final_module_sizes_and_proportions.csv"
  )
)

safe_write_csv(
  module_architecture_summary,
  file.path(
    OUTDIR,
    "tables",
    "modularity",
    "module_architecture_concentration_summary.csv"
  )
)

###############################################################################
# 10) EXACT MODULE BLOCK STRENGTHS FOR ADJACENCY AND TOM
###############################################################################

n_genes <- length(
  gene_ids
)

n_modules <- length(
  modules_use
)

module_index <- match(
  mergedColors,
  modules_use
)

if (anyNA(module_index)) {
  stop(
    "Some genes could not be assigned to the final non-grey modules.",
    call. = FALSE
  )
}

membership_matrix <- matrix(
  0,
  nrow = n_genes,
  ncol = n_modules,
  dimnames = list(
    gene_ids,
    modules_use
  )
)

membership_matrix[
  cbind(
    seq_len(n_genes),
    module_index
  )
] <- 1

adjacency_diag <- diag(
  adjacency_mat
)

tom_diag <- diag(TOM)

cat(
  "Calculating exact adjacency strengths by module...\n"
)

adjacency_strength_by_module <-
  adjacency_mat %*%
  membership_matrix

adjacency_strength_by_module[
  cbind(
    seq_len(n_genes),
    module_index
  )
] <-
  adjacency_strength_by_module[
    cbind(
      seq_len(n_genes),
      module_index
    )
  ] -
  adjacency_diag

cat(
  "Calculating exact TOM strengths by module...\n"
)

tom_strength_by_module <-
  TOM %*%
  membership_matrix

tom_strength_by_module[
  cbind(
    seq_len(n_genes),
    module_index
  )
] <-
  tom_strength_by_module[
    cbind(
      seq_len(n_genes),
      module_index
    )
  ] -
  tom_diag

adjacency_block_sum <- crossprod(
  membership_matrix,
  adjacency_strength_by_module
)

tom_block_sum <- crossprod(
  membership_matrix,
  tom_strength_by_module
)

module_sizes_named <- setNames(
  module_size_tbl$N_genes[
    match(
      modules_use,
      module_size_tbl$Module
    )
  ],
  modules_use
)

make_block_table <- function(
    block_sum,
    matrix_name
) {
  purrr::map_dfr(
    modules_use,
    function(module_row) {
      purrr::map_dfr(
        modules_use,
        function(module_col) {
          n_row <- module_sizes_named[[
            module_row
          ]]

          n_col <- module_sizes_named[[
            module_col
          ]]

          denominator <- if (
            module_row ==
              module_col
          ) {
            n_row *
              (
                n_row - 1
              )
          } else {
            n_row * n_col
          }

          block_value <- block_sum[
            module_row,
            module_col
          ]

          tibble::tibble(
            Matrix = matrix_name,
            Module_row =
              module_row,
            Module_column =
              module_col,
            N_row = n_row,
            N_column = n_col,
            Directed_block_sum =
              as.numeric(
                block_value
              ),
            Possible_directed_edges =
              denominator,
            Mean_edge_weight =
              if (
                denominator > 0
              ) {
                as.numeric(
                  block_value
                ) /
                  denominator
              } else {
                NA_real_
              },
            Same_module =
              module_row ==
              module_col
          )
        }
      )
    }
  )
}

adjacency_block_tbl <- make_block_table(
  adjacency_block_sum,
  "Adjacency"
)

tom_block_tbl <- make_block_table(
  tom_block_sum,
  "TOM"
)

safe_write_csv(
  adjacency_block_tbl,
  file.path(
    OUTDIR,
    "tables",
    "matrix_mixing",
    "adjacency_module_block_mixing.csv"
  )
)

safe_write_csv(
  tom_block_tbl,
  file.path(
    OUTDIR,
    "tables",
    "matrix_mixing",
    "TOM_module_block_mixing.csv"
  )
)

###############################################################################
# 11) EXACT POST HOC WEIGHTED MODULARITY Q
###############################################################################

k_total_adjacency <- rowSums(
  adjacency_strength_by_module
)

two_m <- sum(
  k_total_adjacency
)

if (
  !is.finite(two_m) ||
  two_m <= 0
) {
  stop(
    "Total off-diagonal adjacency strength is not positive.",
    call. = FALSE
  )
}

module_modularity_tbl <- purrr::map_dfr(
  modules_use,
  function(mod) {
    idx <- which(
      mergedColors == mod
    )

    within_directed <- adjacency_block_sum[
      mod,
      mod
    ]

    module_volume <- sum(
      k_total_adjacency[idx]
    )

    within_fraction_of_total <-
      within_directed /
      two_m

    expected_fraction <-
      (
        module_volume /
          two_m
      )^2

    q_contribution <-
      within_fraction_of_total -
      expected_fraction

    external_strength <-
      module_volume -
      within_directed

    conductance <- if (
      module_volume > 0
    ) {
      external_strength /
        module_volume
    } else {
      NA_real_
    }

    tibble::tibble(
      Module = mod,
      N_genes = length(idx),
      Module_volume =
        module_volume,
      Module_strength_share =
        module_volume /
        two_m,
      Within_directed_strength =
        within_directed,
      Within_strength_share_of_network =
        within_fraction_of_total,
      Expected_within_share_under_strength_null =
        expected_fraction,
      Modularity_Q_contribution =
        q_contribution,
      External_strength =
        external_strength,
      Conductance =
        conductance
    )
  }
) %>%
  dplyr::left_join(
    module_reference %>%
      dplyr::select(
        Module,
        Biological_label
      ),
    by = "Module"
  ) %>%
  dplyr::arrange(
    dplyr::desc(
      Modularity_Q_contribution
    )
  )

weighted_modularity_Q <- sum(
  module_modularity_tbl$
    Modularity_Q_contribution
)

global_mean_adjacency <- two_m /
  (
    n_genes *
      (
        n_genes - 1
      )
  )

global_modularity_tbl <- tibble::tibble(
  Metric = c(
    "Posthoc_weighted_modularity_Q",
    "Total_directed_offdiagonal_strength_2m",
    "Global_mean_offdiagonal_adjacency",
    "N_modules",
    "N_genes",
    "Modularity_definition",
    "Interpretation_boundary"
  ),
  Value = c(
    as.character(
      weighted_modularity_Q
    ),
    as.character(two_m),
    as.character(
      global_mean_adjacency
    ),
    as.character(n_modules),
    as.character(n_genes),
    paste(
      "Q = sum_s[e_ss/(2m) - (a_s/(2m))^2],",
      "calculated from the signed weighted adjacency after setting",
      "self-adjacency contributions to zero analytically."
    ),
    paste(
      "Post hoc descriptive diagnostic only;",
      "WGCNA did not optimize Q and no universal biological threshold is assumed."
    )
  )
)

safe_write_csv(
  global_modularity_tbl,
  file.path(
    OUTDIR,
    "tables",
    "modularity",
    "posthoc_weighted_modularity_Q.csv"
  )
)

safe_write_csv(
  module_modularity_tbl,
  file.path(
    OUTDIR,
    "tables",
    "modularity",
    "module_specific_modularity_contributions.csv"
  )
)

###############################################################################
# 12) NODE-LEVEL CONNECTIVITY AND PARTICIPATION
###############################################################################

k_within_adjacency <-
  adjacency_strength_by_module[
    cbind(
      seq_len(n_genes),
      module_index
    )
  ]

k_external_adjacency <-
  k_total_adjacency -
  k_within_adjacency

within_connectivity_fraction <- ifelse(
  k_total_adjacency > 0,
  k_within_adjacency /
    k_total_adjacency,
  NA_real_
)

strength_proportions <- adjacency_strength_by_module /
  k_total_adjacency

strength_proportions[
  !is.finite(
    strength_proportions
  )
] <- NA_real_

participation_coefficient <- 1 -
  rowSums(
    strength_proportions^2,
    na.rm = TRUE
  )

k_total_tom <- rowSums(
  tom_strength_by_module
)

k_within_tom <- tom_strength_by_module[
  cbind(
    seq_len(n_genes),
    module_index
  )
]

k_external_tom <- k_total_tom -
  k_within_tom

node_quality_tbl <- tibble::tibble(
  Gene = gene_ids,
  Module = mergedColors,
  kTotal_adjacency =
    k_total_adjacency,
  kWithin_adjacency =
    k_within_adjacency,
  kExternal_adjacency =
    k_external_adjacency,
  Within_connectivity_fraction =
    within_connectivity_fraction,
  Participation_coefficient =
    participation_coefficient,
  kTotal_TOM =
    k_total_tom,
  kWithin_TOM =
    k_within_tom,
  kExternal_TOM =
    k_external_tom
) %>%
  dplyr::group_by(Module) %>%
  dplyr::mutate(
    kWithin_z =
      if (
        safe_sd(
          kWithin_adjacency
        ) > 0
      ) {
        (
          kWithin_adjacency -
            safe_mean(
              kWithin_adjacency
            )
        ) /
          safe_sd(
            kWithin_adjacency
          )
      } else {
        NA_real_
      },
    kWithin_rank =
      rank(
        -kWithin_adjacency,
        ties.method = "min"
      )
  ) %>%
  dplyr::ungroup()

###############################################################################
# 13) MODULE MEMBERSHIP kME
###############################################################################

me_names <- colnames(
  mergedMEs
)

if (
  is.null(me_names) ||
  length(me_names) == 0
) {
  stop(
    "mergedMEs has no eigengene column names.",
    call. = FALSE
  )
}

me_module_names <- sub(
  "^ME",
  "",
  me_names
)

if (
  !all(
    modules_use %in%
      me_module_names
  )
) {
  stop(
    "Not all modules have a corresponding merged eigengene.",
    call. = FALSE
  )
}

mergedMEs_ordered <- mergedMEs[
  ,
  match(
    modules_use,
    me_module_names
  ),
  drop = FALSE
]

colnames(
  mergedMEs_ordered
) <- modules_use

cat(
  "Calculating gene-module membership correlations...\n"
)

kME_matrix <- stats::cor(
  datExpr_clean,
  mergedMEs_ordered,
  use = "pairwise.complete.obs",
  method = "pearson"
)

assigned_kME <- kME_matrix[
  cbind(
    seq_len(n_genes),
    module_index
  )
]

node_quality_tbl <- node_quality_tbl %>%
  dplyr::mutate(
    Assigned_kME =
      assigned_kME,
    Abs_assigned_kME =
      abs(Assigned_kME)
  ) %>%
  dplyr::group_by(Module) %>%
  dplyr::mutate(
    Abs_kME_rank =
      rank(
        -Abs_assigned_kME,
        ties.method = "min"
      )
  ) %>%
  dplyr::ungroup()

module_connectivity_summary <- node_quality_tbl %>%
  dplyr::group_by(Module) %>%
  dplyr::summarise(
    N_genes = dplyr::n(),
    Mean_kTotal_adjacency =
      safe_mean(
        kTotal_adjacency
      ),
    Median_kTotal_adjacency =
      safe_median(
        kTotal_adjacency
      ),
    Mean_kWithin_adjacency =
      safe_mean(
        kWithin_adjacency
      ),
    Median_kWithin_adjacency =
      safe_median(
        kWithin_adjacency
      ),
    Mean_kExternal_adjacency =
      safe_mean(
        kExternal_adjacency
      ),
    Median_kExternal_adjacency =
      safe_median(
        kExternal_adjacency
      ),
    Mean_within_connectivity_fraction =
      safe_mean(
        Within_connectivity_fraction
      ),
    Median_within_connectivity_fraction =
      safe_median(
        Within_connectivity_fraction
      ),
    Mean_participation_coefficient =
      safe_mean(
        Participation_coefficient
      ),
    Median_participation_coefficient =
      safe_median(
        Participation_coefficient
      ),
    Mean_kWithin_TOM =
      safe_mean(kWithin_TOM),
    Median_kWithin_TOM =
      safe_median(
        kWithin_TOM
      ),
    .groups = "drop"
  )

module_kME_summary <- node_quality_tbl %>%
  dplyr::group_by(Module) %>%
  dplyr::summarise(
    N_genes = dplyr::n(),
    Mean_signed_kME =
      safe_mean(
        Assigned_kME
      ),
    Median_signed_kME =
      safe_median(
        Assigned_kME
      ),
    Mean_abs_kME =
      safe_mean(
        Abs_assigned_kME
      ),
    Median_abs_kME =
      safe_median(
        Abs_assigned_kME
      ),
    Q25_abs_kME =
      safe_quantile(
        Abs_assigned_kME,
        0.25
      ),
    Q75_abs_kME =
      safe_quantile(
        Abs_assigned_kME,
        0.75
      ),
    Maximum_abs_kME =
      safe_max(
        Abs_assigned_kME
      ),
    N_abs_kME_ge_0_50 =
      sum(
        Abs_assigned_kME >=
          0.50,
        na.rm = TRUE
      ),
    N_abs_kME_ge_0_70 =
      sum(
        Abs_assigned_kME >=
          0.70,
        na.rm = TRUE
      ),
    N_abs_kME_ge_0_80 =
      sum(
        Abs_assigned_kME >=
          0.80,
        na.rm = TRUE
      ),
    Proportion_abs_kME_ge_0_50 =
      mean(
        Abs_assigned_kME >=
          0.50,
        na.rm = TRUE
      ),
    kWithin_abs_kME_Pearson =
      safe_cor(
        kWithin_adjacency,
        Abs_assigned_kME,
        method = "pearson"
      ),
    kWithin_abs_kME_Spearman =
      safe_cor(
        kWithin_adjacency,
        Abs_assigned_kME,
        method = "spearman"
      ),
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    module_reference %>%
      dplyr::select(
        Module,
        Biological_label
      ),
    by = "Module"
  ) %>%
  dplyr::arrange(
    dplyr::desc(
      Median_abs_kME
    )
  )

top_hubs_by_kWithin <- node_quality_tbl %>%
  dplyr::group_by(Module) %>%
  dplyr::slice_max(
    order_by =
      kWithin_adjacency,
    n = TOP_HUBS_PER_MODULE,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(
    Module,
    dplyr::desc(
      kWithin_adjacency
    )
  )

top_hubs_by_kME <- node_quality_tbl %>%
  dplyr::group_by(Module) %>%
  dplyr::slice_max(
    order_by =
      Abs_assigned_kME,
    n = TOP_HUBS_PER_MODULE,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(
    Module,
    dplyr::desc(
      Abs_assigned_kME
    )
  )

hub_overlap_summary <- node_quality_tbl %>%
  dplyr::mutate(
    Top_kWithin =
      kWithin_rank <=
      TOP_HUBS_PER_MODULE,
    Top_abs_kME =
      Abs_kME_rank <=
      TOP_HUBS_PER_MODULE
  ) %>%
  dplyr::group_by(Module) %>%
  dplyr::summarise(
    Top_N =
      TOP_HUBS_PER_MODULE,
    N_overlap =
      sum(
        Top_kWithin &
          Top_abs_kME
      ),
    Jaccard_overlap =
      N_overlap /
      (
        sum(Top_kWithin) +
          sum(Top_abs_kME) -
          N_overlap
      ),
    .groups = "drop"
  )

safe_write_csv(
  node_quality_tbl,
  file.path(
    OUTDIR,
    "tables",
    "connectivity",
    "gene_level_network_quality_metrics.csv"
  )
)

safe_write_csv(
  module_connectivity_summary,
  file.path(
    OUTDIR,
    "tables",
    "connectivity",
    "module_connectivity_quality_summary.csv"
  )
)

safe_write_csv(
  module_kME_summary,
  file.path(
    OUTDIR,
    "tables",
    "kME",
    "module_kME_quality_summary.csv"
  )
)

safe_write_csv(
  top_hubs_by_kWithin,
  file.path(
    OUTDIR,
    "tables",
    "kME",
    "top_hubs_by_intramodular_connectivity.csv"
  )
)

safe_write_csv(
  top_hubs_by_kME,
  file.path(
    OUTDIR,
    "tables",
    "kME",
    "top_hubs_by_absolute_kME.csv"
  )
)

safe_write_csv(
  hub_overlap_summary,
  file.path(
    OUTDIR,
    "tables",
    "kME",
    "hub_definition_overlap_summary.csv"
  )
)

###############################################################################
# 14) MODULE-LEVEL ADJACENCY AND TOM SEPARATION
###############################################################################

module_separation_tbl <- purrr::map_dfr(
  modules_use,
  function(mod) {
    idx <- which(
      mergedColors == mod
    )

    n_mod <- length(idx)
    n_out <- n_genes -
      n_mod

    adj_within_sum <-
      adjacency_block_sum[
        mod,
        mod
      ]

    adj_external_sum <-
      sum(
        adjacency_block_sum[
          mod,
          modules_use != mod,
          drop = TRUE
        ]
      )

    tom_within_sum <-
      tom_block_sum[
        mod,
        mod
      ]

    tom_external_sum <-
      sum(
        tom_block_sum[
          mod,
          modules_use != mod,
          drop = TRUE
        ]
      )

    adj_within_mean <- if (
      n_mod > 1
    ) {
      adj_within_sum /
        (
          n_mod *
            (
              n_mod - 1
            )
        )
    } else {
      NA_real_
    }

    adj_external_mean <- if (
      n_mod > 0 &&
      n_out > 0
    ) {
      adj_external_sum /
        (
          n_mod *
            n_out
        )
    } else {
      NA_real_
    }

    tom_within_mean <- if (
      n_mod > 1
    ) {
      tom_within_sum /
        (
          n_mod *
            (
              n_mod - 1
            )
        )
    } else {
      NA_real_
    }

    tom_external_mean <- if (
      n_mod > 0 &&
      n_out > 0
    ) {
      tom_external_sum /
        (
          n_mod *
            n_out
        )
    } else {
      NA_real_
    }

    tibble::tibble(
      Module = mod,
      N_genes = n_mod,
      Mean_within_adjacency =
        adj_within_mean,
      Mean_external_adjacency =
        adj_external_mean,
      Adjacency_within_external_ratio =
        if (
          is.finite(
            adj_external_mean
          ) &&
          adj_external_mean > 0
        ) {
          adj_within_mean /
            adj_external_mean
        } else {
          NA_real_
        },
      Mean_within_TOM =
        tom_within_mean,
      Mean_external_TOM =
        tom_external_mean,
      TOM_within_external_ratio =
        if (
          is.finite(
            tom_external_mean
          ) &&
          tom_external_mean > 0
        ) {
          tom_within_mean /
            tom_external_mean
        } else {
          NA_real_
        }
    )
  }
) %>%
  dplyr::left_join(
    module_modularity_tbl %>%
      dplyr::select(
        Module,
        Modularity_Q_contribution,
        Conductance
      ),
    by = "Module"
  ) %>%
  dplyr::left_join(
    module_connectivity_summary,
    by = c(
      "Module",
      "N_genes"
    )
  ) %>%
  dplyr::left_join(
    module_kME_summary %>%
      dplyr::select(
        Module,
        Median_abs_kME,
        Mean_abs_kME,
        kWithin_abs_kME_Spearman
      ),
    by = "Module"
  ) %>%
  dplyr::left_join(
    module_reference %>%
      dplyr::select(
        Module,
        Biological_label
      ),
    by = "Module"
  ) %>%
  dplyr::arrange(
    dplyr::desc(
      Adjacency_within_external_ratio
    )
  )

safe_write_csv(
  module_separation_tbl,
  file.path(
    OUTDIR,
    "tables",
    "modularity",
    "module_internal_external_separation_summary.csv"
  )
)

###############################################################################
# 15) MODULE EIGENGENE DIAGNOSTICS
###############################################################################

eigengene_cor <- stats::cor(
  mergedMEs_ordered,
  use = "pairwise.complete.obs",
  method = "pearson"
)

eigengene_cor_tbl <- as.data.frame(
  eigengene_cor,
  check.names = FALSE
) %>%
  tibble::rownames_to_column(
    "Module"
  ) %>%
  tibble::as_tibble()

eigengene_pair_tbl <- purrr::map_dfr(
  seq_len(
    n_modules - 1
  ),
  function(i) {
    purrr::map_dfr(
      seq.int(
        i + 1,
        n_modules
      ),
      function(j) {
        tibble::tibble(
          Module_1 =
            modules_use[[i]],
          Module_2 =
            modules_use[[j]],
          Eigengene_correlation =
            eigengene_cor[
              i,
              j
            ],
          Absolute_eigengene_correlation =
            abs(
              eigengene_cor[
                i,
                j
              ]
            ),
          Signed_eigengene_dissimilarity =
            1 -
            eigengene_cor[
              i,
              j
            ],
          Exceeds_merge_correlation_threshold =
            eigengene_cor[
              i,
              j
            ] >
            EIGENGENE_MERGE_COR_THRESHOLD
        )
      }
    )
  }
) %>%
  dplyr::arrange(
    dplyr::desc(
      Eigengene_correlation
    )
  )

highest_signed_pair <-
  eigengene_pair_tbl %>%
  dplyr::slice_max(
    Eigengene_correlation,
    n = 1,
    with_ties = FALSE
  )

highest_absolute_pair <-
  eigengene_pair_tbl %>%
  dplyr::slice_max(
    Absolute_eigengene_correlation,
    n = 1,
    with_ties = FALSE
  )

eigengene_quality_summary <- tibble::tibble(
  Metric = c(
    "N_eigengenes",
    "Maximum_positive_intermodule_correlation",
    "Maximum_positive_pair",
    "Maximum_absolute_intermodule_correlation",
    "Maximum_absolute_pair",
    "Minimum_signed_eigengene_dissimilarity",
    "Merge_correlation_threshold",
    "N_pairs_exceeding_merge_threshold"
  ),
  Value = c(
    as.character(n_modules),
    as.character(
      highest_signed_pair$
        Eigengene_correlation
    ),
    paste(
      highest_signed_pair$
        Module_1,
      highest_signed_pair$
        Module_2,
      sep = " - "
    ),
    as.character(
      highest_absolute_pair$
        Absolute_eigengene_correlation
    ),
    paste(
      highest_absolute_pair$
        Module_1,
      highest_absolute_pair$
        Module_2,
      sep = " - "
    ),
    as.character(
      safe_min(
        eigengene_pair_tbl$
          Signed_eigengene_dissimilarity
      )
    ),
    as.character(
      EIGENGENE_MERGE_COR_THRESHOLD
    ),
    as.character(
      sum(
        eigengene_pair_tbl$
          Exceeds_merge_correlation_threshold,
        na.rm = TRUE
      )
    )
  )
)

safe_write_csv(
  eigengene_cor_tbl,
  file.path(
    OUTDIR,
    "tables",
    "eigengenes",
    "module_eigengene_correlation_matrix.csv"
  )
)

safe_write_csv(
  eigengene_pair_tbl,
  file.path(
    OUTDIR,
    "tables",
    "eigengenes",
    "module_eigengene_pairwise_correlations.csv"
  )
)

safe_write_csv(
  eigengene_quality_summary,
  file.path(
    OUTDIR,
    "tables",
    "eigengenes",
    "module_eigengene_quality_summary.csv"
  )
)

###############################################################################
# 16) REPRODUCIBLE OFF-DIAGONAL EDGE SAMPLE
###############################################################################

set.seed(
  EDGE_SAMPLE_SEED
)

sample_i <- sample.int(
  n_genes,
  size =
    EDGE_SAMPLE_SIZE,
  replace = TRUE
)

sample_j <- sample.int(
  n_genes,
  size =
    EDGE_SAMPLE_SIZE,
  replace = TRUE
)

same_index <- sample_i ==
  sample_j

while (any(same_index)) {
  sample_j[
    same_index
  ] <- sample.int(
    n_genes,
    size =
      sum(same_index),
    replace = TRUE
  )

  same_index <- sample_i ==
    sample_j
}

adjacency_sample <- adjacency_mat[
  cbind(
    sample_i,
    sample_j
  )
]

adjacency_reverse_sample <-
  adjacency_mat[
    cbind(
      sample_j,
      sample_i
    )
  ]

tom_sample <- TOM[
  cbind(
    sample_i,
    sample_j
  )
]

tom_reverse_sample <- TOM[
  cbind(
    sample_j,
    sample_i
  )
]

edge_sample_tbl <- tibble::tibble(
  Gene_i = gene_ids[
    sample_i
  ],
  Gene_j = gene_ids[
    sample_j
  ],
  Module_i = mergedColors[
    sample_i
  ],
  Module_j = mergedColors[
    sample_j
  ],
  Same_module =
    Module_i ==
    Module_j,
  Adjacency =
    adjacency_sample,
  TOM = tom_sample
)

edge_distribution_quantiles <- dplyr::bind_rows(
  matrix_sample_quantiles(
    adjacency_sample,
    "Adjacency"
  ),
  matrix_sample_quantiles(
    tom_sample,
    "TOM"
  )
)

sampled_matrix_audit <- tibble::tibble(
  Metric = c(
    "Edge_sample_seed",
    "N_sampled_directed_offdiagonal_pairs",
    "Maximum_sampled_adjacency_asymmetry",
    "Mean_sampled_adjacency_asymmetry",
    "Maximum_sampled_TOM_asymmetry",
    "Mean_sampled_TOM_asymmetry",
    "Sampled_adjacency_TOM_Pearson",
    "Sampled_adjacency_TOM_Spearman",
    "Sampled_same_module_fraction"
  ),
  Value = c(
    EDGE_SAMPLE_SEED,
    EDGE_SAMPLE_SIZE,
    safe_max(
      abs(
        adjacency_sample -
          adjacency_reverse_sample
      )
    ),
    safe_mean(
      abs(
        adjacency_sample -
          adjacency_reverse_sample
      )
    ),
    safe_max(
      abs(
        tom_sample -
          tom_reverse_sample
      )
    ),
    safe_mean(
      abs(
        tom_sample -
          tom_reverse_sample
      )
    ),
    safe_cor(
      adjacency_sample,
      tom_sample,
      method = "pearson"
    ),
    safe_cor(
      adjacency_sample,
      tom_sample,
      method = "spearman"
    ),
    mean(
      edge_sample_tbl$
        Same_module
    )
  )
)

sampled_edge_group_summary <- edge_sample_tbl %>%
  dplyr::group_by(
    Same_module
  ) %>%
  dplyr::summarise(
    N_sampled_edges =
      dplyr::n(),
    Mean_adjacency =
      safe_mean(Adjacency),
    Median_adjacency =
      safe_median(
        Adjacency
      ),
    Q95_adjacency =
      safe_quantile(
        Adjacency,
        0.95
      ),
    Mean_TOM =
      safe_mean(TOM),
    Median_TOM =
      safe_median(TOM),
    Q95_TOM =
      safe_quantile(
        TOM,
        0.95
      ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Edge_class =
      ifelse(
        Same_module,
        "Within module",
        "Between modules"
      )
  )

safe_write_csv(
  edge_distribution_quantiles,
  file.path(
    OUTDIR,
    "tables",
    "input",
    "sampled_offdiagonal_edge_quantiles.csv"
  )
)

safe_write_csv(
  sampled_matrix_audit,
  file.path(
    OUTDIR,
    "tables",
    "input",
    "sampled_matrix_symmetry_and_concordance_audit.csv"
  )
)

safe_write_csv(
  sampled_edge_group_summary,
  file.path(
    OUTDIR,
    "tables",
    "modularity",
    "sampled_within_between_edge_summary.csv"
  )
)

###############################################################################
# 17) INTEGRATE DEFINITIVE STRUCTURAL PRESERVATION
###############################################################################

integrated_module_quality <- module_size_tbl %>%
  dplyr::select(
    Module,
    N_genes,
    Proportion_of_network,
    Biological_label
  ) %>%
  dplyr::left_join(
    module_separation_tbl %>%
      dplyr::select(
        Module,
        Mean_within_adjacency,
        Mean_external_adjacency,
        Adjacency_within_external_ratio,
        Mean_within_TOM,
        Mean_external_TOM,
        TOM_within_external_ratio,
        Modularity_Q_contribution,
        Conductance,
        Median_within_connectivity_fraction,
        Median_participation_coefficient,
        Median_abs_kME,
        kWithin_abs_kME_Spearman
      ),
    by = "Module"
  ) %>%
  dplyr::left_join(
    hub_overlap_summary,
    by = "Module"
  )

if (preservation_available) {
  integrated_module_quality <-
    integrated_module_quality %>%
    dplyr::left_join(
      preservation_summary %>%
        dplyr::select(
          Module,
          Country_minimum_Zsummary,
          Country_minimum_preservation_class,
          Site_minimum_Zsummary,
          Site_minimum_preservation_class,
          Worst_country,
          Worst_site_country,
          Worst_site_reference,
          Worst_site_test
        ),
      by = "Module"
    )
} else {
  integrated_module_quality <-
    integrated_module_quality %>%
    dplyr::mutate(
      Country_minimum_Zsummary =
        NA_real_,
      Country_minimum_preservation_class =
        NA_character_,
      Site_minimum_Zsummary =
        NA_real_,
      Site_minimum_preservation_class =
        NA_character_,
      Worst_country =
        NA_character_,
      Worst_site_country =
        NA_character_,
      Worst_site_reference =
        NA_character_,
      Worst_site_test =
        NA_character_
    )
}

integrated_module_quality <-
  integrated_module_quality %>%
  dplyr::mutate(
    Quality_interpretation_boundary = paste(
      "Module size, Q contribution, within-between separation, kME and",
      "preservation are complementary descriptive properties and are not",
      "combined into a single validation score."
    )
  ) %>%
  dplyr::arrange(
    dplyr::desc(
      Adjacency_within_external_ratio
    )
  )

safe_write_csv(
  integrated_module_quality,
  file.path(
    OUTDIR,
    "tables",
    "integration",
    "integrated_module_network_quality_summary.csv"
  )
)

###############################################################################
# 18) GLOBAL REVIEWER SUMMARY
###############################################################################

module_with_highest_q_contribution <-
  module_modularity_tbl$Module[[1]]

module_with_lowest_q_contribution <-
  module_modularity_tbl$Module[[
    nrow(module_modularity_tbl)
  ]]

module_with_highest_adj_separation <-
  module_separation_tbl$Module[[1]]

module_with_lowest_adj_separation <-
  module_separation_tbl$Module[[
    nrow(module_separation_tbl)
  ]]

reviewer_summary <- tibble::tibble(
  Item = c(
    "Input level",
    "Samples",
    "Gene-collapsed proteins",
    "Selected beta",
    "Signed scale-free R2 at beta",
    "Mean connectivity at beta",
    "Final modules",
    "Grey genes",
    "Largest module",
    "Largest module fraction",
    "Normalized module-size entropy",
    "Effective number of modules",
    "Post hoc weighted modularity Q",
    "Highest module Q contribution",
    "Lowest module Q contribution",
    "Highest adjacency within/external separation",
    "Lowest adjacency within/external separation",
    "Highest remaining eigengene correlation",
    "Eigengene pairs above merge threshold",
    "Preservation summary available",
    "Recommended interpretation"
  ),
  Value = c(
    "GENE-COLLAPSED, outcome-independent SOMAmer selection",
    as.character(
      nrow(datExpr_clean)
    ),
    as.character(
      ncol(datExpr_clean)
    ),
    as.character(softPower),
    as.character(
      selected_soft_row$
        Signed_scale_free_R2
    ),
    as.character(
      selected_soft_row$
        Mean_connectivity
    ),
    as.character(n_modules),
    as.character(
      grey_gene_count
    ),
    module_size_tbl$Module[[1]],
    as.character(
      module_size_tbl$
        Proportion_of_network[[1]]
    ),
    as.character(
      normalized_module_entropy
    ),
    as.character(
      effective_module_number
    ),
    as.character(
      weighted_modularity_Q
    ),
    module_with_highest_q_contribution,
    module_with_lowest_q_contribution,
    module_with_highest_adj_separation,
    module_with_lowest_adj_separation,
    as.character(
      highest_signed_pair$
        Eigengene_correlation
    ),
    as.character(
      sum(
        eigengene_pair_tbl$
          Exceeds_merge_correlation_threshold
      )
    ),
    as.character(
      preservation_available
    ),
    paste(
      "The network shows coordinated but not disconnected plasma",
      "co-expression structure. Q is descriptive and must be interpreted",
      "with kME, internal/external connectivity, biology and preservation."
    )
  )
)

safe_write_csv(
  reviewer_summary,
  file.path(
    OUTDIR,
    "tables",
    "wgcna_network_quality_reviewer_summary.csv"
  )
)

###############################################################################
# 19) MANUSCRIPT-READY WORDING
###############################################################################

fmt <- function(
    x,
    digits = 3
) {
  if (
    length(x) == 0 ||
    !is.finite(x)
  ) {
    return("NA")
  }

  format(
    round(
      x,
      digits
    ),
    nsmall = digits,
    trim = TRUE
  )
}

fmt1 <- function(x) {
  fmt(x, 1)
}

results_text <- paste0(
  "The outcome-independent signed co-expression network included ",
  ncol(datExpr_clean),
  " gene-collapsed plasma proteins from ",
  nrow(datExpr_clean),
  " participants and resolved ",
  n_modules,
  " modules without an unassigned grey component. ",
  "A soft-thresholding power of beta = ",
  softPower,
  " yielded a signed scale-free topology fit of R2 = ",
  fmt(
    selected_soft_row$
      Signed_scale_free_R2,
    3
  ),
  " while retaining a mean connectivity of ",
  fmt1(
    selected_soft_row$
      Mean_connectivity
  ),
  ". Module sizes ranged from ",
  min(module_size_tbl$N_genes),
  " to ",
  max(module_size_tbl$N_genes),
  " proteins, with the largest module comprising ",
  scales::percent(
    max(
      module_size_tbl$
        Proportion_of_network
    ),
    accuracy = 0.1
  ),
  " of the network. ",
  "As a post hoc descriptive measure, the weighted adjacency yielded ",
  "a Newman-Girvan modularity of Q = ",
  fmt(
    weighted_modularity_Q,
    3
  ),
  ". Because WGCNA did not optimize this quantity and the plasma network ",
  "remained densely connected, Q was interpreted together with ",
  "module-membership, internal-to-external connectivity, biological ",
  "annotation and structural-preservation diagnostics rather than as a ",
  "standalone validation criterion."
)

methods_text <- paste0(
  "Network-quality diagnostics were performed on the complete signed ",
  "gene-collapsed WGCNA adjacency and topological-overlap matrices. ",
  "Soft-threshold selection was summarized using signed scale-free topology ",
  "fit and mean connectivity across candidate powers. Post hoc weighted ",
  "Newman-Girvan modularity was calculated after analytically excluding ",
  "self-adjacency terms as Q = sum_s[e_ss/(2m) - (a_s/(2m))^2], where e_ss ",
  "denotes the directed within-module edge strength and a_s the total strength ",
  "of module s. For each module, mean within-module and external adjacency ",
  "and topological overlap, conductance, intramodular connectivity, ",
  "participation coefficient and eigengene-based module membership were ",
  "summarized. Off-diagonal matrix quantiles and symmetry were evaluated ",
  "using a reproducible sample of ",
  format(
    EDGE_SAMPLE_SIZE,
    big.mark = ",",
    scientific = FALSE
  ),
  " directed gene pairs to avoid materializing complete upper-triangle ",
  "vectors."
)

limitations_text <- paste0(
  "Post hoc modularity Q is a descriptive graph-partition metric and was not ",
  "used to select the WGCNA parameters or optimize module assignments. ",
  "No universal Q threshold was assumed for this dense weighted plasma ",
  "network. Modules should therefore be interpreted as coordinated plasma ",
  "co-expression structures supported jointly by network diagnostics, ",
  "biological enrichment, association analyses and internal structural ",
  "preservation, rather than as disconnected communities or tissue-specific ",
  "mechanisms."
)

manuscript_wording <- tibble::tibble(
  Section = c(
    "Results",
    "Methods",
    "Interpretation limitation"
  ),
  Suggested_text = c(
    results_text,
    methods_text,
    limitations_text
  )
)

safe_write_csv(
  manuscript_wording,
  file.path(
    OUTDIR,
    "tables",
    "script16_manuscript_ready_wording.csv"
  )
)

writeLines(
  results_text,
  file.path(
    OUTDIR,
    "manuscript_ready_WGCNA_quality_Results.txt"
  )
)

writeLines(
  methods_text,
  file.path(
    OUTDIR,
    "manuscript_ready_WGCNA_quality_Methods.txt"
  )
)

writeLines(
  limitations_text,
  file.path(
    OUTDIR,
    "manuscript_ready_WGCNA_quality_Interpretation.txt"
  )
)

###############################################################################
# 20) FIGURES — SOFT THRESHOLD
###############################################################################

p_soft_r2 <- ggplot(
  soft_scan_clean,
  aes(
    x = Power,
    y =
      Signed_scale_free_R2
  )
) +
  geom_hline(
    yintercept = 0.80,
    linetype = 2,
    colour = "grey60"
  ) +
  geom_hline(
    yintercept = 0.90,
    linetype = 3,
    colour = "grey60"
  ) +
  geom_line(
    linewidth = 0.7
  ) +
  geom_point(
    size = 2
  ) +
  geom_vline(
    xintercept = softPower,
    linetype = 2
  ) +
  geom_point(
    data =
      selected_soft_row,
    size = 3.5
  ) +
  labs(
    title =
      "Signed scale-free topology fit",
    subtitle = paste0(
      "Selected beta = ",
      softPower,
      "; signed R2 = ",
      fmt(
        selected_soft_row$
          Signed_scale_free_R2,
        3
      )
    ),
    x = "Soft-threshold power",
    y = "Signed scale-free topology R2"
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
    "soft_threshold",
    "signed_scale_free_topology_fit.pdf"
  ),
  p_soft_r2,
  width = 6,
  height = 4.5
)

ggsave(
  file.path(
    OUTDIR,
    "figures",
    "soft_threshold",
    "signed_scale_free_topology_fit.png"
  ),
  p_soft_r2,
  width = 6,
  height = 4.5,
  dpi = DPI
)

p_soft_connectivity <- ggplot(
  soft_scan_clean,
  aes(
    x = Power,
    y = Mean_connectivity
  )
) +
  geom_line(
    linewidth = 0.7
  ) +
  geom_point(
    size = 2
  ) +
  geom_vline(
    xintercept = softPower,
    linetype = 2
  ) +
  geom_point(
    data =
      selected_soft_row,
    size = 3.5
  ) +
  labs(
    title =
      "Mean network connectivity",
    subtitle = paste0(
      "Selected beta = ",
      softPower,
      "; mean connectivity = ",
      fmt1(
        selected_soft_row$
          Mean_connectivity
      )
    ),
    x = "Soft-threshold power",
    y = "Mean connectivity"
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
    "soft_threshold",
    "mean_connectivity_by_soft_power.pdf"
  ),
  p_soft_connectivity,
  width = 6,
  height = 4.5
)

ggsave(
  file.path(
    OUTDIR,
    "figures",
    "soft_threshold",
    "mean_connectivity_by_soft_power.png"
  ),
  p_soft_connectivity,
  width = 6,
  height = 4.5,
  dpi = DPI
)

###############################################################################
# 21) FIGURES — MODULE ARCHITECTURE AND MODULARITY
###############################################################################

module_order_size <- module_size_tbl %>%
  dplyr::arrange(N_genes) %>%
  dplyr::pull(Module)

p_module_sizes <- module_size_tbl %>%
  dplyr::mutate(
    Module = factor(
      Module,
      levels =
        module_order_size
    )
  ) %>%
  ggplot(
    aes(
      x = Module,
      y = N_genes,
      fill = Module
    )
  ) +
  geom_col(
    width = 0.75
  ) +
  coord_flip() +
  scale_fill_manual(
    values = get_module_colors(
      modules_use
    ),
    guide = "none"
  ) +
  geom_text(
    aes(
      label = paste0(
        scales::comma(
          N_genes
        ),
        " (",
        scales::percent(
          Proportion_of_network,
          accuracy = 0.1
        ),
        ")"
      )
    ),
    hjust = -0.05,
    size = 3.2
  ) +
  scale_y_continuous(
    labels = scales::comma,
    expand = expansion(
      mult = c(
        0,
        0.20
      )
    )
  ) +
  labs(
    title =
      "Final WGCNA module architecture",
    subtitle =
      "Outcome-independent gene-collapsed plasma network",
    x = NULL,
    y = "Number of proteins"
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
    "modules",
    "final_module_sizes_and_proportions.pdf"
  ),
  p_module_sizes,
  width = 7.5,
  height = 5
)

ggsave(
  file.path(
    OUTDIR,
    "figures",
    "modules",
    "final_module_sizes_and_proportions.png"
  ),
  p_module_sizes,
  width = 7.5,
  height = 5,
  dpi = DPI
)

p_q_contribution <- module_modularity_tbl %>%
  dplyr::mutate(
    Module = factor(
      Module,
      levels = rev(
        Module
      )
    )
  ) %>%
  ggplot(
    aes(
      x = Module,
      y =
        Modularity_Q_contribution,
      fill = Module
    )
  ) +
  geom_hline(
    yintercept = 0,
    colour = "grey50"
  ) +
  geom_col(
    width = 0.75
  ) +
  coord_flip() +
  scale_fill_manual(
    values = get_module_colors(
      modules_use
    ),
    guide = "none"
  ) +
  labs(
    title =
      "Module contributions to post hoc weighted modularity",
    subtitle = paste0(
      "Total Q = ",
      fmt(
        weighted_modularity_Q,
        3
      ),
      "; Q was not a WGCNA optimization target"
    ),
    x = NULL,
    y = "Contribution to weighted modularity Q"
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
    "modularity",
    "module_modularity_Q_contributions.pdf"
  ),
  p_q_contribution,
  width = 7,
  height = 5
)

ggsave(
  file.path(
    OUTDIR,
    "figures",
    "modularity",
    "module_modularity_Q_contributions.png"
  ),
  p_q_contribution,
  width = 7,
  height = 5,
  dpi = DPI
)

###############################################################################
# 22) FIGURES — BLOCK MIXING
###############################################################################

p_adjacency_mixing <- adjacency_block_tbl %>%
  dplyr::mutate(
    Module_row = factor(
      Module_row,
      levels = rev(
        modules_use
      )
    ),
    Module_column = factor(
      Module_column,
      levels = modules_use
    )
  ) %>%
  ggplot(
    aes(
      x = Module_column,
      y = Module_row,
      fill = Mean_edge_weight
    )
  ) +
  geom_tile(
    colour = "white",
    linewidth = 0.4
  ) +
  geom_text(
    aes(
      label = sprintf(
        "%.3f",
        Mean_edge_weight
      )
    ),
    size = 2.8
  ) +
  scale_fill_gradient(
    low = "white",
    high = "#1F4E79",
    name = "Mean\nadjacency"
  ) +
  labs(
    title =
      "Weighted adjacency mixing across modules",
    subtitle =
      "Diagonal cells represent within-module mean adjacency",
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
    plot.title = element_text(
      face = "bold"
    )
  )

ggsave(
  file.path(
    OUTDIR,
    "figures",
    "modularity",
    "adjacency_module_mixing_heatmap.pdf"
  ),
  p_adjacency_mixing,
  width = 7.5,
  height = 6.5
)

ggsave(
  file.path(
    OUTDIR,
    "figures",
    "modularity",
    "adjacency_module_mixing_heatmap.png"
  ),
  p_adjacency_mixing,
  width = 7.5,
  height = 6.5,
  dpi = DPI
)

p_tom_mixing <- tom_block_tbl %>%
  dplyr::mutate(
    Module_row = factor(
      Module_row,
      levels = rev(
        modules_use
      )
    ),
    Module_column = factor(
      Module_column,
      levels = modules_use
    )
  ) %>%
  ggplot(
    aes(
      x = Module_column,
      y = Module_row,
      fill = Mean_edge_weight
    )
  ) +
  geom_tile(
    colour = "white",
    linewidth = 0.4
  ) +
  geom_text(
    aes(
      label = sprintf(
        "%.3f",
        Mean_edge_weight
      )
    ),
    size = 2.8
  ) +
  scale_fill_gradient(
    low = "white",
    high = "#7A3E9D",
    name = "Mean\nTOM"
  ) +
  labs(
    title =
      "Topological-overlap mixing across modules",
    subtitle =
      "Diagonal cells represent within-module mean topological overlap",
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
    plot.title = element_text(
      face = "bold"
    )
  )

ggsave(
  file.path(
    OUTDIR,
    "figures",
    "modularity",
    "TOM_module_mixing_heatmap.pdf"
  ),
  p_tom_mixing,
  width = 7.5,
  height = 6.5
)

ggsave(
  file.path(
    OUTDIR,
    "figures",
    "modularity",
    "TOM_module_mixing_heatmap.png"
  ),
  p_tom_mixing,
  width = 7.5,
  height = 6.5,
  dpi = DPI
)

###############################################################################
# 23) FIGURES — CONNECTIVITY AND kME
###############################################################################

p_kme <- node_quality_tbl %>%
  dplyr::mutate(
    Module = factor(
      Module,
      levels = modules_use
    )
  ) %>%
  ggplot(
    aes(
      x = Module,
      y = Abs_assigned_kME,
      fill = Module
    )
  ) +
  geom_boxplot(
    outlier.alpha = 0.15,
    width = 0.70
  ) +
  scale_fill_manual(
    values = get_module_colors(
      modules_use
    ),
    guide = "none"
  ) +
  labs(
    title =
      "Assigned module-membership distributions",
    subtitle =
      "Absolute correlation between each protein and its assigned eigengene",
    x = NULL,
    y = "Absolute assigned kME"
  ) +
  theme_bw(
    base_size = 11
  ) +
  theme(
    axis.text.x = element_text(
      angle = 35,
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
    "connectivity",
    "assigned_absolute_kME_by_module.pdf"
  ),
  p_kme,
  width = 7.5,
  height = 5
)

ggsave(
  file.path(
    OUTDIR,
    "figures",
    "connectivity",
    "assigned_absolute_kME_by_module.png"
  ),
  p_kme,
  width = 7.5,
  height = 5,
  dpi = DPI
)

p_connectivity_kme <- ggplot(
  node_quality_tbl,
  aes(
    x = kWithin_adjacency,
    y = Abs_assigned_kME,
    colour = Module
  )
) +
  geom_point(
    alpha = 0.22,
    size = 0.7
  ) +
  facet_wrap(
    ~ Module,
    scales = "free_x"
  ) +
  scale_colour_manual(
    values = get_module_colors(
      modules_use
    ),
    guide = "none"
  ) +
  labs(
    title =
      "Intramodular connectivity and eigengene membership",
    subtitle =
      "Relationship shown separately for each final module",
    x = "Intramodular adjacency connectivity",
    y = "Absolute assigned kME"
  ) +
  theme_bw(
    base_size = 10
  ) +
  theme(
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
    "connectivity",
    "intramodular_connectivity_vs_absolute_kME.pdf"
  ),
  p_connectivity_kme,
  width = 10,
  height = 7
)

ggsave(
  file.path(
    OUTDIR,
    "figures",
    "connectivity",
    "intramodular_connectivity_vs_absolute_kME.png"
  ),
  p_connectivity_kme,
  width = 10,
  height = 7,
  dpi = DPI
)

separation_plot_tbl <- module_separation_tbl %>%
  dplyr::select(
    Module,
    `Adjacency ratio` =
      Adjacency_within_external_ratio,
    `TOM ratio` =
      TOM_within_external_ratio
  ) %>%
  tidyr::pivot_longer(
    cols = c(
      `Adjacency ratio`,
      `TOM ratio`
    ),
    names_to = "Metric",
    values_to = "Within_external_ratio"
  ) %>%
  dplyr::mutate(
    Module = factor(
      Module,
      levels = modules_use
    )
  )

p_separation <- ggplot(
  separation_plot_tbl,
  aes(
    x = Module,
    y = Within_external_ratio,
    fill = Metric
  )
) +
  geom_hline(
    yintercept = 1,
    linetype = 2,
    colour = "grey50"
  ) +
  geom_col(
    position =
      position_dodge(
        width = 0.8
      ),
    width = 0.72
  ) +
  labs(
    title =
      "Within-module versus external network separation",
    subtitle =
      "Values above one indicate greater average within-module connectivity",
    x = NULL,
    y = "Within/external mean edge-weight ratio",
    fill = NULL
  ) +
  theme_bw(
    base_size = 11
  ) +
  theme(
    axis.text.x = element_text(
      angle = 35,
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
    "modularity",
    "module_within_external_separation_ratios.pdf"
  ),
  p_separation,
  width = 8,
  height = 5
)

ggsave(
  file.path(
    OUTDIR,
    "figures",
    "modularity",
    "module_within_external_separation_ratios.png"
  ),
  p_separation,
  width = 8,
  height = 5,
  dpi = DPI
)

###############################################################################
# 24) FIGURES — EIGENGENE CORRELATION
###############################################################################

eigengene_plot_tbl <- eigengene_cor_tbl %>%
  tidyr::pivot_longer(
    cols = -Module,
    names_to = "Module_2",
    values_to =
      "Correlation"
  ) %>%
  dplyr::mutate(
    Module = factor(
      Module,
      levels = rev(
        modules_use
      )
    ),
    Module_2 = factor(
      Module_2,
      levels = modules_use
    )
  )

p_eigengene <- ggplot(
  eigengene_plot_tbl,
  aes(
    x = Module_2,
    y = Module,
    fill = Correlation
  )
) +
  geom_tile(
    colour = "white",
    linewidth = 0.45
  ) +
  geom_text(
    aes(
      label = sprintf(
        "%.2f",
        Correlation
      )
    ),
    size = 2.9
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1),
    oob = scales::squish,
    name = "Eigengene\ncorrelation"
  ) +
  labs(
    title =
      "Final module-eigengene correlations",
    subtitle = paste0(
      "Merge-correlation reference threshold = ",
      EIGENGENE_MERGE_COR_THRESHOLD
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
    plot.title = element_text(
      face = "bold"
    )
  )

ggsave(
  file.path(
    OUTDIR,
    "figures",
    "eigengenes",
    "module_eigengene_correlation_heatmap.pdf"
  ),
  p_eigengene,
  width = 7.5,
  height = 6.5
)

ggsave(
  file.path(
    OUTDIR,
    "figures",
    "eigengenes",
    "module_eigengene_correlation_heatmap.png"
  ),
  p_eigengene,
  width = 7.5,
  height = 6.5,
  dpi = DPI
)

###############################################################################
# 25) FIGURE — INTEGRATED MODULE QUALITY
###############################################################################

integrated_plot_tbl <- integrated_module_quality %>%
  dplyr::select(
    Module,
    `Adjacency within/external` =
      Adjacency_within_external_ratio,
    `TOM within/external` =
      TOM_within_external_ratio,
    `Median abs kME` =
      Median_abs_kME,
    `Within-connectivity fraction` =
      Median_within_connectivity_fraction,
    `Country minimum Zsummary` =
      Country_minimum_Zsummary,
    `Site minimum Zsummary` =
      Site_minimum_Zsummary
  ) %>%
  tidyr::pivot_longer(
    cols = -Module,
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  dplyr::group_by(Metric) %>%
  dplyr::mutate(
    Metric_z = if (
      safe_sd(Value) > 0
    ) {
      (
        Value -
          safe_mean(Value)
      ) /
        safe_sd(Value)
    } else {
      0
    }
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    Module = factor(
      Module,
      levels = rev(
        modules_use
      )
    )
  )

p_integrated <- ggplot(
  integrated_plot_tbl,
  aes(
    x = Metric,
    y = Module,
    fill = Metric_z
  )
) +
  geom_tile(
    colour = "white",
    linewidth = 0.45
  ) +
  geom_text(
    aes(
      label = ifelse(
        is.finite(Value),
        sprintf(
          "%.2f",
          Value
        ),
        ""
      )
    ),
    size = 2.7
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    name = "Within-metric\nstandardized value"
  ) +
  labs(
    title =
      "Complementary module-quality diagnostics",
    subtitle = paste(
      "Cell labels show raw values; colors are standardized separately",
      "within each metric and do not form a composite score"
    ),
    x = NULL,
    y = NULL
  ) +
  theme_bw(
    base_size = 10
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
    plot.title = element_text(
      face = "bold"
    )
  )

ggsave(
  file.path(
    OUTDIR,
    "figures",
    "integration",
    "integrated_module_quality_diagnostics_heatmap.pdf"
  ),
  p_integrated,
  width = 11,
  height = 6
)

ggsave(
  file.path(
    OUTDIR,
    "figures",
    "integration",
    "integrated_module_quality_diagnostics_heatmap.png"
  ),
  p_integrated,
  width = 11,
  height = 6,
  dpi = DPI
)

###############################################################################
# 26) EXCEL WORKBOOK
###############################################################################

workbook_tables <- list(
  Input_audit =
    input_alignment_audit,
  Reviewer_summary =
    reviewer_summary,
  Soft_threshold =
    soft_threshold_summary,
  Soft_scan =
    soft_scan_clean,
  Module_sizes =
    module_size_tbl,
  Architecture_summary =
    module_architecture_summary,
  Global_modularity_Q =
    global_modularity_tbl,
  Q_by_module =
    module_modularity_tbl,
  Module_separation =
    module_separation_tbl,
  Module_connectivity =
    module_connectivity_summary,
  Module_kME =
    module_kME_summary,
  Hub_overlap =
    hub_overlap_summary,
  Eigengene_summary =
    eigengene_quality_summary,
  Eigengene_pairs =
    eigengene_pair_tbl,
  Edge_quantiles =
    edge_distribution_quantiles,
  Sampled_matrix_audit =
    sampled_matrix_audit,
  Integrated_quality =
    integrated_module_quality,
  Manuscript_wording =
    manuscript_wording
)

write_workbook_safe(
  workbook_tables,
  file.path(
    OUTDIR,
    "WGCNA_16_Network_Quality_Modularity_Diagnostics.xlsx"
  )
)

###############################################################################
# 27) FINAL SUMMARY AND OUTPUT MANIFEST
###############################################################################

script16_summary <- tibble::tibble(
  Metric = c(
    "Base_directory",
    "Output_directory",
    "Input_level",
    "N_samples",
    "N_genes",
    "Soft_power",
    "Signed_scale_free_R2",
    "Mean_connectivity",
    "N_modules",
    "Grey_genes",
    "Largest_module",
    "Largest_module_fraction",
    "Weighted_modularity_Q",
    "Maximum_positive_eigengene_correlation",
    "Eigengene_pairs_above_merge_threshold",
    "Edge_sample_size",
    "Preservation_summary_available",
    "Country_numeric_included"
  ),
  Value = c(
    BASE_DIR,
    OUTDIR,
    "GENE-COLLAPSED, outcome-independent SOMAmer selection",
    as.character(
      nrow(datExpr_clean)
    ),
    as.character(
      ncol(datExpr_clean)
    ),
    as.character(softPower),
    as.character(
      selected_soft_row$
        Signed_scale_free_R2
    ),
    as.character(
      selected_soft_row$
        Mean_connectivity
    ),
    as.character(n_modules),
    as.character(
      grey_gene_count
    ),
    module_size_tbl$Module[[1]],
    as.character(
      module_size_tbl$
        Proportion_of_network[[1]]
    ),
    as.character(
      weighted_modularity_Q
    ),
    as.character(
      highest_signed_pair$
        Eigengene_correlation
    ),
    as.character(
      sum(
        eigengene_pair_tbl$
          Exceeds_merge_correlation_threshold
      )
    ),
    as.character(
      EDGE_SAMPLE_SIZE
    ),
    as.character(
      preservation_available
    ),
    "FALSE"
  )
)

safe_write_csv(
  script16_summary,
  file.path(
    OUTDIR,
    "tables",
    "script16_final_summary.csv"
  )
)

output_manifest <- tibble::tibble(
  Output_file = c(
    "tables/input/script16_input_alignment_audit.csv",
    "tables/input/script11_network_matrix_alignment_audit.csv",
    "tables/soft_threshold/soft_threshold_scan_clean.csv",
    "tables/soft_threshold/soft_threshold_quality_summary.csv",
    "tables/modularity/final_module_sizes_and_proportions.csv",
    "tables/modularity/module_architecture_concentration_summary.csv",
    "tables/modularity/posthoc_weighted_modularity_Q.csv",
    "tables/modularity/module_specific_modularity_contributions.csv",
    "tables/modularity/module_internal_external_separation_summary.csv",
    "tables/matrix_mixing/adjacency_module_block_mixing.csv",
    "tables/matrix_mixing/TOM_module_block_mixing.csv",
    "tables/connectivity/gene_level_network_quality_metrics.csv",
    "tables/connectivity/module_connectivity_quality_summary.csv",
    "tables/kME/module_kME_quality_summary.csv",
    "tables/kME/top_hubs_by_intramodular_connectivity.csv",
    "tables/kME/top_hubs_by_absolute_kME.csv",
    "tables/kME/hub_definition_overlap_summary.csv",
    "tables/eigengenes/module_eigengene_correlation_matrix.csv",
    "tables/eigengenes/module_eigengene_pairwise_correlations.csv",
    "tables/integration/integrated_module_network_quality_summary.csv",
    "tables/wgcna_network_quality_reviewer_summary.csv",
    "tables/script16_manuscript_ready_wording.csv",
    "figures/soft_threshold/signed_scale_free_topology_fit.pdf/png",
    "figures/soft_threshold/mean_connectivity_by_soft_power.pdf/png",
    "figures/modules/final_module_sizes_and_proportions.pdf/png",
    "figures/modularity/module_modularity_Q_contributions.pdf/png",
    "figures/modularity/adjacency_module_mixing_heatmap.pdf/png",
    "figures/modularity/TOM_module_mixing_heatmap.pdf/png",
    "figures/modularity/module_within_external_separation_ratios.pdf/png",
    "figures/connectivity/assigned_absolute_kME_by_module.pdf/png",
    "figures/connectivity/intramodular_connectivity_vs_absolute_kME.pdf/png",
    "figures/eigengenes/module_eigengene_correlation_heatmap.pdf/png",
    "figures/integration/integrated_module_quality_diagnostics_heatmap.pdf/png",
    "WGCNA_16_Network_Quality_Modularity_Diagnostics.xlsx",
    "workspace/script16_network_quality_modularity_workspace.RData",
    "sessionInfo.txt"
  ),
  Description = c(
    "Strict expression, module, eigengene, adjacency and TOM alignment audit.",
    "Network-matrix checks inherited from Script 11.",
    "Clean candidate-power scale-free and connectivity diagnostics.",
    "Selected-power quality summary.",
    "Final module sizes, fractions and biological labels.",
    "Module-size concentration and entropy diagnostics.",
    "Exact post hoc weighted Newman-Girvan modularity Q.",
    "Module-specific contributions to Q and conductance.",
    "Within-module versus external adjacency and TOM separation.",
    "Complete 8 x 8 adjacency mixing table.",
    "Complete 8 x 8 topological-overlap mixing table.",
    "Per-gene connectivity, participation, kME and rank diagnostics.",
    "Connectivity distributions summarized by module.",
    "Assigned kME distributions and connectivity concordance.",
    "Top intramodular-connectivity hubs.",
    "Top eigengene-membership hubs.",
    "Overlap between the two hub definitions.",
    "Final module-eigengene correlation matrix.",
    "Pairwise eigengene similarity and merge-threshold audit.",
    "Integrated network quality and fixed-gene preservation summary.",
    "Reviewer-facing global quality summary.",
    "Methods-, Results- and interpretation-ready wording.",
    "Signed scale-free topology diagnostic.",
    "Mean-connectivity diagnostic.",
    "Module architecture figure.",
    "Module contributions to post hoc Q.",
    "Adjacency module-mixing heatmap.",
    "Topological-overlap module-mixing heatmap.",
    "Within/external network-separation ratios.",
    "Absolute assigned-kME distributions.",
    "Intramodular connectivity versus assigned kME.",
    "Final eigengene-correlation heatmap.",
    "Complementary non-composite module-quality heatmap.",
    "Integrated reviewer-facing workbook.",
    "Reduced workspace excluding dense adjacency and TOM matrices.",
    "R session information."
  )
)

safe_write_csv(
  output_manifest,
  file.path(
    OUTDIR,
    "script16_output_manifest.csv"
  )
)

###############################################################################
# 28) SAVE REDUCED WORKSPACE
###############################################################################

# Dense adjacency and TOM matrices are intentionally excluded.
save(
  soft_scan_clean,
  selected_soft_row,
  soft_threshold_summary,
  module_size_tbl,
  module_architecture_summary,
  adjacency_block_tbl,
  tom_block_tbl,
  global_modularity_tbl,
  module_modularity_tbl,
  node_quality_tbl,
  module_connectivity_summary,
  module_kME_summary,
  top_hubs_by_kWithin,
  top_hubs_by_kME,
  hub_overlap_summary,
  module_separation_tbl,
  eigengene_cor,
  eigengene_cor_tbl,
  eigengene_pair_tbl,
  eigengene_quality_summary,
  edge_distribution_quantiles,
  sampled_matrix_audit,
  sampled_edge_group_summary,
  integrated_module_quality,
  reviewer_summary,
  manuscript_wording,
  script16_summary,
  output_manifest,
  file = file.path(
    OUTDIR,
    "workspace",
    "script16_network_quality_modularity_workspace.RData"
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
# 29) FINAL MESSAGE
###############################################################################

cat("\nScript 16 finished successfully.\n")
cat("Main output directory:\n", OUTDIR, "\n")
cat("Input level: GENE-COLLAPSED, outcome-independent SOMAmer selection.\n")
cat("Samples: ", nrow(datExpr_clean), "\n", sep = "")
cat("Genes: ", ncol(datExpr_clean), "\n", sep = "")
cat("Selected beta: ", softPower, "\n", sep = "")
cat(
  "Signed scale-free R2: ",
  selected_soft_row$
    Signed_scale_free_R2,
  "\n",
  sep = ""
)
cat(
  "Mean connectivity: ",
  selected_soft_row$
    Mean_connectivity,
  "\n",
  sep = ""
)
cat("Final modules: ", n_modules, "\n", sep = "")
cat("Grey genes: ", grey_gene_count, "\n", sep = "")
cat(
  "Largest module: ",
  module_size_tbl$Module[[1]],
  " (",
  module_size_tbl$N_genes[[1]],
  "; ",
  scales::percent(
    module_size_tbl$
      Proportion_of_network[[1]],
    accuracy = 0.1
  ),
  ")\n",
  sep = ""
)
cat(
  "Post hoc weighted modularity Q: ",
  weighted_modularity_Q,
  "\n",
  sep = ""
)
cat(
  "Highest remaining positive eigengene correlation: ",
  highest_signed_pair$
    Eigengene_correlation,
  " (",
  highest_signed_pair$Module_1,
  " - ",
  highest_signed_pair$Module_2,
  ")\n",
  sep = ""
)
cat(
  "Eigengene pairs above merge threshold: ",
  sum(
    eigengene_pair_tbl$
      Exceeds_merge_correlation_threshold
  ),
  "\n",
  sep = ""
)
cat(
  "Definitive Script 15b preservation integrated: ",
  preservation_available,
  "\n",
  sep = ""
)
cat("Country_numeric included: FALSE\n")
cat("\nKey outputs:\n")
cat("- tables/wgcna_network_quality_reviewer_summary.csv\n")
cat("- tables/modularity/posthoc_weighted_modularity_Q.csv\n")
cat("- tables/modularity/module_internal_external_separation_summary.csv\n")
cat("- tables/connectivity/module_connectivity_quality_summary.csv\n")
cat("- tables/kME/module_kME_quality_summary.csv\n")
cat("- tables/eigengenes/module_eigengene_quality_summary.csv\n")
cat("- tables/integration/integrated_module_network_quality_summary.csv\n")
cat("- WGCNA_16_Network_Quality_Modularity_Diagnostics.xlsx\n")
cat("\nDense matrices were not copied to the Script 16 workspace.\n")

###############################################################################
# END
###############################################################################


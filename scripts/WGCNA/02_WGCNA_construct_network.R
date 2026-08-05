###############################################################################
# ReDLat plasma proteomics — WGCNA workflow
# 02. Construct the WGCNA network
# Requires: outputs from Script 01
# Produces: modules, eigengenes, network diagnostics and workspaces
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
  "WGCNA",
  "dynamicTreeCut",
  "dplyr",
  "tibble",
  "readr",
  "ggplot2",
  "stringr"
)

cran_missing <- cran_pkgs[
  !vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(cran_missing) > 0L) {
  stop("Missing required packages: ", paste(cran_missing, collapse = ", "),
       ". Run renv::restore() before this script.", call. = FALSE)
}

invisible(lapply(cran_pkgs, library, character.only = TRUE))

options(stringsAsFactors = FALSE)
options(error = traceback)

WGCNA::allowWGCNAThreads()

###############################################################################
# 2) PATHS
###############################################################################

base_dir <- WGCNA_CONFIG$project_root

input_dir <- file.path(WGCNA_CONFIG$result_root,
  "01_input"
)

outdir <- file.path(WGCNA_CONFIG$result_root,
  "02_network"
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

subdirs <- c(
  "qc",
  "soft_threshold",
  "network",
  "modules",
  "eigengenes",
  "tables",
  "workspace"
)

invisible(lapply(
  file.path(outdir, subdirs),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

expr_file <- file.path(
  input_dir,
  "gene_collapsed_expression_matrix.csv"
)

meta_file <- file.path(
  input_dir,
  "wgcna_sample_metadata.csv"
)

annot_file <- file.path(
  input_dir,
  "gene_level_annotation.csv"
)

manifest_file <- file.path(
  input_dir,
  "wgcna_input_manifest.csv"
)

input_qc_file <- file.path(
  input_dir,
  "wgcna_input_qc_summary.csv"
)

gene_map_audit_file <- file.path(
  input_dir,
  "qc",
  "gene_map_audit.csv"
)

expression_scale_audit_file <- file.path(
  input_dir,
  "qc",
  "expression_scale_audit.csv"
)

required_files <- c(
  expr_file,
  meta_file,
  annot_file,
  input_qc_file,
  gene_map_audit_file,
  expression_scale_audit_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing required files from Script 10:\n",
    paste(missing_files, collapse = "\n"),
    "\n\nRun 01_WGCNA_prepare_input.R first.",
    call. = FALSE
  )
}

###############################################################################
# 3) GLOBAL PARAMETERS
###############################################################################

# Original network parameters retained.
DEFAULT_SOFT_POWER <- 4
MIN_MODULE_SIZE    <- 30
MERGE_CUT_HEIGHT   <- 0.25
DEEPSPLIT          <- 3
NETWORK_TYPE       <- "signed"

# If manual outlier removal is not desired, keep NA.
OUTLIER_CUT_HEIGHT <- NA_real_

# The original analysis used beta = 4 rather than automatic reselection.
AUTO_SELECT_SOFT_POWER <- FALSE
TARGET_SCALE_FREE_R2 <- 0.90

# Strict input audits based on the validated Script 10 output.
EXPECTED_N_SAMPLES <- 639L
EXPECTED_N_GENES <- 9638L
REQUIRE_EXACT_INPUT_DIMENSIONS <- TRUE
REQUIRE_NO_NONFINITE_VALUES <- TRUE
REQUIRE_GSG_TO_PRESERVE_ALL_INPUT <- TRUE

# Optional checkpoint after TOM construction. This can create a very large file.
SAVE_NETWORK_CHECKPOINT_AFTER_TOM <- FALSE

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

make_display_name <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- "Unknown_protein"
  make.unique(x, sep = "_dup")
}

read_optional_csv <- function(file) {
  if (!file.exists(file)) return(NULL)
  readr::read_csv(file, show_col_types = FALSE)
}

safe_pick_value <- function(df, row_filter, candidates) {
  if (is.null(df) || nrow(df) == 0) return(NA_real_)

  hit <- candidates[candidates %in% names(df)][1]
  if (length(hit) == 0 || is.na(hit)) return(NA_real_)

  selected <- df[row_filter, , drop = FALSE]
  if (nrow(selected) == 0) return(NA_real_)

  safe_numeric(selected[[hit]][1])
}

select_soft_power <- function(sft_fit_indices,
                              default_power = 4,
                              target_r2 = 0.90,
                              auto_select = FALSE) {
  if (!auto_select) {
    return(default_power)
  }

  fit_df <- as.data.frame(sft_fit_indices)

  if (!all(c("Power", "SFT.R.sq", "slope") %in% names(fit_df))) {
    return(default_power)
  }

  fit_df <- fit_df %>%
    dplyr::mutate(
      signed_R2 = -sign(slope) * SFT.R.sq
    )

  candidate <- fit_df %>%
    dplyr::filter(
      is.finite(signed_R2),
      signed_R2 >= target_r2
    ) %>%
    dplyr::arrange(Power) %>%
    dplyr::slice(1)

  if (nrow(candidate) == 1) {
    return(candidate$Power[[1]])
  }

  default_power
}

estimate_network_memory <- function(n_genes) {
  bytes_per_dense_matrix <- as.numeric(n_genes)^2 * 8
  gib_per_dense_matrix <- bytes_per_dense_matrix / 1024^3

  tibble::tibble(
    n_genes = n_genes,
    estimated_GiB_per_dense_numeric_matrix = gib_per_dense_matrix,
    estimated_GiB_for_adjacency_TOM_dissTOM = gib_per_dense_matrix * 3,
    note = paste(
      "Approximate lower-bound memory estimate.",
      "R object overhead, temporary copies and clustering objects require",
      "additional memory beyond these values."
    )
  )
}

assert_exact_set <- function(x, y, label_x, label_y) {
  missing_in_y <- setdiff(x, y)
  missing_in_x <- setdiff(y, x)

  if (length(missing_in_y) > 0 || length(missing_in_x) > 0) {
    stop(
      "Feature-set mismatch between ", label_x, " and ", label_y, ".\n",
      "Missing in ", label_y, ": ", length(missing_in_y), "\n",
      "Missing in ", label_x, ": ", length(missing_in_x),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

###############################################################################
# 5) LOAD INPUTS FROM SCRIPT 10
###############################################################################

expr_gene_df <- readr::read_csv(
  expr_file,
  show_col_types = FALSE,
  guess_max = 100000
)

meta_df <- readr::read_csv(
  meta_file,
  show_col_types = FALSE,
  guess_max = 100000
)

gene_annotation <- readr::read_csv(
  annot_file,
  show_col_types = FALSE,
  guess_max = 100000
)

manifest_df <- read_optional_csv(manifest_file)
input_qc_df <- read_optional_csv(input_qc_file)
gene_map_audit_df <- read_optional_csv(gene_map_audit_file)
expression_scale_audit_df <- read_optional_csv(expression_scale_audit_file)

cat("INPUT DIMENSIONS\n")
cat("expr_gene_df    :", dim(expr_gene_df), "\n")
cat("meta_df         :", dim(meta_df), "\n")
cat("gene_annotation :", dim(gene_annotation), "\n\n")

cat("Input directory:\n")
cat(input_dir, "\n\n")

###############################################################################
# 6) STRICT BASIC INPUT CHECKS
###############################################################################

if (!"SampleId" %in% names(expr_gene_df)) {
  stop("expr_gene_df must contain SampleId.", call. = FALSE)
}

if (!"SampleId" %in% names(meta_df)) {
  stop("meta_df must contain SampleId.", call. = FALSE)
}

if (!"EntrezGeneSymbol" %in% names(gene_annotation)) {
  stop("gene_annotation must contain EntrezGeneSymbol.", call. = FALSE)
}

if (anyDuplicated(expr_gene_df$SampleId) > 0) {
  stop("Duplicated SampleId values were found in the expression matrix.", call. = FALSE)
}

if (anyDuplicated(meta_df$SampleId) > 0) {
  stop("Duplicated SampleId values were found in metadata.", call. = FALSE)
}

expr_gene_df <- expr_gene_df %>%
  dplyr::mutate(
    SampleId = as.character(SampleId)
  )

meta_df <- meta_df %>%
  dplyr::mutate(
    SampleId = as.character(SampleId)
  )

gene_annotation <- gene_annotation %>%
  dplyr::mutate(
    EntrezGeneSymbol = clean_text_na(EntrezGeneSymbol),
    Protein_Name = if ("Protein_Name" %in% names(.)) {
      clean_text_na(Protein_Name)
    } else {
      NA_character_
    },
    AptName = if ("AptName" %in% names(.)) {
      clean_text_na(AptName)
    } else {
      NA_character_
    }
  ) %>%
  dplyr::filter(
    !is.na(EntrezGeneSymbol),
    EntrezGeneSymbol != ""
  )

if (anyDuplicated(gene_annotation$EntrezGeneSymbol) > 0) {
  stop("Duplicated EntrezGeneSymbol values were found in gene annotation.", call. = FALSE)
}

if ("AptName" %in% names(gene_annotation) &&
    anyDuplicated(gene_annotation$AptName) > 0) {
  stop("Duplicated AptName values were found in gene annotation.", call. = FALSE)
}

expression_genes <- setdiff(
  names(expr_gene_df),
  "SampleId"
)

if (anyDuplicated(expression_genes) > 0) {
  stop("Duplicated gene columns were found in the expression matrix.", call. = FALSE)
}

assert_exact_set(
  expression_genes,
  gene_annotation$EntrezGeneSymbol,
  "expression matrix",
  "gene annotation"
)

assert_exact_set(
  expr_gene_df$SampleId,
  meta_df$SampleId,
  "expression matrix samples",
  "metadata samples"
)

if (REQUIRE_EXACT_INPUT_DIMENSIONS) {
  if (nrow(expr_gene_df) != EXPECTED_N_SAMPLES) {
    stop(
      "Expected ", EXPECTED_N_SAMPLES,
      " samples, but expression input contains ", nrow(expr_gene_df), ".",
      call. = FALSE
    )
  }

  if (length(expression_genes) != EXPECTED_N_GENES) {
    stop(
      "Expected ", EXPECTED_N_GENES,
      " genes, but expression input contains ", length(expression_genes), ".",
      call. = FALSE
    )
  }

  if (nrow(gene_annotation) != EXPECTED_N_GENES) {
    stop(
      "Expected ", EXPECTED_N_GENES,
      " annotation rows, but found ", nrow(gene_annotation), ".",
      call. = FALSE
    )
  }
}

# Reorder metadata and annotation exactly to the expression input.
meta_df <- meta_df[
  match(expr_gene_df$SampleId, meta_df$SampleId),
  ,
  drop = FALSE
]

gene_annotation <- gene_annotation[
  match(expression_genes, gene_annotation$EntrezGeneSymbol),
  ,
  drop = FALSE
]

if (!all(meta_df$SampleId == expr_gene_df$SampleId)) {
  stop("Sample order could not be aligned exactly.", call. = FALSE)
}

if (!all(gene_annotation$EntrezGeneSymbol == expression_genes)) {
  stop("Gene annotation order could not be aligned exactly.", call. = FALSE)
}

input_alignment_audit <- tibble::tibble(
  check = c(
    "expression_sample_ids_unique",
    "metadata_sample_ids_unique",
    "expression_and_metadata_sample_sets_identical",
    "expression_and_metadata_sample_order_identical",
    "expression_gene_columns_unique",
    "annotation_gene_symbols_unique",
    "expression_and_annotation_gene_sets_identical",
    "expression_and_annotation_gene_order_identical",
    "annotation_AptName_unique",
    "expected_sample_count",
    "expected_gene_count"
  ),
  passed = c(
    anyDuplicated(expr_gene_df$SampleId) == 0,
    anyDuplicated(meta_df$SampleId) == 0,
    setequal(expr_gene_df$SampleId, meta_df$SampleId),
    all(expr_gene_df$SampleId == meta_df$SampleId),
    anyDuplicated(expression_genes) == 0,
    anyDuplicated(gene_annotation$EntrezGeneSymbol) == 0,
    setequal(expression_genes, gene_annotation$EntrezGeneSymbol),
    all(expression_genes == gene_annotation$EntrezGeneSymbol),
    if ("AptName" %in% names(gene_annotation)) {
      anyDuplicated(gene_annotation$AptName) == 0
    } else {
      NA
    },
    nrow(expr_gene_df) == EXPECTED_N_SAMPLES,
    length(expression_genes) == EXPECTED_N_GENES
  ),
  observed = c(
    as.character(dplyr::n_distinct(expr_gene_df$SampleId)),
    as.character(dplyr::n_distinct(meta_df$SampleId)),
    as.character(setequal(expr_gene_df$SampleId, meta_df$SampleId)),
    as.character(all(expr_gene_df$SampleId == meta_df$SampleId)),
    as.character(dplyr::n_distinct(expression_genes)),
    as.character(dplyr::n_distinct(gene_annotation$EntrezGeneSymbol)),
    as.character(setequal(expression_genes, gene_annotation$EntrezGeneSymbol)),
    as.character(all(expression_genes == gene_annotation$EntrezGeneSymbol)),
    if ("AptName" %in% names(gene_annotation)) {
      as.character(dplyr::n_distinct(gene_annotation$AptName))
    } else {
      NA_character_
    },
    as.character(nrow(expr_gene_df)),
    as.character(length(expression_genes))
  )
)

safe_write_csv(
  input_alignment_audit,
  file.path(outdir, "qc", "input_alignment_audit.csv")
)

if (any(input_alignment_audit$passed %in% FALSE, na.rm = TRUE)) {
  stop(
    "One or more input-alignment audits failed. ",
    "See qc/input_alignment_audit.csv.",
    call. = FALSE
  )
}

safe_write_csv(
  gene_annotation,
  file.path(outdir, "tables", "gene_level_annotation_used.csv")
)

###############################################################################
# 7) PREPARE datExpr
###############################################################################

rownames_expr <- expr_gene_df$SampleId

expr_only <- expr_gene_df %>%
  dplyr::select(-SampleId) %>%
  as.data.frame(check.names = FALSE)

expr_only[] <- lapply(
  expr_only,
  safe_numeric
)

datExpr0 <- as.data.frame(
  as.matrix(expr_only),
  check.names = FALSE
)

rownames(datExpr0) <- rownames_expr

nonfinite_input <- sum(
  !is.finite(as.matrix(datExpr0))
)

if (REQUIRE_NO_NONFINITE_VALUES && nonfinite_input > 0) {
  stop(
    "Expression input contains ",
    nonfinite_input,
    " missing or non-finite values.",
    call. = FALSE
  )
}

if (nrow(datExpr0) != nrow(meta_df)) {
  stop("datExpr0 and metadata do not contain the same number of samples.", call. = FALSE)
}

if (!all(rownames(datExpr0) == meta_df$SampleId)) {
  stop("datExpr0 and metadata sample order is not identical.", call. = FALSE)
}

cat("Initial datExpr0 dimensions (samples x genes):\n")
cat(dim(datExpr0), "\n\n")

memory_estimate_tbl <- estimate_network_memory(
  ncol(datExpr0)
)

safe_write_csv(
  memory_estimate_tbl,
  file.path(outdir, "qc", "dense_network_memory_estimate.csv")
)

cat("Approximate dense-network memory estimate:\n")
print(memory_estimate_tbl)
cat("\n")

###############################################################################
# 8) QC — goodSamplesGenes
###############################################################################

gsg <- WGCNA::goodSamplesGenes(
  datExpr0,
  verbose = 3
)

datExpr <- datExpr0[
  gsg$goodSamples,
  gsg$goodGenes,
  drop = FALSE
]

meta_clean <- meta_df[
  gsg$goodSamples,
  ,
  drop = FALSE
]

removed_sample_ids <- rownames(datExpr0)[!gsg$goodSamples]
removed_gene_symbols <- colnames(datExpr0)[!gsg$goodGenes]

safe_write_csv(
  tibble::tibble(SampleId = removed_sample_ids),
  file.path(outdir, "qc", "samples_removed_by_goodSamplesGenes.csv")
)

safe_write_csv(
  tibble::tibble(EntrezGeneSymbol = removed_gene_symbols),
  file.path(outdir, "qc", "genes_removed_by_goodSamplesGenes.csv")
)

qc_tbl <- tibble::tibble(
  metric = c(
    "samples_input",
    "genes_input",
    "samples_after_goodSamplesGenes",
    "genes_after_goodSamplesGenes",
    "samples_removed_by_goodSamplesGenes",
    "genes_removed_by_goodSamplesGenes",
    "nonfinite_values_input",
    "wgcna_input_level",
    "somamer_selection_rule"
  ),
  value = c(
    as.character(nrow(datExpr0)),
    as.character(ncol(datExpr0)),
    as.character(nrow(datExpr)),
    as.character(ncol(datExpr)),
    as.character(length(removed_sample_ids)),
    as.character(length(removed_gene_symbols)),
    as.character(nonfinite_input),
    "GENE-COLLAPSED, not aptamer-level",
    "Outcome-independent: lowest missingness > highest log2 MAD > AptName"
  )
)

safe_write_csv(
  qc_tbl,
  file.path(outdir, "qc", "qc_summary.csv")
)

cat("QC summary:\n")
print(qc_tbl)
cat("\n")

if (
  REQUIRE_GSG_TO_PRESERVE_ALL_INPUT &&
  (
    nrow(datExpr) != nrow(datExpr0) ||
    ncol(datExpr) != ncol(datExpr0)
  )
) {
  stop(
    "goodSamplesGenes removed validated samples or genes. ",
    "Review qc/samples_removed_by_goodSamplesGenes.csv and ",
    "qc/genes_removed_by_goodSamplesGenes.csv.",
    call. = FALSE
  )
}

###############################################################################
# 9) SAMPLE CLUSTERING AND OPTIONAL OUTLIER REMOVAL
###############################################################################

sampleTree <- hclust(
  dist(datExpr),
  method = "average"
)

pdf(
  file.path(outdir, "qc", "sample_clustering.pdf"),
  width = 10,
  height = 6
)

plot(
  sampleTree,
  main = "Sample clustering (outcome-independent gene-collapsed matrix)",
  xlab = "",
  sub = ""
)

if (!is.na(OUTLIER_CUT_HEIGHT)) {
  abline(
    h = OUTLIER_CUT_HEIGHT,
    col = "red"
  )
}

dev.off()

datExpr_clean <- datExpr
meta_final <- meta_clean

if (!is.na(OUTLIER_CUT_HEIGHT)) {
  clust <- WGCNA::cutreeStatic(
    sampleTree,
    cutHeight = OUTLIER_CUT_HEIGHT,
    minSize = 10
  )

  keep_idx <- which(clust == 1)

  if (length(keep_idx) >= 10) {
    datExpr_clean <- datExpr[
      keep_idx,
      ,
      drop = FALSE
    ]

    meta_final <- meta_clean[
      keep_idx,
      ,
      drop = FALSE
    ]
  }
}

if (!all(rownames(datExpr_clean) == meta_final$SampleId)) {
  stop(
    "Sample order became misaligned after outlier handling.",
    call. = FALSE
  )
}

outlier_tbl <- tibble::tibble(
  samples_before_outlier_removal = nrow(datExpr),
  samples_after_outlier_removal = nrow(datExpr_clean),
  n_samples_removed = nrow(datExpr) - nrow(datExpr_clean),
  outlier_cut_height = OUTLIER_CUT_HEIGHT,
  manual_outlier_removal_enabled = !is.na(OUTLIER_CUT_HEIGHT)
)

safe_write_csv(
  outlier_tbl,
  file.path(outdir, "qc", "outlier_summary.csv")
)

cat("Outlier summary:\n")
print(outlier_tbl)
cat("\n")

###############################################################################
# 10) SOFT-THRESHOLD SELECTION
###############################################################################

powers <- c(
  1:10,
  seq(12, 30, by = 2)
)

sft <- WGCNA::pickSoftThreshold(
  datExpr_clean,
  powerVector = powers,
  verbose = 5,
  networkType = NETWORK_TYPE
)

sft_scan <- as.data.frame(
  sft$fitIndices
)

if (all(c("SFT.R.sq", "slope") %in% names(sft_scan))) {
  sft_scan <- sft_scan %>%
    dplyr::mutate(
      signed_R2 = -sign(slope) * SFT.R.sq
    )
}

pdf(
  file.path(outdir, "soft_threshold", "soft_threshold_diagnostics.pdf"),
  width = 12,
  height = 6
)

par(mfrow = c(1, 2))

plot(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  xlab = "Soft Threshold (power)",
  ylab = "Scale Free Topology Model Fit, signed R²",
  type = "n",
  main = "Scale independence"
)

text(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  labels = sft$fitIndices[, 1],
  cex = 0.8,
  col = "darkred"
)

abline(
  h = TARGET_SCALE_FREE_R2,
  col = "blue",
  lty = 2
)

plot(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  xlab = "Soft Threshold (power)",
  ylab = "Mean Connectivity",
  type = "n",
  main = "Mean connectivity"
)

text(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  labels = sft$fitIndices[, 1],
  cex = 0.8,
  col = "darkred"
)

dev.off()

softPower <- select_soft_power(
  sft_fit_indices = sft$fitIndices,
  default_power = DEFAULT_SOFT_POWER,
  target_r2 = TARGET_SCALE_FREE_R2,
  auto_select = AUTO_SELECT_SOFT_POWER
)

selected_power_row <- sft_scan$Power == softPower

selected_signed_R2 <- safe_pick_value(
  sft_scan,
  selected_power_row,
  c("signed_R2")
)

selected_scale_free_R2_raw <- safe_pick_value(
  sft_scan,
  selected_power_row,
  c("SFT.R.sq")
)

selected_mean_connectivity <- safe_pick_value(
  sft_scan,
  selected_power_row,
  c("mean.k.", "mean.k", "mean_connectivity")
)

selected_slope <- safe_pick_value(
  sft_scan,
  selected_power_row,
  c("slope")
)

soft_tbl <- tibble::tibble(
  chosen_soft_power = softPower,
  default_soft_power_if_needed = DEFAULT_SOFT_POWER,
  auto_select_soft_power = AUTO_SELECT_SOFT_POWER,
  target_scale_free_R2 = TARGET_SCALE_FREE_R2,
  network_type = NETWORK_TYPE,
  signed_scale_free_R2_at_chosen_power = selected_signed_R2,
  raw_SFT_R_sq_at_chosen_power = selected_scale_free_R2_raw,
  slope_at_chosen_power = selected_slope,
  mean_connectivity_at_chosen_power = selected_mean_connectivity,
  target_R2_reached_at_chosen_power = (
    is.finite(selected_signed_R2) &&
      selected_signed_R2 >= TARGET_SCALE_FREE_R2
  )
)

safe_write_csv(
  sft_scan,
  file.path(outdir, "soft_threshold", "soft_threshold_scan.csv")
)

safe_write_csv(
  soft_tbl,
  file.path(outdir, "soft_threshold", "soft_threshold_selected.csv")
)

cat("Selected soft power: ", softPower, "\n")
cat("Signed scale-free R² at selected power: ", selected_signed_R2, "\n")
cat("Mean connectivity at selected power: ", selected_mean_connectivity, "\n\n")

if (
  is.finite(selected_signed_R2) &&
  selected_signed_R2 < TARGET_SCALE_FREE_R2
) {
  warning(
    "The retained prespecified beta = ",
    softPower,
    " did not reach the target signed scale-free R² of ",
    TARGET_SCALE_FREE_R2,
    ". Review soft-threshold diagnostics before manuscript reporting."
  )
}

###############################################################################
# 11) NETWORK CONSTRUCTION
###############################################################################

gc(verbose = FALSE)

adjacency_mat <- WGCNA::adjacency(
  datExpr_clean,
  power = softPower,
  type = NETWORK_TYPE
)

dimnames(adjacency_mat) <- list(
  colnames(datExpr_clean),
  colnames(datExpr_clean)
)

if (
  nrow(adjacency_mat) != ncol(datExpr_clean) ||
  ncol(adjacency_mat) != ncol(datExpr_clean)
) {
  stop("Adjacency dimensions do not match the expression feature count.", call. = FALSE)
}

if (!isTRUE(all.equal(adjacency_mat, t(adjacency_mat), tolerance = 1e-10))) {
  stop("Adjacency matrix is not symmetric.", call. = FALSE)
}

gc(verbose = FALSE)

TOM <- WGCNA::TOMsimilarity(
  adjacency_mat,
  TOMType = NETWORK_TYPE
)

dimnames(TOM) <- list(
  colnames(datExpr_clean),
  colnames(datExpr_clean)
)

if (
  nrow(TOM) != ncol(datExpr_clean) ||
  ncol(TOM) != ncol(datExpr_clean)
) {
  stop("TOM dimensions do not match the expression feature count.", call. = FALSE)
}

if (!isTRUE(all.equal(TOM, t(TOM), tolerance = 1e-10))) {
  stop("TOM matrix is not symmetric.", call. = FALSE)
}

dissTOM <- 1 - TOM

dimnames(dissTOM) <- list(
  colnames(datExpr_clean),
  colnames(datExpr_clean)
)

network_matrix_audit <- tibble::tibble(
  object = c(
    "adjacency_mat",
    "TOM",
    "dissTOM"
  ),
  n_rows = c(
    nrow(adjacency_mat),
    nrow(TOM),
    nrow(dissTOM)
  ),
  n_cols = c(
    ncol(adjacency_mat),
    ncol(TOM),
    ncol(dissTOM)
  ),
  has_row_names = c(
    !is.null(rownames(adjacency_mat)),
    !is.null(rownames(TOM)),
    !is.null(rownames(dissTOM))
  ),
  has_col_names = c(
    !is.null(colnames(adjacency_mat)),
    !is.null(colnames(TOM)),
    !is.null(colnames(dissTOM))
  ),
  feature_order_matches_datExpr = c(
    identical(rownames(adjacency_mat), colnames(datExpr_clean)),
    identical(rownames(TOM), colnames(datExpr_clean)),
    identical(rownames(dissTOM), colnames(datExpr_clean))
  ),
  symmetric = c(
    isTRUE(all.equal(adjacency_mat, t(adjacency_mat), tolerance = 1e-10)),
    isTRUE(all.equal(TOM, t(TOM), tolerance = 1e-10)),
    isTRUE(all.equal(dissTOM, t(dissTOM), tolerance = 1e-10))
  )
)

safe_write_csv(
  network_matrix_audit,
  file.path(outdir, "qc", "network_matrix_alignment_audit.csv")
)

if (
  any(!network_matrix_audit$feature_order_matches_datExpr) ||
  any(!network_matrix_audit$symmetric)
) {
  stop(
    "Network-matrix alignment audit failed. ",
    "See qc/network_matrix_alignment_audit.csv.",
    call. = FALSE
  )
}

if (SAVE_NETWORK_CHECKPOINT_AFTER_TOM) {
  save(
    datExpr_clean,
    meta_final,
    gene_annotation,
    softPower,
    adjacency_mat,
    TOM,
    dissTOM,
    file = file.path(
      outdir,
      "workspace",
      "wgcna_network_checkpoint_after_TOM.RData"
    )
  )
}

geneTree <- hclust(
  as.dist(dissTOM),
  method = "average"
)

pdf(
  file.path(outdir, "network", "gene_clustering_dendrogram.pdf"),
  width = 12,
  height = 6
)

plot(
  geneTree,
  xlab = "",
  sub = "",
  main = "Gene clustering on TOM-based dissimilarity",
  labels = FALSE,
  hang = 0.04
)

dev.off()

###############################################################################
# 12) INITIAL MODULE DETECTION
###############################################################################

dynamicMods <- dynamicTreeCut::cutreeDynamic(
  dendro = geneTree,
  distM = dissTOM,
  deepSplit = DEEPSPLIT,
  pamRespectsDendro = FALSE,
  minClusterSize = MIN_MODULE_SIZE
)

names(dynamicMods) <- colnames(datExpr_clean)

dynamicColors <- WGCNA::labels2colors(
  dynamicMods
)

names(dynamicColors) <- colnames(datExpr_clean)

pdf(
  file.path(outdir, "modules", "dynamic_modules_dendrogram.pdf"),
  width = 12,
  height = 6
)

WGCNA::plotDendroAndColors(
  geneTree,
  dynamicColors,
  "Dynamic Tree Cut",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Gene dendrogram and initial module colors"
)

dev.off()

###############################################################################
# 13) MODULE EIGENGENES AND MERGING
###############################################################################

MEList <- WGCNA::moduleEigengenes(
  datExpr_clean,
  colors = dynamicColors
)

MEs <- WGCNA::orderMEs(
  MEList$eigengenes
)

rownames(MEs) <- rownames(datExpr_clean)

MEDiss <- 1 - stats::cor(
  MEs,
  use = "pairwise.complete.obs"
)

METree <- hclust(
  as.dist(MEDiss),
  method = "average"
)

pdf(
  file.path(outdir, "modules", "module_eigengene_clustering.pdf"),
  width = 7,
  height = 6
)

plot(
  METree,
  main = "Clustering of module eigengenes"
)

abline(
  h = MERGE_CUT_HEIGHT,
  col = "red"
)

dev.off()

merge <- WGCNA::mergeCloseModules(
  datExpr_clean,
  dynamicColors,
  cutHeight = MERGE_CUT_HEIGHT,
  verbose = 3
)

mergedColors <- as.character(
  merge$colors
)

names(mergedColors) <- colnames(datExpr_clean)

mergedMEs <- WGCNA::orderMEs(
  merge$newMEs
)

rownames(mergedMEs) <- rownames(datExpr_clean)

if (length(mergedColors) != ncol(datExpr_clean)) {
  stop("Merged module labels do not match the expression feature count.", call. = FALSE)
}

if (!identical(names(mergedColors), colnames(datExpr_clean))) {
  stop("Merged module labels are not aligned to expression genes.", call. = FALSE)
}

if (nrow(mergedMEs) != nrow(datExpr_clean)) {
  stop("Merged eigengenes do not match the expression sample count.", call. = FALSE)
}

if (!identical(rownames(mergedMEs), rownames(datExpr_clean))) {
  stop("Merged eigengenes are not aligned to expression samples.", call. = FALSE)
}

pdf(
  file.path(outdir, "modules", "merged_modules_dendrogram.pdf"),
  width = 12,
  height = 6
)

WGCNA::plotDendroAndColors(
  geneTree,
  cbind(
    dynamicColors,
    mergedColors
  ),
  c(
    "Initial modules",
    "Merged modules"
  ),
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Gene dendrogram and merged module colors"
)

dev.off()

###############################################################################
# 14) MODULE COUNTS
###############################################################################

module_counts <- sort(
  table(mergedColors),
  decreasing = TRUE
)

module_count_tbl <- tibble::tibble(
  Module = names(module_counts),
  N_genes = as.integer(module_counts),
  Prop_total = as.integer(module_counts) / sum(module_counts)
)

safe_write_csv(
  module_count_tbl,
  file.path(outdir, "tables", "module_counts.csv")
)

pdf(
  file.path(outdir, "modules", "module_barplot.pdf"),
  width = 9,
  height = 6
)

par(mar = c(8, 5, 3, 1))

barplot(
  module_counts,
  las = 2,
  col = names(module_counts),
  main = "Number of genes per module",
  ylab = "Number of genes",
  border = NA
)

dev.off()

cat("Module count table:\n")
print(module_count_tbl)
cat("\n")

###############################################################################
# 15) GENE → MODULE ASSIGNMENT TABLE
###############################################################################

if (!"Protein_Name" %in% names(gene_annotation)) {
  gene_annotation$Protein_Name <- gene_annotation$EntrezGeneSymbol
}

if (!"Protein_Display" %in% names(gene_annotation)) {
  gene_annotation$Protein_Display <- gene_annotation$Protein_Name
}

gene_annotation <- gene_annotation %>%
  dplyr::mutate(
    Protein_Name = dplyr::coalesce(
      clean_text_na(Protein_Name),
      EntrezGeneSymbol
    ),
    Protein_Display = dplyr::coalesce(
      clean_text_na(Protein_Display),
      Protein_Name,
      EntrezGeneSymbol
    )
  )

gene_module_assignment <- tibble::tibble(
  EntrezGeneSymbol = colnames(datExpr_clean),
  Module = unname(mergedColors)
) %>%
  dplyr::left_join(
    gene_annotation,
    by = "EntrezGeneSymbol"
  ) %>%
  dplyr::mutate(
    Protein_Name = dplyr::coalesce(
      Protein_Name,
      EntrezGeneSymbol
    ),
    Protein_Display = dplyr::coalesce(
      Protein_Display,
      Protein_Name,
      EntrezGeneSymbol
    )
  ) %>%
  dplyr::left_join(
    module_count_tbl,
    by = "Module"
  ) %>%
  dplyr::arrange(
    Module,
    Protein_Name
  )

if (nrow(gene_module_assignment) != ncol(datExpr_clean)) {
  stop(
    "Gene-module assignment does not contain exactly one row per WGCNA gene.",
    call. = FALSE
  )
}

if (anyDuplicated(gene_module_assignment$EntrezGeneSymbol) > 0) {
  stop(
    "Gene-module assignment contains duplicated genes.",
    call. = FALSE
  )
}

if (any(is.na(gene_module_assignment$Module))) {
  stop(
    "Gene-module assignment contains missing module labels.",
    call. = FALSE
  )
}

safe_write_csv(
  gene_module_assignment,
  file.path(outdir, "tables", "gene_module_assignment.csv")
)

module_assignment_audit <- tibble::tibble(
  metric = c(
    "n_genes_in_datExpr_clean",
    "n_gene_module_assignment_rows",
    "n_unique_assigned_genes",
    "n_missing_module_labels",
    "n_modules_including_grey",
    "n_modules_excluding_grey",
    "n_grey_genes",
    "smallest_non_grey_module",
    "largest_non_grey_module",
    "sum_module_counts"
  ),
  value = c(
    ncol(datExpr_clean),
    nrow(gene_module_assignment),
    dplyr::n_distinct(gene_module_assignment$EntrezGeneSymbol),
    sum(is.na(gene_module_assignment$Module)),
    dplyr::n_distinct(gene_module_assignment$Module),
    dplyr::n_distinct(
      gene_module_assignment$Module[
        gene_module_assignment$Module != "grey"
      ]
    ),
    sum(gene_module_assignment$Module == "grey"),
    suppressWarnings(min(
      module_count_tbl$N_genes[
        module_count_tbl$Module != "grey"
      ],
      na.rm = TRUE
    )),
    suppressWarnings(max(
      module_count_tbl$N_genes[
        module_count_tbl$Module != "grey"
      ],
      na.rm = TRUE
    )),
    sum(module_count_tbl$N_genes)
  )
)

safe_write_csv(
  module_assignment_audit,
  file.path(outdir, "qc", "module_assignment_audit.csv")
)

###############################################################################
# 16) MODULE EIGENGENES PER SAMPLE
###############################################################################

ME_df <- dplyr::bind_cols(
  meta_final %>%
    dplyr::select(SampleId),
  tibble::as_tibble(
    mergedMEs,
    .name_repair = "minimal"
  )
)

if (nrow(ME_df) != nrow(meta_final)) {
  stop(
    "Eigengene table does not match final metadata sample count.",
    call. = FALSE
  )
}

if (!all(ME_df$SampleId == meta_final$SampleId)) {
  stop(
    "Eigengene table is not aligned to final metadata sample order.",
    call. = FALSE
  )
}

safe_write_csv(
  ME_df,
  file.path(outdir, "eigengenes", "module_eigengenes_per_sample.csv")
)

eigengene_alignment_audit <- tibble::tibble(
  metric = c(
    "n_samples_datExpr_clean",
    "n_samples_meta_final",
    "n_rows_eigengene_table",
    "n_eigengene_columns",
    "sample_order_datExpr_vs_meta",
    "sample_order_meta_vs_eigengene_table",
    "all_eigengene_values_finite"
  ),
  value = c(
    nrow(datExpr_clean),
    nrow(meta_final),
    nrow(ME_df),
    ncol(ME_df) - 1,
    all(rownames(datExpr_clean) == meta_final$SampleId),
    all(meta_final$SampleId == ME_df$SampleId),
    all(is.finite(as.matrix(ME_df[, -1, drop = FALSE])))
  )
)

safe_write_csv(
  eigengene_alignment_audit,
  file.path(outdir, "qc", "eigengene_alignment_audit.csv")
)

cat("Eigengene table dimensions:\n")
print(dim(ME_df))
cat("\n")

###############################################################################
# 17) FINAL SUMMARY TABLE
###############################################################################

summary_tbl <- tibble::tibble(
  metric = c(
    "base_dir",
    "input_dir",
    "output_dir",
    "wgcna_input_level",
    "somamer_selection",
    "n_samples_input_script11",
    "n_genes_input_script11",
    "n_samples_final_wgcna",
    "n_genes_final_wgcna",
    "soft_power",
    "signed_scale_free_R2_at_selected_power",
    "mean_connectivity_at_selected_power",
    "target_scale_free_R2",
    "target_scale_free_R2_reached",
    "auto_select_soft_power",
    "network_type",
    "min_module_size",
    "deepSplit",
    "merge_cut_height",
    "manual_outlier_cut_height",
    "n_modules_final_including_grey",
    "n_modules_final_excluding_grey",
    "n_grey_genes",
    "full_workspace_contains_adjacency_TOM_dissTOM"
  ),
  value = c(
    base_dir,
    input_dir,
    outdir,
    "GENE-COLLAPSED, not aptamer-level",
    "Outcome-independent: lowest missingness > highest log2 MAD > AptName",
    as.character(nrow(datExpr0)),
    as.character(ncol(datExpr0)),
    as.character(nrow(datExpr_clean)),
    as.character(ncol(datExpr_clean)),
    as.character(softPower),
    as.character(selected_signed_R2),
    as.character(selected_mean_connectivity),
    as.character(TARGET_SCALE_FREE_R2),
    as.character(
      is.finite(selected_signed_R2) &&
        selected_signed_R2 >= TARGET_SCALE_FREE_R2
    ),
    as.character(AUTO_SELECT_SOFT_POWER),
    NETWORK_TYPE,
    as.character(MIN_MODULE_SIZE),
    as.character(DEEPSPLIT),
    as.character(MERGE_CUT_HEIGHT),
    as.character(OUTLIER_CUT_HEIGHT),
    as.character(dplyr::n_distinct(mergedColors)),
    as.character(dplyr::n_distinct(mergedColors[mergedColors != "grey"])),
    as.character(sum(mergedColors == "grey")),
    "TRUE"
  )
)

safe_write_csv(
  summary_tbl,
  file.path(outdir, "tables", "wgcna_core_summary.csv")
)

cat("Final summary:\n")
print(summary_tbl)
cat("\n")

###############################################################################
# 18) OUTPUT MANIFEST
###############################################################################

output_manifest <- tibble::tibble(
  output_file = c(
    "qc/input_alignment_audit.csv",
    "qc/qc_summary.csv",
    "qc/sample_clustering.pdf",
    "qc/outlier_summary.csv",
    "qc/dense_network_memory_estimate.csv",
    "qc/network_matrix_alignment_audit.csv",
    "qc/module_assignment_audit.csv",
    "qc/eigengene_alignment_audit.csv",
    "soft_threshold/soft_threshold_scan.csv",
    "soft_threshold/soft_threshold_selected.csv",
    "soft_threshold/soft_threshold_diagnostics.pdf",
    "network/gene_clustering_dendrogram.pdf",
    "modules/dynamic_modules_dendrogram.pdf",
    "modules/module_eigengene_clustering.pdf",
    "modules/merged_modules_dendrogram.pdf",
    "modules/module_barplot.pdf",
    "tables/gene_level_annotation_used.csv",
    "tables/gene_module_assignment.csv",
    "tables/module_counts.csv",
    "tables/wgcna_core_summary.csv",
    "eigengenes/module_eigengenes_per_sample.csv",
    "workspace/wgcna_core_collapsed_workspace.RData",
    "workspace/wgcna_core_light_workspace.RData",
    "sessionInfo.txt"
  ),
  description = c(
    "Strict sample and feature alignment audit.",
    "goodSamplesGenes and input-level QC summary.",
    "Hierarchical sample-clustering diagnostic.",
    "Manual outlier-removal summary.",
    "Approximate dense-matrix memory requirement.",
    "Adjacency, TOM and dissTOM dimension/name/symmetry audit.",
    "Gene-to-module assignment audit.",
    "Module-eigengene sample-alignment audit.",
    "Full soft-threshold scan, including signed R2 when available.",
    "Selected power and diagnostics at beta.",
    "Scale-free topology fit and mean-connectivity figure.",
    "TOM-based gene dendrogram.",
    "Initial dynamic-tree-cut modules.",
    "Module eigengene clustering before merging.",
    "Initial and merged module colors.",
    "Final module-size bar plot.",
    "Outcome-independent gene annotation used for network construction.",
    "Final gene-to-module assignment with SOMAmer and DEP annotations.",
    "Final module sizes.",
    "Core WGCNA parameters and results summary.",
    "Module eigengenes aligned to SampleId.",
    "Full workspace preserving the original downstream object contract.",
    "Reduced workspace for biological interpretation and module-trait scripts.",
    "R session information."
  )
)

safe_write_csv(
  output_manifest,
  file.path(outdir, "wgcna_core_output_manifest.csv")
)

###############################################################################
# 19) SAVE WORKSPACES
###############################################################################

# Preserve the full object contract of the original Script 02. This workspace
# intentionally includes adjacency_mat, TOM and dissTOM because later quality
# diagnostics and network figures use them.
save(
  expr_gene_df,
  meta_df,
  meta_clean,
  meta_final,
  gene_annotation,
  manifest_df,
  input_qc_df,
  gene_map_audit_df,
  expression_scale_audit_df,
  datExpr0,
  datExpr,
  datExpr_clean,
  gsg,
  sampleTree,
  sft,
  sft_scan,
  softPower,
  soft_tbl,
  adjacency_mat,
  TOM,
  dissTOM,
  geneTree,
  dynamicMods,
  dynamicColors,
  MEList,
  MEs,
  MEDiss,
  METree,
  merge,
  mergedColors,
  mergedMEs,
  module_count_tbl,
  gene_module_assignment,
  ME_df,
  qc_tbl,
  outlier_tbl,
  memory_estimate_tbl,
  network_matrix_audit,
  module_assignment_audit,
  eigengene_alignment_audit,
  summary_tbl,
  output_manifest,
  file = file.path(
    outdir,
    "workspace",
    "wgcna_core_collapsed_workspace.RData"
  )
)

# Lightweight workspace for scripts that do not need dense network matrices.
save(
  meta_final,
  gene_annotation,
  datExpr_clean,
  softPower,
  mergedColors,
  mergedMEs,
  module_count_tbl,
  gene_module_assignment,
  ME_df,
  summary_tbl,
  file = file.path(
    outdir,
    "workspace",
    "wgcna_core_light_workspace.RData"
  )
)

writeLines(
  capture.output(utils::sessionInfo()),
  con = file.path(outdir, "sessionInfo.txt")
)

###############################################################################
# 20) FINAL MESSAGE
###############################################################################

cat("\nWGCNA core workflow on OUTCOME-INDEPENDENT GENE-COLLAPSED features finished successfully.\n")
cat("Main output directory:\n", outdir, "\n")
cat("Input level: GENE-COLLAPSED, not aptamer-level.\n")
cat("Samples used: ", nrow(datExpr_clean), "\n", sep = "")
cat("Genes used: ", ncol(datExpr_clean), "\n", sep = "")
cat("Selected beta: ", softPower, "\n", sep = "")
cat("Signed scale-free R2 at beta: ", selected_signed_R2, "\n", sep = "")
cat("Mean connectivity at beta: ", selected_mean_connectivity, "\n", sep = "")
cat("Modules including grey: ", dplyr::n_distinct(mergedColors), "\n", sep = "")
cat("Modules excluding grey: ", dplyr::n_distinct(mergedColors[mergedColors != "grey"]), "\n", sep = "")
cat("Grey genes: ", sum(mergedColors == "grey"), "\n", sep = "")
cat("Main outputs for Script 12:\n")
cat(file.path(outdir, "workspace", "wgcna_core_collapsed_workspace.RData"), "\n")
cat(file.path(outdir, "tables", "gene_module_assignment.csv"), "\n")
cat(file.path(outdir, "eigengenes", "module_eigengenes_per_sample.csv"), "\n")

###############################################################################
# END
###############################################################################

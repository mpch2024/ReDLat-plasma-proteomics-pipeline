###############################################################################
# ReDLat plasma proteomics — machine-learning workflow
# 11. Fit country-held-out DEP models
# Requires: Script 10 outputs, reviewer metadata and proteomics
# Produces: country-specific training candidate lists
# Data policy: participant-level inputs and predictions remain local.
###############################################################################

.project_root <- if (nzchar(Sys.getenv("REDLAT_PROJECT_ROOT", unset = ""))) {
  normalizePath(Sys.getenv("REDLAT_PROJECT_ROOT"), winslash = "/", mustWork = TRUE)
} else if (requireNamespace("here", quietly = TRUE)) {
  normalizePath(here::here(), winslash = "/", mustWork = TRUE)
} else {
  stop("Package 'here' is required. Restore the R environment before running the workflow.", call. = FALSE)
}
source(file.path(.project_root, "R", "ml_bootstrap.R"), local = FALSE)
source(file.path(.project_root, "R", "reviewer_data_adapter.R"), local = FALSE)
ML_CONFIG <- ml_load_config(.project_root)
EXCLUDED_SAMPLE_IDS <- ml_read_excluded_ids(ML_CONFIG)

project_root <- ML_CONFIG$project_root
FOLDS_DIR <- file.path(ML_CONFIG$private_root, "loco", "00_folds")
data_dir <- ML_CONFIG$data_dir
outdir <- ML_CONFIG$private_root
csv_file <- ML_CONFIG$metadata_file
adat_file <- ML_CONFIG$adat_file
OUT_DIR <- file.path(ML_CONFIG$private_root, "loco", "01_dep_folds")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
ml_assert_private_path(OUT_DIR, ML_CONFIG)

MAIN_GROUPS <- c("CN", "AD")
MAIN_FDR <- 0.05
STRICT_FDR <- 0.05
LOGFC_THRESHOLD <- 0
SEED_GLOBAL <- 1234

set.seed(SEED_GLOBAL)
options(stringsAsFactors = FALSE)

required_pkgs <- c(
  "dplyr", "tidyr", "purrr", "tibble", "stringr",
  "readr", "limma", "openxlsx"
)

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0L) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
       ". Restore the R environment before running this script.", call. = FALSE)
}
invisible(lapply(required_pkgs, library, character.only = TRUE))

missing_input_files <- c(csv_file, ML_CONFIG$proteomics_file, ML_CONFIG$annotation_file)[!file.exists(c(csv_file, ML_CONFIG$proteomics_file, ML_CONFIG$annotation_file))]
if (length(missing_input_files) > 0) {
  stop("Missing required input files:\n", paste(missing_input_files, collapse = "\n"))
}

###############################################################################
# 01_helpers
###############################################################################

safe_write_csv <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(x, file)
}

required_columns <- function(df, cols, object_name = "data frame") {
  missing_cols <- setdiff(cols, names(df))
  if (length(missing_cols) > 0) {
    stop(object_name, " is missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  invisible(TRUE)
}

clean_text_na <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA", "NaN", "NULL", "null", "N/A")] <- NA_character_
  x
}

norm_char <- function(x) trimws(as.character(x))

pick_col <- function(df, candidates, default = NA_real_) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) > 0) df[[hit[1]]] else rep(default, nrow(df))
}

safe_log2_matrix <- function(mat) {
  mat <- apply(mat, 2, as.numeric)
  mat[mat <= 0] <- NA_real_
  log2(mat)
}

safe_se_from_limma <- function(logFC, t_stat) {
  logFC <- suppressWarnings(as.numeric(logFC))
  t_stat <- suppressWarnings(as.numeric(t_stat))
  se <- rep(NA_real_, length(logFC))
  ok <- is.finite(logFC) & is.finite(t_stat) & abs(t_stat) > .Machine$double.eps
  se[ok] <- abs(logFC[ok] / t_stat[ok])
  se
}

make_protein_name <- function(entrez, target_full, target, apt, raw) {
  entrez <- clean_text_na(entrez)
  target_full <- clean_text_na(target_full)
  target <- clean_text_na(target)
  apt <- clean_text_na(apt)
  raw <- clean_text_na(raw)
  dplyr::coalesce(entrez, target_full, target, apt, raw)
}

validate_design_matrix <- function(design, model_name = "model") {
  rank_design <- qr(design)$rank
  n_col <- ncol(design)
  out <- tibble::tibble(
    model = model_name,
    n_samples = nrow(design),
    n_parameters = n_col,
    design_rank = rank_design,
    full_rank = rank_design == n_col,
    columns = paste(colnames(design), collapse = " | "),
    alias_diagnostic = NA_character_
  )
  if (rank_design < n_col) {
    out$alias_diagnostic <- "Rank deficient design; inspect sparse factors/covariates."
    warning(model_name, " design matrix is rank deficient. See model_design_diagnostics.csv.")
  }
  out
}

build_annotation_table <- function(soma_info_internal, seq_cols) {
  seq_key_tbl <- tibble::tibble(AptName = as.character(seq_cols))
  internal_tbl <- soma_info_internal %>%
    dplyr::mutate(
      AptName = as.character(AptName),
      SeqId = as.character(SeqId),
      SomaId = as.character(SomaId),
      TargetFullName = as.character(TargetFullName),
      Target = as.character(Target),
      EntrezGeneID = suppressWarnings(as.numeric(EntrezGeneID)),
      EntrezGeneSymbol = as.character(EntrezGeneSymbol),
      UniProt = as.character(UniProt),
      Organism = as.character(Organism),
      Type = as.character(Type)
    ) %>%
    dplyr::select(AptName, SeqId, SomaId, TargetFullName, Target, EntrezGeneID, EntrezGeneSymbol, UniProt, Organism, Type) %>%
    dplyr::distinct(AptName, .keep_all = TRUE)

  seq_key_tbl %>%
    dplyr::left_join(internal_tbl, by = "AptName") %>%
    dplyr::mutate(Protein_Name = make_protein_name(EntrezGeneSymbol, TargetFullName, Target, AptName, SeqId))
}

classify_dep <- function(dep_tbl, fdr = MAIN_FDR, logfc = LOGFC_THRESHOLD) {
  dep_tbl %>%
    dplyr::mutate(
      type = dplyr::case_when(
        logFC >  logfc & adj.P.Val < fdr ~ "Up",
        logFC < -logfc & adj.P.Val < fdr ~ "Down",
        TRUE ~ "NS"
      ),
      Direction = dplyr::case_when(
        type == "Up" ~ "Higher in AD",
        type == "Down" ~ "Lower in AD",
        TRUE ~ "Not significant"
      )
    )
}

collapse_dep_to_gene <- function(dep_tbl) {
  dep_tbl %>%
    dplyr::filter(!is.na(EntrezGeneSymbol), EntrezGeneSymbol != "") %>%
    dplyr::arrange(adj.P.Val, dplyr::desc(abs(logFC)), AptName) %>%
    dplyr::distinct(EntrezGeneSymbol, .keep_all = TRUE)
}

export_dep_table <- function(dep_tbl, file) {
  keep_cols <- intersect(
    c(
      "Protein_Name", "EntrezGeneSymbol", "TargetFullName", "Target", "UniProt",
      "AptName", "SeqId", "feature_id_raw", "logFC", "se", "AveExpr",
      "t", "P.Value", "adj.P.Val", "B", "type", "Direction"
    ),
    names(dep_tbl)
  )
  dep_tbl %>%
    dplyr::select(dplyr::all_of(keep_cols)) %>%
    dplyr::arrange(adj.P.Val) %>%
    safe_write_csv(file)
}

summarize_dep_counts <- function(dep_tbl, fdr_values = c(STRICT_FDR, MAIN_FDR), universe_label = "gene_collapsed") {
  purrr::map_dfr(fdr_values, function(fdr_value) {
    tibble::tibble(
      universe = universe_label,
      fdr = fdr_value,
      sig_total = sum(dep_tbl$adj.P.Val < fdr_value, na.rm = TRUE),
      up = sum(dep_tbl$adj.P.Val < fdr_value & dep_tbl$logFC > 0, na.rm = TRUE),
      down = sum(dep_tbl$adj.P.Val < fdr_value & dep_tbl$logFC < 0, na.rm = TRUE)
    )
  })
}

run_limma_dep_model <- function(dat, seq_cols, annot_tbl, formula_str, coef_name = "SampleGroupAD", fdr = MAIN_FDR, model_name = "main_DEP") {
  dat <- dat %>% tibble::as_tibble()
  protein_cols <- intersect(seq_cols, names(dat))
  if (length(protein_cols) == 0) stop("No protein columns found for limma model.")

  metadata <- dat %>% dplyr::select(-dplyr::all_of(protein_cols))
  model_vars <- all.vars(stats::as.formula(formula_str))
  missing_model_vars <- setdiff(model_vars, names(metadata))
  if (length(missing_model_vars) > 0) {
    stop("Missing model variables in ", model_name, ": ", paste(missing_model_vars, collapse = ", "))
  }

  keep <- stats::complete.cases(metadata[, model_vars, drop = FALSE])
  dat <- dat[keep, , drop = FALSE]
  metadata <- metadata[keep, , drop = FALSE]

  expr <- dat %>% dplyr::select(dplyr::all_of(protein_cols)) %>% as.matrix()
  expr <- t(safe_log2_matrix(expr))

  design <- model.matrix(stats::as.formula(formula_str), data = metadata)
  design_diag <- validate_design_matrix(design, model_name = model_name)
  if (!coef_name %in% colnames(design)) {
    stop("Coefficient not found in design matrix for ", model_name, ": ", coef_name)
  }

  fit <- limma::lmFit(expr, design)
  fit <- limma::eBayes(fit)

  tt <- limma::topTable(fit, coef = coef_name, adjust.method = "BH", number = Inf) %>%
    tibble::rownames_to_column(var = "feature_id_raw") %>%
    dplyr::mutate(AptName = as.character(feature_id_raw), se = safe_se_from_limma(logFC, t)) %>%
    dplyr::left_join(annot_tbl, by = "AptName") %>%
    dplyr::mutate(Protein_Name = make_protein_name(EntrezGeneSymbol, TargetFullName, Target, AptName, feature_id_raw)) %>%
    classify_dep(fdr = fdr)

  list(dep = tt, fit = fit, design = design, metadata = metadata, protein_cols = protein_cols, design_diagnostic = design_diag)
}

###############################################################################
# 02_import_metadata_adat_annotation
###############################################################################

message("Importing reviewer metadata and proteomics...")

reviewer_inputs <- reviewer_load_aptamer_data(
  metadata_file = csv_file,
  proteomics_file = ML_CONFIG$proteomics_file,
  annotation_file = ML_CONFIG$annotation_file
)

meta_info_new <- reviewer_inputs$metadata
required_columns(meta_info_new, c("SampleId", "Sex", "Age", "Country", "Education", "ApoE"), "reviewer metadata")

my_adat <- reviewer_inputs$proteomics_raw_compat
soma_info_internal <- reviewer_inputs$annotation

sample_data <- meta_info_new %>%
  dplyr::left_join(my_adat, by = "SampleId")

if (!"SampleType" %in% names(sample_data)) sample_data$SampleType <- "Sample"
if (!"SampleGroup" %in% names(sample_data)) stop("SampleGroup was not found in reviewer metadata.")
if (!"RowCheck" %in% names(sample_data)) sample_data$RowCheck <- NA_character_

# Reviewer release is already post-QC; do not introduce a new RowCheck exclusion.
seq_cols <- grep("^seq[._]", names(sample_data), value = TRUE)
if (length(seq_cols) == 0) stop("No seq.* or seq_* columns found after reviewer-data merge.")

annot_tbl <- reviewer_inputs$annotation
protein_universe <- reviewer_inputs$protein_universe
if (length(protein_universe) == 0) stop("No proteins matched the reviewer annotation universe.")

###############################################################################
# 03_nested_DEP_for_ML_candidate_space
###############################################################################

message("Running nested DEP models for ML candidate space...")

# FOLDS_DIR is defined from ML_CONFIG above; legacy external override removed.
main_formula <- "~ SampleGroup + Age + Sex + Country + Education"
coef_name <- "SampleGroupAD"

dep_df <- sample_data %>%
  dplyr::filter(!is.na(Age),
                !is.na(Sex),
                !is.na(Country),
                !is.na(Education)) %>%
  dplyr::filter(SampleType == "Sample") %>%
  tidyr::drop_na(SampleGroup) %>%
  dplyr::filter(SampleGroup %in% MAIN_GROUPS) %>%
  dplyr::mutate(
    SampleGroup = factor(SampleGroup, levels = MAIN_GROUPS),
    Sex = factor(Sex),
    Country = factor(Country)
  ) %>%
  dplyr::filter(!SampleId %in% EXCLUDED_SAMPLE_IDS) %>%
  tibble::as_tibble()

dep_protein_cols <- intersect(
  protein_universe,
  names(dep_df)
)

if (length(dep_protein_cols) == 0) {
  stop("No protein columns found for DEP.")
}

nested_summary <- list()

countries <- c(
  "Argentina",
  "Chile",
  "Colombia",
  "Mexico",
  "Peru"
)

for (country in countries) {

  message(
    "\n====================================\n",
    "Processing Fold ", country,
    "\n===================================="
  )

  train_ids <- readr::read_csv(
    file.path(
      FOLDS_DIR,
      paste0(country, "_train_ids.csv")
    ),
    show_col_types = FALSE
  )

  dep_df_fold <- dep_df %>%
    dplyr::filter(
      SampleId %in% train_ids$SampleId
    )
  dep_df_fold$Country <- droplevels(dep_df_fold$Country)
  print(table(dep_df_fold$Country))
  print(levels(dep_df_fold$Country))
  fold_dir <- file.path(
    OUT_DIR,
    country
  )

  dir.create(
    fold_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  main_fit <- run_limma_dep_model(
    dat = dep_df_fold,
    seq_cols = dep_protein_cols,
    annot_tbl = annot_tbl,
    formula_str = main_formula,
    coef_name = coef_name,
    fdr = MAIN_FDR,
    model_name = paste0(
      "LOCO_",
      country
    )
  )

  DEP_aptamer <- main_fit$dep %>%
    dplyr::filter(
      AptName %in% dep_protein_cols
    )

  DEP_gene <- collapse_dep_to_gene(
    DEP_aptamer
  )

  DEP_gene_STRICT  <- DEP_gene %>%
    dplyr::filter(
      adj.P.Val < STRICT_FDR
    ) %>%
    dplyr::arrange(adj.P.Val)

  DEP_aptamer_STRICT <- DEP_aptamer %>%
    dplyr::filter(
      adj.P.Val < STRICT_FDR
    ) %>%
    dplyr::arrange(adj.P.Val)

  export_dep_table(
    DEP_gene_STRICT,
    file.path(
      fold_dir,
      "DEP_gene_FDR005.csv"
    )
  )

  export_dep_table(
    DEP_aptamer_STRICT,
    file.path(
      fold_dir,
      "DEP_aptamer_FDR005.csv"
    )
  )

  safe_write_csv(
    tibble::tibble(
      EntrezGeneSymbol =
        unique(
          DEP_gene_STRICT$EntrezGeneSymbol
        )
    ),
    file.path(
      fold_dir,
      "candidate_gene_symbols.csv"
    )
  )

  safe_write_csv(
    tibble::tibble(
      AptName =
        unique(
          DEP_aptamer_STRICT$AptName
        )
    ),
    file.path(
      fold_dir,
      "candidate_aptamers.csv"
    )
  )

  safe_write_csv(
    tibble::tibble(
      SampleId =
        main_fit$metadata$SampleId
    ),
    file.path(
      fold_dir,
      "train_ids_used.csv"
    )
  )

  nested_summary[[country]] <- tibble::tibble(
    Country = country,
    N_train = nrow(dep_df_fold),
    CN = sum(dep_df_fold$SampleGroup == "CN"),
    AD = sum(dep_df_fold$SampleGroup == "AD"),
    DEP_gene_STRICT = nrow(DEP_gene_STRICT),
    DEP_aptamer_STRICT = nrow(DEP_aptamer_STRICT)
  )

  message(
    country,
    " | Train=", nrow(dep_df_fold),
    " | Genes=", nrow(DEP_gene_STRICT),
    " | Aptamers=", nrow(DEP_aptamer_STRICT)
  )
}

nested_summary_df <- dplyr::bind_rows(
  nested_summary
)

safe_write_csv(
  nested_summary_df,
  file.path(
    OUT_DIR,
    "nested_DEP_summary.csv"
  )
)

writeLines(
  capture.output(sessionInfo()),
  file.path(
    OUT_DIR,
    "sessionInfo_nested_DEP.txt"
  )
)

message(
  "\n====================================\n",
  "Nested DEP completed.\n",
  "Outputs written to: ",
  OUT_DIR,
  "\n===================================="
)

warnings()
###############################################################################

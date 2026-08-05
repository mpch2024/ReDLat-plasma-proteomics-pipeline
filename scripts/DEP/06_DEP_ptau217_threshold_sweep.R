###############################################################################
# ReDLat plasma proteomics — DEP workflow
# 06. p-tau217 threshold sensitivity analysis
# Requires: Outputs from Scripts 01–02
# Produces: Threshold-sweep models, summaries and diagnostic plots
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

# =============================================================================
# 0) Paths, packages and user-editable constants
# =============================================================================


PROJECT_ROOT <- project_root
OUTDIR <- project_root

MAIN_FDR <- 0.05
STRICT_FDR <- 0.01
SEED_GLOBAL <- 1234
set.seed(SEED_GLOBAL)

# Cohort-derived threshold sweep. These are NOT clinical positivity cutoffs.
PTAU217_CN_QUANTILES <- c(0.80, 0.90, 0.95)
MAIN_TEXT_PTAU217_QUANTILE <- 0.90

# Optional stringent secondary contrast:
# AD: high p-tau217 and low Aβ42/40
# CN: below p-tau217 threshold and above low Aβ42/40 threshold
# This is expected to reduce sample size and attenuate signal.
RUN_STRINGENT_PTAU217_ABETA_CONTRAST <- TRUE
ABETA_CN_LOW_QUANTILE <- 0.10

# Minimum group size required to fit a contrast.
MIN_N_PER_GROUP <- 20

# Expression handling.
# Use "auto" unless you are certain your workspace stores raw RFU or log2 RFU.
# Options: "auto", "log2", "none".
EXPRESSION_SCALE_MODE <- "auto"
RUN_BASELINE_REPLICATION_CHECK <- TRUE
BASELINE_LOGFC_COR_MIN_WARNING <- 0.995

required_pkgs <- c(
  "dplyr", "tidyr", "purrr", "tibble", "readr", "stringr",
  "limma", "openxlsx", "ggplot2", "scales"
)
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
}
invisible(lapply(required_pkgs, library, character.only = TRUE))
options(stringsAsFactors = FALSE)

first_existing_file <- function(paths) {
  hit <- paths[file.exists(paths)][1]
  if (length(hit) == 0 || is.na(hit)) return(NA_character_)
  hit
}

workspace_file <- first_existing_file(c(
  file.path(analysis_root, "workspace", "analysis_workspace.RData")
))

if (is.na(workspace_file)) {
  stop(
    "Could not find analysis_workspace.RData. Run 01_DEP_primary_analysis.R first.\n",
    "Expected location: result/workspace/analysis_workspace.RData",
    call. = FALSE
  )
}
load(workspace_file)

OUT_ROOT <- file.path(analysis_root, "04_sensitivity", "p_tau217_enrichment")
TABLE_DIR <- file.path(OUT_ROOT, "tables")
FIG_DIR <- file.path(OUT_ROOT, "figures")
LOG_DIR <- file.path(OUT_ROOT, "logs")
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

safe_write_csv <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(x, file)
}

save_xlsx_safe <- function(named_list, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  wb <- openxlsx::createWorkbook()
  for (nm in names(named_list)) {
    sheet <- substr(nm, 1, 31)
    openxlsx::addWorksheet(wb, sheet)
    openxlsx::writeData(wb, sheet, named_list[[nm]])
  }
  openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
}

# =============================================================================
# 1) Helper functions
# =============================================================================

safe_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))

as_numeric_matrix <- function(x) {
  x <- as.data.frame(x, check.names = FALSE)
  x[] <- lapply(x, safe_numeric)
  as.matrix(x)
}

infer_expression_scale <- function(dat, protein_cols, sample_n = 200000) {
  protein_cols <- intersect(protein_cols, names(dat))
  if (length(protein_cols) == 0) stop("No protein columns available for expression-scale inference.", call. = FALSE)

  probe_cols <- protein_cols[seq_len(min(length(protein_cols), 500))]
  vals <- as.vector(as_numeric_matrix(dat[, probe_cols, drop = FALSE]))
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0) stop("Expression matrix contains no finite values.", call. = FALSE)
  if (length(vals) > sample_n) vals <- sample(vals, sample_n)

  qs <- stats::quantile(vals, probs = c(0.01, 0.50, 0.99), na.rm = TRUE, names = FALSE)
  names(qs) <- c("q01", "q50", "q99")

  # Heuristic:
  # Raw SomaScan RFU values are usually far above the log2 scale.
  # log2 RFU commonly has q50 around 6-16 and q99 rarely above ~25.
  decision <- if (EXPRESSION_SCALE_MODE == "log2") {
    "log2"
  } else if (EXPRESSION_SCALE_MODE == "none") {
    "none"
  } else {
    if (is.finite(qs[["q99"]]) && (qs[["q99"]] > 50 || qs[["q50"]] > 25)) "log2" else "none"
  }

  tibble::tibble(
    expression_scale_mode_requested = EXPRESSION_SCALE_MODE,
    expression_transform_applied = decision,
    q01 = qs[["q01"]],
    q50 = qs[["q50"]],
    q99 = qs[["q99"]],
    interpretation = dplyr::case_when(
      decision == "log2" ~ "Values appear to be raw/high-scale RFU; log2 transform will be applied.",
      decision == "none" ~ "Values appear already log2-like or standardized; no additional log2 transform will be applied.",
      TRUE ~ "Unknown"
    )
  )
}

make_expression_matrix <- function(dat, protein_cols, transform_decision) {
  protein_cols <- intersect(protein_cols, names(dat))
  expr <- as_numeric_matrix(dat[, protein_cols, drop = FALSE])

  if (transform_decision == "log2") {
    expr[expr <= 0] <- NA_real_
    expr <- log2(expr)
  } else if (transform_decision == "none") {
    expr <- expr
  } else {
    stop("Invalid transform_decision: ", transform_decision, call. = FALSE)
  }

  # limma expects features x samples.
  t(expr)
}

safe_se_from_limma <- function(logFC, t) {
  out <- abs(safe_numeric(logFC) / safe_numeric(t))
  out[!is.finite(out)] <- NA_real_
  out
}

ensure_annot_cols <- function(annot_tbl) {
  needed <- c("AptName", "EntrezGeneSymbol", "TargetFullName", "Target", "UniProt", "SeqId")
  for (nm in needed) {
    if (!nm %in% names(annot_tbl)) annot_tbl[[nm]] <- NA_character_
  }
  annot_tbl
}

make_protein_name <- function(symbol, target_full, target, apt, seqid = NA) {
  dplyr::coalesce(
    as.character(symbol),
    as.character(target_full),
    as.character(target),
    as.character(seqid),
    as.character(apt)
  )
}

classify_dep <- function(dep_tbl, fdr = MAIN_FDR) {
  dep_tbl %>%
    dplyr::mutate(
      type = dplyr::case_when(
        logFC > 0 & adj.P.Val < fdr ~ "Higher_in_enriched_AD",
        logFC < 0 & adj.P.Val < fdr ~ "Lower_in_enriched_AD",
        TRUE ~ "NS"
      ),
      Direction = dplyr::case_when(
        type == "Higher_in_enriched_AD" ~ "Higher in p-tau217-enriched clinically diagnosed AD",
        type == "Lower_in_enriched_AD" ~ "Lower in p-tau217-enriched clinically diagnosed AD",
        TRUE ~ "Not significant"
      )
    )
}

build_primary_gene_map <- function(primary_gene_tbl) {
  required <- c("EntrezGeneSymbol", "AptName")
  missing_required <- setdiff(required, names(primary_gene_tbl))
  if (length(missing_required) > 0) {
    stop("DEP_gene is missing: ", paste(missing_required, collapse = ", "), call. = FALSE)
  }
  primary_map <- primary_gene_tbl %>%
    dplyr::transmute(
      EntrezGeneSymbol = as.character(EntrezGeneSymbol),
      AptName = as.character(AptName)
    ) %>%
    dplyr::filter(
      !is.na(EntrezGeneSymbol), EntrezGeneSymbol != "",
      !is.na(AptName), AptName != ""
    ) %>%
    dplyr::distinct()

  if (anyDuplicated(primary_map$EntrezGeneSymbol) > 0 ||
      anyDuplicated(primary_map$AptName) > 0 ||
      nrow(primary_map) != nrow(primary_gene_tbl)) {
    stop("Primary DEP_gene does not define a unique one-gene/one-SOMAmer map.", call. = FALSE)
  }
  primary_map
}

apply_primary_somamer_map <- function(dep_tbl, primary_gene_map, model_name = "sensitivity_model") {
  out <- dep_tbl %>%
    dplyr::mutate(
      EntrezGeneSymbol = as.character(EntrezGeneSymbol),
      AptName = as.character(AptName)
    ) %>%
    dplyr::inner_join(primary_gene_map, by = c("EntrezGeneSymbol", "AptName")) %>%
    dplyr::arrange(match(EntrezGeneSymbol, primary_gene_map$EntrezGeneSymbol))

  missing_map <- dplyr::anti_join(
    primary_gene_map,
    out %>% dplyr::select(EntrezGeneSymbol, AptName),
    by = c("EntrezGeneSymbol", "AptName")
  )
  if (nrow(missing_map) > 0) {
    stop(
      model_name, " recovered ", nrow(out), " of ", nrow(primary_gene_map),
      " fixed primary SOMAmers. First missing pairs: ",
      paste(utils::head(paste0(missing_map$EntrezGeneSymbol, "=", missing_map$AptName), 10), collapse = "; "),
      call. = FALSE
    )
  }
  if (anyDuplicated(out$EntrezGeneSymbol) > 0 || anyDuplicated(out$AptName) > 0) {
    stop(model_name, " has duplicated genes or SOMAmers after fixed-map selection.", call. = FALSE)
  }
  out
}

export_dep_table <- function(dep_tbl, file) {
  keep_cols <- intersect(
    c(
      "Protein_Name", "EntrezGeneSymbol", "TargetFullName", "Target", "UniProt", "AptName", "SeqId", "feature_id_raw",
      "logFC", "se", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "type", "Direction"
    ),
    names(dep_tbl)
  )
  dep_tbl %>%
    dplyr::select(dplyr::all_of(keep_cols)) %>%
    dplyr::arrange(adj.P.Val) %>%
    safe_write_csv(file)
}

run_limma_dep_model_local <- function(dat, seq_cols, annot_tbl, formula_str, coef_name, model_name, transform_decision) {
  protein_cols <- intersect(seq_cols, names(dat))
  if (length(protein_cols) == 0) stop("No protein columns found for ", model_name, call. = FALSE)

  metadata <- dat %>% dplyr::select(-dplyr::all_of(protein_cols))
  model_vars <- all.vars(stats::as.formula(formula_str))
  missing_model_vars <- setdiff(model_vars, names(metadata))
  if (length(missing_model_vars) > 0) {
    stop("Missing model variables in ", model_name, ": ", paste(missing_model_vars, collapse = ", "), call. = FALSE)
  }

  keep <- stats::complete.cases(metadata[, model_vars, drop = FALSE])
  dat <- dat[keep, , drop = FALSE]
  metadata <- metadata[keep, , drop = FALSE]

  expr <- make_expression_matrix(dat, protein_cols, transform_decision = transform_decision)
  design <- model.matrix(stats::as.formula(formula_str), data = metadata)

  if (!coef_name %in% colnames(design)) {
    stop(
      "Coefficient not found in design matrix for ", model_name, ": ", coef_name,
      "\nAvailable coefficients: ", paste(colnames(design), collapse = ", "), call. = FALSE
    )
  }

  fit <- limma::lmFit(expr, design)
  fit <- limma::eBayes(fit)

  annot_tbl <- ensure_annot_cols(annot_tbl)

  tt <- limma::topTable(fit, coef = coef_name, adjust.method = "BH", number = Inf) %>%
    tibble::rownames_to_column(var = "feature_id_raw") %>%
    dplyr::mutate(AptName = as.character(feature_id_raw), se = safe_se_from_limma(logFC, t)) %>%
    dplyr::left_join(annot_tbl, by = "AptName") %>%
    dplyr::mutate(Protein_Name = make_protein_name(EntrezGeneSymbol, TargetFullName, Target, AptName, feature_id_raw)) %>%
    classify_dep(fdr = MAIN_FDR)

  effective_counts <- metadata %>%
    dplyr::count(dplyr::across(dplyr::all_of(all.vars(stats::as.formula(formula_str))[1])), name = "n_effective")

  list(
    dep = tt,
    fit = fit,
    design = design,
    metadata = metadata,
    protein_cols = protein_cols,
    effective_counts = effective_counts
  )
}

pick_col <- function(df, candidates) {
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

build_gene_compare_table <- function(primary_gene, secondary_gene, label) {
  primary_std <- primary_gene %>%
    dplyr::select(EntrezGeneSymbol, AptName, Protein_Name, logFC, adj.P.Val, P.Value, type) %>%
    dplyr::rename(
      Protein_Name_primary = Protein_Name,
      logFC_primary = logFC,
      adj.P.Val_primary = adj.P.Val,
      P.Value_primary = P.Value,
      type_primary = type
    )

  secondary_std <- secondary_gene %>%
    dplyr::select(EntrezGeneSymbol, AptName, Protein_Name, logFC, adj.P.Val, P.Value, type) %>%
    dplyr::rename(
      Protein_Name_secondary = Protein_Name,
      logFC_secondary = logFC,
      adj.P.Val_secondary = adj.P.Val,
      P.Value_secondary = P.Value,
      type_secondary = type
    )

  out <- primary_std %>%
    dplyr::inner_join(
      secondary_std,
      by = c("EntrezGeneSymbol", "AptName")
    ) %>%
    dplyr::mutate(
      AptName_primary = AptName,
      AptName_secondary = AptName,
      comparison = label,
      Protein_Name = dplyr::coalesce(Protein_Name_primary, Protein_Name_secondary),
      same_somamer = TRUE,
      same_direction = sign(logFC_primary) == sign(logFC_secondary),
      primary_fdr005 = adj.P.Val_primary < MAIN_FDR,
      secondary_fdr005 = adj.P.Val_secondary < MAIN_FDR,
      primary_fdr001 = adj.P.Val_primary < STRICT_FDR,
      secondary_fdr001 = adj.P.Val_secondary < STRICT_FDR,
      primary_fdr005_preserved_same_direction = primary_fdr005 & secondary_fdr005 & same_direction,
      primary_fdr005_same_direction_any_secondary = primary_fdr005 & same_direction
    )

  if (nrow(out) != nrow(primary_gene_map)) {
    stop(label, " comparison contains ", nrow(out), " rows; expected ",
         nrow(primary_gene_map), ".", call. = FALSE)
  }
  out
}

get_effective_n <- function(fit_obj, group_col = "BioContrast") {
  fit_obj$metadata %>%
    dplyr::count(.data[[group_col]], name = "n_effective") %>%
    tidyr::pivot_wider(names_from = .data[[group_col]], values_from = n_effective, values_fill = 0)
}

summarize_model <- function(dep_gene, compare_tbl, model_name, initial_counts, effective_counts,
                            threshold_family, ptau_q, ptau_cut, abeta_q, abeta_cut, status = "completed") {
  n_init_cn <- initial_counts$n_initial[initial_counts$BioContrast == "CN_below_ptau217_threshold"]
  n_init_ad <- initial_counts$n_initial[initial_counts$BioContrast == "AD_ptau217_enriched"]
  n_eff_cn <- effective_counts$CN_below_ptau217_threshold[1]
  n_eff_ad <- effective_counts$AD_ptau217_enriched[1]

  if (length(n_init_cn) == 0) n_init_cn <- 0
  if (length(n_init_ad) == 0) n_init_ad <- 0
  if (length(n_eff_cn) == 0 || is.na(n_eff_cn)) n_eff_cn <- 0
  if (length(n_eff_ad) == 0 || is.na(n_eff_ad)) n_eff_ad <- 0

  tibble::tibble(
    model = model_name,
    threshold_family = threshold_family,
    p_tau217_CN_quantile = ptau_q,
    p_tau217_cutoff = ptau_cut,
    abeta42_40_CN_low_quantile = abeta_q,
    abeta42_40_low_cutoff = abeta_cut,
    status = status,
    n_CN_initial = as.integer(n_init_cn),
    n_AD_initial = as.integer(n_init_ad),
    n_total_initial = as.integer(n_init_cn + n_init_ad),
    n_CN_effective_model = as.integer(n_eff_cn),
    n_AD_effective_model = as.integer(n_eff_ad),
    n_total_effective_model = as.integer(n_eff_cn + n_eff_ad),
    genes_tested = nrow(dep_gene),
    fdr005_total = sum(dep_gene$adj.P.Val < MAIN_FDR, na.rm = TRUE),
    fdr005_higher = sum(dep_gene$adj.P.Val < MAIN_FDR & dep_gene$logFC > 0, na.rm = TRUE),
    fdr005_lower = sum(dep_gene$adj.P.Val < MAIN_FDR & dep_gene$logFC < 0, na.rm = TRUE),
    fdr001_total = sum(dep_gene$adj.P.Val < STRICT_FDR, na.rm = TRUE),
    fdr001_higher = sum(dep_gene$adj.P.Val < STRICT_FDR & dep_gene$logFC > 0, na.rm = TRUE),
    fdr001_lower = sum(dep_gene$adj.P.Val < STRICT_FDR & dep_gene$logFC < 0, na.rm = TRUE),
    logFC_correlation_with_primary = suppressWarnings(stats::cor(compare_tbl$logFC_primary, compare_tbl$logFC_secondary, use = "complete.obs")),
    direction_consistency_all = mean(compare_tbl$same_direction, na.rm = TRUE),
    direction_consistency_primary_fdr005 = mean(compare_tbl$same_direction[compare_tbl$primary_fdr005], na.rm = TRUE),
    primary_fdr005_same_direction_any_secondary_n = sum(compare_tbl$primary_fdr005_same_direction_any_secondary, na.rm = TRUE),
    primary_fdr005_preserved_same_direction_n = sum(compare_tbl$primary_fdr005_preserved_same_direction, na.rm = TRUE),
    primary_fdr005_n = sum(compare_tbl$primary_fdr005, na.rm = TRUE),
    preservation_fraction_primary_fdr005 = primary_fdr005_preserved_same_direction_n / pmax(primary_fdr005_n, 1),
    interpretation = "Cohort-defined exploratory p-tau217 enrichment; not a diagnostic cutoff or biomarker-confirmed AD classification."
  )
}

make_clean_model_id <- function(family, q) {
  paste0(family, "_CNp", sprintf("%02d", round(q * 100)))
}

format_percent <- function(x, digits = 1) paste0(round(100 * x, digits), "%")

# =============================================================================
# 2) Validate workspace objects and detect biomarker columns
# =============================================================================

required_objects <- c("dep_df", "dep_protein_cols", "annot_tbl", "DEP_gene")
missing_objects <- required_objects[!vapply(required_objects, exists, logical(1), envir = .GlobalEnv)]
if (length(missing_objects) > 0) {
  stop("Workspace is missing required objects: ", paste(missing_objects, collapse = ", "), call. = FALSE)
}

annot_tbl <- ensure_annot_cols(annot_tbl)
primary_gene_map <- build_primary_gene_map(DEP_gene)
safe_write_csv(primary_gene_map, file.path(LOG_DIR, "fixed_primary_SOMAmer_gene_map.csv"))
safe_write_csv(
  tibble::tibble(
    item = c("workspace_file", "primary_map_rows", "fixed_map_rule"),
    value = c(
      workspace_file,
      as.character(nrow(primary_gene_map)),
      "Exact EntrezGeneSymbol--AptName pairs inherited from primary DEP_gene"
    )
  ),
  file.path(LOG_DIR, "fixed_primary_SOMAmer_gene_map_audit.csv")
)

ptau217_col <- pick_col(dep_df, c("p.tau217", "p_tau217", "p-tau217", "ptau217", "Ptau217", "pTau217", "pTau_217", "plasma_p_tau217"))
ratio_col <- pick_col(dep_df, c("ratio.AB42.40", "ratio_AB42_40", "ratio AB42/40", "AB42_40_ratio", "Aβ42/40", "Abeta42_40", "Aβ42_40", "AB42.AB40"))

if (is.na(ptau217_col)) {
  stop("Could not detect p-tau217 column in dep_df. Please inspect column names.", call. = FALSE)
}

biomarker_detection <- tibble::tibble(
  biomarker = c("p_tau217", "Aβ42/40_ratio"),
  detected_column = c(ptau217_col, ratio_col),
  use_in_script = c(TRUE, !is.na(ratio_col) && RUN_STRINGENT_PTAU217_ABETA_CONTRAST)
)
safe_write_csv(biomarker_detection, file.path(LOG_DIR, "biomarker_columns_detected.csv"))

expression_scale_audit <- infer_expression_scale(dep_df, dep_protein_cols)
TRANSFORM_DECISION <- expression_scale_audit$expression_transform_applied[1]
safe_write_csv(expression_scale_audit, file.path(LOG_DIR, "expression_scale_audit.csv"))
message("Expression transform decision: ", TRANSFORM_DECISION)

analysis_df <- dep_df %>%
  dplyr::mutate(
    SampleGroup = factor(SampleGroup, levels = c("CN", "AD")),
    Sex = factor(Sex),
    Country = factor(Country),
    p_tau217_numeric = safe_numeric(.data[[ptau217_col]]),
    ratio_ab42_40_numeric = if (!is.na(ratio_col)) safe_numeric(.data[[ratio_col]]) else NA_real_
  )

biomarker_missingness <- analysis_df %>%
  dplyr::group_by(SampleGroup, Country) %>%
  dplyr::summarise(
    n = dplyr::n(),
    p_tau217_available = sum(!is.na(p_tau217_numeric)),
    p_tau217_missing = sum(is.na(p_tau217_numeric)),
    p_tau217_missing_fraction = p_tau217_missing / n,
    abeta42_40_available = sum(!is.na(ratio_ab42_40_numeric)),
    abeta42_40_missing = sum(is.na(ratio_ab42_40_numeric)),
    abeta42_40_missing_fraction = abeta42_40_missing / n,
    .groups = "drop"
  )
safe_write_csv(biomarker_missingness, file.path(TABLE_DIR, "biomarker_missingness_by_group_country.csv"))

cn_ref <- analysis_df %>% dplyr::filter(SampleGroup == "CN")

ptau_threshold_tbl <- tibble::tibble(
  biomarker = "p_tau217",
  quantile = PTAU217_CN_QUANTILES,
  cutoff = as.numeric(stats::quantile(cn_ref$p_tau217_numeric, probs = PTAU217_CN_QUANTILES, na.rm = TRUE, names = FALSE)),
  reference_group = "CN",
  threshold_role = "cohort_defined_p_tau217_enrichment_threshold",
  diagnostic_interpretation = "Not a clinical positivity cutoff; no PET/CSF reference standard available.",
  analysis_interpretation = "Used only to enrich the clinical AD versus CN contrast and test preservation of the primary proteomic signature."
)

ratio_low_cutoff <- if (!is.na(ratio_col)) {
  as.numeric(stats::quantile(cn_ref$ratio_ab42_40_numeric, probs = ABETA_CN_LOW_QUANTILE, na.rm = TRUE, names = FALSE))
} else {
  NA_real_
}

abeta_threshold_tbl <- tibble::tibble(
  biomarker = "Aβ42/40_ratio",
  quantile = ABETA_CN_LOW_QUANTILE,
  cutoff = ratio_low_cutoff,
  reference_group = "CN",
  threshold_role = "cohort_defined_low_abeta42_40_threshold",
  diagnostic_interpretation = "Not a clinical positivity cutoff; used only for stringent secondary sensitivity.",
  analysis_interpretation = "Lower Aβ42/40 values are treated as more AD-compatible only for exploratory enrichment."
)

threshold_tbl <- dplyr::bind_rows(ptau_threshold_tbl, abeta_threshold_tbl)
safe_write_csv(threshold_tbl, file.path(TABLE_DIR, "cohort_defined_biomarker_thresholds_v3.csv"))

threshold_membership <- purrr::map_dfr(seq_len(nrow(ptau_threshold_tbl)), function(i) {
  q <- ptau_threshold_tbl$quantile[i]
  cut <- ptau_threshold_tbl$cutoff[i]
  analysis_df %>%
    dplyr::mutate(
      p_tau217_threshold_quantile = q,
      p_tau217_cutoff = cut,
      p_tau217_position = dplyr::case_when(
        is.na(p_tau217_numeric) ~ "missing_p_tau217",
        p_tau217_numeric >= cut ~ "at_or_above_threshold",
        p_tau217_numeric < cut ~ "below_threshold",
        TRUE ~ "unclassified"
      )
    ) %>%
    dplyr::count(p_tau217_threshold_quantile, p_tau217_cutoff, SampleGroup, p_tau217_position, name = "n")
})
safe_write_csv(threshold_membership, file.path(TABLE_DIR, "p_tau217_threshold_membership_audit.csv"))

# =============================================================================
# 3) Optional baseline replication check against primary DEP_gene
# =============================================================================

baseline_check <- tibble::tibble()

if (RUN_BASELINE_REPLICATION_CHECK) {
  message("Running baseline replication check against primary DEP_gene...")

  baseline_formula <- "~ SampleGroup + Age + Sex + Country + Education"
  baseline_fit <- run_limma_dep_model_local(
    dat = analysis_df,
    seq_cols = dep_protein_cols,
    annot_tbl = annot_tbl,
    formula_str = baseline_formula,
    coef_name = "SampleGroupAD",
    model_name = "baseline_replication_check",
    transform_decision = TRANSFORM_DECISION
  )

  baseline_dep_gene <- apply_primary_somamer_map(
    baseline_fit$dep %>% dplyr::filter(AptName %in% dep_protein_cols),
    primary_gene_map,
    model_name = "baseline_replication_check"
  )
  baseline_compare <- build_gene_compare_table(DEP_gene, baseline_dep_gene, "baseline_replication_check")

  baseline_check <- tibble::tibble(
    transform_decision = TRANSFORM_DECISION,
    baseline_formula = baseline_formula,
    genes_compared = nrow(baseline_compare),
    logFC_correlation_with_workspace_DEP_gene = suppressWarnings(stats::cor(baseline_compare$logFC_primary, baseline_compare$logFC_secondary, use = "complete.obs")),
    direction_consistency_all = mean(baseline_compare$same_direction, na.rm = TRUE),
    direction_consistency_primary_fdr005 = mean(baseline_compare$same_direction[baseline_compare$primary_fdr005], na.rm = TRUE),
    warning = dplyr::if_else(
      logFC_correlation_with_workspace_DEP_gene < BASELINE_LOGFC_COR_MIN_WARNING,
      paste0("Baseline logFC correlation is below ", BASELINE_LOGFC_COR_MIN_WARNING, ". Check whether expression values were transformed consistently with the primary DEP workflow."),
      "Baseline replication check passed."
    )
  )

  safe_write_csv(baseline_check, file.path(LOG_DIR, "baseline_replication_check.csv"))
  if (baseline_check$logFC_correlation_with_workspace_DEP_gene[1] < BASELINE_LOGFC_COR_MIN_WARNING) {
    warning(baseline_check$warning[1])
  }
}

# =============================================================================
# 4) Define p-tau217-enriched contrasts
# =============================================================================

make_tau_only_contrast <- function(df, ptau_cutoff) {
  df %>%
    dplyr::mutate(
      BioContrast = dplyr::case_when(
        SampleGroup == "AD" & !is.na(p_tau217_numeric) & p_tau217_numeric >= ptau_cutoff ~ "AD_ptau217_enriched",
        SampleGroup == "CN" & !is.na(p_tau217_numeric) & p_tau217_numeric < ptau_cutoff ~ "CN_below_ptau217_threshold",
        TRUE ~ NA_character_
      ),
      excluded_from_enriched_contrast_reason = dplyr::case_when(
        is.na(p_tau217_numeric) ~ "missing_p_tau217",
        SampleGroup == "AD" & p_tau217_numeric < ptau_cutoff ~ "clinically_diagnosed_AD_below_cohort_defined_threshold",
        SampleGroup == "CN" & p_tau217_numeric >= ptau_cutoff ~ "CN_at_or_above_cohort_defined_threshold",
        TRUE ~ NA_character_
      )
    )
}

make_tau_abeta_contrast <- function(df, ptau_cutoff, abeta_cutoff) {
  if (is.na(ratio_col) || !is.finite(abeta_cutoff)) return(NULL)
  df %>%
    dplyr::mutate(
      BioContrast = dplyr::case_when(
        SampleGroup == "AD" &
          !is.na(p_tau217_numeric) & !is.na(ratio_ab42_40_numeric) &
          p_tau217_numeric >= ptau_cutoff & ratio_ab42_40_numeric <= abeta_cutoff ~ "AD_ptau217_enriched",
        SampleGroup == "CN" &
          !is.na(p_tau217_numeric) & !is.na(ratio_ab42_40_numeric) &
          p_tau217_numeric < ptau_cutoff & ratio_ab42_40_numeric > abeta_cutoff ~ "CN_below_ptau217_threshold",
        TRUE ~ NA_character_
      ),
      excluded_from_enriched_contrast_reason = dplyr::case_when(
        is.na(p_tau217_numeric) ~ "missing_p_tau217",
        is.na(ratio_ab42_40_numeric) ~ "missing_abeta42_40_ratio",
        SampleGroup == "AD" & p_tau217_numeric < ptau_cutoff ~ "clinically_diagnosed_AD_below_cohort_defined_p_tau217_threshold",
        SampleGroup == "AD" & ratio_ab42_40_numeric > abeta_cutoff ~ "clinically_diagnosed_AD_not_low_abeta42_40_by_cohort_threshold",
        SampleGroup == "CN" & p_tau217_numeric >= ptau_cutoff ~ "CN_at_or_above_cohort_defined_p_tau217_threshold",
        SampleGroup == "CN" & ratio_ab42_40_numeric <= abeta_cutoff ~ "CN_low_abeta42_40_by_cohort_threshold",
        TRUE ~ NA_character_
      )
    )
}

# =============================================================================
# 5) Run threshold-sweep models
# =============================================================================

all_summary <- list()
all_counts <- list()
all_exclusions <- list()
all_compare <- list()
all_key_dep <- list()

for (ii in seq_len(nrow(ptau_threshold_tbl))) {
  ptau_q <- ptau_threshold_tbl$quantile[ii]
  ptau_cutoff <- ptau_threshold_tbl$cutoff[ii]

  contrast_specs <- list(
    list(
      family = "p_tau217_enriched",
      model_id = make_clean_model_id("p_tau217_enriched", ptau_q),
      dat = make_tau_only_contrast(analysis_df, ptau_cutoff),
      abeta_q = NA_real_,
      abeta_cutoff = NA_real_,
      manuscript_role = ifelse(abs(ptau_q - MAIN_TEXT_PTAU217_QUANTILE) < 1e-9, "main_text_candidate", "threshold_sensitivity")
    )
  )

  if (RUN_STRINGENT_PTAU217_ABETA_CONTRAST && !is.na(ratio_col) && is.finite(ratio_low_cutoff)) {
    contrast_specs[[length(contrast_specs) + 1]] <- list(
      family = "p_tau217_enriched_plus_low_abeta42_40",
      model_id = make_clean_model_id("p_tau217_enriched_low_abeta42_40", ptau_q),
      dat = make_tau_abeta_contrast(analysis_df, ptau_cutoff, ratio_low_cutoff),
      abeta_q = ABETA_CN_LOW_QUANTILE,
      abeta_cutoff = ratio_low_cutoff,
      manuscript_role = "stringent_secondary_sensitivity"
    )
  }

  for (spec in contrast_specs) {
    contrast_name <- spec$model_id
    message("Running exploratory enriched contrast: ", contrast_name)

    dat_model <- spec$dat %>%
      dplyr::filter(!is.na(BioContrast)) %>%
      dplyr::mutate(
        BioContrast = factor(BioContrast, levels = c("CN_below_ptau217_threshold", "AD_ptau217_enriched")),
        Sex = factor(Sex),
        Country = droplevels(factor(Country))
      )

    group_counts <- dat_model %>%
      dplyr::count(BioContrast, name = "n_initial") %>%
      dplyr::mutate(
        model = contrast_name,
        threshold_family = spec$family,
        manuscript_role = spec$manuscript_role,
        p_tau217_CN_quantile = ptau_q,
        p_tau217_cutoff = ptau_cutoff,
        abeta42_40_CN_low_quantile = spec$abeta_q,
        abeta42_40_low_cutoff = spec$abeta_cutoff,
        .before = 1
      )

    exclusion_counts <- spec$dat %>%
      dplyr::filter(is.na(BioContrast)) %>%
      dplyr::count(SampleGroup, excluded_from_enriched_contrast_reason, name = "n_excluded") %>%
      dplyr::mutate(
        model = contrast_name,
        threshold_family = spec$family,
        p_tau217_CN_quantile = ptau_q,
        p_tau217_cutoff = ptau_cutoff,
        .before = 1
      )

    safe_write_csv(group_counts, file.path(TABLE_DIR, paste0(contrast_name, "_initial_sample_counts.csv")))
    safe_write_csv(exclusion_counts, file.path(TABLE_DIR, paste0(contrast_name, "_excluded_sample_counts.csv")))
    all_counts[[contrast_name]] <- group_counts
    all_exclusions[[contrast_name]] <- exclusion_counts

    n_ad <- sum(dat_model$BioContrast == "AD_ptau217_enriched", na.rm = TRUE)
    n_cn <- sum(dat_model$BioContrast == "CN_below_ptau217_threshold", na.rm = TRUE)

    if (n_ad < MIN_N_PER_GROUP || n_cn < MIN_N_PER_GROUP) {
      warning("Skipping ", contrast_name, ": fewer than ", MIN_N_PER_GROUP, " samples in one contrast group.")
      all_summary[[contrast_name]] <- tibble::tibble(
        model = contrast_name,
        threshold_family = spec$family,
        manuscript_role = spec$manuscript_role,
        p_tau217_CN_quantile = ptau_q,
        p_tau217_cutoff = ptau_cutoff,
        abeta42_40_CN_low_quantile = spec$abeta_q,
        abeta42_40_low_cutoff = spec$abeta_cutoff,
        status = "skipped_low_n",
        n_CN_initial = n_cn,
        n_AD_initial = n_ad,
        n_total_initial = n_cn + n_ad
      )
      next
    }

    formula_str <- "~ BioContrast + Age + Sex + Country + Education"
    fit_obj <- run_limma_dep_model_local(
      dat = dat_model,
      seq_cols = dep_protein_cols,
      annot_tbl = annot_tbl,
      formula_str = formula_str,
      coef_name = "BioContrastAD_ptau217_enriched",
      model_name = contrast_name,
      transform_decision = TRANSFORM_DECISION
    )

    dep_apt <- fit_obj$dep %>% dplyr::filter(AptName %in% dep_protein_cols)
    dep_gene <- apply_primary_somamer_map(dep_apt, primary_gene_map, model_name = contrast_name)

    model_dir <- file.path(TABLE_DIR, contrast_name)
    dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

    export_dep_table(dep_apt, file.path(model_dir, paste0(contrast_name, "_limma_results_aptamer_level.csv")))
    export_dep_table(dep_gene, file.path(model_dir, paste0(contrast_name, "_limma_results_gene_collapsed.csv")))
    export_dep_table(dep_gene %>% dplyr::filter(adj.P.Val < MAIN_FDR), file.path(model_dir, paste0(contrast_name, "_DEP_gene_collapsed_FDR005.csv")))
    export_dep_table(dep_gene %>% dplyr::filter(adj.P.Val < STRICT_FDR), file.path(model_dir, paste0(contrast_name, "_DEP_gene_collapsed_FDR001.csv")))

    compare_tbl <- build_gene_compare_table(DEP_gene, dep_gene, contrast_name)
    safe_write_csv(compare_tbl, file.path(model_dir, paste0("primary_vs_", contrast_name, "_gene_comparison.csv")))

    effective_counts <- get_effective_n(fit_obj, group_col = "BioContrast")

    summary_tbl <- summarize_model(
      dep_gene = dep_gene,
      compare_tbl = compare_tbl,
      model_name = contrast_name,
      initial_counts = group_counts,
      effective_counts = effective_counts,
      threshold_family = spec$family,
      ptau_q = ptau_q,
      ptau_cut = ptau_cutoff,
      abeta_q = spec$abeta_q,
      abeta_cut = spec$abeta_cutoff,
      status = "completed"
    ) %>%
      dplyr::mutate(
        manuscript_role = spec$manuscript_role,
        formula = formula_str,
        coefficient = "BioContrastAD_ptau217_enriched",
        expression_transform_applied = TRANSFORM_DECISION,
        .after = threshold_family
      )

    key_dep <- dep_gene %>%
      dplyr::filter(adj.P.Val < MAIN_FDR) %>%
      dplyr::mutate(
        model = contrast_name,
        threshold_family = spec$family,
        manuscript_role = spec$manuscript_role,
        p_tau217_CN_quantile = ptau_q,
        p_tau217_cutoff = ptau_cutoff,
        abeta42_40_CN_low_quantile = spec$abeta_q,
        abeta42_40_low_cutoff = spec$abeta_cutoff,
        .before = 1
      )

    all_summary[[contrast_name]] <- summary_tbl
    all_compare[[contrast_name]] <- compare_tbl
    all_key_dep[[contrast_name]] <- key_dep
  }
}

summary_all <- dplyr::bind_rows(all_summary) %>%
  dplyr::arrange(threshold_family, p_tau217_CN_quantile)
sample_counts_all <- dplyr::bind_rows(all_counts) %>%
  dplyr::arrange(threshold_family, p_tau217_CN_quantile, BioContrast)
exclusions_all <- dplyr::bind_rows(all_exclusions) %>%
  dplyr::arrange(threshold_family, p_tau217_CN_quantile, SampleGroup, excluded_from_enriched_contrast_reason)
compare_all <- if (length(all_compare) > 0) dplyr::bind_rows(all_compare) else tibble::tibble()
key_dep_all <- if (length(all_key_dep) > 0) dplyr::bind_rows(all_key_dep) else tibble::tibble()

reviewer_interpretation <- summary_all %>%
  dplyr::mutate(
    suggested_reporting = dplyr::case_when(
      threshold_family == "p_tau217_enriched" & abs(p_tau217_CN_quantile - MAIN_TEXT_PTAU217_QUANTILE) < 1e-9 ~ "main_text_sensitivity_candidate",
      threshold_family == "p_tau217_enriched" ~ "supplementary_threshold_sensitivity",
      threshold_family == "p_tau217_enriched_plus_low_abeta42_40" ~ "stringent_secondary_sensitivity",
      TRUE ~ "supplementary"
    ),
    reviewer_interpretation = dplyr::case_when(
      status != "completed" ~ "Skipped because of insufficient sample size.",
      threshold_family == "p_tau217_enriched" ~ "Tests whether the primary clinical AD-associated proteomic signature is retained after empirical p-tau217 enrichment. This is not a biomarker-confirmed AD contrast.",
      threshold_family == "p_tau217_enriched_plus_low_abeta42_40" ~ "More stringent p-tau217/Aβ42/40-compatible contrast. Expected to attenuate signal because of lower n and greater biomarker stringency.",
      TRUE ~ "Reviewer sensitivity analysis."
    )
  )

main_text_candidate <- reviewer_interpretation %>%
  dplyr::filter(threshold_family == "p_tau217_enriched", abs(p_tau217_CN_quantile - MAIN_TEXT_PTAU217_QUANTILE) < 1e-9)

safe_write_csv(summary_all, file.path(TABLE_DIR, "p_tau217_enriched_threshold_sweep_summary_v3.csv"))
safe_write_csv(sample_counts_all, file.path(TABLE_DIR, "p_tau217_enriched_threshold_sweep_sample_counts_v3.csv"))
safe_write_csv(exclusions_all, file.path(TABLE_DIR, "p_tau217_enriched_threshold_sweep_excluded_counts_v3.csv"))
safe_write_csv(compare_all, file.path(TABLE_DIR, "primary_vs_p_tau217_enriched_threshold_sweep_all_comparisons_v3.csv"))
safe_write_csv(key_dep_all, file.path(TABLE_DIR, "p_tau217_enriched_threshold_sweep_DEP_FDR005_all_models_v3.csv"))
safe_write_csv(main_text_candidate, file.path(TABLE_DIR, "recommended_main_text_p_tau217_CNp90_result_v3.csv"))
safe_write_csv(reviewer_interpretation, file.path(TABLE_DIR, "reviewer_interpretation_threshold_sweep_v3.csv"))

# =============================================================================
# 6) Manuscript-ready text generated from p90 result
# =============================================================================

if (nrow(main_text_candidate) == 1 && main_text_candidate$status[1] == "completed") {
  mt <- main_text_candidate[1, ]
  manuscript_text <- c(
    "Suggested Results paragraph:",
    "",
    paste0(
      "As an exploratory biomarker-enriched sensitivity analysis, we repeated the differential-abundance analysis using a cohort-defined p-tau217 enrichment strategy. ",
      "Because amyloid PET, tau PET, CSF biomarkers and neuropathological confirmation were unavailable, this analysis did not reclassify participants as biomarker-confirmed AD and did not apply an externally validated diagnostic p-tau217 cutoff. ",
      "Instead, the ", round(mt$p_tau217_CN_quantile * 100), "th percentile of plasma p-tau217 among cognitively normal participants was used as an empirical enrichment threshold to compare p-tau217-low cognitively normal individuals with p-tau217-enriched clinically diagnosed AD participants."
    ),
    "",
    paste0(
      "This cohort-defined p-tau217-enriched contrast included ", mt$n_CN_effective_model, " cognitively normal individuals and ", mt$n_AD_effective_model, " clinically diagnosed AD participants after complete-case covariate filtering. ",
      "The analysis identified ", mt$fdr005_total, " proteins at FDR < 0.05 and ", mt$fdr001_total, " proteins at FDR < 0.01. ",
      "Effect estimates remained concordant with the primary clinically defined model, with a log2 fold-change correlation of ", round(mt$logFC_correlation_with_primary, 3), ", ",
      format_percent(mt$direction_consistency_all, 1), " overall directional consistency and ",
      format_percent(mt$direction_consistency_primary_fdr005, 1), " directional consistency among proteins significant in the primary analysis. ",
      "In total, ", mt$primary_fdr005_preserved_same_direction_n, " of ", mt$primary_fdr005_n, " primary FDR-significant proteins were preserved at FDR < 0.05 in the same direction. ",
      "These findings support substantial retention of the primary clinical AD-associated proteomic architecture under a p-tau217-enriched contrast, while preserving the interpretation of the study as clinically defined rather than biomarker-confirmed AD."
    )
  )
  writeLines(manuscript_text, file.path(OUT_ROOT, "manuscript_ready_results_p_tau217_CNp90.txt"))
}

# =============================================================================
# 7) Supplementary figures
# =============================================================================

# p-tau217 distribution with CN-derived thresholds.
plot_dist_df <- analysis_df %>%
  dplyr::filter(!is.na(p_tau217_numeric), !is.na(SampleGroup))

if (nrow(plot_dist_df) > 0) {
  p_dist <- ggplot(plot_dist_df, aes(x = p_tau217_numeric, linetype = SampleGroup)) +
    geom_density(linewidth = 0.7, na.rm = TRUE) +
    geom_vline(
      data = ptau_threshold_tbl,
      aes(xintercept = cutoff),
      linewidth = 0.4,
      alpha = 0.7
    ) +
    geom_text(
      data = ptau_threshold_tbl,
      aes(x = cutoff, y = Inf, label = paste0("CN p", round(quantile * 100))),
      angle = 90,
      vjust = 1.2,
      hjust = 1.05,
      size = 3,
      inherit.aes = FALSE
    ) +
    labs(
      title = "Cohort-defined p-tau217 enrichment thresholds",
      subtitle = "Thresholds are derived from the CN distribution and are not diagnostic positivity cutoffs",
      x = "Plasma p-tau217",
      y = "Density",
      linetype = "Clinical group"
    ) +
    theme_classic(base_size = 9)

  ggsave(file.path(FIG_DIR, "p_tau217_distribution_CN_thresholds_v3.pdf"), p_dist, width = 7.0, height = 4.2)
  ggsave(file.path(FIG_DIR, "p_tau217_distribution_CN_thresholds_v3.png"), p_dist, width = 7.0, height = 4.2, dpi = 300)
}

summary_plot_df <- summary_all %>%
  dplyr::filter(status == "completed") %>%
  dplyr::mutate(
    p_tau217_percentile = paste0("CN p", round(p_tau217_CN_quantile * 100)),
    threshold_family = factor(
      threshold_family,
      levels = c("p_tau217_enriched", "p_tau217_enriched_plus_low_abeta42_40")
    )
  )

if (nrow(summary_plot_df) > 0) {
  p_fdr <- ggplot(summary_plot_df, aes(x = p_tau217_percentile, y = fdr005_total, group = threshold_family)) +
    geom_line() +
    geom_point(size = 2) +
    facet_wrap(~ threshold_family, scales = "free_y") +
    labs(
      title = "p-tau217-enriched sensitivity across cohort-defined thresholds",
      subtitle = "Exploratory enrichment only; not biomarker-confirmed diagnosis",
      x = "p-tau217 threshold derived from CN distribution",
      y = "Proteins at FDR < 0.05"
    ) +
    theme_classic(base_size = 9)

  p_cor <- ggplot(summary_plot_df, aes(x = p_tau217_percentile, y = logFC_correlation_with_primary, group = threshold_family)) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_line() +
    geom_point(size = 2) +
    facet_wrap(~ threshold_family, scales = "free_y") +
    labs(
      title = "Concordance with primary clinically defined model",
      x = "p-tau217 threshold derived from CN distribution",
      y = "logFC correlation with primary model"
    ) +
    theme_classic(base_size = 9)

  p_pres <- ggplot(summary_plot_df, aes(x = p_tau217_percentile, y = preservation_fraction_primary_fdr005, group = threshold_family)) +
    geom_line() +
    geom_point(size = 2) +
    facet_wrap(~ threshold_family, scales = "free_y") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
      title = "Preservation of primary FDR-significant proteins",
      x = "p-tau217 threshold derived from CN distribution",
      y = "Primary FDR-significant proteins preserved at FDR < 0.05"
    ) +
    theme_classic(base_size = 9)

  ggsave(file.path(FIG_DIR, "threshold_sweep_FDR005_counts_v3.pdf"), p_fdr, width = 7.0, height = 4.2)
  ggsave(file.path(FIG_DIR, "threshold_sweep_FDR005_counts_v3.png"), p_fdr, width = 7.0, height = 4.2, dpi = 300)
  ggsave(file.path(FIG_DIR, "threshold_sweep_logFC_correlation_v3.pdf"), p_cor, width = 7.0, height = 4.2)
  ggsave(file.path(FIG_DIR, "threshold_sweep_logFC_correlation_v3.png"), p_cor, width = 7.0, height = 4.2, dpi = 300)
  ggsave(file.path(FIG_DIR, "threshold_sweep_primary_DEP_preservation_v3.pdf"), p_pres, width = 7.0, height = 4.2)
  ggsave(file.path(FIG_DIR, "threshold_sweep_primary_DEP_preservation_v3.png"), p_pres, width = 7.0, height = 4.2, dpi = 300)
}

# Primary vs p90 scatter.
p90_id <- make_clean_model_id("p_tau217_enriched", MAIN_TEXT_PTAU217_QUANTILE)
if (p90_id %in% names(all_compare)) {
  p90_compare <- all_compare[[p90_id]] %>%
    dplyr::mutate(
      primary_significant = dplyr::if_else(primary_fdr005, "Primary FDR < 0.05", "Not primary FDR < 0.05")
    )

  p_scatter <- ggplot(p90_compare, aes(x = logFC_primary, y = logFC_secondary, shape = primary_significant)) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    geom_point(alpha = 0.45, size = 1.1, na.rm = TRUE) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, na.rm = TRUE) +
    labs(
      title = "Primary model versus p-tau217-enriched p90 sensitivity",
      subtitle = "Cohort-defined p-tau217 enrichment; not biomarker-confirmed diagnosis",
      x = "Primary clinical model log2 fold-change",
      y = "p-tau217-enriched p90 model log2 fold-change",
      shape = NULL
    ) +
    theme_classic(base_size = 9)

  ggsave(file.path(FIG_DIR, "primary_vs_p_tau217_enriched_CNp90_logFC_scatter_v3.pdf"), p_scatter, width = 5.5, height = 5.0)
  ggsave(file.path(FIG_DIR, "primary_vs_p_tau217_enriched_CNp90_logFC_scatter_v3.png"), p_scatter, width = 5.5, height = 5.0, dpi = 300)
}

# =============================================================================
# 8) Workbook and logs
# =============================================================================

save_xlsx_safe(
  list(
    fixed_primary_map = primary_gene_map,
    thresholds = threshold_tbl,
    biomarker_detection = biomarker_detection,
    expression_scale_audit = expression_scale_audit,
    baseline_check = baseline_check,
    biomarker_missingness = biomarker_missingness,
    threshold_membership = threshold_membership,
    summary = summary_all,
    sample_counts = sample_counts_all,
    excluded_counts = exclusions_all,
    reviewer_interpretation = reviewer_interpretation,
    main_text_candidate = main_text_candidate,
    dep_fdr005_all = key_dep_all,
    comparisons = compare_all
  ),
  file.path(OUT_ROOT, "DEP_p_tau217_enrichment_threshold_sweep_REVIEWER_v3.xlsx")
)

script_metadata <- tibble::tibble(
  item = c(
    "script",
    "workspace_file",
    "output_directory",
    "primary_interpretation",
    "main_text_threshold",
    "expression_transform_applied",
    "date_time"
  ),
  value = c(
    "06_DEP_ptau217_threshold_sweep.R",
    workspace_file,
    OUT_ROOT,
    "Exploratory cohort-defined p-tau217 enrichment; not diagnostic cutoff; not biomarker-confirmed AD.",
    paste0("CN p", round(MAIN_TEXT_PTAU217_QUANTILE * 100)),
    TRANSFORM_DECISION,
    as.character(Sys.time())
  )
)
safe_write_csv(script_metadata, file.path(LOG_DIR, "script_metadata.csv"))
writeLines(capture.output(utils::sessionInfo()), file.path(LOG_DIR, "sessionInfo.txt"))

message("Exploratory p-tau217-enriched sensitivity analysis complete.")
message("Main-text candidate result saved to: ", file.path(TABLE_DIR, "recommended_main_text_p_tau217_CNp90_result_v3.csv"))
message("Manuscript-ready text saved to: ", file.path(OUT_ROOT, "manuscript_ready_results_p_tau217_CNp90.txt"))
message("Outputs saved to: ", OUT_ROOT)


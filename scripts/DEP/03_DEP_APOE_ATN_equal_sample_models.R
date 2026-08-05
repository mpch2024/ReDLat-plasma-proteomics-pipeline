###############################################################################
# ReDLat plasma proteomics — DEP workflow
# 03. Equal-sample APOE and AT(N) models
# Requires: Outputs from Scripts 01–02
# Produces: Complete-case baseline and adjusted-model comparisons
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
options(stringsAsFactors = FALSE)
set.seed(1234)

required_pkgs <- c(
  "dplyr", "tidyr", "purrr", "tibble", "readr", "stringr", "limma"
)
missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing required packages: ", paste(missing_pkgs, collapse = ", "),
    "\nRun renv::restore() before running this script.",
    call. = FALSE
  )
}
invisible(lapply(required_pkgs, library, character.only = TRUE))


MAIN_GROUPS <- c("CN", "AD")
MAIN_FDR <- 0.05
STRICT_FDR <- 0.01

first_existing_file <- function(paths) {
  paths <- unique(paths[!is.na(paths) & nzchar(paths)])
  hit <- paths[file.exists(paths)][1]
  if (length(hit) == 0 || is.na(hit)) return(NA_character_)
  normalizePath(hit, winslash = "/", mustWork = TRUE)
}

workspace_file <- first_existing_file(c(
  file.path(analysis_root, "workspace", "analysis_workspace.RData")
))
if (is.na(workspace_file)) {
  stop("Could not find analysis_workspace.RData.", call. = FALSE)
}

ws <- new.env(parent = emptyenv())
load(workspace_file, envir = ws)

required_objects <- c("dep_df", "dep_protein_cols", "annot_tbl", "DEP_gene")
missing_objects <- required_objects[
  !vapply(required_objects, exists, logical(1), envir = ws)
]
if (length(missing_objects) > 0) {
  stop(
    "Workspace is missing required objects: ",
    paste(missing_objects, collapse = ", "),
    call. = FALSE
  )
}

analysis_root <- DEP_CONFIG$result_root
APOE_OUT <- file.path(analysis_root, "04_sensitivity", "apoe", "equal_sample_baseline")
ATN_OUT <- file.path(analysis_root, "04_sensitivity", "atn_adjusted", "equal_sample_baseline")
MANIFEST_OUT <- file.path(analysis_root, "07_manifest", "equal_sample_APOE_ATN")
BACKUP_OUT <- file.path(analysis_root, "workspace", "backups_before_equal_sample_models")
invisible(lapply(c(APOE_OUT, ATN_OUT, MANIFEST_OUT, BACKUP_OUT), dir.create,
                 recursive = TRUE, showWarnings = FALSE))

safe_write_csv <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(tibble::as_tibble(x), file)
}

safe_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))

safe_log2_matrix <- function(mat) {
  mat <- apply(mat, 2, safe_numeric)
  mat[mat <= 0] <- NA_real_
  log2(mat)
}

safe_se_from_limma <- function(logFC, t_stat) {
  logFC <- safe_numeric(logFC)
  t_stat <- safe_numeric(t_stat)
  se <- rep(NA_real_, length(logFC))
  ok <- is.finite(logFC) & is.finite(t_stat) & abs(t_stat) > .Machine$double.eps
  se[ok] <- abs(logFC[ok] / t_stat[ok])
  se
}

clean_text_na <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "NaN", "NULL", "null", "N/A")] <- NA_character_
  x
}

make_protein_name <- function(symbol, target_full, target, apt, seqid = NA_character_) {
  dplyr::coalesce(
    clean_text_na(symbol), clean_text_na(target_full), clean_text_na(target),
    clean_text_na(apt), clean_text_na(seqid)
  )
}

classify_dep <- function(tbl, fdr = MAIN_FDR) {
  tbl %>%
    dplyr::mutate(
      type = dplyr::case_when(
        adj.P.Val < fdr & logFC > 0 ~ "Up",
        adj.P.Val < fdr & logFC < 0 ~ "Down",
        TRUE ~ "NS"
      ),
      Direction = dplyr::case_when(
        type == "Up" ~ "Higher in AD",
        type == "Down" ~ "Lower in AD",
        TRUE ~ "Not significant"
      )
    )
}

standardize_annotation <- function(tbl) {
  tbl <- tibble::as_tibble(tbl)
  needed <- c(
    "AptName", "SeqId", "EntrezGeneID", "EntrezGeneSymbol",
    "TargetFullName", "Target", "UniProt", "Protein_Name"
  )
  for (nm in needed) if (!nm %in% names(tbl)) tbl[[nm]] <- NA
  tbl %>%
    dplyr::mutate(AptName = as.character(AptName)) %>%
    dplyr::distinct(AptName, .keep_all = TRUE)
}

build_primary_map <- function(primary_gene) {
  map <- primary_gene %>%
    dplyr::transmute(
      EntrezGeneSymbol = as.character(EntrezGeneSymbol),
      AptName = as.character(AptName)
    ) %>%
    dplyr::filter(
      !is.na(EntrezGeneSymbol), EntrezGeneSymbol != "",
      !is.na(AptName), AptName != ""
    ) %>%
    dplyr::distinct()

  if (anyDuplicated(map$EntrezGeneSymbol) > 0 || anyDuplicated(map$AptName) > 0) {
    stop("The primary map must contain one unique gene and SOMAmer per row.", call. = FALSE)
  }
  map
}

apply_primary_map <- function(dep_aptamer, primary_map, model_name) {
  fixed <- dep_aptamer %>%
    dplyr::mutate(
      EntrezGeneSymbol = as.character(EntrezGeneSymbol),
      AptName = as.character(AptName)
    ) %>%
    dplyr::inner_join(primary_map, by = c("EntrezGeneSymbol", "AptName")) %>%
    dplyr::arrange(match(EntrezGeneSymbol, primary_map$EntrezGeneSymbol))

  if (nrow(fixed) != nrow(primary_map)) {
    missing <- dplyr::anti_join(
      primary_map,
      fixed %>% dplyr::select(EntrezGeneSymbol, AptName),
      by = c("EntrezGeneSymbol", "AptName")
    )
    stop(
      model_name, " retained ", nrow(fixed), " of ", nrow(primary_map),
      " fixed SOMAmers. Missing examples: ",
      paste(utils::head(paste0(missing$EntrezGeneSymbol, "=", missing$AptName), 10),
            collapse = "; "),
      call. = FALSE
    )
  }
  if (anyDuplicated(fixed$EntrezGeneSymbol) > 0 || anyDuplicated(fixed$AptName) > 0) {
    stop(model_name, " contains duplicated fixed-map entries.", call. = FALSE)
  }
  fixed
}

prepare_core_covariates <- function(dat) {
  dat %>%
    dplyr::mutate(
      SampleGroup = factor(as.character(SampleGroup), levels = MAIN_GROUPS),
      Sex = factor(as.character(Sex)),
      Country = factor(as.character(Country)),
      Age = safe_numeric(Age),
      Education = safe_numeric(Education)
    ) %>%
    dplyr::filter(SampleGroup %in% MAIN_GROUPS)
}

run_limma_model <- function(dat, formula_str, coef_name, model_name,
                            dep_protein_cols, annot_tbl) {
  dat <- prepare_core_covariates(dat)
  protein_cols <- intersect(dep_protein_cols, names(dat))
  if (length(protein_cols) == 0) stop("No protein columns found for ", model_name)

  metadata <- dat %>% dplyr::select(-dplyr::all_of(protein_cols))
  model_vars <- all.vars(stats::as.formula(formula_str))
  missing_vars <- setdiff(model_vars, names(metadata))
  if (length(missing_vars) > 0) {
    stop(model_name, " is missing variables: ", paste(missing_vars, collapse = ", "))
  }

  keep <- stats::complete.cases(metadata[, model_vars, drop = FALSE])
  dat <- dat[keep, , drop = FALSE]
  metadata <- metadata[keep, , drop = FALSE]

  expr <- dat %>%
    dplyr::select(dplyr::all_of(protein_cols)) %>%
    as.matrix() %>%
    safe_log2_matrix() %>%
    t()

  design <- stats::model.matrix(stats::as.formula(formula_str), data = metadata)
  if (qr(design)$rank != ncol(design)) {
    stop(model_name, " has a rank-deficient design matrix.", call. = FALSE)
  }
  if (!coef_name %in% colnames(design)) {
    stop(model_name, " coefficient not found: ", coef_name, call. = FALSE)
  }

  fit <- limma::eBayes(limma::lmFit(expr, design))
  dep <- limma::topTable(
    fit, coef = coef_name, adjust.method = "BH", number = Inf, sort.by = "P"
  ) %>%
    tibble::rownames_to_column("feature_id_raw") %>%
    dplyr::mutate(
      AptName = as.character(feature_id_raw),
      se = safe_se_from_limma(logFC, t)
    ) %>%
    dplyr::left_join(annot_tbl, by = "AptName") %>%
    dplyr::mutate(
      Protein_Name = make_protein_name(
        EntrezGeneSymbol, TargetFullName, Target, AptName, feature_id_raw
      ),
      model_name = model_name,
      formula = formula_str
    ) %>%
    classify_dep()

  list(
    dep = dep,
    metadata = metadata,
    design = design,
    sample_ids = if ("SampleId" %in% names(metadata)) as.character(metadata$SampleId) else character(0)
  )
}

export_dep <- function(tbl, file) {
  keep <- intersect(c(
    "Protein_Name", "EntrezGeneSymbol", "EntrezGeneID", "TargetFullName",
    "Target", "UniProt", "AptName", "SeqId", "feature_id_raw",
    "logFC", "se", "AveExpr", "t", "P.Value", "adj.P.Val", "B",
    "type", "Direction", "model_name", "formula"
  ), names(tbl))
  tbl %>%
    dplyr::select(dplyr::all_of(keep)) %>%
    dplyr::arrange(adj.P.Val, dplyr::desc(abs(logFC)), AptName) %>%
    safe_write_csv(file)
}

count_dep <- function(tbl, label) {
  purrr::map_dfr(c(STRICT_FDR, MAIN_FDR), function(fdr) {
    tibble::tibble(
      model = label,
      fdr = fdr,
      total = sum(tbl$adj.P.Val < fdr, na.rm = TRUE),
      higher_in_AD = sum(tbl$adj.P.Val < fdr & tbl$logFC > 0, na.rm = TRUE),
      lower_in_AD = sum(tbl$adj.P.Val < fdr & tbl$logFC < 0, na.rm = TRUE)
    )
  })
}

build_three_way_comparison <- function(primary_full, subset_baseline, adjusted,
                                       analysis_label) {
  primary_tbl <- primary_full %>%
    dplyr::select(
      EntrezGeneSymbol, AptName, Protein_Name,
      primary_full_logFC = logFC,
      primary_full_P.Value = P.Value,
      primary_full_adj.P.Val = adj.P.Val
    )

  baseline_tbl <- subset_baseline %>%
    dplyr::select(
      EntrezGeneSymbol, AptName,
      subset_baseline_logFC = logFC,
      subset_baseline_P.Value = P.Value,
      subset_baseline_adj.P.Val = adj.P.Val
    )

  adjusted_tbl <- adjusted %>%
    dplyr::select(
      EntrezGeneSymbol, AptName,
      adjusted_logFC = logFC,
      adjusted_P.Value = P.Value,
      adjusted_adj.P.Val = adj.P.Val
    )

  primary_tbl %>%
    dplyr::inner_join(baseline_tbl, by = c("EntrezGeneSymbol", "AptName")) %>%
    dplyr::inner_join(adjusted_tbl, by = c("EntrezGeneSymbol", "AptName")) %>%
    dplyr::mutate(
      analysis = analysis_label,
      same_fixed_somamer = TRUE,
      same_direction_full_vs_subset =
        sign(primary_full_logFC) == sign(subset_baseline_logFC),
      same_direction_subset_vs_adjusted =
        sign(subset_baseline_logFC) == sign(adjusted_logFC),
      primary_full_fdr005 = primary_full_adj.P.Val < MAIN_FDR,
      subset_baseline_fdr005 = subset_baseline_adj.P.Val < MAIN_FDR,
      adjusted_fdr005 = adjusted_adj.P.Val < MAIN_FDR,
      primary_preserved_after_subset_restriction =
        primary_full_fdr005 & subset_baseline_fdr005 & same_direction_full_vs_subset,
      subset_baseline_preserved_after_adjustment =
        subset_baseline_fdr005 & adjusted_fdr005 & same_direction_subset_vs_adjusted,
      sample_restriction_delta_logFC = subset_baseline_logFC - primary_full_logFC,
      covariate_adjustment_delta_logFC = adjusted_logFC - subset_baseline_logFC,
      covariate_adjustment_abs_ratio =
        abs(adjusted_logFC) / pmax(abs(subset_baseline_logFC), 1e-12)
    )
}

summarize_three_way <- function(compare_tbl, n_total, n_cn, n_ad,
                                baseline_formula, adjusted_formula,
                                covariates_added) {
  baseline_sig <- compare_tbl$subset_baseline_fdr005
  primary_sig <- compare_tbl$primary_full_fdr005

  tibble::tibble(
    analysis = unique(compare_tbl$analysis),
    n_complete_case = n_total,
    n_CN = n_cn,
    n_AD = n_ad,
    n_fixed_primary_somamers = nrow(compare_tbl),
    baseline_formula = baseline_formula,
    adjusted_formula = adjusted_formula,
    covariates_added = paste(covariates_added, collapse = ", "),
    primary_full_fdr005 = sum(primary_sig, na.rm = TRUE),
    subset_baseline_fdr005 = sum(baseline_sig, na.rm = TRUE),
    adjusted_fdr005 = sum(compare_tbl$adjusted_fdr005, na.rm = TRUE),
    correlation_full_vs_subset_baseline = suppressWarnings(stats::cor(
      compare_tbl$primary_full_logFC,
      compare_tbl$subset_baseline_logFC,
      use = "complete.obs"
    )),
    correlation_subset_baseline_vs_adjusted = suppressWarnings(stats::cor(
      compare_tbl$subset_baseline_logFC,
      compare_tbl$adjusted_logFC,
      use = "complete.obs"
    )),
    direction_consistency_full_vs_subset = mean(
      compare_tbl$same_direction_full_vs_subset, na.rm = TRUE
    ),
    direction_consistency_subset_vs_adjusted = mean(
      compare_tbl$same_direction_subset_vs_adjusted, na.rm = TRUE
    ),
    proportion_primary_DEPs_preserved_after_subset_restriction = if (any(primary_sig, na.rm = TRUE)) {
      mean(compare_tbl$primary_preserved_after_subset_restriction[primary_sig], na.rm = TRUE)
    } else NA_real_,
    proportion_subset_baseline_DEPs_preserved_after_adjustment = if (any(baseline_sig, na.rm = TRUE)) {
      mean(compare_tbl$subset_baseline_preserved_after_adjustment[baseline_sig], na.rm = TRUE)
    } else NA_real_
  )
}

normalize_apoe4 <- function(x) {
  if (is.numeric(x) || is.integer(x)) {
    out <- safe_numeric(x)
    out[!out %in% c(0, 1)] <- NA_real_
    return(out)
  }
  xx <- tolower(trimws(as.character(x)))
  dplyr::case_when(
    xx %in% c("1", "yes", "y", "true", "carrier", "e4 carrier", "ε4 carrier") ~ 1,
    xx %in% c("0", "no", "n", "false", "non-carrier", "noncarrier", "non-e4") ~ 0,
    TRUE ~ NA_real_
  )
}

create_apoe4_if_needed <- function(dat) {
  if ("APOE4_carrier" %in% names(dat)) {
    dat$APOE4_carrier <- normalize_apoe4(dat$APOE4_carrier)
    return(dat)
  }
  genotype_col <- c("ApoE", "APOE", "apoe_genotype", "APOE_genotype")
  genotype_col <- genotype_col[genotype_col %in% names(dat)][1]
  if (length(genotype_col) == 0 || is.na(genotype_col)) {
    stop("No APOE4_carrier or APOE genotype column was found.", call. = FALSE)
  }
  geno <- tolower(trimws(as.character(dat[[genotype_col]])))
  dat$APOE4_carrier <- dplyr::case_when(
    geno %in% c("e2/e4", "e3/e4", "e4/e4", "2/4", "3/4", "4/4") ~ 1,
    geno %in% c("e2/e2", "e2/e3", "e3/e3", "2/2", "2/3", "3/3") ~ 0,
    TRUE ~ NA_real_
  )
  dat
}

pick_first_existing <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)][1]
  if (length(hit) == 0 || is.na(hit)) NA_character_ else hit
}

detect_atn_covariates <- function(df) {
  detected <- c(
    pick_first_existing(df, c("p.tau217", "p_tau217", "p-tau217", "ptau217")),
    pick_first_existing(df, c("NfL", "NFL", "nfl", "neurofilament_light")),
    pick_first_existing(df, c(
      "ratio.AB42.40", "ratio_AB42_40", "ratio AB42/40",
      "AB42_40_ratio", "Aβ42/40"
    ))
  )
  unique(detected[!is.na(detected)])
}

primary_gene <- tibble::as_tibble(ws$DEP_gene) %>%
  dplyr::mutate(
    EntrezGeneSymbol = as.character(EntrezGeneSymbol),
    AptName = as.character(AptName)
  )
primary_map <- if (exists("primary_gene_map_fixed", envir = ws)) {
  tibble::as_tibble(ws$primary_gene_map_fixed)
} else {
  build_primary_map(primary_gene)
}
primary_map <- build_primary_map(primary_map)

if (nrow(primary_gene) != nrow(primary_map)) {
  stop("DEP_gene and the fixed primary map do not have the same number of rows.")
}

annot_tbl <- standardize_annotation(ws$annot_tbl)
dep_df <- tibble::as_tibble(ws$dep_df)
dep_protein_cols <- as.character(ws$dep_protein_cols)

safe_write_csv(primary_map, file.path(MANIFEST_OUT, "fixed_primary_SOMAmer_gene_map.csv"))

# =============================================================================
# APOE: equal-sample baseline versus APOE-adjusted model
# =============================================================================

apoe_df <- dep_df %>%
  create_apoe4_if_needed() %>%
  prepare_core_covariates() %>%
  dplyr::mutate(APOE4_carrier = factor(APOE4_carrier, levels = c(0, 1))) %>%
  dplyr::filter(stats::complete.cases(
    SampleGroup, Age, Sex, Country, Education, APOE4_carrier
  ))

if (nrow(apoe_df) < 20 || dplyr::n_distinct(apoe_df$APOE4_carrier) < 2) {
  stop("Insufficient complete APOE data for equal-sample models.", call. = FALSE)
}

apoe_baseline_formula <- "~ SampleGroup + Age + Sex + Country + Education"
apoe_adjusted_formula <- paste0(apoe_baseline_formula, " + APOE4_carrier")

apoe_baseline_fit <- run_limma_model(
  apoe_df, apoe_baseline_formula, "SampleGroupAD", "APOE_complete_case_baseline",
  dep_protein_cols, annot_tbl
)
apoe_adjusted_fit <- run_limma_model(
  apoe_df, apoe_adjusted_formula, "SampleGroupAD", "APOE_equal_sample_adjusted",
  dep_protein_cols, annot_tbl
)

apoe_baseline_gene <- apply_primary_map(
  apoe_baseline_fit$dep, primary_map, "APOE_complete_case_baseline"
)
apoe_adjusted_gene <- apply_primary_map(
  apoe_adjusted_fit$dep, primary_map, "APOE_equal_sample_adjusted"
)

apoe_compare <- build_three_way_comparison(
  primary_gene, apoe_baseline_gene, apoe_adjusted_gene,
  "APOE_equal_sample_baseline_vs_adjusted"
)
apoe_summary <- summarize_three_way(
  apoe_compare,
  nrow(apoe_adjusted_fit$metadata),
  sum(apoe_adjusted_fit$metadata$SampleGroup == "CN"),
  sum(apoe_adjusted_fit$metadata$SampleGroup == "AD"),
  apoe_baseline_formula,
  apoe_adjusted_formula,
  "APOE4_carrier"
)

export_dep(apoe_baseline_fit$dep, file.path(
  APOE_OUT, "APOE_complete_case_baseline_full_limma_results_aptamer_level.csv"
))
export_dep(apoe_baseline_gene, file.path(
  APOE_OUT, "APOE_complete_case_baseline_full_limma_results_FIXED_PRIMARY_MAP.csv"
))
export_dep(apoe_adjusted_gene, file.path(
  APOE_OUT, "APOE_equal_sample_adjusted_full_limma_results_FIXED_PRIMARY_MAP.csv"
))
safe_write_csv(apoe_compare, file.path(
  APOE_OUT, "APOE_primary_full_vs_subset_baseline_vs_adjusted_FIXED_PRIMARY_MAP.csv"
))
safe_write_csv(apoe_summary, file.path(
  APOE_OUT, "APOE_equal_sample_model_summary_FIXED_PRIMARY_MAP.csv"
))
safe_write_csv(
  dplyr::bind_rows(
    count_dep(apoe_baseline_gene, "APOE_complete_case_baseline"),
    count_dep(apoe_adjusted_gene, "APOE_equal_sample_adjusted")
  ),
  file.path(APOE_OUT, "APOE_equal_sample_DEP_counts_FIXED_PRIMARY_MAP.csv")
)
safe_write_csv(
  tibble::tibble(SampleId = apoe_adjusted_fit$sample_ids),
  file.path(analysis_root, "private", "participant_sets", "APOE_equal_sample_participant_ids.csv")
)

# =============================================================================
# AT(N): equal-sample baseline versus AT(N)-adjusted model
# =============================================================================

atn_covariates <- detect_atn_covariates(dep_df)
if (length(atn_covariates) < 2) {
  stop(
    "Fewer than two AT(N) covariates were detected: ",
    paste(atn_covariates, collapse = ", "),
    call. = FALSE
  )
}

atn_df <- dep_df %>%
  prepare_core_covariates() %>%
  dplyr::mutate(dplyr::across(
    dplyr::all_of(atn_covariates), safe_numeric
  )) %>%
  dplyr::filter(dplyr::if_all(
    dplyr::all_of(c(
      "SampleGroup", "Age", "Sex", "Country", "Education", atn_covariates
    )),
    ~ !is.na(.x)
  ))

if (nrow(atn_df) < 20) {
  stop("Insufficient complete AT(N) data for equal-sample models.", call. = FALSE)
}

atn_baseline_formula <- "~ SampleGroup + Age + Sex + Country + Education"
atn_adjusted_formula <- paste0(
  atn_baseline_formula, " + ", paste(atn_covariates, collapse = " + ")
)

atn_baseline_fit <- run_limma_model(
  atn_df, atn_baseline_formula, "SampleGroupAD", "ATN_complete_case_baseline",
  dep_protein_cols, annot_tbl
)
atn_adjusted_fit <- run_limma_model(
  atn_df, atn_adjusted_formula, "SampleGroupAD", "ATN_equal_sample_adjusted",
  dep_protein_cols, annot_tbl
)

atn_baseline_gene <- apply_primary_map(
  atn_baseline_fit$dep, primary_map, "ATN_complete_case_baseline"
)
atn_adjusted_gene <- apply_primary_map(
  atn_adjusted_fit$dep, primary_map, "ATN_equal_sample_adjusted"
)

atn_compare <- build_three_way_comparison(
  primary_gene, atn_baseline_gene, atn_adjusted_gene,
  "ATN_equal_sample_baseline_vs_adjusted"
)
atn_summary <- summarize_three_way(
  atn_compare,
  nrow(atn_adjusted_fit$metadata),
  sum(atn_adjusted_fit$metadata$SampleGroup == "CN"),
  sum(atn_adjusted_fit$metadata$SampleGroup == "AD"),
  atn_baseline_formula,
  atn_adjusted_formula,
  atn_covariates
)

export_dep(atn_baseline_fit$dep, file.path(
  ATN_OUT, "ATN_complete_case_baseline_full_limma_results_aptamer_level.csv"
))
export_dep(atn_baseline_gene, file.path(
  ATN_OUT, "ATN_complete_case_baseline_full_limma_results_FIXED_PRIMARY_MAP.csv"
))
export_dep(atn_adjusted_gene, file.path(
  ATN_OUT, "ATN_equal_sample_adjusted_full_limma_results_FIXED_PRIMARY_MAP.csv"
))
safe_write_csv(atn_compare, file.path(
  ATN_OUT, "ATN_primary_full_vs_subset_baseline_vs_adjusted_FIXED_PRIMARY_MAP.csv"
))
safe_write_csv(atn_summary, file.path(
  ATN_OUT, "ATN_equal_sample_model_summary_FIXED_PRIMARY_MAP.csv"
))
safe_write_csv(
  dplyr::bind_rows(
    count_dep(atn_baseline_gene, "ATN_complete_case_baseline"),
    count_dep(atn_adjusted_gene, "ATN_equal_sample_adjusted")
  ),
  file.path(ATN_OUT, "ATN_equal_sample_DEP_counts_FIXED_PRIMARY_MAP.csv")
)
safe_write_csv(
  tibble::tibble(covariate = atn_covariates),
  file.path(ATN_OUT, "ATN_equal_sample_covariates_used.csv")
)
safe_write_csv(
  tibble::tibble(SampleId = atn_adjusted_fit$sample_ids),
  file.path(analysis_root, "private", "participant_sets", "ATN_equal_sample_participant_ids.csv")
)

# =============================================================================
# Workspace update and audit
# =============================================================================

run_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
backup_file <- file.path(
  BACKUP_OUT,
  paste0("analysis_workspace_before_equal_sample_", run_stamp, ".RData")
)
file.copy(workspace_file, backup_file, overwrite = FALSE)

ws$primary_gene_map_fixed <- primary_map
ws$DEP_APOE_complete_case_baseline_gene <- apoe_baseline_gene
ws$DEP_APOE_equal_sample_adjusted_gene <- apoe_adjusted_gene
ws$apoe_equal_sample_compare_tbl <- apoe_compare
ws$apoe_equal_sample_summary <- apoe_summary
ws$DEP_ATN_complete_case_baseline_gene <- atn_baseline_gene
ws$DEP_ATN_equal_sample_adjusted_gene <- atn_adjusted_gene
ws$atn_equal_sample_compare_tbl <- atn_compare
ws$atn_equal_sample_summary <- atn_summary
ws$atn_equal_sample_covariates <- atn_covariates
save(list = ls(ws), file = workspace_file, envir = ws)

run_audit <- tibble::tibble(
  item = c(
    "run_time", "workspace_file", "workspace_backup",
    "fixed_map_rows", "APOE_complete_case_n", "APOE_CN", "APOE_AD",
    "ATN_complete_case_n", "ATN_CN", "ATN_AD", "ATN_covariates",
    "primary_DEP_modified"
  ),
  value = c(
    as.character(Sys.time()), workspace_file, backup_file,
    as.character(nrow(primary_map)),
    as.character(nrow(apoe_adjusted_fit$metadata)),
    as.character(sum(apoe_adjusted_fit$metadata$SampleGroup == "CN")),
    as.character(sum(apoe_adjusted_fit$metadata$SampleGroup == "AD")),
    as.character(nrow(atn_adjusted_fit$metadata)),
    as.character(sum(atn_adjusted_fit$metadata$SampleGroup == "CN")),
    as.character(sum(atn_adjusted_fit$metadata$SampleGroup == "AD")),
    paste(atn_covariates, collapse = ", "),
    "FALSE"
  )
)
safe_write_csv(run_audit, file.path(MANIFEST_OUT, "equal_sample_APOE_ATN_run_audit.csv"))
safe_write_csv(
  dplyr::bind_rows(apoe_summary, atn_summary),
  file.path(MANIFEST_OUT, "equal_sample_APOE_ATN_summary.csv")
)
writeLines(
  capture.output(utils::sessionInfo()),
  file.path(MANIFEST_OUT, "equal_sample_APOE_ATN_sessionInfo.txt")
)

message("Equal-sample APOE and AT(N) models completed.")
message("Primary DEP was not modified.")
message("APOE outputs: ", APOE_OUT)
message("AT(N) outputs: ", ATN_OUT)
message("Updated workspace: ", workspace_file)
###############################################################################
# END
###############################################################################


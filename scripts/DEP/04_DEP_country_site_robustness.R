###############################################################################
# ReDLat plasma proteomics — DEP workflow
# 04. Country and site robustness analyses
# Requires: Outputs from Scripts 01–02
# Produces: LOCO, LOSO, resampling, interaction and meta-analysis results
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
  "dplyr", "tidyr", "purrr", "tibble", "readr", "stringr",
  "rlang", "limma", "metafor"
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
MIN_N_PER_GROUP_COUNTRY <- 10L
MIN_N_PER_SITE_GROUP <- 8L
BALANCED_RESAMPLING_NITER <- 200L
BALANCED_RESAMPLING_MIN_GROUP <- 8L
SEED_GLOBAL <- 1234L
UPDATE_WORKSPACE <- TRUE
set.seed(SEED_GLOBAL)

first_existing_file <- function(paths) {
  paths <- unique(paths[!is.na(paths) & nzchar(paths)])
  hit <- paths[file.exists(paths)][1]
  if (length(hit) == 0 || is.na(hit)) return(NA_character_)
  normalizePath(hit, winslash = "/", mustWork = TRUE)
}

workspace_file <- first_existing_file(c(
  file.path(analysis_root, "workspace", "analysis_workspace.RData")
))
if (is.na(workspace_file)) stop("Could not find analysis_workspace.RData.", call. = FALSE)

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
LOCO_DIR <- file.path(analysis_root, "06_robustness", "country_loco", "tables")
META_DIR <- file.path(analysis_root, "06_robustness", "country_meta", "tables")
INTERACTION_DIR <- file.path(analysis_root, "06_robustness", "country_interaction", "tables")
LOSO_DIR <- file.path(analysis_root, "06_robustness", "site_robustness", "loso", "tables")
BALANCED_DIR <- file.path(analysis_root, "06_robustness", "balanced_country_resampling", "tables")
CLASS_DIR <- file.path(analysis_root, "06_robustness", "formal_classification")
MANIFEST_DIR <- file.path(analysis_root, "07_manifest", "fixed_primary_map_robustness_v2")
BACKUP_DIR <- file.path(analysis_root, "workspace", "backups_before_fixed_robustness_v2")
invisible(lapply(
  c(LOCO_DIR, META_DIR, INTERACTION_DIR, LOSO_DIR, BALANCED_DIR,
    CLASS_DIR, MANIFEST_DIR, BACKUP_DIR),
  dir.create, recursive = TRUE, showWarnings = FALSE
))

safe_write_csv <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(tibble::as_tibble(x), file)
}

safe_file_tag <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
  gsub("^_+|_+$", "", x)
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

make_protein_name <- function(symbol, target_full, target, apt, raw = NA_character_) {
  dplyr::coalesce(
    clean_text_na(symbol), clean_text_na(target_full), clean_text_na(target),
    clean_text_na(apt), clean_text_na(raw)
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
    stop("The fixed primary map contains duplicated genes or SOMAmers.", call. = FALSE)
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
      " fixed pairs. Missing examples: ",
      paste(utils::head(paste0(missing$EntrezGeneSymbol, "=", missing$AptName), 10),
            collapse = "; "),
      call. = FALSE
    )
  }
  if (anyDuplicated(fixed$EntrezGeneSymbol) > 0 || anyDuplicated(fixed$AptName) > 0) {
    stop(model_name, " contains duplicate fixed-map rows.", call. = FALSE)
  }
  fixed
}

prepare_model_data <- function(dat) {
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

run_limma_dep <- function(dat, formula_str, coef_name, model_name,
                          dep_protein_cols, annot_tbl) {
  dat <- prepare_model_data(dat)
  protein_cols <- intersect(dep_protein_cols, names(dat))
  if (length(protein_cols) == 0) stop("No protein columns found for ", model_name)

  metadata <- dat %>% dplyr::select(-dplyr::all_of(protein_cols))
  model_vars <- all.vars(stats::as.formula(formula_str))
  missing_vars <- setdiff(model_vars, names(metadata))
  if (length(missing_vars) > 0) {
    stop(model_name, " is missing model variables: ", paste(missing_vars, collapse = ", "))
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
      model_formula = formula_str
    ) %>%
    classify_dep()

  list(dep = dep, metadata = metadata, design = design, fit = fit)
}

export_dep <- function(tbl, file) {
  keep <- intersect(c(
    "Protein_Name", "EntrezGeneSymbol", "EntrezGeneID", "TargetFullName",
    "Target", "UniProt", "AptName", "SeqId", "feature_id_raw",
    "logFC", "se", "AveExpr", "t", "P.Value", "adj.P.Val", "B",
    "type", "Direction", "model_name", "model_formula",
    "excluded_country", "excluded_site", "Country"
  ), names(tbl))
  tbl %>%
    dplyr::select(dplyr::all_of(keep)) %>%
    dplyr::arrange(adj.P.Val, dplyr::desc(abs(logFC)), AptName) %>%
    safe_write_csv(file)
}

fixed_comparison <- function(main_gene, secondary_gene, context_name, context_value) {
  context_col <- context_name
  main_gene %>%
    dplyr::select(
      EntrezGeneSymbol, AptName, Protein_Name,
      main_logFC = logFC, main_P.Value = P.Value,
      main_adj.P.Val = adj.P.Val, main_type = type
    ) %>%
    dplyr::inner_join(
      secondary_gene %>%
        dplyr::select(
          EntrezGeneSymbol, AptName,
          context_logFC = logFC, context_P.Value = P.Value,
          context_adj.P.Val = adj.P.Val, context_type = type
        ),
      by = c("EntrezGeneSymbol", "AptName")
    ) %>%
    dplyr::mutate(
      !!context_col := as.character(context_value),
      same_fixed_somamer = TRUE,
      same_direction = sign(main_logFC) == sign(context_logFC),
      main_sig_fdr005 = main_adj.P.Val < MAIN_FDR,
      context_sig_fdr005 = context_adj.P.Val < MAIN_FDR,
      main_sig_preserved = main_sig_fdr005 & context_sig_fdr005 & same_direction,
      delta_logFC = context_logFC - main_logFC
    )
}

summarize_comparison <- function(cmp, n_samples, context_name, context_value,
                                 n_fixed_map) {
  main_sig <- cmp$main_sig_fdr005
  slope <- tryCatch(
    unname(stats::coef(stats::lm(cmp$context_logFC ~ cmp$main_logFC))[2]),
    error = function(e) NA_real_
  )
  tibble::tibble(
    !!context_name := as.character(context_value),
    n_samples = n_samples,
    n_fixed_primary_somamers = n_fixed_map,
    logFC_correlation = suppressWarnings(stats::cor(
      cmp$main_logFC, cmp$context_logFC, use = "complete.obs"
    )),
    slope = slope,
    direction_consistency_all = mean(cmp$same_direction, na.rm = TRUE),
    n_main_sig_fdr005 = sum(main_sig, na.rm = TRUE),
    prop_main_sig_preserved = if (any(main_sig, na.rm = TRUE)) {
      mean(cmp$main_sig_preserved[main_sig], na.rm = TRUE)
    } else NA_real_
  )
}

ensure_count_cols <- function(tbl, cols) {
  for (nm in cols) if (!nm %in% names(tbl)) tbl[[nm]] <- 0L
  tbl
}

primary_gene <- tibble::as_tibble(ws$DEP_gene) %>%
  dplyr::mutate(
    EntrezGeneSymbol = as.character(EntrezGeneSymbol),
    AptName = as.character(AptName)
  )
primary_map <- if (exists("primary_gene_map_fixed", envir = ws)) {
  build_primary_map(tibble::as_tibble(ws$primary_gene_map_fixed))
} else {
  build_primary_map(primary_gene)
}
if (nrow(primary_gene) != nrow(primary_map)) {
  stop("DEP_gene and the fixed map differ in size.", call. = FALSE)
}

main_gene <- apply_primary_map(primary_gene, primary_map, "primary_DEP_gene")
annot_tbl <- standardize_annotation(ws$annot_tbl)
dep_df <- prepare_model_data(tibble::as_tibble(ws$dep_df))
dep_protein_cols <- as.character(ws$dep_protein_cols)

safe_write_csv(primary_map, file.path(MANIFEST_DIR, "fixed_primary_SOMAmer_gene_map.csv"))
map_md5 <- unname(tools::md5sum(file.path(MANIFEST_DIR, "fixed_primary_SOMAmer_gene_map.csv")))

# =============================================================================
# Country counts and eligible countries
# =============================================================================

country_counts <- dep_df %>%
  dplyr::count(Country, SampleGroup, name = "n") %>%
  tidyr::pivot_wider(names_from = SampleGroup, values_from = n, values_fill = 0) %>%
  ensure_count_cols(c("CN", "AD")) %>%
  dplyr::mutate(total = CN + AD) %>%
  dplyr::arrange(dplyr::desc(total))

countries_to_test <- country_counts %>%
  dplyr::filter(CN >= MIN_N_PER_GROUP_COUNTRY, AD >= MIN_N_PER_GROUP_COUNTRY) %>%
  dplyr::pull(Country) %>%
  as.character()

if (length(countries_to_test) < 2) {
  stop("Fewer than two countries have adequate CN and AD sample sizes.")
}
safe_write_csv(country_counts, file.path(LOCO_DIR, "country_group_counts.csv"))

# =============================================================================
# LOCO with fixed primary map
# =============================================================================

loco_summary_list <- list()
loco_compare_list <- list()
loco_gene_list <- list()

for (country_i in countries_to_test) {
  message("LOCO fixed-map model excluding: ", country_i)
  dat_loco <- dep_df %>%
    dplyr::filter(as.character(Country) != country_i) %>%
    dplyr::mutate(Country = droplevels(factor(Country)))

  fit_loco <- run_limma_dep(
    dat_loco,
    "~ SampleGroup + Age + Sex + Country + Education",
    "SampleGroupAD",
    paste0("LOCO_excluding_", safe_file_tag(country_i)),
    dep_protein_cols,
    annot_tbl
  )
  gene_loco <- apply_primary_map(
    fit_loco$dep, primary_map, paste0("LOCO_excluding_", country_i)
  ) %>%
    dplyr::mutate(excluded_country = country_i)

  cmp <- fixed_comparison(
    main_gene, gene_loco, "excluded_country", country_i
  )

  export_dep(gene_loco, file.path(
    LOCO_DIR,
    paste0("LOCO_excluding_", safe_file_tag(country_i),
           "_dep_results_FIXED_PRIMARY_MAP.csv")
  ))
  safe_write_csv(cmp, file.path(
    LOCO_DIR,
    paste0("LOCO_excluding_", safe_file_tag(country_i),
           "_comparison_to_main_FIXED_PRIMARY_MAP.csv")
  ))

  loco_summary_list[[country_i]] <- summarize_comparison(
    cmp, nrow(fit_loco$metadata), "excluded_country", country_i,
    nrow(primary_map)
  )
  loco_compare_list[[country_i]] <- cmp
  loco_gene_list[[country_i]] <- gene_loco %>%
    dplyr::select(
      EntrezGeneSymbol, AptName, excluded_country,
      loco_logFC = logFC, loco_P.Value = P.Value,
      loco_adj.P.Val = adj.P.Val
    )
}

loco_summary <- dplyr::bind_rows(loco_summary_list)
loco_all <- dplyr::bind_rows(loco_gene_list)

main_vs_mean_loco <- main_gene %>%
  dplyr::select(
    Protein_Name, EntrezGeneSymbol, AptName,
    main_logFC = logFC, main_adj.P.Val = adj.P.Val,
    main_type = type, Direction
  ) %>%
  dplyr::left_join(
    loco_all %>%
      dplyr::group_by(EntrezGeneSymbol, AptName) %>%
      dplyr::summarise(
        mean_loco_logFC = mean(loco_logFC, na.rm = TRUE),
        sd_loco_logFC = stats::sd(loco_logFC, na.rm = TRUE),
        min_loco_adj.P.Val = min(loco_adj.P.Val, na.rm = TRUE),
        max_loco_adj.P.Val = max(loco_adj.P.Val, na.rm = TRUE),
        prop_same_direction = mean(
          sign(loco_logFC) == sign(main_gene$logFC[
            match(AptName, main_gene$AptName)
          ]),
          na.rm = TRUE
        ),
        prop_sig_fdr005 = mean(loco_adj.P.Val < MAIN_FDR, na.rm = TRUE),
        n_loco_models = dplyr::n(),
        .groups = "drop"
      ),
    by = c("EntrezGeneSymbol", "AptName")
  ) %>%
  dplyr::mutate(
    Protein_Label = Protein_Name,
    loco_delta = mean_loco_logFC - main_logFC,
    all_loco_same_direction = prop_same_direction == 1,
    all_loco_fdr005 = prop_sig_fdr005 == 1,
    all_loco_fdr005_same_direction = all_loco_same_direction & all_loco_fdr005
  ) %>%
  dplyr::select(
    Protein_Label, Protein_Name, EntrezGeneSymbol, AptName,
    main_logFC, main_adj.P.Val, mean_loco_logFC, sd_loco_logFC,
    loco_delta, prop_same_direction, prop_sig_fdr005,
    all_loco_same_direction, all_loco_fdr005,
    all_loco_fdr005_same_direction, n_loco_models, Direction
  )

safe_write_csv(loco_summary, file.path(
  LOCO_DIR, "LOCO_summary_metrics_FIXED_PRIMARY_MAP.csv"
))
safe_write_csv(main_vs_mean_loco, file.path(
  LOCO_DIR, "main_vs_meanLOCO_table_FIXED_PRIMARY_MAP.csv"
))

# =============================================================================
# Country-specific fixed-map models and random-effects meta-analysis
# =============================================================================

country_specific_list <- list()
for (country_i in countries_to_test) {
  message("Country-specific fixed-map model: ", country_i)
  dat_country <- dep_df %>%
    dplyr::filter(as.character(Country) == country_i)

  fit_country <- run_limma_dep(
    dat_country,
    "~ SampleGroup + Age + Sex + Education",
    "SampleGroupAD",
    paste0("country_specific_", safe_file_tag(country_i)),
    dep_protein_cols,
    annot_tbl
  )
  gene_country <- apply_primary_map(
    fit_country$dep, primary_map, paste0("country_specific_", country_i)
  ) %>%
    dplyr::mutate(Country = country_i)

  country_specific_list[[country_i]] <- gene_country %>%
    dplyr::select(
      Country, AptName, Protein_Name, EntrezGeneSymbol,
      logFC, se, P.Value, adj.P.Val, type
    )
}
country_specific_fixed <- dplyr::bind_rows(country_specific_list)

meta_tbl <- country_specific_fixed %>%
  dplyr::group_by(AptName, Protein_Name, EntrezGeneSymbol) %>%
  dplyr::group_modify(function(.x, .y) {
    dd <- .x %>%
      dplyr::filter(is.finite(logFC), is.finite(se), se > 0)
    if (nrow(dd) < 2) {
      return(tibble::tibble(
        meta_logFC = NA_real_, meta_se = NA_real_, meta_p = NA_real_,
        I2 = NA_real_, tau2 = NA_real_, k = nrow(dd)
      ))
    }
    fit <- tryCatch(
      metafor::rma.uni(yi = dd$logFC, sei = dd$se, method = "REML"),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(tibble::tibble(
        meta_logFC = NA_real_, meta_se = NA_real_, meta_p = NA_real_,
        I2 = NA_real_, tau2 = NA_real_, k = nrow(dd)
      ))
    }
    tibble::tibble(
      meta_logFC = as.numeric(fit$b[1]),
      meta_se = as.numeric(fit$se[1]),
      meta_p = as.numeric(fit$pval[1]),
      I2 = as.numeric(fit$I2) / 100,
      tau2 = as.numeric(fit$tau2),
      k = nrow(dd)
    )
  }) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    meta_adj.P.Val = stats::p.adjust(meta_p, method = "BH"),
    meta_type = dplyr::case_when(
      meta_adj.P.Val < MAIN_FDR & meta_logFC > 0 ~ "Up",
      meta_adj.P.Val < MAIN_FDR & meta_logFC < 0 ~ "Down",
      TRUE ~ "NS"
    )
  ) %>%
  dplyr::left_join(
    main_gene %>%
      dplyr::select(
        AptName, EntrezGeneSymbol,
        primary_logFC = logFC,
        primary_adj.P.Val = adj.P.Val
      ),
    by = c("AptName", "EntrezGeneSymbol")
  ) %>%
  dplyr::arrange(meta_adj.P.Val, meta_p)

meta_summary <- tibble::tibble(
  n_fixed_primary_somamers = nrow(meta_tbl),
  n_countries = length(countries_to_test),
  meta_fdr005_total = sum(meta_tbl$meta_adj.P.Val < MAIN_FDR, na.rm = TRUE),
  meta_fdr005_higher_in_AD = sum(
    meta_tbl$meta_adj.P.Val < MAIN_FDR & meta_tbl$meta_logFC > 0,
    na.rm = TRUE
  ),
  meta_fdr005_lower_in_AD = sum(
    meta_tbl$meta_adj.P.Val < MAIN_FDR & meta_tbl$meta_logFC < 0,
    na.rm = TRUE
  ),
  median_I2_among_meta_fdr005 = stats::median(
    meta_tbl$I2[meta_tbl$meta_adj.P.Val < MAIN_FDR], na.rm = TRUE
  ),
  mean_I2_among_meta_fdr005 = mean(
    meta_tbl$I2[meta_tbl$meta_adj.P.Val < MAIN_FDR], na.rm = TRUE
  )
)

safe_write_csv(country_specific_fixed, file.path(
  META_DIR, "country_specific_DEP_results_FIXED_PRIMARY_MAP.csv"
))
safe_write_csv(meta_tbl, file.path(
  META_DIR, "country_meta_analysis_results_FIXED_PRIMARY_MAP.csv"
))
safe_write_csv(meta_summary, file.path(
  META_DIR, "country_meta_summary_FIXED_PRIMARY_MAP.csv"
))

# =============================================================================
# Omnibus diagnosis-by-country interaction on the fixed map
# =============================================================================

interaction_data <- dep_df %>%
  dplyr::filter(as.character(Country) %in% countries_to_test) %>%
  dplyr::mutate(Country = droplevels(factor(Country)))
interaction_formula <- "~ SampleGroup * Country + Age + Sex + Education"
interaction_vars <- all.vars(stats::as.formula(interaction_formula))
interaction_keep <- stats::complete.cases(
  interaction_data[, interaction_vars, drop = FALSE]
)
interaction_data <- interaction_data[interaction_keep, , drop = FALSE]

interaction_expr <- interaction_data %>%
  dplyr::select(dplyr::all_of(intersect(dep_protein_cols, names(interaction_data)))) %>%
  as.matrix() %>%
  safe_log2_matrix() %>%
  t()
interaction_design <- stats::model.matrix(
  stats::as.formula(interaction_formula), data = interaction_data
)
if (qr(interaction_design)$rank != ncol(interaction_design)) {
  stop("Diagnosis-by-country interaction design is rank deficient.", call. = FALSE)
}
interaction_coef <- grep(
  "^SampleGroupAD:Country", colnames(interaction_design), value = TRUE
)
if (length(interaction_coef) == 0) {
  stop("No diagnosis-by-country interaction coefficients were found.")
}
interaction_fit <- limma::eBayes(limma::lmFit(interaction_expr, interaction_design))

interaction_omnibus_aptamer <- limma::topTable(
  interaction_fit,
  coef = interaction_coef,
  number = Inf,
  sort.by = "F",
  adjust.method = "BH"
) %>%
  tibble::rownames_to_column("feature_id_raw") %>%
  dplyr::mutate(AptName = as.character(feature_id_raw)) %>%
  dplyr::left_join(annot_tbl, by = "AptName") %>%
  dplyr::mutate(
    Protein_Name = make_protein_name(
      EntrezGeneSymbol, TargetFullName, Target, AptName, feature_id_raw
    )
  )

interaction_omnibus_fixed <- apply_primary_map(
  interaction_omnibus_aptamer, primary_map,
  "diagnosis_by_country_interaction_omnibus"
) %>%
  dplyr::mutate(
    adj.P.Val_aptamer_universe = adj.P.Val,
    adj.P.Val = stats::p.adjust(P.Value, method = "BH"),
    n_interaction_coefficients = length(interaction_coef),
    interaction_test = "Moderated omnibus F-test across all diagnosis-by-country coefficients"
  ) %>%
  dplyr::arrange(adj.P.Val, P.Value)

interaction_specific <- purrr::map_dfr(interaction_coef, function(coef_i) {
  limma::topTable(
    interaction_fit, coef = coef_i, number = Inf,
    sort.by = "P", adjust.method = "BH"
  ) %>%
    tibble::rownames_to_column("feature_id_raw") %>%
    dplyr::mutate(AptName = as.character(feature_id_raw), interaction_term = coef_i) %>%
    dplyr::left_join(annot_tbl, by = "AptName") %>%
    apply_primary_map(primary_map = primary_map, model_name = coef_i)
})

safe_write_csv(interaction_omnibus_fixed, file.path(
  INTERACTION_DIR,
  "country_interaction_omnibus_FIXED_PRIMARY_MAP.csv"
))
safe_write_csv(interaction_specific, file.path(
  INTERACTION_DIR,
  "country_interaction_specific_terms_FIXED_PRIMARY_MAP.csv"
))

# =============================================================================
# LOSO with fixed primary map
# =============================================================================

site_candidates <- c(
  "Site", "site", "Center", "center", "Cohort", "cohort",
  "RecruitmentSite", "recruitment_site", "site_id", "Site_ID"
)
site_var <- site_candidates[site_candidates %in% names(dep_df)][1]
if (length(site_var) == 0 || is.na(site_var)) site_var <- NA_character_

loso_summary <- tibble::tibble()
loso_gene_all <- tibble::tibble()
site_counts <- tibble::tibble()

if (!is.na(site_var)) {
  site_counts <- dep_df %>%
    dplyr::count(.data[[site_var]], SampleGroup, name = "n") %>%
    dplyr::rename(site = !!rlang::sym(site_var)) %>%
    tidyr::pivot_wider(names_from = SampleGroup, values_from = n, values_fill = 0) %>%
    ensure_count_cols(c("CN", "AD")) %>%
    dplyr::mutate(total = CN + AD) %>%
    dplyr::arrange(dplyr::desc(total))
  safe_write_csv(site_counts, file.path(
    analysis_root, "06_robustness", "site_robustness", "site_group_counts.csv"
  ))

  eligible_sites <- site_counts %>%
    dplyr::filter(CN >= MIN_N_PER_SITE_GROUP, AD >= MIN_N_PER_SITE_GROUP) %>%
    dplyr::pull(site) %>%
    as.character()

  loso_summary_list <- list()
  loso_gene_list <- list()

  for (site_i in eligible_sites) {
    message("LOSO fixed-map model excluding: ", site_i)
    dat_loso <- dep_df %>%
      dplyr::filter(as.character(.data[[site_var]]) != site_i) %>%
      dplyr::mutate(Country = droplevels(factor(Country)))

    fit_loso <- run_limma_dep(
      dat_loso,
      "~ SampleGroup + Age + Sex + Country + Education",
      "SampleGroupAD",
      paste0("LOSO_excluding_", safe_file_tag(site_i)),
      dep_protein_cols,
      annot_tbl
    )
    gene_loso <- apply_primary_map(
      fit_loso$dep, primary_map, paste0("LOSO_excluding_", site_i)
    ) %>%
      dplyr::mutate(excluded_site = site_i)
    cmp <- fixed_comparison(main_gene, gene_loso, "excluded_site", site_i)

    export_dep(gene_loso, file.path(
      LOSO_DIR,
      paste0("LOSO_excluding_", safe_file_tag(site_i),
             "_dep_results_FIXED_PRIMARY_MAP.csv")
    ))
    safe_write_csv(cmp, file.path(
      LOSO_DIR,
      paste0("LOSO_excluding_", safe_file_tag(site_i),
             "_comparison_to_main_FIXED_PRIMARY_MAP.csv")
    ))

    loso_summary_list[[site_i]] <- summarize_comparison(
      cmp, nrow(fit_loso$metadata), "excluded_site", site_i,
      nrow(primary_map)
    )
    loso_gene_list[[site_i]] <- gene_loso %>%
      dplyr::select(
        EntrezGeneSymbol, AptName, excluded_site,
        loso_logFC = logFC, loso_adj.P.Val = adj.P.Val
      )
  }

  loso_summary <- dplyr::bind_rows(loso_summary_list)
  loso_gene_all <- dplyr::bind_rows(loso_gene_list)
  safe_write_csv(loso_summary, file.path(
    LOSO_DIR, "LOSO_summary_metrics_FIXED_PRIMARY_MAP.csv"
  ))
}

# =============================================================================
# Balanced country resampling with fixed primary map
# =============================================================================

eligible_country_tbl <- country_counts %>%
  dplyr::filter(as.character(Country) %in% countries_to_test)
n_per_country_group <- min(
  eligible_country_tbl$CN,
  eligible_country_tbl$AD,
  BALANCED_RESAMPLING_MIN_GROUP,
  na.rm = TRUE
)
if (!is.finite(n_per_country_group) || n_per_country_group < 3) {
  stop("Balanced resampling has fewer than three participants per country/group.")
}

balanced_summary_list <- vector("list", BALANCED_RESAMPLING_NITER)
balanced_protein_list <- vector("list", BALANCED_RESAMPLING_NITER)

for (iter in seq_len(BALANCED_RESAMPLING_NITER)) {
  if (iter %% 20 == 0) message("Balanced fixed-map iteration ", iter, "/", BALANCED_RESAMPLING_NITER)
  set.seed(SEED_GLOBAL + iter)

  sampled_ids <- dep_df %>%
    dplyr::filter(as.character(Country) %in% countries_to_test) %>%
    dplyr::group_by(Country, SampleGroup) %>%
    dplyr::slice_sample(n = n_per_country_group, replace = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::pull(SampleId)

  dat_bal <- dep_df %>%
    dplyr::filter(SampleId %in% sampled_ids) %>%
    dplyr::mutate(Country = droplevels(factor(Country)))

  bal_fit <- tryCatch(
    run_limma_dep(
      dat_bal,
      "~ SampleGroup + Age + Sex + Country + Education",
      "SampleGroupAD",
      paste0("balanced_country_resampling_iter_", iter),
      dep_protein_cols,
      annot_tbl
    ),
    error = function(e) {
      warning("Balanced iteration ", iter, " failed: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(bal_fit)) next

  bal_gene <- apply_primary_map(
    bal_fit$dep, primary_map, paste0("balanced_iteration_", iter)
  )
  cmp <- fixed_comparison(main_gene, bal_gene, "iteration", iter)

  balanced_summary_list[[iter]] <- tibble::tibble(
    iteration = iter,
    n_samples = nrow(bal_fit$metadata),
    n_per_country_group = n_per_country_group,
    n_fixed_primary_somamers = nrow(primary_map),
    n_main_sig_fdr005 = sum(cmp$main_sig_fdr005, na.rm = TRUE),
    logFC_correlation = suppressWarnings(stats::cor(
      cmp$main_logFC, cmp$context_logFC, use = "complete.obs"
    )),
    direction_consistency_all = mean(cmp$same_direction, na.rm = TRUE),
    prop_main_sig_preserved = mean(
      cmp$main_sig_preserved[cmp$main_sig_fdr005], na.rm = TRUE
    )
  )

  balanced_protein_list[[iter]] <- cmp %>%
    dplyr::transmute(
      iteration = iter,
      EntrezGeneSymbol, AptName, Protein_Name,
      bal_logFC = context_logFC,
      bal_adj.P.Val = context_adj.P.Val,
      same_direction,
      main_sig_fdr005,
      main_sig_preserved,
      main_logFC,
      main_adj.P.Val
    )
}

balanced_summary <- dplyr::bind_rows(balanced_summary_list)
balanced_all <- dplyr::bind_rows(balanced_protein_list)
if (nrow(balanced_summary) == 0 || nrow(balanced_all) == 0) {
  stop("All balanced-resampling iterations failed.", call. = FALSE)
}

balanced_protein <- balanced_all %>%
  dplyr::group_by(EntrezGeneSymbol, AptName, Protein_Name) %>%
  dplyr::summarise(
    mean_bal_logFC = mean(bal_logFC, na.rm = TRUE),
    sd_bal_logFC = stats::sd(bal_logFC, na.rm = TRUE),
    prop_same_direction = mean(same_direction, na.rm = TRUE),
    prop_preserved_if_main_sig = mean(main_sig_preserved, na.rm = TRUE),
    prop_sig_fdr005 = mean(bal_adj.P.Val < MAIN_FDR, na.rm = TRUE),
    n_iterations = dplyr::n(),
    main_logFC = dplyr::first(main_logFC),
    main_adj.P.Val = dplyr::first(main_adj.P.Val),
    .groups = "drop"
  ) %>%
  dplyr::arrange(main_adj.P.Val)

safe_write_csv(balanced_summary, file.path(
  BALANCED_DIR, "balanced_resampling_summary_metrics_FIXED_PRIMARY_MAP.csv"
))
safe_write_csv(balanced_protein, file.path(
  BALANCED_DIR, "balanced_resampling_protein_stability_FIXED_PRIMARY_MAP.csv"
))

# =============================================================================
# Fixed-map robustness classification
# =============================================================================

robustness_classification <- main_gene %>%
  dplyr::select(
    AptName, Protein_Name, EntrezGeneSymbol,
    main_logFC = logFC, main_adj.P.Val = adj.P.Val, main_type = type
  ) %>%
  dplyr::mutate(main_sig_fdr005 = main_adj.P.Val < MAIN_FDR) %>%
  dplyr::left_join(
    main_vs_mean_loco %>%
      dplyr::select(
        AptName, EntrezGeneSymbol, mean_loco_logFC,
        prop_same_direction_loco = prop_same_direction,
        prop_sig_fdr005_loco = prop_sig_fdr005
      ),
    by = c("AptName", "EntrezGeneSymbol")
  ) %>%
  dplyr::left_join(
    meta_tbl %>%
      dplyr::select(AptName, EntrezGeneSymbol, meta_logFC, meta_adj.P.Val, I2),
    by = c("AptName", "EntrezGeneSymbol")
  ) %>%
  dplyr::left_join(
    balanced_protein %>%
      dplyr::select(
        AptName, EntrezGeneSymbol, mean_bal_logFC,
        prop_same_direction_balanced = prop_same_direction,
        prop_sig_fdr005_balanced = prop_sig_fdr005
      ),
    by = c("AptName", "EntrezGeneSymbol")
  ) %>%
  dplyr::mutate(
    loco_preserved =
      !is.na(mean_loco_logFC) &
      sign(main_logFC) == sign(mean_loco_logFC) &
      dplyr::coalesce(prop_same_direction_loco, 0) >= 0.80,
    meta_preserved =
      !is.na(meta_logFC) &
      sign(main_logFC) == sign(meta_logFC) &
      dplyr::coalesce(meta_adj.P.Val, 1) < MAIN_FDR,
    balanced_preserved =
      !is.na(mean_bal_logFC) &
      sign(main_logFC) == sign(mean_bal_logFC) &
      dplyr::coalesce(prop_same_direction_balanced, 0) >= 0.80,
    robustness_score =
      as.integer(main_sig_fdr005) +
      as.integer(loco_preserved) +
      as.integer(meta_preserved) +
      as.integer(balanced_preserved),
    robustness_class = dplyr::case_when(
      robustness_score >= 4 ~ "High robustness",
      robustness_score == 3 ~ "Moderate robustness",
      robustness_score == 2 ~ "Partial robustness",
      TRUE ~ "Limited or unclassified"
    ),
    n_fixed_primary_somamers = nrow(primary_map)
  ) %>%
  dplyr::arrange(dplyr::desc(robustness_score), main_adj.P.Val)

robustness_counts <- robustness_classification %>%
  dplyr::count(robustness_class, name = "n") %>%
  dplyr::mutate(n_fixed_primary_somamers = nrow(primary_map)) %>%
  dplyr::arrange(dplyr::desc(n))

safe_write_csv(robustness_classification, file.path(
  CLASS_DIR, "protein_robustness_classification_FIXED_PRIMARY_MAP.csv"
))
safe_write_csv(robustness_counts, file.path(
  CLASS_DIR, "protein_robustness_classification_counts_FIXED_PRIMARY_MAP.csv"
))

# =============================================================================
# Workspace update and hard audits
# =============================================================================

hard_audit <- tibble::tibble(
  output = c(
    "primary_map", "LOCO_protein_table", "country_meta", "LOSO_summary",
    "balanced_protein", "robustness_classification"
  ),
  observed_rows = c(
    nrow(primary_map), nrow(main_vs_mean_loco), nrow(meta_tbl),
    if (nrow(loso_summary) > 0) nrow(loso_summary) else NA_integer_,
    nrow(balanced_protein), nrow(robustness_classification)
  ),
  expected_map_rows = c(
    nrow(primary_map), nrow(primary_map), nrow(primary_map),
    NA_integer_, nrow(primary_map), nrow(primary_map)
  ),
  passed = c(
    nrow(primary_map) == 9638L,
    nrow(main_vs_mean_loco) == nrow(primary_map),
    nrow(meta_tbl) == nrow(primary_map),
    if (nrow(loso_summary) > 0) all(loso_summary$n_fixed_primary_somamers == nrow(primary_map)) else FALSE,
    nrow(balanced_protein) == nrow(primary_map),
    nrow(robustness_classification) == nrow(primary_map)
  )
)
if (!all(hard_audit$passed)) {
  safe_write_csv(hard_audit, file.path(MANIFEST_DIR, "hard_audit_FAILED.csv"))
  stop("One or more fixed-map robustness audits failed.", call. = FALSE)
}

safe_write_csv(hard_audit, file.path(MANIFEST_DIR, "hard_audit.csv"))
run_manifest <- tibble::tibble(
  item = c(
    "run_time", "workspace_file", "primary_map_rows", "primary_map_md5",
    "eligible_countries", "eligible_country_count", "site_variable",
    "eligible_site_count", "balanced_iterations_requested",
    "balanced_iterations_completed", "n_per_country_group",
    "main_DEP_FDR005", "meta_FDR005", "meta_higher", "meta_lower",
    "interpretation"
  ),
  value = c(
    as.character(Sys.time()), workspace_file, as.character(nrow(primary_map)),
    map_md5, paste(countries_to_test, collapse = ", "),
    as.character(length(countries_to_test)),
    ifelse(is.na(site_var), "none", site_var),
    as.character(if (nrow(loso_summary) > 0) nrow(loso_summary) else 0),
    as.character(BALANCED_RESAMPLING_NITER),
    as.character(nrow(balanced_summary)),
    as.character(n_per_country_group),
    as.character(sum(main_gene$adj.P.Val < MAIN_FDR, na.rm = TRUE)),
    as.character(meta_summary$meta_fdr005_total),
    as.character(meta_summary$meta_fdr005_higher_in_AD),
    as.character(meta_summary$meta_fdr005_lower_in_AD),
    "Internal geographic and recruitment-context stability; not external validation"
  )
)
safe_write_csv(run_manifest, file.path(MANIFEST_DIR, "fixed_map_robustness_run_manifest.csv"))
writeLines(
  capture.output(utils::sessionInfo()),
  file.path(MANIFEST_DIR, "fixed_map_robustness_sessionInfo.txt")
)

if (UPDATE_WORKSPACE) {
  run_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  backup_file <- file.path(
    BACKUP_DIR,
    paste0("analysis_workspace_before_fixed_robustness_v2_", run_stamp, ".RData")
  )
  file.copy(workspace_file, backup_file, overwrite = FALSE)

  ws$primary_gene_map_fixed <- primary_map
  ws$loco_summary_metrics_fixed <- loco_summary
  ws$main_vs_loco_mean_fixed <- main_vs_mean_loco
  ws$country_specific_fixed <- country_specific_fixed
  ws$country_meta_fixed <- meta_tbl
  ws$country_meta_summary_fixed <- meta_summary
  ws$country_interaction_omnibus_fixed <- interaction_omnibus_fixed
  ws$loso_summary_metrics_fixed <- loso_summary
  ws$balanced_summary_fixed <- balanced_summary
  ws$balanced_protein_fixed <- balanced_protein
  ws$robustness_classification_fixed <- robustness_classification
  ws$robustness_counts_fixed <- robustness_counts
  save(list = ls(ws), file = workspace_file, envir = ws)
}

message("Fixed-primary-map recruitment-context robustness analyses completed.")
message("Fixed map rows: ", nrow(primary_map))
message("Main DEP FDR < 0.05: ", sum(main_gene$adj.P.Val < MAIN_FDR, na.rm = TRUE))
message("Country meta-analysis FDR < 0.05: ", meta_summary$meta_fdr005_total)
message("LOCO outputs: ", LOCO_DIR)
message("LOSO outputs: ", LOSO_DIR)
message("Balanced outputs: ", BALANCED_DIR)
###############################################################################
# END
###############################################################################


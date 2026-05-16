###############################################################################
# ReDLat plasma proteomics differential analysis
# 01_data_processing_and_differential_analysis.R
#
# Purpose
# - Run the primary statistical analysis layer for the ReDLat SOMAscan plasma
#   proteomics study of Alzheimer's disease (AD) versus cognitively normal (CN)
#   participants.
# - Use the internal ADAT annotation as the primary analyte annotation source.
# - Define the human protein universe from the internal ADAT annotation and use
#   the gene-collapsed universe as the default level for biological interpretation.
# - Export publication-ready source tables, model diagnostics, robustness tables,
#   corrected enrichment outputs, and compatibility aliases for downstream figure
#   and supplementary scripts.
#
# Workflow
# 01 Import metadata and ADAT, harmonize sample identifiers, and build annotation.
# 02 Define the internal human protein universe and sample-tracking tables.
# 03 Generate normalized log2 expression objects and PCA source tables.
# 04 Fit the primary AD versus CN differential expression model.
# 05 Fit sensitivity models: APOE4, within-AD CDR-SB severity, secondary CDR-SB-adjusted diagnostic attenuation, and vascular/metabolic covariates.
# 06 Fit secondary clinical and cognitive severity association models.
# 07 Run corrected enrichment analyses: GO, KEGG, Reactome, Hallmark, and ORA.
# 08 Run internal robustness analyses: LOCO, country meta-analysis, interaction,
#    LOSO, balanced country resampling, and formal robustness classification.
# 09 Export canonical tables, backward-compatible aliases, manifest, diagnostics,
#    session information, and analysis workspace.
#
# Downstream scripts
# - 02_main_figure_generation.R builds manuscript figures from these outputs.
# - 03_supplementary_and_robustness_analyses.R runs exploratory and extended
#   supplementary analyses that are not part of the primary inference layer.
###############################################################################
###############################################################################

###############################################################################
# 00_config_paths_packages
###############################################################################

project_root <- "C:/Users/mnpiz/Desktop/DEPs_Proteomic_Publishable_V2"
outdir <- project_root

csv_file  <- file.path(project_root, "data", "ReDLat_CARD-proteomic_updated_all_data_11_2025.csv")
adat_file <- file.path(project_root, "data", "merged.adat")

MAIN_GROUPS <- c("CN", "AD")
MAIN_FDR <- 0.05
STRICT_FDR <- 0.01
LOGFC_THRESHOLD <- 0
MIN_SAMPLES_COR <- 10
MIN_N_PER_GROUP_COUNTRY <- 10
MIN_N_PER_SITE_GROUP <- 8
BALANCED_RESAMPLING_NITER <- 200
BALANCED_RESAMPLING_MIN_GROUP <- 8
SEED_GLOBAL <- 1234

RUN_FLAGS <- list(
  run_pca = TRUE,
  run_main_dep = TRUE,
  run_apoe_sensitivity = TRUE,
  run_cdrsb_sensitivity = TRUE,
  run_cdrsb_adjusted_diagnostic_sensitivity = TRUE,
  run_vascular_metabolic_sensitivity = TRUE,
  run_atn_adjusted_sensitivity = TRUE,
  run_secondary_clinical_severity = TRUE,
  run_corrected_enrichment = TRUE,
  run_country_loco = TRUE,
  run_country_meta = TRUE,
  run_country_interaction = TRUE,
  run_site_robustness = TRUE,
  run_balanced_country_resampling = TRUE,
  run_formal_robustness_classification = TRUE
)

set.seed(SEED_GLOBAL)
options(stringsAsFactors = FALSE)

SCRIPT_ID <- "01_data_processing_and_differential_analysis"
SCRIPT_ROLE <- "Primary data processing, differential expression, corrected enrichment, and internal robustness analyses"

required_pkgs <- c(
  "SomaDataIO", "dplyr", "tidyr", "purrr", "ggplot2", "ggrepel",
  "tibble", "stringr", "forcats", "readr", "patchwork", "limma",
  "ggpubr", "openxlsx", "metafor", "msigdbr", "rlang",
  "clusterProfiler", "org.Hs.eg.db", "BiocParallel", "ReactomePA", "DOSE"
)

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing required packages: ", paste(missing_pkgs, collapse = ", "),
    "\nInstall them before running this publication pipeline."
  )
}
invisible(lapply(required_pkgs, library, character.only = TRUE))

subdirs <- c(
  "data", "data/raw", "result", "result/00_logs",
  "result/01_qc", "result/01_qc/sample_tracking", "result/01_qc/annotations",
  "result/02_pca",
  "result/03_dep", "result/03_dep/aptamer_level", "result/03_dep/gene_collapsed",
  "result/04_sensitivity", "result/04_sensitivity/apoe", "result/04_sensitivity/cdrsb",
  "result/04_sensitivity/cdrsb/AD_only", "result/04_sensitivity/cdrsb/AD_vs_CN_adjusted",
  "result/04_sensitivity/vascular_metabolic", "result/04_sensitivity/atn_adjusted", "result/04_sensitivity/clinical_severity",
  "result/05_enrichment_corrected", "result/05_enrichment_corrected/gsea", "result/05_enrichment_corrected/ora",
  "result/06_robustness", "result/06_robustness/country_loco/tables",
  "result/06_robustness/country_meta/tables", "result/06_robustness/country_interaction/tables",
  "result/06_robustness/site_robustness", "result/06_robustness/site_robustness/loso/tables",
  "result/06_robustness/balanced_country_resampling/tables",
  "result/06_robustness/formal_classification",
  "result/07_manifest", "result/07_manifest/model_diagnostics", "result/workspace",
  "result/final_source_data"
)
invisible(lapply(file.path(outdir, subdirs), dir.create, recursive = TRUE, showWarnings = FALSE))
setwd(outdir)

missing_input_files <- c(csv_file, adat_file)[!file.exists(c(csv_file, adat_file))]
if (length(missing_input_files) > 0) {
  stop("Missing required input files:\n", paste(missing_input_files, collapse = "\n"))
}

###############################################################################
# 00_helpers
###############################################################################

safe_write_csv <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(x, file)
}

safe_copy <- function(from, to) {
  if (file.exists(from)) {
    dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
    file.copy(from, to, overwrite = TRUE)
  }
}

safe_file_tag <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
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

safe_log2_vector <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[x <= 0] <- NA_real_
  log2(x)
}

safe_log2_matrix <- function(mat) {
  mat <- apply(mat, 2, as.numeric)
  mat[mat <= 0] <- NA_real_
  log2(mat)
}

safe_standardize_vector <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  sdx <- stats::sd(x, na.rm = TRUE)
  mn <- mean(x, na.rm = TRUE)
  if (is.na(sdx) || sdx == 0) return(rep(0, length(x)))
  (x - mn) / sdx
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

analysis_checkpoints <- tibble::tibble(step = character(), n = integer(), value = numeric(), note = character())
add_checkpoint <- function(checkpoint_tbl, step, n = NA_integer_, value = NA_real_, note = NA_character_) {
  dplyr::bind_rows(
    checkpoint_tbl,
    tibble::tibble(step = step, n = as.integer(n), value = as.numeric(value), note = as.character(note))
  )
}

model_formula_manifest <- tibble::tibble(
  model = character(), formula = character(), coefficient = character(),
  n_samples = integer(), n_proteins = integer(), note = character()
)
add_model_formula <- function(manifest, model, formula, coefficient, n_samples, n_proteins, note = NA_character_) {
  dplyr::bind_rows(
    manifest,
    tibble::tibble(
      model = as.character(model), formula = as.character(formula), coefficient = as.character(coefficient),
      n_samples = as.integer(n_samples), n_proteins = as.integer(n_proteins), note = as.character(note)
    )
  )
}

model_design_diagnostics <- tibble::tibble(
  model = character(), n_samples = integer(), n_parameters = integer(), design_rank = integer(),
  full_rank = logical(), columns = character(), alias_diagnostic = character()
)

validate_design_matrix <- function(design, model_name = "model") {
  rank_design <- qr(design)$rank
  n_col <- ncol(design)
  out <- tibble::tibble(
    model = model_name,
    n_samples = nrow(design),
    n_parameters = n_col,
    design_rank = rank_design,
    full_rank = rank_design == n_col,
    columns = paste(colnames(design), collapse = " | ")
  )
  out$alias_diagnostic <- NA_character_
  if (rank_design < n_col) {
    out$alias_diagnostic <- "Rank deficient design; inspect sparse factors/covariates."
    warning(model_name, " design matrix is rank deficient. See model_design_diagnostics.csv.")
  }
  out
}

prepare_optional_model_covariates <- function(df, covariates) {
  out <- df
  for (v in covariates) {
    if (!v %in% names(out)) next
    x <- out[[v]]
    if (is.character(x) || is.factor(x) || is.logical(x)) {
      out[[v]] <- factor(clean_text_na(x))
    } else if (is.numeric(x)) {
      ux <- sort(unique(stats::na.omit(x)))
      if (length(ux) <= 6 && all(ux %in% c(0, 1, 2, 3, 4, 5))) {
        out[[v]] <- factor(x)
      } else {
        out[[v]] <- as.numeric(x)
      }
    } else {
      out[[v]] <- factor(as.character(x))
    }
  }
  out
}

ensure_count_cols <- function(tbl, group_cols = c("CN", "AD")) {
  for (g in group_cols) if (!g %in% names(tbl)) tbl[[g]] <- 0L
  tbl
}

###############################################################################
# Annotation and DEP helpers
###############################################################################

extract_adat_annotation_from_table <- function(adat_file, max_annotation_lines = 500) {
  lines <- readLines(adat_file, warn = FALSE)
  table_begin <- grep("^\\^TABLE_BEGIN", lines)
  if (length(table_begin) == 0) stop("Could not find ^TABLE_BEGIN in the ADAT file.")
  table_begin <- table_begin[1]
  candidate_end <- min(length(lines), table_begin + max_annotation_lines)
  meta_lines <- lines[(table_begin + 1):candidate_end]
  split_lines <- strsplit(meta_lines, "\t", fixed = TRUE)
  find_row <- function(row_name) {
    idx <- which(vapply(split_lines, function(x) length(x) > 0 && row_name %in% x, logical(1)))
    if (length(idx) == 0) return(NULL)
    x <- split_lines[[idx[1]]]
    pos <- which(x == row_name)[1]
    x[(pos + 1):length(x)]
  }
  seq_id <- find_row("SeqId")
  if (is.null(seq_id)) {
    adat_obj <- tryCatch(SomaDataIO::read_adat(adat_file), error = function(e) NULL)
    ai <- tryCatch(SomaDataIO::getAnalyteInfo(adat_obj), error = function(e) NULL)
    if (!is.null(ai) && all(c("AptName", "SeqId") %in% names(ai))) {
      return(tibble::as_tibble(ai) %>% dplyr::mutate(dplyr::across(dplyr::everything(), as.character)))
    }
    stop("Could not find SeqId row in ADAT annotation block and SomaDataIO fallback failed.")
  }
  n <- length(seq_id)
  safe_row <- function(row_name) {
    vals <- find_row(row_name)
    if (is.null(vals)) return(rep(NA_character_, n))
    length(vals) <- n
    vals
  }
  tibble::tibble(
    SeqId = seq_id,
    AptName = paste0("seq.", gsub("-", ".", seq_id)),
    SeqIdVersion = safe_row("SeqIdVersion"),
    SomaId = safe_row("SomaId"),
    TargetFullName = safe_row("TargetFullName"),
    Target = safe_row("Target"),
    UniProt = safe_row("UniProt"),
    EntrezGeneID = safe_row("EntrezGeneID"),
    EntrezGeneSymbol = safe_row("EntrezGeneSymbol"),
    Organism = safe_row("Organism"),
    Units = safe_row("Units"),
    Type = safe_row("Type"),
    Dilution = safe_row("Dilution")
  ) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ trimws(as.character(.x)))) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ dplyr::na_if(.x, "")))
}

build_annotation_table <- function(soma_info_internal, seq_cols) {
  seq_key_tbl <- tibble::tibble(AptName = as.character(seq_cols))
  internal_tbl <- soma_info_internal %>%
    dplyr::mutate(
      AptName = as.character(AptName), SeqId = as.character(SeqId), SomaId = as.character(SomaId),
      TargetFullName = as.character(TargetFullName), Target = as.character(Target),
      EntrezGeneID = suppressWarnings(as.numeric(EntrezGeneID)),
      EntrezGeneSymbol = as.character(EntrezGeneSymbol), UniProt = as.character(UniProt),
      Organism = as.character(Organism), Type = as.character(Type)
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

summarize_dep_counts <- function(dep_tbl, fdr_values = c(STRICT_FDR, MAIN_FDR), universe_label = "aptamer_level") {
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

export_dep_table <- function(dep_tbl, file) {
  keep_cols <- intersect(
    c("Protein_Name", "EntrezGeneSymbol", "TargetFullName", "Target", "UniProt", "AptName", "SeqId", "feature_id_raw",
      "logFC", "se", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "type", "Direction"),
    names(dep_tbl)
  )
  dep_tbl %>% dplyr::select(dplyr::all_of(keep_cols)) %>% dplyr::arrange(adj.P.Val) %>% safe_write_csv(file)
}

export_fdr_specific_dep_tables <- function(dep_tbl, out_folder, file_prefix, fdr_values = c(STRICT_FDR, MAIN_FDR)) {
  # Exports filtered DEP tables for manuscript/reporting convenience.
  # The full limma result remains the canonical source of truth.
  dir.create(out_folder, recursive = TRUE, showWarnings = FALSE)
  purrr::walk(fdr_values, function(fdr_value) {
    tag <- ifelse(abs(fdr_value - 0.01) < 1e-12, "FDR001",
                  ifelse(abs(fdr_value - 0.05) < 1e-12, "FDR005",
                         paste0("FDR", gsub("\\.", "", as.character(fdr_value)))))
    export_dep_table(
      dep_tbl %>% dplyr::filter(adj.P.Val < fdr_value),
      file.path(out_folder, paste0(file_prefix, "_", tag, ".csv"))
    )
  })
  invisible(TRUE)
}

run_limma_dep_model <- function(dat, seq_cols, annot_tbl, formula_str, coef_name = "SampleGroupAD", fdr = MAIN_FDR, model_name = "limma_model") {
  dat <- dat %>% tibble::as_tibble()
  protein_cols <- intersect(seq_cols, names(dat))
  if (length(protein_cols) == 0) stop("No protein columns found for limma model.")
  metadata <- dat %>% dplyr::select(-dplyr::all_of(protein_cols))
  model_vars <- all.vars(stats::as.formula(formula_str))
  missing_model_vars <- setdiff(model_vars, names(metadata))
  if (length(missing_model_vars) > 0) stop("Missing model variables in ", model_name, ": ", paste(missing_model_vars, collapse = ", "))
  keep <- stats::complete.cases(metadata[, model_vars, drop = FALSE])
  dat <- dat[keep, , drop = FALSE]
  metadata <- metadata[keep, , drop = FALSE]
  expr <- dat %>% dplyr::select(dplyr::all_of(protein_cols)) %>% as.matrix()
  expr <- t(safe_log2_matrix(expr))
  design <- model.matrix(stats::as.formula(formula_str), data = metadata)
  design_diag <- validate_design_matrix(design, model_name = model_name)
  if (!coef_name %in% colnames(design)) stop("Coefficient not found in design matrix for ", model_name, ": ", coef_name)
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

build_gene_compare_table <- function(dep_primary_gene, dep_secondary_gene, secondary_label = "secondary") {
  dep_primary_gene %>%
    dplyr::select(EntrezGeneSymbol, AptName, SeqId, Protein_Name, logFC, adj.P.Val, P.Value, type) %>%
    dplyr::rename(
      AptName_primary = AptName, SeqId_primary = SeqId, Protein_Name_primary = Protein_Name,
      logFC_primary = logFC, adj.P.Val_primary = adj.P.Val, P.Value_primary = P.Value, type_primary = type
    ) %>%
    dplyr::inner_join(
      dep_secondary_gene %>%
        dplyr::select(EntrezGeneSymbol, AptName, SeqId, Protein_Name, logFC, adj.P.Val, P.Value, type) %>%
        dplyr::rename(
          AptName_secondary = AptName, SeqId_secondary = SeqId, Protein_Name_secondary = Protein_Name,
          logFC_secondary = logFC, adj.P.Val_secondary = adj.P.Val, P.Value_secondary = P.Value, type_secondary = type
        ),
      by = "EntrezGeneSymbol"
    ) %>%
    dplyr::mutate(
      comparison = secondary_label,
      Protein_Name = dplyr::coalesce(Protein_Name_primary, Protein_Name_secondary),
      AptName = dplyr::coalesce(AptName_primary, AptName_secondary),
      SeqId = dplyr::coalesce(SeqId_primary, SeqId_secondary),
      same_direction = sign(logFC_primary) == sign(logFC_secondary),
      preserved_fdr005 = same_direction & adj.P.Val_secondary < MAIN_FDR,
      delta_logFC = logFC_secondary - logFC_primary,
      abs_ratio = abs(logFC_secondary) / pmax(abs(logFC_primary), 1e-9)
    )
}

make_summary_table <- function(df) {
  df %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ dplyr::na_if(.x, -1))) %>%
    dplyr::summarise(
      N = dplyr::n(),
      Male = sum(Sex == "M", na.rm = TRUE),
      Female = sum(Sex == "F", na.rm = TRUE),
      Age_mean = round(mean(Age, na.rm = TRUE), 1),
      Age_sd = round(stats::sd(Age, na.rm = TRUE), 1),
      Education_mean = round(mean(Education, na.rm = TRUE), 1),
      Education_sd = round(stats::sd(Education, na.rm = TRUE), 1),
      CDRSB_mean = if ("cdr_boxscore" %in% names(df)) round(mean(cdr_boxscore, na.rm = TRUE), 1) else NA_real_,
      CDRSB_sd = if ("cdr_boxscore" %in% names(df)) round(stats::sd(cdr_boxscore, na.rm = TRUE), 1) else NA_real_,
      MMSE_mean = if ("mmse_total" %in% names(df)) round(mean(mmse_total, na.rm = TRUE), 1) else NA_real_,
      MMSE_sd = if ("mmse_total" %in% names(df)) round(stats::sd(mmse_total, na.rm = TRUE), 1) else NA_real_,
      APOE_e4_carrier = sum(trimws(ApoE) %in% c("e2/e4", "e3/e4", "e4/e4"), na.rm = TRUE),
      APOE_non_e4 = sum(trimws(ApoE) %in% c("e2/e2", "e2/e3", "e3/e3"), na.rm = TRUE)
    )
}

detect_vascular_metabolic_covariates <- function(df) {
  candidate_map <- list(
    diabetes = c("contains_diabetes", "diabetes", "Diabetes", "diabetes_status"),
    hypertension = c("hypertension_status", "hypertension", "Hypertension"),
    obesity = c("obesity_status", "obesity", "Obesity"),
    bmi = c("BMI", "bmi", "clinical_bmi"),
    smoking = c("smoking_index", "smoking", "Smoking"),
    alcohol = c("alcohol_status", "Alcoholism_Status", "alcohol")
  )
  detected <- purrr::map_chr(candidate_map, function(cands) {
    hit <- cands[cands %in% names(df)][1]
    ifelse(length(hit) == 0 || is.na(hit), NA_character_, hit)
  })
  unique(as.character(stats::na.omit(detected)))
}

###############################################################################
# Enrichment helpers
###############################################################################


detect_atn_covariates <- function(df) {
  # Detect canonical plasma AT(N) variables after metadata harmonization.
  # p-tau217, NfL and Aβ42/40 are prioritized because they index tau pathology,
  # neurodegeneration and amyloid-related signal while limiting collinearity.
  candidate_map <- list(
    p_tau217 = c("p.tau217", "p_tau217", "p-tau217", "ptau217"),
    NfL = c("NfL", "NFL", "nfl", "Neurofilament_light"),
    ratio_AB42_40 = c("ratio.AB42.40", "ratio_AB42_40", "ratio AB42/40", "AB42_40_ratio", "Aβ42/40"),
    p_tau181 = c("p.tau181", "p_tau181", "p-tau181", "ptau181")
  )
  detected <- purrr::map_chr(candidate_map, function(cands) {
    hit <- cands[cands %in% names(df)][1]
    ifelse(length(hit) == 0 || is.na(hit), NA_character_, hit)
  })
  detected <- detected[!is.na(detected)]
  detected <- detected[vapply(detected, function(v) {
    x <- suppressWarnings(as.numeric(df[[v]]))
    sum(!is.na(x)) >= 20 && is.finite(stats::sd(x, na.rm = TRUE)) && stats::sd(x, na.rm = TRUE) > 0
  }, logical(1))]
  unique(as.character(detected))
}

prepare_ranked_genes <- function(dep_tbl, rank_col = "logFC") {
  dep_tbl %>%
    dplyr::mutate(
      EntrezGeneID = suppressWarnings(as.numeric(EntrezGeneID)),
      rank_metric = suppressWarnings(as.numeric(.data[[rank_col]]))
    ) %>%
    dplyr::filter(!is.na(EntrezGeneID), !is.na(rank_metric), is.finite(rank_metric)) %>%
    dplyr::group_by(EntrezGeneID) %>%
    dplyr::slice_max(order_by = abs(rank_metric), n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(dplyr::desc(rank_metric))
}

run_gsea_go_kegg <- function(dep_tbl, out_prefix, out_folder) {
  ranked_df <- prepare_ranked_genes(dep_tbl)
  if (nrow(ranked_df) < 10) return(list(go = data.frame(), kegg = data.frame(), ranked_df = ranked_df))
  gene_list <- ranked_df$rank_metric
  names(gene_list) <- as.character(ranked_df$EntrezGeneID)
  gene_list <- sort(gene_list, decreasing = TRUE)
  suppressWarnings(BiocParallel::register(BiocParallel::SerialParam(), default = TRUE))
  go <- tryCatch(clusterProfiler::gseGO(geneList = gene_list, ont = "ALL", keyType = "ENTREZID", pvalueCutoff = 1, pAdjustMethod = "BH", OrgDb = org.Hs.eg.db::org.Hs.eg.db, verbose = FALSE), error = function(e) NULL)
  kegg <- tryCatch(clusterProfiler::gseKEGG(geneList = gene_list, organism = "hsa", pvalueCutoff = 1, pAdjustMethod = "BH", verbose = FALSE), error = function(e) NULL)
  go_df <- if (!is.null(go)) as.data.frame(go) else data.frame()
  kegg_df <- if (!is.null(kegg)) as.data.frame(kegg) else data.frame()
  safe_write_csv(go_df, file.path(out_folder, paste0(out_prefix, "_gsea_go_bh.csv")))
  safe_write_csv(kegg_df, file.path(out_folder, paste0(out_prefix, "_gsea_kegg_bh.csv")))
  list(go = go_df, kegg = kegg_df, ranked_df = ranked_df)
}

run_gsea_reactome <- function(dep_tbl, out_prefix, out_folder) {
  ranked_df <- prepare_ranked_genes(dep_tbl)
  if (nrow(ranked_df) < 10) return(list(result = data.frame(), ranked_df = ranked_df))
  gene_list <- ranked_df$rank_metric
  names(gene_list) <- as.character(ranked_df$EntrezGeneID)
  gene_list <- sort(gene_list, decreasing = TRUE)
  suppressWarnings(BiocParallel::register(BiocParallel::SerialParam(), default = TRUE))
  res <- tryCatch(ReactomePA::gsePathway(geneList = gene_list, organism = "human", pvalueCutoff = 1, pAdjustMethod = "BH", verbose = FALSE), error = function(e) NULL)
  df <- if (!is.null(res)) as.data.frame(res) else data.frame()
  safe_write_csv(df, file.path(out_folder, paste0(out_prefix, "_gsea_reactome_bh.csv")))
  list(result = df, ranked_df = ranked_df)
}

run_gsea_hallmark <- function(dep_tbl, out_prefix, out_folder) {
  ranked_df <- prepare_ranked_genes(dep_tbl)
  if (nrow(ranked_df) < 10) return(list(result = data.frame(), ranked_df = ranked_df))
  hallmark_df <- tryCatch(msigdbr::msigdbr(species = "Homo sapiens", category = "H"), error = function(e) NULL)
  if (is.null(hallmark_df) || nrow(hallmark_df) == 0) return(list(result = data.frame(), ranked_df = ranked_df))
  hallmark_t2g <- hallmark_df %>%
    dplyr::select(gs_name, entrez_gene) %>%
    dplyr::distinct() %>%
    dplyr::mutate(entrez_gene = as.character(entrez_gene))
  gene_list <- ranked_df$rank_metric
  names(gene_list) <- as.character(ranked_df$EntrezGeneID)
  gene_list <- sort(gene_list, decreasing = TRUE)
  suppressWarnings(BiocParallel::register(BiocParallel::SerialParam(), default = TRUE))
  res <- tryCatch(clusterProfiler::GSEA(geneList = gene_list, TERM2GENE = hallmark_t2g, pvalueCutoff = 1, pAdjustMethod = "BH", verbose = FALSE), error = function(e) NULL)
  df <- if (!is.null(res)) as.data.frame(res) else data.frame()
  safe_write_csv(df, file.path(out_folder, paste0(out_prefix, "_gsea_hallmark_bh.csv")))
  list(result = df, ranked_df = ranked_df)
}

run_corrected_ora_by_direction <- function(dep_tbl, out_prefix, out_folder, fdr = MAIN_FDR) {
  dir.create(out_folder, recursive = TRUE, showWarnings = FALSE)
  mapped_df <- dep_tbl %>%
    dplyr::mutate(EntrezGeneID = suppressWarnings(as.numeric(EntrezGeneID))) %>%
    dplyr::filter(!is.na(EntrezGeneID)) %>%
    dplyr::group_by(EntrezGeneID) %>%
    dplyr::arrange(adj.P.Val, dplyr::desc(abs(logFC)), .by_group = TRUE) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup()
  universe_ids <- unique(as.character(mapped_df$EntrezGeneID))
  direction_sets <- list(
    higher_in_AD = mapped_df %>% dplyr::filter(logFC > 0, adj.P.Val < fdr) %>% dplyr::pull(EntrezGeneID) %>% as.character() %>% unique(),
    lower_in_AD  = mapped_df %>% dplyr::filter(logFC < 0, adj.P.Val < fdr) %>% dplyr::pull(EntrezGeneID) %>% as.character() %>% unique()
  )
  trace_tbl <- tibble::tibble(
    analysis = out_prefix,
    fdr_used = fdr,
    n_unique_entrez_universe = length(universe_ids),
    n_higher_in_AD = length(direction_sets$higher_in_AD),
    n_lower_in_AD = length(direction_sets$lower_in_AD),
    note = "ORA uses gene-collapsed DEP universe and BH correction; direction sets are defined at MAIN_FDR."
  )
  safe_write_csv(trace_tbl, file.path(out_folder, paste0(out_prefix, "_directional_ORA_trace.csv")))
  run_one <- function(ids, direction) {
    if (length(ids) < 5) {
      skipped <- tibble::tibble(direction = direction, note = paste0("Skipped: fewer than 5 genes at FDR < ", fdr))
      safe_write_csv(skipped, file.path(out_folder, paste0(out_prefix, "_", direction, "_ORA_skipped.csv")))
      return(list(go = data.frame(), kegg = data.frame(), reactome = data.frame()))
    }
    go <- tryCatch(clusterProfiler::enrichGO(gene = ids, universe = universe_ids, OrgDb = org.Hs.eg.db::org.Hs.eg.db, ont = "ALL", keyType = "ENTREZID", pAdjustMethod = "BH", readable = TRUE), error = function(e) NULL)
    kegg <- tryCatch(clusterProfiler::enrichKEGG(gene = ids, universe = universe_ids, organism = "hsa", pAdjustMethod = "BH"), error = function(e) NULL)
    reactome <- tryCatch(ReactomePA::enrichPathway(gene = ids, universe = universe_ids, organism = "human", pAdjustMethod = "BH", readable = TRUE), error = function(e) NULL)
    go_df <- if (!is.null(go)) as.data.frame(go) else data.frame()
    kegg_df <- if (!is.null(kegg)) as.data.frame(kegg) else data.frame()
    reactome_df <- if (!is.null(reactome)) as.data.frame(reactome) else data.frame()
    if (nrow(go_df) > 0) go_df$direction <- direction
    if (nrow(kegg_df) > 0) kegg_df$direction <- direction
    if (nrow(reactome_df) > 0) reactome_df$direction <- direction
    safe_write_csv(go_df, file.path(out_folder, paste0(out_prefix, "_", direction, "_GO_ORA_BH.csv")))
    safe_write_csv(kegg_df, file.path(out_folder, paste0(out_prefix, "_", direction, "_KEGG_ORA_BH.csv")))
    safe_write_csv(reactome_df, file.path(out_folder, paste0(out_prefix, "_", direction, "_Reactome_ORA_BH.csv")))
    list(go = go_df, kegg = kegg_df, reactome = reactome_df)
  }
  higher <- run_one(direction_sets$higher_in_AD, "higher_in_AD")
  lower <- run_one(direction_sets$lower_in_AD, "lower_in_AD")
  combined <- dplyr::bind_rows(
    dplyr::bind_rows(higher$go, lower$go) %>% dplyr::mutate(database = "GO"),
    dplyr::bind_rows(higher$kegg, lower$kegg) %>% dplyr::mutate(database = "KEGG"),
    dplyr::bind_rows(higher$reactome, lower$reactome) %>% dplyr::mutate(database = "Reactome")
  )
  if (nrow(combined) > 0) safe_write_csv(combined, file.path(out_folder, paste0(out_prefix, "_directional_ORA_GO_KEGG_Reactome_BH_combined.csv")))
  invisible(list(trace = trace_tbl, higher_in_AD = higher, lower_in_AD = lower, combined = combined))
}

run_clinical_severity_models <- function(dep_df, seq_cols, annot_tbl, out_folder) {
  candidate_outcomes <- c("cdr_boxscore", "cdr_global", "mmse_total", "udsfaq_total", "NPI", "Mini.SEA", "T.ADLQ")
  candidate_outcomes <- candidate_outcomes[candidate_outcomes %in% names(dep_df)]
  if (length(candidate_outcomes) == 0) return(NULL)
  results <- vector("list", length(candidate_outcomes)); names(results) <- candidate_outcomes
  for (outcome in candidate_outcomes) {
    md <- dep_df %>% dplyr::filter(!is.na(.data[[outcome]]))
    if (nrow(md) < 20) next
    fit_obj <- run_limma_dep_model(
      md, seq_cols, annot_tbl,
      paste0("~ ", outcome, " + SampleGroup + Age + Sex + Country + Education"),
      coef_name = outcome,
      model_name = paste0("clinical_severity_", outcome)
    )
    tt <- fit_obj$dep %>% dplyr::mutate(outcome = outcome) %>% dplyr::arrange(adj.P.Val)
    results[[outcome]] <- tt
    export_dep_table(tt, file.path(out_folder, paste0("severity_assoc_", outcome, ".csv")))
  }
  combined <- dplyr::bind_rows(results)
  if (nrow(combined) > 0) {
    combined <- combined %>% dplyr::group_by(outcome) %>% dplyr::mutate(adj.P.Val_within_outcome = p.adjust(P.Value, method = "BH")) %>% dplyr::ungroup()
    safe_write_csv(combined, file.path(out_folder, "severity_assoc_all_outcomes_combined.csv"))
  }
  combined
}

###############################################################################
# 01_import_annotation_qc
###############################################################################

message("01_import_annotation_qc")

meta_info_new <- read.csv(csv_file, check.names = FALSE, stringsAsFactors = FALSE)
required_columns(meta_info_new, c("SampleId", "Sex", "Age", "Country", "Education", "ApoE"), "metadata CSV")

meta_info_new <- meta_info_new %>%
  dplyr::mutate(
    SampleId = as.character(SampleId),
    Sex = dplyr::recode(as.character(Sex), "1" = "F", "2" = "M", .default = as.character(Sex)),
    Mini.SEA = pick_col(., c("Mini.SEA", "Mini-SEA")),
    T.ADLQ = pick_col(., c("T.ADLQ", "T-ADLQ")),
    p.tau217 = pick_col(., c("p.tau217", "p-tau217", "p_tau217")),
    p.tau181 = pick_col(., c("p.tau181", "p-tau181", "p_tau181")),
    ratio.AB42.40 = pick_col(., c("ratio.AB42.40", "ratio AB42/40", "ratio_AB42_40")),
    ApoE = norm_char(ApoE),
    APOE_group = dplyr::case_when(
      is.na(ApoE) | ApoE == "" ~ "Unknown",
      ApoE %in% c("e2/e4", "e3/e4", "e4/e4") ~ "E4 carrier",
      ApoE %in% c("e2/e2", "e2/e3", "e3/e3") ~ "Non-E4",
      TRUE ~ "Other"
    ),
    APOE4_carrier = dplyr::case_when(
      ApoE %in% c("e2/e4", "e3/e4", "e4/e4") ~ 1,
      ApoE %in% c("e2/e2", "e2/e3", "e3/e3") ~ 0,
      TRUE ~ NA_real_
    )
  )

my_adat <- SomaDataIO::read_adat(adat_file)
my_adat$SampleId <- as.character(my_adat$SampleId)
soma_info_internal <- extract_adat_annotation_from_table(adat_file)

sample_data <- meta_info_new %>% dplyr::left_join(my_adat, by = "SampleId", suffix = c(".csv", ".adat"))

if ("SampleType.csv" %in% names(sample_data) && "SampleType.adat" %in% names(sample_data)) {
  sample_data <- sample_data %>% dplyr::mutate(SampleType = dplyr::coalesce(`SampleType.adat`, `SampleType.csv`))
}
if ("SampleGroup.csv" %in% names(sample_data) && "SampleGroup.adat" %in% names(sample_data)) {
  sample_data <- sample_data %>% dplyr::mutate(SampleGroup = dplyr::coalesce(`SampleGroup.adat`, `SampleGroup.csv`))
}
if ("RowCheck.adat" %in% names(sample_data) && !"RowCheck" %in% names(sample_data)) {
  sample_data <- sample_data %>% dplyr::mutate(RowCheck = `RowCheck.adat`)
}
if ("PlateId.adat" %in% names(sample_data) && !"PlateId" %in% names(sample_data)) {
  sample_data <- sample_data %>% dplyr::mutate(PlateId = `PlateId.adat`)
}
if (!"SampleType" %in% names(sample_data)) stop("SampleType could not be recovered after CSV + ADAT merge.")
if (!"SampleGroup" %in% names(sample_data)) stop("SampleGroup could not be recovered after CSV + ADAT merge.")

sample_data_raw_merged <- sample_data
sample_data <- sample_data %>% dplyr::filter(is.na(RowCheck) | RowCheck == "PASS")

seq_cols <- grep("^seq[._]", names(sample_data), value = TRUE)
if (length(seq_cols) == 0) stop("No seq.* or seq_* columns found after merge.")

annot_tbl <- build_annotation_table(soma_info_internal, seq_cols)
protein_universe <- soma_info_internal %>%
  dplyr::mutate(AptName = as.character(AptName), Organism = as.character(Organism), Type = as.character(Type), EntrezGeneSymbol = as.character(EntrezGeneSymbol)) %>%
  dplyr::filter(Organism == "Human", Type == "Protein", !is.na(EntrezGeneSymbol), EntrezGeneSymbol != "") %>%
  dplyr::pull(AptName) %>%
  unique() %>%
  intersect(seq_cols)
if (length(protein_universe) == 0) stop("No proteins matched the internal ADAT protein universe filtering.")

protein_universe_audit <- tibble::tibble(
  metric = c("seq_columns_detected", "adat_annotation_rows", "human_protein_rows_with_gene_symbol", "protein_universe_intersecting_expression_columns", "annotation_rows_with_missing_gene_symbol"),
  n = c(
    length(seq_cols),
    nrow(soma_info_internal),
    soma_info_internal %>% dplyr::filter(Organism == "Human", Type == "Protein", !is.na(EntrezGeneSymbol), EntrezGeneSymbol != "") %>% nrow(),
    length(protein_universe),
    sum(is.na(annot_tbl$EntrezGeneSymbol) | annot_tbl$EntrezGeneSymbol == "")
  )
)

analysis_checkpoints <- add_checkpoint(analysis_checkpoints, "Raw metadata rows", nrow(meta_info_new))
analysis_checkpoints <- add_checkpoint(analysis_checkpoints, "Merged metadata + ADAT rows", nrow(sample_data_raw_merged))
analysis_checkpoints <- add_checkpoint(analysis_checkpoints, "Rows after RowCheck PASS or missing", nrow(sample_data))
analysis_checkpoints <- add_checkpoint(analysis_checkpoints, "Detected seq columns", length(seq_cols))
analysis_checkpoints <- add_checkpoint(analysis_checkpoints, "Internal ADAT protein universe", length(protein_universe))

safe_write_csv(soma_info_internal, file.path(outdir, "result", "01_qc", "annotations", "soma_info_internal_from_adat.csv"))
safe_write_csv(annot_tbl, file.path(outdir, "result", "01_qc", "annotations", "protein_annotation_dictionary.csv"))
safe_write_csv(tibble::tibble(AptName = protein_universe), file.path(outdir, "result", "01_qc", "sample_tracking", "protein_universe_internal_adat.csv"))
safe_write_csv(protein_universe_audit, file.path(outdir, "result", "01_qc", "sample_tracking", "protein_universe_audit.csv"))

sample_tracking_tbl <- tibble::tibble(
  step = c("Initial merged dataset", "RowCheck PASS or missing", "SampleType == Sample", "Non-missing Age and Sex", "Non-missing Country and Education", "Non-missing SampleGroup", "Restricted to CN/AD for DEP", "Restricted to non-missing APOE4 for APOE sensitivity"),
  n = c(
    nrow(sample_data_raw_merged),
    nrow(sample_data),
    nrow(sample_data %>% dplyr::filter(SampleType == "Sample")),
    nrow(sample_data %>% dplyr::filter(SampleType == "Sample", !is.na(Age), !is.na(Sex))),
    nrow(sample_data %>% dplyr::filter(SampleType == "Sample", !is.na(Age), !is.na(Sex), !is.na(Country), !is.na(Education))),
    nrow(sample_data %>% dplyr::filter(SampleType == "Sample", !is.na(Age), !is.na(Sex), !is.na(Country), !is.na(Education)) %>% tidyr::drop_na(SampleGroup)),
    nrow(sample_data %>% dplyr::filter(SampleType == "Sample", !is.na(Age), !is.na(Sex), !is.na(Country), !is.na(Education)) %>% tidyr::drop_na(SampleGroup) %>% dplyr::filter(SampleGroup %in% MAIN_GROUPS)),
    nrow(sample_data %>% dplyr::filter(SampleType == "Sample", !is.na(Age), !is.na(Sex), !is.na(Country), !is.na(Education)) %>% tidyr::drop_na(SampleGroup) %>% dplyr::filter(SampleGroup %in% MAIN_GROUPS, !is.na(APOE4_carrier)))
  )
)
safe_write_csv(sample_tracking_tbl, file.path(outdir, "result", "01_qc", "sample_tracking", "sample_tracking_workflow.csv"))

summary_tbl <- sample_data %>%
  dplyr::filter(!is.na(Age), !is.na(Sex)) %>%
  tidyr::drop_na(SampleGroup) %>%
  dplyr::select(-dplyr::any_of(seq_cols)) %>%
  dplyr::group_by(SampleGroup) %>%
  dplyr::group_modify(~ make_summary_table(.x)) %>%
  dplyr::ungroup()
na_count_tbl <- sample_data %>%
  dplyr::select(-dplyr::any_of(seq_cols)) %>%
  dplyr::summarise(dplyr::across(dplyr::everything(), ~ sum(is.na(.))))
safe_write_csv(summary_tbl, file.path(outdir, "result", "01_qc", "Table_1_sample_characteristics_by_group.csv"))
safe_write_csv(na_count_tbl, file.path(outdir, "result", "01_qc", "metadata_missingness_count.csv"))

###############################################################################
# 02_expression_matrices_pca
###############################################################################

message("02_expression_matrices_pca")

analysis_df <- sample_data %>%
  dplyr::filter(!is.na(Age), !is.na(Sex), !is.na(Country), !is.na(Education)) %>%
  dplyr::filter(SampleType == "Sample") %>%
  tidyr::drop_na(SampleGroup) %>%
  tibble::as_tibble()

protein_vars_present <- intersect(protein_universe, names(analysis_df))
if (length(protein_vars_present) == 0) stop("No protein_universe columns are present in analysis_df.")

normalized_expr_all <- analysis_df %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(protein_vars_present), safe_log2_vector)) %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(protein_vars_present), safe_standardize_vector))
normalized_expr_CN_AD <- normalized_expr_all %>% dplyr::filter(SampleGroup %in% MAIN_GROUPS)
normalized_expr <- normalized_expr_CN_AD
safe_write_csv(normalized_expr_all, file.path(outdir, "result", "02_pca", "normalized_log2_z_expression_all_samples.csv"))
safe_write_csv(normalized_expr_CN_AD, file.path(outdir, "result", "02_pca", "normalized_log2_z_expression_CN_AD.csv"))
analysis_checkpoints <- add_checkpoint(analysis_checkpoints, "Analysis dataset rows before CN/AD restriction", nrow(normalized_expr_all))
analysis_checkpoints <- add_checkpoint(analysis_checkpoints, "CN/AD normalized expression rows", nrow(normalized_expr_CN_AD))
analysis_checkpoints <- add_checkpoint(analysis_checkpoints, "Protein columns in normalized expression matrix", length(protein_vars_present))

if (RUN_FLAGS$run_pca) {
  pca_input <- normalized_expr_all %>% dplyr::select(dplyr::all_of(protein_vars_present)) %>% as.matrix()
  keep_cols <- apply(pca_input, 2, function(z) stats::sd(z, na.rm = TRUE) > 0)
  pca_input <- pca_input[, keep_cols, drop = FALSE]
  pca_input[!is.finite(pca_input)] <- NA_real_
  for (j in seq_len(ncol(pca_input))) {
    miss <- is.na(pca_input[, j])
    if (any(miss)) pca_input[miss, j] <- stats::median(pca_input[, j], na.rm = TRUE)
  }
  pca <- stats::prcomp(pca_input, center = TRUE, scale. = FALSE)
  pca_var <- as.data.frame(t(summary(pca)$importance))
  colnames(pca_var) <- c("stdv", "percent", "cumulative")
  pca_var$percent <- round(pca_var$percent * 100, 2)
  pca_var$cumulative <- round(pca_var$cumulative * 100, 2)
  pca_var$PC <- rownames(pca_var)
  pca_df <- as.data.frame(pca$x)[, seq_len(min(10, ncol(pca$x))), drop = FALSE] %>%
    dplyr::mutate(SampleId = normalized_expr_all$SampleId, .before = 1) %>%
    dplyr::bind_cols(normalized_expr_all %>% dplyr::select(-dplyr::any_of(c(protein_vars_present, "SampleId"))))
  safe_write_csv(pca_var, file.path(outdir, "result", "02_pca", "pca_variance_explained.csv"))
  safe_write_csv(pca_df, file.path(outdir, "result", "02_pca", "pca_scores_with_metadata.csv"))
}

###############################################################################
# 03_main_DEP_gene_collapsed
###############################################################################

message("03_main_DEP_gene_collapsed")

dep_df <- sample_data %>%
  dplyr::filter(!is.na(Age), !is.na(Sex), !is.na(Country), !is.na(Education)) %>%
  dplyr::filter(SampleType == "Sample") %>%
  tidyr::drop_na(SampleGroup) %>%
  dplyr::filter(SampleGroup %in% MAIN_GROUPS) %>%
  dplyr::mutate(SampleGroup = factor(SampleGroup, levels = MAIN_GROUPS), Sex = factor(Sex), Country = factor(Country)) %>%
  tibble::as_tibble()

dep_protein_cols <- intersect(protein_universe, names(dep_df))
if (length(dep_protein_cols) == 0) stop("No protein columns found for DEP.")
analysis_checkpoints <- add_checkpoint(analysis_checkpoints, "DEP input samples", nrow(dep_df))
analysis_checkpoints <- add_checkpoint(analysis_checkpoints, "DEP input proteins", length(dep_protein_cols))
safe_write_csv(tibble::tibble(SampleId = dep_df$SampleId), file.path(outdir, "result", "01_qc", "sample_tracking", "main_dep_sample_ids.csv"))

main_formula <- "~ SampleGroup + Age + Sex + Country + Education"
main_fit <- run_limma_dep_model(dep_df, dep_protein_cols, annot_tbl, main_formula, "SampleGroupAD", MAIN_FDR, model_name = "main_DEP")
model_design_diagnostics <- dplyr::bind_rows(model_design_diagnostics, main_fit$design_diagnostic)
model_formula_manifest <- add_model_formula(model_formula_manifest, "main_DEP", main_formula, "SampleGroupAD", nrow(main_fit$metadata), length(dep_protein_cols), "Primary DEP model; protein is the outcome; coefficients are adjusted log2 RFU differences.")
fit_main <- main_fit$fit
design_main <- main_fit$design
DEP_aptamer <- main_fit$dep %>% dplyr::filter(AptName %in% dep_protein_cols)
DEP <- DEP_aptamer
DEP_gene <- collapse_dep_to_gene(DEP_aptamer)
DEP_counts <- dplyr::bind_rows(
  summarize_dep_counts(DEP_aptamer, universe_label = "aptamer_level"),
  summarize_dep_counts(DEP_gene, universe_label = "gene_collapsed")
)
analysis_checkpoints <- add_checkpoint(analysis_checkpoints, "Gene-collapsed interpretation universe", nrow(DEP_gene))
analysis_checkpoints <- add_checkpoint(analysis_checkpoints, "Main DEP genes FDR < 0.05", sum(DEP_gene$adj.P.Val < MAIN_FDR, na.rm = TRUE))
analysis_checkpoints <- add_checkpoint(analysis_checkpoints, "Main DEP genes FDR < 0.01", sum(DEP_gene$adj.P.Val < STRICT_FDR, na.rm = TRUE))
export_dep_table(DEP_aptamer, file.path(outdir, "result", "03_dep", "aptamer_level", "AD_vs_CN_full_limma_results_aptamer_level.csv"))
export_dep_table(DEP_gene, file.path(outdir, "result", "03_dep", "gene_collapsed", "AD_vs_CN_full_limma_results_gene_collapsed.csv"))
safe_write_csv(DEP_counts, file.path(outdir, "result", "03_dep", "AD_vs_CN_DEP_counts_by_universe.csv"))

# Main DEP filtered tables for direct FDR 0.05 / FDR 0.01 reporting.
export_fdr_specific_dep_tables(
  DEP_gene,
  file.path(outdir, "result", "03_dep", "gene_collapsed", "FDR_specific"),
  "AD_vs_CN_DEP_gene_collapsed"
)
export_fdr_specific_dep_tables(
  DEP_aptamer,
  file.path(outdir, "result", "03_dep", "aptamer_level", "FDR_specific"),
  "AD_vs_CN_DEP_aptamer_level"
)

###############################################################################
# 04_sensitivity_APOE_CDRSB_vascular
###############################################################################

message("04_sensitivity_APOE_CDRSB_vascular")

DEP_APOE <- NULL; DEP_APOE_gene <- NULL; apoe_compare_tbl <- NULL
if (RUN_FLAGS$run_apoe_sensitivity) {
  dep_df_apoe <- dep_df %>% dplyr::filter(!is.na(APOE4_carrier)) %>% dplyr::mutate(APOE4_carrier = factor(APOE4_carrier))
  if (nrow(dep_df_apoe) >= 20 && length(unique(dep_df_apoe$APOE4_carrier)) >= 2) {
    apoe_formula <- "~ SampleGroup + Age + Sex + Country + Education + APOE4_carrier"
    apoe_fit <- run_limma_dep_model(dep_df_apoe, dep_protein_cols, annot_tbl, apoe_formula, "SampleGroupAD", MAIN_FDR, model_name = "APOE_adjusted_DEP")
    model_design_diagnostics <- dplyr::bind_rows(model_design_diagnostics, apoe_fit$design_diagnostic)
    model_formula_manifest <- add_model_formula(model_formula_manifest, "APOE_adjusted_DEP", apoe_formula, "SampleGroupAD", nrow(apoe_fit$metadata), length(dep_protein_cols), "Sensitivity model adding APOE epsilon 4 carrier status.")
    DEP_APOE <- apoe_fit$dep %>% dplyr::filter(AptName %in% dep_protein_cols)
    DEP_APOE_gene <- collapse_dep_to_gene(DEP_APOE)
    apoe_compare_tbl <- build_gene_compare_table(DEP_gene, DEP_APOE_gene, "APOE-adjusted") %>%
      dplyr::rename(logFC_apoe = logFC_secondary, adj.P.Val_apoe = adj.P.Val_secondary, P.Value_apoe = P.Value_secondary, type_apoe = type_secondary)
    export_dep_table(DEP_APOE, file.path(outdir, "result", "04_sensitivity", "apoe", "AD_vs_CN_APOE_adjusted_full_limma_results_aptamer_level.csv"))
    export_dep_table(DEP_APOE_gene, file.path(outdir, "result", "04_sensitivity", "apoe", "AD_vs_CN_APOE_adjusted_full_limma_results_gene_collapsed.csv"))
    safe_write_csv(apoe_compare_tbl, file.path(outdir, "result", "04_sensitivity", "apoe", "primary_vs_APOE_adjusted_gene_comparison.csv"))
    safe_write_csv(summarize_dep_counts(DEP_APOE_gene, universe_label = "gene_collapsed_APOE_adjusted"), file.path(outdir, "result", "04_sensitivity", "apoe", "APOE_adjusted_DEP_counts_gene_collapsed.csv"))
    export_fdr_specific_dep_tables(
      DEP_APOE_gene,
      file.path(outdir, "result", "04_sensitivity", "apoe", "FDR_specific"),
      "AD_vs_CN_APOE_adjusted_DEP_gene_collapsed"
    )
    analysis_checkpoints <- add_checkpoint(analysis_checkpoints, "APOE sensitivity samples", nrow(apoe_fit$metadata))
  } else {
    analysis_checkpoints <- add_checkpoint(analysis_checkpoints, "APOE sensitivity skipped", nrow(dep_df_apoe), note = "Insufficient APOE4 data or only one APOE4 class.")
  }
}

DEP_CDRSB <- NULL; DEP_CDRSB_gene <- NULL; cdrsb_compare_tbl <- NULL
DEP_CDRSB_ADJ <- NULL; DEP_CDRSB_ADJ_gene <- NULL; cdrsb_adjusted_compare_tbl <- NULL
cdrsb_col <- c("cdr_boxscore", "CDR_SB", "cdr_sb", "cdr_sum_box")[c("cdr_boxscore", "CDR_SB", "cdr_sb", "cdr_sum_box") %in% names(dep_df)][1]
if (RUN_FLAGS$run_cdrsb_sensitivity && !is.na(cdrsb_col)) {
  # Revised reviewer-facing analysis:
  # CDR-SB is almost structurally tied to diagnosis in a CN vs AD case-control design
  # because CN participants have CDR-SB values at or near zero. Therefore, CDR-SB
  # is not added to the diagnostic CN vs AD model as a conventional covariate here.
  # Instead, clinical severity is tested within AD participants only.
  dep_df_cdr <- dep_df %>%
    dplyr::filter(SampleGroup == "AD", !is.na(.data[[cdrsb_col]])) %>%
    dplyr::mutate(
      Sex = factor(Sex),
      Country = factor(Country)
    )

  cdrsb_has_variation <- nrow(dep_df_cdr) >= 20 &&
    stats::sd(suppressWarnings(as.numeric(dep_df_cdr[[cdrsb_col]])), na.rm = TRUE) > 0

  if (cdrsb_has_variation) {
    cdrsb_formula <- paste0("~ ", cdrsb_col, " + Age + Sex + Country + Education")
    cdrsb_fit <- run_limma_dep_model(
      dep_df_cdr,
      dep_protein_cols,
      annot_tbl,
      cdrsb_formula,
      coef_name = cdrsb_col,
      fdr = MAIN_FDR,
      model_name = "AD_only_CDRSB_severity"
    )
    model_design_diagnostics <- dplyr::bind_rows(model_design_diagnostics, cdrsb_fit$design_diagnostic)
    model_formula_manifest <- add_model_formula(
      model_formula_manifest,
      "AD_only_CDRSB_severity",
      cdrsb_formula,
      cdrsb_col,
      nrow(cdrsb_fit$metadata),
      length(dep_protein_cols),
      paste0("Within-AD clinical severity model. Protein is the outcome; coefficient is the association per one-unit increase in ", cdrsb_col, ", adjusted for age, sex, country and education. This avoids collinearity between diagnosis and CDR-SB in CN vs AD models.")
    )

    DEP_CDRSB <- cdrsb_fit$dep %>%
      dplyr::filter(AptName %in% dep_protein_cols) %>%
      dplyr::mutate(
        model_context = "AD_only_CDRSB_severity",
        coefficient_meaning = paste0("Association with ", cdrsb_col, " within AD participants; positive values indicate higher protein abundance with greater clinical severity."),
        severity_variable = cdrsb_col
      )
    DEP_CDRSB_gene <- collapse_dep_to_gene(DEP_CDRSB) %>%
      dplyr::mutate(
        model_context = "AD_only_CDRSB_severity",
        coefficient_meaning = paste0("Association with ", cdrsb_col, " within AD participants; positive values indicate higher protein abundance with greater clinical severity."),
        severity_variable = cdrsb_col
      )

    # Alignment table: this is not a diagnostic attenuation table. It compares
    # primary AD-vs-CN diagnostic effects with within-AD CDR-SB severity slopes
    # only to summarize whether diagnostic direction and severity direction align.
    cdrsb_compare_tbl <- DEP_gene %>%
      dplyr::select(EntrezGeneSymbol, AptName, SeqId, Protein_Name, logFC, adj.P.Val, P.Value, type) %>%
      dplyr::rename(
        AptName_primary = AptName,
        SeqId_primary = SeqId,
        Protein_Name_primary = Protein_Name,
        logFC_primary = logFC,
        adj.P.Val_primary = adj.P.Val,
        P.Value_primary = P.Value,
        type_primary = type
      ) %>%
      dplyr::inner_join(
        DEP_CDRSB_gene %>%
          dplyr::select(EntrezGeneSymbol, AptName, SeqId, Protein_Name, logFC, adj.P.Val, P.Value, type) %>%
          dplyr::rename(
            AptName_severity = AptName,
            SeqId_severity = SeqId,
            Protein_Name_severity = Protein_Name,
            beta_CDRSB_AD_only = logFC,
            adj.P.Val_CDRSB_AD_only = adj.P.Val,
            P.Value_CDRSB_AD_only = P.Value,
            type_CDRSB_AD_only = type
          ),
        by = "EntrezGeneSymbol"
      ) %>%
      dplyr::mutate(
        comparison = "Primary AD-vs-CN effect versus within-AD CDR-SB severity slope",
        Protein_Name = dplyr::coalesce(Protein_Name_primary, Protein_Name_severity),
        AptName = dplyr::coalesce(AptName_primary, AptName_severity),
        SeqId = dplyr::coalesce(SeqId_primary, SeqId_severity),
        same_direction = sign(logFC_primary) == sign(beta_CDRSB_AD_only),
        primary_fdr005 = adj.P.Val_primary < MAIN_FDR,
        severity_fdr005 = adj.P.Val_CDRSB_AD_only < MAIN_FDR,
        primary_fdr005_and_same_direction = primary_fdr005 & same_direction,
        primary_fdr005_and_severity_fdr005_same_direction = primary_fdr005 & severity_fdr005 & same_direction,
        note = "CDR-SB model is AD-only; beta_CDRSB_AD_only is not an adjusted AD-vs-CN logFC."
      )

    export_dep_table(DEP_CDRSB, file.path(outdir, "result", "04_sensitivity", "cdrsb", "AD_only", "AD_only_CDRSB_severity_full_limma_results_aptamer_level.csv"))
    export_dep_table(DEP_CDRSB_gene, file.path(outdir, "result", "04_sensitivity", "cdrsb", "AD_only", "AD_only_CDRSB_severity_full_limma_results_gene_collapsed.csv"))
    safe_write_csv(cdrsb_compare_tbl, file.path(outdir, "result", "04_sensitivity", "cdrsb", "AD_only", "primary_AD_vs_CN_vs_AD_only_CDRSB_severity_alignment.csv"))
    safe_write_csv(summarize_dep_counts(DEP_CDRSB_gene, universe_label = "gene_collapsed_AD_only_CDRSB_severity"), file.path(outdir, "result", "04_sensitivity", "cdrsb", "AD_only", "AD_only_CDRSB_severity_counts_gene_collapsed.csv"))
    export_fdr_specific_dep_tables(
      DEP_CDRSB_gene,
      file.path(outdir, "result", "04_sensitivity", "cdrsb", "AD_only", "FDR_specific"),
      "AD_only_CDRSB_severity_DEP_gene_collapsed"
    )


    analysis_checkpoints <- add_checkpoint(
      analysis_checkpoints,
      "AD-only CDR-SB severity samples",
      nrow(cdrsb_fit$metadata),
      note = paste0("CDR-SB column used: ", cdrsb_col, "; analysis restricted to AD participants to avoid diagnosis/CDR-SB collinearity.")
    )
  } else {
    analysis_checkpoints <- add_checkpoint(
      analysis_checkpoints,
      "AD-only CDR-SB severity skipped",
      nrow(dep_df_cdr),
      note = paste0("Insufficient AD participants with non-missing or variable ", cdrsb_col, ".")
    )
  }
} else if (RUN_FLAGS$run_cdrsb_sensitivity) {
  analysis_checkpoints <- add_checkpoint(analysis_checkpoints, "AD-only CDR-SB severity skipped", NA_integer_, note = "No CDR-SB column detected.")
}

# Secondary diagnostic attenuation model:
# This model intentionally reintroduces CDR-SB into the CN-vs-AD diagnostic model
# as a reviewer-facing sensitivity analysis. It asks how much of the AD-vs-CN
# proteomic contrast is attenuated after accounting for clinical severity. It is
# not interpreted as a within-AD severity model and does not replace the AD-only
# CDR-SB severity analysis above.
if (isTRUE(RUN_FLAGS$run_cdrsb_adjusted_diagnostic_sensitivity) && !is.na(cdrsb_col)) {
  dep_df_cdr_adj <- dep_df %>%
    dplyr::filter(!is.na(.data[[cdrsb_col]])) %>%
    dplyr::mutate(
      Sex = factor(Sex),
      Country = factor(Country),
      SampleGroup = factor(SampleGroup, levels = c("CN", "AD"))
    )

  cdrsb_adj_has_variation <- nrow(dep_df_cdr_adj) >= 20 &&
    length(unique(stats::na.omit(dep_df_cdr_adj$SampleGroup))) >= 2 &&
    stats::sd(suppressWarnings(as.numeric(dep_df_cdr_adj[[cdrsb_col]])), na.rm = TRUE) > 0

  if (cdrsb_adj_has_variation) {
    cdrsb_adj_formula <- paste0("~ SampleGroup + ", cdrsb_col, " + Age + Sex + Country + Education")
    cdrsb_adj_fit <- run_limma_dep_model(
      dep_df_cdr_adj,
      dep_protein_cols,
      annot_tbl,
      cdrsb_adj_formula,
      coef_name = "SampleGroupAD",
      fdr = MAIN_FDR,
      model_name = "CDRSB_adjusted_AD_vs_CN_diagnostic_sensitivity"
    )
    model_design_diagnostics <- dplyr::bind_rows(model_design_diagnostics, cdrsb_adj_fit$design_diagnostic)
    model_formula_manifest <- add_model_formula(
      model_formula_manifest,
      "CDRSB_adjusted_AD_vs_CN_diagnostic_sensitivity",
      cdrsb_adj_formula,
      "SampleGroupAD",
      nrow(cdrsb_adj_fit$metadata),
      length(dep_protein_cols),
      paste0("Secondary diagnostic attenuation model. Protein is the outcome; coefficient is AD versus CN adjusted for age, sex, country, education and ", cdrsb_col, ". Interpreted as diagnostic attenuation, not as within-AD severity biology.")
    )

    DEP_CDRSB_ADJ <- cdrsb_adj_fit$dep %>%
      dplyr::filter(AptName %in% dep_protein_cols) %>%
      dplyr::mutate(
        model_context = "CDRSB_adjusted_AD_vs_CN_diagnostic_sensitivity",
        coefficient_meaning = paste0("AD versus CN diagnostic log2FC after additional adjustment for ", cdrsb_col, "; interpreted as attenuation of the diagnostic contrast."),
        severity_variable = cdrsb_col
      )
    DEP_CDRSB_ADJ_gene <- collapse_dep_to_gene(DEP_CDRSB_ADJ) %>%
      dplyr::mutate(
        model_context = "CDRSB_adjusted_AD_vs_CN_diagnostic_sensitivity",
        coefficient_meaning = paste0("AD versus CN diagnostic log2FC after additional adjustment for ", cdrsb_col, "; interpreted as attenuation of the diagnostic contrast."),
        severity_variable = cdrsb_col
      )

    cdrsb_adjusted_compare_tbl <- build_gene_compare_table(
      DEP_gene,
      DEP_CDRSB_ADJ_gene,
      "CDR-SB-adjusted AD-vs-CN diagnostic sensitivity"
    ) %>%
      dplyr::rename(
        logFC_CDRSB_adjusted = logFC_secondary,
        adj.P.Val_CDRSB_adjusted = adj.P.Val_secondary,
        P.Value_CDRSB_adjusted = P.Value_secondary,
        type_CDRSB_adjusted = type_secondary
      ) %>%
      dplyr::mutate(
        comparison = "Primary AD-vs-CN effect versus CDR-SB-adjusted AD-vs-CN diagnostic effect",
        model_context = "Secondary diagnostic attenuation; not within-AD severity",
        same_direction_CDRSB_adjusted = sign(logFC_primary) == sign(logFC_CDRSB_adjusted),
        preserved_FDR005_CDRSB_adjusted = adj.P.Val_CDRSB_adjusted < MAIN_FDR,
        attenuation_absolute = abs(logFC_primary) - abs(logFC_CDRSB_adjusted),
        attenuation_ratio = dplyr::if_else(abs(logFC_primary) > 0, abs(logFC_CDRSB_adjusted) / abs(logFC_primary), NA_real_),
        note = paste0("CDR-SB-adjusted model includes CN and AD and adjusts the diagnostic coefficient for ", cdrsb_col, "; use only as supplementary attenuation analysis.")
      )

    export_dep_table(
      DEP_CDRSB_ADJ,
      file.path(outdir, "result", "04_sensitivity", "cdrsb", "AD_vs_CN_adjusted", "AD_vs_CN_CDRSB_adjusted_full_limma_results_aptamer_level.csv")
    )
    export_dep_table(
      DEP_CDRSB_ADJ_gene,
      file.path(outdir, "result", "04_sensitivity", "cdrsb", "AD_vs_CN_adjusted", "AD_vs_CN_CDRSB_adjusted_full_limma_results_gene_collapsed.csv")
    )
    safe_write_csv(
      cdrsb_adjusted_compare_tbl,
      file.path(outdir, "result", "04_sensitivity", "cdrsb", "AD_vs_CN_adjusted", "primary_vs_CDRSB_adjusted_AD_vs_CN_gene_comparison.csv")
    )
    safe_write_csv(
      summarize_dep_counts(DEP_CDRSB_ADJ_gene, universe_label = "gene_collapsed_CDRSB_adjusted_AD_vs_CN"),
      file.path(outdir, "result", "04_sensitivity", "cdrsb", "AD_vs_CN_adjusted", "CDRSB_adjusted_AD_vs_CN_DEP_counts_gene_collapsed.csv")
    )
    export_fdr_specific_dep_tables(
      DEP_CDRSB_ADJ_gene,
      file.path(outdir, "result", "04_sensitivity", "cdrsb", "AD_vs_CN_adjusted", "FDR_specific"),
      "AD_vs_CN_CDRSB_adjusted_DEP_gene_collapsed"
    )

    analysis_checkpoints <- add_checkpoint(
      analysis_checkpoints,
      "CDR-SB-adjusted AD-vs-CN sensitivity samples",
      nrow(cdrsb_adj_fit$metadata),
      note = paste0("CDR-SB column used: ", cdrsb_col, "; secondary diagnostic attenuation model, not within-AD severity.")
    )
  } else {
    analysis_checkpoints <- add_checkpoint(
      analysis_checkpoints,
      "CDR-SB-adjusted AD-vs-CN sensitivity skipped",
      nrow(dep_df_cdr_adj),
      note = paste0("Insufficient CN/AD participants with non-missing or variable ", cdrsb_col, ".")
    )
  }
} else if (isTRUE(RUN_FLAGS$run_cdrsb_adjusted_diagnostic_sensitivity)) {
  analysis_checkpoints <- add_checkpoint(analysis_checkpoints, "CDR-SB-adjusted AD-vs-CN sensitivity skipped", NA_integer_, note = "No CDR-SB column detected.")
}

DEP_VASCULAR <- NULL; DEP_VASCULAR_gene <- NULL; vascular_compare_tbl <- NULL; vascular_covariates <- character(0)
if (RUN_FLAGS$run_vascular_metabolic_sensitivity) {
  vascular_covariates <- detect_vascular_metabolic_covariates(dep_df)
  vascular_covariates <- vascular_covariates[vapply(vascular_covariates, function(v) {
    x <- dep_df[[v]]
    if (is.numeric(x)) return(sum(!is.na(x)) >= 20 && stats::sd(x, na.rm = TRUE) > 0)
    length(unique(stats::na.omit(as.character(x)))) >= 2
  }, logical(1))]
  if (length(vascular_covariates) > 0) {
    dep_df_vascular <- dep_df %>%
      dplyr::filter(dplyr::if_all(dplyr::all_of(vascular_covariates), ~ !is.na(.x))) %>%
      prepare_optional_model_covariates(vascular_covariates)
    if (nrow(dep_df_vascular) >= 20) {
      vascular_formula <- paste0("~ SampleGroup + Age + Sex + Country + Education + ", paste(vascular_covariates, collapse = " + "))
      vascular_fit <- run_limma_dep_model(dep_df_vascular, dep_protein_cols, annot_tbl, vascular_formula, "SampleGroupAD", MAIN_FDR, model_name = "vascular_metabolic_adjusted_DEP")
      model_design_diagnostics <- dplyr::bind_rows(model_design_diagnostics, vascular_fit$design_diagnostic)
      model_formula_manifest <- add_model_formula(model_formula_manifest, "vascular_metabolic_adjusted_DEP", vascular_formula, "SampleGroupAD", nrow(vascular_fit$metadata), length(dep_protein_cols), paste("Optional sensitivity model; detected covariates:", paste(vascular_covariates, collapse = ", ")))
      DEP_VASCULAR <- vascular_fit$dep %>% dplyr::filter(AptName %in% dep_protein_cols)
      DEP_VASCULAR_gene <- collapse_dep_to_gene(DEP_VASCULAR)
      vascular_compare_tbl <- build_gene_compare_table(DEP_gene, DEP_VASCULAR_gene, "vascular-metabolic-adjusted") %>%
        dplyr::rename(logFC_vascular = logFC_secondary, adj.P.Val_vascular = adj.P.Val_secondary, P.Value_vascular = P.Value_secondary, type_vascular = type_secondary)
      export_dep_table(DEP_VASCULAR, file.path(outdir, "result", "04_sensitivity", "vascular_metabolic", "AD_vs_CN_vascular_metabolic_adjusted_full_limma_results_aptamer_level.csv"))
      export_dep_table(DEP_VASCULAR_gene, file.path(outdir, "result", "04_sensitivity", "vascular_metabolic", "AD_vs_CN_vascular_metabolic_adjusted_full_limma_results_gene_collapsed.csv"))
      safe_write_csv(vascular_compare_tbl, file.path(outdir, "result", "04_sensitivity", "vascular_metabolic", "main_vs_vascular_metabolic_adjusted_gene_comparison.csv"))
      safe_write_csv(summarize_dep_counts(DEP_VASCULAR_gene, universe_label = "gene_collapsed_vascular_metabolic_adjusted"), file.path(outdir, "result", "04_sensitivity", "vascular_metabolic", "vascular_metabolic_adjusted_DEP_counts_gene_collapsed.csv"))
      export_fdr_specific_dep_tables(
        DEP_VASCULAR_gene,
        file.path(outdir, "result", "04_sensitivity", "vascular_metabolic", "FDR_specific"),
        "AD_vs_CN_vascular_metabolic_adjusted_DEP_gene_collapsed"
      )
      safe_write_csv(tibble::tibble(covariate = vascular_covariates, n_complete_in_model = nrow(vascular_fit$metadata)), file.path(outdir, "result", "04_sensitivity", "vascular_metabolic", "vascular_metabolic_covariates_used.csv"))
      analysis_checkpoints <- add_checkpoint(analysis_checkpoints, "Vascular/metabolic sensitivity samples", nrow(vascular_fit$metadata), note = paste("Covariates:", paste(vascular_covariates, collapse = ", ")))
    }
  } else {
    analysis_checkpoints <- add_checkpoint(analysis_checkpoints, "Vascular/metabolic sensitivity skipped", NA_integer_, note = "No recognized vascular/metabolic covariates were available.")
  }
}


DEP_ATN <- NULL; DEP_ATN_gene <- NULL; atn_compare_tbl <- NULL; atn_covariates <- character(0); preferred_atn_covariates <- character(0)
if (RUN_FLAGS$run_atn_adjusted_sensitivity) {
  atn_covariates <- detect_atn_covariates(dep_df)

  # Primary AT(N)-adjusted sensitivity uses p-tau217, NfL and Aβ42/40 when available.
  # p-tau181 is detected but not added by default to avoid over-adjustment/collinearity
  # with p-tau217; it can be added manually if a separate sensitivity is required.
  preferred_atn_covariates <- intersect(c("p.tau217", "p_tau217", "p-tau217", "ptau217",
                                         "NfL", "NFL", "nfl", "Neurofilament_light",
                                         "ratio.AB42.40", "ratio_AB42_40", "ratio AB42/40", "AB42_40_ratio", "Aβ42/40"), atn_covariates)
  preferred_atn_covariates <- unique(preferred_atn_covariates)

  if (length(preferred_atn_covariates) >= 2) {
    dep_df_atn <- dep_df %>%
      dplyr::mutate(dplyr::across(dplyr::all_of(preferred_atn_covariates), ~ suppressWarnings(as.numeric(.x)))) %>%
      dplyr::filter(dplyr::if_all(dplyr::all_of(preferred_atn_covariates), ~ !is.na(.x)))

    if (nrow(dep_df_atn) >= 20) {
      atn_formula <- paste0("~ SampleGroup + Age + Sex + Country + Education + ", paste(preferred_atn_covariates, collapse = " + "))
      atn_fit <- run_limma_dep_model(
        dep_df_atn,
        dep_protein_cols,
        annot_tbl,
        atn_formula,
        "SampleGroupAD",
        MAIN_FDR,
        model_name = "ATN_adjusted_DEP"
      )
      model_design_diagnostics <- dplyr::bind_rows(model_design_diagnostics, atn_fit$design_diagnostic)
      model_formula_manifest <- add_model_formula(
        model_formula_manifest,
        "ATN_adjusted_DEP",
        atn_formula,
        "SampleGroupAD",
        nrow(atn_fit$metadata),
        length(dep_protein_cols),
        paste("Sensitivity model adjusting for plasma AT(N)-related biomarkers:", paste(preferred_atn_covariates, collapse = ", "))
      )

      DEP_ATN <- atn_fit$dep %>% dplyr::filter(AptName %in% dep_protein_cols)
      DEP_ATN_gene <- collapse_dep_to_gene(DEP_ATN)
      atn_compare_tbl <- build_gene_compare_table(DEP_gene, DEP_ATN_gene, "ATN-adjusted") %>%
        dplyr::rename(
          logFC_atn = logFC_secondary,
          adj.P.Val_atn = adj.P.Val_secondary,
          P.Value_atn = P.Value_secondary,
          type_atn = type_secondary
        ) %>%
        dplyr::mutate(
          primary_fdr005 = adj.P.Val_primary < MAIN_FDR,
          atn_fdr005 = adj.P.Val_atn < MAIN_FDR,
          primary_fdr005_preserved_same_direction = primary_fdr005 & atn_fdr005 & same_direction,
          note = "AT(N)-adjusted model tests whether the AD-vs-CN diagnostic coefficient is preserved after adjustment for plasma p-tau217, NfL and Aβ42/40 when available."
        )

      export_dep_table(DEP_ATN, file.path(outdir, "result", "04_sensitivity", "atn_adjusted", "AD_vs_CN_ATN_adjusted_full_limma_results_aptamer_level.csv"))
      export_dep_table(DEP_ATN_gene, file.path(outdir, "result", "04_sensitivity", "atn_adjusted", "AD_vs_CN_ATN_adjusted_full_limma_results_gene_collapsed.csv"))
      safe_write_csv(atn_compare_tbl, file.path(outdir, "result", "04_sensitivity", "atn_adjusted", "primary_vs_ATN_adjusted_gene_comparison.csv"))
      safe_write_csv(summarize_dep_counts(DEP_ATN_gene, universe_label = "gene_collapsed_ATN_adjusted"), file.path(outdir, "result", "04_sensitivity", "atn_adjusted", "ATN_adjusted_DEP_counts_gene_collapsed.csv"))
      safe_write_csv(tibble::tibble(covariate = preferred_atn_covariates, n_complete_in_model = nrow(atn_fit$metadata)), file.path(outdir, "result", "04_sensitivity", "atn_adjusted", "ATN_covariates_used.csv"))
      export_fdr_specific_dep_tables(
        DEP_ATN_gene,
        file.path(outdir, "result", "04_sensitivity", "atn_adjusted", "FDR_specific"),
        "AD_vs_CN_ATN_adjusted_DEP_gene_collapsed"
      )
      analysis_checkpoints <- add_checkpoint(
        analysis_checkpoints,
        "AT(N)-adjusted sensitivity samples",
        nrow(atn_fit$metadata),
        note = paste("Covariates:", paste(preferred_atn_covariates, collapse = ", "))
      )
    } else {
      analysis_checkpoints <- add_checkpoint(
        analysis_checkpoints,
        "AT(N)-adjusted sensitivity skipped",
        nrow(dep_df_atn),
        note = "Insufficient complete cases for detected AT(N) covariates."
      )
    }
  } else {
    analysis_checkpoints <- add_checkpoint(
      analysis_checkpoints,
      "AT(N)-adjusted sensitivity skipped",
      NA_integer_,
      note = paste0("Fewer than two usable AT(N) covariates detected. Detected: ", paste(atn_covariates, collapse = ", "))
    )
  }
}

###############################################################################
# 05_secondary_clinical_severity
###############################################################################

message("05_secondary_clinical_severity")
severity_results <- NULL
if (RUN_FLAGS$run_secondary_clinical_severity) {
  severity_results <- run_clinical_severity_models(dep_df, dep_protein_cols, annot_tbl, file.path(outdir, "result", "04_sensitivity", "clinical_severity"))
  if (!is.null(severity_results) && nrow(severity_results) > 0) {
    analysis_checkpoints <- add_checkpoint(analysis_checkpoints, "Secondary clinical severity associations rows", nrow(severity_results))
  } else {
    analysis_checkpoints <- add_checkpoint(analysis_checkpoints, "Secondary clinical severity associations skipped or empty")
  }
}

###############################################################################
# 06_enrichment_corrected
###############################################################################

message("06_enrichment_corrected")
enrichment_inventory <- tibble::tibble(analysis = character(), status = character(), file_prefix = character(), note = character())
if (RUN_FLAGS$run_corrected_enrichment) {
  gsea_go_kegg_main <- run_gsea_go_kegg(DEP_gene, "main_dep", file.path(outdir, "result", "05_enrichment_corrected", "gsea"))
  gsea_reactome_main <- run_gsea_reactome(DEP_gene, "main_dep", file.path(outdir, "result", "05_enrichment_corrected", "gsea"))
  gsea_hallmark_main <- run_gsea_hallmark(DEP_gene, "main_dep", file.path(outdir, "result", "05_enrichment_corrected", "gsea"))
  ora_directional_main <- run_corrected_ora_by_direction(DEP_gene, "main_dep", file.path(outdir, "result", "05_enrichment_corrected", "ora"), fdr = MAIN_FDR)
  enrichment_inventory <- dplyr::bind_rows(
    enrichment_inventory,
    tibble::tibble(analysis = "GSEA GO/KEGG", status = "completed", file_prefix = "main_dep", note = "BH-corrected"),
    tibble::tibble(analysis = "GSEA Reactome", status = "completed", file_prefix = "main_dep", note = "BH-corrected"),
    tibble::tibble(analysis = "GSEA Hallmark", status = "completed", file_prefix = "main_dep", note = "BH-corrected"),
    tibble::tibble(analysis = "Directional ORA GO/KEGG/Reactome", status = "completed", file_prefix = "main_dep", note = "BH-corrected; higher and lower in AD separated")
  )
  if (!is.null(DEP_APOE_gene) && nrow(DEP_APOE_gene) > 0) {
    gsea_reactome_apoe <- run_gsea_reactome(DEP_APOE_gene, "apoe_adjusted_dep", file.path(outdir, "result", "05_enrichment_corrected", "gsea"))
    enrichment_inventory <- dplyr::bind_rows(enrichment_inventory, tibble::tibble(analysis = "APOE-adjusted GSEA Reactome", status = "completed", file_prefix = "apoe_adjusted_dep", note = "BH-corrected"))
  }
  if (!is.null(DEP_ATN_gene) && nrow(DEP_ATN_gene) > 0) {
    gsea_reactome_atn <- run_gsea_reactome(DEP_ATN_gene, "atn_adjusted_dep", file.path(outdir, "result", "05_enrichment_corrected", "gsea"))
    enrichment_inventory <- dplyr::bind_rows(enrichment_inventory, tibble::tibble(analysis = "AT(N)-adjusted GSEA Reactome", status = "completed", file_prefix = "atn_adjusted_dep", note = "BH-corrected; diagnostic coefficient after plasma AT(N) adjustment"))
  }
  if (!is.null(DEP_CDRSB_gene) && nrow(DEP_CDRSB_gene) > 0) {
    gsea_reactome_cdrsb <- run_gsea_reactome(DEP_CDRSB_gene, "ad_only_cdrsb_severity", file.path(outdir, "result", "05_enrichment_corrected", "gsea"))
    enrichment_inventory <- dplyr::bind_rows(enrichment_inventory, tibble::tibble(analysis = "AD-only CDR-SB severity GSEA Reactome", status = "completed", file_prefix = "ad_only_cdrsb_severity", note = "BH-corrected; interpreted as within-AD severity association"))
  }
  if (!is.null(DEP_CDRSB_ADJ_gene) && nrow(DEP_CDRSB_ADJ_gene) > 0) {
    gsea_reactome_cdrsb_adjusted <- run_gsea_reactome(DEP_CDRSB_ADJ_gene, "cdrsb_adjusted_ad_vs_cn", file.path(outdir, "result", "05_enrichment_corrected", "gsea"))
    enrichment_inventory <- dplyr::bind_rows(enrichment_inventory, tibble::tibble(analysis = "CDR-SB-adjusted AD-vs-CN GSEA Reactome", status = "completed", file_prefix = "cdrsb_adjusted_ad_vs_cn", note = "BH-corrected; secondary diagnostic attenuation model, not within-AD severity"))
  }
}
safe_write_csv(enrichment_inventory, file.path(outdir, "result", "05_enrichment_corrected", "enrichment_inventory.csv"))

###############################################################################
# 07_country_site_robustness
###############################################################################

message("07_country_site_robustness")

country_counts <- dep_df %>%
  dplyr::count(Country, SampleGroup, name = "n") %>%
  tidyr::pivot_wider(names_from = SampleGroup, values_from = n, values_fill = 0) %>%
  ensure_count_cols(c("CN", "AD")) %>%
  dplyr::mutate(total = CN + AD) %>%
  dplyr::arrange(dplyr::desc(total))
safe_write_csv(country_counts, file.path(outdir, "result", "06_robustness", "country_loco", "tables", "country_group_counts.csv"))
countries_to_test <- country_counts %>% dplyr::filter(CN >= MIN_N_PER_GROUP_COUNTRY, AD >= MIN_N_PER_GROUP_COUNTRY) %>% dplyr::pull(Country)
analysis_checkpoints <- add_checkpoint(analysis_checkpoints, "Countries eligible for country robustness", length(countries_to_test), note = paste(countries_to_test, collapse = ", "))

loco_tables <- list(); loco_summary_metrics <- NULL; main_vs_loco_mean <- NULL
if (RUN_FLAGS$run_country_loco && length(countries_to_test) >= 2) {
  for (ctry in countries_to_test) {
    message("Running LOCO excluding: ", ctry)
    dat_loco <- dep_df %>% dplyr::filter(Country != ctry) %>% dplyr::mutate(Country = droplevels(factor(Country)), SampleGroup = factor(SampleGroup, levels = MAIN_GROUPS))
    formula_loco <- if (dplyr::n_distinct(dat_loco$Country) >= 2) "~ SampleGroup + Age + Sex + Country + Education" else "~ SampleGroup + Age + Sex + Education"
    loco_fit_obj <- run_limma_dep_model(dat_loco, dep_protein_cols, annot_tbl, formula_loco, "SampleGroupAD", MAIN_FDR, model_name = paste0("LOCO_excluding_", safe_file_tag(ctry)))
    model_design_diagnostics <- dplyr::bind_rows(model_design_diagnostics, loco_fit_obj$design_diagnostic)
    model_formula_manifest <- add_model_formula(model_formula_manifest, paste0("LOCO_excluding_", ctry), formula_loco, "SampleGroupAD", nrow(loco_fit_obj$metadata), length(dep_protein_cols), "Leave-one-country-out internal stability model.")
    dep_loco <- loco_fit_obj$dep %>% dplyr::mutate(excluded_country = ctry)
    export_dep_table(dep_loco, file.path(outdir, "result", "06_robustness", "country_loco", "tables", paste0("LOCO_excluding_", safe_file_tag(ctry), "_dep_results.csv")))
    cmp <- DEP %>%
      dplyr::select(AptName, Protein_Name, logFC, adj.P.Val, type) %>%
      dplyr::rename(Protein_Name_main = Protein_Name, main_logFC = logFC, main_adj.P.Val = adj.P.Val, main_type = type) %>%
      dplyr::inner_join(dep_loco %>% dplyr::select(AptName, Protein_Name, logFC, adj.P.Val, type, excluded_country) %>% dplyr::rename(Protein_Name_loco = Protein_Name, loco_logFC = logFC, loco_adj.P.Val = adj.P.Val, loco_type = type), by = "AptName") %>%
      dplyr::mutate(
        Protein_Name = dplyr::coalesce(Protein_Name_main, Protein_Name_loco),
        same_direction = sign(main_logFC) == sign(loco_logFC),
        main_sig_fdr005 = main_adj.P.Val < MAIN_FDR,
        loco_sig_fdr005 = loco_adj.P.Val < MAIN_FDR,
        main_sig_preserved = main_sig_fdr005 & loco_sig_fdr005 & same_direction,
        delta_logFC = loco_logFC - main_logFC
      )
    safe_write_csv(cmp, file.path(outdir, "result", "06_robustness", "country_loco", "tables", paste0("LOCO_excluding_", safe_file_tag(ctry), "_comparison_to_main.csv")))
    loco_tables[[ctry]] <- cmp
    loco_summary_metrics <- dplyr::bind_rows(loco_summary_metrics, tibble::tibble(
      excluded_country = ctry,
      n_samples = nrow(loco_fit_obj$metadata),
      logFC_correlation = suppressWarnings(cor(cmp$main_logFC, cmp$loco_logFC, use = "complete.obs")),
      slope = suppressWarnings(as.numeric(stats::coef(stats::lm(loco_logFC ~ main_logFC, data = cmp))[2])),
      direction_consistency_all = mean(cmp$same_direction, na.rm = TRUE),
      n_main_sig_fdr005 = sum(cmp$main_sig_fdr005, na.rm = TRUE),
      prop_main_sig_preserved = mean(cmp$main_sig_preserved[cmp$main_sig_fdr005], na.rm = TRUE)
    ))
  }
  safe_write_csv(loco_summary_metrics, file.path(outdir, "result", "06_robustness", "country_loco", "tables", "LOCO_summary_metrics.csv"))
  main_vs_loco_mean <- dplyr::bind_rows(loco_tables) %>%
    dplyr::group_by(AptName) %>%
    dplyr::summarise(
      Protein_Name = dplyr::first(Protein_Name),
      main_logFC = dplyr::first(main_logFC),
      main_adj.P.Val = dplyr::first(main_adj.P.Val),
      mean_loco_logFC = mean(loco_logFC, na.rm = TRUE),
      sd_loco_logFC = stats::sd(loco_logFC, na.rm = TRUE),
      n_loco_models = dplyr::n(),
      prop_same_direction = mean(same_direction, na.rm = TRUE),
      prop_main_sig_preserved = mean(main_sig_preserved[main_sig_fdr005], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(loco_delta = mean_loco_logFC - main_logFC)
  safe_write_csv(main_vs_loco_mean, file.path(outdir, "result", "06_robustness", "country_loco", "tables", "main_vs_meanLOCO_table.csv"))
}

country_specific_tbl <- NULL; meta_tbl <- NULL
if (RUN_FLAGS$run_country_meta && length(countries_to_test) >= 2) {
  country_specific_list <- list()
  for (ctry in countries_to_test) {
    message("Running country-specific DEP: ", ctry)
    dat_country <- dep_df %>% dplyr::filter(Country == ctry) %>% dplyr::mutate(SampleGroup = factor(SampleGroup, levels = MAIN_GROUPS), Sex = droplevels(factor(Sex)))
    if (nrow(dat_country) < 20) next
    country_formula <- "~ SampleGroup + Age + Sex + Education"
    country_fit_obj <- run_limma_dep_model(dat_country, dep_protein_cols, annot_tbl, country_formula, "SampleGroupAD", MAIN_FDR, model_name = paste0("country_specific_", safe_file_tag(ctry)))
    model_design_diagnostics <- dplyr::bind_rows(model_design_diagnostics, country_fit_obj$design_diagnostic)
    model_formula_manifest <- add_model_formula(model_formula_manifest, paste0("country_specific_", ctry), country_formula, "SampleGroupAD", nrow(country_fit_obj$metadata), length(dep_protein_cols), "Country-specific model for random-effects meta-analysis.")
    country_specific_list[[ctry]] <- country_fit_obj$dep %>%
      dplyr::mutate(Country = ctry) %>%
      dplyr::select(Country, AptName, Protein_Name, EntrezGeneSymbol, logFC, se, P.Value, adj.P.Val, type)
  }
  country_specific_tbl <- dplyr::bind_rows(country_specific_list)
  safe_write_csv(country_specific_tbl, file.path(outdir, "result", "06_robustness", "country_meta", "tables", "country_specific_DEP_results.csv"))
  if (nrow(country_specific_tbl) > 0) {
    meta_tbl <- country_specific_tbl %>%
      dplyr::filter(!is.na(se), is.finite(se), se > 0) %>%
      dplyr::group_by(AptName) %>%
      dplyr::group_modify(function(.x, .y) {
        if (nrow(.x) < 2) return(tibble::tibble(meta_logFC = NA_real_, meta_se = NA_real_, meta_p = NA_real_, I2 = NA_real_, tau2 = NA_real_, k = nrow(.x)))
        fit <- tryCatch(metafor::rma.uni(yi = .x$logFC, sei = .x$se, method = "REML"), error = function(e) NULL)
        if (is.null(fit)) return(tibble::tibble(meta_logFC = NA_real_, meta_se = NA_real_, meta_p = NA_real_, I2 = NA_real_, tau2 = NA_real_, k = nrow(.x)))
        tibble::tibble(meta_logFC = as.numeric(fit$b[1]), meta_se = fit$se, meta_p = fit$pval, I2 = fit$I2, tau2 = fit$tau2, k = nrow(.x))
      }) %>%
      dplyr::ungroup() %>%
      dplyr::left_join(DEP %>% dplyr::select(AptName, Protein_Name, EntrezGeneSymbol, logFC, adj.P.Val), by = "AptName") %>%
      dplyr::mutate(
        meta_adj.P.Val = p.adjust(meta_p, method = "BH"),
        meta_type = dplyr::case_when(meta_logFC > 0 & meta_adj.P.Val < MAIN_FDR ~ "Up", meta_logFC < 0 & meta_adj.P.Val < MAIN_FDR ~ "Down", TRUE ~ "NS")
      ) %>%
      dplyr::arrange(meta_adj.P.Val)
    safe_write_csv(meta_tbl, file.path(outdir, "result", "06_robustness", "country_meta", "tables", "country_meta_analysis_results.csv"))
  }
}

country_interaction_tbl <- NULL
if (RUN_FLAGS$run_country_interaction && length(countries_to_test) >= 2) {
  dat_interaction <- dep_df %>% dplyr::filter(Country %in% countries_to_test) %>% dplyr::mutate(Country = droplevels(factor(Country)), SampleGroup = factor(SampleGroup, levels = MAIN_GROUPS))
  expr <- dat_interaction %>% dplyr::select(dplyr::all_of(dep_protein_cols)) %>% as.matrix() %>% safe_log2_matrix() %>% t()
  design_interaction <- model.matrix(~ SampleGroup * Country + Age + Sex + Education, data = dat_interaction)
  model_design_diagnostics <- dplyr::bind_rows(model_design_diagnostics, validate_design_matrix(design_interaction, "SampleGroup_by_Country_interaction"))
  fit_interaction <- limma::eBayes(limma::lmFit(expr, design_interaction))
  interaction_coef <- grep("^SampleGroupAD:Country", colnames(design_interaction), value = TRUE)
  if (length(interaction_coef) > 0) {
    interaction_results <- purrr::map_dfr(interaction_coef, function(coef_name) {
      limma::topTable(fit_interaction, coef = coef_name, adjust.method = "BH", number = Inf) %>%
        tibble::rownames_to_column("feature_id_raw") %>%
        dplyr::mutate(AptName = feature_id_raw, interaction_term = coef_name) %>%
        dplyr::left_join(annot_tbl, by = "AptName") %>%
        dplyr::mutate(Protein_Name = make_protein_name(EntrezGeneSymbol, TargetFullName, Target, AptName, feature_id_raw))
    })
    country_interaction_tbl <- interaction_results %>%
      dplyr::group_by(AptName) %>%
      dplyr::summarise(Protein_Name = dplyr::first(Protein_Name), EntrezGeneSymbol = dplyr::first(EntrezGeneSymbol), min_interaction_p = min(P.Value, na.rm = TRUE), min_interaction_adjP = min(adj.P.Val, na.rm = TRUE), n_interaction_terms = dplyr::n(), .groups = "drop") %>%
      dplyr::arrange(min_interaction_adjP)
    safe_write_csv(interaction_results, file.path(outdir, "result", "06_robustness", "country_interaction", "tables", "country_interaction_all_terms.csv"))
    safe_write_csv(country_interaction_tbl, file.path(outdir, "result", "06_robustness", "country_interaction", "tables", "country_interaction_summary_by_protein.csv"))
  }
}

site_var <- c("Site", "site", "Center", "center", "Cohort", "cohort")[c("Site", "site", "Center", "center", "Cohort", "cohort") %in% names(dep_df)][1]
loso_summary_metrics <- NULL
if (RUN_FLAGS$run_site_robustness && !is.na(site_var)) {
  site_counts <- dep_df %>%
    dplyr::count(.data[[site_var]], SampleGroup, name = "n") %>%
    dplyr::rename(site = !!rlang::sym(site_var)) %>%
    tidyr::pivot_wider(names_from = SampleGroup, values_from = n, values_fill = 0) %>%
    ensure_count_cols(c("CN", "AD")) %>%
    dplyr::mutate(total = CN + AD)
  safe_write_csv(site_counts, file.path(outdir, "result", "06_robustness", "site_robustness", "site_group_counts.csv"))
  sites_to_test <- site_counts %>% dplyr::filter(CN >= MIN_N_PER_SITE_GROUP, AD >= MIN_N_PER_SITE_GROUP) %>% dplyr::pull(site)
  if (length(sites_to_test) >= 2) {
    for (site_i in sites_to_test) {
      message("Running LOSO excluding: ", site_i)
      dat_loso <- dep_df %>% dplyr::filter(.data[[site_var]] != site_i) %>% dplyr::mutate(Country = droplevels(factor(Country)), SampleGroup = factor(SampleGroup, levels = MAIN_GROUPS))
      loso_formula <- "~ SampleGroup + Age + Sex + Country + Education"
      loso_fit_obj <- run_limma_dep_model(dat_loso, dep_protein_cols, annot_tbl, loso_formula, "SampleGroupAD", MAIN_FDR, model_name = paste0("LOSO_excluding_", safe_file_tag(site_i)))
      model_design_diagnostics <- dplyr::bind_rows(model_design_diagnostics, loso_fit_obj$design_diagnostic)
      model_formula_manifest <- add_model_formula(model_formula_manifest, paste0("LOSO_excluding_", site_i), loso_formula, "SampleGroupAD", nrow(loso_fit_obj$metadata), length(dep_protein_cols), "Leave-one-site-out internal stability model.")
      dep_loso <- loso_fit_obj$dep %>% dplyr::mutate(excluded_site = site_i)
      export_dep_table(dep_loso, file.path(outdir, "result", "06_robustness", "site_robustness", "loso", "tables", paste0("LOSO_excluding_", safe_file_tag(site_i), "_dep_results.csv")))
      cmp <- DEP %>%
        dplyr::select(AptName, main_logFC = logFC, main_adj.P.Val = adj.P.Val) %>%
        dplyr::inner_join(dep_loso %>% dplyr::select(AptName, loso_logFC = logFC, loso_adj.P.Val = adj.P.Val), by = "AptName") %>%
        dplyr::mutate(same_direction = sign(main_logFC) == sign(loso_logFC), main_sig_fdr005 = main_adj.P.Val < MAIN_FDR, loso_sig_fdr005 = loso_adj.P.Val < MAIN_FDR, main_sig_preserved = main_sig_fdr005 & loso_sig_fdr005 & same_direction)
      loso_summary_metrics <- dplyr::bind_rows(loso_summary_metrics, tibble::tibble(
        excluded_site = site_i,
        n_samples = nrow(loso_fit_obj$metadata),
        logFC_correlation = suppressWarnings(cor(cmp$main_logFC, cmp$loso_logFC, use = "complete.obs")),
        direction_consistency_all = mean(cmp$same_direction, na.rm = TRUE),
        n_main_sig_fdr005 = sum(cmp$main_sig_fdr005, na.rm = TRUE),
        prop_main_sig_preserved = mean(cmp$main_sig_preserved[cmp$main_sig_fdr005], na.rm = TRUE)
      ))
    }
    safe_write_csv(loso_summary_metrics, file.path(outdir, "result", "06_robustness", "site_robustness", "loso", "tables", "LOSO_summary_metrics.csv"))
  }
} else if (RUN_FLAGS$run_site_robustness) {
  analysis_checkpoints <- add_checkpoint(analysis_checkpoints, "Site robustness skipped", NA_integer_, note = "No site/center/cohort variable detected.")
}

###############################################################################
# 08_balanced_resampling_and_formal_classification
###############################################################################

message("08_balanced_resampling_and_formal_classification")

balanced_summary_tbl <- NULL; balanced_protein_tbl <- NULL

# Defensive guard: if this section is run interactively without first running the
# country-robustness block, rebuild country_counts and countries_to_test here.
if (!exists("country_counts") || !exists("countries_to_test")) {
  country_counts <- dep_df %>%
    dplyr::count(Country, SampleGroup, name = "n") %>%
    tidyr::pivot_wider(names_from = SampleGroup, values_from = n, values_fill = 0) %>%
    ensure_count_cols(c("CN", "AD")) %>%
    dplyr::mutate(total = CN + AD) %>%
    dplyr::arrange(dplyr::desc(total))
  countries_to_test <- country_counts %>%
    dplyr::filter(CN >= MIN_N_PER_GROUP_COUNTRY, AD >= MIN_N_PER_GROUP_COUNTRY) %>%
    dplyr::pull(Country)
}

if (RUN_FLAGS$run_balanced_country_resampling && length(countries_to_test) >= 2) {
  eligible_country_tbl <- country_counts %>% dplyr::filter(Country %in% countries_to_test)
  n_per_country_group <- min(eligible_country_tbl$CN, eligible_country_tbl$AD, BALANCED_RESAMPLING_MIN_GROUP, na.rm = TRUE)
  if (is.finite(n_per_country_group) && n_per_country_group >= 3) {
    balanced_summary_list <- list(); balanced_protein_list <- list()
    for (iter in seq_len(BALANCED_RESAMPLING_NITER)) {
      sampled_ids <- dep_df %>%
        dplyr::filter(Country %in% countries_to_test) %>%
        dplyr::group_by(Country, SampleGroup) %>%
        dplyr::slice_sample(n = n_per_country_group, replace = FALSE) %>%
        dplyr::ungroup() %>%
        dplyr::pull(SampleId)
      dat_bal <- dep_df %>% dplyr::filter(SampleId %in% sampled_ids) %>% dplyr::mutate(Country = droplevels(factor(Country)), SampleGroup = factor(SampleGroup, levels = MAIN_GROUPS))
      balanced_formula <- "~ SampleGroup + Age + Sex + Country + Education"
      bal_fit <- tryCatch(run_limma_dep_model(dat_bal, dep_protein_cols, annot_tbl, balanced_formula, "SampleGroupAD", MAIN_FDR, model_name = paste0("balanced_country_resampling_iter_", iter)), error = function(e) NULL)
      if (is.null(bal_fit)) next
      cmp <- DEP %>%
        dplyr::select(AptName, main_logFC = logFC, main_adj.P.Val = adj.P.Val) %>%
        dplyr::inner_join(bal_fit$dep %>% dplyr::select(AptName, bal_logFC = logFC, bal_adj.P.Val = adj.P.Val), by = "AptName") %>%
        dplyr::mutate(iteration = iter, same_direction = sign(main_logFC) == sign(bal_logFC), main_sig_fdr005 = main_adj.P.Val < MAIN_FDR, bal_sig_fdr005 = bal_adj.P.Val < MAIN_FDR, main_sig_preserved = main_sig_fdr005 & bal_sig_fdr005 & same_direction)
      balanced_summary_list[[iter]] <- tibble::tibble(
        iteration = iter,
        n_samples = nrow(bal_fit$metadata),
        n_per_country_group = n_per_country_group,
        logFC_correlation = suppressWarnings(cor(cmp$main_logFC, cmp$bal_logFC, use = "complete.obs")),
        direction_consistency_all = mean(cmp$same_direction, na.rm = TRUE),
        prop_main_sig_preserved = mean(cmp$main_sig_preserved[cmp$main_sig_fdr005], na.rm = TRUE)
      )
      balanced_protein_list[[iter]] <- cmp %>% dplyr::select(iteration, AptName, bal_logFC, bal_adj.P.Val, same_direction, main_sig_preserved)
    }
    balanced_summary_tbl <- dplyr::bind_rows(balanced_summary_list)
    balanced_all_proteins <- dplyr::bind_rows(balanced_protein_list)
    safe_write_csv(balanced_summary_tbl, file.path(outdir, "result", "06_robustness", "balanced_country_resampling", "tables", "balanced_resampling_summary_metrics.csv"))
    if (nrow(balanced_all_proteins) > 0) {
      balanced_protein_tbl <- balanced_all_proteins %>%
        dplyr::group_by(AptName) %>%
        dplyr::summarise(
          mean_bal_logFC = mean(bal_logFC, na.rm = TRUE),
          sd_bal_logFC = stats::sd(bal_logFC, na.rm = TRUE),
          prop_same_direction = mean(same_direction, na.rm = TRUE),
          prop_preserved_if_main_sig = mean(main_sig_preserved, na.rm = TRUE),
          prop_sig_fdr005 = mean(bal_adj.P.Val < MAIN_FDR, na.rm = TRUE),
          n_iterations = dplyr::n(),
          .groups = "drop"
        ) %>%
        dplyr::left_join(DEP %>% dplyr::select(AptName, Protein_Name, EntrezGeneSymbol, main_logFC = logFC, main_adj.P.Val = adj.P.Val), by = "AptName")
      safe_write_csv(balanced_protein_tbl, file.path(outdir, "result", "06_robustness", "balanced_country_resampling", "tables", "balanced_resampling_protein_stability.csv"))
    }
  }
}

robustness_classification_tbl <- NULL; robustness_counts <- NULL
if (RUN_FLAGS$run_formal_robustness_classification) {
  robustness_classification_tbl <- DEP_gene %>%
    dplyr::select(AptName, Protein_Name, EntrezGeneSymbol, main_logFC = logFC, main_adj.P.Val = adj.P.Val, main_type = type) %>%
    dplyr::mutate(main_sig_fdr005 = main_adj.P.Val < MAIN_FDR)
  if (!is.null(main_vs_loco_mean)) robustness_classification_tbl <- robustness_classification_tbl %>% dplyr::left_join(main_vs_loco_mean %>% dplyr::select(AptName, mean_loco_logFC, prop_same_direction_loco = prop_same_direction), by = "AptName")
  if (!is.null(meta_tbl)) robustness_classification_tbl <- robustness_classification_tbl %>% dplyr::left_join(meta_tbl %>% dplyr::select(AptName, meta_logFC, meta_adj.P.Val, I2), by = "AptName")
  if (!is.null(balanced_protein_tbl)) robustness_classification_tbl <- robustness_classification_tbl %>% dplyr::left_join(balanced_protein_tbl %>% dplyr::select(AptName, mean_bal_logFC, prop_same_direction_balanced = prop_same_direction, prop_sig_fdr005_balanced = prop_sig_fdr005), by = "AptName")
  robustness_classification_tbl <- robustness_classification_tbl %>%
    dplyr::mutate(
      loco_preserved = !is.na(mean_loco_logFC) & sign(main_logFC) == sign(mean_loco_logFC) & dplyr::coalesce(prop_same_direction_loco, 0) >= 0.80,
      meta_preserved = !is.na(meta_logFC) & sign(main_logFC) == sign(meta_logFC) & dplyr::coalesce(meta_adj.P.Val, 1) < MAIN_FDR,
      balanced_preserved = !is.na(mean_bal_logFC) & sign(main_logFC) == sign(mean_bal_logFC) & dplyr::coalesce(prop_same_direction_balanced, 0) >= 0.80,
      robustness_score = as.integer(main_sig_fdr005) + as.integer(loco_preserved) + as.integer(meta_preserved) + as.integer(balanced_preserved),
      robustness_class = dplyr::case_when(
        robustness_score >= 4 ~ "High robustness",
        robustness_score == 3 ~ "Moderate robustness",
        robustness_score == 2 ~ "Partial robustness",
        TRUE ~ "Limited or unclassified"
      )
    ) %>%
    dplyr::arrange(dplyr::desc(robustness_score), main_adj.P.Val)
  robustness_counts <- robustness_classification_tbl %>% dplyr::count(robustness_class, name = "n") %>% dplyr::arrange(dplyr::desc(n))
  safe_write_csv(robustness_classification_tbl, file.path(outdir, "result", "06_robustness", "formal_classification", "protein_robustness_classification.csv"))
  safe_write_csv(robustness_counts, file.path(outdir, "result", "06_robustness", "formal_classification", "protein_robustness_classification_counts.csv"))
}

###############################################################################
# 09_manifest_workspace
###############################################################################

message("09_manifest_workspace")

analysis_checkpoints <- add_checkpoint(analysis_checkpoints, "Final DEP input samples", nrow(dep_df), note = "Samples used in the main AD vs CN DEP model")
analysis_checkpoints <- add_checkpoint(analysis_checkpoints, "Final gene-collapsed DEP rows", nrow(DEP_gene), note = "Default interpretation universe after deterministic gene collapse")

analysis_manifest <- tibble::tibble(
  item = c(
    "script", "project_root", "csv_file", "adat_file", "main_fdr", "strict_fdr", "main_groups",
    "n_raw_metadata", "n_merged_rowcheck_pass", "n_dep_samples", "n_dep_CN", "n_dep_AD",
    "n_seq_cols", "n_protein_universe_internal_adat", "n_gene_collapsed_universe",
    "n_main_dep_genes_fdr005", "n_main_dep_genes_fdr001",
    "apoe_sensitivity_available", "ad_only_cdrsb_severity_available", "cdrsb_adjusted_ad_vs_cn_available",
    "vascular_metabolic_sensitivity_available", "vascular_metabolic_covariates",
    "atn_adjusted_sensitivity_available", "atn_adjusted_covariates",
    "directional_ora_available", "loco_available", "country_meta_available", "site_robustness_available",
    "balanced_resampling_available", "formal_robustness_classification_available",
    "n_model_formulas_logged", "all_logged_designs_full_rank"
  ),
  value = c(
    "01_data_processing_and_differential_analysis", project_root, csv_file, adat_file,
    as.character(MAIN_FDR), as.character(STRICT_FDR), paste(MAIN_GROUPS, collapse = ","),
    as.character(nrow(meta_info_new)), as.character(nrow(sample_data)), as.character(nrow(dep_df)),
    as.character(sum(dep_df$SampleGroup == "CN", na.rm = TRUE)),
    as.character(sum(dep_df$SampleGroup == "AD", na.rm = TRUE)),
    as.character(length(seq_cols)), as.character(length(protein_universe)), as.character(nrow(DEP_gene)),
    as.character(sum(DEP_gene$adj.P.Val < MAIN_FDR, na.rm = TRUE)),
    as.character(sum(DEP_gene$adj.P.Val < STRICT_FDR, na.rm = TRUE)),
    as.character(!is.null(DEP_APOE_gene)), as.character(!is.null(DEP_CDRSB_gene)), as.character(!is.null(DEP_CDRSB_ADJ_gene)),
    as.character(!is.null(DEP_VASCULAR_gene)), paste(vascular_covariates, collapse = ","),
    as.character(!is.null(DEP_ATN_gene)), paste(preferred_atn_covariates, collapse = ","),
    as.character(exists("ora_directional_main")),
    as.character(!is.null(main_vs_loco_mean)), as.character(!is.null(meta_tbl)), as.character(!is.null(loso_summary_metrics)),
    as.character(!is.null(balanced_summary_tbl)), as.character(!is.null(robustness_classification_tbl)),
    as.character(nrow(model_formula_manifest)),
    as.character(all(model_design_diagnostics$full_rank, na.rm = TRUE))
  )
)

safe_write_csv(analysis_checkpoints, file.path(outdir, "result", "07_manifest", "analysis_checkpoints.csv"))
safe_write_csv(analysis_manifest, file.path(outdir, "result", "07_manifest", "analysis_manifest.csv"))
safe_write_csv(model_formula_manifest, file.path(outdir, "result", "07_manifest", "model_formula_manifest.csv"))
safe_write_csv(model_design_diagnostics, file.path(outdir, "result", "07_manifest", "model_diagnostics", "model_design_diagnostics.csv"))
writeLines(capture.output(utils::sessionInfo()), con = file.path(outdir, "result", "07_manifest", "sessionInfo.txt"))

workspace_objects <- c(
  "project_root", "outdir", "csv_file", "adat_file", "MAIN_GROUPS", "MAIN_FDR", "STRICT_FDR", "RUN_FLAGS",
  "sample_data_raw_merged", "sample_data", "meta_info_new", "my_adat", "soma_info_internal",
  "annot_tbl", "seq_cols", "protein_universe", "protein_universe_audit",
  "analysis_df", "normalized_expr_all", "normalized_expr_CN_AD", "normalized_expr", "protein_vars_present",
  "pca", "pca_var", "pca_df",
  "dep_df", "dep_protein_cols", "fit_main", "design_main",
  "DEP", "DEP_aptamer", "DEP_gene", "DEP_counts",
  "DEP_APOE", "DEP_APOE_gene", "apoe_compare_tbl",
  "DEP_CDRSB", "DEP_CDRSB_gene", "cdrsb_compare_tbl",
  "DEP_CDRSB_ADJ", "DEP_CDRSB_ADJ_gene", "cdrsb_adjusted_compare_tbl",
  "DEP_VASCULAR", "DEP_VASCULAR_gene", "vascular_compare_tbl", "vascular_covariates",
  "DEP_ATN", "DEP_ATN_gene", "atn_compare_tbl", "atn_covariates", "preferred_atn_covariates",
  "severity_results",
  "gsea_go_kegg_main", "gsea_reactome_main", "gsea_hallmark_main", "gsea_reactome_atn", "gsea_reactome_cdrsb", "gsea_reactome_cdrsb_adjusted", "ora_directional_main",
  "country_counts", "countries_to_test", "loco_tables", "loco_summary_metrics", "main_vs_loco_mean",
  "country_specific_tbl", "meta_tbl", "country_interaction_tbl", "site_var", "loso_summary_metrics",
  "balanced_summary_tbl", "balanced_protein_tbl", "robustness_classification_tbl", "robustness_counts",
  "summary_tbl", "na_count_tbl", "sample_tracking_tbl", "analysis_checkpoints", "analysis_manifest", "model_formula_manifest", "model_design_diagnostics"
)

existing_workspace_objects <- workspace_objects[vapply(workspace_objects, exists, logical(1), envir = .GlobalEnv)]
missing_workspace_objects <- setdiff(workspace_objects, existing_workspace_objects)
safe_write_csv(tibble::tibble(saved_object = existing_workspace_objects), file.path(outdir, "result", "07_manifest", "workspace_saved_objects.csv"))
if (length(missing_workspace_objects) > 0) safe_write_csv(tibble::tibble(missing_object = missing_workspace_objects), file.path(outdir, "result", "07_manifest", "workspace_missing_objects.csv"))

# Canonical workspace used by downstream scripts.
save(list = existing_workspace_objects, file = file.path(outdir, "result", "workspace", "analysis_workspace.RData"), envir = .GlobalEnv)

message("Primary data processing and differential analysis complete.")
message("Canonical workspace saved to: ", file.path(outdir, "result", "workspace", "analysis_workspace.RData"))
message("Manifest saved to: ", file.path(outdir, "result", "07_manifest", "analysis_manifest.csv"))
###############################################################################
# END OF 01_data_processing_and_differential_analysis.R
###############################################################################


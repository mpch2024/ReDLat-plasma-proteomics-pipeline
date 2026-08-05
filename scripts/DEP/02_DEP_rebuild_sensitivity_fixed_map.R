###############################################################################
# ReDLat plasma proteomics — DEP workflow
# 02. Sensitivity results on the primary SOMAmer map
# Requires: Outputs from Script 01
# Produces: Gene-level sensitivity tables and model comparisons
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

# =============================================================================
# 0. Packages and paths
# =============================================================================

required_pkgs <- c("dplyr", "tidyr", "purrr", "tibble", "readr", "stringr")
missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing required packages: ", paste(missing_pkgs, collapse = ", "),
    call. = FALSE
  )
}
invisible(lapply(required_pkgs, library, character.only = TRUE))


MAIN_FDR <- 0.05
STRICT_FDR <- 0.01
RUN_SENSITIVITY_REACTOME_GSEA <- TRUE
STRICT_EXPECTED_MAP_SIZE <- TRUE
EXPECTED_GENE_COLLAPSED_N <- 9638L

first_existing_file <- function(paths) {
  paths <- unique(paths[!is.na(paths) & nzchar(paths)])
  hit <- paths[file.exists(paths)][1]
  if (length(hit) == 0 || is.na(hit)) return(NA_character_)
  normalizePath(hit, winslash = "/", mustWork = TRUE)
}

workspace_candidates <- c(
  file.path(analysis_root, "workspace", "analysis_workspace.RData")
)

workspace_file <- first_existing_file(workspace_candidates)
if (is.na(workspace_file)) {
  stop(
    "Could not find analysis_workspace.RData. Run the primary DEP script first.",
    call. = FALSE
  )
}

analysis_root <- DEP_CONFIG$result_root
SENS_ROOT <- file.path(analysis_root, "04_sensitivity")
ENRICH_ROOT <- file.path(analysis_root, "05_enrichment_corrected", "gsea")
MANIFEST_ROOT <- file.path(analysis_root, "07_manifest", "fixed_primary_SOMAmer_map")
BACKUP_ROOT <- file.path(analysis_root, "workspace", "backups_before_fixed_primary_map")

dir.create(SENS_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(ENRICH_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(MANIFEST_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(BACKUP_ROOT, recursive = TRUE, showWarnings = FALSE)

message("Workspace: ", workspace_file)
message("Analysis root: ", analysis_root)

# Load into a dedicated environment so helper functions are not written into the
# analysis workspace.
ws <- new.env(parent = emptyenv())
load(workspace_file, envir = ws)

required_workspace_objects <- c("dep_df", "dep_protein_cols", "annot_tbl", "DEP_gene")
missing_workspace_objects <- required_workspace_objects[
  !vapply(required_workspace_objects, exists, logical(1), envir = ws)
]
if (length(missing_workspace_objects) > 0) {
  stop(
    "Workspace is missing required objects: ",
    paste(missing_workspace_objects, collapse = ", "),
    call. = FALSE
  )
}

# =============================================================================
# 1. Helpers
# =============================================================================

safe_write_csv <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(tibble::as_tibble(x), file)
}

safe_read_csv <- function(file) {
  if (is.na(file) || !file.exists(file)) return(NULL)
  readr::read_csv(file, show_col_types = FALSE, guess_max = 100000)
}

safe_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))

ensure_annotation_columns <- function(tbl) {
  needed <- c(
    "AptName", "SeqId", "EntrezGeneID", "EntrezGeneSymbol",
    "TargetFullName", "Target", "UniProt", "Protein_Name"
  )
  for (nm in needed) if (!nm %in% names(tbl)) tbl[[nm]] <- NA
  tbl
}

make_protein_name <- function(symbol, target_full, target, apt, seqid = NA_character_) {
  dplyr::coalesce(
    as.character(symbol), as.character(target_full), as.character(target),
    as.character(apt), as.character(seqid)
  )
}

standardize_dep_table <- function(dep_tbl, annot_tbl) {
  dep_tbl <- tibble::as_tibble(dep_tbl)

  if (!"AptName" %in% names(dep_tbl)) {
    candidate <- c("feature_id_raw", "feature_id", "SeqId")
    candidate <- candidate[candidate %in% names(dep_tbl)][1]
    if (length(candidate) == 0 || is.na(candidate)) {
      stop("DEP table has no AptName or compatible feature identifier.", call. = FALSE)
    }
    dep_tbl$AptName <- as.character(dep_tbl[[candidate]])
  }

  dep_tbl <- dep_tbl %>%
    dplyr::mutate(AptName = as.character(AptName))

  annot_tbl <- ensure_annotation_columns(tibble::as_tibble(annot_tbl)) %>%
    dplyr::mutate(AptName = as.character(AptName)) %>%
    dplyr::distinct(AptName, .keep_all = TRUE)

  annotation_fields <- c(
    "SeqId", "EntrezGeneID", "EntrezGeneSymbol", "TargetFullName",
    "Target", "UniProt", "Protein_Name"
  )

  for (nm in annotation_fields) {
    if (!nm %in% names(dep_tbl)) dep_tbl[[nm]] <- NA
  }

  annotation_add <- annot_tbl %>%
    dplyr::select(AptName, dplyr::all_of(annotation_fields)) %>%
    dplyr::rename_with(~ paste0(.x, ".annotation"), -AptName)

  dep_tbl <- dep_tbl %>%
    dplyr::left_join(annotation_add, by = "AptName")

  for (nm in annotation_fields) {
    ann_nm <- paste0(nm, ".annotation")
    if (nm == "EntrezGeneID") {
      dep_tbl[[nm]] <- dplyr::coalesce(
        safe_numeric(dep_tbl[[nm]]),
        safe_numeric(dep_tbl[[ann_nm]])
      )
    } else {
      dep_tbl[[nm]] <- dplyr::coalesce(
        as.character(dep_tbl[[nm]]),
        as.character(dep_tbl[[ann_nm]])
      )
    }
    dep_tbl[[ann_nm]] <- NULL
  }

  required_statistics <- c("logFC", "P.Value", "adj.P.Val")
  missing_statistics <- setdiff(required_statistics, names(dep_tbl))
  if (length(missing_statistics) > 0) {
    stop(
      "DEP table is missing required statistics: ",
      paste(missing_statistics, collapse = ", "),
      call. = FALSE
    )
  }

  dep_tbl %>%
    dplyr::mutate(
      EntrezGeneSymbol = as.character(EntrezGeneSymbol),
      Protein_Name = dplyr::coalesce(
        as.character(Protein_Name),
        make_protein_name(EntrezGeneSymbol, TargetFullName, Target, AptName, SeqId)
      ),
      logFC = safe_numeric(logFC),
      P.Value = safe_numeric(P.Value),
      adj.P.Val = safe_numeric(adj.P.Val),
      type = dplyr::case_when(
        adj.P.Val < MAIN_FDR & logFC > 0 ~ "Up",
        adj.P.Val < MAIN_FDR & logFC < 0 ~ "Down",
        TRUE ~ "NS"
      ),
      Direction = dplyr::case_when(
        type == "Up" ~ "Higher in AD",
        type == "Down" ~ "Lower in AD",
        TRUE ~ "Not significant"
      )
    )
}

build_primary_gene_map <- function(primary_gene_tbl) {
  required <- c("EntrezGeneSymbol", "AptName")
  missing_required <- setdiff(required, names(primary_gene_tbl))
  if (length(missing_required) > 0) {
    stop(
      "Primary DEP_gene is missing: ",
      paste(missing_required, collapse = ", "),
      call. = FALSE
    )
  }

  map <- primary_gene_tbl %>%
    dplyr::transmute(
      EntrezGeneSymbol = as.character(EntrezGeneSymbol),
      AptName = as.character(AptName)
    ) %>%
    dplyr::filter(
      !is.na(EntrezGeneSymbol), EntrezGeneSymbol != "",
      !is.na(AptName), AptName != ""
    ) %>%
    dplyr::distinct()

  if (anyDuplicated(map$EntrezGeneSymbol) > 0) {
    stop("Primary DEP_gene has more than one SOMAmer for at least one gene.", call. = FALSE)
  }
  if (anyDuplicated(map$AptName) > 0) {
    stop("Primary DEP_gene has a duplicated AptName in the fixed map.", call. = FALSE)
  }
  if (nrow(map) != nrow(primary_gene_tbl)) {
    stop(
      "Primary map has ", nrow(map), " rows but DEP_gene has ",
      nrow(primary_gene_tbl), ".", call. = FALSE
    )
  }
  if (STRICT_EXPECTED_MAP_SIZE && nrow(map) != EXPECTED_GENE_COLLAPSED_N) {
    stop(
      "Primary map has ", nrow(map), " rows; expected ",
      EXPECTED_GENE_COLLAPSED_N, ".", call. = FALSE
    )
  }
  map
}

apply_primary_map <- function(dep_aptamer, primary_map, model_name) {
  dep_aptamer <- dep_aptamer %>%
    dplyr::mutate(
      EntrezGeneSymbol = as.character(EntrezGeneSymbol),
      AptName = as.character(AptName)
    )

  fixed <- dep_aptamer %>%
    dplyr::inner_join(primary_map, by = c("EntrezGeneSymbol", "AptName")) %>%
    dplyr::arrange(match(EntrezGeneSymbol, primary_map$EntrezGeneSymbol))

  missing <- dplyr::anti_join(
    primary_map,
    fixed %>% dplyr::select(EntrezGeneSymbol, AptName),
    by = c("EntrezGeneSymbol", "AptName")
  )

  if (nrow(missing) > 0) {
    stop(
      model_name, " contains only ", nrow(fixed), " of ", nrow(primary_map),
      " primary SOMAmers. First missing pairs: ",
      paste(utils::head(paste0(missing$EntrezGeneSymbol, "=", missing$AptName), 10), collapse = "; "),
      call. = FALSE
    )
  }
  if (anyDuplicated(fixed$EntrezGeneSymbol) > 0 || anyDuplicated(fixed$AptName) > 0) {
    stop(model_name, " contains duplicated genes or SOMAmers after fixed-map selection.", call. = FALSE)
  }
  fixed
}

export_dep_table <- function(dep_tbl, file) {
  keep <- intersect(
    c(
      "Protein_Name", "EntrezGeneSymbol", "EntrezGeneID", "TargetFullName",
      "Target", "UniProt", "AptName", "SeqId", "feature_id_raw",
      "logFC", "se", "AveExpr", "t", "P.Value", "adj.P.Val", "B",
      "type", "Direction", "model_context", "coefficient_meaning",
      "severity_variable"
    ),
    names(dep_tbl)
  )
  dep_tbl %>%
    dplyr::select(dplyr::all_of(keep)) %>%
    dplyr::arrange(adj.P.Val, dplyr::desc(abs(logFC)), AptName) %>%
    safe_write_csv(file)
}

summarize_dep_counts <- function(dep_tbl, universe_label) {
  purrr::map_dfr(c(STRICT_FDR, MAIN_FDR), function(fdr) {
    tibble::tibble(
      universe = universe_label,
      fdr = fdr,
      sig_total = sum(dep_tbl$adj.P.Val < fdr, na.rm = TRUE),
      up = sum(dep_tbl$adj.P.Val < fdr & dep_tbl$logFC > 0, na.rm = TRUE),
      down = sum(dep_tbl$adj.P.Val < fdr & dep_tbl$logFC < 0, na.rm = TRUE)
    )
  })
}

export_fdr_specific <- function(dep_tbl, folder, prefix) {
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  export_dep_table(
    dep_tbl %>% dplyr::filter(adj.P.Val < STRICT_FDR),
    file.path(folder, paste0(prefix, "_FDR001.csv"))
  )
  export_dep_table(
    dep_tbl %>% dplyr::filter(adj.P.Val < MAIN_FDR),
    file.path(folder, paste0(prefix, "_FDR005.csv"))
  )
}

build_fixed_compare <- function(primary_gene, secondary_gene, label) {
  primary_std <- primary_gene %>%
    dplyr::select(
      EntrezGeneSymbol, AptName, SeqId, Protein_Name,
      logFC, adj.P.Val, P.Value, type
    ) %>%
    dplyr::rename(
      SeqId_primary = SeqId,
      Protein_Name_primary = Protein_Name,
      logFC_primary = logFC,
      adj.P.Val_primary = adj.P.Val,
      P.Value_primary = P.Value,
      type_primary = type
    )

  secondary_std <- secondary_gene %>%
    dplyr::select(
      EntrezGeneSymbol, AptName, SeqId, Protein_Name,
      logFC, adj.P.Val, P.Value, type
    ) %>%
    dplyr::rename(
      SeqId_secondary = SeqId,
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
      SeqId = dplyr::coalesce(SeqId_primary, SeqId_secondary),
      same_somamer = TRUE,
      same_direction = sign(logFC_primary) == sign(logFC_secondary),
      preserved_fdr005 = same_direction & adj.P.Val_secondary < MAIN_FDR,
      delta_logFC = logFC_secondary - logFC_primary,
      abs_ratio = abs(logFC_secondary) / pmax(abs(logFC_primary), 1e-9)
    )

  if (nrow(out) != nrow(primary_map)) {
    stop(
      label, " comparison contains ", nrow(out), " rows; expected ",
      nrow(primary_map), ".", call. = FALSE
    )
  }
  out
}

load_aptamer_result <- function(object_candidates, file_candidates, model_name) {
  for (nm in object_candidates) {
    if (exists(nm, envir = ws) && is.data.frame(get(nm, envir = ws))) {
      message(model_name, ": using workspace object ", nm)
      return(list(
        table = standardize_dep_table(get(nm, envir = ws), ws$annot_tbl),
        source = paste0("workspace::", nm)
      ))
    }
  }

  path <- first_existing_file(file_candidates)
  if (!is.na(path)) {
    message(model_name, ": reading ", path)
    return(list(
      table = standardize_dep_table(safe_read_csv(path), ws$annot_tbl),
      source = path
    ))
  }

  list(table = NULL, source = NA_character_)
}

run_reactome_gsea_optional <- function(dep_gene, output_file, analysis_label) {
  if (!RUN_SENSITIVITY_REACTOME_GSEA) {
    return(tibble::tibble(analysis = analysis_label, status = "disabled"))
  }
  needed <- c("ReactomePA", "BiocParallel")
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    warning(
      "Skipping Reactome GSEA for ", analysis_label,
      ": missing packages ", paste(missing, collapse = ", ")
    )
    return(tibble::tibble(
      analysis = analysis_label,
      status = "skipped_missing_packages",
      detail = paste(missing, collapse = ", ")
    ))
  }

  ranked <- dep_gene %>%
    dplyr::mutate(
      EntrezGeneID = safe_numeric(EntrezGeneID),
      rank_metric = safe_numeric(logFC)
    ) %>%
    dplyr::filter(
      !is.na(EntrezGeneID), is.finite(EntrezGeneID),
      !is.na(rank_metric), is.finite(rank_metric)
    ) %>%
    dplyr::group_by(EntrezGeneID) %>%
    dplyr::slice_max(abs(rank_metric), n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(dplyr::desc(rank_metric))

  if (nrow(ranked) < 10) {
    return(tibble::tibble(analysis = analysis_label, status = "skipped_too_few_mapped_genes"))
  }

  gene_list <- ranked$rank_metric
  names(gene_list) <- as.character(ranked$EntrezGeneID)
  gene_list <- sort(gene_list, decreasing = TRUE)

  suppressWarnings(
    BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)
  )
  result <- tryCatch(
    ReactomePA::gsePathway(
      geneList = gene_list,
      organism = "human",
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      verbose = FALSE
    ),
    error = function(e) e
  )

  if (inherits(result, "error")) {
    warning("Reactome GSEA failed for ", analysis_label, ": ", conditionMessage(result))
    return(tibble::tibble(
      analysis = analysis_label,
      status = "failed",
      detail = conditionMessage(result)
    ))
  }

  result_df <- as.data.frame(result)
  safe_write_csv(result_df, output_file)
  tibble::tibble(
    analysis = analysis_label,
    status = "completed",
    n_ranked_genes = length(gene_list),
    n_pathways = nrow(result_df),
    output_file = output_file
  )
}

backup_file <- function(path, backup_dir, suffix) {
  if (!file.exists(path)) return(NA_character_)
  dir.create(backup_dir, recursive = TRUE, showWarnings = FALSE)
  dest <- file.path(backup_dir, paste0(basename(path), ".", suffix, ".bak"))
  ok <- file.copy(path, dest, overwrite = TRUE)
  if (!ok) stop("Could not back up: ", path, call. = FALSE)
  dest
}

# =============================================================================
# 2. Freeze and audit the primary map
# =============================================================================

primary_before <- ws$DEP_gene
primary_gene <- standardize_dep_table(ws$DEP_gene, ws$annot_tbl)
primary_map <- build_primary_gene_map(primary_gene)

safe_write_csv(primary_map, file.path(MANIFEST_ROOT, "fixed_primary_SOMAmer_gene_map.csv"))

map_in_protein_universe <- primary_map$AptName %in% as.character(ws$dep_protein_cols)
map_audit <- tibble::tibble(
  item = c(
    "workspace_file",
    "analysis_root",
    "primary_DEP_gene_rows",
    "fixed_primary_map_rows",
    "expected_gene_collapsed_rows",
    "all_fixed_SOMAmers_in_dep_protein_cols",
    "selection_rule"
  ),
  value = c(
    workspace_file,
    analysis_root,
    as.character(nrow(primary_gene)),
    as.character(nrow(primary_map)),
    as.character(EXPECTED_GENE_COLLAPSED_N),
    as.character(all(map_in_protein_universe)),
    "Exact EntrezGeneSymbol--AptName pairs inherited from primary DEP_gene"
  )
)
safe_write_csv(map_audit, file.path(MANIFEST_ROOT, "fixed_primary_SOMAmer_gene_map_audit.csv"))

if (!all(map_in_protein_universe)) {
  missing_universe <- primary_map %>% dplyr::filter(!AptName %in% ws$dep_protein_cols)
  safe_write_csv(
    missing_universe,
    file.path(MANIFEST_ROOT, "fixed_primary_map_missing_from_dep_protein_cols.csv")
  )
  stop("Some primary SOMAmers are absent from dep_protein_cols.", call. = FALSE)
}

# =============================================================================
# 3. Registry and output rebuilding
# =============================================================================

registry <- list(
  APOE = list(
    model_name = "APOE_adjusted_DEP",
    objects = c("DEP_APOE"),
    aptamer_files = c(
      file.path(SENS_ROOT, "apoe", "AD_vs_CN_APOE_adjusted_full_limma_results_aptamer_level.csv")
    ),
    gene_file = file.path(SENS_ROOT, "apoe", "AD_vs_CN_APOE_adjusted_full_limma_results_gene_collapsed.csv"),
    compare_file = file.path(SENS_ROOT, "apoe", "primary_vs_APOE_adjusted_gene_comparison.csv"),
    count_file = file.path(SENS_ROOT, "apoe", "APOE_adjusted_DEP_counts_gene_collapsed.csv"),
    fdr_dir = file.path(SENS_ROOT, "apoe", "FDR_specific"),
    fdr_prefix = "AD_vs_CN_APOE_adjusted_DEP_gene_collapsed",
    reactome_file = file.path(ENRICH_ROOT, "apoe_adjusted_dep_gsea_reactome_bh.csv")
  ),
  CDRSB_AD_ONLY = list(
    model_name = "AD_only_CDRSB_severity",
    objects = c("DEP_CDRSB"),
    aptamer_files = c(
      file.path(SENS_ROOT, "cdrsb", "AD_only", "AD_only_CDRSB_severity_full_limma_results_aptamer_level.csv"),
      file.path(SENS_ROOT, "cdrsb", "AD_only_CDRSB_severity_full_limma_results_aptamer_level.csv")
    ),
    gene_file = file.path(SENS_ROOT, "cdrsb", "AD_only", "AD_only_CDRSB_severity_full_limma_results_gene_collapsed.csv"),
    compare_file = file.path(SENS_ROOT, "cdrsb", "AD_only", "primary_AD_vs_CN_vs_AD_only_CDRSB_severity_alignment.csv"),
    count_file = file.path(SENS_ROOT, "cdrsb", "AD_only", "AD_only_CDRSB_severity_counts_gene_collapsed.csv"),
    fdr_dir = file.path(SENS_ROOT, "cdrsb", "AD_only", "FDR_specific"),
    fdr_prefix = "AD_only_CDRSB_severity_DEP_gene_collapsed",
    reactome_file = file.path(ENRICH_ROOT, "ad_only_cdrsb_severity_gsea_reactome_bh.csv")
  ),
  CDRSB_ADJUSTED = list(
    model_name = "CDRSB_adjusted_AD_vs_CN_diagnostic_sensitivity",
    objects = c("DEP_CDRSB_ADJ", "DEP_CDRSB_adjusted"),
    aptamer_files = c(
      file.path(SENS_ROOT, "cdrsb", "AD_vs_CN_adjusted", "AD_vs_CN_CDRSB_adjusted_full_limma_results_aptamer_level.csv")
    ),
    gene_file = file.path(SENS_ROOT, "cdrsb", "AD_vs_CN_adjusted", "AD_vs_CN_CDRSB_adjusted_full_limma_results_gene_collapsed.csv"),
    compare_file = file.path(SENS_ROOT, "cdrsb", "AD_vs_CN_adjusted", "primary_vs_CDRSB_adjusted_AD_vs_CN_gene_comparison.csv"),
    count_file = file.path(SENS_ROOT, "cdrsb", "AD_vs_CN_adjusted", "CDRSB_adjusted_AD_vs_CN_DEP_counts_gene_collapsed.csv"),
    fdr_dir = file.path(SENS_ROOT, "cdrsb", "AD_vs_CN_adjusted", "FDR_specific"),
    fdr_prefix = "AD_vs_CN_CDRSB_adjusted_DEP_gene_collapsed",
    reactome_file = file.path(ENRICH_ROOT, "cdrsb_adjusted_ad_vs_cn_gsea_reactome_bh.csv")
  ),
  VASCULAR = list(
    model_name = "vascular_metabolic_adjusted_DEP",
    objects = c("DEP_VASCULAR"),
    aptamer_files = c(
      file.path(SENS_ROOT, "vascular_metabolic", "AD_vs_CN_vascular_metabolic_adjusted_full_limma_results_aptamer_level.csv")
    ),
    gene_file = file.path(SENS_ROOT, "vascular_metabolic", "AD_vs_CN_vascular_metabolic_adjusted_full_limma_results_gene_collapsed.csv"),
    compare_file = file.path(SENS_ROOT, "vascular_metabolic", "main_vs_vascular_metabolic_adjusted_gene_comparison.csv"),
    count_file = file.path(SENS_ROOT, "vascular_metabolic", "vascular_metabolic_adjusted_DEP_counts_gene_collapsed.csv"),
    fdr_dir = file.path(SENS_ROOT, "vascular_metabolic", "FDR_specific"),
    fdr_prefix = "AD_vs_CN_vascular_metabolic_adjusted_DEP_gene_collapsed",
    reactome_file = file.path(ENRICH_ROOT, "vascular_metabolic_adjusted_dep_gsea_reactome_bh.csv")
  ),
  ATN = list(
    model_name = "ATN_adjusted_DEP",
    objects = c("DEP_ATN"),
    aptamer_files = c(
      file.path(SENS_ROOT, "atn_adjusted", "AD_vs_CN_ATN_adjusted_full_limma_results_aptamer_level.csv")
    ),
    gene_file = file.path(SENS_ROOT, "atn_adjusted", "AD_vs_CN_ATN_adjusted_full_limma_results_gene_collapsed.csv"),
    compare_file = file.path(SENS_ROOT, "atn_adjusted", "primary_vs_ATN_adjusted_gene_comparison.csv"),
    count_file = file.path(SENS_ROOT, "atn_adjusted", "ATN_adjusted_DEP_counts_gene_collapsed.csv"),
    fdr_dir = file.path(SENS_ROOT, "atn_adjusted", "FDR_specific"),
    fdr_prefix = "AD_vs_CN_ATN_adjusted_DEP_gene_collapsed",
    reactome_file = file.path(ENRICH_ROOT, "atn_adjusted_dep_gsea_reactome_bh.csv")
  )
)

run_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
backup_dir <- file.path(BACKUP_ROOT, run_stamp)
dir.create(backup_dir, recursive = TRUE, showWarnings = FALSE)
workspace_backup <- backup_file(workspace_file, backup_dir, run_stamp)

processing_log <- list()
reactome_log <- list()
rebuilt <- list()

for (registry_name in names(registry)) {
  spec <- registry[[registry_name]]
  loaded <- load_aptamer_result(spec$objects, spec$aptamer_files, spec$model_name)

  if (is.null(loaded$table)) {
    warning("Skipping ", spec$model_name, ": no aptamer-level result was found.")
    processing_log[[registry_name]] <- tibble::tibble(
      analysis = registry_name,
      model = spec$model_name,
      status = "skipped_missing_aptamer_level_result",
      source = NA_character_,
      fixed_gene_rows = NA_integer_
    )
    next
  }

  # Back up outputs that will be replaced.
  files_to_backup <- c(spec$gene_file, spec$compare_file, spec$count_file)
  invisible(lapply(files_to_backup, backup_file, backup_dir = backup_dir, suffix = run_stamp))

  fixed_gene <- apply_primary_map(loaded$table, primary_map, spec$model_name)

  if (registry_name == "CDRSB_AD_ONLY") {
    if (!"severity_variable" %in% names(fixed_gene)) {
      fixed_gene$severity_variable <- NA_character_
    }
    fixed_gene <- fixed_gene %>%
      dplyr::mutate(
        model_context = "AD_only_CDRSB_severity",
        coefficient_meaning = "Within-AD CDR-SB severity slope for the fixed primary SOMAmer.",
        severity_variable = as.character(severity_variable)
      )

    comparison <- primary_gene %>%
      dplyr::select(
        EntrezGeneSymbol, AptName, SeqId, Protein_Name,
        logFC, adj.P.Val, P.Value, type
      ) %>%
      dplyr::rename(
        SeqId_primary = SeqId,
        Protein_Name_primary = Protein_Name,
        logFC_primary = logFC,
        adj.P.Val_primary = adj.P.Val,
        P.Value_primary = P.Value,
        type_primary = type
      ) %>%
      dplyr::inner_join(
        fixed_gene %>%
          dplyr::select(
            EntrezGeneSymbol, AptName, SeqId, Protein_Name,
            logFC, adj.P.Val, P.Value, type
          ) %>%
          dplyr::rename(
            SeqId_severity = SeqId,
            Protein_Name_severity = Protein_Name,
            beta_CDRSB_AD_only = logFC,
            adj.P.Val_CDRSB_AD_only = adj.P.Val,
            P.Value_CDRSB_AD_only = P.Value,
            type_CDRSB_AD_only = type
          ),
        by = c("EntrezGeneSymbol", "AptName")
      ) %>%
      dplyr::mutate(
        AptName_primary = AptName,
        AptName_severity = AptName,
        comparison = "Primary AD-vs-CN effect versus within-AD CDR-SB severity slope",
        Protein_Name = dplyr::coalesce(Protein_Name_primary, Protein_Name_severity),
        SeqId = dplyr::coalesce(SeqId_primary, SeqId_severity),
        same_somamer = TRUE,
        same_direction = sign(logFC_primary) == sign(beta_CDRSB_AD_only),
        primary_fdr005 = adj.P.Val_primary < MAIN_FDR,
        severity_fdr005 = adj.P.Val_CDRSB_AD_only < MAIN_FDR,
        primary_fdr005_and_same_direction = primary_fdr005 & same_direction,
        primary_fdr005_and_severity_fdr005_same_direction =
          primary_fdr005 & severity_fdr005 & same_direction,
        note = paste(
          "CDR-SB model is AD-only; beta_CDRSB_AD_only is not an adjusted",
          "AD-vs-CN logFC. The same primary SOMAmer is used in both models."
        )
      )
  } else {
    comparison <- build_fixed_compare(primary_gene, fixed_gene, spec$model_name)

    if (registry_name == "APOE") {
      comparison <- comparison %>%
        dplyr::rename(
          logFC_apoe = logFC_secondary,
          adj.P.Val_apoe = adj.P.Val_secondary,
          P.Value_apoe = P.Value_secondary,
          type_apoe = type_secondary
        )
    }

    if (registry_name == "CDRSB_ADJUSTED") {
      fixed_gene <- fixed_gene %>%
        dplyr::mutate(
          model_context = "CDRSB_adjusted_AD_vs_CN_diagnostic_sensitivity",
          coefficient_meaning = "AD versus CN diagnostic coefficient after CDR-SB adjustment, using the fixed primary SOMAmer."
        )
      comparison <- comparison %>%
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
          attenuation_ratio = dplyr::if_else(
            abs(logFC_primary) > 0,
            abs(logFC_CDRSB_adjusted) / abs(logFC_primary),
            NA_real_
          ),
          note = "Same primary SOMAmer used in both diagnostic models."
        )
    }

    if (registry_name == "VASCULAR") {
      comparison <- comparison %>%
        dplyr::rename(
          logFC_vascular = logFC_secondary,
          adj.P.Val_vascular = adj.P.Val_secondary,
          P.Value_vascular = P.Value_secondary,
          type_vascular = type_secondary
        )
    }

    if (registry_name == "ATN") {
      comparison <- comparison %>%
        dplyr::rename(
          logFC_atn = logFC_secondary,
          adj.P.Val_atn = adj.P.Val_secondary,
          P.Value_atn = P.Value_secondary,
          type_atn = type_secondary
        ) %>%
        dplyr::mutate(
          primary_fdr005 = adj.P.Val_primary < MAIN_FDR,
          atn_fdr005 = adj.P.Val_atn < MAIN_FDR,
          primary_fdr005_preserved_same_direction =
            primary_fdr005 & atn_fdr005 & same_direction,
          note = paste(
            "AT(N)-adjusted sensitivity uses the exact same primary SOMAmer",
            "for each gene."
          )
        )
    }
  }

  if (nrow(comparison) != nrow(primary_map)) {
    stop(
      spec$model_name, " comparison has ", nrow(comparison),
      " rows; expected ", nrow(primary_map), ".", call. = FALSE
    )
  }

  export_dep_table(fixed_gene, spec$gene_file)
  safe_write_csv(comparison, spec$compare_file)
  safe_write_csv(
    summarize_dep_counts(fixed_gene, paste0("gene_collapsed_fixed_primary_map_", spec$model_name)),
    spec$count_file
  )
  export_fdr_specific(fixed_gene, spec$fdr_dir, spec$fdr_prefix)

  reactome_log[[registry_name]] <- run_reactome_gsea_optional(
    fixed_gene,
    spec$reactome_file,
    spec$model_name
  )

  processing_log[[registry_name]] <- tibble::tibble(
    analysis = registry_name,
    model = spec$model_name,
    status = "completed",
    source = loaded$source,
    fixed_gene_rows = nrow(fixed_gene),
    same_map_rows = sum(
      fixed_gene$EntrezGeneSymbol == primary_map$EntrezGeneSymbol &
        fixed_gene$AptName == primary_map$AptName
    ),
    fdr005_total = sum(fixed_gene$adj.P.Val < MAIN_FDR, na.rm = TRUE),
    fdr001_total = sum(fixed_gene$adj.P.Val < STRICT_FDR, na.rm = TRUE)
  )

  rebuilt[[registry_name]] <- list(
    gene = fixed_gene,
    comparison = comparison
  )
}

# =============================================================================
# 4. Update canonical workspace objects
# =============================================================================

if ("APOE" %in% names(rebuilt)) {
  ws$DEP_APOE_gene <- rebuilt$APOE$gene
  ws$apoe_compare_tbl <- rebuilt$APOE$comparison
}
if ("CDRSB_AD_ONLY" %in% names(rebuilt)) {
  ws$DEP_CDRSB_gene <- rebuilt$CDRSB_AD_ONLY$gene
  ws$cdrsb_compare_tbl <- rebuilt$CDRSB_AD_ONLY$comparison
}
if ("CDRSB_ADJUSTED" %in% names(rebuilt)) {
  ws$DEP_CDRSB_ADJ_gene <- rebuilt$CDRSB_ADJUSTED$gene
  ws$cdrsb_adjusted_compare_tbl <- rebuilt$CDRSB_ADJUSTED$comparison
}
if ("VASCULAR" %in% names(rebuilt)) {
  ws$DEP_VASCULAR_gene <- rebuilt$VASCULAR$gene
  ws$vascular_compare_tbl <- rebuilt$VASCULAR$comparison
}
if ("ATN" %in% names(rebuilt)) {
  ws$DEP_ATN_gene <- rebuilt$ATN$gene
  ws$atn_compare_tbl <- rebuilt$ATN$comparison
}

# Preserve explicit map objects for downstream auditability.
ws$primary_gene_map_fixed <- primary_map
ws$primary_gene_map_rule <- "Exact EntrezGeneSymbol--AptName pairs inherited from primary DEP_gene"

primary_unchanged <- identical(primary_before, ws$DEP_gene)
if (!primary_unchanged) {
  stop("Safety check failed: DEP_gene changed in memory.", call. = FALSE)
}

# Save back to the canonical workspace. The original workspace was backed up.
save(list = ls(ws), file = workspace_file, envir = ws)

# Update compatible aliases only when they already exist.
workspace_aliases <- unique(c(
  file.path(analysis_root, "workspace", "proteomics_master_analysis_workspace.RData"),
  file.path(analysis_root, "workspace", "proteomics_master_reanalysis_workspace.RData"),
  file.path(dirname(analysis_root), "analysis_workspace.RData"),
  file.path(dirname(analysis_root), "proteomics_master_reanalysis_workspace.RData")
))
for (alias in workspace_aliases) {
  if (file.exists(alias)) {
    backup_file(alias, backup_dir, run_stamp)
    save(list = ls(ws), file = alias, envir = ws)
  }
}

processing_log_df <- dplyr::bind_rows(processing_log)
reactome_log_df <- dplyr::bind_rows(reactome_log)

run_audit <- tibble::tibble(
  item = c(
    "run_time",
    "workspace_file",
    "workspace_backup",
    "primary_DEP_unchanged",
    "fixed_map_rows",
    "completed_sensitivity_layers",
    "script_action"
  ),
  value = c(
    as.character(Sys.time()),
    workspace_file,
    workspace_backup,
    as.character(primary_unchanged),
    as.character(nrow(primary_map)),
    paste(processing_log_df$analysis[processing_log_df$status == "completed"], collapse = ", "),
    "Rebuilt gene-level sensitivity outputs from existing aptamer-level models using fixed primary SOMAmer map"
  )
)

safe_write_csv(processing_log_df, file.path(MANIFEST_ROOT, "01B_sensitivity_processing_log.csv"))
safe_write_csv(reactome_log_df, file.path(MANIFEST_ROOT, "01B_reactome_gsea_log.csv"))
safe_write_csv(run_audit, file.path(MANIFEST_ROOT, "01B_run_audit.csv"))
writeLines(capture.output(utils::sessionInfo()), file.path(MANIFEST_ROOT, "01B_sessionInfo.txt"))

message("Fixed-primary-map sensitivity rebuild complete.")
message("Primary DEP was not modified.")
message("Updated workspace: ", workspace_file)
message("Backup workspace: ", workspace_backup)
message("Audit folder: ", MANIFEST_ROOT)
###############################################################################
# END
###############################################################################


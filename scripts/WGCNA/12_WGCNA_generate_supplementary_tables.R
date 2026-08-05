###############################################################################
# ReDLat plasma proteomics — WGCNA workflow
# 12. Generate supplementary tables
# Requires: outputs from Scripts 01–09
# Produces: consolidated WGCNA supplementary tables
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

options(stringsAsFactors = FALSE)
options(error = traceback)

required_pkgs <- c("dplyr", "readr", "tibble", "purrr", "stringr", "openxlsx")
missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs) > 0L) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
       ". Run renv::restore() before this script.", call. = FALSE)
}
invisible(lapply(required_pkgs, library, character.only = TRUE))

run_supplementary_tables <- function() {

###############################################################################
# 00. SETTINGS
###############################################################################

BASE_DIR <- WGCNA_CONFIG$project_root

WGCNA_ROOT <- WGCNA_CONFIG$result_root

S10 <- file.path(WGCNA_ROOT, "01_input")
S11 <- file.path(WGCNA_ROOT, "02_network")
S12 <- file.path(WGCNA_ROOT, "03_modules")
S13 <- file.path(WGCNA_ROOT, "04_module_traits")
S13B <- file.path(WGCNA_ROOT, "05_sensitivity")
S14 <- file.path(WGCNA_ROOT, "06_stability")
S14B <- file.path(WGCNA_ROOT, "07_biomarker_fdr")
S15B <- file.path(WGCNA_ROOT, "08_preservation")
S16 <- file.path(WGCNA_ROOT, "09_network_quality")

OUTDIR <- file.path(WGCNA_CONFIG$publication_root, "supplementary_tables")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

OUTPUT_XLSX <- file.path(
  OUTDIR,
  "Supplementary_Tables_WGCNA_NatureAging_CLEAN.xlsx"
)
OUTPUT_INDEX <- file.path(
  OUTDIR,
  "Supplementary_Tables_WGCNA_NatureAging_Index.csv"
)
OUTPUT_AUDIT <- file.path(
  OUTDIR,
  "WGCNA_Current_Analysis_Audit.csv"
)
OUTPUT_REGISTRY <- file.path(
  OUTDIR,
  "WGCNA_Source_Registry.csv"
)
OUTPUT_MODULE_KEY <- file.path(
  OUTDIR,
  "WGCNA_Module_Key.csv"
)
OUTPUT_MANIFEST <- file.path(
  OUTDIR,
  "script19_final_manifest.csv"
)

EXPECTED <- list(
  genes = 9638L,
  samples = 639L,
  modules = 8L,
  module_trait = 128L,
  continuous_models = 48L,
  diagnosis_models = 8L,
  biomarker_models = 32L,
  fixed_preservation_genes = 4239L,
  soft_power = 4L
)

MODULE_KEY <- tibble::tribble(
  ~Module,   ~Module_ID, ~Module_display, ~Biological_label, ~Figure_role,
  "green",   "M1",       "M1/green",      "Neuronal-connectivity and extracellular-matrix module", "Focal",
  "blue",    "M2",       "M2/blue",       "AD-elevated high-differential-burden module", "Focal",
  "brown",   "M3",       "M3/brown",      "RNA-processing and ribonucleoprotein module", "Focal",
  "black",   "M4",       "M4/black",      "Biomarker-associated covariance module with low DEP burden", "Non-focal",
  "magenta", "M5",       "M5/magenta",    "Broad analytical module without a focal biological label", "Non-focal",
  "red",     "M6",       "M6/red",        "Analytical module without a focal biological label", "Non-focal",
  "pink",    "M7",       "M7/pink",       "Analytical module without a focal biological label", "Non-focal",
  "purple",  "M8",       "M8/purple",     "Analytical module without a focal biological label", "Non-focal"
)

readr::write_csv(MODULE_KEY, OUTPUT_MODULE_KEY)

###############################################################################
# 01. HELPERS
###############################################################################

first_existing_file <- function(paths, required = TRUE, label = "source") {
  paths <- unique(as.character(paths))
  paths <- paths[!is.na(paths) & nzchar(paths)]
  hit <- paths[file.exists(paths)][1]

  if (length(hit) == 0 || is.na(hit)) {
    if (required) {
      stop(
        "Required ", label, " was not found. Checked:\n",
        paste(paths, collapse = "\n"),
        call. = FALSE
      )
    }
    return(NA_character_)
  }

  normalizePath(hit, winslash = "/", mustWork = TRUE)
}

read_csv_safe <- function(path, required = TRUE, label = basename(path)) {
  if (is.na(path) || !file.exists(path)) {
    if (required) {
      stop("Missing required CSV: ", label, "\n", path, call. = FALSE)
    }
    return(tibble::tibble())
  }

  readr::read_csv(
    path,
    show_col_types = FALSE,
    guess_max = 100000,
    name_repair = "unique"
  )
}

sanitize_table <- function(x) {
  x <- as.data.frame(x, check.names = FALSE, stringsAsFactors = FALSE)

  for (i in seq_along(x)) {
    if (is.list(x[[i]])) {
      x[[i]] <- vapply(
        x[[i]],
        function(z) paste(as.character(z), collapse = "; "),
        character(1)
      )
    }
  }

  names(x) <- make.unique(
    ifelse(
      is.na(names(x)) | !nzchar(trimws(names(x))),
      paste0("Column_", seq_along(names(x))),
      names(x)
    ),
    sep = "_dup"
  )

  x
}

find_module_col <- function(tbl) {
  candidates <- c(
    "Module", "module", "Module_list_name", "Assigned_module",
    "assigned_module", "moduleColor", "Module_color"
  )
  hit <- candidates[candidates %in% names(tbl)][1]
  if (length(hit) == 0 || is.na(hit)) return(NA_character_)
  hit
}

add_module_identity <- function(tbl) {
  tbl <- tibble::as_tibble(tbl, .name_repair = "unique")
  module_col <- find_module_col(tbl)

  if (is.na(module_col)) return(tbl)

  # Normalize the analytical module identifier.
  tbl[[module_col]] <- sub("^ME", "", as.character(tbl[[module_col]]))

  # The same table may pass through this helper more than once. Remove all
  # previously attached display fields before joining the definitive key.
  identity_cols <- setdiff(names(MODULE_KEY), "Module")
  tbl <- dplyr::select(tbl, -dplyr::any_of(identity_cols))

  joined <- dplyr::left_join(
    tbl,
    MODULE_KEY,
    by = stats::setNames("Module", module_col)
  )

  dplyr::relocate(
    joined,
    dplyr::any_of(
      c(
        module_col,
        "Module_ID",
        "Module_display",
        "Biological_label",
        "Figure_role"
      )
    )
  )
}

pick_col <- function(tbl, candidates, required = FALSE, label = "column") {
  hit <- candidates[candidates %in% names(tbl)][1]
  if (length(hit) > 0 && !is.na(hit)) return(hit)

  clean <- function(x) tolower(gsub("[^a-z0-9]+", "", x))
  idx <- match(clean(candidates), clean(names(tbl)), nomatch = 0)
  idx <- idx[idx > 0][1]

  if (length(idx) > 0 && !is.na(idx)) return(names(tbl)[idx])

  if (required) {
    stop(
      "Could not identify ", label, ". Available columns: ",
      paste(names(tbl), collapse = ", "),
      call. = FALSE
    )
  }

  NA_character_
}

leading_enrichment_terms <- function(tbl, n_per_module = 5L) {
  tbl <- tibble::as_tibble(tbl, .name_repair = "unique")
  module_col <- find_module_col(tbl)

  if (is.na(module_col)) {
    stop(
      "The enrichment table does not contain a recognized module column.",
      call. = FALSE
    )
  }

  term_col <- pick_col(
    tbl,
    c("Description", "Term", "Pathway"),
    TRUE,
    "enrichment term"
  )
  fdr_col <- pick_col(
    tbl,
    c("p.adjust", "FDR", "adj.P.Val", "qvalue"),
    TRUE,
    "enrichment FDR"
  )

  out <- tbl %>%
    dplyr::mutate(
      .module_clean = sub("^ME", "", as.character(.data[[module_col]])),
      .fdr_numeric = suppressWarnings(as.numeric(.data[[fdr_col]]))
    ) %>%
    dplyr::filter(is.finite(.fdr_numeric), .fdr_numeric < 0.05) %>%
    dplyr::group_by(.module_clean) %>%
    dplyr::arrange(.fdr_numeric, .by_group = TRUE) %>%
    dplyr::distinct(.data[[term_col]], .keep_all = TRUE) %>%
    dplyr::slice_head(n = n_per_module) %>%
    dplyr::ungroup()

  # Replace the values of the existing module column. Never create a second
  # column called Module.
  out[[module_col]] <- out$.module_clean

  out <- dplyr::select(
    out,
    -.module_clean,
    -.fdr_numeric
  )

  add_module_identity(out)
}

relative_path <- function(path) {
  root <- normalizePath(BASE_DIR, winslash = "/", mustWork = FALSE)
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  prefix <- paste0(root, "/")
  if (startsWith(path, prefix)) {
    substring(path, nchar(prefix) + 1L)
  } else {
    path
  }
}

safe_sheet_name <- function(x, used = character()) {
  x <- gsub("[\\[\\]\\*\\?/\\\\:]", "_", x)
  x <- gsub("\\s+", "_", x)
  x <- substr(x, 1, 31)
  candidate <- x
  i <- 1L

  while (candidate %in% used) {
    suffix <- paste0("_", i)
    candidate <- paste0(substr(x, 1, 31 - nchar(suffix)), suffix)
    i <- i + 1L
  }

  candidate
}

###############################################################################
# 02. AUTHORITATIVE SOURCE REGISTRY
###############################################################################

P <- list(
  input_qc = file.path(S10, "wgcna_input_qc_summary.csv"),
  gene_map = file.path(S10, "tables", "outcome_independent_gene_somamer_map.csv"),
  core_summary = file.path(S11, "tables", "wgcna_core_summary.csv"),

  module_biology = file.path(S12, "tables", "integrated_module_biology_summary.csv"),
  dep_burden = file.path(S12, "tables", "module_DEP_burden_dual_definition.csv"),
  dep_overrep = file.path(S12, "tables", "module_DEP_overrepresentation_dual_definition.csv"),
  hub_summary = file.path(S12, "tables", "module_hub_summary.csv"),
  top_hubs = file.path(S12, "tables", "top_hubs_all_modules_combined.csv"),
  enrichment_summary = file.path(S12, "tables", "enrichment_summary_by_module.csv"),
  enrichment_all = file.path(S12, "tables", "enrichment_significant_terms_all_modules.csv"),

  module_trait = file.path(S13, "tables", "correlations", "module_trait_results_long.csv"),
  adjusted_continuous = file.path(S13, "tables", "regression", "adjusted_module_models.csv"),
  country_context = file.path(S13, "tables", "context", "module_recruitment_context_effect_summary.csv"),
  script13_summary = file.path(S13, "tables", "script13_final_summary.csv"),

  diagnosis_hc3 = file.path(S13B, "tables", "diagnosis_robustness", "adjusted_diagnosis_module_models_HC3.csv"),
  nested_site = file.path(S13B, "tables", "context", "corrected_nested_site_all_samples.csv"),

  loco = file.path(S14, "tables", "loco", "loco_country_summary_by_module_trait.csv"),
  wc_loso = file.path(S14, "tables", "within_country_loso", "within_country_loso_summary_by_module_trait.csv"),
  downsampling = file.path(S14, "tables", "downsampling", "balanced_downsampling_summary_by_module_trait.csv"),
  diagnosis_stability = file.path(S14, "tables", "model_stability", "diagnosis", "integrated_diagnosis_HC3_stability_summary.csv"),
  core_module_stability = file.path(S14, "tables", "core_modules", "core_modules_correlation_stability_summary.csv"),

  biomarker_full32 = file.path(S14B, "tables", "full", "full_all_modules_log_HC3_family32.csv"),
  biomarker_focal_stability = file.path(S14B, "tables", "focal_modules", "integrated_focal_log_HC3_stability_family32.csv"),
  biomarker_fdr_audit = file.path(S14B, "tables", "fdr_family_audit", "old_family12_vs_corrected_family32_by_model.csv"),

  fixed_gene_manifest = file.path(S15B, "tables", "input", "fixed_gene_manifest_all_modules.csv"),
  fixed_gene_manifest_summary = file.path(S15B, "tables", "input", "fixed_gene_manifest_summary_by_module.csv"),
  preservation_plan = file.path(S15B, "tables", "input", "reciprocal_within_country_site_run_plan.csv"),
  preservation_definitive = file.path(S15B, "tables", "fixed_geneset_structural_preservation_definitive_summary.csv"),
  preservation_country = file.path(S15B, "tables", "country", "country_fixed_geneset_primary_statistics.csv"),
  preservation_site = file.path(S15B, "tables", "site", "site_fixed_geneset_primary_statistics.csv"),
  script15b_summary = file.path(S15B, "tables", "script15b_final_summary.csv"),

  module_sizes = file.path(S16, "tables", "modularity", "final_module_sizes_and_proportions.csv"),
  weighted_q = file.path(S16, "tables", "modularity", "posthoc_weighted_modularity_Q.csv"),
  modularity_contributions = file.path(S16, "tables", "modularity", "module_specific_modularity_contributions.csv"),
  module_separation = file.path(S16, "tables", "modularity", "module_internal_external_separation_summary.csv"),
  network_quality = file.path(S16, "tables", "integration", "integrated_module_network_quality_summary.csv"),
  kme_quality = file.path(S16, "tables", "kME", "module_kME_quality_summary.csv"),
  eigengene_pairs = file.path(S16, "tables", "eigengenes", "module_eigengene_pairwise_correlations.csv"),
  eigengene_matrix = file.path(S16, "tables", "eigengenes", "module_eigengene_correlation_matrix.csv"),
  soft_threshold = file.path(S16, "tables", "soft_threshold", "soft_threshold_scan_clean.csv"),
  reviewer_summary = file.path(S16, "tables", "wgcna_network_quality_reviewer_summary.csv"),
  script16_summary = file.path(S16, "tables", "script16_final_summary.csv")
)

registry <- tibble::tibble(
  Source_ID = names(P),
  Path = unlist(P, use.names = FALSE),
  Required = TRUE
) %>%
  dplyr::mutate(
    Exists = file.exists(Path),
    Relative_path = vapply(Path, relative_path, character(1))
  )

readr::write_csv(registry, OUTPUT_REGISTRY)

missing_required <- registry %>% dplyr::filter(Required, !Exists)
if (nrow(missing_required) > 0) {
  readr::write_csv(
    missing_required,
    file.path(OUTDIR, "WGCNA_Missing_Required_Sources.csv")
  )
  stop(
    "Required WGCNA sources are missing. See WGCNA_Missing_Required_Sources.csv.",
    call. = FALSE
  )
}

TBL <- purrr::imap(P, ~ read_csv_safe(.x, TRUE, .y))

TOP_ENRICHMENT <- leading_enrichment_terms(
  TBL$enrichment_all,
  n_per_module = 5L
)

###############################################################################
# 03. CURRENT-ANALYSIS AUDIT
###############################################################################

family_col <- pick_col(TBL$biomarker_full32, c("FDR_family_size", "family_size"))
family_fdr_col <- pick_col(TBL$biomarker_full32, c("FDR_family32", "full_FDR_family32"))
fixed_comp_col <- pick_col(
  TBL$preservation_definitive,
  c("Fixed_gene_comparability", "fixed_gene_comparability")
)
power_col <- pick_col(TBL$soft_threshold, c("Power", "power"), TRUE, "soft-threshold power")

module_count <- if ("Module" %in% names(TBL$module_sizes)) {
  dplyr::n_distinct(TBL$module_sizes$Module)
} else {
  NA_integer_
}

fixed_gene_n <- nrow(TBL$fixed_gene_manifest)

audit <- tibble::tribble(
  ~Check, ~Expected, ~Observed, ~Pass,
  "Outcome-independent gene map", as.character(EXPECTED$genes), as.character(nrow(TBL$gene_map)), nrow(TBL$gene_map) == EXPECTED$genes,
  "Final WGCNA modules", as.character(EXPECTED$modules), as.character(module_count), identical(as.integer(module_count), EXPECTED$modules),
  "Complete module-trait family", as.character(EXPECTED$module_trait), as.character(nrow(TBL$module_trait)), nrow(TBL$module_trait) == EXPECTED$module_trait,
  "Adjusted continuous models", as.character(EXPECTED$continuous_models), as.character(nrow(TBL$adjusted_continuous)), nrow(TBL$adjusted_continuous) == EXPECTED$continuous_models,
  "Adjusted diagnosis models", as.character(EXPECTED$diagnosis_models), as.character(nrow(TBL$diagnosis_hc3)), nrow(TBL$diagnosis_hc3) == EXPECTED$diagnosis_models,
  "Corrected biomarker family", as.character(EXPECTED$biomarker_models), as.character(nrow(TBL$biomarker_full32)), nrow(TBL$biomarker_full32) == EXPECTED$biomarker_models,
  "Biomarker FDR column", "FDR_family32", ifelse(is.na(family_fdr_col), "missing", family_fdr_col), !is.na(family_fdr_col),
  "Biomarker family size", "32 in every row", ifelse(is.na(family_col), "missing", paste(unique(TBL$biomarker_full32[[family_col]]), collapse = "; ")), !is.na(family_col) && all(TBL$biomarker_full32[[family_col]] == 32, na.rm = TRUE),
  "Fixed preservation gene set", as.character(EXPECTED$fixed_preservation_genes), as.character(fixed_gene_n), fixed_gene_n == EXPECTED$fixed_preservation_genes,
  "Fixed-gene comparability field", "present", ifelse(is.na(fixed_comp_col), "missing", fixed_comp_col), !is.na(fixed_comp_col),
  "Selected soft power available", as.character(EXPECTED$soft_power), paste(sort(unique(TBL$soft_threshold[[power_col]])), collapse = ", "), EXPECTED$soft_power %in% TBL$soft_threshold[[power_col]],
  "M1-M8 module key", "8 unique identifiers", as.character(nrow(MODULE_KEY)), nrow(MODULE_KEY) == 8L
)

readr::write_csv(audit, OUTPUT_AUDIT)

if (!all(audit$Pass)) {
  stop(
    "Current-analysis audit failed. See:\n",
    OUTPUT_AUDIT,
    call. = FALSE
  )
}

###############################################################################
# 04. SUPPLEMENTARY TABLE ARCHITECTURE
###############################################################################

table_plan <- list(
  S20 = list(
    title = "Supplementary Table 20. WGCNA network construction and global quality diagnostics",
    note = paste(
      "Selected soft-threshold, module sizes, exact post hoc weighted modularity,",
      "module-specific modularity contributions, separation and integrated network-quality metrics.",
      "The modularity Q statistic is descriptive and was not used to optimize the network."
    ),
    blocks = list(
      list(title = "Reviewer-facing global network summary", data = TBL$reviewer_summary),
      list(title = "Soft-threshold scan", data = TBL$soft_threshold),
      list(title = "Final module sizes and proportions", data = TBL$module_sizes),
      list(title = "Exact post hoc weighted modularity", data = TBL$weighted_q),
      list(title = "Module-specific modularity contributions", data = TBL$modularity_contributions),
      list(title = "Integrated module network-quality summary", data = TBL$network_quality),
      list(title = "Module separation diagnostics", data = TBL$module_separation)
    )
  ),
  S21 = list(
    title = "Supplementary Table 21. Module composition, differential-abundance burden and biological summary",
    note = paste(
      "Module-level size, biological annotation and differential-abundance burden.",
      "Analyte-matched and canonical primary DEP definitions are reported separately."
    ),
    blocks = list(
      list(title = "Integrated module biology summary", data = TBL$module_biology),
      list(title = "Dual-definition DEP burden", data = TBL$dep_burden),
      list(title = "Dual-definition DEP overrepresentation", data = TBL$dep_overrep)
    )
  ),
  S22 = list(
    title = "Supplementary Table 22. Hub architecture and functional enrichment",
    note = paste(
      "Module membership and hub summaries together with leading FDR-significant enrichment terms.",
      "Counts of enriched terms are descriptive and are not normalized for module size or pathway redundancy."
    ),
    blocks = list(
      list(title = "Module-level hub summary", data = TBL$hub_summary),
      list(title = "Top hubs across all modules", data = TBL$top_hubs),
      list(title = "Enrichment summary by module and library", data = TBL$enrichment_summary),
      list(title = "Leading enrichment terms by module", data = TOP_ENRICHMENT)
    )
  ),
  S23 = list(
    title = "Supplementary Table 23. Complete module-trait association matrix",
    note = paste(
      "All 128 Spearman module-trait tests.",
      "BH-FDR correction was applied across the complete 8-module by 16-trait family."
    ),
    blocks = list(
      list(title = "Complete module-trait results", data = TBL$module_trait)
    )
  ),
  S24 = list(
    title = "Supplementary Table 24. Covariate-adjusted diagnosis, clinical and plasma biomarker models",
    note = paste(
      "Continuous-outcome models, diagnosis models with HC3 inference and all 32 log-transformed",
      "module-biomarker models. Biomarker FDR values come exclusively from Script 14b."
    ),
    blocks = list(
      list(title = "Forty-eight adjusted continuous-outcome models", data = TBL$adjusted_continuous),
      list(title = "Eight adjusted diagnosis models with HC3 standard errors", data = TBL$diagnosis_hc3),
      list(title = "Complete 32-model plasma biomarker family", data = TBL$biomarker_full32)
    )
  ),
  S25 = list(
    title = "Supplementary Table 25. Recruitment-context effects of country and site within country",
    note = paste(
      "Country and recruitment site were treated as categorical recruitment-context factors.",
      "Nested site estimates represent incremental site-within-country variation."
    ),
    blocks = list(
      list(title = "Country context-effect summary", data = TBL$country_context),
      list(title = "Corrected site-within-country models", data = TBL$nested_site)
    )
  ),
  S26 = list(
    title = "Supplementary Table 26. Stability of module associations across recruitment-context perturbations",
    note = paste(
      "Association stability was evaluated while holding the full-cohort module definitions fixed.",
      "LOCO, valid within-country site deletion and balanced downsampling are internal sensitivity analyses,",
      "not structural-preservation tests or external validation."
    ),
    blocks = list(
      list(title = "Leave-one-country-out module-trait stability", data = TBL$loco),
      list(title = "Corrected within-country leave-one-site-out stability", data = TBL$wc_loso),
      list(title = "Balanced-downsampling module-trait stability", data = TBL$downsampling),
      list(title = "Adjusted diagnosis HC3 stability", data = TBL$diagnosis_stability),
      list(title = "Focal biomarker stability using the corrected family of 32 models", data = TBL$biomarker_focal_stability),
      list(title = "Core-module correlation stability", data = TBL$core_module_stability)
    )
  ),
  S27 = list(
    title = "Supplementary Table 27. Fixed-gene structural module preservation",
    note = paste(
      "Structural preservation was evaluated using the identical fixed set of 4,239 genes in every comparison.",
      "Country-held-out and reciprocal within-country site analyses used 100 permutations.",
      "These results represent internal structural reproducibility rather than external validation."
    ),
    blocks = list(
      list(title = "Definitive preservation summary", data = TBL$preservation_definitive),
      list(title = "Country-held-out primary preservation statistics", data = TBL$preservation_country),
      list(title = "Reciprocal within-country site preservation statistics", data = TBL$preservation_site),
      list(title = "Fixed-gene manifest summary by module", data = TBL$fixed_gene_manifest_summary),
      list(title = "Reciprocal site-comparison plan", data = TBL$preservation_plan)
    )
  )
)

table_index <- purrr::imap_dfr(
  table_plan,
  function(x, id) {
    tibble::tibble(
      Table_ID = id,
      Supplementary_table_number = as.integer(sub("S", "", id)),
      Title = x$title,
      Purpose = x$note,
      Blocks = paste(vapply(x$blocks, `[[`, character(1), "title"), collapse = " | ")
    )
  }
)

readr::write_csv(table_index, OUTPUT_INDEX)

###############################################################################
# 05. WORKBOOK STYLES AND WRITER
###############################################################################

COL_TITLE <- "#CFC7B7"
COL_SECTION <- "#E9E5DC"
COL_HEADER <- "#F4F2EC"
COL_NOTE <- "#FAF9F6"
COL_BORDER <- "#4A4A4A"
COL_TEXT <- "#111111"

styles <- list(
  title = openxlsx::createStyle(
    fontName = "Arial", fontSize = 12, textDecoration = "bold",
    fgFill = COL_TITLE, fontColour = COL_TEXT,
    halign = "left", valign = "center",
    border = "TopBottomLeftRight", borderColour = COL_BORDER
  ),
  note = openxlsx::createStyle(
    fontName = "Arial", fontSize = 9, textDecoration = "italic",
    fgFill = COL_NOTE, fontColour = "#4A4A4A",
    wrapText = TRUE, valign = "top",
    border = "TopBottomLeftRight", borderColour = COL_BORDER
  ),
  section = openxlsx::createStyle(
    fontName = "Arial", fontSize = 10, textDecoration = "bold",
    fgFill = COL_SECTION, fontColour = COL_TEXT,
    halign = "left", valign = "center",
    border = "TopBottomLeftRight", borderColour = COL_BORDER
  ),
  header = openxlsx::createStyle(
    fontName = "Arial", fontSize = 9, textDecoration = "bold",
    fgFill = COL_HEADER, fontColour = COL_TEXT,
    halign = "center", valign = "center", wrapText = TRUE,
    border = "TopBottomLeftRight", borderColour = COL_BORDER
  ),
  body = openxlsx::createStyle(
    fontName = "Arial", fontSize = 8,
    valign = "top", wrapText = TRUE,
    border = "TopBottomLeftRight", borderColour = "#B7B7B7"
  )
)

write_multiblock_sheet <- function(wb, sheet, title, note, blocks) {
  sheet <- safe_sheet_name(sheet, openxlsx::sheets(wb))
  openxlsx::addWorksheet(wb, sheet, gridLines = FALSE)

  max_cols <- max(
    2L,
    vapply(
      blocks,
      function(b) max(1L, ncol(sanitize_table(b$data))),
      integer(1)
    )
  )

  openxlsx::mergeCells(wb, sheet, cols = seq_len(max_cols), rows = 1)
  openxlsx::writeData(wb, sheet, title, startRow = 1, startCol = 1)
  openxlsx::addStyle(wb, sheet, styles$title, rows = 1, cols = seq_len(max_cols), gridExpand = TRUE)
  openxlsx::setRowHeights(wb, sheet, 1, 26)

  openxlsx::mergeCells(wb, sheet, cols = seq_len(max_cols), rows = 2)
  openxlsx::writeData(wb, sheet, note, startRow = 2, startCol = 1)
  openxlsx::addStyle(wb, sheet, styles$note, rows = 2, cols = seq_len(max_cols), gridExpand = TRUE)
  openxlsx::setRowHeights(wb, sheet, 2, 40)

  row_cursor <- 4L

  for (block in blocks) {
    dat <- sanitize_table(add_module_identity(block$data))
    block_cols <- max(2L, ncol(dat))

    openxlsx::mergeCells(wb, sheet, cols = seq_len(block_cols), rows = row_cursor)
    openxlsx::writeData(wb, sheet, block$title, startRow = row_cursor, startCol = 1)
    openxlsx::addStyle(
      wb, sheet, styles$section,
      rows = row_cursor, cols = seq_len(block_cols), gridExpand = TRUE
    )
    openxlsx::setRowHeights(wb, sheet, row_cursor, 22)
    row_cursor <- row_cursor + 1L

    if (ncol(dat) == 0L) {
      dat <- data.frame(Note = "No rows available.", check.names = FALSE)
    }

    openxlsx::writeData(
      wb, sheet, dat,
      startRow = row_cursor, startCol = 1,
      colNames = TRUE, rowNames = FALSE, withFilter = FALSE
    )
    openxlsx::addStyle(
      wb, sheet, styles$header,
      rows = row_cursor, cols = seq_len(ncol(dat)), gridExpand = TRUE, stack = TRUE
    )

    if (nrow(dat) > 0) {
      openxlsx::addStyle(
        wb, sheet, styles$body,
        rows = (row_cursor + 1L):(row_cursor + nrow(dat)),
        cols = seq_len(ncol(dat)), gridExpand = TRUE, stack = TRUE
      )
    }

    numeric_cols <- which(vapply(dat, is.numeric, logical(1)))
    if (length(numeric_cols) > 0) {
      openxlsx::addStyle(
        wb, sheet,
        openxlsx::createStyle(numFmt = "0.000"),
        rows = (row_cursor + 1L):(row_cursor + max(1L, nrow(dat))),
        cols = numeric_cols,
        gridExpand = TRUE,
        stack = TRUE
      )
    }

    row_cursor <- row_cursor + nrow(dat) + 3L
  }

  openxlsx::freezePane(wb, sheet, firstActiveRow = 4)

  openxlsx::setColWidths(wb, sheet, cols = seq_len(max_cols), widths = 14)
  # Wider default columns keep multiblock sheets readable without relying on
  # block-specific column indices, which differ between sections.
  openxlsx::setColWidths(wb, sheet, cols = seq_len(max_cols), widths = 16)

  invisible(sheet)
}

wb <- openxlsx::createWorkbook(creator = "Matías Pizarro")

# Index
openxlsx::addWorksheet(wb, "Index", gridLines = FALSE)
openxlsx::writeData(wb, "Index", table_index, startRow = 1, headerStyle = styles$header)
openxlsx::addStyle(
  wb, "Index", styles$body,
  rows = 2:(nrow(table_index) + 1L),
  cols = seq_len(ncol(table_index)),
  gridExpand = TRUE
)
openxlsx::freezePane(wb, "Index", firstRow = TRUE)
openxlsx::setColWidths(wb, "Index", cols = 1:ncol(table_index), widths = c(10, 12, 42, 60, 70))

# Module key
openxlsx::addWorksheet(wb, "Module_Key", gridLines = FALSE)
openxlsx::writeData(wb, "Module_Key", MODULE_KEY, headerStyle = styles$header)
openxlsx::addStyle(
  wb, "Module_Key", styles$body,
  rows = 2:(nrow(MODULE_KEY) + 1L),
  cols = 1:ncol(MODULE_KEY),
  gridExpand = TRUE
)
openxlsx::freezePane(wb, "Module_Key", firstRow = TRUE)
openxlsx::setColWidths(wb, "Module_Key", cols = 1:ncol(MODULE_KEY), widths = c(12, 10, 16, 48, 14))

# Audit
openxlsx::addWorksheet(wb, "Current_Analysis_Audit", gridLines = FALSE)
openxlsx::writeData(wb, "Current_Analysis_Audit", audit, headerStyle = styles$header)
openxlsx::addStyle(
  wb, "Current_Analysis_Audit", styles$body,
  rows = 2:(nrow(audit) + 1L),
  cols = 1:ncol(audit),
  gridExpand = TRUE
)
openxlsx::setColWidths(wb, "Current_Analysis_Audit", cols = 1:ncol(audit), widths = c(36, 24, 40, 12))

for (id in names(table_plan)) {
  plan <- table_plan[[id]]
  write_multiblock_sheet(
    wb = wb,
    sheet = id,
    title = plan$title,
    note = plan$note,
    blocks = plan$blocks
  )
}

openxlsx::saveWorkbook(wb, OUTPUT_XLSX, overwrite = TRUE)

###############################################################################
# 06. README AND MANIFEST
###############################################################################

readme_lines <- c(
  "WGCNA SUPPLEMENTARY TABLES — NATURE AGING FINAL",
  "",
  "The workbook contains Supplementary Tables 20-27 in manuscript order.",
  "",
  "Authoritative rules:",
  "1. M1-M8 are stable visual identifiers and do not represent a ranking.",
  "2. Biomarker FDR values come only from Script 14b, complete family of 32 models.",
  "3. Structural preservation comes only from Script 15b, fixed 4,239-gene set.",
  "4. Network modularity comes only from Script 16, exact post hoc weighted Q.",
  "5. Association stability and structural preservation are reported separately.",
  "6. Complete 9,638-gene and high-dimensional results belong in Supplementary Data 4-6."
)
writeLines(readme_lines, file.path(OUTDIR, "README_WGCNA_Supplementary_Tables.txt"))

manifest <- tibble::tibble(
  Artifact = c(
    "Supplementary tables workbook",
    "Table index",
    "Current-analysis audit",
    "Source registry",
    "Module key",
    "README"
  ),
  Path = c(
    OUTPUT_XLSX,
    OUTPUT_INDEX,
    OUTPUT_AUDIT,
    OUTPUT_REGISTRY,
    OUTPUT_MODULE_KEY,
    file.path(OUTDIR, "README_WGCNA_Supplementary_Tables.txt")
  )
) %>%
  dplyr::mutate(
    Exists = file.exists(Path),
    Bytes = ifelse(Exists, file.info(Path)$size, NA_real_)
  )

readr::write_csv(manifest, OUTPUT_MANIFEST)
writeLines(capture.output(utils::sessionInfo()), file.path(OUTDIR, "sessionInfo.txt"))

if (!all(manifest$Exists) || any(manifest$Bytes <= 0, na.rm = TRUE)) {
  stop("One or more Script 19 outputs were not created correctly.", call. = FALSE)
}

writeLines(
  c(
    "SCRIPT_19_STATUS=PASS",
    paste0("WORKBOOK=", OUTPUT_XLSX),
    paste0("TIMESTAMP=", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  ),
  file.path(OUTDIR, "SCRIPT_19_SUCCESS.txt")
)

message("Script 19 completed successfully.")
message("Workbook: ", OUTPUT_XLSX)
message("Supplementary Tables: 20-27")
message("All current-analysis checks passed: TRUE")

} # end run_supplementary_tables()

message("============================================================")
message("RUNNING SCRIPT 12: SUPPLEMENTARY TABLES")
message("Duplicate Module protection: ENABLED")
message("============================================================")

run_supplementary_tables()

###############################################################################
# END
###############################################################################


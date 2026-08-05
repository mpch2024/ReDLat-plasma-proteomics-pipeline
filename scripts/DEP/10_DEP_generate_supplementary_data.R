###############################################################################
# ReDLat plasma proteomics — DEP workflow
# 10. DEP Supplementary Data
# Requires: Completed DEP analysis outputs
# Produces: Complete protein, pathway and robustness workbooks
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

required <- c("here", "readr", "dplyr", "tibble", "openxlsx", "stringr")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) stop("Missing packages: ", paste(missing, collapse = ", "))
invisible(lapply(required, library, character.only = TRUE))

out_dir <- file.path(publication_root, "supplementary_data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

COL_TITLE  <- "#CFC7B7"
COL_HEADER <- "#E9E5DC"
COL_NOTE   <- "#F4F2EC"
COL_BORDER <- "#4A4A4A"

styles <- list(
  title = createStyle(
    fontName = "Arial", fontSize = 12, textDecoration = "bold",
    fgFill = COL_TITLE, fontColour = "#111111",
    halign = "center", valign = "center",
    border = "TopBottomLeftRight", borderColour = COL_BORDER
  ),
  note = createStyle(
    fontName = "Arial", fontSize = 9, textDecoration = "italic",
    fgFill = COL_NOTE, fontColour = "#555555", wrapText = TRUE,
    border = "TopBottomLeftRight", borderColour = COL_BORDER
  ),
  header = createStyle(
    fontName = "Arial", fontSize = 9, textDecoration = "bold",
    fgFill = COL_HEADER, fontColour = "#111111",
    halign = "center", valign = "center", wrapText = TRUE,
    border = "TopBottomLeftRight", borderColour = COL_BORDER
  ),
  body = createStyle(
    fontName = "Arial", fontSize = 8, valign = "top", wrapText = TRUE,
    border = "TopBottomLeftRight", borderColour = COL_BORDER
  )
)

first_existing_file <- function(paths) {
  hit <- paths[file.exists(paths)][1]
  if (length(hit) == 0 || is.na(hit)) return(NA_character_)
  normalizePath(hit, winslash = "/", mustWork = TRUE)
}

relative_path <- function(path) {
  root <- normalizePath(project_root, winslash = "/", mustWork = FALSE)
  normalized <- normalizePath(path, winslash = "/", mustWork = FALSE)
  prefix <- paste0(root, "/")
  if (startsWith(normalized, prefix)) substring(normalized, nchar(prefix) + 1L) else basename(normalized)
}

read_required <- function(path, label) {
  path <- first_existing_file(path)
  if (is.na(path)) stop("Required source not found: ", label)
  readr::read_csv(path, show_col_types = FALSE, guess_max = 100000)
}

sanitize_sheet <- function(x) {
  x <- gsub("[\\[\\]\\*\\?/\\\\:]", "_", x)
  substr(x, 1, 31)
}

excel_safe_colnames <- function(x) {
  x <- as.character(x)
  x[is.na(x) | !nzchar(trimws(x))] <- "Column"
  x <- gsub("[\r\n\t]+", "_", x)
  x <- gsub("[^A-Za-z0-9_.]", "_", x)
  x <- gsub("_+", "_", x)
  x <- make.unique(x, sep = "_dup")
  if (anyDuplicated(x)) stop("Column-name repair failed; duplicated names remain.")
  x
}

write_sheet <- function(wb, sheet, title, note, data, table_id) {
  sheet <- sanitize_sheet(sheet)
  data <- as.data.frame(data, check.names = FALSE)
  names(data) <- excel_safe_colnames(names(data))
  ncols <- max(ncol(data), 2L)

  addWorksheet(wb, sheet, gridLines = FALSE)
  mergeCells(wb, sheet, cols = seq_len(ncols), rows = 1)
  writeData(wb, sheet, title, startRow = 1)
  addStyle(wb, sheet, styles$title, rows = 1, cols = seq_len(ncols), gridExpand = TRUE)
  setRowHeights(wb, sheet, rows = 1, heights = 25)

  mergeCells(wb, sheet, cols = seq_len(ncols), rows = 2)
  writeData(wb, sheet, note, startRow = 2)
  addStyle(wb, sheet, styles$note, rows = 2, cols = seq_len(ncols), gridExpand = TRUE)
  setRowHeights(wb, sheet, rows = 2, heights = 34)

  # Each Supplementary Data worksheet contains a single data block. Using
  # writeData(..., withFilter = TRUE) is more robust than writeDataTable for
  # very wide exported CSVs whose original headers may collide after Excel/openxlsx
  # normalization. The visible header, filter and submission formatting are retained.
  writeData(
    wb, sheet, data,
    startRow = 4, startCol = 1,
    colNames = TRUE,
    rowNames = FALSE,
    withFilter = TRUE
  )
  addStyle(wb, sheet, styles$header, rows = 4, cols = seq_len(ncol(data)), gridExpand = TRUE, stack = TRUE)
  if (nrow(data) > 0) {
    addStyle(
      wb, sheet, styles$body,
      rows = 5:(4 + nrow(data)), cols = seq_len(ncol(data)),
      gridExpand = TRUE, stack = TRUE
    )
  }
  freezePane(wb, sheet, firstActiveRow = 5)

  widths <- vapply(names(data), function(nm) {
    if (grepl("description|targetfullname|core_enrichment|leading_edge|note|formula|source", nm, ignore.case = TRUE)) 34 else 18
  }, numeric(1))
  setColWidths(wb, sheet, cols = seq_along(widths), widths = widths)
}

build_workbook <- function(filename, workbook_title, plan) {
  wb <- createWorkbook()

  index <- tibble(
    Sheet = vapply(plan, `[[`, character(1), "sheet"),
    Description = vapply(plan, `[[`, character(1), "title"),
    Source_file = vapply(plan, function(x) relative_path(first_existing_file(x$paths)), character(1))
  )
  write_sheet(
    wb, "Index", workbook_title,
    "Complete high-dimensional outputs supporting the consolidated DEP Supplementary Tables.",
    index, 1L
  )

  for (i in seq_along(plan)) {
    x <- plan[[i]]
    tbl <- read_required(x$paths, x$title)
    message("Writing sheet: ", x$sheet, " (", nrow(tbl), " rows; ", ncol(tbl), " columns)")
    write_sheet(wb, x$sheet, x$title, x$note, tbl, i + 1L)
  }

  path <- file.path(out_dir, filename)
  saveWorkbook(wb, path, overwrite = TRUE)
  message("Generated: ", path)
  invisible(path)
}

# Supplementary Data 1: complete protein-level results --------------------------
sd1 <- list(
  list(
    sheet = "Primary_DEP_full",
    title = "Primary AD-versus-CN differential-abundance results",
    note = "Complete fixed gene–SOMAmer map; one row per protein.",
    paths = c(file.path(analysis_root, "03_dep/gene_collapsed/AD_vs_CN_full_limma_results_gene_collapsed.csv"))
  ),
  list(
    sheet = "APOE_adjusted_full",
    title = "APOE ε4-adjusted differential-abundance results",
    note = "Complete fixed gene–SOMAmer map.",
    paths = c(file.path(analysis_root, "04_sensitivity/apoe/AD_vs_CN_APOE_adjusted_full_limma_results_gene_collapsed.csv"))
  ),
  list(
    sheet = "APOE_equal_sample",
    title = "Primary full, APOE complete-case baseline and APOE-adjusted comparison",
    note = "Separates complete-case restriction from APOE adjustment.",
    paths = c(file.path(analysis_root, "04_sensitivity/apoe/equal_sample_baseline/APOE_primary_full_vs_subset_baseline_vs_adjusted_FIXED_PRIMARY_MAP.csv"))
  ),
  list(
    sheet = "ATN_adjusted_full",
    title = "Plasma AT(N)-adjusted differential-abundance results",
    note = "Complete fixed gene–SOMAmer map.",
    paths = c(file.path(analysis_root, "04_sensitivity/atn_adjusted/AD_vs_CN_ATN_adjusted_full_limma_results_gene_collapsed.csv"))
  ),
  list(
    sheet = "ATN_equal_sample",
    title = "Primary full, AT(N) complete-case baseline and AT(N)-adjusted comparison",
    note = "Separates complete-case restriction from biomarker adjustment.",
    paths = c(file.path(analysis_root, "04_sensitivity/atn_adjusted/equal_sample_baseline/ATN_primary_full_vs_subset_baseline_vs_adjusted_FIXED_PRIMARY_MAP.csv"))
  ),
  list(
    sheet = "AD_only_CDRSB_full",
    title = "Within-AD CDR-SB severity associations",
    note = "Complete fixed gene–SOMAmer map.",
    paths = c(file.path(analysis_root, "04_sensitivity/cdrsb/AD_only/AD_only_CDRSB_severity_full_limma_results_gene_collapsed.csv"))
  ),
  list(
    sheet = "Primary_vs_CDRSB",
    title = "Primary diagnostic effects versus within-AD CDR-SB associations",
    note = "Effect-alignment comparison; not an independent diagnostic validation.",
    paths = c(file.path(analysis_root, "04_sensitivity/cdrsb/AD_only/primary_AD_vs_CN_vs_AD_only_CDRSB_severity_alignment.csv"))
  ),
  list(
    sheet = "pTau217_threshold_sweep",
    title = "Primary versus p-tau217-enriched threshold-sweep comparisons",
    note = "Complete p80, p90, p95 and stringent p90 comparisons.",
    paths = c(file.path(analysis_root, "04_sensitivity/p_tau217_enrichment/tables/primary_vs_p_tau217_enriched_threshold_sweep_all_comparisons_v3.csv"))
  ),
  list(
    sheet = "Matched_DEP_full",
    title = "Propensity-score-matched differential-abundance results",
    note = "Complete matched-sample gene-level output.",
    paths = c(file.path(analysis_root, "08_manuscript_supplementary_sensitivity/matched_DEP_gene_collapsed.csv"))
  ),
  list(
    sheet = "Primary_vs_matched",
    title = "Primary versus propensity-score-matched effects",
    note = "Complete effect comparison.",
    paths = c(file.path(analysis_root, "08_manuscript_supplementary_sensitivity/main_vs_matched_DEP_gene_comparison.csv"))
  ),
  list(
    sheet = "Biomarker_compatible_full",
    title = "Biomarker-compatible complete-case differential-abundance results",
    note = "Exploratory p-tau217 p75 and Aβ42/40 p25 contrast.",
    paths = c(file.path(analysis_root, "08_manuscript_supplementary_sensitivity/biomarker_consistent_DEP_gene_collapsed.csv"))
  ),
  list(
    sheet = "Primary_vs_biomarker",
    title = "Primary versus biomarker-compatible effects",
    note = "Complete effect comparison.",
    paths = c(file.path(analysis_root, "08_manuscript_supplementary_sensitivity/main_vs_biomarker_consistent_DEP_gene_comparison.csv"))
  )
)

build_workbook(
  "Supplementary_Data_1_DEP_Complete_Protein_Results.xlsx",
  "Supplementary Data 1. Complete DEP protein-level results",
  sd1
)

# Supplementary Data 2: complete pathway results -------------------------------
gsea_dir <- file.path(analysis_root, "05_enrichment_corrected", "gsea")
ora_dir  <- file.path(analysis_root, "05_enrichment_corrected", "ora")

sd2 <- list(
  list(sheet = "Primary_GSEA_GO", title = "Primary pre-ranked GO GSEA", note = "Complete BH-adjusted results.", paths = c(file.path(gsea_dir, "main_dep_gsea_go_bh.csv"))),
  list(sheet = "Primary_GSEA_KEGG", title = "Primary pre-ranked KEGG GSEA", note = "Complete BH-adjusted results.", paths = c(file.path(gsea_dir, "main_dep_gsea_kegg_bh.csv"))),
  list(sheet = "Primary_GSEA_Reactome", title = "Primary pre-ranked Reactome GSEA", note = "Complete BH-adjusted results.", paths = c(file.path(gsea_dir, "main_dep_gsea_reactome_bh.csv"))),
  list(sheet = "Primary_GSEA_Hallmark", title = "Primary pre-ranked Hallmark GSEA", note = "Complete BH-adjusted results.", paths = c(file.path(gsea_dir, "main_dep_gsea_hallmark_bh.csv"))),
  list(sheet = "Directional_ORA", title = "Primary direction-specific ORA", note = "Combined GO, KEGG and Reactome results for proteins higher or lower in AD.", paths = c(file.path(ora_dir, "main_dep_directional_ORA_GO_KEGG_Reactome_BH_combined.csv"))),
  list(sheet = "APOE_GSEA_Reactome", title = "APOE-adjusted Reactome GSEA", note = "Complete BH-adjusted results.", paths = c(file.path(gsea_dir, "apoe_adjusted_dep_gsea_reactome_bh.csv"))),
  list(sheet = "ATN_GSEA_Reactome", title = "AT(N)-adjusted Reactome GSEA", note = "Complete BH-adjusted results.", paths = c(file.path(gsea_dir, "atn_adjusted_dep_gsea_reactome_bh.csv"))),
  list(sheet = "CDRSB_GSEA_Reactome", title = "Within-AD CDR-SB Reactome GSEA", note = "Complete BH-adjusted results.", paths = c(file.path(gsea_dir, "ad_only_cdrsb_severity_gsea_reactome_bh.csv")))
)

build_workbook(
  "Supplementary_Data_2_DEP_Complete_Pathway_Results.xlsx",
  "Supplementary Data 2. Complete DEP pathway-enrichment results",
  sd2
)

# Supplementary Data 3: complete recruitment-context robustness ----------------
robust_root <- file.path(analysis_root, "06_robustness")

sd3 <- list(
  list(sheet = "LOCO_summary", title = "Leave-one-country-out summary", note = "Fixed-map model-level stability metrics.", paths = c(file.path(robust_root, "country_loco/tables/LOCO_summary_metrics_FIXED_PRIMARY_MAP.csv"))),
  list(sheet = "LOCO_all_proteins", title = "Primary versus mean LOCO effects", note = "Complete fixed-map protein-level output.", paths = c(file.path(robust_root, "country_loco/tables/main_vs_meanLOCO_table_FIXED_PRIMARY_MAP.csv"))),
  list(sheet = "Country_specific", title = "Country-specific differential-abundance effects", note = "Complete fixed-map country-specific effects.", paths = c(file.path(robust_root, "country_meta/tables/country_specific_DEP_results_FIXED_PRIMARY_MAP.csv"))),
  list(sheet = "Country_meta", title = "Country-level random-effects meta-analysis", note = "Complete fixed-map meta-analytic and heterogeneity statistics.", paths = c(file.path(robust_root, "country_meta/tables/country_meta_analysis_results_FIXED_PRIMARY_MAP.csv"))),
  list(sheet = "Diagnosis_country", title = "Diagnosis-by-country omnibus interaction", note = "Complete fixed-map interaction statistics.", paths = c(file.path(robust_root, "country_interaction/tables/country_interaction_omnibus_FIXED_PRIMARY_MAP.csv"))),
  list(sheet = "LOSO_summary", title = "Leave-one-site-out summary", note = "Fixed-map site-level stability metrics.", paths = c(file.path(robust_root, "site_robustness/loso/tables/LOSO_summary_metrics_FIXED_PRIMARY_MAP.csv"))),
  list(sheet = "Balanced_iterations", title = "Balanced country-resampling iterations", note = "One row per successful iteration.", paths = c(file.path(robust_root, "balanced_country_resampling/tables/balanced_resampling_summary_metrics_FIXED_PRIMARY_MAP.csv"))),
  list(sheet = "Balanced_proteins", title = "Balanced country-resampling protein stability", note = "Complete fixed-map protein-level stability output.", paths = c(file.path(robust_root, "balanced_country_resampling/tables/balanced_resampling_protein_stability_FIXED_PRIMARY_MAP.csv"))),
  list(sheet = "Robustness_class", title = "Formal protein robustness classification", note = "Complete fixed-map classification.", paths = c(file.path(robust_root, "formal_classification/protein_robustness_classification_FIXED_PRIMARY_MAP.csv"))),
  list(sheet = "Robustness_counts", title = "Formal robustness classification counts", note = "Summary counts by class.", paths = c(file.path(robust_root, "formal_classification/protein_robustness_classification_counts_FIXED_PRIMARY_MAP.csv"))),
  list(sheet = "Country_counts", title = "Country diagnostic-group counts", note = "Counts supporting country eligibility and balancing.", paths = c(file.path(robust_root, "country_loco/tables/country_group_counts.csv"))),
  list(sheet = "Site_counts", title = "Site diagnostic-group counts", note = "Counts supporting site eligibility.", paths = c(file.path(robust_root, "site_robustness/site_group_counts.csv")))
)

build_workbook(
  "Supplementary_Data_3_DEP_Complete_Robustness_Results.xlsx",
  "Supplementary Data 3. Complete DEP recruitment-context robustness results",
  sd3
)

message("\nAll three complete DEP Supplementary Data workbooks were generated.")
###############################################################################
# END
###############################################################################


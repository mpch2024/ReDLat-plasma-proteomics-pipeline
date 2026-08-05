###############################################################################
# ReDLat plasma proteomics — DEP workflow
# 08. DEP supplementary tables
# Requires: Completed DEP analysis outputs
# Produces: Consolidated supplementary table workbook
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

required <- c("here", "readr", "dplyr", "purrr", "tibble", "openxlsx", "stringr")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) stop("Missing packages: ", paste(missing, collapse = ", "))
invisible(lapply(required, library, character.only = TRUE))

out_dir <- file.path(publication_root, "supplementary_tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, "Supplementary_Tables_DEP_Consolidated_NatureAging.xlsx")

# Submission palette -----------------------------------------------------------
COL_TITLE   <- "#CFC7B7"
COL_SECTION <- "#D9D4C8"
COL_HEADER  <- "#E9E5DC"
COL_NOTE    <- "#F4F2EC"
COL_BORDER  <- "#4A4A4A"
COL_TEXT    <- "#111111"

st_title <- createStyle(fontName = "Arial", fontSize = 12, textDecoration = "bold",
                        fgFill = COL_TITLE, fontColour = COL_TEXT,
                        halign = "center", valign = "center",
                        border = "TopBottomLeftRight", borderColour = COL_BORDER)
st_section <- createStyle(fontName = "Arial", fontSize = 10, textDecoration = "bold",
                          fgFill = COL_SECTION, fontColour = COL_TEXT,
                          halign = "center", valign = "center",
                          border = "TopBottomLeftRight", borderColour = COL_BORDER)
st_header <- createStyle(fontName = "Arial", fontSize = 9, textDecoration = "bold",
                         fgFill = COL_HEADER, fontColour = COL_TEXT,
                         halign = "center", valign = "center", wrapText = TRUE,
                         border = "TopBottomLeftRight", borderColour = COL_BORDER)
st_body <- createStyle(fontName = "Arial", fontSize = 8, valign = "top", wrapText = TRUE,
                       border = "TopBottomLeftRight", borderColour = COL_BORDER)
st_note <- createStyle(fontName = "Arial", fontSize = 9, textDecoration = "italic",
                       fgFill = COL_NOTE, fontColour = "#555555", wrapText = TRUE,
                       border = "TopBottomLeftRight", borderColour = COL_BORDER)

read_required <- function(path) {
  if (!file.exists(path)) stop("Required file not found: ", path)
  readr::read_csv(path, show_col_types = FALSE, guess_max = 100000)
}

select_existing <- function(x, cols) dplyr::select(x, dplyr::any_of(cols))

# Excel supports only one worksheet-level AutoFilter. Each block is therefore
# written as an independent Excel table so every section retains its own filter.
block_table_counter <- 0L

write_title <- function(wb, sheet, title, note, ncols) {
  mergeCells(wb, sheet, cols = 1:ncols, rows = 1)
  writeData(wb, sheet, title, startRow = 1)
  addStyle(wb, sheet, st_title, rows = 1, cols = 1:ncols, gridExpand = TRUE)
  setRowHeights(wb, sheet, rows = 1, heights = 25)
  mergeCells(wb, sheet, cols = 1:ncols, rows = 2)
  writeData(wb, sheet, note, startRow = 2)
  addStyle(wb, sheet, st_note, rows = 2, cols = 1:ncols, gridExpand = TRUE)
  setRowHeights(wb, sheet, rows = 2, heights = 34)
}

write_block <- function(wb, sheet, start_row, section_title, x, note = NULL) {
  x <- as.data.frame(x, check.names = FALSE)
  ncols <- max(ncol(x), 1)
  mergeCells(wb, sheet, cols = 1:ncols, rows = start_row)
  writeData(wb, sheet, section_title, startRow = start_row)
  addStyle(wb, sheet, st_section, rows = start_row, cols = 1:ncols, gridExpand = TRUE)
  start_row <- start_row + 1
  if (!is.null(note)) {
    mergeCells(wb, sheet, cols = 1:ncols, rows = start_row)
    writeData(wb, sheet, note, startRow = start_row)
    addStyle(wb, sheet, st_note, rows = start_row, cols = 1:ncols, gridExpand = TRUE)
    start_row <- start_row + 1
  }
  if (nrow(x) == 0 || ncol(x) == 0) {
    writeData(wb, sheet, "No data available.", startRow = start_row)
    return(start_row + 3)
  }
  block_table_counter <<- block_table_counter + 1L
  table_stub <- gsub("[^A-Za-z0-9_]+", "_", section_title)
  table_stub <- gsub("^_+|_+$", "", table_stub)
  table_name <- substr(paste0("Tbl_", block_table_counter, "_", table_stub), 1, 250)

  writeDataTable(
    wb, sheet, x,
    startRow = start_row,
    startCol = 1,
    tableName = table_name,
    tableStyle = "TableStyleLight1",
    withFilter = TRUE
  )
  addStyle(wb, sheet, st_header, rows = start_row, cols = 1:ncol(x), gridExpand = TRUE, stack = TRUE)
  addStyle(wb, sheet, st_body,
           rows = (start_row + 1):(start_row + nrow(x)),
           cols = 1:ncol(x), gridExpand = TRUE, stack = TRUE)
  widths <- vapply(names(x), function(nm) {
    if (grepl("description|targetfullname|core_enrichment|leading_edge|note", nm, ignore.case = TRUE)) 34 else 18
  }, numeric(1))
  setColWidths(wb, sheet, cols = seq_along(widths), widths = widths)
  start_row + nrow(x) + 3
}

paths <- list(
  composition = file.path(analysis_root, "08_manuscript_supplementary_sensitivity", "supplementary_table_country_site_composition_consolidated.csv"),
  balance_overall = file.path(analysis_root, "06_robustness", "demographic_balance_reviewer_v3", "tables", "demographic_balance_overall_primary_complete_case.csv"),
  balance_country = file.path(analysis_root, "06_robustness", "demographic_balance_reviewer_v3", "tables", "demographic_balance_country_SMD_long_primary_complete_case.csv"),
  balance_site = file.path(analysis_root, "06_robustness", "demographic_balance_reviewer_v3", "tables", "demographic_balance_site_SMD_long_primary_complete_case.csv"),
  dep_counts = file.path(analysis_root, "03_dep", "AD_vs_CN_DEP_counts_by_universe.csv"),
  dep_fdr005 = file.path(analysis_root, "03_dep", "gene_collapsed", "FDR_specific", "AD_vs_CN_DEP_gene_collapsed_FDR005.csv"),
  apoe_summary = file.path(analysis_root, "04_sensitivity", "apoe", "equal_sample_baseline", "APOE_equal_sample_model_summary_FIXED_PRIMARY_MAP.csv"),
  apoe_fdr005 = file.path(analysis_root, "04_sensitivity", "apoe", "FDR_specific", "AD_vs_CN_APOE_adjusted_DEP_gene_collapsed_FDR005.csv"),
  atn_summary = file.path(analysis_root, "04_sensitivity", "atn_adjusted", "equal_sample_baseline", "ATN_equal_sample_model_summary_FIXED_PRIMARY_MAP.csv"),
  atn_fdr005 = file.path(analysis_root, "04_sensitivity", "atn_adjusted", "FDR_specific", "AD_vs_CN_ATN_adjusted_DEP_gene_collapsed_FDR005.csv"),
  cdr_counts = file.path(analysis_root, "04_sensitivity", "cdrsb", "AD_only", "AD_only_CDRSB_severity_counts_gene_collapsed.csv"),
  ptau_thresholds = file.path(analysis_root, "04_sensitivity", "p_tau217_enrichment", "tables", "cohort_defined_biomarker_thresholds_v3.csv"),
  ptau_summary = file.path(analysis_root, "04_sensitivity", "p_tau217_enrichment", "tables", "p_tau217_enriched_threshold_sweep_summary_v3.csv"),
  ptau_counts = file.path(analysis_root, "04_sensitivity", "p_tau217_enrichment", "tables", "p_tau217_enriched_threshold_sweep_sample_counts_v3.csv"),
  ptau_p90 = file.path(analysis_root, "04_sensitivity", "p_tau217_enrichment", "tables", "p_tau217_enriched_CNp90", "p_tau217_enriched_CNp90_DEP_gene_collapsed_FDR005.csv"),
  ptau_joint_p90 = file.path(analysis_root, "04_sensitivity", "p_tau217_enrichment", "tables", "p_tau217_enriched_low_abeta42_40_CNp90", "p_tau217_enriched_low_abeta42_40_CNp90_DEP_gene_collapsed_FDR005.csv"),
  matched_summary = file.path(analysis_root, "08_manuscript_supplementary_sensitivity", "main_vs_matched_DEP_summary.csv"),
  matched_balance = file.path(analysis_root, "08_manuscript_supplementary_sensitivity", "matched_covariate_balance_before_after.csv"),
  matched_full = file.path(analysis_root, "08_manuscript_supplementary_sensitivity", "matched_DEP_gene_collapsed.csv"),
  biomarker_summary = file.path(analysis_root, "08_manuscript_supplementary_sensitivity", "main_vs_biomarker_consistent_DEP_summary.csv"),
  biomarker_thresholds = file.path(analysis_root, "08_manuscript_supplementary_sensitivity", "biomarker_consistent_thresholds.csv"),
  biomarker_full = file.path(analysis_root, "08_manuscript_supplementary_sensitivity", "biomarker_consistent_DEP_gene_collapsed.csv"),
  loco_summary = file.path(analysis_root, "06_robustness", "country_loco", "tables", "LOCO_summary_metrics_FIXED_PRIMARY_MAP.csv"),
  loco_proteins = file.path(analysis_root, "06_robustness", "country_loco", "tables", "main_vs_meanLOCO_table_FIXED_PRIMARY_MAP.csv"),
  meta_summary = file.path(analysis_root, "06_robustness", "country_meta", "tables", "country_meta_summary_FIXED_PRIMARY_MAP.csv"),
  meta_full = file.path(analysis_root, "06_robustness", "country_meta", "tables", "country_meta_analysis_results_FIXED_PRIMARY_MAP.csv"),
  balanced = file.path(analysis_root, "06_robustness", "balanced_country_resampling", "tables", "balanced_resampling_summary_metrics_FIXED_PRIMARY_MAP.csv"),
  loso = file.path(analysis_root, "06_robustness", "site_robustness", "loso", "tables", "LOSO_summary_metrics_FIXED_PRIMARY_MAP.csv"),
  robustness_counts = file.path(analysis_root, "06_robustness", "formal_classification", "protein_robustness_classification_counts_FIXED_PRIMARY_MAP.csv"),
  interaction = file.path(analysis_root, "06_robustness", "country_interaction", "tables", "country_interaction_omnibus_FIXED_PRIMARY_MAP.csv")
)

protein_cols <- c("Protein_Name", "EntrezGeneSymbol", "TargetFullName", "UniProt", "AptName", "logFC", "se", "P.Value", "adj.P.Val", "Direction")
wb <- createWorkbook()

# Table 1 ----------------------------------------------------------------------
addWorksheet(wb, "Supplementary Table 1", gridLines = FALSE)
write_title(wb, "Supplementary Table 1", "Supplementary Table 1. Cohort composition and demographic balance",
            "CN, cognitively normal; AD, clinically diagnosed Alzheimer’s disease; SMD, standardized mean difference.", 26)
r <- 4
r <- write_block(wb, "Supplementary Table 1", r, "A. Composition overall and by recruitment context", read_required(paths$composition))
r <- write_block(wb, "Supplementary Table 1", r, "B. Overall demographic balance", read_required(paths$balance_overall))
r <- write_block(wb, "Supplementary Table 1", r, "C. Country-level standardized mean differences", read_required(paths$balance_country))
r <- write_block(wb, "Supplementary Table 1", r, "D. Site-level standardized mean differences", read_required(paths$balance_site))
freezePane(wb, "Supplementary Table 1", firstActiveRow = 4)

# Table 2 ----------------------------------------------------------------------
addWorksheet(wb, "Supplementary Table 2", gridLines = FALSE)
write_title(wb, "Supplementary Table 2", "Supplementary Table 2. Primary differential protein abundance",
            "The FDR < 0.05 list is reported once. The FDR < 0.01 subset is directly identifiable from adj.P.Val.", 10)
r <- 4
r <- write_block(wb, "Supplementary Table 2", r, "A. Protein-universe and threshold summary", read_required(paths$dep_counts))
r <- write_block(wb, "Supplementary Table 2", r, "B. Proteins associated with clinically diagnosed AD at FDR < 0.05",
                 select_existing(read_required(paths$dep_fdr005), protein_cols),
                 "Complete statistics for all 9,638 proteins belong in Supplementary Data 1.")
freezePane(wb, "Supplementary Table 2", firstActiveRow = 4)

# Table 3 ----------------------------------------------------------------------
addWorksheet(wb, "Supplementary Table 3", gridLines = FALSE)
write_title(wb, "Supplementary Table 3", "Supplementary Table 3. Pathway enrichment analyses",
            "GSEA is combined across GO, KEGG, Reactome and Hallmark; ORA is combined across database and direction.", 17)
r <- 4
gsea_files <- c(
  GO = file.path(analysis_root, "05_enrichment_corrected", "gsea", "main_dep_gsea_go_bh.csv"),
  KEGG = file.path(analysis_root, "05_enrichment_corrected", "gsea", "main_dep_gsea_kegg_bh.csv"),
  Reactome = file.path(analysis_root, "05_enrichment_corrected", "gsea", "main_dep_gsea_reactome_bh.csv"),
  Hallmark = file.path(analysis_root, "05_enrichment_corrected", "gsea", "main_dep_gsea_hallmark_bh.csv")
)
gsea <- purrr::imap_dfr(gsea_files, ~ read_required(.x) |> filter(p.adjust < 0.05) |> mutate(Database = .y, Model = "Primary AD versus CN", .before = 1))
r <- write_block(wb, "Supplementary Table 3", r, "A. Significant primary pre-ranked GSEA terms", gsea)
ora <- read_required(file.path(analysis_root, "05_enrichment_corrected", "ora", "main_dep_directional_ORA_GO_KEGG_Reactome_BH_combined.csv"))
r <- write_block(wb, "Supplementary Table 3", r, "B. Direction-specific over-representation analysis", ora)
sens_files <- c(
  `APOE-adjusted` = file.path(analysis_root, "05_enrichment_corrected", "gsea", "apoe_adjusted_dep_gsea_reactome_bh.csv"),
  `AT(N)-adjusted` = file.path(analysis_root, "05_enrichment_corrected", "gsea", "atn_adjusted_dep_gsea_reactome_bh.csv"),
  `AD-only CDR-SB` = file.path(analysis_root, "05_enrichment_corrected", "gsea", "ad_only_cdrsb_severity_gsea_reactome_bh.csv")
)
sens <- purrr::imap_dfr(sens_files, ~ read_required(.x) |> filter(p.adjust < 0.05) |> mutate(Database = "Reactome", Model = .y, .before = 1))
r <- write_block(wb, "Supplementary Table 3", r, "C. Significant Reactome GSEA terms in sensitivity models", sens)
freezePane(wb, "Supplementary Table 3", firstActiveRow = 4)

# Table 4 ----------------------------------------------------------------------
addWorksheet(wb, "Supplementary Table 4", gridLines = FALSE)
write_title(wb, "Supplementary Table 4", "Supplementary Table 4. APOE, plasma AT(N) and within-AD CDR-SB sensitivity analyses",
            "Equal-sample baselines distinguish complete-case restriction from covariate adjustment.", 17)
r <- 4
r <- write_block(wb, "Supplementary Table 4", r, "A. APOE ε4 equal-sample summary", read_required(paths$apoe_summary))
r <- write_block(wb, "Supplementary Table 4", r, "B. APOE-adjusted proteins at FDR < 0.05", select_existing(read_required(paths$apoe_fdr005), protein_cols))
r <- write_block(wb, "Supplementary Table 4", r, "C. Plasma AT(N) equal-sample summary", read_required(paths$atn_summary))
r <- write_block(wb, "Supplementary Table 4", r, "D. AT(N)-adjusted proteins at FDR < 0.05", select_existing(read_required(paths$atn_fdr005), protein_cols))
r <- write_block(wb, "Supplementary Table 4", r, "E. Within-AD CDR-SB significance counts", read_required(paths$cdr_counts))
cdr_summary <- tibble(Analysis = "Within-AD CDR-SB severity", Proteins_tested = 9638,
                      FDR_significant = 0, Effect_correlation_with_primary = -0.164,
                      Directional_consistency = 0.514,
                      Interpretation = "Severity-related effects did not recapitulate the primary CN-versus-AD signature.")
r <- write_block(wb, "Supplementary Table 4", r, "F. Within-AD CDR-SB alignment summary", cdr_summary)
freezePane(wb, "Supplementary Table 4", firstActiveRow = 4)

# Table 5 ----------------------------------------------------------------------
addWorksheet(wb, "Supplementary Table 5", gridLines = FALSE)
write_title(wb, "Supplementary Table 5", "Supplementary Table 5. Biomarker-informed and demographic-matching sensitivity analyses",
            "Cohort-derived thresholds are exploratory enrichment anchors and do not redefine clinical diagnosis.", 32)
r <- 4
r <- write_block(wb, "Supplementary Table 5", r, "A. Cohort-defined biomarker thresholds", read_required(paths$ptau_thresholds))
r <- write_block(wb, "Supplementary Table 5", r, "B. p-tau217 threshold-sweep summary", read_required(paths$ptau_summary))
r <- write_block(wb, "Supplementary Table 5", r, "C. p-tau217 threshold-sweep sample counts", read_required(paths$ptau_counts))
r <- write_block(wb, "Supplementary Table 5", r, "D. Primary p-tau217 CN p90 contrast: FDR < 0.05", select_existing(read_required(paths$ptau_p90), protein_cols))
r <- write_block(wb, "Supplementary Table 5", r, "E. Joint p-tau217–Aβ42/40 CN p90 contrast: FDR < 0.05", select_existing(read_required(paths$ptau_joint_p90), protein_cols))
r <- write_block(wb, "Supplementary Table 5", r, "F. Propensity-score matching summary", read_required(paths$matched_summary))
r <- write_block(wb, "Supplementary Table 5", r, "G. Covariate balance before and after matching", read_required(paths$matched_balance))
matched <- read_required(paths$matched_full) |> filter(significant_fdr005)
r <- write_block(wb, "Supplementary Table 5", r, "H. Matched proteins at FDR < 0.05", select_existing(matched, c(setdiff(protein_cols, "Direction"), "type")))
r <- write_block(wb, "Supplementary Table 5", r, "I. Biomarker-compatible analysis summary", read_required(paths$biomarker_summary))
r <- write_block(wb, "Supplementary Table 5", r, "J. Biomarker-compatible thresholds", read_required(paths$biomarker_thresholds))
biomarker <- read_required(paths$biomarker_full) |> filter(significant_fdr005)
r <- write_block(wb, "Supplementary Table 5", r, "K. Biomarker-compatible proteins at FDR < 0.05", select_existing(biomarker, c(setdiff(protein_cols, "Direction"), "type")))
freezePane(wb, "Supplementary Table 5", firstActiveRow = 4)

# Table 6 ----------------------------------------------------------------------
addWorksheet(wb, "Supplementary Table 6", gridLines = FALSE)
write_title(wb, "Supplementary Table 6", "Supplementary Table 6. Multicountry and recruitment-site robustness",
            "Internal stability analyses within ReDLat; not external validation or evidence of country-specific biology.", 16)
r <- 4
r <- write_block(wb, "Supplementary Table 6", r, "A. Leave-one-country-out summary", read_required(paths$loco_summary))
loco_core <- read_required(paths$loco_proteins) |> filter(all_loco_fdr005_same_direction)
r <- write_block(wb, "Supplementary Table 6", r, "B. Proteins preserved across every LOCO model", loco_core)
r <- write_block(wb, "Supplementary Table 6", r, "C. Country-level meta-analysis summary", read_required(paths$meta_summary))
meta_sig <- read_required(paths$meta_full) |> filter(meta_adj.P.Val < 0.05)
r <- write_block(wb, "Supplementary Table 6", r, "D. Multicountry-preserved proteins at meta-analysis FDR < 0.05", meta_sig)
bal <- read_required(paths$balanced)
bal_summary <- tibble(
  Iterations = nrow(bal), Participants_per_iteration = unique(bal$n_samples),
  Participants_per_country_group = unique(bal$n_per_country_group),
  Median_logFC_correlation = median(bal$logFC_correlation),
  IQR_logFC_correlation = paste0(round(quantile(bal$logFC_correlation, .25), 3), "–", round(quantile(bal$logFC_correlation, .75), 3)),
  Median_directional_consistency = median(bal$direction_consistency_all),
  IQR_directional_consistency = paste0(round(quantile(bal$direction_consistency_all, .25), 3), "–", round(quantile(bal$direction_consistency_all, .75), 3)),
  Median_primary_FDR_preservation = median(bal$prop_main_sig_preserved)
)
r <- write_block(wb, "Supplementary Table 6", r, "E. Balanced country-resampling summary", bal_summary)
r <- write_block(wb, "Supplementary Table 6", r, "F. Leave-one-site-out summary", read_required(paths$loso))
r <- write_block(wb, "Supplementary Table 6", r, "G. Formal robustness classification counts", read_required(paths$robustness_counts))
interaction <- read_required(paths$interaction)
interaction_summary <- tibble(Proteins_tested = nrow(interaction),
                              Interaction_FDR005 = sum(interaction$adj.P.Val < 0.05, na.rm = TRUE),
                              Interpretation = "Context-varying diagnostic effects were detectable; complete statistics are retained in Supplementary Data 3.")
r <- write_block(wb, "Supplementary Table 6", r, "H. Diagnosis-by-country interaction summary", interaction_summary)
freezePane(wb, "Supplementary Table 6", firstActiveRow = 4)

# Index ------------------------------------------------------------------------
addWorksheet(wb, "Index", gridLines = FALSE)
index <- tibble(
  Table = paste("Supplementary Table", 1:6),
  Title = c("Cohort composition and demographic balance", "Primary differential protein abundance",
            "Pathway enrichment analyses", "APOE, AT(N) and within-AD CDR-SB sensitivities",
            "Biomarker-informed and demographic-matching sensitivities", "Multicountry and recruitment-site robustness"),
  Complete_data = c("Supplementary Data 1", "Supplementary Data 1", "Supplementary Data 2",
                    "Supplementary Data 1", "Supplementary Data 1", "Supplementary Data 3")
)
write_title(wb, "Index", "DEP supplementary tables — consolidated Nature Aging set",
            "Six numbered tables; complete high-dimensional outputs are supplied separately as Supplementary Data.", 3)
write_block(wb, "Index", 4, "Table index", index)

saveWorkbook(wb, out_file, overwrite = TRUE)
message("Generated: ", out_file)


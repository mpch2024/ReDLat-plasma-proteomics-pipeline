###############################################################################
# ReDLat plasma proteomics
# 99_generate_final_supplementary_tables_MASTER.R
#
# PURPOSE
# Build one final manuscript-facing Supplementary Tables workbook following the
# order in which tables are cited in the manuscript.
#
# PRINCIPLES
# - One Supplementary Table = one worksheet.
# - No "S1", no "1a-1d" internal structure.
# - Numbering follows the manuscript order, not the original script order.
# - Machine-learning tables are intentionally excluded from this master workbook.
# - Source Data files are not mixed into Supplementary Tables.
#
# EXPECTED INPUTS
# Run the DEP and WGCNA scripts first, including:
#   DEP:
#     01_data_processing_and_differential_analysis.R
#     03_robustness_and_supplementary_analyses.R
#     04_DEP_supplementary_tables.R
#     05_additional_tables_sensitivity.R
#   WGCNA:
#     01_define_wgcna_input_from_DEP.R
#     02_wgcna_network_construction.R
#     03_WGCNA_module_biology_hubs_enrichment.R
#     04_WGCNA_module_trait_clinical_integration.R
#     04b_WGCNA_country_categorical_reviewer_analysis.R
#     05_WGCNA_country_site_downsampling_sensitivity.R
#
# MAIN OUTPUTS
#   result/final_supplementary_tables/Supplementary_Tables_FINAL.xlsx
#   result/final_supplementary_tables/Supplementary_Tables_FINAL_Index.csv
#   result/final_supplementary_tables/Supplementary_Tables_FINAL_MissingSources.csv
#
# AUTHOR
# Matias Pizarro + ChatGPT support
###############################################################################

rm(list = ls())

###############################################################################
# 00. Packages
###############################################################################

packages <- c(
  "readr", "openxlsx", "dplyr", "tibble", "purrr", "stringr", "tidyr"
)

installed <- rownames(installed.packages())
for (p in packages) {
  if (!p %in% installed) install.packages(p)
}

suppressPackageStartupMessages({
  library(readr)
  library(openxlsx)
  library(dplyr)
  library(tibble)
  library(purrr)
  library(stringr)
  library(tidyr)
})

###############################################################################
# 01. Paths - EDIT ONLY THESE IF NEEDED
###############################################################################

DEP_PROJECT_ROOT <- file.path(PROJECT_ROOT, "results", "DEP")
WGCNA_PROJECT_ROOT <- file.path(PROJECT_ROOT, "results", "WGCNA")

OUT_DIR <- file.path(DEP_PROJECT_ROOT,
  "final_supplementary_tables"
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUTPUT_XLSX <- file.path(OUT_DIR, "Supplementary_Tables_FINAL.xlsx")
OUTPUT_INDEX <- file.path(OUT_DIR, "Supplementary_Tables_FINAL_Index.csv")
OUTPUT_MISSING <- file.path(OUT_DIR, "Supplementary_Tables_FINAL_MissingSources.csv")

###############################################################################
# 02. Helper functions
###############################################################################

clean_sheet_name <- function(x) {
  x <- gsub("[\\[\\]\\*\\?/\\\\:]", "_", x)
  x <- gsub("\\s+", "_", x)
  substr(x, 1, 31)
}

clean_colnames <- function(df) {
  if (is.null(df)) return(df)
  nm <- names(df)
  nm[is.na(nm) | trimws(nm) == ""] <- paste0("Column_", seq_len(sum(is.na(nm) | trimws(nm) == "")))
  nm <- gsub("[^A-Za-z0-9_\\.]+", "_", nm)
  nm <- gsub("_+", "_", nm)
  nm <- gsub("^_|_$", "", nm)
  nm[nm == ""] <- "Column"
  names(df) <- make.unique(nm, sep = "_")
  df
}

sanitize_df <- function(df) {
  if (is.null(df)) return(NULL)
  df <- as.data.frame(df, check.names = FALSE, stringsAsFactors = FALSE)
  df <- clean_colnames(df)

  for (j in seq_along(df)) {
    if (is.list(df[[j]])) {
      df[[j]] <- vapply(df[[j]], function(z) paste(as.character(z), collapse = "; "), character(1))
    }
  }
  df
}

missing_df <- function(message, attempted_sources = NA_character_) {
  tibble::tibble(
    status = "MISSING_SOURCE",
    message = message,
    attempted_sources = paste(attempted_sources, collapse = " | ")
  )
}

first_existing <- function(paths) {
  paths <- paths[!is.na(paths)]
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

safe_read_csv <- function(path) {
  if (is.na(path) || !file.exists(path)) {
    return(missing_df("CSV source file was not found.", path))
  }
  suppressMessages(
    readr::read_csv(path, show_col_types = FALSE, guess_max = 100000)
  ) %>%
    as.data.frame(check.names = FALSE) %>%
    sanitize_df()
}

safe_read_xlsx_sheet <- function(file, sheet_candidates) {
  if (is.na(file) || !file.exists(file)) {
    return(missing_df("Excel source workbook was not found.", file))
  }

  sheets <- openxlsx::getSheetNames(file)
  sheet <- sheet_candidates[sheet_candidates %in% sheets][1]

  if (length(sheet) == 0 || is.na(sheet)) {
    return(missing_df(
      paste0("None of the requested sheets were found. Available sheets: ", paste(sheets, collapse = ", ")),
      paste0(file, "::", paste(sheet_candidates, collapse = "|"))
    ))
  }

  openxlsx::read.xlsx(file, sheet = sheet, detectDates = TRUE) %>%
    as.data.frame(check.names = FALSE) %>%
    sanitize_df()
}

read_any <- function(source) {
  type <- source$type %||% "csv"

  if (type == "csv") {
    path <- first_existing(source$paths)
    return(safe_read_csv(path))
  }

  if (type == "xlsx_sheet") {
    file <- first_existing(source$paths)
    return(safe_read_xlsx_sheet(file, source$sheets))
  }

  stop("Unknown source type: ", type)
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

# Helper used for combined tables that merge outputs from different analyses.
# Some enrichment files store the same column with different classes across
# databases, for example RichFactor can be numeric in one file and character in
# another. To avoid bind_rows() failures, combined-source tables are coerced to a
# manuscript-safe character representation before binding. This preserves values
# and makes the final Excel stable.
standardize_for_combo_bind <- function(df) {
  df <- sanitize_df(df)
  if (is.null(df)) return(NULL)
  df <- as.data.frame(df, check.names = FALSE, stringsAsFactors = FALSE)
  for (j in seq_along(df)) {
    if (inherits(df[[j]], "Date") || inherits(df[[j]], "POSIXct") || inherits(df[[j]], "POSIXt")) {
      df[[j]] <- as.character(df[[j]])
    } else if (is.factor(df[[j]])) {
      df[[j]] <- as.character(df[[j]])
    } else if (is.list(df[[j]])) {
      df[[j]] <- vapply(df[[j]], function(z) paste(as.character(z), collapse = "; "), character(1))
    } else {
      df[[j]] <- as.character(df[[j]])
    }
  }
  df
}

combine_sources <- function(sources, labels = NULL) {
  dfs <- purrr::map2(
    sources,
    seq_along(sources),
    function(src, i) {
      df <- read_any(src)
      lab <- if (!is.null(labels) && length(labels) >= i) labels[[i]] else paste0("source_", i)

      # Some ORA/GSEA source files can exist but contain no significant rows.
      # Keep an explicit placeholder row instead of failing during mutate/bind.
      if (is.null(df) || nrow(df) == 0 || ncol(df) == 0) {
        df <- tibble::tibble(
          source_section = lab,
          status = "EMPTY_SOURCE",
          message = "Source file exists but has no rows or columns."
        )
      } else {
        df <- df %>% mutate(source_section = lab, .before = 1)
      }

      standardize_for_combo_bind(df)
    }
  )

  # bind_rows can still fail if duplicate names or unusual objects slip through;
  # this final repair ensures all sources share a common schema.
  all_names <- unique(unlist(lapply(dfs, names)))
  dfs <- lapply(dfs, function(df) {
    missing_cols <- setdiff(all_names, names(df))
    for (mc in missing_cols) df[[mc]] <- NA_character_
    df <- df[, all_names, drop = FALSE]
    standardize_for_combo_bind(df)
  })

  dplyr::bind_rows(dfs)
}

write_table_sheet <- function(wb, sheet_name, df, title, note = NULL) {
  sheet_name <- clean_sheet_name(sheet_name)

  if (sheet_name %in% names(wb)) {
    openxlsx::removeWorksheet(wb, sheet_name)
  }

  openxlsx::addWorksheet(wb, sheet_name)

  title_style <- openxlsx::createStyle(
    fontName = "Arial",
    fontSize = 12,
    textDecoration = "bold",
    valign = "top",
    wrapText = TRUE
  )

  note_style <- openxlsx::createStyle(
    fontName = "Arial",
    fontSize = 9,
    fontColour = "#555555",
    textDecoration = "italic",
    valign = "top",
    wrapText = TRUE
  )

  header_style <- openxlsx::createStyle(
    fontName = "Arial",
    fontSize = 10,
    fontColour = "white",
    fgFill = "#1F4E79",
    textDecoration = "bold",
    halign = "center",
    valign = "center",
    border = "Bottom",
    borderColour = "#D9EAF7",
    wrapText = TRUE
  )

  body_style <- openxlsx::createStyle(
    fontName = "Arial",
    fontSize = 9,
    valign = "top",
    wrapText = TRUE
  )

  row_start <- 1
  openxlsx::writeData(wb, sheet_name, title, startRow = row_start, startCol = 1)
  openxlsx::addStyle(wb, sheet_name, title_style, rows = row_start, cols = 1)
  row_start <- row_start + 2

  if (!is.null(note) && !is.na(note) && nchar(note) > 0) {
    openxlsx::writeData(wb, sheet_name, note, startRow = row_start, startCol = 1)
    openxlsx::addStyle(wb, sheet_name, note_style, rows = row_start, cols = 1)
    row_start <- row_start + 2
  }

  df <- sanitize_df(df)
  if (is.null(df) || nrow(df) == 0 || ncol(df) == 0) {
    df <- tibble::tibble(status = "EMPTY_TABLE", message = "No rows or columns available.")
  }

  # Write as a table when possible. Fall back to plain data if Excel rejects a table.
  tryCatch(
    {
      openxlsx::writeDataTable(
        wb,
        sheet = sheet_name,
        x = df,
        startRow = row_start,
        startCol = 1,
        tableStyle = "TableStyleMedium2",
        withFilter = TRUE
      )
    },
    error = function(e) {
      openxlsx::writeData(
        wb,
        sheet = sheet_name,
        x = df,
        startRow = row_start,
        startCol = 1,
        withFilter = TRUE
      )
    }
  )

  openxlsx::addStyle(
    wb,
    sheet_name,
    header_style,
    rows = row_start,
    cols = seq_len(ncol(df)),
    gridExpand = TRUE,
    stack = TRUE
  )

  if (nrow(df) > 0) {
    openxlsx::addStyle(
      wb,
      sheet_name,
      body_style,
      rows = seq(row_start + 1, row_start + nrow(df)),
      cols = seq_len(ncol(df)),
      gridExpand = TRUE,
      stack = TRUE
    )
  }

  openxlsx::freezePane(wb, sheet_name, firstActiveRow = row_start + 1)

  for (j in seq_len(ncol(df))) {
    vals <- as.character(df[[j]])
    max_chars <- suppressWarnings(max(nchar(vals), na.rm = TRUE))
    if (!is.finite(max_chars)) max_chars <- 10
    width <- min(max(max_chars + 2, 10), 45)
    openxlsx::setColWidths(wb, sheet_name, cols = j, widths = width)
  }

  invisible(wb)
}

###############################################################################
# 03. Source paths
###############################################################################

# Manuscript-ready cohort and sensitivity tables generated by:
# 05_additional_tables_sensitivity_MANUSCRIPT_READY.R
#
# Older reviewer-facing versions are retained as fallbacks, but the script now
# searches the manuscript-ready directory first.
MANUSCRIPT_SENS_DIR <- file.path(DEP_PROJECT_ROOT,
  "08_manuscript_supplementary_sensitivity"
)

REVIEWER_AGUSTIN_DIR <- file.path(DEP_PROJECT_ROOT,
  "08_reviewer_agustin"
)

COHORT_XLSX <- c(
  file.path(MANUSCRIPT_SENS_DIR, "Supplementary_Table_1_Diagnostic_demographic_composition_by_recruitment_context.xlsx"),
  file.path(MANUSCRIPT_SENS_DIR, "Supplementary_Table_Cohort_Composition_MANUSCRIPT_READY.xlsx"),
  file.path(MANUSCRIPT_SENS_DIR, "Supplementary_Table_Country_Site_Composition_SUPPORTING.xlsx"),
  file.path(MANUSCRIPT_SENS_DIR, "Supplementary_Table_Country_Site_Composition_CONSOLIDATED.xlsx"),
  file.path(REVIEWER_AGUSTIN_DIR, "Supplementary_Table_1_Diagnostic_demographic_composition_by_recruitment_context.xlsx"),
  file.path(REVIEWER_AGUSTIN_DIR, "Supplementary_Table_Country_Site_Composition_SUPPORTING.xlsx"),
  file.path(REVIEWER_AGUSTIN_DIR, "Supplementary_Table_Country_Site_Composition_CONSOLIDATED.xlsx"),
  file.path(REVIEWER_AGUSTIN_DIR, "Agustin_Table1_and_CountrySite_Composition_MANUSCRIPT_STYLE.xlsx")
)

REVIEWER_SENS_XLSX <- c(
  file.path(MANUSCRIPT_SENS_DIR, "Supplementary_Table_2_Internal_sensitivity_DEP_signature.xlsx"),
  file.path(REVIEWER_AGUSTIN_DIR, "Supplementary_Table_2_Internal_sensitivity_DEP_signature.xlsx"),
  file.path(REVIEWER_AGUSTIN_DIR, "Reviewer_Sensitivity_DEP_Results.xlsx")
)

DEP_ROOT <- DEP_PROJECT_ROOT
DEP_GENE_DIR <- file.path(DEP_ROOT, "03_dep", "gene_collapsed")
DEP_FDR_DIR <- file.path(DEP_GENE_DIR, "FDR_specific")
GSEA_DIR <- file.path(DEP_ROOT, "05_enrichment_corrected", "gsea")
ORA_DIR  <- file.path(DEP_ROOT, "05_enrichment_corrected", "ora")
SENS_DIR <- file.path(DEP_ROOT, "04_sensitivity")
ROBUST_DIR <- file.path(DEP_ROOT, "06_robustness")

WGCNA_SCRIPT1 <- file.path(WGCNA_PROJECT_ROOT, "01_define_wgcna_input_from_DEP")
WGCNA_SCRIPT2 <- file.path(WGCNA_PROJECT_ROOT, "02_WGCNA_core_collapsed_genes")
WGCNA_SCRIPT3 <- file.path(WGCNA_PROJECT_ROOT, "03_WGCNA_module_biology_hubs_enrichment")
WGCNA_SCRIPT4 <- file.path(WGCNA_PROJECT_ROOT, "04_WGCNA_module_trait_clinical_integration")
WGCNA_SCRIPT4B <- file.path(WGCNA_PROJECT_ROOT, "04b_WGCNA_country_categorical_reviewer")
WGCNA_SCRIPT5 <- file.path(WGCNA_PROJECT_ROOT, "05_WGCNA_country_site_downsampling_sensitivity")

###############################################################################
# 04. Table plan - ordered by manuscript appearance
###############################################################################

tbl <- function(number, title, manuscript_section, citation_hint, source, note = NULL) {
  list(
    number = number,
    sheet = paste0("Supplementary_Table_", number),
    title = paste0("Supplementary Table ", number, ". ", title),
    short_title = title,
    manuscript_section = manuscript_section,
    citation_hint = citation_hint,
    source = source,
    note = note
  )
}

xlsx_src <- function(paths, sheets) list(type = "xlsx_sheet", paths = paths, sheets = sheets)
csv_src <- function(...) list(type = "csv", paths = c(...))
combo_src <- function(sources, labels = NULL) list(type = "combo", sources = sources, labels = labels)

# read_any wrapper for combo sources
read_table <- function(src) {
  if (!is.null(src$type) && src$type == "combo") {
    return(combine_sources(src$sources, src$labels))
  }
  read_any(src)
}

table_plan <- list(

  # Study cohort and clinical profile -----------------------------------------
  tbl(
    1,
    "Diagnostic and demographic composition by country.",
    "Study cohort and clinical profile",
    "diagnostic and demographic composition by country",
    xlsx_src(
      COHORT_XLSX,
      c("Suppl_Table_1b_Country", "By_country", "Supp_By_Country", "Country_composition", "Table_1b_Country", "Country_composition")
    )
  ),

  tbl(
    2,
    "Diagnostic and demographic composition by recruitment site.",
    "Study cohort and clinical profile",
    "diagnostic and demographic composition by recruitment site",
    xlsx_src(
      COHORT_XLSX,
      c("Suppl_Table_1c_Site", "By_recruitment_site", "Supp_By_Site", "Site_composition", "Table_1c_Site")
    )
  ),

  tbl(
    3,
    "Diagnostic composition tests by recruitment context.",
    "Study cohort and clinical profile",
    "recruitment-context composition tests",
    xlsx_src(
      COHORT_XLSX,
      c("Suppl_Table_1d_Tests", "Composition_Tests", "Composition_tests", "Table_1d_Tests")
    )
  ),

  # Differential plasma proteomics --------------------------------------------
  tbl(
    4,
    "Full gene-collapsed differential-abundance results for clinically diagnosed AD versus CN.",
    "Differential plasma proteomics",
    "full differential-abundance results",
    csv_src(file.path(DEP_GENE_DIR, "AD_vs_CN_full_limma_results_gene_collapsed.csv"))
  ),

  tbl(
    5,
    "Clinical AD-associated proteins at FDR < 0.05.",
    "Differential plasma proteomics",
    "FDR < 0.05 differentially abundant proteins",
    csv_src(file.path(DEP_FDR_DIR, "AD_vs_CN_DEP_gene_collapsed_FDR005.csv"))
  ),

  tbl(
    6,
    "High-confidence clinical AD-associated proteins at FDR < 0.01.",
    "Differential plasma proteomics",
    "FDR < 0.01 differentially abundant proteins",
    csv_src(file.path(DEP_FDR_DIR, "AD_vs_CN_DEP_gene_collapsed_FDR001.csv"))
  ),

  # Functional enrichment ------------------------------------------------------
  tbl(
    7,
    "Gene Ontology GSEA results for the clinical AD-associated plasma proteomic profile.",
    "Differential plasma proteomics",
    "Gene Ontology GSEA",
    csv_src(file.path(GSEA_DIR, "main_dep_gsea_go_bh.csv"))
  ),

  tbl(
    8,
    "KEGG GSEA results for the clinical AD-associated plasma proteomic profile.",
    "Differential plasma proteomics",
    "KEGG GSEA",
    csv_src(file.path(GSEA_DIR, "main_dep_gsea_kegg_bh.csv"))
  ),

  tbl(
    9,
    "Reactome GSEA results for the clinical AD-associated plasma proteomic profile.",
    "Differential plasma proteomics",
    "Reactome GSEA",
    csv_src(file.path(GSEA_DIR, "main_dep_gsea_reactome_bh.csv"))
  ),

  tbl(
    10,
    "Hallmark GSEA results for the clinical AD-associated plasma proteomic profile.",
    "Differential plasma proteomics",
    "Hallmark GSEA",
    csv_src(file.path(GSEA_DIR, "main_dep_gsea_hallmark_bh.csv"))
  ),

  tbl(
    11,
    "Directional over-representation analysis of proteins with higher or lower abundance in clinically diagnosed AD.",
    "Differential plasma proteomics",
    "directional over-representation analysis",
    combo_src(
      list(
        csv_src(file.path(ORA_DIR, "main_dep_higher_in_AD_GO_ORA_BH.csv")),
        csv_src(file.path(ORA_DIR, "main_dep_higher_in_AD_KEGG_ORA_BH.csv")),
        csv_src(file.path(ORA_DIR, "main_dep_higher_in_AD_Reactome_ORA_BH.csv")),
        csv_src(file.path(ORA_DIR, "main_dep_lower_in_AD_GO_ORA_BH.csv")),
        csv_src(file.path(ORA_DIR, "main_dep_lower_in_AD_KEGG_ORA_BH.csv")),
        csv_src(file.path(ORA_DIR, "main_dep_lower_in_AD_Reactome_ORA_BH.csv")),
        csv_src(file.path(ORA_DIR, "main_dep_directional_ORA_GO_KEGG_Reactome_BH_combined.csv"))
      ),
      labels = c(
        "Higher in AD - GO",
        "Higher in AD - KEGG",
        "Higher in AD - Reactome",
        "Lower in AD - GO",
        "Lower in AD - KEGG",
        "Lower in AD - Reactome",
        "Combined directional ORA"
      )
    )
  ),

  # Sensitivity and multicountry analyses --------------------------------------
  tbl(
    12,
    "AD-only CDR-SB severity model alignment with the primary diagnostic proteomic signature.",
    "Sensitivity and multicountry analyses",
    "AD-only severity alignment",
    csv_src(file.path(SENS_DIR, "cdrsb", "AD_only", "primary_AD_vs_CN_vs_AD_only_CDRSB_severity_alignment.csv"))
  ),

  tbl(
    13,
    "Full AD-only CDR-SB severity model results.",
    "Sensitivity and multicountry analyses",
    "AD-only CDR-SB protein-level model",
    csv_src(file.path(SENS_DIR, "cdrsb", "AD_only", "AD_only_CDRSB_severity_full_limma_results_gene_collapsed.csv"))
  ),

  tbl(
    14,
    "AD-only CDR-SB severity model summary and pathway-level results.",
    "Sensitivity and multicountry analyses",
    "AD-only CDR-SB counts and enrichment",
    combo_src(
      list(
        csv_src(file.path(SENS_DIR, "cdrsb", "AD_only", "AD_only_CDRSB_severity_counts_gene_collapsed.csv")),
        csv_src(file.path(GSEA_DIR, "ad_only_cdrsb_severity_gsea_reactome_bh.csv"))
      ),
      labels = c("AD-only CDR-SB DEP counts", "AD-only CDR-SB Reactome GSEA")
    )
  ),

  tbl(
    15,
    "APOE ε4-adjusted differential-abundance results.",
    "Sensitivity and multicountry analyses",
    "APOE ε4-adjusted protein-level model",
    csv_src(file.path(SENS_DIR, "apoe", "AD_vs_CN_APOE_adjusted_full_limma_results_gene_collapsed.csv"))
  ),

  tbl(
    16,
    "Comparison of primary and APOE ε4-adjusted differential-abundance models.",
    "Sensitivity and multicountry analyses",
    "primary versus APOE-adjusted comparison",
    csv_src(file.path(SENS_DIR, "apoe", "primary_vs_APOE_adjusted_gene_comparison.csv"))
  ),

  tbl(
    17,
    "AT(N)-adjusted differential-abundance results.",
    "Sensitivity and multicountry analyses",
    "AT(N)-adjusted protein-level model",
    csv_src(file.path(SENS_DIR, "atn_adjusted", "AD_vs_CN_ATN_adjusted_full_limma_results_gene_collapsed.csv"))
  ),

  tbl(
    18,
    "Comparison of primary and AT(N)-adjusted differential-abundance models.",
    "Sensitivity and multicountry analyses",
    "primary versus AT(N)-adjusted comparison",
    csv_src(file.path(SENS_DIR, "atn_adjusted", "primary_vs_ATN_adjusted_gene_comparison.csv"))
  ),

  tbl(
    19,
    "AT(N)-adjusted model covariates, DEP counts and Reactome enrichment.",
    "Sensitivity and multicountry analyses",
    "AT(N)-adjusted covariates, counts and enrichment",
    combo_src(
      list(
        csv_src(file.path(SENS_DIR, "atn_adjusted", "ATN_covariates_used.csv")),
        csv_src(file.path(SENS_DIR, "atn_adjusted", "ATN_adjusted_DEP_counts_gene_collapsed.csv")),
        csv_src(file.path(GSEA_DIR, "atn_adjusted_dep_gsea_reactome_bh.csv"))
      ),
      labels = c("AT(N) covariates used", "AT(N)-adjusted DEP counts", "AT(N)-adjusted Reactome GSEA")
    )
  ),

  tbl(
    20,
    "Leave-one-country-out diagnostic model robustness results.",
    "Sensitivity and multicountry analyses",
    "LOCO robustness",
    combo_src(
      list(
        csv_src(file.path(ROBUST_DIR, "country_loco", "tables", "LOCO_summary_metrics.csv")),
        csv_src(file.path(ROBUST_DIR, "country_loco", "tables", "main_vs_meanLOCO_table.csv"))
      ),
      labels = c("LOCO summary metrics", "Primary versus mean LOCO protein-level comparison")
    )
  ),

  tbl(
    21,
    "Country-level meta-analysis of clinical AD-associated proteomic effects.",
    "Sensitivity and multicountry analyses",
    "country-level meta-analysis",
    combo_src(
      list(
        csv_src(file.path(ROBUST_DIR, "country_meta", "tables", "country_meta_analysis_results.csv")),
        csv_src(file.path(ROBUST_DIR, "country_loco", "tables", "country_group_counts.csv"))
      ),
      labels = c("Country-level meta-analysis", "Country group counts")
    )
  ),

  tbl(
    22,
    "Balanced country-resampling sensitivity analysis.",
    "Sensitivity and multicountry analyses",
    "balanced resampling",
    csv_src(file.path(ROBUST_DIR, "balanced_country_resampling", "tables", "balanced_resampling_summary_metrics.csv"))
  ),

  tbl(
    23,
    "Leave-one-site-out recruitment-site sensitivity analysis.",
    "Sensitivity and multicountry analyses",
    "site-level sensitivity",
    csv_src(file.path(ROBUST_DIR, "site_robustness", "loso", "tables", "LOSO_summary_metrics.csv"))
  ),

  tbl(
    24,
    "Propensity-score matched sample characteristics.",
    "Sensitivity and multicountry analyses",
    "propensity-score matched sample",
    xlsx_src(REVIEWER_SENS_XLSX, c("Matched_sample_summary", "Matched_summary"))
  ),

  tbl(
    25,
    "Covariate balance before and after propensity-score matching.",
    "Sensitivity and multicountry analyses",
    "covariate balance before and after matching",
    xlsx_src(REVIEWER_SENS_XLSX, c("Balance_before_after", "Covariate_balance"))
  ),

  tbl(
    26,
    "Comparison of primary and propensity-score matched differential-abundance models.",
    "Sensitivity and multicountry analyses",
    "primary versus matched DEP comparison",
    xlsx_src(REVIEWER_SENS_XLSX, c("Main_vs_matched"))
  ),

  tbl(
    27,
    "Propensity-score matched differential-abundance results.",
    "Sensitivity and multicountry analyses",
    "matched DEP results",
    xlsx_src(REVIEWER_SENS_XLSX, c("Matched_DEP"))
  ),

  tbl(
    28,
    "Biomarker-compatible subgroup thresholds.",
    "Sensitivity and multicountry analyses",
    "biomarker-compatible thresholds",
    xlsx_src(REVIEWER_SENS_XLSX, c("Biomarker_thresholds"))
  ),

  tbl(
    29,
    "Biomarker-compatible subgroup composition.",
    "Sensitivity and multicountry analyses",
    "biomarker-compatible subgroup composition",
    xlsx_src(REVIEWER_SENS_XLSX, c("Biomarker_subset_summary"))
  ),

  tbl(
    30,
    "Comparison of primary and biomarker-compatible differential-abundance models.",
    "Sensitivity and multicountry analyses",
    "primary versus biomarker-compatible DEP comparison",
    xlsx_src(REVIEWER_SENS_XLSX, c("Main_vs_biomarker"))
  ),

  tbl(
    31,
    "Biomarker-compatible differential-abundance results.",
    "Sensitivity and multicountry analyses",
    "biomarker-compatible DEP results",
    xlsx_src(REVIEWER_SENS_XLSX, c("Biomarker_DEP"))
  ),

  # WGCNA ----------------------------------------------------------------------
  tbl(
    32,
    "WGCNA module assignment for gene-collapsed plasma proteins.",
    "Clinical AD-associated plasma proteomic variation shows modular organization",
    "WGCNA module assignments",
    csv_src(file.path(WGCNA_SCRIPT2, "tables", "gene_module_assignment.csv"))
  ),

  tbl(
    33,
    "WGCNA module sizes.",
    "Clinical AD-associated plasma proteomic variation shows modular organization",
    "module sizes",
    csv_src(file.path(WGCNA_SCRIPT2, "tables", "module_counts.csv"))
  ),

  tbl(
    34,
    "Module-trait associations without numerically encoded country.",
    "Clinical AD-associated plasma proteomic variation shows modular organization",
    "module-trait associations",
    csv_src(
      file.path(WGCNA_SCRIPT4B, "tables", "module_trait_results_no_country_numeric.csv"),
      file.path(WGCNA_SCRIPT4, "tables", "module_trait_results_long.csv")
    ),
    note = "Country is treated as a categorical recruitment-context variable in separate analyses and is not interpreted as an ordinal trait."
  ),

  tbl(
    35,
    "Module-level differential-abundance burden.",
    "Clinical AD-associated plasma proteomic variation shows modular organization",
    "module-level DEP burden",
    csv_src(file.path(WGCNA_SCRIPT4, "tables", "module_dep_summary.csv"))
  ),

  tbl(
    36,
    "Integrated WGCNA module prioritization summary.",
    "Clinical AD-associated plasma proteomic variation shows modular organization",
    "module prioritization",
    csv_src(file.path(WGCNA_SCRIPT4, "tables", "final_module_prioritization_table.csv"))
  ),

  tbl(
    37,
    "WGCNA functional enrichment summary by module.",
    "Clinical AD-associated plasma proteomic variation shows modular organization",
    "module enrichment",
    csv_src(file.path(WGCNA_SCRIPT3, "tables", "enrichment_summary_by_module.csv"))
  ),

  tbl(
    38,
    "WGCNA hub protein summary by module.",
    "Clinical AD-associated plasma proteomic variation shows modular organization",
    "module hubs",
    csv_src(file.path(WGCNA_SCRIPT3, "tables", "module_hub_summary.csv"))
  ),

  tbl(
    39,
    "Full kME and module membership table for assigned proteins.",
    "Clinical AD-associated plasma proteomic variation shows modular organization",
    "full module membership",
    csv_src(file.path(WGCNA_SCRIPT3, "tables", "full_kME_assigned_module_long.csv"))
  ),

  tbl(
    40,
    "WGCNA leave-one-country-out module-trait robustness analysis.",
    "Clinical AD-associated plasma proteomic variation shows modular organization",
    "WGCNA LOCO robustness",
    combo_src(
      list(
        csv_src(file.path(WGCNA_SCRIPT5, "tables", "loco", "loco_country_summary_by_module.csv")),
        csv_src(file.path(WGCNA_SCRIPT5, "tables", "loco", "loco_country_summary_by_module_trait.csv"))
      ),
      labels = c("LOCO summary by module", "LOCO summary by module-trait")
    )
  ),

  tbl(
    41,
    "WGCNA recruitment-site and balanced downsampling robustness analyses.",
    "Clinical AD-associated plasma proteomic variation shows modular organization",
    "WGCNA site and downsampling robustness",
    combo_src(
      list(
        csv_src(file.path(WGCNA_SCRIPT5, "tables", "loso", "loso_site_summary_by_module.csv")),
        csv_src(file.path(WGCNA_SCRIPT5, "tables", "downsampling", "balanced_downsampling_summary_by_module_trait.csv"))
      ),
      labels = c("LOSO site summary by module", "Balanced downsampling summary by module-trait")
    )
  ),

  tbl(
    42,
    "Categorical country and recruitment-site analyses of WGCNA module eigengenes.",
    "Clinical AD-associated plasma proteomic variation shows modular organization",
    "categorical recruitment-context analysis",
    combo_src(
      list(
        csv_src(file.path(WGCNA_SCRIPT4B, "tables", "country_categorical_kruskal_by_module.csv")),
        csv_src(file.path(WGCNA_SCRIPT4B, "tables", "country_categorical_adjusted_lm_by_module.csv")),
        csv_src(file.path(WGCNA_SCRIPT4B, "tables", "country_module_descriptives.csv")),
        csv_src(file.path(WGCNA_SCRIPT4B, "tables", "site_categorical_kruskal_by_module.csv")),
        csv_src(file.path(WGCNA_SCRIPT4B, "tables", "site_categorical_adjusted_lm_by_module.csv")),
        csv_src(file.path(WGCNA_SCRIPT4B, "tables", "site_module_descriptives.csv"))
      ),
      labels = c(
        "Country Kruskal-Wallis",
        "Country adjusted linear model",
        "Country descriptives",
        "Site Kruskal-Wallis",
        "Site adjusted linear model",
        "Site descriptives"
      )
    )
  )
)

###############################################################################
# 05. Build workbook
###############################################################################

wb <- openxlsx::createWorkbook()

index_rows <- list()
missing_rows <- list()

for (entry in table_plan) {
  message("Building Supplementary Table ", entry$number, ": ", entry$short_title)

  df <- tryCatch(
    read_table(entry$source),
    error = function(e) {
      missing_df(
        message = paste0("Table failed during read/combine step: ", conditionMessage(e)),
        attempted_sources = paste(capture.output(str(entry$source)), collapse = " ")
      )
    }
  )

  has_missing <- "status" %in% names(df) && any(df$status == "MISSING_SOURCE", na.rm = TRUE)
  if (has_missing) {
    missing_rows[[length(missing_rows) + 1]] <- tibble(
      Supplementary_Table = paste0("Supplementary Table ", entry$number),
      Title = entry$short_title,
      Manuscript_section = entry$manuscript_section,
      Missing_message = paste(unique(df$message), collapse = "; "),
      Attempted_sources = paste(unique(df$attempted_sources), collapse = " | ")
    )
  }

  write_table_sheet(
    wb = wb,
    sheet_name = entry$sheet,
    df = df,
    title = entry$title,
    note = entry$note
  )

  index_rows[[length(index_rows) + 1]] <- tibble(
    Supplementary_Table = paste0("Supplementary Table ", entry$number),
    Worksheet = clean_sheet_name(entry$sheet),
    Title = entry$short_title,
    Manuscript_section = entry$manuscript_section,
    Citation_hint = entry$citation_hint,
    Status = ifelse(has_missing, "MISSING SOURCE - check sheet", "OK")
  )
}

index_df <- bind_rows(index_rows)

# Add index sheet last, then move it first.
write_table_sheet(
  wb = wb,
  sheet_name = "Index",
  df = index_df,
  title = "Index of final Supplementary Tables",
  note = "Tables are ordered by first appearance in the manuscript. Machine-learning outputs are intentionally excluded from this workbook."
)

# Move Index to first position if supported.
try({
  openxlsx::worksheetOrder(wb) <- c(
    which(names(wb) == "Index"),
    setdiff(seq_along(names(wb)), which(names(wb) == "Index"))
  )
}, silent = TRUE)

openxlsx::saveWorkbook(wb, OUTPUT_XLSX, overwrite = TRUE)

readr::write_csv(index_df, OUTPUT_INDEX)

missing_df_final <- if (length(missing_rows) > 0) {
  bind_rows(missing_rows)
} else {
  tibble(
    Supplementary_Table = character(),
    Title = character(),
    Manuscript_section = character(),
    Missing_message = character(),
    Attempted_sources = character()
  )
}

readr::write_csv(missing_df_final, OUTPUT_MISSING)

###############################################################################
# 06. Manuscript citation helper
###############################################################################

citation_helper <- index_df %>%
  group_by(Manuscript_section) %>%
  summarise(
    Suggested_citations = paste(Supplementary_Table, collapse = "; "),
    .groups = "drop"
  )

readr::write_csv(
  citation_helper,
  file.path(OUT_DIR, "Supplementary_Tables_FINAL_Citation_Helper.csv")
)

###############################################################################
# 07. Final messages
###############################################################################

message("\nDone.")
message("Final workbook: ", OUTPUT_XLSX)
message("Index: ", OUTPUT_INDEX)
message("Missing sources report: ", OUTPUT_MISSING)

if (nrow(missing_df_final) > 0) {
  message("\nWARNING: Some tables have missing source files. See: ", OUTPUT_MISSING)
} else {
  message("\nAll source files were found.")
}

###############################################################################
# END
###############################################################################


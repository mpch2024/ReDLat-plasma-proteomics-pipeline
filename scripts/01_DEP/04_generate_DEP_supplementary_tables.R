###############################################################################
# ReDLat plasma proteomics
# 04_DEP_supplementary_tables_NatureAging_FINAL.R
#
# Purpose
# - Compile DEP/Figure 2 supplementary tables directly from pipeline CSV outputs.
# - Generate a clean Nature Aging-style supplementary table workbook.
# - Generate a separate Figure 2 source-data workbook.
# - Archive legacy / non-main tables without mixing them into the main supplement.
#
# Final outputs:
#   1. Supplementary_Tables_DEP_NatureAging_CLEAN.xlsx
#   2. SourceData_Figure2_DEP_NatureAging_CLEAN.xlsx
#   3. Archive_Legacy_DEP_Tables_NOT_FOR_MAIN_SUPPLEMENT.xlsx
#   4. Supplementary_Tables_DEP_NatureAging_Index.csv
###############################################################################

rm(list = ls())

###############################################################################
# 00. Packages
###############################################################################

packages <- c(
  "readr", "openxlsx", "dplyr", "stringr", "purrr", "tibble"
)

installed <- rownames(installed.packages())

for (p in packages) {
  if (!p %in% installed) install.packages(p)
}

suppressPackageStartupMessages({
  library(readr)
  library(openxlsx)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(tibble)
})

###############################################################################
# 01. Paths
###############################################################################

project_root <- "C:/Users/mnpiz/Desktop/DEPs_Proteomic_Publishable_V2"

out_dir <- file.path(
  project_root,
  "result",
  "final_tables",
  "NatureAging_clean"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

output_supp_xlsx <- file.path(
  out_dir,
  "Supplementary_Tables_DEP_NatureAging_CLEAN.xlsx"
)

output_source_xlsx <- file.path(
  out_dir,
  "SourceData_Figure2_DEP_NatureAging_CLEAN.xlsx"
)

output_archive_xlsx <- file.path(
  out_dir,
  "Archive_Legacy_DEP_Tables_NOT_FOR_MAIN_SUPPLEMENT.xlsx"
)

output_index_csv <- file.path(
  out_dir,
  "Supplementary_Tables_DEP_NatureAging_Index.csv"
)

gsea_dir <- file.path(project_root, "result", "05_enrichment_corrected", "gsea")
ora_dir  <- file.path(project_root, "result", "05_enrichment_corrected", "ora")

###############################################################################
# 02. Helper functions
###############################################################################

find_one <- function(dir, pattern) {
  if (!dir.exists(dir)) return(NA_character_)
  
  files <- list.files(
    dir,
    pattern = pattern,
    full.names = TRUE,
    recursive = TRUE
  )
  
  if (length(files) == 0) return(NA_character_)
  
  files[1]
}

first_existing_file <- function(paths) {
  hit <- paths[file.exists(paths)][1]
  
  if (length(hit) == 0 || is.na(hit)) {
    return(NA_character_)
  }
  
  hit
}

safe_read_csv <- function(path) {
  if (is.na(path) || !file.exists(path)) {
    return(NULL)
  }
  
  x <- readr::read_csv(
    path,
    show_col_types = FALSE,
    guess_max = 10000
  )
  
  x <- as.data.frame(x)
  
  return(x)
}

is_empty_table <- function(df) {
  if (is.null(df)) return(TRUE)
  if (nrow(df) == 0) return(TRUE)
  
  non_empty <- df %>%
    mutate(across(everything(), ~ as.character(.x))) %>%
    mutate(across(everything(), ~ ifelse(is.na(.x), "", .x)))
  
  all(apply(non_empty, 1, function(z) all(z == "")))
}

clean_colnames <- function(df) {
  if (is.null(df)) return(df)
  
  new_names <- names(df) %>%
    stringr::str_replace_all("\\s+", "_") %>%
    stringr::str_replace_all("[^A-Za-z0-9_\\.]", "_") %>%
    stringr::str_replace_all("_+", "_") %>%
    stringr::str_replace_all("_$", "")
  
  new_names[new_names == "" | is.na(new_names)] <- "Column"
  
  # openxlsx::writeDataTable() requires unique column names.
  names(df) <- make.unique(new_names, sep = "_dup")
  
  df
}

standardize_table <- function(df) {
  if (is.null(df)) return(NULL)
  
  df <- as.data.frame(df)
  df <- clean_colnames(df)
  
  df
}

as_character_df <- function(df) {
  if (is.null(df)) return(NULL)
  
  df <- as.data.frame(df)
  df <- standardize_table(df)
  
  df[] <- lapply(df, function(x) {
    if (inherits(x, "Date")) {
      as.character(x)
    } else if (inherits(x, "POSIXct") || inherits(x, "POSIXt")) {
      as.character(x)
    } else {
      as.character(x)
    }
  })
  
  df
}

shorten_sheet_name <- function(x) {
  x <- stringr::str_replace_all(x, "[\\[\\]\\*\\?/\\\\:]", "_")
  substr(x, 1, 31)
}

write_df_sheet <- function(wb, sheet_name, df, title = NULL, note = NULL) {
  
  sheet_name <- shorten_sheet_name(sheet_name)
  
  if (sheet_name %in% names(wb)) {
    removeWorksheet(wb, sheet_name)
  }
  
  addWorksheet(wb, sheet_name)
  
  row_start <- 1
  
  if (!is.null(title)) {
    writeData(
      wb,
      sheet_name,
      title,
      startRow = row_start,
      startCol = 1
    )
    
    addStyle(
      wb,
      sheet_name,
      style = createStyle(
        fontName = "Arial",
        textDecoration = "bold",
        fontSize = 12
      ),
      rows = row_start,
      cols = 1
    )
    
    row_start <- row_start + 2
  }
  
  if (!is.null(note)) {
    writeData(
      wb,
      sheet_name,
      note,
      startRow = row_start,
      startCol = 1
    )
    
    addStyle(
      wb,
      sheet_name,
      style = createStyle(
        fontName = "Arial",
        textDecoration = "italic",
        fontSize = 9,
        fontColour = "#555555",
        wrapText = TRUE
      ),
      rows = row_start,
      cols = 1
    )
    
    row_start <- row_start + 2
  }
  
  if (is.null(df) || nrow(df) == 0 || ncol(df) == 0) {
    writeData(
      wb,
      sheet_name,
      "No data available or source file was not found.",
      startRow = row_start,
      startCol = 1
    )
    return(invisible(wb))
  }
  
  df <- standardize_table(df)
  
  if (anyDuplicated(names(df)) > 0) {
    names(df) <- make.unique(names(df), sep = "_dup")
  }
  
  tryCatch(
    {
      writeDataTable(
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
      warning(
        "writeDataTable failed for sheet ",
        sheet_name,
        ". Writing as plain data instead. Error: ",
        conditionMessage(e)
      )
      
      writeData(
        wb,
        sheet = sheet_name,
        x = df,
        startRow = row_start,
        startCol = 1,
        withFilter = TRUE
      )
    }
  )
  
  header_style <- createStyle(
    fontName = "Arial",
    fontSize = 10,
    fontColour = "white",
    fgFill = "#1F4E79",
    halign = "center",
    valign = "center",
    textDecoration = "bold",
    border = "Bottom",
    borderColour = "#D9EAF7"
  )
  
  body_style <- createStyle(
    fontName = "Arial",
    fontSize = 9,
    valign = "top",
    wrapText = TRUE
  )
  
  addStyle(
    wb,
    sheet_name,
    header_style,
    rows = row_start,
    cols = seq_len(ncol(df)),
    gridExpand = TRUE,
    stack = TRUE
  )
  
  if (nrow(df) > 0) {
    addStyle(
      wb,
      sheet_name,
      body_style,
      rows = seq(row_start + 1, row_start + nrow(df)),
      cols = seq_len(ncol(df)),
      gridExpand = TRUE,
      stack = TRUE
    )
  }
  
  freezePane(wb, sheet_name, firstActiveRow = row_start + 1)
  
  for (j in seq_len(ncol(df))) {
    col_values <- as.character(df[[j]])
    max_chars <- suppressWarnings(max(nchar(col_values), na.rm = TRUE))
    
    if (is.infinite(max_chars)) max_chars <- 10
    
    width <- min(max(max_chars + 2, 10), 38)
    
    setColWidths(
      wb,
      sheet_name,
      cols = j,
      widths = width
    )
  }
  
  invisible(wb)
}

write_index_sheet <- function(wb, index_df, sheet_name = "Index") {
  
  sheet_name <- shorten_sheet_name(sheet_name)
  
  if (sheet_name %in% names(wb)) {
    removeWorksheet(wb, sheet_name)
  }
  
  addWorksheet(wb, sheet_name)
  
  writeData(
    wb,
    sheet_name,
    "Supplementary Tables Index",
    startRow = 1,
    startCol = 1
  )
  
  addStyle(
    wb,
    sheet_name,
    createStyle(
      fontName = "Arial",
      textDecoration = "bold",
      fontSize = 14
    ),
    rows = 1,
    cols = 1
  )
  
  index_df <- standardize_table(index_df)
  
  writeDataTable(
    wb,
    sheet = sheet_name,
    x = index_df,
    startRow = 3,
    startCol = 1,
    tableStyle = "TableStyleMedium2",
    withFilter = TRUE
  )
  
  freezePane(wb, sheet_name, firstActiveRow = 4)
  
  setColWidths(
    wb,
    sheet_name,
    cols = 1:ncol(index_df),
    widths = "auto"
  )
  
  invisible(wb)
}

###############################################################################
# 03. Source file registry
###############################################################################

source_registry <- c(
  # Main DEP -------------------------------------------------------------------
  "Table_S1_DEP_main_full" =
    file.path(
      project_root,
      "result/03_dep/gene_collapsed/AD_vs_CN_full_limma_results_gene_collapsed.csv"
    ),
  
  "Table_S2_DEP_main_FDR005" =
    file.path(
      project_root,
      "result/03_dep/gene_collapsed/FDR_specific/AD_vs_CN_DEP_gene_collapsed_FDR005.csv"
    ),
  
  "Table_S3_DEP_main_FDR001" =
    file.path(
      project_root,
      "result/03_dep/gene_collapsed/FDR_specific/AD_vs_CN_DEP_gene_collapsed_FDR001.csv"
    ),
  
  # Enrichment -----------------------------------------------------------------
  "Table_S4_GSEA_GO" =
    file.path(
      gsea_dir,
      "main_dep_gsea_go_bh.csv"
    ),
  
  "Table_S5_GSEA_KEGG" =
    file.path(
      gsea_dir,
      "main_dep_gsea_kegg_bh.csv"
    ),
  
  "Table_S6_GSEA_Reactome" =
    file.path(
      gsea_dir,
      "main_dep_gsea_reactome_bh.csv"
    ),
  
  "Table_S7_GSEA_Hallmark" =
    file.path(
      gsea_dir,
      "main_dep_gsea_hallmark_bh.csv"
    ),
  
  "Table_S8_ORA_Up_GO" =
    find_one(
      ora_dir,
      "higher_in_AD_GO_ORA_BH\\.csv$"
    ),
  
  "Table_S9_ORA_Up_KEGG" =
    find_one(
      ora_dir,
      "higher_in_AD_KEGG_ORA_BH\\.csv$"
    ),
  
  "Table_S10_ORA_Up_Reactome" =
    find_one(
      ora_dir,
      "higher_in_AD_Reactome_ORA_BH\\.csv$"
    ),
  
  "Table_S11_ORA_Down_GO" =
    find_one(
      ora_dir,
      "lower_in_AD_GO_ORA_BH\\.csv$"
    ),
  
  "Table_S12_ORA_Down_KEGG" =
    find_one(
      ora_dir,
      "lower_in_AD_KEGG_ORA_BH\\.csv$"
    ),
  
  "Table_S13_ORA_Down_Reactome" =
    find_one(
      ora_dir,
      "lower_in_AD_Reactome_ORA_BH\\.csv$"
    ),
  
  # AD-only CDR-SB --------------------------------------------------------------
  "Table_S14_AD_only_CDRSB_alignment" =
    file.path(
      project_root,
      "result/04_sensitivity/cdrsb/AD_only/primary_AD_vs_CN_vs_AD_only_CDRSB_severity_alignment.csv"
    ),
  
  "Table_S15_AD_only_CDRSB_full" =
    file.path(
      project_root,
      "result/04_sensitivity/cdrsb/AD_only/AD_only_CDRSB_severity_full_limma_results_gene_collapsed.csv"
    ),
  
  "Table_S16_AD_only_CDRSB_counts" =
    file.path(
      project_root,
      "result/04_sensitivity/cdrsb/AD_only/AD_only_CDRSB_severity_counts_gene_collapsed.csv"
    ),
  
  "Table_S17_AD_only_CDRSB_Reactome" =
    file.path(
      gsea_dir,
      "ad_only_cdrsb_severity_gsea_reactome_bh.csv"
    ),
  
  # Legacy CDR-SB-adjusted AD vs CN --------------------------------------------
  "Table_S18_CDRSBadj_ADvsCN_comparison" =
    file.path(
      project_root,
      "result/04_sensitivity/cdrsb/AD_vs_CN_adjusted/primary_vs_CDRSB_adjusted_AD_vs_CN_gene_comparison.csv"
    ),
  
  "Table_S19_CDRSBadj_ADvsCN_full" =
    file.path(
      project_root,
      "result/04_sensitivity/cdrsb/AD_vs_CN_adjusted/AD_vs_CN_CDRSB_adjusted_full_limma_results_gene_collapsed.csv"
    ),
  
  "Table_S20_CDRSBadj_ADvsCN_count" =
    file.path(
      project_root,
      "result/04_sensitivity/cdrsb/AD_vs_CN_adjusted/CDRSB_adjusted_AD_vs_CN_DEP_counts_gene_collapsed.csv"
    ),
  
  "Table_S21_CDRSBadj_ADvsCN_Reactome" =
    file.path(
      gsea_dir,
      "cdrsb_adjusted_ad_vs_cn_gsea_reactome_bh.csv"
    ),
  
  # Additional sensitivity ------------------------------------------------------
  "Table_S22_Clinical_severity" =
    file.path(
      project_root,
      "result/04_sensitivity/clinical_severity/severity_assoc_all_outcomes_combined.csv"
    ),
  
  "Table_S23_APOE_sensitivity" =
    file.path(
      project_root,
      "result/04_sensitivity/apoe/primary_vs_APOE_adjusted_gene_comparison.csv"
    ),
  
  "Table_S24_APOE_adjusted_Reactome" =
    file.path(
      gsea_dir,
      "apoe_adjusted_dep_gsea_reactome_bh.csv"
    ),
  
  "Table_S25_ATN_adjusted_compare" =
    file.path(
      project_root,
      "result/04_sensitivity/atn_adjusted/primary_vs_ATN_adjusted_gene_comparison.csv"
    ),
  
  "Table_S26_ATN_adjusted_full" =
    file.path(
      project_root,
      "result/04_sensitivity/atn_adjusted/AD_vs_CN_ATN_adjusted_full_limma_results_gene_collapsed.csv"
    ),
  
  "Table_S27_ATN_adjusted_counts" =
    file.path(
      project_root,
      "result/04_sensitivity/atn_adjusted/ATN_adjusted_DEP_counts_gene_collapsed.csv"
    ),
  
  "Table_S28_ATN_covariates_used" =
    file.path(
      project_root,
      "result/04_sensitivity/atn_adjusted/ATN_covariates_used.csv"
    ),
  
  "Table_S29_ATN_adjusted_Reactome" =
    file.path(
      gsea_dir,
      "atn_adjusted_dep_gsea_reactome_bh.csv"
    ),
  
  # Figure 2 source data --------------------------------------------------------
  "Table_S30_Fig2a_volcano_source" =
    file.path(
      project_root,
      "result/final_source_data/Figure2/Figure2a_main_volcano_source_data.csv"
    ),
  
  "Table_S31_Fig2b_Reactome_source" =
    file.path(
      project_root,
      "result/final_source_data/Figure2/Figure2b_reactome_gsea_source_data.csv"
    ),
  
  "Table_S32_Fig2c_heatmap_source" =
    file.path(
      project_root,
      "result/final_source_data/Figure2/Figure2c_representative_heatmap_source_data.csv"
    ),
  
  "Table_S33_Fig2c_selected_proteins" =
    file.path(
      project_root,
      "result/final_source_data/Figure2/Figure2c_representative_heatmap_selected_proteins.csv"
    ),
  
  "Table_S34_Fig2d_ATN_source" =
    file.path(
      project_root,
      "result/final_source_data/Figure2/Figure2d_ATN_adjusted_attenuation_source_data.csv"
    ),
  
  "Table_S35_Fig2e_country_source" =
    file.path(
      project_root,
      "result/final_source_data/Figure2/Figure2e_contextual_country_exclusion_source_data.csv"
    ),
  
  "Table_S36_Fig2f_LOCO_source" =
    file.path(
      project_root,
      "result/final_source_data/Figure2/Figure2f_LOCO_stability_source_data.csv"
    ),
  
  "Table_S37_Supp_ADonly_CDRSB_fig" =
    file.path(
      project_root,
      "result/final_source_data/Figure2/Supplementary_AD_only_CDRSB_severity_alignment_source_data.csv"
    ),
  
  "Table_S38_Supp_CDRSBadj_fig" =
    file.path(
      project_root,
      "result/final_source_data/Supplementary/Supplementary_CDRSB_adjusted_AD_vs_CN_attenuation_source_data.csv"
    ),
  
  # Country/site robustness -----------------------------------------------------
  "Table_S39_LOCO_summary" =
    file.path(
      project_root,
      "result/06_robustness/country_loco/tables/LOCO_summary_metrics.csv"
    ),
  
  "Table_S40_LOCO_proteins" =
    file.path(
      project_root,
      "result/06_robustness/country_loco/tables/main_vs_meanLOCO_table.csv"
    ),
  
  "Table_S41_Country_meta" =
    file.path(
      project_root,
      "result/06_robustness/country_meta/tables/country_meta_analysis_results.csv"
    ),
  
  "Table_S42_Country_group_counts" =
    file.path(
      project_root,
      "result/06_robustness/country_loco/tables/country_group_counts.csv"
    ),
  
  "Table_S43_Balanced_resampling_summary" =
    file.path(
      project_root,
      "result/06_robustness/balanced_country_resampling/tables/balanced_resampling_summary_metrics.csv"
    ),
  
  "Table_S45_LOSO_site" =
    file.path(
      project_root,
      "result/06_robustness/site_robustness/loso/tables/LOSO_summary_metrics.csv"
    ),
  
  # Classification outputs ------------------------------------------------------
  "Table_S49_Robustness_classification" =
    first_existing_file(c(
      file.path(
        project_root,
        "result/supplementary/tables/robustness_classification/main_DEP_FDR005_extended_robustness_classification.csv"
      ),
      file.path(
        project_root,
        "result/06_robustness/formal_classification/protein_robustness_classification.csv"
      )
    )),
  
  "Table_S50_Robustness_counts" =
    first_existing_file(c(
      file.path(
        project_root,
        "result/supplementary/tables/robustness_classification/main_DEP_FDR005_extended_robustness_classification_counts.csv"
      ),
      file.path(
        project_root,
        "result/06_robustness/formal_classification/protein_robustness_classification_counts.csv"
      )
    ))
)

###############################################################################
# 04. File check
###############################################################################

file_check <- tibble::tibble(
  source_table = names(source_registry),
  path = as.character(source_registry),
  exists = file.exists(as.character(source_registry)),
  status = dplyr::if_else(exists, "FOUND", "MISSING")
)

###############################################################################
# 05. Table construction helpers
###############################################################################

read_source <- function(source_id) {
  if (!source_id %in% names(source_registry)) {
    warning("Source id not found in registry: ", source_id)
    return(NULL)
  }
  
  path <- as.character(source_registry[[source_id]])
  
  if (is.na(path) || !file.exists(path)) {
    warning("Source file missing for ", source_id, ": ", path)
    return(NULL)
  }
  
  message("Reading source: ", source_id)
  
  safe_read_csv(path)
}

read_source_character <- function(source_id) {
  df <- read_source(source_id)
  as_character_df(df)
}

combine_ora_tables <- function(source_ids) {
  
  out <- purrr::map_dfr(source_ids, function(source_id) {
    
    df <- read_source_character(source_id)
    
    database <- dplyr::case_when(
      stringr::str_detect(source_id, "GO") ~ "GO",
      stringr::str_detect(source_id, "KEGG") ~ "KEGG",
      stringr::str_detect(source_id, "Reactome") ~ "Reactome",
      TRUE ~ NA_character_
    )
    
    direction <- dplyr::case_when(
      stringr::str_detect(source_id, "_Up_") ~ "Higher in AD",
      stringr::str_detect(source_id, "_Down_") ~ "Lower in AD",
      TRUE ~ NA_character_
    )
    
    if (is.null(df) || is_empty_table(df)) {
      return(
        tibble(
          source_table = source_id,
          database = database,
          direction = direction,
          note = "No significant terms available in this category."
        )
      )
    }
    
    df %>%
      mutate(
        source_table = source_id,
        database = database,
        direction = direction,
        .before = 1
      )
  })
  
  out
}

combine_tagged_tables <- function(source_ids) {
  
  out <- purrr::map_dfr(source_ids, function(source_id) {
    
    df <- read_source_character(source_id)
    
    if (is.null(df) || is_empty_table(df)) {
      return(
        tibble(
          source_table = source_id,
          note = "No data available or source file was not found."
        )
      )
    }
    
    df %>%
      mutate(source_table = source_id, .before = 1)
  })
  
  out
}

build_table <- function(action, source_ids) {
  
  source_ids <- unlist(strsplit(source_ids, ";", fixed = TRUE))
  
  if (action == "copy") {
    return(read_source(source_ids[1]))
  }
  
  if (action == "combine_ora") {
    return(combine_ora_tables(source_ids))
  }
  
  if (action == "combine_tagged") {
    return(combine_tagged_tables(source_ids))
  }
  
  stop("Unknown action: ", action)
}

###############################################################################
# 06. Clean Supplementary Tables plan
###############################################################################

supplementary_plan <- tibble::tribble(
  ~new_id, ~new_sheet, ~title, ~source_ids, ~action,
  
  "S1", "DEP_full",
  "Full gene-collapsed differential abundance results for AD versus CN.",
  "Table_S1_DEP_main_full", "copy",
  
  "S2", "DEP_FDR005",
  "Differentially abundant proteins associated with AD at FDR < 0.05.",
  "Table_S2_DEP_main_FDR005", "copy",
  
  "S3", "DEP_FDR001",
  "High-confidence differentially abundant proteins associated with AD at FDR < 0.01.",
  "Table_S3_DEP_main_FDR001", "copy",
  
  "S4", "GSEA_GO",
  "Pre-ranked Gene Ontology gene set enrichment analysis of the AD versus CN proteomic signature.",
  "Table_S4_GSEA_GO", "copy",
  
  "S5", "GSEA_KEGG",
  "Pre-ranked KEGG gene set enrichment analysis of the AD versus CN proteomic signature.",
  "Table_S5_GSEA_KEGG", "copy",
  
  "S6", "GSEA_Reactome",
  "Pre-ranked Reactome gene set enrichment analysis of the AD versus CN proteomic signature.",
  "Table_S6_GSEA_Reactome", "copy",
  
  "S7", "GSEA_Hallmark",
  "Pre-ranked Hallmark gene set enrichment analysis of the AD versus CN proteomic signature.",
  "Table_S7_GSEA_Hallmark", "copy",
  
  "S8", "Directional_ORA",
  "Direction-specific over-representation analysis of proteins with higher or lower abundance in AD.",
  "Table_S8_ORA_Up_GO;Table_S9_ORA_Up_KEGG;Table_S10_ORA_Up_Reactome;Table_S11_ORA_Down_GO;Table_S12_ORA_Down_KEGG;Table_S13_ORA_Down_Reactome",
  "combine_ora",
  
  "S9", "ADonly_CDRSB_alignment",
  "Alignment between the primary diagnostic AD versus CN proteomic signature and the AD-only CDR-SB severity model.",
  "Table_S14_AD_only_CDRSB_alignment", "copy",
  
  "S10", "ADonly_CDRSB_full",
  "Full AD-only CDR-SB severity model results.",
  "Table_S15_AD_only_CDRSB_full", "copy",
  
  "S11", "ADonly_CDRSB_summary",
  "Summary counts and Reactome enrichment for the AD-only CDR-SB severity analysis.",
  "Table_S16_AD_only_CDRSB_counts;Table_S17_AD_only_CDRSB_Reactome",
  "combine_tagged",
  
  "S12", "APOE_sensitivity",
  "APOE ε4 sensitivity analysis comparing the primary diagnostic model with the APOE-adjusted model.",
  "Table_S23_APOE_sensitivity", "copy",
  
  "S13", "APOE_Reactome",
  "Reactome enrichment analysis for the APOE-adjusted diagnostic model.",
  "Table_S24_APOE_adjusted_Reactome", "copy",
  
  "S14", "ATN_adjusted",
  "AT(N)-adjusted diagnostic model comparison and full protein-level results.",
  "Table_S25_ATN_adjusted_compare;Table_S26_ATN_adjusted_full;Table_S27_ATN_adjusted_counts;Table_S28_ATN_covariates_used",
  "combine_tagged",
  
  "S15", "ATN_Reactome",
  "Reactome enrichment analysis for the AT(N)-adjusted diagnostic model.",
  "Table_S29_ATN_adjusted_Reactome", "copy",
  
  "S16", "LOCO_country",
  "Leave-one-country-out robustness summary and protein-level preservation metrics.",
  "Table_S39_LOCO_summary;Table_S40_LOCO_proteins",
  "combine_tagged",
  
  "S17", "Country_meta",
  "Country-level meta-analysis of AD-associated proteomic effects and sample counts.",
  "Table_S41_Country_meta;Table_S42_Country_group_counts",
  "combine_tagged",
  
  "S18", "Balanced_resampling",
  "Balanced country-resampling robustness summary.",
  "Table_S43_Balanced_resampling_summary", "copy",
  
  "S19", "LOSO_site",
  "Leave-one-site-out site-level robustness analysis.",
  "Table_S45_LOSO_site", "copy"
)

###############################################################################
# 07. Source Data and Archive plans
###############################################################################

source_plan <- tibble::tribble(
  ~new_sheet, ~title, ~source_id,
  
  "SourceData_Fig2a",
  "Source data for Fig. 2a volcano plot.",
  "Table_S30_Fig2a_volcano_source",
  
  "SourceData_Fig2b",
  "Source data for Fig. 2b Reactome enrichment panel.",
  "Table_S31_Fig2b_Reactome_source",
  
  "SourceData_Fig2c",
  "Source data for Fig. 2c protein-trait heatmap.",
  "Table_S32_Fig2c_heatmap_source",
  
  "SourceData_Fig2c_selected",
  "Selected proteins displayed in Fig. 2c.",
  "Table_S33_Fig2c_selected_proteins",
  
  "SourceData_Fig2d",
  "Source data for Fig. 2d AT(N)-adjusted attenuation panel.",
  "Table_S34_Fig2d_ATN_source",
  
  "SourceData_Fig2e",
  "Source data for Fig. 2e country-exclusion contextual stability panel.",
  "Table_S35_Fig2e_country_source",
  
  "SourceData_Fig2f",
  "Source data for Fig. 2f LOCO stability panel.",
  "Table_S36_Fig2f_LOCO_source",
  
  "SourceData_Supp_ADonly_CDRSB",
  "Source data for supplementary AD-only CDR-SB severity figure.",
  "Table_S37_Supp_ADonly_CDRSB_fig",
  
  "SourceData_Supp_CDRSBadj_legacy",
  "Source data for legacy CDR-SB-adjusted AD versus CN figure.",
  "Table_S38_Supp_CDRSBadj_fig"
)

archive_plan <- tibble::tribble(
  ~new_sheet, ~title, ~source_id,
  
  "Archive_CDRSBadj_compare",
  "Legacy CDR-SB-adjusted AD versus CN comparison. Not recommended for the main supplement.",
  "Table_S18_CDRSBadj_ADvsCN_comparison",
  
  "Archive_CDRSBadj_full",
  "Legacy CDR-SB-adjusted AD versus CN full results. Not recommended for the main supplement.",
  "Table_S19_CDRSBadj_ADvsCN_full",
  
  "Archive_CDRSBadj_counts",
  "Legacy CDR-SB-adjusted AD versus CN counts. Not recommended for the main supplement.",
  "Table_S20_CDRSBadj_ADvsCN_count",
  
  "Archive_CDRSBadj_Reactome",
  "Legacy CDR-SB-adjusted AD versus CN Reactome enrichment. Not recommended for the main supplement.",
  "Table_S21_CDRSBadj_ADvsCN_Reactome",
  
  "Archive_Clinical_severity_large",
  "Large clinical severity output. Recommended as archival/source data rather than main supplementary table.",
  "Table_S22_Clinical_severity",
  
  "Archive_Classif_robustness",
  "Classification robustness table. Should be moved to ML/Fig. 4–5 supplement if retained.",
  "Table_S49_Robustness_classification",
  
  "Archive_Classif_counts",
  "Classification robustness counts. Should be moved to ML/Fig. 4–5 supplement if retained.",
  "Table_S50_Robustness_counts"
)

###############################################################################
# 08. Build clean Supplementary Tables workbook
###############################################################################

wb_supp <- createWorkbook()
supp_index_rows <- list()

for (i in seq_len(nrow(supplementary_plan))) {
  
  plan_i <- supplementary_plan[i, ]
  
  message("Building Supplementary Table ", plan_i$new_id, ": ", plan_i$new_sheet)
  
  df <- build_table(
    action = plan_i$action,
    source_ids = plan_i$source_ids
  )
  
  final_sheet_name <- paste0(plan_i$new_id, "_", plan_i$new_sheet)
  
  note_i <- paste0(
    "Generated from: ",
    plan_i$source_ids,
    ". This table is part of the cleaned DEP supplementary table set."
  )
  
  write_df_sheet(
    wb = wb_supp,
    sheet_name = final_sheet_name,
    df = df,
    title = paste0("Supplementary Table ", plan_i$new_id, ". ", plan_i$title),
    note = note_i
  )
  
  supp_index_rows[[i]] <- tibble(
    Supplementary_Table = plan_i$new_id,
    Sheet = shorten_sheet_name(final_sheet_name),
    Title = plan_i$title,
    Source_tables = plan_i$source_ids,
    Status = "Retained in clean Supplementary Tables"
  )
}

supp_index_df <- bind_rows(supp_index_rows)

write_index_sheet(wb_supp, supp_index_df, "Index")

saveWorkbook(
  wb_supp,
  output_supp_xlsx,
  overwrite = TRUE
)

###############################################################################
# 09. Build Source Data workbook
###############################################################################

wb_source <- createWorkbook()

source_index <- source_plan %>%
  mutate(
    Path = as.character(source_registry[source_id]),
    Available = file.exists(Path),
    Status = "Moved to Source Data workbook"
  )

for (i in seq_len(nrow(source_plan))) {
  
  plan_i <- source_plan[i, ]
  
  message("Building source data: ", plan_i$new_sheet)
  
  df <- read_source(plan_i$source_id)
  
  write_df_sheet(
    wb = wb_source,
    sheet_name = plan_i$new_sheet,
    df = df,
    title = plan_i$title,
    note = paste0("Original source table: ", plan_i$source_id)
  )
}

write_index_sheet(wb_source, source_index, "Index")

saveWorkbook(
  wb_source,
  output_source_xlsx,
  overwrite = TRUE
)

###############################################################################
# 10. Build Archive workbook
###############################################################################

wb_archive <- createWorkbook()

archive_index <- archive_plan %>%
  mutate(
    Path = as.character(source_registry[source_id]),
    Available = file.exists(Path),
    Status = "Archived / not recommended for main Supplementary Tables"
  )

for (i in seq_len(nrow(archive_plan))) {
  
  plan_i <- archive_plan[i, ]
  
  message("Building archive table: ", plan_i$new_sheet)
  
  df <- read_source(plan_i$source_id)
  
  write_df_sheet(
    wb = wb_archive,
    sheet_name = plan_i$new_sheet,
    df = df,
    title = plan_i$title,
    note = paste0("Original source table: ", plan_i$source_id)
  )
}

write_index_sheet(wb_archive, archive_index, "Index")

saveWorkbook(
  wb_archive,
  output_archive_xlsx,
  overwrite = TRUE
)

###############################################################################
# 11. Export global audit index
###############################################################################

global_audit <- bind_rows(
  file_check %>%
    transmute(
      Workbook = "Input file registry",
      Table_or_sheet = source_table,
      Title = NA_character_,
      Source = path,
      Status = status
    ),
  supp_index_df %>%
    transmute(
      Workbook = "Supplementary_Tables_DEP_NatureAging_CLEAN.xlsx",
      Table_or_sheet = Sheet,
      Title = Title,
      Source = Source_tables,
      Status = Status
    ),
  source_index %>%
    transmute(
      Workbook = "SourceData_Figure2_DEP_NatureAging_CLEAN.xlsx",
      Table_or_sheet = new_sheet,
      Title = title,
      Source = source_id,
      Status = Status
    ),
  archive_index %>%
    transmute(
      Workbook = "Archive_Legacy_DEP_Tables_NOT_FOR_MAIN_SUPPLEMENT.xlsx",
      Table_or_sheet = new_sheet,
      Title = title,
      Source = source_id,
      Status = Status
    )
)

write.csv(
  global_audit,
  output_index_csv,
  row.names = FALSE
)

###############################################################################
# 12. Console summary
###############################################################################

message("\n====================================================")
message("Clean DEP supplementary table files generated:")
message("1) ", output_supp_xlsx)
message("2) ", output_source_xlsx)
message("3) ", output_archive_xlsx)
message("4) ", output_index_csv)
message("====================================================\n")

message("Recommended manuscript citation ranges after cleaning:")
message("Main DEP: Supplementary Tables S1-S3")
message("Directional enrichment: Supplementary Tables S4-S8")
message("AD-only CDR-SB severity: Supplementary Tables S9-S11")
message("APOE sensitivity: Supplementary Tables S12-S13")
message("AT(N)-adjusted model: Supplementary Tables S14-S15")
message("Country/site robustness: Supplementary Tables S16-S19")

message("\nFile check summary:")
print(file_check %>% count(status, name = "n"))

message("\nMissing source files:")
print(file_check %>% filter(!exists))

###############################################################################
# END
###############################################################################


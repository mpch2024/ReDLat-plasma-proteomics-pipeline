###############################################################################
# ReDLat plasma proteomics
# 05_additional_tables_sensitivity_MANUSCRIPT_READY.R
#
# Purpose
# - Manuscript-ready supplementary tables and internal sensitivity analyses.
# - Generate ONE manuscript-style consolidated supplementary table with CN vs
#   clinically diagnosed AD composition overall, by country, and by site when
#   available.
# - Run internal sensitivity analyses:
#     1) Propensity-score matched/balanced DEP sensitivity.
#     2) Biomarker-consistent DEP sensitivity using plasma p-tau217 and Aβ42/40.
#
# Important interpretation
# - These are internal sensitivity/robustness analyses, not external validation.
# - Country/site are treated as categorical recruitment-context variables.
# - The biomarker-consistent analysis uses data-driven thresholds by default;
#   report them explicitly.
#
# Expected input
# - Run 01_data_processing_and_differential_analysis.R first.
# - This script loads: result/workspace/analysis_workspace.RData
#
# Main outputs
# - result/08_manuscript_supplementary_sensitivity_MATCHIT/Supplementary_Table_1_Diagnostic_demographic_composition_by_recruitment_context.xlsx
# - result/08_manuscript_supplementary_sensitivity_MATCHIT/Supplementary_Table_2_Internal_sensitivity_DEP_signature.xlsx
# - result/08_manuscript_supplementary_sensitivity_MATCHIT/*.csv
###############################################################################

rm(list = ls())

###############################################################################
# 00. Packages and paths
###############################################################################

packages <- c(
  "dplyr", "tidyr", "purrr", "tibble", "stringr", "readr", "openxlsx",
  "limma", "rlang", "stats", "MatchIt", "cobalt"
)

installed <- rownames(installed.packages())
for (p in packages) {
  if (!p %in% installed) install.packages(p)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(stringr)
  library(readr)
  library(openxlsx)
  library(limma)
  library(rlang)
  library(MatchIt)
  library(cobalt)
})

if (requireNamespace("here", quietly = TRUE)) {
  project_root <- here::here()
} else {
  project_root <- getwd()
}

workspace_file <- file.path(project_root, "result", "workspace", "analysis_workspace.RData")
if (!file.exists(workspace_file)) {
  stop(
    "Cannot find analysis workspace: ", workspace_file,
    "\nRun 01_data_processing_and_differential_analysis.R first, or set project_root manually."
  )
}

load(workspace_file)

out_dir <- file.path(project_root, "result", "08_manuscript_supplementary_sensitivity_MATCHIT")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

MAIN_FDR <- if (exists("MAIN_FDR")) MAIN_FDR else 0.05
STRICT_FDR <- if (exists("STRICT_FDR")) STRICT_FDR else 0.01
MAIN_GROUPS <- if (exists("MAIN_GROUPS")) MAIN_GROUPS else c("CN", "AD")
SEED_MANUSCRIPT <- 20260519
set.seed(SEED_MANUSCRIPT)

###############################################################################
# 01. Helper functions
###############################################################################

safe_write_csv <- function(x, file) {
  if (is.null(x)) return(invisible(NULL))
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(as_tibble(x), file)
}

safe_file_tag <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  ifelse(nchar(x) == 0, "NA", x)
}

safe_log2_matrix <- function(mat) {
  mat <- apply(mat, 2, as.numeric)
  mat[mat <= 0] <- NA_real_
  log2(mat)
}

safe_se_from_limma <- function(logFC, t_stat) {
  logFC <- suppressWarnings(as.numeric(logFC))
  t_stat <- suppressWarnings(as.numeric(t_stat))
  se <- rep(NA_real_, length(logFC))
  ok <- is.finite(logFC) & is.finite(t_stat) & abs(t_stat) > .Machine$double.eps
  se[ok] <- abs(logFC[ok] / t_stat[ok])
  se
}

clean_text_na <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA", "NaN", "NULL", "null", "N/A")] <- NA_character_
  x
}

format_p <- function(p) {
  p <- suppressWarnings(as.numeric(p))
  ifelse(
    is.na(p), "-",
    ifelse(p < 0.0001, "< 0.0001",
           ifelse(p < 0.001, "< 0.001", sprintf("%.3f", p)))
  )
}

mean_sd <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  n <- sum(!is.na(x))
  if (n == 0) return("-")
  paste0(sprintf("%.2f", mean(x, na.rm = TRUE)), " ± ", sprintf("%.2f", sd(x, na.rm = TRUE)))
}

n_percent <- function(n, denom) {
  if (is.na(denom) || denom == 0) return("-")
  paste0(n, " (", sprintf("%.1f", 100 * n / denom), ")")
}

n_over_n_percent <- function(n, denom) {
  if (is.na(denom) || denom == 0) return("-")
  paste0(n, "/", denom, " (", sprintf("%.1f", 100 * n / denom), ")")
}

as_numeric_safe <- function(x) suppressWarnings(as.numeric(x))

is_female_vec <- function(x) {
  xx <- tolower(trimws(as.character(x)))
  xx %in% c("f", "female", "woman", "women", "mujer", "2")
}

is_apoe4_carrier_vec <- function(x) {
  if (is.null(x)) return(rep(NA, 0))
  if (is.numeric(x) || is.integer(x)) return(as.numeric(x) == 1)
  xx <- tolower(trimws(as.character(x)))
  xx %in% c("1", "yes", "y", "true", "carrier", "e4", "ε4", "apoe4 carrier")
}

safe_wilcox_p <- function(x, g) {
  x <- as_numeric_safe(x)
  g <- as.character(g)
  ok <- !is.na(x) & !is.na(g) & g %in% c("CN", "AD")
  x <- x[ok]; g <- g[ok]
  if (!all(c("CN", "AD") %in% unique(g))) return(NA_real_)
  if (sum(g == "CN") < 2 || sum(g == "AD") < 2) return(NA_real_)
  tryCatch(wilcox.test(x ~ g)$p.value, error = function(e) NA_real_)
}

safe_p_chisq_or_fisher <- function(tab) {
  tab <- as.matrix(tab)
  if (any(dim(tab) < 2)) return(NA_real_)
  out <- tryCatch({
    cs <- suppressWarnings(chisq.test(tab))
    if (any(cs$expected < 5)) fisher.test(tab, simulate.p.value = TRUE, B = 10000)$p.value else cs$p.value
  }, error = function(e) NA_real_)
  out
}

safe_binary_p <- function(x_binary, g) {
  g <- as.character(g)
  ok <- !is.na(x_binary) & !is.na(g) & g %in% c("CN", "AD")
  if (sum(ok) == 0) return(NA_real_)
  safe_p_chisq_or_fisher(table(x_binary[ok], g[ok]))
}

sanitize_df <- function(df) {
  if (is.null(df)) return(NULL)
  df <- as.data.frame(df, check.names = FALSE, stringsAsFactors = FALSE)
  
  # Robust Excel-safe names. Some outputs can contain duplicated names after joins
  # or names that differ only by punctuation/case. openxlsx::writeDataTable() is
  # strict, so we sanitize names before export and use writeData() below.
  nm <- names(df)
  nm[is.na(nm) | trimws(nm) == ""] <- paste0("column_", which(is.na(nm) | trimws(nm) == ""))
  nm <- gsub("[^A-Za-z0-9_]+", "_", nm)
  nm <- gsub("^_+|_+$", "", nm)
  nm[nm == ""] <- paste0("column_", which(nm == ""))
  nm <- make.unique(nm, sep = "_")
  names(df) <- nm
  
  for (j in seq_along(df)) {
    if (is.list(df[[j]])) {
      df[[j]] <- vapply(df[[j]], function(z) paste(as.character(z), collapse = "; "), character(1))
    }
  }
  df
}

write_xlsx_list <- function(tables, file) {
  wb <- openxlsx::createWorkbook()
  header_style <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#D9EAF7", border = "Bottom")
  note_style <- openxlsx::createStyle(wrapText = TRUE, valign = "top")
  
  used_sheets <- character(0)
  for (nm in names(tables)) {
    df <- sanitize_df(tables[[nm]])
    if (is.null(df) || nrow(df) == 0 || ncol(df) == 0) next
    
    sheet <- substr(gsub("[\\[\\]\\*\\?/\\:]", "_", nm), 1, 31)
    if (sheet == "" || is.na(sheet)) sheet <- "Sheet"
    sheet <- make.unique(c(used_sheets, sheet), sep = "_") |> tail(1)
    used_sheets <- c(used_sheets, sheet)
    
    openxlsx::addWorksheet(wb, sheet)
    # Use writeData instead of writeDataTable to avoid hard failures caused by
    # duplicate/edge-case column names in manuscript diagnostic tables.
    openxlsx::writeData(wb, sheet, df, colNames = TRUE)
    openxlsx::addStyle(wb, sheet, header_style, rows = 1, cols = seq_len(ncol(df)), gridExpand = TRUE)
    openxlsx::addStyle(wb, sheet, note_style, rows = 1:(nrow(df) + 1), cols = seq_len(ncol(df)), gridExpand = TRUE)
    openxlsx::freezePane(wb, sheet, firstRow = TRUE)
    openxlsx::setColWidths(wb, sheet, cols = 1:ncol(df), widths = "auto")
  }
  openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
}

write_manuscript_style_table <- function(table_df, notes_df, file) {
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "Consolidated_Table")
  openxlsx::addWorksheet(wb, "Notes")
  
  title_style <- openxlsx::createStyle(textDecoration = "bold", fontSize = 12, halign = "center")
  header_style <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#EDEBE1", halign = "center", valign = "center", border = "TopBottomLeftRight")
  section_style <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#EDEBE1", halign = "center", border = "TopBottomLeftRight")
  body_style <- openxlsx::createStyle(halign = "center", valign = "center", border = "TopBottomLeftRight")
  var_style <- openxlsx::createStyle(halign = "left", valign = "center", border = "TopBottomLeftRight")
  
  title <- "Supplementary Table 1. Diagnostic and demographic composition of the analytic cohort by recruitment context"
  openxlsx::writeData(wb, "Consolidated_Table", title, startRow = 1, startCol = 1)
  openxlsx::mergeCells(wb, "Consolidated_Table", cols = 1:6, rows = 1)
  openxlsx::addStyle(wb, "Consolidated_Table", title_style, rows = 1, cols = 1:6, gridExpand = TRUE)
  
  df <- table_df %>% select(Stratum, Section, Variable, CN, AD, `p-value CN vs AD`)
  openxlsx::writeData(wb, "Consolidated_Table", df, startRow = 3, startCol = 1, colNames = TRUE)
  openxlsx::addStyle(wb, "Consolidated_Table", header_style, rows = 3, cols = 1:6, gridExpand = TRUE)
  
  data_start <- 4
  for (i in seq_len(nrow(df))) {
    rr <- data_start + i - 1
    openxlsx::addStyle(wb, "Consolidated_Table", body_style, rows = rr, cols = 1:6, gridExpand = TRUE)
    openxlsx::addStyle(wb, "Consolidated_Table", var_style, rows = rr, cols = 1:3, gridExpand = TRUE)
    if (df$Variable[i] == "") {
      openxlsx::mergeCells(wb, "Consolidated_Table", cols = 1:6, rows = rr)
      openxlsx::addStyle(wb, "Consolidated_Table", section_style, rows = rr, cols = 1:6, gridExpand = TRUE)
    }
  }
  openxlsx::setColWidths(wb, "Consolidated_Table", cols = 1:6, widths = c(24, 24, 34, 22, 22, 18))
  openxlsx::freezePane(wb, "Consolidated_Table", firstActiveRow = 4)
  
  openxlsx::writeData(wb, "Notes", sanitize_df(notes_df), colNames = TRUE)
  openxlsx::setColWidths(wb, "Notes", cols = 1:ncol(notes_df), widths = "auto")
  
  openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
}

write_manuscript_supplementary_package <- function(overall_table, country_table, site_table,
                                             composition_tests, notes_df, file) {
  wb <- openxlsx::createWorkbook()
  
  title_style <- openxlsx::createStyle(
    textDecoration = "bold", fontSize = 12,
    halign = "center", valign = "center"
  )
  header_style <- openxlsx::createStyle(
    textDecoration = "bold", fgFill = "#EDEBE1",
    halign = "center", valign = "center",
    border = "TopBottomLeftRight"
  )
  section_style <- openxlsx::createStyle(
    textDecoration = "bold", fgFill = "#EDEBE1",
    halign = "center", valign = "center",
    border = "TopBottomLeftRight"
  )
  body_center <- openxlsx::createStyle(
    halign = "center", valign = "center",
    border = "TopBottomLeftRight"
  )
  body_left <- openxlsx::createStyle(
    halign = "left", valign = "center",
    border = "TopBottomLeftRight"
  )
  note_style <- openxlsx::createStyle(wrapText = TRUE, valign = "top")
  
  # ---------------------------------------------------------------------------
  # Sheet 1: Supplementary Table 1-style overall analytic cohort table.
  # ---------------------------------------------------------------------------
  openxlsx::addWorksheet(wb, "Overall_analytic_cohort")
  
  t1 <- overall_table %>%
    mutate(
      Variables = ifelse(Variable == "", Section, Variable),
      CN = CN,
      `Clinical AD` = AD,
      `p-value CN vs clinical AD` = `p-value CN vs AD`,
      is_section = Variable == ""
    ) %>%
    select(Variables, CN, `Clinical AD`, `p-value CN vs clinical AD`, is_section)
  
  openxlsx::writeData(
    wb, "Overall_analytic_cohort",
    "Supplementary Table 1. Diagnostic and demographic composition of the analytic cohort by recruitment context",
    startRow = 1, startCol = 1
  )
  openxlsx::mergeCells(wb, "Overall_analytic_cohort", cols = 1:4, rows = 1)
  openxlsx::addStyle(wb, "Overall_analytic_cohort", title_style, rows = 1, cols = 1:4, gridExpand = TRUE)
  
  display_t1 <- t1 %>% select(-is_section)
  openxlsx::writeData(wb, "Overall_analytic_cohort", display_t1, startRow = 3, startCol = 1, colNames = TRUE)
  openxlsx::addStyle(wb, "Overall_analytic_cohort", header_style, rows = 3, cols = 1:4, gridExpand = TRUE)
  
  for (i in seq_len(nrow(t1))) {
    rr <- 3 + i
    if (isTRUE(t1$is_section[i])) {
      openxlsx::mergeCells(wb, "Overall_analytic_cohort", cols = 1:4, rows = rr)
      openxlsx::addStyle(wb, "Overall_analytic_cohort", section_style, rows = rr, cols = 1:4, gridExpand = TRUE)
    } else {
      openxlsx::addStyle(wb, "Overall_analytic_cohort", body_left, rows = rr, cols = 1, gridExpand = TRUE)
      openxlsx::addStyle(wb, "Overall_analytic_cohort", body_center, rows = rr, cols = 2:4, gridExpand = TRUE)
    }
  }
  
  foot_row <- nrow(t1) + 6
  footnote <- paste(
    "Values are shown as mean ± SD unless otherwise indicated.",
    "Clinical AD refers to clinically diagnosed Alzheimer's disease assigned using harmonized clinical criteria and supported by neuropsychological assessment.",
    "This table summarizes the analytic cohort with complete covariate information used in the primary differential-abundance models.",
    "APOE ε4 carrier frequency is reported among participants with available genotype data.",
    "P values correspond to CN versus clinically diagnosed AD comparisons."
  )
  openxlsx::writeData(wb, "Overall_analytic_cohort", footnote, startRow = foot_row, startCol = 1)
  openxlsx::mergeCells(wb, "Overall_analytic_cohort", cols = 1:4, rows = foot_row)
  openxlsx::addStyle(wb, "Overall_analytic_cohort", note_style, rows = foot_row, cols = 1:4, gridExpand = TRUE)
  openxlsx::setColWidths(wb, "Overall_analytic_cohort", cols = 1:4, widths = c(34, 22, 22, 24))
  
  # ---------------------------------------------------------------------------
  # Sheet 2-3: compact country/site composition tables.
  # ---------------------------------------------------------------------------
  write_context_sheet <- function(sheet, tbl, title) {
    if (is.null(tbl) || nrow(tbl) == 0) return(invisible(NULL))
    openxlsx::addWorksheet(wb, sheet)
    
    tbl2 <- tbl %>%
      select(-any_of("Context")) %>%
      rename(
        `Recruitment context` = Level,
        `Total N` = N,
        `CN, n` = CN,
        `Clinical AD, n` = AD,
        `% clinical AD` = `% AD`
      )
    
    openxlsx::writeData(wb, sheet, title, startRow = 1, startCol = 1)
    openxlsx::mergeCells(wb, sheet, cols = 1:ncol(tbl2), rows = 1)
    openxlsx::addStyle(wb, sheet, title_style, rows = 1, cols = 1:ncol(tbl2), gridExpand = TRUE)
    openxlsx::writeData(wb, sheet, tbl2, startRow = 3, startCol = 1, colNames = TRUE)
    openxlsx::addStyle(wb, sheet, header_style, rows = 3, cols = 1:ncol(tbl2), gridExpand = TRUE)
    if (nrow(tbl2) > 0) {
      openxlsx::addStyle(wb, sheet, body_center, rows = 4:(nrow(tbl2) + 3), cols = 1:ncol(tbl2), gridExpand = TRUE)
      openxlsx::addStyle(wb, sheet, body_left, rows = 4:(nrow(tbl2) + 3), cols = 1, gridExpand = TRUE)
    }
    openxlsx::setColWidths(wb, sheet, cols = 1:ncol(tbl2), widths = "auto")
    openxlsx::freezePane(wb, sheet, firstActiveRow = 4)
  }
  
  write_context_sheet(
    "By_country",
    country_table,
    "Supplementary Table 1a. Diagnostic and demographic composition by country"
  )
  write_context_sheet(
    "By_recruitment_site",
    site_table,
    "Supplementary Table 1b. Diagnostic and demographic composition by recruitment site"
  )
  
  # ---------------------------------------------------------------------------
  # Sheet 4-5: tests and notes.
  # ---------------------------------------------------------------------------
  openxlsx::addWorksheet(wb, "Composition_Tests")
  tests <- composition_tests %>%
    rename(
      Level = level,
      Variable = variable,
      Test = test,
      `p value` = p_value,
      `Formatted p value` = p_value_formatted
    )
  openxlsx::writeData(wb, "Composition_Tests", tests, colNames = TRUE)
  openxlsx::addStyle(wb, "Composition_Tests", header_style, rows = 1, cols = 1:ncol(tests), gridExpand = TRUE)
  openxlsx::setColWidths(wb, "Composition_Tests", cols = 1:ncol(tests), widths = "auto")
  
  openxlsx::addWorksheet(wb, "Notes")
  openxlsx::writeData(wb, "Notes", sanitize_df(notes_df), colNames = TRUE)
  openxlsx::addStyle(wb, "Notes", header_style, rows = 1, cols = 1:ncol(notes_df), gridExpand = TRUE)
  openxlsx::addStyle(wb, "Notes", note_style, rows = 1:(nrow(notes_df) + 1), cols = 1:ncol(notes_df), gridExpand = TRUE)
  openxlsx::setColWidths(wb, "Notes", cols = 1:ncol(notes_df), widths = c(24, 120))
  
  openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
}

make_protein_name_review <- function(symbol, target_full, target, apt) {
  symbol <- clean_text_na(symbol)
  target_full <- clean_text_na(target_full)
  target <- clean_text_na(target)
  apt <- clean_text_na(apt)
  dplyr::coalesce(symbol, target_full, target, apt)
}

classify_dep_review <- function(df, fdr = MAIN_FDR) {
  df %>%
    mutate(
      type = case_when(
        adj.P.Val < fdr & logFC > 0 ~ "Up",
        adj.P.Val < fdr & logFC < 0 ~ "Down",
        TRUE ~ "NS"
      ),
      significant_fdr005 = adj.P.Val < 0.05,
      significant_fdr001 = adj.P.Val < 0.01
    )
}

collapse_to_gene_review <- function(dep_tbl) {
  dep_tbl %>%
    mutate(
      EntrezGeneSymbol = as.character(EntrezGeneSymbol),
      AptName = as.character(AptName)
    ) %>%
    filter(!is.na(EntrezGeneSymbol), EntrezGeneSymbol != "") %>%
    arrange(adj.P.Val, desc(abs(logFC)), AptName) %>%
    group_by(EntrezGeneSymbol) %>%
    slice(1) %>%
    ungroup()
}

run_limma_dep_review <- function(dat, protein_cols, annot_tbl, formula_str,
                                 coef_name = "SampleGroupAD", model_name = "manuscript_model") {
  dat <- dat %>% as_tibble()
  protein_cols <- intersect(protein_cols, names(dat))
  if (length(protein_cols) == 0) stop("No protein columns available for ", model_name)
  
  metadata <- dat %>% select(-all_of(protein_cols))
  model_vars <- all.vars(as.formula(formula_str))
  missing_model_vars <- setdiff(model_vars, names(metadata))
  if (length(missing_model_vars) > 0) {
    stop("Missing model variables in ", model_name, ": ", paste(missing_model_vars, collapse = ", "))
  }
  
  keep <- complete.cases(metadata[, model_vars, drop = FALSE])
  dat <- dat[keep, , drop = FALSE]
  metadata <- metadata[keep, , drop = FALSE]
  
  expr <- dat %>% select(all_of(protein_cols)) %>% as.matrix()
  expr <- t(safe_log2_matrix(expr))
  design <- model.matrix(as.formula(formula_str), data = metadata)
  
  if (!coef_name %in% colnames(design)) {
    stop("Coefficient not found in design matrix for ", model_name, ": ", coef_name)
  }
  
  fit <- limma::lmFit(expr, design)
  fit <- limma::eBayes(fit)
  
  tt <- limma::topTable(fit, coef = coef_name, adjust.method = "BH", number = Inf) %>%
    tibble::rownames_to_column(var = "feature_id_raw") %>%
    mutate(
      AptName = as.character(feature_id_raw),
      se = safe_se_from_limma(logFC, t)
    ) %>%
    left_join(annot_tbl, by = "AptName") %>%
    mutate(
      Protein_Name = make_protein_name_review(EntrezGeneSymbol, TargetFullName, Target, AptName),
      model = model_name
    ) %>%
    classify_dep_review(fdr = MAIN_FDR)
  
  list(
    dep_aptamer = tt,
    dep_gene = collapse_to_gene_review(tt),
    metadata = metadata,
    design = design
  )
}

compare_to_main_review <- function(main_gene, sensitivity_gene, sensitivity_label) {
  main_gene %>%
    select(
      EntrezGeneSymbol, AptName, Protein_Name, logFC, P.Value, adj.P.Val, type
    ) %>%
    rename(
      main_AptName = AptName,
      main_Protein_Name = Protein_Name,
      main_logFC = logFC,
      main_P.Value = P.Value,
      main_adj.P.Val = adj.P.Val,
      main_type = type
    ) %>%
    inner_join(
      sensitivity_gene %>%
        select(EntrezGeneSymbol, AptName, Protein_Name, logFC, P.Value, adj.P.Val, type) %>%
        rename(
          sensitivity_AptName = AptName,
          sensitivity_Protein_Name = Protein_Name,
          sensitivity_logFC = logFC,
          sensitivity_P.Value = P.Value,
          sensitivity_adj.P.Val = adj.P.Val,
          sensitivity_type = type
        ),
      by = "EntrezGeneSymbol"
    ) %>%
    mutate(
      sensitivity = sensitivity_label,
      Protein_Name = coalesce(main_Protein_Name, sensitivity_Protein_Name),
      same_direction = sign(main_logFC) == sign(sensitivity_logFC),
      main_sig_fdr005 = main_adj.P.Val < 0.05,
      sensitivity_sig_fdr005 = sensitivity_adj.P.Val < 0.05,
      preserved_fdr005 = main_sig_fdr005 & sensitivity_sig_fdr005 & same_direction,
      main_sig_fdr001 = main_adj.P.Val < 0.01,
      sensitivity_sig_fdr001 = sensitivity_adj.P.Val < 0.01,
      preserved_fdr001 = main_sig_fdr001 & sensitivity_sig_fdr001 & same_direction,
      delta_logFC = sensitivity_logFC - main_logFC,
      attenuation_ratio = abs(sensitivity_logFC) / pmax(abs(main_logFC), 1e-9)
    )
}

summarize_comparison_review <- function(compare_tbl, sensitivity_label) {
  sig_main <- compare_tbl %>% filter(main_sig_fdr005)
  tibble(
    sensitivity = sensitivity_label,
    n_common_gene_symbols = nrow(compare_tbl),
    n_main_fdr005 = sum(compare_tbl$main_sig_fdr005, na.rm = TRUE),
    n_sensitivity_fdr005 = sum(compare_tbl$sensitivity_sig_fdr005, na.rm = TRUE),
    n_preserved_fdr005_same_direction = sum(compare_tbl$preserved_fdr005, na.rm = TRUE),
    proportion_main_fdr005_preserved = ifelse(nrow(sig_main) > 0, mean(sig_main$preserved_fdr005, na.rm = TRUE), NA_real_),
    logFC_correlation_all = suppressWarnings(cor(compare_tbl$main_logFC, compare_tbl$sensitivity_logFC, use = "complete.obs")),
    direction_consistency_all = mean(compare_tbl$same_direction, na.rm = TRUE),
    direction_consistency_main_fdr005 = ifelse(nrow(sig_main) > 0, mean(sig_main$same_direction, na.rm = TRUE), NA_real_),
    median_attenuation_ratio_main_fdr005 = ifelse(nrow(sig_main) > 0, median(sig_main$attenuation_ratio, na.rm = TRUE), NA_real_)
  )
}

smd_numeric <- function(x, g) {
  x <- suppressWarnings(as.numeric(x))
  g <- as.character(g)
  if (!all(c("CN", "AD") %in% unique(g))) return(NA_real_)
  x0 <- x[g == "CN"]
  x1 <- x[g == "AD"]
  s_pool <- sqrt((var(x0, na.rm = TRUE) + var(x1, na.rm = TRUE)) / 2)
  if (!is.finite(s_pool) || s_pool == 0) return(NA_real_)
  (mean(x1, na.rm = TRUE) - mean(x0, na.rm = TRUE)) / s_pool
}

smd_binary <- function(x, g, level) {
  g <- as.character(g)
  x <- as.character(x)
  if (!all(c("CN", "AD") %in% unique(g))) return(NA_real_)
  p0 <- mean(x[g == "CN"] == level, na.rm = TRUE)
  p1 <- mean(x[g == "AD"] == level, na.rm = TRUE)
  p <- (p0 + p1) / 2
  denom <- sqrt(p * (1 - p))
  if (!is.finite(denom) || denom == 0) return(NA_real_)
  (p1 - p0) / denom
}

make_balance_table <- function(dat, label, site_var = NA_character_) {
  vars_numeric <- c("Age", "Education")
  vars_numeric <- vars_numeric[vars_numeric %in% names(dat)]
  vars_cat <- c("Sex", "Country")
  if (!is.na(site_var) && site_var %in% names(dat)) vars_cat <- c(vars_cat, site_var)
  vars_cat <- vars_cat[vars_cat %in% names(dat)]
  
  num_tbl <- purrr::map_dfr(vars_numeric, function(v) {
    tibble(
      dataset = label,
      variable = v,
      level = NA_character_,
      type = "numeric",
      smd = smd_numeric(dat[[v]], dat$SampleGroup),
      CN = mean_sd(dat[[v]][dat$SampleGroup == "CN"]),
      AD = mean_sd(dat[[v]][dat$SampleGroup == "AD"])
    )
  })
  
  cat_tbl <- purrr::map_dfr(vars_cat, function(v) {
    lvls <- sort(unique(as.character(dat[[v]])))
    purrr::map_dfr(lvls, function(lv) {
      tibble(
        dataset = label,
        variable = v,
        level = lv,
        type = "categorical_level",
        smd = smd_binary(dat[[v]], dat$SampleGroup, lv),
        CN = paste0(sum(dat$SampleGroup == "CN" & as.character(dat[[v]]) == lv, na.rm = TRUE), " / ", sum(dat$SampleGroup == "CN", na.rm = TRUE)),
        AD = paste0(sum(dat$SampleGroup == "AD" & as.character(dat[[v]]) == lv, na.rm = TRUE), " / ", sum(dat$SampleGroup == "AD", na.rm = TRUE))
      )
    })
  })
  
  bind_rows(num_tbl, cat_tbl)
}

###############################################################################
# 02. Required objects and core analytic data
###############################################################################

required_objects <- c("dep_df", "dep_protein_cols", "annot_tbl", "DEP_gene")
missing_objects <- required_objects[!vapply(required_objects, exists, logical(1))]
if (length(missing_objects) > 0) {
  stop("The workspace is missing required objects: ", paste(missing_objects, collapse = ", "))
}

review_df <- dep_df %>%
  mutate(
    SampleGroup = factor(as.character(SampleGroup), levels = MAIN_GROUPS),
    Sex = factor(Sex),
    Country = factor(Country)
  ) %>%
  filter(SampleGroup %in% MAIN_GROUPS)

# Site detection: add more candidates if your actual site variable has a different name.
site_candidates <- c("Site", "site", "Center", "center", "Cohort", "cohort", "RecruitmentSite", "recruitment_site", "site_id", "Site_ID")
site_var <- site_candidates[site_candidates %in% names(review_df)][1]
if (length(site_var) == 0 || is.na(site_var)) site_var <- NA_character_

# APOE4 detection. If APOE4_carrier does not exist, try common alternatives.
apoe_candidates <- c("APOE4_carrier", "apoe4_carrier", "APOE_e4_carrier", "APOE4", "apoe4")
apoe_var <- apoe_candidates[apoe_candidates %in% names(review_df)][1]
if (length(apoe_var) == 0 || is.na(apoe_var)) apoe_var <- NA_character_

###############################################################################
# 03. Consolidated manuscript-style supplementary table
###############################################################################


# Flexible variable detection for a manuscript-ready Supplementary Table 1.
# This keeps the table stable even if the workspace uses slightly different
# names for neuropsychological or AT(N) variables. Detection is performed in
# three steps: exact match, case-insensitive match, and normalized-name match.
normalize_var_name <- function(x) {
  x <- tolower(as.character(x))
  x <- gsub("β", "b", x)
  x <- gsub("[^a-z0-9]+", "", x)
  x
}

detect_first_col <- function(candidates, data = review_df, regex = NULL) {
  candidates <- candidates[!is.na(candidates)]
  nm <- names(data)

  # 1) exact match
  hit <- candidates[candidates %in% nm][1]
  if (length(hit) > 0 && !is.na(hit)) return(hit)

  # 2) case-insensitive exact match
  nm_lower <- tolower(nm)
  cand_lower <- tolower(candidates)
  idx <- match(cand_lower, nm_lower)
  idx <- idx[!is.na(idx)]
  if (length(idx) > 0) return(nm[idx[1]])

  # 3) normalized match, robust to punctuation such as CDR.SB / CDR-SB / cdr_sb
  nm_norm <- normalize_var_name(nm)
  cand_norm <- normalize_var_name(candidates)
  idx <- match(cand_norm, nm_norm)
  idx <- idx[!is.na(idx)]
  if (length(idx) > 0) return(nm[idx[1]])

  # 4) optional regex fallback
  if (!is.null(regex)) {
    idx <- grep(regex, nm, ignore.case = TRUE, perl = TRUE)
    if (length(idx) > 0) return(nm[idx[1]])
  }

  NA_character_
}

variable_dictionary <- tibble::tribble(
  ~section, ~label, ~var, ~type,
  "Demographic profile", "N", NA_character_, "n",
  "Demographic profile", "Age", detect_first_col(c("Age", "age", "Edad", "edad"), regex = "^age$|^edad$"), "continuous",
  "Demographic profile", "Female, n (%)", detect_first_col(c("Sex", "sex", "Gender", "gender", "sexo"), regex = "^sex$|^gender$|^sexo$"), "female",
  "Demographic profile", "Years of Education", detect_first_col(c("Education", "education", "YearsEducation", "years_education", "educacion", "edu", "Education_years"), regex = "educ|school|years.*educ"), "continuous",
  
  "Genetic profile", "APOE ε4 carriers, n/N (%)", apoe_var, "apoe4",
  
  "Neuropsychological tests", "CDR-SB", detect_first_col(c("CDR.SB", "CDR_SB", "CDR-SB", "cdr_sb", "cdrsb", "CDRSB", "cdr_sum_boxes", "CDRSumOfBoxes", "CDR_Sum_of_Boxes"), regex = "cdr.*(sb|sum|box)|cdrsb"), "continuous",
  "Neuropsychological tests", "MMSE Score", detect_first_col(c("MMSE", "MMSE.Score", "MMSE_score", "mmse_total", "mmse", "MMSE_total"), regex = "mmse"), "continuous",
  "Neuropsychological tests", "PFAQ Score", detect_first_col(c("PFAQ", "PFAQ.Score", "PFAQ_score", "pfaq_total", "pfaq", "FAQ", "FAQ_score", "UDSFAQ", "UDS_FAQ", "functional_activities_questionnaire"), regex = "pfaq|uds.*faq|^faq$|functional.*activ"), "continuous",
  "Neuropsychological tests", "Mini-SEA Score", detect_first_col(c("Mini.Sea", "MiniSea", "Mini_SEA", "MiniSEA", "mini_sea", "mini_sea_total", "MiniSEA_score", "Mini.Sea.Score", "SEA", "social_cognition"), regex = "mini.*sea|sea.*score|social.*cog"), "continuous",
  "Neuropsychological tests", "T-ADLQ score", detect_first_col(c("T.ADLQ", "T_ADLQ", "T-ADLQ", "TADLQ", "t_adlq", "ADLQ", "adlq_total", "TADLQ_score", "Technology_ADLQ"), regex = "t.*adlq|tadlq|adlq|technology.*adl"), "continuous",
  "Neuropsychological tests", "NPI-Q Score", detect_first_col(c("NPI.Q", "NPI_Q", "NPI-Q", "NPIQ", "npi_q", "npi_total", "NPIQ_score", "NPI_total"), regex = "npi"), "continuous",
  
  "Plasma AT(N) biomarkers", "p-tau217", detect_first_col(c("p.tau217", "p-tau217", "p_tau217", "ptau217", "ptau_217", "pTau217"), regex = "p.*tau.*217|ptau217"), "continuous",
  "Plasma AT(N) biomarkers", "p-tau181", detect_first_col(c("p.tau181", "p-tau181", "p_tau181", "ptau181", "ptau_181", "pTau181"), regex = "p.*tau.*181|ptau181"), "continuous",
  "Plasma AT(N) biomarkers", "NfL", detect_first_col(c("NfL", "NFL", "nfl", "Neurofilament_light", "neurofilament_light"), regex = "^nfl$|neurofilament"), "continuous",
  "Plasma AT(N) biomarkers", "Aβ42/40 ratio", detect_first_col(c("ratio.AB42.40", "ratio AB42/40", "ratio_AB42_40", "AB42_40", "abeta42_40", "Aβ42/40", "Abeta42_40", "AB42_AB40_ratio"), regex = "(ab|a).*42.*40|ratio.*42.*40"), "continuous"
)

# Export a variable-detection map so missing fields can be audited explicitly.
variable_detection_map <- variable_dictionary %>%
  transmute(
    section,
    label,
    detected_column = ifelse(is.na(var), "Not detected", var),
    type
  )

make_one_stratum_table <- function(dat, stratum_label) {
  dat <- dat %>% filter(SampleGroup %in% c("CN", "AD"))
  n_cn <- sum(dat$SampleGroup == "CN", na.rm = TRUE)
  n_ad <- sum(dat$SampleGroup == "AD", na.rm = TRUE)
  
  rows <- purrr::pmap_dfr(variable_dictionary, function(section, label, var, type) {
    if (type == "n") {
      return(tibble(
        Stratum = stratum_label,
        Section = section,
        Variable = label,
        CN = as.character(n_cn),
        AD = as.character(n_ad),
        `p-value CN vs AD` = "-",
        Test = "-"
      ))
    }
    
    if (is.na(var) || !var %in% names(dat)) {
      return(tibble(
        Stratum = stratum_label,
        Section = section,
        Variable = label,
        CN = "Not available",
        AD = "Not available",
        `p-value CN vs AD` = "-",
        Test = "Variable not detected"
      ))
    }
    
    if (type == "continuous") {
      p <- safe_wilcox_p(dat[[var]], dat$SampleGroup)
      return(tibble(
        Stratum = stratum_label,
        Section = section,
        Variable = label,
        CN = mean_sd(dat[[var]][dat$SampleGroup == "CN"]),
        AD = mean_sd(dat[[var]][dat$SampleGroup == "AD"]),
        `p-value CN vs AD` = format_p(p),
        Test = "Wilcoxon rank-sum test"
      ))
    }
    
    if (type == "female") {
      fem <- is_female_vec(dat[[var]])
      p <- safe_binary_p(fem, dat$SampleGroup)
      return(tibble(
        Stratum = stratum_label,
        Section = section,
        Variable = label,
        CN = n_percent(sum(fem[dat$SampleGroup == "CN"], na.rm = TRUE), sum(dat$SampleGroup == "CN" & !is.na(fem))),
        AD = n_percent(sum(fem[dat$SampleGroup == "AD"], na.rm = TRUE), sum(dat$SampleGroup == "AD" & !is.na(fem))),
        `p-value CN vs AD` = format_p(p),
        Test = "Chi-square or Fisher exact test"
      ))
    }
    
    if (type == "apoe4") {
      car <- is_apoe4_carrier_vec(dat[[var]])
      p <- safe_binary_p(car, dat$SampleGroup)
      return(tibble(
        Stratum = stratum_label,
        Section = section,
        Variable = label,
        CN = n_over_n_percent(sum(car[dat$SampleGroup == "CN"], na.rm = TRUE), sum(dat$SampleGroup == "CN" & !is.na(car))),
        AD = n_over_n_percent(sum(car[dat$SampleGroup == "AD"], na.rm = TRUE), sum(dat$SampleGroup == "AD" & !is.na(car))),
        `p-value CN vs AD` = format_p(p),
        Test = "Chi-square or Fisher exact test"
      ))
    }
    
    tibble(Stratum = stratum_label, Section = section, Variable = label, CN = "-", AD = "-", `p-value CN vs AD` = "-", Test = "Unknown")
  })
  
  # Add visual section separator rows without group_modify(), because
  # dplyr::group_modify() does not allow returning the grouping variable.
  section_order <- unique(rows$Section)
  purrr::map_dfr(section_order, function(sec) {
    bind_rows(
      tibble(
        Stratum = stratum_label,
        Section = sec,
        Variable = "",
        CN = "",
        AD = "",
        `p-value CN vs AD` = "",
        Test = "Section header"
      ),
      rows %>% filter(Section == sec)
    )
  })
}

overall_table <- make_one_stratum_table(review_df, "Overall")

country_tables <- purrr::map_dfr(levels(droplevels(review_df$Country)), function(cc) {
  make_one_stratum_table(review_df %>% filter(Country == cc), paste0("Country: ", cc))
})

site_tables <- NULL
if (!is.na(site_var) && site_var %in% names(review_df)) {
  site_values <- sort(unique(as.character(review_df[[site_var]])))
  site_tables <- purrr::map_dfr(site_values, function(ss) {
    make_one_stratum_table(review_df %>% filter(as.character(.data[[site_var]]) == ss), paste0("Site: ", ss))
  })
}

composition_consolidated <- bind_rows(overall_table, country_tables, site_tables)

composition_tests <- bind_rows(
  tibble(
    level = "Country",
    variable = "Diagnostic composition across countries",
    test = "Chi-square/Fisher simulation as needed",
    p_value = safe_p_chisq_or_fisher(table(review_df$Country, review_df$SampleGroup)),
    p_value_formatted = format_p(p_value)
  ),
  if (!is.na(site_var) && site_var %in% names(review_df)) {
    tibble(
      level = "Recruitment site",
      variable = "Diagnostic composition across recruitment sites",
      test = "Chi-square/Fisher simulation as needed",
      p_value = safe_p_chisq_or_fisher(table(review_df[[site_var]], review_df$SampleGroup)),
      p_value_formatted = format_p(p_value)
    )
  }
)

composition_notes <- tibble(
  item = c(
    "Analytic sample",
    "Diagnostic definition",
    "Statistical tests",
    "Country handling",
    "Site handling",
    "Interpretation"
  ),
  note = c(
    "This table summarizes the analytic cohort with complete covariate information used in the primary differential-abundance models; counts may therefore differ from the full enrolled cohort.",
    "Clinical AD refers to clinically diagnosed Alzheimer's disease assigned using harmonized clinical criteria and supported by neuropsychological assessment.",
    "Continuous variables are summarized as mean ± SD and compared using Wilcoxon rank-sum tests within each stratum. Categorical variables are summarized as n (%) or n/N (%) and compared using chi-square or Fisher exact tests depending on expected counts.",
    "Country is treated as a categorical recruitment-context variable, not as an ordinal/numeric variable.",
    ifelse(is.na(site_var), "No site/center variable was detected in dep_df.", paste0("Recruitment site variable detected as: ", site_var)),
    "This table documents recruitment-context and covariate composition and should be interpreted as descriptive cohort characterization, not as external validation."
  )
)

# A clean Table 1-style version for direct manuscript/supplementary use.
# This removes the Stratum/Test columns and keeps only the standard table columns.
overall_table_manuscript <- overall_table %>%
  filter(Stratum == "Overall") %>%
  select(Section, Variable, CN, AD, `p-value CN vs AD`)

make_compact_context_table <- function(dat, context_var, context_label) {
  if (is.na(context_var) || !context_var %in% names(dat)) return(NULL)
  
  contexts <- sort(unique(as.character(dat[[context_var]])))
  purrr::map_dfr(contexts, function(ctx) {
    dd <- dat %>% filter(as.character(.data[[context_var]]) == ctx)
    n_total <- nrow(dd)
    n_cn <- sum(dd$SampleGroup == "CN", na.rm = TRUE)
    n_ad <- sum(dd$SampleGroup == "AD", na.rm = TRUE)
    fem <- if ("Sex" %in% names(dd)) is_female_vec(dd$Sex) else rep(NA, nrow(dd))
    apo <- if (!is.na(apoe_var) && apoe_var %in% names(dd)) is_apoe4_carrier_vec(dd[[apoe_var]]) else rep(NA, nrow(dd))
    
    tibble(
      Context = context_label,
      Level = ctx,
      N = n_total,
      CN = n_cn,
      AD = n_ad,
      `% AD` = sprintf("%.1f", 100 * n_ad / pmax(n_total, 1)),
      `Age CN` = mean_sd(dd$Age[dd$SampleGroup == "CN"]),
      `Age AD` = mean_sd(dd$Age[dd$SampleGroup == "AD"]),
      `Education CN` = mean_sd(dd$Education[dd$SampleGroup == "CN"]),
      `Education AD` = mean_sd(dd$Education[dd$SampleGroup == "AD"]),
      `Female CN, n (%)` = n_percent(sum(fem[dd$SampleGroup == "CN"], na.rm = TRUE), sum(dd$SampleGroup == "CN" & !is.na(fem))),
      `Female AD, n (%)` = n_percent(sum(fem[dd$SampleGroup == "AD"], na.rm = TRUE), sum(dd$SampleGroup == "AD" & !is.na(fem))),
      `APOE e4 CN, n/N (%)` = n_over_n_percent(sum(apo[dd$SampleGroup == "CN"], na.rm = TRUE), sum(dd$SampleGroup == "CN" & !is.na(apo))),
      `APOE e4 AD, n/N (%)` = n_over_n_percent(sum(apo[dd$SampleGroup == "AD"], na.rm = TRUE), sum(dd$SampleGroup == "AD" & !is.na(apo)))
    )
  })
}

country_compact_table <- make_compact_context_table(review_df, "Country", "Country")
site_compact_table <- make_compact_context_table(review_df, site_var, "Recruitment site")

safe_write_csv(variable_detection_map, file.path(out_dir, "supplementary_table_1_variable_detection_map.csv"))
safe_write_csv(composition_consolidated, file.path(out_dir, "supplementary_table_country_site_composition_consolidated.csv"))
safe_write_csv(composition_tests, file.path(out_dir, "country_site_diagnostic_composition_tests.csv"))

write_manuscript_style_table(
  table_df = composition_consolidated,
  notes_df = bind_rows(
    composition_notes,
    composition_tests %>% transmute(item = level, note = paste0(variable, "; ", test, "; p = ", p_value_formatted))
  ),
  file = file.path(out_dir, "Supplementary_Table_Country_Site_Composition_CONSOLIDATED.xlsx")
)

# Additional compact workbook: closer to journal Table 1 and easier to inspect.
write_xlsx_list(
  list(
    Table1_style_overall = overall_table_manuscript,
    Country_composition = country_compact_table,
    Site_composition = site_compact_table,
    Composition_tests = composition_tests,
    Notes = composition_notes,
    Variable_detection = variable_detection_map
  ),
  file.path(out_dir, "Supplementary_Table_Cohort_Composition_MANUSCRIPT_READY.xlsx")
)

# Manuscript-ready Supplementary Table 1 workbook.
# This creates:
#   1) an overall analytic cohort table,
#   2) a compact country composition table,
#   3) a compact recruitment-site composition table,
#   4) diagnostic composition tests and notes.
write_manuscript_supplementary_package(
  overall_table = overall_table_manuscript,
  country_table = country_compact_table,
  site_table = site_compact_table,
  composition_tests = composition_tests,
  notes_df = bind_rows(
    composition_notes,
    composition_tests %>%
      transmute(item = level, note = paste0(variable, "; ", test, "; p = ", p_value_formatted))
  ),
  file = file.path(out_dir, "Supplementary_Table_1_Diagnostic_demographic_composition_by_recruitment_context.xlsx")
)

# Compatibility outputs similar to the previous script, but kept secondary.
country_group_counts <- review_df %>%
  count(Country, SampleGroup, name = "n") %>%
  tidyr::pivot_wider(names_from = SampleGroup, values_from = n, values_fill = 0) %>%
  mutate(total = CN + AD, percent_AD = 100 * AD / pmax(total, 1)) %>%
  arrange(desc(total))

site_group_counts <- NULL
if (!is.na(site_var) && site_var %in% names(review_df)) {
  site_group_counts <- review_df %>%
    count(.data[[site_var]], SampleGroup, name = "n") %>%
    rename(site = !!rlang::sym(site_var)) %>%
    tidyr::pivot_wider(names_from = SampleGroup, values_from = n, values_fill = 0) %>%
    mutate(total = CN + AD, percent_AD = 100 * AD / pmax(total, 1)) %>%
    arrange(desc(total))
}

write_xlsx_list(
  list(
    Consolidated_Table = composition_consolidated,
    Notes = composition_notes,
    Composition_tests = composition_tests,
    Country_counts = country_group_counts,
    Site_counts = site_group_counts
  ),
  file.path(out_dir, "Supplementary_Table_Country_Site_Composition_SUPPORTING.xlsx")
)

###############################################################################
# 04. Sensitivity analysis 1: MatchIt propensity-score matched DEP
###############################################################################

matched_df <- tibble()
matched_summary <- tibble()
balance_tbl <- tibble()
matched_dep_gene <- tibble()
matched_compare <- tibble()
matched_compare_summary <- tibble()

# IMPORTANT:
# This matching block intentionally follows the standalone Matching.R logic:
# nearest-neighbor propensity-score matching without replacement using only
# Sex, Age and Education. Country/site are NOT used for matching, but Country is
# retained as an adjustment covariate in the downstream DEP model.
matching_covariates <- c("Sex", "Age", "Education")

match_df <- review_df %>%
  mutate(
    SampleGroup_bin = case_when(
      SampleGroup == "AD" ~ 1,
      SampleGroup == "CN" ~ 0,
      TRUE ~ NA_real_
    ),
    Sex = factor(Sex),
    Age = as.numeric(Age),
    Education = as.numeric(Education),
    Country = factor(Country)
  ) %>%
  filter(
    !is.na(SampleGroup_bin),
    !is.na(Sex),
    !is.na(Age),
    !is.na(Education)
  )

cat("\nMatchIt matching input sample:\n")
cat("Total before matching filter: ", nrow(review_df), "\n", sep = "")
cat("Total after matching filter: ", nrow(match_df), "\n", sep = "")
cat("CN before matching: ", sum(match_df$SampleGroup == "CN"), "\n", sep = "")
cat("AD before matching: ", sum(match_df$SampleGroup == "AD"), "\n", sep = "")

set.seed(1111)

matchit_out <- MatchIt::matchit(
  formula = SampleGroup_bin ~ Sex + Age + Education,
  data = match_df,
  method = "nearest",
  distance = "logit",
  replace = FALSE,
  caliper = 0.20
)

matchit_summary <- summary(matchit_out)

matched_df <- MatchIt::match.data(matchit_out) %>%
  as_tibble() %>%
  mutate(
    SampleGroup = factor(as.character(SampleGroup), levels = MAIN_GROUPS),
    Sex = factor(Sex),
    Country = factor(Country),
    Age = as.numeric(Age),
    Education = as.numeric(Education)
  )

# MatchIt creates "subclass" as the matched-pair identifier.
if ("subclass" %in% names(matched_df)) {
  matched_df <- matched_df %>%
    mutate(match_pair_id = as.character(subclass))
} else {
  matched_df <- matched_df %>%
    mutate(match_pair_id = NA_character_)
}

matched_summary <- tibble(
  sensitivity = "MatchIt_nearest_neighbor_matched_DEP",
  matching_covariates = paste(matching_covariates, collapse = ", "),
  matching_method = "nearest-neighbor propensity-score matching",
  distance = "logit",
  replace = FALSE,
  caliper = 0.20,
  country_in_matching = FALSE,
  site_in_matching = FALSE,
  country_in_DEP_model = TRUE,
  n_before = nrow(match_df),
  n_before_CN = sum(match_df$SampleGroup == "CN"),
  n_before_AD = sum(match_df$SampleGroup == "AD"),
  n_after = nrow(matched_df),
  n_after_CN = sum(matched_df$SampleGroup == "CN"),
  n_after_AD = sum(matched_df$SampleGroup == "AD"),
  n_pairs = length(unique(na.omit(matched_df$match_pair_id)))
)

balance_before <- make_balance_table(match_df, "before_matching", site_var = site_var)
balance_after  <- make_balance_table(matched_df, "after_matching", site_var = site_var)
balance_tbl <- bind_rows(balance_before, balance_after)

safe_write_csv(matched_df, file.path(out_dir, "Matched_Output_MatchIt.csv"))
safe_write_csv(matched_summary, file.path(out_dir, "matched_sample_summary.csv"))
safe_write_csv(balance_tbl, file.path(out_dir, "matched_covariate_balance_before_after.csv"))

# Export MatchIt/cobalt balance diagnostics for reviewer transparency.
capture.output(
  matchit_summary,
  file = file.path(out_dir, "MatchIt_summary.txt")
)

capture.output(
  cobalt::bal.tab(matchit_out, un = TRUE),
  file = file.path(out_dir, "MatchIt_cobalt_balance_table.txt")
)

# MatchIt base diagnostic plots can pause interactively in RStudio depending on
# graphics settings. To keep the manuscript pipeline fully non-interactive, we
# skip those base plots and retain the reproducible cobalt balance outputs:
#   1) MatchIt_summary.txt
#   2) MatchIt_cobalt_balance_table.txt
#   3) MatchIt_love_plot.pdf
# These are sufficient to document matching balance for reviewer inspection.

pdf(file.path(out_dir, "MatchIt_love_plot.pdf"), width = 8, height = 6)
try(
  print(
    cobalt::love.plot(
      matchit_out,
      stat = "mean.diffs",
      threshold = 0.1,
      abs = TRUE,
      var.order = "unadjusted",
      sample.names = c("Unmatched", "Matched")
    )
  ),
  silent = TRUE
)
dev.off()

matched_dep_gene <- NULL
matched_compare <- NULL
matched_compare_summary <- NULL

if (nrow(matched_df) >= 40 && all(table(matched_df$SampleGroup) >= 15)) {
  matched_model_df <- matched_df %>%
    mutate(
      SampleGroup = factor(as.character(SampleGroup), levels = MAIN_GROUPS),
      Sex = factor(Sex),
      Country = factor(Country),
      Age = as.numeric(Age),
      Education = as.numeric(Education)
    )

  matched_fit <- run_limma_dep_review(
    dat = matched_model_df,
    protein_cols = dep_protein_cols,
    annot_tbl = annot_tbl,
    formula_str = "~ SampleGroup + Age + Sex + Country + Education",
    coef_name = "SampleGroupAD",
    model_name = "MatchIt_nearest_neighbor_matched_DEP"
  )

  matched_dep_gene <- matched_fit$dep_gene
  matched_compare <- compare_to_main_review(
    DEP_gene,
    matched_dep_gene,
    "MatchIt_nearest_neighbor_matched_DEP"
  )
  matched_compare_summary <- summarize_comparison_review(
    matched_compare,
    "MatchIt_nearest_neighbor_matched_DEP"
  )

  safe_write_csv(matched_dep_gene, file.path(out_dir, "matched_DEP_gene_collapsed.csv"))
  safe_write_csv(matched_compare, file.path(out_dir, "main_vs_matched_DEP_gene_comparison.csv"))
  safe_write_csv(matched_compare_summary, file.path(out_dir, "main_vs_matched_DEP_summary.csv"))
} else {
  matched_compare_summary <- tibble(
    sensitivity = "MatchIt_nearest_neighbor_matched_DEP",
    warning = "Matched sample too small for stable DEP model using current minimum thresholds."
  )
  safe_write_csv(matched_compare_summary, file.path(out_dir, "main_vs_matched_DEP_summary.csv"))
}

###############################################################################

# 05. Sensitivity analysis 2: biomarker-consistent clinical AD/CN subset
###############################################################################

ptau_candidates <- c("p.tau217", "p-tau217", "p_tau217", "ptau217")
ab_candidates <- c("ratio.AB42.40", "ratio AB42/40", "ratio_AB42_40", "AB42_40")
ptau_col <- ptau_candidates[ptau_candidates %in% names(review_df)][1]
ab_col <- ab_candidates[ab_candidates %in% names(review_df)][1]
if (length(ptau_col) == 0) ptau_col <- NA_character_
if (length(ab_col) == 0) ab_col <- NA_character_

biomarker_subset_df <- NULL
biomarker_thresholds <- NULL
biomarker_dep_gene <- NULL
biomarker_compare <- NULL
biomarker_compare_summary <- NULL
biomarker_subset_summary <- NULL

if (!is.na(ptau_col) && !is.na(ab_col)) {
  cn_ref <- review_df %>% filter(SampleGroup == "CN")
  ptau_threshold <- quantile(as.numeric(cn_ref[[ptau_col]]), 0.75, na.rm = TRUE, names = FALSE)
  ab_threshold <- quantile(as.numeric(cn_ref[[ab_col]]), 0.25, na.rm = TRUE, names = FALSE)
  
  biomarker_thresholds <- tibble(
    ptau_col = ptau_col,
    abeta_ratio_col = ab_col,
    ptau217_high_threshold_source = "CN 75th percentile",
    abeta_ratio_low_threshold_source = "CN 25th percentile",
    ptau217_high_threshold = ptau_threshold,
    abeta_ratio_low_threshold = ab_threshold,
    definition_AD_biomarker_compatible = paste0(ptau_col, " >= CN P75 AND ", ab_col, " <= CN P25"),
    definition_CN_exclusion = paste0("CN excluded if ", ptau_col, " >= CN P75 AND ", ab_col, " <= CN P25")
  )
  
  biomarker_subset_df <- review_df %>%
    mutate(
      ptau_high = as.numeric(.data[[ptau_col]]) >= ptau_threshold,
      abeta_low = as.numeric(.data[[ab_col]]) <= ab_threshold,
      biomarker_AD_like = ptau_high & abeta_low,
      biomarker_consistent_group = case_when(
        SampleGroup == "AD" & biomarker_AD_like ~ TRUE,
        SampleGroup == "CN" & !biomarker_AD_like ~ TRUE,
        TRUE ~ FALSE
      )
    ) %>%
    filter(biomarker_consistent_group) %>%
    filter(complete.cases(Age, Sex, Country, Education)) %>%
    mutate(
      SampleGroup = factor(as.character(SampleGroup), levels = MAIN_GROUPS),
      Sex = factor(Sex),
      Country = factor(Country)
    )
  
  biomarker_subset_summary <- review_df %>%
    mutate(
      ptau_high = as.numeric(.data[[ptau_col]]) >= ptau_threshold,
      abeta_low = as.numeric(.data[[ab_col]]) <= ab_threshold,
      biomarker_AD_like = ptau_high & abeta_low,
      included_biomarker_consistent = case_when(
        SampleGroup == "AD" & biomarker_AD_like ~ TRUE,
        SampleGroup == "CN" & !biomarker_AD_like ~ TRUE,
        TRUE ~ FALSE
      )
    ) %>%
    group_by(SampleGroup, included_biomarker_consistent) %>%
    summarise(n = n(), .groups = "drop")
  
  if (nrow(biomarker_subset_df) >= 40 && all(table(biomarker_subset_df$SampleGroup) >= 15)) {
    biomarker_fit <- run_limma_dep_review(
      dat = biomarker_subset_df,
      protein_cols = dep_protein_cols,
      annot_tbl = annot_tbl,
      formula_str = "~ SampleGroup + Age + Sex + Country + Education",
      coef_name = "SampleGroupAD",
      model_name = "biomarker_consistent_DEP"
    )
    
    biomarker_dep_gene <- biomarker_fit$dep_gene
    biomarker_compare <- compare_to_main_review(DEP_gene, biomarker_dep_gene, "biomarker_consistent_DEP")
    biomarker_compare_summary <- summarize_comparison_review(biomarker_compare, "biomarker_consistent_DEP")
    
    safe_write_csv(biomarker_dep_gene, file.path(out_dir, "biomarker_consistent_DEP_gene_collapsed.csv"))
    safe_write_csv(biomarker_compare, file.path(out_dir, "main_vs_biomarker_consistent_DEP_gene_comparison.csv"))
    safe_write_csv(biomarker_compare_summary, file.path(out_dir, "main_vs_biomarker_consistent_DEP_summary.csv"))
  } else {
    biomarker_compare_summary <- tibble(
      sensitivity = "biomarker_consistent_DEP",
      warning = "Biomarker-consistent subset too small for stable DEP model using current minimum thresholds."
    )
  }
  
  safe_write_csv(biomarker_thresholds, file.path(out_dir, "biomarker_consistent_thresholds.csv"))
  safe_write_csv(biomarker_subset_summary, file.path(out_dir, "biomarker_consistent_sample_summary.csv"))
  
  id_keep <- intersect(c("SampleId", "SampleID", "SubjectID", "ID"), names(biomarker_subset_df))[1]
  sample_cols <- c(id_keep, "SampleGroup", "Country", "Sex", "Age", "Education", ptau_col, ab_col, "biomarker_AD_like")
  sample_cols <- sample_cols[!is.na(sample_cols) & sample_cols %in% names(biomarker_subset_df)]
  safe_write_csv(
    biomarker_subset_df %>% select(all_of(sample_cols)),
    file.path(out_dir, "biomarker_consistent_sample_ids.csv")
  )
} else {
  biomarker_thresholds <- tibble(
    warning = "p-tau217 and/or Aβ42/40 ratio columns were not detected. Biomarker-consistent sensitivity was skipped.",
    detected_ptau_col = ptau_col,
    detected_abeta_ratio_col = ab_col
  )
  biomarker_subset_summary <- biomarker_thresholds
  biomarker_compare_summary <- biomarker_thresholds %>% mutate(sensitivity = "biomarker_consistent_DEP")
  safe_write_csv(biomarker_thresholds, file.path(out_dir, "biomarker_consistent_skipped.csv"))
}

###############################################################################
# 06. Integrated manuscript-ready workbooks and summary
###############################################################################

sensitivity_summary <- bind_rows(
  matched_compare_summary,
  biomarker_compare_summary
)

manuscript_recommendation_notes <- tibble(
  section = c(
    "Cohort composition table",
    "MatchIt matched DEP sensitivity",
    "Biomarker-consistent DEP sensitivity",
    "Recommended Results wording",
    "Recommended Methods wording"
  ),
  text = c(
    "Use Supplementary_Table_1_Diagnostic_demographic_composition_by_recruitment_context.xlsx as the manuscript-ready Supplementary Table 1 for cohort composition by recruitment context.",
    "Report as an internal sensitivity analysis assessing preservation of the clinical AD-associated DEP signature after nearest-neighbor propensity-score matching on age, sex and education. Do not call it validation.",
    "Report as an internal biomarker-enriched sensitivity analysis. Because thresholds are data-driven by default, describe the exact threshold rule used.",
    "Because age, education and recruitment context differed across diagnostic groups, primary and sensitivity models adjusted for these variables and additional internal analyses evaluated robustness to matching and biomarker-compatible subgroup definitions.",
    "State that country/site are categorical recruitment-context variables. Avoid treating country as numeric or ordinal. Continuous descriptive comparisons use Wilcoxon rank-sum tests; categorical comparisons use chi-square/Fisher exact tests."
  )
)

safe_write_csv(sensitivity_summary, file.path(out_dir, "supplementary_table_2_sensitivity_summary.csv"))
safe_write_csv(manuscript_recommendation_notes, file.path(out_dir, "manuscript_recommendation_notes.csv"))

write_xlsx_list(
  list(
    Notes = manuscript_recommendation_notes,
    Sensitivity_summary = sensitivity_summary,
    Matched_sample_summary = matched_summary,
    Balance_before_after = balance_tbl,
    Main_vs_matched = matched_compare,
    Matched_DEP = matched_dep_gene,
    Biomarker_thresholds = biomarker_thresholds,
    Biomarker_subset_summary = biomarker_subset_summary,
    Main_vs_biomarker = biomarker_compare,
    Biomarker_DEP = biomarker_dep_gene
  ),
  file.path(out_dir, "Supplementary_Table_2_Internal_sensitivity_DEP_signature.xlsx")
)

save(
  composition_consolidated,
  composition_tests,
  variable_detection_map,
  country_group_counts,
  site_group_counts,
  matched_df,
  matched_summary,
  balance_tbl,
  matched_dep_gene,
  matched_compare,
  matched_compare_summary,
  biomarker_thresholds,
  biomarker_subset_df,
  biomarker_dep_gene,
  biomarker_compare,
  biomarker_compare_summary,
  sensitivity_summary,
  file = file.path(out_dir, "manuscript_supplementary_sensitivity_workspace.RData")
)

message("Manuscript-ready supplementary tables and sensitivity analyses complete.")
message("Outputs saved to: ", out_dir)
message("Main manuscript-ready Supplementary Table 1: ", file.path(out_dir, "Supplementary_Table_1_Diagnostic_demographic_composition_by_recruitment_context.xlsx"))
message("Supporting composition workbook: ", file.path(out_dir, "Supplementary_Table_Country_Site_Composition_SUPPORTING.xlsx"))
message("Main sensitivity workbook: ", file.path(out_dir, "Supplementary_Table_2_Internal_sensitivity_DEP_signature.xlsx"))
###############################################################################
# END
###############################################################################


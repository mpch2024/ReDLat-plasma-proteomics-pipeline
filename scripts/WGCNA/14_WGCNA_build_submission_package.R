###############################################################################
# ReDLat plasma proteomics — WGCNA workflow
# 14. Build the submission package
# Requires: outputs from Scripts 10–13
# Produces: figures, tables, Source Data and Supplementary Data package
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

required_pkgs <- c(
  "dplyr", "readr", "tibble", "purrr", "stringr", "openxlsx", "zip"
)
missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs) > 0L) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
       ". Run renv::restore() before this script.", call. = FALSE)
}
invisible(lapply(required_pkgs, library, character.only = TRUE))

run_submission_package <- function() {

###############################################################################
# 00. PATHS AND CONSTANTS
###############################################################################

BASE_DIR <- WGCNA_CONFIG$project_root

SCRIPT_DIR <- file.path(WGCNA_CONFIG$project_root, "scripts", "WGCNA")
WGCNA_ROOT <- WGCNA_CONFIG$result_root
FIGURE_ROOT <- file.path(WGCNA_CONFIG$publication_root, "figures")

SCRIPT19_FINAL <- file.path(SCRIPT_DIR, "12_WGCNA_generate_supplementary_tables.R")
SCRIPT20_FINAL <- file.path(SCRIPT_DIR, "13_WGCNA_generate_manuscript_text.R")

S10 <- file.path(WGCNA_ROOT, "01_input")
S11 <- file.path(WGCNA_ROOT, "02_network")
S12 <- file.path(WGCNA_ROOT, "03_modules")
S13 <- file.path(WGCNA_ROOT, "04_module_traits")
S13B <- file.path(WGCNA_ROOT, "05_sensitivity")
S14 <- file.path(WGCNA_ROOT, "06_stability")
S14B <- file.path(WGCNA_ROOT, "07_biomarker_fdr")
S15B <- file.path(WGCNA_ROOT, "08_preservation")
S16 <- file.path(WGCNA_ROOT, "09_network_quality")

S17_FINAL <- file.path(FIGURE_ROOT, "main_figure_3")
S18_FINAL <- file.path(FIGURE_ROOT, "extended_data")
S19_FINAL <- file.path(WGCNA_CONFIG$publication_root, "supplementary_tables")
S20_FINAL <- file.path(WGCNA_CONFIG$publication_root, "manuscript_text")

OUT_ROOT <- file.path(WGCNA_CONFIG$publication_root, "submission_package")
PACKAGE_ROOT <- file.path(OUT_ROOT, "WGCNA_NatureAging_PACKAGE")
FIG_MAIN_DIR <- file.path(PACKAGE_ROOT, "Figures", "Main")
FIG_EXT_DIR <- file.path(PACKAGE_ROOT, "Figures", "Extended_Data")
TABLE_DIR <- file.path(PACKAGE_ROOT, "Tables")
SOURCE_DATA_DIR <- file.path(PACKAGE_ROOT, "Source_Data")
SUPP_DATA_DIR <- file.path(PACKAGE_ROOT, "Supplementary_Data")
TEXT_DIR <- file.path(PACKAGE_ROOT, "Manuscript_Text_Optional")
LEGEND_DIR <- file.path(PACKAGE_ROOT, "Figure_Legends")
AUDIT_DIR <- file.path(PACKAGE_ROOT, "Audit")
LEGACY_DIR <- file.path(PACKAGE_ROOT, "Legacy_Do_Not_Use")

invisible(lapply(
  c(
    OUT_ROOT, PACKAGE_ROOT, FIG_MAIN_DIR, FIG_EXT_DIR, TABLE_DIR,
    SOURCE_DATA_DIR, SUPP_DATA_DIR, TEXT_DIR, LEGEND_DIR, AUDIT_DIR, LEGACY_DIR
  ),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

OUTPUTS <- list(
  supplementary_tables = file.path(TABLE_DIR, "Supplementary_Tables_WGCNA_NatureAging_CLEAN.xlsx"),
  source_fig3 = file.path(SOURCE_DATA_DIR, "SourceData_Fig3_WGCNA_NatureAging.xlsx"),
  source_ed7 = file.path(SOURCE_DATA_DIR, "SourceData_ExtendedData_Figure_07_WGCNA.xlsx"),
  source_ed8 = file.path(SOURCE_DATA_DIR, "SourceData_ExtendedData_Figure_08_WGCNA.xlsx"),
  source_ed9 = file.path(SOURCE_DATA_DIR, "SourceData_ExtendedData_Figure_09_WGCNA.xlsx"),
  source_ed10 = file.path(SOURCE_DATA_DIR, "SourceData_ExtendedData_Figure_10_WGCNA.xlsx"),
  source_ed11 = file.path(SOURCE_DATA_DIR, "SourceData_ExtendedData_Figure_11_WGCNA.xlsx"),
  supp_data4 = file.path(SUPP_DATA_DIR, "Supplementary_Data_4_WGCNA_Complete_Module_Results.xlsx"),
  supp_data5 = file.path(SUPP_DATA_DIR, "Supplementary_Data_5_WGCNA_Complete_Association_and_Context_Results.xlsx"),
  supp_data6 = file.path(SUPP_DATA_DIR, "Supplementary_Data_6_WGCNA_Complete_Preservation_and_Diagnostics.xlsx"),
  global_index_csv = file.path(PACKAGE_ROOT, "WGCNA_Submission_Package_Index.csv"),
  global_index_xlsx = file.path(PACKAGE_ROOT, "WGCNA_Submission_Package_Index.xlsx"),
  source_registry = file.path(AUDIT_DIR, "WGCNA_Source_Registry_FINAL.csv"),
  combined_legends_csv = file.path(LEGEND_DIR, "WGCNA_Final_Figure_Legends.csv"),
  combined_legends_md = file.path(LEGEND_DIR, "WGCNA_Final_Figure_Legends.md"),
  preaudit = file.path(AUDIT_DIR, "WGCNA_PrePackage_Audit.csv"),
  package_manifest = file.path(AUDIT_DIR, "WGCNA_Package_Manifest.csv"),
  readme = file.path(PACKAGE_ROOT, "README_WGCNA_NatureAging_FINAL.txt"),
  session_info = file.path(AUDIT_DIR, "sessionInfo_WGCNA_submission_package.txt"),
  zip = file.path(OUT_ROOT, "WGCNA_NatureAging_PACKAGE.zip")
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

find_one <- function(dir, pattern, required = FALSE, label = "source") {
  if (!dir.exists(dir)) {
    if (required) stop("Directory not found for ", label, ": ", dir, call. = FALSE)
    return(NA_character_)
  }

  hit <- list.files(dir, pattern = pattern, full.names = TRUE, recursive = TRUE)
  if (length(hit) == 0) {
    if (required) stop("No file found for ", label, " using pattern ", pattern, call. = FALSE)
    return(NA_character_)
  }

  normalizePath(hit[1], winslash = "/", mustWork = TRUE)
}

relative_path <- function(path) {
  root <- normalizePath(BASE_DIR, winslash = "/", mustWork = FALSE)
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  prefix <- paste0(root, "/")
  if (startsWith(path, prefix)) substring(path, nchar(prefix) + 1L) else path
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

read_csv_safe <- function(path, required = TRUE, label = basename(path)) {
  if (is.na(path) || !file.exists(path)) {
    if (required) stop("Missing required CSV: ", label, "\n", path, call. = FALSE)
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

  tbl[[module_col]] <- sub("^ME", "", as.character(tbl[[module_col]]))

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

copy_required <- function(from, to_dir, new_name = NULL) {
  if (!file.exists(from)) {
    stop("Required package file is missing:\n", from, call. = FALSE)
  }
  dir.create(to_dir, recursive = TRUE, showWarnings = FALSE)
  target_name <- if (is.null(new_name)) basename(from) else new_name
  target <- file.path(to_dir, target_name)
  ok <- file.copy(from, target, overwrite = TRUE)
  if (!ok || !file.exists(target) || file.info(target)$size <= 0) {
    stop("Failed to copy package file:\n", from, call. = FALSE)
  }
  normalizePath(target, winslash = "/", mustWork = TRUE)
}

run_isolated_script <- function(script_path, label, required = TRUE) {
  if (!file.exists(script_path)) {
    if (required) {
      stop("Missing required script for ", label, ":\n", script_path, call. = FALSE)
    }
    return(list(status = "NOT_FOUND", error = NA_character_))
  }

  script_env <- new.env(parent = globalenv())

  result <- tryCatch(
    {
      sys.source(
        script_path,
        envir = script_env,
        chdir = FALSE,
        keep.source = TRUE
      )
      list(status = "PASS", error = NA_character_)
    },
    error = function(e) {
      if (required) {
        stop(label, " failed:\n", conditionMessage(e), call. = FALSE)
      }
      list(status = "OPTIONAL_FAIL", error = conditionMessage(e))
    }
  )

  result
}

nonempty_file <- function(path) {
  file.exists(path) &&
    is.finite(file.info(path)$size) &&
    file.info(path)$size > 0
}


###############################################################################
# 02. WORKBOOK STYLES AND GENERIC WRITER
###############################################################################

COL_TITLE <- "#CFC7B7"
COL_HEADER <- "#E9E5DC"
COL_NOTE <- "#F4F2EC"
COL_BORDER <- "#4A4A4A"

styles <- list(
  title = openxlsx::createStyle(
    fontName = "Arial", fontSize = 12, textDecoration = "bold",
    fgFill = COL_TITLE, fontColour = "#111111",
    halign = "left", valign = "center",
    border = "TopBottomLeftRight", borderColour = COL_BORDER
  ),
  note = openxlsx::createStyle(
    fontName = "Arial", fontSize = 9, textDecoration = "italic",
    fgFill = COL_NOTE, fontColour = "#555555", wrapText = TRUE,
    border = "TopBottomLeftRight", borderColour = COL_BORDER
  ),
  header = openxlsx::createStyle(
    fontName = "Arial", fontSize = 9, textDecoration = "bold",
    fgFill = COL_HEADER, fontColour = "#111111",
    halign = "center", valign = "center", wrapText = TRUE,
    border = "TopBottomLeftRight", borderColour = COL_BORDER
  ),
  body = openxlsx::createStyle(
    fontName = "Arial", fontSize = 8, valign = "top", wrapText = TRUE,
    border = "TopBottomLeftRight", borderColour = "#B7B7B7"
  )
)

write_sheet <- function(wb, sheet, title, note, data) {
  sheet <- safe_sheet_name(sheet, openxlsx::sheets(wb))
  data <- sanitize_table(add_module_identity(data))
  ncols <- max(2L, ncol(data))

  openxlsx::addWorksheet(wb, sheet, gridLines = FALSE)

  openxlsx::mergeCells(wb, sheet, cols = seq_len(ncols), rows = 1)
  openxlsx::writeData(wb, sheet, title, startRow = 1)
  openxlsx::addStyle(wb, sheet, styles$title, rows = 1, cols = seq_len(ncols), gridExpand = TRUE)
  openxlsx::setRowHeights(wb, sheet, 1, 26)

  openxlsx::mergeCells(wb, sheet, cols = seq_len(ncols), rows = 2)
  openxlsx::writeData(wb, sheet, note, startRow = 2)
  openxlsx::addStyle(wb, sheet, styles$note, rows = 2, cols = seq_len(ncols), gridExpand = TRUE)
  openxlsx::setRowHeights(wb, sheet, 2, 40)

  if (ncol(data) == 0L) data <- data.frame(Note = "No rows available.")

  openxlsx::writeData(
    wb, sheet, data,
    startRow = 4, startCol = 1,
    colNames = TRUE, rowNames = FALSE,
    withFilter = TRUE
  )
  openxlsx::addStyle(wb, sheet, styles$header, rows = 4, cols = seq_len(ncol(data)), gridExpand = TRUE)

  # Avoid applying cell-by-cell body styles to very large Supplementary Data
  # sheets, which can substantially increase memory use and workbook size.
  if (nrow(data) > 0 && nrow(data) <= 20000L) {
    openxlsx::addStyle(
      wb, sheet, styles$body,
      rows = 5:(4 + nrow(data)),
      cols = seq_len(ncol(data)),
      gridExpand = TRUE
    )
  }

  openxlsx::freezePane(wb, sheet, firstActiveRow = 5)
  openxlsx::setColWidths(wb, sheet, cols = seq_len(ncol(data)), widths = 16)

  long_cols <- which(
    grepl(
      "description|term|pathway|label|note|formula|source|interpretation|model",
      names(data),
      ignore.case = TRUE
    )
  )
  if (length(long_cols) > 0) {
    openxlsx::setColWidths(wb, sheet, cols = long_cols, widths = 34)
  }

  invisible(sheet)
}

build_workbook <- function(path, workbook_title, plan, workbook_note) {
  wb <- openxlsx::createWorkbook(creator = "Matías Pizarro")

  index_tbl <- purrr::map_dfr(
    plan,
    function(x) {
      tibble::tibble(
        Sheet = x$sheet,
        Description = x$title,
        Source = if (!is.null(x$path)) relative_path(x$path) else "Derived in Script 22"
      )
    }
  )

  write_sheet(
    wb, "Index", workbook_title, workbook_note, index_tbl
  )
  write_sheet(
    wb, "Module_Key",
    "WGCNA module key",
    "M1-M8 are stable visual identifiers and do not represent a ranking.",
    MODULE_KEY
  )

  for (x in plan) {
    dat <- if (!is.null(x$data)) {
      x$data
    } else {
      read_csv_safe(x$path, required = isTRUE(x$required), label = x$title)
    }

    write_sheet(wb, x$sheet, x$title, x$note, dat)
  }

  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)

  if (!file.exists(path) || file.info(path)$size <= 0) {
    stop("Workbook was not created correctly:\n", path, call. = FALSE)
  }

  invisible(path)
}

###############################################################################
# 03. SELF-HEALING PREFLIGHT
###############################################################################

table_source <- file.path(
  S19_FINAL,
  "Supplementary_Tables_WGCNA_NatureAging_CLEAN.xlsx"
)

script19_status <- "ALREADY_AVAILABLE"
if (!nonempty_file(table_source)) {
  result19 <- run_isolated_script(
    SCRIPT19_FINAL,
    "Script 19 supplementary-table generation",
    required = TRUE
  )
  script19_status <- result19$status
}

if (!nonempty_file(table_source)) {
  stop(
    "Script 19 completed without generating the required supplementary-table workbook:\n",
    table_source,
    call. = FALSE
  )
}

text_expected <- c(
  file.path(S20_FINAL, "WGCNA_Main_Methods_FINAL.txt"),
  file.path(S20_FINAL, "WGCNA_Main_Results_FINAL.txt"),
  file.path(S20_FINAL, "WGCNA_Figure3_Legend_FINAL.txt"),
  file.path(S20_FINAL, "WGCNA_Extended_Data_Legends_FINAL.csv")
)

script21_status <- "ALREADY_AVAILABLE"
script21_error <- NA_character_

if (!all(vapply(text_expected, nonempty_file, logical(1)))) {
  result21 <- run_isolated_script(
    SCRIPT20_FINAL,
    "Script 20 manuscript-text generation",
    required = FALSE
  )
  script21_status <- result21$status
  script21_error <- result21$error
}

preflight_status <- tibble::tibble(
  Component = c("Script 19 tables", "Script 20 editorial text"),
  Status = c(script19_status, script21_status),
  Required_for_package = c(TRUE, FALSE),
  Error = c(NA_character_, script21_error)
)
readr::write_csv(
  preflight_status,
  file.path(AUDIT_DIR, "WGCNA_Package_Preflight_Status.csv")
)

copy_required(
  table_source,
  TABLE_DIR,
  basename(OUTPUTS$supplementary_tables)
)

###############################################################################
# 04. SOURCE DATA FOR FIGURE 3
###############################################################################

fig3_source_dir <- file.path(S17_FINAL, "diagnostics", "source_data")
fig3_csvs <- if (dir.exists(fig3_source_dir)) {
  list.files(fig3_source_dir, pattern = "\\.csv$", full.names = TRUE)
} else {
  character()
}

if (length(fig3_csvs) == 0) {
  stop("No Figure 3 source-data CSVs were found.", call. = FALSE)
}

workspace_path <- file.path(
  S11,
  "workspace",
  "wgcna_core_collapsed_workspace.RData"
)
if (!file.exists(workspace_path)) {
  stop("Required Script 11 workspace is missing:\n", workspace_path, call. = FALSE)
}

workspace_env <- new.env(parent = baseenv())
loaded <- load(workspace_path, envir = workspace_env)
required_objects <- c("geneTree", "dynamicColors", "mergedColors")
missing_objects <- setdiff(required_objects, loaded)
if (length(missing_objects) > 0) {
  stop(
    "Script 11 workspace lacks dendrogram objects: ",
    paste(missing_objects, collapse = ", "),
    call. = FALSE
  )
}

geneTree <- get("geneTree", envir = workspace_env)
dynamicColors <- get("dynamicColors", envir = workspace_env)
mergedColors <- get("mergedColors", envir = workspace_env)

gene_labels <- geneTree$labels
if (is.null(gene_labels) || length(gene_labels) != length(mergedColors)) {
  gene_labels <- names(mergedColors)
}
if (is.null(gene_labels) || length(gene_labels) != length(mergedColors)) {
  gene_labels <- paste0("Gene_", seq_along(mergedColors))
}

dendro_order <- geneTree$order
dendrogram_source <- tibble::tibble(
  Dendrogram_position = seq_along(dendro_order),
  Original_gene_index = dendro_order,
  Gene_label = gene_labels[dendro_order],
  Initial_module = as.character(dynamicColors[dendro_order]),
  Merged_module = as.character(mergedColors[dendro_order])
) %>%
  dplyr::left_join(MODULE_KEY, by = c("Merged_module" = "Module"))

fig3_plan <- c(
  list(
    list(
      sheet = "Panel_a_dendrogram",
      title = "Figure 3a source data: dendrogram order and initial/final module assignments",
      note = "The dendrogram itself is reconstructed from the Script 11 hclust object; this sheet provides the ordered gene labels and module-color tracks.",
      data = dendrogram_source,
      path = NULL,
      required = TRUE
    )
  ),
  purrr::map(
    fig3_csvs,
    function(path) {
      list(
        sheet = tools::file_path_sans_ext(basename(path)),
        title = paste0("Figure 3 source data: ", tools::file_path_sans_ext(basename(path))),
        note = "Underlying numerical data exported by the definitive Script 17 v16n figure workflow.",
        path = path,
        data = NULL,
        required = TRUE
      )
    }
  )
)

build_workbook(
  OUTPUTS$source_fig3,
  "Source Data Fig. 3. WGCNA modular organization",
  fig3_plan,
  "Underlying data for all quantitative panels of the definitive Figure 3. M1-M8 are visual identifiers and not a ranking."
)

###############################################################################
# 05. SOURCE DATA FOR EXTENDED DATA FIGS. 7-11
###############################################################################

ext_source_dir <- file.path(S18_FINAL, "source_data")
if (!dir.exists(ext_source_dir)) {
  stop("Extended Data source-data directory was not found:\n", ext_source_dir, call. = FALSE)
}

ext_outputs <- list(
  `07` = OUTPUTS$source_ed7,
  `08` = OUTPUTS$source_ed8,
  `09` = OUTPUTS$source_ed9,
  `10` = OUTPUTS$source_ed10,
  `11` = OUTPUTS$source_ed11
)

for (fig_no in names(ext_outputs)) {
  files <- list.files(
    ext_source_dir,
    pattern = paste0("^Extended_Data_Figure_", fig_no, ".*\\.csv$"),
    full.names = TRUE
  )

  if (length(files) == 0) {
    stop("No source CSVs found for Extended Data Fig. ", fig_no, ".", call. = FALSE)
  }

  plan <- purrr::map(
    files,
    function(path) {
      list(
        sheet = tools::file_path_sans_ext(basename(path)),
        title = paste0(
          "Extended Data Fig. ", as.integer(fig_no),
          " source data: ", tools::file_path_sans_ext(basename(path))
        ),
        note = "Underlying numerical data exported by the definitive Script 18 v6 figure workflow.",
        path = path,
        data = NULL,
        required = TRUE
      )
    }
  )

  build_workbook(
    ext_outputs[[fig_no]],
    paste0("Source Data Extended Data Fig. ", as.integer(fig_no), ". WGCNA"),
    plan,
    paste0(
      "Underlying data for Extended Data Fig. ", as.integer(fig_no),
      ". M1-M8 are stable visual identifiers and do not represent a ranking."
    )
  )
}

###############################################################################
# 06. SUPPLEMENTARY DATA 4 — COMPLETE MODULE RESULTS
###############################################################################

SD4_PLAN <- list(
  list(
    sheet = "Gene_SOMAmer_map",
    title = "Outcome-independent gene-SOMAmer map",
    note = "One representative SOMAmer per gene selected without diagnosis or other outcome information.",
    path = file.path(S10, "tables", "outcome_independent_gene_somamer_map.csv"),
    data = NULL, required = TRUE
  ),
  list(
    sheet = "Gene_module_annotation",
    title = "Complete gene-level module, kME and DEP annotation",
    note = "Complete 9,638-gene WGCNA annotation with analyte-matched and canonical primary DEP definitions.",
    path = file.path(S12, "tables", "gene_module_kME_DEP_dual_annotation.csv"),
    data = NULL, required = TRUE
  ),
  list(
    sheet = "Full_assigned_kME",
    title = "Complete assigned-module kME results",
    note = "One row per gene with assigned-module membership metrics.",
    path = file.path(S12, "tables", "full_kME_assigned_module_long.csv"),
    data = NULL, required = TRUE
  ),
  list(
    sheet = "Module_biology",
    title = "Integrated module biology summary",
    note = "Module size, DEP burden, hubs and functional annotation.",
    path = file.path(S12, "tables", "integrated_module_biology_summary.csv"),
    data = NULL, required = TRUE
  ),
  list(
    sheet = "DEP_burden",
    title = "Dual-definition module DEP burden",
    note = "Analyte-matched and canonical primary DEP overlap are kept separate.",
    path = file.path(S12, "tables", "module_DEP_burden_dual_definition.csv"),
    data = NULL, required = TRUE
  ),
  list(
    sheet = "DEP_overrepresentation",
    title = "Dual-definition DEP overrepresentation",
    note = "Module-level overrepresentation tests under both DEP definitions.",
    path = file.path(S12, "tables", "module_DEP_overrepresentation_dual_definition.csv"),
    data = NULL, required = TRUE
  ),
  list(
    sheet = "Hub_summary",
    title = "Module hub summary",
    note = "Module-level kME and high-membership hub counts.",
    path = file.path(S12, "tables", "module_hub_summary.csv"),
    data = NULL, required = TRUE
  ),
  list(
    sheet = "Top_hubs",
    title = "Top hubs across all modules",
    note = "Top proteins ranked by absolute module membership.",
    path = file.path(S12, "tables", "top_hubs_all_modules_combined.csv"),
    data = NULL, required = TRUE
  ),
  list(
    sheet = "Enrichment_summary",
    title = "Enrichment summary by module",
    note = "Counts and leading terms across GO BP, KEGG and Reactome.",
    path = file.path(S12, "tables", "enrichment_summary_by_module.csv"),
    data = NULL, required = TRUE
  ),
  list(
    sheet = "Enrichment_all",
    title = "All FDR-significant module enrichment terms",
    note = "Complete significant enrichment results across annotation libraries.",
    path = file.path(S12, "tables", "enrichment_significant_terms_all_modules.csv"),
    data = NULL, required = TRUE
  )
)

build_workbook(
  OUTPUTS$supp_data4,
  "Supplementary Data 4. Complete WGCNA module results",
  SD4_PLAN,
  "Complete high-dimensional module annotation, hub and enrichment results supporting Supplementary Tables 20-22 and Figure 3."
)

###############################################################################
# 07. SUPPLEMENTARY DATA 5 — ASSOCIATIONS AND CONTEXT
###############################################################################

optional_loco_raw <- find_one(
  file.path(S14, "tables", "loco"),
  "loco.*correlation.*long.*\\.csv$",
  required = FALSE
)
optional_loso_raw <- find_one(
  file.path(S14, "tables", "within_country_loso"),
  "within.*country.*loso.*correlation.*long.*\\.csv$",
  required = FALSE
)
optional_downsampling_raw <- find_one(
  file.path(S14, "tables", "downsampling"),
  "downsampling.*correlation.*long.*\\.csv$",
  required = FALSE
)

SD5_PLAN <- list(
  list(
    sheet = "Module_trait_128",
    title = "Complete 128-test module-trait family",
    note = "BH-FDR correction was applied across all eight modules and 16 traits.",
    path = file.path(S13, "tables", "correlations", "module_trait_results_long.csv"),
    data = NULL, required = TRUE
  ),
  list(
    sheet = "Adjusted_continuous_48",
    title = "Complete 48-model adjusted continuous-outcome family",
    note = "Covariate-adjusted continuous clinical and biomarker outcomes.",
    path = file.path(S13, "tables", "regression", "adjusted_module_models.csv"),
    data = NULL, required = TRUE
  ),
  list(
    sheet = "Diagnosis_HC3_8",
    title = "Eight adjusted diagnosis models with HC3 inference",
    note = "Module eigengene differences by diagnosis with robust standard errors.",
    path = file.path(S13B, "tables", "diagnosis_robustness", "adjusted_diagnosis_module_models_HC3.csv"),
    data = NULL, required = TRUE
  ),
  list(
    sheet = "Biomarker_family32",
    title = "Complete corrected 32-model plasma biomarker family",
    note = "Authoritative biomarker inference from Script 14b.",
    path = file.path(S14B, "tables", "full", "full_all_modules_log_HC3_family32.csv"),
    data = NULL, required = TRUE
  ),
  list(
    sheet = "Biomarker_focal_stability",
    title = "Focal-module biomarker stability",
    note = "M1/green, M2/blue and M3/brown extracted only after correction across the complete 32-model family.",
    path = file.path(S14B, "tables", "focal_modules", "integrated_focal_log_HC3_stability_family32.csv"),
    data = NULL, required = TRUE
  ),
  list(
    sheet = "Country_context",
    title = "Categorical country context effects",
    note = "Country contribution to eigengene variation.",
    path = file.path(S13, "tables", "context", "module_recruitment_context_effect_summary.csv"),
    data = NULL, required = TRUE
  ),
  list(
    sheet = "Nested_site_context",
    title = "Corrected site-within-country effects",
    note = "Incremental site-within-country variation after country adjustment.",
    path = file.path(S13B, "tables", "context", "corrected_nested_site_all_samples.csv"),
    data = NULL, required = TRUE
  ),
  list(
    sheet = "LOCO_summary",
    title = "LOCO module-trait stability summary",
    note = "Association stability with fixed module definitions.",
    path = file.path(S14, "tables", "loco", "loco_country_summary_by_module_trait.csv"),
    data = NULL, required = TRUE
  ),
  list(
    sheet = "WC_LOSO_summary",
    title = "Corrected within-country LOSO summary",
    note = "Only valid site deletions that preserve the corresponding country are interpreted.",
    path = file.path(S14, "tables", "within_country_loso", "within_country_loso_summary_by_module_trait.csv"),
    data = NULL, required = TRUE
  ),
  list(
    sheet = "Downsampling_summary",
    title = "Balanced-downsampling module-trait summary",
    note = "Five hundred reproducible country-by-diagnosis balanced samples.",
    path = file.path(S14, "tables", "downsampling", "balanced_downsampling_summary_by_module_trait.csv"),
    data = NULL, required = TRUE
  ),
  list(
    sheet = "Downsampling_iterations",
    title = "Balanced-downsampling iteration summary",
    note = "Iteration-level sample and stability metrics.",
    path = file.path(S14, "tables", "downsampling", "balanced_downsampling_summary_by_iteration.csv"),
    data = NULL, required = TRUE
  ),
  list(
    sheet = "Diagnosis_stability",
    title = "Adjusted diagnosis HC3 stability",
    note = "Full, LOCO, valid LOSO and balanced-downsampling diagnosis-model stability.",
    path = file.path(S14, "tables", "model_stability", "diagnosis", "integrated_diagnosis_HC3_stability_summary.csv"),
    data = NULL, required = TRUE
  ),
  list(
    sheet = "Core_module_stability",
    title = "Core-module association stability",
    note = "Integrated correlation-stability metrics for focal modules.",
    path = file.path(S14, "tables", "core_modules", "core_modules_correlation_stability_summary.csv"),
    data = NULL, required = TRUE
  )
)

if (!is.na(optional_loco_raw)) {
  SD5_PLAN <- c(
    SD5_PLAN,
    list(list(
      sheet = "LOCO_raw_long",
      title = "Complete LOCO module-trait results",
      note = "Optional full long-form LOCO correlation output.",
      path = optional_loco_raw, data = NULL, required = FALSE
    ))
  )
}
if (!is.na(optional_loso_raw)) {
  SD5_PLAN <- c(
    SD5_PLAN,
    list(list(
      sheet = "WC_LOSO_raw_long",
      title = "Complete valid within-country LOSO module-trait results",
      note = "Optional full long-form within-country site-deletion output.",
      path = optional_loso_raw, data = NULL, required = FALSE
    ))
  )
}
if (!is.na(optional_downsampling_raw)) {
  SD5_PLAN <- c(
    SD5_PLAN,
    list(list(
      sheet = "Downsampling_raw_long",
      title = "Complete balanced-downsampling module-trait results",
      note = "Optional full long-form balanced-downsampling output.",
      path = optional_downsampling_raw, data = NULL, required = FALSE
    ))
  )
}

build_workbook(
  OUTPUTS$supp_data5,
  "Supplementary Data 5. Complete WGCNA association and context results",
  SD5_PLAN,
  "Complete module-trait, adjusted-model, recruitment-context and association-stability results. Association stability is distinct from structural preservation."
)

###############################################################################
# 08. SUPPLEMENTARY DATA 6 — PRESERVATION AND DIAGNOSTICS
###############################################################################

SD6_PLAN <- list(
  list(
    sheet = "Fixed_gene_manifest",
    title = "Fixed 4,239-gene preservation manifest",
    note = "The identical gene set was used in every country and site comparison.",
    path = file.path(S15B, "tables", "input", "fixed_gene_manifest_all_modules.csv"),
    data = NULL, required = TRUE
  ),
  list(
    sheet = "Fixed_gene_summary",
    title = "Fixed-gene manifest summary by module",
    note = "Gene retention and reproducible subsampling by module.",
    path = file.path(S15B, "tables", "input", "fixed_gene_manifest_summary_by_module.csv"),
    data = NULL, required = TRUE
  ),
  list(
    sheet = "Site_run_plan",
    title = "Reciprocal within-country site run plan",
    note = "Chile and Colombia reciprocal reference-test comparisons.",
    path = file.path(S15B, "tables", "input", "reciprocal_within_country_site_run_plan.csv"),
    data = NULL, required = TRUE
  ),
  list(
    sheet = "Preservation_summary",
    title = "Definitive fixed-gene preservation summary",
    note = "Country and site minimum Zsummary and preservation counts.",
    path = file.path(S15B, "tables", "fixed_geneset_structural_preservation_definitive_summary.csv"),
    data = NULL, required = TRUE
  ),
  list(
    sheet = "Country_primary",
    title = "Country-held-out primary preservation statistics",
    note = "Primary modulePreservation statistics for all 40 country-module tests.",
    path = file.path(S15B, "tables", "country", "country_fixed_geneset_primary_statistics.csv"),
    data = NULL, required = TRUE
  ),
  list(
    sheet = "Country_all_statistics",
    title = "Country-held-out complete WGCNA preservation statistics",
    note = "Complete raw preservation-statistic output.",
    path = file.path(S15B, "tables", "country", "country_fixed_geneset_all_WGCNA_statistics_raw.csv"),
    data = NULL, required = TRUE
  ),
  list(
    sheet = "Site_primary",
    title = "Reciprocal site primary preservation statistics",
    note = "Primary modulePreservation statistics for all 32 reciprocal site-module tests.",
    path = file.path(S15B, "tables", "site", "site_fixed_geneset_primary_statistics.csv"),
    data = NULL, required = TRUE
  ),
  list(
    sheet = "Site_all_statistics",
    title = "Reciprocal site complete WGCNA preservation statistics",
    note = "Complete raw preservation-statistic output.",
    path = file.path(S15B, "tables", "site", "site_fixed_geneset_all_WGCNA_statistics_raw.csv"),
    data = NULL, required = TRUE
  ),
  list(
    sheet = "Soft_threshold",
    title = "Soft-threshold scan",
    note = "Scale-free topology fit and mean connectivity across candidate powers.",
    path = file.path(S16, "tables", "soft_threshold", "soft_threshold_scan_clean.csv"),
    data = NULL, required = TRUE
  ),
  list(
    sheet = "Gene_network_quality",
    title = "Complete gene-level network-quality metrics",
    note = "Connectivity, intramodular connectivity and assigned-module membership for all genes.",
    path = file.path(S16, "tables", "connectivity", "gene_level_network_quality_metrics.csv"),
    data = NULL, required = TRUE
  ),
  list(
    sheet = "Module_network_quality",
    title = "Integrated module network-quality summary",
    note = "Module-level separation, kME, modularity contribution and preservation.",
    path = file.path(S16, "tables", "integration", "integrated_module_network_quality_summary.csv"),
    data = NULL, required = TRUE
  ),
  list(
    sheet = "Module_separation",
    title = "Within-versus-external module separation",
    note = "Adjacency and TOM separation metrics.",
    path = file.path(S16, "tables", "modularity", "module_internal_external_separation_summary.csv"),
    data = NULL, required = TRUE
  ),
  list(
    sheet = "Weighted_Q",
    title = "Exact post hoc weighted modularity Q",
    note = "Descriptive global modularity metric; not used to optimize WGCNA.",
    path = file.path(S16, "tables", "modularity", "posthoc_weighted_modularity_Q.csv"),
    data = NULL, required = TRUE
  ),
  list(
    sheet = "Q_contributions",
    title = "Module-specific contributions to weighted modularity",
    note = "Module contributions to the exact post hoc weighted Q.",
    path = file.path(S16, "tables", "modularity", "module_specific_modularity_contributions.csv"),
    data = NULL, required = TRUE
  ),
  list(
    sheet = "kME_quality",
    title = "Module kME-quality summary",
    note = "Module membership and intramodular-connectivity concordance.",
    path = file.path(S16, "tables", "kME", "module_kME_quality_summary.csv"),
    data = NULL, required = TRUE
  ),
  list(
    sheet = "Eigengene_pairs",
    title = "Pairwise module eigengene correlations",
    note = "Final eigengene correlations and merge-threshold audit.",
    path = file.path(S16, "tables", "eigengenes", "module_eigengene_pairwise_correlations.csv"),
    data = NULL, required = TRUE
  ),
  list(
    sheet = "Eigengene_matrix",
    title = "Module eigengene correlation matrix",
    note = "Complete final eigengene correlation matrix.",
    path = file.path(S16, "tables", "eigengenes", "module_eigengene_correlation_matrix.csv"),
    data = NULL, required = TRUE
  ),
  list(
    sheet = "Adjacency_mixing",
    title = "Adjacency module-block mixing",
    note = "Mean adjacency across module blocks.",
    path = file.path(S16, "tables", "matrix_mixing", "adjacency_module_block_mixing.csv"),
    data = NULL, required = TRUE
  ),
  list(
    sheet = "TOM_mixing",
    title = "Topological-overlap module-block mixing",
    note = "Mean TOM across module blocks.",
    path = file.path(S16, "tables", "matrix_mixing", "TOM_module_block_mixing.csv"),
    data = NULL, required = TRUE
  )
)

build_workbook(
  OUTPUTS$supp_data6,
  "Supplementary Data 6. Complete WGCNA preservation and network diagnostics",
  SD6_PLAN,
  "Complete fixed-gene preservation and network-quality outputs. Exact weighted modularity is descriptive and is not a validation threshold."
)

###############################################################################
# 09. COPY FINAL FIGURES
###############################################################################

main_base <- "Figure3_WGCNA_COMPACT_CLASSIC_v16n_FINAL"
main_figure_files <- c(
  file.path(S17_FINAL, "Illustrator_ready", "master_vectors", paste0(main_base, ".pdf")),
  file.path(S17_FINAL, "Illustrator_ready", "master_vectors", paste0(main_base, ".svg")),
  file.path(S17_FINAL, "preview", paste0(main_base, ".png")),
  file.path(S17_FINAL, "submission", paste0(main_base, ".tiff"))
)

copied_main <- vapply(
  main_figure_files,
  copy_required,
  character(1),
  to_dir = FIG_MAIN_DIR
)

extended_bases <- c(
  "Extended_Data_Figure_07_WGCNA_network_quality",
  "Extended_Data_Figure_08_WGCNA_module_trait_landscape",
  "Extended_Data_Figure_09_WGCNA_hubs_enrichment",
  "Extended_Data_Figure_10_WGCNA_context_stability",
  "Extended_Data_Figure_11_WGCNA_structural_preservation"
)

copied_extended <- unlist(
  purrr::map(
    extended_bases,
    function(base) {
      files <- c(
        file.path(S18_FINAL, "Illustrator_ready", "master_vectors", paste0(base, ".pdf")),
        file.path(S18_FINAL, "Illustrator_ready", "master_vectors", paste0(base, ".svg")),
        file.path(S18_FINAL, "preview", paste0(base, ".png")),
        file.path(S18_FINAL, "submission", paste0(base, ".tiff"))
      )
      vapply(files, copy_required, character(1), to_dir = FIG_EXT_DIR)
    }
  ),
  use.names = FALSE
)

###############################################################################
# 10. FIGURE LEGENDS, OPTIONAL TEXT AND SUPPORTING AUDITS
###############################################################################

figure3_legend_src <- file.path(
  S17_FINAL,
  "diagnostics",
  "figure_legends.txt"
)
extended_legends_csv_src <- file.path(
  S18_FINAL,
  "diagnostics",
  "WGCNA_Extended_Data_figure_legends.csv"
)
extended_legends_md_src <- file.path(
  S18_FINAL,
  "diagnostics",
  "WGCNA_Extended_Data_figure_legends.md"
)
module_key_src <- file.path(
  S18_FINAL,
  "diagnostics",
  "WGCNA_module_identifier_key.csv"
)

if (!nonempty_file(figure3_legend_src)) {
  stop("Final Figure 3 legend is missing:\n", figure3_legend_src, call. = FALSE)
}
if (!nonempty_file(extended_legends_csv_src)) {
  stop("Final Extended Data legends are missing:\n", extended_legends_csv_src, call. = FALSE)
}

figure3_legend <- paste(
  readLines(figure3_legend_src, warn = FALSE),
  collapse = " "
)

extended_legends <- readr::read_csv(
  extended_legends_csv_src,
  show_col_types = FALSE,
  guess_max = 10000,
  name_repair = "unique"
)

if (!"Figure_ID" %in% names(extended_legends)) {
  names(extended_legends)[1] <- "Figure_ID"
}
if (!"Legend_text" %in% names(extended_legends)) {
  legend_col <- c("Legend", "legend", "Text", "text")
  legend_col <- legend_col[legend_col %in% names(extended_legends)][1]
  if (length(legend_col) == 0 || is.na(legend_col)) {
    stop(
      "Could not identify the legend-text column in:\n",
      extended_legends_csv_src,
      call. = FALSE
    )
  }
  names(extended_legends)[names(extended_legends) == legend_col] <- "Legend_text"
}
if (!"Legend_title" %in% names(extended_legends)) {
  title_col <- c("Title", "title")
  title_col <- title_col[title_col %in% names(extended_legends)][1]
  if (length(title_col) > 0 && !is.na(title_col)) {
    names(extended_legends)[names(extended_legends) == title_col] <- "Legend_title"
  } else {
    extended_legends$Legend_title <- extended_legends$Figure_ID
  }
}

combined_legends <- dplyr::bind_rows(
  tibble::tibble(
    Figure_ID = "Figure 3",
    Legend_title = "Modular organization of the outcome-independent ReDLat plasma proteome",
    Legend_text = figure3_legend
  ),
  extended_legends %>%
    dplyr::select(
      Figure_ID,
      Legend_title,
      Legend_text,
      dplyr::everything()
    )
)

readr::write_csv(
  combined_legends,
  OUTPUTS$combined_legends_csv
)

combined_legend_md <- paste0(
  purrr::pmap_chr(
    combined_legends,
    function(Figure_ID, Legend_title, Legend_text, ...) {
      paste0(
        "## ", Figure_ID, ". ", Legend_title,
        "\n\n", Legend_text, "\n"
      )
    }
  ),
  collapse = "\n"
)
writeLines(
  combined_legend_md,
  OUTPUTS$combined_legends_md,
  useBytes = TRUE
)

if (nonempty_file(extended_legends_md_src)) {
  copy_required(extended_legends_md_src, LEGEND_DIR)
}
if (nonempty_file(module_key_src)) {
  copy_required(module_key_src, LEGEND_DIR)
}
copy_required(figure3_legend_src, LEGEND_DIR, "Figure3_Legend_Source.txt")

text_files <- if (dir.exists(S20_FINAL)) {
  list.files(
    S20_FINAL,
    pattern = "\\.(txt|md|csv)$",
    full.names = TRUE
  )
} else {
  character()
}

copied_text <- if (length(text_files) > 0) {
  vapply(
    text_files,
    copy_required,
    character(1),
    to_dir = TEXT_DIR
  )
} else {
  character()
}

audit_files <- list.files(
  S19_FINAL,
  pattern = "\\.(csv|txt)$",
  full.names = TRUE
)
audit_files <- unique(audit_files)

copied_audit <- if (length(audit_files) > 0) {
  vapply(
    audit_files,
    copy_required,
    character(1),
    to_dir = AUDIT_DIR
  )
} else {
  character()
}

###############################################################################
# 11. GLOBAL INDEX AND README
###############################################################################

package_files <- list.files(
  PACKAGE_ROOT,
  full.names = TRUE,
  recursive = TRUE,
  include.dirs = FALSE
)

index_tbl <- tibble::tibble(
  Artifact = basename(package_files),
  Relative_path = substring(
    normalizePath(package_files, winslash = "/", mustWork = TRUE),
    nchar(normalizePath(PACKAGE_ROOT, winslash = "/", mustWork = TRUE)) + 2L
  ),
  Extension = tolower(tools::file_ext(package_files)),
  Bytes = as.numeric(file.info(package_files)$size)
) %>%
  dplyr::mutate(
    Category = dplyr::case_when(
      stringr::str_starts(Relative_path, "Figures/Main") ~ "Main figure",
      stringr::str_starts(Relative_path, "Figures/Extended_Data") ~ "Extended Data figure",
      stringr::str_starts(Relative_path, "Tables") ~ "Supplementary tables",
      stringr::str_starts(Relative_path, "Source_Data") ~ "Source Data",
      stringr::str_starts(Relative_path, "Supplementary_Data") ~ "Supplementary Data",
      stringr::str_starts(Relative_path, "Manuscript_Text_Optional") ~ "Optional manuscript text",
      stringr::str_starts(Relative_path, "Figure_Legends") ~ "Figure legends",
      stringr::str_starts(Relative_path, "Audit") ~ "Audit",
      stringr::str_starts(Relative_path, "Legacy") ~ "Legacy registry",
      TRUE ~ "Other"
    )
  ) %>%
  dplyr::arrange(Category, Relative_path)

readr::write_csv(index_tbl, OUTPUTS$global_index_csv)

wb_index <- openxlsx::createWorkbook(creator = "Matías Pizarro")
write_sheet(
  wb_index,
  "Package_Index",
  "WGCNA Nature Aging final submission package index",
  "Complete inventory of the definitive WGCNA manuscript-delivery package.",
  index_tbl
)
openxlsx::saveWorkbook(wb_index, OUTPUTS$global_index_xlsx, overwrite = TRUE)

source_registry <- tibble::tibble(
  Source_ID = c(
    "Script17_final_figure",
    "Script18_final_extended_data",
    "Script19_final_tables",
    "Script20_optional_text",
    "Script14b_biomarker",
    "Script15b_preservation",
    "Script16_network_quality"
  ),
  Path = c(
    S17_FINAL,
    S18_FINAL,
    S19_FINAL,
    S20_FINAL,
    S14B,
    S15B,
    S16
  ),
  Exists = file.exists(
    c(
      S17_FINAL,
      S18_FINAL,
      S19_FINAL,
      S20_FINAL,
      S14B,
      S15B,
      S16
    )
  ),
  Required = c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE),
  Rule = c(
    "Definitive Figure 3 v16n",
    "Definitive Extended Data Figs. 7-11 v6",
    "Definitive Supplementary Tables 20-27",
    "Optional manuscript-facing Methods, Results and interpretation text",
    "Only corrected family of 32 biomarker models",
    "Only fixed 4,239-gene structural preservation",
    "Only exact post hoc weighted modularity and current diagnostics"
  )
)
readr::write_csv(source_registry, OUTPUTS$source_registry)

readme <- c(
  "WGCNA NATURE AGING FINAL PACKAGE",
  "",
  "Contents:",
  "- Figure 3 in PDF, SVG, PNG and TIFF.",
  "- Extended Data Figs. 7-11 in PDF, SVG, PNG and TIFF.",
  "- Supplementary Tables 20-27.",
  "- Source Data for Figure 3 and Extended Data Figs. 7-11.",
  "- Supplementary Data 4-6.",
  "- Figure legends generated directly from the definitive figure workflows.",
  "- Optional integrated manuscript Methods and Results when Script 20 succeeds.",
  "- Source registry, pre-package audit and package manifest.",
  "",
  "Authoritative rules:",
  "1. M1-M8 are stable visual identifiers and do not represent a ranking.",
  "2. Biomarker inference uses Script 14b and BH-FDR across 32 models.",
  "3. Structural preservation uses Script 15b and one fixed 4,239-gene set.",
  "4. Exact weighted modularity comes from Script 16 and is descriptive.",
  "5. Association stability and structural preservation are distinct.",
  "6. The package provides internal robustness evidence, not external validation.",
  "7. Script 22 is the post-package audit and is intentionally run after Script 21.",
  "",
  "Do not submit or cite legacy Script 14 family-12 biomarker FDR results,",
  "legacy Script 15 context-specific preservation draws, or previous WGCNA",
  "figure-development versions."
)
writeLines(readme, OUTPUTS$readme, useBytes = TRUE)

writeLines(
  capture.output(utils::sessionInfo()),
  OUTPUTS$session_info
)

###############################################################################
# 12. PACKAGE MANIFEST AND ZIP
###############################################################################

# Regenerate the package index after README, source registry and session files
# have been added, so the index reflects the complete final directory.
package_files_for_index <- list.files(
  PACKAGE_ROOT,
  full.names = TRUE,
  recursive = TRUE,
  include.dirs = FALSE
)

index_tbl <- tibble::tibble(
  Artifact = basename(package_files_for_index),
  Relative_path = substring(
    normalizePath(package_files_for_index, winslash = "/", mustWork = TRUE),
    nchar(normalizePath(PACKAGE_ROOT, winslash = "/", mustWork = TRUE)) + 2L
  ),
  Extension = tolower(tools::file_ext(package_files_for_index)),
  Bytes = as.numeric(file.info(package_files_for_index)$size)
) %>%
  dplyr::mutate(
    Category = dplyr::case_when(
      stringr::str_starts(Relative_path, "Figures/Main") ~ "Main figure",
      stringr::str_starts(Relative_path, "Figures/Extended_Data") ~ "Extended Data figure",
      stringr::str_starts(Relative_path, "Tables") ~ "Supplementary tables",
      stringr::str_starts(Relative_path, "Source_Data") ~ "Source Data",
      stringr::str_starts(Relative_path, "Supplementary_Data") ~ "Supplementary Data",
      stringr::str_starts(Relative_path, "Manuscript_Text_Optional") ~ "Optional manuscript text",
      stringr::str_starts(Relative_path, "Figure_Legends") ~ "Figure legends",
      stringr::str_starts(Relative_path, "Audit") ~ "Audit",
      stringr::str_starts(Relative_path, "Legacy") ~ "Legacy registry",
      TRUE ~ "Other"
    )
  ) %>%
  dplyr::arrange(Category, Relative_path)

readr::write_csv(index_tbl, OUTPUTS$global_index_csv)
wb_index <- openxlsx::createWorkbook(creator = "Matías Pizarro")
write_sheet(
  wb_index,
  "Package_Index",
  "WGCNA Nature Aging final submission package index",
  "Complete inventory of the definitive WGCNA manuscript-delivery package.",
  index_tbl
)
openxlsx::saveWorkbook(wb_index, OUTPUTS$global_index_xlsx, overwrite = TRUE)

required_package_outputs <- c(
  OUTPUTS$supplementary_tables,
  OUTPUTS$source_fig3,
  OUTPUTS$source_ed7,
  OUTPUTS$source_ed8,
  OUTPUTS$source_ed9,
  OUTPUTS$source_ed10,
  OUTPUTS$source_ed11,
  OUTPUTS$supp_data4,
  OUTPUTS$supp_data5,
  OUTPUTS$supp_data6,
  OUTPUTS$combined_legends_csv,
  OUTPUTS$combined_legends_md,
  OUTPUTS$global_index_csv,
  OUTPUTS$global_index_xlsx,
  OUTPUTS$source_registry,
  OUTPUTS$readme
)

required_figure_outputs <- c(
  copied_main,
  copied_extended
)

preaudit <- tibble::tibble(
  Artifact = c(
    basename(required_package_outputs),
    basename(required_figure_outputs)
  ),
  Path = c(
    required_package_outputs,
    required_figure_outputs
  )
) %>%
  dplyr::mutate(
    Exists = file.exists(Path),
    Bytes = ifelse(Exists, file.info(Path)$size, NA_real_),
    Pass = Exists & is.finite(Bytes) & Bytes > 0
  )

readr::write_csv(preaudit, OUTPUTS$preaudit)

if (!all(preaudit$Pass)) {
  failed <- preaudit %>% dplyr::filter(!Pass)
  stop(
    "The package pre-audit failed. Missing or empty artifacts:\n",
    paste(failed$Path, collapse = "\n"),
    call. = FALSE
  )
}

package_files <- list.files(
  PACKAGE_ROOT,
  full.names = TRUE,
  recursive = TRUE,
  include.dirs = FALSE
)

manifest <- tibble::tibble(
  Path = normalizePath(package_files, winslash = "/", mustWork = TRUE),
  Relative_path = substring(
    normalizePath(package_files, winslash = "/", mustWork = TRUE),
    nchar(normalizePath(PACKAGE_ROOT, winslash = "/", mustWork = TRUE)) + 2L
  ),
  Bytes = as.numeric(file.info(package_files)$size),
  Exists = file.exists(package_files)
)

readr::write_csv(manifest, OUTPUTS$package_manifest)

if (!all(manifest$Exists) || any(manifest$Bytes <= 0, na.rm = TRUE)) {
  stop("One or more final package files are missing or empty.", call. = FALSE)
}

if (file.exists(OUTPUTS$zip)) file.remove(OUTPUTS$zip)

wgcna_assert_public_tree(PACKAGE_ROOT)

zip::zipr(
  zipfile = OUTPUTS$zip,
  files = basename(PACKAGE_ROOT),
  root = dirname(PACKAGE_ROOT),
  recurse = TRUE,
  include_directories = TRUE
)

if (!file.exists(OUTPUTS$zip) || file.info(OUTPUTS$zip)$size <= 0) {
  stop("Final WGCNA ZIP package was not created correctly.", call. = FALSE)
}

final_summary <- tibble::tibble(
  Metric = c(
    "Main figures", "Extended Data figures", "Supplementary Tables",
    "Source Data workbooks", "Supplementary Data workbooks",
    "Biomarker FDR family", "Fixed preservation genes",
    "Final package files", "Final ZIP"
  ),
  Value = c(
    "1", "5", "20-27", "6", "3", "32", "4239",
    as.character(nrow(manifest)), OUTPUTS$zip
  )
)
readr::write_csv(
  final_summary,
  file.path(OUT_ROOT, "script22_final_summary.csv")
)

writeLines(
  c(
    "SCRIPT_21_STATUS=PASS",
    paste0("PACKAGE_ROOT=", PACKAGE_ROOT),
    paste0("FINAL_ZIP=", OUTPUTS$zip),
    paste0("TIMESTAMP=", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  ),
  file.path(OUT_ROOT, "SCRIPT_21_SUCCESS.txt")
)

message("Script 21 completed successfully.")
message("Final package directory: ", PACKAGE_ROOT)
message("Final ZIP: ", OUTPUTS$zip)
message("Supplementary Tables: 20-27")
message("Source Data workbooks: Figure 3 and Extended Data Figs. 7-11")
message("Supplementary Data workbooks: 4-6")

} # end run_submission_package()

message("============================================================")
message("RUNNING SCRIPT 14: SUBMISSION PACKAGE")
message("The package is created before the final audit.")
message("Script 20 text is included when available and cannot invalidate numerical outputs.")
message("============================================================")

run_submission_package()

###############################################################################
# END
###############################################################################


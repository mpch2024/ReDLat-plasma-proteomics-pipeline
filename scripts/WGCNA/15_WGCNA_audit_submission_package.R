###############################################################################
# ReDLat plasma proteomics — WGCNA workflow
# 15. Audit the submission package
# Requires: output from Script 14
# Produces: scientific, workbook, privacy and completeness audits
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

while (grDevices::dev.cur() > 1) grDevices::dev.off()
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

run_package_audit <- function() {

###############################################################################
# 00. PATHS
###############################################################################

BASE_DIR <- WGCNA_CONFIG$project_root

SCRIPT_DIR <- file.path(WGCNA_CONFIG$project_root, "scripts", "WGCNA")
WGCNA_ROOT <- WGCNA_CONFIG$result_root

SCRIPT21_FINAL <- file.path(
  SCRIPT_DIR,
  "14_WGCNA_build_submission_package.R"
)

S14B <- file.path(
  WGCNA_ROOT,
  "07_biomarker_fdr"
)
S15B <- file.path(
  WGCNA_ROOT,
  "08_preservation"
)
S16 <- file.path(
  WGCNA_ROOT,
  "09_network_quality"
)
S19 <- file.path(WGCNA_CONFIG$publication_root, "supplementary_tables")

S22 <- file.path(WGCNA_CONFIG$publication_root, "submission_package")
PACKAGE_ROOT <- file.path(WGCNA_CONFIG$publication_root, "submission_package", "WGCNA_NatureAging_PACKAGE")
FINAL_ZIP <- file.path(
  S22,
  "WGCNA_NatureAging_PACKAGE.zip"
)

OUTDIR <- file.path(WGCNA_CONFIG$publication_root, "submission_audit")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

PACKAGE_AUDIT_DIR <- file.path(PACKAGE_ROOT, "Audit")
dir.create(PACKAGE_AUDIT_DIR, recursive = TRUE, showWarnings = FALSE)

###############################################################################
# 01. HELPERS
###############################################################################

nonempty_file <- function(path) {
  file.exists(path) &&
    is.finite(file.info(path)$size) &&
    file.info(path)$size > 0
}

run_isolated_script <- function(script_path, label) {
  if (!file.exists(script_path)) {
    stop("Missing required script for ", label, ":\n", script_path, call. = FALSE)
  }

  script_env <- new.env(parent = globalenv())

  tryCatch(
    {
      sys.source(
        script_path,
        envir = script_env,
        chdir = FALSE,
        keep.source = TRUE
      )
    },
    error = function(e) {
      stop(label, " failed:\n", conditionMessage(e), call. = FALSE)
    }
  )

  invisible(TRUE)
}

file_size_safe <- function(path) {
  if (!file.exists(path)) return(NA_real_)
  as.numeric(file.info(path)$size)
}

read_raw_prefix <- function(path, n = 12L) {
  if (!file.exists(path)) return(raw())
  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  readBin(con, what = "raw", n = n)
}

file_signature_ok <- function(path) {
  if (!nonempty_file(path)) return(FALSE)

  ext <- tolower(tools::file_ext(path))
  prefix <- read_raw_prefix(path, 12L)

  if (ext == "pdf") {
    return(
      length(prefix) >= 4L &&
        rawToChar(prefix[1:4]) == "%PDF"
    )
  }

  if (ext == "png") {
    return(
      length(prefix) >= 4L &&
        identical(
          prefix[1:4],
          as.raw(c(0x89, 0x50, 0x4E, 0x47))
        )
    )
  }

  if (ext %in% c("tif", "tiff")) {
    little <- as.raw(c(0x49, 0x49, 0x2A, 0x00))
    big <- as.raw(c(0x4D, 0x4D, 0x00, 0x2A))
    return(
      length(prefix) >= 4L &&
        (
          identical(prefix[1:4], little) ||
            identical(prefix[1:4], big)
        )
    )
  }

  if (ext == "svg") {
    txt <- paste(readLines(path, n = 20L, warn = FALSE), collapse = " ")
    return(grepl("<svg", txt, fixed = TRUE))
  }

  if (ext %in% c("xlsx", "zip")) {
    return(
      length(prefix) >= 2L &&
        identical(prefix[1:2], as.raw(c(0x50, 0x4B)))
    )
  }

  TRUE
}

relative_to_package <- function(path) {
  root <- normalizePath(PACKAGE_ROOT, winslash = "/", mustWork = FALSE)
  p <- normalizePath(path, winslash = "/", mustWork = FALSE)
  prefix <- paste0(root, "/")

  if (startsWith(p, prefix)) {
    substring(p, nchar(prefix) + 1L)
  } else {
    p
  }
}

read_csv_required <- function(path, label) {
  if (!nonempty_file(path)) {
    stop("Missing required ", label, ":\n", path, call. = FALSE)
  }

  readr::read_csv(
    path,
    show_col_types = FALSE,
    guess_max = 100000,
    name_repair = "unique"
  )
}

###############################################################################
# 02. SELF-HEALING PACKAGE BUILD
###############################################################################

core_package_probe <- c(
  file.path(
    PACKAGE_ROOT,
    "Tables",
    "Supplementary_Tables_WGCNA_NatureAging_CLEAN.xlsx"
  ),
  file.path(
    PACKAGE_ROOT,
    "Source_Data",
    "SourceData_Fig3_WGCNA_NatureAging.xlsx"
  ),
  file.path(
    PACKAGE_ROOT,
    "Supplementary_Data",
    "Supplementary_Data_6_WGCNA_Complete_Preservation_and_Diagnostics.xlsx"
  )
)

if (!all(vapply(core_package_probe, nonempty_file, logical(1)))) {
  message("Final package is incomplete. Running Script 21 automatically.")
  run_isolated_script(
    SCRIPT21_FINAL,
    "Script 21 final package construction"
  )
}

if (!all(vapply(core_package_probe, nonempty_file, logical(1)))) {
  stop(
    "Script 21 completed but the final package remains incomplete:\n",
    paste(
      core_package_probe[
        !vapply(core_package_probe, nonempty_file, logical(1))
      ],
      collapse = "\n"
    ),
    call. = FALSE
  )
}

###############################################################################
# 03. EXPECTED FINAL ARTIFACTS
###############################################################################

main_base <- "Figure3_WGCNA_COMPACT_CLASSIC_v16n_FINAL"
main_expected <- file.path(
  PACKAGE_ROOT,
  "Figures",
  "Main",
  paste0(
    main_base,
    c(".pdf", ".svg", ".png", ".tiff")
  )
)

extended_bases <- c(
  "Extended_Data_Figure_07_WGCNA_network_quality",
  "Extended_Data_Figure_08_WGCNA_module_trait_landscape",
  "Extended_Data_Figure_09_WGCNA_hubs_enrichment",
  "Extended_Data_Figure_10_WGCNA_context_stability",
  "Extended_Data_Figure_11_WGCNA_structural_preservation"
)

extended_expected <- unlist(
  lapply(
    extended_bases,
    function(stem) {
      file.path(
        PACKAGE_ROOT,
        "Figures",
        "Extended_Data",
        paste0(
          stem,
          c(".pdf", ".svg", ".png", ".tiff")
        )
      )
    }
  ),
  use.names = FALSE
)

table_expected <- file.path(
  PACKAGE_ROOT,
  "Tables",
  "Supplementary_Tables_WGCNA_NatureAging_CLEAN.xlsx"
)

source_data_expected <- file.path(
  PACKAGE_ROOT,
  "Source_Data",
  c(
    "SourceData_Fig3_WGCNA_NatureAging.xlsx",
    "SourceData_ExtendedData_Figure_07_WGCNA.xlsx",
    "SourceData_ExtendedData_Figure_08_WGCNA.xlsx",
    "SourceData_ExtendedData_Figure_09_WGCNA.xlsx",
    "SourceData_ExtendedData_Figure_10_WGCNA.xlsx",
    "SourceData_ExtendedData_Figure_11_WGCNA.xlsx"
  )
)

supp_data_expected <- file.path(
  PACKAGE_ROOT,
  "Supplementary_Data",
  c(
    "Supplementary_Data_4_WGCNA_Complete_Module_Results.xlsx",
    "Supplementary_Data_5_WGCNA_Complete_Association_and_Context_Results.xlsx",
    "Supplementary_Data_6_WGCNA_Complete_Preservation_and_Diagnostics.xlsx"
  )
)

legend_expected <- file.path(
  PACKAGE_ROOT,
  "Figure_Legends",
  c(
    "WGCNA_Final_Figure_Legends.csv",
    "WGCNA_Final_Figure_Legends.md"
  )
)

support_expected <- file.path(
  PACKAGE_ROOT,
  c(
    "WGCNA_Submission_Package_Index.csv",
    "WGCNA_Submission_Package_Index.xlsx",
    "README_WGCNA_NatureAging_FINAL.txt"
  )
)

registry_expected <- file.path(
  PACKAGE_ROOT,
  "Audit",
  c(
    "WGCNA_Source_Registry_FINAL.csv",
    "WGCNA_PrePackage_Audit.csv"
  )
)

required_paths <- c(
  main_expected,
  extended_expected,
  table_expected,
  source_data_expected,
  supp_data_expected,
  legend_expected,
  support_expected,
  registry_expected
)

required_inventory <- tibble::tibble(
  Artifact = basename(required_paths),
  Relative_path = vapply(
    required_paths,
    relative_to_package,
    character(1)
  ),
  Path = required_paths
) %>%
  dplyr::mutate(
    Exists = file.exists(Path),
    Bytes = ifelse(Exists, file.info(Path)$size, NA_real_),
    Signature_OK = vapply(Path, file_signature_ok, logical(1)),
    Pass = Exists & is.finite(Bytes) & Bytes > 0 & Signature_OK
  )

###############################################################################
# 04. SCIENTIFIC-SOURCE RULE AUDIT
###############################################################################

biomarker_path <- file.path(
  S14B,
  "tables",
  "full",
  "full_all_modules_log_HC3_family32.csv"
)
preservation_path <- file.path(
  S15B,
  "tables",
  "fixed_geneset_structural_preservation_definitive_summary.csv"
)
fixed_manifest_path <- file.path(
  S15B,
  "tables",
  "input",
  "fixed_gene_manifest_all_modules.csv"
)
weighted_q_path <- file.path(
  S16,
  "tables",
  "modularity",
  "posthoc_weighted_modularity_Q.csv"
)
script19_audit_path <- file.path(
  S19,
  "WGCNA_Current_Analysis_Audit.csv"
)

biomarker_tbl <- read_csv_required(
  biomarker_path,
  "corrected biomarker family"
)
preservation_tbl <- read_csv_required(
  preservation_path,
  "fixed-gene preservation summary"
)
fixed_manifest_tbl <- read_csv_required(
  fixed_manifest_path,
  "fixed-gene manifest"
)
weighted_q_tbl <- read_csv_required(
  weighted_q_path,
  "weighted modularity Q"
)
script19_audit <- read_csv_required(
  script19_audit_path,
  "Script 19 current-analysis audit"
)

family_col <- c("FDR_family_size", "family_size")
family_col <- family_col[family_col %in% names(biomarker_tbl)][1]

fdr_col <- c("FDR_family32", "full_FDR_family32")
fdr_col <- fdr_col[fdr_col %in% names(biomarker_tbl)][1]

comparability_col <- c(
  "Fixed_gene_comparability",
  "fixed_gene_comparability"
)
comparability_col <- comparability_col[
  comparability_col %in% names(preservation_tbl)
][1]

scientific_audit <- tibble::tribble(
  ~Check, ~Expected, ~Observed, ~Pass,
  "Corrected biomarker model rows",
  "32",
  as.character(nrow(biomarker_tbl)),
  nrow(biomarker_tbl) == 32L,

  "Corrected biomarker FDR field",
  "FDR_family32",
  ifelse(length(fdr_col) == 0 || is.na(fdr_col), "missing", fdr_col),
  length(fdr_col) > 0 && !is.na(fdr_col),

  "Biomarker family size",
  "32 in every row",
  ifelse(
    length(family_col) == 0 || is.na(family_col),
    "missing",
    paste(unique(biomarker_tbl[[family_col]]), collapse = "; ")
  ),
  length(family_col) > 0 &&
    !is.na(family_col) &&
    all(biomarker_tbl[[family_col]] == 32, na.rm = TRUE),

  "Fixed preservation gene set",
  "4239",
  as.character(nrow(fixed_manifest_tbl)),
  nrow(fixed_manifest_tbl) == 4239L,

  "Fixed-gene comparability field",
  "present",
  ifelse(
    length(comparability_col) == 0 || is.na(comparability_col),
    "missing",
    comparability_col
  ),
  length(comparability_col) > 0 && !is.na(comparability_col),

  "Exact weighted modularity source",
  "non-empty",
  as.character(nrow(weighted_q_tbl)),
  nrow(weighted_q_tbl) >= 1L,

  "Script 19 current-analysis audit",
  "all PASS",
  paste(unique(script19_audit$Pass), collapse = "; "),
  "Pass" %in% names(script19_audit) &&
    all(script19_audit$Pass)
)

###############################################################################
# 05. WORKBOOK CONTENT AUDIT
###############################################################################

supp_table_sheets <- tryCatch(
  openxlsx::getSheetNames(table_expected),
  error = function(e) character()
)

required_supp_sheets <- c(
  "Index",
  "Module_Key",
  "Current_Analysis_Audit",
  paste0("S", 20:27)
)

workbook_audit <- tibble::tribble(
  ~Check, ~Expected, ~Observed, ~Pass,
  "Supplementary Tables sheets",
  paste(required_supp_sheets, collapse = "; "),
  paste(supp_table_sheets, collapse = "; "),
  all(required_supp_sheets %in% supp_table_sheets),

  "Source Data workbook count",
  "6",
  as.character(sum(vapply(source_data_expected, nonempty_file, logical(1)))),
  sum(vapply(source_data_expected, nonempty_file, logical(1))) == 6L,

  "Supplementary Data workbook count",
  "3",
  as.character(sum(vapply(supp_data_expected, nonempty_file, logical(1)))),
  sum(vapply(supp_data_expected, nonempty_file, logical(1))) == 3L,

  "Main figure format count",
  "4",
  as.character(sum(vapply(main_expected, nonempty_file, logical(1)))),
  sum(vapply(main_expected, nonempty_file, logical(1))) == 4L,

  "Extended Data figure format count",
  "20",
  as.character(sum(vapply(extended_expected, nonempty_file, logical(1)))),
  sum(vapply(extended_expected, nonempty_file, logical(1))) == 20L
)

###############################################################################
# 06. OPTIONAL EDITORIAL TEXT STATUS
###############################################################################

optional_text_dir <- file.path(
  PACKAGE_ROOT,
  "Manuscript_Text_Optional"
)

optional_text_files <- if (dir.exists(optional_text_dir)) {
  list.files(
    optional_text_dir,
    pattern = "\\.(txt|md|csv)$",
    full.names = TRUE
  )
} else {
  character()
}

optional_text_status <- tibble::tibble(
  Component = "Optional manuscript-facing text",
  Required = FALSE,
  Files_present = length(optional_text_files),
  Status = ifelse(
    length(optional_text_files) > 0,
    "AVAILABLE",
    "NOT_AVAILABLE_BUT_NOT_REQUIRED"
  )
)

###############################################################################
# 07. FINAL CHECKLIST
###############################################################################

checklist <- tibble::tribble(
  ~Item, ~Required, ~Status, ~Evidence,
  "All required package artifacts exist and have valid signatures",
  TRUE,
  all(required_inventory$Pass),
  "WGCNA_Final_Required_Artifact_Audit.csv",

  "Corrected biomarker family and fixed-gene preservation rules",
  TRUE,
  all(scientific_audit$Pass),
  "WGCNA_Final_Scientific_Source_Audit.csv",

  "Supplementary Tables and workbook counts",
  TRUE,
  all(workbook_audit$Pass),
  "WGCNA_Final_Workbook_Audit.csv",

  "Optional manuscript text",
  FALSE,
  TRUE,
  "WGCNA_Optional_Text_Status.csv"
)

failed_checks <- checklist %>%
  dplyr::filter(Required, !Status)

failed_artifacts <- required_inventory %>%
  dplyr::filter(!Pass)

failed_science <- scientific_audit %>%
  dplyr::filter(!Pass)

failed_workbooks <- workbook_audit %>%
  dplyr::filter(!Pass)

###############################################################################
# 08. WRITE AUDIT INTO PACKAGE
###############################################################################

audit_outputs <- list(
  required = file.path(
    PACKAGE_AUDIT_DIR,
    "WGCNA_Final_Required_Artifact_Audit.csv"
  ),
  scientific = file.path(
    PACKAGE_AUDIT_DIR,
    "WGCNA_Final_Scientific_Source_Audit.csv"
  ),
  workbook = file.path(
    PACKAGE_AUDIT_DIR,
    "WGCNA_Final_Workbook_Audit.csv"
  ),
  optional_text = file.path(
    PACKAGE_AUDIT_DIR,
    "WGCNA_Optional_Text_Status.csv"
  ),
  checklist = file.path(
    PACKAGE_AUDIT_DIR,
    "WGCNA_FINAL_SUBMISSION_CHECKLIST.csv"
  ),
  failed_checks = file.path(
    PACKAGE_AUDIT_DIR,
    "WGCNA_FAILED_CHECKS.csv"
  ),
  failed_artifacts = file.path(
    PACKAGE_AUDIT_DIR,
    "WGCNA_FAILED_REQUIRED_ARTIFACTS.csv"
  ),
  failed_science = file.path(
    PACKAGE_AUDIT_DIR,
    "WGCNA_FAILED_SCIENTIFIC_RULES.csv"
  ),
  failed_workbooks = file.path(
    PACKAGE_AUDIT_DIR,
    "WGCNA_FAILED_WORKBOOK_CHECKS.csv"
  )
)

readr::write_csv(required_inventory, audit_outputs$required)
readr::write_csv(scientific_audit, audit_outputs$scientific)
readr::write_csv(workbook_audit, audit_outputs$workbook)
readr::write_csv(optional_text_status, audit_outputs$optional_text)
readr::write_csv(checklist, audit_outputs$checklist)
readr::write_csv(failed_checks, audit_outputs$failed_checks)
readr::write_csv(failed_artifacts, audit_outputs$failed_artifacts)
readr::write_csv(failed_science, audit_outputs$failed_science)
readr::write_csv(failed_workbooks, audit_outputs$failed_workbooks)

# Mirror the key audits outside the package for rapid inspection.
readr::write_csv(
  checklist,
  file.path(OUTDIR, "WGCNA_FINAL_SUBMISSION_CHECKLIST.csv")
)
readr::write_csv(
  required_inventory,
  file.path(OUTDIR, "WGCNA_Final_Required_Artifact_Audit.csv")
)
readr::write_csv(
  scientific_audit,
  file.path(OUTDIR, "WGCNA_Final_Scientific_Source_Audit.csv")
)
readr::write_csv(
  workbook_audit,
  file.path(OUTDIR, "WGCNA_Final_Workbook_Audit.csv")
)

if (nrow(failed_checks) > 0) {
  message("")
  message("============================================================")
  message("FINAL WGCNA PACKAGE AUDIT FAILED")
  message("============================================================")

  for (i in seq_len(nrow(failed_checks))) {
    message(
      "- ", failed_checks$Item[i],
      " | Evidence: ", failed_checks$Evidence[i]
    )
  }

  if (nrow(failed_artifacts) > 0) {
    message("")
    message("Missing or invalid artifacts:")
    for (i in seq_len(nrow(failed_artifacts))) {
      message(
        "- ", failed_artifacts$Relative_path[i],
        " | Exists=", failed_artifacts$Exists[i],
        " | Bytes=", failed_artifacts$Bytes[i],
        " | Signature_OK=", failed_artifacts$Signature_OK[i]
      )
    }
  }

  if (nrow(failed_science) > 0) {
    message("")
    message("Scientific-rule failures:")
    for (i in seq_len(nrow(failed_science))) {
      message(
        "- ", failed_science$Check[i],
        " | Expected=", failed_science$Expected[i],
        " | Observed=", failed_science$Observed[i]
      )
    }
  }

  if (nrow(failed_workbooks) > 0) {
    message("")
    message("Workbook failures:")
    for (i in seq_len(nrow(failed_workbooks))) {
      message(
        "- ", failed_workbooks$Check[i],
        " | Expected=", failed_workbooks$Expected[i],
        " | Observed=", failed_workbooks$Observed[i]
      )
    }
  }

  stop(
    "The final package audit failed. Exact failures were written to the package Audit directory.",
    call. = FALSE
  )
}

###############################################################################
# 09. REFRESH FINAL INDEX, MANIFEST AND ZIP AFTER AUDIT
###############################################################################

package_files <- list.files(
  PACKAGE_ROOT,
  full.names = TRUE,
  recursive = TRUE,
  include.dirs = FALSE
)

package_root_norm <- normalizePath(
  PACKAGE_ROOT,
  winslash = "/",
  mustWork = TRUE
)

index_tbl <- tibble::tibble(
  Artifact = basename(package_files),
  Relative_path = substring(
    normalizePath(package_files, winslash = "/", mustWork = TRUE),
    nchar(package_root_norm) + 2L
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
      stringr::str_starts(Relative_path, "Figure_Legends") ~ "Figure legends",
      stringr::str_starts(Relative_path, "Manuscript_Text_Optional") ~ "Optional manuscript text",
      stringr::str_starts(Relative_path, "Audit") ~ "Audit",
      TRUE ~ "Other"
    )
  ) %>%
  dplyr::arrange(Category, Relative_path)

index_csv <- file.path(
  PACKAGE_ROOT,
  "WGCNA_Submission_Package_Index.csv"
)
index_xlsx <- file.path(
  PACKAGE_ROOT,
  "WGCNA_Submission_Package_Index.xlsx"
)
manifest_path <- file.path(
  PACKAGE_AUDIT_DIR,
  "WGCNA_Package_Manifest.csv"
)

readr::write_csv(index_tbl, index_csv)

wb <- openxlsx::createWorkbook(creator = "Matías Pizarro")
openxlsx::addWorksheet(wb, "Package_Index", gridLines = FALSE)
openxlsx::writeData(wb, "Package_Index", index_tbl, withFilter = TRUE)
openxlsx::freezePane(wb, "Package_Index", firstRow = TRUE)
openxlsx::setColWidths(
  wb,
  "Package_Index",
  cols = seq_len(ncol(index_tbl)),
  widths = c(34, 75, 12, 16, 24)
)
openxlsx::saveWorkbook(wb, index_xlsx, overwrite = TRUE)

# Two-pass manifest so the manifest includes the refreshed index and audits.
package_files <- list.files(
  PACKAGE_ROOT,
  full.names = TRUE,
  recursive = TRUE,
  include.dirs = FALSE
)

manifest <- tibble::tibble(
  Path = normalizePath(
    package_files,
    winslash = "/",
    mustWork = TRUE
  ),
  Relative_path = substring(
    normalizePath(package_files, winslash = "/", mustWork = TRUE),
    nchar(package_root_norm) + 2L
  ),
  Bytes = as.numeric(file.info(package_files)$size),
  Signature_OK = vapply(package_files, file_signature_ok, logical(1))
)

readr::write_csv(manifest, manifest_path)

if (file.exists(FINAL_ZIP)) file.remove(FINAL_ZIP)

wgcna_assert_public_tree(PACKAGE_ROOT)

zip::zipr(
  zipfile = FINAL_ZIP,
  files = basename(PACKAGE_ROOT),
  root = dirname(PACKAGE_ROOT),
  recurse = TRUE,
  include_directories = TRUE
)

zip_pass <- nonempty_file(FINAL_ZIP) && file_signature_ok(FINAL_ZIP)

if (!zip_pass) {
  stop("The final audited ZIP could not be created correctly.", call. = FALSE)
}

final_summary <- tibble::tibble(
  Metric = c(
    "Required package checks",
    "Scientific-source checks",
    "Workbook checks",
    "Optional manuscript text files",
    "Final package files",
    "Final ZIP",
    "Audit passed"
  ),
  Value = c(
    paste0(sum(required_inventory$Pass), "/", nrow(required_inventory)),
    paste0(sum(scientific_audit$Pass), "/", nrow(scientific_audit)),
    paste0(sum(workbook_audit$Pass), "/", nrow(workbook_audit)),
    as.character(length(optional_text_files)),
    as.character(nrow(manifest)),
    FINAL_ZIP,
    "TRUE"
  )
)

readr::write_csv(
  final_summary,
  file.path(OUTDIR, "script20_final_summary.csv")
)
readr::write_csv(
  final_summary,
  file.path(PACKAGE_AUDIT_DIR, "script20_final_summary.csv")
)

writeLines(
  c(
    "SCRIPT_22_STATUS=PASS",
    paste0("PACKAGE_ROOT=", PACKAGE_ROOT),
    paste0("FINAL_ZIP=", FINAL_ZIP),
    paste0("TIMESTAMP=", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  ),
  file.path(OUTDIR, "SCRIPT_22_SUCCESS.txt")
)

writeLines(
  capture.output(utils::sessionInfo()),
  file.path(OUTDIR, "sessionInfo.txt")
)

message("============================================================")
message("FINAL WGCNA PACKAGE AUDIT PASSED")
message("============================================================")
message("Final audited ZIP:")
message(FINAL_ZIP)
message("Optional Script 21 text files present: ", length(optional_text_files))

} # end run_package_audit()

message("============================================================")
message("RUNNING SCRIPT 15: SUBMISSION PACKAGE AUDIT")
message("This is the final post-package audit.")
message("Required numerical artifacts are audited separately from editorial text.")
message("============================================================")

run_package_audit()

###############################################################################
# END
###############################################################################


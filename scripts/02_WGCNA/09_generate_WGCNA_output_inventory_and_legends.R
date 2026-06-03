###############################################################################
# 08_wgcna_output_inventory_and_legends_FINAL_v2.R
#
# TITLE:
# Output inventory and figure legends for the WGCNA figure set
# (ReDLat proteomic project)
#
# PURPOSE:
# This script creates a manuscript-facing inventory of the final WGCNA figures,
# supplementary figures, source tables, and ready-to-edit figure legends.
#
# INPUTS:
# - Script 06 outputs:
#   results/06_wgcna_figure3_main/
# - Script 07 outputs:
#   results/07_wgcna_supplementary_figures/
#
# MAIN OUTPUTS:
# - wgcna_final_output_inventory.csv
# - wgcna_final_figure_legends.csv
# - wgcna_final_figure_legends.md
# - wgcna_submission_checklist.csv
#
# AUTHOR:
# Matías Pizarro + ChatGPT support
###############################################################################

rm(list = ls())

###############################################################################
# 1) PACKAGES
###############################################################################

cran_pkgs <- c(
  "dplyr",
  "readr",
  "tibble",
  "stringr",
  "purrr"
)

cran_missing <- cran_pkgs[!sapply(cran_pkgs, requireNamespace, quietly = TRUE)]
if (length(cran_missing) > 0) install.packages(cran_missing)

invisible(lapply(cran_pkgs, library, character.only = TRUE))

options(stringsAsFactors = FALSE)
options(error = traceback)

###############################################################################
# 2) PATHS
###############################################################################

BASE_DIR <- "C:/Users/mnpiz/Desktop/WGCNA_Workflow_april_V2"

SCRIPT06_DIR <- file.path(BASE_DIR, "results", "06_wgcna_figure3_main")
SCRIPT07_DIR <- file.path(BASE_DIR, "results", "07_wgcna_supplementary_figures")

SCRIPT06_FIG_DIR <- file.path(SCRIPT06_DIR, "figures")
SCRIPT06_TAB_DIR <- file.path(SCRIPT06_DIR, "tables")
SCRIPT06_PANEL_DIR <- file.path(SCRIPT06_DIR, "separate_panels")

SCRIPT07_FIG_DIR <- file.path(SCRIPT07_DIR, "figures")
SCRIPT07_TAB_DIR <- file.path(SCRIPT07_DIR, "tables")
SCRIPT07_PANEL_DIR <- file.path(SCRIPT07_DIR, "separate_panels")

OUTDIR <- file.path(BASE_DIR, "results", "08_wgcna_output_inventory_and_legends")
OUT_TAB <- file.path(OUTDIR, "tables")
OUT_TXT <- file.path(OUTDIR, "legends")

invisible(lapply(
  c(OUTDIR, OUT_TAB, OUT_TXT),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

###############################################################################
# 3) HELPERS
###############################################################################

file_exists_safe <- function(path) {
  ifelse(file.exists(path), TRUE, FALSE)
}

make_file_record <- function(figure_id,
                             figure_type,
                             description,
                             file_base,
                             directory,
                             source_table = NA_character_,
                             required = TRUE) {
  pdf_file <- file.path(directory, paste0(file_base, ".pdf"))
  png_file <- file.path(directory, paste0(file_base, ".png"))

  tibble::tibble(
    Figure_ID = figure_id,
    Figure_type = figure_type,
    Description = description,
    File_base = file_base,
    PDF_file = pdf_file,
    PNG_file = png_file,
    PDF_exists = file.exists(pdf_file),
    PNG_exists = file.exists(png_file),
    Source_table = source_table,
    Source_table_exists = ifelse(!is.na(source_table), file.exists(source_table), NA),
    Required = required
  )
}

write_markdown_legends <- function(legend_tbl, file) {
  lines <- purrr::pmap_chr(
    legend_tbl,
    function(Figure_ID, Legend_title, Legend_text, ...) {
      paste0(
        "## ", Figure_ID, ". ", Legend_title, "\n\n",
        Legend_text,
        "\n"
      )
    }
  )
  writeLines(lines, con = file)
}

###############################################################################
# 4) FINAL OUTPUT INVENTORY
###############################################################################

# NOTE:
# Script 06 and Script 07 may generate slightly different folder structures
# depending on the version used. This inventory checks the expected final file
# names and flags missing files without stopping the script.

inventory <- dplyr::bind_rows(

  # ---------------------------------------------------------------------------
  # Main WGCNA figure
  # ---------------------------------------------------------------------------
  make_file_record(
    figure_id = "Figure 3",
    figure_type = "Main figure",
    description = "Integrated WGCNA architecture: module structure, module size, prioritized module-trait associations, and module biological portraits.",
    file_base = "Figure3_WGCNA_main_real_networks_fixed_v2",
    directory = SCRIPT06_FIG_DIR,
    source_table = NA_character_,
    required = TRUE
  ),

  # ---------------------------------------------------------------------------
  # Supplementary Fig. 4
  # ---------------------------------------------------------------------------
  make_file_record(
    figure_id = "Supplementary Fig. 4",
    figure_type = "Supplementary figure",
    description = "Complete WGCNA module-trait association matrix across all non-grey modules and available clinical, cognitive, AT(N), demographic, and genetic variables.",
    file_base = "Supplementary_Fig_4_complete_module_trait_heatmap",
    directory = SCRIPT07_FIG_DIR,
    source_table = file.path(SCRIPT07_TAB_DIR, "Supplementary_Fig_4_complete_module_trait_heatmap_source.csv"),
    required = TRUE
  ),

  # ---------------------------------------------------------------------------
  # Supplementary Fig. 5
  # ---------------------------------------------------------------------------
  make_file_record(
    figure_id = "Supplementary Fig. 5",
    figure_type = "Supplementary figure",
    description = "Soft-thresholding diagnostics supporting WGCNA network construction and selected beta.",
    file_base = "Supplementary_Fig_5_soft_thresholding",
    directory = SCRIPT07_FIG_DIR,
    source_table = NA_character_,
    required = TRUE
  ),

  # ---------------------------------------------------------------------------
  # Supplementary Fig. 6
  # ---------------------------------------------------------------------------
  make_file_record(
    figure_id = "Supplementary Fig. 6",
    figure_type = "Supplementary figure",
    description = "Module size and differential-abundance burden across WGCNA modules.",
    file_base = "Supplementary_Fig_6_module_size_DEP_burden",
    directory = SCRIPT07_FIG_DIR,
    source_table = NA_character_,
    required = TRUE
  ),

  # ---------------------------------------------------------------------------
  # Supplementary Fig. 7
  # ---------------------------------------------------------------------------
  make_file_record(
    figure_id = "Supplementary Fig. 7",
    figure_type = "Supplementary figure",
    description = "Top hub proteins by module, ranked by absolute module membership.",
    file_base = "Supplementary_Fig_7_hub_proteins_by_module",
    directory = SCRIPT07_FIG_DIR,
    source_table = file.path(SCRIPT07_TAB_DIR, "Supplementary_Fig_7_top_hub_proteins_source.csv"),
    required = TRUE
  ),

  # ---------------------------------------------------------------------------
  # Supplementary Fig. 8
  # ---------------------------------------------------------------------------
  make_file_record(
    figure_id = "Supplementary Fig. 8",
    figure_type = "Supplementary figure",
    description = "Module membership and hub-structure overview based on kME distributions and high-membership hub burden.",
    file_base = "Supplementary_Fig_8_enrichment_summary",
    directory = SCRIPT07_FIG_DIR,
    source_table = NA_character_,
    required = TRUE
  ),

  # ---------------------------------------------------------------------------
  # Supplementary Fig. 9
  # ---------------------------------------------------------------------------
  make_file_record(
    figure_id = "Supplementary Fig. 9",
    figure_type = "Supplementary figure",
    description = "Functional enrichment term burden across modules for GO BP, KEGG, and Reactome.",
    file_base = "Supplementary_Fig_9_LOCO_country_robustness",
    directory = SCRIPT07_FIG_DIR,
    source_table = NA_character_,
    required = TRUE
  ),

  # ---------------------------------------------------------------------------
  # Supplementary Fig. 10
  # ---------------------------------------------------------------------------
  make_file_record(
    figure_id = "Supplementary Fig. 10",
    figure_type = "Supplementary figure",
    description = "Leave-one-country-out robustness of module-trait associations.",
    file_base = "Supplementary_Fig_10_downsampling_LOSO_robustness",
    directory = SCRIPT07_FIG_DIR,
    source_table = NA_character_,
    required = TRUE
  ),

  # ---------------------------------------------------------------------------
  # Supplementary Fig. 10
  # There are two possible filenames depending on whether LOSO exists.
  # ---------------------------------------------------------------------------
  make_file_record(
    figure_id = "Supplementary Fig. 10",
    figure_type = "Supplementary figure",
    description = "Balanced downsampling and optional leave-one-site-out robustness of module-trait associations.",
    file_base = "Supplementary_Fig_10_downsampling_LOSO_robustness",
    directory = SCRIPT07_FIG_DIR,
    source_table = NA_character_,
    required = FALSE
  ),

  make_file_record(
    figure_id = "Supplementary Fig. 10",
    figure_type = "Supplementary figure",
    description = "Balanced downsampling robustness of module-trait associations when leave-one-site-out output is unavailable.",
    file_base = "Supplementary_Fig_10_downsampling_robustness",
    directory = SCRIPT07_FIG_DIR,
    source_table = NA_character_,
    required = FALSE
  )
)

# For Supplementary Fig. 10, only one of the two alternative files is required.
supp11_exists <- inventory %>%
  dplyr::filter(Figure_ID == "Supplementary Fig. 10") %>%
  dplyr::summarise(any_exists = any(PDF_exists | PNG_exists, na.rm = TRUE)) %>%
  dplyr::pull(any_exists)

inventory <- inventory %>%
  dplyr::mutate(
    Required_effective = dplyr::case_when(
      Figure_ID == "Supplementary Fig. 10" ~ supp11_exists,
      TRUE ~ Required
    ),
    Missing_expected_output = dplyr::case_when(
      Figure_ID == "Supplementary Fig. 10" ~ !supp11_exists,
      Required ~ !(PDF_exists & PNG_exists),
      TRUE ~ FALSE
    )
  )

readr::write_csv(
  inventory,
  file.path(OUT_TAB, "wgcna_final_output_inventory.csv")
)

###############################################################################
# 5) FINAL FIGURE LEGENDS
###############################################################################

legend_tbl <- tibble::tribble(
  ~Figure_ID, ~Legend_title, ~Legend_text,

  "Figure 3",
  "Weighted co-expression network analysis identifies coordinated plasma proteomic modules associated with Alzheimer’s disease biology.",
  paste0(
    "(a) Gene dendrogram and merged module colors from the signed WGCNA network constructed using the gene-collapsed plasma proteomic universe. ",
    "(b) Distribution of module sizes across non-grey modules, with prioritized modules highlighted. ",
    "(c) Module-trait association heatmap for prioritized modules, showing Spearman correlations with clinical, cognitive, AT(N), demographic, and genetic variables; asterisks denote BH-FDR-corrected significance. ",
    "(d) Integrated module prioritization based on biological richness, clinical association, and hub connectivity. (e-f) Biological portraits of prioritized WGCNA modules, integrating hub proteins, differential-abundance composition, functional enrichment, and TOM-based local network structure. ",
    "Together, these panels summarize the coordinated module-level architecture underlying the AD-associated plasma proteomic signature."
  ),

  "Supplementary Fig. 4",
  "Complete WGCNA module-trait association matrix.",
  paste0(
    "Complete module-trait heatmap showing Spearman associations between all non-grey WGCNA module eigengenes and available clinical, cognitive, AT(N), demographic, and genetic variables. ",
    "Tile color indicates the direction and magnitude of Spearman rho. ",
    "Asterisks indicate BH-FDR-corrected significance across the full tested matrix. ",
    "The side color bar denotes WGCNA module identity."
  ),

  "Supplementary Fig. 5",
  "Soft-thresholding diagnostics for WGCNA network construction.",
  paste0(
    "(a) Scale-free topology model fit across candidate soft-thresholding powers. ",
    "The horizontal dashed line denotes the target criterion for approximate scale-free topology, and the vertical dashed line indicates the selected soft-thresholding power. ",
    "(b) Mean network connectivity across candidate powers, illustrating the trade-off between scale independence and network sparsity."
  ),

  "Supplementary Fig. 6",
  "Module size and differential-abundance burden across WGCNA modules.",
  paste0(
    "(a) Number of proteins assigned to each non-grey WGCNA module. ",
    "(b) Proportion of proteins within each module classified as differentially abundant at FDR < 0.05. ",
    "(c) Absolute number of differentially abundant proteins per module. ",
    "(d) Average absolute log2 fold-change across proteins within each module. ",
    "These panels summarize the extent to which module-level structure captures the differential-abundance signal."
  ),

  "Supplementary Fig. 7",
  "Top hub proteins across WGCNA modules.",
  paste0(
    "Top hub proteins within prioritized WGCNA modules ranked by absolute module membership (|kME|). ",
    "Each panel shows the highest-ranked hub proteins for a given module, with bars colored according to WGCNA module identity. ",
    "This figure complements the module portraits by showing the specific proteins driving module cohesion."
  ),

  "Supplementary Fig. 8",
  "Module membership and hub-structure overview.",
  paste0(
    "(a) Distribution of absolute module membership values (|kME|) across non-grey WGCNA modules. ",
    "(b) Number of high-membership hub proteins per module according to the available kME threshold. ",
    "These analyses characterize the internal cohesion and hub architecture of each module."
  ),

  "Supplementary Fig. 9",
  "Functional enrichment burden across WGCNA modules.",
  paste0(
    "(a) Number of enriched terms per module across GO Biological Process, KEGG, and Reactome libraries. ",
    "(b) Total enrichment burden per module. ",
    "These panels summarize the relative functional annotation density of each WGCNA module."
  ),

  "Supplementary Fig. 10",
  "Leave-one-country-out robustness of module-trait associations.",
  paste0(
    "(a) Mean absolute change in module-trait Spearman rho after excluding one country at a time. ",
    "(b) Maximum absolute change observed across leave-one-country-out iterations. ",
    "(c) Proportion of module-trait associations preserving the same direction as the full-sample estimate. ",
    "These analyses evaluate whether module-trait associations are stable across recruitment countries."
  ),

  "Supplementary Fig. 10",
  "Balanced downsampling and leave-one-site-out robustness of module-trait associations.",
  paste0(
    "(a) Distribution of changes in module-trait Spearman rho after balanced downsampling by country and diagnostic group. ",
    "(b) Direction preservation after downsampling relative to the full-sample estimates. ",
    "(c) Downsampling-related changes stratified by trait. ",
    "(d) If site information is available, leave-one-site-out robustness is summarized as mean absolute change in rho. ",
    "These analyses evaluate whether module-trait patterns are preserved under balanced resampling and site-level sensitivity checks."
  )
)

readr::write_csv(
  legend_tbl,
  file.path(OUT_TAB, "wgcna_final_figure_legends.csv")
)

write_markdown_legends(
  legend_tbl,
  file.path(OUT_TXT, "wgcna_final_figure_legends.md")
)

###############################################################################
# 6) SUBMISSION CHECKLIST
###############################################################################

checklist <- tibble::tibble(
  Check_item = c(
    "Main WGCNA Figure 3 exists as PDF and PNG",
    "Complete module-trait heatmap exists",
    "Soft-thresholding supplementary figure exists",
    "Module DEP burden supplementary figure exists",
    "Top hub proteins by module supplementary figure exists",
    "Hub/kME structure supplementary figure exists",
    "Functional enrichment supplementary figure exists",
    "LOCO robustness supplementary figure exists",
    "Downsampling/LOSO robustness supplementary figure exists",
    "Source table for complete module-trait heatmap exists",
    "Source table for top hub proteins exists",
    "Supplementary figure inventory from Script 07 exists",
    "All expected required PDF/PNG files exist"
  ),
  Status = c(
    any(inventory$Figure_ID == "Figure 3" & inventory$PDF_exists & inventory$PNG_exists),
    any(inventory$Figure_ID == "Supplementary Fig. 4" & inventory$PDF_exists & inventory$PNG_exists),
    any(inventory$Figure_ID == "Supplementary Fig. 5" & inventory$PDF_exists & inventory$PNG_exists),
    any(inventory$Figure_ID == "Supplementary Fig. 6" & inventory$PDF_exists & inventory$PNG_exists),
    any(inventory$Figure_ID == "Supplementary Fig. 7" & inventory$PDF_exists & inventory$PNG_exists),
    any(inventory$Figure_ID == "Supplementary Fig. 8" & inventory$PDF_exists & inventory$PNG_exists),
    any(inventory$Figure_ID == "Supplementary Fig. 9" & inventory$PDF_exists & inventory$PNG_exists),
    any(inventory$Figure_ID == "Supplementary Fig. 10" & inventory$PDF_exists & inventory$PNG_exists),
    any(inventory$Figure_ID == "Supplementary Fig. 10" & (inventory$PDF_exists | inventory$PNG_exists)),
    file.exists(file.path(SCRIPT07_TAB_DIR, "Supplementary_Fig_4_complete_module_trait_heatmap_source.csv")),
    file.exists(file.path(SCRIPT07_TAB_DIR, "Supplementary_Fig_7_top_hub_proteins_source.csv")),
    file.exists(file.path(SCRIPT07_TAB_DIR, "supplementary_figure_inventory.csv")),
    !any(inventory$Missing_expected_output, na.rm = TRUE)
  )
) %>%
  dplyr::mutate(
    Status_label = ifelse(Status, "OK", "CHECK")
  )

readr::write_csv(
  checklist,
  file.path(OUT_TAB, "wgcna_submission_checklist.csv")
)

###############################################################################
# 7) OPTIONAL: READ SCRIPT 07 INVENTORY IF AVAILABLE
###############################################################################

script07_inventory_file <- file.path(SCRIPT07_TAB_DIR, "supplementary_figure_inventory.csv")

if (file.exists(script07_inventory_file)) {
  script07_inventory <- readr::read_csv(script07_inventory_file, show_col_types = FALSE)
  readr::write_csv(
    script07_inventory,
    file.path(OUT_TAB, "script07_supplementary_figure_inventory_copy.csv")
  )
}

###############################################################################
# 8) CONSOLE SUMMARY
###############################################################################

cat("\nDONE.\n")
cat("Output inventory saved to:\n", file.path(OUT_TAB, "wgcna_final_output_inventory.csv"), "\n")
cat("Figure legends saved to:\n", file.path(OUT_TAB, "wgcna_final_figure_legends.csv"), "\n")
cat("Markdown legends saved to:\n", file.path(OUT_TXT, "wgcna_final_figure_legends.md"), "\n")
cat("Checklist saved to:\n", file.path(OUT_TAB, "wgcna_submission_checklist.csv"), "\n")

cat("\nMissing required outputs:\n")
missing_tbl <- inventory %>% dplyr::filter(Missing_expected_output)
if (nrow(missing_tbl) == 0) {
  cat("None.\n")
} else {
  print(missing_tbl %>% dplyr::select(Figure_ID, File_base, PDF_exists, PNG_exists, Required, Required_effective))
}

cat("\nChecklist:\n")
print(checklist %>% dplyr::select(Check_item, Status_label))


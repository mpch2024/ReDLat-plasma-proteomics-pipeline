###############################################################################
# 04_WGCNA_module_trait_clinical_integration_FINAL_v6.R
#
# TITLE
# Integrated module-trait, clinical, DEP burden, robustness and prioritization
# analysis for GENE-COLLAPSED WGCNA plasma modules
#
# PURPOSE
# This script integrates outputs from Scripts 02 and 03 with metadata and DEP
# results to generate:
#   1) module-trait Spearman correlations
#   2) BH-FDR correction across the module-trait matrix
#   3) publication-ready heatmap
#   4) DEP burden by module
#   5) leave-one-country-out robustness
#   6) focalized module-trait scatter plots
#   7) adjusted module-level regression models
#   8) integrated final module prioritization table
#
# IMPORTANT
# This script assumes WGCNA was constructed at the GENE-COLLAPSED level.
# It does NOT use aptamer-level WGCNA modules.
#
# INPUTS
#   Script 02:
#     results/02_WGCNA_core_collapsed_genes/eigengenes/module_eigengenes_per_sample.csv
#     results/02_WGCNA_core_collapsed_genes/tables/gene_module_assignment.csv
#
#   Script 03:
#     results/03_WGCNA_module_biology_hubs_enrichment/tables/module_hub_summary.csv
#     results/03_WGCNA_module_biology_hubs_enrichment/tables/enrichment_summary_by_module.csv
#
#   DEP clean pipeline:
#     C:/Users/mnpiz/Desktop/DEPs_Proteomic_Publishable_V2/result/03_dep/gene_collapsed/
#       AD_vs_CN_full_limma_results_gene_collapsed.csv
#
#   Metadata:
#     results/01_define_wgcna_input_from_DEP/wgcna_sample_metadata.csv
#
# OUTPUTS
#   results/04_WGCNA_module_trait_clinical_integration/
#
# UPDATE v6
# Compatible with Script 03 v3 robust enrichment outputs.
#
# AUTHOR
# Matías Pizarro + ChatGPT support
###############################################################################

rm(list = ls())

###############################################################################
# 1) PACKAGES
###############################################################################

cran_pkgs <- c(
  "dplyr",
  "tidyr",
  "readr",
  "tibble",
  "stringr",
  "purrr",
  "ggplot2",
  "pheatmap",
  "ggrepel",
  "scales",
  "broom",
  "circlize"
)

bioc_pkgs <- c(
  "ComplexHeatmap"
)

cran_missing <- cran_pkgs[
  !sapply(cran_pkgs, requireNamespace, quietly = TRUE)
]

if (length(cran_missing) > 0) {
  install.packages(cran_missing)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_missing <- bioc_pkgs[
  !sapply(bioc_pkgs, requireNamespace, quietly = TRUE)
]

if (length(bioc_missing) > 0) {
  BiocManager::install(bioc_missing, ask = FALSE, update = FALSE)
}

invisible(lapply(c(cran_pkgs, bioc_pkgs), library, character.only = TRUE))

options(stringsAsFactors = FALSE)
options(error = traceback)

###############################################################################
# 2) PATHS
###############################################################################

BASE_DIR <- "C:/Users/mnpiz/Desktop/WGCNA_Workflow_april_V2"
DEP_PROJECT_ROOT <- "C:/Users/mnpiz/Desktop/DEPs_Proteomic_Publishable_V2"

SCRIPT1_DIR <- file.path(
  BASE_DIR,
  "results",
  "01_define_wgcna_input_from_DEP"
)

SCRIPT2_DIR <- file.path(
  BASE_DIR,
  "results",
  "02_WGCNA_core_collapsed_genes"
)

SCRIPT3_DIR <- file.path(
  BASE_DIR,
  "results",
  "03_WGCNA_module_biology_hubs_enrichment"
)

EIGENGENE_FILE <- file.path(
  SCRIPT2_DIR,
  "eigengenes",
  "module_eigengenes_per_sample.csv"
)

GENE_MODULE_ASSIGNMENT_FILE <- file.path(
  SCRIPT2_DIR,
  "tables",
  "gene_module_assignment.csv"
)

META_FILE <- file.path(
  SCRIPT1_DIR,
  "wgcna_sample_metadata.csv"
)

HUB_SUMMARY_FILE <- file.path(
  SCRIPT3_DIR,
  "tables",
  "module_hub_summary.csv"
)

ENRICHMENT_SUMMARY_FILE <- file.path(
  SCRIPT3_DIR,
  "tables",
  "enrichment_summary_by_module.csv"
)

DEP_FILE <- file.path(
  DEP_PROJECT_ROOT,
  "result",
  "03_dep",
  "gene_collapsed",
  "AD_vs_CN_full_limma_results_gene_collapsed.csv"
)

OUTDIR <- file.path(
  BASE_DIR,
  "results",
  "04_WGCNA_module_trait_clinical_integration"
)

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

SUBDIRS <- c(
  "tables",
  "figures",
  "figures/module_trait",
  "figures/dep_burden",
  "figures/loco",
  "figures/scatter",
  "tables/loco",
  "tables/regression",
  "workspace"
)

invisible(lapply(
  file.path(OUTDIR, SUBDIRS),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

REQUIRED_FILES <- c(
  EIGENGENE_FILE,
  GENE_MODULE_ASSIGNMENT_FILE,
  META_FILE,
  HUB_SUMMARY_FILE,
  ENRICHMENT_SUMMARY_FILE,
  DEP_FILE
)

missing_files <- REQUIRED_FILES[!file.exists(REQUIRED_FILES)]

if (length(missing_files) > 0) {
  stop(
    "Faltan archivos requeridos:\n",
    paste(missing_files, collapse = "\n"),
    call. = FALSE
  )
}

###############################################################################
# 3) PARAMETERS
###############################################################################

MIN_N_FOR_CORR <- 6
MIN_N_FOR_REG <- 15

MODULES_OF_INTEREST <- NULL
EXCLUDE_MODULES <- c("grey", "MEgrey")

AD_TRAITS_ORDER <- c(
  "SampleGroup_bin",
  "cdr_global",
  "cdr_boxscore",
  "mmse_total",
  "udsfaq_total",
  "NPI",
  "Mini_SEA",
  "T_ADLQ",
  "p_tau181",
  "p_tau217",
  "NfL",
  "ratio_AB42_40"
)

COV_TRAITS_ORDER <- c(
  "Age",
  "Sex_bin",
  "Education",
  "APOE4_carrier",
  "Country_numeric"
)

MAIN_TRAITS_FOR_INTEGRATION <- c(
  "SampleGroup_bin",
  "cdr_boxscore",
  "mmse_total",
  "udsfaq_total",
  "p_tau181",
  "p_tau217",
  "NfL",
  "ratio_AB42_40",
  "Age"
)

SCATTER_PRIORITY_TRAITS <- c(
  "SampleGroup_bin",
  "cdr_boxscore",
  "mmse_total",
  "udsfaq_total",
  "p_tau181",
  "p_tau217",
  "NfL",
  "ratio_AB42_40",
  "Age"
)

TOP_SCATTERS_PER_MODULE <- 4

REGRESSION_OUTCOMES <- c(
  "p_tau181",
  "p_tau217",
  "NfL",
  "ratio_AB42_40",
  "cdr_boxscore",
  "mmse_total"
)

TRAIT_LABELS <- c(
  SampleGroup_bin = "Diagnosis",
  cdr_global = "CDR global",
  cdr_boxscore = "CDR-SB",
  mmse_total = "MMSE",
  udsfaq_total = "PFAQ",
  NPI = "NPI-Q",
  Mini_SEA = "Mini-SEA",
  T_ADLQ = "T-ADLQ",
  p_tau181 = "p-tau181",
  p_tau217 = "p-tau217",
  NfL = "NfL",
  ratio_AB42_40 = "Aβ42/40",
  Age = "Age",
  Sex_bin = "Sex",
  Education = "Education",
  APOE4_carrier = "APOE ε4",
  Country_numeric = "Country"
)

MODULE_COLORS <- c(
  black       = "#2B2B2B",
  brown       = "#9C6B00",
  yellow      = "#B7D500",
  blue        = "#1535E8",
  green       = "#3E8F4E",
  red         = "#B4334A",
  purple      = "#8A2BE2",
  magenta     = "#C75ACD",
  pink        = "#F3A6D6",
  greenyellow = "#ADFF2F",
  grey        = "#9E9E9E"
)

get_module_colors <- function(modules) {
  modules <- as.character(modules)
  cols <- MODULE_COLORS[modules]
  missing_mods <- modules[is.na(cols)]

  if (length(missing_mods) > 0) {
    extra_cols <- grDevices::rainbow(length(unique(missing_mods)), s = 0.45, v = 0.75)
    names(extra_cols) <- unique(missing_mods)
    cols[is.na(cols)] <- extra_cols[missing_mods]
  }

  names(cols) <- modules
  cols
}

DPI <- 300
HEATMAP_WIDTH <- 12
HEATMAP_HEIGHT <- 7

###############################################################################
# 4) HELPERS
###############################################################################

safe_write_csv <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(x, file)
}

clean_text_na <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA", "NaN", "NULL", "null", "N/A", "nan")] <- NA_character_
  x
}

standardize_module_name <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("^ME", "", x)
  x
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(x))
}

add_sig_stars <- function(q) {
  dplyr::case_when(
    is.na(q) ~ "",
    q < 0.001 ~ "***",
    q < 0.01 ~ "**",
    q < 0.05 ~ "*",
    TRUE ~ ""
  )
}

safe_spearman <- function(df, xcol, ycol, min_n = 6) {
  d <- df %>%
    dplyr::select(dplyr::all_of(c(xcol, ycol))) %>%
    tidyr::drop_na()

  n <- nrow(d)

  if (n < min_n) {
    return(tibble::tibble(rho = NA_real_, p_value = NA_real_, N = n))
  }

  if (length(unique(d[[xcol]])) <= 1 || length(unique(d[[ycol]])) <= 1) {
    return(tibble::tibble(rho = NA_real_, p_value = NA_real_, N = n))
  }

  test <- suppressWarnings(
    stats::cor.test(d[[xcol]], d[[ycol]], method = "spearman", exact = FALSE)
  )

  tibble::tibble(
    rho = unname(test$estimate),
    p_value = test$p.value,
    N = n
  )
}

standardize_metadata_columns <- function(df) {
  rename_map <- c(
    "Mini-SEA" = "Mini_SEA",
    "Mini.SEA" = "Mini_SEA",
    "T-ADLQ" = "T_ADLQ",
    "T.ADLQ" = "T_ADLQ",
    "p-tau181" = "p_tau181",
    "p.tau181" = "p_tau181",
    "p-tau217" = "p_tau217",
    "p.tau217" = "p_tau217",
    "ratio AB42/40" = "ratio_AB42_40",
    "ratio.AB42.40" = "ratio_AB42_40"
  )

  for (old in names(rename_map)) {
    new <- rename_map[[old]]
    if (old %in% names(df) && !new %in% names(df)) {
      names(df)[names(df) == old] <- new
    }
  }

  df
}

encode_basic_variables <- function(df) {
  df <- df %>%
    dplyr::mutate(SampleId = as.character(SampleId))

  if ("SampleGroup" %in% names(df)) {
    df <- df %>%
      dplyr::mutate(
        SampleGroup = as.character(SampleGroup),
        SampleGroup_bin = dplyr::case_when(
          SampleGroup == "CN" ~ 0,
          SampleGroup == "AD" ~ 1,
          TRUE ~ NA_real_
        )
      )
  }

  if ("Sex" %in% names(df)) {
    df <- df %>%
      dplyr::mutate(
        Sex_chr = as.character(Sex),
        Sex_bin = dplyr::case_when(
          Sex_chr %in% c("1", "M", "Male", "male") ~ 0,
          Sex_chr %in% c("2", "F", "Female", "female") ~ 1,
          TRUE ~ safe_numeric(Sex_chr)
        )
      )
  }

  if ("ApoE" %in% names(df) && !"APOE_group" %in% names(df)) {
    df <- df %>%
      dplyr::mutate(
        ApoE = as.character(ApoE),
        APOE_group = dplyr::case_when(
          ApoE %in% c("e2/e4", "e3/e4", "e4/e4") ~ "E4 carrier",
          ApoE %in% c("e2/e2", "e2/e3", "e3/e3") ~ "Non-E4",
          TRUE ~ NA_character_
        )
      )
  }

  if ("APOE_group" %in% names(df) && !"APOE4_carrier" %in% names(df)) {
    df <- df %>%
      dplyr::mutate(
        APOE4_carrier = dplyr::case_when(
          APOE_group == "E4 carrier" ~ 1,
          APOE_group == "Non-E4" ~ 0,
          TRUE ~ NA_real_
        )
      )
  }

  if ("Country" %in% names(df)) {
    df <- df %>%
      dplyr::mutate(
        Country = as.factor(Country),
        Country_numeric = as.numeric(Country)
      )
  }

  numeric_candidates <- c(
    "Age", "Education", "cdr_global", "cdr_boxscore", "mmse_total",
    "udsfaq_total", "NPI", "Mini_SEA", "T_ADLQ", "p_tau181",
    "p_tau217", "NfL", "ratio_AB42_40", "APOE4_carrier"
  )

  for (cc in intersect(numeric_candidates, names(df))) {
    df[[cc]] <- safe_numeric(df[[cc]])
  }

  df
}

standardize_dep_table <- function(dep_df) {
  rename_map <- c(
    "adj.P.Val" = "FDR",
    "adj_p_val" = "FDR",
    "adj_pvalue" = "FDR",
    "adj_p_value" = "FDR",
    "fdr" = "FDR",
    "q_value" = "FDR",
    "logFC" = "logFC",
    "logfc" = "logFC",
    "Gene" = "EntrezGeneSymbol",
    "gene" = "EntrezGeneSymbol",
    "symbol" = "EntrezGeneSymbol"
  )

  for (old in names(rename_map)) {
    new <- rename_map[[old]]
    if (old %in% names(dep_df) && !new %in% names(dep_df)) {
      names(dep_df)[names(dep_df) == old] <- new
    }
  }

  needed <- c("EntrezGeneSymbol", "FDR", "logFC")
  missing <- setdiff(needed, names(dep_df))

  if (length(missing) > 0) {
    stop(
      "Faltan columnas en DEP final: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  dep_df %>%
    dplyr::mutate(
      EntrezGeneSymbol = clean_text_na(EntrezGeneSymbol),
      FDR = safe_numeric(FDR),
      logFC = safe_numeric(logFC)
    ) %>%
    dplyr::filter(!is.na(EntrezGeneSymbol))
}

compute_module_trait <- function(df, modules, traits, label = "full") {
  expand.grid(
    Module = modules,
    Trait = traits,
    stringsAsFactors = FALSE
  ) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
      result = purrr::map2(
        Module,
        Trait,
        ~ safe_spearman(df, .x, .y, min_n = MIN_N_FOR_CORR)
      )
    ) %>%
    tidyr::unnest(result) %>%
    dplyr::mutate(Analysis = label)
}

make_delta_summary <- function(full_tbl, sens_tbl) {
  sens_tbl %>%
    dplyr::left_join(
      full_tbl %>%
        dplyr::select(Module, Trait, rho_full = rho, p_full = p_value, N_full = N),
      by = c("Module", "Trait")
    ) %>%
    dplyr::mutate(
      delta_rho = rho - rho_full,
      abs_delta_rho = abs(delta_rho),
      direction_full = sign(rho_full),
      direction_sensitivity = sign(rho),
      same_direction = direction_full == direction_sensitivity
    )
}

###############################################################################
# 5) LOAD DATA
###############################################################################

eig <- readr::read_csv(EIGENGENE_FILE, show_col_types = FALSE)
meta <- readr::read_csv(META_FILE, show_col_types = FALSE)
gene_module_assignment <- readr::read_csv(GENE_MODULE_ASSIGNMENT_FILE, show_col_types = FALSE)
hub_summary <- readr::read_csv(HUB_SUMMARY_FILE, show_col_types = FALSE)
enrichment_summary <- readr::read_csv(ENRICHMENT_SUMMARY_FILE, show_col_types = FALSE)
dep <- readr::read_csv(DEP_FILE, show_col_types = FALSE)

cat("INPUT DIMENSIONS\n")
cat("eigengenes        :", dim(eig), "\n")
cat("metadata          :", dim(meta), "\n")
cat("gene modules      :", dim(gene_module_assignment), "\n")
cat("hub summary       :", dim(hub_summary), "\n")
cat("enrichment summary:", dim(enrichment_summary), "\n")
cat("DEP table         :", dim(dep), "\n\n")

###############################################################################
# 6) PREPARE METADATA + EIGENGENES
###############################################################################

eig <- eig %>%
  dplyr::mutate(SampleId = as.character(SampleId)) %>%
  dplyr::distinct(SampleId, .keep_all = TRUE)

meta <- meta %>%
  standardize_metadata_columns() %>%
  encode_basic_variables() %>%
  dplyr::distinct(SampleId, .keep_all = TRUE)

eig_cols <- setdiff(names(eig), "SampleId")
eig_cols <- eig_cols[!eig_cols %in% EXCLUDE_MODULES]

# Ensure module eigengene columns use ME prefix and then create clean module names.
eig_renamed <- eig
names(eig_renamed)[names(eig_renamed) %in% eig_cols] <- paste0(
  "ME",
  standardize_module_name(eig_cols)
)

module_cols <- setdiff(names(eig_renamed), "SampleId")
module_cols <- module_cols[!module_cols %in% c("MEgrey", "grey")]

analysis_df <- eig_renamed %>%
  dplyr::inner_join(meta, by = "SampleId")

if (nrow(analysis_df) == 0) {
  stop("No hay muestras compartidas entre eigengenes y metadata.", call. = FALSE)
}

for (mc in module_cols) {
  analysis_df[[mc]] <- safe_numeric(analysis_df[[mc]])
}

module_names_clean <- standardize_module_name(module_cols)
names(analysis_df)[match(module_cols, names(analysis_df))] <- module_names_clean

module_cols_clean <- module_names_clean

if (is.null(MODULES_OF_INTEREST)) {
  modules_use <- setdiff(module_cols_clean, standardize_module_name(EXCLUDE_MODULES))
} else {
  modules_use <- intersect(MODULES_OF_INTEREST, module_cols_clean)
}

traits_available <- intersect(
  c(AD_TRAITS_ORDER, COV_TRAITS_ORDER),
  names(analysis_df)
)

if (length(traits_available) == 0) {
  stop("No hay traits disponibles para correlacionar.", call. = FALSE)
}

safe_write_csv(
  analysis_df,
  file.path(OUTDIR, "module_trait_input_clean.csv")
)

safe_write_csv(
  analysis_df,
  file.path(OUTDIR, "tables", "module_trait_input_clean.csv")
)

cat("Muestras en análisis:", nrow(analysis_df), "\n")
cat("Módulos analizados:", paste(modules_use, collapse = ", "), "\n")
cat("Traits disponibles:", paste(traits_available, collapse = ", "), "\n\n")

###############################################################################
# 7) MODULE-TRAIT CORRELATIONS
###############################################################################

module_trait_long <- compute_module_trait(
  df = analysis_df,
  modules = modules_use,
  traits = traits_available,
  label = "full"
) %>%
  dplyr::mutate(
    FDR = p.adjust(p_value, method = "BH"),
    stars = add_sig_stars(FDR),
    Trait_label = dplyr::recode(Trait, !!!TRAIT_LABELS, .default = Trait)
  )

safe_write_csv(
  module_trait_long,
  file.path(OUTDIR, "module_trait_results_long.csv")
)

safe_write_csv(
  module_trait_long,
  file.path(OUTDIR, "tables", "module_trait_results_long.csv")
)

rho_mat <- module_trait_long %>%
  dplyr::select(Module, Trait, rho) %>%
  tidyr::pivot_wider(names_from = Trait, values_from = rho) %>%
  tibble::column_to_rownames("Module") %>%
  as.matrix()

p_mat <- module_trait_long %>%
  dplyr::select(Module, Trait, p_value) %>%
  tidyr::pivot_wider(names_from = Trait, values_from = p_value) %>%
  tibble::column_to_rownames("Module") %>%
  as.matrix()

fdr_mat <- module_trait_long %>%
  dplyr::select(Module, Trait, FDR) %>%
  tidyr::pivot_wider(names_from = Trait, values_from = FDR) %>%
  tibble::column_to_rownames("Module") %>%
  as.matrix()

annot_mat <- module_trait_long %>%
  dplyr::select(Module, Trait, stars) %>%
  tidyr::pivot_wider(names_from = Trait, values_from = stars) %>%
  tibble::column_to_rownames("Module") %>%
  as.matrix()

safe_write_csv(
  as.data.frame(rho_mat) %>% tibble::rownames_to_column("Module"),
  file.path(OUTDIR, "module_trait_correlations_final.csv")
)

safe_write_csv(
  as.data.frame(p_mat) %>% tibble::rownames_to_column("Module"),
  file.path(OUTDIR, "module_trait_pvalues_final.csv")
)

safe_write_csv(
  as.data.frame(fdr_mat) %>% tibble::rownames_to_column("Module"),
  file.path(OUTDIR, "module_trait_fdr_final.csv")
)

safe_write_csv(
  as.data.frame(annot_mat) %>% tibble::rownames_to_column("Module"),
  file.path(OUTDIR, "module_trait_annotations_final.csv")
)

###############################################################################
# 8) MODULE-TRAIT HEATMAPS
###############################################################################

col_order <- intersect(c(AD_TRAITS_ORDER, COV_TRAITS_ORDER), colnames(rho_mat))

rho_mat_plot <- rho_mat[, col_order, drop = FALSE]
annot_mat_plot <- annot_mat[, col_order, drop = FALSE]

colnames(rho_mat_plot) <- dplyr::recode(
  colnames(rho_mat_plot),
  !!!TRAIT_LABELS,
  .default = colnames(rho_mat_plot)
)

colnames(annot_mat_plot) <- colnames(rho_mat_plot)

# -------------------------------------------------------------------------
# 8A) Classic pheatmap version, kept for compatibility
# -------------------------------------------------------------------------

pdf(
  file.path(OUTDIR, "figures", "module_trait", "module_trait_heatmap_fdr.pdf"),
  width = HEATMAP_WIDTH,
  height = HEATMAP_HEIGHT
)

pheatmap::pheatmap(
  rho_mat_plot,
  color = colorRampPalette(c("#4682b4", "white", "#f46d43"))(101),
  breaks = seq(-1, 1, length.out = 102),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  display_numbers = annot_mat_plot,
  number_color = "black",
  fontsize = 10,
  fontsize_number = 12,
  border_color = "grey90",
  main = "WGCNA module-trait associations\nSpearman rho; stars indicate BH-FDR significance"
)

dev.off()

png(
  file.path(OUTDIR, "figures", "module_trait", "module_trait_heatmap_fdr.png"),
  width = HEATMAP_WIDTH,
  height = HEATMAP_HEIGHT,
  units = "in",
  res = DPI
)

pheatmap::pheatmap(
  rho_mat_plot,
  color = colorRampPalette(c("#4682b4", "white", "#f46d43"))(101),
  breaks = seq(-1, 1, length.out = 102),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  display_numbers = annot_mat_plot,
  number_color = "black",
  fontsize = 10,
  fontsize_number = 12,
  border_color = "grey90",
  main = "WGCNA module-trait associations\nSpearman rho; stars indicate BH-FDR significance"
)

dev.off()

# -------------------------------------------------------------------------
# 8B) Compact manuscript-style heatmap using ComplexHeatmap
#     Visual target: module color bar at left, clean tiles, compact legends.
# -------------------------------------------------------------------------

module_rows <- rownames(rho_mat_plot)

module_palette <- get_module_colors(module_rows)

module_order_tbl <- tibble::tibble(Module = module_rows) %>%
  dplyr::mutate(Module_ID = paste0("M", dplyr::row_number()))

module_row_labels <- paste0(module_order_tbl$Module, "\n", module_order_tbl$Module_ID)

# Soft red/blue palette closer to the example.
# The range is intentionally narrower than -1 to 1 because most module-trait
# Spearman correlations are modest; this improves visible contrast without
# changing the actual rho values.
col_fun <- circlize::colorRamp2(
  c(-0.30, 0, 0.30),
  c("#6FA8DC", "#F7F7F7", "#F4A08B")
)

row_ha <- ComplexHeatmap::rowAnnotation(
  Module = ComplexHeatmap::anno_simple(
    module_rows,
    col = module_palette,
    border = FALSE
  ),
  show_annotation_name = FALSE,
  width = grid::unit(4, "mm")
)

ht <- ComplexHeatmap::Heatmap(
  rho_mat_plot,
  name = "Spearman\nrho",
  col = col_fun,

  cluster_rows = FALSE,
  cluster_columns = FALSE,

  left_annotation = row_ha,

  rect_gp = grid::gpar(
    col = "#E6E6E6",
    lwd = 0.45
  ),

  row_labels = module_row_labels,
  row_names_side = "left",
  row_names_gp = grid::gpar(
    fontsize = 9,
    fontface = "plain",
    col = "black",
    lineheight = 0.85
  ),

  column_names_gp = grid::gpar(
    fontsize = 9,
    fontface = "plain",
    col = "black"
  ),
  column_names_rot = 45,
  column_names_centered = FALSE,

  cell_fun = function(j, i, x, y, width, height, fill) {
    if (!is.na(annot_mat_plot[i, j]) && annot_mat_plot[i, j] != "") {
      grid::grid.text(
        annot_mat_plot[i, j],
        x,
        y,
        gp = grid::gpar(
          fontsize = 7,
          col = "black",
          fontface = "bold"
        )
      )
    }
  },

  heatmap_legend_param = list(
    title = "Spearman\nrho",
    title_gp = grid::gpar(fontsize = 9),
    labels_gp = grid::gpar(fontsize = 8),
    grid_width = grid::unit(3.5, "mm"),
    grid_height = grid::unit(3.5, "mm"),
    legend_height = grid::unit(28, "mm"),
    at = c(-0.2, 0, 0.2),
    labels = c("-0.2", "0.0", "0.2")
  ),

  width = grid::unit(16.5, "cm"),
  height = grid::unit(6.2, "cm")
)

pdf(
  file.path(OUTDIR, "figures", "module_trait", "module_trait_heatmap_manuscript_style.pdf"),
  width = 10.2,
  height = 5.4,
  useDingbats = FALSE,
  bg = "white"
)

grid::grid.newpage()
grid::grid.text(
  "Complete WGCNA module–trait association matrix",
  x = grid::unit(0.50, "npc"),
  y = grid::unit(0.965, "npc"),
  gp = grid::gpar(fontsize = 12, fontface = "plain")
)

ComplexHeatmap::draw(
  ht,
  heatmap_legend_side = "right",
  padding = grid::unit(c(12, 2, 2, 2), "mm")
)

dev.off()

png(
  file.path(OUTDIR, "figures", "module_trait", "module_trait_heatmap_manuscript_style.png"),
  width = 3000,
  height = 1600,
  res = 300,
  bg = "white"
)

grid::grid.newpage()
grid::grid.text(
  "Complete WGCNA module–trait association matrix",
  x = grid::unit(0.50, "npc"),
  y = grid::unit(0.965, "npc"),
  gp = grid::gpar(fontsize = 12, fontface = "plain")
)

ComplexHeatmap::draw(
  ht,
  heatmap_legend_side = "right",
  padding = grid::unit(c(12, 2, 2, 2), "mm")
)

dev.off()

# Keep a copy under the previous Nature-style filename for downstream scripts
# that may expect this name.
file.copy(
  from = file.path(OUTDIR, "figures", "module_trait", "module_trait_heatmap_manuscript_style.pdf"),
  to = file.path(OUTDIR, "figures", "module_trait", "module_trait_heatmap_nature_style.pdf"),
  overwrite = TRUE
)

file.copy(
  from = file.path(OUTDIR, "figures", "module_trait", "module_trait_heatmap_manuscript_style.png"),
  to = file.path(OUTDIR, "figures", "module_trait", "module_trait_heatmap_nature_style.png"),
  overwrite = TRUE
)

###############################################################################
# 9) DEP BURDEN BY MODULE
###############################################################################

dep_std <- standardize_dep_table(dep)

gene_module_small <- gene_module_assignment %>%
  dplyr::mutate(
    EntrezGeneSymbol = clean_text_na(EntrezGeneSymbol),
    Module = standardize_module_name(Module)
  ) %>%
  dplyr::filter(!is.na(EntrezGeneSymbol), !is.na(Module)) %>%
  dplyr::distinct(EntrezGeneSymbol, .keep_all = TRUE)

dep_module <- gene_module_small %>%
  dplyr::select(EntrezGeneSymbol, Module) %>%
  dplyr::left_join(
    dep_std %>%
      dplyr::select(EntrezGeneSymbol, FDR, logFC),
    by = "EntrezGeneSymbol"
  )

module_dep_summary <- dep_module %>%
  dplyr::group_by(Module) %>%
  dplyr::summarise(
    n_genes_module = dplyr::n(),
    n_dep_fdr_0_05 = sum(!is.na(FDR) & FDR < 0.05, na.rm = TRUE),
    n_dep_fdr_0_01 = sum(!is.na(FDR) & FDR < 0.01, na.rm = TRUE),
    prop_dep_fdr_0_05 = n_dep_fdr_0_05 / n_genes_module,
    prop_dep_fdr_0_01 = n_dep_fdr_0_01 / n_genes_module,
    mean_abs_logFC = mean(abs(logFC), na.rm = TRUE),
    median_abs_logFC = median(abs(logFC), na.rm = TRUE),
    n_up_fdr_0_05 = sum(!is.na(FDR) & FDR < 0.05 & logFC > 0, na.rm = TRUE),
    n_down_fdr_0_05 = sum(!is.na(FDR) & FDR < 0.05 & logFC < 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(prop_dep_fdr_0_05), dplyr::desc(n_dep_fdr_0_05))

safe_write_csv(
  module_dep_summary,
  file.path(OUTDIR, "tables", "module_dep_summary.csv")
)

p_dep <- ggplot(module_dep_summary, aes(x = reorder(Module, prop_dep_fdr_0_05), y = prop_dep_fdr_0_05, fill = Module)) +
  geom_col(width = 0.8, show.legend = FALSE) +
  scale_fill_manual(values = get_module_colors(module_dep_summary$Module)) +
  coord_flip() +
  labs(
    title = "DEP burden by WGCNA module",
    x = "Module",
    y = "Proportion of genes with DEP FDR < 0.05"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.major.y = element_blank())

ggsave(
  file.path(OUTDIR, "figures", "dep_burden", "module_dep_burden_prop_fdr005.pdf"),
  p_dep,
  width = 8,
  height = 5
)

ggsave(
  file.path(OUTDIR, "figures", "dep_burden", "module_dep_burden_prop_fdr005.png"),
  p_dep,
  width = 8,
  height = 5,
  dpi = DPI
)

###############################################################################
# 10) LEAVE-ONE-COUNTRY-OUT ROBUSTNESS
###############################################################################

loco_results <- list()

if ("Country" %in% names(analysis_df)) {
  countries <- sort(unique(as.character(analysis_df$Country)))
  countries <- countries[!is.na(countries)]

  full_tbl <- module_trait_long %>%
    dplyr::select(Module, Trait, rho, p_value, N)

  for (cc in countries) {
    df_cc <- analysis_df %>%
      dplyr::filter(as.character(Country) != cc)

    sens_tbl <- compute_module_trait(
      df = df_cc,
      modules = modules_use,
      traits = traits_available,
      label = paste0("leave_out_", cc)
    ) %>%
      dplyr::mutate(Excluded_country = cc)

    loco_results[[cc]] <- sens_tbl
  }

  loco_long <- dplyr::bind_rows(loco_results)

  loco_delta <- make_delta_summary(
    full_tbl = full_tbl,
    sens_tbl = loco_long
  )

  loco_summary_by_module_trait <- loco_delta %>%
    dplyr::group_by(Module, Trait) %>%
    dplyr::summarise(
      mean_abs_delta_rho = mean(abs_delta_rho, na.rm = TRUE),
      max_abs_delta_rho = max(abs_delta_rho, na.rm = TRUE),
      direction_consistency = mean(same_direction, na.rm = TRUE),
      n_sensitivity_runs = dplyr::n(),
      .groups = "drop"
    )

  loco_summary_by_module <- loco_delta %>%
    dplyr::group_by(Module) %>%
    dplyr::summarise(
      mean_abs_delta_rho = mean(abs_delta_rho, na.rm = TRUE),
      max_abs_delta_rho = max(abs_delta_rho, na.rm = TRUE),
      direction_consistency = mean(same_direction, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(mean_abs_delta_rho)

  safe_write_csv(
    loco_long,
    file.path(OUTDIR, "tables", "loco", "loco_country_module_trait_results.csv")
  )

  safe_write_csv(
    loco_delta,
    file.path(OUTDIR, "tables", "loco", "loco_country_delta_results.csv")
  )

  safe_write_csv(
    loco_summary_by_module_trait,
    file.path(OUTDIR, "tables", "loco", "loco_country_summary_by_module_trait.csv")
  )

  safe_write_csv(
    loco_summary_by_module,
    file.path(OUTDIR, "tables", "loco", "loco_country_summary_by_module.csv")
  )

  p_loco <- ggplot(loco_summary_by_module, aes(x = reorder(Module, mean_abs_delta_rho), y = mean_abs_delta_rho, fill = Module)) +
    geom_col(show.legend = FALSE) +
    scale_fill_manual(values = get_module_colors(loco_summary_by_module$Module)) +
    coord_flip() +
    labs(
      title = "Leave-one-country-out robustness",
      x = "Module",
      y = "Mean absolute delta rho"
    ) +
    theme_bw(base_size = 12)

  ggsave(
    file.path(OUTDIR, "figures", "loco", "loco_country_mean_abs_delta_rho.pdf"),
    p_loco,
    width = 7,
    height = 5
  )

  ggsave(
    file.path(OUTDIR, "figures", "loco", "loco_country_mean_abs_delta_rho.png"),
    p_loco,
    width = 7,
    height = 5,
    dpi = DPI
  )

} else {
  loco_summary_by_module <- tibble::tibble(
    Module = modules_use,
    mean_abs_delta_rho = NA_real_,
    max_abs_delta_rho = NA_real_,
    direction_consistency = NA_real_
  )
}

###############################################################################
# 11) FOCALIZED SCATTER PLOTS
###############################################################################

scatter_candidates <- module_trait_long %>%
  dplyr::filter(
    Trait %in% SCATTER_PRIORITY_TRAITS,
    !is.na(rho)
  ) %>%
  dplyr::group_by(Module) %>%
  dplyr::arrange(FDR, dplyr::desc(abs(rho)), .by_group = TRUE) %>%
  dplyr::slice(seq_len(min(TOP_SCATTERS_PER_MODULE, dplyr::n()))) %>%
  dplyr::ungroup()

safe_write_csv(
  scatter_candidates,
  file.path(OUTDIR, "tables", "scatter_candidates_top_module_trait.csv")
)

for (ii in seq_len(nrow(scatter_candidates))) {
  row <- scatter_candidates[ii, ]
  mod <- row$Module
  trait <- row$Trait
  trait_lab <- dplyr::recode(trait, !!!TRAIT_LABELS, .default = trait)

  if (!mod %in% names(analysis_df) || !trait %in% names(analysis_df)) next

  df_plot <- analysis_df %>%
    dplyr::select(SampleId, dplyr::all_of(c(mod, trait))) %>%
    tidyr::drop_na()

  if (nrow(df_plot) < MIN_N_FOR_CORR) next

  p <- ggplot(df_plot, aes(x = .data[[trait]], y = .data[[mod]])) +
    geom_point(alpha = 0.7, size = 2) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.6) +
    labs(
      title = paste0("Module ", mod, " vs ", trait_lab),
      subtitle = paste0("Spearman rho = ", round(row$rho, 3), "; FDR = ", signif(row$FDR, 3)),
      x = trait_lab,
      y = paste0("Module eigengene: ", mod)
    ) +
    theme_bw(base_size = 12)

  file_tag <- paste0("scatter_", mod, "_", trait)

  ggsave(
    file.path(OUTDIR, "figures", "scatter", paste0(file_tag, ".pdf")),
    p,
    width = 5.5,
    height = 4.5
  )

  ggsave(
    file.path(OUTDIR, "figures", "scatter", paste0(file_tag, ".png")),
    p,
    width = 5.5,
    height = 4.5,
    dpi = DPI
  )
}

###############################################################################
# 12) ADJUSTED REGRESSION MODELS
###############################################################################

regression_results <- list()

for (mod in modules_use) {
  for (outcome in REGRESSION_OUTCOMES) {

    if (!mod %in% names(analysis_df) || !outcome %in% names(analysis_df)) next

    candidate_covars <- c("Age", "Sex_bin", "Education", "Country")
    covars <- intersect(candidate_covars, names(analysis_df))

    form_txt <- paste0(
      outcome,
      " ~ ",
      mod,
      if (length(covars) > 0) paste0(" + ", paste(covars, collapse = " + ")) else ""
    )

    model_df <- analysis_df %>%
      dplyr::select(dplyr::all_of(c(outcome, mod, covars))) %>%
      tidyr::drop_na()

    if (nrow(model_df) < MIN_N_FOR_REG) next

    fit <- tryCatch(
      stats::lm(stats::as.formula(form_txt), data = model_df),
      error = function(e) NULL
    )

    if (is.null(fit)) next

    tidy_fit <- broom::tidy(fit) %>%
      dplyr::filter(term == mod) %>%
      dplyr::mutate(
        Module = mod,
        Outcome = outcome,
        N = nrow(model_df),
        Formula = form_txt,
        adj_r_squared = summary(fit)$adj.r.squared
      )

    regression_results[[paste(mod, outcome, sep = "__")]] <- tidy_fit
  }
}

adjusted_module_models <- dplyr::bind_rows(regression_results)

if (nrow(adjusted_module_models) > 0) {
  adjusted_module_models <- adjusted_module_models %>%
    dplyr::mutate(
      FDR = p.adjust(p.value, method = "BH")
    ) %>%
    dplyr::arrange(FDR, p.value)

  safe_write_csv(
    adjusted_module_models,
    file.path(OUTDIR, "tables", "regression", "adjusted_module_models.csv")
  )
}

###############################################################################
# 13) FINAL INTEGRATED PRIORITIZATION
###############################################################################

main_trait_summary <- module_trait_long %>%
  dplyr::filter(Trait %in% MAIN_TRAITS_FOR_INTEGRATION) %>%
  dplyr::group_by(Module) %>%
  dplyr::summarise(
    n_trait_fdr_0_05 = sum(!is.na(FDR) & FDR < 0.05),
    min_trait_FDR = suppressWarnings(min(FDR, na.rm = TRUE)),
    max_abs_trait_rho = suppressWarnings(max(abs(rho), na.rm = TRUE)),
    mean_abs_trait_rho = mean(abs(rho), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    min_trait_FDR = ifelse(is.infinite(min_trait_FDR), NA_real_, min_trait_FDR),
    max_abs_trait_rho = ifelse(is.infinite(max_abs_trait_rho), NA_real_, max_abs_trait_rho)
  )

hub_summary_clean <- hub_summary %>%
  dplyr::mutate(Module = standardize_module_name(Module))

# -------------------------------------------------------------------------
# Enrichment summary from Script 03 v3
# -------------------------------------------------------------------------
# Script 03 v3 exports:
#   GO_BP_terms / KEGG_terms / Reactome_terms              = FDR < 0.05 terms
#   GO_BP_terms_fdr010 / KEGG_terms_fdr010 / Reactome...   = FDR < 0.10 terms
#   GO_BP_terms_total / KEGG_terms_total / Reactome...     = all retrieved terms
#   *_min_FDR and *_top_term                               = diagnostic/interpretive fields
#
# The prioritization score uses FDR < 0.05 terms by default.
# Suggestive and total terms are preserved for reporting/diagnostics.

if (!"Module" %in% names(enrichment_summary)) {
  warning(
    "enrichment_summary_by_module.csv does not contain a Module column. ",
    "Using zero enrichment terms for all modules."
  )

  enrichment_summary_clean <- tibble::tibble(
    Module = modules_use,
    GO_BP_terms = 0,
    KEGG_terms = 0,
    Reactome_terms = 0,
    GO_BP_terms_fdr010 = 0,
    KEGG_terms_fdr010 = 0,
    Reactome_terms_fdr010 = 0,
    GO_BP_terms_total = 0,
    KEGG_terms_total = 0,
    Reactome_terms_total = 0,
    GO_BP_min_FDR = NA_real_,
    KEGG_min_FDR = NA_real_,
    Reactome_min_FDR = NA_real_,
    GO_BP_top_term = NA_character_,
    KEGG_top_term = NA_character_,
    Reactome_top_term = NA_character_,
    total_enrichment_terms = 0,
    total_enrichment_terms_fdr010 = 0,
    total_enrichment_terms_all = 0,
    best_enrichment_FDR = NA_real_,
    top_enrichment_term = NA_character_
  )

} else {

  enrichment_summary_clean <- enrichment_summary %>%
    dplyr::mutate(
      Module = standardize_module_name(Module),

      GO_BP_terms = if ("GO_BP_terms" %in% names(.)) safe_numeric(GO_BP_terms) else 0,
      KEGG_terms = if ("KEGG_terms" %in% names(.)) safe_numeric(KEGG_terms) else 0,
      Reactome_terms = if ("Reactome_terms" %in% names(.)) safe_numeric(Reactome_terms) else 0,

      GO_BP_terms_fdr010 = if ("GO_BP_terms_fdr010" %in% names(.)) safe_numeric(GO_BP_terms_fdr010) else 0,
      KEGG_terms_fdr010 = if ("KEGG_terms_fdr010" %in% names(.)) safe_numeric(KEGG_terms_fdr010) else 0,
      Reactome_terms_fdr010 = if ("Reactome_terms_fdr010" %in% names(.)) safe_numeric(Reactome_terms_fdr010) else 0,

      GO_BP_terms_total = if ("GO_BP_terms_total" %in% names(.)) safe_numeric(GO_BP_terms_total) else GO_BP_terms,
      KEGG_terms_total = if ("KEGG_terms_total" %in% names(.)) safe_numeric(KEGG_terms_total) else KEGG_terms,
      Reactome_terms_total = if ("Reactome_terms_total" %in% names(.)) safe_numeric(Reactome_terms_total) else Reactome_terms,

      GO_BP_min_FDR = if ("GO_BP_min_FDR" %in% names(.)) safe_numeric(GO_BP_min_FDR) else NA_real_,
      KEGG_min_FDR = if ("KEGG_min_FDR" %in% names(.)) safe_numeric(KEGG_min_FDR) else NA_real_,
      Reactome_min_FDR = if ("Reactome_min_FDR" %in% names(.)) safe_numeric(Reactome_min_FDR) else NA_real_,

      GO_BP_top_term = if ("GO_BP_top_term" %in% names(.)) as.character(GO_BP_top_term) else NA_character_,
      KEGG_top_term = if ("KEGG_top_term" %in% names(.)) as.character(KEGG_top_term) else NA_character_,
      Reactome_top_term = if ("Reactome_top_term" %in% names(.)) as.character(Reactome_top_term) else NA_character_,

      total_enrichment_terms = GO_BP_terms + KEGG_terms + Reactome_terms,
      total_enrichment_terms_fdr010 = GO_BP_terms_fdr010 + KEGG_terms_fdr010 + Reactome_terms_fdr010,
      total_enrichment_terms_all = GO_BP_terms_total + KEGG_terms_total + Reactome_terms_total,

      best_enrichment_FDR = pmin(GO_BP_min_FDR, KEGG_min_FDR, Reactome_min_FDR, na.rm = TRUE),
      best_enrichment_FDR = ifelse(is.infinite(best_enrichment_FDR), NA_real_, best_enrichment_FDR),

      top_enrichment_term = dplyr::case_when(
        !is.na(GO_BP_min_FDR) & GO_BP_min_FDR == best_enrichment_FDR ~ GO_BP_top_term,
        !is.na(KEGG_min_FDR) & KEGG_min_FDR == best_enrichment_FDR ~ KEGG_top_term,
        !is.na(Reactome_min_FDR) & Reactome_min_FDR == best_enrichment_FDR ~ Reactome_top_term,
        TRUE ~ NA_character_
      )
    )
}

safe_write_csv(
  enrichment_summary_clean,
  file.path(OUTDIR, "tables", "enrichment_summary_clean_for_prioritization.csv")
)

if (!exists("loco_summary_by_module")) {
  loco_summary_by_module <- tibble::tibble(
    Module = modules_use,
    mean_abs_delta_rho = NA_real_,
    max_abs_delta_rho = NA_real_,
    direction_consistency = NA_real_
  )
}

final_prioritization <- tibble::tibble(Module = modules_use) %>%
  dplyr::left_join(main_trait_summary, by = "Module") %>%
  dplyr::left_join(module_dep_summary, by = "Module") %>%
  dplyr::left_join(hub_summary_clean, by = "Module") %>%
  dplyr::left_join(enrichment_summary_clean, by = "Module") %>%
  dplyr::left_join(loco_summary_by_module, by = "Module") %>%
  dplyr::mutate(
    clinical_score = scales::rescale(replace_na(mean_abs_trait_rho, 0), to = c(0, 1)),
    dep_score = scales::rescale(replace_na(prop_dep_fdr_0_05, 0), to = c(0, 1)),
    hub_score = scales::rescale(replace_na(mean_abs_kME, 0), to = c(0, 1)),
    enrichment_score = scales::rescale(replace_na(total_enrichment_terms, 0), to = c(0, 1)),
    robustness_score = scales::rescale(
      -replace_na(mean_abs_delta_rho, max(mean_abs_delta_rho, na.rm = TRUE)),
      to = c(0, 1)
    ),
    integrated_priority_score = rowMeans(
      dplyr::across(
        dplyr::all_of(c(
          "clinical_score",
          "dep_score",
          "hub_score",
          "enrichment_score",
          "robustness_score"
        ))
      ),
      na.rm = TRUE
    )
  ) %>%
  dplyr::arrange(dplyr::desc(integrated_priority_score))

safe_write_csv(
  final_prioritization,
  file.path(OUTDIR, "tables", "final_module_prioritization_table.csv")
)

p_prio <- ggplot(
  final_prioritization,
  aes(
    x = mean_abs_trait_rho,
    y = total_enrichment_terms,
    size = mean_abs_kME,
    color = Module,
    label = Module
  )
) +
  geom_point(alpha = 0.75) +
  scale_color_manual(values = get_module_colors(final_prioritization$Module)) +
  ggrepel::geom_text_repel(size = 3.5, max.overlaps = Inf) +
  labs(
    title = "Integrated WGCNA module prioritization",
    x = "Mean absolute module-trait rho",
    y = "Total enrichment terms",
    size = "Mean |kME|"
  ) +
  theme_bw(base_size = 12)

ggsave(
  file.path(OUTDIR, "figures", "module_prioritization_scatter.pdf"),
  p_prio,
  width = 7,
  height = 5
)

ggsave(
  file.path(OUTDIR, "figures", "module_prioritization_scatter.png"),
  p_prio,
  width = 7,
  height = 5,
  dpi = DPI
)

###############################################################################
# 14) FINAL SUMMARY AND WORKSPACE
###############################################################################

script4_summary <- tibble::tibble(
  metric = c(
    "base_dir",
    "outdir",
    "wgcna_input_level",
    "n_samples",
    "n_modules",
    "n_traits",
    "n_module_trait_tests",
    "n_module_trait_fdr_0_05",
    "n_modules_prioritized"
  ),
  value = c(
    BASE_DIR,
    OUTDIR,
    "GENE-COLLAPSED, not aptamer-level",
    as.character(nrow(analysis_df)),
    as.character(length(modules_use)),
    as.character(length(traits_available)),
    as.character(nrow(module_trait_long)),
    as.character(sum(module_trait_long$FDR < 0.05, na.rm = TRUE)),
    as.character(nrow(final_prioritization))
  )
)

safe_write_csv(
  script4_summary,
  file.path(OUTDIR, "tables", "script4_final_summary.csv")
)

save(
  analysis_df,
  module_trait_long,
  rho_mat,
  p_mat,
  fdr_mat,
  annot_mat,
  module_dep_summary,
  final_prioritization,
  script4_summary,
  file = file.path(OUTDIR, "workspace", "script4_module_trait_integration_workspace.RData")
)

writeLines(
  capture.output(utils::sessionInfo()),
  con = file.path(OUTDIR, "sessionInfo.txt")
)

cat("\nScript 04 terminado correctamente.\n")
cat("Directorio principal de salida:\n", OUTDIR, "\n")
cat("Input level: GENE-COLLAPSED, not aptamer-level.\n")
cat("\nArchivos clave:\n")
cat("- module_trait_results_long.csv\n")
cat("- figures/module_trait/module_trait_heatmap_fdr.pdf/png\n")
cat("- figures/module_trait/module_trait_heatmap_nature_style.pdf/png\n")
cat("- figures/module_trait/module_trait_heatmap_manuscript_style.pdf/png\n")
cat("- tables/module_dep_summary.csv\n")
cat("- tables/loco/loco_country_summary_by_module.csv\n")
cat("- tables/regression/adjusted_module_models.csv\n")
cat("- tables/enrichment_summary_clean_for_prioritization.csv\n")
cat("- tables/final_module_prioritization_table.csv\n")

###############################################################################
# END
###############################################################################


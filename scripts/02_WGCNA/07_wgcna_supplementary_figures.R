##############################################################################
# 07_wgcna_supplementary_figures_FINAL_v2.R
#
# TITLE:
# WGCNA supplementary figures using the Figure 3 visual identity
# UPDATED for WGCNA_Workflow_april_V2 and Scripts 03–06 current outputs
# (ReDLat proteomic project)
#
# PURPOSE:
# This script generates the final WGCNA supplementary figure set using the same
# publication style used in Script 06 (Figure 3 main WGCNA figure).
#
# INPUTS:
# - Script 02 outputs:
#   results/02_WGCNA_core_collapsed_genes/tables/module_counts.csv
#   results/02_WGCNA_core_collapsed_genes/soft_threshold/soft_threshold_scan.csv
# - Script 03 outputs:
#   results/03_WGCNA_module_biology_hubs_enrichment/tables/full_kME_assigned_module_long.csv
#   results/03_WGCNA_module_biology_hubs_enrichment/tables/module_hub_summary.csv
#   results/03_WGCNA_module_biology_hubs_enrichment/tables/enrichment_summary_by_module.csv
# - Script 04 outputs:
#   results/04_WGCNA_module_trait_clinical_integration/module_trait_results_long.csv
#   results/04_WGCNA_module_trait_clinical_integration/tables/module_dep_summary.csv
# - Script 05 outputs:
#   results/05_WGCNA_country_site_downsampling_sensitivity/tables/loco_country_summary_by_module.csv
#   results/05_WGCNA_country_site_downsampling_sensitivity/tables/loco_country_summary_by_module_trait.csv
#   results/05_WGCNA_country_site_downsampling_sensitivity/tables/loso_site_summary_by_module.csv
#   results/05_WGCNA_country_site_downsampling_sensitivity/tables/balanced_downsampling_summary_by_module_trait.csv
#
# MAIN OUTPUTS:
# - Supplementary_Fig_4_complete_module_trait_heatmap.pdf/png
# - Supplementary_Fig_5_soft_thresholding.pdf/png
# - Supplementary_Fig_6_module_size_DEP_burden.pdf/png
# - Supplementary_Fig_7_hub_proteins_by_module.pdf/png      [if kME/hub files exist]
# - Supplementary_Fig_8_enrichment_summary.pdf/png
# - Supplementary_Fig_9_LOCO_country_robustness.pdf/png
# - Supplementary_Fig_10_downsampling_LOSO_robustness.pdf/png
# - supplementary_figure_inventory.csv
#
# AUTHOR:
# Matías Pizarro + ChatGPT support
##############################################################################

rm(list = ls())

##############################################################################
# 1) PACKAGES
##############################################################################

cran_pkgs <- c(
  "dplyr", "readr", "tidyr", "stringr", "forcats", "purrr",
  "ggplot2", "patchwork", "scales", "tibble"
)

cran_missing <- cran_pkgs[!sapply(cran_pkgs, requireNamespace, quietly = TRUE)]
if (length(cran_missing) > 0) install.packages(cran_missing)

invisible(lapply(cran_pkgs, library, character.only = TRUE))

options(stringsAsFactors = FALSE)
options(error = traceback)

##############################################################################
# 2) FIGURE STYLE — MATCH SCRIPT 06 / FIGURE 3
#    Word-optimized version: larger text for manuscript insertion
##############################################################################

BASE_FONT_FAMILY  <- "sans"
BASE_FONT_SIZE    <- 20
TITLE_SIZE        <- 22
SUBTITLE_SIZE     <- 18
AXIS_TITLE_SIZE   <- 19
AXIS_TEXT_SIZE    <- 17
STRIP_TEXT_SIZE   <- 18
LEGEND_TITLE_SIZE <- 18
LEGEND_TEXT_SIZE  <- 17

LINE_WIDTH_BASE <- 0.85
POINT_SIZE_BASE <- 3.0
PANEL_TAG_SIZE  <- 24

COL_UP      <- "#f46d43"
COL_DOWN    <- "#4682b4"
COL_NEUTRAL <- "grey70"
COL_NS      <- "grey78"
COL_TEXT    <- "black"
COL_BORDER  <- "grey85"
COL_GRID    <- "grey93"

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

TRAIT_LABELS <- c(
  SampleGroup_bin = "Diagnosis",
  cdr_global      = "CDR global",
  cdr_boxscore    = "CDR-SB",
  mmse_total      = "MMSE",
  udsfaq_total    = "PFAQ",
  NPI             = "NPI-Q",
  Mini_SEA        = "Mini-SEA",
  T_ADLQ          = "T-ADLQ",
  p_tau181        = "p-tau181",
  p_tau217        = "p-tau217",
  NfL             = "NfL",
  ratio_AB42_40   = "Aβ42/40",
  Age             = "Age",
  Sex_bin         = "Sex",
  Education       = "Education",
  APOE4_carrier   = "APOE ε4",
  Country_numeric = "Country"
)

TRAIT_ORDER <- c(
  "SampleGroup_bin", "cdr_global", "cdr_boxscore", "mmse_total",
  "udsfaq_total", "NPI", "Mini_SEA", "T_ADLQ",
  "p_tau181", "p_tau217", "NfL", "ratio_AB42_40",
  "Age", "Sex_bin", "Education", "APOE4_carrier", "Country_numeric"
)

theme_pub <- function(base_size = BASE_FONT_SIZE, base_family = BASE_FONT_FAMILY) {
  ggplot2::theme_bw(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(linewidth = LINE_WIDTH_BASE, color = COL_TEXT),
      axis.line = ggplot2::element_line(linewidth = LINE_WIDTH_BASE, color = COL_TEXT),
      axis.ticks = ggplot2::element_line(linewidth = LINE_WIDTH_BASE, color = COL_TEXT),
      axis.title = ggplot2::element_text(size = AXIS_TITLE_SIZE, color = COL_TEXT),
      axis.text = ggplot2::element_text(size = AXIS_TEXT_SIZE, color = COL_TEXT),
      plot.title = ggplot2::element_text(size = TITLE_SIZE, face = "bold", hjust = 0.5, color = COL_TEXT),
      plot.subtitle = ggplot2::element_text(size = SUBTITLE_SIZE, hjust = 0.5, color = "grey35"),
      strip.background = ggplot2::element_rect(fill = "white", color = COL_TEXT, linewidth = LINE_WIDTH_BASE),
      strip.text = ggplot2::element_text(size = STRIP_TEXT_SIZE, face = "bold", color = COL_TEXT),
      legend.title = ggplot2::element_text(size = LEGEND_TITLE_SIZE, color = COL_TEXT),
      legend.text = ggplot2::element_text(size = LEGEND_TEXT_SIZE, color = COL_TEXT),
      legend.key = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(10, 10, 10, 10)
    )
}

theme_card <- function(base_size = 14, base_family = BASE_FONT_FAMILY) {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = COL_BORDER, fill = NA, linewidth = 0.6),
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = 14, color = COL_TEXT),
      axis.text = ggplot2::element_text(size = 13, color = COL_TEXT),
      plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0, color = COL_TEXT),
      plot.subtitle = ggplot2::element_text(size = 13, hjust = 0, color = "grey35"),
      plot.margin = ggplot2::margin(8, 8, 8, 8),
      legend.position = "right"
    )
}

save_pub_pdf <- function(plot_obj, filename, dims = c(8, 6)) {
  ggplot2::ggsave(
    filename = filename,
    plot = plot_obj,
    width = dims[1],
    height = dims[2],
    units = "in",
    device = grDevices::cairo_pdf,
    dpi = 300,
    bg = "white",
    limitsize = FALSE
  )
}

save_pub_png <- function(plot_obj, filename, dims = c(8, 6)) {
  ggplot2::ggsave(
    filename = filename,
    plot = plot_obj,
    width = dims[1],
    height = dims[2],
    units = "in",
    dpi = 600,
    bg = "white",
    limitsize = FALSE
  )
}

save_pub_both <- function(plot_obj, filename_base, dims = c(8, 6)) {
  save_pub_pdf(plot_obj, paste0(filename_base, ".pdf"), dims = dims)
  save_pub_png(plot_obj, paste0(filename_base, ".png"), dims = dims)
}

##############################################################################
# 3) PATHS
##############################################################################

BASE_DIR <- "C:/Users/mnpiz/Desktop/WGCNA_Workflow_april_V2"

SCRIPT2_DIR <- file.path(BASE_DIR, "results", "02_WGCNA_core_collapsed_genes")
SCRIPT3_DIR <- file.path(BASE_DIR, "results", "03_WGCNA_module_biology_hubs_enrichment")
SCRIPT4_DIR <- file.path(BASE_DIR, "results", "04_WGCNA_module_trait_clinical_integration")
SCRIPT5_DIR <- file.path(BASE_DIR, "results", "05_WGCNA_country_site_downsampling_sensitivity")

OUTDIR <- file.path(BASE_DIR, "results", "07_wgcna_supplementary_figures")
OUT_FIG <- file.path(OUTDIR, "figures")
OUT_TAB <- file.path(OUTDIR, "tables")
OUT_PANEL <- file.path(OUTDIR, "separate_panels")

invisible(lapply(c(OUTDIR, OUT_FIG, OUT_TAB, OUT_PANEL), dir.create, recursive = TRUE, showWarnings = FALSE))

MODULE_COUNTS_FILE <- file.path(SCRIPT2_DIR, "tables", "module_counts.csv")

# Different versions of Script 02 may save the soft-threshold table with different names.
SFT_CANDIDATES <- c(
  file.path(SCRIPT2_DIR, "soft_threshold", "soft_threshold_scan.csv"),
  file.path(SCRIPT2_DIR, "soft_threshold", "soft_threshold_fit_indices.csv"),
  file.path(SCRIPT2_DIR, "soft_threshold", "soft_threshold_diagnostics.csv")
)

KME_FILE <- file.path(SCRIPT3_DIR, "tables", "full_kME_assigned_module_long.csv")
HUB_SUMMARY_FILE <- file.path(SCRIPT3_DIR, "tables", "module_hub_summary.csv")
ENRICHMENT_SUMMARY_FILE <- file.path(SCRIPT3_DIR, "tables", "enrichment_summary_by_module.csv")

MODULE_TRAIT_FILE <- file.path(SCRIPT4_DIR, "tables", "module_trait_results_long.csv")
MODULE_DEP_FILE <- file.path(SCRIPT4_DIR, "tables", "module_dep_summary.csv")
PRIORITIZATION_FILE <- file.path(SCRIPT4_DIR, "tables", "final_module_prioritization_table.csv")

LOCO_MODULE_FILE <- file.path(SCRIPT5_DIR, "tables", "loco", "loco_country_summary_by_module.csv")
LOCO_TRAIT_FILE <- file.path(SCRIPT5_DIR, "tables", "loco", "loco_country_summary_by_module_trait.csv")
LOSO_MODULE_FILE <- file.path(SCRIPT5_DIR, "tables", "loso", "loso_site_summary_by_module.csv")
DOWNSAMPLE_FILE <- file.path(SCRIPT5_DIR, "tables", "downsampling", "balanced_downsampling_summary_by_module_trait.csv")

##############################################################################
# 4) HELPERS
##############################################################################

safe_read <- function(file, required = TRUE) {
  if (!file.exists(file)) {
    if (required) stop("No existe archivo requerido:\n", file)
    return(NULL)
  }
  readr::read_csv(file, show_col_types = FALSE)
}

safe_col <- function(df, candidates) {
  if (is.null(df)) return(NA_character_)
  out <- candidates[candidates %in% names(df)]
  if (length(out) == 0) return(NA_character_)
  out[1]
}

standardize_module <- function(x) {
  x <- as.character(x)
  x <- gsub("^ME", "", x)
  x
}

module_fill_scale <- function() {
  ggplot2::scale_fill_manual(values = MODULE_COLORS, na.value = "grey80")
}

module_color_scale <- function() {
  ggplot2::scale_color_manual(values = MODULE_COLORS, na.value = "grey60")
}

label_trait <- function(x) {
  out <- TRAIT_LABELS[x]
  out[is.na(out)] <- x[is.na(out)]
  unname(out)
}

make_module_display_labels <- function(modules, show_internal_name = FALSE) {
  modules <- as.character(modules)
  if (show_internal_name) {
    labels <- paste0("M", seq_along(modules), " / ", modules)
  } else {
    labels <- modules
  }
  stats::setNames(labels, modules)
}

plot_tag_theme <- function() {
  ggplot2::theme(plot.tag = ggplot2::element_text(size = PANEL_TAG_SIZE, face = "bold"))
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

blank_panel <- function(title, subtitle = NULL) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0.5, y = 0.58, label = title, fontface = "bold", size = 5) +
    ggplot2::annotate("text", x = 0.5, y = 0.42, label = ifelse(is.null(subtitle), "", subtitle), size = 4) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", color = NA))
}

##############################################################################
# 5) LOAD DATA
##############################################################################

module_counts <- safe_read(MODULE_COUNTS_FILE)
module_trait  <- safe_read(MODULE_TRAIT_FILE)
module_dep    <- safe_read(MODULE_DEP_FILE)
enrich_summary <- safe_read(ENRICHMENT_SUMMARY_FILE)

kme_tbl <- safe_read(KME_FILE, required = FALSE)
hub_summary <- safe_read(HUB_SUMMARY_FILE, required = FALSE)
prioritization <- safe_read(PRIORITIZATION_FILE, required = FALSE)

loco_module <- safe_read(LOCO_MODULE_FILE, required = FALSE)
loco_trait <- safe_read(LOCO_TRAIT_FILE, required = FALSE)
loso_module <- safe_read(LOSO_MODULE_FILE, required = FALSE)
downsample <- safe_read(DOWNSAMPLE_FILE, required = FALSE)

sft_tbl <- NULL
sft_file_used <- NA_character_
for (candidate in SFT_CANDIDATES) {
  if (file.exists(candidate)) {
    sft_tbl <- safe_read(candidate, required = FALSE)
    sft_file_used <- candidate
    break
  }
}

objs_to_standardize <- c(
  "module_counts", "module_trait", "module_dep", "enrich_summary",
  "kme_tbl", "hub_summary", "prioritization",
  "loco_module", "loco_trait", "loso_module", "downsample"
)

for (obj in objs_to_standardize) {
  x <- get(obj)
  if (!is.null(x) && "Module" %in% names(x)) {
    assign(obj, x %>% dplyr::mutate(Module = standardize_module(Module)))
  }
}

# Compatibility with current Script 04 column names
if (!"rho" %in% names(module_trait)) {
  rho_col <- safe_col(module_trait, c("correlation", "Correlation", "r", "estimate", "Estimate"))
  if (!is.na(rho_col)) module_trait <- module_trait %>% dplyr::rename(rho = dplyr::all_of(rho_col))
}

if (!"FDR" %in% names(module_trait)) {
  fdr_col <- safe_col(module_trait, c("p.adjust", "qvalue", "padj", "adj_p", "adj.P.Val"))
  if (!is.na(fdr_col)) module_trait <- module_trait %>% dplyr::rename(FDR = dplyr::all_of(fdr_col))
}

module_counts <- module_counts %>%
  dplyr::mutate(
    Module = standardize_module(Module),
    N_genes = suppressWarnings(as.numeric(N_genes))
  )

CORE_MODULES <- c("yellow", "brown", "black", "blue", "purple", "pink")
CORE_MODULES <- CORE_MODULES[CORE_MODULES %in% unique(module_counts$Module)]

##############################################################################
# 6) SUPPLEMENTARY FIG. 4 — COMPLETE MODULE–TRAIT HEATMAP
##############################################################################

required_mt_cols <- c("Module", "Trait", "rho", "FDR")
missing_mt_cols <- setdiff(required_mt_cols, names(module_trait))
if (length(missing_mt_cols) > 0) {
  stop("Faltan columnas en module_trait_results_long.csv: ", paste(missing_mt_cols, collapse = ", "))
}

module_order <- module_counts %>%
  dplyr::filter(Module != "grey") %>%
  dplyr::arrange(dplyr::desc(N_genes)) %>%
  dplyr::pull(Module)

trait_order <- TRAIT_ORDER[TRAIT_ORDER %in% unique(module_trait$Trait)]
extra_traits <- setdiff(unique(module_trait$Trait), trait_order)
trait_order <- c(trait_order, extra_traits)

MODULE_DISPLAY_MAP_ALL <- make_module_display_labels(
  modules = module_order,
  show_internal_name = FALSE
)

mt_plot <- module_trait %>%
  dplyr::filter(Module %in% module_order, Trait %in% trait_order) %>%
  dplyr::mutate(
    rho = suppressWarnings(as.numeric(rho)),
    FDR = suppressWarnings(as.numeric(FDR)),
    stars = add_sig_stars(FDR),
    # Match Figure 3: show only BH-FDR stars inside cells, not rho values.
    label = stars,
    Trait_label = label_trait(Trait),
    Module_display = unname(MODULE_DISPLAY_MAP_ALL[Module]),
    Module_display = factor(Module_display, levels = rev(unname(MODULE_DISPLAY_MAP_ALL[module_order]))),
    Trait_label = factor(Trait_label, levels = label_trait(trait_order))
  ) %>%
  dplyr::filter(!is.na(Module_display), !is.na(Trait_label))

rho_lim <- max(abs(mt_plot$rho), na.rm = TRUE)
if (!is.finite(rho_lim)) rho_lim <- 0.5
rho_lim <- max(rho_lim, 0.35)

side_df_s4 <- tibble::tibble(
  Module = module_order,
  Module_display = factor(
    unname(MODULE_DISPLAY_MAP_ALL[module_order]),
    levels = rev(unname(MODULE_DISPLAY_MAP_ALL[module_order]))
  ),
  x = 1
)

p_s4_side <- ggplot(side_df_s4, aes(x = x, y = Module_display, fill = Module)) +
  geom_tile(width = 1, height = 0.98, color = NA) +
  scale_fill_manual(values = MODULE_COLORS, guide = "none", na.value = "grey80") +
  scale_x_continuous(expand = c(0, 0)) +
  theme_void() +
  theme(
    plot.margin = margin(34, 0, 58, 4),
    plot.background = element_rect(fill = "white", color = NA)
  )

p_s4_heat <- ggplot(mt_plot, aes(x = Trait_label, y = Module_display, fill = rho)) +
  geom_tile(color = "grey90", linewidth = 0.55) +
  geom_text(aes(label = label), size = 7.6, color = "black", fontface = "bold") +
  scale_fill_gradient2(
    low = "#6FA8DC",
    mid = "white",
    high = "#E06666",
    midpoint = 0,
    limits = c(-rho_lim, rho_lim),
    name = "Spearman\nrho"
  ) +
  labs(
    title = "Complete WGCNA module–trait association matrix",
    subtitle = "BH-FDR corrected significance overlaid as stars",
    x = NULL,
    y = NULL
  ) +
  theme_pub(base_size = 20) +
  theme(
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_text(face = "bold", color = COL_TEXT, size = 17),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = COL_TEXT, size = 17),
    legend.position = "right",
    plot.margin = margin(10, 12, 10, 0)
  )

p_s4_heatmap <- p_s4_side + p_s4_heat +
  patchwork::plot_layout(widths = c(0.035, 0.965))

save_pub_both(
  p_s4_heatmap,
  file.path(OUT_FIG, "Supplementary_Fig_4_complete_module_trait_heatmap"),
  dims = c(max(18, 1.25 * length(trait_order)), max(10, 0.85 * length(module_order) + 4.0))
)

readr::write_csv(mt_plot, file.path(OUT_TAB, "Supplementary_Fig_4_complete_module_trait_heatmap_source.csv"))

##############################################################################
# 7) SUPPLEMENTARY FIG. 5 — SOFT-THRESHOLDING DIAGNOSTICS
##############################################################################

if (!is.null(sft_tbl)) {

  names(sft_tbl) <- make.names(names(sft_tbl))

  power_col <- safe_col(sft_tbl, c("Power", "power"))
  r2_col <- safe_col(sft_tbl, c("SFT.R.sq", "SFT.R.sq.", "SFT_R_sq", "SFT.R.sq.."))
  mean_k_col <- safe_col(sft_tbl, c("mean.k.", "mean.k", "mean_k", "MeanConnectivity", "mean.k.."))

  if (!is.na(power_col) && !is.na(r2_col) && !is.na(mean_k_col)) {

    sft_tbl <- sft_tbl %>%
      dplyr::mutate(
        Power_plot = suppressWarnings(as.numeric(.data[[power_col]])),
        SFT_R2_plot = suppressWarnings(as.numeric(.data[[r2_col]])),
        Mean_k_plot = suppressWarnings(as.numeric(.data[[mean_k_col]]))
      )

    p_sft_r2 <- ggplot(sft_tbl, aes(x = Power_plot, y = SFT_R2_plot)) +
      geom_line(linewidth = 0.7, color = "black") +
      geom_point(size = 3.2, color = "black") +
      geom_hline(yintercept = 0.90, linetype = "dashed", color = "grey40") +
      geom_vline(xintercept = 4, linetype = "dashed", color = COL_UP) +
      annotate("text", x = 4, y = max(sft_tbl$SFT_R2_plot, na.rm = TRUE),
               label = "β = 4", hjust = -0.1, vjust = 1.1, size = 5.0) +
      labs(
        title = "Scale-free topology fit",
        x = "Soft-thresholding power",
        y = "Signed R²"
      ) +
      theme_pub()

    p_sft_k <- ggplot(sft_tbl, aes(x = Power_plot, y = Mean_k_plot)) +
      geom_line(linewidth = 0.7, color = "black") +
      geom_point(size = 3.2, color = "black") +
      geom_vline(xintercept = 4, linetype = "dashed", color = COL_UP) +
      labs(
        title = "Mean connectivity",
        x = "Soft-thresholding power",
        y = "Mean connectivity"
      ) +
      theme_pub()

    fig_s5 <- p_sft_r2 + p_sft_k +
      patchwork::plot_annotation(tag_levels = "a", theme = plot_tag_theme())

    save_pub_both(fig_s5, file.path(OUT_FIG, "Supplementary_Fig_5_soft_thresholding"), dims = c(13.5, 6.5))

  } else {
    message("No se reconocieron columnas de soft-thresholding en: ", sft_file_used)
  }

} else {
  message("No se encontró tabla de soft-thresholding. Usa el PDF original de Script 02 si es necesario.")
}

##############################################################################
# 8) SUPPLEMENTARY FIG. 6 — MODULE SIZE + DEP BURDEN
##############################################################################

dep_prop_col <- safe_col(module_dep, c("prop_dep_fdr_0_05", "Prop_FDR005", "prop_dep_fdr05", "Prop_DEP_FDR005"))
dep_n_col <- safe_col(module_dep, c("n_dep_fdr_0_05", "N_DEPs_FDR005", "n_dep_fdr05", "N_DEP_FDR005"))
avg_logfc_col <- safe_col(module_dep, c("mean_abs_logFC", "Avg_abs_logFC", "avg_abs_logFC"))

if (is.na(dep_prop_col)) {
  stop("No encuentro columna de proporción DEP FDR<0.05. Revisa module_dep_summary.csv.")
}

p_s6a_size <- module_counts %>%
  dplyr::filter(Module != "grey") %>%
  dplyr::mutate(Module = factor(Module, levels = Module[order(N_genes)])) %>%
  ggplot(aes(x = Module, y = N_genes, fill = as.character(Module))) +
  geom_col(width = 0.75) +
  coord_flip() +
  module_fill_scale() +
  labs(title = "Module size", x = NULL, y = "Number of proteins") +
  theme_pub() +
  theme(legend.position = "none")

p_s6b_dep_prop <- module_dep %>%
  dplyr::filter(Module != "grey") %>%
  dplyr::mutate(Module = factor(Module, levels = Module[order(.data[[dep_prop_col]])])) %>%
  ggplot(aes(x = Module, y = .data[[dep_prop_col]], fill = as.character(Module))) +
  geom_col(width = 0.75) +
  coord_flip() +
  module_fill_scale() +
  labs(title = "DEP density by module", x = NULL, y = "Proportion DEP at FDR < 0.05") +
  theme_pub() +
  theme(legend.position = "none")

if (!is.na(dep_n_col)) {
  p_s6c_dep_n <- module_dep %>%
    dplyr::filter(Module != "grey") %>%
    dplyr::mutate(Module = factor(Module, levels = Module[order(.data[[dep_n_col]])])) %>%
    ggplot(aes(x = Module, y = .data[[dep_n_col]], fill = as.character(Module))) +
    geom_col(width = 0.75) +
    coord_flip() +
    module_fill_scale() +
    labs(title = "Number of DEPs by module", x = NULL, y = "Number of DEPs at FDR < 0.05") +
    theme_pub() +
    theme(legend.position = "none")
} else {
  p_s6c_dep_n <- blank_panel("N_DEP column not available")
}

if (!is.na(avg_logfc_col)) {
  p_s6d_logfc <- module_dep %>%
    dplyr::filter(Module != "grey") %>%
    dplyr::mutate(Module = factor(Module, levels = Module[order(.data[[avg_logfc_col]])])) %>%
    ggplot(aes(x = Module, y = .data[[avg_logfc_col]], fill = as.character(Module))) +
    geom_col(width = 0.75) +
    coord_flip() +
    module_fill_scale() +
    labs(title = "Average absolute logFC", x = NULL, y = "Mean absolute logFC") +
    theme_pub() +
    theme(legend.position = "none")
} else {
  p_s6d_logfc <- blank_panel("Avg_abs_logFC column not available")
}

fig_s6 <- (p_s6a_size + p_s6b_dep_prop) / (p_s6c_dep_n + p_s6d_logfc) +
  patchwork::plot_annotation(tag_levels = "a", theme = plot_tag_theme())

save_pub_both(fig_s6, file.path(OUT_FIG, "Supplementary_Fig_6_module_size_DEP_burden"), dims = c(14.5, 11))

##############################################################################
##############################################################################
# 9) SUPPLEMENTARY FIG. 7 — HUB PROTEINS BY MODULE + kME STRUCTURE
##############################################################################

if (!is.null(kme_tbl) && "abs_kME" %in% names(kme_tbl)) {

  kme_plot <- kme_tbl %>%
    dplyr::filter(!is.na(Module), Module != "grey") %>%
    dplyr::mutate(
      abs_kME = suppressWarnings(as.numeric(abs_kME)),
      Module = standardize_module(Module),
      Module = factor(Module, levels = module_counts$Module[module_counts$Module != "grey"])
    )

  # Panel a: global module-membership distribution
  p_s7a <- kme_plot %>%
    ggplot(aes(x = Module, y = abs_kME, fill = as.character(Module))) +
    geom_boxplot(width = 0.65, outlier.size = 1.1, alpha = 0.9, linewidth = 0.7) +
    module_fill_scale() +
    labs(
      title = "Module membership distribution",
      x = NULL,
      y = "Absolute module membership (|kME|)"
    ) +
    theme_pub() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1, size = 17, face = "bold")
    )

  # Panel b: high-membership hub burden
  if (!is.null(hub_summary)) {
    hub_n_col <- safe_col(hub_summary, c("n_hubs_abs_kME_ge_0_8", "n_hubs_abs_kME_ge_0_7"))

    if (!is.na(hub_n_col)) {
      p_s7b <- hub_summary %>%
        dplyr::filter(!is.na(Module), Module != "grey") %>%
        dplyr::mutate(
          Module = standardize_module(Module),
          Module = factor(Module, levels = Module[order(.data[[hub_n_col]])])
        ) %>%
        ggplot(aes(x = Module, y = .data[[hub_n_col]], fill = as.character(Module))) +
        geom_col(width = 0.75, color = "black", linewidth = 0.35) +
        coord_flip() +
        module_fill_scale() +
        labs(
          title = "High-membership hub burden",
          x = NULL,
          y = ifelse(hub_n_col == "n_hubs_abs_kME_ge_0_8", "Number of hubs with |kME| ≥ 0.80", "Number of hubs with |kME| ≥ 0.70")
        ) +
        theme_pub() +
        theme(legend.position = "none")
    } else {
      p_s7b <- blank_panel("Hub-count column not available")
    }
  } else {
    p_s7b <- blank_panel("Hub summary not available")
  }

  # Panel c: explicit top hub proteins by prioritized module
  HUB_MODULES <- c("yellow", "brown", "black", "blue", "purple", "pink")
  HUB_MODULES <- HUB_MODULES[HUB_MODULES %in% unique(as.character(kme_plot$Module))]
  if (length(HUB_MODULES) == 0) {
    HUB_MODULES <- unique(as.character(kme_plot$Module))[seq_len(min(6, length(unique(as.character(kme_plot$Module)))))]
  }

  hub_label_col <- safe_col(kme_plot, c("Protein_Display", "Protein_Name", "EntrezGeneSymbol", "TargetFullName", "Target"))
  if (is.na(hub_label_col)) hub_label_col <- "EntrezGeneSymbol"

  top_hub_proteins <- kme_plot %>%
    dplyr::filter(as.character(Module) %in% HUB_MODULES) %>%
    dplyr::mutate(
      Hub_label_raw = as.character(.data[[hub_label_col]]),
      Hub_label_raw = dplyr::if_else(is.na(Hub_label_raw) | Hub_label_raw == "", as.character(EntrezGeneSymbol), Hub_label_raw),
      Hub_label_raw = stringr::str_replace_all(Hub_label_raw, "_", " "),
      Hub_label = stringr::str_trunc(Hub_label_raw, width = 34),
      Module_chr = as.character(Module)
    ) %>%
    dplyr::filter(!is.na(abs_kME), !is.na(Hub_label), Hub_label != "") %>%
    dplyr::group_by(Module_chr) %>%
    dplyr::arrange(dplyr::desc(abs_kME), .by_group = TRUE) %>%
    dplyr::slice_head(n = 10) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      Hub_label = forcats::fct_reorder(Hub_label, abs_kME),
      Module_chr = factor(Module_chr, levels = HUB_MODULES)
    )

  readr::write_csv(
    top_hub_proteins,
    file.path(OUT_TAB, "Supplementary_Fig_7_top_hub_proteins_source.csv")
  )

  if (nrow(top_hub_proteins) > 0) {
    p_s7c <- ggplot(top_hub_proteins, aes(x = Hub_label, y = abs_kME, fill = Module_chr)) +
      geom_col(width = 0.72, color = "black", linewidth = 0.25) +
      coord_flip() +
      facet_wrap(~ Module_chr, scales = "free_y", ncol = 2) +
      module_fill_scale() +
      scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.04))) +
      labs(
        title = "Top hub proteins by prioritized WGCNA module",
        subtitle = "Proteins ranked by absolute module membership (|kME|)",
        x = NULL,
        y = "|kME|"
      ) +
      theme_pub(base_size = 19) +
      theme(
        legend.position = "none",
        axis.text.y = element_text(size = 15, color = "black", face = "bold"),
        axis.text.x = element_text(size = 15, color = "black"),
        strip.text = element_text(size = 18, face = "bold"),
        plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 17, hjust = 0.5, color = "grey35")
      )
  } else {
    p_s7c <- blank_panel("Top hub proteins not available")
  }

  fig_s7_overview <- p_s7a + p_s7b +
    patchwork::plot_layout(widths = c(1.25, 1)) +
    patchwork::plot_annotation(tag_levels = "a", theme = plot_tag_theme())

  fig_s7 <- fig_s7_overview / p_s7c +
    patchwork::plot_layout(heights = c(0.9, 1.8)) +
    patchwork::plot_annotation(tag_levels = "a", theme = plot_tag_theme())

  save_pub_both(
    fig_s7,
    file.path(OUT_FIG, "Supplementary_Fig_7_hub_proteins_by_module"),
    dims = c(17.5, 18.5)
  )
}

##############################################################################
# 10) SUPPLEMENTARY FIG. 8 — FUNCTIONAL ENRICHMENT SUMMARY
##############################################################################

needed_enrich_cols <- c("GO_BP_terms", "KEGG_terms", "Reactome_terms")
missing_enrich_cols <- setdiff(needed_enrich_cols, names(enrich_summary))

if (length(missing_enrich_cols) == 0) {

  enrich_long <- enrich_summary %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(needed_enrich_cols),
      names_to = "Library",
      values_to = "N_terms"
    ) %>%
    dplyr::mutate(
      Library = dplyr::recode(
        Library,
        GO_BP_terms = "GO BP",
        KEGG_terms = "KEGG",
        Reactome_terms = "Reactome"
      ),
      Module = factor(Module, levels = module_counts$Module[order(module_counts$N_genes)])
    )

  p_s8a <- enrich_long %>%
    ggplot(aes(x = Module, y = N_terms, fill = Library)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7) +
    coord_flip() +
    labs(title = "Functional enrichment terms by module", x = NULL, y = "Number of enriched terms") +
    scale_fill_manual(values = c("GO BP" = "#8DA0CB", "KEGG" = "#FC8D62", "Reactome" = "#66C2A5")) +
    theme_pub()

  p_s8b <- enrich_long %>%
    dplyr::group_by(Module) %>%
    dplyr::summarise(Total_terms = sum(N_terms, na.rm = TRUE), .groups = "drop") %>%
    ggplot(aes(x = reorder(Module, Total_terms), y = Total_terms, fill = as.character(Module))) +
    geom_col(width = 0.75) +
    coord_flip() +
    module_fill_scale() +
    labs(title = "Total enrichment burden", x = NULL, y = "Total enriched terms") +
    theme_pub() +
    theme(legend.position = "none")

  fig_s8 <- p_s8a + p_s8b +
    patchwork::plot_layout(widths = c(1.25, 1)) +
    patchwork::plot_annotation(tag_levels = "a", theme = plot_tag_theme())

  save_pub_both(fig_s8, file.path(OUT_FIG, "Supplementary_Fig_8_enrichment_summary"), dims = c(14.5, 7.5))

} else {
  message("Faltan columnas de enrichment_summary_by_module.csv: ", paste(missing_enrich_cols, collapse = ", "))
}

##############################################################################
# 11) SUPPLEMENTARY FIG. 9 — LOCO ROBUSTNESS
##############################################################################

if (!is.null(loco_module)) {

  same_dir_col_loco <- safe_col(loco_module, c("direction_consistency", "prop_same_direction", "same_direction"))
  mean_delta_col_loco <- safe_col(loco_module, c("mean_abs_delta_rho"))
  max_delta_col_loco <- safe_col(loco_module, c("max_abs_delta_rho"))

  if (!is.na(mean_delta_col_loco) && !is.na(max_delta_col_loco)) {

    p_s9a <- loco_module %>%
      dplyr::mutate(
        mean_abs_delta_rho = suppressWarnings(as.numeric(.data[[mean_delta_col_loco]])),
        Module = factor(Module, levels = Module[order(mean_abs_delta_rho)])
      ) %>%
      ggplot(aes(x = Module, y = mean_abs_delta_rho, fill = as.character(Module))) +
      geom_col(width = 0.75) +
      coord_flip() +
      module_fill_scale() +
      labs(title = "Leave-one-country-out robustness", x = NULL, y = "Mean absolute delta rho") +
      theme_pub() +
      theme(legend.position = "none")

    p_s9b <- loco_module %>%
      dplyr::mutate(
        max_abs_delta_rho = suppressWarnings(as.numeric(.data[[max_delta_col_loco]])),
        Module = factor(Module, levels = Module[order(max_abs_delta_rho)])
      ) %>%
      ggplot(aes(x = Module, y = max_abs_delta_rho, fill = as.character(Module))) +
      geom_col(width = 0.75) +
      coord_flip() +
      module_fill_scale() +
      labs(title = "Maximum country-exclusion change", x = NULL, y = "Maximum absolute delta rho") +
      theme_pub() +
      theme(legend.position = "none")

    if (!is.na(same_dir_col_loco)) {
      p_s9c <- loco_module %>%
        dplyr::mutate(
          direction_consistency = suppressWarnings(as.numeric(.data[[same_dir_col_loco]])),
          Module = factor(Module, levels = Module[order(direction_consistency)])
        ) %>%
        ggplot(aes(x = Module, y = direction_consistency, fill = as.character(Module))) +
        geom_col(width = 0.75) +
        coord_flip() +
        module_fill_scale() +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
        labs(title = "Direction preservation", x = NULL, y = "Same direction") +
        theme_pub() +
        theme(legend.position = "none")
    } else {
      p_s9c <- blank_panel("Direction consistency column not available")
    }

    fig_s9 <- (p_s9a + p_s9b + p_s9c) +
      patchwork::plot_annotation(tag_levels = "a", theme = plot_tag_theme())

    save_pub_both(fig_s9, file.path(OUT_FIG, "Supplementary_Fig_9_LOCO_country_robustness"), dims = c(17, 6.5))
  }
}

##############################################################################
# 12) SUPPLEMENTARY FIG. 10 — DOWNSAMPLING AND OPTIONAL LOSO
##############################################################################

if (!is.null(downsample)) {

  down_delta_col <- safe_col(downsample, c("mean_abs_delta_rho", "delta_mean_vs_full", "delta_rho", "abs_delta_mean_vs_full"))
  same_dir_col <- safe_col(downsample, c("direction_consistency", "prop_same_direction_as_full", "prop_same_direction"))
  max_delta_col <- safe_col(downsample, c("max_abs_delta_rho"))

  if (!is.na(down_delta_col)) {

    down_plot <- downsample %>%
      dplyr::mutate(
        Trait_label = if ("Trait" %in% names(.)) label_trait(Trait) else NA_character_,
        Module = factor(Module, levels = unique(Module)),
        down_delta_value = suppressWarnings(as.numeric(.data[[down_delta_col]]))
      )

    p_s10a <- down_plot %>%
      ggplot(aes(x = Module, y = down_delta_value, fill = as.character(Module))) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
      geom_boxplot(width = 0.65, outlier.size = 2.0) +
      module_fill_scale() +
      labs(title = "Balanced downsampling stability", x = NULL, y = "Mean absolute delta rho") +
      theme_pub() +
      theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

    if (!is.na(same_dir_col)) {
      p_s10b <- down_plot %>%
        dplyr::mutate(same_dir_value = suppressWarnings(as.numeric(.data[[same_dir_col]]))) %>%
        dplyr::group_by(Module) %>%
        dplyr::summarise(mean_same_direction = mean(same_dir_value, na.rm = TRUE), .groups = "drop") %>%
        ggplot(aes(x = reorder(Module, mean_same_direction), y = mean_same_direction, fill = as.character(Module))) +
        geom_col(width = 0.75) +
        coord_flip() +
        module_fill_scale() +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
        labs(title = "Direction preservation after downsampling", x = NULL, y = "Same direction vs full model") +
        theme_pub() +
        theme(legend.position = "none")
    } else {
      p_s10b <- blank_panel("Direction column not available")
    }

    if ("Trait" %in% names(down_plot)) {
      p_s10c <- down_plot %>%
        ggplot(aes(x = Trait_label, y = down_delta_value, fill = as.character(Module))) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
        geom_boxplot(width = 0.65, outlier.size = 1.8) +
        coord_flip() +
        module_fill_scale() +
        labs(title = "Downsampling robustness by trait", x = NULL, y = "Mean absolute delta rho", fill = "Module") +
        theme_pub()
    } else {
      p_s10c <- blank_panel("Trait column not available")
    }

    if (!is.null(loso_module)) {

      mean_delta_col_loso <- safe_col(loso_module, c("mean_abs_delta_rho"))
      if (!is.na(mean_delta_col_loso)) {
        p_s10d <- loso_module %>%
          dplyr::mutate(
            mean_abs_delta_rho = suppressWarnings(as.numeric(.data[[mean_delta_col_loso]])),
            Module = factor(Module, levels = Module[order(mean_abs_delta_rho)])
          ) %>%
          ggplot(aes(x = Module, y = mean_abs_delta_rho, fill = as.character(Module))) +
          geom_col(width = 0.75) +
          coord_flip() +
          module_fill_scale() +
          labs(title = "Leave-one-site-out robustness", x = NULL, y = "Mean absolute delta rho") +
          theme_pub() +
          theme(legend.position = "none")
      } else {
        p_s10d <- blank_panel("LOSO mean delta column not available")
      }

      fig_s10 <- (p_s10a + p_s10b) / (p_s10c + p_s10d) +
        patchwork::plot_annotation(tag_levels = "a", theme = plot_tag_theme())

      save_pub_both(fig_s10, file.path(OUT_FIG, "Supplementary_Fig_10_downsampling_LOSO_robustness"), dims = c(15.5, 12))

    } else {

      fig_s10 <- (p_s10a + p_s10b) / p_s10c +
        patchwork::plot_annotation(tag_levels = "a", theme = plot_tag_theme())

      save_pub_both(fig_s10, file.path(OUT_FIG, "Supplementary_Fig_10_downsampling_robustness"), dims = c(15.5, 12))
    }
  }
}

##############################################################################
# 13) EXPORT SUMMARY TABLE FOR MANUSCRIPT CHECKING
##############################################################################

supp_fig_inventory <- tibble::tibble(
  Supplementary_Figure = c(
    "Supplementary Fig. 4",
    "Supplementary Fig. 5",
    "Supplementary Fig. 6",
    "Supplementary Fig. 7",
    "Supplementary Fig. 8",
    "Supplementary Fig. 9",
    "Supplementary Fig. 10"
  ),
  Purpose = c(
    "Complete WGCNA module-trait association matrix with BH-FDR annotations",
    "Soft-thresholding diagnostics supporting the selected beta",
    "Module size and DEP burden across WGCNA modules",
    "Top hub proteins by module and module membership structure",
    "GO BP, KEGG, and Reactome enrichment term burden by module",
    "Leave-one-country-out robustness of module-trait associations",
    "Balanced downsampling and optional leave-one-site-out robustness"
  ),
  Output_file_base = c(
    "Supplementary_Fig_4_complete_module_trait_heatmap",
    "Supplementary_Fig_5_soft_thresholding",
    "Supplementary_Fig_6_module_size_DEP_burden",
    "Supplementary_Fig_7_hub_proteins_by_module",
    "Supplementary_Fig_8_enrichment_summary",
    "Supplementary_Fig_9_LOCO_country_robustness",
    "Supplementary_Fig_10_downsampling_LOSO_robustness or Supplementary_Fig_10_downsampling_robustness"
  )
)

readr::write_csv(
  supp_fig_inventory,
  file.path(OUT_TAB, "supplementary_figure_inventory.csv")
)

cat("\nDONE.\n")
cat("Supplementary figures saved in:\n", OUT_FIG, "\n")
cat("Inventory saved in:\n", file.path(OUT_TAB, "supplementary_figure_inventory.csv"), "\n")


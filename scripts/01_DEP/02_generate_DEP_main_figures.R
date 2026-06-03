###############################################################################
# REDLAT PLASMA PROTEOMICS
# SCRIPT 02 — MAIN MANUSCRIPT FIGURE GENERATION
#
# Purpose:
# - Generate publication-ready main and supplementary figures from Script 01
#   analysis outputs.
# - Avoid rerunning statistical models; this script only reads finalized tables
#   and workspace objects.
# - Use the gene-collapsed DEP table as the default interpretation layer for
#   figure-level summaries.
# - Export both final figures and source-data tables for manuscript submission.
# - Read canonical outputs from Script 01 v5, including APOE-, AT(N)- and AD-only CDR-SB sensitivity outputs.
#
# Expected input:
# - Workspace and tables produced by:
#     01_data_processing_and_differential_analysis.R
# - Preferred workspace paths:
#     result/workspace/proteomics_master_analysis_workspace.RData
#     result/workspace/analysis_workspace.RData
#     proteomics_master_reanalysis_workspace.RData
# - Canonical result folders:
#     result/03_dep/
#     result/04_sensitivity/
#     result/05_enrichment_corrected/
#     result/06_robustness/
#
# Main outputs:
# - result/final_figures/Figure2/
# - result/final_figures/Supplementary/
# - result/final_source_data/Figure2/
# - result/final_source_data/Supplementary/
# - result/final_figures/figure_manifest.csv
#
# Manuscript figure logic:
# - Figure 2a: main AD vs CN volcano plot.
# - Figure 2b: directional Reactome GSEA.
# - Figure 2c: representative clinical, cognitive, and plasma AT(N) heatmap.
# - Figure 2d: AT(N)-adjusted attenuation of diagnostic proteomic effects.
# - Figure 2e: contextual country-exclusion summary.
# - Figure 2f: protein-level internal multicountry LOCO stability.
# - AD-only CDR-SB severity alignment and CDR-SB-adjusted AD-vs-CN
#   diagnostic attenuation are exported as supplementary panels.
##############################################################################################################################################################

###############################################################################
# 00_config_packages_paths
###############################################################################

required_cran <- c(
  "dplyr", "tidyr", "purrr", "ggplot2", "ggrepel", "readr", "stringr",
  "forcats", "tibble", "patchwork", "scales", "grid", "rlang", "openxlsx"
)

missing_cran <- required_cran[!vapply(required_cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran) > 0) {
  stop(
    "Missing required CRAN packages: ", paste(missing_cran, collapse = ", "),
    "\nInstall them before running Script 02."
  )
}

required_bioc <- c("ComplexHeatmap", "circlize")
missing_bioc <- required_bioc[!vapply(required_bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc) > 0) {
  stop(
    "Missing required Bioconductor packages: ", paste(missing_bioc, collapse = ", "),
    "\nInstall them with BiocManager before running Script 02."
  )
}

invisible(lapply(c(required_cran, required_bioc), library, character.only = TRUE))
options(stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------
# USER PATH: must match project_root from Script 01.
# -----------------------------------------------------------------------------
project_root <- "C:/Users/mnpiz/Desktop/DEPs_Proteomic_Publishable_V2"
outdir <- project_root
setwd(outdir)

MAIN_FDR <- if (exists("MAIN_FDR")) MAIN_FDR else 0.05
STRICT_FDR <- if (exists("STRICT_FDR")) STRICT_FDR else 0.01

###############################################################################
# 01_helpers
###############################################################################

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

safe_file_tag <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

first_existing_file <- function(paths) {
  hit <- paths[file.exists(paths)][1]
  if (length(hit) == 0 || is.na(hit)) return(NA_character_)
  hit
}

read_csv_if_exists <- function(path, show_message = TRUE) {
  if (file.exists(path)) {
    if (show_message) message("Reading: ", path)
    return(readr::read_csv(path, show_col_types = FALSE))
  }
  if (show_message) message("Missing file: ", path)
  NULL
}

require_table_file <- function(paths, label = NULL) {
  path <- if (length(paths) == 1) paths else first_existing_file(paths)
  if (is.na(path) || !file.exists(path)) {
    if (is.null(label)) label <- paste(paths, collapse = " | ")
    stop("Required table not found: ", label)
  }
  readr::read_csv(path, show_col_types = FALSE)
}

rename_with_candidates <- function(df, mapping_list) {
  if (is.null(df) || !is.data.frame(df)) return(df)
  out <- df
  for (target in names(mapping_list)) {
    if (target %in% names(out)) next
    candidates <- mapping_list[[target]]
    found <- candidates[candidates %in% names(out)][1]
    if (!is.na(found)) {
      out <- out %>% dplyr::rename(!!target := !!rlang::sym(found))
    }
  }
  out
}

safe_numeric_vec <- function(x) {
  if (is.numeric(x)) return(x)
  suppressWarnings(as.numeric(as.character(x)))
}

safe_cor <- function(x, y, method = "pearson") {
  x <- safe_numeric_vec(x)
  y <- safe_numeric_vec(y)
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  suppressWarnings(as.numeric(stats::cor(x[ok], y[ok], method = method)))
}

safe_lm_slope <- function(x, y) {
  x <- safe_numeric_vec(x)
  y <- safe_numeric_vec(y)
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  suppressWarnings(as.numeric(stats::coef(stats::lm(y[ok] ~ x[ok]))[2]))
}

compute_pairwise_spearman_local <- function(data, protein_cols, trait_cols, min_n = 30) {
  out <- vector("list", length(protein_cols) * length(trait_cols))
  k <- 1
  
  for (p in protein_cols) {
    for (t in trait_cols) {
      x <- safe_numeric_vec(data[[p]])
      y <- safe_numeric_vec(data[[t]])
      ok <- is.finite(x) & is.finite(y)
      n_ok <- sum(ok)
      
      if (n_ok >= min_n) {
        ct <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
        out[[k]] <- tibble::tibble(
          AptName = p,
          trait = t,
          n = n_ok,
          rho = unname(ct$estimate),
          p_value = ct$p.value
        )
      } else {
        out[[k]] <- tibble::tibble(
          AptName = p,
          trait = t,
          n = n_ok,
          rho = NA_real_,
          p_value = NA_real_
        )
      }
      k <- k + 1
    }
  }
  
  dplyr::bind_rows(out) %>%
    dplyr::mutate(
      q_value_bh = stats::p.adjust(p_value, method = "BH"),
      sig_bh_005 = !is.na(q_value_bh) & q_value_bh < 0.05,
      sig_bh_001 = !is.na(q_value_bh) & q_value_bh < 0.01
    )
}

capture_complex_heatmap <- function(ht,
                                    heatmap_legend_side = "right",
                                    annotation_legend_side = "right",
                                    merge_legends = FALSE,
                                    padding = grid::unit(c(4, 5, 5, 4), "mm")) {
  grid::grid.grabExpr(
    ComplexHeatmap::draw(
      ht,
      heatmap_legend_side = heatmap_legend_side,
      annotation_legend_side = annotation_legend_side,
      merge_legends = merge_legends,
      padding = padding
    )
  )
}

save_complex_heatmap_pdf_png <- function(ht,
                                         pdf_file,
                                         png_file = NULL,
                                         width_cm = 18,
                                         height_cm = 14,
                                         dpi = 600,
                                         heatmap_legend_side = "right",
                                         annotation_legend_side = "right",
                                         merge_legends = FALSE,
                                         padding = grid::unit(c(5, 5, 5, 5), "mm")) {
  ensure_dir(dirname(pdf_file))
  grDevices::cairo_pdf(pdf_file, width = width_cm / 2.54, height = height_cm / 2.54)
  ComplexHeatmap::draw(
    ht,
    heatmap_legend_side = heatmap_legend_side,
    annotation_legend_side = annotation_legend_side,
    merge_legends = merge_legends,
    padding = padding
  )
  grDevices::dev.off()
  
  if (!is.null(png_file)) {
    ensure_dir(dirname(png_file))
    grDevices::png(png_file, width = width_cm, height = height_cm, units = "cm", res = dpi, bg = "white")
    ComplexHeatmap::draw(
      ht,
      heatmap_legend_side = heatmap_legend_side,
      annotation_legend_side = annotation_legend_side,
      merge_legends = merge_legends,
      padding = padding
    )
    grDevices::dev.off()
  }
}

save_gg_pdf_png_cm <- function(plot_obj,
                               file_base,
                               width_cm,
                               height_cm,
                               dpi = 600,
                               limitsize = FALSE) {
  ensure_dir(dirname(file_base))
  ggplot2::ggsave(
    filename = paste0(file_base, ".pdf"),
    plot = plot_obj,
    width = width_cm,
    height = height_cm,
    units = "cm",
    device = grDevices::cairo_pdf,
    bg = "white",
    limitsize = limitsize
  )
  ggplot2::ggsave(
    filename = paste0(file_base, ".png"),
    plot = plot_obj,
    width = width_cm,
    height = height_cm,
    units = "cm",
    dpi = dpi,
    bg = "white",
    limitsize = limitsize
  )
}

add_manifest <- function(manifest, figure_id, file_base, source_data = NA_character_, note = NA_character_) {
  dplyr::bind_rows(
    manifest,
    tibble::tibble(
      figure_id = figure_id,
      pdf = paste0(file_base, ".pdf"),
      png = paste0(file_base, ".png"),
      source_data = source_data,
      note = note
    )
  )
}

placeholder_plot <- function(label, base_size = 10) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0, label = label, size = 3.5) +
    ggplot2::xlim(-1, 1) +
    ggplot2::ylim(-1, 1) +
    ggplot2::theme_void(base_size = base_size)
}

###############################################################################
# 02_visual_system
###############################################################################

COL_UP        <- "#f46d43"
COL_DOWN      <- "#4682b4"
COL_NEUTRAL   <- "grey70"
COL_NS        <- "grey86"
COL_TEXT      <- "black"
COL_BG        <- "white"
COL_BORDER    <- "black"
COL_AD        <- COL_UP
COL_CN        <- COL_DOWN
COL_PRESERVED <- "#2c7fb8"
COL_DIRONLY   <- "grey70"
COL_NOTPRES   <- "white"
COL_MISSING   <- "grey94"

# Nature Aging typography target:
# - Final figure width: 18 cm
# - Panel titles: Arial regular, 10 pt
# - Axis/legend/tick text: Arial regular, 8 pt
# - Heatmap compacted so it does not dominate the composite.
FONT_FAMILY <- "Arial"
PANEL_TITLE_SIZE <- 10
TEXT_SIZE <- 8
TAG_SIZE <- 10

COMPOSITE_WIDTH_CM  <- 18
COMPOSITE_HEIGHT_CM <- 19
PANEL_REVIEW_W_CM <- 18
PANEL_REVIEW_H_CM <- 12

# Large individual panels.
theme_panel_review <- function(base_size = TEXT_SIZE, base_family = FONT_FAMILY) {
  ggplot2::theme_classic(base_size = TEXT_SIZE, base_family = FONT_FAMILY) +
    ggplot2::theme(
      text = ggplot2::element_text(family = FONT_FAMILY, face = "plain", size = TEXT_SIZE, colour = COL_TEXT),
      axis.text = ggplot2::element_text(family = FONT_FAMILY, face = "plain", size = TEXT_SIZE, colour = COL_TEXT),
      axis.title = ggplot2::element_text(family = FONT_FAMILY, face = "plain", size = TEXT_SIZE, colour = COL_TEXT),
      axis.line = ggplot2::element_line(linewidth = 0.35, colour = COL_BORDER),
      axis.ticks = ggplot2::element_line(linewidth = 0.30, colour = COL_BORDER),
      panel.background = ggplot2::element_rect(fill = COL_BG, colour = NA),
      plot.background = ggplot2::element_rect(fill = COL_BG, colour = NA),
      panel.grid = ggplot2::element_blank(),
      legend.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(family = FONT_FAMILY, face = "plain", size = TEXT_SIZE),
      legend.text = ggplot2::element_text(family = FONT_FAMILY, face = "plain", size = TEXT_SIZE),
      plot.title = ggplot2::element_text(family = FONT_FAMILY, face = "plain", size = PANEL_TITLE_SIZE, hjust = 0, colour = COL_TEXT, margin = ggplot2::margin(b = 3)),
      plot.subtitle = ggplot2::element_text(family = FONT_FAMILY, face = "plain", size = TEXT_SIZE, hjust = 0, colour = COL_TEXT),
      plot.margin = ggplot2::margin(4, 5, 4, 4)
    )
}

# Compact composite panels.
theme_panel_composite <- function(base_size = TEXT_SIZE, base_family = FONT_FAMILY) {
  ggplot2::theme_classic(base_size = TEXT_SIZE, base_family = FONT_FAMILY) +
    ggplot2::theme(
      text = ggplot2::element_text(family = FONT_FAMILY, face = "plain", size = TEXT_SIZE, colour = COL_TEXT),
      axis.text = ggplot2::element_text(family = FONT_FAMILY, face = "plain", size = TEXT_SIZE, colour = COL_TEXT),
      axis.title = ggplot2::element_text(family = FONT_FAMILY, face = "plain", size = TEXT_SIZE, colour = COL_TEXT),
      axis.line = ggplot2::element_line(linewidth = 0.30, colour = COL_BORDER),
      axis.ticks = ggplot2::element_line(linewidth = 0.25, colour = COL_BORDER),
      panel.background = ggplot2::element_rect(fill = COL_BG, colour = NA),
      plot.background = ggplot2::element_rect(fill = COL_BG, colour = NA),
      panel.grid = ggplot2::element_blank(),
      legend.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(family = FONT_FAMILY, face = "plain", size = TEXT_SIZE),
      legend.text = ggplot2::element_text(family = FONT_FAMILY, face = "plain", size = TEXT_SIZE),
      plot.title = ggplot2::element_text(family = FONT_FAMILY, face = "plain", size = PANEL_TITLE_SIZE, hjust = 0, colour = COL_TEXT, margin = ggplot2::margin(b = 2)),
      plot.subtitle = ggplot2::element_text(family = FONT_FAMILY, face = "plain", size = TEXT_SIZE, hjust = 0, colour = COL_TEXT),
      plot.margin = ggplot2::margin(3, 4, 3, 3)
    )
}

###############################################################################
# 03_load_workspace_and_tables
###############################################################################

ensure_dir(file.path(outdir, "result", "final_figures", "Figure2"))
ensure_dir(file.path(outdir, "result", "final_figures", "Supplementary"))
ensure_dir(file.path(outdir, "result", "final_source_data", "Figure2"))
ensure_dir(file.path(outdir, "result", "final_source_data", "Supplementary"))
ensure_dir(file.path(outdir, "result", "final_figures", "logs"))

figure_manifest <- tibble::tibble(
  figure_id = character(),
  pdf = character(),
  png = character(),
  source_data = character(),
  note = character()
)

workspace_file <- first_existing_file(c(
  file.path(outdir, "result", "workspace", "proteomics_master_analysis_workspace.RData"),
  file.path(outdir, "result", "workspace", "analysis_workspace.RData"),
  file.path(outdir, "result", "workspace", "proteomics_master_reanalysis_workspace.RData"),
  file.path(outdir, "proteomics_master_analysis_workspace.RData"),
  file.path(outdir, "analysis_workspace.RData"),
  file.path(outdir, "proteomics_master_reanalysis_workspace.RData"),
  file.path(outdir, "result", "analysis_workspace.RData"),
  file.path(outdir, "result", "proteomics_master_reanalysis_workspace.RData")
))

if (!is.na(workspace_file)) {
  message("Loading workspace: ", workspace_file)
  load(workspace_file)
} else {
  warning("No workspace found. The script will try to read required tables from result/ folders.")
}

MAIN_FDR <- if (exists("MAIN_FDR")) MAIN_FDR else 0.05
STRICT_FDR <- if (exists("STRICT_FDR")) STRICT_FDR else 0.01

# Main DEP: prefer gene-collapsed.
if (exists("DEP_gene") && is.data.frame(DEP_gene)) {
  DEP <- DEP_gene
} else if (exists("DEP") && is.data.frame(DEP)) {
  DEP <- DEP
} else {
  DEP <- require_table_file(
    file.path(outdir, "result", "03_dep", "gene_collapsed", "AD_vs_CN_full_limma_results_gene_collapsed.csv"),
    label = "main DEP gene-collapsed table from result/03_dep"
  )
}

DEP <- rename_with_candidates(
  DEP,
  list(
    Protein_Name     = c("Protein_Name", "ProteinLabel", "Target", "TargetFullName", "EntrezGeneSymbol"),
    EntrezGeneSymbol = c("EntrezGeneSymbol", "Gene", "Symbol"),
    AptName          = c("AptName"),
    SeqId            = c("SeqId", "seqid"),
    logFC            = c("logFC", "main_logFC"),
    adj.P.Val        = c("adj.P.Val", "FDR", "qvalue", "q_value", "main_adj.P.Val"),
    P.Value          = c("P.Value", "p_value", "pvalue"),
    type             = c("type", "Direction", "direction")
  )
)

if (!all(c("AptName", "Protein_Name", "logFC", "adj.P.Val") %in% names(DEP))) {
  stop("DEP table is missing required columns after standardization.")
}

DEP <- DEP %>%
  dplyr::mutate(
    logFC = safe_numeric_vec(logFC),
    adj.P.Val = safe_numeric_vec(adj.P.Val),
    minus_log10_fdr = -log10(pmax(adj.P.Val, 1e-300)),
    Direction = dplyr::case_when(
      adj.P.Val < MAIN_FDR & logFC > 0 ~ "Higher in AD",
      adj.P.Val < MAIN_FDR & logFC < 0 ~ "Lower in AD",
      TRUE ~ "Not significant"
    ),
    Direction = factor(Direction, levels = c("Higher in AD", "Lower in AD", "Not significant")),
    Protein_Label = dplyr::coalesce(Protein_Name, EntrezGeneSymbol, AptName)
  )

# Normalized expression for heatmap and PCA.
if (exists("normalized_expr_CN_AD") && is.data.frame(normalized_expr_CN_AD)) {
  normalized_expr_plot <- normalized_expr_CN_AD
} else if (exists("normalized_expr") && is.data.frame(normalized_expr)) {
  normalized_expr_plot <- normalized_expr
} else {
  stop("Required object 'normalized_expr_CN_AD' or 'normalized_expr' not found. Run Script 01 first.")
}

if (exists("seq_cols")) {
  seq_cols_plot <- intersect(seq_cols, names(normalized_expr_plot))
} else if (exists("protein_universe")) {
  seq_cols_plot <- intersect(protein_universe, names(normalized_expr_plot))
} else {
  seq_cols_plot <- names(normalized_expr_plot)[grepl("^seq[._]", names(normalized_expr_plot), ignore.case = TRUE)]
}
if (length(seq_cols_plot) == 0) stop("Could not infer proteomic columns from normalized_expr.")

# Optional tables, new paths first, aliases second.
main_vs_loco_mean <- if (exists("main_vs_loco_mean") && is.data.frame(main_vs_loco_mean)) {
  main_vs_loco_mean
} else {
  read_csv_if_exists(file.path(outdir, "result", "06_robustness", "country_loco", "tables", "main_vs_meanLOCO_table.csv"), FALSE)
}

apoe_compare_tbl <- if (exists("apoe_compare_tbl") && is.data.frame(apoe_compare_tbl)) {
  apoe_compare_tbl
} else {
  read_csv_if_exists(file.path(outdir, "result", "04_sensitivity", "apoe", "primary_vs_APOE_adjusted_gene_comparison.csv"), FALSE)
}

cdr_compare_tbl <- if (exists("cdrsb_compare_tbl") && is.data.frame(cdrsb_compare_tbl)) {
  cdrsb_compare_tbl
} else if (exists("cdr_compare_tbl") && is.data.frame(cdr_compare_tbl)) {
  cdr_compare_tbl
} else {
  read_csv_if_exists(file.path(outdir, "result", "04_sensitivity", "cdrsb", "AD_only", "primary_AD_vs_CN_vs_AD_only_CDRSB_severity_alignment.csv"), FALSE)
}

# Secondary diagnostic attenuation model: AD vs CN additionally adjusted for CDR-SB.
# This is intentionally kept separate from the AD-only CDR-SB severity model.
cdrsb_adjusted_compare_tbl <- if (exists("cdrsb_adjusted_compare_tbl") && is.data.frame(cdrsb_adjusted_compare_tbl)) {
  cdrsb_adjusted_compare_tbl
} else if (exists("cdrsb_adjusted_ad_vs_cn_compare_tbl") && is.data.frame(cdrsb_adjusted_ad_vs_cn_compare_tbl)) {
  cdrsb_adjusted_ad_vs_cn_compare_tbl
} else {
  read_csv_if_exists(file.path(outdir, "result", "04_sensitivity", "cdrsb", "AD_vs_CN_adjusted", "primary_vs_CDRSB_adjusted_AD_vs_CN_gene_comparison.csv"), FALSE)
}

vascular_compare_tbl <- if (exists("vascular_compare_tbl") && is.data.frame(vascular_compare_tbl)) {
  vascular_compare_tbl
} else {
  read_csv_if_exists(first_existing_file(c(
    file.path(outdir, "result", "04_sensitivity", "vascular_metabolic", "main_vs_vascular_metabolic_adjusted_gene_comparison.csv")
  )), FALSE)
}

atn_compare_tbl <- if (exists("atn_compare_tbl") && is.data.frame(atn_compare_tbl)) {
  atn_compare_tbl
} else if (exists("ATN_compare_tbl") && is.data.frame(ATN_compare_tbl)) {
  ATN_compare_tbl
} else {
  read_csv_if_exists(file.path(outdir, "result", "04_sensitivity", "atn_adjusted", "primary_vs_ATN_adjusted_gene_comparison.csv"), FALSE)
}

# Reactome GSEA: corrected new path first, old path second.
gsea_reactome <- read_csv_if_exists(file.path(outdir, "result", "05_enrichment_corrected", "gsea", "main_dep_gsea_reactome_bh.csv"), FALSE)

# Robustness optional.
robustness_tbl <- if (exists("robustness_classification_tbl") && is.data.frame(robustness_classification_tbl)) {
  robustness_classification_tbl
} else {
  read_csv_if_exists(file.path(outdir, "result", "06_robustness", "formal_classification", "protein_robustness_classification.csv"), FALSE)
}

###############################################################################
# 04_common_plot_builders
###############################################################################

make_volcano_plot <- function(df,
                              label_df = NULL,
                              compact = FALSE,
                              show_legend = TRUE) {
  theme_use <- if (compact) theme_panel_composite() else theme_panel_review()
  label_size <- if (compact) 2.0 else 3.1
  point_sig <- if (compact) 1.25 else 2.0
  point_ns  <- if (compact) 0.75 else 1.4
  
  df <- df %>% dplyr::filter(is.finite(logFC), is.finite(minus_log10_fdr))
  xmax <- max(0.6, ceiling(max(abs(df$logFC), na.rm = TRUE) * 10) / 10)
  ymax <- ceiling(max(df$minus_log10_fdr, na.rm = TRUE) * 10) / 10 + ifelse(compact, 0.15, 0.45)
  
  p <- ggplot() +
    geom_hline(yintercept = -log10(MAIN_FDR), linetype = 2, linewidth = 0.35, colour = "grey45") +
    geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey55") +
    geom_point(
      data = df %>% dplyr::filter(Direction == "Not significant"),
      aes(x = logFC, y = minus_log10_fdr),
      shape = 16,
      size = point_ns,
      colour = "grey88",
      alpha = 0.35
    ) +
    geom_point(
      data = df %>% dplyr::filter(Direction != "Not significant"),
      aes(x = logFC, y = minus_log10_fdr, colour = Direction),
      shape = 16,
      size = point_sig,
      alpha = 0.88
    ) +
    scale_colour_manual(
      values = c("Higher in AD" = COL_UP, "Lower in AD" = COL_DOWN, "Not significant" = COL_NS),
      breaks = c("Higher in AD", "Lower in AD"),
      drop = FALSE
    ) +
    coord_cartesian(xlim = c(-xmax, xmax), ylim = c(0, ymax), expand = FALSE, clip = "off") +
    labs(
      x = if (compact) "log2FC" else "Adjusted log2 fold-change (AD vs CN)",
      y = expression(-log[10](FDR)),
      colour = NULL
    ) +
    theme_use +
    theme(legend.position = if (show_legend) "top" else "none")
  
  if (!is.null(label_df) && nrow(label_df) > 0) {
    p <- p +
      ggrepel::geom_text_repel(
        data = label_df,
        aes(x = logFC, y = minus_log10_fdr, label = Protein_Label),
        size = label_size,
        seed = 123,
        max.overlaps = Inf,
        force = if (compact) 0.7 else 1.1,
        box.padding = if (compact) 0.12 else 0.25,
        point.padding = if (compact) 0.05 else 0.12,
        min.segment.length = 0,
        segment.color = "grey55",
        segment.size = if (compact) 0.18 else 0.25,
        show.legend = FALSE
      )
  }
  
  p
}

make_direction_scatter <- function(df,
                                   xvar,
                                   yvar,
                                   label_df = NULL,
                                   label_var = "Protein_Label",
                                   compact = FALSE,
                                   show_stats = TRUE,
                                   xlab = "",
                                   ylab = "") {
  theme_use <- if (compact) theme_panel_composite() else theme_panel_review()
  label_size <- if (compact) 2.0 else 2.8
  point_size <- if (compact) 1.0 else 1.7
  
  plot_df <- df %>% dplyr::filter(is.finite(.data[[xvar]]), is.finite(.data[[yvar]]))
  stats_lab <- paste0(
    "r=", format(round(safe_cor(plot_df[[xvar]], plot_df[[yvar]]), 2), nsmall = 2),
    if (!compact) "\n" else ", ",
    "slope=", format(round(safe_lm_slope(plot_df[[xvar]], plot_df[[yvar]]), 2), nsmall = 2)
  )
  
  p <- ggplot(plot_df, aes(x = .data[[xvar]], y = .data[[yvar]])) +
    geom_hline(yintercept = 0, colour = COL_NEUTRAL, linewidth = 0.25) +
    geom_vline(xintercept = 0, colour = COL_NEUTRAL, linewidth = 0.25) +
    geom_abline(intercept = 0, slope = 1, linetype = 2, colour = "grey55", linewidth = 0.32) +
    geom_point(aes(colour = Direction), shape = 16, size = point_size, alpha = 0.52) +
    geom_smooth(method = "lm", se = FALSE, colour = COL_TEXT, linewidth = if (compact) 0.35 else 0.55) +
    scale_colour_manual(
      values = c("Higher in AD" = COL_UP, "Lower in AD" = COL_DOWN, "Neutral" = "grey65", "Not significant" = "grey65"),
      guide = "none"
    ) +
    labs(x = xlab, y = ylab) +
    theme_use
  
  if (show_stats) {
    p <- p + annotate("text", x = Inf, y = -Inf, label = stats_lab, hjust = 1.04, vjust = -0.25, size = if (compact) 2.15 else 3.2)
  }
  
  if (!is.null(label_df) && nrow(label_df) > 0) {
    p <- p +
      ggrepel::geom_text_repel(
        data = label_df,
        aes(label = .data[[label_var]]),
        size = label_size,
        seed = 123,
        max.overlaps = Inf,
        box.padding = if (compact) 0.10 else 0.22,
        point.padding = if (compact) 0.05 else 0.10,
        segment.color = "grey55",
        segment.size = if (compact) 0.16 else 0.24,
        show.legend = FALSE
      )
  }
  
  p
}

make_density_scatter <- function(df,
                                 xvar,
                                 yvar,
                                 label_df = NULL,
                                 label_var = "Protein_Label",
                                 compact = FALSE,
                                 xlab = "",
                                 ylab = "") {
  theme_use <- if (compact) theme_panel_composite() else theme_panel_review()
  label_size <- if (compact) 2.0 else 2.8
  
  plot_df <- df %>% dplyr::filter(is.finite(.data[[xvar]]), is.finite(.data[[yvar]]))
  stats_lab <- paste0(
    "r=", format(round(safe_cor(plot_df[[xvar]], plot_df[[yvar]]), 2), nsmall = 2),
    if (!compact) "\n" else ", ",
    "slope=", format(round(safe_lm_slope(plot_df[[xvar]], plot_df[[yvar]]), 2), nsmall = 2)
  )
  
  p <- ggplot(plot_df, aes(x = .data[[xvar]], y = .data[[yvar]])) +
    geom_hline(yintercept = 0, colour = COL_NEUTRAL, linewidth = 0.25) +
    geom_vline(xintercept = 0, colour = COL_NEUTRAL, linewidth = 0.25) +
    geom_abline(intercept = 0, slope = 1, linetype = 2, colour = "grey55", linewidth = 0.32) +
    stat_bin2d(aes(fill = after_stat(count)), bins = if (compact) 25 else 35, alpha = 0.95) +
    scale_fill_gradient(low = "grey94", high = "grey45", guide = if (compact) "none" else guide_colorbar(title = "Density")) +
    geom_smooth(method = "lm", se = FALSE, colour = COL_TEXT, linewidth = if (compact) 0.35 else 0.55) +
    labs(x = xlab, y = ylab) +
    theme_use +
    annotate("text", x = Inf, y = -Inf, label = stats_lab, hjust = 1.04, vjust = -0.25, size = if (compact) 2.15 else 3.2)
  
  if (!is.null(label_df) && nrow(label_df) > 0) {
    p <- p +
      geom_point(
        data = label_df,
        aes(x = .data[[xvar]], y = .data[[yvar]], colour = Direction),
        inherit.aes = FALSE,
        shape = 16,
        size = if (compact) 1.2 else 2.0,
        alpha = 0.95,
        show.legend = FALSE
      ) +
      ggrepel::geom_text_repel(
        data = label_df,
        aes(x = .data[[xvar]], y = .data[[yvar]], label = .data[[label_var]], colour = Direction),
        inherit.aes = FALSE,
        size = label_size,
        seed = 123,
        max.overlaps = Inf,
        box.padding = if (compact) 0.10 else 0.22,
        point.padding = if (compact) 0.05 else 0.10,
        segment.color = "grey55",
        segment.size = if (compact) 0.16 else 0.24,
        show.legend = FALSE
      ) +
      scale_colour_manual(values = c("Higher in AD" = COL_UP, "Lower in AD" = COL_DOWN, "Neutral" = "grey55", "Not significant" = "grey55"))
  }
  
  p
}

###############################################################################
# 05_select_proteins_and_traits
###############################################################################

# Curated labels. If absent, fallback to most significant proteins.
label_targets <- c("SPC25", "LRRN1", "SMOC1", "CPLX2", "ACHE", "IGFBP2", "C3", "AMBRA1", "GDI1", "PHGDH", "PUM2", "CAPN1")

label_df_main <- DEP %>%
  dplyr::filter(Protein_Label %in% label_targets | EntrezGeneSymbol %in% label_targets) %>%
  dplyr::filter(is.finite(logFC), is.finite(minus_log10_fdr)) %>%
  dplyr::arrange(adj.P.Val)

if (nrow(label_df_main) < 6) {
  label_df_main <- DEP %>%
    dplyr::filter(Direction != "Not significant", is.finite(logFC), is.finite(minus_log10_fdr)) %>%
    dplyr::arrange(adj.P.Val, dplyr::desc(abs(logFC))) %>%
    dplyr::slice_head(n = 8)
}

trait_candidates <- list(
  cdr_boxscore = c("cdr_boxscore", "cdr_sb", "cdr_sum_box", "CDR_SB"),
  pfaq = c("udsfaq_total", "pfaq_total", "PFAQ", "pfaq"),
  mmse = c("mmse_total", "MMSE", "mmse"),
  mini_sea = c("Mini.SEA", "mini_sea", "MiniSEA", "mini_sea_total"),
  npi = c("npi_total", "NPI", "NPI_Q", "npi_q", "NPI_Q_total"),
  ptau181 = c("p.tau181", "p_tau181", "ptau181", "p-tau181"),
  ptau217 = c("p.tau217", "p_tau217", "ptau217", "p-tau217"),
  nfl = c("NfL", "nfl", "NFL"),
  gfap = c("GFAP", "gfap"),
  abratio = c("ratio.AB42.40", "ratio_AB42_40", "AB42_40_ratio", "Aβ42.40")
)

trait_vars_detected <- purrr::map_chr(trait_candidates, function(cands) {
  hit <- cands[cands %in% names(normalized_expr_plot)][1]
  ifelse(length(hit) == 0 || is.na(hit), NA_character_, hit)
})
trait_vars <- unique(stats::na.omit(trait_vars_detected))

trait_display <- c(
  cdr_boxscore = "CDR-SB", cdr_sb = "CDR-SB", CDR_SB = "CDR-SB",
  udsfaq_total = "PFAQ", pfaq_total = "PFAQ", PFAQ = "PFAQ", pfaq = "PFAQ",
  mmse_total = "MMSE", MMSE = "MMSE", mmse = "MMSE",
  Mini.SEA = "Mini-SEA", mini_sea = "Mini-SEA", MiniSEA = "Mini-SEA",
  npi_total = "NPI-Q", NPI = "NPI-Q", NPI_Q = "NPI-Q",
  p.tau181 = "p-tau181", p_tau181 = "p-tau181", ptau181 = "p-tau181",
  p.tau217 = "p-tau217", p_tau217 = "p-tau217", ptau217 = "p-tau217",
  NfL = "NfL", nfl = "NfL", GFAP = "GFAP", gfap = "GFAP",
  ratio.AB42.40 = "Aβ42/40", ratio_AB42_40 = "Aβ42/40", AB42_40_ratio = "Aβ42/40"
)

###############################################################################
# 06_figure2a_volcano
###############################################################################

p_volcano_review <- make_volcano_plot(
  DEP,
  label_df = label_df_main,
  compact = FALSE,
  show_legend = TRUE
) +
  labs(title = "AD-associated plasma proteomic remodeling")

p_volcano_compact <- make_volcano_plot(
  DEP,
  label_df = label_df_main %>% dplyr::slice_head(n = 6),
  compact = TRUE,
  show_legend = FALSE
) +
  labs(title = "AD-associated plasma proteomic remodeling")

file_base <- file.path(outdir, "result", "final_figures", "Figure2", "Figure2a_main_volcano")
save_gg_pdf_png_cm(p_volcano_review, file_base, PANEL_REVIEW_W_CM, PANEL_REVIEW_H_CM)

source_file <- file.path(outdir, "result", "final_source_data", "Figure2", "Figure2a_main_volcano_source_data.csv")
readr::write_csv(
  DEP %>% dplyr::select(Protein_Label, Protein_Name, EntrezGeneSymbol, AptName, logFC, adj.P.Val, minus_log10_fdr, Direction),
  source_file
)
figure_manifest <- add_manifest(figure_manifest, "Figure2a_main_volcano", file_base, source_file, "Main gene-collapsed DEP volcano.")

###############################################################################
# 07_figure2b_reactome_gsea
###############################################################################

if (!is.null(gsea_reactome) && nrow(gsea_reactome) > 0) {
  
  gsea_reactome <- rename_with_candidates(
    gsea_reactome,
    list(
      Description = c("Description", "description"),
      NES = c("NES", "nes"),
      p.adjust = c("p.adjust", "qvalue", "q_value", "FDR"),
      setSize = c("setSize", "size")
    )
  )
  
  clean_reactome_label <- function(x) {
    x <- stringr::str_replace_all(x, "^REACTOME_", "")
    x <- stringr::str_replace_all(x, "_", " ")
    x <- stringr::str_to_lower(x)
    x <- stringr::str_replace_all(x, "\\bmrna\\b", "mRNA")
    x <- stringr::str_replace_all(x, "\\brna\\b", "RNA")
    x <- stringr::str_replace_all(x, "\\bpre-mRNA\\b", "pre-mRNA")
    x <- stringr::str_replace_all(x, "\\becm\\b", "ECM")
    x <- stringr::str_replace_all(x, "^([a-z])", toupper)
    x
  }
  
  gsea_plot_tbl <- gsea_reactome %>%
    dplyr::filter(
      !is.na(NES),
      !is.na(p.adjust),
      is.finite(NES),
      is.finite(p.adjust),
      p.adjust < 0.05
    ) %>%
    dplyr::mutate(
      Direction = ifelse(NES > 0, "Higher in AD", "Lower in AD"),
      Description_clean = clean_reactome_label(Description),
      Description_upper = toupper(Description_clean)
    )
  
  EXCLUDE_PATTERNS <- c(
    "DIGESTION",
    "ANTIMICROBIAL"
  )
  
  gsea_plot_tbl <- gsea_plot_tbl %>%
    dplyr::filter(
      !stringr::str_detect(
        toupper(Description_clean),
        paste(EXCLUDE_PATTERNS, collapse="|")
      )
    )
  
  BIO_PRIORITY_HIGHER <- c(
    "EXTRACELLULAR",
    "MATRIX",
    "COLLAGEN",
    "CELL SURFACE",
    "CELL ADHESION",
    "ADHESION",
    "RECEPTOR",
    "CYTOKINE",
    "GROWTH FACTOR",
    "SIGNALING"
  )
  
  BIO_PRIORITY_LOWER <- c(
    "RNA",
    "MRNA",
    "MRNA",
    "SPLIC",
    "POLYADENYL",
    "3'-END",
    "PRE-MRNA",
    "TRANSCRIPTION",
    "POST TRANSCRIPTION"
  )
  
  gsea_plot_tbl <- gsea_plot_tbl %>%
    dplyr::mutate(
      bio_priority_higher = stringr::str_detect(
        Description_upper,
        paste(BIO_PRIORITY_HIGHER, collapse = "|")
      ),
      bio_priority_lower = stringr::str_detect(
        Description_upper,
        paste(BIO_PRIORITY_LOWER, collapse = "|")
      )
    )
  
  top_higher <- gsea_plot_tbl %>%
    dplyr::filter(Direction == "Higher in AD") %>%
    dplyr::arrange(
      dplyr::desc(bio_priority_higher),
      p.adjust,
      dplyr::desc(abs(NES))
    ) %>%
    dplyr::slice_head(n = 6)
  
  top_lower <- gsea_plot_tbl %>%
    dplyr::filter(Direction == "Lower in AD") %>%
    dplyr::arrange(
      dplyr::desc(bio_priority_lower),
      p.adjust,
      dplyr::desc(abs(NES))
    ) %>%
    dplyr::slice_head(n = 6)
  
  gsea_selected_review <- dplyr::bind_rows(top_higher, top_lower) %>%
    dplyr::mutate(
      Description_display = stringr::str_wrap(Description_clean, 36),
      Description_display = factor(
        Description_display,
        levels = Description_display[order(NES)]
      )
    )
  
  gsea_selected_compact <- dplyr::bind_rows(
    top_higher %>% dplyr::slice_head(n = 3),
    top_lower %>% dplyr::slice_head(n = 3)
  ) %>%
    dplyr::mutate(
      Description_display = stringr::str_wrap(Description_clean, 23),
      Description_display = factor(
        Description_display,
        levels = Description_display[order(NES)]
      )
    )
  
  make_gsea_bar <- function(plot_tbl, compact = FALSE) {
    
    theme_use <- if (compact) theme_panel_composite() else theme_panel_review()
    
    xmax_nes <- max(
      2.0,
      ceiling(max(abs(plot_tbl$NES), na.rm = TRUE) * 10) / 10
    )
    
    ggplot(plot_tbl, aes(x = NES, y = Description_display, fill = Direction)) +
      geom_col(
        width = if (compact) 0.58 else 0.65,
        colour = "black",
        linewidth = if (compact) 0.18 else 0.25
      ) +
      geom_vline(xintercept = 0, linewidth = 0.32, colour = "black") +
      scale_fill_manual(
        values = c(
          "Higher in AD" = COL_UP,
          "Lower in AD" = COL_DOWN
        )
      ) +
      coord_cartesian(
        xlim = c(-xmax_nes, xmax_nes),
        expand = FALSE,
        clip = "off"
      ) +
      labs(
        x = "Normalized enrichment score (NES)",
        y = "Representative Reactome pathways",
        fill = NULL,
        title = if (!compact) "Functional asymmetry of the AD-associated plasma proteome" else NULL
      ) +
      theme_use +
      theme(
        legend.position = if (compact) "none" else "top",
        axis.text.y = element_text(size = if (compact) 6.1 else 8.8)
      )
  }
  
  p_gsea_review <- make_gsea_bar(gsea_selected_review, compact = FALSE)
  p_gsea_compact <- make_gsea_bar(gsea_selected_compact, compact = TRUE)
  
  source_file <- file.path(
    outdir,
    "result",
    "final_source_data",
    "Figure2",
    "Figure2b_reactome_gsea_source_data.csv"
  )
  
  readr::write_csv(
    gsea_selected_review %>%
      dplyr::select(
        Description,
        Description_clean,
        NES,
        p.adjust,
        setSize,
        Direction,
        bio_priority_higher,
        bio_priority_lower
      ),
    source_file
  )
  
} else {
  
  p_gsea_review <- placeholder_plot("Reactome GSEA file not found")
  p_gsea_compact <- placeholder_plot("Reactome GSEA not available")
  source_file <- NA_character_
}

file_base <- file.path(
  outdir,
  "result",
  "final_figures",
  "Figure2",
  "Figure2b_reactome_gsea"
)

save_gg_pdf_png_cm(
  p_gsea_review,
  file_base,
  PANEL_REVIEW_W_CM,
  PANEL_REVIEW_H_CM
)

figure_manifest <- add_manifest(
  figure_manifest,
  "Figure2b_reactome_gsea",
  file_base,
  source_file,
  "Directional Reactome GSEA terms; BH-FDR < 0.05; representative pathways prioritized for manuscript-aligned biological interpretation."
)
###############################################################################
# 08_figure2d_ATN_adjusted_attenuation
###############################################################################

if (!is.null(atn_compare_tbl) && nrow(atn_compare_tbl) > 0) {
  atn_compare_tbl <- rename_with_candidates(
    atn_compare_tbl,
    list(
      AptName = c("AptName"),
      Protein_Name = c("Protein_Name", "ProteinLabel", "Target"),
      EntrezGeneSymbol = c("EntrezGeneSymbol", "Gene", "Symbol"),
      logFC_primary = c("logFC_primary", "main_logFC"),
      adj.P.Val_primary = c("adj.P.Val_primary", "main_adj.P.Val"),
      logFC_atn = c("logFC_atn", "logFC_ATN", "logFC_secondary", "logFC_adjusted"),
      adj.P.Val_atn = c("adj.P.Val_atn", "adj.P.Val_ATN", "adj.P.Val_secondary", "adj.P.Val_adjusted")
    )
  )

  atn_plot_tbl <- atn_compare_tbl %>%
    dplyr::mutate(
      logFC_primary = safe_numeric_vec(logFC_primary),
      logFC_atn = safe_numeric_vec(logFC_atn),
      adj.P.Val_primary = safe_numeric_vec(adj.P.Val_primary),
      adj.P.Val_atn = safe_numeric_vec(adj.P.Val_atn),
      Protein_Label = dplyr::coalesce(Protein_Name, EntrezGeneSymbol, AptName),
      Direction = dplyr::case_when(
        logFC_primary > 0 ~ "Higher in AD",
        logFC_primary < 0 ~ "Lower in AD",
        TRUE ~ "Neutral"
      ),
      Direction = factor(Direction, levels = c("Higher in AD", "Lower in AD", "Neutral")),
      primary_fdr005 = !is.na(adj.P.Val_primary) & adj.P.Val_primary < MAIN_FDR,
      atn_fdr005 = !is.na(adj.P.Val_atn) & adj.P.Val_atn < MAIN_FDR,
      same_direction = sign(logFC_primary) == sign(logFC_atn),
      attenuation_ratio = abs(logFC_atn) / pmax(abs(logFC_primary), 1e-9)
    ) %>%
    dplyr::filter(is.finite(logFC_primary), is.finite(logFC_atn))

  label_atn <- atn_plot_tbl %>%
    dplyr::filter(
      Protein_Label %in% label_targets |
        EntrezGeneSymbol %in% label_targets |
        Protein_Label == "GDI1" |
        EntrezGeneSymbol == "GDI1" |
        atn_fdr005
    ) %>%
    dplyr::arrange(dplyr::desc(atn_fdr005), adj.P.Val_primary) %>%
    dplyr::slice_head(n = 7)

  if (nrow(label_atn) == 0) {
    label_atn <- atn_plot_tbl %>%
      dplyr::arrange(adj.P.Val_primary) %>%
      dplyr::slice_head(n = 6)
  }

  n_atn_sig <- sum(atn_plot_tbl$atn_fdr005, na.rm = TRUE)
  dc_atn <- mean(atn_plot_tbl$same_direction, na.rm = TRUE)
  stats_atn <- paste0(
    "r=", format(round(safe_cor(atn_plot_tbl$logFC_primary, atn_plot_tbl$logFC_atn), 2), nsmall = 2),
    "\nDirection consistency=", round(100 * dc_atn, 1), "%",
    "\nFDR<0.05 after AT(N)=", n_atn_sig
  )

  p_atn_review <- make_direction_scatter(
    atn_plot_tbl,
    "logFC_primary", "logFC_atn",
    label_df = label_atn,
    compact = FALSE,
    show_stats = FALSE,
    xlab = "Main-model log2 fold-change (AD vs CN)",
    ylab = "AT(N)-adjusted log2 fold-change"
  ) +
    ggplot2::annotate("text", x = Inf, y = -Inf, label = stats_atn,
                      hjust = 1.04, vjust = -0.20, size = 3.0) +
    ggplot2::labs(title = "AT(N)-adjusted attenuation of diagnostic proteomic effects")

  p_atn_compact <- make_direction_scatter(
    atn_plot_tbl,
    "logFC_primary", "logFC_atn",
    label_df = label_atn %>% dplyr::slice_head(n = 5),
    compact = TRUE,
    show_stats = FALSE,
    xlab = "Main log2FC",
    ylab = "AT(N)-adjusted log2FC"
  ) +
    ggplot2::annotate("text", x = Inf, y = -Inf,
                      label = paste0("r=", format(round(safe_cor(atn_plot_tbl$logFC_primary, atn_plot_tbl$logFC_atn), 2), nsmall = 2),
                                     ", FDR<0.05=", n_atn_sig),
                      hjust = 1.04, vjust = -0.20, size = 2.1)

  source_file <- file.path(outdir, "result", "final_source_data", "Figure2", "Figure2d_ATN_adjusted_attenuation_source_data.csv")
  readr::write_csv(
    atn_plot_tbl %>%
      dplyr::select(Protein_Label, Protein_Name, EntrezGeneSymbol, AptName,
                    logFC_primary, adj.P.Val_primary, logFC_atn, adj.P.Val_atn,
                    primary_fdr005, atn_fdr005, same_direction, attenuation_ratio, Direction),
    source_file
  )
} else {
  p_atn_review <- placeholder_plot("AT(N)-adjusted comparison table not found")
  p_atn_compact <- placeholder_plot("AT(N)-adjusted comparison not available")
  source_file <- NA_character_
}

file_base <- file.path(outdir, "result", "final_figures", "Figure2", "Figure2d_ATN_adjusted_attenuation")
save_gg_pdf_png_cm(p_atn_review, file_base, PANEL_REVIEW_W_CM, PANEL_REVIEW_H_CM)
figure_manifest <- add_manifest(figure_manifest, "Figure2d_ATN_adjusted_attenuation", file_base, source_file, "Main diagnostic DEP effects compared with AT(N)-adjusted effects; interpreted as sensitivity/attenuation, not causal independence.")

###############################################################################
# 08b_supplementary_AD_only_CDRSB_severity_alignment
###############################################################################

if (!is.null(cdr_compare_tbl) && nrow(cdr_compare_tbl) > 0) {
  cdr_compare_tbl <- rename_with_candidates(
    cdr_compare_tbl,
    list(
      AptName = c("AptName"),
      Protein_Name = c("Protein_Name", "ProteinLabel", "Target"),
      EntrezGeneSymbol = c("EntrezGeneSymbol", "Gene", "Symbol"),
      logFC_primary = c("logFC_primary", "main_logFC"),
      adj.P.Val_primary = c("adj.P.Val_primary", "main_adj.P.Val"),
      logFC_severity = c("beta_CDRSB_AD_only", "logFC_severity", "logFC_secondary", "logFC_cdrsb", "logFC_CDRSB"),
      adj.P.Val_severity = c("adj.P.Val_CDRSB_AD_only", "adj.P.Val_severity", "adj.P.Val_secondary", "adj.P.Val_cdrsb", "adj.P.Val_CDRSB")
    )
  )
  
  cdr_plot_tbl <- cdr_compare_tbl %>%
    dplyr::mutate(
      logFC_primary = safe_numeric_vec(logFC_primary),
      logFC_severity = safe_numeric_vec(logFC_severity),
      adj.P.Val_primary = safe_numeric_vec(adj.P.Val_primary),
      adj.P.Val_severity = safe_numeric_vec(adj.P.Val_severity),
      Protein_Label = dplyr::coalesce(Protein_Name, EntrezGeneSymbol, AptName),
      Direction = dplyr::case_when(
        logFC_primary > 0 ~ "Higher in AD",
        logFC_primary < 0 ~ "Lower in AD",
        TRUE ~ "Neutral"
      ),
      Direction = factor(Direction, levels = c("Higher in AD", "Lower in AD", "Neutral")),
      same_direction = sign(logFC_primary) == sign(logFC_severity),
      primary_fdr005 = adj.P.Val_primary < MAIN_FDR,
      severity_fdr005 = adj.P.Val_severity < MAIN_FDR,
      severity_nominal = adj.P.Val_severity < 0.05,
      alignment_class = dplyr::case_when(
        same_direction ~ "Same direction",
        !same_direction ~ "Opposite direction",
        TRUE ~ "Other"
      ),
      delta_slope = logFC_severity - logFC_primary
    ) %>%
    dplyr::filter(is.finite(logFC_primary), is.finite(logFC_severity))
  
  label_cdr <- cdr_plot_tbl %>%
    dplyr::filter(Protein_Label %in% label_targets | EntrezGeneSymbol %in% label_targets) %>%
    dplyr::arrange(adj.P.Val_primary) %>%
    dplyr::slice_head(n = 6)
  if (nrow(label_cdr) == 0) label_cdr <- cdr_plot_tbl %>% dplyr::arrange(adj.P.Val_primary) %>% dplyr::slice_head(n = 6)
  
  p_cdr_review <- make_density_scatter(
    cdr_plot_tbl,
    "logFC_primary", "logFC_severity",
    label_df = label_cdr,
    compact = FALSE,
    xlab = "Main-model log2FC (AD vs CN)",
    ylab = "AD-only CDR-SB severity slope"
  ) + labs(title = "Within-AD clinical severity does not recapitulate the diagnostic proteomic profile")
  
  p_cdr_compact <- make_density_scatter(
    cdr_plot_tbl,
    "logFC_primary", "logFC_severity",
    label_df = label_cdr %>% dplyr::slice_head(n = 4),
    compact = TRUE,
    xlab = "AD-vs-CN log2FC",
    ylab = "AD-only CDR-SB slope"
  )
  
  source_file <- file.path(outdir, "result", "final_source_data", "Figure2", "Supplementary_AD_only_CDRSB_severity_alignment_source_data.csv")
  readr::write_csv(
    cdr_plot_tbl %>% dplyr::select(Protein_Label, Protein_Name, EntrezGeneSymbol, AptName, logFC_primary, adj.P.Val_primary, logFC_severity, adj.P.Val_severity, same_direction, primary_fdr005, severity_fdr005, severity_nominal, alignment_class, delta_slope, Direction),
    source_file
  )
} else {
  p_cdr_review <- placeholder_plot("AD-only CDR-SB severity alignment table not found")
  p_cdr_compact <- placeholder_plot("AD-only CDR-SB severity alignment not available")
  source_file <- NA_character_
}

file_base <- file.path(outdir, "result", "final_figures", "Supplementary", "Supplementary_AD_only_CDRSB_severity_alignment")
save_gg_pdf_png_cm(p_cdr_review, file_base, PANEL_REVIEW_W_CM, PANEL_REVIEW_H_CM)
figure_manifest <- add_manifest(figure_manifest, "Supplementary_AD_only_CDRSB_severity_alignment", file_base, source_file, "Primary AD-vs-CN effects compared with AD-only CDR-SB severity slopes; not a CDR-SB-adjusted diagnostic model.")

###############################################################################
# 08c_supplementary_CDRSB_adjusted_AD_vs_CN_diagnostic_attenuation
###############################################################################

if (!is.null(cdrsb_adjusted_compare_tbl) && nrow(cdrsb_adjusted_compare_tbl) > 0) {
  cdrsb_adjusted_compare_tbl <- rename_with_candidates(
    cdrsb_adjusted_compare_tbl,
    list(
      AptName = c("AptName", "AptName_primary"),
      Protein_Name = c("Protein_Name", "ProteinLabel", "Target", "Protein_Name_primary"),
      EntrezGeneSymbol = c("EntrezGeneSymbol", "Gene", "Symbol"),
      logFC_primary = c("logFC_primary", "main_logFC"),
      adj.P.Val_primary = c("adj.P.Val_primary", "main_adj.P.Val"),
      logFC_CDRSB_adjusted = c("logFC_CDRSB_adjusted", "logFC_cdrsb_adjusted", "logFC_secondary", "logFC_adjusted"),
      adj.P.Val_CDRSB_adjusted = c("adj.P.Val_CDRSB_adjusted", "adj.P.Val_cdrsb_adjusted", "adj.P.Val_secondary", "adj.P.Val_adjusted")
    )
  )

  cdrsb_adjusted_plot_tbl <- cdrsb_adjusted_compare_tbl %>%
    dplyr::mutate(
      logFC_primary = safe_numeric_vec(logFC_primary),
      logFC_CDRSB_adjusted = safe_numeric_vec(logFC_CDRSB_adjusted),
      adj.P.Val_primary = safe_numeric_vec(adj.P.Val_primary),
      adj.P.Val_CDRSB_adjusted = safe_numeric_vec(adj.P.Val_CDRSB_adjusted),
      Protein_Label = dplyr::coalesce(Protein_Name, EntrezGeneSymbol, AptName),
      Direction = dplyr::case_when(
        logFC_primary > 0 ~ "Higher in AD",
        logFC_primary < 0 ~ "Lower in AD",
        TRUE ~ "Neutral"
      ),
      Direction = factor(Direction, levels = c("Higher in AD", "Lower in AD", "Neutral")),
      same_direction = sign(logFC_primary) == sign(logFC_CDRSB_adjusted),
      primary_fdr005 = adj.P.Val_primary < MAIN_FDR,
      cdrsb_adjusted_fdr005 = adj.P.Val_CDRSB_adjusted < MAIN_FDR,
      attenuation_absolute = abs(logFC_primary) - abs(logFC_CDRSB_adjusted),
      attenuation_ratio = dplyr::if_else(abs(logFC_primary) > 0, abs(logFC_CDRSB_adjusted) / abs(logFC_primary), NA_real_)
    ) %>%
    dplyr::filter(is.finite(logFC_primary), is.finite(logFC_CDRSB_adjusted))

  label_cdrsb_adjusted <- cdrsb_adjusted_plot_tbl %>%
    dplyr::filter(Protein_Label %in% label_targets | EntrezGeneSymbol %in% label_targets) %>%
    dplyr::arrange(adj.P.Val_primary) %>%
    dplyr::slice_head(n = 6)
  if (nrow(label_cdrsb_adjusted) == 0) {
    label_cdrsb_adjusted <- cdrsb_adjusted_plot_tbl %>%
      dplyr::arrange(adj.P.Val_primary) %>%
      dplyr::slice_head(n = 6)
  }

  p_cdrsb_adjusted_review <- make_direction_scatter(
    cdrsb_adjusted_plot_tbl,
    "logFC_primary", "logFC_CDRSB_adjusted",
    label_df = label_cdrsb_adjusted,
    compact = FALSE,
    xlab = "Main-model log2FC (AD vs CN)",
    ylab = "CDR-SB-adjusted log2FC (AD vs CN)"
  ) +
    ggplot2::labs(title = "Secondary CDR-SB-adjusted diagnostic attenuation analysis")

  source_file_cdrsb_adjusted <- file.path(outdir, "result", "final_source_data", "Supplementary", "Supplementary_CDRSB_adjusted_AD_vs_CN_attenuation_source_data.csv")
  readr::write_csv(
    cdrsb_adjusted_plot_tbl %>%
      dplyr::select(
        Protein_Label, Protein_Name, EntrezGeneSymbol, AptName,
        logFC_primary, adj.P.Val_primary,
        logFC_CDRSB_adjusted, adj.P.Val_CDRSB_adjusted,
        same_direction, primary_fdr005, cdrsb_adjusted_fdr005,
        attenuation_absolute, attenuation_ratio, Direction
      ),
    source_file_cdrsb_adjusted
  )
} else {
  p_cdrsb_adjusted_review <- placeholder_plot("CDR-SB-adjusted AD-vs-CN diagnostic attenuation table not found")
  source_file_cdrsb_adjusted <- NA_character_
}

file_base <- file.path(outdir, "result", "final_figures", "Supplementary", "Supplementary_CDRSB_adjusted_AD_vs_CN_attenuation")
save_gg_pdf_png_cm(p_cdrsb_adjusted_review, file_base, PANEL_REVIEW_W_CM, PANEL_REVIEW_H_CM)
figure_manifest <- add_manifest(figure_manifest, "Supplementary_CDRSB_adjusted_AD_vs_CN_attenuation", file_base, source_file_cdrsb_adjusted, "Secondary diagnostic attenuation model comparing primary AD-vs-CN effects with AD-vs-CN effects additionally adjusted for CDR-SB; not a within-AD severity model.")

###############################################################################
# 09_figure2f_loco_stability
###############################################################################

if (!is.null(main_vs_loco_mean) && nrow(main_vs_loco_mean) > 0) {
  main_vs_loco_mean <- rename_with_candidates(
    main_vs_loco_mean,
    list(
      AptName = c("AptName"),
      main_logFC = c("main_logFC", "logFC_primary", "logFC_main"),
      mean_loco_logFC = c("mean_loco_logFC", "mean_LOCO_logFC"),
      main_adj.P.Val = c("main_adj.P.Val", "adj.P.Val", "main_FDR")
    )
  )
  
  loco_plot_tbl <- main_vs_loco_mean %>%
    dplyr::left_join(
      DEP %>% dplyr::select(AptName, Protein_Name_DEP = Protein_Name, EntrezGeneSymbol_DEP = EntrezGeneSymbol, Protein_Label_DEP = Protein_Label),
      by = "AptName"
    ) %>%
    dplyr::mutate(
      main_logFC = safe_numeric_vec(main_logFC),
      mean_loco_logFC = safe_numeric_vec(mean_loco_logFC),
      main_adj.P.Val = safe_numeric_vec(main_adj.P.Val),
      Protein_Label = dplyr::coalesce(Protein_Label_DEP, Protein_Name_DEP, EntrezGeneSymbol_DEP, AptName),
      Direction = dplyr::case_when(
        main_logFC > 0 ~ "Higher in AD",
        main_logFC < 0 ~ "Lower in AD",
        TRUE ~ "Neutral"
      ),
      Direction = factor(Direction, levels = c("Higher in AD", "Lower in AD", "Neutral")),
      loco_delta = mean_loco_logFC - main_logFC
    ) %>%
    dplyr::filter(is.finite(main_logFC), is.finite(mean_loco_logFC))
  
  label_loco <- loco_plot_tbl %>%
    dplyr::filter(Protein_Label %in% label_targets | EntrezGeneSymbol_DEP %in% label_targets) %>%
    dplyr::arrange(main_adj.P.Val) %>%
    dplyr::slice_head(n = 6)
  if (nrow(label_loco) == 0) label_loco <- loco_plot_tbl %>% dplyr::arrange(main_adj.P.Val) %>% dplyr::slice_head(n = 6)
  
  p_loco_review <- make_direction_scatter(
    loco_plot_tbl,
    "main_logFC", "mean_loco_logFC",
    label_df = label_loco,
    compact = FALSE,
    xlab = "Main-model log2FC",
    ylab = "Mean leave-one-country-out log2FC"
  ) + labs(title = "Internal multicountry stability of AD proteomic effects")
  
  p_loco_compact <- make_direction_scatter(
    loco_plot_tbl,
    "main_logFC", "mean_loco_logFC",
    label_df = label_loco %>% dplyr::slice_head(n = 4),
    compact = TRUE,
    xlab = "Main-model log2FC",
    ylab = "Mean LOCO log2FC"
  ) +
    labs(title = "Country-exclusion stability of protein effects")
  
  source_file <- file.path(outdir, "result", "final_source_data", "Figure2", "Figure2f_LOCO_stability_source_data.csv")
  readr::write_csv(
    loco_plot_tbl %>% dplyr::select(Protein_Label, AptName, main_logFC, main_adj.P.Val, mean_loco_logFC, loco_delta, Direction),
    source_file
  )
} else {
  p_loco_review <- placeholder_plot("LOCO table not found")
  p_loco_compact <- placeholder_plot("LOCO not available")
  source_file <- NA_character_
}

file_base <- file.path(outdir, "result", "final_figures", "Figure2", "Figure2f_LOCO_stability")
save_gg_pdf_png_cm(p_loco_review, file_base, PANEL_REVIEW_W_CM, PANEL_REVIEW_H_CM)
figure_manifest <- add_manifest(figure_manifest, "Figure2f_LOCO_stability", file_base, source_file, "Protein-level main model versus mean leave-one-country-out effect.")


###############################################################################
# 09b_figure2e_contextual_country_exclusion_summary
###############################################################################

# This panel summarizes country-level exclusion effects using the LOCO summary
# table from Script 01. It is intentionally different from the protein-level
# LOCO scatter: here each point is one excluded country.

loco_summary_main <- read_csv_if_exists(
  file.path(outdir, "result", "06_robustness", "country_loco", "tables", "LOCO_summary_metrics.csv"),
  show_message = FALSE
)

if (!is.null(loco_summary_main) && nrow(loco_summary_main) > 0) {
  loco_summary_main <- rename_with_candidates(
    loco_summary_main,
    list(
      excluded_country = c("excluded_country", "country", "Country"),
      logFC_correlation = c("logFC_correlation", "correlation", "r"),
      direction_consistency = c("direction_consistency", "direction_consistency_all", "direction_consistency_rate")
    )
  )

  if (all(c("excluded_country", "logFC_correlation") %in% names(loco_summary_main))) {
    loco_context_tbl <- loco_summary_main %>%
      dplyr::mutate(
        excluded_country = as.character(excluded_country),
        logFC_correlation = safe_numeric_vec(logFC_correlation),
        direction_consistency = if ("direction_consistency" %in% names(.)) safe_numeric_vec(direction_consistency) else NA_real_,
        stability_class = dplyr::if_else(logFC_correlation >= 0.85, "High concordance", "Context-sensitive"),
        excluded_country = forcats::fct_reorder(excluded_country, logFC_correlation, .desc = FALSE),
        label_corr = sprintf("%.2f", logFC_correlation)
      ) %>%
      dplyr::filter(is.finite(logFC_correlation))

    make_loco_context_plot <- function(plot_tbl, compact = FALSE) {
      theme_use <- if (compact) theme_panel_composite() else theme_panel_review()
      point_size <- if (compact) 1.6 else 2.4
      label_size <- if (compact) 2.2 else 3.0
      axis_size <- if (compact) 6.8 else 9.0
      ref_x <- 0.85

      ggplot2::ggplot(plot_tbl, ggplot2::aes(y = excluded_country, x = logFC_correlation)) +
        ggplot2::geom_segment(
          ggplot2::aes(x = 0, xend = logFC_correlation, yend = excluded_country),
          colour = "grey78", linewidth = if (compact) 0.42 else 0.60
        ) +
        ggplot2::geom_vline(xintercept = ref_x, linetype = 2, linewidth = 0.35, colour = "grey45") +
        ggplot2::geom_point(
          ggplot2::aes(fill = stability_class),
          shape = 21, size = point_size, colour = "black", stroke = if (compact) 0.18 else 0.25
        ) +
        ggplot2::geom_text(
          ggplot2::aes(label = label_corr),
          hjust = -0.35, size = label_size, colour = COL_TEXT
        ) +
        ggplot2::scale_fill_manual(
          values = c("High concordance" = COL_DOWN, "Context-sensitive" = COL_UP),
          breaks = c("High concordance", "Context-sensitive"),
          guide = "none"
        ) +
        ggplot2::coord_cartesian(xlim = c(0, 1.03), expand = FALSE, clip = "off") +
        ggplot2::scale_x_continuous(breaks = c(0, 0.25, 0.50, 0.75, 1.00)) +
        ggplot2::labs(
          title = if (!compact) "Contextual stability across country-specific exclusions" else NULL,
          x = "Correlation with primary-model log2FC",
          y = "Excluded country"
        ) +
        theme_use +
        ggplot2::theme(
          axis.text.y = ggplot2::element_text(size = axis_size, colour = COL_TEXT),
          axis.text.x = ggplot2::element_text(size = axis_size, colour = COL_TEXT),
          plot.margin = ggplot2::margin(4, 10, 4, 4)
        )
    }

    p_loco_context_review <- make_loco_context_plot(loco_context_tbl, compact = FALSE)
    p_loco_context_compact <- make_loco_context_plot(loco_context_tbl, compact = TRUE)

    source_file <- file.path(outdir, "result", "final_source_data", "Figure2", "Figure2e_contextual_country_exclusion_source_data.csv")
    readr::write_csv(loco_context_tbl, source_file)

  } else {
    p_loco_context_review <- placeholder_plot("LOCO summary table missing required columns")
    p_loco_context_compact <- placeholder_plot("LOCO summary unavailable")
    source_file <- NA_character_
  }
} else {
  p_loco_context_review <- placeholder_plot("LOCO summary table not found")
  p_loco_context_compact <- placeholder_plot("LOCO summary not available")
  source_file <- NA_character_
}

file_base <- file.path(outdir, "result", "final_figures", "Figure2", "Figure2e_contextual_country_exclusion")
save_gg_pdf_png_cm(p_loco_context_review, file_base, PANEL_REVIEW_W_CM, PANEL_REVIEW_H_CM)
figure_manifest <- add_manifest(
  figure_manifest,
  "Figure2e_contextual_country_exclusion",
  file_base,
  source_file,
  "Country-level LOCO summary showing effect-size correlation after excluding each country."
)

###############################################################################
# 10_figure2c_heatmap_full_and_compact
###############################################################################

if (length(trait_vars) == 0) {
  warning("No trait variables detected for heatmap. Figure2c will be replaced by a placeholder.")
  p_heatmap_placeholder <- placeholder_plot("No clinical/AT(N) traits detected")
  heatmap_grob_compact <- patchwork::wrap_elements(p_heatmap_placeholder)
  heatmap_source_file <- NA_character_
} else {
  main_tbl <- DEP %>%
    dplyr::filter(Direction %in% c("Higher in AD", "Lower in AD")) %>%
    dplyr::transmute(
      AptName,
      Protein_Name,
      EntrezGeneSymbol,
      Protein_Label,
      main_logFC = logFC,
      main_adjP = adj.P.Val,
      main_type = as.character(Direction)
    )
  
  robustness_heat_tbl <- main_tbl
  
  if (!is.null(apoe_compare_tbl) && nrow(apoe_compare_tbl) > 0) {
    apoe_compare_tbl <- rename_with_candidates(
      apoe_compare_tbl,
      list(
        AptName = c("AptName"),
        logFC_primary = c("logFC_primary", "main_logFC"),
        adj.P.Val_primary = c("adj.P.Val_primary", "main_adj.P.Val"),
        logFC_secondary = c("logFC_secondary", "logFC_apoe"),
        adj.P.Val_secondary = c("adj.P.Val_secondary", "adj.P.Val_apoe")
      )
    )
    apoe_ann <- apoe_compare_tbl %>%
      dplyr::transmute(
        AptName,
        apoe_logFC = safe_numeric_vec(logFC_secondary),
        apoe_adjP = safe_numeric_vec(adj.P.Val_secondary),
        apoe_status = dplyr::case_when(
          is.na(apoe_logFC) ~ "Missing",
          sign(safe_numeric_vec(logFC_primary)) == sign(apoe_logFC) & dplyr::coalesce(apoe_adjP, 1) < MAIN_FDR ~ "Preserved",
          sign(safe_numeric_vec(logFC_primary)) == sign(apoe_logFC) ~ "Direction only",
          TRUE ~ "Not preserved"
        )
      )
    robustness_heat_tbl <- robustness_heat_tbl %>% dplyr::left_join(apoe_ann, by = "AptName")
  } else {
    robustness_heat_tbl$apoe_status <- "Missing"
  }
  
  if (!is.null(cdr_compare_tbl) && nrow(cdr_compare_tbl) > 0) {
    cdr_compare_tbl <- rename_with_candidates(
      cdr_compare_tbl,
      list(
        AptName = c("AptName"),
        logFC_primary = c("logFC_primary", "main_logFC"),
        logFC_severity = c("beta_CDRSB_AD_only", "logFC_severity", "logFC_secondary", "logFC_cdrsb"),
        adj.P.Val_severity = c("adj.P.Val_CDRSB_AD_only", "adj.P.Val_severity", "adj.P.Val_secondary", "adj.P.Val_cdrsb")
      )
    )
    cdr_ann <- cdr_compare_tbl %>%
      dplyr::transmute(
        AptName,
        cdr_logFC = safe_numeric_vec(logFC_severity),
        cdr_adjP = safe_numeric_vec(adj.P.Val_severity),
        cdr_status = dplyr::case_when(
          is.na(cdr_logFC) ~ "Missing",
          sign(safe_numeric_vec(logFC_primary)) == sign(cdr_logFC) ~ "Same direction",
          TRUE ~ "Opposite direction"
        ),
        cdr_abs_ratio = abs(cdr_logFC) / pmax(abs(safe_numeric_vec(logFC_primary)), 1e-9)
      )
    robustness_heat_tbl <- robustness_heat_tbl %>% dplyr::left_join(cdr_ann, by = "AptName")
  } else {
    robustness_heat_tbl$cdr_status <- "Missing"
    robustness_heat_tbl$cdr_abs_ratio <- NA_real_
  }
  
  curated_priority <- c("SPC25", "LRRN1", "SMOC1", "CPLX2", "ACHE", "IGFBP2", "C3", "AMBRA1", "GDI1", "PHGDH", "PUM2", "CAPN1")
  
  selection_tbl <- robustness_heat_tbl %>%
    dplyr::mutate(
      apoe_status = dplyr::coalesce(apoe_status, "Missing"),
      cdr_status = dplyr::coalesce(cdr_status, "Missing"),
      priority_anchor = Protein_Label %in% curated_priority | EntrezGeneSymbol %in% curated_priority,
      robustness_score =
        4 * (apoe_status == "Preserved") +
        2 * (apoe_status == "Direction only") +
        2 * priority_anchor,
      severity_alignment_score = dplyr::case_when(
        cdr_status == "Same direction" ~ 0.5,
        TRUE ~ 0
      ),
      selection_score = robustness_score + severity_alignment_score
    )
  
  select_balanced_proteins <- function(tbl, n_total) {
    n_each <- floor(n_total / 2)
    up <- tbl %>%
      dplyr::filter(main_type == "Higher in AD") %>%
      dplyr::arrange(dplyr::desc(selection_score), main_adjP, dplyr::desc(abs(main_logFC))) %>%
      dplyr::slice_head(n = n_each)
    down <- tbl %>%
      dplyr::filter(main_type == "Lower in AD") %>%
      dplyr::arrange(dplyr::desc(selection_score), main_adjP, dplyr::desc(abs(main_logFC))) %>%
      dplyr::slice_head(n = n_each)
    out <- dplyr::bind_rows(up, down) %>% dplyr::distinct(AptName, .keep_all = TRUE)
    if (nrow(out) < n_total) {
      filler <- tbl %>%
        dplyr::filter(!AptName %in% out$AptName) %>%
        dplyr::arrange(dplyr::desc(selection_score), main_adjP, dplyr::desc(abs(main_logFC))) %>%
        dplyr::slice_head(n = n_total - nrow(out))
      out <- dplyr::bind_rows(out, filler)
    }
    out %>%
      dplyr::mutate(direction_order = ifelse(main_type == "Higher in AD", 1, 2)) %>%
      dplyr::arrange(direction_order, main_adjP, dplyr::desc(selection_score)) %>%
      dplyr::slice_head(n = n_total)
  }
  
  heatmap_map_full <- select_balanced_proteins(selection_tbl, n_total = 12)
  heatmap_map_compact <- heatmap_map_full
  
  make_heatmap_object <- function(protein_map, compact = FALSE, show_legends = TRUE) {
    proteins_present <- intersect(protein_map$AptName, names(normalized_expr_plot))
    traits_present <- intersect(trait_vars, names(normalized_expr_plot))
    
    hm_df <- normalized_expr_plot %>%
      dplyr::select(dplyr::all_of(c(traits_present, proteins_present))) %>%
      dplyr::mutate(dplyr::across(dplyr::everything(), safe_numeric_vec))
    
    valid_traits <- traits_present[
      vapply(traits_present, function(v) {
        x <- hm_df[[v]]
        sum(is.finite(x)) >= 30 && stats::sd(x[is.finite(x)], na.rm = TRUE) > 0
      }, logical(1))
    ]
    valid_proteins <- proteins_present[
      vapply(proteins_present, function(v) {
        x <- hm_df[[v]]
        sum(is.finite(x)) >= 30 && stats::sd(x[is.finite(x)], na.rm = TRUE) > 0
      }, logical(1))
    ]
    
    if (length(valid_traits) == 0 || length(valid_proteins) == 0) {
      stop("No valid trait/protein pairs for heatmap after missingness and variance filters.")
    }
    
    cor_long <- compute_pairwise_spearman_local(hm_df, valid_proteins, valid_traits, min_n = 30) %>%
      dplyr::left_join(
        protein_map %>% dplyr::select(AptName, Protein_Label, EntrezGeneSymbol, main_logFC, main_adjP, main_type, apoe_status, cdr_status, selection_score),
        by = "AptName"
      ) %>%
      dplyr::mutate(
        trait_display = dplyr::recode(trait, !!!trait_display),
        protein_display = dplyr::coalesce(Protein_Label, EntrezGeneSymbol, AptName)
      )
    
    trait_display_order <- unique(dplyr::recode(valid_traits, !!!trait_display))
    
    cor_wide <- cor_long %>%
      dplyr::select(AptName, trait_display, rho) %>%
      dplyr::distinct() %>%
      tidyr::pivot_wider(names_from = trait_display, values_from = rho)
    cor_mat <- as.data.frame(cor_wide)
    rownames(cor_mat) <- cor_mat$AptName
    cor_mat$AptName <- NULL
    cor_mat <- as.matrix(cor_mat[, intersect(trait_display_order, colnames(cor_mat)), drop = FALSE])
    
    q_wide <- cor_long %>%
      dplyr::select(AptName, trait_display, q_value_bh) %>%
      dplyr::distinct() %>%
      tidyr::pivot_wider(names_from = trait_display, values_from = q_value_bh)
    q_mat <- as.data.frame(q_wide)
    rownames(q_mat) <- q_mat$AptName
    q_mat$AptName <- NULL
    q_mat <- as.matrix(q_mat[, colnames(cor_mat), drop = FALSE])
    
    protein_order <- protein_map %>%
      dplyr::filter(AptName %in% rownames(cor_mat)) %>%
      dplyr::mutate(direction_order = ifelse(main_type == "Higher in AD", 1, 2)) %>%
      dplyr::arrange(direction_order, main_adjP, dplyr::desc(selection_score)) %>%
      dplyr::pull(AptName)
    protein_order <- intersect(protein_order, rownames(cor_mat))
    cor_mat <- cor_mat[protein_order, , drop = FALSE]
    q_mat <- q_mat[protein_order, , drop = FALSE]
    
    label_map <- protein_map %>% dplyr::select(AptName, Protein_Label) %>% dplyr::distinct(AptName, .keep_all = TRUE)
    pretty_labels <- label_map$Protein_Label[match(protein_order, label_map$AptName)]
    pretty_labels[is.na(pretty_labels)] <- protein_order[is.na(pretty_labels)]
    pretty_labels <- make.unique(pretty_labels)
    rownames(cor_mat) <- pretty_labels
    rownames(q_mat) <- pretty_labels
    
    cell_labels <- matrix("", nrow = nrow(q_mat), ncol = ncol(q_mat), dimnames = dimnames(q_mat))
    cell_labels[!is.na(q_mat) & q_mat < 0.05] <- "*"
    cell_labels[!is.na(q_mat) & q_mat < 0.01] <- "**"
    
    row_info <- protein_map %>%
      dplyr::filter(AptName %in% protein_order) %>%
      dplyr::mutate(
        Direction = ifelse(main_logFC > 0, "Higher in AD", "Lower in AD"),
        APOE = factor(dplyr::coalesce(apoe_status, "Missing"), levels = c("Preserved", "Direction only", "Not preserved", "Missing")),
        Severity = factor(dplyr::coalesce(cdr_status, "Missing"), levels = c("Same direction", "Opposite direction", "Missing"))
      ) %>% as.data.frame()
    rownames(row_info) <- row_info$AptName
    row_info <- row_info[protein_order, , drop = FALSE]
    rownames(row_info) <- pretty_labels
    
    status_cols <- c("Preserved" = COL_PRESERVED, "Direction only" = COL_DIRONLY, "Not preserved" = COL_NOTPRES, "Missing" = COL_MISSING)
    severity_cols <- c("Same direction" = COL_DIRONLY, "Opposite direction" = COL_NOTPRES, "Missing" = COL_MISSING)
    
    cell_fun_sig <- function(j, i, x, y, width, height, fill) {
      lab <- cell_labels[i, j]
      if (!is.na(lab) && lab != "") {
        grid::grid.text(lab, x = x, y = y, gp = grid::gpar(fontsize = 6, fontfamily = FONT_FAMILY, col = "black", fontface = "plain"))
      }
    }
    
    left_annot <- ComplexHeatmap::rowAnnotation(
      Direction = row_info$Direction,
      APOE = row_info$APOE,
      `AD-only severity` = row_info$Severity,
      col = list(
        Direction = c("Higher in AD" = COL_UP, "Lower in AD" = COL_DOWN),
        APOE = status_cols,
        `AD-only severity` = severity_cols
      ),
      gp = grid::gpar(col = "grey90", lwd = 0.35),
      border = TRUE,
      annotation_name_gp = grid::gpar(fontsize = TEXT_SIZE, fontfamily = FONT_FAMILY, fontface = "plain"),
      annotation_name_rot = 90,
      annotation_name_side = "bottom",
      width = grid::unit(c(1.6, 1.6, 1.6), "mm"),
      show_legend = show_legends
    )
    
    ht <- ComplexHeatmap::Heatmap(
      cor_mat,
      name = "Spearman rho",
      col = circlize::colorRamp2(c(-0.5, 0, 0.5), c(COL_DOWN, "white", COL_UP)),
      rect_gp = grid::gpar(col = "grey92", lwd = 0.25),
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_row_dend = FALSE,
      show_column_dend = FALSE,
      row_names_side = "right",
      row_names_gp = grid::gpar(fontsize = TEXT_SIZE, fontfamily = FONT_FAMILY, fontface = "plain"),
      row_names_max_width = grid::unit(2.6, "cm"),
      column_names_gp = grid::gpar(fontsize = TEXT_SIZE, fontfamily = FONT_FAMILY, fontface = "plain"),
      column_names_rot = 90,
      column_names_side = "bottom",
      column_names_centered = TRUE,
      column_names_max_height = grid::unit(2.1, "cm"),
      cell_fun = cell_fun_sig,
      show_heatmap_legend = show_legends,
      heatmap_legend_param = list(
        title = "Spearman rho",
        title_gp = grid::gpar(fontsize = TEXT_SIZE, fontfamily = FONT_FAMILY, fontface = "plain"),
        labels_gp = grid::gpar(fontsize = TEXT_SIZE, fontfamily = FONT_FAMILY, fontface = "plain"),
        legend_height = grid::unit(2.2, "cm")
      )
    )
    
    list(ht = left_annot + ht, cor_long = cor_long, protein_map = protein_map)
  }
  
  heat_full <- make_heatmap_object(heatmap_map_full, compact = FALSE, show_legends = TRUE)
  heat_compact <- make_heatmap_object(heatmap_map_compact, compact = TRUE, show_legends = TRUE)
  
  save_complex_heatmap_pdf_png(
    heat_full$ht,
    file.path(outdir, "result", "final_figures", "Figure2", "Figure2c_representative_heatmap_full.pdf"),
    file.path(outdir, "result", "final_figures", "Figure2", "Figure2c_representative_heatmap_full.png"),
    width_cm = 18,
    height_cm = 14,
    dpi = 600,
    padding = grid::unit(c(5, 5, 5, 5), "mm")
  )
  
  heatmap_source_file <- file.path(outdir, "result", "final_source_data", "Figure2", "Figure2c_representative_heatmap_source_data.csv")
  readr::write_csv(
    heat_full$cor_long %>%
      dplyr::select(protein_display, AptName, trait, trait_display, n, rho, p_value, q_value_bh, main_logFC, main_adjP, main_type, apoe_status, cdr_status),
    heatmap_source_file
  )
  
  heatmap_selection_file <- file.path(outdir, "result", "final_source_data", "Figure2", "Figure2c_representative_heatmap_selected_proteins.csv")
  readr::write_csv(heat_full$protein_map, heatmap_selection_file)
  
  heatmap_grob_compact <- patchwork::wrap_elements(
    full = capture_complex_heatmap(
      heat_compact$ht,
      heatmap_legend_side = "right",
      annotation_legend_side = "right",
      merge_legends = FALSE,
      padding = grid::unit(c(0.5, 0.5, 0.5, 0.5), "mm")
    )
  )
  
  figure_manifest <- add_manifest(
    figure_manifest,
    "Figure2c_representative_heatmap_full",
    file.path(outdir, "result", "final_figures", "Figure2", "Figure2c_representative_heatmap_full"),
    heatmap_source_file,
    "Full 12-protein heatmap saved separately and used in composite."
  )
}

###############################################################################
# 11_composite_figure2
###############################################################################

figure2_composite <- (
  p_volcano_compact + p_gsea_compact +
    patchwork::plot_layout(ncol = 2, widths = c(1.05, 0.95))
) / (
  patchwork::wrap_elements(heatmap_grob_compact)
) / (
  p_atn_compact + p_loco_context_compact + p_loco_compact +
    patchwork::plot_layout(ncol = 3, widths = c(1.02, 0.82, 1.02))
) +
  patchwork::plot_layout(heights = c(1.05, 0.66, 1.02)) +
  patchwork::plot_annotation(
    tag_levels = "a",
    theme = ggplot2::theme(
      text = ggplot2::element_text(family = FONT_FAMILY),
      plot.tag = ggplot2::element_text(family = FONT_FAMILY, size = TAG_SIZE, face = "bold", colour = COL_TEXT),
      plot.margin = ggplot2::margin(1.5, 1.5, 1.5, 1.5)
    )
  )

file_base <- file.path(outdir, "result", "final_figures", "Figure2", "Figure2_composite_18cm_optimized")
save_gg_pdf_png_cm(figure2_composite, file_base, width_cm = COMPOSITE_WIDTH_CM, height_cm = COMPOSITE_HEIGHT_CM, dpi = 600, limitsize = FALSE)
figure_manifest <- add_manifest(figure_manifest, "Figure2_composite_18cm_optimized", file_base, "See individual panel source-data files", "Composite optimized for final 18-cm width.")

###############################################################################
# 12_supplementary_pca_qc
###############################################################################

make_pca_plot <- function(data, seq_cols, color_var, title_text) {
  if (is.na(color_var) || !color_var %in% names(data)) {
    return(placeholder_plot(paste("Missing variable:", color_var)))
  }
  
  tmp <- data %>% dplyr::select(dplyr::all_of(c(seq_cols, color_var))) %>% dplyr::filter(!is.na(.data[[color_var]]))
  if (nrow(tmp) < 5) return(placeholder_plot(paste("Too few samples for:", title_text)))
  
  xmat <- tmp %>% dplyr::select(dplyr::all_of(seq_cols)) %>% as.matrix()
  storage.mode(xmat) <- "numeric"
  keep <- apply(xmat, 2, function(z) stats::sd(z, na.rm = TRUE) > 0)
  xmat <- xmat[, keep, drop = FALSE]
  if (ncol(xmat) < 2) return(placeholder_plot(paste("Not enough proteins for:", title_text)))
  
  for (j in seq_len(ncol(xmat))) {
    miss <- !is.finite(xmat[, j])
    if (any(miss)) xmat[miss, j] <- stats::median(xmat[, j], na.rm = TRUE)
  }
  
  pca <- stats::prcomp(xmat, center = TRUE, scale. = FALSE)
  pca_df <- tibble::tibble(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = as.factor(tmp[[color_var]]))
  ve <- summary(pca)$importance[2, 1:2] * 100
  
  ggplot(pca_df, aes(PC1, PC2, colour = group)) +
    geom_point(shape = 16, size = 1.7, alpha = 0.85) +
    labs(title = title_text, x = paste0("PC1 (", round(ve[1], 1), "%)"), y = paste0("PC2 (", round(ve[2], 1), "%)"), colour = NULL) +
    theme_panel_review(base_size = 10) +
    theme(legend.position = "right")
}

site_var    <- c("Site", "site", "Center", "center", "Cohort", "cohort")[c("Site", "site", "Center", "center", "Cohort", "cohort") %in% names(normalized_expr_plot)][1]
plate_var   <- c("PlateId", "Plate", "plate")[c("PlateId", "Plate", "plate") %in% names(normalized_expr_plot)][1]
country_var <- c("Country", "country")[c("Country", "country") %in% names(normalized_expr_plot)][1]
group_var   <- c("SampleGroup", "Group", "Diagnosis")[c("SampleGroup", "Group", "Diagnosis") %in% names(normalized_expr_plot)][1]
sex_var     <- c("Sex", "sex")[c("Sex", "sex") %in% names(normalized_expr_plot)][1]
age_var     <- c("Age", "age")[c("Age", "age") %in% names(normalized_expr_plot)][1]

p_s1a <- if (!is.na(plate_var))   make_pca_plot(normalized_expr_plot, seq_cols_plot, plate_var, "PCA colored by plate") else placeholder_plot("Plate variable not available")
p_s1b <- if (!is.na(country_var)) make_pca_plot(normalized_expr_plot, seq_cols_plot, country_var, "PCA colored by country") else placeholder_plot("Country variable not available")
p_s1c <- if (!is.na(site_var))    make_pca_plot(normalized_expr_plot, seq_cols_plot, site_var, "PCA colored by site") else placeholder_plot("Site variable not available")
p_s1d <- if (!is.na(group_var))   make_pca_plot(normalized_expr_plot, seq_cols_plot, group_var, "PCA colored by diagnosis") else placeholder_plot("Diagnosis variable not available")
p_s1e <- if (!is.na(sex_var))     make_pca_plot(normalized_expr_plot, seq_cols_plot, sex_var, "PCA colored by sex") else placeholder_plot("Sex variable not available")
p_s1f <- if (!is.na(age_var))     make_pca_plot(normalized_expr_plot, seq_cols_plot, age_var, "PCA colored by age") else placeholder_plot("Age variable not available")

supp_s1 <- (p_s1a + p_s1b + p_s1c) / (p_s1d + p_s1e + p_s1f) +
  patchwork::plot_annotation(tag_levels = "a", theme = ggplot2::theme(plot.tag = ggplot2::element_text(size = 12, face = "bold")))

file_base <- file.path(outdir, "result", "final_figures", "Supplementary", "Supplementary_Fig_S1_PCA_QC")
save_gg_pdf_png_cm(supp_s1, file_base, width_cm = 18, height_cm = 16, dpi = 600)
figure_manifest <- add_manifest(figure_manifest, "Supplementary_Fig_S1_PCA_QC", file_base, NA_character_, "PCA QC panels colored by available technical and biological variables.")

###############################################################################
# 13_supplementary_sensitivity_plots
###############################################################################

if (!is.null(apoe_compare_tbl) && nrow(apoe_compare_tbl) > 0) {
  apoe_plot_tbl <- apoe_compare_tbl %>%
    rename_with_candidates(
      list(
        AptName = c("AptName"),
        Protein_Name = c("Protein_Name", "ProteinLabel", "Target"),
        EntrezGeneSymbol = c("EntrezGeneSymbol", "Gene", "Symbol"),
        logFC_primary = c("logFC_primary", "main_logFC"),
        adj.P.Val_primary = c("adj.P.Val_primary", "main_adj.P.Val"),
        logFC_secondary = c("logFC_secondary", "logFC_apoe"),
        adj.P.Val_secondary = c("adj.P.Val_secondary", "adj.P.Val_apoe")
      )
    ) %>%
    dplyr::mutate(
      logFC_primary = safe_numeric_vec(logFC_primary),
      logFC_secondary = safe_numeric_vec(logFC_secondary),
      Protein_Label = dplyr::coalesce(Protein_Name, EntrezGeneSymbol, AptName),
      Direction = dplyr::case_when(logFC_primary > 0 ~ "Higher in AD", logFC_primary < 0 ~ "Lower in AD", TRUE ~ "Neutral"),
      Direction = factor(Direction, levels = c("Higher in AD", "Lower in AD", "Neutral"))
    ) %>%
    dplyr::filter(is.finite(logFC_primary), is.finite(logFC_secondary))
  
  label_apoe <- apoe_plot_tbl %>% dplyr::filter(Protein_Label %in% label_targets | EntrezGeneSymbol %in% label_targets) %>% dplyr::slice_head(n = 6)
  if (nrow(label_apoe) == 0) label_apoe <- apoe_plot_tbl %>% dplyr::arrange(adj.P.Val_primary) %>% dplyr::slice_head(n = 6)
  
  p_apoe_supp <- make_direction_scatter(apoe_plot_tbl, "logFC_primary", "logFC_secondary", label_df = label_apoe, compact = FALSE, xlab = "Main-model log2FC", ylab = "APOE-adjusted log2FC") +
    labs(title = "APOE sensitivity analysis")
  
  file_base <- file.path(outdir, "result", "final_figures", "Supplementary", "Supplementary_Fig_S2_APOE_sensitivity")
  save_gg_pdf_png_cm(p_apoe_supp, file_base, width_cm = 18, height_cm = 13, dpi = 600)
  source_file <- file.path(outdir, "result", "final_source_data", "Supplementary", "Supplementary_Fig_S2_APOE_sensitivity_source_data.csv")
  readr::write_csv(apoe_plot_tbl, source_file)
  figure_manifest <- add_manifest(figure_manifest, "Supplementary_Fig_S2_APOE_sensitivity", file_base, source_file, "Main DEP versus APOE-adjusted DEP.")
}

if (!is.null(atn_compare_tbl) && nrow(atn_compare_tbl) > 0) {
  atn_plot_tbl <- atn_compare_tbl %>%
    rename_with_candidates(
      list(
        AptName = c("AptName"),
        Protein_Name = c("Protein_Name", "ProteinLabel", "Target"),
        EntrezGeneSymbol = c("EntrezGeneSymbol", "Gene", "Symbol"),
        logFC_primary = c("logFC_primary", "main_logFC"),
        adj.P.Val_primary = c("adj.P.Val_primary", "main_adj.P.Val"),
        logFC_atn = c("logFC_atn", "logFC_ATN", "logFC_secondary"),
        adj.P.Val_atn = c("adj.P.Val_atn", "adj.P.Val_ATN", "adj.P.Val_secondary")
      )
    ) %>%
    dplyr::mutate(
      logFC_primary = safe_numeric_vec(logFC_primary),
      logFC_atn = safe_numeric_vec(logFC_atn),
      Protein_Label = dplyr::coalesce(Protein_Name, EntrezGeneSymbol, AptName),
      Direction = dplyr::case_when(logFC_primary > 0 ~ "Higher in AD", logFC_primary < 0 ~ "Lower in AD", TRUE ~ "Neutral"),
      Direction = factor(Direction, levels = c("Higher in AD", "Lower in AD", "Neutral"))
    ) %>%
    dplyr::filter(is.finite(logFC_primary), is.finite(logFC_atn))
  
  label_atn <- atn_plot_tbl %>% dplyr::filter(Protein_Label %in% label_targets | EntrezGeneSymbol %in% label_targets) %>% dplyr::slice_head(n = 6)
  if (nrow(label_atn) == 0) label_atn <- atn_plot_tbl %>% dplyr::arrange(adj.P.Val_primary) %>% dplyr::slice_head(n = 6)
  
  p_atn_supp <- make_direction_scatter(atn_plot_tbl, "logFC_primary", "logFC_atn", label_df = label_atn, compact = FALSE, xlab = "Main-model log2FC", ylab = "AT(N)-adjusted log2FC") +
    labs(title = "AT(N)-adjusted sensitivity analysis")
  
  file_base <- file.path(outdir, "result", "final_figures", "Supplementary", "Supplementary_Fig_S3_ATN_adjusted_sensitivity")
  save_gg_pdf_png_cm(p_atn_supp, file_base, width_cm = 18, height_cm = 13, dpi = 600)
  source_file <- file.path(outdir, "result", "final_source_data", "Supplementary", "Supplementary_Fig_S3_ATN_adjusted_sensitivity_source_data.csv")
  readr::write_csv(atn_plot_tbl, source_file)
  figure_manifest <- add_manifest(figure_manifest, "Supplementary_Fig_S3_ATN_adjusted_sensitivity", file_base, source_file, "Main DEP versus AT(N)-adjusted DEP; covariates include p-tau217, NfL and Aβ42/40 when available.")
}

if (!is.null(vascular_compare_tbl) && nrow(vascular_compare_tbl) > 0) {
  vascular_plot_tbl <- vascular_compare_tbl %>%
    rename_with_candidates(
      list(
        AptName = c("AptName"),
        Protein_Name = c("Protein_Name", "ProteinLabel", "Target"),
        EntrezGeneSymbol = c("EntrezGeneSymbol", "Gene", "Symbol"),
        logFC_primary = c("logFC_primary", "main_logFC"),
        adj.P.Val_primary = c("adj.P.Val_primary", "main_adj.P.Val"),
        logFC_vascular = c("logFC_vascular", "logFC_secondary"),
        adj.P.Val_vascular = c("adj.P.Val_vascular", "adj.P.Val_secondary")
      )
    ) %>%
    dplyr::mutate(
      logFC_primary = safe_numeric_vec(logFC_primary),
      logFC_vascular = safe_numeric_vec(logFC_vascular),
      Protein_Label = dplyr::coalesce(Protein_Name, EntrezGeneSymbol, AptName),
      Direction = dplyr::case_when(logFC_primary > 0 ~ "Higher in AD", logFC_primary < 0 ~ "Lower in AD", TRUE ~ "Neutral"),
      Direction = factor(Direction, levels = c("Higher in AD", "Lower in AD", "Neutral"))
    ) %>%
    dplyr::filter(is.finite(logFC_primary), is.finite(logFC_vascular))
  
  label_vascular <- vascular_plot_tbl %>% dplyr::filter(Protein_Label %in% label_targets | EntrezGeneSymbol %in% label_targets) %>% dplyr::slice_head(n = 6)
  if (nrow(label_vascular) == 0) label_vascular <- vascular_plot_tbl %>% dplyr::arrange(adj.P.Val_primary) %>% dplyr::slice_head(n = 6)
  
  p_vascular_supp <- make_direction_scatter(vascular_plot_tbl, "logFC_primary", "logFC_vascular", label_df = label_vascular, compact = FALSE, xlab = "Main-model log2FC", ylab = "Vascular/metabolic-adjusted log2FC") +
    labs(title = "Vascular/metabolic sensitivity analysis")
  
  file_base <- file.path(outdir, "result", "final_figures", "Supplementary", "Supplementary_Fig_S4_vascular_metabolic_sensitivity")
  save_gg_pdf_png_cm(p_vascular_supp, file_base, width_cm = 18, height_cm = 13, dpi = 600)
  source_file <- file.path(outdir, "result", "final_source_data", "Supplementary", "Supplementary_Fig_S4_vascular_metabolic_sensitivity_source_data.csv")
  readr::write_csv(vascular_plot_tbl, source_file)
  figure_manifest <- add_manifest(figure_manifest, "Supplementary_Fig_S4_vascular_metabolic_sensitivity", file_base, source_file, "Main DEP versus vascular/metabolic-adjusted DEP.")
}

if (exists("p_cdr_review")) {
  file_base <- file.path(outdir, "result", "final_figures", "Supplementary", "Supplementary_Fig_S5_AD_only_CDRSB_severity_alignment")
  save_gg_pdf_png_cm(p_cdr_review, file_base, width_cm = 18, height_cm = 13, dpi = 600)
  figure_manifest <- add_manifest(figure_manifest, "Supplementary_Fig_S5_AD_only_CDRSB_severity_alignment", file_base, file.path(outdir, "result", "final_source_data", "Figure2", "Supplementary_AD_only_CDRSB_severity_alignment_source_data.csv"), "Full-size version of the AD-only CDR-SB severity alignment panel.")
}

###############################################################################
# 14_supplementary_country_robustness
###############################################################################

loco_summary <- read_csv_if_exists(file.path(outdir, "result", "06_robustness", "country_loco", "tables", "LOCO_summary_metrics.csv"), show_message = FALSE)

if (!is.null(loco_summary) && nrow(loco_summary) > 0) {
  loco_summary <- rename_with_candidates(loco_summary, list(excluded_country = c("excluded_country", "country", "Country"), logFC_correlation = c("logFC_correlation", "correlation", "r"), direction_consistency = c("direction_consistency", "direction_consistency_all", "direction_consistency_rate")))
  if (all(c("excluded_country", "logFC_correlation") %in% names(loco_summary))) {
    p_loco_summary <- loco_summary %>%
      dplyr::mutate(excluded_country = forcats::fct_reorder(as.factor(excluded_country), safe_numeric_vec(logFC_correlation)), logFC_correlation = safe_numeric_vec(logFC_correlation)) %>%
      ggplot(aes(x = excluded_country, y = logFC_correlation)) +
      geom_col(fill = "grey70", colour = "black", linewidth = 0.25) +
      coord_flip() +
      labs(title = "Leave-one-country-out stability", x = "Excluded country", y = "Correlation with main log2FC") +
      theme_panel_review(base_size = 10.5)
    
    file_base <- file.path(outdir, "result", "final_figures", "Supplementary", "Supplementary_Fig_S6_LOCO_summary")
    save_gg_pdf_png_cm(p_loco_summary, file_base, width_cm = 18, height_cm = 11, dpi = 600)
    source_file <- file.path(outdir, "result", "final_source_data", "Supplementary", "Supplementary_Fig_S6_LOCO_summary_source_data.csv")
    readr::write_csv(loco_summary, source_file)
    figure_manifest <- add_manifest(figure_manifest, "Supplementary_Fig_S6_LOCO_summary", file_base, source_file, "Summary metrics for leave-one-country-out sensitivity.")
  }
}

balanced_summary <- read_csv_if_exists(file.path(outdir, "result", "06_robustness", "balanced_country_resampling", "tables", "balanced_resampling_summary_metrics.csv"), show_message = FALSE)

if (!is.null(balanced_summary) && nrow(balanced_summary) > 0) {
  balanced_summary <- rename_with_candidates(
    balanced_summary,
    list(iteration = c("iteration", "iter"), logFC_correlation = c("logFC_correlation"), direction_consistency_all = c("direction_consistency_all", "direction_consistency"), prop_main_sig_preserved = c("prop_main_sig_preserved"))
  )
  
  plot_cols <- intersect(c("logFC_correlation", "direction_consistency_all", "prop_main_sig_preserved"), names(balanced_summary))
  if (length(plot_cols) > 0) {
    balanced_long <- balanced_summary %>%
      dplyr::select(dplyr::all_of(plot_cols)) %>%
      tidyr::pivot_longer(everything(), names_to = "metric", values_to = "value") %>%
      dplyr::mutate(metric = dplyr::recode(metric,
                                           logFC_correlation = "Effect-size correlation",
                                           direction_consistency_all = "Direction consistency",
                                           prop_main_sig_preserved = "Main DEP preserved"))
    
    p_balanced <- ggplot(balanced_long, aes(x = metric, y = value)) +
      geom_boxplot(outlier.shape = NA, fill = "grey85", color = "black", linewidth = 0.35) +
      geom_jitter(width = 0.12, size = 0.5, alpha = 0.25) +
      coord_cartesian(ylim = c(0, 1)) +
      labs(title = "Balanced country resampling stability", x = NULL, y = "Metric value") +
      theme_panel_review(base_size = 10.5) +
      theme(axis.text.x = element_text(angle = 20, hjust = 1))
    
    file_base <- file.path(outdir, "result", "final_figures", "Supplementary", "Supplementary_Fig_S7_balanced_resampling")
    save_gg_pdf_png_cm(p_balanced, file_base, width_cm = 18, height_cm = 11, dpi = 600)
    source_file <- file.path(outdir, "result", "final_source_data", "Supplementary", "Supplementary_Fig_S7_balanced_resampling_source_data.csv")
    readr::write_csv(balanced_summary, source_file)
    figure_manifest <- add_manifest(figure_manifest, "Supplementary_Fig_S7_balanced_resampling", file_base, source_file, "Distribution of robustness metrics across balanced resampling iterations.")
  }
}

if (!is.null(robustness_tbl) && nrow(robustness_tbl) > 0) {
  robustness_tbl <- rename_with_candidates(robustness_tbl, list(robustness_class = c("robustness_class"), robustness_score = c("robustness_score")))
  if ("robustness_class" %in% names(robustness_tbl)) {
    robust_counts <- robustness_tbl %>% dplyr::count(robustness_class, name = "n") %>% dplyr::arrange(dplyr::desc(n))
    p_robust <- robust_counts %>%
      dplyr::mutate(robustness_class = forcats::fct_reorder(robustness_class, n)) %>%
      ggplot(aes(x = robustness_class, y = n)) +
      geom_col(fill = "grey70", color = "black", linewidth = 0.25) +
      coord_flip() +
      labs(title = "Formal robustness classification", x = NULL, y = "Number of proteins") +
      theme_panel_review(base_size = 10.5)
    
    file_base <- file.path(outdir, "result", "final_figures", "Supplementary", "Supplementary_Fig_S8_robustness_classification")
    save_gg_pdf_png_cm(p_robust, file_base, width_cm = 18, height_cm = 10, dpi = 600)
    source_file <- file.path(outdir, "result", "final_source_data", "Supplementary", "Supplementary_Fig_S8_robustness_classification_source_data.csv")
    readr::write_csv(robustness_tbl, source_file)
    figure_manifest <- add_manifest(figure_manifest, "Supplementary_Fig_S8_robustness_classification", file_base, source_file, "Formal protein-level robustness classification from Script 01.")
  }
}

###############################################################################
# 15_export_manifest_and_session
###############################################################################

manifest_file <- file.path(outdir, "result", "final_figures", "figure_manifest.csv")
readr::write_csv(figure_manifest, manifest_file)

writeLines(capture.output(utils::sessionInfo()), con = file.path(outdir, "result", "final_figures", "logs", "script02_sessionInfo.txt"))

cat("\nSCRIPT 02 COMPLETE.\n")
cat("Figure manifest:\n", manifest_file, "\n")
cat("Main Figure 2 composite:\n", file.path(outdir, "result", "final_figures", "Figure2", "Figure2_composite_18cm_optimized.pdf"), "\n")
###############################################################################
# END SCRIPT 02
###############################################################################


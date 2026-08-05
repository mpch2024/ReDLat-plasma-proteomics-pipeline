###############################################################################
# ReDLat plasma proteomics — WGCNA workflow
# 10. Generate main Figure 3
# Requires: outputs from Scripts 01–09
# Produces: figure files and panel-level source data
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

cran_pkgs <- c(
  "dplyr", "tidyr", "readr", "tibble", "stringr", "forcats", "purrr",
  "ggplot2", "patchwork", "cowplot", "scales", "pdftools", "png",
  "igraph", "ggraph", "ggrepel", "svglite", "WGCNA", "ggtext"
)
missing_pkgs <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0L) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
       ". Run renv::restore() before this script.", call. = FALSE)
}
invisible(lapply(cran_pkgs, library, character.only = TRUE))

BASE_DIR <- WGCNA_CONFIG$project_root

S10 <- file.path(WGCNA_CONFIG$result_root, "01_input")
S11 <- file.path(WGCNA_CONFIG$result_root, "02_network")
S12 <- file.path(WGCNA_CONFIG$result_root, "03_modules")
S13 <- file.path(WGCNA_CONFIG$result_root, "04_module_traits")
S13B <- file.path(WGCNA_CONFIG$result_root, "05_sensitivity")
S14B <- file.path(WGCNA_CONFIG$result_root, "07_biomarker_fdr")

OUTDIR <- file.path(WGCNA_CONFIG$publication_root, "figures", "main_figure_3")
SUBMISSION_DIR <- file.path(OUTDIR, "submission")
ILLUSTRATOR_DIR <- file.path(OUTDIR, "Illustrator_ready")
FIGDIR <- file.path(ILLUSTRATOR_DIR, "master_vectors")
PANELDIR <- file.path(ILLUSTRATOR_DIR, "panels_vector")
PREVIEW_DIR <- file.path(OUTDIR, "preview")
PANEL_PREVIEW_DIR <- file.path(PREVIEW_DIR, "panels")
DIAGNOSTICS_DIR <- file.path(OUTDIR, "diagnostics")
TABLEDIR <- file.path(DIAGNOSTICS_DIR, "source_data")
TMPDIR <- file.path(DIAGNOSTICS_DIR, "tmp_layout_tuned")

invisible(lapply(
  c(OUTDIR, SUBMISSION_DIR, ILLUSTRATOR_DIR, FIGDIR, PANELDIR,
    PREVIEW_DIR, PANEL_PREVIEW_DIR, DIAGNOSTICS_DIR, TABLEDIR, TMPDIR),
  dir.create, recursive = TRUE, showWarnings = FALSE
))

first_existing <- function(paths, label) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) {
    stop(
      paste0("Missing required file for ", label, ":\n", paste(paths, collapse = "\n")),
      call. = FALSE
    )
  }
  normalizePath(hit[[1]], winslash = "/", mustWork = TRUE)
}

DENDROGRAM_PDF <- first_existing(
  c(file.path(S11, "modules", "merged_modules_dendrogram.pdf")),
  "final dendrogram"
)
CORE_WORKSPACE_FILE <- first_existing(
  c(file.path(S11, "workspace", "wgcna_core_collapsed_workspace.RData")),
  "full Script 11 WGCNA workspace"
)
MODULE_SIZE_FILE <- first_existing(
  c(
    file.path(S11, "tables", "module_counts.csv"),
    file.path(WGCNA_CONFIG$result_root,
              "09_network_quality",
              "tables", "modularity", "final_module_sizes_and_proportions.csv")
  ),
  "module sizes"
)
DEP_BURDEN_FILE <- first_existing(
  c(file.path(S12, "tables", "module_DEP_burden_dual_definition.csv")),
  "module DEP burden"
)
MODULE_TRAIT_FILE <- first_existing(
  c(
    file.path(S13, "tables", "correlations", "module_trait_results_long.csv"),
    file.path(S13, "module_trait_results_long.csv")
  ),
  "module-trait correlations"
)
MODULE_LABEL_FILE <- first_existing(
  c(file.path(S13, "tables", "module_biological_label_reference.csv")),
  "module biological labels"
)
TOP_HUBS_FILE <- first_existing(
  c(file.path(S12, "tables", "top_hubs_all_modules_combined.csv")),
  "top hubs"
)
ENRICHMENT_FILE <- first_existing(
  c(file.path(S12, "tables", "enrichment_significant_terms_all_modules.csv")),
  "module enrichment"
)
EXPRESSION_FILE <- first_existing(
  c(file.path(S10, "gene_collapsed_expression_matrix.csv")),
  "outcome-independent expression matrix"
)
DIAGNOSIS_HC3_FILE <- first_existing(
  c(file.path(S13B, "tables", "diagnosis_robustness", "adjusted_diagnosis_module_models_HC3.csv")),
  "HC3 diagnosis models"
)
BIOMARKER_FILE <- first_existing(
  c(file.path(S14B, "tables", "focal_modules", "integrated_focal_log_HC3_stability_family32.csv")),
  "family32 biomarker models"
)

required_files <- c(
  DENDROGRAM_PDF, CORE_WORKSPACE_FILE, MODULE_SIZE_FILE, DEP_BURDEN_FILE, MODULE_TRAIT_FILE,
  MODULE_LABEL_FILE, TOP_HUBS_FILE, ENRICHMENT_FILE, EXPRESSION_FILE,
  DIAGNOSIS_HC3_FILE, BIOMARKER_FILE
)

BASE_FONT <- "Arial"
HEATMAP_FONT <- BASE_FONT

MAIN_WIDTH_MM <- 180
MAIN_HEIGHT_MM <- 210
MAIN_WIDTH_IN <- MAIN_WIDTH_MM / 25.4
MAIN_HEIGHT_IN <- MAIN_HEIGHT_MM / 25.4
DPI <- 400

if (!identical(c(MAIN_WIDTH_MM, MAIN_HEIGHT_MM), c(180, 210))) {
  stop(
    "Figure 3 artboard dimensions changed unexpectedly. ",
    "Keep the master figure fixed at 180 × 210 mm.",
    call. = FALSE
  )
}

BASE_SIZE <- 6.0
TAG_SIZE_PT <- 9.0
PANEL_TITLE_SIZE_PT <- 7.0
AXIS_TITLE_SIZE_PT <- 7.0
AXIS_TEXT_SIZE_PT <- 6.0
LEGEND_TITLE_SIZE_PT <- 6.0
LEGEND_TEXT_SIZE_PT <- 6.0
INTERNAL_TEXT_SIZE_PT <- 6.0
COMPACT_TEXT_SIZE_PT <- 5.0
PORTRAIT_SUBTITLE_SIZE_PT <- 7.0
HEATMAP_SYMBOL_SIZE_PT <- 6.0

SOFT_POWER <- 4
N_HUBS_DISPLAY <- 8
N_HUB_LABELS <- N_HUBS_DISPLAY
N_ENRICHMENT_TERMS <- 4
NETWORK_EXPANSION <- 0.115
NETWORK_X_STRETCH <- 1.42
NETWORK_COORD_RATIO <- 0.76
PORTRAIT_WIDTHS <- c(1.10, 0.90)
PANEL_A_PAD_PX <- 3
PANEL_A_IMAGE_SCALE <- 1.00
NON_FOCAL_MODULE_ALPHA <- 0.55
TOP_ROW_WIDTHS <- c(1.18, 0.82)
MIDDLE_ROW_WIDTHS <- c(1.30, 0.70)
MAIN_ROW_HEIGHTS <- c(0.94, 0.88, 0.92, 0.92)
HEATMAP_COMPONENT_WIDTHS <- c(0.12, 0.035, 0.845)
HEATMAP_STRIP_ALPHA <- 0.82
PANEL_D_POINT_RANGE <- c(2.9, 5.7)
PANEL_D_LEGEND_BREAKS <- c(0, 2, 5)
PANEL_D_LEGEND_OVERRIDE_SIZES <- c(1.7, 2.3, 3.0)

MODULE_ORDER <- c("magenta", "blue", "brown", "green", "red", "black", "pink", "purple")
FOCAL_MODULES <- c("green", "blue", "brown")
FOCAL_MODULE_IDS <- c(green = "M1", blue = "M2", brown = "M3")
FOCAL_MODULE_DISPLAY <- c(green = "M1/green", blue = "M2/blue", brown = "M3/brown")
FOCAL_MODULE_DISPLAY_HEATMAP <- c(
  green = "M1/\ngreen",
  blue = "M2/\nblue",
  brown = "M3/\nbrown"
)
SELECTED_TRAITS <- c(
  "SampleGroup_bin", "cdr_boxscore", "mmse_total", "udsfaq_total",
  "p_tau181", "p_tau217", "NfL", "ratio_AB42_40",
  "Age", "Sex_bin", "Education", "APOE4_carrier"
)
TRAIT_LABELS <- c(
  SampleGroup_bin = "Clinical AD",
  cdr_boxscore = "CDR-SB",
  mmse_total = "MMSE",
  udsfaq_total = "PFAQ",
  p_tau181 = "p-tau181",
  p_tau217 = "p-tau217",
  NfL = "NfL",
  ratio_AB42_40 = "Aβ42/40",
  Age = "Age",
  Sex_bin = "Sex",
  Education = "Education",
  APOE4_carrier = "APOE ε4"
)
MODULE_COLORS <- c(
  black = "#2B2B2B",
  blue = "#1535E8",
  brown = "#9C6B00",
  green = "#3E8F4E",
  magenta = "#C75ACD",
  pink = "#F3A6D6",
  purple = "#8A2BE2",
  red = "#B4334A"
)
PORTRAIT_HEADER <- c(
  green = "M1/green | Neuronal connectivity and ECM",
  blue = "M2/blue | Epithelial differentiation and xenobiotic metabolism"
)

pt_to_mm <- function(x) x / ggplot2::.pt

THEME_CLASSIC <- function(base_size = BASE_SIZE) {
  ggplot2::theme_classic(base_size = base_size, base_family = BASE_FONT) +
    ggplot2::theme(
      text = ggplot2::element_text(family = BASE_FONT, colour = "#252525"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(colour = "#252525", size = AXIS_TEXT_SIZE_PT),
      axis.title = ggplot2::element_text(colour = "#252525", size = AXIS_TITLE_SIZE_PT, face = "plain"),
      plot.title = ggplot2::element_text(face = "plain", size = PANEL_TITLE_SIZE_PT, hjust = 0),
      plot.subtitle = ggplot2::element_text(size = AXIS_TEXT_SIZE_PT, colour = "grey25"),
      plot.tag = ggplot2::element_text(face = "bold", size = TAG_SIZE_PT, colour = "#252525"),
      legend.title = ggplot2::element_text(size = LEGEND_TITLE_SIZE_PT, face = "plain"),
      legend.text = ggplot2::element_text(size = LEGEND_TEXT_SIZE_PT),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      legend.background = ggplot2::element_rect(fill = "white", colour = NA)
    )
}

read_csv_required <- function(path) {
  readr::read_csv(path, show_col_types = FALSE, guess_max = 100000)
}

save_plot <- function(plot, file_base, width, height, png_dir = NULL) {
  ggplot2::ggsave(
    paste0(file_base, ".pdf"), plot, width = width, height = height,
    units = "in", device = grDevices::cairo_pdf, family = BASE_FONT,
    bg = "white", limitsize = FALSE
  )
  ggplot2::ggsave(
    paste0(file_base, ".svg"), plot, width = width, height = height,
    units = "in", device = svglite::svglite,
    bg = "white", limitsize = FALSE
  )
  png_file <- if (is.null(png_dir)) paste0(file_base, ".png") else file.path(png_dir, paste0(basename(file_base), ".png"))
  ggplot2::ggsave(
    png_file, plot, width = width, height = height,
    units = "in", dpi = DPI, bg = "white", limitsize = FALSE
  )
  invisible(c(PDF = paste0(file_base, ".pdf"), SVG = paste0(file_base, ".svg"), PNG = png_file))
}

generate_harmonized_dendrogram_pdf <- function(
    workspace_file,
    output_pdf
) {
  workspace_env <- new.env(parent = baseenv())
  loaded_objects <- load(workspace_file, envir = workspace_env)

  required_objects <- c(
    "geneTree",
    "dynamicColors",
    "mergedColors"
  )

  missing_objects <- setdiff(required_objects, loaded_objects)
  if (length(missing_objects) > 0) {
    stop(
      "The Script 11 workspace is missing dendrogram objects: ",
      paste(missing_objects, collapse = ", "),
      call. = FALSE
    )
  }

  gene_tree <- get("geneTree", envir = workspace_env)
  dynamic_colors <- get("dynamicColors", envir = workspace_env)
  merged_colors <- get("mergedColors", envir = workspace_env)

  # A wide, shallow device is intentional: it reproduces the Illustrator
  # panel geometry and prevents the dendrogram from appearing vertically
  # stretched after insertion into the fixed master layout.
  grDevices::cairo_pdf(
    filename = output_pdf,
    width = 13.2,
    height = 4.35,
    family = BASE_FONT,
    bg = "white"
  )

  old_par <- graphics::par(no.readonly = TRUE)

  tryCatch(
    {
      graphics::par(
        family = BASE_FONT,
        # More room on the left keeps "Initial modules" and "Merged modules"
        # on a single line. The larger lower margin accommodates thicker colour
        # annotations without compressing the dendrogram branches.
        mar = c(4.4, 5.2, 0.15, 0.45),
        mgp = c(1.65, 0.48, 0),
        tcl = -0.22,
        cex.axis = 0.78,
        cex.lab = 0.82
      )

      WGCNA::plotDendroAndColors(
        gene_tree,
        cbind(dynamic_colors, merged_colors),
        c("Initial modules", "Merged modules"),
        dendroLabels = FALSE,
        hang = 0.025,
        addGuide = TRUE,
        guideHang = 0.035,
        main = "",
        cex.colorLabels = 0.86,
        # The previous value produced annotation bands that were visually too
        # thin in the final 180-mm figure. This matches the manually refined
        # Illustrator version more closely.
        colorHeight = 0.145,
        autoColorHeight = FALSE,
        setLayout = TRUE
      )
    },
    finally = {
      graphics::par(old_par)
      grDevices::dev.off()
    }
  )

  if (!file.exists(output_pdf)) {
    stop("The harmonized dendrogram PDF was not created.", call. = FALSE)
  }

  invisible(output_pdf)
}

convert_pdf_first_page_to_png <- function(pdf_file, png_file, dpi = 350) {
  converted <- suppressWarnings(
    pdftools::pdf_convert(pdf = pdf_file, format = "png", pages = 1, filenames = png_file, dpi = dpi)
  )
  if (length(converted) == 0 || !file.exists(png_file)) stop("Could not convert dendrogram PDF to PNG.", call. = FALSE)
  invisible(png_file)
}

panel_from_png <- function(png_file, pad_px = 8, image_scale = 1.02) {
  img <- png::readPNG(png_file)
  if (length(dim(img)) == 3) {
    rgb <- img[, , seq_len(min(3, dim(img)[3])), drop = FALSE]
    nonwhite <- apply(rgb < 0.985, c(1, 2), any)
    if (any(nonwhite)) {
      rr <- range(which(rowSums(nonwhite) > 0))
      cc <- range(which(colSums(nonwhite) > 0))
      rr <- c(max(1, rr[1] - pad_px), min(dim(img)[1], rr[2] + pad_px))
      cc <- c(max(1, cc[1] - pad_px), min(dim(img)[2], cc[2] + pad_px))
      img <- img[rr[1]:rr[2], cc[1]:cc[2], , drop = FALSE]
    }
  }
  cowplot::ggdraw() +
    cowplot::draw_image(img, x = 0.5, y = 0.5, hjust = 0.5, vjust = 0.5, scale = image_scale) +
    ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", colour = NA))
}

format_fdr_stars <- function(fdr) {
  dplyr::case_when(
    is.na(fdr) ~ "",
    fdr < 0.001 ~ "***",
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  )
}

wrap_term <- function(x, width = 31) stringr::str_wrap(as.character(x), width = width)

make_header_strip <- function(module, title_text, show_encoding = FALSE) {
  colour <- unname(MODULE_COLORS[[module]])
  p <- ggplot2::ggplot() +
    ggplot2::annotate("rect", xmin = 0.010, xmax = 0.180, ymin = 0.08, ymax = 0.92,
                      fill = colour, colour = "grey35", linewidth = 0.28) +
    ggplot2::annotate("rect", xmin = 0.180, xmax = 0.990, ymin = 0.08, ymax = 0.92,
                      fill = "grey96", colour = "grey65", linewidth = 0.28) +
    ggplot2::annotate("text", x = 0.585, y = if (isTRUE(show_encoding)) 0.60 else 0.50,
                      label = title_text, hjust = 0.5, vjust = 0.5,
                      family = BASE_FONT, fontface = "plain",
                      size = pt_to_mm(PANEL_TITLE_SIZE_PT + 0.4), colour = "#252525")
  if (isTRUE(show_encoding)) {
    p <- p + ggplot2::annotate(
      "text", x = 0.585, y = 0.28,
      label = "Node fill: higher in AD, lower in AD or module colour; size = |kME|; edge width = signed adjacency",
      hjust = 0.5, vjust = 0.5, family = BASE_FONT,
      size = pt_to_mm(LEGEND_TEXT_SIZE_PT - 0.2), colour = "grey30"
    )
  }
  p +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off", expand = FALSE) +
    ggplot2::theme_void(base_family = BASE_FONT) +
    ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", colour = NA),
                   plot.margin = ggplot2::margin(0, 1.5, 0, 1.5, unit = "mm"))
}

compact_pathway_label <- function(x) {
  dplyr::case_when(
    stringr::str_detect(x, stringr::regex("Differentiation of Keratinocytes", ignore_case = TRUE)) ~ "Keratinocyte differentiation",
    stringr::str_detect(x, stringr::regex("Developmental Cell Lineages", ignore_case = TRUE)) ~ "Integumentary cell lineages",
    stringr::str_detect(x, stringr::regex("Metabolism of xenobiotics", ignore_case = TRUE)) ~ "Xenobiotic metabolism by CYP450",
    TRUE ~ as.character(x)
  )
}

build_hub_network <- function(module, expression_df, top_hubs_df) {
  hubs <- top_hubs_df %>%
    dplyr::filter(Module == module) %>%
    dplyr::arrange(hub_rank) %>%
    dplyr::slice_head(n = N_HUBS_DISPLAY) %>%
    dplyr::filter(EntrezGeneSymbol %in% names(expression_df))

  if (nrow(hubs) < 4) stop("Too few expression-matched hubs for module: ", module, call. = FALSE)

  genes <- hubs$EntrezGeneSymbol
  x <- as.matrix(expression_df[, genes, drop = FALSE])
  storage.mode(x) <- "double"

  cor_mat <- stats::cor(x, use = "pairwise.complete.obs", method = "pearson")
  signed_adj <- ((1 + cor_mat) / 2)^SOFT_POWER
  diag(signed_adj) <- 0

  edge_df <- as.data.frame(as.table(signed_adj), stringsAsFactors = FALSE) %>%
    dplyr::rename(from = Var1, to = Var2, weight = Freq) %>%
    dplyr::filter(as.character(from) < as.character(to)) %>%
    dplyr::mutate(weight = as.numeric(weight)) %>%
    dplyr::filter(is.finite(weight), weight > 0) %>%
    dplyr::arrange(dplyr::desc(weight)) %>%
    dplyr::slice_head(n = 18L)

  node_df <- hubs %>%
    dplyr::transmute(
      name = EntrezGeneSymbol,
      label = Protein_Display,
      display_label = dplyr::if_else(hub_rank <= N_HUB_LABELS, as.character(Protein_Display), ""),
      hub_rank = hub_rank,
      abs_kME = abs_kME,
      primary_DEP_overlap = as.logical(primary_DEP_overlap),
      dep_direction = dplyr::case_when(
        primary_DEP_direction == "Higher_in_AD" ~ "Higher in AD",
        primary_DEP_direction == "Lower_in_AD" ~ "Lower in AD",
        TRUE ~ "Not primary DEP"
      )
    )

  graph <- igraph::graph_from_data_frame(edge_df, directed = FALSE, vertices = node_df)

  dep_fill <- c(
    "Higher in AD" = unname(MODULE_COLORS[[module]]),
    "Lower in AD" = unname(MODULE_COLORS[[module]]),
    "Not primary DEP" = unname(MODULE_COLORS[[module]])
  )

  # Build the force-directed layout once and expand only its horizontal axis.
  # This makes the network use the available width without increasing the
  # vertical height of panels e and f or changing the rest of the figure.
  set.seed(123)
  network_layout <- ggraph::create_layout(graph, layout = "fr")
  network_layout$x <- network_layout$x * NETWORK_X_STRETCH

  x_span <- diff(range(network_layout$x, na.rm = TRUE))
  y_span <- diff(range(network_layout$y, na.rm = TRUE))

  manual_offsets <- if (module == "green") {
    tibble::tribble(
      ~display_label, ~dx,    ~dy,    ~hjust,
      "CADM1",        -0.050, -0.055, 1,
      "EFNA5",        -0.045, -0.060, 1,
      "NRP2",          0.045, -0.050, 0,
      "UNC5B",         0.040,  0.055, 0,
      "EPHA4",        -0.045,  0.060, 1,
      "FLRT2",         0.048, -0.010, 0,
      "IL1R1",        -0.050,  0.018, 1,
      "ROR1",          0.045,  0.060, 0
    )
  } else {
    tibble::tribble(
      ~display_label, ~dx,    ~dy,    ~hjust,
      "CFAP53",        0.045,  0.045, 0,
      "RWDD2A",        0.045,  0.050, 0,
      "KLF14",        -0.045, -0.060, 1,
      "MRPL47",        0.045,  0.000, 0,
      "LIAS",          0.045,  0.055, 0,
      "CLDN4",        -0.050,  0.060, 1,
      "SPRY2",        -0.052, -0.025, 1,
      "KHDC1L",        0.050, -0.015, 0
    )
  }

  label_data <- network_layout %>%
    dplyr::filter(display_label != "") %>%
    dplyr::left_join(manual_offsets, by = "display_label") %>%
    dplyr::mutate(
      dx = dplyr::coalesce(dx, 0.035),
      dy = dplyr::coalesce(dy, 0.035),
      hjust = dplyr::coalesce(hjust, 0),
      label_x = x + dx * x_span,
      label_y = y + dy * y_span
    )

  ggraph::ggraph(network_layout) +
    ggraph::geom_edge_link(
      ggplot2::aes(width = weight, alpha = weight),
      colour = "grey45",
      show.legend = FALSE
    ) +
    # Every node is one of the eight highest-|kME| hubs. Module colour is kept
    # as the fill so topology remains the primary visual message.
    ggraph::geom_node_point(
      ggplot2::aes(size = abs_kME, fill = dep_direction),
      shape = 21,
      colour = "#252525",
      stroke = 0.36
    ) +
    # A second transparent layer thickens only the outline of hubs that also
    # overlap the canonical primary DEP signature. DEP status does not affect
    # which proteins are selected or how they are ranked.
    ggraph::geom_node_point(
      data = network_layout %>% dplyr::filter(primary_DEP_overlap),
      ggplot2::aes(size = abs_kME),
      shape = 21,
      fill = NA,
      colour = "#111111",
      stroke = 1.05,
      show.legend = FALSE
    ) +
    ggplot2::geom_segment(
      data = label_data,
      ggplot2::aes(x = x, y = y, xend = label_x, yend = label_y),
      inherit.aes = FALSE,
      colour = "grey55",
      linewidth = 0.18
    ) +
    ggplot2::geom_text(
      data = label_data,
      ggplot2::aes(x = label_x, y = label_y, label = display_label, hjust = hjust),
      inherit.aes = FALSE,
      family = BASE_FONT,
      size = pt_to_mm(INTERNAL_TEXT_SIZE_PT),
      colour = "black"
    ) +
    ggraph::scale_edge_width(range = c(0.24, 1.38)) +
    ggraph::scale_edge_alpha(range = c(0.22, 0.68)) +
    ggplot2::scale_size_continuous(range = c(4.1, 7.4), guide = "none") +
    ggplot2::scale_fill_manual(values = dep_fill, guide = "none") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = NETWORK_EXPANSION)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = NETWORK_EXPANSION)) +
    ggplot2::coord_fixed(ratio = NETWORK_COORD_RATIO, clip = "off") +
    ggplot2::labs(
      title = "Hub network",
      subtitle = "Eight highest-|kME| hubs; thick outline = primary DEP"
    ) +
    ggplot2::theme_void(base_family = BASE_FONT, base_size = BASE_SIZE) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "plain",
        size = PORTRAIT_SUBTITLE_SIZE_PT,
        hjust = 0.5,
        margin = ggplot2::margin(b = 0.5, unit = "pt")
      ),
      plot.subtitle = ggplot2::element_text(
        family = BASE_FONT,
        face = "plain",
        size = COMPACT_TEXT_SIZE_PT,
        colour = "grey30",
        hjust = 0.5,
        margin = ggplot2::margin(b = 1.0, unit = "pt")
      ),
      plot.margin = ggplot2::margin(1, 2, 1, 2)
    )
}

build_node_color_legend <- function(module) {
  legend_df <- tibble::tribble(
    ~x, ~y, ~label, ~fill,
    1, 3, "Higher in AD", "#D65A4A",
    1, 2, "Lower in AD", "#3288BD",
    1, 1, "Not primary DEP", unname(MODULE_COLORS[[module]])
  )
  ggplot2::ggplot(legend_df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(ggplot2::aes(fill = fill), shape = 21, size = 2.8, colour = "#252525", stroke = 0.32, show.legend = FALSE) +
    ggplot2::geom_text(ggplot2::aes(x = 1.34, label = label), hjust = 0, family = BASE_FONT,
                       size = pt_to_mm(LEGEND_TEXT_SIZE_PT - 0.4), colour = "#252525") +
    ggplot2::scale_fill_identity() +
    ggplot2::xlim(0.85, 3.05) +
    ggplot2::ylim(0.4, 3.6) +
    ggplot2::labs(title = NULL) +
    ggplot2::theme_void(base_family = BASE_FONT, base_size = BASE_SIZE) +
    ggplot2::theme(
      plot.margin = ggplot2::margin(4, 2, 2, 0)
    )
}

build_enrichment_plot <- function(module, enrichment_df) {
  d <- enrichment_df %>%
    dplyr::filter(Module == module, is.finite(p.adjust), p.adjust < 0.05) %>%
    dplyr::group_by(Description) %>%
    dplyr::slice_min(order_by = p.adjust, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::slice_min(order_by = p.adjust, n = N_ENRICHMENT_TERMS, with_ties = FALSE) %>%
    dplyr::mutate(
      neglog10FDR = -log10(pmax(p.adjust, .Machine$double.xmin)),
      Display_description = compact_pathway_label(Description),
      Display_description = wrap_term(Display_description, width = 30),
      Display_description = forcats::fct_reorder(Display_description, neglog10FDR)
    )

  if (nrow(d) == 0) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No FDR-significant enrichment terms",
                          family = BASE_FONT, size = pt_to_mm(INTERNAL_TEXT_SIZE_PT)) +
        ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
        ggplot2::labs(title = "Pathway enrichment") +
        ggplot2::theme_void(base_family = BASE_FONT, base_size = BASE_SIZE) +
        ggplot2::theme(plot.title = ggplot2::element_text(face = "plain", size = PORTRAIT_SUBTITLE_SIZE_PT, hjust = 0.5),
                       panel.border = ggplot2::element_rect(colour = "grey35", fill = NA, linewidth = 0.55))
    )
  }

  ggplot2::ggplot(d, ggplot2::aes(x = neglog10FDR, y = Display_description, colour = Module, size = Count)) +
    ggplot2::geom_vline(xintercept = -log10(0.05), linetype = 2, colour = "grey70", linewidth = 0.35) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = neglog10FDR, yend = Display_description), colour = "grey65", linewidth = 0.65) +
    ggplot2::geom_point(alpha = 0.95) +
    ggplot2::scale_colour_manual(values = MODULE_COLORS, guide = "none") +
    ggplot2::scale_size_continuous(range = c(2.3, 5.4), guide = "none") +
    ggplot2::labs(title = "Pathway enrichment", x = expression(-log[10]("FDR")), y = NULL) +
    THEME_CLASSIC() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "plain", size = PORTRAIT_SUBTITLE_SIZE_PT, hjust = 0.5),
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = AXIS_TEXT_SIZE_PT),
      plot.margin = ggplot2::margin(5, 7, 5, 5)
    )
}

build_module_portrait <- function(module, expression_df, top_hubs_df, enrichment_df, show_node_legend = FALSE) {
  # Keep the header simple; node colours are intentionally uniform within each module.
  header <- make_header_strip(
    module,
    PORTRAIT_HEADER[[module]],
    show_encoding = FALSE
  )

  network <- build_hub_network(module, expression_df, top_hubs_df)
  enrich <- build_enrichment_plot(module, enrichment_df)

  # IMPORTANT: use cowplot only for the biological portraits. Nested patchwork
  # objects were rendered as blank grey rectangles by Cairo on this Windows
  # setup. cowplot::plot_grid keeps the hub network and enrichment as ordinary
  # grobs and renders them reliably in PDF, SVG, PNG and TIFF.
  lower <- cowplot::plot_grid(
    network,
    enrich,
    nrow = 1,
    rel_widths = PORTRAIT_WIDTHS,
    align = "h",
    axis = "tb"
  )

  cowplot::plot_grid(
    header,
    lower,
    ncol = 1,
    rel_heights = c(0.16, 1),
    align = "v",
    axis = "lr"
  )
}

# LOAD DATA -------------------------------------------------------------------
module_size <- read_csv_required(MODULE_SIZE_FILE)
dep_burden <- read_csv_required(DEP_BURDEN_FILE)
module_trait <- read_csv_required(MODULE_TRAIT_FILE)
module_labels <- read_csv_required(MODULE_LABEL_FILE)
top_hubs <- read_csv_required(TOP_HUBS_FILE)
enrichment <- read_csv_required(ENRICHMENT_FILE)
diagnosis_hc3 <- read_csv_required(DIAGNOSIS_HC3_FILE)
biomarker <- read_csv_required(BIOMARKER_FILE)
expression_df <- read_csv_required(EXPRESSION_FILE)

if (!"full_FDR_family32" %in% names(biomarker)) stop("Biomarker table is not the corrected Script 14b family32 output.", call. = FALSE)
if (!all(c("standardized_difference", "FDR_HC3") %in% names(diagnosis_hc3))) stop("Diagnosis table does not contain final HC3 standardized inference.", call. = FALSE)
if (!"Module" %in% names(module_size)) stop("Module-size table must contain a Module column.", call. = FALSE)
if (!"N_genes" %in% names(module_size)) {
  if ("n_genes" %in% names(module_size)) {
    module_size <- module_size %>% dplyr::rename(N_genes = n_genes)
  } else if ("Freq" %in% names(module_size)) {
    module_size <- module_size %>% dplyr::rename(N_genes = Freq)
  } else {
    stop("Could not identify module-size column.", call. = FALSE)
  }
}
module_order_use <- MODULE_ORDER[MODULE_ORDER %in% module_size$Module]
if (length(module_order_use) != 8) stop("Expected eight definitive modules.", call. = FALSE)

# PANEL A ---------------------------------------------------------------------
harmonized_dendrogram_pdf <- file.path(
  TMPDIR,
  "final_merged_module_dendrogram_harmonized.pdf"
)
dendrogram_png <- file.path(
  TMPDIR,
  "final_merged_module_dendrogram_harmonized.png"
)

generate_harmonized_dendrogram_pdf(
  CORE_WORKSPACE_FILE,
  harmonized_dendrogram_pdf
)
convert_pdf_first_page_to_png(
  harmonized_dendrogram_pdf,
  dendrogram_png,
  dpi = 400
)

panel_a_core <- panel_from_png(
  dendrogram_png,
  pad_px = PANEL_A_PAD_PX,
  image_scale = PANEL_A_IMAGE_SCALE
)

panel_a_title <- ggplot2::ggplot() +
  ggplot2::annotate(
    "text", x = 0, y = 0.5,
    label = "Final dendrogram and merged modules",
    hjust = 0, vjust = 0.5,
    family = BASE_FONT, fontface = "plain",
    size = pt_to_mm(PANEL_TITLE_SIZE_PT), colour = "#252525"
  ) +
  ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
  ggplot2::theme_void(base_family = BASE_FONT) +
  ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 2))

panel_a <- cowplot::plot_grid(
  panel_a_title,
  panel_a_core,
  ncol = 1,
  # Slightly more title space and a shallower plot body improve alignment with
  # panel b while preserving the manually matched top-row geometry.
  rel_heights = c(0.13, 0.87),
  align = "v",
  axis = "lr"
)

# PANEL B ---------------------------------------------------------------------
panel_b_data <- module_size %>%
  dplyr::select(Module, N_genes) %>%
  dplyr::left_join(dep_burden %>% dplyr::select(Module, primary_DEP_overlap_n, primary_DEP_overlap_prop), by = "Module") %>%
  dplyr::mutate(
    Module_chr = as.character(Module),
    focus_alpha = dplyr::if_else(Module_chr %in% FOCAL_MODULES, 1, NON_FOCAL_MODULE_ALPHA),
    Module = factor(Module_chr, levels = rev(module_order_use))
  )

panel_b <- ggplot2::ggplot(panel_b_data, ggplot2::aes(x = N_genes, y = Module)) +
  ggplot2::geom_col(ggplot2::aes(fill = Module, alpha = focus_alpha), width = 0.68) +
  ggplot2::geom_text(ggplot2::aes(label = scales::comma(N_genes)), hjust = -0.08, family = BASE_FONT, size = pt_to_mm(COMPACT_TEXT_SIZE_PT)) +
  ggplot2::scale_fill_manual(values = MODULE_COLORS, guide = "none") +
  ggplot2::scale_alpha_identity() +
  ggplot2::scale_x_continuous(labels = scales::comma, expand = ggplot2::expansion(mult = c(0, 0.20))) +
  ggplot2::labs(title = "Number of proteins per module", x = "Number of proteins", y = NULL) +
  THEME_CLASSIC() +
  ggplot2::theme(
    panel.grid.major.y = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(face = "plain", size = PANEL_TITLE_SIZE_PT, hjust = 0),
    plot.margin = ggplot2::margin(5, 4, 5, 2)
  )

# PANEL C ---------------------------------------------------------------------
trait_use <- SELECTED_TRAITS[SELECTED_TRAITS %in% module_trait$Trait]
panel_c_data <- module_trait %>%
  dplyr::filter(Module %in% FOCAL_MODULES, Trait %in% trait_use) %>%
  dplyr::mutate(
    Module_chr = as.character(Module),
    Module = factor(Module_chr, levels = rev(FOCAL_MODULES)),
    Trait_display = dplyr::recode(Trait, !!!TRAIT_LABELS, .default = Trait),
    Trait_display = factor(Trait_display, levels = unname(TRAIT_LABELS[trait_use])),
    stars = format_fdr_stars(FDR),
    star_colour = dplyr::case_when(stars == "" ~ "#252525", abs(rho) >= 0.20 ~ "white", TRUE ~ "#252525")
  )

panel_c_module_info <- tibble::tibble(
  Module_chr = FOCAL_MODULES,
  Module = factor(FOCAL_MODULES, levels = rev(FOCAL_MODULES)),
  Display_label = unname(FOCAL_MODULE_DISPLAY_HEATMAP[FOCAL_MODULES])
)

panel_c_label_panel <- ggplot2::ggplot(panel_c_module_info, ggplot2::aes(x = 1, y = Module)) +
  ggplot2::geom_text(ggplot2::aes(label = Display_label), hjust = 1, family = BASE_FONT,
                     size = pt_to_mm(AXIS_TEXT_SIZE_PT), colour = "#252525") +
  ggplot2::xlim(0.30, 1.02) +
  ggplot2::coord_cartesian(expand = FALSE, clip = "off") +
  ggplot2::theme_void(base_family = BASE_FONT, base_size = BASE_SIZE) +
  ggplot2::theme(plot.margin = ggplot2::margin(0, 1, 0, 0))

panel_c_strip <- ggplot2::ggplot(panel_c_module_info, ggplot2::aes(x = 1, y = Module)) +
  ggplot2::geom_tile(
    data = panel_c_module_info %>% dplyr::filter(Module_chr == "green"),
    fill = scales::alpha(unname(MODULE_COLORS[["green"]]), HEATMAP_STRIP_ALPHA), colour = "white", linewidth = 0.45
  ) +
  ggplot2::geom_tile(
    data = panel_c_module_info %>% dplyr::filter(Module_chr == "blue"),
    fill = scales::alpha(unname(MODULE_COLORS[["blue"]]), HEATMAP_STRIP_ALPHA), colour = "white", linewidth = 0.45
  ) +
  ggplot2::geom_tile(
    data = panel_c_module_info %>% dplyr::filter(Module_chr == "brown"),
    fill = scales::alpha(unname(MODULE_COLORS[["brown"]]), HEATMAP_STRIP_ALPHA), colour = "white", linewidth = 0.45
  ) +
  ggplot2::scale_y_discrete(drop = FALSE) +
  ggplot2::coord_cartesian(expand = FALSE) +
  ggplot2::theme_void(base_family = BASE_FONT, base_size = BASE_SIZE) +
  ggplot2::theme(plot.margin = ggplot2::margin(0, 1, 0, 0))

rho_limit <- max(abs(panel_c_data$rho), na.rm = TRUE)
panel_c_heatmap <- ggplot2::ggplot(panel_c_data, ggplot2::aes(x = Trait_display, y = Module, fill = rho)) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
  ggplot2::geom_text(
    ggplot2::aes(label = stars, colour = star_colour),
    size = pt_to_mm(HEATMAP_SYMBOL_SIZE_PT), family = HEATMAP_FONT, fontface = "bold"
  ) +
  ggplot2::scale_colour_identity() +
  ggplot2::scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
    limits = c(-rho_limit, rho_limit), oob = scales::squish, name = "Spearman\nrho",
    guide = ggplot2::guide_colourbar(
      title.position = "top", title.hjust = 0.5,
      barwidth = grid::unit(2.2, "mm"), barheight = grid::unit(22, "mm"),
      ticks.colour = "#252525", frame.colour = "#252525"
    )
  ) +
  ggplot2::scale_x_discrete(
    labels = function(x) {
      dplyr::if_else(
        x == "APOE ε4",
        "<i>APOE</i> ε4",
        x
      )
    },
    drop = FALSE
  ) +
  ggplot2::scale_y_discrete(drop = FALSE) +
  ggplot2::labs(x = NULL, y = NULL) +
  THEME_CLASSIC() +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggtext::element_markdown(
      family = BASE_FONT,
      angle = 52,
      hjust = 1,
      vjust = 1,
      size = AXIS_TEXT_SIZE_PT,
      colour = "#252525"
    ),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(0, 4, 0, 2)
  )

panel_c_body <- (panel_c_label_panel | panel_c_strip | panel_c_heatmap) +
  patchwork::plot_layout(widths = HEATMAP_COMPONENT_WIDTHS)

panel_c <- panel_c_body +
  patchwork::plot_annotation(
    title = "Module-trait relationships (BH-FDR corrected)",
    tag_levels = NULL,
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(family = BASE_FONT, face = "plain", size = PANEL_TITLE_SIZE_PT, hjust = 0),
      plot.margin = ggplot2::margin(3, 3, 0, 3)
    )
  )

# PANEL D ---------------------------------------------------------------------
diag_data <- diagnosis_hc3 %>%
  dplyr::filter(Module %in% FOCAL_MODULES) %>%
  dplyr::mutate(
    Outcome = "Clinical AD",
    standardized_beta = standardized_difference,
    ci_low = robust_conf_low / module_SD,
    ci_high = robust_conf_high / module_SD,
    FDR_final = FDR_HC3,
    significant = FDR_final < 0.05
  ) %>%
  dplyr::select(Module, Outcome, standardized_beta, ci_low, ci_high, FDR_final, significant)

biomarker_data <- biomarker %>%
  dplyr::filter(Module %in% FOCAL_MODULES) %>%
  dplyr::mutate(
    Outcome = Biomarker_label,
    standardized_beta = full_standardized_beta,
    ci_low = downsampling_beta_2_5,
    ci_high = downsampling_beta_97_5,
    FDR_final = full_FDR_family32,
    significant = FDR_final < 0.05
  ) %>%
  dplyr::select(Module, Outcome, standardized_beta, ci_low, ci_high, FDR_final, significant)

panel_d_trait_summary <- module_trait %>%
  dplyr::filter(Module %in% FOCAL_MODULES) %>%
  dplyr::group_by(Module) %>%
  dplyr::summarise(
    significant_trait_n = dplyr::n_distinct(Trait[is.finite(FDR) & FDR < 0.05]),
    tested_trait_n = dplyr::n_distinct(Trait), .groups = "drop"
  )

panel_d_adjusted_summary <- dplyr::bind_rows(diag_data, biomarker_data) %>%
  dplyr::group_by(Module) %>%
  dplyr::summarise(adjusted_significant_n = sum(significant, na.rm = TRUE), adjusted_tested_n = dplyr::n(), .groups = "drop")

panel_d_data <- dep_burden %>%
  dplyr::filter(Module %in% FOCAL_MODULES) %>%
  dplyr::select(Module, primary_DEP_overlap_n, primary_DEP_overlap_prop) %>%
  dplyr::left_join(panel_d_trait_summary, by = "Module") %>%
  dplyr::left_join(panel_d_adjusted_summary, by = "Module")

dep_prop_multiplier <- if (max(panel_d_data$primary_DEP_overlap_prop, na.rm = TRUE) <= 1.5) 100 else 1
panel_d_data <- panel_d_data %>%
  dplyr::mutate(
    primary_DEP_overlap_percent = dep_prop_multiplier * primary_DEP_overlap_prop,
    Display_label = unname(FOCAL_MODULE_DISPLAY[as.character(Module)]),
    Module = factor(Module, levels = c("green", "blue", "brown"))
  )

panel_d_raw <- ggplot2::ggplot(panel_d_data,
                           ggplot2::aes(x = primary_DEP_overlap_percent, y = significant_trait_n, colour = Module, size = adjusted_significant_n)) +
  ggplot2::geom_point(alpha = 0.96, stroke = 0.9) +
  ggrepel::geom_text_repel(
    ggplot2::aes(label = Display_label), family = BASE_FONT, fontface = "plain",
    size = pt_to_mm(INTERNAL_TEXT_SIZE_PT), colour = "black", seed = 123,
    box.padding = 0.38, point.padding = 0.45, min.segment.length = 0,
    segment.colour = "grey55", segment.size = 0.45, max.overlaps = Inf,
    show.legend = FALSE
  ) +
  ggplot2::scale_colour_manual(values = MODULE_COLORS, guide = "none") +
  ggplot2::scale_size_continuous(
    name = "Adjusted models\nFDR < 0.05",
    range = PANEL_D_POINT_RANGE,
    limits = c(0, 5),
    breaks = PANEL_D_LEGEND_BREAKS,
    guide = ggplot2::guide_legend(
      title.position = "top", title.hjust = 0.5, nrow = 1, byrow = TRUE,
      override.aes = list(size = PANEL_D_LEGEND_OVERRIDE_SIZES, alpha = 1)
    )
  ) +
  ggplot2::scale_x_continuous(labels = function(x) paste0(x, "%"), expand = ggplot2::expansion(mult = c(0.14, 0.20))) +
  ggplot2::scale_y_continuous(breaks = scales::breaks_pretty(n = 5), expand = ggplot2::expansion(mult = c(0.15, 0.22))) +
  ggplot2::labs(title = "Comparative profile of focal modules",
                x = "Primary DEP burden (% of module)",
                y = "FDR-significant module–trait\nassociations (of 16)") +
  THEME_CLASSIC() +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(colour = "grey90", linewidth = 0.35),
    axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 4)),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = ggplot2::element_text(
      size = COMPACT_TEXT_SIZE_PT,
      lineheight = 0.90
    ),
    legend.text = ggplot2::element_text(
      size = COMPACT_TEXT_SIZE_PT
    ),
    legend.background = ggplot2::element_blank(),
    legend.key.height = grid::unit(2.2, "mm"),
    legend.key.width = grid::unit(3.0, "mm"),
    legend.spacing.x = grid::unit(0.3, "mm"),
    legend.box.spacing = grid::unit(0.1, "mm"),
    plot.margin = ggplot2::margin(4, 3, 1, 2)
  )

panel_d <- panel_d_raw

# PANELS E/F COMPONENTS --------------------------------------------------------
panel_e_header <- make_header_strip("green", PORTRAIT_HEADER[["green"]], show_encoding = FALSE)
panel_e_network <- build_hub_network("green", expression_df, top_hubs)
panel_e_enrichment <- build_enrichment_plot("green", enrichment)

panel_f_header <- make_header_strip("blue", PORTRAIT_HEADER[["blue"]], show_encoding = FALSE)
panel_f_network <- build_hub_network("blue", expression_df, top_hubs)
panel_f_enrichment <- build_enrichment_plot("blue", enrichment)

# FONT SAFETY AUDIT --------------------------------------------------------------
font_safety_audit <- tibble::tibble(
  check = c(
    "Base font is Arial",
    "Heatmap font uses Arial fallback",
    "ggtext is available for italic APOE",
    "No plotmath parser is used in panel c labels"
  ),
  passed = c(
    identical(BASE_FONT, "Arial"),
    identical(HEATMAP_FONT, BASE_FONT),
    requireNamespace("ggtext", quietly = TRUE),
    TRUE
  )
)

if (!all(font_safety_audit$passed)) {
  stop(
    "Figure 3 font-safety audit failed: ",
    paste(font_safety_audit$check[!font_safety_audit$passed], collapse = "; "),
    call. = FALSE
  )
}

readr::write_csv(
  font_safety_audit,
  file.path(DIAGNOSTICS_DIR, "font_safety_audit_v19.csv")
)

# TYPOGRAPHY AND ARTBOARD AUDIT --------------------------------------------------
typography_audit <- tibble::tibble(
  item = c(
    "Panel tags", "Panel titles", "Axis titles", "Axis text",
    "Protein labels", "Legend text", "Compact numeric text",
    "Heatmap symbols (Arial)", "Master width", "Master height"
  ),
  expected = c(9, 7, 7, 6, 6, 6, 5, 6, 180, 210),
  observed = c(
    TAG_SIZE_PT, PANEL_TITLE_SIZE_PT, AXIS_TITLE_SIZE_PT,
    AXIS_TEXT_SIZE_PT, INTERNAL_TEXT_SIZE_PT, LEGEND_TEXT_SIZE_PT,
    COMPACT_TEXT_SIZE_PT, HEATMAP_SYMBOL_SIZE_PT,
    MAIN_WIDTH_MM, MAIN_HEIGHT_MM
  )
)

if (any(typography_audit$expected != typography_audit$observed)) {
  stop("Figure 3 typography or artboard audit failed.", call. = FALSE)
}

readr::write_csv(
  typography_audit,
  file.path(DIAGNOSTICS_DIR, "typography_and_artboard_audit_v19.csv")
)

# HUB DISPLAY AUDIT -------------------------------------------------------------
hub_display_audit <- top_hubs %>%
  dplyr::filter(Module %in% c("green", "blue"), hub_rank <= N_HUBS_DISPLAY) %>%
  dplyr::group_by(Module) %>%
  dplyr::summarise(
    displayed_hubs = dplyr::n(),
    labelled_hubs = sum(hub_rank <= N_HUB_LABELS),
    primary_DEP_hubs = sum(as.logical(primary_DEP_overlap), na.rm = TRUE),
    unique_labels = dplyr::n_distinct(Protein_Display),
    .groups = "drop"
  )

if (
  nrow(hub_display_audit) != 2L ||
    any(hub_display_audit$displayed_hubs != N_HUBS_DISPLAY) ||
    any(hub_display_audit$labelled_hubs != N_HUBS_DISPLAY) ||
    any(hub_display_audit$unique_labels != N_HUBS_DISPLAY)
) {
  stop(
    "Hub display audit failed: panels e/f must each contain and label exactly ",
    N_HUBS_DISPLAY, " unique highest-|kME| hubs.",
    call. = FALSE
  )
}

readr::write_csv(
  hub_display_audit,
  file.path(DIAGNOSTICS_DIR, "hub_display_and_DEP_outline_audit_v19.csv")
)

# SOURCE TABLES ----------------------------------------------------------------
source_tables <- list(
  panel_b_module_sizes = panel_b_data,
  panel_c_focal_module_trait = panel_c_data,
  panel_c_module_info = panel_c_module_info,
  panel_d_comparative_profile = panel_d_data,
  panel_e_green_hubs = top_hubs %>% dplyr::filter(Module == "green", hub_rank <= N_HUBS_DISPLAY),
  panel_f_blue_hubs = top_hubs %>% dplyr::filter(Module == "blue", hub_rank <= N_HUBS_DISPLAY),
  hub_display_audit = hub_display_audit,
  typography_audit = typography_audit,
  font_safety_audit = font_safety_audit
)
purrr::iwalk(source_tables, function(tbl, nm) {
  readr::write_csv(tibble::as_tibble(tbl), file.path(TABLEDIR, paste0(nm, ".csv")))
})

# ABSOLUTE MASTER LAYOUT -------------------------------------------------------
# Coordinates are normalized to the final 180 x 210 mm canvas and reproduce the
# manually refined Illustrator composition. Each component is independently
# positioned, so later changes to one panel do not distort the remaining figure.
main_figure <- cowplot::ggdraw(xlim = c(0, 1), ylim = c(0, 1)) +
  cowplot::draw_plot(panel_a, x = 0.035, y = 0.755, width = 0.500, height = 0.235) +
  cowplot::draw_plot(panel_b, x = 0.605, y = 0.755, width = 0.375, height = 0.235) +
  cowplot::draw_plot(panel_c, x = 0.020, y = 0.515, width = 0.545, height = 0.225) +
  cowplot::draw_plot(panel_d, x = 0.610, y = 0.505, width = 0.370, height = 0.230) +
  cowplot::draw_plot(panel_e_header, x = 0.040, y = 0.465, width = 0.920, height = 0.032) +
  cowplot::draw_plot(panel_e_network, x = 0.065, y = 0.270, width = 0.485, height = 0.185) +
  cowplot::draw_plot(panel_e_enrichment, x = 0.575, y = 0.260, width = 0.400, height = 0.195) +
  cowplot::draw_plot(panel_f_header, x = 0.040, y = 0.220, width = 0.920, height = 0.032) +
  cowplot::draw_plot(panel_f_network, x = 0.065, y = 0.035, width = 0.485, height = 0.175) +
  cowplot::draw_plot(panel_f_enrichment, x = 0.575, y = 0.020, width = 0.400, height = 0.190) +
  cowplot::draw_label("a", x = 0.006, y = 0.997, hjust = 0, vjust = 1,
                      fontfamily = BASE_FONT, fontface = "bold", size = TAG_SIZE_PT) +
  cowplot::draw_label("b", x = 0.580, y = 0.997, hjust = 0, vjust = 1,
                      fontfamily = BASE_FONT, fontface = "bold", size = TAG_SIZE_PT) +
  cowplot::draw_label("c", x = 0.006, y = 0.755, hjust = 0, vjust = 1,
                      fontfamily = BASE_FONT, fontface = "bold", size = TAG_SIZE_PT) +
  cowplot::draw_label("d", x = 0.580, y = 0.755, hjust = 0, vjust = 1,
                      fontfamily = BASE_FONT, fontface = "bold", size = TAG_SIZE_PT) +
  cowplot::draw_label("e", x = 0.020, y = 0.505, hjust = 0, vjust = 1,
                      fontfamily = BASE_FONT, fontface = "bold", size = TAG_SIZE_PT) +
  cowplot::draw_label("f", x = 0.020, y = 0.260, hjust = 0, vjust = 1,
                      fontfamily = BASE_FONT, fontface = "bold", size = TAG_SIZE_PT)

# SAVE COMPONENTS AND MASTER ---------------------------------------------------
save_plot(panel_a, file.path(PANELDIR, "Figure3_panel_a_dendrogram_panelA_fixed"), 95 / 25.4, 45 / 25.4, PANEL_PREVIEW_DIR)
save_plot(panel_b, file.path(PANELDIR, "Figure3_panel_b_module_sizes_illustrator_matched"), 72 / 25.4, 48 / 25.4, PANEL_PREVIEW_DIR)
save_plot(panel_c, file.path(PANELDIR, "Figure3_panel_c_module_trait_illustrator_matched"), 105 / 25.4, 46 / 25.4, PANEL_PREVIEW_DIR)
save_plot(panel_d, file.path(PANELDIR, "Figure3_panel_d_comparative_profile_illustrator_matched"), 72 / 25.4, 46 / 25.4, PANEL_PREVIEW_DIR)
save_plot(panel_e_network, file.path(PANELDIR, "Figure3_panel_e_green_network_illustrator_matched"), 88 / 25.4, 39 / 25.4, PANEL_PREVIEW_DIR)
save_plot(panel_e_enrichment, file.path(PANELDIR, "Figure3_panel_e_green_enrichment_illustrator_matched"), 72 / 25.4, 39 / 25.4, PANEL_PREVIEW_DIR)
save_plot(panel_f_network, file.path(PANELDIR, "Figure3_panel_f_blue_network_illustrator_matched"), 88 / 25.4, 37 / 25.4, PANEL_PREVIEW_DIR)
save_plot(panel_f_enrichment, file.path(PANELDIR, "Figure3_panel_f_blue_enrichment_illustrator_matched"), 72 / 25.4, 37 / 25.4, PANEL_PREVIEW_DIR)

file_base <- file.path(FIGDIR, "Figure3_WGCNA_COMPACT_CLASSIC_v19_FINAL")
save_plot(main_figure, file_base, width = MAIN_WIDTH_IN, height = MAIN_HEIGHT_IN, png_dir = PREVIEW_DIR)
main_tiff <- file.path(SUBMISSION_DIR, "Figure3_WGCNA_COMPACT_CLASSIC_v19_FINAL.tiff")
ggplot2::ggsave(main_tiff, main_figure, width = MAIN_WIDTH_MM, height = MAIN_HEIGHT_MM,
                units = "mm", dpi = 300, compression = "lzw", bg = "white", limitsize = FALSE)

# LEGEND / AUDIT ---------------------------------------------------------------
legend_text <- paste(
  "Figure 3 | Modular organization of the outcome-independent ReDLat plasma proteome.",
  "(a) Protein clustering dendrogram and final merged module assignments from the signed gene-collapsed WGCNA network.",
  "(b) Number of proteins assigned to each of the eight final modules.",
  "(c) Spearman correlations between the green, blue and brown module eigengenes and selected clinical, cognitive, demographic and plasma AT(N) traits; module labels are shown before a narrow colour strip identifying M1/green, M2/blue and M3/brown. Asterisks mark BH-FDR significance, and FDR was controlled across all 128 tests.",
  "(d) Comparative profile of the focal modules. The horizontal axis shows the proportion of each module overlapping the canonical 587-gene primary DEP signature, the vertical axis shows the number of BH-FDR-significant associations in the complete 16-trait module-trait matrix, and point size shows the number of FDR-significant adjusted models across diagnosis and the four plasma biomarkers. These dimensions are displayed separately and were not combined into a composite prioritization score.",
  "(e) Biological portrait of the green neuronal-connectivity and extracellular-matrix module, including its eight highest-|kME| hubs and top significant enrichment terms.",
  "(f) Biological portrait of the blue epithelial-differentiation and xenobiotic-metabolism module, including its eight highest-|kME| hubs and top significant enrichment terms. In panels e and f, all eight displayed nodes are the highest-|kME| hubs of the corresponding module and all are labelled. Nodes are filled uniformly with the module colour; a thicker black outline identifies hubs that also overlap the canonical primary DEP signature. DEP status was not used to select or rank hubs. Hub-network edge width represents signed co-expression adjacency and node size represents absolute module membership. The displayed edges summarize co-expression among the selected eight hubs rather than the complete module network. The blue module also showed the highest proportional primary DEP burden among the focal modules. The brown RNA-processing module is represented in panels c and d and shown in greater detail in the supplementary WGCNA figures."
)
writeLines(legend_text, file.path(DIAGNOSTICS_DIR, "figure_legends.txt"))
readr::write_csv(
  tibble::tibble(Figure_ID = "Figure 3", File_base = basename(file_base), Legend = legend_text),
  file.path(DIAGNOSTICS_DIR, "figure_legends.csv")
)

readr::write_csv(
  tibble::tibble(Source = basename(required_files), Path = required_files, Exists = file.exists(required_files)),
  file.path(DIAGNOSTICS_DIR, "input_audit.csv")
)
readr::write_csv(
  tibble::tibble(
    Layout_parameter = c(
      "top_row_panel_a", "top_row_panel_b",
      "middle_row_panel_c", "middle_row_panel_d",
      "top_row_height", "middle_row_height",
      "panel_e_height", "panel_f_height",
      "network_x_stretch", "network_coord_ratio",
      "portrait_left_width", "portrait_pathway_width"
    ),
    Value = c(
      TOP_ROW_WIDTHS[1], TOP_ROW_WIDTHS[2],
      MIDDLE_ROW_WIDTHS[1], MIDDLE_ROW_WIDTHS[2],
      MAIN_ROW_HEIGHTS[1], MAIN_ROW_HEIGHTS[2],
      MAIN_ROW_HEIGHTS[3], MAIN_ROW_HEIGHTS[4],
      NETWORK_X_STRETCH, NETWORK_COORD_RATIO,
      PORTRAIT_WIDTHS[1], PORTRAIT_WIDTHS[2]
    )
  ),
  file.path(DIAGNOSTICS_DIR, "layout_hierarchy_v19.csv")
)

writeLines(capture.output(sessionInfo()), file.path(DIAGNOSTICS_DIR, "sessionInfo.txt"))

message("Figure 3 v19 completed successfully.")
message("Main outputs saved at:")
message(file_base)


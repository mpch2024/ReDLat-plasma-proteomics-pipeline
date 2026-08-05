###############################################################################
# ReDLat plasma proteomics — WGCNA workflow
# 11. Generate Extended Data figures
# Requires: outputs from Scripts 01–09
# Produces: Extended Data figures, legends and source data
# Data policy: participant-level inputs and intermediate outputs remain local.
###############################################################################

local({

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

###############################################################################
# 1) PACKAGES
###############################################################################

cran_pkgs <- c(
  "dplyr", "tidyr", "readr", "tibble", "stringr", "forcats", "purrr",
  "ggplot2", "cowplot", "scales", "svglite", "ggtext", "rsvg", "systemfonts"
)

missing_pkgs <- cran_pkgs[
  !vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0L) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
       ". Run renv::restore() before this script.", call. = FALSE)
}

invisible(lapply(cran_pkgs, library, character.only = TRUE))

###############################################################################
# 2) PATHS
###############################################################################

BASE_DIR <- WGCNA_CONFIG$project_root


S12 <- file.path(WGCNA_CONFIG$result_root,
  "03_modules"
)

S13 <- file.path(WGCNA_CONFIG$result_root,
  "04_module_traits"
)

S13B <- file.path(WGCNA_CONFIG$result_root,
  "05_sensitivity"
)

S14 <- file.path(WGCNA_CONFIG$result_root,
  "06_stability"
)

S14B <- file.path(WGCNA_CONFIG$result_root,
  "07_biomarker_fdr"
)

S15B <- file.path(WGCNA_CONFIG$result_root,
  "08_preservation"
)

S16 <- file.path(WGCNA_CONFIG$result_root,
  "09_network_quality"
)

first_existing <- function(paths, label) {
  paths <- unique(paths[!is.na(paths) & nzchar(paths)])
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) {
    stop(
      "Missing required source for ", label, ":\n",
      paste(paths, collapse = "\n"),
      call. = FALSE
    )
  }
  normalizePath(hit[[1]], winslash = "/", mustWork = TRUE)
}

PATHS <- list(
  # Extended Data Fig. 7 -------------------------------------------------------
  soft_scan = first_existing(
    file.path(S16, "tables", "soft_threshold", "soft_threshold_scan_clean.csv"),
    "soft-threshold scan"
  ),
  module_separation = first_existing(
    file.path(S16, "tables", "modularity", "module_internal_external_separation_summary.csv"),
    "module separation summary"
  ),
  eigengene = first_existing(
    file.path(S16, "tables", "eigengenes", "module_eigengene_correlation_matrix.csv"),
    "module eigengene correlations"
  ),
  modularity_global = first_existing(
    file.path(S16, "tables", "modularity", "posthoc_weighted_modularity_Q.csv"),
    "global weighted modularity"
  ),
  modularity_module = first_existing(
    file.path(S16, "tables", "modularity", "module_specific_modularity_contributions.csv"),
    "module modularity contributions"
  ),

  # Extended Data Fig. 8 -------------------------------------------------------
  module_trait = first_existing(
    file.path(S13, "tables", "correlations", "module_trait_results_long.csv"),
    "complete module-trait matrix"
  ),
  adjusted_continuous = first_existing(
    file.path(S13, "tables", "regression", "adjusted_module_models.csv"),
    "48 adjusted continuous-outcome models"
  ),
  diagnosis_hc3 = first_existing(
    file.path(S13B, "tables", "diagnosis_robustness", "adjusted_diagnosis_module_models_HC3.csv"),
    "adjusted diagnosis HC3 models"
  ),
  biomarker_full32 = first_existing(
    file.path(S14B, "tables", "full", "full_all_modules_log_HC3_family32.csv"),
    "complete 32-model biomarker family"
  ),

  # Extended Data Fig. 9 -------------------------------------------------------
  kme_assigned = first_existing(
    file.path(S12, "tables", "full_kME_assigned_module_long.csv"),
    "assigned-module kME table"
  ),
  hub_summary = first_existing(
    file.path(S12, "tables", "module_hub_summary.csv"),
    "module hub summary"
  ),
  top_hubs = first_existing(
    file.path(S12, "tables", "top_hubs_all_modules_combined.csv"),
    "top hubs"
  ),
  enrichment_all = first_existing(
    file.path(S12, "tables", "enrichment_significant_terms_all_modules.csv"),
    "module enrichment terms"
  ),
  enrichment_summary = first_existing(
    file.path(S12, "tables", "enrichment_summary_by_module.csv"),
    "module enrichment summary"
  ),

  # Extended Data Fig. 10 ------------------------------------------------------
  context_country = first_existing(
    file.path(S13, "tables", "context", "module_recruitment_context_effect_summary.csv"),
    "country context effects"
  ),
  context_nested = first_existing(
    file.path(S13B, "tables", "context", "corrected_nested_site_all_samples.csv"),
    "nested site-within-country effects"
  ),
  loco = first_existing(
    file.path(S14, "tables", "loco", "loco_country_summary_by_module_trait.csv"),
    "LOCO association stability"
  ),
  downsampling = first_existing(
    file.path(S14, "tables", "downsampling", "balanced_downsampling_summary_by_module_trait.csv"),
    "balanced-downsampling association stability"
  ),
  diagnosis_stability = first_existing(
    file.path(S14, "tables", "model_stability", "diagnosis", "integrated_diagnosis_HC3_stability_summary.csv"),
    "adjusted diagnosis stability"
  ),

  # Extended Data Fig. 11 ------------------------------------------------------
  country_preservation = first_existing(
    file.path(S15B, "tables", "country", "country_fixed_geneset_primary_statistics.csv"),
    "country structural preservation"
  ),
  site_preservation = first_existing(
    file.path(S15B, "tables", "site", "site_fixed_geneset_primary_statistics.csv"),
    "site structural preservation"
  )
)

OUTDIR <- file.path(WGCNA_CONFIG$publication_root, "figures", "extended_data")

SUBMISSION_DIR <- file.path(OUTDIR, "submission")
MASTER_VECTOR_DIR <- file.path(OUTDIR, "Illustrator_ready", "master_vectors")
PANEL_VECTOR_DIR <- file.path(OUTDIR, "Illustrator_ready", "panels_vector")
PREVIEW_DIR <- file.path(OUTDIR, "preview")
SOURCE_DIR <- file.path(OUTDIR, "source_data")
DIAGNOSTICS_DIR <- file.path(OUTDIR, "diagnostics")

invisible(lapply(
  c(
    OUTDIR, SUBMISSION_DIR, MASTER_VECTOR_DIR, PANEL_VECTOR_DIR,
    PREVIEW_DIR, SOURCE_DIR, DIAGNOSTICS_DIR
  ),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

###############################################################################
# 3) STYLE, ORDERS AND HELPERS
###############################################################################

BASE_FONT <- "Arial"
HEATMAP_FONT <- BASE_FONT
PLOTMATH_LABELS_DISABLED <- TRUE

BASE_SIZE <- 6.0
TITLE_SIZE_PT <- 7.0
AXIS_TITLE_SIZE_PT <- 7.0
AXIS_TEXT_SIZE_PT <- 6.0
LEGEND_TITLE_SIZE_PT <- 6.0
LEGEND_TEXT_SIZE_PT <- 6.0
TAG_SIZE_PT <- 9.0
INTERNAL_TEXT_SIZE_PT <- 6.0
COMPACT_TEXT_SIZE_PT <- 5.0
HEATMAP_SYMBOL_SIZE_PT <- 6.0
DPI <- 400

FIG_WIDTH_MM <- 180

if (!identical(FIG_WIDTH_MM, 180)) {
  stop("Extended Data figure width must remain fixed at 180 mm.", call. = FALSE)
}

MODULE_ORDER <- c(
  "green", "blue", "brown", "black",
  "magenta", "red", "pink", "purple"
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

MODULE_LABELS <- c(
  green = "M1/green",
  blue = "M2/blue",
  brown = "M3/brown",
  black = "M4/black",
  magenta = "M5/magenta",
  red = "M6/red",
  pink = "M7/pink",
  purple = "M8/purple"
)

FOCAL_MODULES <- c("green", "blue", "brown")

TRAIT_LABELS <- c(
  SampleGroup_bin = "Clinical AD",
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
  APOE4_carrier = "APOE ε4"
)

TRAIT_ORDER <- names(TRAIT_LABELS)

format_trait_markdown <- function(x) {
  dplyr::if_else(
    as.character(x) == "APOE ε4",
    "<i>APOE</i> ε4",
    as.character(x)
  )
}

SELECTED_STABILITY_TRAITS <- c(
  "SampleGroup_bin", "cdr_boxscore", "mmse_total", "udsfaq_total",
  "p_tau181", "p_tau217", "NfL", "ratio_AB42_40", "Age"
)

BIOMARKER_LABELS <- c(
  p_tau181 = "p-tau181",
  p_tau217 = "p-tau217",
  NfL = "NfL",
  ratio_AB42_40 = "Aβ42/40"
)

pt_to_mm <- function(x) x / ggplot2::.pt

THEME_NATURE <- function(base_size = BASE_SIZE) {
  ggplot2::theme_bw(base_size = base_size, base_family = BASE_FONT) +
    ggplot2::theme(
      text = ggplot2::element_text(family = BASE_FONT, colour = "#252525"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "grey90", linewidth = 0.24),
      panel.border = ggplot2::element_rect(colour = "grey35", linewidth = 0.34),
      axis.text = ggplot2::element_text(colour = "#252525", size = AXIS_TEXT_SIZE_PT),
      axis.title = ggplot2::element_text(colour = "#252525", size = AXIS_TITLE_SIZE_PT),
      plot.title = ggplot2::element_text(face = "plain", size = TITLE_SIZE_PT, hjust = 0),
      plot.subtitle = ggplot2::element_text(size = AXIS_TEXT_SIZE_PT, colour = "grey30"),
      legend.title = ggplot2::element_text(size = LEGEND_TITLE_SIZE_PT),
      legend.text = ggplot2::element_text(size = LEGEND_TEXT_SIZE_PT),
      strip.background = ggplot2::element_rect(fill = "grey96", colour = "grey60", linewidth = 0.3),
      strip.text = ggplot2::element_text(size = AXIS_TEXT_SIZE_PT, face = "plain"),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      legend.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin = ggplot2::margin(2.2, 2.6, 2.2, 2.6)
    )
}

read_required <- function(path) {
  readr::read_csv(path, show_col_types = FALSE, guess_max = 100000)
}

safe_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))

safe_positive_limit <- function(x, fallback = 1) {
  x <- abs(safe_numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(fallback)
  out <- max(x)
  if (!is.finite(out) || out <= 0) fallback else out
}

safe_minimum <- function(x, fallback = NA_real_) {
  x <- safe_numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(fallback)
  min(x)
}

pick_col <- function(df, candidates, required = TRUE, label = "column") {
  exact <- candidates[candidates %in% names(df)][1]
  if (length(exact) > 0 && !is.na(exact)) return(exact)

  clean <- function(x) tolower(gsub("[^a-z0-9]+", "", x))
  idx <- match(clean(candidates), clean(names(df)), nomatch = 0)
  idx <- idx[idx > 0][1]
  if (length(idx) > 0 && !is.na(idx)) return(names(df)[idx])

  if (isTRUE(required)) {
    stop(
      "Could not identify ", label, ". Available columns: ",
      paste(names(df), collapse = ", "),
      call. = FALSE
    )
  }
  NA_character_
}

fdr_stars <- function(x) {
  dplyr::case_when(
    is.na(x) ~ "",
    x < 0.001 ~ "***",
    x < 0.01 ~ "**",
    x < 0.05 ~ "*",
    TRUE ~ ""
  )
}

wrap_term <- function(x, width = 34) {
  stringr::str_wrap(as.character(x), width = width)
}

module_display_factor <- function(x, reverse = TRUE) {
  levels_use <- if (reverse) rev(MODULE_ORDER) else MODULE_ORDER
  factor(as.character(x), levels = levels_use)
}

module_label <- function(x) {
  unname(MODULE_LABELS[as.character(x)])
}

save_master_figure <- function(plot, filename, width_mm, height_mm) {
  if (!is.finite(width_mm) || !is.finite(height_mm) || width_mm <= 0 || height_mm <= 0) {
    stop("Invalid master-figure dimensions for ", filename, call. = FALSE)
  }

  pdf_file <- file.path(MASTER_VECTOR_DIR, paste0(filename, ".pdf"))
  svg_file <- file.path(MASTER_VECTOR_DIR, paste0(filename, ".svg"))
  png_file <- file.path(PREVIEW_DIR, paste0(filename, ".png"))
  tiff_file <- file.path(SUBMISSION_DIR, paste0(filename, ".tiff"))

  tryCatch(
    ggplot2::ggsave(
      svg_file, plot,
      width = width_mm, height = height_mm, units = "mm",
      device = svglite::svglite,
      bg = "white", limitsize = FALSE
    ),
    error = function(e) {
      stop(
        "Master SVG could not be rendered: ", filename,
        "\n", conditionMessage(e),
        call. = FALSE
      )
    }
  )

  if (!file.exists(svg_file) || file.info(svg_file)$size <= 0) {
    stop("Master SVG was not created or is empty: ", filename, call. = FALSE)
  }

  tryCatch(
    rsvg::rsvg_pdf(svg_file, pdf_file),
    error = function(e) {
      stop(
        "SVG-to-PDF conversion failed: ", filename,
        "\n", conditionMessage(e),
        call. = FALSE
      )
    }
  )

  ggplot2::ggsave(
    png_file, plot,
    width = width_mm, height = height_mm, units = "mm",
    dpi = DPI, bg = "white", limitsize = FALSE
  )

  ggplot2::ggsave(
    tiff_file, plot,
    width = width_mm, height = height_mm, units = "mm",
    dpi = 300, compression = "lzw",
    bg = "white", limitsize = FALSE
  )

  invisible(c(PDF = pdf_file, SVG = svg_file, PNG = png_file, TIFF = tiff_file))
}

save_panel <- function(plot, figure_id, panel_id, width_mm, height_mm) {
  base <- file.path(PANEL_VECTOR_DIR, paste0(figure_id, "_panel_", panel_id))
  svg_file <- paste0(base, ".svg")
  pdf_file <- paste0(base, ".pdf")

  ggplot2::ggsave(
    svg_file, plot,
    width = width_mm, height = height_mm, units = "mm",
    device = svglite::svglite,
    bg = "white", limitsize = FALSE
  )

  if (!file.exists(svg_file) || file.info(svg_file)$size <= 0) {
    stop("Panel SVG was not created or is empty: ", figure_id, panel_id, call. = FALSE)
  }

  rsvg::rsvg_pdf(svg_file, pdf_file)
}

write_source <- function(tbl, filename) {
  readr::write_csv(
    tibble::as_tibble(tbl),
    file.path(SOURCE_DIR, paste0(filename, ".csv"))
  )
}

make_tagged_grid <- function(..., labels, rel_widths = NULL, rel_heights = NULL,
                             nrow = NULL, ncol = NULL, align = "none", axis = "none") {
  plots <- list(...)
  args <- c(
    plots,
    list(
      labels = labels,
      label_fontfamily = BASE_FONT,
      label_fontface = "bold",
      label_size = TAG_SIZE_PT,
      label_x = 0.008,
      label_y = 1.010,
      hjust = 0,
      vjust = 1,
      align = align,
      axis = axis
    )
  )
  if (!is.null(rel_widths)) args$rel_widths <- rel_widths
  if (!is.null(rel_heights)) args$rel_heights <- rel_heights
  if (!is.null(nrow)) args$nrow <- nrow
  if (!is.null(ncol)) args$ncol <- ncol
  do.call(cowplot::plot_grid, args)
}


make_module_label_panel <- function(module_levels, label_size = AXIS_TEXT_SIZE_PT) {
  d <- tibble::tibble(
    Module = factor(module_levels, levels = module_levels),
    Label = module_label(module_levels)
  )

  ggplot2::ggplot(d, ggplot2::aes(x = 1, y = Module, label = Label)) +
    ggplot2::geom_text(
      hjust = 1,
      family = BASE_FONT,
      size = pt_to_mm(label_size),
      colour = "#252525"
    ) +
    ggplot2::xlim(0.15, 1.02) +
    ggplot2::scale_y_discrete(drop = FALSE) +
    ggplot2::coord_cartesian(expand = FALSE, clip = "off") +
    ggplot2::theme_void(base_family = BASE_FONT) +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 1, 0, 0))
}

make_module_strip_panel <- function(module_levels) {
  d <- tibble::tibble(
    Module = factor(module_levels, levels = module_levels),
    x = 1,
    Module_chr = as.character(module_levels)
  )

  ggplot2::ggplot(d, ggplot2::aes(x = x, y = Module, fill = Module_chr)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.35) +
    ggplot2::scale_fill_manual(values = MODULE_COLORS, guide = "none") +
    ggplot2::scale_y_discrete(drop = FALSE) +
    ggplot2::coord_cartesian(expand = FALSE) +
    ggplot2::theme_void(base_family = BASE_FONT) +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 1, 0, 0))
}

wrap_module_heatmap <- function(heatmap_plot, module_levels,
                                label_width = 0.13, strip_width = 0.025) {
  label_plot <- make_module_label_panel(module_levels)
  strip_plot <- make_module_strip_panel(module_levels)

  cowplot::plot_grid(
    label_plot,
    strip_plot,
    heatmap_plot,
    nrow = 1,
    rel_widths = c(label_width, strip_width, 1 - label_width - strip_width),
    align = "h",
    axis = "tb"
  )
}

trim_png_whitespace <- function(input_path, output_path, threshold = 0.992, padding = 18L) {
  img <- png::readPNG(input_path)
  if (length(dim(img)) < 3L || dim(img)[3] < 3L) {
    png::writePNG(img, output_path)
    return(output_path)
  }

  rgb <- img[, , 1:3, drop = FALSE]
  nonwhite <- apply(rgb, c(1, 2), function(v) any(v < threshold))
  if (dim(img)[3] >= 4L) nonwhite <- nonwhite & img[, , 4] > 0.01

  rows <- which(rowSums(nonwhite) > 0)
  cols <- which(colSums(nonwhite) > 0)
  if (length(rows) == 0L || length(cols) == 0L) {
    png::writePNG(img, output_path)
    return(output_path)
  }

  r1 <- max(1L, min(rows) - as.integer(padding))
  r2 <- min(dim(img)[1], max(rows) + as.integer(padding))
  c1 <- max(1L, min(cols) - as.integer(padding))
  c2 <- min(dim(img)[2], max(cols) + as.integer(padding))

  cropped <- img[r1:r2, c1:c2, , drop = FALSE]
  png::writePNG(cropped, output_path)
  output_path
}

render_pdf_panel <- function(
  pdf_path,
  png_path,
  dpi = 300,
  crop_top = 0,
  crop_bottom = 0,
  crop_left = 0,
  crop_right = 0
) {
  if (!file.exists(pdf_path)) {
    stop("PDF panel source was not found: ", pdf_path, call. = FALSE)
  }

  crop_values <- c(crop_top, crop_bottom, crop_left, crop_right)
  if (any(!is.finite(crop_values)) || any(crop_values < 0) ||
      crop_top + crop_bottom >= 0.90 || crop_left + crop_right >= 0.90) {
    stop("Invalid fractional crop supplied to render_pdf_panel().", call. = FALSE)
  }

  needs_render <- !file.exists(png_path) ||
    file.info(png_path)$mtime < file.info(pdf_path)$mtime

  if (isTRUE(needs_render)) {
    pdftools::pdf_convert(
      pdf = pdf_path,
      format = "png",
      pages = 1,
      filenames = png_path,
      dpi = dpi,
      verbose = FALSE
    )
  }

  trimmed_path <- sub("\\.png$", "_trimmed.png", png_path, ignore.case = TRUE)
  needs_trim <- !file.exists(trimmed_path) ||
    file.info(trimmed_path)$mtime < file.info(png_path)$mtime
  if (isTRUE(needs_trim)) {
    trim_png_whitespace(png_path, trimmed_path, threshold = 0.992, padding = 18L)
  }

  image_array <- png::readPNG(trimmed_path)
  nr <- dim(image_array)[1]
  nc <- dim(image_array)[2]

  r1 <- max(1L, floor(1 + nr * crop_top))
  r2 <- min(nr, ceiling(nr * (1 - crop_bottom)))
  c1 <- max(1L, floor(1 + nc * crop_left))
  c2 <- min(nc, ceiling(nc * (1 - crop_right)))

  if (r2 <= r1 || c2 <= c1) {
    stop("Fractional crop removed the entire rendered PDF panel.", call. = FALSE)
  }

  image_array <- image_array[r1:r2, c1:c2, , drop = FALSE]
  image_grob <- grid::rasterGrob(image_array, interpolate = TRUE)

  cowplot::ggdraw() +
    cowplot::draw_grob(image_grob, x = 0, y = 0, width = 1, height = 1)
}

add_absolute_panel <- function(canvas, plot, tag, x, y, width, height,
                               tag_x = x - 0.027, tag_y = y + height + 0.012) {
  canvas +
    cowplot::draw_plot(plot, x = x, y = y, width = width, height = height) +
    cowplot::draw_label(
      tag,
      x = tag_x,
      y = tag_y,
      hjust = 0,
      vjust = 1,
      fontfamily = BASE_FONT,
      fontface = "bold",
      size = TAG_SIZE_PT,
      colour = "#111111"
    )
}

add_absolute_grob <- function(canvas, grob, x, y, width, height) {
  canvas + cowplot::draw_grob(grob, x = x, y = y, width = width, height = height)
}

extract_horizontal_legend <- function(plot, barwidth = 34, barheight = 4.0, nrow = 1) {
  cowplot::get_legend(
    plot +
      ggplot2::guides(
        fill = ggplot2::guide_colorbar(
          direction = "horizontal",
          title.position = "top",
          title.hjust = 0.5,
          barwidth = grid::unit(barwidth, "mm"),
          barheight = grid::unit(barheight, "mm")
        ),
        colour = ggplot2::guide_legend(nrow = nrow, byrow = TRUE),
        shape = ggplot2::guide_legend(nrow = nrow, byrow = TRUE),
        size = ggplot2::guide_legend(nrow = nrow, byrow = TRUE)
      ) +
      ggplot2::theme(
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.margin = ggplot2::margin(0, 0, 0, 0),
        legend.box.margin = ggplot2::margin(0, 0, 0, 0)
      )
  )
}

module_scale_x <- function(...) {
  ggplot2::scale_x_discrete(labels = MODULE_LABELS, ...)
}

module_scale_y <- function(...) {
  ggplot2::scale_y_discrete(labels = MODULE_LABELS, ...)
}

###############################################################################
# 4) LOAD DEFINITIVE TABLES
###############################################################################

soft_scan <- read_required(PATHS$soft_scan)
module_separation <- read_required(PATHS$module_separation)
eigengene_wide <- read_required(PATHS$eigengene)
modularity_global <- read_required(PATHS$modularity_global)
modularity_module <- read_required(PATHS$modularity_module)

module_trait <- read_required(PATHS$module_trait)
adjusted_continuous <- read_required(PATHS$adjusted_continuous)
diagnosis_hc3 <- read_required(PATHS$diagnosis_hc3)
biomarker_full32 <- read_required(PATHS$biomarker_full32)

kme_assigned <- read_required(PATHS$kme_assigned)
hub_summary <- read_required(PATHS$hub_summary)
top_hubs <- read_required(PATHS$top_hubs)
enrichment_all <- read_required(PATHS$enrichment_all)
enrichment_summary <- read_required(PATHS$enrichment_summary)

context_country <- read_required(PATHS$context_country)
context_nested <- read_required(PATHS$context_nested)
loco <- read_required(PATHS$loco)
downsampling <- read_required(PATHS$downsampling)
diagnosis_stability <- read_required(PATHS$diagnosis_stability)

country_preservation <- read_required(PATHS$country_preservation)
site_preservation <- read_required(PATHS$site_preservation)

if (!"FDR_family32" %in% names(biomarker_full32)) {
  stop("The full biomarker source is not the corrected 32-model family.", call. = FALSE)
}

if (nrow(biomarker_full32) != 32L) {
  stop(
    "Expected exactly 32 full-sample module–biomarker models; observed ",
    nrow(biomarker_full32), ".",
    call. = FALSE
  )
}

###############################################################################
# 5) EXTENDED DATA FIGURE 7
# Network construction, quality and global organization
#
# FINAL COMPACT-CLASSIC ARCHITECTURE
# The six-panel 3 × 2 structure intentionally follows the original balanced
# composition preferred during visual review. The participant dendrogram remains
# available as the original Script 11 QC output but is not repeated in this master
# figure because its 639 terminal labels reduce legibility without adding inference.
###############################################################################

SELECTED_POWER <- 4
selected_soft <- soft_scan %>%
  dplyr::filter(Power == SELECTED_POWER) %>%
  dplyr::slice(1)

if (nrow(selected_soft) != 1L) {
  stop("Soft-threshold scan does not contain power 4.", call. = FALSE)
}

# a, signed scale-free topology fit ---------------------------------------------

p7a <- ggplot2::ggplot(
  soft_scan,
  ggplot2::aes(x = Power, y = Signed_scale_free_R2)
) +
  ggplot2::geom_hline(
    yintercept = 0.90,
    linetype = 2,
    colour = "grey55",
    linewidth = 0.35
  ) +
  ggplot2::geom_line(linewidth = 0.62, colour = "grey35") +
  ggplot2::geom_point(size = 1.65, colour = "grey35") +
  ggplot2::geom_point(
    data = selected_soft,
    size = 3.0,
    shape = 21,
    fill = "white",
    colour = "black",
    stroke = 0.60
  ) +
  ggplot2::annotate(
    "label",
    x = 6.5,
    y = 0.82,
    label = paste0(
      "β = 4; R² = ",
      sprintf("%.3f", selected_soft$Signed_scale_free_R2)
    ),
    family = BASE_FONT,
    size = pt_to_mm(COMPACT_TEXT_SIZE_PT),
    hjust = 0,
    label.size = 0,
    fill = scales::alpha("white", 0.86)
  ) +
  ggplot2::labs(
    title = "Signed scale-free topology fit",
    x = "Soft-threshold power",
    y = "Signed R²"
  ) +
  THEME_NATURE()

# b, mean connectivity ----------------------------------------------------------

p7b <- ggplot2::ggplot(
  soft_scan,
  ggplot2::aes(x = Power, y = Mean_connectivity)
) +
  ggplot2::geom_line(linewidth = 0.62, colour = "grey35") +
  ggplot2::geom_point(size = 1.65, colour = "grey35") +
  ggplot2::geom_point(
    data = selected_soft,
    size = 3.0,
    shape = 21,
    fill = "white",
    colour = "black",
    stroke = 0.60
  ) +
  ggplot2::annotate(
    "label",
    x = 6.0,
    y = selected_soft$Mean_connectivity * 1.08,
    label = paste0(
      "β = 4; connectivity = ",
      scales::comma(round(selected_soft$Mean_connectivity, 1))
    ),
    family = BASE_FONT,
    size = pt_to_mm(COMPACT_TEXT_SIZE_PT),
    hjust = 0,
    label.size = 0,
    fill = scales::alpha("white", 0.86)
  ) +
  ggplot2::labs(
    title = "Mean network connectivity",
    x = "Soft-threshold power",
    y = "Mean connectivity"
  ) +
  THEME_NATURE()

# c, informative module-membership distributions -------------------------------

kme_module_col7 <- pick_col(
  kme_assigned,
  c("Module", "module"),
  TRUE,
  "module in assigned-kME table"
)
kme_abs_col7 <- pick_col(
  kme_assigned,
  c("abs_kME", "Abs_kME", "Abs_assigned_kME"),
  TRUE,
  "absolute assigned kME"
)

ed7_kme <- kme_assigned %>%
  dplyr::transmute(
    Module = as.character(.data[[kme_module_col7]]),
    abs_kME = safe_numeric(.data[[kme_abs_col7]])
  ) %>%
  dplyr::filter(Module %in% MODULE_ORDER, is.finite(abs_kME)) %>%
  dplyr::mutate(Module = factor(Module, levels = rev(MODULE_ORDER)))

ed7_kme_median <- ed7_kme %>%
  dplyr::group_by(Module) %>%
  dplyr::summarise(
    abs_kME = stats::median(abs_kME, na.rm = TRUE),
    .groups = "drop"
  )

p7c <- ggplot2::ggplot(
  ed7_kme,
  ggplot2::aes(x = abs_kME, y = Module, fill = Module)
) +
  ggplot2::geom_violin(
    trim = FALSE,
    scale = "width",
    orientation = "y",
    alpha = 0.52,
    linewidth = 0.30
  ) +
  ggplot2::geom_boxplot(
    width = 0.14,
    orientation = "y",
    outlier.shape = NA,
    fill = "white",
    colour = "grey25",
    linewidth = 0.32
  ) +
  ggplot2::geom_point(
    data = ed7_kme_median,
    ggplot2::aes(x = abs_kME, y = Module),
    inherit.aes = FALSE,
    shape = 21,
    fill = "white",
    colour = "black",
    size = 1.65,
    stroke = 0.42
  ) +
  ggplot2::scale_fill_manual(values = MODULE_COLORS, guide = "none") +
  module_scale_y(drop = FALSE) +
  ggplot2::labs(
    title = "Assigned module membership",
    subtitle = "Violin density, boxplot and median",
    x = "Absolute module membership (|kME|)",
    y = NULL
  ) +
  THEME_NATURE() +
  ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())

# d, adjacency and TOM separation on a shared log scale -------------------------

ed7_separation <- module_separation %>%
  dplyr::filter(Module %in% MODULE_ORDER) %>%
  dplyr::transmute(
    Module = as.character(Module),
    `Adjacency ratio` = safe_numeric(Adjacency_within_external_ratio),
    `TOM ratio` = safe_numeric(TOM_within_external_ratio)
  ) %>%
  tidyr::pivot_longer(
    cols = c(`Adjacency ratio`, `TOM ratio`),
    names_to = "Metric",
    values_to = "Ratio"
  ) %>%
  dplyr::filter(is.finite(Ratio), Ratio > 0) %>%
  dplyr::mutate(
    Module = factor(Module, levels = rev(MODULE_ORDER)),
    Metric = factor(Metric, levels = c("Adjacency ratio", "TOM ratio"))
  )

p7d <- ggplot2::ggplot(
  ed7_separation,
  ggplot2::aes(x = Ratio, y = Module, colour = Metric, shape = Metric)
) +
  ggplot2::geom_vline(
    xintercept = 1,
    linetype = 2,
    colour = "grey55",
    linewidth = 0.35
  ) +
  ggplot2::geom_line(
    ggplot2::aes(group = Module),
    colour = "grey78",
    linewidth = 0.55
  ) +
  ggplot2::geom_point(size = 2.45) +
  ggplot2::scale_x_log10(
    breaks = scales::log_breaks(n = 5),
    labels = scales::label_number(accuracy = 0.01)
  ) +
  ggplot2::scale_colour_manual(
    values = c(
      "Adjacency ratio" = "#355C7D",
      "TOM ratio" = "#9A6FB0"
    )
  ) +
  ggplot2::scale_shape_manual(
    values = c(
      "Adjacency ratio" = 16,
      "TOM ratio" = 17
    )
  ) +
  module_scale_y(drop = FALSE) +
  ggplot2::labs(
    title = "Within/external network separation",
    subtitle = "Shared log scale; dashed line marks ratio = 1",
    x = "Within/external mean-weight ratio (log scale)",
    y = NULL,
    colour = NULL,
    shape = NULL
  ) +
  THEME_NATURE() +
  ggplot2::theme(
    panel.grid.major.y = ggplot2::element_blank(),
    legend.position = "bottom",
    legend.box = "horizontal"
  ) +
  ggplot2::guides(
    colour = ggplot2::guide_legend(nrow = 1, byrow = TRUE),
    shape = ggplot2::guide_legend(nrow = 1, byrow = TRUE)
  )

# e, module eigengene correlations ---------------------------------------------

ed7_eig <- eigengene_wide %>%
  tidyr::pivot_longer(
    cols = -Module,
    names_to = "Module_2",
    values_to = "Correlation"
  ) %>%
  dplyr::filter(Module %in% MODULE_ORDER, Module_2 %in% MODULE_ORDER) %>%
  dplyr::mutate(
    Module = factor(Module, levels = rev(MODULE_ORDER)),
    Module_2 = factor(Module_2, levels = MODULE_ORDER)
  )

p7e <- ggplot2::ggplot(
  ed7_eig,
  ggplot2::aes(x = Module_2, y = Module, fill = Correlation)
) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.35) +
  ggplot2::geom_text(
    ggplot2::aes(label = sprintf("%.2f", Correlation)),
    family = BASE_FONT,
    size = pt_to_mm(COMPACT_TEXT_SIZE_PT)
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Eigengene
correlation"
  ) +
  module_scale_x(drop = FALSE) +
  module_scale_y(drop = FALSE) +
  ggplot2::labs(
    title = "Final module eigengene correlations",
    x = NULL,
    y = NULL
  ) +
  THEME_NATURE() +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 35, hjust = 1)
  )

# f, post hoc weighted modularity -----------------------------------------------

q_row <- modularity_global %>%
  dplyr::filter(Metric == "Posthoc_weighted_modularity_Q") %>%
  dplyr::slice(1)

if (nrow(q_row) != 1L) {
  stop(
    "Global modularity table does not contain Posthoc_weighted_modularity_Q.",
    call. = FALSE
  )
}

global_q <- safe_numeric(q_row$Value)

ed7_q <- modularity_module %>%
  dplyr::filter(Module %in% MODULE_ORDER) %>%
  dplyr::mutate(
    Module = factor(Module, levels = rev(MODULE_ORDER)),
    Module_chr = as.character(Module)
  )

p7f <- ggplot2::ggplot(
  ed7_q,
  ggplot2::aes(
    x = Modularity_Q_contribution,
    y = Module,
    fill = Module_chr
  )
) +
  ggplot2::geom_vline(
    xintercept = 0,
    colour = "grey55",
    linewidth = 0.35
  ) +
  ggplot2::geom_col(width = 0.68) +
  ggplot2::scale_fill_manual(values = MODULE_COLORS, guide = "none") +
  module_scale_y(drop = FALSE) +
  ggplot2::labs(
    title = "Weighted modularity contributions",
    subtitle = paste0(
      "Post hoc descriptive metric; global Q = ",
      sprintf("%.3f", global_q)
    ),
    x = "Contribution to weighted modularity Q",
    y = NULL
  ) +
  THEME_NATURE() +
  ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())

# Original balanced 3 × 2 composition ------------------------------------------

ed7_row1 <- make_tagged_grid(
  p7a, p7b,
  labels = c("a", "b"),
  nrow = 1,
  rel_widths = c(1, 1)
)

ed7_row2 <- make_tagged_grid(
  p7c, p7d,
  labels = c("c", "d"),
  nrow = 1,
  rel_widths = c(1, 1)
)

ed7_row3 <- make_tagged_grid(
  p7e, p7f,
  labels = c("e", "f"),
  nrow = 1,
  rel_widths = c(1.07, 0.93)
)

ed7_figure <- cowplot::plot_grid(
  ed7_row1,
  ed7_row2,
  ed7_row3,
  ncol = 1,
  rel_heights = c(0.93, 1.00, 1.08)
)

save_panel(p7a, "Extended_Data_Figure_07", "a", 86, 54)
save_panel(p7b, "Extended_Data_Figure_07", "b", 86, 54)
save_panel(p7c, "Extended_Data_Figure_07", "c", 86, 58)
save_panel(p7d, "Extended_Data_Figure_07", "d", 86, 58)
save_panel(p7e, "Extended_Data_Figure_07", "e", 94, 64)
save_panel(p7f, "Extended_Data_Figure_07", "f", 78, 64)

save_master_figure(
  ed7_figure,
  "Extended_Data_Figure_07_WGCNA_network_quality",
  FIG_WIDTH_MM,
  182
)

write_source(
  soft_scan,
  "Extended_Data_Figure_07a_b_soft_threshold"
)
write_source(
  ed7_kme,
  "Extended_Data_Figure_07c_module_membership"
)
write_source(
  ed7_separation,
  "Extended_Data_Figure_07d_module_separation"
)
write_source(
  ed7_eig,
  "Extended_Data_Figure_07e_eigengene_correlations"
)
write_source(
  ed7_q,
  "Extended_Data_Figure_07f_modularity_contributions"
)
write_source(
  modularity_global,
  "Extended_Data_Figure_07_global_modularity"
)

###############################################################################
# 6) EXTENDED DATA FIGURE 8
# Complete module–trait and adjusted association landscape
###############################################################################

trait_order_use <- TRAIT_ORDER[TRAIT_ORDER %in% module_trait$Trait]
module_levels_heatmap <- rev(MODULE_ORDER)

ed8_trait <- module_trait %>%
  dplyr::filter(Module %in% MODULE_ORDER, Trait %in% trait_order_use) %>%
  dplyr::mutate(
    Module = factor(Module, levels = module_levels_heatmap),
    Trait_label = dplyr::recode(Trait, !!!TRAIT_LABELS, .default = Trait),
    Trait_label = factor(Trait_label, levels = unname(TRAIT_LABELS[trait_order_use])),
    stars = fdr_stars(FDR),
    star_colour = dplyr::if_else(abs(rho) >= 0.18, "white", "black")
  )

n_trait_sig <- sum(ed8_trait$FDR < 0.05, na.rm = TRUE)
rho_limit <- safe_positive_limit(ed8_trait$rho)

p8a_heat <- ggplot2::ggplot(
  ed8_trait,
  ggplot2::aes(x = Trait_label, y = Module, fill = rho)
) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.42) +
  ggplot2::geom_text(
    ggplot2::aes(label = stars, colour = star_colour),
    family = HEATMAP_FONT,
    fontface = "bold",
    size = pt_to_mm(HEATMAP_SYMBOL_SIZE_PT)
  ) +
  ggplot2::scale_colour_identity() +
  ggplot2::scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = c(-rho_limit, rho_limit), oob = scales::squish,
    name = "Spearman rho"
  ) +
  ggplot2::scale_x_discrete(
    labels = format_trait_markdown,
    drop = FALSE
  ) +
  ggplot2::labs(
    title = "Complete module–trait correlation matrix",
    subtitle = paste0(n_trait_sig, " of ", nrow(ed8_trait), " associations reached BH-FDR < 0.05"),
    x = NULL,
    y = NULL
  ) +
  THEME_NATURE() +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggtext::element_markdown(
      family = BASE_FONT, angle = 48, hjust = 1,
      size = AXIS_TEXT_SIZE_PT, colour = "#252525"
    ),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    legend.position = "right"
  )

p8a <- wrap_module_heatmap(p8a_heat, module_levels_heatmap, 0.14, 0.025)

ed8_cont <- adjusted_continuous %>%
  dplyr::filter(Module %in% MODULE_ORDER, model_status == "ok") %>%
  dplyr::mutate(
    Module = factor(Module, levels = module_levels_heatmap),
    Outcome_label = dplyr::recode(Outcome, !!!TRAIT_LABELS, .default = Outcome),
    stars = fdr_stars(FDR),
    signed_t = safe_numeric(statistic),
    star_colour = dplyr::if_else(abs(signed_t) >= 2.6, "white", "black")
  )

continuous_order <- unique(ed8_cont$Outcome_label)
ed8_cont$Outcome_label <- factor(ed8_cont$Outcome_label, levels = continuous_order)
n_cont_sig <- sum(ed8_cont$FDR < 0.05, na.rm = TRUE)
t_limit <- safe_positive_limit(ed8_cont$signed_t)

p8b_heat <- ggplot2::ggplot(
  ed8_cont,
  ggplot2::aes(x = Outcome_label, y = Module, fill = signed_t)
) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.42) +
  ggplot2::geom_text(
    ggplot2::aes(label = stars, colour = star_colour),
    family = HEATMAP_FONT,
    fontface = "bold",
    size = pt_to_mm(HEATMAP_SYMBOL_SIZE_PT)
  ) +
  ggplot2::scale_colour_identity() +
  ggplot2::scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = c(-t_limit, t_limit), oob = scales::squish,
    name = "Signed t statistic"
  ) +
  ggplot2::scale_x_discrete(
    labels = format_trait_markdown,
    drop = FALSE
  ) +
  ggplot2::labs(
    title = "Covariate-adjusted continuous-outcome models",
    subtitle = paste0(n_cont_sig, " of ", nrow(ed8_cont), " models reached BH-FDR < 0.05"),
    x = NULL,
    y = NULL
  ) +
  THEME_NATURE() +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggtext::element_markdown(
      family = BASE_FONT, angle = 45, hjust = 1,
      size = AXIS_TEXT_SIZE_PT, colour = "#252525"
    ),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank()
  )

p8b <- wrap_module_heatmap(p8b_heat, module_levels_heatmap, 0.14, 0.025)

ed8_diag <- diagnosis_hc3 %>%
  dplyr::filter(Module %in% MODULE_ORDER) %>%
  dplyr::mutate(
    Module = factor(Module, levels = rev(MODULE_ORDER)),
    Module_chr = as.character(Module),
    standardized_low = robust_conf_low / module_SD,
    standardized_high = robust_conf_high / module_SD,
    Significant = FDR_HC3 < 0.05,
    FDR_label = paste0("FDR ", sprintf("%.3f", FDR_HC3))
  )

diag_min <- min(c(ed8_diag$standardized_low, ed8_diag$standardized_difference), na.rm = TRUE)
diag_max <- max(c(ed8_diag$standardized_high, ed8_diag$standardized_difference), na.rm = TRUE)
diag_span <- max(diag_max - diag_min, 0.20)
diag_text_x <- diag_max + 0.08 * diag_span

p8c <- ggplot2::ggplot(
  ed8_diag,
  ggplot2::aes(y = Module, colour = Module_chr)
) +
  ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey55", linewidth = 0.35) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = standardized_low,
      xend = standardized_high,
      yend = Module
    ),
    linewidth = 0.75
  ) +
  ggplot2::geom_point(
    data = dplyr::filter(ed8_diag, !Significant),
    ggplot2::aes(x = standardized_difference),
    shape = 21,
    fill = "white",
    size = 3.0,
    stroke = 0.75
  ) +
  ggplot2::geom_point(
    data = dplyr::filter(ed8_diag, Significant),
    ggplot2::aes(x = standardized_difference, fill = Module_chr),
    shape = 21,
    size = 3.6,
    stroke = 0.75
  ) +
  ggplot2::geom_text(
    ggplot2::aes(x = diag_text_x, label = FDR_label),
    hjust = 0,
    colour = "#252525",
    family = BASE_FONT,
    size = pt_to_mm(COMPACT_TEXT_SIZE_PT)
  ) +
  ggplot2::scale_colour_manual(values = MODULE_COLORS, guide = "none") +
  ggplot2::scale_fill_manual(values = MODULE_COLORS, guide = "none") +
  module_scale_y(drop = FALSE) +
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.32))) +
  ggplot2::labs(
    title = "Diagnosis associations (HC3)",
    subtitle = "Filled points: BH-FDR < 0.05; open points: BH-FDR ≥ 0.05",
    x = "Standardized module difference: clinical AD versus CN",
    y = NULL
  ) +
  THEME_NATURE() +
  ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())

ed8_bio <- biomarker_full32 %>%
  dplyr::filter(Module %in% MODULE_ORDER, Biomarker %in% names(BIOMARKER_LABELS)) %>%
  dplyr::mutate(
    Module = factor(Module, levels = module_levels_heatmap),
    Biomarker_label = dplyr::recode(Biomarker, !!!BIOMARKER_LABELS, .default = Biomarker),
    Biomarker_label = factor(Biomarker_label, levels = unname(BIOMARKER_LABELS)),
    stars = fdr_stars(FDR_family32),
    star_colour = dplyr::if_else(abs(standardized_beta) >= 0.18, "white", "black")
  )

bio_limit <- safe_positive_limit(ed8_bio$standardized_beta)

p8d_heat <- ggplot2::ggplot(
  ed8_bio,
  ggplot2::aes(x = Biomarker_label, y = Module, fill = standardized_beta)
) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.42) +
  ggplot2::geom_text(
    ggplot2::aes(label = stars, colour = star_colour),
    family = HEATMAP_FONT,
    fontface = "bold",
    size = pt_to_mm(HEATMAP_SYMBOL_SIZE_PT)
  ) +
  ggplot2::scale_colour_identity() +
  ggplot2::scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = c(-bio_limit, bio_limit), oob = scales::squish,
    name = "Standardized beta"
  ) +
  ggplot2::labs(
    title = "Plasma biomarkers (log-HC3)",
    subtitle = "BH-FDR across the complete 8-module × 4-biomarker family (n = 32)",
    x = NULL,
    y = NULL
  ) +
  THEME_NATURE() +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank()
  )

p8d <- wrap_module_heatmap(p8d_heat, module_levels_heatmap, 0.18, 0.032)

# Compact original-style composition: full-width heatmaps followed by a paired
# diagnosis/biomarker row. Legends remain attached to their corresponding panels.
ed8_top <- make_tagged_grid(
  p8a,
  labels = "a",
  nrow = 1
)

ed8_mid <- make_tagged_grid(
  p8b,
  labels = "b",
  nrow = 1
)

ed8_bottom <- make_tagged_grid(
  p8c, p8d,
  labels = c("c", "d"),
  nrow = 1,
  rel_widths = c(1.08, 0.92)
)

ed8_figure <- cowplot::plot_grid(
  ed8_top,
  ed8_mid,
  ed8_bottom,
  ncol = 1,
  rel_heights = c(1.08, 0.86, 0.92)
)

save_panel(p8a, "Extended_Data_Figure_08", "a", 178, 70)
save_panel(p8b, "Extended_Data_Figure_08", "b", 178, 58)
save_panel(p8c, "Extended_Data_Figure_08", "c", 96, 62)
save_panel(p8d, "Extended_Data_Figure_08", "d", 82, 62)

save_master_figure(
  ed8_figure,
  "Extended_Data_Figure_08_WGCNA_module_trait_landscape",
  FIG_WIDTH_MM,
  188
)

write_source(ed8_trait, "Extended_Data_Figure_08a_complete_module_trait")
write_source(ed8_cont, "Extended_Data_Figure_08b_adjusted_continuous_models")
write_source(ed8_diag, "Extended_Data_Figure_08c_diagnosis_HC3")
write_source(ed8_bio, "Extended_Data_Figure_08d_biomarker_family32")

###############################################################################

# 7) EXTENDED DATA FIGURE 9
# Hub architecture and functional enrichment
###############################################################################

kme_module_col <- pick_col(kme_assigned, c("Module", "module"), TRUE, "module in kME table")
kme_abs_col <- pick_col(kme_assigned, c("abs_kME", "Abs_kME", "Abs_assigned_kME"), TRUE, "absolute kME")

ed9_kme <- kme_assigned %>%
  dplyr::transmute(
    Module = as.character(.data[[kme_module_col]]),
    abs_kME = safe_numeric(.data[[kme_abs_col]])
  ) %>%
  dplyr::filter(Module %in% MODULE_ORDER, is.finite(abs_kME)) %>%
  dplyr::mutate(Module = factor(Module, levels = rev(MODULE_ORDER)))

ed9_kme_n <- ed9_kme %>%
  dplyr::count(Module, name = "n_proteins") %>%
  dplyr::mutate(Module_chr = as.character(Module))

ed9_kme_median <- ed9_kme %>%
  dplyr::group_by(Module) %>%
  dplyr::summarise(abs_kME = stats::median(abs_kME, na.rm = TRUE), .groups = "drop")

kme_axis_labels <- setNames(
  vapply(rev(MODULE_ORDER), function(m) {
    n_m <- ed9_kme_n$n_proteins[ed9_kme_n$Module_chr == m]
    if (length(n_m) == 0) n_m <- NA_integer_
    paste0(MODULE_LABELS[[m]], "  (n = ", scales::comma(n_m), ")")
  }, character(1)),
  rev(MODULE_ORDER)
)

p9a <- ggplot2::ggplot(
  ed9_kme,
  ggplot2::aes(x = abs_kME, y = Module, fill = Module)
) +
  ggplot2::geom_violin(
    trim = FALSE,
    scale = "width",
    orientation = "y",
    alpha = 0.55,
    linewidth = 0.32
  ) +
  ggplot2::geom_boxplot(
    width = 0.16,
    orientation = "y",
    outlier.shape = NA,
    fill = "white",
    colour = "grey25",
    linewidth = 0.34
  ) +
  ggplot2::geom_point(
    data = ed9_kme_median,
    ggplot2::aes(x = abs_kME, y = Module),
    inherit.aes = FALSE,
    shape = 21,
    fill = "white",
    colour = "black",
    size = 1.8,
    stroke = 0.45
  ) +
  ggplot2::scale_fill_manual(values = MODULE_COLORS, guide = "none") +
  ggplot2::scale_y_discrete(labels = kme_axis_labels, drop = FALSE) +
  ggplot2::labs(
    title = "Module membership distributions",
    subtitle = "Violin density with embedded boxplot and median",
    x = "Absolute module membership (|kME|)",
    y = NULL
  ) +
  THEME_NATURE() +
  ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())

hub_count_col <- pick_col(
  hub_summary,
  c("n_hubs_abs_kME_ge_0_80", "n_hubs_kME_ge_0_80"),
  TRUE,
  "high-membership hub count"
)
hub_prop_col <- pick_col(
  hub_summary,
  c("prop_hubs_abs_kME_ge_0_80", "proportion_hubs_abs_kME_ge_0_80"),
  TRUE,
  "high-membership hub proportion"
)
hub_size_col <- pick_col(
  hub_summary,
  c("n_genes", "module_size", "n_proteins"),
  TRUE,
  "module size in hub summary"
)

ed9_hub_summary <- hub_summary %>%
  dplyr::transmute(
    Module = as.character(Module),
    N_hubs = safe_numeric(.data[[hub_count_col]]),
    Module_size = safe_numeric(.data[[hub_size_col]]),
    Prop_hubs = safe_numeric(.data[[hub_prop_col]])
  ) %>%
  dplyr::filter(Module %in% MODULE_ORDER) %>%
  dplyr::mutate(
    Module = factor(Module, levels = rev(MODULE_ORDER)),
    Module_chr = as.character(Module),
    Count_label = paste0(scales::comma(N_hubs), " / ", scales::comma(Module_size))
  )

p9b <- ggplot2::ggplot(
  ed9_hub_summary,
  ggplot2::aes(x = Prop_hubs, y = Module, colour = Module_chr)
) +
  ggplot2::geom_segment(
    ggplot2::aes(x = 0, xend = Prop_hubs, yend = Module),
    colour = "grey75",
    linewidth = 0.65
  ) +
  ggplot2::geom_point(size = 3.0) +
  ggplot2::geom_text(
    ggplot2::aes(label = Count_label),
    hjust = -0.10,
    colour = "#252525",
    family = BASE_FONT,
    size = pt_to_mm(INTERNAL_TEXT_SIZE_PT)
  ) +
  ggplot2::scale_colour_manual(values = MODULE_COLORS, guide = "none") +
  module_scale_y(drop = FALSE) +
  ggplot2::scale_x_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = ggplot2::expansion(mult = c(0, 0.26))
  ) +
  ggplot2::labs(
    title = "High-membership hub proportion",
    subtitle = "Labels show |kME| ≥ 0.80 proteins / total module proteins",
    x = "Proportion with |kME| ≥ 0.80",
    y = NULL
  ) +
  THEME_NATURE() +
  ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())

top_module_col <- pick_col(top_hubs, c("Module", "module"), TRUE, "module in top-hub table")
top_label_col <- pick_col(
  top_hubs,
  c("Protein_Display", "EntrezGeneSymbol", "GeneSymbol", "Gene", "label"),
  TRUE,
  "hub label"
)
top_abs_col <- pick_col(top_hubs, c("abs_kME", "Abs_kME", "kME_abs"), TRUE, "hub absolute kME")
top_rank_col <- pick_col(top_hubs, c("hub_rank", "Rank", "rank"), FALSE, "hub rank")

ed9_top_hubs <- top_hubs %>%
  dplyr::transmute(
    Module = as.character(.data[[top_module_col]]),
    Hub = as.character(.data[[top_label_col]]),
    abs_kME = safe_numeric(.data[[top_abs_col]]),
    hub_rank = if (!is.na(top_rank_col)) safe_numeric(.data[[top_rank_col]]) else NA_real_
  ) %>%
  dplyr::filter(Module %in% FOCAL_MODULES, is.finite(abs_kME)) %>%
  dplyr::group_by(Module) %>%
  dplyr::arrange(dplyr::desc(abs_kME), .by_group = TRUE) %>%
  dplyr::slice_head(n = 8) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    Module = factor(Module, levels = FOCAL_MODULES),
    Module_label = factor(module_label(Module), levels = module_label(FOCAL_MODULES))
  ) %>%
  dplyr::arrange(Module, abs_kME) %>%
  dplyr::mutate(
    Hub_key = paste(Hub, Module_label, sep = "___"),
    Hub_key = factor(Hub_key, levels = unique(Hub_key))
  )

hub_x_start <- max(0, min(ed9_top_hubs$abs_kME, na.rm = TRUE) - 0.025)

p9c <- ggplot2::ggplot(
  ed9_top_hubs,
  ggplot2::aes(x = abs_kME, y = Hub_key, colour = Module)
) +
  ggplot2::geom_segment(
    ggplot2::aes(x = hub_x_start, xend = abs_kME, yend = Hub_key),
    colour = "grey72",
    linewidth = 0.65
  ) +
  ggplot2::geom_point(size = 2.5) +
  ggplot2::geom_text(
    ggplot2::aes(label = sprintf("%.3f", abs_kME)),
    nudge_x = 0.004,
    hjust = 0,
    colour = "#252525",
    family = BASE_FONT,
    size = pt_to_mm(COMPACT_TEXT_SIZE_PT)
  ) +
  ggplot2::scale_colour_manual(values = MODULE_COLORS, guide = "none") +
  ggplot2::scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
  ggplot2::facet_grid(
    Module_label ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  ggplot2::coord_cartesian(xlim = c(hub_x_start, 1.015), clip = "off") +
  ggplot2::labs(
    title = "Highest-membership hubs in focal modules",
    subtitle = "Eight proteins per module; common |kME| scale",
    x = "Absolute module membership (|kME|)",
    y = NULL
  ) +
  THEME_NATURE() +
  ggplot2::theme(
    strip.placement = "outside",
    strip.text.y.left = ggplot2::element_text(angle = 0),
    panel.grid.major.y = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(3, 12, 3, 3)
  )

enrich_module_col <- pick_col(enrichment_all, c("Module", "module"), TRUE, "module in enrichment table")
enrich_library_col <- pick_col(enrichment_all, c("Library", "library", "category"), TRUE, "library")
enrich_term_col <- pick_col(enrichment_all, c("Description", "Term", "Pathway"), TRUE, "enrichment term")
enrich_fdr_col <- pick_col(enrichment_all, c("p.adjust", "FDR", "adj.P.Val", "qvalue"), TRUE, "enrichment FDR")
enrich_count_col <- pick_col(enrichment_all, c("Count", "Gene_count", "n_genes", "setSize"), FALSE, "enrichment gene count")

ed9_enrichment <- enrichment_all %>%
  dplyr::transmute(
    Module = as.character(.data[[enrich_module_col]]),
    Library = as.character(.data[[enrich_library_col]]),
    Term = as.character(.data[[enrich_term_col]]),
    FDR = safe_numeric(.data[[enrich_fdr_col]]),
    Count = if (!is.na(enrich_count_col)) safe_numeric(.data[[enrich_count_col]]) else 1
  ) %>%
  dplyr::filter(Module %in% MODULE_ORDER, is.finite(FDR), FDR < 0.05)

ed9_counts <- enrichment_summary %>%
  dplyr::filter(Module %in% MODULE_ORDER) %>%
  dplyr::select(Module, GO_BP_terms, KEGG_terms, Reactome_terms) %>%
  tidyr::pivot_longer(
    cols = -Module,
    names_to = "Library",
    values_to = "N_significant_terms"
  ) %>%
  dplyr::mutate(
    Library = dplyr::recode(
      Library,
      GO_BP_terms = "GO BP",
      KEGG_terms = "KEGG",
      Reactome_terms = "Reactome"
    ),
    Module = factor(Module, levels = rev(MODULE_ORDER)),
    Library = factor(Library, levels = c("GO BP", "KEGG", "Reactome"))
  )

library_palette <- c(
  "GO BP" = "#86A7C4",
  "GO" = "#86A7C4",
  "KEGG" = "#E7A76F",
  "Reactome" = "#74B8A6"
)

p9d <- ggplot2::ggplot(
  ed9_counts,
  ggplot2::aes(x = N_significant_terms, y = Module, fill = Library)
) +
  ggplot2::geom_col(width = 0.68) +
  ggplot2::scale_fill_manual(values = library_palette) +
  module_scale_y(drop = FALSE) +
  ggplot2::labs(
    title = "Functional enrichment burden",
    subtitle = "Descriptive BH-FDR-significant term count; not normalized biological strength",
    x = "FDR-significant terms",
    y = NULL,
    fill = NULL
  ) +
  THEME_NATURE() +
  ggplot2::theme(
    panel.grid.major.y = ggplot2::element_blank(),
    legend.position = "bottom",
    plot.subtitle = ggplot2::element_text(size = AXIS_TEXT_SIZE_PT - 0.2)
  )

ed9_top_terms <- ed9_enrichment %>%
  dplyr::filter(Module %in% FOCAL_MODULES) %>%
  dplyr::group_by(Module) %>%
  dplyr::arrange(FDR, .by_group = TRUE) %>%
  dplyr::distinct(Term, .keep_all = TRUE) %>%
  dplyr::slice_head(n = 3) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    Module = factor(Module, levels = FOCAL_MODULES),
    Module_label = factor(module_label(Module), levels = module_label(FOCAL_MODULES)),
    Neg_log10_FDR = -log10(pmax(FDR, .Machine$double.xmin)),
    Term_wrapped = wrap_term(Term, 27)
  ) %>%
  dplyr::arrange(Module, Neg_log10_FDR) %>%
  dplyr::mutate(
    Term_key = paste(Term_wrapped, Module_label, sep = "___"),
    Term_key = factor(Term_key, levels = unique(Term_key))
  )

p9e <- ggplot2::ggplot(
  ed9_top_terms,
  ggplot2::aes(x = Neg_log10_FDR, y = Term_key, colour = Module, size = Count)
) +
  ggplot2::geom_vline(
    xintercept = -log10(0.05),
    linetype = 2,
    colour = "grey55",
    linewidth = 0.35
  ) +
  ggplot2::geom_segment(
    ggplot2::aes(x = 0, xend = Neg_log10_FDR, yend = Term_key),
    colour = "grey70",
    linewidth = 0.55
  ) +
  ggplot2::geom_point(alpha = 0.95) +
  ggplot2::scale_colour_manual(values = MODULE_COLORS, guide = "none") +
  ggplot2::scale_size_continuous(range = c(2.0, 5.2), name = "Genes in term") +
  ggplot2::scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
  ggplot2::facet_grid(
    Module_label ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  ggplot2::labs(
    title = "Leading enrichment terms in focal modules",
    subtitle = "Three terms per module; dashed line marks BH-FDR = 0.05",
    x = "−log10(FDR)",
    y = NULL
  ) +
  THEME_NATURE() +
  ggplot2::theme(
    strip.placement = "outside",
    strip.text.y.left = ggplot2::element_text(angle = 0),
    panel.grid.major.y = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_text(size = AXIS_TEXT_SIZE_PT - 0.4, lineheight = 0.90),
    legend.position = "bottom",
    plot.margin = ggplot2::margin(3, 5, 3, 3)
  )

# Compact original-style composition:
# top row = distribution and hub proportion;
# lower block = focal hubs at left and enrichment summaries stacked at right.
ed9_top <- make_tagged_grid(
  p9a, p9b,
  labels = c("a", "b"),
  nrow = 1,
  rel_widths = c(1.04, 0.96)
)

ed9_right <- make_tagged_grid(
  p9d, p9e,
  labels = c("d", "e"),
  ncol = 1,
  rel_heights = c(0.72, 1.28)
)

ed9_bottom <- make_tagged_grid(
  p9c, ed9_right,
  labels = c("c", ""),
  nrow = 1,
  rel_widths = c(0.94, 1.06)
)

ed9_figure <- cowplot::plot_grid(
  ed9_top,
  ed9_bottom,
  ncol = 1,
  rel_heights = c(0.76, 1.55)
)

save_panel(p9a, "Extended_Data_Figure_09", "a", 92, 58)
save_panel(p9b, "Extended_Data_Figure_09", "b", 84, 58)
save_panel(p9c, "Extended_Data_Figure_09", "c", 84, 124)
save_panel(p9d, "Extended_Data_Figure_09", "d", 92, 52)
save_panel(p9e, "Extended_Data_Figure_09", "e", 92, 84)

save_master_figure(
  ed9_figure,
  "Extended_Data_Figure_09_WGCNA_hubs_enrichment",
  FIG_WIDTH_MM,
  208
)

write_source(ed9_kme, "Extended_Data_Figure_09a_kME_distribution")
write_source(ed9_hub_summary, "Extended_Data_Figure_09b_high_membership_hub_proportion")
write_source(ed9_top_hubs, "Extended_Data_Figure_09c_top_hubs")
write_source(ed9_counts, "Extended_Data_Figure_09d_enrichment_counts")
write_source(ed9_top_terms, "Extended_Data_Figure_09e_top_enrichment_terms")

###############################################################################

# 8) EXTENDED DATA FIGURE 10
# Recruitment-context dependence and internal association stability
###############################################################################

ed10_context <- context_country %>%
  dplyr::select(Module, Country_partial_R2 = country_partial_R2) %>%
  dplyr::left_join(
    context_nested %>%
      dplyr::select(Module, Nested_site_partial_R2 = partial_R2_site_within_country),
    by = "Module"
  ) %>%
  dplyr::filter(Module %in% MODULE_ORDER) %>%
  tidyr::pivot_longer(
    cols = c(Country_partial_R2, Nested_site_partial_R2),
    names_to = "Context",
    values_to = "Partial_R2"
  ) %>%
  dplyr::mutate(
    Context = dplyr::recode(
      Context,
      Country_partial_R2 = "Country",
      Nested_site_partial_R2 = "Site within country"
    ),
    Module = factor(Module, levels = rev(MODULE_ORDER)),
    R2_label = scales::percent(Partial_R2, accuracy = 0.1)
  )

p10a <- ggplot2::ggplot(
  ed10_context,
  ggplot2::aes(x = Partial_R2, y = Module, colour = Context, shape = Context)
) +
  ggplot2::geom_line(ggplot2::aes(group = Module), colour = "grey75", linewidth = 0.45) +
  ggplot2::geom_point(size = 2.8) +
  ggplot2::geom_text(
    data = dplyr::filter(ed10_context, Context == "Country"),
    ggplot2::aes(label = R2_label),
    nudge_x = 0.007,
    nudge_y = 0.13,
    hjust = 0,
    family = BASE_FONT,
    size = pt_to_mm(COMPACT_TEXT_SIZE_PT),
    show.legend = FALSE
  ) +
  ggplot2::geom_text(
    data = dplyr::filter(ed10_context, Context == "Site within country"),
    ggplot2::aes(label = R2_label),
    nudge_x = -0.007,
    nudge_y = -0.13,
    hjust = 1,
    family = BASE_FONT,
    size = pt_to_mm(COMPACT_TEXT_SIZE_PT),
    show.legend = FALSE
  ) +
  ggplot2::scale_x_continuous(
    labels = scales::percent,
    expand = ggplot2::expansion(mult = c(0.05, 0.10))
  ) +
  ggplot2::scale_colour_manual(values = c("Country" = "#355C7D", "Site within country" = "#C06C84")) +
  module_scale_y(drop = FALSE) +
  ggplot2::labs(
    title = "Recruitment-context contribution to eigengene variation",
    subtitle = "Adjusted for diagnosis, age, sex and education",
    x = "Partial R²",
    y = NULL,
    colour = NULL,
    shape = NULL
  ) +
  THEME_NATURE() +
  ggplot2::theme(
    panel.grid.major.y = ggplot2::element_blank(),
    legend.position = "bottom"
  )

stability_trait_use <- SELECTED_STABILITY_TRAITS[
  SELECTED_STABILITY_TRAITS %in% loco$Trait &
    SELECTED_STABILITY_TRAITS %in% downsampling$Trait
]

make_stability_heatmap <- function(tbl, value_col, title, fill_type = c("rate", "delta"), subtitle = NULL) {
  fill_type <- match.arg(fill_type)

  d <- tbl %>%
    dplyr::filter(Module %in% MODULE_ORDER, Trait %in% stability_trait_use) %>%
    dplyr::mutate(
      Module = factor(Module, levels = rev(MODULE_ORDER)),
      Trait_label = dplyr::recode(Trait, !!!TRAIT_LABELS, .default = Trait),
      Trait_label = factor(Trait_label, levels = unname(TRAIT_LABELS[stability_trait_use])),
      Value = safe_numeric(.data[[value_col]])
    )

  p <- ggplot2::ggplot(d, ggplot2::aes(x = Trait_label, y = Module, fill = Value)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.38) +
    module_scale_y(drop = FALSE) +
    ggplot2::scale_x_discrete(labels = format_trait_markdown, drop = FALSE) +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    THEME_NATURE() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggtext::element_markdown(
        family = BASE_FONT, angle = 45, hjust = 1,
        size = AXIS_TEXT_SIZE_PT, colour = "#252525"
      )
    )

  if (fill_type == "rate") {
    p <- p + ggplot2::scale_fill_gradient(
      low = "#F7FBFF", high = "#2166AC",
      limits = c(0, 1), oob = scales::squish,
      labels = scales::percent,
      name = "Rate"
    )
  } else {
    upper <- safe_positive_limit(d$Value)
    p <- p + ggplot2::scale_fill_gradient(
      low = "white", high = "#B2182B",
      limits = c(0, upper), oob = scales::squish,
      name = "Mean |Δρ|"
    )
  }

  list(plot = p, data = d)
}

h10b <- make_stability_heatmap(
  loco,
  "direction_consistency",
  "LOCO direction consistency",
  "rate"
)

h10c <- make_stability_heatmap(
  loco,
  "mean_abs_delta_rho",
  "LOCO mean absolute effect change",
  "delta"
)

h10d <- make_stability_heatmap(
  downsampling,
  "direction_consistency",
  "Balanced-downsampling direction consistency",
  "rate"
)

h10e <- make_stability_heatmap(
  downsampling,
  "FDR_significance_rate",
  "Balanced-downsampling FDR-significance rate",
  "rate",
  "Descriptive significance-retention frequency"
)

p10b <- h10b$plot
p10c <- h10c$plot
p10d <- h10d$plot
p10e <- h10e$plot

ed10_diag_base <- diagnosis_stability %>%
  dplyr::filter(Module %in% MODULE_ORDER) %>%
  dplyr::transmute(
    Module,
    Full = full_standardized_difference,
    Full_FDR = full_FDR,
    LOCO_low = loco_effect_2_5,
    LOCO_high = loco_effect_97_5,
    WC_low = wc_loso_effect_2_5,
    WC_high = wc_loso_effect_97_5,
    Down_low = downsampling_effect_2_5,
    Down_high = downsampling_effect_97_5
  )

ed10_diag <- dplyr::bind_rows(
  ed10_diag_base %>%
    dplyr::transmute(Module, Analysis = "Full sample", Estimate = Full, Low = Full, High = Full),
  ed10_diag_base %>%
    dplyr::transmute(Module, Analysis = "LOCO", Estimate = (LOCO_low + LOCO_high) / 2, Low = LOCO_low, High = LOCO_high),
  ed10_diag_base %>%
    dplyr::transmute(Module, Analysis = "Within-country LOSO", Estimate = (WC_low + WC_high) / 2, Low = WC_low, High = WC_high),
  ed10_diag_base %>%
    dplyr::transmute(Module, Analysis = "Balanced downsampling", Estimate = (Down_low + Down_high) / 2, Low = Down_low, High = Down_high)
) %>%
  dplyr::mutate(
    Module = factor(Module, levels = rev(MODULE_ORDER)),
    Analysis = factor(
      Analysis,
      levels = c("Full sample", "LOCO", "Within-country LOSO", "Balanced downsampling")
    )
  )

analysis_cols <- c(
  "Full sample" = "#252525",
  "LOCO" = "#3B78A3",
  "Within-country LOSO" = "#C06C84",
  "Balanced downsampling" = "#7B6D9C"
)
analysis_shapes <- c(
  "Full sample" = 18,
  "LOCO" = 16,
  "Within-country LOSO" = 17,
  "Balanced downsampling" = 15
)

p10f <- ggplot2::ggplot(
  ed10_diag,
  ggplot2::aes(x = Estimate, y = Module, colour = Analysis, shape = Analysis)
) +
  ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey55", linewidth = 0.35) +
  ggplot2::geom_errorbarh(
    ggplot2::aes(xmin = Low, xmax = High),
    height = 0,
    position = ggplot2::position_dodge(width = 0.62),
    linewidth = 0.62
  ) +
  ggplot2::geom_point(
    position = ggplot2::position_dodge(width = 0.62),
    size = 2.3
  ) +
  ggplot2::scale_colour_manual(values = analysis_cols) +
  ggplot2::scale_shape_manual(values = analysis_shapes) +
  module_scale_y(drop = FALSE) +
  ggplot2::labs(
    title = "Adjusted diagnosis association stability",
    subtitle = "Intervals summarize empirical sensitivity ranges; full-sample estimates are shown as points",
    x = "Standardized module difference: clinical AD versus CN",
    y = NULL,
    colour = NULL,
    shape = NULL
  ) +
  THEME_NATURE() +
  ggplot2::theme(
    panel.grid.major.y = ggplot2::element_blank(),
    legend.position = "bottom"
  )

# Compact paired 3 × 2 composition, matching the preferred original geometry.
p10f_compact <- p10f +
  ggplot2::guides(
    colour = ggplot2::guide_legend(nrow = 2, byrow = TRUE),
    shape = ggplot2::guide_legend(nrow = 2, byrow = TRUE)
  ) +
  ggplot2::theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal"
  )

ed10_row1 <- make_tagged_grid(
  p10a, p10b,
  labels = c("a", "b"),
  nrow = 1,
  rel_widths = c(0.96, 1.04)
)

ed10_row2 <- make_tagged_grid(
  p10c, p10d,
  labels = c("c", "d"),
  nrow = 1,
  rel_widths = c(1, 1)
)

ed10_row3 <- make_tagged_grid(
  p10e, p10f_compact,
  labels = c("e", "f"),
  nrow = 1,
  rel_widths = c(0.96, 1.04)
)

ed10_figure <- cowplot::plot_grid(
  ed10_row1,
  ed10_row2,
  ed10_row3,
  ncol = 1,
  rel_heights = c(0.94, 1.02, 1.07)
)

save_panel(p10a, "Extended_Data_Figure_10", "a", 86, 64)
save_panel(p10b, "Extended_Data_Figure_10", "b", 92, 68)
save_panel(p10c, "Extended_Data_Figure_10", "c", 88, 68)
save_panel(p10d, "Extended_Data_Figure_10", "d", 88, 68)
save_panel(p10e, "Extended_Data_Figure_10", "e", 88, 72)
save_panel(p10f_compact, "Extended_Data_Figure_10", "f", 92, 72)

save_master_figure(
  ed10_figure,
  "Extended_Data_Figure_10_WGCNA_context_stability",
  FIG_WIDTH_MM,
  202
)

write_source(ed10_context, "Extended_Data_Figure_10a_context_partial_R2")
write_source(h10b$data, "Extended_Data_Figure_10b_LOCO_direction")
write_source(h10c$data, "Extended_Data_Figure_10c_LOCO_delta")
write_source(h10d$data, "Extended_Data_Figure_10d_downsampling_direction")
write_source(h10e$data, "Extended_Data_Figure_10e_downsampling_FDR_rate")
write_source(ed10_diag, "Extended_Data_Figure_10f_diagnosis_stability")

###############################################################################

# 9) EXTENDED DATA FIGURE 11
# Fixed-gene structural module preservation
###############################################################################

PRESERVATION_COLORS <- c(
  "No evidence (<2)" = "#B2182B",
  "Weak–moderate (2–<10)" = "#E6B84C",
  "Strong (≥10)" = "#2166AC"
)

classify_preservation <- function(z) {
  dplyr::case_when(
    is.na(z) ~ NA_character_,
    z < 2 ~ "No evidence (<2)",
    z < 10 ~ "Weak–moderate (2–<10)",
    TRUE ~ "Strong (≥10)"
  )
}

ed11_country <- country_preservation %>%
  dplyr::filter(Module %in% MODULE_ORDER) %>%
  dplyr::mutate(
    Module = factor(Module, levels = rev(MODULE_ORDER)),
    Held_out_country = factor(
      Held_out_country,
      levels = c("Argentina", "Chile", "Colombia", "Mexico", "Peru")
    ),
    Zsummary = safe_numeric(.data[["Zsummary.pres"]]),
    Preservation = factor(
      classify_preservation(Zsummary),
      levels = names(PRESERVATION_COLORS)
    ),
    Text_colour = dplyr::if_else(Preservation == "Weak–moderate (2–<10)", "black", "white")
  )

p11a <- ggplot2::ggplot(
  ed11_country,
  ggplot2::aes(x = Held_out_country, y = Module, fill = Preservation)
) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.42) +
  ggplot2::geom_text(
    ggplot2::aes(label = sprintf("%.1f", Zsummary), colour = Text_colour),
    family = BASE_FONT,
    size = pt_to_mm(COMPACT_TEXT_SIZE_PT)
  ) +
  ggplot2::scale_fill_manual(values = PRESERVATION_COLORS, drop = FALSE, name = "Preservation") +
  ggplot2::scale_colour_identity() +
  module_scale_y(drop = FALSE) +
  ggplot2::labs(
    title = "Country-held-out structural preservation",
    subtitle = "Fixed 4,239-gene set; 100 permutations; internal structural reproducibility",
    x = "Held-out test country",
    y = NULL
  ) +
  THEME_NATURE() +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 25, hjust = 1),
    legend.position = "bottom"
  )

site_lookup <- site_preservation %>%
  dplyr::select(Country, Reference_site, Test_site) %>%
  tidyr::pivot_longer(
    cols = c(Reference_site, Test_site),
    names_to = "Position",
    values_to = "Site"
  ) %>%
  dplyr::distinct(Country, Site) %>%
  dplyr::group_by(Country) %>%
  dplyr::arrange(Site, .by_group = TRUE) %>%
  dplyr::mutate(Site_number = dplyr::row_number()) %>%
  dplyr::ungroup()

ed11_site <- site_preservation %>%
  dplyr::filter(Module %in% MODULE_ORDER) %>%
  dplyr::left_join(
    site_lookup %>% dplyr::rename(Reference_site = Site, Reference_number = Site_number),
    by = c("Country", "Reference_site")
  ) %>%
  dplyr::left_join(
    site_lookup %>% dplyr::rename(Test_site = Site, Test_number = Site_number),
    by = c("Country", "Test_site")
  ) %>%
  dplyr::mutate(
    Module = factor(Module, levels = rev(MODULE_ORDER)),
    Site_label = paste0(Country, " site ", Reference_number, " → site ", Test_number),
    Zsummary = safe_numeric(.data[["Zsummary.pres"]]),
    Preservation = factor(
      classify_preservation(Zsummary),
      levels = names(PRESERVATION_COLORS)
    ),
    Text_colour = dplyr::if_else(Preservation == "Weak–moderate (2–<10)", "black", "white")
  )

p11b <- ggplot2::ggplot(
  ed11_site,
  ggplot2::aes(x = Site_label, y = Module, fill = Preservation)
) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.42) +
  ggplot2::geom_text(
    ggplot2::aes(label = sprintf("%.1f", Zsummary), colour = Text_colour),
    family = BASE_FONT,
    size = pt_to_mm(COMPACT_TEXT_SIZE_PT)
  ) +
  ggplot2::scale_fill_manual(values = PRESERVATION_COLORS, drop = FALSE, name = "Preservation") +
  ggplot2::scale_colour_identity() +
  module_scale_y(drop = FALSE) +
  ggplot2::labs(
    title = "Reciprocal within-country site preservation",
    subtitle = "Same fixed gene set and permutation design as panel a",
    x = NULL,
    y = NULL
  ) +
  THEME_NATURE() +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 25, hjust = 1),
    legend.position = "none"
  )

ed11_country_summary <- ed11_country %>%
  dplyr::group_by(Module) %>%
  dplyr::summarise(
    Minimum_Zsummary = safe_minimum(Zsummary),
    Strong_tests = sum(Zsummary >= 10, na.rm = TRUE),
    Moderate_or_strong_tests = sum(Zsummary >= 2, na.rm = TRUE),
    Total_tests = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Module = factor(as.character(Module), levels = rev(MODULE_ORDER)),
    Summary_label = paste0(Strong_tests, "/", Total_tests, " strong")
  )

ed11_site_summary <- ed11_site %>%
  dplyr::group_by(Module) %>%
  dplyr::summarise(
    Minimum_Zsummary = safe_minimum(Zsummary),
    Strong_tests = sum(Zsummary >= 10, na.rm = TRUE),
    Moderate_or_strong_tests = sum(Zsummary >= 2, na.rm = TRUE),
    Total_tests = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Module = factor(as.character(Module), levels = rev(MODULE_ORDER)),
    Summary_label = paste0(
      Strong_tests, "/", Total_tests, " strong; ",
      Moderate_or_strong_tests, "/", Total_tests, " ≥2"
    )
  )

make_preservation_summary <- function(tbl, title) {
  x_max <- max(tbl$Minimum_Zsummary, na.rm = TRUE)
  x_span <- max(x_max, 10)
  text_x <- x_max + 0.08 * x_span
  plot_xmax <- text_x + 0.38 * x_span

  ggplot2::ggplot(
    tbl,
    ggplot2::aes(x = Minimum_Zsummary, y = Module, colour = Module)
  ) +
    ggplot2::geom_vline(xintercept = 2, linetype = 3, colour = "grey55", linewidth = 0.35) +
    ggplot2::geom_vline(xintercept = 10, linetype = 2, colour = "grey45", linewidth = 0.35) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = Minimum_Zsummary, yend = Module),
      colour = "grey78",
      linewidth = 0.55
    ) +
    ggplot2::geom_point(size = 3.0) +
    ggplot2::geom_text(
      ggplot2::aes(x = text_x, label = Summary_label),
      hjust = 0,
      colour = "#252525",
      family = BASE_FONT,
      size = pt_to_mm(COMPACT_TEXT_SIZE_PT)
    ) +
    ggplot2::scale_colour_manual(values = MODULE_COLORS, guide = "none") +
    module_scale_y(drop = FALSE) +
    ggplot2::scale_x_continuous(limits = c(0, plot_xmax), expand = ggplot2::expansion(mult = c(0, 0.01))) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(
      title = title,
      subtitle = "Thresholds: <2 no evidence; 2–<10 weak–moderate; ≥10 strong",
      x = "Minimum Zsummary across comparisons",
      y = NULL
    ) +
    THEME_NATURE() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(3, 10, 3, 3)
    )
}

p11c <- make_preservation_summary(
  ed11_country_summary,
  "Minimum country-held-out preservation by module"
)

p11d <- make_preservation_summary(
  ed11_site_summary,
  "Minimum reciprocal site preservation by module"
)

# Compact original-style composition: two full-width heatmaps and paired
# preservation summaries below.
ed11_top <- make_tagged_grid(
  p11a,
  labels = "a",
  nrow = 1
)

ed11_mid <- make_tagged_grid(
  p11b,
  labels = "b",
  nrow = 1
)

ed11_bottom <- make_tagged_grid(
  p11c, p11d,
  labels = c("c", "d"),
  nrow = 1,
  rel_widths = c(1, 1)
)

ed11_figure <- cowplot::plot_grid(
  ed11_top,
  ed11_mid,
  ed11_bottom,
  ncol = 1,
  rel_heights = c(0.94, 0.88, 0.90)
)

save_panel(p11a, "Extended_Data_Figure_11", "a", 178, 64)
save_panel(p11b, "Extended_Data_Figure_11", "b", 178, 58)
save_panel(p11c, "Extended_Data_Figure_11", "c", 90, 60)
save_panel(p11d, "Extended_Data_Figure_11", "d", 90, 60)

save_master_figure(
  ed11_figure,
  "Extended_Data_Figure_11_WGCNA_structural_preservation",
  FIG_WIDTH_MM,
  186
)

write_source(ed11_country, "Extended_Data_Figure_11a_country_preservation")
write_source(ed11_site, "Extended_Data_Figure_11b_site_preservation")
write_source(ed11_country_summary, "Extended_Data_Figure_11c_country_summary")
write_source(ed11_site_summary, "Extended_Data_Figure_11d_site_summary")

###############################################################################

# 10) TYPOGRAPHY AND ARTBOARD AUDIT
###############################################################################

typography_audit <- tibble::tibble(
  item = c(
    "Panel tags",
    "Panel titles",
    "Axis titles",
    "Axis text",
    "Legend titles",
    "Legend text",
    "Internal/protein labels",
    "Compact numeric text",
    "Heatmap symbols",
    "Master width"
  ),
  expected = c(9, 7, 7, 6, 6, 6, 6, 5, 6, 180),
  observed = c(
    TAG_SIZE_PT,
    TITLE_SIZE_PT,
    AXIS_TITLE_SIZE_PT,
    AXIS_TEXT_SIZE_PT,
    LEGEND_TITLE_SIZE_PT,
    LEGEND_TEXT_SIZE_PT,
    INTERNAL_TEXT_SIZE_PT,
    COMPACT_TEXT_SIZE_PT,
    HEATMAP_SYMBOL_SIZE_PT,
    FIG_WIDTH_MM
  )
)

if (any(typography_audit$expected != typography_audit$observed)) {
  stop("WGCNA Extended Data typography audit failed.", call. = FALSE)
}

readr::write_csv(
  typography_audit,
  file.path(DIAGNOSTICS_DIR, "WGCNA_Extended_Data_typography_audit_v10.csv")
)

font_safety_audit <- tibble::tibble(
  check = c(
    "Base font is Arial",
    "Heatmap font uses Arial",
    "ggtext is available",
    "APOE markdown helper is active",
    "Arial is found by systemfonts",
    "rsvg PDF conversion is available",
    "Plotmath labels are disabled"
  ),
  passed = c(
    identical(BASE_FONT, "Arial"),
    identical(HEATMAP_FONT, BASE_FONT),
    requireNamespace("ggtext", quietly = TRUE),
    identical(format_trait_markdown("APOE ε4"), "<i>APOE</i> ε4"),
    nzchar(systemfonts::match_font("Arial")$path),
    is.function(rsvg::rsvg_pdf),
    isTRUE(PLOTMATH_LABELS_DISABLED)
  )
)

if (!all(font_safety_audit$passed)) {
  stop(
    "WGCNA Extended Data font-safety audit failed: ",
    paste(font_safety_audit$check[!font_safety_audit$passed], collapse = "; "),
    call. = FALSE
  )
}

readr::write_csv(
  font_safety_audit,
  file.path(DIAGNOSTICS_DIR, "WGCNA_Extended_Data_font_safety_audit_v10.csv")
)

###############################################################################
# 11) FIGURE LEGENDS, INVENTORY AND MANUSCRIPT-REFERENCE AUDIT
###############################################################################

figure_legends <- tibble::tribble(
  ~Figure_ID, ~Title, ~Legend,
  "Extended Data Fig. 7",
  "Network construction, quality and global organization of the outcome-independent WGCNA network.",
  paste(
    "a, Signed scale-free topology fit across candidate soft-thresholding powers. b, Mean network connectivity across candidate powers; the selected power was β = 4.",
    "c, Distribution of absolute assigned module-membership values (|kME|), shown as violin densities with embedded boxplots and medians.",
    "d, Mean within-module versus external adjacency and topological-overlap ratios on a common logarithmic axis; the dashed line marks a ratio of one.",
    "e, Correlations among final module eigengenes. f, Module-specific contributions to the exact post hoc weighted Newman–Girvan modularity Q.",
    "The original sample-clustering diagnostic remains available among the Script 11 quality-control outputs. Modularity was calculated descriptively after module definition and was not used to optimize the WGCNA solution. M1–M8 are stable visual identifiers and do not represent a ranking."
  ),

  "Extended Data Fig. 8",
  "Complete module–trait landscape and covariate-adjusted clinical and plasma biomarker associations.",
  paste(
    "a, Spearman correlations between all eight module eigengenes and the 16 prespecified clinical, cognitive, demographic, genetic and plasma biomarker traits. Asterisks denote BH-FDR significance across the complete 128-test matrix.",
    "b, Covariate-adjusted continuous-outcome models across all modules and six outcomes; colour indicates the signed t statistic and asterisks denote BH-FDR significance across 48 models.",
    "c, Adjusted diagnosis-related module differences using HC3 robust inference; points show standardized differences and horizontal lines show 95% confidence intervals. Filled points reached BH-FDR < 0.05 and open points did not; exact FDR values are displayed.",
    "d, Adjusted log-transformed plasma biomarker models with HC3 inference. Colour indicates standardized beta and asterisks denote BH-FDR significance across the complete 32-model family.",
    "M1–M8 are stable visual identifiers and do not represent a ranking."
  ),

  "Extended Data Fig. 9",
  "Hub architecture and functional enrichment across WGCNA modules.",
  paste(
    "a, Horizontal violin distributions of absolute module membership values with embedded boxplots and medians; labels report the number of assigned proteins.",
    "b, Proportion of proteins with |kME| ≥ 0.80 in each module; labels give the corresponding count and module denominator.",
    "c, Eight highest-membership hubs in the green, blue and brown focal modules on a common |kME| scale.",
    "d, Descriptive number of BH-FDR-significant enrichment terms across GO Biological Process, KEGG and Reactome; term counts are not normalized measures of biological strength.",
    "e, Three leading BH-FDR-significant enrichment terms per focal module; point size represents the number of mapped genes and the dashed line marks BH-FDR = 0.05.",
    "M1–M8 are stable visual identifiers and do not represent a ranking."
  ),

  "Extended Data Fig. 10",
  "Recruitment-context dependence and internal stability of WGCNA associations.",
  paste(
    "a, Partial R² associated with country and with recruitment site nested within country after adjustment for diagnosis, age, sex and education.",
    "b,c, Leave-one-country-out directional consistency and mean absolute change in selected module–trait Spearman correlations.",
    "d,e, Directional consistency and descriptive FDR-significance-retention rates across 500 country-by-diagnosis balanced-downsampling iterations; significance-retention frequencies are sensitive to sample size and power.",
    "f, Full-sample adjusted diagnosis effects and empirical ranges across leave-one-country-out, corrected within-country site-deletion and balanced-downsampling analyses. Shape and colour distinguish analysis layers.",
    "These analyses assess internal association stability while holding module definitions fixed and are distinct from structural module preservation. M1–M8 are stable visual identifiers and do not represent a ranking."
  ),

  "Extended Data Fig. 11",
  "Structural preservation of WGCNA modules across countries and recruitment sites using a fixed gene set.",
  paste(
    "a, WGCNA Zsummary statistics for each module when one country was held out as the test dataset. b, Reciprocal structural preservation between the two Chilean recruitment sites and between the two Colombian recruitment sites.",
    "The same fixed set of 4,239 genes and 100 permutations was used in every comparison. Tile colours encode the prespecified interpretive categories and exact Zsummary values are printed within cells.",
    "c,d, Minimum Zsummary and the number of strong or weak-to-moderate-or-strong comparisons by module across country and site analyses, respectively.",
    "Zsummary ≥ 10 was interpreted descriptively as strong preservation, 2 to <10 as weak-to-moderate preservation and <2 as no evidence of preservation. These analyses represent internal structural reproducibility rather than external validation. M1–M8 are stable visual identifiers and do not represent a ranking."
  )
)
readr::write_csv(
  figure_legends,
  file.path(DIAGNOSTICS_DIR, "WGCNA_Extended_Data_figure_legends.csv")
)

module_identifier_key <- tibble::tibble(
  Module_ID = paste0("M", seq_along(MODULE_ORDER)),
  WGCNA_module = MODULE_ORDER,
  Display_label = unname(MODULE_LABELS[MODULE_ORDER]),
  Hex_colour = unname(MODULE_COLORS[MODULE_ORDER]),
  Interpretation = "Stable visual identifier only; M1–M8 do not represent a ranking."
)

readr::write_csv(
  module_identifier_key,
  file.path(DIAGNOSTICS_DIR, "WGCNA_module_identifier_key.csv")
)

legend_md <- purrr::pmap_chr(
  figure_legends,
  function(Figure_ID, Title, Legend) {
    paste0("## ", Figure_ID, ". ", Title, "\n\n", Legend, "\n")
  }
)

writeLines(
  legend_md,
  file.path(DIAGNOSTICS_DIR, "WGCNA_Extended_Data_figure_legends.md")
)

figure_inventory <- tibble::tribble(
  ~Figure_ID, ~File_stem, ~Width_mm, ~Height_mm, ~Manuscript_order,
  "Extended Data Fig. 7", "Extended_Data_Figure_07_WGCNA_network_quality", 180, 182, 1,
  "Extended Data Fig. 8", "Extended_Data_Figure_08_WGCNA_module_trait_landscape", 180, 188, 2,
  "Extended Data Fig. 9", "Extended_Data_Figure_09_WGCNA_hubs_enrichment", 180, 208, 3,
  "Extended Data Fig. 10", "Extended_Data_Figure_10_WGCNA_context_stability", 180, 202, 4,
  "Extended Data Fig. 11", "Extended_Data_Figure_11_WGCNA_structural_preservation", 180, 186, 5
) %>%
  dplyr::mutate(
    PDF = file.path(MASTER_VECTOR_DIR, paste0(File_stem, ".pdf")),
    SVG = file.path(MASTER_VECTOR_DIR, paste0(File_stem, ".svg")),
    PNG = file.path(PREVIEW_DIR, paste0(File_stem, ".png")),
    TIFF = file.path(SUBMISSION_DIR, paste0(File_stem, ".tiff")),
    All_exist = file.exists(PDF) & file.exists(SVG) & file.exists(PNG) & file.exists(TIFF),
    PDF_bytes = ifelse(file.exists(PDF), file.info(PDF)$size, NA_real_),
    SVG_bytes = ifelse(file.exists(SVG), file.info(SVG)$size, NA_real_),
    PNG_bytes = ifelse(file.exists(PNG), file.info(PNG)$size, NA_real_),
    TIFF_bytes = ifelse(file.exists(TIFF), file.info(TIFF)$size, NA_real_),
    Nonempty = All_exist & PDF_bytes > 0 & SVG_bytes > 0 & PNG_bytes > 0 & TIFF_bytes > 0
  )

readr::write_csv(
  figure_inventory,
  file.path(DIAGNOSTICS_DIR, "WGCNA_Extended_Data_figure_inventory.csv")
)

input_audit <- tibble::tibble(
  Source_ID = names(PATHS),
  Path = unlist(PATHS, use.names = FALSE),
  Exists = file.exists(unlist(PATHS, use.names = FALSE))
)

readr::write_csv(
  input_audit,
  file.path(DIAGNOSTICS_DIR, "WGCNA_Extended_Data_input_audit.csv")
)

manuscript_reference_audit <- tibble::tribble(
  ~Location, ~Current_reference_or_wording, ~Required_update,
  "Results — structural preservation paragraph",
  "Fig. 3f",
  "Replace with Extended Data Fig. 11 because final Fig. 3f is the blue-module biological portrait.",
  "Results — blue adjusted-association paragraph",
  "Fig. 3e",
  "Verify against final Fig. 3 panel order; the blue-module portrait is Fig. 3f, whereas adjusted associations are displayed in Extended Data Fig. 8.",
  "Results — WGCNA quality paragraph",
  "Extended Data Fig. 7",
  "Retain: this is the first WGCNA Extended Data figure in manuscript order.",
  "Results — complete module–trait paragraph",
  "Extended Data Fig. 8",
  "Retain and use for the complete 128-test matrix and adjusted model families.",
  "Results — focal module biology paragraphs",
  "Extended Data Fig. 9",
  "Retain for hub architecture and enrichment across focal modules.",
  "Results — recruitment-context stability paragraph",
  "Extended Data Fig. 10",
  "Retain for country/site context and internal association stability."
)

readr::write_csv(
  manuscript_reference_audit,
  file.path(DIAGNOSTICS_DIR, "WGCNA_manuscript_reference_audit.csv")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(DIAGNOSTICS_DIR, "sessionInfo.txt")
)

if (!all(figure_inventory$Nonempty)) {
  stop(
    "One or more final Extended Data outputs were not generated. Inspect: ",
    file.path(DIAGNOSTICS_DIR, "WGCNA_Extended_Data_figure_inventory.csv"),
    call. = FALSE
  )
}

message("============================================================")
message("WGCNA Extended Data figures were generated successfully.")
message("Output directory: ", OUTDIR)
message("Manuscript order: Extended Data Figs. 7–11")
message("============================================================")

}) # end local render-safe execution


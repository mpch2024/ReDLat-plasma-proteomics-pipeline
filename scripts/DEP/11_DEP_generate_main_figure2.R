###############################################################################
# ReDLat plasma proteomics — DEP workflow
# 11. Main Figure 2
# Requires: Approved DEP Source Data and supplementary outputs
# Produces: Editable and submission-ready Figure 2 files
# Data policy: participant-level data and intermediate outputs remain local.
###############################################################################

rm(list = ls())

# -----------------------------------------------------------------------------
# Project setup
# -----------------------------------------------------------------------------
.project_root_env <- Sys.getenv("REDLAT_PROJECT_ROOT", unset = "")
if (nzchar(.project_root_env)) {
  project_root <- normalizePath(.project_root_env, winslash = "/", mustWork = TRUE)
} else if (requireNamespace("here", quietly = TRUE)) {
  project_root <- normalizePath(here::here(), winslash = "/", mustWork = TRUE)
} else {
  stop("Package 'here' is required. Restore the project environment with renv::restore().", call. = FALSE)
}
source(file.path(project_root, "R", "dep_bootstrap.R"), local = FALSE)
DEP_CONFIG <- dep_load_config(project_root)
analysis_root <- DEP_CONFIG$result_root
publication_root <- DEP_CONFIG$publication_root
options(stringsAsFactors = FALSE)

message(
  "Run this file as one complete script with source(..., encoding = 'UTF-8'). ",
  "Do not continue line by line after an error."
)

required_pkgs <- c(
  "dplyr", "tidyr", "tibble", "stringr", "forcats",
  "ggplot2", "cowplot", "openxlsx", "scales", "svglite", "readr", "ggtext"
)

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0L) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
       ". Run renv::restore() before this script.", call. = FALSE)
}

invisible(lapply(required_pkgs, library, character.only = TRUE))

message("cowplot version: ", as.character(utils::packageVersion("cowplot")))

###############################################################################
# 1. Paths
###############################################################################


first_existing_dir <- function(paths) {
  paths <- unique(paths[!is.na(paths) & nzchar(paths)])
  hit <- paths[dir.exists(paths)][1]
  if (length(hit) == 0 || is.na(hit)) return(NA_character_)
  normalizePath(hit, winslash = "/", mustWork = TRUE)
}

source_root <- first_existing_dir(c(
  file.path(publication_root, "source_data")
))

tables_root <- first_existing_dir(c(
  file.path(publication_root, "supplementary_data")
))

if (is.na(source_root)) {
  stop(
    "Could not locate publication_candidate/DEP/source_data.\n",
    "Set DEP_SOURCE_DATA_DIR explicitly.",
    call. = FALSE
  )
}

if (is.na(tables_root)) {
  stop(
    "Could not locate publication_candidate/DEP/supplementary_data.\n",
    "Set DEP_FINAL_TABLES_DIR explicitly.",
    call. = FALSE
  )
}

out_root <- file.path(publication_root, "figures", "main_figure_2")

illustrator_root <- file.path(out_root, "Illustrator_ready")
master_vector_root <- file.path(illustrator_root, "master_vector")
panel_root <- file.path(illustrator_root, "panels_vector")
preview_root <- file.path(out_root, "preview")
preview_panel_root <- file.path(preview_root, "panels")

for (d in c(
  master_vector_root,
  panel_root,
  preview_root,
  preview_panel_root
)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

fig2_xlsx <- file.path(source_root, "SourceData_Fig2_DEP.xlsx")
suppdata1_xlsx <- file.path(
  tables_root, "Supplementary_Data_1_DEP_Complete_Protein_Results.xlsx"
)
representative_xlsx <- file.path(
  source_root, "SourceData_DEP_RepresentativeProteins.xlsx"
)

for (f in c(fig2_xlsx, suppdata1_xlsx, representative_xlsx)) {
  if (!file.exists(f)) stop("Required workbook not found: ", f, call. = FALSE)
}

###############################################################################
# 2. Constants and helpers
###############################################################################

COL_ORANGE <- "#F46D43"
COL_BLUE   <- "#3288BD"
COL_DARK   <- "#252525"
COL_MID    <- "#8C8C8C"
COL_LIGHT  <- "#D9D9D9"
COL_PALE   <- "#F2F2F2"
COL_WHITE  <- "#FFFFFF"

###############################################################################
# Fixed typography and artboard
#
# IMPORTANT:
# - ggplot2 theme text sizes are specified directly in points.
# - geom_text()/annotate("text") sizes are specified in mm-like units, so all
#   layer text uses pt_to_geom() for reproducible point-equivalent sizing.
# - The final artboard remains fixed at 180 mm × 170 mm. Typography changes
#   must never be compensated by increasing figure dimensions.
###############################################################################

FONT_ARIAL <- "Arial"
FONT_MYRIAD <- "Myriad Pro"

TEXT_PT <- list(
  panel_tag = 9,
  panel_title = 7,
  axis_title = 7,
  axis_tick = 6,
  protein = 6,
  legend = 6,
  legend_title = 6,
  heatmap_colourbar_tick = 5,
  compact_numeric = 5,
  micro = 4.08,
  myriad_heatmap = 6,
  myriad_exponent = 3.5
)

pt_to_geom <- function(pt) pt / ggplot2::.pt

BASE_SIZE <- TEXT_PT$axis_tick
TAG_SIZE <- TEXT_PT$panel_tag
PANEL_TITLE_SIZE <- TEXT_PT$panel_title
LINE_MM <- 0.28

MAIN_WIDTH_MM <- 180
MAIN_HEIGHT_MM <- 170

FOCAL_PROTEINS_HIGH <- c("SPC25", "LRRN1", "GJB3", "CPLX2", "SMOC1")
FOCAL_PROTEINS_LOW  <- c("C3", "PLCH1", "PUM2", "GDI1", "AMBRA1")
FOCAL_PROTEINS <- c(FOCAL_PROTEINS_HIGH, FOCAL_PROTEINS_LOW)

as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
as_chr <- function(x) trimws(as.character(x))
as_bool <- function(x) {
  if (is.logical(x)) return(replace(x, is.na(x), FALSE))
  tolower(as_chr(x)) %in% c("true", "t", "1", "yes", "y")
}
neglog10 <- function(x) -log10(pmax(as_num(x), .Machine$double.xmin))

read_sheet <- function(path, sheet) {
  available <- openxlsx::getSheetNames(path)
  if (!sheet %in% available) {
    stop(
      "Sheet '", sheet, "' not found in ", basename(path),
      "\nAvailable sheets: ", paste(available, collapse = ", "),
      call. = FALSE
    )
  }

  openxlsx::read.xlsx(
    path,
    sheet = sheet,
    startRow = 4,
    colNames = TRUE,
    check.names = FALSE,
    skipEmptyRows = TRUE,
    skipEmptyCols = TRUE
  ) |>
    tibble::as_tibble()
}

assert_cols <- function(x, required, label) {
  missing <- setdiff(required, names(x))
  if (length(missing) > 0) {
    stop(
      label, " is missing columns: ", paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

nature_theme <- function(base_size = BASE_SIZE) {
  ggplot2::theme_classic(base_size = base_size, base_family = FONT_ARIAL) +
    ggplot2::theme(
      text = ggplot2::element_text(family = FONT_ARIAL, colour = COL_DARK),
      axis.title = ggplot2::element_text(
        family = FONT_ARIAL,
        size = TEXT_PT$axis_title,
        face = "plain",
        colour = COL_DARK
      ),
      axis.text = ggplot2::element_text(
        family = FONT_ARIAL,
        size = TEXT_PT$axis_tick,
        face = "plain",
        colour = COL_DARK
      ),
      axis.line = ggplot2::element_line(linewidth = LINE_MM, colour = COL_DARK),
      axis.ticks = ggplot2::element_line(linewidth = LINE_MM, colour = COL_DARK),
      axis.ticks.length = grid::unit(1.1, "mm"),
      legend.title = ggplot2::element_text(
        family = FONT_ARIAL,
        size = TEXT_PT$legend_title,
        face = "plain"
      ),
      legend.text = ggplot2::element_text(
        family = FONT_ARIAL,
        size = TEXT_PT$legend,
        face = "plain"
      ),
      legend.key.height = grid::unit(2.5, "mm"),
      legend.key.width = grid::unit(2.8, "mm"),
      panel.grid = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      legend.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin = ggplot2::margin(0.7, 1.5, 0.9, 1.5, unit = "mm")
    )
}

# Panel titles are overlaid only AFTER the underlying plots have been aligned.
# This avoids the v10 problem in which patchwork aligned the external boxes
# rather than the actual plotting regions.
make_titled_panel <- function(
  aligned_plot,
  tag,
  title,
  legend_grob = NULL,
  plot_height = 0.950,
  side_note = NULL,
  side_note_x = 0.79,
  side_note_y = 0.16,
  side_note_size = TEXT_PT$compact_numeric
) {
  g <- cowplot::ggdraw() +
    cowplot::draw_plot(
      aligned_plot,
      x = 0,
      y = 0,
      width = 1,
      height = plot_height
    ) +
    cowplot::draw_label(
      tag,
      x = 0.010,
      y = 1.004,
      hjust = 0,
      vjust = 1,
      fontfamily = FONT_ARIAL,
      fontface = "bold",
      size = TEXT_PT$panel_tag,
      colour = COL_DARK
    ) +
    cowplot::draw_label(
      title,
      x = 0.072,
      y = 0.995,
      hjust = 0,
      vjust = 1,
      fontfamily = FONT_ARIAL,
      fontface = "plain",
      size = TEXT_PT$panel_title,
      colour = COL_DARK
    )

  if (!is.null(legend_grob)) {
    # Place the real ggplot legend below the title and centered over the
    # Reactome bar area rather than floating against the top border.
    g <- g + cowplot::draw_grob(
      legend_grob,
      x = 0.515,
      y = 0.882,
      width = 0.46,
      height = 0.075
    )
  }

  if (!is.null(side_note) && nzchar(side_note)) {
    # Compact annotation key positioned inside the right-side legend region
    # of panel c, beneath the Spearman-rho colour bar.
    g <- g + cowplot::draw_label(
      side_note,
      x = side_note_x,
      y = side_note_y,
      hjust = 0,
      vjust = 0.5,
      fontfamily = FONT_ARIAL,
      fontface = "plain",
      size = side_note_size,
      colour = COL_DARK,
      lineheight = 0.95
    )
  }

  g
}

make_direction_legend_grob <- function() {
  dummy <- tibble::tibble(
    x = c(1, 2),
    y = c(1, 1),
    direction = factor(
      c("Higher in AD", "Lower in AD"),
      levels = c("Higher in AD", "Lower in AD")
    )
  )

  legend_plot <- ggplot2::ggplot(
    dummy,
    ggplot2::aes(x = x, y = y, fill = direction)
  ) +
    ggplot2::geom_point(shape = 22, size = 2.8, stroke = 0.28, colour = COL_DARK) +
    ggplot2::scale_fill_manual(
      values = c("Higher in AD" = COL_ORANGE, "Lower in AD" = COL_BLUE)
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(shape = 22, size = 4.0, colour = COL_DARK)
      )
    ) +
    ggplot2::theme_void(base_family = FONT_ARIAL) +
    ggplot2::theme(
      legend.position = "top",
      legend.direction = "horizontal",
      legend.justification = "center",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(
        family = FONT_ARIAL, size = TEXT_PT$legend, face = "plain", colour = COL_DARK
      ),
      legend.key.width = grid::unit(4.2, "mm"),
      legend.key.height = grid::unit(3.8, "mm"),
      legend.spacing.x = grid::unit(1.4, "mm"),
      legend.margin = ggplot2::margin(0, 0, 0, 0, "mm"),
      legend.box.margin = ggplot2::margin(0, 0, 0, 0, "mm")
    )

  cowplot::get_legend(legend_plot)
}

align_row_pair <- function(left_plot, right_plot, row_name) {
  aligned <- tryCatch(
    cowplot::align_plots(
      left_plot,
      right_plot,
      align = "h",
      axis = "tb"
    ),
    error = function(e) {
      stop(
        "Cowplot could not align ", row_name, ": ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  if (length(aligned) != 2L) {
    stop("Unexpected alignment result for ", row_name, call. = FALSE)
  }

  aligned
}


###############################################################################
# 3. Safe export
###############################################################################

validate_file <- function(path, min_bytes) {
  if (!file.exists(path)) stop("Output not created: ", path, call. = FALSE)
  sz <- file.info(path)$size
  if (!is.finite(sz) || sz < min_bytes) {
    stop("Output appears blank or truncated (", sz, " bytes): ", path, call. = FALSE)
  }
  invisible(path)
}

# SVG files should be validated by structure rather than a fixed 10-KB rule.
# A relatively simple vector panel, such as the Reactome bar chart, can be
# fully valid at 8–9 KB.
validate_svg_file <- function(path, min_bytes = 1000L) {
  if (!file.exists(path)) {
    stop("SVG output not created: ", path, call. = FALSE)
  }

  sz <- file.info(path)$size
  if (!is.finite(sz) || sz < min_bytes) {
    stop(
      "SVG appears blank or truncated (", sz, " bytes): ",
      path,
      call. = FALSE
    )
  }

  svg_text <- paste(
    readLines(path, warn = FALSE, encoding = "UTF-8"),
    collapse = "\n"
  )

  has_opening_svg <- grepl(
    "<svg(?:\\s|>)",
    svg_text,
    perl = TRUE,
    ignore.case = TRUE
  )
  has_closing_svg <- grepl(
    "</svg>",
    tolower(svg_text),
    fixed = TRUE
  )
  has_dimensions <- grepl(
    "(viewBox=|width=)",
    svg_text,
    perl = TRUE,
    ignore.case = TRUE
  )
  has_graphics <- grepl(
    "<(path|text|rect|circle|line|polyline|polygon|use|g)(?:\\s|>)",
    svg_text,
    perl = TRUE,
    ignore.case = TRUE
  )

  if (!(has_opening_svg && has_closing_svg && has_dimensions && has_graphics)) {
    stop(
      "SVG failed structural validation despite having ",
      sz, " bytes: ", path,
      call. = FALSE
    )
  }

  invisible(path)
}

safe_pdf <- function(plot, path, width_mm, height_mm) {
  tmp <- paste0(path, ".tmp.pdf")
  if (file.exists(tmp)) unlink(tmp)

  grDevices::cairo_pdf(
    filename = tmp,
    width = width_mm / 25.4,
    height = height_mm / 25.4,
    family = FONT_ARIAL,
    bg = "white",
    onefile = TRUE
  )
  err <- NULL
  tryCatch(print(plot), error = function(e) err <<- e)
  grDevices::dev.off()

  if (!is.null(err)) {
    if (file.exists(tmp)) unlink(tmp)
    stop("PDF drawing failed: ", conditionMessage(err), call. = FALSE)
  }

  validate_file(tmp, 10000L)
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) {
    if (!file.copy(tmp, path, overwrite = TRUE)) {
      stop("Could not finalize PDF: ", path, call. = FALSE)
    }
    unlink(tmp)
  }
  invisible(path)
}

safe_svg <- function(plot, path, width_mm, height_mm) {
  tmp <- paste0(path, ".tmp.svg")
  if (file.exists(tmp)) unlink(tmp)

  svglite::svglite(
    file = tmp,
    width = width_mm / 25.4,
    height = height_mm / 25.4,
    bg = "white"
  )
  err <- NULL
  tryCatch(print(plot), error = function(e) err <<- e)
  grDevices::dev.off()

  if (!is.null(err)) {
    if (file.exists(tmp)) unlink(tmp)
    stop("SVG drawing failed: ", conditionMessage(err), call. = FALSE)
  }

  validate_svg_file(tmp, min_bytes = 1000L)
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) {
    if (!file.copy(tmp, path, overwrite = TRUE)) {
      stop("Could not finalize SVG: ", path, call. = FALSE)
    }
    unlink(tmp)
  }

  # Validate the finalized file once more after rename/copy.
  validate_svg_file(path, min_bytes = 1000L)
  invisible(path)
}

safe_png <- function(plot, path, width_mm, height_mm, dpi = 600L) {
  width_px <- round(width_mm / 25.4 * dpi)
  height_px <- round(height_mm / 25.4 * dpi)
  tmp <- paste0(path, ".tmp.png")
  if (file.exists(tmp)) unlink(tmp)

  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(
      tmp,
      width = width_px,
      height = height_px,
      units = "px",
      res = dpi,
      background = "white"
    )
  } else {
    grDevices::png(
      tmp,
      width = width_px,
      height = height_px,
      res = dpi,
      bg = "white",
      type = "cairo"
    )
  }

  err <- NULL
  tryCatch(print(plot), error = function(e) err <<- e)
  grDevices::dev.off()

  if (!is.null(err)) {
    if (file.exists(tmp)) unlink(tmp)
    stop("PNG drawing failed: ", conditionMessage(err), call. = FALSE)
  }

  validate_file(tmp, 20000L)
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) {
    if (!file.copy(tmp, path, overwrite = TRUE)) {
      stop("Could not finalize PNG: ", path, call. = FALSE)
    }
    unlink(tmp)
  }
  invisible(path)
}

save_panel <- function(plot, stem, width_mm = 90, height_mm = 70) {
  # Illustrator-editable vector panel.
  safe_pdf(
    plot,
    file.path(panel_root, paste0(stem, ".pdf")),
    width_mm, height_mm
  )
  safe_svg(
    plot,
    file.path(panel_root, paste0(stem, ".svg")),
    width_mm, height_mm
  )

  # Raster preview only; do not use this PNG for Illustrator editing.
  safe_png(
    plot,
    file.path(preview_panel_root, paste0(stem, "_preview.png")),
    width_mm, height_mm, dpi = 300L
  )

  svg_panel_path <- file.path(panel_root, paste0(stem, ".svg"))
  message(
    "Illustrator SVG panel created: ",
    basename(svg_panel_path),
    " (",
    file.info(svg_panel_path)$size,
    " bytes)"
  )
}

###############################################################################
# 4. Read the approved data
###############################################################################

fig2a <- read_sheet(fig2_xlsx, "Fig2a_volcano")
fig2b <- read_sheet(fig2_xlsx, "Fig2b_displayed_paths")
fig2e <- read_sheet(fig2_xlsx, "Fig2e_LOCO_summary")
fig2f <- read_sheet(fig2_xlsx, "Fig2f_mean_LOCO")

# The main Figure 2 workbook contains correlations only for the previously
# displayed protein set. The representative-protein Source Data workbook
# contains the broader 50-protein correlation landscape and supplies GJB3
# and AMBRA1 without recomputing correlations from participant-level data.
fig2c_corr_full <- read_sheet(representative_xlsx, "Heatmap_correlations")

# Rebuild annotation columns from canonical complete protein-level results.
primary_dep <- read_sheet(suppdata1_xlsx, "Primary_DEP_full")
apoe_equal <- read_sheet(suppdata1_xlsx, "APOE_equal_sample")
primary_vs_cdrsb <- read_sheet(suppdata1_xlsx, "Primary_vs_CDRSB")

# Crucial correction: same-sample baseline versus AT(N)-adjusted model.
fig2d <- read_sheet(suppdata1_xlsx, "ATN_equal_sample")

assert_cols(fig2a, c("EntrezGeneSymbol", "logFC", "adj.P.Val"), "Fig. 2a")
assert_cols(fig2b, c("Display_label", "NES", "p.adjust"), "Fig. 2b")
assert_cols(
  fig2c_corr_full,
  c("Protein", "Trait", "rho", "P.Value"),
  "Representative-protein correlations"
)
assert_cols(
  primary_dep,
  c("EntrezGeneSymbol", "AptName", "logFC", "adj.P.Val"),
  "Primary DEP results"
)
assert_cols(
  apoe_equal,
  c("EntrezGeneSymbol", "adjusted_fdr005"),
  "APOE equal-sample results"
)
assert_cols(
  primary_vs_cdrsb,
  c("EntrezGeneSymbol", "same_direction"),
  "Primary-versus-CDR-SB results"
)
assert_cols(
  fig2d,
  c(
    "EntrezGeneSymbol", "primary_full_logFC", "primary_full_fdr005",
    "subset_baseline_logFC", "subset_baseline_adj.P.Val",
    "adjusted_logFC", "adjusted_adj.P.Val"
  ),
  "Fig. 2d equal-sample AT(N)"
)
assert_cols(fig2e, c("excluded_country", "logFC_correlation"), "Fig. 2e")
assert_cols(
  fig2f,
  c("EntrezGeneSymbol", "main_logFC", "mean_loco_logFC", "main_adj.P.Val"),
  "Fig. 2f"
)

###############################################################################
# 5. Deterministic label helper
###############################################################################

make_manual_labels <- function(
  data, gene_col, x_col, y_col, offset_table,
  x_fraction = TRUE, y_fraction = TRUE
) {
  x <- as_num(data[[x_col]])
  y <- as_num(data[[y_col]])

  x_span <- diff(range(x, finite = TRUE))
  y_span <- diff(range(y, finite = TRUE))
  if (!is.finite(x_span) || x_span == 0) x_span <- 1
  if (!is.finite(y_span) || y_span == 0) y_span <- 1

  d <- data |>
    dplyr::transmute(
      gene = as_chr(.data[[gene_col]]),
      x_point = as_num(.data[[x_col]]),
      y_point = as_num(.data[[y_col]]),
      class = if ("class" %in% names(data)) {
        as_chr(.data[["class"]])
      } else {
        "Not significant"
      }
    ) |>
    dplyr::inner_join(offset_table, by = "gene") |>
    dplyr::mutate(
      label_x = x_point + if (x_fraction) dx * x_span else dx,
      label_y = y_point + if (y_fraction) dy * y_span else dy,
      hjust = dplyr::if_else(dx >= 0, 0, 1)
    )

  d
}

# Protein labels use a white text halo rather than rectangular boxes.
# The highlighted point retains its original orange/blue direction.
add_manual_labels <- function(
  p, labels,
  label_size = pt_to_geom(TEXT_PT$protein),
  point_size = 1.38
) {
  if (nrow(labels) == 0) return(p)

  p +
    ggplot2::geom_point(
      data = labels,
      ggplot2::aes(x = x_point, y = y_point, fill = class),
      inherit.aes = FALSE,
      shape = 21,
      size = point_size,
      stroke = 0.28,
      colour = COL_DARK
    ) +
    ggplot2::geom_segment(
      data = labels,
      ggplot2::aes(
        x = x_point, y = y_point,
        xend = label_x, yend = label_y
      ),
      inherit.aes = FALSE,
      linewidth = 0.20,
      colour = COL_MID
    ) +
    # White underprint creates a clean halo without a visible rectangle.
    ggplot2::geom_text(
      data = labels,
      ggplot2::aes(
        x = label_x, y = label_y,
        label = gene, hjust = hjust
      ),
      inherit.aes = FALSE,
      family = FONT_ARIAL,
      size = label_size + pt_to_geom(1.6),
      fontface = "plain",
      colour = "white"
    ) +
    ggplot2::geom_text(
      data = labels,
      ggplot2::aes(
        x = label_x, y = label_y,
        label = gene, hjust = hjust
      ),
      inherit.aes = FALSE,
      family = FONT_ARIAL,
      size = label_size,
      fontface = "plain",
      colour = COL_DARK
    )
}

###############################################################################
# 6. Panel a — primary volcano
###############################################################################

volcano_data <- fig2a |>
  dplyr::transmute(
    gene = as_chr(EntrezGeneSymbol),
    effect = as_num(logFC),
    fdr = as_num(adj.P.Val),
    y = neglog10(adj.P.Val),
    class = dplyr::case_when(
      as_num(adj.P.Val) < 0.05 & as_num(logFC) > 0 ~ "Higher in AD",
      as_num(adj.P.Val) < 0.05 & as_num(logFC) < 0 ~ "Lower in AD",
      TRUE ~ "Not significant"
    )
  ) |>
  dplyr::filter(is.finite(effect), is.finite(y))

volcano_offsets <- tibble::tribble(
  ~gene,    ~dx,    ~dy,
  "C3",     -0.060,  0.030,
  "PLCH1",  -0.095,  0.055,
  "PUM2",   -0.035,  0.115,
  "GDI1",    0.038, -0.085,
  "AMBRA1", -0.075,  0.075,
  "SMOC1",   0.040,  0.060,
  "CPLX2",   0.050,  0.040,
  "GJB3",    0.040,  0.075,
  "LRRN1",   0.030,  0.050,
  "SPC25",   0.030,  0.020
)

volcano_labels <- make_manual_labels(
  volcano_data, "gene", "effect", "y", volcano_offsets
) |>
  dplyr::filter(gene %in% FOCAL_PROTEINS)

missing_focal_volcano <- setdiff(FOCAL_PROTEINS, volcano_data$gene)
if (length(missing_focal_volcano) > 0) {
  stop(
    "The following focal proteins are absent from the volcano source data: ",
    paste(missing_focal_volcano, collapse = ", "),
    call. = FALSE
  )
}

n_high <- sum(volcano_data$class == "Higher in AD")
n_low <- sum(volcano_data$class == "Lower in AD")

p2a <- ggplot2::ggplot(volcano_data, ggplot2::aes(effect, y)) +
  ggplot2::geom_point(
    data = dplyr::filter(volcano_data, class == "Not significant"),
    colour = COL_MID,
    size = 0.48,
    alpha = 0.42
  ) +
  ggplot2::geom_point(
    data = dplyr::filter(volcano_data, class != "Not significant"),
    ggplot2::aes(fill = class),
    shape = 21,
    size = 0.72,
    stroke = 0.18,
    colour = COL_DARK,
    alpha = 0.88
  ) +
  ggplot2::geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    linewidth = LINE_MM,
    colour = COL_MID
  ) +
  ggplot2::geom_vline(
    xintercept = 0,
    linewidth = LINE_MM,
    colour = COL_DARK
  ) +
  ggplot2::scale_fill_manual(
    values = c("Higher in AD" = COL_ORANGE, "Lower in AD" = COL_BLUE),
    guide = "none"
  ) +
  ggplot2::annotate(
    "text",
    x = min(volcano_data$effect, na.rm = TRUE) +
      0.12 * diff(range(volcano_data$effect, na.rm = TRUE)),
    y = max(volcano_data$y, na.rm = TRUE) * 0.91,
    label = paste0(n_low, " lower"),
    family = FONT_ARIAL,
    fontface = "bold",
    size = pt_to_geom(TEXT_PT$protein),
    hjust = 0,
    colour = COL_BLUE
  ) +
  ggplot2::annotate(
    "text",
    x = max(volcano_data$effect, na.rm = TRUE) -
      0.12 * diff(range(volcano_data$effect, na.rm = TRUE)),
    y = max(volcano_data$y, na.rm = TRUE) * 0.91,
    label = paste0(n_high, " higher"),
    family = FONT_ARIAL,
    fontface = "bold",
    size = pt_to_geom(TEXT_PT$protein),
    hjust = 1,
    colour = COL_ORANGE
  ) +
  ggplot2::labs(
    x = expression(Adjusted~log[2]~fold-change~("AD versus CN")),
    y = expression(-log[10]~("BH-FDR"))
  ) +
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.08, 0.11))) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.12))) +
  ggplot2::coord_cartesian(clip = "off") +
  nature_theme()

p2a <- add_manual_labels(p2a, volcano_labels)

###############################################################################
# 7. Panel b — displayed Reactome pathways
###############################################################################

path_data <- fig2b |>
  dplyr::transmute(
    label = stringr::str_wrap(as_chr(Display_label), width = 23),
    NES = as_num(NES),
    fdr = as_num(p.adjust),
    direction = dplyr::if_else(NES >= 0, "Higher in AD", "Lower in AD")
  ) |>
  dplyr::filter(is.finite(NES)) |>
  dplyr::mutate(label = forcats::fct_reorder(label, NES))

p2b <- ggplot2::ggplot(
  path_data,
  ggplot2::aes(NES, label, fill = direction)
) +
  ggplot2::geom_col(
    width = 0.62,
    colour = COL_DARK,
    linewidth = 0.22
  ) +
  ggplot2::geom_vline(
    xintercept = 0,
    linewidth = LINE_MM,
    colour = COL_DARK
  ) +
  ggplot2::scale_fill_manual(
    values = c("Higher in AD" = COL_ORANGE, "Lower in AD" = COL_BLUE)
  ) +
  ggplot2::labs(
    x = "Normalized enrichment score (NES)",
    y = NULL,
    fill = NULL
  ) +
  nature_theme() +
  ggplot2::theme(
    legend.position = "none",
    axis.text.y = ggplot2::element_text(
      size = TEXT_PT$protein,
      hjust = 1,
      lineheight = 0.90,
      margin = ggplot2::margin(r = 1.5, unit = "pt")
    ),
    plot.margin = ggplot2::margin(1.2, 1.8, 1.2, 1.2, "mm")
  )

###############################################################################
# 8. Panel c — ONE integrated annotation + correlation heatmap
###############################################################################

protein_order <- FOCAL_PROTEINS

missing_corr_proteins <- setdiff(
  protein_order,
  unique(as_chr(fig2c_corr_full$Protein))
)
if (length(missing_corr_proteins) > 0) {
  stop(
    "Focal proteins absent from the complete representative-protein correlation source: ",
    paste(missing_corr_proteins, collapse = ", "),
    call. = FALSE
  )
}

trait_order <- unique(as_chr(fig2c_corr_full$Trait))
trait_order <- trait_order[trait_order != "GFAP"]

fig2c_corr <- fig2c_corr_full |>
  dplyr::filter(
    as_chr(Protein) %in% protein_order,
    as_chr(Trait) %in% trait_order
  ) |>
  dplyr::mutate(
    Protein = as_chr(Protein),
    Trait = as_chr(Trait),
    rho = as_num(rho),
    P.Value = as_num(P.Value)
  )

# Recalculate BH-FDR across exactly the protein-trait pairs displayed in the
# definitive 10-protein heatmap.
fig2c_corr$adj.P.Val <- stats::p.adjust(fig2c_corr$P.Value, method = "BH")

expected_corr_rows <- length(protein_order) * length(trait_order)
if (nrow(fig2c_corr) != expected_corr_rows) {
  stop(
    "Incomplete focal heatmap correlation matrix: expected ",
    expected_corr_rows, " rows but found ", nrow(fig2c_corr), ".",
    call. = FALSE
  )
}

fig2c_ann <- primary_dep |>
  dplyr::filter(as_chr(EntrezGeneSymbol) %in% protein_order) |>
  dplyr::transmute(
    EntrezGeneSymbol = as_chr(EntrezGeneSymbol),
    AptName = as_chr(AptName),
    primary_logFC = as_num(logFC),
    primary_FDR = as_num(adj.P.Val)
  ) |>
  dplyr::left_join(
    apoe_equal |>
      dplyr::transmute(
        EntrezGeneSymbol = as_chr(EntrezGeneSymbol),
        APOE_preserved_FDR005 = as_bool(adjusted_fdr005)
      ),
    by = "EntrezGeneSymbol"
  ) |>
  dplyr::left_join(
    primary_vs_cdrsb |>
      dplyr::transmute(
        EntrezGeneSymbol = as_chr(EntrezGeneSymbol),
        CDRSB_same_direction = as_bool(same_direction)
      ),
    by = "EntrezGeneSymbol"
  )

missing_ann_proteins <- setdiff(
  protein_order,
  unique(as_chr(fig2c_ann$EntrezGeneSymbol))
)
if (length(missing_ann_proteins) > 0) {
  stop(
    "Could not rebuild focal annotations for: ",
    paste(missing_ann_proteins, collapse = ", "),
    call. = FALSE
  )
}

# Internal IDs keep annotation columns distinct from traits with the same
# displayed name (especially CDR-SB). This avoids duplicated factor levels.
annotation_ids <- c("ANN_DEP", "ANN_APOE", "ANN_CDRSB")
column_order_ids <- c(annotation_ids, paste0("TRAIT__", trait_order))
column_labels <- c(
  ANN_DEP = "DEP dir.",
  ANN_APOE = "APOE pres.",
  ANN_CDRSB = "CDR-SB conc.",
  stats::setNames(trait_order, paste0("TRAIT__", trait_order))
)

# Rich-text labels italicize only the gene symbol APOE without invoking
# plotmath. This avoids the Windows PDF CID-font error that can occur when
# parsed expressions are measured during cowplot legend extraction.
column_labels_rich <- unname(column_labels[column_order_ids])
column_labels_rich[column_order_ids == "ANN_APOE"] <- "<i>APOE</i> pres."
names(column_labels_rich) <- column_order_ids

trait_cells <- fig2c_corr |>
  dplyr::filter(
    as_chr(Trait) != "GFAP",
    as_chr(Protein) %in% protein_order
  ) |>
  dplyr::transmute(
    Protein = as_chr(Protein),
    Column_ID = paste0("TRAIT__", as_chr(Trait)),
    rho = as_num(rho),
    stars = dplyr::case_when(
      as_num(adj.P.Val) < 0.001 ~ "***",
      as_num(adj.P.Val) < 0.01 ~ "**",
      as_num(adj.P.Val) < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

annotation_cells <- fig2c_ann |>
  dplyr::filter(as_chr(EntrezGeneSymbol) %in% protein_order) |>
  dplyr::transmute(
    Protein = as_chr(EntrezGeneSymbol),
    ANN_DEP = dplyr::case_when(
      as_num(primary_logFC) > 0 ~ "Higher in AD",
      as_num(primary_logFC) < 0 ~ "Lower in AD",
      TRUE ~ "Neutral"
    ),
    ANN_APOE = dplyr::if_else(
      as_bool(APOE_preserved_FDR005), "Preserved", "Not preserved"
    ),
    ANN_CDRSB = dplyr::if_else(
      as_bool(CDRSB_same_direction), "Same direction", "Opposite direction"
    )
  ) |>
  tidyr::pivot_longer(
    cols = -Protein,
    names_to = "Column_ID",
    values_to = "annotation_state"
  ) |>
  dplyr::mutate(
    annotation_colour = dplyr::case_when(
      Column_ID == "ANN_DEP" & annotation_state == "Higher in AD" ~ COL_ORANGE,
      Column_ID == "ANN_DEP" & annotation_state == "Lower in AD" ~ COL_BLUE,
      Column_ID %in% c("ANN_APOE", "ANN_CDRSB") ~ "#F1F1F1",
      TRUE ~ COL_PALE
    ),
    annotation_symbol = dplyr::case_when(
      Column_ID == "ANN_APOE" & annotation_state == "Preserved" ~ "●",
      Column_ID == "ANN_APOE" & annotation_state == "Not preserved" ~ "○",
      Column_ID == "ANN_CDRSB" & annotation_state == "Same direction" ~ "●",
      Column_ID == "ANN_CDRSB" & annotation_state == "Opposite direction" ~ "×",
      TRUE ~ ""
    )
  )

trait_cells <- trait_cells |>
  dplyr::mutate(
    Protein = factor(Protein, levels = rev(protein_order)),
    Column_ID = factor(Column_ID, levels = column_order_ids)
  )

annotation_cells <- annotation_cells |>
  dplyr::mutate(
    Protein = factor(Protein, levels = rev(protein_order)),
    Column_ID = factor(Column_ID, levels = column_order_ids)
  )

if (anyDuplicated(column_order_ids)) {
  stop("Internal heatmap column IDs are duplicated.", call. = FALSE)
}

separator_x <- 3.5

p2c_base <- ggplot2::ggplot() +
  # Annotation tiles use fixed colours and therefore do not enter the rho legend.
  ggplot2::geom_tile(
    data = annotation_cells,
    ggplot2::aes(Column_ID, Protein),
    fill = annotation_cells$annotation_colour,
    colour = "white",
    linewidth = 0.17
  ) +
  ggplot2::geom_text(
    data = dplyr::filter(annotation_cells, annotation_symbol != ""),
    ggplot2::aes(Column_ID, Protein, label = annotation_symbol),
    family = FONT_ARIAL,
    size = pt_to_geom(TEXT_PT$myriad_heatmap),
    fontface = "bold",
    colour = COL_DARK
  ) +
  # Trait tiles alone define the continuous Spearman-rho scale.
  ggplot2::geom_tile(
    data = trait_cells,
    ggplot2::aes(Column_ID, Protein, fill = rho),
    colour = "white",
    linewidth = 0.17
  ) +
  ggplot2::geom_text(
    data = trait_cells,
    ggplot2::aes(Column_ID, Protein, label = stars),
    family = FONT_MYRIAD,
    size = pt_to_geom(TEXT_PT$myriad_heatmap),
    colour = COL_DARK
  ) +
  ggplot2::geom_vline(
    xintercept = separator_x,
    colour = COL_DARK,
    linewidth = 0.30
  ) +
  ggplot2::scale_x_discrete(
    limits = column_order_ids,
    labels = column_labels_rich,
    drop = FALSE
  ) +
  ggplot2::scale_fill_gradient2(
    low = COL_BLUE,
    mid = "white",
    high = COL_ORANGE,
    midpoint = 0,
    limits = c(-0.50, 0.50),
    oob = scales::squish,
    breaks = c(-0.50, -0.25, 0, 0.25, 0.50),
    labels = c("-0.50", "-0.25", "0", "0.25", "0.50"),
    name = "Spearman\nrho"
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_colourbar(
      order = 1,
      title.position = "top",
      title.hjust = 0.5,
      barheight = grid::unit(15, "mm"),
      barwidth = grid::unit(2.4, "mm")
    )
  ) +
  ggplot2::labs(
    x = NULL,
    y = NULL
  ) +
  nature_theme(TEXT_PT$axis_tick) +
  ggplot2::theme(
    axis.text.x = ggtext::element_markdown(
      family = FONT_ARIAL,
      angle = 50, hjust = 1, vjust = 1,
      size = TEXT_PT$axis_tick,
      face = "plain",
      colour = COL_DARK
    ),
    axis.text.y = ggplot2::element_text(
      family = FONT_ARIAL,
      size = TEXT_PT$protein,
      face = "plain",
      hjust = 1,
      margin = ggplot2::margin(r = 2.0, unit = "pt")
    ),
    axis.ticks = ggplot2::element_blank(),
    axis.line = ggplot2::element_blank(),
    legend.position = "right",
    legend.title = ggplot2::element_text(
      family = FONT_ARIAL, size = TEXT_PT$legend_title, face = "plain"
    ),
    legend.text = ggplot2::element_text(
      family = FONT_ARIAL, size = TEXT_PT$heatmap_colourbar_tick, face = "plain"
    ),
    plot.margin = ggplot2::margin(1.0, 1.2, 2.0, 3.5, "mm")
  )


# Extract the continuous Spearman-rho legend before suppressing it in the
# matrix. This creates a dedicated legend column and prevents any overlap with
# Aβ42/40 or the angled x-axis labels.
p2c_rho_legend <- cowplot::get_legend(
  p2c_base +
    ggplot2::theme(
      legend.position = "right",
      legend.margin = ggplot2::margin(0, 0, 0, 0, "mm"),
      legend.box.margin = ggplot2::margin(0, 0, 0, 0, "mm")
    )
)

if (is.null(p2c_rho_legend)) {
  stop("Could not extract the Spearman-rho legend for panel c.", call. = FALSE)
}

p2c_matrix <- p2c_base +
  ggplot2::theme(
    legend.position = "none",
    plot.margin = ggplot2::margin(0.5, 0.6, 1.1, 4.5, "mm")
  )

# A fixed annotation key occupies the lower part of the right-side legend
# column. The key is intentionally stacked because this column is narrow.
p2c_annotation_key <- cowplot::ggdraw() +
  cowplot::draw_label(
    "Annotations",
    x = 0.015, y = 0.98,
    hjust = 0, vjust = 1,
    fontfamily = FONT_ARIAL,
    fontface = "bold",
    size = TEXT_PT$legend_title,
    colour = COL_DARK
  ) +
  cowplot::draw_label(
    "APOE",
    x = 0.015, y = 0.79,
    hjust = 0, vjust = 1,
    fontfamily = FONT_ARIAL,
    fontface = "bold.italic",
    size = TEXT_PT$legend,
    colour = COL_DARK
  ) +
  cowplot::draw_label(
    "● preserved",
    x = 0.025, y = 0.61,
    hjust = 0, vjust = 1,
    fontfamily = FONT_ARIAL,
    fontface = "plain",
    size = TEXT_PT$legend,
    colour = COL_DARK
  ) +
  cowplot::draw_label(
    "○ not preserved",
    x = 0.025, y = 0.44,
    hjust = 0, vjust = 1,
    fontfamily = FONT_ARIAL,
    fontface = "plain",
    size = TEXT_PT$legend,
    colour = COL_DARK
  ) +
  cowplot::draw_label(
    "CDR-SB",
    x = 0.015, y = 0.27,
    hjust = 0, vjust = 1,
    fontfamily = FONT_ARIAL,
    fontface = "bold",
    size = TEXT_PT$legend,
    colour = COL_DARK
  ) +
  cowplot::draw_label(
    "● concordant",
    x = 0.025, y = 0.11,
    hjust = 0, vjust = 1,
    fontfamily = FONT_ARIAL,
    fontface = "plain",
    size = TEXT_PT$legend,
    colour = COL_DARK
  ) +
  cowplot::draw_label(
    "× opposite",
    x = 0.025, y = -0.05,
    hjust = 0, vjust = 1,
    fontfamily = FONT_ARIAL,
    fontface = "plain",
    size = TEXT_PT$legend,
    colour = COL_DARK
  )

p2c_legend_spacer <- cowplot::ggdraw()

p2c_legend_column <- cowplot::plot_grid(
  p2c_rho_legend,
  p2c_legend_spacer,
  p2c_annotation_key,
  ncol = 1,
  rel_heights = c(0.41, 0.07, 0.52)
)

# The matrix retains most of the width, while the legend column has a fixed
# editorial role. No legend text is drawn over the data region.
p2c <- cowplot::plot_grid(
  p2c_matrix,
  p2c_legend_column,
  nrow = 1,
  rel_widths = c(0.85, 0.15),
  align = "h",
  axis = "tb"
)

###############################################################################
# 9. Panel d — SAME-SAMPLE baseline versus AT(N)-adjusted
###############################################################################

atn_data <- fig2d |>
  dplyr::transmute(
    gene = as_chr(EntrezGeneSymbol),
    x = as_num(subset_baseline_logFC),
    y = as_num(adjusted_logFC),
    primary_effect = as_num(primary_full_logFC),
    primary_sig = as_bool(primary_full_fdr005),
    baseline_fdr = as_num(subset_baseline_adj.P.Val),
    adjusted_fdr = as_num(adjusted_adj.P.Val),
    class = dplyr::case_when(
      primary_sig & primary_effect > 0 ~ "Higher in AD",
      primary_sig & primary_effect < 0 ~ "Lower in AD",
      TRUE ~ "Not significant"
    )
  ) |>
  dplyr::filter(is.finite(x), is.finite(y))

atn_r <- stats::cor(atn_data$x, atn_data$y, method = "pearson")
atn_direction <- mean(sign(atn_data$x) == sign(atn_data$y))
atn_secondary_n <- sum(atn_data$adjusted_fdr < 0.05, na.rm = TRUE)

atn_offsets <- tibble::tribble(
  ~gene,   ~dx,    ~dy,
  "SPC25",  0.040,  0.060,
  "LRRN1",  0.050,  0.045,
  "GJB3",   0.055, -0.045,
  "PUM2",  -0.060,  0.060,
  "PLCH1", -0.110,  0.040,
  "C3",     0.030, -0.035,
  "GDI1",   0.035, -0.135
)

atn_labels <- make_manual_labels(
  atn_data, "gene", "x", "y", atn_offsets
)

atn_lim <- max(abs(c(atn_data$x, atn_data$y)), na.rm = TRUE) * 1.14

p2d <- ggplot2::ggplot(atn_data, ggplot2::aes(x, y)) +
  ggplot2::geom_point(
    data = dplyr::filter(atn_data, class == "Not significant"),
    colour = COL_MID,
    size = 0.50,
    alpha = 0.50
  ) +
  ggplot2::geom_point(
    data = dplyr::filter(atn_data, class != "Not significant"),
    ggplot2::aes(fill = class),
    shape = 21,
    size = 0.72,
    stroke = 0.17,
    colour = COL_DARK,
    alpha = 0.85
  ) +
  ggplot2::geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    linewidth = LINE_MM,
    colour = COL_MID
  ) +
  ggplot2::geom_hline(yintercept = 0, linewidth = LINE_MM, colour = COL_DARK) +
  ggplot2::geom_vline(xintercept = 0, linewidth = LINE_MM, colour = COL_DARK) +
  ggplot2::geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    linewidth = 0.40,
    colour = COL_DARK
  ) +
  ggplot2::scale_fill_manual(
    values = c("Higher in AD" = COL_ORANGE, "Lower in AD" = COL_BLUE),
    guide = "none"
  ) +
  ggplot2::annotate(
    "text",
    x = -atn_lim * 0.94,
    y = atn_lim * 0.92,
    hjust = 0,
    vjust = 1,
    label = paste0(
      "Equal-sample n = 343\n",
      "Adjusted FDR < 0.05: ", atn_secondary_n
    ),
    family = FONT_ARIAL,
    size = pt_to_geom(TEXT_PT$protein),
    colour = COL_DARK
  ) +
  ggplot2::annotate(
    "text",
    x = atn_lim * 0.94,
    y = -atn_lim * 0.92,
    hjust = 1,
    vjust = 0,
    label = sprintf(
      "r = %.2f\nDirection = %.1f%%",
      atn_r, 100 * atn_direction
    ),
    family = FONT_ARIAL,
    size = pt_to_geom(TEXT_PT$protein),
    colour = COL_DARK
  ) +
  ggplot2::coord_cartesian(
    xlim = c(-atn_lim, atn_lim),
    ylim = c(-atn_lim, atn_lim),
    clip = "off"
  ) +
  ggplot2::labs(
    x = "Same-sample baseline log2 fold-change",
    y = "AT(N)-adjusted log2 fold-change"
  ) +
  nature_theme() +
  ggplot2::theme(
    axis.title.x = ggplot2::element_text(
      margin = ggplot2::margin(t = 2.0, unit = "pt")
    ),
    plot.margin = ggplot2::margin(1.5, 2.5, 2.2, 2.5, "mm")
  )

p2d <- add_manual_labels(p2d, atn_labels)

###############################################################################
# 10. Panel e — LOCO summary
###############################################################################

loco_data <- fig2e |>
  dplyr::transmute(
    country = as_chr(excluded_country),
    r = as_num(logFC_correlation),
    highlight = stringr::str_detect(tolower(country), "colombia")
  ) |>
  dplyr::filter(is.finite(r)) |>
  dplyr::arrange(r) |>
  dplyr::mutate(country = factor(country, levels = country))

LOCO_REFERENCE_R <- 0.90

p2e <- ggplot2::ggplot(loco_data, ggplot2::aes(r, country)) +
  ggplot2::geom_vline(
    xintercept = LOCO_REFERENCE_R,
    linetype = "dashed",
    linewidth = 0.30,
    colour = COL_MID
  ) +
  ggplot2::geom_segment(
    ggplot2::aes(x = 0, xend = r, yend = country),
    colour = COL_LIGHT,
    linewidth = 0.45
  ) +
  ggplot2::geom_point(
    ggplot2::aes(fill = highlight),
    shape = 21,
    size = 1.55,
    stroke = 0.25,
    colour = COL_DARK
  ) +
  ggplot2::geom_text(
    ggplot2::aes(x = r + 0.018, label = sprintf("%.2f", r)),
    family = FONT_ARIAL,
    size = pt_to_geom(TEXT_PT$protein),
    hjust = 0,
    colour = COL_DARK
  ) +
  ggplot2::scale_fill_manual(
    values = c(`TRUE` = COL_ORANGE, `FALSE` = COL_BLUE),
    guide = "none"
  ) +
  ggplot2::scale_x_continuous(
    limits = c(0, 1.08),
    breaks = c(0, 0.3, 0.6, 0.9),
    expand = ggplot2::expansion(mult = c(0.01, 0.01))
  ) +
  ggplot2::labs(
    x = "Correlation with primary model",
    y = NULL
  ) +
  nature_theme() +
  ggplot2::theme(axis.text.y = ggplot2::element_text(
    family = FONT_ARIAL, size = TEXT_PT$protein, face = "plain"
  ))

###############################################################################
# 11. Panel f — primary versus mean LOCO
###############################################################################

mean_loco <- fig2f |>
  dplyr::transmute(
    gene = as_chr(EntrezGeneSymbol),
    x = as_num(main_logFC),
    y = as_num(mean_loco_logFC),
    fdr = as_num(main_adj.P.Val),
    class = dplyr::case_when(
      fdr < 0.05 & x > 0 ~ "Higher in AD",
      fdr < 0.05 & x < 0 ~ "Lower in AD",
      TRUE ~ "Not significant"
    )
  ) |>
  dplyr::filter(is.finite(x), is.finite(y))

loco_r <- stats::cor(mean_loco$x, mean_loco$y, method = "pearson")
loco_direction <- mean(sign(mean_loco$x) == sign(mean_loco$y))
loco_lim <- max(abs(c(mean_loco$x, mean_loco$y)), na.rm = TRUE) * 1.14

loco_offsets <- tibble::tribble(
  ~gene,   ~dx,    ~dy,
  "SPC25", -0.070, -0.010,
  "LRRN1",  0.060,  0.075,
  "CPLX2",  0.075, -0.100,
  "SMOC1", -0.045,  0.075,
  "GDI1",  -0.030, -0.060,
  "C3",     0.030,  0.040
)

loco_labels <- make_manual_labels(
  mean_loco, "gene", "x", "y", loco_offsets
)

p2f <- ggplot2::ggplot(mean_loco, ggplot2::aes(x, y)) +
  ggplot2::geom_point(
    data = dplyr::filter(mean_loco, class == "Not significant"),
    colour = COL_MID,
    size = 0.46,
    alpha = 0.50
  ) +
  ggplot2::geom_point(
    data = dplyr::filter(mean_loco, class != "Not significant"),
    ggplot2::aes(fill = class),
    shape = 21,
    size = 0.68,
    stroke = 0.16,
    colour = COL_DARK,
    alpha = 0.85
  ) +
  ggplot2::geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    linewidth = 0.35,
    colour = COL_MID
  ) +
  ggplot2::geom_hline(yintercept = 0, linewidth = LINE_MM, colour = COL_DARK) +
  ggplot2::geom_vline(xintercept = 0, linewidth = LINE_MM, colour = COL_DARK) +
  ggplot2::scale_fill_manual(
    values = c("Higher in AD" = COL_ORANGE, "Lower in AD" = COL_BLUE),
    guide = "none"
  ) +
  ggplot2::annotate(
    "text",
    x = loco_lim * 0.94,
    y = -loco_lim * 0.90,
    hjust = 1,
    vjust = 0,
    label = sprintf(
      "r = %.3f\nDirection = %.1f%%",
      loco_r, 100 * loco_direction
    ),
    family = FONT_ARIAL,
    size = pt_to_geom(TEXT_PT$protein),
    colour = COL_DARK
  ) +
  ggplot2::coord_cartesian(
    xlim = c(-loco_lim, loco_lim),
    ylim = c(-loco_lim, loco_lim),
    clip = "off"
  ) +
  ggplot2::labs(
    x = "Primary-model log2 fold-change",
    y = "Mean LOCO log2 fold-change"
  ) +
  nature_theme() +
  ggplot2::theme(
    axis.title.x = ggplot2::element_text(
      margin = ggplot2::margin(t = 2.0, unit = "pt")
    ),
    plot.margin = ggplot2::margin(1.5, 2.8, 2.2, 2.5, "mm")
  )

p2f <- add_manual_labels(p2f, loco_labels)

###############################################################################
# 12. Cowplot-aligned three-row editorial layout
###############################################################################

# First align the REAL plotting regions in each row. Titles are added only
# afterwards, so they cannot distort panel alignment.
aligned_row1 <- align_row_pair(p2a, p2b, "row 1 (a–b)")
aligned_row2 <- align_row_pair(p2c, p2d, "row 2 (c–d)")
aligned_row3 <- align_row_pair(p2e, p2f, "row 3 (e–f)")

direction_legend_grob <- make_direction_legend_grob()
m2a <- make_titled_panel(
  aligned_row1[[1]], "a", "Differential abundance"
)
m2b <- make_titled_panel(
  aligned_row1[[2]], "b", "Reactome enrichment",
  legend_grob = direction_legend_grob
)

m2c <- make_titled_panel(
  aligned_row2[[1]], "c", "Protein–trait associations",
  plot_height = 0.985
)
m2d <- make_titled_panel(
  aligned_row2[[2]], "d", "AT(N) sensitivity"
)

m2e <- make_titled_panel(
  aligned_row3[[1]], "e", "Country-exclusion stability"
)
m2f <- make_titled_panel(
  aligned_row3[[2]], "f", "Mean LOCO preservation"
)

# A single vertical division is used across all rows.
# The narrow spacer corresponds to approximately 2–3 mm at 183 mm width.
column_gap <- cowplot::ggdraw()

row1 <- cowplot::plot_grid(
  m2a, column_gap, m2b,
  nrow = 1,
  rel_widths = c(0.50, 0.012, 0.50),
  align = "h",
  axis = "tb"
)

row2 <- cowplot::plot_grid(
  m2c, column_gap, m2d,
  nrow = 1,
  rel_widths = c(0.50, 0.012, 0.50),
  align = "h",
  axis = "tb"
)

row3 <- cowplot::plot_grid(
  m2e, column_gap, m2f,
  nrow = 1,
  rel_widths = c(0.50, 0.012, 0.50),
  align = "h",
  axis = "tb"
)

# The bottom row remains modestly smaller, but no longer compressed as in v10.
figure2 <- cowplot::plot_grid(
  row1, row2, row3,
  ncol = 1,
  rel_heights = c(1, 1, 1),
  align = "v",
  axis = "lr"
)

###############################################################################
# 13. Export
###############################################################################

# Hard-stop audit: publication dimensions are fixed and must not drift when
# typography or panel content changes.
if (!identical(MAIN_WIDTH_MM, 180) || !identical(MAIN_HEIGHT_MM, 170)) {
  stop(
    "Publication artboard changed unexpectedly. Required: 180 mm × 170 mm; ",
    "observed: ", MAIN_WIDTH_MM, " mm × ", MAIN_HEIGHT_MM, " mm.",
    call. = FALSE
  )
}

message(
  "Fixed Nature Aging artboard confirmed: ",
  MAIN_WIDTH_MM, " mm × ", MAIN_HEIGHT_MM, " mm."
)
message(
  "Focal proteins in volcano and heatmap: ",
  paste(FOCAL_PROTEINS, collapse = ", ")
)

save_panel(m2a, "Fig2a_primary_volcano_ILLUSTRATOR_v29", 90, 58)
save_panel(m2b, "Fig2b_Reactome_ILLUSTRATOR_v29", 90, 58)
save_panel(m2c, "Fig2c_integrated_heatmap_ILLUSTRATOR_v29", 90, 58)
save_panel(m2d, "Fig2d_ATN_equal_sample_ILLUSTRATOR_v29", 90, 58)
save_panel(m2e, "Fig2e_LOCO_summary_ILLUSTRATOR_v29", 90, 58)
save_panel(m2f, "Fig2f_mean_LOCO_ILLUSTRATOR_v29", 90, 58)

pdf_path <- file.path(
  master_vector_root,
  "Figure2_DEP_NatureAging_ILLUSTRATOR_v29.pdf"
)
svg_path <- file.path(
  master_vector_root,
  "Figure2_DEP_NatureAging_ILLUSTRATOR_v29.svg"
)
png_path <- file.path(
  preview_root,
  "Figure2_DEP_NatureAging_ILLUSTRATOR_v29_preview.png"
)

# Open the SVG in Adobe Illustrator for the highest editability.
# The PDF is a vector backup and submission-quality reference.
safe_pdf(figure2, pdf_path, MAIN_WIDTH_MM, MAIN_HEIGHT_MM)
safe_svg(figure2, svg_path, MAIN_WIDTH_MM, MAIN_HEIGHT_MM)
safe_png(figure2, png_path, MAIN_WIDTH_MM, MAIN_HEIGHT_MM, dpi = 300L)

###############################################################################
# 14. Review metrics and temporary corrected panel-d source
###############################################################################

metrics <- tibble::tibble(
  panel = c("Fig2d_ATN_equal_sample", "Fig2f_mean_LOCO"),
  correlation = c(atn_r, loco_r),
  direction_consistency = c(atn_direction, loco_direction),
  significant_secondary_fdr005 = c(atn_secondary_n, NA_integer_)
)

readme_patch <- data.frame(
  Field = c(
    "Purpose",
    "Important",
    "Final integration"
  ),
  Content = c(
    "Temporary plot data for the corrected Figure 2d equal-sample AT(N) comparison.",
    "The x axis is the baseline model refitted in the same 343 participants used by the AT(N)-adjusted model.",
    "After visual approval, this sheet must replace the old Fig2d source within the definitive SourceData_Fig2_DEP workbook."
  )
)

source_patch_path <- file.path(
  out_root, "Figure2d_ATN_equal_sample_SOURCE_ILLUSTRATOR_v29.xlsx"
)

openxlsx::write.xlsx(
  list(
    README = readme_patch,
    Fig2d_ATN_equal_sample = as.data.frame(fig2d),
    Review_metrics = as.data.frame(metrics)
  ),
  file = source_patch_path,
  overwrite = TRUE
)

readr::write_csv(
  metrics,
  file.path(out_root, "Figure2_ILLUSTRATOR_v29_metrics.csv")
)


illustrator_guide_path <- file.path(
  illustrator_root,
  "README_OPEN_IN_ADOBE_ILLUSTRATOR.txt"
)

writeLines(
  c(
    "FIGURE 2 — ADOBE ILLUSTRATOR WORKFLOW",
    "",
    "1. Open master_vector/Figure2_DEP_NatureAging_ILLUSTRATOR_v29.svg",
    "   in Adobe Illustrator. The SVG is the preferred editable master.",
    "2. Confirm the artboard remains exactly 180 mm × 170 mm.",
    "3. Confirm Arial and Myriad Pro are available before opening, so text",
    "   remains editable and is not substituted.",
    "4. Save immediately as an Adobe Illustrator .ai working file.",
    "5. Use the individual files in panels_vector/ when panel-level editing",
    "   or replacement is easier than editing the full composite.",
    "6. The PDF is a vector backup and visual reference.",
    "7. Files in preview/ are raster previews only and should not be edited",
    "   or submitted as the vector master.",
    "",
    "Illustrator notes:",
    "- Some objects may open inside clipping groups. Ungroup only as needed.",
    "- Do not release every clipping mask globally; this can alter axes.",
    "- Keep text as live Arial text whenever possible.",
    "- Keep the orange #F46D43 and blue #3288BD unchanged.",
    "- Keep the 180 mm × 170 mm artboard fixed; do not enlarge the figure",
    "  to accommodate text. Reposition objects within the existing artboard.",
    "- Export the final submission copy from Illustrator as PDF, with fonts",
    "  embedded and no downsampling of vector elements.",
    "",
    "R cannot create a native .ai file directly. SVG and vector PDF are the",
    "appropriate editable interchange formats; save the opened SVG as .ai."
  ),
  con = illustrator_guide_path,
  useBytes = TRUE
)

vector_files <- c(
  pdf_path,
  svg_path,
  list.files(
    panel_root,
    pattern = "\\.(pdf|svg)$",
    full.names = TRUE,
    ignore.case = TRUE
  )
)

vector_manifest <- data.frame(
  file = basename(vector_files),
  type = ifelse(
    grepl("\\.svg$", vector_files, ignore.case = TRUE),
    "SVG editable vector",
    "PDF editable vector"
  ),
  full_path = normalizePath(
    vector_files,
    winslash = "/",
    mustWork = FALSE
  ),
  exists = file.exists(vector_files),
  size_bytes = ifelse(
    file.exists(vector_files),
    file.info(vector_files)$size,
    NA_real_
  ),
  stringsAsFactors = FALSE
)

vector_manifest_path <- file.path(
  illustrator_root,
  "Figure2_ILLUSTRATOR_v29_vector_manifest.csv"
)
readr::write_csv(vector_manifest, vector_manifest_path)

if (!all(vector_manifest$exists)) {
  stop(
    "Illustrator vector export is incomplete. Missing files: ",
    paste(vector_manifest$file[!vector_manifest$exists], collapse = ", "),
    call. = FALSE
  )
}

# Final audit is based on the files actually written to disk.
# The previous object-name audit used exists(..., inherits = FALSE) inside
# vapply(), which checked the temporary callback environment and falsely
# reported that every object was missing even after successful export.
required_final_files <- c(
  pdf_path, svg_path, png_path,
  source_patch_path, vector_manifest_path
)
missing_final_files <- required_final_files[!file.exists(required_final_files)]
if (length(missing_final_files) > 0) {
  stop(
    "Figure 2 was NOT completed. Missing final files: ",
    paste(basename(missing_final_files), collapse = ", "),
    call. = FALSE
  )
}

# Confirm that the assembled figure and six titled panels are valid plot/grob
# objects before reporting success. These direct checks avoid environment-scope
# ambiguity while still detecting an interrupted build.
plot_objects <- list(m2a, m2b, m2c, m2d, m2e, m2f, figure2)
if (!all(vapply(plot_objects, function(x) inherits(x, c("gg", "ggplot", "grob", "gTree")), logical(1)))) {
  stop("Figure 2 was NOT completed: one or more assembled plot objects are invalid.", call. = FALSE)
}

message("============================================================")
message("Figure 2 Illustrator-ready v28 build completed successfully.")
message("Illustrator SVG master: ", svg_path)
message("Vector PDF master:     ", pdf_path)
message("PNG preview:           ", png_path)
message("Editable vector panels: ", panel_root)
message("Illustrator guide:      ", illustrator_guide_path)
message("Vector manifest:        ", vector_manifest_path)
message("AT(N) equal-sample r: ", sprintf("%.3f", atn_r))
message("AT(N) equal-sample direction: ", sprintf("%.1f%%", 100 * atn_direction))
message("Mean LOCO r: ", sprintf("%.4f", loco_r))
message("Illustrator output directory: ", illustrator_root)
legend_note_path <- file.path(
  out_root, "Figure2_ILLUSTRATOR_v29_legend_note.txt"
)
writeLines(
  c(
    "Panel c contains a dedicated right-side legend column.",
    "APOE: a filled circle indicates preservation after APOE adjustment; an open circle indicates non-preservation.",
    "CDR-SB: a filled circle indicates concordant direction; a multiplication sign indicates opposite direction."
  ),
  con = legend_note_path,
  useBytes = TRUE
)


warning_log_path <- file.path(
  out_root, "Figure2_ILLUSTRATOR_v29_warnings.txt"
)

current_warnings <- warnings()

if (is.null(current_warnings)) {
  writeLines(
    "No warnings were recorded during this R session.",
    con = warning_log_path,
    useBytes = TRUE
  )
} else {
  capture.output(
    print(current_warnings),
    file = warning_log_path
  )
}

message("============================================================")


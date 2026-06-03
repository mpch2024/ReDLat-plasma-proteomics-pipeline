###############################################################################
# ReDLat plasma proteomics — supplementary and robustness analyses
#
# This script generates supplementary analyses that support the primary
# differential plasma proteomics workflow. It is designed to be run after
# 01_data_processing_and_differential_analysis.R and, when needed, alongside
# 02_figure_generation_main_manuscript_Figure2_NatureAging_typography_v3.R.
#
# Scope:
# - Compare primary FDR thresholds (FDR < 0.05 and FDR < 0.01).
# - Summarize overlap and preservation across sensitivity models, including AT(N)-adjusted DEP.
# - Run supplementary corrected and nominal enrichment analyses.
# - Export disease, tissue, and cell-type enrichR analyses when available.
# - Generate extended cross-country summaries, intersections, Jaccard indices,
#   Fisher/Stouffer meta-analytic summaries, and country-level exports.
# - Provide interpretation-oriented robustness classification tables.
#
# Interpretation:
# - FDR < 0.05 is the primary discovery threshold used in the main analysis.
# - FDR < 0.01 is used as a stricter high-confidence subset.
# - Nominal P < 0.05 enrichment is exploratory and hypothesis-generating only.
# - LOCO/LOSO/country analyses are internal stability analyses, not external
#   replication.
###############################################################################

###############################################################################
# 00_setup
###############################################################################

project_root <- "C:/Users/mnpiz/Desktop/DEPs_Proteomic_Publishable_V2"
outdir <- project_root
setwd(outdir)

MAIN_FDR <- 0.05
STRICT_FDR <- 0.01
NOMINAL_P <- 0.05
MIN_GENES_ORA <- 5
MIN_SIGNAL_COUNTRY_DEPS <- 20
K_PERCENT_COMMON <- 0.67
SEED_GLOBAL <- 1234
set.seed(SEED_GLOBAL)

required_cran <- c(
  "dplyr", "tidyr", "purrr", "readr", "tibble", "stringr", "forcats",
  "rlang", "ggplot2", "ggrepel", "openxlsx", "ComplexUpset", "eulerr", "scales"
)

required_bioc <- c(
  "clusterProfiler", "org.Hs.eg.db", "ReactomePA", "BiocParallel"
)

optional_cran <- c("enrichR")

missing_cran <- required_cran[!vapply(required_cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran) > 0) {
  stop("Missing required CRAN packages: ", paste(missing_cran, collapse = ", "),
       "\nInstall them before running Script 03.")
}

missing_bioc <- required_bioc[!vapply(required_bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc) > 0) {
  stop("Missing required Bioconductor packages: ", paste(missing_bioc, collapse = ", "),
       "\nInstall them with BiocManager before running Script 03.")
}

invisible(lapply(c(required_cran, required_bioc), library, character.only = TRUE))
options(stringsAsFactors = FALSE)

###############################################################################
# 00_helpers
###############################################################################

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

first_existing_file <- function(paths) {
  hit <- paths[file.exists(paths)][1]
  if (length(hit) == 0 || is.na(hit)) return(NA_character_)
  hit
}

safe_file_tag <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

safe_write_csv <- function(x, file) {
  ensure_dir(dirname(file))
  readr::write_csv(x, file)
}

safe_read_csv <- function(file) {
  if (is.na(file) || !file.exists(file)) return(NULL)
  suppressMessages(readr::read_csv(file, show_col_types = FALSE))
}

save_xlsx_safe <- function(df, path, sheet = "Sheet1") {
  ensure_dir(dirname(path))
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheet)
  openxlsx::writeData(wb, sheet, df)
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
}

safe_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))

clip_p <- function(p) {
  p <- safe_numeric(p)
  pmin(pmax(p, 1e-300), 1 - 1e-16)
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

add_inventory <- function(inventory, section, output, note = NA_character_) {
  dplyr::bind_rows(
    inventory,
    tibble::tibble(
      section = as.character(section),
      output = as.character(output),
      note = as.character(note)
    )
  )
}

###############################################################################
# 00_visual_system_consistent_with_script02
###############################################################################

# Visual system copied from the canonical Script 02 v3 so that supplementary
# robustness figures use the same typography, colors, margins, and export logic
# as the main manuscript figure panels.
COL_UP        <- "#f46d43"
COL_DOWN      <- "#4682b4"
COL_NEUTRAL   <- "grey70"
COL_NS        <- "grey86"
COL_TEXT      <- "black"
COL_BG        <- "white"
COL_BORDER    <- "black"
COL_AD        <- COL_UP
COL_CN        <- COL_DOWN
COL_MISSING   <- "grey94"

FONT_FAMILY <- "Arial"
PANEL_TITLE_SIZE <- 10
TEXT_SIZE <- 8
TAG_SIZE <- 10

SUPP_PANEL_W_CM <- 18
SUPP_PANEL_H_CM <- 12
SUPP_SMALL_W_CM <- 18
SUPP_SMALL_H_CM <- 10

save_gg_pdf_png_cm <- function(plot_obj,
                               file_base,
                               width_cm = SUPP_PANEL_W_CM,
                               height_cm = SUPP_PANEL_H_CM,
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

standardize_dep_tbl <- function(dep_tbl) {
  dep_tbl <- rename_with_candidates(
    dep_tbl,
    list(
      Protein_Name     = c("Protein_Name", "ProteinLabel", "Target", "TargetFullName"),
      EntrezGeneSymbol = c("EntrezGeneSymbol", "Gene", "Symbol"),
      EntrezGeneID     = c("EntrezGeneID", "ENTREZID", "entrez_gene"),
      AptName          = c("AptName", "feature_id_raw", "feature_id"),
      SeqId            = c("SeqId", "seqid"),
      logFC            = c("logFC", "main_logFC", "estimate"),
      P.Value          = c("P.Value", "p_value", "pvalue", "P"),
      adj.P.Val        = c("adj.P.Val", "FDR", "qvalue", "q_value", "main_adj.P.Val"),
      type             = c("type", "Direction", "direction")
    )
  )
  
  if (!"Protein_Name" %in% names(dep_tbl)) dep_tbl$Protein_Name <- NA_character_
  if (!"EntrezGeneSymbol" %in% names(dep_tbl)) dep_tbl$EntrezGeneSymbol <- NA_character_
  if (!"EntrezGeneID" %in% names(dep_tbl)) dep_tbl$EntrezGeneID <- NA_real_
  if (!"SeqId" %in% names(dep_tbl)) dep_tbl$SeqId <- NA_character_
  if (!"P.Value" %in% names(dep_tbl)) dep_tbl$P.Value <- NA_real_
  if (!"type" %in% names(dep_tbl)) dep_tbl$type <- NA_character_
  
  required <- c("AptName", "logFC", "adj.P.Val")
  missing_required <- setdiff(required, names(dep_tbl))
  if (length(missing_required) > 0) {
    stop("DEP table is missing required columns after standardization: ",
         paste(missing_required, collapse = ", "))
  }
  
  dep_tbl %>%
    dplyr::mutate(
      logFC = safe_numeric(logFC),
      P.Value = safe_numeric(P.Value),
      adj.P.Val = safe_numeric(adj.P.Val),
      EntrezGeneID = safe_numeric(EntrezGeneID),
      Protein_Label = dplyr::coalesce(Protein_Name, EntrezGeneSymbol, AptName),
      Direction_FDR005 = dplyr::case_when(
        adj.P.Val < MAIN_FDR & logFC > 0 ~ "Higher in AD",
        adj.P.Val < MAIN_FDR & logFC < 0 ~ "Lower in AD",
        TRUE ~ "Not significant"
      ),
      Direction_FDR001 = dplyr::case_when(
        adj.P.Val < STRICT_FDR & logFC > 0 ~ "Higher in AD",
        adj.P.Val < STRICT_FDR & logFC < 0 ~ "Lower in AD",
        TRUE ~ "Not significant"
      ),
      type = dplyr::case_when(
        adj.P.Val < MAIN_FDR & logFC > 0 ~ "Up",
        adj.P.Val < MAIN_FDR & logFC < 0 ~ "Down",
        TRUE ~ "NS"
      )
    )
}

collapse_dep_to_gene_local <- function(dep_tbl) {
  dep_tbl %>%
    dplyr::filter(!is.na(EntrezGeneSymbol), EntrezGeneSymbol != "") %>%
    dplyr::arrange(adj.P.Val, dplyr::desc(abs(logFC)), AptName) %>%
    dplyr::distinct(EntrezGeneSymbol, .keep_all = TRUE)
}

make_down_up_annot <- function(df, score_col = "logFC", id_col = "AptName") {
  if (is.null(df) || nrow(df) == 0) {
    return(tibble::tibble())
  }
  
  df2 <- df %>% dplyr::arrange(dplyr::desc(abs(.data[[score_col]])))
  
  up <- df2 %>%
    dplyr::filter(.data[[score_col]] > 0) %>%
    dplyr::mutate(rank = dplyr::row_number()) %>%
    dplyr::select(
      rank,
      Up_AptName = dplyr::all_of(id_col),
      Up_Protein = Protein_Label,
      Up_Gene = EntrezGeneSymbol,
      Up_logFC = dplyr::all_of(score_col),
      Up_adjP = adj.P.Val,
      Up_P = P.Value
    )
  
  down <- df2 %>%
    dplyr::filter(.data[[score_col]] < 0) %>%
    dplyr::mutate(rank = dplyr::row_number()) %>%
    dplyr::select(
      rank,
      Down_AptName = dplyr::all_of(id_col),
      Down_Protein = Protein_Label,
      Down_Gene = EntrezGeneSymbol,
      Down_logFC = dplyr::all_of(score_col),
      Down_adjP = adj.P.Val,
      Down_P = P.Value
    )
  
  n <- max(nrow(up), nrow(down), 1L)
  up2 <- up %>% tidyr::complete(rank = 1:n)
  down2 <- down %>% tidyr::complete(rank = 1:n)
  
  dplyr::full_join(down2, up2, by = "rank") %>%
    dplyr::arrange(rank) %>%
    dplyr::select(-rank)
}

pairwise_set_matrices <- function(sets) {
  cn <- names(sets)
  inter <- matrix(0, nrow = length(cn), ncol = length(cn), dimnames = list(cn, cn))
  jacc <- matrix(NA_real_, nrow = length(cn), ncol = length(cn), dimnames = list(cn, cn))
  
  for (i in seq_along(cn)) {
    for (j in seq_along(cn)) {
      a <- sets[[i]]
      b <- sets[[j]]
      inter[i, j] <- length(intersect(a, b))
      u <- length(union(a, b))
      jacc[i, j] <- ifelse(u == 0, NA_real_, inter[i, j] / u)
    }
  }
  
  list(
    intersection = as.data.frame(inter) %>% tibble::rownames_to_column("Country"),
    jaccard = as.data.frame(jacc) %>% tibble::rownames_to_column("Country")
  )
}

compute_K <- function(m, pct = K_PERCENT_COMMON) {
  if (m <= 1) return(1)
  max(2, ceiling(pct * m))
}

prepare_entrez_sets <- function(dep_tbl, threshold_type = c("fdr", "nominal"), threshold = 0.05) {
  threshold_type <- match.arg(threshold_type)
  p_col <- if (threshold_type == "fdr") "adj.P.Val" else "P.Value"
  
  dep_tbl <- dep_tbl %>%
    dplyr::mutate(
      EntrezGeneID = safe_numeric(EntrezGeneID),
      logFC = safe_numeric(logFC),
      P.Value = safe_numeric(P.Value),
      adj.P.Val = safe_numeric(adj.P.Val)
    ) %>%
    dplyr::filter(!is.na(EntrezGeneID), !is.na(logFC), is.finite(logFC), !is.na(.data[[p_col]])) %>%
    dplyr::group_by(EntrezGeneID) %>%
    dplyr::arrange(.data[[p_col]], dplyr::desc(abs(logFC)), .by_group = TRUE) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup()
  
  list(
    mapped_df = dep_tbl,
    universe_ids = unique(as.character(dep_tbl$EntrezGeneID)),
    up_ids = dep_tbl %>% dplyr::filter(.data[[p_col]] < threshold, logFC > 0) %>% dplyr::pull(EntrezGeneID) %>% as.character() %>% unique(),
    down_ids = dep_tbl %>% dplyr::filter(.data[[p_col]] < threshold, logFC < 0) %>% dplyr::pull(EntrezGeneID) %>% as.character() %>% unique(),
    p_col = p_col,
    threshold = threshold,
    threshold_type = threshold_type
  )
}

run_directional_ora <- function(dep_tbl,
                                out_prefix,
                                out_folder,
                                threshold_type = c("fdr", "nominal"),
                                threshold = 0.05,
                                min_genes = MIN_GENES_ORA) {
  threshold_type <- match.arg(threshold_type)
  ensure_dir(out_folder)
  
  sets <- prepare_entrez_sets(dep_tbl, threshold_type = threshold_type, threshold = threshold)
  
  trace_tbl <- tibble::tibble(
    out_prefix = out_prefix,
    threshold_type = threshold_type,
    threshold = threshold,
    n_unique_entrez_universe = length(sets$universe_ids),
    n_higher_in_AD = length(sets$up_ids),
    n_lower_in_AD = length(sets$down_ids),
    note = ifelse(threshold_type == "nominal", "Exploratory only; not primary inference", "BH-FDR corrected input set")
  )
  safe_write_csv(trace_tbl, file.path(out_folder, paste0(out_prefix, "_ORA_trace.csv")))
  
  run_one <- function(ids, direction) {
    if (length(ids) < min_genes) {
      skipped <- tibble::tibble(direction = direction, n_genes = length(ids), note = "Skipped: too few genes")
      safe_write_csv(skipped, file.path(out_folder, paste0(out_prefix, "_", direction, "_too_few_genes.csv")))
      return(data.frame())
    }
    
    go <- tryCatch(
      clusterProfiler::enrichGO(
        gene = ids,
        universe = sets$universe_ids,
        OrgDb = org.Hs.eg.db::org.Hs.eg.db,
        ont = "ALL",
        keyType = "ENTREZID",
        pAdjustMethod = "BH",
        readable = TRUE
      ),
      error = function(e) NULL
    )
    
    kegg <- tryCatch(
      clusterProfiler::enrichKEGG(
        gene = ids,
        universe = sets$universe_ids,
        organism = "hsa",
        pAdjustMethod = "BH"
      ),
      error = function(e) NULL
    )
    
    reactome <- tryCatch(
      ReactomePA::enrichPathway(
        gene = ids,
        universe = sets$universe_ids,
        organism = "human",
        pAdjustMethod = "BH",
        readable = TRUE
      ),
      error = function(e) NULL
    )
    
    go_df <- if (!is.null(go)) as.data.frame(go) else data.frame()
    kegg_df <- if (!is.null(kegg)) as.data.frame(kegg) else data.frame()
    reactome_df <- if (!is.null(reactome)) as.data.frame(reactome) else data.frame()
    
    if (nrow(go_df) > 0) go_df$database <- "GO"
    if (nrow(kegg_df) > 0) kegg_df$database <- "KEGG"
    if (nrow(reactome_df) > 0) reactome_df$database <- "Reactome"
    
    safe_write_csv(go_df, file.path(out_folder, paste0(out_prefix, "_", direction, "_GO.csv")))
    safe_write_csv(kegg_df, file.path(out_folder, paste0(out_prefix, "_", direction, "_KEGG.csv")))
    safe_write_csv(reactome_df, file.path(out_folder, paste0(out_prefix, "_", direction, "_Reactome.csv")))
    
    dplyr::bind_rows(go_df, kegg_df, reactome_df) %>%
      dplyr::mutate(direction = direction, threshold_type = threshold_type, threshold = threshold)
  }
  
  up_combined <- run_one(sets$up_ids, "higher_in_AD")
  down_combined <- run_one(sets$down_ids, "lower_in_AD")
  combined <- dplyr::bind_rows(up_combined, down_combined)
  
  if (nrow(combined) > 0) {
    safe_write_csv(combined, file.path(out_folder, paste0(out_prefix, "_directional_ORA_combined.csv")))
  }
  
  invisible(list(trace = trace_tbl, combined = combined))
}

plot_threshold_overlap <- function(dep_tbl, out_file_base) {
  ensure_dir(dirname(out_file_base))
  
  plot_df <- dep_tbl %>%
    dplyr::mutate(
      status = dplyr::case_when(
        adj.P.Val < STRICT_FDR & logFC > 0 ~ "FDR<0.01 higher",
        adj.P.Val < MAIN_FDR & logFC > 0 ~ "FDR<0.05 higher only",
        adj.P.Val < STRICT_FDR & logFC < 0 ~ "FDR<0.01 lower",
        adj.P.Val < MAIN_FDR & logFC < 0 ~ "FDR<0.05 lower only",
        TRUE ~ "Not significant"
      ),
      status = factor(status, levels = c(
        "FDR<0.01 higher", "FDR<0.05 higher only",
        "FDR<0.01 lower", "FDR<0.05 lower only", "Not significant"
      ))
    )
  
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = logFC, y = -log10(pmax(adj.P.Val, 1e-300)), color = status)) +
    ggplot2::geom_point(size = 1.15, alpha = 0.75, stroke = 0) +
    ggplot2::geom_hline(yintercept = -log10(MAIN_FDR), linetype = 2, color = "grey55", linewidth = 0.30) +
    ggplot2::geom_hline(yintercept = -log10(STRICT_FDR), linetype = 3, color = "grey35", linewidth = 0.30) +
    ggplot2::geom_vline(xintercept = 0, color = "grey70", linewidth = 0.25) +
    ggplot2::scale_color_manual(values = c(
      "FDR<0.01 higher" = COL_UP,
      "FDR<0.05 higher only" = scales::alpha(COL_UP, 0.55),
      "FDR<0.01 lower" = COL_DOWN,
      "FDR<0.05 lower only" = scales::alpha(COL_DOWN, 0.55),
      "Not significant" = COL_NS
    )) +
    ggplot2::labs(
      x = "Adjusted log2 fold-change",
      y = expression(-log[10](FDR)),
      color = NULL,
      title = "FDR threshold structure of the AD plasma proteome"
    ) +
    theme_panel_review() +
    ggplot2::theme(legend.position = "right")
  
  save_gg_pdf_png_cm(p, out_file_base, width_cm = SUPP_PANEL_W_CM, height_cm = SUPP_PANEL_H_CM, dpi = 600)
}

run_enrichr_directional <- function(dep_tbl, out_prefix, out_folder, fdr = MAIN_FDR) {
  ensure_dir(out_folder)
  if (!requireNamespace("enrichR", quietly = TRUE)) {
    safe_write_csv(tibble::tibble(status = "skipped", reason = "enrichR package not installed"),
                   file.path(out_folder, paste0(out_prefix, "_enrichR_skipped.csv")))
    return(invisible(NULL))
  }
  
  available <- tryCatch(enrichR::listEnrichrDbs(), error = function(e) NULL)
  if (is.null(available) || nrow(available) == 0) {
    safe_write_csv(tibble::tibble(status = "skipped", reason = "Could not retrieve enrichR database list"),
                   file.path(out_folder, paste0(out_prefix, "_enrichR_skipped.csv")))
    return(invisible(NULL))
  }
  
  db_names <- available$libraryName
  disease_db_candidates <- c("DisGeNET", "DisGeNET_2023", "Jensen_DISEASES", "Disease_Perturbations_from_GEO_up", "Disease_Perturbations_from_GEO_down")
  tissue_db_candidates <- c("ARCHS4_Tissues", "Jensen_TISSUES", "Human_Gene_Atlas")
  cell_db_candidates <- c("PanglaoDB_Augmented_2021", "CellMarker_Augmented_2021", "Azimuth_Cell_Types_2021")
  
  chosen_dbs <- unique(stats::na.omit(c(
    disease_db_candidates[disease_db_candidates %in% db_names][1],
    tissue_db_candidates[tissue_db_candidates %in% db_names][1],
    cell_db_candidates[cell_db_candidates %in% db_names][1]
  )))
  
  if (length(chosen_dbs) == 0) {
    safe_write_csv(tibble::tibble(status = "skipped", reason = "No compatible enrichR disease/tissue/cell database found"),
                   file.path(out_folder, paste0(out_prefix, "_enrichR_skipped.csv")))
    return(invisible(NULL))
  }
  
  chosen_tbl <- tibble::tibble(database = chosen_dbs)
  safe_write_csv(chosen_tbl, file.path(out_folder, paste0(out_prefix, "_enrichR_databases_used.csv")))
  
  run_one <- function(symbols, direction) {
    symbols <- unique(stats::na.omit(as.character(symbols)))
    symbols <- symbols[symbols != ""]
    if (length(symbols) < 5) {
      safe_write_csv(tibble::tibble(direction = direction, n_symbols = length(symbols), status = "too_few_symbols"),
                     file.path(out_folder, paste0(out_prefix, "_", direction, "_enrichR_too_few_symbols.csv")))
      return(NULL)
    }
    
    enr <- tryCatch(enrichR::enrichr(symbols, databases = chosen_dbs), error = function(e) NULL)
    if (is.null(enr) || length(enr) == 0) return(NULL)
    
    out_all <- list()
    for (db in names(enr)) {
      tbl <- tibble::as_tibble(enr[[db]])
      if (nrow(tbl) == 0) next
      tbl <- tbl %>% dplyr::mutate(database = db, direction = direction)
      out_all[[db]] <- tbl
      safe_write_csv(tbl, file.path(out_folder, paste0(out_prefix, "_", direction, "_", safe_file_tag(db), "_enrichR.csv")))
    }
    dplyr::bind_rows(out_all)
  }
  
  up_symbols <- dep_tbl %>%
    dplyr::filter(adj.P.Val < fdr, logFC > 0) %>%
    dplyr::pull(EntrezGeneSymbol)
  
  down_symbols <- dep_tbl %>%
    dplyr::filter(adj.P.Val < fdr, logFC < 0) %>%
    dplyr::pull(EntrezGeneSymbol)
  
  up_res <- run_one(up_symbols, "higher_in_AD")
  down_res <- run_one(down_symbols, "lower_in_AD")
  combined <- dplyr::bind_rows(up_res, down_res)
  
  if (nrow(combined) > 0) {
    safe_write_csv(combined, file.path(out_folder, paste0(out_prefix, "_enrichR_combined.csv")))
  }
  
  invisible(combined)
}

###############################################################################
# 01_load_workspace_and_tables
###############################################################################

message("01_load_workspace_and_tables")

supp_root <- file.path(outdir, "result", "supplementary")
ensure_dir(supp_root)
ensure_dir(file.path(supp_root, "tables"))
ensure_dir(file.path(supp_root, "figures"))
ensure_dir(file.path(supp_root, "logs"))

supplementary_inventory <- tibble::tibble(section = character(), output = character(), note = character())

workspace_file <- first_existing_file(c(
  file.path(outdir, "result", "workspace", "analysis_workspace.RData"),
  file.path(outdir, "result", "workspace", "proteomics_master_analysis_workspace.RData"),
  file.path(outdir, "result", "workspace", "proteomics_master_reanalysis_workspace.RData"),
  file.path(outdir, "result", "analysis_workspace.RData"),
  file.path(outdir, "proteomics_master_reanalysis_workspace.RData"),
  file.path(outdir, "proteomics_master_analysis_workspace.RData")
))

if (is.na(workspace_file)) {
  stop("Workspace not found. Run Script 01 first. Expected result/workspace/analysis_workspace.RData or compatible alias.")
}

message("Loading workspace: ", workspace_file)
load(workspace_file)

# Main DEP gene-collapsed
if (exists("DEP_gene") && is.data.frame(DEP_gene)) {
  DEP_main_gene <- DEP_gene
} else if (exists("DEP") && is.data.frame(DEP)) {
  DEP_main_gene <- collapse_dep_to_gene_local(DEP)
} else {
  dep_file <- first_existing_file(c(
    file.path(outdir, "result", "03_dep", "gene_collapsed", "AD_vs_CN_full_limma_results_gene_collapsed.csv"),
    file.path(outdir, "result", "dep", "AD_vs_CN_full_limma_results_gene_collapsed.csv"),
    file.path(outdir, "result", "dep", "AD_vs_CN_full_limma_results.csv")
  ))
  if (is.na(dep_file)) stop("Main DEP gene-collapsed table not found.")
  DEP_main_gene <- readr::read_csv(dep_file, show_col_types = FALSE)
}
DEP_main_gene <- standardize_dep_tbl(DEP_main_gene)

# Main DEP aptamer-level
DEP_main_aptamer <- NULL
if (exists("DEP_aptamer") && is.data.frame(DEP_aptamer)) {
  DEP_main_aptamer <- standardize_dep_tbl(DEP_aptamer)
} else if (exists("DEP") && is.data.frame(DEP)) {
  DEP_main_aptamer <- standardize_dep_tbl(DEP)
} else {
  apt_file <- first_existing_file(c(
    file.path(outdir, "result", "03_dep", "aptamer_level", "AD_vs_CN_full_limma_results_aptamer_level.csv"),
    file.path(outdir, "result", "dep", "AD_vs_CN_full_limma_results_aptamer_level.csv")
  ))
  if (!is.na(apt_file)) DEP_main_aptamer <- standardize_dep_tbl(readr::read_csv(apt_file, show_col_types = FALSE))
}

# Sensitivity DEP tables
load_sensitivity_dep <- function(object_names, file_candidates) {
  for (nm in object_names) {
    if (exists(nm, envir = .GlobalEnv) && is.data.frame(get(nm, envir = .GlobalEnv))) {
      return(standardize_dep_tbl(get(nm, envir = .GlobalEnv)))
    }
  }
  f <- first_existing_file(file_candidates)
  if (!is.na(f)) return(standardize_dep_tbl(readr::read_csv(f, show_col_types = FALSE)))
  NULL
}

DEP_APOE_gene_local <- load_sensitivity_dep(
  c("DEP_APOE_gene"),
  c(
    file.path(outdir, "result", "04_sensitivity", "apoe", "AD_vs_CN_APOE_adjusted_full_limma_results_gene_collapsed.csv"),
    file.path(outdir, "result", "dep_apoe", "AD_vs_CN_APOE_adjusted_full_limma_results.csv")
  )
)

DEP_CDRSB_gene_local <- load_sensitivity_dep(
  c("DEP_CDRSB_gene"),
  c(
    file.path(outdir, "result", "04_sensitivity", "cdrsb", "AD_only", "AD_only_CDRSB_severity_full_limma_results_gene_collapsed.csv"),
    file.path(outdir, "result", "04_sensitivity", "cdrsb", "AD_only_CDRSB_severity_full_limma_results_gene_collapsed.csv")
  )
)

DEP_CDRSB_ADJ_gene_local <- load_sensitivity_dep(
  c("DEP_CDRSB_ADJ_gene", "DEP_CDRSB_adjusted_gene", "DEP_CDRSB_ADJUSTED_gene"),
  c(
    file.path(outdir, "result", "04_sensitivity", "cdrsb", "AD_vs_CN_adjusted", "AD_vs_CN_CDRSB_adjusted_full_limma_results_gene_collapsed.csv"),
    file.path(outdir, "result", "04_sensitivity", "cdrsb", "AD_vs_CN_adjusted", "CDRSB_adjusted_AD_vs_CN_full_limma_results_gene_collapsed.csv")
  )
)

DEP_VASCULAR_gene_local <- load_sensitivity_dep(
  c("DEP_VASCULAR_gene", "DEP_METABOLIC_gene"),
  c(
    file.path(outdir, "result", "04_sensitivity", "vascular_metabolic", "AD_vs_CN_vascular_metabolic_adjusted_full_limma_results_gene_collapsed.csv")
  )
)

DEP_ATN_gene_local <- load_sensitivity_dep(
  c("DEP_ATN_gene", "DEP_ATN_adjusted_gene", "DEP_ATN_ADJUSTED_gene"),
  c(
    file.path(outdir, "result", "04_sensitivity", "atn_adjusted", "AD_vs_CN_ATN_adjusted_full_limma_results_gene_collapsed.csv")
  )
)

supplementary_inventory <- add_inventory(
  supplementary_inventory,
  "Workspace",
  workspace_file,
  "Workspace loaded successfully."
)

###############################################################################
# 02_FDR_thresholds_and_strict_outputs
###############################################################################

message("02_FDR_thresholds_and_strict_outputs")

strict_dir <- file.path(supp_root, "tables", "thresholds")
ensure_dir(strict_dir)

DEP_main_gene <- DEP_main_gene %>%
  dplyr::mutate(
    sig_fdr_005 = adj.P.Val < MAIN_FDR,
    sig_fdr_001 = adj.P.Val < STRICT_FDR,
    direction_fdr_005 = dplyr::case_when(
      sig_fdr_005 & logFC > 0 ~ "Higher in AD",
      sig_fdr_005 & logFC < 0 ~ "Lower in AD",
      TRUE ~ "Not significant"
    ),
    direction_fdr_001 = dplyr::case_when(
      sig_fdr_001 & logFC > 0 ~ "Higher in AD",
      sig_fdr_001 & logFC < 0 ~ "Lower in AD",
      TRUE ~ "Not significant"
    )
  )

safe_write_csv(DEP_main_gene, file.path(strict_dir, "main_DEP_gene_collapsed_with_FDR005_FDR001_flags.csv"))
safe_write_csv(DEP_main_gene %>% dplyr::filter(sig_fdr_001), file.path(strict_dir, "main_DEP_gene_collapsed_FDR001_strict.csv"))
safe_write_csv(DEP_main_gene %>% dplyr::filter(sig_fdr_005 & !sig_fdr_001), file.path(strict_dir, "main_DEP_gene_collapsed_FDR005_only_not_FDR001.csv"))

threshold_summary_gene <- tibble::tibble(
  universe = "gene_collapsed",
  threshold = c("FDR < 0.05", "FDR < 0.01"),
  total = c(sum(DEP_main_gene$adj.P.Val < MAIN_FDR, na.rm = TRUE), sum(DEP_main_gene$adj.P.Val < STRICT_FDR, na.rm = TRUE)),
  higher_in_AD = c(sum(DEP_main_gene$adj.P.Val < MAIN_FDR & DEP_main_gene$logFC > 0, na.rm = TRUE), sum(DEP_main_gene$adj.P.Val < STRICT_FDR & DEP_main_gene$logFC > 0, na.rm = TRUE)),
  lower_in_AD = c(sum(DEP_main_gene$adj.P.Val < MAIN_FDR & DEP_main_gene$logFC < 0, na.rm = TRUE), sum(DEP_main_gene$adj.P.Val < STRICT_FDR & DEP_main_gene$logFC < 0, na.rm = TRUE))
)

threshold_summary_aptamer <- NULL
if (!is.null(DEP_main_aptamer)) {
  threshold_summary_aptamer <- tibble::tibble(
    universe = "aptamer_level",
    threshold = c("FDR < 0.05", "FDR < 0.01"),
    total = c(sum(DEP_main_aptamer$adj.P.Val < MAIN_FDR, na.rm = TRUE), sum(DEP_main_aptamer$adj.P.Val < STRICT_FDR, na.rm = TRUE)),
    higher_in_AD = c(sum(DEP_main_aptamer$adj.P.Val < MAIN_FDR & DEP_main_aptamer$logFC > 0, na.rm = TRUE), sum(DEP_main_aptamer$adj.P.Val < STRICT_FDR & DEP_main_aptamer$logFC > 0, na.rm = TRUE)),
    lower_in_AD = c(sum(DEP_main_aptamer$adj.P.Val < MAIN_FDR & DEP_main_aptamer$logFC < 0, na.rm = TRUE), sum(DEP_main_aptamer$adj.P.Val < STRICT_FDR & DEP_main_aptamer$logFC < 0, na.rm = TRUE))
  )
}

threshold_summary <- dplyr::bind_rows(threshold_summary_gene, threshold_summary_aptamer)
safe_write_csv(threshold_summary, file.path(strict_dir, "FDR005_vs_FDR001_counts_summary.csv"))
plot_threshold_overlap(DEP_main_gene, file.path(supp_root, "figures", "main_DEP_FDR005_vs_FDR001_volcano"))

supplementary_inventory <- add_inventory(
  supplementary_inventory,
  "FDR thresholds",
  file.path(strict_dir, "FDR005_vs_FDR001_counts_summary.csv"),
  "Counts and tables comparing the main FDR 0.05 threshold with the stricter FDR 0.01 threshold."
)

###############################################################################
# 03_model_overlap_and_preservation_tables
###############################################################################

message("03_model_overlap_and_preservation_tables")

overlap_dir <- file.path(supp_root, "tables", "model_overlap")
ensure_dir(overlap_dir)

model_list <- list(
  main = DEP_main_gene,
  apoe = DEP_APOE_gene_local,
  atn_adjusted = DEP_ATN_gene_local,
  cdrsb_adjusted_ad_vs_cn = DEP_CDRSB_ADJ_gene_local,
  vascular_metabolic = DEP_VASCULAR_gene_local
)
model_list <- model_list[!vapply(model_list, is.null, logical(1))]

make_model_sig_tbl <- function(tbl, model_name, fdr = MAIN_FDR) {
  tbl %>%
    dplyr::transmute(
      EntrezGeneSymbol,
      AptName,
      Protein_Label,
      !!paste0(model_name, "_logFC") := logFC,
      !!paste0(model_name, "_P") := P.Value,
      !!paste0(model_name, "_adjP") := adj.P.Val,
      !!paste0(model_name, "_sig") := adj.P.Val < fdr,
      !!paste0(model_name, "_direction") := dplyr::case_when(
        adj.P.Val < fdr & logFC > 0 ~ "Higher in AD",
        adj.P.Val < fdr & logFC < 0 ~ "Lower in AD",
        TRUE ~ "Not significant"
      )
    )
}

overlap_tbl <- NULL
for (nm in names(model_list)) {
  tmp <- make_model_sig_tbl(model_list[[nm]], nm, fdr = MAIN_FDR)
  if (is.null(overlap_tbl)) overlap_tbl <- tmp else overlap_tbl <- overlap_tbl %>% dplyr::full_join(tmp, by = c("EntrezGeneSymbol", "AptName", "Protein_Label"))
}

if (!is.null(overlap_tbl)) {
  safe_write_csv(overlap_tbl, file.path(overlap_dir, "model_overlap_main_APOE_ATN_CDRSBadjusted_vascular_FDR005.csv"))
  save_xlsx_safe(overlap_tbl, file.path(overlap_dir, "model_overlap_main_APOE_ATN_CDRSBadjusted_vascular_FDR005.xlsx"))
  
  sig_cols <- grep("_sig$", names(overlap_tbl), value = TRUE)
  overlap_summary <- tibble::tibble(
    model = gsub("_sig$", "", sig_cols),
    n_sig_fdr005 = vapply(sig_cols, function(x) sum(overlap_tbl[[x]], na.rm = TRUE), numeric(1))
  )
  safe_write_csv(overlap_summary, file.path(overlap_dir, "model_overlap_counts_by_model_FDR005.csv"))
  
  if (all(c("main_sig", "apoe_sig") %in% names(overlap_tbl))) {
    main_apoe_preserved <- overlap_tbl %>%
      dplyr::filter(main_sig) %>%
      dplyr::mutate(
        apoe_preservation = dplyr::case_when(
          apoe_sig & sign(main_logFC) == sign(apoe_logFC) ~ "FDR-preserved same direction",
          !apoe_sig & sign(main_logFC) == sign(apoe_logFC) ~ "Direction-preserved only",
          !is.na(apoe_logFC) & sign(main_logFC) != sign(apoe_logFC) ~ "Direction changed",
          TRUE ~ "Missing"
        ),
        apoe_abs_ratio = abs(apoe_logFC) / pmax(abs(main_logFC), 1e-9)
      )
    safe_write_csv(main_apoe_preserved, file.path(overlap_dir, "main_DEP_APOE_preservation_FDR005.csv"))
  }
  
  if (all(c("main_sig", "cdrsb_adjusted_ad_vs_cn_sig") %in% names(overlap_tbl))) {
    main_cdrsb_adjusted_preserved <- overlap_tbl %>%
      dplyr::filter(main_sig) %>%
      dplyr::mutate(
        cdrsb_adjusted_ad_vs_cn_preservation = dplyr::case_when(
          cdrsb_adjusted_ad_vs_cn_sig & sign(main_logFC) == sign(cdrsb_adjusted_ad_vs_cn_logFC) ~ "FDR-preserved same direction",
          !cdrsb_adjusted_ad_vs_cn_sig & sign(main_logFC) == sign(cdrsb_adjusted_ad_vs_cn_logFC) ~ "Direction-preserved only",
          !is.na(cdrsb_adjusted_ad_vs_cn_logFC) & sign(main_logFC) != sign(cdrsb_adjusted_ad_vs_cn_logFC) ~ "Direction changed",
          TRUE ~ "Missing"
        ),
        cdrsb_adjusted_ad_vs_cn_abs_ratio = abs(cdrsb_adjusted_ad_vs_cn_logFC) / pmax(abs(main_logFC), 1e-9)
      )
    safe_write_csv(main_cdrsb_adjusted_preserved, file.path(overlap_dir, "main_DEP_CDRSB_adjusted_AD_vs_CN_preservation_FDR005.csv"))
  }
  
  if (all(c("main_sig", "atn_adjusted_sig") %in% names(overlap_tbl))) {
    main_atn_preserved <- overlap_tbl %>%
      dplyr::filter(main_sig) %>%
      dplyr::mutate(
        atn_adjusted_preservation = dplyr::case_when(
          atn_adjusted_sig & sign(main_logFC) == sign(atn_adjusted_logFC) ~ "FDR-preserved same direction",
          !atn_adjusted_sig & sign(main_logFC) == sign(atn_adjusted_logFC) ~ "Direction-preserved only",
          !is.na(atn_adjusted_logFC) & sign(main_logFC) != sign(atn_adjusted_logFC) ~ "Direction changed",
          TRUE ~ "Missing"
        ),
        atn_adjusted_abs_ratio = abs(atn_adjusted_logFC) / pmax(abs(main_logFC), 1e-9)
      )
    safe_write_csv(main_atn_preserved, file.path(overlap_dir, "main_DEP_ATN_adjusted_preservation_FDR005.csv"))
  }
  
  if (all(c("main_sig", "vascular_metabolic_sig") %in% names(overlap_tbl))) {
    main_vascular_preserved <- overlap_tbl %>%
      dplyr::filter(main_sig) %>%
      dplyr::mutate(
        vascular_metabolic_preservation = dplyr::case_when(
          vascular_metabolic_sig & sign(main_logFC) == sign(vascular_metabolic_logFC) ~ "FDR-preserved same direction",
          !vascular_metabolic_sig & sign(main_logFC) == sign(vascular_metabolic_logFC) ~ "Direction-preserved only",
          !is.na(vascular_metabolic_logFC) & sign(main_logFC) != sign(vascular_metabolic_logFC) ~ "Direction changed",
          TRUE ~ "Missing"
        ),
        vascular_metabolic_abs_ratio = abs(vascular_metabolic_logFC) / pmax(abs(main_logFC), 1e-9)
      )
    safe_write_csv(main_vascular_preserved, file.path(overlap_dir, "main_DEP_vascular_metabolic_preservation_FDR005.csv"))
  }
  
  supplementary_inventory <- add_inventory(
    supplementary_inventory,
    "Model overlap",
    file.path(overlap_dir, "model_overlap_main_APOE_ATN_CDRSBadjusted_vascular_FDR005.csv"),
    "Overlap and preservation of main DEP signals across APOE, AT(N)-adjusted, CDR-SB-adjusted AD-vs-CN attenuation, and vascular/metabolic sensitivity models. AD-only CDR-SB remains handled separately as within-AD severity alignment."
  )
}

###############################################################################
# 04_exploratory_nominal_and_corrected_ORA
###############################################################################

message("04_exploratory_nominal_and_corrected_ORA")

ora_root <- file.path(supp_root, "exploratory_nominal_enrichment")
nominal_lists_dir <- file.path(ora_root, "nominal_DEP_lists")
ensure_dir(nominal_lists_dir)

export_nominal_list <- function(tbl, label) {
  if (is.null(tbl)) return(NULL)
  out <- tbl %>%
    dplyr::filter(!is.na(P.Value), P.Value < NOMINAL_P) %>%
    dplyr::arrange(P.Value)
  file <- file.path(nominal_lists_dir, paste0(label, "_nominal_P005_gene_collapsed.csv"))
  safe_write_csv(out, file)
  tibble::tibble(
    model = label,
    n_nominal = nrow(out),
    n_nominal_higher_in_AD = sum(out$logFC > 0, na.rm = TRUE),
    n_nominal_lower_in_AD = sum(out$logFC < 0, na.rm = TRUE),
    file = file
  )
}

nominal_dep_summary <- dplyr::bind_rows(
  export_nominal_list(DEP_main_gene, "main"),
  export_nominal_list(DEP_APOE_gene_local, "APOE_adjusted"),
  export_nominal_list(DEP_CDRSB_gene_local, "AD_only_CDRSB_severity"),
  export_nominal_list(DEP_CDRSB_ADJ_gene_local, "CDRSB_adjusted_AD_vs_CN"),
  export_nominal_list(DEP_VASCULAR_gene_local, "vascular_metabolic_adjusted"),
  export_nominal_list(DEP_ATN_gene_local, "ATN_adjusted")
)
safe_write_csv(nominal_dep_summary, file.path(nominal_lists_dir, "nominal_DEP_lists_summary.csv"))

ora_runs <- list(
  list(tbl = DEP_main_gene, label = "main_FDR005", type = "fdr", threshold = MAIN_FDR),
  list(tbl = DEP_main_gene, label = "main_FDR001", type = "fdr", threshold = STRICT_FDR),
  list(tbl = DEP_main_gene, label = "main_nominal_P005", type = "nominal", threshold = NOMINAL_P),
  list(tbl = DEP_APOE_gene_local, label = "APOE_FDR005", type = "fdr", threshold = MAIN_FDR),
  list(tbl = DEP_APOE_gene_local, label = "APOE_FDR001", type = "fdr", threshold = STRICT_FDR),
  list(tbl = DEP_APOE_gene_local, label = "APOE_nominal_P005", type = "nominal", threshold = NOMINAL_P),
  list(tbl = DEP_ATN_gene_local, label = "ATN_adjusted_FDR005", type = "fdr", threshold = MAIN_FDR),
  list(tbl = DEP_ATN_gene_local, label = "ATN_adjusted_FDR001", type = "fdr", threshold = STRICT_FDR),
  list(tbl = DEP_ATN_gene_local, label = "ATN_adjusted_nominal_P005", type = "nominal", threshold = NOMINAL_P),
  list(tbl = DEP_CDRSB_gene_local, label = "AD_only_CDRSB_severity_FDR005", type = "fdr", threshold = MAIN_FDR),
  list(tbl = DEP_CDRSB_gene_local, label = "AD_only_CDRSB_severity_FDR001", type = "fdr", threshold = STRICT_FDR),
  list(tbl = DEP_CDRSB_gene_local, label = "AD_only_CDRSB_severity_nominal_P005", type = "nominal", threshold = NOMINAL_P),
  list(tbl = DEP_CDRSB_ADJ_gene_local, label = "CDRSB_adjusted_AD_vs_CN_FDR005", type = "fdr", threshold = MAIN_FDR),
  list(tbl = DEP_CDRSB_ADJ_gene_local, label = "CDRSB_adjusted_AD_vs_CN_FDR001", type = "fdr", threshold = STRICT_FDR),
  list(tbl = DEP_CDRSB_ADJ_gene_local, label = "CDRSB_adjusted_AD_vs_CN_nominal_P005", type = "nominal", threshold = NOMINAL_P),
  list(tbl = DEP_VASCULAR_gene_local, label = "vascular_metabolic_FDR005", type = "fdr", threshold = MAIN_FDR),
  list(tbl = DEP_VASCULAR_gene_local, label = "vascular_metabolic_nominal_P005", type = "nominal", threshold = NOMINAL_P)
)

ora_inventory <- list()
for (run in ora_runs) {
  if (is.null(run$tbl)) next
  out_folder <- file.path(ora_root, run$label)
  message("Running supplementary ORA: ", run$label)
  res <- tryCatch(
    run_directional_ora(
      dep_tbl = run$tbl,
      out_prefix = run$label,
      out_folder = out_folder,
      threshold_type = run$type,
      threshold = run$threshold
    ),
    error = function(e) {
      warning("ORA failed for ", run$label, ": ", e$message)
      NULL
    }
  )
  ora_inventory[[run$label]] <- tibble::tibble(
    label = run$label,
    threshold_type = run$type,
    threshold = run$threshold,
    output_folder = out_folder,
    status = ifelse(is.null(res), "failed", "completed")
  )
}

ora_inventory_tbl <- dplyr::bind_rows(ora_inventory)
safe_write_csv(ora_inventory_tbl, file.path(ora_root, "supplementary_ORA_inventory.csv"))

writeLines(c(
  "Supplementary ORA interpretation note:",
  "",
  "FDR-based ORA outputs use BH-corrected DEP sets and may support the main biological interpretation.",
  "Nominal P<0.05 ORA is exploratory only and should not be used as primary inference.",
  "Direction-specific sets are separated into higher-in-AD and lower-in-AD proteins."
), file.path(ora_root, "supplementary_ORA_interpretive_note.txt"))

supplementary_inventory <- add_inventory(
  supplementary_inventory,
  "Supplementary ORA",
  file.path(ora_root, "supplementary_ORA_inventory.csv"),
  "Includes corrected and nominal ORA. Nominal ORA is explicitly exploratory."
)

###############################################################################
# 05_enrichR_disease_tissue_celltype
###############################################################################

message("05_enrichR_disease_tissue_celltype")

enrichr_dir <- file.path(supp_root, "enrichr_disease_tissue_celltype")
ensure_dir(enrichr_dir)

enrichr_main <- tryCatch(
  run_enrichr_directional(DEP_main_gene, "main_DEP_FDR005", enrichr_dir, fdr = MAIN_FDR),
  error = function(e) {
    warning("enrichR main DEP failed: ", e$message)
    NULL
  }
)

enrichr_strict <- tryCatch(
  run_enrichr_directional(DEP_main_gene, "main_DEP_FDR001", enrichr_dir, fdr = STRICT_FDR),
  error = function(e) {
    warning("enrichR strict DEP failed: ", e$message)
    NULL
  }
)

supplementary_inventory <- add_inventory(
  supplementary_inventory,
  "enrichR disease/tissue/cell-type",
  enrichr_dir,
  "Restored from the original master script as supplementary/exploratory biological context. Requires internet access through enrichR."
)

###############################################################################
# 06_cross_country_extended_exports
###############################################################################

message("06_cross_country_extended_exports")

country_root <- file.path(supp_root, "cross_country_extended")
dir_qc <- file.path(country_root, "01_qc")
dir_by_country <- file.path(country_root, "02_by_country")
dir_overlap <- file.path(country_root, "03_overlap")
dir_meta <- file.path(country_root, "04_fisher_stouffer_meta")
dir_common <- file.path(country_root, "05_common_consistency")
ensure_dir(dir_qc); ensure_dir(dir_by_country); ensure_dir(dir_overlap); ensure_dir(dir_meta); ensure_dir(dir_common)

country_specific_tbl <- NULL
country_specific_file <- first_existing_file(c(
  file.path(outdir, "result", "06_robustness", "country_meta", "tables", "country_specific_DEP_results.csv"),
  file.path(outdir, "result", "06_robustness", "country_meta", "tables", "country_specific_effects_all.csv"),
  file.path(outdir, "result", "step3_sensitivity", "meta_analysis", "tables", "country_specific_effects_all.csv")
))

if (!is.na(country_specific_file)) {
  country_specific_tbl <- readr::read_csv(country_specific_file, show_col_types = FALSE)
} else if (exists("country_specific_tbl") && is.data.frame(country_specific_tbl)) {
  country_specific_tbl <- country_specific_tbl
} else if (exists("country_effects_tbl") && is.data.frame(country_effects_tbl)) {
  country_specific_tbl <- country_effects_tbl
}

if (!is.null(country_specific_tbl) && nrow(country_specific_tbl) > 0) {
  country_specific_tbl <- rename_with_candidates(
    country_specific_tbl,
    list(
      Country = c("Country", "country"),
      AptName = c("AptName", "feature_id_raw"),
      Protein_Name = c("Protein_Name", "ProteinLabel", "Target"),
      EntrezGeneSymbol = c("EntrezGeneSymbol", "Gene", "Symbol"),
      logFC = c("logFC", "country_logFC"),
      P.Value = c("P.Value", "p_value", "P"),
      adj.P.Val = c("adj.P.Val", "FDR", "qvalue"),
      se = c("se", "SE", "std_error")
    )
  ) %>%
    standardize_dep_tbl() %>%
    dplyr::mutate(
      Country = as.character(Country),
      P.Value = clip_p(P.Value),
      adj.P.Val = clip_p(adj.P.Val),
      is_sig_main = adj.P.Val < MAIN_FDR,
      is_sig_strict = adj.P.Val < STRICT_FDR,
      Direction_country = dplyr::case_when(
        logFC > 0 ~ "Higher in AD",
        logFC < 0 ~ "Lower in AD",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::distinct(Country, AptName, .keep_all = TRUE)
  
  countries <- sort(unique(country_specific_tbl$Country))
  n_countries <- length(countries)
  
  qc_counts <- country_specific_tbl %>%
    dplyr::group_by(Country) %>%
    dplyr::summarise(
      n_total = dplyr::n(),
      n_sig_fdr005 = sum(is_sig_main, na.rm = TRUE),
      n_sig_fdr001 = sum(is_sig_strict, na.rm = TRUE),
      n_higher_fdr005 = sum(is_sig_main & logFC > 0, na.rm = TRUE),
      n_lower_fdr005 = sum(is_sig_main & logFC < 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(n_sig_fdr005))
  
  safe_write_csv(qc_counts, file.path(dir_qc, "country_specific_DEP_QC_counts.csv"))
  save_xlsx_safe(qc_counts, file.path(dir_qc, "country_specific_DEP_QC_counts.xlsx"))
  
  signal_countries <- qc_counts %>%
    dplyr::filter(n_sig_fdr005 >= MIN_SIGNAL_COUNTRY_DEPS) %>%
    dplyr::pull(Country)
  
  safe_write_csv(tibble::tibble(signal_country = signal_countries), file.path(dir_qc, "signal_countries_FDR005.csv"))
  
  # By-country exports
  for (cty in countries) {
    cdir <- file.path(dir_by_country, safe_file_tag(cty))
    ensure_dir(cdir)
    
    df_ct <- country_specific_tbl %>%
      dplyr::filter(Country == cty) %>%
      dplyr::arrange(adj.P.Val, dplyr::desc(abs(logFC)))
    
    safe_write_csv(df_ct, file.path(cdir, "country_DEP_all.csv"))
    safe_write_csv(df_ct %>% dplyr::filter(is_sig_main), file.path(cdir, "country_DEP_FDR005.csv"))
    safe_write_csv(df_ct %>% dplyr::filter(is_sig_strict), file.path(cdir, "country_DEP_FDR001.csv"))
    save_xlsx_safe(make_down_up_annot(df_ct, "logFC", "AptName"), file.path(cdir, "country_DEP_all_DownUp.xlsx"))
    save_xlsx_safe(make_down_up_annot(df_ct %>% dplyr::filter(is_sig_main), "logFC", "AptName"), file.path(cdir, "country_DEP_FDR005_DownUp.xlsx"))
  }
  
  # Overlap sets and Jaccard matrices
  dep_sets_fdr005 <- country_specific_tbl %>%
    dplyr::filter(is_sig_main) %>%
    dplyr::group_by(Country) %>%
    dplyr::summarise(set = list(unique(AptName)), .groups = "drop") %>%
    tibble::deframe()
  
  dep_sets_fdr001 <- country_specific_tbl %>%
    dplyr::filter(is_sig_strict) %>%
    dplyr::group_by(Country) %>%
    dplyr::summarise(set = list(unique(AptName)), .groups = "drop") %>%
    tibble::deframe()
  
  for (ct in countries) {
    if (!ct %in% names(dep_sets_fdr005)) dep_sets_fdr005[[ct]] <- character(0)
    if (!ct %in% names(dep_sets_fdr001)) dep_sets_fdr001[[ct]] <- character(0)
  }
  dep_sets_fdr005 <- dep_sets_fdr005[countries]
  dep_sets_fdr001 <- dep_sets_fdr001[countries]
  
  mat_005 <- pairwise_set_matrices(dep_sets_fdr005)
  mat_001 <- pairwise_set_matrices(dep_sets_fdr001)
  
  safe_write_csv(mat_005$intersection, file.path(dir_overlap, "pairwise_intersections_FDR005.csv"))
  safe_write_csv(mat_005$jaccard, file.path(dir_overlap, "pairwise_jaccard_FDR005.csv"))
  safe_write_csv(mat_001$intersection, file.path(dir_overlap, "pairwise_intersections_FDR001.csv"))
  safe_write_csv(mat_001$jaccard, file.path(dir_overlap, "pairwise_jaccard_FDR001.csv"))
  
  global_intersection_005 <- if (length(dep_sets_fdr005) > 0) Reduce(intersect, dep_sets_fdr005) else character(0)
  global_intersection_001 <- if (length(dep_sets_fdr001) > 0) Reduce(intersect, dep_sets_fdr001) else character(0)
  
  safe_write_csv(tibble::tibble(AptName = global_intersection_005), file.path(dir_overlap, "global_intersection_all_countries_FDR005.csv"))
  safe_write_csv(tibble::tibble(AptName = global_intersection_001), file.path(dir_overlap, "global_intersection_all_countries_FDR001.csv"))
  
  # Fisher and signed Stouffer meta-analysis based on country-level P values
  meta_input <- country_specific_tbl %>%
    dplyr::transmute(
      Country,
      AptName,
      Protein_Label,
      EntrezGeneSymbol,
      P = clip_p(P.Value),
      logFC = logFC,
      Direction_country = Direction_country,
      is_sig_main = is_sig_main
    )
  
  fisher_tbl <- meta_input %>%
    dplyr::group_by(AptName) %>%
    dplyr::summarise(
      Protein_Label = dplyr::first(Protein_Label),
      EntrezGeneSymbol = dplyr::first(EntrezGeneSymbol),
      n_present = dplyr::n(),
      fisher_stat = -2 * sum(log(P), na.rm = TRUE),
      fisher_df = 2 * sum(!is.na(P)),
      p_fisher = ifelse(fisher_df > 0, stats::pchisq(fisher_stat, df = fisher_df, lower.tail = FALSE), NA_real_),
      .groups = "drop"
    ) %>%
    dplyr::mutate(p_fisher_adj = p.adjust(p_fisher, method = "BH"))
  
  stouffer_tbl <- meta_input %>%
    dplyr::mutate(
      z_abs = stats::qnorm(P / 2, lower.tail = FALSE),
      z_signed = z_abs * sign(logFC)
    ) %>%
    dplyr::group_by(AptName) %>%
    dplyr::summarise(
      n_present = dplyr::n(),
      n_z = sum(!is.na(z_signed)),
      z_stouffer = ifelse(n_z > 0, sum(z_signed, na.rm = TRUE) / sqrt(n_z), NA_real_),
      p_stouffer = ifelse(!is.na(z_stouffer), 2 * stats::pnorm(abs(z_stouffer), lower.tail = FALSE), NA_real_),
      direction_stouffer = dplyr::case_when(
        is.na(z_stouffer) ~ NA_character_,
        z_stouffer > 0 ~ "Higher in AD",
        z_stouffer < 0 ~ "Lower in AD",
        TRUE ~ NA_character_
      ),
      .groups = "drop"
    ) %>%
    dplyr::mutate(p_stouffer_adj = p.adjust(p_stouffer, method = "BH"))
  
  effect_summary <- meta_input %>%
    dplyr::group_by(AptName) %>%
    dplyr::summarise(
      median_logFC_all = stats::median(logFC, na.rm = TRUE),
      mean_logFC_all = mean(logFC, na.rm = TRUE),
      n_sig_fdr005 = sum(is_sig_main, na.rm = TRUE),
      n_higher_sig_fdr005 = sum(is_sig_main & logFC > 0, na.rm = TRUE),
      n_lower_sig_fdr005 = sum(is_sig_main & logFC < 0, na.rm = TRUE),
      sign_concordance_fraction = max(sum(logFC > 0, na.rm = TRUE), sum(logFC < 0, na.rm = TRUE)) / dplyr::n(),
      .groups = "drop"
    )
  
  fisher_stouffer_meta <- effect_summary %>%
    dplyr::left_join(fisher_tbl, by = "AptName") %>%
    dplyr::left_join(stouffer_tbl, by = "AptName") %>%
    dplyr::arrange(p_fisher_adj, p_stouffer_adj)
  
  safe_write_csv(fisher_stouffer_meta, file.path(dir_meta, "country_level_Fisher_Stouffer_meta_analysis.csv"))
  save_xlsx_safe(fisher_stouffer_meta, file.path(dir_meta, "country_level_Fisher_Stouffer_meta_analysis.xlsx"))
  
  # Common consistency sets
  K_global <- max(n_countries - 1, 2)
  common_global <- effect_summary %>%
    dplyr::mutate(
      common_direction = dplyr::case_when(
        n_higher_sig_fdr005 >= K_global ~ "Higher in AD",
        n_lower_sig_fdr005 >= K_global ~ "Lower in AD",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(common_direction)) %>%
    dplyr::left_join(DEP_main_gene %>% dplyr::select(AptName, Protein_Label, EntrezGeneSymbol, main_logFC = logFC, main_adjP = adj.P.Val), by = "AptName") %>%
    dplyr::arrange(dplyr::desc(abs(median_logFC_all)))
  
  safe_write_csv(common_global, file.path(dir_common, paste0("common_global_K", K_global, "_FDR005.csv")))
  save_xlsx_safe(common_global, file.path(dir_common, paste0("common_global_K", K_global, "_FDR005.xlsx")))
  
  if (length(signal_countries) >= 2) {
    K_signal <- compute_K(length(signal_countries))
    signal_data <- country_specific_tbl %>% dplyr::filter(Country %in% signal_countries)
    
    signal_summary <- signal_data %>%
      dplyr::group_by(AptName) %>%
      dplyr::summarise(
        n_signal_countries = dplyr::n_distinct(Country),
        n_higher_sig_fdr005 = sum(is_sig_main & logFC > 0, na.rm = TRUE),
        n_lower_sig_fdr005 = sum(is_sig_main & logFC < 0, na.rm = TRUE),
        median_logFC_signal = stats::median(logFC, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        common_direction_signal = dplyr::case_when(
          n_higher_sig_fdr005 >= K_signal ~ "Higher in AD",
          n_lower_sig_fdr005 >= K_signal ~ "Lower in AD",
          TRUE ~ NA_character_
        )
      ) %>%
      dplyr::filter(!is.na(common_direction_signal)) %>%
      dplyr::left_join(DEP_main_gene %>% dplyr::select(AptName, Protein_Label, EntrezGeneSymbol, main_logFC = logFC, main_adjP = adj.P.Val), by = "AptName") %>%
      dplyr::arrange(dplyr::desc(abs(median_logFC_signal)))
    
    safe_write_csv(signal_summary, file.path(dir_common, paste0("common_signal_countries_K", K_signal, "_FDR005.csv")))
    save_xlsx_safe(signal_summary, file.path(dir_common, paste0("common_signal_countries_K", K_signal, "_FDR005.xlsx")))
  }
  
  supplementary_inventory <- add_inventory(
    supplementary_inventory,
    "Cross-country extended",
    country_root,
    "Restores extended country exports, Jaccard/intersections, Fisher/Stouffer meta-analysis, and common-country consistency tables."
  )
  
} else {
  safe_write_csv(
    tibble::tibble(status = "not_available", note = "Country-specific DEP table was not found. Run country meta-analysis in Script 01 first."),
    file.path(country_root, "cross_country_extended_not_available.csv")
  )
  
  supplementary_inventory <- add_inventory(
    supplementary_inventory,
    "Cross-country extended",
    file.path(country_root, "cross_country_extended_not_available.csv"),
    "Skipped because country-specific DEP table was not available."
  )
}

###############################################################################
# 07_interpretive_robustness_classification
###############################################################################

message("07_interpretive_robustness_classification")

robust_dir <- file.path(supp_root, "robustness_classification")
ensure_dir(robust_dir)

main_vs_loco_mean_local <- NULL
if (exists("main_vs_loco_mean") && is.data.frame(main_vs_loco_mean)) {
  main_vs_loco_mean_local <- main_vs_loco_mean
} else {
  f <- first_existing_file(c(
    file.path(outdir, "result", "06_robustness", "country_loco", "tables", "main_vs_meanLOCO_table.csv"),
    file.path(outdir, "result", "step3_sensitivity", "country_loco", "tables", "main_vs_meanLOCO_table.csv")
  ))
  if (!is.na(f)) main_vs_loco_mean_local <- readr::read_csv(f, show_col_types = FALSE)
}

if (!is.null(main_vs_loco_mean_local)) {
  main_vs_loco_mean_local <- rename_with_candidates(
    main_vs_loco_mean_local,
    list(
      AptName = c("AptName"),
      main_logFC = c("main_logFC", "logFC_primary", "logFC_main"),
      mean_loco_logFC = c("mean_loco_logFC", "mean_LOCO_logFC"),
      main_adj.P.Val = c("main_adj.P.Val", "adj.P.Val", "main_FDR"),
      prop_same_direction = c("prop_same_direction", "prop_same_direction_loco")
    )
  ) %>%
    dplyr::mutate(
      main_logFC = safe_numeric(main_logFC),
      mean_loco_logFC = safe_numeric(mean_loco_logFC),
      loco_same_direction = sign(main_logFC) == sign(mean_loco_logFC),
      loco_abs_ratio = abs(mean_loco_logFC) / pmax(abs(main_logFC), 1e-9)
    )
}

meta_tbl_local <- NULL
if (exists("meta_tbl") && is.data.frame(meta_tbl)) {
  meta_tbl_local <- meta_tbl
} else {
  f <- first_existing_file(c(
    file.path(outdir, "result", "06_robustness", "country_meta", "tables", "country_meta_analysis_results.csv"),
    file.path(outdir, "result", "step3_sensitivity", "meta_analysis", "tables", "country_meta_analysis_results.csv")
  ))
  if (!is.na(f)) meta_tbl_local <- readr::read_csv(f, show_col_types = FALSE)
}

if (!is.null(meta_tbl_local)) {
  meta_tbl_local <- rename_with_candidates(
    meta_tbl_local,
    list(
      AptName = c("AptName"),
      meta_logFC = c("meta_logFC"),
      meta_adj.P.Val = c("meta_adj.P.Val", "meta_adjP", "meta_FDR"),
      I2 = c("I2")
    )
  ) %>%
    dplyr::mutate(
      meta_logFC = safe_numeric(meta_logFC),
      meta_adj.P.Val = safe_numeric(meta_adj.P.Val),
      I2 = safe_numeric(I2)
    )
}

balanced_tbl_local <- NULL
if (exists("balanced_protein_tbl") && is.data.frame(balanced_protein_tbl)) {
  balanced_tbl_local <- balanced_protein_tbl
} else {
  f <- first_existing_file(c(
    file.path(outdir, "result", "06_robustness", "balanced_country_resampling", "tables", "balanced_resampling_protein_stability.csv"),
    file.path(outdir, "result", "step4_balanced_country_resampling", "tables", "balanced_resampling_protein_stability.csv")
  ))
  if (!is.na(f)) balanced_tbl_local <- readr::read_csv(f, show_col_types = FALSE)
}

if (!is.null(balanced_tbl_local)) {
  balanced_tbl_local <- rename_with_candidates(
    balanced_tbl_local,
    list(
      AptName = c("AptName"),
      mean_bal_logFC = c("mean_bal_logFC", "balanced_logFC"),
      prop_same_direction = c("prop_same_direction", "prop_same_direction_balanced"),
      prop_sig_fdr005 = c("prop_sig_fdr005", "prop_sig_fdr005_balanced")
    )
  ) %>%
    dplyr::mutate(
      mean_bal_logFC = safe_numeric(mean_bal_logFC),
      prop_same_direction = safe_numeric(prop_same_direction),
      prop_sig_fdr005 = safe_numeric(prop_sig_fdr005)
    )
}

robust_tbl <- DEP_main_gene %>%
  dplyr::filter(adj.P.Val < MAIN_FDR) %>%
  dplyr::select(Protein_Label, Protein_Name, EntrezGeneSymbol, AptName, main_logFC = logFC, main_adjP = adj.P.Val, main_direction = Direction_FDR005)

if (!is.null(DEP_APOE_gene_local)) {
  robust_tbl <- robust_tbl %>%
    dplyr::left_join(DEP_APOE_gene_local %>% dplyr::select(AptName, apoe_logFC = logFC, apoe_adjP = adj.P.Val), by = "AptName") %>%
    dplyr::mutate(
      apoe_same_direction = sign(main_logFC) == sign(apoe_logFC),
      apoe_fdr_preserved = apoe_same_direction & apoe_adjP < MAIN_FDR,
      apoe_abs_ratio = abs(apoe_logFC) / pmax(abs(main_logFC), 1e-9)
    )
} else {
  robust_tbl <- robust_tbl %>%
    dplyr::mutate(apoe_logFC = NA_real_, apoe_adjP = NA_real_, apoe_same_direction = NA, apoe_fdr_preserved = NA, apoe_abs_ratio = NA_real_)
}

robust_tbl <- robust_tbl %>%
  dplyr::mutate(
    cdrsb_logFC = NA_real_,
    cdrsb_adjP = NA_real_,
    cdrsb_same_direction = NA,
    cdrsb_fdr_preserved = NA,
    cdrsb_abs_ratio = NA_real_,
    cdrsb_note = "CDR-SB is AD-only severity alignment and is not included in diagnostic robustness scoring."
  )

if (!is.null(DEP_VASCULAR_gene_local)) {
  robust_tbl <- robust_tbl %>%
    dplyr::left_join(DEP_VASCULAR_gene_local %>% dplyr::select(AptName, vascular_logFC = logFC, vascular_adjP = adj.P.Val), by = "AptName") %>%
    dplyr::mutate(
      vascular_same_direction = sign(main_logFC) == sign(vascular_logFC),
      vascular_fdr_preserved = vascular_same_direction & vascular_adjP < MAIN_FDR,
      vascular_abs_ratio = abs(vascular_logFC) / pmax(abs(main_logFC), 1e-9)
    )
} else {
  robust_tbl <- robust_tbl %>%
    dplyr::mutate(vascular_logFC = NA_real_, vascular_adjP = NA_real_, vascular_same_direction = NA, vascular_fdr_preserved = NA, vascular_abs_ratio = NA_real_)
}

if (!is.null(main_vs_loco_mean_local)) {
  robust_tbl <- robust_tbl %>%
    dplyr::left_join(main_vs_loco_mean_local %>% dplyr::select(AptName, mean_loco_logFC, loco_same_direction, loco_abs_ratio, dplyr::any_of("prop_same_direction")), by = "AptName")
} else {
  robust_tbl <- robust_tbl %>% dplyr::mutate(mean_loco_logFC = NA_real_, loco_same_direction = NA, loco_abs_ratio = NA_real_, prop_same_direction = NA_real_)
}

if (!is.null(meta_tbl_local)) {
  robust_tbl <- robust_tbl %>%
    dplyr::left_join(meta_tbl_local %>% dplyr::select(AptName, meta_logFC, meta_adj.P.Val, I2), by = "AptName") %>%
    dplyr::mutate(meta_same_direction = sign(main_logFC) == sign(meta_logFC), meta_fdr_preserved = meta_same_direction & meta_adj.P.Val < MAIN_FDR)
} else {
  robust_tbl <- robust_tbl %>% dplyr::mutate(meta_logFC = NA_real_, meta_adj.P.Val = NA_real_, I2 = NA_real_, meta_same_direction = NA, meta_fdr_preserved = NA)
}

if (!is.null(balanced_tbl_local)) {
  robust_tbl <- robust_tbl %>%
    dplyr::left_join(balanced_tbl_local %>% dplyr::select(AptName, mean_bal_logFC, prop_same_direction_balanced = prop_same_direction, prop_sig_fdr005_balanced = prop_sig_fdr005), by = "AptName") %>%
    dplyr::mutate(balanced_same_direction = sign(main_logFC) == sign(mean_bal_logFC))
} else {
  robust_tbl <- robust_tbl %>% dplyr::mutate(mean_bal_logFC = NA_real_, prop_same_direction_balanced = NA_real_, prop_sig_fdr005_balanced = NA_real_, balanced_same_direction = NA)
}

robust_tbl <- robust_tbl %>%
  dplyr::mutate(
    apoe_component = dplyr::case_when(apoe_fdr_preserved %in% TRUE ~ 2, apoe_same_direction %in% TRUE ~ 1, is.na(apoe_same_direction) ~ 0, TRUE ~ -1),
    cdrsb_component = dplyr::case_when(cdrsb_fdr_preserved %in% TRUE ~ 2, cdrsb_same_direction %in% TRUE ~ 1, is.na(cdrsb_same_direction) ~ 0, TRUE ~ -1),
    vascular_component = dplyr::case_when(vascular_fdr_preserved %in% TRUE ~ 2, vascular_same_direction %in% TRUE ~ 1, is.na(vascular_same_direction) ~ 0, TRUE ~ -1),
    loco_component = dplyr::case_when(loco_same_direction %in% TRUE & !is.na(loco_abs_ratio) & loco_abs_ratio >= 0.50 ~ 2, loco_same_direction %in% TRUE ~ 1, is.na(loco_same_direction) ~ 0, TRUE ~ -1),
    meta_component = dplyr::case_when(meta_fdr_preserved %in% TRUE ~ 2, meta_same_direction %in% TRUE ~ 1, is.na(meta_same_direction) ~ 0, TRUE ~ -1),
    balanced_component = dplyr::case_when(balanced_same_direction %in% TRUE & dplyr::coalesce(prop_same_direction_balanced, 0) >= 0.80 ~ 2, balanced_same_direction %in% TRUE ~ 1, is.na(balanced_same_direction) ~ 0, TRUE ~ -1),
    robustness_score = apoe_component + cdrsb_component + vascular_component + loco_component + meta_component + balanced_component,
    robustness_class = dplyr::case_when(
      robustness_score >= 9 ~ "High robustness across sensitivity layers",
      robustness_score >= 6 ~ "Moderate robustness",
      robustness_score >= 3 ~ "Partial/context-sensitive support",
      TRUE ~ "Limited or inconsistent support"
    )
  ) %>%
  dplyr::arrange(dplyr::desc(robustness_score), main_adjP)

safe_write_csv(robust_tbl, file.path(robust_dir, "main_DEP_FDR005_extended_robustness_classification.csv"))
save_xlsx_safe(robust_tbl, file.path(robust_dir, "main_DEP_FDR005_extended_robustness_classification.xlsx"))

robust_counts <- robust_tbl %>%
  dplyr::count(robustness_class, name = "n") %>%
  dplyr::arrange(dplyr::desc(n))
safe_write_csv(robust_counts, file.path(robust_dir, "main_DEP_FDR005_extended_robustness_classification_counts.csv"))

p_robust <- robust_counts %>%
  dplyr::mutate(robustness_class = forcats::fct_reorder(robustness_class, n)) %>%
  ggplot2::ggplot(ggplot2::aes(x = robustness_class, y = n)) +
  ggplot2::geom_col(fill = COL_NEUTRAL, color = COL_BORDER, linewidth = 0.25) +
  ggplot2::coord_flip() +
  ggplot2::labs(
    title = "Supplementary robustness classification of main FDR-significant proteins",
    x = NULL,
    y = "Number of proteins"
  ) +
  theme_panel_review()

save_gg_pdf_png_cm(
  p_robust,
  file.path(robust_dir, "main_DEP_FDR005_extended_robustness_classification_counts"),
  width_cm = SUPP_PANEL_W_CM,
  height_cm = SUPP_SMALL_H_CM,
  dpi = 600
)

writeLines(c(
  "Supplementary robustness classification note:",
  "",
  "This classification is intended for interpretation and reviewer support, not as a separate primary inferential layer.",
  "LOCO, country meta-analysis, and balanced resampling are internal stability analyses, not external replication.",
  "APOE, AT(N), vascular/metabolic and CDR-SB-adjusted AD-vs-CN models are diagnostic sensitivity/attenuation analyses. AD-only CDR-SB is analyzed separately as within-disease severity alignment and should not be interpreted as diagnostic attenuation or causal mediation.",
  "Proteins in stronger classes preserve direction and/or statistical support across available sensitivity layers."
), file.path(robust_dir, "robustness_classification_interpretive_note.txt"))

supplementary_inventory <- add_inventory(
  supplementary_inventory,
  "Extended robustness classification",
  file.path(robust_dir, "main_DEP_FDR005_extended_robustness_classification.csv"),
  "Classifies main FDR-significant proteins across APOE, vascular/metabolic, LOCO, meta-analysis, and balanced resampling layers when available. AD-only CDR-SB is excluded from diagnostic robustness scoring; CDR-SB-adjusted AD-vs-CN attenuation is exported in the model-overlap layer."
)

###############################################################################
# 08_secondary_clinical_severity_inventory
###############################################################################

message("08_secondary_clinical_severity_inventory")

clinical_dir <- file.path(supp_root, "tables", "clinical_severity_secondary")
ensure_dir(clinical_dir)

severity_combined_file <- first_existing_file(c(
  file.path(outdir, "result", "04_sensitivity", "clinical_severity", "severity_assoc_all_outcomes_combined.csv"),
  file.path(outdir, "result", "clinical_severity", "severity_assoc_all_outcomes_combined.csv")
))

if (!is.na(severity_combined_file)) {
  severity_combined <- readr::read_csv(severity_combined_file, show_col_types = FALSE)
  safe_write_csv(severity_combined, file.path(clinical_dir, "severity_assoc_all_outcomes_combined_copy.csv"))
  
  severity_summary <- severity_combined %>%
    dplyr::mutate(
      adj.P.Val = safe_numeric(adj.P.Val),
      adj.P.Val_within_outcome = if ("adj.P.Val_within_outcome" %in% names(.)) safe_numeric(adj.P.Val_within_outcome) else adj.P.Val
    ) %>%
    dplyr::group_by(outcome) %>%
    dplyr::summarise(
      n_tests = dplyr::n(),
      n_fdr005_global = sum(adj.P.Val < MAIN_FDR, na.rm = TRUE),
      n_fdr001_global = sum(adj.P.Val < STRICT_FDR, na.rm = TRUE),
      n_fdr005_within_outcome = sum(adj.P.Val_within_outcome < MAIN_FDR, na.rm = TRUE),
      n_fdr001_within_outcome = sum(adj.P.Val_within_outcome < STRICT_FDR, na.rm = TRUE),
      .groups = "drop"
    )
  
  safe_write_csv(severity_summary, file.path(clinical_dir, "clinical_severity_secondary_models_summary.csv"))
  
  supplementary_inventory <- add_inventory(
    supplementary_inventory,
    "Secondary clinical severity models",
    file.path(clinical_dir, "clinical_severity_secondary_models_summary.csv"),
    "Inventory of secondary protein-outcome associations from Script 01."
  )
} else {
  severity_files <- c(
    list.files(file.path(outdir, "result", "04_sensitivity", "clinical_severity"), pattern = "^severity_assoc_.*\\.csv$", full.names = TRUE),
    list.files(file.path(outdir, "result", "clinical_severity"), pattern = "^severity_assoc_.*\\.csv$", full.names = TRUE)
  )
  if (length(severity_files) > 0) {
    severity_inventory <- tibble::tibble(file = severity_files, outcome = gsub("^severity_assoc_|\\.csv$", "", basename(severity_files)))
    safe_write_csv(severity_inventory, file.path(clinical_dir, "clinical_severity_secondary_models_file_inventory.csv"))
  } else {
    safe_write_csv(tibble::tibble(status = "not_found", note = "No secondary clinical severity outputs were found."), file.path(clinical_dir, "clinical_severity_secondary_models_not_found.csv"))
  }
}

###############################################################################
# 09_protein_universe_summary_and_final_tables
###############################################################################

message("09_protein_universe_summary_and_final_tables")

summary_dir <- file.path(supp_root, "tables", "protein_universe_summary")
ensure_dir(summary_dir)

protein_summary <- tibble::tibble(
  metric = c(
    "seq_columns_detected",
    "internal_ADAT_protein_universe",
    "gene_collapsed_DEP_universe",
    "main_DEP_FDR005",
    "main_DEP_FDR001",
    "main_DEP_FDR005_higher_in_AD",
    "main_DEP_FDR005_lower_in_AD",
    "main_DEP_FDR001_higher_in_AD",
    "main_DEP_FDR001_lower_in_AD"
  ),
  value = c(
    if (exists("seq_cols")) length(seq_cols) else NA_integer_,
    if (exists("protein_universe")) length(protein_universe) else NA_integer_,
    nrow(DEP_main_gene),
    sum(DEP_main_gene$adj.P.Val < MAIN_FDR, na.rm = TRUE),
    sum(DEP_main_gene$adj.P.Val < STRICT_FDR, na.rm = TRUE),
    sum(DEP_main_gene$adj.P.Val < MAIN_FDR & DEP_main_gene$logFC > 0, na.rm = TRUE),
    sum(DEP_main_gene$adj.P.Val < MAIN_FDR & DEP_main_gene$logFC < 0, na.rm = TRUE),
    sum(DEP_main_gene$adj.P.Val < STRICT_FDR & DEP_main_gene$logFC > 0, na.rm = TRUE),
    sum(DEP_main_gene$adj.P.Val < STRICT_FDR & DEP_main_gene$logFC < 0, na.rm = TRUE)
  )
)

safe_write_csv(protein_summary, file.path(summary_dir, "protein_universe_and_DEP_summary.csv"))
safe_write_csv(DEP_main_gene, file.path(summary_dir, "main_DEP_gene_collapsed_complete_supplementary_table.csv"))
if (!is.null(DEP_main_aptamer)) safe_write_csv(DEP_main_aptamer, file.path(summary_dir, "main_DEP_aptamer_level_complete_supplementary_table.csv"))
save_xlsx_safe(DEP_main_gene, file.path(summary_dir, "main_DEP_gene_collapsed_complete_supplementary_table.xlsx"))

supplementary_inventory <- add_inventory(
  supplementary_inventory,
  "Protein universe summary",
  file.path(summary_dir, "protein_universe_and_DEP_summary.csv"),
  "Traceable counts for protein universe, gene-collapsed universe, and DEP thresholds."
)

###############################################################################
# 10_export_inventory_session_and_workspace
###############################################################################

message("10_export_inventory_session_and_workspace")

inventory_file <- file.path(supp_root, "supplementary_analysis_inventory.csv")
safe_write_csv(supplementary_inventory, inventory_file)

writeLines(capture.output(utils::sessionInfo()), file.path(supp_root, "logs", "script03_sessionInfo.txt"))

supp_workspace_file <- file.path(supp_root, "script03_supplementary_workspace.RData")
save(
  DEP_main_gene,
  DEP_main_aptamer,
  DEP_APOE_gene_local,
  DEP_CDRSB_gene_local,
  DEP_VASCULAR_gene_local,
  threshold_summary,
  overlap_tbl,
  nominal_dep_summary,
  ora_inventory_tbl,
  robust_tbl,
  robust_counts,
  protein_summary,
  supplementary_inventory,
  file = supp_workspace_file
)

cat("\nSCRIPT 03 SUPPLEMENTARY / EXPLORATORY 10/10 COMPLETE.\n")
cat("Supplementary inventory:\n", inventory_file, "\n")
cat("Supplementary workspace:\n", supp_workspace_file, "\n")
cat("Main supplementary folder:\n", supp_root, "\n")


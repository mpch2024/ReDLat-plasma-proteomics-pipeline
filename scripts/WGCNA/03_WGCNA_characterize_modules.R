###############################################################################
# ReDLat plasma proteomics — WGCNA workflow
# 03. Characterize modules and hubs
# Requires: outputs from Script 02
# Produces: kME, hubs, DEP burden and enrichment results
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

###############################################################################
# 1) PACKAGES
###############################################################################

cran_pkgs <- c(
  "dplyr",
  "tibble",
  "readr",
  "stringr",
  "ggplot2",
  "purrr",
  "tidyr",
  "forcats"
)

bioc_pkgs <- c(
  "WGCNA",
  "clusterProfiler",
  "ReactomePA",
  "org.Hs.eg.db",
  "enrichplot"
)

cran_missing <- cran_pkgs[
  !vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(cran_missing) > 0L) {
  stop("Missing required packages: ", paste(cran_missing, collapse = ", "),
       ". Run renv::restore() before this script.", call. = FALSE)
}

bioc_missing <- bioc_pkgs[
  !vapply(bioc_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(bioc_missing) > 0L) {
  stop("Missing required Bioconductor packages: ", paste(bioc_missing, collapse = ", "),
       ". Run renv::restore() before this script.", call. = FALSE)
}

invisible(lapply(
  c(cran_pkgs, bioc_pkgs),
  library,
  character.only = TRUE
))

options(stringsAsFactors = FALSE)
options(error = traceback)

###############################################################################
# 2) PATHS
###############################################################################

base_dir <- WGCNA_CONFIG$project_root

core_dir <- file.path(WGCNA_CONFIG$result_root,
  "02_network"
)

workspace_candidates <- c(
  file.path(
    core_dir,
    "workspace",
    "wgcna_core_light_workspace.RData"
  ),
  file.path(
    core_dir,
    "workspace",
    "wgcna_core_collapsed_workspace.RData"
  )
)

workspace_file <- workspace_candidates[
  file.exists(workspace_candidates)
][1]

gene_module_file <- file.path(
  core_dir,
  "tables",
  "gene_module_assignment.csv"
)

outdir <- file.path(WGCNA_CONFIG$result_root,
  "03_modules"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

subdirs <- c(
  "tables",
  "hubs",
  "enrichment",
  "enrichment/GO_BP",
  "enrichment/KEGG",
  "enrichment/Reactome",
  "enrichment/significant",
  "enrichment/top_terms",
  "figures",
  "figures/hubs",
  "figures/enrichment",
  "network_optional",
  "qc",
  "workspace"
)

invisible(lapply(
  file.path(outdir, subdirs),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

if (
  length(workspace_file) == 0 ||
  is.na(workspace_file) ||
  !file.exists(workspace_file)
) {
  stop(
    "No Script 11 workspace was found.\nAttempted:\n",
    paste(workspace_candidates, collapse = "\n"),
    "\n\nRun 02_WGCNA_construct_network.R first.",
    call. = FALSE
  )
}

if (!file.exists(gene_module_file)) {
  stop(
    "gene_module_assignment.csv from Script 11 was not found:\n",
    gene_module_file,
    call. = FALSE
  )
}

###############################################################################
# 3) PARAMETERS
###############################################################################

# NULL = use all automatically detected modules except excluded modules.
modules_of_interest <- NULL
exclude_modules <- c("grey")

top_hubs_per_module <- 30
top_plot_hubs <- 15

# Enrichment strategy:
# retrieve terms permissively, then classify by adjusted P value.
enrich_retrieve_p_cutoff <- 1
enrich_retrieve_q_cutoff <- 1

enrich_sig_fdr_cutoff <- 0.05
enrich_suggestive_fdr_cutoff <- 0.10

min_genes_for_enrichment <- 10
show_terms_per_plot <- 12

# Validated dimensions and DEP counts.
EXPECTED_N_SAMPLES <- 639L
EXPECTED_N_GENES <- 9638L
EXPECTED_N_MODULES_EXCLUDING_GREY <- 8L
EXPECTED_PRIMARY_DEP_GENES <- 587L
EXPECTED_ANALYTE_MATCHED_DEP_GENES <- 544L

REQUIRE_EXPECTED_DIMENSIONS <- TRUE
REQUIRE_EXPECTED_DEP_COUNTS <- TRUE
REQUIRE_ONE_ASSIGNED_KME_PER_GENE <- TRUE

DEP_FDR_CUTOFF <- 0.05

# A module occupying at least this fraction is highlighted for additional audit.
LARGE_MODULE_PROP_THRESHOLD <- 0.50

###############################################################################
# 4) COLORS AND HELPERS
###############################################################################

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
  turquoise   = "#40E0D0",
  cyan        = "#00BFC4",
  tan         = "#D2B48C",
  salmon      = "#FA8072",
  midnightblue = "#191970",
  lightcyan   = "#E0FFFF",
  grey        = "#9E9E9E"
)

get_module_color <- function(module_name) {
  module_name <- as.character(module_name)
  cols <- MODULE_COLORS[module_name]
  cols[is.na(cols)] <- "#7F7F7F"
  unname(cols)
}

safe_write_csv <- function(x, file) {
  dir.create(
    dirname(file),
    recursive = TRUE,
    showWarnings = FALSE
  )
  readr::write_csv(x, file)
}

safe_file_tag <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9_\\-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

clean_text_na <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c(
    "",
    "NA",
    "NaN",
    "NULL",
    "null",
    "N/A",
    "nan"
  )] <- NA_character_
  x
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(x))
}

ensure_column <- function(df, col, default = NA_character_) {
  if (!col %in% names(df)) {
    df[[col]] <- default
  }
  df
}

first_existing_value <- function(x) {
  x <- clean_text_na(x)
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  sort(unique(x))[[1]]
}

normalize_entrez_id <- function(x) {
  x <- clean_text_na(x)
  x <- gsub("\\.0$", "", x)
  x[!grepl("^[0-9]+$", x)] <- NA_character_
  x
}

safe_quantile <- function(x, prob) {
  x <- safe_numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  unname(stats::quantile(
    x,
    probs = prob,
    na.rm = TRUE
  ))
}

safe_fisher_module <- function(
    module_vector,
    dep_flag,
    target_module,
    definition_label
) {
  module_vector <- as.character(module_vector)
  dep_flag <- as.logical(dep_flag)

  keep <- !is.na(module_vector) & !is.na(dep_flag)
  module_vector <- module_vector[keep]
  dep_flag <- dep_flag[keep]

  in_module <- module_vector == target_module

  a <- sum(in_module & dep_flag)
  b <- sum(in_module & !dep_flag)
  c <- sum(!in_module & dep_flag)
  d <- sum(!in_module & !dep_flag)

  mat <- matrix(
    c(a, b, c, d),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      Module_membership = c("In_module", "Outside_module"),
      DEP_status = c("DEP", "Not_DEP")
    )
  )

  ft <- tryCatch(
    stats::fisher.test(
      mat,
      alternative = "greater"
    ),
    error = function(e) NULL
  )

  tibble::tibble(
    Module = target_module,
    DEP_definition = definition_label,
    n_module_DEP = a,
    n_module_not_DEP = b,
    n_outside_DEP = c,
    n_outside_not_DEP = d,
    odds_ratio = if (is.null(ft)) {
      NA_real_
    } else {
      unname(ft$estimate)
    },
    p_value = if (is.null(ft)) {
      NA_real_
    } else {
      ft$p.value
    }
  )
}

plot_hub_bar <- function(
    df_hubs,
    module_name,
    n_plot = 15,
    outdir_fig
) {
  if (is.null(df_hubs) || nrow(df_hubs) == 0) {
    return(NULL)
  }

  module_col <- get_module_color(module_name)

  df_plot <- df_hubs %>%
    dplyr::slice(
      seq_len(min(n_plot, nrow(df_hubs)))
    ) %>%
    dplyr::mutate(
      Protein_Label_raw = dplyr::coalesce(
        clean_text_na(Protein_Display),
        clean_text_na(Protein_Name),
        clean_text_na(EntrezGeneSymbol)
      ),
      Protein_Label = make.unique(
        as.character(Protein_Label_raw),
        sep = "_"
      ),
      Protein_Label = factor(
        Protein_Label,
        levels = rev(Protein_Label)
      )
    )

  p <- ggplot(
    df_plot,
    aes(
      x = Protein_Label,
      y = abs_kME
    )
  ) +
    geom_col(
      fill = module_col,
      width = 0.8
    ) +
    coord_flip() +
    labs(
      title = paste0(
        "Top hubs — module ",
        module_name
      ),
      x = NULL,
      y = "|kME|"
    ) +
    theme_bw(
      base_size = 12
    ) +
    theme(
      plot.title = element_text(
        face = "bold"
      ),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )

  ggsave(
    filename = file.path(
      outdir_fig,
      paste0(
        "top_hubs_barplot_",
        module_name,
        ".pdf"
      )
    ),
    plot = p,
    width = 7,
    height = 5
  )

  ggsave(
    filename = file.path(
      outdir_fig,
      paste0(
        "top_hubs_barplot_",
        module_name,
        ".png"
      )
    ),
    plot = p,
    width = 7,
    height = 5,
    dpi = 300
  )

  p
}

plot_enrichment_bar <- function(
    res,
    module_name,
    analysis_name,
    outdir_fig,
    n_show = 12
) {
  if (is.null(res) || nrow(res) == 0) {
    return(NULL)
  }

  module_col <- get_module_color(module_name)

  top_res <- res %>%
    dplyr::arrange(
      p.adjust,
      pvalue
    ) %>%
    dplyr::slice(
      seq_len(min(n_show, nrow(res)))
    ) %>%
    dplyr::mutate(
      Description_raw = stringr::str_wrap(
        as.character(Description),
        width = 45
      ),
      Description_plot = make.unique(
        Description_raw,
        sep = "_"
      ),
      Description_plot = factor(
        Description_plot,
        levels = rev(Description_plot)
      ),
      minus_log10_padj = -log10(
        pmax(
          p.adjust,
          .Machine$double.xmin
        )
      )
    )

  p <- ggplot(
    top_res,
    aes(
      x = Description_plot,
      y = minus_log10_padj
    )
  ) +
    geom_col(
      fill = module_col,
      width = 0.8
    ) +
    coord_flip() +
    labs(
      title = paste0(
        analysis_name,
        " enrichment — ",
        module_name
      ),
      x = NULL,
      y = expression(-log[10]("FDR"))
    ) +
    theme_bw(
      base_size = 12
    ) +
    theme(
      plot.title = element_text(
        face = "bold"
      ),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )

  ggsave(
    filename = file.path(
      outdir_fig,
      paste0(
        "barplot_",
        analysis_name,
        "_",
        module_name,
        ".pdf"
      )
    ),
    plot = p,
    width = 8,
    height = 5.5
  )

  ggsave(
    filename = file.path(
      outdir_fig,
      paste0(
        "barplot_",
        analysis_name,
        "_",
        module_name,
        ".png"
      )
    ),
    plot = p,
    width = 8,
    height = 5.5,
    dpi = 300
  )

  p
}

empty_enrichment_table <- function() {
  tibble::tibble(
    ID = character(),
    Description = character(),
    GeneRatio = character(),
    BgRatio = character(),
    pvalue = numeric(),
    p.adjust = numeric(),
    qvalue = numeric(),
    geneID = character(),
    Count = integer()
  )
}

save_enrichment_tables <- function(
    res,
    module_name,
    analysis_name
) {
  if (is.null(res)) {
    res_df <- empty_enrichment_table()
  } else {
    res_df <- as.data.frame(res)

    if (nrow(res_df) == 0) {
      res_df <- empty_enrichment_table()
    }
  }

  if (nrow(res_df) > 0) {
    res_df <- res_df %>%
      dplyr::mutate(
        pvalue = safe_numeric(pvalue),
        p.adjust = safe_numeric(p.adjust)
      ) %>%
      dplyr::arrange(
        p.adjust,
        pvalue
      )
  }

  all_file <- file.path(
    outdir,
    "enrichment",
    analysis_name,
    paste0(
      "enrichment_",
      analysis_name,
      "_",
      module_name,
      "_all_terms.csv"
    )
  )

  safe_write_csv(
    res_df,
    all_file
  )

  if (nrow(res_df) > 0) {
    sig <- res_df %>%
      dplyr::filter(
        !is.na(p.adjust),
        p.adjust < enrich_sig_fdr_cutoff
      )

    top <- res_df %>%
      dplyr::slice(
        seq_len(min(
          show_terms_per_plot,
          nrow(res_df)
        ))
      )
  } else {
    sig <- res_df
    top <- res_df
  }

  safe_write_csv(
    sig,
    file.path(
      outdir,
      "enrichment",
      "significant",
      paste0(
        "enrichment_",
        analysis_name,
        "_",
        module_name,
        "_FDR005.csv"
      )
    )
  )

  safe_write_csv(
    top,
    file.path(
      outdir,
      "enrichment",
      "top_terms",
      paste0(
        "enrichment_",
        analysis_name,
        "_",
        module_name,
        "_top_terms.csv"
      )
    )
  )

  if (nrow(top) > 0) {
    plot_enrichment_bar(
      res = top,
      module_name = module_name,
      analysis_name = analysis_name,
      outdir_fig = file.path(
        outdir,
        "figures",
        "enrichment"
      ),
      n_show = show_terms_per_plot
    )
  }

  res_df
}

make_enrich_count <- function(res) {
  if (is.null(res) || nrow(res) == 0) {
    return(tibble::tibble(
      n_terms_total = 0L,
      n_terms_fdr_0_05 = 0L,
      n_terms_fdr_0_10 = 0L,
      min_FDR = NA_real_,
      top_term = NA_character_
    ))
  }

  res <- as.data.frame(res) %>%
    dplyr::arrange(
      p.adjust,
      pvalue
    )

  min_fdr <- suppressWarnings(
    min(
      res$p.adjust,
      na.rm = TRUE
    )
  )

  if (!is.finite(min_fdr)) {
    min_fdr <- NA_real_
  }

  tibble::tibble(
    n_terms_total = nrow(res),
    n_terms_fdr_0_05 = sum(
      !is.na(res$p.adjust) &
        res$p.adjust < enrich_sig_fdr_cutoff
    ),
    n_terms_fdr_0_10 = sum(
      !is.na(res$p.adjust) &
        res$p.adjust < enrich_suggestive_fdr_cutoff
    ),
    min_FDR = min_fdr,
    top_term = if (nrow(res) > 0) {
      as.character(res$Description[[1]])
    } else {
      NA_character_
    }
  )
}

###############################################################################
# 5) LOAD WORKSPACE AND TABLES
###############################################################################

load(workspace_file)

gene_module_assignment_tbl <- readr::read_csv(
  gene_module_file,
  show_col_types = FALSE,
  guess_max = 100000
)

cat(
  "Workspace loaded from Script 11:\n",
  workspace_file,
  "\n\n"
)

cat("Objects available:\n")
cat(" - datExpr_clean:", exists("datExpr_clean"), "\n")
cat(" - mergedMEs:", exists("mergedMEs"), "\n")
cat(" - mergedColors:", exists("mergedColors"), "\n")
cat(
  " - gene_module_assignment:",
  exists("gene_module_assignment"),
  "\n\n"
)

if (!exists("datExpr_clean")) {
  stop(
    "datExpr_clean is absent from the Script 11 workspace.",
    call. = FALSE
  )
}

if (!exists("mergedMEs")) {
  stop(
    "mergedMEs is absent from the Script 11 workspace.",
    call. = FALSE
  )
}

if (!exists("mergedColors")) {
  stop(
    "mergedColors is absent from the Script 11 workspace.",
    call. = FALSE
  )
}

required_module_cols <- c(
  "EntrezGeneSymbol",
  "Module"
)

missing_module_cols <- setdiff(
  required_module_cols,
  names(gene_module_assignment_tbl)
)

if (length(missing_module_cols) > 0) {
  stop(
    "gene_module_assignment.csv lacks required columns: ",
    paste(
      missing_module_cols,
      collapse = ", "
    ),
    call. = FALSE
  )
}

###############################################################################
# 6) BASIC CHECKS AND ANNOTATION HARMONIZATION
###############################################################################

datExpr_clean <- as.data.frame(
  datExpr_clean,
  check.names = FALSE
)

mergedMEs <- as.data.frame(
  mergedMEs,
  check.names = FALSE
)

mergedColors <- as.character(
  mergedColors
)

if (is.null(names(mergedColors))) {
  names(mergedColors) <- colnames(datExpr_clean)
}

if (nrow(datExpr_clean) != nrow(mergedMEs)) {
  stop(
    "datExpr_clean and mergedMEs do not contain the same number of samples.",
    call. = FALSE
  )
}

if (!setequal(
  rownames(datExpr_clean),
  rownames(mergedMEs)
)) {
  stop(
    "datExpr_clean and mergedMEs do not contain the same sample set.",
    call. = FALSE
  )
}

if (!identical(
  rownames(datExpr_clean),
  rownames(mergedMEs)
)) {
  mergedMEs <- mergedMEs[
    rownames(datExpr_clean),
    ,
    drop = FALSE
  ]
}

if (!identical(
  names(mergedColors),
  colnames(datExpr_clean)
)) {
  if (!setequal(
    names(mergedColors),
    colnames(datExpr_clean)
  )) {
    stop(
      "mergedColors and datExpr_clean do not contain the same gene set.",
      call. = FALSE
    )
  }

  mergedColors <- mergedColors[
    colnames(datExpr_clean)
  ]
}

gene_module_assignment_tbl <- gene_module_assignment_tbl %>%
  dplyr::mutate(
    EntrezGeneSymbol = clean_text_na(
      EntrezGeneSymbol
    ),
    Module = clean_text_na(
      Module
    )
  ) %>%
  dplyr::filter(
    !is.na(EntrezGeneSymbol),
    !is.na(Module)
  )

if (
  anyDuplicated(
    gene_module_assignment_tbl$EntrezGeneSymbol
  ) > 0
) {
  stop(
    "gene_module_assignment.csv contains duplicated genes.",
    call. = FALSE
  )
}

if (!setequal(
  gene_module_assignment_tbl$EntrezGeneSymbol,
  colnames(datExpr_clean)
)) {
  stop(
    "gene_module_assignment.csv and datExpr_clean do not contain the same gene set.",
    call. = FALSE
  )
}

gene_module_assignment_tbl <- gene_module_assignment_tbl[
  match(
    colnames(datExpr_clean),
    gene_module_assignment_tbl$EntrezGeneSymbol
  ),
  ,
  drop = FALSE
]

if (!all(
  gene_module_assignment_tbl$EntrezGeneSymbol ==
    colnames(datExpr_clean)
)) {
  stop(
    "Gene annotation could not be aligned to datExpr_clean.",
    call. = FALSE
  )
}

if (!all(
  gene_module_assignment_tbl$Module ==
    unname(mergedColors)
)) {
  stop(
    "Module labels in gene_module_assignment.csv do not match mergedColors.",
    call. = FALSE
  )
}

# Preserve original annotations while ensuring expected columns exist.
default_character_cols <- c(
  "Protein_Name",
  "Protein_Display",
  "AptName",
  "Representative_AptName",
  "SeqId",
  "EntrezGeneID",
  "TargetFullName",
  "Target",
  "UniProt",
  "Primary_DEP_AptName",
  "Primary_DEP_SeqId",
  "Primary_DEP_Protein_Name",
  "type",
  "Direction",
  "Primary_DEP_type",
  "Primary_DEP_Direction"
)

for (cc in default_character_cols) {
  gene_module_assignment_tbl <- ensure_column(
    gene_module_assignment_tbl,
    cc,
    NA_character_
  )
}

default_numeric_cols <- c(
  "N_genes",
  "Prop_total",
  "logFC",
  "P.Value",
  "adj.P.Val",
  "Primary_DEP_logFC",
  "Primary_DEP_P.Value",
  "Primary_DEP_adj.P.Val"
)

for (cc in default_numeric_cols) {
  gene_module_assignment_tbl <- ensure_column(
    gene_module_assignment_tbl,
    cc,
    NA_real_
  )

  gene_module_assignment_tbl[[cc]] <- safe_numeric(
    gene_module_assignment_tbl[[cc]]
  )
}

if (all(is.na(
  gene_module_assignment_tbl$Representative_AptName
))) {
  gene_module_assignment_tbl$Representative_AptName <-
    gene_module_assignment_tbl$AptName
}

gene_module_assignment_tbl <- gene_module_assignment_tbl %>%
  dplyr::mutate(
    Protein_Name = dplyr::coalesce(
      clean_text_na(Protein_Name),
      EntrezGeneSymbol
    ),
    Protein_Display = dplyr::coalesce(
      clean_text_na(Protein_Display),
      Protein_Name,
      EntrezGeneSymbol
    ),
    Representative_AptName = dplyr::coalesce(
      clean_text_na(Representative_AptName),
      clean_text_na(AptName)
    ),
    analyte_matched_DEP = (
      !is.na(adj.P.Val) &
        adj.P.Val < DEP_FDR_CUTOFF
    ),
    primary_DEP_overlap = (
      !is.na(Primary_DEP_adj.P.Val) &
        Primary_DEP_adj.P.Val < DEP_FDR_CUTOFF
    ),
    analyte_matched_direction = dplyr::case_when(
      analyte_matched_DEP & logFC > 0 ~ "Higher_in_AD",
      analyte_matched_DEP & logFC < 0 ~ "Lower_in_AD",
      TRUE ~ "Not_DEP"
    ),
    primary_DEP_direction = dplyr::case_when(
      primary_DEP_overlap &
        Primary_DEP_logFC > 0 ~ "Higher_in_AD",
      primary_DEP_overlap &
        Primary_DEP_logFC < 0 ~ "Lower_in_AD",
      TRUE ~ "Not_DEP"
    )
  )

me_cols <- colnames(mergedMEs)
module_colors_from_MEs <- gsub(
  "^ME",
  "",
  me_cols
)

module_colors_from_assignment <- unique(
  gene_module_assignment_tbl$Module
)

if (!setequal(
  module_colors_from_MEs,
  module_colors_from_assignment
)) {
  stop(
    "Module labels derived from mergedMEs do not match gene-module assignments.",
    call. = FALSE
  )
}

cat("Modules detected in mergedMEs:\n")
print(module_colors_from_MEs)

if (is.null(modules_of_interest)) {
  modules_of_interest <- setdiff(
    module_colors_from_MEs,
    exclude_modules
  )
} else {
  modules_of_interest <- modules_of_interest[
    modules_of_interest %in%
      module_colors_from_MEs
  ]

  modules_of_interest <- setdiff(
    modules_of_interest,
    exclude_modules
  )
}

if (length(modules_of_interest) == 0) {
  stop(
    "No valid modules remained after exclusions.",
    call. = FALSE
  )
}

n_primary_dep_observed <- sum(
  gene_module_assignment_tbl$primary_DEP_overlap,
  na.rm = TRUE
)

n_analyte_dep_observed <- sum(
  gene_module_assignment_tbl$analyte_matched_DEP,
  na.rm = TRUE
)

if (REQUIRE_EXPECTED_DIMENSIONS) {
  if (nrow(datExpr_clean) != EXPECTED_N_SAMPLES) {
    stop(
      "Expected ",
      EXPECTED_N_SAMPLES,
      " samples, but found ",
      nrow(datExpr_clean),
      ".",
      call. = FALSE
    )
  }

  if (ncol(datExpr_clean) != EXPECTED_N_GENES) {
    stop(
      "Expected ",
      EXPECTED_N_GENES,
      " genes, but found ",
      ncol(datExpr_clean),
      ".",
      call. = FALSE
    )
  }

  if (
    length(modules_of_interest) !=
      EXPECTED_N_MODULES_EXCLUDING_GREY
  ) {
    stop(
      "Expected ",
      EXPECTED_N_MODULES_EXCLUDING_GREY,
      " non-grey modules, but found ",
      length(modules_of_interest),
      ".",
      call. = FALSE
    )
  }
}

if (REQUIRE_EXPECTED_DEP_COUNTS) {
  if (
    n_primary_dep_observed !=
      EXPECTED_PRIMARY_DEP_GENES
  ) {
    stop(
      "Expected ",
      EXPECTED_PRIMARY_DEP_GENES,
      " primary DEP genes, but observed ",
      n_primary_dep_observed,
      ".",
      call. = FALSE
    )
  }

  if (
    n_analyte_dep_observed !=
      EXPECTED_ANALYTE_MATCHED_DEP_GENES
  ) {
    stop(
      "Expected ",
      EXPECTED_ANALYTE_MATCHED_DEP_GENES,
      " analyte-matched DEP genes, but observed ",
      n_analyte_dep_observed,
      ".",
      call. = FALSE
    )
  }
}

cat("\nModules analyzed in Script 12:\n")
print(modules_of_interest)
cat("\n")

input_audit <- tibble::tibble(
  metric = c(
    "workspace_file",
    "n_samples_datExpr_clean",
    "n_genes_datExpr_clean",
    "n_module_eigengenes",
    "n_modules_excluding_grey",
    "n_gene_module_rows",
    "n_unique_genes",
    "n_primary_DEP_genes",
    "n_analyte_matched_DEP_genes",
    "gene_order_aligned",
    "module_labels_aligned",
    "sample_order_aligned",
    "wgcna_input_level"
  ),
  value = c(
    workspace_file,
    as.character(nrow(datExpr_clean)),
    as.character(ncol(datExpr_clean)),
    as.character(ncol(mergedMEs)),
    as.character(length(modules_of_interest)),
    as.character(nrow(gene_module_assignment_tbl)),
    as.character(dplyr::n_distinct(
      gene_module_assignment_tbl$EntrezGeneSymbol
    )),
    as.character(n_primary_dep_observed),
    as.character(n_analyte_dep_observed),
    as.character(all(
      gene_module_assignment_tbl$EntrezGeneSymbol ==
        colnames(datExpr_clean)
    )),
    as.character(all(
      gene_module_assignment_tbl$Module ==
        unname(mergedColors)
    )),
    as.character(all(
      rownames(datExpr_clean) ==
        rownames(mergedMEs)
    )),
    "GENE-COLLAPSED, outcome-independent SOMAmer selection"
  )
)

safe_write_csv(
  input_audit,
  file.path(
    outdir,
    "qc",
    "script12_input_alignment_audit.csv"
  )
)

###############################################################################
# 7) kME / MODULE MEMBERSHIP
###############################################################################

kME_mat <- WGCNA::signedKME(
  datExpr = datExpr_clean,
  datME = mergedMEs,
  exprWeights = NULL,
  MEWeights = NULL,
  outputColumnName = "kME",
  corFnc = "cor",
  corOptions = "use = 'pairwise.complete.obs'"
)

if (is.null(rownames(kME_mat))) {
  rownames(kME_mat) <- colnames(datExpr_clean)
}

if (!identical(
  rownames(kME_mat),
  colnames(datExpr_clean)
)) {
  if (!setequal(
    rownames(kME_mat),
    colnames(datExpr_clean)
  )) {
    stop(
      "signedKME output does not contain the WGCNA gene set.",
      call. = FALSE
    )
  }

  kME_mat <- kME_mat[
    colnames(datExpr_clean),
    ,
    drop = FALSE
  ]
}

kME_tbl <- tibble::tibble(
  EntrezGeneSymbol = rownames(kME_mat)
) %>%
  dplyr::bind_cols(
    as.data.frame(
      kME_mat,
      check.names = FALSE
    )
  )

n_samples_kme <- nrow(datExpr_clean)

kME_long <- kME_tbl %>%
  tidyr::pivot_longer(
    cols = starts_with("kME"),
    names_to = "ME_col",
    values_to = "kME"
  ) %>%
  dplyr::mutate(
    Module_from_ME = gsub(
      "^kME",
      "",
      ME_col
    ),
    abs_kME = abs(kME)
  )

kME_annot <- kME_long %>%
  dplyr::left_join(
    gene_module_assignment_tbl,
    by = "EntrezGeneSymbol"
  ) %>%
  dplyr::mutate(
    is_assigned_module = (
      Module == Module_from_ME
    )
  )

kME_assigned <- kME_annot %>%
  dplyr::filter(
    is_assigned_module
  ) %>%
  dplyr::mutate(
    t_stat = kME *
      sqrt(
        (n_samples_kme - 2) /
          pmax(
            1e-12,
            1 - kME^2
          )
      ),
    p_value = 2 *
      stats::pt(
        -abs(t_stat),
        df = n_samples_kme - 2
      )
  ) %>%
  dplyr::group_by(
    Module
  ) %>%
  dplyr::mutate(
    FDR_within_module = p.adjust(
      p_value,
      method = "BH"
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(
    Module,
    dplyr::desc(abs_kME),
    p_value
  )

if (
  REQUIRE_ONE_ASSIGNED_KME_PER_GENE &&
  nrow(kME_assigned) != ncol(datExpr_clean)
) {
  stop(
    "Expected one assigned-module kME row per gene, but found ",
    nrow(kME_assigned),
    " rows for ",
    ncol(datExpr_clean),
    " genes.",
    call. = FALSE
  )
}

if (
  anyDuplicated(
    kME_assigned$EntrezGeneSymbol
  ) > 0
) {
  stop(
    "Assigned-module kME table contains duplicated genes.",
    call. = FALSE
  )
}

if (any(!is.finite(kME_assigned$kME))) {
  stop(
    "Assigned-module kME table contains non-finite kME values.",
    call. = FALSE
  )
}

safe_write_csv(
  kME_tbl,
  file.path(
    outdir,
    "tables",
    "full_kME_table_wide.csv"
  )
)

safe_write_csv(
  kME_assigned,
  file.path(
    outdir,
    "tables",
    "full_kME_assigned_module_long.csv"
  )
)

###############################################################################
# 8) HUB SUMMARY
###############################################################################

hub_summary <- kME_assigned %>%
  dplyr::group_by(
    Module
  ) %>%
  dplyr::summarise(
    n_genes = dplyr::n(),
    prop_total_genes = dplyr::n() /
      ncol(datExpr_clean),
    mean_kME = mean(
      kME,
      na.rm = TRUE
    ),
    mean_abs_kME = mean(
      abs_kME,
      na.rm = TRUE
    ),
    median_abs_kME = median(
      abs_kME,
      na.rm = TRUE
    ),
    q25_abs_kME = safe_quantile(
      abs_kME,
      0.25
    ),
    q75_abs_kME = safe_quantile(
      abs_kME,
      0.75
    ),
    max_abs_kME = max(
      abs_kME,
      na.rm = TRUE
    ),
    n_hubs_abs_kME_ge_0_90 = sum(
      abs_kME >= 0.90,
      na.rm = TRUE
    ),
    n_hubs_abs_kME_ge_0_80 = sum(
      abs_kME >= 0.80,
      na.rm = TRUE
    ),
    n_hubs_abs_kME_ge_0_70 = sum(
      abs_kME >= 0.70,
      na.rm = TRUE
    ),
    prop_hubs_abs_kME_ge_0_80 = mean(
      abs_kME >= 0.80,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    dplyr::desc(n_genes)
  )

safe_write_csv(
  hub_summary,
  file.path(
    outdir,
    "tables",
    "module_hub_summary.csv"
  )
)

###############################################################################
# 9) DUAL DEP BURDEN AND OVERREPRESENTATION
###############################################################################

module_dep_burden <- kME_assigned %>%
  dplyr::group_by(
    Module
  ) %>%
  dplyr::summarise(
    n_genes = dplyr::n(),

    analyte_matched_DEP_n = sum(
      analyte_matched_DEP,
      na.rm = TRUE
    ),
    analyte_matched_DEP_prop = mean(
      analyte_matched_DEP,
      na.rm = TRUE
    ),
    analyte_matched_higher_in_AD_n = sum(
      analyte_matched_DEP &
        logFC > 0,
      na.rm = TRUE
    ),
    analyte_matched_lower_in_AD_n = sum(
      analyte_matched_DEP &
        logFC < 0,
      na.rm = TRUE
    ),
    analyte_matched_mean_abs_logFC = mean(
      abs(logFC),
      na.rm = TRUE
    ),

    primary_DEP_overlap_n = sum(
      primary_DEP_overlap,
      na.rm = TRUE
    ),
    primary_DEP_overlap_prop = mean(
      primary_DEP_overlap,
      na.rm = TRUE
    ),
    primary_DEP_higher_in_AD_n = sum(
      primary_DEP_overlap &
        Primary_DEP_logFC > 0,
      na.rm = TRUE
    ),
    primary_DEP_lower_in_AD_n = sum(
      primary_DEP_overlap &
        Primary_DEP_logFC < 0,
      na.rm = TRUE
    ),
    primary_DEP_mean_abs_logFC = mean(
      abs(Primary_DEP_logFC),
      na.rm = TRUE
    ),

    same_SOMAmer_as_primary_DEP_n = sum(
      same_somamer_as_primary_DEP %in% TRUE,
      na.rm = TRUE
    ),
    different_SOMAmer_from_primary_DEP_n = sum(
      same_somamer_as_primary_DEP %in% FALSE,
      na.rm = TRUE
    ),

    mean_abs_kME = mean(
      abs_kME,
      na.rm = TRUE
    ),
    median_abs_kME = median(
      abs_kME,
      na.rm = TRUE
    ),

    .groups = "drop"
  ) %>%
  dplyr::mutate(
    analyte_matched_DEP_definition = paste0(
      "Exact WGCNA SOMAmer; FDR < ",
      DEP_FDR_CUTOFF
    ),
    primary_DEP_overlap_definition = paste0(
      "Canonical primary DEP gene map; FDR < ",
      DEP_FDR_CUTOFF
    )
  ) %>%
  dplyr::arrange(
    dplyr::desc(primary_DEP_overlap_prop),
    dplyr::desc(analyte_matched_DEP_prop)
  )

safe_write_csv(
  module_dep_burden,
  file.path(
    outdir,
    "tables",
    "module_DEP_burden_dual_definition.csv"
  )
)

safe_write_csv(
  kME_assigned %>%
    dplyr::select(
      EntrezGeneSymbol,
      Module,
      Protein_Name,
      Protein_Display,
      AptName,
      Representative_AptName,
      Primary_DEP_AptName,
      same_somamer_as_primary_DEP,
      kME,
      abs_kME,
      analyte_matched_DEP,
      logFC,
      P.Value,
      adj.P.Val,
      analyte_matched_direction,
      primary_DEP_overlap,
      Primary_DEP_logFC,
      Primary_DEP_P.Value,
      Primary_DEP_adj.P.Val,
      primary_DEP_direction,
      dplyr::everything()
    ),
  file.path(
    outdir,
    "tables",
    "gene_module_kME_DEP_dual_annotation.csv"
  )
)

dep_overrepresentation <- purrr::map_dfr(
  modules_of_interest,
  function(mod) {
    dplyr::bind_rows(
      safe_fisher_module(
        module_vector = kME_assigned$Module,
        dep_flag = kME_assigned$analyte_matched_DEP,
        target_module = mod,
        definition_label = "Analyte-matched exact WGCNA SOMAmer"
      ),
      safe_fisher_module(
        module_vector = kME_assigned$Module,
        dep_flag = kME_assigned$primary_DEP_overlap,
        target_module = mod,
        definition_label = "Primary 587-gene DEP overlap"
      )
    )
  }
) %>%
  dplyr::group_by(
    DEP_definition
  ) %>%
  dplyr::mutate(
    FDR_across_modules = p.adjust(
      p_value,
      method = "BH"
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(
    DEP_definition,
    FDR_across_modules,
    dplyr::desc(odds_ratio)
  )

safe_write_csv(
  dep_overrepresentation,
  file.path(
    outdir,
    "tables",
    "module_DEP_overrepresentation_dual_definition.csv"
  )
)

###############################################################################
# 10) TOP HUBS BY MODULE
###############################################################################

hub_tables <- list()

for (mod in modules_of_interest) {
  df_mod <- kME_assigned %>%
    dplyr::filter(
      Module == mod
    ) %>%
    dplyr::arrange(
      dplyr::desc(abs_kME),
      p_value
    ) %>%
    dplyr::mutate(
      hub_rank = dplyr::row_number()
    )

  if (nrow(df_mod) == 0) {
    message(
      "Module without assigned genes in kME_assigned: ",
      mod
    )
    next
  }

  n_top <- min(
    top_hubs_per_module,
    nrow(df_mod)
  )

  top_mod <- df_mod %>%
    dplyr::slice(
      seq_len(n_top)
    )

  hub_tables[[mod]] <- top_mod

  safe_write_csv(
    df_mod,
    file.path(
      outdir,
      "hubs",
      paste0(
        "all_hubs_ranked_",
        safe_file_tag(mod),
        ".csv"
      )
    )
  )

  safe_write_csv(
    top_mod,
    file.path(
      outdir,
      "hubs",
      paste0(
        "top_hubs_",
        safe_file_tag(mod),
        ".csv"
      )
    )
  )

  plot_hub_bar(
    df_hubs = top_mod,
    module_name = mod,
    n_plot = top_plot_hubs,
    outdir_fig = file.path(
      outdir,
      "figures",
      "hubs"
    )
  )
}

top_hubs_all_modules <- dplyr::bind_rows(
  hub_tables,
  .id = "Module_list_name"
)

safe_write_csv(
  top_hubs_all_modules,
  file.path(
    outdir,
    "tables",
    "top_hubs_all_modules_combined.csv"
  )
)

###############################################################################
# 11) LARGEST-MODULE COHESION AUDIT
###############################################################################

largest_module_name <- hub_summary %>%
  dplyr::arrange(
    dplyr::desc(n_genes)
  ) %>%
  dplyr::slice(1) %>%
  dplyr::pull(Module)

largest_module_summary <- hub_summary %>%
  dplyr::filter(
    Module == largest_module_name
  ) %>%
  dplyr::left_join(
    module_dep_burden %>%
      dplyr::select(
        -n_genes,
        -mean_abs_kME,
        -median_abs_kME
      ),
    by = "Module"
  ) %>%
  dplyr::mutate(
    large_module_threshold = LARGE_MODULE_PROP_THRESHOLD,
    exceeds_large_module_threshold = (
      prop_total_genes >=
        LARGE_MODULE_PROP_THRESHOLD
    ),
    interpretation_note = paste(
      "Descriptive audit only.",
      "A large module is not automatically invalid;",
      "its cohesion, enrichment and trait associations",
      "must be interpreted jointly."
    )
  )

safe_write_csv(
  largest_module_summary,
  file.path(
    outdir,
    "tables",
    "largest_module_cohesion_audit.csv"
  )
)

safe_write_csv(
  kME_assigned %>%
    dplyr::filter(
      Module == largest_module_name
    ) %>%
    dplyr::arrange(
      dplyr::desc(abs_kME)
    ),
  file.path(
    outdir,
    "tables",
    paste0(
      "largest_module_",
      safe_file_tag(largest_module_name),
      "_gene_level_audit.csv"
    )
  )
)

###############################################################################
# 12) GENE ID MAPPING FOR ENRICHMENT
###############################################################################

gene_universe <- gene_module_assignment_tbl %>%
  dplyr::transmute(
    EntrezGeneSymbol = clean_text_na(
      EntrezGeneSymbol
    ),
    EntrezGeneID_native = normalize_entrez_id(
      EntrezGeneID
    ),
    Module = Module
  ) %>%
  dplyr::filter(
    !is.na(EntrezGeneSymbol),
    EntrezGeneSymbol != ""
  )

if (nrow(gene_universe) != EXPECTED_N_GENES) {
  stop(
    "The enrichment source universe does not contain exactly ",
    EXPECTED_N_GENES,
    " genes.",
    call. = FALSE
  )
}

symbol_to_entrez_raw <- tryCatch(
  clusterProfiler::bitr(
    gene_universe$EntrezGeneSymbol,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  ),
  error = function(e) {
    message(
      "SYMBOL -> ENTREZID mapping failed: ",
      e$message
    )
    NULL
  }
)

if (
  is.null(symbol_to_entrez_raw) ||
  nrow(symbol_to_entrez_raw) == 0
) {
  symbol_to_entrez <- tibble::tibble(
    EntrezGeneSymbol = character(),
    ENTREZID_map = character()
  )
} else {
  symbol_to_entrez <- symbol_to_entrez_raw %>%
    dplyr::transmute(
      EntrezGeneSymbol = clean_text_na(
        SYMBOL
      ),
      ENTREZID_map = normalize_entrez_id(
        ENTREZID
      )
    ) %>%
    dplyr::filter(
      !is.na(EntrezGeneSymbol),
      !is.na(ENTREZID_map)
    ) %>%
    dplyr::group_by(
      EntrezGeneSymbol
    ) %>%
    dplyr::summarise(
      ENTREZID_map = first_existing_value(
        ENTREZID_map
      ),
      .groups = "drop"
    )
}

gene_universe_mapped <- gene_universe %>%
  dplyr::left_join(
    symbol_to_entrez,
    by = "EntrezGeneSymbol"
  ) %>%
  dplyr::mutate(
    ENTREZID_final = dplyr::coalesce(
      EntrezGeneID_native,
      ENTREZID_map
    ),
    mapping_source = dplyr::case_when(
      !is.na(EntrezGeneID_native) ~ "native_annotation",
      is.na(EntrezGeneID_native) &
        !is.na(ENTREZID_map) ~ "org.Hs.eg.db_SYMBOL_mapping",
      TRUE ~ "unmapped"
    )
  )

background_entrez <- unique(
  gene_universe_mapped$ENTREZID_final[
    !is.na(
      gene_universe_mapped$ENTREZID_final
    )
  ]
)

if (
  length(background_entrez) <
    min_genes_for_enrichment
) {
  stop(
    "The mapped enrichment background is too small.",
    call. = FALSE
  )
}

mapping_audit <- tibble::tibble(
  metric = c(
    "n_WGCNA_source_genes",
    "n_WGCNA_source_modules",
    "n_genes_with_native_Entrez",
    "n_genes_mapped_by_org_Hs_eg_db",
    "n_genes_with_final_Entrez",
    "n_genes_without_final_Entrez",
    "mapping_proportion",
    "n_unique_background_Entrez",
    "background_definition"
  ),
  value = c(
    as.character(nrow(gene_universe_mapped)),
    as.character(dplyr::n_distinct(
      gene_universe_mapped$Module
    )),
    as.character(sum(
      gene_universe_mapped$mapping_source ==
        "native_annotation"
    )),
    as.character(sum(
      gene_universe_mapped$mapping_source ==
        "org.Hs.eg.db_SYMBOL_mapping"
    )),
    as.character(sum(
      !is.na(
        gene_universe_mapped$ENTREZID_final
      )
    )),
    as.character(sum(
      is.na(
        gene_universe_mapped$ENTREZID_final
      )
    )),
    as.character(mean(
      !is.na(
        gene_universe_mapped$ENTREZID_final
      )
    )),
    as.character(length(background_entrez)),
    paste(
      "All 9,638 outcome-independent WGCNA genes;",
      "enrichment functions use the uniquely mapped ENTREZ subset."
    )
  )
)

safe_write_csv(
  gene_universe_mapped,
  file.path(
    outdir,
    "tables",
    "gene_universe_symbol_to_entrez.csv"
  )
)

safe_write_csv(
  mapping_audit,
  file.path(
    outdir,
    "qc",
    "enrichment_background_mapping_audit.csv"
  )
)

safe_write_csv(
  gene_universe_mapped %>%
    dplyr::filter(
      is.na(ENTREZID_final)
    ),
  file.path(
    outdir,
    "qc",
    "enrichment_unmapped_WGCNA_genes.csv"
  )
)

###############################################################################
# 13) ENRICHMENT BY MODULE
###############################################################################

enrichment_summary <- list()
enrichment_errors <- list()
all_significant_terms <- list()

for (mod in modules_of_interest) {
  cat("\n============================\n")
  cat("Processing module:", mod, "\n")
  cat("============================\n")

  module_genes <- gene_module_assignment_tbl %>%
    dplyr::filter(
      Module == mod
    ) %>%
    dplyr::transmute(
      EntrezGeneSymbol = clean_text_na(
        EntrezGeneSymbol
      ),
      Protein_Name = Protein_Name,
      AptName = AptName,
      analyte_matched_DEP = analyte_matched_DEP,
      primary_DEP_overlap = primary_DEP_overlap
    ) %>%
    dplyr::filter(
      !is.na(EntrezGeneSymbol),
      EntrezGeneSymbol != ""
    ) %>%
    dplyr::distinct(
      EntrezGeneSymbol,
      .keep_all = TRUE
    )

  module_genes_mapped <- module_genes %>%
    dplyr::left_join(
      gene_universe_mapped %>%
        dplyr::select(
          EntrezGeneSymbol,
          ENTREZID_final,
          mapping_source
        ),
      by = "EntrezGeneSymbol"
    )

  safe_write_csv(
    module_genes_mapped,
    file.path(
      outdir,
      "tables",
      paste0(
        "genes_for_enrichment_",
        safe_file_tag(mod),
        ".csv"
      )
    )
  )

  module_entrez <- unique(
    module_genes_mapped$ENTREZID_final[
      !is.na(
        module_genes_mapped$ENTREZID_final
      )
    ]
  )

  enrich_info <- tibble::tibble(
    Module = mod,
    n_symbols = nrow(module_genes),
    n_mapped_genes = sum(
      !is.na(
        module_genes_mapped$ENTREZID_final
      )
    ),
    n_unique_mapped_entrez = length(module_entrez),
    mapping_proportion = mean(
      !is.na(
        module_genes_mapped$ENTREZID_final
      )
    ),
    n_background_entrez = length(background_entrez),
    background_source_genes = EXPECTED_N_GENES
  )

  safe_write_csv(
    enrich_info,
    file.path(
      outdir,
      "tables",
      paste0(
        "enrichment_input_summary_",
        safe_file_tag(mod),
        ".csv"
      )
    )
  )

  if (
    length(module_entrez) <
      min_genes_for_enrichment
  ) {
    message(
      "Module ",
      mod,
      " skipped because too few genes mapped to ENTREZID."
    )

    enrichment_summary[[mod]] <- tibble::tibble(
      Module = mod,
      n_symbols = nrow(module_genes),
      n_mapped_entrez = length(module_entrez),
      n_background = length(background_entrez),
      GO_BP_terms = 0L,
      GO_BP_terms_fdr010 = 0L,
      GO_BP_terms_total = 0L,
      GO_BP_min_FDR = NA_real_,
      GO_BP_top_term = NA_character_,
      KEGG_terms = 0L,
      KEGG_terms_fdr010 = 0L,
      KEGG_terms_total = 0L,
      KEGG_min_FDR = NA_real_,
      KEGG_top_term = NA_character_,
      Reactome_terms = 0L,
      Reactome_terms_fdr010 = 0L,
      Reactome_terms_total = 0L,
      Reactome_min_FDR = NA_real_,
      Reactome_top_term = NA_character_
    )

    next
  }

  # -------------------------------------------------------------------------
  # GO Biological Process
  # -------------------------------------------------------------------------

  ego_bp <- tryCatch(
    clusterProfiler::enrichGO(
      gene = module_entrez,
      universe = background_entrez,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = enrich_retrieve_p_cutoff,
      qvalueCutoff = enrich_retrieve_q_cutoff,
      readable = TRUE
    ),
    error = function(e) {
      enrichment_errors[[
        paste(
          mod,
          "GO_BP",
          sep = "_"
        )
      ]] <<- tibble::tibble(
        Module = mod,
        Library = "GO_BP",
        Error = e$message
      )

      message(
        "GO BP enrichment failed for module ",
        mod,
        ": ",
        e$message
      )

      NULL
    }
  )

  res_go <- save_enrichment_tables(
    res = ego_bp,
    module_name = mod,
    analysis_name = "GO_BP"
  )

  count_go <- make_enrich_count(
    res_go
  )

  if (nrow(res_go) > 0) {
    all_significant_terms[[
      paste(
        mod,
        "GO_BP",
        sep = "_"
      )
    ]] <- res_go %>%
      dplyr::filter(
        !is.na(p.adjust),
        p.adjust < enrich_sig_fdr_cutoff
      ) %>%
      dplyr::mutate(
        Module = mod,
        Library = "GO_BP",
        .before = 1
      )
  }

  # -------------------------------------------------------------------------
  # KEGG
  # -------------------------------------------------------------------------

  ekegg <- tryCatch(
    clusterProfiler::enrichKEGG(
      gene = module_entrez,
      universe = background_entrez,
      organism = "hsa",
      pAdjustMethod = "BH",
      pvalueCutoff = enrich_retrieve_p_cutoff,
      qvalueCutoff = enrich_retrieve_q_cutoff
    ),
    error = function(e) {
      enrichment_errors[[
        paste(
          mod,
          "KEGG",
          sep = "_"
        )
      ]] <<- tibble::tibble(
        Module = mod,
        Library = "KEGG",
        Error = e$message
      )

      message(
        "KEGG enrichment failed for module ",
        mod,
        ": ",
        e$message
      )

      NULL
    }
  )

  res_kegg <- save_enrichment_tables(
    res = ekegg,
    module_name = mod,
    analysis_name = "KEGG"
  )

  count_kegg <- make_enrich_count(
    res_kegg
  )

  if (nrow(res_kegg) > 0) {
    all_significant_terms[[
      paste(
        mod,
        "KEGG",
        sep = "_"
      )
    ]] <- res_kegg %>%
      dplyr::filter(
        !is.na(p.adjust),
        p.adjust < enrich_sig_fdr_cutoff
      ) %>%
      dplyr::mutate(
        Module = mod,
        Library = "KEGG",
        .before = 1
      )
  }

  # -------------------------------------------------------------------------
  # Reactome
  # -------------------------------------------------------------------------

  ereact <- tryCatch(
    ReactomePA::enrichPathway(
      gene = module_entrez,
      universe = background_entrez,
      organism = "human",
      pvalueCutoff = enrich_retrieve_p_cutoff,
      pAdjustMethod = "BH",
      qvalueCutoff = enrich_retrieve_q_cutoff,
      readable = TRUE
    ),
    error = function(e) {
      enrichment_errors[[
        paste(
          mod,
          "Reactome",
          sep = "_"
        )
      ]] <<- tibble::tibble(
        Module = mod,
        Library = "Reactome",
        Error = e$message
      )

      message(
        "Reactome enrichment failed for module ",
        mod,
        ": ",
        e$message
      )

      NULL
    }
  )

  res_react <- save_enrichment_tables(
    res = ereact,
    module_name = mod,
    analysis_name = "Reactome"
  )

  count_react <- make_enrich_count(
    res_react
  )

  if (nrow(res_react) > 0) {
    all_significant_terms[[
      paste(
        mod,
        "Reactome",
        sep = "_"
      )
    ]] <- res_react %>%
      dplyr::filter(
        !is.na(p.adjust),
        p.adjust < enrich_sig_fdr_cutoff
      ) %>%
      dplyr::mutate(
        Module = mod,
        Library = "Reactome",
        .before = 1
      )
  }

  enrichment_summary[[mod]] <- tibble::tibble(
    Module = mod,
    n_symbols = nrow(module_genes),
    n_mapped_entrez = length(module_entrez),
    n_background = length(background_entrez),

    GO_BP_terms = count_go$n_terms_fdr_0_05,
    GO_BP_terms_fdr010 = count_go$n_terms_fdr_0_10,
    GO_BP_terms_total = count_go$n_terms_total,
    GO_BP_min_FDR = count_go$min_FDR,
    GO_BP_top_term = count_go$top_term,

    KEGG_terms = count_kegg$n_terms_fdr_0_05,
    KEGG_terms_fdr010 = count_kegg$n_terms_fdr_0_10,
    KEGG_terms_total = count_kegg$n_terms_total,
    KEGG_min_FDR = count_kegg$min_FDR,
    KEGG_top_term = count_kegg$top_term,

    Reactome_terms = count_react$n_terms_fdr_0_05,
    Reactome_terms_fdr010 = count_react$n_terms_fdr_0_10,
    Reactome_terms_total = count_react$n_terms_total,
    Reactome_min_FDR = count_react$min_FDR,
    Reactome_top_term = count_react$top_term
  )
}

enrichment_summary_tbl <- dplyr::bind_rows(
  enrichment_summary
)

if (
  nrow(enrichment_summary_tbl) !=
    length(modules_of_interest)
) {
  stop(
    "Enrichment summary does not contain one row per analyzed module.",
    call. = FALSE
  )
}

safe_write_csv(
  enrichment_summary_tbl,
  file.path(
    outdir,
    "tables",
    "enrichment_summary_by_module.csv"
  )
)

if (length(enrichment_errors) > 0) {
  enrichment_errors_tbl <- dplyr::bind_rows(
    enrichment_errors
  )
} else {
  enrichment_errors_tbl <- tibble::tibble(
    Module = character(),
    Library = character(),
    Error = character()
  )
}

safe_write_csv(
  enrichment_errors_tbl,
  file.path(
    outdir,
    "tables",
    "enrichment_errors_by_module.csv"
  )
)

if (length(all_significant_terms) > 0) {
  enrichment_significant_all <- dplyr::bind_rows(
    all_significant_terms
  )
} else {
  enrichment_significant_all <- tibble::tibble(
    Module = character(),
    Library = character(),
    ID = character(),
    Description = character(),
    pvalue = numeric(),
    p.adjust = numeric()
  )
}

safe_write_csv(
  enrichment_significant_all,
  file.path(
    outdir,
    "tables",
    "enrichment_significant_terms_all_modules.csv"
  )
)

###############################################################################
# 14) OPTIONAL NETWORK TABLES
###############################################################################

for (mod in modules_of_interest) {
  if (!mod %in% names(hub_tables)) {
    next
  }

  node_tbl <- hub_tables[[mod]] %>%
    dplyr::mutate(
      node = EntrezGeneSymbol,
      label = dplyr::coalesce(
        clean_text_na(Protein_Display),
        clean_text_na(Protein_Name),
        clean_text_na(EntrezGeneSymbol)
      ),
      module = Module,
      module_color = get_module_color(
        Module
      )
    ) %>%
    dplyr::select(
      node,
      label,
      module,
      module_color,
      kME,
      abs_kME,
      hub_rank,
      analyte_matched_DEP,
      logFC,
      adj.P.Val,
      primary_DEP_overlap,
      Primary_DEP_logFC,
      Primary_DEP_adj.P.Val,
      AptName,
      Primary_DEP_AptName,
      same_somamer_as_primary_DEP,
      dplyr::everything()
    )

  safe_write_csv(
    node_tbl,
    file.path(
      outdir,
      "network_optional",
      paste0(
        "nodes_top_hubs_",
        safe_file_tag(mod),
        ".csv"
      )
    )
  )
}

###############################################################################
# 15) FINAL INTEGRATED MODULE SUMMARY
###############################################################################

integrated_module_summary <- hub_summary %>%
  dplyr::left_join(
    module_dep_burden %>%
      dplyr::select(
        -n_genes,
        -mean_abs_kME,
        -median_abs_kME
      ),
    by = "Module"
  ) %>%
  dplyr::left_join(
    enrichment_summary_tbl,
    by = "Module"
  ) %>%
  dplyr::left_join(
    dep_overrepresentation %>%
      dplyr::filter(
        DEP_definition ==
          "Analyte-matched exact WGCNA SOMAmer"
      ) %>%
      dplyr::select(
        Module,
        analyte_matched_DEP_OR = odds_ratio,
        analyte_matched_DEP_overrep_P = p_value,
        analyte_matched_DEP_overrep_FDR =
          FDR_across_modules
      ),
    by = "Module"
  ) %>%
  dplyr::left_join(
    dep_overrepresentation %>%
      dplyr::filter(
        DEP_definition ==
          "Primary 587-gene DEP overlap"
      ) %>%
      dplyr::select(
        Module,
        primary_DEP_overlap_OR = odds_ratio,
        primary_DEP_overlap_overrep_P = p_value,
        primary_DEP_overlap_overrep_FDR =
          FDR_across_modules
      ),
    by = "Module"
  ) %>%
  dplyr::mutate(
    total_significant_enrichment_terms =
      dplyr::coalesce(
        GO_BP_terms,
        0L
      ) +
      dplyr::coalesce(
        KEGG_terms,
        0L
      ) +
      dplyr::coalesce(
        Reactome_terms,
        0L
      ),
    is_largest_module = (
      Module == largest_module_name
    ),
    large_module_flag = (
      prop_total_genes >=
        LARGE_MODULE_PROP_THRESHOLD
    )
  ) %>%
  dplyr::arrange(
    dplyr::desc(
      primary_DEP_overlap_prop
    ),
    dplyr::desc(
      total_significant_enrichment_terms
    ),
    dplyr::desc(
      median_abs_kME
    )
  )

safe_write_csv(
  integrated_module_summary,
  file.path(
    outdir,
    "tables",
    "integrated_module_biology_summary.csv"
  )
)

###############################################################################
# 16) FINAL SUMMARY AND OUTPUT MANIFEST
###############################################################################

final_summary <- tibble::tibble(
  metric = c(
    "base_dir",
    "core_dir",
    "workspace_file",
    "outdir",
    "wgcna_input_level",
    "n_samples_used_for_kME",
    "n_genes_in_datExpr_clean",
    "n_module_eigengenes",
    "n_modules_analyzed",
    "module_names",
    "largest_module",
    "largest_module_n_genes",
    "largest_module_prop_total",
    "top_hubs_per_module_exported",
    "primary_DEP_gene_count",
    "analyte_matched_DEP_gene_count",
    "min_genes_for_enrichment",
    "n_source_background_genes",
    "n_unique_background_entrez",
    "enrichment_retrieve_p_cutoff",
    "enrichment_retrieve_q_cutoff",
    "enrichment_sig_fdr_cutoff",
    "enrichment_suggestive_fdr_cutoff",
    "n_enrichment_errors"
  ),
  value = c(
    base_dir,
    core_dir,
    workspace_file,
    outdir,
    "GENE-COLLAPSED, outcome-independent SOMAmer selection",
    as.character(nrow(datExpr_clean)),
    as.character(ncol(datExpr_clean)),
    as.character(ncol(mergedMEs)),
    as.character(length(modules_of_interest)),
    paste(
      modules_of_interest,
      collapse = ", "
    ),
    largest_module_name,
    as.character(
      hub_summary$n_genes[
        hub_summary$Module ==
          largest_module_name
      ]
    ),
    as.character(
      hub_summary$prop_total_genes[
        hub_summary$Module ==
          largest_module_name
      ]
    ),
    as.character(top_hubs_per_module),
    as.character(n_primary_dep_observed),
    as.character(n_analyte_dep_observed),
    as.character(min_genes_for_enrichment),
    as.character(nrow(gene_universe_mapped)),
    as.character(length(background_entrez)),
    as.character(enrich_retrieve_p_cutoff),
    as.character(enrich_retrieve_q_cutoff),
    as.character(enrich_sig_fdr_cutoff),
    as.character(enrich_suggestive_fdr_cutoff),
    as.character(nrow(enrichment_errors_tbl))
  )
)

safe_write_csv(
  final_summary,
  file.path(
    outdir,
    "tables",
    "script12_final_summary.csv"
  )
)

output_manifest <- tibble::tibble(
  output_file = c(
    "qc/script12_input_alignment_audit.csv",
    "qc/enrichment_background_mapping_audit.csv",
    "qc/enrichment_unmapped_WGCNA_genes.csv",
    "tables/full_kME_table_wide.csv",
    "tables/full_kME_assigned_module_long.csv",
    "tables/module_hub_summary.csv",
    "tables/module_DEP_burden_dual_definition.csv",
    "tables/module_DEP_overrepresentation_dual_definition.csv",
    "tables/gene_module_kME_DEP_dual_annotation.csv",
    "tables/largest_module_cohesion_audit.csv",
    "tables/gene_universe_symbol_to_entrez.csv",
    "tables/enrichment_summary_by_module.csv",
    "tables/enrichment_errors_by_module.csv",
    "tables/enrichment_significant_terms_all_modules.csv",
    "tables/integrated_module_biology_summary.csv",
    "tables/top_hubs_all_modules_combined.csv",
    "tables/script12_final_summary.csv",
    "workspace/script12_module_biology_workspace.RData",
    "sessionInfo.txt"
  ),
  description = c(
    "Strict Script 11 sample, gene, module and DEP-count audit.",
    "Mapping of the fixed 9,638-gene WGCNA universe to ENTREZ identifiers.",
    "WGCNA genes without a usable ENTREZ identifier.",
    "Wide kME matrix for every gene against every module eigengene.",
    "One assigned-module kME row per WGCNA gene, with complete annotation.",
    "Module size and kME-cohesion summary.",
    "Analyte-matched DEP burden and primary 587-gene DEP overlap by module.",
    "One-sided Fisher overrepresentation for both DEP definitions.",
    "Gene-level module, kME and dual DEP annotations.",
    "Additional descriptive audit of the largest module.",
    "Gene-symbol to final ENTREZ mapping for the WGCNA universe.",
    "GO BP, KEGG and Reactome term counts by module.",
    "Enrichment errors retained rather than silently discarded.",
    "All enrichment terms significant at FDR < 0.05.",
    "Integrated hubs, DEP burden and enrichment summary by module.",
    "Top hubs from all automatically detected modules.",
    "Script-level provenance and result summary.",
    "Workspace for subsequent module-trait integration and figures.",
    "R session information."
  )
)

safe_write_csv(
  output_manifest,
  file.path(
    outdir,
    "script12_output_manifest.csv"
  )
)

###############################################################################
# 17) SAVE WORKSPACE
###############################################################################

save(
  datExpr_clean,
  mergedMEs,
  mergedColors,
  gene_module_assignment_tbl,
  modules_of_interest,
  kME_mat,
  kME_tbl,
  kME_assigned,
  hub_summary,
  hub_tables,
  module_dep_burden,
  dep_overrepresentation,
  largest_module_name,
  largest_module_summary,
  gene_universe_mapped,
  background_entrez,
  mapping_audit,
  enrichment_summary_tbl,
  enrichment_errors_tbl,
  enrichment_significant_all,
  integrated_module_summary,
  final_summary,
  output_manifest,
  file = file.path(
    outdir,
    "workspace",
    "script12_module_biology_workspace.RData"
  )
)

writeLines(
  capture.output(
    utils::sessionInfo()
  ),
  con = file.path(
    outdir,
    "sessionInfo.txt"
  )
)

###############################################################################
# 18) FINAL MESSAGE
###############################################################################

cat("\nScript 12 finished successfully.\n")
cat("Main output directory:\n", outdir, "\n")
cat("Input level: GENE-COLLAPSED, outcome-independent SOMAmer selection.\n")
cat("Samples used for kME: ", nrow(datExpr_clean), "\n", sep = "")
cat("Genes used: ", ncol(datExpr_clean), "\n", sep = "")
cat("Modules analyzed: ", length(modules_of_interest), "\n", sep = "")
cat("Module labels: ", paste(modules_of_interest, collapse = ", "), "\n", sep = "")
cat("Largest module: ", largest_module_name, "\n", sep = "")
cat(
  "Largest module size: ",
  hub_summary$n_genes[
    hub_summary$Module ==
      largest_module_name
  ],
  "\n",
  sep = ""
)
cat("Primary DEP overlap count: ", n_primary_dep_observed, "\n", sep = "")
cat("Analyte-matched DEP count: ", n_analyte_dep_observed, "\n", sep = "")
cat("Mapped enrichment background ENTREZ IDs: ", length(background_entrez), "\n", sep = "")
cat("Enrichment errors recorded: ", nrow(enrichment_errors_tbl), "\n", sep = "")
cat("\nKey outputs:\n")
cat("- tables/full_kME_assigned_module_long.csv\n")
cat("- tables/module_hub_summary.csv\n")
cat("- tables/module_DEP_burden_dual_definition.csv\n")
cat("- tables/module_DEP_overrepresentation_dual_definition.csv\n")
cat("- tables/enrichment_summary_by_module.csv\n")
cat("- tables/integrated_module_biology_summary.csv\n")
cat("- tables/largest_module_cohesion_audit.csv\n")
cat("- hubs/top_hubs_<module>.csv\n")
cat("- enrichment/GO_BP/*_all_terms.csv\n")
cat("- enrichment/KEGG/*_all_terms.csv\n")
cat("- enrichment/Reactome/*_all_terms.csv\n")

###############################################################################
# END
###############################################################################


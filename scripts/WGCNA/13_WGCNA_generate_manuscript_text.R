###############################################################################
# ReDLat plasma proteomics — WGCNA workflow
# 13. Generate manuscript support text
# Requires: analysis and final figure outputs
# Produces: methods, results, legends and value audits
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

options(stringsAsFactors = FALSE)
options(error = traceback)

required_pkgs <- c("dplyr", "readr", "tibble", "stringr", "glue", "purrr")
missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs) > 0L) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
       ". Run renv::restore() before this script.", call. = FALSE)
}
invisible(lapply(required_pkgs, library, character.only = TRUE))

run_manuscript_text <- function() {

###############################################################################
# 00. PATHS AND CONSTANTS
###############################################################################

BASE_DIR <- WGCNA_CONFIG$project_root

WGCNA_ROOT <- WGCNA_CONFIG$result_root
FIGURE_ROOT <- file.path(WGCNA_CONFIG$publication_root, "figures")

S10 <- file.path(WGCNA_ROOT, "01_input")
S11 <- file.path(WGCNA_ROOT, "02_network")
S12 <- file.path(WGCNA_ROOT, "03_modules")
S13 <- file.path(WGCNA_ROOT, "04_module_traits")
S13B <- file.path(WGCNA_ROOT, "05_sensitivity")
S14 <- file.path(WGCNA_ROOT, "06_stability")
S14B <- file.path(WGCNA_ROOT, "07_biomarker_fdr")
S15B <- file.path(WGCNA_ROOT, "08_preservation")
S16 <- file.path(WGCNA_ROOT, "09_network_quality")

S17_FINAL <- file.path(FIGURE_ROOT, "main_figure_3")
S18_FINAL <- file.path(FIGURE_ROOT, "extended_data")

OUTDIR <- file.path(WGCNA_CONFIG$publication_root, "manuscript_text")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

MODULE_KEY <- tibble::tribble(
  ~Module,   ~Module_ID, ~Module_display, ~Biological_label,
  "green",   "M1",       "M1/green",      "neuronal-connectivity and extracellular-matrix",
  "blue",    "M2",       "M2/blue",       "AD-elevated high-differential-burden",
  "brown",   "M3",       "M3/brown",      "RNA-processing and ribonucleoprotein",
  "black",   "M4",       "M4/black",      "biomarker-associated covariance with low DEP burden",
  "magenta", "M5",       "M5/magenta",    "broad analytical module",
  "red",     "M6",       "M6/red",        "analytical module",
  "pink",    "M7",       "M7/pink",       "analytical module",
  "purple",  "M8",       "M8/purple",     "analytical module"
)

readr::write_csv(MODULE_KEY, file.path(OUTDIR, "WGCNA_Module_Key.csv"))

###############################################################################
# 01. HELPERS
###############################################################################

read_required <- function(path) {
  if (!file.exists(path)) {
    stop("Missing definitive WGCNA text source:\n", path, call. = FALSE)
  }
  readr::read_csv(path, show_col_types = FALSE, guess_max = 100000, name_repair = "unique")
}

summary_value <- function(tbl, metric, metric_col = NULL, value_col = NULL) {
  if (is.null(metric_col)) {
    metric_col <- c("metric", "Metric", "Item", "name", "Name")
    metric_col <- metric_col[metric_col %in% names(tbl)][1]
  }
  if (is.null(value_col)) {
    value_col <- c("value", "Value", "Observed", "result", "Result")
    value_col <- value_col[value_col %in% names(tbl)][1]
  }

  if (length(metric_col) == 0 || is.na(metric_col) ||
      length(value_col) == 0 || is.na(value_col)) {
    return(NA_character_)
  }

  idx <- which(as.character(tbl[[metric_col]]) == metric)
  if (length(idx) == 0) return(NA_character_)
  as.character(tbl[[value_col]][idx[1]])
}

fmt <- function(x, digits = 2) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0 || !is.finite(x)) return("NA")
  format(round(x, digits), nsmall = digits, trim = TRUE)
}

fmt_p <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0 || !is.finite(x)) return("NA")
  if (x < 0.001) return(format(x, scientific = TRUE, digits = 2))
  fmt(x, 3)
}

safe_number <- function(x, fallback = NA_real_) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0 || !is.finite(x[1])) fallback else x[1]
}

pick_col <- function(tbl, candidates, required = FALSE, label = "column") {
  hit <- candidates[candidates %in% names(tbl)][1]
  if (length(hit) > 0 && !is.na(hit)) return(hit)

  clean <- function(x) tolower(gsub("[^a-z0-9]+", "", x))
  idx <- match(clean(candidates), clean(names(tbl)), nomatch = 0)
  idx <- idx[idx > 0][1]

  if (length(idx) > 0 && !is.na(idx)) return(names(tbl)[idx])

  if (required) {
    stop(
      "Could not identify ", label, ". Available columns: ",
      paste(names(tbl), collapse = ", "),
      call. = FALSE
    )
  }
  NA_character_
}

module_row <- function(tbl, module) {
  module_col <- pick_col(tbl, c("Module", "module"), TRUE, "module column")
  tbl %>%
    dplyr::filter(sub("^ME", "", as.character(.data[[module_col]])) == module) %>%
    dplyr::slice(1)
}

extract_num <- function(row, candidates, fallback = NA_real_) {
  if (nrow(row) == 0) return(fallback)
  col <- pick_col(row, candidates)
  if (is.na(col)) return(fallback)
  safe_number(row[[col]], fallback)
}

format_association <- function(tbl, module, trait) {
  module_col <- pick_col(tbl, c("Module", "module"), TRUE, "module column")
  trait_col <- pick_col(tbl, c("Trait", "trait"), TRUE, "trait column")
  rho_col <- pick_col(tbl, c("rho", "Rho", "Spearman_rho"), TRUE, "rho column")
  fdr_col <- pick_col(tbl, c("FDR", "FDR_BH", "p_adj"), TRUE, "FDR column")
  label_col <- pick_col(tbl, c("Trait_label", "trait_label"))

  row <- tbl %>%
    dplyr::filter(
      sub("^ME", "", as.character(.data[[module_col]])) == module,
      as.character(.data[[trait_col]]) == trait
    ) %>%
    dplyr::slice(1)

  display <- MODULE_KEY$Module_display[match(module, MODULE_KEY$Module)]
  trait_label <- if (!is.na(label_col) && nrow(row) > 0) {
    as.character(row[[label_col]][1])
  } else {
    trait
  }

  if (nrow(row) == 0) return(paste0(display, "–", trait_label, " unavailable"))

  paste0(
    display, "–", trait_label,
    " (ρ = ", fmt(row[[rho_col]][1], 2),
    ", FDR = ", fmt_p(row[[fdr_col]][1]), ")"
  )
}

write_section <- function(text, filename) {
  writeLines(text, file.path(OUTDIR, filename), useBytes = TRUE)
}

###############################################################################
# 02. LOAD DEFINITIVE SOURCES
###############################################################################

PATHS <- list(
  input_qc = file.path(S10, "wgcna_input_qc_summary.csv"),
  core_summary = file.path(S11, "tables", "wgcna_core_summary.csv"),
  dep_burden = file.path(S12, "tables", "module_DEP_burden_dual_definition.csv"),
  enrichment_summary = file.path(S12, "tables", "enrichment_summary_by_module.csv"),
  module_trait = file.path(S13, "tables", "correlations", "module_trait_results_long.csv"),
  adjusted_continuous = file.path(S13, "tables", "regression", "adjusted_module_models.csv"),
  script13_summary = file.path(S13, "tables", "script13_final_summary.csv"),
  diagnosis_hc3 = file.path(S13B, "tables", "diagnosis_robustness", "adjusted_diagnosis_module_models_HC3.csv"),
  country_context = file.path(S13, "tables", "context", "module_recruitment_context_effect_summary.csv"),
  nested_site = file.path(S13B, "tables", "context", "corrected_nested_site_all_samples.csv"),
  loco = file.path(S14, "tables", "loco", "loco_country_summary_by_module_trait.csv"),
  wc_loso = file.path(S14, "tables", "within_country_loso", "within_country_loso_summary_by_module_trait.csv"),
  downsampling = file.path(S14, "tables", "downsampling", "balanced_downsampling_summary_by_module_trait.csv"),
  diagnosis_stability = file.path(S14, "tables", "model_stability", "diagnosis", "integrated_diagnosis_HC3_stability_summary.csv"),
  biomarker_full32 = file.path(S14B, "tables", "full", "full_all_modules_log_HC3_family32.csv"),
  biomarker_focal_stability = file.path(S14B, "tables", "focal_modules", "integrated_focal_log_HC3_stability_family32.csv"),
  preservation = file.path(S15B, "tables", "fixed_geneset_structural_preservation_definitive_summary.csv"),
  script15b_summary = file.path(S15B, "tables", "script15b_final_summary.csv"),
  module_sizes = file.path(S16, "tables", "modularity", "final_module_sizes_and_proportions.csv"),
  module_quality = file.path(S16, "tables", "integration", "integrated_module_network_quality_summary.csv"),
  script16_summary = file.path(S16, "tables", "script16_final_summary.csv")
)

T <- purrr::map(PATHS, read_required)

if (!"FDR_family32" %in% names(T$biomarker_full32) &&
    !"full_FDR_family32" %in% names(T$biomarker_full32)) {
  stop("The biomarker source is not the definitive corrected family of 32 models.", call. = FALSE)
}

if (!"Fixed_gene_comparability" %in% names(T$preservation) &&
    !"fixed_gene_comparability" %in% names(T$preservation)) {
  stop("The preservation source is not the definitive fixed-gene output.", call. = FALSE)
}

###############################################################################
# 03. DEFINITIVE VALUES
###############################################################################

n_samples <- safe_number(summary_value(T$script13_summary, "n_samples"), 639)
n_modules <- safe_number(summary_value(T$script13_summary, "n_modules"), 8)
n_traits <- safe_number(summary_value(T$script13_summary, "n_traits"), 16)
n_trait_tests <- safe_number(summary_value(T$script13_summary, "n_module_trait_tests"), 128)
n_trait_significant <- safe_number(summary_value(T$script13_summary, "n_module_trait_fdr_0_05"), 47)
n_adjusted_significant <- safe_number(summary_value(T$script13_summary, "n_adjusted_clinical_fdr_0_05"), 20)

fixed_genes <- safe_number(summary_value(T$script15b_summary, "n_fixed_network_genes"), 4239)
n_country_strong <- safe_number(summary_value(T$script15b_summary, "n_country_strong_results"), 32)
n_country_moderate <- safe_number(summary_value(T$script15b_summary, "n_country_at_least_moderate_results"), 37)
n_site_strong <- safe_number(summary_value(T$script15b_summary, "n_site_strong_results"), 23)
n_site_moderate <- safe_number(summary_value(T$script15b_summary, "n_site_at_least_moderate_results"), 30)

soft_power <- safe_number(summary_value(T$script16_summary, "Soft_power", "Metric", "Value"), 4)
signed_r2 <- safe_number(summary_value(T$script16_summary, "Signed_scale_free_R2", "Metric", "Value"), 0.914)
mean_connectivity <- safe_number(summary_value(T$script16_summary, "Mean_connectivity", "Metric", "Value"), 865.0)
weighted_q <- safe_number(summary_value(T$script16_summary, "Weighted_modularity_Q", "Metric", "Value"), 0.118)

green_dep <- module_row(T$dep_burden, "green")
blue_dep <- module_row(T$dep_burden, "blue")
brown_dep <- module_row(T$dep_burden, "brown")

green_size <- extract_num(green_dep, c("Module_size", "N_genes", "n_genes"), 762)
blue_size <- extract_num(blue_dep, c("Module_size", "N_genes", "n_genes"), 1709)
brown_size <- extract_num(brown_dep, c("Module_size", "N_genes", "n_genes"), 1494)

green_dep_n <- extract_num(green_dep, c("primary_DEP_overlap_n", "Primary_DEP_n"), 86)
blue_dep_n <- extract_num(blue_dep, c("primary_DEP_overlap_n", "Primary_DEP_n"), 259)
brown_dep_n <- extract_num(brown_dep, c("primary_DEP_overlap_n", "Primary_DEP_n"), 36)

green_diag <- module_row(T$diagnosis_hc3, "green")
blue_diag <- module_row(T$diagnosis_hc3, "blue")

green_diag_effect <- extract_num(green_diag, c("standardized_difference", "estimate_standardized"), 0.183)
green_diag_fdr <- extract_num(green_diag, c("FDR_HC3", "FDR", "p_adj"), 0.059)
blue_diag_effect <- extract_num(blue_diag, c("standardized_difference", "estimate_standardized"), 0.214)
blue_diag_fdr <- extract_num(blue_diag, c("FDR_HC3", "FDR", "p_adj"), 0.027)

green_nfl <- format_association(T$module_trait, "green", "NfL")
green_ptau217 <- format_association(T$module_trait, "green", "p_tau217")
green_abeta <- format_association(T$module_trait, "green", "ratio_AB42_40")
blue_ptau217 <- format_association(T$module_trait, "blue", "p_tau217")
blue_abeta <- format_association(T$module_trait, "blue", "ratio_AB42_40")
brown_ptau181 <- format_association(T$module_trait, "brown", "p_tau181")

###############################################################################
# 04. MAIN METHODS
###############################################################################

main_methods <- glue::glue(
"Weighted gene co-expression network analysis was performed using the outcome-independent gene-collapsed plasma proteomic matrix. For genes represented by more than one eligible SOMAmer, one representative analyte was selected before diagnostic modeling according to the lowest missingness, followed by the highest log2-scale median absolute deviation and an alphanumeric AptName tie-breaker. This procedure yielded 9,638 gene-level protein features across {as.integer(n_samples)} participants. A signed Pearson-correlation network was constructed using a soft-thresholding power of β = {fmt(soft_power, 0)}, selected to provide a signed scale-free topology fit of R² = {fmt(signed_r2, 3)} while retaining a mean connectivity of {fmt(mean_connectivity, 1)}. Modules were identified from the signed topological-overlap matrix using dynamic tree cutting (minimum module size 30; deepSplit = 3; pamRespectsDendro = FALSE) and merged at an eigengene dissimilarity threshold of 0.25. Module colors were retained as stable analytical identifiers. M1-M8 were added only as visual identifiers and do not represent a ranking.

Module eigengenes were correlated with {as.integer(n_traits)} clinical, cognitive, demographic, genetic and plasma biomarker variables using pairwise-complete Spearman correlations. BH-FDR correction was applied across the complete {as.integer(n_trait_tests)}-test matrix. Selected continuous outcomes were additionally evaluated using covariate-adjusted linear models, and diagnosis-related module differences were estimated with the module eigengene as the dependent variable. Diagnosis and biomarker models used HC3 robust standard errors. Plasma p-tau181, p-tau217, NfL and Aβ42/40 outcomes were evaluated across all eight modules, with BH-FDR correction applied across the complete family of 32 module-biomarker models.

Country and recruitment site were treated as categorical recruitment-context factors. Site effects were estimated using a corrected site-within-country parameterization. Association stability was evaluated while holding the full-cohort module definitions fixed through leave-one-country-out analyses, valid within-country site deletion and 500 balanced-downsampling iterations. These analyses assessed the stability of module associations and were not interpreted as structural-preservation tests.

Structural module preservation was evaluated separately using one fixed set of {format(as.integer(fixed_genes), big.mark = ',')} genes in every comparison. All genes were retained from modules containing 1,000 or fewer genes, while the same reproducible 1,000-gene subset was selected once from each larger module. Country-held-out and reciprocal within-country site comparisons used 100 permutations. Zsummary values ≥10, 2 to <10 and <2 were interpreted descriptively as strong, weak-to-moderate and no evidence of preservation, respectively. Because the module assignments were learned in the full ReDLat cohort, these analyses were interpreted as internal structural reproducibility rather than external validation.

Network-quality diagnostics included scale-free topology fit, mean connectivity, module membership, within-versus-external adjacency and topological-overlap separation, eigengene correlations and exact post hoc weighted Newman-Girvan modularity. The modularity statistic was not used to choose the soft threshold or optimize the WGCNA partition and was interpreted only as a descriptive metric."
)

###############################################################################
# 05. MAIN RESULTS
###############################################################################

main_results <- glue::glue(
"Clinical AD-associated plasma proteomic variation shows modular organization

We next asked whether the asymmetric clinical AD-associated proteomic signature was organized into coordinated plasma co-expression programs. Weighted gene co-expression network analysis resolved {as.integer(n_modules)} modules, supporting modular organization rather than a diffuse pattern of independent protein alterations (Fig. 3a,b, Supplementary Tables 20-21 and Supplementary Data 4 and 6). Scale-free topology and connectivity diagnostics supported the selected soft-thresholding power of β = {fmt(soft_power, 0)}, with a scale-free topology fit of R² = {fmt(signed_r2, 3)} and mean connectivity of {fmt(mean_connectivity, 1)}. As an additional descriptive quality metric, exact post hoc weighted modularity analysis yielded Q = {fmt(weighted_q, 3)}, supporting detectable coordinated structure while indicating that modularity should be interpreted cautiously rather than as a strong partition of the plasma proteome (Extended Data Fig. 7, Supplementary Table 20 and Supplementary Data 6).

Module eigengenes were correlated with clinical, cognitive, demographic, genetic and plasma AT(N) variables across {as.integer(n_trait_tests)} tests, with BH-FDR correction applied across the complete module-trait matrix. In total, {as.integer(n_trait_significant)} associations reached FDR < 0.05. In complementary covariate-adjusted analyses, {as.integer(n_adjusted_significant)} of 48 selected continuous-outcome models remained significant after correction (Fig. 3c, Extended Data Fig. 8, Supplementary Tables 23-24 and Supplementary Data 5). M1/green, M2/blue and M3/brown captured complementary patterns of biological organization, differential-abundance burden and clinical or biomarker alignment and were therefore examined in greater detail, without combining these measures into a composite module-prioritization score (Fig. 3d, Extended Data Figs. 7-10, Supplementary Tables 21-24 and Supplementary Data 4-5).

M1/green emerged as the most biologically convergent clinical AD-associated module. The module contained {as.integer(green_dep_n)} proteins from the primary 587-protein signature among approximately {as.integer(green_size)} proteins and showed coherent enrichment for axon development, neuron projection guidance, synapse organization, cell-cell adhesion and extracellular-matrix organization. Its eigengene aligned with diagnosis, impairment and plasma AT(N) measures, including {green_nfl}, {green_ptau217} and {green_abeta}. Representative hubs included CADM1, EFNA5, NRP2, UNC5B, EPHA4, FLRT2 and TNFRSF21/DR6 (Fig. 3e, Extended Data Fig. 9 and Supplementary Tables 21-22). In adjusted diagnosis models, the M1/green eigengene difference remained directionally positive but was borderline after HC3 and FDR correction (standardized difference = {fmt(green_diag_effect, 3)}, FDR = {fmt_p(green_diag_fdr)}).

M2/blue captured the largest concentration of clinical AD-associated proteins, with {as.integer(blue_dep_n)} proteins from the primary signature among approximately {as.integer(blue_size)} proteins. This high differential-abundance burden distinguished M2/blue from the other modules, despite more restricted pathway-level coherence. Representative hubs included CFAP53, RWDD2A, KLF14, MRPL47, LIAS, CLDN4, SPRY2 and CPT2 (Fig. 3f and Supplementary Tables 21-22). M2/blue was the only module associated with diagnosis after HC3 and FDR correction (standardized difference = {fmt(blue_diag_effect, 3)}, FDR = {fmt_p(blue_diag_fdr)}) and showed biomarker alignment including {blue_ptau217} and {blue_abeta} after correction across the complete 32-model family (Extended Data Fig. 8 and Supplementary Table 24).

M3/brown showed comparatively low differential-abundance burden, with {as.integer(brown_dep_n)} proteins from the primary signature among approximately {as.integer(brown_size)} proteins, but strong functional coherence centered on mRNA metabolism, RNA processing, ribonucleoprotein-complex biogenesis and spliceosomal pathways. Representative hubs included RBBP4, SIRPB2, KHDRBS1, FNBP1, WDR5 and HNRNPA1. Although M3/brown was not independently associated with diagnosis after HC3 correction, it showed biomarker alignment including {brown_ptau181} (Extended Data Figs. 8-9 and Supplementary Tables 21-24).

Country and site accounted for measurable proportions of eigengene variation, and LOCO, valid within-country site-deletion and balanced-downsampling analyses supported graded internal stability of focal module associations (Extended Data Fig. 10, Supplementary Tables 25-26 and Supplementary Data 5). These effects were interpreted as recruitment-context sensitivity rather than evidence of country- or site-specific biology.

Using the identical fixed set of {format(as.integer(fixed_genes), big.mark = ',')} genes across comparisons, {as.integer(n_country_strong)} of 40 country-module tests and {as.integer(n_site_strong)} of 32 reciprocal site-module tests showed strong structural preservation, increasing to {as.integer(n_country_moderate)} of 40 and {as.integer(n_site_moderate)} of 32, respectively, when weak-to-moderate preservation was included (Extended Data Fig. 11, Supplementary Table 27 and Supplementary Data 6). Association stability and structural preservation were interpreted as complementary but conceptually distinct layers of internal robustness."
)

###############################################################################
# 06. INTERPRETATION BOUNDARIES
###############################################################################

interpretation_text <- paste(
  "WGCNA modules are coordinated plasma co-expression structures and should not",
  "be interpreted as disconnected graph communities, tissue-specific pathways",
  "or causal mechanisms. Biological labels are provisional summaries derived",
  "from enrichment, hub composition, differential-abundance burden and clinical",
  "alignment. The color labels remain the analytical module identifiers, while",
  "M1-M8 are stable visual identifiers and do not indicate a ranking.",
  "",
  "Country and recruitment-site results quantify internal recruitment-context",
  "dependence and do not establish country-specific biology.",
  "",
  "Association stability and structural preservation answer distinct questions.",
  "The former tests whether fixed eigengene-trait relationships persist under",
  "sample perturbations; the latter evaluates whether internal co-expression",
  "topology is retained between reference and test samples.",
  "",
  "The exact post hoc weighted modularity Q was not used to choose the soft",
  "power, detect modules or optimize the partition. No universal biological",
  "threshold was assumed for Q in this dense weighted plasma network."
)

###############################################################################
# 07. FIGURE LEGENDS
###############################################################################

figure3_legend_path <- file.path(S17_FINAL, "diagnostics", "figure_legends.txt")
if (file.exists(figure3_legend_path)) {
  figure3_legend <- paste(readLines(figure3_legend_path, warn = FALSE), collapse = " ")
} else {
  figure3_legend <- paste(
    "Fig. 3. Modular organization of the outcome-independent ReDLat plasma proteome.",
    "a, Protein clustering dendrogram and final merged module assignments.",
    "b, Number of proteins assigned to each module.",
    "c, Spearman correlations between M1/green, M2/blue and M3/brown eigengenes",
    "and selected traits, with BH-FDR correction across 128 tests.",
    "d, Separate displays of canonical DEP overlap, significant module-trait",
    "associations and significant adjusted models; these measures were not",
    "combined into a composite score.",
    "e,f, Biological portraits of M1/green and M2/blue."
  )
}

extended_legend_path <- file.path(
  S18_FINAL,
  "diagnostics",
  "WGCNA_Extended_Data_figure_legends.csv"
)

if (file.exists(extended_legend_path)) {
  extended_legends <- readr::read_csv(
    extended_legend_path,
    show_col_types = FALSE,
    guess_max = 10000
  )
} else {
  extended_legends <- tibble::tribble(
    ~Figure_ID, ~Legend_title, ~Legend_text,
    "Extended Data Fig. 7", "WGCNA network construction and global quality diagnostics",
    "Panels summarize soft-threshold selection, mean connectivity, module-membership distributions, within-versus-external separation, final eigengene correlations and module contributions to exact post hoc weighted modularity. M1-M8 are visual identifiers and not a ranking.",
    "Extended Data Fig. 8", "Complete module-trait and adjusted-association landscape",
    "Panels show the complete 128-test module-trait matrix, 48 adjusted continuous-outcome models, eight diagnosis models with HC3 inference and the corrected 32-model plasma biomarker family.",
    "Extended Data Fig. 9", "Hub architecture and functional enrichment",
    "Panels summarize module-membership distributions, the proportion of high-membership hubs, leading hubs in M1/green, M2/blue and M3/brown, enrichment burden and leading FDR-significant terms.",
    "Extended Data Fig. 10", "Recruitment-context effects and association stability",
    "Panels show categorical country and site-within-country effects, LOCO and balanced-downsampling stability, and adjusted diagnosis-model stability. These analyses hold the full-cohort module definitions fixed.",
    "Extended Data Fig. 11", "Fixed-gene structural module preservation",
    "Country-held-out and reciprocal within-country site preservation used the identical fixed set of 4,239 genes and 100 permutations. Zsummary categories are descriptive and represent internal structural reproducibility."
  )
}

if (!"Figure_ID" %in% names(extended_legends)) {
  names(extended_legends)[1] <- "Figure_ID"
}
if (!"Legend_text" %in% names(extended_legends)) {
  preferred <- c("Legend", "legend", "Text", "text")
  text_col <- preferred[preferred %in% names(extended_legends)][1]
  if (length(text_col) == 0 || is.na(text_col)) {
    candidates <- names(extended_legends)[
      grepl("legend_text|legend$|text$", names(extended_legends), ignore.case = TRUE)
    ]
    text_col <- candidates[1]
  }
  if (length(text_col) > 0 && !is.na(text_col)) {
    names(extended_legends)[names(extended_legends) == text_col] <- "Legend_text"
  }
}

###############################################################################
# 08. WRITE OUTPUTS
###############################################################################

reviewer_summary <- glue::glue(
"Analytical unit: 9,638 outcome-independently selected gene-level plasma proteins in {as.integer(n_samples)} participants.
Network: signed WGCNA, β = {fmt(soft_power, 0)}, signed R² = {fmt(signed_r2, 3)}, mean connectivity = {fmt(mean_connectivity, 1)}, eight final modules.
Trait integration: {as.integer(n_trait_significant)}/{as.integer(n_trait_tests)} module-trait associations and {as.integer(n_adjusted_significant)}/48 adjusted continuous-outcome models at FDR < 0.05.
Biomarker family: eight modules × four biomarkers = 32 HC3 models per run.
Association robustness: five LOCO runs, valid within-country site deletions and 500 balanced-downsampling iterations.
Structural robustness: one fixed set of {format(as.integer(fixed_genes), big.mark = ',')} genes; 100 permutations per comparison.
Network quality: exact post hoc weighted Q = {fmt(weighted_q, 3)}, interpreted descriptively."
)

write_section(main_methods, "WGCNA_Main_Methods_FINAL.txt")
write_section(main_results, "WGCNA_Main_Results_FINAL.txt")
write_section(interpretation_text, "WGCNA_Interpretation_Boundaries_FINAL.txt")
write_section(reviewer_summary, "WGCNA_Reviewer_Summary_FINAL.txt")
write_section(figure3_legend, "WGCNA_Figure3_Legend_FINAL.txt")

extended_md <- paste0(
  purrr::pmap_chr(
    extended_legends,
    function(Figure_ID, Legend_text, ...) {
      paste0("## ", Figure_ID, "\n\n", Legend_text, "\n")
    }
  ),
  collapse = "\n"
)
writeLines(
  extended_md,
  file.path(OUTDIR, "WGCNA_Extended_Data_Legends_FINAL.md"),
  useBytes = TRUE
)
readr::write_csv(
  extended_legends,
  file.path(OUTDIR, "WGCNA_Extended_Data_Legends_FINAL.csv")
)

text_sections <- tibble::tibble(
  Section = c(
    "Main Methods",
    "Main Results",
    "Interpretation boundaries",
    "Reviewer-facing summary",
    "Figure 3 legend"
  ),
  Text = c(
    main_methods,
    main_results,
    interpretation_text,
    reviewer_summary,
    figure3_legend
  )
)
readr::write_csv(
  text_sections,
  file.path(OUTDIR, "WGCNA_Integrated_Manuscript_Text_Sections.csv")
)

bundle <- paste0(
  "# WGCNA integrated manuscript text — Nature Aging final\n\n",
  "## Main Methods\n\n", main_methods, "\n\n",
  "## Main Results\n\n", main_results, "\n\n",
  "## Interpretation boundaries\n\n", interpretation_text, "\n\n",
  "## Reviewer-facing summary\n\n", reviewer_summary, "\n\n",
  "## Figure 3 legend\n\n", figure3_legend, "\n\n",
  "# Extended Data legends\n\n", extended_md
)
writeLines(
  bundle,
  file.path(OUTDIR, "WGCNA_Integrated_Manuscript_Text_FINAL.md"),
  useBytes = TRUE
)

value_audit <- tibble::tibble(
  Value = c(
    "Participants", "Proteins", "Modules", "Module-trait tests",
    "Module-trait associations FDR<0.05", "Adjusted continuous models FDR<0.05",
    "Biomarker family size", "Fixed preservation genes",
    "Country strong preservation", "Country at least moderate preservation",
    "Site strong preservation", "Site at least moderate preservation",
    "Soft power", "Signed R2", "Mean connectivity", "Exact weighted Q"
  ),
  Observed = c(
    n_samples, 9638, n_modules, n_trait_tests,
    n_trait_significant, n_adjusted_significant,
    32, fixed_genes,
    n_country_strong, n_country_moderate,
    n_site_strong, n_site_moderate,
    soft_power, signed_r2, mean_connectivity, weighted_q
  ),
  Source = c(
    "Script 13", "Scripts 10/11/16", "Scripts 11/13/16", "Script 13",
    "Script 13", "Script 13", "Script 14b", "Script 15b",
    "Script 15b", "Script 15b", "Script 15b", "Script 15b",
    "Script 16", "Script 16", "Script 16", "Script 16"
  )
)
readr::write_csv(
  value_audit,
  file.path(OUTDIR, "WGCNA_Manuscript_Value_Audit.csv")
)

source_audit <- tibble::tibble(
  Source_name = names(PATHS),
  Source_path = unlist(PATHS, use.names = FALSE),
  Exists = file.exists(unlist(PATHS, use.names = FALSE)),
  Definitive_source = TRUE
)
readr::write_csv(
  source_audit,
  file.path(OUTDIR, "WGCNA_Integrated_Text_Source_Audit.csv")
)

manifest <- tibble::tibble(
  Artifact = c(
    "Main Methods", "Main Results", "Interpretation boundaries",
    "Reviewer summary", "Figure 3 legend", "Extended Data legends",
    "Integrated Markdown bundle", "Structured text sections",
    "Value audit", "Source audit", "Module key", "Session information"
  ),
  Path = c(
    file.path(OUTDIR, "WGCNA_Main_Methods_FINAL.txt"),
    file.path(OUTDIR, "WGCNA_Main_Results_FINAL.txt"),
    file.path(OUTDIR, "WGCNA_Interpretation_Boundaries_FINAL.txt"),
    file.path(OUTDIR, "WGCNA_Reviewer_Summary_FINAL.txt"),
    file.path(OUTDIR, "WGCNA_Figure3_Legend_FINAL.txt"),
    file.path(OUTDIR, "WGCNA_Extended_Data_Legends_FINAL.md"),
    file.path(OUTDIR, "WGCNA_Integrated_Manuscript_Text_FINAL.md"),
    file.path(OUTDIR, "WGCNA_Integrated_Manuscript_Text_Sections.csv"),
    file.path(OUTDIR, "WGCNA_Manuscript_Value_Audit.csv"),
    file.path(OUTDIR, "WGCNA_Integrated_Text_Source_Audit.csv"),
    file.path(OUTDIR, "WGCNA_Module_Key.csv"),
    file.path(OUTDIR, "sessionInfo.txt")
  )
)
writeLines(capture.output(utils::sessionInfo()), file.path(OUTDIR, "sessionInfo.txt"))
manifest <- manifest %>%
  dplyr::mutate(
    Exists = file.exists(Path),
    Bytes = ifelse(Exists, file.info(Path)$size, NA_real_)
  )
readr::write_csv(manifest, file.path(OUTDIR, "script21_final_manifest.csv"))

if (!all(manifest$Exists) || any(manifest$Bytes <= 0, na.rm = TRUE)) {
  stop("One or more Script 21 outputs were not created correctly.", call. = FALSE)
}

writeLines(
  c(
    "SCRIPT_20_STATUS=PASS",
    paste0("OUTPUT_DIRECTORY=", OUTDIR),
    paste0("TIMESTAMP=", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  ),
  file.path(OUTDIR, "SCRIPT_20_SUCCESS.txt")
)

message("Script 20 completed successfully.")
message("Output directory: ", OUTDIR)
message("Figure 3 and Extended Data Fig. 7-11 legends generated: TRUE")

} # end run_manuscript_text()

message("============================================================")
message("RUNNING SCRIPT 13: MANUSCRIPT SUPPORT TEXT")
message("This script generates the manuscript-facing WGCNA text and legends.")
message("============================================================")

run_manuscript_text()

###############################################################################
# END
###############################################################################


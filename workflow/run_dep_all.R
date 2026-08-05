# ReDLat DEP workflow — complete runner

args <- commandArgs(trailingOnly = TRUE)

project_root <- Sys.getenv(
  "REDLAT_PROJECT_ROOT",
  unset = normalizePath(getwd(), winslash = "/", mustWork = TRUE)
)

rscript <- file.path(R.home("bin"), "Rscript")

run_optional_simulation <- "--with-simulation" %in% args

analysis_scripts <- c(
  "scripts/DEP/01_DEP_primary_analysis.R",
  "scripts/DEP/02_DEP_rebuild_sensitivity_fixed_map.R",
  "scripts/DEP/03_DEP_APOE_ATN_equal_sample_models.R",
  "scripts/DEP/04_DEP_country_site_robustness.R",
  "scripts/DEP/05_DEP_matching_biomarker_sensitivity.R",
  "scripts/DEP/06_DEP_ptau217_threshold_sweep.R",
  "scripts/DEP/07_DEP_demographic_balance.R"
)

if (run_optional_simulation) {
  analysis_scripts <- c(
    analysis_scripts,
    "scripts/DEP/13_DEP_simulation_based_detectable_effect_analysis.R"
  )
} else {
  message(
    "Skipping optional simulation-based sensitivity analysis. ",
    "Use --with-simulation to run script 13."
  )
}

reporting_scripts <- c(
  "scripts/DEP/08_DEP_generate_supplementary_tables.R",
  "scripts/DEP/09_DEP_generate_source_data.R",
  "scripts/DEP/10_DEP_generate_supplementary_data.R",
  "scripts/DEP/11_DEP_generate_main_figure2.R",
  "scripts/DEP/12_DEP_generate_extended_data_figures.R"
)

scripts <- c(
  analysis_scripts,
  reporting_scripts
)

for (script in scripts) {
  path <- file.path(project_root, script)
  
  if (!file.exists(path)) {
    stop(
      "Required script not found: ",
      path,
      call. = FALSE
    )
  }
  
  message("\n=== Running ", script, " ===")
  
  status <- system2(
    rscript,
    c("--vanilla", shQuote(path)),
    env = paste0(
      "REDLAT_PROJECT_ROOT=",
      shQuote(project_root)
    )
  )
  
  if (!identical(status, 0L)) {
    stop(
      "Pipeline stopped after failure in ",
      script,
      call. = FALSE
    )
  }
}

message("Completed: Complete DEP pipeline")
# Run the WGCNA analysis workflow
rm(list = ls())
options(stringsAsFactors = FALSE)

project_root <- if (nzchar(Sys.getenv("REDLAT_PROJECT_ROOT", unset = ""))) {
  normalizePath(Sys.getenv("REDLAT_PROJECT_ROOT"), winslash = "/", mustWork = TRUE)
} else if (requireNamespace("here", quietly = TRUE)) {
  normalizePath(here::here(), winslash = "/", mustWork = TRUE)
} else {
  stop("Package 'here' is required. Run renv::restore().", call. = FALSE)
}

scripts <- c(
  "01_WGCNA_prepare_input.R",
  "02_WGCNA_construct_network.R",
  "03_WGCNA_characterize_modules.R",
  "04_WGCNA_module_trait_models.R",
  "05_WGCNA_site_biomarker_sensitivity.R",
  "06_WGCNA_association_stability.R",
  "07_WGCNA_biomarker_fdr_correction.R",
  "08_WGCNA_structural_preservation.R",
  "09_WGCNA_network_quality.R"
)
log_dir <- file.path(project_root, "result", "WGCNA", "workflow_logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

for (script in scripts) {
  path <- file.path(project_root, "scripts", "WGCNA", script)
  if (!file.exists(path)) stop("Missing workflow script: ", path, call. = FALSE)
  message("Running ", script)
  started <- Sys.time()
  tryCatch(
    sys.source(path, envir = new.env(parent = globalenv()), encoding = "UTF-8"),
    error = function(e) {
      writeLines(c(
        paste0("SCRIPT=", script),
        paste0("STARTED=", started),
        paste0("FAILED=", Sys.time()),
        paste0("ERROR=", conditionMessage(e))
      ), file.path(log_dir, paste0(script, ".failed.txt")))
      stop("WGCNA workflow stopped in ", script, ": ", conditionMessage(e), call. = FALSE)
    }
  )
  writeLines(c(
    paste0("SCRIPT=", script),
    paste0("STARTED=", started),
    paste0("COMPLETED=", Sys.time())
  ), file.path(log_dir, paste0(script, ".completed.txt")))
}

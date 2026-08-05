# ReDLat DEP workflow — reporting runner
args <- commandArgs(trailingOnly = TRUE)
project_root <- Sys.getenv("REDLAT_PROJECT_ROOT", unset = normalizePath(getwd(), winslash = "/", mustWork = TRUE))
rscript <- file.path(R.home("bin"), "Rscript")
scripts <- c(
  "scripts/DEP/08_DEP_generate_supplementary_tables.R",
  "scripts/DEP/09_DEP_generate_source_data.R",
  "scripts/DEP/10_DEP_generate_supplementary_data.R",
  "scripts/DEP/11_DEP_generate_main_figure2.R",
  "scripts/DEP/12_DEP_generate_extended_data_figures.R"
)
for (script in scripts) {
  path <- file.path(project_root, script)
  message("\n=== Running ", script, " ===")
  status <- system2(rscript, c("--vanilla", shQuote(path)), env = paste0("REDLAT_PROJECT_ROOT=", shQuote(project_root)))
  if (!identical(status, 0L)) stop("Pipeline stopped after failure in ", script, call. = FALSE)
}
message("Completed: DEP reporting pipeline")

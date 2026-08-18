# Project-relative configuration only; no personal absolute paths.
if (file.exists(file.path(getwd(), ".redlat-root"))) {
  Sys.setenv(
    REDLAT_PROJECT_ROOT = normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  )
}
options(stringsAsFactors = FALSE)

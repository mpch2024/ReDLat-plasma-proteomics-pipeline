# Audit the WGCNA publication candidate for identifiers and private paths.
rm(list = ls())
project_root <- if (nzchar(Sys.getenv("REDLAT_PROJECT_ROOT", unset = ""))) normalizePath(Sys.getenv("REDLAT_PROJECT_ROOT"), winslash = "/", mustWork = TRUE) else here::here()
source(file.path(project_root, "R", "wgcna_bootstrap.R"), local = FALSE)
config <- wgcna_load_config(project_root)
wgcna_assert_public_tree(config$publication_root)
message("WGCNA publication-candidate privacy audit passed.")

###############################################################################
# ReDLat plasma proteomics — DEP workflow
# 13. Simulation-based sensitivity and detectable-effect analysis
# Requires: Complete-case metadata, quality-controlled SOMAscan data and DEP inputs
# Produces: Empirical power, FDR and detectable-effect summaries
# Interpretation: Retrospective sensitivity analysis; not an a priori sample-size
# calculation or observed post hoc power.
# Data policy: participant-level data and intermediate outputs remain local.
###############################################################################

rm(list = ls())
options(stringsAsFactors = FALSE, warn = 1)

# -----------------------------------------------------------------------------
# Project setup
# -----------------------------------------------------------------------------
.project_root_env <- Sys.getenv("REDLAT_PROJECT_ROOT", unset = "")
if (nzchar(.project_root_env)) {
  project_root <- normalizePath(.project_root_env, winslash = "/", mustWork = TRUE)
} else if (requireNamespace("here", quietly = TRUE)) {
  project_root <- normalizePath(here::here(), winslash = "/", mustWork = TRUE)
} else {
  stop(
    "Package 'here' is required. Restore the project environment with renv::restore().",
    call. = FALSE
  )
}

source(file.path(project_root, "R", "dep_bootstrap.R"), local = FALSE)
DEP_CONFIG <- dep_load_config(project_root)
analysis_root <- DEP_CONFIG$result_root
publication_root <- DEP_CONFIG$publication_root

PROJECT_ROOT <- project_root
WORKSPACE_FILE <- file.path(
  analysis_root,
  "workspace",
  "analysis_workspace.RData"
)

###############################################################################
# 01. USER SETTINGS
###############################################################################

# The canonical DEP workspace used in the existing ReDLat pipeline stores the
# protein columns before the model-level log2 transformation. Keep FALSE unless
# you confirm that dep_df already contains log2 abundances.
EXPRESSION_IS_LOG2 <- FALSE

# Exact settings used by the primary manuscript pipeline.
FDR_TARGET <- 0.05
EBAYES_TREND <- FALSE
EBAYES_ROBUST <- FALSE

# Execution profile. Use "validation" to confirm that the script runs and that
# the null model is calibrated. Use "final" for manuscript-quality estimates.
RUN_PROFILE <- "final"
if (!RUN_PROFILE %in% c("validation", "final")) {
  stop("RUN_PROFILE must be 'validation' or 'final'.", call. = FALSE)
}

# Standardized diagnosis effects, measured in each analyte's empirical residual
# SD. The fine grid around 0.30-0.42 follows the preliminary ReDLat run, in which
# 0.40 was the smallest evaluated effect reaching >=80% average power.
EFFECT_SIZES_SD_PRIMARY <- c(
  0, 0.15, 0.20, 0.25, 0.30,
  0.32, 0.34, 0.36, 0.38, 0.40, 0.42,
  0.45, 0.50
)

# Focused bootstrap grid used as a sensitivity analysis around the expected
# transition to 80% average power. It is intentionally smaller to limit
# computational burden while directly testing the threshold region.
EFFECT_SIZES_SD_BOOTSTRAP <- c(
  0, 0.30, 0.32, 0.34, 0.36, 0.38, 0.40, 0.42, 0.45
)

# Assumed fractions of truly associated analytes. These are prespecified
# sensitivity scenarios and are not estimated from the 587 significant genes.
NON_NULL_PROPORTIONS <- c(0.01, 0.05, 0.10)

# Primary residual permutation and secondary residual bootstrap.
RUN_BOOTSTRAP_SENSITIVITY <- TRUE
N_SIM_PERMUTATION <- if (RUN_PROFILE == "final") 500L else 100L
N_SIM_BOOTSTRAP <- if (RUN_PROFILE == "final") 250L else 50L

# By default, evaluate only the observed complete-case cohort because the
# manuscript question concerns the operating characteristics supported by the
# available sample. Set TRUE only for a separate prospective planning analysis.
EVALUATE_SAMPLE_SIZE_CURVE <- FALSE
SAMPLE_SIZES <- NULL

TARGET_AVERAGE_POWER <- 0.80
SEED_GLOBAL <- 20260805L
SAVE_ITERATION_LEVEL_RESULTS <- TRUE

###############################################################################
# 02. PACKAGES
###############################################################################

dep_require_packages("limma")

###############################################################################
# 03. HELPERS
###############################################################################

first_existing_file <- function(paths) {
  paths <- unique(paths[!is.na(paths) & nzchar(paths)])
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0L) return(NA_character_)
  normalizePath(hit[[1]], winslash = "/", mustWork = TRUE)
}

find_workspace <- function(root, configured_analysis_root = NULL) {
  direct <- first_existing_file(c(
    if (!is.null(configured_analysis_root)) {
      file.path(configured_analysis_root, "workspace", "analysis_workspace.RData")
    },
    file.path(root, "result", "workspace", "analysis_workspace.RData"),
    file.path(root, "results", "DEP", "workspace", "analysis_workspace.RData"),
    file.path(root, "results", "DEP", "result", "workspace", "analysis_workspace.RData"),
    file.path(root, "workspace", "analysis_workspace.RData")
  ))
  if (!is.na(direct)) return(direct)

  hits <- list.files(
    root,
    pattern = "^analysis_workspace\\.RData$",
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(hits) == 0L) {
    stop(
      "Could not find analysis_workspace.RData under:\n", root,
      "\nSet WORKSPACE_FILE manually in section 01.",
      call. = FALSE
    )
  }
  normalizePath(hits[[1]], winslash = "/", mustWork = TRUE)
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

safe_log2_matrix <- function(mat) {
  storage.mode(mat) <- "double"
  if (any(!is.finite(mat))) {
    stop("The expression matrix contains non-finite values before log2.", call. = FALSE)
  }
  if (any(mat <= 0)) {
    stop(
      "Non-positive abundance values were found. Review preprocessing before log2.",
      call. = FALSE
    )
  }
  log2(mat)
}

clamp <- function(x, lo = 0, hi = 1) {
  pmax(lo, pmin(hi, x))
}

mean_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else mean(x)
}

sd_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1L) NA_real_ else stats::sd(x)
}

mc_summary <- function(x, bounded_01 = FALSE) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    return(c(mean = NA_real_, mc_se = NA_real_, lo = NA_real_, hi = NA_real_))
  }
  m <- mean(x)
  se <- if (length(x) > 1L) stats::sd(x) / sqrt(length(x)) else NA_real_
  lo <- if (is.finite(se)) m - 1.96 * se else NA_real_
  hi <- if (is.finite(se)) m + 1.96 * se else NA_real_
  if (bounded_01) {
    lo <- clamp(lo)
    hi <- clamp(hi)
  }
  c(mean = m, mc_se = se, lo = lo, hi = hi)
}

allocate_stratified_counts <- function(strata, n_target) {
  strata <- droplevels(factor(strata))
  available <- table(strata)
  n_available <- sum(available)

  if (n_target > n_available) {
    stop("Candidate sample size exceeds observed complete-case sample.", call. = FALSE)
  }

  raw <- n_target * as.numeric(available) / n_available
  allocated <- floor(raw)
  names(allocated) <- names(available)

  # Preserve every observed country-by-diagnosis stratum when feasible.
  if (n_target >= length(available)) {
    zero_cells <- which(allocated == 0L & available > 0L)
    if (length(zero_cells) > 0L) allocated[zero_cells] <- 1L
  }

  # If the minimum-one adjustment overshoots, remove from the largest cells.
  while (sum(allocated) > n_target) {
    removable <- which(allocated > 1L)
    if (length(removable) == 0L) removable <- which(allocated > 0L)
    j <- removable[which.max(allocated[removable] - raw[removable])]
    allocated[j] <- allocated[j] - 1L
  }

  fractional <- raw - floor(raw)
  while (sum(allocated) < n_target) {
    eligible <- which(allocated < as.numeric(available))
    if (length(eligible) == 0L) {
      stop("Unable to allocate the requested stratified sample.", call. = FALSE)
    }
    score <- fractional[eligible] + runif(length(eligible), 0, 1e-9)
    j <- eligible[which.max(score)]
    allocated[j] <- allocated[j] + 1L
    fractional[j] <- -Inf
  }

  allocated
}

draw_stratified_indices <- function(strata, n_target) {
  strata <- droplevels(factor(strata))
  allocated <- allocate_stratified_counts(strata, n_target)
  out <- integer(0)

  for (s in names(allocated)) {
    pool <- which(as.character(strata) == s)
    k <- allocated[[s]]
    if (k > 0L) {
      out <- c(out, sample(pool, size = k, replace = FALSE))
    }
  }

  sample(out, length(out), replace = FALSE)
}

draw_residual_donors_by_country <- function(
    full_country,
    selected_country,
    resampling_method = c("permutation", "bootstrap")
) {
  resampling_method <- match.arg(resampling_method)
  full_country <- as.character(full_country)
  selected_country <- as.character(selected_country)
  donors <- integer(length(selected_country))
  use_replacement <- identical(resampling_method, "bootstrap")

  for (country_i in unique(selected_country)) {
    positions <- which(selected_country == country_i)
    pool <- which(full_country == country_i)

    if (!use_replacement && length(pool) < length(positions)) {
      stop(
        "Residual donor pool is smaller than the selected sample in country: ",
        country_i, call. = FALSE
      )
    }

    donors[positions] <- sample(
      pool,
      size = length(positions),
      replace = use_replacement
    )
  }

  donors
}

make_design <- function(meta) {
  meta <- as.data.frame(meta)
  meta$SampleGroup <- factor(as.character(meta$SampleGroup), levels = c("CN", "AD"))
  meta$Sex <- droplevels(factor(as.character(meta$Sex)))
  meta$Country <- droplevels(factor(as.character(meta$Country)))
  meta$Age <- safe_numeric(meta$Age)
  meta$Education <- safe_numeric(meta$Education)

  design <- stats::model.matrix(
    ~ SampleGroup + Age + Sex + Country + Education,
    data = meta
  )

  if (qr(design)$rank != ncol(design)) {
    stop("The simulated design matrix is rank deficient.", call. = FALSE)
  }
  if (!"SampleGroupAD" %in% colnames(design)) {
    stop("SampleGroupAD coefficient was not found.", call. = FALSE)
  }
  design
}

power_two_sided_nct <- function(effect_sd, alpha, df, variance_multiplier) {
  crit <- stats::qt(1 - alpha / 2, df = df)
  ncp <- effect_sd / sqrt(variance_multiplier)
  stats::pt(-crit, df = df, ncp = ncp) +
    (1 - stats::pt(crit, df = df, ncp = ncp))
}

solve_detectable_effect <- function(alpha, target_power, df, variance_multiplier) {
  f <- function(d) {
    power_two_sided_nct(d, alpha, df, variance_multiplier) - target_power
  }
  upper <- 0.25
  while (f(upper) < 0 && upper < 10) upper <- upper * 2
  if (upper >= 10 && f(upper) < 0) return(NA_real_)
  stats::uniroot(f, lower = 0, upper = upper, tol = 1e-8)$root
}

write_session_info <- function(path) {
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(capture.output(sessionInfo()), con)
}

###############################################################################
# 04. LOAD THE CANONICAL ReDLat WORKSPACE
###############################################################################

if (!file.exists(WORKSPACE_FILE)) {
  WORKSPACE_FILE <- find_workspace(
    PROJECT_ROOT,
    configured_analysis_root = analysis_root
  )
} else {
  WORKSPACE_FILE <- normalizePath(
    WORKSPACE_FILE,
    winslash = "/",
    mustWork = TRUE
  )
}

message("Using workspace: ", WORKSPACE_FILE)

ws <- new.env(parent = baseenv())
load(WORKSPACE_FILE, envir = ws)

required_objects <- c("dep_df", "dep_protein_cols")
missing_objects <- required_objects[
  !vapply(required_objects, exists, logical(1), envir = ws, inherits = FALSE)
]
if (length(missing_objects) > 0L) {
  stop(
    "Workspace is missing required objects: ",
    paste(missing_objects, collapse = ", "),
    call. = FALSE
  )
}

dep_df <- as.data.frame(get("dep_df", envir = ws))
protein_cols <- as.character(get("dep_protein_cols", envir = ws))
protein_cols <- intersect(protein_cols, names(dep_df))

required_metadata <- c("SampleGroup", "Age", "Sex", "Country", "Education")
missing_metadata <- setdiff(required_metadata, names(dep_df))
if (length(missing_metadata) > 0L) {
  stop(
    "dep_df is missing required metadata: ",
    paste(missing_metadata, collapse = ", "),
    call. = FALSE
  )
}
if (length(protein_cols) == 0L) {
  stop("No dep_protein_cols were found in dep_df.", call. = FALSE)
}

meta0 <- dep_df[, required_metadata, drop = FALSE]
meta0$SampleGroup <- factor(as.character(meta0$SampleGroup), levels = c("CN", "AD"))
meta0$Sex <- factor(as.character(meta0$Sex))
meta0$Country <- factor(as.character(meta0$Country))
meta0$Age <- safe_numeric(meta0$Age)
meta0$Education <- safe_numeric(meta0$Education)

keep <- stats::complete.cases(meta0) &
  meta0$SampleGroup %in% c("CN", "AD")

meta <- droplevels(meta0[keep, , drop = FALSE])
expr_samples_by_features <- dep_df[keep, protein_cols, drop = FALSE]
expr <- t(as.matrix(data.frame(
  lapply(expr_samples_by_features, safe_numeric),
  check.names = FALSE
)))
rownames(expr) <- protein_cols

if (!EXPRESSION_IS_LOG2) {
  expr <- safe_log2_matrix(expr)
} else {
  storage.mode(expr) <- "double"
}

if (any(!is.finite(expr))) {
  stop(
    "The complete-case expression matrix contains missing/non-finite values. ",
    "Use the same filtering/imputation rules as the primary DEP pipeline.",
    call. = FALSE
  )
}

n_observed <- ncol(expr)
m_tests <- nrow(expr)

if (is.null(SAMPLE_SIZES)) {
  SAMPLE_SIZES <- if (isTRUE(EVALUATE_SAMPLE_SIZE_CURVE)) {
    unique(as.integer(round(c(
      max(100, 0.40 * n_observed),
      0.60 * n_observed,
      0.80 * n_observed,
      n_observed
    ))))
  } else {
    n_observed
  }
}
SAMPLE_SIZES <- sort(unique(SAMPLE_SIZES))
n_observed_strata <- nlevels(droplevels(interaction(
  meta$SampleGroup,
  meta$Country,
  drop = TRUE
)))
SAMPLE_SIZES <- SAMPLE_SIZES[
  SAMPLE_SIZES >= max(30L, n_observed_strata) &
    SAMPLE_SIZES <= n_observed
]
if (length(SAMPLE_SIZES) == 0L) {
  stop("No valid SAMPLE_SIZES remain after checking the observed sample.", call. = FALSE)
}

message("Complete-case N: ", n_observed)
message("Tested analytes: ", m_tests)
message("Candidate N: ", paste(SAMPLE_SIZES, collapse = ", "))

###############################################################################
# 05. EMPIRICAL MODEL USED TO GENERATE SIMULATED DATA
###############################################################################

design_full <- make_design(meta)

fit_raw <- limma::lmFit(expr, design_full)

if (any(!is.finite(fit_raw$coefficients))) {
  stop("Non-finite coefficients found in the empirical limma fit.", call. = FALSE)
}

fitted_full <- fit_raw$coefficients %*% t(design_full)
residuals_full <- expr - fitted_full
residuals_full <- sweep(
  residuals_full,
  1,
  rowMeans(residuals_full),
  FUN = "-"
)

# Use the feature-specific residual SD estimated by the same linear model.
residual_sd <- as.numeric(fit_raw$sigma)
if (any(!is.finite(residual_sd) | residual_sd <= 0)) {
  stop("Invalid empirical residual SD for one or more analytes.", call. = FALSE)
}

# Counterfactual null mean: retain fitted demographic/country effects while
# setting the diagnosis coefficient to zero.
design_null <- design_full
design_null[, "SampleGroupAD"] <- 0
mu_null <- fit_raw$coefficients %*% t(design_null)

# Residual vectors are centered within country. Resampling entire participant
# columns preserves cross-protein dependence while avoiding diagnosis leakage.
residuals_country_centered <- residuals_full
for (country_i in levels(meta$Country)) {
  idx_country <- which(meta$Country == country_i)
  residuals_country_centered[, idx_country] <- sweep(
    residuals_full[, idx_country, drop = FALSE],
    1,
    rowMeans(residuals_full[, idx_country, drop = FALSE]),
    FUN = "-"
  )
}

###############################################################################
# 06. ANALYTICAL CROSS-CHECK FOR THE OBSERVED DESIGN
###############################################################################

xtx_inv <- solve(crossprod(design_full))
variance_multiplier <- xtx_inv["SampleGroupAD", "SampleGroupAD"]
df_residual <- nrow(design_full) - ncol(design_full)

analytical_crosscheck <- data.frame(
  method = c(
    "Nominal two-sided alpha 0.05",
    "Bonferroni across tested analytes"
  ),
  alpha = c(FDR_TARGET, FDR_TARGET / m_tests),
  target_power = TARGET_AVERAGE_POWER,
  detectable_effect_residual_SD = c(
    solve_detectable_effect(
      FDR_TARGET,
      TARGET_AVERAGE_POWER,
      df_residual,
      variance_multiplier
    ),
    solve_detectable_effect(
      FDR_TARGET / m_tests,
      TARGET_AVERAGE_POWER,
      df_residual,
      variance_multiplier
    )
  ),
  stringsAsFactors = FALSE
)

###############################################################################
# 07. SCENARIO GRID
###############################################################################

make_scenario_grid <- function(
    sample_sizes,
    effect_sizes,
    non_null_proportions,
    resampling_method,
    n_sim_requested
) {
  positive_effects <- sort(unique(effect_sizes[effect_sizes > 0]))

  nonnull <- expand.grid(
    sample_size = sample_sizes,
    effect_sd = positive_effects,
    pi1 = non_null_proportions,
    stringsAsFactors = FALSE
  )

  null <- data.frame(
    sample_size = sample_sizes,
    effect_sd = 0,
    pi1 = 0,
    stringsAsFactors = FALSE
  )

  out <- rbind(null, nonnull)
  out$resampling_method <- resampling_method
  out$n_sim_requested <- as.integer(n_sim_requested)
  out
}

scenario_grid <- make_scenario_grid(
  sample_sizes = SAMPLE_SIZES,
  effect_sizes = EFFECT_SIZES_SD_PRIMARY,
  non_null_proportions = NON_NULL_PROPORTIONS,
  resampling_method = "permutation",
  n_sim_requested = N_SIM_PERMUTATION
)

if (isTRUE(RUN_BOOTSTRAP_SENSITIVITY)) {
  scenario_grid <- rbind(
    scenario_grid,
    make_scenario_grid(
      sample_sizes = SAMPLE_SIZES,
      effect_sizes = EFFECT_SIZES_SD_BOOTSTRAP,
      non_null_proportions = NON_NULL_PROPORTIONS,
      resampling_method = "bootstrap",
      n_sim_requested = N_SIM_BOOTSTRAP
    )
  )
}

scenario_grid$scenario_id <- seq_len(nrow(scenario_grid))
scenario_grid <- scenario_grid[
  order(
    scenario_grid$resampling_method,
    scenario_grid$sample_size,
    scenario_grid$pi1,
    scenario_grid$effect_sd
  ),
]

###############################################################################
# 08. ONE SIMULATION
###############################################################################

run_one_simulation <- function(
    sample_size, effect_sd, pi1, resampling_method, iteration, scenario_id
) {
  set.seed(SEED_GLOBAL + scenario_id * 100000L + iteration)

  strata <- interaction(
    meta$Country,
    meta$SampleGroup,
    drop = TRUE,
    lex.order = TRUE
  )
  selected <- draw_stratified_indices(strata, sample_size)
  meta_sim <- droplevels(meta[selected, , drop = FALSE])

  design_sim <- make_design(meta_sim)
  mu_sim <- mu_null[, selected, drop = FALSE]

  residual_donors <- draw_residual_donors_by_country(
    full_country = meta$Country,
    selected_country = meta_sim$Country,
    resampling_method = resampling_method
  )
  y_sim <- mu_sim + residuals_country_centered[, residual_donors, drop = FALSE]

  n_nonnull <- if (effect_sd > 0 && pi1 > 0) {
    max(1L, as.integer(round(pi1 * m_tests)))
  } else {
    0L
  }

  truth <- rep(FALSE, m_tests)
  median_abs_log2fc <- 0

  if (n_nonnull > 0L) {
    signal_idx <- sample.int(m_tests, size = n_nonnull, replace = FALSE)
    truth[signal_idx] <- TRUE

    directions <- sample(c(-1, 1), size = n_nonnull, replace = TRUE)
    beta_signal <- effect_sd * residual_sd[signal_idx] * directions
    ad_indicator <- as.numeric(meta_sim$SampleGroup == "AD")

    y_sim[signal_idx, ] <- y_sim[signal_idx, , drop = FALSE] +
      beta_signal %o% ad_indicator

    median_abs_log2fc <- stats::median(abs(beta_signal))
  }

  fit_sim <- limma::eBayes(
    limma::lmFit(y_sim, design_sim),
    trend = EBAYES_TREND,
    robust = EBAYES_ROBUST
  )

  p_values <- fit_sim$p.value[, "SampleGroupAD"]
  q_values <- stats::p.adjust(p_values, method = "BH")
  discovered <- is.finite(q_values) & q_values < FDR_TARGET

  tp <- sum(discovered & truth)
  fp <- sum(discovered & !truth)
  total_discoveries <- sum(discovered)

  data.frame(
    scenario_id = scenario_id,
    iteration = iteration,
    resampling_method = resampling_method,
    sample_size = sample_size,
    effect_sd = effect_sd,
    pi1 = pi1,
    n_nonnull = n_nonnull,
    median_abs_log2fc = median_abs_log2fc,
    median_fold_change = 2^median_abs_log2fc,
    median_percent_difference = (2^median_abs_log2fc - 1) * 100,
    true_positives = tp,
    false_positives = fp,
    total_discoveries = total_discoveries,
    average_power = if (n_nonnull > 0L) tp / n_nonnull else NA_real_,
    false_discovery_proportion = if (total_discoveries > 0L) {
      fp / total_discoveries
    } else {
      0
    },
    any_discovery = total_discoveries > 0L,
    stringsAsFactors = FALSE
  )
}

###############################################################################
# 09. RUN SIMULATIONS
###############################################################################

set.seed(SEED_GLOBAL)
iteration_results <- vector("list", nrow(scenario_grid))

for (i in seq_len(nrow(scenario_grid))) {
  sc <- scenario_grid[i, ]
  message(
    "[", i, "/", nrow(scenario_grid), "] method=", sc$resampling_method,
    "; N=", sc$sample_size,
    "; effect SD=", sc$effect_sd,
    "; pi1=", sc$pi1,
    "; simulations=", sc$n_sim_requested
  )

  reps <- vector("list", sc$n_sim_requested)
  for (b in seq_len(sc$n_sim_requested)) {
    reps[[b]] <- run_one_simulation(
      sample_size = sc$sample_size,
      effect_sd = sc$effect_sd,
      pi1 = sc$pi1,
      resampling_method = sc$resampling_method,
      iteration = b,
      scenario_id = sc$scenario_id
    )
  }
  iteration_results[[i]] <- do.call(rbind, reps)
}

iteration_results <- do.call(rbind, iteration_results)

###############################################################################
# 10. SUMMARIZE OPERATING CHARACTERISTICS
###############################################################################

split_key <- interaction(
  iteration_results$resampling_method,
  iteration_results$sample_size,
  iteration_results$effect_sd,
  iteration_results$pi1,
  drop = TRUE,
  lex.order = TRUE
)

summary_list <- lapply(split(iteration_results, split_key), function(x) {
  pwr <- mc_summary(x$average_power, bounded_01 = TRUE)
  fdr <- mc_summary(x$false_discovery_proportion, bounded_01 = TRUE)
  disc <- mc_summary(x$total_discoveries, bounded_01 = FALSE)

  data.frame(
    resampling_method = x$resampling_method[[1]],
    sample_size = x$sample_size[[1]],
    effect_sd = x$effect_sd[[1]],
    pi1 = x$pi1[[1]],
    n_sim = nrow(x),
    median_abs_log2fc = stats::median(x$median_abs_log2fc),
    median_fold_change = stats::median(x$median_fold_change),
    median_percent_difference = stats::median(x$median_percent_difference),
    mean_average_power = pwr[["mean"]],
    power_mc_se = pwr[["mc_se"]],
    power_mc95_low = pwr[["lo"]],
    power_mc95_high = pwr[["hi"]],
    mean_fdr = fdr[["mean"]],
    fdr_mc_se = fdr[["mc_se"]],
    fdr_mc95_low = fdr[["lo"]],
    fdr_mc95_high = fdr[["hi"]],
    probability_fdp_le_target = mean(
      x$false_discovery_proportion <= FDR_TARGET
    ),
    mean_discoveries = disc[["mean"]],
    median_discoveries = stats::median(x$total_discoveries),
    probability_any_discovery = mean(x$any_discovery),
    stringsAsFactors = FALSE
  )
})

simulation_summary <- do.call(rbind, summary_list)
simulation_summary <- simulation_summary[
  order(
    simulation_summary$resampling_method,
    simulation_summary$sample_size,
    simulation_summary$pi1,
    simulation_summary$effect_sd
  ),
]

# Minimum EVALUATED effect reaching the prespecified operating target.
nonnull_summary <- simulation_summary[
  simulation_summary$effect_sd > 0 &
    simulation_summary$pi1 > 0,
]

mde_list <- lapply(
  split(
    nonnull_summary,
    interaction(
      nonnull_summary$resampling_method,
      nonnull_summary$sample_size,
      nonnull_summary$pi1,
      drop = TRUE,
      lex.order = TRUE
    )
  ),
  function(x) {
    x <- x[order(x$effect_sd), ]
    eligible <- x[
      is.finite(x$mean_average_power) &
        x$mean_average_power >= TARGET_AVERAGE_POWER &
        is.finite(x$mean_fdr) &
        x$mean_fdr <= FDR_TARGET,
    ]

    data.frame(
      resampling_method = x$resampling_method[[1]],
      sample_size = x$sample_size[[1]],
      pi1 = x$pi1[[1]],
      minimum_evaluated_effect_sd_for_target = if (nrow(eligible) > 0L) {
        eligible$effect_sd[[1]]
      } else {
        NA_real_
      },
      corresponding_median_abs_log2fc = if (nrow(eligible) > 0L) {
        eligible$median_abs_log2fc[[1]]
      } else {
        NA_real_
      },
      corresponding_median_fold_change = if (nrow(eligible) > 0L) {
        eligible$median_fold_change[[1]]
      } else {
        NA_real_
      },
      corresponding_median_percent_difference = if (nrow(eligible) > 0L) {
        eligible$median_percent_difference[[1]]
      } else {
        NA_real_
      },
      target_average_power = TARGET_AVERAGE_POWER,
      target_fdr = FDR_TARGET,
      note = if (nrow(eligible) > 0L) {
        paste(
          "Smallest effect on the prespecified grid meeting both targets;",
          "this is a grid-based threshold, not an exact continuous MDE."
        )
      } else {
        "No evaluated effect met both targets; extend the effect-size grid."
      },
      stringsAsFactors = FALSE
    )
  }
)

mde_table <- do.call(rbind, mde_list)
mde_table <- mde_table[
  order(
    mde_table$resampling_method,
    mde_table$sample_size,
    mde_table$pi1
  ),
]

# Exploratory linear interpolation of the 80% power crossing. This output is
# useful for planning the next simulation grid but should not be the principal
# reported result.
interpolate_power_crossing <- function(x, target_power = 0.80) {
  x <- x[
    is.finite(x$effect_sd) &
      is.finite(x$mean_average_power) &
      x$effect_sd > 0,
  ]
  x <- x[order(x$effect_sd), ]

  below <- which(x$mean_average_power < target_power)
  above <- which(x$mean_average_power >= target_power)

  if (length(below) == 0L || length(above) == 0L) return(NA_real_)

  i_hi <- min(above)
  i_lo_candidates <- below[below < i_hi]
  if (length(i_lo_candidates) == 0L) return(NA_real_)
  i_lo <- max(i_lo_candidates)

  x1 <- x$effect_sd[i_lo]
  x2 <- x$effect_sd[i_hi]
  y1 <- x$mean_average_power[i_lo]
  y2 <- x$mean_average_power[i_hi]

  if (!is.finite(y2 - y1) || y2 == y1) return(NA_real_)
  x1 + (target_power - y1) * (x2 - x1) / (y2 - y1)
}

interpolated_mde_list <- lapply(
  split(
    nonnull_summary,
    interaction(
      nonnull_summary$resampling_method,
      nonnull_summary$sample_size,
      nonnull_summary$pi1,
      drop = TRUE,
      lex.order = TRUE
    )
  ),
  function(x) {
    data.frame(
      resampling_method = x$resampling_method[[1]],
      sample_size = x$sample_size[[1]],
      pi1 = x$pi1[[1]],
      interpolated_effect_sd_at_80pct_power = interpolate_power_crossing(
        x,
        target_power = TARGET_AVERAGE_POWER
      ),
      note = paste(
        "Exploratory linear interpolation between adjacent simulated effects;",
        "do not report as the primary MDE without direct simulation."
      ),
      stringsAsFactors = FALSE
    )
  }
)
interpolated_mde_table <- do.call(rbind, interpolated_mde_list)

# Global-null audit. Under the complete null, every discovery is false, so the
# empirical FDR equals the probability of at least one discovery.
null_iterations <- iteration_results[
  iteration_results$effect_sd == 0 &
    iteration_results$pi1 == 0,
]

null_split_key <- interaction(
  null_iterations$resampling_method,
  null_iterations$sample_size,
  drop = TRUE,
  lex.order = TRUE
)

null_diagnostics_list <- lapply(
  split(null_iterations, null_split_key),
  function(x) {
    data.frame(
      resampling_method = x$resampling_method[[1]],
      sample_size = x$sample_size[[1]],
      n_sim = nrow(x),
      empirical_fdr = mean(x$false_discovery_proportion),
      probability_any_discovery = mean(x$any_discovery),
      mean_discoveries = mean(x$total_discoveries),
      median_discoveries = stats::median(x$total_discoveries),
      maximum_discoveries = max(x$total_discoveries),
      iterations_with_any_discovery = sum(x$any_discovery),
      stringsAsFactors = FALSE
    )
  }
)
null_diagnostics <- do.call(rbind, null_diagnostics_list)
null_diagnostics <- null_diagnostics[
  order(
    null_diagnostics$resampling_method,
    null_diagnostics$sample_size
  ),
]

null_events <- null_iterations[
  null_iterations$total_discoveries > 0,
  c(
    "resampling_method", "sample_size", "scenario_id", "iteration",
    "false_positives", "total_discoveries",
    "false_discovery_proportion"
  ),
  drop = FALSE
]

# Direct comparison at effects evaluated by both methods.
permutation_summary <- nonnull_summary[
  nonnull_summary$resampling_method == "permutation",
]
bootstrap_summary <- nonnull_summary[
  nonnull_summary$resampling_method == "bootstrap",
]

resampling_comparison <- merge(
  permutation_summary,
  bootstrap_summary,
  by = c("sample_size", "effect_sd", "pi1"),
  suffixes = c("_permutation", "_bootstrap")
)

if (nrow(resampling_comparison) > 0L) {
  resampling_comparison$power_difference_bootstrap_minus_permutation <-
    resampling_comparison$mean_average_power_bootstrap -
    resampling_comparison$mean_average_power_permutation
  resampling_comparison$fdr_difference_bootstrap_minus_permutation <-
    resampling_comparison$mean_fdr_bootstrap -
    resampling_comparison$mean_fdr_permutation
}

###############################################################################
# 11. OUTPUTS
###############################################################################

OUTDIR <- file.path(analysis_root, "08_sample_size_power_simulation")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

write.csv(
  scenario_grid,
  file.path(OUTDIR, "simulation_scenario_grid.csv"),
  row.names = FALSE
)
write.csv(
  simulation_summary,
  file.path(OUTDIR, "simulation_operating_characteristics_summary.csv"),
  row.names = FALSE
)
write.csv(
  mde_table,
  file.path(OUTDIR, "minimum_detectable_effect_grid_80pct_power.csv"),
  row.names = FALSE
)
write.csv(
  analytical_crosscheck,
  file.path(OUTDIR, "analytical_design_based_crosscheck.csv"),
  row.names = FALSE
)

write.csv(
  interpolated_mde_table,
  file.path(OUTDIR, "interpolated_effect_at_80pct_power_EXPLORATORY.csv"),
  row.names = FALSE
)
write.csv(
  null_diagnostics,
  file.path(OUTDIR, "global_null_calibration_summary.csv"),
  row.names = FALSE
)
write.csv(
  null_events,
  file.path(OUTDIR, "global_null_nonzero_discovery_events.csv"),
  row.names = FALSE
)
write.csv(
  resampling_comparison,
  file.path(OUTDIR, "permutation_vs_bootstrap_comparison.csv"),
  row.names = FALSE
)

if (SAVE_ITERATION_LEVEL_RESULTS) {
  write.csv(
    iteration_results,
    file.path(OUTDIR, "simulation_iteration_level_results.csv"),
    row.names = FALSE
  )
}

sample_audit <- as.data.frame(table(meta$Country, meta$SampleGroup))
names(sample_audit) <- c("Country", "SampleGroup", "n")
write.csv(
  sample_audit,
  file.path(OUTDIR, "complete_case_country_diagnosis_counts.csv"),
  row.names = FALSE
)

config <- data.frame(
  parameter = c(
    "workspace_file", "complete_case_n", "tested_analytes",
    "expression_is_log2", "fdr_target", "target_average_power",
    "ebayes_trend", "ebayes_robust", "run_profile",
    "n_sim_permutation", "n_sim_bootstrap",
    "run_bootstrap_sensitivity", "seed",
    "effect_sizes_sd_primary", "effect_sizes_sd_bootstrap",
    "non_null_proportions", "evaluate_sample_size_curve",
    "candidate_sample_sizes"
  ),
  value = c(
    WORKSPACE_FILE, n_observed, m_tests,
    EXPRESSION_IS_LOG2, FDR_TARGET, TARGET_AVERAGE_POWER,
    EBAYES_TREND, EBAYES_ROBUST, RUN_PROFILE,
    N_SIM_PERMUTATION, N_SIM_BOOTSTRAP,
    RUN_BOOTSTRAP_SENSITIVITY, SEED_GLOBAL,
    paste(EFFECT_SIZES_SD_PRIMARY, collapse = ";"),
    paste(EFFECT_SIZES_SD_BOOTSTRAP, collapse = ";"),
    paste(NON_NULL_PROPORTIONS, collapse = ";"),
    EVALUATE_SAMPLE_SIZE_CURVE,
    paste(SAMPLE_SIZES, collapse = ";")
  ),
  stringsAsFactors = FALSE
)
write.csv(config, file.path(OUTDIR, "analysis_configuration.csv"), row.names = FALSE)

# Plot power and empirical FDR separately for each residual-resampling method.
for (method_i in unique(nonnull_summary$resampling_method)) {
  current_n_tbl <- nonnull_summary[
    nonnull_summary$sample_size == n_observed &
      nonnull_summary$resampling_method == method_i,
  ]

  if (nrow(current_n_tbl) == 0L) next

  pi_values <- sort(unique(current_n_tbl$pi1))
  effect_values <- sort(unique(current_n_tbl$effect_sd))

  ymat <- sapply(pi_values, function(p_i) {
    z <- current_n_tbl[current_n_tbl$pi1 == p_i, ]
    z$mean_average_power[match(effect_values, z$effect_sd)]
  })

  png(
    file.path(
      OUTDIR,
      paste0("power_curves_observed_sample_", method_i, ".png")
    ),
    width = 1800, height = 1200, res = 180
  )
  matplot(
    effect_values, ymat,
    type = "b", pch = seq_along(pi_values),
    lty = seq_along(pi_values),
    xlab = "Prespecified effect (residual SD units)",
    ylab = "Mean average power",
    ylim = c(0, 1),
    main = paste0(
      "ReDLat limma simulation: ", method_i,
      ", complete-case N = ", n_observed
    )
  )
  abline(h = TARGET_AVERAGE_POWER, lty = 2)
  legend(
    "bottomright",
    legend = paste0("Non-null proportion = ", pi_values),
    pch = seq_along(pi_values),
    lty = seq_along(pi_values),
    bty = "n"
  )
  dev.off()

  fdr_mat <- sapply(pi_values, function(p_i) {
    z <- current_n_tbl[current_n_tbl$pi1 == p_i, ]
    z$mean_fdr[match(effect_values, z$effect_sd)]
  })

  finite_fdr <- fdr_mat[is.finite(fdr_mat)]
  fdr_ymax <- if (length(finite_fdr) > 0L) {
    min(1, max(FDR_TARGET * 1.5, finite_fdr))
  } else {
    FDR_TARGET * 1.5
  }

  png(
    file.path(
      OUTDIR,
      paste0("fdr_curves_observed_sample_", method_i, ".png")
    ),
    width = 1800, height = 1200, res = 180
  )
  matplot(
    effect_values, fdr_mat,
    type = "b", pch = seq_along(pi_values),
    lty = seq_along(pi_values),
    xlab = "Prespecified effect (residual SD units)",
    ylab = "Mean empirical FDR",
    ylim = c(0, fdr_ymax),
    main = paste0(
      "ReDLat limma simulation: ", method_i,
      ", empirical FDR at N = ", n_observed
    )
  )
  abline(h = FDR_TARGET, lty = 2)
  legend(
    "topright",
    legend = paste0("Non-null proportion = ", pi_values),
    pch = seq_along(pi_values),
    lty = seq_along(pi_values),
    bty = "n"
  )
  dev.off()
}

guardrails <- c(
  "INTERPRETATION GUARDRAILS",
  "",
  "1. Describe this as a retrospective simulation-based sensitivity analysis,",
  "   not as an a priori sample-size determination.",
  "2. Do not report observed post hoc power calculated from significant effects.",
  "3. The primary result is the residual-permutation analysis. The residual",
  "   bootstrap is a sensitivity analysis and must be reported separately.",
  "4. Report the smallest DIRECTLY EVALUATED effect meeting average power >=0.80",
  "   and empirical mean FDR <=0.05. Treat interpolation as exploratory.",
  "5. Report all prespecified non-null-proportion scenarios; do not select only",
  "   the most favorable scenario.",
  "6. Confirm global-null calibration for both resampling methods.",
  "7. The multiplicity universe is the aptamer-level feature set submitted to",
  "   limma and BH, not the later gene-collapsed reporting universe.",
  "8. Effects are in empirical residual-SD units. The median log2-RFU, fold-change",
  "   and percentage translations summarize the imposed signals across analytes;",
  "   they are not universal per-protein thresholds.",
  "9. Do not claim that this analysis establishes power for WGCNA module recovery",
  "   or machine-learning AUC precision; those require separate simulations.",
  "10. PROPER and MSstatsSampleSize informed the design, but were not executed.",
  "    Describe this as custom ReDLat code."
)
writeLines(
  guardrails,
  con = file.path(OUTDIR, "interpretation_guardrails.txt")
)

methods_template <- c(
  "SUGGESTED SUPPLEMENTARY METHODS",
  "",
  paste0(
    "No statistical methods were used to predetermine the sample size. All ",
    "eligible participants with available plasma proteomic measurements passing ",
    "predefined quality-control criteria were included. To contextualize the ",
    "operating characteristics supported by the available cohort, we performed ",
    "a retrospective simulation-based sensitivity analysis reproducing the ",
    "primary aptamer-level differential-abundance pipeline. The analysis used ",
    n_observed, " complete-case participants and ", m_tests,
    " aptamer-level features, and refitted the same covariate-adjusted limma ",
    "model with empirical-Bayes moderation and Benjamini-Hochberg false-",
    "discovery-rate control. We retained the observed diagnosis-by-country ",
    "composition, fitted age, sex, education and country effects, feature-",
    "specific residual variation and cross-feature dependence by resampling ",
    "complete participant-level residual vectors within country. Residual ",
    "permutation was prespecified as the primary analysis, and residual ",
    "bootstrap resampling was used as a sensitivity analysis. Diagnosis effects ",
    "were introduced under prespecified standardized effect sizes and assumed ",
    "non-null proportions, without using the effects of the observed significant ",
    "features as ground truth. Average power and empirical false discovery rate ",
    "were estimated across ", N_SIM_PERMUTATION,
    " primary simulations per scenario; the bootstrap sensitivity used ",
    N_SIM_BOOTSTRAP, " simulations per scenario."
  ),
  "",
  "Insert numerical results only after auditing the final output tables.",
  "",
  "Do not state that the sample size was prospectively powered."
)
writeLines(
  methods_template,
  con = file.path(OUTDIR, "supplementary_methods_template.txt")
)

references_and_provenance <- c(
  "METHODS, REFERENCES AND PROVENANCE",
  "",
  "STATUS",
  "This is custom ReDLat simulation code. It directly executes limma and the",
  "Benjamini-Hochberg adjustment, but it does not execute PROPER,",
  "FDRsamplesize2 or MSstatsSampleSize.",
  "",
  "RATIONALE",
  "High-dimensional power depends on the fraction and distribution of true",
  "signals, empirical variation, feature dependence, the statistical model and",
  "multiplicity correction. The custom simulation therefore reruns the exact",
  "primary ReDLat model under prespecified effects and non-null proportions.",
  "",
  "REFERENCES",
  "",
  "1. Ritchie ME, Phipson B, Wu D, et al. limma powers differential expression",
  "   analyses for RNA-sequencing and microarray studies. Nucleic Acids Research.",
  "   2015;43(7):e47. doi:10.1093/nar/gkv007.",
  "   Supports the empirical-Bayes moderated linear-model framework.",
  "",
  "2. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical",
  "   and powerful approach to multiple testing. Journal of the Royal Statistical",
  "   Society Series B. 1995;57(1):289-300.",
  "   doi:10.1111/j.2517-6161.1995.tb02031.x.",
  "   Supports the BH false-discovery-rate procedure reproduced in every run.",
  "",
  "3. Wu H, Wang C, Wu Z. PROPER: comprehensive power evaluation for differential",
  "   expression using RNA-seq. Bioinformatics. 2015;31(2):233-241.",
  "   doi:10.1093/bioinformatics/btu640.",
  "   Supports semi-parametric simulation based on real high-dimensional data,",
  "   explicit non-null proportions and direct evaluation of power and errors.",
  "",
  "4. Ni Y, Eames Seffernick A, Onar-Thomas A. Computing power and sample size",
  "   for the false discovery rate in multiple applications. Genes.",
  "   2024;15(3):344. doi:10.3390/genes15030344.",
  "   Supports average-power calculations under FDR-based multiple testing and",
  "   motivates the analytical cross-check.",
  "",
  "5. Huang T, Choi M, Vitek O. MSstatsSampleSize: simulation tool for optimal",
  "   design of high-dimensional MS-based proteomics experiments.",
  "   Bioconductor DOI:10.18129/B9.bioc.MSstatsSampleSize.",
  "   GitHub: https://github.com/Vitek-Lab/MSstatsSampleSize",
  "   Provides a proteomics-specific precedent for using empirical abundance",
  "   variation and repeated simulation. It was not directly used because its",
  "   MS-specific model does not reproduce the ReDLat SOMAscan covariate-adjusted",
  "   limma analysis.",
  "",
  "REPOSITORIES",
  "limma: https://bioconductor.org/packages/limma",
  "PROPER: https://bioconductor.org/packages/PROPER",
  "MSstatsSampleSize: https://github.com/Vitek-Lab/MSstatsSampleSize"
)
writeLines(
  references_and_provenance,
  con = file.path(OUTDIR, "METHODS_AND_REFERENCES.txt")
)

# Compact automatic audit report.
audit_lines <- c(
  "AUTOMATIC AUDIT SUMMARY",
  "",
  paste0("Run profile: ", RUN_PROFILE),
  paste0("Complete-case N: ", n_observed),
  paste0("Tested aptamer-level features: ", m_tests),
  paste0("Primary simulations per scenario: ", N_SIM_PERMUTATION),
  paste0("Bootstrap simulations per scenario: ", N_SIM_BOOTSTRAP),
  "",
  "GLOBAL-NULL CALIBRATION",
  capture.output(print(null_diagnostics, row.names = FALSE)),
  "",
  "GRID-BASED MINIMUM EVALUATED EFFECTS",
  capture.output(print(mde_table, row.names = FALSE)),
  "",
  "Interpret numerical findings only after manual review of all output tables."
)
writeLines(
  audit_lines,
  con = file.path(OUTDIR, "automatic_audit_summary.txt")
)

write_session_info(file.path(OUTDIR, "sessionInfo.txt"))

###############################################################################
# 12. COMPACT CONSOLE SUMMARY
###############################################################################

cat("\n============================================================\n")
cat("SIMULATION-BASED SENSITIVITY ANALYSIS\n")
cat("============================================================\n")
cat("Complete-case N: ", n_observed, "\n", sep = "")
cat("Tested aptamer-level features: ", m_tests, "\n", sep = "")
cat("Run profile: ", RUN_PROFILE, "\n", sep = "")

cat("\nGLOBAL-NULL CALIBRATION\n")
print(null_diagnostics, row.names = FALSE)

cat("\nMINIMUM DIRECTLY EVALUATED EFFECTS MEETING BOTH TARGETS\n")
print(mde_table, row.names = FALSE)

cat("\nANALYTICAL DESIGN-BASED CROSS-CHECK\n")
print(analytical_crosscheck, row.names = FALSE)

message("Completed.")
message("Outputs written to: ", OUTDIR)

###############################################################################
# END
###############################################################################

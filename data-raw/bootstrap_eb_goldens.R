# Bootstrap script for Sub-phase 8.6.B EB goldens (2026-05-25 closure
# session). Parallel to data-raw/bootstrap_K2_goldens.R and similar
# 8.3.X bootstrap scripts; mirrors the env-var gating pattern of
# tests/testthat/test-golden_regression_8_3_9.R.
#
# Usage:
#
#   GDPAR_BOOTSTRAP_EB_GOLDENS=1 Rscript data-raw/bootstrap_eb_goldens.R
#
# When the env var is set, the script:
#   1. Bootstraps the 5 fitable goldens of .gdpar_golden_roster_8_6_B()
#      with seed = 20260525L and saves them under
#      tests/testthat/data/golden_eb_<family>_<regime>.rds.
#   2. Bootstraps the 12 guard-rejected goldens as metadata-only RDS
#      files (no Stan; capture only the expected condition class so the
#      structural tests can verify the guard fires as documented).
#   3. Appends one row per golden to tests/testthat/data/golden_manifest.csv
#      with sub_phase = "8.6.B" (16 canonical columns of
#      .gdpar_golden_manifest_columns()).
#
# Without the env var, the script is a no-op informational message.

if (!nzchar(Sys.getenv("GDPAR_BOOTSTRAP_EB_GOLDENS"))) {
  message("GDPAR_BOOTSTRAP_EB_GOLDENS not set; bootstrap is a no-op.")
  quit(save = "no", status = 0L)
}

suppressMessages({
  devtools::load_all(".", quiet = TRUE)
})

stopifnot(requireNamespace("cmdstanr", quietly = TRUE))

set.seed(20260525L)
data_dir <- file.path("tests", "testthat", "data")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
manifest_path <- file.path(data_dir, "golden_manifest.csv")
gdpar:::.gdpar_golden_manifest_init(manifest_path)

R_version_str <- paste0(R.version$major, ".", R.version$minor)
cmdstan_version_str <- gdpar:::.gdpar_cmdstan_version()
bootstrap_date <- format(Sys.Date(), "%Y-%m-%d")

.simulate_data <- function(family, n = 80L, seed = 20260525L) {
  set.seed(seed)
  df <- data.frame(x1 = rnorm(n))
  if (family == "gaussian") {
    df$y <- 0.5 + 0.3 * df$x1 + rnorm(n, sd = 0.4)
  } else if (family == "poisson") {
    df$y <- rpois(n, lambda = exp(0.6 + 0.2 * df$x1))
  } else if (family == "neg_binomial_2") {
    df$y <- rnbinom(n, size = 5, mu = exp(0.5 + 0.2 * df$x1))
  } else if (family == "bernoulli") {
    df$y <- rbinom(n, size = 1L,
                   prob = plogis(0.3 + 0.4 * df$x1))
  } else {
    stop(sprintf("No data generator for family '%s'.", family))
  }
  df
}

.bootstrap_one_fitable <- function(cfg, seed = 20260525L) {
  family <- cfg$family
  regime <- cfg$regime
  fname <- gdpar:::.gdpar_golden_eb_rds_name(family, regime)
  rds_path_existing <- file.path(data_dir, fname)
  if (file.exists(rds_path_existing) &&
      !nzchar(Sys.getenv("GDPAR_BOOTSTRAP_EB_FORCE"))) {
    message(sprintf("[EB-bootstrap] %s/%s already present; skipping (set GDPAR_BOOTSTRAP_EB_FORCE=1 to refit).",
                    family, regime))
    return(invisible(rds_path_existing))
  }
  message(sprintf("[EB-bootstrap] fitting %s / %s ...", family, regime))
  df <- .simulate_data(family, n = 80L, seed = seed)
  amm_obj <- if (regime == "bspline") {
    amm_spec(a = ~ x1,
             W = W_basis(type = "bspline", degree = 3L,
                         boundary_knots = c(-3, 3), df = 4L))
  } else {
    amm_spec(a = ~ x1)
  }
  # B-spline diagnostic identifiability test point sits at the anchor;
  # for the bootstrap we trust the structural template and skip the
  # data-driven check (same convention used in the 8.3.8 fit-smoke
  # bootstraps).
  skip_id <- (regime == "bspline")
  fit_eb <- gdpar_eb(
    formula       = y ~ x1,
    family        = gdpar_family(family),
    amm           = amm_obj,
    data          = df,
    iter_warmup   = 200L, iter_sampling = 200L, chains = 2L,
    refresh       = 0L, seed = seed,
    skip_id_check = skip_id,
    laplace_control = list(multi_start_M = 2L, laplace_draws = 200L),
    verbose       = FALSE
  )
  draws <- fit_eb$conditional_fit$draws()
  draws_matrix <- posterior::as_draws_matrix(draws)
  snapshot <- list(
    family             = family,
    regime             = regime,
    status             = "fitable",
    sub_phase          = "8.6.B",
    bootstrap_date     = bootstrap_date,
    cmdstan_version    = cmdstan_version_str,
    R_version          = R_version_str,
    seed               = seed,
    n_chains           = 2L,
    n_iter_warmup      = 200L,
    n_iter_sampling    = 200L,
    n_obs              = nrow(df),
    theta_ref_hat      = fit_eb$theta_ref_hat,
    theta_ref_se       = fit_eb$theta_ref_se,
    correction_applied = fit_eb$correction_applied,
    eb_correction_constant = fit_eb$eb_correction_constant,
    diagnostics_numerical  = fit_eb$diagnostics_numerical,
    draws_matrix       = draws_matrix
  )
  rds_path <- file.path(data_dir, fname)
  saveRDS(snapshot, rds_path)
  gdpar:::.gdpar_golden_manifest_add(
    list(
      family = family, K = 1L, basis_type = regime, p = 1L,
      sub_phase = "8.6.B", bootstrap_date = bootstrap_date,
      cmdstan_version = cmdstan_version_str, R_version = R_version_str,
      DHARMa_version = gdpar:::.gdpar_dharma_version(),
      fit_code_hash = gdpar:::.gdpar_golden_fit_code_hash(.bootstrap_one_fitable),
      seed = seed, n_chains = 2L, n_iter_warmup = 200L,
      n_iter_sampling = 200L,
      rds_path = file.path("data", fname), n_obs = nrow(df)
    ),
    path = manifest_path
  )
  invisible(rds_path)
}

.bootstrap_one_guard <- function(cfg) {
  family <- cfg$family
  regime <- cfg$regime
  fname <- gdpar:::.gdpar_golden_eb_rds_name(family, regime)
  message(sprintf("[EB-bootstrap] metadata-only guard %s / %s ...",
                  family, regime))
  snapshot <- list(
    family          = family,
    regime          = regime,
    status          = "guard",
    sub_phase       = "8.6.B",
    bootstrap_date  = bootstrap_date,
    cmdstan_version = cmdstan_version_str,
    R_version       = R_version_str,
    guard_class     = cfg$guard_class,
    rationale       = cfg$rationale
  )
  rds_path <- file.path(data_dir, fname)
  saveRDS(snapshot, rds_path)
  gdpar:::.gdpar_golden_manifest_add(
    list(
      family = family, K = 1L, basis_type = regime, p = 1L,
      sub_phase = "8.6.B", bootstrap_date = bootstrap_date,
      cmdstan_version = cmdstan_version_str, R_version = R_version_str,
      DHARMa_version = gdpar:::.gdpar_dharma_version(),
      fit_code_hash = gdpar:::.gdpar_golden_fit_code_hash(.bootstrap_one_guard),
      seed = NA_integer_, n_chains = NA_integer_,
      n_iter_warmup = NA_integer_, n_iter_sampling = NA_integer_,
      rds_path = file.path("data", fname), n_obs = NA_integer_
    ),
    path = manifest_path
  )
  invisible(rds_path)
}

roster <- gdpar:::.gdpar_golden_roster_8_6_B()
for (cfg in roster) {
  if (cfg$status == "fitable") {
    .bootstrap_one_fitable(cfg)
  } else {
    .bootstrap_one_guard(cfg)
  }
}

message(sprintf("[EB-bootstrap] DONE. Wrote %d EB goldens + manifest rows.",
                length(roster)))

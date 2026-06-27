# Bootstrap script for Sub-phase 8.6.D EB goldens (Session 13c closure,
# 2026-05-26). Parallel to data-raw/bootstrap_eb_goldens_C.R (8.6.C),
# mirrors the env-var gating pattern but writes Path C golden
# snapshots under the .gdpar_golden_eb_rds_name_8_6_D() convention so
# the on-disk namespaces of 8.6.B / 8.6.C / 8.6.D never collide.
#
# Usage:
#
#   GDPAR_BOOTSTRAP_EB_GOLDENS=1 \
#   GDPAR_BOOTSTRAP_EB_GOLDENS_D=1 \
#     Rscript data-raw/bootstrap_eb_goldens_D.R
#
# Both env vars must be set: the outer flag GDPAR_BOOTSTRAP_EB_GOLDENS
# matches the 8.6.B / 8.6.C convention (a single switch for any EB
# bootstrap), and the inner flag GDPAR_BOOTSTRAP_EB_GOLDENS_D scopes
# the run to the 8.6.D roster so a user can rebuild only one sub-phase
# selectively. Without both flags the script is a no-op informational
# message.
#
# When both env vars are set, the script:
#   1. Bootstraps the 4 fitable Path C goldens of
#      .gdpar_golden_roster_8_6_D() with seed = 20260526L and saves
#      them under
#      tests/testthat/data/golden_eb_<family>_pathC_<regime>.rds.
#   2. Bootstraps the 4 guard-rejected Path C goldens as metadata-only
#      RDS files (no Stan; capture only the expected condition class
#      so the structural tests can verify the guard fires as
#      documented per D40').
#   3. Appends one row per golden to
#      tests/testthat/data/golden_manifest.csv with
#      sub_phase = "8.6.D" (16 canonical columns).

if (!nzchar(Sys.getenv("GDPAR_BOOTSTRAP_EB_GOLDENS")) ||
    !nzchar(Sys.getenv("GDPAR_BOOTSTRAP_EB_GOLDENS_D"))) {
  message("GDPAR_BOOTSTRAP_EB_GOLDENS / GDPAR_BOOTSTRAP_EB_GOLDENS_D not set; bootstrap_eb_goldens_D is a no-op.")
  quit(save = "no", status = 0L)
}

suppressMessages({
  devtools::load_all(".", quiet = TRUE)
})

stopifnot(requireNamespace("cmdstanr", quietly = TRUE))
stopifnot(requireNamespace("posterior", quietly = TRUE))

set.seed(20260526L)
data_dir <- file.path("tests", "testthat", "data")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
manifest_path <- file.path(data_dir, "golden_manifest.csv")
gdpar:::.gdpar_golden_manifest_init(manifest_path)

R_version_str <- paste0(R.version$major, ".", R.version$minor)
cmdstan_version_str <- gdpar:::.gdpar_cmdstan_version()
bootstrap_date <- format(Sys.Date(), "%Y-%m-%d")

# ---- Path C K x p fit-code closures ------------------------------------

.fit_path_C_gaussian_polyP2 <- function(seed) {
  set.seed(seed)
  n <- 80L
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  y_mat <- matrix(NA_real_, n, 2L)
  for (j in 1:2) {
    mu_j <- c(1.0, -0.5)[j] + c(0.5, -0.3)[j] * d$x1
    ls_j <- c(-0.2, 0.0)[j] + c(0.2, 0.1)[j] * d$x2
    y_mat[, j] <- rnorm(n, mean = mu_j, sd = exp(ls_j))
  }
  d$y <- y_mat
  gdpar_eb(
    formula = y ~ x1 + x2,
    family  = gdpar_family("gaussian"),
    amm     = list(
      mu    = amm_spec(p = 2L, dims = dimwise(a = ~ x1)),
      sigma = amm_spec(p = 2L, dims = dimwise(a = ~ x2))
    ),
    data    = d,
    iter_warmup = 200L, iter_sampling = 200L, chains = 2L,
    refresh = 0L, seed = seed, skip_id_check = TRUE,
    laplace_control = list(multi_start_M = 2L, laplace_draws = 200L),
    verbose = FALSE
  )
}

.fit_path_C_gaussian_polyP3 <- function(seed) {
  set.seed(seed)
  n <- 100L
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  y_mat <- matrix(NA_real_, n, 3L)
  for (j in 1:3) {
    mu_j <- c(1.0, -0.5, 0.2)[j] + c(0.5, -0.3, 0.1)[j] * d$x1
    ls_j <- c(-0.2, 0.0, -0.1)[j] + c(0.2, 0.1, 0.0)[j] * d$x2
    y_mat[, j] <- rnorm(n, mean = mu_j, sd = exp(ls_j))
  }
  d$y <- y_mat
  gdpar_eb(
    formula = y ~ x1 + x2,
    family  = gdpar_family("gaussian"),
    amm     = list(
      mu    = amm_spec(p = 3L, dims = dimwise(a = ~ x1)),
      sigma = amm_spec(p = 3L, dims = dimwise(a = ~ x2))
    ),
    data    = d,
    iter_warmup = 200L, iter_sampling = 200L, chains = 2L,
    refresh = 0L, seed = seed, skip_id_check = TRUE,
    laplace_control = list(multi_start_M = 2L, laplace_draws = 200L),
    verbose = FALSE
  )
}

.fit_path_C_neg_binomial_2_polyP2 <- function(seed) {
  set.seed(seed)
  n <- 100L
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  y_mat <- matrix(NA_integer_, n, 2L)
  for (j in 1:2) {
    eta_mu <- 1.2 + c(0.2, -0.1)[j] * d$x1
    eta_phi <- 0.5 + c(0.1, 0.0)[j] * d$x2
    mu_ij <- exp(eta_mu)
    phi_ij <- exp(eta_phi)
    y_mat[, j] <- rnbinom(n, mu = mu_ij, size = phi_ij)
  }
  d$y <- y_mat
  gdpar_eb(
    formula = y ~ x1 + x2,
    family  = gdpar_family("neg_binomial_2"),
    amm     = list(
      mu  = amm_spec(p = 2L, dims = dimwise(a = ~ x1)),
      phi = amm_spec(p = 2L, dims = dimwise(a = ~ x2))
    ),
    data    = d,
    iter_warmup = 200L, iter_sampling = 200L, chains = 2L,
    refresh = 0L, seed = seed, skip_id_check = TRUE,
    laplace_control = list(multi_start_M = 2L, laplace_draws = 200L),
    verbose = FALSE
  )
}

.fit_path_C_neg_binomial_2_polyP3 <- function(seed) {
  set.seed(seed)
  n <- 120L
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  y_mat <- matrix(NA_integer_, n, 3L)
  for (j in 1:3) {
    eta_mu <- 1.0 + c(0.2, -0.1, 0.05)[j] * d$x1
    eta_phi <- 0.5 + c(0.1, 0.0, -0.05)[j] * d$x2
    mu_ij <- exp(eta_mu)
    phi_ij <- exp(eta_phi)
    y_mat[, j] <- rnbinom(n, mu = mu_ij, size = phi_ij)
  }
  d$y <- y_mat
  gdpar_eb(
    formula = y ~ x1 + x2,
    family  = gdpar_family("neg_binomial_2"),
    amm     = list(
      mu  = amm_spec(p = 3L, dims = dimwise(a = ~ x1)),
      phi = amm_spec(p = 3L, dims = dimwise(a = ~ x2))
    ),
    data    = d,
    iter_warmup = 200L, iter_sampling = 200L, chains = 2L,
    refresh = 0L, seed = seed, skip_id_check = TRUE,
    laplace_control = list(multi_start_M = 2L, laplace_draws = 200L),
    verbose = FALSE
  )
}

# ---- Dispatcher: roster entry -> fit-code closure -----------------------

.dispatch_fit_code_D <- function(cfg) {
  key <- paste(cfg$family, cfg$regime, sep = "/")
  switch(
    key,
    "gaussian/polyP2"         = .fit_path_C_gaussian_polyP2,
    "gaussian/polyP3"         = .fit_path_C_gaussian_polyP3,
    "neg_binomial_2/polyP2"   = .fit_path_C_neg_binomial_2_polyP2,
    "neg_binomial_2/polyP3"   = .fit_path_C_neg_binomial_2_polyP3,
    stop(sprintf("No fit-code dispatch for C/%s", key))
  )
}

.bootstrap_one_fitable_D <- function(cfg, seed = 20260526L) {
  family <- cfg$family
  fname <- gdpar:::.gdpar_golden_eb_rds_name_8_6_D(family, cfg$regime)
  rds_path_existing <- file.path(data_dir, fname)
  if (file.exists(rds_path_existing) &&
      !nzchar(Sys.getenv("GDPAR_BOOTSTRAP_EB_FORCE"))) {
    message(sprintf("[EB-8.6.D-bootstrap] C/%s/%s already present; skipping.",
                    family, cfg$regime))
    return(invisible(rds_path_existing))
  }
  message(sprintf("[EB-8.6.D-bootstrap] fitting C / %s / %s ...",
                  family, cfg$regime))
  fit_code_fn <- .dispatch_fit_code_D(cfg)
  fit_eb <- fit_code_fn(seed)
  draws <- fit_eb$conditional_fit$draws()
  draws_matrix <- posterior::as_draws_matrix(draws)
  snapshot <- list(
    family                       = family,
    regime                       = cfg$regime,
    path                         = "C",
    K                            = cfg$K,
    p                            = cfg$p,
    status                       = "fitable",
    sub_phase                    = "8.6.D",
    bootstrap_date               = bootstrap_date,
    cmdstan_version              = cmdstan_version_str,
    R_version                    = R_version_str,
    seed                         = seed,
    n_chains                     = 2L,
    n_iter_warmup                = 200L,
    n_iter_sampling              = 200L,
    theta_ref_kp_hat             = fit_eb$theta_ref_kp_hat,
    theta_ref_kp_se              = fit_eb$theta_ref_kp_se,
    theta_ref_kp_cov_per_slot    = fit_eb$theta_ref_kp_cov_per_slot,
    correction_applied           = fit_eb$correction_applied,
    correction_tensor_constant   = fit_eb$correction_tensor_constant,
    correction_tensor_dispositions = fit_eb$correction_tensor_dispositions,
    diagnostics_numerical        = fit_eb$diagnostics_numerical,
    draws_matrix                 = draws_matrix
  )
  rds_path <- file.path(data_dir, fname)
  saveRDS(snapshot, rds_path)
  gdpar:::.gdpar_golden_manifest_add(
    list(
      family = family, K = cfg$K, basis_type = cfg$regime, p = cfg$p,
      sub_phase = "8.6.D", bootstrap_date = bootstrap_date,
      cmdstan_version = cmdstan_version_str, R_version = R_version_str,
      DHARMa_version = gdpar:::.gdpar_dharma_version(),
      fit_code_hash = gdpar:::.gdpar_golden_fit_code_hash(fit_code_fn),
      seed = seed, n_chains = 2L, n_iter_warmup = 200L,
      n_iter_sampling = 200L,
      rds_path = file.path("data", fname),
      n_obs = NA_integer_
    ),
    path = manifest_path
  )
  invisible(rds_path)
}

.bootstrap_one_guard_D <- function(cfg) {
  family <- cfg$family
  fname <- gdpar:::.gdpar_golden_eb_rds_name_8_6_D(family, cfg$regime)
  message(sprintf("[EB-8.6.D-bootstrap] metadata-only guard C / %s / %s ...",
                  family, cfg$regime))
  snapshot <- list(
    family          = family,
    regime          = cfg$regime,
    path            = "C",
    K               = cfg$K,
    p               = cfg$p,
    status          = "guard",
    sub_phase       = "8.6.D",
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
      family = family, K = cfg$K, basis_type = cfg$regime, p = cfg$p,
      sub_phase = "8.6.D", bootstrap_date = bootstrap_date,
      cmdstan_version = cmdstan_version_str, R_version = R_version_str,
      DHARMa_version = gdpar:::.gdpar_dharma_version(),
      fit_code_hash = gdpar:::.gdpar_golden_fit_code_hash(.bootstrap_one_guard_D),
      seed = NA_integer_, n_chains = NA_integer_,
      n_iter_warmup = NA_integer_, n_iter_sampling = NA_integer_,
      rds_path = file.path("data", fname),
      n_obs = NA_integer_
    ),
    path = manifest_path
  )
  invisible(rds_path)
}

roster <- gdpar:::.gdpar_golden_roster_8_6_D()
for (cfg in roster) {
  if (cfg$status == "fitable") {
    .bootstrap_one_fitable_D(cfg)
  } else {
    .bootstrap_one_guard_D(cfg)
  }
}

message(sprintf("[EB-8.6.D-bootstrap] DONE. Wrote %d EB goldens + manifest rows.",
                length(roster)))

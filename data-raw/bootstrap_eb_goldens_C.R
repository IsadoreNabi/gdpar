# Bootstrap script for Sub-phase 8.6.C EB goldens (Session 12 closure
# of 8.6.C, 2026-05-25). Parallel to data-raw/bootstrap_eb_goldens.R
# (8.6.B), mirrors the env-var gating pattern of the 8.6.B script
# but writes Path A / Path B golden snapshots under the
# .gdpar_golden_eb_rds_name_8_6_C() convention so the on-disk
# namespaces of 8.6.B and 8.6.C never collide.
#
# Usage:
#
#   GDPAR_BOOTSTRAP_EB_GOLDENS=1 \
#   GDPAR_BOOTSTRAP_EB_GOLDENS_C=1 \
#     Rscript data-raw/bootstrap_eb_goldens_C.R
#
# Both env vars must be set: the outer flag GDPAR_BOOTSTRAP_EB_GOLDENS
# matches the 8.6.B convention (a single switch for any EB bootstrap),
# and the inner flag GDPAR_BOOTSTRAP_EB_GOLDENS_C scopes the run to
# the 8.6.C roster so a user can rebuild only 8.6.B or only 8.6.C
# selectively. Without both flags the script is a no-op informational
# message.
#
# When both env vars are set, the script:
#   1. Bootstraps the 13 fitable goldens of .gdpar_golden_roster_8_6_C()
#      with seed = 20260525L and saves them under
#      tests/testthat/data/golden_eb_<family>_path<A|B>_<regime>.rds.
#   2. Bootstraps the 4 guard-rejected Path A goldens as metadata-only
#      RDS files (no Stan; capture only the expected condition class
#      so the structural tests can verify the guard fires as
#      documented).
#   3. Appends one row per golden to
#      tests/testthat/data/golden_manifest.csv with
#      sub_phase = "8.6.C" (16 canonical columns).

if (!nzchar(Sys.getenv("GDPAR_BOOTSTRAP_EB_GOLDENS")) ||
    !nzchar(Sys.getenv("GDPAR_BOOTSTRAP_EB_GOLDENS_C"))) {
  message("GDPAR_BOOTSTRAP_EB_GOLDENS / GDPAR_BOOTSTRAP_EB_GOLDENS_C not set; bootstrap_eb_goldens_C is a no-op.")
  quit(save = "no", status = 0L)
}

suppressMessages({
  devtools::load_all(".", quiet = TRUE)
})

stopifnot(requireNamespace("cmdstanr", quietly = TRUE))
stopifnot(requireNamespace("posterior", quietly = TRUE))

set.seed(20260525L)
data_dir <- file.path("tests", "testthat", "data")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
manifest_path <- file.path(data_dir, "golden_manifest.csv")
gdpar:::.gdpar_golden_manifest_init(manifest_path)

R_version_str <- paste0(R.version$major, ".", R.version$minor)
cmdstan_version_str <- gdpar:::.gdpar_cmdstan_version()
bootstrap_date <- format(Sys.Date(), "%Y-%m-%d")

# ---- Path A K=1 p>1 fit-code closures -----------------------------------

.fit_path_A_gaussian_polyP2 <- function(seed) {
  set.seed(seed)
  n <- 80L
  df <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  df$y <- matrix(rnorm(n * 2L, sd = 0.4), n, 2L)
  spec <- amm_spec(p = 2L, dims = dimwise(a = ~ x1 + x2, b = NULL))
  gdpar_eb(
    formula = y ~ x1 + x2,
    family  = gdpar_family_multi("gaussian", p = 2L),
    amm     = spec, data = df,
    iter_warmup = 200L, iter_sampling = 200L, chains = 2L,
    refresh = 0L, seed = seed,
    laplace_control = list(multi_start_M = 2L, laplace_draws = 200L),
    verbose = FALSE
  )
}

.fit_path_A_gaussian_bsplineP2 <- function(seed) {
  set.seed(seed)
  n <- 100L
  df <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  df$y <- matrix(rnorm(n * 2L, sd = 0.4), n, 2L)
  Wspec <- W_basis(type = "bspline", degree = 3L,
                   boundary_knots = c(-3, 3), df = 4L, p = 2L)
  spec <- amm_spec(p = 2L, dims = dimwise(a = ~ x1 + x2, b = NULL),
                   W = Wspec)
  gdpar_eb(
    formula = y ~ x1 + x2,
    family  = gdpar_family_multi("gaussian", p = 2L),
    amm     = spec, data = df,
    iter_warmup = 200L, iter_sampling = 200L, chains = 2L,
    refresh = 0L, seed = seed, skip_id_check = TRUE,
    laplace_control = list(multi_start_M = 2L, laplace_draws = 200L),
    verbose = FALSE
  )
}

.fit_path_A_gaussian_polyP3 <- function(seed) {
  set.seed(seed)
  n <- 80L
  df <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  df$y <- matrix(rnorm(n * 3L, sd = 0.4), n, 3L)
  spec <- amm_spec(p = 3L,
                   dims = dimwise(a = ~ x1 + x2 + x3, b = NULL))
  gdpar_eb(
    formula = y ~ x1 + x2 + x3,
    family  = gdpar_family_multi("gaussian", p = 3L),
    amm     = spec, data = df,
    iter_warmup = 200L, iter_sampling = 200L, chains = 2L,
    refresh = 0L, seed = seed,
    laplace_control = list(multi_start_M = 2L, laplace_draws = 200L),
    verbose = FALSE
  )
}

.fit_path_A_poisson_polyP2 <- function(seed) {
  set.seed(seed)
  n <- 80L
  df <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  df$y <- matrix(rpois(n * 2L, lambda = exp(0.5)), n, 2L)
  spec <- amm_spec(p = 2L, dims = dimwise(a = ~ x1 + x2, b = NULL))
  gdpar_eb(
    formula = y ~ x1 + x2,
    family  = gdpar_family_multi("poisson", p = 2L),
    amm     = spec, data = df,
    iter_warmup = 200L, iter_sampling = 200L, chains = 2L,
    refresh = 0L, seed = seed,
    laplace_control = list(multi_start_M = 2L, laplace_draws = 200L),
    verbose = FALSE
  )
}

.fit_path_A_neg_binomial_2_polyP2 <- function(seed) {
  set.seed(seed)
  n <- 80L
  df <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  df$y <- matrix(rnbinom(n * 2L, size = 5, mu = exp(0.5)), n, 2L)
  spec <- amm_spec(p = 2L, dims = dimwise(a = ~ x1 + x2, b = NULL))
  gdpar_eb(
    formula = y ~ x1 + x2,
    family  = gdpar_family_multi("neg_binomial_2", p = 2L),
    amm     = spec, data = df,
    iter_warmup = 200L, iter_sampling = 200L, chains = 2L,
    refresh = 0L, seed = seed,
    laplace_control = list(multi_start_M = 2L, laplace_draws = 200L),
    verbose = FALSE
  )
}

.fit_path_A_bernoulli_polyP2 <- function(seed) {
  set.seed(seed)
  n <- 120L
  df <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  df$y <- matrix(rbinom(n * 2L, 1L, prob = plogis(0.2)), n, 2L)
  spec <- amm_spec(p = 2L, dims = dimwise(a = ~ x1 + x2, b = NULL))
  gdpar_eb(
    formula = y ~ x1 + x2,
    family  = gdpar_family_multi("bernoulli", p = 2L),
    amm     = spec, data = df,
    iter_warmup = 200L, iter_sampling = 200L, chains = 2L,
    refresh = 0L, seed = seed,
    laplace_control = list(multi_start_M = 2L, laplace_draws = 200L),
    verbose = FALSE
  )
}

# ---- Path B K>1 p=1 fit-code closures -----------------------------------

.fit_path_B_beta_K2 <- function(seed) {
  set.seed(seed)
  n <- 200L
  x <- rnorm(n); z <- rnorm(n)
  eta_mu <- plogis(0.4 + 0.3 * (x - mean(x)))
  y <- stats::rbeta(n, shape1 = 2 * eta_mu * 5,
                    shape2 = 2 * (1 - eta_mu) * 5)
  d <- data.frame(y = y, x = x, z = z)
  fs <- gdpar_bf(y ~ a(x), phi ~ a(z))
  gdpar_eb(
    formula = fs, family = gdpar_family("beta"), data = d,
    iter_warmup = 200L, iter_sampling = 200L, chains = 2L,
    refresh = 0L, seed = seed, skip_id_check = TRUE,
    laplace_control = list(multi_start_M = 2L, laplace_draws = 200L),
    verbose = FALSE
  )
}

.fit_path_B_gamma_K2 <- function(seed) {
  set.seed(seed)
  n <- 200L
  x <- rnorm(n); z <- rnorm(n)
  y <- stats::rgamma(n, shape = 2, rate = 1)
  d <- data.frame(y = y, x = x, z = z)
  fs <- gdpar_bf(y ~ a(x), shape ~ a(z))
  gdpar_eb(
    formula = fs, family = gdpar_family("gamma"), data = d,
    iter_warmup = 200L, iter_sampling = 200L, chains = 2L,
    refresh = 0L, seed = seed, skip_id_check = TRUE,
    laplace_control = list(multi_start_M = 2L, laplace_draws = 200L),
    verbose = FALSE
  )
}

.fit_path_B_student_t_K3 <- function(seed) {
  set.seed(seed)
  n <- 250L
  x <- rnorm(n); z <- rnorm(n); w <- rnorm(n)
  y <- 0.2 + 0.4 * (x - mean(x)) + rt(n, df = 8)
  d <- data.frame(y = y, x = x, z = z, w = w)
  fs <- gdpar_bf(y ~ a(x), sigma ~ a(z), nu ~ a(w))
  gdpar_eb(
    formula = fs, family = gdpar_family("student_t"), data = d,
    iter_warmup = 200L, iter_sampling = 200L, chains = 2L,
    refresh = 0L, seed = seed, skip_id_check = TRUE,
    laplace_control = list(multi_start_M = 2L, laplace_draws = 200L),
    verbose = FALSE
  )
}

.fit_path_B_tweedie_K3 <- function(seed) {
  rtweedie <- function(n, mu, phi, p) {
    lambda <- mu^(2 - p) / (phi * (2 - p))
    shape  <- (2 - p) / (p - 1)
    rate   <- 1 / (phi * (p - 1) * mu^(p - 1))
    N <- rpois(n, lambda)
    out <- numeric(n)
    pos <- N > 0L
    if (any(pos)) {
      out[pos] <- vapply(which(pos), function(i) {
        stats::rgamma(1L, shape = N[i] * shape, rate = rate[i])
      }, numeric(1L))
    }
    out
  }
  set.seed(seed)
  n <- 250L
  x <- rnorm(n); z <- rnorm(n); w <- rnorm(n)
  y <- rtweedie(n, mu = exp(1.0), phi = exp(0.0), p = rep(1.5, n))
  d <- data.frame(y = y, x = x, z = z, w = w)
  fs <- gdpar_bf(y ~ a(x), phi ~ a(z), p ~ a(w))
  gdpar_eb(
    formula = fs, family = gdpar_family("tweedie"), data = d,
    iter_warmup = 200L, iter_sampling = 200L, chains = 2L,
    refresh = 0L, seed = seed, skip_id_check = TRUE,
    laplace_control = list(multi_start_M = 2L, laplace_draws = 200L),
    verbose = FALSE
  )
}

.fit_path_B_lognormal_loc_scale_K2 <- function(seed) {
  set.seed(seed)
  n <- 200L
  x <- rnorm(n); z <- rnorm(n)
  y <- rlnorm(n, meanlog = 0.2, sdlog = 0.8)
  d <- data.frame(y = y, x = x, z = z)
  fam_ln <- suppressMessages(gdpar_family_custom_K(
    name = "boot_ln_8_6_C", stan_lpdf_id = "lognormal_loc_scale",
    did_holds = TRUE, did_reference = "test"
  ))
  fs <- gdpar_bf(y ~ a(x), sigma ~ a(z))
  gdpar_eb(
    formula = fs, family = fam_ln, data = d,
    iter_warmup = 200L, iter_sampling = 200L, chains = 2L,
    refresh = 0L, seed = seed, skip_id_check = TRUE,
    laplace_control = list(multi_start_M = 2L, laplace_draws = 200L),
    verbose = FALSE
  )
}

.fit_path_B_zip_K2 <- function(seed) {
  set.seed(seed)
  n <- 250L
  x <- rnorm(n); z <- rnorm(n)
  pi_vec <- plogis(-1.0 + 0.3 * z)
  is_zero <- rbinom(n, 1L, pi_vec) == 1L
  y <- ifelse(is_zero, 0L, rpois(n, lambda = exp(1.2 + 0.4 * x)))
  d <- data.frame(y = y, x = x, z = z)
  fs <- gdpar_bf(y ~ a(x), pi ~ a(z))
  gdpar_eb(
    formula = fs, family = gdpar_family("zip"), data = d,
    iter_warmup = 200L, iter_sampling = 200L, chains = 2L,
    refresh = 0L, seed = seed, skip_id_check = TRUE,
    laplace_control = list(multi_start_M = 2L, laplace_draws = 200L),
    verbose = FALSE
  )
}

.fit_path_B_het_gauss_beta_K2 <- function(seed) {
  set.seed(seed)
  n <- 250L
  x <- rnorm(n); z <- rnorm(n)
  sigma_vec <- plogis(-0.8 + 0.3 * z)
  y <- rnorm(n, mean = 0.4 + 0.4 * x, sd = sigma_vec)
  d <- data.frame(y = y, x = x, z = z)
  fs <- gdpar_bf(y ~ a(x), sigma ~ a(z))
  gdpar_eb(
    formula = fs,
    family  = list(mu = gdpar_family("gaussian"),
                   sigma = gdpar_family("beta")),
    data    = d,
    iter_warmup = 200L, iter_sampling = 200L, chains = 2L,
    refresh = 0L, seed = seed, skip_id_check = TRUE,
    laplace_control = list(multi_start_M = 2L, laplace_draws = 200L),
    verbose = FALSE
  )
}

# ---- Dispatcher: roster entry -> fit-code closure -----------------------

.dispatch_fit_code <- function(cfg) {
  key <- paste(cfg$path, cfg$family, cfg$regime, sep = "/")
  switch(
    key,
    "A/gaussian/polyP2"            = .fit_path_A_gaussian_polyP2,
    "A/gaussian/bsplineP2"         = .fit_path_A_gaussian_bsplineP2,
    "A/gaussian/polyP3"            = .fit_path_A_gaussian_polyP3,
    "A/poisson/polyP2"             = .fit_path_A_poisson_polyP2,
    "A/neg_binomial_2/polyP2"      = .fit_path_A_neg_binomial_2_polyP2,
    "A/bernoulli/polyP2"           = .fit_path_A_bernoulli_polyP2,
    "B/beta/K2"                    = .fit_path_B_beta_K2,
    "B/gamma/K2"                   = .fit_path_B_gamma_K2,
    "B/student_t/K3"               = .fit_path_B_student_t_K3,
    "B/tweedie/K3"                 = .fit_path_B_tweedie_K3,
    "B/lognormal_loc_scale/K2"     = .fit_path_B_lognormal_loc_scale_K2,
    "B/zip/K2"                     = .fit_path_B_zip_K2,
    "B/het_gauss_beta_K2/K2"       = .fit_path_B_het_gauss_beta_K2,
    stop(sprintf("No fit-code dispatch for %s", key))
  )
}

.bootstrap_one_fitable_C <- function(cfg, seed = 20260525L) {
  family <- cfg$family
  fname <- gdpar:::.gdpar_golden_eb_rds_name_8_6_C(family, cfg$path,
                                                   cfg$regime)
  rds_path_existing <- file.path(data_dir, fname)
  if (file.exists(rds_path_existing) &&
      !nzchar(Sys.getenv("GDPAR_BOOTSTRAP_EB_FORCE"))) {
    message(sprintf("[EB-8.6.C-bootstrap] %s/%s/%s already present; skipping.",
                    cfg$path, family, cfg$regime))
    return(invisible(rds_path_existing))
  }
  message(sprintf("[EB-8.6.C-bootstrap] fitting %s / %s / %s ...",
                  cfg$path, family, cfg$regime))
  fit_code_fn <- .dispatch_fit_code(cfg)
  fit_eb <- fit_code_fn(seed)
  draws <- fit_eb$conditional_fit$draws()
  draws_matrix <- posterior::as_draws_matrix(draws)
  snapshot <- list(
    family             = family,
    regime             = cfg$regime,
    path               = cfg$path,
    K                  = cfg$K,
    p                  = cfg$p,
    status             = "fitable",
    sub_phase          = "8.6.C",
    bootstrap_date     = bootstrap_date,
    cmdstan_version    = cmdstan_version_str,
    R_version          = R_version_str,
    seed               = seed,
    n_chains           = 2L,
    n_iter_warmup      = 200L,
    n_iter_sampling    = 200L,
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
      family = family, K = cfg$K, basis_type = cfg$regime, p = cfg$p,
      sub_phase = "8.6.C", bootstrap_date = bootstrap_date,
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

.bootstrap_one_guard_C <- function(cfg) {
  family <- cfg$family
  fname <- gdpar:::.gdpar_golden_eb_rds_name_8_6_C(family, cfg$path,
                                                   cfg$regime)
  message(sprintf("[EB-8.6.C-bootstrap] metadata-only guard %s / %s / %s ...",
                  cfg$path, family, cfg$regime))
  snapshot <- list(
    family          = family,
    regime          = cfg$regime,
    path            = cfg$path,
    K               = cfg$K,
    p               = cfg$p,
    status          = "guard",
    sub_phase       = "8.6.C",
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
      sub_phase = "8.6.C", bootstrap_date = bootstrap_date,
      cmdstan_version = cmdstan_version_str, R_version = R_version_str,
      DHARMa_version = gdpar:::.gdpar_dharma_version(),
      fit_code_hash = gdpar:::.gdpar_golden_fit_code_hash(.bootstrap_one_guard_C),
      seed = NA_integer_, n_chains = NA_integer_,
      n_iter_warmup = NA_integer_, n_iter_sampling = NA_integer_,
      rds_path = file.path("data", fname),
      n_obs = NA_integer_
    ),
    path = manifest_path
  )
  invisible(rds_path)
}

roster <- gdpar:::.gdpar_golden_roster_8_6_C()
for (cfg in roster) {
  if (cfg$status == "fitable") {
    .bootstrap_one_fitable_C(cfg)
  } else {
    .bootstrap_one_guard_C(cfg)
  }
}

message(sprintf("[EB-8.6.C-bootstrap] DONE. Wrote %d EB goldens + manifest rows.",
                length(roster)))

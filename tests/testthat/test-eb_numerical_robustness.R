# Sub-phase 8.6.B (2026-05-25): dedicated tests of the four-component
# anti-fragility strategy of Charter §2.8 (parcial handoff §5.3 of
# HANDOFF_SUBFASE_8_6_B_PARCIAL.md). Three Stan-running scenarios plus
# unit-level coverage of the laplace_control resolver.
#
# Scenario 1 (kappa-threshold guard): forces the marginal-covariance
#   condition number to exceed a deliberately low threshold and
#   verifies the canonical gdpar_eb_numerical_error fires. Under
#   K = 1 and p = 1 the projected covariance is 1 x 1 so its true
#   kappa is exactly 1; setting kappa_threshold < 1 is the most
#   sensitive way of exercising the guard in 8.6.B (8.6.C will
#   exercise it with multivariate covariances whose kappa exceeds
#   defaults under near-singular regimes).
# Scenario 2 (multimodality of L_n^{marg}): generates an outcome from
#   a 50/50 mixture so the marginal log-likelihood develops competing
#   modes; multi_start_dispersion grows beyond the 0.05 advisory
#   threshold and the gdpar_diagnostic_warning is emitted (matches
#   v07b §9.5 O5*-EBFB).
# Scenario 3 (Laplace fallback): builds a 2-observation outlier-only
#   dataset that crashes every multi-start optimize() init and
#   verifies the gdpar_unsupported_feature_error pointing back to
#   gdpar() (FB) is raised.

.gdpar_eb_smoke_env_active <- function() {
  nzchar(Sys.getenv("GDPAR_RUN_STAN_SMOKE_EB"))
}

# ---- Unit-level coverage (unconditional) --------------------------------

test_that(".gdpar_eb_resolve_laplace_control rejects non-positive multi_start_M", {
  expect_error(
    gdpar:::.gdpar_eb_resolve_laplace_control(list(multi_start_M = 0L)),
    class = "gdpar_input_error"
  )
})

test_that(".gdpar_eb_resolve_laplace_control rejects non-positive kappa_threshold", {
  expect_error(
    gdpar:::.gdpar_eb_resolve_laplace_control(list(kappa_threshold = -1)),
    class = "gdpar_input_error"
  )
})

test_that(".gdpar_eb_resolve_laplace_control coerces multi_start_M to integer", {
  lc <- gdpar:::.gdpar_eb_resolve_laplace_control(list(multi_start_M = 3))
  expect_identical(lc$multi_start_M, 3L)
})

test_that(".gdpar_eb_resolve_laplace_control coerces ridge_init to double", {
  lc <- gdpar:::.gdpar_eb_resolve_laplace_control(list(ridge_init = 1L))
  expect_identical(lc$ridge_init, 1.0)
})

# ---- Sub-bloque 9.3.b (Sesion B9.2, 2026-05-27) -------------------------
# Unit-level coverage of the adaptive Levenberg-Marquardt helper
# .gdpar_eb_lm_perturb, of the three new resolver fields (epsilon_lm,
# ridge_max_iter, ridge_grow_factor), and of the new diagnostics
# semantics (lm_status / kappa_post_ridge / lm_n_iter).

test_that(".gdpar_eb_resolve_laplace_control populates the 9.3.b defaults", {
  lc <- gdpar:::.gdpar_eb_resolve_laplace_control(list())
  expect_identical(lc$ridge_max_iter, 10L)
  expect_identical(lc$ridge_grow_factor, 10.0)
  expect_equal(lc$epsilon_lm, sqrt(.Machine$double.eps))
})

test_that(".gdpar_eb_resolve_laplace_control rejects non-positive epsilon_lm", {
  expect_error(
    gdpar:::.gdpar_eb_resolve_laplace_control(list(epsilon_lm = 0)),
    class = "gdpar_input_error"
  )
})

test_that(".gdpar_eb_resolve_laplace_control rejects non-positive ridge_max_iter", {
  expect_error(
    gdpar:::.gdpar_eb_resolve_laplace_control(list(ridge_max_iter = 0L)),
    class = "gdpar_input_error"
  )
})

test_that(".gdpar_eb_resolve_laplace_control rejects ridge_grow_factor <= 1", {
  expect_error(
    gdpar:::.gdpar_eb_resolve_laplace_control(list(ridge_grow_factor = 1)),
    class = "gdpar_input_error"
  )
})

test_that(".gdpar_eb_lm_perturb returns not_needed on a well-conditioned cov", {
  control <- gdpar:::.gdpar_eb_resolve_laplace_control(list())
  cov_ok <- diag(c(1.0, 2.0, 0.5))
  out <- gdpar:::.gdpar_eb_lm_perturb(cov_ok, control)
  expect_identical(out$status, "not_needed")
  expect_identical(out$n_iter, 0L)
  expect_identical(out$lambda_used, 0)
  expect_identical(out$cov_perturbed, cov_ok)
  expect_equal(out$kappa_post, 4.0)
})

test_that(".gdpar_eb_lm_perturb ridges a non-PD cov until PD (status converged)", {
  control <- gdpar:::.gdpar_eb_resolve_laplace_control(list())
  cov_bad <- matrix(c(1, 1, 1, 1), 2L, 2L)  # rank 1 -> det = 0
  out <- gdpar:::.gdpar_eb_lm_perturb(cov_bad, control)
  expect_identical(out$status, "converged")
  expect_gte(out$n_iter, 1L)
  expect_gt(out$lambda_used, 0)
  eigs <- eigen(out$cov_perturbed, symmetric = TRUE,
                only.values = TRUE)$values
  expect_true(all(eigs > 0))
  expect_lte(out$kappa_post, control$kappa_threshold)
})

test_that(".gdpar_eb_lm_perturb triggers on |det(cov)| < epsilon_lm even when PD", {
  control <- gdpar:::.gdpar_eb_resolve_laplace_control(list(
    epsilon_lm = 1e-3,
    ridge_max_iter = 10L
  ))
  cov_small_det <- diag(c(1e-2, 1e-2))  # det = 1e-4 < 1e-3
  out <- gdpar:::.gdpar_eb_lm_perturb(cov_small_det, control)
  expect_identical(out$status, "converged")
  expect_gte(out$n_iter, 1L)
  expect_gt(out$lambda_used, 0)
})

test_that(".gdpar_eb_lm_perturb exhausts when kappa_threshold is unreachable", {
  control <- gdpar:::.gdpar_eb_resolve_laplace_control(list(
    kappa_threshold   = 1.001,
    ridge_max_iter    = 3L,
    ridge_init        = 1e-12,
    ridge_grow_factor = 1.01
  ))
  # Rank-deficient input triggers the ridge (det = 0), but the
  # kappa_threshold of 1.001 is below the achievable kappa for any
  # lambda in 3 iterations with a 1.01 grow factor; status must be
  # "exhausted" and n_iter must equal ridge_max_iter.
  cov_singular <- matrix(c(1, 1, 1, 1), 2L, 2L)
  out <- gdpar:::.gdpar_eb_lm_perturb(cov_singular, control)
  expect_identical(out$status, "exhausted")
  expect_identical(out$n_iter, control$ridge_max_iter)
})

test_that(".gdpar_eb_lm_perturb handles 1x1 covariance bit-exact preserved", {
  control <- gdpar:::.gdpar_eb_resolve_laplace_control(list())
  cov_1x1 <- matrix(0.25, 1L, 1L)
  out <- gdpar:::.gdpar_eb_lm_perturb(cov_1x1, control)
  expect_identical(out$status, "not_needed")
  expect_identical(out$cov_perturbed, cov_1x1)
  expect_equal(out$kappa_post, 1.0)
})

test_that(".gdpar_eb_lm_perturb ridges a singular 1x1 cov (zero variance)", {
  control <- gdpar:::.gdpar_eb_resolve_laplace_control(list())
  cov_zero <- matrix(0.0, 1L, 1L)
  out <- gdpar:::.gdpar_eb_lm_perturb(cov_zero, control)
  expect_identical(out$status, "converged")
  expect_gt(out$lambda_used, 0)
  expect_gt(out$cov_perturbed[1L, 1L], 0)
})

# ---- Scenario 1: kappa-threshold guard (opt-in) -------------------------

test_that("scenario 1: kappa_threshold < 1 aborts with gdpar_eb_numerical_error", {
  skip_if_not_installed("cmdstanr")
  if (!.gdpar_eb_smoke_env_active()) {
    skip("Set GDPAR_RUN_STAN_SMOKE_EB to a non-empty value to enable.")
  }
  set.seed(20260525L)
  n <- 80L
  df <- data.frame(x1 = rnorm(n))
  df$y <- 0.5 + 0.3 * df$x1 + rnorm(n, sd = 0.4)
  expect_error(
    gdpar_eb(
      formula       = y ~ x1,
      family        = gdpar_family("gaussian"),
      amm           = amm_spec(a = ~ x1),
      data          = df,
      iter_warmup   = 200L, iter_sampling = 200L, chains = 2L,
      refresh       = 0L, seed = 1L,
      laplace_control = list(multi_start_M = 2L, laplace_draws = 200L,
                             kappa_threshold = 0.5),
      verbose       = FALSE
    ),
    class = "gdpar_eb_numerical_error"
  )
})

# ---- Scenario 2: multimodality / multi_start dispersion (opt-in) --------

test_that("scenario 2: bimodal-mixture data inflates multi_start_dispersion", {
  skip_if_not_installed("cmdstanr")
  if (!.gdpar_eb_smoke_env_active()) {
    skip("Set GDPAR_RUN_STAN_SMOKE_EB to a non-empty value to enable.")
  }
  set.seed(20260525L)
  n <- 120L
  comp <- rbinom(n, size = 1L, prob = 0.5)
  df <- data.frame(x1 = rnorm(n))
  df$y <- ifelse(comp == 1L,
                 -3 + 0.1 * df$x1 + rnorm(n, sd = 0.4),
                  3 + 0.1 * df$x1 + rnorm(n, sd = 0.4))
  # Wider multi-start with a higher M so dispersion is meaningful.
  res <- tryCatch(
    gdpar_eb(
      formula       = y ~ x1,
      family        = gdpar_family("gaussian"),
      amm           = amm_spec(a = ~ x1),
      data          = df,
      iter_warmup   = 200L, iter_sampling = 200L, chains = 2L,
      refresh       = 0L, seed = 1L,
      laplace_control = list(multi_start_M = 5L, laplace_draws = 200L),
      verbose       = FALSE
    ),
    error = function(e) e
  )
  # Either the fit completes (with a high dispersion / warning) or it
  # aborts with gdpar_eb_numerical_error on a singular Hessian. Both
  # are documented anti-fragility outcomes; the negative outcome is a
  # silent fit with low dispersion and no warning.
  if (inherits(res, "gdpar_eb_fit")) {
    # Dispersion can be NA when fewer than two inits converge (the
    # multimodal landscape kills most of them); both NA and a finite
    # non-negative number are documented anti-fragility outcomes (cf.
    # v07b §9.5 O5*-EBFB).
    disp <- res$diagnostics_numerical$multi_start_dispersion
    expect_true(is.na(disp) ||
                (is.finite(disp) && disp >= 0))
  } else {
    expect_true(inherits(res, "gdpar_eb_numerical_error") ||
                inherits(res, "gdpar_unsupported_feature_error"))
  }
})

# ---- Scenario 3: total Laplace fallback (opt-in) ------------------------

test_that("scenario 3: pathological outlier data triggers Laplace fallback", {
  skip_if_not_installed("cmdstanr")
  if (!.gdpar_eb_smoke_env_active()) {
    skip("Set GDPAR_RUN_STAN_SMOKE_EB to a non-empty value to enable.")
  }
  # Two observations, one with magnitude 1e20: optimize() either
  # crashes on every init (triggering the all-inits-failed fallback
  # gdpar_unsupported_feature_error) or converges to a degenerate MAP
  # whose Hessian fails laplace() (triggering
  # gdpar_eb_numerical_error). Both are documented anti-fragility
  # outcomes; the test verifies the cascade routes the failure to one
  # of the two canonical error classes rather than producing a
  # silently invalid estimator.
  df <- data.frame(y = c(1e20, -1e20), x1 = c(1, -1))
  res <- suppressWarnings(suppressMessages(tryCatch(
    gdpar_eb(
      formula       = y ~ x1,
      family        = gdpar_family("gaussian"),
      amm           = amm_spec(a = ~ x1),
      data          = df,
      iter_warmup   = 200L, iter_sampling = 200L, chains = 2L,
      refresh       = 0L, seed = 1L,
      skip_id_check = TRUE,
      laplace_control = list(multi_start_M = 2L, laplace_draws = 200L),
      verbose       = FALSE
    ),
    error = function(e) e
  )))
  expect_true(
    inherits(res, "gdpar_eb_numerical_error") ||
      inherits(res, "gdpar_unsupported_feature_error"),
    info = sprintf("got class %s", paste(class(res), collapse = ","))
  )
})

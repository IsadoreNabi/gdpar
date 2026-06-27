# Sub-phase 8.6.B (2026-05-25): first-order EB vs FB coherence test
# (Theorem 7A of v07 §5; the scalar case of Theorem 7A* of v07b §4).
#
# Under (EB-MARG-ID), (PRIOR-FB-WEAK) and (HIER-COMPLEX) with K = 1
# and p = 1, the EB plug-in estimator $\widehat\theta_{ref}^{EB}$
# concentrates around the FB posterior mean $E_{FB}[\theta_{ref}]$ at
# rate $O(n^{-1/2})$. The test fits the same Gaussian canonical model
# twice (once with gdpar() FB, once with gdpar_eb()) on a large
# synthetic sample (n = 500) with a weak prior, then verifies that
# the EB point estimate matches the FB posterior mean within
# epsilon = 5 / sqrt(n) (a slack-by-5x version of the corollary
# bound; tight enough to fail under regression but loose enough to
# survive MC variation across the seed-fixed runs).
#
# Gated by env var GDPAR_RUN_STAN_SMOKE_EB and
# skip_if_not_installed("cmdstanr"); runtime ~2 minutes (two compiled
# Stan models, each samples 200 warmup + 200 iter x 2 chains).

.gdpar_eb_smoke_env_active <- function() {
  nzchar(Sys.getenv("GDPAR_RUN_STAN_SMOKE_EB"))
}

test_that("EB plug-in matches FB posterior mean of theta_ref up to O(1/sqrt(n))", {
  skip_if_not_installed("cmdstanr")
  if (!.gdpar_eb_smoke_env_active()) {
    skip("Set GDPAR_RUN_STAN_SMOKE_EB to a non-empty value to enable.")
  }
  set.seed(20260525L)
  n <- 500L
  df <- data.frame(x1 = rnorm(n))
  df$y <- 0.7 + 0.25 * df$x1 + rnorm(n, sd = 0.5)

  fit_fb <- gdpar(
    formula       = y ~ x1,
    family        = gdpar_family("gaussian"),
    amm           = amm_spec(a = ~ x1),
    data          = df,
    iter_warmup   = 200L, iter_sampling = 200L, chains = 2L,
    refresh       = 0L, seed = 9999L,
    verbose       = FALSE
  )
  fit_eb <- gdpar_eb(
    formula       = y ~ x1,
    family        = gdpar_family("gaussian"),
    amm           = amm_spec(a = ~ x1),
    data          = df,
    iter_warmup   = 200L, iter_sampling = 200L, chains = 2L,
    refresh       = 0L, seed = 9999L,
    laplace_control = list(multi_start_M = 2L, laplace_draws = 200L),
    verbose       = FALSE
  )

  theta_ref_eb <- fit_eb$theta_ref_hat[[1L]]
  fb_draws <- posterior::as_draws_matrix(fit_fb$fit$draws(variables = "theta_ref"))
  theta_ref_fb_mean <- mean(fb_draws[, 1L])
  epsilon <- 5 / sqrt(n)  # 5x slack on the asymptotic 1 / sqrt(n) bound
  diff <- abs(theta_ref_eb - theta_ref_fb_mean)
  expect_true(
    diff < epsilon,
    info = sprintf(
      "EB - FB |diff| = %.4f exceeds epsilon = %.4f (5 / sqrt(n))",
      diff, epsilon
    )
  )
})

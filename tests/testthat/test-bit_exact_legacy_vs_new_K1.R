# Sub-phase 8.3.3 Unit 3 (2026-05-20): bit-exact regression net before
# the K > 1 codegen lands. The new dispatch in gdpar() unwraps K = 1
# inputs (gdpar_formula_set, named-list amm_spec, classic formula with
# AMM wrappers) to the legacy K = 1 path. This file asserts that the
# legacy invocation
#
#   gdpar(y ~ x, amm = amm_spec(a = ~ x), ...)
#
# and the new-API invocation
#
#   gdpar(gdpar_bf(y ~ a(x)), ...)
#
# produce draws bit-identical to one another (same Stan source, same
# seed, same data, same sampler config). The test is gated by
# GDPAR_RUN_STAN_SMOKE because it requires cmdstan compilation and
# adds 10-20 seconds to the suite when active.

.skip_if_no_stan_smoke <- function() {
  if (!nzchar(Sys.getenv("GDPAR_RUN_STAN_SMOKE"))) {
    testthat::skip("GDPAR_RUN_STAN_SMOKE not set; skipping Stan compilation")
  }
  testthat::skip_if_not_installed("cmdstanr")
}

test_that("legacy K=1 path and new-API K=1 path produce bit-identical draws", {
  .skip_if_no_stan_smoke()
  set.seed(101L)
  n <- 30L
  d <- data.frame(y = rnorm(n), x = rnorm(n))
  cfg <- list(chains = 1L, iter_warmup = 200L, iter_sampling = 200L,
              seed = 1234L, refresh = 0L, verbose = FALSE,
              skip_id_check = TRUE)
  fit_legacy <- do.call(gdpar, c(
    list(formula = y ~ x, amm = amm_spec(a = ~ x), data = d,
         family = gdpar_family("gaussian")),
    cfg
  ))
  fit_new <- do.call(gdpar, c(
    list(formula = gdpar_bf(y ~ a(x)), data = d,
         family = gdpar_family("gaussian")),
    cfg
  ))
  expect_equal(fit_legacy$fit$draws(),
               fit_new$fit$draws())
})

test_that("legacy K=1 path and named-list K=1 path produce bit-identical draws", {
  .skip_if_no_stan_smoke()
  set.seed(202L)
  n <- 30L
  d <- data.frame(y = rnorm(n), x = rnorm(n))
  cfg <- list(chains = 1L, iter_warmup = 200L, iter_sampling = 200L,
              seed = 7777L, refresh = 0L, verbose = FALSE,
              skip_id_check = TRUE)
  fit_legacy <- do.call(gdpar, c(
    list(formula = y ~ x, amm = amm_spec(a = ~ x), data = d,
         family = gdpar_family("gaussian")),
    cfg
  ))
  fit_list <- do.call(gdpar, c(
    list(formula = y ~ x,
         amm     = list(mu = amm_spec(a = ~ x)),
         data    = d,
         family  = gdpar_family("gaussian")),
    cfg
  ))
  expect_equal(fit_legacy$fit$draws(),
               fit_list$fit$draws())
})

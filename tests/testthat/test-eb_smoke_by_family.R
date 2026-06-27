# Tests for gdpar_eb (Sub-phase 8.6.B closure session, 2026-05-25):
# Charter §3.2 / parcial handoff §5.1 — full 17-case smoke matrix.
#
# The matrix exercises all family / setup combinations canonized for
# Sub-phase 8.6.B and verifies, for each combination, that the
# orchestrator either fits cleanly (4 fitable univariate families +
# B-spline W variant of Gaussian) or aborts with the canonical
# gdpar_unsupported_feature_error class (12 deferred combinations:
# 5 canonically K >= 2 families per Decision D33 of this session,
# 4 mixtures / hurdle families per family$min_K guard, 3 heterogeneous
# K > 1 specifications per the amm_spec class guard).
#
# Fitable smokes are gated by env var GDPAR_RUN_STAN_SMOKE_EB and
# skip_if_not_installed("cmdstanr"), consistent with the per-feature
# opt-in pattern of the suite (matches test-causal_bridge.R, the
# 8.5.A precedent). Guard-reject smokes run unconditionally and need
# no Stan toolchain.

.gdpar_eb_smoke_env_active <- function() {
  nzchar(Sys.getenv("GDPAR_RUN_STAN_SMOKE_EB"))
}

# ---- Fitable smokes (env-gated) -----------------------------------------

test_that("smoke 1/17: gaussian K=1 polynomial fits end-to-end (opt-in)", {
  skip_if_not_installed("cmdstanr")
  if (!.gdpar_eb_smoke_env_active()) {
    skip("Set GDPAR_RUN_STAN_SMOKE_EB to a non-empty value to enable.")
  }
  set.seed(20260525L)
  n <- 80L
  df <- data.frame(x1 = rnorm(n))
  df$y <- 0.5 + 0.3 * df$x1 + rnorm(n, sd = 0.4)
  fit_eb <- gdpar_eb(
    formula       = y ~ x1,
    family        = gdpar_family("gaussian"),
    amm           = amm_spec(a = ~ x1),
    data          = df,
    iter_warmup   = 200L, iter_sampling = 200L, chains = 2L,
    refresh       = 0L, seed = 1001L,
    laplace_control = list(multi_start_M = 2L, laplace_draws = 200L),
    verbose       = FALSE
  )
  expect_s3_class(fit_eb, "gdpar_eb_fit")
  expect_true(isTRUE(is.finite(fit_eb$theta_ref_hat)))
  expect_true(isTRUE(is.finite(fit_eb$diagnostics_numerical$kappa)))
  expect_identical(fit_eb$family$name, "gaussian")
  expect_true(fit_eb$correction_applied)
})

test_that("smoke 2/17: poisson K=1 log link fits end-to-end (opt-in)", {
  skip_if_not_installed("cmdstanr")
  if (!.gdpar_eb_smoke_env_active()) {
    skip("Set GDPAR_RUN_STAN_SMOKE_EB to a non-empty value to enable.")
  }
  set.seed(20260525L)
  n <- 80L
  df <- data.frame(x1 = rnorm(n))
  df$y <- rpois(n, lambda = exp(0.6 + 0.2 * df$x1))
  fit_eb <- gdpar_eb(
    formula       = y ~ x1,
    family        = gdpar_family("poisson"),
    amm           = amm_spec(a = ~ x1),
    data          = df,
    iter_warmup   = 200L, iter_sampling = 200L, chains = 2L,
    refresh       = 0L, seed = 1002L,
    laplace_control = list(multi_start_M = 2L, laplace_draws = 200L),
    verbose       = FALSE
  )
  expect_s3_class(fit_eb, "gdpar_eb_fit")
  expect_true(isTRUE(is.finite(fit_eb$theta_ref_hat)))
  expect_identical(fit_eb$family$name, "poisson")
})

test_that("smoke 3/17: neg_binomial_2 K=1 fits end-to-end (opt-in)", {
  skip_if_not_installed("cmdstanr")
  if (!.gdpar_eb_smoke_env_active()) {
    skip("Set GDPAR_RUN_STAN_SMOKE_EB to a non-empty value to enable.")
  }
  set.seed(20260525L)
  n <- 80L
  df <- data.frame(x1 = rnorm(n))
  mu <- exp(0.5 + 0.2 * df$x1)
  df$y <- rnbinom(n, size = 5, mu = mu)
  fit_eb <- gdpar_eb(
    formula       = y ~ x1,
    family        = gdpar_family("neg_binomial_2"),
    amm           = amm_spec(a = ~ x1),
    data          = df,
    iter_warmup   = 200L, iter_sampling = 200L, chains = 2L,
    refresh       = 0L, seed = 1003L,
    laplace_control = list(multi_start_M = 2L, laplace_draws = 200L),
    verbose       = FALSE
  )
  expect_s3_class(fit_eb, "gdpar_eb_fit")
  expect_true(isTRUE(is.finite(fit_eb$theta_ref_hat)))
  expect_identical(fit_eb$family$name, "neg_binomial_2")
})

test_that("smoke 4/17: bernoulli K=1 logit link fits end-to-end (opt-in)", {
  skip_if_not_installed("cmdstanr")
  if (!.gdpar_eb_smoke_env_active()) {
    skip("Set GDPAR_RUN_STAN_SMOKE_EB to a non-empty value to enable.")
  }
  set.seed(20260525L)
  n <- 120L
  df <- data.frame(x1 = rnorm(n))
  p <- plogis(0.3 + 0.4 * df$x1)
  df$y <- rbinom(n, size = 1L, prob = p)
  fit_eb <- gdpar_eb(
    formula       = y ~ x1,
    family        = gdpar_family("bernoulli"),
    amm           = amm_spec(a = ~ x1),
    data          = df,
    iter_warmup   = 200L, iter_sampling = 200L, chains = 2L,
    refresh       = 0L, seed = 1004L,
    laplace_control = list(multi_start_M = 2L, laplace_draws = 200L),
    verbose       = FALSE
  )
  expect_s3_class(fit_eb, "gdpar_eb_fit")
  expect_true(isTRUE(is.finite(fit_eb$theta_ref_hat)))
  expect_identical(fit_eb$family$name, "bernoulli")
})

test_that("smoke 5/17: gaussian K=1 with B-spline W fits end-to-end (opt-in)", {
  skip_if_not_installed("cmdstanr")
  if (!.gdpar_eb_smoke_env_active()) {
    skip("Set GDPAR_RUN_STAN_SMOKE_EB to a non-empty value to enable.")
  }
  set.seed(20260525L)
  n <- 100L
  df <- data.frame(x1 = rnorm(n))
  df$y <- 0.5 + 0.3 * df$x1 + rnorm(n, sd = 0.4)
  Wspec <- W_basis(type = "bspline", degree = 3L,
                   boundary_knots = c(-3, 3), df = 4L)
  amm_W <- amm_spec(a = ~ x1, W = Wspec)
  fit_eb <- gdpar_eb(
    formula       = y ~ x1,
    family        = gdpar_family("gaussian"),
    amm           = amm_W,
    data          = df,
    iter_warmup   = 200L, iter_sampling = 200L, chains = 2L,
    refresh       = 0L, seed = 1005L,
    laplace_control = list(multi_start_M = 2L, laplace_draws = 200L),
    verbose       = FALSE
  )
  expect_s3_class(fit_eb, "gdpar_eb_fit")
  expect_true(isTRUE(is.finite(fit_eb$theta_ref_hat)))
})

# ---- Guard-reject smokes (unconditional, no Stan toolchain needed) ------

# Decision D33 (2026-05-25): family$stan_id outside {1,2,3,4} rejected.

test_that("smoke 6/17: beta family rejected by stan_id guard (D33)", {
  expect_error(
    gdpar_eb(
      formula = y ~ x1,
      family  = gdpar_family("beta"),
      amm     = amm_spec(a = ~ x1),
      data    = data.frame(y = stats::rbeta(20, 2, 5), x1 = rnorm(20))
    ),
    class = "gdpar_unsupported_feature_error"
  )
})

test_that("smoke 7/17: gamma family rejected by stan_id guard (D33)", {
  expect_error(
    gdpar_eb(
      formula = y ~ x1,
      family  = gdpar_family("gamma"),
      amm     = amm_spec(a = ~ x1),
      data    = data.frame(y = stats::rgamma(20, 2, 1), x1 = rnorm(20))
    ),
    class = "gdpar_unsupported_feature_error"
  )
})

test_that("smoke 8/17: student_t family rejected by stan_id guard (D33)", {
  expect_error(
    gdpar_eb(
      formula = y ~ x1,
      family  = gdpar_family("student_t"),
      amm     = amm_spec(a = ~ x1),
      data    = data.frame(y = rt(20, df = 5), x1 = rnorm(20))
    ),
    class = "gdpar_unsupported_feature_error"
  )
})

test_that("smoke 9/17: tweedie family rejected by stan_id guard (D33)", {
  expect_error(
    gdpar_eb(
      formula = y ~ x1,
      family  = gdpar_family("tweedie"),
      amm     = amm_spec(a = ~ x1),
      data    = data.frame(y = abs(rnorm(20)) + 0.01, x1 = rnorm(20))
    ),
    class = "gdpar_unsupported_feature_error"
  )
})

test_that("smoke 10/17: lognormal_loc_scale not registered as K=1 family", {
  # lognormal-loc-scale lives only at K >= 2 via gdpar_family_custom_K()
  # and is therefore unreachable through gdpar_family() in 8.6.B; the
  # match.arg() set rejects the name. We exercise the failure mode that
  # an 8.6.B user would actually hit when trying to wire lognormal-EB.
  expect_error(
    gdpar_family("lognormal_loc_scale")
  )
})

test_that("smoke 11/17: ZIP rejected by family$min_K guard", {
  zip_fake <- list(name = "zip", link = "log", K = 1L, min_K = 2L,
                   stan_id = 10L)
  class(zip_fake) <- c("gdpar_family", "list")
  expect_error(
    gdpar_eb(
      formula = y ~ x1,
      family  = zip_fake,
      amm     = amm_spec(a = ~ x1),
      data    = data.frame(y = rpois(20, 2), x1 = rnorm(20))
    ),
    class = "gdpar_unsupported_feature_error"
  )
})

test_that("smoke 12/17: ZINB rejected by family$min_K guard", {
  zinb_fake <- list(name = "zinb", link = "log", K = 1L, min_K = 3L,
                    stan_id = 11L)
  class(zinb_fake) <- c("gdpar_family", "list")
  expect_error(
    gdpar_eb(
      formula = y ~ x1,
      family  = zinb_fake,
      amm     = amm_spec(a = ~ x1),
      data    = data.frame(y = rpois(20, 2), x1 = rnorm(20))
    ),
    class = "gdpar_unsupported_feature_error"
  )
})

test_that("smoke 13/17: hurdle_poisson rejected by family$min_K guard", {
  hp_fake <- list(name = "hurdle_poisson", link = "log",
                  K = 1L, min_K = 2L, stan_id = 12L)
  class(hp_fake) <- c("gdpar_family", "list")
  expect_error(
    gdpar_eb(
      formula = y ~ x1,
      family  = hp_fake,
      amm     = amm_spec(a = ~ x1),
      data    = data.frame(y = rpois(20, 2), x1 = rnorm(20))
    ),
    class = "gdpar_unsupported_feature_error"
  )
})

test_that("smoke 14/17: hurdle_neg_binomial_2 rejected by family$min_K guard", {
  hnb_fake <- list(name = "hurdle_neg_binomial_2", link = "log",
                   K = 1L, min_K = 3L, stan_id = 13L)
  class(hnb_fake) <- c("gdpar_family", "list")
  expect_error(
    gdpar_eb(
      formula = y ~ x1,
      family  = hnb_fake,
      amm     = amm_spec(a = ~ x1),
      data    = data.frame(y = rpois(20, 2), x1 = rnorm(20))
    ),
    class = "gdpar_unsupported_feature_error"
  )
})

# Heterogeneous K > 1 specs under Sub-phase 8.6.C (decision D34):
# named-list amm input is no longer rejected by a scope guard. The
# rejection now comes from the family-promotion layer when the slots
# named in the amm list (phi, p_power, nu, ...) are not eligible in
# the supplied family (here gaussian, whose param_specs are mu/sigma).
# The error class changes from gdpar_unsupported_feature_error (8.6.B
# scope guard) to gdpar_input_error (semantic mismatch between amm
# slot names and family eligibility).

test_that("smoke 15/17: heterogeneous amm with phi slot rejected for gaussian family (D34)", {
  het_amm <- list(mu = amm_spec(a = ~ x1), phi = amm_spec(a = ~ x1))
  expect_error(
    gdpar_eb(
      formula = y ~ x1,
      family  = gdpar_family("gaussian"),
      amm     = het_amm,
      data    = data.frame(y = rnorm(20), x1 = rnorm(20))
    ),
    class = "gdpar_input_error"
  )
})

test_that("smoke 16/17: heterogeneous amm with p_power slot rejected for gaussian family (D34)", {
  het_amm <- list(mu = amm_spec(a = ~ x1),
                  phi = amm_spec(a = ~ x1),
                  p_power = amm_spec(a = ~ x1))
  expect_error(
    gdpar_eb(
      formula = y ~ x1,
      family  = gdpar_family("gaussian"),
      amm     = het_amm,
      data    = data.frame(y = rnorm(20), x1 = rnorm(20))
    ),
    class = "gdpar_input_error"
  )
})

test_that("smoke 17/17: heterogeneous amm with nu slot rejected for gaussian family (D34)", {
  het_amm <- list(mu = amm_spec(a = ~ x1),
                  sigma = amm_spec(a = ~ x1),
                  nu = amm_spec(a = ~ x1))
  expect_error(
    gdpar_eb(
      formula = y ~ x1,
      family  = gdpar_family("gaussian"),
      amm     = het_amm,
      data    = data.frame(y = rnorm(20), x1 = rnorm(20))
    ),
    class = "gdpar_input_error"
  )
})

# ===========================================================================
# Sub-phase 8.6.C smokes (D34 Path A K=1+p>1 + Path B K>1+p=1)
# ===========================================================================
# Smokes 18-27 cover Path A (K=1+p>1) via amm_eb_marginal_multi.stan and
# amm_eb_conditional_multi.stan; the supported stan_id set for Path A
# inherits amm_distrib_multi.stan (Gaussian/Poisson/NB/Bernoulli only),
# so the matched FB family scope of {1,2,3,4} carries over to Path A.
# Smokes 28-34 cover Path B (K>1+p=1) via amm_eb_marginal_K.stan and
# amm_eb_conditional_K.stan; the supported stan_id set for Path B
# inherits amm_distrib_K.stan {1, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13}.

# ---- Path A K=1 p>1 fitable smokes (env-gated) --------------------------

test_that("smoke 18/34: Path A gaussian p=2 polynomial fits (D34, opt-in)", {
  skip_if_not_installed("cmdstanr")
  if (!.gdpar_eb_smoke_env_active()) {
    skip("Set GDPAR_RUN_STAN_SMOKE_EB to a non-empty value to enable.")
  }
  set.seed(20260525L)
  n <- 80L
  df <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  df$y <- matrix(rnorm(n * 2L, sd = 0.4), n, 2L)
  spec <- amm_spec(p = 2L, dims = dimwise(a = ~ x1 + x2, b = NULL))
  fit_eb <- gdpar_eb(
    formula = y ~ x1 + x2,
    family  = gdpar_family_multi("gaussian", p = 2L),
    amm     = spec, data = df,
    iter_warmup = 200L, iter_sampling = 200L, chains = 2L,
    refresh = 0L, seed = 2001L,
    laplace_control = list(multi_start_M = 2L, laplace_draws = 200L),
    verbose = FALSE
  )
  expect_s3_class(fit_eb, "gdpar_eb_fit")
  expect_true(all(is.finite(fit_eb$theta_ref_hat)))
})

test_that("smoke 19/34: Path A gaussian p=2 B-spline fits (D34, opt-in)", {
  skip_if_not_installed("cmdstanr")
  if (!.gdpar_eb_smoke_env_active()) {
    skip("Set GDPAR_RUN_STAN_SMOKE_EB to a non-empty value to enable.")
  }
  set.seed(20260525L)
  n <- 100L
  df <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  df$y <- matrix(rnorm(n * 2L, sd = 0.4), n, 2L)
  Wspec <- W_basis(type = "bspline", degree = 3L,
                   boundary_knots = c(-3, 3), df = 4L, p = 2L)
  spec <- amm_spec(p = 2L, dims = dimwise(a = ~ x1 + x2, b = NULL),
                   W = Wspec)
  fit_eb <- gdpar_eb(
    formula = y ~ x1 + x2,
    family  = gdpar_family_multi("gaussian", p = 2L),
    amm     = spec, data = df,
    iter_warmup = 200L, iter_sampling = 200L, chains = 2L,
    refresh = 0L, seed = 2002L, skip_id_check = TRUE,
    laplace_control = list(multi_start_M = 2L, laplace_draws = 200L),
    verbose = FALSE
  )
  expect_s3_class(fit_eb, "gdpar_eb_fit")
})

test_that("smoke 20/34: Path A gaussian p=3 polynomial fits (D34, opt-in)", {
  skip_if_not_installed("cmdstanr")
  if (!.gdpar_eb_smoke_env_active()) {
    skip("Set GDPAR_RUN_STAN_SMOKE_EB to a non-empty value to enable.")
  }
  set.seed(20260525L)
  n <- 80L
  df <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  df$y <- matrix(rnorm(n * 3L, sd = 0.4), n, 3L)
  spec <- amm_spec(p = 3L,
                   dims = dimwise(a = ~ x1 + x2 + x3, b = NULL))
  fit_eb <- gdpar_eb(
    formula = y ~ x1 + x2 + x3,
    family  = gdpar_family_multi("gaussian", p = 3L),
    amm     = spec, data = df,
    iter_warmup = 200L, iter_sampling = 200L, chains = 2L,
    refresh = 0L, seed = 2003L,
    laplace_control = list(multi_start_M = 2L, laplace_draws = 200L),
    verbose = FALSE
  )
  expect_s3_class(fit_eb, "gdpar_eb_fit")
})

test_that("smoke 21/34: Path A poisson p=2 polynomial fits (D34, opt-in)", {
  skip_if_not_installed("cmdstanr")
  if (!.gdpar_eb_smoke_env_active()) {
    skip("Set GDPAR_RUN_STAN_SMOKE_EB to a non-empty value to enable.")
  }
  set.seed(20260525L)
  n <- 80L
  df <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  df$y <- matrix(rpois(n * 2L, lambda = exp(0.5)), n, 2L)
  spec <- amm_spec(p = 2L, dims = dimwise(a = ~ x1 + x2, b = NULL))
  fit_eb <- gdpar_eb(
    formula = y ~ x1 + x2,
    family  = gdpar_family_multi("poisson", p = 2L),
    amm     = spec, data = df,
    iter_warmup = 200L, iter_sampling = 200L, chains = 2L,
    refresh = 0L, seed = 2004L,
    laplace_control = list(multi_start_M = 2L, laplace_draws = 200L),
    verbose = FALSE
  )
  expect_s3_class(fit_eb, "gdpar_eb_fit")
})

test_that("smoke 22/34: Path A neg_binomial_2 p=2 fits (D34, opt-in)", {
  skip_if_not_installed("cmdstanr")
  if (!.gdpar_eb_smoke_env_active()) {
    skip("Set GDPAR_RUN_STAN_SMOKE_EB to a non-empty value to enable.")
  }
  set.seed(20260525L)
  n <- 80L
  df <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  df$y <- matrix(rnbinom(n * 2L, size = 5, mu = exp(0.5)), n, 2L)
  spec <- amm_spec(p = 2L, dims = dimwise(a = ~ x1 + x2, b = NULL))
  fit_eb <- gdpar_eb(
    formula = y ~ x1 + x2,
    family  = gdpar_family_multi("neg_binomial_2", p = 2L),
    amm     = spec, data = df,
    iter_warmup = 200L, iter_sampling = 200L, chains = 2L,
    refresh = 0L, seed = 2005L,
    laplace_control = list(multi_start_M = 2L, laplace_draws = 200L),
    verbose = FALSE
  )
  expect_s3_class(fit_eb, "gdpar_eb_fit")
})

test_that("smoke 23/34: Path A bernoulli p=2 logit link fits (D34, opt-in)", {
  skip_if_not_installed("cmdstanr")
  if (!.gdpar_eb_smoke_env_active()) {
    skip("Set GDPAR_RUN_STAN_SMOKE_EB to a non-empty value to enable.")
  }
  set.seed(20260525L)
  n <- 120L
  df <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  df$y <- matrix(rbinom(n * 2L, 1L, prob = plogis(0.2)), n, 2L)
  spec <- amm_spec(p = 2L, dims = dimwise(a = ~ x1 + x2, b = NULL))
  fit_eb <- gdpar_eb(
    formula = y ~ x1 + x2,
    family  = gdpar_family_multi("bernoulli", p = 2L),
    amm     = spec, data = df,
    iter_warmup = 200L, iter_sampling = 200L, chains = 2L,
    refresh = 0L, seed = 2006L,
    laplace_control = list(multi_start_M = 2L, laplace_draws = 200L),
    verbose = FALSE
  )
  expect_s3_class(fit_eb, "gdpar_eb_fit")
})

# ---- Path A K=1 p>1 guard smokes (D34 Path A supported {1,2,3,4}) -------

test_that("smoke 24/34: Path A beta p=2 rejected by check_stan_id_for_path", {
  spec <- amm_spec(p = 2L, dims = dimwise(a = ~ x1 + x2, b = NULL))
  df <- data.frame(x1 = rnorm(20), x2 = rnorm(20),
                   y = matrix(stats::rbeta(40, 2, 5), 20, 2L))
  fam <- gdpar_family_multi("beta", p = 2L)
  expect_error(
    gdpar_eb(formula = y ~ x1 + x2, family = fam,
             amm = spec, data = df, verbose = FALSE),
    class = "gdpar_unsupported_feature_error"
  )
})

test_that("smoke 25/34: Path A gamma p=2 rejected by check_stan_id_for_path", {
  spec <- amm_spec(p = 2L, dims = dimwise(a = ~ x1 + x2, b = NULL))
  df <- data.frame(x1 = rnorm(20), x2 = rnorm(20),
                   y = matrix(stats::rgamma(40, 2, 1), 20, 2L))
  fam <- gdpar_family_multi("gamma", p = 2L)
  expect_error(
    gdpar_eb(formula = y ~ x1 + x2, family = fam,
             amm = spec, data = df, verbose = FALSE),
    class = "gdpar_unsupported_feature_error"
  )
})

test_that("smoke 26/34: Path A student_t p=2 rejected by check_stan_id_for_path", {
  spec <- amm_spec(p = 2L, dims = dimwise(a = ~ x1 + x2, b = NULL))
  df <- data.frame(x1 = rnorm(20), x2 = rnorm(20),
                   y = matrix(rt(40, df = 5), 20, 2L))
  fam <- gdpar_family_multi("student_t", p = 2L)
  expect_error(
    gdpar_eb(formula = y ~ x1 + x2, family = fam,
             amm = spec, data = df, verbose = FALSE),
    class = "gdpar_unsupported_feature_error"
  )
})

test_that("smoke 27/34: Path A tweedie p=2 rejected by check_stan_id_for_path", {
  spec <- amm_spec(p = 2L, dims = dimwise(a = ~ x1 + x2, b = NULL))
  df <- data.frame(x1 = rnorm(20), x2 = rnorm(20),
                   y = matrix(abs(rnorm(40)) + 0.01, 20, 2L))
  fam <- gdpar_family_multi("tweedie", p = 2L)
  expect_error(
    gdpar_eb(formula = y ~ x1 + x2, family = fam,
             amm = spec, data = df, verbose = FALSE),
    class = "gdpar_unsupported_feature_error"
  )
})

# ---- Path B K>1 p=1 fitable smokes (env-gated) --------------------------

test_that("smoke 28/34: Path B beta K=2 polynomial fits (D34, opt-in)", {
  skip_if_not_installed("cmdstanr")
  if (!.gdpar_eb_smoke_env_active()) {
    skip("Set GDPAR_RUN_STAN_SMOKE_EB to a non-empty value to enable.")
  }
  set.seed(20260525L)
  n <- 200L
  x <- rnorm(n); z <- rnorm(n)
  eta_mu <- plogis(0.4 + 0.3 * (x - mean(x)))
  y <- stats::rbeta(n, shape1 = 2 * eta_mu * 5,
                    shape2 = 2 * (1 - eta_mu) * 5)
  d <- data.frame(y = y, x = x, z = z)
  fs <- gdpar_bf(y ~ a(x), phi ~ a(z))
  fit_eb <- gdpar_eb(
    formula = fs, family = gdpar_family("beta"), data = d,
    iter_warmup = 200L, iter_sampling = 200L, chains = 2L,
    refresh = 0L, seed = 2007L, skip_id_check = TRUE,
    laplace_control = list(multi_start_M = 2L, laplace_draws = 200L),
    verbose = FALSE
  )
  expect_s3_class(fit_eb, "gdpar_eb_fit")
  expect_equal(fit_eb$K, 2L)
})

test_that("smoke 29/34: Path B gamma K=2 polynomial fits (D34, opt-in)", {
  skip_if_not_installed("cmdstanr")
  if (!.gdpar_eb_smoke_env_active()) {
    skip("Set GDPAR_RUN_STAN_SMOKE_EB to a non-empty value to enable.")
  }
  set.seed(20260525L)
  n <- 200L
  x <- rnorm(n); z <- rnorm(n)
  y <- stats::rgamma(n, shape = 2, rate = 1)
  d <- data.frame(y = y, x = x, z = z)
  fs <- gdpar_bf(y ~ a(x), shape ~ a(z))
  fit_eb <- gdpar_eb(
    formula = fs, family = gdpar_family("gamma"), data = d,
    iter_warmup = 200L, iter_sampling = 200L, chains = 2L,
    refresh = 0L, seed = 2008L, skip_id_check = TRUE,
    laplace_control = list(multi_start_M = 2L, laplace_draws = 200L),
    verbose = FALSE
  )
  expect_s3_class(fit_eb, "gdpar_eb_fit")
  expect_equal(fit_eb$K, 2L)
})

test_that("smoke 30/34: Path B student_t K=3 polynomial fits (D34, opt-in)", {
  skip_if_not_installed("cmdstanr")
  if (!.gdpar_eb_smoke_env_active()) {
    skip("Set GDPAR_RUN_STAN_SMOKE_EB to a non-empty value to enable.")
  }
  set.seed(20260525L)
  n <- 250L
  x <- rnorm(n); z <- rnorm(n); w <- rnorm(n)
  y <- 0.2 + 0.4 * (x - mean(x)) + rt(n, df = 8)
  d <- data.frame(y = y, x = x, z = z, w = w)
  fs <- gdpar_bf(y ~ a(x), sigma ~ a(z), nu ~ a(w))
  fit_eb <- gdpar_eb(
    formula = fs, family = gdpar_family("student_t"), data = d,
    iter_warmup = 200L, iter_sampling = 200L, chains = 2L,
    refresh = 0L, seed = 2009L, skip_id_check = TRUE,
    laplace_control = list(multi_start_M = 2L, laplace_draws = 200L),
    verbose = FALSE
  )
  expect_s3_class(fit_eb, "gdpar_eb_fit")
  expect_equal(fit_eb$K, 3L)
})

test_that("smoke 31/34: Path B tweedie K=3 polynomial fits (D34, opt-in)", {
  skip_if_not_installed("cmdstanr")
  if (!.gdpar_eb_smoke_env_active()) {
    skip("Set GDPAR_RUN_STAN_SMOKE_EB to a non-empty value to enable.")
  }
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
  set.seed(20260525L)
  n <- 250L
  x <- rnorm(n); z <- rnorm(n); w <- rnorm(n)
  y <- rtweedie(n, mu = exp(1.0), phi = exp(0.0), p = rep(1.5, n))
  d <- data.frame(y = y, x = x, z = z, w = w)
  fs <- gdpar_bf(y ~ a(x), phi ~ a(z), p ~ a(w))
  fit_eb <- gdpar_eb(
    formula = fs, family = gdpar_family("tweedie"), data = d,
    iter_warmup = 200L, iter_sampling = 200L, chains = 2L,
    refresh = 0L, seed = 2010L, skip_id_check = TRUE,
    laplace_control = list(multi_start_M = 2L, laplace_draws = 200L),
    verbose = FALSE
  )
  expect_s3_class(fit_eb, "gdpar_eb_fit")
  expect_equal(fit_eb$K, 3L)
})

test_that("smoke 32/34: Path B lognormal_loc_scale K=2 polynomial fits (D34, opt-in)", {
  skip_if_not_installed("cmdstanr")
  if (!.gdpar_eb_smoke_env_active()) {
    skip("Set GDPAR_RUN_STAN_SMOKE_EB to a non-empty value to enable.")
  }
  set.seed(20260525L)
  n <- 200L
  x <- rnorm(n); z <- rnorm(n)
  y <- rlnorm(n, meanlog = 0.2, sdlog = 0.8)
  d <- data.frame(y = y, x = x, z = z)
  fam_ln <- suppressMessages(gdpar_family_custom_K(
    name = "smoke_ln_8_6_C", stan_lpdf_id = "lognormal_loc_scale",
    did_holds = TRUE, did_reference = "test"
  ))
  fs <- gdpar_bf(y ~ a(x), sigma ~ a(z))
  fit_eb <- gdpar_eb(
    formula = fs, family = fam_ln, data = d,
    iter_warmup = 200L, iter_sampling = 200L, chains = 2L,
    refresh = 0L, seed = 2011L, skip_id_check = TRUE,
    laplace_control = list(multi_start_M = 2L, laplace_draws = 200L),
    verbose = FALSE
  )
  expect_s3_class(fit_eb, "gdpar_eb_fit")
  expect_equal(fit_eb$K, 2L)
})

test_that("smoke 33/34: Path B zip K=2 polynomial fits (D34, opt-in)", {
  skip_if_not_installed("cmdstanr")
  if (!.gdpar_eb_smoke_env_active()) {
    skip("Set GDPAR_RUN_STAN_SMOKE_EB to a non-empty value to enable.")
  }
  set.seed(20260525L)
  n <- 250L
  x <- rnorm(n); z <- rnorm(n)
  pi_vec <- plogis(-1.0 + 0.3 * z)
  is_zero <- rbinom(n, 1L, pi_vec) == 1L
  y <- ifelse(is_zero, 0L, rpois(n, lambda = exp(1.2 + 0.4 * x)))
  d <- data.frame(y = y, x = x, z = z)
  fs <- gdpar_bf(y ~ a(x), pi ~ a(z))
  fit_eb <- gdpar_eb(
    formula = fs, family = gdpar_family("zip"), data = d,
    iter_warmup = 200L, iter_sampling = 200L, chains = 2L,
    refresh = 0L, seed = 2012L, skip_id_check = TRUE,
    laplace_control = list(multi_start_M = 2L, laplace_draws = 200L),
    verbose = FALSE
  )
  expect_s3_class(fit_eb, "gdpar_eb_fit")
  expect_equal(fit_eb$K, 2L)
})

test_that("smoke 34/34: Path B heterogeneous gauss+beta K=2 fits (D34, opt-in)", {
  skip_if_not_installed("cmdstanr")
  if (!.gdpar_eb_smoke_env_active()) {
    skip("Set GDPAR_RUN_STAN_SMOKE_EB to a non-empty value to enable.")
  }
  set.seed(20260525L)
  n <- 250L
  x <- rnorm(n); z <- rnorm(n)
  sigma_vec <- plogis(-0.8 + 0.3 * z)
  y <- rnorm(n, mean = 0.4 + 0.4 * x, sd = sigma_vec)
  d <- data.frame(y = y, x = x, z = z)
  fs <- gdpar_bf(y ~ a(x), sigma ~ a(z))
  fit_eb <- gdpar_eb(
    formula = fs,
    family  = list(mu = gdpar_family("gaussian"),
                   sigma = gdpar_family("beta")),
    data    = d,
    iter_warmup = 200L, iter_sampling = 200L, chains = 2L,
    refresh = 0L, seed = 2013L, skip_id_check = TRUE,
    laplace_control = list(multi_start_M = 2L, laplace_draws = 200L),
    verbose = FALSE
  )
  expect_s3_class(fit_eb, "gdpar_eb_fit")
  expect_equal(fit_eb$K, 2L)
})

# -----------------------------------------------------------------------------
# Sub-phase 8.6.D Path C K x p smokes (35..42)
# -----------------------------------------------------------------------------
# Coverage initial: D40' = stan_id in {1, 3} (Gaussian K=2 + NB K=2),
# crossed with p in {2, 3}. Remaining Path B set deferred per opening
# Section 6.1 numerical caveat -- guarded entries below validate the
# guard message rather than fit end-to-end.

test_that("smoke 35/42: Path C Gaussian K=2 p=2 polynomial fits (D38'', opt-in)", {
  skip_if_not_installed("cmdstanr")
  if (!.gdpar_eb_smoke_env_active()) {
    skip("Set GDPAR_RUN_STAN_SMOKE_EB to a non-empty value to enable.")
  }
  set.seed(20260526L)
  n <- 80L
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  y_mat <- matrix(NA_real_, n, 2L)
  for (j in 1:2) {
    mu_j <- c(1.0, -0.5)[j] + c(0.5, -0.3)[j] * d$x1
    ls_j <- c(-0.2, 0.0)[j] + c(0.2, 0.1)[j] * d$x2
    y_mat[, j] <- rnorm(n, mean = mu_j, sd = exp(ls_j))
  }
  d$y <- y_mat
  fit_eb <- gdpar_eb(
    formula = y ~ x1 + x2,
    family  = gdpar_family("gaussian"),
    amm     = list(
      mu    = amm_spec(p = 2L, dims = dimwise(a = ~ x1)),
      sigma = amm_spec(p = 2L, dims = dimwise(a = ~ x2))
    ),
    data    = d,
    iter_warmup = 200L, iter_sampling = 200L, chains = 2L,
    refresh = 0L, seed = 35L, skip_id_check = TRUE,
    laplace_control = list(multi_start_M = 2L, laplace_draws = 200L),
    verbose = FALSE
  )
  expect_s3_class(fit_eb, "gdpar_eb_fit")
  expect_identical(fit_eb$path, "eb_KxP")
  expect_equal(fit_eb$K, 2L)
  expect_equal(fit_eb$p, 2L)
  expect_equal(dim(fit_eb$theta_ref_kp_hat), c(1L, 2L, 2L))
  expect_true(all(is.finite(fit_eb$theta_ref_kp_hat)))
  expect_true(isTRUE(fit_eb$correction_applied))
  expect_equal(dim(fit_eb$correction_tensor_constant), c(2L, 2L, 2L))
})

test_that("smoke 36/42: Path C Gaussian K=2 p=3 polynomial fits (D38'', opt-in)", {
  skip_if_not_installed("cmdstanr")
  if (!.gdpar_eb_smoke_env_active()) {
    skip("Set GDPAR_RUN_STAN_SMOKE_EB to a non-empty value to enable.")
  }
  set.seed(20260526L)
  n <- 100L
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  y_mat <- matrix(NA_real_, n, 3L)
  for (j in 1:3) {
    mu_j <- c(1.0, -0.5, 0.2)[j] + c(0.5, -0.3, 0.1)[j] * d$x1
    ls_j <- c(-0.2, 0.0, -0.1)[j] + c(0.2, 0.1, 0.0)[j] * d$x2
    y_mat[, j] <- rnorm(n, mean = mu_j, sd = exp(ls_j))
  }
  d$y <- y_mat
  fit_eb <- gdpar_eb(
    formula = y ~ x1 + x2,
    family  = gdpar_family("gaussian"),
    amm     = list(
      mu    = amm_spec(p = 3L, dims = dimwise(a = ~ x1)),
      sigma = amm_spec(p = 3L, dims = dimwise(a = ~ x2))
    ),
    data    = d,
    iter_warmup = 200L, iter_sampling = 200L, chains = 2L,
    refresh = 0L, seed = 36L, skip_id_check = TRUE,
    laplace_control = list(multi_start_M = 2L, laplace_draws = 200L),
    verbose = FALSE
  )
  expect_s3_class(fit_eb, "gdpar_eb_fit")
  expect_identical(fit_eb$path, "eb_KxP")
  expect_equal(fit_eb$K, 2L)
  expect_equal(fit_eb$p, 3L)
  expect_equal(dim(fit_eb$theta_ref_kp_hat), c(1L, 2L, 3L))
  expect_equal(dim(fit_eb$correction_tensor_constant), c(2L, 3L, 3L))
})

test_that("smoke 37/42: Path C NB K=2 p=2 polynomial fits (D40', opt-in)", {
  skip_if_not_installed("cmdstanr")
  if (!.gdpar_eb_smoke_env_active()) {
    skip("Set GDPAR_RUN_STAN_SMOKE_EB to a non-empty value to enable.")
  }
  set.seed(20260526L)
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
  fit_eb <- gdpar_eb(
    formula = y ~ x1 + x2,
    family  = gdpar_family("neg_binomial_2"),
    amm     = list(
      mu  = amm_spec(p = 2L, dims = dimwise(a = ~ x1)),
      phi = amm_spec(p = 2L, dims = dimwise(a = ~ x2))
    ),
    data    = d,
    iter_warmup = 200L, iter_sampling = 200L, chains = 2L,
    refresh = 0L, seed = 37L, skip_id_check = TRUE,
    laplace_control = list(multi_start_M = 2L, laplace_draws = 200L),
    verbose = FALSE
  )
  expect_s3_class(fit_eb, "gdpar_eb_fit")
  expect_identical(fit_eb$path, "eb_KxP")
  expect_equal(fit_eb$K, 2L)
  expect_equal(fit_eb$p, 2L)
})

test_that("smoke 38/42: Path C NB K=2 p=3 polynomial fits (D40', opt-in)", {
  skip_if_not_installed("cmdstanr")
  if (!.gdpar_eb_smoke_env_active()) {
    skip("Set GDPAR_RUN_STAN_SMOKE_EB to a non-empty value to enable.")
  }
  set.seed(20260526L)
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
  fit_eb <- gdpar_eb(
    formula = y ~ x1 + x2,
    family  = gdpar_family("neg_binomial_2"),
    amm     = list(
      mu  = amm_spec(p = 3L, dims = dimwise(a = ~ x1)),
      phi = amm_spec(p = 3L, dims = dimwise(a = ~ x2))
    ),
    data    = d,
    iter_warmup = 200L, iter_sampling = 200L, chains = 2L,
    refresh = 0L, seed = 38L, skip_id_check = TRUE,
    laplace_control = list(multi_start_M = 2L, laplace_draws = 200L),
    verbose = FALSE
  )
  expect_s3_class(fit_eb, "gdpar_eb_fit")
  expect_identical(fit_eb$path, "eb_KxP")
  expect_equal(fit_eb$K, 2L)
  expect_equal(fit_eb$p, 3L)
})

test_that("smoke 39/42: Path C Beta K=2 p=2 rejected by D40' guard", {
  set.seed(20260526L)
  n <- 50L
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n),
                  y = matrix(runif(n * 2L, 0.1, 0.9), n, 2L))
  expect_error(
    gdpar_eb(
      formula = y ~ x1 + x2,
      family  = gdpar_family("beta"),
      amm     = list(
        mu  = amm_spec(p = 2L, dims = dimwise(a = ~ x1)),
        phi = amm_spec(p = 2L, dims = dimwise(a = ~ x2))
      ),
      data    = d,
      iter_warmup = 50L, iter_sampling = 50L, chains = 1L,
      refresh = 0L, seed = 39L, skip_id_check = TRUE,
      verbose = FALSE
    ),
    class = "gdpar_unsupported_feature_error"
  )
})

test_that("smoke 40/42: Path C Gamma K=2 p=2 rejected by D40' guard", {
  set.seed(20260526L)
  n <- 50L
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n),
                  y = matrix(rgamma(n * 2L, 2, 1), n, 2L))
  expect_error(
    gdpar_eb(
      formula = y ~ x1 + x2,
      family  = gdpar_family("gamma"),
      amm     = list(
        mu    = amm_spec(p = 2L, dims = dimwise(a = ~ x1)),
        shape = amm_spec(p = 2L, dims = dimwise(a = ~ x2))
      ),
      data    = d,
      iter_warmup = 50L, iter_sampling = 50L, chains = 1L,
      refresh = 0L, seed = 40L, skip_id_check = TRUE,
      verbose = FALSE
    ),
    class = "gdpar_unsupported_feature_error"
  )
})

test_that("smoke 41/42: Path C Student-t K=3 p=2 rejected by D40' guard", {
  set.seed(20260526L)
  n <- 50L
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n),
                  y = matrix(rt(n * 2L, df = 5), n, 2L))
  expect_error(
    gdpar_eb(
      formula = y ~ x1 + x2,
      family  = gdpar_family("student_t"),
      amm     = list(
        mu    = amm_spec(p = 2L, dims = dimwise(a = ~ x1)),
        sigma = amm_spec(p = 2L, dims = dimwise(a = ~ x2)),
        nu    = amm_spec(p = 2L, dims = dimwise(a = ~ x1))
      ),
      data    = d,
      iter_warmup = 50L, iter_sampling = 50L, chains = 1L,
      refresh = 0L, seed = 41L, skip_id_check = TRUE,
      verbose = FALSE
    ),
    class = "gdpar_unsupported_feature_error"
  )
})

test_that("smoke 42/42: Path C ZIP K=2 p=2 rejected by D40' guard", {
  set.seed(20260526L)
  n <- 50L
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n),
                  y = matrix(rpois(n * 2L, lambda = 1.5), n, 2L))
  expect_error(
    gdpar_eb(
      formula = y ~ x1 + x2,
      family  = gdpar_family("zip"),
      amm     = list(
        mu = amm_spec(p = 2L, dims = dimwise(a = ~ x1)),
        pi = amm_spec(p = 2L, dims = dimwise(a = ~ x2))
      ),
      data    = d,
      iter_warmup = 50L, iter_sampling = 50L, chains = 1L,
      refresh = 0L, seed = 42L, skip_id_check = TRUE,
      verbose = FALSE
    ),
    class = "gdpar_unsupported_feature_error"
  )
})

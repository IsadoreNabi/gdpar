# Sub-phase 8.3.9 (2026-05-22): structural tests for G1/G2/G3 residuals,
# gdpar_posterior_predict, pp_check.gdpar_fit, gdpar_dharma_object.
#
# Structural tests (no Stan): mock gdpar_fit objects with synthetic
# draws of y_pred suffice to exercise the path-dispatching, family
# classification, and per-type residual computations.
#
# Gated smokes against real fits are factored into
# test-fit_K_smoke.R and test-fit_multi_smoke.R via env vars.

# ---- Helpers ------------------------------------------------------------

.mock_scalar_fit <- function(family_name = "gaussian",
                             y_obs = NULL,
                             y_pred_mat = NULL,
                             S = 100L,
                             n = 20L,
                             seed = 42L) {
  set.seed(seed)
  if (is.null(y_obs)) {
    y_obs <- rnorm(n)
  } else {
    n <- length(y_obs)
  }
  if (is.null(y_pred_mat)) {
    y_pred_mat <- matrix(rnorm(S * n, mean = mean(y_obs), sd = sd(y_obs)),
                         nrow = S)
  } else {
    S <- nrow(y_pred_mat)
    n <- ncol(y_pred_mat)
  }
  fam <- gdpar_family(family_name)
  needs_int <- fam$stan_id %in% .gdpar_integer_family_stan_ids()
  stan_data <- list(
    y_real = if (needs_int) numeric(n) else as.numeric(y_obs),
    y_int  = if (needs_int) as.integer(y_obs) else integer(n)
  )
  # Build a minimal fake cmdstanr-like fit that has $draws(variables, format)
  # returning a draws_matrix compatible with posterior::as_draws_matrix.
  draws_arr <- posterior::as_draws_matrix(
    matrix(as.numeric(y_pred_mat),
           nrow = S, ncol = n,
           dimnames = list(NULL, paste0("y_pred[", seq_len(n), "]")))
  )
  fake_fit <- list()
  fake_fit$draws <- function(variables, format = "draws_matrix") {
    draws_arr
  }
  obj <- list(
    fit       = fake_fit,
    family    = fam,
    stan_data = stan_data,
    path      = "bayes"
  )
  class(obj) <- c("gdpar_fit", "list")
  obj
}

.mock_K_fit <- function(K = 2L, family_name = "gaussian",
                        y_obs = NULL, S = 80L, n = 15L, seed = 11L) {
  fit_obj <- .mock_scalar_fit(family_name = family_name, y_obs = y_obs,
                              S = S, n = n, seed = seed)
  fit_obj$K <- K
  fit_obj$slot_names <- if (K == 2L) c("mu", "sigma") else paste0("slot_", seq_len(K))
  fit_obj
}

.mock_multivariate_fit <- function(p = 2L, S = 60L, n = 12L, seed = 7L) {
  set.seed(seed)
  y_obs <- matrix(rnorm(n * p), nrow = n, ncol = p)
  y_pred_arr <- array(rnorm(S * n * p), dim = c(S, n, p))
  var_names <- character(0)
  flat_cols <- numeric(0)
  for (k in seq_len(p)) {
    for (i in seq_len(n)) {
      var_names <- c(var_names, sprintf("y_pred[%d,%d]", i, k))
      flat_cols <- c(flat_cols, list(y_pred_arr[, i, k]))
    }
  }
  mat <- do.call(cbind, flat_cols)
  colnames(mat) <- var_names
  draws_arr <- posterior::as_draws_matrix(mat)
  fake_fit <- list()
  fake_fit$draws <- function(variables, format = "draws_matrix") draws_arr
  fams <- replicate(p, gdpar_family("gaussian"), simplify = FALSE)
  multi_fam <- list(
    name = "gaussian_multi",
    families = fams,
    p = p
  )
  class(multi_fam) <- c("gdpar_family_multi", "list")
  stan_data <- list(y_real = y_obs, y_int = matrix(0L, nrow = n, ncol = p))
  obj <- list(
    fit       = fake_fit,
    family    = multi_fam,
    stan_data = stan_data,
    p         = p,
    path      = "bayes"
  )
  class(obj) <- c("gdpar_fit", "list")
  obj
}

# ---- Section 1: path classification + integer-family table -------------

test_that("`.gdpar_fit_path_class` returns 'scalar' for K=1+p=1 fits", {
  obj <- .mock_scalar_fit()
  expect_identical(.gdpar_fit_path_class(obj), "scalar")
})

test_that("`.gdpar_fit_path_class` returns 'K_individual' for K>1 fits", {
  obj <- .mock_K_fit(K = 2L)
  expect_identical(.gdpar_fit_path_class(obj), "K_individual")
  obj3 <- .mock_K_fit(K = 3L, family_name = "student_t")
  expect_identical(.gdpar_fit_path_class(obj3), "K_individual")
})

test_that("`.gdpar_fit_path_class` returns 'multivariate' for p>1 fits", {
  obj <- .mock_multivariate_fit(p = 2L)
  expect_identical(.gdpar_fit_path_class(obj), "multivariate")
})

test_that(".gdpar_integer_family_stan_ids covers Poisson/NB/Bernoulli/mixtures", {
  ids <- .gdpar_integer_family_stan_ids()
  expect_true(2L %in% ids)
  expect_true(3L %in% ids)
  expect_true(4L %in% ids)
  expect_true(10L %in% ids)
  expect_true(11L %in% ids)
  expect_true(12L %in% ids)
  expect_true(13L %in% ids)
  expect_false(1L %in% ids)
  expect_false(5L %in% ids)
  expect_false(6L %in% ids)
  expect_false(8L %in% ids)
  expect_false(9L %in% ids)
})

# ---- Section 2: family discreteness classifier --------------------------

test_that(".gdpar_family_is_discrete classifies all families correctly", {
  expect_false(.gdpar_family_is_discrete("gaussian"))
  expect_false(.gdpar_family_is_discrete("beta"))
  expect_false(.gdpar_family_is_discrete("gamma"))
  expect_false(.gdpar_family_is_discrete("student_t"))
  expect_false(.gdpar_family_is_discrete("tweedie"))
  expect_false(.gdpar_family_is_discrete("lognormal_loc_scale"))
  expect_true(.gdpar_family_is_discrete("poisson"))
  expect_true(.gdpar_family_is_discrete("neg_binomial_2"))
  expect_true(.gdpar_family_is_discrete("bernoulli"))
  expect_true(.gdpar_family_is_discrete("zip"))
  expect_true(.gdpar_family_is_discrete("zinb"))
  expect_true(.gdpar_family_is_discrete("hurdle_poisson"))
  expect_true(.gdpar_family_is_discrete("hurdle_neg_binomial_2"))
})

# ---- Section 3: deviance residual analytic checks -----------------------

test_that(".gdpar_deviance_residual returns sign(y - mu) * |y - mu| / sigma for gaussian", {
  y  <- c(1, 2, 3, 5)
  mu <- c(0, 2, 5, 4)
  d  <- .gdpar_deviance_residual(y, mu, "gaussian")
  expect_equal(d, c(1, 0, -2, 1))
})

test_that(".gdpar_deviance_residual returns 0 for y = mu in all families", {
  for (fam in c("gaussian", "poisson", "neg_binomial_2", "bernoulli",
                "beta", "gamma", "lognormal_loc_scale", "student_t",
                "tweedie", "zip", "zinb", "hurdle_poisson",
                "hurdle_neg_binomial_2")) {
    y  <- c(1, 2, 3)
    mu <- c(1, 2, 3)
    d  <- .gdpar_deviance_residual(y, mu, fam)
    expect_equal(d, c(0, 0, 0), tolerance = 1e-8,
                 info = sprintf("family=%s", fam))
  }
})

test_that(".gdpar_deviance_residual sign convention is y > mu => positive", {
  for (fam in c("gaussian", "poisson", "neg_binomial_2",
                "beta", "gamma", "lognormal_loc_scale")) {
    y  <- c(2, 1)
    mu <- c(1, 2)
    d  <- .gdpar_deviance_residual(y, mu, fam)
    expect_true(d[1L] >= 0, info = sprintf("family=%s y>mu", fam))
    expect_true(d[2L] <= 0, info = sprintf("family=%s y<mu", fam))
  }
})

test_that(".gdpar_deviance_residual returns numeric vector of length n", {
  for (fam in c("gaussian", "poisson", "beta", "gamma", "student_t",
                "tweedie", "zip", "zinb", "hurdle_poisson")) {
    n <- 20L
    y <- if (fam %in% c("poisson", "zip", "hurdle_poisson", "neg_binomial_2",
                        "zinb", "hurdle_neg_binomial_2")) {
      sample(0:10, n, replace = TRUE)
    } else if (fam == "beta") {
      runif(n, 0.1, 0.9)
    } else if (fam == "gamma") {
      rgamma(n, shape = 2, rate = 1)
    } else if (fam == "bernoulli") {
      sample(0:1, n, replace = TRUE)
    } else {
      rnorm(n)
    }
    mu <- if (fam == "beta") {
      runif(n, 0.1, 0.9)
    } else if (fam == "gamma") {
      rgamma(n, shape = 2, rate = 1)
    } else {
      abs(rnorm(n)) + 0.5
    }
    d <- .gdpar_deviance_residual(y, mu, fam)
    expect_length(d, n)
    expect_true(is.numeric(d))
    expect_true(all(is.finite(d)),
                info = sprintf("family=%s has non-finite residuals", fam))
  }
})

# ---- Section 4: quantile residual Bayesian -------------------------------

test_that(".gdpar_quantile_residuals_bayesian: well-specified Gaussian gives ~N(0,1)", {
  set.seed(1234L)
  n <- 200L
  S <- 800L
  mu_true <- rnorm(n)
  sigma_true <- 1
  y_obs <- rnorm(n, mean = mu_true, sd = sigma_true)
  y_pred_mat <- matrix(NA_real_, nrow = S, ncol = n)
  for (i in seq_len(n)) {
    y_pred_mat[, i] <- rnorm(S, mean = mu_true[i], sd = sigma_true)
  }
  qr <- .gdpar_quantile_residuals_bayesian(y_obs, y_pred_mat, "gaussian")
  expect_length(qr, n)
  expect_true(all(is.finite(qr)))
  ks <- suppressWarnings(stats::ks.test(qr, "pnorm"))
  expect_gt(ks$p.value, 0.001)
})

test_that(".gdpar_quantile_residuals_bayesian: discrete family randomization runs", {
  set.seed(2345L)
  n <- 30L
  S <- 200L
  lambda <- 3
  y_obs <- rpois(n, lambda)
  y_pred_mat <- matrix(rpois(S * n, lambda), nrow = S)
  qr <- .gdpar_quantile_residuals_bayesian(
    y_obs, y_pred_mat, "poisson", randomize_seed = 7L
  )
  expect_length(qr, n)
  expect_true(all(is.finite(qr)))
})

test_that(".gdpar_quantile_residuals_bayesian: randomize_seed gives reproducible output", {
  set.seed(3456L)
  n <- 20L; S <- 100L
  y_obs <- rpois(n, 2)
  y_pred_mat <- matrix(rpois(S * n, 2), nrow = S)
  qr1 <- .gdpar_quantile_residuals_bayesian(y_obs, y_pred_mat, "poisson",
                                            randomize_seed = 1L)
  qr2 <- .gdpar_quantile_residuals_bayesian(y_obs, y_pred_mat, "poisson",
                                            randomize_seed = 1L)
  expect_identical(qr1, qr2)
})

# ---- Section 5: residuals.gdpar_fit on mock fits ------------------------

test_that("residuals.gdpar_fit returns vector of length n in scalar path (all types)", {
  obj <- .mock_scalar_fit(n = 30L, S = 200L)
  for (t in c("response", "pearson", "deviance", "quantile")) {
    r <- residuals(obj, type = t)
    expect_length(r, 30L)
    expect_true(is.numeric(r))
    expect_true(all(is.finite(r)))
  }
})

test_that("residuals.gdpar_fit returns vector of length n in K_individual path", {
  obj <- .mock_K_fit(K = 2L, family_name = "gaussian", n = 20L, S = 150L)
  for (t in c("response", "pearson", "deviance", "quantile")) {
    r <- residuals(obj, type = t)
    expect_length(r, 20L)
    expect_true(is.numeric(r))
  }
})

test_that("residuals.gdpar_fit returns matrix n*p in multivariate path", {
  obj <- .mock_multivariate_fit(p = 2L, n = 15L, S = 80L)
  for (t in c("response", "pearson", "deviance", "quantile")) {
    r <- residuals(obj, type = t)
    expect_true(is.matrix(r))
    expect_equal(dim(r), c(15L, 2L))
    expect_identical(colnames(r), c("dim_1", "dim_2"))
  }
})

test_that("residuals.gdpar_fit with coord=k returns vector in multivariate path", {
  obj <- .mock_multivariate_fit(p = 3L, n = 12L, S = 80L)
  r1 <- residuals(obj, type = "quantile", coord = 1L)
  r2 <- residuals(obj, type = "quantile", coord = 2L)
  expect_length(r1, 12L)
  expect_length(r2, 12L)
  expect_false(identical(r1, r2))
})

test_that("residuals.gdpar_fit rejects coord out of range in multivariate path", {
  obj <- .mock_multivariate_fit(p = 2L, n = 10L, S = 50L)
  expect_error(
    residuals(obj, type = "quantile", coord = 3L),
    class = "gdpar_input_error"
  )
  expect_error(
    residuals(obj, type = "quantile", coord = 0L),
    class = "gdpar_input_error"
  )
})

test_that("residuals.gdpar_fit rejects unknown type", {
  obj <- .mock_scalar_fit()
  expect_error(residuals(obj, type = "foo"))
})

test_that("residuals.gdpar_fit randomize_seed yields reproducible output", {
  obj <- .mock_scalar_fit(family_name = "poisson", n = 15L, S = 80L,
                          y_obs = rpois(15L, 3))
  r1 <- residuals(obj, type = "quantile", randomize_seed = 11L)
  r2 <- residuals(obj, type = "quantile", randomize_seed = 11L)
  expect_identical(r1, r2)
})

# ---- Section 6: gdpar_posterior_predict on mock fits --------------------

test_that("gdpar_posterior_predict returns S*n matrix in scalar path", {
  obj <- .mock_scalar_fit(n = 25L, S = 120L)
  pp <- gdpar_posterior_predict(obj)
  expect_true(is.matrix(pp))
  expect_equal(dim(pp), c(120L, 25L))
})

test_that("gdpar_posterior_predict returns S*n*p array in multivariate path", {
  obj <- .mock_multivariate_fit(p = 3L, n = 10L, S = 50L)
  pp <- gdpar_posterior_predict(obj)
  expect_true(is.array(pp))
  expect_equal(dim(pp), c(50L, 10L, 3L))
})

test_that("gdpar_posterior_predict subsamples by ndraws", {
  obj <- .mock_scalar_fit(n = 10L, S = 100L)
  pp_full <- gdpar_posterior_predict(obj)
  pp_sub  <- gdpar_posterior_predict(obj, ndraws = 30L)
  expect_equal(nrow(pp_sub), 30L)
  expect_equal(pp_sub[1:5, ], pp_full[1:5, ])
})

test_that("gdpar_posterior_predict rejects invalid ndraws", {
  obj <- .mock_scalar_fit(n = 10L, S = 50L)
  expect_error(gdpar_posterior_predict(obj, ndraws = -1L),
               class = "gdpar_input_error")
  expect_error(gdpar_posterior_predict(obj, ndraws = NA_integer_),
               class = "gdpar_input_error")
})

test_that("gdpar_posterior_predict rejects non-gdpar_fit input", {
  expect_error(gdpar_posterior_predict(list()),
               class = "gdpar_input_error")
})

# ---- Section 7: gdpar_dharma_object detection ----------------------------

test_that("gdpar_dharma_object errors out when DHARMa is not installed", {
  if (.gdpar_has_dharma()) {
    skip("DHARMa is installed; cannot test the not-installed branch")
  }
  obj <- .mock_scalar_fit(n = 20L, S = 100L)
  expect_error(gdpar_dharma_object(obj),
               class = "gdpar_input_error")
})

test_that("gdpar_dharma_object works when DHARMa is installed (gated)", {
  skip_if_not_installed("DHARMa")
  obj <- .mock_scalar_fit(n = 30L, S = 150L)
  d <- gdpar_dharma_object(obj)
  expect_true(inherits(d, "DHARMa"))
})

test_that("gdpar_dharma_object multivariate requires coord", {
  skip_if_not_installed("DHARMa")
  obj <- .mock_multivariate_fit(p = 2L, n = 15L, S = 80L)
  expect_error(gdpar_dharma_object(obj),
               class = "gdpar_input_error")
  d <- gdpar_dharma_object(obj, coord = 1L)
  expect_true(inherits(d, "DHARMa"))
})

# ---- Section 8: pp_check.gdpar_fit dispatch ------------------------------

test_that("pp_check.gdpar_fit returns a ggplot when bayesplot is installed", {
  skip_if_not_installed("bayesplot")
  skip_if_not_installed("ggplot2")
  obj <- .mock_scalar_fit(n = 30L, S = 100L)
  g <- suppressWarnings(pp_check.gdpar_fit(obj, type = "dens_overlay",
                                           ndraws = 20L))
  expect_true(inherits(g, "ggplot"))
})

test_that("pp_check.gdpar_fit multivariate requires coord", {
  skip_if_not_installed("bayesplot")
  obj <- .mock_multivariate_fit(p = 2L, n = 15L, S = 60L)
  expect_error(pp_check.gdpar_fit(obj, type = "dens_overlay"),
               class = "gdpar_input_error")
})

test_that("pp_check.gdpar_fit rejects unknown type", {
  obj <- .mock_scalar_fit()
  expect_error(pp_check.gdpar_fit(obj, type = "foo"))
})

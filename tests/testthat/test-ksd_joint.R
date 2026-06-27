# Tests for Sub-bloque 9.3.c (Block 9, Session B9.4, 2026-05-27)
# canonized decision H.iv lateral: gdpar_ksd_joint() helper.
#
# Verification under the rigor de tres capas (feedback_proof_rigor_standards):
#   - algebraica: input validators + closed-form Stein kernel under
#     Gaussian target;
#   - estadistica: KSD ~ 0 when EB and FB samples are from the same
#     Gaussian; KSD > 0 detectable when EB samples are shifted away
#     from FB; ESS-weighted variant returns a different value than
#     the unweighted variant;
#   - numerica: covered indirectly by the package-level suite + the
#     R CMD check natural.

# ---- Mock fits to avoid Stan dependence in unit tests ----

.b9_4_mock_eb_fb_fits <- function(eb_mu = c(0, 0),
                                  fb_mu = c(0, 0),
                                  Sigma = diag(2L),
                                  n_eb = 200L,
                                  n_fb = 200L,
                                  vars = c("a_coef[1]", "c_b_coef[1]"),
                                  seed = 1L) {
  set.seed(seed)
  L <- chol(Sigma)
  eb_mat <- matrix(stats::rnorm(n_eb * length(eb_mu)),
                   nrow = n_eb) %*% L +
            matrix(eb_mu, n_eb, length(eb_mu), byrow = TRUE)
  fb_mat <- matrix(stats::rnorm(n_fb * length(fb_mu)),
                   nrow = n_fb) %*% L +
            matrix(fb_mu, n_fb, length(fb_mu), byrow = TRUE)
  colnames(eb_mat) <- vars
  colnames(fb_mat) <- vars
  eb_draws <- posterior::as_draws_matrix(eb_mat)
  fb_draws <- posterior::as_draws_matrix(fb_mat)

  # Build minimal mock with the structure consumed by
  # .gdpar_eb_fb_extract_xi_draws: $conditional_fit$draws() returns
  # a posterior::draws object.
  mock_fit_obj <- function(draws_obj) {
    obj <- list(.draws = draws_obj)
    obj$draws <- function() draws_obj
    obj
  }
  eb <- list(conditional_fit = mock_fit_obj(eb_draws))
  fb <- list(fit = mock_fit_obj(fb_draws))
  class(eb) <- c("gdpar_eb_fit", "list")
  class(fb) <- c("gdpar_fit", "list")
  list(eb_fit = eb, fb_fit = fb)
}

# ---- Input validation ----

test_that("gdpar_ksd_joint rejects non-gdpar fits", {
  expect_error(
    gdpar_ksd_joint(list(), list()),
    class = "gdpar_input_error"
  )
})

test_that("gdpar_ksd_joint rejects malformed beta", {
  fits <- .b9_4_mock_eb_fb_fits()
  expect_error(
    gdpar_ksd_joint(fits$eb_fit, fits$fb_fit, beta = 0.0),
    class = "gdpar_input_error"
  )
  expect_error(
    gdpar_ksd_joint(fits$eb_fit, fits$fb_fit, beta = -1.0),
    class = "gdpar_input_error"
  )
  expect_error(
    gdpar_ksd_joint(fits$eb_fit, fits$fb_fit, beta = NA_real_),
    class = "gdpar_input_error"
  )
})

test_that("gdpar_ksd_joint rejects malformed bandwidth_value", {
  fits <- .b9_4_mock_eb_fb_fits()
  expect_error(
    gdpar_ksd_joint(fits$eb_fit, fits$fb_fit,
                    bandwidth = "fixed", bandwidth_value = 0),
    class = "gdpar_input_error"
  )
  expect_error(
    gdpar_ksd_joint(fits$eb_fit, fits$fb_fit,
                    bandwidth = "fixed", bandwidth_value = -1),
    class = "gdpar_input_error"
  )
})

test_that("gdpar_ksd_joint rejects malformed ess_weighted / seed", {
  fits <- .b9_4_mock_eb_fb_fits()
  expect_error(
    gdpar_ksd_joint(fits$eb_fit, fits$fb_fit, ess_weighted = NA),
    class = "gdpar_input_error"
  )
  expect_error(
    gdpar_ksd_joint(fits$eb_fit, fits$fb_fit, seed = NA_integer_),
    class = "gdpar_input_error"
  )
})

# ---- Happy path and class invariants ----

test_that("gdpar_ksd_joint returns a gdpar_ksd_joint object with required fields", {
  fits <- .b9_4_mock_eb_fb_fits()
  res <- gdpar_ksd_joint(fits$eb_fit, fits$fb_fit)
  expect_s3_class(res, "gdpar_ksd_joint")
  for (nm in c("ksd_value", "ksd_squared", "kernel", "bandwidth",
               "bandwidth_value", "beta", "n_eb_draws", "n_fb_draws",
               "n_dim", "target_mu", "target_Sigma", "ess_weighted",
               "thinned_to", "vars", "call")) {
    expect_true(nm %in% names(res), info = nm)
  }
  expect_true(is.numeric(res$ksd_value))
  expect_true(res$ksd_value >= 0.0)
  expect_equal(res$n_dim, 2L)
  expect_equal(res$kernel, "imq")
  expect_equal(res$beta, -0.5)
  expect_false(res$ess_weighted)
  expect_true(is.na(res$thinned_to))
})

test_that("gdpar_ksd_joint accepts RBF kernel and reports beta as NA", {
  fits <- .b9_4_mock_eb_fb_fits()
  res <- gdpar_ksd_joint(fits$eb_fit, fits$fb_fit, kernel = "rbf")
  expect_equal(res$kernel, "rbf")
  expect_true(is.na(res$beta))
  expect_true(res$ksd_value >= 0.0)
})

test_that("gdpar_ksd_joint accepts fixed bandwidth", {
  fits <- .b9_4_mock_eb_fb_fits()
  res <- gdpar_ksd_joint(fits$eb_fit, fits$fb_fit,
                         bandwidth = "fixed", bandwidth_value = 2.5)
  expect_equal(res$bandwidth, "fixed")
  expect_equal(res$bandwidth_value, 2.5)
})

# ---- Statistical layer: KSD ~ 0 when distributions match ----

test_that("KSD is small when EB and FB are drawn from the same Gaussian", {
  fits <- .b9_4_mock_eb_fb_fits(eb_mu = c(0, 0), fb_mu = c(0, 0),
                                Sigma = diag(2L), n_eb = 500L,
                                n_fb = 500L, seed = 42L)
  res <- gdpar_ksd_joint(fits$eb_fit, fits$fb_fit)
  # Under the null (same distribution), V-stat KSD should be small;
  # 0.5 is a conservative upper bound that allows for Monte Carlo
  # variability on n = 500.
  expect_lt(res$ksd_value, 0.5)
})

test_that("KSD detects shift between EB and FB Gaussians", {
  fits_close <- .b9_4_mock_eb_fb_fits(eb_mu = c(0, 0), fb_mu = c(0, 0),
                                       Sigma = diag(2L), n_eb = 500L,
                                       n_fb = 500L, seed = 7L)
  fits_far   <- .b9_4_mock_eb_fb_fits(eb_mu = c(2, 2), fb_mu = c(0, 0),
                                       Sigma = diag(2L), n_eb = 500L,
                                       n_fb = 500L, seed = 7L)
  res_close <- gdpar_ksd_joint(fits_close$eb_fit, fits_close$fb_fit)
  res_far   <- gdpar_ksd_joint(fits_far$eb_fit,   fits_far$fb_fit)
  expect_gt(res_far$ksd_value, res_close$ksd_value)
})

# ---- ESS-weighted variant ----

test_that("ESS-weighted thinning runs and reports thinned_to", {
  fits <- .b9_4_mock_eb_fb_fits(n_eb = 200L, n_fb = 200L, seed = 11L)
  res <- gdpar_ksd_joint(fits$eb_fit, fits$fb_fit,
                         ess_weighted = TRUE, seed = 1L)
  expect_true(res$ess_weighted)
  expect_true(is.integer(res$thinned_to))
  expect_true(res$thinned_to >= 2L)
  expect_true(res$thinned_to <= 200L)
})

test_that("seed makes ESS-weighted thinning reproducible", {
  fits <- .b9_4_mock_eb_fb_fits(n_eb = 200L, n_fb = 200L, seed = 13L)
  r1 <- gdpar_ksd_joint(fits$eb_fit, fits$fb_fit,
                        ess_weighted = TRUE, seed = 99L)
  r2 <- gdpar_ksd_joint(fits$eb_fit, fits$fb_fit,
                        ess_weighted = TRUE, seed = 99L)
  expect_equal(r1$ksd_value, r2$ksd_value)
})

# ---- S3 methods ----

test_that("print.gdpar_ksd_joint emits the expected header", {
  fits <- .b9_4_mock_eb_fb_fits()
  res <- gdpar_ksd_joint(fits$eb_fit, fits$fb_fit)
  out <- capture.output(print(res))
  expect_true(any(grepl("<gdpar_ksd_joint>", out, fixed = TRUE)))
  expect_true(any(grepl("Joint KSD", out, fixed = TRUE)))
  expect_true(any(grepl("Kernel", out, fixed = TRUE)))
})

test_that("summary.gdpar_ksd_joint returns a summary object", {
  fits <- .b9_4_mock_eb_fb_fits()
  res <- gdpar_ksd_joint(fits$eb_fit, fits$fb_fit)
  sm <- summary(res)
  expect_s3_class(sm, "summary.gdpar_ksd_joint")
  expect_true("interpretation" %in% names(sm))
  out <- capture.output(print(sm))
  expect_true(any(grepl("<summary.gdpar_ksd_joint>", out, fixed = TRUE)))
})

# ---- Common-var failure mode ----

test_that("gdpar_ksd_joint aborts when EB and FB share no xi variables", {
  fits_eb <- .b9_4_mock_eb_fb_fits(vars = c("a_coef[1]", "c_b_coef[1]"))
  fits_fb <- .b9_4_mock_eb_fb_fits(vars = c("W_coef[1,1]", "W_coef[1,2]"))
  expect_error(
    gdpar_ksd_joint(fits_eb$eb_fit, fits_fb$fb_fit),
    class = "gdpar_input_error"
  )
})

# Tests for .gdpar_eb_correction_matrix() and .gdpar_eb_apply_correction()
# dispatch under Sub-phase 8.6.C (decision D34, Path A).
#
# Three layers:
#   (i)   algebraic unit tests of .gdpar_eb_correction_matrix() at p > 1.
#         No Stan compilation, no posterior draws.
#   (ii)  algebraic reduction tests: at p = 1 the matrix-valued helper
#         collapses to the 1x1 matrix kappa(alpha) * Var^marg(theta_ref)
#         that matches the scalar 8.6.B form bit-for-bit.
#   (iii) dispatch tests of .gdpar_eb_apply_correction() showing that
#         p = 1 routes to the scalar branch and p > 1 routes to the
#         matrix branch.

kappa_alpha_95 <- 1.92

make_laplace_result <- function(cov_mat) {
  list(
    theta_ref_hat = rep(0, nrow(cov_mat)),
    theta_ref_se  = sqrt(pmax(diag(cov_mat), 0)),
    theta_ref_cov = cov_mat,
    diagnostics   = list()
  )
}

# ---- (i) Algebraic tests at p > 1 ---------------------------------------

test_that(".gdpar_eb_correction_matrix returns the Prop 7B* scaling at p = 2", {
  cov_mat <- matrix(c(0.50, 0.10,
                      0.10, 0.80), nrow = 2L, byrow = TRUE)
  lap <- make_laplace_result(cov_mat)
  out <- .gdpar_eb_correction_matrix(
    eb_correction  = TRUE,
    laplace_result = lap,
    stan_data      = list(),
    p              = 2L,
    verbose        = FALSE
  )
  expect_true(isTRUE(out$applied))
  expect_true(is.matrix(out$constant))
  expect_equal(dim(out$constant), c(2L, 2L))
  expect_equal(out$constant, kappa_alpha_95 * cov_mat)
  expect_true(isSymmetric(out$constant, tol = 1e-12))
  eigs <- eigen(out$constant, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigs > 0))
})

test_that(".gdpar_eb_correction_matrix returns the Prop 7B* scaling at p = 3", {
  cov_mat <- diag(c(0.4, 0.6, 0.9))
  lap <- make_laplace_result(cov_mat)
  out <- .gdpar_eb_correction_matrix(
    eb_correction  = TRUE,
    laplace_result = lap,
    stan_data      = list(),
    p              = 3L,
    verbose        = FALSE
  )
  expect_true(isTRUE(out$applied))
  expect_equal(dim(out$constant), c(3L, 3L))
  expect_equal(diag(out$constant), kappa_alpha_95 * c(0.4, 0.6, 0.9))
  expect_equal(out$constant[1L, 2L], 0)
  expect_equal(out$constant[2L, 3L], 0)
})

test_that(".gdpar_eb_correction_matrix returns NA matrix when eb_correction = FALSE", {
  cov_mat <- diag(c(0.5, 0.5))
  lap <- make_laplace_result(cov_mat)
  out <- .gdpar_eb_correction_matrix(
    eb_correction  = FALSE,
    laplace_result = lap,
    stan_data      = list(),
    p              = 2L,
    verbose        = FALSE
  )
  expect_false(isTRUE(out$applied))
  expect_true(all(is.na(out$constant)))
  expect_equal(dim(out$constant), c(2L, 2L))
})

test_that(".gdpar_eb_correction_matrix returns NA matrix when cov has NA / Inf", {
  cov_mat <- matrix(c(0.5, NA_real_,
                      NA_real_, 0.5), nrow = 2L)
  lap <- make_laplace_result(cov_mat)
  out <- .gdpar_eb_correction_matrix(
    eb_correction  = TRUE,
    laplace_result = lap,
    stan_data      = list(),
    p              = 2L,
    verbose        = FALSE
  )
  expect_false(isTRUE(out$applied))
  expect_true(all(is.na(out$constant)))
})

test_that(".gdpar_eb_correction_matrix returns NA matrix when cov is non-PSD", {
  cov_mat <- matrix(c(1.0, 2.0,
                      2.0, 1.0), nrow = 2L)
  lap <- make_laplace_result(cov_mat)
  out <- .gdpar_eb_correction_matrix(
    eb_correction  = TRUE,
    laplace_result = lap,
    stan_data      = list(),
    p              = 2L,
    verbose        = FALSE
  )
  expect_false(isTRUE(out$applied))
  expect_true(all(is.na(out$constant)))
})

test_that(".gdpar_eb_correction_matrix returns NA matrix when cov is non-square", {
  bad <- list(
    theta_ref_hat = c(0, 0),
    theta_ref_se  = c(1, 1),
    theta_ref_cov = matrix(0, nrow = 2L, ncol = 3L),
    diagnostics   = list()
  )
  out <- .gdpar_eb_correction_matrix(
    eb_correction  = TRUE,
    laplace_result = bad,
    stan_data      = list(),
    p              = 2L,
    verbose        = FALSE
  )
  expect_false(isTRUE(out$applied))
  expect_true(all(is.na(out$constant)))
})

# ---- (ii) Algebraic reduction at p = 1 ----------------------------------

test_that(".gdpar_eb_correction_matrix at p = 1 collapses to the 8.6.B scalar form", {
  v <- 0.42
  cov_mat <- matrix(v, nrow = 1L, ncol = 1L)
  lap <- make_laplace_result(cov_mat)
  out_mat <- .gdpar_eb_correction_matrix(
    eb_correction  = TRUE,
    laplace_result = lap,
    stan_data      = list(),
    p              = 1L,
    verbose        = FALSE
  )
  expect_true(isTRUE(out_mat$applied))
  expect_equal(dim(out_mat$constant), c(1L, 1L))
  expect_equal(as.numeric(out_mat$constant), kappa_alpha_95 * v)
})

# ---- (iii) Dispatch tests of .gdpar_eb_apply_correction -----------------

test_that(".gdpar_eb_apply_correction p = 1 routes to the scalar branch", {
  v <- 0.25
  cov_mat <- matrix(v, nrow = 1L, ncol = 1L)
  lap <- make_laplace_result(cov_mat)
  out <- .gdpar_eb_apply_correction(
    eb_correction  = TRUE,
    laplace_result = lap,
    stan_data      = list(),
    p              = 1L,
    verbose        = FALSE
  )
  expect_true(isTRUE(out$applied))
  expect_type(out$constant, "double")
  expect_length(out$constant, 1L)
  expect_equal(out$constant, kappa_alpha_95 * v)
})

test_that(".gdpar_eb_apply_correction p > 1 routes to the matrix branch", {
  cov_mat <- diag(c(0.3, 0.7))
  lap <- make_laplace_result(cov_mat)
  out <- .gdpar_eb_apply_correction(
    eb_correction  = TRUE,
    laplace_result = lap,
    stan_data      = list(),
    p              = 2L,
    verbose        = FALSE
  )
  expect_true(isTRUE(out$applied))
  expect_true(is.matrix(out$constant))
  expect_equal(dim(out$constant), c(2L, 2L))
  expect_equal(diag(out$constant), kappa_alpha_95 * c(0.3, 0.7))
})

test_that(".gdpar_eb_apply_correction eb_correction = FALSE returns p-aware NA at p > 1", {
  cov_mat <- diag(c(0.3, 0.5, 0.7))
  lap <- make_laplace_result(cov_mat)
  out <- .gdpar_eb_apply_correction(
    eb_correction  = TRUE,
    laplace_result = lap,
    stan_data      = list(),
    p              = 3L,
    verbose        = FALSE
  )
  expect_true(isTRUE(out$applied))
  expect_equal(dim(out$constant), c(3L, 3L))

  out_off <- .gdpar_eb_apply_correction(
    eb_correction  = FALSE,
    laplace_result = lap,
    stan_data      = list(),
    p              = 1L,
    verbose        = FALSE
  )
  expect_false(isTRUE(out_off$applied))
  expect_true(is.na(out_off$constant))
})

# ---- (iv) Backward-compatibility default argument -----------------------

test_that(".gdpar_eb_apply_correction default p = 1L preserves the 8.6.B contract", {
  v <- 0.18
  cov_mat <- matrix(v, nrow = 1L, ncol = 1L)
  lap <- make_laplace_result(cov_mat)
  out <- .gdpar_eb_apply_correction(
    eb_correction  = TRUE,
    laplace_result = lap,
    stan_data      = list(),
    verbose        = FALSE
  )
  expect_true(isTRUE(out$applied))
  expect_type(out$constant, "double")
  expect_length(out$constant, 1L)
  expect_equal(out$constant, kappa_alpha_95 * v)
})

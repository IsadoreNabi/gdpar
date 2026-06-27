# Tests for sub-phase 8.3.8 (2026-05-22): B-splines W basis (D-D3).
# Structural coverage (no Stan fit, no env-var gating) for:
#   1. W_basis() constructor extension with boundary_knots + default
#      degree = 3 for bspline (decision D3 = a + b + c).
#   2. Internal helpers .gdpar_resolve_bspline_knots,
#      .gdpar_bspline_knots_full, .gdpar_estimate_theta_ref_proxy_range,
#      .gdpar_validate_bspline_boundary_range,
#      .gdpar_resolve_W_stan_data (decisions D2 + D4 + D5).
#   3. stan_data shape for bspline across the three regimes (K = 1 + p
#      = 1; K = 1 + p > 1; K > 1 + p = 1).
#   4. Generated Stan source contains apply_W_basis_diff +
#      bspline_basis_eval + W_type_id (decision D5 = 3 helper-based
#      dispatch).
#
# Smoke fits with Stan compilation are gated and live in
# tests/testthat/test-fit_K_smoke.R.

# ---------------------------------------------------------------------
# Section 1: W_basis() constructor (decision D3)
# ---------------------------------------------------------------------
test_that("W_basis(type = 'bspline') defaults degree to 3", {
  wb <- W_basis(type = "bspline", knots = c(0.3, 0.7),
                boundary_knots = c(0, 1))
  expect_equal(wb$degree, 3L)
})

test_that("W_basis(type = 'polynomial') still defaults degree to 1", {
  wb <- W_basis(type = "polynomial")
  expect_equal(wb$degree, 1L)
})

test_that("W_basis(type = 'bspline') stores boundary_knots field", {
  wb <- W_basis(type = "bspline", df = 5, boundary_knots = c(-2, 2))
  expect_equal(wb$boundary_knots, c(-2, 2))
})

test_that("W_basis(type = 'bspline') rejects malformed boundary_knots", {
  expect_error(
    W_basis(type = "bspline", knots = c(0.5),
            boundary_knots = c(1, 0)),
    class = "gdpar_input_error"
  )
  expect_error(
    W_basis(type = "bspline", knots = c(0.5),
            boundary_knots = c(0)),
    class = "gdpar_input_error"
  )
  expect_error(
    W_basis(type = "bspline", knots = c(0.5),
            boundary_knots = c(NA, 1)),
    class = "gdpar_input_error"
  )
})

test_that("W_basis(type = 'bspline') rejects knots outside boundary", {
  expect_error(
    W_basis(type = "bspline", knots = c(-0.5, 0.5),
            boundary_knots = c(0, 1)),
    class = "gdpar_input_error"
  )
  expect_error(
    W_basis(type = "bspline", knots = c(0.5, 1.5),
            boundary_knots = c(0, 1)),
    class = "gdpar_input_error"
  )
})

test_that("W_basis(type = 'bspline') df = K + degree matches dim", {
  wb <- W_basis(type = "bspline", degree = 3, knots = c(0.3, 0.7),
                boundary_knots = c(0, 1), p = 1L)
  expect_equal(wb$dim, 5L)
})

# ---------------------------------------------------------------------
# Section 2: knot resolver + augmented knot vector
# ---------------------------------------------------------------------
test_that(".gdpar_resolve_bspline_knots is passthrough when knots set", {
  wb <- W_basis(type = "bspline", degree = 3, knots = c(0.2, 0.5, 0.8),
                boundary_knots = c(0, 1))
  out <- .gdpar_resolve_bspline_knots(wb)
  expect_equal(out$knots, c(0.2, 0.5, 0.8))
})

test_that(".gdpar_resolve_bspline_knots derives equally-spaced knots from df", {
  wb <- W_basis(type = "bspline", degree = 3, df = 5,
                boundary_knots = c(0, 1))
  out <- .gdpar_resolve_bspline_knots(wb)
  expect_equal(length(out$knots), 2L)
  expect_equal(out$knots, c(1 / 3, 2 / 3))
})

test_that(".gdpar_resolve_bspline_knots returns numeric(0) when n_int = 0", {
  wb <- W_basis(type = "bspline", degree = 3, df = 3,
                boundary_knots = c(0, 1))
  out <- .gdpar_resolve_bspline_knots(wb)
  expect_equal(length(out$knots), 0L)
})

test_that(".gdpar_resolve_bspline_knots errors on df < degree", {
  wb <- W_basis(type = "bspline", degree = 3, df = 2,
                boundary_knots = c(0, 1))
  expect_error(.gdpar_resolve_bspline_knots(wb),
               class = "gdpar_input_error")
})

test_that(".gdpar_bspline_knots_full constructs augmented vector", {
  wb <- W_basis(type = "bspline", degree = 3, knots = c(0.3, 0.7),
                boundary_knots = c(0, 1))
  out <- .gdpar_bspline_knots_full(wb)
  expect_equal(out, c(0, 0, 0, 0, 0.3, 0.7, 1, 1, 1, 1))
})

test_that(".gdpar_bspline_knots_full passes through df resolution", {
  wb <- W_basis(type = "bspline", degree = 3, df = 5,
                boundary_knots = c(-2, 2))
  out <- .gdpar_bspline_knots_full(wb)
  expect_equal(out, c(-2, -2, -2, -2, -2 / 3, 2 / 3, 2, 2, 2, 2))
})

test_that(".gdpar_bspline_knots_full returns numeric(0) for polynomial", {
  wb <- W_basis(type = "polynomial", degree = 2)
  expect_equal(length(.gdpar_bspline_knots_full(wb)), 0L)
})

# ---------------------------------------------------------------------
# Section 3: proxy range + boundary validator (decision D4)
# ---------------------------------------------------------------------
test_that(".gdpar_estimate_theta_ref_proxy_range for gaussian identity", {
  fam <- gdpar_family("gaussian")
  rng <- .gdpar_estimate_theta_ref_proxy_range(c(0, 1, 2, 3, 4), fam)
  expect_equal(length(rng), 2L)
  expect_true(rng[1L] < 0 && rng[2L] > 4)
})

test_that(".gdpar_estimate_theta_ref_proxy_range for beta logit", {
  fam <- gdpar_family("beta")
  rng <- .gdpar_estimate_theta_ref_proxy_range(c(0.1, 0.5, 0.9), fam)
  expect_equal(length(rng), 2L)
  expect_true(all(is.finite(rng)))
})

test_that(".gdpar_estimate_theta_ref_proxy_range for gamma log", {
  fam <- gdpar_family("gamma")
  rng <- .gdpar_estimate_theta_ref_proxy_range(c(0.5, 1, 2, 4), fam)
  expect_equal(length(rng), 2L)
  expect_true(all(is.finite(rng)))
})

test_that(".gdpar_validate_bspline_boundary_range accepts contained range", {
  wb <- W_basis(type = "bspline", degree = 3, df = 5,
                boundary_knots = c(-10, 10))
  expect_silent(.gdpar_validate_bspline_boundary_range(wb, c(-3, 3)))
})

test_that(".gdpar_validate_bspline_boundary_range rejects exceeding range", {
  wb <- W_basis(type = "bspline", degree = 3, knots = c(0.5),
                boundary_knots = c(0, 1))
  expect_error(
    .gdpar_validate_bspline_boundary_range(wb, c(-0.5, 1.5)),
    class = "gdpar_input_error"
  )
  expect_error(
    .gdpar_validate_bspline_boundary_range(wb, c(0, 1.5)),
    class = "gdpar_input_error"
  )
  expect_error(
    .gdpar_validate_bspline_boundary_range(wb, c(-0.5, 1)),
    class = "gdpar_input_error"
  )
})

test_that(".gdpar_validate_bspline_boundary_range passthrough for polynomial", {
  wb <- W_basis(type = "polynomial", degree = 2)
  expect_silent(.gdpar_validate_bspline_boundary_range(wb, c(-100, 100)))
})

# ---------------------------------------------------------------------
# Section 4: W stan_data resolver (decision D5)
# ---------------------------------------------------------------------
test_that(".gdpar_resolve_W_stan_data: off path emits empty fields", {
  fam <- gdpar_family("gaussian")
  out <- .gdpar_resolve_W_stan_data(NULL, 0L, rnorm(10), fam)
  expect_equal(out$W_type_id, 0L)
  expect_equal(out$W_n_knots_full, 0L)
  expect_equal(length(out$W_knots_full), 0L)
  expect_equal(out$W_degree, 0L)
})

test_that(".gdpar_resolve_W_stan_data: polynomial path emits id = 1", {
  fam <- gdpar_family("gaussian")
  wb <- W_basis(type = "polynomial", degree = 2)
  out <- .gdpar_resolve_W_stan_data(wb, 1L, rnorm(10), fam)
  expect_equal(out$W_type_id, 1L)
  expect_equal(out$W_n_knots_full, 0L)
  expect_equal(out$W_degree, 2L)
})

test_that(".gdpar_resolve_W_stan_data: bspline path emits id = 2 + knots_full", {
  fam <- gdpar_family("gaussian")
  wb <- W_basis(type = "bspline", degree = 3, df = 5,
                boundary_knots = c(-10, 10))
  out <- .gdpar_resolve_W_stan_data(wb, 1L, rnorm(10), fam)
  expect_equal(out$W_type_id, 2L)
  expect_equal(out$W_n_knots_full, 10L)
  expect_equal(out$W_degree, 3L)
  expect_equal(out$W_knots_full[1L], -10)
  expect_equal(out$W_knots_full[10L], 10)
})

test_that(".gdpar_resolve_W_stan_data: bspline without boundary_knots aborts", {
  fam <- gdpar_family("gaussian")
  wb <- W_basis(type = "bspline", degree = 3, knots = c(0.5))
  expect_error(.gdpar_resolve_W_stan_data(wb, 1L, rnorm(10), fam),
               class = "gdpar_input_error")
})

test_that(".gdpar_resolve_W_stan_data: bspline triggers boundary validator", {
  fam <- gdpar_family("gaussian")
  wb <- W_basis(type = "bspline", degree = 3, df = 5,
                boundary_knots = c(0, 1))
  expect_error(
    .gdpar_resolve_W_stan_data(wb, 1L, c(-5, 5, 10), fam),
    class = "gdpar_input_error"
  )
})

# ---------------------------------------------------------------------
# Section 5: stan_data shape for bspline across the three regimes
# ---------------------------------------------------------------------
test_that("assemble_stan_data K = 1 + p = 1 with bspline populates W_* fields", {
  set.seed(1L)
  n <- 30L
  d_df <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  amm <- amm_spec(
    a = ~ x1,
    W = W_basis(type = "bspline", degree = 3, df = 5,
                boundary_knots = c(-10, 10))
  )
  amm$W <- materialize_W_basis(amm$W, 1L)
  design <- list(
    Z_a = matrix(rnorm(n * 1L), n, 1L),
    Z_b = matrix(0, n, 0L),
    X   = matrix(rnorm(n * 2L), n, 2L)
  )
  y <- rnorm(n)
  fam <- gdpar_family("gaussian")
  sd <- assemble_stan_data(design, fam, amm, y, theta_anchor = 0)
  expect_equal(sd$W_type_id, 2L)
  expect_equal(sd$W_degree, 3L)
  expect_equal(sd$W_n_knots_full, 10L)
  expect_equal(length(sd$W_knots_full), 10L)
})

test_that("assemble_stan_data K = 1 + p = 1 with polynomial keeps W_type_id = 1", {
  set.seed(1L)
  n <- 30L
  amm <- amm_spec(a = ~ x1, W = W_basis(type = "polynomial", degree = 2))
  amm$W <- materialize_W_basis(amm$W, 1L)
  design <- list(
    Z_a = matrix(rnorm(n * 1L), n, 1L),
    Z_b = matrix(0, n, 0L),
    X   = matrix(rnorm(n * 1L), n, 1L)
  )
  fam <- gdpar_family("gaussian")
  sd <- assemble_stan_data(design, fam, amm, rnorm(n), theta_anchor = 0)
  expect_equal(sd$W_type_id, 1L)
  expect_equal(sd$W_degree, 2L)
  expect_equal(sd$W_n_knots_full, 0L)
})

# ---------------------------------------------------------------------
# Section 6: generated Stan source includes the new helper + data fields
# ---------------------------------------------------------------------
test_that("generate_stan_code emits apply_W_basis_diff + bspline_basis_eval + W_type_id", {
  prior <- gdpar_prior()
  src <- generate_stan_code(prior)
  expect_true(grepl("apply_W_basis_diff", src, fixed = TRUE))
  expect_true(grepl("bspline_basis_eval", src, fixed = TRUE))
  expect_true(grepl("W_type_id", src, fixed = TRUE))
  expect_true(grepl("W_knots_full", src, fixed = TRUE))
  expect_false(grepl("W_is_polynomial", src, fixed = TRUE))
  expect_false(grepl("{{", src, fixed = TRUE))
})

test_that("generate_stan_code_multi emits the same helpers + fields", {
  prior <- gdpar_prior()
  src <- generate_stan_code_multi(prior)
  expect_true(grepl("apply_W_basis_diff", src, fixed = TRUE))
  expect_true(grepl("bspline_basis_eval", src, fixed = TRUE))
  expect_true(grepl("W_type_id", src, fixed = TRUE))
  expect_true(grepl("W_knots_full", src, fixed = TRUE))
  expect_false(grepl("W_is_polynomial", src, fixed = TRUE))
})

test_that("generate_stan_code_K emits the helpers in functions block", {
  prior <- gdpar_prior()
  src <- generate_stan_code_K(prior, cp_a_per_K = rep(FALSE, 2L))
  expect_true(grepl("apply_W_basis_diff", src, fixed = TRUE))
  expect_true(grepl("bspline_basis_eval", src, fixed = TRUE))
  expect_true(grepl("apply_inv_link_by_id", src, fixed = TRUE))
  expect_true(grepl("W_type_id", src, fixed = TRUE))
  expect_false(grepl("W_is_polynomial", src, fixed = TRUE))
})

test_that("amm_spec records the implied AMM level", {
  expect_equal(amm_spec()$level, 0L)
  expect_equal(amm_spec(a = ~ x1)$level, 1L)
  expect_equal(amm_spec(a = ~ x1, b = ~ x1)$level, 2L)
  expect_equal(amm_spec(W = W_basis("polynomial", 1))$level, 2L)
})

test_that("amm_spec rejects two-sided formulas", {
  expect_error(
    amm_spec(a = y ~ x1),
    class = "gdpar_input_error"
  )
})

test_that("amm_spec accepts NULL components", {
  spec <- amm_spec(a = NULL, b = NULL, W = NULL)
  expect_null(spec$a)
  expect_null(spec$b)
  expect_null(spec$W)
})

test_that("build_amm_design centers Z_a columns", {
  set.seed(NULL)
  df <- make_gaussian_data(n = 100, seed = 42)
  spec <- amm_spec(a = ~ x1 + x2)
  design <- gdpar:::build_amm_design(spec, df, formula_rhs = ~ x1 + x2)
  expect_equal(colMeans(design$Z_a),
               stats::setNames(c(0, 0), c("x1", "x2")),
               tolerance = 1e-12)
})

test_that("build_amm_design centers Z_b and standardizes X", {
  df <- make_gaussian_data(n = 100, seed = 42)
  spec <- amm_spec(a = ~ x1, b = ~ x1, W = W_basis("polynomial", 1))
  spec$W <- gdpar:::materialize_W_basis(spec$W, p = 1L)
  design <- gdpar:::build_amm_design(spec, df, formula_rhs = ~ x1 + x2)
  expect_equal(unname(colMeans(design$Z_b)), 0, tolerance = 1e-12)
  expect_equal(unname(apply(design$X, 2L, stats::sd)),
               c(1, 1), tolerance = 1e-12)
})

test_that("build_amm_design rejects constant covariates in W", {
  df <- make_gaussian_data(n = 50, seed = 42)
  df$x_const <- 1
  spec <- amm_spec(W = W_basis("polynomial", 1), x_vars = "x_const")
  spec$W <- gdpar:::materialize_W_basis(spec$W, p = 1L)
  expect_error(
    gdpar:::build_amm_design(spec, df, formula_rhs = ~ x_const),
    class = "gdpar_input_error"
  )
})

test_that("print method for amm_spec runs without error", {
  expect_output(print(amm_spec(a = ~ x1)), "amm_spec")
})

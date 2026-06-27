.make_mv_data <- function(n = 20L, seed = 1L) {
  set.seed(seed)
  data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n),
    x3 = rnorm(n)
  )
}

test_that("build_amm_design() multi: uniform dimwise builds equal-shape Z_a_list", {
  d <- .make_mv_data()
  spec <- amm_spec(p = 2L, dims = dimwise(a = ~ x1 + x2, b = ~ x1))
  out <- build_amm_design(spec, d, formula_rhs = ~ x1 + x2)
  expect_named(out, c("p", "Z_a_list", "Z_b_list", "X",
                       "Z_a_means_list", "Z_b_means_list",
                       "X_means", "X_sds",
                       "Z_a_names_list", "Z_b_names_list", "X_names"))
  expect_equal(out$p, 2L)
  expect_length(out$Z_a_list, 2L)
  expect_length(out$Z_b_list, 2L)
  for (k in seq_len(2L)) {
    expect_equal(dim(out$Z_a_list[[k]]), c(20L, 2L))
    expect_equal(dim(out$Z_b_list[[k]]), c(20L, 1L))
  }
})

test_that("build_amm_design() multi: per-k overrides yield ragged column counts", {
  d <- .make_mv_data()
  dims <- override(
    dimwise(a = ~ x1, b = NULL),
    k = 2L, a = ~ x1 + x2 + x3, b = ~ x1 + x2
  )
  spec <- amm_spec(p = 2L, dims = dims)
  out <- build_amm_design(spec, d, formula_rhs = ~ x1 + x2 + x3)
  expect_equal(ncol(out$Z_a_list[[1L]]), 1L)
  expect_equal(ncol(out$Z_a_list[[2L]]), 3L)
  expect_equal(ncol(out$Z_b_list[[1L]]), 0L)
  expect_equal(ncol(out$Z_b_list[[2L]]), 2L)
})

test_that("build_amm_design() multi: NULL components yield empty matrices", {
  d <- .make_mv_data()
  spec <- amm_spec(p = 3L, dims = dimwise(a = NULL, b = NULL),
                   W = W_basis(type = "polynomial", degree = 2),
                   x_vars = c("x1", "x2"))
  out <- build_amm_design(spec, d, formula_rhs = ~ x1 + x2)
  expect_length(out$Z_a_list, 3L)
  for (k in seq_len(3L)) {
    expect_equal(ncol(out$Z_a_list[[k]]), 0L)
    expect_equal(ncol(out$Z_b_list[[k]]), 0L)
    expect_equal(nrow(out$Z_a_list[[k]]), 20L)
  }
})

test_that("build_amm_design() multi: Z_a_list columns are centered (colMeans == 0)", {
  d <- .make_mv_data()
  spec <- amm_spec(p = 2L, dims = dimwise(a = ~ x1 + x2, b = ~ x1))
  out <- build_amm_design(spec, d, formula_rhs = ~ x1 + x2)
  for (k in seq_len(2L)) {
    cm <- colMeans(out$Z_a_list[[k]])
    expect_equal(unname(cm), rep(0, length(cm)), tolerance = 1e-12)
    cm_b <- colMeans(out$Z_b_list[[k]])
    expect_equal(unname(cm_b), rep(0, length(cm_b)), tolerance = 1e-12)
  }
})

test_that("build_amm_design() multi: means and names are recorded per coordinate", {
  d <- .make_mv_data()
  spec <- amm_spec(p = 2L, dims = dimwise(a = ~ x1 + x2, b = ~ x1))
  out <- build_amm_design(spec, d, formula_rhs = ~ x1 + x2)
  expect_length(out$Z_a_means_list, 2L)
  expect_length(out$Z_a_names_list, 2L)
  for (k in seq_len(2L)) {
    expect_length(out$Z_a_means_list[[k]], 2L)
    expect_equal(out$Z_a_names_list[[k]], c("x1", "x2"))
    expect_equal(out$Z_b_names_list[[k]], "x1")
  }
})

test_that("build_amm_design() multi: X is shared across coordinates and centered+scaled", {
  d <- .make_mv_data()
  spec <- amm_spec(p = 2L, dims = dimwise(a = ~ x1, b = NULL),
                   W = W_basis(type = "polynomial", degree = 2),
                   x_vars = c("x1", "x2"))
  out <- build_amm_design(spec, d, formula_rhs = ~ x1 + x2)
  expect_equal(dim(out$X), c(20L, 2L))
  expect_equal(colMeans(out$X), c(x1 = 0, x2 = 0), tolerance = 1e-12)
  expect_equal(apply(out$X, 2L, sd), c(x1 = 1, x2 = 1), tolerance = 1e-12)
  expect_equal(out$X_names, c("x1", "x2"))
})

test_that("build_amm_design() multi: no W -> empty X with zero columns", {
  d <- .make_mv_data()
  spec <- amm_spec(p = 2L, dims = dimwise(a = ~ x1, b = NULL))
  out <- build_amm_design(spec, d, formula_rhs = ~ x1)
  expect_equal(dim(out$X), c(20L, 0L))
  expect_length(out$X_means, 0L)
  expect_length(out$X_sds, 0L)
})

test_that("build_amm_design() multi: NA in needed covariate aborts", {
  d <- .make_mv_data()
  d$x1[3L] <- NA
  spec <- amm_spec(p = 2L, dims = dimwise(a = ~ x1, b = NULL))
  expect_error(
    build_amm_design(spec, d, formula_rhs = ~ x1),
    class = "gdpar_input_error"
  )
})

test_that("build_amm_design() multi: missing x_vars in data aborts when W is active", {
  d <- .make_mv_data()
  spec <- amm_spec(p = 2L, dims = dimwise(a = ~ x1, b = NULL),
                   W = W_basis(type = "polynomial", degree = 2),
                   x_vars = c("x1", "x_missing"))
  expect_error(
    build_amm_design(spec, d, formula_rhs = ~ x1),
    class = "gdpar_input_error"
  )
})

test_that("build_amm_design() multi: constant x_var aborts when W is active", {
  d <- .make_mv_data()
  d$x_const <- 1
  spec <- amm_spec(p = 2L, dims = dimwise(a = ~ x1, b = NULL),
                   W = W_basis(type = "polynomial", degree = 2),
                   x_vars = c("x1", "x_const"))
  expect_error(
    build_amm_design(spec, d, formula_rhs = ~ x1 + x_const),
    class = "gdpar_input_error"
  )
})

test_that("build_amm_design() multi: dispatch happens only for p > 1L (scalar untouched)", {
  d <- .make_mv_data()
  spec <- amm_spec(a = ~ x1, b = NULL)
  out <- build_amm_design(spec, d, formula_rhs = ~ x1)
  expect_named(out, c("Z_a", "Z_b", "X",
                       "Z_a_means", "Z_b_means",
                       "X_means", "X_sds",
                       "Z_a_names", "Z_b_names", "X_names"))
  expect_false("Z_a_list" %in% names(out))
})

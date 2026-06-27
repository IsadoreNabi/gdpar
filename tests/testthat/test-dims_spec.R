test_that("dimwise returns a dims_spec with base components", {
  ds <- dimwise(a = ~ x1 + x2, b = ~ x1)
  expect_s3_class(ds, "dims_spec")
  expect_equal(deparse(ds$base$a), "~x1 + x2")
  expect_equal(deparse(ds$base$b), "~x1")
  expect_length(ds$overrides, 0L)
})

test_that("dimwise accepts NULL components", {
  ds <- dimwise(a = NULL, b = NULL)
  expect_s3_class(ds, "dims_spec")
  expect_null(ds$base$a)
  expect_null(ds$base$b)
})

test_that("dimwise rejects two-sided formulas", {
  expect_error(
    dimwise(a = y ~ x1),
    class = "gdpar_input_error"
  )
  expect_error(
    dimwise(b = y ~ x1),
    class = "gdpar_input_error"
  )
})

test_that("override registers a per-dimension override on dims_spec", {
  base <- dimwise(a = ~ x1 + x2, b = ~ x1)
  ov <- override(base, k = 2L, a = ~ x1)
  expect_s3_class(ov, "dims_spec")
  expect_length(ov$overrides, 1L)
  expect_true("2" %in% names(ov$overrides))
  expect_equal(deparse(ov$overrides[["2"]]$a), "~x1")
  expect_true(isTRUE(ov$overrides[["2"]]$a_set))
  expect_false(isTRUE(ov$overrides[["2"]]$b_set))
})

test_that("override distinguishes missing argument from explicit NULL", {
  base <- dimwise(a = ~ x1 + x2, b = ~ x1)
  ov_inherit <- override(base, k = 2L, a = ~ x1)
  expect_false(isTRUE(ov_inherit$overrides[["2"]]$b_set))
  ov_disable <- override(base, k = 2L, b = NULL)
  expect_true(isTRUE(ov_disable$overrides[["2"]]$b_set))
  expect_null(ov_disable$overrides[["2"]]$b)
})

test_that("override calls compose and overwrite per-k", {
  base <- dimwise(a = ~ x1 + x2, b = ~ x1)
  v1 <- override(base, k = 2L, a = ~ x1)
  v2 <- override(v1, k = 2L, a = ~ x2)
  expect_length(v2$overrides, 1L)
  expect_equal(deparse(v2$overrides[["2"]]$a), "~x2")
  v3 <- override(v2, k = 3L, b = NULL)
  expect_length(v3$overrides, 2L)
  expect_true("3" %in% names(v3$overrides))
})

test_that("override rejects calls with neither a nor b", {
  base <- dimwise(a = ~ x1, b = ~ x1)
  expect_error(
    override(base, k = 2L),
    class = "gdpar_input_error"
  )
})

test_that("override rejects non-positive integer k", {
  base <- dimwise(a = ~ x1, b = ~ x1)
  expect_error(override(base, k = 0L, a = ~ x2), class = "gdpar_input_error")
  expect_error(override(base, k = -1L, a = ~ x2), class = "gdpar_input_error")
  expect_error(override(base, k = 1.5, a = ~ x2), class = "gdpar_input_error")
})

test_that("override rejects non-dims_spec input", {
  expect_error(
    override(list(), k = 1L, a = ~ x1),
    class = "gdpar_input_error"
  )
})

test_that("resolve_dims_spec produces correct per-k list", {
  base <- dimwise(a = ~ x1 + x2, b = ~ x1)
  resolved <- gdpar:::resolve_dims_spec(base, p = 3L)
  expect_length(resolved, 3L)
  for (k in seq_len(3L)) {
    expect_equal(deparse(resolved[[k]]$a), "~x1 + x2")
    expect_equal(deparse(resolved[[k]]$b), "~x1")
  }
})

test_that("resolve_dims_spec applies overrides on top of base", {
  base <- dimwise(a = ~ x1 + x2, b = ~ x1)
  spec <- override(base, k = 2L, a = ~ x1)
  spec <- override(spec, k = 3L, b = NULL)
  resolved <- gdpar:::resolve_dims_spec(spec, p = 3L)
  expect_equal(deparse(resolved[[1L]]$a), "~x1 + x2")
  expect_equal(deparse(resolved[[1L]]$b), "~x1")
  expect_equal(deparse(resolved[[2L]]$a), "~x1")
  expect_equal(deparse(resolved[[2L]]$b), "~x1")
  expect_equal(deparse(resolved[[3L]]$a), "~x1 + x2")
  expect_null(resolved[[3L]]$b)
})

test_that("resolve_dims_spec aborts on out-of-range override", {
  base <- dimwise(a = ~ x1, b = ~ x1)
  spec <- override(base, k = 5L, a = ~ x2)
  expect_error(
    gdpar:::resolve_dims_spec(spec, p = 3L),
    class = "gdpar_input_error"
  )
})

test_that("resolve_dims_spec accepts overrides at boundary k = p", {
  base <- dimwise(a = ~ x1, b = ~ x1)
  spec <- override(base, k = 3L, a = ~ x2)
  resolved <- gdpar:::resolve_dims_spec(spec, p = 3L)
  expect_equal(deparse(resolved[[3L]]$a), "~x2")
})

test_that("print method for dims_spec runs without error", {
  base <- dimwise(a = ~ x1, b = ~ x1)
  expect_output(print(base), "dims_spec")
  expect_output(print(base), "overrides: <none>")
  spec <- override(base, k = 2L, a = ~ x1 + x2)
  expect_output(print(spec), "k = 2")
})

test_that("dimwise and override are pipe-friendly", {
  spec <- dimwise(a = ~ x1, b = ~ x1) |>
    override(k = 2L, a = ~ x1 + x2) |>
    override(k = 3L, b = NULL)
  expect_s3_class(spec, "dims_spec")
  expect_length(spec$overrides, 2L)
})

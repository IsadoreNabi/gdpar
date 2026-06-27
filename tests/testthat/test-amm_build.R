test_that("amm_build initialises with default p = 1 and empty state", {
  b <- amm_build()
  expect_s3_class(b, "amm_builder")
  expect_identical(b$p, 1L)
  expect_s3_class(b$dims, "dims_spec")
  expect_null(b$dims$base$a)
  expect_null(b$dims$base$b)
  expect_length(b$dims$overrides, 0L)
  expect_null(b$W)
  expect_null(b$x_vars)
})

test_that("amm_build accepts a custom positive integer p", {
  b <- amm_build(p = 3L)
  expect_identical(b$p, 3L)
})

test_that("amm_build rejects non-positive or non-integer p", {
  expect_error(amm_build(p = 0L), class = "gdpar_input_error")
  expect_error(amm_build(p = -1L), class = "gdpar_input_error")
  expect_error(amm_build(p = 1.5), class = "gdpar_input_error")
  expect_error(amm_build(p = "two"), class = "gdpar_input_error")
})

test_that("amm_set_a_uniform mutates base$a and preserves overrides", {
  b <- amm_build(p = 3L) |>
    amm_set_a(k = 2L, ~ x1) |>
    amm_set_a_uniform(~ x1 + x2)
  expect_equal(deparse(b$dims$base$a), "~x1 + x2")
  expect_equal(deparse(b$dims$overrides[["2"]][["a"]]), "~x1")
  expect_true(isTRUE(b$dims$overrides[["2"]]$a_set))
})

test_that("amm_set_a_uniform accepts NULL to clear the base", {
  b <- amm_build() |>
    amm_set_a_uniform(~ x1) |>
    amm_set_a_uniform(NULL)
  expect_null(b$dims$base$a)
})

test_that("amm_set_a_uniform rejects two-sided formulas", {
  expect_error(
    amm_build() |> amm_set_a_uniform(y ~ x1),
    class = "gdpar_input_error"
  )
})

test_that("amm_set_a_uniform rejects non-builder first argument", {
  expect_error(
    amm_set_a_uniform(list(), ~ x1),
    class = "gdpar_input_error"
  )
})

test_that("amm_set_b_uniform mutates base$b and preserves overrides", {
  b <- amm_build(p = 2L) |>
    amm_set_b(k = 1L, ~ x2) |>
    amm_set_b_uniform(~ x1)
  expect_equal(deparse(b$dims$base$b), "~x1")
  expect_equal(deparse(b$dims$overrides[["1"]][["b"]]), "~x2")
})

test_that("amm_set_a registers a per-k override", {
  b <- amm_build(p = 3L) |>
    amm_set_a(k = 2L, ~ x1)
  ov <- b$dims$overrides[["2"]]
  expect_true(isTRUE(ov$a_set))
  expect_equal(deparse(ov[["a"]]), "~x1")
})

test_that("amm_set_a with NULL disables additive component for k only", {
  b <- amm_build(p = 3L) |>
    amm_set_a_uniform(~ x1 + x2) |>
    amm_set_a(k = 2L, NULL)
  resolved <- gdpar:::resolve_dims_spec(b$dims, 3L)
  expect_equal(deparse(resolved[[1L]]$a), "~x1 + x2")
  expect_null(resolved[[2L]]$a)
  expect_equal(deparse(resolved[[3L]]$a), "~x1 + x2")
})

test_that("amm_set_a replaces a previous override on the same k", {
  b <- amm_build(p = 2L) |>
    amm_set_a(k = 1L, ~ x1) |>
    amm_set_a(k = 1L, ~ x2)
  expect_equal(deparse(b$dims$overrides[["1"]][["a"]]), "~x2")
})

test_that("amm_set_a rejects k out of range and non-positive k", {
  expect_error(
    amm_build(p = 2L) |> amm_set_a(k = 3L, ~ x1),
    class = "gdpar_input_error"
  )
  expect_error(
    amm_build(p = 2L) |> amm_set_a(k = 0L, ~ x1),
    class = "gdpar_input_error"
  )
})

test_that("amm_set_b registers per-k override and supports NULL", {
  b <- amm_build(p = 2L) |>
    amm_set_b_uniform(~ x1) |>
    amm_set_b(k = 2L, NULL)
  resolved <- gdpar:::resolve_dims_spec(b$dims, 2L)
  expect_equal(deparse(resolved[[1L]]$b), "~x1")
  expect_null(resolved[[2L]]$b)
})

test_that("amm_set_W stores a W_basis and accepts NULL to clear", {
  wb <- W_basis(type = "polynomial", degree = 2)
  b <- amm_build() |> amm_set_W(wb)
  expect_s3_class(b$W, "W_basis")
  b2 <- amm_set_W(b, NULL)
  expect_null(b2$W)
})

test_that("amm_set_W rejects non-W_basis input", {
  expect_error(
    amm_build() |> amm_set_W(list()),
    class = "gdpar_input_error"
  )
})

test_that("amm_set_x_vars stores a character vector and accepts NULL", {
  b <- amm_build() |> amm_set_x_vars(c("x1", "x2"))
  expect_identical(b$x_vars, c("x1", "x2"))
  b2 <- amm_set_x_vars(b, NULL)
  expect_null(b2$x_vars)
})

test_that("amm_set_x_vars rejects empty or non-character input", {
  expect_error(
    amm_build() |> amm_set_x_vars(character(0)),
    class = "gdpar_input_error"
  )
  expect_error(
    amm_build() |> amm_set_x_vars(1:3),
    class = "gdpar_input_error"
  )
})

test_that("as_amm_spec on empty builder returns Level 0 scalar spec", {
  spec <- amm_build() |> as_amm_spec()
  expect_s3_class(spec, "amm_spec")
  expect_equal(spec$level, 0L)
  expect_identical(spec$p, 1L)
  expect_null(spec$a)
  expect_null(spec$b)
  expect_null(spec$W)
  expect_null(spec$dims)
})

test_that("as_amm_spec on scalar builder with a only yields Level 1", {
  spec <- amm_build() |>
    amm_set_a_uniform(~ x1 + x2) |>
    as_amm_spec()
  expect_equal(spec$level, 1L)
  expect_equal(deparse(spec$a), "~x1 + x2")
  expect_null(spec$b)
})

test_that("as_amm_spec on scalar builder with a + b yields Level 2", {
  spec <- amm_build() |>
    amm_set_a_uniform(~ x1) |>
    amm_set_b_uniform(~ x1) |>
    as_amm_spec()
  expect_equal(spec$level, 2L)
})

test_that("as_amm_spec on scalar builder with W only yields Level 2", {
  spec <- amm_build() |>
    amm_set_W(W_basis(type = "polynomial", degree = 1)) |>
    as_amm_spec()
  expect_equal(spec$level, 2L)
})

test_that("as_amm_spec on scalar builder applies k=1 override over base", {
  spec <- amm_build() |>
    amm_set_a_uniform(~ x1) |>
    amm_set_a(k = 1L, ~ x2) |>
    as_amm_spec()
  expect_equal(deparse(spec$a), "~x2")
})

test_that("as_amm_spec on multivariate builder broadcasts uniform a and b", {
  spec <- amm_build(p = 3L) |>
    amm_set_a_uniform(~ x1 + x2) |>
    amm_set_b_uniform(~ x1) |>
    as_amm_spec()
  expect_identical(spec$p, 3L)
  expect_length(spec$dims, 3L)
  for (k in seq_len(3L)) {
    expect_equal(deparse(spec$dims[[k]]$a), "~x1 + x2")
    expect_equal(deparse(spec$dims[[k]]$b), "~x1")
  }
})

test_that("as_amm_spec on multivariate builder applies per-k override", {
  spec <- amm_build(p = 3L) |>
    amm_set_a_uniform(~ x1 + x2) |>
    amm_set_a(k = 2L, ~ x1) |>
    as_amm_spec()
  expect_equal(deparse(spec$dims[[1L]]$a), "~x1 + x2")
  expect_equal(deparse(spec$dims[[2L]]$a), "~x1")
  expect_equal(deparse(spec$dims[[3L]]$a), "~x1 + x2")
})

test_that("as_amm_spec carries W and x_vars through to amm_spec", {
  wb <- W_basis(type = "polynomial", degree = 2)
  spec <- amm_build(p = 2L) |>
    amm_set_a_uniform(~ x1) |>
    amm_set_W(wb) |>
    amm_set_x_vars(c("x1", "x2")) |>
    as_amm_spec()
  expect_s3_class(spec$W, "W_basis")
  expect_identical(spec$x_vars, c("x1", "x2"))
})

test_that("as_amm_spec produces a spec equivalent to direct amm_spec call (scalar)", {
  spec_b <- amm_build() |>
    amm_set_a_uniform(~ x1 + x2) |>
    amm_set_b_uniform(~ x1) |>
    amm_set_W(W_basis(type = "polynomial", degree = 1)) |>
    as_amm_spec()
  spec_d <- amm_spec(
    a = ~ x1 + x2,
    b = ~ x1,
    W = W_basis(type = "polynomial", degree = 1)
  )
  expect_equal(spec_b$level, spec_d$level)
  expect_equal(deparse(spec_b$a), deparse(spec_d$a))
  expect_equal(deparse(spec_b$b), deparse(spec_d$b))
  expect_equal(spec_b$p, spec_d$p)
})

test_that("as_amm_spec produces a spec equivalent to direct amm_spec call (multivariate)", {
  spec_b <- amm_build(p = 2L) |>
    amm_set_a_uniform(~ x1 + x2) |>
    amm_set_a(k = 2L, ~ x1) |>
    amm_set_b_uniform(~ x1) |>
    as_amm_spec()
  spec_d <- amm_spec(
    p = 2L,
    dims = override(dimwise(a = ~ x1 + x2, b = ~ x1), k = 2L, a = ~ x1)
  )
  expect_equal(spec_b$level, spec_d$level)
  expect_equal(spec_b$p, spec_d$p)
  expect_equal(length(spec_b$dims), length(spec_d$dims))
  for (k in seq_along(spec_b$dims)) {
    expect_equal(deparse(spec_b$dims[[k]]$a), deparse(spec_d$dims[[k]]$a))
    expect_equal(deparse(spec_b$dims[[k]]$b), deparse(spec_d$dims[[k]]$b))
  }
})

test_that("as_amm_spec rejects non-builder input", {
  expect_error(as_amm_spec(list()), class = "gdpar_input_error")
})

test_that("print.amm_builder runs without error for scalar and multivariate", {
  expect_output(
    print(amm_build()),
    "amm_builder"
  )
  b_mv <- amm_build(p = 3L) |>
    amm_set_a_uniform(~ x1) |>
    amm_set_a(k = 2L, NULL) |>
    amm_set_W(W_basis(type = "polynomial", degree = 1)) |>
    amm_set_x_vars(c("x1", "x2"))
  expect_output(print(b_mv), "amm_builder")
  expect_output(print(b_mv), "overrides")
  expect_output(print(b_mv), "x_vars")
})

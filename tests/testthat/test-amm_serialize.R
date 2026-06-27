test_that("round-trip preserves an empty scalar spec (Level 0)", {
  spec <- amm_spec()
  path <- tempfile(fileext = ".gdpar")
  amm_save_spec(spec, path)
  spec2 <- amm_load_spec(path)
  expect_s3_class(spec2, "amm_spec")
  expect_equal(spec2$level, 0L)
  expect_identical(spec2$p, 1L)
  expect_null(spec2$a)
  expect_null(spec2$b)
  expect_null(spec2$W)
  expect_null(spec2$x_vars)
  expect_null(spec2$dims)
})

test_that("round-trip preserves a scalar spec with a and b", {
  spec <- amm_spec(a = ~ x1 + x2, b = ~ x1)
  path <- tempfile(fileext = ".gdpar")
  amm_save_spec(spec, path)
  spec2 <- amm_load_spec(path)
  expect_equal(spec2$level, spec$level)
  expect_equal(deparse(spec2$a), deparse(spec$a))
  expect_equal(deparse(spec2$b), deparse(spec$b))
})

test_that("round-trip preserves a scalar spec with W polynomial", {
  spec <- amm_spec(W = W_basis(type = "polynomial", degree = 3))
  path <- tempfile(fileext = ".gdpar")
  amm_save_spec(spec, path)
  spec2 <- amm_load_spec(path)
  expect_s3_class(spec2$W, "W_basis")
  expect_equal(spec2$W$type, "polynomial")
  expect_equal(spec2$W$degree, 3L)
})

test_that("round-trip preserves a scalar spec with W bspline (knots variant)", {
  spec <- amm_spec(W = W_basis(type = "bspline", degree = 3,
                               knots = c(0.25, 0.5, 0.75)))
  path <- tempfile(fileext = ".gdpar")
  amm_save_spec(spec, path)
  spec2 <- amm_load_spec(path)
  expect_equal(spec2$W$type, "bspline")
  expect_equal(spec2$W$degree, 3L)
  expect_equal(spec2$W$knots, c(0.25, 0.5, 0.75))
  expect_null(spec2$W$df)
})

test_that("round-trip preserves a scalar spec with W bspline (df variant)", {
  spec <- amm_spec(W = W_basis(type = "bspline", degree = 3, df = 5))
  path <- tempfile(fileext = ".gdpar")
  amm_save_spec(spec, path)
  spec2 <- amm_load_spec(path)
  expect_equal(spec2$W$type, "bspline")
  expect_equal(spec2$W$degree, 3L)
  expect_null(spec2$W$knots)
  expect_equal(spec2$W$df, 5L)
})

test_that("round-trip preserves x_vars on a scalar spec", {
  spec <- amm_spec(a = ~ x1, W = W_basis(type = "polynomial", degree = 1),
                   x_vars = c("x1", "x2"))
  path <- tempfile(fileext = ".gdpar")
  amm_save_spec(spec, path)
  spec2 <- amm_load_spec(path)
  expect_identical(spec2$x_vars, c("x1", "x2"))
})

test_that("round-trip preserves a multivariate spec with uniform a and b", {
  spec <- amm_spec(p = 3L,
                   dims = dimwise(a = ~ x1 + x2, b = ~ x1),
                   W = W_basis(type = "polynomial", degree = 2))
  path <- tempfile(fileext = ".gdpar")
  amm_save_spec(spec, path)
  spec2 <- amm_load_spec(path)
  expect_identical(spec2$p, 3L)
  expect_length(spec2$dims, 3L)
  for (k in seq_len(3L)) {
    expect_equal(deparse(spec2$dims[[k]]$a), "~x1 + x2")
    expect_equal(deparse(spec2$dims[[k]]$b), "~x1")
  }
  expect_equal(spec2$W$type, "polynomial")
  expect_equal(spec2$W$degree, 2L)
})

test_that("round-trip preserves a multivariate spec with overrides and NULLs", {
  spec <- amm_spec(p = 3L,
                   dims = override(
                     override(dimwise(a = ~ x1 + x2, b = ~ x1),
                              k = 2L, a = ~ x1),
                     k = 3L, b = NULL))
  path <- tempfile(fileext = ".gdpar")
  amm_save_spec(spec, path)
  spec2 <- amm_load_spec(path)
  expect_equal(deparse(spec2$dims[[1L]]$a), "~x1 + x2")
  expect_equal(deparse(spec2$dims[[2L]]$a), "~x1")
  expect_equal(deparse(spec2$dims[[3L]]$a), "~x1 + x2")
  expect_equal(deparse(spec2$dims[[1L]]$b), "~x1")
  expect_equal(deparse(spec2$dims[[2L]]$b), "~x1")
  expect_null(spec2$dims[[3L]]$b)
})

test_that("re-saving a loaded spec produces a bit-exact identical file", {
  spec <- amm_spec(p = 2L,
                   dims = dimwise(a = ~ x1 + x2, b = ~ x1),
                   W = W_basis(type = "bspline", degree = 2,
                               knots = c(0.3, 0.7)),
                   x_vars = c("x1", "x2"))
  path1 <- tempfile(fileext = ".gdpar")
  path2 <- tempfile(fileext = ".gdpar")
  amm_save_spec(spec, path1)
  spec2 <- amm_load_spec(path1)
  amm_save_spec(spec2, path2)
  expect_identical(readLines(path1), readLines(path2))
})

test_that("amm_save_spec rejects a user-defined W_basis", {
  wb_user <- W_basis(type = "user",
                     basis_fn = function(theta) as.numeric(theta),
                     dim = 1L)
  spec <- amm_spec(W = wb_user)
  path <- tempfile(fileext = ".gdpar")
  expect_error(amm_save_spec(spec, path), class = "gdpar_input_error")
})

test_that("amm_load_spec rejects a file declaring W.type = user", {
  pkg_ver <- as.character(utils::packageVersion("gdpar"))
  path <- tempfile(fileext = ".gdpar")
  writeLines(c(
    paste0("# gdpar_spec_version: ", pkg_ver),
    "p: 1",
    "a: NULL",
    "b: NULL",
    "x_vars: NULL",
    "W.type: user"
  ), path)
  expect_error(amm_load_spec(path), class = "gdpar_input_error")
})

test_that("amm_save_spec rejects non-amm_spec input", {
  path <- tempfile(fileext = ".gdpar")
  expect_error(amm_save_spec(list(), path), class = "gdpar_input_error")
})

test_that("amm_save_spec rejects an empty path argument", {
  spec <- amm_spec()
  expect_error(amm_save_spec(spec, ""), class = "gdpar_input_error")
  expect_error(amm_save_spec(spec, NULL), class = "gdpar_input_error")
  expect_error(amm_save_spec(spec, c("a", "b")), class = "gdpar_input_error")
})

test_that("amm_load_spec aborts on a non-existent file", {
  expect_error(
    amm_load_spec(tempfile(fileext = ".gdpar")),
    class = "gdpar_input_error"
  )
})

test_that("amm_load_spec aborts when the version header is missing", {
  path <- tempfile(fileext = ".gdpar")
  writeLines(c("p: 1", "a: NULL", "b: NULL"), path)
  expect_error(amm_load_spec(path), class = "gdpar_input_error")
})

test_that("amm_load_spec aborts on a version mismatch", {
  path <- tempfile(fileext = ".gdpar")
  writeLines(c(
    "# gdpar_spec_version: 99.99.99",
    "p: 1",
    "a: NULL",
    "b: NULL"
  ), path)
  expect_error(amm_load_spec(path), class = "gdpar_input_error")
})

test_that("amm_load_spec aborts on a malformed record line (no colon)", {
  pkg_ver <- as.character(utils::packageVersion("gdpar"))
  path <- tempfile(fileext = ".gdpar")
  writeLines(c(
    paste0("# gdpar_spec_version: ", pkg_ver),
    "p 1",
    "a: NULL",
    "b: NULL"
  ), path)
  expect_error(amm_load_spec(path), class = "gdpar_input_error")
})

test_that("amm_load_spec aborts on an unknown key", {
  pkg_ver <- as.character(utils::packageVersion("gdpar"))
  path <- tempfile(fileext = ".gdpar")
  writeLines(c(
    paste0("# gdpar_spec_version: ", pkg_ver),
    "p: 1",
    "a: NULL",
    "b: NULL",
    "rogue.key: 42"
  ), path)
  expect_error(amm_load_spec(path), class = "gdpar_input_error")
})

test_that("amm_load_spec aborts on a duplicate key", {
  pkg_ver <- as.character(utils::packageVersion("gdpar"))
  path <- tempfile(fileext = ".gdpar")
  writeLines(c(
    paste0("# gdpar_spec_version: ", pkg_ver),
    "p: 1",
    "p: 2",
    "a: NULL",
    "b: NULL"
  ), path)
  expect_error(amm_load_spec(path), class = "gdpar_input_error")
})

test_that("amm_load_spec aborts when required key p is missing", {
  pkg_ver <- as.character(utils::packageVersion("gdpar"))
  path <- tempfile(fileext = ".gdpar")
  writeLines(c(
    paste0("# gdpar_spec_version: ", pkg_ver),
    "a: NULL",
    "b: NULL"
  ), path)
  expect_error(amm_load_spec(path), class = "gdpar_input_error")
})

test_that("amm_load_spec aborts on an invalid p value", {
  pkg_ver <- as.character(utils::packageVersion("gdpar"))
  path <- tempfile(fileext = ".gdpar")
  writeLines(c(
    paste0("# gdpar_spec_version: ", pkg_ver),
    "p: not_an_int",
    "a: NULL",
    "b: NULL"
  ), path)
  expect_error(amm_load_spec(path), class = "gdpar_input_error")
})

test_that("amm_load_spec aborts when a/b are non-NULL while p > 1", {
  pkg_ver <- as.character(utils::packageVersion("gdpar"))
  path <- tempfile(fileext = ".gdpar")
  writeLines(c(
    paste0("# gdpar_spec_version: ", pkg_ver),
    "p: 2",
    "a: ~ x1",
    "b: NULL",
    "dims.1.a: NULL",
    "dims.1.b: NULL",
    "dims.2.a: NULL",
    "dims.2.b: NULL"
  ), path)
  expect_error(amm_load_spec(path), class = "gdpar_input_error")
})

test_that("amm_load_spec aborts when dims.K.* records are missing for p > 1", {
  pkg_ver <- as.character(utils::packageVersion("gdpar"))
  path <- tempfile(fileext = ".gdpar")
  writeLines(c(
    paste0("# gdpar_spec_version: ", pkg_ver),
    "p: 3",
    "a: NULL",
    "b: NULL",
    "dims.1.a: ~ x1",
    "dims.1.b: NULL",
    "dims.2.a: NULL",
    "dims.2.b: NULL"
  ), path)
  expect_error(amm_load_spec(path), class = "gdpar_input_error")
})

test_that("amm_load_spec aborts when scalar p = 1 file contains dims records", {
  pkg_ver <- as.character(utils::packageVersion("gdpar"))
  path <- tempfile(fileext = ".gdpar")
  writeLines(c(
    paste0("# gdpar_spec_version: ", pkg_ver),
    "p: 1",
    "a: NULL",
    "b: NULL",
    "dims.1.a: ~ x1",
    "dims.1.b: NULL"
  ), path)
  expect_error(amm_load_spec(path), class = "gdpar_input_error")
})

test_that("amm_load_spec aborts on a malformed formula value", {
  pkg_ver <- as.character(utils::packageVersion("gdpar"))
  path <- tempfile(fileext = ".gdpar")
  writeLines(c(
    paste0("# gdpar_spec_version: ", pkg_ver),
    "p: 1",
    "a: not_a_formula_or_null",
    "b: NULL"
  ), path)
  expect_error(amm_load_spec(path), class = "gdpar_input_error")
})

test_that("amm_load_spec aborts on a two-sided formula value", {
  pkg_ver <- as.character(utils::packageVersion("gdpar"))
  path <- tempfile(fileext = ".gdpar")
  writeLines(c(
    paste0("# gdpar_spec_version: ", pkg_ver),
    "p: 1",
    "a: y ~ x1",
    "b: NULL"
  ), path)
  expect_error(amm_load_spec(path), class = "gdpar_input_error")
})

test_that("amm_load_spec aborts on a malformed x_vars value", {
  pkg_ver <- as.character(utils::packageVersion("gdpar"))
  path <- tempfile(fileext = ".gdpar")
  writeLines(c(
    paste0("# gdpar_spec_version: ", pkg_ver),
    "p: 1",
    "a: NULL",
    "b: NULL",
    "x_vars: [x1, x2]"
  ), path)
  expect_error(amm_load_spec(path), class = "gdpar_input_error")
})

test_that("amm_load_spec aborts on bspline without knots or df", {
  pkg_ver <- as.character(utils::packageVersion("gdpar"))
  path <- tempfile(fileext = ".gdpar")
  writeLines(c(
    paste0("# gdpar_spec_version: ", pkg_ver),
    "p: 1",
    "a: NULL",
    "b: NULL",
    "W.type: bspline",
    "W.degree: 3"
  ), path)
  expect_error(amm_load_spec(path), class = "gdpar_input_error")
})

test_that("amm_load_spec aborts when both knots and df are present", {
  pkg_ver <- as.character(utils::packageVersion("gdpar"))
  path <- tempfile(fileext = ".gdpar")
  writeLines(c(
    paste0("# gdpar_spec_version: ", pkg_ver),
    "p: 1",
    "a: NULL",
    "b: NULL",
    "W.type: bspline",
    "W.degree: 3",
    "W.knots: c(0.5)",
    "W.df: 4"
  ), path)
  expect_error(amm_load_spec(path), class = "gdpar_input_error")
})

test_that("amm_save_spec writes a stable header and body in canonical order", {
  spec <- amm_spec(a = ~ x1 + x2, b = ~ x1,
                   W = W_basis(type = "polynomial", degree = 2),
                   x_vars = c("x1", "x2"))
  path <- tempfile(fileext = ".gdpar")
  amm_save_spec(spec, path)
  lines <- readLines(path)
  pkg_ver <- as.character(utils::packageVersion("gdpar"))
  expect_equal(lines[[1L]], paste0("# gdpar_spec_version: ", pkg_ver))
  expect_equal(lines[[2L]], "p: 1")
  expect_equal(lines[[3L]], "a: ~x1 + x2")
  expect_equal(lines[[4L]], "b: ~x1")
  expect_equal(lines[[5L]], "x_vars: c(\"x1\", \"x2\")")
  expect_equal(lines[[6L]], "W.type: polynomial")
  expect_equal(lines[[7L]], "W.degree: 2")
})

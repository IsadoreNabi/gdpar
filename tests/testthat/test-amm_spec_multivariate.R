test_that("amm_spec defaults preserve p = 1L scalar path", {
  spec <- amm_spec()
  expect_equal(spec$p, 1L)
  expect_null(spec$dims)
  expect_null(spec$a)
  expect_null(spec$b)
  expect_null(spec$W)
})

test_that("amm_spec(p = 1L) backward compat: structure unchanged", {
  spec <- amm_spec(a = ~ x1 + x2, b = ~ x1,
                   W = W_basis("polynomial", 2))
  expect_equal(spec$p, 1L)
  expect_null(spec$dims)
  expect_equal(deparse(spec$a), "~x1 + x2")
  expect_equal(deparse(spec$b), "~x1")
  expect_s3_class(spec$W, "W_basis")
  expect_equal(spec$level, 2L)
})

test_that("amm_spec rejects dims when p = 1L", {
  expect_error(
    amm_spec(p = 1L, dims = dimwise(a = ~ x1)),
    class = "gdpar_input_error"
  )
})

test_that("amm_spec(p > 1L) requires dims", {
  expect_error(
    amm_spec(p = 2L),
    class = "gdpar_input_error"
  )
})

test_that("amm_spec(p > 1L) rejects scalar a or b arguments", {
  expect_error(
    amm_spec(p = 2L, a = ~ x1, dims = dimwise(b = ~ x1)),
    class = "gdpar_input_error"
  )
  expect_error(
    amm_spec(p = 2L, b = ~ x1, dims = dimwise(a = ~ x1)),
    class = "gdpar_input_error"
  )
})

test_that("amm_spec accepts dimwise() as dims for p > 1L", {
  spec <- amm_spec(
    p    = 2L,
    dims = dimwise(a = ~ x1 + x2, b = ~ x1),
    W    = W_basis("polynomial", 1)
  )
  expect_equal(spec$p, 2L)
  expect_null(spec$a)
  expect_null(spec$b)
  expect_length(spec$dims, 2L)
  for (k in seq_len(2L)) {
    expect_equal(deparse(spec$dims[[k]]$a), "~x1 + x2")
    expect_equal(deparse(spec$dims[[k]]$b), "~x1")
  }
  expect_equal(spec$level, 2L)
})

test_that("amm_spec applies override() composed into dimwise()", {
  spec <- amm_spec(
    p    = 3L,
    dims = override(dimwise(a = ~ x1, b = ~ x1), k = 2L, a = ~ x1 + x2)
  )
  expect_equal(deparse(spec$dims[[1L]]$a), "~x1")
  expect_equal(deparse(spec$dims[[2L]]$a), "~x1 + x2")
  expect_equal(deparse(spec$dims[[3L]]$a), "~x1")
})

test_that("amm_spec applies override(k, b = NULL) to disable per-dim", {
  spec <- amm_spec(
    p    = 3L,
    dims = override(dimwise(a = ~ x1, b = ~ x1), k = 3L, b = NULL)
  )
  expect_equal(deparse(spec$dims[[1L]]$b), "~x1")
  expect_equal(deparse(spec$dims[[2L]]$b), "~x1")
  expect_null(spec$dims[[3L]]$b)
})

test_that("amm_spec accepts plain list of length p for dims", {
  spec <- amm_spec(
    p    = 2L,
    dims = list(
      list(a = ~ x1 + x2, b = ~ x1),
      list(a = ~ x1,      b = ~ x2)
    )
  )
  expect_equal(spec$p, 2L)
  expect_equal(deparse(spec$dims[[1L]]$a), "~x1 + x2")
  expect_equal(deparse(spec$dims[[2L]]$b), "~x2")
})

test_that("amm_spec rejects plain list of wrong length", {
  expect_error(
    amm_spec(
      p    = 3L,
      dims = list(list(a = ~ x1), list(a = ~ x1))
    ),
    class = "gdpar_input_error"
  )
})

test_that("amm_spec rejects bare formula in dims when p > 1L", {
  expect_error(
    amm_spec(p = 2L, dims = ~ x1),
    class = "gdpar_input_error"
  )
})

test_that("amm_spec rejects list entries that are not lists", {
  expect_error(
    amm_spec(p = 2L, dims = list(~ x1, ~ x2)),
    class = "gdpar_input_error"
  )
})

test_that("amm_spec rejects two-sided formula inside plain dims list", {
  expect_error(
    amm_spec(
      p    = 2L,
      dims = list(list(a = y ~ x1), list(a = ~ x1))
    ),
    class = "gdpar_input_error"
  )
})

test_that("amm_spec rejects unknown dims class", {
  expect_error(
    amm_spec(p = 2L, dims = 42L),
    class = "gdpar_input_error"
  )
})

test_that("amm_spec records level correctly for p > 1 multi-dim cases", {
  spec0 <- amm_spec(p = 2L, dims = dimwise())
  expect_equal(spec0$level, 0L)
  spec1 <- amm_spec(p = 2L, dims = dimwise(a = ~ x1))
  expect_equal(spec1$level, 1L)
  spec2 <- amm_spec(p = 2L, dims = dimwise(a = ~ x1, b = ~ x1))
  expect_equal(spec2$level, 2L)
  spec_W <- amm_spec(p = 2L, dims = dimwise(),
                     W = W_basis("polynomial", 1))
  expect_equal(spec_W$level, 2L)
})

test_that("amm_spec level for partial per-dim overrides reflects union", {
  spec <- amm_spec(
    p    = 3L,
    dims = override(dimwise(), k = 2L, a = ~ x1)
  )
  expect_equal(spec$level, 1L)
})

test_that("amm_spec assigns class amm_spec", {
  spec <- amm_spec(p = 2L, dims = dimwise(a = ~ x1, b = ~ x1))
  expect_s3_class(spec, "amm_spec")
})

test_that("print method for amm_spec(p > 1L) prints per-k entries", {
  spec <- amm_spec(p = 2L, dims = dimwise(a = ~ x1, b = ~ x1))
  expect_output(print(spec), "p \\(dim theta_i\\)\\s+: 2")
  expect_output(print(spec), "dims \\(per-k a, b\\)")
  expect_output(print(spec), "k = 1")
  expect_output(print(spec), "k = 2")
})

test_that("print method for amm_spec(p = 1L) keeps scalar fields", {
  spec <- amm_spec(a = ~ x1)
  expect_output(print(spec), "p \\(dim theta_i\\)\\s+: 1")
  expect_output(print(spec), "a \\(additive\\)")
})

test_that("gdpar rejects amm_spec with p > 1L with informative error", {
  df <- make_gaussian_data(n = 50, seed = 1)
  spec <- amm_spec(p = 2L, dims = dimwise(a = ~ x1, b = ~ x1))
  expect_error(
    gdpar(y ~ x1, amm = spec, data = df,
          chains = 1L, iter_warmup = 50L, iter_sampling = 50L,
          verbose = FALSE),
    class = "gdpar_input_error"
  )
})

test_that("amm_spec rejects non-positive integer p", {
  expect_error(amm_spec(p = 0L), class = "gdpar_input_error")
  expect_error(amm_spec(p = -1L), class = "gdpar_input_error")
  expect_error(amm_spec(p = 1.5), class = "gdpar_input_error")
})

test_that("amm_spec(p > 1L) preserves x_vars argument", {
  spec <- amm_spec(
    p      = 2L,
    dims   = dimwise(a = ~ x1, b = ~ x1),
    W      = W_basis("polynomial", 1),
    x_vars = c("x1", "x2")
  )
  expect_equal(spec$x_vars, c("x1", "x2"))
})

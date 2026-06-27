test_that("polynomial W_basis at p = 1 keeps scalar behavior unchanged", {
  wb <- W_basis(type = "polynomial", degree = 2)
  expect_null(wb$p)
  expect_null(wb$block_indices)
  expect_equal(wb$evaluator(0.5), c(0.5, 0.25))
  expect_equal(wb$evaluator(2), c(2, 4))
})

test_that("polynomial W_basis at p > 1 uses blocked-by-coordinate ordering", {
  wb <- W_basis(type = "polynomial", degree = 2, p = 2L)
  out <- wb$evaluator(c(0.5, 0.3))
  expect_equal(out, c(0.5, 0.25, 0.3, 0.09))
  expect_equal(wb$dim, 4L)
  expect_equal(wb$p, 2L)
})

test_that("polynomial W_basis populates block_indices when p is supplied", {
  wb <- W_basis(type = "polynomial", degree = 2, p = 3L)
  expect_length(wb$block_indices, 3L)
  expect_equal(wb$block_indices[[1L]], 1:2)
  expect_equal(wb$block_indices[[2L]], 3:4)
  expect_equal(wb$block_indices[[3L]], 5:6)
})

test_that("polynomial W_basis with degree 3 and p 2 has correct dims", {
  wb <- W_basis(type = "polynomial", degree = 3, p = 2L)
  expect_equal(wb$dim, 6L)
  expect_equal(wb$block_indices[[1L]], 1:3)
  expect_equal(wb$block_indices[[2L]], 4:6)
})

test_that("materialize_W_basis works when p was not given at construction", {
  wb <- W_basis(type = "polynomial", degree = 2)
  expect_null(wb$p)
  wb_mat <- gdpar:::materialize_W_basis(wb, p = 2L)
  expect_equal(wb_mat$p, 2L)
  expect_equal(wb_mat$dim, 4L)
  expect_length(wb_mat$block_indices, 2L)
})

test_that("materialize_W_basis errors on p mismatch with construction-time p", {
  wb <- W_basis(type = "polynomial", degree = 2, p = 2L)
  expect_error(
    gdpar:::materialize_W_basis(wb, p = 3L),
    class = "gdpar_input_error"
  )
})

test_that("materialize_W_basis with same p is idempotent", {
  wb <- W_basis(type = "polynomial", degree = 2, p = 2L)
  wb2 <- gdpar:::materialize_W_basis(wb, p = 2L)
  expect_equal(wb2$dim, wb$dim)
  expect_equal(wb2$p, wb$p)
  expect_equal(wb2$block_indices, wb$block_indices)
})

test_that("bspline W_basis with df and p > 1 yields blocked output", {
  wb <- W_basis(type = "bspline", df = 3, degree = 2, p = 2L)
  expect_equal(wb$p, 2L)
  expect_length(wb$block_indices, 2L)
  expect_equal(length(wb$block_indices[[1L]]), 3L)
  expect_equal(length(wb$block_indices[[2L]]), 3L)
  expect_equal(wb$dim, 6L)
})

test_that("user W_basis at p = 2 validates dim against evaluator output", {
  wb <- W_basis(type = "user",
                basis_fn = function(theta) c(sum(theta), prod(theta), 1),
                dim = 3,
                p = 2L)
  expect_equal(wb$p, 2L)
  expect_equal(wb$dim, 3L)
  expect_null(wb$block_indices)
})

test_that("user W_basis dim mismatch is detected at construction with p", {
  expect_error(
    W_basis(type = "user",
            basis_fn = function(theta) c(sum(theta), prod(theta)),
            dim = 3,
            p = 2L),
    class = "gdpar_input_error"
  )
})

test_that("as_per_k splits polynomial into p univariate bases", {
  wb <- W_basis(type = "polynomial", degree = 2, p = 3L)
  subs <- as_per_k(wb)
  expect_length(subs, 3L)
  for (sub in subs) {
    expect_s3_class(sub, "W_basis")
    expect_equal(sub$type, "polynomial")
    expect_equal(sub$p, 1L)
    expect_equal(sub$dim, 2L)
    expect_equal(sub$degree, 2L)
  }
})

test_that("as_per_k splits bspline into p univariate bases", {
  wb <- W_basis(type = "bspline", df = 4, degree = 2, p = 2L)
  subs <- as_per_k(wb)
  expect_length(subs, 2L)
  expect_equal(subs[[1L]]$type, "bspline")
  expect_equal(subs[[1L]]$dim, 4L)
  expect_equal(subs[[1L]]$p, 1L)
})

test_that("as_per_k warns and returns NULL for user bases", {
  wb <- W_basis(type = "user",
                basis_fn = function(theta) c(sum(theta), 1),
                dim = 2,
                p = 2L)
  expect_warning(out <- as_per_k(wb))
  expect_null(out)
})

test_that("as_per_k errors when wb is not materialized", {
  wb <- W_basis(type = "polynomial", degree = 2)
  expect_error(as_per_k(wb), class = "gdpar_input_error")
})

test_that("as_per_k preserves degree across split sub-bases", {
  wb <- W_basis(type = "polynomial", degree = 3, p = 2L)
  subs <- as_per_k(wb)
  for (sub in subs) {
    expect_equal(sub$degree, 3L)
  }
})

test_that("print method shows p and block_indices when present", {
  wb <- W_basis(type = "polynomial", degree = 2, p = 2L)
  out <- paste(capture.output(print(wb)), collapse = "\n")
  expect_match(out, "p\\s+: 2")
  expect_match(out, "multivariate")
  expect_match(out, "block_indices")
})

test_that("print method for p = 1 shows univariate label", {
  wb <- W_basis(type = "polynomial", degree = 2, p = 1L)
  out <- paste(capture.output(print(wb)), collapse = "\n")
  expect_match(out, "univariate")
})

test_that("W_basis rejects non-positive p at construction", {
  expect_error(
    W_basis(type = "polynomial", degree = 2, p = 0L),
    class = "gdpar_input_error"
  )
  expect_error(
    W_basis(type = "polynomial", degree = 2, p = -1L),
    class = "gdpar_input_error"
  )
})

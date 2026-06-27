test_that("polynomial W_basis evaluates correctly", {
  wb <- W_basis(type = "polynomial", degree = 2)
  expect_s3_class(wb, "W_basis")
  expect_equal(wb$evaluator(0.5), c(0.5, 0.25))
  expect_equal(wb$evaluator(2), c(2, 4))
})

test_that("polynomial degree must be a positive integer", {
  expect_error(
    W_basis(type = "polynomial", degree = -1),
    class = "gdpar_input_error"
  )
  expect_error(
    W_basis(type = "polynomial", degree = 0),
    class = "gdpar_input_error"
  )
})

test_that("user W_basis enforces dim consistency", {
  wb <- W_basis(type = "user",
                basis_fn = function(theta) c(theta, theta^2, sin(theta)),
                dim = 3)
  expect_equal(wb$evaluator(0.5), c(0.5, 0.25, sin(0.5)))
})

test_that("user W_basis requires basis_fn", {
  expect_error(
    W_basis(type = "user", dim = 3),
    class = "gdpar_input_error"
  )
})

test_that("bspline W_basis requires either knots or df", {
  expect_error(
    W_basis(type = "bspline", degree = 3),
    class = "gdpar_input_error"
  )
  expect_error(
    W_basis(type = "bspline", degree = 3, knots = c(0.3), df = 5),
    class = "gdpar_input_error"
  )
})

test_that("materialize_W_basis sets dim from probe evaluation", {
  wb <- W_basis(type = "polynomial", degree = 3)
  wb_mat <- gdpar:::materialize_W_basis(wb, p = 1L)
  expect_equal(wb_mat$dim, 3L)
})

test_that("user W_basis dim mismatch is detected", {
  wb <- W_basis(type = "user",
                basis_fn = function(theta) c(theta, theta^2),
                dim = 5)
  expect_error(
    gdpar:::materialize_W_basis(wb, p = 1L),
    class = "gdpar_input_error"
  )
})

test_that("print method for W_basis runs without error", {
  expect_output(print(W_basis(type = "polynomial", degree = 2)),
                "W_basis")
})

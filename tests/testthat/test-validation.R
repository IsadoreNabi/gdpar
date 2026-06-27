test_that("assert_inherits accepts a matching class", {
  x <- structure(list(), class = "my_class")
  expect_silent(gdpar:::assert_inherits(x, "my_class", "x"))
})

test_that("assert_inherits rejects non-matching class", {
  expect_error(
    gdpar:::assert_inherits(1L, "my_class", "x"),
    class = "gdpar_input_error"
  )
})

test_that("assert_numeric_scalar enforces range", {
  expect_silent(gdpar:::assert_numeric_scalar(0.5, "x", 0, 1))
  expect_error(
    gdpar:::assert_numeric_scalar(2, "x", 0, 1),
    class = "gdpar_input_error"
  )
  expect_error(
    gdpar:::assert_numeric_scalar(c(0.5, 0.6), "x", 0, 1),
    class = "gdpar_input_error"
  )
  expect_error(
    gdpar:::assert_numeric_scalar(NULL, "x"),
    class = "gdpar_input_error"
  )
  expect_silent(
    gdpar:::assert_numeric_scalar(NULL, "x", allow_null = TRUE)
  )
})

test_that("assert_count enforces positive integer", {
  expect_silent(gdpar:::assert_count(3L, "x"))
  expect_silent(gdpar:::assert_count(3, "x"))
  expect_error(
    gdpar:::assert_count(0, "x"),
    class = "gdpar_input_error"
  )
  expect_error(
    gdpar:::assert_count(-1, "x"),
    class = "gdpar_input_error"
  )
  expect_error(
    gdpar:::assert_count(2.5, "x"),
    class = "gdpar_input_error"
  )
})

test_that("assert_one_sided_formula accepts ~ x and rejects y ~ x", {
  expect_silent(gdpar:::assert_one_sided_formula(~ x, "f"))
  expect_error(
    gdpar:::assert_one_sided_formula(y ~ x, "f"),
    class = "gdpar_input_error"
  )
  expect_silent(gdpar:::assert_one_sided_formula(NULL, "f"))
})

test_that("assert_data_frame checks required variables", {
  df <- data.frame(x = 1, y = 2)
  expect_silent(gdpar:::assert_data_frame(df, "df", c("x", "y")))
  expect_error(
    gdpar:::assert_data_frame(df, "df", c("z")),
    class = "gdpar_input_error"
  )
})

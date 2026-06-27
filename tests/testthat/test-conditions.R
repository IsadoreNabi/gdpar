test_that("gdpar_abort signals a structured error with custom class", {
  expect_error(
    gdpar:::gdpar_abort("test error", class = "my_test_error"),
    class = "my_test_error"
  )
  expect_error(
    gdpar:::gdpar_abort("test error", class = "my_test_error"),
    class = "gdpar_error"
  )
})

test_that("gdpar_warn signals a structured warning with custom class", {
  expect_warning(
    gdpar:::gdpar_warn("test warning", class = "my_test_warning"),
    class = "my_test_warning"
  )
  expect_warning(
    gdpar:::gdpar_warn("test warning", class = "my_test_warning"),
    class = "gdpar_warning"
  )
})

test_that("gdpar_inform signals a structured message with custom class", {
  expect_message(
    gdpar:::gdpar_inform("test message", class = "my_test_message"),
    class = "my_test_message"
  )
})

test_that("require_suggested aborts when package missing", {
  expect_error(
    gdpar:::require_suggested("package_that_does_not_exist_12345",
                              "perform a fictional task"),
    class = "gdpar_missing_dependency_error"
  )
})

test_that("require_suggested returns invisibly when package present", {
  expect_silent(
    gdpar:::require_suggested("stats", "use base statistics")
  )
})

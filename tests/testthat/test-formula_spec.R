test_that("gdpar_formula_set() builds with single named formula (K=1)", {
  fs <- gdpar_formula_set(mu = y ~ x1 + x2)
  expect_s3_class(fs, "gdpar_formula_set")
  expect_equal(fs$outcome, "y")
  expect_equal(names(fs), "mu")
  expect_length(fs, 1L)
  expect_true(inherits(fs[["mu"]], "formula"))
})

test_that("gdpar_formula_set() builds with K=2 mu + sigma", {
  fs <- gdpar_formula_set(mu = y ~ x1, sigma = ~ x2)
  expect_s3_class(fs, "gdpar_formula_set")
  expect_equal(fs$outcome, "y")
  expect_equal(names(fs), c("mu", "sigma"))
  expect_length(fs, 2L)
})

test_that("gdpar_formula_set() preserves declaration order for K>=3", {
  fs <- gdpar_formula_set(
    mu    = y ~ x1,
    sigma = ~ x2,
    nu    = ~ x3
  )
  expect_equal(names(fs), c("mu", "sigma", "nu"))
  expect_length(fs, 3L)
})

test_that("gdpar_formula_set() rejects empty call", {
  expect_error(
    gdpar_formula_set(),
    class = "gdpar_input_error"
  )
})

test_that("gdpar_formula_set() rejects unnamed args", {
  expect_error(
    gdpar_formula_set(y ~ x1, ~ x2),
    class = "gdpar_input_error"
  )
})

test_that("gdpar_formula_set() rejects partially unnamed args", {
  expect_error(
    gdpar_formula_set(mu = y ~ x1, ~ x2),
    class = "gdpar_input_error"
  )
})

test_that("gdpar_formula_set() rejects first slot 1-sided", {
  expect_error(
    gdpar_formula_set(mu = ~ x1),
    class = "gdpar_input_error"
  )
})

test_that("gdpar_formula_set() rejects subsequent slots 2-sided", {
  expect_error(
    gdpar_formula_set(mu = y ~ x1, sigma = y ~ x2),
    class = "gdpar_input_error"
  )
})

test_that("gdpar_formula_set() rejects intercept suppression via -1 in first slot", {
  expect_error(
    gdpar_formula_set(mu = y ~ x1 - 1),
    class = "gdpar_input_error"
  )
})

test_that("gdpar_formula_set() rejects intercept suppression via +0 in subsequent slot", {
  expect_error(
    gdpar_formula_set(mu = y ~ x1, sigma = ~ x2 + 0),
    class = "gdpar_input_error"
  )
})

test_that("gdpar_formula_set() error message references theta_ref anchor canonization", {
  err <- tryCatch(
    gdpar_formula_set(mu = y ~ x1 - 1),
    gdpar_input_error = function(e) e
  )
  expect_match(conditionMessage(err), "theta_ref")
  expect_match(conditionMessage(err), "anchor")
})

test_that("gdpar_formula_set() rejects duplicate slot names", {
  # R coerces duplicate-name argument syntax in `list(...)`; we provoke
  # the same condition by constructing the call via do.call to bypass
  # R's argument-name unique-ifier.
  args <- list(mu = y ~ x1, mu = y ~ x2)
  names(args) <- c("mu", "mu")
  expect_error(
    do.call(gdpar_formula_set, args),
    class = "gdpar_input_error"
  )
})

test_that("gdpar_formula_set() rejects non-formula slot", {
  expect_error(
    gdpar_formula_set(mu = "y ~ x"),
    class = "gdpar_input_error"
  )
})

test_that("gdpar_formula_set() rejects non-symbol LHS in first slot", {
  expect_error(
    gdpar_formula_set(mu = log(y) ~ x1),
    class = "gdpar_input_error"
  )
})

test_that("gdpar_formula_set [[ accessor returns formula", {
  fs <- gdpar_formula_set(mu = y ~ x1, sigma = ~ x2)
  expect_true(inherits(fs[["mu"]], "formula"))
  expect_true(inherits(fs[["sigma"]], "formula"))
})

test_that("gdpar_formula_set [ accessor returns named list of formulas", {
  fs <- gdpar_formula_set(mu = y ~ x1, sigma = ~ x2)
  sub <- fs["mu"]
  expect_type(sub, "list")
  expect_length(sub, 1L)
  expect_equal(names(sub), "mu")
  expect_true(inherits(sub[[1L]], "formula"))
})

test_that("gdpar_formula_set print works and shows K + slots", {
  fs <- gdpar_formula_set(mu = y ~ x1, sigma = ~ x2)
  expect_output(print(fs), "gdpar_formula_set")
  expect_output(print(fs), "K = 2")
  expect_output(print(fs), "outcome")
  expect_output(print(fs), "y")
  expect_output(print(fs), "mu")
  expect_output(print(fs), "sigma")
})

test_that("gdpar_formula_set stores environment of the first formula", {
  e <- new.env()
  e$marker <- "test_env_marker"
  f <- y ~ x1
  environment(f) <- e
  fs <- gdpar_formula_set(mu = f)
  expect_identical(fs$env, e)
})

# ---- gdpar_bf() ----

test_that("gdpar_bf() with K=1 defaults first slot to 'mu'", {
  fs <- gdpar_bf(y ~ x1)
  expect_s3_class(fs, "gdpar_formula_set")
  expect_equal(names(fs), "mu")
  expect_equal(fs$outcome, "y")
})

test_that("gdpar_bf() with K=2 derives slot names from LHS of subsequent formulas", {
  fs <- gdpar_bf(y ~ x1, sigma ~ x2)
  expect_equal(names(fs), c("mu", "sigma"))
  expect_equal(fs$outcome, "y")
})

test_that("gdpar_bf() with K=3 preserves order", {
  fs <- gdpar_bf(y ~ x1, sigma ~ x2, nu ~ x3)
  expect_equal(names(fs), c("mu", "sigma", "nu"))
})

test_that("gdpar_bf() preserves explicit first-slot name when provided", {
  fs <- gdpar_bf(theta = y ~ x1)
  expect_equal(names(fs), "theta")
  expect_equal(fs$outcome, "y")
})

test_that("gdpar_bf() rejects empty call", {
  expect_error(
    gdpar_bf(),
    class = "gdpar_input_error"
  )
})

test_that("gdpar_bf() rejects non-formula args", {
  expect_error(
    gdpar_bf("y ~ x"),
    class = "gdpar_input_error"
  )
})

test_that("gdpar_bf() rejects first formula not two-sided", {
  expect_error(
    gdpar_bf(~ x1),
    class = "gdpar_input_error"
  )
})

test_that("gdpar_bf() rejects subsequent formulas one-sided (must carry LHS param name)", {
  expect_error(
    gdpar_bf(y ~ x1, ~ x2),
    class = "gdpar_input_error"
  )
})

test_that("gdpar_bf() rejects subsequent formulas with non-symbol LHS", {
  expect_error(
    gdpar_bf(y ~ x1, log(sigma) ~ x2),
    class = "gdpar_input_error"
  )
})

test_that("gdpar_bf() rejects intercept suppression in any formula", {
  expect_error(
    gdpar_bf(y ~ x1 - 1),
    class = "gdpar_input_error"
  )
  expect_error(
    gdpar_bf(y ~ x1, sigma ~ x2 + 0),
    class = "gdpar_input_error"
  )
})

test_that("gdpar_bf() rejects duplicate parameter names across formulas", {
  expect_error(
    gdpar_bf(y ~ x1, sigma ~ x2, sigma ~ x3),
    class = "gdpar_input_error"
  )
})

test_that("gdpar_bf() rejects explicit arg name not matching LHS", {
  expect_error(
    gdpar_bf(y ~ x1, log_sigma = sigma ~ x2),
    class = "gdpar_input_error"
  )
})

test_that("gdpar_bf() produces equivalent object to gdpar_formula_set() canonical form", {
  fs_bf  <- gdpar_bf(y ~ x1, sigma ~ x2)
  fs_can <- gdpar_formula_set(mu = y ~ x1, sigma = ~ x2)
  expect_equal(names(fs_bf), names(fs_can))
  expect_equal(fs_bf$outcome, fs_can$outcome)
  expect_equal(length(fs_bf), length(fs_can))
})

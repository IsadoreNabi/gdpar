test_that("gdpar_prior() returns class gdpar_prior with defaults", {
  pr <- gdpar_prior()
  expect_s3_class(pr, "gdpar_prior")
  expect_true(is.character(pr$theta_ref))
})

test_that("gdpar_prior rejects empty strings", {
  expect_error(
    gdpar_prior(theta_ref = ""),
    class = "gdpar_input_error"
  )
})

test_that("gdpar_prior rejects non-character", {
  expect_error(
    gdpar_prior(theta_ref = 1),
    class = "gdpar_input_error"
  )
})

test_that("print method for gdpar_prior runs without error", {
  expect_output(print(gdpar_prior()), "gdpar_prior")
})

test_that("priors_by_kind defaults to empty list and is preserved", {
  pr <- gdpar_prior()
  expect_true(!is.null(pr$priors_by_kind))
  expect_length(pr$priors_by_kind, 0L)
})

test_that("priors_by_kind accepts valid named overrides", {
  pr <- gdpar_prior(priors_by_kind = list(
    log_sigma = "normal(0, 0.5)",
    log_phi   = "student_t(3, 0, 0.5)"
  ))
  expect_equal(pr$priors_by_kind$log_sigma, "normal(0, 0.5)")
  expect_equal(pr$priors_by_kind$log_phi, "student_t(3, 0, 0.5)")
})

test_that("priors_by_kind rejects unrecognized kinds", {
  expect_error(
    gdpar_prior(priors_by_kind = list(foo_bar = "normal(0, 1)")),
    class = "gdpar_input_error"
  )
})

test_that("priors_by_kind rejects unnamed lists", {
  expect_error(
    gdpar_prior(priors_by_kind = list("normal(0, 1)")),
    class = "gdpar_input_error"
  )
})

test_that("priors_by_kind rejects empty strings", {
  expect_error(
    gdpar_prior(priors_by_kind = list(mu = "")),
    class = "gdpar_input_error"
  )
})

test_that(".gdpar_prior_for_kind returns override when present, canonical otherwise", {
  pr <- gdpar_prior(priors_by_kind = list(log_sigma = "normal(0, 0.5)"))
  expect_equal(gdpar:::.gdpar_prior_for_kind(pr, "log_sigma"),
               "normal(0, 0.5)")
  expect_equal(gdpar:::.gdpar_prior_for_kind(pr, "mu"),
               "normal(0, 2.5)")
})

test_that(".gdpar_canonical_prior_for_kind registers all known kinds", {
  kinds <- gdpar:::.gdpar_known_canonical_kinds()
  for (k in kinds) {
    p <- gdpar:::.gdpar_canonical_prior_for_kind(k)
    expect_true(is.character(p) && nchar(p) > 0L)
  }
})

test_that("print method shows priors_by_kind block when populated", {
  pr <- gdpar_prior(priors_by_kind = list(log_sigma = "normal(0, 0.5)"))
  out <- capture.output(print(pr))
  expect_true(any(grepl("priors_by_kind", out)))
  expect_true(any(grepl("log_sigma", out)))
})

test_that("backward compat: legacy args still accepted unchanged", {
  pr <- gdpar_prior(
    theta_ref       = "normal(0, 1)",
    sigma_theta_ref = "student_t(3, 0, 0.5)",
    sigma_a         = "student_t(3, 0, 0.5)",
    sigma_b         = "student_t(3, 0, 0.5)",
    sigma_W         = "student_t(3, 0, 0.5)",
    sigma_y         = "student_t(3, 0, 1)",
    phi             = "gamma(1, 0.1)"
  )
  expect_equal(pr$theta_ref, "normal(0, 1)")
  expect_equal(pr$phi, "gamma(1, 0.1)")
  expect_length(pr$priors_by_kind, 0L)
})

# Tests for gdpar_contraction_diagnostic (Block 4 Theorem 4B).
#
# Coverage:
#   (i)  argument-validation guards run unconditionally (no Stan).
#   (ii) end-to-end smoke test that exercises the full refit-curve pipeline.
#        Gated only by `skip_if_not_installed("cmdstanr")` (and `posterior`);
#        runs by default whenever the Stan toolchain is available
#        (D16 ampliado, Sesion 17 Etapa 2 de Sesion 8.4).

# ---- (i) Argument-validation guards ------------------------------------

test_that("gdpar_contraction_diagnostic rejects non gdpar_fit `fit`", {
  expect_error(
    gdpar_contraction_diagnostic(
      fit = list(), data = data.frame(x = 1, y = 2),
      replicates = 1L, verbose = FALSE
    ),
    class = "gdpar_input_error"
  )
})

test_that("gdpar_contraction_diagnostic rejects non data.frame `data`", {
  fit_mock <- structure(list(), class = "gdpar_fit")
  expect_error(
    gdpar_contraction_diagnostic(
      fit = fit_mock, data = "not a data frame",
      replicates = 1L, verbose = FALSE
    ),
    class = "gdpar_input_error"
  )
})

test_that("gdpar_contraction_diagnostic rejects level outside (0, 1)", {
  fit_mock <- structure(list(), class = "gdpar_fit")
  expect_error(
    gdpar_contraction_diagnostic(
      fit = fit_mock, data = data.frame(x = 1, y = 2),
      replicates = 1L, level = 1.5, verbose = FALSE
    ),
    class = "gdpar_input_error"
  )
})

# ---- (ii) End-to-end smoke (gated only by cmdstanr availability) -------

test_that("gdpar_contraction_diagnostic end-to-end on Gaussian (smoke)", {
  testthat::skip_if_not_installed("cmdstanr")
  testthat::skip_if_not_installed("posterior")
  df <- make_gaussian_data(n = 400L, seed = 19L)
  fit <- gdpar(
    formula       = y ~ x1,
    family        = gdpar_family("gaussian"),
    amm           = amm_spec(a = ~ x1),
    data          = df,
    chains        = 1L,
    iter_warmup   = 200L,
    iter_sampling = 200L,
    refresh       = 0L,
    verbose       = FALSE
  )
  report <- gdpar_contraction_diagnostic(
    fit, data = df, replicates = 1L,
    iter_warmup = 200L, iter_sampling = 200L, chains = 1L,
    verbose = FALSE
  )
  expect_s3_class(report, "gdpar_contraction_report")
  # D48 (Sesion 17 Etapa 2 de Sesion 8.4): the canonical print method emits
  # a structured report (header + table + slope estimate + verdict) by
  # design; the original expect_silent assumption was semantically
  # incorrect (it only passed because the smoke was opt-in and never
  # executed). We validate (i) that the print returns invisibly the
  # report object and (ii) that the output carries the canonical header.
  expect_output(print(report), "<gdpar_contraction_report>", fixed = TRUE)
  expect_invisible(print(report))
  # D49 (Sesion 18 Etapa 2 de Sesion 8.4): the report carries a structured
  # `$warnings` slot that records per-refit fallback notifications. Empty
  # in the happy path; populated when any refit failed (each entry is the
  # human-readable message that was simultaneously emitted via gdpar_warn).
  expect_true("warnings" %in% names(report))
  expect_type(report$warnings, "character")
})

# ---- (iii) D49 (Sesion 18 Etapa 2): print method surfaces $warnings

test_that("print.gdpar_contraction_report displays $warnings when populated", {
  # Mock report structure (no Stan invocation). Verifies that when the
  # diagnostic populates $warnings with refit-failure messages, the
  # canonical print method surfaces them under a `Warnings:` section,
  # matching the D48 norm (notification at runtime + persistent record
  # in object + display in print method).
  mock_report <- structure(
    list(
      table = data.frame(n = c(100L, 200L, 400L),
                          replicate = c(1L, 1L, 1L),
                          median_width = c(0.5, NA_real_, 0.2)),
      slope_estimate = -0.5,
      slope_se = 0.05,
      slope_ci_lower = -0.6,
      slope_ci_upper = -0.4,
      verdict = "Consistent with parametric n^{-1/2} rate.",
      level = 0.95,
      warnings = c("Refit at n = 200 (replicate 1) failed: timeout.")
    ),
    class = c("gdpar_contraction_report", "list")
  )
  expect_output(print(mock_report), "Warnings:", fixed = TRUE)
  expect_output(print(mock_report), "Refit at n = 200", fixed = TRUE)
  expect_invisible(print(mock_report))
})

test_that("print.gdpar_contraction_report omits Warnings: block when empty", {
  mock_report <- structure(
    list(
      table = data.frame(n = c(100L, 200L, 400L),
                          replicate = c(1L, 1L, 1L),
                          median_width = c(0.5, 0.35, 0.2)),
      slope_estimate = -0.5,
      slope_se = 0.05,
      slope_ci_lower = -0.6,
      slope_ci_upper = -0.4,
      verdict = "Consistent with parametric n^{-1/2} rate.",
      level = 0.95,
      warnings = character(0L)
    ),
    class = c("gdpar_contraction_report", "list")
  )
  out <- capture.output(print(mock_report))
  expect_false(any(grepl("^Warnings:", out)))
})

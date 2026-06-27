# Tests for gdpar_bvm_check (Block 4 Theorem 4C; refactored to MLE+Hessian
# in Phase E of Audit 8.4, D11 closed).
#
# Layered coverage:
#   (i)   unit tests of the codegen helper `.strip_priors_block` and of
#         `generate_stan_code(mle = TRUE)`. Run unconditionally.
#   (ii)  guard test that hierarchical fits (use_groups == 1) abort with
#         `gdpar_unsupported_feature_error`. Built with a minimal mock
#         object; does not invoke cmdstan and runs unconditionally.
#   (iii) end-to-end smoke test that fits a small Path 1 Gaussian model and
#         invokes the full bvm_check pipeline. Gated only by
#         `skip_if_not_installed("cmdstanr")` (and `posterior`); runs by
#         default whenever the Stan toolchain is available (D16 ampliado,
#         Sesion 17 Etapa 2 de Sesion 8.4).

# ---- (i) Codegen helpers: strip-priors and generate_stan_code(mle = TRUE) ----

test_that(".strip_priors_block removes the marked block inclusive of markers", {
  src <- paste(
    "model {",
    "  // BEGIN PRIORS",
    "  theta_ref[1] ~ normal(0, 1);",
    "  sigma_y[1] ~ exponential(1);",
    "  // END PRIORS",
    "  y_real ~ normal(eta, sigma_y[1]);",
    "}",
    sep = "\n"
  )
  out <- gdpar:::.strip_priors_block(src)
  expect_false(grepl("BEGIN PRIORS", out, fixed = TRUE))
  expect_false(grepl("END PRIORS",   out, fixed = TRUE))
  expect_false(grepl("theta_ref[1] ~", out, fixed = TRUE))
  expect_false(grepl("sigma_y[1] ~ exponential", out, fixed = TRUE))
  expect_true(grepl("y_real ~ normal(eta, sigma_y[1])", out, fixed = TRUE))
})

test_that(".strip_priors_block aborts when markers are absent or malformed", {
  expect_error(
    gdpar:::.strip_priors_block("model { y ~ normal(0, 1); }"),
    class = "gdpar_internal_error"
  )
  expect_error(
    gdpar:::.strip_priors_block(
      paste("// END PRIORS", "// BEGIN PRIORS", sep = "\n")
    ),
    class = "gdpar_internal_error"
  )
})

test_that("generate_stan_code(mle = TRUE) yields a prior-stripped CP variant", {
  prior <- gdpar_prior()
  src <- gdpar:::generate_stan_code(prior, mle = TRUE)
  expect_false(grepl("BEGIN PRIORS", src, fixed = TRUE))
  expect_false(grepl("END PRIORS",   src, fixed = TRUE))
  expect_false(grepl("theta_ref[1] ~", src, fixed = TRUE))
  expect_false(grepl("sigma_y[1] ~",   src, fixed = TRUE))
  expect_false(grepl("a_raw ~",        src, fixed = TRUE))
  expect_false(grepl("c_b_raw ~",      src, fixed = TRUE))
  expect_false(grepl("to_vector(W_raw) ~", src, fixed = TRUE))
  expect_true(grepl("y_real ~ normal(eta, sigma_y[1])", src, fixed = TRUE))
  # MLE mode forces CP for a and W: the assignments `a_coef[j] = a_raw[j];`
  # and `W_diff_x[jj] += basis_diff[k] * W_raw[k, jj];` (no `* sigma_a[1]`,
  # no `* sigma_W[1]`) must appear in the transformed-parameters block.
  # The template carries a documentation comment that lists the NCP form
  # verbatim, so we anchor the assertion on the assignment string itself.
  expect_true(grepl("a_coef[j] = a_raw[j];", src, fixed = TRUE))
  expect_true(grepl("basis_diff[k] * W_raw[k, jj];", src, fixed = TRUE))
})

# ---- (ii) Guard: bvm_check aborts on hierarchical fits (use_groups == 1) ----

test_that("gdpar_bvm_check aborts on hierarchical (use_groups == 1) fit", {
  fit_mock <- structure(
    list(stan_data = list(use_groups = 1L)),
    class = "gdpar_fit"
  )
  expect_error(
    gdpar_bvm_check(fit_mock, verbose = FALSE),
    class = "gdpar_unsupported_feature_error"
  )
})

# ---- (iii) End-to-end smoke (gated only by cmdstanr availability) -------

test_that("gdpar_bvm_check end-to-end on a Gaussian Path 1 fit (smoke)", {
  testthat::skip_if_not_installed("cmdstanr")
  testthat::skip_if_not_installed("posterior")
  df <- make_gaussian_data(n = 200L, seed = 17L)
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
  report <- gdpar_bvm_check(fit, verbose = FALSE)
  expect_s3_class(report, "gdpar_bvm_report")
  expect_true(is.data.frame(report$table))
  expect_true(all(c("variable", "bayes_mean", "bayes_lower", "bayes_upper",
                    "asymp_mean", "asymp_lower", "asymp_upper",
                    "width_ratio") %in% names(report$table)))
  expect_true("theta_ref[1]" %in% report$table$variable)
  expect_equal(report$level, 0.95)
  # D48 (Sesion 17 Etapa 2 de Sesion 8.4): the canonical print method emits
  # a structured report (header + table + any fallback warnings) by design;
  # the original expect_silent assumption was semantically incorrect (it
  # only passed because the smoke was opt-in and never executed). We
  # validate (i) that the print returns invisibly the report object and
  # (ii) that the output carries the canonical header.
  expect_output(print(report), "<gdpar_bvm_report>", fixed = TRUE)
  expect_invisible(print(report))
})

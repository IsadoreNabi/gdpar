# Tests for gdpar_eb (Sub-phase 8.6.B; Empirical-Bayes orchestrator).
#
# Three layers:
#   (i)   unit tests of the input-scope guards of .gdpar_eb_validate_inputs()
#         and .gdpar_eb_resolve_laplace_control(). No Stan compilation, no
#         posterior draws.
#   (ii)  unit test of the EB-rendering of amm_eb_marginal.stan and
#         amm_eb_conditional.stan from a default gdpar_prior() (placeholder
#         substitution + structural invariants). No compilation.
#   (iii) opt-in smoke test gated by env var GDPAR_RUN_STAN_SMOKE_EB and
#         skip_if_not_installed("cmdstanr"), consistent with the per-feature
#         opt-in pattern of the suite. The expanded smoke-by-family matrix
#         (17 cases of Charter §3.2 Sub-phase 8.6.B), the EB-vs-FB
#         first-order coherence test (Theorem 7A scalar), and the 17
#         golden snapshots are listed as TODO in
#         HANDOFF_SUBFASE_8_6_B_PARCIAL.md §5 and are deferred to the
#         next clean session.

# ---- (i) Input scope guards (unconditional) -----------------------------

test_that("gdpar_eb rejects formula that is not a two-sided formula", {
  expect_error(
    gdpar_eb(formula = ~ x1, data = data.frame(x1 = 1)),
    class = "gdpar_input_error"
  )
})

test_that("gdpar_eb rejects data that is not a data.frame", {
  expect_error(
    gdpar_eb(formula = y ~ x1, data = list(y = 1, x1 = 1)),
    class = "gdpar_input_error"
  )
})

test_that("gdpar_eb under D34 (8.6.C) admits amm$p > 1 (Path A; no scope guard)", {
  # Sub-phase 8.6.C decision D34 (canonized Session 10): Path A admits
  # K = 1, p > 1. The validate-inputs scope guard against p > 1 is
  # removed; failures downstream (e.g. inconsistent family class) raise
  # gdpar_input_error rather than gdpar_unsupported_feature_error.
  amm_p2 <- amm_spec(a = ~ x1)
  amm_p2$p <- 2L
  err <- tryCatch(
    gdpar_eb(
      formula = y ~ x1,
      amm     = amm_p2,
      data    = data.frame(y = rnorm(10), x1 = rnorm(10))
    ),
    error = identity
  )
  expect_true(inherits(err, "error"))
  expect_false(inherits(err, "gdpar_unsupported_feature_error"))
})

test_that("gdpar_eb rejects eb_correction that is not a scalar logical", {
  expect_error(
    gdpar_eb(
      formula = y ~ x1,
      data = data.frame(y = rnorm(5), x1 = rnorm(5)),
      eb_correction = "yes"
    ),
    class = "gdpar_input_error"
  )
})

test_that("gdpar_eb rejects laplace_control with un-named entries", {
  expect_error(
    gdpar_eb(
      formula = y ~ x1,
      data = data.frame(y = rnorm(5), x1 = rnorm(5)),
      laplace_control = list(5L)
    ),
    class = "gdpar_input_error"
  )
})

test_that("gdpar_eb K=1 path rejects stan_id outside {1,2,3,4} (D34 check_stan_id_for_path)", {
  # Sub-phase 8.6.C (D34): for K = 1 the EB template body inherits
  # amm_main.stan or amm_distrib_multi.stan, both of which only branch
  # on family_id in {1, 2, 3, 4}. Canonical K >= 2 families (Beta,
  # Gamma, Student-t, Tweedie, mixtures, hurdle) must enter via the
  # Path B input (gdpar_formula_set or named list of amm_spec).
  fake_family <- list(
    name    = "beta",
    link    = "logit",
    K       = 1L,
    stan_id = 5L
  )
  class(fake_family) <- c("gdpar_family", "list")
  expect_error(
    gdpar_eb(
      formula = y ~ x1,
      family  = fake_family,
      data    = data.frame(y = runif(5, 0.1, 0.9), x1 = rnorm(5))
    ),
    class = "gdpar_unsupported_feature_error"
  )
})

# ---- (ii) Stan template rendering (unconditional) -----------------------

test_that(".gdpar_eb_generate_stan_marginal substitutes every placeholder", {
  prior <- gdpar_prior()
  src <- gdpar:::.gdpar_eb_generate_stan_marginal(prior)
  expect_type(src, "character")
  expect_length(src, 1L)
  expect_false(grepl("{{", src, fixed = TRUE))
  expect_true(grepl("vector\\[J_groups\\] theta_ref", src))
  expect_false(grepl("vector\\[J_groups\\] theta_ref_data", src))
  expect_true(grepl("K_slots", src, fixed = TRUE))
  expect_true(grepl("p_dim", src, fixed = TRUE))
})

test_that(".gdpar_eb_generate_stan_conditional substitutes every placeholder", {
  prior <- gdpar_prior()
  src <- gdpar:::.gdpar_eb_generate_stan_conditional(prior)
  expect_type(src, "character")
  expect_length(src, 1L)
  expect_false(grepl("{{", src, fixed = TRUE))
  expect_true(grepl("vector\\[J_groups\\] theta_ref_data", src))
  expect_false(grepl("vector\\[J_groups\\] theta_ref;", src))
  # mu_theta_ref / sigma_theta_ref appear only in human-readable
  # comments that document why the conditional template drops them; the
  # actual parameter declarations and sampling statements must be
  # absent. Strip line-comments before asserting.
  code_only <- paste(
    gsub("//.*$", "", strsplit(src, "\n", fixed = TRUE)[[1L]]),
    collapse = "\n"
  )
  expect_false(grepl("mu_theta_ref", code_only, fixed = TRUE))
  expect_false(grepl("sigma_theta_ref", code_only, fixed = TRUE))
})

test_that(".gdpar_eb_resolve_laplace_control merges defaults correctly", {
  lc <- gdpar:::.gdpar_eb_resolve_laplace_control(list())
  expect_equal(lc$multi_start_M, 5L)
  expect_equal(lc$kappa_threshold, 1e10)
  expect_equal(lc$ridge_init, 1e-6)
  expect_equal(lc$laplace_draws, 1000L)

  lc2 <- gdpar:::.gdpar_eb_resolve_laplace_control(
    list(multi_start_M = 3L, kappa_threshold = 1e8)
  )
  expect_equal(lc2$multi_start_M, 3L)
  expect_equal(lc2$kappa_threshold, 1e8)
  expect_equal(lc2$ridge_init, 1e-6)
})

test_that(".gdpar_eb_resolve_laplace_control warns on unknown entries", {
  expect_warning(
    gdpar:::.gdpar_eb_resolve_laplace_control(list(unknown_arg = 1)),
    class = "gdpar_diagnostic_warning"
  )
})

# ---- (iii) Opt-in smoke test (Gaussian, single canonical setup) ----------

test_that("gdpar_eb runs end-to-end on a Gaussian model (opt-in)", {
  skip_if_not_installed("cmdstanr")
  if (!nzchar(Sys.getenv("GDPAR_RUN_STAN_SMOKE_EB"))) {
    skip("Set GDPAR_RUN_STAN_SMOKE_EB to a non-empty value to enable.")
  }
  set.seed(20260524L)
  n <- 80L
  df <- data.frame(x1 = rnorm(n))
  df$y <- 0.5 + 0.3 * df$x1 + rnorm(n, sd = 0.4)
  fit_eb <- gdpar_eb(
    formula       = y ~ x1,
    family        = gdpar_family("gaussian"),
    amm           = amm_spec(a = ~ x1),
    data          = df,
    iter_warmup   = 200L,
    iter_sampling = 200L,
    chains        = 2L,
    refresh       = 0L,
    seed          = 12345L,
    laplace_control = list(multi_start_M = 2L, laplace_draws = 200L),
    verbose       = FALSE
  )
  expect_s3_class(fit_eb, "gdpar_eb_fit")
  expect_length(fit_eb$theta_ref_hat, 1L)
  expect_true(isTRUE(is.finite(fit_eb$theta_ref_hat)))
  expect_true(isTRUE(is.finite(fit_eb$theta_ref_se)))
  expect_true(isTRUE(is.finite(fit_eb$diagnostics_numerical$kappa)))
  expect_identical(fit_eb$path, "eb")
  expect_true(fit_eb$correction_applied)
})

test_that("print + summary + coef methods on gdpar_eb_fit work (opt-in)", {
  skip_if_not_installed("cmdstanr")
  if (!nzchar(Sys.getenv("GDPAR_RUN_STAN_SMOKE_EB"))) {
    skip("Set GDPAR_RUN_STAN_SMOKE_EB to a non-empty value to enable.")
  }
  set.seed(20260524L)
  n <- 60L
  df <- data.frame(x1 = rnorm(n))
  df$y <- 1 + 0.3 * df$x1 + rnorm(n, sd = 0.4)
  fit_eb <- gdpar_eb(
    formula       = y ~ x1,
    amm           = amm_spec(a = ~ x1),
    data          = df,
    iter_warmup   = 150L, iter_sampling = 150L, chains = 2L,
    refresh       = 0L, seed = 12345L,
    laplace_control = list(multi_start_M = 2L, laplace_draws = 200L),
    verbose       = FALSE
  )
  expect_invisible(print(fit_eb))
  sm <- summary(fit_eb)
  expect_s3_class(sm, "summary.gdpar_eb_fit")
  expect_true("estimate" %in% colnames(sm$theta_table))
  cf <- coef(fit_eb)
  expect_s3_class(cf, "gdpar_coef_eb")
  expect_identical(cf$theta_ref$method, "EB")
})

# Tests for the multivariate preflight (Phase H.1): resolver
# behavior under the no-preflight branches, all explicit error
# paths, and decision_to_logical helper. End-to-end Stan fits are
# env-gated under GDPAR_F6_FIT=1 (heavy compile + sample).

# ---- decision_to_logical ----------------------------------------

test_that("decision_to_logical maps CP/NCP/absent/per_k correctly", {
  expect_identical(
    gdpar:::decision_to_logical(c("CP", "NCP", "absent", "per_k")),
    c(TRUE, FALSE, FALSE, NA)
  )
  expect_identical(gdpar:::decision_to_logical(character(0)),
                   logical(0))
})

test_that("decision_to_logical aborts on unknown values", {
  expect_error(gdpar:::decision_to_logical("MAYBE"),
               class = "gdpar_internal_error")
})


# ---- resolver: explicit branches skip preflight -----------------

# Note: under the no-preflight branches we pass prior = NULL,
# stan_data = NULL, amm = list(p = 2L). The resolver does not touch
# them because no preflight is invoked.

mk_dummy_amm <- function(p = 2L) {
  structure(list(p = as.integer(p)), class = "amm_spec")
}

test_that("resolver: parametrization = 'ncp' skips preflight (multi)", {
  res <- gdpar:::resolve_parametrization_multi(
    parametrization             = "ncp",
    parametrization_a           = NULL,
    parametrization_W           = NULL,
    parametrization_aggregation = NULL,
    prior = NULL, stan_data = NULL, amm = mk_dummy_amm(2L),
    preflight_seed = NULL, verbose = FALSE
  )
  expect_false(res$cp_a)
  expect_false(res$cp_W)
  expect_identical(res$cp_a_per_k, c(FALSE, FALSE))
  expect_identical(res$cp_W_per_k, c(FALSE, FALSE))
  expect_null(res$report)
})

test_that("resolver: parametrization = 'cp' skips preflight (multi)", {
  # Even though cp on a multi is normally rejected, the resolver's
  # only abort path is on the per-component arg parametrization_a;
  # parametrization = 'cp' alone is not the same: the H.1 contract
  # is that 'cp' translates to (cp_a = TRUE, cp_W = TRUE) via the
  # globals, BUT cp_a = TRUE for multi is still unwired (H.2). To
  # keep tests behavior-locked, we exercise only the cp_W branch
  # via parametrization_W = 'cp'; for the full global 'cp' we
  # accept the resolver returns TRUE and let downstream codegen
  # decide.
  res <- gdpar:::resolve_parametrization_multi(
    parametrization             = "auto",
    parametrization_a           = "ncp",
    parametrization_W           = "cp",
    parametrization_aggregation = NULL,
    prior = NULL, stan_data = NULL, amm = mk_dummy_amm(3L),
    preflight_seed = NULL, verbose = FALSE
  )
  expect_false(res$cp_a)
  expect_true(res$cp_W)
  expect_identical(res$cp_a_per_k, rep(FALSE, 3L))
  expect_identical(res$cp_W_per_k, rep(TRUE,  3L))
  expect_null(res$report)
})

test_that("resolver honors cp_a = TRUE from global 'cp' (Phase H.2 no clamp)", {
  expect_silent(
    res <- gdpar:::resolve_parametrization_multi(
      parametrization             = "cp",
      parametrization_a           = NULL,
      parametrization_W           = NULL,
      parametrization_aggregation = NULL,
      prior = NULL, stan_data = NULL, amm = mk_dummy_amm(2L),
      preflight_seed = NULL, verbose = FALSE
    )
  )
  expect_true(res$cp_a)
  expect_true(res$cp_W)
  expect_identical(res$cp_a_per_k, c(TRUE, TRUE))
  expect_identical(res$cp_W_per_k, c(TRUE, TRUE))
  expect_null(res$report)
})

test_that("resolver honors parametrization_a = 'cp' (Phase H.2 unblock)", {
  expect_silent(
    res <- gdpar:::resolve_parametrization_multi(
      parametrization             = "auto",
      parametrization_a           = "cp",
      parametrization_W           = "ncp",
      parametrization_aggregation = NULL,
      prior = NULL, stan_data = NULL, amm = mk_dummy_amm(3L),
      preflight_seed = NULL, verbose = FALSE
    )
  )
  expect_true(res$cp_a)
  expect_false(res$cp_W)
  expect_identical(res$cp_a_per_k, c(TRUE, TRUE, TRUE))
  expect_identical(res$cp_W_per_k, c(FALSE, FALSE, FALSE))
})

test_that("resolver accepts aggregation = 'per_k' (Phase H.2 unblock)", {
  expect_silent(
    res <- gdpar:::resolve_parametrization_multi(
      parametrization             = "ncp",
      parametrization_a           = NULL,
      parametrization_W           = NULL,
      parametrization_aggregation = "per_k",
      prior = NULL, stan_data = NULL, amm = mk_dummy_amm(2L),
      preflight_seed = NULL, verbose = FALSE
    )
  )
  expect_identical(res$aggregation, "per_k")
})

test_that("resolver: both per-component explicit skips preflight (multi)", {
  res <- gdpar:::resolve_parametrization_multi(
    parametrization             = "auto",
    parametrization_a           = "ncp",
    parametrization_W           = "ncp",
    parametrization_aggregation = NULL,
    prior = NULL, stan_data = NULL, amm = mk_dummy_amm(2L),
    preflight_seed = NULL, verbose = FALSE
  )
  expect_false(res$cp_a)
  expect_false(res$cp_W)
  expect_null(res$report)
})


# ---- resolver: error paths --------------------------------------

test_that("resolver: unknown aggregation aborts with input_error", {
  expect_error(
    gdpar:::resolve_parametrization_multi(
      parametrization             = "auto",
      parametrization_a           = NULL,
      parametrization_W           = NULL,
      parametrization_aggregation = "weird",
      prior = NULL, stan_data = NULL, amm = mk_dummy_amm(2L),
      preflight_seed = NULL, verbose = FALSE
    ),
    class = "gdpar_input_error"
  )
})


# ---- preflight_parametrization_multi: argument validation -------

test_that("preflight_parametrization_multi rejects amm without p", {
  bad_amm <- structure(list(p = NULL), class = "amm_spec")
  expect_error(
    gdpar:::preflight_parametrization_multi(
      gdpar_prior(), stan_data = list(), amm = bad_amm
    ),
    class = "gdpar_internal_error"
  )
})

test_that("preflight_parametrization_multi rejects unknown aggregation", {
  expect_error(
    gdpar:::preflight_parametrization_multi(
      gdpar_prior(), stan_data = list(), amm = mk_dummy_amm(2L),
      aggregation = "weird"
    ),
    class = "gdpar_input_error"
  )
})

test_that("preflight_parametrization_multi short-circuits when no component needs decision", {
  # use_a = 0, use_W = 0 => nothing to decide => no Stan, NCP defaults.
  empty_stan_data <- list(
    use_a = 0L, use_W = 0L,
    J_a_per_k = rep(0L, 2L), W_per_k_dim = 0L, d = 0L,
    theta_anchor = c(0, 0)
  )
  res <- gdpar:::preflight_parametrization_multi(
    gdpar_prior(), empty_stan_data, mk_dummy_amm(2L)
  )
  expect_false(res$cp_a)
  expect_false(res$cp_W)
  expect_identical(res$cp_a_per_k, c(FALSE, FALSE))
  expect_identical(res$cp_W_per_k, c(FALSE, FALSE))
  expect_null(res$cs_model_ncp)
  expect_s3_class(res$report, "gdpar_preflight_report")
  expect_true(all(res$report$per_dim$decision == "absent"))
  expect_identical(res$report$global$global_decision,
                   c("absent", "absent"))
})


# ---- gdpar() exposes the new arg ---------------------------------

test_that("gdpar() rejects unknown parametrization_aggregation values", {
  df <- make_gaussian_data(n = 50, seed = 1)
  df$y <- cbind(df$y, df$y)
  spec <- amm_spec(p = 2L, dims = dimwise(a = ~ x1 + x2, b = NULL))
  fam <- gdpar_family_multi("gaussian", p = 2L)
  expect_error(
    gdpar(
      formula = y ~ x1 + x2, family = fam, amm = spec, data = df,
      parametrization_aggregation = "weird"
    )
  )
})


# ---- env-gated: end-to-end multi preflight + gdpar() -------------

test_that("preflight_parametrization_multi runs end-to-end p=2 NCP", {
  testthat::skip_if(Sys.getenv("GDPAR_F6_FIT") != "1",
    "Set GDPAR_F6_FIT=1 to run the multi preflight smoke (compiles a Stan model; takes minutes).")
  skip_if_no_cmdstan()
  set.seed(42)
  n <- 80L
  x1 <- rnorm(n); x2 <- rnorm(n)
  x1_c <- x1 - mean(x1); x2_c <- x2 - mean(x2)
  y1 <- 0.5 + 0.3 * x1_c - 0.2 * x2_c + rnorm(n, sd = 0.5)
  y2 <- 0.4 + 0.2 * x1_c - 0.1 * x2_c + rnorm(n, sd = 0.5)
  df <- data.frame(x1 = x1, x2 = x2)
  y_mat <- cbind(y1, y2)
  spec <- amm_spec(p = 2L, dims = dimwise(a = ~ x1 + x2, b = NULL))
  fam <- gdpar_family_multi("gaussian", p = 2L)
  design <- gdpar:::build_amm_design(
    spec, df, formula_rhs = ~ x1 + x2
  )
  stan_data <- gdpar:::assemble_stan_data(
    design, fam, spec, y_mat, theta_anchor = 0
  )
  res <- gdpar:::preflight_parametrization_multi(
    gdpar_prior(), stan_data, spec,
    preflight_seed = 42L, verbose = FALSE,
    aggregation = "any_ncp"
  )
  expect_s3_class(res$report, "gdpar_preflight_report")
  expect_identical(nrow(res$report$per_dim), 4L)
  expect_true(all(res$report$per_dim$component %in% c("a", "W")))
  expect_true(all(res$report$per_dim$dim %in% c(1L, 2L)))
  expect_true(all(res$report$per_dim$decision[
    res$report$per_dim$component == "W"] == "absent"))
  expect_identical(length(res$cp_a_per_k), 2L)
  expect_identical(length(res$cp_W_per_k), 2L)
  expect_false(res$cp_W)
})

test_that("gdpar() multi fit attaches preflight_report when auto active", {
  testthat::skip_if(Sys.getenv("GDPAR_F6_FIT") != "1",
    "Set GDPAR_F6_FIT=1 to run the multi end-to-end gdpar() preflight (compiles a Stan model; takes minutes).")
  skip_if_no_cmdstan()
  set.seed(42)
  n <- 80L
  x1 <- rnorm(n); x2 <- rnorm(n)
  x1_c <- x1 - mean(x1); x2_c <- x2 - mean(x2)
  y1 <- 0.5 + 0.3 * x1_c - 0.2 * x2_c + rnorm(n, sd = 0.5)
  y2 <- 0.4 + 0.2 * x1_c - 0.1 * x2_c + rnorm(n, sd = 0.5)
  df <- data.frame(x1 = x1, x2 = x2)
  df$y <- cbind(y1, y2)
  spec <- amm_spec(p = 2L, dims = dimwise(a = ~ x1 + x2, b = NULL))
  fam <- gdpar_family_multi("gaussian", p = 2L)
  fit <- suppressMessages(suppressWarnings(
    gdpar(
      formula = y ~ x1 + x2, family = fam, amm = spec, data = df,
      chains = 2L, iter_warmup = 200L, iter_sampling = 200L,
      refresh = 0L, verbose = FALSE, seed = 42L,
      parametrization = "auto", skip_id_check = TRUE
    )
  ))
  expect_s3_class(fit, "gdpar_fit")
  expect_s3_class(fit$parametrization$report,
                  "gdpar_preflight_report")
  expect_identical(length(fit$parametrization$cp_a_per_k), 2L)
  expect_identical(length(fit$parametrization$cp_W_per_k), 2L)
  expect_false(fit$parametrization$cp_W)
  # H.2 contract: cp_a follows the preflight aggregation; the cp_a
  # field is now logical (no longer auto-clamped to FALSE).
  expect_type(fit$parametrization$cp_a, "logical")
  expect_identical(fit$parametrization$meta$aggregation, "any_ncp")
})

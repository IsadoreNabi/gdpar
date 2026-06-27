# Block RG, RG.6 part ii (session B9.31): the integration layer that wires the
# geometry-adaptive orchestrator to the package fit engine.
#   * Ungated (algebraic / logic, deterministic, no cmdstan): the K front-end
#     resolver, the bridge-core assembly, input validation, the print methods.
#   * Gated (GDPAR_RUN_GEOMETRY_PILOTS): the end-to-end wrapper and bridge over a
#     real compiled cmdstan model.
# The refactor's bit-identity (.gdpar_K_build extracted from .gdpar_K) is proven
# separately by the K-path goldens (test-golden_regression_K2.R, the 8.3.9
# compare-path), which go through gdpar() -> .gdpar_K() -> .gdpar_K_build().

.skip_unless_geometry_pilots <- function() {
  if (!nzchar(Sys.getenv("GDPAR_RUN_GEOMETRY_PILOTS"))) {
    testthat::skip("Set GDPAR_RUN_GEOMETRY_PILOTS=1 to enable the bridge pilots.")
  }
  testthat::skip_if_not_installed("cmdstanr")
  testthat::skip_if_not_installed("posterior")
}

# --- ungated: K front-end resolver ----------------------------------------

test_that(".gdpar_geom_fit_resolve_K resolves the formula-set (gdpar_bf) form", {
  fs <- gdpar_bf(y ~ a(x), sigma ~ a(z))
  res <- .gdpar_geom_fit_resolve_K(fs, gdpar_family("gaussian"), amm_spec(),
                                   NULL, data.frame())
  expect_equal(res$K, 2L)
  expect_equal(res$outcome_name, "y")
  expect_named(res$amm_list_canonical, c("mu", "sigma"))
  expect_s3_class(res$family, "gdpar_family")
  expect_null(res$family_id_k_vector)
})

test_that(".gdpar_geom_fit_resolve_K resolves a heterogeneous named-family list", {
  fs <- gdpar_bf(y ~ a(x), sigma ~ a(z))
  fam <- list(mu = gdpar_family("gaussian"), sigma = gdpar_family("gaussian"))
  res <- .gdpar_geom_fit_resolve_K(fs, fam, amm_spec(), NULL, data.frame())
  expect_equal(res$K, 2L)
  expect_false(is.null(res$family_id_k_vector))
})

test_that(".gdpar_geom_fit_resolve_K rejects K = 1 and unsupported forms", {
  expect_error(
    .gdpar_geom_fit_resolve_K(gdpar_bf(y ~ a(x)), gdpar_family("gaussian"),
                              amm_spec(), NULL, data.frame()),
    class = "gdpar_input_error")
  # A plain two-sided formula with no AMM wrappers and a default amm is not a
  # K-individual spec.
  expect_error(
    .gdpar_geom_fit_resolve_K(y ~ x, gdpar_family("gaussian"), amm_spec(),
                              NULL, data.frame()),
    class = "gdpar_input_error")
})

# --- ungated: bridge-core assembly + validation ---------------------------

test_that(".gdpar_geom_bridge_core assembles a bridge from a closure engine", {
  inst <- list(log_prob = function(t) -0.5 * sum(t^2),
               grad_log_prob = function(t) -t, dim = 3L)
  br <- .gdpar_geom_bridge_core(model = "dummy_model",
                               stan_data = list(N = 1L), dim = 3L,
                               reference = c(0, 0, 0), engine_fit = inst)
  expect_s3_class(br, "gdpar_geom_bridge")
  expect_equal(br$dim, 3L)
  expect_s3_class(br$geom_target, "gdpar_geom_target")
  expect_equal(br$geom_target$dim, 3L)
  # the target carries the (re-samplable) model + dim + data for the diagnostic
  expect_identical(br$target$model, "dummy_model")
  expect_equal(br$target$dim, 3L)
  expect_equal(br$geom_target$grad_log_prob(c(1, 2, 3)), c(-1, -2, -3))
})

test_that(".gdpar_geom_bridge_core rejects a non-function fisher", {
  inst <- list(log_prob = function(t) -0.5 * sum(t^2),
               grad_log_prob = function(t) -t, dim = 2L)
  expect_error(
    .gdpar_geom_bridge_core(model = "m", stan_data = NULL, dim = 2L,
                            fisher = "not_a_function", engine_fit = inst),
    class = "gdpar_input_error")
})

test_that("gdpar_geom_bridge validates its input", {
  expect_error(gdpar_geom_bridge(list(a = 1)), class = "gdpar_input_error")
  fake <- structure(list(fit = NULL, stan_data = list()), class = "gdpar_fit")
  expect_error(gdpar_geom_bridge(fake), class = "gdpar_input_error")
})

test_that("print methods are well-behaved", {
  inst <- list(log_prob = function(t) -0.5 * sum(t^2),
               grad_log_prob = function(t) -t, dim = 2L)
  br <- .gdpar_geom_bridge_core(model = "m", stan_data = NULL, dim = 2L,
                               engine_fit = inst)
  expect_output(print(br), "gdpar_geom_bridge")
  ff <- structure(list(status = "resolved", level = "euclidean_diagonal",
                       draws = matrix(0, 4, 2), K = 2L),
                  class = "gdpar_geom_fit")
  expect_output(print(ff), "gdpar_geom_fit")
  expect_invisible(print(br))
})

# --- gated: end-to-end over a real compiled cmdstan model -----------------

.build_bridge_K2_data <- function(n = 80L, seed = 909L) {
  set.seed(seed)
  x <- rnorm(n); z <- rnorm(n)
  mu_true <- 0.5 + 0.8 * (x - mean(x))
  log_sigma <- -0.2 + 0.4 * (z - mean(z))
  y <- rnorm(n, mu_true, exp(log_sigma))
  data.frame(y = y, x = x, z = z)
}

test_that("gdpar_geom_bridge bridges a fitted K=2 model end to end", {
  .skip_unless_geometry_pilots()
  d <- .build_bridge_K2_data()
  fit <- suppressMessages(suppressWarnings(
    gdpar(gdpar_bf(y ~ a(x), sigma ~ a(z)), data = d,
          family = gdpar_family("gaussian"), chains = 1L,
          iter_warmup = 200L, iter_sampling = 200L, refresh = 0L,
          seed = 313L, skip_id_check = TRUE, verbose = FALSE)))
  br <- gdpar_geom_bridge(fit)
  expect_s3_class(br, "gdpar_geom_bridge")
  expect_gt(br$dim, 0L)
  expect_length(br$reference, br$dim)
  expect_true(inherits(br$model, "CmdStanModel"))
  # the engine target evaluates the unconstrained log-density and its gradient
  lp <- br$geom_target$log_prob(br$reference)
  gl <- br$geom_target$grad_log_prob(br$reference)
  expect_true(is.finite(lp))
  expect_length(gl, br$dim)
})

test_that("gdpar_geom_fit runs the orchestrator on a fresh K=2 build", {
  .skip_unless_geometry_pilots()
  d <- .build_bridge_K2_data()
  b <- gdpar_geom_orchestrate_budget()
  b$tune_epsilon <- FALSE
  b$probe_warmup <- 80L; b$probe_iter <- 80L
  b$full_warmup <- 100L; b$full_iter <- 100L
  b$max_rounds <- 2L
  res <- suppressMessages(suppressWarnings(
    gdpar_geom_fit(gdpar_bf(y ~ a(x), sigma ~ a(z)), data = d,
                   family = gdpar_family("gaussian"), skip_id_check = TRUE,
                   budget = b, n_grid = 1, seed = 42L, verbose = FALSE)))
  expect_s3_class(res, "gdpar_geom_fit")
  expect_equal(res$K, 2L)
  expect_true(res$status %in% c("resolved", "certified_limit", "out_of_scope"))
  expect_s3_class(res$orchestration, "gdpar_geom_orchestration")
})

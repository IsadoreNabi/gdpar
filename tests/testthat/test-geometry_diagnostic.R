# Block RG, RG.1: posterior-geometry diagnostic + synthetic suite.
#
# Default tests are Stan-free: suite contract, analytic-gradient cross-check
# (numerical layer), deterministic classifier branches, input validation. The
# end-to-end cmdstan smoke is gated by GDPAR_RUN_GEOMETRY_PILOTS.

# --- helpers ---------------------------------------------------------------

fd_grad_geom <- function(f, theta, h = 1e-5) {
  vapply(seq_along(theta), function(j) {
    tp <- theta; tp[j] <- tp[j] + h
    tm <- theta; tm[j] <- tm[j] - h
    (f(tp) - f(tm)) / (2 * h)
  }, numeric(1))
}

# A benign baseline signal set; override fields to exercise each branch.
benign_signals <- function(...) {
  base <- list(
    divergent_rate = 0, ebfmi_min = 1.0, treedepth_sat_rate = 0,
    condition_number = 2, lambda_max_cov = 2, step_scale_ratio = 1.0,
    mean_leapfrog = 8, multimodality = 0.2, heavy_kurtosis = 0.1,
    boundary_proximity = 0, elapsed = 1, n_sampling = 600
  )
  modifyList(base, list(...))
}

# --- suite contract --------------------------------------------------------

test_that("gdpar_geometry_suite returns the eight taxonomy targets", {
  suite <- gdpar_geometry_suite()
  expect_length(suite, 8L)
  expect_setequal(
    names(suite),
    c("G0_isotropic", "G1_anisotropic", "G2_funnel", "G3_heavy_tails",
      "G4_quasi_deterministic", "G5_multimodal", "G6_boundary",
      "G7_flat_direction")
  )
  for (tg in suite) {
    expect_s3_class(tg, "gdpar_geometry_target")
    expect_true(all(c("id", "pathology", "geometry_remedy", "culprit",
                      "difficulty_scales_with_n", "make") %in% names(tg)))
    inst <- tg$make(tg$default_n, tg$default_difficulty)
    expect_true(all(c("stan_code", "stan_data", "log_prob", "grad_log_prob",
                      "dim", "param_names") %in% names(inst)))
    expect_length(inst$param_names, inst$dim)
    expect_true(is.function(inst$log_prob))
    expect_true(is.function(inst$grad_log_prob))
  }
})

test_that("ground-truth labels match the taxonomy", {
  suite <- gdpar_geometry_suite()
  expect_identical(suite$G0_isotropic$pathology, "isotropic")
  expect_identical(suite$G2_funnel$pathology, "funnel")
  expect_identical(suite$G4_quasi_deterministic$pathology,
                   "quasi_deterministic")
  expect_identical(suite$G7_flat_direction$pathology, "flat_direction")
  # Only the quasi-deterministic case grows with n.
  grows <- vapply(suite, function(tg) tg$difficulty_scales_with_n, logical(1))
  expect_identical(names(grows)[grows], "G4_quasi_deterministic")
})

test_that("which= subsets and rejects invalid ids", {
  sub <- gdpar_geometry_suite(which = c("G2_funnel", "G4_quasi_deterministic"))
  expect_length(sub, 2L)
  expect_error(gdpar_geometry_suite(which = "nope"),
               class = "gdpar_input_error")
})

# --- numerical layer: analytic gradients --------------------------------

test_that("analytic gradients match finite differences across the suite", {
  suite <- gdpar_geometry_suite()
  set.seed(202)
  for (nm in names(suite)) {
    tg <- suite[[nm]]
    inst <- tg$make(tg$default_n, tg$default_difficulty)
    theta <- stats::rnorm(inst$dim, sd = 0.5)
    ga <- inst$grad_log_prob(theta)
    gn <- fd_grad_geom(inst$log_prob, theta)
    expect_lt(max(abs(ga - gn)), 1e-4)
  }
})

# --- deterministic classifier branches (Stan-free) -------------------------

test_that("classifier recovers each pathology from crafted signals", {
  th <- gdpar_geometry_thresholds()
  flat_n <- list(slope = 0, grows_with_n = FALSE)
  grow_n <- list(slope = 1, grows_with_n = TRUE)
  empty <- gdpar:::.gdpar_geom_empty_culprit()
  cl <- function(sig, nc) {
    gdpar:::.gdpar_geom_classify(sig, nc, empty, th)$pathology
  }
  expect_identical(cl(benign_signals(), flat_n), "isotropic")
  expect_identical(
    cl(benign_signals(condition_number = 200, lambda_max_cov = 200), flat_n),
    "anisotropic")
  # funnel, energy path: high kurtosis with a very low energy fraction.
  expect_identical(
    cl(benign_signals(ebfmi_min = 0.13, heavy_kurtosis = 20,
                      condition_number = 4), flat_n),
    "funnel")
  # funnel, divergence path: a divergence neck with tiny steps.
  expect_identical(
    cl(benign_signals(divergent_rate = 0.05, step_scale_ratio = 0.02,
                      condition_number = 4), flat_n),
    "funnel")
  # heavy tails: high kurtosis but moderate energy and low condition number.
  expect_identical(
    cl(benign_signals(ebfmi_min = 0.5, heavy_kurtosis = 10,
                      condition_number = 2.6), flat_n),
    "heavy_tails")
  expect_identical(
    cl(benign_signals(condition_number = 500, lambda_max_cov = 1), grow_n),
    "quasi_deterministic")
  expect_identical(cl(benign_signals(multimodality = 5), flat_n), "multimodal")
  expect_identical(cl(benign_signals(boundary_proximity = 0.5), flat_n),
                   "boundary")
  expect_identical(
    cl(benign_signals(lambda_max_cov = 1e4, condition_number = 2e4), flat_n),
    "flat_direction")
})

test_that("remedy lookup maps the count case to sub-Riemannian", {
  expect_identical(
    gdpar:::.gdpar_geom_remedy_for("quasi_deterministic")$remedy,
    "sub_riemannian")
  expect_identical(
    gdpar:::.gdpar_geom_remedy_for("flat_direction")$remedy,
    "reparam_eliminate")
  expect_identical(gdpar:::.gdpar_geom_remedy_for("funnel")$remedy,
                   "riemannian")
})

test_that("difficulty curve flags n-growth only when present", {
  th <- gdpar_geometry_thresholds()
  flat <- data.frame(n = c(10, 100, 1000),
                     condition_number = c(50, 52, 49), failed = FALSE)
  grow <- data.frame(n = c(10, 100, 1000),
                     condition_number = c(10, 100, 1000), failed = FALSE)
  expect_false(gdpar:::.gdpar_geom_difficulty_curve(flat, th)$grows_with_n)
  expect_true(gdpar:::.gdpar_geom_difficulty_curve(grow, th)$grows_with_n)
  one <- data.frame(n = 100, condition_number = 50, failed = FALSE)
  expect_true(is.na(gdpar:::.gdpar_geom_difficulty_curve(one, th)$grows_with_n))
})

# --- input validation ------------------------------------------------------

test_that("gdpar_geometry_diagnostic validates arguments", {
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("posterior")
  expect_error(
    gdpar_geometry_diagnostic(list(foo = 1), verbose = FALSE),
    class = "gdpar_input_error")
  suite <- gdpar_geometry_suite()
  expect_error(
    gdpar_geometry_diagnostic(suite$G0_isotropic, adapt_delta = 2,
                              verbose = FALSE),
    class = "gdpar_input_error")
})

test_that("thresholds are a named numeric list", {
  th <- gdpar_geometry_thresholds()
  expect_true(is.list(th))
  expect_true(all(vapply(th, is.numeric, logical(1))))
  expect_true("flat_var_high" %in% names(th))
})

# --- gated end-to-end cmdstan smoke ----------------------------------------

test_that("[gated] diagnostic classifies G0 and G2 through real pilots", {
  if (!nzchar(Sys.getenv("GDPAR_RUN_GEOMETRY_PILOTS"))) {
    skip("Set GDPAR_RUN_GEOMETRY_PILOTS=1 to enable cmdstan pilots.")
  }
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("posterior")
  has_cmdstan <- tryCatch({ cmdstanr::cmdstan_version(); TRUE },
                          error = function(e) FALSE)
  if (!has_cmdstan) skip("cmdstan toolchain not available.")

  suite <- gdpar_geometry_suite()
  d0 <- gdpar_geometry_diagnostic(suite$G0_isotropic, verbose = FALSE)
  expect_s3_class(d0, "gdpar_geometry_diagnostic")
  expect_identical(d0$pathology, "isotropic")
  expect_true(d0$correct)

  d2 <- gdpar_geometry_diagnostic(suite$G2_funnel, verbose = FALSE)
  expect_identical(d2$pathology, "funnel")
  expect_true(d2$correct)
})

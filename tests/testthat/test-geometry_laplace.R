# Block RG, RG.7 step 4 (D99): the automated Laplace fallback.
#
# Default (ungated) tests are pure R (no cmdstan): gdpar_geom_laplace() on cheap
# closure targets with analytic / finite-difference curvature (the mode, the
# positive-definite precision, the draw covariance ~ M^{-1}, the fidelity fields
# and the scalar label), the saddle guard, input validation, and the
# orchestrator's attach-fallback finalizer (status relabelling, the out-of-scope
# PD gate, additivity = the bit-identical default). The end-to-end layer (a real
# orchestration ending in a certified limit, with the fallback firing) is gated
# by GDPAR_RUN_GEOMETRY_PILOTS, the same gate as the other RG pilots.

# --- helpers ---------------------------------------------------------------

# A correlated Gaussian: the Laplace approximation is EXACT, so the fidelity
# label is "good". Constant analytic Hessian -H = A (positive-definite).
gauss_target <- function(mu = c(1, -0.5), A = matrix(c(2, 0.8, 0.8, 1), 2, 2),
                         with_hessian = TRUE) {
  lp <- function(th) -0.5 * as.numeric(t(th - mu) %*% A %*% (th - mu))
  gr <- function(th) -as.numeric(A %*% (th - mu))
  he <- if (with_hessian) function(th) -A else NULL
  gdpar_geom_target(log_prob = lp, grad_log_prob = gr, dim = length(mu),
                    hessian = he)
}

# An independent Student-t with small nu: heavy tails the Laplace Gaussian cannot
# match, so the importance weights p/q blow up in the tails and the label
# degrades. Mode at 0; analytic Hessian.
studentt_target <- function(d = 3L, nu = 2.5) {
  lp <- function(th) -((nu + 1) / 2) * sum(log1p(th^2 / nu))
  gr <- function(th) -((nu + 1) / 2) * (2 * th / nu) / (1 + th^2 / nu)
  he <- function(th) {
    u <- 1 + th^2 / nu
    diag(-((nu + 1) / nu) * (1 - 2 * th^2 / nu / u) / u, d)
  }
  gdpar_geom_target(log_prob = lp, grad_log_prob = gr, dim = d, hessian = he)
}

# --- gdpar_geom_laplace on the exact-Gaussian case -------------------------

test_that("the Laplace of a correlated Gaussian is exact and labelled 'good'", {
  mu <- c(1, -0.5); A <- matrix(c(2, 0.8, 0.8, 1), 2, 2)
  lap <- gdpar_geom_laplace(gauss_target(mu, A), draws = 2000L, seed = 1L)
  expect_s3_class(lap, "gdpar_geom_laplace")
  # mode hits the analytic centre; the precision is A; the covariance is A^{-1}.
  expect_lt(max(abs(lap$mode - mu)), 1e-6)
  expect_lt(max(abs(lap$M - A)), 1e-6)
  expect_lt(max(abs(lap$cov - solve(A))), 1e-6)
  expect_true(lap$all_pos)
  expect_identical(lap$method, "exact_hessian")
  expect_true(is.finite(lap$cond) && is.finite(lap$cond_unfloored))
  # iid draws reproduce M^{-1}; fidelity is perfect, the label is "good".
  expect_lt(max(abs(stats::cov(lap$draws) - solve(A))), 0.15)
  fq <- lap$fit_quality
  expect_true(is.finite(fq$ess_is) && is.finite(fq$logdrop_mean))
  expect_gt(fq$ess_is_frac, 0.9)
  expect_identical(lap$fit_quality_label, "good")
})

test_that("the draws argument controls the returned draw matrix", {
  tgt <- gauss_target()
  lap0 <- gdpar_geom_laplace(tgt, draws = 0L, seed = 2L)
  expect_identical(nrow(lap0$draws), 0L)
  expect_identical(ncol(lap0$draws), 2L)
  # the label is still computed (on the internal fidelity sample).
  expect_true(lap0$fit_quality_label %in% c("good", "poor", "very_poor"))
  lap50 <- gdpar_geom_laplace(tgt, draws = 50L, seed = 2L)
  expect_identical(nrow(lap50$draws), 50L)
})

test_that("a Laplace with no exact Hessian uses finite differences", {
  mu <- c(1, -0.5); A <- matrix(c(2, 0.8, 0.8, 1), 2, 2)
  lap <- gdpar_geom_laplace(gauss_target(mu, A, with_hessian = FALSE),
                            draws = 0L, seed = 3L)
  expect_identical(lap$method, "finite_difference")
  expect_lt(max(abs(lap$mode - mu)), 1e-5)
  expect_lt(max(abs(lap$M - A)), 1e-4)
  expect_true(lap$all_pos)
})

# --- the fidelity label degrades on a genuinely non-Gaussian target --------

test_that("a heavy-tailed posterior gets a non-'good' fidelity label", {
  lap <- gdpar_geom_laplace(studentt_target(d = 3L, nu = 2.5),
                            draws = 0L, seed = 4L)
  expect_lt(max(abs(lap$mode)), 1e-3)          # mode at the origin.
  expect_true(lap$all_pos)                      # PD curvature at the mode...
  expect_false(identical(lap$fit_quality_label, "good"))  # ...but a crude Gaussian.
})

# --- the saddle guard ------------------------------------------------------

test_that("a non-positive-definite curvature is flagged 'very_poor' with a warning", {
  # log p = 0.5 t1^2 - 0.5 t2^2: a saddle at the origin (grad = 0), -Hessian
  # indefinite. With climb = FALSE the origin is read directly.
  saddle <- gdpar_geom_target(
    log_prob = function(th) 0.5 * th[1]^2 - 0.5 * th[2]^2,
    grad_log_prob = function(th) c(th[1], -th[2]),
    hessian = function(th) diag(c(1, -1)), dim = 2L)
  expect_warning(
    lap <- gdpar_geom_laplace(saddle, reference = c(0, 0), climb = FALSE,
                              seed = 5L),
    class = "gdpar_geometry_warning")
  expect_false(lap$all_pos)
  expect_identical(lap$fit_quality_label, "very_poor")
})

# --- input validation ------------------------------------------------------

test_that("gdpar_geom_laplace validates draws and reference", {
  tgt <- gauss_target()
  expect_error(gdpar_geom_laplace(tgt, draws = -1L), class = "gdpar_input_error")
  expect_error(gdpar_geom_laplace(tgt, reference = c(1, 2, 3)),
               class = "gdpar_input_error")
})

# --- print method ----------------------------------------------------------

test_that("the print method reports the label and the fidelity", {
  lap <- gdpar_geom_laplace(gauss_target(), seed = 6L)
  expect_output(print(lap), "gdpar_geom_laplace")
  expect_output(print(lap), "fit-quality")
})

# --- the orchestrator attach-fallback finalizer ----------------------------

quiet_inform <- function(msg) invisible(NULL)

test_that("the finalizer relabels a certified limit and attaches the Laplace", {
  obj0 <- structure(
    list(status = "certified_limit", certificate = list(verdict = "x"),
         ledger = list(), diagnosis = NULL, best = NULL,
         budget_spent = list(fits = 1L, seconds = 0)),
    class = c("gdpar_geom_orchestration", "list"))
  obj1 <- gdpar:::.gdpar_orch_attach_laplace(
    obj0, gauss_target(), reference = c(0, 0), draws = 0L, seed = 7L,
    inform = quiet_inform)
  expect_identical(obj1$status, "certified_limit_laplace")
  expect_s3_class(obj1$laplace, "gdpar_geom_laplace")
  # purely additive: everything else is untouched (=> default FALSE is
  # bit-identical, since the finalizer is the only writer).
  expect_identical(obj1[setdiff(names(obj0), "status")],
                   obj0[setdiff(names(obj0), "status")])
})

test_that("the out-of-scope path attaches only when the curvature is PD", {
  pd_obj <- structure(list(status = "out_of_scope"),
                      class = c("gdpar_geom_orchestration", "list"))
  out_pd <- gdpar:::.gdpar_orch_attach_laplace(
    pd_obj, gauss_target(), reference = c(0, 0), draws = 0L, seed = 8L,
    inform = quiet_inform)
  expect_identical(out_pd$status, "out_of_scope_laplace")
  expect_s3_class(out_pd$laplace, "gdpar_geom_laplace")
  # a non-PD (saddle) curvature: no Laplace attached, status unchanged.
  saddle <- gdpar_geom_target(
    log_prob = function(th) 0.5 * th[1]^2 - 0.5 * th[2]^2,
    grad_log_prob = function(th) c(th[1], -th[2]),
    hessian = function(th) diag(c(1, -1)), dim = 2L)
  np_obj <- structure(list(status = "out_of_scope"),
                      class = c("gdpar_geom_orchestration", "list"))
  out_np <- suppressWarnings(gdpar:::.gdpar_orch_attach_laplace(
    np_obj, saddle, reference = c(0, 0), draws = 0L, seed = 8L,
    inform = quiet_inform))
  expect_identical(out_np$status, "out_of_scope")
  expect_null(out_np$laplace)
})

# --- gated end-to-end: a real certified limit firing the fallback ----------

test_that("laplace_fallback turns a certified limit into certified_limit_laplace", {
  if (!nzchar(Sys.getenv("GDPAR_RUN_GEOMETRY_PILOTS"))) {
    skip("Set GDPAR_RUN_GEOMETRY_PILOTS=1 to enable the orchestration pilot.")
  }
  suite <- gdpar_geometry_suite()
  b <- gdpar_geom_orchestrate_budget()
  b$probe_warmup <- 80L; b$probe_iter <- 80L
  b$full_warmup <- 120L; b$full_iter <- 120L
  b$tune_iter <- 40L; b$L <- 18L; b$max_rounds <- 4L
  # The quasi-deterministic canyon WITHOUT the Fisher certifies a limit; with the
  # fallback it returns a labelled Laplace instead of a bare certificate.
  r <- gdpar_geom_orchestrate(suite$G4_quasi_deterministic,
         n_grid = c(10, 100, 1000), budget = b, laplace_fallback = TRUE,
         laplace_draws = 100L, pilot_warmup = 80L, pilot_sampling = 80L,
         verbose = FALSE)
  expect_identical(r$status, "certified_limit_laplace")
  expect_s3_class(r$laplace, "gdpar_geom_laplace")
  expect_true(r$laplace$fit_quality_label %in% c("good", "poor", "very_poor"))
  expect_identical(nrow(r$laplace$draws), 100L)
  # the default (no fallback) keeps the bare certified limit.
  r0 <- gdpar_geom_orchestrate(suite$G4_quasi_deterministic,
          n_grid = c(10, 100, 1000), budget = b, pilot_warmup = 80L,
          pilot_sampling = 80L, verbose = FALSE)
  expect_identical(r0$status, "certified_limit")
  expect_null(r0$laplace)
})

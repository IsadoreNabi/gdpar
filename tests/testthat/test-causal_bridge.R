# Tests for gdpar_causal_bridge (Sub-phase 8.5.A; T-learner AMM-side bridge).
#
# Layered coverage:
#   (i)   unit tests of the compatibility validators run unconditionally on
#         mock gdpar_fit objects (no Stan compilation, no posterior draws).
#   (ii)  unit tests of the internal helpers `.align_bridge_draws` and
#         `.summarize_cate` covering matrix and 3D-array shapes.
#   (iii) opt-in smoke + algorithmic-equivalence test gated by env var
#         GDPAR_RUN_STAN_SMOKE_CAUSAL_BRIDGE and
#         skip_if_not_installed("cmdstanr"), consistent with the per-feature
#         opt-in pattern of the suite (cf. test-bvm_check.R Section iii).

# ---- (i) Compatibility validators (mock fits, unconditional) ------------

make_mock_bridge_fit <- function(family_name = "gaussian",
                                  link = "identity",
                                  amm = NULL,
                                  anchor = 0,
                                  use_groups = 0L,
                                  path = "bayes",
                                  K = NULL,
                                  p = NULL,
                                  param_specs = NULL) {
  if (is.null(amm)) {
    amm <- amm_spec(a = ~ x1)
  }
  fit <- list(
    path      = path,
    family    = list(name = family_name, link = link,
                     param_specs = param_specs),
    amm       = amm,
    anchor    = anchor,
    stan_data = list(use_groups = as.integer(use_groups)),
    identifiability_report = list(passed = TRUE),
    call      = call("gdpar")
  )
  if (!is.null(K)) fit$K <- as.integer(K)
  if (!is.null(p)) fit$p <- as.integer(p)
  class(fit) <- c("gdpar_fit", "list")
  fit
}

test_that("gdpar_causal_bridge rejects non-gdpar_fit inputs", {
  expect_error(
    gdpar_causal_bridge(list(), make_mock_bridge_fit()),
    class = "gdpar_input_error"
  )
  expect_error(
    gdpar_causal_bridge(make_mock_bridge_fit(), list()),
    class = "gdpar_input_error"
  )
})

test_that("gdpar_causal_bridge rejects level out of range", {
  fit_t <- make_mock_bridge_fit()
  fit_c <- make_mock_bridge_fit()
  expect_error(
    gdpar_causal_bridge(fit_t, fit_c, level = 0),
    class = "gdpar_input_error"
  )
  expect_error(
    gdpar_causal_bridge(fit_t, fit_c, level = 1),
    class = "gdpar_input_error"
  )
  expect_error(
    gdpar_causal_bridge(fit_t, fit_c, level = -0.1),
    class = "gdpar_input_error"
  )
})

test_that("gdpar_causal_bridge aborts when one fit is non-bayes path", {
  fit_t <- make_mock_bridge_fit(path = "bayes")
  fit_c <- make_mock_bridge_fit(path = "vcm")
  expect_error(
    gdpar_causal_bridge(fit_t, fit_c),
    class = "gdpar_unsupported_feature_error"
  )
})

test_that("gdpar_causal_bridge aborts when either fit is hierarchical", {
  fit_t <- make_mock_bridge_fit(use_groups = 1L)
  fit_c <- make_mock_bridge_fit(use_groups = 0L)
  expect_error(
    gdpar_causal_bridge(fit_t, fit_c),
    class = "gdpar_unsupported_feature_error"
  )
  fit_t <- make_mock_bridge_fit(use_groups = 0L)
  fit_c <- make_mock_bridge_fit(use_groups = 1L)
  expect_error(
    gdpar_causal_bridge(fit_t, fit_c),
    class = "gdpar_unsupported_feature_error"
  )
})

test_that("gdpar_causal_bridge aborts on mismatched family name or link", {
  fit_t <- make_mock_bridge_fit(family_name = "gaussian", link = "identity")
  fit_c <- make_mock_bridge_fit(family_name = "poisson", link = "log")
  expect_error(
    gdpar_causal_bridge(fit_t, fit_c),
    class = "gdpar_unsupported_feature_error"
  )
  fit_t <- make_mock_bridge_fit(family_name = "gaussian", link = "identity")
  fit_c <- make_mock_bridge_fit(family_name = "gaussian", link = "log")
  expect_error(
    gdpar_causal_bridge(fit_t, fit_c),
    class = "gdpar_unsupported_feature_error"
  )
})

test_that("gdpar_causal_bridge aborts on mismatched K or p", {
  fit_t <- make_mock_bridge_fit(K = 2L)
  fit_c <- make_mock_bridge_fit(K = 3L)
  expect_error(
    gdpar_causal_bridge(fit_t, fit_c),
    class = "gdpar_unsupported_feature_error"
  )
  fit_t <- make_mock_bridge_fit(p = 2L)
  fit_c <- make_mock_bridge_fit(p = 1L)
  expect_error(
    gdpar_causal_bridge(fit_t, fit_c),
    class = "gdpar_unsupported_feature_error"
  )
})

test_that("gdpar_causal_bridge aborts on mismatched anchor length or value", {
  fit_t <- make_mock_bridge_fit(anchor = c(0, 1))
  fit_c <- make_mock_bridge_fit(anchor = 0)
  expect_error(
    gdpar_causal_bridge(fit_t, fit_c),
    class = "gdpar_unsupported_feature_error"
  )
  fit_t <- make_mock_bridge_fit(anchor = 0.0)
  fit_c <- make_mock_bridge_fit(anchor = 0.5)
  expect_error(
    gdpar_causal_bridge(fit_t, fit_c),
    class = "gdpar_unsupported_feature_error"
  )
})

test_that("gdpar_causal_bridge aborts on mismatched AMM level", {
  fit_t <- make_mock_bridge_fit(amm = amm_spec(a = ~ x1))
  fit_c <- make_mock_bridge_fit(amm = amm_spec())
  expect_error(
    gdpar_causal_bridge(fit_t, fit_c),
    class = "gdpar_unsupported_feature_error"
  )
})

test_that("gdpar_causal_bridge aborts on mismatched covariate columns", {
  fit_t <- make_mock_bridge_fit(amm = amm_spec(a = ~ x1))
  fit_c <- make_mock_bridge_fit(amm = amm_spec(a = ~ x2))
  expect_error(
    gdpar_causal_bridge(fit_t, fit_c),
    class = "gdpar_unsupported_feature_error"
  )
})

test_that("gdpar_causal_bridge with NULL newdata aborts when no recovery", {
  fit_t <- make_mock_bridge_fit()
  fit_c <- make_mock_bridge_fit()
  fit_t$call <- call("gdpar", data = quote(.nonexistent_data_treat))
  fit_c$call <- call("gdpar", data = quote(.nonexistent_data_ctrl))
  expect_error(
    gdpar_causal_bridge(fit_t, fit_c),
    class = "gdpar_input_error"
  )
})

# ---- (ii) Internal helpers: shape handling ------------------------------

test_that(".align_bridge_draws preserves equal-sized inputs untouched", {
  pred_t <- matrix(seq_len(20L), nrow = 4L, ncol = 5L)
  pred_c <- matrix(seq_len(20L) + 100L, nrow = 4L, ncol = 5L)
  out <- expect_silent(gdpar:::.align_bridge_draws(pred_t, pred_c))
  expect_equal(out$S, 4L)
  expect_equal(out$treat, pred_t)
  expect_equal(out$ctrl, pred_c)
  # D50 (Sesion 18 Etapa 2 de Sesion 8.4): equal-sized inputs yield NA
  # in the `warning` field (no trimming occurred).
  expect_true(is.na(out$warning))
})

test_that(".align_bridge_draws trims to the shorter input with warning", {
  pred_t <- matrix(seq_len(15L), nrow = 5L, ncol = 3L)
  pred_c <- matrix(seq_len(9L), nrow = 3L, ncol = 3L)
  expect_warning(
    out <- gdpar:::.align_bridge_draws(pred_t, pred_c),
    class = "gdpar_diagnostic_warning"
  )
  expect_equal(out$S, 3L)
  expect_equal(nrow(out$treat), 3L)
  expect_equal(nrow(out$ctrl), 3L)
  # D50: trimming populates the `warning` field with the human-readable
  # notification message that was simultaneously emitted via gdpar_warn.
  expect_type(out$warning, "character")
  expect_false(is.na(out$warning))
  expect_match(out$warning, "trimmed posterior draws", fixed = TRUE)
})

test_that(".align_bridge_draws supports 3D arrays", {
  pred_t <- array(seq_len(60L), dim = c(6L, 5L, 2L))
  pred_c <- array(seq_len(40L), dim = c(4L, 5L, 2L))
  expect_warning(
    out <- gdpar:::.align_bridge_draws(pred_t, pred_c),
    class = "gdpar_diagnostic_warning"
  )
  expect_equal(out$S, 4L)
  expect_equal(dim(out$treat), c(4L, 5L, 2L))
  expect_equal(dim(out$ctrl),  c(4L, 5L, 2L))
  expect_false(is.na(out$warning))
})

# ---- D50 (Sesion 18 Etapa 2): print.gdpar_causal_bridge surfaces $warnings

test_that("print.gdpar_causal_bridge shows $warnings when populated", {
  mock_bridge <- structure(
    list(
      cate_draws = matrix(0.0, nrow = 4L, ncol = 3L),
      cate_mean  = c(0.0, 0.0, 0.0),
      cate_ci    = matrix(0.0, nrow = 3L, ncol = 2L,
                           dimnames = list(NULL, c("lower", "upper"))),
      newdata    = data.frame(x1 = seq_len(3L)),
      id_check   = list(treat = NULL, ctrl = NULL),
      fits       = list(
        treat = structure(list(family = list(name = "gaussian",
                                              link = "identity"),
                                  anchor = 0.0),
                           class = "gdpar_fit"),
        ctrl  = structure(list(family = list(name = "gaussian",
                                              link = "identity"),
                                  anchor = 0.0),
                           class = "gdpar_fit")
      ),
      type       = "response",
      level      = 0.95,
      n_draws    = 4L,
      n_obs      = 3L,
      call       = quote(gdpar_causal_bridge()),
      warnings   = "gdpar_causal_bridge trimmed posterior draws to S = 4 (treat had 10, ctrl had 4).",
      meta       = list(dim_kind = "scalar", dim_size = 1L,
                          dim_names = NULL,
                          newdata_source = "user_provided")
    ),
    class = c("gdpar_causal_bridge", "list")
  )
  expect_output(print(mock_bridge), "Warnings:", fixed = TRUE)
  expect_output(print(mock_bridge), "trimmed posterior draws", fixed = TRUE)
  expect_invisible(print(mock_bridge))
})

test_that("print.gdpar_causal_bridge omits Warnings: block when empty", {
  mock_bridge <- structure(
    list(
      cate_draws = matrix(0.0, nrow = 4L, ncol = 3L),
      cate_mean  = c(0.0, 0.0, 0.0),
      cate_ci    = matrix(0.0, nrow = 3L, ncol = 2L,
                           dimnames = list(NULL, c("lower", "upper"))),
      newdata    = data.frame(x1 = seq_len(3L)),
      id_check   = list(treat = NULL, ctrl = NULL),
      fits       = list(
        treat = structure(list(family = list(name = "gaussian",
                                              link = "identity"),
                                  anchor = 0.0),
                           class = "gdpar_fit"),
        ctrl  = structure(list(family = list(name = "gaussian",
                                              link = "identity"),
                                  anchor = 0.0),
                           class = "gdpar_fit")
      ),
      type       = "response",
      level      = 0.95,
      n_draws    = 4L,
      n_obs      = 3L,
      call       = quote(gdpar_causal_bridge()),
      warnings   = character(0L),
      meta       = list(dim_kind = "scalar", dim_size = 1L,
                          dim_names = NULL,
                          newdata_source = "user_provided")
    ),
    class = c("gdpar_causal_bridge", "list")
  )
  out <- capture.output(print(mock_bridge))
  expect_false(any(grepl("^Warnings:", out)))
})

test_that(".summarize_cate (matrix) produces vector mean and n x 2 CI", {
  set.seed(11L)
  draws <- matrix(rnorm(1000L * 3L, mean = 1.5, sd = 0.2),
                   nrow = 1000L, ncol = 3L)
  s <- gdpar:::.summarize_cate(draws, ql = 0.025, qu = 0.975)
  expect_equal(s$dim_kind, "scalar")
  expect_equal(length(s$mean), 3L)
  expect_equal(dim(s$ci), c(3L, 2L))
  expect_equal(colnames(s$ci), c("lower", "upper"))
  expect_true(all(s$ci[, 1L] < s$mean))
  expect_true(all(s$mean < s$ci[, 2L]))
  expect_lt(max(abs(s$mean - 1.5)), 0.05)
})

test_that(".summarize_cate (3D array) produces n x K mean and n x K x 2 CI", {
  set.seed(13L)
  draws <- array(rnorm(800L * 4L * 2L, mean = 0.5, sd = 0.3),
                  dim = c(800L, 4L, 2L),
                  dimnames = list(NULL, NULL, c("dim_1", "dim_2")))
  s <- gdpar:::.summarize_cate(draws, ql = 0.025, qu = 0.975)
  expect_true(s$dim_kind %in% c("multi", "K_individual"))
  expect_equal(dim(s$mean), c(4L, 2L))
  expect_equal(dim(s$ci),   c(4L, 2L, 2L))
  expect_equal(dimnames(s$ci)[[3L]], c("lower", "upper"))
  expect_lt(max(abs(s$mean - 0.5)), 0.10)
})

# ---- (iii) Opt-in end-to-end smoke + algorithmic equivalence (gated) ----

test_that("gdpar_causal_bridge end-to-end Gaussian T-learner (smoke + equivalence, opt-in)", {
  testthat::skip_if(Sys.getenv("GDPAR_RUN_STAN_SMOKE_CAUSAL_BRIDGE") != "1",
                    "GDPAR_RUN_STAN_SMOKE_CAUSAL_BRIDGE not set; skipping smoke")
  skip_if_no_cmdstan()
  testthat::skip_if_not_installed("posterior")
  set.seed(202605231L)
  n_per_arm <- 300L
  beta0 <- 0.2
  beta1 <- 0.8
  tau0  <- 1.0
  tau1  <- 0.5
  x1_t <- rnorm(n_per_arm)
  x1_c <- rnorm(n_per_arm)
  y_t  <- (beta0 + tau0) + (beta1 + tau1) * x1_t + rnorm(n_per_arm, sd = 0.4)
  y_c  <- beta0 + beta1 * x1_c + rnorm(n_per_arm, sd = 0.4)
  df_t <- data.frame(x1 = x1_t, y = y_t)
  df_c <- data.frame(x1 = x1_c, y = y_c)
  fit_t <- gdpar(
    formula       = y ~ x1,
    family        = gdpar_family("gaussian"),
    amm           = amm_spec(a = ~ x1),
    data          = df_t,
    chains        = 1L,
    iter_warmup   = 300L,
    iter_sampling = 300L,
    refresh       = 0L,
    verbose       = FALSE
  )
  fit_c <- gdpar(
    formula       = y ~ x1,
    family        = gdpar_family("gaussian"),
    amm           = amm_spec(a = ~ x1),
    data          = df_c,
    chains        = 1L,
    iter_warmup   = 300L,
    iter_sampling = 300L,
    refresh       = 0L,
    verbose       = FALSE
  )
  grid <- data.frame(x1 = seq(-2, 2, length.out = 11L))
  bridge <- gdpar_causal_bridge(fit_t, fit_c, newdata = grid)
  expect_s3_class(bridge, "gdpar_causal_bridge")
  expect_equal(bridge$n_obs, nrow(grid))
  expect_true(is.matrix(bridge$cate_draws))
  expect_equal(ncol(bridge$cate_draws), nrow(grid))
  expect_equal(bridge$type, "response")
  expect_equal(bridge$level, 0.95)
  # D48 (Sesion 17 Etapa 2 de Sesion 8.4): the canonical print method emits
  # a structured object summary by design; the original expect_silent
  # assumption was semantically incorrect. We validate (i) that print
  # returns invisibly the object and (ii) that the output carries the
  # canonical header.
  expect_output(print(bridge), "<gdpar_causal_bridge>", fixed = TRUE)
  expect_invisible(print(bridge))
  # D50 (Sesion 18 Etapa 2 de Sesion 8.4): the bridge object carries a
  # structured `$warnings` slot. Empty when both fits had equal draw
  # counts (the smoke configures both arms identically); populated
  # otherwise (the draw-trimming notification is persisted there).
  expect_true("warnings" %in% names(bridge))
  expect_type(bridge$warnings, "character")
  sm <- summary(bridge)
  expect_s3_class(sm, "summary.gdpar_causal_bridge")
  expect_true(is.data.frame(sm$table))
  expect_true(all(c("observation", "cate_mean", "cate_lower",
                    "cate_upper") %in% names(sm$table)))
  expect_output(print(sm), "<summary.gdpar_causal_bridge>", fixed = TRUE)
  expect_invisible(print(sm))
  # Algorithmic equivalence: the true CATE in the DGP is
  # tau(x) = tau0 + tau1 * x. The bridge estimate must be within
  # a permissive Monte Carlo tolerance per evaluation point.
  truth <- tau0 + tau1 * grid$x1
  expect_lt(mean(abs(bridge$cate_mean - truth)), 0.15)
  # predict.gdpar_causal_bridge on a fresh grid returns the same shape.
  grid2 <- data.frame(x1 = seq(-1, 1, length.out = 5L))
  re <- predict(bridge, newdata = grid2)
  expect_true(is.matrix(re$cate_draws))
  expect_equal(ncol(re$cate_draws), nrow(grid2))
  expect_equal(length(re$cate_mean), nrow(grid2))
})

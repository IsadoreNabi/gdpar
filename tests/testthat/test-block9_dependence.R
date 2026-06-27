# Block 9, Sub-block 9.1, Session B9.8 (D75): dependence-robust inference
# (Axis 2). Layers 1-2 (resampler algebra + API guards) run by default and
# are Stan-free. Layer 3 (end-to-end fit-based behaviour) is gated behind
# GDPAR_RUN_BLOCK9_DEP_FITS=1 because each assertion refits the EB model.

# ---------------------------------------------------------------------------
# Layer 1: block-bootstrap resampler algebra (no Stan)
# ---------------------------------------------------------------------------

test_that("moving block resampler returns a valid length-n index vector", {
  set.seed(1)
  idx <- .gdpar_block_bootstrap_data_indices(50L, 5L, "moving")
  expect_length(idx, 50L)
  expect_true(all(idx >= 1L & idx <= 50L))
  expect_true(is.integer(idx))
})

test_that("moving block with block_length = n is the identity order", {
  # n_starts = 1, so the single contiguous block is deterministically 1:n.
  expect_identical(
    .gdpar_block_bootstrap_data_indices(20L, 20L, "moving"),
    seq_len(20L)
  )
})

test_that("circular block resampler wraps and stays in range", {
  set.seed(2)
  idx <- .gdpar_block_bootstrap_data_indices(30L, 7L, "circular")
  expect_length(idx, 30L)
  expect_true(all(idx >= 1L & idx <= 30L))
})

test_that("circular block with block_length = n is a single rotation", {
  set.seed(3)
  idx <- .gdpar_block_bootstrap_data_indices(12L, 12L, "circular")
  expect_setequal(idx, seq_len(12L))  # a cyclic permutation: all distinct
})

test_that("block resampler is deterministic under a fixed seed", {
  set.seed(99)
  a <- .gdpar_block_bootstrap_data_indices(40L, 4L, "moving")
  set.seed(99)
  b <- .gdpar_block_bootstrap_data_indices(40L, 4L, "moving")
  expect_identical(a, b)
})

test_that("GOLDEN: block resampler indices are frozen bit-exact (Stan-free)", {
  # Deterministic RNG golden (D76): locks the resampling algorithm itself,
  # including contiguity of blocks and the circular wrap. The HMC-dependent
  # bootstrap standard errors are deliberately NOT frozen bit-exact, because
  # HMC is not reproducible across Stan toolchains; freezing them would be a
  # fragile, masking golden contrary to the no-check-masking regime.
  set.seed(20260528L)
  expect_identical(
    .gdpar_block_bootstrap_data_indices(20L, 4L, "moving"),
    c(14L, 15L, 16L, 17L, 16L, 17L, 18L, 19L, 12L, 13L, 14L, 15L,
      11L, 12L, 13L, 14L, 8L, 9L, 10L, 11L)
  )
  set.seed(20260528L)
  expect_identical(
    .gdpar_block_bootstrap_data_indices(16L, 4L, "circular"),
    c(2L, 3L, 4L, 5L, 14L, 15L, 16L, 1L, 16L, 1L, 2L, 3L, 13L, 14L,
      15L, 16L)
  )
})

test_that("block resampler rejects out-of-range block lengths", {
  expect_error(
    .gdpar_block_bootstrap_data_indices(10L, 0L, "moving"),
    class = "gdpar_input_error"
  )
  expect_error(
    .gdpar_block_bootstrap_data_indices(10L, 11L, "moving"),
    class = "gdpar_input_error"
  )
})

test_that("default block length follows the n^(1/3) rate", {
  expect_identical(.gdpar_default_block_length(27L), 3L)
  expect_identical(.gdpar_default_block_length(1000L), 10L)
  expect_true(.gdpar_default_block_length(1L) >= 1L)
})

# ---------------------------------------------------------------------------
# Layer 1b: data-driven block length, Politis-White (D101; no Stan)
# ---------------------------------------------------------------------------

test_that("flat-top lag window has the exact Politis (2003) shape", {
  expect_equal(.gdpar_flat_top_window(c(0, 0.5, 0.75, 1, 1.5, -0.75)),
               c(1, 1, 0.5, 0, 0, 0.5))
})

test_that("Politis (2003) bandwidth picks the first run of Kn negligible lags", {
  Kn <- 5L; crit <- 0.1
  # Significant lags 1..3, then 5 negligible -> run starts at j = 4.
  rho <- c(0.5, 0.4, 0.3, rep(0.01, 5L), 0.2)
  expect_identical(.gdpar_pw_mhat(rho, Kn, crit), 4L)
  # Every lag negligible -> 1; no negligible run -> largest significant lag.
  expect_identical(.gdpar_pw_mhat(rep(0.01, 10L), Kn, crit), 1L)
  expect_identical(.gdpar_pw_mhat(rep(0.5, 10L), Kn, crit), 10L)
})

test_that("is_auto recognises only the literal string \"auto\"", {
  expect_true(.gdpar_is_auto("auto"))
  expect_false(.gdpar_is_auto(NULL))
  expect_false(.gdpar_is_auto(5L))
  expect_false(.gdpar_is_auto(NA_character_))
  expect_false(.gdpar_is_auto(c("auto", "auto")))
  expect_false(.gdpar_is_auto("AUTO"))
})

test_that("NUMERICAL: Politis-White matches an independent reference (no Stan)", {
  # An independent, loop-based re-implementation of the overlapping (circular-
  # block) Politis-White / Patton-Politis-White selector. Agreement of two
  # independent code paths is the numerical-layer check; no external 'np'
  # dependency is incurred and the result is deterministic given the input.
  pw_ref <- function(resid, c_thresh = stats::qnorm(0.975)) {
    n <- length(resid)
    Kn <- max(5L, as.integer(ceiling(log10(n))))
    Mmax <- min(as.integer(ceiling(sqrt(n)) + Kn), n - 1L)
    r <- resid - mean(resid); v <- sum(r^2) / n
    rho <- vapply(seq_len(Mmax),
                  function(k) sum(r[1:(n - k)] * r[(k + 1):n]) / n / v, numeric(1))
    crit <- c_thresh * sqrt(log10(n) / n)
    insig <- abs(rho) < crit
    mhat <- NA_integer_
    for (j in seq_len(Mmax - Kn + 1L)) {
      if (all(insig[j:(j + Kn - 1L)])) { mhat <- j; break }
    }
    if (is.na(mhat)) {
      sig <- which(!insig)
      mhat <- if (length(sig)) sig[length(sig)] else 1L
    }
    M <- min(2L * mhat, Mmax)
    Rk <- vapply(0:M,
                 function(k) sum(r[1:(n - k)] * r[(k + 1):n]) / n, numeric(1))
    lam <- .gdpar_flat_top_window((0:M) / M)
    spec <- Rk[1] * lam[1] + 2 * sum(lam[-1] * Rk[-1])
    ghat <- 2 * sum(lam[-1] * (1:M) * Rk[-1])
    b <- (2 * ghat^2 / ((4 / 3) * spec^2))^(1 / 3) * n^(1 / 3)
    Bmax <- ceiling(min(3 * sqrt(n), n / 3))
    max(1L, min(as.integer(round(b)), as.integer(Bmax)))
  }
  set.seed(2026L)
  for (phi in c(0.2, 0.6, -0.4)) {
    x <- as.numeric(stats::arima.sim(list(ar = phi), 200L))
    expect_identical(.gdpar_politis_white_block_length(x)$block_length,
                     pw_ref(x))
  }
})

test_that("Politis-White: longer block for stronger dependence; 1 for white noise", {
  set.seed(7L)
  med_block <- function(phi) {
    stats::median(vapply(seq_len(15L), function(i) {
      .gdpar_politis_white_block_length(
        as.numeric(stats::arima.sim(list(ar = phi), 300L)))$block_length
    }, numeric(1)))
  }
  expect_gt(med_block(0.8), med_block(0.2))
  set.seed(11L)
  wn <- .gdpar_politis_white_block_length(rnorm(400L))
  expect_identical(wn$method, "auto")
  expect_identical(wn$block_length, 1L)        # no dependence -> unit blocks
})

test_that("Politis-White falls back to the rate on a degenerate series", {
  expect_identical(.gdpar_politis_white_block_length(rnorm(5L))$method, "rate")
  expect_identical(.gdpar_politis_white_block_length(rep(3, 50L))$method, "rate")
})

# ---------------------------------------------------------------------------
# Layer 2: public API guards (no Stan; fail before any refit)
# ---------------------------------------------------------------------------

test_that("dependence functions validate the supplied fit (EB or FB, D102)", {
  # A non-fit object is rejected outright.
  expect_error(
    gdpar_dependence_diagnostic(list(1, 2, 3)),
    class = "gdpar_input_error"
  )
  # A scalar full-Bayes gdpar_fit is now accepted as a class (D102 closes the
  # EB/FB asymmetry), but a malformed one lacking its HMC fit is caught.
  no_fit <- structure(list(), class = c("gdpar_fit", "list"))
  expect_error(
    gdpar_dependence_diagnostic(no_fit),
    class = "gdpar_internal_error"
  )
  # Multivariate / K-individual full-Bayes fits remain deferred on both paths.
  multi_fb <- structure(list(fit = list(), p = 2L),
                        class = c("gdpar_fit", "list"))
  expect_error(
    gdpar_dependence_robust(multi_fb, data = data.frame(x = 1)),
    class = "gdpar_unsupported_feature_error"
  )
  k_fb <- structure(list(fit = list(), K = 2L),
                    class = c("gdpar_fit", "list"))
  expect_error(
    gdpar_dependence_diagnostic(k_fb),
    class = "gdpar_unsupported_feature_error"
  )
})

test_that("gdpar_dependence_robust rejects a non-data-frame data argument", {
  fake_eb <- structure(
    list(conditional_fit = list(), family = list(name = "gaussian")),
    class = c("gdpar_eb_fit", "list")
  )
  expect_error(
    gdpar_dependence_robust(fake_eb, data = matrix(1:4, 2)),
    class = "gdpar_input_error"
  )
})

test_that("gdpar_dependence_robust rejects a non-\"auto\" string block_length", {
  fake_eb <- structure(
    list(conditional_fit = list(), family = list(name = "gaussian")),
    class = c("gdpar_eb_fit", "list")
  )
  # Caught Stan-free, before any require_suggested / refit.
  expect_error(
    gdpar_dependence_robust(fake_eb, data = data.frame(x = 1, y = 1),
                            block_length = "bogus"),
    class = "gdpar_input_error"
  )
})

# ---------------------------------------------------------------------------
# Layer 2b: full-Bayes (gdpar_fit) touchpoints, Stan-free (D102)
#
# The FB extension dispatches three EB-specific touchpoints by object class.
# These tests exercise the dispatch + the FB extractors against a mock
# gdpar_fit carrying synthetic posterior draws, with no Stan involved.
# ---------------------------------------------------------------------------

# Build a mock gdpar_fit whose $fit$draws() returns a synthetic draws object
# over the requested coefficient columns, plus the slots the touchpoints read.
.make_fake_fb_fit <- function(coef_mat, amm = list(a = TRUE),
                              y_pred_mat = NULL, y_obs = NULL,
                              stan_id = 1L,
                              diagnostics = list(rhat_max = 1.0,
                                                 ess_bulk_min = 500,
                                                 divergent_count = 0L)) {
  cols <- coef_mat
  if (!is.null(y_pred_mat)) {
    colnames(y_pred_mat) <- paste0("y_pred[", seq_len(ncol(y_pred_mat)), "]")
    cols <- cbind(cols, y_pred_mat)
  }
  dr <- posterior::as_draws_array(posterior::as_draws_matrix(cols))
  stan_data <- list()
  if (!is.null(y_obs)) {
    if (stan_id %in% c(2L, 3L, 4L)) stan_data$y_int <- y_obs else
      stan_data$y_real <- y_obs
  }
  structure(
    list(
      fit         = list(draws = function(...) {
        dots <- list(...)
        out <- dr
        if (!is.null(dots$variables)) {
          out <- posterior::subset_draws(out, variable = dots$variables)
        }
        if (identical(dots$format, "draws_matrix")) {
          out <- posterior::as_draws_matrix(out)
        }
        out
      }),
      amm         = amm,
      family      = list(name = "gaussian", stan_id = stan_id),
      stan_data   = stan_data,
      diagnostics = diagnostics
    ),
    class = c("gdpar_fit", "list")
  )
}

test_that("the scalar gate accepts a well-formed FB fit and defers K>1/p>1", {
  skip_if_not_installed("posterior")
  ok_fb <- structure(list(fit = list()), class = c("gdpar_fit", "list"))
  expect_invisible(gdpar:::.gdpar_assert_scalar_dep(ok_fb))
  # No HMC fit -> internal error.
  expect_error(
    gdpar:::.gdpar_assert_scalar_dep(
      structure(list(), class = c("gdpar_fit", "list"))),
    class = "gdpar_internal_error"
  )
  # Multivariate / K-individual deferred.
  expect_error(
    gdpar:::.gdpar_assert_scalar_dep(
      structure(list(fit = list(), p = 2L), class = c("gdpar_fit", "list"))),
    class = "gdpar_unsupported_feature_error"
  )
  expect_error(
    gdpar:::.gdpar_assert_scalar_dep(
      structure(list(fit = list(), K = 3L), class = c("gdpar_fit", "list"))),
    class = "gdpar_unsupported_feature_error"
  )
})

test_that("FB estimate/SE extractors are posterior mean / SD, name-aligned", {
  skip_if_not_installed("posterior")
  set.seed(11)
  S <- 300L
  m <- cbind(
    theta_ref   = stats::rnorm(S, 2.0, 0.50),
    `a_coef[1]` = stats::rnorm(S, -1.0, 0.30),
    `a_coef[2]` = stats::rnorm(S, 0.7, 0.20)
  )
  fb <- .make_fake_fb_fit(m, amm = list(a = TRUE))
  est <- gdpar:::.gdpar_dep_estimate_vector(fb)
  se  <- gdpar:::.gdpar_dep_model_se_vector(fb)
  expect_equal(est, colMeans(m), tolerance = 1e-10)
  expect_equal(se, apply(m, 2L, stats::sd), tolerance = 1e-10)
  expect_identical(names(est), c("theta_ref", "a_coef[1]", "a_coef[2]"))
  # Hyperparameters are excluded for EB/FB parity: add a mu_theta_ref column and
  # confirm it does NOT enter the coefficient vector.
  m2 <- cbind(m, mu_theta_ref = stats::rnorm(S))
  fb2 <- .make_fake_fb_fit(m2, amm = list(a = TRUE))
  expect_false("mu_theta_ref" %in% names(gdpar:::.gdpar_dep_estimate_vector(fb2)))
})

test_that("FB residual touchpoint delegates to residuals.gdpar_fit", {
  skip_if_not_installed("posterior")
  set.seed(12)
  n <- 8L
  S <- 50L
  y_obs <- stats::rnorm(n, 3, 1)
  y_pred <- matrix(stats::rnorm(S * n, rep(y_obs, each = S), 0.5), S, n)
  m <- cbind(theta_ref = stats::rnorm(S))
  fb <- .make_fake_fb_fit(m, amm = list(), y_pred_mat = y_pred,
                          y_obs = y_obs, stan_id = 1L)
  r_touch <- gdpar:::.gdpar_dependence_residuals(fb, "response", NULL)
  r_meth  <- as.numeric(residuals(fb, type = "response"))
  expect_equal(r_touch, r_meth, tolerance = 1e-12)
  # Response residual is exactly y_obs - colMeans(y_pred).
  expect_equal(r_touch, y_obs - colMeans(y_pred), tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# Layer 3: end-to-end behaviour (gated; refits the EB model)
# ---------------------------------------------------------------------------

test_that("diagnostic detects AR(1) dependence and not iid noise", {
  skip_if(Sys.getenv("GDPAR_RUN_BLOCK9_DEP_FITS") != "1",
          "set GDPAR_RUN_BLOCK9_DEP_FITS=1 to run EB fit-based tests")
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("posterior")
  set.seed(440)
  n <- 120L
  x <- rnorm(n)
  y_ar  <- 1 + 0.5 * x + as.numeric(stats::arima.sim(list(ar = 0.7), n))
  y_iid <- 1 + 0.5 * x + rnorm(n)
  df_ar  <- data.frame(x = x, y = y_ar,  t = seq_len(n))
  df_iid <- data.frame(x = x, y = y_iid, t = seq_len(n))
  fit_ar <- gdpar_eb(y ~ x, amm = amm_spec(a = ~ x), data = df_ar,
                     chains = 2, iter_warmup = 150, iter_sampling = 150,
                     verbose = FALSE, refresh = 0, seed = 41)
  fit_iid <- gdpar_eb(y ~ x, amm = amm_spec(a = ~ x), data = df_iid,
                      chains = 2, iter_warmup = 150, iter_sampling = 150,
                      verbose = FALSE, refresh = 0, seed = 42)
  d_ar <- gdpar_dependence_diagnostic(fit_ar, index = df_ar$t,
                                      residual_type = "response")
  d_iid <- gdpar_dependence_diagnostic(fit_iid, index = df_iid$t,
                                       residual_type = "response")
  expect_s3_class(d_ar, "gdpar_dependence_diagnostic")
  expect_lt(d_ar$ljung_box_p_value, 0.05)
  expect_gt(d_ar$lag1_autocorr, 0.3)
  expect_gt(d_iid$ljung_box_p_value, 0.05)
})

test_that("block bootstrap returns a robust-SE table and is seed-deterministic", {
  skip_if(Sys.getenv("GDPAR_RUN_BLOCK9_DEP_FITS") != "1",
          "set GDPAR_RUN_BLOCK9_DEP_FITS=1 to run EB fit-based tests")
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("posterior")
  set.seed(450)
  n <- 100L
  x <- rnorm(n)
  y <- 1 + 0.5 * x + as.numeric(stats::arima.sim(list(ar = 0.6), n))
  df <- data.frame(x = x, y = y, t = seq_len(n))
  fit <- gdpar_eb(y ~ x, amm = amm_spec(a = ~ x), data = df,
                  chains = 2, iter_warmup = 150, iter_sampling = 150,
                  verbose = FALSE, refresh = 0, seed = 45)
  rb <- gdpar_dependence_robust(fit, data = df, index = df$t, B = 15,
                                seed = 7, verbose = FALSE,
                                iter_warmup = 150, iter_sampling = 150,
                                chains = 2)
  expect_s3_class(rb, "gdpar_dependence_robust")
  expect_true(all(c("estimate", "model_se", "robust_se", "se_ratio",
                    "ci_lower", "ci_upper") %in% names(rb$table)))
  expect_true(all(is.finite(rb$table$robust_se)))
  expect_gte(rb$B_ok, 2L)
  expect_true(all(rb$table$ci_lower <= rb$table$ci_upper))
  rb2 <- gdpar_dependence_robust(fit, data = df, index = df$t, B = 15,
                                 seed = 7, verbose = FALSE,
                                 iter_warmup = 150, iter_sampling = 150,
                                 chains = 2)
  expect_equal(rb$table$robust_se, rb2$table$robust_se)
})

test_that("block_length = \"auto\" selects a data-driven length end-to-end", {
  skip_if(Sys.getenv("GDPAR_RUN_BLOCK9_DEP_FITS") != "1",
          "set GDPAR_RUN_BLOCK9_DEP_FITS=1 to run EB fit-based tests")
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("posterior")
  set.seed(451)
  n <- 120L
  x <- rnorm(n)
  y <- 1 + 0.5 * x + as.numeric(stats::arima.sim(list(ar = 0.6), n))
  df <- data.frame(x = x, y = y, t = seq_len(n))
  fit <- gdpar_eb(y ~ x, amm = amm_spec(a = ~ x), data = df,
                  chains = 2, iter_warmup = 150, iter_sampling = 150,
                  verbose = FALSE, refresh = 0, seed = 46)
  rb <- gdpar_dependence_robust(fit, data = df, index = df$t,
                                block_length = "auto", B = 12, seed = 7,
                                verbose = FALSE, iter_warmup = 150,
                                iter_sampling = 150, chains = 2)
  expect_s3_class(rb, "gdpar_dependence_robust")
  expect_true(rb$block_length_method %in% c("auto", "rate"))
  expect_gte(rb$block_length, 1L)
  expect_lte(rb$block_length, n)
  expect_true(all(is.finite(rb$table$robust_se)))
  # The auto length is a deterministic function of the (seeded) residuals.
  rb2 <- gdpar_dependence_robust(fit, data = df, index = df$t,
                                 block_length = "auto", B = 12, seed = 7,
                                 verbose = FALSE, iter_warmup = 150,
                                 iter_sampling = 150, chains = 2)
  expect_identical(rb$block_length, rb2$block_length)
})

# ---------------------------------------------------------------------------
# Layer 3b: full-Bayes (gdpar_fit) end-to-end (gated; refits the FB model)
# ---------------------------------------------------------------------------

test_that("FB path: diagnostic + temporal robust SE end-to-end (D102)", {
  skip_if(Sys.getenv("GDPAR_RUN_BLOCK9_DEP_FITS") != "1",
          "set GDPAR_RUN_BLOCK9_DEP_FITS=1 to run fit-based tests")
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("posterior")
  set.seed(460)
  n <- 90L
  x <- rnorm(n)
  y <- 1 + 0.5 * x + as.numeric(stats::arima.sim(list(ar = 0.6), n))
  df <- data.frame(x = x, y = y, t = seq_len(n))
  fb <- gdpar(y ~ x, amm = amm_spec(a = ~ x), data = df,
              chains = 2, iter_warmup = 150, iter_sampling = 150,
              verbose = FALSE, refresh = 0, seed = 46)
  expect_s3_class(fb, "gdpar_fit")
  # The same diagnostic accepts a gdpar_fit (D102).
  d <- gdpar_dependence_diagnostic(fb, index = df$t, residual_type = "response")
  expect_s3_class(d, "gdpar_dependence_diagnostic")
  expect_lt(d$ljung_box_p_value, 0.05)
  # Robust SE by full-Bayes refit (full HMC). B kept tiny for the gated suite.
  rb <- gdpar_dependence_robust(fb, data = df, index = df$t, B = 4, seed = 7,
                                verbose = FALSE, iter_warmup = 150,
                                iter_sampling = 150, chains = 2)
  expect_s3_class(rb, "gdpar_dependence_robust")
  expect_true(all(c("estimate", "model_se", "robust_se", "se_ratio",
                    "ci_lower", "ci_upper") %in% names(rb$table)))
  expect_true(all(is.finite(rb$table$robust_se)))
  expect_true(all(is.finite(rb$table$model_se)))
  expect_gte(rb$B_ok, 2L)
  expect_true(all(rb$table$ci_lower <= rb$table$ci_upper))
  # The FB-specific refit diagnostics are populated.
  expect_false(is.null(rb$refit_diagnostics))
  expect_true(is.finite(rb$refit_diagnostics$max_rhat))
  # Same-session seed reproducibility.
  rb2 <- gdpar_dependence_robust(fb, data = df, index = df$t, B = 4, seed = 7,
                                 verbose = FALSE, iter_warmup = 150,
                                 iter_sampling = 150, chains = 2)
  expect_equal(rb$table$robust_se, rb2$table$robust_se)
})

test_that("FB path: temporal block_length = \"auto\" end-to-end (D102)", {
  skip_if(Sys.getenv("GDPAR_RUN_BLOCK9_DEP_FITS") != "1",
          "set GDPAR_RUN_BLOCK9_DEP_FITS=1 to run fit-based tests")
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("posterior")
  set.seed(461)
  n <- 90L
  x <- rnorm(n)
  y <- 1 + 0.5 * x + as.numeric(stats::arima.sim(list(ar = 0.6), n))
  df <- data.frame(x = x, y = y, t = seq_len(n))
  fb <- gdpar(y ~ x, amm = amm_spec(a = ~ x), data = df,
              chains = 2, iter_warmup = 150, iter_sampling = 150,
              verbose = FALSE, refresh = 0, seed = 47)
  # The Politis-White selector runs on FB Dunn-Smyth residuals (path-agnostic).
  rb <- gdpar_dependence_robust(fb, data = df, index = df$t,
                                block_length = "auto", B = 4, seed = 7,
                                verbose = FALSE, iter_warmup = 150,
                                iter_sampling = 150, chains = 2)
  expect_true(rb$block_length_method %in% c("auto", "rate"))
  expect_gte(rb$block_length, 1L)
  expect_true(all(is.finite(rb$table$robust_se)))
})

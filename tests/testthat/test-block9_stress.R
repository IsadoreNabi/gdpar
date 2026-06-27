# Block 9.1 -- Internal synthetic adversarial re-validation.
# Layer 3 (Decision C = C.iv, Session B9.7, 2026-05-28): STRESS tests for
# extreme outliers in the outcome and ROBUSTNESS to NA / missing data.
#
#   Sections A-B  -- run by default. NA / missing handling (clean aborts)
#                    and extreme-but-finite covariate outliers (numerical
#                    stability of the identifiability algebra). No Stan
#                    compilation; the NA-in-outcome aborts are guarded by
#                    skip_if_not_installed("cmdstanr") because they enter
#                    gdpar_eb().
#   Section C     -- GATED extreme-outlier OUTCOME stress over each of the
#                    12 canonical families. Substantive execution is the
#                    nocturnal B9.8 sub-unit under Decision C.iv (each
#                    family requires a real fit). The contract is: a fit
#                    on outlier-laden, in-support data must either return
#                    FINITE estimates or abort with a CANONICAL error --
#                    never silently produce NaN/Inf estimates.

# ---- Section A: NA / missing robustness ---------------------------------

test_that("STR-A1: NA in an AMM covariate aborts with gdpar_input_error (pre-fit)", {
  df <- make_gaussian_data(n = 80, seed = 2)
  df$x1[5] <- NA
  expect_error(
    gdpar_check_identifiability(amm_spec(a = ~ x1), df, formula_rhs = ~ x1),
    class = "gdpar_input_error"
  )
})

test_that("STR-A2: NA in the outcome aborts cleanly with gdpar_input_error", {
  skip_if_not_installed("cmdstanr")
  df <- make_gaussian_data(n = 80, seed = 3)
  df$y[7] <- NA
  expect_error(
    suppressMessages(gdpar_eb(y ~ x1, data = df, family = gdpar_family("gaussian"))),
    class = "gdpar_input_error"
  )
})

test_that("STR-A3: NaN in the outcome is treated as missing and aborts cleanly", {
  skip_if_not_installed("cmdstanr")
  df <- make_gaussian_data(n = 80, seed = 4)
  df$y[9] <- NaN
  expect_error(
    suppressMessages(gdpar_eb(y ~ x1, data = df, family = gdpar_family("gaussian"))),
    class = "gdpar_input_error"
  )
})

test_that("STR-A4: Inf in the outcome aborts pre-fit with gdpar_input_error (D74)", {
  # Before D74 an infinite outcome passed the is.na() guard (is.na(Inf) is
  # FALSE) and surfaced late as an opaque gdpar_eb_numerical_error. The
  # broadened !is.finite() guard now rejects it cleanly at the door.
  skip_if_not_installed("cmdstanr")
  df <- make_gaussian_data(n = 80, seed = 5)
  df$y[11] <- Inf
  expect_error(
    suppressMessages(gdpar_eb(y ~ x1, data = df, family = gdpar_family("gaussian"))),
    class = "gdpar_input_error"
  )
})

# ---- Section B: extreme-but-finite covariate outliers -------------------

test_that("STR-B1: 99.9% / 0.1% quantile covariate outliers keep the conditioning finite", {
  set.seed(515)
  n <- 400
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  # Inject symmetric extreme outliers at the tails.
  x1[1] <- 50 * stats::quantile(abs(x1), 0.999)
  x1[2] <- -50 * stats::quantile(abs(x1), 0.999)
  df <- data.frame(x1 = x1, x2 = x2, y = rnorm(n))
  rep <- gdpar_check_identifiability(amm_spec(a = ~ x1 + x2), df,
                                     formula_rhs = ~ x1 + x2)
  expect_true(is.finite(rep$condition_number))
  expect_true(rep$passed)
})

test_that("STR-B2: a covariate on a 1e6 scale does not destroy the conditioning (normalization)", {
  set.seed(616)
  n <- 300
  df <- data.frame(x1 = 1e6 * rnorm(n), x2 = rnorm(n), y = rnorm(n))
  rep <- gdpar_check_identifiability(amm_spec(a = ~ x1 + x2), df,
                                     formula_rhs = ~ x1 + x2)
  expect_true(is.finite(rep$condition_number))
  expect_true(rep$passed)
})

# ---- Section C: outcome-outlier stress per family (GATED; Stan-bound) ----
# Substantive execution is the nocturnal B9.8 sub-unit under Decision C.iv.
# Enable with GDPAR_RUN_BLOCK9_STRESS_FITS=1 and a configured cmdstan.

.skip_if_no_block9_stress_fits <- function() {
  if (!nzchar(Sys.getenv("GDPAR_RUN_BLOCK9_STRESS_FITS"))) {
    skip("GDPAR_RUN_BLOCK9_STRESS_FITS not set; skipping Block 9.1 stress fits")
  }
  skip_if_not_installed("cmdstanr")
}

# Support-aware synthetic data with extreme tail outliers in y_obs, kept
# strictly inside the family's support. Returns a data.frame with x1, x2, y.
.block9_make_outlier_data <- function(family, n = 250, seed = 1L) {
  set.seed(seed)
  x1 <- rnorm(n); x2 <- rnorm(n)
  eta <- 0.3 * x1 - 0.2 * x2
  lo <- max(1L, floor(0.001 * n)); hi <- seq_len(lo)
  hi2 <- (n - lo + 1L):n
  y <- switch(
    family,
    gaussian              = { v <- 1 + eta + rnorm(n, sd = 0.5); v[hi] <- 50; v[hi2] <- -50; v },
    student_t             = { v <- 1 + eta + rt(n, df = 3);      v[hi] <- 80; v[hi2] <- -80; v },
    poisson               = { v <- rpois(n, exp(1 + eta));       v[hi] <- 1e4L; v },
    zip                   = { v <- rpois(n, exp(1 + eta));       v[hi2] <- 0L; v[hi] <- 1e4L; v },
    hurdle_poisson        = { v <- rpois(n, exp(1 + eta)) + 1L;  v[hi] <- 1e4L; v },
    neg_binomial_2        = { v <- rnbinom(n, mu = exp(1 + eta), size = 5); v[hi] <- 1e5L; v },
    zinb                  = { v <- rnbinom(n, mu = exp(1 + eta), size = 5); v[hi2] <- 0L; v[hi] <- 1e5L; v },
    hurdle_neg_binomial_2 = { v <- rnbinom(n, mu = exp(1 + eta), size = 5) + 1L; v[hi] <- 1e5L; v },
    gamma                 = { v <- rgamma(n, shape = 2, rate = exp(-(1 + eta))); v[hi] <- 1e6; v[hi2] <- 1e-6; v },
    tweedie               = { v <- rgamma(n, shape = 2, rate = exp(-(1 + eta))); v[hi] <- 1e6; v[hi2] <- 0; v },
    beta                  = { v <- rbeta(n, 2, 2); v[hi] <- 1 - 1e-6; v[hi2] <- 1e-6; v },
    bernoulli             = { v <- rbinom(n, 1, plogis(eta)); v[] <- 0L; v[hi] <- 1L; v },
    stop("unknown family in .block9_make_outlier_data")
  )
  data.frame(x1 = x1, x2 = x2, y = y)
}

.BLOCK9_CANONICAL_FAMILIES <- c(
  "gaussian", "poisson", "neg_binomial_2", "bernoulli", "beta", "gamma",
  "student_t", "tweedie", "zip", "zinb", "hurdle_poisson",
  "hurdle_neg_binomial_2"
)

test_that("STR-C (gated): outcome outliers yield finite estimates or a canonical abort, never silent NaN -- all 12 families", {
  .skip_if_no_block9_stress_fits()
  for (fam_name in .BLOCK9_CANONICAL_FAMILIES) {
    df <- .block9_make_outlier_data(fam_name, n = 200, seed = 909L)
    res <- tryCatch(
      suppressMessages(gdpar_eb(y ~ x1, data = df,
                                family = gdpar_family(fam_name), seed = 1L)),
      error = function(e) e
    )
    if (inherits(res, "error")) {
      # A clean canonical abort is an acceptable outcome.
      expect_true(
        inherits(res, c("gdpar_input_error", "gdpar_eb_numerical_error",
                        "gdpar_unsupported_feature_error")),
        info = sprintf("%s: non-canonical error class %s",
                       fam_name, class(res)[1])
      )
    } else {
      # Otherwise the point estimate must be finite (no silent NaN/Inf).
      est <- as.numeric(res$theta_ref_hat)
      expect_true(all(is.finite(est)),
                  info = sprintf("%s: non-finite theta_ref estimate", fam_name))
    }
  }
})

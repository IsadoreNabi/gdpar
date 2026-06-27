# Block 9.1 -- Internal synthetic adversarial re-validation.
# Layer 2 (Decision C = C.iv, Session B9.7, 2026-05-28): PROPERTY-BASED
# invariances. A correct estimator must respect the equivariances of the
# model under group actions on the data. We test them at two levels:
#
#   Sections A-B  -- DETERMINISTIC algebraic invariances of the
#                    identifiability eigenstructure (row permutation,
#                    covariate rescaling). These hold *exactly* (up to
#                    floating-point summation order) and run by default;
#                    no Stan compilation.
#   Section C     -- STATISTICAL invariances of the fitted EB estimates
#                    (permutation invariance of the likelihood, continuity
#                    under small noise injection). These require Stan and
#                    are gated by GDPAR_RUN_BLOCK9_FIT_FUZZ; substantive
#                    execution is the nocturnal B9.8 sub-unit under C.iv.

# ---- Section A: permutation invariance (deterministic) ------------------

test_that("INV-A1: row permutation leaves the identifiability eigenstructure invariant", {
  df <- make_gaussian_data(n = 250, seed = 42)
  spec <- amm_spec(a = ~ x1 + x2)
  rep0 <- gdpar_check_identifiability(spec, df, formula_rhs = ~ x1 + x2)

  set.seed(99)
  perm <- sample(nrow(df))
  rep1 <- gdpar_check_identifiability(spec, df[perm, , drop = FALSE],
                                      formula_rhs = ~ x1 + x2)

  expect_identical(rep0$passed, rep1$passed)
  # The Gram matrix Z'Z is a sum over rows, hence permutation-invariant;
  # eigenvalues match up to floating-point summation order.
  expect_equal(rep1$condition_number, rep0$condition_number,
               tolerance = 1e-6)
  expect_equal(rep1$lambda_min, rep0$lambda_min, tolerance = 1e-8)
  expect_equal(rep1$lambda_max, rep0$lambda_max, tolerance = 1e-8)
})

test_that("INV-A2: per-coord (C4-bis) verdict is permutation-invariant for p > 1", {
  set.seed(11)
  n <- 220
  df <- data.frame(x1 = rnorm(n), x2 = rnorm(n),
                   z1 = rnorm(n), z2 = rnorm(n))
  spec <- amm_spec(
    p = 2L,
    dims = dimwise(a = ~ x1 + x2, b = NULL),
    W = W_basis(type = "polynomial", degree = 2),
    x_vars = c("z1", "z2")
  )
  args <- list(amm = spec, theta_ref_init = c(0.4, -0.3),
               formula_rhs = ~ x1 + x2, rigor = "full")
  rep0 <- do.call(gdpar_check_identifiability, c(list(data = df), args))
  set.seed(123)
  perm <- sample(n)
  rep1 <- do.call(gdpar_check_identifiability,
                  c(list(data = df[perm, , drop = FALSE]), args))
  expect_identical(rep0$passed, rep1$passed)
  expect_identical(length(rep0$c4_bis$per_k), length(rep1$c4_bis$per_k))
  for (k in seq_along(rep0$c4_bis$per_k)) {
    expect_identical(rep0$c4_bis$per_k[[k]]$passed,
                     rep1$c4_bis$per_k[[k]]$passed)
  }
})

# ---- Section B: scale invariance (deterministic) ------------------------

test_that("INV-B1: rescaling a covariate leaves the identifiability verdict invariant", {
  df <- make_gaussian_data(n = 250, seed = 7)
  spec <- amm_spec(a = ~ x1 + x2)
  rep0 <- gdpar_check_identifiability(spec, df, formula_rhs = ~ x1 + x2)

  df_scaled <- df
  df_scaled$x2 <- 1e4 * df_scaled$x2
  rep1 <- gdpar_check_identifiability(spec, df_scaled, formula_rhs = ~ x1 + x2)

  # The check normalizes columns to unit norm, so a positive rescaling of
  # any single covariate must not change pass/fail nor the conditioning.
  expect_identical(rep0$passed, rep1$passed)
  expect_equal(rep1$condition_number, rep0$condition_number,
               tolerance = 1e-4)
})

test_that("INV-B2: sign flip of a covariate leaves the eigenstructure invariant", {
  df <- make_gaussian_data(n = 250, seed = 21)
  spec <- amm_spec(a = ~ x1 + x2)
  rep0 <- gdpar_check_identifiability(spec, df, formula_rhs = ~ x1 + x2)
  df_flip <- df
  df_flip$x1 <- -df_flip$x1
  rep1 <- gdpar_check_identifiability(spec, df_flip, formula_rhs = ~ x1 + x2)
  expect_identical(rep0$passed, rep1$passed)
  expect_equal(rep1$condition_number, rep0$condition_number,
               tolerance = 1e-6)
})

# ---- Section C: fitted-estimate invariances (GATED; Stan-bound) ---------
# Substantive execution is the nocturnal B9.8 sub-unit under Decision C.iv.
# Enable with GDPAR_RUN_BLOCK9_FIT_FUZZ=1 and a configured cmdstan.

.skip_if_no_block9_fit_fuzz <- function() {
  if (!nzchar(Sys.getenv("GDPAR_RUN_BLOCK9_FIT_FUZZ"))) {
    skip("GDPAR_RUN_BLOCK9_FIT_FUZZ not set; skipping Block 9.1 fit fuzz")
  }
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("posterior")
}

.block9_eb_location <- function(fit) {
  # Robust scalar summary: theta_ref estimate plus the additive-block
  # posterior-mean coefficients (NULL-safe).
  est <- as.numeric(fit$theta_ref_hat)
  co <- coef(fit)
  if (!is.null(co$a) && !is.null(co$a$estimate)) {
    est <- c(est, as.numeric(co$a$estimate))
  }
  est
}

test_that("INV-C1 (gated): EB location estimate is permutation-invariant up to MC error", {
  .skip_if_no_block9_fit_fuzz()
  df <- make_gaussian_data(n = 200, seed = 4242)
  fam <- gdpar_family("gaussian")
  fit0 <- suppressMessages(gdpar_eb(y ~ x1, data = df, family = fam, seed = 1L))
  set.seed(7)
  perm <- sample(nrow(df))
  fit1 <- suppressMessages(
    gdpar_eb(y ~ x1, data = df[perm, , drop = FALSE], family = fam, seed = 1L)
  )
  # The likelihood is exchangeable in the rows; posterior means agree up
  # to Monte-Carlo error.
  expect_equal(.block9_eb_location(fit1), .block9_eb_location(fit0),
               tolerance = 0.1)
})

test_that("INV-C2 (gated): EB location estimate is continuous under small noise injection", {
  .skip_if_no_block9_fit_fuzz()
  df <- make_gaussian_data(n = 200, seed = 909)
  fam <- gdpar_family("gaussian")
  fit0 <- suppressMessages(gdpar_eb(y ~ x1, data = df, family = fam, seed = 1L))
  df_noisy <- df
  set.seed(31)
  df_noisy$y <- df_noisy$y + rnorm(nrow(df), sd = 1e-3)
  fit1 <- suppressMessages(
    gdpar_eb(y ~ x1, data = df_noisy, family = fam, seed = 1L)
  )
  # An O(1e-3) perturbation of the outcome must produce an O(1e-2) change
  # in the location estimate (Lipschitz continuity of the estimator).
  expect_equal(.block9_eb_location(fit1), .block9_eb_location(fit0),
               tolerance = 0.1)
})

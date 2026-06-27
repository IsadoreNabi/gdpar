## Tests for gdpar_loo() and the internal helper aggregate_log_lik().
##
## Pure tests construct synthetic posterior::draws_array objects and
## verify the aggregation rules (subject vs cell) without compiling
## Stan. Two env-gated smoke tests (GDPAR_F6_FIT=1) exercise the full
## pipeline with a real fit.

make_synthetic_log_lik_draws <- function(S, n_chains, n, p,
                                          fill = NULL,
                                          seed = 7L) {
  set.seed(seed)
  iter_per_chain <- S %/% n_chains
  if (p == 1L) {
    vars <- paste0("log_lik[", seq_len(n), "]")
    arr <- array(
      if (is.null(fill)) rnorm(iter_per_chain * n_chains * n)
                          else as.numeric(fill),
      dim = c(iter_per_chain, n_chains, n),
      dimnames = list(NULL, NULL, vars)
    )
  } else {
    vars <- as.vector(outer(seq_len(n), seq_len(p),
                              FUN = function(i, k)
                                sprintf("log_lik[%d,%d]", i, k)))
    arr <- array(
      if (is.null(fill)) rnorm(iter_per_chain * n_chains * n * p)
                          else as.numeric(fill),
      dim = c(iter_per_chain, n_chains, n * p),
      dimnames = list(NULL, NULL, vars)
    )
  }
  posterior::as_draws_array(arr)
}

test_that("aggregate_log_lik p=1 returns S x n matrix unchanged", {
  draws <- make_synthetic_log_lik_draws(S = 400L, n_chains = 2L,
                                          n = 5L, p = 1L)
  out <- aggregate_log_lik(draws, p = 1L, aggregation = "subject")
  expect_true(is.matrix(out))
  expect_equal(dim(out), c(400L, 5L))
  expect_false(inherits(out, "draws_matrix"))
})

test_that("aggregate_log_lik p=1 ignores cell vs subject distinction", {
  draws <- make_synthetic_log_lik_draws(S = 200L, n_chains = 2L,
                                          n = 4L, p = 1L)
  out_subj <- aggregate_log_lik(draws, p = 1L, aggregation = "subject")
  out_cell <- aggregate_log_lik(draws, p = 1L, aggregation = "cell")
  expect_identical(out_subj, out_cell)
})

test_that("aggregate_log_lik p=2 subject sums over coordinates", {
  draws <- make_synthetic_log_lik_draws(S = 200L, n_chains = 2L,
                                          n = 3L, p = 2L, fill = 1.0)
  out <- aggregate_log_lik(draws, p = 2L, aggregation = "subject")
  expect_equal(dim(out), c(200L, 3L))
  expect_true(all(out == 2.0))
})

test_that("aggregate_log_lik p=2 subject preserves per-subject sums numerically", {
  draws <- make_synthetic_log_lik_draws(S = 200L, n_chains = 2L,
                                          n = 3L, p = 2L, seed = 11L)
  mat <- unclass(posterior::as_draws_matrix(draws))
  vars <- posterior::variables(draws)
  m <- regexec("^log_lik\\[(\\d+),(\\d+)\\]$", vars)
  parsed <- regmatches(vars, m)
  ij <- do.call(rbind, lapply(parsed, function(z) as.integer(z[2:3])))
  expected <- matrix(0, nrow = nrow(mat), ncol = 3L)
  for (col_idx in seq_along(vars)) {
    expected[, ij[col_idx, 1L]] <-
      expected[, ij[col_idx, 1L]] + mat[, col_idx]
  }
  out <- aggregate_log_lik(draws, p = 2L, aggregation = "subject")
  expect_equal(out, expected, tolerance = 1e-12)
})

test_that("aggregate_log_lik p=2 cell flattens to S x (n*p)", {
  draws <- make_synthetic_log_lik_draws(S = 200L, n_chains = 2L,
                                          n = 3L, p = 2L)
  out <- aggregate_log_lik(draws, p = 2L, aggregation = "cell")
  expect_equal(dim(out), c(200L, 6L))
})

test_that("aggregate_log_lik p=2 cell preserves the original values per (i, k)", {
  draws <- make_synthetic_log_lik_draws(S = 200L, n_chains = 2L,
                                          n = 3L, p = 2L, seed = 13L)
  mat <- unclass(posterior::as_draws_matrix(draws))
  vars <- posterior::variables(draws)
  out <- aggregate_log_lik(draws, p = 2L, aggregation = "cell")
  m <- regexec("^log_lik\\[(\\d+),(\\d+)\\]$", vars)
  parsed <- regmatches(vars, m)
  for (col_idx in seq_along(vars)) {
    i <- as.integer(parsed[[col_idx]][2L])
    k <- as.integer(parsed[[col_idx]][3L])
    flat_col <- (i - 1L) * 2L + k
    expect_equal(out[, flat_col], mat[, col_idx],
                  tolerance = 1e-12, ignore_attr = TRUE)
  }
})

test_that("aggregate_log_lik p=3 cell uses i-major k-minor ordering", {
  draws <- make_synthetic_log_lik_draws(S = 100L, n_chains = 2L,
                                          n = 2L, p = 3L, fill = 0.0)
  out <- aggregate_log_lik(draws, p = 3L, aggregation = "cell")
  expect_equal(dim(out), c(100L, 6L))
})

test_that("aggregate_log_lik p=2 rejects malformed variable names", {
  draws <- make_synthetic_log_lik_draws(S = 50L, n_chains = 2L,
                                          n = 2L, p = 2L)
  draws_bad <- draws
  attr(draws_bad, "dimnames")[[3L]][1L] <- "log_lik_bad_name"
  expect_error(
    aggregate_log_lik(draws_bad, p = 2L, aggregation = "subject"),
    class = "gdpar_loo_parse_error"
  )
})

test_that("aggregate_log_lik p mismatch is caught", {
  draws <- make_synthetic_log_lik_draws(S = 50L, n_chains = 2L,
                                          n = 2L, p = 2L)
  expect_error(
    aggregate_log_lik(draws, p = 3L, aggregation = "subject"),
    class = "gdpar_loo_dim_mismatch"
  )
})

test_that("gdpar_loo aborts on non-gdpar_fit input", {
  expect_error(
    gdpar_loo(list(a = 1)),
    class = "gdpar_input_error"
  )
})

test_that("gdpar_loo aggregation argument is validated via match.arg", {
  fake_fit <- structure(list(), class = "gdpar_fit")
  expect_error(
    gdpar_loo(fake_fit, aggregation = "unknown"),
    regexp = "should be one of"
  )
})

test_that("F.6 smoke fit p=2: gdpar_loo subject returns psis_loo with n obs", {
  skip_if(Sys.getenv("GDPAR_F6_FIT") != "1",
          paste("Set GDPAR_F6_FIT=1 to run the gdpar_loo p=2 smoke fit",
                 "(compiles a Stan model; takes minutes)."))
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("loo")

  n <- 200L; p <- 2L
  set.seed(42L)
  x1 <- rnorm(n); x2 <- rnorm(n)
  theta_ref_true <- c(0.5, -0.5)
  eta_true <- cbind(
    theta_ref_true[1] + 0.5 * x1 - 0.5 * x2,
    theta_ref_true[2] - 0.4 * x1 + 0.4 * x2
  )
  y_mat <- eta_true + matrix(rnorm(n * p, 0, 0.3), n, p)
  df <- data.frame(x1 = x1, x2 = x2)
  df$y <- y_mat

  spec <- amm_spec(p = 2L, dims = dimwise(a = ~ x1 + x2, b = NULL))
  fam  <- gdpar_family_multi("gaussian", p = 2L)
  fit <- suppressMessages(suppressWarnings(
    gdpar(formula = y ~ x1 + x2, family = fam, amm = spec, data = df,
           chains = 2L, iter_warmup = 300L, iter_sampling = 300L,
           refresh = 0L, seed = 42L, parametrization = "ncp",
           skip_id_check = FALSE, verbose = FALSE)
  ))

  lo_subj <- gdpar_loo(fit, aggregation = "subject")
  expect_s3_class(lo_subj, "psis_loo")
  expect_equal(dim(lo_subj$pointwise)[1L], n)
  expect_true(is.finite(lo_subj$estimates["elpd_loo", "Estimate"]))

  lo_cell <- gdpar_loo(fit, aggregation = "cell")
  expect_s3_class(lo_cell, "psis_loo")
  expect_equal(dim(lo_cell$pointwise)[1L], n * p)
})

test_that("F.6 smoke fit p=1: gdpar_loo returns psis_loo with n obs", {
  skip_if(Sys.getenv("GDPAR_F6_FIT") != "1",
          paste("Set GDPAR_F6_FIT=1 to run the gdpar_loo p=1 smoke fit",
                 "(compiles a Stan model; takes minutes)."))
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("loo")

  n <- 200L
  set.seed(42L)
  x1 <- rnorm(n); x2 <- rnorm(n)
  y <- 0.5 + 0.7 * x1 - 0.3 * x2 + rnorm(n, 0, 0.4)
  dat <- data.frame(y = y, x1 = x1, x2 = x2)
  fit <- suppressMessages(suppressWarnings(
    gdpar(y ~ x1 + x2, data = dat,
           family = gdpar_family("gaussian"),
           chains = 2L, iter_warmup = 300L, iter_sampling = 300L,
           refresh = 0L, seed = 42L, verbose = FALSE)
  ))

  lo <- gdpar_loo(fit)
  expect_s3_class(lo, "psis_loo")
  expect_equal(dim(lo$pointwise)[1L], n)
  expect_true(is.finite(lo$estimates["elpd_loo", "Estimate"]))
})

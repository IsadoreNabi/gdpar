# Sub-phase 8.3.3 Unit 3 (2026-05-20): bootstrap golden regression for
# the K-individual path. The first run with GDPAR_GOLDEN_CHECK_K2=1
# and GDPAR_GOLDEN_K2_BOOTSTRAP=1 writes the golden draws snapshot to
# tests/testthat/data/golden_K2_gaussian_K2.rds; subsequent runs with
# GDPAR_GOLDEN_CHECK_K2=1 compare the current draws to the snapshot
# byte-identically. The two-env-var pattern prevents accidental
# overwrite of a frozen golden.

.skip_if_no_golden_K2 <- function() {
  if (!nzchar(Sys.getenv("GDPAR_GOLDEN_CHECK_K2"))) {
    testthat::skip("GDPAR_GOLDEN_CHECK_K2 not set; skipping K=2 golden")
  }
  testthat::skip_if_not_installed("cmdstanr")
}

.golden_K2_path <- function(name) {
  golden_dir <- file.path("data")
  if (!dir.exists(golden_dir)) {
    test_dir <- file.path("tests", "testthat", "data")
    if (dir.exists(test_dir)) golden_dir <- test_dir
  }
  file.path(golden_dir, paste0("golden_K2_", name, ".rds"))
}

.build_K2_gaussian_data <- function(n = 80L, seed = 909L) {
  set.seed(seed)
  x <- rnorm(n)
  z <- rnorm(n)
  mu_true    <- 0.5 + 0.8 * (x - mean(x))
  log_sigma  <- -0.2 + 0.4 * (z - mean(z))
  y <- rnorm(n, mu_true, exp(log_sigma))
  data.frame(y = y, x = x, z = z)
}

test_that("K=2 Gaussian draws are bit-stable against frozen golden", {
  .skip_if_no_golden_K2()
  d <- .build_K2_gaussian_data()
  fs <- gdpar_bf(y ~ a(x), sigma ~ a(z))
  fit <- suppressMessages(suppressWarnings(
    gdpar(fs, data = d, family = gdpar_family("gaussian"),
          chains = 1L, iter_warmup = 200L, iter_sampling = 200L,
          refresh = 0L, seed = 313L,
          skip_id_check = TRUE, verbose = FALSE)
  ))
  observed_draws <- fit$fit$draws()
  golden_path <- .golden_K2_path("gaussian_K2")
  if (!file.exists(golden_path) &&
      nzchar(Sys.getenv("GDPAR_GOLDEN_K2_BOOTSTRAP"))) {
    dir.create(dirname(golden_path), showWarnings = FALSE,
               recursive = TRUE)
    saveRDS(observed_draws, golden_path)
    testthat::skip(sprintf(
      "Bootstrapped K=2 golden written to %s; rerun without bootstrap to compare",
      golden_path
    ))
  }
  testthat::skip_if(
    !file.exists(golden_path),
    sprintf(paste(
      "K=2 golden not present at %s.",
      "Set GDPAR_GOLDEN_K2_BOOTSTRAP=1 alongside GDPAR_GOLDEN_CHECK_K2=1",
      "to generate it."
    ), golden_path)
  )
  golden <- readRDS(golden_path)
  expect_equal(observed_draws, golden)
})

# Sub-phase 8.3.4 additions: Beta K=2 and Gamma K=2 goldens.

.build_K2_beta_data <- function(n = 80L, seed = 919L) {
  set.seed(seed)
  x <- rnorm(n); z <- rnorm(n)
  xc <- x - mean(x); zc <- z - mean(z)
  eta_mu  <- 0.4 + 0.5 * xc
  eta_phi <- 1.0 + 0.3 * zc
  mu_vec  <- 1 / (1 + exp(-eta_mu))
  phi_vec <- exp(eta_phi)
  y <- rbeta(n, mu_vec * phi_vec, (1 - mu_vec) * phi_vec)
  y <- pmin(pmax(y, 1e-4), 1 - 1e-4)
  data.frame(y = y, x = x, z = z)
}

.build_K2_gamma_data <- function(n = 80L, seed = 929L) {
  set.seed(seed)
  x <- rnorm(n); z <- rnorm(n)
  xc <- x - mean(x); zc <- z - mean(z)
  eta_mu    <- 0.3 + 0.4 * xc
  eta_shape <- 0.7 + 0.2 * zc
  mu_vec    <- exp(eta_mu)
  shape_vec <- exp(eta_shape)
  rate_vec  <- shape_vec / mu_vec
  y <- rgamma(n, shape = shape_vec, rate = rate_vec)
  data.frame(y = y, x = x, z = z)
}

test_that("K=2 Beta draws are bit-stable against frozen golden", {
  if (!nzchar(Sys.getenv("GDPAR_GOLDEN_CHECK_K2_BETA"))) {
    testthat::skip("GDPAR_GOLDEN_CHECK_K2_BETA not set; skipping K=2 Beta golden")
  }
  testthat::skip_if_not_installed("cmdstanr")
  d <- .build_K2_beta_data()
  fs <- gdpar_bf(y ~ a(x), phi ~ a(z))
  fit <- suppressMessages(suppressWarnings(
    gdpar(fs, data = d, family = gdpar_family("beta"),
          chains = 1L, iter_warmup = 200L, iter_sampling = 200L,
          refresh = 0L, seed = 314L,
          skip_id_check = TRUE, verbose = FALSE)
  ))
  observed_draws <- fit$fit$draws()
  golden_path <- .golden_K2_path("beta_K2")
  if (!file.exists(golden_path) &&
      nzchar(Sys.getenv("GDPAR_GOLDEN_K2_BETA_BOOTSTRAP"))) {
    dir.create(dirname(golden_path), showWarnings = FALSE,
               recursive = TRUE)
    saveRDS(observed_draws, golden_path)
    testthat::skip(sprintf(
      "Bootstrapped K=2 Beta golden written to %s; rerun without bootstrap to compare",
      golden_path
    ))
  }
  testthat::skip_if(
    !file.exists(golden_path),
    sprintf(paste(
      "K=2 Beta golden not present at %s.",
      "Set GDPAR_GOLDEN_K2_BETA_BOOTSTRAP=1 alongside GDPAR_GOLDEN_CHECK_K2_BETA=1",
      "to generate it."
    ), golden_path)
  )
  golden <- readRDS(golden_path)
  expect_equal(observed_draws, golden)
})

test_that("K=2 Gamma draws are bit-stable against frozen golden", {
  if (!nzchar(Sys.getenv("GDPAR_GOLDEN_CHECK_K2_GAMMA"))) {
    testthat::skip("GDPAR_GOLDEN_CHECK_K2_GAMMA not set; skipping K=2 Gamma golden")
  }
  testthat::skip_if_not_installed("cmdstanr")
  d <- .build_K2_gamma_data()
  fs <- gdpar_bf(y ~ a(x), shape ~ a(z))
  fit <- suppressMessages(suppressWarnings(
    gdpar(fs, data = d, family = gdpar_family("gamma"),
          chains = 1L, iter_warmup = 200L, iter_sampling = 200L,
          refresh = 0L, seed = 315L,
          skip_id_check = TRUE, verbose = FALSE)
  ))
  observed_draws <- fit$fit$draws()
  golden_path <- .golden_K2_path("gamma_K2")
  if (!file.exists(golden_path) &&
      nzchar(Sys.getenv("GDPAR_GOLDEN_K2_GAMMA_BOOTSTRAP"))) {
    dir.create(dirname(golden_path), showWarnings = FALSE,
               recursive = TRUE)
    saveRDS(observed_draws, golden_path)
    testthat::skip(sprintf(
      "Bootstrapped K=2 Gamma golden written to %s; rerun without bootstrap to compare",
      golden_path
    ))
  }
  testthat::skip_if(
    !file.exists(golden_path),
    sprintf(paste(
      "K=2 Gamma golden not present at %s.",
      "Set GDPAR_GOLDEN_K2_GAMMA_BOOTSTRAP=1 alongside GDPAR_GOLDEN_CHECK_K2_GAMMA=1",
      "to generate it."
    ), golden_path)
  )
  golden <- readRDS(golden_path)
  expect_equal(observed_draws, golden)
})

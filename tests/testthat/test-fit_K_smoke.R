# Sub-phase 8.3.3 Unit 3 (2026-05-20): Stan compilation + smoke fits
# for the K-individual path. Mirror of test-fit_multi_smoke.R adapted
# to the K > 1 + p = 1 regime that routes through
# inst/stan/amm_distrib_K.stan.

test_that("K=2 Gaussian brms-style compiles, samples, and recovers theta_ref", {
  testthat::skip_if(Sys.getenv("GDPAR_RUN_STAN_SMOKE") != "1",
                     "Set GDPAR_RUN_STAN_SMOKE=1 to run the K=2 smoke fit (compiles a Stan model; takes minutes).")
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("cmdstanr")

  set.seed(42L)
  n <- 200L
  x <- rnorm(n)
  z <- rnorm(n)
  x_c <- x - mean(x)
  z_c <- z - mean(z)
  theta_ref_mu_true    <- 0.5
  theta_ref_sigma_true <- -0.3
  beta_a_mu_true       <- 0.8
  beta_a_sigma_true    <- 0.4
  eta_mu    <- theta_ref_mu_true    + beta_a_mu_true    * x_c
  eta_sigma <- theta_ref_sigma_true + beta_a_sigma_true * z_c
  y <- rnorm(n, eta_mu, exp(eta_sigma))
  d <- data.frame(y = y, x = x, z = z)

  fs <- gdpar_bf(y ~ a(x), sigma ~ a(z))
  fit <- suppressMessages(suppressWarnings(
    gdpar(fs, data = d, family = gdpar_family("gaussian"),
          chains = 2L, iter_warmup = 500L, iter_sampling = 500L,
          refresh = 0L, seed = 31415L,
          skip_id_check = TRUE, verbose = FALSE)
  ))
  expect_s3_class(fit, "gdpar_fit")
  expect_equal(fit$K, 2L)
  expect_equal(fit$slot_names, c("mu", "sigma"))
  expect_named(fit$stan_data,
               c("n", "K", "family_id_k", "inv_link_id_per_slot",
                 "use_a_k", "use_b_k", "use_W",
                 "J_a_per_k", "J_b_per_k", "J_a_max", "J_b_max",
                 "Z_a_k", "Z_b_k", "dim_W", "d", "X",
                 "y_real", "y_int", "theta_anchor_K",
                 "use_dispersion_y_k", "use_dispersion_phi_k",
                 "W_type_id", "W_n_knots_full", "W_knots_full",
                 "W_degree", "use_groups", "J_groups", "group_id"))
  # Sanity recovery of the per-slot anchor: posterior median within 0.3
  # of the data-generating value. This is a smoke threshold, not a
  # calibration claim; tighter coverage is checked in the formal
  # validation of Bloque 9.
  draws <- fit$fit$draws(format = "df")
  theta_mu_med    <- median(draws[["theta_ref_k[1,1]"]])
  theta_sigma_med <- median(draws[["theta_ref_k[1,2]"]])
  expect_lt(abs(theta_mu_med - theta_ref_mu_true), 0.3)
  expect_lt(abs(theta_sigma_med - theta_ref_sigma_true), 0.3)
})

# Sub-phase 8.3.4 (2026-05-20 / 21): smoke fits for Beta K=2 and
# Gamma K=2 bi-parametric families plus the descriptor-based custom
# lognormal_loc_scale (D-A3.B). Each gated by its own env var so the
# user can opt into one family at a time.

test_that("K=2 Beta brms-style compiles, samples, and recovers anchors", {
  testthat::skip_if(Sys.getenv("GDPAR_RUN_STAN_SMOKE_BETA") != "1",
                     "Set GDPAR_RUN_STAN_SMOKE_BETA=1 to run the K=2 Beta smoke fit.")
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("cmdstanr")

  set.seed(43L)
  n <- 250L
  x <- rnorm(n); z <- rnorm(n)
  xc <- x - mean(x); zc <- z - mean(z)
  theta_mu_true  <- 0.4
  theta_phi_true <- 1.0
  beta_a_mu_true  <- 0.5
  beta_a_phi_true <- 0.3
  eta_mu  <- theta_mu_true  + beta_a_mu_true  * xc
  eta_phi <- theta_phi_true + beta_a_phi_true * zc
  mu_vec  <- 1 / (1 + exp(-eta_mu))
  phi_vec <- exp(eta_phi)
  y <- rbeta(n, mu_vec * phi_vec, (1 - mu_vec) * phi_vec)
  y <- pmin(pmax(y, 1e-4), 1 - 1e-4)
  d <- data.frame(y = y, x = x, z = z)
  fs <- gdpar_bf(y ~ a(x), phi ~ a(z))
  fit <- suppressMessages(suppressWarnings(
    gdpar(fs, data = d, family = gdpar_family("beta"),
          chains = 2L, iter_warmup = 500L, iter_sampling = 500L,
          refresh = 0L, seed = 31416L,
          skip_id_check = TRUE, verbose = FALSE)
  ))
  expect_s3_class(fit, "gdpar_fit")
  expect_equal(fit$K, 2L)
  expect_equal(fit$slot_names, c("mu", "phi"))
  expect_equal(fit$stan_data$family_id_k, c(5L, 5L))
  draws <- fit$fit$draws(format = "df")
  theta_mu_med  <- median(draws[["theta_ref_k[1,1]"]])
  theta_phi_med <- median(draws[["theta_ref_k[1,2]"]])
  expect_lt(abs(theta_mu_med  - theta_mu_true),  0.4)
  expect_lt(abs(theta_phi_med - theta_phi_true), 0.4)
})

test_that("K=2 Gamma brms-style compiles, samples, and recovers anchors", {
  testthat::skip_if(Sys.getenv("GDPAR_RUN_STAN_SMOKE_GAMMA") != "1",
                     "Set GDPAR_RUN_STAN_SMOKE_GAMMA=1 to run the K=2 Gamma smoke fit.")
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("cmdstanr")

  set.seed(44L)
  n <- 250L
  x <- rnorm(n); z <- rnorm(n)
  xc <- x - mean(x); zc <- z - mean(z)
  theta_mu_true    <- 0.3   # log mu around exp(0.3) ~ 1.35
  theta_shape_true <- 0.7   # log shape around 2.0
  beta_a_mu_true    <- 0.4
  beta_a_shape_true <- 0.2
  eta_mu    <- theta_mu_true    + beta_a_mu_true    * xc
  eta_shape <- theta_shape_true + beta_a_shape_true * zc
  mu_vec    <- exp(eta_mu)
  shape_vec <- exp(eta_shape)
  rate_vec  <- shape_vec / mu_vec
  y <- rgamma(n, shape = shape_vec, rate = rate_vec)
  d <- data.frame(y = y, x = x, z = z)
  fs <- gdpar_bf(y ~ a(x), shape ~ a(z))
  fit <- suppressMessages(suppressWarnings(
    gdpar(fs, data = d, family = gdpar_family("gamma"),
          chains = 2L, iter_warmup = 500L, iter_sampling = 500L,
          refresh = 0L, seed = 31417L,
          skip_id_check = TRUE, verbose = FALSE)
  ))
  expect_s3_class(fit, "gdpar_fit")
  expect_equal(fit$K, 2L)
  expect_equal(fit$slot_names, c("mu", "shape"))
  expect_equal(fit$stan_data$family_id_k, c(6L, 6L))
  draws <- fit$fit$draws(format = "df")
  theta_mu_med    <- median(draws[["theta_ref_k[1,1]"]])
  theta_shape_med <- median(draws[["theta_ref_k[1,2]"]])
  expect_lt(abs(theta_mu_med    - theta_mu_true),    0.4)
  expect_lt(abs(theta_shape_med - theta_shape_true), 0.4)
})

test_that("K=2 custom lognormal_loc_scale compiles, samples, and recovers anchors", {
  testthat::skip_if(Sys.getenv("GDPAR_RUN_STAN_SMOKE_CUSTOM") != "1",
                     "Set GDPAR_RUN_STAN_SMOKE_CUSTOM=1 to run the K=2 custom lognormal smoke fit.")
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("cmdstanr")

  set.seed(45L)
  n <- 250L
  x <- rnorm(n); z <- rnorm(n)
  xc <- x - mean(x); zc <- z - mean(z)
  theta_mu_true    <- 0.2
  theta_sigma_true <- -0.3
  beta_a_mu_true    <- 0.4
  beta_a_sigma_true <- 0.2
  eta_mu    <- theta_mu_true    + beta_a_mu_true    * xc
  eta_sigma <- theta_sigma_true + beta_a_sigma_true * zc
  y <- rlnorm(n, eta_mu, exp(eta_sigma))
  d <- data.frame(y = y, x = x, z = z)
  fam_ln <- suppressMessages(gdpar_family_custom_K(
    name = "smoke_ln", stan_lpdf_id = "lognormal_loc_scale",
    did_holds = TRUE, did_reference = "test"
  ))
  fs <- gdpar_bf(y ~ a(x), sigma ~ a(z))
  fit <- suppressMessages(suppressWarnings(
    gdpar(fs, data = d, family = fam_ln,
          chains = 2L, iter_warmup = 500L, iter_sampling = 500L,
          refresh = 0L, seed = 31418L,
          skip_id_check = TRUE, verbose = FALSE)
  ))
  expect_s3_class(fit, "gdpar_fit")
  expect_equal(fit$K, 2L)
  expect_equal(fit$slot_names, c("mu", "sigma"))
  expect_equal(fit$stan_data$family_id_k, c(7L, 7L))
  draws <- fit$fit$draws(format = "df")
  theta_mu_med    <- median(draws[["theta_ref_k[1,1]"]])
  theta_sigma_med <- median(draws[["theta_ref_k[1,2]"]])
  expect_lt(abs(theta_mu_med    - theta_mu_true),    0.4)
  expect_lt(abs(theta_sigma_med - theta_sigma_true), 0.4)
})

# Sub-phase 8.3.5a (2026-05-21): smoke fit for Student-t K=3 wired in
# amm_distrib_K.stan branch family_id_k == 8 with the canonical
# (mu identity, sigma log, nu log) parametrization. Gated by a
# dedicated env var so the user can opt into the Student-t path
# without compiling Stan models for the K=2 paths.

test_that("K=3 Student-t brms-style compiles, samples, and recovers anchors", {
  testthat::skip_if(Sys.getenv("GDPAR_RUN_STAN_SMOKE_STUDENT_T") != "1",
                     "Set GDPAR_RUN_STAN_SMOKE_STUDENT_T=1 to run the K=3 Student-t smoke fit.")
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("cmdstanr")

  set.seed(46L)
  n <- 300L
  x <- rnorm(n); z <- rnorm(n); w <- rnorm(n)
  xc <- x - mean(x); zc <- z - mean(z); wc <- w - mean(w)
  theta_mu_true    <- 0.2
  theta_sigma_true <- 0.0
  theta_nu_true    <- log(15)
  beta_a_mu_true    <- 0.4
  beta_a_sigma_true <- 0.2
  beta_a_nu_true    <- 0.1
  eta_mu    <- theta_mu_true    + beta_a_mu_true    * xc
  eta_sigma <- theta_sigma_true + beta_a_sigma_true * zc
  eta_nu    <- theta_nu_true    + beta_a_nu_true    * wc
  sigma_vec <- exp(eta_sigma)
  nu_vec    <- exp(eta_nu)
  y <- eta_mu + sigma_vec * rt(n, df = nu_vec)
  d <- data.frame(y = y, x = x, z = z, w = w)
  fs <- gdpar_bf(y ~ a(x), sigma ~ a(z), nu ~ a(w))
  fit <- suppressMessages(suppressWarnings(
    gdpar(fs, data = d, family = gdpar_family("student_t"),
          chains = 2L, iter_warmup = 500L, iter_sampling = 500L,
          refresh = 0L, seed = 31419L,
          skip_id_check = TRUE, verbose = FALSE)
  ))
  expect_s3_class(fit, "gdpar_fit")
  expect_equal(fit$K, 3L)
  expect_equal(fit$slot_names, c("mu", "sigma", "nu"))
  expect_equal(fit$stan_data$family_id_k, c(8L, 8L, 8L))
  draws <- fit$fit$draws(format = "df")
  theta_mu_med    <- median(draws[["theta_ref_k[1,1]"]])
  theta_sigma_med <- median(draws[["theta_ref_k[1,2]"]])
  theta_nu_med    <- median(draws[["theta_ref_k[1,3]"]])
  expect_lt(abs(theta_mu_med    - theta_mu_true),    0.4)
  expect_lt(abs(theta_sigma_med - theta_sigma_true), 0.4)
  expect_lt(abs(theta_nu_med    - theta_nu_true),    1.0)
})

# Sub-phase 8.3.5b (2026-05-21): smoke fit for Tweedie K=3 wired in
# amm_distrib_K.stan branch family_id_k == 9 with the canonical
# (mu log, phi log, p identity) parametrization. The hybrid lpdf
# uses the Dunn--Smyth series in |p - 1.5| < 0.4 and saddlepoint
# elsewhere. Gated by GDPAR_RUN_STAN_SMOKE_TWEEDIE; the slot p
# anchor is sampled with uniform(1.01, 1.99) prior per decision E7.

test_that("K=3 Tweedie brms-style compiles, samples, and recovers anchors", {
  testthat::skip_if(Sys.getenv("GDPAR_RUN_STAN_SMOKE_TWEEDIE") != "1",
                     "Set GDPAR_RUN_STAN_SMOKE_TWEEDIE=1 to run the K=3 Tweedie smoke fit.")
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("cmdstanr")

  # Compound Poisson-gamma generator: Y = sum_{i=1}^N G_i with
  # N ~ Poisson(lambda) and G_i ~ Gamma(shape, rate); the resulting
  # marginal is exact Tweedie for 1 < p < 2.
  rtweedie <- function(n, mu, phi, p) {
    lambda <- mu^(2 - p) / (phi * (2 - p))
    shape  <- (2 - p) / (p - 1)
    rate   <- 1 / (phi * (p - 1) * mu^(p - 1))
    N <- rpois(n, lambda)
    out <- numeric(n)
    pos <- N > 0L
    if (any(pos)) {
      out[pos] <- vapply(which(pos), function(i) {
        rgamma(1L, shape = N[i] * shape, rate = rate[i])
      }, numeric(1L))
    }
    out
  }

  set.seed(47L)
  n <- 300L
  x <- rnorm(n); z <- rnorm(n); w <- rnorm(n)
  xc <- x - mean(x); zc <- z - mean(z); wc <- w - mean(w)
  theta_mu_true    <- 1.0       # log scale; mu = exp(eta_mu)
  theta_phi_true   <- 0.0       # log scale; phi = exp(eta_phi)
  theta_p_true     <- 1.5       # identity; constrained to (1.01, 1.99)
  beta_a_mu_true   <- 0.3
  beta_a_phi_true  <- 0.2
  beta_a_p_true    <- 0.05      # small slope to keep p well inside the range
  eta_mu   <- theta_mu_true   + beta_a_mu_true  * xc
  eta_phi  <- theta_phi_true  + beta_a_phi_true * zc
  eta_p    <- theta_p_true    + beta_a_p_true   * wc
  mu_vec   <- exp(eta_mu)
  phi_vec  <- exp(eta_phi)
  p_vec    <- pmin(pmax(eta_p, 1.05), 1.95)
  y <- rtweedie(n, mu = mu_vec, phi = phi_vec, p = p_vec)
  d <- data.frame(y = y, x = x, z = z, w = w)
  fs <- gdpar_bf(y ~ a(x), phi ~ a(z), p ~ a(w))
  fit <- suppressMessages(suppressWarnings(
    gdpar(fs, data = d, family = gdpar_family("tweedie"),
          chains = 2L, iter_warmup = 500L, iter_sampling = 500L,
          refresh = 0L, seed = 31420L,
          skip_id_check = TRUE, verbose = FALSE)
  ))
  expect_s3_class(fit, "gdpar_fit")
  expect_equal(fit$K, 3L)
  expect_equal(fit$slot_names, c("mu", "phi", "p"))
  expect_equal(fit$stan_data$family_id_k, c(9L, 9L, 9L))
  draws <- fit$fit$draws(format = "df")
  theta_mu_med  <- median(draws[["theta_ref_k[1,1]"]])
  theta_phi_med <- median(draws[["theta_ref_k[1,2]"]])
  theta_p_med   <- median(draws[["theta_ref_k[1,3]"]])
  expect_lt(abs(theta_mu_med  - theta_mu_true),  0.5)
  expect_lt(abs(theta_phi_med - theta_phi_true), 0.6)
  expect_lt(abs(theta_p_med   - theta_p_true),   0.15)
})

# Sub-phase 8.3.6 (2026-05-21): smoke fits for the four count mixture
# families wired in amm_distrib_K.stan: ZIP K=2 (stan_id 10), ZINB K=3
# (stan_id 11), Hurdle-Poisson K=2 (stan_id 12), Hurdle-NB K=3
# (stan_id 13). Each smoke generates synthetic data with the canonical
# mixture mechanism on the R side using only built-in random
# generators (rpois, rnbinom, rbinom) per gotcha 79; rejection
# sampling realizes the truncated-at-one positive branch of the
# hurdle families. Each smoke is gated by its dedicated env var.

test_that("K=2 ZIP brms-style compiles, samples, and recovers anchors", {
  testthat::skip_if(Sys.getenv("GDPAR_RUN_STAN_SMOKE_ZIP") != "1",
                     "Set GDPAR_RUN_STAN_SMOKE_ZIP=1 to run the K=2 ZIP smoke fit.")
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("cmdstanr")

  set.seed(51L)
  n <- 300L
  x <- rnorm(n); z <- rnorm(n)
  xc <- x - mean(x); zc <- z - mean(z)
  theta_mu_true   <- 1.2          # log scale; mu = exp(eta_mu)
  theta_pi_true   <- -1.0         # logit scale; pi = inv_logit(eta_pi)
  beta_a_mu_true  <- 0.4
  beta_a_pi_true  <- 0.3
  eta_mu  <- theta_mu_true + beta_a_mu_true * xc
  eta_pi  <- theta_pi_true + beta_a_pi_true * zc
  mu_vec  <- exp(eta_mu)
  pi_vec  <- 1 / (1 + exp(-eta_pi))
  is_structural_zero <- rbinom(n, size = 1L, prob = pi_vec) == 1L
  y <- ifelse(is_structural_zero, 0L, rpois(n, lambda = mu_vec))
  d <- data.frame(y = y, x = x, z = z)
  fs <- gdpar_bf(y ~ a(x), pi ~ a(z))
  fit <- suppressMessages(suppressWarnings(
    gdpar(fs, data = d, family = gdpar_family("zip"),
          chains = 2L, iter_warmup = 500L, iter_sampling = 500L,
          refresh = 0L, seed = 41510L,
          skip_id_check = TRUE, verbose = FALSE)
  ))
  expect_s3_class(fit, "gdpar_fit")
  expect_equal(fit$K, 2L)
  expect_equal(fit$slot_names, c("mu", "pi"))
  expect_equal(fit$stan_data$family_id_k, c(10L, 10L))
  draws <- fit$fit$draws(format = "df")
  theta_mu_med <- median(draws[["theta_ref_k[1,1]"]])
  theta_pi_med <- median(draws[["theta_ref_k[1,2]"]])
  expect_lt(abs(theta_mu_med - theta_mu_true), 0.5)
  expect_lt(abs(theta_pi_med - theta_pi_true), 0.6)
})

test_that("K=3 ZINB brms-style compiles, samples, and recovers anchors", {
  testthat::skip_if(Sys.getenv("GDPAR_RUN_STAN_SMOKE_ZINB") != "1",
                     "Set GDPAR_RUN_STAN_SMOKE_ZINB=1 to run the K=3 ZINB smoke fit.")
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("cmdstanr")

  set.seed(52L)
  n <- 300L
  x <- rnorm(n); z <- rnorm(n); w <- rnorm(n)
  xc <- x - mean(x); zc <- z - mean(z); wc <- w - mean(w)
  theta_mu_true    <- 1.0         # log; mu = exp(eta_mu)
  theta_phi_true   <- 0.5         # log; phi = exp(eta_phi)
  theta_pi_true    <- -1.2        # logit; pi = inv_logit(eta_pi)
  beta_a_mu_true   <- 0.3
  beta_a_phi_true  <- 0.2
  beta_a_pi_true   <- 0.3
  eta_mu   <- theta_mu_true  + beta_a_mu_true  * xc
  eta_phi  <- theta_phi_true + beta_a_phi_true * zc
  eta_pi   <- theta_pi_true  + beta_a_pi_true  * wc
  mu_vec   <- exp(eta_mu)
  phi_vec  <- exp(eta_phi)
  pi_vec   <- 1 / (1 + exp(-eta_pi))
  is_structural_zero <- rbinom(n, size = 1L, prob = pi_vec) == 1L
  y <- ifelse(is_structural_zero, 0L,
              rnbinom(n, mu = mu_vec, size = phi_vec))
  d <- data.frame(y = y, x = x, z = z, w = w)
  fs <- gdpar_bf(y ~ a(x), phi ~ a(z), pi ~ a(w))
  fit <- suppressMessages(suppressWarnings(
    gdpar(fs, data = d, family = gdpar_family("zinb"),
          chains = 2L, iter_warmup = 500L, iter_sampling = 500L,
          refresh = 0L, seed = 41511L,
          skip_id_check = TRUE, verbose = FALSE)
  ))
  expect_s3_class(fit, "gdpar_fit")
  expect_equal(fit$K, 3L)
  expect_equal(fit$slot_names, c("mu", "phi", "pi"))
  expect_equal(fit$stan_data$family_id_k, c(11L, 11L, 11L))
  draws <- fit$fit$draws(format = "df")
  theta_mu_med  <- median(draws[["theta_ref_k[1,1]"]])
  theta_phi_med <- median(draws[["theta_ref_k[1,2]"]])
  theta_pi_med  <- median(draws[["theta_ref_k[1,3]"]])
  expect_lt(abs(theta_mu_med  - theta_mu_true),  0.5)
  expect_lt(abs(theta_phi_med - theta_phi_true), 0.7)
  expect_lt(abs(theta_pi_med  - theta_pi_true),  0.7)
})

test_that("K=2 Hurdle-Poisson brms-style compiles, samples, and recovers anchors", {
  testthat::skip_if(Sys.getenv("GDPAR_RUN_STAN_SMOKE_HURDLE_P") != "1",
                     "Set GDPAR_RUN_STAN_SMOKE_HURDLE_P=1 to run the K=2 Hurdle-Poisson smoke fit.")
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("cmdstanr")

  # Zero-truncated Poisson via rejection sampling. Robust when mu is
  # not pathologically small (smoke uses mu around exp(1) = 2.7).
  rtpois <- function(n, mu) {
    out <- integer(n)
    for (i in seq_len(n)) {
      val <- 0L
      iter <- 0L
      while (val == 0L && iter < 10000L) {
        val <- rpois(1L, lambda = mu[i])
        iter <- iter + 1L
      }
      out[i] <- val
    }
    out
  }

  set.seed(53L)
  n <- 300L
  x <- rnorm(n); z <- rnorm(n)
  xc <- x - mean(x); zc <- z - mean(z)
  theta_mu_true   <- 1.0          # log
  theta_pi_true   <- -0.5         # logit
  beta_a_mu_true  <- 0.4
  beta_a_pi_true  <- 0.3
  eta_mu <- theta_mu_true + beta_a_mu_true * xc
  eta_pi <- theta_pi_true + beta_a_pi_true * zc
  mu_vec <- exp(eta_mu)
  pi_vec <- 1 / (1 + exp(-eta_pi))
  is_zero <- rbinom(n, size = 1L, prob = pi_vec) == 1L
  y <- integer(n)
  y[is_zero]  <- 0L
  y[!is_zero] <- rtpois(sum(!is_zero), mu = mu_vec[!is_zero])
  d <- data.frame(y = y, x = x, z = z)
  fs <- gdpar_bf(y ~ a(x), pi ~ a(z))
  fit <- suppressMessages(suppressWarnings(
    gdpar(fs, data = d, family = gdpar_family("hurdle_poisson"),
          chains = 2L, iter_warmup = 500L, iter_sampling = 500L,
          refresh = 0L, seed = 41512L,
          skip_id_check = TRUE, verbose = FALSE)
  ))
  expect_s3_class(fit, "gdpar_fit")
  expect_equal(fit$K, 2L)
  expect_equal(fit$stan_data$family_id_k, c(12L, 12L))
  draws <- fit$fit$draws(format = "df")
  theta_mu_med <- median(draws[["theta_ref_k[1,1]"]])
  theta_pi_med <- median(draws[["theta_ref_k[1,2]"]])
  expect_lt(abs(theta_mu_med - theta_mu_true), 0.5)
  expect_lt(abs(theta_pi_med - theta_pi_true), 0.6)
})

test_that("K=3 Hurdle-NB brms-style compiles, samples, and recovers anchors", {
  testthat::skip_if(Sys.getenv("GDPAR_RUN_STAN_SMOKE_HURDLE_NB") != "1",
                     "Set GDPAR_RUN_STAN_SMOKE_HURDLE_NB=1 to run the K=3 Hurdle-NB smoke fit.")
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("cmdstanr")

  rtnbinom <- function(n, mu, size) {
    out <- integer(n)
    for (i in seq_len(n)) {
      val <- 0L
      iter <- 0L
      while (val == 0L && iter < 10000L) {
        val <- rnbinom(1L, mu = mu[i], size = size[i])
        iter <- iter + 1L
      }
      out[i] <- val
    }
    out
  }

  set.seed(54L)
  n <- 300L
  x <- rnorm(n); z <- rnorm(n); w <- rnorm(n)
  xc <- x - mean(x); zc <- z - mean(z); wc <- w - mean(w)
  theta_mu_true    <- 1.0         # log
  theta_phi_true   <- 0.5         # log
  theta_pi_true    <- -0.5        # logit
  beta_a_mu_true   <- 0.3
  beta_a_phi_true  <- 0.2
  beta_a_pi_true   <- 0.3
  eta_mu  <- theta_mu_true  + beta_a_mu_true  * xc
  eta_phi <- theta_phi_true + beta_a_phi_true * zc
  eta_pi  <- theta_pi_true  + beta_a_pi_true  * wc
  mu_vec  <- exp(eta_mu)
  phi_vec <- exp(eta_phi)
  pi_vec  <- 1 / (1 + exp(-eta_pi))
  is_zero <- rbinom(n, size = 1L, prob = pi_vec) == 1L
  y <- integer(n)
  y[is_zero]  <- 0L
  y[!is_zero] <- rtnbinom(sum(!is_zero),
                          mu   = mu_vec[!is_zero],
                          size = phi_vec[!is_zero])
  d <- data.frame(y = y, x = x, z = z, w = w)
  fs <- gdpar_bf(y ~ a(x), phi ~ a(z), pi ~ a(w))
  fit <- suppressMessages(suppressWarnings(
    gdpar(fs, data = d,
          family = gdpar_family("hurdle_neg_binomial_2"),
          chains = 2L, iter_warmup = 500L, iter_sampling = 500L,
          refresh = 0L, seed = 41513L,
          skip_id_check = TRUE, verbose = FALSE)
  ))
  expect_s3_class(fit, "gdpar_fit")
  expect_equal(fit$K, 3L)
  expect_equal(fit$slot_names, c("mu", "phi", "pi"))
  expect_equal(fit$stan_data$family_id_k, c(13L, 13L, 13L))
  draws <- fit$fit$draws(format = "df")
  theta_mu_med  <- median(draws[["theta_ref_k[1,1]"]])
  theta_phi_med <- median(draws[["theta_ref_k[1,2]"]])
  theta_pi_med  <- median(draws[["theta_ref_k[1,3]"]])
  expect_lt(abs(theta_mu_med  - theta_mu_true),  0.5)
  expect_lt(abs(theta_phi_med - theta_phi_true), 0.7)
  expect_lt(abs(theta_pi_med  - theta_pi_true),  0.7)
})

# Sub-phase 8.3.7 (2026-05-21): smoke fits for heterogeneous families
# per slot under the K = 2 location families. Three smokes covering
# the validator's coherence space (D6 = a2):
#   - Gauss + Beta scale: sigma forced into (0, 1) via Beta slot 1
#     (inv_logit). Demonstrates effective heterogeneity (link change
#     log -> logit).
#   - Gauss + Gamma scale: sigma stays in positive_real via Gamma
#     slot 1 (exp). "Virtual" heterogeneity (same inv_link as
#     homogeneous Gaussian K = 2 slot 2) that validates the API
#     dispatcher without exercising a link change.
#   - NB + Beta dispersion: phi forced into (0, 1) via Beta slot 1
#     (inv_logit). Effective heterogeneity for a non-Gaussian location.
# Each smoke gated by its dedicated env var.

test_that("K=2 het Gauss+Beta scale recovers anchors with sigma in (0, 1)", {
  testthat::skip_if(Sys.getenv("GDPAR_RUN_STAN_SMOKE_HET_GAUSS_BETA") != "1",
                     "Set GDPAR_RUN_STAN_SMOKE_HET_GAUSS_BETA=1 to run the K=2 het Gauss+Beta smoke fit.")
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("cmdstanr")

  set.seed(61L)
  n <- 300L
  x <- rnorm(n); z <- rnorm(n)
  xc <- x - mean(x); zc <- z - mean(z)
  theta_mu_true     <- 0.4                       # identity
  theta_sigma_true  <- -0.8                      # logit -> inv_logit(-0.8) ~ 0.31
  beta_a_mu_true    <- 0.4
  beta_a_sigma_true <- 0.3
  eta_mu    <- theta_mu_true    + beta_a_mu_true    * xc
  eta_sigma <- theta_sigma_true + beta_a_sigma_true * zc
  sigma_vec <- 1 / (1 + exp(-eta_sigma))         # heterogeneous Beta link
  y <- rnorm(n, eta_mu, sigma_vec)
  d <- data.frame(y = y, x = x, z = z)
  fs <- gdpar_bf(y ~ a(x), sigma ~ a(z))
  fit <- suppressMessages(suppressWarnings(
    gdpar(fs, data = d,
          family = list(mu = gdpar_family("gaussian"),
                        sigma = gdpar_family("beta")),
          chains = 2L, iter_warmup = 500L, iter_sampling = 500L,
          refresh = 0L, seed = 41521L,
          skip_id_check = TRUE, verbose = FALSE)
  ))
  expect_s3_class(fit, "gdpar_fit")
  expect_equal(fit$K, 2L)
  expect_equal(fit$slot_names, c("mu", "sigma"))
  expect_equal(fit$stan_data$family_id_k, c(1L, 5L))
  expect_equal(fit$stan_data$inv_link_id_per_slot, c(0L, 1L))
  draws <- fit$fit$draws(format = "df")
  theta_mu_med    <- median(draws[["theta_ref_k[1,1]"]])
  theta_sigma_med <- median(draws[["theta_ref_k[1,2]"]])
  expect_lt(abs(theta_mu_med    - theta_mu_true),    0.5)
  expect_lt(abs(theta_sigma_med - theta_sigma_true), 0.7)
})

test_that("K=2 het Gauss+Gamma scale recovers anchors (virtual heterogeneity)", {
  testthat::skip_if(Sys.getenv("GDPAR_RUN_STAN_SMOKE_HET_GAUSS_GAMMA") != "1",
                     "Set GDPAR_RUN_STAN_SMOKE_HET_GAUSS_GAMMA=1 to run the K=2 het Gauss+Gamma smoke fit.")
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("cmdstanr")

  set.seed(62L)
  n <- 300L
  x <- rnorm(n); z <- rnorm(n)
  xc <- x - mean(x); zc <- z - mean(z)
  theta_mu_true     <- 0.5                       # identity
  theta_sigma_true  <- -0.3                      # log -> exp(-0.3) ~ 0.74
  beta_a_mu_true    <- 0.4
  beta_a_sigma_true <- 0.2
  eta_mu    <- theta_mu_true    + beta_a_mu_true    * xc
  eta_sigma <- theta_sigma_true + beta_a_sigma_true * zc
  sigma_vec <- exp(eta_sigma)                    # Gamma slot 1 = exp; identical to
                                                 # homogeneous Gaussian slot 2 (log)
  y <- rnorm(n, eta_mu, sigma_vec)
  d <- data.frame(y = y, x = x, z = z)
  fs <- gdpar_bf(y ~ a(x), sigma ~ a(z))
  fit <- suppressMessages(suppressWarnings(
    gdpar(fs, data = d,
          family = list(mu = gdpar_family("gaussian"),
                        sigma = gdpar_family("gamma")),
          chains = 2L, iter_warmup = 500L, iter_sampling = 500L,
          refresh = 0L, seed = 41522L,
          skip_id_check = TRUE, verbose = FALSE)
  ))
  expect_s3_class(fit, "gdpar_fit")
  expect_equal(fit$K, 2L)
  expect_equal(fit$slot_names, c("mu", "sigma"))
  expect_equal(fit$stan_data$family_id_k, c(1L, 6L))
  expect_equal(fit$stan_data$inv_link_id_per_slot, c(0L, 2L))
  draws <- fit$fit$draws(format = "df")
  theta_mu_med    <- median(draws[["theta_ref_k[1,1]"]])
  theta_sigma_med <- median(draws[["theta_ref_k[1,2]"]])
  expect_lt(abs(theta_mu_med    - theta_mu_true),    0.5)
  expect_lt(abs(theta_sigma_med - theta_sigma_true), 0.7)
})

test_that("K=2 het NB+Beta dispersion recovers anchors with phi in (0, 1)", {
  testthat::skip_if(Sys.getenv("GDPAR_RUN_STAN_SMOKE_HET_NB_BETA") != "1",
                     "Set GDPAR_RUN_STAN_SMOKE_HET_NB_BETA=1 to run the K=2 het NB+Beta smoke fit.")
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("cmdstanr")

  set.seed(63L)
  n <- 300L
  x <- rnorm(n); z <- rnorm(n)
  xc <- x - mean(x); zc <- z - mean(z)
  theta_mu_true    <- 0.8                        # log -> mu around exp(0.8) ~ 2.2
  theta_phi_true   <- 0.3                        # logit -> inv_logit(0.3) ~ 0.57
  beta_a_mu_true   <- 0.3
  beta_a_phi_true  <- 0.2
  eta_mu  <- theta_mu_true  + beta_a_mu_true  * xc
  eta_phi <- theta_phi_true + beta_a_phi_true * zc
  mu_vec  <- exp(eta_mu)
  phi_vec <- 1 / (1 + exp(-eta_phi))             # heterogeneous Beta link
  y <- rnbinom(n, mu = mu_vec, size = phi_vec)
  d <- data.frame(y = y, x = x, z = z)
  fs <- gdpar_bf(y ~ a(x), phi ~ a(z))
  fit <- suppressMessages(suppressWarnings(
    gdpar(fs, data = d,
          family = list(mu = gdpar_family("neg_binomial_2"),
                        phi = gdpar_family("beta")),
          chains = 2L, iter_warmup = 500L, iter_sampling = 500L,
          refresh = 0L, seed = 41523L,
          skip_id_check = TRUE, verbose = FALSE)
  ))
  expect_s3_class(fit, "gdpar_fit")
  expect_equal(fit$K, 2L)
  expect_equal(fit$slot_names, c("mu", "phi"))
  expect_equal(fit$stan_data$family_id_k, c(3L, 5L))
  expect_equal(fit$stan_data$inv_link_id_per_slot, c(2L, 1L))
  draws <- fit$fit$draws(format = "df")
  theta_mu_med  <- median(draws[["theta_ref_k[1,1]"]])
  theta_phi_med <- median(draws[["theta_ref_k[1,2]"]])
  expect_lt(abs(theta_mu_med  - theta_mu_true),  0.5)
  expect_lt(abs(theta_phi_med - theta_phi_true), 0.7)
})

# Sub-phase 8.3.8 (2026-05-22): smoke fits exercising the B-spline W
# basis end-to-end through the three Stan templates. Each gated by its
# own env var. The first one drives K = 1 + p = 1 via amm_main.stan
# directly (low-level path); the second combines heterogeneous family
# (8.3.7) with B-spline W via the named-list family API.

test_that("K=1+p=1 bspline W (Gaussian) compiles, samples, and recovers anchor", {
  testthat::skip_if(Sys.getenv("GDPAR_RUN_STAN_SMOKE_BSPLINE_K1") != "1",
                     "Set GDPAR_RUN_STAN_SMOKE_BSPLINE_K1=1 to run the K=1+p=1 bspline smoke fit.")
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("cmdstanr")

  set.seed(81L)
  n <- 250L
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x1_c <- x1 - mean(x1)
  x2_c <- x2 - mean(x2)
  theta_ref_true <- 0.4
  beta_a_true <- 0.3
  sigma_true <- 0.4
  # Modulating term: W(theta_ref) - W(theta_anchor=0). With bspline the
  # contribution per observation depends on theta_ref (a parameter, not
  # on the covariate); we keep its true contribution flat at zero by
  # zeroing W coefficients in the data-generating model. The smoke
  # therefore validates the wiring (compile + sample + recover anchor)
  # rather than W recovery per se, which belongs in the goldens of
  # 8.3.9.
  eta_true <- theta_ref_true + beta_a_true * x1_c
  y <- rnorm(n, eta_true, sigma_true)
  data_df <- data.frame(x1 = x1, x2 = x2)

  spec <- amm_spec(
    a = ~ x1,
    W = W_basis(type = "bspline", degree = 3, df = 5,
                boundary_knots = c(-10, 10))
  )
  spec$W <- materialize_W_basis(spec$W, 1L)
  fam <- gdpar_family("gaussian")
  prior <- gdpar_prior()

  src <- generate_stan_code(prior, cp_a = FALSE, cp_W = FALSE)
  stan_file <- tempfile(fileext = ".stan")
  writeLines(src, stan_file)
  mod <- cmdstanr::cmdstan_model(stan_file)
  expect_s3_class(mod, "CmdStanModel")

  design <- build_amm_design(spec, data_df, formula_rhs = ~ x2)
  stan_data <- assemble_stan_data(design, fam, spec, y, theta_anchor = 0)
  expect_equal(stan_data$W_type_id, 2L)
  expect_equal(stan_data$W_degree, 3L)
  expect_equal(stan_data$W_n_knots_full, 10L)
  expect_equal(stan_data$dim_W, 5L)

  fit <- mod$sample(
    data            = stan_data,
    chains          = 2L,
    parallel_chains = 2L,
    iter_warmup     = 300L,
    iter_sampling   = 300L,
    seed            = 81L,
    refresh         = 0L,
    show_messages   = FALSE,
    show_exceptions = FALSE
  )
  expect_s3_class(fit, "CmdStanMCMC")

  theta_post <- fit$summary("theta_ref")
  expect_lt(abs(theta_post$mean[1] - theta_ref_true), 0.3)
  sigma_post <- fit$summary("sigma_y")
  expect_lt(abs(sigma_post$mean[1] - sigma_true), 0.2)
})

test_that("K=2 het Gauss+Beta scale + bspline W recovers anchors", {
  testthat::skip_if(Sys.getenv("GDPAR_RUN_STAN_SMOKE_BSPLINE_HET") != "1",
                     "Set GDPAR_RUN_STAN_SMOKE_BSPLINE_HET=1 to run the K=2 het + bspline W smoke fit.")
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("cmdstanr")

  set.seed(82L)
  n <- 300L
  x <- rnorm(n); z <- rnorm(n)
  xc <- x - mean(x); zc <- z - mean(z)
  theta_mu_true    <- 0.3   # identity (Gaussian)
  theta_sigma_true <- 0.2   # logit (Beta) -> sigma ~ inv_logit(0.2) ~ 0.55
  beta_a_mu_true    <- 0.25
  beta_a_sigma_true <- 0.15
  eta_mu    <- theta_mu_true    + beta_a_mu_true    * xc
  eta_sigma <- theta_sigma_true + beta_a_sigma_true * zc
  sigma_vec <- 1 / (1 + exp(-eta_sigma))
  y <- rnorm(n, eta_mu, sigma_vec)
  d <- data.frame(y = y, x = x, z = z)

  fs <- gdpar_bf(y ~ a(x) + W(), sigma ~ a(z))
  fit <- suppressMessages(suppressWarnings(
    gdpar(fs, data = d,
          family = list(mu = gdpar_family("gaussian"),
                        sigma = gdpar_family("beta")),
          W = W_basis(type = "bspline", degree = 3, df = 5,
                      boundary_knots = c(-5, 5)),
          chains = 2L, iter_warmup = 400L, iter_sampling = 400L,
          refresh = 0L, seed = 41581L,
          skip_id_check = TRUE, verbose = FALSE)
  ))
  expect_s3_class(fit, "gdpar_fit")
  expect_equal(fit$K, 2L)
  expect_equal(fit$slot_names, c("mu", "sigma"))
  expect_equal(fit$stan_data$family_id_k, c(1L, 5L))
  expect_equal(fit$stan_data$inv_link_id_per_slot, c(0L, 1L))
  expect_equal(fit$stan_data$W_type_id, 2L)
  expect_equal(fit$stan_data$W_degree, 3L)
  expect_equal(fit$stan_data$W_n_knots_full, 10L)

  draws <- fit$fit$draws(format = "df")
  theta_mu_med    <- median(draws[["theta_ref_k[1,1]"]])
  theta_sigma_med <- median(draws[["theta_ref_k[1,2]"]])
  expect_lt(abs(theta_mu_med    - theta_mu_true),    0.5)
  expect_lt(abs(theta_sigma_med - theta_sigma_true), 0.7)
})

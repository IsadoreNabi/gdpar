# Sub-phase 8.3.9 (2026-05-22): unified bootstrap script for the
# 14 golden RDS snapshots called for by Decision D1=(1D) and the
# manifest CSV called for by Decisions D2=(2B) + E3=(E3.A).
#
# Run with:
#   NOT_CRAN=true GDPAR_BOOTSTRAP_8_3_9=1 \
#     /usr/bin/Rscript inst/extdata/scripts/bootstrap_8_3_9_goldens.R
#
# Each config is a self-contained closure that returns
#   list(rds_data, manifest_row, fit_code_hash)
# where rds_data is the named list saved to .rds (draws_matrix + meta).
#
# Configs covered (14):
#   1. lognormal_loc_scale K=2 (8.3.4 atrasado)
#   2. student_t K=3          (8.3.5a)
#   3. tweedie K=3            (8.3.5b)
#   4. zip K=2                (8.3.6)
#   5. zinb K=3               (8.3.6)
#   6. hurdle_poisson K=2     (8.3.6)
#   7. hurdle_neg_binomial_2 K=3 (8.3.6)
#   8. het_gauss_beta K=2     (8.3.7)
#   9. het_gauss_gamma K=2    (8.3.7)
#  10. het_nb_beta K=2        (8.3.7)
#  11. gaussian K=1 bspline p=1 (8.3.8)
#  12. het_gauss_beta K=2 bspline (8.3.8)
#  13. gaussian K=1 bspline p=2 (8.3.8)
#  14. gaussian K=1 polynomial p=2 baseline (8.3.9 new)

.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages(pkgload::load_all(".", quiet = TRUE))
stopifnot(requireNamespace("cmdstanr", quietly = TRUE))
stopifnot(requireNamespace("posterior", quietly = TRUE))

# ---- Common parameters --------------------------------------------------

CHAINS         <- 1L
ITER_WARMUP    <- 300L
ITER_SAMPLING  <- 300L
REFRESH        <- 0L
DATA_DIR       <- file.path("tests", "testthat", "data")

# Sub-phase 8.3.10 (2026-05-22): emit the startup banner only when
# the script is executed top-level (Rscript) with GDPAR_BOOTSTRAP_8_3_9
# set. When `source()`-ed from a test file that wants to extract the
# closures + configs list (e.g., the bit-exact compare-path tests),
# stays silent and skips the full run. Closures themselves are
# unchanged so that the SHA256 fit_code_hash recorded in the manifest
# remains bit-exact stable.
.is_bootstrap_main_run <- function() {
  sys.nframe() == 0L && nzchar(Sys.getenv("GDPAR_BOOTSTRAP_8_3_9"))
}
if (.is_bootstrap_main_run()) {
  if (!dir.exists(DATA_DIR)) {
    stop(sprintf("Data directory not found: %s", DATA_DIR))
  }
  cat("Bootstrap 8.3.9 starting at ", format(Sys.time()), "\n", sep = "")
  cat("cmdstan version: ", .gdpar_cmdstan_version(), "\n", sep = "")
  cat("DHARMa version : ", .gdpar_dharma_version(), "\n", sep = "")
  cat("DATA_DIR       : ", DATA_DIR, "\n", sep = "")
}

# ---- Helper: run one config, save RDS, add manifest row -----------------

run_config <- function(label, fit_code_fn, family, K, basis_type, p,
                       sub_phase, seed) {
  cat(sprintf("[%s] starting (seed=%d) ... ",
              label, seed))
  t0 <- Sys.time()
  fit <- tryCatch(fit_code_fn(),
                  error = function(e) {
                    cat("ERROR: ", conditionMessage(e), "\n", sep = "")
                    return(NULL)
                  })
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (is.null(fit)) {
    cat(sprintf("FAILED after %.1fs\n", elapsed))
    return(invisible(NULL))
  }
  draws_mat <- suppressWarnings(
    posterior::as_draws_matrix(fit$fit$draws())
  )
  n_obs <- if (!is.null(fit$stan_data$n)) fit$stan_data$n else NA_integer_
  rds_name <- .gdpar_golden_rds_name(family, K, basis_type, p)
  rds_path <- file.path(DATA_DIR, rds_name)
  rds_data <- list(
    draws_matrix    = draws_mat,
    family          = family,
    K               = K,
    basis_type      = basis_type,
    p               = p,
    sub_phase       = sub_phase,
    seed            = seed,
    bootstrap_date  = as.character(Sys.Date()),
    cmdstan_version = .gdpar_cmdstan_version(),
    R_version       = R.version.string,
    n_chains        = CHAINS,
    n_iter_warmup   = ITER_WARMUP,
    n_iter_sampling = ITER_SAMPLING,
    n_obs           = n_obs
  )
  saveRDS(rds_data, rds_path)
  .gdpar_golden_manifest_add(list(
    family          = family,
    K               = K,
    basis_type      = basis_type,
    p               = p,
    sub_phase       = sub_phase,
    bootstrap_date  = as.character(Sys.Date()),
    cmdstan_version = .gdpar_cmdstan_version(),
    R_version       = R.version.string,
    DHARMa_version  = .gdpar_dharma_version(),
    fit_code_hash   = .gdpar_golden_fit_code_hash(fit_code_fn),
    seed            = seed,
    n_chains        = CHAINS,
    n_iter_warmup   = ITER_WARMUP,
    n_iter_sampling = ITER_SAMPLING,
    rds_path        = rds_path,
    n_obs           = n_obs
  ))
  cat(sprintf("OK %.1fs -> %s\n", elapsed, rds_path))
  invisible(rds_path)
}

# ---- 14 Config closures --------------------------------------------------

cfg_01_lognormal_K2 <- function() {
  set.seed(45L)
  n <- 250L
  x <- rnorm(n); z <- rnorm(n)
  xc <- x - mean(x); zc <- z - mean(z)
  eta_mu    <- 0.2 + 0.4 * xc
  eta_sigma <- -0.3 + 0.2 * zc
  y <- rlnorm(n, eta_mu, exp(eta_sigma))
  d <- data.frame(y = y, x = x, z = z)
  fam_ln <- suppressMessages(gdpar_family_custom_K(
    name = "smoke_ln", stan_lpdf_id = "lognormal_loc_scale",
    did_holds = TRUE, did_reference = "test"
  ))
  fs <- gdpar_bf(y ~ a(x), sigma ~ a(z))
  suppressMessages(suppressWarnings(
    gdpar(fs, data = d, family = fam_ln,
          chains = CHAINS, iter_warmup = ITER_WARMUP,
          iter_sampling = ITER_SAMPLING, refresh = REFRESH,
          seed = 31418L, skip_id_check = TRUE, verbose = FALSE)
  ))
}

cfg_02_student_t_K3 <- function() {
  set.seed(46L)
  n <- 300L
  x <- rnorm(n); z <- rnorm(n); w <- rnorm(n)
  xc <- x - mean(x); zc <- z - mean(z); wc <- w - mean(w)
  eta_mu    <- 0.2 + 0.4 * xc
  eta_sigma <- 0.0 + 0.2 * zc
  eta_nu    <- log(15) + 0.1 * wc
  y <- eta_mu + exp(eta_sigma) * rt(n, df = exp(eta_nu))
  d <- data.frame(y = y, x = x, z = z, w = w)
  fs <- gdpar_bf(y ~ a(x), sigma ~ a(z), nu ~ a(w))
  suppressMessages(suppressWarnings(
    gdpar(fs, data = d, family = gdpar_family("student_t"),
          chains = CHAINS, iter_warmup = ITER_WARMUP,
          iter_sampling = ITER_SAMPLING, refresh = REFRESH,
          seed = 31419L, skip_id_check = TRUE, verbose = FALSE)
  ))
}

cfg_03_tweedie_K3 <- function() {
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
  eta_mu  <- 1.0 + 0.3 * xc
  eta_phi <- 0.0 + 0.2 * zc
  eta_p   <- 1.5 + 0.05 * wc
  p_vec   <- pmin(pmax(eta_p, 1.05), 1.95)
  y <- rtweedie(n, mu = exp(eta_mu), phi = exp(eta_phi), p = p_vec)
  d <- data.frame(y = y, x = x, z = z, w = w)
  fs <- gdpar_bf(y ~ a(x), phi ~ a(z), p ~ a(w))
  suppressMessages(suppressWarnings(
    gdpar(fs, data = d, family = gdpar_family("tweedie"),
          chains = CHAINS, iter_warmup = ITER_WARMUP,
          iter_sampling = ITER_SAMPLING, refresh = REFRESH,
          seed = 31420L, skip_id_check = TRUE, verbose = FALSE)
  ))
}

cfg_04_zip_K2 <- function() {
  set.seed(51L)
  n <- 300L
  x <- rnorm(n); z <- rnorm(n)
  xc <- x - mean(x); zc <- z - mean(z)
  eta_mu <- 1.2 + 0.4 * xc
  eta_pi <- -1.0 + 0.3 * zc
  pi_vec <- 1 / (1 + exp(-eta_pi))
  is_zero <- rbinom(n, size = 1L, prob = pi_vec) == 1L
  y <- ifelse(is_zero, 0L, rpois(n, lambda = exp(eta_mu)))
  d <- data.frame(y = y, x = x, z = z)
  fs <- gdpar_bf(y ~ a(x), pi ~ a(z))
  suppressMessages(suppressWarnings(
    gdpar(fs, data = d, family = gdpar_family("zip"),
          chains = CHAINS, iter_warmup = ITER_WARMUP,
          iter_sampling = ITER_SAMPLING, refresh = REFRESH,
          seed = 41510L, skip_id_check = TRUE, verbose = FALSE)
  ))
}

cfg_05_zinb_K3 <- function() {
  set.seed(52L)
  n <- 300L
  x <- rnorm(n); z <- rnorm(n); w <- rnorm(n)
  xc <- x - mean(x); zc <- z - mean(z); wc <- w - mean(w)
  eta_mu  <- 1.0 + 0.3 * xc
  eta_phi <- 0.5 + 0.2 * zc
  eta_pi  <- -1.2 + 0.3 * wc
  pi_vec  <- 1 / (1 + exp(-eta_pi))
  is_zero <- rbinom(n, size = 1L, prob = pi_vec) == 1L
  y <- ifelse(is_zero, 0L,
              rnbinom(n, mu = exp(eta_mu), size = exp(eta_phi)))
  d <- data.frame(y = y, x = x, z = z, w = w)
  fs <- gdpar_bf(y ~ a(x), phi ~ a(z), pi ~ a(w))
  suppressMessages(suppressWarnings(
    gdpar(fs, data = d, family = gdpar_family("zinb"),
          chains = CHAINS, iter_warmup = ITER_WARMUP,
          iter_sampling = ITER_SAMPLING, refresh = REFRESH,
          seed = 41511L, skip_id_check = TRUE, verbose = FALSE)
  ))
}

cfg_06_hurdle_p_K2 <- function() {
  rtpois <- function(n, mu) {
    out <- integer(n)
    for (i in seq_len(n)) {
      val <- 0L; iter <- 0L
      while (val == 0L && iter < 10000L) {
        val <- rpois(1L, lambda = mu[i]); iter <- iter + 1L
      }
      out[i] <- val
    }
    out
  }
  set.seed(53L)
  n <- 300L
  x <- rnorm(n); z <- rnorm(n)
  xc <- x - mean(x); zc <- z - mean(z)
  eta_mu <- 1.0 + 0.4 * xc
  eta_pi <- -0.5 + 0.3 * zc
  mu_vec <- exp(eta_mu); pi_vec <- 1 / (1 + exp(-eta_pi))
  is_zero <- rbinom(n, size = 1L, prob = pi_vec) == 1L
  y <- integer(n); y[is_zero] <- 0L
  y[!is_zero] <- rtpois(sum(!is_zero), mu = mu_vec[!is_zero])
  d <- data.frame(y = y, x = x, z = z)
  fs <- gdpar_bf(y ~ a(x), pi ~ a(z))
  suppressMessages(suppressWarnings(
    gdpar(fs, data = d, family = gdpar_family("hurdle_poisson"),
          chains = CHAINS, iter_warmup = ITER_WARMUP,
          iter_sampling = ITER_SAMPLING, refresh = REFRESH,
          seed = 41512L, skip_id_check = TRUE, verbose = FALSE)
  ))
}

cfg_07_hurdle_nb_K3 <- function() {
  rtnbinom <- function(n, mu, size) {
    out <- integer(n)
    for (i in seq_len(n)) {
      val <- 0L; iter <- 0L
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
  eta_mu  <- 1.0 + 0.3 * xc
  eta_phi <- 0.5 + 0.2 * zc
  eta_pi  <- -0.5 + 0.3 * wc
  pi_vec  <- 1 / (1 + exp(-eta_pi))
  is_zero <- rbinom(n, size = 1L, prob = pi_vec) == 1L
  y <- integer(n); y[is_zero] <- 0L
  y[!is_zero] <- rtnbinom(sum(!is_zero),
                          mu   = exp(eta_mu[!is_zero]),
                          size = exp(eta_phi[!is_zero]))
  d <- data.frame(y = y, x = x, z = z, w = w)
  fs <- gdpar_bf(y ~ a(x), phi ~ a(z), pi ~ a(w))
  suppressMessages(suppressWarnings(
    gdpar(fs, data = d,
          family = gdpar_family("hurdle_neg_binomial_2"),
          chains = CHAINS, iter_warmup = ITER_WARMUP,
          iter_sampling = ITER_SAMPLING, refresh = REFRESH,
          seed = 41513L, skip_id_check = TRUE, verbose = FALSE)
  ))
}

cfg_08_het_gauss_beta_K2 <- function() {
  set.seed(61L)
  n <- 300L
  x <- rnorm(n); z <- rnorm(n)
  xc <- x - mean(x); zc <- z - mean(z)
  eta_mu    <- 0.4 + 0.4 * xc
  eta_sigma <- -0.8 + 0.3 * zc
  sigma_vec <- 1 / (1 + exp(-eta_sigma))
  y <- rnorm(n, eta_mu, sigma_vec)
  d <- data.frame(y = y, x = x, z = z)
  fs <- gdpar_bf(y ~ a(x), sigma ~ a(z))
  suppressMessages(suppressWarnings(
    gdpar(fs, data = d,
          family = list(mu = gdpar_family("gaussian"),
                        sigma = gdpar_family("beta")),
          chains = CHAINS, iter_warmup = ITER_WARMUP,
          iter_sampling = ITER_SAMPLING, refresh = REFRESH,
          seed = 41521L, skip_id_check = TRUE, verbose = FALSE)
  ))
}

cfg_09_het_gauss_gamma_K2 <- function() {
  set.seed(62L)
  n <- 300L
  x <- rnorm(n); z <- rnorm(n)
  xc <- x - mean(x); zc <- z - mean(z)
  eta_mu    <- 0.5 + 0.4 * xc
  eta_sigma <- -0.3 + 0.2 * zc
  sigma_vec <- exp(eta_sigma)
  y <- rnorm(n, eta_mu, sigma_vec)
  d <- data.frame(y = y, x = x, z = z)
  fs <- gdpar_bf(y ~ a(x), sigma ~ a(z))
  suppressMessages(suppressWarnings(
    gdpar(fs, data = d,
          family = list(mu = gdpar_family("gaussian"),
                        sigma = gdpar_family("gamma")),
          chains = CHAINS, iter_warmup = ITER_WARMUP,
          iter_sampling = ITER_SAMPLING, refresh = REFRESH,
          seed = 41522L, skip_id_check = TRUE, verbose = FALSE)
  ))
}

cfg_10_het_nb_beta_K2 <- function() {
  set.seed(63L)
  n <- 300L
  x <- rnorm(n); z <- rnorm(n)
  xc <- x - mean(x); zc <- z - mean(z)
  eta_mu  <- 0.8 + 0.3 * xc
  eta_phi <- 0.3 + 0.2 * zc
  mu_vec  <- exp(eta_mu)
  phi_vec <- 1 / (1 + exp(-eta_phi))
  y <- rnbinom(n, mu = mu_vec, size = phi_vec)
  d <- data.frame(y = y, x = x, z = z)
  fs <- gdpar_bf(y ~ a(x), phi ~ a(z))
  suppressMessages(suppressWarnings(
    gdpar(fs, data = d,
          family = list(mu = gdpar_family("neg_binomial_2"),
                        phi = gdpar_family("beta")),
          chains = CHAINS, iter_warmup = ITER_WARMUP,
          iter_sampling = ITER_SAMPLING, refresh = REFRESH,
          seed = 41523L, skip_id_check = TRUE, verbose = FALSE)
  ))
}

cfg_11_gaussian_K1_bspline_p1 <- function() {
  # Mirror of test-fit_K_smoke.R BSPLINE_K1 smoke (low-level path
  # through assemble_stan_data + amm_main.stan; brms-style formula
  # path is incompatible with explicit 'amm' arg).
  set.seed(81L)
  n <- 250L
  x1 <- rnorm(n); x2 <- rnorm(n)
  x1c <- x1 - mean(x1)
  y <- rnorm(n, 0.4 + 0.3 * x1c, 0.4)
  data_df <- data.frame(x1 = x1, x2 = x2)
  spec <- amm_spec(
    a = ~ x1,
    W = W_basis(type = "bspline", degree = 3, df = 5,
                boundary_knots = c(-10, 10))
  )
  spec$W <- materialize_W_basis(spec$W, 1L)
  fam   <- gdpar_family("gaussian")
  prior <- gdpar_prior()
  src <- generate_stan_code(prior, cp_a = FALSE, cp_W = FALSE)
  stan_file <- tempfile(fileext = ".stan")
  writeLines(src, stan_file)
  mod <- cmdstanr::cmdstan_model(stan_file)
  design <- build_amm_design(spec, data_df, formula_rhs = ~ x2)
  stan_data <- assemble_stan_data(design, fam, spec, y, theta_anchor = 0)
  csm <- mod$sample(
    data = stan_data, chains = CHAINS,
    iter_warmup = ITER_WARMUP, iter_sampling = ITER_SAMPLING,
    refresh = REFRESH, seed = 81L,
    show_messages = FALSE, show_exceptions = FALSE
  )
  list(fit = csm, stan_data = stan_data, amm = spec, family = fam)
}

cfg_12_het_gauss_beta_bspline_K2 <- function() {
  set.seed(82L)
  n <- 300L
  x <- rnorm(n); z <- rnorm(n)
  xc <- x - mean(x); zc <- z - mean(z)
  eta_mu    <- 0.3 + 0.25 * xc
  eta_sigma <- 0.2 + 0.15 * zc
  sigma_vec <- 1 / (1 + exp(-eta_sigma))
  y <- rnorm(n, eta_mu, sigma_vec)
  d <- data.frame(y = y, x = x, z = z)
  fs <- gdpar_bf(y ~ a(x) + W(), sigma ~ a(z))
  suppressMessages(suppressWarnings(
    gdpar(fs, data = d,
          family = list(mu = gdpar_family("gaussian"),
                        sigma = gdpar_family("beta")),
          W = W_basis(type = "bspline", degree = 3, df = 5,
                      boundary_knots = c(-5, 5)),
          chains = CHAINS, iter_warmup = ITER_WARMUP,
          iter_sampling = ITER_SAMPLING, refresh = REFRESH,
          seed = 41581L, skip_id_check = TRUE, verbose = FALSE)
  ))
}

cfg_13_gaussian_K1_bspline_p2 <- function() {
  set.seed(83L)
  n <- 300L; p <- 2L
  x1 <- rnorm(n); x2 <- rnorm(n)
  x1c <- x1 - mean(x1); x2c <- x2 - mean(x2)
  theta_ref_true <- c(0.3, -0.4); sigma_true <- 0.4
  beta_a_true <- list(c(0.4, -0.3), c(-0.2, 0.3))
  eta_true <- cbind(
    theta_ref_true[1] + beta_a_true[[1]][1] * x1c + beta_a_true[[1]][2] * x2c,
    theta_ref_true[2] + beta_a_true[[2]][1] * x1c + beta_a_true[[2]][2] * x2c
  )
  y <- eta_true + matrix(rnorm(n * p, 0, sigma_true), n, p)
  d <- data.frame(x1 = x1, x2 = x2)
  spec <- amm_spec(
    p = 2L,
    dims = dimwise(a = ~ x1 + x2, b = NULL),
    W = W_basis(type = "bspline", degree = 3, df = 5,
                boundary_knots = c(-5, 5), p = 2L)
  )
  fam   <- gdpar_family_multi("gaussian", p = 2L)
  prior <- gdpar_prior()
  src <- generate_stan_code_multi(prior, cp_a = FALSE, cp_W = FALSE)
  stan_file <- tempfile(fileext = ".stan")
  writeLines(src, stan_file)
  mod <- cmdstanr::cmdstan_model(stan_file)
  design <- build_amm_design(spec, d, formula_rhs = ~ x1 + x2)
  stan_data <- assemble_stan_data(design, fam, spec, y, theta_anchor = 0)
  csm <- mod$sample(
    data = stan_data, chains = CHAINS,
    iter_warmup = ITER_WARMUP, iter_sampling = ITER_SAMPLING,
    refresh = REFRESH, seed = 83L,
    show_messages = FALSE, show_exceptions = FALSE
  )
  # Wrap to a gdpar_fit-like object so the manifest helper can read $stan_data
  list(fit = csm, stan_data = stan_data, amm = spec, family = fam, p = 2L)
}

cfg_14_gaussian_K1_polynomial_p2 <- function() {
  set.seed(84L)
  n <- 250L; p <- 2L
  x1 <- rnorm(n); x2 <- rnorm(n)
  x1c <- x1 - mean(x1); x2c <- x2 - mean(x2)
  theta_ref_true <- c(0.2, -0.3); sigma_true <- 0.35
  beta_a_true <- list(c(0.4, -0.2), c(-0.3, 0.25))
  eta_true <- cbind(
    theta_ref_true[1] + beta_a_true[[1]][1] * x1c + beta_a_true[[1]][2] * x2c,
    theta_ref_true[2] + beta_a_true[[2]][1] * x1c + beta_a_true[[2]][2] * x2c
  )
  y <- eta_true + matrix(rnorm(n * p, 0, sigma_true), n, p)
  d <- data.frame(x1 = x1, x2 = x2)
  spec <- amm_spec(p = 2L, dims = dimwise(a = ~ x1 + x2, b = NULL))
  fam   <- gdpar_family_multi("gaussian", p = 2L)
  prior <- gdpar_prior()
  src <- generate_stan_code_multi(prior, cp_a = FALSE, cp_W = FALSE)
  stan_file <- tempfile(fileext = ".stan")
  writeLines(src, stan_file)
  mod <- cmdstanr::cmdstan_model(stan_file)
  design <- build_amm_design(spec, d, formula_rhs = ~ x1 + x2)
  stan_data <- assemble_stan_data(design, fam, spec, y, theta_anchor = 0)
  csm <- mod$sample(
    data = stan_data, chains = CHAINS,
    iter_warmup = ITER_WARMUP, iter_sampling = ITER_SAMPLING,
    refresh = REFRESH, seed = 84L,
    show_messages = FALSE, show_exceptions = FALSE
  )
  list(fit = csm, stan_data = stan_data, amm = spec, family = fam, p = 2L)
}

# ---- Run all 14 ---------------------------------------------------------

configs <- list(
  list("01_lognormal_K2",          cfg_01_lognormal_K2,
       "lognormal_loc_scale", 2L, "polynomial", 1L, "8.3.4", 31418L),
  list("02_student_t_K3",          cfg_02_student_t_K3,
       "student_t",           3L, "polynomial", 1L, "8.3.5a", 31419L),
  list("03_tweedie_K3",            cfg_03_tweedie_K3,
       "tweedie",             3L, "polynomial", 1L, "8.3.5b", 31420L),
  list("04_zip_K2",                cfg_04_zip_K2,
       "zip",                 2L, "polynomial", 1L, "8.3.6", 41510L),
  list("05_zinb_K3",               cfg_05_zinb_K3,
       "zinb",                3L, "polynomial", 1L, "8.3.6", 41511L),
  list("06_hurdle_p_K2",           cfg_06_hurdle_p_K2,
       "hurdle_poisson",      2L, "polynomial", 1L, "8.3.6", 41512L),
  list("07_hurdle_nb_K3",          cfg_07_hurdle_nb_K3,
       "hurdle_neg_binomial_2", 3L, "polynomial", 1L, "8.3.6", 41513L),
  list("08_het_gauss_beta_K2",     cfg_08_het_gauss_beta_K2,
       "het_gauss_beta",      2L, "polynomial", 1L, "8.3.7", 41521L),
  list("09_het_gauss_gamma_K2",    cfg_09_het_gauss_gamma_K2,
       "het_gauss_gamma",     2L, "polynomial", 1L, "8.3.7", 41522L),
  list("10_het_nb_beta_K2",        cfg_10_het_nb_beta_K2,
       "het_nb_beta",         2L, "polynomial", 1L, "8.3.7", 41523L),
  list("11_gaussian_K1_bspline_p1", cfg_11_gaussian_K1_bspline_p1,
       "gaussian",            1L, "bspline",    1L, "8.3.8", 81L),
  list("12_het_gauss_beta_bspline_K2", cfg_12_het_gauss_beta_bspline_K2,
       "het_gauss_beta",      2L, "bspline",    1L, "8.3.8", 41581L),
  list("13_gaussian_K1_bspline_p2", cfg_13_gaussian_K1_bspline_p2,
       "gaussian",            1L, "bspline",    2L, "8.3.8", 83L),
  list("14_gaussian_K1_polynomial_p2", cfg_14_gaussian_K1_polynomial_p2,
       "gaussian",            1L, "polynomial", 2L, "8.3.9", 84L)
)

if (.is_bootstrap_main_run()) {
  skip_set <- strsplit(Sys.getenv("GDPAR_BOOTSTRAP_8_3_9_SKIP", ""), ",")[[1L]]
  skip_set <- trimws(skip_set)

  results <- list()
  for (cfg in configs) {
    label <- cfg[[1L]]
    if (label %in% skip_set) {
      cat(sprintf("[%s] SKIPPED by GDPAR_BOOTSTRAP_8_3_9_SKIP\n", label))
      next
    }
    rds_path <- run_config(
      label       = label,
      fit_code_fn = cfg[[2L]],
      family      = cfg[[3L]],
      K           = cfg[[4L]],
      basis_type  = cfg[[5L]],
      p           = cfg[[6L]],
      sub_phase   = cfg[[7L]],
      seed        = cfg[[8L]]
    )
    results[[label]] <- rds_path
  }

  cat("\nBootstrap 8.3.9 completed at ", format(Sys.time()), "\n", sep = "")
  cat("Generated ", sum(!vapply(results, is.null, logical(1L))), " RDS files\n",
      sep = "")
  cat("Manifest at: ",
      .gdpar_golden_manifest_path(), "\n", sep = "")
}

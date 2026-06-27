test_that("F.6 smoke fit p=2 gaussian homogeneous compiles, samples, and recovers theta_ref", {
  skip_if(Sys.getenv("GDPAR_F6_FIT") != "1",
          "Set GDPAR_F6_FIT=1 to run the F.6 smoke fit (compiles a Stan model; takes minutes).")
  skip_on_cran()
  skip_if_not_installed("cmdstanr")

  n <- 300L
  p <- 2L
  set.seed(42L)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x1_c <- x1 - mean(x1)
  x2_c <- x2 - mean(x2)
  data_df <- data.frame(x1 = x1, x2 = x2)

  theta_ref_true <- c(0.5, -0.5)
  beta_a_true <- list(c(0.5, -0.5), c(-0.4, 0.4))
  sigma_true <- 0.3
  eta_true <- cbind(
    theta_ref_true[1] + beta_a_true[[1]][1] * x1_c + beta_a_true[[1]][2] * x2_c,
    theta_ref_true[2] + beta_a_true[[2]][1] * x1_c + beta_a_true[[2]][2] * x2_c
  )
  y <- eta_true + matrix(rnorm(n * p, 0, sigma_true), n, p)

  spec  <- amm_spec(p = 2L, dims = dimwise(a = ~ x1 + x2, b = NULL))
  fam   <- gdpar_family_multi("gaussian", p = 2L)
  prior <- gdpar_prior()

  src <- generate_stan_code_multi(prior, cp_a = FALSE, cp_W = FALSE)
  stan_file <- write_stan_to_tempfile(src)
  mod <- cmdstanr::cmdstan_model(stan_file)
  expect_s3_class(mod, "CmdStanModel")

  design <- build_amm_design(spec, data_df, formula_rhs = ~ x1 + x2)
  stan_data <- assemble_stan_data(design, fam, spec, y, theta_anchor = 0)
  expect_equal(stan_data$p, 2L)
  expect_equal(stan_data$family_id, 1L)
  expect_equal(dim(stan_data$y_real), c(n, p))

  fit <- mod$sample(
    data            = stan_data,
    chains          = 2L,
    parallel_chains = 2L,
    iter_warmup     = 200L,
    iter_sampling   = 200L,
    seed            = 42L,
    refresh         = 0L,
    show_messages   = FALSE,
    show_exceptions = FALSE
  )
  expect_s3_class(fit, "CmdStanMCMC")

  diag_summary <- fit$summary(c("theta_ref", "sigma_y"))
  expect_true(all(diag_summary$rhat < 1.10),
              info = paste("max Rhat:",
                           max(diag_summary$rhat, na.rm = TRUE)))
  expect_true(all(diag_summary$ess_bulk > 100),
              info = paste("min ESS bulk:",
                           min(diag_summary$ess_bulk, na.rm = TRUE)))

  theta_post <- fit$summary("theta_ref")
  expect_lt(abs(theta_post$mean[1] - theta_ref_true[1]), 0.2)
  expect_lt(abs(theta_post$mean[2] - theta_ref_true[2]), 0.2)

  sigma_post <- fit$summary("sigma_y")
  expect_lt(abs(sigma_post$mean[1] - sigma_true), 0.1)
  expect_lt(abs(sigma_post$mean[2] - sigma_true), 0.1)
})

# Sub-phase 8.3.8 (2026-05-22): smoke fit exercising the B-spline W
# basis under the multivariate K = 1 + p = 2 regime that routes through
# amm_distrib_multi.stan. The basis is separable (one block per
# coordinate of theta_ref); dim_W = df * p = 5 * 2 = 10.

test_that("F.6 smoke fit p=2 Gaussian + bspline W compiles and recovers theta_ref", {
  skip_if(Sys.getenv("GDPAR_RUN_STAN_SMOKE_BSPLINE_K1_P2") != "1",
          "Set GDPAR_RUN_STAN_SMOKE_BSPLINE_K1_P2=1 to run the multi-coord bspline smoke fit.")
  skip_on_cran()
  skip_if_not_installed("cmdstanr")

  n <- 300L
  p <- 2L
  set.seed(83L)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x1_c <- x1 - mean(x1)
  x2_c <- x2 - mean(x2)
  data_df <- data.frame(x1 = x1, x2 = x2)

  theta_ref_true <- c(0.3, -0.4)
  beta_a_true <- list(c(0.4, -0.3), c(-0.2, 0.3))
  sigma_true <- 0.4
  eta_true <- cbind(
    theta_ref_true[1] + beta_a_true[[1]][1] * x1_c + beta_a_true[[1]][2] * x2_c,
    theta_ref_true[2] + beta_a_true[[2]][1] * x1_c + beta_a_true[[2]][2] * x2_c
  )
  y <- eta_true + matrix(rnorm(n * p, 0, sigma_true), n, p)

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
  expect_s3_class(mod, "CmdStanModel")

  design <- build_amm_design(spec, data_df, formula_rhs = ~ x1 + x2)
  stan_data <- assemble_stan_data(design, fam, spec, y, theta_anchor = 0)
  expect_equal(stan_data$p, 2L)
  expect_equal(stan_data$W_type_id, 2L)
  expect_equal(stan_data$W_degree, 3L)
  expect_equal(stan_data$dim_W, 10L)
  expect_equal(stan_data$W_per_k_dim, 5L)

  fit <- mod$sample(
    data            = stan_data,
    chains          = 2L,
    parallel_chains = 2L,
    iter_warmup     = 300L,
    iter_sampling   = 300L,
    seed            = 83L,
    refresh         = 0L,
    show_messages   = FALSE,
    show_exceptions = FALSE
  )
  expect_s3_class(fit, "CmdStanMCMC")

  theta_post <- fit$summary("theta_ref")
  expect_lt(abs(theta_post$mean[1] - theta_ref_true[1]), 0.3)
  expect_lt(abs(theta_post$mean[2] - theta_ref_true[2]), 0.3)
})

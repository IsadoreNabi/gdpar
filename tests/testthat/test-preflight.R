test_that("preflight_attribution_score is small under H0 (random divergences)", {
  set.seed(42)
  n <- 400L
  log_sigma <- rnorm(n)
  raw_norm <- abs(rnorm(n, mean = 1, sd = 0.3))
  divergent <- rep(FALSE, n)
  divergent[sample.int(n, 40L)] <- TRUE
  score <- gdpar:::preflight_attribution_score(
    log_sigma, raw_norm, divergent
  )
  expect_false(is.na(score))
  expect_lt(abs(score), 3)
})

test_that("preflight_attribution_score is large positive under funnel concentration", {
  set.seed(42)
  n_nondiv <- 380L
  n_div <- 20L
  log_sigma <- c(rnorm(n_nondiv, mean = 0, sd = 1),
                 rnorm(n_div, mean = -2, sd = 0.3))
  raw_norm <- c(abs(rnorm(n_nondiv, mean = 1, sd = 0.4)),
                abs(rnorm(n_div, mean = 3, sd = 0.3)))
  divergent <- c(rep(FALSE, n_nondiv), rep(TRUE, n_div))
  score <- gdpar:::preflight_attribution_score(
    log_sigma, raw_norm, divergent
  )
  expect_false(is.na(score))
  expect_gt(score, 5)
})

test_that("preflight_attribution_score is negative when divergent set is positively correlated", {
  set.seed(42)
  n_nondiv <- 380L
  n_div <- 20L
  log_sigma <- c(rnorm(n_nondiv, mean = 0, sd = 1),
                 rnorm(n_div, mean = -2, sd = 0.3))
  raw_norm <- c(abs(rnorm(n_nondiv, mean = 2, sd = 0.5)),
                abs(rnorm(n_div, mean = 0.3, sd = 0.05)))
  divergent <- c(rep(FALSE, n_nondiv), rep(TRUE, n_div))
  score <- gdpar:::preflight_attribution_score(
    log_sigma, raw_norm, divergent
  )
  expect_false(is.na(score))
  expect_lt(score, -2)
})

test_that("preflight_attribution_score returns NA on degenerate input", {
  expect_true(is.na(gdpar:::preflight_attribution_score(
    log_sigma = rep(1, 100), raw_norm = runif(100, 1, 2),
    divergent = c(rep(TRUE, 10), rep(FALSE, 90))
  )))
  expect_true(is.na(gdpar:::preflight_attribution_score(
    log_sigma = rnorm(100), raw_norm = rnorm(50),
    divergent = rep(c(TRUE, FALSE), 50)
  )))
  expect_true(is.na(gdpar:::preflight_attribution_score(
    log_sigma = rnorm(100), raw_norm = abs(rnorm(100)),
    divergent = rep(FALSE, 100)
  )))
})

test_that("block_bootstrap_indices respects per-chain block structure", {
  set.seed(123)
  idx <- gdpar:::block_bootstrap_indices(
    n_draws = 400L, n_chains = 2L, block_size = 10L
  )
  expect_length(idx, 400L)
  expect_true(all(idx >= 1L & idx <= 400L))
  expect_true(all(idx[1:200] <= 200L))
  expect_true(all(idx[201:400] > 200L))
  set.seed(123)
  idx2 <- gdpar:::block_bootstrap_indices(
    n_draws = 400L, n_chains = 2L, block_size = 10L
  )
  expect_identical(idx, idx2)
})

test_that("block_bootstrap_indices errors on invalid arguments", {
  # D19 (Sesion 17 Etapa 2 de Sesion 8.4): the two guards in
  # block_bootstrap_indices now signal class = "gdpar_input_error" via
  # gdpar_abort() (previously plain stop()).
  expect_error(
    gdpar:::block_bootstrap_indices(
      n_draws = 401L, n_chains = 2L, block_size = 10L
    ),
    class = "gdpar_input_error"
  )
  expect_error(
    gdpar:::block_bootstrap_indices(
      n_draws = 400L, n_chains = 2L, block_size = 0L
    ),
    class = "gdpar_input_error"
  )
  expect_error(
    gdpar:::block_bootstrap_indices(
      n_draws = 400L, n_chains = 2L, block_size = 300L
    ),
    class = "gdpar_input_error"
  )
})

test_that("preflight_info_ratio_t flags high-info regime via t_cp (Path B')", {
  set.seed(42)
  n_draws <- 400L
  n_coords <- 6L
  effective_coef <- matrix(rnorm(n_draws * n_coords, sd = 0.05),
                           n_draws, n_coords)
  reference_scale_per_draw <- rep(1.0, n_draws)
  res <- gdpar:::preflight_info_ratio_t(
    effective_coef = effective_coef,
    reference_scale_per_draw = reference_scale_per_draw,
    n_chains = 2L,
    tau_cp = 10, tau_ncp = 3
  )
  expect_false(is.na(res$t_cp))
  expect_gt(res$t_cp, 1.96)
  expect_gt(res$t_ncp, res$t_cp)
})

test_that("preflight_info_ratio_t flags low-info regime via t_ncp (Path B')", {
  set.seed(42)
  n_draws <- 400L
  n_coords <- 6L
  effective_coef <- matrix(rnorm(n_draws * n_coords, sd = 0.95),
                           n_draws, n_coords)
  reference_scale_per_draw <- rep(1.0, n_draws)
  res <- gdpar:::preflight_info_ratio_t(
    effective_coef = effective_coef,
    reference_scale_per_draw = reference_scale_per_draw,
    n_chains = 2L,
    tau_cp = 10, tau_ncp = 3
  )
  expect_false(is.na(res$t_ncp))
  expect_lt(res$t_ncp, -1.96)
  expect_lt(res$t_cp, res$t_ncp)
})

test_that("preflight_info_ratio_t handles single-coord input (Path B')", {
  set.seed(42)
  n_draws <- 400L
  effective_coef <- matrix(rnorm(n_draws, sd = 0.05),
                           n_draws, 1L)
  reference_scale_per_draw <- rep(1.0, n_draws)
  res <- gdpar:::preflight_info_ratio_t(
    effective_coef = effective_coef,
    reference_scale_per_draw = reference_scale_per_draw,
    n_chains = 2L,
    tau_cp = 10, tau_ncp = 3
  )
  expect_false(is.na(res$t_cp))
  expect_false(is.na(res$t_ncp))
  expect_gt(res$t_cp, 1.96)
})

test_that("preflight_info_ratio_t returns NA on degenerate input (Path B')", {
  set.seed(42)
  n_draws <- 400L
  expect_true(is.na(gdpar:::preflight_info_ratio_t(
    effective_coef = matrix(0, n_draws, 5L),
    reference_scale_per_draw = rep(1.0, n_draws),
    n_chains = 2L
  )$t_cp))
  expect_true(is.na(gdpar:::preflight_info_ratio_t(
    effective_coef = matrix(rnorm(n_draws * 5L), n_draws, 5L),
    reference_scale_per_draw = rep(0.0, n_draws),
    n_chains = 2L
  )$t_cp))
  expect_true(is.na(gdpar:::preflight_info_ratio_t(
    effective_coef = matrix(rnorm((n_draws - 1L) * 5L),
                            n_draws - 1L, 5L),
    reference_scale_per_draw = rep(1.0, n_draws),
    n_chains = 2L
  )$t_cp))
  expect_true(is.na(gdpar:::preflight_info_ratio_t(
    effective_coef = matrix(rnorm(n_draws * 5L), n_draws, 5L),
    reference_scale_per_draw = c(NA_real_, rep(1.0, n_draws - 1L)),
    n_chains = 2L
  )$t_cp))
  expect_true(is.na(gdpar:::preflight_info_ratio_t(
    effective_coef = matrix(rnorm(401L * 5L), 401L, 5L),
    reference_scale_per_draw = rep(1.0, 401L),
    n_chains = 2L
  )$t_cp))
})

test_that("resolve_parametrization with global 'ncp' skips pre-flight", {
  res <- gdpar:::resolve_parametrization(
    parametrization = "ncp",
    parametrization_a = NULL, parametrization_W = NULL,
    prior = NULL, stan_data = NULL, amm = NULL,
    preflight_seed = NULL, verbose = FALSE
  )
  expect_false(res$cp_a)
  expect_false(res$cp_W)
  expect_false(res$meta$used_preflight)
  expect_equal(res$meta$decision_reason_a, "user_global")
  expect_equal(res$meta$decision_reason_W, "user_global")
})

test_that("resolve_parametrization with global 'cp' skips pre-flight", {
  res <- gdpar:::resolve_parametrization(
    parametrization = "cp",
    parametrization_a = NULL, parametrization_W = NULL,
    prior = NULL, stan_data = NULL, amm = NULL,
    preflight_seed = NULL, verbose = FALSE
  )
  expect_true(res$cp_a)
  expect_true(res$cp_W)
  expect_false(res$meta$used_preflight)
})

test_that("resolve_parametrization with both per-component explicit skips pre-flight", {
  res <- gdpar:::resolve_parametrization(
    parametrization = "auto",
    parametrization_a = "cp", parametrization_W = "ncp",
    prior = NULL, stan_data = NULL, amm = NULL,
    preflight_seed = NULL, verbose = FALSE
  )
  expect_true(res$cp_a)
  expect_false(res$cp_W)
  expect_false(res$meta$used_preflight)
  expect_equal(res$meta$decision_reason_a, "user_explicit_a")
  expect_equal(res$meta$decision_reason_W, "user_explicit_W")
})

test_that("resolve_parametrization mixed: global 'ncp' plus explicit 'cp' for a", {
  res <- gdpar:::resolve_parametrization(
    parametrization = "ncp",
    parametrization_a = "cp", parametrization_W = NULL,
    prior = NULL, stan_data = NULL, amm = NULL,
    preflight_seed = NULL, verbose = FALSE
  )
  expect_true(res$cp_a)
  expect_false(res$cp_W)
  expect_false(res$meta$used_preflight)
  expect_equal(res$meta$decision_reason_a, "user_explicit_a")
  expect_equal(res$meta$decision_reason_W, "user_global")
})

test_that("preflight picks CP for a in high-info regime (default thresholds)", {
  skip_if_no_cmdstan()
  set.seed(42)
  n <- 500L
  df <- data.frame(x1 = rnorm(n), x2 = rnorm(n),
                   x3 = rnorm(n), x4 = rnorm(n))
  df$y <- 1 + 1.0*df$x1 - 0.8*df$x2 + 0.6*df$x3 - 0.4*df$x4 +
          rnorm(n, sd = 0.1)
  spec <- amm_spec(a = ~ x1 + x2 + x3 + x4)
  design <- gdpar:::build_amm_design(
    spec, df, formula_rhs = ~ x1 + x2 + x3 + x4
  )
  fam <- gdpar_family("gaussian")
  stan_data <- gdpar:::assemble_stan_data(
    design, fam, spec, df$y, theta_anchor = 0
  )
  res <- gdpar:::preflight_parametrization(
    gdpar_prior(), stan_data, spec,
    preflight_seed = 42L, verbose = FALSE
  )
  expect_true(res$cp_a)
  expect_match(res$meta$decision_reason_a, "filter_info_high")
})

test_that("preflight picks NCP for a in low-info regime (default thresholds)", {
  skip_if_no_cmdstan()
  set.seed(42)
  n <- 40L
  df <- data.frame(x1 = rnorm(n), x2 = rnorm(n),
                   x3 = rnorm(n), x4 = rnorm(n))
  df$y <- 1 + 0.02*df$x1 - 0.02*df$x2 + 0.01*df$x3 - 0.01*df$x4 +
          rnorm(n, sd = 3.0)
  spec <- amm_spec(a = ~ x1 + x2 + x3 + x4)
  design <- gdpar:::build_amm_design(
    spec, df, formula_rhs = ~ x1 + x2 + x3 + x4
  )
  fam <- gdpar_family("gaussian")
  stan_data <- gdpar:::assemble_stan_data(
    design, fam, spec, df$y, theta_anchor = 0
  )
  res <- gdpar:::preflight_parametrization(
    gdpar_prior(), stan_data, spec,
    preflight_seed = 42L, verbose = FALSE
  )
  expect_false(res$cp_a)
  expect_match(res$meta$decision_reason_a,
               "filter_info_(low|ambiguous_ncp)")
})

test_that("preflight picks NCP for W in low-info regime (default thresholds)", {
  skip_if_no_cmdstan()
  set.seed(42)
  n <- 40L
  x1 <- rnorm(n); x2 <- rnorm(n)
  c1 <- 0.3 * 0.05 + 0.09 * (-0.03)
  c2 <- 0.3 * (-0.04) + 0.09 * 0.02
  y <- 0.3 + c1*x1 + c2*x2 + rnorm(n, sd = 2.0)
  df <- data.frame(x1 = x1, x2 = x2, y = y)
  spec <- amm_spec(
    W = W_basis(type = "polynomial", degree = 2),
    x_vars = c("x1", "x2")
  )
  spec$W <- gdpar:::materialize_W_basis(spec$W, p = 1L)
  design <- gdpar:::build_amm_design(
    spec, df, formula_rhs = ~ x1 + x2
  )
  fam <- gdpar_family("gaussian")
  stan_data <- gdpar:::assemble_stan_data(
    design, fam, spec, df$y, theta_anchor = 0
  )
  res <- gdpar:::preflight_parametrization(
    gdpar_prior(), stan_data, spec,
    preflight_seed = 42L, verbose = FALSE
  )
  expect_false(res$cp_W)
  expect_match(res$meta$decision_reason_W,
               "filter_info_(low|ambiguous_ncp)")
})

test_that("preflight picks CP for W in high-info regime (Path B', default thresholds)", {
  skip_if_no_cmdstan()
  set.seed(42)
  n <- 500L
  x1 <- rnorm(n); x2 <- rnorm(n)
  theta_ref_true <- 0.3
  c1 <- theta_ref_true * 2.0 + theta_ref_true^2 * (-1.0)
  c2 <- theta_ref_true * (-1.5) + theta_ref_true^2 * 0.8
  y <- theta_ref_true + c1*x1 + c2*x2 + rnorm(n, sd = 0.05)
  df <- data.frame(x1 = x1, x2 = x2, y = y)
  spec <- amm_spec(
    W = W_basis(type = "polynomial", degree = 2),
    x_vars = c("x1", "x2")
  )
  spec$W <- gdpar:::materialize_W_basis(spec$W, p = 1L)
  design <- gdpar:::build_amm_design(
    spec, df, formula_rhs = ~ x1 + x2
  )
  fam <- gdpar_family("gaussian")
  stan_data <- gdpar:::assemble_stan_data(
    design, fam, spec, df$y, theta_anchor = 0
  )
  res <- gdpar:::preflight_parametrization(
    gdpar_prior(), stan_data, spec,
    preflight_seed = 42L, verbose = FALSE
  )
  expect_true(res$cp_W)
  expect_match(res$meta$decision_reason_W, "filter_info_high")
})

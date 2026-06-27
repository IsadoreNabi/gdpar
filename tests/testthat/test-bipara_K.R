# Sub-phase 8.3.4 (2026-05-20 / 21): structural tests for K = 2
# bi-parametric families (Beta, Gamma), descriptor-based custom-family
# wiring (D-A3.B), per-slot rank check (D-B3), cross-slot extended
# Gram check (D-B2), and post-fit information contraction (D-B1).
# Companion gated smoke tests live in test-fit_K_smoke.R for Stan
# compilation paths.

# -----------------------------------------------------------------
# 1. Beta family construction (decision 1: (mu, log phi))
# -----------------------------------------------------------------

test_that("gdpar_family('beta') constructs (mu logit, phi log) param_specs", {
  fam <- gdpar_family("beta")
  expect_s3_class(fam, "gdpar_family")
  expect_equal(fam$name, "beta")
  expect_equal(fam$link, "logit")
  expect_equal(fam$stan_id, 5L)
  expect_true(fam$has_dispersion)
  expect_length(fam$param_specs, 2L)
  expect_equal(fam$param_specs[[1L]]$name, "mu")
  expect_equal(fam$param_specs[[1L]]$link, "logit")
  expect_equal(fam$param_specs[[1L]]$prior_canonical_kind, "logit_p")
  expect_equal(fam$param_specs[[1L]]$scope, "per_observation")
  expect_equal(fam$param_specs[[1L]]$family_role, "location")
  expect_equal(fam$param_specs[[2L]]$name, "phi")
  expect_equal(fam$param_specs[[2L]]$link, "log")
  expect_equal(fam$param_specs[[2L]]$prior_canonical_kind, "log_phi")
  expect_equal(fam$param_specs[[2L]]$scope, "population")
  expect_equal(fam$param_specs[[2L]]$family_role, "scale")
  expect_equal(fam$param_specs[[2L]]$support, "positive_real")
})

test_that("gdpar_family('beta') rejects non-canonical link", {
  expect_error(
    gdpar_family("beta", link = "identity"),
    class = "gdpar_input_error"
  )
  expect_error(
    gdpar_family("beta", link = "log"),
    class = "gdpar_input_error"
  )
})

# -----------------------------------------------------------------
# 2. Gamma family construction (decision 2: (mu, log shape))
# -----------------------------------------------------------------

test_that("gdpar_family('gamma') constructs (mu log, shape log) param_specs", {
  fam <- gdpar_family("gamma")
  expect_s3_class(fam, "gdpar_family")
  expect_equal(fam$name, "gamma")
  expect_equal(fam$link, "log")
  expect_equal(fam$stan_id, 6L)
  expect_true(fam$has_dispersion)
  expect_length(fam$param_specs, 2L)
  expect_equal(fam$param_specs[[1L]]$name, "mu")
  expect_equal(fam$param_specs[[1L]]$link, "log")
  expect_equal(fam$param_specs[[1L]]$prior_canonical_kind, "mu")
  expect_equal(fam$param_specs[[1L]]$scope, "per_observation")
  expect_equal(fam$param_specs[[1L]]$family_role, "location")
  expect_equal(fam$param_specs[[1L]]$support, "positive_real")
  expect_equal(fam$param_specs[[2L]]$name, "shape")
  expect_equal(fam$param_specs[[2L]]$link, "log")
  expect_equal(fam$param_specs[[2L]]$prior_canonical_kind, "log_shape")
  expect_equal(fam$param_specs[[2L]]$scope, "population")
  expect_equal(fam$param_specs[[2L]]$family_role, "shape")
})

test_that("gdpar_family('gamma') rejects non-canonical link", {
  expect_error(
    gdpar_family("gamma", link = "identity"),
    class = "gdpar_input_error"
  )
  expect_error(
    gdpar_family("gamma", link = "logit"),
    class = "gdpar_input_error"
  )
})

# -----------------------------------------------------------------
# 3. K = 1 guard for bi-parametric families
# -----------------------------------------------------------------

test_that(".gdpar_guard_K_below_family_min fires for Beta with K=1 (population phi)", {
  fam <- gdpar_family("beta")
  expect_equal(gdpar:::.gdpar_n_params_individual(fam), 1L)
  expect_error(
    gdpar:::.gdpar_guard_K_below_family_min(fam),
    class = "gdpar_unsupported_feature_error"
  )
})

test_that(".gdpar_guard_K_below_family_min fires for Gamma with K=1 (population shape)", {
  fam <- gdpar_family("gamma")
  expect_equal(gdpar:::.gdpar_n_params_individual(fam), 1L)
  expect_error(
    gdpar:::.gdpar_guard_K_below_family_min(fam),
    class = "gdpar_unsupported_feature_error"
  )
})

test_that(".gdpar_guard_K_below_family_min passes for non-bipara built-ins (gaussian K=1)", {
  fam <- gdpar_family("gaussian")
  expect_silent(gdpar:::.gdpar_guard_K_below_family_min(fam))
})

test_that(".gdpar_guard_K_below_family_min passes when Beta is elevated to K=2 via promotion", {
  fam <- gdpar_family("beta")
  fam_promoted <- gdpar:::.gdpar_promote_scope_per_observation(
    fam, c("mu", "phi")
  )
  expect_equal(gdpar:::.gdpar_n_params_individual(fam_promoted), 2L)
  expect_silent(gdpar:::.gdpar_guard_K_below_family_min(fam_promoted))
})

# -----------------------------------------------------------------
# 4. assemble_stan_data_K accepts stan_id 5/6/7
# -----------------------------------------------------------------

test_that(".assemble_stan_data_K accepts Beta K=2 (stan_id 5) and emits y_real", {
  set.seed(1L)
  n <- 30L
  df <- data.frame(y = pmin(pmax(rbeta(n, 2, 5), 0.05), 0.95),
                   x = rnorm(n), z = rnorm(n))
  fam <- gdpar_family("beta")
  fam_p <- gdpar:::.gdpar_promote_scope_per_observation(fam, c("mu", "phi"))
  fs <- gdpar_bf(y ~ a(x), phi ~ a(z))
  amm_list <- gdpar:::.gdpar_formula_set_to_amm_spec_list(fs, NULL)
  rhs <- y ~ x + z + 0
  design_K <- gdpar:::.build_amm_design_K(amm_list, df, formula_rhs = rhs)
  stan_data <- gdpar:::.assemble_stan_data_K(
    design_K = design_K, family = fam_p, amm_list_canonical = amm_list,
    y = df$y, theta_anchor_K = c(mu = 0, phi = 0)
  )
  expect_equal(stan_data$family_id_k, c(5L, 5L))
  expect_equal(length(stan_data$y_real), n)
  expect_true(all(stan_data$y_int == 0L))
})

test_that(".assemble_stan_data_K accepts Gamma K=2 (stan_id 6) and emits y_real", {
  set.seed(2L)
  n <- 30L
  df <- data.frame(y = rgamma(n, 2, 1), x = rnorm(n), z = rnorm(n))
  fam <- gdpar_family("gamma")
  fam_p <- gdpar:::.gdpar_promote_scope_per_observation(fam, c("mu", "shape"))
  fs <- gdpar_bf(y ~ a(x), shape ~ a(z))
  amm_list <- gdpar:::.gdpar_formula_set_to_amm_spec_list(fs, NULL)
  rhs <- y ~ x + z + 0
  design_K <- gdpar:::.build_amm_design_K(amm_list, df, formula_rhs = rhs)
  stan_data <- gdpar:::.assemble_stan_data_K(
    design_K = design_K, family = fam_p, amm_list_canonical = amm_list,
    y = df$y, theta_anchor_K = c(mu = 0, shape = 0)
  )
  expect_equal(stan_data$family_id_k, c(6L, 6L))
  expect_equal(length(stan_data$y_real), n)
})

test_that(".assemble_stan_data_K rejects Beta y outside (0,1)", {
  set.seed(3L)
  n <- 30L
  df <- data.frame(y = c(0, runif(n - 1)), x = rnorm(n), z = rnorm(n))
  fam <- gdpar_family("beta")
  fam_p <- gdpar:::.gdpar_promote_scope_per_observation(fam, c("mu", "phi"))
  fs <- gdpar_bf(y ~ a(x), phi ~ a(z))
  amm_list <- gdpar:::.gdpar_formula_set_to_amm_spec_list(fs, NULL)
  rhs <- y ~ x + z + 0
  design_K <- gdpar:::.build_amm_design_K(amm_list, df, formula_rhs = rhs)
  expect_error(
    gdpar:::.assemble_stan_data_K(
      design_K = design_K, family = fam_p, amm_list_canonical = amm_list,
      y = df$y, theta_anchor_K = c(mu = 0, phi = 0)
    ),
    class = "gdpar_input_error"
  )
})

test_that(".assemble_stan_data_K rejects Gamma y <= 0", {
  set.seed(4L)
  n <- 30L
  df <- data.frame(y = c(-1, rgamma(n - 1, 2, 1)),
                   x = rnorm(n), z = rnorm(n))
  fam <- gdpar_family("gamma")
  fam_p <- gdpar:::.gdpar_promote_scope_per_observation(fam, c("mu", "shape"))
  fs <- gdpar_bf(y ~ a(x), shape ~ a(z))
  amm_list <- gdpar:::.gdpar_formula_set_to_amm_spec_list(fs, NULL)
  rhs <- y ~ x + z + 0
  design_K <- gdpar:::.build_amm_design_K(amm_list, df, formula_rhs = rhs)
  expect_error(
    gdpar:::.assemble_stan_data_K(
      design_K = design_K, family = fam_p, amm_list_canonical = amm_list,
      y = df$y, theta_anchor_K = c(mu = 0, shape = 0)
    ),
    class = "gdpar_input_error"
  )
})

# -----------------------------------------------------------------
# 5. D-A3.B custom-family wiring (descriptor-based, whitelist)
# -----------------------------------------------------------------

test_that("gdpar_family_custom_K builds a lognormal_loc_scale family", {
  fam <- suppressMessages(gdpar_family_custom_K(
    name          = "my_ln",
    stan_lpdf_id  = "lognormal_loc_scale",
    did_holds     = TRUE,
    did_reference = "test"
  ))
  expect_s3_class(fam, "gdpar_family")
  expect_equal(fam$name, "my_ln")
  expect_equal(fam$stan_id, 7L)
  expect_true(fam$has_dispersion)
  expect_true(isTRUE(fam$is_custom))
  expect_equal(fam$stan_lpdf_id, "lognormal_loc_scale")
  expect_length(fam$param_specs, 2L)
  expect_equal(fam$param_specs[[1L]]$name, "mu")
  expect_equal(fam$param_specs[[1L]]$link, "identity")
  expect_equal(fam$param_specs[[2L]]$name, "sigma")
  expect_equal(fam$param_specs[[2L]]$link, "log")
})

test_that("gdpar_family_custom_K rejects builtin name collision", {
  expect_error(
    suppressMessages(gdpar_family_custom_K(
      name          = "gaussian",
      stan_lpdf_id  = "lognormal_loc_scale",
      did_holds     = TRUE,
      did_reference = "test"
    )),
    class = "gdpar_input_error"
  )
})

test_that("gdpar_family_custom_K rejects unknown stan_lpdf_id", {
  expect_error(
    gdpar_family_custom_K(
      name          = "foo",
      stan_lpdf_id  = "nope_not_in_registry",
      did_holds     = TRUE,
      did_reference = "test"
    ),
    class = "gdpar_input_error"
  )
})

test_that(".gdpar_K_custom_patterns registry shape is stable", {
  reg <- gdpar:::.gdpar_K_custom_patterns()
  expect_true("lognormal_loc_scale" %in% names(reg))
  ln <- reg$lognormal_loc_scale
  expect_equal(ln$stan_id, 7L)
  expect_equal(ln$y_type, "real")
  expect_length(ln$slot_specs, 2L)
  expect_equal(ln$slot_specs[[1L]]$name, "mu")
  expect_equal(ln$slot_specs[[2L]]$name, "sigma")
})

test_that("custom_K family routes through K = 2 assembler with stan_id 7", {
  set.seed(5L)
  n <- 30L
  df <- data.frame(y = rlnorm(n, 0.2, 0.5), x = rnorm(n), z = rnorm(n))
  fam <- suppressMessages(gdpar_family_custom_K(
    name = "my_ln", stan_lpdf_id = "lognormal_loc_scale",
    did_holds = TRUE, did_reference = "test"
  ))
  fam_p <- gdpar:::.gdpar_promote_scope_per_observation(fam, c("mu", "sigma"))
  fs <- gdpar_bf(y ~ a(x), sigma ~ a(z))
  amm_list <- gdpar:::.gdpar_formula_set_to_amm_spec_list(fs, NULL)
  rhs <- y ~ x + z + 0
  design_K <- gdpar:::.build_amm_design_K(amm_list, df, formula_rhs = rhs)
  stan_data <- gdpar:::.assemble_stan_data_K(
    design_K = design_K, family = fam_p, amm_list_canonical = amm_list,
    y = df$y, theta_anchor_K = c(mu = 0, sigma = 0)
  )
  expect_equal(stan_data$family_id_k, c(7L, 7L))
  expect_equal(length(stan_data$y_real), n)
})

# -----------------------------------------------------------------
# 6. D-B3 per-slot Z_a^{(k)} rank check
# -----------------------------------------------------------------

mock_design_K <- function(K = 2L, n = 50L, J = 2L,
                          slot_names = c("mu", "sigma"),
                          seed = 11L) {
  set.seed(seed)
  Z_list <- replicate(K, matrix(rnorm(n * J), n, J), simplify = FALSE)
  names_list <- replicate(K, paste0("x", seq_len(J)), simplify = FALSE)
  list(
    Z_a_k_list       = Z_list,
    Z_b_k_list       = replicate(K, matrix(0, n, 0L), simplify = FALSE),
    Z_a_k_names_list = names_list,
    slot_names       = slot_names,
    K                = as.integer(K),
    X                = matrix(0, n, 0L)
  )
}

test_that(".check_Z_a_K_per_slot passes for well-conditioned per-slot design", {
  d <- mock_design_K()
  rep <- gdpar:::.check_Z_a_K_per_slot(d, rigor = "full", tol = 1e-8)
  expect_true(rep$passed)
  expect_equal(rep$rigor, "full")
  expect_length(rep$per_slot, 2L)
  expect_true(all(vapply(rep$per_slot,
                         function(s) isTRUE(s$passed), logical(1L))))
})

test_that(".check_Z_a_K_per_slot detects per-slot collinearity (slot 1 only)", {
  d <- mock_design_K()
  d$Z_a_k_list[[1L]] <- cbind(d$Z_a_k_list[[1L]][, 1L],
                              2 * d$Z_a_k_list[[1L]][, 1L])
  rep <- gdpar:::.check_Z_a_K_per_slot(d, rigor = "full", tol = 1e-8)
  expect_false(rep$passed)
  expect_false(rep$per_slot[[1L]]$passed)
  expect_true(rep$per_slot[[2L]]$passed)
  expect_lt(rep$per_slot[[1L]]$condition_number^-1, 1e-8)
})

test_that(".check_Z_a_K_per_slot detects zero-norm column", {
  d <- mock_design_K()
  d$Z_a_k_list[[2L]][, 2L] <- 0
  rep <- gdpar:::.check_Z_a_K_per_slot(d, rigor = "full", tol = 1e-8)
  expect_false(rep$passed)
  expect_false(rep$per_slot[[2L]]$passed)
  expect_identical(rep$per_slot[[2L]]$collinear_columns, "x2")
})

test_that(".check_Z_a_K_per_slot with rigor='fast' skips eigendecomposition", {
  d <- mock_design_K()
  rep <- gdpar:::.check_Z_a_K_per_slot(d, rigor = "fast", tol = 1e-8)
  expect_true(rep$passed)
  for (s in rep$per_slot) {
    expect_equal(s$rigor, "fast")
    expect_true(is.na(s$condition_number))
  }
})

test_that(".check_Z_a_K_per_slot passes trivially on empty Z_a", {
  d <- mock_design_K()
  d$Z_a_k_list[[1L]] <- matrix(0, nrow = 50L, ncol = 0L)
  d$Z_a_k_names_list[[1L]] <- character(0L)
  rep <- gdpar:::.check_Z_a_K_per_slot(d, rigor = "full", tol = 1e-8)
  expect_true(rep$passed)
  expect_true(rep$per_slot[[1L]]$passed)
  expect_match(rep$per_slot[[1L]]$message, "empty")
})

# -----------------------------------------------------------------
# 7. D-B2 cross-slot extended Gram check
# -----------------------------------------------------------------

test_that(".check_C4_bis_K_cross_slot passes for non-overlapping slots", {
  d <- mock_design_K()
  rep <- gdpar:::.check_C4_bis_K_cross_slot(d, rigor = "full", tol = 1e-8)
  expect_true(rep$passed)
  expect_equal(rep$total_columns, 4L)
})

test_that(".check_C4_bis_K_cross_slot detects cross-slot collinearity", {
  d <- mock_design_K()
  d$Z_a_k_list[[2L]] <- d$Z_a_k_list[[1L]]
  d$Z_a_k_names_list[[2L]] <- d$Z_a_k_names_list[[1L]]
  rep <- gdpar:::.check_C4_bis_K_cross_slot(d, rigor = "full", tol = 1e-8)
  expect_false(rep$passed)
  expect_false(is.null(rep$collinear_directions))
})

test_that(".check_C4_bis_K_cross_slot fast mode skips eigendecomp", {
  d <- mock_design_K()
  rep <- gdpar:::.check_C4_bis_K_cross_slot(d, rigor = "fast", tol = 1e-8)
  expect_true(rep$passed)
  expect_true(is.na(rep$condition_number))
})

test_that(".check_C4_bis_K_cross_slot handles all-empty Z_a (vacuous pass)", {
  d <- mock_design_K()
  for (k in seq_along(d$Z_a_k_list)) {
    d$Z_a_k_list[[k]] <- matrix(0, nrow = 50L, ncol = 0L)
    d$Z_a_k_names_list[[k]] <- character(0L)
  }
  rep <- gdpar:::.check_C4_bis_K_cross_slot(d, rigor = "full", tol = 1e-8)
  expect_true(rep$passed)
  expect_equal(rep$total_columns, 0L)
})

# -----------------------------------------------------------------
# 8. .check_identifiability_K orchestrator
# -----------------------------------------------------------------

test_that(".check_identifiability_K combines D-B2 and D-B3 layers", {
  d <- mock_design_K()
  rep <- gdpar:::.check_identifiability_K(d, rigor = "full")
  expect_true(rep$passed)
  expect_equal(rep$K, 2L)
  expect_named(rep, c("passed", "rigor", "tol", "K",
                      "per_slot_rank", "cross_slot_gram"))
})

test_that(".check_identifiability_K fails when D-B3 fails", {
  d <- mock_design_K()
  d$Z_a_k_list[[1L]][, 2L] <- 0
  rep <- gdpar:::.check_identifiability_K(d, rigor = "full")
  expect_false(rep$passed)
  expect_false(rep$per_slot_rank$passed)
})

test_that(".check_identifiability_K fails when D-B2 fails (clean per-slot)", {
  d <- mock_design_K(seed = 13L)
  d$Z_a_k_list[[2L]] <- d$Z_a_k_list[[1L]]
  d$Z_a_k_names_list[[2L]] <- d$Z_a_k_names_list[[1L]]
  rep <- gdpar:::.check_identifiability_K(d, rigor = "full")
  expect_false(rep$passed)
  expect_true(rep$per_slot_rank$passed)
  expect_false(rep$cross_slot_gram$passed)
})

# -----------------------------------------------------------------
# 9. D-B1 prior canonical variance lookup
# -----------------------------------------------------------------

test_that(".gdpar_canonical_prior_variance covers all registered kinds", {
  expect_equal(gdpar:::.gdpar_canonical_prior_variance("mu"), 6.25)
  expect_equal(gdpar:::.gdpar_canonical_prior_variance("log_sigma"), 1.0)
  expect_equal(gdpar:::.gdpar_canonical_prior_variance("log_phi"), 1.0)
  expect_equal(gdpar:::.gdpar_canonical_prior_variance("logit_p"), 6.25)
  expect_equal(gdpar:::.gdpar_canonical_prior_variance("log_shape"), 1.0)
  expect_equal(gdpar:::.gdpar_canonical_prior_variance("log_nu"), 1.0)
  expect_equal(gdpar:::.gdpar_canonical_prior_variance("logit_pi"), 6.25)
  expect_equal(gdpar:::.gdpar_canonical_prior_variance("power_p"),
               (1.99 - 1.01)^2 / 12)
  expect_true(is.na(gdpar:::.gdpar_canonical_prior_variance("unknown_kind")))
})

# -----------------------------------------------------------------
# 10. D-B1 .compute_info_ratio_K classification
# -----------------------------------------------------------------

mock_fit_with_draws <- function(draws_list) {
  draws_mat <- do.call(cbind, draws_list)
  storage.mode(draws_mat) <- "double"
  list(
    draws = function(format = "draws_matrix") draws_mat
  )
}

test_that(".compute_info_ratio_K classifies pass/warn/info_error by contraction", {
  set.seed(91L)
  fam <- gdpar_family("gaussian")
  fam_p <- gdpar:::.gdpar_promote_scope_per_observation(fam, c("mu", "sigma"))
  # Slot 1 (mu, kind 'mu', prior var 6.25): post var = 0.5 -> contraction 0.92 -> pass
  # Slot 2 (sigma, kind 'log_sigma', prior var 1): post var = 0.95 -> contraction 0.05 -> info_error
  draws <- list(
    `theta_ref_k[1,1]` = rnorm(800, 0, sqrt(0.5)),
    `theta_ref_k[1,2]` = rnorm(800, 0, sqrt(0.95))
  )
  fit_mock <- mock_fit_with_draws(draws)
  rep <- gdpar:::.compute_info_ratio_K(
    fit_mock, fam_p, slot_names = c("mu", "sigma"),
    use_groups = 0L, prior = NULL
  )
  expect_equal(rep$per_slot[[1L]]$status, "pass")
  expect_equal(rep$per_slot[[2L]]$status, "information_error")
  expect_true(rep$any_info_error)
  expect_false(rep$passed)
})

test_that(".compute_info_ratio_K returns 'warn' for contraction in [0.1, 0.5)", {
  set.seed(92L)
  fam <- gdpar_family("gaussian")
  fam_p <- gdpar:::.gdpar_promote_scope_per_observation(fam, c("mu", "sigma"))
  # Slot 1: post var = 4.0 / prior 6.25 -> contraction = 0.36 -> warn
  draws <- list(
    `theta_ref_k[1,1]` = rnorm(800, 0, sqrt(4.0)),
    `theta_ref_k[1,2]` = rnorm(800, 0, sqrt(0.1))
  )
  fit_mock <- mock_fit_with_draws(draws)
  rep <- gdpar:::.compute_info_ratio_K(
    fit_mock, fam_p, slot_names = c("mu", "sigma"),
    use_groups = 0L, prior = NULL
  )
  expect_equal(rep$per_slot[[1L]]$status, "warn")
  expect_equal(rep$per_slot[[2L]]$status, "pass")
  expect_true(rep$any_warn)
  expect_false(rep$any_info_error)
})

test_that(".compute_info_ratio_K respects use_groups = 1 (mu_theta_ref_k)", {
  set.seed(93L)
  fam <- gdpar_family("gaussian")
  fam_p <- gdpar:::.gdpar_promote_scope_per_observation(fam, c("mu", "sigma"))
  draws <- list(
    `mu_theta_ref_k[1,1]` = rnorm(800, 0, sqrt(0.5)),
    `mu_theta_ref_k[1,2]` = rnorm(800, 0, sqrt(0.5))
  )
  fit_mock <- mock_fit_with_draws(draws)
  rep <- gdpar:::.compute_info_ratio_K(
    fit_mock, fam_p, slot_names = c("mu", "sigma"),
    use_groups = 1L, prior = NULL
  )
  expect_equal(rep$per_slot[[1L]]$status, "pass")
  expect_equal(rep$per_slot[[2L]]$status, "pass")
})

test_that(".compute_info_ratio_K skips slots when draws are missing", {
  fam <- gdpar_family("gaussian")
  fam_p <- gdpar:::.gdpar_promote_scope_per_observation(fam, c("mu", "sigma"))
  fit_mock <- mock_fit_with_draws(list(
    `theta_ref_k[1,1]` = rnorm(100)
    # slot 2 draws intentionally missing
  ))
  rep <- gdpar:::.compute_info_ratio_K(
    fit_mock, fam_p, slot_names = c("mu", "sigma"),
    use_groups = 0L, prior = NULL
  )
  expect_equal(rep$per_slot[[2L]]$status, "skipped")
})

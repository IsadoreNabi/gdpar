# Sub-phase 8.3.5b (2026-05-21): structural tests for K = 3 Tweedie
# family (stan_id 9, mu log + phi log + p identity bounded on
# (1.01, 1.99)) plus the THETA_REF_PRIOR_BLOCK builder that enforces
# the per-slot prior dispatch per decision E7 of D5. The K = 3 smoke
# fit lives in test-fit_K_smoke.R gated by
# GDPAR_RUN_STAN_SMOKE_TWEEDIE.

# -----------------------------------------------------------------
# 1. Tweedie family construction
# -----------------------------------------------------------------

test_that("gdpar_family('tweedie') constructs (mu log, phi log, p identity) param_specs", {
  fam <- gdpar_family("tweedie")
  expect_s3_class(fam, "gdpar_family")
  expect_equal(fam$name, "tweedie")
  expect_equal(fam$link, "log")
  expect_equal(fam$stan_id, 9L)
  expect_true(fam$has_dispersion)
  expect_length(fam$param_specs, 3L)
  expect_equal(fam$param_specs[[1L]]$name, "mu")
  expect_equal(fam$param_specs[[1L]]$link, "log")
  expect_equal(fam$param_specs[[1L]]$prior_canonical_kind, "mu")
  expect_equal(fam$param_specs[[1L]]$scope, "per_observation")
  expect_equal(fam$param_specs[[1L]]$family_role, "location")
  expect_equal(fam$param_specs[[1L]]$support, "positive_real")
  expect_equal(fam$param_specs[[2L]]$name, "phi")
  expect_equal(fam$param_specs[[2L]]$link, "log")
  expect_equal(fam$param_specs[[2L]]$prior_canonical_kind, "log_phi")
  expect_equal(fam$param_specs[[2L]]$scope, "population")
  expect_equal(fam$param_specs[[2L]]$family_role, "scale")
  expect_equal(fam$param_specs[[2L]]$support, "positive_real")
  expect_equal(fam$param_specs[[3L]]$name, "p")
  expect_equal(fam$param_specs[[3L]]$link, "identity")
  expect_equal(fam$param_specs[[3L]]$prior_canonical_kind, "power_p")
  expect_equal(fam$param_specs[[3L]]$scope, "population")
  expect_equal(fam$param_specs[[3L]]$family_role, "power")
  expect_equal(fam$param_specs[[3L]]$support, "bounded_open")
})

test_that("gdpar_family('tweedie') rejects non-canonical links", {
  expect_error(
    gdpar_family("tweedie", link = "identity"),
    class = "gdpar_input_error"
  )
  expect_error(
    gdpar_family("tweedie", link = "logit"),
    class = "gdpar_input_error"
  )
})

# -----------------------------------------------------------------
# 2. K < 3 guard for Tweedie (stan_id 9 -> min_K 3)
# -----------------------------------------------------------------

test_that(".gdpar_guard_K_below_family_min fires for Tweedie with K=1", {
  fam <- gdpar_family("tweedie")
  expect_equal(gdpar:::.gdpar_n_params_individual(fam), 1L)
  expect_error(
    gdpar:::.gdpar_guard_K_below_family_min(fam),
    class = "gdpar_unsupported_feature_error"
  )
})

test_that(".gdpar_guard_K_below_family_min fires for Tweedie with K=2 (mu + phi promoted)", {
  fam <- gdpar_family("tweedie")
  fam_p <- gdpar:::.gdpar_promote_scope_per_observation(fam,
                                                       c("mu", "phi"))
  expect_equal(gdpar:::.gdpar_n_params_individual(fam_p), 2L)
  expect_error(
    gdpar:::.gdpar_guard_K_below_family_min(fam_p),
    class = "gdpar_unsupported_feature_error"
  )
})

test_that(".gdpar_guard_K_below_family_min passes when Tweedie is elevated to K=3", {
  fam <- gdpar_family("tweedie")
  fam_p <- gdpar:::.gdpar_promote_scope_per_observation(fam,
                                                       c("mu", "phi", "p"))
  expect_equal(gdpar:::.gdpar_n_params_individual(fam_p), 3L)
  expect_silent(gdpar:::.gdpar_guard_K_below_family_min(fam_p))
})

# -----------------------------------------------------------------
# 3. assemble_stan_data_K accepts Tweedie K=3 (stan_id 9)
# -----------------------------------------------------------------

test_that(".assemble_stan_data_K accepts Tweedie K=3 (stan_id 9) and emits y_real", {
  set.seed(91L)
  n <- 30L
  df <- data.frame(y = c(rep(0, 5L), exp(rnorm(n - 5L))),
                   x = rnorm(n), z = rnorm(n), w = rnorm(n))
  fam <- gdpar_family("tweedie")
  fam_p <- gdpar:::.gdpar_promote_scope_per_observation(fam,
                                                       c("mu", "phi", "p"))
  fs <- gdpar_bf(y ~ a(x), phi ~ a(z), p ~ a(w))
  amm_list <- gdpar:::.gdpar_formula_set_to_amm_spec_list(fs, NULL)
  rhs <- y ~ x + z + w + 0
  design_K <- gdpar:::.build_amm_design_K(amm_list, df, formula_rhs = rhs)
  stan_data <- gdpar:::.assemble_stan_data_K(
    design_K = design_K, family = fam_p, amm_list_canonical = amm_list,
    y = df$y, theta_anchor_K = c(mu = 0, phi = 0, p = 1.5)
  )
  expect_equal(stan_data$family_id_k, c(9L, 9L, 9L))
  expect_equal(stan_data$K, 3L)
  expect_equal(length(stan_data$y_real), n)
  expect_true(all(stan_data$y_int == 0L))
})

test_that(".assemble_stan_data_K accepts Tweedie with y = 0 (point mass)", {
  set.seed(92L)
  n <- 30L
  df <- data.frame(y = c(rep(0, 10L), exp(rnorm(n - 10L))),
                   x = rnorm(n), z = rnorm(n), w = rnorm(n))
  fam <- gdpar_family("tweedie")
  fam_p <- gdpar:::.gdpar_promote_scope_per_observation(fam,
                                                       c("mu", "phi", "p"))
  fs <- gdpar_bf(y ~ a(x), phi ~ a(z), p ~ a(w))
  amm_list <- gdpar:::.gdpar_formula_set_to_amm_spec_list(fs, NULL)
  rhs <- y ~ x + z + w + 0
  design_K <- gdpar:::.build_amm_design_K(amm_list, df, formula_rhs = rhs)
  expect_silent(
    gdpar:::.assemble_stan_data_K(
      design_K = design_K, family = fam_p, amm_list_canonical = amm_list,
      y = df$y, theta_anchor_K = c(mu = 0, phi = 0, p = 1.5)
    )
  )
})

test_that(".assemble_stan_data_K rejects Tweedie with y < 0", {
  set.seed(93L)
  n <- 30L
  df <- data.frame(y = c(-1, exp(rnorm(n - 1L))),
                   x = rnorm(n), z = rnorm(n), w = rnorm(n))
  fam <- gdpar_family("tweedie")
  fam_p <- gdpar:::.gdpar_promote_scope_per_observation(fam,
                                                       c("mu", "phi", "p"))
  fs <- gdpar_bf(y ~ a(x), phi ~ a(z), p ~ a(w))
  amm_list <- gdpar:::.gdpar_formula_set_to_amm_spec_list(fs, NULL)
  rhs <- y ~ x + z + w + 0
  design_K <- gdpar:::.build_amm_design_K(amm_list, df, formula_rhs = rhs)
  expect_error(
    gdpar:::.assemble_stan_data_K(
      design_K = design_K, family = fam_p, amm_list_canonical = amm_list,
      y = df$y, theta_anchor_K = c(mu = 0, phi = 0, p = 1.5)
    ),
    class = "gdpar_input_error"
  )
})

# -----------------------------------------------------------------
# 4. Custom family blacklist: tweedie collision
# -----------------------------------------------------------------

test_that("gdpar_family_custom rejects 'tweedie' as builtin name collision", {
  expect_error(
    gdpar_family_custom(
      name               = "tweedie",
      link               = "log",
      did_holds          = TRUE,
      did_condition      = NA_character_,
      stan_loglik_block  = "target += normal_lpdf(y_real[i] | eta[i], 1);",
      stan_log_lik_block = "log_lik[i] = normal_lpdf(y_real[i] | eta[i], 1);",
      stan_y_pred_block  = "y_pred[i] = normal_rng(eta[i], 1);",
      y_type             = "real",
      did_reference      = "test"
    ),
    class = "gdpar_input_error"
  )
})

test_that("gdpar_family_custom_K rejects 'tweedie' as builtin name collision", {
  expect_error(
    gdpar_family_custom_K(
      name         = "tweedie",
      stan_lpdf_id = "lognormal_loc_scale",
      did_holds    = TRUE
    ),
    class = "gdpar_input_error"
  )
})

# -----------------------------------------------------------------
# 5. .gdpar_canonical_prior_variance covers power_p (already registered)
# -----------------------------------------------------------------

test_that(".gdpar_canonical_prior_variance(power_p) returns (1.99-1.01)^2/12", {
  expect_equal(gdpar:::.gdpar_canonical_prior_variance("power_p"),
               (1.99 - 1.01)^2 / 12)
})

test_that(".gdpar_canonical_prior_for_kind(power_p) returns uniform(1.01, 1.99)", {
  expect_equal(gdpar:::.gdpar_canonical_prior_for_kind("power_p"),
               "uniform(1.01, 1.99)")
})

# -----------------------------------------------------------------
# 6. THETA_REF_PRIOR_BLOCK builder (E7 of D5)
# -----------------------------------------------------------------

test_that(".gdpar_build_theta_ref_prior_block(NULL) returns the default vectorized block", {
  block <- gdpar:::.gdpar_build_theta_ref_prior_block(NULL)
  expect_match(block, "theta_ref_k\\[1\\] ~ \\{\\{PRIOR_THETA_REF\\}\\};",
               fixed = FALSE)
  expect_false(grepl("uniform(1.01, 1.99)", block, fixed = TRUE))
})

test_that(".gdpar_build_theta_ref_prior_block(Gaussian) returns the default block", {
  fam <- gdpar_family("gaussian")
  block <- gdpar:::.gdpar_build_theta_ref_prior_block(fam)
  expect_match(block, "theta_ref_k\\[1\\] ~ \\{\\{PRIOR_THETA_REF\\}\\};",
               fixed = FALSE)
  expect_false(grepl("uniform(1.01, 1.99)", block, fixed = TRUE))
})

test_that(".gdpar_build_theta_ref_prior_block(Student-t) returns the default block", {
  fam <- gdpar_family("student_t")
  block <- gdpar:::.gdpar_build_theta_ref_prior_block(fam)
  expect_match(block, "theta_ref_k\\[1\\] ~ \\{\\{PRIOR_THETA_REF\\}\\};",
               fixed = FALSE)
  expect_false(grepl("uniform(1.01, 1.99)", block, fixed = TRUE))
})

test_that(".gdpar_build_theta_ref_prior_block(Tweedie) returns the slice + uniform block", {
  fam <- gdpar_family("tweedie")
  block <- gdpar:::.gdpar_build_theta_ref_prior_block(fam)
  expect_true(grepl("theta_ref_k[1][1:2] ~ {{PRIOR_THETA_REF}}", block,
                    fixed = TRUE))
  expect_true(grepl("theta_ref_k[1][3] ~ uniform(1.01, 1.99)", block,
                    fixed = TRUE))
  expect_true(grepl("mu_theta_ref_k[1][1:2] ~ {{PRIOR_THETA_REF}}", block,
                    fixed = TRUE))
  expect_true(grepl("mu_theta_ref_k[1][3] ~ uniform(1.01, 1.99)", block,
                    fixed = TRUE))
})

test_that("generate_stan_code_K leaves no unsubstituted placeholders for Tweedie", {
  fam <- gdpar_family("tweedie")
  src <- gdpar:::generate_stan_code_K(gdpar_prior(), family = fam)
  expect_false(grepl("{{", src, fixed = TRUE))
  expect_true(grepl("theta_ref_k[1][1:2] ~", src, fixed = TRUE))
  expect_true(grepl("theta_ref_k[1][3] ~ uniform(1.01, 1.99)", src,
                    fixed = TRUE))
  expect_true(grepl("tweedie_lpdf", src, fixed = TRUE))
  expect_true(grepl("tweedie_rng", src, fixed = TRUE))
})

test_that("generate_stan_code_K with default family preserves vectorized block", {
  src <- gdpar:::generate_stan_code_K(gdpar_prior())
  expect_false(grepl("{{", src, fixed = TRUE))
  expect_true(grepl("theta_ref_k[1] ~ normal", src, fixed = TRUE))
  expect_false(grepl("uniform(1.01, 1.99)", src, fixed = TRUE))
})

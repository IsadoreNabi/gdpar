# Sub-phase 8.3.6 (2026-05-21): structural tests for the four count
# mixture families wired in the K-individual distributional regression
# path: ZIP K = 2 (stan_id 10), ZINB K = 3 (stan_id 11), Hurdle-Poisson
# K = 2 (stan_id 12), Hurdle-NB K = 3 (stan_id 13). All four share the
# canonical (mu log, [phi log,] pi logit) parametrization and the
# logit_pi prior_canonical_kind already registered in prior_spec.R.
# Smoke fits live in test-fit_K_smoke.R gated by family-specific env
# vars.

# -----------------------------------------------------------------
# 1. Family construction
# -----------------------------------------------------------------

test_that("gdpar_family('zip') constructs (mu log, pi logit) param_specs", {
  fam <- gdpar_family("zip")
  expect_s3_class(fam, "gdpar_family")
  expect_equal(fam$name, "zip")
  expect_equal(fam$link, "log")
  expect_equal(fam$stan_id, 10L)
  expect_true(fam$has_dispersion)
  expect_length(fam$param_specs, 2L)
  expect_equal(fam$param_specs[[1L]]$name, "mu")
  expect_equal(fam$param_specs[[1L]]$link, "log")
  expect_equal(fam$param_specs[[1L]]$prior_canonical_kind, "mu")
  expect_equal(fam$param_specs[[1L]]$scope, "per_observation")
  expect_equal(fam$param_specs[[1L]]$family_role, "location")
  expect_equal(fam$param_specs[[1L]]$support, "positive_real")
  expect_equal(fam$param_specs[[2L]]$name, "pi")
  expect_equal(fam$param_specs[[2L]]$link, "logit")
  expect_equal(fam$param_specs[[2L]]$prior_canonical_kind, "logit_pi")
  expect_equal(fam$param_specs[[2L]]$scope, "population")
  expect_equal(fam$param_specs[[2L]]$family_role, "mixture_pi")
  expect_equal(fam$param_specs[[2L]]$support, "unit_interval")
})

test_that("gdpar_family('zinb') constructs (mu log, phi log, pi logit) param_specs", {
  fam <- gdpar_family("zinb")
  expect_s3_class(fam, "gdpar_family")
  expect_equal(fam$name, "zinb")
  expect_equal(fam$stan_id, 11L)
  expect_length(fam$param_specs, 3L)
  expect_equal(fam$param_specs[[1L]]$name, "mu")
  expect_equal(fam$param_specs[[1L]]$family_role, "location")
  expect_equal(fam$param_specs[[2L]]$name, "phi")
  expect_equal(fam$param_specs[[2L]]$link, "log")
  expect_equal(fam$param_specs[[2L]]$prior_canonical_kind, "log_phi")
  expect_equal(fam$param_specs[[2L]]$family_role, "scale")
  expect_equal(fam$param_specs[[3L]]$name, "pi")
  expect_equal(fam$param_specs[[3L]]$link, "logit")
  expect_equal(fam$param_specs[[3L]]$prior_canonical_kind, "logit_pi")
  expect_equal(fam$param_specs[[3L]]$family_role, "mixture_pi")
})

test_that("gdpar_family('hurdle_poisson') constructs (mu log, pi logit) param_specs", {
  fam <- gdpar_family("hurdle_poisson")
  expect_equal(fam$stan_id, 12L)
  expect_length(fam$param_specs, 2L)
  expect_equal(fam$param_specs[[1L]]$name, "mu")
  expect_equal(fam$param_specs[[1L]]$link, "log")
  expect_equal(fam$param_specs[[2L]]$name, "pi")
  expect_equal(fam$param_specs[[2L]]$link, "logit")
  expect_equal(fam$param_specs[[2L]]$family_role, "mixture_pi")
})

test_that("gdpar_family('hurdle_neg_binomial_2') constructs (mu log, phi log, pi logit) param_specs", {
  fam <- gdpar_family("hurdle_neg_binomial_2")
  expect_equal(fam$stan_id, 13L)
  expect_length(fam$param_specs, 3L)
  expect_equal(fam$param_specs[[1L]]$name, "mu")
  expect_equal(fam$param_specs[[2L]]$name, "phi")
  expect_equal(fam$param_specs[[2L]]$link, "log")
  expect_equal(fam$param_specs[[2L]]$family_role, "scale")
  expect_equal(fam$param_specs[[3L]]$name, "pi")
  expect_equal(fam$param_specs[[3L]]$link, "logit")
})

test_that("mixture families reject non-canonical links", {
  for (nm in c("zip", "zinb", "hurdle_poisson",
               "hurdle_neg_binomial_2")) {
    expect_error(
      gdpar_family(nm, link = "identity"),
      class = "gdpar_input_error",
      info = sprintf("family '%s' should reject link='identity'", nm)
    )
    expect_error(
      gdpar_family(nm, link = "logit"),
      class = "gdpar_input_error",
      info = sprintf("family '%s' should reject link='logit' on mu", nm)
    )
  }
})

# -----------------------------------------------------------------
# 2. K < min_K guard fires per family
# -----------------------------------------------------------------

test_that(".gdpar_guard_K_below_family_min fires for ZIP K=1", {
  fam <- gdpar_family("zip")
  expect_equal(gdpar:::.gdpar_n_params_individual(fam), 1L)
  expect_error(
    gdpar:::.gdpar_guard_K_below_family_min(fam),
    class = "gdpar_unsupported_feature_error"
  )
})

test_that(".gdpar_guard_K_below_family_min passes when ZIP is elevated to K=2", {
  fam <- gdpar_family("zip")
  fam_p <- gdpar:::.gdpar_promote_scope_per_observation(fam,
                                                       c("mu", "pi"))
  expect_equal(gdpar:::.gdpar_n_params_individual(fam_p), 2L)
  expect_silent(gdpar:::.gdpar_guard_K_below_family_min(fam_p))
})

test_that(".gdpar_guard_K_below_family_min fires for ZINB K=1 and K=2", {
  fam <- gdpar_family("zinb")
  expect_error(
    gdpar:::.gdpar_guard_K_below_family_min(fam),
    class = "gdpar_unsupported_feature_error"
  )
  fam_p2 <- gdpar:::.gdpar_promote_scope_per_observation(fam,
                                                        c("mu", "phi"))
  expect_error(
    gdpar:::.gdpar_guard_K_below_family_min(fam_p2),
    class = "gdpar_unsupported_feature_error"
  )
})

test_that(".gdpar_guard_K_below_family_min passes when ZINB is elevated to K=3", {
  fam <- gdpar_family("zinb")
  fam_p <- gdpar:::.gdpar_promote_scope_per_observation(fam,
                                                       c("mu", "phi", "pi"))
  expect_equal(gdpar:::.gdpar_n_params_individual(fam_p), 3L)
  expect_silent(gdpar:::.gdpar_guard_K_below_family_min(fam_p))
})

test_that(".gdpar_guard_K_below_family_min fires for hurdle_poisson K=1", {
  fam <- gdpar_family("hurdle_poisson")
  expect_error(
    gdpar:::.gdpar_guard_K_below_family_min(fam),
    class = "gdpar_unsupported_feature_error"
  )
  fam_p <- gdpar:::.gdpar_promote_scope_per_observation(fam,
                                                       c("mu", "pi"))
  expect_silent(gdpar:::.gdpar_guard_K_below_family_min(fam_p))
})

test_that(".gdpar_guard_K_below_family_min fires for hurdle_neg_binomial_2 K=1 and K=2", {
  fam <- gdpar_family("hurdle_neg_binomial_2")
  expect_error(
    gdpar:::.gdpar_guard_K_below_family_min(fam),
    class = "gdpar_unsupported_feature_error"
  )
  fam_p2 <- gdpar:::.gdpar_promote_scope_per_observation(fam,
                                                        c("mu", "phi"))
  expect_error(
    gdpar:::.gdpar_guard_K_below_family_min(fam_p2),
    class = "gdpar_unsupported_feature_error"
  )
  fam_p3 <- gdpar:::.gdpar_promote_scope_per_observation(fam,
                                                        c("mu", "phi", "pi"))
  expect_silent(gdpar:::.gdpar_guard_K_below_family_min(fam_p3))
})

# -----------------------------------------------------------------
# 3. .assemble_stan_data_K accepts each mixture family
# -----------------------------------------------------------------

test_that(".assemble_stan_data_K accepts ZIP K=2 and emits y_int", {
  set.seed(101L)
  n <- 40L
  df <- data.frame(
    y = c(rep(0L, 10L), rpois(n - 10L, lambda = 3)),
    x = rnorm(n), z = rnorm(n)
  )
  fam <- gdpar_family("zip")
  fam_p <- gdpar:::.gdpar_promote_scope_per_observation(fam,
                                                       c("mu", "pi"))
  fs <- gdpar_bf(y ~ a(x), pi ~ a(z))
  amm_list <- gdpar:::.gdpar_formula_set_to_amm_spec_list(fs, NULL)
  rhs <- y ~ x + z + 0
  design_K <- gdpar:::.build_amm_design_K(amm_list, df, formula_rhs = rhs)
  stan_data <- gdpar:::.assemble_stan_data_K(
    design_K = design_K, family = fam_p, amm_list_canonical = amm_list,
    y = df$y, theta_anchor_K = c(mu = 0, pi = 0)
  )
  expect_equal(stan_data$family_id_k, c(10L, 10L))
  expect_equal(stan_data$K, 2L)
  expect_equal(length(stan_data$y_int), n)
  expect_true(all(stan_data$y_real == 0))
  expect_true(all(stan_data$y_int >= 0L))
})

test_that(".assemble_stan_data_K accepts ZINB K=3 and emits y_int", {
  set.seed(102L)
  n <- 40L
  df <- data.frame(
    y = c(rep(0L, 8L), rnbinom(n - 8L, mu = 4, size = 3)),
    x = rnorm(n), z = rnorm(n), w = rnorm(n)
  )
  fam <- gdpar_family("zinb")
  fam_p <- gdpar:::.gdpar_promote_scope_per_observation(fam,
                                                       c("mu", "phi", "pi"))
  fs <- gdpar_bf(y ~ a(x), phi ~ a(z), pi ~ a(w))
  amm_list <- gdpar:::.gdpar_formula_set_to_amm_spec_list(fs, NULL)
  rhs <- y ~ x + z + w + 0
  design_K <- gdpar:::.build_amm_design_K(amm_list, df, formula_rhs = rhs)
  stan_data <- gdpar:::.assemble_stan_data_K(
    design_K = design_K, family = fam_p, amm_list_canonical = amm_list,
    y = df$y, theta_anchor_K = c(mu = 0, phi = 0, pi = 0)
  )
  expect_equal(stan_data$family_id_k, c(11L, 11L, 11L))
  expect_equal(stan_data$K, 3L)
  expect_true(all(stan_data$y_real == 0))
})

test_that(".assemble_stan_data_K accepts hurdle_poisson K=2", {
  set.seed(103L)
  n <- 40L
  df <- data.frame(
    y = c(rep(0L, 12L), 1L + rpois(n - 12L, lambda = 2)),
    x = rnorm(n), z = rnorm(n)
  )
  fam <- gdpar_family("hurdle_poisson")
  fam_p <- gdpar:::.gdpar_promote_scope_per_observation(fam,
                                                       c("mu", "pi"))
  fs <- gdpar_bf(y ~ a(x), pi ~ a(z))
  amm_list <- gdpar:::.gdpar_formula_set_to_amm_spec_list(fs, NULL)
  rhs <- y ~ x + z + 0
  design_K <- gdpar:::.build_amm_design_K(amm_list, df, formula_rhs = rhs)
  stan_data <- gdpar:::.assemble_stan_data_K(
    design_K = design_K, family = fam_p, amm_list_canonical = amm_list,
    y = df$y, theta_anchor_K = c(mu = 0, pi = 0)
  )
  expect_equal(stan_data$family_id_k, c(12L, 12L))
})

test_that(".assemble_stan_data_K accepts hurdle_neg_binomial_2 K=3", {
  set.seed(104L)
  n <- 40L
  df <- data.frame(
    y = c(rep(0L, 10L), 1L + rnbinom(n - 10L, mu = 3, size = 2)),
    x = rnorm(n), z = rnorm(n), w = rnorm(n)
  )
  fam <- gdpar_family("hurdle_neg_binomial_2")
  fam_p <- gdpar:::.gdpar_promote_scope_per_observation(fam,
                                                       c("mu", "phi", "pi"))
  fs <- gdpar_bf(y ~ a(x), phi ~ a(z), pi ~ a(w))
  amm_list <- gdpar:::.gdpar_formula_set_to_amm_spec_list(fs, NULL)
  rhs <- y ~ x + z + w + 0
  design_K <- gdpar:::.build_amm_design_K(amm_list, df, formula_rhs = rhs)
  stan_data <- gdpar:::.assemble_stan_data_K(
    design_K = design_K, family = fam_p, amm_list_canonical = amm_list,
    y = df$y, theta_anchor_K = c(mu = 0, phi = 0, pi = 0)
  )
  expect_equal(stan_data$family_id_k, c(13L, 13L, 13L))
})

test_that(".assemble_stan_data_K rejects mixture y < 0", {
  set.seed(105L)
  n <- 30L
  df <- data.frame(
    y = c(-1L, rpois(n - 1L, lambda = 2)),
    x = rnorm(n), z = rnorm(n)
  )
  fam <- gdpar_family("zip")
  fam_p <- gdpar:::.gdpar_promote_scope_per_observation(fam,
                                                       c("mu", "pi"))
  fs <- gdpar_bf(y ~ a(x), pi ~ a(z))
  amm_list <- gdpar:::.gdpar_formula_set_to_amm_spec_list(fs, NULL)
  rhs <- y ~ x + z + 0
  design_K <- gdpar:::.build_amm_design_K(amm_list, df, formula_rhs = rhs)
  expect_error(
    gdpar:::.assemble_stan_data_K(
      design_K = design_K, family = fam_p, amm_list_canonical = amm_list,
      y = df$y, theta_anchor_K = c(mu = 0, pi = 0)
    ),
    class = "gdpar_input_error"
  )
})

test_that(".assemble_stan_data_K rejects mixture y with non-integer values", {
  set.seed(106L)
  n <- 30L
  df <- data.frame(
    y = c(0.5, rpois(n - 1L, lambda = 2)),
    x = rnorm(n), z = rnorm(n)
  )
  fam <- gdpar_family("hurdle_poisson")
  fam_p <- gdpar:::.gdpar_promote_scope_per_observation(fam,
                                                       c("mu", "pi"))
  fs <- gdpar_bf(y ~ a(x), pi ~ a(z))
  amm_list <- gdpar:::.gdpar_formula_set_to_amm_spec_list(fs, NULL)
  rhs <- y ~ x + z + 0
  design_K <- gdpar:::.build_amm_design_K(amm_list, df, formula_rhs = rhs)
  expect_error(
    gdpar:::.assemble_stan_data_K(
      design_K = design_K, family = fam_p, amm_list_canonical = amm_list,
      y = df$y, theta_anchor_K = c(mu = 0, pi = 0)
    ),
    class = "gdpar_input_error"
  )
})

# -----------------------------------------------------------------
# 4. Custom-family blacklist: mixture names collide with built-ins
# -----------------------------------------------------------------

test_that("gdpar_family_custom rejects each mixture name as builtin collision", {
  for (nm in c("zip", "zinb", "hurdle_poisson",
               "hurdle_neg_binomial_2")) {
    expect_error(
      gdpar_family_custom(
        name               = nm,
        link               = "log",
        did_holds          = TRUE,
        did_condition      = NA_character_,
        stan_loglik_block  = "target += normal_lpdf(y_real[i] | eta[i], 1);",
        stan_log_lik_block = "log_lik[i] = normal_lpdf(y_real[i] | eta[i], 1);",
        stan_y_pred_block  = "y_pred[i] = normal_rng(eta[i], 1);",
        y_type             = "real",
        did_reference      = "test"
      ),
      class = "gdpar_input_error",
      info = sprintf("custom should reject '%s' as builtin", nm)
    )
  }
})

test_that("gdpar_family_custom_K rejects each mixture name as builtin collision", {
  for (nm in c("zip", "zinb", "hurdle_poisson",
               "hurdle_neg_binomial_2")) {
    expect_error(
      gdpar_family_custom_K(
        name         = nm,
        stan_lpdf_id = "lognormal_loc_scale",
        did_holds    = TRUE
      ),
      class = "gdpar_input_error",
      info = sprintf("custom_K should reject '%s' as builtin", nm)
    )
  }
})

# -----------------------------------------------------------------
# 5. Canonical prior for logit_pi (already registered)
# -----------------------------------------------------------------

test_that(".gdpar_canonical_prior_for_kind('logit_pi') returns normal(0, 2.5)", {
  expect_equal(gdpar:::.gdpar_canonical_prior_for_kind("logit_pi"),
               "normal(0, 2.5)")
})

test_that(".gdpar_canonical_prior_variance('logit_pi') returns 6.25", {
  expect_equal(gdpar:::.gdpar_canonical_prior_variance("logit_pi"),
               6.25)
})

# -----------------------------------------------------------------
# 6. generate_stan_code_K wires each mixture branch with no
#    unsubstituted placeholders. Decision D6 of 8.3.6 (a): the default
#    vectorized THETA_REF_PRIOR_BLOCK is used for mixtures (the logit
#    space of pi is covered by normal(0, 2.5) without structural
#    bounding); generalization of the helper is deferred to Session
#    8.4 per gotcha 55.
# -----------------------------------------------------------------

test_that("generate_stan_code_K leaves no placeholders for ZIP and contains branch 10", {
  fam <- gdpar_family("zip")
  src <- gdpar:::generate_stan_code_K(gdpar_prior(), family = fam)
  expect_false(grepl("{{", src, fixed = TRUE))
  expect_true(grepl("family_id_k[1] == 10", src, fixed = TRUE))
  expect_true(grepl("bernoulli_logit_lpmf", src, fixed = TRUE))
  expect_true(grepl("poisson_log_lpmf", src, fixed = TRUE))
  expect_true(grepl("theta_ref_k[1] ~ normal", src, fixed = TRUE))
})

test_that("generate_stan_code_K leaves no placeholders for ZINB and contains branch 11", {
  fam <- gdpar_family("zinb")
  src <- gdpar:::generate_stan_code_K(gdpar_prior(), family = fam)
  expect_false(grepl("{{", src, fixed = TRUE))
  expect_true(grepl("family_id_k[1] == 11", src, fixed = TRUE))
  expect_true(grepl("neg_binomial_2_log_lpmf", src, fixed = TRUE))
})

test_that("generate_stan_code_K leaves no placeholders for hurdle_poisson and contains branch 12", {
  fam <- gdpar_family("hurdle_poisson")
  src <- gdpar:::generate_stan_code_K(gdpar_prior(), family = fam)
  expect_false(grepl("{{", src, fixed = TRUE))
  expect_true(grepl("family_id_k[1] == 12", src, fixed = TRUE))
  expect_true(grepl("log1m_exp(-exp(eta_mu_i))", src, fixed = TRUE))
})

test_that("generate_stan_code_K leaves no placeholders for hurdle_neg_binomial_2 and contains branch 13", {
  fam <- gdpar_family("hurdle_neg_binomial_2")
  src <- gdpar:::generate_stan_code_K(gdpar_prior(), family = fam)
  expect_false(grepl("{{", src, fixed = TRUE))
  expect_true(grepl("family_id_k[1] == 13", src, fixed = TRUE))
  expect_true(grepl("log1m_exp(neg_binomial_2_log_lpmf(0", src,
                    fixed = TRUE))
})

test_that(".gdpar_build_theta_ref_prior_block(mixtures) returns the default block", {
  for (nm in c("zip", "zinb", "hurdle_poisson",
               "hurdle_neg_binomial_2")) {
    fam <- gdpar_family(nm)
    block <- gdpar:::.gdpar_build_theta_ref_prior_block(fam)
    expect_match(block, "theta_ref_k\\[1\\] ~ \\{\\{PRIOR_THETA_REF\\}\\};",
                 fixed = FALSE,
                 info = sprintf("mixture '%s' should use default block", nm))
    expect_false(grepl("uniform(1.01, 1.99)", block, fixed = TRUE),
                 info = sprintf("mixture '%s' must not carry tweedie slot prior",
                                nm))
  }
})

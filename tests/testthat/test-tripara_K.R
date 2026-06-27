# Sub-phase 8.3.5a (2026-05-21): structural tests for K = 3 Student-t
# family + did_override plasticity argument (D2). The K = 3 smoke fit
# lives in test-fit_K_smoke.R gated by GDPAR_RUN_STAN_SMOKE_STUDENT_T.

# -----------------------------------------------------------------
# 1. Student-t family construction
# -----------------------------------------------------------------

test_that("gdpar_family('student_t') constructs (mu identity, sigma log, nu log) param_specs", {
  fam <- gdpar_family("student_t")
  expect_s3_class(fam, "gdpar_family")
  expect_equal(fam$name, "student_t")
  expect_equal(fam$link, "identity")
  expect_equal(fam$stan_id, 8L)
  expect_true(fam$has_dispersion)
  expect_length(fam$param_specs, 3L)
  expect_equal(fam$param_specs[[1L]]$name, "mu")
  expect_equal(fam$param_specs[[1L]]$link, "identity")
  expect_equal(fam$param_specs[[1L]]$prior_canonical_kind, "mu")
  expect_equal(fam$param_specs[[1L]]$scope, "per_observation")
  expect_equal(fam$param_specs[[1L]]$family_role, "location")
  expect_equal(fam$param_specs[[1L]]$support, "real_line")
  expect_equal(fam$param_specs[[2L]]$name, "sigma")
  expect_equal(fam$param_specs[[2L]]$link, "log")
  expect_equal(fam$param_specs[[2L]]$prior_canonical_kind, "log_sigma")
  expect_equal(fam$param_specs[[2L]]$scope, "population")
  expect_equal(fam$param_specs[[2L]]$family_role, "scale")
  expect_equal(fam$param_specs[[2L]]$support, "positive_real")
  expect_equal(fam$param_specs[[3L]]$name, "nu")
  expect_equal(fam$param_specs[[3L]]$link, "log")
  expect_equal(fam$param_specs[[3L]]$prior_canonical_kind, "log_nu")
  expect_equal(fam$param_specs[[3L]]$scope, "population")
  expect_equal(fam$param_specs[[3L]]$family_role, "df")
  expect_equal(fam$param_specs[[3L]]$support, "positive_real")
})

test_that("gdpar_family('student_t') accepts log link override on mu", {
  fam <- gdpar_family("student_t", link = "log")
  expect_equal(fam$link, "log")
  expect_equal(fam$param_specs[[1L]]$link, "log")
  expect_equal(fam$param_specs[[2L]]$link, "log")
  expect_equal(fam$param_specs[[3L]]$link, "log")
})

test_that("gdpar_family('student_t') rejects non-canonical link", {
  expect_error(
    gdpar_family("student_t", link = "logit"),
    class = "gdpar_input_error"
  )
})

# -----------------------------------------------------------------
# 2. K < 3 guard for Student-t
# -----------------------------------------------------------------

test_that(".gdpar_guard_K_below_family_min fires for Student-t with K=1", {
  fam <- gdpar_family("student_t")
  expect_equal(gdpar:::.gdpar_n_params_individual(fam), 1L)
  expect_error(
    gdpar:::.gdpar_guard_K_below_family_min(fam),
    class = "gdpar_unsupported_feature_error"
  )
})

test_that(".gdpar_guard_K_below_family_min fires for Student-t with K=2 (mu + sigma promoted)", {
  fam <- gdpar_family("student_t")
  fam_p <- gdpar:::.gdpar_promote_scope_per_observation(fam,
                                                       c("mu", "sigma"))
  expect_equal(gdpar:::.gdpar_n_params_individual(fam_p), 2L)
  expect_error(
    gdpar:::.gdpar_guard_K_below_family_min(fam_p),
    class = "gdpar_unsupported_feature_error"
  )
})

test_that(".gdpar_guard_K_below_family_min passes when Student-t is elevated to K=3", {
  fam <- gdpar_family("student_t")
  fam_p <- gdpar:::.gdpar_promote_scope_per_observation(fam,
                                                       c("mu", "sigma", "nu"))
  expect_equal(gdpar:::.gdpar_n_params_individual(fam_p), 3L)
  expect_silent(gdpar:::.gdpar_guard_K_below_family_min(fam_p))
})

# -----------------------------------------------------------------
# 3. assemble_stan_data_K accepts Student-t K=3 (stan_id 8)
# -----------------------------------------------------------------

test_that(".assemble_stan_data_K accepts Student-t K=3 (stan_id 8) and emits y_real", {
  set.seed(11L)
  n <- 30L
  df <- data.frame(y = rnorm(n), x = rnorm(n), z = rnorm(n), w = rnorm(n))
  fam <- gdpar_family("student_t")
  fam_p <- gdpar:::.gdpar_promote_scope_per_observation(fam,
                                                       c("mu", "sigma", "nu"))
  fs <- gdpar_bf(y ~ a(x), sigma ~ a(z), nu ~ a(w))
  amm_list <- gdpar:::.gdpar_formula_set_to_amm_spec_list(fs, NULL)
  rhs <- y ~ x + z + w + 0
  design_K <- gdpar:::.build_amm_design_K(amm_list, df, formula_rhs = rhs)
  stan_data <- gdpar:::.assemble_stan_data_K(
    design_K = design_K, family = fam_p, amm_list_canonical = amm_list,
    y = df$y, theta_anchor_K = c(mu = 0, sigma = 0, nu = log(10))
  )
  expect_equal(stan_data$family_id_k, c(8L, 8L, 8L))
  expect_equal(stan_data$K, 3L)
  expect_equal(length(stan_data$y_real), n)
  expect_true(all(stan_data$y_int == 0L))
})

test_that(".assemble_stan_data_K imposes no y-range constraint for Student-t", {
  set.seed(12L)
  n <- 30L
  df <- data.frame(y = c(rnorm(n - 1L), -1e6),
                   x = rnorm(n), z = rnorm(n), w = rnorm(n))
  fam <- gdpar_family("student_t")
  fam_p <- gdpar:::.gdpar_promote_scope_per_observation(fam,
                                                       c("mu", "sigma", "nu"))
  fs <- gdpar_bf(y ~ a(x), sigma ~ a(z), nu ~ a(w))
  amm_list <- gdpar:::.gdpar_formula_set_to_amm_spec_list(fs, NULL)
  rhs <- y ~ x + z + w + 0
  design_K <- gdpar:::.build_amm_design_K(amm_list, df, formula_rhs = rhs)
  expect_silent(
    gdpar:::.assemble_stan_data_K(
      design_K = design_K, family = fam_p, amm_list_canonical = amm_list,
      y = df$y, theta_anchor_K = c(mu = 0, sigma = 0, nu = log(10))
    )
  )
})

# -----------------------------------------------------------------
# 4. Custom family blacklist: student_t collision
# -----------------------------------------------------------------

test_that("gdpar_family_custom rejects 'student_t' as builtin name collision", {
  expect_error(
    gdpar_family_custom(
      name               = "student_t",
      link               = "identity",
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

test_that("gdpar_family_custom_K rejects 'student_t' as builtin name collision", {
  expect_error(
    gdpar_family_custom_K(
      name         = "student_t",
      stan_lpdf_id = "lognormal_loc_scale",
      did_holds    = TRUE
    ),
    class = "gdpar_input_error"
  )
})

# -----------------------------------------------------------------
# 5. .gdpar_canonical_prior_variance covers log_nu (already registered)
# -----------------------------------------------------------------

test_that(".gdpar_canonical_prior_variance(log_nu) returns 1.0", {
  expect_equal(gdpar:::.gdpar_canonical_prior_variance("log_nu"), 1.0)
})

# -----------------------------------------------------------------
# 6. did_override plasticity (D2 of 8.3.5a)
# -----------------------------------------------------------------

test_that("did_override = NULL preserves canonical D-ID", {
  fam <- gdpar_family("gaussian")
  fam_default <- gdpar_family("gaussian", did_override = NULL)
  expect_equal(fam_default$param_specs[[1L]]$did_status,
               fam$param_specs[[1L]]$did_status)
  expect_equal(fam_default$param_specs[[2L]]$did_status,
               fam$param_specs[[2L]]$did_status)
  expect_equal(fam_default$param_specs[[2L]]$did_reference,
               fam$param_specs[[2L]]$did_reference)
})

test_that("did_override applies to a single named slot (Gaussian sigma)", {
  fam <- gdpar_family("gaussian", did_override = list(
    sigma = list(did_status    = "user_responsible",
                 did_condition = NA_character_,
                 did_reference = "User declaration")
  ))
  expect_equal(fam$param_specs[[1L]]$did_status, "holds")
  expect_equal(fam$param_specs[[2L]]$did_status, "user_responsible")
  expect_true(is.na(fam$param_specs[[2L]]$did_condition))
  expect_equal(fam$param_specs[[2L]]$did_reference, "User declaration")
})

test_that("did_override on the first slot propagates to family$did_status", {
  fam <- gdpar_family("student_t", did_override = list(
    mu = list(did_status    = "holds_under_condition",
              did_condition = "design-specific anchor identification")
  ))
  expect_equal(fam$param_specs[[1L]]$did_status, "holds_under_condition")
  expect_equal(fam$param_specs[[1L]]$did_condition,
               "design-specific anchor identification")
  expect_equal(fam$did_status, "holds_under_condition")
  expect_equal(fam$did_condition,
               "design-specific anchor identification")
})

test_that("did_override applies to multiple slots independently", {
  fam <- gdpar_family("student_t", did_override = list(
    sigma = list(did_status = "user_responsible"),
    nu    = list(did_reference = "Custom citation")
  ))
  expect_equal(fam$param_specs[[2L]]$did_status, "user_responsible")
  expect_equal(fam$param_specs[[3L]]$did_status, "holds_under_condition")
  expect_equal(fam$param_specs[[3L]]$did_reference, "Custom citation")
})

test_that("did_override rejects unknown slot names", {
  expect_error(
    gdpar_family("gaussian", did_override = list(
      not_a_slot = list(did_status = "holds")
    )),
    class = "gdpar_input_error"
  )
})

test_that("did_override rejects unknown fields", {
  expect_error(
    gdpar_family("gaussian", did_override = list(
      mu = list(weird_field = "x")
    )),
    class = "gdpar_input_error"
  )
})

test_that("did_override rejects invalid did_status values", {
  expect_error(
    gdpar_family("gaussian", did_override = list(
      mu = list(did_status = "definitely_not_valid")
    )),
    class = "gdpar_input_error"
  )
})

test_that("did_override rejects non-list top-level argument", {
  expect_error(
    gdpar_family("gaussian", did_override = c(mu = "x")),
    class = "gdpar_input_error"
  )
})

test_that("did_override rejects unnamed entries at top level", {
  expect_error(
    gdpar_family("gaussian",
                 did_override = list(list(did_status = "holds"))),
    class = "gdpar_input_error"
  )
})

test_that("did_override entry must be a named list (not a character)", {
  expect_error(
    gdpar_family("gaussian", did_override = list(mu = "holds")),
    class = "gdpar_input_error"
  )
})

test_that("did_override covers all three Student-t slots", {
  fam <- gdpar_family("student_t", did_override = list(
    mu    = list(did_reference = "Citation A"),
    sigma = list(did_status    = "user_responsible",
                 did_condition = NA_character_,
                 did_reference = "Citation B"),
    nu    = list(did_condition = "tail thickness assumed identified")
  ))
  expect_equal(fam$param_specs[[1L]]$did_reference, "Citation A")
  expect_equal(fam$param_specs[[2L]]$did_status, "user_responsible")
  expect_true(is.na(fam$param_specs[[2L]]$did_condition))
  expect_equal(fam$param_specs[[2L]]$did_reference, "Citation B")
  expect_equal(fam$param_specs[[3L]]$did_condition,
               "tail thickness assumed identified")
})

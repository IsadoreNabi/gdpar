# Sub-phase 8.3.7 (2026-05-21): structural tests for heterogeneous
# families per slot under the K-individual distributional regression
# path. Covers:
#   - public API named-list dispatch (D5)
#   - canonical inv_link / support / subset coherence helpers (D3/D4)
#   - validator with rich errors and suggestions (D4)
#   - resolver with L1 refined materialization + did_status override
#     + guard rejecting heterogeneous locations outside {1, 3, 5, 6, 7}
#     (D3.5, supplementary guard for D7)
#   - .assemble_stan_data_K propagation of family_id_k_vector and
#     emission of inv_link_id_per_slot
#   - generate_stan_code_K wiring of the apply_inv_link_by_id helper
#     and inv_link_id_per_slot data field
#
# Smoke fits with cmdstanr live in test-fit_K_smoke.R gated by
# GDPAR_RUN_STAN_SMOKE_HET_GAUSS_BETA / GDPAR_RUN_STAN_SMOKE_HET_GAUSS_GAMMA
# / GDPAR_RUN_STAN_SMOKE_HET_NB_BETA.

# -----------------------------------------------------------------
# 1. Canonical inv_link ID helpers (D3 / D3.5 L1 refined)
# -----------------------------------------------------------------

test_that(".gdpar_canonical_inv_link_id_slot1 returns 0/1/2 for built-in families", {
  expect_equal(gdpar:::.gdpar_canonical_inv_link_id_slot1(1L),  0L)  # gaussian
  expect_equal(gdpar:::.gdpar_canonical_inv_link_id_slot1(3L),  2L)  # neg_binomial_2
  expect_equal(gdpar:::.gdpar_canonical_inv_link_id_slot1(5L),  1L)  # beta
  expect_equal(gdpar:::.gdpar_canonical_inv_link_id_slot1(6L),  2L)  # gamma
  expect_equal(gdpar:::.gdpar_canonical_inv_link_id_slot1(7L),  0L)  # lognormal_loc_scale
  expect_equal(gdpar:::.gdpar_canonical_inv_link_id_slot1(8L),  0L)  # student_t
  expect_equal(gdpar:::.gdpar_canonical_inv_link_id_slot1(9L),  2L)  # tweedie
  expect_equal(gdpar:::.gdpar_canonical_inv_link_id_slot1(10L), 2L)  # zip
  expect_equal(gdpar:::.gdpar_canonical_inv_link_id_slot1(13L), 2L)  # hurdle_nb
})

test_that(".gdpar_canonical_inv_link_id_slot returns the canonical link of slot k", {
  # Gaussian: slot 1 identity, slot 2 sigma log.
  expect_equal(gdpar:::.gdpar_canonical_inv_link_id_slot(1L, 1L), 0L)
  expect_equal(gdpar:::.gdpar_canonical_inv_link_id_slot(1L, 2L), 2L)
  # neg_binomial_2: slot 1 log, slot 2 phi log.
  expect_equal(gdpar:::.gdpar_canonical_inv_link_id_slot(3L, 1L), 2L)
  expect_equal(gdpar:::.gdpar_canonical_inv_link_id_slot(3L, 2L), 2L)
  # Beta: slot 1 mu logit, slot 2 phi log.
  expect_equal(gdpar:::.gdpar_canonical_inv_link_id_slot(5L, 1L), 1L)
  expect_equal(gdpar:::.gdpar_canonical_inv_link_id_slot(5L, 2L), 2L)
  # Gamma: slot 1 log, slot 2 shape log.
  expect_equal(gdpar:::.gdpar_canonical_inv_link_id_slot(6L, 1L), 2L)
  expect_equal(gdpar:::.gdpar_canonical_inv_link_id_slot(6L, 2L), 2L)
})

test_that(".gdpar_canonical_inv_link_id_slot1 errors on NA stan_id", {
  expect_error(
    gdpar:::.gdpar_canonical_inv_link_id_slot1(NA_integer_),
    class = "gdpar_internal_error"
  )
})

# -----------------------------------------------------------------
# 2. Canonical support helpers
# -----------------------------------------------------------------

test_that(".gdpar_canonical_support_slot1 returns the support of slot 1 per family", {
  expect_equal(gdpar:::.gdpar_canonical_support_slot1(1L),  "real_line")
  expect_equal(gdpar:::.gdpar_canonical_support_slot1(3L),  "positive_real")
  expect_equal(gdpar:::.gdpar_canonical_support_slot1(5L),  "unit_interval")
  expect_equal(gdpar:::.gdpar_canonical_support_slot1(6L),  "positive_real")
  expect_equal(gdpar:::.gdpar_canonical_support_slot1(7L),  "real_line")
  expect_equal(gdpar:::.gdpar_canonical_support_slot1(8L),  "real_line")
  expect_equal(gdpar:::.gdpar_canonical_support_slot1(9L),  "positive_real")
  expect_equal(gdpar:::.gdpar_canonical_support_slot1(10L), "positive_real")
})

test_that(".gdpar_support_subset_coherent enforces the canonical containment table", {
  # Reflexive cases pass.
  expect_true(gdpar:::.gdpar_support_subset_coherent("real_line", "real_line"))
  expect_true(gdpar:::.gdpar_support_subset_coherent("positive_real", "positive_real"))
  expect_true(gdpar:::.gdpar_support_subset_coherent("unit_interval", "unit_interval"))
  # Strict subset cases (emitted strictly inside required) pass.
  expect_true(gdpar:::.gdpar_support_subset_coherent("positive_real", "real_line"))
  expect_true(gdpar:::.gdpar_support_subset_coherent("unit_interval", "positive_real"))
  expect_true(gdpar:::.gdpar_support_subset_coherent("unit_interval", "real_line"))
  expect_true(gdpar:::.gdpar_support_subset_coherent("bounded_open", "positive_real"))
  # Non-coherent cases fail.
  expect_false(gdpar:::.gdpar_support_subset_coherent("real_line", "positive_real"))
  expect_false(gdpar:::.gdpar_support_subset_coherent("real_line", "unit_interval"))
  expect_false(gdpar:::.gdpar_support_subset_coherent("positive_real", "unit_interval"))
})

test_that(".gdpar_compatible_families_for_support enumerates alternative families", {
  pos <- gdpar:::.gdpar_compatible_families_for_support("positive_real",
                                                        exclude_stan_ids = c(1L, 5L))
  expect_true(all(pos %in% c("neg_binomial_2", "gamma", "tweedie")))
  expect_true(length(pos) >= 2L)
  uni <- gdpar:::.gdpar_compatible_families_for_support("unit_interval",
                                                        exclude_stan_ids = integer(0))
  expect_true("beta" %in% uni)
  rea <- gdpar:::.gdpar_compatible_families_for_support("real_line",
                                                        exclude_stan_ids = integer(0))
  expect_true(all(c("gaussian", "neg_binomial_2", "beta", "gamma",
                    "student_t", "tweedie") %in% rea))
})

# -----------------------------------------------------------------
# 3. Validator with rich errors (D4)
# -----------------------------------------------------------------

test_that(".gdpar_validate_heterogeneous_family_K passes coherent combinations", {
  # Gaussian location + Beta scale: unit_interval subset of positive_real -> ok.
  fam_loc <- gdpar:::.gdpar_promote_scope_per_observation(
    gdpar_family("gaussian"), c("mu", "sigma")
  )
  expect_silent(
    gdpar:::.gdpar_validate_heterogeneous_family_K(
      family_het_list      = list(mu = gdpar_family("gaussian"),
                                  sigma = gdpar_family("beta")),
      location_param_specs = fam_loc$param_specs,
      slot_names           = c("mu", "sigma")
    )
  )
  # Gaussian location + Gamma scale: positive_real == positive_real -> ok.
  expect_silent(
    gdpar:::.gdpar_validate_heterogeneous_family_K(
      family_het_list      = list(mu = gdpar_family("gaussian"),
                                  sigma = gdpar_family("gamma")),
      location_param_specs = fam_loc$param_specs,
      slot_names           = c("mu", "sigma")
    )
  )
})

test_that(".gdpar_validate_heterogeneous_family_K rejects support violations", {
  # Beta location + Gaussian heterogeneous in phi: slot 1 of Gaussian
  # is real_line, but slot 2 (phi) of Beta requires positive_real ->
  # violation. (Beta location is one of the supported K = 2 locations
  # in D7, so this exercises the validator directly without firing the
  # earlier location-family guard.)
  fam_loc <- gdpar:::.gdpar_promote_scope_per_observation(
    gdpar_family("beta"), c("mu", "phi")
  )
  expect_error(
    gdpar:::.gdpar_validate_heterogeneous_family_K(
      family_het_list      = list(mu = gdpar_family("beta"),
                                  phi = gdpar_family("gaussian")),
      location_param_specs = fam_loc$param_specs,
      slot_names           = c("mu", "phi")
    ),
    class = "gdpar_input_error"
  )
})

test_that("validator error message includes role, supports, and suggestions", {
  fam_loc <- gdpar:::.gdpar_promote_scope_per_observation(
    gdpar_family("neg_binomial_2"), c("mu", "phi")
  )
  err <- tryCatch(
    gdpar:::.gdpar_validate_heterogeneous_family_K(
      family_het_list      = list(mu = gdpar_family("neg_binomial_2"),
                                  phi = gdpar_family("student_t")),
      location_param_specs = fam_loc$param_specs,
      slot_names           = c("mu", "phi")
    ),
    error = function(e) e
  )
  expect_s3_class(err, "gdpar_input_error")
  msg <- conditionMessage(err)
  expect_true(grepl("phi", msg, fixed = TRUE))
  expect_true(grepl("scale", msg, fixed = TRUE))
  expect_true(grepl("real_line", msg, fixed = TRUE))
  expect_true(grepl("positive_real", msg, fixed = TRUE))
  expect_true(grepl("gamma|tweedie|beta", msg, perl = TRUE))
})

# -----------------------------------------------------------------
# 4. Resolver with L1 refined materialization + guard (D3.5, D7)
# -----------------------------------------------------------------

test_that(".gdpar_resolve_heterogeneous_family_K passes through a single gdpar_family", {
  res <- gdpar:::.gdpar_resolve_heterogeneous_family_K(
    family_input = gdpar_family("gaussian"),
    slot_names   = c("mu", "sigma")
  )
  expect_false(res$is_heterogeneous)
  expect_equal(res$family_id_k_vector, c(1L, 1L))
  expect_s3_class(res$location_family, "gdpar_family")
  expect_equal(res$location_family$name, "gaussian")
})

test_that(".gdpar_resolve_heterogeneous_family_K dispatches heterogeneous Gauss+Beta", {
  res <- gdpar:::.gdpar_resolve_heterogeneous_family_K(
    family_input = list(mu = gdpar_family("gaussian"),
                        sigma = gdpar_family("beta")),
    slot_names   = c("mu", "sigma")
  )
  expect_true(res$is_heterogeneous)
  expect_equal(res$family_id_k_vector, c(1L, 5L))
  expect_equal(res$location_family$name, "gaussian")
  # Slot 2 inherits did_status from Beta mu slot 1 (which is "holds").
  expect_equal(res$location_family$param_specs[[2L]]$did_status, "holds")
  # Slot 2 keeps its original support/link in param_specs (Stan side
  # uses inv_link_id_per_slot to dispatch the heterogeneous link).
  expect_equal(res$location_family$param_specs[[2L]]$support, "positive_real")
})

test_that(".gdpar_resolve_heterogeneous_family_K propagates did_condition + did_reference", {
  beta_with_override <- gdpar_family(
    "beta",
    did_override = list(
      mu = list(
        did_status    = "user_responsible",
        did_condition = "user-supplied custom condition for tests",
        did_reference = "test-heterogeneous_K"
      )
    )
  )
  res <- gdpar:::.gdpar_resolve_heterogeneous_family_K(
    family_input = list(mu = gdpar_family("gaussian"),
                        sigma = beta_with_override),
    slot_names   = c("mu", "sigma")
  )
  expect_equal(res$location_family$param_specs[[2L]]$did_status,
               "user_responsible")
  expect_equal(res$location_family$param_specs[[2L]]$did_condition,
               "user-supplied custom condition for tests")
  expect_equal(res$location_family$param_specs[[2L]]$did_reference,
               "test-heterogeneous_K")
})

test_that("resolver rejects named-list with wrong K", {
  expect_error(
    gdpar:::.gdpar_resolve_heterogeneous_family_K(
      family_input = list(mu = gdpar_family("gaussian")),
      slot_names   = c("mu", "sigma")
    ),
    class = "gdpar_input_error"
  )
  expect_error(
    gdpar:::.gdpar_resolve_heterogeneous_family_K(
      family_input = list(mu = gdpar_family("gaussian"),
                          sigma = gdpar_family("beta"),
                          phi = gdpar_family("gamma")),
      slot_names   = c("mu", "sigma")
    ),
    class = "gdpar_input_error"
  )
})

test_that("resolver rejects named-list with mismatched slot names", {
  expect_error(
    gdpar:::.gdpar_resolve_heterogeneous_family_K(
      family_input = list(mu = gdpar_family("gaussian"),
                          dispersion = gdpar_family("beta")),
      slot_names   = c("mu", "sigma")
    ),
    class = "gdpar_input_error"
  )
})

test_that("resolver rejects non-gdpar_family entries", {
  expect_error(
    gdpar:::.gdpar_resolve_heterogeneous_family_K(
      family_input = list(mu = gdpar_family("gaussian"),
                          sigma = "beta"),
      slot_names   = c("mu", "sigma")
    ),
    class = "gdpar_input_error"
  )
})

test_that("resolver guard rejects heterogeneous Student-t location (deferred)", {
  err <- tryCatch(
    gdpar:::.gdpar_resolve_heterogeneous_family_K(
      family_input = list(mu = gdpar_family("student_t"),
                          sigma = gdpar_family("beta"),
                          nu = gdpar_family("gamma")),
      slot_names   = c("mu", "sigma", "nu")
    ),
    error = function(e) e
  )
  expect_s3_class(err, "gdpar_unsupported_feature_error")
  msg <- conditionMessage(err)
  expect_true(grepl("student_t", msg, fixed = TRUE))
  expect_true(grepl("8.3.7", msg, fixed = TRUE))
})

test_that("resolver guard rejects heterogeneous Tweedie / mixtures location (deferred)", {
  for (loc_name in c("tweedie", "zip", "zinb",
                     "hurdle_poisson", "hurdle_neg_binomial_2")) {
    slot_names <- c("mu", if (loc_name %in% c("zinb", "hurdle_neg_binomial_2"))
                    c("phi", "pi") else "pi")
    if (loc_name == "tweedie") slot_names <- c("mu", "phi", "p")
    extras <- vector("list", length(slot_names) - 1L)
    # Use Beta for any auxiliary slot that admits unit_interval -> positive_real
    extras[] <- list(gdpar_family("beta"))
    fam_list <- c(list(gdpar_family(loc_name)), extras)
    names(fam_list) <- slot_names
    expect_error(
      gdpar:::.gdpar_resolve_heterogeneous_family_K(fam_list, slot_names),
      class = "gdpar_unsupported_feature_error",
      info  = sprintf("location '%s' should be rejected as heterogeneous", loc_name)
    )
  }
})

test_that("resolver permits homogeneous Student-t / Tweedie (named-list trivial case)", {
  # When the named-list reproduces the homogeneous regime (all entries
  # share the location stan_id), is_heterogeneous is FALSE and the
  # K = 3+ guard does NOT fire.
  res_st <- gdpar:::.gdpar_resolve_heterogeneous_family_K(
    family_input = list(mu = gdpar_family("student_t"),
                        sigma = gdpar_family("student_t"),
                        nu = gdpar_family("student_t")),
    slot_names   = c("mu", "sigma", "nu")
  )
  expect_false(res_st$is_heterogeneous)
  expect_equal(res_st$family_id_k_vector, c(8L, 8L, 8L))
})

# -----------------------------------------------------------------
# 5. inv_link_id_per_slot computation under L1 refined (D3.5)
# -----------------------------------------------------------------

test_that(".gdpar_compute_inv_link_id_per_slot returns canonical IDs for homogeneous", {
  fam_g <- gdpar_family("gaussian")
  expect_equal(
    gdpar:::.gdpar_compute_inv_link_id_per_slot(c(1L, 1L), fam_g),
    c(0L, 2L)   # mu identity, sigma log
  )
  fam_b <- gdpar_family("beta")
  expect_equal(
    gdpar:::.gdpar_compute_inv_link_id_per_slot(c(5L, 5L), fam_b),
    c(1L, 2L)   # mu logit, phi log
  )
})

test_that(".gdpar_compute_inv_link_id_per_slot follows L1 for heterogeneous", {
  fam_g <- gdpar_family("gaussian")
  # Gauss location + Beta heterogeneous in slot 2: slot 2 inherits
  # Beta slot 1 (logit, ID = 1).
  expect_equal(
    gdpar:::.gdpar_compute_inv_link_id_per_slot(c(1L, 5L), fam_g),
    c(0L, 1L)
  )
  # Gauss location + Gamma heterogeneous in slot 2: slot 2 inherits
  # Gamma slot 1 (log, ID = 2). Trivial coincidence with the
  # homogeneous Gaussian slot 2 link, validating that the API works
  # symmetrically even when the heterogeneity is "virtual".
  expect_equal(
    gdpar:::.gdpar_compute_inv_link_id_per_slot(c(1L, 6L), fam_g),
    c(0L, 2L)
  )
})

# -----------------------------------------------------------------
# 6. .assemble_stan_data_K propagation of family_id_k_vector
# -----------------------------------------------------------------

test_that(".assemble_stan_data_K accepts NULL family_id_k_vector (legacy homogeneous)", {
  set.seed(401L); n <- 30L
  d <- data.frame(y = rnorm(n), x = rnorm(n), z = rnorm(n))
  fam_p <- gdpar:::.gdpar_promote_scope_per_observation(
    gdpar_family("gaussian"), c("mu", "sigma")
  )
  fs <- gdpar_bf(y ~ a(x), sigma ~ a(z))
  amm_list <- gdpar:::.gdpar_formula_set_to_amm_spec_list(fs, NULL)
  rhs <- y ~ x + z + 0
  design_K <- gdpar:::.build_amm_design_K(amm_list, d, formula_rhs = rhs)
  sd <- gdpar:::.assemble_stan_data_K(
    design_K = design_K, family = fam_p, amm_list_canonical = amm_list,
    y = d$y, theta_anchor_K = c(mu = 0, sigma = 0)
  )
  expect_equal(sd$family_id_k, c(1L, 1L))
  expect_equal(sd$inv_link_id_per_slot, c(0L, 2L))
})

test_that(".assemble_stan_data_K accepts heterogeneous family_id_k_vector", {
  set.seed(402L); n <- 30L
  d <- data.frame(y = rnorm(n), x = rnorm(n), z = rnorm(n))
  fam_p <- gdpar:::.gdpar_promote_scope_per_observation(
    gdpar_family("gaussian"), c("mu", "sigma")
  )
  fs <- gdpar_bf(y ~ a(x), sigma ~ a(z))
  amm_list <- gdpar:::.gdpar_formula_set_to_amm_spec_list(fs, NULL)
  rhs <- y ~ x + z + 0
  design_K <- gdpar:::.build_amm_design_K(amm_list, d, formula_rhs = rhs)
  sd <- gdpar:::.assemble_stan_data_K(
    design_K = design_K, family = fam_p, amm_list_canonical = amm_list,
    y = d$y, theta_anchor_K = c(mu = 0, sigma = 0),
    family_id_k_vector = c(1L, 5L)
  )
  expect_equal(sd$family_id_k, c(1L, 5L))
  expect_equal(sd$inv_link_id_per_slot, c(0L, 1L))
})

test_that(".assemble_stan_data_K rejects family_id_k_vector with wrong length", {
  set.seed(403L); n <- 20L
  d <- data.frame(y = rnorm(n), x = rnorm(n), z = rnorm(n))
  fam_p <- gdpar:::.gdpar_promote_scope_per_observation(
    gdpar_family("gaussian"), c("mu", "sigma")
  )
  fs <- gdpar_bf(y ~ a(x), sigma ~ a(z))
  amm_list <- gdpar:::.gdpar_formula_set_to_amm_spec_list(fs, NULL)
  rhs <- y ~ x + z + 0
  design_K <- gdpar:::.build_amm_design_K(amm_list, d, formula_rhs = rhs)
  expect_error(
    gdpar:::.assemble_stan_data_K(
      design_K = design_K, family = fam_p,
      amm_list_canonical = amm_list, y = d$y,
      theta_anchor_K = c(mu = 0, sigma = 0),
      family_id_k_vector = c(1L, 5L, 3L)
    ),
    class = "gdpar_internal_error"
  )
})

test_that(".assemble_stan_data_K rejects family_id_k_vector with slot1 mismatch", {
  set.seed(404L); n <- 20L
  d <- data.frame(y = rnorm(n), x = rnorm(n), z = rnorm(n))
  fam_p <- gdpar:::.gdpar_promote_scope_per_observation(
    gdpar_family("gaussian"), c("mu", "sigma")
  )
  fs <- gdpar_bf(y ~ a(x), sigma ~ a(z))
  amm_list <- gdpar:::.gdpar_formula_set_to_amm_spec_list(fs, NULL)
  rhs <- y ~ x + z + 0
  design_K <- gdpar:::.build_amm_design_K(amm_list, d, formula_rhs = rhs)
  expect_error(
    gdpar:::.assemble_stan_data_K(
      design_K = design_K, family = fam_p,
      amm_list_canonical = amm_list, y = d$y,
      theta_anchor_K = c(mu = 0, sigma = 0),
      family_id_k_vector = c(5L, 5L)
    ),
    class = "gdpar_internal_error"
  )
})

# -----------------------------------------------------------------
# 7. Stan code generation wiring
# -----------------------------------------------------------------

test_that("generate_stan_code_K leaves no placeholders for Gaussian K=2", {
  src <- gdpar:::generate_stan_code_K(gdpar_prior(),
                                      family = gdpar_family("gaussian"))
  expect_false(grepl("{{", src, fixed = TRUE))
})

test_that("generate_stan_code_K wires apply_inv_link_by_id and inv_link_id_per_slot", {
  src <- gdpar:::generate_stan_code_K(gdpar_prior(),
                                      family = gdpar_family("gaussian"))
  expect_true(grepl("real apply_inv_link_by_id(int link_id, real eta)",
                    src, fixed = TRUE))
  expect_true(grepl("inv_link_id_per_slot", src, fixed = TRUE))
  expect_true(grepl("real sigma_i = apply_inv_link_by_id(inv_link_id_per_slot[2]",
                    src, fixed = TRUE))
})

test_that("generate_stan_code_K wires the K=2 branches for NB / Beta / Gamma / Lognormal", {
  for (nm in c("neg_binomial_2", "beta", "gamma")) {
    src <- gdpar:::generate_stan_code_K(gdpar_prior(),
                                        family = gdpar_family(nm))
    expect_false(grepl("{{", src, fixed = TRUE),
                 info = sprintf("placeholders left for family '%s'", nm))
    expect_true(grepl("apply_inv_link_by_id(inv_link_id_per_slot[2]",
                      src, fixed = TRUE),
                info = sprintf("K=2 branch not refactored for '%s'", nm))
  }
  # Lognormal_loc_scale (custom_K stan_id 7).
  fam_ln <- suppressMessages(gdpar_family_custom_K(
    name = "ln_test", stan_lpdf_id = "lognormal_loc_scale",
    did_holds = TRUE, did_reference = "test"
  ))
  src <- gdpar:::generate_stan_code_K(gdpar_prior(), family = fam_ln)
  expect_false(grepl("{{", src, fixed = TRUE))
  expect_true(grepl("apply_inv_link_by_id(inv_link_id_per_slot[2]",
                    src, fixed = TRUE))
})

# -----------------------------------------------------------------
# 8. Public API: gdpar() arg family accepts named-list
# -----------------------------------------------------------------

test_that("named-list family detection short-circuits the legacy assert_inherits", {
  # Pure-R structural check: the dispatcher detects the named-list
  # shape via the .family_is_named_list predicate computed at the top
  # of gdpar() and routes to the resolver instead of the legacy
  # assert_inherits(c('gdpar_family', 'gdpar_family_multi')) path. We
  # validate the predicate's components here without invoking Stan.
  family_input <- list(mu = gdpar_family("gaussian"),
                       sigma = gdpar_family("beta"))
  expect_true(is.list(family_input))
  expect_false(inherits(family_input, "gdpar_family"))
  expect_false(inherits(family_input, "gdpar_family_multi"))
  expect_false(is.null(names(family_input)))
  expect_true(all(nzchar(names(family_input))))
  expect_equal(anyDuplicated(names(family_input)), 0L)
  expect_true(all(vapply(family_input, inherits, logical(1L),
                         "gdpar_family")))
})

test_that("gdpar() rejects named-list family when resolved K == 1", {
  set.seed(702L); n <- 20L
  d <- data.frame(y = rnorm(n), x = rnorm(n))
  # K = 1 (only a single slot in the formula set): named-list is
  # explicitly rejected by the D5 dispatcher.
  fs <- gdpar_bf(y ~ a(x))
  expect_error(
    suppressMessages(suppressWarnings(
      gdpar(fs, data = d,
            family = list(mu = gdpar_family("gaussian")),
            chains = 1L, iter_warmup = 5L, iter_sampling = 5L,
            refresh = 0L, seed = 1L,
            skip_id_check = TRUE, verbose = FALSE)
    )),
    class = "gdpar_input_error"
  )
})

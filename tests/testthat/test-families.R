test_that("gdpar_family constructs each built-in family", {
  for (name in c("gaussian", "poisson", "neg_binomial_2", "bernoulli")) {
    fam <- gdpar_family(name)
    expect_s3_class(fam, "gdpar_family")
    expect_equal(fam$name, name)
    expect_true(fam$did_status %in% c("holds", "holds_under_condition"))
  }
})

test_that("gdpar_family rejects unsupported link", {
  expect_error(
    gdpar_family("poisson", link = "identity"),
    class = "gdpar_input_error"
  )
})

test_that("inverse link is correct for built-in families", {
  expect_equal(gdpar_family("gaussian")$inv_link(2), 2)
  expect_equal(gdpar_family("poisson")$inv_link(0), 1)
  expect_equal(gdpar_family("bernoulli")$inv_link(0), 0.5)
})

.custom_args_min <- function(name = "my_test",
                              stan_loglik_block  =
                                "target += normal_lpdf(y_real[i] | eta[i], 1);",
                              stan_log_lik_block =
                                "log_lik[i] = normal_lpdf(y_real[i] | eta[i], 1);",
                              stan_y_pred_block  =
                                "y_pred[i] = normal_rng(eta[i], 1);",
                              y_type = "real") {
  list(
    name               = name,
    link               = "log",
    did_holds          = TRUE,
    did_condition      = NA_character_,
    stan_loglik_block  = stan_loglik_block,
    stan_log_lik_block = stan_log_lik_block,
    stan_y_pred_block  = stan_y_pred_block,
    y_type             = y_type,
    did_reference      = "test"
  )
}

test_that("gdpar_family_custom requires explicit did_holds", {
  args <- .custom_args_min(name = "x")
  args$did_holds <- NA
  expect_error(
    do.call(gdpar_family_custom, args),
    class = "gdpar_input_error"
  )
})

test_that("gdpar_family_custom rejects names of built-ins", {
  args <- .custom_args_min(name = "gaussian")
  expect_error(
    do.call(gdpar_family_custom, args),
    class = "gdpar_input_error"
  )
})

test_that("gdpar_family_custom emits a did message on construction", {
  args <- .custom_args_min()
  expect_message(
    do.call(gdpar_family_custom, args),
    class = "gdpar_did_message"
  )
})

test_that("gdpar_family_custom rejects legacy 'y[i]' placeholder", {
  args <- .custom_args_min(
    stan_loglik_block = "target += normal_lpdf(y[i] | eta[i], 1);"
  )
  expect_error(
    do.call(gdpar_family_custom, args),
    class = "gdpar_input_error"
  )
})

test_that("gdpar_family_custom rejects invalid y_type", {
  args <- .custom_args_min(y_type = "double")
  expect_error(
    do.call(gdpar_family_custom, args),
    class = "gdpar_input_error"
  )
})

test_that("gdpar_family_custom stores all three Stan snippets and y_type", {
  args <- .custom_args_min()
  fam <- suppressMessages(do.call(gdpar_family_custom, args))
  expect_equal(fam$stan_loglik_block, args$stan_loglik_block)
  expect_equal(fam$stan_log_lik_block, args$stan_log_lik_block)
  expect_equal(fam$stan_y_pred_block, args$stan_y_pred_block)
  expect_equal(fam$y_type, "real")
  expect_true(fam$is_custom)
})

test_that("gdpar_family_custom accepts y_type='integer' with y_int snippets", {
  args <- .custom_args_min(
    stan_loglik_block  = "target += poisson_log_lpmf(y_int[i] | eta[i]);",
    stan_log_lik_block = "log_lik[i] = poisson_log_lpmf(y_int[i] | eta[i]);",
    stan_y_pred_block  = "y_pred[i] = poisson_log_rng(eta[i]);",
    y_type             = "integer"
  )
  fam <- suppressMessages(do.call(gdpar_family_custom, args))
  expect_equal(fam$y_type, "integer")
})

test_that("param_specs of built-in families carry family_role", {
  expected_roles <- list(
    gaussian       = c("location", "scale"),
    poisson        = "location",
    neg_binomial_2 = c("location", "scale"),
    bernoulli      = "location"
  )
  for (name in names(expected_roles)) {
    fam <- gdpar_family(name)
    roles <- vapply(fam$param_specs,
                    function(s) s$family_role, character(1L))
    expect_equal(roles, expected_roles[[name]],
                 info = sprintf("family = %s", name))
  }
})

test_that("custom family param_spec carries family_role = 'location'", {
  args <- .custom_args_min()
  fam <- suppressMessages(do.call(gdpar_family_custom, args))
  expect_equal(fam$param_specs[[1L]]$family_role, "location")
})

test_that(".gdpar_known_family_roles enumerates the canonical roster", {
  expect_setequal(
    gdpar:::.gdpar_known_family_roles(),
    c("location", "scale", "shape", "df", "mixture_pi", "power")
  )
})

test_that("print.gdpar_param_spec shows family_role", {
  fam <- gdpar_family("gaussian")
  out <- capture.output(print(fam$param_specs[[1L]]))
  expect_true(any(grepl("family_role", out)))
  expect_true(any(grepl("location", out)))
})

test_that("print method for gdpar_family runs without error", {
  expect_output(print(gdpar_family("poisson")), "gdpar_family")
})

test_that("gdpar_family carries a param_specs list (decision 1C)", {
  fam <- gdpar_family("gaussian")
  expect_true(!is.null(fam$param_specs))
  expect_true(is.list(fam$param_specs))
  expect_length(fam$param_specs, 2L)
  for (spec in fam$param_specs) {
    expect_s3_class(spec, "gdpar_param_spec")
    expect_true(all(c("name", "link", "linkfun", "inv_link",
                      "did_status", "did_condition", "did_reference",
                      "support", "prior_canonical_kind", "scope") %in%
                    names(spec)))
  }
})

test_that("param_specs length matches family structural cardinality", {
  expected <- list(
    gaussian       = 2L,
    poisson        = 1L,
    neg_binomial_2 = 2L,
    bernoulli      = 1L
  )
  for (name in names(expected)) {
    fam <- gdpar_family(name)
    expect_length(fam$param_specs, expected[[name]])
  }
})

test_that("first param_spec is per-observation in every built-in", {
  for (name in c("gaussian", "poisson", "neg_binomial_2", "bernoulli")) {
    fam <- gdpar_family(name)
    expect_equal(fam$param_specs[[1L]]$scope, "per_observation")
    expect_equal(fam$param_specs[[1L]]$name, "mu")
  }
})

test_that("auxiliary param_specs are population-scoped in K=1 default", {
  gauss <- gdpar_family("gaussian")
  expect_equal(gauss$param_specs[[2L]]$name, "sigma")
  expect_equal(gauss$param_specs[[2L]]$scope, "population")
  expect_equal(gauss$param_specs[[2L]]$did_status, "holds_under_condition")
  expect_equal(gauss$param_specs[[2L]]$prior_canonical_kind, "log_sigma")
  nb <- gdpar_family("neg_binomial_2")
  expect_equal(nb$param_specs[[2L]]$name, "phi")
  expect_equal(nb$param_specs[[2L]]$scope, "population")
  expect_equal(nb$param_specs[[2L]]$prior_canonical_kind, "log_phi")
})

test_that("param_spec inv_link matches family inv_link for the first entry", {
  for (name in c("gaussian", "poisson", "neg_binomial_2", "bernoulli")) {
    fam <- gdpar_family(name)
    eta <- 0.7
    expect_equal(fam$inv_link(eta), fam$param_specs[[1L]]$inv_link(eta))
  }
})

test_that("gdpar_family_custom exposes a single param_spec", {
  fam <- suppressMessages(gdpar_family_custom(
    name              = "my_log_normal",
    link               = "log",
    did_holds          = TRUE,
    did_condition      = NA_character_,
    stan_loglik_block  =
      "target += normal_lpdf(log(y_real[i]) | eta[i], 1);",
    stan_log_lik_block =
      "log_lik[i] = normal_lpdf(log(y_real[i]) | eta[i], 1);",
    stan_y_pred_block  =
      "y_pred[i] = exp(normal_rng(eta[i], 1));",
    y_type             = "real",
    did_reference      = "test"
  ))
  expect_true(!is.null(fam$param_specs))
  expect_length(fam$param_specs, 1L)
  expect_s3_class(fam$param_specs[[1L]], "gdpar_param_spec")
  expect_equal(fam$param_specs[[1L]]$scope, "per_observation")
  expect_equal(fam$param_specs[[1L]]$support, "custom")
})

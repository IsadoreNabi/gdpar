## Tests for the AMM RHS parser and supporting helpers introduced in
## Unit 2 of sub-phase 8.3.3. Cover decisions P1-P5 of the parser
## scoping canonized in project_gdpar_block_8_3_extended_plan
## (2026-05-20).

test_that(".gdpar_rhs_has_amm_calls detects a(), b(), W() at top level", {
  expect_true(gdpar:::.gdpar_rhs_has_amm_calls(y ~ a(x)))
  expect_true(gdpar:::.gdpar_rhs_has_amm_calls(y ~ b(x)))
  expect_true(gdpar:::.gdpar_rhs_has_amm_calls(y ~ W()))
  expect_true(gdpar:::.gdpar_rhs_has_amm_calls(y ~ a(x) + b(z)))
  expect_true(gdpar:::.gdpar_rhs_has_amm_calls(y ~ a(x) + W()))
})

test_that(".gdpar_rhs_has_amm_calls returns FALSE on pure-symbol RHS", {
  expect_false(gdpar:::.gdpar_rhs_has_amm_calls(y ~ x1 + x2))
  expect_false(gdpar:::.gdpar_rhs_has_amm_calls(y ~ 1))
  expect_false(gdpar:::.gdpar_rhs_has_amm_calls(y ~ poly(x, 2)))
})

test_that(".gdpar_rhs_has_amm_calls is non-throwing on malformed RHS", {
  expect_false(gdpar:::.gdpar_rhs_has_amm_calls(y ~ a(x) - 1))
  expect_silent(gdpar:::.gdpar_rhs_has_amm_calls(y ~ a(x) * b(z)))
})

test_that(".gdpar_is_default_amm_spec detects the empty Level 0 spec", {
  expect_true(gdpar:::.gdpar_is_default_amm_spec(amm_spec()))
  expect_false(gdpar:::.gdpar_is_default_amm_spec(amm_spec(a = ~ x)))
  expect_false(gdpar:::.gdpar_is_default_amm_spec(amm_spec(b = ~ x)))
  expect_false(gdpar:::.gdpar_is_default_amm_spec(NULL))
  expect_false(gdpar:::.gdpar_is_default_amm_spec(list(a = NULL)))
})

test_that("parser extracts a(), b(), W() from a single-slot formula set", {
  fs <- gdpar_formula_set(mu = y ~ a(x1) + b(x2) + W())
  W <- W_basis(type = "polynomial", degree = 2L)
  out <- gdpar:::.gdpar_formula_set_to_amm_spec_list(fs, W_basis_arg = W)
  expect_named(out, "mu")
  expect_s3_class(out[["mu"]], "amm_spec")
  expect_equal(all.vars(out[["mu"]]$a), "x1")
  expect_equal(all.vars(out[["mu"]]$b), "x2")
  expect_s3_class(out[["mu"]]$W, "W_basis")
  expect_equal(out[["mu"]]$level, 2L)
})

test_that("parser tolerates expression-style RHS inside a() and b()", {
  fs <- gdpar_formula_set(mu = y ~ a(x1 + x2) + b(z1 + z1:z2))
  out <- gdpar:::.gdpar_formula_set_to_amm_spec_list(fs, NULL)
  expect_setequal(all.vars(out[["mu"]]$a), c("x1", "x2"))
  expect_setequal(all.vars(out[["mu"]]$b), c("z1", "z2"))
})

test_that("parser accepts ~ 1 as Level 0 AMM (anchor only)", {
  fs <- gdpar_formula_set(mu = y ~ 1)
  out <- gdpar:::.gdpar_formula_set_to_amm_spec_list(fs, NULL)
  expect_null(out[["mu"]]$a)
  expect_null(out[["mu"]]$b)
  expect_null(out[["mu"]]$W)
  expect_equal(out[["mu"]]$level, 0L)
})

test_that("parser rejects unknown functions in the RHS (P1)", {
  fs <- gdpar_formula_set(mu = y ~ a(x) + s(z))
  err <- expect_error(
    gdpar:::.gdpar_formula_set_to_amm_spec_list(fs, NULL),
    class = "gdpar_input_error"
  )
  expect_match(conditionMessage(err), "unknown function call `s\\(\\)`")
})

test_that("parser rejects multiple occurrences of the same wrapper (P2)", {
  fs_a <- gdpar_formula_set(mu = y ~ a(x1) + a(x2))
  expect_error(
    gdpar:::.gdpar_formula_set_to_amm_spec_list(fs_a, NULL),
    class = "gdpar_input_error"
  )
  fs_b <- gdpar_formula_set(mu = y ~ b(x1) + b(x2))
  expect_error(
    gdpar:::.gdpar_formula_set_to_amm_spec_list(fs_b, NULL),
    class = "gdpar_input_error"
  )
  fs_W <- gdpar_formula_set(mu = y ~ a(x) + W() + W())
  W <- W_basis(type = "polynomial", degree = 2L)
  expect_error(
    gdpar:::.gdpar_formula_set_to_amm_spec_list(fs_W, W),
    class = "gdpar_input_error"
  )
})

test_that("parser rejects W() with arguments (P3)", {
  fs <- gdpar_formula_set(mu = y ~ a(x) + W(degree = 2))
  err <- expect_error(
    gdpar:::.gdpar_formula_set_to_amm_spec_list(fs, NULL),
    class = "gdpar_input_error"
  )
  expect_match(conditionMessage(err), "zero")
})

test_that("parser errors when W() appears but external W arg is NULL", {
  fs <- gdpar_formula_set(mu = y ~ a(x) + W())
  err <- expect_error(
    gdpar:::.gdpar_formula_set_to_amm_spec_list(fs, W_basis_arg = NULL),
    class = "gdpar_input_error"
  )
  expect_match(conditionMessage(err), "external")
})

test_that("parser rejects a non-W_basis external W argument", {
  fs <- gdpar_formula_set(mu = y ~ a(x) + W())
  expect_error(
    gdpar:::.gdpar_formula_set_to_amm_spec_list(fs, W_basis_arg = list(type = "poly")),
    class = "gdpar_input_error"
  )
})

test_that("parser rejects bare symbols in the RHS (P5)", {
  fs <- gdpar_formula_set(mu = y ~ a(x) + extra_term)
  err <- expect_error(
    gdpar:::.gdpar_formula_set_to_amm_spec_list(fs, NULL),
    class = "gdpar_input_error"
  )
  expect_match(conditionMessage(err), "bare term")
})

test_that("parser rejects non-additive top-level operators", {
  fs_minus <- gdpar_formula_set(mu = y ~ a(x) - b(z))
  expect_error(
    gdpar:::.gdpar_formula_set_to_amm_spec_list(fs_minus, NULL),
    class = "gdpar_input_error"
  )
  fs_times <- gdpar_formula_set(mu = y ~ a(x) * b(z))
  expect_error(
    gdpar:::.gdpar_formula_set_to_amm_spec_list(fs_times, NULL),
    class = "gdpar_input_error"
  )
  fs_colon <- gdpar_formula_set(mu = y ~ a(x) + b(z):W())
  expect_error(
    gdpar:::.gdpar_formula_set_to_amm_spec_list(fs_colon, NULL),
    class = "gdpar_input_error"
  )
})

test_that("parser rejects a() or b() with zero or more than one argument", {
  fs_a0 <- gdpar_formula_set(mu = y ~ a() + b(z))
  expect_error(
    gdpar:::.gdpar_formula_set_to_amm_spec_list(fs_a0, NULL),
    class = "gdpar_input_error"
  )
  fs_a2 <- gdpar_formula_set(mu = y ~ a(x1, x2))
  expect_error(
    gdpar:::.gdpar_formula_set_to_amm_spec_list(fs_a2, NULL),
    class = "gdpar_input_error"
  )
})

test_that("parser propagates intercept-suppression rejection through wrappers", {
  fs <- gdpar_formula_set(mu = y ~ a(x - 1))
  expect_error(
    gdpar:::.gdpar_formula_set_to_amm_spec_list(fs, NULL),
    class = "gdpar_input_error"
  )
})

test_that("parser preserves the environment of the slot formula in sub-formulas", {
  marker <- new.env()
  assign("marker_var", 42L, envir = marker)
  f_with_env <- y ~ a(x)
  environment(f_with_env) <- marker
  fs <- gdpar_formula_set(mu = f_with_env)
  out <- gdpar:::.gdpar_formula_set_to_amm_spec_list(fs, NULL)
  expect_identical(environment(out[["mu"]]$a), marker)
})

test_that(".gdpar_promote_scope_per_observation promotes named params", {
  fam <- gdpar_family("gaussian")
  expect_equal(fam$param_specs[[1]]$scope, "per_observation")
  expect_equal(fam$param_specs[[2]]$scope, "population")
  fam2 <- gdpar:::.gdpar_promote_scope_per_observation(fam, "sigma")
  expect_equal(fam2$param_specs[[1]]$scope, "per_observation")
  expect_equal(fam2$param_specs[[2]]$scope, "per_observation")
})

test_that(".gdpar_promote_scope_per_observation is a no-op for empty k_names", {
  fam <- gdpar_family("gaussian")
  expect_identical(
    gdpar:::.gdpar_promote_scope_per_observation(fam, character(0)),
    fam
  )
  expect_identical(
    gdpar:::.gdpar_promote_scope_per_observation(fam, NULL),
    fam
  )
})

test_that(".gdpar_promote_scope_per_observation errors on non-eligible names", {
  fam <- gdpar_family("gaussian")
  err <- expect_error(
    gdpar:::.gdpar_promote_scope_per_observation(fam, "phi"),
    class = "gdpar_input_error"
  )
  expect_match(conditionMessage(err), "not eligible")
  expect_match(conditionMessage(err), "Eligible parameters")
})

test_that(".gdpar_promote_scope_per_observation does not mutate the input", {
  fam <- gdpar_family("gaussian")
  scopes_before <- vapply(fam$param_specs, function(s) s$scope, character(1L))
  invisible(
    gdpar:::.gdpar_promote_scope_per_observation(fam, "sigma")
  )
  scopes_after <- vapply(fam$param_specs, function(s) s$scope, character(1L))
  expect_identical(scopes_before, scopes_after)
})

test_that("scope promotion updates K reported by .gdpar_n_params_individual", {
  fam <- gdpar_family("gaussian")
  expect_equal(gdpar:::.gdpar_n_params_individual(fam), 1L)
  fam2 <- gdpar:::.gdpar_promote_scope_per_observation(fam, "sigma")
  expect_equal(gdpar:::.gdpar_n_params_individual(fam2), 2L)
})

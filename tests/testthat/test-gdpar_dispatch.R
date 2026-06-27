## Tests for the K-individual dispatch wired into gdpar() in Unit 2
## of sub-phase 8.3.3. Decision P-dispatch canonized in memory entry
## project_gdpar_block_8_3_extended_plan (2026-05-20): classic formula
## without AMM wrappers routes through the legacy single-amm_spec
## path bit-exactly; the three new K-individual entries
## (gdpar_formula_set, named list of amm_spec, classic formula with
## AMM wrappers in the RHS) canonize at the top of gdpar() and pass
## downstream through the same single-amm_spec body when K = 1.
## K > 1 fit-time support landed in Unit 3 via .gdpar_K() and the
## amm_distrib_K.stan template; the dispatch + fit assertions live
## in tests/testthat/test-gdpar_K.R.

.toy_data <- function() {
  data.frame(y = 1:5, x = 1:5, z = 1:5)
}

test_that("formula_set + non-default amm errors", {
  fs <- gdpar_formula_set(mu = y ~ a(x))
  err <- expect_error(
    gdpar(fs, amm = amm_spec(a = ~ x), data = .toy_data()),
    class = "gdpar_input_error"
  )
  expect_match(conditionMessage(err), "formula set")
})

test_that("classic formula with AMM wrappers + non-default amm errors", {
  err <- expect_error(
    gdpar(y ~ a(x), amm = amm_spec(a = ~ x), data = .toy_data()),
    class = "gdpar_input_error"
  )
  expect_match(conditionMessage(err), "AMM wrapper")
})

test_that("legacy classic formula + non-NULL W argument errors", {
  err <- expect_error(
    gdpar(y ~ x, amm = amm_spec(a = ~ x),
          W = W_basis(type = "polynomial", degree = 2L),
          data = .toy_data()),
    class = "gdpar_input_error"
  )
  expect_match(conditionMessage(err), "reserved for the K-individual")
})

test_that("named-list amm with non-amm_spec entry errors", {
  err <- expect_error(
    gdpar(y ~ x,
          amm  = list(mu = list(a = ~ x)),
          data = .toy_data()),
    class = "gdpar_input_error"
  )
  expect_match(conditionMessage(err), "must be an object of class amm_spec")
})

test_that("named-list amm with empty slot name errors", {
  amm_bad <- list(amm_spec(a = ~ x))
  names(amm_bad) <- ""
  expect_error(
    gdpar(y ~ x, amm = amm_bad, data = .toy_data()),
    class = "gdpar_input_error"
  )
})

test_that("named-list amm with duplicated slot names errors", {
  amm_dup <- list(mu = amm_spec(a = ~ x), mu = amm_spec(a = ~ z))
  expect_error(
    gdpar(y ~ x, amm = amm_dup, data = .toy_data()),
    class = "gdpar_input_error"
  )
})

test_that("formula_set declares W() but no external W arg supplied errors", {
  fs <- gdpar_formula_set(mu = y ~ a(x) + W())
  expect_error(
    gdpar(fs, data = .toy_data()),
    class = "gdpar_input_error"
  )
})

test_that("classic + AMM wrappers in RHS with W() and no W arg errors", {
  expect_error(
    gdpar(y ~ a(x) + W(), data = .toy_data()),
    class = "gdpar_input_error"
  )
})

test_that("classic legacy path is left untouched (no dispatch interception)", {
  err <- expect_error(
    gdpar(y_missing ~ x, amm = amm_spec(a = ~ x), data = .toy_data())
  )
  expect_s3_class(err, "gdpar_input_error")
  expect_match(conditionMessage(err), "Outcome")
})

test_that("formula_set is canonized to a classic formula reaching outcome check", {
  fs <- gdpar_formula_set(mu = y_missing ~ a(x))
  err <- expect_error(
    gdpar(fs, data = .toy_data())
  )
  expect_match(conditionMessage(err), "Outcome")
})

test_that("named-list amm with K=1 canonizes to outcome-not-in-data point", {
  err <- expect_error(
    gdpar(y_missing ~ x,
          amm  = list(mu = amm_spec(a = ~ x)),
          data = .toy_data())
  )
  expect_match(conditionMessage(err), "Outcome")
})

test_that("classic + AMM wrappers in RHS canonizes outcome to LHS symbol", {
  err <- expect_error(
    gdpar(y_missing ~ a(x), data = .toy_data())
  )
  expect_match(conditionMessage(err), "Outcome")
})

test_that("named-list amm without names falls through to legacy validator", {
  err <- expect_error(
    gdpar(y ~ x,
          amm  = list(amm_spec(a = ~ x)),
          data = .toy_data())
  )
  ## An unnamed list is NOT a K-individual named-list dispatch (the
  ## dispatch check requires names); the legacy assert_inherits(amm,
  ## "amm_spec") catches it.
  expect_match(conditionMessage(err), "amm_spec|amm")
})

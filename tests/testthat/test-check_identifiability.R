test_that("Level 0 identifiability check passes trivially", {
  df <- data.frame(x = rnorm(50))
  spec <- amm_spec()
  rep <- gdpar_check_identifiability(spec, df)
  expect_true(rep$passed)
})

test_that("Independent additive basis passes identifiability check", {
  df <- make_gaussian_data(n = 200, seed = 42)
  spec <- amm_spec(a = ~ x1 + x2)
  rep <- gdpar_check_identifiability(spec, df,
                                      formula_rhs = ~ x1 + x2)
  expect_true(rep$passed)
  expect_lt(rep$condition_number, 1e6)
})

test_that("Perfectly collinear design fails identifiability check", {
  df <- make_gaussian_data(n = 200, seed = 42)
  df$x_dup <- df$x1
  spec <- amm_spec(a = ~ x1 + x_dup)
  rep <- gdpar_check_identifiability(spec, df,
                                      formula_rhs = ~ x1 + x_dup)
  expect_false(rep$passed)
  expect_true(!is.null(rep$collinear_directions))
})

test_that("Identifiability report has all required components", {
  df <- make_gaussian_data(n = 100, seed = 42)
  spec <- amm_spec(a = ~ x1 + x2)
  rep <- gdpar_check_identifiability(spec, df,
                                      formula_rhs = ~ x1 + x2)
  expect_s3_class(rep, "gdpar_identifiability_report")
  expect_true(all(c("passed", "lambda_min", "lambda_max",
                    "condition_number", "collinear_directions",
                    "theta_ref_used", "tol_used", "column_labels",
                    "message") %in% names(rep)))
})

test_that("Identifiability check uses normalized columns", {
  set.seed(NULL)
  n <- 200
  df <- data.frame(
    x1 = rnorm(n, sd = 1),
    x2 = rnorm(n, sd = 1000),
    y = rnorm(n)
  )
  spec <- amm_spec(a = ~ x1 + x2)
  rep <- gdpar_check_identifiability(spec, df,
                                      formula_rhs = ~ x1 + x2)
  expect_true(rep$passed)
})

test_that("AMM Level 2 with collinear b and W blocks fails", {
  set.seed(NULL)
  n <- 200
  x1 <- rnorm(n)
  df <- data.frame(x1 = x1, y = rnorm(n))
  spec <- amm_spec(
    b = ~ x1,
    W = W_basis("polynomial", 1),
    x_vars = "x1"
  )
  rep <- gdpar_check_identifiability(spec, df,
                                      theta_ref_init = 1,
                                      formula_rhs = ~ x1)
  expect_false(rep$passed)
})

test_that("print method for identifiability report runs without error", {
  df <- make_gaussian_data(n = 50, seed = 42)
  spec <- amm_spec(a = ~ x1)
  rep <- gdpar_check_identifiability(spec, df,
                                      formula_rhs = ~ x1)
  expect_output(print(rep), "gdpar_identifiability_report")
})


# ---- C4-bis (multi p > 1) tests --------------------------------

test_that("C4-bis: p=1 returns c4_bis = NULL (single-coord case skips)", {
  df <- make_gaussian_data(n = 100, seed = 42)
  spec <- amm_spec(a = ~ x1 + x2)
  rep <- gdpar_check_identifiability(spec, df,
                                     formula_rhs = ~ x1 + x2)
  expect_null(rep$c4_bis)
})

test_that("C4-bis: multi spec without W still runs per-k check on Z_a[k] alone", {
  df <- make_gaussian_data(n = 100, seed = 42)
  spec <- amm_spec(p = 2L, dims = dimwise(a = ~ x1 + x2, b = NULL))
  rep <- gdpar_check_identifiability(spec, df,
                                     theta_ref_init = c(0.3, -0.2),
                                     formula_rhs = ~ x1 + x2)
  expect_false(is.null(rep$c4_bis))
  expect_identical(length(rep$c4_bis$per_k), 2L)
  for (pk in rep$c4_bis$per_k) {
    expect_true(pk$passed)
    expect_identical(pk$shared_cols, character(0))
  }
})

test_that("C4-bis full: clean multi spec passes (no overlap between Z_a[k] and X)", {
  set.seed(42)
  n <- 200
  df <- data.frame(x1 = rnorm(n), x2 = rnorm(n),
                   z1 = rnorm(n), z2 = rnorm(n))
  spec <- amm_spec(
    p = 2L,
    dims = dimwise(a = ~ x1 + x2, b = NULL),
    W = W_basis(type = "polynomial", degree = 2),
    x_vars = c("z1", "z2")
  )
  rep <- gdpar_check_identifiability(
    spec, df, theta_ref_init = c(0.4, -0.3),
    formula_rhs = ~ x1 + x2,
    rigor = "full"
  )
  expect_true(rep$passed)
  expect_false(is.null(rep$c4_bis))
  expect_identical(length(rep$c4_bis$per_k), 2L)
  for (pk in rep$c4_bis$per_k) {
    expect_true(pk$passed)
    expect_identical(pk$rigor, "full")
    expect_identical(pk$shared_cols, character(0))
  }
})

test_that("C4-bis fast: structural overlap emits warning but does not fail", {
  set.seed(42)
  n <- 200
  df <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  spec <- amm_spec(
    p = 2L,
    dims = dimwise(a = ~ x1 + x2, b = NULL),
    W = W_basis(type = "polynomial", degree = 2),
    x_vars = c("x1")
  )
  expect_warning(
    rep <- gdpar_check_identifiability(
      spec, df, theta_ref_init = c(0.4, -0.3),
      formula_rhs = ~ x1 + x2,
      rigor = "fast"
    ),
    class = "gdpar_c4bis_overlap_warning"
  )
  expect_true(rep$passed)
  for (pk in rep$c4_bis$per_k) {
    expect_identical(pk$rigor, "fast")
    expect_true(pk$passed)
    expect_true("x1" %in% pk$shared_cols)
  }
})

test_that("C4-bis full: degree-1 W with overlap and theta != 0 fails (deterministic colinearity)", {
  set.seed(42)
  n <- 300
  df <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  spec <- amm_spec(
    p = 2L,
    dims = dimwise(a = ~ x1 + x2, b = NULL),
    W = W_basis(type = "polynomial", degree = 1),
    x_vars = c("x1", "x2")
  )
  rep <- gdpar_check_identifiability(
    spec, df, theta_ref_init = c(0.5, -0.3),
    formula_rhs = ~ x1 + x2,
    rigor = "full"
  )
  # With degree 1 W, the W*x_l columns become (theta_ref[k] * x_l).
  # Combined with Z_a[k] = (x_l), these are linearly dependent
  # (scalar multiple). C4-bis should fail per-k.
  expect_false(rep$passed)
  for (pk in rep$c4_bis$per_k) {
    expect_false(pk$passed)
    expect_false(is.null(pk$collinear_directions))
  }
})

test_that("gdpar_check_identifiability rejects invalid rigor argument", {
  df <- make_gaussian_data(n = 50, seed = 1)
  spec <- amm_spec(a = ~ x1 + x2)
  expect_error(
    gdpar_check_identifiability(spec, df, rigor = "weird"),
    regexp = "should be one of"
  )
})

test_that("C4-bis report renders in print without errors", {
  set.seed(42)
  n <- 200
  df <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  spec <- amm_spec(
    p = 2L,
    dims = dimwise(a = ~ x1 + x2, b = NULL),
    W = W_basis(type = "polynomial", degree = 2),
    x_vars = c("x1")
  )
  rep <- gdpar_check_identifiability(
    spec, df, theta_ref_init = c(0.4, -0.3),
    formula_rhs = ~ x1 + x2,
    rigor = "full"
  )
  out <- capture.output(print(rep))
  expect_true(any(grepl("C4-bis", out)))
  expect_true(any(grepl("rigor", out)))
  expect_true(any(grepl("coord", out)))
})

test_that("D-ID pre-fit slot is NULL by default (backward compatible)", {
  set.seed(7)
  df <- data.frame(x1 = rnorm(40))
  spec <- amm_spec(a = ~ x1, b = NULL, W = NULL)
  rep <- gdpar_check_identifiability(spec, df)
  expect_null(rep$did_pre_fit)
})

test_that("D-ID pre-fit slot populated when family is supplied (K=1 default)", {
  set.seed(8)
  df <- data.frame(x1 = rnorm(40))
  spec <- amm_spec(a = ~ x1, b = NULL, W = NULL)
  fam <- gdpar_family("gaussian")
  rep <- gdpar_check_identifiability(spec, df, family = fam)
  expect_true(!is.null(rep$did_pre_fit))
  expect_equal(rep$did_pre_fit$K, 1L)
  expect_equal(rep$did_pre_fit$per_param[[1L]]$name, "mu")
  expect_equal(rep$did_pre_fit$per_param[[1L]]$scope, "per_observation")
  expect_true(rep$did_pre_fit$passed)
  expect_null(rep$did_pre_fit$symbolic_separability)
})

test_that("D-ID pre-fit symbolic separability is checked when K >= 2 with rigor='full'", {
  fake_K2 <- gdpar_family("gaussian")
  fake_K2$param_specs[[2L]]$scope <- "per_group"
  set.seed(9)
  df <- data.frame(x1 = rnorm(40))
  spec <- amm_spec(a = ~ x1, b = NULL, W = NULL)
  rep <- gdpar_check_identifiability(spec, df, family = fake_K2,
                                     rigor = "full")
  expect_equal(rep$did_pre_fit$K, 2L)
  expect_true(!is.null(rep$did_pre_fit$symbolic_separability))
  expect_true(rep$did_pre_fit$symbolic_separability$passed)
})

test_that("D-ID pre-fit symbolic separability flags duplicate prior_canonical_kind", {
  fake_collision <- gdpar_family("gaussian")
  fake_collision$param_specs[[2L]]$scope <- "per_group"
  fake_collision$param_specs[[2L]]$prior_canonical_kind <- "mu"
  set.seed(10)
  df <- data.frame(x1 = rnorm(40))
  spec <- amm_spec(a = ~ x1, b = NULL, W = NULL)
  rep <- gdpar_check_identifiability(spec, df, family = fake_collision,
                                     rigor = "full")
  expect_false(rep$did_pre_fit$symbolic_separability$passed)
  expect_false(rep$did_pre_fit$passed)
  expect_false(rep$passed)
})

test_that("D-ID pre-fit per-param echoes did_status from family's param_specs", {
  set.seed(11)
  df <- data.frame(x1 = rnorm(40))
  spec <- amm_spec(a = ~ x1, b = NULL, W = NULL)
  fam <- gdpar_family("gaussian")
  rep <- gdpar_check_identifiability(spec, df, family = fam)
  expect_equal(rep$did_pre_fit$per_param[[1L]]$did_status, "holds")
  expect_equal(rep$did_pre_fit$per_param[[1L]]$prior_canonical_kind, "mu")
})

test_that("print() of identifiability report shows D-ID block when populated", {
  set.seed(12)
  df <- data.frame(x1 = rnorm(40))
  spec <- amm_spec(a = ~ x1, b = NULL, W = NULL)
  fam <- gdpar_family("gaussian")
  rep <- gdpar_check_identifiability(spec, df, family = fam)
  out <- capture.output(print(rep))
  expect_true(any(grepl("D-ID pre-fit", out)))
  expect_true(any(grepl("K \\(individual params\\)", out)))
  expect_true(any(grepl("param 'mu'", out)))
})

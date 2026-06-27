toggle_combinations <- list(
  list(a = "ncp", W = "ncp"),
  list(a = "ncp", W = "cp"),
  list(a = "cp",  W = "ncp"),
  list(a = "cp",  W = "cp")
)

fit_with_toggles <- function(family_name, df, tg) {
  spec <- amm_spec(
    a = ~ x1 + x2,
    W = W_basis(type = "polynomial", degree = 2),
    x_vars = c("x1", "x2")
  )
  gdpar(
    formula           = y ~ x1 + x2,
    family            = gdpar_family(family_name),
    amm               = spec,
    data              = df,
    iter_warmup       = 200L,
    iter_sampling     = 200L,
    chains            = 2L,
    parametrization_a = tg$a,
    parametrization_W = tg$W,
    skip_id_check     = TRUE,
    verbose           = FALSE,
    refresh           = 0L,
    seed              = 42L
  )
}

test_that("gdpar fits all 4 CP/NCP toggle combinations for gaussian family", {
  skip_if_no_cmdstan()
  df <- make_gaussian_data(n = 100, seed = 42)
  for (tg in toggle_combinations) {
    fit <- fit_with_toggles("gaussian", df, tg)
    expect_true(inherits(fit, "gdpar_fit"))
    expect_equal(fit$parametrization$cp_a, tg$a == "cp")
    expect_equal(fit$parametrization$cp_W, tg$W == "cp")
    expect_false(fit$parametrization$meta$used_preflight)
  }
})

test_that("gdpar fits all 4 CP/NCP toggle combinations for poisson family", {
  skip_if_no_cmdstan()
  df <- make_poisson_data(n = 100, seed = 42)
  for (tg in toggle_combinations) {
    fit <- fit_with_toggles("poisson", df, tg)
    expect_true(inherits(fit, "gdpar_fit"))
    expect_equal(fit$parametrization$cp_a, tg$a == "cp")
    expect_equal(fit$parametrization$cp_W, tg$W == "cp")
    expect_false(fit$parametrization$meta$used_preflight)
  }
})

test_that("gdpar fits all 4 CP/NCP toggle combinations for bernoulli family", {
  skip_if_no_cmdstan()
  df <- make_bernoulli_data(n = 200, seed = 42)
  for (tg in toggle_combinations) {
    fit <- fit_with_toggles("bernoulli", df, tg)
    expect_true(inherits(fit, "gdpar_fit"))
    expect_equal(fit$parametrization$cp_a, tg$a == "cp")
    expect_equal(fit$parametrization$cp_W, tg$W == "cp")
    expect_false(fit$parametrization$meta$used_preflight)
  }
})

test_that("gdpar fits all 4 CP/NCP toggle combinations for neg_binomial_2 family", {
  skip_if_no_cmdstan()
  df <- make_neg_binomial_data(n = 100, seed = 42)
  for (tg in toggle_combinations) {
    fit <- fit_with_toggles("neg_binomial_2", df, tg)
    expect_true(inherits(fit, "gdpar_fit"))
    expect_equal(fit$parametrization$cp_a, tg$a == "cp")
    expect_equal(fit$parametrization$cp_W, tg$W == "cp")
    expect_false(fit$parametrization$meta$used_preflight)
  }
})

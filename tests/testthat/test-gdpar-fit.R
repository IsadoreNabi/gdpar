test_that("gdpar rejects unsupported paths", {
  df <- make_gaussian_data(n = 50, seed = 42)
  expect_error(
    gdpar(y ~ x1, data = df, path = "vcm"),
    class = "gdpar_unsupported_feature_error"
  )
  expect_error(
    gdpar(y ~ x1, data = df, path = "hyper"),
    class = "gdpar_unsupported_feature_error"
  )
})

test_that("gdpar rejects malformed inputs", {
  df <- make_gaussian_data(n = 50, seed = 42)
  expect_error(
    gdpar("not a formula", data = df),
    class = "gdpar_input_error"
  )
  expect_error(
    gdpar(y ~ x1, data = "not a data frame"),
    class = "gdpar_input_error"
  )
  expect_error(
    gdpar(y ~ x1, data = df, chains = 0),
    class = "gdpar_input_error"
  )
  expect_error(
    gdpar(y ~ x1, data = df, adapt_delta = 1.5),
    class = "gdpar_input_error"
  )
})

test_that("gdpar fits a Level 0 Gaussian model end-to-end", {
  skip_on_cran()
  skip_if_no_cmdstan()
  df <- make_gaussian_data(n = 100, seed = 42)
  fit <- gdpar(
    formula       = y ~ x1,
    family        = gdpar_family("gaussian"),
    amm           = amm_spec(),
    data          = df,
    iter_warmup   = 200,
    iter_sampling = 200,
    chains        = 2,
    refresh       = 0,
    verbose       = FALSE
  )
  expect_s3_class(fit, "gdpar_fit")
  expect_s3_class(fit$diagnostics, "gdpar_diagnostics")
})

test_that("gdpar fits an AMM Level 1 Gaussian model end-to-end", {
  skip_on_cran()
  skip_if_no_cmdstan()
  df <- make_gaussian_data(n = 200, seed = 42)
  fit <- gdpar(
    formula       = y ~ x1 + x2,
    family        = gdpar_family("gaussian"),
    amm           = amm_spec(a = ~ x1 + x2),
    data          = df,
    iter_warmup   = 300,
    iter_sampling = 300,
    chains        = 2,
    refresh       = 0,
    verbose       = FALSE
  )
  expect_s3_class(fit, "gdpar_fit")
  cf <- coef(fit)
  expect_s3_class(cf, "gdpar_coef")
  expect_equal(cf$p, 1L)
  expect_true(!is.null(cf$a))
  expect_length(cf$a, 1L)
  expect_s3_class(cf$a[[1L]], "data.frame")
  expect_equal(cf$a[[1L]]$term, c("x1", "x2"))
  expect_named(cf$a[[1L]],
               c("term", "mean", "q05", "q50", "q95"))
})

test_that("predict on training data returns a draws matrix", {
  skip_on_cran()
  skip_if_no_cmdstan()
  df <- make_gaussian_data(n = 100, seed = 42)
  fit <- gdpar(
    formula       = y ~ x1,
    amm           = amm_spec(a = ~ x1),
    data          = df,
    iter_warmup   = 200,
    iter_sampling = 200,
    chains        = 2,
    refresh       = 0,
    verbose       = FALSE
  )
  pr <- predict(fit, summary = "draws")
  expect_true(is.matrix(pr))
  expect_equal(ncol(pr), nrow(df))
})

test_that("predict with mean_se returns one row per observation", {
  skip_on_cran()
  skip_if_no_cmdstan()
  df <- make_gaussian_data(n = 100, seed = 42)
  fit <- gdpar(
    formula       = y ~ x1,
    amm           = amm_spec(a = ~ x1),
    data          = df,
    iter_warmup   = 200,
    iter_sampling = 200,
    chains        = 2,
    refresh       = 0,
    verbose       = FALSE
  )
  pr <- predict(fit, summary = "mean_se")
  expect_s3_class(pr, "data.frame")
  expect_equal(nrow(pr), nrow(df))
  expect_named(pr, c("mean", "se"))
})

test_that("predict on newdata produces n_new columns", {
  skip_on_cran()
  skip_if_no_cmdstan()
  df <- make_gaussian_data(n = 100, seed = 42)
  newdf <- make_gaussian_data(n = 30, seed = 7)
  fit <- gdpar(
    formula       = y ~ x1,
    amm           = amm_spec(a = ~ x1),
    data          = df,
    iter_warmup   = 200,
    iter_sampling = 200,
    chains        = 2,
    refresh       = 0,
    verbose       = FALSE
  )
  pr <- predict(fit, newdata = newdf, summary = "mean_se")
  expect_equal(nrow(pr), nrow(newdf))
})

test_that("predict on newdata with J_a >= 2 returns full draws matrix", {
  skip_on_cran()
  skip_if_no_cmdstan()
  df <- make_gaussian_data(n = 100, seed = 42)
  newdf <- make_gaussian_data(n = 30, seed = 7)
  fit <- gdpar(
    formula       = y ~ x1 + x2,
    amm           = amm_spec(a = ~ x1 + x2),
    data          = df,
    iter_warmup   = 200,
    iter_sampling = 200,
    chains        = 2,
    refresh       = 0,
    verbose       = FALSE
  )
  pr <- predict(fit, newdata = newdf, summary = "draws")
  expect_true(is.matrix(pr))
  expect_equal(ncol(pr), nrow(newdf))
  expect_equal(nrow(pr), 400L)
  expect_true(all(is.finite(pr)))
})

test_that("identifiability failure aborts the fit", {
  skip_on_cran()
  skip_if_no_cmdstan()
  df <- make_gaussian_data(n = 100, seed = 42)
  df$x_dup <- df$x1
  expect_error(
    gdpar(
      formula = y ~ x1 + x_dup,
      amm     = amm_spec(a = ~ x1 + x_dup),
      data    = df,
      iter_warmup = 200, iter_sampling = 200, chains = 2,
      refresh = 0, verbose = FALSE
    ),
    class = "gdpar_identifiability_error"
  )
})

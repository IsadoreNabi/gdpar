.make_multi_dataset <- function(n = 20L, p = 2L, seed = 1L) {
  set.seed(seed)
  df <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  df$y <- matrix(rnorm(n * p), nrow = n, ncol = p)
  df
}

test_that("gdpar() dispatches to .gdpar_multi when amm$p > 1L (not the old guard)", {
  skip_if_not_installed("cmdstanr")
  df <- .make_multi_dataset()
  spec <- amm_spec(p = 2L, dims = dimwise(a = ~ x1 + x2, b = NULL))
  df$y <- rnorm(nrow(df))
  expect_error(
    gdpar(
      formula = y ~ x1 + x2,
      family  = gdpar_family_multi("gaussian", p = 2L),
      amm     = spec,
      data    = df,
      verbose = FALSE
    ),
    class = "gdpar_input_error",
    regexp = "matrix of dimensions n by p"
  )
})

test_that("gdpar() multi rejects outcome with wrong number of columns", {
  skip_if_not_installed("cmdstanr")
  df <- .make_multi_dataset(p = 2L)
  df$y <- matrix(rnorm(nrow(df) * 3L), nrow = nrow(df), ncol = 3L)
  spec <- amm_spec(p = 2L, dims = dimwise(a = ~ x1 + x2, b = NULL))
  expect_error(
    gdpar(
      formula = y ~ x1 + x2,
      family  = gdpar_family_multi("gaussian", p = 2L),
      amm     = spec,
      data    = df,
      verbose = FALSE
    ),
    class = "gdpar_input_error",
    regexp = "3 columns but amm\\$p = 2"
  )
})

test_that("gdpar() multi rejects outcome with NAs", {
  skip_if_not_installed("cmdstanr")
  df <- .make_multi_dataset(p = 2L)
  df$y[1L, 1L] <- NA_real_
  spec <- amm_spec(p = 2L, dims = dimwise(a = ~ x1 + x2, b = NULL))
  expect_error(
    gdpar(
      formula = y ~ x1 + x2,
      family  = gdpar_family_multi("gaussian", p = 2L),
      amm     = spec,
      data    = df,
      verbose = FALSE
    ),
    class = "gdpar_input_error",
    regexp = "contains 1 non-finite"
  )
})

test_that("gdpar() multi auto-promotes a univariate gdpar_family to gdpar_family_multi", {
  skip_if_not_installed("cmdstanr")
  df <- .make_multi_dataset(p = 2L)
  spec <- amm_spec(p = 2L, dims = dimwise(a = ~ x1 + x2, b = NULL))
  # H.2 contract: auto-promote emits the message and the fit
  # completes (cp_a multi is now wired). We assert only on the
  # promotion message; the diagnostic warnings from a short fit are
  # suppressed so they do not bury the test assertion.
  expect_message(
    suppressWarnings(
      gdpar(
        formula           = y ~ x1 + x2,
        family            = gdpar_family("gaussian"),
        amm               = spec,
        data              = df,
        parametrization   = "ncp",
        chains            = 2L,
        iter_warmup       = 50L,
        iter_sampling     = 50L,
        refresh           = 0L,
        seed              = 42L,
        verbose           = TRUE
      )
    ),
    regexp = "Auto-promoted univariate family 'gaussian'.*p = 2"
  )
})

test_that("generate_stan_code_multi accepts cp_a = TRUE (Phase H.2 wiring)", {
  src_ncp <- gdpar:::generate_stan_code_multi(
    gdpar_prior(), cp_a = FALSE, cp_W = FALSE
  )
  src_cp <- gdpar:::generate_stan_code_multi(
    gdpar_prior(), cp_a = TRUE, cp_W = FALSE
  )
  expect_type(src_cp, "character")
  expect_match(src_cp, "segment\\(a_raw")
  expect_false(grepl("segment\\(a_raw", src_ncp))
  # NCP-only via legacy scalar arg is byte-identical to NCP via
  # cp_a_per_k explicit (bit-exact paridad with H.1).
  src_ncp_explicit <- gdpar:::generate_stan_code_multi(
    gdpar_prior(), cp_a_per_k = c(FALSE, FALSE)
  )
  expect_identical(src_ncp, src_ncp_explicit)
})

test_that("gdpar() multi rejects family with mismatched p", {
  skip_if_not_installed("cmdstanr")
  df <- .make_multi_dataset(p = 2L)
  spec <- amm_spec(p = 2L, dims = dimwise(a = ~ x1 + x2, b = NULL))
  expect_error(
    gdpar(
      formula = y ~ x1 + x2,
      family  = gdpar_family_multi("gaussian", p = 3L),
      amm     = spec,
      data    = df,
      verbose = FALSE
    ),
    class = "gdpar_input_error",
    regexp = "Family p = 3 does not match amm\\$p = 2"
  )
})

test_that("gdpar() scalar path: amm_spec default assigns amm$p = 1L (no dispatch to multi)", {
  spec <- amm_spec(a = ~ x1)
  expect_equal(spec$p, 1L)
})

test_that("gdpar() scalar path: amm_spec(p = 1L, a = ...) keeps amm$p = 1L (no dispatch to multi)", {
  spec <- amm_spec(p = 1L, a = ~ x1, b = NULL)
  expect_equal(spec$p, 1L)
})

test_that("resolve_anchor_multi: scalar broadcasts to length p", {
  fam <- gdpar_family_multi("gaussian", p = 3L)
  y <- matrix(rnorm(30L), nrow = 10L, ncol = 3L)
  out <- resolve_anchor_multi(0.7, fam, y, gdpar_prior(), p = 3L,
                              verbose = FALSE)
  expect_equal(out, rep(0.7, 3L))
})

test_that("resolve_anchor_multi: vector of length p is preserved", {
  fam <- gdpar_family_multi("gaussian", p = 2L)
  y <- matrix(rnorm(20L), nrow = 10L, ncol = 2L)
  out <- resolve_anchor_multi(c(1.5, -2.5), fam, y, gdpar_prior(), p = 2L,
                              verbose = FALSE)
  expect_equal(out, c(1.5, -2.5))
})

test_that("resolve_anchor_multi: 'prior_mean' returns zeros", {
  fam <- gdpar_family_multi("gaussian", p = 4L)
  y <- matrix(rnorm(40L), nrow = 10L, ncol = 4L)
  out <- resolve_anchor_multi("prior_mean", fam, y, gdpar_prior(), p = 4L,
                              verbose = FALSE)
  expect_equal(out, rep(0.0, 4L))
})

test_that("resolve_anchor_multi: 'empirical_y' applies link per coordinate", {
  fam <- gdpar_family_multi("gaussian", p = 2L)
  y <- cbind(rep(3.0, 10L), rep(-1.0, 10L))
  out <- resolve_anchor_multi("empirical_y", fam, y, gdpar_prior(), p = 2L,
                              verbose = FALSE)
  expect_equal(out, c(3.0, -1.0))
})

test_that("resolve_anchor_multi: wrong vector length errors", {
  fam <- gdpar_family_multi("gaussian", p = 2L)
  y <- matrix(rnorm(20L), nrow = 10L, ncol = 2L)
  expect_error(
    resolve_anchor_multi(c(1, 2, 3), fam, y, gdpar_prior(), p = 2L,
                         verbose = FALSE),
    class = "gdpar_input_error",
    regexp = "scalar.*vector of length p"
  )
})

test_that("resolve_anchor_multi: invalid character anchor errors", {
  fam <- gdpar_family_multi("gaussian", p = 2L)
  y <- matrix(rnorm(20L), nrow = 10L, ncol = 2L)
  expect_error(
    resolve_anchor_multi("nonsense", fam, y, gdpar_prior(), p = 2L,
                         verbose = FALSE),
    class = "gdpar_input_error"
  )
})

test_that("print.gdpar_fit() multi shows p (theta_ref dim) and observations", {
  obj <- structure(
    list(
      p = 2L,
      path = "bayes",
      family = list(name = "gaussian", link = "identity"),
      amm = structure(list(level = 1, p = 2L),
                       class = c("amm_spec", "list")),
      anchor = c(0.5, -0.3),
      design = list(
        Z_a_list = list(matrix(0, 15L, 2L), matrix(0, 15L, 2L)),
        Z_b_list = list(matrix(0, 15L, 0L), matrix(0, 15L, 0L)),
        X = matrix(0, 15L, 0L)
      ),
      diagnostics = list(
        converged       = TRUE,
        rhat_max        = 1.002,
        ess_bulk_min    = 500,
        divergent_count = 0L
      )
    ),
    class = c("gdpar_fit", "list")
  )
  out <- capture.output(print(obj))
  expect_true(any(grepl("p \\(theta_ref dim\\)\\s*:\\s*2", out)))
  expect_true(any(grepl("observations\\s*:\\s*15", out)))
  anchor_lines <- out[grepl("^\\s*anchor\\s*:", out)]
  expect_length(anchor_lines, 1L)
  expect_match(anchor_lines, "0\\.5")
  expect_match(anchor_lines, "-0\\.3")
  expect_match(anchor_lines, "\\[")
  expect_match(anchor_lines, "\\]")
})

test_that("coef.gdpar_fit() multi no longer aborts with deferred-feature error (Phase H.4 unblock)", {
  # Phase H.4 wired coef() multi to a real gdpar_coef builder. The
  # deferred-feature abort path was removed. Without a real fit object
  # (no $fit slot with $draws()), the function still fails downstream,
  # but it must not produce the gdpar_unsupported_feature_error.
  obj <- structure(list(p = 2L), class = c("gdpar_fit", "list"))
  err <- tryCatch(coef(obj), error = function(e) e)
  expect_s3_class(err, "error")
  expect_false(inherits(err, "gdpar_unsupported_feature_error"))
})

test_that("coef.gdpar_fit() scalar path also does not abort with the multi deferred-feature error", {
  obj <- structure(list(p = NULL), class = c("gdpar_fit", "list"))
  err <- tryCatch(coef(obj), error = function(e) e)
  expect_s3_class(err, "error")
  expect_false(inherits(err, "gdpar_unsupported_feature_error"))
})

test_that("gdpar() rejects gdpar_family_multi when amm$p is NULL (scalar amm)", {
  spec <- amm_spec(a = ~ x1)
  df <- data.frame(x1 = rnorm(10L), y = rnorm(10L))
  expect_error(
    gdpar(
      formula = y ~ x1,
      family  = gdpar_family_multi("gaussian", p = 2L),
      amm     = spec,
      data    = df,
      verbose = FALSE
    ),
    class = "gdpar_input_error",
    regexp = "multivariate family.*amm\\$p is NULL or 1"
  )
})

test_that("gdpar() rejects gdpar_family_multi when amm$p = 1L", {
  spec <- amm_spec(p = 1L, a = ~ x1, b = NULL)
  df <- data.frame(x1 = rnorm(10L), y = rnorm(10L))
  expect_error(
    gdpar(
      formula = y ~ x1,
      family  = gdpar_family_multi("gaussian", p = 1L),
      amm     = spec,
      data    = df,
      verbose = FALSE
    ),
    class = "gdpar_input_error",
    regexp = "multivariate family.*amm\\$p is NULL or 1"
  )
})

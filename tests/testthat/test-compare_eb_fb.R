# Sub-phase 8.6.E (Session 14 closure, 2026-05-26): unit tests for
# gdpar_compare_eb_fb() and its S3 methods. Mirror of
# test-compare_meta_learners.R (8.5.B) for the EB-vs-FB comparator.
#
# Section 1: input validation guards (do not require cmdstanr).
# Section 2: structural assertions on the returned object using a
#   synthetic Gaussian K=1+p=1 pair (cmdstanr-gated, opt-in via
#   GDPAR_RUN_STAN_SMOKE_EB).
# Section 3: internal helpers
#   .gdpar_eb_fb_tv_table / .gdpar_eb_fb_coverage_table / etc.

# ---- Section 1: input validation guards ---------------------------------

test_that("gdpar_compare_eb_fb rejects non-gdpar_eb_fit eb_fit", {
  expect_error(
    gdpar_compare_eb_fb(eb_fit = list(), fb_fit = list()),
    class = "gdpar_input_error"
  )
})

test_that("gdpar_compare_eb_fb rejects non-gdpar_fit fb_fit", {
  fake_eb <- structure(list(), class = c("gdpar_eb_fit", "list"))
  expect_error(
    gdpar_compare_eb_fb(eb_fit = fake_eb, fb_fit = list()),
    class = "gdpar_input_error"
  )
})

test_that("gdpar_compare_eb_fb rejects level outside (0, 1)", {
  fake_eb <- structure(list(), class = c("gdpar_eb_fit", "list"))
  fake_fb <- structure(list(), class = c("gdpar_fit", "list"))
  expect_error(
    gdpar_compare_eb_fb(fake_eb, fake_fb, level = 0),
    class = "gdpar_input_error"
  )
  expect_error(
    gdpar_compare_eb_fb(fake_eb, fake_fb, level = 1),
    class = "gdpar_input_error"
  )
  expect_error(
    gdpar_compare_eb_fb(fake_eb, fake_fb, level = c(0.5, 0.95)),
    class = "gdpar_input_error"
  )
})

test_that("gdpar_compare_eb_fb rejects tv_bins below the minimum 5", {
  fake_eb <- structure(list(), class = c("gdpar_eb_fit", "list"))
  fake_fb <- structure(list(), class = c("gdpar_fit", "list"))
  expect_error(
    gdpar_compare_eb_fb(fake_eb, fake_fb, tv_bins = 4L),
    class = "gdpar_input_error"
  )
})

# ---- Section 2: structural assertions on a synthetic K=1+p=1 fit -------

.compare_eb_fb_smoke_env_active <- function() {
  nzchar(Sys.getenv("GDPAR_RUN_STAN_SMOKE_EB"))
}

test_that("gdpar_compare_eb_fb on Gaussian K=1+p=1 returns the expected shape (opt-in)", {
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("posterior")
  if (!.compare_eb_fb_smoke_env_active()) {
    skip("Set GDPAR_RUN_STAN_SMOKE_EB to a non-empty value to enable.")
  }
  set.seed(20260526L)
  n <- 120L
  df <- data.frame(x = stats::rnorm(n))
  df$y <- 0.5 + 0.4 * df$x + stats::rnorm(n, sd = 0.3)
  spec <- amm_spec(a = ~ x)
  fit_eb <- gdpar_eb(
    formula = y ~ x, family = gdpar_family("gaussian"),
    amm = spec, data = df,
    iter_warmup = 200L, iter_sampling = 200L, chains = 2L,
    refresh = 0L, verbose = FALSE, seed = 100L
  )
  fit_fb <- gdpar(
    formula = y ~ x, family = gdpar_family("gaussian"),
    amm = spec, data = df,
    iter_warmup = 200L, iter_sampling = 200L, chains = 2L,
    refresh = 0L, verbose = FALSE, seed = 100L
  )
  cmp <- gdpar_compare_eb_fb(fit_eb, fit_fb, level = 0.95,
                              tv_bins = 25L)
  expect_s3_class(cmp, "gdpar_eb_fb_comparison")
  expect_equal(cmp$level, 0.95)
  expect_equal(cmp$tv_bins, 25L)
  expect_true(is.data.frame(cmp$theta_diff_table))
  expect_true("eb_estimate" %in% colnames(cmp$theta_diff_table))
  expect_true("fb_mean" %in% colnames(cmp$theta_diff_table))
  expect_true("diff" %in% colnames(cmp$theta_diff_table))
  expect_true(is.data.frame(cmp$tv_table))
  expect_true("tv" %in% colnames(cmp$tv_table))
  expect_true(all(is.finite(cmp$tv_table$tv) | is.na(cmp$tv_table$tv)))
  if (any(is.finite(cmp$tv_table$tv))) {
    tv_vals <- cmp$tv_table$tv[is.finite(cmp$tv_table$tv)]
    expect_true(all(tv_vals >= 0 - 1e-8 & tv_vals <= 1 + 1e-8))
  }
  expect_true(is.data.frame(cmp$coverage_table))
  expect_true("width_ratio" %in% colnames(cmp$coverage_table))
  s_cmp <- summary(cmp)
  expect_s3_class(s_cmp, "summary.gdpar_eb_fb_comparison")
  expect_invisible(print(cmp))
  expect_invisible(print(s_cmp))
})

# ---- Section 3: internal helpers ----------------------------------------

test_that(".gdpar_eb_fb_tv_table returns NULL when one input is NULL", {
  m <- matrix(stats::rnorm(50L), 25L, 2L,
              dimnames = list(NULL, c("a", "b")))
  expect_null(.gdpar_eb_fb_tv_table(NULL, m, tv_bins = 10L))
  expect_null(.gdpar_eb_fb_tv_table(m, NULL, tv_bins = 10L))
})

test_that(".gdpar_eb_fb_tv_table on identical samples returns near-zero TV", {
  set.seed(1L)
  m <- matrix(stats::rnorm(2000L), 500L, 4L,
              dimnames = list(NULL, c("a", "b", "c", "d")))
  res <- .gdpar_eb_fb_tv_table(m, m, tv_bins = 30L)
  expect_true(is.data.frame(res))
  expect_equal(nrow(res), 4L)
  expect_true(all(res$tv == 0 | abs(res$tv) < 1e-8))
})

test_that(".gdpar_eb_fb_tv_table on disjoint samples returns TV near 1", {
  set.seed(2L)
  m1 <- matrix(stats::rnorm(500L, mean = -10), 500L, 1L,
               dimnames = list(NULL, "a"))
  m2 <- matrix(stats::rnorm(500L, mean =  10), 500L, 1L,
               dimnames = list(NULL, "a"))
  res <- .gdpar_eb_fb_tv_table(m1, m2, tv_bins = 30L)
  expect_true(is.data.frame(res))
  expect_equal(nrow(res), 1L)
  expect_true(isTRUE(res$tv > 0.9))
})

test_that(".gdpar_eb_fb_extract_xi_draws returns NULL for NULL fit", {
  expect_null(.gdpar_eb_fb_extract_xi_draws(NULL))
})

# ---- Section 4: D51 (Sesion 18 Etapa 2) print method surfaces $warnings -

test_that("print.gdpar_eb_fb_comparison shows $warnings when populated", {
  mock_cmp <- structure(
    list(
      theta_diff_table = data.frame(
        cell = 1L, eb_estimate = 0.5, eb_se = 0.1,
        fb_mean = NA_real_, fb_se = NA_real_,
        diff = NA_real_, diff_rel = NA_real_,
        stringsAsFactors = FALSE
      ),
      tv_table = NULL,
      coverage_table = data.frame(
        cell = 1L, eb_width = 0.4, fb_width = NA_real_,
        width_ratio = NA_real_, inflation = 1,
        stringsAsFactors = FALSE
      ),
      level = 0.95,
      tv_bins = 30L,
      n_common_params = 0L,
      path_eb = "eb",
      path_fb = "fb",
      family_eb = "gaussian",
      family_fb = "gaussian",
      call = quote(gdpar_compare_eb_fb()),
      warnings = c(
        "Empirical-Bayes vs Fully-Bayes anchor diff table: every FB cell is NA.",
        "Fully-Bayes xi draws unavailable from fb_fit$conditional_fit / fb_fit$fit."
      ),
      meta = list(mode = "compare_eb_fb", note = "n/a")
    ),
    class = c("gdpar_eb_fb_comparison", "list")
  )
  expect_output(print(mock_cmp), "Warnings:", fixed = TRUE)
  expect_output(print(mock_cmp), "anchor diff table", fixed = TRUE)
  expect_invisible(print(mock_cmp))
})

test_that("print.gdpar_eb_fb_comparison omits Warnings: block when empty", {
  mock_cmp <- structure(
    list(
      theta_diff_table = data.frame(
        cell = 1L, eb_estimate = 0.5, eb_se = 0.1,
        fb_mean = 0.45, fb_se = 0.12,
        diff = 0.05, diff_rel = 0.4,
        stringsAsFactors = FALSE
      ),
      tv_table = data.frame(parameter = "a", tv = 0.1, n_eb = 100L,
                              n_fb = 100L, mean_eb = 0.5, mean_fb = 0.5,
                              stringsAsFactors = FALSE),
      coverage_table = data.frame(cell = 1L, eb_width = 0.4,
                                    fb_width = 0.42, width_ratio = 0.95,
                                    inflation = 1,
                                    stringsAsFactors = FALSE),
      level = 0.95,
      tv_bins = 30L,
      n_common_params = 1L,
      path_eb = "eb",
      path_fb = "fb",
      family_eb = "gaussian",
      family_fb = "gaussian",
      call = quote(gdpar_compare_eb_fb()),
      warnings = character(0L),
      meta = list(mode = "compare_eb_fb", note = "n/a")
    ),
    class = c("gdpar_eb_fb_comparison", "list")
  )
  out <- capture.output(print(mock_cmp))
  expect_false(any(grepl("^Warnings:", out)))
})

test_that("summary.gdpar_eb_fb_comparison propagates $warnings", {
  mock_cmp <- structure(
    list(
      theta_diff_table = data.frame(
        cell = 1L, eb_estimate = 0.5, eb_se = 0.1,
        fb_mean = NA_real_, fb_se = NA_real_,
        diff = NA_real_, diff_rel = NA_real_,
        stringsAsFactors = FALSE
      ),
      tv_table = NULL,
      coverage_table = NULL,
      level = 0.95,
      tv_bins = 30L,
      n_common_params = 0L,
      path_eb = "eb",
      path_fb = "fb",
      family_eb = "gaussian",
      family_fb = "gaussian",
      call = quote(gdpar_compare_eb_fb()),
      warnings = "Empirical-Bayes vs Fully-Bayes anchor diff table: every FB cell is NA.",
      meta = list(mode = "compare_eb_fb", note = "n/a")
    ),
    class = c("gdpar_eb_fb_comparison", "list")
  )
  s <- summary(mock_cmp)
  expect_true("warnings" %in% names(s))
  expect_identical(s$warnings, mock_cmp$warnings)
  expect_output(print(s), "Warnings:", fixed = TRUE)
})

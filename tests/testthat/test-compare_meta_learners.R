# Tests for gdpar_compare_meta_learners (Sub-phase 8.5.B; external
# meta-learner comparator with pluggable adapter contract).
#
# Layered coverage:
#   (i)    constructor + predicate + S3 print of gdpar_meta_learner_adapter
#          (unit; no external dependencies).
#   (ii)   unit tests of the internal helpers .compute_comparison_metrics,
#          .extract_covariates, .guard_scalar_outcome, .validate_adapter_output,
#          .check_adapter_requirements (unit; no external dependencies).
#   (iii)  orchestrator tests with mock bridge + mock adapter (no Stan,
#          no grf, no econml).
#   (iv)   smoke + algorithmic-equivalence tests of gdpar_adapter_grf gated
#          by requireNamespace("grf").
#   (v)    smoke test of gdpar_adapter_econml gated by env-var
#          GDPAR_RUN_PYTHON_SMOKE_ECONML and reticulate::py_module_available.
#   (vi)   S3 method tests (print, summary, predict).

# ---- helpers ------------------------------------------------------------

make_mock_compare_bridge <- function(newdata, cate_mean,
                                      cate_ci = NULL,
                                      level = 0.95,
                                      dim_kind = "scalar",
                                      data_t = NULL,
                                      data_c = NULL,
                                      outcome_name = "y",
                                      formula = y ~ x1) {
  if (is.null(cate_ci)) {
    cate_ci <- cbind(lower = cate_mean - 0.1,
                     upper = cate_mean + 0.1)
  }
  fit_skeleton <- function(d) {
    cl <- as.call(list(
      as.name("gdpar"),
      formula = formula,
      data    = d
    ))
    obj <- list(
      family = list(name = "gaussian", link = "identity"),
      stan_data = list(use_groups = 0L),
      anchor = 0,
      call = cl
    )
    class(obj) <- c("gdpar_fit", "list")
    obj
  }
  bridge <- list(
    cate_draws = NULL,
    cate_mean  = cate_mean,
    cate_ci    = cate_ci,
    newdata    = newdata,
    fits       = list(treat = fit_skeleton(data_t),
                      ctrl  = fit_skeleton(data_c)),
    type       = "response",
    level      = level,
    n_draws    = 200L,
    n_obs      = nrow(newdata),
    call       = call("gdpar_causal_bridge"),
    meta       = list(dim_kind = dim_kind,
                      dim_size = 1L,
                      dim_names = NULL,
                      newdata_source = "user")
  )
  class(bridge) <- c("gdpar_causal_bridge", "list")
  bridge
}

make_mock_adapter <- function(name = "mock",
                              cate_value = 0.5,
                              has_predict = TRUE,
                              native_ci = TRUE) {
  fp <- function(X, Y, T, X_newdata, level, seed_run) {
    n <- nrow(X_newdata)
    m <- rep(cate_value, n)
    list(cate_mean = m,
         cate_ci   = if (native_ci) cbind(lower = m - 0.05,
                                           upper = m + 0.05) else NULL,
         state     = list(value = cate_value),
         notes     = character(0L))
  }
  pf <- if (has_predict) {
    function(state, X_newdata, level) {
      n <- nrow(X_newdata)
      list(cate_mean = rep(state$value, n),
           cate_ci   = if (native_ci) cbind(lower = rep(state$value - 0.05, n),
                                             upper = rep(state$value + 0.05, n))
                       else NULL)
    }
  } else NULL
  gdpar_meta_learner_adapter(
    name = name,
    fit_predict_fun = fp,
    predict_fun = pf,
    native_ci = native_ci,
    description = sprintf("Mock adapter '%s' returning constant CATE", name)
  )
}

# ---- (i) adapter constructor / predicate / S3 print ---------------------

test_that("gdpar_meta_learner_adapter validates inputs", {
  expect_error(
    gdpar_meta_learner_adapter(name = "", fit_predict_fun = identity),
    class = "gdpar_input_error"
  )
  expect_error(
    gdpar_meta_learner_adapter(name = "x", fit_predict_fun = 42),
    class = "gdpar_input_error"
  )
  expect_error(
    gdpar_meta_learner_adapter(name = "x", fit_predict_fun = identity,
                                predict_fun = 1),
    class = "gdpar_input_error"
  )
  expect_error(
    gdpar_meta_learner_adapter(name = "x", fit_predict_fun = identity,
                                requires_r = 1),
    class = "gdpar_input_error"
  )
  expect_error(
    gdpar_meta_learner_adapter(name = "x", fit_predict_fun = identity,
                                native_ci = NA),
    class = "gdpar_input_error"
  )
  expect_error(
    gdpar_meta_learner_adapter(name = "x", fit_predict_fun = identity,
                                description = 42),
    class = "gdpar_input_error"
  )
})

test_that("gdpar_meta_learner_adapter constructs a valid object", {
  a <- make_mock_adapter()
  expect_true(is_gdpar_meta_learner_adapter(a))
  expect_false(is_gdpar_meta_learner_adapter(list()))
  expect_s3_class(a, "gdpar_meta_learner_adapter")
  expect_identical(a$name, "mock")
  expect_true(a$native_ci)
  expect_true(is.function(a$predict_fun))
})

test_that("print.gdpar_meta_learner_adapter does not error", {
  a <- make_mock_adapter()
  expect_output(print(a), "gdpar_meta_learner_adapter")
  expect_output(print(a), "mock")
})

# ---- (ii) internal helpers ----------------------------------------------

test_that(".compute_comparison_metrics returns symmetric matrices", {
  set.seed(1L)
  cate_list <- list(
    bridge = rnorm(50),
    grf    = rnorm(50),
    econml = rnorm(50)
  )
  m <- gdpar:::.compute_comparison_metrics(cate_list)
  expect_named(m, c("rmse", "pearson", "mad"))
  for (mat_name in names(m)) {
    expect_equal(dim(m[[mat_name]]), c(3L, 3L))
    expect_equal(rownames(m[[mat_name]]), c("bridge", "grf", "econml"))
  }
  expect_equal(diag(m$rmse), c(bridge = 0, grf = 0, econml = 0))
  expect_equal(diag(m$mad),  c(bridge = 0, grf = 0, econml = 0))
  expect_equal(diag(m$pearson), c(bridge = 1, grf = 1, econml = 1))
  expect_equal(m$rmse, t(m$rmse))
  expect_equal(m$mad,  t(m$mad))
  expect_equal(m$pearson, t(m$pearson))
})

test_that(".extract_covariates drops the outcome column", {
  df <- data.frame(y = 1:3, x1 = c(0.1, 0.2, 0.3), x2 = c(1L, 2L, 3L))
  out <- gdpar:::.extract_covariates(df, "y")
  expect_equal(colnames(out), c("x1", "x2"))
  expect_equal(nrow(out), 3L)
})

test_that(".guard_scalar_outcome rejects non-scalar bridges", {
  bridge <- list(meta = list(dim_kind = "multi"))
  class(bridge) <- "gdpar_causal_bridge"
  expect_error(gdpar:::.guard_scalar_outcome(bridge),
               class = "gdpar_unsupported_feature_error")
  bridge2 <- list(meta = list(dim_kind = "K_individual"))
  class(bridge2) <- "gdpar_causal_bridge"
  expect_error(gdpar:::.guard_scalar_outcome(bridge2),
               class = "gdpar_unsupported_feature_error")
  bridge3 <- list(meta = list(dim_kind = "scalar"))
  class(bridge3) <- "gdpar_causal_bridge"
  expect_silent(gdpar:::.guard_scalar_outcome(bridge3))
})

test_that(".validate_adapter_output enforces shape", {
  expect_error(
    gdpar:::.validate_adapter_output(list(), 5L, "x"),
    class = "gdpar_internal_error"
  )
  expect_error(
    gdpar:::.validate_adapter_output(list(cate_mean = 1:3), 5L, "x"),
    class = "gdpar_internal_error"
  )
  expect_error(
    gdpar:::.validate_adapter_output(
      list(cate_mean = 1:5,
           cate_ci = matrix(0, nrow = 5L, ncol = 3L)),
      5L, "x"),
    class = "gdpar_internal_error"
  )
  expect_silent(
    gdpar:::.validate_adapter_output(
      list(cate_mean = 1:5,
           cate_ci   = matrix(0, nrow = 5L, ncol = 2L)),
      5L, "x"
    )
  )
})

test_that(".check_adapter_requirements detects missing R packages", {
  a_bogus <- gdpar_meta_learner_adapter(
    name = "bogus",
    fit_predict_fun = function(X, Y, T, X_newdata, level, seed_run) {
      list(cate_mean = numeric(0L), cate_ci = NULL,
           state = NULL, notes = character(0L))
    },
    requires_r = "ThisPackageDoesNotExist____ZZ"
  )
  req <- gdpar:::.check_adapter_requirements(a_bogus)
  expect_false(req$ok)
  expect_true("ThisPackageDoesNotExist____ZZ" %in% req$missing_r)
})

# ---- (iii) orchestrator with mock bridge + mock adapter ------------------

test_that("gdpar_compare_meta_learners rejects non-bridge inputs", {
  expect_error(
    gdpar_compare_meta_learners(list(), list(make_mock_adapter())),
    class = "gdpar_input_error"
  )
})

test_that("gdpar_compare_meta_learners rejects empty methods list", {
  bridge <- make_mock_compare_bridge(
    newdata = data.frame(x1 = seq(-1, 1, length.out = 5L)),
    cate_mean = rep(0.5, 5L),
    data_t = data.frame(y = rnorm(10), x1 = rnorm(10)),
    data_c = data.frame(y = rnorm(10), x1 = rnorm(10))
  )
  expect_error(
    gdpar_compare_meta_learners(bridge, methods = list()),
    class = "gdpar_input_error"
  )
})

test_that("gdpar_compare_meta_learners rejects non-adapter methods", {
  bridge <- make_mock_compare_bridge(
    newdata = data.frame(x1 = seq(-1, 1, length.out = 5L)),
    cate_mean = rep(0.5, 5L),
    data_t = data.frame(y = rnorm(10), x1 = rnorm(10)),
    data_c = data.frame(y = rnorm(10), x1 = rnorm(10))
  )
  expect_error(
    gdpar_compare_meta_learners(bridge,
                                 methods = list(make_mock_adapter(), list())),
    class = "gdpar_input_error"
  )
})

test_that("gdpar_compare_meta_learners rejects duplicate method names", {
  bridge <- make_mock_compare_bridge(
    newdata = data.frame(x1 = seq(-1, 1, length.out = 5L)),
    cate_mean = rep(0.5, 5L),
    data_t = data.frame(y = rnorm(10), x1 = rnorm(10)),
    data_c = data.frame(y = rnorm(10), x1 = rnorm(10))
  )
  expect_error(
    gdpar_compare_meta_learners(bridge,
                                 methods = list(make_mock_adapter("a"),
                                                make_mock_adapter("a"))),
    class = "gdpar_input_error"
  )
})

test_that("gdpar_compare_meta_learners returns a valid comparison", {
  set.seed(42L)
  n <- 30L
  data_t <- data.frame(y = rnorm(n), x1 = rnorm(n))
  data_c <- data.frame(y = rnorm(n), x1 = rnorm(n))
  newdata <- data.frame(x1 = seq(-1, 1, length.out = 5L))
  bridge <- make_mock_compare_bridge(
    newdata = newdata,
    cate_mean = rep(0.5, 5L),
    data_t = data_t, data_c = data_c
  )
  cmp <- gdpar_compare_meta_learners(
    bridge,
    methods = list(make_mock_adapter("m1", cate_value = 0.4),
                   make_mock_adapter("m2", cate_value = 0.6))
  )
  expect_s3_class(cmp, "gdpar_meta_learner_comparison")
  expect_equal(cmp$n_obs, 5L)
  expect_equal(cmp$n_methods, 2L)
  expect_equal(names(cmp$external), c("m1", "m2"))
  expect_equal(rownames(cmp$comparison$rmse),
               c("bridge", "m1", "m2"))
  expect_equal(cmp$comparison$rmse["bridge", "m1"], 0.1, tolerance = 1e-10)
  expect_equal(cmp$comparison$rmse["bridge", "m2"], 0.1, tolerance = 1e-10)
  expect_equal(cmp$comparison$rmse["m1", "m2"],     0.2, tolerance = 1e-10)
  expect_equal(cmp$comparison$mad["m1", "m2"],      0.2, tolerance = 1e-10)
})

test_that("gdpar_compare_meta_learners uses adapter$name when methods unnamed", {
  set.seed(42L)
  n <- 20L
  bridge <- make_mock_compare_bridge(
    newdata = data.frame(x1 = seq(-1, 1, length.out = 4L)),
    cate_mean = rep(0.5, 4L),
    data_t = data.frame(y = rnorm(n), x1 = rnorm(n)),
    data_c = data.frame(y = rnorm(n), x1 = rnorm(n))
  )
  cmp <- gdpar_compare_meta_learners(
    bridge,
    methods = list(make_mock_adapter("alpha"))
  )
  expect_equal(names(cmp$external), "alpha")
})

test_that("gdpar_compare_meta_learners accepts explicit data argument", {
  n <- 20L
  bridge <- make_mock_compare_bridge(
    newdata = data.frame(x1 = seq(-1, 1, length.out = 4L)),
    cate_mean = rep(0.5, 4L),
    data_t = data.frame(y = rnorm(n), x1 = rnorm(n)),
    data_c = data.frame(y = rnorm(n), x1 = rnorm(n))
  )
  user_data <- list(
    X = data.frame(x1 = rnorm(2L * n)),
    T = c(rep(1L, n), rep(0L, n)),
    Y = rnorm(2L * n)
  )
  cmp <- gdpar_compare_meta_learners(
    bridge,
    methods = list(make_mock_adapter("m1")),
    data = user_data
  )
  expect_s3_class(cmp, "gdpar_meta_learner_comparison")
})

test_that("gdpar_compare_meta_learners rejects malformed user data", {
  n <- 20L
  bridge <- make_mock_compare_bridge(
    newdata = data.frame(x1 = seq(-1, 1, length.out = 4L)),
    cate_mean = rep(0.5, 4L),
    data_t = data.frame(y = rnorm(n), x1 = rnorm(n)),
    data_c = data.frame(y = rnorm(n), x1 = rnorm(n))
  )
  expect_error(
    gdpar_compare_meta_learners(
      bridge,
      methods = list(make_mock_adapter("m1")),
      data = list(X = data.frame(x1 = 1:3))
    ),
    class = "gdpar_input_error"
  )
  expect_error(
    gdpar_compare_meta_learners(
      bridge,
      methods = list(make_mock_adapter("m1")),
      data = list(X = data.frame(x1 = 1:3), T = c(0L, 1L, 2L),
                  Y = rep(0, 3L))
    ),
    class = "gdpar_input_error"
  )
})

test_that("gdpar_compare_meta_learners aborts when bridge is non-scalar", {
  n <- 20L
  bridge <- make_mock_compare_bridge(
    newdata = data.frame(x1 = seq(-1, 1, length.out = 4L)),
    cate_mean = rep(0.5, 4L),
    data_t = data.frame(y = rnorm(n), x1 = rnorm(n)),
    data_c = data.frame(y = rnorm(n), x1 = rnorm(n)),
    dim_kind = "multi"
  )
  expect_error(
    gdpar_compare_meta_learners(bridge,
                                 methods = list(make_mock_adapter("m1"))),
    class = "gdpar_unsupported_feature_error"
  )
})

test_that("gdpar_compare_meta_learners aborts on missing adapter requirements", {
  n <- 20L
  bridge <- make_mock_compare_bridge(
    newdata = data.frame(x1 = seq(-1, 1, length.out = 4L)),
    cate_mean = rep(0.5, 4L),
    data_t = data.frame(y = rnorm(n), x1 = rnorm(n)),
    data_c = data.frame(y = rnorm(n), x1 = rnorm(n))
  )
  bogus <- gdpar_meta_learner_adapter(
    name = "bogus",
    fit_predict_fun = function(X, Y, T, X_newdata, level, seed_run) {
      list(cate_mean = rep(0, nrow(X_newdata)), cate_ci = NULL,
           state = NULL, notes = character(0L))
    },
    requires_r = "ThisPackageDoesNotExist____ZZ"
  )
  expect_error(
    gdpar_compare_meta_learners(bridge, methods = list(bogus)),
    class = "gdpar_missing_dependency_error"
  )
})

# ---- (vi) S3 methods on a real (mock) comparison ------------------------

test_that("print.gdpar_meta_learner_comparison emits a summary block", {
  n <- 20L
  bridge <- make_mock_compare_bridge(
    newdata = data.frame(x1 = seq(-1, 1, length.out = 4L)),
    cate_mean = rep(0.5, 4L),
    data_t = data.frame(y = rnorm(n), x1 = rnorm(n)),
    data_c = data.frame(y = rnorm(n), x1 = rnorm(n))
  )
  cmp <- gdpar_compare_meta_learners(
    bridge,
    methods = list(make_mock_adapter("m1"), make_mock_adapter("m2"))
  )
  expect_output(print(cmp), "gdpar_meta_learner_comparison")
  expect_output(print(cmp), "RMSE")
  expect_output(print(cmp), "Pearson")
})

test_that("summary.gdpar_meta_learner_comparison returns structured object", {
  n <- 20L
  bridge <- make_mock_compare_bridge(
    newdata = data.frame(x1 = seq(-1, 1, length.out = 4L)),
    cate_mean = rep(0.5, 4L),
    data_t = data.frame(y = rnorm(n), x1 = rnorm(n)),
    data_c = data.frame(y = rnorm(n), x1 = rnorm(n))
  )
  cmp <- gdpar_compare_meta_learners(
    bridge,
    methods = list(make_mock_adapter("m1"), make_mock_adapter("m2"))
  )
  s <- summary(cmp)
  expect_s3_class(s, "summary.gdpar_meta_learner_comparison")
  expect_named(s, c("ate_table", "metrics", "timing", "level",
                    "n_obs", "n_methods"))
  expect_equal(nrow(s$ate_table), 3L)
  expect_equal(nrow(s$metrics),
               3L * 2L)
  expect_output(print(s), "ATE per method")
})

# ---- (iv) grf reference adapter (gated by Suggests:grf) ----------------

test_that("gdpar_adapter_grf constructs a valid adapter when grf installed", {
  skip_if_not_installed("grf")
  a <- gdpar_adapter_grf(num_trees = 200L)
  expect_s3_class(a, "gdpar_meta_learner_adapter")
  expect_identical(a$requires_r, "grf")
  expect_true(a$native_ci)
  expect_true(is.function(a$predict_fun))
})

test_that("gdpar_adapter_grf rejects malformed hyperparameters", {
  skip_if_not_installed("grf")
  expect_error(gdpar_adapter_grf(num_trees = -1),
               class = "gdpar_input_error")
  expect_error(gdpar_adapter_grf(sample_fraction = 0.9),
               class = "gdpar_input_error")
  expect_error(gdpar_adapter_grf(honesty = NA),
               class = "gdpar_input_error")
})

test_that("gdpar_adapter_grf smoke on synthetic data (algorithmic equivalence)", {
  skip_if_not_installed("grf")
  set.seed(2026L)
  n <- 250L
  beta_x_t <- 1.2; beta_x_c <- 0.6
  x_t <- rnorm(n); x_c <- rnorm(n)
  y_t <- 0.5 + beta_x_t * x_t + rnorm(n, sd = 0.3)
  y_c <- 0.0 + beta_x_c * x_c + rnorm(n, sd = 0.3)
  data_t <- data.frame(y = y_t, x1 = x_t)
  data_c <- data.frame(y = y_c, x1 = x_c)
  newdata <- data.frame(x1 = seq(-2, 2, length.out = 25L))
  cate_pop <- (0.5 + beta_x_t * newdata$x1) -
              (0.0 + beta_x_c * newdata$x1)
  bridge <- make_mock_compare_bridge(
    newdata = newdata,
    cate_mean = cate_pop + rnorm(length(cate_pop), sd = 0.05),
    data_t = data_t, data_c = data_c
  )
  cmp <- gdpar_compare_meta_learners(
    bridge,
    methods = list(grf = gdpar_adapter_grf(num_trees = 500L,
                                            seed = 2026L))
  )
  expect_s3_class(cmp, "gdpar_meta_learner_comparison")
  expect_true(!is.null(cmp$external$grf$cate_ci))
  rmse_grf_pop <- sqrt(mean((cmp$external$grf$cate_mean - cate_pop)^2))
  scale_y <- sd(c(y_t, y_c))
  expect_lt(rmse_grf_pop, 0.5 * scale_y)
})

test_that("predict.gdpar_meta_learner_comparison reuses grf state without refit", {
  skip_if_not_installed("grf")
  set.seed(2027L)
  n <- 200L
  x_t <- rnorm(n); x_c <- rnorm(n)
  y_t <- 0.3 + 0.8 * x_t + rnorm(n, sd = 0.3)
  y_c <- 0.0 + 0.4 * x_c + rnorm(n, sd = 0.3)
  bridge <- make_mock_compare_bridge(
    newdata = data.frame(x1 = seq(-2, 2, length.out = 11L)),
    cate_mean = rep(0.3, 11L),
    data_t = data.frame(y = y_t, x1 = x_t),
    data_c = data.frame(y = y_c, x1 = x_c)
  )
  cmp <- gdpar_compare_meta_learners(
    bridge,
    methods = list(grf = gdpar_adapter_grf(num_trees = 300L, seed = 2027L))
  )
  newdata2 <- data.frame(x1 = seq(-1.5, 1.5, length.out = 7L))
  pred <- expect_silent(
    predict(cmp, newdata = newdata2)
  )
  expect_s3_class(pred, "predict.gdpar_meta_learner_comparison")
  expect_equal(length(pred$external$grf$cate_mean), 7L)
  expect_true(!is.null(pred$external$grf$cate_ci))
})

test_that("predict.gdpar_meta_learner_comparison warns on missing predict_fun", {
  set.seed(2028L)
  n <- 20L
  bridge <- make_mock_compare_bridge(
    newdata = data.frame(x1 = seq(-1, 1, length.out = 4L)),
    cate_mean = rep(0.5, 4L),
    data_t = data.frame(y = rnorm(n), x1 = rnorm(n)),
    data_c = data.frame(y = rnorm(n), x1 = rnorm(n))
  )
  cmp <- gdpar_compare_meta_learners(
    bridge,
    methods = list(make_mock_adapter("m_nopred", has_predict = FALSE))
  )
  newdata2 <- data.frame(x1 = seq(-0.5, 0.5, length.out = 3L))
  expect_warning(
    predict(cmp, newdata = newdata2),
    class = "gdpar_diagnostic_warning"
  )
})

# ---- (v) EconML smoke (env-var gated) -----------------------------------

test_that("gdpar_adapter_econml constructs a valid adapter when reticulate installed", {
  skip_if_not_installed("reticulate")
  a <- gdpar_adapter_econml(n_estimators = 100L)
  expect_s3_class(a, "gdpar_meta_learner_adapter")
  expect_identical(a$requires_r, "reticulate")
  expect_identical(a$requires_py, "econml")
  expect_true(a$native_ci)
  expect_true(is.function(a$predict_fun))
})

test_that("gdpar_adapter_econml rejects unsupported estimators", {
  expect_error(
    gdpar_adapter_econml(estimator = "Frankenstein"),
    class = "gdpar_unsupported_feature_error"
  )
  expect_error(
    gdpar_adapter_econml(estimator = ""),
    class = "gdpar_input_error"
  )
})

test_that("gdpar_adapter_econml smoke (env-var gated)", {
  if (!nzchar(Sys.getenv("GDPAR_RUN_PYTHON_SMOKE_ECONML"))) {
    skip("Set GDPAR_RUN_PYTHON_SMOKE_ECONML=1 to run the EconML smoke test.")
  }
  skip_if_not_installed("reticulate")
  if (exists("py_require", envir = asNamespace("reticulate"),
             inherits = FALSE)) {
    tryCatch(reticulate::py_require("econml"),
             error = function(e) NULL)
  }
  if (!isTRUE(tryCatch(reticulate::py_module_available("econml"),
                       error = function(e) FALSE))) {
    skip("Python module 'econml' not available.")
  }
  set.seed(2029L)
  n <- 150L
  x_t <- rnorm(n); x_c <- rnorm(n)
  y_t <- 0.4 + 0.9 * x_t + rnorm(n, sd = 0.3)
  y_c <- 0.0 + 0.5 * x_c + rnorm(n, sd = 0.3)
  bridge <- make_mock_compare_bridge(
    newdata = data.frame(x1 = seq(-1.5, 1.5, length.out = 9L)),
    cate_mean = rep(0.4, 9L),
    data_t = data.frame(y = y_t, x1 = x_t),
    data_c = data.frame(y = y_c, x1 = x_c)
  )
  cmp <- gdpar_compare_meta_learners(
    bridge,
    methods = list(econml = gdpar_adapter_econml(n_estimators = 200L,
                                                  seed = 2029L))
  )
  expect_s3_class(cmp, "gdpar_meta_learner_comparison")
  expect_equal(length(cmp$external$econml$cate_mean), 9L)
})

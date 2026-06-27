# ============================================================
# gdpar / inst/benchmarks/benchmark_competitors.R
# Comparison of gdpar against canonical R competitors.
# Two distinct comparison groups, reported separately:
#
# Group A (structural): mgcv::gam, rstanarm::stan_glmer, brms::brm
#   if installed. These methods have conceptual correspondence with
#   the AMM canonical form (additive/varying-coefficient structure)
#   and are therefore informative comparisons.
#
# Group B (predictive agnostic): ranger::ranger (random forest),
#   xgboost::xgboost (gradient boosting). These methods have no
#   structural correspondence with AMM but are widely used as
#   black-box predictors. The comparison here is pure predictive
#   performance under the same metrics; the methodological reading
#   is "structural interpretable model versus black-box predictor".
#
# Designed to take ~30-60 minutes total on a modern multi-core CPU
# depending on which competitors are installed.
# ============================================================

if (!exists("bench_section", mode = "function")) {
  helpers_path <- local({
    candidates <- character(0)
    if (requireNamespace("rstudioapi", quietly = TRUE) &&
        rstudioapi::isAvailable()) {
      ctx <- tryCatch(rstudioapi::getSourceEditorContext()$path,
                      error = function(e) "")
      if (nzchar(ctx) && file.exists(ctx)) {
        candidates <- c(candidates,
                        file.path(dirname(ctx), "helpers.R"))
      }
    }
    sourced <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
    if (!is.null(sourced) && nzchar(sourced)) {
      candidates <- c(candidates,
                      file.path(dirname(normalizePath(sourced,
                                                     mustWork = FALSE)),
                                "helpers.R"))
    }
    candidates <- c(candidates,
                    file.path(getwd(), "helpers.R"),
                    file.path(getwd(), "inst", "benchmarks",
                              "helpers.R"))
    inst_path <- system.file("benchmarks", "helpers.R",
                             package = "gdpar")
    if (nzchar(inst_path)) candidates <- c(candidates, inst_path)
    found <- ""
    for (p in candidates) {
      if (nzchar(p) && file.exists(p)) { found <- p; break }
    }
    if (!nzchar(found)) {
      stop("Could not locate helpers.R. Set the working directory ",
           "to the gdpar package root or to its inst/benchmarks ",
           "subdirectory, or install gdpar so that system.file() ",
           "resolves it.")
    }
    found
  })
  source(helpers_path)
}

bench_section("gdpar versus competitors benchmark")

if (!bench_require("gdpar")) stop("gdpar must be installed first.")
library(gdpar)
have_cmdstan <- bench_have_cmdstan()
have_mgcv <- bench_require("mgcv")
have_rstanarm <- bench_require("rstanarm")
have_brms <- bench_require("brms")
have_ranger <- bench_require("ranger")
have_xgboost <- bench_require("xgboost")

if (!have_cmdstan) {
  bench_skip("All competitor benchmarks require cmdstan for gdpar; aborting.")
  invisible(NULL)
} else {

# ------------------------------------------------------------
# Common train/test split
# ------------------------------------------------------------

set.seed(2026L)
N <- 800L
train_frac <- 0.75
n_train <- floor(train_frac * N)

# ------------------------------------------------------------
# Dataset 1: Gaussian, well-conditioned
# ------------------------------------------------------------

bench_subsection("Dataset 1: Gaussian, well-conditioned (n = 800, 75/25 split)")

df_g <- make_gauss_clean(n = N, seed = 11L)
idx_train <- sample.int(N, n_train)
df_g_train <- df_g[idx_train, ]
df_g_test  <- df_g[-idx_train, ]

results_g <- list()

# gdpar
fit_g_gdpar <- gdpar(
  formula       = y ~ x1 + x2,
  family        = gdpar_family("gaussian"),
  amm           = amm_spec(a = ~ x1 + x2),
  data          = df_g_train,
  iter_warmup   = 500,
  iter_sampling = 500,
  chains        = 4,
  refresh       = 0,
  verbose       = FALSE
)
pp_gdpar <- bench_posterior_predict(fit_g_gdpar, df_g_test, "gaussian")
y_pred_mean_gdpar <- colMeans(pp_gdpar)
results_g[[length(results_g) + 1L]] <- build_compare_df(
  "gdpar(bayes, AMM L1)",
  list(
    rmse = rmse(df_g_test$y, y_pred_mean_gdpar),
    mae = mae(df_g_test$y, y_pred_mean_gdpar),
    crps = mean(crps_ensemble(df_g_test$y, pp_gdpar))
  )
)

# Group A: mgcv::gam
if (have_mgcv) {
  fit_g_gam <- mgcv::gam(y ~ s(x1) + s(x2), data = df_g_train,
                          method = "REML")
  y_pred_gam <- predict(fit_g_gam, newdata = df_g_test)
  results_g[[length(results_g) + 1L]] <- build_compare_df(
    "mgcv::gam (s(x1) + s(x2))",
    list(
      rmse = rmse(df_g_test$y, y_pred_gam),
      mae = mae(df_g_test$y, y_pred_gam),
      crps = NA_real_
    )
  )
}

# Group A: rstanarm
if (have_rstanarm) {
  fit_g_rsa <- rstanarm::stan_glm(
    y ~ x1 + x2, data = df_g_train, family = gaussian(),
    chains = 4, iter = 1000, refresh = 0
  )
  pp_rsa <- rstanarm::posterior_predict(fit_g_rsa, newdata = df_g_test)
  y_pred_rsa <- colMeans(pp_rsa)
  results_g[[length(results_g) + 1L]] <- build_compare_df(
    "rstanarm::stan_glm",
    list(
      rmse = rmse(df_g_test$y, y_pred_rsa),
      mae = mae(df_g_test$y, y_pred_rsa),
      crps = mean(crps_ensemble(df_g_test$y, pp_rsa))
    )
  )
}

# Group A: brms
if (have_brms) {
  fit_g_brms <- brms::brm(
    y ~ x1 + x2, data = df_g_train, family = gaussian(),
    chains = 4, iter = 1000, refresh = 0, silent = 2
  )
  pp_brms <- brms::posterior_predict(fit_g_brms, newdata = df_g_test)
  y_pred_brms <- colMeans(pp_brms)
  results_g[[length(results_g) + 1L]] <- build_compare_df(
    "brms::brm",
    list(
      rmse = rmse(df_g_test$y, y_pred_brms),
      mae = mae(df_g_test$y, y_pred_brms),
      crps = mean(crps_ensemble(df_g_test$y, pp_brms))
    )
  )
}

# Group B: ranger
if (have_ranger) {
  fit_g_ranger <- ranger::ranger(y ~ x1 + x2, data = df_g_train,
                                   num.trees = 500, quantreg = TRUE)
  pred_ranger <- predict(fit_g_ranger, data = df_g_test)$predictions
  qpred_ranger <- predict(fit_g_ranger, data = df_g_test,
                            type = "quantiles",
                            quantiles = seq(0.05, 0.95, by = 0.05))
  results_g[[length(results_g) + 1L]] <- build_compare_df(
    "ranger (RF, 500 trees)",
    list(
      rmse = rmse(df_g_test$y, pred_ranger),
      mae = mae(df_g_test$y, pred_ranger),
      crps = NA_real_
    )
  )
}

# Group B: xgboost
if (have_xgboost) {
  X_train <- as.matrix(df_g_train[, c("x1", "x2")])
  X_test  <- as.matrix(df_g_test[, c("x1", "x2")])
  dtrain <- xgboost::xgb.DMatrix(data = X_train, label = df_g_train$y)
  dtest  <- xgboost::xgb.DMatrix(data = X_test)
  fit_g_xgb <- xgboost::xgb.train(
    data = dtrain,
    params = list(objective = "reg:squarederror"),
    nrounds = 200,
    verbose = 0
  )
  pred_xgb <- predict(fit_g_xgb, dtest)
  results_g[[length(results_g) + 1L]] <- build_compare_df(
    "xgboost (200 rounds)",
    list(
      rmse = rmse(df_g_test$y, pred_xgb),
      mae = mae(df_g_test$y, pred_xgb),
      crps = NA_real_
    )
  )
}

bench_subsection("Dataset 1 results")
df_results_g <- do.call(rbind, results_g)
bench_table(df_results_g, digits = 4)

# ------------------------------------------------------------
# Dataset 2: Negative binomial, well-conditioned
# ------------------------------------------------------------

bench_subsection("Dataset 2: Negative binomial, well-conditioned (n = 800, 75/25 split)")

df_n <- make_negbin_clean(n = N, phi = 5, seed = 22L)
df_n_train <- df_n[idx_train, ]
df_n_test  <- df_n[-idx_train, ]

results_n <- list()

# gdpar
fit_n_gdpar <- gdpar(
  formula       = y ~ x1 + x2,
  family        = gdpar_family("neg_binomial_2"),
  amm           = amm_spec(a = ~ x1 + x2),
  data          = df_n_train,
  iter_warmup   = 500,
  iter_sampling = 500,
  chains        = 4,
  refresh       = 0,
  verbose       = FALSE
)
pp_n_gdpar <- bench_posterior_predict(fit_n_gdpar, df_n_test,
                                       "neg_binomial_2")
y_pred_n_gdpar <- colMeans(pp_n_gdpar)
results_n[[length(results_n) + 1L]] <- build_compare_df(
  "gdpar(bayes, NB AMM L1)",
  list(
    rmse = rmse(df_n_test$y, y_pred_n_gdpar),
    mae = mae(df_n_test$y, y_pred_n_gdpar),
    crps = mean(crps_ensemble(df_n_test$y, pp_n_gdpar))
  )
)

# Group A: mgcv (NB family)
if (have_mgcv) {
  fit_n_gam <- mgcv::gam(y ~ s(x1) + s(x2), data = df_n_train,
                          family = mgcv::nb(), method = "REML")
  pred_n_gam <- predict(fit_n_gam, newdata = df_n_test, type = "response")
  results_n[[length(results_n) + 1L]] <- build_compare_df(
    "mgcv::gam (NB, s(x1) + s(x2))",
    list(
      rmse = rmse(df_n_test$y, pred_n_gam),
      mae = mae(df_n_test$y, pred_n_gam),
      crps = NA_real_
    )
  )
}

# Group A: rstanarm (NB family)
if (have_rstanarm) {
  fit_n_rsa <- rstanarm::stan_glm(
    y ~ x1 + x2, data = df_n_train,
    family = rstanarm::neg_binomial_2(),
    chains = 4, iter = 1000, refresh = 0
  )
  pp_n_rsa <- rstanarm::posterior_predict(fit_n_rsa, newdata = df_n_test)
  y_pred_n_rsa <- colMeans(pp_n_rsa)
  results_n[[length(results_n) + 1L]] <- build_compare_df(
    "rstanarm::stan_glm (NB)",
    list(
      rmse = rmse(df_n_test$y, y_pred_n_rsa),
      mae = mae(df_n_test$y, y_pred_n_rsa),
      crps = mean(crps_ensemble(df_n_test$y, pp_n_rsa))
    )
  )
}

# Group B: ranger
if (have_ranger) {
  fit_n_ranger <- ranger::ranger(y ~ x1 + x2, data = df_n_train,
                                   num.trees = 500)
  pred_n_ranger <- predict(fit_n_ranger, data = df_n_test)$predictions
  results_n[[length(results_n) + 1L]] <- build_compare_df(
    "ranger (RF, 500 trees)",
    list(
      rmse = rmse(df_n_test$y, pred_n_ranger),
      mae = mae(df_n_test$y, pred_n_ranger),
      crps = NA_real_
    )
  )
}

# Group B: xgboost (Poisson objective)
if (have_xgboost) {
  X_train <- as.matrix(df_n_train[, c("x1", "x2")])
  X_test  <- as.matrix(df_n_test[, c("x1", "x2")])
  dtrain <- xgboost::xgb.DMatrix(data = X_train, label = df_n_train$y)
  dtest  <- xgboost::xgb.DMatrix(data = X_test)
  fit_n_xgb <- xgboost::xgb.train(
    data = dtrain,
    params = list(objective = "count:poisson"),
    nrounds = 200,
    verbose = 0
  )
  pred_n_xgb <- predict(fit_n_xgb, dtest)
  results_n[[length(results_n) + 1L]] <- build_compare_df(
    "xgboost (Poisson, 200 rounds)",
    list(
      rmse = rmse(df_n_test$y, pred_n_xgb),
      mae = mae(df_n_test$y, pred_n_xgb),
      crps = NA_real_
    )
  )
}

bench_subsection("Dataset 2 results")
df_results_n <- do.call(rbind, results_n)
bench_table(df_results_n, digits = 4)

# ------------------------------------------------------------
# Dataset 3: Bernoulli (binary classification)
# ------------------------------------------------------------

bench_subsection("Dataset 3: Bernoulli (n = 800, 75/25 split)")

df_b <- make_bernoulli_clean(n = N, seed = 33L)
df_b_train <- df_b[idx_train, ]
df_b_test  <- df_b[-idx_train, ]

results_b <- list()

# gdpar
fit_b_gdpar <- gdpar(
  formula       = y ~ x1 + x2,
  family        = gdpar_family("bernoulli"),
  amm           = amm_spec(a = ~ x1 + x2),
  data          = df_b_train,
  iter_warmup   = 500,
  iter_sampling = 500,
  chains        = 4,
  refresh       = 0,
  verbose       = FALSE
)
pp_b_gdpar <- predict(fit_b_gdpar, newdata = df_b_test,
                        type = "response", summary = "draws")
p_pred_gdpar <- colMeans(pp_b_gdpar)
results_b[[length(results_b) + 1L]] <- build_compare_df(
  "gdpar(bayes, Bernoulli AMM L1)",
  list(
    auc = auc_roc(df_b_test$y, p_pred_gdpar),
    log_loss = log_loss(df_b_test$y, p_pred_gdpar),
    brier = brier_score(df_b_test$y, p_pred_gdpar)
  )
)

if (have_mgcv) {
  fit_b_gam <- mgcv::gam(y ~ s(x1) + s(x2), data = df_b_train,
                          family = binomial(), method = "REML")
  pred_b_gam <- predict(fit_b_gam, newdata = df_b_test,
                          type = "response")
  results_b[[length(results_b) + 1L]] <- build_compare_df(
    "mgcv::gam (binomial, s(x1) + s(x2))",
    list(
      auc = auc_roc(df_b_test$y, pred_b_gam),
      log_loss = log_loss(df_b_test$y, pred_b_gam),
      brier = brier_score(df_b_test$y, pred_b_gam)
    )
  )
}

if (have_rstanarm) {
  fit_b_rsa <- rstanarm::stan_glm(
    y ~ x1 + x2, data = df_b_train, family = binomial(),
    chains = 4, iter = 1000, refresh = 0
  )
  pred_b_rsa <- rstanarm::posterior_epred(fit_b_rsa,
                                            newdata = df_b_test)
  p_pred_rsa <- colMeans(pred_b_rsa)
  results_b[[length(results_b) + 1L]] <- build_compare_df(
    "rstanarm::stan_glm (binomial)",
    list(
      auc = auc_roc(df_b_test$y, p_pred_rsa),
      log_loss = log_loss(df_b_test$y, p_pred_rsa),
      brier = brier_score(df_b_test$y, p_pred_rsa)
    )
  )
}

if (have_ranger) {
  fit_b_ranger <- ranger::ranger(
    factor(y) ~ x1 + x2, data = df_b_train,
    num.trees = 500, probability = TRUE
  )
  p_pred_ranger <- predict(fit_b_ranger, data = df_b_test)$predictions[, "1"]
  results_b[[length(results_b) + 1L]] <- build_compare_df(
    "ranger (RF, 500 trees)",
    list(
      auc = auc_roc(df_b_test$y, p_pred_ranger),
      log_loss = log_loss(df_b_test$y, p_pred_ranger),
      brier = brier_score(df_b_test$y, p_pred_ranger)
    )
  )
}

if (have_xgboost) {
  X_train <- as.matrix(df_b_train[, c("x1", "x2")])
  X_test  <- as.matrix(df_b_test[, c("x1", "x2")])
  dtrain <- xgboost::xgb.DMatrix(data = X_train, label = df_b_train$y)
  dtest  <- xgboost::xgb.DMatrix(data = X_test)
  fit_b_xgb <- xgboost::xgb.train(
    data = dtrain,
    params = list(objective = "binary:logistic"),
    nrounds = 200,
    verbose = 0
  )
  p_pred_xgb <- predict(fit_b_xgb, dtest)
  results_b[[length(results_b) + 1L]] <- build_compare_df(
    "xgboost (binary:logistic, 200 rounds)",
    list(
      auc = auc_roc(df_b_test$y, p_pred_xgb),
      log_loss = log_loss(df_b_test$y, p_pred_xgb),
      brier = brier_score(df_b_test$y, p_pred_xgb)
    )
  )
}

bench_subsection("Dataset 3 results")
df_results_b <- do.call(rbind, results_b)
bench_table(df_results_b, digits = 4)

# ------------------------------------------------------------
# Final consolidated report
# ------------------------------------------------------------

bench_section("Competitor benchmark summary")

cat("Dataset 1 (Gaussian):\n")
bench_table(df_results_g, digits = 4)
cat("\nDataset 2 (Negative Binomial):\n")
bench_table(df_results_n, digits = 4)
cat("\nDataset 3 (Bernoulli):\n")
bench_table(df_results_b, digits = 4)

cat("\nMethodological notes on the comparison:\n")
cat("\nGroup A (structural: mgcv, rstanarm, brms):\n")
cat("- These methods share the structural framing of gdpar (additive / varying-coefficient on the linear-predictor scale of a generalized linear model). Differences in predictive metrics reflect parametrization choices, prior structure, and basis specifics.\n")
cat("- A reasonable expectation is that gdpar performs comparably to rstanarm on Level 1 specifications, and that mgcv with smooth bases can outperform any parametric method when the truth is non-linear in covariates.\n")
cat("\nGroup B (predictive agnostic: ranger, xgboost):\n")
cat("- These methods have no structural correspondence with the AMM canonical form. They are flexible function approximators that may fit the conditional mean very accurately at the cost of interpretability and uncertainty quantification.\n")
cat("- A reasonable expectation is that ranger and xgboost can win on point-prediction metrics (RMSE, MAE) on datasets with strong non-linearities, while gdpar wins on probabilistic metrics (CRPS, ELPD-LOO) and on inferential outputs (calibrated credible intervals, parameter-level interpretation).\n")
cat("- The comparison is best read as: 'when do you need structural interpretation and uncertainty quantification, and when do you only need accurate point predictions?'\n")

assign("competitor_results", list(
  gaussian = df_results_g,
  negbin = df_results_n,
  bernoulli = df_results_b
), envir = .GlobalEnv)
invisible(NULL)
}

# ============================================================
# gdpar / inst/benchmarks/helpers.R
# Shared utilities for the validation and benchmarking scripts.
# This file is a standalone helper, not part of the package API.
# Sourced by quick_smoke.R, synthetic_recovery.R, synthetic_hard.R,
# benchmark_competitors.R, ebird_validation.R, and full_validation.R.
# ============================================================

# ------------------------------------------------------------
# Section / status reporting
# ------------------------------------------------------------

bench_section <- function(title) {
  bar <- paste(rep("=", 70), collapse = "")
  cat("\n", bar, "\n", title, "\n", bar, "\n", sep = "")
}

bench_subsection <- function(title) {
  cat("\n--- ", title, " ---\n", sep = "")
}

bench_status <- function(...) {
  cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n", sep = "")
}

bench_skip <- function(reason) {
  cat("[skip] ", reason, "\n", sep = "")
}

bench_pass <- function(msg) {
  cat("[ ok ] ", msg, "\n", sep = "")
}

bench_fail <- function(msg) {
  cat("[FAIL] ", msg, "\n", sep = "")
}

# ------------------------------------------------------------
# Conditional package loading
# ------------------------------------------------------------

bench_require <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    bench_skip(sprintf("package '%s' not installed; sections that require it will be skipped.", pkg))
    return(FALSE)
  }
  TRUE
}

bench_have_cmdstan <- function() {
  if (!bench_require("cmdstanr")) return(FALSE)
  cs_path <- tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL)
  if (is.null(cs_path) || !nzchar(cs_path)) {
    bench_skip("cmdstan binary not configured (run cmdstanr::install_cmdstan() first).")
    return(FALSE)
  }
  TRUE
}

# ------------------------------------------------------------
# Synthetic data generators with controlled difficulty
# ------------------------------------------------------------

# Standard well-conditioned Gaussian AMM Level 1 dataset.
# Coefficients on x1 and x2 are anti-symmetric by historical
# construction (legacy from the prior Stan implementation that imposed
# a sum-to-zero reparametrisation on basis coefficients). After the
# (C2)/(C3) implementation fix --- centering of Z_a/Z_b in R is
# sufficient to satisfy E_mu[a(X)] = E_mu[b(X)] = 0 --- the model
# recovers any a_coef in R^{J_a}; the anti-symmetry of the truths is
# no longer required for recoverability and is retained only for
# reproducibility of legacy benchmarks.
make_gauss_clean <- function(n = 500, seed = 1L) {
  set.seed(seed)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  eta <- 1 + 0.4 * x1 - 0.4 * x2
  y <- eta + rnorm(n, sd = 0.5)
  data.frame(x1 = x1, x2 = x2, y = y)
}

# Poisson AMM Level 1 dataset, well-conditioned. Sum-to-zero truths.
make_poisson_clean <- function(n = 500, seed = 2L) {
  set.seed(seed)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  eta <- 1 + 0.25 * x1 - 0.25 * x2
  y <- rpois(n, exp(eta))
  data.frame(x1 = x1, x2 = x2, y = y)
}

# Negative binomial AMM Level 1 dataset. Sum-to-zero truths.
make_negbin_clean <- function(n = 500, phi = 5, seed = 3L) {
  set.seed(seed)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  eta <- 1 + 0.25 * x1 - 0.25 * x2
  mu <- exp(eta)
  y <- rnbinom(n, size = phi, mu = mu)
  data.frame(x1 = x1, x2 = x2, y = y)
}

# Bernoulli AMM Level 1 dataset. Sum-to-zero truths.
make_bernoulli_clean <- function(n = 500, seed = 4L) {
  set.seed(seed)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  eta <- 0 + 0.55 * x1 - 0.55 * x2
  p <- 1 / (1 + exp(-eta))
  y <- rbinom(n, 1, p)
  data.frame(x1 = x1, x2 = x2, y = y)
}

# AMM Level 2 with active multiplicative component.
# a acts on (x1, x2); b acts on (x3, x4); W acts on (x1, x2).
# Disjoint covariate sets for a and b yield a block-diagonal Gram
# matrix at any theta_ref, so the basis-restricted identifiability
# check passes. Covariates x3, x4 are independent of x1, x2.
make_gauss_amm_level2 <- function(n = 500, theta_ref = 1.5,
                                  a_coef = c(0.4, -0.4),
                                  b_coef = c(0.3, -0.3),
                                  W_coef = c(0.3, -0.2),
                                  sd_y = 0.5, seed = 5L) {
  set.seed(seed)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x3 <- rnorm(n)
  x4 <- rnorm(n)
  Z_a <- cbind(x1 - mean(x1), x2 - mean(x2))
  Z_b <- cbind(x3 - mean(x3), x4 - mean(x4))
  X_W <- cbind(scale(x1)[, 1], scale(x2)[, 1])
  W_diff <- (theta_ref - 0)
  eta <- theta_ref +
    Z_a %*% a_coef +
    (Z_b %*% b_coef) * theta_ref +
    W_diff * (X_W %*% W_coef)
  y <- as.numeric(eta) + rnorm(n, sd = sd_y)
  data.frame(x1 = x1, x2 = x2, x3 = x3, x4 = x4, y = y)
}

# Hard synthetic: latent strata violating homogeneity (HOM).
# Two latent groups with different theta_ref values; the model
# assumes a single theta_ref. This breaks (HOM) of Block 2.
# Coefficients are anti-symmetric by legacy construction; no longer
# required by the model after the (C2)/(C3) implementation fix.
make_hom_violated <- function(n = 600, theta_a = 1.0, theta_b = 3.0,
                              prop_a = 0.5, sd_y = 0.4, seed = 6L) {
  set.seed(seed)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  group <- rbinom(n, 1, prop_a)
  theta_i <- ifelse(group == 1, theta_a, theta_b)
  eta <- theta_i + 0.3 * x1 - 0.3 * x2
  y <- eta + rnorm(n, sd = sd_y)
  data.frame(x1 = x1, x2 = x2, y = y, group = group)
}

# Hard synthetic: near-collinearity between two basis terms.
# Coefficients are anti-symmetric by legacy construction; no longer
# required by the model after the (C2)/(C3) implementation fix.
make_near_collinear <- function(n = 400, rho = 0.99, sd_y = 0.4,
                                seed = 7L) {
  set.seed(seed)
  x1 <- rnorm(n)
  x2 <- rho * x1 + sqrt(1 - rho^2) * rnorm(n)
  eta <- 1 + 0.4 * x1 - 0.4 * x2
  y <- eta + rnorm(n, sd = sd_y)
  data.frame(x1 = x1, x2 = x2, y = y)
}

# Hard synthetic: heavy-tailed errors (Cauchy contamination).
# Coefficients are anti-symmetric by legacy construction; no longer
# required by the model after the (C2)/(C3) implementation fix.
make_heavy_tailed <- function(n = 400, contamination = 0.1, sd_y = 0.4,
                              seed = 8L) {
  set.seed(seed)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  eta <- 1 + 0.4 * x1 - 0.4 * x2
  is_outlier <- rbinom(n, 1, contamination)
  err <- ifelse(is_outlier == 1,
                rcauchy(n, scale = 5),
                rnorm(n, sd = sd_y))
  y <- eta + err
  data.frame(x1 = x1, x2 = x2, y = y, is_outlier = is_outlier)
}

# Hard synthetic: zero-inflated counts.
# Coefficients are anti-symmetric by legacy construction; no longer
# required by the model after the (C2)/(C3) implementation fix.
make_zero_inflated <- function(n = 500, pi_zero = 0.3, seed = 9L) {
  set.seed(seed)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  eta <- 1 + 0.3 * x1 - 0.3 * x2
  is_struct_zero <- rbinom(n, 1, pi_zero)
  y_count <- rpois(n, exp(eta))
  y <- ifelse(is_struct_zero == 1, 0L, y_count)
  data.frame(x1 = x1, x2 = x2, y = y)
}

# ------------------------------------------------------------
# Posterior-predictive sampling for gdpar fits on new data
# ------------------------------------------------------------

# Generate posterior-predictive samples y_pred[s, i] for newdata, with
# observation-level noise drawn from the appropriate family. The
# returned matrix is comparable in semantics to rstanarm::posterior_predict
# and brms::posterior_predict, enabling fair CRPS / ELPD comparisons in
# benchmark_competitors.R.
#
# The Stan template ships with a generated_quantities y_pred for the
# training data; this helper reconstructs the equivalent draws on
# arbitrary newdata using the predictive linear-predictor draws and
# the family-specific dispersion / probability transform.
#
# fit: a gdpar_fit object.
# newdata: data.frame on which to predict.
# family: character scalar in c("gaussian", "neg_binomial_2",
#         "poisson", "bernoulli").
# Returns: numeric matrix of shape (n_draws, n_obs) with predictive
#   samples on the response scale.
bench_posterior_predict <- function(fit, newdata, family) {
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("Package 'posterior' is required for bench_posterior_predict().")
  }
  eta_draws <- predict(fit, newdata = newdata, summary = "draws")
  n_draws <- nrow(eta_draws)
  n_obs   <- ncol(eta_draws)
  if (family == "gaussian") {
    sigma_y_draws <- as.numeric(posterior::as_draws_matrix(
      posterior::subset_draws(fit$fit$draws(), variable = "sigma_y")))
    y_pred <- eta_draws + matrix(rnorm(n_draws * n_obs, 0, sigma_y_draws),
                                 n_draws, n_obs)
  } else if (family == "neg_binomial_2") {
    phi_draws <- as.numeric(posterior::as_draws_matrix(
      posterior::subset_draws(fit$fit$draws(), variable = "phi")))
    mu <- exp(eta_draws)
    y_pred <- matrix(rnbinom(n_draws * n_obs,
                             size = phi_draws,
                             mu = as.vector(mu)),
                     n_draws, n_obs)
  } else if (family == "poisson") {
    lambda <- exp(eta_draws)
    y_pred <- matrix(rpois(n_draws * n_obs, lambda = as.vector(lambda)),
                     n_draws, n_obs)
  } else if (family == "bernoulli") {
    p <- plogis(eta_draws)
    y_pred <- matrix(rbinom(n_draws * n_obs, size = 1,
                            prob = as.vector(p)),
                     n_draws, n_obs)
  } else {
    stop(sprintf("Unsupported family in bench_posterior_predict: %s",
                  family))
  }
  y_pred
}

# ------------------------------------------------------------
# Predictive performance metrics
# ------------------------------------------------------------

# Continuous Ranked Probability Score for an ensemble of posterior
# predictive samples. Generic implementation: works for any
# (continuous or discrete) response when an ensemble of draws is
# available.
#
# y: numeric vector of length n_obs (observed outcomes)
# y_pred: numeric matrix of shape n_draws x n_obs (posterior
#   predictive samples)
# Returns: numeric vector of length n_obs (CRPS per observation).
crps_ensemble <- function(y, y_pred) {
  if (!is.matrix(y_pred)) {
    stop("y_pred must be a numeric matrix of shape n_draws x n_obs.")
  }
  n_obs <- length(y)
  if (ncol(y_pred) != n_obs) {
    stop("ncol(y_pred) must equal length(y).")
  }
  if (requireNamespace("scoringRules", quietly = TRUE)) {
    scoringRules::crps_sample(y = y, dat = t(y_pred))
  } else {
    out <- numeric(n_obs)
    for (i in seq_len(n_obs)) {
      x <- y_pred[, i]
      n_draws <- length(x)
      term1 <- mean(abs(x - y[i]))
      term2 <- 0
      for (a in seq_len(n_draws)) {
        term2 <- term2 + sum(abs(x[a] - x))
      }
      term2 <- term2 / (2 * n_draws^2)
      out[i] <- term1 - term2
    }
    out
  }
}

# Brier score for binary outcomes.
brier_score <- function(y, p_pred) {
  if (length(y) != length(p_pred)) {
    stop("y and p_pred must have the same length.")
  }
  mean((p_pred - y)^2)
}

# Logarithmic loss for binary outcomes.
log_loss <- function(y, p_pred, eps = 1e-15) {
  p <- pmin(pmax(p_pred, eps), 1 - eps)
  -mean(y * log(p) + (1 - y) * log(1 - p))
}

# Area under the receiver operating characteristic curve.
auc_roc <- function(y, p_pred) {
  if (length(unique(y)) < 2) {
    return(NA_real_)
  }
  ord <- order(p_pred, decreasing = TRUE)
  y_ord <- y[ord]
  n_pos <- sum(y == 1)
  n_neg <- sum(y == 0)
  if (n_pos == 0 || n_neg == 0) return(NA_real_)
  cum_tp <- cumsum(y_ord == 1) / n_pos
  cum_fp <- cumsum(y_ord == 0) / n_neg
  sum(diff(c(0, cum_fp)) * (cum_tp + c(0, cum_tp[-length(cum_tp)])) / 2)
}

# Calibration table for binary outcomes via decile bins.
calibration_table <- function(y, p_pred, n_bins = 10) {
  bin_id <- cut(p_pred,
                breaks = seq(0, 1, length.out = n_bins + 1),
                include.lowest = TRUE, labels = FALSE)
  out <- data.frame(
    bin = seq_len(n_bins),
    n_obs = tapply(y, bin_id, length)[as.character(seq_len(n_bins))],
    mean_predicted = tapply(p_pred, bin_id, mean)[as.character(seq_len(n_bins))],
    mean_observed = tapply(y, bin_id, mean)[as.character(seq_len(n_bins))]
  )
  out
}

# Probability Integral Transform values for an ensemble of posterior
# predictive samples. For continuous outcomes only; returns the
# empirical quantile of y_i within the predictive distribution.
pit_values <- function(y, y_pred) {
  if (!is.matrix(y_pred)) {
    stop("y_pred must be a numeric matrix of shape n_draws x n_obs.")
  }
  vapply(seq_along(y), function(i) {
    mean(y_pred[, i] <= y[i])
  }, numeric(1L))
}

# RMSE.
rmse <- function(y, y_pred_mean) {
  sqrt(mean((y - y_pred_mean)^2))
}

# Mean absolute error.
mae <- function(y, y_pred_mean) {
  mean(abs(y - y_pred_mean))
}

# ELPD via PSIS-LOO from a log-likelihood matrix (n_draws x n_obs).
elpd_loo_psis <- function(log_lik) {
  if (!bench_require("loo")) return(NA_real_)
  loo_obj <- loo::loo(log_lik)
  list(
    elpd_loo = loo_obj$estimates["elpd_loo", "Estimate"],
    elpd_loo_se = loo_obj$estimates["elpd_loo", "SE"],
    p_loo = loo_obj$estimates["p_loo", "Estimate"],
    pareto_k_summary = summary(loo_obj$diagnostics$pareto_k)
  )
}

# ------------------------------------------------------------
# Coverage of credible / confidence intervals
# ------------------------------------------------------------

# Empirical coverage of nominal-level intervals.
# y_lower, y_upper: numeric vectors of length n_obs with the
#   nominal interval bounds.
# y: observed outcomes.
# Returns the empirical proportion of observations inside the interval.
empirical_coverage <- function(y, y_lower, y_upper) {
  mean(y >= y_lower & y <= y_upper)
}

# ------------------------------------------------------------
# Common compare-frame builder
# ------------------------------------------------------------

# Build a tidy comparison data.frame with one row per (method, metric).
build_compare_df <- function(method, metrics_named_list) {
  data.frame(
    method = method,
    metric = names(metrics_named_list),
    value = unlist(metrics_named_list),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

# ------------------------------------------------------------
# Reporting
# ------------------------------------------------------------

bench_table <- function(df, digits = 4) {
  if (nrow(df) == 0L) {
    cat("(empty table)\n")
    return(invisible(NULL))
  }
  num_cols <- vapply(df, is.numeric, logical(1L))
  for (j in which(num_cols)) {
    df[[j]] <- format(df[[j]], digits = digits, nsmall = digits,
                      scientific = FALSE)
  }
  print(df, row.names = FALSE)
  invisible(NULL)
}

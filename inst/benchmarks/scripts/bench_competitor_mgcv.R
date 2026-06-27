## mgcv wrapper for the Block 6 multi adversarial bench.
##
## Implements PSIS-LOO via the Bayesian interpretation of REML
## (Wood 2017, "Generalized Additive Models in R", section 6.10):
## the covariance Vp = vcov(gam_obj) defines a multivariate Gaussian
## approximation of the posterior over coefficients, conditional on
## the smoothing parameters. We sample beta_s ~ MVN(coef(gam_obj),
## Vp) for s = 1..S, recompute the linear predictor X_pred %*% beta_s,
## map through the inverse link, evaluate log p(y_ik | eta_s_i_k,
## sigma_s) per family, sum over k for the subject aggregation, and
## pass the resulting S x n matrix to loo::loo with chain_id = 1.
##
## Wood reports that this Bayesian interpretation reproduces the
## frequentist confidence intervals of gam under mild regularity, so
## PSIS-LOO on these draws is comparable in scale to fully Bayesian
## fits (gdpar, brms, rstanarm).

bench_mgcv_fit <- function(data, family_spec, scenario_meta, ...) {
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("Package 'mgcv' is required.")
  }
  p <- as.integer(scenario_meta$p)
  fams <- family_spec$mgcv_family
  formulas <- scenario_meta$mgcv_formula_per_coord
  fits <- vector("list", p)
  t0 <- Sys.time()
  for (k in seq_len(p)) {
    fits[[k]] <- mgcv::gam(formulas[[k]], data = data,
                            family = fams[[k]], method = "REML")
  }
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  list(fits = fits, wall_seconds = wall, p = p,
       family_spec = family_spec, formulas = formulas)
}

bench_mgcv_predict <- function(fit_obj, newdata) {
  p <- fit_obj$p
  fam_name <- fit_obj$family_spec$brms_family
  mean_mat <- matrix(0, nrow = nrow(newdata), ncol = p)
  sd_mat   <- matrix(0, nrow = nrow(newdata), ncol = p)
  prob_mat <- if (fam_name == "bernoulli")
    matrix(0, nrow = nrow(newdata), ncol = p) else NULL
  for (k in seq_len(p)) {
    fit_k <- fit_obj$fits[[k]]
    pred_lp <- stats::predict(fit_k, newdata = newdata, type = "link",
                                se.fit = TRUE)
    inv_link <- fit_k$family$linkinv
    mean_mat[, k] <- inv_link(pred_lp$fit)
    sd_mat[, k] <- pred_lp$se.fit
    if (fam_name == "bernoulli") {
      prob_mat[, k] <- inv_link(pred_lp$fit)
    }
  }
  list(mean = if (p == 1L) mean_mat[, 1L] else mean_mat,
       sd   = if (p == 1L) sd_mat[, 1L]   else sd_mat,
       prob = if (is.null(prob_mat)) NULL
              else if (p == 1L) prob_mat[, 1L] else prob_mat,
       draws = NULL)
}

bench_mgcv_loo <- function(fit_obj, S = 4000L, seed = 42L) {
  if (!requireNamespace("loo", quietly = TRUE) ||
      !requireNamespace("mvnfast", quietly = TRUE)) {
    return(NULL)
  }
  p <- fit_obj$p
  n <- nrow(fit_obj$fits[[1L]]$model)
  log_lik_subject <- matrix(0, nrow = S, ncol = n)
  fam_name <- fit_obj$family_spec$brms_family
  for (k in seq_len(p)) {
    gam_k <- fit_obj$fits[[k]]
    Vp <- gam_k$Vp
    bhat <- stats::coef(gam_k)
    set.seed(seed + k)
    beta_draws <- mvnfast::rmvn(n = S, mu = bhat, sigma = Vp)
    X_train <- stats::predict(gam_k, type = "lpmatrix")
    eta_draws <- beta_draws %*% t(X_train)
    inv_link <- gam_k$family$linkinv
    y_k <- gam_k$model[[1L]]
    if (fam_name == "gaussian") {
      sigma_hat <- sqrt(gam_k$sig2)
      ll_k <- t(apply(eta_draws, 1L, function(eta_s) {
        stats::dnorm(y_k, mean = inv_link(eta_s),
                      sd = sigma_hat, log = TRUE)
      }))
    } else if (fam_name == "bernoulli") {
      ll_k <- t(apply(eta_draws, 1L, function(eta_s) {
        prob_s <- inv_link(eta_s)
        prob_s <- pmin(pmax(prob_s, 1e-15), 1 - 1e-15)
        y_k * log(prob_s) + (1 - y_k) * log(1 - prob_s)
      }))
    } else if (fam_name == "poisson") {
      ll_k <- t(apply(eta_draws, 1L, function(eta_s) {
        stats::dpois(y_k, lambda = inv_link(eta_s), log = TRUE)
      }))
    } else {
      return(NULL)
    }
    log_lik_subject <- log_lik_subject + ll_k
  }
  r_eff <- loo::relative_eff(exp(log_lik_subject),
                              chain_id = rep(1L, S))
  tryCatch(
    loo::loo(log_lik_subject, r_eff = r_eff),
    error = function(e) {
      warning("mgcv loo failed: ", conditionMessage(e))
      NULL
    }
  )
}

bench_mgcv_metrics <- function(fit_obj, test, scenario_meta) {
  preds <- bench_mgcv_predict(fit_obj, test)
  p <- fit_obj$p
  fam_name <- fit_obj$family_spec$brms_family
  classification <- isTRUE(scenario_meta$classification_metrics) ||
                     identical(fam_name, "bernoulli")
  rows <- list()
  for (k in seq_len(p)) {
    y_col <- paste0("y", k)
    y_true <- test[[y_col]]
    mean_k <- if (p == 1L) preds$mean else preds$mean[, k]
    sd_k   <- if (p == 1L) preds$sd else preds$sd[, k]
    prob_k <- if (!is.null(preds$prob)) {
      if (p == 1L) preds$prob else preds$prob[, k]
    } else NULL
    if (classification) {
      y_int <- as.integer(y_true)
      rows[[length(rows) + 1L]] <- bench_record(
        scenario_meta$scenario_id, p, "mgcv", fam_name,
        "auc", as.character(k), bench_auc(y_int, prob_k))
      rows[[length(rows) + 1L]] <- bench_record(
        scenario_meta$scenario_id, p, "mgcv", fam_name,
        "log_loss", as.character(k), bench_logloss(y_int, prob_k))
      rows[[length(rows) + 1L]] <- bench_record(
        scenario_meta$scenario_id, p, "mgcv", fam_name,
        "brier", as.character(k), bench_brier(y_int, prob_k))
    } else {
      rows[[length(rows) + 1L]] <- bench_record(
        scenario_meta$scenario_id, p, "mgcv", fam_name,
        "rmse", as.character(k), bench_rmse(y_true, mean_k))
      rows[[length(rows) + 1L]] <- bench_record(
        scenario_meta$scenario_id, p, "mgcv", fam_name,
        "mae", as.character(k), bench_mae(y_true, mean_k))
      sigma_hat <- sqrt(fit_obj$fits[[k]]$sig2)
      rows[[length(rows) + 1L]] <- bench_record(
        scenario_meta$scenario_id, p, "mgcv", fam_name,
        "crps", as.character(k),
        bench_crps_gaussian(y_true, mean_k, sigma_hat))
    }
  }
  lo <- bench_mgcv_loo(fit_obj)
  if (!is.null(lo)) {
    pk <- bench_pareto_k_summary(lo)
    rows[[length(rows) + 1L]] <- bench_record(
      scenario_meta$scenario_id, p, "mgcv", fam_name,
      "elpd_loo", "agg",
      lo$estimates["elpd_loo", "Estimate"],
      se_value = lo$estimates["elpd_loo", "SE"],
      convergence_flag = bench_convergence_flag_from_loo(lo),
      max_pareto_k = pk$max_pareto_k,
      n_pareto_k_above_07 = pk$n_above_07,
      extra = "wood_2017_laplace_approx")
  } else {
    rows[[length(rows) + 1L]] <- bench_record(
      scenario_meta$scenario_id, p, "mgcv", fam_name,
      "elpd_loo", "agg", NA_real_,
      convergence_flag = "loo_failed")
  }
  aic_total <- sum(vapply(fit_obj$fits, AIC, numeric(1L)))
  rows[[length(rows) + 1L]] <- bench_record(
    scenario_meta$scenario_id, p, "mgcv", fam_name,
    "aic", "agg", aic_total,
    extra = "mgcv_builtin_cross_check")
  rows[[length(rows) + 1L]] <- bench_record(
    scenario_meta$scenario_id, p, "mgcv", fam_name,
    "wall", "fit", fit_obj$wall_seconds,
    wall_seconds = fit_obj$wall_seconds)
  do.call(rbind, rows)
}

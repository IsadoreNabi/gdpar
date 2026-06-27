## INLA wrapper for the Block 6 multi adversarial bench.
##
## Per-coordinate INLA fits (one inla() call per k). PSIS-LOO via
## inla.posterior.sample(S, fit) + reconstructed log_lik + loo::loo()
## as the principal ELPD metric (comparable to gdpar/brms/rstanarm/mgcv).
## CPO sum is persisted as a secondary cross-check; if PSIS-LOO and
## sum(log(cpo)) diverge by more than 2 * se_elpd the convergence_flag
## is set to "inla_cpo_loo_divergent".

bench_inla_fit <- function(data, family_spec, scenario_meta, ...) {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop("Package 'INLA' is required.")
  }
  p <- as.integer(scenario_meta$p)
  fams <- family_spec$inla_family
  formulas <- scenario_meta$inla_formula_per_coord
  fits <- vector("list", p)
  cpo_list <- vector("list", p)
  t0 <- Sys.time()
  for (k in seq_len(p)) {
    fits[[k]] <- INLA::inla(
      formulas[[k]], family = fams[k], data = data,
      control.compute = list(cpo = TRUE, config = TRUE,
                              return.marginals = TRUE),
      control.predictor = list(compute = TRUE, link = 1L),
      verbose = FALSE
    )
    cpo_list[[k]] <- fits[[k]]$cpo$cpo
  }
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  list(fits = fits, cpo = cpo_list, wall_seconds = wall, p = p,
       family_spec = family_spec, formulas = formulas)
}

bench_inla_predict <- function(fit_obj, newdata) {
  p <- fit_obj$p
  fam_name <- fit_obj$family_spec$brms_family
  mean_mat <- matrix(0, nrow = nrow(newdata), ncol = p)
  sd_mat   <- matrix(0, nrow = nrow(newdata), ncol = p)
  prob_mat <- if (fam_name == "bernoulli")
    matrix(0, nrow = nrow(newdata), ncol = p) else NULL
  for (k in seq_len(p)) {
    fit_k <- fit_obj$fits[[k]]
    formula_k <- fit_obj$formulas[[k]]
    fam_k <- fit_obj$family_spec$inla_family[k]
    train_data <- fit_k$.args$data
    augmented <- newdata
    augmented[[paste0("y", k)]] <- NA_real_
    keep_cols <- setdiff(intersect(colnames(train_data),
                                       colnames(augmented)),
                          character(0))
    augmented <- augmented[, keep_cols, drop = FALSE]
    for (cn in setdiff(colnames(train_data), keep_cols)) {
      augmented[[cn]] <- NA
    }
    augmented <- augmented[, colnames(train_data), drop = FALSE]
    full_data <- rbind(train_data, augmented)
    refit <- tryCatch(
      INLA::inla(formula_k, family = fam_k, data = full_data,
                  control.predictor = list(compute = TRUE, link = 1L),
                  control.compute = list(config = TRUE),
                  verbose = FALSE),
      error = function(e) NULL
    )
    if (is.null(refit)) {
      mean_mat[, k] <- NA_real_; sd_mat[, k] <- NA_real_
      if (!is.null(prob_mat)) prob_mat[, k] <- NA_real_
      next
    }
    n_train <- nrow(train_data)
    pred_rows <- (n_train + 1L):(n_train + nrow(newdata))
    mean_mat[, k] <- refit$summary.fitted.values$mean[pred_rows]
    sd_eta_k <- refit$summary.fitted.values$sd[pred_rows]
    if (fam_name == "gaussian") {
      prec_name <- "Precision for the Gaussian observations"
      tau_y <- if (!is.null(refit$summary.hyperpar) &&
                   prec_name %in% rownames(refit$summary.hyperpar))
                 refit$summary.hyperpar[prec_name, "mean"]
               else NA_real_
      sigma_y_k <- if (is.finite(tau_y) && tau_y > 0)
                     sqrt(1 / as.numeric(tau_y)) else 0
      sd_mat[, k] <- sqrt(sd_eta_k^2 + sigma_y_k^2)
    } else {
      sd_mat[, k] <- sd_eta_k
    }
    if (fam_name == "bernoulli") {
      prob_mat[, k] <- refit$summary.fitted.values$mean[pred_rows]
    }
  }
  list(mean = if (p == 1L) mean_mat[, 1L] else mean_mat,
       sd   = if (p == 1L) sd_mat[, 1L]   else sd_mat,
       prob = if (is.null(prob_mat)) NULL
              else if (p == 1L) prob_mat[, 1L] else prob_mat,
       draws = NULL)
}

bench_inla_loo <- function(fit_obj, S = 2000L) {
  if (!requireNamespace("INLA", quietly = TRUE) ||
      !requireNamespace("loo", quietly = TRUE)) return(NULL)
  p <- fit_obj$p
  fam_name <- fit_obj$family_spec$brms_family
  fits <- fit_obj$fits
  n <- length(fits[[1L]]$cpo$cpo)
  log_lik_subject <- matrix(0, nrow = S, ncol = n)
  for (k in seq_len(p)) {
    samples <- tryCatch(
      INLA::inla.posterior.sample(S, fits[[k]]),
      error = function(e) NULL
    )
    if (is.null(samples)) return(NULL)
    y_k <- fits[[k]]$.args$data[[paste0("y", k)]]
    if (is.null(y_k)) y_k <- fits[[k]]$.args$data[["y"]]
    pred_pattern <- "^Predictor:"
    pred_idx <- grep(pred_pattern,
                       rownames(samples[[1L]]$latent))
    ll_k <- matrix(0, nrow = S, ncol = n)
    for (s in seq_len(S)) {
      eta_s <- samples[[s]]$latent[pred_idx, 1L][seq_len(n)]
      if (fam_name == "gaussian") {
        prec_name <- "Precision for the Gaussian observations"
        tau_s <- samples[[s]]$hyperpar[prec_name]
        sigma_s <- if (is.null(tau_s) || is.na(tau_s)) 1
                    else sqrt(1 / as.numeric(tau_s))
        ll_k[s, ] <- stats::dnorm(y_k, mean = eta_s, sd = sigma_s,
                                    log = TRUE)
      } else if (fam_name == "bernoulli") {
        prob_s <- 1 / (1 + exp(-eta_s))
        prob_s <- pmin(pmax(prob_s, 1e-15), 1 - 1e-15)
        ll_k[s, ] <- y_k * log(prob_s) + (1 - y_k) * log(1 - prob_s)
      } else if (fam_name == "poisson") {
        ll_k[s, ] <- stats::dpois(y_k, lambda = exp(eta_s),
                                    log = TRUE)
      } else {
        return(NULL)
      }
    }
    log_lik_subject <- log_lik_subject + ll_k
  }
  r_eff <- loo::relative_eff(exp(log_lik_subject),
                              chain_id = rep(1L, S))
  tryCatch(loo::loo(log_lik_subject, r_eff = r_eff),
            error = function(e) {
              warning("inla loo failed: ", conditionMessage(e))
              NULL
            })
}

bench_inla_metrics <- function(fit_obj, test, scenario_meta) {
  preds <- bench_inla_predict(fit_obj, test)
  p <- fit_obj$p
  fam_name <- fit_obj$family_spec$brms_family
  classification <- isTRUE(scenario_meta$classification_metrics) ||
                     identical(fam_name, "bernoulli")
  rows <- list()
  for (k in seq_len(p)) {
    y_col <- paste0("y", k); y_true <- test[[y_col]]
    mean_k <- if (p == 1L) preds$mean else preds$mean[, k]
    sd_k   <- if (p == 1L) preds$sd else preds$sd[, k]
    prob_k <- if (!is.null(preds$prob)) {
      if (p == 1L) preds$prob else preds$prob[, k]
    } else NULL
    if (classification) {
      y_int <- as.integer(y_true)
      rows[[length(rows) + 1L]] <- bench_record(
        scenario_meta$scenario_id, p, "inla", fam_name,
        "auc", as.character(k), bench_auc(y_int, prob_k))
      rows[[length(rows) + 1L]] <- bench_record(
        scenario_meta$scenario_id, p, "inla", fam_name,
        "log_loss", as.character(k), bench_logloss(y_int, prob_k))
      rows[[length(rows) + 1L]] <- bench_record(
        scenario_meta$scenario_id, p, "inla", fam_name,
        "brier", as.character(k), bench_brier(y_int, prob_k))
    } else {
      rows[[length(rows) + 1L]] <- bench_record(
        scenario_meta$scenario_id, p, "inla", fam_name,
        "rmse", as.character(k), bench_rmse(y_true, mean_k))
      rows[[length(rows) + 1L]] <- bench_record(
        scenario_meta$scenario_id, p, "inla", fam_name,
        "mae", as.character(k), bench_mae(y_true, mean_k))
      rows[[length(rows) + 1L]] <- bench_record(
        scenario_meta$scenario_id, p, "inla", fam_name,
        "crps", as.character(k),
        bench_crps_gaussian(y_true, mean_k, sd_k))
    }
  }
  lo <- bench_inla_loo(fit_obj)
  cpo_total <- sum(vapply(fit_obj$cpo, function(cv)
                            sum(log(pmax(cv, 1e-15))),
                          numeric(1L)))
  conv_flag <- "ok"
  if (!is.null(lo)) {
    elpd <- lo$estimates["elpd_loo", "Estimate"]
    se   <- lo$estimates["elpd_loo", "SE"]
    if (is.finite(elpd) && is.finite(se) &&
        abs(cpo_total - elpd) > 2 * se) {
      conv_flag <- "inla_cpo_loo_divergent"
    } else {
      conv_flag <- bench_convergence_flag_from_loo(lo)
    }
    pk <- bench_pareto_k_summary(lo)
    rows[[length(rows) + 1L]] <- bench_record(
      scenario_meta$scenario_id, p, "inla", fam_name,
      "elpd_loo", "agg", elpd, se_value = se,
      convergence_flag = conv_flag,
      max_pareto_k = pk$max_pareto_k,
      n_pareto_k_above_07 = pk$n_above_07)
  } else {
    rows[[length(rows) + 1L]] <- bench_record(
      scenario_meta$scenario_id, p, "inla", fam_name,
      "elpd_loo", "agg", NA_real_,
      convergence_flag = "loo_failed")
  }
  rows[[length(rows) + 1L]] <- bench_record(
    scenario_meta$scenario_id, p, "inla", fam_name,
    "elpd_cpo", "agg", cpo_total,
    extra = "inla_builtin_cross_check")
  rows[[length(rows) + 1L]] <- bench_record(
    scenario_meta$scenario_id, p, "inla", fam_name,
    "wall", "fit", fit_obj$wall_seconds,
    wall_seconds = fit_obj$wall_seconds)
  do.call(rbind, rows)
}

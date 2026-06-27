## brms wrapper for the Block 6 multi adversarial bench.
##
## Policy: set_rescor(FALSE) in all scenarios for consistency with
## gdpar coord-wise factorization. S8 additionally runs a second fit
## with set_rescor(TRUE) under the method label "brms_rescor" to
## verify recovery of the true residual correlation rho.

bench_brms_fit <- function(data, family_spec, scenario_meta,
                            rescor = FALSE, ...) {
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("Package 'brms' is required.")
  }
  iter <- bench_default_iter()
  p <- as.integer(scenario_meta$p)
  rhs <- scenario_meta$brms_rhs
  fam_name <- family_spec$brms_family
  fam <- switch(fam_name,
                  gaussian  = brms::brmsfamily("gaussian"),
                  bernoulli = brms::brmsfamily("bernoulli"),
                  poisson   = brms::brmsfamily("poisson"),
                  negbinomial = brms::brmsfamily("negbinomial"),
                  stop(sprintf("Unsupported brms family: %s", fam_name)))
  formula_obj <- if (p == 1L) {
    brms::bf(stats::as.formula(sprintf("y1 ~ %s", rhs)))
  } else {
    y_terms <- paste0("y", seq_len(p))
    response <- sprintf("mvbind(%s)", paste(y_terms, collapse = ", "))
    f <- brms::bf(stats::as.formula(sprintf("%s ~ %s", response, rhs)))
    if (!rescor) f <- f + brms::set_rescor(FALSE)
    else         f <- f + brms::set_rescor(TRUE)
    f
  }
  timed <- bench_time_fit(
    suppressMessages(suppressWarnings(
      brms::brm(formula_obj, data = data, family = fam,
                  chains = iter$chains,
                  warmup = iter$iter_warmup,
                  iter = iter$iter_warmup + iter$iter_sampling,
                  seed = iter$seed, refresh = iter$refresh,
                  backend = "rstan", silent = 2L,
                  cores = iter$chains)
    ))
  )
  list(fit = timed$value, wall_seconds = timed$wall_seconds,
       p = p, family_spec = family_spec, rescor = rescor)
}

bench_brms_predict <- function(fit_obj, newdata) {
  fit <- fit_obj$fit
  p <- fit_obj$p
  fam_name <- fit_obj$family_spec$brms_family
  if (fam_name == "bernoulli") {
    draws_prob <- brms::posterior_epred(fit, newdata = newdata)
  } else {
    draws_prob <- NULL
  }
  draws_mean <- brms::posterior_epred(fit, newdata = newdata)
  draws_pred <- brms::posterior_predict(fit, newdata = newdata)
  if (p == 1L) {
    mean_pred <- colMeans(draws_mean)
    sd_pred   <- apply(draws_pred, 2L, stats::sd)
    return(list(mean = mean_pred, sd = sd_pred,
                 prob = if (!is.null(draws_prob)) colMeans(draws_prob)
                         else NULL,
                 draws = draws_pred))
  }
  if (length(dim(draws_mean)) == 3L) {
    mean_mat <- apply(draws_mean, c(2L, 3L), mean)
    sd_mat   <- apply(draws_pred, c(2L, 3L), stats::sd)
    prob_mat <- if (!is.null(draws_prob))
      apply(draws_prob, c(2L, 3L), mean) else NULL
    return(list(mean = mean_mat, sd = sd_mat, prob = prob_mat,
                 draws = draws_pred))
  }
  list(mean = colMeans(draws_mean),
       sd = apply(draws_pred, 2L, stats::sd),
       prob = if (!is.null(draws_prob)) colMeans(draws_prob) else NULL,
       draws = draws_pred)
}

bench_brms_loo <- function(fit_obj) {
  tryCatch(
    suppressWarnings(brms::loo(fit_obj$fit)),
    error = function(e) {
      warning("brms::loo failed: ", conditionMessage(e))
      NULL
    }
  )
}

bench_brms_recovered_rescor <- function(fit_obj) {
  if (!isTRUE(fit_obj$rescor)) return(NA_real_)
  tryCatch({
    summ <- summary(fit_obj$fit)
    rescor_row <- summ$rescor_pars
    if (is.null(rescor_row) || nrow(rescor_row) == 0L) return(NA_real_)
    as.numeric(rescor_row$Estimate[1L])
  }, error = function(e) NA_real_)
}

bench_brms_metrics <- function(fit_obj, test, scenario_meta,
                                 method_label = "brms") {
  preds <- bench_brms_predict(fit_obj, test)
  p <- fit_obj$p
  fam_name <- fit_obj$family_spec$brms_family
  classification <- isTRUE(scenario_meta$classification_metrics) ||
                     identical(fam_name, "bernoulli")
  rows <- list()
  for (k in seq_len(p)) {
    y_col <- if (p == 1L) "y1" else paste0("y", k)
    y_true <- test[[y_col]]
    mean_k <- if (p == 1L) preds$mean else preds$mean[, k]
    sd_k   <- if (p == 1L) preds$sd else preds$sd[, k]
    prob_k <- if (!is.null(preds$prob)) {
      if (p == 1L) preds$prob else preds$prob[, k]
    } else NULL
    draws_k <- if (length(dim(preds$draws)) == 3L) preds$draws[, , k]
                else preds$draws
    if (classification) {
      y_int <- as.integer(y_true)
      rows[[length(rows) + 1L]] <- bench_record(
        scenario_meta$scenario_id, p, method_label, fam_name,
        "auc", as.character(k), bench_auc(y_int, prob_k))
      rows[[length(rows) + 1L]] <- bench_record(
        scenario_meta$scenario_id, p, method_label, fam_name,
        "log_loss", as.character(k), bench_logloss(y_int, prob_k))
      rows[[length(rows) + 1L]] <- bench_record(
        scenario_meta$scenario_id, p, method_label, fam_name,
        "brier", as.character(k), bench_brier(y_int, prob_k))
    } else {
      rows[[length(rows) + 1L]] <- bench_record(
        scenario_meta$scenario_id, p, method_label, fam_name,
        "rmse", as.character(k), bench_rmse(y_true, mean_k))
      rows[[length(rows) + 1L]] <- bench_record(
        scenario_meta$scenario_id, p, method_label, fam_name,
        "mae", as.character(k), bench_mae(y_true, mean_k))
      rows[[length(rows) + 1L]] <- bench_record(
        scenario_meta$scenario_id, p, method_label, fam_name,
        "crps", as.character(k),
        if (is.matrix(draws_k)) bench_crps_sample(y_true, draws_k)
                                 else bench_crps_gaussian(y_true, mean_k, sd_k))
    }
  }
  lo <- bench_brms_loo(fit_obj)
  if (!is.null(lo)) {
    pk <- bench_pareto_k_summary(lo)
    rows[[length(rows) + 1L]] <- bench_record(
      scenario_meta$scenario_id, p, method_label, fam_name,
      "elpd_loo", "agg",
      lo$estimates["elpd_loo", "Estimate"],
      se_value = lo$estimates["elpd_loo", "SE"],
      convergence_flag = bench_convergence_flag_from_loo(lo),
      max_pareto_k = pk$max_pareto_k,
      n_pareto_k_above_07 = pk$n_above_07,
      recovered_rescor = bench_brms_recovered_rescor(fit_obj))
  } else {
    rows[[length(rows) + 1L]] <- bench_record(
      scenario_meta$scenario_id, p, method_label, fam_name,
      "elpd_loo", "agg", NA_real_,
      convergence_flag = "loo_failed",
      recovered_rescor = bench_brms_recovered_rescor(fit_obj))
  }
  rows[[length(rows) + 1L]] <- bench_record(
    scenario_meta$scenario_id, p, method_label, fam_name,
    "wall", "fit", fit_obj$wall_seconds,
    wall_seconds = fit_obj$wall_seconds,
    recovered_rescor = bench_brms_recovered_rescor(fit_obj))
  do.call(rbind, rows)
}

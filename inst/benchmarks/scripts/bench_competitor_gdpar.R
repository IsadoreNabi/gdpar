## gdpar wrapper for the Block 6 multi adversarial bench.
##
## Exposes:
##   bench_gdpar_fit(data, family_spec, scenario_meta) -> list
##   bench_gdpar_predict(fit_obj, newdata) -> list with mean / sd /
##     prob (depending on family) and posterior draws
##   bench_gdpar_loo(fit_obj) -> psis_loo
##   bench_gdpar_metrics(fit_obj, test, scenario_meta) -> data.frame
##
## All wall_seconds are sampling time only (compile cached across
## scenarios via the cmdstanr default).

bench_gdpar_fit <- function(data, family_spec, scenario_meta,
                              ...) {
  iter <- bench_default_iter()
  spec <- scenario_meta$gdpar_amm
  fam  <- family_spec$gdpar_family
  rhs  <- paste(scenario_meta$gdpar_rhs, collapse = " + ")
  p <- as.integer(scenario_meta$p)
  if (p == 1L) {
    fm <- stats::as.formula(sprintf("y ~ %s", rhs))
    if (!"y" %in% colnames(data)) {
      data$y <- data$y1
    }
  } else {
    fm <- stats::as.formula(sprintf("y ~ %s", rhs))
    y_cols <- paste0("y", seq_len(p))
    if (!all(y_cols %in% colnames(data))) {
      stop("Expected outcome columns y1..yp in data.")
    }
    data$y <- as.matrix(data[, y_cols, drop = FALSE])
  }
  rigor <- if (!is.null(scenario_meta$rigor)) scenario_meta$rigor
            else "full"
  call_args <- list(
    formula         = fm,
    family          = fam,
    amm             = spec,
    data            = data,
    chains          = iter$chains,
    iter_warmup     = iter$iter_warmup,
    iter_sampling   = iter$iter_sampling,
    seed            = iter$seed,
    refresh         = iter$refresh,
    parametrization = "auto",
    skip_id_check   = FALSE,
    verbose         = FALSE
  )
  if (!is.null(scenario_meta$gdpar_group_formula)) {
    call_args$group <- scenario_meta$gdpar_group_formula
  }
  if (!is.null(scenario_meta$rigor) &&
      scenario_meta$rigor == "fast") {
    call_args$id_check_rigor <- "fast"
  }
  timed <- bench_time_fit(
    suppressMessages(suppressWarnings(
      do.call(gdpar, call_args)
    ))
  )
  list(fit = timed$value, wall_seconds = timed$wall_seconds,
       p = p, family_spec = family_spec)
}

bench_gdpar_predict <- function(fit_obj, newdata) {
  fit <- fit_obj$fit
  p <- fit_obj$p
  family_name <- fit_obj$family_spec$gdpar_family$name
  sigma_y_post <- if (family_name == "gaussian") {
    sy_draws <- fit$fit$draws(variables = "sigma_y",
                                 format = "draws_matrix")
    sy_mat <- unclass(sy_draws)
    apply(sy_mat, 2L, function(z) sqrt(mean(z^2)))
  } else NULL
  if (p == 1L) {
    eta_draws <- stats::predict(fit, newdata = newdata,
                                  type = "linear_predictor",
                                  draws = TRUE)
    mean_pred <- colMeans(eta_draws)
    sd_eta <- apply(eta_draws, 2L, stats::sd)
    sd_pred <- if (family_name == "gaussian") {
      sigma_y_scalar <- if (length(sigma_y_post) >= 1L)
        as.numeric(sigma_y_post[1L]) else 0
      sqrt(sd_eta^2 + sigma_y_scalar^2)
    } else sd_eta
    if (family_name == "gaussian") {
      n_test <- ncol(eta_draws)
      noise <- matrix(stats::rnorm(nrow(eta_draws) * n_test, 0,
                                      sigma_y_scalar),
                       nrow = nrow(eta_draws), ncol = n_test)
      y_pred_draws <- eta_draws + noise
    } else if (family_name == "bernoulli") {
      y_pred_draws <- eta_draws
    } else {
      y_pred_draws <- eta_draws
    }
    prob_pred <- if (family_name == "bernoulli")
      colMeans(1 / (1 + exp(-eta_draws)))
    else NULL
    return(list(mean = mean_pred, sd = sd_pred, prob = prob_pred,
                 draws = y_pred_draws))
  }
  arr <- stats::predict(fit, newdata = newdata,
                         type = "linear_predictor")
  S <- dim(arr)[1L]; n <- dim(arr)[2L]
  mean_mat <- apply(arr, c(2L, 3L), mean)
  sd_eta_mat <- apply(arr, c(2L, 3L), stats::sd)
  sd_mat <- if (family_name == "gaussian") {
    sigma_y_full <- rep(sigma_y_post, each = n)
    sigma_y_mat <- matrix(sigma_y_full, nrow = n, ncol = p, byrow = FALSE)
    sqrt(sd_eta_mat^2 + sigma_y_mat^2)
  } else sd_eta_mat
  if (family_name == "gaussian") {
    y_pred_arr <- array(0, dim = c(S, n, p))
    for (k in seq_len(p)) {
      sy_k <- as.numeric(sigma_y_post[k])
      noise_k <- matrix(stats::rnorm(S * n, 0, sy_k),
                          nrow = S, ncol = n)
      y_pred_arr[, , k] <- arr[, , k] + noise_k
    }
  } else {
    y_pred_arr <- arr
  }
  prob_mat <- if (family_name == "bernoulli") {
    apply(arr, c(2L, 3L), function(z) mean(1 / (1 + exp(-z))))
  } else NULL
  list(mean = mean_mat, sd = sd_mat, prob = prob_mat,
       draws = y_pred_arr)
}

bench_gdpar_loo <- function(fit_obj) {
  tryCatch(
    gdpar_loo(fit_obj$fit, aggregation = "subject"),
    error = function(e) {
      warning("gdpar_loo failed: ", conditionMessage(e))
      NULL
    }
  )
}

bench_gdpar_metrics <- function(fit_obj, test, scenario_meta) {
  preds <- bench_gdpar_predict(fit_obj, test)
  p <- fit_obj$p
  family_name <- fit_obj$family_spec$gdpar_family$name
  classification <- isTRUE(scenario_meta$classification_metrics) ||
                     identical(family_name, "bernoulli")
  rows <- list()
  for (k in seq_len(p)) {
    y_col <- if (p == 1L) "y1" else paste0("y", k)
    if (!y_col %in% colnames(test)) y_col <- "y"
    y_true <- test[[y_col]]
    if (p == 1L) {
      mean_k <- preds$mean; sd_k <- preds$sd
      draws_k <- preds$draws
      prob_k <- preds$prob
    } else {
      mean_k <- preds$mean[, k]; sd_k <- preds$sd[, k]
      draws_k <- preds$draws[, , k]
      prob_k <- if (!is.null(preds$prob)) preds$prob[, k] else NULL
    }
    if (classification) {
      y_int <- as.integer(y_true)
      rows[[length(rows) + 1L]] <- bench_record(
        scenario_meta$scenario_id, p, "gdpar", family_name,
        "auc", as.character(k),
        bench_auc(y_int, prob_k))
      rows[[length(rows) + 1L]] <- bench_record(
        scenario_meta$scenario_id, p, "gdpar", family_name,
        "log_loss", as.character(k),
        bench_logloss(y_int, prob_k))
      rows[[length(rows) + 1L]] <- bench_record(
        scenario_meta$scenario_id, p, "gdpar", family_name,
        "brier", as.character(k),
        bench_brier(y_int, prob_k))
    } else {
      rows[[length(rows) + 1L]] <- bench_record(
        scenario_meta$scenario_id, p, "gdpar", family_name,
        "rmse", as.character(k), bench_rmse(y_true, mean_k))
      rows[[length(rows) + 1L]] <- bench_record(
        scenario_meta$scenario_id, p, "gdpar", family_name,
        "mae", as.character(k), bench_mae(y_true, mean_k))
      rows[[length(rows) + 1L]] <- bench_record(
        scenario_meta$scenario_id, p, "gdpar", family_name,
        "crps", as.character(k),
        bench_crps_gaussian(y_true, mean_k, sd_k))
    }
  }
  lo <- bench_gdpar_loo(fit_obj)
  if (!is.null(lo)) {
    pk <- bench_pareto_k_summary(lo)
    rows[[length(rows) + 1L]] <- bench_record(
      scenario_meta$scenario_id, p, "gdpar", family_name,
      "elpd_loo", "agg",
      lo$estimates["elpd_loo", "Estimate"],
      se_value = lo$estimates["elpd_loo", "SE"],
      convergence_flag = bench_convergence_flag_from_loo(lo),
      max_pareto_k = pk$max_pareto_k,
      n_pareto_k_above_07 = pk$n_above_07)
  } else {
    rows[[length(rows) + 1L]] <- bench_record(
      scenario_meta$scenario_id, p, "gdpar", family_name,
      "elpd_loo", "agg", NA_real_,
      convergence_flag = "loo_failed")
  }
  rows[[length(rows) + 1L]] <- bench_record(
    scenario_meta$scenario_id, p, "gdpar", family_name,
    "wall", "fit", fit_obj$wall_seconds,
    wall_seconds = fit_obj$wall_seconds)
  do.call(rbind, rows)
}

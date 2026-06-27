## rstanarm wrapper for the Block 6 multi adversarial bench.
##
## Implementation policy: per-coordinate fits. Routing per formula:
##   - if the RHS contains a (...|...) random-effects term (lme4
##     syntax), switch to rstanarm::stan_glmer for that coord
##   - otherwise fall back to rstanarm::stan_glm
## stan_mvmer flagged non-deterministic across versions for
## heterogeneous families; per-coord routing is robust and consistent
## with the coord-wise factorization. PSIS-LOO via rstanarm::loo()
## per coord summed over k for subject aggregation.

.bench_rstanarm_has_random_effect <- function(formula_k) {
  rhs_txt <- paste(deparse(formula_k), collapse = " ")
  grepl("\\([^)]*\\|", rhs_txt)
}

bench_rstanarm_fit <- function(data, family_spec, scenario_meta, ...) {
  if (!requireNamespace("rstanarm", quietly = TRUE)) {
    stop("Package 'rstanarm' is required.")
  }
  iter <- bench_default_iter()
  p <- as.integer(scenario_meta$p)
  fams <- family_spec$rstanarm_family
  formulas <- scenario_meta$rstanarm_formula_per_coord
  fits <- vector("list", p)
  t0 <- Sys.time()
  for (k in seq_len(p)) {
    fit_fn <- if (.bench_rstanarm_has_random_effect(formulas[[k]])) {
      rstanarm::stan_glmer
    } else {
      rstanarm::stan_glm
    }
    fits[[k]] <- suppressMessages(suppressWarnings(
      fit_fn(
        formula = formulas[[k]], data = data, family = fams[[k]],
        chains = iter$chains,
        iter = iter$iter_warmup + iter$iter_sampling,
        warmup = iter$iter_warmup, seed = iter$seed + k,
        refresh = iter$refresh, cores = iter$chains
      )
    ))
  }
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  list(fits = fits, wall_seconds = wall, p = p,
       family_spec = family_spec, formulas = formulas)
}

bench_rstanarm_predict <- function(fit_obj, newdata) {
  p <- fit_obj$p
  fam_name <- fit_obj$family_spec$brms_family
  mean_mat <- matrix(0, nrow = nrow(newdata), ncol = p)
  sd_mat   <- matrix(0, nrow = nrow(newdata), ncol = p)
  prob_mat <- if (fam_name == "bernoulli")
    matrix(0, nrow = nrow(newdata), ncol = p) else NULL
  draws_list <- vector("list", p)
  for (k in seq_len(p)) {
    fit_k <- fit_obj$fits[[k]]
    draws_pred <- rstanarm::posterior_predict(fit_k, newdata = newdata)
    draws_mean <- rstanarm::posterior_epred(fit_k, newdata = newdata)
    mean_mat[, k] <- colMeans(draws_mean)
    sd_mat[, k]   <- apply(draws_pred, 2L, stats::sd)
    if (fam_name == "bernoulli") prob_mat[, k] <- colMeans(draws_mean)
    draws_list[[k]] <- draws_pred
  }
  list(mean = if (p == 1L) mean_mat[, 1L] else mean_mat,
       sd   = if (p == 1L) sd_mat[, 1L]   else sd_mat,
       prob = if (is.null(prob_mat)) NULL
              else if (p == 1L) prob_mat[, 1L] else prob_mat,
       draws = draws_list)
}

bench_rstanarm_loo <- function(fit_obj) {
  if (!requireNamespace("loo", quietly = TRUE)) return(NULL)
  p <- fit_obj$p
  ll_list <- lapply(fit_obj$fits, function(fit_k) {
    rstanarm::log_lik(fit_k)
  })
  S <- nrow(ll_list[[1L]])
  n <- ncol(ll_list[[1L]])
  log_lik_subject <- matrix(0, nrow = S, ncol = n)
  for (k in seq_len(p)) log_lik_subject <- log_lik_subject + ll_list[[k]]
  iter <- bench_default_iter()
  chain_id <- rep(seq_len(iter$chains), each = S / iter$chains)
  r_eff <- loo::relative_eff(exp(log_lik_subject), chain_id = chain_id)
  tryCatch(loo::loo(log_lik_subject, r_eff = r_eff),
            error = function(e) {
              warning("rstanarm loo failed: ", conditionMessage(e))
              NULL
            })
}

bench_rstanarm_metrics <- function(fit_obj, test, scenario_meta) {
  preds <- bench_rstanarm_predict(fit_obj, test)
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
    draws_k <- preds$draws[[k]]
    if (classification) {
      y_int <- as.integer(y_true)
      rows[[length(rows) + 1L]] <- bench_record(
        scenario_meta$scenario_id, p, "rstanarm", fam_name,
        "auc", as.character(k), bench_auc(y_int, prob_k))
      rows[[length(rows) + 1L]] <- bench_record(
        scenario_meta$scenario_id, p, "rstanarm", fam_name,
        "log_loss", as.character(k), bench_logloss(y_int, prob_k))
      rows[[length(rows) + 1L]] <- bench_record(
        scenario_meta$scenario_id, p, "rstanarm", fam_name,
        "brier", as.character(k), bench_brier(y_int, prob_k))
    } else {
      rows[[length(rows) + 1L]] <- bench_record(
        scenario_meta$scenario_id, p, "rstanarm", fam_name,
        "rmse", as.character(k), bench_rmse(y_true, mean_k))
      rows[[length(rows) + 1L]] <- bench_record(
        scenario_meta$scenario_id, p, "rstanarm", fam_name,
        "mae", as.character(k), bench_mae(y_true, mean_k))
      rows[[length(rows) + 1L]] <- bench_record(
        scenario_meta$scenario_id, p, "rstanarm", fam_name,
        "crps", as.character(k), bench_crps_sample(y_true, draws_k))
    }
  }
  lo <- bench_rstanarm_loo(fit_obj)
  if (!is.null(lo)) {
    pk <- bench_pareto_k_summary(lo)
    rows[[length(rows) + 1L]] <- bench_record(
      scenario_meta$scenario_id, p, "rstanarm", fam_name,
      "elpd_loo", "agg",
      lo$estimates["elpd_loo", "Estimate"],
      se_value = lo$estimates["elpd_loo", "SE"],
      convergence_flag = bench_convergence_flag_from_loo(lo),
      max_pareto_k = pk$max_pareto_k,
      n_pareto_k_above_07 = pk$n_above_07)
  } else {
    rows[[length(rows) + 1L]] <- bench_record(
      scenario_meta$scenario_id, p, "rstanarm", fam_name,
      "elpd_loo", "agg", NA_real_,
      convergence_flag = "loo_failed")
  }
  rows[[length(rows) + 1L]] <- bench_record(
    scenario_meta$scenario_id, p, "rstanarm", fam_name,
    "wall", "fit", fit_obj$wall_seconds,
    wall_seconds = fit_obj$wall_seconds)
  do.call(rbind, rows)
}

## Shared helpers for the Block 6 multi-method bench harness.
## Sourced from each bench_competitor_<method>.R and from the runner.
##
## All metric functions accept y_true (numeric vector or 0/1 integer
## vector for classification) and either:
##   - point predictions:   y_pred (numeric vector of same length)
##   - probabilistic:       mean_pred + sd_pred (gaussian) OR
##                           draws_mat (n_draws x n) sample-based, OR
##                           prob_pred (vector in [0, 1]).
## The runner decides which to call based on family_spec.

bench_record <- function(scenario_id, p, method, family_name,
                          metric, coord, value,
                          se_value = NA_real_,
                          wall_seconds = NA_real_,
                          convergence_flag = "ok",
                          max_pareto_k = NA_real_,
                          n_pareto_k_above_07 = NA_integer_,
                          recovered_rescor = NA_real_,
                          extra = NA_character_) {
  data.frame(
    scenario_id        = scenario_id,
    p                  = as.integer(p),
    method             = method,
    family_name        = family_name,
    metric             = metric,
    coord              = coord,
    value              = as.numeric(value),
    se_value           = as.numeric(se_value),
    wall_seconds       = as.numeric(wall_seconds),
    convergence_flag   = as.character(convergence_flag),
    max_pareto_k       = as.numeric(max_pareto_k),
    n_pareto_k_above_07 = as.integer(n_pareto_k_above_07),
    recovered_rescor   = as.numeric(recovered_rescor),
    extra              = as.character(extra),
    stringsAsFactors   = FALSE
  )
}

bench_rmse <- function(y_true, y_pred) {
  sqrt(mean((as.numeric(y_true) - as.numeric(y_pred))^2))
}

bench_mae <- function(y_true, y_pred) {
  mean(abs(as.numeric(y_true) - as.numeric(y_pred)))
}

bench_brier <- function(y_true, prob_pred) {
  mean((as.numeric(y_true) - as.numeric(prob_pred))^2)
}

bench_logloss <- function(y_true, prob_pred) {
  eps <- 1e-15
  p <- pmin(pmax(as.numeric(prob_pred), eps), 1 - eps)
  y <- as.numeric(y_true)
  -mean(y * log(p) + (1 - y) * log(1 - p))
}

bench_auc <- function(y_true, prob_pred) {
  y <- as.integer(y_true)
  s <- as.numeric(prob_pred)
  n_pos <- sum(y == 1L); n_neg <- sum(y == 0L)
  if (n_pos == 0L || n_neg == 0L) return(NA_real_)
  ranks <- rank(s)
  (sum(ranks[y == 1L]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

bench_crps_gaussian <- function(y_true, mean_pred, sd_pred) {
  if (!requireNamespace("scoringRules", quietly = TRUE)) {
    return(NA_real_)
  }
  mean(scoringRules::crps_norm(y = as.numeric(y_true),
                                mean = as.numeric(mean_pred),
                                sd = as.numeric(sd_pred)))
}

bench_crps_sample <- function(y_true, draws_mat) {
  if (!requireNamespace("scoringRules", quietly = TRUE)) {
    return(NA_real_)
  }
  if (is.null(draws_mat) || !is.matrix(draws_mat)) return(NA_real_)
  n <- length(y_true)
  if (nrow(draws_mat) == n && ncol(draws_mat) != n) {
    dat <- draws_mat
  } else if (ncol(draws_mat) == n && nrow(draws_mat) != n) {
    dat <- t(draws_mat)
  } else if (nrow(draws_mat) == n && ncol(draws_mat) == n) {
    dat <- draws_mat
  } else {
    return(NA_real_)
  }
  mean(scoringRules::crps_sample(y = as.numeric(y_true), dat = dat))
}

bench_pareto_k_summary <- function(loo_obj) {
  if (is.null(loo_obj)) return(list(max_pareto_k = NA_real_,
                                       n_above_07 = NA_integer_))
  pk <- loo_obj$diagnostics$pareto_k
  if (is.null(pk)) return(list(max_pareto_k = NA_real_,
                                  n_above_07 = NA_integer_))
  list(max_pareto_k = as.numeric(max(pk, na.rm = TRUE)),
       n_above_07   = as.integer(sum(pk > 0.7, na.rm = TRUE)))
}

bench_time_fit <- function(expr) {
  t0 <- Sys.time()
  out <- force(expr)
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  list(value = out, wall_seconds = wall)
}

bench_convergence_flag_from_loo <- function(loo_obj) {
  if (is.null(loo_obj)) return("loo_failed")
  pk_n <- bench_pareto_k_summary(loo_obj)$n_above_07
  if (is.na(pk_n)) return("loo_failed")
  if (pk_n > 0L)   return("pareto_k_high")
  "ok"
}

bench_default_iter <- function() {
  list(chains = 2L, iter_warmup = 500L, iter_sampling = 1000L,
       seed = 42L, refresh = 0L)
}

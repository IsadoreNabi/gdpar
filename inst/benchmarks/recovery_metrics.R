## Block 9.2.S: recovery metric module (DESIGN_9_2_EXTERNAL section 2.4).
##
## Operates on a fixed evaluation grid with KNOWN ground truth. Each
## adapter (recovery_adapters.R) returns a standard prediction object:
##   $status      "ok" | "unsupported" | "pending_impl" | "fit_failed" |
##                 "metrics_failed" | "canonical_abort"
##   $eyx_draws   S x n_grid matrix of posterior draws of E[y|x] on the
##                grid (for the mean function + its 90% interval).
##   $sd_hat      length-n_grid posterior-mean conditional SD on grid.
##   $tail_q_hat  length-n_grid plug-in 0.95 predictive quantile on grid.
##   $loo         psis_loo object on test (elpd cross-check) or NULL.
##   $diag        list(rhat_max, ess_bulk_min, n_divergences).
##   $wall_seconds, $error_class, $gdpar_family, $K.
##
## The metrics are computed against truth_grid = list(mu, sd, tail_q095)
## from recovery_realize_cell(). Strata interpretation (criterion 2.S,
## DESIGN 2.5): rmse_sd grades heteroscedasticity; tail_metric grades
## heavy tails; rmse_mean/coverage90 are global; on autocorrelation the
## only requirement is finiteness (robustness, not leadership).

## --- core scalar metrics -------------------------------------------

rec_rmse <- function(est, truth) {
  sqrt(mean((as.numeric(est) - as.numeric(truth))^2))
}
rec_bias <- function(est, truth) {
  mean(as.numeric(est) - as.numeric(truth))
}

## Coverage of the true mean function by the empirical 90% credible
## interval built from the draws of E[y|x].
rec_coverage90 <- function(eyx_draws, mu_truth) {
  if (is.null(eyx_draws) || !is.matrix(eyx_draws)) return(NA_real_)
  qs <- apply(eyx_draws, 2L, stats::quantile, probs = c(0.05, 0.95),
              names = FALSE, na.rm = TRUE)
  lo <- qs[1L, ]; hi <- qs[2L, ]
  mean(mu_truth >= lo & mu_truth <= hi, na.rm = TRUE)
}

## Build the long-tidy metric rows for one (scenario, method, rep) cell.
## Returns a one-row data.frame of the manifest metric columns. A failed
## or unsupported prediction yields finite-NA metrics but a recorded
## status (the coverage gate counts status, not the metric values).
recovery_compute_metrics <- function(pred, truth_grid) {
  status <- pred$status %||% "fit_failed"
  out <- list(
    status        = status,
    wall_seconds  = pred$wall_seconds %||% NA_real_,
    rmse_mean     = NA_real_, bias_mean = NA_real_,
    rmse_sd       = NA_real_, tail_metric = NA_real_,
    coverage90    = NA_real_, elpd_loo = NA_real_,
    pareto_k_max  = NA_real_, rhat_max = NA_real_,
    ess_bulk_min  = NA_real_, n_divergences = NA_real_,
    error_class   = pred$error_class %||% NA_character_
  )
  if (!identical(status, "ok")) {
    return(as.data.frame(out, stringsAsFactors = FALSE))
  }
  mu_hat <- if (!is.null(pred$eyx_draws)) colMeans(pred$eyx_draws)
            else pred$mu_hat
  out$rmse_mean   <- rec_rmse(mu_hat, truth_grid$mu)
  out$bias_mean   <- rec_bias(mu_hat, truth_grid$mu)
  out$coverage90  <- rec_coverage90(pred$eyx_draws, truth_grid$mu)
  if (!is.null(pred$sd_hat)) {
    out$rmse_sd <- rec_rmse(pred$sd_hat, truth_grid$sd)
  }
  if (!is.null(pred$tail_q_hat)) {
    out$tail_metric <- rec_rmse(pred$tail_q_hat, truth_grid$tail_q095)
  }
  if (!is.null(pred$loo)) {
    est <- tryCatch(pred$loo$estimates["elpd_loo", "Estimate"],
                    error = function(e) NA_real_)
    out$elpd_loo <- as.numeric(est)
    pk <- pred$loo$diagnostics$pareto_k
    out$pareto_k_max <- if (!is.null(pk)) max(pk, na.rm = TRUE) else NA_real_
  }
  if (!is.null(pred$diag)) {
    out$rhat_max      <- pred$diag$rhat_max %||% NA_real_
    out$ess_bulk_min  <- pred$diag$ess_bulk_min %||% NA_real_
    out$n_divergences <- pred$diag$n_divergences %||% NA_real_
  }
  as.data.frame(out, stringsAsFactors = FALSE)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

## --- shared family-moment helpers ----------------------------------
## Map fitted slot parameters (on the grid) to (eyx_draws, sd_hat,
## tail_q_hat). Adapters call these so the family algebra lives in ONE
## place. loc/scale/etc. are matrices [S, n_grid] of posterior draws
## (response scale); scalars are recycled.

## Continuous Gaussian / Student-t. loc = mean draws; sd_mat = scale
## (SD) draws; nu = degrees of freedom (scalar or draws; Inf -> normal).
rec_moments_continuous <- function(loc_mat, sd_mat, nu = Inf) {
  eyx <- loc_mat                                  ## E[y|x] = location
  if (all(is.infinite(nu))) {
    sd_eff <- sd_mat
    tail_mat <- loc_mat + sd_mat * stats::qnorm(0.95)
  } else {
    nu_v <- if (length(nu) == 1L) nu else nu
    nu_safe <- pmax(nu_v, 2.01)
    sd_eff <- sd_mat * sqrt(nu_safe / (nu_safe - 2))
    tail_mat <- loc_mat + sd_mat * stats::qt(0.95, df = nu_safe)
  }
  list(eyx_draws = eyx,
       sd_hat = colMeans(sd_eff),
       tail_q_hat = colMeans(tail_mat))
}

## Count Poisson / Neg-Binomial with optional zero-inflation. lambda_mat
## = rate draws [S, n_grid]; size = nb size (Inf -> Poisson); pi_mat =
## structural-zero prob draws (0 -> no inflation). tail_q computed at the
## posterior-mean parameters (plug-in) to avoid per-draw CDF inversion.
rec_moments_count <- function(lambda_mat, size = Inf, pi_mat = 0) {
  if (length(pi_mat) == 1L) {
    pi_mat <- matrix(pi_mat, nrow = nrow(lambda_mat),
                     ncol = ncol(lambda_mat))
  }
  eyx <- (1 - pi_mat) * lambda_mat
  base_var <- if (is.finite(size)) lambda_mat + lambda_mat^2 / size
              else lambda_mat
  ey2 <- (1 - pi_mat) * (base_var + lambda_mat^2)
  sd_mat <- sqrt(pmax(ey2 - eyx^2, 0))
  lam_bar <- colMeans(lambda_mat)
  pi_bar  <- colMeans(pi_mat)
  tail_q <- rec_count_tail_plugin(lam_bar, pi_bar, size, prob = 0.95)
  list(eyx_draws = eyx,
       sd_hat = colMeans(sd_mat),
       tail_q_hat = tail_q)
}

## Plug-in 0.95 quantile of the (zero-inflated) count predictive.
rec_count_tail_plugin <- function(lambda, pi, size, prob = 0.95) {
  n <- length(lambda)
  out <- numeric(n)
  for (i in seq_len(n)) {
    lam_i <- max(lambda[i], 1e-8); pi_i <- pi[i]
    k <- 0L
    repeat {
      f0 <- if (is.finite(size)) stats::pnbinom(k, size = size, mu = lam_i)
            else stats::ppois(k, lambda = lam_i)
      f <- pi_i + (1 - pi_i) * f0
      if (f >= prob || k > 1e6L) break
      k <- k + 1L
    }
    out[i] <- k
  }
  out
}

## Build a standard prediction object skeleton for non-ok outcomes.
rec_pred_status <- function(status, error_class = NA_character_,
                            wall_seconds = NA_real_,
                            gdpar_family = NA_character_, K = NA_integer_) {
  list(status = status, error_class = error_class,
       wall_seconds = wall_seconds, gdpar_family = gdpar_family, K = K,
       eyx_draws = NULL, sd_hat = NULL, tail_q_hat = NULL,
       loo = NULL, diag = NULL)
}

## Block RG, sub-phase RG.7 -- the CERTIFIED-LIMIT endpoint for the four
## near-deterministic Tweedie cells of benchmark 9.2.O
## (N__N1__gdpar__{nef,mac,app,grl}).
##
## VERDICT (B9.33). The eBird count Tweedie posterior (condition number of the
## observed information ~ 3-5e8 at the mode) is NOT sampleable at feasible cost.
## Default NUTS, constant-mass Laplace HMC, AND whitened NUTS (dense metric =
## the Laplace precision, init at the mode, 1000 warmup + 4 x 1000) all FAIL the
## multi-chain gate: whitened NUTS returns R-hat = 1.117, ESS-bulk = 22.5, and
## 193 / 4000 = 5% divergences (~/rg33_whiten_full.out). Divergences in the
## ALREADY-WHITENED coordinates prove genuinely non-Gaussian geometry -- the
## canyon curves; the 90-degree Fisher / Hessian misalignment (B9.32) foretold
## it -- so a frozen Gaussian metric cannot track it and a per-step Riemannian
## metric is infeasible (minutes per exact Hessian). Both lineages (GLM-5.2 and
## MiMo-2.5-Pro) converged on, and the user pre-approved, the honest endpoint:
## CERTIFY the limit and report the LAPLACE / plug-in ELPD. This is exactly
## competitor parity -- mgcv/REML and INLA report a Laplace/plug-in predictive,
## O(d/n)-exact -- accompanied by a sensitivity check. No further samplers (the
## stopping rule of both lineages). The "RESOLVED" of B9.32 was a false positive
## of a single-chain acceptance gate; the multi-chain R-hat/ESS gate (a GLM/MiMo
## contribution, now a permanent diagnostic) caught it.
##
## What this module does, per Tweedie cell, by the CLEAN route:
##   (1) find the posterior MODE on the unconstrained scale (L-BFGS warm start
##       from the best on-ridge draw + a modified-Newton polish, rg7_climb...);
##   (2) the EXACT Hessian at the mode (cmdstan autodiff) -> precision
##       M = -Hessian = observed information (sign!), eigen-floored to PD;
##   (3) draw theta_s ~ N(mode, M^{-1}) = mode + M^{-1/2} z_s (the Laplace
##       posterior approximation);
##   (4) CONSTRAIN each theta_s through the compiled model
##       (constrain_variables with the transformed parameters) and INJECT the
##       resulting constrained draws as the posterior of a gdpar_fit, so the
##       PUBLIC predict() / .ebird_gdpar_param_draws() produce the predictive
##       (mu, phi, p) draws on the test set and the SHARED Tweedie density +
##       Jacobian machinery (ebird_metrics_9_2_O.R) scores the ELPD identically
##       to every competitor;
##   (5) a SENSITIVITY check: the Laplace-averaged ELPD vs the plug-in (mode)
##       ELPD, and the per-point spread of the log predictive density across the
##       Laplace draws. A small spread / delta means the posterior-uncertainty
##       contribution (the only place the Gaussian approximation acts) is small,
##       so the Laplace ELPD is close to the (intractable) exact-posterior ELPD.
##
## Nothing here touches the package R/ sources (the automated in-library Laplace
## fallback is a separate, later step) -> the R CMD check baseline is unchanged.
## All curvature helpers (rg7_best_reference, rg7_climb_to_mode,
## rg7_observed_information) are reused from rg7_geometry_helpers.R.

`%||%` <- function(a, b) if (is.null(a)) b else a

## ------------------------------------------------------------------
## (A) Mode + Laplace precision at the mode.
##     `cached` may carry a $mode (an unconstrained vector from an earlier fit of
##     the SAME cell): the unconstrained coordinate system is a deterministic
##     function of the Stan model + data, independent of the sampler budget/seed,
##     so a cached mode is a valid WARM START here -- but we always re-polish
##     (cheap) so the reported mode is reproducible from this run alone.
## ------------------------------------------------------------------
rg7_laplace_mode_precision <- function(bridge, fit, cached = NULL,
                                       eigen_floor_rel = 1e-10,
                                       climb_steps = 300L, verbose = TRUE) {
  gt <- bridge$geom_target
  d <- as.integer(gt$dim)
  flog <- function(...) if (verbose) cat(format(Sys.time(), "%H:%M:%S"),
                                         sprintf(...), "\n")

  start <- NULL; mode_src <- "fresh"
  if (!is.null(cached) && !is.null(cached$mode)) {
    cm <- as.numeric(cached$mode)
    if (length(cm) == d) {
      gn <- tryCatch(sqrt(sum(gt$grad_log_prob(cm)^2)),
                     error = function(e) Inf)
      if (is.finite(gn)) { start <- cm; mode_src <- sprintf("cached|grad=%.3g", gn) }
    }
  }
  if (is.null(start)) {
    ref0 <- rg7_best_reference(bridge, fit)
    neg_lp <- function(th) {
      v <- tryCatch(-gt$log_prob(th), error = function(e) Inf)
      if (!is.finite(v)) 1e18 else v
    }
    neg_grd <- function(th) {
      g <- tryCatch(-gt$grad_log_prob(th), error = function(e) rep(0, length(th)))
      g[!is.finite(g)] <- 0; g
    }
    opt <- tryCatch(stats::optim(ref0, neg_lp, neg_grd, method = "L-BFGS-B",
                                 control = list(maxit = 5000L, factr = 1e1)),
                    error = function(e) NULL)
    start <- if (!is.null(opt)) opt$par else ref0
  }
  cl <- rg7_climb_to_mode(gt, start, n_steps = climb_steps, tol = 1e-3)
  mode <- as.numeric(cl$mode)
  gnorm <- sqrt(sum(gt$grad_log_prob(mode)^2))
  logp <- gt$log_prob(mode)
  flog("mode: src=%s |grad|=%.4g logp=%.2f dim=%d conv=%s",
       mode_src, gnorm, logp, d, cl$converged)

  ## precision = observed information = -Hessian(log p); exact cmdstan Hessian
  ## preferred (rg7_observed_information records the method).
  M <- rg7_observed_information(gt, mode)
  M <- 0.5 * (M + t(M))
  eM <- eigen(M, symmetric = TRUE)
  all_pos <- all(eM$values > 0)
  ## A genuine maximum has -Hessian positive-definite. A negative raw eigenvalue
  ## means the climb stalled at a saddle, not the mode; the Laplace Gaussian
  ## would then be ill-defined (GLM/MiMo review, MINOR-5). The eigen-floor below
  ## keeps the draw machinery numerically alive, but we WARN loudly and record
  ## the flag rather than silently certify a saddle.
  if (!all_pos)
    warning(sprintf(paste("rg7 Laplace: -Hessian at the mode has %d non-positive",
                          "eigenvalue(s) (min=%.3g) -- possible saddle; Laplace",
                          "approximation flagged."),
                    sum(eM$values <= 0), min(eM$values)))
  floor_eig <- max(eM$values) * eigen_floor_rel
  lam <- pmax(eM$values, floor_eig)
  logdetM <- sum(log(lam))
  Hinv  <- eM$vectors %*% (t(eM$vectors) / lam)         # M^{-1} (Laplace covariance)
  Lhalf <- eM$vectors %*% (t(eM$vectors) / sqrt(lam))   # M^{-1/2} (symmetric)
  ## Mode-accuracy diagnostic (GLM/MiMo review, MAJOR-1): the residual gradient
  ## norm is not directly interpretable under cond ~1e8, but the Newton
  ## displacement M^{-1} g, measured in posterior SDs, is. |M^{-1} g| <=
  ## |g| / lambda_min bounds the worst-direction offset of the reported mode
  ## from the true mode in SD units; report it honestly.
  g_mode <- gt$grad_log_prob(mode)
  newton_step <- as.numeric(Hinv %*% g_mode)            # M^{-1} g
  mode_offset_sd <- sqrt(max(sum(g_mode * newton_step), 0))  # Newton decrement = sqrt(g^T M^{-1} g) = ||offset||_M in SD
  flog("precision(-H): cond=%.3e min_eig=%.4g all_pos=%s method=%s mode_offset=%.3g SD",
       max(lam) / min(lam), min(eM$values), all_pos,
       attr(M, "method") %||% "unknown", mode_offset_sd)

  list(mode = mode, M = M, Hinv = Hinv, Lhalf = Lhalf, logdetM = logdetM,
       grad_norm = gnorm, logp = logp, eig = eM$values, lam = lam,
       eig_min = min(eM$values), eig_max = max(eM$values),
       cond = max(lam) / min(lam), all_pos = all_pos,
       mode_offset_sd = mode_offset_sd,
       method = attr(M, "method") %||% "unknown", mode_source = mode_src,
       converged = isTRUE(cl$converged), dim = d)
}

## ------------------------------------------------------------------
## (B) Laplace draws on the UNCONSTRAINED scale: theta_s = mode + M^{-1/2} z_s.
##     iid (no autocorrelation), so S of them carry more information than S
##     correlated NUTS draws.
## ------------------------------------------------------------------
rg7_laplace_unconstrained_draws <- function(mp, S, seed) {
  ## local RNG: save + restore the global seed so the extractor leaves the
  ## caller's RNG stream untouched (GLM/MiMo review, MINOR-7).
  old_seed <- if (exists(".Random.seed", envir = globalenv())) {
    get(".Random.seed", envir = globalenv())
  } else NULL
  on.exit(if (!is.null(old_seed))
    assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
  set.seed(seed)
  d <- mp$dim
  Z <- matrix(stats::rnorm(S * d), nrow = S, ncol = d)
  U <- Z %*% mp$Lhalf                 # row s = (Lhalf z_s)^T (Lhalf symmetric)
  sweep(U, 2L, mp$mode, "+")          # + mode  => cov(row) = Lhalf Lhalf = M^{-1}
}

## ------------------------------------------------------------------
## (B2) FIT-QUALITY diagnostic of the Laplace Gaussian against the TRUE
##      unconstrained posterior (GLM/MiMo review, MAJOR-4). The per-point sd of
##      the log predictive density (the sensitivity check below) measures
##      WITHIN-Laplace spread, NOT the gap between the Laplace ELPD and the
##      intractable exact-posterior ELPD: a skewed posterior (the canyon curves)
##      can leave that gap large even when the within-Laplace spread is moderate.
##      The right diagnostic is how well the Laplace Gaussian q approximates the
##      true posterior p over the SAME draws:
##        - the self-normalised importance-sampling ESS with log-weights
##          w_s = log p(theta_s|y) - log q(theta_s); ESS << S flags a poor
##          Gaussian match (and, via PSIS, a Pareto-k tail index);
##        - the mean / max log-density drop  log p(mode|y) - log p(theta_s|y),
##          which is ~ chi^2_d / 2 (mean d/2) under a perfect Gaussian; a much
##          larger mean flags non-Gaussian (lighter/heavier) tails.
##      Cheap: one gt$log_prob() per draw (the quadratic log q is vectorised).
## ------------------------------------------------------------------
rg7_laplace_fit_quality <- function(gt, mp, U) {
  S <- nrow(U); d <- mp$dim
  cen <- sweep(U, 2L, mp$mode, "-")
  quad <- rowSums((cen %*% mp$M) * cen)                  # (theta-mode)^T M (theta-mode)
  log_q <- -0.5 * d * log(2 * pi) + 0.5 * mp$logdetM - 0.5 * quad
  log_p <- vapply(seq_len(S), function(s)
    tryCatch(gt$log_prob(U[s, ]), error = function(e) NA_real_), numeric(1L))
  fin <- is.finite(log_p) & is.finite(log_q)
  lw <- log_p[fin] - log_q[fin]                          # log importance weights
  lw <- lw - max(lw)
  w <- exp(lw)
  ess <- (sum(w)^2) / sum(w^2)
  ## PSIS Pareto-k of the importance ratios = the canonical reliability index of
  ## the Laplace proposal vs the true posterior (k < 0.7 reliable).
  pk <- tryCatch(
    if (requireNamespace("loo", quietly = TRUE) && sum(fin) > 10L)
      suppressWarnings(loo::psis(matrix(log_p[fin] - log_q[fin], ncol = 1L),
                                 r_eff = NA)$diagnostics$pareto_k)
    else NA_real_,
    error = function(e) NA_real_)
  drop <- mp$logp - log_p[fin]                            # log p(mode) - log p(theta_s)
  list(ess_is = ess, ess_is_frac = ess / sum(fin), pareto_k = pk,
       logdrop_mean = mean(drop), logdrop_max = max(drop),
       logdrop_expected = d / 2, n_finite = sum(fin), n_total = S)
}

## ------------------------------------------------------------------
## (C) Flatten a constrain_variables() named list into a single named numeric
##     vector with posterior bracket names in COLUMN-MAJOR order (the convention
##     cmdstanr's $draws() and posterior use, hence what predict's regexes
##     parse: theta_ref_k[1,k], a_coef_k[k,j], ...). Verified against the model's
##     own draws by rg7_laplace_verify_flatten() before any cell is trusted.
## ------------------------------------------------------------------
.rg7_flatten_constrained <- function(cv) {
  out <- vector("list", length(cv))
  nms <- names(cv)
  for (i in seq_along(cv)) {
    v <- cv[[i]]; nm <- nms[i]; dd <- dim(v)
    if (length(v) == 0L) next        # drop empty params (e.g. mu_theta_ref_k
                                     # declared (0,3) under Option A/D96); they
                                     # carry no values and predict never reads them
    if (is.null(dd)) {
      out[[i]] <- if (length(v) == 1L) stats::setNames(as.numeric(v), nm)
                  else stats::setNames(as.numeric(v),
                                       paste0(nm, "[", seq_along(v), "]"))
    } else {
      ai <- arrayInd(seq_len(prod(dd)), dd)               # column-major indices
      lbl <- apply(ai, 1L, paste, collapse = ",")
      out[[i]] <- stats::setNames(as.numeric(v), paste0(nm, "[", lbl, "]"))
    }
  }
  unlist(out)
}

## The exhaustive set of draw variables predict_from_newdata_K() reads (R/
## methods.R): the per-slot reference theta_ref_k, the additive a_coef_k, the
## multiplicative c_b_k, and the globally-shared modulating W_raw / sigma_W.
## Everything else returned by constrain_variables (raw parameters, the
## generated quantity theta_i_k over the training rows) is irrelevant to
## newdata prediction and is dropped to keep the injected draws small.
.RG7_PREDICT_VARS <- c("theta_ref_k", "a_coef_k", "c_b_k", "W_raw", "sigma_W")

## ------------------------------------------------------------------
## (D) GROUND-TRUTH self-check of the flatten + naming. The model's own
##     posterior draws ARE the constrained images of its unconstrained draws, so
##     constrain_variables() on the unconstrained draws must reproduce the draws
##     exactly. This de-risks the only bespoke logic (the column-major bracket
##     naming of array[K] vector[J] containers): it certifies, against the
##     model's ground truth, that the injected draws carry the right values
##     under the names predict() parses. Returns the max abs discrepancy over
##     the predict-relevant columns; the driver asserts it is ~0.
## ------------------------------------------------------------------
rg7_laplace_verify_flatten <- function(fobj, fit, n_check = 4L) {
  um <- as.matrix(posterior::as_draws_matrix(fit$fit$unconstrain_draws()))
  dm <- as.matrix(posterior::as_draws_matrix(fit$fit$draws()))
  S0 <- nrow(um); n_check <- min(n_check, S0)
  rows <- unique(round(seq(1, S0, length.out = n_check)))
  base_of <- sub("\\[.*$", "", colnames(dm))
  cols <- colnames(dm)[base_of %in% .RG7_PREDICT_VARS]
  if (!length(cols)) {
    return(list(ok = FALSE, max_abs_diff = NA_real_,
                reason = "no predict-relevant variables in real draws"))
  }
  worst <- 0; missing <- character(0)
  for (s in rows) {
    cv <- fobj$constrain_variables(unconstrained_variables = um[s, ],
                                   transformed_parameters = TRUE,
                                   generated_quantities = FALSE)
    fl <- .rg7_flatten_constrained(cv)
    miss <- setdiff(cols, names(fl)); missing <- union(missing, miss)
    common <- intersect(cols, names(fl))
    worst <- max(worst, max(abs(fl[common] - dm[s, common])))
  }
  list(ok = (length(missing) == 0L) && is.finite(worst) && worst < 1e-6,
       max_abs_diff = worst, missing = missing, n_vars = length(cols),
       rows_checked = length(rows))
}

## ------------------------------------------------------------------
## (E) Build the injected synthetic-posterior draws from the Laplace
##     unconstrained draws U (S x d): constrain each, flatten, keep the
##     predict-relevant columns, stack into a posterior draws_matrix.
## ------------------------------------------------------------------
rg7_laplace_synth_draws <- function(fobj, U, verbose = TRUE) {
  S <- nrow(U)
  rows <- vector("list", S); n_fail <- 0L
  for (s in seq_len(S)) {
    cv <- tryCatch(
      fobj$constrain_variables(unconstrained_variables = U[s, ],
                               transformed_parameters = TRUE,
                               generated_quantities = FALSE),
      error = function(e) NULL)
    if (is.null(cv)) { n_fail <- n_fail + 1L; next }
    rows[[s]] <- .rg7_flatten_constrained(cv)
  }
  ok <- !vapply(rows, is.null, TRUE)
  if (!any(ok)) stop("rg7 Laplace: every constrain_variables() call failed.")
  ref_names <- names(rows[[which(ok)[1L]]])
  base_of <- sub("\\[.*$", "", ref_names)
  keep <- ref_names[base_of %in% .RG7_PREDICT_VARS]
  mat <- matrix(NA_real_, nrow = sum(ok), ncol = length(keep),
                dimnames = list(NULL, keep))
  r <- 0L
  for (s in which(ok)) { r <- r + 1L; mat[r, ] <- rows[[s]][keep] }
  ## Non-random draw dropping bias (GLM/MiMo review, MAJOR-2): extreme draws in
  ## the stiffest eigen-directions can overflow bounded transforms and fail to
  ## constrain; because failures correlate with the Gaussian tails, the survivors
  ## are a TRUNCATED Laplace posterior, which can bias the ELPD. The eigen-floor
  ## (well below the smallest real eigenvalue at cond ~1e8 -> the floor never
  ## binds) keeps the tails finite, so the rate is expected ~0; we GATE on it and
  ## warn loudly if non-trivial so the bias is never silent.
  fail_rate <- n_fail / S
  if (fail_rate > 0.005)
    warning(sprintf(paste("rg7 Laplace: %d/%d (%.2f%%) draws failed to",
                          "constrain -- the surviving Laplace draws are",
                          "truncated; the ELPD may be biased. Inspect."),
                    n_fail, S, 100 * fail_rate))
  if (verbose) cat(sprintf("  [synth] draws ok=%d fail=%d (%.3f%%) vars=%d\n",
                           sum(ok), n_fail, 100 * fail_rate, length(keep)))
  list(draws = posterior::as_draws_matrix(mat), n_ok = sum(ok),
       n_fail = n_fail, fail_rate = fail_rate, kept_vars = keep)
}

## ------------------------------------------------------------------
## (F) Predictive (mu, phi, p) draws on the test set from the injected Laplace
##     posterior. predict.gdpar_fit() touches object$fit ONLY through
##     object$fit$draws() (no args; R/methods.R line 154); a minimal stand-in
##     whose $draws() returns the synthetic Laplace draws makes the PUBLIC
##     predict path consume them with ZERO changes to gdpar. All other fit
##     metadata (K, slot_names, family$param_specs, amm_list_canonical,
##     design_K, anchor, parametrization) is the live refit's, so the linear
##     predictor reconstruction + per-slot inverse links are gdpar's own.
## ------------------------------------------------------------------
rg7_laplace_param_draws <- function(fit, fobj, U, test, level = "Nt",
                                    verbose = TRUE) {
  syn <- rg7_laplace_synth_draws(fobj, U, verbose = verbose)
  fit_laplace <- fit
  fit_laplace$fit <- list(draws = function(...) syn$draws)
  pd <- .ebird_gdpar_param_draws(fit_laplace, test, level)
  ## predict()'s S x n x K array drops to a length-n vector per slot when S = 1
  ## (the plug-in / mode draw); the shared metric functions expect S x n
  ## matrices, so re-attach the single-draw row dimension.
  pd <- lapply(pd, function(x) if (is.matrix(x)) x else matrix(x, nrow = 1L))
  list(param_draws = pd, n_ok = syn$n_ok, n_fail = syn$n_fail,
       fail_rate = syn$fail_rate, kept_vars = syn$kept_vars)
}

## ------------------------------------------------------------------
## (G) Sensitivity of the Laplace ELPD to posterior uncertainty: the
##     Laplace-averaged ELPD vs the plug-in (mode) ELPD on the original measure,
##     and the per-test-point spread of the log predictive density across the
##     Laplace draws. coord/track default to the Tweedie count coordinate
##     (coord 2, track "N": Jacobian 0, density already on the original scale).
## ------------------------------------------------------------------
rg7_laplace_sensitivity <- function(param_draws, mode_param_draws, y_orig,
                                    coord = 2L, track = "N",
                                    family_tag = "tweedie") {
  ll <- ebird_loglik_orig_from_params(family_tag, param_draws, y_orig,
                                      coord, track)          # S x n
  pw <- apply(ll, 2L, .ebird_logmeanexp_col)
  fin <- is.finite(pw)
  elpd_laplace <- sum(pw[fin])
  ## Monte Carlo SE of the ELPD across test points (GLM/MiMo review, MINOR-6):
  ## the pointwise contributions are independent across test rows, so the total
  ## SE is sqrt(n) * sd(pw). With S=2000 iid draws this is ~0.5 ELPD units.
  se_elpd <- if (sum(fin) > 1L) sqrt(sum(fin)) * stats::sd(pw[fin]) else NA_real_
  ll0 <- ebird_loglik_orig_from_params(family_tag, mode_param_draws, y_orig,
                                       coord, track)          # 1 x n
  elpd_plugin <- sum(ll0[is.finite(ll0)])
  sd_pt <- apply(ll, 2L, function(c) {
    g <- c[is.finite(c)]; if (length(g) > 1L) stats::sd(g) else NA_real_
  })
  list(elpd_laplace = elpd_laplace, elpd_plugin = elpd_plugin,
       delta_elpd = elpd_laplace - elpd_plugin, se_elpd = se_elpd,
       sd_logdens_median = stats::median(sd_pt, na.rm = TRUE),
       sd_logdens_mean = mean(sd_pt, na.rm = TRUE),
       sd_logdens_max = max(sd_pt, na.rm = TRUE),
       n_test = length(y_orig))
}

## ------------------------------------------------------------------
## (H) End-to-end Laplace ELPD for ONE Tweedie cell. Refits the Tweedie scalar
##     (Nt) to obtain a LIVE fit + bridge (the cell stores only param draws),
##     finds the mode + precision, draws S Laplace draws, injects them, predicts
##     (mu, phi, p) on the test set, and returns the param draws + the mode
##     (plug-in) param draws + the full Laplace certification metadata. The
##     caller turns param_draws into the coord-2 result via ebird_finalize_coord
##     and pairs it with the coord-1 beta fit.
## ------------------------------------------------------------------
rg7_laplace_tweedie_coord <- function(train, test, fit_seed, cached = NULL,
                                      S = 2000L, refit_iters = NULL,
                                      refit_treedepth = 8L, verbose = TRUE) {
  flog <- function(...) if (verbose) cat(format(Sys.time(), "%H:%M:%S"),
                                         sprintf(...), "\n")
  t0 <- Sys.time()
  ri <- refit_iters %||% list(chains = 2L, iter_warmup = 300L,
                              iter_sampling = 300L)
  flog("refit tweedie (live fit for the bridge + predict metadata) ...")
  rf <- rg7_refit_tweedie(train, iters = ri, fit_seed = fit_seed,
                          max_treedepth = refit_treedepth)
  fit <- rf$fit; tr <- rf$train
  fobj <- fit$fit
  flog("bridge (hessian=TRUE) ...")
  bridge <- gdpar_geom_bridge(fit, hessian = TRUE)

  flog("verify flatten/naming against the model's own draws ...")
  vf <- rg7_laplace_verify_flatten(fobj, fit)
  if (!isTRUE(vf$ok)) {
    stop(sprintf(paste("rg7 Laplace flatten self-check FAILED",
                       "(max_abs_diff=%.3g, missing={%s}). The injected draws",
                       "would not reproduce the model -- aborting this cell."),
                 vf$max_abs_diff %||% NA, paste(vf$missing, collapse = ",")))
  }
  flog("  flatten OK: max_abs_diff=%.3g over %d vars, %d rows",
       vf$max_abs_diff, vf$n_vars, vf$rows_checked)

  mp <- rg7_laplace_mode_precision(bridge, fit, cached = cached,
                                   verbose = verbose)

  U <- rg7_laplace_unconstrained_draws(mp, S = S, seed = fit_seed + 700L)
  ## Laplace-vs-true-posterior fit quality on the SAME draws (the honest bound on
  ## the Laplace-to-exact ELPD gap; the per-point sd alone does not bound it).
  flog("Laplace fit-quality (IS-ESS + log-density drop) ...")
  fq <- rg7_laplace_fit_quality(bridge$geom_target, mp, U)
  flog("  fit-quality: ESS_is=%.1f (%.1f%% of %d) pareto_k=%.3g logdrop mean=%.2f (exp d/2=%.1f) max=%.2f",
       fq$ess_is, 100 * fq$ess_is_frac, fq$n_finite, fq$pareto_k,
       fq$logdrop_mean, fq$logdrop_expected, fq$logdrop_max)

  flog("Laplace draws S=%d -> constrain + inject + predict ...", S)
  pdr <- rg7_laplace_param_draws(fit, fobj, U, test, level = "Nt",
                                 verbose = verbose)
  ## plug-in (mode) predictive: a single "draw" at the mode.
  mpd <- rg7_laplace_param_draws(fit, fobj, matrix(mp$mode, nrow = 1L),
                                 test, level = "Nt", verbose = FALSE)

  sens <- rg7_laplace_sensitivity(pdr$param_draws, mpd$param_draws, test$cnt,
                                  coord = 2L, track = "N",
                                  family_tag = "tweedie")
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  flog("Laplace ELPD=%.2f +/- %.2f (plug-in=%.2f, delta=%.2f) sd_logdens med=%.3g max=%.3g | cond=%.2e |grad|=%.3g offset=%.3gSD fail=%.2f%% wall=%.0fs",
       sens$elpd_laplace, sens$se_elpd %||% NA, sens$elpd_plugin,
       sens$delta_elpd, sens$sd_logdens_median, sens$sd_logdens_max,
       mp$cond, mp$grad_norm, mp$mode_offset_sd,
       100 * (pdr$fail_rate %||% NA), wall)

  list(param_draws = pdr$param_draws, mode_param_draws = mpd$param_draws,
       sensitivity = sens, fit_quality = fq, mode_precision = mp,
       verify_flatten = vf, n_draw_ok = pdr$n_ok, n_draw_fail = pdr$n_fail,
       fail_rate = pdr$fail_rate, kept_vars = pdr$kept_vars, wall_seconds = wall)
}

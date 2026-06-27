## Block RG, sub-phase RG.7 -- applying the geometry-adaptive capability to
## the REAL Tweedie count of benchmark 9.2.O (the four N__N1__gdpar cells).
##
## These helpers re-fit a Tweedie count cell, bridge it to the geometry engine,
## furnish the local curvature at the reference position by the user's chosen
## CROSS-CHECK route (observed information by finite differences of the gradient
## + expected Fisher by simulation through the real compiled model), confirm the
## two agree on the wall/floor split, and run the geometry-adaptive orchestrator
## (and a forced sub-Riemannian pass) on the genuine posterior.
##
## Nothing here touches the package R/ sources: the whole RG.7 application lives
## under inst/benchmarks/, so the R CMD check baseline is unchanged. The engine,
## bridge, orchestrator, sub-Riemannian metric and Fisher simulator it calls are
## all already in the package (RG.1-RG.6).
##
## Honesty (charter section 5; the inherited caveat): the eBird count is
## near-deterministic (condition number ~ 1e6; RG.4 validated cond ~ 90). Both a
## RESOLVED outcome (an honest ELPD) and a CERTIFIED reproducible LIMIT are
## first-class results; neither is advertised from partial coverage.

## ------------------------------------------------------------------
## (1) Re-fit the Tweedie scalar (count) of one sub-region -> live fit.
##     Mirrors .ebird_fit_gdpar_scalar()'s level "Nt" config EXACTLY so the
##     posterior is the same object the 9.2.O cell would have produced
##     (the cell stores only param draws, not the fit; we need a live fit for
##     the bridge). `iters` overrides the sampler budget for smoke vs full.
## ------------------------------------------------------------------
## `max_treedepth`: the full 9.2.O cell uses 10; RG.7 only needs the fit to
## supply the bridge ingredients (compiled model + a representative reference
## on the posterior ridge + Stan data), since the GEOMETRY sampler -- not NUTS
## -- produces the final draws, so a shallower, faster fit on the same posterior
## suffices and is honest (the reference is a point on the same canyon ridge).
rg7_refit_tweedie <- function(train, iters = NULL, fit_seed = 91079L,
                              parallel_chains = NULL, max_treedepth = 9L) {
  rhs <- paste(.ebird_9_2_O_covariates, collapse = " + ")
  has_g <- nlevels(droplevels(factor(train$group))) >= 2L
  loc_rhs <- if (has_g) paste(rhs, "+ group") else rhs
  loc <- stats::as.formula(sprintf("y ~ a(%s)", loc_rhs))
  fs <- gdpar_bf(loc, phi ~ 1, p ~ 1)               # level "Nt" (tweedie K=3)
  tr <- train; tr$y <- tr$cnt
  it <- iters %||% list(chains = 2L, iter_warmup = 500L, iter_sampling = 1000L)
  pc <- parallel_chains %||% it$chains
  fit <- gdpar(fs, family = gdpar_family("tweedie"), data = tr,
               chains = it$chains, parallel_chains = pc,
               iter_warmup = it$iter_warmup, iter_sampling = it$iter_sampling,
               adapt_delta = 0.95, max_treedepth = max_treedepth, refresh = 0L,
               seed = fit_seed, skip_id_check = TRUE, verbose = FALSE,
               parametrization = "auto")
  list(fit = fit, train = tr)
}

## ------------------------------------------------------------------
## (1b) A robust reference position: the highest-log-density posterior draw
##      (the most "on-ridge" point), not the posterior mean of a possibly
##      unconverged multi-chain fit (whose mean can sit OFF the canyon ridge,
##      between chains, where the local curvature read is noise). The wall/floor
##      split is a global property of the ridge, so any on-ridge point reads it
##      consistently; the mode-like best-density draw is the natural choice.
## ------------------------------------------------------------------
rg7_best_reference <- function(bridge, fit) {
  um <- tryCatch(
    as.matrix(posterior::as_draws_matrix(fit$fit$unconstrain_draws())),
    error = function(e) NULL)
  if (is.null(um) || nrow(um) == 0L) return(bridge$reference)
  lp <- bridge$geom_target$log_prob
  lps <- apply(um, 1L, function(th) tryCatch(lp(th), error = function(e) -Inf))
  if (all(!is.finite(lps))) return(bridge$reference)
  unname(um[which.max(lps), ])
}

## ------------------------------------------------------------------
## (1c) Climb from the best-density draw to the posterior MODE by damped Newton
##      steps using the engine gradient + exact Hessian. RG.7 found that the
##      sub-Riemannian remedy must read the curvature at the genuine mode: a
##      reference left off-mode (the mean / a draw of a stuck NUTS pilot) makes
##      the expected Fisher and the observed information point in different
##      directions, and a metric built from the wrong directions preconditions
##      phantom walls (acceptance collapses to zero). The optimiser traverses
##      the canyon ridge the way the penalised-likelihood competitors do, so a
##      few Newton steps reach the mode where the exact Hessian is the true,
##      positive-definite local curvature. Damped + backtracking for safety.
## ------------------------------------------------------------------
rg7_climb_to_mode <- function(geom_target, ref, n_steps = 150L, ridge = 1e-6,
                              tol = 1e-4, verbose = FALSE, max_step_norm = 5) {
  lp <- geom_target$log_prob
  gl <- geom_target$grad_log_prob
  he <- geom_target$hessian
  if (!is.function(he)) return(list(mode = ref, converged = FALSE,
                                    reason = "no exact Hessian"))
  d <- length(ref)
  th <- ref; f0 <- lp(th)
  for (it in seq_len(n_steps)) {
    g <- gl(th); H <- he(th)
    if (any(!is.finite(g))) break       # gradient unusable -> stop at best so far.
    # Modified-Newton ascent of log p when the exact Hessian is FINITE: form
    # A = -H (pos-def at a maximum) and FLOOR its eigenvalues to a positive
    # minimum, so an indefinite Hessian (saddle / not-yet-mode) still yields an
    # ascent direction A^{-1} g. More robust than a fixed ridge when the
    # curvature spans ~10^6. A GRADIENT-ONLY damped step when the Hessian is
    # NON-FINITE: the stiff off-mode region can return a non-finite Hessian,
    # which must NOT abort the climb (eigen() would throw) -- the climb has to
    # TRAVERSE that region to reach the benign mode where the Hessian is finite
    # and positive-definite. So eigen() is never fed a non-finite matrix.
    if (all(is.finite(H))) {
      H <- 0.5 * (H + t(H))
      e <- eigen(-H, symmetric = TRUE)
      lam <- pmax(e$values, ridge * max(abs(e$values), 1))
      step <- as.numeric(e$vectors %*% ((crossprod(e$vectors, g)) / lam))
    } else {
      step <- g                         # steepest ascent fallback.
    }
    if (any(!is.finite(step))) break
    # trust-region cap on the step length.
    sn <- sqrt(sum(step^2))
    if (sn > max_step_norm) step <- step * (max_step_norm / sn)
    # backtracking line search on log p.
    t <- 1; ok <- FALSE
    for (bt in seq_len(20L)) {
      cand <- th + t * step
      fc <- tryCatch(lp(cand), error = function(e) -Inf)
      if (is.finite(fc) && fc >= f0 + 1e-4 * t * sum(g * step)) { ok <- TRUE; break }
      t <- t / 2
    }
    if (!ok) break
    gnorm <- sqrt(sum(g^2))
    if (verbose) cat(sprintf("  [mode] it=%d logp=%.4f |grad|=%.3g step=%.3g\n",
                             it, fc, gnorm, t))
    improve <- fc - f0
    th <- cand; f0 <- fc
    if (gnorm < tol || improve < tol * 1e-2) {
      return(list(mode = th, converged = TRUE, iters = it, grad_norm = gnorm,
                  logp = f0))
    }
  }
  list(mode = th, converged = FALSE, logp = f0,
       grad_norm = sqrt(sum(gl(th)^2)))
}

## ------------------------------------------------------------------
## (1d) Multi-start mode climb: climb from the top-K highest-density posterior
##      draws and keep the cleanest mode (smallest gradient norm). Records
##      whether the starts AGREE on one basin (max pairwise distance between the
##      climbed modes, relative to the posterior scale) -- a frozen Hessian read
##      at a non-unique mode would precondition the wrong basin (GLM/MiMo review,
##      rank 3-5). A few extra Newton climbs are cheap relative to the HMC run.
## ------------------------------------------------------------------
rg7_climb_multistart <- function(bridge, fit, geom_target, K = 3L, ...) {
  um <- tryCatch(as.matrix(posterior::as_draws_matrix(fit$fit$unconstrain_draws())),
                 error = function(e) NULL)
  lp <- geom_target$log_prob
  refs <- list()
  if (!is.null(um) && nrow(um) > 0L) {
    lps <- apply(um, 1L, function(th) tryCatch(lp(th), error = function(e) -Inf))
    ord <- order(lps, decreasing = TRUE)
    ord <- ord[is.finite(lps[ord])]
    for (i in head(ord, K)) refs[[length(refs) + 1L]] <- unname(um[i, ])
  }
  if (!length(refs)) refs <- list(bridge$reference)
  climbs <- lapply(refs, function(r) rg7_climb_to_mode(geom_target, r, ...))
  gnorms <- vapply(climbs, function(c) c$grad_norm %||% Inf, 0)
  best <- climbs[[which.min(gnorms)]]
  ## basin agreement: max pairwise distance between climbed modes.
  basin_spread <- NA_real_
  if (length(climbs) >= 2L) {
    M <- do.call(rbind, lapply(climbs, `[[`, "mode"))
    dm <- as.matrix(stats::dist(M))
    basin_spread <- max(dm)
  }
  list(mode = best$mode, converged = best$converged,
       grad_norm = best$grad_norm %||% NA_real_, logp = best$logp %||% NA_real_,
       n_starts = length(climbs), grad_norms = gnorms,
       basin_spread = basin_spread)
}

## ------------------------------------------------------------------
## (1e) Multi-chain convergence diagnostics for the geometric draws. A frozen
##      mass at a non-converged mode can accept ~1 while only exploring a slice
##      of the canyon: acceptance + divergences alone DO NOT prove resolution
##      (GLM/MiMo review, rank 1). R-hat / ESS across independent chains is the
##      real proof. `draws_list` = list of (iter x d) matrices, one per chain.
## ------------------------------------------------------------------
rg7_multichain_diag <- function(draws_list) {
  draws_list <- Filter(function(x) !is.null(x) && is.matrix(x) && nrow(x) > 1L,
                       draws_list)
  if (length(draws_list) < 2L || !requireNamespace("posterior", quietly = TRUE))
    return(list(max_rhat = NA_real_, min_ess_bulk = NA_real_,
                n_chains = length(draws_list)))
  d <- ncol(draws_list[[1L]])
  niter <- min(vapply(draws_list, nrow, 0L))
  m <- length(draws_list)
  arr <- array(NA_real_, dim = c(niter, m, d))
  for (ci in seq_len(m)) arr[, ci, ] <- draws_list[[ci]][seq_len(niter), , drop = FALSE]
  rh <- apply(arr, 3L, function(x) tryCatch(posterior::rhat(x),
                                            error = function(e) NA_real_))
  eb <- apply(arr, 3L, function(x) tryCatch(posterior::ess_bulk(x),
                                            error = function(e) NA_real_))
  list(max_rhat = suppressWarnings(max(rh, na.rm = TRUE)),
       min_ess_bulk = suppressWarnings(min(eb, na.rm = TRUE)),
       n_chains = m, n_iter = niter, rhat = rh, ess_bulk = eb)
}

## ------------------------------------------------------------------
## (2) Observed information at the reference: -Hessian of the log-posterior,
##     symmetrised. Prefers the EXACT cmdstan Hessian (higher-order autodiff,
##     which -- found in RG.7 -- DOES compile for this Tweedie model) over
##     finite differences: on a near-deterministic posterior the curvature spans
##     ~10 orders of magnitude, where a fixed finite-difference step cannot
##     resolve the stiff wall directions (its condition number blows up to ~1e17
##     = numerical noise, not the geometry). The exact Hessian avoids that.
##     Finite differences remain the fallback when no Hessian is available, with
##     a Richardson-extrapolated step. Records which route was used.
## ------------------------------------------------------------------
rg7_observed_information <- function(geom_target, reference, h = 1e-4) {
  d <- geom_target$dim
  if (is.function(geom_target$hessian)) {
    H <- tryCatch(geom_target$hessian(reference), error = function(e) NULL)
    if (!is.null(H) && all(is.finite(H))) {
      obs <- -0.5 * (H + t(H))
      attr(obs, "method") <- "exact_hessian"
      return(obs)
    }
  }
  gl <- geom_target$grad_log_prob
  J <- matrix(0, d, d)
  for (j in seq_len(d)) {
    ep <- reference; ep[j] <- ep[j] + h
    em <- reference; em[j] <- em[j] - h
    J[, j] <- (gl(ep) - gl(em)) / (2 * h)
  }
  obs <- -0.5 * (J + t(J))           # -Hessian(log post), symmetrised.
  attr(obs, "method") <- if (all(is.finite(obs))) "finite_difference" else
                         "finite_difference_nonfinite"
  obs
}

## ------------------------------------------------------------------
## (3) Expected Fisher by simulation THROUGH THE REAL MODEL.
##     simulate(theta_ref): draw a synthetic count from the model's Tweedie
##       likelihood at the posterior-mean predictive (mu_i, phi, p).
##     score(theta_ref, y_sim): re-instantiate the SAME compiled model with the
##       synthetic outcome (only y_real swapped; the prior anchor stays fixed)
##       and read grad_log_prob at theta_ref -> the per-data-set posterior score
##       s_post = s_lik + g_prior.
##     The expected-likelihood Fisher is the CENTRED second moment of the scores
##     (E[s_lik]=0 under the model, so the mean score estimates the fixed prior
##     gradient g_prior, and the covariance removes the rank-one prior term):
##       I_sim = mean(s s^T) - sbar sbar^T.
##     This is the gdpar_geom_fisher_simulator estimator (PSD by construction,
##     unbiased) made robust to the heavy real-model score (non-finite scores
##     skipped). Evaluated once, at the reference.
## ------------------------------------------------------------------
rg7_simulated_fisher <- function(bridge, fit, train, n_sim = 64L, seed = 9L,
                                 verbose = TRUE, reference = NULL) {
  ref <- reference %||% bridge$reference
  d <- bridge$dim
  model <- bridge$model
  base_data <- bridge$stan_data
  stopifnot(!is.null(base_data$y_real))            # Tweedie is continuous.
  n <- length(base_data$y_real)

  ## predictive (mu_i, phi, p) at the posterior mean, from the fitted model.
  arr <- stats::predict(fit, newdata = train, type = "response",
                        summary = "draws")          # S x n x 3
  mu_i <- colMeans(arr[, , 1L])
  phi  <- mean(arr[, , 2L])
  p    <- mean(arr[, , 3L])
  ## Root-cause instrumentation (GLM/MiMo review, rank 2): a degenerate
  ## predictive (mu/phi/p non-finite or pathological) is a model-side signal,
  ## NOT merely a fragile witness -- it would also threaten the metric path.
  pred_finite <- all(is.finite(mu_i)) && is.finite(phi) && is.finite(p) &&
                 phi > 0 && p > 1 && p < 2
  if (!requireNamespace("tweedie", quietly = TRUE)) {
    stop("package 'tweedie' is required to simulate the Tweedie count.")
  }
  simulate_count <- function() {
    y <- tweedie::rtweedie(n, mu = mu_i, phi = phi, power = p)
    y[!is.finite(y)] <- 0
    pmax(y, 0)
  }
  ## re-instantiation score: 1-iteration sample with the swapped outcome
  ## exposes grad_log_prob at theta_ref for the synthetic data set. Returns the
  ## score (ok) or a status string ("sample_fail"/"grad_fail") for attribution.
  score_at_ref <- function(y_sim) {
    sd2 <- base_data; sd2$y_real <- as.double(y_sim)
    f <- tryCatch(model$sample(data = sd2, chains = 1L, iter_warmup = 1L,
                               iter_sampling = 1L, refresh = 0L,
                               show_messages = FALSE, show_exceptions = FALSE),
                  error = function(e) NULL)
    if (is.null(f)) return("sample_fail")
    g <- tryCatch({
      f$init_model_methods(seed = 1L, verbose = FALSE, hessian = FALSE)
      as.numeric(f$grad_log_prob(unconstrained_variables = ref))
    }, error = function(e) NULL)
    if (is.null(g) || length(g) != d || any(!is.finite(g))) return("grad_fail")
    g
  }

  set.seed(seed)
  S <- matrix(NA_real_, n_sim, d)
  n_ok <- 0L; n_fail <- 0L; n_sample_fail <- 0L; n_grad_fail <- 0L
  for (i in seq_len(n_sim)) {
    s <- score_at_ref(simulate_count())
    if (is.character(s)) {
      n_fail <- n_fail + 1L
      if (identical(s, "sample_fail")) n_sample_fail <- n_sample_fail + 1L
      else n_grad_fail <- n_grad_fail + 1L
      next
    }
    n_ok <- n_ok + 1L; S[n_ok, ] <- s
    if (verbose && (i %% 8L == 0L))
      cat(sprintf("  [fisher-sim] %d/%d (ok=%d fail=%d)\n", i, n_sim, n_ok, n_fail))
  }
  ## Degenerate witness: too few valid scores (the heavy real-model scores can
  ## come back non-finite for the stiffest sub-regions) -> the centred second
  ## moment would be a NaN matrix. Return I_centred = NULL (a clean degenerate
  ## marker), NOT a NaN matrix: the simulated Fisher is a cross-check WITNESS,
  ## never the metric, so a degenerate witness must NEVER abort a resolvable
  ## cell. The caller proceeds on the (finite, all-positive) observed Hessian.
  if (n_ok < d + 1L) {
    warning(sprintf("simulated Fisher degenerate: %d valid scores < dim %d + 1; witness skipped.",
                    n_ok, d))
    return(list(I_centred = NULL, I_raw = NULL, g_prior = NULL,
                n_ok = n_ok, n_fail = n_fail, degenerate = TRUE,
                n_sample_fail = n_sample_fail, n_grad_fail = n_grad_fail,
                pred_finite = pred_finite,
                predictive = list(mu_range = range(mu_i), phi = phi, p = p)))
  }
  Sok <- S[seq_len(n_ok), , drop = FALSE]
  sbar <- colMeans(Sok)
  raw <- crossprod(Sok) / n_ok                       # mean(s s^T)
  centred <- raw - tcrossprod(sbar)                  # expected-likelihood Fisher
  centred <- 0.5 * (centred + t(centred))
  if (!all(is.finite(centred))) {                    # belt-and-suspenders finitude guard.
    warning("simulated Fisher centred matrix non-finite; witness skipped.")
    return(list(I_centred = NULL, I_raw = NULL, g_prior = sbar,
                n_ok = n_ok, n_fail = n_fail, degenerate = TRUE,
                n_sample_fail = n_sample_fail, n_grad_fail = n_grad_fail,
                pred_finite = pred_finite,
                predictive = list(mu_range = range(mu_i), phi = phi, p = p)))
  }
  list(I_centred = centred, I_raw = 0.5 * (raw + t(raw)),
       g_prior = sbar, n_ok = n_ok, n_fail = n_fail, degenerate = FALSE,
       n_sample_fail = n_sample_fail, n_grad_fail = n_grad_fail,
       pred_finite = pred_finite,
       predictive = list(mu_range = range(mu_i), phi = phi, p = p))
}

## ------------------------------------------------------------------
## (4) Cross-check the two curvature estimates: do they agree on the wall/floor
##     split? Reports each spectrum, the condition numbers, and the principal
##     angle between their dominant (wall) eigenspaces -- in the near-
##     deterministic regime observed information and expected Fisher coincide
##     (Bernstein-von Mises rigidity), so a small angle confirms the metric is
##     reading a real, stable wall/floor structure.
## ------------------------------------------------------------------
rg7_crosscheck <- function(J_obs, I_sim, floor = 1e-8, n_wall = NULL) {
  ## Defense in depth: the witness may be degenerate (I_sim = NULL or non-finite
  ## from too few valid simulated scores). A non-informative cross-check is a
  ## valid outcome; an error is NOT -- the cross-check never feeds the metric.
  if (is.null(J_obs) || is.null(I_sim) ||
      !all(is.finite(J_obs)) || !all(is.finite(I_sim))) {
    return(list(obs_spectrum = NA_real_, sim_spectrum = NA_real_,
                cond_obs = NA_real_, cond_sim = NA_real_, k_wall = NA_integer_,
                max_principal_angle_deg = NA_real_, top_eig_log10_rmsd = NA_real_,
                degenerate = TRUE))
  }
  spec <- function(M) {
    e <- eigen(0.5 * (M + t(M)), symmetric = TRUE)
    list(values = e$values, vectors = e$vectors,
         pos = e$values[e$values > floor],
         cond = {
           p <- e$values[e$values > floor]
           if (length(p)) max(p) / min(p) else NA_real_
         })
  }
  so <- spec(J_obs); si <- spec(I_sim)
  d <- nrow(J_obs)
  ## dominant-subspace agreement: principal angles between the top-k eigvecs.
  k <- n_wall %||% max(1L, sum(so$values > stats::median(so$values[so$values > floor]) ))
  k <- min(k, d - 1L, length(si$values))
  Uo <- so$vectors[, seq_len(k), drop = FALSE]
  Ui <- si$vectors[, seq_len(k), drop = FALSE]
  sv <- svd(crossprod(Uo, Ui))$d
  sv <- pmin(pmax(sv, -1), 1)
  max_angle_deg <- acos(min(sv)) * 180 / pi
  ## eigenvalue agreement on the log scale (top-k).
  lo <- log10(pmax(sort(so$values, decreasing = TRUE)[seq_len(k)], floor))
  li <- log10(pmax(sort(si$values, decreasing = TRUE)[seq_len(k)], floor))
  list(obs_spectrum = sort(so$values), sim_spectrum = sort(si$values),
       cond_obs = so$cond, cond_sim = si$cond,
       k_wall = k, max_principal_angle_deg = max_angle_deg,
       top_eig_log10_rmsd = sqrt(mean((lo - li)^2)), degenerate = FALSE)
}

## ------------------------------------------------------------------
## (5) Apply the geometry-adaptive capability to the genuine posterior, using
##     the WINNING recipe discovered in RG.7 (B9.32 step-size sweep): the real
##     eBird tweedie IS sampleable. At the Newton-climbed mode the exact Hessian
##     is all-positive (Bernstein-von Mises holds: the posterior is locally
##     near-Gaussian), so a CONSTANT dense Euclidean mass equal to that Hessian
##     (Laplace preconditioning = the Riemannian level with a frozen metric)
##     resolves it cleanly at a small step (eps ~ 0.01-0.03: accept ~1, 0
##     divergences, E-BFMI ~1.3). The sub-Riemannian level works too but needs a
##     smaller step and mixes worse (a constant metric suffices here -- no need
##     for a position-dependent one); the identity-mass Euclidean never works
##     (condition ~5e8 unpreconditioned). This function therefore:
##       (5a) PRIMARY: dense Euclidean (mass = exact Hessian at the mode), with a
##            short epsilon mini-tune over the resolved band -> the resolving
##            draws.
##       (5b) RECORD: the sub-Riemannian pass (for the comparison / the ladder).
##       (5c) RECORD: the identity-Euclidean baseline (the canyon-wrecks-it
##            contrast).
##     `curvature` is the PSD curvature feeding the metrics (the exact Hessian at
##     the mode = the observed information; the simulated Fisher is the
##     cross-check witness, not the metric source -- RG.7 found it points the
##     wrong way off-mode). `reference` must be the climbed mode.
## ------------------------------------------------------------------
rg7_apply_geometry <- function(bridge, curvature, reference, geom = NULL,
                               seed = 20260604L, verbose = TRUE,
                               eps_grid = c(0.03, 0.02, 0.01),
                               run_subriemannian = TRUE,
                               n_chains = 4L, init_jitter = 2) {
  ## init_jitter = OVERDISPERSION factor applied to the Laplace-covariance
  ## chain inits (metric-aware), not a flat absolute jitter (see Minv_half).
  gt <- bridge$geom_target
  g <- geom %||% list(n_warmup = 500L, n_iter = 1000L, L = 24L,
                      tune_iter = 80L)
  ## The curvature (exact Hessian at the mode) feeds the Laplace metric; if it
  ## is unusable, return a structured non-resolved result instead of letting the
  ## eigen() below throw -- a certified limit is a first-class outcome, a crash
  ## is not.
  if (is.null(curvature) || !all(is.finite(curvature)))
    return(list(curvature_cond = NA_real_, tune = list(), best_eps = NA_real_,
                dense_euclidean = NULL, sub_riemannian = NULL, sub_meta = NULL,
                identity_euclidean = NULL, resolved = FALSE,
                error = "non-finite curvature: cannot build the Laplace metric"))
  e_cur <- eigen(0.5 * (curvature + t(curvature)), symmetric = TRUE)
  lam_cur <- pmax(e_cur$values, 1e-8)
  Jpsd <- e_cur$vectors %*% (lam_cur * t(e_cur$vectors))
  ## M^{-1/2} for OVERDISPERSED, metric-aware chain inits: a flat jitter is
  ## catastrophic here (cond ~1e8 -> a flat 0.01 is ~hundreds of posterior SDs
  ## in the stiff directions, so jittered chains start deep in the wall and
  ## diverge, inflating R-hat artificially). Disperse along the Laplace
  ## covariance instead: big in soft directions, tiny in stiff ones.
  Minv_half <- e_cur$vectors %*% (t(e_cur$vectors) / sqrt(lam_cur))

  ## (5a) PRIMARY: dense Euclidean (Laplace). Mini-tune eps for a healthy band.
  m_dense <- tryCatch(gdpar_geom_metric_euclidean(M = Jpsd),
                      error = function(e) NULL)
  tune <- list()
  best <- list(eps = NA_real_, accept = -1, fit = NULL)
  if (!is.null(m_dense)) {
    for (eps in eps_grid) {
      ft <- tryCatch(gdpar_geom_hmc(gt, metric = m_dense, epsilon = eps,
                                    L = g$L, n_iter = g$tune_iter,
                                    n_warmup = g$tune_iter, init = reference,
                                    seed = seed),
                     error = function(e) NULL)
      acc <- if (is.null(ft)) NA_real_ else ft$accept_rate
      dv  <- if (is.null(ft)) NA_real_ else ft$n_divergent
      tune[[length(tune) + 1L]] <- list(eps = eps, accept = acc, div = dv)
      if (verbose) cat(sprintf("  [tune] eps=%.4g accept=%.3f div=%s\n",
                               eps, acc %||% NA, dv %||% NA))
      ## a healthy resolved fit: high acceptance, no divergences.
      score <- if (is.null(ft) || !is.finite(acc)) -Inf else
        acc - 5 * (dv / (2 * g$tune_iter))
      if (score > best$accept) best <- list(eps = eps, accept = score, fit = ft)
    }
  }
  ## Final resolution as n_chains INDEPENDENT chains: chain 1 from the climbed
  ## mode, the rest from jittered inits, so R-hat / ESS across chains can prove
  ## the sampler explores the whole canyon (acceptance + divergences alone do
  ## NOT -- a too-conservative frozen mass accepts ~1 while crawling a slice;
  ## GLM/MiMo review, rank 1). dense_full = the primary (mode-started) chain.
  dense_full <- NULL; dense_chains <- NULL; dense_diag <- NULL
  if (!is.na(best$eps)) {
    chains <- vector("list", n_chains)
    for (ci in seq_len(n_chains)) {
      init_ci <- if (ci == 1L) reference else
        reference + init_jitter *
          as.numeric(Minv_half %*% stats::rnorm(length(reference)))
      chains[[ci]] <- tryCatch(
        gdpar_geom_hmc(gt, metric = m_dense, epsilon = best$eps, L = g$L,
                       n_iter = g$n_iter, n_warmup = g$n_warmup, init = init_ci,
                       seed = seed + 5L + ci),
        error = function(e) structure(list(error = conditionMessage(e)),
                                      class = "rg7_hmc_error"))
      if (verbose && !is.null(chains[[ci]]$accept_rate))
        cat(sprintf("  [chain %d/%d] accept=%.3f div=%s ebfmi=%.3f\n", ci, n_chains,
                    chains[[ci]]$accept_rate, chains[[ci]]$n_divergent,
                    chains[[ci]]$ebfmi %||% NA))
    }
    dense_chains <- chains
    dense_full <- chains[[1L]]
    ok_draws <- Filter(function(c) is.null(c$error) && !is.null(c$draws),
                       chains)
    if (length(ok_draws) >= 2L)
      dense_diag <- rg7_multichain_diag(lapply(ok_draws, `[[`, "draws"))
    if (verbose && !is.null(dense_diag))
      cat(sprintf("  [multichain] chains=%d max_rhat=%.4g min_ess=%.4g\n",
                  dense_diag$n_chains, dense_diag$max_rhat, dense_diag$min_ess_bulk))
  }

  ## (5b) RECORD the sub-Riemannian pass (smaller step, the ladder's level 5).
  sub <- NULL; sub_meta <- NULL
  if (isTRUE(run_subriemannian)) {
    m_sub <- tryCatch(gdpar_geom_metric_subriemannian(
      gt, fisher = function(th) Jpsd, reference = reference),
      error = function(e) NULL)
    if (!is.null(m_sub)) {
      sub_meta <- list(n_walls = m_sub$n_walls, tau = m_sub$tau,
                       suggested_epsilon = m_sub$suggested_epsilon)
      sub <- tryCatch(
        gdpar_geom_hmc(gt, metric = m_sub, epsilon = 5e-4, L = g$L,
                       n_iter = min(g$n_iter, 400L),
                       n_warmup = min(g$n_warmup, 200L), init = reference,
                       seed = seed + 9L),
        error = function(e) structure(list(error = conditionMessage(e)),
                                      class = "rg7_hmc_error"))
    }
  }

  ## (5c) RECORD the identity-Euclidean baseline (the canyon contrast).
  euc_id <- tryCatch(
    gdpar_geom_hmc(gt, epsilon = best$eps %||% 0.01, L = g$L, n_iter = 100L,
                   n_warmup = 60L, init = reference, seed = seed + 13L),
    error = function(e) structure(list(error = conditionMessage(e)),
                                  class = "rg7_hmc_error"))

  ## RESOLVED requires a HEALTHY sampler on ALL diagnostics, not just acceptance
  ## + divergences: E-BFMI >= 0.3 (mass not catastrophically scaled) AND
  ## multi-chain R-hat < 1.05 + ESS-bulk >= 100 (chains agree on the whole
  ## posterior, ruling out a frozen-mass random walk on a canyon slice).
  acc_ok   <- !is.null(dense_full) && is.null(dense_full$error) &&
              is.finite(dense_full$accept_rate) && dense_full$accept_rate > 0.6
  div_ok   <- !is.null(dense_full) && is.finite(dense_full$n_divergent) &&
              dense_full$n_divergent / (g$n_warmup + g$n_iter) < 0.02
  ebfmi_ok <- !is.null(dense_full) && is.finite(dense_full$ebfmi %||% NA_real_) &&
              (dense_full$ebfmi %||% 0) >= 0.3
  conv_ok  <- !is.null(dense_diag) && is.finite(dense_diag$max_rhat) &&
              dense_diag$max_rhat < 1.05 &&
              is.finite(dense_diag$min_ess_bulk) && dense_diag$min_ess_bulk >= 100
  list(curvature_cond = { ev <- eigen(Jpsd, symmetric = TRUE, only.values = TRUE)$values
                          max(ev) / min(ev) },
       tune = tune, best_eps = best$eps,
       dense_euclidean = dense_full, dense_chains = dense_chains,
       dense_diag = dense_diag, sub_riemannian = sub,
       sub_meta = sub_meta, identity_euclidean = euc_id,
       resolved_components = list(acc = acc_ok, div = div_ok, ebfmi = ebfmi_ok,
                                  convergence = conv_ok),
       resolved = acc_ok && div_ok && ebfmi_ok && conv_ok)
}

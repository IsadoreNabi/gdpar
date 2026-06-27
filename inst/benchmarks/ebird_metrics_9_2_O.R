## Block 9.2.O: predictive-metric module for the eBird external
## re-validation (DESIGN_9_2_EXTERNAL section 4). Unlike 9.2.S there is
## NO ground truth (model-vs-model S&T outcomes), so every metric is a
## held-out PREDICTIVE quantity on the test set.
##
## THE primary metric is the held-out test log predictive density on the
## ORIGINAL response measure (occurrence in [0,1], count in [0, Inf)).
## This is the single comparability countermeasure (B9.19, criterion 2.O,
## three-layer rigor / statistical layer): a Gaussian fit on a TRANSFORMED
## response (logit-occurrence / log1p-count, "track G") and a natural
## support family (beta-occurrence / tweedie-count, "track N") live on
## DIFFERENT base measures, so their raw log densities are NOT comparable.
## Reporting every method's log density on the same original measure via
## the change-of-variables Jacobian makes (a) all methods within a track
## and (b) the two tracks against each other directly comparable.
##
## Jacobian (per observation, method- and draw-independent):
##   coord 1 occurrence, track G: z1 = logit(y1c); log|dz1/dy1| =
##     -log(y1c) - log1p(-y1c)   [y1c = y1 clipped to (eps, 1-eps)]
##   coord 2 count,      track G: z2 = log1p(y2);  log|dz2/dy2| =
##     -log1p(y2)
##   track N (beta / tweedie): density already on the original measure ->
##     Jacobian 0.
## The Jacobian offset is a constant added equally to every track-G method
## (cancels in within-track-G ranking) but is essential to compare track G
## against track N on the common original measure.
##
## Each adapter (ebird_adapters_9_2_O.R) returns, per coordinate, EITHER
##   (a) $param_draws: named list of S x n_test posterior draw matrices of
##       the predictive-distribution parameters on the family's native
##       scale (mu on the linear-predictor inverse-link scale, etc.), with
##       a $family_tag in {gaussian_z, student_t_z, beta, tweedie}; OR
##   (b) $loglik_native: an S x n_test matrix of log p(y_test | theta) on
##       the family's NATIVE response scale (brms/rstanarm native log_lik),
##       plus $family_tag to decide the Jacobian.
## ebird_coord_elpd() turns either into the elpd on the original measure.

`%||%` <- function(a, b) if (is.null(a)) b else a

.ebird_eps_logit <- 1e-4

## --- numerically stable per-point log-mean-exp over draws ------------
.ebird_logmeanexp_col <- function(col) {
  good <- col[is.finite(col)]
  if (length(good) == 0L) return(-Inf)
  m <- max(good)
  m + log(mean(exp(good - m)))
}

## elpd = sum over test points of logmeanexp_s log p(y_i | theta_s).
## Returns list(elpd, n_pts, n_nonfinite_pts).
ebird_elpd_from_loglik <- function(loglik) {
  if (is.null(loglik) || !is.matrix(loglik)) {
    return(list(elpd = NA_real_, n_pts = NA_integer_,
                n_nonfinite_pts = NA_integer_))
  }
  pw <- apply(loglik, 2L, .ebird_logmeanexp_col)
  nbad <- sum(!is.finite(pw))
  finite_pw <- pw[is.finite(pw)]
  list(elpd = if (length(finite_pw) > 0L) sum(finite_pw) else NA_real_,
       n_pts = ncol(loglik), n_nonfinite_pts = as.integer(nbad))
}

## --- change-of-variables Jacobian on the original measure -----------
## coord in {1, 2}; track in {"G", "N"}. y_orig is the observed test
## outcome on the ORIGINAL scale (occurrence in [0,1] / count in [0,Inf)).
ebird_logjac <- function(coord, track, y_orig) {
  if (identical(track, "N")) return(rep(0, length(y_orig)))
  if (coord == 1L) {
    y1c <- pmin(pmax(y_orig, .ebird_eps_logit), 1 - .ebird_eps_logit)
    return(-log(y1c) - log1p(-y1c))
  }
  -log1p(y_orig)             ## coord 2, track G: log1p
}

## Transform an original-scale outcome to the track-G modelling scale.
ebird_to_z <- function(coord, y_orig) {
  if (coord == 1L) {
    return(stats::qlogis(pmin(pmax(y_orig, .ebird_eps_logit),
                              1 - .ebird_eps_logit)))
  }
  log1p(y_orig)
}

## --- family log densities (S x n) on the family's native scale ------
## params: named list of S x n matrices (or length-S vectors recycled to
## columns). Returns S x n log density. The caller adds the Jacobian.

.ebird_recycle <- function(v, S, n) {
  if (is.matrix(v)) return(v)
  if (length(v) == 1L) return(matrix(v, S, n))
  if (length(v) == S)  return(matrix(v, S, n, byrow = FALSE))
  if (length(v) == n)  return(matrix(v, S, n, byrow = TRUE))
  stop("ebird: cannot recycle parameter to S x n")
}

## Gaussian density of z (the transformed response) given mu, sigma draws.
.ebird_logdens_gaussian_z <- function(z, mu, sigma) {
  S <- nrow(mu); n <- ncol(mu)
  sig <- .ebird_recycle(sigma, S, n)
  zmat <- matrix(z, S, n, byrow = TRUE)
  stats::dnorm(zmat, mean = mu, sd = pmax(sig, 1e-8), log = TRUE)
}

## Scaled Student-t density of z given mu, sigma, nu draws.
.ebird_logdens_student_z <- function(z, mu, sigma, nu) {
  S <- nrow(mu); n <- ncol(mu)
  sig <- pmax(.ebird_recycle(sigma, S, n), 1e-8)
  nuu <- pmax(.ebird_recycle(nu, S, n), 2.01)
  zmat <- matrix(z, S, n, byrow = TRUE)
  stats::dt((zmat - mu) / sig, df = nuu, log = TRUE) - log(sig)
}

## Beta density on (0,1) with mean-precision parameterization mu, phi.
.ebird_logdens_beta <- function(y, mu, phi) {
  S <- nrow(mu); n <- ncol(mu)
  ph <- pmax(.ebird_recycle(phi, S, n), 1e-8)
  muc <- pmin(pmax(mu, 1e-6), 1 - 1e-6)
  yc <- pmin(pmax(y, .ebird_eps_logit), 1 - .ebird_eps_logit)
  ymat <- matrix(yc, S, n, byrow = TRUE)
  stats::dbeta(ymat, shape1 = muc * ph, shape2 = (1 - muc) * ph, log = TRUE)
}

## Tweedie density on [0, Inf) with mu, phi, power p draws. tweedie::dtweedie
## requires a SCALAR power; p is a population parameter (p ~ 1, negligible
## posterior spread), so we use the per-column posterior-mean power as the
## scalar and vectorize over the S draws of mu/phi (documented approximation,
## three-layer rigor: algebraic density vs the tweedie-package scalar-power
## constraint).
.ebird_logdens_tweedie <- function(y, mu, phi, power) {
  if (!requireNamespace("tweedie", quietly = TRUE)) {
    stop("Package 'tweedie' is required for the tweedie density.")
  }
  S <- nrow(mu); n <- ncol(mu)
  ph <- pmax(.ebird_recycle(phi, S, n), 1e-8)
  pw <- pmin(pmax(.ebird_recycle(power, S, n), 1.01), 1.99)
  out <- matrix(-Inf, S, n)
  for (i in seq_len(n)) {
    yi <- max(y[i], 0)
    p_i <- mean(pw[, i])
    out[, i] <- log(pmax(tweedie::dtweedie(rep(yi, S), mu = pmax(mu[, i], 1e-8),
                                           phi = ph[, i], power = p_i),
                         1e-300))
  }
  out
}

## Dispatch: log density on ORIGINAL measure for a param-draw adapter.
## family_tag in {gaussian_z, student_t_z, beta, tweedie}.
ebird_loglik_orig_from_params <- function(family_tag, param_draws,
                                          y_orig, coord, track) {
  pd <- param_draws
  ll_native <- switch(family_tag,
    gaussian_z  = .ebird_logdens_gaussian_z(ebird_to_z(coord, y_orig),
                                            pd$mu, pd$sigma),
    student_t_z = .ebird_logdens_student_z(ebird_to_z(coord, y_orig),
                                           pd$mu, pd$sigma, pd$nu),
    beta        = .ebird_logdens_beta(y_orig, pd$mu, pd$phi),
    tweedie     = .ebird_logdens_tweedie(y_orig, pd$mu, pd$phi, pd$p),
    stop(sprintf("ebird: unknown family_tag '%s'", family_tag)))
  jac <- ebird_logjac(coord, track, y_orig)
  sweep(ll_native, 2L, jac, `+`)
}

## Dispatch for a native-log_lik adapter (brms/rstanarm): add the Jacobian.
ebird_loglik_orig_from_native <- function(loglik_native, y_orig,
                                          coord, track) {
  jac <- ebird_logjac(coord, track, y_orig)
  sweep(loglik_native, 2L, jac, `+`)
}

## --- in-sample PSIS-LOO on the ORIGINAL measure (track B) -----------
## Secondary cross-check (DESIGN 4.2). The native loo object gives elpd on
## the method's modelling response scale over the TRAINING points; adding
## the per-training-point Jacobian offset (constant across methods within a
## track) puts it on the original measure, comparable across tracks. Only
## the fully-Bayesian per-coordinate fits supply a loo object (gdpar/brms/
## rstanarm); mgcv/INLA loo is skipped for cost (DESIGN 2.4), as in 9.2.S.
ebird_loo_orig <- function(loo_obj, y_train_orig, coord, track) {
  if (is.null(loo_obj)) {
    return(list(elpd_loo = NA_real_, pareto_k_max = NA_real_))
  }
  est <- tryCatch(loo_obj$estimates["elpd_loo", "Estimate"],
                  error = function(e) NA_real_)
  jac <- sum(ebird_logjac(coord, track, y_train_orig))
  pk <- tryCatch(loo_obj$diagnostics$pareto_k, error = function(e) NULL)
  list(elpd_loo = as.numeric(est) + jac,
       pareto_k_max = if (!is.null(pk)) max(pk, na.rm = TRUE) else NA_real_)
}

## --- coverage of the 90% posterior predictive interval --------------
## ypred_draws: S x n predictive draws on the ORIGINAL scale. Returns the
## empirical fraction of test outcomes inside [q05, q95].
ebird_coverage90 <- function(ypred_draws, y_orig) {
  if (is.null(ypred_draws) || !is.matrix(ypred_draws)) return(NA_real_)
  qs <- apply(ypred_draws, 2L, stats::quantile, probs = c(0.05, 0.95),
              names = FALSE, na.rm = TRUE)
  mean(y_orig >= qs[1L, ] & y_orig <= qs[2L, ], na.rm = TRUE)
}

## --- posterior predictive simulation on the ORIGINAL scale ----------
## Used to compute coverage uniformly from param draws (one y draw per
## param draw). For track G the simulation is on z then back-transformed.
ebird_sim_ypred_orig <- function(family_tag, param_draws, coord, track,
                                 seed = 1L) {
  set.seed(seed)
  pd <- param_draws
  S <- nrow(pd$mu); n <- ncol(pd$mu)
  z <- switch(family_tag,
    gaussian_z  = matrix(stats::rnorm(S * n, pd$mu,
                                      .ebird_recycle(pd$sigma, S, n)), S, n),
    student_t_z = pd$mu + .ebird_recycle(pd$sigma, S, n) *
                  matrix(stats::rt(S * n, df = as.vector(
                    pmax(.ebird_recycle(pd$nu, S, n), 2.01))), S, n),
    beta        = {
      ph <- .ebird_recycle(pd$phi, S, n); muc <- pmin(pmax(pd$mu, 1e-6), 1 - 1e-6)
      matrix(stats::rbeta(S * n, muc * ph, (1 - muc) * ph), S, n)
    },
    tweedie     = {
      ph <- .ebird_recycle(pd$phi, S, n); pw <- .ebird_recycle(pd$p, S, n)
      m <- matrix(0, S, n)
      for (j in seq_len(n)) {
        m[, j] <- tweedie::rtweedie(S, mu = pmax(pd$mu[, j], 1e-8),
                                    phi = ph[, j],
                                    power = mean(pmin(pmax(pw[, j], 1.01), 1.99)))
      }
      m
    },
    stop(sprintf("ebird: unknown family_tag '%s'", family_tag)))
  if (identical(family_tag, "beta") || identical(family_tag, "tweedie")) {
    return(z)                                  ## already original scale
  }
  if (coord == 1L) stats::plogis(z) else expm1(z)
}

# Automated Laplace fallback for the geometry-adaptive orchestrator (Block RG,
# RG.7 step 4, session B9.35; canonised D99).
#
# RG.7 certified that the near-deterministic Tweedie canyon of benchmark 9.2.O
# (the `quasi_deterministic` pathology; observed-information condition number
# ~9e7-4.4e8 at the mode; unconstrained dimension 14) is GENUINELY non-Gaussian
# and not sampleable at feasible cost: default NUTS, a constant-mass Laplace HMC
# and a whitened NUTS (dense metric = the Laplace precision, every chain
# initialised at the mode) all FAIL the multi-chain R-hat / ESS gate -- the
# divergences persist in the already-whitened coordinates, which can only mean
# the canyon curves. The honest endpoint (both lineages and the user agreed) is
# the LAPLACE / plug-in predictive: exactly the regime of mgcv/REML and
# INLA/Laplace, accompanied by a fidelity diagnostic. The stopping rule is "no
# more samplers".
#
# This file PROMOTES the Laplace machinery validated in the benchmark extractor
# (inst/benchmarks/scripts/rg7_laplace_elpd.R) into R/ as a first-class,
# exported capability, gdpar_geom_laplace(), and wires it into the orchestrator
# (R/geometry_orchestrator.R) and the one-call fit (R/geometry_bridge.R) behind
# an opt-in laplace_fallback = FALSE flag: the default fit path is byte-for-byte
# unchanged (the orchestrator output carries no $laplace and the status set is
# unchanged), so the goldens stay bit-identical. When opted in and a run ends in
# a certified limit (the non-Gaussian canyon), the orchestrator attaches the
# Laplace approximation and relabels the status `certified_limit_laplace`.
#
# Honesty (charter section 5; ORPHEUS section 16.3; [[feedback_no_overreach]]).
# The Laplace approximation is NEVER advertised as exact. It carries its own
# fidelity diagnostics against the true posterior over the SAME draws -- the
# self-normalised importance-sampling effective sample size, the PSIS Pareto-k of
# the weights p/q, and the mean log-density drop log p(mode) - log p(theta)
# against its Gaussian expectation d/2 -- distilled into a single, unambiguous
# scalar label (good / poor / very_poor) so a caller can never mis-read raw
# diagnostics. On the real 9.2.O canyon this label is `very_poor`, which is the
# scientific point: a measurement of HOW non-Gaussian the canyon is.
#
# Design feedback adjudicated with MiMo-2.5-Pro (B9.34, all accepted; recorded in
# DESIGN_RG7_LAPLACE_FALLBACK.md): hook only on the certified-limit path (and the
# out-of-scope path when the curvature is genuinely positive-definite), never
# pre-empting the sub-Riemannian climb (the climb is the EVIDENCE for the
# certificate); expose the flag on gdpar_geom_fit() too; default draws = 0L (the
# mode + precision Gaussian IS the approximation); add the scalar fit-quality
# label; store the UN-floored condition number and warn when the eigen-floor
# could be hiding ill-conditioning; keep the mode-climb reading the same target
# and gradient as the sampler.
#
# Anchors (gdpar integrates established patterns, it does not claim to invent
# them): Tierney & Kadane 1986 (the Laplace approximation to posterior
# expectations); Vehtari, Gelman & Gabry 2017 and Vehtari et al. 2024 (the PSIS
# Pareto-k reliability index and ELPD); Bernstein--von Mises (the local-Gaussian
# regime, and where it breaks down); Rue, Martino & Chopin 2009 (INLA, the
# nested-Laplace predictive that this fallback matches).

# ---------------------------------------------------------------------------
# Internal helpers (the promoted extractor core).
# ---------------------------------------------------------------------------

# Modified-Newton climb to the posterior mode using the engine gradient and the
# exact Hessian (the promotion of rg7_climb_to_mode). At each step it forms
# A = -Hessian, floors A's eigenvalues to a positive minimum so an indefinite
# Hessian (a saddle / not-yet-mode) still yields an ascent direction A^{-1} g,
# caps the step in a trust region and backtracks on log p. When the Hessian is
# non-finite in a stiff off-mode region it falls back to a steepest-ascent step
# so the climb can TRAVERSE that region to the benign mode (eigen() is never fed
# a non-finite matrix). Requires geom_target$hessian; the caller checks first.
.gdpar_geom_newton_climb <- function(geom_target, ref, n_steps = 300L,
                                     ridge = 1e-6, tol = 1e-3,
                                     max_step_norm = 5) {
  lp <- geom_target$log_prob
  gl <- geom_target$grad_log_prob
  he <- geom_target$hessian
  if (!is.function(he)) {
    return(list(mode = ref, converged = FALSE, reason = "no exact Hessian",
                grad_norm = sqrt(sum(gl(ref)^2)), logp = as.numeric(lp(ref))))
  }
  th <- ref
  f0 <- as.numeric(lp(th))
  for (it in seq_len(n_steps)) {
    g <- gl(th)
    H <- he(th)
    if (any(!is.finite(g))) break
    if (all(is.finite(H))) {
      H <- 0.5 * (H + t(H))
      e <- eigen(-H, symmetric = TRUE)
      lam <- pmax(e$values, ridge * max(abs(e$values), 1))
      step <- as.numeric(e$vectors %*% ((crossprod(e$vectors, g)) / lam))
    } else {
      step <- g                                   # steepest-ascent fallback.
    }
    if (any(!is.finite(step))) break
    sn <- sqrt(sum(step^2))
    if (sn > max_step_norm) step <- step * (max_step_norm / sn)
    t <- 1
    ok <- FALSE
    for (bt in seq_len(20L)) {
      cand <- th + t * step
      fc <- tryCatch(as.numeric(lp(cand)), error = function(e) -Inf)
      if (is.finite(fc) && fc >= f0 + 1e-4 * t * sum(g * step)) {
        ok <- TRUE
        break
      }
      t <- t / 2
    }
    if (!ok) break
    gnorm <- sqrt(sum(g^2))
    improve <- fc - f0
    th <- cand
    f0 <- fc
    if (gnorm < tol || improve < tol * 1e-2) {
      return(list(mode = th, converged = TRUE, iters = it, grad_norm = gnorm,
                  logp = f0))
    }
  }
  list(mode = th, converged = FALSE, logp = f0,
       grad_norm = sqrt(sum(gl(th)^2)))
}

# Reach the mode from a reference: an L-BFGS-B warm start on -log p (robust when
# the reference is far) followed by the modified-Newton polish when an exact
# Hessian is available; L-BFGS alone otherwise. Both stages read the SAME target
# and gradient the sampler uses (the geom_target).
.gdpar_geom_laplace_climb <- function(geom_target, start, climb_steps) {
  neg_lp <- function(th) {
    v <- tryCatch(-as.numeric(geom_target$log_prob(th)), error = function(e) Inf)
    if (!is.finite(v)) 1e18 else v
  }
  neg_gr <- function(th) {
    g <- tryCatch(-geom_target$grad_log_prob(th),
                  error = function(e) rep(0, length(th)))
    g[!is.finite(g)] <- 0
    g
  }
  opt <- tryCatch(
    stats::optim(start, neg_lp, neg_gr, method = "L-BFGS-B",
                 control = list(maxit = 5000L, factr = 1e1)),
    error = function(e) NULL)
  warm <- if (!is.null(opt)) opt$par else start
  if (is.function(geom_target$hessian)) {
    return(.gdpar_geom_newton_climb(geom_target, warm, n_steps = climb_steps))
  }
  list(mode = warm,
       converged = !is.null(opt) && isTRUE(opt$convergence == 0L),
       grad_norm = sqrt(sum(geom_target$grad_log_prob(warm)^2)),
       logp = as.numeric(geom_target$log_prob(warm)))
}

# Observed information = -Hessian(log p), symmetrised (the promotion of
# rg7_observed_information). Prefers the exact cmdstan Hessian; on a
# near-deterministic posterior a fixed finite-difference step cannot resolve the
# stiff wall directions, so finite differences are only the fallback when no
# Hessian is exposed. Records the route in attr(., "method").
.gdpar_geom_observed_information <- function(geom_target, reference, h = 1e-4) {
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
  obs <- -0.5 * (J + t(J))
  attr(obs, "method") <- if (all(is.finite(obs))) "finite_difference" else
    "finite_difference_nonfinite"
  obs
}

# iid Laplace draws on the unconstrained scale: theta_s = mode + M^{-1/2} z_s,
# with cov(theta_s) = M^{-1}. Uses a LOCAL RNG (saves + restores the global seed)
# so the call leaves the caller's RNG stream untouched.
.gdpar_geom_laplace_draws_unconstrained <- function(mode, Lhalf, S, seed) {
  old_seed <- if (exists(".Random.seed", envir = globalenv()))
    get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (!is.null(old_seed))
    assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
  set.seed(seed)
  d <- length(mode)
  Z <- matrix(stats::rnorm(S * d), nrow = S, ncol = d)
  U <- Z %*% Lhalf                              # Lhalf symmetric => row cov = M^{-1}.
  sweep(U, 2L, mode, "+")
}

# Fit quality of the Laplace Gaussian q against the TRUE posterior p over the
# SAME draws (the promotion of rg7_laplace_fit_quality). The per-point spread of
# a downstream predictive does NOT bound the Laplace-to-exact gap; the right
# diagnostics are (i) the self-normalised IS effective sample size with
# log-weights log p - log q, plus its PSIS Pareto-k tail index, and (ii) the
# mean / max log-density drop log p(mode) - log p(theta), which is ~ chi^2_d / 2
# (mean d/2) under a perfect Gaussian. Cheap: one log_prob per draw.
.gdpar_geom_laplace_fit_quality <- function(geom_target, mode, M, logdetM, U) {
  S <- nrow(U)
  d <- length(mode)
  cen <- sweep(U, 2L, mode, "-")
  quad <- rowSums((cen %*% M) * cen)
  log_q <- -0.5 * d * log(2 * pi) + 0.5 * logdetM - 0.5 * quad
  log_p <- vapply(seq_len(S), function(s)
    tryCatch(as.numeric(geom_target$log_prob(U[s, ])),
             error = function(e) NA_real_), numeric(1L))
  logp_mode <- tryCatch(as.numeric(geom_target$log_prob(mode)),
                        error = function(e) NA_real_)
  fin <- is.finite(log_p) & is.finite(log_q)
  if (sum(fin) < 2L) {
    return(list(ess_is = NA_real_, ess_is_frac = NA_real_, pareto_k = NA_real_,
                logdrop_mean = NA_real_, logdrop_max = NA_real_,
                logdrop_expected = d / 2, n_finite = sum(fin), n_total = S))
  }
  lw <- log_p[fin] - log_q[fin]
  lw <- lw - max(lw)
  w <- exp(lw)
  ess <- (sum(w)^2) / sum(w^2)
  pk <- tryCatch(
    if (requireNamespace("loo", quietly = TRUE) && sum(fin) > 10L)
      suppressWarnings(loo::psis(matrix(log_p[fin] - log_q[fin], ncol = 1L),
                                 r_eff = NA)$diagnostics$pareto_k)
    else NA_real_,
    error = function(e) NA_real_)
  drop <- logp_mode - log_p[fin]
  list(ess_is = ess, ess_is_frac = ess / sum(fin), pareto_k = as.numeric(pk),
       logdrop_mean = mean(drop), logdrop_max = max(drop),
       logdrop_expected = d / 2, n_finite = sum(fin), n_total = S)
}

# Distil the fidelity diagnostics into a single, unambiguous scalar label so a
# caller can never mis-read raw numbers (MiMo design feedback, accepted B9.34).
# Data-driven thresholds ([[feedback_data_driven_knobs]]): the IS-ESS fraction is
# the STABLE primary signal (no tail fit) -- floor 0.5 says the Gaussian proposal
# retains half its nominal information, and a near-perfect ess ~ 1 with a benign
# log-density drop is `good` regardless of the Pareto-k (which, on numerically
# constant importance weights, is meaningless); PSIS Pareto-k < 0.7 is the
# reliability cut of Vehtari et al., used to refine the middle ground and damning
# only when the ESS is also weak; the mean log-density drop is d/2 under a perfect
# Gaussian, so <= d is benign and >= 2d is severe. A non-positive-definite
# curvature at the "mode" (a saddle) is `very_poor` outright: the Laplace Gaussian
# is ill-defined there.
.gdpar_geom_laplace_label <- function(fq, all_pos, d, mode_offset_sd = 0) {
  if (!isTRUE(all_pos)) return("very_poor")
  pk <- fq$pareto_k
  ess <- fq$ess_is_frac
  ld <- fq$logdrop_mean
  if (!is.finite(ess) || !is.finite(ld)) return("very_poor")
  # Mode-location certification (GLM design review, B9.35): a Laplace centred on a
  # point that is NOT the mode is confidently wrong. The Newton decrement
  # mode_offset_sd bounds the offset in posterior SDs; a large offset is
  # `very_poor` outright, and a "good" label requires a well-centred mode.
  off <- if (is.finite(mode_offset_sd)) mode_offset_sd else Inf
  if (off > 2) return("very_poor")
  well_centred <- off <= 0.5
  # ESS is the STABLE primary signal (no tail fit). A near-perfect importance
  # fit (ess ~ 1, log-drop ~ d/2) means q matches p essentially exactly; the PSIS
  # Pareto-k is then estimated on (numerically) constant weights and is NOT
  # informative, so it must not override a clearly excellent fit. Pareto-k only
  # refines the middle ground, and is damning only when the ESS is also weak.
  if (well_centred && ess >= 0.9 && ld <= d) return("good")
  pk_bad  <- isTRUE(pk >= 1.0)
  if (ess < 0.1 || ld >= 2 * d || (pk_bad && ess < 0.5)) return("very_poor")
  pk_ok <- is.na(pk) || pk < 0.7
  if (well_centred && ess >= 0.5 && ld <= d && pk_ok) return("good")
  "poor"
}

# ---------------------------------------------------------------------------
# gdpar_geom_laplace(): the exported Laplace approximation.
# ---------------------------------------------------------------------------

#' Laplace approximation of a posterior at its mode
#'
#' Compute the Laplace (mode + curvature) Gaussian approximation
#' \eqn{N(\hat\theta, M^{-1})} of a posterior exposed as a
#' \code{\link{gdpar_geom_target}}, on the unconstrained scale: climb to the mode
#' \eqn{\hat\theta}, form the precision \eqn{M = -\nabla^2 \log p(\hat\theta)}
#' (the observed information), and report its covariance, optional draws, and a
#' \strong{fidelity diagnostic} of the Gaussian against the true posterior. This
#' is the first-class, honest endpoint for a posterior the geometry-adaptive
#' orchestrator certifies as a non-sampleable, genuinely non-Gaussian canyon
#' (Block RG, RG.7): exactly the regime of \pkg{mgcv}/REML and \pkg{INLA}/Laplace
#' competitors, accompanied by a measurement of how good the Gaussian is.
#'
#' @details
#' The mode is reached by an L-BFGS-B warm start on \eqn{-\log p} followed, when
#' the target exposes an exact Hessian, by a modified-Newton polish; both stages
#' read the \emph{same} target and gradient the sampler uses. The precision
#' \eqn{M} is the symmetrised \eqn{-}Hessian (the exact cmdstan Hessian when
#' available, central finite differences otherwise), eigen-floored to be
#' positive-definite so the draw machinery stays numerically alive; the
#' \strong{un-floored} condition number is reported separately and a warning is
#' raised when the floor could be masking ill-conditioning. A non-positive-definite
#' raw curvature (a saddle, not a maximum) is flagged loudly.
#'
#' \strong{Fidelity, never blind trust.} The Gaussian \eqn{q} is scored against
#' the true posterior \eqn{p} over the same draws: the self-normalised
#' importance-sampling effective sample size with log-weights
#' \eqn{\log p - \log q}, its PSIS Pareto-\eqn{k} tail index (when \pkg{loo} is
#' installed), and the mean / max log-density drop
#' \eqn{\log p(\hat\theta) - \log p(\theta)} against its Gaussian expectation
#' \eqn{d/2}. These are distilled into a single scalar label,
#' \code{"good"} / \code{"poor"} / \code{"very_poor"}, so the approximation is
#' never mistaken for exact MCMC. On a curved canyon (the real 9.2.O Tweedie
#' count) the label is \code{"very_poor"} -- the scientific finding, not a defect.
#'
#' @param geom_target A \code{\link{gdpar_geom_target}} (or an object coercible to
#'   one) exposing \code{log_prob}, \code{grad_log_prob}, \code{dim} and,
#'   ideally, \code{hessian}.
#' @param reference Optional unconstrained warm-start position for the mode climb;
#'   defaults to the origin. When called by the orchestrator this is the best
#'   on-ridge position found during sampling.
#' @param draws Number of iid Laplace draws \eqn{\hat\theta + M^{-1/2} z} to
#'   return (default \code{0L}: the mode plus precision Gaussian \emph{is} the
#'   approximation; the caller asks for draws explicitly when it needs them for a
#'   downstream predictive).
#' @param climb Logical; whether to climb to the mode from \code{reference}
#'   (default \code{TRUE}). \code{FALSE} treats \code{reference} as the mode and
#'   only reads the curvature there.
#' @param seed Integer seed for the (local, stream-preserving) draw RNG; defaults
#'   to a fixed value so the result is reproducible.
#' @param fit_quality_draws Number of internal iid draws used to assess fidelity
#'   when \code{draws} is small (default \code{256L}); the fidelity label is
#'   always computed, independently of how many user-facing \code{draws} are
#'   requested.
#' @param eigen_floor_rel Relative eigenvalue floor applied to \eqn{M} for
#'   positive-definiteness (default \code{1e-10} of the largest eigenvalue).
#' @param climb_steps Maximum modified-Newton steps in the mode climb (default
#'   \code{300L}).
#' @param cond_warn Condition-number threshold above which the un-floored
#'   curvature triggers an ill-conditioning warning (default \code{1e12}).
#'
#' @return An object of class \code{gdpar_geom_laplace}: a list with \code{mode},
#'   the precision \code{M}, the covariance \code{cov} (\eqn{M^{-1}}), the
#'   symmetric square root \code{Lhalf} (\eqn{M^{-1/2}}), \code{logdet} (of
#'   \eqn{M}), the eigenvalues \code{eig}, the floored and un-floored condition
#'   numbers \code{cond} / \code{cond_unfloored}, the count \code{n_floored} of
#'   eigen-floored directions and the \code{floor_value} used, \code{all_pos}
#'   (was the raw curvature positive-definite), \code{mode_offset_sd} (the
#'   Newton-decrement bound on the mode offset in posterior SDs),
#'   \code{grad_norm}, \code{logp}, \code{converged}, \code{method} (the Hessian
#'   route), \code{dim}, \code{draws} (a \code{draws} \eqn{\times d} matrix,
#'   possibly zero rows, carrying \code{attr(., "approximation") = "laplace"} so
#'   it is never mistaken for exact MCMC draws downstream),
#'   \code{fit_quality} (the fidelity diagnostics) and
#'   \code{fit_quality_label}.
#'
#' @seealso \code{\link{gdpar_geom_orchestrate}} (the \code{laplace_fallback}
#'   opt-in), \code{\link{gdpar_geom_fit}}, \code{\link{gdpar_geom_target}}.
#'
#' @examples
#' # A correlated Gaussian: the Laplace approximation is exact, so the fidelity
#' # label is "good". The target carries an analytic (constant) Hessian.
#' A <- matrix(c(2, 0.8, 0.8, 1), 2, 2)
#' mu <- c(1, -0.5)
#' tgt <- gdpar_geom_target(
#'   log_prob = function(th) -0.5 * as.numeric(t(th - mu) %*% A %*% (th - mu)),
#'   grad_log_prob = function(th) -as.numeric(A %*% (th - mu)),
#'   hessian = function(th) -A, dim = 2L)
#' lap <- gdpar_geom_laplace(tgt, draws = 200L, seed = 1L)
#' lap$fit_quality_label
#' max(abs(lap$mode - mu)) < 1e-6
#'
#' @export
gdpar_geom_laplace <- function(geom_target, reference = NULL, draws = 0L,
                               climb = TRUE, seed = NULL,
                               fit_quality_draws = 256L,
                               eigen_floor_rel = 1e-10, climb_steps = 300L,
                               cond_warn = 1e12) {
  if (!inherits(geom_target, "gdpar_geom_target")) {
    geom_target <- gdpar_geom_target(geom_target)
  }
  d <- as.integer(geom_target$dim)
  draws <- as.integer(draws)
  if (is.na(draws) || draws < 0L) {
    gdpar_abort("'draws' must be a non-negative integer.",
                class = "gdpar_input_error")
  }
  fit_quality_draws <- max(as.integer(fit_quality_draws), 2L)
  if (is.null(seed)) seed <- 20260625L
  assert_count(seed, "seed")
  if (is.null(reference)) {
    reference <- rep(0, d)
  } else {
    reference <- as.numeric(reference)
    if (length(reference) != d) {
      gdpar_abort(sprintf("'reference' must have length %d (the unconstrained dimension).", d),
                  class = "gdpar_input_error")
    }
  }

  # (1) mode.
  if (isTRUE(climb)) {
    cl <- .gdpar_geom_laplace_climb(geom_target, reference, climb_steps)
  } else {
    cl <- list(mode = reference, converged = NA,
               grad_norm = sqrt(sum(geom_target$grad_log_prob(reference)^2)),
               logp = as.numeric(geom_target$log_prob(reference)))
  }
  mode <- as.numeric(cl$mode)
  logp <- as.numeric(cl$logp %||% geom_target$log_prob(mode))

  # (2) precision M = -Hessian, eigen-floored to PD; covariance + sqrt + logdet.
  M <- .gdpar_geom_observed_information(geom_target, mode)
  M <- 0.5 * (M + t(M))
  eM <- eigen(M, symmetric = TRUE)
  all_pos <- all(eM$values > 0)
  if (!all_pos) {
    gdpar_warn(sprintf(paste(
      "gdpar_geom_laplace: -Hessian at the mode has %d non-positive eigenvalue(s)",
      "(min = %.3g) -- possible saddle, not a maximum; the Laplace Gaussian is",
      "flagged 'very_poor'."), sum(eM$values <= 0), min(eM$values)),
      class = "gdpar_geometry_warning")
  }
  cond_unfloored <- {
    av <- abs(eM$values)
    if (min(av) > 0) max(av) / min(av) else Inf
  }
  if (is.finite(cond_warn) && cond_unfloored > cond_warn) {
    gdpar_warn(sprintf(paste(
      "gdpar_geom_laplace: the un-floored curvature condition number is %.3g",
      "(> %.3g); the eigen-floor keeps the covariance defined but may be masking",
      "ill-conditioning -- treat the stiffest directions with care."),
      cond_unfloored, cond_warn), class = "gdpar_geometry_warning")
  }
  floor_eig <- max(eM$values) * eigen_floor_rel
  lam <- pmax(eM$values, floor_eig)
  # The eigen-floor IS part of the approximation (GLM design review, B9.35): a
  # floored eigenvalue means a huge posterior variance the Gaussian invents in
  # that direction. Report how many directions were floored and the floor value
  # so a silent floor is never a silent prior.
  n_floored <- sum(eM$values < floor_eig)
  logdetM <- sum(log(lam))
  Hinv  <- eM$vectors %*% (t(eM$vectors) / lam)          # M^{-1}.
  Lhalf <- eM$vectors %*% (t(eM$vectors) / sqrt(lam))    # M^{-1/2} (symmetric).
  g_mode <- geom_target$grad_log_prob(mode)
  newton_step <- as.numeric(Hinv %*% g_mode)
  mode_offset_sd <- sqrt(max(sum(g_mode * newton_step), 0))   # Newton decrement.

  # (3) draws + fidelity (the label is always computed on an internal sample).
  S_fit <- max(draws, fit_quality_draws)
  U_full <- .gdpar_geom_laplace_draws_unconstrained(mode, Lhalf, S_fit, seed)
  fq <- .gdpar_geom_laplace_fit_quality(geom_target, mode, M, logdetM, U_full)
  label <- .gdpar_geom_laplace_label(fq, all_pos, d, mode_offset_sd)
  out_draws <- if (draws > 0L) U_full[seq_len(draws), , drop = FALSE] else
    U_full[integer(0), , drop = FALSE]
  # Tag the draws so they can never be silently consumed as exact MCMC draws by a
  # downstream method (GLM design review, B9.35): the attribute travels with the
  # matrix even when extracted from the object.
  attr(out_draws, "approximation") <- "laplace"

  obj <- list(
    mode = mode, M = M, cov = Hinv, Lhalf = Lhalf, logdet = logdetM,
    eig = eM$values, cond = max(lam) / min(lam),
    cond_unfloored = cond_unfloored, n_floored = n_floored,
    floor_value = floor_eig, all_pos = all_pos,
    mode_offset_sd = mode_offset_sd, grad_norm = cl$grad_norm %||% NA_real_,
    logp = logp, converged = isTRUE(cl$converged),
    method = attr(M, "method") %||% "unknown", dim = d,
    draws = out_draws, fit_quality = fq, fit_quality_label = label)
  class(obj) <- c("gdpar_geom_laplace", "list")
  obj
}

#' Print method for gdpar_geom_laplace objects
#'
#' @param x A \code{gdpar_geom_laplace}.
#' @param ... Unused.
#' @return Invisibly returns \code{x}.
#' @export
print.gdpar_geom_laplace <- function(x, ...) {
  cat("<gdpar_geom_laplace> dim ", x$dim, " | fit-quality: ",
      x$fit_quality_label, "\n", sep = "")
  cat("  mode: |grad| ", format(x$grad_norm, digits = 3), " | offset ",
      format(x$mode_offset_sd, digits = 3), " SD | PD ", x$all_pos,
      " | method ", x$method, "\n", sep = "")
  cat("  condition: floored ", format(x$cond, digits = 3), " | un-floored ",
      format(x$cond_unfloored, digits = 3), " | eigen-floored dirs ",
      x$n_floored, "/", x$dim, "\n", sep = "")
  fq <- x$fit_quality
  is_unreliable <- isTRUE(fq$pareto_k >= 0.7)
  cat("  fidelity: IS-ESS ", format(fq$ess_is, digits = 3), " (",
      format(100 * fq$ess_is_frac, digits = 2), "%)",
      if (is_unreliable) " [unreliable: Pareto-k >= 0.7]" else "",
      " | Pareto-k ", format(fq$pareto_k, digits = 3), " | log-drop ",
      format(fq$logdrop_mean, digits = 3), " (E[d/2] = ", fq$logdrop_expected,
      ")\n", sep = "")
  if (!identical(x$fit_quality_label, "good")) {
    cat("  NOTE: a '", x$fit_quality_label, "' Laplace approximation -- a ",
        "plug-in / INLA-style fallback, NOT exact MCMC.\n", sep = "")
  }
  cat("  draws: ", nrow(x$draws), " x ", ncol(x$draws), "\n", sep = "")
  invisible(x)
}

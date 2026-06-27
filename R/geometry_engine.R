# R-native geometric sampling engine (Block RG, RG.2).
#
# Decision A (canonised D88): a geometric integrator written in pure R that
# delegates the log-density, its gradient and (for the Riemannian level of
# RG.3) its Hessian to a compiled backend -- a cmdstan model built with
# compile_model_methods = TRUE exposing $log_prob / $grad_log_prob / $hessian,
# or an R closure such as the dual targets of gdpar_geometry_suite(). The
# engine is the scaffolding the higher levels of the Block RG geometry
# hierarchy plug into: this file delivers the abstract pieces (a position-aware
# metric, a kinetic energy, a symplectic integrator) and the Euclidean level
# (constant metric + standard leapfrog + static Hamiltonian Monte Carlo),
# validated against the no-U-turn sampler of cmdstan on the easy targets. The
# Riemannian metric (Fisher / SoftAbs) and its implicit generalised leapfrog
# (Girolami & Calderhead 2011), and the Finsler / relativistic and
# sub-Riemannian levels, extend these same interfaces in RG.3 / RG.4.
#
# Conventions. The engine works on the UNCONSTRAINED scale, the scale on which
# Hamiltonian dynamics act and on which cmdstan's $grad_log_prob is defined.
# The potential energy is U(theta) = -log_prob(theta); the Hamiltonian is
# H(theta, p) = U(theta) + K(theta, p). With a Gaussian kinetic energy and
# mass matrix M, Hamilton's equations are dtheta/dt = M^{-1} p and
# dp/dt = grad_log_prob(theta) - grad_theta K (the second term vanishing for a
# constant, Euclidean metric). This is an opt-in, standalone capability; no
# existing fit path is touched, so the default branch stays bit-identical.

#' Build a target for the geometric sampling engine
#'
#' Wrap a posterior so the Block RG geometric engine
#' (\code{\link{gdpar_geom_hmc}}) can evaluate its log-density and gradient on
#' the unconstrained scale. This is the three-way adapter of decision A: the
#' integrator runs in R while the density work is delegated to a compiled
#' cmdstan model, to an R closure, or to a suite target.
#'
#' @param object Optional. Either a \pkg{cmdstanr} fit/model compiled with
#'   \code{compile_model_methods = TRUE} (exposing \code{$log_prob},
#'   \code{$grad_log_prob} and optionally \code{$hessian}), or a suite instance
#'   from \code{\link{gdpar_geometry_suite}}'s \code{make()} (a list carrying
#'   \code{log_prob}, \code{grad_log_prob} and \code{dim}). If supplied, the
#'   closure arguments are ignored.
#' @param log_prob,grad_log_prob Functions of an unconstrained parameter vector
#'   returning the log-density and its gradient. Used when \code{object} is
#'   \code{NULL}.
#' @param dim Integer dimension of the unconstrained parameter vector. Required
#'   for the closure form and for a cmdstan \code{object} (whose unconstrained
#'   dimension cannot always be inferred).
#' @param hessian Optional function returning the Hessian of \code{log_prob}
#'   (used by the Riemannian level of RG.3). For a cmdstan \code{object} the
#'   built-in \code{$hessian} is used when available.
#' @param param_names Optional character vector of parameter names.
#' @param data Optional data list when \code{object} is an uninstantiated
#'   cmdstan model rather than a fit (forwarded to a one-iteration fit used only
#'   to expose the methods).
#' @param simulate Optional generative function \code{simulate(theta)} returning
#'   one synthetic data set drawn from the model at \code{theta} (any object the
#'   companion \code{score} understands). Supplying \code{simulate} and
#'   \code{score} turns the target into a \emph{generative} target, the input the
#'   simulation-based expected-Fisher estimator
#'   \code{\link{gdpar_geom_fisher_simulator}} consumes.
#' @param score Optional function \code{score(theta, y)} returning the gradient
#'   of the model's log-\emph{likelihood} of a data set \code{y} at \code{theta},
#'   on the unconstrained scale (the per-data-set score whose outer product is
#'   the expected Fisher). Distinct from \code{grad_log_prob}, which is the
#'   gradient of the log-\emph{posterior} at the fixed observed data.
#'
#' @return A list of class \code{gdpar_geom_target} with elements
#'   \code{log_prob}, \code{grad_log_prob}, \code{hessian} (or \code{NULL}),
#'   \code{dim}, \code{param_names}, \code{backend}
#'   (\code{"closure"} or \code{"cmdstan"}), and the optional generative pieces
#'   \code{simulate} and \code{score}.
#'
#' @seealso \code{\link{gdpar_geom_hmc}}, \code{\link{gdpar_geom_metric_euclidean}}.
#'
#' @examples
#' # A two-dimensional standard normal via an R closure.
#' tgt <- gdpar_geom_target(
#'   log_prob = function(theta) -0.5 * sum(theta^2),
#'   grad_log_prob = function(theta) -theta,
#'   dim = 2
#' )
#' tgt$log_prob(c(0, 0))
#'
#' @export
gdpar_geom_target <- function(object = NULL, log_prob = NULL,
                              grad_log_prob = NULL, dim = NULL, hessian = NULL,
                              param_names = NULL, data = NULL,
                              simulate = NULL, score = NULL) {
  if (!is.null(simulate) && !is.function(simulate)) {
    gdpar_abort("'simulate' must be a function of theta or NULL.",
                class = "gdpar_input_error")
  }
  if (!is.null(score) && !is.function(score)) {
    gdpar_abort("'score' must be a function of (theta, y) or NULL.",
                class = "gdpar_input_error")
  }
  if (!is.null(object)) {
    # Suite instance: a plain list carrying the closures.
    if (is.list(object) && is.function(object$log_prob) &&
        is.function(object$grad_log_prob)) {
      d <- dim %||% object$dim
      return(.gdpar_geom_target_obj(
        log_prob = object$log_prob, grad_log_prob = object$grad_log_prob,
        hessian = hessian, dim = d,
        param_names = param_names %||% object$param_names,
        backend = "closure",
        simulate = simulate %||% object$simulate,
        score = score %||% object$score))
    }
    # cmdstan fit/model compiled with methods.
    if (is.function(object$grad_log_prob) ||
        inherits(object, c("CmdStanFit", "CmdStanMCMC", "CmdStanModel"))) {
      return(.gdpar_geom_target_cmdstan(object, dim, data, param_names,
                                        simulate, score))
    }
    gdpar_abort(
      paste("Unrecognised 'object'. Supply a cmdstanr fit/model compiled with",
            "compile_model_methods = TRUE, or a gdpar_geometry_suite make()",
            "instance, or use the log_prob/grad_log_prob/dim arguments."),
      class = "gdpar_input_error")
  }
  if (!is.function(log_prob) || !is.function(grad_log_prob)) {
    gdpar_abort("Both 'log_prob' and 'grad_log_prob' must be functions.",
                class = "gdpar_input_error")
  }
  if (is.null(dim)) {
    gdpar_abort("Argument 'dim' is required for a closure target.",
                class = "gdpar_input_error")
  }
  assert_count(dim, "dim")
  .gdpar_geom_target_obj(log_prob, grad_log_prob, hessian, as.integer(dim),
                         param_names, "closure", simulate, score)
}

.gdpar_geom_target_obj <- function(log_prob, grad_log_prob, hessian, dim,
                                   param_names, backend, simulate = NULL,
                                   score = NULL) {
  obj <- list(log_prob = log_prob, grad_log_prob = grad_log_prob,
              hessian = hessian, dim = as.integer(dim),
              param_names = param_names %||%
                paste0("theta[", seq_len(dim), "]"),
              backend = backend, simulate = simulate, score = score)
  class(obj) <- c("gdpar_geom_target", "list")
  obj
}

# Wrap a cmdstan fit/model: $log_prob / $grad_log_prob / $hessian act on the
# unconstrained scale and embody decision A (R integrator over Stan backend).
.gdpar_geom_target_cmdstan <- function(object, dim, data, param_names,
                                       simulate = NULL, score = NULL) {
  fit <- object
  if (inherits(object, "CmdStanModel")) {
    require_suggested("cmdstanr", "expose log_prob / grad_log_prob methods")
    fit <- object$sample(data = data, chains = 1, iter_warmup = 1,
                         iter_sampling = 1, refresh = 0, show_messages = FALSE,
                         show_exceptions = FALSE)
  }
  if (!is.function(fit$grad_log_prob)) {
    gdpar_abort(
      paste("The cmdstan object does not expose grad_log_prob; compile the",
            "model with compile_model_methods = TRUE."),
      class = "gdpar_input_error")
  }
  if (is.null(dim)) {
    gdpar_abort("Argument 'dim' is required for a cmdstan target.",
                class = "gdpar_input_error")
  }
  lp <- function(theta) fit$log_prob(unconstrained_variables = theta)
  gl <- function(theta) as.numeric(fit$grad_log_prob(
    unconstrained_variables = theta))
  he <- if (is.function(fit$hessian)) {
    function(theta) fit$hessian(unconstrained_variables = theta)$hessian
  } else NULL
  .gdpar_geom_target_obj(lp, gl, he, as.integer(dim), param_names, "cmdstan",
                         simulate, score)
}

#' Print method for gdpar_geom_target objects
#'
#' @param x A \code{gdpar_geom_target}.
#' @param ... Unused.
#' @return Invisibly returns \code{x}.
#' @export
print.gdpar_geom_target <- function(x, ...) {
  cat("<gdpar_geom_target> backend: ", x$backend, " | dim: ", x$dim, "\n",
      sep = "")
  cat("  hessian available: ", !is.null(x$hessian),
      " | generative: ", !is.null(x$simulate) && !is.null(x$score), "\n",
      sep = "")
  invisible(x)
}

#' Euclidean (constant) metric for the geometric sampling engine
#'
#' Build the level-0/1 metric of the Block RG geometry hierarchy: a
#' position-independent mass matrix. With the identity it is the default
#' diagonal Euclidean metric; with a supplied symmetric positive-definite
#' matrix (or a positive vector of variances) it is the dense Euclidean metric
#' (a constant linear preconditioner), the remedy for a straight anisotropic
#' canyon. The Riemannian level of RG.3 replaces this with a position-dependent
#' metric implementing the same interface.
#'
#' @param dim Integer dimension. Required when \code{M} is \code{NULL}.
#' @param M Optional mass matrix: a \code{dim x dim} symmetric positive-definite
#'   matrix, or a length-\code{dim} positive vector interpreted as its diagonal.
#'   Defaults to the identity.
#'
#' @return A list of class \code{gdpar_geom_metric} with
#'   \code{position_dependent = FALSE} and functions \code{mass(theta)},
#'   \code{inv_mass(theta)}, \code{chol_mass(theta)} (lower Cholesky factor, for
#'   drawing momenta) and \code{logdet(theta)}, each ignoring \code{theta}.
#'
#' @seealso \code{\link{gdpar_geom_hmc}}.
#'
#' @examples
#' m <- gdpar_geom_metric_euclidean(dim = 3)
#' m$mass(c(0, 0, 0))
#'
#' @export
gdpar_geom_metric_euclidean <- function(dim = NULL, M = NULL) {
  if (is.null(M)) {
    if (is.null(dim)) {
      gdpar_abort("Supply either 'dim' or a mass matrix 'M'.",
                  class = "gdpar_input_error")
    }
    assert_count(dim, "dim")
    M <- diag(as.integer(dim))
  } else if (is.vector(M) && !is.matrix(M)) {
    if (any(!is.finite(M)) || any(M <= 0)) {
      gdpar_abort("Diagonal mass vector 'M' must be finite and positive.",
                  class = "gdpar_input_error")
    }
    M <- diag(M, nrow = length(M))
  }
  M <- as.matrix(M)
  if (nrow(M) != ncol(M)) {
    gdpar_abort("Mass matrix 'M' must be square.", class = "gdpar_input_error")
  }
  ch <- tryCatch(chol(M), error = function(e) NULL)
  if (is.null(ch)) {
    gdpar_abort("Mass matrix 'M' must be symmetric positive-definite.",
                class = "gdpar_input_error")
  }
  d <- nrow(M)
  Minv <- chol2inv(ch)
  L <- t(ch)                      # lower factor: M = L L^T, momenta p = L z.
  ld <- 2 * sum(log(diag(ch)))    # log det M.
  obj <- list(
    position_dependent = FALSE, dim = d,
    mass = function(theta) M,
    inv_mass = function(theta) Minv,
    chol_mass = function(theta) L,
    logdet = function(theta) ld)
  class(obj) <- c("gdpar_geom_metric", "list")
  obj
}

#' Riemannian (position-dependent) metric for the geometric sampling engine
#'
#' Build the level-3 metric of the Block RG geometry hierarchy: a
#' position-dependent mass matrix \eqn{M(\theta)} that adapts the sampler's
#' local notion of distance to the curvature of the log-posterior, the remedy
#' for a funnel (variable curvature). Two curvature sources are offered, matching
#' the two ways to obtain the local bending of the density.
#'
#' \describe{
#'   \item{\code{curvature = "fisher"}}{The expected Fisher information (the
#'     natural metric of the statistical manifold; Rao--Amari). It is
#'     positive-definite by construction wherever the model is identifiable, so
#'     no eigenvalue surgery is needed. It is model-specific and must be supplied
#'     as a function \code{fisher(theta)} returning a symmetric positive-definite
#'     matrix (optionally with its derivative \code{dfisher}). This is the
#'     primary, maximally robust choice; its fully general, learned amortisation
#'     across families is a separate Block RG sub-phase.}
#'   \item{\code{curvature = "softabs"}}{The SoftAbs regularisation of the
#'     observed Hessian (Betancourt 2013): the eigenvalues \eqn{\lambda} of the
#'     Hessian of \eqn{-\log\pi} are mapped to \eqn{\lambda\coth(\alpha\lambda)},
#'     turning any bending into a sensible positive mass and flooring nearly flat
#'     directions at \eqn{1/\alpha}. It needs only the Hessian (taken from the
#'     target's \code{$hessian} when available, otherwise finite-differenced from
#'     the gradient), so it is fully general and serves as the cold-start and
#'     extrapolation fallback for the Fisher metric.}
#' }
#'
#' The metric exposes, beyond the Euclidean interface, the spatial derivatives
#' \code{dmass(theta)} (a length-\code{dim} list of \eqn{\partial M/\partial
#' \theta_k}) that the generalised implicit leapfrog of \code{\link{gdpar_geom_hmc}}
#' requires. The Riemannian sampler stays exact regardless of how crude the
#' metric is: the metric is a preconditioner, not part of the target, so the
#' Metropolis correction with the exact log-density is the corrector.
#'
#' @param target A \code{\link{gdpar_geom_target}} (or an object accepted by it).
#'   Used for the Hessian / gradient of the \code{"softabs"} source and for the
#'   dimension.
#' @param curvature Curvature source: \code{"fisher"} (expected Fisher, supplied)
#'   or \code{"softabs"} (SoftAbs of the observed Hessian).
#' @param fisher,dfisher For \code{curvature = "fisher"}: a function
#'   \code{fisher(theta)} returning the \code{dim x dim} expected Fisher matrix,
#'   and optionally \code{dfisher(theta)} returning a length-\code{dim} list of
#'   its partial derivatives. When \code{dfisher} is \code{NULL} the derivatives
#'   are finite-differenced from \code{fisher}.
#' @param alpha SoftAbs softening parameter (\code{curvature = "softabs"}).
#'   Larger values track the true curvature more faithfully (\eqn{\to|\lambda|})
#'   but make the metric vary faster; the default \code{1e6} is the near-absolute
#'   limit.
#' @param floor Minimum eigenvalue imposed on \eqn{M(\theta)} to keep it strictly
#'   positive-definite at genuinely flat directions (non-identification).
#' @param fd_step Finite-difference step for the Hessian (softabs) or for the
#'   metric derivative when an analytic one is unavailable.
#'
#' @return A list of class \code{gdpar_geom_metric} with
#'   \code{position_dependent = TRUE} and functions \code{mass(theta)},
#'   \code{inv_mass(theta)}, \code{chol_mass(theta)}, \code{logdet(theta)} and
#'   \code{dmass(theta)}.
#'
#' @seealso \code{\link{gdpar_geom_hmc}}, \code{\link{gdpar_geom_metric_euclidean}}.
#'
#' @examples
#' # SoftAbs Riemannian metric for a two-dimensional standard normal.
#' tgt <- gdpar_geom_target(
#'   log_prob = function(theta) -0.5 * sum(theta^2),
#'   grad_log_prob = function(theta) -theta, dim = 2)
#' m <- gdpar_geom_metric_riemannian(tgt, curvature = "softabs")
#' m$mass(c(0, 0))
#'
#' @export
gdpar_geom_metric_riemannian <- function(target, curvature = c("fisher", "softabs"),
                                         fisher = NULL, dfisher = NULL,
                                         alpha = 1e6, floor = 1e-8,
                                         fd_step = 1e-4) {
  curvature <- match.arg(curvature)
  if (!inherits(target, "gdpar_geom_target")) target <- gdpar_geom_target(target)
  d <- target$dim
  assert_numeric_scalar(alpha, "alpha", lower = 0)
  assert_numeric_scalar(floor, "floor", lower = 0)
  assert_numeric_scalar(fd_step, "fd_step", lower = 0)
  if (identical(curvature, "fisher")) {
    if (!is.function(fisher)) {
      gdpar_abort(
        "curvature = 'fisher' requires a 'fisher' function returning the expected Fisher matrix.",
        class = "gdpar_input_error")
    }
    mass_fn <- function(theta) .gdpar_geom_floor_spd(as.matrix(fisher(theta)), floor)
    dmass_fn <- if (is.function(dfisher)) {
      function(theta) lapply(dfisher(theta), as.matrix)
    } else {
      function(theta) .gdpar_geom_fd_dmass(mass_fn, theta, fd_step)
    }
  } else {
    hess_lp <- .gdpar_geom_hessian_fn(target, fd_step)
    mass_fn <- function(theta) .gdpar_geom_softabs_mass(-hess_lp(theta), alpha, floor)
    dmass_fn <- function(theta)
      .gdpar_geom_softabs_dmass(theta, hess_lp, alpha, fd_step)
  }
  .gdpar_geom_metric_from_mass(mass_fn, dmass_fn, d, curvature, alpha)
}

# Assemble a gdpar_geom_metric from a mass function and its derivative, caching
# the (mass, Cholesky, inverse) of the last queried position so the implicit
# leapfrog does not redo the factorisation across its fixed-point iterations.
.gdpar_geom_metric_from_mass <- function(mass_fn, dmass_fn, d, kind, alpha) {
  cache <- new.env(parent = emptyenv())
  cache$theta <- NULL
  ensure <- function(theta) {
    if (!is.null(cache$theta) && identical(cache$theta, theta)) return(invisible())
    M <- mass_fn(theta)
    ch <- .gdpar_geom_chol_spd(M)
    cache$theta <- theta; cache$M <- M; cache$ch <- ch
    cache$Minv <- chol2inv(ch); cache$ld <- 2 * sum(log(diag(ch)))
    invisible()
  }
  obj <- list(
    position_dependent = TRUE, dim = d, metric_kind = kind, alpha = alpha,
    mass = function(theta) { ensure(theta); cache$M },
    inv_mass = function(theta) { ensure(theta); cache$Minv },
    chol_mass = function(theta) { ensure(theta); t(cache$ch) },
    logdet = function(theta) { ensure(theta); cache$ld },
    dmass = function(theta) dmass_fn(theta))
  class(obj) <- c("gdpar_geom_metric", "list")
  obj
}

# Cholesky of a symmetric matrix, symmetrising and adding a tiny jitter if the
# raw factorisation fails (a numerically near-singular but valid metric).
.gdpar_geom_chol_spd <- function(M) {
  M <- (M + t(M)) / 2
  ch <- tryCatch(chol(M), error = function(e) NULL)
  if (!is.null(ch)) return(ch)
  jit <- 1e-8 * max(abs(diag(M)), 1)
  chol(M + diag(jit, nrow(M)))
}

# Project a symmetric matrix onto the SPD cone by flooring its eigenvalues.
.gdpar_geom_floor_spd <- function(M, floor) {
  M <- (M + t(M)) / 2
  ev <- eigen(M, symmetric = TRUE)
  lam <- pmax(ev$values, floor)
  ev$vectors %*% (lam * t(ev$vectors))
}

# SoftAbs transform of eigenvalues, lambda * coth(alpha * lambda), with a stable
# series near zero (where the limit is 1 / alpha).
.gdpar_geom_softabs_vals <- function(lambda, alpha) {
  u <- alpha * lambda
  ucoth <- ifelse(abs(u) < 1e-4, 1 + u^2 / 3, u / tanh(u))
  ucoth / alpha
}

# Derivative d/dlambda [lambda * coth(alpha * lambda)] = coth(u) - u csch^2(u),
# u = alpha * lambda, with the near-zero limit (2/3) u.
.gdpar_geom_softabs_deriv <- function(lambda, alpha) {
  u <- alpha * lambda
  small <- abs(u) < 1e-3
  out <- numeric(length(u))
  un <- u[!small]
  out[!small] <- 1 / tanh(un) - un / sinh(un)^2
  out[small] <- (2 / 3) * u[small]
  out
}

# SoftAbs mass M = Q softabs(Lambda) Q^T from the Hessian HU of -log_prob.
.gdpar_geom_softabs_mass <- function(HU, alpha, floor) {
  HU <- (HU + t(HU)) / 2
  ev <- eigen(HU, symmetric = TRUE)
  sav <- pmax(.gdpar_geom_softabs_vals(ev$values, alpha), floor)
  ev$vectors %*% (sav * t(ev$vectors))
}

# Derivative of the SoftAbs mass via the Daleckii--Krein (Loewner) formula:
# dM/dtheta_k = Q ( R .* (Q^T dHU_k Q) ) Q^T, with R the matrix of divided
# differences of softabs over the eigenvalues, and dHU_k the third-derivative
# tensor obtained by central finite differences of the Hessian.
.gdpar_geom_softabs_dmass <- function(theta, hess_lp_fn, alpha, h) {
  HU <- -hess_lp_fn(theta); HU <- (HU + t(HU)) / 2
  ev <- eigen(HU, symmetric = TRUE)
  Q <- ev$vectors; lam <- ev$values; sav <- .gdpar_geom_softabs_vals(lam, alpha)
  d <- length(lam)
  R <- matrix(0, d, d)
  for (i in seq_len(d)) for (j in seq_len(d)) {
    R[i, j] <- if (abs(lam[i] - lam[j]) > 1e-8) {
      (sav[i] - sav[j]) / (lam[i] - lam[j])
    } else {
      .gdpar_geom_softabs_deriv((lam[i] + lam[j]) / 2, alpha)
    }
  }
  lapply(seq_len(d), function(k) {
    e <- numeric(d); e[k] <- h
    HUp <- -hess_lp_fn(theta + e); HUm <- -hess_lp_fn(theta - e)
    dHU <- ((HUp + t(HUp)) - (HUm + t(HUm))) / (4 * h)
    A <- crossprod(Q, dHU %*% Q)
    Q %*% (R * A) %*% t(Q)
  })
}

# A Hessian-of-log_prob function: the target's own when exposed (cmdstan
# $hessian or an analytic closure), else central finite differences of the
# gradient (symmetrised).
.gdpar_geom_hessian_fn <- function(target, h) {
  if (is.function(target$hessian)) return(target$hessian)
  function(theta) {
    d <- length(theta)
    H <- matrix(0, d, d)
    for (k in seq_len(d)) {
      e <- numeric(d); e[k] <- h
      H[, k] <- (target$grad_log_prob(theta + e) -
                   target$grad_log_prob(theta - e)) / (2 * h)
    }
    (H + t(H)) / 2
  }
}

# Central finite-difference derivative of a matrix-valued mass function.
.gdpar_geom_fd_dmass <- function(mass_fn, theta, h) {
  d <- length(theta)
  lapply(seq_len(d), function(k) {
    e <- numeric(d); e[k] <- h
    (mass_fn(theta + e) - mass_fn(theta - e)) / (2 * h)
  })
}

# ---------------------------------------------------------------------------
# RG.3 Capa 2: learned Gaussian-process surrogate of the expected Fisher.
#
# The Riemannian metric M(theta) = L(theta) L(theta)^T is built from a smooth,
# deterministic, differentiable surrogate of the expected Fisher information --
# the general realisation of the natural Rao--Amari metric where no closed form
# exists. Following the organic-critical reading of ORPHEUS-PIMC-A (sections 8
# and 16), the surrogate is a Gaussian process whose MEAN FUNCTION is the
# SoftAbs curvature of cold-start Capa 1 and which learns only the RESIDUAL to
# the (closed-form or simulated) expected Fisher, in the log-Cholesky
# parametrisation. This single object unifies three ORPHEUS components without
# inheriting their drawbacks:
#   * surrogate + graceful fallback: far from the reservoir the kernel decays to
#     zero, the posterior mean returns to the SoftAbs mean, and the metric
#     degrades CONTINUOUSLY to the always-available SoftAbs -- there is no hard
#     metric switch that would break the implicit leapfrog's reversibility;
#   * principled novelty (section 12.2.5): the predictive variance is the
#     epistemic-uncertainty extrapolation detector;
#   * smoothness as a structural necessity (D90): a raw Monte-Carlo Fisher is
#     non-smooth and would break reversibility, whereas the GP posterior mean is
#     smooth and analytically differentiable.
# The sampler stays EXACT for any surrogate quality: the metric is a
# preconditioner, not part of the target, so the Metropolis correction with the
# exact density is the corrector -- no delayed acceptance is needed (the
# improvement over ORPHEUS, where the surrogate enters the acceptance).
#
# Log-Cholesky parametrisation. A symmetric positive-definite M is encoded by
# the lower Cholesky factor L (M = L L^T) through psi: the LOG of the diagonal
# entries and the raw strict-lower entries, in column-major lower-triangular
# order. The map psi -> M is smooth and SPD for any psi, so a smooth surrogate
# on psi yields a smooth SPD metric BY CONSTRUCTION.

# Column-major lower-triangular index (including the diagonal); the diagonal
# entries (j, j) are the first of each column block, in order j = 1..d.
.gdpar_geom_tri_index <- function(d) {
  which(lower.tri(matrix(0, d, d), diag = TRUE), arr.ind = TRUE)
}

# psi (and the lower factor L) from an SPD matrix M.
.gdpar_geom_logchol <- function(M, idx) {
  L <- t(chol((M + t(M)) / 2))
  diagp <- idx[, 1] == idx[, 2]
  psi <- numeric(nrow(idx))
  psi[diagp] <- log(diag(L))
  psi[!diagp] <- L[idx[!diagp, , drop = FALSE]]
  list(L = L, psi = psi)
}

# Lower factor L from psi (diagonal through exp, so L_ii > 0 always -> SPD).
.gdpar_geom_L_from_psi <- function(psi, idx, d) {
  L <- matrix(0, d, d)
  diagp <- idx[, 1] == idx[, 2]
  L[idx[diagp, , drop = FALSE]] <- exp(psi[diagp])
  L[idx[!diagp, , drop = FALSE]] <- psi[!diagp]
  L
}

# dM from dpsi: dL_ii = L_ii dpsi_ii (chain through exp), dL_ij = dpsi_ij;
# dM = dL L^T + L dL^T = (dL L^T) + (dL L^T)^T.
.gdpar_geom_dM_from_dpsi <- function(L, dpsi, idx, d) {
  dL <- matrix(0, d, d)
  diagp <- idx[, 1] == idx[, 2]
  id <- idx[diagp, , drop = FALSE]
  dL[id] <- L[id] * dpsi[diagp]
  dL[idx[!diagp, , drop = FALSE]] <- dpsi[!diagp]
  A <- dL %*% t(L)
  A + t(A)
}

# dpsi from dM (the Cholesky differential): P = L^{-1} dM L^{-T} (symmetric);
# Phi(P) = strict-lower + half-diagonal; dL = L Phi(P); then dpsi_ii = dL_ii /
# L_ii and dpsi_ij = dL_ij. Inverse map of .gdpar_geom_dM_from_dpsi.
.gdpar_geom_dpsi_from_dM <- function(L, dM, idx, d) {
  X <- forwardsolve(L, (dM + t(dM)) / 2)
  P <- t(forwardsolve(L, t(X)))
  Phi <- P
  Phi[upper.tri(Phi)] <- 0
  diag(Phi) <- diag(P) / 2
  dL <- L %*% Phi
  diagp <- idx[, 1] == idx[, 2]
  dpsi <- numeric(nrow(idx))
  dpsi[diagp] <- diag(dL) / diag(L)
  dpsi[!diagp] <- dL[idx[!diagp, , drop = FALSE]]
  dpsi
}

#' Learned Gaussian-process Riemannian metric (expected Fisher surrogate)
#'
#' Build the general, learned realisation of the level-3 Riemannian metric of
#' the Block RG hierarchy: a position-dependent mass matrix
#' \eqn{M(\theta) = L(\theta) L(\theta)^\top} whose log-Cholesky factor is a
#' Gaussian-process surrogate of the expected Fisher information (the natural
#' Rao--Amari metric), for use where the Fisher has no closed form. Where the
#' Fisher \emph{is} closed (the funnel, a Gaussian, generalised-linear-model
#' slots) \code{\link{gdpar_geom_metric_riemannian}} is exact and preferable;
#' this surrogate covers the general case and is validated against those closed
#' forms.
#'
#' @details
#' The surrogate's \emph{mean function} is the SoftAbs curvature of the observed
#' Hessian (the cold-start metric, always positive-definite); the Gaussian
#' process learns only the smooth residual to the expected Fisher at the
#' reservoir sites, in the log-Cholesky parametrisation. Three properties follow
#' by construction:
#' \itemize{
#'   \item \strong{Positive-definite always}: \eqn{M = L L^\top} with
#'     \eqn{L_{ii} = \exp(\cdot) > 0}.
#'   \item \strong{Graceful degradation}: far from the reservoir the kernel
#'     decays to zero and the posterior mean returns to the SoftAbs mean, so the
#'     metric degrades \emph{continuously} to the always-available SoftAbs --
#'     there is no hard metric switch that would break the reversibility of the
#'     implicit generalised leapfrog of \code{\link{gdpar_geom_hmc}}.
#'   \item \strong{Exactness independent of surrogate quality}: the metric is a
#'     preconditioner, not part of the target, so the Metropolis correction with
#'     the exact density keeps the sampler exact for any surrogate; delayed
#'     acceptance is unnecessary.
#' }
#' The spatial derivative \code{dmass(theta)} is closed form: the analytic
#' kernel derivative for the learned residual plus the Daleckii--Krein
#' derivative of the SoftAbs mean (reused from Capa 1), pushed through the
#' log-Cholesky map. The predictive standard deviation is exposed as
#' \code{novelty(theta)}, the epistemic-uncertainty extrapolation detector.
#'
#' @param target A \code{\link{gdpar_geom_target}} (or an object accepted by it),
#'   used for the SoftAbs mean (its Hessian / gradient) and for the dimension.
#' @param fisher A function \code{fisher(theta)} returning the expected Fisher
#'   information at \code{theta} as a \code{dim x dim} symmetric positive-definite
#'   matrix. The training targets at the reservoir sites are evaluated through
#'   it; a simulation-based estimator (the score outer product over simulated
#'   data sets) plugs into this same slot.
#' @param sites A numeric \code{m x dim} matrix of reservoir positions (one row
#'   per site) at which the Fisher is evaluated to train the surrogate. See
#'   \code{\link{gdpar_geom_reservoir}} for collecting sites from a warmup run.
#' @param weights Optional positive vector of length \code{m} of importance
#'   weights (e.g. \eqn{1/Q}) reweighting a reservoir biased towards rare
#'   positions back to the typical set; entered as per-site noise scaling.
#'   Defaults to equal weights.
#' @param lengthscale Optional positive radial-basis-function length-scale on the
#'   standardised inputs. Defaults to the median pairwise-distance heuristic.
#' @param nugget Non-negative kernel nugget (observation-noise variance). Small
#'   values interpolate the Fisher at the sites; larger values smooth it.
#' @param alpha,floor,fd_step SoftAbs softening, eigenvalue floor and
#'   finite-difference step governing the SoftAbs mean function (passed to the
#'   Capa 1 machinery).
#'
#' @return A list of class \code{gdpar_geom_metric} with
#'   \code{position_dependent = TRUE}, \code{metric_kind = "gp_fisher"}, the
#'   functions \code{mass}, \code{inv_mass}, \code{chol_mass}, \code{logdet} and
#'   \code{dmass} of the metric interface, plus \code{novelty(theta)} (the
#'   predictive standard deviation, higher = more out-of-distribution) and the
#'   fields \code{n_sites} and \code{lengthscale}.
#'
#' @seealso \code{\link{gdpar_geom_metric_riemannian}},
#'   \code{\link{gdpar_geom_reservoir}}, \code{\link{gdpar_geom_hmc}}.
#'
#' @examples
#' # A two-dimensional target whose expected Fisher is the identity; the learned
#' # metric recovers it from a small reservoir.
#' tgt <- gdpar_geom_target(
#'   log_prob = function(theta) -0.5 * sum(theta^2),
#'   grad_log_prob = function(theta) -theta, dim = 2)
#' sites <- matrix(stats::rnorm(40), ncol = 2)
#' m <- gdpar_geom_metric_gp_fisher(tgt, fisher = function(theta) diag(2),
#'                                  sites = sites)
#' round(m$mass(c(0, 0)), 3)
#'
#' @export
gdpar_geom_metric_gp_fisher <- function(target, fisher, sites, weights = NULL,
                                        lengthscale = NULL, nugget = 1e-6,
                                        alpha = 1e6, floor = 1e-8,
                                        fd_step = 1e-4) {
  if (!inherits(target, "gdpar_geom_target")) target <- gdpar_geom_target(target)
  d <- target$dim
  if (!is.function(fisher)) {
    gdpar_abort(
      "'fisher' must be a function returning the expected Fisher matrix at a position.",
      class = "gdpar_input_error")
  }
  if (is.numeric(sites) && is.null(dim(sites))) sites <- matrix(sites, ncol = d)
  if (!is.matrix(sites) || !is.numeric(sites) || ncol(sites) != d) {
    gdpar_abort(
      sprintf("'sites' must be a numeric matrix with %d columns (one row per reservoir site).", d),
      class = "gdpar_input_error")
  }
  m <- nrow(sites)
  if (m < 2L) {
    gdpar_abort("'sites' needs at least two reservoir sites.",
                class = "gdpar_input_error")
  }
  assert_numeric_scalar(nugget, "nugget", lower = 0)
  assert_numeric_scalar(alpha, "alpha", lower = 0)
  assert_numeric_scalar(floor, "floor", lower = 0)
  assert_numeric_scalar(fd_step, "fd_step", lower = 0)
  idx <- .gdpar_geom_tri_index(d)
  q <- nrow(idx)

  # SoftAbs mean machinery (the cold-start curvature, Capa 1).
  hess_lp <- .gdpar_geom_hessian_fn(target, fd_step)
  softabs_mass <- function(theta)
    .gdpar_geom_softabs_mass(-hess_lp(theta), alpha, floor)

  # Standardise the inputs; RBF length-scale by the median-distance heuristic.
  center <- colMeans(sites)
  scale <- apply(sites, 2, stats::sd)
  scale[!is.finite(scale) | scale <= 0] <- 1
  Zr <- sweep(sweep(sites, 2, center, "-"), 2, scale, "/")
  if (is.null(lengthscale)) {
    dd <- stats::dist(Zr)
    lengthscale <- stats::median(dd[dd > 0])
    if (!is.finite(lengthscale) || lengthscale <= 0) lengthscale <- sqrt(d)
  }
  assert_numeric_scalar(lengthscale, "lengthscale", lower = 0)
  ell2 <- lengthscale^2

  # Reservoir kernel, importance noise, and residual targets (Fisher minus the
  # SoftAbs mean, both in log-Cholesky coordinates).
  sq <- rowSums(Zr^2)
  K <- exp(-0.5 * pmax(outer(sq, sq, "+") - 2 * tcrossprod(Zr), 0) / ell2)
  if (is.null(weights)) {
    noise <- rep(nugget, m)
  } else {
    if (!is.numeric(weights) || length(weights) != m ||
        any(!is.finite(weights)) || any(weights <= 0)) {
      gdpar_abort(
        "'weights' must be a positive finite vector with one entry per site.",
        class = "gdpar_input_error")
    }
    noise <- nugget / (weights / mean(weights))   # 1/Q importance weighting.
  }
  Kn <- K + diag(noise, m)
  Kn_inv <- chol2inv(.gdpar_geom_chol_spd(Kn))
  Yres <- matrix(0, m, q)
  for (r in seq_len(m)) {
    sr <- sites[r, ]
    psiF <- .gdpar_geom_logchol(as.matrix(fisher(sr)), idx)$psi
    psi0 <- .gdpar_geom_logchol(softabs_mass(sr), idx)$psi
    Yres[r, ] <- psiF - psi0
  }
  coefs <- Kn_inv %*% Yres                         # m x q dual coefficients.

  # Per-position state shared by mass / inv_mass / chol_mass / logdet / dmass /
  # novelty so a single query computes the SoftAbs mean and the GP factor once.
  st <- new.env(parent = emptyenv()); st$theta <- NULL
  ensure <- function(theta) {
    if (!is.null(st$theta) && identical(st$theta, theta)) return(invisible())
    z <- (theta - center) / scale
    diffs <- sweep(Zr, 2, z, "-")
    kstar <- exp(-0.5 * rowSums(diffs^2) / ell2)
    lc0 <- .gdpar_geom_logchol(softabs_mass(theta), idx)
    psi <- lc0$psi + as.numeric(crossprod(kstar, coefs))
    st$theta <- theta; st$z <- z; st$kstar <- kstar
    st$L0 <- lc0$L; st$L <- .gdpar_geom_L_from_psi(psi, idx, d); st$dM0 <- NULL
    invisible()
  }
  obj <- list(
    position_dependent = TRUE, dim = d, metric_kind = "gp_fisher",
    alpha = alpha, n_sites = m, lengthscale = lengthscale,
    mass = function(theta) { ensure(theta); tcrossprod(st$L) },
    inv_mass = function(theta) { ensure(theta); chol2inv(t(st$L)) },
    chol_mass = function(theta) { ensure(theta); st$L },
    logdet = function(theta) { ensure(theta); 2 * sum(log(diag(st$L))) },
    dmass = function(theta) {
      ensure(theta)
      if (is.null(st$dM0)) {
        st$dM0 <- .gdpar_geom_softabs_dmass(theta, hess_lp, alpha, fd_step)
      }
      lapply(seq_len(d), function(j) {
        dpsi0 <- .gdpar_geom_dpsi_from_dM(st$L0, st$dM0[[j]], idx, d)
        dk <- st$kstar * (-(st$z[j] - Zr[, j]) / ell2) / scale[j]
        .gdpar_geom_dM_from_dpsi(st$L, dpsi0 + as.numeric(crossprod(dk, coefs)),
                                 idx, d)
      })
    },
    novelty = function(theta) {
      ensure(theta)
      sqrt(max(0, 1 - as.numeric(crossprod(st$kstar, Kn_inv %*% st$kstar))))
    })
  class(obj) <- c("gdpar_geom_metric", "list")
  obj
}

#' Collect a reservoir of positions for the learned Gaussian-process metric
#'
#' Phase one of the two-phase, decoupled-archivist design (ORPHEUS-PIMC-A
#' section 16): run a short warmup of the geometric sampler and return the
#' retained positions as a reservoir of sites at which the expected Fisher is
#' later evaluated to train \code{\link{gdpar_geom_metric_gp_fisher}}. For an
#' easy or moderately curved target the default Euclidean warmup suffices; for a
#' funnel-like target pass a SoftAbs metric (or a hand-designed sweep of the
#' curvature axis) instead.
#'
#' @param target A \code{\link{gdpar_geom_target}} (or an object accepted by it).
#' @param n_sites Integer number of reservoir sites to return (the retained
#'   warmup draws).
#' @param metric Optional metric for the warmup run (defaults to the Euclidean
#'   identity); a \code{\link{gdpar_geom_metric_riemannian}} SoftAbs metric is a
#'   robust choice for curved targets.
#' @param epsilon,L,n_warmup Leapfrog step, trajectory length and discarded
#'   warmup length for the collecting run.
#' @param init Optional initial position (defaults to zeros).
#' @param seed Optional integer seed for reproducibility.
#'
#' @return A numeric \code{n_sites x dim} matrix of reservoir positions.
#'
#' @seealso \code{\link{gdpar_geom_metric_gp_fisher}}, \code{\link{gdpar_geom_hmc}}.
#'
#' @examples
#' tgt <- gdpar_geom_target(
#'   log_prob = function(theta) -0.5 * sum(theta^2),
#'   grad_log_prob = function(theta) -theta, dim = 2)
#' sites <- gdpar_geom_reservoir(tgt, n_sites = 30, n_warmup = 100, seed = 1)
#' dim(sites)
#'
#' @export
gdpar_geom_reservoir <- function(target, n_sites = 50L, metric = NULL,
                                 epsilon = 0.25, L = 15L, n_warmup = 200L,
                                 init = NULL, seed = NULL) {
  if (!inherits(target, "gdpar_geom_target")) target <- gdpar_geom_target(target)
  assert_count(n_sites, "n_sites")
  fit <- gdpar_geom_hmc(target, metric = metric, epsilon = epsilon, L = L,
                        n_iter = n_sites, n_warmup = n_warmup, init = init,
                        seed = seed)
  fit$draws
}

# A small deterministic non-negative integer derived from a position, so a
# simulation-based estimate evaluated at theta is a reproducible function of
# theta regardless of call order (each evaluation reseeds from this key).
.gdpar_geom_seed_from_theta <- function(theta) {
  v <- as.numeric(theta)
  v[!is.finite(v)] <- 0
  key <- sum(floor(abs(v) * 1009 + 1) * (seq_along(v) * 7919))
  as.integer(key %% 2147483646) + 1L
}

#' Simulation-based estimator of the expected Fisher information
#'
#' Build the \code{fisher(theta)} function that feeds the learned metric
#' \code{\link{gdpar_geom_metric_gp_fisher}} where the expected Fisher
#' information has no closed form -- the general realisation completing the
#' Riemannian level of the Block RG hierarchy. The expected Fisher is
#' \eqn{I(\theta) = \mathbb{E}_{y\sim p(\cdot\mid\theta)}\!\left[s(\theta, y)\,
#' s(\theta, y)^\top\right]} with \eqn{s} the score (the gradient of the
#' log-likelihood). The estimator is the average outer product of the scores
#' over \code{n_sim} data sets simulated from the model at \eqn{\theta}: it is
#' \strong{positive semi-definite by construction} (a sum of rank-one outer
#' products) and unbiased for the expected Fisher.
#'
#' @details
#' The estimate is a \emph{deterministic} function of \eqn{\theta}: each
#' evaluation reseeds the RNG from a key derived from \eqn{\theta} and the base
#' \code{seed} (and the RNG state is saved and restored around the call), so the
#' surrogate trained on it is reproducible and independent of call order. The
#' simulation noise is then absorbed by the Gaussian-process surrogate: its
#' SoftAbs mean carries the bulk of the curvature deterministically and the
#' process learns only the smooth residual, so the SoftAbs acts as a structural
#' control variate and few replicates per site suffice (the deferred half of
#' decision D91, now closed). Antithetic simulation is deliberately \emph{not}
#' used: for location families the score is odd in the centred data, so the
#' antithetic partner has score \eqn{-s} and \eqn{(-s)(-s)^\top = s s^\top}
#' leaves the outer product unchanged -- no variance reduction.
#'
#' For a well-conditioned estimate take \code{n_sim} comfortably larger than the
#' dimension (the average of fewer than \code{dim} rank-one terms is singular and
#' is only made strictly positive-definite by the eigenvalue \code{floor}).
#'
#' @param target A generative \code{\link{gdpar_geom_target}} carrying both
#'   \code{simulate(theta)} (one synthetic data set drawn from the model at
#'   \code{theta}) and \code{score(theta, y)} (the gradient of the
#'   log-likelihood of \code{y} at \code{theta}, on the unconstrained scale).
#' @param n_sim Integer number of simulated data sets averaged per evaluation.
#' @param seed Integer base seed combined with the position key.
#' @param floor Minimum eigenvalue imposed on the returned matrix so it is
#'   strictly positive-definite (needed by the log-Cholesky map of the surrogate).
#'
#' @return A function \code{fisher(theta)} returning the \code{dim x dim}
#'   estimated expected Fisher matrix (symmetric positive-definite), with the
#'   number of simulations recorded in its \code{"n_sim"} attribute. Pass it as
#'   the \code{fisher} argument of \code{\link{gdpar_geom_metric_gp_fisher}}.
#'
#' @seealso \code{\link{gdpar_geom_metric_gp_fisher}},
#'   \code{\link{gdpar_geom_target}}, \code{\link{gdpar_geom_rmhmc_adaptive}}.
#'
#' @examples
#' # A bivariate normal location model y ~ N(theta, Sigma0): the expected Fisher
#' # is the constant precision Sigma0^{-1}, which the estimator recovers.
#' Sigma0 <- matrix(c(1, 0.3, 0.3, 2), 2, 2)
#' P0 <- solve(Sigma0)
#' tgt <- gdpar_geom_target(
#'   log_prob = function(theta) -0.5 * sum(theta^2),
#'   grad_log_prob = function(theta) -theta, dim = 2,
#'   simulate = function(theta)
#'     as.numeric(theta + t(chol(Sigma0)) %*% stats::rnorm(2)),
#'   score = function(theta, y) as.numeric(P0 %*% (y - theta)))
#' fisher <- gdpar_geom_fisher_simulator(tgt, n_sim = 4000, seed = 1)
#' round(fisher(c(0, 0)), 2)   # close to solve(Sigma0)
#'
#' @export
gdpar_geom_fisher_simulator <- function(target, n_sim = 64L, seed = 1L,
                                        floor = 1e-8) {
  if (!inherits(target, "gdpar_geom_target")) target <- gdpar_geom_target(target)
  if (!is.function(target$simulate) || !is.function(target$score)) {
    gdpar_abort(
      paste("gdpar_geom_fisher_simulator() needs a generative target: supply",
            "'simulate' and 'score' to gdpar_geom_target()."),
      class = "gdpar_input_error")
  }
  assert_count(n_sim, "n_sim")
  assert_count(seed, "seed")
  assert_numeric_scalar(floor, "floor", lower = 0)
  d <- target$dim
  base_seed <- as.integer(seed)
  function(theta) {
    sub <- (base_seed + .gdpar_geom_seed_from_theta(theta)) %% 2147483646 + 1L
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      old <- get(".Random.seed", envir = .GlobalEnv)
      on.exit(assign(".Random.seed", old, envir = .GlobalEnv), add = TRUE)
    } else {
      on.exit(rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
    }
    set.seed(sub)
    acc <- matrix(0, d, d)
    for (i in seq_len(n_sim)) {
      s <- as.numeric(target$score(theta, target$simulate(theta)))
      if (length(s) != d || any(!is.finite(s))) {
        gdpar_abort(
          "The generative 'score' must return a finite gradient of length dim.",
          class = "gdpar_input_error")
      }
      acc <- acc + tcrossprod(s)
    }
    I <- .gdpar_geom_floor_spd(acc / n_sim, floor)
    attr(I, "n_sim") <- as.integer(n_sim)
    I
  }
}

# Exact harmonic (Gaussian) flow of the reference quadratic, the wall part of the
# sub-Riemannian splitting. In the Fisher eigenbasis U each mode i evolves under
# H_i = 0.5 q_i^2 + 0.5 omega_i^2 z_i^2 (z = U^T (theta - ref), q = U^T p). The
# flow is the closed-form symplectic rotation [[cos, sin/omega], [-omega sin,
# cos]]; as omega -> 0 it degrades continuously to the free drift [[1, t],[0, 1]]
# (the canyon floor). Exact and reversible for any step t, with no step-size
# penalty from the stiff (high-omega) wall modes.
.gdpar_geom_subriemann_flow <- function(theta, p, ref, U, omega, t) {
  z <- as.numeric(crossprod(U, theta - ref))
  q <- as.numeric(crossprod(U, p))
  wt <- omega * t
  c_wt <- cos(wt)
  drift <- omega <= 1e-8
  s_over_w <- ifelse(drift, t, sin(wt) / omega)   # sin(omega t)/omega, limit t.
  w_sin <- ifelse(drift, 0, omega * sin(wt))       # omega sin(omega t), limit 0.
  z_new <- c_wt * z + s_over_w * q
  q_new <- -w_sin * z + c_wt * q
  list(theta = as.numeric(ref + U %*% z_new),
       p = as.numeric(U %*% q_new))
}

#' Sub-Riemannian metric and integrator for the geometric sampling engine
#'
#' Build the level-5 geometry of the Block RG hierarchy: the remedy for a
#' quasi-deterministic posterior (the eBird count / tweedie case), where the
#' typical set contracts onto a lower-dimensional manifold and the expected
#' Fisher information grows without bound along the stiff "wall" directions while
#' the soft "floor" directions still carry genuine variation. A sub-Riemannian
#' structure equips only a distribution \eqn{D_\theta \subsetneq T_\theta M} of
#' accessible directions with an inner product (Montgomery 2002): the sampler
#' glides along the floor instead of fighting the walls with a vanishing step
#' size.
#'
#' @details
#' The accessible distribution is read from the \strong{near-null space of the
#' expected Fisher}. At a reference position the Fisher is eigendecomposed,
#' \eqn{I = U \Lambda U^\top}, and a \strong{continuous} verticality filter
#' \eqn{w_i = \sigma((\log\lambda_i - \log\tau)/s)} assigns each eigendirection a
#' weight in \eqn{(0, 1)} -- one for a wall (large \eqn{\lambda}), zero for the
#' floor (small \eqn{\lambda}), with no hard cut, so a borderline direction is
#' blended smoothly. The wall curvature \eqn{A = U\,\mathrm{diag}(w_i\lambda_i)\,
#' U^\top} (frequencies \eqn{\omega_i = \sqrt{w_i\lambda_i}}) defines a fixed
#' reference quadratic.
#'
#' Sampling uses a Strang splitting with a Euclidean (identity) kinetic energy:
#' each step is a half momentum kick by \eqn{\nabla U_{\mathrm{rest}} =
#' -\nabla\log p(\theta) - A(\theta - \theta^*)}, an \strong{exact} harmonic flow
#' of the reference quadratic (\code{.gdpar_geom_subriemann_flow}; a closed-form
#' symplectic rotation per mode), and a second half kick. The stiff walls are
#' integrated exactly regardless of the step size, so the step is limited only by
#' the gentle floor; the soft directions follow the free-drift limit. The scheme
#' is symplectic and time-reversible by construction (split HMC with a Gaussian
#' part; Shahbaba et al. 2014). Because the reference quadratic is a
#' preconditioner inside the integrator and \strong{not} part of the target, the
#' Metropolis correction with the exact log-density of \code{\link{gdpar_geom_hmc}}
#' keeps the sampler exact however coarse the Gaussian approximation of the walls
#' is: only efficiency, never correctness, depends on it.
#'
#' @param target A \code{\link{gdpar_geom_target}} (or an object accepted by it).
#' @param fisher A function \code{fisher(theta)} returning the \code{dim x dim}
#'   expected Fisher information matrix (symmetric positive semi-definite); for
#'   models without a closed form pass
#'   \code{\link{gdpar_geom_fisher_simulator}}. Evaluated once at
#'   \code{reference}.
#' @param reference Optional numeric vector at which the Fisher and the reference
#'   quadratic are evaluated; for a real model use a warmup mode or mean.
#'   Defaults to zeros.
#' @param tau Eigenvalue threshold of the verticality filter (directions with
#'   \eqn{\lambda \gg \tau} are walls). Defaults to the smallest Fisher
#'   eigenvalue (the floor scale): with the rational filter the leapfrog residual
#'   curvature is then capped at the floor, so the step is never limited by the
#'   stiff walls.
#' @param softness Positive logistic width of the filter in log-eigenvalue. Small
#'   values approach a hard floor/wall split; the default \code{1} keeps the
#'   blend smooth.
#' @param floor Minimum eigenvalue imposed on the Fisher before the filter, so a
#'   genuinely flat direction stays well defined.
#'
#' @return A list of class \code{gdpar_geom_metric} with
#'   \code{position_dependent = FALSE}, the identity kinetic interface
#'   (\code{mass}, \code{inv_mass}, \code{chol_mass}, \code{logdet}), an
#'   \code{integrator} closure consumed by \code{\link{gdpar_geom_hmc}}, and the
#'   diagnostic fields \code{metric_kind = "sub_riemannian"}, \code{reference},
#'   \code{eigenvalues}, \code{verticality} (the \eqn{w_i}), \code{frequencies}
#'   (the \eqn{\omega_i}), \code{n_walls}, \code{tau}, \code{softness} and
#'   \code{suggested_epsilon} (a leapfrog step matched to the floor scale, since
#'   the exact walls impose no step-size limit -- the source of the speed-up
#'   over a Euclidean sampler, whose step is bottlenecked by the stiffest wall).
#'
#' @references
#' Montgomery, R. (2002) \emph{A Tour of Subriemannian Geometries, Their
#' Geodesics and Applications}. AMS.
#'
#' Shahbaba, B., Lan, S., Johnson, W. O. and Neal, R. M. (2014) Split
#' Hamiltonian Monte Carlo. \emph{Statistics and Computing} 24, 339--349.
#'
#' @seealso \code{\link{gdpar_geom_hmc}}, \code{\link{gdpar_geom_fisher_simulator}},
#'   \code{\link{gdpar_geom_metric_riemannian}}.
#'
#' @examples
#' # A two-dimensional quasi-deterministic canyon: a soft floor (variance one)
#' # and a stiff wall (variance one hundredth). The expected Fisher is constant.
#' tgt <- gdpar_geom_target(
#'   log_prob = function(theta) -0.5 * (theta[1]^2 + 100 * theta[2]^2),
#'   grad_log_prob = function(theta) -c(theta[1], 100 * theta[2]), dim = 2)
#' metric <- gdpar_geom_metric_subriemannian(
#'   tgt, fisher = function(theta) diag(c(1, 100)))
#' fit <- gdpar_geom_hmc(tgt, metric = metric, epsilon = 0.5, L = 10,
#'                       n_iter = 200, n_warmup = 100, seed = 1)
#' fit$ebfmi
#'
#' @export
gdpar_geom_metric_subriemannian <- function(target, fisher, reference = NULL,
                                            tau = NULL, softness = 1,
                                            floor = 1e-8) {
  if (!inherits(target, "gdpar_geom_target")) target <- gdpar_geom_target(target)
  d <- target$dim
  if (!is.function(fisher)) {
    gdpar_abort(
      paste("'fisher' must be a function returning the expected Fisher matrix",
            "(e.g. gdpar_geom_fisher_simulator())."),
      class = "gdpar_input_error")
  }
  assert_numeric_scalar(softness, "softness", lower = 0)
  assert_numeric_scalar(floor, "floor", lower = 0)
  ref <- if (is.null(reference)) rep(0, d) else as.numeric(reference)
  if (length(ref) != d) {
    gdpar_abort(sprintf("'reference' must have length %d.", d),
                class = "gdpar_input_error")
  }
  I_ref <- .gdpar_geom_floor_spd(as.matrix(fisher(ref)), floor)
  if (nrow(I_ref) != d || ncol(I_ref) != d) {
    gdpar_abort(sprintf("'fisher' must return a %d x %d matrix.", d, d),
                class = "gdpar_input_error")
  }
  ev <- eigen(I_ref, symmetric = TRUE)
  lambda <- pmax(ev$values, floor)
  U <- ev$vectors
  # Default threshold at the floor scale (the softest direction's curvature):
  # with the rational filter the leapfrog residual curvature is tau*lambda /
  # (lambda + tau) <= tau for every direction, so no residual is ever stiffer
  # than the floor and the step size is governed by the floor alone, never the
  # walls. (A threshold at the spectral mean would leave stiff residuals and
  # forfeit the stiffness removal.)
  if (is.null(tau)) tau <- min(lambda)
  assert_numeric_scalar(tau, "tau", lower = 0)
  w <- 1 / (1 + exp(-(log(lambda) - log(tau)) / softness))   # verticality.
  a <- w * lambda                                            # wall curvature.
  omega <- sqrt(pmax(a, 0))

  # Residual potential gradient: the full -grad log p minus the reference
  # quadratic A (theta - ref) that the exact flow already integrates.
  A <- U %*% (a * t(U))
  grad_U_rest <- function(theta)
    -target$grad_log_prob(theta) - as.numeric(A %*% (theta - ref))

  # Strang trajectory: half kick, exact harmonic flow, half kick, L times. Any
  # non-finite iterate is reported as non-convergence so the proposal is
  # rejected rather than crashing the run.
  integrator <- function(theta, p, target, eps, L,
                         fp_tol = NULL, fp_max = NULL) {
    ok <- TRUE
    for (i in seq_len(L)) {
      p <- p - (eps / 2) * grad_U_rest(theta)
      st <- .gdpar_geom_subriemann_flow(theta, p, ref, U, omega, eps)
      theta <- st$theta; p <- st$p
      p <- p - (eps / 2) * grad_U_rest(theta)
      if (any(!is.finite(theta)) || any(!is.finite(p))) { ok <- FALSE; break }
    }
    list(theta = theta, p = p, converged = ok)
  }

  # The leapfrog integrates only the residual curvature (1 - w) * lambda; the
  # walls are exact and impose no step-size limit. A suggested step that keeps
  # the residual leapfrog stable (well inside the omega*eps < 2 limit) is the
  # floor scale, not the wall scale -- the whole point of the geometry. This is
  # a hint exposed for the caller, not enforced.
  resid_max <- max((1 - w) * lambda, 0)
  suggested_epsilon <- if (resid_max <= 0) 2 else min(2, 0.7 / sqrt(resid_max))
  Id <- diag(d)
  obj <- list(
    position_dependent = FALSE, dim = d, metric_kind = "sub_riemannian",
    mass = function(theta) Id, inv_mass = function(theta) Id,
    chol_mass = function(theta) Id, logdet = function(theta) 0,
    integrator = integrator,
    reference = ref, eigenvalues = lambda, verticality = w,
    frequencies = omega, n_walls = sum(w > 0.5), tau = tau, softness = softness,
    suggested_epsilon = suggested_epsilon)
  class(obj) <- c("gdpar_geom_metric", "list")
  obj
}

# ---------------------------------------------------------------------------
# RG.4: Finsler / relativistic geometry -- a bounded, non-Gaussian kinetic
# energy coupled to the position-dependent Riemannian metric of RG.3. The remedy
# for heavy tails and directional anisotropy (the G3 target), NOT the count
# (that is the sub-Riemannian level of RG.4 above).
#
# Why a bounded kinetic energy. A Gaussian kinetic energy gives the unbounded
# velocity grad_p K = M^{-1} p, so a large momentum drawn in a heavy tail
# produces an arbitrarily large position step that a fixed-step integrator
# cannot follow (overshoot, divergences, a step size that must be tuned to the
# stiffest region). The relativistic kinetic energy caps the velocity at a
# finite speed c (Lu et al. 2017; Livingstone, Faulkner & Roberts 2019): however
# large the momentum, the step stays bounded, taming the tails and the
# ill-conditioning.
#
# Coupled to the Riemannian metric M(theta) of RG.3 the kinetic energy is
#   K(theta, p) = c * sqrt( p^T M(theta)^{-1} p + m^2 c^2 ) + 0.5 * log det M,
# the relativistic energy of a particle of rest mass m moving on the statistical
# manifold whose local metric is M(theta) (the expected Fisher / SoftAbs of
# RG.3). The 0.5 log det M term is the SAME normaliser as the Gaussian RMHMC
# kinetic: after the change of variables p = L(theta) q (M = L L^T) the momentum
# marginal at fixed theta is the theta-independent integral of
# exp(-c sqrt(q^T q + m^2 c^2)), so the theta-marginal of the joint exp(-H) is
# exactly exp(log_prob). The sampler is therefore EXACT for any c and m -- they
# govern only efficiency, never correctness. As c -> infinity the kinetic energy
# reduces to the Gaussian RMHMC kinetic of RG.3 (the non-relativistic limit), so
# this level strictly generalises the Riemannian one.
#
# Finsler reading (the user's document section 12.3, organic-critical stance). A
# relativistic kinetic energy is the Legendre dual of a Finsler norm on
# velocities: F(theta, v) is NOT induced by an inner product, so the cost of
# motion is direction- and magnitude-dependent and the bounded speed c is the
# Finsler unit ball. The Randers form F = sqrt(g(v, v)) + beta(v) adds an
# asymmetric 1-form beta -- an irreversible / "magnetic" drift -- which makes the
# kinetic energy ODD in p and would break the p -> -p reversibility that the
# Metropolis correction relies on. It is deferred deliberately: an asymmetric
# beta models irreversible DYNAMICS (the document's ontological interest, the
# cost of raising complexity differing from lowering it), a different goal than
# exact SAMPLING of a fixed posterior, for which reversibility is required. The
# even, bounded relativistic kinetic realises exactly the part of the Finsler
# insight that serves robust, exact sampling.
#
# Non-separable Hamiltonian -> dedicated integrator. Because K depends on both
# theta (through M) and p, H is non-separable: the explicit leapfrog does not
# apply, and neither does the existing implicit leapfrog, whose drift sub-step
# hardwires the Gaussian velocity M^{-1} p. The dedicated integrator below is the
# generalised implicit leapfrog of Girolami & Calderhead (2011) for a GENERAL
# kinetic energy -- the same three reversible, volume-preserving sub-steps with
# the relativistic velocity grad_p K -- carried in the metric's `integrator`
# slot (the opt-in routing of gdpar_geom_hmc), leaving the default leapfrog
# bit-identical.

# Inverse-CDF sampler of the relativistic momentum radius. Under p = L q the
# radius r = ||q|| has density proportional to r^(d-1) exp(-c sqrt(r^2 + m^2 c^2))
# (the direction is uniform), independent of theta, so the sampler is built once.
.gdpar_geom_relativistic_radial <- function(d, speed, mass, n_grid = 8192L) {
  mc2 <- (mass * speed)^2
  logg <- function(r) (d - 1) * log(pmax(r, 1e-300)) - speed * sqrt(r^2 + mc2)
  opt <- stats::optimize(logg, lower = 1e-8, upper = 1e6, maximum = TRUE)
  peak <- opt$objective
  r_hi <- max(opt$maximum * 4, 1)
  while (logg(r_hi) > peak - 60 && r_hi < 1e9) r_hi <- r_hi * 2
  grid <- seq(0, r_hi, length.out = n_grid)
  lg <- logg(grid)
  lg[1] <- if (d > 1) -Inf else lg[2]            # r = 0: r^(d-1) vanishes for d > 1.
  w <- exp(lg - peak); w[!is.finite(w)] <- 0
  cdf <- cumsum(c(0, (w[-1] + w[-n_grid]) / 2 * diff(grid)))
  cdf <- cdf / cdf[n_grid]
  keep <- c(TRUE, diff(cdf) > 0)                 # strictly increasing for approx().
  inv <- stats::approxfun(cdf[keep], grid[keep], rule = 2, ties = "ordered")
  function(n) inv(stats::runif(n))
}

# Relativistic kinetic energy, its gradients and momentum sampler over a
# (position-dependent) mass metric M(theta) with rest mass m and speed c.
# s(theta, p) = sqrt(p^T M^{-1} p + m^2 c^2);  K = c s + 0.5 log det M.
.gdpar_geom_kinetic_relativistic <- function(metric, speed, mass) {
  mc2 <- (mass * speed)^2
  d <- metric$dim
  radial <- .gdpar_geom_relativistic_radial(d, speed, mass)
  list(
    value = function(theta, p) {
      Minv <- metric$inv_mass(theta)
      s <- sqrt(as.numeric(crossprod(p, Minv %*% p)) + mc2)
      speed * s + 0.5 * metric$logdet(theta)
    },
    # The velocity grad_p K = c M^{-1} p / s; its M-norm sqrt(v^T M v) = c
    # sqrt(p^T M^{-1} p) / s < c is bounded -- the tail-taming property.
    grad_p = function(theta, p) {
      Minv <- metric$inv_mass(theta)
      Minv_p <- as.numeric(Minv %*% p)
      s <- sqrt(as.numeric(crossprod(p, Minv_p)) + mc2)
      speed * Minv_p / s
    },
    # d K / d theta_i = -c (p^T M^{-1} dM_i M^{-1} p) / (2 s) + 0.5 tr(M^{-1} dM_i).
    grad_theta = function(theta, p) {
      if (!metric$position_dependent) return(rep(0, length(theta)))
      Minv <- metric$inv_mass(theta)
      Minv_p <- as.numeric(Minv %*% p)
      s <- sqrt(as.numeric(crossprod(p, Minv_p)) + mc2)
      dM <- metric$dmass(theta)
      vapply(seq_along(theta), function(i) {
        dMi <- dM[[i]]
        -speed * as.numeric(crossprod(Minv_p, dMi %*% Minv_p)) / (2 * s) +
          0.5 * sum(Minv * dMi)
      }, numeric(1))
    },
    # p = L q with q = r u, r from the relativistic radial law, u uniform on the
    # sphere -> p | theta ~ exp(-K(theta, .)) exactly (p^T M^{-1} p = r^2).
    draw_momentum = function(theta) {
      L <- metric$chol_mass(theta)
      z <- stats::rnorm(d)
      u <- z / sqrt(sum(z^2))
      as.numeric(L %*% (radial(1L) * u))
    }
  )
}

# Generalised implicit leapfrog (Girolami & Calderhead 2011) for the NON-separable
# relativistic Hamiltonian H = -log_prob + K(theta, p): three reversible sub-steps
# (implicit half kick in p, implicit drift in theta with the relativistic velocity,
# explicit half kick in p). Time-reversible and volume-preserving to the
# fixed-point tolerance; a non-converged solve or non-finite iterate is reported
# so the proposal is rejected rather than crashing the run.
.gdpar_geom_relativistic_integrator <- function(kinetic) {
  function(theta, p, target, eps, L, fp_tol = 1e-9, fp_max = 100L) {
    dH_dtheta <- function(th, pp)
      -target$grad_log_prob(th) + kinetic$grad_theta(th, pp)
    tryCatch({
      conv <- TRUE
      for (i in seq_len(L)) {
        # Sub-step 1: implicit half kick in p (fixed point in p_half).
        p_half <- p; ok1 <- FALSE
        for (it in seq_len(fp_max)) {
          p_new <- p - (eps / 2) * dH_dtheta(theta, p_half)
          delta <- max(abs(p_new - p_half))
          if (!is.finite(delta)) break
          p_half <- p_new
          if (delta < fp_tol) { ok1 <- TRUE; break }
        }
        # Sub-step 2: implicit drift in theta with the relativistic velocity,
        # theta' = theta + (eps/2)(grad_p K(theta, p_half) + grad_p K(theta', p_half)).
        v0 <- kinetic$grad_p(theta, p_half); ok2 <- FALSE
        theta_new <- theta
        for (it in seq_len(fp_max)) {
          th_new <- theta + (eps / 2) * (v0 + kinetic$grad_p(theta_new, p_half))
          delta <- max(abs(th_new - theta_new))
          if (!is.finite(delta)) break
          theta_new <- th_new
          if (delta < fp_tol) { ok2 <- TRUE; break }
        }
        # Sub-step 3: explicit half kick in p at the new position.
        p <- p_half - (eps / 2) * dH_dtheta(theta_new, p_half)
        theta <- theta_new
        if (!(ok1 && ok2) || any(!is.finite(theta)) || any(!is.finite(p))) {
          conv <- FALSE; break
        }
      }
      list(theta = theta, p = p, converged = conv)
    }, error = function(e) list(theta = theta, p = p, converged = FALSE))
  }
}

#' Finsler / relativistic metric for the geometric sampling engine
#'
#' Build the level-4 geometry of the Block RG hierarchy: a bounded,
#' non-Gaussian (relativistic) kinetic energy coupled to the position-dependent
#' Riemannian metric of \code{\link{gdpar_geom_metric_riemannian}}, the remedy
#' for heavy tails and directional anisotropy (the \code{G3_heavy_tails} target).
#' A Gaussian kinetic energy lets the velocity grow without bound, so a large
#' momentum in a heavy tail overshoots and the integrator must be tuned to the
#' stiffest region; the relativistic kinetic energy caps the velocity at a finite
#' speed, taming the tails and the ill-conditioning while staying exact.
#'
#' @details
#' The kinetic energy is the relativistic energy of a particle of rest mass
#' \eqn{m} on the statistical manifold with local metric \eqn{M(\theta)} (the
#' expected Fisher or the SoftAbs of the observed Hessian, supplied through the
#' RG.3 machinery):
#' \deqn{K(\theta, p) = c\,\sqrt{p^\top M(\theta)^{-1} p + m^2 c^2} +
#'   \tfrac{1}{2}\log\det M(\theta).}
#' The velocity is \eqn{\nabla_p K = c\,M^{-1}p/\sqrt{p^\top M^{-1}p + m^2 c^2}},
#' whose \eqn{M}-norm is strictly below the speed \eqn{c} for every momentum --
#' the bounded-velocity property that tames heavy tails. Three properties hold by
#' construction:
#' \itemize{
#'   \item \strong{Exactness independent of \eqn{c} and \eqn{m}}: the
#'     \eqn{\tfrac{1}{2}\log\det M} term is the same normaliser as the Gaussian
#'     Riemannian kinetic, so the \eqn{\theta}-marginal of the joint
#'     \eqn{e^{-H}} is exactly \eqn{e^{\log p(\theta)}}; the kinetic energy is a
#'     preconditioner, not part of the target, and the Metropolis correction with
#'     the exact density keeps the sampler exact for any \eqn{c} and \eqn{m}
#'     (they govern only efficiency).
#'   \item \strong{Non-relativistic limit}: as \eqn{c \to \infty} the kinetic
#'     energy reduces to the Gaussian Riemannian kinetic of
#'     \code{\link{gdpar_geom_metric_riemannian}}, so this level strictly
#'     generalises the Riemannian one. Larger \code{speed} approaches that limit
#'     (less tail-taming, faster bulk mixing); smaller \code{speed} caps the
#'     velocity sooner (more robust tails, slower bulk).
#'   \item \strong{Finsler structure} (the user's document section 12.3): a
#'     relativistic kinetic energy is the Legendre dual of a Finsler norm on
#'     velocities, a norm not induced by an inner product, so the cost of motion
#'     depends on direction as well as magnitude (the speed \eqn{c} is the Finsler
#'     unit ball). The asymmetric Randers extension
#'     \eqn{F = \sqrt{g(v, v)} + \beta(v)} (an irreversible drift) is deliberately
#'     not included: it makes the kinetic energy odd in \eqn{p} and would break
#'     the reversibility the Metropolis correction relies on; it models
#'     irreversible dynamics rather than exact sampling of a fixed target.
#' }
#' Because \eqn{K} depends on both \eqn{\theta} (through \eqn{M}) and \eqn{p}, the
#' Hamiltonian is non-separable. Sampling therefore uses a \strong{dedicated}
#' generalised implicit leapfrog (Girolami & Calderhead 2011) carried in the
#' metric's \code{integrator} slot -- three reversible, volume-preserving
#' sub-steps with the relativistic velocity -- which \code{\link{gdpar_geom_hmc}}
#' runs in place of the default leapfrog, leaving the default branch bit-identical.
#' The momentum is refreshed from the exact relativistic momentum law (an
#' inverse-CDF radial sampler under \eqn{p = L q}), so every momentum draw is
#' from \eqn{e^{-K(\theta, \cdot)}}.
#'
#' @param target A \code{\link{gdpar_geom_target}} (or an object accepted by it).
#'   Supplies the dimension and, for \code{curvature = "softabs"}, the Hessian /
#'   gradient of the underlying Riemannian metric.
#' @param curvature Curvature source of the underlying Riemannian mass:
#'   \code{"fisher"} (expected Fisher, supplied through \code{fisher}) or
#'   \code{"softabs"} (SoftAbs of the observed Hessian). See
#'   \code{\link{gdpar_geom_metric_riemannian}}.
#' @param fisher,dfisher For \code{curvature = "fisher"}: the expected Fisher
#'   matrix function \code{fisher(theta)} and optionally its derivative
#'   \code{dfisher(theta)} (a length-\code{dim} list of partials); finite
#'   differenced from \code{fisher} when \code{dfisher} is \code{NULL}.
#' @param speed Positive speed of light \eqn{c} bounding the velocity. Larger
#'   values approach the Gaussian Riemannian kinetic (less tail-taming); smaller
#'   values cap the velocity sooner (more robust). Tune to the momentum scale.
#' @param rest_mass Positive rest mass \eqn{m}; with \code{speed} it sets the
#'   bulk curvature of the kinetic energy (the non-relativistic limit is
#'   \code{speed} or \code{rest_mass} to infinity).
#' @param alpha,floor,fd_step SoftAbs softening, eigenvalue floor and
#'   finite-difference step of the underlying Riemannian metric.
#'
#' @return A list of class \code{gdpar_geom_metric} with
#'   \code{position_dependent = TRUE}, \code{metric_kind = "relativistic"}, the
#'   metric interface (\code{mass}, \code{inv_mass}, \code{chol_mass},
#'   \code{logdet}, \code{dmass}) of the underlying Riemannian mass, a
#'   \code{kinetic} object (the relativistic \code{value}, \code{grad_p},
#'   \code{grad_theta} and \code{draw_momentum}) and an \code{integrator} closure,
#'   both consumed by \code{\link{gdpar_geom_hmc}}, plus the fields \code{speed},
#'   \code{rest_mass} and \code{curvature}.
#'
#' @references
#' Lu, X., Perrone, V., Hasenclever, L., Teh, Y. W. and Vollmer, S. (2017)
#' Relativistic Monte Carlo. \emph{AISTATS} 54, 1236--1245.
#'
#' Livingstone, S., Faulkner, M. F. and Roberts, G. O. (2019) Kinetic energy
#' choice in Hamiltonian/hybrid Monte Carlo. \emph{Biometrika} 106, 303--319.
#'
#' Girolami, M. and Calderhead, B. (2011) Riemann manifold Langevin and
#' Hamiltonian Monte Carlo methods. \emph{JRSS-B} 73, 123--214.
#'
#' Randers, G. (1941) On an asymmetrical metric in the four-space of general
#' relativity. \emph{Physical Review} 59, 195--199.
#'
#' @seealso \code{\link{gdpar_geom_metric_riemannian}},
#'   \code{\link{gdpar_geom_metric_subriemannian}}, \code{\link{gdpar_geom_hmc}}.
#'
#' @examples
#' # A two-dimensional Student-t (heavy tails); the expected Fisher of the
#' # independent t is diagonal. The bounded kinetic energy samples the tails
#' # without the overshoot of a Gaussian kinetic.
#' nu <- 2
#' tgt <- gdpar_geom_target(
#'   log_prob = function(theta) -((nu + 1) / 2) * sum(log1p(theta^2 / nu)),
#'   grad_log_prob = function(theta) -(nu + 1) * theta / (nu + theta^2), dim = 2)
#' fisher <- function(theta) diag((nu + 1) / (nu + 3), 2)   # expected Fisher of t.
#' metric <- gdpar_geom_metric_relativistic(tgt, fisher = fisher, speed = 5)
#' fit <- gdpar_geom_hmc(tgt, metric = metric, epsilon = 0.4, L = 10,
#'                       n_iter = 200, n_warmup = 100, seed = 1)
#' fit$ebfmi
#'
#' @export
gdpar_geom_metric_relativistic <- function(target,
                                           curvature = c("fisher", "softabs"),
                                           fisher = NULL, dfisher = NULL,
                                           speed = 10, rest_mass = 1,
                                           alpha = 1e6, floor = 1e-8,
                                           fd_step = 1e-4) {
  curvature <- match.arg(curvature)
  if (!inherits(target, "gdpar_geom_target")) target <- gdpar_geom_target(target)
  d <- target$dim
  assert_numeric_scalar(speed, "speed", lower = 0)
  assert_numeric_scalar(rest_mass, "rest_mass", lower = 0)
  if (speed <= 0 || rest_mass <= 0) {
    gdpar_abort("'speed' and 'rest_mass' must be strictly positive.",
                class = "gdpar_input_error")
  }
  base <- gdpar_geom_metric_riemannian(target, curvature = curvature,
                                       fisher = fisher, dfisher = dfisher,
                                       alpha = alpha, floor = floor,
                                       fd_step = fd_step)
  kinetic <- .gdpar_geom_kinetic_relativistic(base, speed, rest_mass)
  integrator <- .gdpar_geom_relativistic_integrator(kinetic)
  obj <- list(
    position_dependent = TRUE, dim = d, metric_kind = "relativistic",
    speed = speed, rest_mass = rest_mass, curvature = curvature,
    mass = base$mass, inv_mass = base$inv_mass, chol_mass = base$chol_mass,
    logdet = base$logdet, dmass = base$dmass,
    kinetic = kinetic, integrator = integrator)
  class(obj) <- c("gdpar_geom_metric", "list")
  obj
}

# Gaussian kinetic energy K(theta, p) = 0.5 p^T M(theta)^{-1} p + 0.5 logdet M.
# For a constant Euclidean metric grad_theta K = 0; the position-dependent term
# (nonzero for the Riemannian level) is supplied by RG.3.
.gdpar_geom_kinetic_gaussian <- function(metric) {
  list(
    value = function(theta, p) {
      Minv <- metric$inv_mass(theta)
      0.5 * as.numeric(crossprod(p, Minv %*% p)) + 0.5 * metric$logdet(theta)
    },
    grad_p = function(theta, p) as.numeric(metric$inv_mass(theta) %*% p),
    grad_theta = function(theta, p) {
      if (!metric$position_dependent) return(rep(0, length(theta)))
      # d K / d theta_i = 0.5 tr(M^{-1} dM_i) - 0.5 p^T M^{-1} dM_i M^{-1} p.
      Minv <- metric$inv_mass(theta)
      dM <- metric$dmass(theta)
      Minv_p <- as.numeric(Minv %*% p)
      vapply(seq_along(theta), function(i) {
        dMi <- dM[[i]]
        0.5 * sum(Minv * dMi) -
          0.5 * as.numeric(crossprod(Minv_p, dMi %*% Minv_p))
      }, numeric(1))
    },
    draw_momentum = function(theta) {
      L <- metric$chol_mass(theta)
      as.numeric(L %*% stats::rnorm(nrow(L)))
    }
  )
}

# Gradient of the Hamiltonian in theta: dH/dtheta = -grad_log_prob + grad_theta K
# (the second term is zero for a constant Euclidean metric).
.gdpar_geom_dH_dtheta <- function(target, kinetic, theta, p) {
  -target$grad_log_prob(theta) + kinetic$grad_theta(theta, p)
}

# One leapfrog step. For a constant Euclidean metric this is the explicit,
# symplectic, reversible, second-order Stoermer--Verlet step. For a
# position-dependent (Riemannian) metric it is the generalised implicit leapfrog
# of Girolami & Calderhead (2011): two fixed-point sub-steps (in the momentum,
# then in the position) bracket an explicit momentum sub-step. The integrator is
# exactly time-reversible and volume-preserving up to the fixed-point tolerance;
# non-convergence of either solve is reported so the proposal is rejected.
.gdpar_geom_leapfrog_step <- function(theta, p, target, metric, kinetic, eps,
                                      fp_tol = 1e-9, fp_max = 100L) {
  if (!metric$position_dependent) {
    p <- p + (eps / 2) * target$grad_log_prob(theta)
    theta <- theta + eps * kinetic$grad_p(theta, p)
    p <- p + (eps / 2) * target$grad_log_prob(theta)
    return(list(theta = theta, p = p))
  }
  # The three implicit sub-steps. Any numerical failure -- a step too large for
  # the local curvature: a non-finite iterate, a Hessian that overflows, a
  # metric that is not positive-definite -- is caught and reported as
  # non-convergence so the proposal is rejected rather than crashing the run.
  tryCatch({
    # Sub-step 1: implicit half kick in p (fixed point in p_half).
    p_half <- p; ok1 <- FALSE
    for (it in seq_len(fp_max)) {
      dH <- .gdpar_geom_dH_dtheta(target, kinetic, theta, p_half)
      p_new <- p - (eps / 2) * dH
      delta <- max(abs(p_new - p_half))
      if (!is.finite(delta)) break
      p_half <- p_new
      if (delta < fp_tol) { ok1 <- TRUE; break }
    }
    # Sub-step 2: implicit full drift in theta (fixed point in theta_new),
    # theta_new = theta + (eps/2)(M(theta)^{-1} + M(theta_new)^{-1}) p_half.
    drift0 <- as.numeric(metric$inv_mass(theta) %*% p_half); ok2 <- FALSE
    theta_new <- theta
    for (it in seq_len(fp_max)) {
      th_new <- theta + (eps / 2) *
        (drift0 + as.numeric(metric$inv_mass(theta_new) %*% p_half))
      delta <- max(abs(th_new - theta_new))
      if (!is.finite(delta)) break
      theta_new <- th_new
      if (delta < fp_tol) { ok2 <- TRUE; break }
    }
    # Sub-step 3: explicit half kick in p at the new position.
    p_final <- p_half - (eps / 2) *
      .gdpar_geom_dH_dtheta(target, kinetic, theta_new, p_half)
    list(theta = theta_new, p = p_final, converged = ok1 && ok2)
  }, error = function(e) list(theta = theta, p = p, converged = FALSE))
}

# Integrate L leapfrog steps; return the endpoint state and whether every
# implicit solve converged (always TRUE for the explicit Euclidean step).
.gdpar_geom_leapfrog_traj <- function(theta, p, target, metric, kinetic,
                                      eps, L, fp_tol = 1e-9, fp_max = 100L) {
  conv <- TRUE
  for (i in seq_len(L)) {
    st <- .gdpar_geom_leapfrog_step(theta, p, target, metric, kinetic, eps,
                                    fp_tol, fp_max)
    theta <- st$theta; p <- st$p
    if (isFALSE(st$converged)) { conv <- FALSE; break }
  }
  list(theta = theta, p = p, converged = conv)
}

# Hamiltonian H = U + K = -log_prob + K.
.gdpar_geom_hamiltonian <- function(target, kinetic, theta, p) {
  -target$log_prob(theta) + kinetic$value(theta, p)
}

#' Static Hamiltonian Monte Carlo with a pluggable geometry (Block RG engine)
#'
#' Sample a target with the R-native geometric engine of decision A: a fixed
#' step-size, fixed trajectory-length Hamiltonian Monte Carlo with a
#' Metropolis correction, over a pluggable metric. With the default Euclidean
#' metric this is textbook HMC; the higher levels of the Block RG hierarchy
#' (Riemannian, Finsler / relativistic, sub-Riemannian) reuse the same loop
#' with a richer metric and integrator. This is the validated Euclidean
#' scaffolding for those levels, not a replacement for the package's cmdstan
#' fit path.
#'
#' @param target A \code{\link{gdpar_geom_target}}, or a list/closure accepted
#'   by it.
#' @param metric A \code{\link{gdpar_geom_metric_euclidean}} (or compatible)
#'   metric. Defaults to the identity Euclidean metric of the target dimension.
#' @param epsilon Numeric leapfrog step size.
#' @param L Integer number of leapfrog steps per proposal.
#' @param n_iter Integer number of retained iterations.
#' @param n_warmup Integer number of warmup iterations discarded from the
#'   returned draws (no adaptation is performed; warmup only burns in).
#' @param init Optional numeric vector of length \code{target$dim} giving the
#'   initial position. Defaults to zeros.
#' @param seed Optional integer seed. If supplied, the RNG state is set and
#'   restored around the run.
#' @param fp_tol,fp_max Convergence tolerance and maximum iteration count for the
#'   fixed-point solves of the implicit generalised leapfrog (used only when the
#'   metric is position-dependent). A proposal whose solves do not converge is
#'   counted as divergent and rejected.
#'
#' @return A list of class \code{gdpar_geom_hmc} with \code{draws} (an
#'   \code{n_iter x dim} matrix), \code{accept_rate}, \code{n_divergent}
#'   (proposals with a non-finite Hamiltonian, a large energy error, or a
#'   non-converged implicit solve), \code{energy} (the per-iteration Hamiltonian
#'   energy trace), \code{ebfmi} (the energy Bayesian fraction of missing
#'   information; higher is better, low values flag a metric ill-matched to the
#'   geometry), \code{epsilon}, \code{L} and \code{metric_type}.
#'
#' @seealso \code{\link{gdpar_geom_target}},
#'   \code{\link{gdpar_geom_metric_euclidean}},
#'   \code{\link{gdpar_geometry_diagnostic}}.
#'
#' @examples
#' tgt <- gdpar_geom_target(
#'   log_prob = function(theta) -0.5 * sum(theta^2),
#'   grad_log_prob = function(theta) -theta, dim = 2)
#' fit <- gdpar_geom_hmc(tgt, epsilon = 0.3, L = 12, n_iter = 200,
#'                       n_warmup = 100, seed = 1)
#' colMeans(fit$draws)
#'
#' @export
gdpar_geom_hmc <- function(target, metric = NULL, epsilon = 0.1, L = 20L,
                           n_iter = 1000L, n_warmup = 500L, init = NULL,
                           seed = NULL, fp_tol = 1e-9, fp_max = 100L) {
  if (!inherits(target, "gdpar_geom_target")) target <- gdpar_geom_target(target)
  assert_numeric_scalar(epsilon, "epsilon", lower = 0)
  assert_count(L, "L")
  assert_count(n_iter, "n_iter")
  if (!is.numeric(n_warmup) || length(n_warmup) != 1L || n_warmup < 0) {
    gdpar_abort("Argument 'n_warmup' must be a non-negative integer scalar.",
                class = "gdpar_input_error")
  }
  d <- target$dim
  if (is.null(metric)) metric <- gdpar_geom_metric_euclidean(dim = d)
  if (!is.null(seed)) {
    assert_count(seed, "seed")
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv)
      on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
    }
    set.seed(seed)
  }
  # A metric may carry its own kinetic energy (the relativistic kinetic of RG.4);
  # otherwise the Gaussian kinetic is used. Every existing metric leaves
  # metric$kinetic NULL, so the default branch stays bit-identical.
  kinetic <- metric$kinetic %||% .gdpar_geom_kinetic_gaussian(metric)

  theta <- if (is.null(init)) rep(0, d) else as.numeric(init)
  if (length(theta) != d) {
    gdpar_abort(sprintf("'init' must have length %d.", d),
                class = "gdpar_input_error")
  }
  total <- as.integer(n_warmup) + as.integer(n_iter)
  draws <- matrix(NA_real_, nrow = n_iter, ncol = d)
  colnames(draws) <- target$param_names
  energy <- numeric(n_iter)
  n_accept <- 0L; n_div <- 0L; kept <- 0L
  for (it in seq_len(total)) {
    p0 <- kinetic$draw_momentum(theta)
    H0 <- .gdpar_geom_hamiltonian(target, kinetic, theta, p0)
    # A proposal that drives the target undefined (a cmdstan log-density that
    # throws on a non-finite unconstrained value, an overflowing gradient) is
    # caught and rejected rather than crashing the run -- the same hardening as
    # the implicit leapfrog. Stable proposals never trigger the handler, so the
    # default branch stays bit-identical.
    prop <- tryCatch(
      if (is.function(metric$integrator)) {
        # A metric carrying its own integrator (the sub-Riemannian splitting or
        # the relativistic generalised leapfrog of RG.4) overrides the default
        # leapfrog; the Metropolis correction with the exact Hamiltonian below --
        # using the metric's own kinetic when present -- keeps the sampler exact.
        metric$integrator(theta, p0, target, epsilon, L, fp_tol, fp_max)
      } else {
        .gdpar_geom_leapfrog_traj(theta, p0, target, metric, kinetic,
                                  epsilon, L, fp_tol, fp_max)
      },
      error = function(e) list(theta = theta, p = p0, converged = FALSE))
    H1 <- tryCatch(
      .gdpar_geom_hamiltonian(target, kinetic, prop$theta, prop$p),
      error = function(e) Inf)
    dH <- H1 - H0
    divergent <- isFALSE(prop$converged) || !is.finite(dH) || abs(dH) > 1000
    if (divergent) n_div <- n_div + 1L
    accept <- !divergent && (log(stats::runif(1)) < -dH)
    if (accept) { theta <- prop$theta; n_accept <- n_accept + 1L }
    if (it > n_warmup) {
      kept <- kept + 1L; draws[kept, ] <- theta; energy[kept] <- H0
    }
  }
  obj <- list(
    draws = draws, accept_rate = n_accept / total, n_divergent = n_div,
    energy = energy, ebfmi = .gdpar_geom_ebfmi(energy),
    epsilon = epsilon, L = as.integer(L),
    metric_type = if (metric$position_dependent) "position_dependent"
                  else "euclidean_constant")
  class(obj) <- c("gdpar_geom_hmc", "list")
  obj
}

# Energy Bayesian fraction of missing information: the ratio of the mean squared
# successive energy difference to the energy variance. Low values warn that the
# metric resolves the geometry poorly (Betancourt 2016).
.gdpar_geom_ebfmi <- function(energy) {
  energy <- energy[is.finite(energy)]
  if (length(energy) < 2L) return(NA_real_)
  den <- sum((energy - mean(energy))^2)
  if (!is.finite(den) || den <= 0) return(NA_real_)
  sum(diff(energy)^2) / den
}

#' Print method for gdpar_geom_hmc objects
#'
#' @param x A \code{gdpar_geom_hmc}.
#' @param ... Unused.
#' @return Invisibly returns \code{x}.
#' @export
print.gdpar_geom_hmc <- function(x, ...) {
  cat("<gdpar_geom_hmc> ", nrow(x$draws), " draws x ", ncol(x$draws),
      " dims\n", sep = "")
  cat("  metric: ", x$metric_type, " | epsilon: ", x$epsilon, " | L: ", x$L,
      "\n", sep = "")
  cat("  accept rate: ", format(x$accept_rate, digits = 3),
      " | divergent: ", x$n_divergent,
      " | E-BFMI: ", format(x$ebfmi, digits = 3), "\n", sep = "")
  invisible(x)
}

#' Adaptive Riemannian HMC with online novelty-driven active learning
#'
#' Drive the geometric sampler with the learned expected-Fisher metric
#' (\code{\link{gdpar_geom_metric_gp_fisher}}) while refining that metric online
#' between trajectories -- the full novelty / active-learning loop that completes
#' the Riemannian level of the Block RG hierarchy. It realises the two-phase,
#' decoupled-archivist design of ORPHEUS-PIMC-A (section 16) with a correctness
#' improvement specific to Riemannian sampling.
#'
#' @details
#' The metric is held \strong{fixed within each trajectory} (preserving the
#' reversibility of the implicit generalised leapfrog) and is re-learned only
#' \strong{between rounds} from a growing reservoir; the number of new sites
#' admitted shrinks geometrically per round (a decreasing, Robbins--Monro-style
#' schedule) so the metric sequence settles. A final sampling phase with the
#' \strong{frozen} metric produces the returned draws.
#'
#' One adaptation round:
#' \enumerate{
#'   \item sample a batch of trajectories with the current frozen metric
#'     (continuing the chain from the previous round);
#'   \item score the kept positions by the surrogate's epistemic
#'     \code{novelty} (predictive standard deviation); positions above
#'     \code{novelty_tol} are candidate new reservoir sites;
#'   \item admit at most \code{ceiling(n_add * decay^(round - 1))} of the most
#'     novel candidates and re-learn the surrogate on the augmented reservoir.
#' }
#' The loop stops when the batch's maximum novelty falls below \code{novelty_tol}
#' (the reservoir covers the typical set) or \code{max_rounds} is reached.
#'
#' \strong{Correctness vs efficiency} (the honesty convention of ORPHEUS-PIMC-A
#' section 16.3): the sampler is \emph{exact} in every phase regardless of the
#' metric -- the metric is a preconditioner, not part of the target, so the
#' Metropolis correction with the exact density is the corrector, and no delayed
#' acceptance is needed (the improvement over ORPHEUS, where the surrogate enters
#' the acceptance). What the active learning buys is \emph{efficiency}, which is
#' measured (E-BFMI, acceptance, novelty trace), not asserted.
#'
#' @param target A generative \code{\link{gdpar_geom_target}} carrying
#'   \code{simulate} and \code{score}.
#' @param fisher Optional \code{fisher(theta)} function; defaults to the
#'   simulation-based estimator \code{\link{gdpar_geom_fisher_simulator}} built
#'   from \code{target} with \code{n_sim} simulations.
#' @param n_sim Simulations per Fisher evaluation for the default estimator.
#' @param n_sites_init Initial reservoir size (cold-start warmup draws).
#' @param max_rounds Maximum number of adaptation rounds.
#' @param batch Trajectories sampled per adaptation round.
#' @param n_add Maximum new sites admitted in the first round.
#' @param decay Geometric decay in \eqn{(0, 1]} of admitted sites per round (the
#'   decreasing-adaptation schedule).
#' @param novelty_tol Novelty (predictive standard deviation) threshold, used
#'   both to admit new sites and to stop the loop.
#' @param warmup_metric Optional metric for the cold-start reservoir run;
#'   defaults to a SoftAbs Riemannian metric (good neck coverage on curved
#'   targets). Pass \code{gdpar_geom_metric_euclidean(dim = target$dim)} for a
#'   cheaper warmup.
#' @param epsilon,L Leapfrog step and trajectory length for all phases.
#' @param n_iter,n_warmup Retained and discarded iterations of the final phase.
#' @param init Optional initial position (defaults to zeros).
#' @param seed Optional integer seed for the whole loop.
#' @param alpha,nugget,lengthscale Surrogate hyperparameters forwarded to
#'   \code{\link{gdpar_geom_metric_gp_fisher}}.
#'
#' @return A list of class \code{gdpar_geom_rmhmc_adaptive} with \code{draws}
#'   (final phase), the final \code{metric}, the \code{reservoir} sites,
#'   \code{n_sites_trace}, \code{novelty_trace}, \code{accept_rate},
#'   \code{ebfmi}, \code{n_divergent}, \code{n_rounds}, \code{epsilon} and
#'   \code{L}.
#'
#' @seealso \code{\link{gdpar_geom_fisher_simulator}},
#'   \code{\link{gdpar_geom_metric_gp_fisher}}, \code{\link{gdpar_geom_hmc}}.
#'
#' @examples
#' # Bivariate normal location model; the loop learns the constant expected
#' # Fisher and samples with the frozen metric.
#' Sigma0 <- diag(c(1, 4))
#' P0 <- solve(Sigma0)
#' tgt <- gdpar_geom_target(
#'   log_prob = function(theta) -0.5 * sum((theta^2) / diag(Sigma0)),
#'   grad_log_prob = function(theta) -theta / diag(Sigma0),
#'   hessian = function(theta) -P0, dim = 2,
#'   simulate = function(theta)
#'     as.numeric(theta + sqrt(diag(Sigma0)) * stats::rnorm(2)),
#'   score = function(theta, y) as.numeric(P0 %*% (y - theta)))
#' fit <- gdpar_geom_rmhmc_adaptive(tgt, n_sim = 80, n_sites_init = 8,
#'                                  max_rounds = 1, batch = 10, L = 8,
#'                                  n_iter = 20, n_warmup = 10, seed = 1)
#' fit$n_rounds
#'
#' @export
gdpar_geom_rmhmc_adaptive <- function(target, fisher = NULL, n_sim = 64L,
                                      n_sites_init = 30L, max_rounds = 5L,
                                      batch = 60L, n_add = 20L, decay = 0.5,
                                      novelty_tol = 0.2, warmup_metric = NULL,
                                      epsilon = 0.05, L = 25L, n_iter = 500L,
                                      n_warmup = 200L, init = NULL, seed = NULL,
                                      alpha = 1e6, nugget = 1e-4,
                                      lengthscale = NULL) {
  if (!inherits(target, "gdpar_geom_target")) target <- gdpar_geom_target(target)
  d <- target$dim
  assert_count(n_sites_init, "n_sites_init")
  assert_count(max_rounds, "max_rounds")
  assert_count(batch, "batch")
  assert_count(n_add, "n_add")
  assert_numeric_scalar(decay, "decay", lower = 0)
  assert_numeric_scalar(novelty_tol, "novelty_tol", lower = 0)
  if (n_sites_init < 2L) {
    gdpar_abort("'n_sites_init' must be at least two.",
                class = "gdpar_input_error")
  }
  if (is.null(fisher)) {
    fisher <- gdpar_geom_fisher_simulator(target, n_sim = n_sim,
                                          seed = (seed %||% 1L))
  } else if (!is.function(fisher)) {
    gdpar_abort("'fisher' must be a function of theta or NULL.",
                class = "gdpar_input_error")
  }
  if (!is.null(seed)) {
    assert_count(seed, "seed")
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv)
      on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
    }
    set.seed(seed)
  }
  if (is.null(warmup_metric)) {
    warmup_metric <- gdpar_geom_metric_riemannian(target, curvature = "softabs",
                                                  alpha = alpha)
  }
  build_metric <- function(sites) {
    gdpar_geom_metric_gp_fisher(target, fisher = fisher, sites = sites,
                                alpha = alpha, nugget = nugget,
                                lengthscale = lengthscale)
  }

  # Phase 1: cold-start reservoir + first surrogate.
  sites <- gdpar_geom_reservoir(target, n_sites = n_sites_init,
                                metric = warmup_metric, epsilon = epsilon, L = L,
                                n_warmup = n_warmup, init = init)
  metric <- build_metric(sites)
  n_sites_trace <- nrow(sites)
  novelty_trace <- numeric(0)
  theta0 <- if (is.null(init)) rep(0, d) else as.numeric(init)
  rounds_done <- 0L

  # Phase 2: novelty-driven adaptation with decreasing admission.
  for (r in seq_len(max_rounds)) {
    fit_b <- gdpar_geom_hmc(target, metric = metric, epsilon = epsilon, L = L,
                            n_iter = batch, n_warmup = 0L, init = theta0)
    draws_b <- fit_b$draws
    last <- draws_b[nrow(draws_b), ]
    if (all(is.finite(last))) theta0 <- last
    nov <- apply(draws_b, 1L, metric$novelty)
    max_nov <- max(nov[is.finite(nov)], 0)
    novelty_trace <- c(novelty_trace, max_nov)
    rounds_done <- r
    if (max_nov < novelty_tol) break
    n_admit <- ceiling(n_add * decay^(r - 1))
    eligible <- which(nov >= novelty_tol)
    if (length(eligible) == 0L) break
    take <- eligible[order(nov[eligible], decreasing = TRUE)]
    take <- take[seq_len(min(n_admit, length(take)))]
    sites <- rbind(sites, draws_b[take, , drop = FALSE])
    metric <- build_metric(sites)
    n_sites_trace <- c(n_sites_trace, nrow(sites))
  }

  # Phase 3: final sampling with the frozen metric.
  final <- gdpar_geom_hmc(target, metric = metric, epsilon = epsilon, L = L,
                          n_iter = n_iter, n_warmup = n_warmup, init = theta0)
  obj <- list(
    draws = final$draws, metric = metric, reservoir = sites,
    n_sites_trace = n_sites_trace, novelty_trace = novelty_trace,
    accept_rate = final$accept_rate, ebfmi = final$ebfmi,
    n_divergent = final$n_divergent, n_rounds = rounds_done,
    epsilon = epsilon, L = as.integer(L))
  class(obj) <- c("gdpar_geom_rmhmc_adaptive", "list")
  obj
}

#' Print method for gdpar_geom_rmhmc_adaptive objects
#'
#' @param x A \code{gdpar_geom_rmhmc_adaptive}.
#' @param ... Unused.
#' @return Invisibly returns \code{x}.
#' @export
print.gdpar_geom_rmhmc_adaptive <- function(x, ...) {
  nt <- x$novelty_trace
  cat("<gdpar_geom_rmhmc_adaptive> ", nrow(x$draws), " draws x ",
      ncol(x$draws), " dims\n", sep = "")
  cat("  rounds: ", x$n_rounds, " | reservoir sites: ", nrow(x$reservoir),
      " | epsilon: ", x$epsilon, " | L: ", x$L, "\n", sep = "")
  cat("  accept rate: ", format(x$accept_rate, digits = 3),
      " | divergent: ", x$n_divergent,
      " | E-BFMI: ", format(x$ebfmi, digits = 3),
      " | final novelty: ",
      format(if (length(nt)) nt[length(nt)] else NA_real_, digits = 3),
      "\n", sep = "")
  invisible(x)
}

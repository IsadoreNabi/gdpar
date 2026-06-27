#' Joint kernel Stein discrepancy between EB and FB posteriors
#'
#' Sub-bloque 9.3.c (Block 9, Session B9.4, 2026-05-27) under
#' canonized decision H.iv lateral. Operationalizes the open question
#' recorded in the Roxygen of \code{\link{gdpar_compare_eb_fb}} that
#' the marginal TV reported by the latter is only a coarse proxy of
#' the distributional discrepancy between the EB and FB posteriors,
#' and that the joint discrepancy over \eqn{\xi = (a, b, W,
#' \mathrm{dispersion})} deserves a density-free spectral metric. The
#' canonical choice is the kernel Stein discrepancy (KSD) of Gorham
#' and Mackey (2017) "Measuring Sample Quality with Kernels", JMLR
#' 18(196):1-72; Liu, Lee, and Jordan (2016) "A Kernelized Stein
#' Discrepancy for Goodness-of-Fit Tests", ICML.
#'
#' \strong{Target choice in this iteration.} The KSD requires the
#' target distribution's score function \eqn{s_p(x) = \nabla \log
#' p(x)}. For tractability and atomicity within Session B9.4, the
#' implementation uses an empirical Gaussian target derived from the
#' FB posterior draws via the empirical mean \eqn{\hat\mu} and
#' covariance \eqn{\hat\Sigma} (a Laplace approximation of the FB
#' target): \eqn{s(x) = -\hat\Sigma^{-1} (x - \hat\mu)}. The
#' full-KSD variant that uses the FB Stan model's exact gradient via
#' \code{cmdstanr}'s \code{grad_log_prob()} (true one-sample KSD of
#' EB samples against the actual FB target) is a documented
#' extension for B9.x.
#'
#' \strong{Base kernel.} Inverse multi-quadric (IMQ) of Gorham-Mackey
#' \eqn{k(x, y) = (h + \|x - y\|^2)^\beta}, \eqn{\beta \in (-1, 0)},
#' default \eqn{\beta = -1/2}; the bandwidth \eqn{h} (in squared
#' units of \eqn{x}) defaults to the \emph{median heuristic}: the
#' median of squared pairwise distances between FB draws, which is
#' dimension-adaptive. The RBF (Gaussian) kernel
#' \eqn{k(x, y) = \exp(-\|x - y\|^2 / (2 h))} is provided as a
#' textbook alternative (Liu, Lee, Jordan 2016).
#'
#' \strong{Stein kernel under a Gaussian target.} With \eqn{s(x) =
#' -\hat\Sigma^{-1}(x - \hat\mu)}, the canonical Stein kernel is
#' \deqn{k_p(x, y) = \langle s(x), s(y) \rangle k(x, y) +
#'                   \langle s(x), \nabla_y k(x, y) \rangle +
#'                   \langle s(y), \nabla_x k(x, y) \rangle +
#'                   \mathrm{tr}(\nabla_x \nabla_y^\top k(x, y)),}
#' and the KSD V-statistic is
#' \deqn{\mathrm{KSD}(Q, P) = \sqrt{\max\{0,\,
#'   n^{-2} \sum_{i, j = 1}^{n} k_p(x_i, x_j)\}},
#'   \qquad x_i \sim Q.}
#' A value close to zero indicates that the EB posterior \eqn{Q}
#' matches the empirical Gaussian fit of the FB target \eqn{P}; a
#' positive value indicates joint distributional deviation. The
#' marginal TV reported by \code{\link{gdpar_compare_eb_fb}} may be
#' small even when the joint KSD detects a deviation in the joint
#' dependence structure; the two diagnostics are complementary.
#'
#' \strong{ESS-weighted variant.} When \code{ess_weighted = TRUE},
#' both posteriors are thinned to a common count
#' \eqn{\min\{ \widehat{\mathrm{ESS}}_{EB},
#' \widehat{\mathrm{ESS}}_{FB}\}}
#' via uniform subsampling. This mitigates the autocorrelation bias
#' on the V-statistic when the MCMC chains are short relative to the
#' integrated autocorrelation time; per-variable basic ESS is
#' computed via \code{posterior::ess_basic} and the minimum across
#' retained \eqn{\xi}-variables is used.
#'
#' \strong{Cross-reference.} Pedagogical interpretation of the Stein
#' operator and of the asymmetry between the score-based KSD and the
#' density-free TV is documented in the vignette
#' \code{vignette("v07b_eb_multivariate", package = "gdpar")},
#' Section 11.
#'
#' @param eb_fit An object of class \code{gdpar_eb_fit} produced by
#'   \code{\link{gdpar_eb}}.
#' @param fb_fit An object of class \code{gdpar_fit} produced by
#'   \code{\link{gdpar}}. Must have been fitted on the same dataset
#'   as \code{eb_fit} (same outcome, same covariates, same \eqn{K /
#'   p} regime).
#' @param kernel Character scalar; one of \code{"imq"} (default) or
#'   \code{"rbf"}.
#' @param bandwidth Character scalar; one of \code{"median"} (default,
#'   classic median heuristic on FB squared pairwise distances) or
#'   \code{"fixed"} (use \code{bandwidth_value}).
#' @param bandwidth_value Numeric scalar; bandwidth value (squared
#'   units of \eqn{x}) when \code{bandwidth = "fixed"}. Defaults to
#'   1.0 if \code{NULL}.
#' @param beta Numeric scalar in (-1, 0); IMQ exponent. Defaults to
#'   -0.5 (Gorham-Mackey canonical choice). Ignored when
#'   \code{kernel = "rbf"}.
#' @param ess_weighted Logical scalar; if \code{TRUE}, thins both EB
#'   and FB draws to \eqn{\min(\widehat{\mathrm{ESS}}_{EB},
#'   \widehat{\mathrm{ESS}}_{FB})} common rows. Defaults to
#'   \code{FALSE}.
#' @param seed Integer scalar; RNG seed for ess_weighted thinning (no
#'   effect if \code{ess_weighted = FALSE}). Defaults to \code{NULL}.
#' @param ... Reserved for future arguments; currently unused.
#'
#' @return An object of class \code{gdpar_ksd_joint} with components:
#'   \describe{
#'     \item{\code{ksd_value}}{Numeric scalar; the KSD V-statistic
#'       (square root, clamped to \eqn{\ge 0}).}
#'     \item{\code{ksd_squared}}{Numeric scalar; the raw V-statistic
#'       before clamping (may be slightly negative under numerical
#'       noise; a negative value of small magnitude is consistent
#'       with a true KSD of zero).}
#'     \item{\code{kernel}}{Character scalar; \code{"imq"} or
#'       \code{"rbf"}.}
#'     \item{\code{bandwidth}}{Character scalar; \code{"median"} or
#'       \code{"fixed"}.}
#'     \item{\code{bandwidth_value}}{Numeric scalar; the bandwidth
#'       used (median heuristic result or supplied fixed value).}
#'     \item{\code{beta}}{Numeric scalar; IMQ exponent (or \code{NA}
#'       for RBF).}
#'     \item{\code{n_eb_draws}, \code{n_fb_draws}}{Integer scalars;
#'       row counts after optional thinning.}
#'     \item{\code{n_dim}}{Integer scalar; dimension of the common
#'       \eqn{\xi} vector.}
#'     \item{\code{target_mu}, \code{target_Sigma}}{Empirical mean
#'       and covariance of the FB draws over the common \eqn{\xi}
#'       variables; the Gaussian target.}
#'     \item{\code{ess_weighted}, \code{thinned_to}}{Logical and
#'       integer; thinning configuration.}
#'     \item{\code{vars}}{Character vector; common parameter names
#'       between EB and FB used for the computation.}
#'     \item{\code{call}}{The matched call.}
#'   }
#'   See \code{\link{print.gdpar_ksd_joint}} and
#'   \code{\link{summary.gdpar_ksd_joint}}.
#'
#' @seealso \code{\link{gdpar_compare_eb_fb}} (marginal TV
#'   comparator); \code{\link{gdpar_eb}}; \code{\link{gdpar}}.
#'
#' @references
#' Gorham, J., Mackey, L. (2017). Measuring Sample Quality with
#'   Kernels. JMLR 18(196):1-72.
#'
#' Liu, Q., Lee, J., Jordan, M. (2016). A Kernelized Stein
#'   Discrepancy for Goodness-of-Fit Tests. ICML.
#'
#' @examples
#' \dontrun{
#' # Fit EB and FB on the same data set
#' fit_fb <- gdpar(formula = y ~ x, data = my_data,
#'                 prior = gdpar_prior(), refresh = 0)
#' fit_eb <- gdpar_eb(formula = y ~ x, data = my_data,
#'                    prior = gdpar_prior(), refresh = 0)
#'
#' # Joint KSD on the xi posterior under empirical Gaussian target
#' ksd <- gdpar_ksd_joint(fit_eb, fit_fb)
#' print(ksd)
#'
#' # ESS-weighted variant
#' ksd_ess <- gdpar_ksd_joint(fit_eb, fit_fb,
#'                            ess_weighted = TRUE, seed = 1L)
#' summary(ksd_ess)
#' }
#'
#' @export
gdpar_ksd_joint <- function(eb_fit, fb_fit,
                            kernel = c("imq", "rbf"),
                            bandwidth = c("median", "fixed"),
                            bandwidth_value = NULL,
                            beta = -0.5,
                            ess_weighted = FALSE,
                            seed = NULL,
                            ...) {
  call <- match.call()
  kernel <- match.arg(kernel)
  bandwidth <- match.arg(bandwidth)
  .gdpar_ksd_validate_inputs(eb_fit, fb_fit, kernel, bandwidth,
                             bandwidth_value, beta, ess_weighted,
                             seed)

  fit_obj_eb <- eb_fit$conditional_fit %||% eb_fit$fit
  fit_obj_fb <- fb_fit$conditional_fit %||% fb_fit$fit
  draws_eb <- .gdpar_eb_fb_extract_xi_draws(fit_obj_eb)
  draws_fb <- .gdpar_eb_fb_extract_xi_draws(fit_obj_fb)
  if (is.null(draws_eb) || is.null(draws_fb)) {
    gdpar_abort(
      "Unable to extract xi draws from one of the fits; both must contain a non-empty draws array.",
      class = "gdpar_input_error"
    )
  }

  vars_common <- intersect(colnames(draws_eb), colnames(draws_fb))
  if (length(vars_common) == 0L) {
    gdpar_abort(
      "EB and FB fits share no common xi-variables; cannot compute joint KSD.",
      class = "gdpar_input_error"
    )
  }
  eb_mat <- as.matrix(draws_eb[, vars_common, drop = FALSE])
  fb_mat <- as.matrix(draws_fb[, vars_common, drop = FALSE])

  thinned_to <- NA_integer_
  if (isTRUE(ess_weighted)) {
    ess_eb <- .gdpar_ksd_safe_min_ess(draws_eb[, vars_common, drop = FALSE])
    ess_fb <- .gdpar_ksd_safe_min_ess(draws_fb[, vars_common, drop = FALSE])
    n_eff <- as.integer(floor(min(ess_eb, ess_fb,
                                   nrow(eb_mat), nrow(fb_mat))))
    if (n_eff < 2L) {
      gdpar_abort(
        sprintf(
          "ESS-weighted thinning yields n_eff = %d < 2; refusing to compute KSD on degenerate samples.",
          n_eff
        ),
        class = "gdpar_input_error"
      )
    }
    if (!is.null(seed)) set.seed(as.integer(seed))
    eb_idx <- sample(seq_len(nrow(eb_mat)), n_eff)
    fb_idx <- sample(seq_len(nrow(fb_mat)), n_eff)
    eb_mat <- eb_mat[eb_idx, , drop = FALSE]
    fb_mat <- fb_mat[fb_idx, , drop = FALSE]
    thinned_to <- n_eff
  }

  target_mu <- colMeans(fb_mat)
  target_Sigma <- stats::cov(fb_mat)
  Lambda <- .gdpar_ksd_safe_inverse(target_Sigma)

  bw_used <- if (bandwidth == "median") {
    .gdpar_ksd_median_squared_distance(fb_mat)
  } else if (is.null(bandwidth_value)) {
    1.0
  } else {
    as.numeric(bandwidth_value)
  }

  kp_sum <- .gdpar_ksd_stein_kernel_sum(eb_mat, target_mu, Lambda,
                                        kernel, bw_used, beta)
  n_eb <- nrow(eb_mat)
  ksd_sq <- kp_sum / (n_eb * n_eb)
  ksd_val <- sqrt(max(0.0, ksd_sq))

  out <- list(
    ksd_value = ksd_val,
    ksd_squared = ksd_sq,
    kernel = kernel,
    bandwidth = bandwidth,
    bandwidth_value = bw_used,
    beta = if (kernel == "imq") beta else NA_real_,
    n_eb_draws = n_eb,
    n_fb_draws = nrow(fb_mat),
    n_dim = length(vars_common),
    target_mu = target_mu,
    target_Sigma = target_Sigma,
    ess_weighted = isTRUE(ess_weighted),
    thinned_to = thinned_to,
    vars = vars_common,
    call = call
  )
  class(out) <- c("gdpar_ksd_joint", "list")
  out
}

#' @keywords internal
#' @noRd
.gdpar_ksd_validate_inputs <- function(eb_fit, fb_fit, kernel,
                                       bandwidth, bandwidth_value,
                                       beta, ess_weighted, seed) {
  if (!inherits(eb_fit, "gdpar_eb_fit")) {
    gdpar_abort(
      "Argument 'eb_fit' must be an object of class 'gdpar_eb_fit'.",
      class = "gdpar_input_error",
      data = list(received_class = class(eb_fit))
    )
  }
  if (!inherits(fb_fit, "gdpar_fit")) {
    gdpar_abort(
      "Argument 'fb_fit' must be an object of class 'gdpar_fit'.",
      class = "gdpar_input_error",
      data = list(received_class = class(fb_fit))
    )
  }
  if (!is.numeric(beta) || length(beta) != 1L || is.na(beta) ||
      beta <= -1.0 || beta >= 0.0) {
    gdpar_abort(
      "Argument 'beta' must be a numeric scalar in the open interval (-1, 0).",
      class = "gdpar_input_error",
      data = list(received = beta)
    )
  }
  if (bandwidth == "fixed" && !is.null(bandwidth_value)) {
    if (!is.numeric(bandwidth_value) || length(bandwidth_value) != 1L ||
        is.na(bandwidth_value) || bandwidth_value <= 0.0) {
      gdpar_abort(
        "Argument 'bandwidth_value' must be a positive numeric scalar when supplied.",
        class = "gdpar_input_error",
        data = list(received = bandwidth_value)
      )
    }
  }
  if (!is.logical(ess_weighted) || length(ess_weighted) != 1L ||
      is.na(ess_weighted)) {
    gdpar_abort(
      "Argument 'ess_weighted' must be a logical scalar (TRUE or FALSE).",
      class = "gdpar_input_error",
      data = list(received = ess_weighted)
    )
  }
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L ||
                          is.na(seed))) {
    gdpar_abort(
      "Argument 'seed' must be a numeric scalar or NULL.",
      class = "gdpar_input_error",
      data = list(received = seed)
    )
  }
  invisible(NULL)
}

#' @keywords internal
#' @noRd
.gdpar_ksd_median_squared_distance <- function(X) {
  n <- nrow(X)
  if (n < 2L) return(1.0)
  ss <- rowSums(X * X)
  D2 <- outer(ss, ss, "+") - 2.0 * tcrossprod(X)
  D2[D2 < 0.0] <- 0.0
  upper <- D2[upper.tri(D2)]
  med <- stats::median(upper)
  if (!is.finite(med) || med <= 0.0) 1.0 else med
}

#' @keywords internal
#' @noRd
.gdpar_ksd_safe_inverse <- function(Sigma) {
  d <- ncol(Sigma)
  diag_eps <- max(1e-8, mean(diag(Sigma)) * 1e-10)
  inverted <- tryCatch(chol2inv(chol(Sigma)), error = function(e) NULL)
  if (is.null(inverted)) {
    Sigma_ridged <- Sigma + diag(diag_eps, d, d)
    inverted <- tryCatch(chol2inv(chol(Sigma_ridged)),
                         error = function(e) NULL)
  }
  if (is.null(inverted)) {
    gdpar_abort(
      "Empirical Sigma is numerically singular even after ridging; cannot compute Gaussian-target score.",
      class = "gdpar_internal_error"
    )
  }
  inverted
}

#' @keywords internal
#' @noRd
.gdpar_ksd_safe_min_ess <- function(draws_mat) {
  # B9.4 G.iv/H.iv post-fix: use the canonical posterior::summarise_draws
  # API to obtain per-variable ess_basic; calling posterior::ess_basic()
  # directly on a draws_matrix is method-dispatched ambiguously across
  # posterior versions (returns either a numeric scalar or a tibble).
  ess_res <- tryCatch(
    posterior::summarise_draws(draws_mat,
                               ess_basic = posterior::ess_basic),
    error = function(e) NULL
  )
  if (is.null(ess_res)) return(nrow(draws_mat))
  vals <- if (is.data.frame(ess_res) && "ess_basic" %in% names(ess_res)) {
    ess_res$ess_basic
  } else if (is.numeric(ess_res)) {
    as.numeric(ess_res)
  } else {
    return(nrow(draws_mat))
  }
  vals <- vals[is.finite(vals) & vals > 0.0]
  if (length(vals) == 0L) return(nrow(draws_mat))
  min(vals)
}

#' Internal: Stein-kernel V-statistic numerator under Gaussian target
#'
#' Computes \eqn{\sum_{i,j} k_p(x_i, x_j)} where \eqn{k_p} is the
#' Stein kernel with target score \eqn{s(x) = -\Lambda (x - \mu)}
#' and base kernel either IMQ
#' \eqn{k(x, y) = (h + \|x - y\|^2)^\beta} or RBF
#' \eqn{k(x, y) = \exp(-\|x - y\|^2 / (2 h))}. Vectorized in R; cost
#' \eqn{O(n^2 d^2)} dominated by the score-score outer products.
#'
#' @keywords internal
#' @noRd
.gdpar_ksd_stein_kernel_sum <- function(X, mu, Lambda, kernel, h,
                                        beta) {
  n <- nrow(X)
  d <- ncol(X)
  Xc <- X - matrix(mu, n, d, byrow = TRUE)
  S <- -(Xc %*% Lambda)
  ss <- rowSums(X * X)
  D2 <- outer(ss, ss, "+") - 2.0 * tcrossprod(X)
  D2[D2 < 0.0] <- 0.0
  SX <- rowSums(S * X)
  S_X_t <- S %*% t(X)
  X_S_t <- X %*% t(S)
  SS_t <- tcrossprod(S)

  if (kernel == "imq") {
    u <- h + D2
    K  <- u ^ beta
    K1 <- u ^ (beta - 1.0)
    K2 <- u ^ (beta - 2.0)
    T1 <- SS_t * K
    T2 <- -2.0 * beta * K1 * (matrix(SX, n, n) - S_X_t)
    T3 <-  2.0 * beta * K1 * (X_S_t - matrix(SX, n, n, byrow = TRUE))
    T4 <- -4.0 * beta * (beta - 1.0) * K2 * D2 -
           2.0 * beta * d * K1
  } else {
    K  <- exp(-D2 / (2.0 * h))
    coef1 <- -1.0 / h
    T1 <- SS_t * K
    T2 <- coef1 * K * (-(matrix(SX, n, n) - S_X_t))
    T3 <- coef1 * K * ( (X_S_t - matrix(SX, n, n, byrow = TRUE)))
    T4 <- K * (d / h - D2 / (h * h))
  }
  sum(T1 + T2 + T3 + T4)
}

#' Print method for gdpar_ksd_joint
#'
#' @param x A \code{gdpar_ksd_joint} object.
#' @param digits Integer scalar; significant digits for the KSD value.
#'   Defaults to 4.
#' @param ... Reserved.
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.gdpar_ksd_joint <- function(x, digits = 4L, ...) {
  cat("<gdpar_ksd_joint>\n")
  cat(sprintf("  Joint KSD (sqrt V-statistic) : %s\n",
              format(x$ksd_value, digits = digits)))
  cat(sprintf("  Kernel                       : %s",
              x$kernel))
  if (x$kernel == "imq") {
    cat(sprintf(" (beta = %s)\n", format(x$beta, digits = digits)))
  } else {
    cat("\n")
  }
  cat(sprintf("  Bandwidth (squared units)    : %s (%s heuristic)\n",
              format(x$bandwidth_value, digits = digits),
              x$bandwidth))
  cat(sprintf("  Common xi dimension          : %d\n", x$n_dim))
  cat(sprintf("  EB draws / FB draws used     : %d / %d\n",
              x$n_eb_draws, x$n_fb_draws))
  if (isTRUE(x$ess_weighted)) {
    cat(sprintf("  ESS-weighted thinning        : yes (n_eff = %d)\n",
                x$thinned_to))
  } else {
    cat("  ESS-weighted thinning        : no\n")
  }
  invisible(x)
}

#' Summary method for gdpar_ksd_joint
#'
#' @param object A \code{gdpar_ksd_joint} object.
#' @param ... Reserved.
#'
#' @return An object of class \code{summary.gdpar_ksd_joint} with
#'   components \code{ksd_value}, \code{ksd_squared}, \code{kernel},
#'   \code{bandwidth_value}, \code{n_dim}, \code{vars}, and
#'   \code{interpretation} (character scalar). Use \code{print()} on
#'   the returned object for a formatted display.
#'
#' @export
summary.gdpar_ksd_joint <- function(object, ...) {
  out <- list(
    ksd_value = object$ksd_value,
    ksd_squared = object$ksd_squared,
    kernel = object$kernel,
    bandwidth_value = object$bandwidth_value,
    n_dim = object$n_dim,
    vars = object$vars,
    interpretation = paste0(
      "A value close to zero indicates EB matches the empirical ",
      "Gaussian fit of FB on the joint xi posterior; positive ",
      "values indicate joint distributional deviation. ",
      "Cross-check with gdpar_compare_eb_fb()'s marginal TV table."
    )
  )
  class(out) <- c("summary.gdpar_ksd_joint", "list")
  out
}

#' @export
print.summary.gdpar_ksd_joint <- function(x, digits = 4L, ...) {
  cat("<summary.gdpar_ksd_joint>\n")
  cat(sprintf("  KSD value                : %s (V-stat sqrt)\n",
              format(x$ksd_value, digits = digits)))
  cat(sprintf("  KSD squared              : %s\n",
              format(x$ksd_squared, digits = digits)))
  cat(sprintf("  Kernel                   : %s\n", x$kernel))
  cat(sprintf("  Bandwidth (sq. units)    : %s\n",
              format(x$bandwidth_value, digits = digits)))
  cat(sprintf("  Common xi dimension      : %d\n", x$n_dim))
  cat(sprintf("  Common xi variables      : %s\n",
              paste(x$vars, collapse = ", ")))
  cat("\n  Interpretation:\n  ")
  cat(strwrap(x$interpretation, width = 70L, indent = 0L, exdent = 2L),
      sep = "\n  ")
  cat("\n")
  invisible(x)
}

# Block 9, Sub-block 9.1, Session B9.8 (decision D75): dependence-robust
# inference (Axis 2). gdpar does NOT model residual correlated noise or
# spatial random effects (that is Axis 1 / Block 10, deferred and
# evidence-gated by the external validation 9.2). What this file adds is
# the *inferential* counterpart: (i) a residual dependence diagnostic that
# turns the iid-violation risk from invisible-theoretical into a measured
# quantity, and (ii) a temporal block-bootstrap-by-refit that delivers
# dependence-robust standard errors and intervals in the spirit of the
# working-independence + robust-variance estimator of Liang & Zeger (1986).
#
# Honest scope, stated without ornament: the block bootstrap yields robust
# *variance / intervals*, NOT better point estimates. The estimates remain
# those of a working-independence model: consistent when the mean structure
# is correctly specified (the dynamic-parameter + B-spline flexibility of
# the AMM helps achieve this), but neither efficient nor exploiting the
# correlation.
#
# Coverage evolution (Axis 2): the diagnostic and block-bootstrap machinery
# below started on the scalar Empirical-Bayes path (D75), then grew the
# spatial analogue (Moran's I + spatial block bootstrap; D100, B9.37), the
# data-driven Politis-White / spatial block length (D101, B9.38) and finally
# the full-Bayes analogue (this file's class-dispatched touchpoints; D102,
# B9.39), which closes the EB/FB asymmetry by accepting a scalar gdpar_fit
# wherever a scalar gdpar_eb_fit was accepted. The shared engine, the
# temporal/spatial diagnostics and the block-length selectors are
# path-agnostic; the FB extension dispatches only THREE EB-specific
# touchpoints by object class (the point-estimate vector, the model-SE
# vector and the Dunn-Smyth residuals). K > 1 / p > 1 (both EB and FB)
# remain deferred.

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' Internal: observed scalar outcome of a scalar Empirical-Bayes fit.
#'
#' Reads the outcome that was fed to Stan from \code{object$stan_data}
#' (\code{y_real} for continuous families, \code{y_int} for count /
#' Bernoulli families). Aborts for non-scalar outcomes, which are deferred.
#'
#' @keywords internal
#' @noRd
.gdpar_eb_scalar_y_obs <- function(object) {
  sd <- object$stan_data
  y_raw <- if (!is.null(sd$y_real)) sd$y_real else sd$y_int
  if (is.null(y_raw)) {
    gdpar_abort(
      "Internal error: gdpar_eb_fit lacks y_real / y_int in stan_data.",
      class = "gdpar_internal_error"
    )
  }
  if (is.matrix(y_raw) && ncol(y_raw) > 1L) {
    gdpar_abort(
      paste(
        "Dependence-robust inference is canonized for the scalar",
        "Empirical-Bayes path (K = 1, p = 1) in this sub-block; multivariate",
        "(p > 1) and multi-slot (K > 1) outcomes are deferred."
      ),
      class = "gdpar_unsupported_feature_error"
    )
  }
  as.numeric(y_raw)
}

#' Internal: assert that a fit is a scalar Empirical-Bayes fit.
#'
#' @keywords internal
#' @noRd
.gdpar_assert_scalar_eb <- function(object, arg_name = "object") {
  assert_inherits(object, "gdpar_eb_fit", arg_name)
  if (!is.null(object$family$families)) {
    gdpar_abort(
      paste(
        "Dependence-robust inference is canonized for the scalar",
        "Empirical-Bayes path in this sub-block; heterogeneous K > 1",
        "families are deferred."
      ),
      class = "gdpar_unsupported_feature_error"
    )
  }
  if (is.null(object$conditional_fit)) {
    gdpar_abort(
      "The gdpar_eb_fit object lacks its conditional HMC fit.",
      class = "gdpar_internal_error"
    )
  }
  invisible(object)
}

#' Internal: assert a scalar fit on either the EB or the full-Bayes path.
#'
#' The Axis 2 gate (decision D102): dependence-robust inference accepts a
#' scalar Empirical-Bayes fit (\code{gdpar_eb_fit}) \emph{or} a scalar
#' full-Bayes fit (\code{gdpar_fit}, \code{K = 1}, \code{p = 1}). For the EB
#' branch the existing \code{.gdpar_assert_scalar_eb()} is reused verbatim so
#' that the EB path stays byte-identical; the FB branch checks the path class
#' (\code{.gdpar_fit_path_class}) and the presence of the HMC fit. Heterogeneous
#' / multivariate / K-individual fits are deferred on both paths.
#'
#' @keywords internal
#' @noRd
.gdpar_assert_scalar_dep <- function(object, arg_name = "object") {
  if (inherits(object, "gdpar_eb_fit")) {
    return(.gdpar_assert_scalar_eb(object, arg_name))
  }
  if (inherits(object, "gdpar_fit")) {
    path <- .gdpar_fit_path_class(object)
    if (!identical(path, "scalar")) {
      gdpar_abort(
        paste(
          "Dependence-robust inference is canonized for the scalar path",
          "(K = 1, p = 1); multivariate (p > 1) and K-individual (K > 1)",
          "full-Bayes fits are deferred."
        ),
        class = "gdpar_unsupported_feature_error"
      )
    }
    if (is.null(object$fit)) {
      gdpar_abort(
        "The gdpar_fit object lacks its HMC fit.",
        class = "gdpar_internal_error"
      )
    }
    return(invisible(object))
  }
  gdpar_abort(
    sprintf("Argument '%s' must be a gdpar_eb_fit or a gdpar_fit.", arg_name),
    class = "gdpar_input_error"
  )
}

#' Internal: flatten the EB coefficient list into a named numeric vector.
#'
#' Concatenates, in a fixed order, the EB point estimates of
#' \code{theta_ref}, \code{a_coef}, \code{b_coef} and \code{W_raw} as
#' returned by \code{coef.gdpar_eb_fit()}. The names and order are stable
#' across refits of the same model specification, which the block
#' bootstrap relies on for column alignment.
#'
#' @keywords internal
#' @noRd
.gdpar_eb_estimate_vector <- function(fit) {
  co <- stats::coef(fit)
  parts <- list()
  tr <- co$theta_ref$estimate
  if (!is.null(tr)) {
    nm <- names(tr)
    if (is.null(nm)) {
      nm <- if (length(tr) == 1L) "theta_ref" else
        paste0("theta_ref[", seq_along(tr), "]")
    }
    parts[["theta_ref"]] <- stats::setNames(as.numeric(tr), nm)
  }
  for (comp in c("a", "b", "W")) {
    est <- co[[comp]]$estimate
    if (!is.null(est)) {
      nm <- names(est)
      if (is.null(nm)) nm <- paste0(comp, "[", seq_along(est), "]")
      parts[[comp]] <- stats::setNames(as.numeric(est), nm)
    }
  }
  if (length(parts) == 0L) {
    gdpar_abort(
      "Could not extract any EB point estimates from the fit.",
      class = "gdpar_internal_error"
    )
  }
  do.call(c, unname(parts))
}

#' Internal: model-based (Laplace / posterior) standard errors as a vector.
#'
#' Mirrors \code{.gdpar_eb_estimate_vector()} but reads the \code{se}
#' field of each component, so the two vectors are name-aligned.
#'
#' @keywords internal
#' @noRd
.gdpar_eb_model_se_vector <- function(fit) {
  co <- stats::coef(fit)
  parts <- list()
  tr_se <- co$theta_ref$se
  tr_est <- co$theta_ref$estimate
  if (!is.null(tr_est)) {
    nm <- names(tr_est)
    if (is.null(nm)) {
      nm <- if (length(tr_est) == 1L) "theta_ref" else
        paste0("theta_ref[", seq_along(tr_est), "]")
    }
    val <- if (is.null(tr_se)) rep(NA_real_, length(tr_est)) else
      as.numeric(tr_se)
    parts[["theta_ref"]] <- stats::setNames(val, nm)
  }
  for (comp in c("a", "b", "W")) {
    est <- co[[comp]]$estimate
    if (!is.null(est)) {
      nm <- names(est)
      if (is.null(nm)) nm <- paste0(comp, "[", seq_along(est), "]")
      se <- co[[comp]]$se
      val <- if (is.null(se)) rep(NA_real_, length(est)) else as.numeric(se)
      parts[[comp]] <- stats::setNames(val, nm)
    }
  }
  do.call(c, unname(parts))
}

#' Internal: posterior draws matrix of the AMM coefficients of a full-Bayes fit.
#'
#' The full-Bayes counterpart of the EB coefficient extraction. Reads the
#' posterior draws of the AMM coefficients of a scalar \code{gdpar_fit} --
#' \code{theta_ref} (incl. its per-group entries under grouping), and, when the
#' corresponding AMM component is active, \code{a_coef}, \code{b_coef} and the
#' raw basis weights \code{W_raw} -- as a single \code{S x P} matrix whose
#' columns carry the Stan variable names (\code{theta_ref[1]}, \code{a_coef[1]},
#' ...). Two deliberate parity choices with \code{.gdpar_eb_estimate_vector()}
#' (decision D102): (i) the raw \code{W_raw} draws are used, NOT the
#' \code{sigma_W}-scaled effective weights of \code{coef.gdpar_fit()}, exactly as
#' the EB extractor reads the raw conditional \code{W_raw} draws -- the
#' dependence-robust SE only needs a coefficient parametrization that is
#' consistent between the original fit and its refits, and the raw scale matches
#' the EB contract; (ii) the \code{theta_ref} hyperparameters
#' (\code{mu_theta_ref}, \code{sigma_theta_ref}) are excluded, again for EB/FB
#' parity and the "AMM coefficient" framing. Column names are stable across
#' refits of the same specification, which the block bootstrap relies on for
#' alignment.
#'
#' @keywords internal
#' @noRd
.gdpar_fb_coef_draws_matrix <- function(object) {
  require_suggested("posterior", "extract posterior draws")
  draws <- object$fit$draws()
  vars <- "theta_ref"
  if (!is.null(object$amm$a)) vars <- c(vars, "a_coef")
  if (!is.null(object$amm$b)) vars <- c(vars, "b_coef")
  if (!is.null(object$amm$W)) vars <- c(vars, "W_raw")
  mat <- unclass(posterior::as_draws_matrix(
    posterior::subset_draws(draws, variable = vars)
  ))
  if (is.null(mat) || ncol(mat) == 0L) {
    gdpar_abort(
      "Could not extract any AMM coefficient draws from the full-Bayes fit.",
      class = "gdpar_internal_error"
    )
  }
  mat
}

#' Internal: full-Bayes point-estimate vector (posterior mean per coefficient).
#'
#' @keywords internal
#' @noRd
.gdpar_fb_estimate_vector <- function(object) {
  mat <- .gdpar_fb_coef_draws_matrix(object)
  stats::setNames(colMeans(mat), colnames(mat))
}

#' Internal: full-Bayes model-SE vector (posterior SD per coefficient).
#'
#' @keywords internal
#' @noRd
.gdpar_fb_model_se_vector <- function(object) {
  mat <- .gdpar_fb_coef_draws_matrix(object)
  stats::setNames(apply(mat, 2L, stats::sd), colnames(mat))
}

#' Internal: class-dispatched point-estimate vector (EB or full-Bayes).
#'
#' The shared block-bootstrap engine's first touchpoint. For a
#' \code{gdpar_eb_fit} it calls \code{.gdpar_eb_estimate_vector()} verbatim
#' (so the EB path is byte-identical and the EB regression gate is preserved);
#' for a \code{gdpar_fit} it returns the full-Bayes posterior means.
#'
#' @keywords internal
#' @noRd
.gdpar_dep_estimate_vector <- function(object) {
  if (inherits(object, "gdpar_eb_fit")) {
    return(.gdpar_eb_estimate_vector(object))
  }
  .gdpar_fb_estimate_vector(object)
}

#' Internal: class-dispatched model-SE vector (EB or full-Bayes).
#'
#' The engine's second touchpoint. EB: Laplace / conditional posterior SD
#' (verbatim). Full-Bayes: posterior SD of each coefficient. In both cases the
#' "model SE" is a within-model (posterior / Laplace) standard deviation, so
#' \code{se_ratio = robust_se / model_se} is a like-for-like (SD-vs-SD) ratio.
#'
#' @keywords internal
#' @noRd
.gdpar_dep_model_se_vector <- function(object) {
  if (inherits(object, "gdpar_eb_fit")) {
    return(.gdpar_eb_model_se_vector(object))
  }
  .gdpar_fb_model_se_vector(object)
}

#' Internal: rate-optimal default block length for the moving block bootstrap.
#'
#' Returns \code{max(1, round(n^(1/3)))}, the optimal growth rate for the
#' moving block bootstrap variance estimator (Kuensch 1989; Hall, Horowitz
#' & Jing 1995). The data-driven constant of Politis & White (2004) is a
#' deferred refinement; this default fixes only the rate.
#'
#' @keywords internal
#' @noRd
.gdpar_default_block_length <- function(n) {
  max(1L, as.integer(round(n^(1 / 3))))
}

#' Internal: is a block-size argument the literal string \code{"auto"}?
#'
#' Distinguishes the opt-in data-driven path (\code{block_length = "auto"} /
#' \code{block_size = "auto"}) from the rate default (\code{NULL}) and a fixed
#' integer. Shared by the temporal and spatial robust estimators.
#'
#' @keywords internal
#' @noRd
.gdpar_is_auto <- function(x) {
  is.character(x) && length(x) == 1L && !is.na(x) && identical(x, "auto")
}

#' Internal: flat-top lag window of Politis (2003) / Politis & White (2004).
#'
#' \eqn{\lambda(s) = 1} for \eqn{|s| \le 1/2}, \eqn{2(1 - |s|)} for
#' \eqn{1/2 < |s| \le 1}, and 0 otherwise. Vectorised over \code{s}.
#'
#' @keywords internal
#' @noRd
.gdpar_flat_top_window <- function(s) {
  a <- abs(s)
  ifelse(a <= 0.5, 1, ifelse(a <= 1, 2 * (1 - a), 0))
}

#' Internal: adaptive bandwidth \eqn{\hat m} of Politis (2003).
#'
#' Given the lag-1.. autocorrelations \code{rho} and the critical value
#' \code{crit}, returns the smallest \eqn{m} such that
#' \eqn{|\hat\rho(m + j)| < \mathrm{crit}} for \code{Kn} consecutive lags
#' \eqn{j = 1, \dots, K_N} (the "first run of \code{Kn} negligible lags"
#' rule). If no such run exists, returns the largest significant lag (so the
#' window keeps the dependence it can see); if every lag is negligible,
#' returns 1. Factored out for direct unit testing.
#'
#' @keywords internal
#' @noRd
.gdpar_pw_mhat <- function(rho, Kn, crit) {
  insig <- abs(rho) < crit
  L <- length(rho)
  if (L >= Kn) {
    for (j in seq_len(L - Kn + 1L)) {
      if (all(insig[j:(j + Kn - 1L)])) return(as.integer(j))
    }
  }
  sig <- which(!insig)
  if (length(sig) > 0L) as.integer(sig[length(sig)]) else 1L
}

#' Internal: data-driven block length via the Politis & White (2004) selector.
#'
#' Implements the automatic block-length rule of Politis & White (2004) with
#' the Patton, Politis & White (2009) correction, hand-rolled in base R (no
#' \pkg{np} dependency) and matching \code{np::b.star}. It operates on the
#' \emph{residuals} of the fitted working-independence model (cheap, no
#' refit): the adaptive flat-top spectral estimates \eqn{\hat g} and
#' \eqn{\widehat{spec}} of the lag-window of the autocovariances give
#' \eqn{b_{opt} = (2\hat g^2 / D)^{1/3} n^{1/3}}, capped to
#' \eqn{[1, \lceil\min(3\sqrt n, n/3)\rceil]}. gdpar's moving / circular block
#' bootstrap is \emph{overlapping}, so the overlapping (circular-block)
#' constant \eqn{D = (4/3)\,\widehat{spec}^2} is used (Lahiri 2003: the moving
#' and circular block bootstrap share the leading variance constant), not the
#' stationary-bootstrap \eqn{2\,\widehat{spec}^2}.
#'
#' Degenerate inputs fall back to the \code{n^(1/3)} rate with a stated reason:
#' \code{n} too small for the lag scan, non-positive residual variance, or a
#' non-positive flat-top spectral estimate. A near-white residual gives
#' \eqn{\hat g \approx 0 \Rightarrow b_{opt} \to 0}, floored at 1 (unit
#' blocks), which is the honest data-driven answer, not a fallback.
#'
#' @param resid Numeric residual vector, already in the bootstrap (index)
#'   order.
#' @param c_thresh Critical-value constant for the bandwidth test; default
#'   \code{qnorm(0.975)} (matches \code{np::b.star}).
#' @return A list with \code{block_length} (integer), \code{method}
#'   (\code{"auto"} or \code{"rate"}), and a human-readable \code{reason}.
#' @keywords internal
#' @noRd
.gdpar_politis_white_block_length <- function(resid, c_thresh = stats::qnorm(0.975)) {
  n <- length(resid)
  rate <- .gdpar_default_block_length(n)
  fb <- function(reason) list(block_length = rate, method = "rate", reason = reason)
  if (n < 8L) {
    return(fb(sprintf("n = %d too small for Politis-White; used n^(1/3) rate.", n)))
  }
  if (!is.finite(stats::var(resid)) || stats::var(resid) <= 0) {
    return(fb("Residuals have non-positive variance; used n^(1/3) rate."))
  }
  Kn <- max(5L, as.integer(ceiling(log10(n))))
  Mmax <- min(as.integer(ceiling(sqrt(n)) + Kn), n - 1L)
  if (Mmax < Kn + 1L) {
    return(fb("Series too short for the Politis-White lag scan; used rate."))
  }
  rho <- stats::acf(resid, lag.max = Mmax, type = "correlation",
                    plot = FALSE, demean = TRUE)$acf[-1L]
  crit <- c_thresh * sqrt(log10(n) / n)
  mhat <- .gdpar_pw_mhat(rho, Kn, crit)
  M <- min(2L * mhat, Mmax)
  acov <- stats::acf(resid, lag.max = M, type = "covariance",
                     plot = FALSE, demean = TRUE)$acf
  R0 <- acov[1L]
  Rk <- acov[-1L]
  k <- seq_len(M)
  lam <- .gdpar_flat_top_window(k / M)
  spec <- R0 + 2 * sum(lam * Rk)        # sum_{|k|<=M} lambda(k/M) Rhat(k)
  ghat <- 2 * sum(lam * k * Rk)         # sum_{|k|<=M} lambda(k/M) |k| Rhat(k)
  if (!is.finite(spec) || spec <= 0 || !is.finite(ghat)) {
    return(fb("Non-positive flat-top spectral estimate; used n^(1/3) rate."))
  }
  D <- (4 / 3) * spec^2                  # overlapping (circular/moving) block
  b <- (2 * ghat^2 / D)^(1 / 3) * n^(1 / 3)
  Bmax <- ceiling(min(3 * sqrt(n), n / 3))
  b_int <- max(1L, min(as.integer(round(b)), as.integer(Bmax)))
  list(block_length = b_int, method = "auto",
       reason = sprintf("Politis-White (m_hat = %d, M = %d, b = %.2f, capped to [1, %d]).",
                        mhat, M, b, as.integer(Bmax)))
}

#' Internal: resampled row indices for a temporal block bootstrap.
#'
#' Forms a length-\code{n} index vector by concatenating
#' \code{ceiling(n / block_length)} contiguous blocks of size
#' \code{block_length} drawn with replacement, then truncating to
#' \code{n}. The \code{"moving"} scheme draws block starts in
#' \code{[1, n - block_length + 1]} (Kuensch 1989); the \code{"circular"}
#' scheme draws starts in \code{[1, n]} and wraps around, giving every
#' observation equal resampling weight (Politis & Romano 1992). Sibling of
#' the MCMC-draw resampler \code{block_bootstrap_indices()} with a single
#' chain.
#'
#' @keywords internal
#' @noRd
.gdpar_block_bootstrap_data_indices <- function(n, block_length,
                                                type = c("moving",
                                                         "circular")) {
  type <- match.arg(type)
  n <- as.integer(n)
  block_length <- as.integer(block_length)
  if (block_length < 1L || block_length > n) {
    gdpar_abort(
      sprintf("block_length must be in [1, %d].", n),
      class = "gdpar_input_error",
      data = list(block_length = block_length, n = n)
    )
  }
  n_blocks <- as.integer(ceiling(n / block_length))
  if (type == "moving") {
    n_starts <- n - block_length + 1L
    starts <- sample.int(n_starts, n_blocks, replace = TRUE)
    out <- integer(n_blocks * block_length)
    for (b in seq_len(n_blocks)) {
      pos <- (b - 1L) * block_length + 1L
      out[pos:(pos + block_length - 1L)] <-
        starts[b] + seq_len(block_length) - 1L
    }
  } else {
    starts <- sample.int(n, n_blocks, replace = TRUE)
    out <- integer(n_blocks * block_length)
    for (b in seq_len(n_blocks)) {
      pos <- (b - 1L) * block_length + 1L
      raw <- starts[b] + seq_len(block_length) - 1L
      out[pos:(pos + block_length - 1L)] <- ((raw - 1L) %% n) + 1L
    }
  }
  out[seq_len(n)]
}

#' Internal: Dunn-Smyth (or requested) residuals of a scalar fit (EB or FB).
#'
#' Shared by the temporal (\code{gdpar_dependence_diagnostic}) and spatial
#' (\code{gdpar_spatial_dependence_diagnostic}) dependence diagnostics: reads
#' the scalar outcome, the posterior predictive draws and dispatches to the
#' requested residual type. Extracted verbatim from the temporal diagnostic so
#' both paths share a single residual definition (decision D100). The third
#' EB/FB touchpoint (decision D102): for a full-Bayes \code{gdpar_fit} the
#' canonical \code{residuals.gdpar_fit()} method already computes the same
#' Bayesian Dunn-Smyth (randomized quantile) residuals from the posterior
#' predictive \code{y_pred} draws via the identical \code{.gdpar_residuals_dispatch()}
#' worker, so the FB branch simply delegates to it; the EB branch is unchanged.
#'
#' @keywords internal
#' @noRd
.gdpar_dependence_residuals <- function(object, residual_type, randomize_seed) {
  if (inherits(object, "gdpar_fit") && !inherits(object, "gdpar_eb_fit")) {
    return(as.numeric(stats::residuals(
      object, type = residual_type, randomize_seed = randomize_seed
    )))
  }
  y_obs <- .gdpar_eb_scalar_y_obs(object)
  pred <- stats::predict(object, type = "response")
  family_name <- object$family$name
  .gdpar_residuals_dispatch(
    y_obs, pred$draws, residual_type, family_name, randomize_seed
  )
}

# ---------------------------------------------------------------------------
# Residual dependence diagnostic
# ---------------------------------------------------------------------------

#' Residual dependence diagnostic for a scalar Empirical-Bayes fit
#'
#' Quantifies serial (temporal) dependence in the residuals of a fitted
#' scalar Path 1 Empirical-Bayes model. gdpar assumes conditional
#' independence; under temporal (or spatial) autocorrelation that
#' assumption is violated and the model-based (posterior / Laplace)
#' uncertainty is too narrow. This diagnostic makes the violation
#' *visible and measurable* before any remedy is applied, and is the
#' natural gate for \code{\link{gdpar_dependence_robust}}.
#'
#' @param object A scalar Path 1 fit (\code{K = 1}, \code{p = 1}): either a
#'   \code{gdpar_eb_fit} (Empirical Bayes, from \code{\link{gdpar_eb}}) or a
#'   \code{gdpar_fit} (full Bayes, from \code{\link{gdpar}}).
#' @param index Optional vector of length \code{n} giving the temporal (or
#'   any one-dimensional) ordering of the observations. Residuals are
#'   sorted by \code{order(index)} before the autocorrelation statistics
#'   are computed. When \code{NULL} (default) the natural row order of the
#'   training data is used.
#' @param residual_type One of \code{"quantile"} (default; randomized
#'   quantile / Dunn-Smyth residuals, the gdpar canonical choice),
#'   \code{"response"}, \code{"pearson"} or \code{"deviance"}. Forwarded to
#'   the internal residual dispatcher.
#' @param max_lag Integer scalar; the maximum lag for the Ljung-Box test.
#'   Defaults to \code{min(floor(10 * log10(n)), n - 1)}.
#' @param level Numeric scalar in (0, 1); the confidence level used to turn
#'   the test p-values into the verdict. Defaults to 0.95 (i.e. dependence
#'   is flagged when a p-value falls below \code{1 - level}).
#' @param randomize_seed Optional integer seed used by the randomized
#'   quantile residuals for discrete families; ignored otherwise. Pass a
#'   value for reproducibility.
#' @param ... Unused; present for signature stability.
#'
#' @return A list of class \code{gdpar_dependence_diagnostic} with
#'   components \code{residual_type}, \code{n}, \code{max_lag},
#'   \code{lag1_autocorr}, \code{lag1_p_value} (normal approximation
#'   \eqn{\sqrt{n}\,\hat{\rho}_1 \sim N(0, 1)}), \code{durbin_watson},
#'   \code{ljung_box_statistic}, \code{ljung_box_df},
#'   \code{ljung_box_p_value}, \code{level}, \code{index_supplied} and
#'   \code{verdict}. A \code{print} method provides a human-readable
#'   summary.
#'
#' @details
#' The Durbin-Watson statistic is reported descriptively as
#' \eqn{DW = \sum_{t=2}^{n}(r_t - r_{t-1})^2 / \sum_{t=1}^{n} r_t^2}
#' (\eqn{DW \approx 2(1 - \hat{\rho}_1)}); values near 2 indicate no
#' first-order autocorrelation. The Ljung-Box test (\code{stats::Box.test})
#' provides the omnibus p-value across lags and drives the verdict. The
#' Ljung-Box degrees of freedom are not reduced by the number of estimated
#' AMM coefficients (\code{fitdf = 0}); for residuals of a fitted model
#' this makes the test mildly optimistic, a caveat stated honestly rather
#' than masked.
#'
#' Spatial dependence (Moran's I and a spatial weight structure) is handled
#' by the sibling \code{\link{gdpar_spatial_dependence_diagnostic}}. Both the
#' scalar Empirical-Bayes (\code{gdpar_eb_fit}) and the scalar full-Bayes
#' (\code{gdpar_fit}) paths are supported (decision D102); the residuals are the
#' Bayesian Dunn-Smyth residuals of whichever fit is supplied. The K > 1 /
#' p > 1 paths remain deferred.
#'
#' @section Methodological note:
#' A flagged verdict says the model-based uncertainty is not trustworthy
#' under the detected dependence, not that the point estimates are wrong.
#' The companion remedy \code{\link{gdpar_dependence_robust}} re-estimates
#' the uncertainty by a temporal block bootstrap.
#'
#' @references
#' Ljung, G. M. & Box, G. E. P. (1978). On a measure of lack of fit in
#' time series models. \emph{Biometrika} 65(2), 297-303.
#'
#' Durbin, J. & Watson, G. S. (1950). Testing for serial correlation in
#' least squares regression. I. \emph{Biometrika} 37(3/4), 409-428.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("cmdstanr", quietly = TRUE) &&
#'     requireNamespace("posterior", quietly = TRUE)) {
#'   n <- 100
#'   x <- rnorm(n)
#'   y <- 1 + 0.5 * x + as.numeric(stats::arima.sim(list(ar = 0.6), n))
#'   df <- data.frame(x = x, y = y, t = seq_len(n))
#'   fit <- gdpar_eb(y ~ x, amm = amm_spec(a = ~ x), data = df,
#'                   chains = 2, iter_warmup = 100, iter_sampling = 100)
#'   gdpar_dependence_diagnostic(fit, index = df$t)
#' }
#' }
#'
#' @seealso \code{\link{gdpar_dependence_robust}}, \code{\link{gdpar_eb}}
#' @export
gdpar_dependence_diagnostic <- function(object, index = NULL,
                                        residual_type = c("quantile",
                                                          "response",
                                                          "pearson",
                                                          "deviance"),
                                        max_lag = NULL, level = 0.95,
                                        randomize_seed = NULL, ...) {
  .gdpar_assert_scalar_dep(object, "object")
  residual_type <- match.arg(residual_type)
  assert_numeric_scalar(level, "level", lower = 0, upper = 1)
  require_suggested("posterior", "extract posterior draws for residuals")

  resid <- .gdpar_dependence_residuals(object, residual_type, randomize_seed)
  n <- length(resid)

  if (!is.null(index)) {
    if (length(index) != n) {
      gdpar_abort(
        sprintf("Argument 'index' must have length n = %d.", n),
        class = "gdpar_input_error"
      )
    }
    resid <- resid[order(index)]
  }

  if (is.null(max_lag)) {
    max_lag <- min(as.integer(floor(10 * log10(n))), n - 1L)
  }
  max_lag <- as.integer(max_lag)
  if (max_lag < 1L || max_lag > n - 1L) {
    gdpar_abort(
      sprintf("Argument 'max_lag' must be in [1, %d].", n - 1L),
      class = "gdpar_input_error"
    )
  }

  r_centered <- resid - mean(resid)
  denom <- sum(r_centered^2)
  if (denom <= 0) {
    gdpar_abort(
      "Residuals have zero variance; cannot assess dependence.",
      class = "gdpar_diagnostic_error"
    )
  }
  rho1 <- sum(r_centered[-1L] * r_centered[-n]) / denom
  dw <- sum(diff(resid)^2) / sum(resid^2)
  lag1_p <- 2 * stats::pnorm(-abs(sqrt(n) * rho1))
  lb <- stats::Box.test(resid, lag = max_lag, type = "Ljung-Box", fitdf = 0L)

  alpha <- 1 - level
  verdict <- if (lb$p.value < alpha) {
    sprintf(
      paste(
        "Dependence detected (Ljung-Box p = %.3g < %.3g): the model-based",
        "uncertainty is not robust to it. Use gdpar_dependence_robust()."
      ),
      lb$p.value, alpha
    )
  } else {
    sprintf(
      "No evidence against independence at level %.2f (Ljung-Box p = %.3g).",
      level, lb$p.value
    )
  }

  obj <- list(
    residual_type       = residual_type,
    n                   = n,
    max_lag             = max_lag,
    lag1_autocorr       = unname(rho1),
    lag1_p_value        = unname(lag1_p),
    durbin_watson       = unname(dw),
    ljung_box_statistic = unname(lb$statistic),
    ljung_box_df        = unname(lb$parameter),
    ljung_box_p_value   = unname(lb$p.value),
    level               = level,
    index_supplied      = !is.null(index),
    verdict             = verdict
  )
  class(obj) <- c("gdpar_dependence_diagnostic", "list")
  obj
}

#' Print method for gdpar_dependence_diagnostic objects
#'
#' @param x A \code{gdpar_dependence_diagnostic} object.
#' @param digits Integer; significant digits for the printed statistics.
#' @param ... Unused; present for S3 generic compatibility.
#' @return Invisibly returns \code{x}.
#' @export
print.gdpar_dependence_diagnostic <- function(x, digits = 3L, ...) {
  cat("<gdpar_dependence_diagnostic>\n", sep = "")
  cat("  residual type   : ", x$residual_type,
      if (x$index_supplied) " (ordered by supplied index)" else
        " (natural row order)", "\n", sep = "")
  cat("  n               : ", x$n, "\n", sep = "")
  cat("  lag-1 autocorr  : ", format(x$lag1_autocorr, digits = digits),
      "  (p = ", format(x$lag1_p_value, digits = digits), ")\n", sep = "")
  cat("  Durbin-Watson   : ", format(x$durbin_watson, digits = digits),
      "  (2 = no autocorrelation)\n", sep = "")
  cat("  Ljung-Box       : Q(", x$ljung_box_df, ") = ",
      format(x$ljung_box_statistic, digits = digits),
      "  (p = ", format(x$ljung_box_p_value, digits = digits), ")\n", sep = "")
  cat("  Verdict         : ", x$verdict, "\n", sep = "")
  invisible(x)
}

# ---------------------------------------------------------------------------
# Shared block-bootstrap-by-refit engine
# ---------------------------------------------------------------------------

#' Internal: shared block-bootstrap-by-refit engine.
#'
#' Factors the refit loop, seed management, bootstrap-SE / percentile-interval
#' assembly and failure accounting shared by the temporal
#' (\code{gdpar_dependence_robust}) and spatial
#' (\code{gdpar_spatial_dependence_robust}) dependence-robust estimators. The
#' two public functions differ only in their \code{resample_fun} (a nullary
#' closure returning a length-\code{n} integer index vector into \code{data})
#' and in the descriptive metadata they attach to the result; everything
#' downstream of the resample is identical, so it lives here exactly once.
#'
#' Since decision D102 the point-estimate and model-SE extraction are
#' class-dispatched (\code{.gdpar_dep_estimate_vector} /
#' \code{.gdpar_dep_model_se_vector}) and the refit \code{call} is whatever
#' \code{object$call} encodes (\code{gdpar_eb()} or \code{gdpar()}), so the
#' engine serves the Empirical-Bayes and full-Bayes paths identically. For the
#' EB path the dispatch resolves to the original helpers and the refit is a
#' \code{gdpar_eb()} call, so the engine's behaviour is byte-identical to the
#' pre-D102 code -- the EB regression gate is preserved.
#'
#' The RNG-consumption order is contractually frozen: \code{set.seed(seed)} (when
#' supplied), then the \code{B} per-refit Stan seeds via
#' \code{sample.int(.Machine$integer.max, B)}, then one \code{resample_fun()}
#' call per bootstrap iteration. This reproduces, bit for bit, the pre-refactor
#' temporal path when \code{resample_fun} wraps
#' \code{.gdpar_block_bootstrap_data_indices()} (regression-gated, decision
#' D100).
#'
#' @param resample_fun A nullary function returning an integer vector of length
#'   \code{nrow(data)} indexing the rows of \code{data} to refit on.
#' @param verbose_msg Pre-formatted opt-in cost message, emitted once when
#'   \code{verbose} is \code{TRUE}.
#' @param caller_env The environment in which each refit \code{call} is
#'   evaluated (the public function's own \code{parent.frame()}), so model
#'   symbols resolve exactly as for a direct \code{gdpar_eb()} call.
#'
#' @return A list with \code{table} (one row per coefficient), \code{B_ok},
#'   \code{seed}, \code{warnings} and \code{refit_diagnostics} (aggregate
#'   per-refit convergence: \code{max_rhat}, \code{min_ess_bulk},
#'   \code{n_divergent_refits}, \code{n_high_rhat_refits}, \code{rhat_threshold});
#'   the public wrappers prepend their own descriptive fields and set the class.
#' @keywords internal
#' @noRd
.gdpar_dependence_robust_engine <- function(object, data, resample_fun, B, level,
                                            seed, iter_warmup, iter_sampling,
                                            chains, verbose, verbose_msg,
                                            caller_env, ...) {
  estimate <- .gdpar_dep_estimate_vector(object)
  model_se <- .gdpar_dep_model_se_vector(object)
  param_names <- names(estimate)
  model_se <- model_se[param_names]

  if (verbose && !is.null(verbose_msg)) {
    gdpar_inform(verbose_msg, class = "gdpar_optin_message")
  }

  if (!is.null(seed)) set.seed(as.integer(seed))
  refit_seeds <- sample.int(.Machine$integer.max, B)

  call <- object$call
  boot <- matrix(NA_real_, nrow = B, ncol = length(param_names),
                 dimnames = list(NULL, param_names))
  warnings_msg <- character(0L)
  # Per-refit convergence accounting (decision D102; cross-lineage-endorsed).
  # Informational only: refits are NEVER excluded or down-weighted -- excluding
  # under-converged refits is non-random (it removes the very data
  # configurations the bootstrap is meant to probe) and would bias the SE.
  # Both gdpar_eb_fit and gdpar_fit carry a $diagnostics slot with the same
  # rhat_max / ess_bulk_min / divergent_count fields, so this is path-agnostic
  # and does not change any number (the EB regression gate is preserved).
  refit_rhat <- rep(NA_real_, B)
  refit_ess  <- rep(NA_real_, B)
  refit_div  <- rep(NA_integer_, B)

  for (b in seq_len(B)) {
    idx <- resample_fun()
    sub <- data[idx, , drop = FALSE]
    sub_call <- call
    sub_call$data          <- quote(sub)
    sub_call$iter_warmup   <- iter_warmup
    sub_call$iter_sampling <- iter_sampling
    sub_call$chains        <- chains
    sub_call$verbose       <- FALSE
    sub_call$refresh       <- 0L
    sub_call$skip_id_check <- TRUE
    sub_call$seed          <- refit_seeds[b]
    env <- new.env(parent = caller_env)
    env$sub <- sub
    refit_failure_msg <- NULL
    vec <- tryCatch({
      fit_b <- eval(sub_call, envir = env)
      v <- .gdpar_dep_estimate_vector(fit_b)
      dg <- fit_b$diagnostics
      if (!is.null(dg)) {
        refit_rhat[b] <- dg$rhat_max %||% NA_real_
        refit_ess[b]  <- dg$ess_bulk_min %||% NA_real_
        refit_div[b]  <- as.integer(dg$divergent_count %||% NA_integer_)
      }
      v
    }, error = function(e) {
      refit_failure_msg <<- sprintf(
        "Refit %d failed: %s.", b, conditionMessage(e)
      )
      NULL
    })
    if (is.null(vec)) {
      if (!is.null(refit_failure_msg)) {
        warnings_msg <- c(warnings_msg, refit_failure_msg)
      }
      next
    }
    common <- intersect(param_names, names(vec))
    boot[b, common] <- vec[common]
  }

  ok <- stats::complete.cases(boot)
  B_ok <- sum(ok)
  if (B_ok < 2L) {
    gdpar_abort(
      "Fewer than two successful bootstrap refits; cannot estimate robust SE.",
      class = "gdpar_diagnostic_error",
      data = list(B = B, B_ok = B_ok)
    )
  }
  boot_ok <- boot[ok, , drop = FALSE]
  robust_se <- apply(boot_ok, 2L, stats::sd)
  alpha <- 1 - level
  ci <- apply(boot_ok, 2L, stats::quantile,
              probs = c(alpha / 2, 1 - alpha / 2), names = FALSE)

  table <- data.frame(
    parameter = param_names,
    estimate  = unname(estimate),
    model_se  = unname(model_se),
    robust_se = unname(robust_se),
    se_ratio  = unname(robust_se / model_se),
    ci_lower  = ci[1L, ],
    ci_upper  = ci[2L, ],
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  # Aggregate the per-refit convergence diagnostics over the successful refits.
  # The R-hat threshold (1.05) is the classical Gelman-Rubin "clearly
  # non-converged" line, declared as a heuristic; only a genuine R-hat breach
  # triggers a warning (a handful of divergent refits is common with short
  # bootstrap chains and the bootstrap SD is robust to them, so the divergence
  # count is reported but not, on its own, warned on).
  rhat_threshold <- 1.05
  refit_diagnostics <- list(
    max_rhat           = if (any(!is.na(refit_rhat)))
                           max(refit_rhat, na.rm = TRUE) else NA_real_,
    min_ess_bulk       = if (any(!is.na(refit_ess)))
                           min(refit_ess, na.rm = TRUE) else NA_real_,
    n_divergent_refits = as.integer(sum(refit_div > 0L, na.rm = TRUE)),
    n_high_rhat_refits = as.integer(sum(refit_rhat > rhat_threshold,
                                        na.rm = TRUE)),
    rhat_threshold     = rhat_threshold
  )
  if (isTRUE(refit_diagnostics$max_rhat > rhat_threshold)) {
    warnings_msg <- c(warnings_msg, sprintf(
      paste0("Refit convergence: max R-hat across refits = %s (%d refit(s) ",
             "with R-hat > %.2f, %d with divergences). The robust SE may be ",
             "affected by under-converged refits; consider raising ",
             "iter_warmup / iter_sampling."),
      format(refit_diagnostics$max_rhat, digits = 3L),
      refit_diagnostics$n_high_rhat_refits, rhat_threshold,
      refit_diagnostics$n_divergent_refits))
  }

  list(
    table             = table,
    B_ok              = as.integer(B_ok),
    seed              = if (is.null(seed)) NA_integer_ else as.integer(seed),
    warnings          = warnings_msg,
    refit_diagnostics = refit_diagnostics
  )
}

# ---------------------------------------------------------------------------
# Temporal block-bootstrap-by-refit robust inference
# ---------------------------------------------------------------------------

#' Dependence-robust standard errors via a temporal block bootstrap
#'
#' Re-estimates the uncertainty of a scalar Path 1 Empirical-Bayes fit so
#' that it is robust to temporal (serial) dependence in the data, without
#' modelling that dependence. It refits the model on \code{B} moving (or
#' circular) block bootstrap resamples of the data ordered by
#' \code{index}, and reports the bootstrap standard deviation and
#' percentile intervals of each AMM coefficient alongside the model-based
#' (Laplace / posterior) standard errors. This is the working-independence
#' + robust-variance stance of Liang & Zeger (1986): the point estimates
#' are unchanged (consistent when the mean structure is correct, not
#' efficient), only the reported uncertainty is made dependence-robust.
#'
#' @param object A scalar Path 1 fit (\code{K = 1}, \code{p = 1}): either a
#'   \code{gdpar_eb_fit} (Empirical Bayes, from \code{\link{gdpar_eb}}) or a
#'   \code{gdpar_fit} (full Bayes, from \code{\link{gdpar}}).
#' @param data The data frame originally passed to the fitting function
#'   (\code{\link{gdpar_eb}} or \code{\link{gdpar}}; the fit object deliberately
#'   does not store the data, to stay lightweight). It is resampled by
#'   contiguous blocks and the model is refit on each resample; the model
#'   specification (formula, AMM, family, prior) is recovered from
#'   \code{object$call}.
#' @param index Optional vector of length \code{n} giving the temporal
#'   ordering of the rows of \code{data}. The data are sorted by
#'   \code{order(index)} so that contiguous blocks correspond to
#'   contiguous time. When \code{NULL} (default) the natural row order is
#'   assumed to be the temporal order.
#' @param block_length The block size, one of three forms: \code{NULL}
#'   (default) uses the rate-optimal \code{max(1, round(n^(1/3)))} (Kuensch
#'   1989; Hall, Horowitz & Jing 1995); a positive integer fixes it
#'   manually; or the string \code{"auto"} selects it data-drivenly by the
#'   Politis & White (2004) automatic rule (with the Patton, Politis & White
#'   2009 correction), computed from the fitted residuals (no extra refit),
#'   falling back to the rate on a degenerate series. The chosen value and
#'   method are reported in the result (\code{block_length},
#'   \code{block_length_method}).
#' @param residual_type One of \code{"quantile"} (default; Dunn-Smyth
#'   randomized quantile residuals), \code{"response"}, \code{"pearson"} or
#'   \code{"deviance"}. Used \strong{only} when \code{block_length = "auto"},
#'   to feed the Politis-White selector; ignored otherwise.
#' @param randomize_seed Optional integer seed for the randomized quantile
#'   residuals of discrete families; used only by the \code{"auto"} selector,
#'   for a reproducible block-length choice. Ignored otherwise.
#' @param type One of \code{"moving"} (default) or \code{"circular"} block
#'   bootstrap.
#' @param B Integer scalar; the number of bootstrap refits. Defaults to
#'   199.
#' @param level Numeric scalar in (0, 1); the percentile-interval level.
#'   Defaults to 0.95.
#' @param seed Optional integer seed controlling both the block resampling
#'   and (deterministically derived) the per-refit Stan seeds, for full
#'   reproducibility.
#' @param iter_warmup,iter_sampling,chains Integer scalars controlling each
#'   refit's conditional HMC. Defaults (500, 500, 2) keep the refits short.
#' @param verbose Logical scalar; when TRUE, prints an opt-in cost message.
#' @param ... Additional arguments forwarded to \code{\link{gdpar_eb}} for
#'   every refit.
#'
#' @return A list of class \code{gdpar_dependence_robust} with components
#'   \code{table} (data frame with one row per coefficient and columns
#'   \code{estimate}, \code{model_se}, \code{robust_se}, \code{se_ratio},
#'   \code{ci_lower}, \code{ci_upper}), \code{block_length},
#'   \code{block_length_method} (\code{"rate"}, \code{"fixed"} or
#'   \code{"auto"}; \code{"rate"} also flags an \code{"auto"} request that
#'   fell back), \code{type}, \code{B}, \code{B_ok} (successful refits),
#'   \code{level}, \code{index_supplied}, \code{seed}, \code{warnings} and
#'   \code{refit_diagnostics} (aggregate per-refit convergence: \code{max_rhat},
#'   \code{min_ess_bulk}, \code{n_divergent_refits}, \code{n_high_rhat_refits}).
#'   A \code{print} method provides a human-readable summary.
#'
#' @section Honest scope:
#' The bootstrap delivers robust variance, not better point estimates, and
#' is valid for weak / short-range dependence relative to
#' \code{block_length}; it does not rescue long-memory or unit-root
#' processes. gdpar does not model the dependence (that is deferred to a
#' future block); here it only makes its inference robust to it.
#'
#' @section Empirical-Bayes and full-Bayes paths:
#' Both the scalar Empirical-Bayes (\code{gdpar_eb_fit}) and the scalar
#' full-Bayes (\code{gdpar_fit}) paths are supported (decision D102), through a
#' single shared engine; only the per-fit extraction of the point estimate, the
#' model SE and the residuals is class-dispatched. On the EB path the point
#' estimate / model SE are the Laplace / conditional-posterior mean and SD; on
#' the full-Bayes path they are the \strong{posterior mean and posterior SD} of
#' each AMM coefficient (\code{theta_ref}, \code{a_coef}, \code{b_coef},
#' \code{W_raw}, the latter on its raw scale, for parity with the EB extractor).
#' In both cases \code{robust_se} is the block-bootstrap SD of the per-refit
#' point estimate and \code{se_ratio = robust_se / model_se} is a like-for-like
#' SD-vs-SD ratio: it contrasts the dependence-robust sampling variability of
#' the point estimate against the within-model (posterior / Laplace) SD, and
#' exceeds 1 when the latter understates the former. The posterior mean / SD
#' choice (rather than median / IQR) preserves this parity and avoids an
#' undeclared normal-scaling constant. (For a strongly skewed coefficient
#' posterior -- \code{W_raw} under sparse data or strong shrinkage is the usual
#' culprit -- the posterior mean can sit off the posterior mode; inspect
#' \code{object$fit$draws()} directly in that case. The reported point estimate
#' is always the posterior mean.)
#'
#' Three honest full-Bayes caveats. (1) Each refit re-runs the \strong{full} HMC
#' (markedly more costly than an EB refit). (2) A finite-iteration refit carries
#' Monte-Carlo error in its posterior mean that slightly and
#' \emph{conservatively} inflates \code{robust_se} (reducible by a larger
#' \code{iter_sampling}; the aggregate refit ESS is reported in
#' \code{refit_diagnostics}). (3) The \code{se_ratio} has a subtly different
#' reading across paths: under an \strong{informative prior} the full-Bayes
#' posterior SD can be smaller than the bootstrap SD \emph{even under correct
#' independent specification}, because the prior concentrates the posterior
#' beyond what the data alone support, giving \code{se_ratio < 1}. That is benign
#' prior regularization, not the model SE overstating uncertainty; only
#' \code{se_ratio} clearly \strong{above} 1 signals dependence / misspecification
#' (the analytic conjugate-Gaussian check gives
#' \eqn{\mathrm{se\_ratio}^2 = n\tau / (n\tau + \tau_0) < 1} with prior precision
#' \eqn{\tau_0}). Note too that the EB and full-Bayes \code{theta_ref} point
#' estimates are different estimands -- the Laplace \strong{mode} (EB) versus the
#' posterior \strong{mean} (full Bayes) -- which coincide asymptotically
#' (Bernstein-von Mises) but may differ in finite samples; their \code{se_ratio}
#' values are therefore not expected to match to the last digit when the same
#' data are run through both paths. A widened / bagged posterior (BayesBag;
#' Huggins & Miller 2019) is a different object -- a re-architected estimator
#' rather than a robust variance for the same one -- and is a documented deferred
#' lateral, not adopted here.
#'
#' @section Dependencies:
#' Uses \pkg{cmdstanr} for the refits and \pkg{posterior} to extract the
#' coefficient estimates (Empirical-Bayes or full-Bayes).
#'
#' @references
#' Liang, K.-Y. & Zeger, S. L. (1986). Longitudinal data analysis using
#' generalized linear models. \emph{Biometrika} 73(1), 13-22.
#'
#' Kuensch, H. R. (1989). The jackknife and the bootstrap for general
#' stationary observations. \emph{Annals of Statistics} 17(3), 1217-1241.
#'
#' Politis, D. N. & White, H. (2004). Automatic block-length selection for
#' the dependent bootstrap. \emph{Econometric Reviews} 23(1), 53-70.
#'
#' Patton, A., Politis, D. N. & White, H. (2009). Correction to "Automatic
#' block-length selection for the dependent bootstrap". \emph{Econometric
#' Reviews} 28(4), 372-375.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("cmdstanr", quietly = TRUE) &&
#'     requireNamespace("posterior", quietly = TRUE)) {
#'   n <- 100
#'   x <- rnorm(n)
#'   y <- 1 + 0.5 * x + as.numeric(stats::arima.sim(list(ar = 0.6), n))
#'   df <- data.frame(x = x, y = y, t = seq_len(n))
#'   fit <- gdpar_eb(y ~ x, amm = amm_spec(a = ~ x), data = df,
#'                   chains = 2, iter_warmup = 100, iter_sampling = 100)
#'   # B is kept small here for a fast example; use B >= 199 in practice.
#'   gdpar_dependence_robust(fit, data = df, index = df$t, B = 10,
#'                           seed = 1, iter_warmup = 100,
#'                           iter_sampling = 100, chains = 2)
#'   # Data-driven block length (Politis-White), opt-in:
#'   gdpar_dependence_robust(fit, data = df, index = df$t,
#'                           block_length = "auto", B = 10, seed = 1,
#'                           iter_warmup = 100, iter_sampling = 100, chains = 2)
#' }
#' }
#'
#' @seealso \code{\link{gdpar_dependence_diagnostic}}, \code{\link{gdpar_eb}}
#' @export
gdpar_dependence_robust <- function(object, data, index = NULL,
                                    block_length = NULL,
                                    residual_type = c("quantile", "response",
                                                      "pearson", "deviance"),
                                    randomize_seed = NULL,
                                    type = c("moving", "circular"),
                                    B = 199L, level = 0.95, seed = NULL,
                                    iter_warmup = 500L, iter_sampling = 500L,
                                    chains = 2L, verbose = TRUE, ...) {
  .gdpar_assert_scalar_dep(object, "object")
  assert_data_frame(data, "data")
  residual_type <- match.arg(residual_type)
  type <- match.arg(type)
  assert_count(B, "B")
  assert_numeric_scalar(level, "level", lower = 0, upper = 1)
  assert_count(iter_warmup, "iter_warmup")
  assert_count(iter_sampling, "iter_sampling")
  assert_count(chains, "chains")
  if (!is.logical(verbose) || length(verbose) != 1L) {
    gdpar_abort("Argument 'verbose' must be a logical scalar.",
                class = "gdpar_input_error")
  }
  if (!is.null(block_length) && is.character(block_length) &&
      !.gdpar_is_auto(block_length)) {
    gdpar_abort(
      "Argument 'block_length' must be NULL, a positive integer, or \"auto\".",
      class = "gdpar_input_error"
    )
  }
  require_suggested("cmdstanr", "refit Path 1 models")
  require_suggested("posterior", "extract posterior draws")

  n <- nrow(data)
  ord <- NULL
  if (!is.null(index)) {
    if (length(index) != n) {
      gdpar_abort(
        sprintf("Argument 'index' must have length n = %d.", n),
        class = "gdpar_input_error"
      )
    }
    ord <- order(index)
    data <- data[ord, , drop = FALSE]
  }
  block_length_method <- "fixed"
  if (is.null(block_length)) {
    block_length <- .gdpar_default_block_length(n)
    block_length_method <- "rate"
  } else if (.gdpar_is_auto(block_length)) {
    resid <- .gdpar_dependence_residuals(object, residual_type, randomize_seed)
    if (!is.null(ord)) resid <- resid[ord]
    pw <- .gdpar_politis_white_block_length(resid)
    block_length <- pw$block_length
    block_length_method <- pw$method
    if (verbose) {
      gdpar_inform(paste0("block_length = \"auto\": ", pw$reason),
                   class = "gdpar_optin_message")
    }
  } else {
    block_length <- as.integer(block_length)
  }
  if (block_length < 1L || block_length > n) {
    gdpar_abort(
      sprintf("Argument 'block_length' must be in [1, %d].", n),
      class = "gdpar_input_error"
    )
  }

  is_fb <- inherits(object, "gdpar_fit") && !inherits(object, "gdpar_eb_fit")
  verbose_msg <- sprintf(
    paste(
      "gdpar_dependence_robust is an opt-in, computationally expensive",
      "audit: it will perform %d block-bootstrap refit(s) (block length",
      "%d, %s scheme)%s. It does not change the point estimates; it",
      "re-estimates their uncertainty to be robust to temporal",
      "dependence."
    ),
    B, block_length, type,
    if (is_fb) ", each a full-Bayes refit (full HMC; markedly more expensive than the Empirical-Bayes path)" else ""
  )

  resample_fun <- function() {
    .gdpar_block_bootstrap_data_indices(n, block_length, type)
  }

  caller_env <- parent.frame()
  core <- .gdpar_dependence_robust_engine(
    object, data, resample_fun, B = B, level = level, seed = seed,
    iter_warmup = iter_warmup, iter_sampling = iter_sampling,
    chains = chains, verbose = verbose, verbose_msg = verbose_msg,
    caller_env = caller_env, ...
  )

  obj <- list(
    table               = core$table,
    block_length        = block_length,
    block_length_method = block_length_method,
    type                = type,
    B                   = as.integer(B),
    B_ok                = core$B_ok,
    level               = level,
    index_supplied      = !is.null(index),
    seed                = core$seed,
    warnings            = core$warnings,
    refit_diagnostics   = core$refit_diagnostics
  )
  class(obj) <- c("gdpar_dependence_robust", "list")
  obj
}

#' Print method for gdpar_dependence_robust objects
#'
#' @param x A \code{gdpar_dependence_robust} object.
#' @param digits Integer; significant digits for the printed table.
#' @param ... Unused; present for S3 generic compatibility.
#' @return Invisibly returns \code{x}.
#' @export
print.gdpar_dependence_robust <- function(x, digits = 3L, ...) {
  cat("<gdpar_dependence_robust> ", x$type, " block bootstrap\n", sep = "")
  bl_method <- x$block_length_method %||% "fixed"
  bl_label <- switch(bl_method,
    auto  = " (auto: Politis-White)",
    rate  = " (rate: n^(1/3))",
    fixed = "",
    "")
  cat("  block length    : ", x$block_length, bl_label,
      "   B = ", x$B, " (", x$B_ok, " ok)",
      if (x$index_supplied) "   [index supplied]" else
        "   [natural row order]", "\n", sep = "")
  cat("  level           : ", x$level, "\n\n", sep = "")
  tab <- x$table
  num <- vapply(tab, is.numeric, logical(1L))
  tab[num] <- lapply(tab[num], function(col) format(col, digits = digits))
  print(tab, row.names = FALSE)
  cat("\n  se_ratio = robust_se / model_se; > 1 means the model-based\n",
      "  uncertainty understates the dependence-robust uncertainty.\n",
      sep = "")
  .gdpar_print_refit_diagnostics(x$refit_diagnostics, digits)
  if (length(x$warnings) > 0L) {
    cat("\nNotes (", length(x$warnings), "):\n", sep = "")
    for (w in utils::head(x$warnings, 5L)) cat("  - ", w, "\n", sep = "")
    if (length(x$warnings) > 5L) {
      cat("  ... and ", length(x$warnings) - 5L, " more.\n", sep = "")
    }
  }
  invisible(x)
}

#' Internal: print the aggregate per-refit convergence diagnostics, if any.
#'
#' @keywords internal
#' @noRd
.gdpar_print_refit_diagnostics <- function(rd, digits = 3L) {
  if (is.null(rd)) return(invisible(NULL))
  mr <- rd$max_rhat %||% NA_real_
  if (!is.finite(mr)) return(invisible(NULL))
  cat(sprintf(
    paste0("  refit convergence: max R-hat = %s, min ESS = %s; ",
           "%d divergent, %d R-hat > %.2f refit(s)\n"),
    format(mr, digits = digits),
    format(rd$min_ess_bulk %||% NA_real_, digits = digits),
    rd$n_divergent_refits %||% 0L, rd$n_high_rhat_refits %||% 0L,
    rd$rhat_threshold %||% 1.05))
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Spatial dependence: internal helpers (base R; no spdep / sf / geosphere)
# ---------------------------------------------------------------------------

#' Internal: variance-optimal default cells-per-axis for the spatial bootstrap.
#'
#' Returns \code{max(2, round(n^(1/4)))}, the \eqn{d = 2} specialisation of the
#' variance-MSE-optimal block rate \eqn{M \sim n^{d/(d+2)}} points per block (so
#' \eqn{g^2 \sim n^{1/2}} cells, \eqn{M \sim n^{1/2}} points per cell). See the
#' \code{@details} of \code{\link{gdpar_spatial_dependence_robust}} for the
#' derivation, the reduction to the \eqn{n^{1/3}} temporal rate at \eqn{d = 1},
#' and the registered dissent (decision D100).
#'
#' @keywords internal
#' @noRd
.gdpar_spatial_default_g <- function(n) {
  max(2L, as.integer(round(n^(1 / 4))))
}

#' Internal: validate and coerce a 2-column numeric coordinate matrix.
#'
#' @keywords internal
#' @noRd
.gdpar_validate_coords <- function(coords, n, arg = "coords") {
  if (is.data.frame(coords)) coords <- as.matrix(coords)
  if (!is.matrix(coords) || !is.numeric(coords)) {
    gdpar_abort(
      sprintf("Argument '%s' must be a numeric matrix or data frame.", arg),
      class = "gdpar_input_error"
    )
  }
  if (ncol(coords) != 2L) {
    gdpar_abort(
      sprintf("Argument '%s' must have exactly 2 columns (x, y).", arg),
      class = "gdpar_input_error"
    )
  }
  if (nrow(coords) != n) {
    gdpar_abort(
      sprintf("Argument '%s' must have n = %d rows (one per observation).",
              arg, n),
      class = "gdpar_input_error"
    )
  }
  if (any(!is.finite(coords))) {
    gdpar_abort(
      sprintf("Argument '%s' contains non-finite values.", arg),
      class = "gdpar_input_error"
    )
  }
  coords
}

#' Internal: binary k-nearest-neighbour adjacency (Euclidean; ties broken by
#' index, so duplicate locations are well-defined but spatially degenerate).
#'
#' @keywords internal
#' @noRd
.gdpar_knn_adjacency <- function(coords, k) {
  n <- nrow(coords)
  D <- as.matrix(stats::dist(coords))
  W <- matrix(0, n, n)
  for (i in seq_len(n)) {
    di <- D[i, ]
    di[i] <- Inf
    nn <- order(di)[seq_len(k)]
    W[i, nn] <- 1
  }
  W
}

#' Internal: binary distance-band adjacency. The band is the smallest threshold
#' that leaves no location isolated (the max over points of the
#' nearest-neighbour distance), a declared data-driven heuristic.
#'
#' @keywords internal
#' @noRd
.gdpar_distance_band_adjacency <- function(coords) {
  n <- nrow(coords)
  D <- as.matrix(stats::dist(coords))
  diag(D) <- Inf
  threshold <- max(apply(D, 1L, min))
  W <- (D <= threshold) * 1
  diag(W) <- 0
  W
}

#' Internal: Moran's I for residuals \code{resid} under (possibly asymmetric)
#' weights \code{W}. \eqn{I = (n / S_0) \sum_{ij} w_{ij} z_i z_j / \sum_i z_i^2}
#' with \eqn{z = r - \bar r} and \eqn{S_0 = \sum_{ij} w_{ij}}.
#'
#' @keywords internal
#' @noRd
.gdpar_morans_i <- function(resid, W, S0 = sum(W)) {
  n <- length(resid)
  z <- resid - mean(resid)
  (n / S0) * (sum(z * as.numeric(W %*% z)) / sum(z^2))
}

#' Internal: length-n resampled row indices for a spatial block bootstrap.
#'
#' \code{"tiled"} partitions the bounding box into a \code{g x g} grid (cell side
#' = range / g per axis); with \code{random_origin} the grid is shifted by a
#' \eqn{U(0, 1) \times} cell-side draw per replicate (Politis-Romano-Lahiri
#' randomized partition), then non-empty cells are sampled with replacement and
#' concatenated until \eqn{\geq n} indices accrue, truncated to \code{n}.
#' \code{"moving"} draws overlapping square blocks (side = cell side) anchored to
#' cover a sampled observation, so every block is non-empty. Sibling of
#' \code{.gdpar_block_bootstrap_data_indices()} for 2-D coordinates.
#'
#' @keywords internal
#' @noRd
.gdpar_spatial_block_indices <- function(coords, g, scheme, random_origin,
                                         mins, ranges) {
  n <- nrow(coords)
  L <- ranges / g
  if (scheme == "tiled") {
    origin <- if (random_origin) mins - stats::runif(2L) * L else mins
    cx <- floor((coords[, 1L] - origin[1L]) / L[1L])
    cy <- floor((coords[, 2L] - origin[2L]) / L[2L])
    groups <- split(seq_len(n), paste(cx, cy, sep = "_"))
    n_cells <- length(groups)
    out <- integer(0L)
    while (length(out) < n) {
      out <- c(out, groups[[sample.int(n_cells, 1L)]])
    }
  } else {
    out <- integer(0L)
    while (length(out) < n) {
      seed_pt <- coords[sample.int(n, 1L), ]
      origin <- seed_pt - stats::runif(2L) * L
      inblock <- which(
        coords[, 1L] >= origin[1L] & coords[, 1L] < origin[1L] + L[1L] &
          coords[, 2L] >= origin[2L] & coords[, 2L] < origin[2L] + L[2L]
      )
      out <- c(out, inblock)
    }
  }
  out[seq_len(n)]
}

#' Internal: data-driven spatial block size by calibration over the cell grid.
#'
#' The spatial counterpart of \code{\link{.gdpar_politis_white_block_length}}.
#' Politis & White (2004) has \strong{no} established spatial plug-in (the
#' flat-top spectral-density-at-zero estimator does not extend cleanly to a
#' random field in R^2), so the cells-per-axis \eqn{g} is chosen by a
#' \emph{data-driven calibration over a grid of \eqn{g}}, declared as a
#' heuristic (decision D101; cross-lineage-adjudicated). The construction:
#'
#' \enumerate{
#'   \item \strong{Surrogate} (the influence directions, not the raw residual
#'     mean): for each candidate \eqn{g}, draw \code{B0} spatial block
#'     resamples and form the design-weighted residual functionals
#'     \eqn{t = (1/n)\,D^\top z} where \eqn{D = [1,\ \tilde{gx},\ \tilde{gy}]}
#'     (intercept plus standardized coordinates) and \eqn{z} the centered
#'     residuals. \eqn{V(g) = \mathrm{tr}\,\widehat{\mathrm{Cov}}_{boot}(t)} is
#'     the bootstrap variance the robust SE actually depends on (to first order
#'     the score / influence function of the coefficient, so its MSE-optimal
#'     \eqn{g} matches the coefficient's).
#'   \item \strong{Criterion} (direct empirical MSE): smooth \eqn{V(g)} with a
#'     width-3 running median; bias is anchored at \strong{the largest blocks},
#'     \eqn{g_{ref} = g_{\min}}, which are the \emph{least} biased for variance
#'     estimation (the dependence-breaking bias of the block-bootstrap variance
#'     estimator grows like \eqn{g/\sqrt n}); the variance term uses the
#'     \emph{sampling} variance of the estimator, whose relative variance scales
#'     like the inverse number of blocks (Lahiri 2003), \eqn{\propto V(g)^2 /
#'     n_{tiles}(g)}; \eqn{g^\ast = \arg\min_g [(\tilde V(g) - \tilde
#'     V(g_{\min}))^2 + c\,V(g)^2 / n_{tiles}(g)]}.
#' }
#'
#' \code{var_const} (\eqn{c}) is a declared calibration constant fixing the
#' bias/variance operating point. Two corrections to the cross-lineage proposal,
#' both validated empirically (decision D101): (i) the bias anchor is at
#' \eqn{g_{\min}} (largest blocks), not \eqn{g_{\max}} -- anchoring at
#' \eqn{g_{\max}} steers toward the smallest blocks and anticonservative SEs (the
#' audited blocker); (ii) the variance term is the \eqn{1/n_{tiles}} sampling
#' variance, not the leave-one-out jackknife over the \code{B0} replicates -- the
#' jackknife estimates the finite-\code{B0} Monte-Carlo noise, which vanishes as
#' \code{B0} grows and leaves the bias term to collapse the choice onto
#' \eqn{g_{\min}} (the selector then never adapts). Fall back to the
#' \eqn{n^{1/4}} rate for tiny \code{n}, fewer than three valid grid points, a
#' non-finite criterion, or an argmin at the smallest-block (anticonservative)
#' grid boundary.
#'
#' @param coords Validated \code{n x 2} coordinate matrix.
#' @param resid Residual vector, row-aligned with \code{coords}.
#' @param scheme,random_origin,mins,ranges As in
#'   \code{\link{.gdpar_spatial_block_indices}} (the calibration resamples with
#'   the same scheme as the main bootstrap).
#' @param B0 Integer number of cheap (no-refit) resamples per candidate.
#' @param var_const Declared constant \eqn{c} of the \eqn{1/n_{tiles}} variance
#'   term (the bias/variance operating point).
#' @param seed Optional integer seed for the calibration resamples.
#' @return A list with \code{block_size} (integer \eqn{g}), \code{method}
#'   (\code{"auto"} or \code{"rate"}) and a \code{reason}.
#' @keywords internal
#' @noRd
.gdpar_spatial_block_length_auto <- function(coords, resid, scheme,
                                             random_origin, mins, ranges,
                                             B0 = 200L, var_const = 1,
                                             seed = NULL) {
  n <- length(resid)
  g_def <- .gdpar_spatial_default_g(n)
  fb <- function(reason) list(block_size = g_def, method = "rate", reason = reason)
  if (n < 25L) {
    return(fb(sprintf(
      "n = %d < 25 too small for spatial auto-selection; used n^(1/4) rate.", n)))
  }
  z <- resid - mean(resid)
  gx <- coords[, 1L]; gy <- coords[, 2L]
  sx <- stats::sd(gx); sy <- stats::sd(gy)
  if (!is.finite(sx) || sx <= 0 || !is.finite(sy) || sy <= 0) {
    return(fb("Degenerate coordinate spread; used n^(1/4) rate."))
  }
  Dsurr <- cbind(1, (gx - mean(gx)) / sx, (gy - mean(gy)) / sy)
  g_lo <- max(2L, as.integer(floor(0.5 * g_def)))
  # Upper cap is the cell-occupancy bound (>= ~3 points per cell), the binding
  # validity constraint for within-cell resampling, not the tighter n^(1/3).
  g_hi <- min(as.integer(floor(3 * g_def)), as.integer(floor(sqrt(n / 3))))
  if (g_hi <= g_lo) g_hi <- g_lo + 1L
  grid <- unique(round(exp(seq(log(g_lo), log(g_hi), length.out = 6L))))
  grid <- as.integer(grid[grid >= 2L & (n / grid^2) >= 3])
  if (length(grid) < 3L) {
    return(fb("Fewer than three valid grid points; used n^(1/4) rate."))
  }
  if (!is.null(seed)) set.seed(as.integer(seed))
  Vg <- numeric(length(grid))
  ntiles <- numeric(length(grid))
  for (gi in seq_along(grid)) {
    g <- grid[gi]
    Tmat <- matrix(NA_real_, B0, 3L)
    for (b in seq_len(B0)) {
      idx <- .gdpar_spatial_block_indices(coords, g, scheme, random_origin,
                                          mins, ranges)
      Tmat[b, ] <- crossprod(Dsurr[idx, , drop = FALSE], z[idx]) / n
    }
    Vg[gi] <- sum(apply(Tmat, 2L, stats::var))
    L0 <- ranges / g
    cell <- paste(pmin(floor((coords[, 1L] - mins[1L]) / L0[1L]), g - 1L),
                  pmin(floor((coords[, 2L] - mins[2L]) / L0[2L]), g - 1L),
                  sep = "_")
    ntiles[gi] <- length(unique(cell))
  }
  if (any(!is.finite(Vg)) || all(Vg <= 0) || any(ntiles < 1)) {
    return(fb("Non-finite or zero bootstrap variance; used n^(1/4) rate."))
  }
  Vt <- as.numeric(stats::runmed(Vg, k = 3L))
  bias2 <- (Vt - Vt[1L])^2                 # anchor g_ref = g_min (largest blocks)
  # Variance of the block-bootstrap variance estimator scales like the inverse
  # number of blocks (relative variance ~ 1/n_blocks; Lahiri 2003), NOT the
  # finite-B0 Monte-Carlo noise: the latter vanishes as B0 grows and would leave
  # the bias term to collapse the choice onto g_min. var_const is declared.
  var_term <- var_const * Vt^2 / ntiles
  mse <- bias2 + var_term
  if (any(!is.finite(mse))) {
    return(fb("Non-finite MSE criterion; used n^(1/4) rate."))
  }
  jstar <- which.min(mse)
  if (jstar == length(grid)) {
    return(fb(paste("MSE minimised at the smallest-block grid boundary",
                    "(anticonservative risk); used n^(1/4) rate.")))
  }
  g_star <- as.integer(grid[jstar])
  list(block_size = g_star, method = "auto",
       reason = sprintf(
         "Data-driven MSE over g (grid %s; anchor g_min; B0 = %d, c = %g): g* = %d (rate g = %d).",
         paste(grid, collapse = ","), as.integer(B0), var_const, g_star, g_def))
}

# ---------------------------------------------------------------------------
# Spatial residual dependence diagnostic
# ---------------------------------------------------------------------------

#' Spatial residual dependence diagnostic for a scalar Empirical-Bayes fit
#'
#' Quantifies spatial autocorrelation in the residuals of a fitted scalar Path 1
#' Empirical-Bayes model via Moran's I over a spatial weight structure. gdpar
#' assumes conditional independence; under spatial dependence that assumption is
#' violated and the model-based (posterior / Laplace) uncertainty is too narrow.
#' This is the spatial sibling of \code{\link{gdpar_dependence_diagnostic}} and
#' the natural gate for \code{\link{gdpar_spatial_dependence_robust}}.
#'
#' @param object A scalar Path 1 fit (\code{K = 1}, \code{p = 1}): either a
#'   \code{gdpar_eb_fit} (Empirical Bayes, from \code{\link{gdpar_eb}}) or a
#'   \code{gdpar_fit} (full Bayes, from \code{\link{gdpar}}).
#' @param coords A numeric \code{n x 2} matrix or data frame of spatial
#'   coordinates, row-aligned with the training data.
#' @param W Optional user-supplied \code{n x n} spatial weight matrix. When
#'   given it overrides \code{weights}/\code{k} and is row-standardized
#'   internally (its diagonal is zeroed). Supplying \code{W} is the right choice
#'   when domain knowledge (adjacency, flow, network connectivity) defines the
#'   neighbourhood structure.
#' @param weights One of \code{"knn"} (default; \code{k}-nearest-neighbour
#'   adjacency, robust to irregular spacing) or \code{"distance"} (a
#'   distance-band whose threshold is the smallest that isolates no location).
#'   Ignored when \code{W} is supplied. Both are row-standardized.
#' @param k Integer number of neighbours for \code{weights = "knn"}. When
#'   \code{NULL} (default) the declared heuristic \code{max(4, min(round(log n),
#'   n - 1))} is used.
#' @param residual_type One of \code{"quantile"} (default; randomized quantile /
#'   Dunn-Smyth residuals), \code{"response"}, \code{"pearson"} or
#'   \code{"deviance"}.
#' @param test One of \code{"permutation"} (default; \code{n_perm}
#'   location-relabelling permutations, two-sided via \eqn{|I - E[I]|}, robust to
#'   non-normal residuals and asymmetric \code{W}) or \code{"analytic"} (the
#'   Cliff-Ord normal approximation, cheaper but assuming a symmetric \code{W} --
#'   a warning is emitted otherwise).
#' @param n_perm Integer number of permutations for the permutation test
#'   (default 999; capped below \code{n!} for tiny \code{n}).
#' @param level Numeric scalar in (0, 1); the confidence level used to turn the
#'   p-value into the verdict. Defaults to 0.95.
#' @param randomize_seed Optional integer seed for the randomized quantile
#'   residuals of discrete families; ignored otherwise.
#' @param seed Optional integer seed for the permutation test, for
#'   reproducibility.
#' @param ... Unused; present for signature stability.
#'
#' @return A list of class \code{gdpar_spatial_dependence_diagnostic} with
#'   components \code{residual_type}, \code{n}, \code{weights}, \code{k},
#'   \code{style}, \code{n_zero_weight}, \code{morans_i}, \code{expected_i}
#'   (\eqn{-1/(n-1)}), \code{var_i} (analytic, else \code{NA}), \code{z}
#'   (analytic, else \code{NA}), \code{p_value}, \code{test}, \code{n_perm},
#'   \code{level} and \code{verdict}. A \code{print} method is provided.
#'
#' @details
#' Moran's I is hand-rolled in base R (no \pkg{spdep} / \pkg{sf} dependency).
#' With row-standardized weights \eqn{S_0 = n} and
#' \eqn{E[I] = -1/(n-1)} under the null of spatial exchangeability.
#'
#' \strong{Guards and caveats.} Locations with zero total weight (isolated under
#' a supplied \code{W} or a too-small \code{k}) make Moran's I undefined: a
#' warning is emitted and \code{morans_i} is returned as \code{NA} (kNN with
#' \eqn{k \geq 1} never isolates a point). Duplicate coordinates are permitted
#' (kNN ties broken by index; spatially degenerate but well-defined). For
#' \code{n < 20} a hard and for \code{n < 50} a soft small-sample warning is
#' issued. Coordinates are treated as \strong{Euclidean}: lon/lat data should be
#' projected first (e.g. UTM), or the neighbour graph is distorted, severely so
#' at high latitudes -- great-circle distance is deliberately not supported to
#' avoid a heavy \pkg{geosphere}/\pkg{sf} dependency. Finally, a significant
#' Moran's I may reflect \strong{either} true spatial dependence \strong{or}
#' model misspecification (e.g. an omitted nonlinear covariate effect); the
#' diagnostic tests residual spatial exchangeability, not its cause.
#'
#' @section Honest scope:
#' gdpar does \strong{not} model the spatial dependence (that is Axis 1 / a
#' future block, deferred and evidence-gated). This diagnostic only makes the
#' violation visible; the companion remedy
#' \code{\link{gdpar_spatial_dependence_robust}} re-estimates the uncertainty by
#' a spatial block bootstrap.
#'
#' @references
#' Moran, P. A. P. (1950). Notes on continuous stochastic phenomena.
#' \emph{Biometrika} 37(1/2), 17-23.
#'
#' Cliff, A. D. & Ord, J. K. (1981). \emph{Spatial Processes: Models and
#' Applications}. Pion, London.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("cmdstanr", quietly = TRUE) &&
#'     requireNamespace("posterior", quietly = TRUE)) {
#'   n <- 100
#'   gx <- runif(n); gy <- runif(n)
#'   x <- rnorm(n)
#'   y <- 1 + 0.5 * x + (gx + gy) + rnorm(n)  # smooth spatial trend in residuals
#'   df <- data.frame(x = x, y = y)
#'   fit <- gdpar_eb(y ~ x, amm = amm_spec(a = ~ x), data = df,
#'                   chains = 2, iter_warmup = 100, iter_sampling = 100)
#'   gdpar_spatial_dependence_diagnostic(fit, coords = cbind(gx, gy), seed = 1)
#' }
#' }
#'
#' @seealso \code{\link{gdpar_spatial_dependence_robust}},
#'   \code{\link{gdpar_dependence_diagnostic}}, \code{\link{gdpar_eb}}
#' @export
gdpar_spatial_dependence_diagnostic <- function(object, coords, W = NULL,
                                                weights = c("knn", "distance"),
                                                k = NULL,
                                                residual_type = c("quantile",
                                                                  "response",
                                                                  "pearson",
                                                                  "deviance"),
                                                test = c("permutation",
                                                         "analytic"),
                                                n_perm = 999L, level = 0.95,
                                                randomize_seed = NULL,
                                                seed = NULL, ...) {
  .gdpar_assert_scalar_dep(object, "object")
  weights <- match.arg(weights)
  residual_type <- match.arg(residual_type)
  test <- match.arg(test)
  assert_count(n_perm, "n_perm")
  assert_numeric_scalar(level, "level", lower = 0, upper = 1)
  require_suggested("posterior", "extract posterior draws for residuals")

  resid <- .gdpar_dependence_residuals(object, residual_type, randomize_seed)
  n <- length(resid)
  coords <- .gdpar_validate_coords(coords, n)
  if (stats::var(resid) <= 0) {
    gdpar_abort("Residuals have zero variance; cannot assess dependence.",
                class = "gdpar_diagnostic_error")
  }

  w_supplied <- !is.null(W)
  if (w_supplied) {
    if (is.data.frame(W)) W <- as.matrix(W)
    if (!is.matrix(W) || !is.numeric(W) || nrow(W) != n || ncol(W) != n) {
      gdpar_abort(
        sprintf("Supplied 'W' must be a numeric n x n matrix (n = %d).", n),
        class = "gdpar_input_error"
      )
    }
    if (any(!is.finite(W))) {
      gdpar_abort("Supplied 'W' contains non-finite values.",
                  class = "gdpar_input_error")
    }
    W_raw <- W
    diag(W_raw) <- 0
    k_used <- NA_integer_
  } else if (weights == "knn") {
    if (is.null(k)) k <- max(4L, min(as.integer(round(log(n))), n - 1L))
    k <- as.integer(k)
    if (k < 1L || k > n - 1L) {
      gdpar_abort(sprintf("Argument 'k' must be in [1, %d].", n - 1L),
                  class = "gdpar_input_error")
    }
    W_raw <- .gdpar_knn_adjacency(coords, k)
    k_used <- k
  } else {
    W_raw <- .gdpar_distance_band_adjacency(coords)
    k_used <- NA_integer_
  }

  rs <- rowSums(W_raw)
  n_zero <- sum(rs == 0)
  morans_i <- expected_i <- var_i <- z_val <- p_value <- NA_real_
  n_perm_eff <- NA_integer_
  verdict <- NA_character_

  if (n_zero > 0L) {
    warning(sprintf(
      paste("%d location(s) have zero spatial weight; Moran's I is undefined.",
            "Increase 'k', widen the band, or revise the supplied 'W'."),
      n_zero
    ), call. = FALSE)
    verdict <- sprintf("Undefined: %d location(s) have no spatial neighbours.",
                       n_zero)
  } else {
    Wn <- W_raw / rs
    S0 <- sum(Wn)
    expected_i <- -1 / (n - 1)
    morans_i <- .gdpar_morans_i(resid, Wn, S0 = S0)

    if (test == "permutation") {
      max_distinct <- if (n <= 10L) factorial(n) else Inf
      n_perm_eff <- as.integer(min(as.integer(n_perm),
                                   max(1, max_distinct - 1)))
      if (n < 20L) {
        warning(sprintf(
          paste("n = %d is very small for a spatial permutation test;",
                "treat the p-value as indicative only."), n
        ), call. = FALSE)
      } else if (n < 50L) {
        warning(sprintf(
          "n = %d is small; the spatial permutation p-value is approximate.", n
        ), call. = FALSE)
      }
      if (!is.null(seed)) set.seed(as.integer(seed))
      dev_obs <- abs(morans_i - expected_i)
      ge <- 0L
      for (p in seq_len(n_perm_eff)) {
        i_perm <- .gdpar_morans_i(resid[sample.int(n)], Wn, S0 = S0)
        if (abs(i_perm - expected_i) >= dev_obs) ge <- ge + 1L
      }
      p_value <- (1 + ge) / (n_perm_eff + 1)
    } else {
      if (!isSymmetric(unname(Wn))) {
        warning(paste(
          "Analytic Moran's I variance assumes a symmetric W; the",
          "row-standardized kNN / band weights are asymmetric, so the normal",
          "approximation is only indicative. Prefer test = \"permutation\"."
        ), call. = FALSE)
      }
      S1 <- 0.5 * sum((Wn + t(Wn))^2)
      S2 <- sum((rowSums(Wn) + colSums(Wn))^2)
      var_i <- (n^2 * S1 - n * S2 + 3 * S0^2) / (S0^2 * (n^2 - 1)) -
        expected_i^2
      if (!is.finite(var_i) || var_i <= 0) {
        warning(paste("Analytic Moran's I variance is non-positive;",
                      "use test = \"permutation\"."), call. = FALSE)
        var_i <- NA_real_
      } else {
        z_val <- (morans_i - expected_i) / sqrt(var_i)
        p_value <- 2 * stats::pnorm(-abs(z_val))
      }
    }

    alpha <- 1 - level
    if (is.na(p_value)) {
      verdict <- "Spatial dependence test undefined for the chosen weights/test."
    } else if (p_value < alpha) {
      verdict <- sprintf(
        paste("Spatial dependence detected (Moran's I = %.3g, p = %.3g <",
              "%.3g): the model-based uncertainty is not robust to it. Use",
              "gdpar_spatial_dependence_robust()."),
        morans_i, p_value, alpha
      )
    } else {
      verdict <- sprintf(
        paste("No evidence against spatial independence at level %.2f",
              "(Moran's I = %.3g, p = %.3g)."),
        level, morans_i, p_value
      )
    }
  }

  obj <- list(
    residual_type = residual_type,
    n             = n,
    weights       = if (w_supplied) "user" else weights,
    k             = k_used,
    style         = "W",
    n_zero_weight = as.integer(n_zero),
    morans_i      = morans_i,
    expected_i    = expected_i,
    var_i         = var_i,
    z             = z_val,
    p_value       = p_value,
    test          = test,
    n_perm        = n_perm_eff,
    level         = level,
    verdict       = verdict
  )
  class(obj) <- c("gdpar_spatial_dependence_diagnostic", "list")
  obj
}

#' Print method for gdpar_spatial_dependence_diagnostic objects
#'
#' @param x A \code{gdpar_spatial_dependence_diagnostic} object.
#' @param digits Integer; significant digits for the printed statistics.
#' @param ... Unused; present for S3 generic compatibility.
#' @return Invisibly returns \code{x}.
#' @export
print.gdpar_spatial_dependence_diagnostic <- function(x, digits = 3L, ...) {
  cat("<gdpar_spatial_dependence_diagnostic>\n", sep = "")
  cat("  residual type   : ", x$residual_type, "\n", sep = "")
  cat("  n               : ", x$n, "\n", sep = "")
  wlab <- if (identical(x$weights, "knn")) {
    sprintf("knn (k = %d), row-standardized", x$k)
  } else if (identical(x$weights, "distance")) {
    "distance band, row-standardized"
  } else {
    "user-supplied, row-standardized"
  }
  cat("  weights         : ", wlab, "\n", sep = "")
  if (is.na(x$morans_i)) {
    cat("  Moran's I       : undefined\n", sep = "")
  } else {
    cat("  Moran's I       : ", format(x$morans_i, digits = digits),
        "  (E[I] = ", format(x$expected_i, digits = digits), ")\n", sep = "")
    if (identical(x$test, "analytic")) {
      cat("  analytic z      : ", format(x$z, digits = digits),
          "  (p = ", format(x$p_value, digits = digits), ")\n", sep = "")
    } else {
      cat("  permutation     : ", x$n_perm, " draws  (p = ",
          format(x$p_value, digits = digits), ")\n", sep = "")
    }
  }
  cat("  Verdict         : ", x$verdict, "\n", sep = "")
  invisible(x)
}

# ---------------------------------------------------------------------------
# Spatial block-bootstrap-by-refit robust inference
# ---------------------------------------------------------------------------

#' Dependence-robust standard errors via a spatial block bootstrap
#'
#' Re-estimates the uncertainty of a scalar Path 1 Empirical-Bayes fit so that
#' it is robust to spatial dependence in the data, without modelling that
#' dependence. It refits the model on \code{B} spatial block-bootstrap resamples
#' (tiled or moving blocks over \code{coords}) and reports the bootstrap standard
#' deviation and percentile intervals of each AMM coefficient alongside the
#' model-based (Laplace / posterior) standard errors. As in
#' \code{\link{gdpar_dependence_robust}} (its temporal sibling, with which it
#' shares one refit engine), this is the working-independence + robust-variance
#' stance of Liang & Zeger (1986): the point estimates are unchanged, only the
#' reported uncertainty is made dependence-robust.
#'
#' @param object A scalar Path 1 fit (\code{K = 1}, \code{p = 1}): either a
#'   \code{gdpar_eb_fit} (Empirical Bayes, from \code{\link{gdpar_eb}}) or a
#'   \code{gdpar_fit} (full Bayes, from \code{\link{gdpar}}).
#' @param data The data frame originally passed to the fitting function
#'   (\code{\link{gdpar_eb}} or \code{\link{gdpar}}), row-aligned with
#'   \code{coords}. It is resampled by spatial blocks and the model is refit on
#'   each resample; the specification is recovered from \code{object$call}.
#' @param coords A numeric \code{n x 2} matrix or data frame of spatial
#'   coordinates, row-aligned with \code{data}.
#' @param block_size The number of grid cells per axis \code{g} (the block is a
#'   square of side range/\code{g}), one of three forms: \code{NULL} (default)
#'   uses the variance-optimal rate \code{max(2, round(n^(1/4)))} (see Details);
#'   a positive integer fixes it manually; or the string \code{"auto"} selects
#'   it data-drivenly by a calibration over a grid of \code{g} (no spatial
#'   plug-in exists; decision D101), computed from the fitted residuals (no
#'   extra refit), falling back to the rate when the calibration degenerates.
#'   The chosen value and method are reported (\code{block_size},
#'   \code{block_size_method}).
#' @param residual_type One of \code{"quantile"} (default; Dunn-Smyth
#'   randomized quantile residuals), \code{"response"}, \code{"pearson"} or
#'   \code{"deviance"}. Used \strong{only} when \code{block_size = "auto"}, to
#'   feed the data-driven selector; ignored otherwise.
#' @param randomize_seed Optional integer seed for the randomized quantile
#'   residuals of discrete families; used only by the \code{"auto"} selector,
#'   for a reproducible block-size choice. Ignored otherwise.
#' @param scheme One of \code{"tiled"} (default; non-overlapping cells) or
#'   \code{"moving"} (overlapping square blocks anchored on sampled points).
#' @param random_origin Logical; when \code{TRUE} (default) and
#'   \code{scheme = "tiled"}, the grid origin is randomized per replicate
#'   (Politis-Romano-Lahiri), breaking the deterministic boundary artifact at
#'   the cost of one extra random draw per refit.
#' @param B Integer number of bootstrap refits. Defaults to 199.
#' @param level Numeric scalar in (0, 1); the percentile-interval level.
#'   Defaults to 0.95.
#' @param seed Optional integer seed controlling the block resampling and the
#'   per-refit Stan seeds, for reproducibility.
#' @param iter_warmup,iter_sampling,chains Integer scalars controlling each
#'   refit's conditional HMC. Defaults (500, 500, 2) keep the refits short.
#' @param verbose Logical scalar; when \code{TRUE}, prints an opt-in cost
#'   message.
#' @param ... Additional arguments (currently absorbed; reserved for forward
#'   compatibility, matching the temporal sibling).
#'
#' @return A list of class \code{gdpar_spatial_dependence_robust} with
#'   components \code{table} (one row per coefficient with \code{estimate},
#'   \code{model_se}, \code{robust_se}, \code{se_ratio}, \code{ci_lower},
#'   \code{ci_upper}), \code{block_size}, \code{block_size_method}
#'   (\code{"rate"}, \code{"fixed"} or \code{"auto"}; \code{"rate"} also flags
#'   an \code{"auto"} request that fell back), \code{scheme},
#'   \code{random_origin}, \code{n_tiles}, \code{B}, \code{B_ok}, \code{level},
#'   \code{seed}, \code{warnings} and \code{refit_diagnostics} (aggregate
#'   per-refit convergence, as in \code{\link{gdpar_dependence_robust}}). A
#'   \code{print} method is provided.
#'
#' @details
#' \strong{Default block-size rate (decision D100).} The block side per axis is
#' \code{g = max(2, round(n^(1/4)))}. This is the \eqn{d = 2} case of the rate
#' that minimises the mean-squared error of the block-bootstrap \emph{variance}
#' estimator. Writing \eqn{M} for the number of points per block (linear extent
#' \eqn{M^{1/d}} per axis), the first-order bias from dependence broken at block
#' edges is \eqn{O(M^{-1/d})} (Kuensch 1989; Hall, Horowitz & Jing 1995) and the
#' estimator variance is \eqn{O(M/n)}, so
#' \eqn{\mathrm{MSE}(M) \sim M^{-2/d} + M/n} is minimised at
#' \eqn{M \sim n^{d/(d+2)}}. At \eqn{d = 1} this gives \eqn{M \sim n^{1/3}}
#' points per block, \strong{exactly} the \eqn{n^{1/3}} block length of the
#' temporal default (\code{\link{gdpar_dependence_robust}}); at \eqn{d = 2} it
#' gives \eqn{M \sim n^{1/2}} points per block, i.e. \eqn{g^2 = n/M \sim n^{1/2}}
#' cells, hence \eqn{g \sim n^{1/4}} cells per axis. The exponent is therefore
#' the variance-optimal rate that reduces correctly to the canonical temporal
#' rate; \code{block_size} is user-overridable, and the data-driven
#' \emph{constant} (a spatial analogue of Politis & White 2004, which has no
#' established plug-in form) is available opt-in via \code{block_size = "auto"}
#' (the calibration over \code{g} described below; decision D101). A
#' decorrelating cross-lineage review
#' argued for the \eqn{n^{1/(d+4)}} rate (\eqn{n^{1/6}} at \eqn{d = 2}); that
#' rate governs a different estimand -- the second-order bias / two-sided
#' distribution-function coverage, which gives \eqn{n^{1/5}} at \eqn{d = 1} and
#' so does \emph{not} reduce to the variance default's \eqn{n^{1/3}} -- and is
#' recorded here as a registered dissent rather than adopted.
#'
#' \strong{Resampling.} \code{"tiled"} samples non-empty cells with replacement
#' and truncates to \code{n} (negative bias \eqn{O(1/n)}, negligible);
#' \code{random_origin} draws a fresh sub-cell grid shift per replicate.
#' \code{"moving"} draws overlapping square blocks anchored to cover a sampled
#' observation (never empty).
#'
#' \strong{Guards.} Collinear coordinates (zero range on an axis) abort. If all
#' locations fall in one cell, a warning is emitted and the bootstrap SE
#' collapses toward zero. Coordinates are Euclidean (project lon/lat first).
#'
#' \strong{Data-driven block size (decision D101).} With \code{block_size =
#' "auto"} the cells-per-axis \eqn{g} is chosen by a calibration over a grid of
#' \eqn{g}, because Politis & White (2004) has \strong{no} established spatial
#' plug-in. For each candidate \eqn{g}, \eqn{B_0} cheap (no-refit) spatial block
#' resamples give the bootstrap variance \eqn{V(g)} of the design-weighted
#' residual functionals \eqn{(1/n)\,[1,\tilde{gx},\tilde{gy}]^\top z} (the
#' influence directions of the coefficient, so their MSE-optimal \eqn{g} matches
#' the coefficient's, not merely the residual mean's); \eqn{g} is then chosen to
#' minimise an empirical mean-squared error, the squared bias (anchored at the
#' \emph{largest} blocks, which are the least biased because the
#' dependence-breaking bias grows like \eqn{g/\sqrt n}) plus a leave-one-out
#' jackknife variance, with the \eqn{n^{1/4}} rate as the fallback. \eqn{B_0} is
#' a declared calibration constant. A single isotropic \eqn{g} is used; strongly
#' anisotropic residual dependence is a documented limitation (the minimal fix,
#' two independent coordinate-wise calibrations, is deferred). The bias anchor
#' is the corrected form of a cross-lineage proposal that anchored at the
#' smallest blocks, which would have biased the selector toward anticonservative
#' standard errors.
#'
#' @section Honest scope:
#' The bootstrap delivers robust variance, not better point estimates, and is
#' valid for weak / short-range spatial dependence relative to the block size;
#' it does not rescue strong long-range dependence. gdpar does not model the
#' dependence; here it only makes its inference robust to it.
#'
#' @section Empirical-Bayes and full-Bayes paths:
#' Like its temporal sibling, this function accepts both a scalar
#' \code{gdpar_eb_fit} and a scalar \code{gdpar_fit} (decision D102); see the
#' identically named section of \code{\link{gdpar_dependence_robust}} for the
#' full-Bayes point-estimate / model-SE convention (posterior mean / SD), the
#' \code{se_ratio} interpretation and the full-Bayes cost / Monte-Carlo caveats.
#'
#' @section Dependencies:
#' Uses \pkg{cmdstanr} for the refits and \pkg{posterior} to extract the
#' coefficient estimates (Empirical-Bayes or full-Bayes). The spatial weights
#' and blocks are hand-rolled in base R (no \pkg{spdep} / \pkg{sf} /
#' \pkg{geosphere}).
#'
#' @references
#' Liang, K.-Y. & Zeger, S. L. (1986). Longitudinal data analysis using
#' generalized linear models. \emph{Biometrika} 73(1), 13-22.
#'
#' Kuensch, H. R. (1989). The jackknife and the bootstrap for general
#' stationary observations. \emph{Annals of Statistics} 17(3), 1217-1241.
#'
#' Hall, P., Horowitz, J. L. & Jing, B.-Y. (1995). On blocking rules for the
#' bootstrap with dependent data. \emph{Biometrika} 82(3), 561-574.
#'
#' Politis, D. N. & Romano, J. P. (1992). A circular block resampling procedure
#' for stationary data. In \emph{Exploring the Limits of Bootstrap}, 263-270.
#' Wiley, New York.
#'
#' Lahiri, S. N. (2003). \emph{Resampling Methods for Dependent Data}. Springer,
#' New York.
#'
#' Nordman, D. J. & Lahiri, S. N. (2004). On optimal spatial subsample size for
#' variance estimation. \emph{Annals of Statistics} 32(5), 1981-2027.
#'
#' Politis, D. N. & White, H. (2004). Automatic block-length selection for the
#' dependent bootstrap. \emph{Econometric Reviews} 23(1), 53-70 (with the
#' Patton, Politis & White 2009 correction; the temporal plug-in whose spatial
#' analogue is the \code{"auto"} calibration).
#'
#' @examples
#' \donttest{
#' if (requireNamespace("cmdstanr", quietly = TRUE) &&
#'     requireNamespace("posterior", quietly = TRUE)) {
#'   n <- 100
#'   gx <- runif(n); gy <- runif(n)
#'   x <- rnorm(n)
#'   y <- 1 + 0.5 * x + (gx + gy) + rnorm(n)
#'   df <- data.frame(x = x, y = y)
#'   fit <- gdpar_eb(y ~ x, amm = amm_spec(a = ~ x), data = df,
#'                   chains = 2, iter_warmup = 100, iter_sampling = 100)
#'   # B kept small for a fast example; use B >= 199 in practice.
#'   gdpar_spatial_dependence_robust(fit, data = df, coords = cbind(gx, gy),
#'                                   B = 10, seed = 1, iter_warmup = 100,
#'                                   iter_sampling = 100, chains = 2)
#'   # Data-driven block size (calibration over the cell grid), opt-in:
#'   gdpar_spatial_dependence_robust(fit, data = df, coords = cbind(gx, gy),
#'                                   block_size = "auto", B = 10, seed = 1,
#'                                   iter_warmup = 100, iter_sampling = 100,
#'                                   chains = 2)
#' }
#' }
#'
#' @seealso \code{\link{gdpar_spatial_dependence_diagnostic}},
#'   \code{\link{gdpar_dependence_robust}}, \code{\link{gdpar_eb}}
#' @export
gdpar_spatial_dependence_robust <- function(object, data, coords,
                                            block_size = NULL,
                                            residual_type = c("quantile",
                                                              "response",
                                                              "pearson",
                                                              "deviance"),
                                            randomize_seed = NULL,
                                            scheme = c("tiled", "moving"),
                                            random_origin = TRUE,
                                            B = 199L, level = 0.95, seed = NULL,
                                            iter_warmup = 500L,
                                            iter_sampling = 500L, chains = 2L,
                                            verbose = TRUE, ...) {
  .gdpar_assert_scalar_dep(object, "object")
  assert_data_frame(data, "data")
  residual_type <- match.arg(residual_type)
  scheme <- match.arg(scheme)
  assert_count(B, "B")
  assert_numeric_scalar(level, "level", lower = 0, upper = 1)
  assert_count(iter_warmup, "iter_warmup")
  assert_count(iter_sampling, "iter_sampling")
  assert_count(chains, "chains")
  if (!is.logical(random_origin) || length(random_origin) != 1L) {
    gdpar_abort("Argument 'random_origin' must be a logical scalar.",
                class = "gdpar_input_error")
  }
  if (!is.logical(verbose) || length(verbose) != 1L) {
    gdpar_abort("Argument 'verbose' must be a logical scalar.",
                class = "gdpar_input_error")
  }
  if (!is.null(block_size) && is.character(block_size) &&
      !.gdpar_is_auto(block_size)) {
    gdpar_abort(
      "Argument 'block_size' must be NULL, a positive integer, or \"auto\".",
      class = "gdpar_input_error"
    )
  }

  n <- nrow(data)
  coords <- .gdpar_validate_coords(coords, n)
  ranges <- apply(coords, 2L, function(z) diff(range(z)))
  if (any(ranges < .Machine$double.eps)) {
    gdpar_abort(
      paste("Coordinates are collinear (zero range on an axis); the spatial",
            "block bootstrap requires 2-D variation."),
      class = "gdpar_input_error"
    )
  }
  mins <- apply(coords, 2L, min)

  block_size_method <- "fixed"
  if (is.null(block_size)) {
    block_size <- .gdpar_spatial_default_g(n)
    block_size_method <- "rate"
  } else if (.gdpar_is_auto(block_size)) {
    require_suggested("posterior", "extract posterior draws for residuals")
    resid <- .gdpar_dependence_residuals(object, residual_type, randomize_seed)
    sel <- .gdpar_spatial_block_length_auto(
      coords, resid, scheme = scheme, random_origin = random_origin,
      mins = mins, ranges = ranges, seed = seed
    )
    block_size <- sel$block_size
    block_size_method <- sel$method
    if (verbose) {
      gdpar_inform(paste0("block_size = \"auto\": ", sel$reason),
                   class = "gdpar_optin_message")
    }
  }
  block_size <- as.integer(block_size)
  if (block_size < 1L) {
    gdpar_abort("Argument 'block_size' (cells per axis) must be >= 1.",
                class = "gdpar_input_error")
  }
  g <- block_size

  L0 <- ranges / g
  base_cells <- paste(
    pmin(floor((coords[, 1L] - mins[1L]) / L0[1L]), g - 1L),
    pmin(floor((coords[, 2L] - mins[2L]) / L0[2L]), g - 1L),
    sep = "_"
  )
  n_tiles <- length(unique(base_cells))
  warnings_pre <- character(0L)
  if (n_tiles <= 1L) {
    msg <- paste("All locations fall in a single spatial cell at this",
                 "resolution; the bootstrap SE will collapse toward 0.",
                 "Increase 'block_size' or check the coordinates.")
    warning(msg, call. = FALSE)
    warnings_pre <- msg
  }

  is_fb <- inherits(object, "gdpar_fit") && !inherits(object, "gdpar_eb_fit")
  verbose_msg <- sprintf(
    paste("gdpar_spatial_dependence_robust is an opt-in, computationally",
          "expensive audit: it will perform %d spatial block-bootstrap",
          "refit(s) (%d x %d grid, %s scheme%s)%s. It does not change the point",
          "estimates; it re-estimates their uncertainty to be robust to",
          "spatial dependence."),
    B, g, g, scheme,
    if (scheme == "tiled" && random_origin) ", randomized origin" else "",
    if (is_fb) ", each a full-Bayes refit (full HMC; markedly more expensive than the Empirical-Bayes path)" else ""
  )

  require_suggested("cmdstanr", "refit Path 1 models")
  require_suggested("posterior", "extract posterior draws")

  resample_fun <- function() {
    .gdpar_spatial_block_indices(coords, g, scheme, random_origin, mins, ranges)
  }

  caller_env <- parent.frame()
  core <- .gdpar_dependence_robust_engine(
    object, data, resample_fun, B = B, level = level, seed = seed,
    iter_warmup = iter_warmup, iter_sampling = iter_sampling,
    chains = chains, verbose = verbose, verbose_msg = verbose_msg,
    caller_env = caller_env, ...
  )

  obj <- list(
    table             = core$table,
    block_size        = g,
    block_size_method = block_size_method,
    scheme            = scheme,
    random_origin     = random_origin,
    n_tiles           = as.integer(n_tiles),
    B                 = as.integer(B),
    B_ok              = core$B_ok,
    level             = level,
    seed              = core$seed,
    warnings          = c(warnings_pre, core$warnings),
    refit_diagnostics = core$refit_diagnostics
  )
  class(obj) <- c("gdpar_spatial_dependence_robust", "list")
  obj
}

#' Print method for gdpar_spatial_dependence_robust objects
#'
#' @param x A \code{gdpar_spatial_dependence_robust} object.
#' @param digits Integer; significant digits for the printed table.
#' @param ... Unused; present for S3 generic compatibility.
#' @return Invisibly returns \code{x}.
#' @export
print.gdpar_spatial_dependence_robust <- function(x, digits = 3L, ...) {
  cat("<gdpar_spatial_dependence_robust> ", x$scheme,
      " spatial block bootstrap\n", sep = "")
  bs_method <- x$block_size_method %||% "fixed"
  bs_label <- switch(bs_method,
    auto  = " (auto: data-driven calibration)",
    rate  = " (rate: n^(1/4))",
    fixed = "",
    "")
  cat("  grid            : ", x$block_size, " x ", x$block_size, " cells",
      bs_label,
      if (isTRUE(x$random_origin) && identical(x$scheme, "tiled"))
        " (randomized origin)" else "",
      "   (", x$n_tiles, " non-empty)\n", sep = "")
  cat("  refits          : B = ", x$B, " (", x$B_ok, " ok)\n", sep = "")
  cat("  level           : ", x$level, "\n\n", sep = "")
  tab <- x$table
  num <- vapply(tab, is.numeric, logical(1L))
  tab[num] <- lapply(tab[num], function(col) format(col, digits = digits))
  print(tab, row.names = FALSE)
  cat("\n  se_ratio = robust_se / model_se; > 1 means the model-based\n",
      "  uncertainty understates the spatial-dependence-robust uncertainty.\n",
      sep = "")
  .gdpar_print_refit_diagnostics(x$refit_diagnostics, digits)
  if (length(x$warnings) > 0L) {
    cat("\nNotes (", length(x$warnings), "):\n", sep = "")
    for (w in utils::head(x$warnings, 5L)) cat("  - ", w, "\n", sep = "")
    if (length(x$warnings) > 5L) {
      cat("  ... and ", length(x$warnings) - 5L, " more.\n", sep = "")
    }
  }
  invisible(x)
}

#' Approximate leave-one-out cross-validation for a gdpar fit
#'
#' Computes PSIS-LOO ("Pareto smoothed importance sampling
#' leave-one-out") approximate cross-validation via the \pkg{loo}
#' package, using the per-observation log-likelihood persisted by the
#' Stan model in the generated quantity \code{log_lik}.
#'
#' For univariate fits (\code{p == 1L}) the Stan model emits
#' \code{log_lik} as a vector of length \code{n}; observations are the
#' \code{n} rows of the input data.
#'
#' For multivariate fits (\code{p > 1L}) the Stan model emits
#' \code{log_lik} as an \code{n} by \code{p} matrix with
#' \code{log_lik[i, k]} equal to
#' \eqn{\log p(y_{ik} \mid \theta_i[k])}. Two aggregations are
#' available, selected by \code{aggregation}:
#' \describe{
#'   \item{\code{"subject"} (default)}{
#'     The natural observational unit is the row (subject) of the
#'     input data. Following the coord-wise factorization
#'     \eqn{p(y_i \mid \theta_i) = \prod_{k=1}^p D_k(y_{ik} \mid
#'     \theta_i[k])}, the per-subject log-likelihood is the sum over
#'     coordinates,
#'     \eqn{\log p(y_i \mid \theta_i) = \sum_k \log p(y_{ik} \mid
#'     \theta_i[k])}. This aggregation matches the convention used by
#'     \pkg{brms} multivariate fits with
#'     \code{set_rescor(FALSE)} and yields ELPD values directly
#'     comparable to per-coordinate competitors aggregated identically.
#'   }
#'   \item{\code{"cell"}}{
#'     Each pair \eqn{(i, k)} is treated as an independent
#'     observation, yielding PSIS-LOO over \eqn{n \cdot p} cells. This
#'     is useful for per-coordinate diagnostics (Pareto-k mass
#'     concentrated in a specific coordinate is a signal of a
#'     marginally identified component for that dimension), but
#'     breaks the leave-one-subject-out interpretation of classical
#'     ELPD: the implicit assumption is that the cells are exchangeable
#'     given \eqn{\theta_i}, which is technically true under the
#'     coord-wise factorization but conflates subject-level and
#'     coordinate-level cross-validation. Use for diagnostics, not for
#'     reporting comparable ELPD values across methods.
#'   }
#' }
#'
#' @param fit A \code{gdpar_fit} object.
#' @param aggregation Character scalar, one of \code{"subject"}
#'   (default) or \code{"cell"}. Ignored for univariate fits. See
#'   Details.
#' @param r_eff Optional numeric vector of relative effective sample
#'   sizes per observation. If \code{NULL} (default) it is computed
#'   from the draws via \code{loo::relative_eff()} with the actual
#'   chain identifiers extracted from the fit. Supplying it explicitly
#'   is useful when reusing the result across multiple LOO calls on
#'   the same fit.
#' @param cores Number of cores for the LOO computation; passed to
#'   \code{loo::loo()}. Defaults to 1 (sequential) to avoid
#'   non-determinism.
#' @param ... Additional arguments forwarded to \code{loo::loo()}.
#'
#' @return A \code{loo} object (S3 class \code{"psis_loo"}) with the
#'   ELPD estimate, the elpd_loo standard error, the Pareto-k
#'   diagnostics, and pointwise contributions.
#'
#' @section Status:
#' This function is flagged \code{@keywords experimental}. The
#' aggregation rule (sum over k for multivariate, default
#' \code{"subject"}) is stable and documented above; the signature
#' may gain additional arguments in future versions (e.g.
#' \code{integrand} for non-pointwise predictive quantities).
#' Pareto-k diagnostics with \code{k > 0.7} signal that the PSIS
#' approximation is unreliable for the affected observations; consider
#' \code{loo::loo_moment_match()} or \code{loo::reloo()} as
#' refinements.
#'
#' @seealso \code{\link{gdpar}}, \code{\link[loo]{loo}}.
#'
#' @examples
#' \donttest{
#'   if (requireNamespace("cmdstanr", quietly = TRUE) &&
#'       requireNamespace("loo", quietly = TRUE)) {
#'     n <- 200
#'     set.seed(42)
#'     x1 <- rnorm(n); x2 <- rnorm(n)
#'     y <- 0.5 + 0.7 * x1 - 0.3 * x2 + rnorm(n, sd = 0.5)
#'     dat <- data.frame(y = y, x1 = x1, x2 = x2)
#'     fit <- gdpar(y ~ x1 + x2, data = dat,
#'                  family = gdpar_family("gaussian"),
#'                  chains = 2, iter_warmup = 500,
#'                  iter_sampling = 500, refresh = 0)
#'     lo <- gdpar_loo(fit)
#'     print(lo)
#'   }
#' }
#'
#' @keywords experimental
#' @export
gdpar_loo <- function(fit,
                       aggregation = c("subject", "cell"),
                       r_eff = NULL,
                       cores = 1L,
                       ...) {
  assert_inherits(fit, "gdpar_fit", "fit")
  aggregation <- match.arg(aggregation)
  require_suggested("loo", "compute PSIS-LOO approximate cross-validation")

  draws_arr <- fit$fit$draws(variables = "log_lik", format = "draws_array")
  p <- if (!is.null(fit$amm$p)) as.integer(fit$amm$p) else 1L

  log_lik_mat <- aggregate_log_lik(draws_arr, p, aggregation)

  if (is.null(r_eff)) {
    n_iter <- dim(draws_arr)[1L]
    n_chain <- dim(draws_arr)[2L]
    chain_id <- rep(seq_len(n_chain), each = n_iter)
    r_eff <- loo::relative_eff(exp(log_lik_mat), chain_id = chain_id,
                                cores = cores)
  }

  loo::loo(log_lik_mat, r_eff = r_eff, cores = cores, ...)
}

#' Internal: aggregate log_lik draws to the requested scale
#'
#' Converts a \code{draws_array} of \code{log_lik} into a matrix of
#' shape \code{S x n} (subject aggregation) or \code{S x (n*p)} (cell
#' aggregation), suitable for \code{loo::loo()}.
#'
#' For \code{p == 1L} the input is already a vector \code{log_lik[i]}
#' and both aggregations return the same \code{S x n} matrix.
#'
#' For \code{p > 1L}, the variable names follow the Stan convention
#' \code{"log_lik[i,k]"}; we parse them and either sum over k
#' (subject) or concatenate (cell). The order of columns in the
#' resulting matrix is i-major, k-minor for \code{"cell"} (i.e.,
#' columns are c((1,1), (1,2), ..., (1,p), (2,1), ...)).
#'
#' @keywords internal
#' @noRd
aggregate_log_lik <- function(draws_arr, p, aggregation) {
  vars <- posterior::variables(draws_arr)
  mat <- unclass(posterior::as_draws_matrix(draws_arr))

  if (p == 1L) {
    return(mat)
  }

  m <- regexec("^log_lik\\[(\\d+),(\\d+)\\]$", vars)
  parsed <- regmatches(vars, m)
  bad <- vapply(parsed, function(z) length(z) == 0L, logical(1L))
  if (any(bad)) {
    gdpar_abort(
      sprintf(
        "Unexpected log_lik variable names in multivariate fit: %s.",
        paste(sQuote(utils::head(vars[bad], 5L)), collapse = ", ")
      ),
      class = "gdpar_loo_parse_error",
      data = list(unparsed = vars[bad])
    )
  }
  ij <- do.call(rbind, lapply(parsed, function(z) as.integer(z[2:3])))
  colnames(ij) <- c("i", "k")
  n <- max(ij[, "i"])
  p_in <- max(ij[, "k"])
  if (p_in != p) {
    gdpar_abort(
      sprintf(
        "log_lik dimension mismatch: parsed p = %d from draws; fit declares p = %d.",
        p_in, p
      ),
      class = "gdpar_loo_dim_mismatch",
      data = list(p_parsed = p_in, p_fit = p)
    )
  }

  S <- nrow(mat)
  if (aggregation == "subject") {
    out <- matrix(0, nrow = S, ncol = n)
    for (col_idx in seq_along(vars)) {
      i <- ij[col_idx, "i"]
      out[, i] <- out[, i] + mat[, col_idx]
    }
    return(out)
  }

  out <- matrix(0, nrow = S, ncol = n * p)
  for (col_idx in seq_along(vars)) {
    i <- ij[col_idx, "i"]
    k <- ij[col_idx, "k"]
    flat_col <- (i - 1L) * p + k
    out[, flat_col] <- mat[, col_idx]
  }
  out
}

# Sub-phase 8.3.9 (2026-05-22): G1 deviance / Pearson / response
# residuals, G2 randomized quantile residuals (Dunn-Smyth 1996) and
# G3 posterior predictive machinery for gdpar_fit objects.
#
# DHARMa is consumed as Suggests with detect-and-use (decision E1.A
# of 8.3.9 scoping). When DHARMa is available the user can build a
# DHARMa object via .gdpar_build_dharma_object(); otherwise the
# Bayesian quantile residuals implemented here (averaged across
# posterior draws of y_pred) cover the same metodological role.
#
# Bayesplot is Suggests; pp_check.gdpar_fit dispatches off the
# bayesplot generic, so users that do not load bayesplot cannot call
# pp_check directly but can still call posterior_predict and feed the
# result to ggplot or to base R graphics.

# ---- Internal helpers --------------------------------------------------

#' Internal: detect the fit "path" classification for residuals dispatch.
#'
#' Returns one of "scalar" (K = 1 + p = 1), "multivariate" (K = 1 +
#' p > 1) or "K_individual" (K > 1 + p = 1). Mirrors the three paths
#' in gdpar.R (assemble_stan_data, .assemble_stan_data_multi,
#' .assemble_stan_data_K) and the three Stan templates.
#'
#' @keywords internal
#' @noRd
.gdpar_fit_path_class <- function(object) {
  K_val <- object[["K"]]
  p_val <- object[["p"]]
  if (!is.null(K_val) && K_val > 1L) {
    return("K_individual")
  }
  if (!is.null(p_val) && p_val > 1L) {
    return("multivariate")
  }
  "scalar"
}

#' Internal: integer-valued family stan_ids (Poisson, NB, Bernoulli,
#' mixtures and Hurdle).
#'
#' @keywords internal
#' @noRd
.gdpar_integer_family_stan_ids <- function() {
  c(2L, 3L, 4L, 10L, 11L, 12L, 13L)
}

#' Internal: extract observed y from a gdpar_fit object via stan_data.
#'
#' Returns a numeric vector for the scalar / K-individual paths and
#' a numeric matrix \code{n x p} for the multivariate path. The
#' choice between \code{stan_data$y_real} and \code{stan_data$y_int}
#' is driven by the family's \code{stan_id} (integer families use
#' \code{y_int}; real families use \code{y_real}). For the
#' multivariate path the matrix is returned as-is.
#'
#' @keywords internal
#' @noRd
.gdpar_get_y_obs <- function(object) {
  sd <- object$stan_data
  path <- .gdpar_fit_path_class(object)
  if (path == "multivariate") {
    if (!is.null(sd$y_real) && is.matrix(sd$y_real) && nrow(sd$y_real) > 0L) {
      return(sd$y_real)
    }
    if (!is.null(sd$y_int) && is.matrix(sd$y_int) && nrow(sd$y_int) > 0L) {
      return(sd$y_int)
    }
    gdpar_abort(
      "Internal error: multivariate fit lacks both y_real and y_int in stan_data.",
      class = "gdpar_internal_error"
    )
  }
  stan_id <- if (path == "K_individual") {
    object$family$stan_id
  } else {
    object$family$stan_id
  }
  if (is.null(stan_id)) {
    gdpar_abort(
      "Internal error: family lacks stan_id; cannot determine y_real vs y_int.",
      class = "gdpar_internal_error"
    )
  }
  if (stan_id %in% .gdpar_integer_family_stan_ids()) {
    return(as.numeric(sd$y_int))
  }
  as.numeric(sd$y_real)
}

#' Internal: extract posterior draws of y_pred as a matrix / array.
#'
#' Returns a numeric matrix of shape \code{S x n} for the scalar
#' and K-individual paths, and a numeric array of shape
#' \code{S x n x p} for the multivariate path.
#'
#' @keywords internal
#' @noRd
.gdpar_get_y_pred_draws <- function(object) {
  require_suggested("posterior", "extract posterior y_pred draws")
  draws <- object$fit$draws(variables = "y_pred", format = "draws_matrix")
  mat <- unclass(draws)
  path <- .gdpar_fit_path_class(object)
  var_names <- colnames(mat)
  if (path == "multivariate") {
    matches <- regmatches(
      var_names,
      regexec("^y_pred\\[(\\d+),(\\d+)\\]$", var_names)
    )
    idx_i <- as.integer(vapply(matches, function(m) m[2L], character(1L)))
    idx_k <- as.integer(vapply(matches, function(m) m[3L], character(1L)))
    if (any(is.na(idx_i)) || any(is.na(idx_k))) {
      gdpar_abort(
        "Internal error: failed to parse multivariate y_pred variable names.",
        class = "gdpar_internal_error"
      )
    }
    n <- max(idx_i)
    p <- max(idx_k)
    S <- nrow(mat)
    arr <- array(NA_real_, dim = c(S, n, p))
    for (col in seq_len(ncol(mat))) {
      arr[, idx_i[col], idx_k[col]] <- mat[, col]
    }
    return(arr)
  }
  matches <- regmatches(
    var_names,
    regexec("^y_pred\\[(\\d+)\\]$", var_names)
  )
  idx_i <- as.integer(vapply(matches, function(m) m[2L], character(1L)))
  if (any(is.na(idx_i))) {
    gdpar_abort(
      "Internal error: failed to parse scalar y_pred variable names.",
      class = "gdpar_internal_error"
    )
  }
  out <- matrix(NA_real_, nrow = nrow(mat), ncol = max(idx_i))
  for (col in seq_len(ncol(mat))) {
    out[, idx_i[col]] <- mat[, col]
  }
  out
}

#' Internal: family-name classifier used by quantile and deviance helpers.
#'
#' For the K-individual path the relevant family is the slot-1
#' (location) family; for the multivariate path it is one of the
#' per-coord families. The helper returns a single character.
#'
#' @keywords internal
#' @noRd
.gdpar_family_name_for_residuals <- function(object, coord = NULL) {
  path <- .gdpar_fit_path_class(object)
  if (path == "multivariate") {
    if (is.null(coord)) {
      gdpar_abort(
        "Internal error: multivariate path requires explicit coord index.",
        class = "gdpar_internal_error"
      )
    }
    fam <- object$family$families[[coord]]
    return(fam$name)
  }
  object$family$name
}

#' Internal: detect discrete families (count or Bernoulli).
#'
#' Used to switch on randomization in quantile residuals (Dunn-Smyth).
#'
#' @keywords internal
#' @noRd
.gdpar_family_is_discrete <- function(name) {
  name %in% c("poisson", "neg_binomial_2", "bernoulli",
              "zip", "zinb", "hurdle_poisson", "hurdle_neg_binomial_2")
}

#' Internal: family-specific deviance residual on (y, mu) pairs.
#'
#' Uses the canonical deviance contribution for the gaussian, poisson,
#' neg_binomial_2 (approximated by Poisson deviance with shape inflation
#' factor, conservative), bernoulli, beta, gamma, lognormal, student_t,
#' tweedie families. For mixtures and Hurdle the helper falls back to
#' a sign(y - mu) * sqrt(2 * |y - mu|) approximation (a Pearson-like
#' surrogate). The contribution is then signed by sign(y - mu).
#'
#' @keywords internal
#' @noRd
.gdpar_deviance_residual <- function(y, mu, family_name, dispersion = NULL) {
  eps <- .Machine$double.eps^(1 / 3)
  mu_safe <- pmax(mu, eps)
  if (family_name == "gaussian") {
    sigma <- if (is.null(dispersion)) 1 else max(dispersion, eps)
    return(sign(y - mu) * abs(y - mu) / sigma)
  }
  if (family_name == "poisson" || family_name == "zip" ||
      family_name == "hurdle_poisson") {
    term1 <- ifelse(y > 0, y * log(y / mu_safe), 0)
    d_i   <- 2 * (term1 - (y - mu))
    return(sign(y - mu) * sqrt(pmax(d_i, 0)))
  }
  if (family_name == "neg_binomial_2" || family_name == "zinb" ||
      family_name == "hurdle_neg_binomial_2") {
    phi <- if (is.null(dispersion)) 1 else max(dispersion, eps)
    term1 <- ifelse(y > 0, y * log(y / mu_safe), 0)
    term2 <- (y + phi) * log((mu_safe + phi) / (y + phi))
    d_i <- 2 * (term1 + term2)
    return(sign(y - mu) * sqrt(pmax(d_i, 0)))
  }
  if (family_name == "bernoulli") {
    mu_clip <- pmin(pmax(mu_safe, eps), 1 - eps)
    term <- y * log(mu_clip) + (1 - y) * log(1 - mu_clip)
    d_i <- -2 * term
    return(sign(y - mu) * sqrt(pmax(d_i, 0)))
  }
  if (family_name == "beta") {
    mu_clip <- pmin(pmax(mu_safe, eps), 1 - eps)
    y_clip <- pmin(pmax(y, eps), 1 - eps)
    term <- y_clip * log(y_clip / mu_clip) +
      (1 - y_clip) * log((1 - y_clip) / (1 - mu_clip))
    d_i <- 2 * term
    return(sign(y - mu) * sqrt(pmax(d_i, 0)))
  }
  if (family_name == "gamma") {
    y_safe <- pmax(y, eps)
    d_i <- 2 * (-log(y_safe / mu_safe) + (y_safe - mu_safe) / mu_safe)
    return(sign(y - mu) * sqrt(pmax(d_i, 0)))
  }
  if (family_name == "lognormal_loc_scale") {
    y_safe <- pmax(y, eps)
    return(sign(y - mu) * abs(log(y_safe) - log(mu_safe)))
  }
  if (family_name == "student_t") {
    sigma <- if (is.null(dispersion)) 1 else max(dispersion, eps)
    return(sign(y - mu) * abs(y - mu) / sigma)
  }
  if (family_name == "tweedie") {
    return(sign(y - mu) * sqrt(pmax(2 * abs(y - mu), 0)))
  }
  sign(y - mu) * sqrt(pmax(2 * abs(y - mu), 0))
}

#' Internal: Bayesian randomized quantile residuals (Dunn-Smyth 1996).
#'
#' For each observation \eqn{y_i} the helper averages
#' \eqn{P(Y \le y_i \mid \theta_i^{(s)})} over posterior draws s of
#' \eqn{\theta_i}, using the empirical ECDF of the y_pred draws as a
#' nonparametric estimator of the family-conditional CDF. Discrete
#' families add a uniform jitter on the equality mass to avoid
#' degenerate u_i, then map via qnorm.
#'
#' @keywords internal
#' @noRd
.gdpar_quantile_residuals_bayesian <- function(y_obs, y_pred_mat,
                                               family_name,
                                               randomize_seed = NULL) {
  if (!is.null(randomize_seed)) {
    if (exists(".Random.seed", envir = globalenv())) {
      old_seed <- get(".Random.seed", envir = globalenv())
      on.exit(assign(".Random.seed", old_seed, envir = globalenv()),
              add = TRUE)
    }
    set.seed(randomize_seed)
  }
  n <- length(y_obs)
  S <- nrow(y_pred_mat)
  discrete <- .gdpar_family_is_discrete(family_name)
  u <- numeric(n)
  for (i in seq_len(n)) {
    yp <- y_pred_mat[, i]
    if (discrete) {
      p_lt <- mean(yp < y_obs[i])
      p_le <- mean(yp <= y_obs[i])
      v <- stats::runif(1L)
      u[i] <- p_lt + v * (p_le - p_lt)
    } else {
      u[i] <- mean(yp <= y_obs[i])
    }
  }
  eps <- .Machine$double.eps^(1 / 3)
  u <- pmin(pmax(u, eps), 1 - eps)
  stats::qnorm(u)
}

#' Internal: detect DHARMa availability (E1.A Suggests with detect-and-use).
#'
#' @keywords internal
#' @noRd
.gdpar_has_dharma <- function() {
  requireNamespace("DHARMa", quietly = TRUE)
}

#' Internal: build a DHARMa object from a fitted gdpar model.
#'
#' Lifts the y_obs and y_pred draws into a DHARMa simulationOutput
#' object so the user can call DHARMa::testResiduals,
#' DHARMa::testZeroInflation, DHARMa::plotResiduals on the result.
#' Returns NULL if DHARMa is not installed (the user-facing path
#' surfaces a clear error in that case).
#'
#' @keywords internal
#' @noRd
.gdpar_build_dharma_object <- function(object, coord = NULL) {
  if (!.gdpar_has_dharma()) {
    gdpar_abort(
      paste(
        "DHARMa is not installed. Install via install.packages('DHARMa') to",
        "use gdpar_dharma_object(); or call residuals(., type = 'quantile')",
        "for the built-in Bayesian quantile residuals fallback."
      ),
      class = "gdpar_input_error"
    )
  }
  path <- .gdpar_fit_path_class(object)
  if (path == "multivariate" && is.null(coord)) {
    gdpar_abort(
      paste(
        "For multivariate fits (p > 1), pass coord = k to",
        "gdpar_dharma_object() to build a DHARMa object for coordinate k."
      ),
      class = "gdpar_input_error"
    )
  }
  y_obs <- .gdpar_get_y_obs(object)
  y_pred <- .gdpar_get_y_pred_draws(object)
  family_name <- .gdpar_family_name_for_residuals(object, coord = coord)
  integer_family <- .gdpar_family_is_discrete(family_name)
  if (path == "multivariate") {
    y_pred_mat <- y_pred[, , coord]
    y_obs_vec <- as.numeric(y_obs[, coord])
  } else {
    y_pred_mat <- y_pred
    y_obs_vec <- as.numeric(y_obs)
  }
  fitted_pred <- colMeans(y_pred_mat)
  DHARMa::createDHARMa(
    simulatedResponse = t(y_pred_mat),
    observedResponse  = y_obs_vec,
    fittedPredictedResponse = fitted_pred,
    integerResponse   = integer_family
  )
}

# ---- Exported S3 methods ----------------------------------------------

#' Residuals for a fitted gdpar model
#'
#' Posterior-predictive residuals for a \code{gdpar_fit} object,
#' covering response, Pearson, deviance and randomized quantile
#' (Dunn-Smyth 1996) residual types. The residuals are computed from
#' the \code{y_pred} draws emitted by the Stan templates and are
#' returned as a numeric vector (scalar and K-individual paths) or as
#' a numeric matrix \code{n} by \code{p} (multivariate path).
#'
#' @param object An object of class \code{gdpar_fit}.
#' @param type Character scalar in \code{"response"} (raw residual
#'   \eqn{y_i - \bar y_{i, \text{pred}}}), \code{"pearson"} (response
#'   divided by the posterior-predictive standard deviation per
#'   observation), \code{"deviance"} (family-specific deviance
#'   contribution; falls back to a Pearson-like surrogate for
#'   mixtures and Hurdle), or \code{"quantile"} (Bayesian randomized
#'   quantile residuals; the canonical residual under the model).
#'   Default is \code{"quantile"}.
#' @param coord Integer scalar between 1 and p; only used for
#'   multivariate fits. When \code{NULL} (default) the function
#'   returns a matrix \code{n} by \code{p} with one column per
#'   coordinate.
#' @param randomize_seed Optional integer scalar to set the RNG seed
#'   used by the randomized quantile residual under discrete families.
#'   When \code{NULL}, the residuals use the global RNG state.
#' @param ... Unused; present for S3 generic compatibility.
#'
#' @return A numeric vector for scalar / K-individual paths; a
#'   numeric matrix \code{n} by \code{p} (or numeric vector of length
#'   \code{n} if \code{coord} is specified) for multivariate paths.
#'
#' @section Methodology:
#' The Bayesian randomized quantile residual averages the
#' nonparametric ECDF of \code{y_pred} draws at \eqn{y_i} across
#' posterior draws, adding a uniform jitter on the equality mass when
#' the family is discrete (Dunn and Smyth 1996). Under a correctly
#' specified model the residuals are marginally \eqn{\mathcal N(0,
#' 1)} regardless of the family. Deviance and Pearson residuals are
#' provided for parity with frequentist diagnostics; their
#' distribution under the null model is approximate for non-Gaussian
#' families. For mixtures (ZIP/ZINB) and Hurdle families the deviance
#' residual is approximated by a Pearson-like surrogate; the
#' quantile residual remains canonical and is recommended.
#'
#' @references
#' Dunn, P. K. and Smyth, G. K. (1996). Randomized Quantile
#' Residuals. \emph{Journal of Computational and Graphical
#' Statistics}, 5(3), 236-244.
#'
#' @export
residuals.gdpar_fit <- function(object,
                                type = c("quantile", "response",
                                         "pearson", "deviance"),
                                coord = NULL,
                                randomize_seed = NULL,
                                ...) {
  assert_inherits(object, "gdpar_fit", "object")
  type <- match.arg(type)
  path <- .gdpar_fit_path_class(object)
  y_obs <- .gdpar_get_y_obs(object)
  y_pred <- .gdpar_get_y_pred_draws(object)
  if (path == "multivariate") {
    if (!is.null(coord)) {
      coord <- as.integer(coord)
      p <- object[["p"]]
      if (coord < 1L || coord > p) {
        gdpar_abort(
          sprintf(
            "Argument 'coord' must be an integer in [1, %d]; received %s.",
            p, coord
          ),
          class = "gdpar_input_error"
        )
      }
      y_pred_mat <- y_pred[, , coord]
      y_obs_vec <- as.numeric(y_obs[, coord])
      family_name <- .gdpar_family_name_for_residuals(object, coord = coord)
      return(.gdpar_residuals_dispatch(
        y_obs_vec, y_pred_mat, type, family_name, randomize_seed
      ))
    }
    p <- object[["p"]]
    out <- matrix(NA_real_, nrow = nrow(y_obs), ncol = p)
    colnames(out) <- paste0("dim_", seq_len(p))
    for (k in seq_len(p)) {
      family_name <- .gdpar_family_name_for_residuals(object, coord = k)
      out[, k] <- .gdpar_residuals_dispatch(
        as.numeric(y_obs[, k]), y_pred[, , k], type, family_name,
        randomize_seed
      )
    }
    return(out)
  }
  family_name <- .gdpar_family_name_for_residuals(object)
  .gdpar_residuals_dispatch(
    as.numeric(y_obs), y_pred, type, family_name, randomize_seed
  )
}

#' Internal: type-dispatcher for a (y_obs, y_pred_mat) pair.
#'
#' @keywords internal
#' @noRd
.gdpar_residuals_dispatch <- function(y_obs, y_pred_mat, type, family_name,
                                      randomize_seed) {
  if (type == "response") {
    return(y_obs - colMeans(y_pred_mat))
  }
  if (type == "pearson") {
    mu_hat <- colMeans(y_pred_mat)
    sd_hat <- apply(y_pred_mat, 2L, stats::sd)
    eps <- .Machine$double.eps^(1 / 3)
    sd_safe <- pmax(sd_hat, eps)
    return((y_obs - mu_hat) / sd_safe)
  }
  if (type == "deviance") {
    mu_hat <- colMeans(y_pred_mat)
    return(.gdpar_deviance_residual(y_obs, mu_hat, family_name))
  }
  .gdpar_quantile_residuals_bayesian(
    y_obs, y_pred_mat, family_name, randomize_seed
  )
}

#' Posterior predictive draws for a fitted gdpar model
#'
#' Extracts the posterior predictive draws \code{y_pred} emitted by
#' the Stan templates. Returns a matrix of dimensions S by n (scalar
#' and K-individual paths) or an array of dimensions S by n by p
#' (multivariate path).
#'
#' @param object An object of class \code{gdpar_fit}.
#' @param ndraws Optional integer scalar. When given, subsamples the
#'   first \code{ndraws} posterior draws. When NULL (default), all
#'   draws are returned.
#' @param ... Unused; present for S3 generic compatibility.
#'
#' @return Numeric matrix S by n (scalar/K-individual paths) or
#'   numeric array S by n by p (multivariate path).
#'
#' @details
#' This function is independent of the \pkg{rstantools} generic; it
#' is exported under the name \code{gdpar_posterior_predict} to avoid
#' a hard dependency on \pkg{rstantools}. Users that load
#' \pkg{rstantools} or \pkg{brms} will see \code{posterior_predict()}
#' route to \code{gdpar_posterior_predict()} via the S3 method
#' registered in NAMESPACE.
#'
#' @export
gdpar_posterior_predict <- function(object, ndraws = NULL, ...) {
  assert_inherits(object, "gdpar_fit", "object")
  out <- .gdpar_get_y_pred_draws(object)
  if (is.null(ndraws)) {
    return(out)
  }
  ndraws <- as.integer(ndraws)
  if (length(ndraws) != 1L || is.na(ndraws) || ndraws < 1L) {
    gdpar_abort(
      "Argument 'ndraws' must be a positive integer scalar or NULL.",
      class = "gdpar_input_error"
    )
  }
  if (is.matrix(out)) {
    n_total <- nrow(out)
    keep <- seq_len(min(ndraws, n_total))
    return(out[keep, , drop = FALSE])
  }
  n_total <- dim(out)[1L]
  keep <- seq_len(min(ndraws, n_total))
  out[keep, , , drop = FALSE]
}

#' Build a DHARMa simulation object from a fitted gdpar model
#'
#' If \pkg{DHARMa} is available (Suggests), returns a
#' \code{DHARMa::DHARMa} object built from the posterior predictive
#' draws of \code{y_pred} and the observed response \code{y}. The
#' user can then call \code{DHARMa::testResiduals()},
#' \code{DHARMa::testZeroInflation()}, \code{DHARMa::plotResiduals()}
#' and other DHARMa tests on the returned object.
#'
#' @param object An object of class \code{gdpar_fit}.
#' @param coord Integer scalar between 1 and p; required for
#'   multivariate fits. Ignored for scalar / K-individual paths.
#'
#' @return A \code{DHARMa} simulation object.
#'
#' @section Dependencies:
#' Requires \pkg{DHARMa} (Suggests). If the package is not installed
#' the function raises \code{gdpar_input_error} pointing to
#' \code{residuals(., type = "quantile")} as the built-in fallback.
#'
#' @export
gdpar_dharma_object <- function(object, coord = NULL) {
  assert_inherits(object, "gdpar_fit", "object")
  .gdpar_build_dharma_object(object, coord = coord)
}

#' Posterior predictive check for a fitted gdpar model
#'
#' S3 method off the \pkg{bayesplot} generic. Forwards
#' \code{y_pred} draws and the observed response \code{y} to one of
#' the \code{bayesplot::ppc_*} family of functions selected by
#' \code{type}.
#'
#' @param object An object of class \code{gdpar_fit}.
#' @param type Character scalar selecting the bayesplot ppc plot:
#'   \code{"dens_overlay"} (default), \code{"hist"},
#'   \code{"ecdf_overlay"}, \code{"stat"} or \code{"intervals"}.
#' @param coord Integer scalar between 1 and p; required for
#'   multivariate fits.
#' @param ndraws Integer scalar; subsamples the first \code{ndraws}
#'   posterior draws before plotting. Defaults to 50.
#' @param ... Additional arguments forwarded to the bayesplot
#'   function.
#'
#' @return A \code{ggplot} object produced by bayesplot.
#'
#' @section Dependencies:
#' Requires \pkg{bayesplot} (Suggests). If not installed, raises a
#' clear error.
#'
#' @exportS3Method bayesplot::pp_check
pp_check.gdpar_fit <- function(object,
                               type = c("dens_overlay", "hist",
                                        "ecdf_overlay", "stat",
                                        "intervals"),
                               coord = NULL,
                               ndraws = 50L,
                               ...) {
  assert_inherits(object, "gdpar_fit", "object")
  type <- match.arg(type)
  require_suggested("bayesplot", "produce a posterior predictive check plot")
  path <- .gdpar_fit_path_class(object)
  y_obs <- .gdpar_get_y_obs(object)
  y_pred <- .gdpar_get_y_pred_draws(object)
  if (path == "multivariate") {
    if (is.null(coord)) {
      gdpar_abort(
        paste(
          "For multivariate fits (p > 1), pass coord = k to",
          "pp_check() to plot a posterior predictive check for",
          "coordinate k."
        ),
        class = "gdpar_input_error"
      )
    }
    coord <- as.integer(coord)
    p <- object[["p"]]
    if (coord < 1L || coord > p) {
      gdpar_abort(
        sprintf(
          "Argument 'coord' must be an integer in [1, %d]; received %s.",
          p, coord
        ),
        class = "gdpar_input_error"
      )
    }
    y_obs_vec <- as.numeric(y_obs[, coord])
    yrep_mat <- y_pred[, , coord]
  } else {
    y_obs_vec <- as.numeric(y_obs)
    yrep_mat <- y_pred
  }
  ndraws <- as.integer(ndraws)
  if (length(ndraws) != 1L || is.na(ndraws) || ndraws < 1L) {
    gdpar_abort(
      "Argument 'ndraws' must be a positive integer scalar.",
      class = "gdpar_input_error"
    )
  }
  keep <- seq_len(min(ndraws, nrow(yrep_mat)))
  yrep_mat <- yrep_mat[keep, , drop = FALSE]
  fn <- switch(type,
    dens_overlay  = bayesplot::ppc_dens_overlay,
    hist          = bayesplot::ppc_hist,
    ecdf_overlay  = bayesplot::ppc_ecdf_overlay,
    stat          = bayesplot::ppc_stat,
    intervals     = bayesplot::ppc_intervals
  )
  fn(y = y_obs_vec, yrep = yrep_mat, ...)
}

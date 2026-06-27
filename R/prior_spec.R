#' Internal: canonical default prior per parameter kind
#'
#' Block 8 Session 1 decision (2b-iii) canonizes prior specification by
#' kind of parameter, in distinction from the legacy specification by
#' hyperparameter name (\code{theta_ref}, \code{sigma_a}, etc., which
#' remain backward-compatible for K = 1). The kind classifies the role
#' of an individual-scope parameter under the multi-parametric
#' extension: \code{mu} (location on the link scale), \code{log_sigma}
#' (Gaussian residual scale on log scale), \code{log_phi} (NB
#' dispersion on log scale), \code{logit_p} (Bernoulli probability on
#' logit scale), \code{log_shape} (Gamma shape on log scale),
#' \code{log_nu} (Student-t degrees of freedom on log scale),
#' \code{logit_pi} (mixture probability for ZIP/Hurdle on logit scale),
#' \code{power_p} (Tweedie power parameter on the natural (1, 2)
#' interval).
#'
#' Defaults are weakly informative on the canonical scale of each
#' kind. They are consumed by \code{generate_stan_code*()} when a
#' \code{param_spec} carries a \code{prior_canonical_kind} that has
#' not been overridden by the user via the \code{priors_by_kind}
#' argument of \code{\link{gdpar_prior}}.
#'
#' @param kind Character scalar identifying the parameter kind.
#' @return Character scalar with the Stan-syntax prior expression.
#' @keywords internal
#' @noRd
.gdpar_canonical_prior_for_kind <- function(kind) {
  switch(kind,
    mu        = "normal(0, 2.5)",
    log_sigma = "normal(0, 1)",
    log_phi   = "normal(0, 1)",
    logit_p   = "normal(0, 2.5)",
    log_shape = "normal(0, 1)",
    log_nu    = "normal(log(10), 1)",
    logit_pi  = "normal(0, 2.5)",
    power_p   = "uniform(1.01, 1.99)",
    gdpar_abort(
      sprintf("Internal error: no canonical prior for kind '%s'.", kind),
      class = "gdpar_internal_error"
    )
  )
}

#' Internal: canonical prior variance per parameter kind
#'
#' Returns the variance of the canonical default prior for a given
#' \code{prior_canonical_kind}, used by the post-fit information ratio
#' diagnostic (D-B1) introduced in sub-phase 8.3.4 of Block 8. The
#' contraction metric
#' \eqn{1 - \mathrm{var}_{post} / \mathrm{var}_{prior}} compares each
#' K-individual slot's posterior variance against the canonical prior
#' variance returned here.
#'
#' Values map to the priors emitted by
#' \code{.gdpar_canonical_prior_for_kind()}:
#' \itemize{
#'   \item \code{mu}, \code{logit_p}, \code{logit_pi}: \code{normal(0, 2.5)} -> 6.25
#'   \item \code{log_sigma}, \code{log_phi}, \code{log_shape}: \code{normal(0, 1)} -> 1.0
#'   \item \code{log_nu}: \code{normal(log(10), 1)} -> 1.0 (location-shifted, scale identical)
#'   \item \code{power_p}: \code{uniform(1.01, 1.99)} -> 0.08 ((b-a)^2/12)
#' }
#'
#' @param kind Character scalar identifying the parameter kind.
#' @return Numeric scalar with the canonical prior variance.
#' @keywords internal
#' @noRd
.gdpar_canonical_prior_variance <- function(kind) {
  switch(kind,
    mu        = 6.25,
    log_sigma = 1.0,
    log_phi   = 1.0,
    logit_p   = 6.25,
    log_shape = 1.0,
    log_nu    = 1.0,
    logit_pi  = 6.25,
    power_p   = (1.99 - 1.01)^2 / 12,
    NA_real_
  )
}

#' Internal: known canonical kinds
#'
#' Returns the character vector of all parameter kinds for which a
#' canonical prior is registered. Used by the validator of
#' \code{priors_by_kind} in \code{\link{gdpar_prior}}.
#'
#' @return Character vector.
#' @keywords internal
#' @noRd
.gdpar_known_canonical_kinds <- function() {
  c("mu", "log_sigma", "log_phi", "logit_p",
    "log_shape", "log_nu", "logit_pi", "power_p")
}

#' Specify the priors for the AMM hierarchical Bayesian model
#'
#' Build a prior specification consumed by \code{\link{gdpar}} when
#' \code{path = "bayes"}. All defaults are weakly informative on the
#' linear-predictor scale of the family, calibrated for covariates
#' standardized to unit variance and outcomes on the natural scale of
#' the family. Each component can be overridden individually.
#'
#' @param theta_ref Character scalar with the prior on the population
#'   reference theta_ref on the linear-predictor scale, in Stan
#'   syntax. Default is \code{"normal(0, 2.5)"}. When grouping is
#'   active (the \code{group} argument of \code{\link{gdpar}} is not
#'   \code{NULL}), this prior is applied to the hyperparameter
#'   \code{mu_theta_ref} (the population mean of \code{theta_ref[g]});
#'   when grouping is inactive it is applied directly to the single
#'   \code{theta_ref}, preserving the Block 6 semantics.
#' @param sigma_theta_ref Character scalar with the prior on the
#'   hierarchical scale of \code{theta_ref[g]} across groups, in Stan
#'   syntax. Default is \code{"student_t(3, 0, 1)"} (truncated to
#'   positive values). Used only when grouping is active; ignored in
#'   the single-anchor regime.
#' @param sigma_a Character scalar with the prior on the hierarchical
#'   scale of the additive component coefficients, in Stan syntax.
#'   Default is \code{"student_t(3, 0, 1)"} (truncated to positive
#'   values inside the Stan code).
#' @param sigma_b Character scalar with the prior on the hierarchical
#'   scale of the multiplicative contribution to the linear predictor,
#'   in Stan syntax. Default is \code{"student_t(3, 0, 1)"} (truncated
#'   to positive values). Internally the model samples
#'   \code{c_b = theta_ref * b_coef} (the linearly identified quantity
#'   in \code{eta = theta_ref + Z_a * a + Z_b * c_b + ruido}); the
#'   prior \code{sigma_b} is applied as the scale of \code{c_b}, which
#'   coincides with the scale of \code{b_coef} only when
#'   \code{|theta_ref|} is close to 1. For weakly informative priors on
#'   covariates standardized to unit variance and outcomes on the
#'   linear-predictor scale, \code{student_t(3, 0, 1)} remains a
#'   reasonable default. See the Methodological notes section.
#' @param sigma_W Character scalar with the prior on the hierarchical
#'   scale of the modulating component coefficients, in Stan syntax.
#'   Default is \code{"student_t(3, 0, 1)"} (truncated to positive
#'   values).
#' @param sigma_y Character scalar with the prior on the residual
#'   standard deviation for Gaussian families, in Stan syntax. Default
#'   is \code{"student_t(3, 0, 2.5)"} (truncated to positive values).
#' @param phi Character scalar with the prior on the negative-binomial
#'   dispersion phi (Stan parametrization neg_binomial_2: variance =
#'   mu + mu^2 / phi). Default is \code{"gamma(2, 0.1)"}.
#' @param priors_by_kind Optional named list of Stan-syntax prior
#'   strings, indexed by \code{prior_canonical_kind} (decision 2b-iii of
#'   Block 8 Session 1). Each entry overrides the canonical default for
#'   that kind, used by the codegen when K > 1 (multi-parametric
#'   extension) introduces individual-scope parameters beyond the
#'   primary one. Recognized kinds: \code{mu}, \code{log_sigma},
#'   \code{log_phi}, \code{logit_p}, \code{log_shape}, \code{log_nu},
#'   \code{logit_pi}, \code{power_p}. Unrecognized kinds raise an error.
#'   Defaults to \code{NULL} (empty override list; canonical priors
#'   apply when K > 1). In Block 8.0 K = 1 in every built-in family, so
#'   the slot is inert; it becomes active under the multi-parametric
#'   extension queued for Block 8.1+.
#'
#' @return An object of class \code{gdpar_prior} containing the
#'   seven legacy character snippets above plus the
#'   \code{priors_by_kind} override list.
#'
#' @details
#' The defaults follow the standard weakly informative recommendations
#' of the Stan team and are calibrated for problems in which (i) the
#' covariates entering the additive and multiplicative bases are
#' centered and scaled to unit variance and (ii) theta_ref is on the
#' linear-predictor scale of the family. The package standardizes the
#' covariates internally and reports posterior summaries on both the
#' standardized scale (used during sampling) and the user's original
#' scale (after back-transformation).
#'
#' All scale parameters are declared on the positive real line in
#' Stan; the prior strings supplied here are interpreted as
#' positive-truncated when the parameter has a lower bound of zero in
#' the Stan model.
#'
#' @section Methodological notes:
#' For finite-dimensional parametric AMM specifications the conditions
#' (PRIOR-KL) and (PRIOR-THICK) of Block 4 are satisfied automatically
#' when the prior on the parameter is absolutely continuous with
#' positive density at the true parameter. Each of the defaults above
#' satisfies this property unconditionally on the relevant parameter
#' space. For non-parametric extensions the matching of prior
#' smoothness becomes a substantive question (see Block 4, Sections
#' 6.1 and 7); such extensions are deferred to a future version of the
#' package.
#'
#' Internal sampling parametrization for the multiplicative component:
#' the AMM canonical form \code{theta_i = theta_ref + a(x_i) +
#' b(x_i) * theta_ref + W_term} is preserved at the user-facing level,
#' but the Stan model samples \code{c_b = theta_ref * b_coef} as the
#' free parameter and reports \code{b_coef = c_b / theta_ref} as a
#' derived quantity. Centering condition (C3) is enforced empirically
#' by column-wise centering of Z_b in the AMM design constructor, not
#' by any constraint on c_b.
#' This linear reparametrization yields a strictly log-concave
#' (Gaussian-conditional) posterior in \code{(theta_ref, a, c_b)},
#' eliminating the bimodality that the non-linear parametrization
#' \code{(theta_ref, a, b_coef)} admits as an artefact of the term
#' \code{theta_ref * (Z_b * b_coef)}. The prior \code{sigma_b} is the
#' hierarchical scale of \code{c_b}; users supplying custom priors on
#' \code{sigma_b} should interpret it accordingly.
#'
#' The user is free to supply any Stan-syntax prior string. The package
#' performs a syntactic check at generation time but does not attempt
#' to verify the prior's mathematical properties.
#'
#' @section Dependencies:
#' None at construction time. The prior strings are inserted into the
#' Stan model by the code generator and parsed by Stan at compile
#' time.
#'
#' @references
#' See \code{vignette("v04_asymptotics_path1_bayesian", package = "gdpar")},
#' Section 10.1, for the role of (PRIOR-KL) and (PRIOR-THICK) in
#' posterior consistency and contraction.
#'
#' Stan Development Team (2024). Stan User's Guide, version 2.35,
#' Section "Prior Choice Recommendations".
#'
#' @examples
#' pr <- gdpar_prior()
#' print(pr)
#'
#' pr2 <- gdpar_prior(theta_ref = "normal(0, 5)")
#' print(pr2)
#'
#' @seealso \code{\link{gdpar}}, \code{\link{amm_spec}}
#' @export
gdpar_prior <- function(theta_ref       = "normal(0, 2.5)",
                        sigma_theta_ref = "student_t(3, 0, 1)",
                        sigma_a         = "student_t(3, 0, 1)",
                        sigma_b         = "student_t(3, 0, 1)",
                        sigma_W         = "student_t(3, 0, 1)",
                        sigma_y         = "student_t(3, 0, 2.5)",
                        phi             = "gamma(2, 0.1)",
                        priors_by_kind  = NULL) {
  for (nm in c("theta_ref", "sigma_theta_ref", "sigma_a", "sigma_b",
               "sigma_W", "sigma_y", "phi")) {
    val <- get(nm)
    if (!is.character(val) || length(val) != 1L || nchar(val) == 0L) {
      gdpar_abort(
        sprintf(
          "Argument '%s' must be a non-empty character scalar (Stan prior syntax).",
          nm
        ),
        class = "gdpar_input_error",
        data = list(argument = nm, received = val)
      )
    }
  }
  if (is.null(priors_by_kind)) {
    priors_by_kind <- list()
  } else {
    if (!is.list(priors_by_kind) || is.null(names(priors_by_kind)) ||
        any(names(priors_by_kind) == "")) {
      gdpar_abort(
        "Argument 'priors_by_kind' must be a named list of Stan-syntax prior strings.",
        class = "gdpar_input_error"
      )
    }
    known <- .gdpar_known_canonical_kinds()
    unknown <- setdiff(names(priors_by_kind), known)
    if (length(unknown) > 0L) {
      gdpar_abort(
        sprintf(
          "Argument 'priors_by_kind' contains unrecognized kind(s): %s. Allowed: %s.",
          paste(sQuote(unknown), collapse = ", "),
          paste(sQuote(known), collapse = ", ")
        ),
        class = "gdpar_input_error",
        data = list(unknown = unknown, allowed = known)
      )
    }
    for (k in names(priors_by_kind)) {
      v <- priors_by_kind[[k]]
      if (!is.character(v) || length(v) != 1L || nchar(v) == 0L) {
        gdpar_abort(
          sprintf(
            "Entry 'priors_by_kind$%s' must be a non-empty character scalar (Stan prior syntax).",
            k
          ),
          class = "gdpar_input_error",
          data = list(kind = k, received = v)
        )
      }
    }
  }
  obj <- list(
    theta_ref       = theta_ref,
    sigma_theta_ref = sigma_theta_ref,
    sigma_a         = sigma_a,
    sigma_b         = sigma_b,
    sigma_W         = sigma_W,
    sigma_y         = sigma_y,
    phi             = phi,
    priors_by_kind  = priors_by_kind
  )
  class(obj) <- c("gdpar_prior", "list")
  obj
}

#' Internal: resolve the effective prior for a given canonical kind
#'
#' Returns the user-supplied override from \code{prior$priors_by_kind}
#' if present, otherwise the canonical default registered in
#' \code{.gdpar_canonical_prior_for_kind()}. Used by the codegen and by
#' tests of the canonical-kind override mechanism.
#'
#' @param prior A \code{gdpar_prior} object.
#' @param kind Character scalar.
#' @return Character scalar with the Stan-syntax prior expression.
#' @keywords internal
#' @noRd
.gdpar_prior_for_kind <- function(prior, kind) {
  override <- prior$priors_by_kind[[kind]]
  if (!is.null(override)) {
    return(override)
  }
  .gdpar_canonical_prior_for_kind(kind)
}

#' Print method for gdpar_prior objects
#'
#' @param x An object of class \code{gdpar_prior}.
#' @param ... Unused; present for S3 generic compatibility.
#' @return Invisibly returns \code{x}.
#' @export
print.gdpar_prior <- function(x, ...) {
  cat("<gdpar_prior>\n")
  for (nm in c("theta_ref", "sigma_theta_ref", "sigma_a", "sigma_b",
               "sigma_W", "sigma_y", "phi")) {
    cat("  ", format(nm, width = 15), " : ", x[[nm]], "\n", sep = "")
  }
  if (!is.null(x$priors_by_kind) && length(x$priors_by_kind) > 0L) {
    cat("  priors_by_kind  :\n")
    for (k in names(x$priors_by_kind)) {
      cat("    ", format(k, width = 13), " : ", x$priors_by_kind[[k]],
          "\n", sep = "")
    }
  }
  invisible(x)
}

#' Bernstein-von Mises calibration check (opt-in, costly)
#'
#' Verify numerically the conclusion of the Bernstein-von Mises
#' theorem (Theorem 4C of Block 4) for a fitted Path 1 model: that the
#' posterior is asymptotically Gaussian around the maximum likelihood
#' estimator with covariance equal to the inverse Fisher information
#' matrix divided by n. The function refits the model by maximum
#' likelihood (MLE) using a derived Stan model in which the prior
#' block is stripped (see \code{generate_stan_code(mle = TRUE)} and the
#' \code{// BEGIN PRIORS} / \code{// END PRIORS} markers in
#' \code{inst/stan/amm_main.stan}); it computes a Hessian-based
#' covariance estimate via a Laplace approximation around the MLE,
#' both on the unconstrained scale (the Stan optimizer is invoked with
#' \code{jacobian = FALSE}); and it compares the resulting interval
#' coverage with the Bayesian posterior intervals reported by the
#' fitted model.
#'
#' This function is opt-in and computationally expensive: it refits
#' the model in MLE mode and inverts a Hessian matrix. It is intended
#' as a methodological audit, not as part of the standard inference
#' flow. The conclusions of \code{\link{gdpar}} are not affected by
#' calling this function.
#'
#' @param fit An object of class \code{gdpar_fit} produced by
#'   \code{\link{gdpar}} with \code{path = "bayes"}.
#' @param parameters Optional character vector of parameter names to
#'   include in the comparison. Defaults to the user-facing parameters
#'   that the prior-stripped likelihood identifies: \code{theta_ref},
#'   \code{sigma_y} (when present), and \code{phi} (when present).
#' @param level Numeric scalar in (0, 1) with the nominal credible /
#'   confidence level. Defaults to 0.95.
#' @param verbose Logical scalar; when TRUE, prints an estimated cost
#'   message before starting. Defaults to TRUE.
#'
#' @return A list of class \code{gdpar_bvm_report} with components
#'   \code{table} (data frame comparing posterior intervals with
#'   Hessian-based intervals per parameter), \code{discrepancy}
#'   (numeric vector of relative interval-width differences),
#'   \code{level} and \code{warnings}. A \code{print} method provides
#'   a human-readable summary.
#'
#' @details
#' Theorem 4C of Block 4 establishes that, for finite-dimensional
#' parametric AMM specifications under the (LAN) condition with
#' non-singular Fisher information at the true parameter, the
#' posterior distribution converges in total variation to the
#' asymptotic-Gaussian distribution of the maximum likelihood
#' estimator. Empirically, this entails that posterior credible
#' intervals at any nominal level should agree with the
#' Hessian-based asymptotic confidence intervals as the sample size
#' grows.
#'
#' This function performs the empirical comparison at the observed
#' sample size. Substantial discrepancy between the two interval
#' families at large n signals either (i) the limit has not yet been
#' approached, requiring more data; (ii) the (LAN) condition fails
#' (e.g., singular Fisher information at the true parameter); or
#' (iii) the model is misspecified (Block 7).
#'
#' Applies only to fits with finite-dimensional parametric AMM
#' specifications. Non-parametric components (planned for a future
#' version) are outside the scope of Theorem 4C; the function will
#' abort if invoked on a non-parametric fit. The hierarchical regime
#' activated by the \code{group} argument of \code{gdpar()} (Block 6.5)
#' is likewise out of scope, because the classical asymptotic theory
#' that underwrites Theorem 4C assumes a fixed-dimension parameter
#' vector while the hierarchical case introduces an increasing number
#' of random anchors; the function aborts with
#' \code{gdpar_unsupported_feature_error} when invoked on a grouped
#' fit.
#'
#' @section Methodological notes:
#' This function is part of the methodological audit toolkit, not of
#' the standard inference flow. It does not modify the \code{fit}
#' object in any way; its output is informational. Users running
#' large simulations are advised to call this function selectively
#' rather than after every fit.
#'
#' The MLE estimate uses Stan's \code{optimize} method with the LBFGS
#' algorithm on the prior-stripped variant of the model produced by
#' \code{generate_stan_code(mle = TRUE)}; the optimizer is invoked
#' with \code{jacobian = FALSE} so the maximum is taken on the
#' constrained (natural) scale of the parameters rather than on the
#' unconstrained scale that would carry the Jacobian of the
#' transformation. The Hessian-based interval estimate comes from
#' \code{cmdstanr::laplace} applied around the MLE on the same
#' prior-stripped model and with the same \code{jacobian = FALSE}
#' convention.
#'
#' The comparison is restricted to the user-facing parameters that
#' the likelihood identifies in MLE mode without the prior anchoring:
#' the global reference parameter \code{theta_ref}, the response-level
#' scale \code{sigma_y} (Gaussian families), and the dispersion
#' \code{phi} (Negative Binomial). The hierarchical scales
#' \code{sigma_a}, \code{sigma_b}, \code{sigma_W} are not identified
#' by the likelihood alone (they enter only through the random-effect
#' priors that have been stripped) and are excluded from the table.
#'
#' @section Dependencies:
#' Uses \pkg{cmdstanr} for the optimization run and \pkg{posterior}
#' to extract the Bayesian intervals.
#'
#' @references
#' See \code{vignette("v04_asymptotics_path1_bayesian", package = "gdpar")},
#' Section 7 (Theorem 4C).
#'
#' @examples
#' \donttest{
#' if (requireNamespace("cmdstanr", quietly = TRUE)) {
#'   df <- data.frame(x1 = rnorm(200), y = rnorm(200))
#'   fit <- gdpar(y ~ x1, amm = amm_spec(a = ~ x1), data = df,
#'                iter_warmup = 200, iter_sampling = 200, chains = 2)
#'   gdpar_bvm_check(fit)
#' }
#' }
#'
#' @seealso \code{\link{gdpar}}, \code{\link{gdpar_contraction_diagnostic}}
#' @export
gdpar_bvm_check <- function(fit, parameters = NULL, level = 0.95,
                            verbose = TRUE) {
  assert_inherits(fit, "gdpar_fit", "fit")
  assert_numeric_scalar(level, "level", lower = 0, upper = 1)
  if (!is.logical(verbose) || length(verbose) != 1L) {
    gdpar_abort("Argument 'verbose' must be a logical scalar.",
                class = "gdpar_input_error")
  }
  use_groups_in_fit <- isTRUE(
    !is.null(fit$stan_data$use_groups) &&
      as.integer(fit$stan_data$use_groups) == 1L
  )
  if (use_groups_in_fit) {
    gdpar_abort(
      paste(
        "gdpar_bvm_check does not apply to hierarchical fits (Block 6.5).",
        "The Bernstein-von Mises theorem of Block 4 (Theorem 4C) is",
        "stated for finite-dimensional parametric specifications;",
        "per-group anchors introduce a random-effect dimension that",
        "is outside its scope. Refit without the 'group' argument to",
        "audit a non-hierarchical fit."
      ),
      class = "gdpar_unsupported_feature_error"
    )
  }
  require_suggested("cmdstanr",
                    "fit Path 1 hierarchical Bayesian models")
  require_suggested("posterior", "extract posterior draws")
  if (verbose) {
    gdpar_inform(
      paste(
        "gdpar_bvm_check is an opt-in methodological audit.",
        "It refits the model by maximum likelihood (prior block",
        "stripped) and inverts a Hessian via Laplace approximation.",
        "The output does not modify the fit; it complements it with",
        "a calibration comparison."
      ),
      class = "gdpar_optin_message"
    )
  }
  draws <- fit$fit$draws()
  vars <- posterior::variables(draws)
  ignore_pat <- paste0(
    "^(eta|log_lik|y_pred|theta_i|a_coef|b_coef|a_raw|b_raw|W_raw|",
    "c_b|c_b_raw|mu_theta_ref|sigma_theta_ref|sigma_a|sigma_b|sigma_W)"
  )
  candidate_vars <- vars[!grepl(ignore_pat, vars)]
  if (is.null(parameters)) {
    parameters <- candidate_vars
  } else {
    missing_p <- setdiff(parameters, candidate_vars)
    if (length(missing_p) > 0L) {
      gdpar_abort(
        sprintf(
          "Parameter(s) not available in the fit: %s.",
          paste(sQuote(missing_p), collapse = ", ")
        ),
        class = "gdpar_input_error"
      )
    }
  }
  alpha <- 1 - level
  ql_val <- alpha / 2
  qu_val <- 1 - alpha / 2
  fn_q_lower <- function(x) stats::quantile(x, ql_val, names = FALSE)
  fn_q_upper <- function(x) stats::quantile(x, qu_val, names = FALSE)
  bayes_summary <- posterior::summarise_draws(
    posterior::subset_draws(draws, variable = parameters),
    mean = function(x) mean(x),
    sd = function(x) stats::sd(x),
    q_lower = fn_q_lower,
    q_upper = fn_q_upper
  )
  stan_data <- fit$stan_data
  src <- generate_stan_code(fit$prior, mle = TRUE)
  stan_path <- write_stan_to_tempfile(src)
  cs_model <- cmdstanr::cmdstan_model(stan_path)
  opt <- cs_model$optimize(
    data = stan_data,
    algorithm = "lbfgs",
    refresh = 0,
    history_size = 5,
    init_alpha = 0.001,
    iter = 2000,
    jacobian = FALSE
  )
  hessian_table <- tryCatch({
    laplace <- cs_model$laplace(
      mode = opt,
      data = stan_data,
      draws = 1000,
      refresh = 0,
      jacobian = FALSE
    )
    laplace_summary <- posterior::summarise_draws(
      posterior::subset_draws(laplace$draws(), variable = parameters),
      mean = function(x) mean(x),
      sd = function(x) stats::sd(x),
      q_lower = fn_q_lower,
      q_upper = fn_q_upper
    )
    laplace_summary
  }, error = function(e) {
    gdpar_warn(
      sprintf(
        "Laplace approximation failed: %s. Falling back to MLE-only summary.",
        conditionMessage(e)
      ),
      class = "gdpar_diagnostic_warning"
    )
    NULL
  })
  if (is.null(hessian_table)) {
    table <- data.frame(
      variable = parameters,
      bayes_mean = bayes_summary$mean,
      bayes_lower = bayes_summary$q_lower,
      bayes_upper = bayes_summary$q_upper,
      bayes_width = bayes_summary$q_upper - bayes_summary$q_lower,
      asymp_mean = NA_real_,
      asymp_lower = NA_real_,
      asymp_upper = NA_real_,
      asymp_width = NA_real_,
      width_ratio = NA_real_
    )
    discrepancy <- rep(NA_real_, length(parameters))
    warnings_msg <- "Laplace approximation failed; asymptotic comparison unavailable."
  } else {
    table <- data.frame(
      variable = parameters,
      bayes_mean = bayes_summary$mean,
      bayes_lower = bayes_summary$q_lower,
      bayes_upper = bayes_summary$q_upper,
      bayes_width = bayes_summary$q_upper - bayes_summary$q_lower,
      asymp_mean = hessian_table$mean,
      asymp_lower = hessian_table$q_lower,
      asymp_upper = hessian_table$q_upper,
      asymp_width = hessian_table$q_upper - hessian_table$q_lower
    )
    table$width_ratio <- table$bayes_width / table$asymp_width
    discrepancy <- abs(log(table$width_ratio))
    warnings_msg <- character(0)
    bad <- table$width_ratio < 0.5 | table$width_ratio > 2
    if (any(bad, na.rm = TRUE)) {
      warnings_msg <- c(warnings_msg,
                        sprintf(
                          "Width ratio outside [0.5, 2] for: %s.",
                          paste(table$variable[bad], collapse = ", ")
                        ))
    }
  }
  obj <- list(
    table = table,
    discrepancy = discrepancy,
    level = level,
    warnings = warnings_msg
  )
  class(obj) <- c("gdpar_bvm_report", "list")
  obj
}

#' Print method for gdpar_bvm_report objects
#'
#' @param x An object of class \code{gdpar_bvm_report}.
#' @param ... Unused; present for S3 generic compatibility.
#' @return Invisibly returns \code{x}.
#' @export
print.gdpar_bvm_report <- function(x, ...) {
  cat("<gdpar_bvm_report> level = ", x$level, "\n", sep = "")
  print(x$table, row.names = FALSE)
  if (length(x$warnings) > 0L) {
    cat("\nWarnings:\n")
    for (w in x$warnings) cat("  - ", w, "\n", sep = "")
  }
  invisible(x)
}

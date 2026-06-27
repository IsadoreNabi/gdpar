#' Empirical posterior contraction rate diagnostic (opt-in, costly)
#'
#' Verify numerically the predicted posterior contraction rate
#' (Theorem 4B of Block 4) for a fitted Path 1 model. Refits the model
#' at multiple subsample sizes, records the median posterior credible
#' interval width across the user-facing parameters, and fits the
#' regression \eqn{\log(\text{width}) = \alpha + \beta \log(n)}
#' across the subsample sizes. The slope estimate is consistent with
#' a parametric \eqn{n^{-1/2}} contraction rate when its value is in
#' the interval (-0.6, -0.4).
#'
#' This function is opt-in and computationally expensive: it refits
#' the model \code{length(sizes) * replicates} times. A cost message
#' is printed at the start. The conclusions of \code{\link{gdpar}}
#' are not affected by calling this function.
#'
#' @param fit An object of class \code{gdpar_fit} produced by
#'   \code{\link{gdpar}} with \code{path = "bayes"}.
#' @param data The data frame originally passed to
#'   \code{\link{gdpar}} (or another data frame compatible with the
#'   AMM specification of \code{fit}).
#' @param sizes Integer vector with the subsample sizes at which to
#'   refit. Defaults to a length-five geometric sequence between
#'   \code{ceiling(n / 8)} and \code{n}.
#' @param replicates Integer scalar with the number of independent
#'   subsamples per size. Defaults to 1; a higher value reduces
#'   Monte Carlo variance of the curve at additional cost.
#' @param parameters Optional character vector of parameter names to
#'   include in the credible-width calculation. Defaults to the
#'   user-facing parameters.
#' @param level Numeric scalar in (0, 1) with the nominal credible
#'   level used for the width calculation. Defaults to 0.95.
#' @param iter_warmup Integer scalar; warmup iterations for each
#'   refit. Defaults to 500.
#' @param iter_sampling Integer scalar; sampling iterations for each
#'   refit. Defaults to 500.
#' @param chains Integer scalar; chains per refit. Defaults to 2.
#' @param verbose Logical scalar; when TRUE, prints an estimated cost
#'   message before starting. Defaults to TRUE.
#' @param ... Additional arguments forwarded to \code{\link{gdpar}}.
#'
#' @return A list of class \code{gdpar_contraction_report} with
#'   components \code{table} (data frame with columns \code{n},
#'   \code{replicate}, \code{median_width}), \code{slope_estimate},
#'   \code{slope_se}, \code{slope_ci_lower}, \code{slope_ci_upper},
#'   \code{verdict} (a character indicating whether the empirical
#'   slope is consistent with the parametric \eqn{n^{-1/2}} rate),
#'   and \code{warnings} (character vector recording per-refit
#'   fallback notifications; empty when every refit succeeded).
#'   A \code{print} method provides a human-readable summary.
#'
#' @details
#' Theorem 4B of Block 4 establishes that, under the conditions of
#' Theorem 4A plus the prior thickness condition (PRIOR-THICK) and
#' the sieve condition (SIEVE), the posterior contracts at rate
#' \eqn{\varepsilon_n} with \eqn{n \varepsilon_n^2 \to \infty}. For
#' finite-dimensional parametric AMM specifications, the rate is
#' \eqn{n^{-1/2}}. This diagnostic checks the empirical slope of the
#' log-width against log-n.
#'
#' Deviations from the predicted slope can indicate (i) prior
#' misspecification (the prior fails (PRIOR-THICK) at the true
#' parameter), (ii) failure of the homogeneity (HOM) or regularity
#' (REG) conditions of Block 2, or (iii) non-parametric components
#' whose smoothness assumption does not match the truth. The
#' diagnostic flags the discrepancy without diagnosing the cause.
#'
#' @section Methodological notes:
#' This function is part of the methodological audit toolkit, not of
#' the standard inference flow. It is computationally expensive: each
#' subsample size requires a full refit. The default settings keep
#' the refits short (500 + 500 iterations, 2 chains) to make the
#' diagnostic affordable on moderate datasets; users with more
#' computational budget should pass higher values via the relevant
#' arguments.
#'
#' Subsamples are drawn without replacement and stratified by row
#' order; users with structured data (time series, clustered
#' observations) should pass an explicit \code{sizes} vector that
#' respects the structure if random subsampling is inappropriate.
#'
#' @section Dependencies:
#' Uses \pkg{cmdstanr} for the refits and \pkg{posterior} to extract
#' credible-interval widths.
#'
#' @references
#' See \code{vignette("v04_asymptotics_path1_bayesian", package = "gdpar")},
#' Section 6.2 (numerical verification of contraction).
#'
#' @examples
#' \donttest{
#' if (requireNamespace("cmdstanr", quietly = TRUE)) {
#'   df <- data.frame(x1 = rnorm(400), y = rnorm(400))
#'   fit <- gdpar(y ~ x1, amm = amm_spec(a = ~ x1), data = df,
#'                iter_warmup = 200, iter_sampling = 200, chains = 2)
#'   gdpar_contraction_diagnostic(fit, data = df, replicates = 1)
#' }
#' }
#'
#' @seealso \code{\link{gdpar}}, \code{\link{gdpar_bvm_check}}
#' @export
gdpar_contraction_diagnostic <- function(fit, data, sizes = NULL,
                                          replicates = 1L,
                                          parameters = NULL,
                                          level = 0.95,
                                          iter_warmup = 500L,
                                          iter_sampling = 500L,
                                          chains = 2L,
                                          verbose = TRUE,
                                          ...) {
  assert_inherits(fit, "gdpar_fit", "fit")
  assert_data_frame(data, "data")
  assert_count(replicates, "replicates")
  assert_numeric_scalar(level, "level", lower = 0, upper = 1)
  assert_count(iter_warmup, "iter_warmup")
  assert_count(iter_sampling, "iter_sampling")
  assert_count(chains, "chains")
  if (!is.logical(verbose) || length(verbose) != 1L) {
    gdpar_abort("Argument 'verbose' must be a logical scalar.",
                class = "gdpar_input_error")
  }
  require_suggested("cmdstanr",
                    "fit Path 1 hierarchical Bayesian models")
  require_suggested("posterior", "extract posterior draws")
  n <- nrow(data)
  if (is.null(sizes)) {
    sizes <- unique(round(exp(seq(log(max(20, ceiling(n / 8))),
                                  log(n), length.out = 5))))
  }
  if (!is.numeric(sizes) || any(sizes < 5) || any(sizes > n)) {
    gdpar_abort(
      sprintf(
        "Argument 'sizes' must be a numeric vector with entries in [5, %d].",
        n
      ),
      class = "gdpar_input_error"
    )
  }
  if (verbose) {
    n_refits <- length(sizes) * replicates
    gdpar_inform(
      sprintf(
        paste(
          "gdpar_contraction_diagnostic is an opt-in methodological audit.",
          "It will perform %d refit(s) (%d sizes x %d replicates).",
          "The output does not modify the fit; it complements it with",
          "an empirical contraction-rate estimate."
        ),
        n_refits, length(sizes), replicates
      ),
      class = "gdpar_optin_message"
    )
  }
  call <- fit$call
  call_data_arg_name <- "data"
  alpha <- 1 - level
  ql_val <- alpha / 2
  qu_val <- 1 - alpha / 2
  fn_q_lower <- function(x) stats::quantile(x, ql_val, names = FALSE)
  fn_q_upper <- function(x) stats::quantile(x, qu_val, names = FALSE)
  rows <- list()
  warnings_msg <- character(0L)
  for (s in seq_along(sizes)) {
    sz <- sizes[s]
    for (r in seq_len(replicates)) {
      idx <- sample.int(n, size = sz)
      sub <- data[idx, , drop = FALSE]
      sub_call <- call
      sub_call$data <- quote(sub)
      sub_call$iter_warmup <- iter_warmup
      sub_call$iter_sampling <- iter_sampling
      sub_call$chains <- chains
      sub_call$verbose <- FALSE
      sub_call$refresh <- 0L
      sub_call$skip_id_check <- TRUE
      env <- new.env(parent = parent.frame())
      env$sub <- sub
      refit_failure_msg <- NULL
      sub_fit <- tryCatch(
        eval(sub_call, envir = env),
        error = function(e) {
          refit_failure_msg <<- sprintf(
            "Refit at n = %d (replicate %d) failed: %s.",
            sz, r, conditionMessage(e)
          )
          gdpar_warn(
            sprintf(
              paste(
                "%s Falling back to median_width = NA for this cell; the",
                "contraction slope will be estimated from the remaining",
                "successful refits (see $warnings of the returned report)."
              ),
              refit_failure_msg
            ),
            class = "gdpar_diagnostic_warning"
          )
          NULL
        }
      )
      if (is.null(sub_fit)) {
        if (!is.null(refit_failure_msg)) {
          warnings_msg <- c(warnings_msg, refit_failure_msg)
        }
        rows[[length(rows) + 1L]] <- data.frame(
          n = sz, replicate = r, median_width = NA_real_
        )
        next
      }
      draws <- sub_fit$fit$draws()
      vars <- posterior::variables(draws)
      ignore_pat <- "^(eta|log_lik|y_pred|theta_i|a_coef|b_coef|a_raw|b_raw|W_raw)"
      candidate_vars <- vars[!grepl(ignore_pat, vars)]
      use_vars <- if (is.null(parameters)) candidate_vars else parameters
      sub_summary <- posterior::summarise_draws(
        posterior::subset_draws(draws, variable = use_vars),
        q_lower = fn_q_lower,
        q_upper = fn_q_upper
      )
      widths <- sub_summary$q_upper - sub_summary$q_lower
      rows[[length(rows) + 1L]] <- data.frame(
        n = sz, replicate = r,
        median_width = stats::median(widths)
      )
    }
  }
  table <- do.call(rbind, rows)
  table_ok <- table[!is.na(table$median_width), , drop = FALSE]
  if (nrow(table_ok) < 3L) {
    gdpar_abort(
      "Not enough successful refits to estimate the contraction slope.",
      class = "gdpar_diagnostic_error"
    )
  }
  log_n <- log(table_ok$n)
  log_w <- log(table_ok$median_width)
  reg <- stats::lm(log_w ~ log_n)
  cf <- stats::coef(reg)
  se <- summary(reg)$coefficients[, "Std. Error"]
  slope <- cf["log_n"]
  slope_se <- se["log_n"]
  slope_lower <- slope - 1.96 * slope_se
  slope_upper <- slope + 1.96 * slope_se
  verdict <- if (slope_upper >= -0.6 && slope_lower <= -0.4) {
    "Consistent with parametric n^{-1/2} rate."
  } else if (slope_upper < -0.6) {
    "Faster than n^{-1/2}; check for spurious artefacts."
  } else {
    "Slower than n^{-1/2}; check for prior misspecification or model misspecification."
  }
  obj <- list(
    table = table,
    slope_estimate = unname(slope),
    slope_se = unname(slope_se),
    slope_ci_lower = unname(slope_lower),
    slope_ci_upper = unname(slope_upper),
    verdict = verdict,
    level = level,
    warnings = warnings_msg
  )
  class(obj) <- c("gdpar_contraction_report", "list")
  obj
}

#' Print method for gdpar_contraction_report objects
#'
#' @param x An object of class \code{gdpar_contraction_report}.
#' @param ... Unused; present for S3 generic compatibility.
#' @return Invisibly returns \code{x}.
#' @export
print.gdpar_contraction_report <- function(x, ...) {
  cat("<gdpar_contraction_report> level = ", x$level, "\n", sep = "")
  print(x$table, row.names = FALSE)
  cat("\nSlope estimate (log_width ~ log_n): ",
      format(x$slope_estimate, digits = 3),
      " (SE = ", format(x$slope_se, digits = 3), ")\n", sep = "")
  cat("95% CI: [",
      format(x$slope_ci_lower, digits = 3), ", ",
      format(x$slope_ci_upper, digits = 3), "]\n", sep = "")
  cat("Verdict: ", x$verdict, "\n", sep = "")
  if (length(x$warnings) > 0L) {
    cat("\nWarnings:\n")
    for (w in x$warnings) cat("  - ", w, "\n", sep = "")
  }
  invisible(x)
}

#' Print method for gdpar_eb_fit
#'
#' Concise console summary of an Empirical-Bayes fit: family, link, AMM
#' level, EB point estimate(s) of theta_ref, numerical diagnostics of
#' the Step (i) Laplace approximation, and a one-line summary of the
#' conditional HMC fit.
#'
#' @param x A \code{gdpar_eb_fit} object.
#' @param digits Integer scalar passed to \code{format()}; defaults to
#'   3.
#' @param ... Unused.
#' @return The object \code{x} invisibly.
#' @export
print.gdpar_eb_fit <- function(x, digits = 3L, ...) {
  cat("<gdpar_eb_fit>\n")
  fam <- x$family
  if (!is.null(fam)) {
    cat(sprintf("  family       : %s (link = %s)\n",
                fam$name %||% "?", fam$link %||% "?"))
  }
  is_path_C <- identical(x$path, "eb_KxP")
  if (is_path_C) {
    K <- x$K %||% (x$stan_data$K_slots %||% 1L)
    p <- x$p %||% (x$stan_data$p_dim %||% 1L)
    J <- dim(x$theta_ref_kp_hat)[1L] %||% 1L
    cat(sprintf("  path         : Path C (K = %d, p = %d, K x p regime)\n",
                K, p))
    cat(sprintf("  J_groups     : %d\n", J))
    cat(sprintf("  slot_names   : %s\n",
                paste(x$slot_names, collapse = ", ")))
    cat("  theta_ref_kp_hat (group 1, per slot per coord):\n")
    for (k in seq_len(K)) {
      cat(sprintf("    %s : %s\n",
                  x$slot_names[k],
                  paste(format(x$theta_ref_kp_hat[1L, k, ],
                                digits = digits),
                        collapse = " ")))
    }
    cat("  theta_ref_kp_se (group 1, per slot per coord):\n")
    for (k in seq_len(K)) {
      cat(sprintf("    %s : %s\n",
                  x$slot_names[k],
                  paste(format(x$theta_ref_kp_se[1L, k, ],
                                digits = digits),
                        collapse = " ")))
    }
    diag_num <- x$diagnostics_numerical %||% list()
    if (!is.null(diag_num$kappa_per_slot)) {
      cat(sprintf("  kappa(H) per slot: %s\n",
                  paste(format(diag_num$kappa_per_slot, digits = digits),
                        collapse = ", ")))
    }
    cat(sprintf("  multi_start_dispersion: %s\n",
                format(diag_num$multi_start_dispersion %||% NA_real_,
                       digits = digits)))
    cat(sprintf("  correction   : %s\n",
                if (isTRUE(x$correction_applied))
                  "Proposition 7B* (tensor full K x p x p, D37 = (i))"
                else "off"))
    if (!is.null(x$correction_tensor_dispositions)) {
      cat(sprintf("  slot dispositions: %s\n",
                  paste(x$correction_tensor_dispositions, collapse = ", ")))
    }
    diag_hmc <- x$diagnostics
    if (!is.null(diag_hmc)) {
      cat(sprintf(
        "  conditional HMC: R-hat <= %s, min ESS = %s, divergent = %d\n",
        format(diag_hmc$rhat_max, digits = digits),
        format(diag_hmc$ess_bulk_min, digits = digits),
        diag_hmc$divergent_count %||% 0L
      ))
    }
    return(invisible(x))
  }
  J <- length(x$theta_ref_hat)
  cat(sprintf("  J_groups     : %d\n", J))
  cat(sprintf("  K_slots      : %d\n", x$stan_data$K_slots %||% 1L))
  cat(sprintf("  p_dim        : %d\n", x$stan_data$p_dim %||% 1L))
  cat(sprintf("  anchor       : %s\n",
              format(x$anchor, digits = digits)))
  cat("  theta_ref_hat:")
  cat("", format(x$theta_ref_hat, digits = digits), "\n")
  cat("  theta_ref_se :")
  cat("", format(x$theta_ref_se, digits = digits), "\n")
  diag_num <- x$diagnostics_numerical %||% list()
  cat(sprintf("  kappa(H)     : %s\n",
              format(diag_num$kappa, digits = digits)))
  cat(sprintf("  multi_start_dispersion: %s\n",
              format(diag_num$multi_start_dispersion, digits = digits)))
  cat(sprintf("  correction   : %s\n",
              if (isTRUE(x$correction_applied)) "Proposition 7B (scalar)"
              else "off"))
  diag_hmc <- x$diagnostics
  if (!is.null(diag_hmc)) {
    cat(sprintf("  conditional HMC: R-hat <= %s, min ESS = %s, divergent = %d\n",
                format(diag_hmc$rhat_max, digits = digits),
                format(diag_hmc$ess_bulk_min, digits = digits),
                diag_hmc$divergent_count %||% 0L))
  }
  invisible(x)
}

#' Summary method for gdpar_eb_fit
#'
#' Returns a structured summary suitable for programmatic access and
#' for the canonical \code{print.summary.gdpar_eb_fit} method. Inflates
#' the conditional credible intervals by the Proposition 7B scalar
#' correction when \code{eb_correction = TRUE} was requested at fit
#' time.
#'
#' @param object A \code{gdpar_eb_fit} object.
#' @param level Numeric scalar in (0, 1); credible-interval level.
#'   Defaults to 0.95.
#' @param ... Unused.
#' @return An object of class \code{summary.gdpar_eb_fit}.
#' @export
summary.gdpar_eb_fit <- function(object, level = 0.95, ...) {
  if (!is.numeric(level) || length(level) != 1L ||
      level <= 0 || level >= 1) {
    gdpar_abort(
      "Argument 'level' must be a single numeric value in (0, 1).",
      class = "gdpar_input_error"
    )
  }
  alpha <- 1 - level
  z <- stats::qnorm(1 - alpha / 2)
  is_path_C <- identical(object$path, "eb_KxP")
  if (is_path_C) {
    K <- object$K %||% (object$stan_data$K_slots %||% 1L)
    p <- object$p %||% (object$stan_data$p_dim %||% 1L)
    J <- dim(object$theta_ref_kp_hat)[1L] %||% 1L
    rows <- vector("list", 0L)
    tensor <- object$correction_tensor_constant
    correction_on <- isTRUE(object$correction_applied)
    for (g in seq_len(J)) {
      for (k in seq_len(K)) {
        for (c in seq_len(p)) {
          est_gkc <- object$theta_ref_kp_hat[g, k, c]
          se_gkc <- object$theta_ref_kp_se[g, k, c]
          # Per-slot inflation: diagonal entry of the slot block of the
          # correction tensor (kappa(alpha) * Sigma^marg_{theta_ref,k}).
          inflate_kc <- if (correction_on && !is.null(tensor) &&
                            all(is.finite(tensor[k, c, c]))) {
                          sqrt(1 + tensor[k, c, c] / max(1, J))
                        } else {
                          1
                        }
          lo <- est_gkc - z * se_gkc * inflate_kc
          hi <- est_gkc + z * se_gkc * inflate_kc
          rows[[length(rows) + 1L]] <- data.frame(
            group    = g,
            slot     = object$slot_names[k],
            coord    = c,
            estimate = est_gkc,
            se       = se_gkc,
            lower    = lo,
            upper    = hi,
            inflation = inflate_kc,
            stringsAsFactors = FALSE
          )
        }
      }
    }
    theta_table <- do.call(rbind, rows)
    conditional_summary <- tryCatch(
      {
        fit_cond <- object$conditional_fit
        if (!is.null(fit_cond) &&
            requireNamespace("posterior", quietly = TRUE)) {
          draws <- fit_cond$draws()
          vars <- posterior::variables(draws)
          keep <- vars[!grepl(
            "^(eta|log_lik|y_pred|theta_i|a_raw|c_b_raw|c_b_kp_raw|W_raw|a_coef_kp|c_b_kp|b_coef_kp|eta_kp)",
            vars
          )]
          if (length(keep) > 0L) {
            posterior::summarise_draws(
              posterior::subset_draws(draws, variable = keep),
              posterior::default_summary_measures(),
              posterior::default_convergence_measures()
            )
          } else NULL
        } else NULL
      },
      error = function(e) NULL
    )
    out <- list(
      theta_table              = theta_table,
      conditional_summary      = conditional_summary,
      correction_applied       = correction_on,
      correction_tensor        = tensor,
      slot_dispositions        = object$correction_tensor_dispositions,
      level                    = level,
      family                   = object$family$name %||% NA_character_,
      link                     = object$family$link %||% NA_character_,
      J_groups                 = J,
      K_slots                  = K,
      p_dim                    = p,
      slot_names               = object$slot_names,
      diagnostics_numerical    = object$diagnostics_numerical,
      diagnostics_hmc          = object$diagnostics,
      path                     = "eb_KxP",
      call                     = object$call
    )
    class(out) <- c("summary.gdpar_eb_fit", "list")
    return(out)
  }
  inflate <- if (isTRUE(object$correction_applied)) {
    sqrt(1 + (object$eb_correction_constant %||% 0) /
                max(1, length(object$theta_ref_hat)))
  } else {
    1
  }
  theta_table <- data.frame(
    estimate = object$theta_ref_hat,
    se       = object$theta_ref_se,
    lower    = object$theta_ref_hat - z * object$theta_ref_se * inflate,
    upper    = object$theta_ref_hat + z * object$theta_ref_se * inflate,
    stringsAsFactors = FALSE
  )
  if (!is.null(object$group_info)) {
    rownames(theta_table) <- object$group_info$levels
  }
  conditional_summary <- tryCatch(
    {
      fit_cond <- object$conditional_fit
      if (!is.null(fit_cond) && requireNamespace("posterior", quietly = TRUE)) {
        draws <- fit_cond$draws()
        vars <- posterior::variables(draws)
        keep <- vars[!grepl("^(eta|log_lik|y_pred|theta_i|a_raw|c_b_raw|W_raw)",
                            vars)]
        if (length(keep) > 0L) {
          posterior::summarise_draws(
            posterior::subset_draws(draws, variable = keep),
            posterior::default_summary_measures(),
            posterior::default_convergence_measures()
          )
        } else NULL
      } else NULL
    },
    error = function(e) NULL
  )
  out <- list(
    theta_table         = theta_table,
    conditional_summary = conditional_summary,
    correction_applied  = isTRUE(object$correction_applied),
    correction_constant = object$eb_correction_constant,
    inflation_factor    = inflate,
    level               = level,
    family              = object$family$name %||% NA_character_,
    link                = object$family$link %||% NA_character_,
    J_groups            = length(object$theta_ref_hat),
    K_slots             = object$stan_data$K_slots %||% 1L,
    p_dim               = object$stan_data$p_dim %||% 1L,
    diagnostics_numerical = object$diagnostics_numerical,
    diagnostics_hmc       = object$diagnostics,
    call                  = object$call
  )
  class(out) <- c("summary.gdpar_eb_fit", "list")
  out
}

#' Print method for summary.gdpar_eb_fit
#'
#' @param x A \code{summary.gdpar_eb_fit} object.
#' @param digits Integer scalar passed to \code{format()}; defaults to
#'   3.
#' @param ... Unused.
#' @return The object \code{x} invisibly.
#' @export
print.summary.gdpar_eb_fit <- function(x, digits = 3L, ...) {
  cat("<gdpar_eb_fit summary>\n")
  cat(sprintf("  family       : %s (link = %s)\n",
              x$family, x$link))
  cat(sprintf("  J_groups     : %d ; K_slots = %d ; p_dim = %d\n",
              x$J_groups, x$K_slots, x$p_dim))
  cat(sprintf("  level        : %s\n", format(x$level, digits = digits)))
  is_path_C <- identical(x$path, "eb_KxP")
  if (is_path_C) {
    cat(sprintf("  path         : Path C (K x p regime)\n"))
    cat(sprintf("  correction   : %s\n",
                if (isTRUE(x$correction_applied)) sprintf(
                  "Proposition 7B* tensor full K x p x p (D37 = (i))"
                ) else "off"))
    if (!is.null(x$slot_dispositions)) {
      cat(sprintf("  slot dispositions: %s\n",
                  paste(x$slot_dispositions, collapse = ", ")))
    }
    cat("\n  theta_ref_kp estimates (per group, per slot, per coord):\n")
    print(format(x$theta_table, digits = digits))
    if (!is.null(x$conditional_summary)) {
      cat("\n  Conditional posterior summary (first 8 parameters):\n")
      print(utils::head(x$conditional_summary, 8L))
    }
    return(invisible(x))
  }
  cat(sprintf("  correction   : %s\n",
              if (x$correction_applied) sprintf(
                "Proposition 7B (scalar), C = %s, inflation = %s",
                format(x$correction_constant, digits = digits),
                format(x$inflation_factor, digits = digits))
              else "off"))
  cat("\n  theta_ref estimates:\n")
  print(format(x$theta_table, digits = digits))
  if (!is.null(x$conditional_summary)) {
    cat("\n  Conditional posterior summary (first 8 parameters):\n")
    print(utils::head(x$conditional_summary, 8L))
  }
  invisible(x)
}

#' Coefficient extraction for gdpar_eb_fit
#'
#' Returns an object analogous to \code{coef.gdpar_fit}: a
#' \code{gdpar_coef} list with components \code{theta_ref}, \code{a},
#' \code{b}, \code{W}. The EB version reports
#' \code{theta_ref$method == "EB"} and includes a flag
#' \code{theta_ref$eb_correction_applied} so downstream consumers can
#' tell that the credible intervals are EB-corrected (or that the
#' correction was disabled).
#'
#' @param object A \code{gdpar_eb_fit} object.
#' @param ... Forwarded to the underlying \code{coef.gdpar_fit}-like
#'   accessor on the conditional fit (none recognized in v0).
#' @return A list with class \code{c("gdpar_coef_eb", "gdpar_coef", "list")}.
#' @export
coef.gdpar_eb_fit <- function(object, ...) {
  fit_cond <- object$conditional_fit
  out <- list(
    theta_ref = list(
      method                    = "EB",
      estimate                  = object$theta_ref_hat,
      se                        = object$theta_ref_se,
      cov                       = object$diagnostics_numerical$theta_ref_cov,
      eb_correction_applied     = isTRUE(object$correction_applied),
      eb_correction_constant    = object$eb_correction_constant
    )
  )
  if (!is.null(fit_cond) &&
      requireNamespace("posterior", quietly = TRUE)) {
    draws <- fit_cond$draws()
    vars <- posterior::variables(draws)
    pick <- function(pat) {
      v <- grep(pat, vars, value = TRUE)
      if (length(v) == 0L) return(NULL)
      mat <- posterior::as_draws_matrix(
        posterior::subset_draws(draws, variable = v)
      )
      list(estimate = colMeans(mat),
            se      = apply(mat, 2L, stats::sd))
    }
    out$a <- pick("^a_coef")
    out$b <- pick("^b_coef")
    out$W <- pick("^W_raw")
  }
  class(out) <- c("gdpar_coef_eb", "gdpar_coef", "list")
  out
}

#' Prediction for gdpar_eb_fit
#'
#' Computes posterior predictions from the conditional HMC draws at
#' the plug-in EB estimate \eqn{\widehat{\theta}_{ref}^{EB}}. Supports
#' both in-sample and out-of-sample prediction via the \code{newdata}
#' argument; the multivariate \code{newdata} path is deferred to
#' Sub-phase 8.6.C together with the rest of \code{p > 1}.
#'
#' @param object A \code{gdpar_eb_fit} object.
#' @param newdata Optional data frame with the same variables as the
#'   training data. When \code{NULL} (default), in-sample predictions
#'   are returned.
#' @param type One of \code{"response"} (default; on the y scale via
#'   the family's inverse link) or \code{"linear_predictor"} (on the
#'   eta scale).
#' @param level Numeric scalar in (0, 1); credible-interval level.
#'   Defaults to 0.95.
#' @param ... Unused.
#' @return A list with components \code{mean}, \code{lower},
#'   \code{upper}, \code{draws} (matrix S x n).
#' @export
predict.gdpar_eb_fit <- function(object, newdata = NULL,
                                  type = c("response",
                                           "linear_predictor"),
                                  level = 0.95, ...) {
  type <- match.arg(type)
  if (!is.null(newdata)) {
    gdpar_abort(
      paste(
        "predict.gdpar_eb_fit(newdata = ...) is canonized for",
        "Sub-phase 8.6.C; in 8.6.B only in-sample predictions",
        "(newdata = NULL) are supported."
      ),
      class = "gdpar_unsupported_feature_error"
    )
  }
  fit_cond <- object$conditional_fit
  if (is.null(fit_cond)) {
    gdpar_abort(
      "Conditional HMC fit is missing from the gdpar_eb_fit object.",
      class = "gdpar_internal_error"
    )
  }
  require_suggested("posterior",
                    "extract posterior draws for prediction")
  draws <- fit_cond$draws()
  vars <- posterior::variables(draws)
  pick_var <- if (type == "linear_predictor") "^eta\\[" else "^y_pred\\["
  matches <- grep(pick_var, vars, value = TRUE)
  if (length(matches) == 0L) {
    gdpar_abort(
      sprintf(
        "Expected '%s' variables in the conditional fit draws, found none.",
        pick_var
      ),
      class = "gdpar_internal_error"
    )
  }
  mat <- posterior::as_draws_matrix(
    posterior::subset_draws(draws, variable = matches)
  )
  alpha <- 1 - level
  qs <- apply(mat, 2L, stats::quantile,
              probs = c(alpha / 2, 1 - alpha / 2),
              names = FALSE)
  list(
    mean  = colMeans(mat),
    lower = qs[1L, ],
    upper = qs[2L, ],
    draws = mat,
    level = level,
    type  = type
  )
}

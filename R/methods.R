#' Print method for gdpar_fit objects
#'
#' Concise summary of the fitted model: AMM specification, family,
#' anchor, sampler dimensions and convergence verdict.
#'
#' @param x An object of class \code{gdpar_fit}.
#' @param ... Unused; present for S3 generic compatibility.
#' @return Invisibly returns \code{x}.
#' @export
print.gdpar_fit <- function(x, ...) {
  cat("<gdpar_fit>\n")
  cat("  path                 : ", x$path, "\n", sep = "")
  cat("  family               : ", x$family$name,
      " (link = ", x$family$link, ")\n", sep = "")
  cat("  AMM Level            : ", x$amm$level, "\n", sep = "")
  multi <- !is.null(x$p) && x$p > 1L
  if (multi) {
    cat("  p (theta_ref dim)    : ", x$p, "\n", sep = "")
    cat("  anchor               : [",
        paste(format(x$anchor, digits = 4), collapse = ", "),
        "]\n", sep = "")
  } else {
    cat("  anchor               : ", format(x$anchor, digits = 4),
        "\n", sep = "")
  }
  n_obs <- if (multi) {
    if (length(x$design$Z_a_list) > 0L &&
        nrow(x$design$Z_a_list[[1L]]) > 0L) {
      nrow(x$design$Z_a_list[[1L]])
    } else if (length(x$design$Z_b_list) > 0L &&
               nrow(x$design$Z_b_list[[1L]]) > 0L) {
      nrow(x$design$Z_b_list[[1L]])
    } else if (!is.null(x$design$X)) {
      nrow(x$design$X)
    } else {
      NA_integer_
    }
  } else {
    n_obs_tmp <- nrow(x$design$Z_a)
    if (is.null(n_obs_tmp) || n_obs_tmp == 0L) {
      n_obs_tmp <- nrow(x$design$Z_b)
    }
    if (is.null(n_obs_tmp) || n_obs_tmp == 0L) {
      n_obs_tmp <- nrow(x$design$X)
    }
    n_obs_tmp
  }
  cat("  observations         : ", n_obs, "\n", sep = "")
  if (!is.null(x$identifiability_report)) {
    cat("  identifiability_pass : ", x$identifiability_report$passed,
        "\n", sep = "")
  }
  cat("  converged            : ", x$diagnostics$converged,
      "\n", sep = "")
  if (!is.na(x$diagnostics$rhat_max)) {
    cat("  rhat_max             : ",
        format(x$diagnostics$rhat_max, digits = 4), "\n", sep = "")
    cat("  ess_bulk_min         : ",
        format(x$diagnostics$ess_bulk_min, digits = 4), "\n", sep = "")
    cat("  divergent_count      : ", x$diagnostics$divergent_count,
        "\n", sep = "")
  }
  invisible(x)
}

#' Summary method for gdpar_fit objects
#'
#' Returns the posterior summary table for the user-facing parameters
#' (theta_ref, hierarchical scales, family-specific dispersion).
#'
#' @param object An object of class \code{gdpar_fit}.
#' @param ... Unused; present for S3 generic compatibility.
#' @return A data frame of posterior summaries (mean, median, standard
#'   deviation, 5\% and 95\% quantiles, R-hat, ESS bulk and tail).
#'
#' @section Dependencies:
#' Posterior summaries are computed by \pkg{posterior}.
#' @export
summary.gdpar_fit <- function(object, ...) {
  assert_inherits(object, "gdpar_fit", "object")
  require_suggested("posterior", "summarize posterior draws")
  draws <- object$fit$draws()
  vars <- posterior::variables(draws)
  ignore_pat <- "^(eta|log_lik|y_pred|theta_i|a_coef|b_coef|a_raw|b_raw|W_raw)"
  keep_vars <- vars[!grepl(ignore_pat, vars)]
  if (length(keep_vars) == 0L) {
    keep_vars <- intersect(vars, "theta_ref")
    if (length(keep_vars) == 0L) keep_vars <- vars[1L]
  }
  do.call(
    posterior::summarise_draws,
    c(list(posterior::subset_draws(draws, variable = keep_vars)),
      posterior::default_summary_measures(),
      posterior::default_convergence_measures())
  )
}

#' Posterior draws of theta_i for a fitted gdpar model
#'
#' Returns the posterior draws of the individual parameters
#' \eqn{\theta_i = \theta_{\text{ref}} + \Delta(x_i, \theta_{\text{ref}})}
#' on the linear-predictor scale by default, or on the response scale
#' after applying the family's inverse link.
#'
#' @param object An object of class \code{gdpar_fit}.
#' @param newdata Optional data frame on which to compute predictions.
#'   When NULL (default), predictions are computed on the training
#'   data and returned from the Stan-side \code{theta_i} draws,
#'   without re-evaluation.
#' @param type Character scalar: \code{"theta_i"} (default; the linear
#'   predictor of the individual parameter), \code{"linear_predictor"}
#'   (synonym), or \code{"response"} (the inverse link of the linear
#'   predictor).
#' @param summary Character scalar: \code{"draws"} (default; full
#'   draws as a numeric matrix of shape draws-by-observations),
#'   \code{"mean_se"} (posterior mean and standard error per
#'   observation), or \code{"quantiles"} (posterior 5\%, 50\% and 95\%
#'   quantiles per observation).
#' @param ... Unused; present for S3 generic compatibility.
#'
#' @return The shape depends on \code{summary}: a numeric matrix when
#'   \code{summary = "draws"}; a data frame when \code{summary} is
#'   \code{"mean_se"} or \code{"quantiles"}.
#'
#' @details
#' Predictions on \code{newdata} require evaluating the centered
#' bases at the new covariate values. The package uses the centering
#' parameters (column means of Z_a and Z_b, column means and standard
#' deviations of X) recorded at fit time so that the transformation is
#' identical to that applied during training.
#'
#' For \code{type = "response"}, the inverse link is applied at the
#' draw level. The reported quantiles are therefore the quantiles of
#' the response-scale distribution and are not the inverse link of the
#' linear-predictor quantiles unless the link is the identity.
#'
#' @section Dependencies:
#' Uses \pkg{posterior} to extract draws.
#' @export
predict.gdpar_fit <- function(object, newdata = NULL,
                              type = c("theta_i", "linear_predictor",
                                       "response"),
                              summary = c("draws", "mean_se",
                                          "quantiles"),
                              ...) {
  type <- match.arg(type)
  summary <- match.arg(summary)
  assert_inherits(object, "gdpar_fit", "object")
  require_suggested("posterior", "extract posterior draws")
  K_val <- object[["K"]]
  p_val <- object[["p"]]
  is_K_individual <- !is.null(K_val) && K_val > 1L
  multi <- !is.null(p_val) && p_val > 1L
  draws <- object$fit$draws()
  if (is_K_individual) {
    return(predict_gdpar_fit_K(object, newdata, type, summary, draws))
  }
  if (multi) {
    return(predict_gdpar_fit_multi(object, newdata, type, summary, draws))
  }
  if (is.null(newdata)) {
    theta_draws <- posterior::subset_draws(draws, variable = "theta_i")
    mat <- posterior::as_draws_matrix(theta_draws)
  } else if (!is.null(object$group_info)) {
    mat <- predict_from_newdata_grouped(object, newdata, draws)
  } else {
    mat <- predict_from_newdata(object, newdata, draws)
  }
  if (type == "response") {
    mat <- object$family$inv_link(mat)
  }
  if (summary == "draws") {
    return(mat)
  }
  if (summary == "mean_se") {
    return(data.frame(
      mean = colMeans(mat),
      se   = apply(mat, 2L, stats::sd)
    ))
  }
  q <- apply(mat, 2L, stats::quantile, probs = c(0.05, 0.5, 0.95),
              names = FALSE)
  data.frame(q05 = q[1L, ], q50 = q[2L, ], q95 = q[3L, ])
}

#' Internal: multivariate path of predict.gdpar_fit
#'
#' Returns predictions for fits with \code{p > 1}. When \code{newdata}
#' is NULL the function reads the generated quantity \code{theta_i}
#' (\code{matrix[n, p]} in the multivariate Stan template) and reshapes
#' it to an array of dimensions \code{(S, n, p)}. When \code{newdata}
#' is supplied the function rebuilds \eqn{\eta_{i,k}} per coordinate
#' from the posterior draws of \code{theta_ref}, \code{a_coef},
#' \code{c_b}, \code{W_raw} and \code{sigma_W} via
#' \code{predict_from_newdata_multi}.
#'
#' @keywords internal
#' @noRd
predict_gdpar_fit_multi <- function(object, newdata, type, summary,
                                    draws) {
  p <- object$p
  if (is.null(newdata)) {
    theta_draws <- posterior::subset_draws(draws, variable = "theta_i")
    mat <- posterior::as_draws_matrix(theta_draws)
    var_names <- posterior::variables(theta_draws)
    matches <- regmatches(
      var_names,
      regexec("^theta_i\\[(\\d+),(\\d+)\\]$", var_names)
    )
    idx_i <- as.integer(vapply(matches, function(m) m[2L], character(1L)))
    idx_k <- as.integer(vapply(matches, function(m) m[3L], character(1L)))
    if (any(is.na(idx_i)) || any(is.na(idx_k))) {
      gdpar_abort(
        "Internal error: failed to parse theta_i variable names from draws.",
        class = "gdpar_internal_error",
        data = list(var_names = utils::head(var_names, 5L))
      )
    }
    n <- max(idx_i)
    S <- nrow(mat)
    arr <- array(NA_real_, dim = c(S, n, p))
    for (col in seq_len(ncol(mat))) {
      arr[, idx_i[col], idx_k[col]] <- mat[, col]
    }
    row_names <- NULL
  } else if (!is.null(object$group_info)) {
    arr <- predict_from_newdata_grouped_multi(object, newdata, draws)
    row_names <- rownames(newdata)
    n <- dim(arr)[2L]
    S <- dim(arr)[1L]
  } else {
    arr <- predict_from_newdata_multi(object, newdata, draws)
    row_names <- rownames(newdata)
    n <- dim(arr)[2L]
    S <- dim(arr)[1L]
  }
  if (type == "response") {
    fam_k_inv <- lapply(object$family$families, function(f) f$inv_link)
    for (k in seq_len(p)) {
      arr[, , k] <- fam_k_inv[[k]](arr[, , k])
    }
  }
  dimnames(arr) <- list(NULL, row_names, paste0("dim_", seq_len(p)))
  if (summary == "draws") {
    return(arr)
  }
  if (summary == "mean_se") {
    out <- lapply(seq_len(p), function(k) {
      mat_k <- arr[, , k, drop = TRUE]
      data.frame(
        mean = colMeans(mat_k),
        se   = apply(mat_k, 2L, stats::sd)
      )
    })
    names(out) <- paste0("dim_", seq_len(p))
    return(out)
  }
  out <- lapply(seq_len(p), function(k) {
    mat_k <- arr[, , k, drop = TRUE]
    q <- apply(mat_k, 2L, stats::quantile, probs = c(0.05, 0.5, 0.95),
                names = FALSE)
    data.frame(q05 = q[1L, ], q50 = q[2L, ], q95 = q[3L, ])
  })
  names(out) <- paste0("dim_", seq_len(p))
  out
}

#' Internal: K-individual path of predict.gdpar_fit (Sub-phase 8.3.9 Ruta B)
#'
#' Returns predictions for fits with \code{K > 1} (distributional
#' regression on the K-individual Stan template
#' \code{amm_distrib_K.stan}). The Stan-side generated quantity is
#' \code{theta_i_k} (\code{matrix[n, K]}), not \code{theta_i}; this
#' helper parses the \code{theta_i_k[i, k]} variable names and
#' reshapes the draws into an array of dimensions \code{(S, n, K)}.
#'
#' For \code{type = "response"}, each slot's canonical inverse link
#' is applied via \code{object$family$param_specs[[k]]$inv_link}.
#' Slot k = 1 is the location (uses the family's primary link); slot
#' k > 1 are dispersion / shape / mixture-weight parameters with
#' their own canonical link (log for sigma, log for phi, logit for
#' pi, etc.). Applying the location link to all slots indiscriminately
#' would be incorrect.
#'
#' When \code{newdata} is non-NULL the function dispatches to
#' \code{predict_from_newdata_K}, which mirrors
#' \code{predict_from_newdata_multi} but iterates the K slots: for
#' each slot k the per-slot \eqn{Z_{a,k}}/\eqn{Z_{b,k}} design
#' matrices are rebuilt from \code{amm_list_canonical[[k]]} and the
#' per-slot \eqn{\eta_{i,k}} is composed via \code{theta_ref_k[k]},
#' \code{a_coef_k[k]}, \code{c_b_k[k]}, and the globally-shared
#' \code{W_raw} / \code{sigma_W} (with the basis-difference vector
#' evaluated at the slot-specific \code{theta_ref_k[k]} and
#' \code{theta_anchor_K[k]}). K > 1 with grouping (\code{J_groups > 1})
#' and K > 1 with B-spline W bases on \code{newdata} are queued for
#' Session 8.4 and raise \code{gdpar_unsupported_feature_error}.
#'
#' @param object A \code{gdpar_fit} object with \code{K > 1}.
#' @param newdata Either \code{NULL} (in-sample prediction reading the
#'   generated quantity \code{theta_i_k}) or a data frame with the
#'   covariates required by every slot's AMM specification.
#' @param type One of \code{"theta_i"}, \code{"linear_predictor"},
#'   \code{"response"}.
#' @param summary One of \code{"draws"}, \code{"mean_se"},
#'   \code{"quantiles"}.
#' @param draws Posterior draws array.
#' @return Depends on \code{summary}: array \code{(S, n, K)} when
#'   \code{summary = "draws"}; named list of length \code{K} with
#'   per-slot summary data frames otherwise.
#' @keywords internal
#' @noRd
predict_gdpar_fit_K <- function(object, newdata, type, summary, draws) {
  require_suggested("posterior", "extract posterior draws")
  K <- as.integer(object[["K"]])
  slot_names <- object[["slot_names"]]
  if (is.null(slot_names) || length(slot_names) != K) {
    slot_names <- paste0("slot_", seq_len(K))
  }
  if (is.null(newdata)) {
    theta_draws <- posterior::subset_draws(draws, variable = "theta_i_k")
    mat <- unclass(posterior::as_draws_matrix(theta_draws))
    var_names <- colnames(mat)
    matches <- regmatches(
      var_names,
      regexec("^theta_i_k\\[(\\d+),(\\d+)\\]$", var_names)
    )
    idx_i <- as.integer(vapply(matches, function(m) m[2L], character(1L)))
    idx_k <- as.integer(vapply(matches, function(m) m[3L], character(1L)))
    if (any(is.na(idx_i)) || any(is.na(idx_k))) {
      gdpar_abort(
        "Internal error: failed to parse theta_i_k variable names.",
        class = "gdpar_internal_error",
        data = list(var_names = utils::head(var_names, 5L))
      )
    }
    n <- max(idx_i)
    S <- nrow(mat)
    arr <- array(NA_real_, dim = c(S, n, K))
    for (col in seq_len(ncol(mat))) {
      arr[, idx_i[col], idx_k[col]] <- mat[, col]
    }
  } else {
    arr <- predict_from_newdata_K(object, newdata, draws)
  }
  dimnames(arr) <- list(NULL, NULL, slot_names)
  if (type == "response") {
    param_specs <- object$family$param_specs
    if (is.null(param_specs) || length(param_specs) != K) {
      gdpar_abort(
        paste(
          "Internal error: K-individual fit lacks per-slot param_specs",
          "for type = 'response'. Use type = 'linear_predictor' instead."
        ),
        class = "gdpar_internal_error"
      )
    }
    for (k in seq_len(K)) {
      arr[, , k] <- param_specs[[k]]$inv_link(arr[, , k])
    }
  }
  if (summary == "draws") {
    return(arr)
  }
  if (summary == "mean_se") {
    out <- lapply(seq_len(K), function(k) {
      mat_k <- arr[, , k, drop = TRUE]
      data.frame(
        mean = colMeans(mat_k),
        se   = apply(mat_k, 2L, stats::sd)
      )
    })
    names(out) <- slot_names
    return(out)
  }
  out <- lapply(seq_len(K), function(k) {
    mat_k <- arr[, , k, drop = TRUE]
    q <- apply(mat_k, 2L, stats::quantile, probs = c(0.05, 0.5, 0.95),
                names = FALSE)
    data.frame(q05 = q[1L, ], q50 = q[2L, ], q95 = q[3L, ])
  })
  names(out) <- slot_names
  out
}

#' Internal: compute theta_i posterior draws on new data
#'
#' Reconstructs the linear predictor on a new data frame using the
#' centering parameters recorded at fit time, the posterior draws of
#' theta_ref, the basis coefficients a_coef and b_coef, and the
#' modulating coefficients W_raw.
#'
#' @param object A \code{gdpar_fit} object.
#' @param newdata Data frame with the covariates required by the AMM
#'   specification.
#' @param draws Posterior draws array.
#' @return A numeric matrix of draws-by-newdata-rows.
#' @keywords internal
#' @noRd
predict_from_newdata <- function(object, newdata, draws) {
  require_suggested("posterior", "extract posterior draws")
  amm <- object$amm
  design_train <- object$design
  new_design <- list()
  if (!is.null(amm$a)) {
    Za_full <- stats::model.matrix(stats::update(amm$a, ~ . + 0), newdata)
    Za_full <- as.matrix(Za_full)
    new_design$Z_a <- sweep(Za_full, 2L, design_train$Z_a_means, "-")
  } else {
    new_design$Z_a <- matrix(numeric(0), nrow = nrow(newdata), ncol = 0)
  }
  if (!is.null(amm$b)) {
    Zb_full <- stats::model.matrix(stats::update(amm$b, ~ . + 0), newdata)
    Zb_full <- as.matrix(Zb_full)
    new_design$Z_b <- sweep(Zb_full, 2L, design_train$Z_b_means, "-")
  } else {
    new_design$Z_b <- matrix(numeric(0), nrow = nrow(newdata), ncol = 0)
  }
  if (!is.null(amm$W) && length(design_train$X_names) > 0L) {
    X_full <- as.matrix(newdata[, design_train$X_names, drop = FALSE])
    storage.mode(X_full) <- "double"
    X_centered <- sweep(X_full, 2L, design_train$X_means, "-")
    new_design$X <- sweep(X_centered, 2L, design_train$X_sds, "/")
  } else {
    new_design$X <- matrix(numeric(0), nrow = nrow(newdata), ncol = 0)
  }
  theta_ref_dr <- unclass(posterior::as_draws_matrix(
    posterior::subset_draws(draws, variable = "theta_ref")
  ))
  S <- nrow(theta_ref_dr)
  n_new <- nrow(newdata)
  out <- matrix(NA_real_, nrow = S, ncol = n_new)
  has_a <- ncol(new_design$Z_a) > 0L
  has_b <- ncol(new_design$Z_b) > 0L
  has_W <- ncol(new_design$X) > 0L &&
    !is.null(amm$W) && amm$W$type == "polynomial"
  a_coef_dr <- if (has_a) {
    unclass(posterior::as_draws_matrix(
      posterior::subset_draws(draws, variable = "a_coef")
    ))
  } else {
    NULL
  }
  b_coef_dr <- if (has_b) {
    unclass(posterior::as_draws_matrix(
      posterior::subset_draws(draws, variable = "b_coef")
    ))
  } else {
    NULL
  }
  W_raw_dr <- if (has_W) {
    unclass(posterior::as_draws_matrix(
      posterior::subset_draws(draws, variable = "W_raw")
    ))
  } else {
    NULL
  }
  sigma_W_dr <- if (has_W) {
    unclass(posterior::as_draws_matrix(
      posterior::subset_draws(draws, variable = "sigma_W")
    ))
  } else {
    NULL
  }
  d <- ncol(new_design$X)
  dim_W <- if (has_W) amm$W$dim else 0L
  theta_anchor <- object$anchor
  for (s in seq_len(S)) {
    eta <- rep(theta_ref_dr[s, 1L], n_new)
    if (has_a) {
      eta <- eta + as.numeric(new_design$Z_a %*% a_coef_dr[s, ])
    }
    if (has_b) {
      eta <- eta + as.numeric(new_design$Z_b %*% b_coef_dr[s, ]) *
        theta_ref_dr[s, 1L]
    }
    if (has_W) {
      W_mat <- matrix(W_raw_dr[s, ], nrow = dim_W, ncol = d)
      basis_diff <- vapply(seq_len(dim_W), function(k) {
        theta_ref_dr[s, 1L]^k - theta_anchor^k
      }, numeric(1L))
      W_diff_x <- as.numeric(crossprod(W_mat, basis_diff)) *
        sigma_W_dr[s, 1L]
      eta <- eta + as.numeric(new_design$X %*% W_diff_x)
    }
    out[s, ] <- eta
  }
  out
}

#' Internal: compute theta_i posterior draws on new data for p > 1
#'
#' Multivariate counterpart of \code{predict_from_newdata}. Rebuilds
#' \eqn{\eta_{i,k}} per coordinate \eqn{k = 1, \dots, p} from posterior
#' draws of \code{theta_ref}, \code{a_coef}, \code{c_b}, \code{W_raw}
#' and \code{sigma_W}, applying the per-coordinate centering recorded
#' at fit time and the shared standardization of the modulating
#' covariates.
#'
#' The reconstruction mirrors exactly the formula encoded in the
#' multivariate Stan template (\code{inst/stan/amm_distrib_multi.stan}):
#' \deqn{\eta_{i,k} = \theta_{\mathrm{ref},k} + Z_{a,k}[i,]\cdot a_{\mathrm{coef},k}
#'        + Z_{b,k}[i,]\cdot c_{b,k}
#'        + \sum_{j=1}^{\dim_W/p}
#'         (\theta_{\mathrm{ref},k}^j - \theta_{\mathrm{anchor},k}^j)\,
#'         W_{\mathrm{raw}}[r_{k,j},\cdot]\,\sigma_W\,X[i,]^{\top}}
#' with \eqn{r_{k,j} = (k-1) W_{\mathrm{per\_k\_dim}} + j} and the
#' \eqn{\sigma_W} multiplier present only when the modulating component
#' was sampled in the non-centered parametrization.
#'
#' @param object A \code{gdpar_fit} object with \code{p > 1}.
#' @param newdata Data frame with the covariates required by the AMM
#'   specification.
#' @param draws Posterior draws array.
#' @return A three-dimensional array of dimensions \code{(S, n_new, p)}.
#' @keywords internal
#' @noRd
predict_from_newdata_multi <- function(object, newdata, draws) {
  require_suggested("posterior", "extract posterior draws")
  amm <- object$amm
  design_train <- object$design
  p <- object$p
  n_new <- nrow(newdata)
  Z_a_list_new <- vector("list", p)
  Z_b_list_new <- vector("list", p)
  for (k in seq_len(p)) {
    a_k <- amm$dims[[k]]$a
    b_k <- amm$dims[[k]]$b
    if (!is.null(a_k)) {
      Z_full <- stats::model.matrix(stats::update(a_k, ~ . + 0), newdata)
      Z_full <- as.matrix(Z_full)
      Z_a_list_new[[k]] <- sweep(
        Z_full, 2L, design_train$Z_a_means_list[[k]], "-"
      )
    } else {
      Z_a_list_new[[k]] <- matrix(numeric(0), nrow = n_new, ncol = 0L)
    }
    if (!is.null(b_k)) {
      Z_full <- stats::model.matrix(stats::update(b_k, ~ . + 0), newdata)
      Z_full <- as.matrix(Z_full)
      Z_b_list_new[[k]] <- sweep(
        Z_full, 2L, design_train$Z_b_means_list[[k]], "-"
      )
    } else {
      Z_b_list_new[[k]] <- matrix(numeric(0), nrow = n_new, ncol = 0L)
    }
  }
  has_W <- !is.null(amm$W) && length(design_train$X_names) > 0L &&
    amm$W$type == "polynomial"
  if (has_W) {
    X_full <- as.matrix(newdata[, design_train$X_names, drop = FALSE])
    storage.mode(X_full) <- "double"
    X_centered <- sweep(X_full, 2L, design_train$X_means, "-")
    X_new <- sweep(X_centered, 2L, design_train$X_sds, "/")
  } else {
    X_new <- matrix(numeric(0), nrow = n_new, ncol = 0L)
  }
  theta_ref_dr <- .extract_theta_ref_multi_flat(draws, J_groups = 1L,
                                                p = p)
  S <- nrow(theta_ref_dr)
  has_a_any <- any(vapply(Z_a_list_new, function(z) ncol(z) > 0L,
                          logical(1L)))
  has_b_any <- any(vapply(Z_b_list_new, function(z) ncol(z) > 0L,
                          logical(1L)))
  a_coef_dr_list <- if (has_a_any) {
    posterior_array_var_to_list(draws, "a_coef", p,
                                ncols = vapply(Z_a_list_new, ncol,
                                                integer(1L)))
  } else {
    NULL
  }
  c_b_dr_list <- if (has_b_any) {
    posterior_array_var_to_list(draws, "c_b", p,
                                ncols = vapply(Z_b_list_new, ncol,
                                                integer(1L)))
  } else {
    NULL
  }
  if (has_W) {
    dim_W <- amm$W$dim
    W_per_k_dim <- as.integer(dim_W / p)
    d <- ncol(X_new)
    W_raw_dr <- posterior_var_to_matrix_3d(draws, "W_raw",
                                           rows = dim_W, cols = d)
    cp_W <- isTRUE(object$parametrization$cp_W)
    sigma_W_dr <- if (cp_W) {
      NULL
    } else {
      posterior_var_to_array_2d(draws, "sigma_W", 1L)
    }
  } else {
    dim_W <- 0L
    W_per_k_dim <- 0L
    d <- 0L
    W_raw_dr <- NULL
    sigma_W_dr <- NULL
    cp_W <- FALSE
  }
  theta_anchor <- object$anchor
  out <- array(NA_real_, dim = c(S, n_new, p))
  for (s in seq_len(S)) {
    for (k in seq_len(p)) {
      eta_k <- rep(theta_ref_dr[s, k], n_new)
      if (!is.null(a_coef_dr_list) && ncol(Z_a_list_new[[k]]) > 0L) {
        eta_k <- eta_k +
          as.numeric(Z_a_list_new[[k]] %*% a_coef_dr_list[[k]][s, ])
      }
      if (!is.null(c_b_dr_list) && ncol(Z_b_list_new[[k]]) > 0L) {
        eta_k <- eta_k +
          as.numeric(Z_b_list_new[[k]] %*% c_b_dr_list[[k]][s, ])
      }
      if (has_W && d > 0L && W_per_k_dim > 0L) {
        row_start <- (k - 1L) * W_per_k_dim + 1L
        row_end <- k * W_per_k_dim
        W_block <- W_raw_dr[s, row_start:row_end, , drop = FALSE]
        dim(W_block) <- c(W_per_k_dim, d)
        basis_diff <- vapply(seq_len(W_per_k_dim), function(jj) {
          theta_ref_dr[s, k]^jj - theta_anchor[k]^jj
        }, numeric(1L))
        W_diff_x_k <- as.numeric(crossprod(W_block, basis_diff))
        if (!cp_W) {
          W_diff_x_k <- W_diff_x_k * sigma_W_dr[s, 1L]
        }
        eta_k <- eta_k + as.numeric(X_new %*% W_diff_x_k)
      }
      out[s, , k] <- eta_k
    }
  }
  out
}

#' Internal: compute theta_i_k posterior draws on new data for K > 1
#'
#' K-individual counterpart of \code{predict_from_newdata} /
#' \code{predict_from_newdata_multi}. Rebuilds \eqn{\eta_{i,k}} for
#' \eqn{k = 1, \dots, K} from posterior draws of \code{theta_ref_k},
#' \code{a_coef_k}, \code{c_b_k}, the globally-shared \code{W_raw},
#' and \code{sigma_W}, mirroring the formula encoded in the
#' K-individual Stan template (\code{inst/stan/amm_distrib_K.stan},
#' lines 483-526):
#' \deqn{\eta_{i,k} = \theta_{\mathrm{ref},k} + Z_{a,k}[i,]\cdot a_{\mathrm{coef},k}
#'        + Z_{b,k}[i,]\cdot c_{b,k}
#'        + \sum_{j=1}^{\dim_W}
#'         (\theta_{\mathrm{ref},k}^j - \theta_{\mathrm{anchor},k}^j)\,
#'         W_{\mathrm{raw}}[j,\cdot]\,\sigma_W\,X[i,]^{\top}}
#' Centering of \eqn{Z_{a,k}}, \eqn{Z_{b,k}} and \eqn{X} uses the
#' per-slot means and global standard deviations recorded at fit time
#' in \code{object$design_K}.
#'
#' Bit-exact compatibility with the in-sample \code{theta_i_k}
#' generated quantity is preserved only when \code{newdata} matches
#' the training data row-for-row (same covariates, same anchoring).
#' Otherwise the function returns the posterior of the linear
#' predictor evaluated on the user's new covariates.
#'
#' B-spline W bases on new data are queued for Session 8.4 and aborted
#' here with \code{gdpar_unsupported_feature_error} (the in-sample
#' path supports both polynomial and B-spline via the helper
#' \code{apply_W_basis_diff} in Stan; the R-side reconstruction
#' presently mirrors only the polynomial branch).
#'
#' Grouping (\code{J_groups > 1}) is also queued for Session 8.4.
#'
#' @param object A \code{gdpar_fit} object with \code{K > 1}.
#' @param newdata Data frame with the covariates required by every
#'   slot's AMM specification.
#' @param draws Posterior draws array.
#' @return A three-dimensional array of dimensions \code{(S, n_new, K)}.
#' @keywords internal
#' @noRd
predict_from_newdata_K <- function(object, newdata, draws) {
  require_suggested("posterior", "extract posterior draws")
  K <- as.integer(object[["K"]])
  if (!is.null(object$group_info)) {
    gdpar_abort(
      paste(
        "predict() on newdata for K-individual fits with grouping",
        "(K > 1 and J_groups > 1) is not yet implemented and is",
        "queued for Session 8.4."
      ),
      class = "gdpar_unsupported_feature_error",
      data = list(K = K,
                  J_groups = length(object$group_info$levels))
    )
  }
  amm_list <- object$amm_list_canonical
  design_K <- object$design_K
  W_active_per_slot <- vapply(amm_list, function(a) {
    !is.null(a$W)
  }, logical(1L))
  has_W <- any(W_active_per_slot) && length(design_K$X_names) > 0L
  if (has_W) {
    first_active <- which(W_active_per_slot)[1L]
    W_spec <- amm_list[[first_active]]$W
    if (!identical(W_spec$type, "polynomial")) {
      gdpar_abort(
        paste(
          "predict() on newdata for K-individual fits with a",
          "non-polynomial W basis (e.g., B-spline) is not yet",
          "implemented and is queued for Session 8.4. The in-sample",
          "path (newdata = NULL) supports both polynomial and",
          "B-spline via the Stan-side apply_W_basis_diff() helper."
        ),
        class = "gdpar_unsupported_feature_error",
        data = list(K = K, W_type = W_spec$type)
      )
    }
  }
  n_new <- nrow(newdata)
  Z_a_k_new <- vector("list", K)
  Z_b_k_new <- vector("list", K)
  for (k in seq_len(K)) {
    a_k <- amm_list[[k]]$a
    b_k <- amm_list[[k]]$b
    if (!is.null(a_k)) {
      Z_full <- stats::model.matrix(stats::update(a_k, ~ . + 0), newdata)
      Z_full <- as.matrix(Z_full)
      Z_a_k_new[[k]] <- sweep(
        Z_full, 2L, design_K$Z_a_k_means_list[[k]], "-"
      )
    } else {
      Z_a_k_new[[k]] <- matrix(numeric(0), nrow = n_new, ncol = 0L)
    }
    if (!is.null(b_k)) {
      Z_full <- stats::model.matrix(stats::update(b_k, ~ . + 0), newdata)
      Z_full <- as.matrix(Z_full)
      Z_b_k_new[[k]] <- sweep(
        Z_full, 2L, design_K$Z_b_k_means_list[[k]], "-"
      )
    } else {
      Z_b_k_new[[k]] <- matrix(numeric(0), nrow = n_new, ncol = 0L)
    }
  }
  if (has_W) {
    X_full <- as.matrix(newdata[, design_K$X_names, drop = FALSE])
    storage.mode(X_full) <- "double"
    X_centered <- sweep(X_full, 2L, design_K$X_means, "-")
    X_new <- sweep(X_centered, 2L, design_K$X_sds, "/")
  } else {
    X_new <- matrix(numeric(0), nrow = n_new, ncol = 0L)
  }
  theta_ref_arr <- posterior_var_to_matrix_3d(draws, "theta_ref_k",
                                              rows = 1L, cols = K)
  S <- dim(theta_ref_arr)[1L]
  theta_ref_dr <- matrix(theta_ref_arr[, 1L, ], nrow = S, ncol = K)
  ncols_a <- vapply(seq_len(K), function(k) ncol(Z_a_k_new[[k]]),
                    integer(1L))
  ncols_b <- vapply(seq_len(K), function(k) ncol(Z_b_k_new[[k]]),
                    integer(1L))
  has_a_any <- any(ncols_a > 0L)
  has_b_any <- any(ncols_b > 0L)
  a_dr_list <- if (has_a_any) {
    posterior_array_var_to_list(draws, "a_coef_k", K, ncols = ncols_a)
  } else {
    NULL
  }
  c_b_dr_list <- if (has_b_any) {
    posterior_array_var_to_list(draws, "c_b_k", K, ncols = ncols_b)
  } else {
    NULL
  }
  if (has_W) {
    dim_W <- as.integer(W_spec$dim)
    d <- ncol(X_new)
    W_raw_dr <- posterior_var_to_matrix_3d(draws, "W_raw",
                                           rows = dim_W, cols = d)
    cp_W <- isTRUE(object$parametrization$cp_W)
    sigma_W_dr <- if (cp_W) {
      NULL
    } else {
      posterior_var_to_array_2d(draws, "sigma_W", 1L)
    }
  } else {
    dim_W <- 0L
    d <- 0L
    W_raw_dr <- NULL
    sigma_W_dr <- NULL
    cp_W <- FALSE
  }
  theta_anchor <- object$anchor
  if (length(theta_anchor) < K) {
    theta_anchor <- rep_len(theta_anchor, K)
  }
  out <- array(NA_real_, dim = c(S, n_new, K))
  for (s in seq_len(S)) {
    for (k in seq_len(K)) {
      eta_k <- rep(theta_ref_dr[s, k], n_new)
      if (!is.null(a_dr_list) && ncols_a[k] > 0L) {
        eta_k <- eta_k +
          as.numeric(Z_a_k_new[[k]] %*% a_dr_list[[k]][s, ])
      }
      if (!is.null(c_b_dr_list) && ncols_b[k] > 0L) {
        eta_k <- eta_k +
          as.numeric(Z_b_k_new[[k]] %*% c_b_dr_list[[k]][s, ])
      }
      if (has_W && d > 0L && dim_W > 0L) {
        basis_diff_k <- vapply(seq_len(dim_W), function(jj) {
          theta_ref_dr[s, k]^jj - theta_anchor[k]^jj
        }, numeric(1L))
        W_diff_x_k <- as.numeric(crossprod(W_raw_dr[s, , ],
                                            basis_diff_k))
        if (!cp_W) {
          W_diff_x_k <- W_diff_x_k * sigma_W_dr[s, 1L]
        }
        eta_k <- eta_k + as.numeric(X_new %*% W_diff_x_k)
      }
      out[s, , k] <- eta_k
    }
  }
  out
}

#' Internal: extract a one-dimensional Stan variable as an (S, p) matrix
#'
#' Variables of the form \code{var[k]} (e.g., \code{theta_ref[k]},
#' \code{sigma_W[1]}) are parsed and reshaped into a matrix with
#' \code{S} rows (draws) and \code{p} columns (indexed lexicographically
#' by the position parsed from the variable name). Robust against
#' implementation-specific ordering in \pkg{posterior}.
#'
#' @keywords internal
#' @noRd
posterior_var_to_array_2d <- function(draws, var_name, p) {
  sub <- posterior::subset_draws(draws, variable = var_name)
  mat <- unclass(posterior::as_draws_matrix(sub))
  var_names <- posterior::variables(sub)
  pat <- sprintf("^%s\\[(\\d+)\\]$", var_name)
  matches <- regmatches(var_names, regexec(pat, var_names))
  idx <- as.integer(vapply(matches, function(m) m[2L], character(1L)))
  if (any(is.na(idx))) {
    gdpar_abort(
      sprintf(
        "Internal error: failed to parse '%s' variable names from draws.",
        var_name
      ),
      class = "gdpar_internal_error"
    )
  }
  if (length(idx) != p) {
    gdpar_abort(
      sprintf(
        "Internal error: expected %d entries for '%s', found %d.",
        p, var_name, length(idx)
      ),
      class = "gdpar_internal_error"
    )
  }
  S <- nrow(mat)
  out <- matrix(NA_real_, nrow = S, ncol = p)
  for (col in seq_len(ncol(mat))) {
    out[, idx[col]] <- mat[, col]
  }
  out
}

#' Internal: extract a two-dimensional Stan variable as an (S, rows, cols) array
#'
#' For matrices declared in Stan (e.g., \code{matrix[dim_W, d] W_raw}),
#' parse the \code{var[r, c]} variable names and reshape into an array
#' of dimensions \code{(S, rows, cols)}.
#'
#' @keywords internal
#' @noRd
posterior_var_to_matrix_3d <- function(draws, var_name, rows, cols) {
  sub <- posterior::subset_draws(draws, variable = var_name)
  mat <- unclass(posterior::as_draws_matrix(sub))
  var_names <- posterior::variables(sub)
  pat <- sprintf("^%s\\[(\\d+),(\\d+)\\]$", var_name)
  matches <- regmatches(var_names, regexec(pat, var_names))
  idx_r <- as.integer(vapply(matches, function(m) m[2L], character(1L)))
  idx_c <- as.integer(vapply(matches, function(m) m[3L], character(1L)))
  if (any(is.na(idx_r)) || any(is.na(idx_c))) {
    gdpar_abort(
      sprintf(
        "Internal error: failed to parse '%s' variable names from draws.",
        var_name
      ),
      class = "gdpar_internal_error"
    )
  }
  S <- nrow(mat)
  out <- array(NA_real_, dim = c(S, rows, cols))
  for (col in seq_len(ncol(mat))) {
    out[, idx_r[col], idx_c[col]] <- mat[, col]
  }
  out
}

#' Internal: extract an array Stan variable into a per-k list of (S, ncols[k]) matrices
#'
#' For variables declared as \code{array[p] vector[J_max] var} (e.g.,
#' \code{a_coef[k][j]}, \code{c_b[k][j]}), parse the \code{var[k, j]}
#' variable names, restrict to the leading \code{ncols[k]} columns for
#' each coordinate (the trailing columns are padded zeros in the
#' template and must be skipped), and return a list of length \code{p}
#' whose \eqn{k}-th entry is an \code{(S, ncols[k])} matrix.
#'
#' @keywords internal
#' @noRd
posterior_array_var_to_list <- function(draws, var_name, p, ncols) {
  sub <- posterior::subset_draws(draws, variable = var_name)
  mat <- unclass(posterior::as_draws_matrix(sub))
  var_names <- posterior::variables(sub)
  pat <- sprintf("^%s\\[(\\d+),(\\d+)\\]$", var_name)
  matches <- regmatches(var_names, regexec(pat, var_names))
  idx_k <- as.integer(vapply(matches, function(m) m[2L], character(1L)))
  idx_j <- as.integer(vapply(matches, function(m) m[3L], character(1L)))
  if (any(is.na(idx_k)) || any(is.na(idx_j))) {
    gdpar_abort(
      sprintf(
        "Internal error: failed to parse '%s' variable names from draws.",
        var_name
      ),
      class = "gdpar_internal_error"
    )
  }
  S <- nrow(mat)
  out <- vector("list", p)
  for (k in seq_len(p)) {
    out[[k]] <- matrix(NA_real_, nrow = S, ncol = ncols[k])
  }
  for (col in seq_len(ncol(mat))) {
    k <- idx_k[col]
    j <- idx_j[col]
    if (j <= ncols[k]) {
      out[[k]][, j] <- mat[, col]
    }
  }
  out
}

#' Coefficients of a fitted gdpar model
#'
#' Returns posterior summaries (mean and 5\%/50\%/95\% quantiles) of
#' \code{theta_ref} and the basis coefficients (additive
#' \code{a_coef}, multiplicative \code{b_coef}/\code{c_b}, modulating
#' \code{W_raw}) as an object of class \code{gdpar_coef}. The format
#' is consistent across \code{p = 1} (scalar) and \code{p > 1} (multi)
#' fits.
#'
#' For multi fits, the modulating slot uses the effective
#' per-sample product \code{W_raw * sigma_W} (when \code{cp_W} is
#' FALSE) so that the coefficients are reported on the natural
#' modulating scale. The per-coordinate \code{b} slot draws from the
#' multi parameter \code{c_b} (already on the
#' theta_ref-multiplied scale; see Recovery 2, handoff 4).
#'
#' For K-individual fits (\code{K > 1}, distributional regression on
#' the per-slot template \code{amm_distrib_K.stan}) the function
#' returns a named list of length \code{K} whose entries are
#' \code{gdpar_coef} objects (each with \code{p = 1L}), one per slot
#' (decision E4.A of sub-phase 8.3.10). The modulating block
#' \code{W_raw} is globally shared across slots in the K-individual
#' template and its contribution enters \code{eta_k} for every slot
#' (see \code{amm_distrib_K.stan}, lines 500-523, and the canonical
#' decision "Scope of W: global" in the likelihood K > 1 decision of
#' handoff 28); accordingly, when any slot declared a non-NULL
#' \code{W} the resulting \code{W} component is attached to every
#' slot's \code{gdpar_coef} (and the per-sample product
#' \code{W_raw * sigma_W} is reported when \code{cp_W} is FALSE,
#' identically to the scalar path). K-individual fits with grouping
#' (\code{J_groups > 1}) are queued for Session 8.4 and raise
#' \code{gdpar_unsupported_feature_error}.
#'
#' @param object An object of class \code{gdpar_fit}.
#' @param ... Unused; present for S3 generic compatibility.
#' @return For scalar (\code{p = 1}) and multi (\code{p > 1}) fits an
#'   object of class \code{gdpar_coef}; for K-individual fits
#'   (\code{K > 1}) a named list of length \code{K} of \code{gdpar_coef}
#'   objects. Use \code{as.data.frame()} on a single \code{gdpar_coef}
#'   to obtain a long-tidy table; for K-individual fits use
#'   \code{lapply(coef(fit), as.data.frame)}.
#' @export
coef.gdpar_fit <- function(object, ...) {
  assert_inherits(object, "gdpar_fit", "object")
  require_suggested("posterior", "extract posterior draws")
  K_val <- object[["K"]]
  if (!is.null(K_val) && K_val > 1L) {
    return(coef_K_to_gdpar_coef_list(object))
  }
  p_val <- object[["p"]]
  p <- if (is.null(p_val)) 1L else as.integer(p_val)
  draws <- object$fit$draws()
  if (p == 1L) {
    coef_scalar_to_gdpar_coef(object, draws)
  } else {
    coef_multi_to_gdpar_coef(object, draws)
  }
}

#' Internal: build a named list of gdpar_coef from a K-individual fit
#'
#' Returns a named list of length \code{K} whose entries are
#' \code{gdpar_coef} objects (each with \code{p = 1L}), one per slot of
#' the K-individual distributional regression. Decision E4.A of sub-phase
#' 8.3.10: each slot reuses the existing S3 contract verbatim and the
#' top-level container is a plain named list (dispatch print/summary
#' on each entry via \code{lapply}).
#'
#' Stan-side variables consumed: \code{theta_ref_k[g, k]} (array per
#' \code{J_groups} of length-K vector; the K-individual path with
#' grouping is queued for Session 8.4 and aborts here),
#' \code{a_coef_k[k, j]} and \code{c_b_k[k, j]} (array[K] vector
#' padded to the maximum across slots; padding columns are silently
#' skipped via \code{posterior_array_var_to_list}), and the
#' globally-shared \code{W_raw[r, c]} plus \code{sigma_W[1]} (which
#' contribute to every slot's \code{eta_k} in the K-individual
#' template, lines 500-523 of \code{amm_distrib_K.stan}; the
#' resulting \code{W} block is therefore replicated across all
#' slots' \code{gdpar_coef} when any slot declared a non-NULL
#' \code{W}).
#'
#' The \code{b} slot is built from the linear-reparametrization variable
#' \code{c_b_k} (not \code{b_coef_k}), exactly as in the multivariate
#' path of \code{coef.gdpar_fit} for \code{p > 1} (Recovery 2, handoff 4
#' addendum).
#'
#' @param object A \code{gdpar_fit} object with \code{K > 1}.
#' @return A named list of \code{gdpar_coef} objects of length \code{K}.
#' @keywords internal
#' @noRd
coef_K_to_gdpar_coef_list <- function(object) {
  require_suggested("posterior", "extract posterior draws")
  K <- as.integer(object[["K"]])
  if (!is.null(object$group_info)) {
    gdpar_abort(
      paste(
        "coef() on K-individual fits with grouping",
        "(K > 1 and J_groups > 1) is not yet implemented and is",
        "queued for Session 8.4. The K-individual Stan template",
        "supports per-group theta_ref_k, but the R-side per-slot",
        "extraction needs to thread group_levels through every gdpar_coef",
        "entry; this requires an extended contract on gdpar_coef and is",
        "deferred until the Session 8.4 audit."
      ),
      class = "gdpar_unsupported_feature_error",
      data = list(K = K,
                  J_groups = length(object$group_info$levels))
    )
  }
  design_K <- object$design_K
  amm_list <- object$amm_list_canonical
  slot_names <- object[["slot_names"]]
  if (is.null(slot_names) || length(slot_names) != K) {
    slot_names <- paste0("slot_", seq_len(K))
  }
  draws <- object$fit$draws()
  theta_ref_arr <- posterior_var_to_matrix_3d(draws, "theta_ref_k",
                                              rows = 1L, cols = K)
  S <- dim(theta_ref_arr)[1L]
  theta_ref_mat <- matrix(theta_ref_arr[, 1L, ], nrow = S, ncol = K)
  ncols_a <- vapply(seq_len(K), function(k) {
    length(design_K$Z_a_k_names_list[[k]])
  }, integer(1L))
  ncols_b <- vapply(seq_len(K), function(k) {
    length(design_K$Z_b_k_names_list[[k]])
  }, integer(1L))
  has_a_any <- any(ncols_a > 0L)
  has_b_any <- any(ncols_b > 0L)
  a_dr_list <- if (has_a_any) {
    posterior_array_var_to_list(draws, "a_coef_k", K, ncols = ncols_a)
  } else {
    NULL
  }
  b_dr_list <- if (has_b_any) {
    posterior_array_var_to_list(draws, "c_b_k", K, ncols = ncols_b)
  } else {
    NULL
  }
  W_active_per_slot <- vapply(amm_list, function(a) {
    !is.null(a$W)
  }, logical(1L))
  W_global_df <- NULL
  if (any(W_active_per_slot) && length(design_K$X_names) > 0L) {
    first_active <- which(W_active_per_slot)[1L]
    W_spec <- amm_list[[first_active]]$W
    if (!is.null(W_spec) && identical(W_spec$type, "polynomial") &&
        W_spec$dim > 0L) {
      dim_W <- as.integer(W_spec$dim)
      d_x <- length(design_K$X_names)
      W_dr_arr <- posterior_var_to_matrix_3d(draws, "W_raw",
                                             rows = dim_W, cols = d_x)
      cp_W <- isTRUE(object$parametrization$cp_W)
      if (!cp_W) {
        sigma_W_dr <- posterior_var_to_array_2d(draws, "sigma_W", 1L)
        for (s in seq_len(dim(W_dr_arr)[1L])) {
          W_dr_arr[s, , ] <- W_dr_arr[s, , ] * sigma_W_dr[s, 1L]
        }
      }
      W_flat <- matrix(W_dr_arr, nrow = dim(W_dr_arr)[1L],
                       ncol = dim_W * d_x)
      W_global_df <- build_coef_W_df(W_flat, basis_dim = dim_W,
                                     x_names = design_K$X_names)
    }
  }
  out <- vector("list", K)
  for (k in seq_len(K)) {
    theta_ref_df_k <- build_coef_theta_ref_df(
      matrix(theta_ref_mat[, k], nrow = S, ncol = 1L), p = 1L
    )
    a_list_k <- NULL
    if (ncols_a[k] > 0L) {
      a_list_k <- list(build_coef_term_df(
        a_dr_list[[k]], term_names = design_K$Z_a_k_names_list[[k]]
      ))
    }
    b_list_k <- NULL
    if (ncols_b[k] > 0L) {
      b_list_k <- list(build_coef_term_df(
        b_dr_list[[k]], term_names = design_K$Z_b_k_names_list[[k]]
      ))
    }
    W_list_k <- NULL
    if (!is.null(W_global_df)) {
      W_list_k <- list(W_global_df)
    }
    out[[k]] <- new_gdpar_coef(
      theta_ref       = theta_ref_df_k,
      a               = a_list_k,
      b               = b_list_k,
      W               = W_list_k,
      p               = 1L,
      mu_theta_ref    = NULL,
      sigma_theta_ref = NULL,
      J_groups        = 1L,
      group_levels    = NULL
    )
  }
  names(out) <- slot_names
  out
}

#' Internal: build a gdpar_coef object from a scalar (p=1) fit
#'
#' @keywords internal
#' @noRd
coef_scalar_to_gdpar_coef <- function(object, draws) {
  J_groups <- if (!is.null(object$group_info))
    length(object$group_info$levels) else 1L
  group_levels <- if (!is.null(object$group_info))
    object$group_info$levels else NULL
  theta_ref_arr <- .extract_theta_ref_uni_grouped(draws, J_groups)
  theta_ref_df <- build_coef_theta_ref_df_grouped(theta_ref_arr,
                                                  J_groups, p = 1L)
  mu_df <- NULL
  sigma_df <- NULL
  if (!is.null(object$group_info)) {
    mu_dr <- .extract_mu_sigma_theta_ref(draws, "mu_theta_ref",
                                          p = 1L)
    sigma_dr <- .extract_mu_sigma_theta_ref(draws, "sigma_theta_ref",
                                             p = 1L)
    mu_df <- build_coef_hyper_df(mu_dr, p = 1L)
    sigma_df <- build_coef_hyper_df(sigma_dr, p = 1L)
  }
  a_list <- NULL
  if (!is.null(object$amm$a)) {
    a_dr <- unclass(posterior::as_draws_matrix(
      posterior::subset_draws(draws, variable = "a_coef")
    ))
    a_df <- build_coef_term_df(a_dr, object$design$Z_a_names)
    a_list <- list(a_df)
  }
  b_list <- NULL
  if (!is.null(object$amm$b)) {
    b_dr <- unclass(posterior::as_draws_matrix(
      posterior::subset_draws(draws, variable = "b_coef")
    ))
    b_df <- build_coef_term_df(b_dr, object$design$Z_b_names)
    b_list <- list(b_df)
  }
  W_list <- NULL
  if (!is.null(object$amm$W)) {
    cp_W <- isTRUE(object$parametrization$cp_W)
    W_dr <- unclass(posterior::as_draws_matrix(
      posterior::subset_draws(draws, variable = "W_raw")
    ))
    if (!cp_W) {
      sigma_W_dr <- unclass(posterior::as_draws_matrix(
        posterior::subset_draws(draws, variable = "sigma_W")
      ))
      W_eff <- W_dr * as.numeric(sigma_W_dr[, 1L])
    } else {
      W_eff <- W_dr
    }
    W_df <- build_coef_W_df(
      W_eff,
      basis_dim = object$amm$W$dim,
      x_names   = object$design$X_names
    )
    W_list <- list(W_df)
  }
  new_gdpar_coef(
    theta_ref       = theta_ref_df,
    a               = a_list,
    b               = b_list,
    W               = W_list,
    p               = 1L,
    mu_theta_ref    = mu_df,
    sigma_theta_ref = sigma_df,
    J_groups        = J_groups,
    group_levels    = group_levels
  )
}

#' Internal: build a gdpar_coef object from a multi (p>1) fit
#'
#' @keywords internal
#' @noRd
coef_multi_to_gdpar_coef <- function(object, draws) {
  amm <- object$amm
  design <- object$design
  p <- as.integer(object$p)
  J_groups <- if (!is.null(object$group_info))
    length(object$group_info$levels) else 1L
  group_levels <- if (!is.null(object$group_info))
    object$group_info$levels else NULL
  theta_ref_arr <- .extract_theta_ref_multi_grouped(draws, J_groups, p)
  theta_ref_df <- build_coef_theta_ref_df_grouped(theta_ref_arr,
                                                  J_groups, p)
  mu_df <- NULL
  sigma_df <- NULL
  if (!is.null(object$group_info)) {
    mu_dr <- .extract_mu_sigma_theta_ref(draws, "mu_theta_ref",
                                          p = p)
    sigma_dr <- .extract_mu_sigma_theta_ref(draws, "sigma_theta_ref",
                                             p = p)
    mu_df <- build_coef_hyper_df(mu_dr, p = p)
    sigma_df <- build_coef_hyper_df(sigma_dr, p = p)
  }
  a_list <- NULL
  has_a_global <- any(vapply(amm$dims, function(d) !is.null(d$a),
                             logical(1L)))
  if (has_a_global) {
    ncols_a <- vapply(seq_len(p), function(k) {
      length(design$Z_a_names_list[[k]])
    }, integer(1L))
    a_dr_list <- posterior_array_var_to_list(draws, "a_coef", p,
                                             ncols = ncols_a)
    a_list <- lapply(seq_len(p), function(k) {
      if (ncols_a[k] == 0L) return(NULL)
      build_coef_term_df(a_dr_list[[k]],
                         term_names = design$Z_a_names_list[[k]])
    })
  }
  b_list <- NULL
  has_b_global <- any(vapply(amm$dims, function(d) !is.null(d$b),
                             logical(1L)))
  if (has_b_global) {
    ncols_b <- vapply(seq_len(p), function(k) {
      length(design$Z_b_names_list[[k]])
    }, integer(1L))
    b_dr_list <- posterior_array_var_to_list(draws, "c_b", p,
                                             ncols = ncols_b)
    b_list <- lapply(seq_len(p), function(k) {
      if (ncols_b[k] == 0L) return(NULL)
      build_coef_term_df(b_dr_list[[k]],
                         term_names = design$Z_b_names_list[[k]])
    })
  }
  W_list <- NULL
  has_W <- !is.null(amm$W) && length(design$X_names) > 0L &&
    amm$W$type == "polynomial"
  if (has_W) {
    dim_W <- amm$W$dim
    W_per_k_dim <- as.integer(dim_W / p)
    d_x <- length(design$X_names)
    W_raw_dr <- posterior_var_to_matrix_3d(draws, "W_raw",
                                           rows = dim_W, cols = d_x)
    cp_W <- isTRUE(object$parametrization$cp_W)
    sigma_W_dr <- if (cp_W) NULL else
      posterior_var_to_array_2d(draws, "sigma_W", 1L)
    W_list <- lapply(seq_len(p), function(k) {
      if (W_per_k_dim == 0L) return(NULL)
      row_start <- (k - 1L) * W_per_k_dim + 1L
      row_end <- k * W_per_k_dim
      W_block_dr <- W_raw_dr[, row_start:row_end, , drop = FALSE]
      if (!cp_W) {
        for (s in seq_len(dim(W_block_dr)[1L])) {
          W_block_dr[s, , ] <- W_block_dr[s, , ] * sigma_W_dr[s, 1L]
        }
      }
      W_block_flat <- matrix(W_block_dr,
                             nrow = dim(W_block_dr)[1L],
                             ncol = W_per_k_dim * d_x)
      build_coef_W_df(W_block_flat,
                      basis_dim = W_per_k_dim,
                      x_names   = design$X_names)
    })
  }
  new_gdpar_coef(
    theta_ref       = theta_ref_df,
    a               = a_list,
    b               = b_list,
    W               = W_list,
    p               = p,
    mu_theta_ref    = mu_df,
    sigma_theta_ref = sigma_df,
    J_groups        = J_groups,
    group_levels    = group_levels
  )
}

#' Print method for gdpar_diagnostics objects
#'
#' @param x An object of class \code{gdpar_diagnostics}.
#' @param ... Unused; present for S3 generic compatibility.
#' @return Invisibly returns \code{x}.
#' @export
print.gdpar_diagnostics <- function(x, ...) {
  cat("<gdpar_diagnostics>\n")
  cat("  converged           : ", x$converged, "\n", sep = "")
  cat("  rhat_max            : ",
      if (is.na(x$rhat_max)) "NA" else format(x$rhat_max, digits = 4),
      "\n", sep = "")
  cat("  ess_bulk_min        : ",
      if (is.na(x$ess_bulk_min)) "NA" else format(x$ess_bulk_min, digits = 4),
      "\n", sep = "")
  cat("  ess_tail_min        : ",
      if (is.na(x$ess_tail_min)) "NA" else format(x$ess_tail_min, digits = 4),
      "\n", sep = "")
  cat("  divergent_count     : ", x$divergent_count, "\n", sep = "")
  cat("  treedepth_saturated : ", x$treedepth_saturated, "\n", sep = "")
  cat("  efmi_min            : ",
      if (is.na(x$efmi_min)) "NA" else format(x$efmi_min, digits = 4),
      "\n", sep = "")
  if (length(x$messages) > 0L) {
    cat("  messages:\n")
    for (m in x$messages) cat("    - ", m, "\n", sep = "")
  }
  invisible(x)
}

#' Internal: extract theta_ref draws (univariate Stan template)
#'
#' Parses the variable theta_ref of the univariate template, which is
#' declared as vector[J_groups] theta_ref and therefore emits posterior
#' draws under names theta_ref[g] for g in 1..J_groups. Returns a matrix
#' of dimensions (S, J_groups).
#'
#' @keywords internal
#' @noRd
.extract_theta_ref_uni_grouped <- function(draws, J_groups) {
  sub <- posterior::subset_draws(draws, variable = "theta_ref")
  mat <- unclass(posterior::as_draws_matrix(sub))
  var_names <- posterior::variables(sub)
  pat <- "^theta_ref\\[(\\d+)\\]$"
  matches <- regmatches(var_names, regexec(pat, var_names))
  idx <- as.integer(vapply(matches, function(m) m[2L], character(1L)))
  if (any(is.na(idx)) || length(idx) != J_groups) {
    gdpar_abort(
      sprintf(
        "Internal error: expected %d entries for 'theta_ref', found %d.",
        J_groups, length(idx)
      ),
      class = "gdpar_internal_error"
    )
  }
  S <- nrow(mat)
  out <- matrix(NA_real_, nrow = S, ncol = J_groups)
  for (col in seq_len(ncol(mat))) {
    out[, idx[col]] <- mat[, col]
  }
  out
}

#' Internal: extract theta_ref draws (multivariate Stan template)
#'
#' Parses the variable theta_ref of the multivariate template, which is
#' declared as array[J_groups] vector[p] theta_ref and therefore emits
#' posterior draws under names theta_ref[g, k]. Returns an array of
#' dimensions (S, J_groups, p).
#'
#' @keywords internal
#' @noRd
.extract_theta_ref_multi_grouped <- function(draws, J_groups, p) {
  sub <- posterior::subset_draws(draws, variable = "theta_ref")
  mat <- unclass(posterior::as_draws_matrix(sub))
  var_names <- posterior::variables(sub)
  pat <- "^theta_ref\\[(\\d+),\\s*(\\d+)\\]$"
  matches <- regmatches(var_names, regexec(pat, var_names))
  idx_g <- as.integer(vapply(matches, function(m) m[2L], character(1L)))
  idx_k <- as.integer(vapply(matches, function(m) m[3L], character(1L)))
  expected <- J_groups * p
  if (any(is.na(idx_g)) || any(is.na(idx_k)) ||
      length(idx_g) != expected) {
    gdpar_abort(
      sprintf(
        "Internal error: expected %d entries for 'theta_ref', found %d.",
        expected, length(idx_g)
      ),
      class = "gdpar_internal_error"
    )
  }
  S <- nrow(mat)
  out <- array(NA_real_, dim = c(S, J_groups, p))
  for (col in seq_len(ncol(mat))) {
    out[, idx_g[col], idx_k[col]] <- mat[, col]
  }
  out
}

#' Internal: flat (S, p) matrix view of theta_ref for the J_groups = 1 case
#'
#' Convenience wrapper used by callers that do not need the per-group
#' axis (the no-grouping branch under the multivariate Stan template):
#' collapses the (S, 1, p) array to an (S, p) matrix.
#'
#' @keywords internal
#' @noRd
.extract_theta_ref_multi_flat <- function(draws, J_groups, p) {
  arr <- .extract_theta_ref_multi_grouped(draws, J_groups, p)
  if (J_groups != 1L) {
    gdpar_abort(
      "Internal error: .extract_theta_ref_multi_flat called with J_groups != 1.",
      class = "gdpar_internal_error"
    )
  }
  matrix(arr[, 1L, ], nrow = dim(arr)[1L], ncol = p)
}

#' Internal: extract hierarchical hyperparameters mu_theta_ref and sigma_theta_ref
#'
#' Both are emitted by Stan only when use_groups = 1, declared as
#' array[1] real (univariate) or array[1] vector[p] (multivariate). The
#' draws variable names are therefore mu_theta_ref[1] (univariate) or
#' mu_theta_ref[1, k] (multivariate). The helper handles both shapes and
#' returns a matrix (S, p), collapsing the leading singleton axis.
#'
#' @keywords internal
#' @noRd
.extract_mu_sigma_theta_ref <- function(draws, var_name, p) {
  sub <- posterior::subset_draws(draws, variable = var_name)
  mat <- unclass(posterior::as_draws_matrix(sub))
  var_names <- posterior::variables(sub)
  pat_2idx <- sprintf("^%s\\[1,\\s*(\\d+)\\]$", var_name)
  matches_2idx <- regmatches(var_names, regexec(pat_2idx, var_names))
  idx_k_2idx <- as.integer(vapply(matches_2idx,
                                   function(m) if (length(m) >= 2L)
                                     m[2L] else NA_character_,
                                   character(1L)))
  if (!any(is.na(idx_k_2idx))) {
    if (length(idx_k_2idx) != p) {
      gdpar_abort(
        sprintf(
          "Internal error: expected %d entries for '%s', found %d.",
          p, var_name, length(idx_k_2idx)
        ),
        class = "gdpar_internal_error"
      )
    }
    S <- nrow(mat)
    out <- matrix(NA_real_, nrow = S, ncol = p)
    for (col in seq_len(ncol(mat))) {
      out[, idx_k_2idx[col]] <- mat[, col]
    }
    return(out)
  }
  pat_1idx <- sprintf("^%s\\[(\\d+)\\]$", var_name)
  matches_1idx <- regmatches(var_names, regexec(pat_1idx, var_names))
  idx_1idx <- as.integer(vapply(matches_1idx,
                                 function(m) if (length(m) >= 2L)
                                   m[2L] else NA_character_,
                                 character(1L)))
  if (any(is.na(idx_1idx)) || length(idx_1idx) != 1L) {
    gdpar_abort(
      sprintf(
        "Internal error: failed to parse '%s' variable names from draws.",
        var_name
      ),
      class = "gdpar_internal_error"
    )
  }
  S <- nrow(mat)
  out <- matrix(NA_real_, nrow = S, ncol = 1L)
  out[, 1L] <- mat[, 1L]
  out
}

#' Internal: predict theta_i on newdata under the univariate grouped fit
#'
#' Mirrors predict_from_newdata but accepts the per-group anchor
#' structure of Block 6.5. For each row i of newdata, the function:
#' (a) resolves the grouping variable against the levels seen at fit
#' time; (b) for known levels reuses theta_ref[s, g_i] from the
#' posterior draws; (c) for unseen levels draws from the marginal prior
#' predictive Normal(mu_theta_ref[s], sigma_theta_ref[s]) and emits a
#' gdpar_predict_unseen_group_warning naming the offending levels. The
#' b-block uses c_b directly (not b_coef), since the derived b_coef is
#' not reported under grouping (decision D4 of Block 6.5).
#'
#' @keywords internal
#' @noRd
predict_from_newdata_grouped <- function(object, newdata, draws) {
  require_suggested("posterior", "extract posterior draws")
  amm <- object$amm
  design_train <- object$design
  group_info <- object$group_info
  var_name <- group_info$var_name
  if (!var_name %in% colnames(newdata)) {
    gdpar_abort(
      sprintf(
        paste(
          "Grouping variable '%s' not found in newdata. predict() under",
          "a grouped fit requires the same grouping variable to be",
          "present in newdata."
        ),
        var_name
      ),
      class = "gdpar_input_error",
      data = list(missing_variable = var_name)
    )
  }
  raw_new <- newdata[[var_name]]
  fac_new <- factor(raw_new, levels = group_info$levels)
  group_id_new <- as.integer(fac_new)
  unseen <- is.na(group_id_new)
  if (any(unseen)) {
    unseen_levels <- unique(as.character(raw_new[unseen]))
    gdpar_warn(
      sprintf(
        paste(
          "newdata contains %d row(s) with %d unseen level(s) of '%s'",
          "(%s). For those rows theta_ref is drawn from the marginal",
          "prior predictive Normal(mu_theta_ref, sigma_theta_ref);",
          "interpret those predictions with caution."
        ),
        sum(unseen), length(unseen_levels), var_name,
        paste(unseen_levels, collapse = ", ")
      ),
      class = "gdpar_predict_unseen_group_warning",
      data = list(variable = var_name, unseen_levels = unseen_levels)
    )
  }
  J_groups <- length(group_info$levels)
  theta_ref_dr <- .extract_theta_ref_uni_grouped(draws, J_groups)
  S <- nrow(theta_ref_dr)
  mu_dr <- .extract_mu_sigma_theta_ref(draws, "mu_theta_ref", p = 1L)
  sigma_dr <- .extract_mu_sigma_theta_ref(draws, "sigma_theta_ref",
                                          p = 1L)
  new_design <- list()
  if (!is.null(amm$a)) {
    Za_full <- stats::model.matrix(stats::update(amm$a, ~ . + 0),
                                    newdata)
    Za_full <- as.matrix(Za_full)
    new_design$Z_a <- sweep(Za_full, 2L, design_train$Z_a_means, "-")
  } else {
    new_design$Z_a <- matrix(numeric(0), nrow = nrow(newdata), ncol = 0)
  }
  if (!is.null(amm$b)) {
    Zb_full <- stats::model.matrix(stats::update(amm$b, ~ . + 0),
                                    newdata)
    Zb_full <- as.matrix(Zb_full)
    new_design$Z_b <- sweep(Zb_full, 2L, design_train$Z_b_means, "-")
  } else {
    new_design$Z_b <- matrix(numeric(0), nrow = nrow(newdata), ncol = 0)
  }
  if (!is.null(amm$W) && length(design_train$X_names) > 0L) {
    X_full <- as.matrix(newdata[, design_train$X_names, drop = FALSE])
    storage.mode(X_full) <- "double"
    X_centered <- sweep(X_full, 2L, design_train$X_means, "-")
    new_design$X <- sweep(X_centered, 2L, design_train$X_sds, "/")
  } else {
    new_design$X <- matrix(numeric(0), nrow = nrow(newdata), ncol = 0)
  }
  has_a <- ncol(new_design$Z_a) > 0L
  has_b <- ncol(new_design$Z_b) > 0L
  has_W <- ncol(new_design$X) > 0L &&
    !is.null(amm$W) && amm$W$type == "polynomial"
  a_coef_dr <- if (has_a) {
    unclass(posterior::as_draws_matrix(
      posterior::subset_draws(draws, variable = "a_coef")
    ))
  } else {
    NULL
  }
  c_b_dr <- if (has_b) {
    unclass(posterior::as_draws_matrix(
      posterior::subset_draws(draws, variable = "c_b")
    ))
  } else {
    NULL
  }
  W_raw_dr <- if (has_W) {
    unclass(posterior::as_draws_matrix(
      posterior::subset_draws(draws, variable = "W_raw")
    ))
  } else {
    NULL
  }
  sigma_W_dr <- if (has_W) {
    unclass(posterior::as_draws_matrix(
      posterior::subset_draws(draws, variable = "sigma_W")
    ))
  } else {
    NULL
  }
  d <- ncol(new_design$X)
  dim_W <- if (has_W) amm$W$dim else 0L
  theta_anchor <- object$anchor
  n_new <- nrow(newdata)
  out <- matrix(NA_real_, nrow = S, ncol = n_new)
  for (s in seq_len(S)) {
    tr_per_row <- numeric(n_new)
    if (any(!unseen)) {
      tr_per_row[!unseen] <- theta_ref_dr[s, group_id_new[!unseen]]
    }
    if (any(unseen)) {
      tr_per_row[unseen] <- stats::rnorm(sum(unseen),
                                          mean = mu_dr[s, 1L],
                                          sd   = sigma_dr[s, 1L])
    }
    eta <- tr_per_row
    if (has_a) {
      eta <- eta + as.numeric(new_design$Z_a %*% a_coef_dr[s, ])
    }
    if (has_b) {
      eta <- eta + as.numeric(new_design$Z_b %*% c_b_dr[s, ])
    }
    if (has_W) {
      W_mat <- matrix(W_raw_dr[s, ], nrow = dim_W, ncol = d)
      for (i in seq_len(n_new)) {
        basis_diff <- vapply(seq_len(dim_W), function(k) {
          tr_per_row[i]^k - theta_anchor^k
        }, numeric(1L))
        W_diff_x <- as.numeric(crossprod(W_mat, basis_diff)) *
          sigma_W_dr[s, 1L]
        eta[i] <- eta[i] + sum(new_design$X[i, ] * W_diff_x)
      }
    }
    out[s, ] <- eta
  }
  out
}

#' Internal: predict theta_i on newdata under the multivariate grouped fit
#'
#' Mirrors predict_from_newdata_multi with the per-group anchor of Block
#' 6.5. Same dispatch as the univariate grouped helper: known levels
#' read theta_ref[s, g_i, k]; unseen levels draw from the marginal prior
#' predictive Normal(mu_theta_ref[s, k], sigma_theta_ref[s, k]) per
#' coord and emit a single warning naming the offending levels.
#'
#' @keywords internal
#' @noRd
predict_from_newdata_grouped_multi <- function(object, newdata, draws) {
  require_suggested("posterior", "extract posterior draws")
  amm <- object$amm
  design_train <- object$design
  p <- object$p
  group_info <- object$group_info
  var_name <- group_info$var_name
  if (!var_name %in% colnames(newdata)) {
    gdpar_abort(
      sprintf(
        paste(
          "Grouping variable '%s' not found in newdata. predict() under",
          "a grouped fit requires the same grouping variable to be",
          "present in newdata."
        ),
        var_name
      ),
      class = "gdpar_input_error",
      data = list(missing_variable = var_name)
    )
  }
  raw_new <- newdata[[var_name]]
  fac_new <- factor(raw_new, levels = group_info$levels)
  group_id_new <- as.integer(fac_new)
  unseen <- is.na(group_id_new)
  if (any(unseen)) {
    unseen_levels <- unique(as.character(raw_new[unseen]))
    gdpar_warn(
      sprintf(
        paste(
          "newdata contains %d row(s) with %d unseen level(s) of '%s'",
          "(%s). For those rows theta_ref is drawn coord-wise from the",
          "marginal prior predictive Normal(mu_theta_ref, sigma_theta_ref);",
          "interpret those predictions with caution."
        ),
        sum(unseen), length(unseen_levels), var_name,
        paste(unseen_levels, collapse = ", ")
      ),
      class = "gdpar_predict_unseen_group_warning",
      data = list(variable = var_name, unseen_levels = unseen_levels)
    )
  }
  J_groups <- length(group_info$levels)
  theta_ref_arr <- .extract_theta_ref_multi_grouped(draws, J_groups, p)
  S <- dim(theta_ref_arr)[1L]
  mu_dr <- .extract_mu_sigma_theta_ref(draws, "mu_theta_ref", p = p)
  sigma_dr <- .extract_mu_sigma_theta_ref(draws, "sigma_theta_ref",
                                          p = p)
  n_new <- nrow(newdata)
  Z_a_list_new <- vector("list", p)
  Z_b_list_new <- vector("list", p)
  for (k in seq_len(p)) {
    a_k <- amm$dims[[k]]$a
    b_k <- amm$dims[[k]]$b
    if (!is.null(a_k)) {
      Z_full <- stats::model.matrix(stats::update(a_k, ~ . + 0),
                                     newdata)
      Z_full <- as.matrix(Z_full)
      Z_a_list_new[[k]] <- sweep(
        Z_full, 2L, design_train$Z_a_means_list[[k]], "-"
      )
    } else {
      Z_a_list_new[[k]] <- matrix(numeric(0), nrow = n_new, ncol = 0L)
    }
    if (!is.null(b_k)) {
      Z_full <- stats::model.matrix(stats::update(b_k, ~ . + 0),
                                     newdata)
      Z_full <- as.matrix(Z_full)
      Z_b_list_new[[k]] <- sweep(
        Z_full, 2L, design_train$Z_b_means_list[[k]], "-"
      )
    } else {
      Z_b_list_new[[k]] <- matrix(numeric(0), nrow = n_new, ncol = 0L)
    }
  }
  has_W <- !is.null(amm$W) && length(design_train$X_names) > 0L &&
    amm$W$type == "polynomial"
  if (has_W) {
    X_full <- as.matrix(newdata[, design_train$X_names, drop = FALSE])
    storage.mode(X_full) <- "double"
    X_centered <- sweep(X_full, 2L, design_train$X_means, "-")
    X_new <- sweep(X_centered, 2L, design_train$X_sds, "/")
  } else {
    X_new <- matrix(numeric(0), nrow = n_new, ncol = 0L)
  }
  has_a_any <- any(vapply(Z_a_list_new, function(z) ncol(z) > 0L,
                          logical(1L)))
  has_b_any <- any(vapply(Z_b_list_new, function(z) ncol(z) > 0L,
                          logical(1L)))
  a_coef_dr_list <- if (has_a_any) {
    posterior_array_var_to_list(draws, "a_coef", p,
                                ncols = vapply(Z_a_list_new, ncol,
                                                integer(1L)))
  } else {
    NULL
  }
  c_b_dr_list <- if (has_b_any) {
    posterior_array_var_to_list(draws, "c_b", p,
                                ncols = vapply(Z_b_list_new, ncol,
                                                integer(1L)))
  } else {
    NULL
  }
  if (has_W) {
    dim_W <- amm$W$dim
    W_per_k_dim <- as.integer(dim_W / p)
    d <- ncol(X_new)
    W_raw_dr <- posterior_var_to_matrix_3d(draws, "W_raw",
                                            rows = dim_W, cols = d)
    cp_W <- isTRUE(object$parametrization$cp_W)
    sigma_W_dr <- if (cp_W) {
      NULL
    } else {
      posterior_var_to_array_2d(draws, "sigma_W", 1L)
    }
  } else {
    dim_W <- 0L
    W_per_k_dim <- 0L
    d <- 0L
    W_raw_dr <- NULL
    sigma_W_dr <- NULL
    cp_W <- FALSE
  }
  theta_anchor <- object$anchor
  out <- array(NA_real_, dim = c(S, n_new, p))
  for (s in seq_len(S)) {
    for (k in seq_len(p)) {
      tr_per_row <- numeric(n_new)
      if (any(!unseen)) {
        tr_per_row[!unseen] <-
          theta_ref_arr[s, group_id_new[!unseen], k]
      }
      if (any(unseen)) {
        tr_per_row[unseen] <- stats::rnorm(sum(unseen),
                                            mean = mu_dr[s, k],
                                            sd   = sigma_dr[s, k])
      }
      eta_k <- tr_per_row
      if (!is.null(a_coef_dr_list) && ncol(Z_a_list_new[[k]]) > 0L) {
        eta_k <- eta_k +
          as.numeric(Z_a_list_new[[k]] %*% a_coef_dr_list[[k]][s, ])
      }
      if (!is.null(c_b_dr_list) && ncol(Z_b_list_new[[k]]) > 0L) {
        eta_k <- eta_k +
          as.numeric(Z_b_list_new[[k]] %*% c_b_dr_list[[k]][s, ])
      }
      if (has_W && d > 0L && W_per_k_dim > 0L) {
        row_start <- (k - 1L) * W_per_k_dim + 1L
        row_end <- k * W_per_k_dim
        W_block <- W_raw_dr[s, row_start:row_end, , drop = FALSE]
        dim(W_block) <- c(W_per_k_dim, d)
        for (i in seq_len(n_new)) {
          basis_diff <- vapply(seq_len(W_per_k_dim), function(jj) {
            tr_per_row[i]^jj - theta_anchor[k]^jj
          }, numeric(1L))
          W_diff_x_k <- as.numeric(crossprod(W_block, basis_diff))
          if (!cp_W) {
            W_diff_x_k <- W_diff_x_k * sigma_W_dr[s, 1L]
          }
          eta_k[i] <- eta_k[i] + sum(X_new[i, ] * W_diff_x_k)
        }
      }
      out[s, , k] <- eta_k
    }
  }
  out
}

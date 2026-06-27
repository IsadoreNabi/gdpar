#' Print method for gdpar_meta_learner_comparison objects
#'
#' Concise summary of the comparison: bridge identifier, number of
#' observations and methods, per-method timing and CI availability,
#' and a head of the three concordance matrices.
#'
#' @param x A \code{gdpar_meta_learner_comparison} object.
#' @param ... Unused; present for S3 generic compatibility.
#' @return Invisibly returns \code{x}.
#' @export
print.gdpar_meta_learner_comparison <- function(x, ...) {
  cat("<gdpar_meta_learner_comparison>\n")
  cat("  n_obs                : ", x$n_obs, "\n", sep = "")
  cat("  n_methods (external) : ", x$n_methods, "\n", sep = "")
  cat("  level                : ", x$level, "\n", sep = "")
  cat("  methods              :\n")
  for (nm in names(x$external)) {
    e <- x$external[[nm]]
    cat(sprintf(
      "    - %-12s  native_ci = %s  time = %.3f s  notes = %d  predict = %s\n",
      nm, e$native_ci, e$time_sec, length(e$notes), e$has_predict_fun
    ))
  }
  cat("\n  concordance matrices (m-by-m, m = 1 + n_methods):\n")
  cat("    RMSE:\n");    print(round(x$comparison$rmse, 4L))
  cat("\n    Pearson:\n"); print(round(x$comparison$pearson, 4L))
  cat("\n    MAD:\n");     print(round(x$comparison$mad, 4L))
  invisible(x)
}

#' Summary method for gdpar_meta_learner_comparison objects
#'
#' Returns a structured summary object with the three concordance
#' matrices in long format, per-method ATE (mean of \code{cate_mean}),
#' per-method ATE CI bounds (when the adapter exposes native
#' per-observation CIs, the bounds are the mean of the per-observation
#' bounds; otherwise \code{NA}), and per-method timing.
#'
#' @param object A \code{gdpar_meta_learner_comparison} object.
#' @param ... Unused; present for S3 generic compatibility.
#' @return A list of class \code{summary.gdpar_meta_learner_comparison}.
#' @export
summary.gdpar_meta_learner_comparison <- function(object, ...) {
  assert_inherits(object, "gdpar_meta_learner_comparison", "object")
  method_names <- c("bridge", names(object$external))
  ate_vec <- stats::setNames(numeric(length(method_names)), method_names)
  ate_lower <- stats::setNames(rep(NA_real_, length(method_names)),
                                method_names)
  ate_upper <- stats::setNames(rep(NA_real_, length(method_names)),
                                method_names)
  ate_vec["bridge"] <- mean(object$bridge_cate$cate_mean)
  if (!is.null(object$bridge_cate$cate_ci)) {
    ate_lower["bridge"] <- mean(object$bridge_cate$cate_ci[, 1L])
    ate_upper["bridge"] <- mean(object$bridge_cate$cate_ci[, 2L])
  }
  for (nm in names(object$external)) {
    e <- object$external[[nm]]
    ate_vec[nm] <- mean(e$cate_mean)
    if (!is.null(e$cate_ci)) {
      ate_lower[nm] <- mean(e$cate_ci[, 1L])
      ate_upper[nm] <- mean(e$cate_ci[, 2L])
    }
  }
  ate_table <- data.frame(
    method = method_names,
    ate    = unname(ate_vec),
    ate_lower = unname(ate_lower),
    ate_upper = unname(ate_upper),
    stringsAsFactors = FALSE
  )
  long_metrics <- .comparison_long(object$comparison)
  timing <- data.frame(
    method = names(object$external),
    time_sec = vapply(object$external,
                       function(e) e$time_sec, numeric(1L)),
    native_ci = vapply(object$external,
                       function(e) e$native_ci, logical(1L)),
    stringsAsFactors = FALSE
  )
  obj <- list(
    ate_table   = ate_table,
    metrics     = long_metrics,
    timing      = timing,
    level       = object$level,
    n_obs       = object$n_obs,
    n_methods   = object$n_methods
  )
  class(obj) <- c("summary.gdpar_meta_learner_comparison", "list")
  obj
}

#' Print method for summary.gdpar_meta_learner_comparison objects
#'
#' @param x A \code{summary.gdpar_meta_learner_comparison} object.
#' @param ... Unused; present for S3 generic compatibility.
#' @return Invisibly returns \code{x}.
#' @export
print.summary.gdpar_meta_learner_comparison <- function(x, ...) {
  cat("<summary.gdpar_meta_learner_comparison>\n")
  cat("  level    : ", x$level, "\n", sep = "")
  cat("  n_obs    : ", x$n_obs, "\n", sep = "")
  cat("  n_methods: ", x$n_methods, "\n", sep = "")
  cat("\nATE per method (CI bounds derived from per-observation native CIs when available):\n")
  print(x$ate_table, row.names = FALSE)
  cat("\nTiming and CI availability:\n")
  print(x$timing, row.names = FALSE)
  cat("\nPairwise concordance metrics (long format):\n")
  print(utils::head(x$metrics, 20L), row.names = FALSE)
  if (nrow(x$metrics) > 20L) {
    cat(sprintf("  ... (%d more rows)\n", nrow(x$metrics) - 20L))
  }
  invisible(x)
}

#' Internal: pivot the three concordance matrices into a long data frame
#' @keywords internal
#' @noRd
.comparison_long <- function(comparison) {
  rmse <- comparison$rmse
  pearson <- comparison$pearson
  mad <- comparison$mad
  nms <- rownames(rmse)
  out_rows <- list()
  k <- 1L
  for (i in seq_along(nms)) {
    for (j in seq_along(nms)) {
      if (i == j) next
      out_rows[[k]] <- data.frame(
        method_i = nms[i],
        method_j = nms[j],
        rmse    = rmse[i, j],
        pearson = pearson[i, j],
        mad     = mad[i, j],
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }
  do.call(rbind, out_rows)
}

#' Predict method for gdpar_meta_learner_comparison objects
#'
#' Re-evaluate the CATE on a new grid for every method in the
#' comparison. Adapters that expose \code{predict_fun} reuse the
#' cached fitted state without a refit; adapters that do not are
#' invoked through \code{fit_predict_fun} (full refit) and a
#' \code{gdpar_diagnostic_warning} is emitted. The bridge component is
#' re-evaluated via \code{predict.gdpar_causal_bridge}.
#'
#' @param object A \code{gdpar_meta_learner_comparison} object.
#' @param newdata Data frame with the new evaluation grid. Required.
#' @param level Optional numeric scalar in \code{(0, 1)} overriding
#'   the credible level used when constructing the original
#'   comparison. Default \code{NULL} reuses \code{object$level}.
#' @param bridge Optional \code{gdpar_causal_bridge} object to use
#'   instead of the bridge cached inside \code{object}. Default
#'   \code{NULL} reuses \code{object$bridge}; supply this argument
#'   when the cached bridge was stripped (e.g. after a \code{saveRDS}
#'   round-trip that lost the two fits).
#' @param data Optional list with components \code{X}, \code{T},
#'   \code{Y} (and optionally \code{X_newdata}) for the case of a
#'   forced re-fit. Default \code{NULL} attempts to recover the
#'   training data of each arm from the bridge's stored fits (same
#'   convention as \code{\link{gdpar_compare_meta_learners}}).
#' @param ... Reserved for future arguments; currently unused.
#' @return A list of class
#'   \code{predict.gdpar_meta_learner_comparison} with components
#'   \code{bridge}, \code{external}, \code{comparison}, and the new
#'   \code{newdata}. The structure mirrors a comparison object but
#'   without the cached state.
#' @export
predict.gdpar_meta_learner_comparison <- function(object, newdata,
                                                  level = NULL,
                                                  bridge = NULL,
                                                  data = NULL, ...) {
  assert_inherits(object, "gdpar_meta_learner_comparison", "object")
  assert_data_frame(newdata, "newdata")
  if (is.null(level)) {
    level <- object$level
  } else {
    assert_numeric_scalar(level, "level",
                          lower = 1e-3, upper = 1 - 1e-3)
  }
  if (is.null(bridge)) {
    bridge_obj <- object$bridge
  } else {
    bridge_obj <- bridge
  }
  if (!inherits(bridge_obj, "gdpar_causal_bridge")) {
    gdpar_abort(
      paste(
        "predict.gdpar_meta_learner_comparison cannot re-evaluate the",
        "bridge: the comparison does not embed a valid",
        "gdpar_causal_bridge object. Pass one via the 'bridge' argument."
      ),
      class = "gdpar_input_error"
    )
  }
  outcome_name <- .bridge_outcome_name(bridge_obj$fits$treat,
                                       bridge_obj$fits$ctrl)
  has_real_fits <- !is.null(bridge_obj$fits$treat$fit) &&
                   !is.null(bridge_obj$fits$ctrl$fit)
  if (has_real_fits) {
    bridge_pred <- stats::predict(bridge_obj, newdata = newdata,
                                  level = level, summary = "mean_ci")
  } else {
    bridge_pred <- list(
      cate_mean = if (nrow(newdata) == bridge_obj$n_obs)
                    bridge_obj$cate_mean else rep(NA_real_, nrow(newdata)),
      cate_ci   = if (nrow(newdata) == bridge_obj$n_obs)
                    bridge_obj$cate_ci   else NULL
    )
  }
  X_newdata <- .extract_covariates(newdata, outcome_name)
  external <- vector("list", length(object$external))
  names(external) <- names(object$external)
  needs_refit <- character(0L)
  for (nm in names(object$external)) {
    e <- object$external[[nm]]
    pf <- e$predict_fun
    if (is.function(pf)) {
      out <- tryCatch(
        pf(state = e$state, X_newdata = X_newdata, level = level),
        error = function(err) err
      )
      if (inherits(out, "error")) {
        needs_refit <- c(needs_refit, nm)
        external[[nm]] <- list(
          cate_mean = rep(NA_real_, nrow(newdata)),
          cate_ci   = NULL,
          method    = nm,
          native_ci = FALSE,
          time_sec  = NA_real_,
          notes     = c(e$notes,
                         sprintf("predict_fun failed: %s",
                                 conditionMessage(out)))
        )
        next
      }
      external[[nm]] <- list(
        cate_mean = as.numeric(out$cate_mean),
        cate_ci   = out$cate_ci,
        method    = nm,
        native_ci = e$native_ci && !is.null(out$cate_ci),
        time_sec  = NA_real_,
        notes     = c(e$notes, "reused cached state via predict_fun")
      )
      next
    }
    needs_refit <- c(needs_refit, nm)
    external[[nm]] <- list(
      cate_mean = rep(NA_real_, nrow(newdata)),
      cate_ci   = NULL,
      method    = nm,
      native_ci = FALSE,
      time_sec  = NA_real_,
      notes     = c(e$notes,
                     "predict_fun unavailable; a full refit would be required")
    )
  }
  if (length(needs_refit) > 0L) {
    gdpar_warn(
      sprintf(
        paste(
          "predict.gdpar_meta_learner_comparison: adapters (%s) did",
          "not expose a usable predict_fun; the corresponding entries",
          "are returned with NA cate_mean. Rebuild the comparison",
          "with gdpar_compare_meta_learners() to obtain fresh fits."
        ),
        paste(needs_refit, collapse = ", ")
      ),
      class = "gdpar_diagnostic_warning",
      data = list(needs_refit = needs_refit)
    )
  }
  cate_list <- c(
    list(bridge = as.numeric(bridge_pred$cate_mean)),
    lapply(external, function(e) e$cate_mean)
  )
  comparison <- .compute_comparison_metrics(cate_list)
  obj <- list(
    bridge     = list(cate_mean = bridge_pred$cate_mean,
                      cate_ci   = bridge_pred$cate_ci),
    external   = external,
    comparison = comparison,
    newdata    = newdata,
    level      = level
  )
  class(obj) <- c("predict.gdpar_meta_learner_comparison", "list")
  obj
}


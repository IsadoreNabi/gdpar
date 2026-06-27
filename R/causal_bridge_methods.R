#' Print method for gdpar_causal_bridge objects
#'
#' Concise summary of the bridge: structural compatibility (family,
#' AMM level, anchor, K, p), number of posterior draws, number of
#' evaluation observations, and the credible level used for the
#' intervals.
#'
#' @param x An object of class \code{gdpar_causal_bridge}.
#' @param ... Unused; present for S3 generic compatibility.
#' @return Invisibly returns \code{x}.
#' @export
print.gdpar_causal_bridge <- function(x, ...) {
  cat("<gdpar_causal_bridge>\n")
  fam <- x$fits$treat$family
  cat("  family               : ", fam$name,
      " (link = ", fam$link, ")\n", sep = "")
  K_val <- x$fits$treat[["K"]] %||% 1L
  p_val <- x$fits$treat[["p"]] %||% 1L
  if (K_val > 1L) {
    cat("  K                    : ", K_val, "\n", sep = "")
  } else if (p_val > 1L) {
    cat("  p                    : ", p_val, "\n", sep = "")
  }
  anchor <- x$fits$treat$anchor
  cat("  anchor               : [",
      paste(format(anchor, digits = 4), collapse = ", "),
      "]\n", sep = "")
  cat("  type                 : ", x$type, "\n", sep = "")
  cat("  level                : ", x$level, "\n", sep = "")
  cat("  n_obs                : ", x$n_obs, "\n", sep = "")
  cat("  n_draws              : ", x$n_draws, "\n", sep = "")
  cat("  newdata source       : ",
      x$meta$newdata_source %||% "<unknown>",
      "\n", sep = "")
  cat("  cate_mean (head)     : ",
      .bridge_format_head(x$cate_mean),
      "\n", sep = "")
  if (length(x$warnings) > 0L) {
    cat("\nWarnings:\n")
    for (w in x$warnings) cat("  - ", w, "\n", sep = "")
  }
  invisible(x)
}

#' Internal: format the head of cate_mean for the print method
#' @keywords internal
#' @noRd
.bridge_format_head <- function(cate_mean, n_show = 6L) {
  if (is.matrix(cate_mean)) {
    n <- nrow(cate_mean)
    n_show <- min(n_show, n)
    rows <- vapply(seq_len(n_show), function(i) {
      sprintf("[%s]",
              paste(format(cate_mean[i, ], digits = 3), collapse = ", "))
    }, character(1L))
    paste(rows, collapse = "; ")
  } else {
    n <- length(cate_mean)
    n_show <- min(n_show, n)
    paste(format(cate_mean[seq_len(n_show)], digits = 3), collapse = ", ")
  }
}

#' Summary method for gdpar_causal_bridge objects
#'
#' Returns a structured summary object with a per-observation table of
#' the posterior CATE (mean, lower and upper credible bounds), the
#' marginal average treatment effect (ATE) computed as the mean of the
#' per-observation CATE, and the credible level used.
#'
#' For scalar bridges (\code{K = 1, p = 1}) the table has one row per
#' observation in \code{newdata}. For multivariate (\code{p > 1}) or
#' K-individual (\code{K > 1}) bridges, the table has one row per
#' (observation, dim/slot) pair and includes a \code{slot} column.
#'
#' @param object An object of class \code{gdpar_causal_bridge}.
#' @param ... Unused; present for S3 generic compatibility.
#' @return A list of class \code{summary.gdpar_causal_bridge} with
#'   components \code{table} (data frame), \code{ate} (named vector of
#'   marginal ATE per slot), \code{ate_ci} (matrix of marginal ATE
#'   credible bounds per slot), \code{level}, \code{type},
#'   \code{n_draws}, \code{n_obs}. The companion \code{print} method
#'   formats the object.
#' @export
summary.gdpar_causal_bridge <- function(object, ...) {
  assert_inherits(object, "gdpar_causal_bridge", "object")
  level <- object$level
  alpha <- 1 - level
  ql <- alpha / 2
  qu <- 1 - alpha / 2
  cate_draws <- object$cate_draws
  if (is.matrix(cate_draws)) {
    n <- ncol(cate_draws)
    table <- data.frame(
      observation = seq_len(n),
      cate_mean   = object$cate_mean,
      cate_lower  = object$cate_ci[, 1L],
      cate_upper  = object$cate_ci[, 2L]
    )
    ate_draws <- rowMeans(cate_draws)
    ate <- mean(ate_draws)
    ate_ci <- stats::quantile(ate_draws, probs = c(ql, qu),
                              names = FALSE)
    ate_vec <- c(ate = ate)
    ate_mat <- matrix(ate_ci, nrow = 1L, ncol = 2L,
                       dimnames = list("ate", c("lower", "upper")))
  } else {
    n <- dim(cate_draws)[2L]
    K <- dim(cate_draws)[3L]
    slot_names <- object$meta$dim_names
    if (is.null(slot_names)) {
      slot_names <- paste0("dim_", seq_len(K))
    }
    table_rows <- list()
    for (k in seq_len(K)) {
      table_rows[[k]] <- data.frame(
        observation = seq_len(n),
        slot        = slot_names[k],
        cate_mean   = object$cate_mean[, k],
        cate_lower  = object$cate_ci[, k, 1L],
        cate_upper  = object$cate_ci[, k, 2L]
      )
    }
    table <- do.call(rbind, table_rows)
    ate_vec <- vapply(seq_len(K), function(k) mean(cate_draws[, , k]),
                       numeric(1L))
    names(ate_vec) <- slot_names
    ate_mat <- matrix(NA_real_, nrow = K, ncol = 2L,
                       dimnames = list(slot_names, c("lower", "upper")))
    for (k in seq_len(K)) {
      mean_per_draw <- apply(cate_draws[, , k, drop = FALSE], 1L, mean)
      qs <- stats::quantile(mean_per_draw, probs = c(ql, qu),
                            names = FALSE)
      ate_mat[k, ] <- qs
    }
  }
  obj <- list(
    table   = table,
    ate     = ate_vec,
    ate_ci  = ate_mat,
    level   = level,
    type    = object$type,
    n_draws = object$n_draws,
    n_obs   = object$n_obs
  )
  class(obj) <- c("summary.gdpar_causal_bridge", "list")
  obj
}

#' Print method for summary.gdpar_causal_bridge objects
#'
#' @param x An object of class \code{summary.gdpar_causal_bridge}.
#' @param ... Unused; present for S3 generic compatibility.
#' @return Invisibly returns \code{x}.
#' @export
print.summary.gdpar_causal_bridge <- function(x, ...) {
  cat("<summary.gdpar_causal_bridge>\n")
  cat("  type                 : ", x$type, "\n", sep = "")
  cat("  level                : ", x$level, "\n", sep = "")
  cat("  n_obs                : ", x$n_obs, "\n", sep = "")
  cat("  n_draws              : ", x$n_draws, "\n", sep = "")
  cat("\nMarginal ATE:\n")
  ate_df <- data.frame(
    slot  = names(x$ate),
    mean  = unname(x$ate),
    lower = x$ate_ci[, 1L],
    upper = x$ate_ci[, 2L],
    stringsAsFactors = FALSE
  )
  print(ate_df, row.names = FALSE)
  cat("\nPer-observation CATE (head):\n")
  print(utils::head(x$table, 10L), row.names = FALSE)
  invisible(x)
}

#' Predict method for gdpar_causal_bridge objects
#'
#' Recompute the per-observation CATE on a new evaluation grid using
#' the two underlying fits stored in the bridge. The structural
#' compatibility of the two fits was validated when the bridge was
#' constructed and is not re-checked.
#'
#' @param object An object of class \code{gdpar_causal_bridge}.
#' @param newdata Data frame on which to evaluate the CATE. Required.
#' @param level Numeric scalar in (0, 1) with the credible level for
#'   the new intervals. Defaults to the level recorded on
#'   \code{object}.
#' @param summary Character scalar selecting the output form:
#'   \code{"all"} (default; returns the same structure as
#'   \code{gdpar_causal_bridge}'s \code{cate_*} slots),
#'   \code{"draws"} (returns the raw \code{cate_draws} object),
#'   \code{"mean_ci"} (returns a list with \code{cate_mean} and
#'   \code{cate_ci}).
#' @param ... Unused; present for S3 generic compatibility.
#' @return Depends on \code{summary}: a list with components
#'   \code{cate_draws}, \code{cate_mean}, \code{cate_ci},
#'   \code{n_draws}, \code{n_obs} for \code{summary = "all"}; the
#'   draws array for \code{summary = "draws"}; or a list with the two
#'   summary slots for \code{summary = "mean_ci"}.
#' @export
predict.gdpar_causal_bridge <- function(object, newdata,
                                        level = NULL,
                                        summary = c("all", "draws",
                                                    "mean_ci"),
                                        ...) {
  assert_inherits(object, "gdpar_causal_bridge", "object")
  assert_data_frame(newdata, "newdata")
  summary <- match.arg(summary)
  if (is.null(level)) {
    level <- object$level
  } else {
    assert_numeric_scalar(level, "level", lower = 1e-3, upper = 1 - 1e-3)
  }
  type <- object$type
  fit_t <- object$fits$treat
  fit_c <- object$fits$ctrl
  pred_t <- stats::predict(fit_t, newdata = newdata, type = type,
                           summary = "draws")
  pred_c <- stats::predict(fit_c, newdata = newdata, type = type,
                           summary = "draws")
  aligned <- .align_bridge_draws(pred_t, pred_c)
  cate_draws <- aligned$treat - aligned$ctrl
  if (summary == "draws") {
    return(cate_draws)
  }
  alpha <- 1 - level
  ql <- alpha / 2
  qu <- 1 - alpha / 2
  summaries <- .summarize_cate(cate_draws, ql = ql, qu = qu)
  if (summary == "mean_ci") {
    return(list(cate_mean = summaries$mean, cate_ci = summaries$ci))
  }
  list(
    cate_draws = cate_draws,
    cate_mean  = summaries$mean,
    cate_ci    = summaries$ci,
    n_draws    = aligned$S,
    n_obs      = nrow(newdata)
  )
}

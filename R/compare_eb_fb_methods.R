#' Print method for gdpar_eb_fb_comparison
#'
#' Concise console summary of an Empirical-Bayes vs Fully-Bayes
#' comparison: paths involved, number of common xi parameters with TV
#' values, summary statistics of the TV distribution and width-ratio
#' distribution, and the first six rows of the per-anchor diff table.
#'
#' @param x A \code{gdpar_eb_fb_comparison} object.
#' @param digits Integer scalar passed to \code{format()}; defaults to
#'   3.
#' @param ... Unused.
#' @return The object \code{x} invisibly.
#' @export
print.gdpar_eb_fb_comparison <- function(x, digits = 3L, ...) {
  cat("<gdpar_eb_fb_comparison>\n")
  cat(sprintf("  family_eb / family_fb : %s / %s\n",
              x$family_eb, x$family_fb))
  cat(sprintf("  path_eb / path_fb     : %s / %s\n",
              x$path_eb, x$path_fb))
  cat(sprintf("  level                 : %s\n",
              format(x$level, digits = digits)))
  cat(sprintf("  tv_bins               : %d\n", x$tv_bins))
  cat(sprintf("  n common xi params    : %d\n",
              x$n_common_params))
  if (!is.null(x$tv_table) && nrow(x$tv_table) > 0L) {
    tv_vals <- x$tv_table$tv[is.finite(x$tv_table$tv)]
    if (length(tv_vals) > 0L) {
      cat(sprintf("  marginal TV (min / median / max): %s / %s / %s\n",
                  format(min(tv_vals), digits = digits),
                  format(stats::median(tv_vals), digits = digits),
                  format(max(tv_vals), digits = digits)))
    }
  }
  if (!is.null(x$coverage_table) && nrow(x$coverage_table) > 0L) {
    wr <- x$coverage_table$width_ratio
    wr <- wr[is.finite(wr)]
    if (length(wr) > 0L) {
      cat(sprintf("  width_ratio EB/FB (min / median / max): %s / %s / %s\n",
                  format(min(wr), digits = digits),
                  format(stats::median(wr), digits = digits),
                  format(max(wr), digits = digits)))
    }
  }
  if (!is.null(x$theta_diff_table) && nrow(x$theta_diff_table) > 0L) {
    cat("\n  theta_diff (first 6 rows):\n")
    print(format(utils::head(x$theta_diff_table, 6L), digits = digits))
  }
  if (length(x$warnings) > 0L) {
    cat("\nWarnings:\n")
    for (w in x$warnings) cat("  - ", w, "\n", sep = "")
  }
  invisible(x)
}

#' Summary method for gdpar_eb_fb_comparison
#'
#' Returns a structured summary suitable for programmatic access and
#' for the canonical \code{print.summary.gdpar_eb_fb_comparison} method.
#' Aggregates the TV table (mean / median / max / quartiles) and the
#' coverage table (mean width_ratio per slot under Path C, or overall
#' under the other regimes).
#'
#' @param object A \code{gdpar_eb_fb_comparison} object.
#' @param ... Unused.
#' @return An object of class \code{summary.gdpar_eb_fb_comparison}.
#' @export
summary.gdpar_eb_fb_comparison <- function(object, ...) {
  tv_tbl <- object$tv_table
  tv_summary <- if (!is.null(tv_tbl) && nrow(tv_tbl) > 0L) {
    tv_vals <- tv_tbl$tv[is.finite(tv_tbl$tv)]
    if (length(tv_vals) > 0L) {
      list(
        n     = length(tv_vals),
        min   = min(tv_vals),
        q25   = unname(stats::quantile(tv_vals, 0.25)),
        median = stats::median(tv_vals),
        q75   = unname(stats::quantile(tv_vals, 0.75)),
        max   = max(tv_vals),
        mean  = mean(tv_vals)
      )
    } else NULL
  } else NULL
  cov_tbl <- object$coverage_table
  cov_summary <- if (!is.null(cov_tbl) && nrow(cov_tbl) > 0L) {
    wr <- cov_tbl$width_ratio[is.finite(cov_tbl$width_ratio)]
    if (length(wr) > 0L) {
      list(
        n     = length(wr),
        min   = min(wr),
        q25   = unname(stats::quantile(wr, 0.25)),
        median = stats::median(wr),
        q75   = unname(stats::quantile(wr, 0.75)),
        max   = max(wr),
        mean  = mean(wr)
      )
    } else NULL
  } else NULL
  out <- list(
    family_eb        = object$family_eb,
    family_fb        = object$family_fb,
    path_eb          = object$path_eb,
    path_fb          = object$path_fb,
    level            = object$level,
    tv_bins          = object$tv_bins,
    n_common_params  = object$n_common_params,
    n_anchor_cells   = if (is.null(object$theta_diff_table)) 0L
                       else nrow(object$theta_diff_table),
    tv_summary       = tv_summary,
    coverage_summary = cov_summary,
    theta_diff_table = object$theta_diff_table,
    tv_table         = tv_tbl,
    coverage_table   = cov_tbl,
    warnings         = object$warnings %||% character(0L),
    call             = object$call
  )
  class(out) <- c("summary.gdpar_eb_fb_comparison", "list")
  out
}

#' Print method for summary.gdpar_eb_fb_comparison
#'
#' @param x A \code{summary.gdpar_eb_fb_comparison} object.
#' @param digits Integer scalar passed to \code{format()}; defaults to
#'   3.
#' @param ... Unused.
#' @return The object \code{x} invisibly.
#' @export
print.summary.gdpar_eb_fb_comparison <- function(x, digits = 3L, ...) {
  cat("<gdpar_eb_fb_comparison summary>\n")
  cat(sprintf("  family_eb / family_fb : %s / %s\n",
              x$family_eb, x$family_fb))
  cat(sprintf("  path_eb / path_fb     : %s / %s\n",
              x$path_eb, x$path_fb))
  cat(sprintf("  level                 : %s\n",
              format(x$level, digits = digits)))
  cat(sprintf("  n common xi params    : %d\n", x$n_common_params))
  cat(sprintf("  n anchor cells        : %d\n", x$n_anchor_cells))
  if (!is.null(x$tv_summary)) {
    cat(sprintf("\n  Marginal TV summary (n = %d):\n", x$tv_summary$n))
    cat(sprintf("    min / q25 / median / q75 / max : %s / %s / %s / %s / %s\n",
                format(x$tv_summary$min, digits = digits),
                format(x$tv_summary$q25, digits = digits),
                format(x$tv_summary$median, digits = digits),
                format(x$tv_summary$q75, digits = digits),
                format(x$tv_summary$max, digits = digits)))
    cat(sprintf("    mean : %s\n",
                format(x$tv_summary$mean, digits = digits)))
  }
  if (!is.null(x$coverage_summary)) {
    cat(sprintf("\n  Width ratio EB/FB summary (n = %d):\n",
                x$coverage_summary$n))
    cat(sprintf("    min / q25 / median / q75 / max : %s / %s / %s / %s / %s\n",
                format(x$coverage_summary$min, digits = digits),
                format(x$coverage_summary$q25, digits = digits),
                format(x$coverage_summary$median, digits = digits),
                format(x$coverage_summary$q75, digits = digits),
                format(x$coverage_summary$max, digits = digits)))
    cat(sprintf("    mean : %s\n",
                format(x$coverage_summary$mean, digits = digits)))
  }
  if (!is.null(x$theta_diff_table) && nrow(x$theta_diff_table) > 0L) {
    cat("\n  theta_diff (full table):\n")
    print(format(x$theta_diff_table, digits = digits))
  }
  if (length(x$warnings) > 0L) {
    cat("\nWarnings:\n")
    for (w in x$warnings) cat("  - ", w, "\n", sep = "")
  }
  invisible(x)
}

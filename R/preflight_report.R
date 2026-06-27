#' Preflight report for multivariate (`p > 1`) CP/NCP decisions
#'
#' S3 class that bundles the per-coordinate, per-component CP/NCP
#' decisions produced by the multivariate preflight (Path B'
#' applied per coordinate) together with an aggregated per-component
#' summary controlled by an explicit aggregation strategy. The
#' constructor is internal; the report is created by
#' \code{preflight_parametrization()} when invoked on a spec with
#' \code{p > 1}.
#'
#' Aggregation strategies:
#' \itemize{
#'   \item \code{"any_ncp"} (default, conservative): the component is
#'         CP at the global level only if every coordinate's
#'         per-dim decision is CP. Any single NCP decision flips the
#'         component-wide decision to NCP. Reflects the asymmetry
#'         that NCP geometry is monotonically safer in the worst
#'         case (more warmup steps in high-info regime, fewer
#'         divergences in low-info regime).
#'   \item \code{"majority"}: component is CP if the strict majority
#'         of coordinates with an effective decision votes CP. Ties
#'         break toward NCP (conservative tie-breaking).
#'   \item \code{"per_k"}: no aggregation; the global slot records a
#'         sentinel \code{"per_k"} value indicating that the Stan
#'         template should consume the per-coordinate decision
#'         vector as-is (handled by the codegen in Phase H.2).
#' }
#'
#' The per-component \code{agreement} column reports the share of
#' coordinates whose per-dim decision matches the aggregated
#' \code{global_decision}. For the \code{"per_k"} strategy, it
#' reports the frequency of the modal per-dim decision (a summary
#' of whether the per-dim decisions are uniform or mixed).
#'
#' @section Inspection:
#' Use \code{preflight_per_dim()} and \code{preflight_global_decision()}
#' rather than touching attributes directly. The S3 methods
#' \code{print}, \code{summary}, and \code{as.data.frame} provide the
#' usual human/programmatic views.
#'
#' @name gdpar_preflight_report
#' @keywords internal
NULL


# ----- internal columns required by the schema --------------------

.preflight_per_dim_cols <- c(
  "component", "dim", "decision", "decision_reason",
  "n_divergent", "div_pct", "ebfmi_min",
  "t_attribution", "t_info_cp", "t_info_ncp"
)

.preflight_global_cols <- c(
  "component", "global_decision", "agreement", "method"
)

.preflight_decision_levels <- c("CP", "NCP", "absent")
.preflight_global_decision_levels <- c("CP", "NCP", "per_k", "absent")
.preflight_aggregation_methods <- c("any_ncp", "majority", "per_k")


#' Internal: validate the per-dim data frame schema
#'
#' Enforces presence and types of the required columns plus the
#' invariants on \code{decision} (one of the allowed levels) and on
#' the absence of negative \code{dim} or \code{n_divergent} values.
#' Side effect only; signals \code{gdpar_input_error} on failure and
#' returns invisibly on success.
#'
#' @keywords internal
#' @noRd
validate_preflight_per_dim <- function(df) {
  if (!is.data.frame(df)) {
    gdpar_abort(
      "Argument 'per_dim' must be a data frame.",
      class = "gdpar_input_error",
      data = list(received = class(df))
    )
  }
  missing_cols <- setdiff(.preflight_per_dim_cols, colnames(df))
  if (length(missing_cols) > 0L) {
    gdpar_abort(
      sprintf(
        "Data frame 'per_dim' is missing the following column(s): %s.",
        paste(sQuote(missing_cols), collapse = ", ")
      ),
      class = "gdpar_input_error",
      data = list(missing = missing_cols)
    )
  }
  if (!is.character(df$component) ||
      !all(df$component %in% c("a", "W"))) {
    gdpar_abort(
      "Column 'component' must be character with values in {'a', 'W'}.",
      class = "gdpar_input_error"
    )
  }
  if (!is.integer(df$dim) || any(df$dim < 1L) || any(is.na(df$dim))) {
    gdpar_abort(
      "Column 'dim' must be a positive integer with no NA.",
      class = "gdpar_input_error"
    )
  }
  if (!is.character(df$decision) ||
      !all(df$decision %in% .preflight_decision_levels)) {
    gdpar_abort(
      sprintf(
        "Column 'decision' must be character with values in {%s}.",
        paste(sQuote(.preflight_decision_levels), collapse = ", ")
      ),
      class = "gdpar_input_error"
    )
  }
  if (!is.character(df$decision_reason)) {
    gdpar_abort(
      "Column 'decision_reason' must be character.",
      class = "gdpar_input_error"
    )
  }
  if (!is.integer(df$n_divergent)) {
    gdpar_abort(
      "Column 'n_divergent' must be integer.",
      class = "gdpar_input_error"
    )
  }
  for (nm in c("div_pct", "ebfmi_min", "t_attribution",
               "t_info_cp", "t_info_ncp")) {
    if (!is.numeric(df[[nm]])) {
      gdpar_abort(
        sprintf("Column '%s' must be numeric.", nm),
        class = "gdpar_input_error"
      )
    }
  }
  dup_idx <- duplicated(df[, c("component", "dim")])
  if (any(dup_idx)) {
    gdpar_abort(
      "Duplicate (component, dim) pairs in 'per_dim'.",
      class = "gdpar_input_error"
    )
  }
  invisible(TRUE)
}


#' Internal: validate the global data frame schema
#'
#' @keywords internal
#' @noRd
validate_preflight_global <- function(df) {
  if (!is.data.frame(df)) {
    gdpar_abort(
      "Argument 'global' must be a data frame.",
      class = "gdpar_input_error",
      data = list(received = class(df))
    )
  }
  missing_cols <- setdiff(.preflight_global_cols, colnames(df))
  if (length(missing_cols) > 0L) {
    gdpar_abort(
      sprintf(
        "Data frame 'global' is missing the following column(s): %s.",
        paste(sQuote(missing_cols), collapse = ", ")
      ),
      class = "gdpar_input_error",
      data = list(missing = missing_cols)
    )
  }
  if (!is.character(df$component) ||
      !all(df$component %in% c("a", "W"))) {
    gdpar_abort(
      "Column 'component' must be character with values in {'a', 'W'}.",
      class = "gdpar_input_error"
    )
  }
  if (!is.character(df$global_decision) ||
      !all(df$global_decision %in% .preflight_global_decision_levels)) {
    gdpar_abort(
      sprintf(
        "Column 'global_decision' must be character with values in {%s}.",
        paste(sQuote(.preflight_global_decision_levels),
              collapse = ", ")
      ),
      class = "gdpar_input_error"
    )
  }
  if (!is.numeric(df$agreement)) {
    gdpar_abort(
      "Column 'agreement' must be numeric.",
      class = "gdpar_input_error"
    )
  }
  finite_agreement <- df$agreement[is.finite(df$agreement)]
  if (length(finite_agreement) > 0L &&
      (any(finite_agreement < 0) || any(finite_agreement > 1))) {
    gdpar_abort(
      "Column 'agreement' must lie in [0, 1] or be NA.",
      class = "gdpar_input_error"
    )
  }
  if (!is.character(df$method) ||
      !all(df$method %in% .preflight_aggregation_methods)) {
    gdpar_abort(
      sprintf(
        "Column 'method' must be character with values in {%s}.",
        paste(sQuote(.preflight_aggregation_methods),
              collapse = ", ")
      ),
      class = "gdpar_input_error"
    )
  }
  if (any(duplicated(df$component))) {
    gdpar_abort(
      "Duplicate 'component' entries in 'global'.",
      class = "gdpar_input_error"
    )
  }
  invisible(TRUE)
}


#' Internal: aggregate per-dim decisions for one component
#'
#' Returns a length-1 character with the aggregated decision.
#' \code{"any_ncp"}: CP only if all effective decisions are CP.
#' \code{"majority"}: CP if strict majority votes CP; ties to NCP.
#' \code{"per_k"}: returns the sentinel "per_k" when there is at
#' least one effective decision; "absent" otherwise.
#'
#' @keywords internal
#' @noRd
aggregate_preflight_one_component <- function(decisions, method) {
  effective <- decisions[decisions != "absent"]
  if (length(effective) == 0L) return("absent")
  switch(method,
    any_ncp  = if (all(effective == "CP")) "CP" else "NCP",
    majority = {
      n_cp <- sum(effective == "CP")
      n_ncp <- sum(effective == "NCP")
      if (n_cp > n_ncp) "CP" else "NCP"
    },
    per_k = "per_k",
    gdpar_abort(
      sprintf("Unknown aggregation method '%s'.", method),
      class = "gdpar_input_error",
      data = list(method = method)
    )
  )
}


#' Internal: compute the agreement score for one component
#'
#' Share of effective per-dim decisions that match the
#' \code{global_decision}. For \code{global_decision == "per_k"},
#' returns the frequency of the modal per-dim decision instead
#' (since there is no single global value to compare to). Returns
#' \code{NA_real_} when there is no effective decision.
#'
#' @keywords internal
#' @noRd
compute_preflight_agreement <- function(decisions, global_decision) {
  effective <- decisions[decisions != "absent"]
  if (length(effective) == 0L) return(NA_real_)
  if (identical(global_decision, "per_k")) {
    n_cp <- sum(effective == "CP")
    n_ncp <- sum(effective == "NCP")
    return(max(n_cp, n_ncp) / length(effective))
  }
  if (identical(global_decision, "absent")) return(NA_real_)
  mean(effective == global_decision)
}


#' Internal: build the global data frame from per-dim and method
#'
#' Iterates over components present in \code{per_dim}, applies the
#' aggregator, and computes agreement. Returns a data frame with
#' one row per component, ordered as in \code{per_dim} (which is
#' canonically \code{c("a", "W")} when both are present).
#'
#' @keywords internal
#' @noRd
build_preflight_global <- function(per_dim, method) {
  if (!method %in% .preflight_aggregation_methods) {
    gdpar_abort(
      sprintf("Unknown aggregation method '%s'.", method),
      class = "gdpar_input_error",
      data = list(method = method)
    )
  }
  comps <- unique(per_dim$component)
  global <- data.frame(
    component       = comps,
    global_decision = character(length(comps)),
    agreement       = numeric(length(comps)),
    method          = rep(method, length(comps)),
    stringsAsFactors = FALSE
  )
  for (i in seq_along(comps)) {
    rows <- per_dim$decision[per_dim$component == comps[i]]
    g <- aggregate_preflight_one_component(rows, method)
    global$global_decision[i] <- g
    global$agreement[i] <- compute_preflight_agreement(rows, g)
  }
  global
}


#' Internal: constructor for the S3 class
#'
#' @param per_dim Data frame with the per-coordinate, per-component
#'   decisions (schema enforced by \code{validate_preflight_per_dim}).
#' @param method Character scalar in \code{c("any_ncp", "majority",
#'   "per_k")}. Default \code{"any_ncp"}.
#' @param settings Named list with the preflight settings actually
#'   used (\code{tau_cp}, \code{tau_ncp}, \code{n_chains},
#'   \code{n_warmup}, \code{n_sampling}, \code{adapt_delta},
#'   \code{max_treedepth}). Stored for introspection.
#'
#' @return An object of class \code{gdpar_preflight_report}, which
#'   is a list with elements \code{per_dim} (data frame), \code{global}
#'   (data frame), \code{method} (character scalar), \code{settings}
#'   (list).
#' @keywords internal
#' @noRd
new_gdpar_preflight_report <- function(per_dim,
                                       method = "any_ncp",
                                       settings = list()) {
  if (!is.character(method) || length(method) != 1L ||
      !method %in% .preflight_aggregation_methods) {
    gdpar_abort(
      sprintf(
        "Argument 'method' must be one of {%s}.",
        paste(sQuote(.preflight_aggregation_methods),
              collapse = ", ")
      ),
      class = "gdpar_input_error",
      data = list(received = method)
    )
  }
  if (!is.list(settings)) {
    gdpar_abort(
      "Argument 'settings' must be a list.",
      class = "gdpar_input_error"
    )
  }
  validate_preflight_per_dim(per_dim)
  global <- build_preflight_global(per_dim, method)
  validate_preflight_global(global)
  structure(
    list(
      per_dim  = per_dim,
      global   = global,
      method   = method,
      settings = settings
    ),
    class = "gdpar_preflight_report"
  )
}


#' Accessor: per-dimension preflight decisions
#'
#' Returns the per-coordinate, per-component data frame that drives
#' the report. Use this for \code{dplyr}/\code{ggplot2} analysis.
#'
#' @param report An object of class \code{gdpar_preflight_report}.
#' @return Data frame with one row per (component, dim) pair.
#' @export
preflight_per_dim <- function(report) {
  assert_inherits(report, "gdpar_preflight_report", "report")
  report$per_dim
}


#' Accessor: aggregated per-component preflight decisions
#'
#' Returns the per-component summary (one row per component) with
#' the aggregated decision, the agreement share, and the
#' aggregation method used.
#'
#' @param report An object of class \code{gdpar_preflight_report}.
#' @return Data frame with one row per component.
#' @export
preflight_global_decision <- function(report) {
  assert_inherits(report, "gdpar_preflight_report", "report")
  report$global
}


#' Print method for gdpar_preflight_report
#'
#' Prints a compact summary by default (\code{level = "global"}). Use
#' \code{level = "dim"} to print only the per-coordinate table, or
#' \code{level = "both"} to print both.
#'
#' @param x An object of class \code{gdpar_preflight_report}.
#' @param level Character scalar: \code{"global"} (default),
#'   \code{"dim"}, or \code{"both"}.
#' @param ... Unused; present for S3 generic compatibility.
#' @return Invisibly returns \code{x}.
#' @export
print.gdpar_preflight_report <- function(x,
                                         level = c("global", "dim", "both"),
                                         ...) {
  level <- match.arg(level)
  cat("<gdpar_preflight_report>\n")
  n_dims <- length(unique(x$per_dim$dim))
  cat("  p (theta_ref dim) : ", n_dims, "\n", sep = "")
  cat("  aggregation       : ", x$method, "\n", sep = "")
  cat("  components        : ",
      paste(unique(x$per_dim$component), collapse = ", "),
      "\n", sep = "")
  if (level %in% c("global", "both")) {
    cat("\n  global decisions:\n")
    print_preflight_global_block(x$global)
  }
  if (level %in% c("dim", "both")) {
    cat("\n  per-coordinate decisions:\n")
    print_preflight_per_dim_block(x$per_dim)
  }
  if (level == "global") {
    cat("\n  use print(x, level = \"dim\") for per-coordinate detail",
        " or as.data.frame(x) for the tidy table.\n", sep = "")
  }
  invisible(x)
}


#' @keywords internal
#' @noRd
format_preflight_table <- function(cols, indent = "    ") {
  if (length(cols) == 0L) return(character())
  n_rows <- length(cols[[1L]])
  if (n_rows == 0L) return(character())
  headers <- names(cols)
  body <- vapply(seq_along(cols), function(j) {
    val <- cols[[j]]
    chr <- if (is.numeric(val)) formatC(val, format = "g", digits = 4)
           else as.character(val)
    chr[is.na(val)] <- "NA"
    chr
  }, character(n_rows))
  if (n_rows == 1L) body <- matrix(body, nrow = 1L)
  if (!is.matrix(body)) body <- matrix(body, nrow = n_rows)
  widths <- vapply(seq_along(headers), function(j) {
    max(nchar(c(headers[j], body[, j])))
  }, integer(1L))
  pad <- function(x, w) formatC(x, width = w, flag = "-")
  header_line <- paste(mapply(pad, headers, widths,
                              SIMPLIFY = TRUE), collapse = "  ")
  body_lines <- vapply(seq_len(n_rows), function(i) {
    paste(mapply(pad, body[i, ], widths, SIMPLIFY = TRUE),
          collapse = "  ")
  }, character(1L))
  paste0(indent, c(header_line, body_lines))
}


#' @keywords internal
#' @noRd
print_preflight_global_block <- function(global) {
  if (nrow(global) == 0L) {
    cat("    (no components)\n")
    return(invisible(NULL))
  }
  flag <- ifelse(
    is.na(global$agreement) | !is.finite(global$agreement), "",
    ifelse(global$agreement >= 1 - 1e-12, "uniform",
           ifelse(global$agreement <= 0.5 + 1e-12, "split",
                  "mixed"))
  )
  cols <- list(
    component       = global$component,
    global_decision = global$global_decision,
    agreement       = formatC(global$agreement, format = "f",
                              digits = 2),
    flag            = flag
  )
  lines <- format_preflight_table(cols)
  cat(lines, sep = "\n")
  cat("\n")
}


#' @keywords internal
#' @noRd
print_preflight_per_dim_block <- function(per_dim) {
  if (nrow(per_dim) == 0L) {
    cat("    (no per-dim rows)\n")
    return(invisible(NULL))
  }
  cols <- list(
    component       = per_dim$component,
    dim             = per_dim$dim,
    decision        = per_dim$decision,
    decision_reason = per_dim$decision_reason,
    t_info_cp       = formatC(per_dim$t_info_cp, format = "g",
                              digits = 3),
    t_info_ncp      = formatC(per_dim$t_info_ncp, format = "g",
                              digits = 3)
  )
  lines <- format_preflight_table(cols)
  cat(lines, sep = "\n")
  cat("\n")
}


#' Summary method for gdpar_preflight_report
#'
#' Returns a list with the per-component aggregated table, the
#' overall agreement (mean of the per-component agreement, NA
#' values dropped), the number of components, the number of
#' coordinates, and a count of per-dim CP/NCP/absent decisions.
#'
#' @param object An object of class \code{gdpar_preflight_report}.
#' @param ... Unused; present for S3 generic compatibility.
#' @return Named list as described above.
#' @export
summary.gdpar_preflight_report <- function(object, ...) {
  assert_inherits(object, "gdpar_preflight_report", "object")
  finite_agree <- object$global$agreement[
    is.finite(object$global$agreement)
  ]
  overall_agreement <- if (length(finite_agree) > 0L) {
    mean(finite_agree)
  } else {
    NA_real_
  }
  list(
    n_components      = nrow(object$global),
    n_dims            = length(unique(object$per_dim$dim)),
    method            = object$method,
    per_component     = object$global,
    overall_agreement = overall_agreement,
    per_dim_counts    = c(
      CP     = sum(object$per_dim$decision == "CP"),
      NCP    = sum(object$per_dim$decision == "NCP"),
      absent = sum(object$per_dim$decision == "absent")
    ),
    settings          = object$settings
  )
}


#' as.data.frame method for gdpar_preflight_report
#'
#' Returns the tidy per-dim data frame. Suitable for use with
#' \code{subset}, \code{aggregate}, \code{dplyr::filter}, etc.
#'
#' @param x An object of class \code{gdpar_preflight_report}.
#' @param row.names \code{NULL} or character vector forwarded as
#'   row.names of the returned data frame.
#' @param optional Ignored; present for S3 generic compatibility.
#' @param ... Unused; present for S3 generic compatibility.
#' @return Data frame with one row per (component, dim).
#' @export
as.data.frame.gdpar_preflight_report <- function(x,
                                                 row.names = NULL,
                                                 optional = FALSE,
                                                 ...) {
  assert_inherits(x, "gdpar_preflight_report", "x")
  df <- x$per_dim
  if (!is.null(row.names)) rownames(df) <- row.names
  df
}


#' Format method for gdpar_preflight_report
#'
#' One-line representation of the report. Used by \code{format} in
#' contexts where a print is not desired.
#'
#' @param x An object of class \code{gdpar_preflight_report}.
#' @param ... Unused; present for S3 generic compatibility.
#' @return Character scalar.
#' @export
format.gdpar_preflight_report <- function(x, ...) {
  n_dims <- length(unique(x$per_dim$dim))
  sprintf(
    "<gdpar_preflight_report: p=%d, aggregation=%s, components=%s>",
    n_dims,
    x$method,
    paste(unique(x$per_dim$component), collapse = ",")
  )
}

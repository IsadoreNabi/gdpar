#' Broadcast a uniform per-component specification to multiple dimensions
#'
#' Construct a \code{dims_spec} object that records the same additive and
#' multiplicative formulas for every dimension \eqn{k = 1, \ldots, p} of a
#' multivariate \code{theta_i} vector. Used as the value of the \code{dims}
#' argument of \code{\link{amm_spec}} when all dimensions share the same
#' per-component specification. Per-dimension overrides can be layered on
#' top via \code{\link{override}}.
#'
#' @param a One-sided formula declaring the additive basis applied
#'   uniformly to every dimension of \code{theta_i}. Use \code{NULL} to
#'   disable the additive component for all dimensions.
#' @param b One-sided formula declaring the multiplicative basis applied
#'   uniformly to every dimension of \code{theta_i}. Use \code{NULL} to
#'   disable the multiplicative component for all dimensions.
#'
#' @return An object of class \code{dims_spec} with components \code{base}
#'   (the uniform template, a list with \code{a} and \code{b}) and
#'   \code{overrides} (an initially empty named list of per-dimension
#'   overrides keyed by character form of the integer index \code{k}).
#'
#' @details
#' This helper enables a low-verbosity declaration of \code{p > 1} specs
#' for the common case where all dimensions of the multivariate
#' \code{theta_i} share the same \code{a_k} and \code{b_k} structures.
#' Broadcasting is explicit: the user must wrap the per-component
#' formulas in \code{dimwise()} to opt in. Bare formulas passed directly
#' to the \code{dims} argument of \code{\link{amm_spec}} when
#' \code{p > 1} are rejected, to avoid silent recycling masking bugs.
#'
#' The dimension \code{p} is not stored in the \code{dims_spec}; it is
#' resolved at the point of consumption by \code{\link{amm_spec}}, which
#' takes \code{p} as an explicit argument and validates coherence across
#' the spec, the multivariate W basis, and any overrides.
#'
#' The modulating component \code{W} is not part of \code{dims_spec}
#' because \code{W} couples all dimensions of \code{theta_ref} in the
#' canonical AMM form and therefore cannot be declared per-dimension.
#' \code{W} is supplied to \code{\link{amm_spec}} as a separate top-level
#' argument via a multivariate basis.
#'
#' @section Methodological notes:
#' The asymmetry between per-dimension components (\code{a}, \code{b})
#' and the cross-dimension component (\code{W}) reflects the canonical
#' AMM form
#' \deqn{\theta_i[k] = \theta_{\mathrm{ref}}[k] + a_k(x_i)
#'       + b_k(x_i) \theta_{\mathrm{ref}}[k]
#'       + (W_k(\theta_{\mathrm{ref}}) - W_k(\theta_{\mathrm{anchor}})) x_i,
#'       \quad k = 1, \ldots, p,}
#' in which \eqn{a_k} and \eqn{b_k} depend only on the covariates
#' \eqn{x_i} (and therefore factor per \eqn{k}) while each \eqn{W_k}
#' depends on the full vector \eqn{\theta_{\mathrm{ref}}} (and therefore
#' couples the dimensions). Representing \code{W} as a list of
#' independent per-dimension bases would silently restrict the model
#' class to the separable sub-class. \code{dims_spec} therefore excludes
#' \code{W} by construction.
#'
#' @examples
#' base <- dimwise(a = ~ x1 + x2, b = ~ x1)
#' print(base)
#'
#' with_override <- override(base, k = 2L, a = ~ x1)
#' print(with_override)
#'
#' @seealso \code{\link{override}}, \code{\link{amm_spec}}
#' @export
dimwise <- function(a = NULL, b = NULL) {
  assert_one_sided_formula(a, "a", allow_null = TRUE)
  assert_one_sided_formula(b, "b", allow_null = TRUE)
  obj <- list(
    base = list(a = a, b = b),
    overrides = list()
  )
  class(obj) <- c("dims_spec", "list")
  obj
}

#' Layer a per-dimension override on a dims_spec
#'
#' Given a \code{dims_spec} produced by \code{\link{dimwise}}, attach a
#' per-dimension override that replaces the additive and/or multiplicative
#' formula for a specific dimension index \code{k}.
#'
#' @param dims A \code{dims_spec} object produced by \code{\link{dimwise}}.
#' @param k Integer scalar with the dimension index to override. Must be
#'   a positive integer. Coherence with the global dimension \code{p} is
#'   checked at the point of consumption by \code{\link{amm_spec}}.
#' @param a Optional one-sided formula replacing the additive basis for
#'   dimension \code{k}. Pass \code{NULL} explicitly to disable the
#'   additive component for that dimension while keeping it active
#'   elsewhere. Omit to leave it unchanged from the base.
#' @param b Optional one-sided formula replacing the multiplicative
#'   basis for dimension \code{k}. Same semantics as \code{a}.
#'
#' @return A new \code{dims_spec} with the override registered. Multiple
#'   calls to \code{override} compose; calling \code{override} twice with
#'   the same \code{k} replaces the previous override for that index.
#'
#' @details
#' Overrides are recorded by integer index, not by position in a list,
#' so the order of \code{override} calls does not matter beyond the
#' overwrite semantics noted above.
#'
#' The semantics of "unchanged" versus "disabled" requires distinguishing
#' between an argument that is missing from the call and an argument
#' that is explicitly \code{NULL}. The function uses
#' \code{\link[base]{missing}} for this distinction: omit the argument
#' to inherit from the base; pass \code{NULL} to disable for this
#' dimension only.
#'
#' At least one of \code{a} or \code{b} must be supplied; calling
#' \code{override} without any change is treated as a user error and
#' aborts with an informative message.
#'
#' @examples
#' base <- dimwise(a = ~ x1 + x2, b = ~ x1)
#' v1 <- override(base, k = 2L, a = ~ x1)
#' v2 <- override(v1, k = 3L, b = NULL)
#' print(v2)
#'
#' @seealso \code{\link{dimwise}}, \code{\link{amm_spec}}
#' @export
override <- function(dims, k, a, b) {
  assert_inherits(dims, "dims_spec", "dims")
  assert_count(k, "k")
  if (missing(a) && missing(b)) {
    gdpar_abort(
      paste(
        "override(): at least one of 'a' or 'b' must be supplied.",
        "To leave a dimension unchanged, do not call override() for it."
      ),
      class = "gdpar_input_error"
    )
  }
  k_int <- as.integer(k)
  key <- as.character(k_int)
  ov <- dims$overrides[[key]]
  if (is.null(ov)) {
    ov <- list(a = NULL, b = NULL, a_set = FALSE, b_set = FALSE)
  }
  if (!missing(a)) {
    assert_one_sided_formula(a, "a", allow_null = TRUE)
    ov["a"] <- list(a)
    ov$a_set <- TRUE
  }
  if (!missing(b)) {
    assert_one_sided_formula(b, "b", allow_null = TRUE)
    ov["b"] <- list(b)
    ov$b_set <- TRUE
  }
  dims$overrides[[key]] <- ov
  dims
}

#' Internal: resolve a dims_spec into a flat per-dimension list
#'
#' Given a \code{dims_spec} and an integer \code{p}, produce the canonical
#' per-dimension representation: a list of length \code{p}, each element
#' a list with \code{a} and \code{b} components. Overrides are applied on
#' top of the base.
#'
#' @param dims A \code{dims_spec} object.
#' @param p Positive integer with the global dimension.
#'
#' @return A list of length \code{p}. Each entry is a list with two
#'   components, \code{a} and \code{b}, each either a one-sided formula
#'   or \code{NULL}.
#' @keywords internal
#' @noRd
resolve_dims_spec <- function(dims, p) {
  assert_inherits(dims, "dims_spec", "dims")
  assert_count(p, "p")
  p <- as.integer(p)
  ovs <- dims$overrides
  if (length(ovs) > 0L) {
    bad <- character(0)
    for (key in names(ovs)) {
      k_idx <- suppressWarnings(as.integer(key))
      if (is.na(k_idx) || k_idx < 1L || k_idx > p) {
        bad <- c(bad, key)
      }
    }
    if (length(bad) > 0L) {
      gdpar_abort(
        sprintf(
          paste(
            "Override(s) declared for out-of-range dimension index/indices:",
            "%s. Valid range is 1:%d."
          ),
          paste(bad, collapse = ", "), p
        ),
        class = "gdpar_input_error",
        data = list(bad_keys = bad, p = p)
      )
    }
  }
  out <- vector("list", p)
  for (k in seq_len(p)) {
    a_k <- dims$base$a
    b_k <- dims$base$b
    ov <- ovs[[as.character(k)]]
    if (!is.null(ov)) {
      if (isTRUE(ov$a_set)) a_k <- ov[["a"]]
      if (isTRUE(ov$b_set)) b_k <- ov[["b"]]
    }
    out[[k]] <- list(a = a_k, b = b_k)
  }
  out
}

#' Print method for dims_spec objects
#'
#' @param x A \code{dims_spec} object.
#' @param ... Unused; present for S3 generic compatibility.
#' @return Invisibly returns \code{x}.
#' @export
print.dims_spec <- function(x, ...) {
  cat("<dims_spec>\n")
  cat("  base:\n")
  cat("    a : ",
      if (is.null(x$base$a)) "NULL" else deparse(x$base$a),
      "\n", sep = "")
  cat("    b : ",
      if (is.null(x$base$b)) "NULL" else deparse(x$base$b),
      "\n", sep = "")
  if (length(x$overrides) > 0L) {
    cat("  overrides:\n")
    keys <- sort(as.integer(names(x$overrides)))
    for (k in keys) {
      ov <- x$overrides[[as.character(k)]]
      parts <- character(0)
      if (isTRUE(ov$a_set)) {
        a_val <- ov[["a"]]
        parts <- c(parts, paste0("a = ",
                                 if (is.null(a_val)) "NULL" else deparse(a_val)))
      }
      if (isTRUE(ov$b_set)) {
        b_val <- ov[["b"]]
        parts <- c(parts, paste0("b = ",
                                 if (is.null(b_val)) "NULL" else deparse(b_val)))
      }
      cat("    k = ", k, " : ", paste(parts, collapse = "; "),
          "\n", sep = "")
    }
  } else {
    cat("  overrides: <none>\n")
  }
  invisible(x)
}
